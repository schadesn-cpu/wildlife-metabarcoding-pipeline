#!/usr/bin/env python3
"""diagnose_analysis_table_join.py

I measure what putting build_analysis_table.py on the shared library lookup
would change, before anything is patched. I write nothing.

WHY
---
build_analysis_table.py builds its own cwmap from name_16S. It is the most
careful of the private resolvers, and its comment says so: it casefolds and
strips the marker tag by regex, after the join "has been got wrong three times:
2026-08-05 x2, 2026-08-06". But it is still private, and it registers only
three forms per library name: the name itself, the marker-tag-stripped stem,
and stem + '-16s'.

That does not cover the replicate suffix. MARKER_TAG is

    -(16s|18s|cytob|cb)([A-Za-z]*\\d*)(-\\d+)?$

which strips '-16s-2' but leaves a bare '-2' alone, so 'EQ02441291-2' stays
'EQ02441291-2' and only matches if the crosswalk carries that exact string.
Section 12 of the ground truth records 122 libraries in this very table with
no crosswalk match, which are "NOT controls, they are the replicate-suffix
shape". This is the script that count came from.

WHAT I REPORT
-------------
For every sample in the decontam 16S feature table, whether cwmap resolves it,
whether the lookup resolves it, and where they disagree. I split the current
'no_host_or_pool' exclusion into its two very different causes:

    join failure        the name resolved to nothing   -> a bug
    incomplete record   it resolved, but the crosswalk has no host or no
                        tick count for that specimen   -> a data gap, and the
                        script is right to exclude it

Those are lumped together today, so a join failure is invisible inside a count
that also contains legitimate exclusions.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 diagnose_analysis_table_join.py
"""
import glob
import os
import re
import sys
import tempfile
import zipfile

import pandas as pd

sys.path.insert(0, "scripts")
try:
    import tick_lookup
except ImportError:
    sys.exit("could not import scripts/tick_lookup.py. Run this from the "
             "project root.")

try:
    import biom
except ImportError:
    sys.exit("biom is not importable. Run this inside the qiime2 environment, "
             "the same one build_analysis_table.py runs in.")

TABLE = "results/decontam/16S/table_decontam_16S.qza"
CW = "results/metadata/crosswalk_4marker.tsv"
THRESH = 1000
CONTROL_RE = re.compile(r"(XB|NTC|POS|PC-WELL|EXTRACTIONBL|XB-OR|BLANK)", re.I)
MARKER_TAG = re.compile(r"-(16s|18s|cytob|cb)([A-Za-z]*\d*)(-\d+)?$", re.I)


def load_table(qza):
    """I unpack a QIIME artifact and return its biom table, as the source does."""
    workdir = tempfile.mkdtemp()
    with zipfile.ZipFile(qza) as archive:
        archive.extractall(workdir)
    hits = glob.glob(workdir + "/*/data/feature-table.biom")
    if not hits:
        sys.exit("no feature-table.biom inside " + qza)
    return biom.load_table(hits[0])


def main():
    for path in (TABLE, CW):
        if not os.path.exists(path):
            sys.exit("not found: " + path)

    crosswalk = pd.read_csv(CW, sep="\t")
    lookup = tick_lookup.load()

    # cwmap, copied from build_analysis_table.py so this measures that script
    # and not my idea of it.
    cwmap = {}
    for _, row in crosswalk.iterrows():
        record = {"Barcode": row["sample_id"], "Host": row["host_spp"],
                  "NumTicks": row["num_ticks"]}
        for name in str(row["name_16S"]).split(";"):
            name = name.strip()
            if not name or name.lower() == "nan":
                continue
            stem = MARKER_TAG.sub("", name)
            for form in (name, stem, stem + "-16s"):
                cwmap[form.lower()] = record
    print("crosswalk: %d specimens -> %d cwmap library keys"
          % (len(crosswalk), len(cwmap)))

    by_specimen = crosswalk.set_index("sample_id")

    table = load_table(TABLE)
    samples = list(table.ids("sample"))
    print("feature table: %d samples" % len(samples))
    print()

    counts = {"control": 0, "fungal_trial": 0, "below_thresh": 0,
              "qualifying": 0}
    old_hit = new_hit = 0
    rescued, disagree, incomplete = [], [], []

    for sample in samples:
        base = str(sample).split(".")[0]
        if CONTROL_RE.search(base):
            counts["control"] += 1
            continue
        if "VT-12" in base:
            counts["fungal_trial"] += 1
            continue
        if table.data(sample, "sample").sum() < THRESH:
            counts["below_thresh"] += 1
            continue
        counts["qualifying"] += 1

        by_script = (cwmap.get(str(sample).lower())
                     or cwmap.get(str(base).lower()))
        by_lookup = (tick_lookup.resolve(lookup, "16S", sample)
                     or tick_lookup.resolve(lookup, "16S", base)
                     or tick_lookup.resolve(lookup, "specimen", base))

        if by_script is not None:
            old_hit += 1
        if by_lookup is not None:
            new_hit += 1

        if by_script is None and by_lookup is not None:
            rescued.append((base, by_lookup))
        elif (by_script is not None and by_lookup is not None
              and by_script["Barcode"] != by_lookup):
            disagree.append((base, by_script["Barcode"], by_lookup))

        # of the ones that DO resolve, which are excluded for a real data gap
        specimen = by_lookup or (by_script["Barcode"] if by_script else None)
        if specimen is not None and specimen in by_specimen.index:
            row = by_specimen.loc[specimen]
            if pd.isna(row["host_spp"]) or pd.isna(row["num_ticks"]):
                incomplete.append((base, specimen,
                                   "no host" if pd.isna(row["host_spp"])
                                   else "no tick count"))

    print("=== sample disposition, before any join ===")
    for key in ("control", "fungal_trial", "below_thresh", "qualifying"):
        print("  %-14s %d" % (key, counts[key]))
    print()
    print("=== of the %d qualifying samples ===" % counts["qualifying"])
    print("  resolved by build_analysis_table.py's cwmap: %d" % old_hit)
    print("  resolved by the shared library lookup:       %d" % new_hit)
    print()

    if rescued:
        print("RESCUED by the lookup: %d sample(s) the cwmap drops" % len(rescued))
        specimens = sorted({s for _, s in rescued})
        print("  %d distinct specimens" % len(specimens))
        flags = crosswalk[crosswalk["sample_id"].isin(specimens)]
        single = flags["is_single_tick"].astype(str).str.lower()
        print("  of which flagged is_single_tick: %d"
              % int(single.isin(["true", "1"]).sum()))
        print()
        print("  %-42s %s" % ("library", "specimen"))
        for lib, spec in rescued[:40]:
            print("  %-42s %s" % (lib, spec))
        if len(rescued) > 40:
            print("  ... and %d more" % (len(rescued) - 40))
        print()
    else:
        print("NOTHING RESCUED. The cwmap resolves everything the lookup does.")
        print()

    if disagree:
        print("DISAGREEMENTS: %d sample(s) resolve to DIFFERENT specimens"
              % len(disagree))
        for lib, a, b in disagree[:40]:
            print("  %-42s cwmap %s  lookup %s" % (lib, a, b))
        print("  These matter more than the rescues. Do not patch anything")
        print("  until each one is understood.")
        print()

    if incomplete:
        specimens = sorted({s for _, s, _ in incomplete})
        print("EXCLUDED FOR A REAL DATA GAP, not a join failure: %d sample(s), "
              "%d specimen(s)" % (len(incomplete), len(specimens)))
        print("  These resolve fine; the crosswalk simply has no host or no")
        print("  tick count for them, and the script is right to drop them.")
        print("  Today they are counted in 'no_host_or_pool' alongside genuine")
        print("  join failures, which is why a join bug can hide there.")
        for lib, spec, why in incomplete[:20]:
            print("  %-42s %-28s %s" % (lib, spec, why))
        if len(incomplete) > 20:
            print("  ... and %d more" % (len(incomplete) - 20))
        print()

    print("=== verdict ===")
    if not rescued and not disagree:
        print("The join is already complete. Putting this script on the lookup")
        print("is tidying, and analysis_table_16S.tsv will not change.")
    else:
        print("The join is NOT complete. analysis_table_16S.tsv is short by")
        print("%d sample(s), and every pooled analysis downstream of it is too:"
              % len(rescued))
        print("PERMANOVA, PERMDISP, pathogen_cooccurrence.py, fig1a, fig1b.")


if __name__ == "__main__":
    main()
