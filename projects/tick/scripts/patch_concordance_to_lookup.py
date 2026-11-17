#!/usr/bin/env python3
"""patch_concordance_to_lookup.py

I put build_concordance.py on the shared library lookup, which takes its
cohort from 83 to 89.

THE BUG BEING REMOVED
---------------------
build_concordance.py's norm() stripped hyphens along with spaces and
underscores:

    'ALHT WTD_09'   -> 'alhtwtd09'
    'alht-deer-09'  -> 'alhtdeer09'

Two keys, one specimen. Six single ticks carry the 'ALHT WTD_nn' shape while
the qPCR sheet keys them as 'alht-deer-nn', so all six failed to pair and
the script reported 83 of 93 screened by both methods with no warning that
the ten it dropped included six it should have kept. ALHT WTD_65 is one of
the six, and it is one of the three ticks whose Anaplasma call changed in
the genus anchoring, so it lands in a cell that just moved.

WHAT REPLACES IT
----------------
tick_lookup, which reads results/metadata/library_lookup.tsv. No parsing in
this script at all: a name either is in the table or it is not, and a miss
returns None rather than a partial match.

I also record the cohort split, because 89 and not 93 is now a deliberate
number rather than an accident. The four single ticks that cannot pair are
USGSTICKS-4103, 4108, 4118 and 7309, whose only qPCR rows are pools of 2, 4,
3 and 2 subsamples. A pooled qPCR result cannot be compared per tick.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_concordance_to_lookup.py
    python3 scripts/build_concordance.py

I refuse to run twice and refuse to guess if any anchor has moved.
"""
import os
import sys

TARGET = "scripts/build_concordance.py"
NEEDS = "scripts/tick_lookup.py"

OLD_NORM = '''def norm(v):
    """Normalise an identifier for joining: lowercase, drop spaces/underscores/hyphens."""
    return re.sub(r'[\\s_\\-]', '', str(v).strip().lower())
'''

NEW_NORM = '''# Name resolution is NOT done here. It is done once, by
# build_4marker_crosswalk.py, and read from results/metadata/library_lookup.tsv.
# The previous norm() in this slot stripped hyphens as well as spaces, so
# 'ALHT WTD_09' and the qPCR sheet's 'alht-deer-09' became two different keys
# and six single ticks silently failed to pair. Do not reintroduce a local
# normaliser; if a name does not resolve, fix the crosswalk.
import tick_lookup

LOOKUP = tick_lookup.load()


def specimen_of(value):
    """I resolve one accessioning ID, barcode or qPCR key to a specimen."""
    return tick_lookup.resolve(LOOKUP, "specimen", value)
'''

OLD_JOIN = '''# qPCR has BOTH a Barcode column (sometimes blank) and host_id. Key on either.
q1["k1"] = q1["Barcode"].map(norm)
q1["k2"] = q1["host_id"].map(norm)
qmap = {}
for _, r in q1.iterrows():
    for k in (r["k1"], r["k2"]):
        if k and k != "nan":
            qmap[k] = r

st["k"] = st["Barcode"].map(norm)
st["_q"] = st["k"].map(lambda k: qmap.get(k))
'''

NEW_JOIN = '''# qPCR has BOTH a Barcode column (sometimes blank) and host_id. Either one
# resolves through the lookup, so both are tried and the result is keyed on the
# SPECIMEN rather than on whichever string that row happened to carry.
qmap = {}
qclash = []
for _, r in q1.iterrows():
    hit = specimen_of(r["Barcode"]) or specimen_of(r["host_id"])
    if hit is None:
        continue
    if hit in qmap:
        qclash.append(hit)
    else:
        qmap[hit] = r

say(f"qPCR single-subsample rows resolving to a specimen: {len(qmap)}")
if qclash:
    say(f"  WARN {len(qclash)} specimen(s) have more than one single-subsample "
        f"qPCR row; the first was kept: {sorted(set(qclash))}")

st["_specimen"] = st["Barcode"].map(specimen_of)
_unresolved = st.loc[st["_specimen"].isna(), "Barcode"].astype(str).tolist()
if _unresolved:
    say(f"  WARN {len(_unresolved)} sequencing specimen(s) are not in the "
        f"lookup at all: {_unresolved}")
    say("       That is a crosswalk problem, not a join problem. Rebuild it.")
st["_q"] = st["_specimen"].map(lambda s: qmap.get(s) if s is not None else None)
'''

OLD_UNPAIRED = '''unpaired = sorted(st.loc[st["_q"].isna(), "Barcode"].astype(str))
say(f"in the sequencing set but NOT qPCR-screened: {len(unpaired)}")
'''

NEW_UNPAIRED = '''unpaired = sorted(st.loc[st["_q"].isna(), "Barcode"].astype(str))
say(f"in the sequencing set but with NO single-subsample qPCR row: {len(unpaired)}")
say("  These specimens ARE in the qPCR sheet. Their only rows are pools, and a")
say("  pooled qPCR result cannot be compared against a single tick. This is the")
say("  reason the cohort is smaller than the sequencing set, and it is the only")
say("  reason: every single tick resolves against the qPCR sheet.")
'''

OLD_DOC = ('''The qPCR-vs-metabarcoding concordance existed ONLY as a PNG (78 / 90 / 86%, n=93 as
then reported; the current values are 79.5 / 90.4 / 85.5% on n=83).''')

NEW_DOC = ('''The qPCR-vs-metabarcoding concordance existed ONLY as a PNG (78 / 90 / 86%, n=93
as then reported). The values were later recomputed as 79.5 / 90.4 / 85.5% on
n=83, and that n was wrong: a local key normaliser stripped hyphens and lost the
six 'ALHT WTD_nn' specimens. On the shared lookup the cohort is n=89, being the
93 sequenced single ticks minus the 4 whose only qPCR rows are pools.''')

OLD_READS = '''Reads:
    results/single_tick_v3/single_tick_v3_full.tsv
    metadata/tick_full_metadata.xlsx  sheet 'qPCR Results\''''

NEW_READS = '''Reads:
    results/metadata/library_lookup.tsv   (via tick_lookup; the ONLY name resolver)
    results/single_tick_v3/single_tick_v3_full.tsv
    metadata/tick_full_metadata.xlsx  sheet 'qPCR Results\''''

REPLACEMENTS = [
    (OLD_DOC, NEW_DOC, "the docstring cohort note"),
    (OLD_READS, NEW_READS, "the Reads: block"),
    (OLD_NORM, NEW_NORM, "the norm() definition"),
    (OLD_JOIN, NEW_JOIN, "the qPCR join block"),
    (OLD_UNPAIRED, NEW_UNPAIRED, "the unpaired report"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    if not os.path.isfile(NEEDS):
        sys.exit("not found: " + NEEDS
                 + "\nWrite tick_lookup.py first; this patch imports it.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "import tick_lookup" in source:
        sys.exit("already patched: " + TARGET + " imports tick_lookup. "
                 "Nothing changed.")

    for old, new, label in REPLACEMENTS:
        count = source.count(old)
        if count != 1:
            sys.exit("expected exactly one occurrence of " + label
                     + ", found " + str(count)
                     + ". The script has changed since this patch was written; "
                       "patch by hand rather than trusting this.")
        source = source.replace(old, new, 1)

    compile(source, TARGET, "exec")
    with open(TARGET, "w", encoding="utf-8") as handle:
        handle.write(source)
    print("patched " + TARGET)
    print("now run: python3 " + TARGET)
    print("expect: 89 specimens screened by both methods")


if __name__ == "__main__":
    main()
