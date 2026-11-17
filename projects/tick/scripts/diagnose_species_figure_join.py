#!/usr/bin/env python3
"""diagnose_species_figure_join.py

I check whether build_species_figure.py's private library-to-specimen map is
losing real specimens. I change nothing and write nothing.

WHY
---
build_species_figure.py builds its own lib2spec by splitting name_18S on
semicolons and swapping a literal "-18s" suffix, then reports every miss as
"(controls and excluded specimens — expected)". That is a hardcoded claim
about data nobody has looked at, and it is word for word the claim that was
removed from make_specimen_taxon_table.py on 2026-09-09 after it turned out to
be hiding 15 real specimens, 2 of them single ticks.

TICK_GROUND_TRUTH section 3 also records that 199 18S library names carry a
variable-region suffix such as -18sV11, which a literal "-18s" swap does not
produce.

So the question is whether 316, 265 and 20 are complete. This script answers it
by resolving the same libraries through results/metadata/library_lookup.tsv,
which build_4marker_crosswalk.py now emits, and printing the difference.

WHAT TO DO WITH THE ANSWER
--------------------------
If "libraries the lookup resolves and lib2spec does not" is 0, the species
numbers are complete and swapping the resolver is tidying that can wait.

If it is more than 0, those are specimens missing from the species figure, the
locked counts move, and it needs its own change rather than being folded into
a cosmetic patch.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 diagnose_species_figure_join.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
try:
    import tick_lookup
except ImportError:
    sys.exit("could not import scripts/tick_lookup.py. Run this from the "
             "project root, with tick_lookup.py in scripts/.")

TG = "results/18S/all/taxonomy_verified/tick_genus_per_sample_corrected.tsv"
CW = "results/metadata/crosswalk_4marker.tsv"


def main():
    for path in (TG, CW):
        if not os.path.exists(path):
            sys.exit("not found: " + path)

    tick_genus = pd.read_csv(TG, sep="\t")
    crosswalk = pd.read_csv(CW, sep="\t")
    lookup = tick_lookup.load()

    # lib2spec, copied verbatim from build_species_figure.py so this measures
    # that script and not my idea of it.
    lib2spec = {}
    for _, row in crosswalk.iterrows():
        for name in str(row["name_18S"]).split(";"):
            name = name.strip()
            if not name or name.lower() == "nan":
                continue
            for form in {name,
                         name[:-4] if name.endswith("-18s") else name + "-18s"}:
                lib2spec[form] = row["sample_id"]

    libs = (tick_genus["sample_id"].astype(str)
            .str.replace(r"\..*$", "", regex=True))

    old_hit, new_hit, rows = 0, 0, []
    for lib in libs:
        by_script = lib2spec.get(lib)
        by_lookup = (tick_lookup.resolve(lookup, "18S", lib)
                     or tick_lookup.resolve(lookup, "specimen", lib))
        if by_script is not None:
            old_hit += 1
        if by_lookup is not None:
            new_hit += 1
        if by_script is None and by_lookup is not None:
            rows.append((lib, by_lookup, "MISSED by lib2spec"))
        elif (by_script is not None and by_lookup is not None
              and by_script != by_lookup):
            rows.append((lib, by_script + " vs " + by_lookup, "DISAGREE"))

    print("18S tick-genus libraries: " + str(len(libs)))
    print("  resolved by build_species_figure.py's lib2spec: " + str(old_hit))
    print("  resolved by the shared library lookup:          " + str(new_hit))
    print()

    if not rows:
        print("NO DIFFERENCE. Every library the lookup resolves, lib2spec also")
        print("resolves, and they never disagree. 316 / 265 / 20 are complete;")
        print("swapping the resolver is tidying, not a correction.")
        return

    print("DIFFERENCES: " + str(len(rows)))
    print("These are libraries the species figure is dropping or mis-assigning.")
    print()
    print("%-42s %-34s %s" % ("library", "specimen", "problem"))
    for lib, spec, why in rows[:60]:
        print("%-42s %-34s %s" % (lib, spec, why))
    if len(rows) > 60:
        print("... and " + str(len(rows) - 60) + " more")

    missed = sorted({spec for lib, spec, why in rows if why.startswith("MISSED")})
    if missed:
        single = crosswalk[crosswalk["sample_id"].isin(missed)]
        flags = single["is_single_tick"].astype(str).str.lower()
        print()
        print("distinct specimens missed: " + str(len(missed)))
        print("  of which flagged is_single_tick: "
              + str(int(flags.isin(["true", "1"]).sum())))
        print()
        print("The species counts are NOT complete. Do not treat 316, 265 or 20")
        print("as locked until this is resolved.")


if __name__ == "__main__":
    main()
