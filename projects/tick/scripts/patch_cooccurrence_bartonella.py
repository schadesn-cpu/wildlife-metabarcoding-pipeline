#!/usr/bin/env python3
"""patch_cooccurrence_bartonella.py

I bring pathogen_cooccurrence.py into line with the rest of the project: it
still counts Bartonella as a pathogen, and its closing caveat still says the
single-tick cohort is n=91.

WHAT CHANGES
------------
1. BARTONELLA leaves the pathogen list. Ground truth 1.16: Section 7 reports
   these as ruminant-associated and ked-associated organisms outside the
   pathogen results, and 6 of the 10 Bartonella-positive single ticks carry
   deer ked reads (odds ratio 23.4, Fisher p = 0.0001).
   build_single_tick_set.py and build_analysis_table.py were both changed on
   2026-09-18; this script was the third place the old classification lived.

   Bartonella_pos is still WRITTEN to analysis_table_16S.tsv, so it is still
   available; it is simply not tested for co-occurrence with pathogens.

2. THE LIST IS NOW CHECKED AGAINST THE TABLE. This script hardcodes a pathogen
   list that must match PATHOGENS in build_analysis_table.py, and nothing
   enforced that. It now asserts that every name it expects has a _pos column,
   and that Bartonella_pos is present but deliberately unused. If either
   changes, the script stops rather than quietly testing a different set.

3. THE CLOSING CAVEAT says single-tick co-infection is n=91. It is 93, and has
   been since 2026-08-07. It also now names the current figure and says
   Bartonella is excluded there too, so the two co-occurrence analyses read as
   one decision rather than two.

WHAT TO EXPECT
--------------
The sample count moves from 558 to 602 (ground truth 1.20), so every cell
changes. The pairs drop from 6 to 3, because removing Bartonella removes the
three pairs it was in and Ehrlichia is all-zero and already dropped.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_cooccurrence_bartonella.py
    python3 scripts/pathogen_cooccurrence.py
"""
import os
import sys

TARGET = "scripts/pathogen_cooccurrence.py"

OLD_LIST = """d = pd.read_csv('results/analysis/analysis_table_16S.tsv', sep='\\t')
pathogens = ['Anaplasma','Borreliella','Borrelia','Bartonella','Ehrlichia']
poscols = {p: p+'_pos' for p in pathogens if p+'_pos' in d.columns}"""

NEW_LIST = """d = pd.read_csv('results/analysis/analysis_table_16S.tsv', sep='\\t')

# This list MUST match PATHOGENS in build_analysis_table.py. Nothing enforced
# that until 2026-09-18, which is how Bartonella stayed a pathogen here after
# it had been reclassified in build_single_tick_set.py.
pathogens = ['Anaplasma','Borreliella','Borrelia','Ehrlichia']

# Bartonella is a CO-DETECTION, not a pathogen. Ground truth 1.16: these are
# ruminant-associated and ked-associated organisms reported outside the
# pathogen results, and 6 of the 10 Bartonella-positive single ticks carry deer
# ked reads (OR 23.4, Fisher p = 0.0001). The column is still written by
# build_analysis_table.py; it is simply not tested here.
co_detections = ['Bartonella']

missing = [p for p in pathogens if p + '_pos' not in d.columns]
if missing:
    raise SystemExit(
        "analysis_table_16S.tsv has no _pos column for: " + ", ".join(missing)
        + "\\nThis list must match PATHOGENS in build_analysis_table.py. "
          "Rebuild the table or fix the list; do not silently test a "
          "different set.")
for c in co_detections:
    if c + '_pos' not in d.columns:
        print("note: " + c + "_pos is not in the table. It is excluded here by "
              "classification, not by absence; check build_analysis_table.py.")
    else:
        print("excluded as a co-detection, not a pathogen: " + c
              + " (" + str(int(d[c + '_pos'].sum())) + "/" + str(len(d)) + ")")

poscols = {p: p+'_pos' for p in pathogens if p+'_pos' in d.columns}"""

OLD_CAVEAT = """print("\\nCAVEAT: co-occurrence in a POOL = both detected in the same pool, NOT proof of co-infection")
print("in one tick. Single-tick co-infection (n=91, single_tick_v3) is the direct co-infection measure.")"""

NEW_CAVEAT = """print("\\nCAVEAT: co-occurrence in a POOL = both detected in the same pool, NOT proof of co-infection")
print("in one tick. Single-tick co-infection is the direct measure: n = 93, from")
print("results/single_tick_v3/single_tick_v3_full.tsv, computed by")
print("scripts/recompute_cooccurrence.py. NOT n=91; the cohort has been 93 since the")
print("2026-08-07 exclusion rule (ground truth 1.6).")
print("Bartonella is excluded from BOTH analyses, here and there, for the same reason")
print("(ground truth 1.16). Neither is a pathogen co-occurrence result if it is in.")"""

REPLACEMENTS = [
    (OLD_LIST, NEW_LIST, "the pathogen list"),
    (OLD_CAVEAT, NEW_CAVEAT, "the closing caveat"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "co_detections" in source:
        sys.exit("already patched: " + TARGET + " defines co_detections. "
                 "Nothing changed.")

    for old, new, label in REPLACEMENTS:
        count = source.count(old)
        if count != 1:
            sys.exit("expected exactly one occurrence of " + label
                     + ", found " + str(count)
                     + ". Patch by hand rather than trusting this.")
        source = source.replace(old, new, 1)

    compile(source, TARGET, "exec")
    with open(TARGET, "w", encoding="utf-8") as handle:
        handle.write(source)
    print("patched " + TARGET)
    print()
    print("Expect: samples 558 -> 602, and 3 pairs rather than 6")
    print("(Bartonella's three pairs go; Ehrlichia was already all-zero).")
    print()
    print("run: python3 " + TARGET)


if __name__ == "__main__":
    main()
