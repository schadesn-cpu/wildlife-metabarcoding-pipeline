#!/usr/bin/env python3
"""patch_species_figure_to_lookup.py

I put build_species_figure.py on the shared library lookup. This CHANGES
LOCKED NUMBERS. Read the whole docstring before running it.

WHAT WAS WRONG
--------------
The script built its own lib2spec from name_18S by splitting on semicolons and
swapping a literal "-18s" suffix, then discarded every miss under the claim
"(controls and excluded specimens — expected)".

diagnose_species_figure_join.py, 2026-09-18:

    18S tick-genus libraries: 907
      resolved by lib2spec:               498
      resolved by the shared lookup:      853
      difference:                         355 libraries
      distinct specimens missed:          185
      of which flagged is_single_tick:     33

The 355 fall into the three name shapes TICK_GROUND_TRUTH section 3 already
documents, none of which lib2spec handles:

    bare accession        EQ02440612          (no marker tag to swap)
    replicate/subsample   EQ02441291-2 .. -10 (suffix after the accession)
    ALHT alias            ALHT-Deer-103       -> ALHT WTD_103

WHAT THIS CHANGES
-----------------
Everything downstream of the join: the informative-specimen count, the Ixodes
and Dermacentor counts, the per-host proportions, and the Dermacentor detail
table. Before this patch the figure reported 316 informative, 265 Ixodes, 20
Dermacentor, and zero Dermacentor from deer.

THE ZERO-DEER RESULT IS THE ONE TO WATCH. It is what justifies restricting the
diversity statistics to the deer subset, and the script asserts it at the end.
If it survives on the full join, it is a much stronger result than it was on
two-thirds of the data. If it does not, several downstream decisions rest on
it and it needs raising before anything is written.

Summing replicate libraries into one specimen is correct here and is the
script's own stated RULE 1: reads are summed across a specimen's libraries
before thresholding, because a specimen split into two tubes must not fail a
threshold a single tube would pass.

I also replace the hardcoded "expected" claim with a real breakdown, written
to results/analysis/species_figure_unmatched.tsv, so the next person does not
have to take anyone's word for what was dropped.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_species_figure_to_lookup.py
    python3 scripts/build_species_figure.py
"""
import os
import sys

TARGET = "scripts/build_species_figure.py"
NEEDS = "scripts/tick_lookup.py"

OLD_JOIN = '''lib2spec = {}
for _, r in cw.iterrows():
    for nm in str(r["name_18S"]).split(";"):
        nm = nm.strip()
        if not nm or nm.lower() == "nan":
            continue
        for form in {nm, nm[:-4] if nm.endswith("-18s") else nm + "-18s"}:
            lib2spec[form] = r["sample_id"]
say(f"library keys:     {len(lib2spec)}")

tg["lib"] = tg["sample_id"].astype(str).str.replace(r"\\..*$", "", regex=True)
tg["specimen"] = tg["lib"].map(lib2spec)

unmatched = tg[tg["specimen"].isna()]
say(f"\\nlibraries with no crosswalk row: {len(unmatched)} "
    f"(controls and excluded specimens — expected)")
'''

NEW_JOIN = '''# Name resolution is NOT done here. It is done once, by
# build_4marker_crosswalk.py, and read from results/metadata/library_lookup.tsv.
#
# What used to sit here built a private lib2spec from name_18S by splitting on
# semicolons and swapping a literal "-18s" suffix, then discarded every miss as
# "(controls and excluded specimens — expected)". Measured 2026-09-18: it
# resolved 498 of 907 libraries where the lookup resolves 853, missing 185
# distinct specimens, 33 of them single ticks. The misses were bare accessions
# (EQ02440612), replicate suffixes (EQ02441291-2 .. -10) and the ALHT alias
# (ALHT-Deer-103 -> ALHT WTD_103): the three shapes in TICK_GROUND_TRUTH
# section 3.
#
# Do not reintroduce a local resolver. If a library does not resolve, that is a
# crosswalk problem and it is now named below rather than assumed away.
sys.path.insert(0, "scripts")
import tick_lookup

LOOKUP = tick_lookup.load()

tg["lib"] = tg["sample_id"].astype(str).str.replace(r"\\..*$", "", regex=True)
tg["specimen"] = [
    tick_lookup.resolve(LOOKUP, "18S", lib)
    or tick_lookup.resolve(LOOKUP, "specimen", lib)
    for lib in tg["lib"]]

unmatched = tg[tg["specimen"].isna()]
say(f"\\nlibraries with no crosswalk specimen: {len(unmatched)} of {len(tg)}")

# No claim about what these are. The list is written out so it can be checked.
_ctrl_markers = ("NTC", "NEG", "BLK", "POS", "XB", "EXTRACTIONBLANK",
                 "PC-WELL", "BLANK", "CONTROL")
if len(unmatched):
    _u = unmatched[["sample_id", "lib"]].copy()
    _u["looks_like_control"] = _u["lib"].astype(str).str.upper().apply(
        lambda s: any(m in s for m in _ctrl_markers))
    _u.to_csv(f"{ODIR}/species_figure_unmatched.tsv", sep="\\t", index=False)
    _n_ctrl = int(_u["looks_like_control"].sum())
    say(f"  {_n_ctrl} carry a control prefix; {len(_u) - _n_ctrl} do not")
    say(f"  named in {ODIR}/species_figure_unmatched.tsv")
    if len(_u) - _n_ctrl:
        say("  The non-control ones are libraries with no crosswalk row. That is")
        say("  a crosswalk gap, not an expected loss. Check them before quoting")
        say("  any number from this figure.")
'''

OLD_DUP = '''say("  - categories OVERLAP: these are pooled samples and a pool off one animal")
say("    can contain both genera, so a specimen may be positive for both")
'''
NEW_DUP = ''

OLD_EXPECT = '''WHAT TO EXPECT (locked 2026-08-07)
----------------------------------
18 Dermacentor specimens: 9 moose (NH 5, VT 4), 5 Virginia opossum (ME 4, NH 1),
4 porcupine (ME), and ZERO deer. The zero-deer result is what justifies restricting
the diversity statistics to the deer subset: it controls host AND tick species at once.'''

NEW_EXPECT = '''WHAT TO EXPECT
--------------
The counts are NOT restated here. A hardcoded expectation in a docstring goes
stale silently, which is what happened to this one: it said 18 Dermacentor
specimens with 9 moose while the script had been printing 20 and 11 for weeks,
and nothing noticed because prose cannot fail. The same failure is recorded for
GOLDEN in build_single_tick_v2.py.

The current counts are in results/analysis/tick_species_report.txt, written by
this script on every run.

The result that MATTERS, and that this script does assert rather than describe,
is ZERO Dermacentor from deer. That is what justifies restricting the diversity
statistics to the deer subset: it controls host and tick species at once. The
script prints a WARNING instead if it ever stops being true.

2026-09-18: the library join moved from a private resolver to the shared
library lookup, which added 185 specimens the old join was silently dropping.
Any count quoted from a run before that date is on two-thirds of the data.'''

REPLACEMENTS = [
    (OLD_EXPECT, NEW_EXPECT, "the WHAT TO EXPECT docstring block"),
    (OLD_JOIN, NEW_JOIN, "the lib2spec join block"),
    (OLD_DUP, NEW_DUP, "the duplicated overlap caption lines"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    if not os.path.isfile(NEEDS):
        sys.exit("not found: " + NEEDS + "\nWrite tick_lookup.py first.")
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
                     + ". Patch by hand rather than trusting this.")
        source = source.replace(old, new, 1)

    compile(source, TARGET, "exec")
    with open(TARGET, "w", encoding="utf-8") as handle:
        handle.write(source)
    print("patched " + TARGET)
    print()
    print("THIS WILL MOVE LOCKED NUMBERS. Before this patch the figure said:")
    print("    informative specimens   316")
    print("    Ixodes-positive         265   (83.9%)")
    print("    Dermacentor-positive     20   (6.3%)")
    print("    Dermacentor from deer     0")
    print("    moose                  11/11 Dermacentor")
    print()
    print("Write the new values down next to those, and watch the zero-deer")
    print("line: the script WARNs instead of asserting if it stops holding.")
    print()
    print("now run: python3 " + TARGET)


if __name__ == "__main__":
    main()
