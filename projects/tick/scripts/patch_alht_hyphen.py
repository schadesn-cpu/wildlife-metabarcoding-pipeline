#!/usr/bin/env python3
"""patch_alht_hyphen.py

I close the missing-hyphen name shape in build_4marker_crosswalk.py. Two
changes, both in the crosswalk, because that is where name resolution lives.

THE FINDING
-----------
Rows 262 and 263 of 'All Run Data - w Metadata' are the SAME specimen:

    Num 262  Barcode USGSTICKS-LHT-Deer-143  27 ticks  ALHT-Deer-143-1-16s
    Num 263  Barcode USGSTICKS-LHT-Deer-143  27 ticks  ALHT-Deer143-2-16s

Identical in every metadata field: barcode, accession date 2025-02-18,
collection date 2024-11-17, host, town, tick count, and all four
run batches. The only difference is a missing hyphen in two of the four
library names.

The ALHT/LHT rule keys on ALHT_LIB_RE = ^ALHT-Deer-(\\d+), so it re-attributed
row 262 to ALHT WTD_143 and left row 263 on USGSTICKS-LHT-Deer-143. One
animal's libraries, split across two specimens with different hosts, different
towns and different tick counts (10 against 27).

WHY BOTH BELONG TO ALHT WTD_143
-------------------------------
This is the Num 247 case exactly, and Lawrence Gordon's dating argument
applies unchanged. Both rows carry 16S run batch AMP-120924-USGS-TICKS-P3-5,
sequenced 2024-12-09, on rows accessioned 2025-02-18. The libraries predate
the accession whose barcode was backfilled over them, which is the documented
basis of the rule. USGSTICKS-LHT-Deer-143 survives as its own specimen with
its own library, USGSTICKS-LHT-Deer-143-01, exactly as
USGSTICKS-LHT-Deer-091 does.

WHAT I CHANGE
-------------
1. ALHT_LIB_RE becomes ^ALHT-?Deer-?(\\d+), so the rule sees row 263. Measured
   against the workbook: exactly ONE row changes.

2. lookup_key_forms emits a hyphen-normalised variant, so alhtdeer-02 and
   alht-deer-02 are one key. 138 library names carry the missing-hyphen shape
   (ALHTDeer-02, ALHTDeer-103-1, ALHTDeer-91-4 and so on). Today they resolve
   only when the feature table happens to spell a name the same way the
   crosswalk does. ALHT-Deer-143-2 resolved at all only because that
   specimen's 18S name IS hyphenated and library names are registered into the
   'specimen' namespace too. That was luck.

I also add an assertion that ALHT WTD_143 ends up holding both libraries, so
this cannot silently regress.

AFTER RUNNING THIS
------------------
Rebuild the crosswalk, then rerun the three scripts already on the lookup and
confirm nothing else moves:

    python3 scripts/build_4marker_crosswalk.py
    python3 scripts/build_concordance.py        # expect n = 89, 16/15/3/55
    python3 scripts/build_qpcr_figure.py        # expect n = 89
    python3 scripts/build_species_figure.py     # expect 476 / 406 / 20, 0 of 424

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_alht_hyphen.py
"""
import os
import sys

TARGET = "scripts/build_4marker_crosswalk.py"

OLD_RE = "ALHT_LIB_RE = re.compile(r'^ALHT-Deer-(\\d+)')"
NEW_RE = ("# 2026-09-18: hyphen-tolerant. 138 library names in the workbook are\n"
          "# written ALHTDeer-NN or ALHT-DeerNN rather than ALHT-Deer-NN, and the\n"
          "# strict form missed row Num 263 (ALHT-Deer143-2-16s), leaving one\n"
          "# specimen's two libraries on two different specimens. Exactly one row\n"
          "# changes; verified against the workbook 2026-09-18.\n"
          "ALHT_LIB_RE = re.compile(r'^ALHT-?Deer-?(\\d+)')")

OLD_FORMS = """    expanded = set(forms)
    for form in forms:
        alias = LOOKUP_ALHT_ALIAS.match(form)
        if alias:
            expanded.add('alht-wtd-' + alias.group(1))
            expanded.add('alht-deer-' + alias.group(1))
    return expanded"""

NEW_FORMS = """    # Missing-hyphen shape first: alhtdeer-02 and alht-deer02 become
    # alht-deer-02, so all three spellings are one key. 138 library names carry
    # it. Then the WTD/Deer alias, which needs the hyphenated form to match.
    expanded = set(forms)
    for form in list(expanded):
        shape = LOOKUP_ALHT_SHAPE.match(form)
        if shape:
            expanded.add('alht-deer-' + shape.group(1))
    for form in list(expanded):
        alias = LOOKUP_ALHT_ALIAS.match(form)
        if alias:
            expanded.add('alht-wtd-' + alias.group(1))
            expanded.add('alht-deer-' + alias.group(1))
    return expanded"""

OLD_ALIAS_DEF = "LOOKUP_ALHT_ALIAS = re.compile(r'^alht-(?:wtd|deer)-(\\d+)$')"
NEW_ALIAS_DEF = ("LOOKUP_ALHT_ALIAS = re.compile(r'^alht-(?:wtd|deer)-(\\d+)$')\n"
                 "LOOKUP_ALHT_SHAPE = re.compile(r'^alht-?deer-?(\\d.*)$')")

OLD_CHECK = """for _bc in ('ALHT WTD_09','USGSTICKS-2213'):"""
NEW_CHECK = """# 2026-09-18: one specimen's libraries must not split across two specimens.
# Rows Num 262 and 263 are the same animal and differ only by a missing hyphen
# in two library names; before the hyphen-tolerant ALHT_LIB_RE, row 263 stayed
# on USGSTICKS-LHT-Deer-143 while row 262 moved to ALHT WTD_143.
_143 = x[x['sample_id'] == 'ALHT WTD_143']
if len(_143):
    _names = str(_143['name_16S'].iloc[0])
    check('143-1' in _names and '143-2' in _names,
          "ALHT WTD_143 holds BOTH 16S libraries (-143-1 and -143-2); got: "
          + _names)
else:
    say("  WARN ALHT WTD_143 is not in the crosswalk; skipping the "
        "split-specimen check")

for _bc in ('ALHT WTD_09','USGSTICKS-2213'):"""

REPLACEMENTS = [
    (OLD_RE, NEW_RE, "the ALHT_LIB_RE definition"),
    (OLD_ALIAS_DEF, NEW_ALIAS_DEF, "the lookup alias regex"),
    (OLD_FORMS, NEW_FORMS, "the key-form expansion"),
    (OLD_CHECK, NEW_CHECK, "the specimen assertions block"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "LOOKUP_ALHT_SHAPE" in source:
        sys.exit("already patched: " + TARGET + " defines LOOKUP_ALHT_SHAPE. "
                 "Nothing changed.")
    if "library_lookup.tsv" not in source:
        sys.exit(TARGET + " has no library lookup yet. Run "
                 "patch_add_library_lookup.py first.")

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
    print("Expect on the next crosswalk build:")
    print("  - the ALHT/LHT line reports 26 library rows across 16 specimens")
    print("    (was 25 across 16)")
    print("  - a new OK line: ALHT WTD_143 holds BOTH 16S libraries")
    print("  - the lookup grows by the hyphen-normalised keys")
    print("  - 563 specimens, 3,939 ticks and 93 single ticks UNCHANGED")
    print()
    print("run: python3 scripts/build_4marker_crosswalk.py")


if __name__ == "__main__":
    main()
