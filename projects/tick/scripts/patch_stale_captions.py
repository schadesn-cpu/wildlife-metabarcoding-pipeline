#!/usr/bin/env python3
"""patch_stale_captions.py

I fix three caption lines that still quote superseded numbers, in the two
figure scripts that print them. Small, but this is the exact class of defect
this project keeps paying for: a number hardcoded in prose, in a caption that
travels with the figure, which cannot fail when it goes stale.

WHAT IS STALE
-------------
build_species_figure.py
    "moose stayed 11/11 and deer 0-1/222-273"
    That was the 2026-08-10 check, run before the script was put on the shared
    lookup. The deer denominator is now 424. The claim still holds; it was
    rerun on the full join on 2026-09-18 by check_threshold_independence.py,
    which also named the single detection. The caption should say the numbers
    it was actually verified at.

build_qpcr_figure.py
    "Powassan is 1 of 100"
    The cohort is 89. This is a literal in a caption about why Wilson
    intervals are used at small n, and it names the wrong n.

    "a DIFFERENT population from the metabarcoding figure (93 sequenced single
     ticks, overlapping in 89)"
    This is no longer true and is now misleading rather than conservative.
    Since the cohort selects on is_single_tick, these 89 ARE a subset of the
    93: the same 89 as the concordance cohort. Ground truth 1.15.

WHAT I DO ABOUT IT
------------------
Where a number can be computed at run time, I compute it, so the caption
cannot go stale again. The deer denominator and the tick count both come from
data already in scope. The one thing I leave as a literal is the ALHT WTD_64
detail, because it is a specific observation with a date, not a derived count.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_stale_captions.py
    python3 scripts/build_species_figure.py
    python3 scripts/build_qpcr_figure.py
"""
import os
import sys

SPECIES = "scripts/build_species_figure.py"
QPCR = "scripts/build_qpcr_figure.py"

OLD_THRESH = '''say("  - the result is threshold-independent: across five rules (>=1, >=10, >=100")
say("    reads; >=1%, >=10% relative) moose stayed 11/11 and deer 0-1/222-273")'''

NEW_THRESH = '''_n_moose = int((df["host_bin"] == "Moose").sum())
_n_deer = int((df["host_bin"] == "Deer").sum())
say("  - the result is threshold-independent, rerun 2026-09-18 on the full")
say("    library join: across five rules (>=1, >=10, >=100 reads; >=1%, >=10%")
say(f"    relative) moose stayed {_n_moose}/{_n_moose} and deer 0-1 of {_n_deer}.")
say("    The single deer detection is ALHT WTD_64, 15 Dermacentor reads against")
say("    53,119 Ixodes (0.03%), inside the blank ceiling of 28 reads.")
say("    Source: results/analysis/threshold_independence.tsv")'''

OLD_WILSON = '''say("  - Wilson 95% intervals; correct for proportions at small n (Powassan is")
say("    1 of 100)")'''

NEW_WILSON = '''say(f"  - Wilson 95% intervals; correct for proportions at small n (Powassan")
say(f"    is 1 of {n_ticks})")'''

OLD_POP = '''say("  - a DIFFERENT population from the metabarcoding figure (93 sequenced single")
say("    ticks, overlapping in 89). The two panel-B co-infection figures are not")
say("    the same ticks and are not directly comparable")'''

NEW_POP = '''say("  - these ticks are a SUBSET of the sequenced single ticks, not a")
say("    different population. The cohort selects on is_single_tick, the same")
say("    rule the sequencing cohort uses, so these are exactly the same ticks")
say("    as the qPCR/metabarcoding concordance cohort. The sequenced specimens")
say("    absent here have only pooled qPCR rows. See ground truth 1.15.")
say("  - the two panel-B co-infection figures ARE the same ticks, but they")
say("    count different things; see the next point.")'''


def patch(path, replacements, sentinel):
    if not os.path.isfile(path):
        sys.exit("not found: " + path + "\nRun this from ~/meedlab/tick_project.")
    with open(path, encoding="utf-8") as handle:
        source = handle.read()
    if sentinel in source:
        print("already patched: " + path)
        return False
    for old, new, label in replacements:
        count = source.count(old)
        if count != 1:
            sys.exit("in " + path + ": expected exactly one occurrence of "
                     + label + ", found " + str(count)
                     + ". Patch by hand rather than trusting this.")
        source = source.replace(old, new, 1)
    compile(source, path, "exec")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(source)
    print("patched " + path)
    return True


def main():
    changed = patch(SPECIES,
                    [(OLD_THRESH, NEW_THRESH, "the threshold-independence caption")],
                    "threshold_independence.tsv")
    changed |= patch(QPCR,
                     [(OLD_WILSON, NEW_WILSON, "the Wilson caption"),
                      (OLD_POP, NEW_POP, "the different-population caption")],
                     "SUBSET of the sequenced single ticks")
    if changed:
        print()
        print("rerun both; no plotted value should move, only caption text:")
        print("  python3 scripts/build_species_figure.py")
        print("  python3 scripts/build_qpcr_figure.py")


if __name__ == "__main__":
    main()
