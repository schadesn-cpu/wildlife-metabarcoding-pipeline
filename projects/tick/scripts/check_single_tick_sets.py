#!/usr/bin/env python3
"""check_single_tick_sets.py

I say how many single ticks have 16S, how many have qPCR, and how many have
both, so the 93 and the 89 in two different figures can be named. I write
nothing.

WHY
---
The prevalence figure says n = 93, the concordance figure says n = 89, and the
qPCR figure says 75 deer out of 89 single ticks. If 93 is every single tick
with 16S and 89 is the ones that also have qPCR, then neither figure is wrong
and each only has to say which set it used. If the 89 is not exactly the
overlap, one of them is drawing on a cohort the paper never defines.

I print the four ticks that would separate the two sets by name, because if it
is four specific animals the answer is settled in one glance.

WHAT I READ
-----------
    results/metadata/crosswalk_4marker.tsv        is_single_tick, name_16S
    metadata/tick_full_metadata.xlsx, 'qPCR Results'
    results/decontam/16S/table_decontam_16S.tsv   which libraries survived

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_single_tick_sets.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
BOOK = "metadata/tick_full_metadata.xlsx"
SHEET = "qPCR Results"
COUNTS = "results/decontam/16S/table_decontam_16S.tsv"
TRUE_WORDS = ("1", "1.0", "true", "yes", "y")
SHOW = 10


def truthy(value):
    return str(value).strip().lower() in TRUE_WORDS


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s)." % problem)
    return tick_lookup.load(), tick_lookup


def resolve(table, module, value):
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    hit = module.resolve(table, "specimen", text)
    return "" if hit is None else str(hit)


def decontam_libraries():
    """I return the library names that survived decontam, or None."""
    if not os.path.isfile(COUNTS):
        print("no %s, so I cannot say which 16S libraries survived" % COUNTS)
        return None
    frame = pd.read_csv(COUNTS, sep="\t", index_col=0, nrows=2,
                        low_memory=False)
    names = set(str(c) for c in frame.columns)
    head = pd.read_csv(COUNTS, sep="\t", usecols=[0], low_memory=False)
    rows = set(str(v) for v in head[head.columns[0]])
    return names if len(names) > len(rows) else rows


def main():
    for path in (CROSSWALK, BOOK):
        if not os.path.isfile(path):
            sys.exit("missing %s. Run me from the project root." % path)
    table, module = lookup()

    crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    for column in ("sample_id", "barcode", "is_single_tick", "name_16S"):
        if column not in crosswalk.columns:
            sys.exit("%s has no '%s' column" % (CROSSWALK, column))

    crosswalk["animal"] = [resolve(table, module, b) or resolve(table, module, s)
                           for b, s in zip(crosswalk["barcode"],
                                           crosswalk["sample_id"])]
    flagged = crosswalk[[truthy(v) for v in crosswalk["is_single_tick"]]]
    singles = set(a for a in flagged["animal"] if a)
    print("single ticks by the crosswalk flag:            %d" % len(singles))

    named = set(a for a, n in zip(flagged["animal"], flagged["name_16S"])
                if a and pd.notna(n) and str(n).strip()
                and str(n).strip().lower() != "nan")
    print("  of those, with a 16S library named:          %d" % len(named))

    survived = decontam_libraries()
    if survived is not None:
        kept = set()
        for animal, name in zip(flagged["animal"], flagged["name_16S"]):
            if animal and str(name).strip() in survived:
                kept.add(animal)
        print("  of those, whose library survived decontam:   %d" % len(kept))
        sixteen = kept or named
    else:
        sixteen = named

    book = pd.read_excel(BOOK, sheet_name=SHEET, dtype=str)
    columns = {str(c).strip().lower(): c for c in book.columns}
    if "subsample_count" not in columns:
        sys.exit("the '%s' sheet has no subsample_count column; it has %s"
                 % (SHEET, ", ".join(map(str, book.columns))))
    tested = set()
    for _, row in book.iterrows():
        subs = str(row[columns["subsample_count"]]).strip()
        if subs not in ("1", "1.0"):
            continue
        for key in ("barcode", "host_id"):
            if key in columns:
                animal = resolve(table, module, row[columns[key]])
                if animal:
                    tested.add(animal)
                    break
    qpcr = tested & singles
    print("single ticks with a qPCR result:               %d" % len(qpcr))
    print()

    both = sixteen & qpcr
    print("=" * 64)
    print("THE THREE SETS")
    print("=" * 64)
    print("  16S only:    %4d" % len(sixteen - qpcr))
    print("  qPCR only:   %4d" % len(qpcr - sixteen))
    print("  both:        %4d" % len(both))
    print("  either:      %4d" % len(sixteen | qpcr))
    print()

    print("  Is 89 the overlap? %s" % ("yes" if len(both) == 89 else
                                       "no, the overlap is %d" % len(both)))
    print("  Is 93 the 16S set? %s" % ("yes" if len(sixteen) == 93 else
                                       "no, the 16S set is %d" % len(sixteen)))
    print()

    for label, group in (("in 16S but not qPCR", sixteen - qpcr),
                         ("in qPCR but not 16S", qpcr - sixteen)):
        if group:
            print("  %s (%d): %s" % (label, len(group),
                                     ", ".join(sorted(group)[:SHOW])))
    print()
    print("  If the overlap is 89 and the 16S set is 93, both figures are")
    print("  right and each needs its denominator named in the caption.")
    print()


if __name__ == "__main__":
    main()
