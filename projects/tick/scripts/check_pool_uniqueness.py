#!/usr/bin/env python3
"""check_pool_uniqueness.py

I ask whether the 298 pooled qPCR rows are 298 different animals, or whether
any animal appears twice. I write nothing.

WHY
---
This is the one genuinely blocking item on the verification list, and it does
not need a person. Section 3.2, Table 2 and Figure 2 all report pooled
detection per animal. That only holds if one pooled row is one animal. If
subsample_count counts subsamples rather than tubes, an animal could appear in
two rows and every pooled percentage would be inflated.

The same question was settled for the single-tick cohort by resolving each row
to an animal through tick_lookup and counting the distinct ones, which gave 89
from 89. I do the same here, and I do it for both cohorts so the single-tick
answer serves as a control on the method.

WHAT WOULD BE BAD
-----------------
Any animal appearing in more than one row of a cohort. I print every one of
them by name, with how many rows it holds, rather than a count, because the
next question would immediately be which ones.

WHAT I READ
-----------
    metadata/tick_full_metadata.xlsx, the qPCR sheet
    results/metadata/crosswalk_4marker.tsv, through check_qpcr_by_host

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_pool_uniqueness.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

BOOK = "metadata/tick_full_metadata.xlsx"
SHOW = 20


def qpcr_module():
    """I return check_qpcr_by_host, whose cohort rules I reuse."""
    try:
        import check_qpcr_by_host as module
    except ImportError as problem:
        sys.exit("cannot import check_qpcr_by_host (%s). It must sit in "
                 "scripts/ beside me." % problem)
    for wanted in ("read_crosswalk", "specimen_of", "QPCR_SHEET"):
        if not hasattr(module, wanted):
            sys.exit("check_qpcr_by_host has no %s; tell me what it calls "
                     "that and I will use it instead." % wanted)
    return module


def truthy(value):
    return str(value).strip().lower() in ("1", "1.0", "true", "yes", "y")


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s)." % problem)
    return tick_lookup.load(), tick_lookup


def report(label, rows, expected):
    """I print one cohort's row count, animal count and any repeats."""
    print("=" * 66)
    print(label.upper())
    print("=" * 66)
    named = [a for a in rows if a]
    blank = len(rows) - len(named)
    counted = {}
    for animal in named:
        counted[animal] = counted.get(animal, 0) + 1
    repeats = {a: n for a, n in counted.items() if n > 1}
    print("  rows:                    %d" % len(rows))
    print("  rows resolving to an animal: %d" % len(named))
    if blank:
        print("  rows resolving to nothing:   %d" % blank)
    print("  distinct animals:        %d" % len(counted))
    if expected is not None:
        print("  the paper reports:       %d" % expected)
    print()
    if not repeats:
        print("  No animal appears twice. One row is one animal.")
    else:
        print("  %d animals appear in more than one row, covering %d rows:"
              % (len(repeats), sum(repeats.values())))
        for animal, n in sorted(repeats.items(), key=lambda kv: -kv[1])[:SHOW]:
            print("    %-30s %d rows" % (str(animal)[:30], n))
        print()
        print("  Every pooled percentage in Section 3.2, Table 2 and Figure 2")
        print("  is computed per row, so repeats mean those numbers are wrong.")
    print()


def main():
    module = qpcr_module()
    table, finder = lookup()
    if not os.path.isfile(BOOK):
        sys.exit("missing %s. Run me from the project root." % BOOK)

    counts, flags, hosts = module.read_crosswalk()
    book = pd.read_excel(BOOK, sheet_name=module.QPCR_SHEET, dtype=str)
    if "subsample_count" not in book.columns:
        sys.exit("the %r sheet has no subsample_count column; it has %s"
                 % (module.QPCR_SHEET, ", ".join(map(str, book.columns))))
    print("qPCR sheet: %d rows" % len(book))

    single, pooled, other = [], [], 0
    for _, row in book.iterrows():
        if str(row["subsample_count"]).strip() not in ("1", "1.0"):
            other += 1
            continue
        specimen = module.specimen_of(table, finder, row)
        ticks = pd.to_numeric(pd.Series([counts.get(specimen)]),
                              errors="coerce").iloc[0]
        if truthy(flags.get(specimen)):
            single.append(specimen)
        elif pd.notna(ticks) and ticks > 1:
            pooled.append(specimen)
    print("rows with subsample_count other than 1: %d" % other)
    print()

    report("single ticks, as a control on the method", single, 89)
    report("pools, the question", pooled, 298)

    print("  If the pool block says 298 rows and 298 distinct animals, the")
    print("  blocking item on Kellen's list is answered and comes off.")
    print()


if __name__ == "__main__":
    main()
