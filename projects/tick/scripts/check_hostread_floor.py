#!/usr/bin/env python3
"""check_hostread_floor.py

I find where the minimum host-read count belongs, by showing the specimens
whose collection host is not their dominant blood-meal group against the
distribution of host reads for everything else. I write nothing.

WHY
---
Section 3.1 ends with a PENDING: four specimens had the collection host
outranked by something else, all of them with few host reads in total, so a
floor is applied before a blood meal is called, and nobody has said what the
floor is. It is not a number to pick by taste. If those four sit in the bottom
tail of the host-read distribution, the floor is wherever the tail ends, and
the data says where that is.

WHAT I REPORT
-------------
    the distribution of total host reads across deer specimens
    every specimen whose dominant group is not its collection host
    what each candidate floor would keep and what it would still let through

A floor that removes the four but also removes a quarter of the cohort is not
worth having, so the cost is printed beside the benefit.

I take the loading, grouping and cohort definitions from the existing
blood-meal scripts rather than reimplementing them, so the numbers here are
the same numbers Section 3.1 reports.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_hostread_floor.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

UNNAMED = "Unassigned or other"
FLOORS = (0, 100, 250, 500, 1000, 2000, 2500, 3000, 5000, 10000)


def deer_module():
    try:
        import check_deer_bloodmeal as module
    except ImportError as problem:
        sys.exit("cannot import check_deer_bloodmeal (%s). It must sit in "
                 "scripts/ beside me." % problem)
    return module


def main():
    deer_side = deer_module()
    base = deer_side.helpers()
    module = base.bloodmeal()
    table, finder = module.lookup()
    counts = module.read_table()
    path, covered, mapping = module.read_taxonomy(list(counts.columns))
    print("cytb table: %d libraries x %d ASVs" % counts.shape)
    print("taxonomy:   %s (%d of %d ASVs)" % (path, covered, counts.shape[1]))

    groups = pd.Series({f: module.group_of(mapping.get(f, ""))
                        for f in counts.columns})
    by_group = counts.T.groupby(groups).sum().T
    named = [c for c in by_group.columns if c != UNNAMED]
    animals = pd.Series([module.animal_of(table, finder, s)
                         for s in counts.index], index=counts.index)
    hosts = base.read_hosts()

    deer = deer_side.deer_specimens(module, table, finder)
    rows = []
    for animal in sorted(deer):
        libraries = animals.index[animals == animal]
        if not len(libraries):
            continue
        totals = by_group.loc[libraries, named].sum(axis=0)
        total = float(totals.sum())
        if total <= 0:
            continue
        ordered = totals.sort_values(ascending=False)
        expected = base.collected_group(hosts.get(animal, ""))
        rows.append({"specimen": animal,
                     "expected": expected or "(none recorded)",
                     "dominant": ordered.index[0],
                     "dominant_reads": int(ordered.iloc[0]),
                     "total": int(total),
                     "matches": bool(expected) and ordered.index[0] == expected})
    frame = pd.DataFrame(rows)
    if frame.empty:
        sys.exit("no deer specimen has cytb host reads. Something is wrong "
                 "upstream, not here.")
    print("deer specimens with cytb host reads: %d" % len(frame))
    print()

    print("=" * 72)
    print("TOTAL HOST READS PER SPECIMEN")
    print("=" * 72)
    for label, value in (("minimum", frame["total"].min()),
                         ("5th percentile", frame["total"].quantile(0.05)),
                         ("10th percentile", frame["total"].quantile(0.10)),
                         ("25th percentile", frame["total"].quantile(0.25)),
                         ("median", frame["total"].median()),
                         ("75th percentile", frame["total"].quantile(0.75)),
                         ("maximum", frame["total"].max())):
        print("  %-18s %12s" % (label, "{:,}".format(int(value))))
    print()

    wrong = frame[~frame["matches"]].sort_values("total")
    print("=" * 72)
    print("SPECIMENS WHOSE COLLECTION HOST IS NOT THE DOMINANT GROUP")
    print("=" * 72)
    if wrong.empty:
        print("  none.")
    else:
        print("  %-22s %-16s %-16s %10s %10s"
              % ("specimen", "collected from", "dominant", "its reads",
                 "total"))
        for _, row in wrong.iterrows():
            print("  %-22s %-16s %-16s %10s %10s"
                  % (str(row["specimen"])[:22], str(row["expected"])[:16],
                     str(row["dominant"])[:16],
                     "{:,}".format(row["dominant_reads"]),
                     "{:,}".format(row["total"])))
    print()
    print("  Section 3.1 says four. If this says four, the paragraph and this")
    print("  check agree and the floor can be set from the totals column.")
    print()

    print("=" * 72)
    print("WHAT EACH CANDIDATE FLOOR WOULD DO")
    print("=" * 72)
    print("  %8s %10s %9s %14s" % ("floor", "specimens", "of %d" % len(frame),
                                   "still wrong"))
    for floor in FLOORS:
        kept = frame[frame["total"] >= floor]
        remaining = int((~kept["matches"]).sum())
        print("  %8s %10d %8.0f%% %14d"
              % ("{:,}".format(floor), len(kept),
                 100.0 * len(kept) / len(frame), remaining))
    print()

    if not wrong.empty:
        clean = [f for f in FLOORS
                 if int((~frame[frame["total"] >= f]["matches"]).sum()) == 0]
        if clean:
            floor = min(clean)
            kept = frame[frame["total"] >= floor]
            print("  The lowest floor that removes every mismatch is %s "
                  "reads," % "{:,}".format(floor))
            print("  which keeps %d of %d specimens (%.0f%%)."
                  % (len(kept), len(frame), 100.0 * len(kept) / len(frame)))
        else:
            print("  No floor in the list removes every mismatch, which means")
            print("  at least one of them has plenty of host reads and is not")
            print("  a depth problem. That one needs looking at on its own.")
    print()


if __name__ == "__main__":
    main()
