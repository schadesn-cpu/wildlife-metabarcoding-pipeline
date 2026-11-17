#!/usr/bin/env python3
"""check_single_tick_coinfection.py

I count how many of the individually screened ticks carried more than one
pathogen, and whether any pair of agents occurs together more often than chance.
I write nothing.

WHY
---
Co-infection is the question the Introduction opens on, and pooled screening
cannot answer it: a pool positive for two agents may hold two ticks with one
each. Only a tube holding one tick can say that one tick carried both. The
single-tick qPCR cohort is the one place in this study where the question is
answerable, so this is the whole of the co-infection result.

This is qPCR, not 16S. An earlier 16S count over 93 single ticks exists in
results/single_tick_v3/cooccurrence_n93.tsv and reports genus-level
co-occurrence. qPCR is the paper's primary pathogen assay, resolves Borreliella
to species, and covers the 89 animals whose cohort rule is already stated, so
the two should not be mixed.

WHAT I REPORT
-------------
    how many targets each animal carried, as a distribution
    the co-infection rate, over all animals and over positive animals only
    every combination actually observed, by name
    pairwise association for the targets with enough positives to test

WHAT WOULD BE MISLEADING
------------------------
Pairwise tests on a target with two positives. I test only pairs where both
agents reach MIN_POSITIVE animals, and I correct across the pairs I do test,
because twenty-one uncorrected pairs on 89 animals will manufacture a result.

An association here is not transmission biology. Two agents sharing a reservoir
host, or a tick feeding where both circulate, produce co-occurrence without any
interaction inside the tick.

WHAT I READ
-----------
    metadata/tick_full_metadata.xlsx and the crosswalk, through
        check_qpcr_by_host

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_single_tick_coinfection.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

MIN_POSITIVE = 5
PRETTY = {
    "a-phagocytophilum": "A. phagocytophilum",
    "b-burgdorferi": "B. burgdorferi",
    "b-microti": "B. microti",
    "b-miyamotoi": "B. miyamotoi",
    "powassan-virus": "Powassan virus",
    "b-mayonii": "B. mayonii",
    "e-muris": "E. muris",
}


def qpcr_module():
    """I return check_qpcr_by_host, whose cohort rules I reuse."""
    try:
        import check_qpcr_by_host as module
    except ImportError as problem:
        sys.exit("cannot import check_qpcr_by_host (%s). It must sit in "
                 "scripts/ beside me." % problem)
    for wanted in ("read_qpcr", "read_crosswalk", "specimen_of", "TARGETS",
                   "host_group"):
        if not hasattr(module, wanted):
            sys.exit("check_qpcr_by_host has no %s; tell me what it calls that "
                     "and I will use it instead." % wanted)
    return module


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s)." % problem)
    return tick_lookup.load(), tick_lookup


def fisher():
    """I return scipy's fisher_exact, or None if scipy is unavailable."""
    try:
        from scipy.stats import fisher_exact
    except ImportError:
        return None
    return fisher_exact


def benjamini_hochberg(pvalues):
    """I return BH-adjusted q values in the order the p values came in."""
    indexed = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    total = len(pvalues)
    adjusted = [0.0] * total
    running = 1.0
    for rank, position in enumerate(reversed(indexed), start=1):
        value = pvalues[position] * total / (total - rank + 1)
        running = min(running, value)
        adjusted[position] = min(1.0, running)
    return adjusted


def single_tick_calls(module, table, finder):
    """I return one row per single-tick animal, each target as 0 or 1."""
    counts, flags, hosts = module.read_crosswalk()
    book = module.read_qpcr()
    rows = []
    for _, row in book.iterrows():
        if str(row["subsample_count"]).strip() not in ("1", "1.0"):
            continue
        specimen = module.specimen_of(table, finder, row)
        if not specimen or not flags.get(specimen):
            continue
        entry = {"animal": specimen, "host": hosts.get(specimen, "")}
        for target in module.TARGETS:
            value = pd.to_numeric(pd.Series([row.get(target)]),
                                  errors="coerce").fillna(0).iloc[0]
            entry[target] = int(value > 0)
        rows.append(entry)
    frame = pd.DataFrame(rows)
    if frame.empty:
        sys.exit("no single-tick row resolved to a specimen. Check the sheet "
                 "name and the crosswalk.")
    before = len(frame)
    frame = frame.drop_duplicates(subset="animal", keep="first")
    if len(frame) < before:
        print("dropped %d repeated animal rows, keeping the first of each"
              % (before - len(frame)))
    return frame


def burden_table(frame, targets):
    """I print how many agents each animal carried."""
    burden = frame[targets].sum(axis=1)
    total = len(frame)
    print("=" * 70)
    print("HOW MANY AGENTS EACH TICK CARRIED")
    print("=" * 70)
    print("  %-18s %8s %9s" % ("agents", "ticks", "of %d" % total))
    for number in range(0, int(burden.max()) + 1):
        count = int((burden == number).sum())
        if count or number <= 2:
            print("  %-18d %8d %8.0f%%" % (number, count,
                                           100.0 * count / total))
    positive = int((burden >= 1).sum())
    multiple = int((burden >= 2).sum())
    print()
    print("  positive for at least one: %d of %d (%.0f%%)"
          % (positive, total, 100.0 * positive / total))
    print("  carrying two or more:      %d of %d (%.0f%%) of all ticks"
          % (multiple, total, 100.0 * multiple / total))
    if positive:
        print("  carrying two or more:      %d of %d (%.0f%%) of positive ticks"
              % (multiple, positive, 100.0 * multiple / positive))
    print()
    return burden


def combinations(frame, targets, burden):
    """I print every combination of two or more agents actually seen."""
    print("=" * 70)
    print("COMBINATIONS OBSERVED")
    print("=" * 70)
    seen = {}
    for _, row in frame[burden >= 2].iterrows():
        names = tuple(sorted(PRETTY.get(t, t) for t in targets if row[t]))
        seen[names] = seen.get(names, 0) + 1
    if not seen:
        print("  none. No tick carried more than one agent.")
        print()
        return
    for names, count in sorted(seen.items(), key=lambda kv: -kv[1]):
        print("  %-52s %d" % (" + ".join(names), count))
    print()


def pairwise(frame, targets, exact):
    """I print pairwise association for targets with enough positives."""
    print("=" * 70)
    print("PAIRWISE ASSOCIATION")
    print("=" * 70)
    usable = [t for t in targets if int(frame[t].sum()) >= MIN_POSITIVE]
    skipped = [t for t in targets if t not in usable]
    if skipped:
        print("  not tested, fewer than %d positives: %s"
              % (MIN_POSITIVE, ", ".join(PRETTY.get(t, t) for t in skipped)))
        print()
    if len(usable) < 2:
        print("  fewer than two targets have enough positives to pair.")
        print()
        return
    if exact is None:
        print("  scipy is unavailable, so I report counts without a test.")

    rows, pvalues = [], []
    for i, first in enumerate(usable):
        for second in usable[i + 1:]:
            both = int(((frame[first] == 1) & (frame[second] == 1)).sum())
            only_first = int(((frame[first] == 1) & (frame[second] == 0)).sum())
            only_second = int(((frame[first] == 0) & (frame[second] == 1)).sum())
            neither = int(((frame[first] == 0) & (frame[second] == 0)).sum())
            expected = (float(frame[first].sum()) * float(frame[second].sum())
                        / len(frame))
            if exact is None:
                odds, p = float("nan"), float("nan")
            else:
                odds, p = exact([[both, only_first], [only_second, neither]])
            rows.append({"pair": "%s x %s" % (PRETTY.get(first, first),
                                              PRETTY.get(second, second)),
                         "both": both, "expected": expected,
                         "only_first": only_first, "only_second": only_second,
                         "neither": neither, "odds": odds, "p": p})
            pvalues.append(p if p == p else 1.0)

    qvalues = benjamini_hochberg(pvalues)
    print("  %-34s %6s %9s %8s %9s %8s"
          % ("pair", "both", "expected", "odds", "p", "q"))
    for row, q in zip(rows, qvalues):
        print("  %-34s %6d %9.1f %8s %9s %8s"
              % (row["pair"], row["both"], row["expected"],
                 "%.2f" % row["odds"] if row["odds"] == row["odds"] else "-",
                 "%.3f" % row["p"] if row["p"] == row["p"] else "-",
                 "%.3f" % q))
    print()
    print("  'expected' is what independence would give: the two prevalences")
    print("  multiplied by the number of ticks. A pair near its expectation is")
    print("  co-occurring by chance, which is the answer the paper needs if no")
    print("  q falls below 0.05.")
    print()


def by_host(frame, targets, module, burden):
    """I print the burden split by collection host, deer being most of it."""
    print("=" * 70)
    print("BY COLLECTION HOST")
    print("=" * 70)
    groups = pd.Series([module.host_group(h) for h in frame["host"]],
                       index=frame.index)
    print("  %-26s %7s %10s %12s" % ("host", "ticks", "positive", "two or more"))
    for name in sorted(set(groups)):
        rows = groups.index[groups == name]
        part = burden.loc[rows]
        print("  %-26s %7d %10d %12d"
              % (str(name)[:26], len(rows), int((part >= 1).sum()),
                 int((part >= 2).sum())))
    print()


def main():
    module = qpcr_module()
    table, finder = lookup()
    frame = single_tick_calls(module, table, finder)
    targets = list(module.TARGETS)
    print("single-tick animals with a qPCR result: %d" % len(frame))
    print("targets: %s" % ", ".join(PRETTY.get(t, t) for t in targets))
    print()

    burden = burden_table(frame, targets)
    combinations(frame, targets, burden)
    pairwise(frame, targets, fisher())
    by_host(frame, targets, module, burden)

    print("=" * 70)
    print("READING IT")
    print("=" * 70)
    print("  These are single ticks, so a tick carrying two agents carried")
    print("  both itself. That cannot be said of any pooled result.")
    print()
    print("  Co-occurrence is not interaction. Two agents sharing a reservoir")
    print("  host will appear together without anything happening inside the")
    print("  tick, and this design cannot separate those.")
    print()


if __name__ == "__main__":
    main()
