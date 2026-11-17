#!/usr/bin/env python3
"""check_rarefaction_floor.py

I read the alpha rarefaction visualisation and say whether the 1,000-read floor
used throughout the paper is where diversity saturates. I write nothing.

WHY
---
The Methods state that a library needs 1,000 reads to enter the analysis. At
the moment that is a convention. The rarefaction curves already computed can
turn it into a statement about these data: if the median diversity is within a
few percent of its plateau by 1,000 reads, the floor is justified, and if it is
still climbing steeply there, the floor is too low and the Methods should say
what that costs.

A .qzv is a zip. The alpha rarefaction one holds one CSV per metric, with a
column per depth and iteration, so I read the values directly rather than
asking for the curves to be rebuilt.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_rarefaction_floor.py
    python3 scripts/check_rarefaction_floor.py path/to/other.qzv
"""
import os
import re
import shutil
import sys
import tempfile
import zipfile

import pandas as pd

DEFAULT = "results/16S/rarefaction/alpha_rarefaction_16S.qzv"
MARKS = (0.90, 0.95, 0.99)
FLOORS = (500, 1000, 2000, 5000)
COLUMN = re.compile(r"depth-(\d+)_iter-(\d+)")


def read_tables(path):
    """I return {metric name: frame} for every CSV inside the visualisation."""
    if not os.path.isfile(path):
        sys.exit("missing %s. Run me from the project root." % path)
    workdir = tempfile.mkdtemp(prefix="qzv_")
    tables = {}
    try:
        with zipfile.ZipFile(path) as archive:
            wanted = [n for n in archive.namelist()
                      if n.endswith(".csv") and "/data/" in n]
            if not wanted:
                sys.exit("no CSV inside %s. Is it an alpha rarefaction "
                         "visualisation?" % path)
            for name in wanted:
                archive.extract(name, workdir)
                frame = pd.read_csv(os.path.join(workdir, name), index_col=0)
                if any(COLUMN.match(str(c)) for c in frame.columns):
                    metric = os.path.splitext(os.path.basename(name))[0]
                    tables[metric] = frame
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    if not tables:
        sys.exit("no rarefaction columns found. The CSVs inside %s do not "
                 "look like depth-N_iter-M tables." % path)
    return tables


def by_depth(frame):
    """I return {depth: (median across samples, samples retained)}."""
    groups = {}
    for column in frame.columns:
        hit = COLUMN.match(str(column))
        if hit:
            groups.setdefault(int(hit.group(1)), []).append(column)
    result = {}
    for depth, columns in sorted(groups.items()):
        block = frame[columns].apply(pd.to_numeric, errors="coerce")
        per_sample = block.mean(axis=1)
        kept = per_sample.notna().sum()
        if kept:
            result[depth] = (float(per_sample.median()), int(kept))
    return result


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    tables = read_tables(path)
    print("read %s" % path)
    print("metrics: %s" % ", ".join(sorted(tables)))
    print()

    for metric in sorted(tables):
        curve = by_depth(tables[metric])
        if not curve:
            continue
        depths = sorted(curve)
        best = max(curve[d][0] for d in depths)
        print("=" * 70)
        print(metric.upper())
        print("=" * 70)
        print("  %10s %10s %14s %10s" % ("depth", "samples", "median", "of max"))
        for depth in depths:
            median, kept = curve[depth]
            print("  %10d %10d %14.2f %9.1f%%"
                  % (depth, kept, median, 100 * median / best if best else 0))
        print()
        for mark in MARKS:
            reached = [d for d in depths if curve[d][0] >= mark * best]
            if reached:
                print("  median reaches %.0f%% of its plateau by %d reads"
                      % (100 * mark, min(reached)))
            else:
                print("  median never reaches %.0f%% of its plateau"
                      % (100 * mark))
        print()
        print("  samples retained at each candidate floor:")
        for floor in FLOORS:
            at_or_below = [d for d in depths if d <= floor]
            if not at_or_below:
                continue
            depth = max(at_or_below)
            median, kept = curve[depth]
            print("    %6d reads: %d samples, median %.2f (%.0f%% of plateau)"
                  % (floor, kept, median, 100 * median / best if best else 0))
        print()

    print("  A floor is defensible when the median is near its plateau there.")
    print("  If diversity is still climbing steeply at 1,000 reads, the floor")
    print("  is a practical compromise rather than a saturation point, and the")
    print("  Methods should say so.")
    print()


if __name__ == "__main__":
    main()
