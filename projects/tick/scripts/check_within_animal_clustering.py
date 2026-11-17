#!/usr/bin/env python3
"""check_within_animal_clustering.py

I ask whether two ticks off the same animal agree on infection status more
often than two ticks off different animals. That is what decides between the
two readings of the flat pool-size curve. I write nothing and change nothing.

WHY THIS EXISTS
---------------
check_burden_empirical_null.py, 2026-09-18, on 311 single-library animals:

    Borreliella    1 tick 38.9% (28/72)    2 ticks 43.9% (29/66)
                   27-62 ticks 36.0% (9/25)
    Bartonella     1 tick 11.1% (8/72)     27-62 ticks 12.0% (3/25)
    Anaplasma      1 tick 26.4% (19/72)    27-62 ticks 44.0% (11/25)

A pool of forty ticks detects at the same rate as one tick. The binomial pool
model predicts 87% at that top band for Borreliella and it is 36%.

Two things explain that and they lead opposite ways:

  BIOLOGY.  Infection is a property of the ANIMAL. Anaplasma and Borreliella
            are acquired by ticks feeding on an infected host, so ticks off one
            deer share its exposure and are near-perfectly correlated. A pool
            of forty then tells you what one tick tells you: what the deer had.
            That is a real finding.

  METADATA. NumTicks does not describe what went into the tube, and the
            pooling did not happen the way the run sheet says. Then the curve
            is flat because every library is really about one tick.

Both kill PoolTestR for this data, because its estimator assumes ticks within a
pool are independent, and the flat curve says pool size carries no information
either way. Only the first is publishable.

HOW TO TELL THEM APART
----------------------
Find animals from which more than one tick was measured SEPARATELY, and ask
whether those ticks agree. Under the biology reading they agree far more than
chance. Under the metadata reading they agree at the rate two random ticks do.

I look for per-tick tables under results/, resolve each tick to its animal
through the library lookup, and report what is available before testing
anything, because the answer may be that the data cannot support this and I
would rather say so than produce a number from four animals.

THE TEST
--------
Among animals contributing two or more ticks, I take the share of within-animal
tick PAIRS that agree, and compare it against a null in which the same tick
statuses are shuffled across animals. That holds the overall prevalence and the
number of ticks per animal fixed, so the only thing that can move the statistic
is whether ticks cluster by animal.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_within_animal_clustering.py
"""
import glob
import itertools
import os
import re
import sys

import numpy as np
import pandas as pd

GENERA = ["Anaplasma", "Borreliella", "Borrelia", "Ehrlichia", "Bartonella",
          "Babesia", "Powassan", "microti", "miyamotoi", "phagocytophilum",
          "burgdorferi", "mayonii", "muris"]
SEARCH = ["results/single_tick*/*.tsv", "results/single_tick*.tsv",
          "results/analysis/single_tick*.tsv", "results/qpcr/*.tsv",
          "results/analysis/*qpcr*.tsv", "results/*/*single*tick*.tsv"]
ID_HINTS = ["tick", "library", "sample", "barcode", "specimen", "name", "id"]
DRAWS = 4000
SEED = 20260918
MIN_PAIRED_ANIMALS = 8

RUN_SUFFIX = re.compile(r"\..*$")


def blank(value):
    """I return '' for anything meaning 'no value', and the text otherwise."""
    text = str(value).strip()
    return "" if text.lower() in ("nan", "none", "<na>", "") else text


def binary_columns(frame):
    """I return columns that name an organism and hold only 0/1-like values."""
    found = []
    for column in frame.columns:
        if not any(genus.lower() in str(column).lower() for genus in GENERA):
            continue
        values = pd.to_numeric(frame[column], errors="coerce").dropna()
        if not len(values):
            continue
        if set(np.unique(values)).issubset({0.0, 1.0}) and values.sum() > 0:
            found.append(column)
    return found


def load_lookup():
    """I return the lookup table and resolver, or (None, None)."""
    try:
        import tick_lookup
    except ImportError as problem:
        print("could not import tick_lookup (%s)." % problem)
        return None, None
    if not hasattr(tick_lookup, "resolve") or not hasattr(tick_lookup, "load"):
        print("tick_lookup has no load()/resolve().")
        return None, None
    return tick_lookup.load(), tick_lookup


def resolve_column(frame, column, table, module):
    """I resolve one column of names to animals, and report the hit rate."""
    labels = []
    for value in frame[column]:
        stem = RUN_SUFFIX.sub("", str(value).strip())
        hit = None
        for marker, candidate in (("16S", value), ("16S", stem),
                                  ("specimen", stem), ("specimen", value)):
            hit = module.resolve(table, marker, candidate)
            if hit is not None:
                break
        labels.append(blank(hit))
    return labels


def survey(table, module):
    """I find per-tick tables and report how many animals have several ticks."""
    print("=" * 72)
    print("1. WHAT PER-TICK DATA EXISTS, AND DOES ANY ANIMAL APPEAR TWICE?")
    print("=" * 72)
    paths = []
    for pattern in SEARCH:
        paths.extend(glob.glob(pattern))
    paths = sorted(set(p for p in paths if os.path.isfile(p)))
    if not paths:
        print("  no candidate per-tick tables found under results/.")
        print()
        return []
    usable = []
    for path in paths:
        try:
            frame = pd.read_csv(path, sep="\t", dtype=str)
        except (ValueError, OSError, pd.errors.ParserError) as problem:
            print("  %-58s unreadable (%s)" % (path[:58], problem))
            continue
        organisms = binary_columns(frame)
        if not organisms:
            print("  %-58s %4d rows, no binary organism column"
                  % (path[:58], len(frame)))
            continue
        best = None
        for column in frame.columns:
            if not any(hint in str(column).lower() for hint in ID_HINTS):
                continue
            labels = resolve_column(frame, column, table, module)
            resolved = sum(1 for label in labels if label)
            if best is None or resolved > best[1]:
                best = (column, resolved, labels)
        if best is None or not best[1]:
            print("  %-58s %4d rows, no column resolves to an animal"
                  % (path[:58], len(frame)))
            continue
        column, resolved, labels = best
        work = frame.assign(_animal=labels)
        work = work[work["_animal"] != ""]
        counts = work.groupby("_animal").size()
        repeated = counts[counts > 1]
        print("  %-58s %4d rows" % (path[:58], len(frame)))
        print("      id column '%s' resolves %d of %d rows to %d animals"
              % (column, resolved, len(frame), work["_animal"].nunique()))
        print("      animals contributing more than one tick: %d,"
              " covering %d ticks" % (len(repeated), int(repeated.sum())))
        print("      organism columns: %s" % ", ".join(organisms[:8]))
        if len(repeated) >= MIN_PAIRED_ANIMALS:
            usable.append((path, work, organisms, repeated.index))
        print()
    if not usable:
        print("  No table has at least %d animals with more than one tick, so"
              % MIN_PAIRED_ANIMALS)
        print("  this question cannot be answered from the files on disk. What")
        print("  it needs is per-tick results where several ticks from one")
        print("  animal were measured separately. If the qPCR was run per tick")
        print("  and the per-tick results still exist outside these tables,")
        print("  that is the file to find.")
        print()
    return usable


def pair_indices(groups):
    """I return the two arrays of positions making up every within-animal pair.

    Built once. The permutation test evaluates the statistic thousands of
    times, and rebuilding the pair list inside that loop is the difference
    between seconds and minutes.
    """
    left, right = [], []
    for positions in groups:
        for first, second in itertools.combinations(positions, 2):
            left.append(first)
            right.append(second)
    return np.asarray(left, dtype=np.int64), np.asarray(right, dtype=np.int64)


def agreement(statuses, left, right):
    """I return the share of within-animal tick pairs that agree."""
    if not len(left):
        return float("nan")
    return float(np.mean(statuses[left] == statuses[right]))


def test_clustering(path, work, organisms, repeated):
    """I compare within-animal agreement against a shuffled null."""
    print("=" * 72)
    print("2. DO TICKS FROM ONE ANIMAL AGREE? %s" % path)
    print("=" * 72)
    subset = work[work["_animal"].isin(repeated)].reset_index(drop=True)
    groups = [np.asarray(positions) for positions
              in subset.groupby("_animal").indices.values()]
    left, right = pair_indices(groups)
    print("  %d ticks from %d animals, %d within-animal pairs."
          % (len(subset), len(groups), len(left)))
    if not len(left):
        print("  no pairs to compare.")
        print()
        return
    print()
    print("  %-26s %8s %10s %12s %s"
          % ("organism", "positive", "observed", "if unclustered", "verdict"))
    generator = np.random.default_rng(SEED)
    for organism in organisms:
        statuses = pd.to_numeric(subset[organism],
                                 errors="coerce").fillna(0).values.astype(int)
        if statuses.sum() in (0, len(statuses)):
            print("  %-26s all ticks the same; nothing to test" % organism[:26])
            continue
        observed = agreement(statuses, left, right)
        shuffled = statuses.copy()
        null = np.empty(DRAWS)
        for index in range(DRAWS):
            generator.shuffle(shuffled)
            null[index] = agreement(shuffled, left, right)
        share = float((null >= observed).mean())
        if share <= 0.05:
            verdict = "CLUSTERED by animal, p = %.4f" % share
        else:
            verdict = "no more than chance, p = %.3f" % share
        print("  %-26s %7d %9.1f%% %11.1f%% %s"
              % (organism[:26], int(statuses.sum()), 100 * observed,
                 100 * float(np.median(null)), verdict))
    print()


def main():
    if not os.path.isdir("results"):
        sys.exit("no results/ directory here. Run this from the project root.")
    table, module = load_lookup()
    if table is None:
        sys.exit("the library lookup is needed to tell which animal a tick "
                 "came from. Build it:\n"
                 "    python3 scripts/build_4marker_crosswalk.py")
    usable = survey(table, module)
    for path, work, organisms, repeated in usable:
        test_clustering(path, work, organisms, repeated)

    print("=" * 72)
    print("HOW TO READ THIS")
    print("=" * 72)
    print("  CLUSTERED means ticks off one animal agree far more than two")
    print("  random ticks do. Infection is then an animal-level property, the")
    print("  flat pool-size curve is biology rather than a metadata fault, and")
    print("  it is reportable as a finding. PoolTestR is still invalid here,")
    print("  because its per-tick estimator assumes the opposite.")
    print()
    print("  NO MORE THAN CHANCE means ticks off one animal are as independent")
    print("  as any others. The flat curve then has no biological explanation")
    print("  and points at the pooling records instead: NumTicks would not be")
    print("  describing what went into the tube.")
    print()
    print("  A null result here is weak if the number of paired animals is")
    print("  small. Check the pair count above before reading anything into a")
    print("  non-rejection.")


if __name__ == "__main__":
    main()
