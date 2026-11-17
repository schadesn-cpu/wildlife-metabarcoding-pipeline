#!/usr/bin/env python3
"""check_burden_empirical_null.py

I redo the tick-burden question without assuming anything about how pool size
turns into a positive result. I measure that relationship instead, and build
everything on the measurement. I write nothing and change nothing.

WHY THE PREVIOUS VERSION WAS WRONG
----------------------------------
check_pool_structure_and_burden.py assumed a pool of n ticks is positive with
probability 1-(1-p)^n, fitted p, and used it to work out the odds ratio tick
burden alone would produce. It came back with a burden null of 4.24 against an
observed 1.97 and declared the association explained.

Section 3 of the same run refuted the model, on single-library animals where
the pool size is not in dispute:

    Anaplasma, per-row positivity by the animal's tick count
        2-3 ticks   28.6%  (28/98)
        4-6         53.5%  (23/43)
        7-10        58.8%  (20/34)
        11-20       50.0%  (13/26)
        21+         42.1%  (16/38)

It rises to about ten ticks and then FALLS. The fitted model predicts roughly
93% at forty ticks. So the null was built on a curve the data rejects, in
exactly the region that drives it, and the 4.24 was an artifact of my
assumption rather than a property of the design.

The likely mechanism is dilution: one positive tick in a pool of fifty-four
does not clear detection. If that is right it also means PoolTestR is biased on
the large pools, because its estimator assumes detection does not depend on
pool size.

WHAT I DO INSTEAD
-----------------
1. MEASURE the relationship. Positivity against tick count on the
   single-library animals, with Wilson intervals, printed beside what the
   binomial pool model would predict, so the disagreement is on the page
   rather than in my argument.

2. BUILD THE NULL FROM THAT MEASUREMENT. For each animal I take the observed
   positivity of single-library animals carrying a similar number of ticks,
   and draw the two organisms INDEPENDENTLY at those rates. Whatever shape the
   burden relationship really has, flat or rising or falling, the null
   inherits it. If burden induces association, this null shows it. If burden
   induces nothing, this null sits at 1 and the association has to be
   explained some other way.

3. REDO ALIQUOTS VERSUS DIVIDED the same way. A multi-library row is compared
   against the observed rate at n ticks (aliquots) and at n/L ticks (divided),
   with both rates read off the measured curve instead of a formula.

WHAT THIS STILL CANNOT DO
-------------------------
The rates are estimated from single-library animals and applied to
multi-library ones. If those two groups differ in some way other than tick
count, that difference rides along. And nothing here separates a batch effect
from an animal effect, because 530 of 551 animals sit entirely within one
batch by design.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_burden_empirical_null.py
"""
import itertools
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

TABLE = "results/analysis/analysis_table_16S.tsv"
PATHOGENS = ["Anaplasma", "Borreliella", "Borrelia", "Ehrlichia"]
CO_DETECTIONS = ["Bartonella"]
DRAWS = 2000
SEED = 20260918
MIN_POSITIVES = 30
MIN_PER_BAND = 25
SUSPECT_TOKEN = "241209_P3-5"

RUN_SUFFIX = re.compile(r"\..*$")
RUN_TOKEN = re.compile(r"\.(.+)$")


def blank(value):
    """I return '' for anything meaning 'no value', and the text otherwise."""
    text = str(value).strip()
    return "" if text.lower() in ("nan", "none", "<na>", "") else text


def odds_ratio(a, b, c, d):
    """I return the Haldane-corrected odds ratio so a zero cell does not kill it."""
    if 0 in (a, b, c, d):
        return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return (a * d) / (b * c)


def two_by_two(x, y):
    """I return the four cells of a 2x2 from two boolean arrays."""
    return (int(np.count_nonzero(x & y)), int(np.count_nonzero(x & ~y)),
            int(np.count_nonzero(~x & y)), int(np.count_nonzero(~x & ~y)))


def resolve_animals(frame):
    """I return a specimen label per row, and the scheme that produced it."""
    try:
        import tick_lookup
    except ImportError:
        tick_lookup = None
    if tick_lookup is not None and hasattr(tick_lookup, "resolve"):
        table = tick_lookup.load()
        labels = []
        for sample in frame["sample"]:
            stem = RUN_SUFFIX.sub("", str(sample).strip())
            hit = None
            for marker, value in (("16S", sample), ("16S", stem),
                                  ("specimen", stem)):
                hit = tick_lookup.resolve(table, marker, value)
                if hit is not None:
                    break
            labels.append(blank(hit))
        if sum(1 for label in labels if label):
            return labels, "library lookup"
    if "Barcode" in frame.columns:
        labels = [blank(v) for v in frame["Barcode"]]
        if sum(1 for label in labels if label):
            return labels, "Barcode column"
    sys.exit("cannot identify animals in " + TABLE)


def build_bands(ticks, target=MIN_PER_BAND):
    """I cut tick counts into bands with at least `target` animals each.

    I grow a band until it holds enough animals AND the next tick count is
    different, so a single very common count cannot be split across two bands
    with different rates, which would be an artefact of the cut rather than of
    the data.
    """
    ordered = np.sort(np.asarray(ticks, dtype=float))
    edges = []
    start = 0
    while start < len(ordered):
        stop = min(start + target, len(ordered))
        while stop < len(ordered) and ordered[stop] == ordered[stop - 1]:
            stop += 1
        if len(ordered) - stop < target:
            stop = len(ordered)
        edges.append(ordered[stop - 1])
        start = stop
    edges[-1] = np.inf
    return np.asarray(edges, dtype=float)


def band_of(values, edges):
    """I return the band index for each value."""
    return np.searchsorted(edges, np.asarray(values, dtype=float),
                           side="left")


def wilson(count, total):
    """I return the Wilson interval as percentages."""
    if total == 0:
        return (0.0, 0.0)
    low, high = stats.binomtest(int(count), int(total)).proportion_ci(
        method="wilson")
    return (100 * low, 100 * high)


def fit_binomial_p(ticks, positives):
    """I fit 1-(1-p)^n by grid MLE, only so I can show it failing."""
    ticks = np.asarray(ticks, dtype=float)
    positives = np.asarray(positives, dtype=bool)
    keep = np.isfinite(ticks) & (ticks > 0)
    ticks, positives = ticks[keep], positives[keep]
    if not len(ticks) or not positives.any():
        return None
    best, best_p = -np.inf, None
    for candidate in np.linspace(1e-5, 0.95, 4000):
        miss = (1.0 - candidate) ** ticks
        hit = np.clip(1.0 - miss, 1e-12, 1.0)
        total = float(np.sum(np.log(hit[positives]))
                      + np.sum(np.log(np.clip(miss[~positives], 1e-12, 1.0))))
        if total > best:
            best, best_p = total, candidate
    return best_p


def measure_curve(single, present, edges):
    """I print observed positivity by tick band beside the model's prediction."""
    print("=" * 72)
    print("1. WHAT POSITIVITY ACTUALLY DOES AS TICK COUNT RISES")
    print("=" * 72)
    print("  %d single-library animals, where the pool size is the animal's"
          % len(single))
    print("  tick count and nothing is in dispute.")
    print()
    bands = band_of(single["_ticks"], edges)
    rates = {}
    for pathogen in present:
        column = pathogen + "_pos"
        fitted = fit_binomial_p(single["_ticks"], single[column] == 1)
        print("  %s   (binomial-pool fit: p = %s per tick)"
              % (pathogen, "%.4f" % fitted if fitted else "n/a"))
        print("    %-14s %6s %-22s %-16s %s"
              % ("ticks", "n", "observed", "95% CI", "model says"))
        table = []
        for index in range(len(edges)):
            mask = bands == index
            total = int(np.count_nonzero(mask))
            if not total:
                table.append(0.0)
                continue
            values = single.loc[mask.nonzero()[0], column] == 1
            count = int(values.sum())
            rate = count / float(total)
            table.append(rate)
            low, high = wilson(count, total)
            span = single.loc[mask.nonzero()[0], "_ticks"]
            predicted = ("%.1f%%" % (100 * (1 - (1 - fitted) ** span.median()))
                         if fitted else "n/a")
            print("    %-14s %6d %-22s %-16s %s"
                  % ("%g-%g" % (span.min(), span.max()), total,
                     "%.1f%% (%d/%d)" % (100 * rate, count, total),
                     "%.1f-%.1f" % (low, high), predicted))
        rates[pathogen] = np.asarray(table)
        print()
    print("  If 'observed' stops rising while 'model says' keeps climbing, the")
    print("  binomial pool model does not describe this assay, and any null")
    print("  built on it is fiction. Dilution is the obvious cause: one")
    print("  positive tick in a pool of fifty does not clear detection.")
    print()
    print("  That also bears on PoolTestR, whose estimator assumes detection")
    print("  does not depend on pool size.")
    print()
    return rates


def aliquot_or_divided(work, present, edges, rates):
    """I compare multi-library rows against the MEASURED rate at n and at n/L."""
    print("=" * 72)
    print("2. ALIQUOTS OR DIVIDED, ON THE MEASURED CURVE")
    print("=" * 72)
    many = work[work["_libraries"] >= 2]
    if len(many) < 20:
        print("  too few multi-library rows to compare.")
        print()
        return
    print("  %d rows from %d multi-library animals. Under ALIQUOTS each row"
          % (len(many), many["_animal"].nunique()))
    print("  saw the whole pool, so it should behave like an animal with n")
    print("  ticks. Under DIVIDED it saw n/L of them.")
    print()
    sizes = many["_ticks"].values.astype(float)
    libraries = many["_libraries"].values.astype(float)
    as_aliquot_band = band_of(sizes, edges)
    as_divided_band = band_of(np.maximum(sizes / libraries, 1.0), edges)
    print("  %-14s %10s %12s %12s %s"
          % ("", "observed", "if aliquot", "if divided", "verdict"))
    for pathogen in present:
        curve = rates[pathogen]
        observed = float((many[pathogen + "_pos"] == 1).mean())
        aliquot = float(np.mean(curve[as_aliquot_band]))
        divided = float(np.mean(curve[as_divided_band]))
        if abs(aliquot - divided) < 0.02:
            verdict = "cannot tell, predictions too close"
        elif abs(observed - aliquot) < abs(observed - divided):
            verdict = "closer to ALIQUOT"
        else:
            verdict = "closer to DIVIDED"
        print("  %-14s %9.1f%% %11.1f%% %11.1f%% %s"
              % (pathogen, 100 * observed, 100 * aliquot, 100 * divided,
                 verdict))
    print()
    print("  Where the two predictions sit close together the test has nothing")
    print("  to say, and that happens wherever the measured curve is flat.")
    print()


def null_report(work, present, edges, rates, label, mask=None):
    """I compare the observed association against a null built on the curve."""
    animals = work.groupby("_animal")
    names, sizes, positions = [], [], []
    for name, group in animals:
        names.append(name)
        sizes.append(float(group["_ticks"].max()))
        positions.append(np.asarray(group.index))
    sizes = np.asarray(sizes)
    keep = np.ones(len(names), dtype=bool)
    if mask is not None:
        keep = np.asarray([mask(name, work.loc[positions[i]])
                           for i, name in enumerate(names)])
    if int(np.count_nonzero(keep)) < 40:
        print("  %s: fewer than 40 animals, skipping." % label)
        print()
        return
    kept_positions = [positions[i] for i in range(len(names)) if keep[i]]
    kept_sizes = sizes[keep]
    bands = band_of(kept_sizes, edges)
    lengths = np.array([len(p) for p in kept_positions])
    padded = np.zeros((len(kept_positions), int(lengths.max())), dtype=np.int64)
    for i, p in enumerate(kept_positions):
        padded[i, :len(p)] = p
    rows = np.arange(len(kept_positions))
    positives = {p: (work[p + "_pos"].values == 1) for p in present}
    generator = np.random.default_rng(SEED)

    print("  %s: %d animals" % (label, len(kept_positions)))
    print("    %-26s %10s %11s %13s %s"
          % ("pair", "observed", "burden null", "null 97.5th", "verdict"))
    for first, second in itertools.combinations(present, 2):
        observed = []
        for _ in range(DRAWS // 4):
            choice = (generator.random(len(rows)) * lengths).astype(np.int64)
            picked = padded[rows, choice]
            observed.append(odds_ratio(*two_by_two(positives[first][picked],
                                                   positives[second][picked])))
        observed_median = float(np.median([v for v in observed
                                           if np.isfinite(v)]))
        p_first = rates[first][bands]
        p_second = rates[second][bands]
        null = []
        for _ in range(DRAWS):
            x = generator.random(len(bands)) < p_first
            y = generator.random(len(bands)) < p_second
            null.append(odds_ratio(*two_by_two(x, y)))
        null = np.asarray([v for v in null if np.isfinite(v)])
        share = float((null >= observed_median).mean())
        if share <= 0.05:
            verdict = "ABOVE burden, p = %.3f" % share
        elif share >= 0.95:
            verdict = "BELOW burden, p = %.3f" % share
        else:
            verdict = "within burden, p = %.3f" % share
        print("    %-26s %10.2f %11.2f %13.2f %s"
              % ("%s x %s" % (first, second), observed_median,
                 float(np.median(null)), float(np.percentile(null, 97.5)),
                 verdict))
    print()


def main():
    if not os.path.isfile(TABLE):
        sys.exit("not found: " + TABLE
                 + "\nRun this from the project root, after "
                   "python3 scripts/build_analysis_table.py")
    frame = pd.read_csv(TABLE, sep="\t")
    for column in ("sample", "NumTicks"):
        if column not in frame.columns:
            sys.exit(TABLE + " has no '" + column + "' column.")
    frame = frame.copy()
    tokens = []
    for sample in frame["sample"]:
        found = RUN_TOKEN.search(str(sample).strip())
        tokens.append(found.group(1) if found else "(no token)")
    frame["_token"] = tokens
    labels, scheme = resolve_animals(frame)
    frame["_animal"] = labels
    frame["_ticks"] = pd.to_numeric(frame["NumTicks"], errors="coerce")
    work = frame[(frame["_animal"] != "") & frame["_ticks"].notna()
                 & (frame["_ticks"] > 0)].reset_index(drop=True)
    work["_libraries"] = work.groupby("_animal")["sample"].transform("size")

    present = [p for p in PATHOGENS + CO_DETECTIONS
               if p + "_pos" in work.columns
               and int((work[p + "_pos"] == 1).sum()) >= MIN_POSITIVES]
    if len(present) < 2:
        sys.exit("fewer than two organisms with %d positives." % MIN_POSITIVES)

    single = work[work["_libraries"] == 1].reset_index(drop=True)
    if len(single) < 60:
        sys.exit("only %d single-library animals; too few to measure the "
                 "curve on." % len(single))
    edges = build_bands(single["_ticks"], MIN_PER_BAND)

    print("rows: %d" % len(work))
    print("animals: %d, identified by the %s"
          % (work["_animal"].nunique(), scheme))
    print("organisms: %s" % ", ".join(present))
    print("tick bands: %d, at least %d single-library animals each"
          % (len(edges), MIN_PER_BAND))
    print()

    rates = measure_curve(single, present, edges)
    aliquot_or_divided(work, present, edges, rates)

    print("=" * 72)
    print("3. OBSERVED ASSOCIATION AGAINST A NULL BUILT ON THE MEASURED CURVE")
    print("=" * 72)
    print("  Each animal's two organisms are drawn INDEPENDENTLY at the")
    print("  observed rate for its tick band. Whatever burden really does,")
    print("  including nothing, the null inherits it.")
    print()
    null_report(work, present, edges, rates, "all animals")
    null_report(work, present, edges, rates,
                "excluding " + SUSPECT_TOKEN,
                mask=lambda name, group: SUSPECT_TOKEN not in set(
                    group["_token"]))
    null_report(work, present, edges, rates,
                "only " + SUSPECT_TOKEN,
                mask=lambda name, group: SUSPECT_TOKEN in set(group["_token"]))

    print("=" * 72)
    print("HOW TO READ THIS")
    print("=" * 72)
    print("  If the null sits near 1, tick burden induces nothing in this")
    print("  dataset, the argument I have been making since this morning is")
    print("  empty, and an observed 1.97 is an association that needs another")
    print("  explanation.")
    print()
    print("  If the null sits well above 1, burden does induce association and")
    print("  the observed value has to clear it to mean anything.")
    print()
    print("  Comparing the three blocks says whether P3-5 is different once")
    print("  burden is accounted for on the measured curve rather than an")
    print("  assumed one.")
    print()
    print("  Ticks on one animal are not independent, so the true null is")
    print("  wider than this one. A result that only just clears the bar has")
    print("  not cleared much.")


if __name__ == "__main__":
    main()
