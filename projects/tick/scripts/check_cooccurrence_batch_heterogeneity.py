#!/usr/bin/env python3
"""check_cooccurrence_batch_heterogeneity.py

I ask whether the pooled pathogen co-occurrence is the same association in
every sequencing batch, or one or two batches carrying the rest. I write
nothing and change nothing.

WHY THIS EXISTS
---------------
check_pooled_cooccurrence_confounding.py, 2026-09-18, on 602 pooled samples,
Anaplasma x Borreliella:

    crude                              OR 2.03
    Mantel-Haenszel over pool size     OR 2.05
    Mantel-Haenszel over read depth    OR 1.89
    Mantel-Haenszel over run batch     OR 1.79

The batch summary looks like survival. The strata behind it do not:

    AMP-021925-VDL-Tick-P7-16s            n=77   OR 0.98
    AMP-032725-USGS-ticks-batch-08-16s    n=81   OR 1.16
    AMP-032725-USGS-ticks-batch-09-16s    n=72   OR 1.78
    AMP-043024-DN-USGSTICKS-2023-2024     n=61   OR 0.99
    AMP-060925-USGS-Ticks-Batch-10-16s    n=41   OR 2.13
    AMP-060925-USGS-Ticks-Batch-11-16s    n=74   OR 2.22
    AMP-072924-DN-USGS-Tick-...Batch-02   n=43   OR 1.20
    AMP-111524-USGS-ticks-batch-06        n=14   OR 81.00   (one co-positive)
    AMP-120924-USGS-TICKS-P3-5            n=96   OR 5.86
    (no batch recorded)                   n=43   OR 8.23

Mantel-Haenszel assumes a COMMON odds ratio across strata and returns the
weighted average of them. When the strata disagree, that average summarises
nothing: 1.79 is a null in six batches averaged against a 5.86 in one. So the
earlier reading rule, "collapses toward 1 or holds near crude", was a false
dichotomy, and this script covers the case it left out.

WHAT I DO
---------
1. HOMOGENEITY. Breslow-Day with Tarone's correction, per pair, across
   batches. This is the standard test of "is one odds ratio enough to describe
   all the strata". A small p-value says no, and then the Mantel-Haenszel
   number should not be quoted at all.

2. LEAVE ONE BATCH OUT. For each batch, the crude and Mantel-Haenszel odds
   ratios recomputed with that batch dropped. If removing a single batch takes
   the association to 1, the association is that batch.

3. THE UNBATCHED 43. They carry the largest odds ratio on the board (8.23) and
   we do not know what they are. I print what the analysis table knows about
   them so they can be chased in the run sheet.

WHAT I DO NOT DO
----------------
I do not decide anything. The leave-one-out table is the part to read with your
own eyes.

WHAT THE TEST CAN AND CANNOT SEE, AT THESE SIZES
------------------------------------------------
Checked by simulation before this script was written, at the batch sizes above
(77, 81, 72, 61, 41, 74, 43, 14, 96):

    nine batches that TRULY share an odds ratio of 1.8 still produce observed
    per-batch ratios spanning about 0.7 to 4.5, half the time

which is close to the spread in the real output. So the scatter in that column
is not by itself evidence of anything, and my earlier reading of it as
heterogeneity was wrong.

Against a pattern like the one that column suggests (most batches near 1, one
near 6), both the Breslow-Day chi2 and the bootstrap below reject about two
times in three. So a rejection here means something and a non-rejection does
not settle it. Both are calibrated: under a true common odds ratio they reject
5% and 4% of the time respectively, against a nominal 5%.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_cooccurrence_batch_heterogeneity.py
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

TABLE = "results/analysis/analysis_table_16S.tsv"
PATHOGENS = ["Anaplasma", "Borreliella", "Borrelia", "Ehrlichia"]
BATCH_COLUMN = "RunBatch"
MIN_STRATUM = 10
SIM_DRAWS = 4000
SIM_SEED = 20260918


def blank(value):
    """I return '' for anything that means 'no value', and the text otherwise.

    Not .astype(str).replace(): pandas turns a replacement of the string 'nan'
    back into a real NaN. This project has been bitten by the string-versus-NA
    distinction twice.
    """
    text = str(value).strip()
    return "" if text.lower() in ("nan", "none", "<na>", "") else text


def odds_ratio(a, b, c, d):
    """I return the Haldane-corrected odds ratio so a zero cell does not kill it."""
    if 0 in (a, b, c, d):
        return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return (a * d) / (b * c)


def two_by_two(frame, first, second):
    """I return (both, first only, second only, neither) for one pair."""
    x = frame[first + "_pos"] == 1
    y = frame[second + "_pos"] == 1
    return (int((x & y).sum()), int((x & ~y).sum()),
            int((~x & y).sum()), int((~x & ~y).sum()))


def mantel_haenszel(tables):
    """I return the Mantel-Haenszel summary odds ratio over a list of 2x2s."""
    numerator = denominator = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            continue
        numerator += (a * d) / n
        denominator += (b * c) / n
    if denominator == 0:
        return None
    return numerator / denominator


def expected_under_common_or(a, b, c, d, psi):
    """I solve for the expected 'both' cell under a common odds ratio.

    Breslow-Day fixes the margins of each stratum and asks what the top-left
    cell would be if every stratum shared the odds ratio psi. That is a
    quadratic; I take the root inside the range the margins allow.
    """
    row_one = a + b
    col_one = a + c
    total = a + b + c + d
    if row_one == 0 or col_one == 0 or row_one == total or col_one == total:
        return None
    low = max(0.0, float(row_one + col_one - total))
    high = float(min(row_one, col_one))
    if abs(psi - 1.0) < 1e-9:
        return row_one * col_one / float(total)
    alpha = 1.0 - psi
    beta = total - row_one - col_one + psi * (row_one + col_one)
    gamma = -psi * row_one * col_one
    discriminant = beta * beta - 4.0 * alpha * gamma
    if discriminant < 0:
        return None
    root = discriminant ** 0.5
    for candidate in ((-beta + root) / (2.0 * alpha),
                      (-beta - root) / (2.0 * alpha)):
        if low - 1e-9 <= candidate <= high + 1e-9:
            return min(max(candidate, low), high)
    return None


def breslow_day(tables, psi):
    """I return (statistic, degrees of freedom, p) for Tarone's corrected test.

    The null is that one odds ratio describes every stratum. A small p means it
    does not, and the Mantel-Haenszel summary should not be quoted.
    """
    statistic = 0.0
    residual_sum = 0.0
    variance_sum = 0.0
    used = 0
    for a, b, c, d in tables:
        expected = expected_under_common_or(a, b, c, d, psi)
        if expected is None:
            continue
        row_one = a + b
        col_one = a + c
        total = a + b + c + d
        cells = (expected, row_one - expected, col_one - expected,
                 total - row_one - col_one + expected)
        if min(cells) <= 1e-9:
            continue
        variance = 1.0 / sum(1.0 / cell for cell in cells)
        statistic += (a - expected) ** 2 / variance
        residual_sum += a - expected
        variance_sum += variance
        used += 1
    if used < 2:
        return None
    if variance_sum > 0:
        statistic -= residual_sum ** 2 / variance_sum      # Tarone
    degrees = used - 1
    return statistic, degrees, float(stats.chi2.sf(statistic, degrees))


def simulate_common_or(tables, psi, draws=SIM_DRAWS, seed=SIM_SEED):
    """I simulate strata that truly share one odds ratio, at the observed sizes.

    WHY THIS IS HERE
    ----------------
    Reading a column of per-stratum odds ratios and calling the scatter
    "heterogeneity" does not work at these sample sizes. Nine strata of 14 to
    96 samples that genuinely share an odds ratio of 1.8 still produce observed
    ratios spanning roughly 0.7 to 4.5 half the time. So a spread has to be
    compared against the spread chance alone would give.

    For each stratum I hold the size and the number of samples positive for the
    first pathogen fixed, take that stratum's own baseline rate for the second
    pathogen from its unexposed samples, impose the common odds ratio psi, and
    redraw. I return the null distributions of the Breslow-Day statistic and of
    the max-to-min odds ratio span.

    The seed is fixed so this returns the same number twice. A p-value that
    moves between runs is not a p-value anyone can put in a manuscript.
    """
    generator = np.random.default_rng(seed)
    setup = []
    for a, b, c, d in tables:
        exposed = a + b
        unexposed = c + d
        if exposed == 0 or unexposed == 0:
            continue
        baseline = c / float(unexposed)
        baseline = min(max(baseline, 1e-6), 1 - 1e-6)
        odds = baseline / (1 - baseline) * psi
        setup.append((exposed, unexposed, odds / (1 + odds), baseline))
    if len(setup) < 2:
        return None
    statistics = []
    spans = []
    for _ in range(draws):
        drawn = []
        for exposed, unexposed, p_exposed, p_unexposed in setup:
            a = int(generator.binomial(exposed, p_exposed))
            c = int(generator.binomial(unexposed, p_unexposed))
            drawn.append((a, exposed - a, c, unexposed - c))
        summary = mantel_haenszel(drawn)
        if summary is None:
            continue
        test = breslow_day(drawn, summary)
        if test is not None:
            statistics.append(test[0])
        ratios = [odds_ratio(*table) for table in drawn
                  if table[1] > 0 and table[2] > 0]
        if len(ratios) >= 2:
            spans.append(max(ratios) / min(ratios))
    if not statistics or not spans:
        return None
    return np.array(statistics), np.array(spans)


def strata_for(frame, first, second, min_n=MIN_STRATUM):
    """I return [(batch, sub-frame)] for batches big enough to inform."""
    kept = []
    for level, sub in frame.groupby("_batch", observed=True):
        if level == "" or len(sub) < min_n:
            continue
        kept.append((str(level), sub))
    return kept


def report_pair(frame, first, second):
    """I print homogeneity and leave-one-out for one pathogen pair."""
    print("%s x %s" % (first, second))
    a, b, c, d = two_by_two(frame, first, second)
    crude = odds_ratio(a, b, c, d)
    print("  crude over all %d samples: both=%d  OR=%.2f" % (len(frame), a, crude))

    kept = strata_for(frame, first, second)
    if len(kept) < 2:
        print("  fewer than two batches large enough to compare. Nothing to test.")
        print()
        return
    tables = [two_by_two(sub, first, second) for _, sub in kept]
    summary = mantel_haenszel(tables)
    informative = sum(1 for _, bb, cc, _ in tables if bb > 0 and cc > 0)
    print("  %d batches of at least %d samples, covering %d samples"
          % (len(kept), MIN_STRATUM, sum(len(sub) for _, sub in kept)))
    print("  informative batches (both discordant cells non-zero): %d of %d"
          % (informative, len(tables)))
    if summary is None:
        print("  no Mantel-Haenszel summary is defined here.")
        print()
        return
    print("  Mantel-Haenszel summary OR: %.2f" % summary)

    ratios = [odds_ratio(*table) for table in tables
              if table[1] > 0 and table[2] > 0]
    observed_span = max(ratios) / min(ratios) if len(ratios) >= 2 else None

    test = breslow_day(tables, summary)
    print()
    print("  HOMOGENEITY (Breslow-Day with Tarone's correction)")
    if test is None:
        print("    too few usable batches to test. The summary rests on nothing.")
    else:
        statistic, degrees, pvalue = test
        print("    chi2 = %.2f on %d df, p = %.4g" % (statistic, degrees, pvalue))

    print()
    print("  IS THE SPREAD BIGGER THAN CHANCE?")
    print("    %d simulations of batches that truly share OR %.2f, at these"
          % (SIM_DRAWS, summary))
    print("    batch sizes and each batch's own baseline rate. Seed %d."
          % SIM_SEED)
    null = simulate_common_or(tables, summary)
    if null is None or observed_span is None:
        print("    not enough usable batches to simulate.")
    else:
        null_statistics, null_spans = null
        print("    observed per-batch OR span (max/min over informative"
              " batches): %.1f-fold" % observed_span)
        print("    under one common OR, that span is %.1f-fold at the median,"
              % float(np.median(null_spans)))
        print("    %.1f-fold at the 95th percentile"
              % float(np.percentile(null_spans, 95)))
        span_p = float((null_spans >= observed_span).mean())
        print("    -> a span this wide or wider arises by chance %.1f%% of the"
              " time" % (100 * span_p))
        if test is not None:
            empirical = float((null_statistics >= test[0]).mean())
            print("    Breslow-Day empirical p (not the chi2 approximation):"
                  " %.4g" % empirical)
            verdict = empirical
        else:
            verdict = span_p
        print()
        if verdict < 0.05:
            print("    The batches do NOT share one odds ratio. The")
            print("    Mantel-Haenszel summary above averages strata that")
            print("    disagree, and must not be quoted as an adjusted")
            print("    association.")
        else:
            print("    The spread across batches is within what chance gives at")
            print("    these sizes. That is NOT evidence the batches agree; the")
            print("    test only has about 2 chances in 3 of catching a real")
            print("    difference here. It means the scatter alone cannot be")
            print("    read as a batch effect, and the leave-one-out table")
            print("    below is the thing to look at instead.")

    print()
    print("  LEAVE ONE BATCH OUT")
    print("    %-46s %7s %7s %7s" % ("batch removed", "n left", "crude", "MH"))
    print("    %-46s %7d %7.2f %7.2f"
          % ("(none removed)", len(frame), crude, summary))
    rows = []
    for name, sub in kept:
        rest = frame[frame["_batch"] != name]
        ra, rb, rc, rd = two_by_two(rest, first, second)
        rest_tables = [two_by_two(other, first, second)
                       for other_name, other in kept if other_name != name]
        rest_summary = mantel_haenszel(rest_tables)
        rows.append((name, len(rest), odds_ratio(ra, rb, rc, rd), rest_summary))
    for name, n_left, crude_left, mh_left in sorted(rows, key=lambda r: r[2]):
        mh_text = "  n/a" if mh_left is None else "%7.2f" % mh_left
        print("    %-46s %7d %7.2f %s"
              % (name[:46], n_left, crude_left, mh_text))
    print()


def describe_unbatched(frame):
    """I print what the analysis table knows about samples carrying no batch."""
    missing = frame[frame["_batch"] == ""]
    print("=" * 72)
    print("SAMPLES WITH NO RUN BATCH: %d of %d" % (len(missing), len(frame)))
    print("=" * 72)
    if not len(missing):
        print("  none. Nothing to chase.")
        return
    print("  These are excluded from every summary above because they are not")
    print("  one stratum, they are an unknown number of them. They also carry")
    print("  the largest odds ratio on the board, so it matters what they are.")
    print()
    for column in ("Host_group", "Host_clean", "Host"):
        if column in missing.columns:
            counts = missing[column].fillna("(blank)").value_counts()
            print("  %s:" % column)
            for name, n in counts.items():
                print("      %-40s %4d" % (str(name)[:40], n))
            break
    if "NumTicks" in missing.columns:
        ticks = pd.to_numeric(missing["NumTicks"], errors="coerce")
        print("  pool size: median %g, range %g to %g, total ticks %g"
              % (ticks.median(), ticks.min(), ticks.max(), ticks.sum()))
    if "total_reads" in missing.columns:
        depth = pd.to_numeric(missing["total_reads"], errors="coerce")
        print("  read depth: median %g, range %g to %g"
              % (depth.median(), depth.min(), depth.max()))
    for pathogen in PATHOGENS:
        column = pathogen + "_pos"
        if column in missing.columns:
            print("  %-14s positive in %d of %d"
                  % (pathogen, int((missing[column] == 1).sum()), len(missing)))
    print()
    print("  sample ids, to grep in the run sheet:")
    for name in sorted(str(s) for s in missing["sample"]):
        print("      " + name)


def main():
    if not os.path.isfile(TABLE):
        sys.exit("not found: " + TABLE
                 + "\nRun this from the project root, after "
                   "python3 scripts/build_analysis_table.py")
    frame = pd.read_csv(TABLE, sep="\t")
    if BATCH_COLUMN not in frame.columns:
        sys.exit(TABLE + " has no '" + BATCH_COLUMN + "' column. Rebuild it:\n"
                 "    python3 scripts/build_analysis_table.py")
    if "sample" not in frame.columns:
        sys.exit(TABLE + " has no 'sample' column. Something is wrong with it.")
    frame = frame.copy()
    frame["_batch"] = [blank(value) for value in frame[BATCH_COLUMN]]

    present = [p for p in PATHOGENS
               if p + "_pos" in frame.columns and frame[p + "_pos"].sum() > 0]
    if len(present) < 2:
        sys.exit("fewer than two pathogens with detections in " + TABLE)

    print("samples: %d" % len(frame))
    print("pathogens with detections: %s" % ", ".join(present))
    print("distinct run batches: %d, plus %d samples with none"
          % (len(set(b for b in frame["_batch"] if b)),
             int((frame["_batch"] == "").sum())))
    print()

    for first, second in itertools.combinations(present, 2):
        report_pair(frame, first, second)

    describe_unbatched(frame)

    print()
    print("=" * 72)
    print("HOW TO READ THIS")
    print("=" * 72)
    print("  Breslow-Day small  ->  the batches disagree, and no single adjusted")
    print("                         odds ratio describes the data. The pooled")
    print("                         association is then a property of particular")
    print("                         runs, not of the pathogens, and the per-tick")
    print("                         analysis (n = 93) is the only co-occurrence")
    print("                         result that can be reported.")
    print()
    print("  Leave-one-out      ->  if dropping ONE batch takes the crude odds")
    print("                         ratio to about 1, the association is that")
    print("                         batch. If every row stays near the full-data")
    print("                         value, no single run is carrying it.")
    print()
    print("  Neither test can rescue this into co-infection. A pool is a tube,")
    print("  and 530 of 551 animals sit entirely within one batch (ground truth")
    print("  2b, lens 2), so animal and batch cannot be separated in this design.")


if __name__ == "__main__":
    main()
