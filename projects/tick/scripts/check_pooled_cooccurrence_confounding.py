#!/usr/bin/env python3
"""check_pooled_cooccurrence_confounding.py

I test whether the pooled pathogen co-occurrence signal survives conditioning
on pool size and sequencing depth. I write nothing and change nothing.

WHY
---
pathogen_cooccurrence.py, 2026-09-18, on 602 pooled samples:

    Anaplasma x Borreliella   104 observed vs 80.9 expected   OR 2.03  q = 0.0002
    Anaplasma x Borrelia       15 observed vs  8.9 expected   OR 3.08  q = 0.022

The same two pairs, measured per tick on the 93 single ticks by
recompute_cooccurrence.py the same day:

    Borreliella x Anaplasma    8 of 93   OR 1.31   p = 0.608   q = 0.912
    Anaplasma x Borrelia       0 of 93             p = 1.000

A pool holding twelve ticks is more likely to contain ANY given pathogen than
a pool holding one. That alone induces positive association between any two
reasonably prevalent taxa, with no biology involved. Sequencing depth does the
same: ground truth section 7 records Anaplasma positivity rising from 26% to
59% across depth quartiles.

So the pooled result is exactly what a confounded analysis looks like, and the
per-tick result is exactly what the unconfounded version looks like. This
script checks whether that is what is happening.

WHAT I DO
---------
For each pair I compute the crude odds ratio, then the odds ratio WITHIN
strata of pool size and within strata of read depth, and then the
Mantel-Haenszel summary odds ratio across those strata. Mantel-Haenszel is the
standard way to ask "is this association still there once I hold the
confounder constant", and it needs no modelling assumptions beyond the strata.

HOW TO READ IT
--------------
If the stratified odds ratios collapse toward 1 while the crude one stays at 2
or 3, the association is a pooling artifact and must not be reported as
co-occurrence. If they hold near the crude value inside every stratum, the
association is not explained by pool size or depth and deserves a closer look.

Either answer is worth having before anything is written. I do not decide
which; I print both and say what each would mean.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 check_pooled_cooccurrence_confounding.py
"""
import itertools
import os
import sys

import pandas as pd

TABLE = "results/analysis/analysis_table_16S.tsv"
PATHOGENS = ["Anaplasma", "Borreliella", "Borrelia", "Ehrlichia"]


def odds_ratio(a, b, c, d):
    """I return the Haldane-corrected odds ratio so a zero cell does not kill it."""
    if 0 in (a, b, c, d):
        return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return (a * d) / (b * c)


def two_by_two(frame, first, second):
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


def strata(frame, column, labels):
    """I cut a column into quartiles, dropping strata that cannot inform."""
    values = pd.to_numeric(frame[column], errors="coerce")
    try:
        binned = pd.qcut(values, 4, labels=labels, duplicates="drop")
    except ValueError:
        return None
    return binned


def report(frame, pair, column, label):
    first, second = pair
    binned = strata(frame, column, ["Q1", "Q2", "Q3", "Q4"])
    if binned is None:
        print("    %-14s could not be split into strata" % label)
        return
    tables = []
    print("    by %s:" % label)
    for level in binned.cat.categories:
        sub = frame[binned == level]
        if len(sub) < 10:
            continue
        a, b, c, d = two_by_two(sub, first, second)
        tables.append((a, b, c, d))
        lo = pd.to_numeric(sub[column], errors="coerce").min()
        hi = pd.to_numeric(sub[column], errors="coerce").max()
        print("      %-4s n=%-4d %-18s both=%-4d OR=%.2f"
              % (level, len(sub), "[%g, %g]" % (lo, hi), a,
                 odds_ratio(a, b, c, d)))
    summary = mantel_haenszel(tables)
    if summary is not None:
        print("      Mantel-Haenszel summary OR, holding %s constant: %.2f"
              % (label, summary))


def report_categorical(frame, pair, column, label, min_n=10):
    """I stratify on a categorical column, grouping by value rather than binning.

    pd.qcut is for numbers. Pushing a categorical through it would either fail
    or order the strata by string, which is exactly the sort of quiet wrongness
    this script exists to catch.

    Samples with no value get their own line but do NOT enter the summary. They
    are not one stratum, they are an unknown number of strata pooled together,
    so including them would put the very confounding I am trying to remove back
    into the answer.
    """
    first, second = pair
    if column not in frame.columns:
        print("    %-14s column not present; rebuild the analysis table"
              % label)
        return
    work = frame.copy()
    # Not .astype(str).replace(): pandas turns the replacement of the string
    # "nan" back into a real NaN, so the blanks vanish and the missing group
    # silently rejoins the strata. This project has been bitten by the
    # string-versus-NA distinction twice already. A comprehension has no
    # opinions.
    work["_stratum"] = ["" if str(v).strip().lower() in ("nan", "none",
                                                         "<na>", "")
                        else str(v).strip() for v in work[column]]
    missing = work[work["_stratum"] == ""]
    work = work[work["_stratum"] != ""]

    tables, skipped, informative = [], [], 0
    print("    by %s:" % label)
    for level, sub in work.groupby("_stratum", observed=True):
        if len(sub) < min_n:
            skipped.append((str(level), len(sub)))
            continue
        a, b, c, d = two_by_two(sub, first, second)
        tables.append((a, b, c, d))
        if b > 0 and c > 0:
            informative += 1
        print("      %-46s n=%-4d both=%-4d OR=%.2f"
              % (str(level)[:46], len(sub), a, odds_ratio(a, b, c, d)))
    for level, n in skipped:
        print("      %-46s n=%-4d SKIPPED, below %d"
              % (str(level)[:46], n, min_n))
    if len(missing):
        a, b, c, d = two_by_two(missing, first, second)
        print("      %-46s n=%-4d both=%-4d OR=%.2f  NOT in the summary"
              % ("(no " + label + " recorded)", len(missing), a,
                 odds_ratio(a, b, c, d)))
    summary = mantel_haenszel(tables)
    if summary is None:
        print("      no stratum large enough to summarise")
        return
    print("      Mantel-Haenszel summary OR, holding %s constant: %.2f"
          % (label, summary))
    print("      strata in the summary: %d, covering %d of %d samples"
          % (len(tables), sum(sum(t) for t in tables), len(frame)))
    print("      informative strata (both discordant cells non-zero): "
          "%d of %d" % (informative, len(tables)))
    if informative < 3:
        print("      CAUTION: fewer than three informative strata. This summary"
              " rests on very little.")


def main():
    if not os.path.isfile(TABLE):
        sys.exit("not found: " + TABLE
                 + "\nRun build_analysis_table.py first.")
    frame = pd.read_csv(TABLE, sep="\t")
    for column in ("NumTicks", "total_reads"):
        if column not in frame.columns:
            sys.exit("column '" + column + "' missing from " + TABLE)

    present = [p for p in PATHOGENS
               if p + "_pos" in frame.columns
               and frame[p + "_pos"].sum() > 0]
    print("samples: %d" % len(frame))
    print("pathogens with detections: %s" % ", ".join(present))
    print("pool size: median %g, range %g to %g"
          % (frame["NumTicks"].median(), frame["NumTicks"].min(),
             frame["NumTicks"].max()))
    print("read depth: median %g, range %g to %g"
          % (frame["total_reads"].median(), frame["total_reads"].min(),
             frame["total_reads"].max()))
    print()

    for first, second in itertools.combinations(present, 2):
        a, b, c, d = two_by_two(frame, first, second)
        crude = odds_ratio(a, b, c, d)
        print("%s x %s" % (first, second))
        print("    crude: both=%d  %s only=%d  %s only=%d  neither=%d  OR=%.2f"
              % (a, first, b, second, c, d, crude))
        report(frame, (first, second), "NumTicks", "pool size")
        report(frame, (first, second), "total_reads", "read depth")
        report_categorical(frame, (first, second), "RunBatch", "run batch")
        print()

    print("HOW TO READ THIS")
    print("  If the stratified and Mantel-Haenszel odds ratios collapse toward 1")
    print("  while the crude one stays at 2 or 3, the pooled association is an")
    print("  artifact of pool size or depth and must NOT be reported as pathogen")
    print("  co-occurrence. The per-tick analysis (n = 93,")
    print("  recompute_cooccurrence.py) is then the only co-occurrence result,")
    print("  and it finds no pair above chance.")
    print()
    print("  If they hold near the crude value inside every stratum, the")
    print("  association is not explained by these two confounders. That is a")
    print("  finding, but it is still a POOL-level co-detection and not")
    print("  co-infection of one tick.")
    print()
    print("  BATCH is now tested too. It is the largest effect in this dataset")
    print("  (R2 10.9%, p = 0.001) and 530 of 551 animals sit entirely within one")
    print("  batch (ground truth 2b, lens 2), so animal and batch cannot be")
    print("  separated here. If the association survives batch as well, the three")
    print("  obvious artifacts are ruled out and it is a real pool-level")
    print("  co-detection. It is STILL not co-infection of one tick, because a")
    print("  pool is a tube, and the animal-batch confound means it still cannot")
    print("  be attributed to tick biology rather than to something shared by")
    print("  ticks off the same animal.")


if __name__ == "__main__":
    main()
