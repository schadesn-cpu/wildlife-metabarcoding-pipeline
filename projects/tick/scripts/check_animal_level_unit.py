#!/usr/bin/env python3
"""check_animal_level_unit.py

I measure prevalence and co-occurrence with the ANIMAL as the unit, four ways,
and I show what NumTicks means now that we know it is repeated. I write
nothing and change nothing.

WHY THIS EXISTS
---------------
check_specimen_replication.py, 2026-09-18: the 602 rows of
analysis_table_16S.tsv are 413 animals. 102 animals hold more than one row,
accounting for 291 of 602 rows, or 48.3%. Barcode and the library lookup agree
exactly on that grouping, so it is not a resolver artifact.

The analysis table emits one row per 16S library because nothing ever told it
otherwise. 602 was never chosen as an n. It is the number of rows in a file.

That leaves two analyses on the table and NEITHER is clean:

    per row (n=602)       291 rows are not independent of another row
    per animal, any-pos   an animal with ten rows has ten chances to be called
                          (n=413)   positive and an animal with one row has one

The second is why collapsing raised prevalence 6 to 7 points and took
Anaplasma x Borrelia from 3.08 to 8.00. It replaced pseudoreplication with a
detection-opportunity artifact. So the association has still not been measured
once, and that is what this script is for.

WHAT I DO
---------
1. WHAT NumTicks MEANS. Every repeated animal seems to carry one tick count
   copied onto each of its rows. I check that across ALL repeated animals, not
   the first forty, and I total the column both ways so the size of the
   overcount is visible. If NumTicks is the animal's total repeated, then any
   analysis weighted by pool size is weighted by a number that does not mean
   what it says, and that includes PoolTestR, whose estimator IS pool size.

2. DETECTION OPPORTUNITY. Positivity against rows-per-animal. If animals with
   more rows are called positive more often, the any-positive collapse is
   biased and by how much.

3. THE CLEAN SUBSET. Animals holding exactly one row. No clustering, no
   collapse, no opportunity bias. Smaller, and honest.

4. ONE ROW PER ANIMAL, DRAWN AT RANDOM, many times. This uses all 413 animals,
   gives each exactly one vote, and the spread across draws is the uncertainty
   that choosing a row introduces. Seed fixed so the numbers repeat.

Then all four side by side.

WHAT I DO NOT DO
----------------
I do not pick the unit and I do not write a corrected table. Which rows of one
animal to combine, and how, is a decision about the experiment, not about this
file.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_animal_level_unit.py
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
# The crosswalk's tick total, for comparison against the NumTicks column.
# From build_4marker_crosswalk.py, 2026-09-18: 563 specimens, 3,939 ticks.
CROSSWALK_TICKS = 3939

RUN_SUFFIX = re.compile(r"\..*$")


def blank(value):
    """I return '' for anything meaning 'no value', and the text otherwise."""
    text = str(value).strip()
    return "" if text.lower() in ("nan", "none", "<na>", "") else text


def resolve_animals(frame):
    """I return a specimen label per row, and the name of the scheme I used.

    I prefer the library lookup, then the Barcode column. I do NOT fall back to
    a name prefix here: stripping a trailing index reads EQ02441291-10
    correctly and reads ALHT-Deer-10 as ALHT-Deer, which merges two animals.
    check_specimen_replication.py reports that comparison; this script needs an
    identity it can trust, so if neither is available I stop.
    """
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
    sys.exit("cannot identify animals in " + TABLE + ": the library lookup "
             "did not resolve and there is no usable Barcode column.")


def numticks_report(frame, key):
    """I report what the NumTicks column means across repeated animals."""
    print("=" * 72)
    print("1. WHAT DOES NumTicks MEAN?")
    print("=" * 72)
    if "NumTicks" not in frame.columns:
        print("  no NumTicks column. Nothing to check.")
        print()
        return
    ticks = pd.to_numeric(frame["NumTicks"], errors="coerce")
    work = frame.assign(_ticks=ticks)
    usable = work[(work[key] != "") & work["_ticks"].notna()]
    grouped = usable.groupby(key)["_ticks"]
    identical = differing = 0
    for _, values in grouped:
        if values.nunique() == 1:
            identical += 1
        else:
            differing += 1
    repeated = grouped.size()
    repeated_names = repeated[repeated > 1].index
    identical_repeated = sum(
        1 for name in repeated_names
        if usable.loc[usable[key] == name, "_ticks"].nunique() == 1)
    print("  across ALL %d animals (not the first forty):" % len(repeated))
    print("    %d carry one NumTicks on every row" % identical)
    print("    %d carry different NumTicks between rows" % differing)
    print("  among the %d animals with MORE than one row:" % len(repeated_names))
    print("    %d carry one NumTicks on every row" % identical_repeated)
    print("    %d carry different NumTicks between rows"
          % (len(repeated_names) - identical_repeated))
    print()
    row_total = float(ticks.sum())
    animal_total = float(usable.groupby(key)["_ticks"].max().sum())
    print("  NumTicks summed over %d rows:    %.0f" % (len(frame), row_total))
    print("  NumTicks summed over %d animals: %.0f  (max per animal)"
          % (usable[key].nunique(), animal_total))
    print("  crosswalk total ticks:            %d" % CROSSWALK_TICKS)
    print()
    if row_total > animal_total * 1.05:
        print("  The row total exceeds the animal total by %.0f ticks (%.0f%%)."
              % (row_total - animal_total,
                 100.0 * (row_total - animal_total) / max(animal_total, 1)))
        print("  So NumTicks is the ANIMAL's tick count repeated onto each of")
        print("  its rows, not the count of ticks in that particular tube.")
        print("  Anything weighted by pool size is therefore weighted by a")
        print("  number that does not mean what it says. That includes the")
        print("  pool-size strata in the confounding check and PoolTestR,")
        print("  whose estimator is pool size.")
    else:
        print("  The two totals are close, so NumTicks varies between an")
        print("  animal's rows and is plausibly the count for each tube.")
    print()


def opportunity_report(frame, key, present):
    """I show whether more rows makes an animal more likely to be positive."""
    print("=" * 72)
    print("2. DOES HAVING MORE ROWS MAKE AN ANIMAL 'POSITIVE' MORE OFTEN?")
    print("=" * 72)
    print("  Collapsing with ANY-positive gives an animal with ten rows ten")
    print("  chances and an animal with one row one chance. If these columns")
    print("  rise from left to right, that is the artifact, measured.")
    print()
    usable = frame[frame[key] != ""]
    columns = [p + "_pos" for p in present]
    collapsed = usable.groupby(key)[columns].max()
    sizes = usable.groupby(key).size()
    buckets = [(1, 1, "1 row"), (2, 2, "2 rows"), (3, 3, "3 rows"),
               (4, 10 ** 6, "4+ rows")]
    print("  %-14s %s" % ("", "  ".join("%-12s" % label
                                        for _, _, label in buckets)))
    print("  %-14s %s" % ("animals", "  ".join(
        "%-12d" % int(((sizes >= low) & (sizes <= high)).sum())
        for low, high, _ in buckets)))
    for pathogen in present:
        cells = []
        for low, high, _ in buckets:
            names = sizes[(sizes >= low) & (sizes <= high)].index
            subset = collapsed.loc[names, pathogen + "_pos"]
            if not len(subset):
                cells.append("%-12s" % "n/a")
            else:
                cells.append("%-12s" % ("%.1f%% (%d)"
                                        % (100.0 * (subset == 1).mean(),
                                           int((subset == 1).sum()))))
        print("  %-14s %s" % (pathogen, "  ".join(cells)))
    print()


def pair_table(source, present, label):
    """I return {pair: (odds ratio, p, both)} for one already-collapsed frame."""
    out = {}
    for first, second in itertools.combinations(present, 2):
        x = source[first + "_pos"] == 1
        y = source[second + "_pos"] == 1
        a = int((x & y).sum())
        b = int((x & ~y).sum())
        c = int((~x & y).sum())
        d = int((~x & ~y).sum())
        ratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
        out[(first, second)] = (ratio, pvalue, a)
    return out


def singleton_report(frame, key, present):
    """I analyse only the animals holding exactly one row."""
    print("=" * 72)
    print("3. THE CLEAN SUBSET: ANIMALS WITH EXACTLY ONE ROW")
    print("=" * 72)
    usable = frame[frame[key] != ""]
    sizes = usable.groupby(key).size()
    singles = usable[usable[key].isin(sizes[sizes == 1].index)]
    print("  %d animals, each contributing one row. No clustering, no"
          % len(singles))
    print("  collapsing, no detection-opportunity bias. Just smaller.")
    print()
    for pathogen in present:
        column = pathogen + "_pos"
        k = int((singles[column] == 1).sum())
        low, high = stats.binomtest(k, len(singles)).proportion_ci(
            method="wilson")
        print("    %-14s %.1f%% (%d/%d)  95%% CI %.1f-%.1f"
              % (pathogen, 100.0 * k / len(singles), k, len(singles),
                 100 * low, 100 * high))
    print()
    return singles, pair_table(singles, present, "singletons")


def resample_report(frame, key, present):
    """I draw one row per animal, many times, and report the spread.

    Every animal gets exactly one vote, nothing is discarded, and the variation
    across draws is the uncertainty that choosing a row introduces. That
    uncertainty is real and belongs in the answer.
    """
    print("=" * 72)
    print("4. ONE ROW PER ANIMAL, DRAWN AT RANDOM, %d TIMES" % DRAWS)
    print("=" * 72)
    usable = frame[frame[key] != ""].reset_index(drop=True)
    groups = [np.asarray(positions)
              for positions in usable.groupby(key).indices.values()]
    lengths = np.array([len(g) for g in groups])
    widest = int(lengths.max())
    padded = np.zeros((len(groups), widest), dtype=np.int64)
    for i, g in enumerate(groups):
        padded[i, :len(g)] = g
    positives = {p: (usable[p + "_pos"].values == 1) for p in present}
    generator = np.random.default_rng(SEED)
    rows = np.arange(len(groups))

    prevalence = {p: [] for p in present}
    ratios = {pair: [] for pair in itertools.combinations(present, 2)}
    significant = {pair: 0 for pair in ratios}
    for _ in range(DRAWS):
        choice = (generator.random(len(groups)) * lengths).astype(np.int64)
        picked = padded[rows, choice]
        drawn = {p: positives[p][picked] for p in present}
        for p in present:
            prevalence[p].append(100.0 * drawn[p].mean())
        for first, second in ratios:
            x, y = drawn[first], drawn[second]
            a = int(np.count_nonzero(x & y))
            b = int(np.count_nonzero(x & ~y))
            c = int(np.count_nonzero(~x & y))
            d = int(np.count_nonzero(~x & ~y))
            ratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
            ratios[(first, second)].append(ratio)
            if pvalue < 0.05:
                significant[(first, second)] += 1

    print("  %d animals, one row each, %d draws. Seed %d."
          % (len(groups), DRAWS, SEED))
    print()
    print("  prevalence, median across draws (2.5th to 97.5th percentile):")
    for p in present:
        values = np.array(prevalence[p])
        print("    %-14s %.1f%%  (%.1f to %.1f)"
              % (p, np.median(values), np.percentile(values, 2.5),
                 np.percentile(values, 97.5)))
    print()
    print("  odds ratios, median across draws:")
    summary = {}
    for pair in ratios:
        values = np.array([r for r in ratios[pair] if np.isfinite(r)])
        if not len(values):
            continue
        share = 100.0 * significant[pair] / DRAWS
        summary[pair] = (float(np.median(values)),
                         float(np.percentile(values, 2.5)),
                         float(np.percentile(values, 97.5)), share)
        print("    %-26s %.2f  (%.2f to %.2f)   p<0.05 in %.0f%% of draws"
              % ("%s x %s" % pair, summary[pair][0], summary[pair][1],
                 summary[pair][2], share))
    print()
    return summary


def comparison(frame, key, present, singles_pairs, resampled):
    """I print all four estimates of each pair beside each other."""
    usable = frame[frame[key] != ""]
    columns = [p + "_pos" for p in present]
    collapsed = usable.groupby(key)[columns].max()
    rows_pairs = pair_table(frame, present, "rows")
    anypos_pairs = pair_table(collapsed, present, "any-positive")
    print("=" * 72)
    print("ALL FOUR, SIDE BY SIDE")
    print("=" * 72)
    print("  %-26s %-13s %-13s %-13s %s"
          % ("pair", "per row", "any-positive", "one-row only", "resampled"))
    print("  %-26s %-13s %-13s %-13s %s"
          % ("", "n=%d" % len(frame), "n=%d" % len(collapsed),
             "n=%d" % int((usable.groupby(key).size() == 1).sum()),
             "n=%d" % usable[key].nunique()))
    for pair in itertools.combinations(present, 2):
        pieces = []
        for source in (rows_pairs, anypos_pairs, singles_pairs):
            ratio, pvalue, _ = source[pair]
            pieces.append("%-13s" % ("%.2f p=%.3g" % (ratio, pvalue)))
        if pair in resampled:
            median, low, high, share = resampled[pair]
            pieces.append("%.2f (%.2f-%.2f)" % (median, low, high))
        else:
            pieces.append("n/a")
        print("  %-26s %s" % ("%s x %s" % pair, " ".join(pieces)))
    print()
    print("  per row       pseudoreplicated: 291 of 602 rows share an animal")
    print("  any-positive  biased upward by unequal detection opportunity")
    print("  one-row only  clean, and the smallest")
    print("  resampled     clean, uses every animal, interval covers the")
    print("                choice of row")
    print()
    print("  A pair that is strong in the first two and gone in the last two")
    print("  was an artifact of the unit. A pair that holds in all four is")
    print("  not explained by the unit, though pool size, depth and batch")
    print("  still apply and the batches did not share one odds ratio.")


def main():
    if not os.path.isfile(TABLE):
        sys.exit("not found: " + TABLE
                 + "\nRun this from the project root, after "
                   "python3 scripts/build_analysis_table.py")
    frame = pd.read_csv(TABLE, sep="\t")
    if "sample" not in frame.columns:
        sys.exit(TABLE + " has no 'sample' column.")
    frame = frame.copy()
    labels, scheme = resolve_animals(frame)
    frame["_animal"] = labels
    key = "_animal"

    present = [p for p in PATHOGENS + CO_DETECTIONS
               if p + "_pos" in frame.columns and frame[p + "_pos"].sum() > 0]
    if len(present) < 2:
        sys.exit("fewer than two organisms with detections in " + TABLE)

    unresolved = int((frame[key] == "").sum())
    print("rows: %d" % len(frame))
    print("animals: %d, identified by the %s"
          % (frame.loc[frame[key] != "", key].nunique(), scheme))
    if unresolved:
        print("rows with no animal: %d. They are left out of everything below."
              % unresolved)
    print("organisms with detections: %s" % ", ".join(present))
    print()

    numticks_report(frame, key)
    opportunity_report(frame, key, present)
    singles, singles_pairs = singleton_report(frame, key, present)
    resampled = resample_report(frame, key, present)
    comparison(frame, key, present, singles_pairs, resampled)


if __name__ == "__main__":
    main()
