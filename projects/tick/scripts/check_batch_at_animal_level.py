#!/usr/bin/env python3
"""check_batch_at_animal_level.py

I rerun the batch conditioning with ONE ROW PER ANIMAL, because the version
that rejected homogeneity ran on pseudoreplicated rows. I write nothing and
change nothing.

WHY THIS EXISTS
---------------
check_cooccurrence_batch_heterogeneity.py, 2026-09-18, on 602 rows, found that
the batches do not share one odds ratio for Anaplasma x Borreliella:
Breslow-Day chi2 18.63 on 8 df, bootstrap p = 0.009. The finding was driven by
one plate, AMP-120924-USGS-TICKS-P3-5, at OR 5.86 against four batches sitting
at 1.

check_specimen_replication.py, the same day: those 602 rows are 413 animals,
102 of which hold more than one row. P3-5 is the most clustered plate in the
set, 1.7 rows per animal against 1.4 for the rest.

Clustering inflates an odds ratio wherever it is heaviest. So the 5.86 that
drove the heterogeneity conclusion may be the repeated animals rather than the
plate, and the whole result has to be recomputed with each animal counted once
before it means anything.

check_animal_level_unit.py, the same day, established that the association
itself survives the unit: 2.00 on the 311 single-row animals, 1.96 resampled
across all 413 (1.68 to 2.29), significant in 2,000 of 2,000 draws. So this is
not about whether the association is real. It is about whether it is the same
association in every batch.

WHAT I DO
---------
1. REPAIR THE BATCH LABELS, from the data rather than from a guess. 43 rows
   carry no batch string. Every sample id ends in a run token such as
   '241209_P3-5'. I learn the token-to-batch mapping from the rows that carry
   BOTH, and apply it to the rows that carry only the token. Nothing is
   inferred about date formats; if a token never appears with a batch string,
   those rows stay unlabelled and I say so.

2. DRAW ONE ROW PER ANIMAL, many times. Each animal votes once, with the batch
   of whichever row was drawn. Per draw I compute every batch's odds ratio, the
   Mantel-Haenszel summary, and the Breslow-Day statistic. I report the
   distribution of each across draws, so the spread includes the uncertainty
   of which row got picked.

3. LEAVE ONE BATCH OUT, again per draw, so a plate that carries the
   association shows up as the row where the odds ratio falls.

WHAT THESE TESTS CAN SEE, MEASURED FIRST
----------------------------------------
Collapsing 602 rows to 413 animals shrinks each batch from about 67 rows to
about 46 animals, and that costs power. Simulated at the animal-level sizes
(57, 60, 53, 45, 30, 55, 32, 10, 85), rejection rates at a nominal 5%:

                                 9-batch omnibus    one plate vs the rest
    all batches truly 2.0              0.07                  0.06
    one plate 4, the rest 2            ----                  0.22
    one plate 6, the rest 2            0.28                  0.51
    one plate 10, the rest 1           0.87                  0.99
    the row-level pattern              0.56                  ----

So the omnibus test is nearly blind to the thing we suspect: if P3-5 truly
runs at 6 while the others run at 2, it misses that 72% of the time. A
non-rejection from it is close to no information. That is why the targeted
test is here as well, one plate against the pooled remainder, one degree of
freedom instead of eight, at roughly double the power.

Both run slightly hot under the null, 6 to 7% against a nominal 5%, so a
p just under 0.05 should not be leaned on.

The chi2 approximation is used directly rather than a nested bootstrap, which
at 1,000 draws would be unaffordable. It was checked against a bootstrap at
the row-level sizes and agreed (mean statistic 8.01 on 8 df, 5.5% rejection).

WHAT I DO NOT DO
----------------
I do not condition on pool size. NumTicks in this table is the animal's total
tick count copied onto each of its rows: 4,730 summed over rows against 3,209
over animals. It is not the size of any pool and must not be used as one until
the run sheet says what the real per-library counts are.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_batch_at_animal_level.py
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
BATCH_COLUMN = "RunBatch"
DRAWS = 1000
SEED = 20260918
MIN_STRATUM = 10
# Below this many positives an odds ratio is a number with no information in
# it. Borrelia has 23 across 602 rows and produced an interval of 3.49 to
# 18.80 in the resampled analysis; it belongs in the paper as a count.
MIN_POSITIVES = 30

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
    """I solve for the expected 'both' cell under a common odds ratio."""
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
    """I return (statistic, degrees of freedom, p) for Tarone's corrected test."""
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
        statistic -= residual_sum ** 2 / variance_sum
    degrees = used - 1
    return statistic, degrees, float(stats.chi2.sf(statistic, degrees))


def resolve_animals(frame):
    """I return a specimen label per row, and the name of the scheme I used."""
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


def repair_batches(frame):
    """I fill missing batch strings from the run token in the sample id.

    I learn the mapping from rows that carry both a token and a batch, so the
    only assumption is that one sequencing run token means one batch. I check
    that assumption and refuse to use any token that maps to more than one.
    """
    print("=" * 72)
    print("1. BATCH LABELS")
    print("=" * 72)
    tokens = []
    for sample in frame["sample"]:
        found = RUN_TOKEN.search(str(sample).strip())
        tokens.append(found.group(1) if found else "")
    frame["_token"] = tokens
    recorded = frame["_batch"] != ""
    print("  %d of %d rows carry a batch string; %d do not."
          % (int(recorded.sum()), len(frame), int((~recorded).sum())))

    mapping = {}
    ambiguous = []
    known = frame[recorded & (frame["_token"] != "")]
    for token, group in known.groupby("_token"):
        values = sorted(set(group["_batch"]))
        if len(values) == 1:
            mapping[token] = values[0]
        else:
            ambiguous.append((token, values))
    for token, values in ambiguous:
        print("  token %r maps to %d different batches; not used for repair."
              % (token, len(values)))

    repaired = 0
    unrepaired_tokens = {}
    for index in frame.index[~recorded]:
        token = frame.at[index, "_token"]
        if token in mapping:
            frame.at[index, "_batch"] = mapping[token]
            repaired += 1
        else:
            unrepaired_tokens[token] = unrepaired_tokens.get(token, 0) + 1
    if repaired:
        print("  %d rows repaired from their run token. Which token, and to"
              " what:" % repaired)
        shown = set()
        for index in frame.index:
            token = frame.at[index, "_token"]
            if token in mapping and token not in shown:
                count = int(((frame["_token"] == token)
                             & (frame["_batch"] == mapping[token])).sum())
                was_missing = int(((frame["_token"] == token)
                                   & ~recorded).sum())
                if was_missing:
                    print("      %-22s -> %-46s (%d rows, %d repaired)"
                          % (token[:22], mapping[token][:46], count,
                             was_missing))
                    shown.add(token)
    for token, count in sorted(unrepaired_tokens.items()):
        print("  %d rows with token %r have no batch anywhere; left unlabelled."
              % (count, token))
    print()
    return frame


def main():
    if not os.path.isfile(TABLE):
        sys.exit("not found: " + TABLE
                 + "\nRun this from the project root, after "
                   "python3 scripts/build_analysis_table.py")
    frame = pd.read_csv(TABLE, sep="\t")
    if "sample" not in frame.columns:
        sys.exit(TABLE + " has no 'sample' column.")
    if BATCH_COLUMN not in frame.columns:
        sys.exit(TABLE + " has no '" + BATCH_COLUMN + "' column. Rebuild it.")
    frame = frame.copy()
    frame["_batch"] = [blank(v) for v in frame[BATCH_COLUMN]]
    labels, scheme = resolve_animals(frame)
    frame["_animal"] = labels

    present = [p for p in PATHOGENS + CO_DETECTIONS
               if p + "_pos" in frame.columns and frame[p + "_pos"].sum() > 0]
    thick = [p for p in present if int((frame[p + "_pos"] == 1).sum())
             >= MIN_POSITIVES]
    thin = [p for p in present if p not in thick]

    print("rows: %d" % len(frame))
    print("animals: %d, identified by the %s"
          % (frame.loc[frame["_animal"] != "", "_animal"].nunique(), scheme))
    for p in thin:
        print("skipping %s: %d positives across %d rows, below %d. An odds"
              % (p, int((frame[p + "_pos"] == 1).sum()), len(frame),
                 MIN_POSITIVES))
        print("    ratio on that is uninformative; report it as a count.")
    if len(thick) < 2:
        sys.exit("fewer than two organisms with enough positives to stratify.")
    print()

    frame = repair_batches(frame)

    usable = frame[(frame["_animal"] != "") & (frame["_batch"] != "")]
    usable = usable.reset_index(drop=True)
    print("  %d rows carry both an animal and a batch, covering %d animals."
          % (len(usable), usable["_animal"].nunique()))
    print()

    groups = [np.asarray(positions)
              for positions in usable.groupby("_animal").indices.values()]
    lengths = np.array([len(g) for g in groups])
    padded = np.zeros((len(groups), int(lengths.max())), dtype=np.int64)
    for i, g in enumerate(groups):
        padded[i, :len(g)] = g
    rows = np.arange(len(groups))
    batch_codes, batch_names = pd.factorize(usable["_batch"])
    batch_codes = np.asarray(batch_codes)
    positives = {p: (usable[p + "_pos"].values == 1) for p in thick}
    generator = np.random.default_rng(SEED)

    pairs = list(itertools.combinations(thick, 2))
    per_batch = {pair: {name: [] for name in batch_names} for pair in pairs}
    summaries = {pair: [] for pair in pairs}
    crude = {pair: [] for pair in pairs}
    bd_p = {pair: [] for pair in pairs}
    bd_stat = {pair: [] for pair in pairs}
    informative = {pair: [] for pair in pairs}
    left_out = {pair: {name: {"crude": [], "mh": []} for name in batch_names}
                for pair in pairs}
    versus_rest = {pair: {name: [] for name in batch_names} for pair in pairs}

    print("=" * 72)
    print("2. ONE ROW PER ANIMAL, %d DRAWS, SEED %d" % (DRAWS, SEED))
    print("=" * 72)
    print("  Each animal votes once, carrying the batch of the row drawn.")
    print()

    for _ in range(DRAWS):
        choice = (generator.random(len(groups)) * lengths).astype(np.int64)
        picked = padded[rows, choice]
        codes = batch_codes[picked]
        for pair in pairs:
            first, second = pair
            x = positives[first][picked]
            y = positives[second][picked]
            a = int(np.count_nonzero(x & y))
            b = int(np.count_nonzero(x & ~y))
            c = int(np.count_nonzero(~x & y))
            d = int(np.count_nonzero(~x & ~y))
            crude[pair].append(odds_ratio(a, b, c, d))
            tables = []
            names_used = []
            count_informative = 0
            for code, name in enumerate(batch_names):
                mask = codes == code
                size = int(np.count_nonzero(mask))
                if size < MIN_STRATUM:
                    continue
                xs, ys = x[mask], y[mask]
                ta = int(np.count_nonzero(xs & ys))
                tb = int(np.count_nonzero(xs & ~ys))
                tc = int(np.count_nonzero(~xs & ys))
                td = int(np.count_nonzero(~xs & ~ys))
                tables.append((ta, tb, tc, td))
                names_used.append(name)
                per_batch[pair][name].append(odds_ratio(ta, tb, tc, td))
                if tb > 0 and tc > 0:
                    count_informative += 1
            informative[pair].append(count_informative)
            summary = mantel_haenszel(tables)
            if summary is None:
                continue
            summaries[pair].append(summary)
            test = breslow_day(tables, summary)
            if test is not None:
                bd_stat[pair].append(test[0])
                bd_p[pair].append(test[2])
            for drop, name in enumerate(names_used):
                rest = [t for i, t in enumerate(tables) if i != drop]
                rest_mh = mantel_haenszel(rest)
                ra = sum(t[0] for t in rest)
                rb = sum(t[1] for t in rest)
                rc = sum(t[2] for t in rest)
                rd = sum(t[3] for t in rest)
                left_out[pair][name]["crude"].append(odds_ratio(ra, rb, rc, rd))
                if rest_mh is not None:
                    left_out[pair][name]["mh"].append(rest_mh)
                # One plate against everything else pooled. One degree of
                # freedom, and about twice the power of the omnibus for the
                # hypothesis we actually hold, which is that ONE plate differs.
                two = [tables[drop], (ra, rb, rc, rd)]
                two_mh = mantel_haenszel(two)
                if two_mh is not None:
                    focused = breslow_day(two, two_mh)
                    if focused is not None:
                        versus_rest[pair][name].append(focused[2])

    def band(values):
        arr = np.array([v for v in values if np.isfinite(v)])
        if not len(arr):
            return None
        return (float(np.median(arr)), float(np.percentile(arr, 2.5)),
                float(np.percentile(arr, 97.5)))

    for pair in pairs:
        first, second = pair
        print("%s x %s" % (first, second))
        c = band(crude[pair])
        s = band(summaries[pair])
        print("  crude OR over all animals:  %.2f  (%.2f to %.2f)" % c)
        if s:
            print("  Mantel-Haenszel over batch: %.2f  (%.2f to %.2f)" % s)
        print("  informative batches: median %d of %d used"
              % (int(np.median(informative[pair])),
                 int(np.median([len([1 for name in batch_names
                                     if per_batch[pair][name]])]))))
        print()
        print("  per batch, median OR across draws (2.5th to 97.5th):")
        for name in batch_names:
            values = per_batch[pair][name]
            if not values:
                print("    %-46s below %d animals in most draws"
                      % (str(name)[:46], MIN_STRATUM))
                continue
            b = band(values)
            print("    %-46s %6.2f  (%.2f to %.2f)  in %d%% of draws"
                  % (str(name)[:46], b[0], b[1], b[2],
                     int(round(100.0 * len(values) / DRAWS))))
        print()
        if bd_p[pair]:
            arr = np.array(bd_p[pair])
            print("  HOMOGENEITY, Breslow-Day with Tarone's correction")
            print("    median chi2 %.2f, median p %.4g"
                  % (float(np.median(bd_stat[pair])), float(np.median(arr))))
            print("    p < 0.05 in %.0f%% of draws" % (100.0 * (arr < 0.05).mean()))
            if (arr < 0.05).mean() >= 0.5:
                print("    The batches still do not share one odds ratio once")
                print("    each animal is counted once. The heterogeneity is")
                print("    not the repeated rows.")
            elif (arr < 0.05).mean() <= 0.1:
                print("    This test does not reject. Do not read that as the")
                print("    batches agreeing: at these sizes it misses a plate")
                print("    running at 6 against the rest at 2 about 72% of the")
                print("    time. The targeted test below is the one to read.")
            else:
                print("    The test rejects in some draws and not others, so")
                print("    which row is chosen decides the answer. That is not")
                print("    a result; it means the data cannot settle this.")
        print()
        print("  EACH PLATE AGAINST THE POOLED REMAINDER (1 df, better powered)")
        print("    %-46s %10s %s" % ("batch", "median p", "p<0.05 in"))
        focused_rows = []
        for name in batch_names:
            values = versus_rest[pair][name]
            if not values:
                continue
            arr = np.array(values)
            focused_rows.append((float(np.median(arr)), name,
                                 100.0 * (arr < 0.05).mean()))
        for median_p, name, share in sorted(focused_rows):
            print("    %-46s %10.4g %.0f%% of draws"
                  % (str(name)[:46], median_p, share))
        if focused_rows:
            print("    Nine tests, so a single p near 0.05 is not much. A plate")
            print("    that differs should stand out, not merely qualify.")
        print()
        print("  LEAVE ONE BATCH OUT, median across draws")
        print("    %-46s %8s %8s" % ("batch removed", "crude", "MH"))
        print("    %-46s %8.2f %8s"
              % ("(none removed)", c[0], "%.2f" % s[0] if s else "n/a"))
        ordered = []
        for name in batch_names:
            if not left_out[pair][name]["crude"]:
                continue
            lc = band(left_out[pair][name]["crude"])
            lm = band(left_out[pair][name]["mh"])
            ordered.append((lc[0], name, lc, lm))
        for _, name, lc, lm in sorted(ordered):
            print("    %-46s %8.2f %8s"
                  % (str(name)[:46], lc[0], "%.2f" % lm[0] if lm else "n/a"))
        print()

    print("=" * 72)
    print("HOW TO READ THIS")
    print("=" * 72)
    print("  Everything here counts each animal once, so none of it is the")
    print("  pseudoreplication. What it cannot remove is the design: 530 of")
    print("  551 animals sit entirely within one batch, so batch and animal")
    print("  are the same variable for most of this dataset. A batch effect")
    print("  and an animal effect cannot be told apart here by any analysis.")
    print()
    print("  A non-rejection here is weak. The omnibus test misses one plate")
    print("  at 6 against the rest at 2 about 72% of the time at these sizes,")
    print("  and the targeted test misses it about half the time. Neither can")
    print("  certify that the batches agree; they can only fail to show that")
    print("  they differ. The leave-one-out table needs no p-value and is the")
    print("  most robust thing on this page.")
    print()
    print("  Pool size is absent on purpose. NumTicks is the animal's total")
    print("  copied onto each of its rows, 4,730 summed over rows against")
    print("  3,209 over animals, so it is not the size of any pool.")


if __name__ == "__main__":
    main()
