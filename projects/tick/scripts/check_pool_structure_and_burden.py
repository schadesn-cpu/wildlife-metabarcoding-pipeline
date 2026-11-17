#!/usr/bin/env python3
"""check_pool_structure_and_burden.py

I ask two things the analysis has never asked. Whether an animal's several
libraries are ALIQUOTS of one pooled extract or the animal's ticks DIVIDED
between tubes, which decides what a pool size means. And whether the pooled
co-occurrence survives conditioning on how many ticks the animal carried,
which is the confounder we have never tested because we were conditioning on
the wrong unit. I write nothing and change nothing.

WHERE THIS SITS
---------------
The design is a pooled design in which the POOL IS THE ANIMAL: 413 animals,
pool size NumTicks, and the 602 rows of analysis_table_16S.tsv are repeat
measurements of those 413 pools. Everything before today treated the 602 rows
as 602 pools.

At the animal level the association is real and survives the unit:

    Anaplasma x Borreliella    1.96 resampled over 413 animals (1.68 to 2.29)
                               2.00 on the 311 single-library animals

but it is not the same in every plate. AMP-120924-USGS-TICKS-P3-5 runs at
6.80 against about 1.5 everywhere else, and Breslow-Day rejects homogeneity in
85% of draws even after each animal is counted once.

TICK BURDEN IS THE OBVIOUS CONFOUNDER AND WE NEVER TESTED IT
------------------------------------------------------------
An animal carrying 54 ticks has many more chances to hold at least one
Anaplasma-positive tick AND at least one Borreliella-positive tick than an
animal carrying one. That produces a positive association between any two
reasonably prevalent organisms with no biology in it whatsoever.

This is the same argument as the pool-size conditioning from this morning. It
failed then because it was applied to 602 rows in which NumTicks was the
animal's total repeated, so the strata were not pool sizes. On 413 animals it
is the correct variable, and P3-5 spans 1 to 54 ticks per animal.

ALIQUOTS OR DIVIDED? A TEST THAT DOES NOT NEED THE RUN SHEET
------------------------------------------------------------
EQ02441291 has 11 ticks and 10 subsample libraries. Either the 11 ticks went
into one extract that was split ten ways, in which case every row is a real
11-tick pool and the rows are technical replicates; or the ticks were divided,
in which case each row is about one tick and NumTicks is wrong on every row.

If they were DIVIDED, a library from a ten-library animal carries a tenth of
the ticks, so it should be positive far less often than a single-library
animal carrying the same number of ticks. A pool of n ticks is positive with
probability 1-(1-p)^n, so at p = 0.10 an 11-tick pool is positive 69% of the
time and a 1-tick pool 10% of the time. That is a large, visible difference.

If they were ALIQUOTS, positivity should be about the same in both.

So I compare per-row positivity between single-library and many-library
animals, holding NumTicks constant. That does not need Kellen to answer, and
it will not be ambiguous.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_pool_structure_and_burden.py
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
MIN_STRATUM = 15
MIN_POSITIVES = 30
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
    row_one, col_one, total = a + b, a + c, a + b + c + d
    if row_one in (0, total) or col_one in (0, total):
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
    """I return (statistic, df, p) for Tarone's corrected homogeneity test."""
    statistic = residual_sum = variance_sum = 0.0
    used = 0
    for a, b, c, d in tables:
        expected = expected_under_common_or(a, b, c, d, psi)
        if expected is None:
            continue
        row_one, col_one, total = a + b, a + c, a + b + c + d
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
    return statistic, used - 1, float(stats.chi2.sf(statistic, used - 1))


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


def token_report(frame):
    """I use the run token as the plate, and show where RunBatch disagrees."""
    print("=" * 72)
    print("1. THE RUN TOKEN AS THE PLATE")
    print("=" * 72)
    print("  RunBatch is missing on some rows and wrong on others. The token at")
    print("  the end of every sample id is present on all of them, so I")
    print("  stratify on that and use RunBatch only to check it.")
    print()
    disagreements = 0
    for token, group in frame.groupby("_token"):
        labels = sorted(set(b for b in group["_batch"] if b))
        missing = int((group["_batch"] == "").sum())
        note = ""
        if len(labels) > 1:
            note = "  <- %d DIFFERENT RunBatch strings" % len(labels)
            disagreements += 1
        print("  %-18s %4d rows, %4d animals, %d RunBatch string(s),"
              " %d blank%s"
              % (str(token)[:18], len(group), group["_animal"].nunique(),
                 len(labels), missing, note))
        if len(labels) > 1:
            counts = group["_batch"].replace("", "(blank)").value_counts()
            for name, n in counts.items():
                print("        %-52s %4d" % (str(name)[:52], n))
            minority = counts[counts.index != ""].idxmin()
            rows = group[group["_batch"] == minority]
            print("        the minority rows, to chase in the run sheet:")
            for name in sorted(str(s) for s in rows["sample"]):
                print("          " + name)
    if not disagreements:
        print("  no token carries more than one RunBatch string.")
    print()


def plate_profile(frame, present):
    """I profile each plate: rows per animal, host mix, tick burden."""
    print("=" * 72)
    print("2. WHAT IS DIFFERENT ABOUT EACH PLATE")
    print("=" * 72)
    host_column = next((c for c in ("Host_group", "Host_clean", "Host")
                        if c in frame.columns), None)
    animals = frame.groupby(["_token", "_animal"]).first().reset_index()
    print("  %-18s %7s %7s %7s %-22s %s"
          % ("plate", "rows", "animals", "r/a", "ticks per animal",
             "hosts" if host_column else ""))
    for token, group in animals.groupby("_token"):
        rows = int((frame["_token"] == token).sum())
        ticks = pd.to_numeric(group["NumTicks"], errors="coerce").dropna()
        burden = ("med %g, max %g" % (ticks.median(), ticks.max())
                  if len(ticks) else "n/a")
        hosts = ""
        if host_column:
            counts = group[host_column].fillna("(blank)").value_counts()
            hosts = ", ".join("%s %d" % (str(name)[:12], n)
                              for name, n in list(counts.items())[:4])
        print("  %-18s %7d %7d %7.2f %-22s %s"
              % (str(token)[:18], rows, len(group), rows / float(len(group)),
                 burden, hosts))
    print()
    print("  A plate that mixes host species is mixing populations, and mixing")
    print("  two populations with different prevalences of BOTH organisms")
    print("  produces a positive association on its own.")
    print()


def fit_per_tick_prevalence(ticks, positives):
    """I fit the per-tick prevalence p from pools of known size, by grid MLE.

    A pool of n ticks is positive with probability 1-(1-p)^n. Given pool sizes
    and outcomes, that is a one-parameter likelihood, so a grid is plenty and
    it cannot fail to converge the way an optimiser can.
    """
    ticks = np.asarray(ticks, dtype=float)
    positives = np.asarray(positives, dtype=bool)
    keep = np.isfinite(ticks) & (ticks > 0)
    ticks, positives = ticks[keep], positives[keep]
    if not len(ticks) or not positives.any():
        return None
    grid = np.linspace(1e-5, 0.95, 4000)
    best, best_p = -np.inf, None
    for candidate in grid:
        miss = (1.0 - candidate) ** ticks
        hit = np.clip(1.0 - miss, 1e-12, 1.0)
        total = float(np.sum(np.log(hit[positives]))
                      + np.sum(np.log(np.clip(miss[~positives], 1e-12, 1.0))))
        if total > best:
            best, best_p = total, candidate
    return best_p


def aliquot_or_divided(frame, present):
    """I test whether extra libraries are aliquots of one pool or split ticks.

    I fit the per-tick prevalence on the SINGLE-library animals, where the pool
    size is certainly the animal's tick count and there is nothing to argue
    about. Then for the multi-library rows I predict positivity two ways:

        aliquots  every library saw the whole pool   1-(1-p)^n
        divided   each library saw n/L of the ticks  1-(1-p)^(n/L)

    and compare both predictions against what those rows actually do. This
    turns the question into two numbers to sit between rather than a column to
    squint at.
    """
    print("=" * 72)
    print("3. ALIQUOTS OF ONE POOL, OR THE ANIMAL'S TICKS DIVIDED?")
    print("=" * 72)
    ticks = pd.to_numeric(frame["NumTicks"], errors="coerce")
    work = frame.assign(_ticks=ticks)
    work = work[(work["_animal"] != "") & work["_ticks"].notna()]
    work = work.assign(
        _libraries=work.groupby("_animal")["sample"].transform("size"))
    single = work[work["_libraries"] == 1]
    many = work[work["_libraries"] >= 2]
    print("  fitting per-tick prevalence on %d single-library animals, then"
          % len(single))
    print("  predicting the %d rows from %d multi-library animals two ways."
          % (len(many), many["_animal"].nunique()))
    print()
    if len(single) < 30 or len(many) < 20:
        print("  not enough of one group to compare. Nothing to conclude.")
        print()
        return
    print("  %-14s %8s %11s %11s %11s %s"
          % ("", "p/tick", "observed", "if aliquot", "if divided", "verdict"))
    for pathogen in present:
        column = pathogen + "_pos"
        rate = fit_per_tick_prevalence(single["_ticks"],
                                       single[column] == 1)
        if rate is None:
            print("  %-14s could not be fitted" % pathogen)
            continue
        sizes = many["_ticks"].values.astype(float)
        libraries = many["_libraries"].values.astype(float)
        as_aliquot = float(np.mean(1.0 - (1.0 - rate) ** sizes))
        as_divided = float(np.mean(
            1.0 - (1.0 - rate) ** np.maximum(sizes / libraries, 1.0)))
        observed = float((many[column] == 1).mean())
        if abs(as_aliquot - as_divided) < 0.02:
            verdict = "cannot tell, predictions too close"
        elif abs(observed - as_aliquot) < abs(observed - as_divided):
            verdict = "closer to ALIQUOT"
        else:
            verdict = "closer to DIVIDED"
        print("  %-14s %8.3f %10.1f%% %10.1f%% %10.1f%% %s"
              % (pathogen, rate, 100 * observed, 100 * as_aliquot,
                 100 * as_divided, verdict))
    print()
    print("  The fit assumes every tick is independent and equally likely to")
    print("  carry the organism, which is the same assumption PoolTestR makes.")
    print("  It is not true, and it does not need to be: both predictions use")
    print("  it, so the comparison between them survives it even where the")
    print("  absolute numbers do not.")
    print()
    buckets = [(1, 1, "1"), (2, 3, "2-3"), (4, 6, "4-6"), (7, 10, "7-10"),
               (11, 20, "11-20"), (21, 10 ** 6, "21+")]
    print("  and the same thing without any model, per-row positivity by the")
    print("  animal's tick count:")
    for pathogen in present:
        print("    %s:" % pathogen)
        print("      %-8s %-22s %-22s %s"
              % ("ticks", "1-library", "3+-library", "ratio"))
        for low, high, label in buckets:
            band = work[(work["_ticks"] >= low) & (work["_ticks"] <= high)]
            one = band[band["_libraries"] == 1]
            lots = band[band["_libraries"] >= 3]
            if len(one) < 10 or len(lots) < 10:
                print("      %-8s %-22s n=%d, too few"
                      % (label, "n=%d" % len(one), len(lots)))
                continue
            one_rate = float((one[pathogen + "_pos"] == 1).mean())
            lots_rate = float((lots[pathogen + "_pos"] == 1).mean())
            print("      %-8s %-22s %-22s %s"
                  % (label,
                     "%.1f%% (%d/%d)" % (100 * one_rate,
                                         int((one[pathogen + "_pos"] == 1).sum()),
                                         len(one)),
                     "%.1f%% (%d/%d)" % (100 * lots_rate,
                                         int((lots[pathogen + "_pos"] == 1).sum()),
                                         len(lots)),
                     "%.2f" % (lots_rate / one_rate) if one_rate > 0 else "n/a"))
        print()
    print("  Ratios near 1 mean ALIQUOTS: every library saw the whole pool, the")
    print("  rows are technical replicates, and NumTicks is the pool size on")
    print("  every row.")
    print("  Ratios well below 1 mean DIVIDED: each library saw a fraction of")
    print("  the ticks, NumTicks is wrong on every row of a multi-library")
    print("  animal, and much of P3-5 is closer to single ticks than to pools.")
    print("  The run sheet still decides it. This says which answer the data")
    print("  is consistent with.")
    print()


def burden_report(frame, present):
    """I condition the animal-level co-occurrence on the animal's tick count."""
    print("=" * 72)
    print("4. CO-OCCURRENCE CONDITIONED ON TICK BURDEN, ONE ROW PER ANIMAL")
    print("=" * 72)
    work = frame[frame["_animal"] != ""].reset_index(drop=True)
    ticks = pd.to_numeric(work["NumTicks"], errors="coerce")
    animal_ticks = work.assign(_ticks=ticks).groupby("_animal")["_ticks"].max()
    try:
        bands = pd.qcut(animal_ticks, 4, labels=["Q1", "Q2", "Q3", "Q4"],
                        duplicates="drop")
    except ValueError:
        print("  tick counts cannot be split into quartiles.")
        print()
        return
    band_of = bands.to_dict()
    ranges = {}
    for level in bands.cat.categories:
        values = animal_ticks[bands == level]
        ranges[level] = (values.min(), values.max(), len(values))

    groups, names = [], []
    for name, positions in work.groupby("_animal").indices.items():
        groups.append(np.asarray(positions))
        names.append(name)
    lengths = np.array([len(g) for g in groups])
    padded = np.zeros((len(groups), int(lengths.max())), dtype=np.int64)
    for i, g in enumerate(groups):
        padded[i, :len(g)] = g
    rows = np.arange(len(groups))
    band_codes = np.array([list(bands.cat.categories).index(band_of[n])
                           for n in names])
    suspect = np.array([SUSPECT_TOKEN in set(work["_token"].iloc[g])
                        for g in groups])
    positives = {p: (work[p + "_pos"].values == 1) for p in present}
    generator = np.random.default_rng(SEED)

    pairs = list(itertools.combinations(present, 2))
    collected = {pair: {"crude": [], "mh": [], "bd": [],
                        "strata": {level: [] for level in ranges},
                        "suspect": [], "rest": [], "rest_mh": []}
                 for pair in pairs}

    for _ in range(DRAWS):
        choice = (generator.random(len(groups)) * lengths).astype(np.int64)
        picked = padded[rows, choice]
        for pair in pairs:
            first, second = pair
            x = positives[first][picked]
            y = positives[second][picked]
            store = collected[pair]
            store["crude"].append(odds_ratio(
                int(np.count_nonzero(x & y)), int(np.count_nonzero(x & ~y)),
                int(np.count_nonzero(~x & y)), int(np.count_nonzero(~x & ~y))))
            tables = []
            for code, level in enumerate(ranges):
                mask = band_codes == code
                if int(np.count_nonzero(mask)) < MIN_STRATUM:
                    continue
                xs, ys = x[mask], y[mask]
                table = (int(np.count_nonzero(xs & ys)),
                         int(np.count_nonzero(xs & ~ys)),
                         int(np.count_nonzero(~xs & ys)),
                         int(np.count_nonzero(~xs & ~ys)))
                tables.append(table)
                store["strata"][level].append(odds_ratio(*table))
            summary = mantel_haenszel(tables)
            if summary is not None:
                store["mh"].append(summary)
                test = breslow_day(tables, summary)
                if test is not None:
                    store["bd"].append(test[2])
            for label, mask in (("suspect", suspect), ("rest", ~suspect)):
                if int(np.count_nonzero(mask)) < MIN_STRATUM:
                    continue
                xs, ys = x[mask], y[mask]
                store[label].append(odds_ratio(
                    int(np.count_nonzero(xs & ys)),
                    int(np.count_nonzero(xs & ~ys)),
                    int(np.count_nonzero(~xs & ys)),
                    int(np.count_nonzero(~xs & ~ys))))
            inner = []
            for code, level in enumerate(ranges):
                mask = (~suspect) & (band_codes == code)
                if int(np.count_nonzero(mask)) < MIN_STRATUM:
                    continue
                xs, ys = x[mask], y[mask]
                inner.append((int(np.count_nonzero(xs & ys)),
                              int(np.count_nonzero(xs & ~ys)),
                              int(np.count_nonzero(~xs & ys)),
                              int(np.count_nonzero(~xs & ~ys))))
            rest_mh = mantel_haenszel(inner)
            if rest_mh is not None:
                store["rest_mh"].append(rest_mh)

    observed_crude = {}
    for pair in pairs:
        values = [v for v in collected[pair]["crude"] if np.isfinite(v)]
        if values:
            observed_crude[pair] = float(np.median(values))
    animal_ticks_array = [animal_ticks[name] for name in names]

    def band_of_values(values):
        arr = np.array([v for v in values if np.isfinite(v)])
        if not len(arr):
            return None
        return (float(np.median(arr)), float(np.percentile(arr, 2.5)),
                float(np.percentile(arr, 97.5)))

    print("  %d animals, %d draws, seed %d. Strata are quartiles of the"
          % (len(groups), DRAWS, SEED))
    print("  animal's own tick count:")
    for level, (low, high, n) in ranges.items():
        print("    %-4s %g to %g ticks, %d animals" % (level, low, high, n))
    print()
    for pair in pairs:
        store = collected[pair]
        crude = band_of_values(store["crude"])
        summary = band_of_values(store["mh"])
        print("  %s x %s" % pair)
        print("    crude over animals:            %.2f  (%.2f to %.2f)" % crude)
        if summary:
            print("    Mantel-Haenszel over burden:   %.2f  (%.2f to %.2f)"
                  % summary)
        for level in ranges:
            values = band_of_values(store["strata"][level])
            if values:
                print("      %-4s %g-%g ticks   OR %.2f  (%.2f to %.2f)"
                      % (level, ranges[level][0], ranges[level][1],
                         values[0], values[1], values[2]))
        if store["bd"]:
            arr = np.array(store["bd"])
            print("    homogeneity across burden: median p %.4g,"
                  " rejects in %.0f%% of draws"
                  % (float(np.median(arr)), 100.0 * (arr < 0.05).mean()))
        inside = band_of_values(store["suspect"])
        outside = band_of_values(store["rest"])
        outside_adjusted = band_of_values(store["rest_mh"])
        print("    on %s:                %s"
              % (SUSPECT_TOKEN,
                 "%.2f  (%.2f to %.2f)" % inside if inside else "too few"))
        print("    everywhere else:               %s"
              % ("%.2f  (%.2f to %.2f)" % outside if outside else "too few"))
        print("    everywhere else, burden-adjusted: %s"
              % ("%.2f  (%.2f to %.2f)" % outside_adjusted
                 if outside_adjusted else "too few"))
        print()
    print("  A quartile-adjusted number still carries most of the confounding:")
    print("  on simulated animals whose organisms are independent tick by tick,")
    print("  the crude came out 3.32 and this adjustment left 1.50. Section 5")
    print("  is the one that settles it.")
    print()
    return observed_crude, animal_ticks_array


def burden_null_report(frame, present, animal_status, animal_ticks_array):
    """I ask what odds ratio TICK BURDEN ALONE would produce, and compare.

    WHY STRATIFYING IS NOT ENOUGH
    -----------------------------
    Checked by simulation before this was written: on animals whose two
    organisms are independent tick by tick, so the only mechanism is burden,
    the crude odds ratio came out 3.32 and the Mantel-Haenszel over burden
    QUARTILES came out 1.50. Quartiles are too coarse. Inside a stratum running
    from 17 to 54 ticks the burden still varies four-fold, and that residual
    confounding is most of what a quartile-adjusted number leaves behind.

    So a quartile-adjusted 1.45 is not evidence of anything. It is roughly what
    burden alone gives.

    WHAT I DO INSTEAD
    -----------------
    I fit each organism's per-tick prevalence on the single-library animals,
    then simulate every animal's two outcomes INDEPENDENTLY, each with
    probability 1-(1-p)^n at that animal's own tick count. That is the null in
    which burden is the only thing linking the two organisms. The odds ratio
    it produces is the number to beat, and I report where the observed one
    falls in that distribution.
    """
    print("=" * 72)
    print("5. WHAT ODDS RATIO WOULD TICK BURDEN ALONE PRODUCE?")
    print("=" * 72)
    ticks = pd.to_numeric(frame["NumTicks"], errors="coerce")
    work = frame.assign(_ticks=ticks)
    work = work[(work["_animal"] != "") & work["_ticks"].notna()]
    work = work.assign(
        _libraries=work.groupby("_animal")["sample"].transform("size"))
    single = work[work["_libraries"] == 1]
    rates = {}
    for pathogen in present:
        rate = fit_per_tick_prevalence(single["_ticks"],
                                       single[pathogen + "_pos"] == 1)
        if rate is not None:
            rates[pathogen] = rate
    if len(rates) < 2:
        print("  could not fit per-tick prevalence for two organisms.")
        print()
        return
    print("  per-tick prevalence fitted on %d single-library animals:"
          % len(single))
    for pathogen, rate in rates.items():
        print("    %-14s %.4f" % (pathogen, rate))
    print()
    sizes = np.asarray(animal_ticks_array, dtype=float)
    sizes = np.where(np.isfinite(sizes) & (sizes > 0), sizes, 1.0)
    generator = np.random.default_rng(SEED + 1)
    print("  %-26s %10s %10s %12s %s"
          % ("pair", "observed", "burden null", "null 97.5th", "verdict"))
    for first, second in itertools.combinations(
            [p for p in present if p in rates], 2):
        observed = animal_status.get((first, second))
        if observed is None:
            continue
        p_first = 1.0 - (1.0 - rates[first]) ** sizes
        p_second = 1.0 - (1.0 - rates[second]) ** sizes
        null = []
        for _ in range(DRAWS):
            x = generator.random(len(sizes)) < p_first
            y = generator.random(len(sizes)) < p_second
            null.append(odds_ratio(
                int(np.count_nonzero(x & y)), int(np.count_nonzero(x & ~y)),
                int(np.count_nonzero(~x & y)), int(np.count_nonzero(~x & ~y))))
        null = np.array(null)
        share = float((null >= observed).mean())
        if share > 0.05:
            verdict = "within burden, p = %.3f" % share
        else:
            verdict = "ABOVE burden, p = %.3f" % share
        print("  %-26s %10.2f %10.2f %12.2f %s"
              % ("%s x %s" % (first, second), observed,
                 float(np.median(null)), float(np.percentile(null, 97.5)),
                 verdict))
    print()
    print("  'burden null' is the odds ratio these animals would show if the")
    print("  two organisms were independent in every tick and the only thing")
    print("  connecting them were how many ticks each animal carried.")
    print()
    print("  An observed value inside that distribution is not evidence of")
    print("  co-occurrence at all. One above it is the association that tick")
    print("  burden does not explain.")
    print()
    print("  This assumes ticks within an animal are independent, which they")
    print("  are not: ticks off one deer share its exposure. Real clustering")
    print("  would make the true null WIDER than this one, so a result that")
    print("  only just clears this bar has not cleared much.")
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
    frame["_batch"] = ([blank(v) for v in frame[BATCH_COLUMN]]
                       if BATCH_COLUMN in frame.columns else "")
    tokens = []
    for sample in frame["sample"]:
        found = RUN_TOKEN.search(str(sample).strip())
        tokens.append(found.group(1) if found else "(no token)")
    frame["_token"] = tokens
    labels, scheme = resolve_animals(frame)
    frame["_animal"] = labels

    present = [p for p in PATHOGENS + CO_DETECTIONS
               if p + "_pos" in frame.columns and frame[p + "_pos"].sum() > 0]
    thick = [p for p in present
             if int((frame[p + "_pos"] == 1).sum()) >= MIN_POSITIVES]

    print("rows: %d" % len(frame))
    print("animals: %d, identified by the %s"
          % (frame.loc[frame["_animal"] != "", "_animal"].nunique(), scheme))
    print()

    token_report(frame)
    plate_profile(frame, present)
    aliquot_or_divided(frame, present)
    if len(thick) >= 2:
        result = burden_report(frame, thick)
        if result is not None:
            burden_null_report(frame, thick, result[0], result[1])
    else:
        print("fewer than two organisms with %d positives; skipping section 4."
              % MIN_POSITIVES)

    print("=" * 72)
    print("HOW TO READ THIS")
    print("=" * 72)
    print("  Section 3 decides what a row is. Until that is settled, section 4")
    print("  assumes NumTicks is the animal's pool size, which is right if the")
    print("  libraries are aliquots and wrong if the ticks were divided.")
    print()
    print("  In section 4, if the burden-adjusted odds ratio falls toward 1")
    print("  while the crude one holds near 2, the association is animals")
    print("  carrying more ticks having more chances to show both organisms,")
    print("  and it is not co-infection in any sense.")
    print()
    print("  If it holds after burden AND the plate difference survives, then")
    print("  P3-5 differs for a reason burden does not explain, and the")
    print("  reportable number is the one from the other plates.")
    print()
    print("  Nothing here can separate a batch effect from an animal effect:")
    print("  530 of 551 animals sit entirely within one batch by design.")


if __name__ == "__main__":
    main()
