#!/usr/bin/env python3
"""I contrast the Ixodes-on-deer microbiome against three groupings.

WHY
---
Table 3 reports one microbiome for 330 deer. Three things recorded about
those tubes could divide it, and none has ever been tested:

  a deer ked in the tube, which the bench counted with real zeros in
  batches 07 to 10; whether the tube held one tick or a pool, which the
  paper already treats as two cohorts for qPCR; and the sex composition of
  the pool, which the accessioning notes give as engorged females and males
  for 96 animals.

The third carries a prediction. Rickettsia buchneri is maternally
inherited, so a pool that is mostly males should carry less of it. If that
holds it validates the notes as data; if it does not, that is worth knowing
before anyone leans on them.

WHAT I RESTRICT, AND WHY
------------------------
The ked comparison uses only animals whose ked record came from the bench,
because batches 01 to 06 recorded keds in prose and never recorded a zero,
so a blank there means unknown rather than none. Mixing the two regimes
would compare ascertainment rather than biology.

I report every restriction as it is applied, so the group sizes are
traceable rather than appearing from nowhere.

WHAT I DO NOT DO
----------------
I do not correct for multiple comparisons across the genera, and I say so
beside the numbers. Five genera against three groupings is fifteen tests
and one will reach p < 0.05 by chance about half the time. The sex
composition test is the only one specified in advance, and it is the only
one I would read on its own.

I write nothing. Everything goes to stdout.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import tick_lookup
import build_core_microbiome as core
import check_batch_arms as arms

COVARIATES = "results/metadata/accessioning_covariates.tsv"
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
WATCH = ("Rickettsia", "Bartonella", "Arsenophonus", "Anaplasma", "Borreliella")
DROP = ("Mitochondria", "Chloroplast", "Unassigned")


def stats():
    """I return the three tests, or None where scipy is absent."""
    try:
        from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
        return fisher_exact, mannwhitneyu, spearmanr
    except ImportError:
        return None, None, None


def compare(frame, mask_a, mask_b, label_a, label_b, title, mannwhitneyu,
            fisher_exact):
    """I compare genus prevalence and abundance between two animal groups."""
    a, b = frame[mask_a], frame[mask_b]
    print("=" * 76)
    print(title)
    print("=" * 76)
    print("  %s n = %d, %s n = %d" % (label_a, len(a), label_b, len(b)))
    if len(a) < 3 or len(b) < 3:
        print("  too few animals in one group to compare.\n")
        return
    print()
    print("  %-14s %-26s %-26s %9s %9s"
          % ("genus", label_a[:26], label_b[:26], "p abund", "p prev"))
    for name in WATCH:
        if name not in frame.columns:
            continue
        pa, pb = a[name], b[name]
        prev_a, prev_b = int((pa > 0).sum()), int((pb > 0).sum())
        p_abund = p_prev = float("nan")
        if mannwhitneyu is not None:
            try:
                p_abund = float(mannwhitneyu(pa, pb, alternative="two-sided")[1])
            except ValueError:
                p_abund = float("nan")
            p_prev = float(fisher_exact([[prev_a, len(a) - prev_a],
                                         [prev_b, len(b) - prev_b]])[1])
        print("  %-14s %2d/%-3d median %7.2f%%  %2d/%-3d median %7.2f%%  %9.3g %9.3g"
              % (name, prev_a, len(a), pa.median(),
                 prev_b, len(b), pb.median(), p_abund, p_prev))
    print()


def main():
    if not os.path.isfile(COVARIATES):
        sys.exit("missing %s. Run build_accessioning_covariates.py first."
                 % COVARIATES)
    fisher_exact, mannwhitneyu, spearmanr = stats()
    if mannwhitneyu is None:
        print("scipy is not importable here, so no tests are run; the counts "
              "and medians below still stand.\n")

    table, module = core.lookup()
    counts = core.read_counts(table, module)
    animals = [core.animal_of(table, module, s) for s in counts.index]
    _, _, mapping = arms.read_taxonomy(list(counts.columns))
    host_like = [f for f in counts.columns
                 if any(word in str(mapping.get(f, "")).lower()
                        for word in ("mitochondri", "chloroplast", "eukaryota"))]
    counts = counts.drop(columns=host_like)
    genus = pd.Series({f: arms.genus_label(mapping.get(f, ""))
                       for f in counts.columns})
    deer, moose, species = core.group_members(table, module)
    chosen = core.one_per_animal(counts, animals)
    deer_chosen = chosen[chosen["animal"].isin(deer)]
    print("cohort: %d deer with a library above the floor\n" % len(deer_chosen))

    part = counts.loc[deer_chosen.index]
    by_genus = part.T.groupby(genus).sum().T
    kept = by_genus.drop(columns=[c for c in by_genus.columns
                                  if c in DROP or c.endswith("(domain)")])
    frac = 100 * kept.div(kept.sum(axis=1), axis=0)
    frac["animal"] = deer_chosen["animal"].values
    frac = frac.set_index("animal")

    cov = pd.read_csv(COVARIATES, sep="\t", dtype=str)
    cov["animal"] = [tick_lookup.resolve(table, "specimen", a)
                     for a in cov["animal"]]
    cov = cov[cov["animal"].notna()].copy()
    cov["animal"] = cov["animal"].astype(str)
    cov = cov.drop_duplicates("animal").set_index("animal")
    shared = frac.index.intersection(cov.index)
    print("covariate rows %d, of which also in the cohort %d\n"
          % (len(cov), len(shared)))
    frame = frac.loc[shared].join(cov.loc[shared], rsuffix="_cov")

    # --- 1. ked, bench regime only ----------------------------------------
    bench = frame["ked_source"] == "bench"
    ked_yes = bench & (pd.to_numeric(frame["ked_count"], errors="coerce") > 0)
    ked_no = bench & (frame["ked_scored_zero"].astype(str).str.lower() == "true")
    print("ked: %d animals in the cohort have a bench ked record, %d with a "
          "ked and %d scored zero" % (int(bench.sum()), int(ked_yes.sum()),
                                      int(ked_no.sum())))
    compare(frame, ked_yes, ked_no, "ked present", "ked scored zero",
            "1. KED IN THE TUBE, BATCHES 07 TO 10 ONLY",
            mannwhitneyu, fisher_exact)

    # --- 2. single tick against pool ---------------------------------------
    cross = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    single = {}
    for _, row in cross.iterrows():
        animal = tick_lookup.resolve(table, "specimen", row["sample_id"])
        if animal is None:
            continue
        flag = str(row.get("is_single_tick")).strip().lower() in ("true", "1", "yes")
        ticks = pd.to_numeric(pd.Series([row.get("num_ticks")]),
                              errors="coerce").iloc[0]
        single[str(animal)] = (flag, ticks)
    frame["is_single"] = [single.get(a, (None, None))[0] for a in frame.index]
    frame["n_ticks"] = [single.get(a, (None, None))[1] for a in frame.index]
    is_single = frame["is_single"] == True            # noqa: E712
    is_pool = (frame["is_single"] == False) & (frame["n_ticks"] > 1)  # noqa: E712
    compare(frame, is_single, is_pool, "single tick", "pool",
            "2. ONE TICK AGAINST A POOL", mannwhitneyu, fisher_exact)

    # --- 3. sex composition, the specified test ----------------------------
    engorged = pd.to_numeric(frame["engorged_count"], errors="coerce")
    males = pd.to_numeric(frame["male_count"], errors="coerce")
    complete = frame["sex_census_complete"].astype(str).str.lower() == "true"
    total = engorged + males
    usable = complete & total.notna() & (total > 0)
    print("=" * 76)
    print("3. RICKETTSIA AGAINST THE FEMALE FRACTION, THE SPECIFIED TEST")
    print("=" * 76)
    print("  animals with a complete sex census and a library: %d"
          % int(usable.sum()))
    if int(usable.sum()) < 5 or "Rickettsia" not in frame.columns:
        print("  too few to test.\n")
        return
    female_fraction = (engorged[usable] / total[usable])
    rick = frame.loc[usable, "Rickettsia"]
    print("  female fraction: median %.2f, range %.2f to %.2f"
          % (female_fraction.median(), female_fraction.min(),
             female_fraction.max()))
    print("  Rickettsia share: median %.1f%%" % rick.median())
    if spearmanr is not None:
        rho, p = spearmanr(female_fraction, rick)
        print("  Spearman rho = %.3f, p = %.4g" % (rho, p))
        print()
        print("  Rickettsia buchneri is maternally inherited, so a positive")
        print("  rho is the direction predicted. This is the one test specified")
        print("  before looking; the genus tables above are not corrected for")
        print("  multiple comparison and a single result there is not evidence")
        print("  on its own.")
    print()


if __name__ == "__main__":
    main()
