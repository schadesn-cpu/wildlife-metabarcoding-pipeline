#!/usr/bin/env python3
"""check_rickettsia_qpcr.py

I ask whether the Rickettsia share of a tick's microbiome differs between
animals whose qPCR was positive and animals whose qPCR was negative. I write
nothing.

WHY THIS AND NOT check_endosymbiont_suppression.py
--------------------------------------------------
That script asks the same question with 16S pathogen reads as the outcome, and
has to spend three sections defending against the obvious objection: a library
that is 99% Rickettsia has 1% of its reads left for everything else, so a
pathogen at the same absolute level is harder to see. Detection and predictor
are measured on the same compositional scale.

qPCR does not have that problem. A qPCR call is made on the extraction, not on
the 16S library, so it cannot be depressed by Rickettsia crowding the reads.
That leaves compositionality only on the predictor side, which rarefaction
handles, and it makes the outcome the one the paper actually reports.

WHAT PREDICTS WHAT
------------------
Cull and colleagues (2022) found Rickettsia buchneri inhibited Anaplasma
phagocytophilum and pathogenic Rickettsia in tick cells. If that operates in
field ticks, animals positive for A. phagocytophilum should carry a lower
Rickettsia share. A. phagocytophilum is therefore the primary test, named
before looking. The other targets are secondary and are read as such: with
seven of them, one will reach p < 0.05 by chance about a third of the time.

COHORTS
-------
Single ticks are the clean comparison, because one tick gives one microbiome
and one qPCR call. Pools are reported second and are harder to read: a pool is
positive if any of its ticks is, so larger pools are more likely positive,
while the Rickettsia share is an average across all of them. I print the tick
count of both groups so an imbalance is visible rather than assumed away.

WHAT WOULD MAKE ME DISTRUST A DIFFERENCE
----------------------------------------
Sequencing depth or plate differing between the positive and negative groups.
I print both. A difference that travels with the plate is a plate effect.

WHAT I READ
-----------
    results/decontam/16S/table_decontam_16S.tsv, through
        check_endosymbiont_suppression
    metadata/tick_full_metadata.xlsx and the crosswalk, through
        check_qpcr_by_host

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_rickettsia_qpcr.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PRIMARY = "a-phagocytophilum"
SEED = 20260928
MIN_GROUP = 5


def suppression_module():
    """I return check_endosymbiont_suppression, whose cohort rules I reuse."""
    try:
        import check_endosymbiont_suppression as module
    except ImportError as problem:
        sys.exit("cannot import check_endosymbiont_suppression (%s). It must "
                 "sit in scripts/ beside me." % problem)
    for wanted in ("arms", "lookup", "animal_of", "read_counts",
                   "deer_with_ixodes", "rarefy", "permutation_test", "run_of",
                   "MIN_READS", "RAREFY_TO", "HOST_WORDS"):
        if not hasattr(module, wanted):
            sys.exit("check_endosymbiont_suppression has no %s; it has changed "
                     "since I was written and I will not guess." % wanted)
    return module


def qpcr_module():
    """I return check_qpcr_by_host, whose cohort rules I reuse."""
    try:
        import check_qpcr_by_host as module
    except ImportError as problem:
        sys.exit("cannot import check_qpcr_by_host (%s). It must sit in "
                 "scripts/ beside me." % problem)
    for wanted in ("read_qpcr", "read_crosswalk", "specimen_of", "TARGETS"):
        if not hasattr(module, wanted):
            sys.exit("check_qpcr_by_host has no %s; tell me what it calls that "
                     "and I will use it instead." % wanted)
    return module


def rickettsia_share(suppression, table, finder):
    """I return the rarefied Rickettsia share, depth and run, per deer animal.

    The cohort and the rarefaction are the ones check_endosymbiont_suppression
    already uses, so this is the same set of animals as the core microbiome,
    not a new one.
    """
    arms = suppression.arms()
    for wanted in ("read_taxonomy", "genus_label"):
        if not hasattr(arms, wanted):
            sys.exit("check_batch_arms has no %s, so I cannot label genera the "
                     "way the rest of the project does." % wanted)
    counts = suppression.read_counts(table, finder)
    path, covered, mapping = arms.read_taxonomy(list(counts.columns))
    print("16S table: %d libraries x %d features" % counts.shape)
    print("taxonomy:  %s (%d of %d features)" % (path, covered, counts.shape[1]))

    host_like = [f for f in counts.columns
                 if any(word in str(mapping.get(f, "")).lower()
                        for word in suppression.HOST_WORDS)]
    counts = counts.drop(columns=host_like)
    genus = pd.Series({f: arms.genus_label(mapping.get(f, ""))
                       for f in counts.columns})

    animals = pd.Series([suppression.animal_of(table, finder, s)
                         for s in counts.index], index=counts.index)
    depth = counts.sum(axis=1)
    frame = pd.DataFrame({"animal": animals, "depth": depth})
    frame = frame[(frame["animal"] != "")
                  & (frame["depth"] >= suppression.MIN_READS)]
    frame = frame.sort_values("depth", ascending=False).drop_duplicates("animal")
    deer = suppression.deer_with_ixodes(table, finder)
    chosen = frame[frame["animal"].isin(deer)]
    if chosen.empty:
        sys.exit("no deer library qualifies. Something is wrong upstream of me.")

    genus_counts = counts.loc[chosen.index].T.groupby(genus).sum().T
    rng = np.random.default_rng(SEED)
    rarefied = suppression.rarefy(genus_counts, suppression.RAREFY_TO, rng)
    if "Rickettsia" not in rarefied.columns:
        sys.exit("no Rickettsia column after collapsing to genus. The taxonomy "
                 "labels have changed and I will not guess which column is it.")
    shares = rarefied.div(rarefied.sum(axis=1), axis=0)

    kept = chosen.loc[rarefied.index]
    out = pd.DataFrame({
        "animal": kept["animal"].values,
        "share": shares["Rickettsia"].values,
        "depth": kept["depth"].values,
        "run": [suppression.run_of(s) for s in kept.index],
    })
    print("deer with Ixodes, one library each, rarefied to %d reads: %d"
          % (suppression.RAREFY_TO, len(out)))
    print()
    return out.set_index("animal")


def qpcr_calls(qpcr, table, finder):
    """I return one row per animal per cohort, with each target as 0 or 1."""
    counts, flags, hosts = qpcr.read_crosswalk()
    book = qpcr.read_qpcr()
    rows = []
    for _, row in book.iterrows():
        if str(row["subsample_count"]).strip() not in ("1", "1.0"):
            continue
        specimen = qpcr.specimen_of(table, finder, row)
        if not specimen:
            continue
        ticks = pd.to_numeric(pd.Series([counts.get(specimen)]),
                              errors="coerce").iloc[0]
        if flags.get(specimen):
            cohort = "single ticks"
        elif pd.notna(ticks) and ticks > 1:
            cohort = "pools"
        else:
            continue
        entry = {"animal": specimen, "cohort": cohort,
                 "ticks": float(ticks) if pd.notna(ticks) else 1.0}
        for target in qpcr.TARGETS:
            value = pd.to_numeric(pd.Series([row.get(target)]),
                                  errors="coerce").fillna(0).iloc[0]
            entry[target] = int(value > 0)
        rows.append(entry)
    frame = pd.DataFrame(rows)
    if frame.empty:
        sys.exit("no qPCR row resolved to a cohort. Check the sheet name and "
                 "the crosswalk.")
    # One animal, one result, the same rule check_qpcr_by_host applies.
    before = len(frame)
    frame = frame.drop_duplicates(subset=["cohort", "animal"], keep="first")
    if len(frame) < before:
        print("dropped %d repeated animal rows, keeping the first of each"
              % (before - len(frame)))
    return frame


def compare(shares, part, target, rng, suppression):
    """I return one printable row comparing positive against negative."""
    joined = part.join(shares, on="animal", how="inner")
    positive = joined[target] == 1
    if int(positive.sum()) < MIN_GROUP or int((~positive).sum()) < MIN_GROUP:
        return None
    values = joined["share"]
    _, p_mean = suppression.permutation_test(values, positive, rng)
    _, p_rank = suppression.permutation_test(values.rank(), positive, rng)
    return {
        "target": target,
        "pos": int(positive.sum()),
        "neg": int((~positive).sum()),
        "mean_pos": values[positive].mean(),
        "mean_neg": values[~positive].mean(),
        "median_pos": values[positive].median(),
        "median_neg": values[~positive].median(),
        "p_mean": p_mean,
        "p_rank": p_rank,
    }


def show(results, label):
    """I print one cohort's comparison table."""
    print("=" * 78)
    print(label.upper())
    print("=" * 78)
    if not results:
        print("  no target has at least %d animals on both sides, so this "
              "cohort cannot answer the question." % MIN_GROUP)
        print()
        return
    print("  %-18s %5s %5s %19s %19s %8s %8s"
          % ("target", "pos", "neg", "mean share pos/neg",
             "median share pos/neg", "p mean", "p rank"))
    for row in results:
        marker = " <- primary" if row["target"] == PRIMARY else ""
        print("  %-18s %5d %5d   %7.3f %7.3f   %7.3f %7.3f %8.4f %8.4f%s"
              % (row["target"], row["pos"], row["neg"],
                 row["mean_pos"], row["mean_neg"],
                 row["median_pos"], row["median_neg"],
                 row["p_mean"], row["p_rank"], marker))
    print()
    print("  A negative difference in the primary row, meaning positives carry")
    print("  less Rickettsia, is the direction the cell culture work predicts.")
    print()


def confounders(shares, part, target, label):
    """I print depth, plate and tick count for the two groups of one target."""
    joined = part.join(shares, on="animal", how="inner")
    positive = joined[target] == 1
    if int(positive.sum()) < MIN_GROUP or int((~positive).sum()) < MIN_GROUP:
        return
    print("-" * 78)
    print("WOULD I BELIEVE IT? %s, %s" % (target, label))
    print("-" * 78)
    print("  %-22s %14s %14s" % ("", "positive", "negative"))
    print("  %-22s %14d %14d"
          % ("animals", int(positive.sum()), int((~positive).sum())))
    print("  %-22s %14.0f %14.0f"
          % ("median 16S depth", joined.loc[positive, "depth"].median(),
             joined.loc[~positive, "depth"].median()))
    print("  %-22s %14.1f %14.1f"
          % ("median ticks in tube", joined.loc[positive, "ticks"].median(),
             joined.loc[~positive, "ticks"].median()))
    crossed = pd.crosstab(joined["run"], positive)
    crossed.columns = ["negative", "positive"][:crossed.shape[1]]
    print()
    print("  by sequencing run:")
    print("\n".join("    " + line for line in crossed.to_string().split("\n")))
    print()


def main():
    suppression = suppression_module()
    qpcr = qpcr_module()
    table, finder = suppression.lookup()

    shares = rickettsia_share(suppression, table, finder)
    calls = qpcr_calls(qpcr, table, finder)

    rng = np.random.default_rng(SEED)
    ordered = [PRIMARY] + [t for t in qpcr.TARGETS if t != PRIMARY]
    for cohort in ("single ticks", "pools"):
        part = calls[calls["cohort"] == cohort]
        overlap = part.join(shares, on="animal", how="inner")
        print("%s: %d animals screened, %d of them with a 16S library"
              % (cohort, len(part), len(overlap)))
        results = []
        for target in ordered:
            row = compare(shares, part, target, rng, suppression)
            if row is not None:
                results.append(row)
        show(results, "%s, Rickettsia share by qPCR result" % cohort)
        confounders(shares, part, PRIMARY, cohort)

    print("=" * 78)
    print("READING IT")
    print("=" * 78)
    print("  The primary test is %s. A secondary row reaching p < 0.05 on its"
          % PRIMARY)
    print("  own is not a finding: seven targets are tested and one will do")
    print("  that by chance about a third of the time.")
    print()
    print("  These are observed abundances, not manipulated ones, so neither")
    print("  a difference nor its absence establishes a causal direction.")
    print()


if __name__ == "__main__":
    main()
