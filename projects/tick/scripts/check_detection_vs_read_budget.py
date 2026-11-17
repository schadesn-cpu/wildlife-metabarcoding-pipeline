#!/usr/bin/env python3
"""I ask whether pathogen detection falls with Rickettsia because of biology
or because there are no reads left to find anything in.

WHY
---
The permutation test in check_endosymbiont_suppression.py rarefies every
animal to 1,000 reads. Where Rickettsia is 98% of classified reads that
leaves about 20 reads for every other organism combined, and a genus sitting
at a few percent of the non-Rickettsia fraction cannot be found in 20 reads
at all. So the rarefaction, which was added so that sequencing effort could
not manufacture an association, removes most heavily the reads the outcome
depends on, and removes them hardest in exactly the animals with the most
Rickettsia. A negative association between Rickettsia share and pathogen
detection is the expected result of that design whether or not anything
biological is happening.

The question I can ask instead is: given the same amount of looking at the
part of the library that is not Rickettsia, is a pathogen still less likely
to be found when the endosymbiont dominates?

HOW
---
I rarefy the non-Rickettsia reads rather than the whole library. Every
animal is subsampled to the same number of non-Rickettsia reads, so looking
effort is equalised where detection actually happens. I then repeat the same
permutation on the resulting detections.

I do this at several depths rather than one, because the choice of depth is
arbitrary and a single number hides whether the answer depends on it. An
association that is biology should persist, and sharpen, as depth rises and
detection improves. One that is a read-budget artifact should weaken once
every animal is given the same number of non-Rickettsia reads to be found
in.

I report the unrarefied result first, so the comparison has a baseline, and
I report the correlation between Rickettsia share and non-Rickettsia budget,
which is the mechanism of the artifact in one number.

WHAT I DO NOT CLAIM
-------------------
Equalising non-Rickettsia reads does not make this an absolute-abundance
measurement. An animal can have a small non-Rickettsia fraction because the
endosymbiont crowded the library or because there was little else in the
tick, and 16S alone cannot separate those. This tests whether the
association survives equal looking, which is narrower than whether the
endosymbiont suppresses anything.

Raising the depth also changes which animals are in the test, because an
animal is only usable at a depth it reaches. I print the n at every depth so
that a p value which moves because the cohort shrank is visible as such.

I write nothing. Everything goes to stdout.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "scripts")
import build_core_microbiome as core
import check_batch_arms as arms

PATHOGENS = ("Anaplasma", "Borreliella", "Borrelia")
DEPTHS = (100, 250, 500, 1000, 2000)
PERMUTATIONS = 9999
SEED = 20260924
DROP = ("Mitochondria", "Chloroplast", "Unassigned")


def permutation_p(share, detected, rng):
    """I permute the detection labels and compare mean Rickettsia share.

    This is the test check_endosymbiont_suppression.py runs, on a different
    detection call, so the two are directly comparable.
    """
    share = np.asarray(share, dtype=float)
    detected = np.asarray(detected, dtype=bool)
    if detected.sum() < 2 or (~detected).sum() < 2:
        return float("nan"), float("nan"), float("nan")
    observed = share[detected].mean() - share[~detected].mean()
    labels = detected.copy()
    bigger = 0
    for _ in range(PERMUTATIONS):
        rng.shuffle(labels)
        if abs(share[labels].mean() - share[~labels].mean()) >= abs(observed):
            bigger += 1
    return (share[detected].mean(), share[~detected].mean(),
            (bigger + 1.0) / (PERMUTATIONS + 1.0))


def rarefy(matrix, depth, rng):
    """I subsample each row to `depth` reads without replacement.

    My caller guarantees every row holds at least `depth` reads. If one does
    not I stop, because silently returning a shallower row would put animals
    with different looking effort into a test whose whole purpose is to
    equalise it.
    """
    out = np.zeros(matrix.shape, dtype=np.int64)
    for i in range(matrix.shape[0]):
        row = matrix[i].astype(np.int64)
        total = int(row.sum())
        if total < depth:
            sys.exit("row %d holds %d reads and I was asked for %d."
                     % (i, total, depth))
        if total == depth:
            out[i] = row
            continue
        out[i] = rng.multivariate_hypergeometric(row, depth)
    return out


def main():
    if not os.path.isdir("results/analysis"):
        sys.exit("run me from the project root; there is no results/analysis here.")
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
    part = counts.loc[deer_chosen.index]
    by_genus = part.T.groupby(genus).sum().T
    kept = by_genus.drop(columns=[c for c in by_genus.columns
                                  if c in DROP or c.endswith("(domain)")])
    if "Rickettsia" not in kept.columns:
        sys.exit("no Rickettsia column after collapsing to genus; there is "
                 "nothing to test against.")
    share = (kept["Rickettsia"] / kept.sum(axis=1)).values
    nonrick = kept.drop(columns=["Rickettsia"])
    budget = nonrick.sum(axis=1).values.astype(np.int64)

    print("cohort %d deer carrying Ixodes and no Dermacentor" % len(deer_chosen))
    print()
    print("THE BUDGET DETECTION ACTUALLY HAS")
    print("  non-Rickettsia classified reads per animal")
    for q in (0, 5, 10, 25, 50, 75, 100):
        print("     %3dth percentile %10d" % (q, int(np.percentile(budget, q))))
    print("  animals left with fewer than 1,000 reads once Rickettsia is")
    print("  removed, though the Section 2.6 test rarefies to 1,000: %d of %d"
          % (int((budget < 1000).sum()), len(budget)))
    try:
        from scipy.stats import spearmanr
        rho, p = spearmanr(share, budget)
        print("  Rickettsia share against that budget: Spearman rho = %.3f, "
              "p = %.3g" % (rho, p))
        print("  A negative rho is the artifact itself: the animals the test")
        print("  calls endosymbiont-dominated are the animals it gives the")
        print("  fewest reads to find a pathogen in.")
    except ImportError:
        print("  scipy is not importable, so no correlation was computed.")
    print()

    print("=" * 76)
    print("A. UNRAREFIED: detection on every read the library holds")
    print("=" * 76)
    rng = np.random.default_rng(SEED)
    print("  %-13s %9s %13s %13s %10s"
          % ("pathogen", "positive", "Rick if pos", "Rick if neg", "p"))
    for name in PATHOGENS:
        if name not in nonrick.columns:
            print("  %-13s not present in this cohort" % name)
            continue
        detected = nonrick[name].values > 0
        pos, neg, pval = permutation_p(share, detected, rng)
        print("  %-13s %9d %13.3f %13.3f %10.4f"
              % (name, int(detected.sum()), pos, neg, pval))
    print()

    print("=" * 76)
    print("B. NON-RICKETTSIA READS EQUALISED, AT SEVERAL DEPTHS")
    print("=" * 76)
    matrix = nonrick.values
    columns = list(nonrick.columns)
    for depth in DEPTHS:
        usable = budget >= depth
        if int(usable.sum()) < 20:
            print("  depth %-5d only %d animals reach it; skipped"
                  % (depth, int(usable.sum())))
            print()
            continue
        rng = np.random.default_rng(SEED)
        sub = pd.DataFrame(rarefy(matrix[usable], depth, rng), columns=columns)
        sub_share = share[usable]
        richness = (sub > 0).sum(axis=1)
        print("  depth %-5d animals %3d of %d, Rickettsia share median %.3f, "
              "genera seen median %d"
              % (depth, int(usable.sum()), len(budget),
                 float(np.median(sub_share)), int(richness.median())))
        for name in PATHOGENS:
            if name not in sub.columns:
                continue
            detected = sub[name].values > 0
            pos, neg, pval = permutation_p(sub_share, detected, rng)
            if np.isnan(pval):
                print("     %-13s positive %4d, too few to permute"
                      % (name, int(detected.sum())))
                continue
            print("     %-13s positive %4d   Rickettsia %.3f against %.3f   "
                  "p = %.4f"
                  % (name, int(detected.sum()), pos, neg, pval))
        print()

    print("HOW TO READ THIS")
    print("  Compare each pathogen's p down the depths, not across pathogens.")
    print("  Persisting and sharpening as depth rises is what biology looks")
    print("  like. Weakening as the budget is equalised is what the artifact")
    print("  looks like. A p that moves while n drops sharply is the cohort")
    print("  changing, not the answer.")
    print("  %d permutations, seed %d, matching Section 2.6."
          % (PERMUTATIONS, SEED))


if __name__ == "__main__":
    main()
