#!/usr/bin/env python3
"""check_library_choice.py

I ask whether taking the deepest library per animal changes the core
microbiome numbers, by comparing it against drawing one library per animal at
random. I write nothing.

WHY
---
The Methods say the deepest library is kept for each animal. Depth is not
random: a library sequenced deeper may also have been handled differently, so
the rule could in principle select a biased subset. The honest check is to
draw at random instead and see whether the estimates move.

WHAT I COMPARE
--------------
For the same animals, under both rules:

    prevalence of Rickettsia, Anaplasma, Borreliella, Borrelia, Bartonella
    and Arsenophonus, as the share of animals with at least one read
    the median and mean Rickettsia share of bacterial reads

Only animals with more than one qualifying library can move, so I say how many
those are first. If almost every animal has one library, the rules agree by
construction and the check says little.

COHORT
------
The same as build_core_microbiome.py: deer with Ixodes and without
Dermacentor, libraries of at least 1,000 reads after host, chloroplast and
eukaryotic features are removed.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_library_choice.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DRAWS = 500
SEED = 20260924
GENERA = ("Rickettsia", "Anaplasma", "Borreliella", "Borrelia", "Bartonella",
          "Arsenophonus")


def suppression():
    try:
        import check_endosymbiont_suppression as module
    except ImportError as problem:
        sys.exit("cannot import check_endosymbiont_suppression (%s). It must "
                 "be in scripts/ beside me." % problem)
    return module


def statistics(genus_counts, animals, chosen):
    """I return the prevalence of each genus plus the Rickettsia share."""
    part = genus_counts.loc[chosen]
    result = {}
    for name in GENERA:
        if name in part.columns:
            result[name] = 100.0 * float((part[name] > 0).mean())
        else:
            result[name] = float("nan")
    shares = part.div(part.sum(axis=1), axis=0)
    rick = shares["Rickettsia"] if "Rickettsia" in shares.columns \
        else shares.iloc[:, 0] * 0
    result["Rickettsia median share"] = 100.0 * float(rick.median())
    result["Rickettsia mean share"] = 100.0 * float(rick.mean())
    return result


def main():
    module = suppression()
    arms = module.arms()
    table, lookup_module = module.lookup()
    counts = module.read_counts(table, lookup_module)
    path, covered, mapping = arms.read_taxonomy(list(counts.columns))
    print("counts: %d libraries x %d features" % counts.shape)
    print("taxonomy: %s (%d of %d)" % (path, covered, counts.shape[1]))
    host_like = [f for f in counts.columns
                 if any(word in str(mapping.get(f, "")).lower()
                        for word in module.HOST_WORDS)]
    counts = counts.drop(columns=host_like)
    genus = pd.Series({f: arms.genus_label(mapping.get(f, ""))
                       for f in counts.columns})
    genus_counts = counts.T.groupby(genus).sum().T

    animals = pd.Series([module.animal_of(table, lookup_module, s)
                         for s in counts.index], index=counts.index)
    depth = genus_counts.sum(axis=1)
    frame = pd.DataFrame({"animal": animals, "depth": depth})
    keep = set(module.deer_with_ixodes(table, lookup_module))
    frame = frame[(frame["animal"].isin(keep)) & (frame["depth"] >= module.MIN_READS)]
    if frame.empty:
        sys.exit("no deer library qualifies. Check the species table.")

    by_animal = {}
    for library, row in frame.iterrows():
        by_animal.setdefault(row["animal"], []).append(library)
    several = [a for a, v in by_animal.items() if len(v) > 1]
    print("deer with at least one qualifying library: %d" % len(by_animal))
    print("of those, with more than one: %d (the only ones that can move)"
          % len(several))
    print("qualifying libraries in total: %d" % len(frame))
    print()

    deepest = frame.sort_values("depth", ascending=False) \
                   .drop_duplicates("animal").index
    fixed = statistics(genus_counts, animals, list(deepest))

    rng = np.random.default_rng(SEED)
    drawn = []
    ordered = sorted(by_animal)
    for _ in range(DRAWS):
        picks = [by_animal[a][int(rng.integers(len(by_animal[a])))]
                 for a in ordered]
        drawn.append(statistics(genus_counts, animals, picks))
    spread = pd.DataFrame(drawn)

    print("=" * 76)
    print("DEEPEST LIBRARY AGAINST %d RANDOM DRAWS (seed %d)" % (DRAWS, SEED))
    print("=" * 76)
    print("  %-26s %10s %10s %18s" % ("statistic", "deepest", "median",
                                      "2.5 to 97.5 pct"))
    for name in list(GENERA) + ["Rickettsia median share",
                                "Rickettsia mean share"]:
        values = spread[name]
        low, high = np.percentile(values, [2.5, 97.5])
        flag = "" if low - 0.05 <= fixed[name] <= high + 0.05 else "   <-- outside"
        print("  %-26s %9.1f%% %9.1f%% %8.1f to %.1f%s"
              % (name, fixed[name], values.median(), low, high, flag))
    print()
    print("  A deepest-library value inside the interval means the rule is not")
    print("  selecting a biased set and the Methods sentence stands.")
    print()


if __name__ == "__main__":
    main()
