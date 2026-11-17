#!/usr/bin/env python3
"""check_endosymbiont_suppression.py

I ask whether pathogen detection in Ixodes from white-tailed deer falls as the
Rickettsia endosymbiont takes up more of the microbiome. I write nothing.

WHY
---
Rickettsia buchneri, the Ixodes scapularis endosymbiont, inhibits Anaplasma
phagocytophilum and pathogenic Rickettsia by three to five orders of magnitude
in tick cells (Kurtti and colleagues, 2021). If that happens in real ticks, deer
whose ticks are nearly all Rickettsia should carry fewer pathogen reads.

THE ALTERNATIVE EXPLANATION, WHICH I TEST RATHER THAN IGNORE
------------------------------------------------------------
Relative abundance is compositional. A library that is 99% Rickettsia has only
1% of its reads left for everything else, so a pathogen present at the same
absolute level is harder to detect. That alone would produce the pattern.

So I report three things, and only the third is evidence:

  1. detection by Rickettsia band, with the median depth of each band, so a
     depth difference between bands is visible
  2. the same after rarefying every animal to a common depth, which removes
     the compositional advantage of a shallow library
  3. a permutation test on the rarefied data

I also print the batch composition of each band, because the low-Rickettsia
deer concentrate on the P3-5 plate and a plate effect would masquerade as a
biological one.

COHORT
------
The same one as build_core_microbiome.py: deer with Ixodes and without
Dermacentor, one library per animal, the deepest, with at least 1,000 reads.

WHAT I READ
-----------
    results/decontam/16S/table_decontam_16S.tsv
    the 16S taxonomy, through check_batch_arms
    results/analysis/tick_species_by_host.tsv
    scripts/tick_lookup.py

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_endosymbiont_suppression.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

COUNTS = "results/decontam/16S/table_decontam_16S.tsv"
SPECIES = "results/analysis/tick_species_by_host.tsv"
MIN_READS = 1000
RAREFY_TO = 1000
PATHOGENS = ("Anaplasma", "Borreliella", "Borrelia")
HOST_WORDS = ("mitochondri", "chloroplast", "eukaryota")
BANDS = [(-0.01, 0.5, "under 50%"), (0.5, 0.9, "50 to 90%"),
         (0.9, 0.99, "90 to 99%"), (0.99, 1.01, "99% and over")]
PERMUTATIONS = 9999
SEED = 20260924
TRUE_WORDS = ("1", "1.0", "true")


def arms():
    try:
        import check_batch_arms as module
    except ImportError as problem:
        sys.exit("cannot import check_batch_arms (%s). Run me from the project "
                 "root." % problem)
    return module


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s)." % problem)
    return tick_lookup.load(), tick_lookup


def animal_of(table, module, name):
    """I turn a library name into its animal, or ''."""
    text = str(name).strip()
    stem = text.split(".")[0]
    for marker, value in (("16S", text), ("16S", stem), ("specimen", stem)):
        hit = module.resolve(table, marker, value)
        if hit is not None:
            return str(hit)
    return ""


def truthy(value):
    return str(value).strip().lower() in TRUE_WORDS


def run_of(library):
    text = str(library)
    return text.rsplit(".", 1)[1] if "." in text else "unknown"


def read_counts(table, module):
    """I return the decontaminated counts as libraries x features."""
    if not os.path.isfile(COUNTS):
        sys.exit("missing %s. Run me from the project root." % COUNTS)
    frame = pd.read_csv(COUNTS, sep="\t", index_col=0, low_memory=False)
    frame = frame.loc[:, [c for c in frame.columns
                          if pd.to_numeric(frame[c], errors="coerce").notna().all()]]
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0)
    rows = sum(1 for name in frame.index if animal_of(table, module, name))
    columns = sum(1 for name in frame.columns if animal_of(table, module, name))
    if columns > rows:
        frame = frame.T
    return frame


def deer_with_ixodes(table, module):
    """I return the set of deer animals carrying Ixodes and not Dermacentor."""
    if not os.path.isfile(SPECIES):
        sys.exit("missing %s" % SPECIES)
    species = pd.read_csv(SPECIES, sep="\t", dtype=str)
    for column in ("specimen", "host_bin", "Ixodes_pos", "Dermacentor_pos"):
        if column not in species.columns:
            sys.exit("%s has no column '%s'" % (SPECIES, column))
    species["animal"] = [animal_of(table, module, s) for s in species["specimen"]]
    keep = set()
    for animal, host, ixodes, derm in zip(species["animal"], species["host_bin"],
                                          species["Ixodes_pos"],
                                          species["Dermacentor_pos"]):
        if not animal:
            continue
        lowered = str(host).lower()
        if ("deer" in lowered or "wtd" in lowered) and truthy(ixodes) \
                and not truthy(derm):
            keep.add(animal)
    # The exclusion list lives in build_core_microbiome so there is one copy.
    # Two scripts define this cohort and a copied list is how they drifted.
    try:
        from build_core_microbiome import EXCLUDED_NO_TICK_MATERIAL
    except ImportError as problem:
        sys.exit("cannot import the exclusion list from build_core_microbiome "
                 "(%s). Both definitions of the deer cohort must use the same "
                 "list." % problem)
    dropped = set()
    for listed in EXCLUDED_NO_TICK_MATERIAL:
        animal = animal_of(table, module, listed)
        if not animal:
            sys.exit("the exclusion list names %s and tick_lookup does not "
                     "resolve it. I will not drop an animal I cannot name."
                     % listed)
        if animal in keep:
            keep.discard(animal)
            dropped.add(animal)
    print("excluded %d of the %d listed libraries holding no tick material%s"
          % (len(dropped), len(EXCLUDED_NO_TICK_MATERIAL),
             (": " + ", ".join(sorted(dropped))) if dropped else ""))
    return keep


def band_of(share):
    for low, high, name in BANDS:
        if low < share <= high:
            return name
    return BANDS[0][2]


def detection_table(genus_counts, label, depth_note):
    """I print detection by Rickettsia band for one count table."""
    shares = genus_counts.div(genus_counts.sum(axis=1), axis=0)
    rick = shares["Rickettsia"] if "Rickettsia" in shares.columns else 0 * shares.iloc[:, 0]
    bands = rick.map(band_of)
    print("  %-14s %5s %9s %11s %12s %8s %8s"
          % ("Rickettsia", "n", depth_note, "Anaplasma", "Borreliella",
             "either", "both"))
    for _, _, name in BANDS:
        rows = genus_counts.index[bands == name]
        if not len(rows):
            continue
        part = genus_counts.loc[rows]
        present = {}
        for pathogen in PATHOGENS:
            present[pathogen] = (part[pathogen] > 0) if pathogen in part.columns \
                else pd.Series(False, index=rows)
        either = present["Anaplasma"] | present["Borreliella"]
        both = present["Anaplasma"] & present["Borreliella"]
        print("  %-14s %5d %9.0f %10.0f%% %11.0f%% %7.0f%% %7.0f%%"
              % (name, len(rows), part.sum(axis=1).median(),
                 100 * present["Anaplasma"].mean(),
                 100 * present["Borreliella"].mean(),
                 100 * either.mean(), 100 * both.mean()))
    print()
    return rick, bands


def rarefy(genus_counts, depth, rng):
    """I subsample every row to the same depth, dropping rows that are short."""
    keep = genus_counts.index[genus_counts.sum(axis=1) >= depth]
    values = genus_counts.loc[keep].values.astype(np.int64)
    try:
        drawn = np.vstack([rng.multivariate_hypergeometric(row, depth)
                           for row in values])
    except (AttributeError, ValueError) as problem:
        sys.exit("cannot rarefy (%s). numpy 1.18 or newer is needed for "
                 "multivariate_hypergeometric." % problem)
    return pd.DataFrame(drawn, index=keep, columns=genus_counts.columns)


def permutation_test(rick, positive, rng):
    """I return (difference in mean Rickettsia share, p) for positive vs not."""
    values = rick.values.astype(float)
    flags = positive.values.astype(bool)
    if flags.all() or not flags.any():
        return float("nan"), float("nan")
    observed = values[flags].mean() - values[~flags].mean()
    count = 0
    for _ in range(PERMUTATIONS):
        shuffled = rng.permutation(flags)
        if abs(values[shuffled].mean() - values[~shuffled].mean()) >= abs(observed):
            count += 1
    return observed, (count + 1.0) / (PERMUTATIONS + 1.0)


def main():
    module = arms()
    table, lookup_module = lookup()
    counts = read_counts(table, lookup_module)
    print("counts: %d libraries x %d features" % counts.shape)

    path, covered, mapping = module.read_taxonomy(list(counts.columns))
    print("taxonomy: %s (%d of %d features)" % (path, covered, counts.shape[1]))
    host_like = [f for f in counts.columns
                 if any(word in str(mapping.get(f, "")).lower()
                        for word in HOST_WORDS)]
    counts = counts.drop(columns=host_like)
    print("dropped %d host, chloroplast or eukaryotic features" % len(host_like))

    genus = pd.Series({f: module.genus_label(mapping.get(f, ""))
                       for f in counts.columns})
    animals = pd.Series([animal_of(table, lookup_module, s) for s in counts.index],
                        index=counts.index)
    depth = counts.sum(axis=1)
    frame = pd.DataFrame({"animal": animals, "depth": depth})
    frame = frame[(frame["animal"] != "") & (frame["depth"] >= MIN_READS)]
    frame = frame.sort_values("depth", ascending=False).drop_duplicates("animal")
    keep_animals = deer_with_ixodes(table, lookup_module)
    chosen = frame[frame["animal"].isin(keep_animals)]
    if chosen.empty:
        sys.exit("no deer library qualifies. Check %s." % SPECIES)
    genus_counts = counts.loc[chosen.index].T.groupby(genus).sum().T
    print("deer with Ixodes, one library each: %d (median depth %.0f reads)"
          % (len(chosen), chosen["depth"].median()))
    print()

    print("=" * 78)
    print("A. DETECTION BY RICKETTSIA BAND, AS SEQUENCED")
    print("=" * 78)
    rick, bands = detection_table(genus_counts, "as sequenced", "depth")

    print("=" * 78)
    print("B. THE SAME AFTER RAREFYING EVERY ANIMAL TO %d READS" % RAREFY_TO)
    print("=" * 78)
    rng = np.random.default_rng(SEED)
    rarefied = rarefy(genus_counts, RAREFY_TO, rng)
    print("  %d of %d animals reach %d reads" % (len(rarefied), len(chosen),
                                                 RAREFY_TO))
    detection_table(rarefied, "rarefied", "depth")

    print("=" * 78)
    print("C. PERMUTATION TEST ON THE RAREFIED DATA (%d permutations, seed %d)"
          % (PERMUTATIONS, SEED))
    print("=" * 78)
    shares = rarefied.div(rarefied.sum(axis=1), axis=0)
    rick_rare = shares["Rickettsia"] if "Rickettsia" in shares.columns \
        else 0 * shares.iloc[:, 0]
    print("  %-14s %8s %28s %9s" % ("pathogen", "positive",
                                    "mean Rickettsia share", "p"))
    for pathogen in PATHOGENS:
        if pathogen not in rarefied.columns:
            continue
        positive = rarefied[pathogen] > 0
        difference, p = permutation_test(rick_rare, positive, rng)
        if pd.isna(difference):
            continue
        print("  %-14s %8d  positive %.3f, negative %.3f  %8.4f"
              % (pathogen, int(positive.sum()), rick_rare[positive].mean(),
                 rick_rare[~positive].mean(), p))
    print()
    print("  A negative difference means pathogen-positive animals carry less")
    print("  Rickettsia, which is the direction the cell culture work predicts.")
    print()

    print("=" * 78)
    print("D. IS THIS A PLATE EFFECT? BAND BY SEQUENCING RUN")
    print("=" * 78)
    runs = pd.Series([run_of(s) for s in chosen.index], index=chosen.index)
    crossed = pd.crosstab(bands.reindex(chosen.index), runs)
    print(crossed.to_string())
    print()
    biggest = runs.value_counts().index[0]
    inside = chosen.index[runs == biggest]
    print("  within the largest run alone (%s, %d animals):"
          % (biggest, len(inside)))
    detection_table(genus_counts.loc[inside], "one run", "depth")


if __name__ == "__main__":
    main()
