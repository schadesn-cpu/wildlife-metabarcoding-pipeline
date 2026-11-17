#!/usr/bin/env python3
"""check_batch_sources.py

I test three explanations for the deer batch effect. I write nothing.

A. ASV SPLITTING: how Rickettsia reads divide among its top ASVs per batch, and
   batch R2 (one-factor PERMANOVA, 999 permutations, seed 20260922) on
   Bray-Curtis at ASV level versus after collapsing to genus. The ASV column
   must reproduce QIIME (about 13.8% for all deer, 33.5% for the main cluster).
B. KED CARRY-OVER: 18S ked reads (Lipoptena, Hippoboscidae) per ordination group.
C. P3-5 MECHANISM: Anaplasma and Borreliella detection by Rickettsia share and
   by plate.

USAGE
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_batch_sources.py
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_batch_arms as arms  # noqa: E402

ASV_18S = "results/tables/asv_table_18S.tsv"
TAX_18S = "results/18S/all/taxonomy_verified/taxonomy_18S_blast_verified.tsv"
KED_WORDS = ("lipoptena", "hippobosc", "melophagus", "neolipoptena")
PERMUTATIONS = 999
SEED = 20260922


def permanova(distances, labels):
    """I return (R2, p) for a one-factor PERMANOVA on a square matrix."""
    d2 = distances ** 2
    n = len(labels)
    labels = np.asarray(labels)
    total = d2[np.triu_indices(n, 1)].sum() / n

    def within(lab):
        ss = 0.0
        for g in np.unique(lab):
            idx = np.where(lab == g)[0]
            if len(idx) > 1:
                block = d2[np.ix_(idx, idx)]
                ss += block[np.triu_indices(len(idx), 1)].sum() / len(idx)
        return ss

    k = len(np.unique(labels))
    if total <= 0 or k < 2:
        return float("nan"), float("nan")
    observed_w = within(labels)
    f_obs = ((total - observed_w) / (k - 1)) / (observed_w / (n - k))
    rng = np.random.default_rng(SEED)
    exceed = 0
    for _ in range(PERMUTATIONS):
        w = within(rng.permutation(labels))
        exceed += ((total - w) / (k - 1)) / (w / (n - k)) >= f_obs
    return 1 - observed_w / total, (exceed + 1) / (PERMUTATIONS + 1)


def bray(frame):
    return squareform(pdist(frame.values, "braycurtis"))


def load():
    distances = arms.read_matrix()
    samples = list(distances.index)
    coords = arms.principal_coordinates(distances)
    meta = pd.read_csv(arms.META, sep="\t", dtype=str)
    meta = meta[meta["sample-id"] != "#q2:types"].set_index("sample-id")
    work = pd.DataFrame(index=samples)
    work["group"] = [arms.group_of(x, y) for x, y in coords]
    work["batch"] = meta.loc[samples, "batch_clean"].values
    hits = arms.plate_hits(samples)
    work["plate"] = ["P3-5" if s in hits else "other" for s in samples]
    shares = arms.read_rarefied().loc[samples]
    shares = shares.loc[:, shares.sum(axis=0) > 0]
    _, _, mapping = arms.read_taxonomy(list(shares.columns))
    genus = pd.Series({f: arms.genus_label(mapping.get(f, "")) for f in shares.columns})
    return work, shares, genus


def section_a(work, shares, genus):
    print("=" * 72)
    print("A. IS THE BATCH EFFECT IN TYPICAL TICKS AN ASV-SPLITTING ARTIFACT?")
    print("=" * 72)
    rick = [f for f in shares.columns if genus[f] == "Rickettsia"]
    main = work.index[work["group"] == "main"]
    if rick:
        within = shares.loc[main, rick]
        top = within.sum().sort_values(ascending=False).index[:3]
        print("  Rickettsia ASVs: %d. Share of Rickettsia reads in the top three,"
              % len(rick))
        print("  per batch (main-cluster deer):")
        print("  %-22s %4s %8s %8s %8s" % ("batch", "n", "ASV 1", "ASV 2", "ASV 3"))
        for batch, rows in work.loc[main].groupby("batch"):
            part = within.loc[rows.index]
            split = part[top].div(part.sum(axis=1).replace(0, np.nan), axis=0).mean()
            values = list(100 * split.reindex(top).fillna(0).values) + [0.0, 0.0, 0.0]
            print("  %-22s %4d %7.0f%% %7.0f%% %7.0f%%"
                  % (batch[:22], len(rows), values[0], values[1], values[2]))
        print()
    genus_table = shares.T.groupby(genus).sum().T
    print("  batch R2, one-factor PERMANOVA (%d permutations):" % PERMUTATIONS)
    print("  %-26s %14s %14s" % ("set", "ASV level", "genus level"))
    for label, index in (("all deer", work.index), ("main cluster only", main)):
        labels = work.loc[index, "batch"].values
        asv_r2, asv_p = permanova(bray(shares.loc[index]), labels)
        gen_r2, gen_p = permanova(bray(genus_table.loc[index]), labels)
        print("  %-26s %5.1f%% p=%.3f %5.1f%% p=%.3f"
              % ("%s (n = %d)" % (label, len(index)), 100 * asv_r2, asv_p,
                 100 * gen_r2, gen_p))
    print("  The ASV column should read about 13.8% and 33.5% (the QIIME runs).")
    print()


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s)." % problem)
    return tick_lookup.load(), tick_lookup


def resolve(table, module, sample, marker):
    stem = str(sample).split(".")[0]
    for m, value in ((marker, sample), (marker, stem), ("specimen", stem)):
        hit = module.resolve(table, m, value)
        if hit is not None:
            return str(hit)
    return ""


def section_b(work):
    print("=" * 72)
    print("B. IS THE LOWER ARM KED DNA?")
    print("=" * 72)
    if not (os.path.isfile(ASV_18S) and os.path.isfile(TAX_18S)):
        print("  missing %s or %s\n" % (ASV_18S, TAX_18S))
        return
    tax = pd.read_csv(TAX_18S, sep="\t", dtype=str)
    tax_col = next((c for c in tax.columns if c.lower() == "taxon"), tax.columns[1])
    ked = set(tax.loc[tax[tax_col].str.lower().str.contains("|".join(KED_WORDS), na=False),
                      tax.columns[0]])
    if not ked:
        print("  no 18S ASV is assigned to a ked (%s).\n" % ", ".join(KED_WORDS))
        return
    counts = pd.read_csv(ASV_18S, sep="\t", index_col=0, low_memory=False)
    if not set(counts.index) & ked:
        counts = counts.T
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    per_sample = counts.loc[[i for i in counts.index if i in ked]].sum(axis=0)
    table, module = lookup()
    per_animal = {}
    for sample, reads in per_sample.items():
        animal = resolve(table, module, sample, "18S")
        if animal:
            per_animal[animal] = max(per_animal.get(animal, 0), float(reads))
    print("  ked 18S ASVs: %d; 18S libraries matched to an animal: %d"
          % (len(ked & set(counts.index)), len(per_animal)))
    animals = [resolve(table, module, s, "16S") for s in work.index]
    work = work.assign(ked=[per_animal.get(a, np.nan) for a in animals])
    print("  %-6s %4s %10s %16s" % ("group", "n", "with 18S", "carrying ked"))
    for name in ("main", "upper", "lower", "other"):
        rows = work[work["group"] == name]
        if rows.empty:
            continue
        known = int(rows["ked"].notna().sum())
        carrying = int((rows["ked"] > 0).sum())
        print("  %-6s %4d %10d %9d (%.0f%%)" % (name, len(rows), known, carrying,
                                               100.0 * carrying / max(known, 1)))
    print()


def section_c(work, shares, genus):
    print("=" * 72)
    print("C. DOES LOW RICKETTSIA DRIVE MULTI-PATHOGEN DETECTION ON P3-5?")
    print("=" * 72)
    if not os.path.isfile(arms.ANALYSIS):
        print("  missing %s\n" % arms.ANALYSIS)
        return
    calls = pd.read_csv(arms.ANALYSIS, sep="\t", dtype=str)
    calls["stem"] = calls["sample"].astype(str).str.split(".").str[0]
    by_stem = calls.drop_duplicates("stem").set_index("stem")
    rick = shares.loc[:, [f for f in shares.columns if genus[f] == "Rickettsia"]].sum(axis=1)
    rows = []
    for sample in work.index:
        stem = str(sample).split(".")[0]
        if stem in by_stem.index:
            a = pd.to_numeric(by_stem.loc[stem, "Anaplasma_pos"], errors="coerce")
            b = pd.to_numeric(by_stem.loc[stem, "Borreliella_pos"], errors="coerce")
            rows.append((sample, a == 1, b == 1))
    if not rows:
        print("  no deer library matched the analysis table by name.\n")
        return
    frame = pd.DataFrame(rows, columns=["sample", "ana", "bor"]).set_index("sample")
    frame["rick"] = rick.reindex(frame.index).values
    frame["plate"] = work.loc[frame.index, "plate"].values
    frame["band"] = pd.cut(frame["rick"], [-0.01, 0.5, 0.9, 1.01],
                           labels=["<50%", "50-90%", ">=90%"])
    print("  %d of %d deer matched to pathogen calls." % (len(frame), len(work)))
    print("  %-7s %-7s %4s %11s %12s %8s" % ("plate", "Rick.", "n", "Anaplasma",
                                          "Borreliella", "both"))
    for (plate, band), part in frame.groupby(["plate", "band"], observed=True):
        print("  %-7s %-7s %4d %10.0f%% %11.0f%% %7.0f%%"
              % (plate, band, len(part), 100 * part["ana"].mean(),
                 100 * part["bor"].mean(), 100 * (part["ana"] & part["bor"]).mean()))
    print()


def main():
    work, shares, genus = load()
    section_a(work, shares, genus)
    section_b(work)
    section_c(work, shares, genus)


if __name__ == "__main__":
    main()
