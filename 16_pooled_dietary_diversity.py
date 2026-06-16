#!/usr/bin/env python3
"""
16_pooled_dietary_diversity.py
==============================
Pooled multi-marker descriptive dietary analysis combining MiFish 12S,
cytochrome b, and 18S V9 invertebrate detections into a coherent
feeding-ecology synthesis.

Outputs
───────
  • Pooled per-bird × prey-species presence/absence matrix
  • Per-species × per-season detection counts and frequencies
  • Habitat-category × per-season composition (table + stacked figure)
  • Per-bird marine prey ratio + Kruskal-Wallis test by season
  • Per-bird marine ratio strip plot

This script is descriptive. It is intentionally not a hypothesis test for
DvT or COD structure on the pooled dietary matrix — Jaccard PERMANOVA at
species-level resolution is underpowered at this sample size and produces
results that do not add to the descriptive synthesis. Marine ratio is
retained as a habitat-summary metric with a single Kruskal-Wallis test
across seasons (validation that the pooled approach recovers known
freshwater↔marine ecology).

DESIGN DECISIONS (5/18 planning session)
─────────────────────────────────────────
  • Inclusion: any-marker passing read-depth threshold.
  • Taxonomic unit: species-level. Cytb 'Unclassified prey' bulk dropped.
  • 18S inverts: PA called by absolute read count (default ≥10), not
    relabund, to avoid host-dominance bias.
  • Marine group treated as ecological reference (kept in ratio analysis).
  • Pacifastacus gambelii reclassified to native Faxonius (biogeographic).

18S INVERTEBRATE LOOKUP
───────────────────────
  Nine verified invertebrate taxa per BLAST + manual triage (5/18 session).

Usage
─────
    python 16_pooled_dietary_diversity.py \\
        --mifish-counts results/MiFish/all/taxonomy_annotated/taxonomy_counts_annotated_MiFish.tsv \\
        --mifish-annot  results/MiFish/all/taxonomy_annotated/annotation_table.tsv \\
        --cytb-counts   results/cytb/all/taxonomy_annotated/taxonomy_counts_annotated_cytb.tsv \\
        --cytb-annot    results/cytb/all/taxonomy_annotated/annotation_table.tsv \\
        --inv18s-table  18S_table_export/feature-table.tsv \\
        --metadata      metadata/qiime/metadata_MiFish.tsv \\
        --outdir        results/multimarker/pooled/
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu, kruskal

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# 18S invertebrate ASV lookup (9 verified taxa, post-triage 5/18)
# ─────────────────────────────────────────────────────────────────────────────
INVERT_LOOKUP: Dict[str, Tuple[str, str, str]] = {
    # ASV_ID (full hash)                  : (species,             habitat,      common_group)
    "d78b1982a1470fc9aeb1513d58335db2"    : ("Daphnia galeata",     "Freshwater", "Freshwater zooplankton"),
    "712b4c2302427aad4f22d841bb398884"    : ("Daphnia pulicaria",   "Freshwater", "Freshwater zooplankton"),
    "f929f19ebd2d16c7744753165c423572"    : ("Faxonius sp.",        "Freshwater", "Freshwater crustacean"),   # reclassified from Pacifastacus gambelii
    "da916174c2bdaa6eccded6f4d897855c"    : ("Byblis gaimardi",     "Marine",     "Marine amphipod"),
    "6282de4bbcfaae277d6eab13f0ff7702"    : ("Erpobdella costata",  "Freshwater", "Freshwater leech"),
    "33c3cffda5914c6cee2617e721edfea4"    : ("Caulleriella parva",  "Marine",     "Marine polychaete"),
    "f3d498a518847cc868d0686ef514af5b"    : ("Centropages typicus", "Marine",     "Marine copepod"),
    "9392f1f2f01d352158fba7411a674dfb"    : ("Neohela monstrosa",   "Marine",     "Marine amphipod"),
    "e6eae451c70477fa95453a8105a921de"    : ("Hatschekia sp.",      "Marine",     "Fish parasite (via prey)"),
}

# Habitat category groupings for the marine-ratio metric
MARINE_HABITATS     = {"Marine", "Estuarine", "Marine/Estuary", "Marine/Anadromous"}
FRESHWATER_HABITATS = {"Freshwater"}
ANADROMOUS_HABITATS = {"Anadromous"}

SEASON_ORDER = ["Breeding", "Freshwater_Nonbreeding", "Saltwater"]

# Habitat category buckets used for the stacked composition plot and table.
# Order here = stacking order (bottom → top) and table column order.
HABITAT_CATEGORIES = ["Marine", "Estuarine", "Anadromous", "Freshwater", "Other"]

HABITAT_COLORS = {
    "Marine":     "#0072B2",  # deep blue
    "Estuarine":  "#56B4E9",  # sky blue
    "Anadromous": "#E69F00",  # orange (transition)
    "Freshwater": "#009E73",  # green
    "Other":      "#999999",  # grey
}

def collapse_habitat(h: str) -> str:
    """Map detailed habitat label → broad category for the composition plot/table."""
    if h in {"Marine", "Marine/Estuary"}:
        return "Marine"
    if h == "Estuarine":
        return "Estuarine"
    if h in {"Anadromous", "Marine/Anadromous"}:
        return "Anadromous"
    if h == "Freshwater":
        return "Freshwater"
    return "Other"


# ─────────────────────────────────────────────────────────────────────────────
# Species-name cleanup and biogeographic reclassifications
# ─────────────────────────────────────────────────────────────────────────────
# Reclassifications follow the manuscript Methods (BLAST-verified, biogeographic):
#   • Sprattus sprattus  → Clupea harengus     (Atlantic herring; 99.5% BLAST)
#   • Sardina pilchardus → Clupea harengus     (89–91% BLAST; lower confidence)
#   • Brevoortia patronus → Brevoortia tyrannus (Atlantic menhaden, gulf reference)
#   • Trachurus trachurus → Trachurus lathami  (rough scad; BLAST 5/18)
# T. lathami reclassification applies to BOTH Trachurus ASVs (genus-only and
# species-level classifier calls) — the underlying BLAST result is the same.

RECLASSIFICATIONS: Dict[str, str] = {
    "Sprattus sprattus":   "Clupea harengus",
    "Sardina pilchardus":  "Clupea harengus",
    "Brevoortia patronus": "Brevoortia tyrannus",
    "Trachurus trachurus": "Trachurus lathami",
    "Trachurus sp.":       "Trachurus lathami",  # catches genus-only classifier calls
    "Trachurus":           "Trachurus lathami",  # catches bare-genus rows
}


def extract_species_name(taxonomy_string: str) -> str:
    """
    Extract a display species name from a full taxonomy path.

    Examples
    --------
    >>> extract_species_name("Eukaryota;Chordata;...;Brevoortia;Brevoortia patronus")
    'Brevoortia patronus'
    >>> extract_species_name("Eukaryota;Chordata;...;Micropterus;Unclassified")
    'Micropterus sp.'
    >>> extract_species_name("k__Metazoa;...;g__Alosa;s__pseudoharengus")
    'Alosa pseudoharengus'
    >>> extract_species_name("k__Metazoa;...;g__Brevoortia;Unclassified")
    'Brevoortia sp.'
    >>> extract_species_name("Faxonius sp.")    # already-clean (INVERT_LOOKUP)
    'Faxonius sp.'
    """
    if ";" not in taxonomy_string:
        return taxonomy_string

    parts = [p.strip() for p in taxonomy_string.split(";")]
    parts = [p.split("__", 1)[1] if "__" in p else p for p in parts]
    parts = [p for p in parts if p]
    if not parts:
        return taxonomy_string

    last = parts[-1]
    prev = next((p for p in reversed(parts[:-1])
                 if p and p.lower() != "unclassified"), None)

    if last.lower() == "unclassified":
        return f"{prev} sp." if prev else "Unclassified"
    if " " in last:                                           # already "Genus species"
        return last
    if last[:1].islower() and prev:                           # species epithet → combine
        return f"{prev} {last}"
    return last


def clean_and_reclassify(
    pa: pd.DataFrame,
    info: pd.DataFrame,
    marker: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Map full taxonomy strings → display species names; apply biogeographic
    reclassifications; collapse rows that now share a display name.
    Logs each reclassification that fires so they're auditable.
    """
    if pa.empty:
        return pa, info

    extracted = {row: extract_species_name(row) for row in pa.index}
    final_names = {row: RECLASSIFICATIONS.get(name, name)
                   for row, name in extracted.items()}

    # Log reclassifications that fired
    seen = set()
    for row in pa.index:
        e, r = extracted[row], final_names[row]
        if e != r and (e, r) not in seen:
            log.info("  %s: reclassified %s → %s", marker, e, r)
            seen.add((e, r))

    pa = pa.copy()
    pa.index = pa.index.map(final_names)
    pa = pa.groupby(level=0).max()                            # collapse merged rows

    if not info.empty:
        info = info.copy()
        info.index = info.index.map(final_names)
        info = info[~info.index.duplicated(keep="first")]     # keep first habitat for merged rows

    return pa, info


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def shorten_sample_id(sample_id: str) -> str:
    """TV230007-GI-MiFish_S123 → TV230007"""
    return sample_id.split("-")[0] if "-" in sample_id else sample_id


def load_metadata(path: Path) -> pd.DataFrame:
    """Load QIIME2-style metadata; drop the #q2:types row."""
    df = pd.read_csv(path, sep="\t", index_col=0, dtype=str)
    df = df[~df.index.str.startswith("#", na=False)]
    df.index = [shorten_sample_id(i) for i in df.index]
    df = df[~df.index.duplicated(keep="first")]
    return df


def dedup_columns_max(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicate columns by taking the maximum across replicates."""
    return df.T.groupby(level=0).max().T


def load_annotated_marker(
    counts_path: Path,
    annot_path: Path,
    marker: str,
    min_reads: int,
    min_relabund: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load an annotated count table (output of 11c_annotate_diet_ecology.py).

    Order of operations matters:
      1. Compute per-sample read total on the FULL marker table (incl. unresolved).
      2. Apply min_reads filter on that total. This keeps birds with adequate
         marker sequencing even if cytb is mostly Chordata Unclassified.
      3. Then filter rows to species-resolved taxa.
      4. Recompute within-sample relabund over the resolved subset (rows that
         remain), and convert to PA at min_relabund.

    Returns species-level PA + taxon info DataFrame.
    """
    log.info("Loading %s annotated counts: %s", marker, counts_path)
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    annot = pd.read_csv(annot_path, sep="\t", index_col=0)
    annot = annot[~annot.index.astype(str).str.startswith("#", na=False)]

    counts_num = counts.select_dtypes(include=[np.number])

    # Step 1+2: sample threshold on FULL marker reads
    full_totals = counts_num.sum(axis=0)
    keep_samples = full_totals[full_totals >= min_reads].index
    counts_num = counts_num[keep_samples]
    log.info("  %s: %d samples pass ≥%d-read threshold (on FULL marker total)",
             marker, len(keep_samples), min_reads)
    if counts_num.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Step 3: drop unresolved / junk taxa
    if "common_group" in annot.columns:
        annot = annot[annot["common_group"] != "Unclassified prey"]
        annot = annot[annot["common_group"].notna()]

    resolved = counts_num.loc[counts_num.index.isin(annot.index)]
    if resolved.empty:
        log.warning("  %s: no rows after filtering to resolved taxa.", marker)
        return pd.DataFrame(), pd.DataFrame()

    # Step 4: within-sample relabund computed over resolved subset (so 1%
    # means 1% of the species-resolved signal — a defensible per-marker PA call)
    resolved_totals = resolved.sum(axis=0).replace(0, np.nan)
    relabund = resolved.div(resolved_totals, axis=1).fillna(0)
    pa = (relabund >= min_relabund).astype(int)

    # Drop birds that have zero resolved-taxon signal even if their full
    # marker total was above threshold (e.g., cytb bird with 5000 reads
    # all of which were Chordata Unclassified)
    nonzero_birds = pa.columns[pa.sum(axis=0) > 0]
    if len(nonzero_birds) < pa.shape[1]:
        log.info("  %s: %d birds dropped (had marker reads but no resolved-taxon detections)",
                 marker, pa.shape[1] - len(nonzero_birds))
    pa = pa[nonzero_birds]

    pa.columns = [shorten_sample_id(c) for c in pa.columns]
    pa = dedup_columns_max(pa)

    cols = [c for c in ["habitat", "common_group"] if c in annot.columns]
    info = annot.loc[pa.index, cols].copy() if cols else pd.DataFrame(index=pa.index)
    info["marker"] = marker

    # Strip taxonomy strings to display species names and apply manuscript
    # biogeographic reclassifications (Sprattus → Clupea, B. patronus →
    # B. tyrannus, Trachurus → T. lathami).
    pa, info = clean_and_reclassify(pa, info, marker)

    return pa, info


def load_18s_inverts(
    table_path: Path,
    min_reads: int,
    min_invert_reads: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load 18S feature-table.tsv (QIIME2 export). First line is
    '# Constructed from biom file'; second line is the header starting
    with '#OTU ID'. Skip the first line and let pandas parse the rest.

    For inverts specifically: PA is called by ABSOLUTE READ COUNT, not
    relative abundance. Inverts are a minority signal against a host-
    dominated eukaryotic pool, and a within-sample relabund threshold
    creates a bias toward samples with low eukaryotic background. The
    standard for dietary metabarcoding occurrence data is absolute read
    count (default: ≥10 reads of the ASV in a sample = detected).
    """
    log.info("Loading 18S feature table: %s", table_path)
    table = pd.read_csv(table_path, sep="\t", index_col=0, skiprows=1)
    table = table.select_dtypes(include=[np.number])

    sample_totals = table.sum(axis=0)
    keep_samples = sample_totals[sample_totals >= min_reads].index
    table = table[keep_samples]
    log.info("  18S: %d samples pass ≥%d-read threshold (on FULL eukaryotic table)",
             len(keep_samples), min_reads)

    invert_asvs = [a for a in INVERT_LOOKUP if a in table.index]
    missing = [a for a in INVERT_LOOKUP if a not in table.index]
    if missing:
        log.warning("  18S: %d invert ASVs NOT found in table: %s",
                    len(missing), [a[:8] for a in missing])
    if not invert_asvs:
        return pd.DataFrame(), pd.DataFrame()

    inv = table.loc[invert_asvs]
    # Absolute-count PA call (not relabund)
    pa_asv = (inv >= min_invert_reads).astype(int)
    log.info("  18S inverts: ≥%d-read absolute threshold; %d ASVs × %d samples detection landscape",
             min_invert_reads, pa_asv.shape[0], pa_asv.shape[1])

    species_map = {asv: INVERT_LOOKUP[asv][0] for asv in invert_asvs}
    pa_asv.index = pa_asv.index.map(species_map)
    pa = pa_asv.groupby(level=0).max()

    pa.columns = [shorten_sample_id(c) for c in pa.columns]
    pa = dedup_columns_max(pa)

    # Drop birds with no invert detections from the per-marker matrix
    # (they still appear via MiFish/cytb if those markers caught them)
    nonzero = pa.columns[pa.sum(axis=0) > 0]
    log.info("  18S inverts: %d birds with at least one invert detection", len(nonzero))
    pa = pa[nonzero]

    info_rows = []
    for asv in invert_asvs:
        sp, hab, grp = INVERT_LOOKUP[asv]
        info_rows.append({"species": sp, "habitat": hab, "common_group": grp, "marker": "18S_invert"})
    info = pd.DataFrame(info_rows).drop_duplicates("species").set_index("species")
    info = info.loc[info.index.isin(pa.index)]

    # Apply the same cleanup/reclassification pass for consistency (no-op for
    # current invert names but covers future additions)
    pa, info = clean_and_reclassify(pa, info, "18S_invert")

    return pa, info


def pool_markers(pa_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Union across markers. Rows = species, cols = birds. Value=1 if any marker detected."""
    all_birds   = sorted({b for pa in pa_dict.values() for b in pa.columns})
    all_species = sorted({s for pa in pa_dict.values() for s in pa.index})
    pooled = pd.DataFrame(0, index=all_species, columns=all_birds, dtype=int)
    for marker, pa in pa_dict.items():
        sp = pa.index.intersection(pooled.index)
        bd = pa.columns.intersection(pooled.columns)
        pooled.loc[sp, bd] = np.maximum(
            pooled.loc[sp, bd].values.astype(int),
            pa.loc[sp, bd].values.astype(int),
        )
    return pooled


def per_bird_marker_count(pa_dict: Dict[str, pd.DataFrame]) -> pd.Series:
    """Series: bird → number of markers in which the bird passed read threshold."""
    all_birds = sorted({b for pa in pa_dict.values() for b in pa.columns})
    return pd.Series(
        {b: sum(1 for pa in pa_dict.values() if b in pa.columns) for b in all_birds},
        name="markers_passing",
    )


def freshwater_marine_ratio(pa: pd.DataFrame, habitat_lookup: pd.Series) -> pd.DataFrame:
    """
    Per bird: marine, freshwater, anadromous detection counts + marine ratio.
    Ratio = marine / (marine + freshwater); anadromous excluded from denominator
    to keep the metric clean (anadromous species span both habitats by life
    stage). Anadromous counts are reported separately for context.
    """
    h = habitat_lookup.reindex(pa.index).fillna("Unknown")
    marine_taxa     = h[h.isin(MARINE_HABITATS)].index
    freshwater_taxa = h[h.isin(FRESHWATER_HABITATS)].index
    anadromous_taxa = h[h.isin(ANADROMOUS_HABITATS)].index

    marine     = pa.loc[pa.index.isin(marine_taxa)].sum(axis=0)
    freshwater = pa.loc[pa.index.isin(freshwater_taxa)].sum(axis=0)
    anadromous = pa.loc[pa.index.isin(anadromous_taxa)].sum(axis=0)
    denom = marine + freshwater
    ratio = (marine / denom).where(denom > 0, np.nan)

    return pd.DataFrame({
        "marine_taxa":     marine,
        "freshwater_taxa": freshwater,
        "anadromous_taxa": anadromous,
        "marine_plus_freshwater": denom,
        "marine_ratio":    ratio,
    })


def plot_marine_ratio_by_season(ratio_df: pd.DataFrame, outdir: Path):
    """Strip plot: per-bird marine ratio across seasons, colored by Group."""
    plot_df = ratio_df.dropna(subset=["marine_ratio", "Season"]) if "Season" in ratio_df else None
    if plot_df is None or plot_df.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor("white")
    colors = {"Diseased": "#D55E00", "Trauma": "#0072B2", "Marine": "#009E73"}
    for i, season in enumerate(SEASON_ORDER):
        sub = plot_df[plot_df["Season"] == season]
        for g, color in colors.items():
            grp = sub[sub.get("Group", pd.Series()) == g] if "Group" in sub.columns else pd.DataFrame()
            if grp.empty:
                continue
            jitter = np.random.uniform(-0.12, 0.12, size=len(grp))
            ax.scatter(np.full(len(grp), i) + jitter, grp["marine_ratio"],
                       color=color, s=60, alpha=0.75, edgecolor="white", linewidth=0.8,
                       label=g if i == 0 else None)
    ax.set_xticks(range(len(SEASON_ORDER)))
    ax.set_xticklabels(SEASON_ORDER)
    ax.set_ylabel("Marine prey ratio\n(marine / [marine + freshwater] detections)")
    ax.set_ylim(-0.05, 1.05)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    ax.set_title("Per-bird marine prey ratio by season — pooled multi-marker dietary signal")
    fig.tight_layout()
    fig.savefig(outdir / "marine_ratio_by_season.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "marine_ratio_by_season.svg", bbox_inches="tight")
    plt.close(fig)


def build_taxonomy_by_season_table(
    pa: pd.DataFrame,
    taxon_info: pd.DataFrame,
    season: pd.Series,
) -> pd.DataFrame:
    """
    Per-species × per-season detection counts and frequencies.

    Returns one row per detected taxon, with columns:
      species, habitat, marker, total_n,
      <Season>_n  (count of birds in that season with a detection),
      <Season>_n_cohort  (count of birds in that season total),
      <Season>_freq  (count/cohort).

    Sorted by habitat then total_n (desc) so freshwater taxa appear together,
    marine together, etc., and the most-frequently-detected lead each group.
    """
    season_cohort = {s: int((season == s).sum()) for s in SEASON_ORDER}
    rows = []
    for sp in pa.index:
        info_row = taxon_info.loc[sp] if sp in taxon_info.index else None
        habitat = info_row["habitat"] if info_row is not None and "habitat" in info_row.index else "Unknown"
        marker = info_row["marker"] if info_row is not None and "marker" in info_row.index else "Unknown"

        row = {"species": sp, "habitat": habitat, "marker": marker}
        total_n = 0
        for s in SEASON_ORDER:
            birds_in_s = season[season == s].index
            cols = [b for b in pa.columns if b in birds_in_s]
            n_det = int(pa.loc[sp, cols].sum()) if cols else 0
            row[f"{s}_n"] = n_det
            row[f"{s}_n_cohort"] = season_cohort[s]
            row[f"{s}_freq"] = round(n_det / season_cohort[s], 3) if season_cohort[s] else 0.0
            total_n += n_det
        row["total_n"] = total_n
        rows.append(row)

    df = pd.DataFrame(rows)
    # Drop taxa with zero detections across all seasons (should be rare but possible)
    df = df[df["total_n"] > 0]
    df = df.sort_values(["habitat", "total_n"], ascending=[True, False]).reset_index(drop=True)
    return df


def build_habitat_by_season_table(
    pa: pd.DataFrame,
    taxon_info: pd.DataFrame,
    season: pd.Series,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Per-bird habitat composition, averaged within each season.

    For each bird:
      1. Map each detected species → broad habitat (Marine/Estuarine/Anadromous/Freshwater/Other).
      2. Compute proportion of detections in each habitat category.
    For each season:
      3. Average the per-bird proportions across birds with ≥1 detection.

    Returns (long_form, wide_form):
      long_form: rows = (season, habitat), cols = mean_proportion, n_birds_with_detection,
                 n_total_cohort
      wide_form: rows = season, cols = habitat (for direct plotting)
    """
    h = taxon_info["habitat"].reindex(pa.index).fillna("Unknown")
    h_broad = h.map(collapse_habitat)

    # Per-bird per-habitat raw detection counts
    bird_habitat = pd.DataFrame(0, index=pa.columns, columns=HABITAT_CATEGORIES, dtype=int)
    for cat in HABITAT_CATEGORIES:
        species_in_cat = h_broad[h_broad == cat].index
        if len(species_in_cat):
            common = pa.index.intersection(species_in_cat)
            bird_habitat[cat] = pa.loc[common].sum(axis=0)

    # Per-bird proportions (rows summing to 1 where bird had any detection)
    totals = bird_habitat.sum(axis=1)
    bird_props = bird_habitat.div(totals.replace(0, np.nan), axis=0)

    # Average within season
    long_rows = []
    for s in SEASON_ORDER:
        birds_in_s = season[season == s].index
        n_cohort = len(birds_in_s)
        # Birds with detection AND in this season
        cohort_with_det = bird_props.loc[bird_props.index.isin(birds_in_s) &
                                          bird_props.notna().any(axis=1)]
        n_with_det = len(cohort_with_det)
        for cat in HABITAT_CATEGORIES:
            mean_prop = cohort_with_det[cat].mean() if n_with_det else 0.0
            long_rows.append({
                "season": s,
                "habitat": cat,
                "mean_proportion": round(float(mean_prop) if pd.notna(mean_prop) else 0.0, 3),
                "n_birds_with_detection": n_with_det,
                "n_total_cohort": n_cohort,
            })

    long_df = pd.DataFrame(long_rows)
    wide_df = long_df.pivot(index="season", columns="habitat",
                             values="mean_proportion").reindex(SEASON_ORDER)
    wide_df = wide_df[[c for c in HABITAT_CATEGORIES if c in wide_df.columns]]
    return long_df, wide_df


def plot_habitat_composition_stacked(
    wide_df: pd.DataFrame,
    long_df: pd.DataFrame,
    outdir: Path,
):
    """
    Stacked bar: mean per-bird habitat composition by season.
    Each bar sums to ~1; segments = habitat categories.
    """
    if wide_df.empty:
        log.warning("habitat-by-season wide_df is empty; skipping stacked plot")
        return

    fig, ax = plt.subplots(figsize=(9, 5.5))
    fig.patch.set_facecolor("white")

    seasons = wide_df.index.tolist()
    bottom = np.zeros(len(seasons))
    handles = []
    for cat in HABITAT_CATEGORIES:
        if cat not in wide_df.columns:
            continue
        vals = wide_df[cat].fillna(0).values
        if vals.sum() == 0:
            continue
        bars = ax.bar(range(len(seasons)), vals, bottom=bottom,
                      color=HABITAT_COLORS[cat], label=cat,
                      edgecolor="white", linewidth=0.6)
        bottom += vals
        handles.append(bars)

    # n labels under each bar
    n_lookup = long_df.drop_duplicates("season").set_index("season")[
        ["n_birds_with_detection", "n_total_cohort"]
    ]
    xticklabels = [
        f"{s}\n(n={int(n_lookup.loc[s, 'n_birds_with_detection'])}"
        f" / {int(n_lookup.loc[s, 'n_total_cohort'])})"
        for s in seasons
    ]
    ax.set_xticks(range(len(seasons)))
    ax.set_xticklabels(xticklabels, fontsize=10)
    ax.set_ylabel("Mean per-bird prey composition\n(proportion of detections)")
    ax.set_ylim(0, 1.05)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              title="Prey habitat", fontsize=9, title_fontsize=10)
    ax.set_title(
        "Pooled multi-marker dietary composition by ecological season\n"
        "MiFish 12S + cytochrome b + 18S V9 invertebrates  ·  per-bird "
        "detection proportions averaged within season"
    )
    fig.tight_layout()
    for ext in ("png", "svg"):
        out = outdir / f"habitat_composition_by_season_stacked.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        log.info("Saved: %s", out)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--mifish-counts", type=Path, required=True)
    p.add_argument("--mifish-annot",  type=Path, required=True)
    p.add_argument("--cytb-counts",   type=Path, required=True)
    p.add_argument("--cytb-annot",    type=Path, required=True)
    p.add_argument("--inv18s-table",  type=Path, required=True,
                   help="18S feature-table.tsv (export of feature-table.qza)")
    p.add_argument("--metadata",      type=Path, required=True)
    p.add_argument("--outdir",        type=Path, required=True)
    p.add_argument("--mifish-min-reads", type=int, default=1000,
                   help="MiFish per-sample threshold (default 1000 — fixed from old ≥500 bug)")
    p.add_argument("--cytb-min-reads",   type=int, default=50,
                   help="cytb per-sample threshold (default 50; cytb has shallower depth)")
    p.add_argument("--inv18s-min-reads", type=int, default=1000,
                   help="18S per-sample threshold on FULL eukaryotic reads (default 1000; matches rarefaction depth)")
    p.add_argument("--inv18s-min-invert-reads", type=int, default=10,
                   help="18S per-invert-ASV absolute read threshold for PA call (default 10). "
                        "Inverts are a minority signal vs host DNA; absolute count avoids the "
                        "bias of relabund thresholds against host-dominated samples.")
    p.add_argument("--min-relabund",  type=float, default=0.01,
                   help="Within-sample relabund threshold for fish markers (default 0.01)")
    p.add_argument("--seed",          type=int, default=42)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    np.random.seed(args.seed)

    # ── Load data ────────────────────────────────────────────────────────────
    meta = load_metadata(args.metadata)
    log.info("Metadata: %d birds, columns = %s", len(meta), list(meta.columns))

    mf_pa, mf_info = load_annotated_marker(args.mifish_counts, args.mifish_annot,
                                            "MiFish", args.mifish_min_reads, args.min_relabund)
    cy_pa, cy_info = load_annotated_marker(args.cytb_counts, args.cytb_annot,
                                            "cytb",   args.cytb_min_reads,   args.min_relabund)
    iv_pa, iv_info = load_18s_inverts(args.inv18s_table, args.inv18s_min_reads, args.inv18s_min_invert_reads)

    pa_dict = {m: pa for m, pa in
               [("MiFish", mf_pa), ("cytb", cy_pa), ("18S_invert", iv_pa)]
               if not pa.empty}
    if not pa_dict:
        log.error("All markers empty after filtering. Aborting.")
        return 1

    # ── Pool ──────────────────────────────────────────────────────────────────
    pooled = pool_markers(pa_dict)
    log.info("Pooled matrix: %d taxa × %d birds", *pooled.shape)
    pooled.to_csv(args.outdir / "pooled_pa_matrix.tsv", sep="\t")

    habitat_lookup = pd.concat(
        [df["habitat"] for df in (mf_info, cy_info, iv_info) if "habitat" in df.columns],
    )
    habitat_lookup = habitat_lookup[~habitat_lookup.index.duplicated(keep="first")]

    mp = per_bird_marker_count(pa_dict)
    mp.to_csv(args.outdir / "markers_passing_per_bird.tsv", sep="\t")

    # Combined taxon_info across markers (for taxonomy table & habitat lookup)
    combined_info = pd.concat(
        [df for df in (mf_info, cy_info, iv_info) if not df.empty],
    )
    combined_info = combined_info[~combined_info.index.duplicated(keep="first")]

    cohort = [b for b in pooled.columns if b in meta.index]

    # ── Season vector for the cohort ─────────────────────────────────────────
    season_vec = pd.Series(dtype=str)
    if "Season" in meta.columns:
        season_vec = meta.loc[cohort, "Season"]
        season_vec = season_vec[season_vec.isin(SEASON_ORDER)]

    # ── Per-species × per-season detection table ─────────────────────────────
    if not season_vec.empty:
        tax_season = build_taxonomy_by_season_table(
            pooled[[b for b in pooled.columns if b in season_vec.index]],
            combined_info, season_vec,
        )
        tax_season.to_csv(args.outdir / "taxonomy_by_season.tsv", sep="\t", index=False)
        log.info("Saved: %s (%d species × seasons)",
                 args.outdir / "taxonomy_by_season.tsv", len(tax_season))

        # ── Per-habitat × per-season composition table + stacked plot ────────
        habitat_long, habitat_wide = build_habitat_by_season_table(
            pooled[[b for b in pooled.columns if b in season_vec.index]],
            combined_info, season_vec,
        )
        habitat_long.to_csv(args.outdir / "habitat_by_season.tsv", sep="\t", index=False)
        log.info("Saved: %s", args.outdir / "habitat_by_season.tsv")

        try:
            plot_habitat_composition_stacked(habitat_wide, habitat_long, args.outdir)
        except Exception as e:
            log.warning("Stacked composition plot failed (non-fatal): %s", e)
    else:
        tax_season = pd.DataFrame()
        habitat_long = pd.DataFrame()
        habitat_wide = pd.DataFrame()

    # ── Analysis C: per-bird marine prey ratio ───────────────────────────────
    ratio_df = freshwater_marine_ratio(pooled, habitat_lookup)
    meta_cols = [c for c in ["Season", "Group", "COD_broad"] if c in meta.columns]
    ratio_df = ratio_df.join(meta[meta_cols], how="left")
    ratio_df.to_csv(args.outdir / "per_bird_marine_ratio.tsv", sep="\t")

    results_C: dict = {}
    # KW across seasons
    if "Season" in ratio_df.columns:
        seasonal = {s: ratio_df.loc[(ratio_df["Season"] == s) &
                                    ratio_df["marine_ratio"].notna(), "marine_ratio"].values
                    for s in SEASON_ORDER}
        seasonal = {s: v for s, v in seasonal.items() if len(v) >= 2}
        if len(seasonal) >= 2:
            stat, pval = kruskal(*seasonal.values())
            results_C["KW_Season"] = {
                "label":      "Marine-ratio across seasons (Kruskal-Wallis)",
                "statistic":  float(stat),
                "p_value":    float(pval),
                "group_ns":    {s: int(len(v)) for s, v in seasonal.items()},
                "group_means": {s: float(np.mean(v)) for s, v in seasonal.items()},
            }

    # MWU DvT (Marine excluded)
    if "Group" in ratio_df.columns:
        dvt_r = ratio_df[ratio_df["Group"].isin(["Diseased", "Trauma"]) &
                          ratio_df["marine_ratio"].notna()]
        d = dvt_r.loc[dvt_r["Group"] == "Diseased", "marine_ratio"].values
        t = dvt_r.loc[dvt_r["Group"] == "Trauma",   "marine_ratio"].values
        if len(d) >= 2 and len(t) >= 2:
            stat, pval = mannwhitneyu(d, t, alternative="two-sided")
            results_C["MWU_DvT"] = {
                "label":         "Marine-ratio Diseased vs Trauma (Mann-Whitney U)",
                "statistic":     float(stat),
                "p_value":       float(pval),
                "diseased_n":    int(len(d)),
                "diseased_mean": float(np.mean(d)),
                "trauma_n":      int(len(t)),
                "trauma_mean":   float(np.mean(t)),
            }

    # Plot
    try:
        plot_marine_ratio_by_season(ratio_df, args.outdir)
    except Exception as e:
        log.warning("Plot failed (non-fatal): %s", e)

    # ── Summary report ──────────────────────────────────────────────────────
    rpt = args.outdir / "pooled_dietary_summary.txt"
    with open(rpt, "w") as f:
        f.write("=" * 72 + "\n")
        f.write("POOLED MULTI-MARKER DIETARY DIVERSITY — SUMMARY\n")
        f.write("=" * 72 + "\n\n")
        f.write("Markers (taxa × birds after filtering):\n")
        for m, pa in [("MiFish 12S", mf_pa), ("cytochrome b", cy_pa), ("18S inverts", iv_pa)]:
            if not pa.empty:
                f.write(f"  {m:14s}: {pa.shape[0]:3d} × {pa.shape[1]:3d}\n")
            else:
                f.write(f"  {m:14s}: empty\n")
        f.write(f"\nPooled matrix: {pooled.shape[0]} taxa × {pooled.shape[1]} birds\n\n")

        f.write("Bird marker coverage:\n")
        for k, v in mp.value_counts().sort_index().items():
            f.write(f"  {k} marker(s) passing: {v} birds\n")
        f.write("\n")

        f.write("─" * 72 + "\n")
        f.write("HABITAT COMPOSITION BY SEASON (per-bird proportions, averaged)\n")
        f.write("─" * 72 + "\n\n")
        if not habitat_wide.empty:
            f.write(habitat_wide.to_string(float_format="%.3f"))
            f.write("\n\n")
        else:
            f.write("  (no season information available)\n\n")

        f.write("─" * 72 + "\n")
        f.write("MARINE PREY RATIO BY SEASON\n")
        f.write("─" * 72 + "\n\n")
        for key, r in results_C.items():
            f.write(f"  {r['label']}\n")
            f.write(f"    statistic = {r['statistic']:.4f}\n")
            f.write(f"    p-value   = {r['p_value']:.4f}\n")
            if "group_ns" in r:
                for s in SEASON_ORDER:
                    if s in r["group_ns"]:
                        f.write(f"      {s:22s}: n={r['group_ns'][s]:2d}, mean ratio = {r['group_means'][s]:.3f}\n")
            else:
                f.write(f"      Diseased: n={r['diseased_n']:2d}, mean = {r['diseased_mean']:.3f}\n")
                f.write(f"      Trauma:   n={r['trauma_n']:2d}, mean = {r['trauma_mean']:.3f}\n")
            f.write("\n")

        f.write("─" * 72 + "\n")
        f.write("TOP DETECTED TAXA (per-season frequency, full table in taxonomy_by_season.tsv)\n")
        f.write("─" * 72 + "\n\n")
        if not tax_season.empty:
            # Show top 5 per habitat to keep the summary readable
            for habitat in tax_season["habitat"].unique():
                sub = tax_season[tax_season["habitat"] == habitat].head(5)
                if sub.empty:
                    continue
                f.write(f"  {habitat}:\n")
                for _, row in sub.iterrows():
                    parts = []
                    for s in SEASON_ORDER:
                        n, ncoh = row[f"{s}_n"], row[f"{s}_n_cohort"]
                        if n > 0:
                            parts.append(f"{s}={n}/{ncoh}")
                    detail = ", ".join(parts) if parts else "no detections in any season"
                    f.write(f"    {row['species']:30s} [{row['marker']:10s}]  {detail}\n")
                f.write("\n")
        else:
            f.write("  (no taxonomy-by-season data)\n\n")

        f.write("─" * 72 + "\n")
        f.write("Files written:\n")
        f.write("  pooled_pa_matrix.tsv                         — taxa × birds, species-level PA\n")
        f.write("  markers_passing_per_bird.tsv                 — bird → markers contributing\n")
        f.write("  per_bird_marine_ratio.tsv                    — bird-level marine ratio + metadata\n")
        f.write("  taxonomy_by_season.tsv                       — species × season detection counts/frequencies\n")
        f.write("  habitat_by_season.tsv                        — habitat category × season composition\n")
        f.write("  habitat_composition_by_season_stacked.{png,svg} — stacked composition figure\n")
        f.write("  marine_ratio_by_season.{png,svg}             — per-bird ratio strip plot\n")
        f.write("  pooled_dietary_summary.txt                   — this file\n")

    # Echo to stdout
    print(rpt.read_text())
    log.info("Done. Summary at %s", rpt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
