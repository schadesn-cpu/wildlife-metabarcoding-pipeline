#!/usr/bin/env python3
"""
plot_combined_habitat_season.py
================================
Combined MiFish 12S + cytochrome b dietary detection figure.

Collapses species-level detections into prey habitat categories
(Marine forage fish, Freshwater prey fish, etc.) and shows detection
frequency per ecological season. Designed for general audiences who
don't know individual fish species.

Each bar = proportion of samples in that season where at least one
species from that habitat category was detected (presence/absence,
not relative abundance).

Usage:
    python plot_combined_habitat_season.py \
        --mifish-counts  results/MiFish/all/taxonomy_cleaned/taxonomy_counts_cleaned_MiFish.tsv \
        --mifish-annot   results/MiFish/all/taxonomy_annotated/annotation_table.tsv \
        --cytb-counts    results/cytb/all/taxonomy_cleaned/taxonomy_relabund_L7_cytb_cleaned.tsv \
        --cytb-annot     results/cytb/all/taxonomy_annotated/annotation_table.tsv \
        --metadata       metadata/qiime/metadata_MiFish_ecoseason.tsv \
        --outdir         results/multimarker/figures/
"""

from __future__ import annotations
import argparse
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Wong 2011 palette ─────────────────────────────────────────────────────────
WONG = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # green
    "#56B4E9",  # sky blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#F0E442",  # yellow
    "#999999",  # grey
]

# Habitat category order — most ecologically meaningful first
CATEGORY_ORDER = [
    "Marine forage fish",
    "Marine bottom fish",
    "Anadromous fish",
    "Marine or anadromous fish",
    "Estuarine fish",
    "Estuarine forage fish",
    "Sand lance",
    "Freshwater prey fish",
    "Freshwater predators",
    # "Unclassified prey" intentionally excluded — cytb reference database
    # does not resolve many sequences below class level, so this bucket
    # absorbed ~100% of detections in every season and obscured the real
    # habitat signature. Report in Methods as a known limitation of the
    # current cytb reference; revisit after targeted NE-coastal-fish
    # reference rebuild.
]

# Simplified display labels for lay audiences
CATEGORY_LABELS = {
    "Marine forage fish":       "Marine forage fish",
    "Marine bottom fish":       "Marine bottom fish",
    "Anadromous fish":          "Anadromous fish",
    "Marine or anadromous fish":"Marine / anadromous fish",
    "Estuarine fish":           "Estuarine fish",
    "Estuarine forage fish":    "Estuarine forage fish",
    "Sand lance":               "Sand lance",
    "Freshwater prey fish":     "Freshwater prey fish",
    "Freshwater predators":     "Freshwater predators",
    "Unclassified prey":        "Unclassified prey",
}

SEASON_ORDER  = ["Breeding", "Freshwater_Nonbreeding", "Saltwater"]
SEASON_LABELS = {
    "Breeding":               "Breeding\n(May–Aug)",
    "Freshwater_Nonbreeding": "Freshwater\nNon-breeding\n(Apr+Sep)",
    "Saltwater":              "Saltwater\n(Oct–Mar)",
}

# ── Data loading ──────────────────────────────────────────────────────────────

def load_counts(path: Path) -> pd.DataFrame:
    """Load a taxonomy count table. Rows = taxa, columns = samples."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    # Drop any non-numeric columns
    return df.select_dtypes(include=[np.number])


def load_annot(path: Path) -> pd.Series:
    """Load annotation table and return species → common_group mapping."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    # Skip the q2:types row if present
    df = df[~df.index.str.startswith("#", na=False)]
    if "common_group" not in df.columns:
        raise ValueError(f"'common_group' column not found in {path}")
    return df["common_group"].dropna()


def load_metadata(path: Path, season_col: str = "Season") -> pd.Series:
    """Load metadata and return sample → season mapping."""
    df = pd.read_csv(path, sep="\t", index_col=0, dtype=str)
    df = df[~df.index.str.startswith("#", na=False)]
    if season_col not in df.columns:
        raise ValueError(f"Column '{season_col}' not found in metadata. "
                         f"Available: {list(df.columns)}")
    season = df[season_col].dropna()
    season = season[season.isin(SEASON_ORDER)]
    return season


def shorten_sample_id(sample_id: str) -> str:
    """TV230007-GI-MiFish_S123 → TV230007"""
    return sample_id.split("-")[0] if "-" in sample_id else sample_id


# ── Core computation ──────────────────────────────────────────────────────────

def counts_to_pa(counts: pd.DataFrame,
                 min_sample_reads: int = 500,
                 min_relabund: float = 0.01) -> pd.DataFrame:
    """
    Convert count table to presence/absence.
    - Drop samples below min_sample_reads
    - Zero out taxa below min_relabund within each sample
    - Convert remainder to 0/1
    """
    totals = counts.sum(axis=0)
    keep = totals[totals >= min_sample_reads].index
    counts = counts[keep]
    if counts.empty:
        return counts

    relabund = counts.div(counts.sum(axis=0), axis=1)
    pa = (relabund >= min_relabund).astype(int)
    return pa


def collapse_to_categories(pa: pd.DataFrame,
                            annot: pd.Series) -> pd.DataFrame:
    """
    Map taxa to common_group categories and collapse.
    Returns DataFrame: rows = categories, cols = samples.
    Value = 1 if any taxon in that category was detected.
    """
    # Map taxa to categories
    mapped = annot.reindex(pa.index).dropna()
    rows = []
    for cat in CATEGORY_ORDER:
        taxa_in_cat = mapped[mapped == cat].index.tolist()
        if not taxa_in_cat:
            continue
        taxa_present = [t for t in taxa_in_cat if t in pa.index]
        if not taxa_present:
            continue
        # Any detection in category = 1
        detected = pa.loc[taxa_present].sum(axis=0).gt(0).astype(int)
        detected.name = cat
        rows.append(detected)

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def compute_detection_freq(cat_pa: pd.DataFrame,
                           season: pd.Series) -> pd.DataFrame:
    """
    For each season and category, compute proportion of samples
    where the category was detected.
    Returns DataFrame: rows = categories, cols = seasons.
    """
    results = {}
    for s in SEASON_ORDER:
        samples_in_season = season[season == s].index.tolist()
        # Shorten IDs to match PA table columns
        short_to_orig = {}
        for col in cat_pa.columns:
            short = shorten_sample_id(col)
            short_to_orig[short] = col

        matched = []
        for sid in samples_in_season:
            short = shorten_sample_id(sid)
            if short in short_to_orig:
                matched.append(short_to_orig[short])
            elif sid in cat_pa.columns:
                matched.append(sid)

        if not matched:
            results[s] = pd.Series(0.0, index=cat_pa.index)
            continue

        sub = cat_pa[matched]
        freq = sub.mean(axis=1)
        results[s] = freq

    return pd.DataFrame(results)


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_figure(freq: pd.DataFrame,
                season_ns: dict,
                outdir: Path,
                palette: str = "wong") -> None:
    """
    Plot grouped bar chart: x = seasons, grouped bars = habitat categories.
    Each bar = detection frequency (0–100%).
    """
    # Filter to categories that have any detections
    freq = freq.loc[(freq > 0).any(axis=1)]
    if freq.empty:
        log.error("No detections to plot.")
        return

    categories = [c for c in CATEGORY_ORDER if c in freq.index]
    seasons    = [s for s in SEASON_ORDER if s in freq.columns]

    n_cats    = len(categories)
    n_seasons = len(seasons)

    colors = WONG[:n_cats]

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor("white")

    bar_width  = 0.12
    group_gap  = 0.15
    x_positions = np.arange(n_seasons) * (n_cats * bar_width + group_gap)

    handles = []
    for i, cat in enumerate(categories):
        offsets = x_positions + i * bar_width
        vals    = [freq.loc[cat, s] * 100 if s in freq.columns else 0
                   for s in seasons]
        bars = ax.bar(offsets, vals, bar_width,
                      color=colors[i], label=CATEGORY_LABELS.get(cat, cat),
                      edgecolor="white", linewidth=0.5)
        handles.append(mpatches.Patch(color=colors[i],
                                      label=CATEGORY_LABELS.get(cat, cat)))

    # X-axis ticks at center of each group
    group_centers = x_positions + (n_cats * bar_width) / 2 - bar_width / 2
    ax.set_xticks(group_centers)
    ax.set_xticklabels(
        [f"{SEASON_LABELS.get(s, s)}\n(n={season_ns.get(s, 0)})"
         for s in seasons],
        fontsize=11,
    )

    ax.set_ylabel("Detection frequency (%)", fontsize=12)
    ax.set_ylim(0, 105)
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0f}%")
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5, zorder=0)

    ax.set_title(
        "MiFish 12S + Cytochrome b — Prey habitat by ecological season\n"
        "(presence/absence, ≥500 reads/sample, ≥1% relative abundance)",
        fontsize=12, pad=12,
    )

    ax.legend(
        handles=handles,
        title="Prey category",
        title_fontsize=9,
        fontsize=8.5,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        framealpha=0.9,
    )

    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)

    for ext in ("png", "svg"):
        out = outdir / f"combined_habitat_season_wong.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        log.info("Saved: %s", out)

    plt.close(fig)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mifish-counts", required=True, type=Path)
    p.add_argument("--mifish-annot",  required=True, type=Path)
    p.add_argument("--cytb-counts",   required=True, type=Path)
    p.add_argument("--cytb-annot",    required=True, type=Path)
    p.add_argument("--metadata",      required=True, type=Path)
    p.add_argument("--outdir",        required=True, type=Path)
    p.add_argument("--season-col",    default="Season")
    p.add_argument("--min-sample-reads", type=int, default=500)
    p.add_argument("--min-relabund",  type=float, default=0.01)
    p.add_argument("--palette",       default="wong")
    args = p.parse_args()

    # Load
    log.info("Loading MiFish counts...")
    mf_counts = load_counts(args.mifish_counts)
    mf_annot  = load_annot(args.mifish_annot)

    log.info("Loading cytb counts...")
    cy_counts = load_counts(args.cytb_counts)
    cy_annot  = load_annot(args.cytb_annot)

    log.info("Loading metadata...")
    season = load_metadata(args.metadata, args.season_col)
    log.info("Season breakdown: %s",
             season.value_counts().to_dict())

    # PA conversion
    log.info("Converting to presence/absence...")
    mf_pa = counts_to_pa(mf_counts,
                         min_sample_reads=args.min_sample_reads,
                         min_relabund=args.min_relabund)
    cy_pa = counts_to_pa(cy_counts,
                         min_sample_reads=50,   # cytb has lower depth
                         min_relabund=args.min_relabund)

    log.info("MiFish: %d taxa x %d samples after PA filter",
             *mf_pa.shape)
    log.info("cytb:   %d taxa x %d samples after PA filter",
             *cy_pa.shape)

    # Collapse to habitat categories
    mf_cat = collapse_to_categories(mf_pa, mf_annot)
    cy_cat = collapse_to_categories(cy_pa, cy_annot)

    # Combine — union of samples, union of categories, OR logic
    if mf_cat.empty and cy_cat.empty:
        log.error("No detections in either marker.")
        return 1

    all_cats = sorted(set(list(mf_cat.index) + list(cy_cat.index)),
                      key=lambda x: CATEGORY_ORDER.index(x)
                      if x in CATEGORY_ORDER else 99)

    # Build combined PA: for each category, a sample is positive if
    # either MiFish OR cytb detected it
    all_samples = sorted(set(
        [shorten_sample_id(c) for c in mf_cat.columns] +
        [shorten_sample_id(c) for c in cy_cat.columns]
    ))

    combined = {}
    for cat in all_cats:
        row = {}
        for s in all_samples:
            mf_val = 0
            cy_val = 0
            # Check MiFish
            if cat in mf_cat.index:
                for col in mf_cat.columns:
                    if shorten_sample_id(col) == s:
                        mf_val = mf_cat.loc[cat, col]
                        break
            # Check cytb
            if cat in cy_cat.index:
                for col in cy_cat.columns:
                    if shorten_sample_id(col) == s:
                        cy_val = cy_cat.loc[cat, col]
                        break
            row[s] = int(mf_val > 0 or cy_val > 0)
        combined[cat] = row

    combined_df = pd.DataFrame(combined).T  # rows=cats, cols=samples

    # Shorten season index too
    season.index = [shorten_sample_id(i) for i in season.index]

    # Detection frequency by season
    freq = compute_detection_freq(combined_df, season)

    # Season sample sizes
    season_ns = season.value_counts().to_dict()

    log.info("Detection frequencies:\n%s", (freq * 100).round(1))

    # Plot
    plot_figure(freq, season_ns, args.outdir, args.palette)
    log.info("Done.")
    return 0


if __name__ == "__main__":
    import matplotlib.ticker
    raise SystemExit(main())
