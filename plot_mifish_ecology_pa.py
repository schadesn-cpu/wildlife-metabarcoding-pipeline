#!/usr/bin/env python3
"""
plot_mifish_ecology_pa.py
=========================
Two-panel grouped-bar detection-frequency figure for MiFish presence/absence
data.

Panel A: Detection frequency by ecological season
Panel B: Detection frequency by cause of death

Each bar = proportion of samples in that grouping level where at least one
species from a given common_group prey category was detected. Because a
single loon can be positive for multiple prey categories simultaneously
(marine forage fish AND freshwater prey fish, etc.), these per-category
proportions are NOT mutually exclusive and therefore cannot be stacked.
Bars are grouped side-by-side within each season / COD level.

Bars colored by common_group using the Wong 2011 colorblind-safe palette.

Usage:
    python plot_mifish_ecology_pa.py \\
        --pa       results/MiFish/all/presence_absence/presence_absence_MiFish.tsv \\
        --annot    results/MiFish/all/taxonomy_annotated/annotation_table.tsv \\
        --meta     metadata/qiime/metadata_MiFish.tsv \\
        --cod-meta metadata/qiime/metadata_MiFish_cod.tsv \\
        --outdir   results/MiFish/all/figures/ecology_pa/
"""

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S")
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
    "#000000",  # black
    "#999999",  # grey
]

# Preferred prey-category order. Colors assigned positionally from WONG.
GROUP_ORDER = [
    "Marine bottom fish",
    "Marine forage fish",
    "Anadromous fish",
    "Estuarine fish",
    "Estuarine forage fish",
    "Freshwater prey fish",
    "Freshwater predators",
    "Sand lance",
    "Unclassified prey",
]
GROUP_COLORS = {g: WONG[i % len(WONG)] for i, g in enumerate(GROUP_ORDER)}


def detection_freq_by_group(pa: pd.DataFrame,
                            annot: pd.DataFrame,
                            grouping: pd.Series,
                            group_order: list) -> pd.DataFrame:
    """
    For each grouping level and each common_group category, compute the
    proportion of samples in which at least one species from that category
    was detected.

    Returns DataFrame: rows = common_group (GROUP_ORDER), cols = grouping levels.
    Values are independent per category (NOT mutually exclusive) and will
    exceed 1.0 if summed across categories.
    """
    results = {}
    for grp in group_order:
        grp_samples = grouping[grouping == grp].index.tolist()
        grp_samples = [s for s in grp_samples if s in pa.columns]
        if not grp_samples:
            results[grp] = pd.Series(0.0, index=GROUP_ORDER)
            continue

        pa_sub = pa[grp_samples]
        freqs = {}
        for cat in GROUP_ORDER:
            taxa_in_cat = annot[annot["common_group"] == cat].index.tolist()
            taxa_present = [t for t in taxa_in_cat if t in pa_sub.index]
            if not taxa_present:
                freqs[cat] = 0.0
            else:
                detected = pa_sub.loc[taxa_present].sum(axis=0) > 0
                freqs[cat] = detected.mean()
        results[grp] = pd.Series(freqs)

    return pd.DataFrame(results)


def make_grouped_bar(ax, freq_df: pd.DataFrame,
                     title: str, n_per_group: dict) -> list:
    """
    Draw grouped (side-by-side) bar chart.

    freq_df:  rows = prey categories, cols = grouping levels.
    Each prey category gets its own colored bar within each grouping level.
    Y-axis is detection frequency (%) capped at 0-105.

    Returns a list of legend handles (one per plotted category).
    """
    groups = freq_df.columns.tolist()
    categories = [g for g in GROUP_ORDER if freq_df.loc[g].sum() > 0]

    if not categories:
        ax.text(0.5, 0.5, "No detections",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=12, color="#666666")
        ax.set_title(title, fontsize=11, fontweight="bold", loc="left", pad=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        return []

    n_groups = len(groups)
    n_cats   = len(categories)

    cluster_width = 0.85
    bar_width = cluster_width / n_cats
    x = np.arange(n_groups)

    handles = []
    for i, cat in enumerate(categories):
        offsets = x + (i - (n_cats - 1) / 2) * bar_width
        vals = freq_df.loc[cat].values.astype(float) * 100
        ax.bar(offsets, vals, bar_width,
               color=GROUP_COLORS[cat], label=cat,
               edgecolor="white", linewidth=0.5, zorder=3)
        handles.append(mpatches.Patch(color=GROUP_COLORS[cat], label=cat))

    ax.set_ylim(0, 105)
    labels = [f"{g}\n(n={n_per_group.get(g, 0)})" for g in groups]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Detection frequency (%)", fontsize=10)
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0f}%")
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5, zorder=0)
    ax.set_title(title, fontsize=11, fontweight="bold", loc="left", pad=8)

    return handles


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pa",       required=True)
    parser.add_argument("--annot",    required=True)
    parser.add_argument("--meta",     required=True)
    parser.add_argument("--cod-meta", required=True, dest="cod_meta")
    parser.add_argument("--outdir",   required=True)
    parser.add_argument("--palette",  default="wong")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────────────
    log.info("Loading PA table...")
    pa = pd.read_csv(args.pa, sep="\t", index_col=0)
    log.info("  %d taxa x %d samples", *pa.shape)

    log.info("Loading annotation table...")
    annot = pd.read_csv(args.annot, sep="\t", index_col=0)

    log.info("Loading metadata...")
    meta     = pd.read_csv(args.meta,     sep="\t", index_col=0, skiprows=[1])
    cod_meta = pd.read_csv(args.cod_meta, sep="\t", index_col=0, skiprows=[1])

    samples = pa.columns.tolist()

    # ── Season grouping ───────────────────────────────────────────────────────
    season_col   = "Season"
    season       = meta.loc[meta.index.isin(samples), season_col].dropna()
    season       = season[season != ""]
    season_order = ["Breeding", "Freshwater_Nonbreeding", "Saltwater"]
    season_n     = season.value_counts().to_dict()
    log.info("Season breakdown: %s", season_n)

    # ── COD grouping ──────────────────────────────────────────────────────────
    cod_col = "COD_broad"
    cod = (cod_meta.loc[cod_meta.index.isin(samples), cod_col]
           if cod_col in cod_meta.columns else pd.Series(dtype=str))
    group = meta.loc[meta.index.isin(samples), "Group"]
    cod_combined = cod.copy()
    for s in samples:
        if s not in cod_combined.index or pd.isna(cod_combined.get(s, None)):
            if group.get(s) == "Marine":
                cod_combined[s] = "Marine"

    cod_order_pref = ["Lead", "Parasitic_Infectious", "Trauma", "Marine"]
    cod_order      = [g for g in cod_order_pref if g in cod_combined.values]
    cod_n          = cod_combined.value_counts().to_dict()
    log.info("COD breakdown: %s", cod_n)

    # ── Compute detection frequencies ────────────────────────────────────────
    log.info("Computing detection frequencies...")
    season_freq = detection_freq_by_group(pa, annot, season,       season_order)
    cod_freq    = detection_freq_by_group(pa, annot, cod_combined, cod_order)

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    plt.subplots_adjust(wspace=0.28)

    make_grouped_bar(axes[0], season_freq,
                     "A   Detection frequency by ecological season",
                     season_n)
    make_grouped_bar(axes[1], cod_freq,
                     "B   Detection frequency by cause of death",
                     cod_n)

    # Shared legend — union of categories seen across both panels, GROUP_ORDER sorted
    seen = set()
    union_cats = []
    for cat in GROUP_ORDER:
        has_data = ((cat in season_freq.index and season_freq.loc[cat].sum() > 0) or
                    (cat in cod_freq.index    and cod_freq.loc[cat].sum() > 0))
        if has_data and cat not in seen:
            union_cats.append(cat)
            seen.add(cat)
    union_handles = [mpatches.Patch(color=GROUP_COLORS[c], label=c) for c in union_cats]

    fig.legend(
        handles=union_handles,
        title="Prey category",
        title_fontsize=9,
        fontsize=8.5,
        loc="lower center",
        ncol=min(len(union_handles), 5),
        bbox_to_anchor=(0.5, -0.10),
        framealpha=0.9,
    )

    fig.suptitle(
        "MiFish 12S \u2014 Dietary detection frequency by ecological group\n"
        "(presence/absence, \u2265500 reads/sample, \u22651% relative abundance)",
        fontsize=11, y=1.02,
    )

    # ── Save ──────────────────────────────────────────────────────────────────
    for ext in ("png", "svg"):
        out = outdir / f"MiFish_ecology_pa_detection.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        log.info("Saved: %s", out)

    plt.close(fig)
    log.info("Done.")


if __name__ == "__main__":
    main()
