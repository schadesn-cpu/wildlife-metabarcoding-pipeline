#!/usr/bin/env python3
"""
plot_mifish_ecology_pa.py
=========================
Two-panel stacked detection frequency figure for MiFish presence/absence data.

Panel A: Detection frequency by ecological season
Panel B: Detection frequency by COD group

Each bar = proportion of samples in that group where each common_group category
was detected (at least one species from that category). Bars colored by
common_group using the Wong 2011 colorblind-safe palette.

Usage:
    python plot_mifish_ecology_pa.py \
        --pa       results/MiFish/all/presence_absence/presence_absence_MiFish.tsv \
        --annot    results/MiFish/all/taxonomy_annotated/annotation_table.tsv \
        --meta     metadata/qiime/metadata_MiFish.tsv \
        --cod-meta metadata/qiime/metadata_MiFish_cod.tsv \
        --outdir   results/MiFish/all/figures/ecology_pa/
"""

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

# Preferred group order and color assignment
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

# Strip colors for the group-identity bar drawn below the x-axis
SEASON_STRIP_COLORS = {
    "Breeding":               "#0072B2",  # wong blue
    "Freshwater_Nonbreeding": "#E69F00",  # wong orange
    "Saltwater":              "#009E73",  # wong green
}

COD_STRIP_COLORS = {
    "Lead":                  "#0072B2",  # wong blue  (Diseased subset)
    "Parasitic_Infectious":  "#CC79A7",  # wong purple (Diseased subset)
    "Trauma":                "#E69F00",  # wong orange
    "Marine":                "#009E73",  # wong green
}


def detection_freq_by_group(pa: pd.DataFrame,
                             annot: pd.DataFrame,
                             grouping: pd.Series,
                             group_order: list) -> pd.DataFrame:
    """
    For each grouping level and each common_group category, compute the
    proportion of samples where at least one species from that category
    was detected.

    Returns DataFrame: rows = common_group, cols = grouping levels.
    Values = detection frequency (0-1).
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
            # Taxa belonging to this common_group
            taxa_in_cat = annot[annot["common_group"] == cat].index.tolist()
            taxa_present = [t for t in taxa_in_cat if t in pa_sub.index]
            if not taxa_present:
                freqs[cat] = 0.0
            else:
                # Detection = at least one taxon in category detected in sample
                detected = pa_sub.loc[taxa_present].sum(axis=0) > 0
                freqs[cat] = detected.mean()
        results[grp] = pd.Series(freqs)

    return pd.DataFrame(results)


def make_stacked_bar(ax, freq_df: pd.DataFrame,
                     title: str, n_per_group: dict,
                     strip_colors: dict = None) -> None:
    """
    Draw stacked horizontal bar chart.
    freq_df: rows = common_group categories, cols = grouping levels.
    strip_colors: optional dict mapping group name → color for the
                  identity strip drawn below the x-axis.
    """
    groups = freq_df.columns.tolist()
    categories = [g for g in GROUP_ORDER if freq_df.loc[g].sum() > 0]

    x = np.arange(len(groups))
    bar_width = 0.55

    bottoms = np.zeros(len(groups))
    handles = []
    for cat in categories:
        vals = freq_df.loc[cat].values.astype(float)
        bars = ax.bar(x, vals * 100, bar_width,
                      bottom=bottoms,
                      color=GROUP_COLORS[cat],
                      label=cat,
                      edgecolor="white", linewidth=0.5)
        bottoms += vals * 100
        handles.append(mpatches.Patch(color=GROUP_COLORS[cat], label=cat))

    # ── Group identity strip below x-axis ────────────────────────────────────
    STRIP_H   = 4.5   # height in data-% units
    STRIP_GAP = 2.0   # gap between strip top and the zero line
    if strip_colors:
        for i, grp in enumerate(groups):
            color = strip_colors.get(grp, "#999999")
            ax.bar(i, STRIP_H, bar_width,
                   bottom=-(STRIP_H + STRIP_GAP),
                   color=color, clip_on=False, zorder=3,
                   edgecolor="white", linewidth=0.5)
        ax.set_ylim(-(STRIP_H + STRIP_GAP + 0.5), 105)
        # Keep the spine anchored at y=0 so it doesn't drop into the strip
        ax.spines["bottom"].set_position(("data", 0))
    else:
        ax.set_ylim(0, 105)

    # X-axis labels with n
    labels = [f"{g}\n(n={n_per_group.get(g, 0)})" for g in groups]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Detection frequency (%)", fontsize=10)
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0f}%" if v >= 0 else "")
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5, zorder=0)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)

    return handles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
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
    season_col = "Season"
    season = meta.loc[meta.index.isin(samples), season_col].dropna()
    season = season[season != ""]
    season_order = ["Breeding", "Freshwater_Nonbreeding", "Saltwater"]
    season_n = season.value_counts().to_dict()

    log.info("Season breakdown: %s", season_n)

    # ── COD grouping ──────────────────────────────────────────────────────────
    cod_col = "COD_broad"
    # Combine main metadata Group (for Marine) with COD metadata
    cod = cod_meta.loc[cod_meta.index.isin(samples), cod_col] if cod_col in cod_meta.columns else pd.Series()
    # Fill in Marine from Group column
    group = meta.loc[meta.index.isin(samples), "Group"]
    cod_combined = cod.copy()
    for s in samples:
        if s not in cod_combined.index or pd.isna(cod_combined.get(s, None)):
            if group.get(s) == "Marine":
                cod_combined[s] = "Marine"
    # Order by prevalence, keeping meaningful groups
    cod_order_pref = ["Lead", "Parasitic_Infectious", "Trauma", "Marine"]
    cod_order = [g for g in cod_order_pref if g in cod_combined.values]
    cod_n = cod_combined.value_counts().to_dict()

    log.info("COD breakdown: %s", cod_n)

    # ── Compute detection frequencies ────────────────────────────────────────
    log.info("Computing detection frequencies...")
    season_freq = detection_freq_by_group(pa, annot, season, season_order)
    cod_freq    = detection_freq_by_group(pa, annot, cod_combined, cod_order)

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    plt.subplots_adjust(wspace=0.35)

    handles_a = make_stacked_bar(
        axes[0], season_freq,
        "A   Detection frequency by ecological season",
        season_n,
        strip_colors=SEASON_STRIP_COLORS)

    handles_b = make_stacked_bar(
        axes[1], cod_freq,
        "B   Detection frequency by cause of death",
        cod_n,
        strip_colors=COD_STRIP_COLORS)

    # Shared legend — use handles from whichever panel has more categories
    all_handles = handles_a if len(handles_a) >= len(handles_b) else handles_b
    fig.legend(
        handles=all_handles,
        title="Prey category",
        title_fontsize=9,
        fontsize=8.5,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.12),
        framealpha=0.9,
    )

    fig.suptitle(
        "MiFish 12S \u2014 Dietary detection frequency by ecological group\n"
        "(classified prey only, \u2265500 reads/sample, \u22651% relative abundance)",
        fontsize=11, y=1.02
    )

    # ── Save ──────────────────────────────────────────────────────────────────
    for ext in ("png", "svg"):
        out = outdir / f"MiFish_ecology_pa_detection.{ext}"
        fig.savefig(out, dpi=180, bbox_inches="tight")
        log.info("Saved: %s", out)

    plt.close(fig)
    log.info("Done.")


if __name__ == "__main__":
    import matplotlib.ticker
    main()
