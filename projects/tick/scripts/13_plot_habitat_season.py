#!/usr/bin/env python3
"""
13_plot_habitat_season.py
=========================
Generate a publication-quality stacked barplot showing fish prey habitat
composition (marine / freshwater / anadromous / estuarine) grouped by
ecological season and optionally colored by mortality group.

This figure serves two purposes:
  1. Seasonal validation — confirms the pipeline detects known loon
     migratory ecology (marine prey in Saltwater season, freshwater
     prey in Breeding season).
  2. Habitat context for dietary differences — shows whether the DvT
     dietary signal reflects habitat-specific prey availability or
     something beyond simple seasonal ecology.

Requires:
  - Cleaned + annotated taxonomy count TSV (from 10b_annotate_diet_ecology.py)
  - QIIME2 metadata TSV with 'Season' and 'Group' columns

ADAPT FOR YOUR STUDY:
  - Change HABITAT_COLORS to match your prey categories
  - Change SEASON_ORDER to match your study species' seasonal ecology
  - The --habitat-col, --season-col, --group-col arguments let you
    rename columns without editing this file

Usage
-----
  # Basic: habitat by season, all samples
  python scripts/12_plot_habitat_season.py \\
      --counts   results/MiFish/all/taxonomy_annotated/taxonomy_counts_annotated_MiFish.tsv \\
      --annot    results/MiFish/all/taxonomy_annotated/annotation_table.tsv \\
      --metadata metadata/qiime/metadata_MiFish.tsv \\
      --marker   MiFish \\
      --outdir   results/MiFish/all/figures/habitat_season/

  # DvT only (exclude Marine birds)
  python scripts/12_plot_habitat_season.py \\
      --counts   results/MiFish/all/taxonomy_annotated/taxonomy_counts_annotated_MiFish.tsv \\
      --annot    results/MiFish/all/taxonomy_annotated/annotation_table.tsv \\
      --metadata metadata/qiime/metadata_MiFish.tsv \\
      --marker   MiFish \\
      --groups   Diseased Trauma \\
      --outdir   results/MiFish/all/figures/habitat_season/

  # Both palettes
  for palette in purple wong; do
    python scripts/12_plot_habitat_season.py \\
        --counts results/MiFish/all/taxonomy_annotated/taxonomy_counts_annotated_MiFish.tsv \\
        --annot  results/MiFish/all/taxonomy_annotated/annotation_table.tsv \\
        --metadata metadata/qiime/metadata_MiFish.tsv \\
        --marker MiFish --palette $palette \\
        --outdir results/MiFish/all/figures/habitat_season/
  done
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

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

# ---------------------------------------------------------------------------
# Habitat group ordering and colors
# Ordered from marine → estuarine → anadromous → freshwater → unclassified
# ADAPT: change these to match your annotation categories and color scheme
# ---------------------------------------------------------------------------
HABITAT_ORDER = [
    "Marine forage fish",
    "Marine bottom fish",
    "Sand lance",
    "Estuarine fish",
    "Estuarine forage fish",
    "Anadromous fish",
    "Marine or anadromous fish",
    "Freshwater prey fish",
    "Freshwater predators",
    "Unclassified prey",
]

# Wong colorblind-safe palette for habitat groups
PALETTES = {
    "wong": {
        "Marine forage fish":      "#0072B2",   # deep blue
        "Marine bottom fish":      "#56B4E9",   # sky blue
        "Sand lance":              "#009E73",   # green
        "Estuarine fish":          "#F0E442",   # yellow
        "Estuarine forage fish":   "#E69F00",   # orange
        "Anadromous fish":         "#CC79A7",   # pink
        "Marine or anadromous fish": "#D55E00", # vermillion
        "Freshwater prey fish":    "#9DC209",   # lime green
        "Freshwater predators":    "#375623",   # dark green
        "Unclassified prey":       "#AAAAAA",   # grey
    },
    "purple": {
        "Marine forage fish":      "#1F4E79",   # dark navy
        "Marine bottom fish":      "#2E75B6",   # mid blue
        "Sand lance":              "#7B2D8B",   # dark purple
        "Estuarine fish":          "#9966CC",   # medium purple
        "Estuarine forage fish":   "#C19FD8",   # lavender
        "Anadromous fish":         "#4B1369",   # deep purple
        "Marine or anadromous fish": "#6B3FA0", # purple
        "Freshwater prey fish":    "#B8D4E8",   # pale blue
        "Freshwater predators":    "#375623",   # dark green
        "Unclassified prey":       "#CCCCCC",   # grey
    }
}

# Ecological season display order and labels
# ADAPT: change these for your study species
SEASON_ORDER   = ["Breeding", "Freshwater_Nonbreeding", "Saltwater"]
SEASON_LABELS  = {
    "Breeding":              "Breeding\n(May–Aug)",
    "Freshwater_Nonbreeding": "Freshwater\nNon-breeding\n(Apr+Sep)",
    "Saltwater":             "Saltwater\n(Oct–Mar)",
}

# Group display colors for the group indicator bar at the bottom
GROUP_COLORS = {
    "Diseased":  "#0072B2",   # Wong blue  — consistent with all other figures
    "Trauma":    "#E69F00",   # Wong orange — consistent with all other figures
    "Marine":    "#009E73",   # Wong green  — already correct
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_counts(counts_path: Path) -> pd.DataFrame:
    """Load annotated taxonomy counts TSV. Index = taxon string."""
    df = pd.read_csv(counts_path, sep="\t", index_col=0)
    return df


def load_annotations(annot_path: Path) -> pd.DataFrame:
    """Load annotation table (Species, habitat, trophic_role, common_group)."""
    df = pd.read_csv(annot_path, sep="\t", index_col=0)
    return df


def load_metadata(meta_path: Path, season_col: str, group_col: str) -> pd.DataFrame:
    """Load QIIME2 metadata, use full sample-id for joining."""
    df = pd.read_csv(meta_path, sep="\t", dtype=str)
    df = df[~df.iloc[:, 0].str.startswith("#", na=False)].reset_index(drop=True)
    # Use the full sample-id as the join key — counts table uses TV220031-GI format
    df["TV"] = df.iloc[:, 0].str.strip()
    cols = ["TV", season_col, group_col]
    return df[[c for c in cols if c in df.columns]].dropna(subset=["TV"])


def build_plot_data(
    counts: pd.DataFrame,
    annot: pd.DataFrame,
    meta: pd.DataFrame,
    habitat_col: str,
    season_col: str,
    group_col: str,
    groups_filter: Optional[List[str]],
    min_reads: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build per-sample relative abundance by habitat group.

    Returns:
        plot_df: DataFrame with samples as columns, habitat groups as rows,
                 values = relative abundance (0-1)
        meta_df: metadata aligned to plot_df columns
    """
    # Sample columns in counts are short TV IDs (TV230084-GI format)
    # Extract TV short ID
    sample_cols = [c for c in counts.columns if c.startswith("TV")]

    # Join counts with habitat annotation
    if habitat_col not in annot.columns:
        raise ValueError(f"Habitat column '{habitat_col}' not found in annotation table. "
                         f"Available: {list(annot.columns)}")

    hab = annot[[habitat_col]].copy()
    counts_hab = counts[sample_cols].join(hab, how="left")
    counts_hab[habitat_col] = counts_hab[habitat_col].fillna("Unclassified prey")

    # Sum reads by habitat group per sample
    grouped = counts_hab.groupby(habitat_col)[sample_cols].sum()

    # Convert to relative abundance
    totals = grouped.sum(axis=0)
    relabund = grouped.div(totals, axis=1).fillna(0)

    # Filter samples by total reads
    keep = totals[totals >= min_reads].index.tolist()
    relabund = relabund[keep]
    totals = totals[keep]

    # Build metadata aligned to samples
    # Sample cols are "TV230084-GI" — use full ID for join
    sample_tv = pd.Series(
        {s: s for s in keep},
        name="TV"
    )
    meta_aligned = sample_tv.reset_index().rename(columns={"index": "sample_id"})
    meta_aligned = meta_aligned.merge(meta, on="TV", how="left")
    meta_aligned = meta_aligned.set_index("sample_id")

    # Filter by groups if requested
    if groups_filter and group_col in meta_aligned.columns:
        keep_samples = meta_aligned[
            meta_aligned[group_col].isin(groups_filter)
        ].index.tolist()
        relabund = relabund[[s for s in keep_samples if s in relabund.columns]]
        meta_aligned = meta_aligned.loc[[s for s in keep_samples if s in meta_aligned.index]]

    # Order samples by season then group
    if season_col in meta_aligned.columns:
        season_map = {s: i for i, s in enumerate(SEASON_ORDER)}
        meta_aligned["_season_order"] = meta_aligned[season_col].map(season_map).fillna(99)
        group_map = {"Diseased": 0, "Trauma": 1, "Marine": 2}
        if group_col in meta_aligned.columns:
            meta_aligned["_group_order"] = meta_aligned[group_col].map(group_map).fillna(99)
            meta_aligned = meta_aligned.sort_values(["_season_order", "_group_order"])
        else:
            meta_aligned = meta_aligned.sort_values("_season_order")
        relabund = relabund[[s for s in meta_aligned.index if s in relabund.columns]]

    return relabund, meta_aligned


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_habitat_season(
    relabund: pd.DataFrame,
    meta: pd.DataFrame,
    palette: Dict[str, str],
    season_col: str,
    group_col: str,
    marker: str,
    title: str,
    outdir: Path,
    output_stem: str,
) -> None:
    """
    Generate the stacked barplot grouped by ecological season.

    Layout:
      - Main panel: stacked bars (relative abundance by habitat group)
      - Bottom strip: mortality group color indicator
      - Season dividers: vertical dashed lines between seasons
      - Season labels: centered below each group
    """
    n_samples = relabund.shape[1]
    fig_width = max(12, n_samples * 0.45)

    fig = plt.figure(figsize=(fig_width, 7), facecolor="white")

    # Main axes + group strip
    ax_main  = fig.add_axes([0.08, 0.18, 0.88, 0.72])
    ax_strip = fig.add_axes([0.08, 0.11, 0.88, 0.05])

    # Order habitat groups
    hab_order = [h for h in HABITAT_ORDER if h in relabund.index]
    # Add any extra groups not in HABITAT_ORDER
    extras = [h for h in relabund.index if h not in hab_order]
    hab_order = hab_order + extras

    # Stack bars
    bottoms = np.zeros(n_samples)
    x = np.arange(n_samples)
    bar_width = 0.85

    for hab in hab_order:
        if hab not in relabund.index:
            continue
        vals = relabund.loc[hab].values
        color = palette.get(hab, "#CCCCCC")
        ax_main.bar(x, vals, bottom=bottoms, width=bar_width,
                    color=color, linewidth=0, label=hab)
        bottoms += vals

    # Season dividers and labels
    seasons_seen = []
    season_starts = {}
    season_ends = {}

    samples = list(relabund.columns)
    for i, samp in enumerate(samples):
        season = meta.loc[samp, season_col] if samp in meta.index else ""
        if season not in season_starts:
            season_starts[season] = i
        season_ends[season] = i
        if season and season not in seasons_seen:
            seasons_seen.append(season)

    for season in seasons_seen[1:]:
        div_x = season_starts[season] - 0.5
        ax_main.axvline(div_x, color="#333333", linewidth=1.2,
                        linestyle="--", alpha=0.6, zorder=5)
        ax_strip.axvline(div_x, color="#333333", linewidth=1.2,
                         linestyle="--", alpha=0.6, zorder=5)

    # Season labels below strip
    for season in seasons_seen:
        mid = (season_starts[season] + season_ends[season]) / 2
        label = SEASON_LABELS.get(season, season)
        ax_strip.text(mid, -1.2, label,
                      ha="center", va="top", fontsize=8.5,
                      fontfamily="Arial", color="#222222")

    # Group indicator strip
    for i, samp in enumerate(samples):
        group = meta.loc[samp, group_col] if (samp in meta.index and group_col in meta.columns) else ""
        color = GROUP_COLORS.get(group, "#DDDDDD")
        ax_strip.barh(0, 1, left=i - 0.5, height=1,
                      color=color, linewidth=0)

    ax_strip.set_xlim(-0.5, n_samples - 0.5)
    ax_strip.set_ylim(-0.5, 0.5)
    ax_strip.axis("off")

    # Main axes formatting
    ax_main.set_xlim(-0.5, n_samples - 0.5)
    ax_main.set_ylim(0, 1.0)
    ax_main.set_xticks([])
    ax_main.set_ylabel("Relative abundance", fontsize=11, fontfamily="Arial")
    # Title removed for publication — add via figure legend or caption
    # ax_main.set_title(title, fontsize=12, fontweight="bold",
    #                   fontfamily="Arial", pad=10)
    ax_main.spines["top"].set_visible(False)
    ax_main.spines["right"].set_visible(False)
    ax_main.tick_params(axis="y", labelsize=9)

    # Legend — habitat groups
    legend_handles = []
    for hab in hab_order:
        if hab in relabund.index and relabund.loc[hab].sum() > 0:
            color = palette.get(hab, "#CCCCCC")
            legend_handles.append(
                mpatches.Patch(facecolor=color, label=hab, linewidth=0)
            )
    # Group legend
    for grp, color in GROUP_COLORS.items():
        legend_handles.append(
            mpatches.Patch(facecolor=color, label=f"● {grp}", linewidth=0)
        )

    ax_main.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize=7.5,
        framealpha=0.9,
        edgecolor="#CCCCCC",
        ncol=2,
        handlelength=1.2,
        handleheight=0.9,
    )

    # Save
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        path = outdir / f"{output_stem}.{ext}"
        fig.savefig(path, dpi=300 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
        log.info("Saved: %s", path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="12_plot_habitat_season.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    req = p.add_argument_group("required")
    req.add_argument("--counts",   required=True, type=Path,
                     help="Annotated taxonomy counts TSV (from 10b_annotate_diet_ecology.py)")
    req.add_argument("--annot",    required=True, type=Path,
                     help="Annotation table TSV (from 10b_annotate_diet_ecology.py)")
    req.add_argument("--metadata", required=True, type=Path,
                     help="QIIME2 metadata TSV with Season and Group columns")
    req.add_argument("--marker",   required=True,
                     help="Marker name for output labeling (e.g. MiFish)")
    req.add_argument("--outdir",   required=True, type=Path,
                     help="Output directory for figures")

    opt = p.add_argument_group("optional")
    opt.add_argument("--palette",     default="wong", choices=list(PALETTES.keys()),
                     help="Color palette. Default: wong")
    opt.add_argument("--groups",      nargs="+", default=None,
                     help="Restrict to these Group values (e.g. Diseased Trauma). "
                          "Default: all groups")
    opt.add_argument("--habitat-col", default="common_group",
                     help="Annotation column to use for habitat grouping. "
                          "Default: common_group")
    opt.add_argument("--season-col",  default="Season",
                     help="Metadata column for ecological season. Default: Season")
    opt.add_argument("--group-col",   default="Group",
                     help="Metadata column for mortality group. Default: Group")
    opt.add_argument("--min-reads",   type=int, default=500,
                     help="Minimum total reads per sample. Default: 500")
    opt.add_argument("--output-stem", default=None,
                     help="Output filename stem. Default: auto-generated")
    opt.add_argument("--title",       default=None,
                     help="Figure title. Default: auto-generated")
    opt.add_argument("--exclude-unclassified", action="store_true", default=False,
                     help="Exclude 'Unclassified prey' from the figure and renormalize "
                          "to 100%% of identified taxa only. Generates a cleaner figure "
                          "but note in caption that unclassified reads are excluded.")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    for path, name in [(args.counts, "--counts"), (args.annot, "--annot"),
                       (args.metadata, "--metadata")]:
        if not path.exists():
            log.error("%s not found: %s", name, path)
            return 1

    log.info("Loading counts: %s", args.counts)
    counts = load_counts(args.counts)

    log.info("Loading annotations: %s", args.annot)
    annot = load_annotations(args.annot)

    log.info("Loading metadata: %s", args.metadata)
    meta = load_metadata(args.metadata, args.season_col, args.group_col)

    try:
        relabund, meta_aligned = build_plot_data(
            counts, annot, meta,
            habitat_col=args.habitat_col,
            season_col=args.season_col,
            group_col=args.group_col,
            groups_filter=args.groups,
            min_reads=args.min_reads,
        )
    except ValueError as e:
        log.error("%s", e)
        return 1

    if relabund.empty:
        log.error("No samples remaining after filters. Check --groups and --min-reads.")
        return 1

    # Optionally exclude unclassified prey and renormalize
    if args.exclude_unclassified:
        unclass_rows = [r for r in relabund.index
                        if "Unclassified" in r or "unclassified" in r.lower()]
        if unclass_rows:
            log.info("Excluding %d unclassified row(s): %s", len(unclass_rows), unclass_rows)
            relabund = relabund.drop(index=unclass_rows)
        # Renormalize to identified reads only
        totals = relabund.sum(axis=0)
        # Drop samples with no identified taxa at all
        has_signal = totals[totals > 0].index.tolist()
        if len(has_signal) < relabund.shape[1]:
            log.warning("Dropping %d sample(s) with no identified taxa after exclusion",
                        relabund.shape[1] - len(has_signal))
        relabund = relabund[has_signal]
        meta_aligned = meta_aligned.loc[
            [s for s in has_signal if s in meta_aligned.index]
        ]
        totals = relabund.sum(axis=0)
        relabund = relabund.div(totals, axis=1).fillna(0)
        log.info("Renormalized to identified taxa only (%d samples with signal)",
                 len(has_signal))

    # Drop samples with no season assigned (missing date)
    if args.season_col in meta_aligned.columns:
        no_season = meta_aligned[
            meta_aligned[args.season_col].isna() |
            (meta_aligned[args.season_col] == "")
        ].index.tolist()
        if no_season:
            log.warning("Dropping %d sample(s) with no season assigned: %s",
                        len(no_season), no_season)
            keep_s = [s for s in relabund.columns if s not in no_season]
            relabund = relabund[keep_s]
            meta_aligned = meta_aligned.loc[[s for s in keep_s if s in meta_aligned.index]]

    log.info("Plotting %d samples, %d habitat groups",
             relabund.shape[1], relabund.shape[0])

    palette = PALETTES[args.palette]
    groups_str = "_".join(args.groups) if args.groups else "all"
    excl_tag = "_identified_only" if args.exclude_unclassified else ""
    stem = args.output_stem or f"{args.marker}_habitat_season_{groups_str}{excl_tag}_{args.palette}"
    excl_note = " (identified taxa only)" if args.exclude_unclassified else ""
    title = args.title or (
        f"{args.marker} fish prey habitat composition by ecological season{excl_note}"
        + (f" — {', '.join(args.groups)}" if args.groups else "")
    )

    plot_habitat_season(
        relabund=relabund,
        meta=meta_aligned,
        palette=palette,
        season_col=args.season_col,
        group_col=args.group_col,
        marker=args.marker,
        title=title,
        outdir=args.outdir,
        output_stem=stem,
    )

    log.info("=== Done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
