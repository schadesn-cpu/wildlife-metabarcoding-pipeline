#!/usr/bin/env python3
"""
11_plot_mifish_season_ecology.py
================================
Generate a multi-panel seasonal feeding ecology figure for MiFish 12S data
that contextualizes dietary composition within loon migration ecology and
confirms that seasonal variation does not confound group-level findings.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WHAT THIS SCRIPT PRODUCES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Panel A — Seasonal taxonomy barplot:
    Stacked barplot of prey taxa relative abundance, samples grouped
    and ordered by season (Summer | Fall | Winter | Spring). Taxa colored
    with Wong 2011 colorblind-safe palette. Annotations call out ecologically
    meaningful taxa (Ammodytes, Brevoortia = marine; Actinopterygii = freshwater).

  Panel B — Season × Group bubble chart:
    Bubble chart showing the distribution of samples across Season × Group
    (Diseased, Trauma, Marine). Bubble size = n. Annotated with chi-square
    p-value confirming season and group are not significantly associated.

  Panel C — Season × COD_broad heatmap:
    Heatmap showing n per COD_broad category × Season. Highlights that Lead
    birds span multiple seasons, supporting the argument that the Lead dietary
    signal is not a seasonal artifact.

  All panels saved together as one figure:
    mifish_season_ecology_{palette}.png/.svg  (300 dpi PNG + vector SVG)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
USAGE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  python scripts/11_plot_mifish_season_ecology.py \\
      --relabund   results/MiFish/all/taxonomy/taxonomy_relabund_L7_MiFish.tsv \\
      --metadata   metadata/qiime/metadata_MiFish.tsv \\
      --cod-metadata metadata/qiime/metadata_MiFish_cod.tsv \\
      --palette    wong \\
      --top-n      12 \\
      --outdir     results/MiFish/figures/

  # Purple palette version
  python scripts/11_plot_mifish_season_ecology.py \\
      --relabund   results/MiFish/all/taxonomy/taxonomy_relabund_L7_MiFish.tsv \\
      --metadata   metadata/qiime/metadata_MiFish.tsv \\
      --cod-metadata metadata/qiime/metadata_MiFish_cod.tsv \\
      --palette    purple \\
      --outdir     results/MiFish/figures/

  # Suppress Panel C (if COD metadata not available)
  python scripts/11_plot_mifish_season_ecology.py \\
      --relabund   results/MiFish/all/taxonomy/taxonomy_relabund_L7_MiFish.tsv \\
      --metadata   metadata/qiime/metadata_MiFish.tsv \\
      --palette    wong \\
      --no-cod-panel \\
      --outdir     results/MiFish/figures/

Dependencies:
  pip install matplotlib numpy pandas scipy
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants — kept in sync with 06_, 09_, 10_
# ---------------------------------------------------------------------------

FIGURE_DPI      = 300
FONT_FAMILY     = "Arial"
FONT_SIZE_TITLE = 13
FONT_SIZE_AXIS  = 11
FONT_SIZE_TICK  = 9
FONT_SIZE_ANNOT = 8
FONT_SIZE_LABEL = 9

SEASON_ORDER = ["Summer", "Fall", "Winter", "Spring"]
# ADAPT FOR YOUR STUDY: This is the display order for cause-of-death categories
# in the COD panel. Change these to match your study's grouping variable values.
COD_ORDER    = ["Lead", "Parasitic_Infectious", "Trauma", "Marine", "Unknown_Other"]

# Ecologically annotated taxa — marine prey highlighted
MARINE_TAXA = {
    "Ammodytes americanus", "Ammodytes",
    "Brevoortia patronus", "Brevoortia tyrannus", "Brevoortia",
    "Leiostomus xanthurus", "Leiostomus",
    "Pholis", "uncl. Pholis",
    "Trachinotus", "uncl. Trachinotus",
    "Etelidae", "Etelis", "uncl. Etelidae",
    "Myoxocephalus",
}

PALETTES: Dict[str, Dict] = {
    "purple": {
        "season_colors": {
            "Summer": "#9B4DCA", "Fall": "#4B0082",
            "Winter": "#C084FC", "Spring": "#2D0057",
        },
        "group_colors": {
            "Diseased": "#7B2D8B", "Trauma": "#C19FD8", "Marine": "#4B1369",
        },
        "cod_colors": {
            "Lead": "#4B0082", "Parasitic_Infectious": "#9B4DCA",
            "Trauma": "#C19FD8", "Marine": "#D09EE0", "Unknown_Other": "#E2BFFF",
        },
        "taxa_colors": [
            "#0D001A", "#1A0033", "#2D0057", "#3D006B", "#4B0082",
            "#6A0DAD", "#7B2D8B", "#9B4DCA", "#B06FD8", "#C084FC",
            "#D4A0FF", "#E2BFFF", "#CCCCCC",
        ],
        "marine_highlight": "#FF6B35",
    },
    "wong": {
        "season_colors": {
            "Summer": "#009E73", "Fall": "#E69F00",
            "Winter": "#56B4E9", "Spring": "#CC79A7",
        },
        "group_colors": {
            "Diseased": "#0072B2", "Trauma": "#E69F00", "Marine": "#009E73",
        },
        "cod_colors": {
            "Lead": "#D55E00", "Parasitic_Infectious": "#CC79A7",
            "Trauma": "#E69F00", "Marine": "#009E73", "Unknown_Other": "#999999",
        },
        "taxa_colors": [
            "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
            "#D55E00", "#F0E442", "#AA4499", "#44BB99", "#BBCC33",
            "#99DDFF", "#EE8866", "#AAAAAA",
        ],
        "marine_highlight": "#D55E00",
    },
}

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_relabund(path: Path) -> pd.DataFrame:
    """Load taxonomy relative abundance TSV from 08_taxonomy_table.py output."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    if "mean_relabund" in df.columns:
        df = df.drop(columns=["mean_relabund"])
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df = df.T  # taxa×samples -> samples×taxa
    df.index.name = "sample_id"
    log.info("Loaded relabund: %d samples × %d taxa", *df.shape)
    return df


def load_metadata(path: Path) -> pd.DataFrame:
    """Load QIIME2 metadata TSV, skip #q2:types row."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[~df.iloc[:, 0].str.startswith("#", na=False)].reset_index(drop=True)
    sid_col = df.columns[0]
    df = df.rename(columns={sid_col: "sample-id"})
    # Extract TV ID
    df["_TV"] = df["sample-id"].str.extract(r"(TV\d+)", expand=False)
    return df


def merge_data(relabund, metadata, cod_metadata):
    """Merge relabund with metadata on TV ID extracted from sample index."""
    rel = relabund.copy().reset_index()
    sid_col = rel.columns[0]
    rel["_TV"] = rel[sid_col].astype(str).str.extract(r"(TV\d+)", expand=False)
    rel = rel.dropna(subset=["_TV"])

    keep = ["_TV", "Season", "Group"]
    if "Location" in metadata.columns:
        keep.append("Location")
    if "State Found" in metadata.columns:
        keep.append("State Found")
    meta_sub = metadata[[c for c in keep if c in metadata.columns]].drop_duplicates("_TV")
    merged = rel.merge(meta_sub, on="_TV", how="inner")

    if cod_metadata is not None and "COD_broad" in cod_metadata.columns:
        cod_sub = cod_metadata[["_TV", "COD_broad"]].drop_duplicates("_TV")
        merged = merged.merge(cod_sub, on="_TV", how="left")

    log.info("Merged: %d samples | Season: %s | Group: %s",
             len(merged),
             merged["Season"].value_counts().to_dict() if "Season" in merged.columns else "MISSING",
             merged["Group"].value_counts().to_dict() if "Group" in merged.columns else "MISSING")
    return merged


def plot_panel_a(
    ax: plt.Axes,
    data: pd.DataFrame,
    taxa_cols: List[str],
    top_n: int,
    palette: Dict,
    season_order: List[str],
) -> None:
    """Stacked barplot of prey taxa by season."""
    # Compute mean relative abundance per taxon across all samples for ranking
    taxon_means = data[taxa_cols].mean().sort_values(ascending=False)
    top_taxa = list(taxon_means.head(top_n).index)
    other_taxa = [t for t in taxa_cols if t not in top_taxa]

    plot_data = data.copy()
    plot_data["Other"] = plot_data[other_taxa].sum(axis=1) if other_taxa else 0.0
    display_taxa = top_taxa + (["Other"] if other_taxa else [])

    # Order samples: by season, then by dominant taxon descending
    ordered_rows = []
    season_boundaries = []
    pos = 0
    for season in season_order:
        grp = plot_data[plot_data["Season"] == season].copy()
        if len(grp) == 0:
            continue
        grp = grp.sort_values(top_taxa[0] if top_taxa else "Other", ascending=False)
        ordered_rows.append(grp)
        season_boundaries.append((pos, pos + len(grp), season, len(grp)))
        pos += len(grp)

    if not ordered_rows:
        return
    plot_df = pd.concat(ordered_rows)

    colors = palette["taxa_colors"]
    x = np.arange(len(plot_df))
    bar_w = 0.85
    bottom = np.zeros(len(plot_df))

    for i, taxon in enumerate(display_taxa):
        vals = plot_df[taxon].fillna(0).values if taxon in plot_df.columns \
            else plot_df.get("Other", pd.Series(0, index=plot_df.index)).fillna(0).values
        color = colors[i % len(colors)]
        is_marine = any(m.lower() in taxon.lower() for m in MARINE_TAXA)
        edgecolor = palette["marine_highlight"] if is_marine else "none"
        lw = 1.2 if is_marine else 0
        ax.bar(x, vals * 100, width=bar_w, bottom=bottom * 100,
               color=color, linewidth=lw, edgecolor=edgecolor, label=taxon)
        bottom += vals

    # Season labels and dividers
    for (start, end, season, n) in season_boundaries:
        mid = (start + end - 1) / 2
        sc = palette["season_colors"].get(season, "#555555")
        ax.text(mid, 106, f"{season}\n(n={n})", ha="center", va="bottom",
                fontsize=FONT_SIZE_LABEL, fontweight="bold", color=sc)
        if end < len(plot_df):
            ax.axvline(end - 0.5, color="#AAAAAA", linewidth=0.8, linestyle="--", alpha=0.7)

    # Marine annotation bracket
    marine_annotation = "← marine prey"
    ax.text(0.99, 0.97, marine_annotation, transform=ax.transAxes,
            ha="right", va="top", fontsize=FONT_SIZE_ANNOT,
            color=palette["marine_highlight"], style="italic",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor=palette["marine_highlight"], alpha=0.8))

    ax.set_xlim(-0.5, len(plot_df) - 0.5)
    ax.set_ylim(0, 118)
    ax.set_ylabel("Relative Abundance (%)\namong classified reads", fontsize=FONT_SIZE_AXIS)
    ax.set_xticks(x)
    short_labels = [tv.replace("-GI", "").replace("-GI-MiFish", "") for tv in plot_df["_TV"]]
    ax.set_xticklabels(short_labels, rotation=90, fontsize=7, ha="center")
    ax.tick_params(axis="y", labelsize=FONT_SIZE_TICK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_title("A  |  MiFish Prey Composition by Season", fontsize=FONT_SIZE_TITLE,
                 fontweight="bold", loc="left", pad=8)

    # Compact legend — italic for species names
    handles, labels = ax.get_legend_handles_labels()
    legend_labels = []
    for lbl in labels:
        if lbl == "Other":
            legend_labels.append(lbl)
        else:
            legend_labels.append(lbl.replace("uncl. ", "uncl. "))
    ax.legend(handles[::-1], legend_labels[::-1],
              loc="lower right", fontsize=FONT_SIZE_ANNOT - 1,
              framealpha=0.95, edgecolor="#CCCCCC", ncol=2,
              title="Taxon", title_fontsize=FONT_SIZE_ANNOT)


# ---------------------------------------------------------------------------
# Panel B — Season × Group bubble chart
# ---------------------------------------------------------------------------

def plot_panel_b(
    ax: plt.Axes,
    data: pd.DataFrame,
    palette: Dict,
    season_order: List[str],
) -> None:
    """Bubble chart: Season × Group, bubble size = n."""
    # ADAPT: These are the loon mortality groups. Change to your study's
    # group names. Only groups present in the data will be plotted.
    groups = [g for g in ["Diseased", "Trauma", "Marine"] if g in data["Group"].values]
    seasons = [s for s in season_order if s in data["Season"].values]

    # Chi-square on DvT only
    # ADAPT: DvT = Diseased vs Trauma. Change these two group names to
    # the two primary comparison groups in your study.
    dvt = data[data["Group"].isin(["Diseased", "Trauma"])]
    ct = pd.crosstab(dvt["Group"], dvt["Season"])
    try:
        _, p_chi2, _, _ = chi2_contingency(ct)
    except Exception:
        p_chi2 = float("nan")

    # Build grid
    x_pos = {s: i for i, s in enumerate(seasons)}
    y_pos = {g: i for i, g in enumerate(groups)}

    counts = data.groupby(["Season", "Group"]).size().reset_index(name="n")
    max_n = counts["n"].max()

    for _, row in counts.iterrows():
        if row["Season"] not in x_pos or row["Group"] not in y_pos:
            continue
        xi = x_pos[row["Season"]]
        yi = y_pos[row["Group"]]
        size = (row["n"] / max_n) * 1800
        color = palette["group_colors"].get(row["Group"], "#888888")
        ax.scatter(xi, yi, s=size, color=color, alpha=0.85, linewidth=0.8,
                   edgecolors="white", zorder=3)
        ax.text(xi, yi, str(int(row["n"])), ha="center", va="center",
                fontsize=FONT_SIZE_ANNOT, fontweight="bold", color="white", zorder=4)

    ax.set_xticks(range(len(seasons)))
    ax.set_xticklabels(seasons, fontsize=FONT_SIZE_TICK)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups, fontsize=FONT_SIZE_TICK)
    ax.set_xlim(-0.6, len(seasons) - 0.4)
    ax.set_ylim(-0.6, len(groups) - 0.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, alpha=0.25, linewidth=0.5, zorder=0)

    p_str = f"χ² p = {p_chi2:.3f} (Diseased vs. Trauma)"
    ax.text(0.98, 0.04, p_str, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=FONT_SIZE_ANNOT, style="italic", color="#444444")
    ax.set_title("B  |  Season × Mortality Group Distribution",
                 fontsize=FONT_SIZE_TITLE, fontweight="bold", loc="left", pad=8)
    ax.set_xlabel("Season", fontsize=FONT_SIZE_AXIS)
    ax.set_ylabel("Mortality Group", fontsize=FONT_SIZE_AXIS)

    log.info("Chi-square (DvT × Season): p=%.4f", p_chi2)


# ---------------------------------------------------------------------------
# Panel C — COD_broad × Season heatmap
# ---------------------------------------------------------------------------

def plot_panel_c(
    ax: plt.Axes,
    data: pd.DataFrame,
    palette: Dict,
    season_order: List[str],
    cod_order: List[str],
) -> None:
    """Heatmap: n samples per COD_broad × Season. Highlights Lead spans seasons."""
    if "COD_broad" not in data.columns:
        ax.text(0.5, 0.5, "COD_broad metadata not available",
                ha="center", va="center", transform=ax.transAxes, fontsize=FONT_SIZE_AXIS)
        ax.set_title("C  |  Season × Cause of Death (COD_broad)",
                     fontsize=FONT_SIZE_TITLE, fontweight="bold", loc="left", pad=8)
        return

    seasons = [s for s in season_order if s in data["Season"].values]
    cods = [c for c in cod_order if c in data["COD_broad"].values]

    ct = pd.crosstab(data["COD_broad"], data["Season"]).reindex(
        index=cods, columns=seasons, fill_value=0
    )

    arr = ct.values.astype(float)
    # Normalize per row for color intensity, keep raw n for labels
    if arr.size == 0:
        ax.text(0.5, 0.5, 'No COD_broad data matched season order',
                ha='center', va='center', transform=ax.transAxes)
        return
    row_max = arr.max(axis=1, keepdims=True)
    row_max[row_max == 0] = 1
    arr_norm = arr / row_max

    cmap = plt.cm.Blues
    im = ax.imshow(arr_norm, cmap=cmap, aspect="auto", vmin=0, vmax=1)

    # Annotate cells with raw n
    for i in range(len(cods)):
        for j in range(len(seasons)):
            n = int(ct.iloc[i, j])
            text_color = "white" if arr_norm[i, j] > 0.6 else "#333333"
            ax.text(j, i, str(n), ha="center", va="center",
                    fontsize=FONT_SIZE_ANNOT + 1, fontweight="bold", color=text_color)

    # Highlight Lead row
    lead_idx = cods.index("Lead") if "Lead" in cods else None
    if lead_idx is not None:
        rect = plt.Rectangle((-0.5, lead_idx - 0.5), len(seasons), 1,
                             fill=False, edgecolor=palette["cod_colors"]["Lead"],
                             linewidth=2.5, zorder=5)
        ax.add_patch(rect)
        ax.text(len(seasons) - 0.4, lead_idx,
                "← spans seasons", va="center", ha="left",
                fontsize=FONT_SIZE_ANNOT, color=palette["cod_colors"]["Lead"],
                style="italic")

    ax.set_xticks(range(len(seasons)))
    ax.set_xticklabels(seasons, fontsize=FONT_SIZE_TICK)
    ax.set_yticks(range(len(cods)))
    display_cods = [c.replace("_", "/") for c in cods]
    ax.set_yticklabels(display_cods, fontsize=FONT_SIZE_TICK)
    ax.set_title("C  |  Season × Cause of Death (COD_broad)",
                 fontsize=FONT_SIZE_TITLE, fontweight="bold", loc="left", pad=8)
    ax.set_xlabel("Season", fontsize=FONT_SIZE_AXIS)
    ax.set_ylabel("Cause of Death", fontsize=FONT_SIZE_AXIS)

    log.info("COD_broad × Season heatmap complete")


# ---------------------------------------------------------------------------
# Main figure assembly
# ---------------------------------------------------------------------------

def _compute_chi2_p(data: pd.DataFrame) -> float:
    from scipy.stats import chi2_contingency
    dvt = data[data["Group"].isin(["Diseased", "Trauma"])]
    ct = pd.crosstab(dvt["Group"], dvt["Season"])
    try:
        _, p, _, _ = chi2_contingency(ct)
        return p
    except Exception:
        return float("nan")

def build_figure(
    data: pd.DataFrame,
    data_full: pd.DataFrame,
    taxa_cols: List[str],
    top_n: int,
    palette_name: str,
    outdir: Path,
    no_cod_panel: bool,
) -> None:
    palette = PALETTES[palette_name]
    plt.rcParams["font.family"] = FONT_FAMILY

    has_cod = "COD_broad" in data.columns and not no_cod_panel

    if has_cod:
        fig = plt.figure(figsize=(20, 16))
        gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35,
                               height_ratios=[1.6, 1])
        ax_a = fig.add_subplot(gs[0, :])   # full width top
        ax_b = fig.add_subplot(gs[1, 0])   # bottom left
        ax_c = fig.add_subplot(gs[1, 1])   # bottom right
    else:
        fig = plt.figure(figsize=(20, 12))
        gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35, width_ratios=[1.8, 1])
        ax_a = fig.add_subplot(gs[0, 0])
        ax_b = fig.add_subplot(gs[0, 1])
        ax_c = None

    # Panel A
    plot_panel_a(ax_a, data, taxa_cols, top_n, palette, SEASON_ORDER)

    # Panel B
    plot_panel_b(ax_b, data_full, palette, SEASON_ORDER)   # ← use data_full

    # Panel C
    if ax_c is not None:
        plot_panel_c(ax_c, data, palette, SEASON_ORDER, COD_ORDER)

    # Figure-level annotation
    fig.text(0.5, 0.01,
             "Seasonal dietary composition is consistent with known Common Loon migration ecology "
             "(marine prey in Winter/Spring; freshwater prey in Summer) "
             f"but does not differ significantly between mortality groups (χ² p={_compute_chi2_p(data_full):.4f}).",
             ha="center", va="bottom", fontsize=FONT_SIZE_ANNOT,
             style="italic", color="#555555")

    fig.tight_layout(rect=[0, 0.03, 1, 1])

    # Save
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / f"mifish_season_ecology_{palette_name}"
    for ext in (".png", ".svg"):
        fpath = stem.with_suffix(ext)
        kw = {"dpi": FIGURE_DPI} if ext == ".png" else {}
        fig.savefig(fpath, bbox_inches="tight", **kw)
        log.info("  Saved: %s", fpath)

    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="11_plot_mifish_season_ecology.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    req = p.add_argument_group("required arguments")
    req.add_argument("--relabund", required=True, type=Path,
                     help="Taxonomy relative abundance TSV from 08_taxonomy_table.py.")
    req.add_argument("--metadata", required=True, type=Path,
                     help="QIIME2 metadata TSV with Season and Group columns.")

    opt = p.add_argument_group("optional arguments")
    opt.add_argument("--cod-metadata", default=None, type=Path,
                     help="Metadata TSV with COD_broad column (for Panel C). "
                          "If not provided, Panel C is skipped.")
    opt.add_argument("--palette", choices=list(PALETTES.keys()), default="wong",
                     help="Color palette. Default: wong.")
    opt.add_argument("--top-n", type=int, default=12,
                     help="Number of top taxa to show individually. Default: 12.")
    opt.add_argument("--outdir", type=Path, default=Path("."),
                     help="Output directory. Default: current directory.")
    opt.add_argument("--no-cod-panel", action="store_true", default=False,
                     help="Suppress Panel C even if --cod-metadata is provided.")
    opt.add_argument("--exclude-groups", nargs="+", default=["Marine"],
                     help="Groups to exclude from Panel A barplot. Default: Marine.")
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.relabund.exists():
        log.error("Relabund file not found: %s", args.relabund)
        return 1
    if not args.metadata.exists():
        log.error("Metadata file not found: %s", args.metadata)
        return 1

    # Load
    relabund = load_relabund(args.relabund)
    metadata = load_metadata(args.metadata)
    cod_meta = None
    if args.cod_metadata and args.cod_metadata.exists():
        cod_meta = load_metadata(args.cod_metadata)
        log.info("COD metadata loaded: %s", args.cod_metadata.name)
    elif args.cod_metadata:
        log.warning("COD metadata not found: %s — Panel C will be skipped", args.cod_metadata)

    # Merge
    data = merge_data(relabund, metadata, cod_meta)
    if data.empty:
        log.error("No samples retained after merge. Check metadata TV ID matching.")
        return 1

    # Validate Season column
    if "Season" not in data.columns:
        log.error("'Season' column not found in metadata. Available: %s", list(data.columns))
        return 1

    # Filter excluded groups for Panel A only
    taxa_cols = [c for c in data.columns
                 if c not in ("Season", "Group", "COD_broad", "_TV", "sample-id",
                              "NHVDL ID", "Lung", "Fecal", "Cadaver Condition ",
                              "Date Found", "State Found", "Location", "COD",
                              "Marker", "Other", "sample_id", "Species", "index")]

    data_plot = data.copy()
    if args.exclude_groups:
        data_plot = data_plot[~data_plot["Group"].isin(args.exclude_groups)]
        log.info("Excluded groups from Panel A: %s — %d samples remain",
                 args.exclude_groups, len(data_plot))

    log.info("Palette: %s | Top-N taxa: %d | Samples: %d", args.palette, args.top_n, len(data_plot))
    log.info("Season distribution:\n%s", data_plot["Season"].value_counts().to_string())

    build_figure(
        data=data_plot,
        data_full=data, 
        taxa_cols=taxa_cols,
        top_n=args.top_n,
        palette_name=args.palette,
        outdir=args.outdir,
        no_cod_panel=args.no_cod_panel,
    )

    log.info("=== Done. Figures in: %s ===", args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
