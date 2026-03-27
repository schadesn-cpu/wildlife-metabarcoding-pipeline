#!/usr/bin/env python3
"""
09_plot_taxonomy.py
===================
Generate publication-quality stacked taxonomy barplots from the relative
abundance TSV produced by 08_taxonomy_table.py.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PIPELINE POSITION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  08_taxonomy_table.py  →  taxonomy_relabund_L{N}_{marker}.tsv
                                          ↓
                          09_plot_taxonomy.py   ← metadata TSV
                                          ↓
                          barplot_{marker}_{group}_{palette}.png/.svg

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WHAT THIS SCRIPT PRODUCES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  - Stacked barplot of relative abundance (%) per sample
  - Top-N taxa shown individually, remainder collapsed to "Other"
  - Samples ordered within groups by dominant taxon
  - Group labels above bars with n= counts
  - Vertical dashed dividers between groups
  - Italic taxon names in legend (non-genus labels like "Other" left upright)
  - Both PNG (300 dpi, slides/SharePoint) and SVG (Illustrator-editable)

Relative abundance is calculated among classified reads only (output of 08_).
Methods note: "Relative abundance was calculated among reads assigned at
[level] level following [database] classification."

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PALETTES — kept in sync with 06_plot_diversity.py and 07_visualize_diversity.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  --palette purple   Dark-to-light purple gradient for taxa bars.
                     Group headers colored by palette group colors.
                     Best for single-marker manuscripts.

  --palette redblue  Red-to-blue gradient for taxa bars.
                     Best for multi-marker panels or colorblind reviewers.

  --palette wong     Wong 2011 8-color colorblind-safe palette for taxa bars.
                     Use when reviewers request a citable colorblind palette.
                     Suitable for ≤8 taxa; extra taxa fall back to grey.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
USAGE EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # 16S Diseased vs Trauma, purple palette
  python 09_plot_taxonomy.py \\
      --relabund  results/16S/DvT/taxonomy/taxonomy_relabund_L6_16S.tsv \\
      --metadata  metadata/qiime/metadata_16S.tsv \\
      --group-by  Group \\
      --marker    16S \\
      --palette   purple \\
      --outdir    results/16S/DvT/figures/taxonomy/

  # Same data, red-blue palette
  python 09_plot_taxonomy.py \\
      --relabund  results/16S/DvT/taxonomy/taxonomy_relabund_L6_16S.tsv \\
      --metadata  metadata/qiime/metadata_16S.tsv \\
      --group-by  Group \\
      --marker    16S \\
      --palette   redblue \\
      --outdir    results/16S/DvT/figures/taxonomy/

  # Seasonal grouping, explicit order, top 20 taxa
  python 09_plot_taxonomy.py \\
      --relabund    results/16S/DvT/taxonomy/taxonomy_relabund_L6_16S.tsv \\
      --metadata    metadata/qiime/metadata_16S.tsv \\
      --group-by    Season \\
      --group-order Spring Summer Fall Winter \\
      --top-n       20 \\
      --marker      16S \\
      --palette     purple \\
      --outdir      results/16S/DvT/figures/taxonomy/

  # MiFish species level
  python 09_plot_taxonomy.py \\
      --relabund  results/MiFish/all/taxonomy/taxonomy_relabund_L7_MiFish.tsv \\
      --metadata  metadata/qiime/metadata_MiFish.tsv \\
      --group-by  Group \\
      --marker    MiFish \\
      --palette   redblue \\
      --outdir    results/MiFish/all/figures/taxonomy/

  # List metadata columns available for --group-by
  python 09_plot_taxonomy.py \\
      --relabund  results/16S/DvT/taxonomy/taxonomy_relabund_L6_16S.tsv \\
      --metadata  metadata/qiime/metadata_16S.tsv \\
      --list-columns

Dependencies:
  pip install matplotlib numpy pandas
"""
from __future__ import annotations

import argparse
import io
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Style constants — match 06_ / 07_
# ---------------------------------------------------------------------------

FIGURE_DPI    = 300
FONT_SIZE_TITLE  = 13
FONT_SIZE_AXIS   = 11
FONT_SIZE_TICK   = 9
FONT_SIZE_LEGEND = 8.5
FONT_SIZE_GROUP  = 11

# ---------------------------------------------------------------------------
# Palettes
# ---------------------------------------------------------------------------
# Group colors: used for group header labels above bars.
# Kept in sync with 06_plot_diversity.py and 07_visualize_diversity.py.
#
# Taxa colors: used for the stacked bar segments.
# 16-step gradients ensure visual separation for up to 15 named taxa + Other.
# The final color in each list is always "Other".

PALETTES: Dict[str, Dict] = {
    "purple": {
        # Group header label colors (up to 8 groups)
        "group_colors": [
            "#7B2D8B",  # dark purple   — Diseased / Group 1
            "#C19FD8",  # lavender      — Trauma / Group 2
            "#4B1369",  # deep purple   — Group 3
            "#D09EE0",  # light purple  — Group 4
            "#2D0A40", "#E0BAEC", "#9870B0", "#F0D6F5",
        ],
        # 16-step dark-to-light purple gradient for taxa bars
        "taxa_colors": [
            "#0D001A", "#1A0033", "#2D0057", "#3D006B", "#4B0082",
            "#6A0DAD", "#7B2D8B", "#9B4DCA", "#B06FD8", "#C084FC",
            "#D4A0FF", "#E2BFFF", "#EEDAFF", "#F7EEFF", "#FDF6FF",
            "#CCCCCC",  # Other
        ],
    },
    "redblue": {
        "group_colors": [
            "#B22222",  # red           — Diseased / Group 1
            "#2E86C1",  # blue          — Trauma / Group 2
            "#E74C3C",  # light red     — Group 3
            "#7FB3D3",  # light blue    — Group 4
            "#7B241C", "#1A5276", "#F1948A", "#2980B9",
        ],
        "taxa_colors": [
            "#7B0000", "#A00000", "#C0392B", "#E74C3C", "#F1948A",
            "#1A3A5C", "#1E5F8A", "#2471A3", "#2E86C1", "#7FB3D3",
            "#AED6F1", "#D6EAF8", "#1A5276", "#117A65", "#1E8449",
            "#AAAAAA",  # Other
        ],
    },
    "wong": {
        # Wong 2011 + Paul Tol 'Muted' extension — both colorblind-safe and citable
        # Group header label colors
        "group_colors": [
            "#0072B2",  # blue          — Diseased / Group 1
            "#E69F00",  # orange        — Trauma / Group 2
            "#009E73",  # green         — Marine / Group 3
            "#CC79A7",  # pink          — Group 4
            "#56B4E9",  # sky blue      — Group 5
            "#D55E00",  # vermillion    — Group 6
            "#F0E442",  # yellow        — Group 7
            "#999999",  # grey          — Group 8
        ],
        # 15-taxon colorblind-safe palette:
        # Colors 1-7: Wong 2011 (no black — replaced with Tol Indigo)
        # Colors 8-15: Paul Tol 'Muted' qualitative palette (doi:10.5281/zenodo.3381072)
        "taxa_colors": [
            "#0072B2",  # Wong blue         — taxon 1  (most abundant)
            "#E69F00",  # Wong orange       — taxon 2
            "#009E73",  # Wong green        — taxon 3
            "#CC79A7",  # Wong pink         — taxon 4
            "#56B4E9",  # Wong sky blue     — taxon 5
            "#D55E00",  # Wong vermillion   — taxon 6
            "#F0E442",  # Wong yellow       — taxon 7
            "#AA4499",  # Tol indigo        — taxon 8  (replaces black)
            "#44BB99",  # Tol teal          — taxon 9
            "#BBCC33",  # Tol olive         — taxon 10
            "#99DDFF",  # Tol pale blue     — taxon 11
            "#EE8866",  # Tol salmon        — taxon 12
            "#FFAABB",  # Tol rose          — taxon 13
            "#DDDDDD",  # Tol pale grey     — taxon 14
            "#994F00",  # Tol brown         — taxon 15
            "#AAAAAA",  # Other
        ],
    },
}

# ---------------------------------------------------------------------------
# Metadata loading — matches pattern used in 06_ and 07_
# ---------------------------------------------------------------------------

def load_metadata(metadata_path: Path, group_column: str) -> pd.Series:
    """
    Load a QIIME 2 metadata TSV or standard CSV/TSV.
    Returns a Series {sample_id: group_value}.

    Handles:
    - QIIME 2 TSV (#q2:types row, first col = sample-id)
    - Standard TSV/CSV with TV column
    - BOM characters from Excel exports
    """
    with metadata_path.open(encoding="utf-8-sig") as f:
        lines = f.readlines()

    # Strip optional #q2:types line
    header = lines[0]
    data_lines = [l for l in lines[1:] if not l.startswith("#")]
    content = header + "".join(data_lines)

    sep = "\t" if "\t" in header else ","
    df = pd.read_csv(io.StringIO(content), sep=sep, index_col=0, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df.index = df.index.str.strip()

    if group_column not in df.columns:
        log.error(
            "Column '%s' not found in metadata. Available columns: %s",
            group_column, list(df.columns),
        )
        sys.exit(1)

    series = df[group_column].dropna()
    series.name = group_column
    return series


def list_metadata_columns(metadata_path: Path) -> None:
    """Print available columns and their unique values."""
    with metadata_path.open(encoding="utf-8-sig") as f:
        lines = f.readlines()
    header = lines[0]
    data_lines = [l for l in lines[1:] if not l.startswith("#")]
    content = header + "".join(data_lines)
    sep = "\t" if "\t" in header else ","
    df = pd.read_csv(io.StringIO(content), sep=sep, index_col=0, dtype=str)
    df.columns = [c.strip() for c in df.columns]

    print(f"\nMetadata: {metadata_path}")
    print(f"Samples : {len(df)}\n")
    print(f"{'Column':<30}  {'Unique values'}")
    print("-" * 70)
    for col in df.columns:
        unique = sorted(df[col].dropna().unique(), key=str)
        if len(unique) <= 10:
            print(f"  {col:<28}  {unique}")
        else:
            print(f"  {col:<28}  ({len(unique)} unique values)")


# ---------------------------------------------------------------------------
# Sample ID matching — same TV#### pattern as 06_ / 07_
# ---------------------------------------------------------------------------

def match_samples(
    relabund_cols: pd.Index,
    group_series: pd.Series,
) -> Dict[str, Optional[str]]:
    """
    Map each sample ID in the relabund table to a group label.

    Handles the common pattern where relabund columns look like
    'TV230007-GI-16S_S1483' but metadata uses 'TV230007-GI-16S_S1483'
    (direct match from QIIME metadata) or 'TV230007' (short key).
    """
    mapping: Dict[str, Optional[str]] = {}
    meta_ids = group_series.index.tolist()

    for sid in relabund_cols:
        # 1. Direct match
        if sid in group_series.index:
            mapping[sid] = group_series[sid]
            continue
        # 2. Metadata ID is a prefix of the sample ID (TV230007 in TV230007-GI-16S_S1483)
        matched = [m for m in meta_ids if sid.startswith(m)]
        if matched:
            best = max(matched, key=len)
            mapping[sid] = group_series[best]
            continue
        # 3. Sample ID is contained in the metadata ID (reverse of above)
        matched = [m for m in meta_ids if m.startswith(sid)]
        if matched:
            best = max(matched, key=len)
            mapping[sid] = group_series[best]
            continue
        # No match
        mapping[sid] = None

    unmatched = [k for k, v in mapping.items() if v is None]
    if unmatched:
        log.warning(
            "%d sample(s) in relabund table not matched to metadata "
            "(will be excluded from plot): %s%s",
            len(unmatched),
            unmatched[:5],
            " ..." if len(unmatched) > 5 else "",
        )

    return mapping


# ---------------------------------------------------------------------------
# Relabund loading and top-N collapsing
# ---------------------------------------------------------------------------

def load_relabund(tsv_path: Path) -> pd.DataFrame:
    """
    Load a taxonomy_relabund_*.tsv produced by 08_taxonomy_table.py.
    Returns DataFrame: taxa (rows) × samples (cols), values 0–1.
    """
    df = pd.read_csv(tsv_path, sep="\t", index_col=0)

    # Drop the mean_relabund column if present (from top-N table)
    if "mean_relabund" in df.columns:
        df = df.drop(columns=["mean_relabund"])

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    log.info(
        "Loaded relabund table: %d taxa × %d samples from %s",
        len(df), len(df.columns), tsv_path.name,
    )
    return df


def collapse_to_top_n(
    relabund: pd.DataFrame,
    top_n: int,
    samples: Optional[List[str]] = None,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Keep the top_n taxa by mean relative abundance across `samples`
    (or all samples if None). Remaining taxa are summed into 'Other'.

    Returns (collapsed_df, display_taxa_list) where display_taxa_list
    includes 'Other' at the end.
    """
    subset = relabund[samples] if samples else relabund
    means = subset.mean(axis=1).sort_values(ascending=False)
    top_taxa = means.head(top_n).index.tolist()

    result = relabund.reindex(top_taxa).copy()
    other_mask = ~relabund.index.isin(top_taxa)
    if other_mask.any():
        result.loc["Other"] = relabund.loc[other_mask].sum(axis=0)
    else:
        result.loc["Other"] = 0.0

    # Re-normalize to 100% per sample (relabund values are 0-1, convert to %)
    result = result * 100

    display_taxa = top_taxa + ["Other"]
    return result, display_taxa


# ---------------------------------------------------------------------------
# Barplot
# ---------------------------------------------------------------------------

def _is_plain_name(name: str) -> bool:
    """
    Returns True if the taxon name should NOT be italicized in the legend.
    Non-italic: Other, uncl.* labels, fully unclassified rows.
    """
    if name == "Other":
        return True
    if name.startswith("uncl."):
        return True
    if name.lower() in ("unclassified", "uncultured bacterium"):
        return True
    return False


def _sort_group_samples(
    samples: List[str],
    relabund_pct: pd.DataFrame,
    display_taxa: List[str],
) -> List[str]:
    """Sort samples within a group by the abundance of the dominant taxon."""
    top_non_other = [t for t in display_taxa if t != "Other"]
    if not top_non_other:
        return samples
    dominant = top_non_other[0]
    if dominant in relabund_pct.index:
        order = relabund_pct.loc[dominant, samples].sort_values(ascending=False)
        return list(order.index)
    return samples


def plot_barplot(
    relabund: pd.DataFrame,
    sample_groups: Dict[str, Optional[str]],
    group_order: List[str],
    top_n: int,
    palette_name: str,
    marker: str,
    group_column: str,
    title: Optional[str],
    outpath_stem: Path,
    show_title: bool = True,
) -> None:
    """
    Generate and save a stacked barplot (PNG + SVG).

    Parameters
    ----------
    relabund      : taxa × samples DataFrame (values 0–1, from 08_)
    sample_groups : {sample_id: group_label} mapping
    group_order   : ordered list of group labels to plot
    top_n         : number of top taxa to show before collapsing to Other
    palette_name  : 'purple', 'redblue', or 'wong'
    marker        : marker label for figure title (e.g. '16S')
    group_column  : metadata column name used for grouping (for display)
    title         : override figure title (None = auto-generate)
    outpath_stem  : output path without extension (will save .png and .svg)
    """
    palette = PALETTES[palette_name]
    group_colors = palette["group_colors"]
    taxa_colors  = palette["taxa_colors"]

    # ── Build ordered sample list ─────────────────────────────────────────
    ordered_samples: List[str] = []
    group_spans: List[Tuple[str, int, int, str]] = []  # (label, start_idx, end_idx, color)

    for gi, grp in enumerate(group_order):
        grp_samples = [s for s in relabund.columns
                       if sample_groups.get(s) == grp]
        if not grp_samples:
            log.warning("Group '%s' has no samples in the relabund table — skipping.", grp)
            continue

        # Collapse top-N using only this group's samples for ranking
        rel_pct_all, display_taxa = collapse_to_top_n(relabund, top_n)
        sorted_samps = _sort_group_samples(grp_samples, rel_pct_all, display_taxa)

        start = len(ordered_samples)
        ordered_samples.extend(sorted_samps)
        end = len(ordered_samples) - 1
        color = group_colors[gi % len(group_colors)]
        group_spans.append((grp, start, end, color))

    if not ordered_samples:
        log.error("No samples matched any group in %s. Check --group-by and metadata.", group_order)
        sys.exit(1)

    # ── Final top-N collapse across all ordered samples ───────────────────
    rel_pct, display_taxa = collapse_to_top_n(relabund, top_n, ordered_samples)

    # Assign taxa colors
    tax_color_map = {
        t: taxa_colors[i] if i < len(taxa_colors) else "#CCCCCC"
        for i, t in enumerate(display_taxa)
    }

    # ── Figure ────────────────────────────────────────────────────────────
    n = len(ordered_samples)
    fig_w = max(14, n * 0.62)
    fig, ax = plt.subplots(figsize=(fig_w, 6.5))
    fig.patch.set_facecolor("white")

    xs = np.arange(n)
    bottom = np.zeros(n)

    for taxon in display_taxa:
        if taxon not in rel_pct.index:
            continue
        vals = rel_pct.reindex([taxon])[ordered_samples].fillna(0).values[0]
        ax.bar(
            xs, vals, bottom=bottom,
            color=tax_color_map[taxon],
            width=0.85, linewidth=0, zorder=2,
        )
        bottom += vals

    # ── Group dividers and header labels ──────────────────────────────────
    for i, (grp, start, end, color) in enumerate(group_spans):
        mid = (start + end) / 2
        n_grp = end - start + 1
        ax.text(
            mid, 101.5,
            f"{grp}  (n={n_grp})",
            ha="center", va="bottom",
            fontsize=FONT_SIZE_GROUP, fontweight="bold", color=color,
        )
        if i > 0:
            ax.axvline(
                start - 0.5,
                color="#999999", lw=1.2, linestyle="--", zorder=5,
            )

    # ── X-axis tick labels ────────────────────────────────────────────────
    short_labels = [
        re.sub(r"-GI-.*", "", s).replace("TV", "") for s in ordered_samples
    ]
    ax.set_xticks(xs)
    ax.set_xticklabels(short_labels, rotation=45, ha="right", fontsize=FONT_SIZE_TICK)

    # ── Axes formatting ───────────────────────────────────────────────────
    level_label = "Relative Abundance (%)\namong classified reads"
    ax.set_ylabel(level_label, fontsize=FONT_SIZE_AXIS)
    ax.set_ylim(0, 108)
    ax.set_xlim(-0.6, n - 0.4)
    ax.yaxis.grid(True, alpha=0.3, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for sp in ["top", "right", "bottom"]:
        ax.spines[sp].set_visible(False)

    # ── Title ─────────────────────────────────────────────────────────────
    if show_title:
        fig_title = title or (
            f"Genus-Level Relative Abundance — {marker}\n"
            f"Grouped by {group_column}"
        )
        ax.set_title(fig_title, fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=10)

    # ── Legend ────────────────────────────────────────────────────────────
    legend_handles = [
        mpatches.Patch(
            facecolor=tax_color_map[t],
            edgecolor="#aaaaaa",
            linewidth=0.3,
            label=t,
        )
        for t in reversed(display_taxa)
        if t in tax_color_map
    ]
    leg = ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        fontsize=FONT_SIZE_LEGEND,
        framealpha=0.9,
        edgecolor="#cccccc",
        title="Taxon",
        title_fontsize=FONT_SIZE_LEGEND,
    )
    for text in leg.get_texts():
        if not _is_plain_name(text.get_text()):
            text.set_fontstyle("italic")

    # ── Save ──────────────────────────────────────────────────────────────
    plt.tight_layout()
    outpath_stem.parent.mkdir(parents=True, exist_ok=True)

    for ext in (".png", ".svg"):
        fpath = outpath_stem.with_suffix(ext)
        fig.savefig(
            fpath, dpi=FIGURE_DPI if ext == ".png" else None,
            bbox_inches="tight", facecolor="white",
            format="svg" if ext == ".svg" else None,
        )
        log.info("  Saved: %s", fpath)

    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """
Build and return the argument parser for 09_plot_taxonomy.py.

    Key arguments: --relabund (TSV from 08_taxonomy_table.py), --marker,
    --group-by (metadata column), --palette, --top-n, --outdir.
    A --list-columns flag is also available to inspect available metadata
    columns without producing any figures.
    """
    p = argparse.ArgumentParser(
        prog="09_plot_taxonomy.py",
        description=(
            "Generate stacked taxonomy barplots from 08_taxonomy_table.py output.\n"
            "Produces PNG (300 dpi) and SVG per group × palette combination."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    req = p.add_argument_group("required arguments")
    req.add_argument(
        "--relabund", required=True, type=Path,
        help=(
            "Relative abundance TSV produced by 08_taxonomy_table.py "
            "(taxonomy_relabund_L{N}_{marker}.tsv). "
            "Rows = taxa, columns = samples, values = 0–1."
        ),
    )
    req.add_argument(
        "--metadata", required=True, type=Path,
        help=(
            "QIIME 2 metadata TSV or source metadata CSV/TSV. "
            "Must contain a sample-id column (or TV column) and the --group-by column."
        ),
    )
    req.add_argument(
        "--group-by", dest="group_by", default=None,
        help="Metadata column to group samples by (e.g. Group, Season).",
    )

    opt = p.add_argument_group("plot options")
    opt.add_argument(
        "--group-order", dest="group_order", nargs="*", default=None,
        help=(
            "Explicit group order left-to-right (e.g. --group-order Diseased Trauma). "
            "Default: alphabetical."
        ),
    )
    opt.add_argument(
        "--top-n", dest="top_n", type=int, default=15,
        help="Number of top taxa to show individually. Rest collapsed to 'Other'. Default: 15.",
    )
    opt.add_argument(
        "--palette",
        choices=list(PALETTES.keys()),
        default="purple",
        help=(
            "Color palette. "
            "purple = dark-to-light purple gradient (single-marker manuscript). "
            "redblue = red-to-blue gradient (multi-marker or colorblind). "
            "wong = Wong 2011 citable 8-color colorblind-safe palette. "
            "Default: purple."
        ),
    )
    opt.add_argument(
        "--marker", default="",
        help="Marker label for figure titles and output filenames (e.g. 16S, MiFish).",
    )
    opt.add_argument(
        "--title", default=None,
        help="Override the auto-generated figure title.",
    )
    opt.add_argument(
        "--outdir", type=Path, default=None,
        help=(
            "Output directory. "
            "Default: same directory as --relabund."
        ),
    )
    opt.add_argument(
        "--no-title", action="store_true", dest="no_title",
        help="Suppress the figure title (use for journal submission figures).",
    )
    opt.add_argument(
        "--output-stem", default=None, dest="output_stem",
        help=(
            "Base filename stem (no extension). "
            "Default: barplot_{marker}_{group_by}_{palette}. "
            "Example: --output-stem 16S_unrarefied_DvT_barplot_purple"
        ),
    )
    opt.add_argument(
        "--suffix", default="",
        help="Optional extra suffix added to output filenames.",
    )

    util = p.add_argument_group("utility")
    util.add_argument(
        "--list-columns", action="store_true",
        help="Print available metadata columns and exit (requires --metadata).",
    )
    util.add_argument(
        "--dry-run", action="store_true",
        help="Print planned outputs without generating figures.",
    )

    return p


def main(argv=None) -> int:
    """
Parse arguments and generate stacked taxonomy barplots.

    Reads the relative abundance TSV produced by 08_taxonomy_table.py, groups
    samples by the specified metadata column, collapses low-abundance taxa into
    'Other', and writes one PNG (300 dpi) + SVG pair per palette to --outdir.
    Returns 0 on success, 2 if required input files are missing, 1 on error.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    # ── Utility: list columns ─────────────────────────────────────────────
    if args.list_columns:
        if not args.metadata or not Path(args.metadata).exists():
            log.error("--list-columns requires --metadata.")
            return 2
        list_metadata_columns(Path(args.metadata))
        return 0

    # ── Validate required args ────────────────────────────────────────────
    if not args.group_by:
        log.error("--group-by is required. Use --list-columns to see available columns.")
        return 2
    if not args.relabund.exists():
        log.error("Relabund file not found: %s", args.relabund)
        return 2
    if not args.metadata.exists():
        log.error("Metadata file not found: %s", args.metadata)
        return 2

    # ── Load data ─────────────────────────────────────────────────────────
    relabund = load_relabund(args.relabund)
    group_series = load_metadata(args.metadata, args.group_by)

    # ── Match samples to groups ───────────────────────────────────────────
    sample_groups = match_samples(relabund.columns, group_series)

    matched = {k: v for k, v in sample_groups.items() if v is not None}
    if not matched:
        log.error(
            "No samples in the relabund table matched any sample ID in the metadata.\n"
            "  Relabund columns (first 3): %s\n"
            "  Metadata IDs (first 3): %s",
            list(relabund.columns[:3]),
            list(group_series.index[:3]),
        )
        return 1

    log.info("Matched %d / %d samples to metadata", len(matched), len(relabund.columns))

    # ── Determine group order ─────────────────────────────────────────────
    all_groups = sorted({v for v in sample_groups.values() if v is not None})
    if args.group_order:
        group_order = args.group_order
        missing = [g for g in group_order if g not in all_groups]
        if missing:
            log.warning("Groups in --group-order not found in metadata: %s", missing)
    else:
        group_order = all_groups

    log.info("Groups: %s", group_order)
    for g in group_order:
        n = sum(1 for v in sample_groups.values() if v == g)
        log.info("  %s: n=%d", g, n)

    # ── Build output path stem ────────────────────────────────────────────
    outdir = args.outdir or args.relabund.parent
    if args.output_stem:
        stem_name = args.output_stem
    else:
        marker_str   = f"_{args.marker}"  if args.marker   else ""
        group_str    = f"_{args.group_by}"
        suffix_str   = f"_{args.suffix}"  if args.suffix   else ""
        palette_str  = f"_{args.palette}"
        stem_name    = f"barplot{marker_str}{group_str}{suffix_str}{palette_str}"
    outpath_stem = outdir / stem_name

    log.info("Output stem : %s", outpath_stem)
    log.info("Top-N taxa  : %d", args.top_n)
    log.info("Palette     : %s", args.palette)

    if args.dry_run:
        log.info("DRY RUN — no files written.")
        for ext in (".png", ".svg"):
            log.info("  Would write: %s", outpath_stem.with_suffix(ext))
        return 0

    # ── Plot ──────────────────────────────────────────────────────────────
    plot_barplot(
        relabund     = relabund,
        sample_groups= sample_groups,
        group_order  = group_order,
        top_n        = args.top_n,
        palette_name = args.palette,
        marker       = args.marker,
        group_column = args.group_by,
        title        = args.title,
        outpath_stem = outpath_stem,
        show_title   = not args.no_title,
    )

    log.info("=== Done. Figures in: %s ===", outdir)
    log.info(
        "\nAll taxonomy figures complete. Suggested next steps:\n"
        "  - Run ANCOM-BC for differential abundance:\n"
        "      python 10_ancombc.py --table <feature_table.qza> \\\n"
        "        --metadata <metadata.tsv> --group-col Group \\\n"
        "        --outdir results/<marker>/ancombc/\n"
        "  - Or review figures and proceed to manuscript assembly."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
