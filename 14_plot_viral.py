#!/usr/bin/env python3
"""
14_plot_viral.py
================
Generate publication-quality viral detection figures from amplicon sequencing
data classified against NCBI nt via BLAST. Designed for the TGF-IYG
pan-herpesvirus primer dataset but generalisable to any marker where
"no-hit" reads represent a putative novel viral signal.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PIPELINE POSITION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  loon_amplicon_analysis.xlsx (herpes sheet)
              ↓
  10_plot_viral.py   ← metadata TSV (with Group, COD_broad, Season)
              ↓
  herpes_relabund_{palette}.png/.svg     relative abundance barplot
  herpes_presence_{palette}.png/.svg     presence/absence barplot

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WHAT THIS SCRIPT PRODUCES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Relative abundance barplot:
    - Per-sample stacked bar: no-hit reads vs. classified reads vs. other
    - Grouped by --group-by column (default: Group)
    - Group labels with n= counts and vertical dashed dividers
    - Relative abundance calculated as no-hit reads / total reads

  Presence/absence barplot:
    - Binary detection per sample above --threshold (default: 0.01 = 1%)
    - Filled bar = detected; empty/hatched bar = not detected
    - Group-level detection rate annotated above each group
    - Fisher's exact test p-value shown between groups (2-group only)

  Both figures:
    - PNG (300 dpi) + SVG (vector)
    - Purple palette (dark-to-light) or Wong 2011 colorblind-safe palette
    - No title (--no-title default) for journal submission
    - Samples ordered within groups by no-hit relative abundance

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SAMPLE QUALITY FILTERING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Before detection is called, samples with fewer than --min-sample-reads
total reads are dropped. This prevents near-empty samples from producing
spurious detections or non-detections purely from relative abundance math.

  --min-sample-reads  N   Drop samples with < N total reads (default: 500)

Why this matters: a sample with 10 total reads and 1 no-hit read is
called "detected" at 10% relative abundance, while a sample with 10,000
reads and 99 no-hit reads is called "not detected" at 0.99%. The relative
abundance threshold is only meaningful when applied to samples with
sufficient read depth. This filter is applied before any detection
calling, matching the approach used in 11b_presence_absence.py.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PALETTES — kept in sync with 09_plot_diversity.py and 10_plot_taxonomy.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  --palette purple   Dark-to-light purple for group bars; grey for classified.
  --palette wong     Wong 2011 8-color colorblind-safe palette.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
USAGE EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # Adenovirus — absolute read count threshold (>=10 reads)
  python scripts/10_plot_viral.py \\
      --xlsx   loon_amplicon_analysis.xlsx \\
      --sheet  adeno \\
      --metadata metadata/qiime/metadata_16S_updated.tsv \\
      --group-by Group \\
      --group-order Diseased Trauma \\
      --palette wong \\
      --no-title \\
      --min-reads 10 \\
      --output-stem adeno_Group_wong \\
      --outdir results/adenovirus/figures/

  # Herpesvirus — default relative abundance threshold, drop low-depth samples
  python scripts/10_plot_viral.py \\
      --xlsx   loon_amplicon_analysis.xlsx \\
      --sheet  herpes \\
      --metadata metadata/qiime/metadata_16S_updated.tsv \\
      --group-by Group \\
      --group-order Diseased Trauma \\
      --palette purple \\
      --no-title \\
      --min-sample-reads 500 \\
      --outdir results/herpes/figures/

  # Wong palette
  python scripts/10_plot_viral.py \\
      --xlsx   loon_amplicon_analysis.xlsx \\
      --sheet  herpes \\
      --metadata metadata/qiime/metadata_16S_updated.tsv \\
      --group-by Group \\
      --group-order Diseased Trauma \\
      --palette wong \\
      --no-title \\
      --outdir results/herpes/figures/

  # COD_broad grouping (Lead / Parasitic_Infectious / Trauma)
  python scripts/10_plot_viral.py \\
      --xlsx   loon_amplicon_analysis.xlsx \\
      --sheet  herpes \\
      --metadata metadata/qiime/metadata_16S_updated.tsv \\
      --group-by COD_broad \\
      --group-order Lead Parasitic_Infectious Trauma \\
      --palette wong \\
      --no-title \\
      --outdir results/herpes/figures/

  # Custom threshold and output stem
  python scripts/10_plot_viral.py \\
      --xlsx   loon_amplicon_analysis.xlsx \\
      --sheet  herpes \\
      --metadata metadata/qiime/metadata_16S_updated.tsv \\
      --group-by Group \\
      --threshold 0.50 \\
      --output-stem herpes_threshold50 \\
      --palette purple \\
      --outdir results/herpes/figures/

  # List available metadata columns
  python scripts/10_plot_viral.py \\
      --xlsx loon_amplicon_analysis.xlsx --sheet herpes \\
      --metadata metadata/qiime/metadata_16S_updated.tsv \\
      --list-columns

  # Dry run — preview per-sample stats and filter results without writing figures
  python scripts/10_plot_viral.py \\
      --xlsx   loon_amplicon_analysis.xlsx \\
      --sheet  herpes \\
      --metadata metadata/qiime/metadata_16S_updated.tsv \\
      --group-by Group \\
      --dry-run

Dependencies:
  pip install matplotlib numpy pandas scipy openpyxl
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
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants — kept in sync with 09_plot_diversity.py and 10_plot_taxonomy.py
# ---------------------------------------------------------------------------

FIGURE_DPI        = 300
FONT_FAMILY       = "Arial"
FONT_SIZE_AXIS    = 11
FONT_SIZE_TICK    = 9
FONT_SIZE_LABEL   = 9
FONT_SIZE_ANNOT   = 8
FONT_SIZE_TITLE   = 12

BLAST_COL  = "BLAST nt"
NO_HIT_TAG = "no-hit"
TV_REGEX   = r"(TV\d+)"

# ---------------------------------------------------------------------------
# Palettes — identical structure to 10_plot_taxonomy.py
# ---------------------------------------------------------------------------

PALETTES: Dict[str, Dict] = {
    "purple": {
        "group_colors": [
            "#7B2D8B",  # dark purple  — Diseased / Group 1
            "#C19FD8",  # lavender     — Trauma / Group 2
            "#4B1369",  # deep purple  — Group 3
            "#D09EE0",  # light purple — Group 4
            "#2D0A40", "#E0BAEC", "#9870B0", "#F0D6F5",
        ],
        "nohit_color":      "#4B0082",
        "classified_color": "#C9B8D8",
        "absent_color":     "#E8E0EE",
        "absent_hatch":     "///",
    },
    "wong": {
        "group_colors": [
            "#0072B2",  # blue    — Diseased / Group 1
            "#E69F00",  # orange  — Trauma / Group 2
            "#009E73",  # green   — Group 3
            "#CC79A7",  # pink    — Group 4
            "#56B4E9", "#D55E00", "#F0E442", "#999999",
        ],
        "nohit_color":      "#0072B2",
        "classified_color": "#BBBBBB",
        "absent_color":     "#EEEEEE",
        "absent_hatch":     "///",
    },
}

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_xlsx(xlsx_path: Path, sheet: str) -> pd.DataFrame:
    """Load the BLAST amplicon sheet from the Excel workbook."""
    log.info("Loading %s [sheet: %s]", xlsx_path.name, sheet)
    df = pd.read_excel(xlsx_path, sheet_name=sheet)
    if BLAST_COL not in df.columns:
        raise ValueError(
            f"Column '{BLAST_COL}' not found in sheet '{sheet}'. "
            f"Available: {list(df.columns[:10])}"
        )
    return df


def load_metadata(metadata_path: Path, group_column: str) -> pd.DataFrame:
    """
    Load QIIME2 metadata TSV. Skips the #q2:types row if present.
    Returns DataFrame indexed by TV ID with the group column included.
    """
    df = pd.read_csv(metadata_path, sep="\t", dtype=str)
    df = df[~df.iloc[:, 0].str.startswith("#", na=False)].copy()
    df = df.reset_index(drop=True)

    sid_col = df.columns[0]
    df["_TV"] = df[sid_col].str.extract(TV_REGEX, expand=False)
    df = df.dropna(subset=["_TV"])
    df = df.set_index("_TV")

    if group_column not in df.columns:
        raise ValueError(
            f"Column '{group_column}' not in metadata. "
            f"Available: {list(df.columns)}"
        )
    return df


def compute_per_sample(df: pd.DataFrame, taxon_filter: str = None) -> pd.DataFrame:
    """
    From the raw amplicon sheet, compute per-sample:
      - total_reads
      - nohit_reads   (signal reads -- either no-hit or taxon-filtered)
      - classified_reads
      - nohit_relabund  (signal / total)

    taxon_filter: if provided, use rows whose BLAST nt contains this string
                  as the signal (e.g. 'Aviadenovirus' for the adeno sheet).
                  If None, uses no-hit rows as signal (default, for herpesvirus).

    Returns a DataFrame indexed by TV ID.
    """
    sample_cols = [c for c in df.columns if str(c).startswith("TV")]

    if taxon_filter:
        signal_mask = df[BLAST_COL].str.contains(taxon_filter, case=False, na=False)
        log.info("Taxon filter '%s': %d / %d OTUs match", taxon_filter,
                 signal_mask.sum(), len(df))
    else:
        signal_mask = df[BLAST_COL].str.startswith(NO_HIT_TAG, na=False)

    total      = df[sample_cols].sum()
    nohit      = df.loc[signal_mask, sample_cols].sum()
    classified = total - nohit

    result = pd.DataFrame({
        "total_reads":      total,
        "nohit_reads":      nohit,
        "classified_reads": classified,
        "nohit_relabund":   nohit / total.replace(0, np.nan),
    })

    result.index = pd.Series(result.index).str.extract(TV_REGEX, expand=False).values
    result.index.name = "TV"

    log.info(
        "Loaded %d samples | total reads: %d–%d | no-hit relabund: %.1f%%–%.1f%%",
        len(result),
        int(result["total_reads"].min()), int(result["total_reads"].max()),
        result["nohit_relabund"].min() * 100,
        result["nohit_relabund"].max() * 100,
    )
    return result


def filter_low_depth_samples(
    per_sample: pd.DataFrame,
    min_sample_reads: int,
) -> pd.DataFrame:
    """
    Drop samples with fewer than min_sample_reads total reads before any
    detection calling.

    Low-depth samples produce unreliable relative abundance values — a
    sample with 10 reads and 1 no-hit read registers 10% no-hit relabund
    and is called 'detected' under the default 1% threshold, even though
    1 read is not meaningful evidence. This filter mirrors the approach in
    11b_presence_absence.py and should always be applied before threshold-
    based detection.

    Returns the filtered DataFrame and logs dropped sample IDs.
    """
    if min_sample_reads <= 0:
        return per_sample

    mask_keep   = per_sample["total_reads"] >= min_sample_reads
    mask_drop   = ~mask_keep
    n_dropped   = int(mask_drop.sum())

    if n_dropped > 0:
        dropped_ids = per_sample.index[mask_drop].tolist()
        log.warning(
            "Dropping %d sample(s) with < %d total reads: %s",
            n_dropped, min_sample_reads,
            dropped_ids[:10] + (["..."] if len(dropped_ids) > 10 else []),
        )
    else:
        log.info(
            "Sample depth filter (>= %d reads): all %d samples retained",
            min_sample_reads, len(per_sample),
        )

    filtered = per_sample[mask_keep].copy()
    log.info(
        "Samples after depth filter: %d / %d",
        len(filtered), len(per_sample),
    )
    return filtered


def merge_with_metadata(
    per_sample: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Join per-sample stats with metadata group assignments."""
    merged = per_sample.join(metadata[[group_column]], how="inner")
    n_lost = len(per_sample) - len(merged)
    if n_lost:
        log.warning("%d sample(s) in xlsx not matched to metadata (excluded)", n_lost)
    log.info("Matched %d / %d samples to metadata", len(merged), len(per_sample))
    merged = merged.rename(columns={group_column: "_group"})
    merged = merged[merged["_group"].notna() & (merged["_group"] != "missing")]
    return merged


# ---------------------------------------------------------------------------
# Relative abundance barplot
# ---------------------------------------------------------------------------

def plot_relabund(
    data: pd.DataFrame,
    group_order: List[str],
    palette_name: str,
    title: Optional[str],
    outpath_stem: Path,
    nohit_label: str = "No-hit reads",
) -> None:
    """Stacked per-sample barplot: no-hit vs classified reads (relative %)."""

    palette = PALETTES[palette_name]
    groups  = [g for g in group_order if g in data["_group"].values]

    ordered_data = []
    for g in groups:
        grp = data[data["_group"] == g].sort_values("nohit_relabund", ascending=False)
        ordered_data.append(grp)
    plot_df = pd.concat(ordered_data)

    fig, ax = plt.subplots(figsize=(max(10, len(plot_df) * 0.38), 5))
    plt.rcParams["font.family"] = FONT_FAMILY

    x     = np.arange(len(plot_df))
    bar_w = 0.85

    nohit_pct      = plot_df["nohit_relabund"].fillna(0).values * 100
    classified_pct = 100 - nohit_pct

    ax.bar(x, classified_pct, width=bar_w,
           color=palette["classified_color"], label="Classified reads", linewidth=0)
    ax.bar(x, nohit_pct, width=bar_w, bottom=classified_pct,
           color=palette["nohit_color"], label=nohit_label, linewidth=0)

    ax.set_xlim(-0.5, len(plot_df) - 0.5)
    ax.set_ylim(0, 115)
    ax.set_ylabel("Relative Abundance (%)\namong total reads", fontsize=FONT_SIZE_AXIS)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [tv.replace("-lung-TGF-IYG", "").replace("-lung-TGF-IYG-NX", "")
         for tv in plot_df.index],
        rotation=90, fontsize=FONT_SIZE_TICK - 1, ha="center",
    )
    ax.tick_params(axis="y", labelsize=FONT_SIZE_TICK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)

    pos = 0
    for gi, g in enumerate(groups):
        n   = (plot_df["_group"] == g).sum()
        mid = pos + (n - 1) / 2
        gc  = palette["group_colors"][gi % len(palette["group_colors"])]
        ax.text(mid, 108, f"{g}  (n={n})", ha="center", va="bottom",
                fontsize=FONT_SIZE_LABEL, fontweight="bold", color=gc)
        if gi < len(groups) - 1:
            ax.axvline(pos + n - 0.5, color="#999999", linewidth=0.8,
                       linestyle="--", alpha=0.7)
        pos += n

    if title:
        ax.set_title(title, fontsize=FONT_SIZE_TITLE, pad=8)

    ax.legend(loc="lower left", fontsize=FONT_SIZE_ANNOT, framealpha=0.9,
              edgecolor="#CCCCCC")

    fig.tight_layout()
    _save(fig, outpath_stem)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Presence / absence barplot
# ---------------------------------------------------------------------------

def plot_presence_absence(
    data: pd.DataFrame,
    group_order: List[str],
    palette_name: str,
    threshold: float,
    title: Optional[str],
    outpath_stem: Path,
    min_reads: Optional[int] = None,
) -> None:
    """
    Binary detection barplot. Bars reach 1.0 if detected above threshold,
    show hatched empty bar if not detected. Group-level detection rate
    annotated above each group. Fisher's exact test for 2-group comparisons.

    If min_reads is set, detection is based on absolute read count >= min_reads
    instead of relative abundance >= threshold.
    """
    palette = PALETTES[palette_name]
    groups  = [g for g in group_order if g in data["_group"].values]

    data = data.copy()
    if min_reads is not None:
        data["detected"] = data["nohit_reads"] >= min_reads
        thresh_label = f">={min_reads} reads"
    else:
        data["detected"] = data["nohit_relabund"] >= threshold
        thresh_label = f">={threshold*100:.0f}%"

    ordered_data = []
    for g in groups:
        grp = data[data["_group"] == g].copy()
        grp = grp.sort_values(["detected", "nohit_relabund"], ascending=[False, False])
        ordered_data.append(grp)
    plot_df = pd.concat(ordered_data)

    fig, ax = plt.subplots(figsize=(max(10, len(plot_df) * 0.38), 4))
    plt.rcParams["font.family"] = FONT_FAMILY

    x     = np.arange(len(plot_df))
    bar_w = 0.85

    for gi, g in enumerate(groups):
        mask = plot_df["_group"] == g
        xi   = x[mask.values]
        det  = plot_df.loc[mask, "detected"].values
        gc   = palette["group_colors"][gi % len(palette["group_colors"])]

        for xi_val, d in zip(xi, det):
            if d:
                ax.bar(xi_val, 1.0, width=bar_w, color=gc,
                       linewidth=0.5, edgecolor="#555555")
            else:
                ax.bar(xi_val, 1.0, width=bar_w,
                       color=palette["absent_color"],
                       hatch=palette["absent_hatch"],
                       linewidth=0.5, edgecolor="#AAAAAA")

    ax.axhline(threshold, color="#CC0000", linewidth=0.8, linestyle=":",
               label=f"Detection threshold ({thresh_label})")

    ax.set_xlim(-0.5, len(plot_df) - 0.5)
    ax.set_ylim(0, 1.35)
    ax.set_yticks([0, 0.25, 0.50, 0.75, 1.0])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=FONT_SIZE_TICK)
    ax.set_ylabel("No-hit relative abundance", fontsize=FONT_SIZE_AXIS)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [tv.replace("-lung-TGF-IYG", "").replace("-lung-TGF-IYG-NX", "")
         for tv in plot_df.index],
        rotation=90, fontsize=FONT_SIZE_TICK - 1, ha="center",
    )
    ax.tick_params(axis="y", labelsize=FONT_SIZE_TICK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)

    group_det_counts: List[int] = []
    group_ns:         List[int] = []
    pos = 0
    for gi, g in enumerate(groups):
        n     = (plot_df["_group"] == g).sum()
        n_det = int(plot_df.loc[plot_df["_group"] == g, "detected"].sum())
        mid   = pos + (n - 1) / 2
        gc    = palette["group_colors"][gi % len(palette["group_colors"])]
        ax.text(mid, 1.27, f"{g}  (n={n})", ha="center", va="bottom",
                fontsize=FONT_SIZE_LABEL, fontweight="bold", color=gc)
        ax.text(mid, 1.12, f"{n_det}/{n} detected",
                ha="center", va="bottom", fontsize=FONT_SIZE_ANNOT, color=gc)
        if gi < len(groups) - 1:
            ax.axvline(pos + n - 0.5, color="#999999", linewidth=0.8,
                       linestyle="--", alpha=0.7)
        group_det_counts.append(n_det)
        group_ns.append(n)
        pos += n

    if len(groups) == 2:
        a, b   = group_det_counts
        na, nb = group_ns
        contingency = [[a, na - a], [b, nb - b]]
        _, pval = fisher_exact(contingency)
        pstr = f"Fisher's exact  p = {pval:.3f}"
        if pval < 0.05:
            pstr += " *"
        ax.text(0.98, 0.03, pstr, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=FONT_SIZE_ANNOT,
                style="italic", color="#444444")
        log.info("Fisher's exact: %s", pstr)

    present_patch = mpatches.Patch(color=palette["group_colors"][0],
                                   label=f"Detected ({thresh_label})")
    absent_patch  = mpatches.Patch(facecolor=palette["absent_color"],
                                   hatch=palette["absent_hatch"],
                                   edgecolor="#AAAAAA",
                                   label="Not detected (below threshold)")
    ax.legend(handles=[present_patch, absent_patch],
              loc="lower left", fontsize=FONT_SIZE_ANNOT,
              framealpha=0.9, edgecolor="#CCCCCC")

    if title:
        ax.set_title(title, fontsize=FONT_SIZE_TITLE, pad=8)

    fig.tight_layout()
    _save(fig, outpath_stem)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Save helper
# ---------------------------------------------------------------------------

def _save(fig: plt.Figure, stem: Path) -> None:
    for ext in (".png", ".svg"):
        fpath = stem.with_suffix(ext)
        fpath.parent.mkdir(parents=True, exist_ok=True)
        kw = {"dpi": FIGURE_DPI} if ext == ".png" else {}
        fig.savefig(fpath, bbox_inches="tight", **kw)
        log.info("  Saved: %s", fpath)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """
    Build and return the argument parser for 10_plot_viral.py.

    Key arguments: --xlsx (BLAST amplicon Excel workbook), --metadata
    (QIIME2 TSV), --group-by, --palette, --threshold, --min-sample-reads,
    --outdir, --output-stem. A --list-columns flag inspects available
    metadata columns without producing any figures.
    """
    p = argparse.ArgumentParser(
        prog="10_plot_viral.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    req = p.add_argument_group("required arguments")
    req.add_argument("--xlsx",     required=True,
                     help="Path to Excel workbook (e.g. loon_amplicon_analysis.xlsx).")
    req.add_argument("--sheet",    required=True,
                     help="Sheet name containing BLAST nt amplicon data (e.g. herpes).")
    req.add_argument("--metadata", required=True,
                     help="QIIME2 metadata TSV with Group, COD_broad, etc.")

    plot = p.add_argument_group("plot options")
    plot.add_argument("--group-by",    default="Group",
                      help="Metadata column to group samples by. Default: Group.")
    plot.add_argument("--group-order", nargs="+", default=None,
                      help="Explicit group order left-to-right. Default: alphabetical.")
    plot.add_argument("--threshold",   type=float, default=0.01,
                      help="No-hit relative abundance threshold for presence/absence (0–1). "
                           "Default: 0.01 (1%%). Ignored if --min-reads is set.")
    plot.add_argument("--taxon-filter", default=None, metavar="TAXON",
                      help="Filter signal to OTUs whose BLAST nt contains TAXON "
                           "(e.g. 'Aviadenovirus'). Use for sheets where the signal "
                           "is a classified taxon rather than no-hit reads. "
                           "Default: None (uses no-hit reads as signal, correct for herpesvirus).")
    plot.add_argument("--min-reads",   type=int, default=None,
                      help="Absolute no-hit read count threshold for presence/absence. "
                           "If set, overrides --threshold. Use for adenovirus (e.g. --min-reads 10).")
    plot.add_argument("--palette",
                      choices=list(PALETTES.keys()), default="purple",
                      help="Color palette. Default: purple.")
    plot.add_argument("--title",       default=None,
                      help="Override figure title.")
    plot.add_argument("--no-title",    action="store_true", default=False,
                      help="Suppress figure title (use for journal submission).")
    plot.add_argument("--output-stem", default=None,
                      help="Base filename stem (no extension). "
                           "Default: {sheet}_{group_by}_{palette}.")
    plot.add_argument("--outdir",      default=".",
                      help="Output directory. Default: current directory.")

    filt = p.add_argument_group("sample quality filter")
    filt.add_argument(
        "--min-sample-reads", type=int, default=500, metavar="N",
        help=(
            "Drop samples with fewer than N total reads before detection "
            "calling. Prevents near-empty samples from producing spurious "
            "relative-abundance detections. Default: 500. "
            "Set to 0 to disable."
        ),
    )

    util = p.add_argument_group("utility")
    util.add_argument("--list-columns", action="store_true",
                      help="Print available metadata columns and exit.")
    util.add_argument("--dry-run", action="store_true",
                      help="Compute per-sample stats and print summary; do not save figures.")

    return p


def main(argv: Optional[List[str]] = None) -> int:
    """
    Parse arguments and generate viral detection barplots.

    Loads the BLAST amplicon sheet from the Excel workbook, computes
    per-sample no-hit read counts, drops low-depth samples
    (--min-sample-reads), joins with metadata group assignments, and
    produces two figure types (relative abundance and presence/absence).
    Returns 0 on success, 2 if required input files are missing, 1 on error.
    """
    parser = build_parser()
    args   = parser.parse_args(argv)

    xlsx_path = Path(args.xlsx)
    meta_path = Path(args.metadata)
    outdir    = Path(args.outdir)

    if not xlsx_path.exists():
        log.error("Excel file not found: %s", xlsx_path)
        return 1
    if not meta_path.exists():
        log.error("Metadata file not found: %s", meta_path)
        return 1

    # --list-columns
    if args.list_columns:
        meta = pd.read_csv(meta_path, sep="\t", dtype=str, nrows=3)
        meta = meta[~meta.iloc[:, 0].str.startswith("#", na=False)]
        print("Available metadata columns:")
        for c in meta.columns:
            print(f"  {c}")
        return 0

    # Load + compute
    df         = load_xlsx(xlsx_path, args.sheet)
    per_sample = compute_per_sample(df, taxon_filter=args.taxon_filter)

    # ── Sample quality filter (NEW) ───────────────────────────────────────
    if args.min_sample_reads > 0:
        per_sample = filter_low_depth_samples(per_sample, args.min_sample_reads)
    else:
        log.info("Sample depth filter disabled (--min-sample-reads 0)")

    if per_sample.empty:
        log.error(
            "No samples remain after depth filter "
            "(--min-sample-reads %d). Lower the threshold or check your data.",
            args.min_sample_reads,
        )
        return 1

    metadata = load_metadata(meta_path, args.group_by)
    data     = merge_with_metadata(per_sample, metadata, args.group_by)

    if data.empty:
        log.error(
            "No samples remain after merging with metadata. "
            "Check that sample IDs in the xlsx match the metadata."
        )
        return 1

    # Group order
    if args.group_order:
        group_order = args.group_order
    else:
        group_order = sorted(data["_group"].dropna().unique().tolist())

    missing_groups = [g for g in group_order if g not in data["_group"].values]
    if missing_groups:
        log.warning("Groups in --group-order not found in data: %s", missing_groups)

    log.info("Groups (after depth filter): %s", group_order)
    for g in group_order:
        n = (data["_group"] == g).sum()
        if args.min_reads is not None:
            n_det = int((data.loc[data["_group"] == g, "nohit_reads"] >= args.min_reads).sum())
            log.info("  %s: n=%d  detected=%d/%d (>=%d reads)",
                     g, n, n_det, n, args.min_reads)
        else:
            n_det = int((data.loc[data["_group"] == g, "nohit_relabund"] >= args.threshold).sum())
            log.info("  %s: n=%d  detected=%d/%d (>=%.0f%%)",
                     g, n, n_det, n, args.threshold * 100)

    if args.dry_run:
        log.info("DRY RUN — no figures written.")
        return 0

    # Output stem
    stem_base = args.output_stem or f"{args.sheet}_{args.group_by}_{args.palette}"
    outdir.mkdir(parents=True, exist_ok=True)

    title = None if args.no_title else (args.title or None)

    sheet_labels = {"herpes": "herpesvirus", "adeno": "Aviadenovirus reads"}
    nohit_label  = f"No-hit reads ({sheet_labels.get(args.sheet, args.sheet)})"

    log.info("Palette:   %s", args.palette)
    if args.min_reads is not None:
        log.info("Threshold (presence/absence): >=%d reads (absolute)", args.min_reads)
    else:
        log.info("Threshold (presence/absence): %.0f%%", args.threshold * 100)
    if args.min_sample_reads > 0:
        log.info("Sample depth filter applied : >=%d total reads", args.min_sample_reads)

    relabund_stem = outdir / f"{stem_base}_relabund"
    plot_relabund(
        data         = data,
        group_order  = group_order,
        palette_name = args.palette,
        title        = title,
        outpath_stem = relabund_stem,
        nohit_label  = nohit_label,
    )

    pa_stem = outdir / f"{stem_base}_presence"
    plot_presence_absence(
        data         = data,
        group_order  = group_order,
        palette_name = args.palette,
        threshold    = args.threshold,
        title        = title,
        outpath_stem = pa_stem,
        min_reads    = args.min_reads,
    )

    log.info("=== Done. Figures in: %s ===", outdir)
    log.info("")
    log.info("Suggested next steps:")
    log.info("  - Run Fisher's exact in R or scipy for multi-group comparisons")
    log.info("  - Extract no-hit OTU sequences for BLASTn confirmation:")
    log.info("      grep no-hit sequences from xlsx -> FASTA -> blastn -db nt")
    log.info("  - Consider relative abundance >50%% threshold for stringent detection")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
