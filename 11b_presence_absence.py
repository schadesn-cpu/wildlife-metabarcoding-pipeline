#!/usr/bin/env python3
"""
11b_presence_absence.py
=======================
Convert a taxonomy count table to presence/absence, compute detection
frequencies, and generate detection barplots. Designed as a universal
presence/absence framework for any amplicon metabarcoding marker where
relative read abundance is not a reliable proxy for biological quantity.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHILOSOPHY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Amplicon metabarcoding count data is not a standard community diversity
dataset in many study designs. Each sample is treated as an independent
detection unit. The key questions are:

  - Was taxon X detected in sample Y?          (presence/absence)
  - Across all samples, how frequently was X detected?  (detection freq)

This framework is appropriate whenever relative read abundance cannot be
defended as a proxy for biological quantity. Common cases include:

  - Dietary metabarcoding (MiFish, cytb, COI): PCR efficiency differs
    across prey taxa, mitochondrial copy number varies by tissue and
    species, and digestion state affects DNA yield independently of
    consumption amount. Deagle et al. (2019, Mol. Ecol.) is the
    standard reference for why presence/absence is the preferred unit.

  - Blood meal metabarcoding (COI, 12S, 16S): blood meal DNA is
    degraded and variably amplified; a single sample may contain DNA
    from multiple feeding events at different times. Framework
    recommended by Borland & Kading (2021, DOI: 10.3390/insects12010037)
    and Balasubramanian et al. (2024, DOI: 10.1002/edn3.522).

  - Any marker with low rarefaction depth: at <=200 reads per sample,
    relative abundance estimates have variance too large to be
    meaningfully interpreted. Presence/absence is the only statistically
    honest choice.

The approach: apply read-count quality filters, convert to binary 0/1,
then summarize as detection frequency per group.

Use --sample-label to replace the generic "sample" with a study-specific
term in all output text (e.g. "loon", "tick", "individual").

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PIPELINE POSITION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  07_taxonomy_table.py  ->  taxonomy_counts_L{N}_{marker}.tsv
                                        |
                        11b_presence_absence.py  <- metadata TSV (optional)
                                        |
              presence_absence_L{N}_{marker}.tsv   <- input to 08_run_diversity_stats.py
              detection_freq_{marker}.tsv
              detection_freq_by_{group}_{marker}.tsv
              detection_summary_{marker}.txt
              detection_barplot_{marker}_{group}.png/.svg

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FILTERING LOGIC
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Applied in order before presence/absence conversion:

  1. --min-sample-reads   Drop samples with fewer total reads than this.
                          Avoids calling detections in near-empty samples.
                          Default: 500. Set lower for low-depth markers
                          (e.g. --min-sample-reads 50 for cytb at depth 200).

  2. --min-taxon-reads    Drop taxa with fewer reads across the whole
                          dataset than this. Removes globally rare taxa
                          that are likely amplification artifacts.
                          Default: 50

  3. --min-relabund       Within each sample, zero out any taxon whose
                          relative read abundance is below this threshold
                          before converting to 0/1. Optional; useful for
                          removing low-level cross-contamination.
                          Default: 0.0 (disabled)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
OUTPUTS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  presence_absence_L{N}_{marker}.tsv
      Binary 0/1 table: taxa (rows) x samples (columns).
      Direct input to 08_run_diversity_stats.py for Jaccard + PERMANOVA.

  detection_freq_{marker}.tsv
      Per-taxon detection frequency across all retained samples.
      Columns: taxon, n_detected, n_samples, detection_freq, pct_detected

  detection_freq_by_{group}_{marker}.tsv
      Same as above but computed per group (requires --metadata --group-by).
      Wide format: taxon x group, values = detection frequency (0-1).

  detection_summary_{marker}.txt
      Human-readable summary: how many samples had at least one detection,
      per-taxon detection frequencies, and per-group breakdowns.

  detection_barplot_{marker}.png / .svg
      Horizontal bar chart of detection frequency (%) per taxon.
      Sorted by overall detection frequency descending.
      If --group-by is set, grouped bars are drawn per group.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
USAGE EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # MiFish diet — group comparison, standard thresholds
  python 11b_presence_absence.py \\
      --counts            results/MiFish/all/taxonomy/taxonomy_counts_L7_MiFish.tsv \\
      --metadata          metadata/qiime/metadata_MiFish.tsv \\
      --marker            MiFish \\
      --group-by          Group \\
      --min-sample-reads  10000 \\
      --min-taxon-reads   10 \\
      --min-relabund      0.01 \\
      --sample-label      loon \\
      --outdir            results/MiFish/all/presence_absence/

  # cytb — low depth marker, thresholds adjusted to match reality
  # At rarefaction depth 200, --min-sample-reads must be <= 200.
  # --min-relabund 0.01 requires >= 2 reads for a detection call at depth 200,
  # which is the practical minimum for distinguishing signal from noise.
  python 11b_presence_absence.py \\
      --counts            results/cytb/all/taxonomy/taxonomy_counts_L7_cytb.tsv \\
      --metadata          metadata/qiime/metadata_cytb.tsv \\
      --marker            cytb \\
      --group-by          Group \\
      --min-sample-reads  50 \\
      --min-taxon-reads   5 \\
      --min-relabund      0.01 \\
      --sample-label      loon \\
      --outdir            results/cytb/all/presence_absence/

  # Blood meal tick study — original use case still fully supported
  python 11b_presence_absence.py \\
      --counts            results/COI/all/taxonomy/taxonomy_counts_L7_COI.tsv \\
      --metadata          metadata/qiime/metadata_COI.tsv \\
      --marker            COI \\
      --group-by          Site \\
      --min-sample-reads  1000 \\
      --min-taxon-reads   100 \\
      --min-relabund      0.01 \\
      --sample-label      tick \\
      --outdir            results/COI/all/presence_absence/

  # Dry run -- see what would happen without writing files
  python 11b_presence_absence.py \\
      --counts   results/cytb/all/taxonomy/taxonomy_counts_L7_cytb.tsv \\
      --marker   cytb \\
      --outdir   results/cytb/all/presence_absence/ \\
      --dry-run

Dependencies:
  pip install matplotlib numpy pandas
"""
from __future__ import annotations

import argparse
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

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Style constants — kept in sync with 09_plot_diversity.py / 10_plot_taxonomy.py
# ---------------------------------------------------------------------------

FIGURE_DPI       = 300
FONT_SIZE_TITLE  = 13
FONT_SIZE_AXIS   = 11
FONT_SIZE_TICK   = 9
FONT_SIZE_LEGEND = 8.5

# Group colors: match 10_plot_taxonomy.py purple palette exactly
GROUP_COLORS = [
    "#7B2D8B",  # dark purple   — Group 1
    "#C19FD8",  # lavender      — Group 2
    "#4B1369",  # deep purple   — Group 3
    "#D09EE0",  # light purple  — Group 4
    "#2D0A40",  # very dark     — Group 5
    "#E0BAEC",  # pale lavender — Group 6
    "#9870B0",  # mid purple    — Group 7
    "#F0D6F5",  # near-white    — Group 8
]

# Single-bar color for overall (no grouping) plots
OVERALL_BAR_COLOR = "#7B2D8B"


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def eprint(*args, **kwargs) -> None:
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def safe_mkdir(path: Path) -> None:
    """Create path and any missing parents; no-op if it already exists."""
    path.mkdir(parents=True, exist_ok=True)


def save_tsv(df: pd.DataFrame, path: Path, dry_run: bool = False) -> None:
    """
    Write a DataFrame to a tab-separated file, creating parent dirs as needed.

    Logs the output path and dimensions. In dry-run mode the file is not
    written but the log message is still printed so intended output can
    be audited before committing to a full run.
    """
    log.info("  Writing: %s  (%d rows × %d cols)", path, len(df), len(df.columns))
    if dry_run:
        return
    safe_mkdir(path.parent)
    df.to_csv(path, sep="\t", index=True)


def save_text(text: str, path: Path, dry_run: bool = False) -> None:
    """Write a plain-text string to path. No-op in dry-run mode."""
    log.info("  Writing: %s", path)
    if dry_run:
        return
    safe_mkdir(path.parent)
    path.write_text(text, encoding="utf-8")


def save_figure(fig: plt.Figure, path: Path, dry_run: bool = False) -> None:
    """Save a matplotlib figure as both PNG and SVG. No-op in dry-run mode."""
    png_path = path.with_suffix(".png")
    svg_path = path.with_suffix(".svg")
    log.info("  Writing: %s", png_path)
    log.info("  Writing: %s", svg_path)
    if dry_run:
        plt.close(fig)
        return
    safe_mkdir(path.parent)
    fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Metadata loading
# ---------------------------------------------------------------------------

def load_metadata(metadata_path: Path, group_by: str) -> pd.Series:
    """
    Load a QIIME2-format metadata TSV and return a Series mapping
    sample-id → group label for the requested column.

    Handles the optional #q2:types second row automatically. Raises
    ValueError if the requested column is not found.
    """
    df = pd.read_csv(metadata_path, sep="\t", dtype=str, index_col=0)
    df.index.name = "sample-id"

    # Drop the optional QIIME2 directive row (#q2:types)
    if df.index[0].startswith("#"):
        df = df.iloc[1:]

    if group_by not in df.columns:
        available = ", ".join(df.columns.tolist())
        raise ValueError(
            f"Column '{group_by}' not found in metadata.\n"
            f"Available columns: {available}"
        )

    log.info(
        "  Metadata: %d samples, using column '%s' (%d unique groups)",
        len(df),
        group_by,
        df[group_by].nunique(),
    )
    return df[group_by].rename("group")


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def apply_filters(
    counts: pd.DataFrame,
    min_sample_reads: int,
    min_taxon_reads: int,
    min_relabund: float,
) -> pd.DataFrame:
    """
    Apply quality filters to a count table (taxa × samples) before
    presence/absence conversion.

    Steps (in order):
      1. Drop samples with fewer than min_sample_reads total reads.
      2. Drop taxa with fewer than min_taxon_reads reads across the dataset.
      3. If min_relabund > 0, zero out per-sample entries whose within-sample
         relative abundance is below the threshold.

    Returns the filtered count DataFrame.
    """
    n_samples_start = len(counts.columns)
    n_taxa_start    = len(counts)

    # ── 1. Sample-level read filter ───────────────────────────────────────
    sample_totals = counts.sum(axis=0)
    keep_samples  = sample_totals[sample_totals >= min_sample_reads].index
    dropped_samples = sample_totals[sample_totals < min_sample_reads].index.tolist()

    if dropped_samples:
        log.warning(
            "  Dropping %d sample(s) with <%d total reads: %s",
            len(dropped_samples), min_sample_reads, dropped_samples,
        )
    counts = counts[keep_samples]
    log.info(
        "  Sample filter (≥%d reads): %d → %d samples",
        min_sample_reads, n_samples_start, len(counts.columns),
    )

    # ── 2. Taxon-level read filter ────────────────────────────────────────
    taxon_totals = counts.sum(axis=1)
    keep_taxa    = taxon_totals[taxon_totals >= min_taxon_reads].index
    dropped_taxa = taxon_totals[taxon_totals < min_taxon_reads].index.tolist()

    if dropped_taxa:
        log.info(
            "  Dropping %d taxon/taxa with <%d total reads across dataset: %s",
            len(dropped_taxa), min_taxon_reads, dropped_taxa,
        )
    counts = counts.loc[keep_taxa]
    log.info(
        "  Taxon filter  (≥%d reads): %d → %d taxa",
        min_taxon_reads, n_taxa_start, len(counts),
    )

    # ── 3. Within-sample relative abundance filter ────────────────────────
    if min_relabund > 0.0:
        col_sums = counts.sum(axis=0).replace(0, 1)
        relabund = counts.div(col_sums, axis=1)
        mask     = relabund < min_relabund
        n_zeroed = int(mask.values.sum())
        counts   = counts.where(~mask, other=0)
        log.info(
            "  Relative abundance filter (<%s%%): zeroed %d entries",
            f"{min_relabund * 100:.1f}", n_zeroed,
        )

    if counts.empty:
        raise ValueError(
            "No data remains after filtering. "
            "Try lowering --min-sample-reads, --min-taxon-reads, "
            "or --min-relabund."
        )

    return counts


# ---------------------------------------------------------------------------
# Presence / absence conversion
# ---------------------------------------------------------------------------

def to_presence_absence(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Convert a count table (taxa × samples) to binary presence/absence (0/1).

    Values > 0 become 1; all other values become 0.
    Returns a DataFrame of the same shape with integer dtype.
    """
    pa = (counts > 0).astype(int)
    log.info(
        "  Presence/absence table: %d taxa × %d samples  "
        "(%d detections total)",
        len(pa), len(pa.columns), int(pa.values.sum()),
    )
    return pa


# ---------------------------------------------------------------------------
# Detection frequency
# ---------------------------------------------------------------------------

def compute_detection_freq(
    pa: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute per-taxon detection frequency across all samples.

    Returns a DataFrame with columns:
      taxon, n_detected, n_samples, detection_freq, pct_detected
    Sorted by detection_freq descending.
    """
    n_samples   = len(pa.columns)
    n_detected  = pa.sum(axis=1).rename("n_detected")
    freq        = (n_detected / n_samples).rename("detection_freq")
    pct         = (freq * 100).round(1).rename("pct_detected")

    df = pd.concat([n_detected, freq, pct], axis=1)
    df.insert(1, "n_samples", n_samples)
    df = df.sort_values("detection_freq", ascending=False)
    df.index.name = "taxon"
    return df


def compute_detection_freq_by_group(
    pa: pd.DataFrame,
    groups: pd.Series,
) -> pd.DataFrame:
    """
    Compute per-taxon detection frequency separately for each group.

    Returns a wide DataFrame: taxa (rows) × groups (columns),
    values = detection frequency (0–1). Sorted by overall freq descending.
    """
    group_labels = groups.reindex(pa.columns).dropna()
    unique_groups = sorted(group_labels.unique())

    result: Dict[str, pd.Series] = {}
    for grp in unique_groups:
        samples_in_grp = group_labels[group_labels == grp].index
        samples_in_grp = samples_in_grp.intersection(pa.columns)
        if len(samples_in_grp) == 0:
            log.warning("  Group '%s': no samples found in presence/absence table", grp)
            continue
        sub = pa[samples_in_grp]
        freq = sub.sum(axis=1) / len(samples_in_grp)
        result[f"{grp} (n={len(samples_in_grp)})"] = freq

    if not result:
        raise ValueError("No group-level detection frequencies could be computed.")

    df = pd.DataFrame(result)
    df.index.name = "taxon"

    # Sort by overall detection frequency (mean across groups)
    df["_sort"] = df.mean(axis=1)
    df = df.sort_values("_sort", ascending=False).drop(columns="_sort")
    return df


# ---------------------------------------------------------------------------
# Summary text (Lawrence Gordon's reporting format)
# ---------------------------------------------------------------------------

def build_summary_text(
    pa: pd.DataFrame,
    freq_df: pd.DataFrame,
    marker: str,
    groups: Optional[pd.Series] = None,
    group_by: Optional[str] = None,
    sample_label: str = "sample",
) -> str:
    """
    Build a human-readable detection summary.

    sample_label controls the unit word used throughout (e.g. "loon",
    "tick", "individual"). Defaults to "sample" for generic use.

    Reports: total samples analyzed, how many had at least one detection,
    per-taxon detection frequencies, and per-group breakdowns when groups
    are provided.
    """
    # Capitalised version for sentence starts
    sample_label_cap = sample_label.capitalize()

    n_total     = len(pa.columns)
    has_det     = (pa.sum(axis=0) > 0)
    n_with_det  = int(has_det.sum())
    pct_det     = n_with_det / n_total * 100 if n_total > 0 else 0.0

    # Detection frequencies among samples WITH at least one detection
    pa_positive = pa[has_det[has_det].index]
    n_pos = len(pa_positive.columns)

    freq_among_positive = (
        pa_positive.sum(axis=1) / n_pos * 100
        if n_pos > 0 else pa_positive.sum(axis=1) * 0
    ).sort_values(ascending=False)

    lines: List[str] = []
    lines.append(f"=== Detection Summary: {marker} ===\n")
    lines.append(f"Total {sample_label}s analyzed        : {n_total}")
    lines.append(
        f"{sample_label_cap}s with detection       : {n_with_det} / {n_total} "
        f"({pct_det:.1f}%)"
    )
    lines.append(
        f"{sample_label_cap}s with no detection    : {n_total - n_with_det} / {n_total} "
        f"(includes failed amplification and true negatives — "
        f"these are indistinguishable without additional data)\n"
    )

    lines.append(
        f"Taxon detection among {sample_label}s with at least one detection "
        f"(n={n_pos}):"
    )
    for taxon, pct in freq_among_positive.items():
        n_det = int(pa_positive.loc[taxon].sum())
        lines.append(f"  {taxon:<45s}  {n_det:>4d} / {n_pos}  ({pct:.1f}%)")

    # ── Per-group breakdown ───────────────────────────────────────────────
    if groups is not None and group_by is not None:
        lines.append(f"\nPer-group breakdown (column: '{group_by}'):")
        group_labels = groups.reindex(pa.columns).dropna()
        for grp in sorted(group_labels.unique()):
            samps   = group_labels[group_labels == grp].index
            samps   = samps.intersection(pa.columns)
            sub     = pa[samps]
            n_g     = len(samps)
            n_g_pos = int((sub.sum(axis=0) > 0).sum())
            pct_g   = n_g_pos / n_g * 100 if n_g > 0 else 0.0
            lines.append(
                f"\n  {grp} (n={n_g}): "
                f"{n_g_pos} {sample_label}s with detection ({pct_g:.1f}%)"
            )
            sub_pos = sub.loc[:, sub.sum(axis=0) > 0]
            if len(sub_pos.columns) > 0:
                gfreq = (
                    sub_pos.sum(axis=1) / len(sub_pos.columns) * 100
                ).sort_values(ascending=False)
                for taxon, pct in gfreq.items():
                    if pct > 0:
                        n_det = int(sub_pos.loc[taxon].sum())
                        lines.append(
                            f"    {taxon:<43s}  "
                            f"{n_det:>4d} / {len(sub_pos.columns)}  ({pct:.1f}%)"
                        )

    lines.append(
        f"\nMethods note: \"Presence/absence of taxa was determined from "
        f"amplicon metabarcoding data. {sample_label_cap}s with fewer than "
        f"[min_sample_reads] reads and taxa with fewer than [min_taxon_reads] "
        f"reads across the dataset were excluded prior to conversion. "
        f"Detection frequency was calculated as the proportion of {sample_label}s "
        f"in which each taxon was identified.\""
    )

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_detection_freq(
    freq_df: pd.DataFrame,
    marker: str,
    outpath: Path,
    group_freq_df: Optional[pd.DataFrame] = None,
    group_by: Optional[str] = None,
    top_n: int = 20,
    sample_label: str = "sample",
    dry_run: bool = False,
) -> None:
    """
    Generate a horizontal bar chart of detection frequency (%) per taxon.

    If group_freq_df is provided, draws grouped bars (one per group) using
    the pipeline's purple palette. Otherwise draws a single-color overall plot.

    Taxa are sorted by overall detection frequency, descending.
    Only the top_n taxa by detection frequency are shown.

    sample_label controls the unit word used in the figure subtitle
    (e.g. "loon", "tick", "sample").
    """
    # ── Prepare data ──────────────────────────────────────────────────────
    plot_data = freq_df.head(top_n).copy()
    pct_vals  = (plot_data["detection_freq"] * 100).values
    taxa      = plot_data.index.tolist()
    n_taxa    = len(taxa)

    if n_taxa == 0:
        log.warning("  No taxa to plot — skipping detection barplot")
        return

    # ── Figure geometry ───────────────────────────────────────────────────
    fig_height = max(4, n_taxa * 0.45 + 1.5)
    fig, ax    = plt.subplots(figsize=(9, fig_height))

    if group_freq_df is not None and group_by is not None:
        # ── Grouped bars ──────────────────────────────────────────────────
        groups     = group_freq_df.columns.tolist()
        n_groups   = len(groups)
        bar_height = 0.7 / n_groups
        offsets    = np.linspace(
            -(n_groups - 1) * bar_height / 2,
             (n_groups - 1) * bar_height / 2,
            n_groups,
        )
        taxa_order = [t for t in group_freq_df.index if t in taxa]

        for gi, (grp, offset) in enumerate(zip(groups, offsets)):
            color  = GROUP_COLORS[gi % len(GROUP_COLORS)]
            vals   = [
                group_freq_df.loc[t, grp] * 100 if t in group_freq_df.index else 0.0
                for t in taxa_order
            ]
            y_pos  = np.arange(len(taxa_order)) + offset
            ax.barh(
                y_pos, vals,
                height=bar_height * 0.9,
                color=color, alpha=0.88,
                label=grp,
            )

        ax.set_yticks(np.arange(len(taxa_order)))
        ax.set_yticklabels(
            [f"$\\it{{{t}}}$" if _is_binomial(t) else t for t in taxa_order],
            fontsize=FONT_SIZE_TICK,
        )
        legend = ax.legend(
            title=group_by,
            fontsize=FONT_SIZE_LEGEND,
            title_fontsize=FONT_SIZE_LEGEND,
            loc="lower right",
            framealpha=0.85,
        )
        title_suffix = f"by {group_by}"

    else:
        # ── Single overall bars ───────────────────────────────────────────
        y_pos = np.arange(n_taxa)
        ax.barh(
            y_pos, pct_vals,
            height=0.65,
            color=OVERALL_BAR_COLOR, alpha=0.88,
        )
        ax.set_yticks(y_pos)
        ax.set_yticklabels(
            [f"$\\it{{{t}}}$" if _is_binomial(t) else t for t in taxa],
            fontsize=FONT_SIZE_TICK,
        )
        title_suffix = "overall"

    # ── Axes formatting ───────────────────────────────────────────────────
    ax.invert_yaxis()
    ax.set_xlabel("Detection frequency (%)", fontsize=FONT_SIZE_AXIS)
    ax.set_xlim(0, min(105, ax.get_xlim()[1] * 1.08))
    ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f"{x:.0f}%")
    )
    ax.tick_params(axis="x", labelsize=FONT_SIZE_TICK)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5, zorder=0)

    ax.set_title(
        f"Detection frequency — {marker} ({title_suffix})",
        fontsize=FONT_SIZE_TITLE,
        pad=10,
    )

    n_shown = len(group_freq_df) if group_freq_df is not None else n_taxa
    fig.text(
        0.5, -0.02,
        f"Detection frequency among {sample_label}s with at least one detection. "
        f"Showing top {min(top_n, n_shown)} taxa.",
        ha="center", fontsize=7.5, color="#555555",
    )

    plt.tight_layout()

    # ── Save ─────────────────────────────────────────────────────────────
    group_slug = re.sub(r"[^\w]", "_", group_by).lower() if group_by else "overall"
    out_base   = outpath / f"detection_barplot_{marker}_{group_slug}"
    save_figure(fig, out_base, dry_run=dry_run)


def _is_binomial(name: str) -> bool:
    """
    Return True if name looks like a binomial species name
    (two words, first capitalised, second lower-case).
    Used to italicize species labels in plots.
    """
    parts = name.strip().split()
    if len(parts) != 2:
        return False
    return parts[0][0].isupper() and parts[1][0].islower()


# ---------------------------------------------------------------------------
# Infer taxonomic level from filename
# ---------------------------------------------------------------------------

def infer_level_from_path(path: Path) -> Optional[int]:
    """
    Try to extract the taxonomic level N from a filename like
    'taxonomy_counts_L7_cytb.tsv'. Returns None if not found.
    """
    m = re.search(r"_L(\d+)_", path.name)
    if m:
        return int(m.group(1))
    return None


# ---------------------------------------------------------------------------
# Main pipeline function
# ---------------------------------------------------------------------------

def run_presence_absence(
    counts_tsv: Path,
    marker: str,
    outdir: Path,
    metadata_path: Optional[Path],
    group_by: Optional[str],
    min_sample_reads: int,
    min_taxon_reads: int,
    min_relabund: float,
    top_n: int,
    force: bool,
    dry_run: bool,
    sample_label: str = "sample",
) -> None:
    """
    Full presence/absence pipeline:

      load counts -> filter -> convert to 0/1 -> compute detection
      frequencies -> write tables -> write summary -> write barplot.

    Args:
        counts_tsv:        Count TSV from 07_taxonomy_table.py.
        marker:            Marker name (for output naming only, e.g. 'cytb').
        outdir:            Directory to write all outputs.
        metadata_path:     QIIME2 metadata TSV (required if group_by is set).
        group_by:          Metadata column for group-level summaries/plots.
        min_sample_reads:  Drop samples with fewer total reads.
        min_taxon_reads:   Drop taxa with fewer reads across dataset.
        min_relabund:      Zero within-sample entries below this rel. abundance.
        top_n:             Max taxa to show in barplots.
        force:             Overwrite existing outputs.
        dry_run:           Log intended actions without writing files.
        sample_label:      Unit word for output text (e.g. "loon", "tick").
                           Defaults to "sample".
    """
    safe_mkdir(outdir)

    log.info("=== 11b_presence_absence: %s ===", marker)
    log.info("Input counts     : %s", counts_tsv)
    log.info("Output dir       : %s", outdir.resolve())
    log.info("Min sample reads : %d", min_sample_reads)
    log.info("Min taxon reads  : %d", min_taxon_reads)
    log.info("Min rel. abund   : %s",
             f"{min_relabund:.3f}" if min_relabund > 0 else "disabled")
    if group_by:
        log.info("Grouping column  : %s", group_by)
    if dry_run:
        log.info("DRY RUN — no files will be written")

    # ── Infer level from filename ─────────────────────────────────────────
    level = infer_level_from_path(counts_tsv)
    level_tag = f"L{level}_" if level is not None else ""
    log.info("Taxonomic level  : %s", str(level) if level else "unknown (not parsed from filename)")

    # ── 1. Load count table ───────────────────────────────────────────────
    log.info("── Loading count table ───────────────────────────────────────")
    counts = pd.read_csv(counts_tsv, sep="\t", index_col=0)
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    log.info(
        "  Loaded: %d taxa × %d samples  (%s total reads)",
        len(counts), len(counts.columns),
        f"{int(counts.values.sum()):,}",
    )

    # ── 2. Load metadata (optional) ───────────────────────────────────────
    groups: Optional[pd.Series] = None
    if metadata_path is not None and group_by is not None:
        log.info("── Loading metadata ──────────────────────────────────────────")
        groups = load_metadata(metadata_path, group_by)

    # ── 3. Filter ─────────────────────────────────────────────────────────
    log.info("── Applying filters ──────────────────────────────────────────")
    counts = apply_filters(counts, min_sample_reads, min_taxon_reads, min_relabund)

    # ── 4. Presence / absence ─────────────────────────────────────────────
    log.info("── Converting to presence/absence ────────────────────────────")
    pa = to_presence_absence(counts)

    # ── 5. Detection frequency ────────────────────────────────────────────
    log.info("── Computing detection frequencies ───────────────────────────")
    freq_df = compute_detection_freq(pa)

    group_freq_df: Optional[pd.DataFrame] = None
    if groups is not None:
        group_freq_df = compute_detection_freq_by_group(pa, groups)

    # ── 6. Write outputs ──────────────────────────────────────────────────
    log.info("── Writing outputs ───────────────────────────────────────────")

    # Binary presence/absence table
    pa_path = outdir / f"presence_absence_{level_tag}{marker}.tsv"
    save_tsv(pa, pa_path, dry_run)

    # Overall detection frequency table
    freq_path = outdir / f"detection_freq_{marker}.tsv"
    save_tsv(freq_df, freq_path, dry_run)

    # Per-group detection frequency table
    if group_freq_df is not None and group_by is not None:
        group_slug = re.sub(r"[^\w]", "_", group_by).lower()
        gfreq_path = outdir / f"detection_freq_by_{group_slug}_{marker}.tsv"
        save_tsv(group_freq_df, gfreq_path, dry_run)

    # Human-readable summary
    summary_text = build_summary_text(pa, freq_df, marker, groups, group_by,
                                      sample_label=sample_label)
    summary_path = outdir / f"detection_summary_{marker}.txt"
    save_text(summary_text, summary_path, dry_run)

    # Detection barplot
    log.info("── Generating detection barplot ──────────────────────────────")
    plot_detection_freq(
        freq_df       = freq_df,
        marker        = marker,
        outpath       = outdir,
        group_freq_df = group_freq_df,
        group_by      = group_by,
        top_n         = top_n,
        sample_label  = sample_label,
        dry_run       = dry_run,
    )

    # ── 7. Console summary ────────────────────────────────────────────────
    log.info("")
    log.info("=== Summary ===")
    n_total    = len(pa.columns)
    n_with_det = int((pa.sum(axis=0) > 0).sum())
    log.info(
        "%ss with detection : %d / %d  (%.1f%%)",
        sample_label.capitalize(),
        n_with_det, n_total,
        n_with_det / n_total * 100 if n_total > 0 else 0.0,
    )
    log.info("Top 5 taxa by detection frequency:")
    for taxon, row in freq_df.head(5).iterrows():
        log.info(
            "  %-45s  %d / %d  (%.1f%%)",
            taxon, int(row["n_detected"]), int(row["n_samples"]),
            float(row["pct_detected"]),
        )

    log.info("")
    log.info("Output files:")
    for f in sorted(outdir.iterdir()) if not dry_run else []:
        log.info("  %s", f)

    log.info("")
    log.info(
        "Next step — run Jaccard + PERMANOVA with 08_run_diversity_stats.py:\n"
        "  python 08_run_diversity_stats.py \\\n"
        "    --marker       %s \\\n"
        "    --dataset      all \\\n"
        "    --metadata     metadata/qiime/metadata_%s.tsv \\\n"
        "    --metrics-dir  <your core-metrics directory> \\\n"
        "    --group-column <your grouping column>\n\n"
        "  Note: 08_run_diversity_stats.py operates on QIIME2 core-metrics output.\n"
        "  The presence_absence TSV written here can be imported into R (vegan)\n"
        "  for Jaccard + PERMANOVA if preferred:\n\n"
        "    library(vegan)\n"
        "    pa <- read.table('%s', sep='\\t', header=TRUE, row.names=1)\n"
        "    jac <- vegdist(t(pa), method='jaccard', binary=TRUE)\n"
        "    meta <- read.table('metadata.tsv', sep='\\t', header=TRUE)\n"
        "    adonis2(jac ~ Group, data=meta)",
        marker, marker,
        outdir / f"presence_absence_{level_tag}{marker}.tsv",
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """
    Build and return the argument parser for 11b_presence_absence.py.

    Key arguments: --counts (TSV from 07_taxonomy_table.py), --marker,
    --outdir. Metadata and grouping are optional but required together.
    Filter thresholds have sensible defaults matching Zeb Antonioli's
    recommended approach for blood meal data.
    """
    p = argparse.ArgumentParser(
        prog="11b_presence_absence.py",
        description=(
            "Convert taxonomy count table to presence/absence and compute "
            "detection frequencies for amplicon metabarcoding data.\n"
            "Designed to run after 07_taxonomy_table.py.\n"
            "See module docstring for full rationale and usage examples."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required
    p.add_argument(
        "--counts", required=True, type=Path,
        help=(
            "Taxonomy count TSV from 07_taxonomy_table.py "
            "(e.g. taxonomy_counts_L7_cytb.tsv). "
            "Taxa as rows, samples as columns."
        ),
    )
    p.add_argument(
        "--marker", required=True,
        help="Marker name used for output file naming (e.g. cytb, COI).",
    )
    p.add_argument(
        "--outdir", required=True, type=Path,
        help="Directory to write all outputs.",
    )

    # Metadata / grouping (optional but must be used together)
    p.add_argument(
        "--metadata", default=None, type=Path,
        help=(
            "QIIME2-format metadata TSV. "
            "Required if --group-by is set."
        ),
    )
    p.add_argument(
        "--group-by", default=None,
        help=(
            "Metadata column to use for group-level detection summaries "
            "and grouped barplots (e.g. Site, TickSpecies)."
        ),
    )

    # Filter thresholds
    p.add_argument(
        "--min-sample-reads", type=int, default=500,
        metavar="N",
        help=(
            "Drop samples with fewer than N total reads before conversion. "
            "Default: 500"
        ),
    )
    p.add_argument(
        "--min-taxon-reads", type=int, default=50,
        metavar="N",
        help=(
            "Drop taxa with fewer than N reads across the whole dataset. "
            "Default: 50"
        ),
    )
    p.add_argument(
        "--min-relabund", type=float, default=0.0,
        metavar="F",
        help=(
            "Within each sample, zero out taxa whose relative read abundance "
            "is below F before converting to 0/1. "
            "Value is a fraction (e.g. 0.01 = 1%%). "
            "Default: 0.0 (disabled)"
        ),
    )

    # Plot options
    p.add_argument(
        "--top-n", type=int, default=20,
        metavar="N",
        help="Maximum number of taxa to show in barplots. Default: 20",
    )
    p.add_argument(
        "--sample-label", default="sample",
        metavar="WORD",
        help=(
            "Unit word used in summary text and figure subtitles. "
            "Set to match your study design, e.g. 'loon', 'tick', 'individual'. "
            "Default: sample"
        ),
    )

    # Execution flags
    p.add_argument(
        "--force", action="store_true", default=False,
        help="Overwrite existing outputs without prompting.",
    )
    p.add_argument(
        "--dry-run", action="store_true", default=False,
        help="Print planned actions without writing any files.",
    )

    return p


def main() -> int:
    """
    Parse arguments and run the presence/absence pipeline.

    Validates that --counts exists, that --metadata is provided when
    --group-by is set, and that min_relabund is in [0, 1]. Returns 0
    on success, 2 for missing/invalid inputs, 1 for any other error.
    """
    parser = build_parser()
    args   = parser.parse_args()

    # ── Input validation ──────────────────────────────────────────────────
    if not args.counts.exists():
        log.error("Count TSV not found: %s", args.counts)
        return 2

    if args.group_by is not None and args.metadata is None:
        log.error("--group-by requires --metadata to also be set.")
        return 2

    if args.metadata is not None and not args.metadata.exists():
        log.error("Metadata file not found: %s", args.metadata)
        return 2

    if not 0.0 <= args.min_relabund < 1.0:
        log.error(
            "--min-relabund must be between 0.0 and 1.0 (got %s)",
            args.min_relabund,
        )
        return 2

    # ── Run ───────────────────────────────────────────────────────────────
    try:
        run_presence_absence(
            counts_tsv       = args.counts,
            marker           = args.marker,
            outdir           = args.outdir,
            metadata_path    = args.metadata,
            group_by         = args.group_by,
            min_sample_reads = args.min_sample_reads,
            min_taxon_reads  = args.min_taxon_reads,
            min_relabund     = args.min_relabund,
            top_n            = args.top_n,
            force            = args.force,
            dry_run          = args.dry_run,
            sample_label     = args.sample_label,
        )
    except ValueError as e:
        log.error("%s", e)
        return 2
    except Exception as e:
        eprint(f"\n[ERROR] {e}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
