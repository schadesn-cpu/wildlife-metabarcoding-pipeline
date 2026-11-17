#!/usr/bin/env python3
"""
14c_plot_herpes_cutadapt.py
===========================
Generate herpesvirus detection figures from cutadapt-confirmed primer reads.

SCIENTIFIC RATIONALE
--------------------
This script replaces the earlier 10_plot_viral.py approach, which used DADA2
"no-hit" reads as a proxy for herpesvirus signal. That approach was invalid
because DADA2 discards reads that do not match the expected amplicon length
and quality, meaning "no-hit" reads were an artifact of the DADA2 pipeline,
not a direct measure of primer-confirmed herpesvirus amplicons.

This script uses cutadapt (v4.9) output directly. Cutadapt was run with
--discard-untrimmed, meaning only read pairs where BOTH forward (TGV) and
reverse (IYG) primers were found are retained. This is the methodologically
correct way to confirm that a read pair originated from a genuine herpesvirus
DPOL amplicon rather than from non-specific amplification or contamination.

BLAST validation on a subset of reads confirmed 91-97% nucleotide identity
to Gallid alphaherpesvirus 1 DPOL, confirming that primer-confirmed reads
are genuine herpesvirus amplicons.

THRESHOLD RATIONALE
-------------------
Detection is called at >= 1,000 primer-confirmed reads per sample. This
threshold was chosen because:
  1. The minimum positive sample has 5,097 reads -- well above threshold
  2. The distribution is bimodal: all 40 samples are clearly positive
  3. At any threshold from 1 to 5,097 reads, all 40/40 samples are positive
  4. The result is therefore threshold-independent for this dataset

PRIMERS
-------
Forward: TGV  5'-TGYAACTCGGTGTAYGGNTTYACNGGNGT-3' (VanDevanter et al. 1996)
Reverse: IYG  5'-CACMGAGTCCGTRTCNCCRTADAT-3'       (VanDevanter et al. 1996)
Target: Herpesvirus DNA polymerase (DPOL) gene, ~231 bp amplicon

USAGE
-----
    python scripts/10b_plot_herpes_cutadapt.py \
        --cutadapt results/herpesvirus/cutadapt/cutadapt_summary.txt \
        --metadata metadata/qiime/metadata_16S_updated.tsv \
        --group-by Group \
        --group-order Diseased Trauma \
        --threshold 1000 \
        --palette wong \
        --outdir results/herpesvirus/figures/

OUTPUT
------
    herpes_cutadapt_{group}_{palette}_presence.png/.svg  -- binary detection
    herpes_cutadapt_{group}_{palette}_relabund.png/.svg  -- read proportions

Author: Samantha Schade | MEED Lab, UNH | 2026-04-04
"""

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# ── Color palettes ─────────────────────────────────────────────────────────────
# Wong (2011) colorblind-safe palette -- used throughout this project for
# consistency across all figures. Colors chosen to be distinguishable by
# individuals with the most common forms of color vision deficiency.
WONG_COLORS = [
    "#0072B2",  # blue       -- Diseased primary
    "#E69F00",  # orange     -- Trauma primary
    "#009E73",  # green      -- Marine / third group
    "#CC79A7",  # pink
    "#56B4E9",  # light blue
    "#D55E00",  # vermillion
    "#F0E442",  # yellow
    "#000000",  # black
]

# Grey used for non-primer-confirmed reads in the relative abundance figure.
# This is scientifically meaningful: grey = reads that did NOT contain both
# TGV and IYG primers and therefore cannot be attributed to herpesvirus DPOL.
NON_PRIMER_COLOR = "#CCCCCC"


# ── Cutadapt parsing ───────────────────────────────────────────────────────────

def parse_cutadapt_summary(summary_path: Path) -> dict:
    """
    Parse a cutadapt summary file and return per-sample read counts.

    Cutadapt writes one '=== Summary ===' block per sample in the order
    the samples were processed. We match each block to its sample ID by
    looking at the R1 fastq files in the same directory, sorted
    alphabetically -- which is how cutadapt was called in this pipeline.

    Parameters
    ----------
    summary_path : Path
        Path to cutadapt_summary.txt

    Returns
    -------
    dict mapping TV-ID -> {'total': int, 'kept': int, 'pct': float}
        total = total read pairs input to cutadapt
        kept  = read pairs where both primers were found (--discard-untrimmed)
        pct   = kept / total * 100

    Raises
    ------
    FileNotFoundError : if summary_path does not exist
    ValueError : if no samples could be parsed from the summary file
    ValueError : if R1 fastq files and summary blocks are mismatched in count
    """
    summary_path = Path(summary_path)

    if not summary_path.exists():
        raise FileNotFoundError(
            f"Cutadapt summary not found: {summary_path}\n"
            f"Expected output from: cutadapt --discard-untrimmed"
        )

    # Read the full summary file
    try:
        raw = summary_path.read_text()
    except OSError as e:
        raise OSError(f"Could not read cutadapt summary file: {e}") from e

    # Split into per-sample blocks on the summary header
    blocks = raw.split("=== Summary ===")
    data_blocks = blocks[1:]  # first element is the file header, skip it

    if not data_blocks:
        raise ValueError(
            f"No '=== Summary ===' blocks found in {summary_path}. "
            f"Is this a valid cutadapt summary file?"
        )

    # Find R1 fastq files in the same directory as the summary to get
    # sample IDs in the same order cutadapt processed them.
    cutadapt_dir = summary_path.parent
    r1_files = sorted(cutadapt_dir.glob("*_R1.fastq.gz"))

    if not r1_files:
        raise FileNotFoundError(
            f"No *_R1.fastq.gz files found in {cutadapt_dir}. "
            f"Cannot determine sample order for parsing."
        )

    if len(r1_files) != len(data_blocks):
        raise ValueError(
            f"Mismatch: found {len(r1_files)} R1 fastq files but "
            f"{len(data_blocks)} summary blocks in {summary_path}. "
            f"Summary file may be incomplete or from a different run."
        )

    # Extract TV-IDs from filenames (first hyphen-delimited field)
    # e.g. TV230007-lung-TGF-IYG_S1581_R1.fastq.gz -> TV230007
    sample_ids = [f.name.split("-")[0] for f in r1_files]

    # Parse each block for total and retained read counts
    results = {}
    for i, block in enumerate(data_blocks):
        tv = sample_ids[i]

        # Regex patterns match cutadapt's standard summary output format
        total_match = re.search(
            r"Total read pairs processed:\s+([\d,]+)", block
        )
        kept_match = re.search(
            r"Pairs written \(passing filters\):\s+([\d,]+)", block
        )

        if not total_match or not kept_match:
            # A missing match means the block is malformed -- report which
            # sample so the user can inspect the raw summary file.
            raise ValueError(
                f"Could not parse read counts for sample {tv} "
                f"(block {i+1} of {len(data_blocks)}). "
                f"Check {summary_path} for truncation or formatting errors."
            )

        total = int(total_match.group(1).replace(",", ""))
        kept = int(kept_match.group(1).replace(",", ""))

        if total == 0:
            # A sample with zero total reads indicates a failed library --
            # flag it explicitly rather than producing a division-by-zero.
            print(
                f"  WARNING: {tv} has 0 total reads in cutadapt summary. "
                f"Library may have failed. Reporting pct=0.",
                file=sys.stderr,
            )
            pct = 0.0
        else:
            pct = kept / total * 100

        results[tv] = {"total": total, "kept": kept, "pct": pct}

    return results


# ── Metadata loading ───────────────────────────────────────────────────────────

def load_group_metadata(meta_path: Path, group_col: str, group_order: list) -> dict:
    """
    Load QIIME2-format metadata and return a TV-ID to group mapping.

    QIIME2 metadata files have a comment line starting with '#q2:types'
    as the second row, which must be skipped. The first column is the
    sample-id; we extract the TV-ID (e.g. TV230007) from it because
    the cutadapt output uses TV-IDs while QIIME2 metadata uses full
    sample IDs (e.g. TV230007-GI-16S).

    Parameters
    ----------
    meta_path : Path
        Path to QIIME2 metadata TSV
    group_col : str
        Column name to use for grouping (e.g. 'Group', 'COD_broad')
    group_order : list of str
        Groups to include; samples in other groups are excluded

    Returns
    -------
    dict mapping TV-ID -> group label (str)

    Raises
    ------
    FileNotFoundError : if meta_path does not exist
    KeyError : if group_col is not found in the metadata
    ValueError : if no samples match any group in group_order
    """
    meta_path = Path(meta_path)

    if not meta_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {meta_path}")

    try:
        # comment='#' skips the QIIME2 '#q2:types' directive line
        df = pd.read_csv(meta_path, sep="\t", comment="#")
    except Exception as e:
        raise OSError(f"Could not read metadata file {meta_path}: {e}") from e

    if group_col not in df.columns:
        available = ", ".join(df.columns.tolist())
        raise KeyError(
            f"Column '{group_col}' not found in {meta_path}. "
            f"Available columns: {available}"
        )

    id_col = df.columns[0]  # first column is always the sample-id

    # Extract TV-ID from sample-id. The TV-ID is always the first
    # hyphen-delimited field (e.g. 'TV230007-GI-16S' -> 'TV230007').
    # Use a regex to be robust to varying sample-id formats.
    df["TV"] = df[id_col].str.extract(r"(TV\d+)", expand=False)

    # Drop rows where TV-ID or group could not be extracted
    n_before = len(df)
    df = df.dropna(subset=["TV", group_col])
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        print(
            f"  NOTE: Dropped {n_dropped} metadata rows with missing TV-ID "
            f"or '{group_col}' value.",
            file=sys.stderr,
        )

    # Filter to requested groups only
    df = df[df[group_col].isin(group_order)]

    if df.empty:
        raise ValueError(
            f"No samples found in metadata matching groups: {group_order}. "
            f"Check that '{group_col}' column contains these values."
        )

    # Build TV -> group dict. If the same TV appears multiple times
    # (e.g. same bird in multiple marker tables), keep the first occurrence.
    mapping = df.drop_duplicates(subset="TV").set_index("TV")[group_col].to_dict()

    return mapping


# ── Data assembly ──────────────────────────────────────────────────────────────

def build_sample_table(
    cutadapt_counts: dict,
    group_mapping: dict,
    group_order: list,
    threshold: int,
) -> pd.DataFrame:
    """
    Join cutadapt counts with group metadata and apply detection threshold.

    Samples that appear in cutadapt output but not in metadata are reported
    as warnings and excluded from figures. Samples in metadata but not in
    cutadapt output (e.g. TV240083 which had 0 reads in the sequencing run)
    are also reported.

    Parameters
    ----------
    cutadapt_counts : dict
        Output of parse_cutadapt_summary()
    group_mapping : dict
        Output of load_group_metadata()
    group_order : list of str
        Groups to include, in display order
    threshold : int
        Minimum primer-confirmed reads to call a sample positive

    Returns
    -------
    pd.DataFrame with columns:
        tv, group, total, kept, pct, detected
    Sorted by group (per group_order) then by kept reads descending.
    """
    rows = []
    unmatched_cutadapt = []  # in cutadapt but not metadata
    unmatched_meta = []       # in metadata but not cutadapt

    for tv, vals in cutadapt_counts.items():
        grp = group_mapping.get(tv)
        if grp is None:
            unmatched_cutadapt.append(tv)
            continue
        if grp not in group_order:
            # Sample is in metadata but not in the requested group subset --
            # e.g. Marine birds excluded from DvT analysis. Not a warning.
            continue
        rows.append({
            "tv": tv,
            "group": grp,
            "total": vals["total"],
            "kept": vals["kept"],
            "pct": vals["pct"],
            "detected": vals["kept"] >= threshold,
        })

    # Check for metadata samples with no cutadapt output
    for tv in group_mapping:
        if tv not in cutadapt_counts and group_mapping[tv] in group_order:
            unmatched_meta.append(tv)

    if unmatched_cutadapt:
        print(
            f"  WARNING: {len(unmatched_cutadapt)} sample(s) in cutadapt output "
            f"have no metadata match and will be excluded: "
            f"{', '.join(sorted(unmatched_cutadapt))}",
            file=sys.stderr,
        )
    if unmatched_meta:
        print(
            f"  WARNING: {len(unmatched_meta)} sample(s) in metadata have no "
            f"cutadapt output (may have failed sequencing): "
            f"{', '.join(sorted(unmatched_meta))}",
            file=sys.stderr,
        )

    if not rows:
        raise ValueError(
            "No samples remained after joining cutadapt output with metadata. "
            "Check that TV-IDs match between the two sources."
        )

    df = pd.DataFrame(rows)

    # Sort by group order, then by kept reads descending within each group
    # (highest-read samples on the left within each group)
    group_order_map = {g: i for i, g in enumerate(group_order)}
    df["group_rank"] = df["group"].map(group_order_map)
    df = df.sort_values(["group_rank", "kept"], ascending=[True, False])
    df = df.drop(columns="group_rank").reset_index(drop=True)

    return df


# ── Plotting ───────────────────────────────────────────────────────────────────

def plot_presence_absence(
    df: pd.DataFrame,
    group_order: list,
    color_map: dict,
    threshold: int,
    outdir: Path,
    stem: str,
) -> None:
    """
    Presence/absence barplot: filled = detected, hatched = not detected.

    Each bar represents one bird. Bars are grouped by mortality group with
    a dashed vertical divider. Group-level detection rates (n/N) are
    annotated above each group. For two-group comparisons, Fisher's exact
    test p-value is shown in the bottom-right corner.

    The threshold is displayed on the figure so readers know exactly what
    criterion was used to call a sample positive.

    Parameters
    ----------
    df : pd.DataFrame
        Output of build_sample_table()
    group_order : list of str
        Groups in display order
    color_map : dict
        Maps group label -> hex color string
    threshold : int
        Reads threshold used for detection (display only, already applied)
    outdir : Path
        Output directory
    stem : str
        Output filename stem (no extension)
    """
    n_samples = len(df)
    fig_width = max(8, n_samples * 0.5 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    x_positions = []
    x_labels = []
    group_x_positions = {g: [] for g in group_order}

    x = 0
    for grp in group_order:
        grp_df = df[df["group"] == grp]

        for _, row in grp_df.iterrows():
            if row["detected"]:
                # Detected: solid fill in group color
                ax.bar(
                    x, 1,
                    color=color_map[grp],
                    edgecolor="white",
                    linewidth=0.5,
                    zorder=2,
                )
            else:
                # Not detected: white fill with hatching and colored border
                # Using hatching rather than grey fill makes the not-detected
                # status unambiguous even in greyscale printing.
                ax.bar(
                    x, 1,
                    color="white",
                    edgecolor=color_map[grp],
                    linewidth=1.5,
                    hatch="///",
                    zorder=2,
                )

            x_positions.append(x)
            # Display just the numeric portion of the TV-ID for readability
            x_labels.append(row["tv"].replace("TV", ""))
            group_x_positions[grp].append(x)
            x += 1

        # Add gap between groups for visual separation
        x += 0.8

    # Annotate each group with detection rate and a divider line
    for i, grp in enumerate(group_order):
        positions = group_x_positions[grp]
        if not positions:
            continue

        midpoint = np.mean(positions)
        n_detected = df[df["group"] == grp]["detected"].sum()
        n_total = len(positions)

        ax.text(
            midpoint, 1.06,
            f"{grp}\n{n_detected}/{n_total} positive",
            ha="center", va="bottom",
            fontsize=10, fontweight="bold",
            color=color_map[grp],
        )

        # Dashed divider after all groups except the last
        if i < len(group_order) - 1:
            divider_x = max(positions) + 0.4
            ax.axvline(
                divider_x,
                color="#cccccc",
                linestyle="--",
                linewidth=1,
                zorder=1,
            )

    # Fisher's exact test for exactly two groups.
    # We only run this for two-group comparisons because the test is designed
    # for a 2x2 contingency table. For three or more groups, a chi-square
    # or Fisher-Freeman-Halton test would be needed instead.
    if len(group_order) == 2:
        g1_df = df[df["group"] == group_order[0]]
        g2_df = df[df["group"] == group_order[1]]
        g1_pos = g1_df["detected"].sum()
        g2_pos = g2_df["detected"].sum()
        g1_neg = len(g1_df) - g1_pos
        g2_neg = len(g2_df) - g2_pos

        _, p_val = fisher_exact([[g1_pos, g1_neg], [g2_pos, g2_neg]])

        star = " *" if p_val <= 0.05 else " ns"
        annot_color = "#444444"

        ax.annotate(
            f"Fisher's exact  p={p_val:.3f}{star}",
            xy=(0.98, 0.04),
            xycoords="axes fraction",
            ha="right", va="bottom",
            fontsize=9, style="italic",
            color=annot_color,
        )
    elif len(group_order) > 2:
        # Remind the reader that pairwise tests need FDR correction when
        # there are more than two groups
        ax.annotate(
            "Multiple groups: use pairwise Fisher's exact with FDR correction",
            xy=(0.98, 0.04),
            xycoords="axes fraction",
            ha="right", va="bottom",
            fontsize=8, style="italic",
            color="#666666",
        )

    # Threshold annotation in bottom-left so readers know the detection criterion
    ax.annotate(
        f"Threshold: \u2265{threshold:,} primer-confirmed reads (TGV/IYG; VanDevanter et al. 1996)",
        xy=(0.02, 0.04),
        xycoords="axes fraction",
        ha="left", va="bottom",
        fontsize=8, style="italic",
        color="#666666",
    )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=90, fontsize=7)
    ax.set_ylim(0, 1.3)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Not detected", "Detected"], fontsize=10)
    ax.set_ylabel(
        "Herpesvirus DPOL detection\n(cutadapt primer-confirmed reads)",
        fontsize=10,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    try:
        out_png = outdir / f"{stem}_presence.png"
        out_svg = outdir / f"{stem}_presence.svg"
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        fig.savefig(out_svg, bbox_inches="tight")
        print(f"  \u2713 Saved: {out_png}")
        print(f"  \u2713 Saved: {out_svg}")
    except OSError as e:
        raise OSError(f"Could not save presence/absence figure: {e}") from e
    finally:
        plt.close(fig)


def plot_relative_abundance(
    df: pd.DataFrame,
    group_order: list,
    color_map: dict,
    outdir: Path,
    stem: str,
) -> None:
    """
    Relative abundance barplot: primer-confirmed reads vs. total reads.

    Each bar represents one bird. The colored portion = primer-confirmed reads
    (i.e. reads containing both TGV and IYG primers, attributed to herpesvirus
    DPOL). The grey portion = all other reads (non-primer reads that did not
    survive cutadapt --discard-untrimmed filtering).

    This figure shows the MAGNITUDE of herpesvirus signal per bird, not just
    the binary positive/negative. High primer retention (>70%) across all
    samples indicates the primers were highly specific to this amplicon.

    Parameters
    ----------
    df : pd.DataFrame
        Output of build_sample_table()
    group_order : list of str
        Groups in display order
    color_map : dict
        Maps group label -> hex color string
    outdir : Path
        Output directory
    stem : str
        Output filename stem (no extension)
    """
    n_samples = len(df)
    fig_width = max(8, n_samples * 0.5 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    x_positions = []
    x_labels = []
    group_x_positions = {g: [] for g in group_order}

    x = 0
    for grp in group_order:
        grp_df = df[df["group"] == grp]

        for _, row in grp_df.iterrows():
            # Proportion of reads that are primer-confirmed herpesvirus DPOL
            prop_herpes = row["kept"] / row["total"] if row["total"] > 0 else 0
            prop_other = 1.0 - prop_herpes

            # Herpesvirus portion (bottom of stack)
            ax.bar(
                x, prop_herpes,
                color=color_map[grp],
                edgecolor="white",
                linewidth=0.3,
                zorder=2,
                label=grp if x == list(df[df["group"]==grp].index)[0] else "",
            )
            # Non-primer reads (top of stack, grey)
            ax.bar(
                x, prop_other,
                bottom=prop_herpes,
                color=NON_PRIMER_COLOR,
                edgecolor="white",
                linewidth=0.3,
                zorder=2,
            )

            x_positions.append(x)
            x_labels.append(row["tv"].replace("TV", ""))
            group_x_positions[grp].append(x)
            x += 1

        x += 0.8  # gap between groups

    # Group labels and dividers
    for i, grp in enumerate(group_order):
        positions = group_x_positions[grp]
        if not positions:
            continue

        midpoint = np.mean(positions)
        n_total = len(positions)
        grp_df = df[df["group"] == grp]
        median_pct = grp_df["pct"].median()

        ax.text(
            midpoint, 1.02,
            f"{grp} (n={n_total})\nmedian {median_pct:.0f}% primer-confirmed",
            ha="center", va="bottom",
            fontsize=9, fontweight="bold",
            color=color_map[grp],
        )

        if i < len(group_order) - 1:
            ax.axvline(
                max(positions) + 0.4,
                color="#cccccc",
                linestyle="--",
                linewidth=1,
                zorder=1,
            )

    # Legend distinguishing herpesvirus reads from non-primer reads
    legend_handles = [
        mpatches.Patch(color=color_map[g], label=f"{g} — herpesvirus DPOL reads")
        for g in group_order
    ] + [
        mpatches.Patch(
            color=NON_PRIMER_COLOR,
            label="Non-primer reads (not herpesvirus-attributed)",
        )
    ]
    ax.legend(
        handles=legend_handles,
        loc="lower right",
        fontsize=8,
        framealpha=0.85,
        edgecolor="#cccccc",
    )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=90, fontsize=7)
    ax.set_ylim(0, 1.25)
    ax.set_ylabel(
        "Proportion of reads\n(primer-confirmed / total input reads)",
        fontsize=10,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotate primer source in bottom-left for reproducibility
    ax.annotate(
        "Primers: TGV/IYG (VanDevanter et al. 1996) | cutadapt --discard-untrimmed",
        xy=(0.02, 0.04),
        xycoords="axes fraction",
        ha="left", va="bottom",
        fontsize=7.5, style="italic",
        color="#666666",
    )

    plt.tight_layout()

    try:
        out_png = outdir / f"{stem}_relabund.png"
        out_svg = outdir / f"{stem}_relabund.svg"
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        fig.savefig(out_svg, bbox_inches="tight")
        print(f"  \u2713 Saved: {out_png}")
        print(f"  \u2713 Saved: {out_svg}")
    except OSError as e:
        raise OSError(f"Could not save relative abundance figure: {e}") from e
    finally:
        plt.close(fig)


# ── Entry point ────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate herpesvirus detection figures from cutadapt-confirmed "
            "primer reads. See script header for full scientific rationale."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--cutadapt", required=True,
        help="Path to cutadapt_summary.txt from the herpesvirus cutadapt run",
    )
    parser.add_argument(
        "--metadata", required=True,
        help="Path to QIIME2-format metadata TSV with group information",
    )
    parser.add_argument(
        "--group-by", default="Group",
        help="Metadata column to use for grouping (default: Group)",
    )
    parser.add_argument(
        "--group-order", nargs="+", default=["Diseased", "Trauma"],
        help="Groups to display, in order (default: Diseased Trauma)",
    )
    parser.add_argument(
        "--threshold", type=int, default=1000,
        help=(
            "Minimum primer-confirmed reads to call a sample positive "
            "(default: 1000). For this dataset all 40 samples are positive "
            "at any threshold >= 1 read."
        ),
    )
    parser.add_argument(
        "--palette", choices=["wong"], default="wong",
        help="Color palette (default: wong -- Wong 2011 colorblind-safe)",
    )
    parser.add_argument(
        "--outdir", default="results/herpesvirus/figures/",
        help="Output directory for figures (default: results/herpesvirus/figures/)",
    )
    parser.add_argument(
        "--stem", default="herpes_cutadapt_Group_wong",
        help="Output filename stem, no extension (default: herpes_cutadapt_Group_wong)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(f"ERROR: Could not create output directory {outdir}: {e}", file=sys.stderr)
        sys.exit(1)

    # Assign colors to groups in the order they were specified
    color_map = {
        grp: WONG_COLORS[i % len(WONG_COLORS)]
        for i, grp in enumerate(args.group_order)
    }

    # ── Step 1: Parse cutadapt output ────────────────────────────────────────
    print(f"\nStep 1: Parsing cutadapt summary")
    print(f"  File: {args.cutadapt}")
    try:
        counts = parse_cutadapt_summary(Path(args.cutadapt))
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"  Parsed {len(counts)} samples")
    kept_vals = [v["kept"] for v in counts.values()]
    print(f"  Primer-confirmed reads range: {min(kept_vals):,} -- {max(kept_vals):,}")
    print(f"  Median: {int(np.median(kept_vals)):,}")

    # ── Step 2: Load metadata ─────────────────────────────────────────────────
    print(f"\nStep 2: Loading metadata")
    print(f"  File: {args.metadata}")
    print(f"  Grouping by: {args.group_by}")
    try:
        group_mapping = load_group_metadata(
            Path(args.metadata), args.group_by, args.group_order
        )
    except (FileNotFoundError, KeyError, ValueError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"  Found {len(group_mapping)} samples in requested groups")

    # ── Step 3: Build sample table ────────────────────────────────────────────
    print(f"\nStep 3: Joining data (threshold >= {args.threshold:,} reads)")
    try:
        df = build_sample_table(
            counts, group_mapping, args.group_order, args.threshold
        )
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"  Final sample count: {len(df)}")
    for grp in args.group_order:
        grp_df = df[df["group"] == grp]
        n_pos = grp_df["detected"].sum()
        print(f"    {grp}: {n_pos}/{len(grp_df)} positive")

    # ── Step 4: Generate figures ──────────────────────────────────────────────
    print(f"\nStep 4: Generating figures")
    print(f"  Output directory: {outdir}")

    try:
        plot_presence_absence(
            df, args.group_order, color_map,
            args.threshold, outdir, args.stem,
        )
    except OSError as e:
        print(f"ERROR generating presence/absence figure: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        plot_relative_abundance(
            df, args.group_order, color_map,
            outdir, args.stem,
        )
    except OSError as e:
        print(f"ERROR generating relative abundance figure: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\nDone.")


if __name__ == "__main__":
    main()
