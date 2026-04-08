#!/usr/bin/env python3
"""
10c_plot_adeno_cutadapt.py
==========================
Generate aviadenovirus detection figures from cutadapt-confirmed primer reads.

SCIENTIFIC RATIONALE
--------------------
This script generates presence/absence and relative abundance figures for
aviadenovirus detection in Common Loon lung tissue, using reads confirmed
by cutadapt primer trimming rather than DADA2 classification.

The adenovirus pipeline uses polF/polR consensus DPOL primers (Wellehan et al.
2004), which target the conserved DNA polymerase gene of adenoviruses. Cutadapt
was run with --discard-untrimmed, retaining only read pairs where BOTH forward
(polF) and reverse (polR) primers were found. This is the methodologically
correct approach because it confirms that retained reads originated from a
genuine adenovirus DPOL amplicon rather than non-specific amplification.

BLAST analysis of the two recovered OTUs showed 76.8% (OTU1, 9,123 reads) and
76.6% (OTU2, 128 reads) nucleotide identity to the closest known sequences
(Aviadenovirus sp. YN06/YN10), below the ~85% species delineation threshold,
supporting a putative novel aviadenovirus species in Common Loons.

THRESHOLD RATIONALE
-------------------
Unlike herpesvirus (where all 40 samples are positive at any threshold),
adenovirus detection is genuinely threshold-dependent. The read count
distribution shows a continuous range from 0 to 10,508 reads with no
obvious bimodal gap. Two thresholds are scientifically defensible:

  >= 500 reads:  Diseased 7/19, Trauma 2/16, Fisher p=0.135 ns
                 Parasitic_Infectious 6/11 vs Trauma 2/16, Fisher p=0.033 *
                 -- more sensitive, captures more positives

  >= 1000 reads: Diseased 5/19, Trauma 1/16, Fisher p=0.187 ns
                 Parasitic_Infectious 4/11 vs Trauma 1/16, Fisher p=0.125 ns
                 -- more conservative, result is no longer significant

Because the PI-vs-Trauma finding (p=0.033) is threshold-sensitive, figures are
generated at both thresholds. The threshold-sensitivity is a known limitation
and is explicitly noted in the figure annotations and manuscript.

This script generates four figures:
  (1) DvT presence/absence       -- Diseased vs Trauma
  (2) DvT relative abundance     -- Diseased vs Trauma
  (3) COD presence/absence       -- Lead, Parasitic_Infectious, Trauma
  (4) COD relative abundance     -- Lead, Parasitic_Infectious, Trauma

PRIMERS
-------
Forward: polF  5'-GTNTWYGAYATHTGYGGHATGTAYGC-3'  (Wellehan et al. 2004)
Reverse: polR  5'-CCANCCBCDRTTRTGNARNGTRA-3'      (Wellehan et al. 2004)
Target: Aviadenovirus DNA polymerase (DPOL) gene

USAGE
-----
    python scripts/10c_plot_adeno_cutadapt.py \
        --cutadapt results/adenovirus/cutadapt/cutadapt_summary.txt \
        --metadata-dvt metadata/qiime/metadata_16S_updated.tsv \
        --metadata-cod metadata/qiime/metadata_16S_updated.tsv \
        --group-by-dvt Group \
        --group-by-cod COD_broad \
        --threshold 500 \
        --palette wong \
        --outdir results/adenovirus/figures/

OUTPUT
------
    adeno_cutadapt_DvT_wong_{threshold}reads_presence.png/.svg
    adeno_cutadapt_DvT_wong_{threshold}reads_relabund.png/.svg
    adeno_cutadapt_COD_wong_{threshold}reads_presence.png/.svg
    adeno_cutadapt_COD_wong_{threshold}reads_relabund.png/.svg

Author: Samantha Schade | MEED Lab, UNH | 2026-04-04
"""

import argparse
import glob
import re
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# ── Color palettes ─────────────────────────────────────────────────────────────
# Wong (2011) colorblind-safe palette used throughout this project.
# Color assignments follow the convention established in 06_plot_diversity.py:
#   DvT:  Diseased = blue (#0072B2), Trauma = orange (#E69F00)
#   COD:  Lead = blue, Parasitic_Infectious = orange, Trauma = green
WONG_COLORS = [
    "#0072B2",  # blue       -- Diseased / Lead
    "#E69F00",  # orange     -- Trauma (DvT) / Parasitic_Infectious (COD)
    "#009E73",  # green      -- Trauma (COD)
    "#CC79A7",  # pink
    "#56B4E9",  # light blue
    "#D55E00",  # vermillion
    "#F0E442",  # yellow
    "#000000",  # black
]

# Grey for reads that did not survive --discard-untrimmed filtering.
# Scientifically: these reads did not contain both polF and polR primers
# and therefore cannot be attributed to adenovirus DPOL amplification.
NON_PRIMER_COLOR = "#CCCCCC"

# Purple used for significant Fisher's exact results, consistent with
# PERMANOVA annotation color used in 06_plot_diversity.py.
SIG_COLOR = "#6A0DAD"
NS_COLOR = "#444444"


# ── Cutadapt parsing ───────────────────────────────────────────────────────────

def parse_cutadapt_summary(summary_path: Path) -> dict:
    """
    Parse a cutadapt summary file and return per-sample primer-confirmed
    read counts.

    Cutadapt writes one '=== Summary ===' block per sample in the order
    samples were processed. We recover sample IDs by matching blocks to
    R1 fastq files in the same directory, sorted alphabetically -- which
    is how cutadapt was called in this pipeline.

    Parameters
    ----------
    summary_path : Path
        Path to cutadapt_summary.txt from the adenovirus cutadapt run.

    Returns
    -------
    dict mapping TV-ID (str) -> {
        'total': int,   # total read pairs input to cutadapt
        'kept':  int,   # read pairs where both polF and polR primers found
        'pct':   float  # kept / total * 100
    }

    Raises
    ------
    FileNotFoundError
        If summary_path does not exist.
    FileNotFoundError
        If no R1 fastq files are found in the same directory as the summary.
    ValueError
        If the number of R1 files does not match the number of summary blocks,
        indicating a truncated or mismatched summary file.
    ValueError
        If a summary block is missing the expected read count fields,
        indicating file corruption or an unexpected cutadapt version.
    """
    summary_path = Path(summary_path)

    if not summary_path.exists():
        raise FileNotFoundError(
            f"Cutadapt summary not found: {summary_path}\n"
            f"Expected output from: cutadapt --discard-untrimmed "
            f"using polF/polR primers (Wellehan et al. 2004)"
        )

    try:
        raw = summary_path.read_text()
    except OSError as e:
        raise OSError(
            f"Could not read cutadapt summary {summary_path}: {e}"
        ) from e

    # Split on the cutadapt version header rather than "=== Summary ==="
    # because the Command line parameters line appears BEFORE each Summary
    # block, meaning splitting on Summary puts the Command line at the end
    # of the previous block rather than the start of the current one.
    # Splitting on "This is cutadapt" keeps Command line and Summary together.
    blocks = raw.split("This is cutadapt")
    data_blocks = blocks[1:]  # first element is file content before first run

    if not data_blocks:
        raise ValueError(
            f"No '=== Summary ===' blocks found in {summary_path}. "
            f"Is this a valid cutadapt output file?"
        )

    # Match each block to its sample by extracting the TV-ID from the -o
    # argument in the "Command line parameters" line of each block.
    # This is more robust than position-based matching because cutadapt
    # skips entirely empty input files (e.g. TV240083, whose fastq was
    # 20 bytes / effectively empty) and produces no summary block for them.
    # Position-based matching would silently shift all subsequent sample
    # assignments by one, corrupting all downstream results without error.
    #
    # Example Command line parameters line:
    #   -o results/adenovirus/cutadapt/TV230007-lung-..._R1.fastq.gz
    # Regex captures the TV-ID directly from the output path.

    results = {}
    for i, block in enumerate(data_blocks):
        cmd_match = re.search(r"-o\s+\S*?(TV\d+)", block)
        if not cmd_match:
            raise ValueError(
                f"Could not extract TV-ID from summary block {i + 1} of "
                f"{len(data_blocks)} in {summary_path}. "
                f"Expected a '-o .../TV######-...' argument in the "
                f"'Command line parameters' line of each block. "
                f"Check that this is a valid cutadapt summary file."
            )
        tv = cmd_match.group(1)

        total_match = re.search(
            r"Total read pairs processed:\s+([\d,]+)", block
        )
        kept_match = re.search(
            r"Pairs written \(passing filters\):\s+([\d,]+)", block
        )

        if not total_match or not kept_match:
            # No Summary block means cutadapt ran but found no reads --
            # "No reads processed!" case. Treat as 0 reads rather than
            # raising an error, since this is a known library failure mode.
            if "No reads processed" in block:
                print(
                    f"  NOTE: {tv} produced 'No reads processed' in cutadapt "
                    f"(empty input fastq). Treating as 0 reads / not detected.",
                    file=sys.stderr,
                )
                results[tv] = {"total": 0, "kept": 0, "pct": 0.0}
                continue
            raise ValueError(
                f"Could not parse read counts for sample {tv} "
                f"(block {i + 1} of {len(data_blocks)}). "
                f"Check {summary_path} for truncation or formatting errors. "
                f"Expected lines matching:\n"
                f"  'Total read pairs processed: N'\n"
                f"  'Pairs written (passing filters): N'"
            )

        total = int(total_match.group(1).replace(",", ""))
        kept = int(kept_match.group(1).replace(",", ""))

        if total == 0:
            # Zero total reads from a valid Summary block -- different from
            # "No reads processed" but also indicates a failed library.
            print(
                f"  WARNING: {tv} has 0 total reads in cutadapt summary. "
                f"Library may have failed. This sample will appear as "
                f"not detected regardless of threshold.",
                file=sys.stderr,
            )
            pct = 0.0
        else:
            pct = kept / total * 100

        results[tv] = {"total": total, "kept": kept, "pct": pct}

    # Report R1 files that had no summary block -- these were skipped by
    # cutadapt because the input fastq was empty (e.g. TV240083, 20 bytes).
    # We add them as 0-read entries so downstream code sees them as
    # not-detected rather than silently missing from the analysis.
    cutadapt_dir = summary_path.parent
    r1_tvids = set(
        f.name.split("-")[0]
        for f in cutadapt_dir.glob("*_R1.fastq.gz")
    )
    missing = r1_tvids - set(results.keys())
    if missing:
        print(
            f"  NOTE: {len(missing)} sample(s) had R1 files but no cutadapt "
            f"summary block (likely empty input fastq, skipped by cutadapt): "
            f"{', '.join(sorted(missing))}. "
            f"Treating as 0 reads / not detected.",
            file=sys.stderr,
        )
        for tv in sorted(missing):
            results[tv] = {"total": 0, "kept": 0, "pct": 0.0}

    return results


# ── Metadata loading ───────────────────────────────────────────────────────────

def load_group_metadata(
    meta_path: Path,
    group_col: str,
    group_order: list,
) -> dict:
    """
    Load QIIME2-format metadata and return a TV-ID to group label mapping.

    QIIME2 metadata files may contain a '#q2:types' directive line as the
    second row, which is skipped via comment='#'. The first column is always
    the sample-id; we extract the TV-ID from it because cutadapt output uses
    TV-IDs while QIIME2 sample-ids include marker suffixes
    (e.g. 'TV230007-GI-16S' -> 'TV230007').

    Parameters
    ----------
    meta_path : Path
        Path to QIIME2-format metadata TSV.
    group_col : str
        Column name to group samples by (e.g. 'Group', 'COD_broad').
    group_order : list of str
        Which groups to include. Samples in other groups are silently excluded
        (e.g. Marine birds excluded from DvT comparison).

    Returns
    -------
    dict mapping TV-ID (str) -> group label (str)
        If the same TV-ID appears in multiple rows (e.g. same bird in multiple
        marker tables), the first occurrence is used.

    Raises
    ------
    FileNotFoundError
        If meta_path does not exist.
    OSError
        If the file exists but cannot be read.
    KeyError
        If group_col is not found among the metadata columns.
    ValueError
        If no samples remain after filtering to group_order.
    """
    meta_path = Path(meta_path)

    if not meta_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {meta_path}")

    try:
        df = pd.read_csv(meta_path, sep="\t", comment="#")
    except Exception as e:
        raise OSError(
            f"Could not read metadata file {meta_path}: {e}"
        ) from e

    if group_col not in df.columns:
        available = ", ".join(df.columns.tolist())
        raise KeyError(
            f"Column '{group_col}' not found in {meta_path}.\n"
            f"Available columns: {available}"
        )

    id_col = df.columns[0]  # QIIME2 convention: first column is sample-id

    # Extract TV-ID using a regex that matches the TV prefix followed by digits.
    # This is more robust than splitting on '-' because some sample IDs may
    # have different delimiter conventions across metadata files.
    df["TV"] = df[id_col].str.extract(r"(TV\d+)", expand=False)

    # Report how many rows were dropped due to missing TV or group values.
    n_before = len(df)
    df = df.dropna(subset=["TV", group_col])
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        print(
            f"  NOTE: Dropped {n_dropped} metadata rows with missing "
            f"TV-ID or '{group_col}' value.",
            file=sys.stderr,
        )

    # Filter to requested groups. Groups outside group_order are excluded
    # silently -- e.g. Marine, Unknown_Other, or unrelated markers.
    df_filtered = df[df[group_col].isin(group_order)]

    if df_filtered.empty:
        raise ValueError(
            f"No samples found in '{group_col}' matching groups: {group_order}.\n"
            f"Values present in this column: "
            f"{df[group_col].dropna().unique().tolist()}"
        )

    # Keep first occurrence of each TV-ID to avoid duplicates from
    # multi-marker metadata files.
    mapping = (
        df_filtered
        .drop_duplicates(subset="TV")
        .set_index("TV")[group_col]
        .to_dict()
    )

    return mapping


# ── Data assembly ──────────────────────────────────────────────────────────────

def build_sample_table(
    cutadapt_counts: dict,
    group_mapping: dict,
    group_order: list,
    threshold: int,
) -> pd.DataFrame:
    """
    Join cutadapt read counts with group metadata and apply detection threshold.

    Detection is defined as >= threshold primer-confirmed reads. This threshold
    represents a minimum signal level above which we consider adenovirus DPOL
    amplification to have occurred in the sample. The choice of threshold is
    scientifically important for adenovirus (unlike herpesvirus) because read
    counts range from 0 to 10,508 with no clear bimodal gap.

    Samples in cutadapt output but absent from metadata are reported as
    warnings. These are typically Marine birds or new birds not yet in the
    metadata file. Samples in metadata but absent from cutadapt output
    (e.g. TV240083, which had 0 input reads) are also reported.

    Parameters
    ----------
    cutadapt_counts : dict
        Output of parse_cutadapt_summary().
    group_mapping : dict
        Output of load_group_metadata().
    group_order : list of str
        Groups to include, in display order left to right.
    threshold : int
        Minimum primer-confirmed reads to call a sample positive.

    Returns
    -------
    pd.DataFrame with columns:
        tv (str), group (str), total (int), kept (int),
        pct (float), detected (bool)
    Sorted by group order then by kept reads descending within each group,
    so the highest-read positive samples appear on the left.

    Raises
    ------
    ValueError
        If no samples remain after joining cutadapt output with metadata.
    """
    rows = []
    unmatched_in_cutadapt = []  # cutadapt has them, metadata does not
    unmatched_in_metadata = []  # metadata has them, cutadapt does not

    for tv, vals in cutadapt_counts.items():
        grp = group_mapping.get(tv)
        if grp is None:
            # This TV-ID is not in the metadata for the requested groups.
            # Could be a Marine bird, a new sample, or a metadata gap.
            unmatched_in_cutadapt.append(tv)
            continue
        if grp not in group_order:
            # Sample is in metadata but belongs to an excluded group
            # (e.g. Unknown_Other excluded from COD analysis). Not a warning.
            continue
        rows.append({
            "tv": tv,
            "group": grp,
            "total": vals["total"],
            "kept": vals["kept"],
            "pct": vals["pct"],
            "detected": vals["kept"] >= threshold,
        })

    # Check for metadata samples with no cutadapt output at all
    for tv, grp in group_mapping.items():
        if tv not in cutadapt_counts and grp in group_order:
            unmatched_in_metadata.append(tv)

    if unmatched_in_cutadapt:
        print(
            f"  NOTE: {len(unmatched_in_cutadapt)} sample(s) in cutadapt "
            f"output have no metadata match and are excluded: "
            f"{', '.join(sorted(unmatched_in_cutadapt))}",
            file=sys.stderr,
        )

    if unmatched_in_metadata:
        print(
            f"  WARNING: {len(unmatched_in_metadata)} sample(s) in metadata "
            f"have no cutadapt output (may have failed sequencing or were "
            f"not in this run): {', '.join(sorted(unmatched_in_metadata))}",
            file=sys.stderr,
        )

    if not rows:
        raise ValueError(
            "No samples remained after joining cutadapt output with metadata. "
            "Check that TV-IDs match between cutadapt filenames and metadata."
        )

    df = pd.DataFrame(rows)

    # Sort: group order first (left-to-right), then by kept reads descending
    # so highest-signal positive samples appear at the left of each group.
    group_rank = {g: i for i, g in enumerate(group_order)}
    df["_group_rank"] = df["group"].map(group_rank)
    df = (
        df.sort_values(["_group_rank", "kept"], ascending=[True, False])
        .drop(columns="_group_rank")
        .reset_index(drop=True)
    )

    return df


# ── Fisher's exact test helper ─────────────────────────────────────────────────

def run_fisher(df: pd.DataFrame, group_a: str, group_b: str) -> tuple:
    """
    Run Fisher's exact test comparing detection rates between two groups.

    Uses a one-sided 2x2 contingency table:
        [[pos_a, neg_a],
         [pos_b, neg_b]]

    This is the appropriate test when sample sizes are small and we are
    comparing proportions (positive/negative counts). The two-sided p-value
    is reported because we do not have a strong prior on the direction of
    the effect for adenovirus.

    Parameters
    ----------
    df : pd.DataFrame
        Sample table with 'group' and 'detected' columns.
    group_a : str
        First group label.
    group_b : str
        Second group label.

    Returns
    -------
    tuple of (odds_ratio: float, p_value: float)

    Raises
    ------
    ValueError
        If either group has no samples in df.
    """
    a = df[df["group"] == group_a]
    b = df[df["group"] == group_b]

    if a.empty:
        raise ValueError(
            f"Group '{group_a}' has no samples in the data table. "
            f"Cannot run Fisher's exact test."
        )
    if b.empty:
        raise ValueError(
            f"Group '{group_b}' has no samples in the data table. "
            f"Cannot run Fisher's exact test."
        )

    pos_a = int(a["detected"].sum())
    neg_a = len(a) - pos_a
    pos_b = int(b["detected"].sum())
    neg_b = len(b) - pos_b

    odds, p = fisher_exact([[pos_a, neg_a], [pos_b, neg_b]])
    return odds, p


# ── Shared plot utilities ──────────────────────────────────────────────────────

def _add_group_dividers_and_labels(
    ax: plt.Axes,
    df: pd.DataFrame,
    group_order: list,
    group_x_positions: dict,
    color_map: dict,
    y_label_pos: float,
    label_template: str = "{group}\n{n_pos}/{n_tot} positive",
) -> None:
    """
    Add group labels above bars and dashed vertical dividers between groups.

    Parameters
    ----------
    ax : plt.Axes
        The axes to annotate.
    df : pd.DataFrame
        Sample table.
    group_order : list of str
        Groups in display order.
    group_x_positions : dict
        Maps group label -> list of x-positions for bars in that group.
    color_map : dict
        Maps group label -> hex color string.
    y_label_pos : float
        Y-coordinate for group label text (in data coordinates).
    label_template : str
        Format string for group labels. Available variables:
        {group}, {n_pos}, {n_tot}.
    """
    for i, grp in enumerate(group_order):
        positions = group_x_positions[grp]
        if not positions:
            continue

        midpoint = float(np.mean(positions))
        n_pos = int(df[df["group"] == grp]["detected"].sum())
        n_tot = len(positions)

        label = label_template.format(group=grp, n_pos=n_pos, n_tot=n_tot)
        ax.text(
            midpoint, y_label_pos,
            label,
            ha="center", va="bottom",
            fontsize=10, fontweight="bold",
            color=color_map[grp],
        )

        # Dashed divider after every group except the last
        if i < len(group_order) - 1:
            divider_x = max(positions) + 0.4
            ax.axvline(
                divider_x,
                color="#cccccc",
                linestyle="--",
                linewidth=1,
                zorder=1,
            )


def _save_figure(fig: plt.Figure, out_png: Path, out_svg: Path) -> None:
    """
    Save a figure to PNG (300 dpi) and SVG, then close it.

    Wraps both save calls in a try/finally so the figure is always closed
    even if saving fails, preventing matplotlib memory leaks.

    Parameters
    ----------
    fig : plt.Figure
        The figure to save.
    out_png : Path
        Output path for the PNG file.
    out_svg : Path
        Output path for the SVG file.

    Raises
    ------
    OSError
        If either file cannot be written.
    """
    try:
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        fig.savefig(out_svg, bbox_inches="tight")
        print(f"  \u2713 Saved: {out_png}")
        print(f"  \u2713 Saved: {out_svg}")
    except OSError as e:
        raise OSError(f"Could not save figure: {e}") from e
    finally:
        plt.close(fig)


# ── Presence / absence figure ──────────────────────────────────────────────────

def plot_presence_absence(
    df: pd.DataFrame,
    group_order: list,
    color_map: dict,
    threshold: int,
    fisher_pairs: list,
    primer_label: str,
    outdir: Path,
    stem: str,
) -> None:
    """
    Presence/absence barplot for adenovirus detection.

    Each bar represents one bird. Filled bars = detected (>= threshold
    primer-confirmed reads). Hatched bars = not detected. Bars are sorted
    within each group by read count descending so the strongest positives
    appear on the left.

    For each pair in fisher_pairs, a Fisher's exact test p-value is
    annotated. This supports reporting the DvT overall comparison AND the
    COD subgroup comparison (e.g. PI vs Trauma) in the same figure.

    The threshold is annotated on the figure because the adenovirus result
    is threshold-sensitive (unlike herpesvirus). This makes the detection
    criterion explicit to reviewers without requiring them to find it in
    the methods section.

    Parameters
    ----------
    df : pd.DataFrame
        Output of build_sample_table().
    group_order : list of str
        Groups in display order.
    color_map : dict
        Maps group label -> hex color string.
    threshold : int
        Reads threshold used for detection (for annotation only).
    fisher_pairs : list of (str, str)
        Pairs of group labels for Fisher's exact test annotation.
        e.g. [('Diseased', 'Trauma')] or [('Parasitic_Infectious', 'Trauma')]
    primer_label : str
        Short primer description for figure annotation.
        e.g. 'polF/polR (Wellehan et al. 2004)'
    outdir : Path
        Output directory.
    stem : str
        Output filename stem.
    """
    n_samples = len(df)
    fig_width = max(9, n_samples * 0.55 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    x = 0
    x_positions = []
    x_labels = []
    group_x_positions = {g: [] for g in group_order}

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
                # Not detected: white fill with hatching and colored border.
                # Hatching ensures the not-detected state is distinguishable
                # in greyscale printing and for color-blind readers.
                ax.bar(
                    x, 1,
                    color="white",
                    edgecolor=color_map[grp],
                    linewidth=1.5,
                    hatch="///",
                    zorder=2,
                )

            x_positions.append(x)
            x_labels.append(row["tv"].replace("TV", ""))
            group_x_positions[grp].append(x)
            x += 1

        x += 0.8  # visual gap between groups

    _add_group_dividers_and_labels(
        ax, df, group_order, group_x_positions, color_map,
        y_label_pos=1.06,
    )

    # Annotate Fisher's exact test results for each requested pair.
    # Stack multiple annotations vertically if more than one pair is tested.
    annotation_y = 0.04
    for pair in fisher_pairs:
        if len(pair) != 2:
            print(
                f"  WARNING: fisher_pairs entry {pair} does not have exactly "
                f"2 elements. Skipping this comparison.",
                file=sys.stderr,
            )
            continue
        group_a, group_b = pair
        try:
            _, p_val = run_fisher(df, group_a, group_b)
        except ValueError as e:
            print(f"  WARNING: Could not run Fisher's test: {e}", file=sys.stderr)
            continue

        star = " \u2605" if p_val <= 0.05 else " ns"
        color = SIG_COLOR if p_val <= 0.05 else NS_COLOR
        label = f"{group_a} vs {group_b}:  Fisher\u2019s exact  p={p_val:.3f}{star}"

        ax.annotate(
            label,
            xy=(0.98, annotation_y),
            xycoords="axes fraction",
            ha="right", va="bottom",
            fontsize=9, style="italic",
            color=color,
        )
        annotation_y += 0.10  # move next annotation up

    # Threshold annotation in bottom-left. Critical for this figure because
    # the adenovirus result changes with threshold choice.
    ax.annotate(
        f"Threshold: \u2265{threshold:,} primer-confirmed reads "
        f"({primer_label}) | result is threshold-sensitive",
        xy=(0.02, 0.04),
        xycoords="axes fraction",
        ha="left", va="bottom",
        fontsize=7.5, style="italic",
        color="#666666",
    )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=90, fontsize=7)
    ax.set_ylim(0, 1.35)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Not detected", "Detected"], fontsize=10)
    ax.set_ylabel(
        "Aviadenovirus DPOL detection\n(cutadapt primer-confirmed reads)",
        fontsize=10,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    _save_figure(
        fig,
        outdir / f"{stem}_presence.png",
        outdir / f"{stem}_presence.svg",
    )


# ── Relative abundance figure ──────────────────────────────────────────────────

def plot_relative_abundance(
    df: pd.DataFrame,
    group_order: list,
    color_map: dict,
    threshold: int,
    primer_label: str,
    outdir: Path,
    stem: str,
) -> None:
    """
    Relative abundance barplot for adenovirus detection.

    Each bar represents one bird. The colored portion = primer-confirmed reads
    (both polF and polR primers found by cutadapt). The grey portion = all
    other reads that did not survive --discard-untrimmed filtering.

    Unlike the herpesvirus relabund figure (where all bars are nearly full),
    the adenovirus figure shows genuine variation -- most birds have very
    few or no adenovirus reads, while a subset have substantial signal.
    This figure conveys the MAGNITUDE of detection, not just binary
    presence/absence. Both figures together tell the complete story.

    The threshold line is drawn as a horizontal reference at the proportion
    corresponding to the detection threshold, calculated using the median
    total read count across all samples. This gives the reader a visual
    sense of where the detection cutoff falls in the read distribution.

    Parameters
    ----------
    df : pd.DataFrame
        Output of build_sample_table().
    group_order : list of str
        Groups in display order.
    color_map : dict
        Maps group label -> hex color string.
    threshold : int
        Detection threshold in absolute read counts (for reference line).
    primer_label : str
        Short primer description for figure annotation.
    outdir : Path
        Output directory.
    stem : str
        Output filename stem.
    """
    n_samples = len(df)
    fig_width = max(9, n_samples * 0.55 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    x = 0
    x_positions = []
    x_labels = []
    group_x_positions = {g: [] for g in group_order}

    for grp in group_order:
        grp_df = df[df["group"] == grp]

        for _, row in grp_df.iterrows():
            prop_adeno = (
                row["kept"] / row["total"] if row["total"] > 0 else 0.0
            )
            prop_other = 1.0 - prop_adeno

            # Adenovirus primer-confirmed reads (bottom of stack)
            ax.bar(
                x, prop_adeno,
                color=color_map[grp],
                edgecolor="white",
                linewidth=0.3,
                zorder=2,
            )
            # Non-primer reads (top of stack, grey)
            ax.bar(
                x, prop_other,
                bottom=prop_adeno,
                color=NON_PRIMER_COLOR,
                edgecolor="white",
                linewidth=0.3,
                zorder=2,
            )

            x_positions.append(x)
            x_labels.append(row["tv"].replace("TV", ""))
            group_x_positions[grp].append(x)
            x += 1

        x += 0.8

    _add_group_dividers_and_labels(
        ax, df, group_order, group_x_positions, color_map,
        y_label_pos=1.02,
        label_template="{group} (n={n_tot})\n{n_pos} detected",
    )

    # Draw a horizontal reference line at the threshold proportion.
    # We use the median total reads across samples to convert the absolute
    # threshold to a proportion. This is approximate -- individual samples
    # with very different total read counts will have different absolute
    # thresholds -- but gives the reader a useful visual reference.
    median_total = df["total"].median()
    if median_total > 0:
        threshold_prop = threshold / median_total
        ax.axhline(
            threshold_prop,
            color="#444444",
            linestyle=":",
            linewidth=1,
            zorder=3,
            label=f"Detection threshold (~{threshold:,} reads at median depth)",
        )
        ax.legend(loc="upper right", fontsize=8, framealpha=0.85,
                  edgecolor="#cccccc")

    # Legend distinguishing adenovirus reads from non-primer reads
    legend_handles = [
        mpatches.Patch(
            color=color_map[g],
            label=f"{g} — aviadenovirus DPOL reads",
        )
        for g in group_order
    ] + [
        mpatches.Patch(
            color=NON_PRIMER_COLOR,
            label="Non-primer reads (not adenovirus-attributed)",
        )
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize=8,
        framealpha=0.85,
        edgecolor="#cccccc",
    )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=90, fontsize=7)
    ax.set_ylim(0, 1.2)
    ax.set_ylabel(
        "Proportion of reads\n(primer-confirmed / total input reads)",
        fontsize=10,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.annotate(
        f"Primers: {primer_label} | cutadapt --discard-untrimmed",
        xy=(0.02, 0.04),
        xycoords="axes fraction",
        ha="left", va="bottom",
        fontsize=7.5, style="italic",
        color="#666666",
    )

    plt.tight_layout()

    _save_figure(
        fig,
        outdir / f"{stem}_relabund.png",
        outdir / f"{stem}_relabund.svg",
    )


# ── Entry point ────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate aviadenovirus detection figures from cutadapt-confirmed "
            "polF/polR primer reads. Produces DvT and COD subgroup figures. "
            "See script header for full scientific rationale."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--cutadapt", required=True,
        help="Path to cutadapt_summary.txt from the adenovirus cutadapt run",
    )
    parser.add_argument(
        "--metadata-dvt", required=True,
        help="Metadata TSV for DvT comparison (needs 'Group' column)",
    )
    parser.add_argument(
        "--metadata-cod", required=True,
        help=(
            "Metadata TSV for COD subgroup comparison "
            "(needs COD_broad or equivalent column)"
        ),
    )
    parser.add_argument(
        "--group-by-dvt", default="Group",
        help="Column for DvT grouping (default: Group)",
    )
    parser.add_argument(
        "--group-by-cod", default="COD_broad",
        help="Column for COD grouping (default: COD_broad)",
    )
    parser.add_argument(
        "--dvt-order", nargs="+", default=["Diseased", "Trauma"],
        help="DvT group order (default: Diseased Trauma)",
    )
    parser.add_argument(
        "--cod-order", nargs="+",
        default=["Lead", "Parasitic_Infectious", "Trauma"],
        help="COD group order (default: Lead Parasitic_Infectious Trauma)",
    )
    parser.add_argument(
        "--threshold", type=int, default=500,
        help=(
            "Minimum primer-confirmed reads to call positive (default: 500). "
            "Key result: PI vs Trauma p=0.033 at >=500; p=0.125 at >=1000. "
            "Result is threshold-sensitive -- run at both thresholds."
        ),
    )
    parser.add_argument(
        "--palette", choices=["wong"], default="wong",
        help="Color palette (default: wong)",
    )
    parser.add_argument(
        "--outdir", default="results/adenovirus/figures/",
        help="Output directory (default: results/adenovirus/figures/)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(
            f"ERROR: Could not create output directory {outdir}: {e}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Stem includes threshold so files from different threshold runs
    # are clearly distinguished and do not overwrite each other.
    dvt_stem = f"adeno_cutadapt_DvT_wong_{args.threshold}reads"
    cod_stem = f"adeno_cutadapt_COD_wong_{args.threshold}reads"
    primer_label = "polF/polR (Wellehan et al. 2004)"

    # Color maps: DvT follows project convention (Diseased=blue, Trauma=orange)
    # COD follows convention (Lead=blue, PI=orange, Trauma=green)
    dvt_colors = {
        g: WONG_COLORS[i % len(WONG_COLORS)]
        for i, g in enumerate(args.dvt_order)
    }
    cod_colors = {
        g: WONG_COLORS[i % len(WONG_COLORS)]
        for i, g in enumerate(args.cod_order)
    }

    # ── Step 1: Parse cutadapt output ────────────────────────────────────────
    print(f"\nStep 1: Parsing cutadapt summary")
    print(f"  File: {args.cutadapt}")
    try:
        counts = parse_cutadapt_summary(Path(args.cutadapt))
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    n_samples = len(counts)
    kept_vals = [v["kept"] for v in counts.values()]
    n_zero = sum(1 for v in kept_vals if v == 0)
    print(f"  Parsed {n_samples} samples")
    print(f"  Primer-confirmed reads range: {min(kept_vals):,} -- {max(kept_vals):,}")
    print(f"  Median: {int(np.median(kept_vals)):,}")
    print(f"  Samples with 0 primer-confirmed reads: {n_zero}")
    print(f"  Samples positive at >= {args.threshold:,} reads: "
          f"{sum(1 for v in kept_vals if v >= args.threshold)}/{n_samples}")

    # ── Step 2: Load metadata (DvT) ──────────────────────────────────────────
    print(f"\nStep 2a: Loading DvT metadata")
    print(f"  File: {args.metadata_dvt}")
    try:
        dvt_groups = load_group_metadata(
            Path(args.metadata_dvt), args.group_by_dvt, args.dvt_order
        )
    except (FileNotFoundError, KeyError, ValueError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    print(f"  Found {len(dvt_groups)} samples in groups: {args.dvt_order}")

    # ── Step 3: Load metadata (COD) ──────────────────────────────────────────
    print(f"\nStep 2b: Loading COD metadata")
    print(f"  File: {args.metadata_cod}")
    try:
        cod_groups = load_group_metadata(
            Path(args.metadata_cod), args.group_by_cod, args.cod_order
        )
    except (FileNotFoundError, KeyError, ValueError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    print(f"  Found {len(cod_groups)} samples in groups: {args.cod_order}")

    # ── Step 4: Build sample tables ───────────────────────────────────────────
    print(f"\nStep 3: Building sample tables (threshold >= {args.threshold:,} reads)")

    print(f"  DvT:")
    try:
        dvt_df = build_sample_table(
            counts, dvt_groups, args.dvt_order, args.threshold
        )
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    for grp in args.dvt_order:
        grp_df = dvt_df[dvt_df["group"] == grp]
        n_pos = grp_df["detected"].sum()
        print(f"    {grp}: {n_pos}/{len(grp_df)} positive")

    print(f"  COD:")
    try:
        cod_df = build_sample_table(
            counts, cod_groups, args.cod_order, args.threshold
        )
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    for grp in args.cod_order:
        grp_df = cod_df[cod_df["group"] == grp]
        n_pos = grp_df["detected"].sum()
        print(f"    {grp}: {n_pos}/{len(grp_df)} positive")

    # ── Step 5: Generate figures ──────────────────────────────────────────────
    print(f"\nStep 4: Generating figures")
    print(f"  Output: {outdir}")

    # DvT presence/absence
    # Fisher's pair: Diseased vs Trauma (the primary DvT comparison)
    print(f"\n  DvT presence/absence:")
    try:
        plot_presence_absence(
            dvt_df, args.dvt_order, dvt_colors,
            args.threshold,
            fisher_pairs=[("Diseased", "Trauma")],
            primer_label=primer_label,
            outdir=outdir,
            stem=dvt_stem,
        )
    except OSError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # DvT relative abundance
    print(f"\n  DvT relative abundance:")
    try:
        plot_relative_abundance(
            dvt_df, args.dvt_order, dvt_colors,
            args.threshold, primer_label,
            outdir=outdir, stem=dvt_stem,
        )
    except OSError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # COD presence/absence
    # Fisher's pairs: PI vs Trauma is the key finding (p=0.033 at >=500 reads)
    # Lead vs Trauma is included for completeness
    # We do NOT apply FDR correction here because this is exploratory and
    # we are reporting all pairwise comparisons transparently.
    print(f"\n  COD presence/absence:")
    try:
        plot_presence_absence(
            cod_df, args.cod_order, cod_colors,
            args.threshold,
            fisher_pairs=[
                ("Parasitic_Infectious", "Trauma"),
                ("Lead", "Trauma"),
            ],
            primer_label=primer_label,
            outdir=outdir,
            stem=cod_stem,
        )
    except OSError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # COD relative abundance
    print(f"\n  COD relative abundance:")
    try:
        plot_relative_abundance(
            cod_df, args.cod_order, cod_colors,
            args.threshold, primer_label,
            outdir=outdir, stem=cod_stem,
        )
    except OSError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\nDone. Four figures written to {outdir}")


if __name__ == "__main__":
    main()
