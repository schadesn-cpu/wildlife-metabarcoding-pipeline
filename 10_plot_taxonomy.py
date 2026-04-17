#!/usr/bin/env python3
"""
10_plot_taxonomy.py
===================
Generate publication-quality stacked taxonomy barplots from the relative
abundance TSV produced by 07_taxonomy_table.py.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PIPELINE POSITION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  07_taxonomy_table.py  →  taxonomy_relabund_L{N}_{marker}.tsv
                                          ↓
                          10_plot_taxonomy.py   ← metadata TSV
                                          ↓
                          barplot_{marker}_{group}_{palette}.png/.svg

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WHAT THIS SCRIPT PRODUCES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  - Stacked barplot of relative abundance (%) per sample
  - Top-N taxa shown individually, remainder collapsed to "Other"
  - Samples ordered within groups by dominant taxon (default sort mode)
    OR
  - Samples clustered by their own dominant taxon across all samples,
    for unstructured baseline figures (--sort-mode dominant-taxon)
  - Group labels above bars with n= counts (group mode only)
  - Vertical dashed dividers between groups (group mode only)
  - Optional thin dividers and small italic labels between dominant-taxon
    clusters (dominant-taxon mode with --cluster-labels)
  - Italic taxon names in legend (non-genus labels like "Other" left upright)
  - Both PNG (300 dpi, slides/SharePoint) and SVG (Illustrator-editable)

Relative abundance is calculated among classified reads only (output of 07_).
Methods note: "Relative abundance was calculated among reads assigned at
[level] level following [database] classification."

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SORT MODES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  --sort-mode group           (default, existing behavior)
    Samples grouped by --group-by column, ordered left-to-right by
    --group-order, sorted within each group by the overall-dominant taxon.
    Produces grouped bar charts with colored headers and dashed dividers.

  --sort-mode dominant-taxon  (new, for baseline figures)
    All samples placed in a single waterfall, clustered by each sample's
    own dominant taxon. Clusters ordered left-to-right by the number of
    samples with that taxon as dominant. Within each cluster, samples
    ordered by that taxon's relative abundance descending. No group
    dividers or headers — use for unstructured baseline community views.
    When this mode is active, --group-by is optional; if provided, it is
    ignored for layout but still checked for sample matching.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PALETTES — kept in sync with 09_plot_diversity.py and 09c_visualize_diversity.py
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

  # 16S Diseased vs Trauma, purple palette (group mode, default)
  python 10_plot_taxonomy.py \\
      --relabund  results/16S/DvT/taxonomy/taxonomy_relabund_L6_16S.tsv \\
      --metadata  metadata/qiime/metadata_16S.tsv \\
      --group-by  Group \\
      --marker    16S \\
      --palette   purple \\
      --outdir    results/16S/DvT/figures/taxonomy/

  # Unstructured baseline across all samples (NEW — for slide 8 / Figure 2)
  python 10_plot_taxonomy.py \\
      --relabund     results/16S/all/taxonomy/taxonomy_relabund_L6_16S.tsv \\
      --metadata     metadata/qiime/metadata_16S_ecoseason.tsv \\
      --sort-mode    dominant-taxon \\
      --cluster-labels \\
      --marker       16S \\
      --palette      wong \\
      --output-stem  16S_baseline_dominant_taxon_wong \\
      --outdir       results/16S/all/figures/taxonomy/

  # Seasonal grouping, explicit order, top 20 taxa (group mode)
  python 10_plot_taxonomy.py \\
      --relabund    results/16S/DvT/taxonomy/taxonomy_relabund_L6_16S.tsv \\
      --metadata    metadata/qiime/metadata_16S.tsv \\
      --group-by    Season \\
      --group-order Spring Summer Fall Winter \\
      --top-n       20 \\
      --marker      16S \\
      --palette     purple \\
      --outdir      results/16S/DvT/figures/taxonomy/

  # MiFish species level
  python 10_plot_taxonomy.py \\
      --relabund  results/MiFish/all/taxonomy/taxonomy_relabund_L7_MiFish.tsv \\
      --metadata  metadata/qiime/metadata_MiFish.tsv \\
      --group-by  Group \\
      --marker    MiFish \\
      --palette   redblue \\
      --outdir    results/MiFish/all/figures/taxonomy/

  # List metadata columns available for --group-by
  python 10_plot_taxonomy.py \\
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
from collections import Counter
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
FONT_SIZE_CLUSTER = 9

# ---------------------------------------------------------------------------
# Palettes
# ---------------------------------------------------------------------------
# Group colors: used for group header labels above bars.
# Kept in sync with 09_plot_diversity.py and 09c_visualize_diversity.py.
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
        # 25-taxon colorblind-safe palette:
        # Colors 1-7:   Wong 2011 (no black — replaced with Tol indigo in pos 8)
        # Colors 8-15:  Paul Tol 'Muted' qualitative palette (doi:10.5281/zenodo.3381072)
        # Colors 16-25: Additional Tol Muted + Tol Bright colors for larger top_n
        # Positions 1-15 are byte-identical to the prior palette — existing
        # wong plots with top-n <= 15 are unchanged.
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
            "#332288",  # Tol Muted indigo  — taxon 16  (deep blue-purple)
            "#117733",  # Tol Muted forest  — taxon 17
            "#882255",  # Tol Muted wine    — taxon 18
            "#DDCC77",  # Tol Muted sand    — taxon 19
            "#EE6677",  # Tol Bright red    — taxon 20
            "#228833",  # Tol Bright green  — taxon 21
            "#4477AA",  # Tol Bright blue   — taxon 22
            "#CCBB44",  # Tol Bright yellow — taxon 23
            "#AA3377",  # Tol Bright magenta — taxon 24
            "#66CCEE",  # Tol Bright cyan   — taxon 25
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


def load_all_samples_from_metadata(metadata_path: Path) -> List[str]:
    """
    Load the sample-id index from a metadata TSV without requiring a specific
    grouping column. Used by dominant-taxon sort mode when --group-by is omitted.
    """
    with metadata_path.open(encoding="utf-8-sig") as f:
        lines = f.readlines()
    header = lines[0]
    data_lines = [l for l in lines[1:] if not l.startswith("#")]
    content = header + "".join(data_lines)
    sep = "\t" if "\t" in header else ","
    df = pd.read_csv(io.StringIO(content), sep=sep, index_col=0, dtype=str)
    df.index = df.index.str.strip()
    return df.index.tolist()


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
    Load a taxonomy_relabund_*.tsv produced by 07_taxonomy_table.py.
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

    # Clean taxonomy labels for display
    label_map = {t: _clean_taxon_label(t) for t in top_taxa}
    result.index = [label_map.get(i, i) for i in result.index]
    display_taxa = [label_map.get(t, t) for t in top_taxa] + ["Other"]
    return result, display_taxa

# ---------------------------------------------------------------------------
# Barplot
# ---------------------------------------------------------------------------

def _clean_taxon_label(name: str) -> str:
    """
    Convert a full QIIME2 taxonomy string to a short readable label.

    e.g. 'k__Metazoa;p__Chordata;c__Actinopteri;o__Clupeiformes;f__Clupeidae;g__Alosa;Unclassified'
         → 'Alosa Unclassified'
         'k__Metazoa;p__Chordata;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified'
         → 'uncl. Chordata'
    """
    if ";" not in name:
        return name  # Already a short name (e.g. 'Other', 'Trematoda')

    parts = [p.strip() for p in name.split(";")]
    # Strip rank prefixes: k__, p__, c__, o__, f__, g__, s__
    cleaned = [re.sub(r"^[kpcofgs]__", "", p) for p in parts]

    # Separate classified from unclassified levels
    UNCLASSIFIED = {"Unclassified", "uncultured", "uncl.", "", "X"}
    classified = [(i, v) for i, v in enumerate(cleaned)
                  if v not in UNCLASSIFIED and not v.startswith("uncl")]
    unclassified = [(i, v) for i, v in enumerate(cleaned)
                    if v in UNCLASSIFIED or v.startswith("uncl")]

    if not classified:
        return "uncl. " + cleaned[0] if cleaned else name

    last_idx, last_val = classified[-1]

    # If the level immediately after is unclassified, note it
    if last_idx + 1 < len(cleaned) and cleaned[last_idx + 1] in UNCLASSIFIED:
        return f"{last_val} X" if last_val.endswith("ae") or last_val.endswith("a") \
               else f"{last_val} Unclassified"

    return last_val


def _is_plain_name(name: str) -> bool:
    if name in ("Other",): return True
    if name.startswith("uncl."): return True
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


def _sort_by_dominant_taxon(
    samples: List[str],
    relabund_pct: pd.DataFrame,
    display_taxa: List[str],
) -> Tuple[List[str], List[Tuple[str, int, int]]]:
    """
    Sort samples by their own dominant taxon, producing a waterfall-style
    layout for unstructured baseline figures.

    For each sample, identify the single most abundant named taxon (from
    display_taxa excluding 'Other'). Cluster samples by that dominant taxon.
    Order clusters left-to-right by the number of samples with that taxon as
    dominant (most common first). Within each cluster, order samples by that
    taxon's relative abundance descending.

    Samples whose largest contribution is 'Other' (or whose named-taxon
    abundances all sum to zero) fall into a trailing 'Other' cluster.

    Returns
    -------
    ordered_samples : list of sample IDs in the new order
    cluster_spans   : list of (dominant_taxon, start_idx, end_idx) tuples
                      suitable for annotating cluster boundaries
    """
    non_other = [t for t in display_taxa if t != "Other"]
    if not non_other or not samples:
        return list(samples), []

    # Identify each sample's dominant named taxon. If a sample has no named-taxon
    # signal (all zero), assign it to the 'Other' cluster at the end.
    sample_dominant: Dict[str, str] = {}
    for s in samples:
        col = relabund_pct[s].reindex(non_other).fillna(0.0)
        if col.sum() == 0:
            sample_dominant[s] = "Other"
        else:
            sample_dominant[s] = col.idxmax()

    # Count samples per dominant taxon, then order clusters.
    # Primary sort: descending sample count.
    # Tie-break: overall taxon rank in display_taxa (most abundant overall first).
    counts = Counter(sample_dominant.values())
    named_clusters = sorted(
        [t for t in counts if t != "Other"],
        key=lambda t: (-counts[t], non_other.index(t) if t in non_other else 999),
    )
    cluster_order = named_clusters + (["Other"] if "Other" in counts else [])

    ordered_samples: List[str] = []
    cluster_spans: List[Tuple[str, int, int]] = []

    for taxon in cluster_order:
        members = [s for s in samples if sample_dominant[s] == taxon]
        if not members:
            continue
        if taxon in relabund_pct.index:
            order = relabund_pct.loc[taxon, members].sort_values(ascending=False)
            members = list(order.index)
        start = len(ordered_samples)
        ordered_samples.extend(members)
        end = len(ordered_samples) - 1
        cluster_spans.append((taxon, start, end))

    return ordered_samples, cluster_spans


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
    sort_mode: str = "group",
    cluster_labels: bool = False,
    min_cluster_label_size: int = 2,
) -> None:
    """
    Generate and save a stacked barplot (PNG + SVG).

    Parameters
    ----------
    relabund       : taxa × samples DataFrame (values 0–1, from 07_)
    sample_groups  : {sample_id: group_label} mapping. In dominant-taxon mode
                     this can contain Nones; samples with None are still
                     plotted (they pass sample-matching at the caller level).
    group_order    : ordered list of group labels to plot (group mode only).
                     In dominant-taxon mode this is ignored for layout.
    top_n          : number of top taxa to show before collapsing to Other
    palette_name   : 'purple', 'redblue', or 'wong'
    marker         : marker label for figure title (e.g. '16S')
    group_column   : metadata column name used for grouping (for display)
    title          : override figure title (None = auto-generate)
    outpath_stem   : output path without extension (will save .png and .svg)
    show_title     : whether to draw the figure title
    sort_mode      : 'group' (default, grouped bar chart with headers) or
                     'dominant-taxon' (unstructured waterfall baseline)
    cluster_labels : when sort_mode='dominant-taxon', draw small italic labels
                     above each dominant-taxon cluster
    min_cluster_label_size : when cluster_labels is True, only label clusters
                     with at least this many samples. Single-sample clusters
                     remain unlabeled to prevent overlap. All taxa are still
                     visible in the legend. Default: 2.
    """
    palette = PALETTES[palette_name]
    group_colors = palette["group_colors"]
    taxa_colors  = palette["taxa_colors"]

    # ── Build ordered sample list ─────────────────────────────────────────
    ordered_samples: List[str] = []
    group_spans: List[Tuple[str, int, int, str]] = []  # (label, start, end, color)
    cluster_spans: List[Tuple[str, int, int]] = []

    if sort_mode == "dominant-taxon":
        # Unstructured baseline: pool all samples, cluster by dominant taxon.
        all_samples = [s for s in relabund.columns]
        if not all_samples:
            log.error("No samples in relabund table. Cannot plot.")
            sys.exit(1)

        # Rank top-N taxa across the pooled sample set, then use that ranking
        # to determine each sample's dominant taxon.
        rel_pct_all, display_taxa = collapse_to_top_n(relabund, top_n)
        ordered_samples, cluster_spans = _sort_by_dominant_taxon(
            all_samples, rel_pct_all, display_taxa,
        )
    else:
        # Existing grouped behavior.
        for gi, grp in enumerate(group_order):
            grp_samples = [s for s in relabund.columns
                           if sample_groups.get(s) == grp]
            if not grp_samples:
                log.warning(
                    "Group '%s' has no samples in the relabund table — skipping.", grp,
                )
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
            log.error(
                "No samples matched any group in %s. Check --group-by and metadata.",
                group_order,
            )
            sys.exit(1)

    # ── Final top-N collapse across the ordered sample set ────────────────
    rel_pct, display_taxa = collapse_to_top_n(relabund, top_n, ordered_samples)

    # Assign taxa colors. For top_n small enough to fit the hardcoded palette
    # (≤15 named taxa), behavior is unchanged — existing figures stay
    # byte-identical. For larger top_n, sample a matplotlib colormap so every
    # named taxon gets a distinguishable color instead of falling back to grey.
    # 'Other' always keeps the palette's designated grey.
    _other_color = taxa_colors[-1]
    _named_taxa = [t for t in display_taxa if t != "Other"]
    _n_named = len(_named_taxa)
    _n_palette_named = len(taxa_colors) - 1  # last slot reserved for Other

    if _n_named <= _n_palette_named:
        _named_colors = list(taxa_colors[:_n_named])
    else:
        _cmap_name = {
            "purple":  "Purples",
            "redblue": "RdBu_r",
            "wong":    "tab20",
        }.get(palette_name, "viridis")
        _cmap = plt.get_cmap(_cmap_name)
        if _cmap_name == "Purples":
            # Skip the palest end so light bars remain visible on white.
            _xs = np.linspace(0.35, 0.95, _n_named)
        elif _cmap_name == "tab20":
            # Qualitative cycle through tab20's 20 distinct colors.
            _xs = [(i % 20) / 20 + 0.025 for i in range(_n_named)]
        else:
            _xs = np.linspace(0, 1, _n_named)
        _named_colors = [_cmap(x) for x in _xs]

    tax_color_map = {t: c for t, c in zip(_named_taxa, _named_colors)}
    if "Other" in display_taxa:
        tax_color_map["Other"] = _other_color
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

    # ── Group dividers and header labels (group mode only) ────────────────
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

    # ── Cluster dividers and labels (dominant-taxon mode only) ────────────
    if sort_mode == "dominant-taxon" and cluster_spans:
        n_suppressed = 0
        for i, (taxon, start, end) in enumerate(cluster_spans):
            if i > 0:
                ax.axvline(
                    start - 0.5,
                    color="#cccccc", lw=0.8, linestyle=":", zorder=5,
                )
            if cluster_labels:
                n_cl = end - start + 1
                if n_cl < min_cluster_label_size:
                    # Suppress label for small clusters to prevent overlap.
                    # Taxa remain visible in the legend.
                    n_suppressed += 1
                    continue
                mid = (start + end) / 2
                label_text = f"{taxon}  (n={n_cl})"
                txt = ax.text(
                    mid, 101.5,
                    label_text,
                    ha="center", va="bottom",
                    fontsize=FONT_SIZE_CLUSTER,
                    color="#444444",
                )
                # Italicize the taxon portion if it's a genus-style name.
                # Use a simple approach: italicize the whole label unless it's
                # a plain name like 'Other' or 'uncl. Chordata'.
                if not _is_plain_name(taxon):
                    txt.set_fontstyle("italic")
        if cluster_labels and n_suppressed > 0:
            log.info(
                "Suppressed %d cluster label(s) below --min-cluster-label-size=%d "
                "(taxa still visible in legend).",
                n_suppressed, min_cluster_label_size,
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
        if title:
            fig_title = title
        elif sort_mode == "dominant-taxon":
            fig_title = (
                f"Genus-Level Relative Abundance — {marker}\n"
                f"All samples, clustered by dominant taxon"
            )
        else:
            fig_title = (
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
Build and return the argument parser for 10_plot_taxonomy.py.

    Key arguments: --relabund (TSV from 07_taxonomy_table.py), --marker,
    --group-by (metadata column), --palette, --top-n, --outdir, --sort-mode.
    A --list-columns flag is also available to inspect available metadata
    columns without producing any figures.
    """
    p = argparse.ArgumentParser(
        prog="10_plot_taxonomy.py",
        description=(
            "Generate stacked taxonomy barplots from 07_taxonomy_table.py output.\n"
            "Produces PNG (300 dpi) and SVG per group × palette combination."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    req = p.add_argument_group("required arguments")
    req.add_argument(
        "--relabund", required=True, type=Path,
        help=(
            "Relative abundance TSV produced by 07_taxonomy_table.py "
            "(taxonomy_relabund_L{N}_{marker}.tsv). "
            "Rows = taxa, columns = samples, values = 0–1."
        ),
    )
    req.add_argument(
        "--metadata", required=True, type=Path,
        help=(
            "QIIME 2 metadata TSV or source metadata CSV/TSV. "
            "Must contain a sample-id column (or TV column) and the --group-by column "
            "(unless --sort-mode dominant-taxon, which does not require --group-by)."
        ),
    )
    req.add_argument(
        "--group-by", dest="group_by", default=None,
        help=(
            "Metadata column to group samples by (e.g. Group, Season). "
            "Required for --sort-mode group (default). "
            "Optional for --sort-mode dominant-taxon; if provided in that mode, "
            "it is ignored for layout."
        ),
    )

    opt = p.add_argument_group("plot options")
    opt.add_argument(
        "--group-order", dest="group_order", nargs="*", default=None,
        help=(
            "Explicit group order left-to-right (e.g. --group-order Diseased Trauma). "
            "Default: alphabetical. Ignored when --sort-mode dominant-taxon."
        ),
    )
    opt.add_argument(
        "--sort-mode",
        choices=["group", "dominant-taxon"],
        default="group",
        dest="sort_mode",
        help=(
            "How to order samples on the x-axis. "
            "'group' (default) = group samples by --group-by, sort within each "
            "group by the overall-dominant taxon. Preserves existing behavior. "
            "'dominant-taxon' = unstructured baseline view. Samples clustered by "
            "their own dominant taxon, clusters ordered by prevalence, samples "
            "within each cluster ordered by that taxon's abundance descending. "
            "Use for baseline figures without group faceting."
        ),
    )
    opt.add_argument(
        "--cluster-labels", action="store_true", dest="cluster_labels",
        help=(
            "When --sort-mode dominant-taxon, add small italic labels above "
            "each dominant-taxon cluster showing taxon name and sample count. "
            "Default: off (clean waterfall with no annotations)."
        ),
    )
    opt.add_argument(
        "--min-cluster-label-size", type=int, default=2,
        dest="min_cluster_label_size",
        help=(
            "When --cluster-labels is set, only label clusters with at least "
            "this many samples. Single-sample clusters remain unlabeled to "
            "prevent overlap. All taxa remain visible in the legend regardless. "
            "Set to 1 to label every cluster. Default: 2."
        ),
    )
    opt.add_argument(
        "--include-controls", action="store_true", dest="include_controls",
        help=(
            "When --sort-mode dominant-taxon, include process controls "
            "(NTC-, PAC-, XB- prefixes) in the plot. Default: excluded. "
            "Use for QC visualization of contamination signal; omit for "
            "clean baseline figures of biological samples only."
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

    Reads the relative abundance TSV produced by 07_taxonomy_table.py, groups
    samples by the specified metadata column (group mode) or clusters them by
    their own dominant taxon (dominant-taxon mode), collapses low-abundance
    taxa into 'Other', and writes one PNG (300 dpi) + SVG pair per palette to
    --outdir. Returns 0 on success, 2 if required input files are missing,
    1 on error.
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
    if args.sort_mode == "group":
        if not args.group_by:
            log.error(
                "--group-by is required in --sort-mode group (default). "
                "Use --list-columns to see available columns, or pass "
                "--sort-mode dominant-taxon for an unstructured baseline plot."
            )
            return 2
    else:
        # dominant-taxon mode
        if args.group_by:
            log.warning(
                "--group-by '%s' is ignored in --sort-mode dominant-taxon "
                "(baseline layout does not use group faceting).",
                args.group_by,
            )
        if args.group_order:
            log.warning(
                "--group-order is ignored in --sort-mode dominant-taxon.",
            )
        if args.cluster_labels:
            log.info("Cluster labels enabled above each dominant-taxon cluster.")

    if not args.relabund.exists():
        log.error("Relabund file not found: %s", args.relabund)
        return 2
    if not args.metadata.exists():
        log.error("Metadata file not found: %s", args.metadata)
        return 2

    # ── Load data ─────────────────────────────────────────────────────────
    relabund = load_relabund(args.relabund)

    # ── Match samples to groups (or build a passthrough map) ──────────────
    if args.sort_mode == "group":
        group_series = load_metadata(args.metadata, args.group_by)
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
        log.info("Matched %d / %d samples to metadata",
                 len(matched), len(relabund.columns))

        # Determine group order
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
    else:
        # dominant-taxon mode: filter controls, then build passthrough sample map.
        #
        # Controls (NTC-, PAC-, XB-) must be removed explicitly here because
        # dominant-taxon mode does not use --group-by, so the usual group-based
        # exclusion does not apply. Pass --include-controls to override (e.g.
        # for QC visualization of contamination signal in process controls).
        if not args.include_controls:
            CONTROL_PREFIXES = ("NTC-", "PAC-", "XB-")
            biological_samples = [
                s for s in relabund.columns
                if not any(s.startswith(p) for p in CONTROL_PREFIXES)
            ]
            n_controls = len(relabund.columns) - len(biological_samples)
            if n_controls > 0:
                excluded = [s for s in relabund.columns if s not in biological_samples]
                log.info(
                    "Excluded %d control sample(s) from baseline plot: %s",
                    n_controls, excluded,
                )
                log.info("Pass --include-controls to keep controls in the plot.")
                relabund = relabund[biological_samples]
        else:
            log.info("--include-controls set: controls will be plotted alongside biological samples.")

        # metadata is used only for sanity-checking sample IDs in this mode
        try:
            meta_samples = load_all_samples_from_metadata(args.metadata)
            matched_count = sum(
                1 for s in relabund.columns
                if s in meta_samples
                or any(s.startswith(m) for m in meta_samples)
                or any(m.startswith(s) for m in meta_samples)
            )
            log.info(
                "Metadata sanity check: %d / %d relabund samples have a "
                "matching entry in metadata",
                matched_count, len(relabund.columns),
            )
        except Exception as e:
            log.warning("Could not sanity-check sample IDs against metadata: %s", e)

        sample_groups = {s: None for s in relabund.columns}
        group_order = []
        log.info(
            "Sort mode: dominant-taxon (%d samples will be plotted as one pool)",
            len(relabund.columns),
        )

    # ── Build output path stem ────────────────────────────────────────────
    outdir = args.outdir or args.relabund.parent
    if args.output_stem:
        stem_name = args.output_stem
    else:
        marker_str   = f"_{args.marker}"  if args.marker   else ""
        if args.sort_mode == "dominant-taxon":
            group_str = "_baseline_dominant_taxon"
        else:
            group_str = f"_{args.group_by}"
        suffix_str   = f"_{args.suffix}"  if args.suffix   else ""
        palette_str  = f"_{args.palette}"
        stem_name    = f"barplot{marker_str}{group_str}{suffix_str}{palette_str}"
    outpath_stem = outdir / stem_name

    log.info("Output stem : %s", outpath_stem)
    log.info("Top-N taxa  : %d", args.top_n)
    log.info("Palette     : %s", args.palette)
    log.info("Sort mode   : %s", args.sort_mode)

    if args.dry_run:
        log.info("DRY RUN — no files written.")
        for ext in (".png", ".svg"):
            log.info("  Would write: %s", outpath_stem.with_suffix(ext))
        return 0

    # ── Plot ──────────────────────────────────────────────────────────────
    plot_barplot(
        relabund                = relabund,
        sample_groups           = sample_groups,
        group_order             = group_order,
        top_n                   = args.top_n,
        palette_name            = args.palette,
        marker                  = args.marker,
        group_column            = args.group_by or "",
        title                   = args.title,
        outpath_stem            = outpath_stem,
        show_title              = not args.no_title,
        sort_mode               = args.sort_mode,
        cluster_labels          = args.cluster_labels,
        min_cluster_label_size  = args.min_cluster_label_size,
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
