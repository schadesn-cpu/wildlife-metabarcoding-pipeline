#!/usr/bin/env python3
"""
09_plot_diversity.py
====================
Standalone diversity visualization script for metabarcoding projects.
Reads QIIME2 artifacts directly (.qza files) — no QIIME2 installation needed.

Generates:
  - 2D PCoA plots (beta diversity) — one per distance matrix artifact
  - Alpha diversity strip/box plots — one per alpha vector artifact

Output formats: PNG (300 dpi) + SVG (Illustrator-editable) per figure.

Usage:
  # Single beta metric, colored by Group
  python 09_plot_diversity.py pcoa \\
    --artifact weighted_unifrac_pcoa_results.qza \\
    --metadata metadata/qiime/metadata_16S.tsv \\
    --color-by Group \\
    --output-dir results/figures/

  # All beta metrics at once
  python 09_plot_diversity.py pcoa \\
    --artifact weighted_unifrac_pcoa_results.qza unweighted_unifrac_pcoa_results.qza \\
                bray_curtis_pcoa_results.qza jaccard_pcoa_results.qza \\
    --metadata metadata/qiime/metadata_16S.tsv \\
    --color-by Group \\
    --output-dir results/figures/

  # Alpha diversity
  python 09_plot_diversity.py alpha \\
    --artifact faith_pd_vector.qza shannon_vector.qza \\
               observed_features_vector.qza evenness_vector.qza \\
    --metadata metadata/qiime/metadata_16S.tsv \\
    --group-by Group \\
    --output-dir results/figures/

  # Panel of all four beta metrics in one figure
  python 09_plot_diversity.py pcoa \\
    --artifact weighted_unifrac_pcoa_results.qza unweighted_unifrac_pcoa_results.qza \\
                bray_curtis_pcoa_results.qza jaccard_pcoa_results.qza \\
    --metadata metadata/qiime/metadata_16S.tsv \\
    --color-by Group --panel \\
    --output-dir results/figures/
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import re
import sys
import traceback
import warnings
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Suppress matplotlib layout warnings that fire spuriously on tight_layout.
# Scoped to matplotlib only — other UserWarnings (e.g. from pandas or scipy)
# remain visible so genuine issues are not silently swallowed.
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Style constants — edit these to match your publication style
# ---------------------------------------------------------------------------

FIGURE_DPI      = 300
POINT_SIZE      = 80
POINT_ALPHA     = 0.85
FONT_FAMILY     = "Arial"
FONT_SIZE_TITLE = 13
FONT_SIZE_AXIS  = 11
FONT_SIZE_TICK  = 9
FONT_SIZE_LEGEND = 9
FONT_SIZE_ANNOT  = 8

# Transparency for the filled interior of confidence ellipses.
# The ellipse edge is drawn at full opacity; the fill uses this alpha.
ELLIPSE_FILL_ALPHA = 0.12
ELLIPSE_EDGE_ALPHA = 0.70

GRID_ALPHA = 0.25

# Colorblind-friendly palette (Wong 2011, doi:10.1038/nmeth.1618)
# Use this palette when reviewers request a citable colorblind-safe scheme.
PALETTE_DEFAULT = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # green
    "#CC79A7",  # pink
    "#56B4E9",  # sky blue
    "#D55E00",  # vermillion
    "#F0E442",  # yellow
    "#000000",  # black
]

MARKER_CYCLE = ["o", "s", "^", "D", "v", "P", "X", "*"]

# ---------------------------------------------------------------------------
# Named palettes — edit here to change figures globally.
# Use --palette purple / --palette redblue / --palette wong on the CLI,
# or pass raw hex codes: --palette "#B22222,#2E86C1"
# ---------------------------------------------------------------------------
PALETTES = {
    "purple": [
        "#7B2D8B",  # dark purple  — Group 1 / Diseased
        "#C19FD8",  # lavender     — Group 2 / Trauma
        "#4B1369",  # deep purple  — Group 3
        "#D09EE0",  # light purple — Group 4
        "#2D0A40", "#E0BAEC", "#9870B0", "#F0D6F5",
    ],
    "redblue": [
        "#B22222",  # red      — recommended for multi-marker or colorblind
        "#2E86C1",  # blue
        "#E74C3C",  # light red    — Group 3
        "#7FB3D3",  # light blue   — Group 4
        "#7B241C", "#1A5276", "#F1948A", "#2980B9",
    ],
    "wong": [
        "#0072B2",  # blue     — Wong 2011 colorblind-safe (citable)
        "#E69F00",  # orange
        "#009E73",  # green
        "#CC79A7",  # pink
        "#56B4E9",  # sky blue
        "#D55E00",  # vermillion
        "#F0E442",  # yellow
        "#000000",  # black
    ],
}

# ---------------------------------------------------------------------------
# Module-level metric name lookup — used in cmd_pcoa to find PERMANOVA QZVs.
# Maps pretty display names back to the snake_case used in QIIME2 filenames.
# ---------------------------------------------------------------------------
_PRETTY_TO_SNAKE: Dict[str, str] = {
    "Weighted UniFrac":   "weighted_unifrac",
    "Unweighted UniFrac": "unweighted_unifrac",
    "Bray-Curtis":        "bray_curtis",
    "Jaccard":            "jaccard",
}


def _resolve_palette(palette_arg: Optional[str]) -> List[str]:
    """
    Resolve --palette argument to a list of hex color strings.

    Accepts either a named palette key ('purple', 'redblue', 'wong')
    or a comma-separated list of raw hex codes ('#B22222,#2E86C1').
    Falls back to PALETTE_DEFAULT (Wong 2011) if None.
    """
    if palette_arg is None:
        return PALETTE_DEFAULT
    if palette_arg in PALETTES:
        return PALETTES[palette_arg]
    # Raw hex string: "#B22222,#2E86C1,..."
    return [c.strip() for c in palette_arg.split(",") if c.strip()]


# ===========================================================================
# QIIME2 artifact reading (no QIIME2 needed — artifacts are just zip files)
# ===========================================================================

def _find_data_file(zf: zipfile.ZipFile, filename: str) -> Optional[str]:
    """Find a file by name anywhere inside the QZA zip archive."""
    for name in zf.namelist():
        if name.endswith(f"/data/{filename}") or name.endswith(f"/{filename}"):
            return name
    return None


def read_pcoa_artifact(path: Path) -> Tuple[pd.DataFrame, List[float], str]:
    """
    Read a QIIME2 PCoA results artifact (.qza).

    Returns
    -------
    coords_df   : DataFrame with sample IDs as index, PC1/PC2/... columns
    prop_expl   : List of proportion of variance explained per axis (0–1 scale)
    metric_name : Human-readable metric name guessed from the filename
    """
    if not path.exists():
        raise FileNotFoundError(f"Artifact not found: {path}")
    if path.suffix != ".qza":
        raise ValueError(f"Expected a .qza file, got: {path}")

    try:
        with zipfile.ZipFile(path, "r") as zf:
            ord_path = _find_data_file(zf, "ordination.txt")
            if not ord_path:
                raise ValueError(
                    f"Could not find ordination.txt inside {path.name}\n"
                    f"  Is this a PCoA results artifact? Check it was generated by "
                    f"qiime diversity pcoa or core-metrics-phylogenetic."
                )
            with zf.open(ord_path) as f:
                content = f.read().decode("utf-8")
    except zipfile.BadZipFile as exc:
        raise ValueError(
            f"Could not open {path.name} as a zip file: {exc}\n"
            f"  The file may be incomplete or corrupted."
        ) from exc

    coords, prop_expl = _parse_ordination(content)
    metric_name = _guess_metric_name(path.name)
    return coords, prop_expl, metric_name


def _parse_ordination(content: str) -> Tuple[pd.DataFrame, List[float]]:
    """
    Parse QIIME2 ordination.txt format into a coordinates DataFrame.

    The ordination.txt file is divided into blank-line-separated sections.
    We extract the 'Proportion explained' and 'Site' sections, which contain
    the per-axis variance fractions and the per-sample PC coordinates.
    """
    sections = content.strip().split("\n\n")
    prop_expl: List[float] = []
    coords: Dict[str, List[float]] = {}
    n_axes = 0

    for section in sections:
        lines = section.strip().split("\n")
        if not lines:
            continue
        header = lines[0]

        if header.startswith("Proportion explained"):
            try:
                n_axes = int(header.split("\t")[1])
                prop_expl = [float(x) for x in lines[1].split("\t")]
            except (IndexError, ValueError) as exc:
                raise ValueError(
                    f"Could not parse 'Proportion explained' section: {exc}"
                ) from exc

        elif header.startswith("Site"):
            try:
                parts = header.split("\t")
                n_samples = int(parts[1])
                n_axes    = int(parts[2])
            except (IndexError, ValueError) as exc:
                raise ValueError(
                    f"Could not parse 'Site' header: {exc}"
                ) from exc

            for line in lines[1:]:
                if not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                sample_id = parts[0]
                try:
                    values = [float(x) for x in parts[1:n_axes + 1]]
                    coords[sample_id] = values
                except ValueError as exc:
                    raise ValueError(
                        f"Could not parse coordinates for sample '{sample_id}': {exc}"
                    ) from exc

    if not coords:
        raise ValueError("No sample coordinates found in ordination.txt")
    if not prop_expl:
        raise ValueError("No 'Proportion explained' section found in ordination.txt")

    n_pc = min(len(v) for v in coords.values())
    col_names = [f"PC{i+1}" for i in range(n_pc)]
    df = pd.DataFrame.from_dict(coords, orient="index", columns=col_names)
    df.index.name = "SampleID"
    return df, prop_expl[:n_pc]


def read_alpha_artifact(path: Path) -> Tuple[pd.Series, str]:
    """
    Read a QIIME2 alpha diversity vector artifact (.qza).

    Returns
    -------
    series      : Series with sample IDs as index and diversity values
    metric_name : Human-readable metric name from the TSV column header
    """
    if not path.exists():
        raise FileNotFoundError(f"Artifact not found: {path}")
    if path.suffix != ".qza":
        raise ValueError(f"Expected a .qza file, got: {path}")

    try:
        with zipfile.ZipFile(path, "r") as zf:
            alpha_path = _find_data_file(zf, "alpha-diversity.tsv")
            if not alpha_path:
                raise ValueError(
                    f"Could not find alpha-diversity.tsv inside {path.name}\n"
                    f"  Is this an alpha diversity vector artifact?"
                )
            with zf.open(alpha_path) as f:
                content = f.read().decode("utf-8")
    except zipfile.BadZipFile as exc:
        raise ValueError(
            f"Could not open {path.name} as a zip file: {exc}\n"
            f"  The file may be incomplete or corrupted."
        ) from exc

    try:
        df = pd.read_csv(io.StringIO(content), sep="\t", index_col=0)
    except Exception as exc:
        raise ValueError(
            f"Could not parse alpha-diversity.tsv from {path.name}: {exc}"
        ) from exc

    if df.empty or df.shape[1] == 0:
        raise ValueError(f"Alpha diversity TSV appears empty in {path.name}")

    metric_col  = df.columns[0]
    series      = df[metric_col].copy()
    series.index.name = "SampleID"
    metric_name = _pretty_metric_name(metric_col)
    return series, metric_name


def read_permanova_qzv(qzv_path: Path) -> Optional[Dict]:
    """
    Parse a QIIME2 beta-group-significance QZV (PERMANOVA) for figure annotation.

    Returns a dict with keys 'pseudo_f' and 'p_value', or None if parsing fails.

    QIIME2 renders PERMANOVA results as an HTML table. The relevant fields are:
      - test_statistic  → pseudo-F statistic (Anderson 2001)
      - p-value         → permutation p-value (default 999 permutations)
      - sample_size     → n (total samples tested — NOT the F-statistic)

    IMPORTANT — PERMANOVA assumption:
      PERMANOVA tests whether group centroids differ in multivariate space but
      assumes equal within-group dispersion (homogeneity of dispersions).
      If PERMDISP is also significant, interpret PERMANOVA results cautiously —
      the signal may reflect dispersion differences rather than centroid shifts.
      Always report PERMDISP alongside PERMANOVA in the methods.

    Parsing is best-effort: QIIME2's HTML schema has changed across versions.
    On failure, logs a warning and returns None so the figure is still produced
    without stats rather than crashing the entire run.
    """
    try:
        with zipfile.ZipFile(qzv_path) as zf:
            html_files = [n for n in zf.namelist() if n.endswith("index.html")]
            if not html_files:
                log.warning(
                    "No index.html found inside %s — cannot parse PERMANOVA stats.",
                    qzv_path.name,
                )
                return None
            with zf.open(html_files[0]) as f:
                html = f.read().decode("utf-8", errors="replace")

        # Parse all <th>...</th> <td>...</td> pairs from the overview table.
        # QIIME2 uses this structure for its stats summary section.
        rows = re.findall(
            r'<th[^>]*>(.*?)</th>\s*<td[^>]*>(.*?)</td>',
            html, re.DOTALL
        )
        data = {}
        for k, v in rows:
            key = re.sub(r'<[^>]+>', '', k).strip().lower().replace(' ', '_')
            val = re.sub(r'<[^>]+>', '', v).strip()
            data[key] = val

        pseudo_f = data.get("test_statistic", "")
        p_value  = data.get("p-value", data.get("p_value", ""))

        if pseudo_f and p_value:
            return {
                "pseudo_f": float(pseudo_f),
                "p_value":  float(p_value),
            }

        # QZV parsed but expected keys not found — QIIME2 may have changed
        # its HTML schema. Log what keys were present for debugging.
        log.warning(
            "Could not find test_statistic/p-value in %s. "
            "Keys found: %s. QIIME2 HTML schema may have changed.",
            qzv_path.name, list(data.keys()),
        )
        return None

    except zipfile.BadZipFile as exc:
        # File exists on disk but is not a valid ZIP/QZV — likely a partial write.
        log.warning(
            "Could not open %s as a ZIP archive: %s", qzv_path.name, exc
        )
        return None

    except Exception as exc:
        # Catch-all for unexpected parsing failures (malformed HTML, encoding
        # errors, etc.). Log the full error so it can be diagnosed without
        # silently producing a figure with missing statistics.
        log.warning(
            "Unexpected error parsing PERMANOVA stats from %s: %s",
            qzv_path.name, exc,
        )
        return None


# ===========================================================================
# Metadata loading
# ===========================================================================

def load_metadata(path: Path) -> pd.DataFrame:
    """
    Load a sample metadata file (CSV or TSV).

    Handles:
      - BOM characters introduced by Excel exports (utf-8-sig encoding)
      - QIIME2 TSV format with '#SampleID' as first column
      - '#q2:types' directive rows (skipped)
      - Leading/trailing whitespace in column names and values
    """
    if not path.exists():
        raise FileNotFoundError(f"Metadata file not found: {path}")

    suffix = path.suffix.lower()
    sep = "\t" if suffix in (".tsv", ".txt") else ","

    try:
        df = pd.read_csv(path, sep=sep, dtype=str, encoding="utf-8-sig")
    except Exception as exc:
        raise ValueError(f"Could not read metadata file {path}: {exc}") from exc

    # Drop QIIME2 '#q2:types' directive rows if present
    df = df[~df.iloc[:, 0].str.startswith("#q2:types", na=False)].reset_index(drop=True)

    # Normalize column names: strip whitespace, remove BOM, remove leading '#'
    df.columns = [c.strip().lstrip("\ufeff").lstrip("#") for c in df.columns]

    # Find the sample ID column — try common conventions in priority order
    id_col = None
    for candidate in ["SampleID", "sample-id", "sample_id", "TV", "id"]:
        if candidate in df.columns:
            id_col = candidate
            break
    if id_col is None:
        id_col = df.columns[0]

    df = df.set_index(id_col)
    df.index = df.index.str.strip()
    df.index.name = "SampleID"

    # Strip whitespace from all string columns
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.strip()

    return df


def match_sample_ids(
    artifact_ids: pd.Index,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """
    Match QIIME2 artifact sample IDs to metadata rows.

    Handles the common loon project pattern where artifact IDs look like
    'TV230007-GI-16S_S1483_L002' but metadata uses short IDs like 'TV230007'.

    Matching strategy (in priority order):
      1. Direct match (exact string equality)
      2. Prefix match — artifact ID starts with the metadata ID
      3. Substring match — metadata ID is contained anywhere in artifact ID

    Unmatched samples are reported as a warning and excluded from figures,
    rather than crashing, to allow partial runs on incomplete metadata.
    """
    # Fast path: all IDs match directly
    direct = artifact_ids.intersection(metadata.index)
    if len(direct) == len(artifact_ids):
        return metadata.loc[artifact_ids]

    mapping: Dict[str, Optional[str]] = {}
    meta_ids = metadata.index.tolist()

    for art_id in artifact_ids:
        if art_id in metadata.index:
            mapping[art_id] = art_id
            continue
        # Prefix: 'TV230007-GI-16S_S1483' starts with 'TV230007'
        matched = [m for m in meta_ids if art_id.startswith(m)]
        if matched:
            mapping[art_id] = max(matched, key=len)  # longest match wins
            continue
        # Substring: metadata ID contained anywhere in artifact ID
        matched = [m for m in meta_ids if m in art_id]
        if matched:
            mapping[art_id] = max(matched, key=len)
            continue
        mapping[art_id] = None  # no match found

    unmatched = [k for k, v in mapping.items() if v is None]
    if unmatched:
        log.warning(
            "%d sample(s) in artifact not found in metadata (will be excluded):\n%s%s",
            len(unmatched),
            "\n".join(f"  {u}" for u in unmatched[:5]),
            f"\n  ... and {len(unmatched)-5} more" if len(unmatched) > 5 else "",
        )

    matched_art_ids  = [a for a in artifact_ids if mapping[a] is not None]
    matched_meta_ids = [mapping[a] for a in matched_art_ids]

    result = metadata.loc[matched_meta_ids].copy()
    result.index = matched_art_ids
    result.index.name = "SampleID"
    return result


# ===========================================================================
# Utility helpers
# ===========================================================================

def _guess_metric_name(filename: str) -> str:
    """Extract a human-readable metric name from a QIIME2 artifact filename."""
    name = filename.replace(".qza", "").replace(".qzv", "")
    for suffix in ["_pcoa_results", "_distance_matrix", "_emperor", "_vector"]:
        name = name.replace(suffix, "")
    return _pretty_metric_name(name)


def _pretty_metric_name(name: str) -> str:
    """Convert snake_case QIIME2 metric names to publication-ready labels."""
    mapping = {
        "weighted_unifrac":   "Weighted UniFrac",
        "unweighted_unifrac": "Unweighted UniFrac",
        "bray_curtis":        "Bray-Curtis",
        "jaccard":            "Jaccard",
        "faith_pd":           "Faith's PD",
        "faith_pd_vector":    "Faith's PD",
        "shannon":            "Shannon Diversity",
        "shannon_entropy":    "Shannon Diversity",
        "observed_features":  "Observed Features",
        "pielou_evenness":    "Pielou's Evenness",
        "evenness":           "Pielou's Evenness",
    }
    return mapping.get(name, name.replace("_", " ").title())


def _get_color_map(
    values: pd.Series,
    palette: Optional[List[str]] = None,
) -> Tuple[Dict[str, str], List[str]]:
    """
    Build a color mapping for unique values in a Series.

    Returns
    -------
    color_dict        : dict mapping each unique value → hex color string
    ordered_categories: list of unique values in sorted order
    """
    if palette is None:
        palette = PALETTE_DEFAULT

    unique_vals = sorted(values.dropna().unique(), key=str)

    if len(unique_vals) > len(palette):
        # Fall back to a matplotlib colormap for large numbers of categories
        cmap = plt.cm.get_cmap("tab20", len(unique_vals))
        palette = [matplotlib.colors.to_hex(cmap(i)) for i in range(len(unique_vals))]

    color_dict = {str(v): palette[i % len(palette)] for i, v in enumerate(unique_vals)}
    return color_dict, [str(v) for v in unique_vals]


def _confidence_ellipse(
    x: np.ndarray,
    y: np.ndarray,
    ax: plt.Axes,
    n_std: float = 2.448,
    facecolor: str = "none",
    **kwargs,
) -> None:
    """
    Draw a covariance-based confidence ellipse around a set of 2D points.

    The ellipse represents the region containing approximately 95% of the
    probability mass for a bivariate normal distribution. The correct
    multiplier for 95% in 2D is n_std = sqrt(chi2.ppf(0.95, df=2)) ≈ 2.448,
    NOT 2.0. Using n_std=2.0 would produce an ~86% ellipse — a common error.

    Reference: Friendly, Monette & Fox (2013), Statistical Science.
               scipy.stats.chi2.ppf(0.95, df=2) = 5.991 → sqrt = 2.448

    Parameters
    ----------
    n_std : float
        Number of standard deviations for the ellipse radius.
        Default 2.448 gives the true 95% confidence ellipse for bivariate normal data.
    """
    if len(x) < 3:
        return  # Cannot compute a meaningful covariance with fewer than 3 points

    cov = np.cov(x, y)
    if np.any(np.isnan(cov)) or np.linalg.det(cov) == 0:
        return  # Degenerate covariance (e.g. all points collinear) — skip

    pearson    = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)

    ellipse = mpatches.Ellipse(
        (0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs,
    )

    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_x, mean_y = np.mean(x), np.mean(y)

    transform = (
        matplotlib.transforms.Affine2D()
        .rotate_deg(45)
        .scale(scale_x, scale_y)
        .translate(mean_x, mean_y)
    )
    ellipse.set_transform(transform + ax.transData)
    ax.add_patch(ellipse)


def _save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    """
    Save a figure as both PNG (300 dpi, for slides/print) and SVG
    (vector, for Illustrator editing).

    Raises a descriptive OSError if the output directory is not writable,
    rather than letting matplotlib produce a cryptic permission error.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = [
        (out_dir / f"{stem}.png", "png"),
        (out_dir / f"{stem}.svg", "svg"),
    ]
    for path, fmt in outputs:
        try:
            fig.savefig(
                path,
                dpi=FIGURE_DPI if fmt == "png" else None,
                format=fmt,
                bbox_inches="tight",
                facecolor="white",
            )
            print(f"  ✓ Saved: {path.name}")
        except OSError as exc:
            raise OSError(
                f"Could not write figure to {path}.\n"
                f"  Check that {out_dir} exists and is writable.\n"
                f"  Original error: {exc}"
            ) from exc


def _p_to_stars(p: float) -> str:
    """Convert a p-value to APA-style significance stars."""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return ""


# ===========================================================================
# PCoA plotting
# ===========================================================================

def plot_pcoa(
    coords: pd.DataFrame,
    prop_expl: List[float],
    metadata_col: pd.Series,
    metric_name: str,
    color_col: str,
    pc_axes: Tuple[int, int] = (1, 2),
    ellipses: bool = True,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    palette: Optional[List[str]] = None,
    permanova_result: Optional[Dict] = None,
    show_title: bool = True,
    group_order: Optional[List[str]] = None,
) -> plt.Figure:
    """
    Draw a 2D PCoA scatter plot with optional 95% confidence ellipses.

    Parameters
    ----------
    coords           : DataFrame of PC coordinates (index = SampleID)
    prop_expl        : Proportion of variance explained per axis (0–1 scale)
    metadata_col     : Series mapping SampleID → group label
    metric_name      : e.g. 'Weighted UniFrac' — used for axis and title labels
    color_col        : Name of the metadata column being plotted (for legend title)
    pc_axes          : Which PC axes to plot, 1-indexed (default: PC1 vs PC2)
    ellipses         : If True, draw 95% confidence ellipses per group (n >= 3)
    ax               : Existing axes to draw into; if None, a new figure is created
    permanova_result : Optional dict {'pseudo_f': float, 'p_value': float} to annotate.
                       Note: PERMANOVA assumes homogeneous within-group dispersion;
                       always report PERMDISP alongside PERMANOVA in the manuscript.
    group_order      : Explicit ordering of groups for legend and color assignment.
                       Groups not present in data are silently skipped.
    show_title       : If False, suppress the subplot title (manuscript mode)
    """
    pc_x = f"PC{pc_axes[0]}"
    pc_y = f"PC{pc_axes[1]}"

    # Align metadata with ordination coordinates
    common = coords.index.intersection(metadata_col.index)
    if len(common) == 0:
        raise ValueError(
            f"No overlapping sample IDs between PCoA coords and metadata "
            f"column '{color_col}'.\n"
            f"  PCoA IDs (first 3): {list(coords.index[:3])}\n"
            f"  Metadata IDs (first 3): {list(metadata_col.index[:3])}"
        )
    coords_plot = coords.loc[common]
    meta_plot   = metadata_col.loc[common]

    # Drop samples with no group assignment (e.g. controls with blank Group)
    valid = meta_plot.notna() & (meta_plot.astype(str).str.strip() != "")
    coords_plot = coords_plot.loc[valid]
    meta_plot   = meta_plot.loc[valid].astype(str)

    color_dict, categories = _get_color_map(meta_plot, palette)

    # Apply group ordering if specified — controls legend order and color assignment
    if group_order:
        categories = [c for c in group_order if c in categories]

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.get_figure()

    # Plot each group's points and optional confidence ellipse
    for i, group in enumerate(categories):
        mask  = meta_plot == group
        if not mask.any():
            continue
        x      = coords_plot.loc[mask, pc_x].values
        y      = coords_plot.loc[mask, pc_y].values
        color  = color_dict[group]
        marker = MARKER_CYCLE[i % len(MARKER_CYCLE)]

        ax.scatter(
            x, y,
            c=color,
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            marker=marker,
            label=group,
            edgecolors="white",
            linewidths=0.5,
            zorder=3,
        )

        if ellipses and len(x) >= 3:
            # Filled ellipse (transparent)
            _confidence_ellipse(
                x, y, ax,
                facecolor=color,
                alpha=ELLIPSE_FILL_ALPHA,
                edgecolor="none",
                zorder=2,
            )
            # Ellipse border (more opaque)
            _confidence_ellipse(
                x, y, ax,
                facecolor="none",
                edgecolor=color,
                alpha=ELLIPSE_EDGE_ALPHA,
                linewidth=1.2,
                linestyle="--",
                zorder=2,
            )

    # Axis labels: include the proportion of variance explained per axis
    pct_x = prop_expl[pc_axes[0] - 1] * 100 if len(prop_expl) >= pc_axes[0] else 0
    pct_y = prop_expl[pc_axes[1] - 1] * 100 if len(prop_expl) >= pc_axes[1] else 0
    ax.set_xlabel(f"{pc_x} ({pct_x:.1f}%)", fontsize=FONT_SIZE_AXIS)
    ax.set_ylabel(f"{pc_y} ({pct_y:.1f}%)", fontsize=FONT_SIZE_AXIS)

    # Title — starred and colored purple if PERMANOVA significant (p ≤ 0.05)
    sig         = permanova_result and permanova_result.get("p_value", 1) <= 0.05
    plot_title  = title or metric_name
    title_color = "#6A0DAD" if sig else "black"
    if show_title:
        ax.set_title(
            f"{plot_title} *" if sig else plot_title,
            fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=8, color=title_color,
        )

    # PERMANOVA annotation in bottom-right corner
    if permanova_result:
        pf = permanova_result.get("pseudo_f")
        pv = permanova_result.get("p_value")
        if pf is not None and pv is not None:
            star         = " *" if pv <= 0.05 else ""
            annot        = f"PERMANOVA  F={pf:.3f}, p={pv:.3f}{star}"
            annot_color  = "#6A0DAD" if pv <= 0.05 else "#444444"
            ax.annotate(
                annot,
                xy=(0.98, 0.02),
                xycoords="axes fraction",
                fontsize=FONT_SIZE_ANNOT,
                color=annot_color,
                style="italic",
                ha="right",
                va="bottom",
            )

    # Crosshairs, grid, and spine cleanup
    ax.axhline(0, color="grey", linewidth=0.5, alpha=0.4, zorder=1)
    ax.axvline(0, color="grey", linewidth=0.5, alpha=0.4, zorder=1)
    ax.grid(True, alpha=GRID_ALPHA, linewidth=0.5)
    ax.tick_params(labelsize=FONT_SIZE_TICK)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    if standalone:
        ax.legend(
            title=color_col,
            fontsize=FONT_SIZE_LEGEND,
            title_fontsize=FONT_SIZE_LEGEND,
            framealpha=0.8,
            edgecolor="#cccccc",
        )
        fig.tight_layout()

    return fig


def cmd_pcoa(args: argparse.Namespace) -> None:
    """Generate 2D PCoA plots from one or more PCoA artifacts."""
    out_dir   = Path(args.output_dir)
    metadata  = load_metadata(Path(args.metadata))
    palette   = _resolve_palette(args.palette)
    stats_dir = Path(args.stats_dir) if args.stats_dir else None

    if args.color_by not in metadata.columns:
        raise ValueError(
            f"--color-by column '{args.color_by}' not found in metadata.\n"
            f"  Available columns: {', '.join(metadata.columns.tolist())}"
        )

    pc_axes = tuple(int(x) for x in args.pc_axes.split(","))
    if len(pc_axes) != 2:
        raise ValueError(
            f"--pc-axes must be two comma-separated integers, e.g. '1,2', got: {args.pc_axes}"
        )

    artifacts = [Path(a) for a in args.artifact]
    for a in artifacts:
        if not a.exists():
            raise FileNotFoundError(f"Artifact not found: {a}")

    # Load all artifacts first so we fail fast on missing files before plotting
    all_data = []
    for artifact_path in artifacts:
        print(f"\n  Reading: {artifact_path.name}")
        coords, prop_expl, metric_name = read_pcoa_artifact(artifact_path)
        matched_meta = match_sample_ids(coords.index, metadata)

        if args.color_by not in matched_meta.columns:
            raise ValueError(
                f"Column '{args.color_by}' not found after metadata matching."
            )
        color_col = matched_meta[args.color_by]

        # Load PERMANOVA result if a stats directory was provided.
        # Uses module-level _PRETTY_TO_SNAKE to convert display name → filename.
        permanova_result = None
        if stats_dir:
            snake_metric = _PRETTY_TO_SNAKE.get(
                metric_name,
                metric_name.lower().replace(" ", "_").replace("-", "_"),
            )
            group_col = args.color_by
            candidates = [
                stats_dir / f"{snake_metric}_permanova_{group_col}.qzv",
                stats_dir / f"{snake_metric}_permanova.qzv",
                stats_dir / f"{metric_name}_permanova_{group_col}.qzv",
                stats_dir / f"{metric_name}_permanova.qzv",
            ]
            for qzv in candidates:
                if qzv.exists():
                    permanova_result = read_permanova_qzv(qzv)
                    if permanova_result:
                        print(
                            f"  PERMANOVA loaded: F={permanova_result['pseudo_f']:.3f}, "
                            f"p={permanova_result['p_value']:.3f}"
                        )
                    break

        all_data.append((
            coords, prop_expl, metric_name, color_col, artifact_path.stem, permanova_result
        ))

    group_order = getattr(args, "group_order", None)

    if args.panel and len(all_data) > 1:
        n     = len(all_data)
        ncols = min(2, n)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
        axes_flat = np.array(axes).flatten()

        for i, (coords, prop_expl, metric_name, color_col, stem, permanova_result) in enumerate(all_data):
            ax    = axes_flat[i]
            label = "ABCDEFGHIJ"[i]
            # Panel letter in top-left corner with a white background box for legibility
            ax.text(
                0.02, 0.98, label, transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="top", ha="left", zorder=10,
                bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.85),
            )
            plot_pcoa(
                coords, prop_expl, color_col, metric_name,
                args.color_by, pc_axes=pc_axes,
                ellipses=not args.no_ellipses,
                palette=palette,
                permanova_result=permanova_result,
                show_title=not args.no_title,
                group_order=group_order,
                ax=ax,
            )
            # Show legend only on the last subplot to avoid redundancy
            if i == len(all_data) - 1:
                ax.legend(
                    title=args.color_by,
                    fontsize=FONT_SIZE_LEGEND,
                    title_fontsize=FONT_SIZE_LEGEND,
                    framealpha=0.8,
                    edgecolor="#cccccc",
                    bbox_to_anchor=(1.02, 1),
                    loc="upper left",
                )

        for j in range(len(all_data), len(axes_flat)):
            axes_flat[j].set_visible(False)

        if not args.no_title:
            fig.suptitle(
                f"Beta Diversity PCoA — colored by {args.color_by}",
                fontsize=FONT_SIZE_TITLE + 1, fontweight="bold", y=1.01,
            )
        fig.tight_layout()
        stem = args.output_stem or f"pcoa_panel_{args.color_by}"
        _save_figure(fig, out_dir, stem)
        plt.close(fig)

    else:
        for coords, prop_expl, metric_name, color_col, stem, permanova_result in all_data:
            fig = plot_pcoa(
                coords, prop_expl, color_col, metric_name,
                args.color_by, pc_axes=pc_axes,
                ellipses=not args.no_ellipses,
                palette=palette,
                permanova_result=permanova_result,
                show_title=not args.no_title,
                group_order=group_order,
            )
            out_stem = args.output_stem or f"pcoa_{stem}_{args.color_by}"
            _save_figure(fig, out_dir, out_stem)
            plt.close(fig)

    print(f"\n  All PCoA figures saved to: {out_dir.resolve()}")


# ===========================================================================
# Alpha diversity plotting
# ===========================================================================

def plot_alpha(
    series: pd.Series,
    metadata_col: pd.Series,
    metric_name: str,
    group_col: str,
    ax: Optional[plt.Axes] = None,
    palette: Optional[List[str]] = None,
    show_stats: bool = True,
    show_title: bool = True,
    group_order: Optional[List[str]] = None,
) -> plt.Figure:
    """
    Draw an alpha diversity strip + box plot with per-sample jittered points.

    Statistical annotation:
      Non-parametric tests are used throughout because alpha diversity metrics
      are not normally distributed, particularly at small sample sizes (n < 30).
      - 2 groups:  Mann-Whitney U (two-sided) — tests whether one group
                   tends to have higher diversity than the other.
      - 3+ groups: Kruskal-Wallis (omnibus) — tests whether at least one group
                   differs. No post-hoc correction is applied here; report as
                   exploratory if making pairwise comparisons after a KW test.

    Parameters
    ----------
    series       : Alpha diversity values indexed by SampleID
    metadata_col : Series mapping SampleID → group label
    metric_name  : Display name for the y-axis label and title
    group_col    : Metadata column name being plotted (for legend title)
    show_stats   : If True, annotate with test name and p-value
    show_title   : If False, suppress the subplot title (manuscript mode)
    group_order  : Explicit group ordering for x-axis. Groups absent from
                   data are silently skipped.
    """
    # Align diversity values with metadata
    common = series.index.intersection(metadata_col.index)
    if len(common) == 0:
        raise ValueError(
            f"No overlapping sample IDs between alpha diversity and "
            f"metadata column '{group_col}'."
        )
    series_plot = series.loc[common]
    meta_plot   = metadata_col.loc[common]

    # Drop samples with no group assignment (e.g. controls with blank Group)
    valid       = meta_plot.notna() & (meta_plot.astype(str).str.strip() != "")
    series_plot = series_plot.loc[valid]
    meta_plot   = meta_plot.loc[valid].astype(str)

    # Apply group order filter before building color map
    if group_order:
        mask        = meta_plot.isin(group_order)
        meta_plot   = meta_plot[mask]
        series_plot = series_plot[mask]

    color_dict, categories = _get_color_map(meta_plot, palette)

    # Enforce group ordering for x-axis position and color assignment
    if group_order:
        categories = [c for c in group_order if c in categories]

    group_data   = []
    group_labels = []
    for cat in categories:
        vals = series_plot[meta_plot == cat].dropna().values
        if len(vals) > 0:
            group_data.append(vals)
            group_labels.append(cat)

    if not group_data:
        raise ValueError(f"No data to plot for metric '{metric_name}'")

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(max(3, len(group_labels) * 1.6), 5))
    else:
        fig = ax.get_figure()

    positions = list(range(1, len(group_labels) + 1))

    # Box plot — fliers suppressed because individual points are shown below
    bp = ax.boxplot(
        group_data,
        positions=positions,
        widths=0.45,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=2),
        whiskerprops=dict(linewidth=1.2, color="#555555"),
        capprops=dict(linewidth=1.2, color="#555555"),
        boxprops=dict(linewidth=1.2),
        zorder=2,
    )
    for patch, cat in zip(bp["boxes"], group_labels):
        color = color_dict[cat]
        patch.set_facecolor(color)
        patch.set_alpha(0.35)
        patch.set_edgecolor(color)

    # Jittered strip plot — shows every individual sample for transparency
    rng = np.random.default_rng(seed=42)  # fixed seed for reproducibility
    for i, (vals, cat) in enumerate(zip(group_data, group_labels)):
        color  = color_dict[cat]
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(
            positions[i] + jitter,
            vals,
            c=color,
            s=POINT_SIZE * 0.7,
            alpha=0.8,
            zorder=4,
            edgecolors="white",
            linewidths=0.4,
        )

    # Statistical annotation
    # Non-parametric tests used — see docstring for justification.
    if show_stats and len(group_data) >= 2:
        if len(group_data) == 2:
            # Mann-Whitney U: tests whether one distribution is stochastically
            # greater than the other. Equivalent to Wilcoxon rank-sum test.
            _stat, p_val = stats.mannwhitneyu(
                group_data[0], group_data[1], alternative="two-sided"
            )
            stat_label = (
                f"Mann-Whitney U\np < 0.001" if p_val < 0.001
                else f"Mann-Whitney U\np = {p_val:.3f}"
            )
        else:
            # Kruskal-Wallis: omnibus test across all groups.
            # Significant result indicates at least one group differs;
            # does not identify which pairs differ (no post-hoc here).
            _stat, p_val = stats.kruskal(*group_data)
            stat_label = (
                f"Kruskal-Wallis\np < 0.001" if p_val < 0.001
                else f"Kruskal-Wallis\np = {p_val:.3f}"
            )

        ax.annotate(
            stat_label,
            xy=(0.02, 0.88),
            xycoords="axes fraction",
            fontsize=FONT_SIZE_ANNOT,
            color="#333333",
            style="italic",
            ha="left",
            va="top",
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="white",
                edgecolor="none",
                alpha=0.85,
            ),
        )

        # Significance bracket for 2-group comparison only
        if len(group_data) == 2:
            sig_str = _p_to_stars(p_val)
            if sig_str:
                y_max     = max(np.max(v) for v in group_data)
                y_range   = y_max - min(np.min(v) for v in group_data)
                y_bracket = y_max + y_range * 0.08
                ax.plot(
                    [1, 1, 2, 2],
                    [y_bracket, y_bracket + y_range * 0.03,
                     y_bracket + y_range * 0.03, y_bracket],
                    lw=1.2, color="#333333",
                )
                ax.text(
                    1.5, y_bracket + y_range * 0.04,
                    sig_str, ha="center", va="bottom",
                    fontsize=FONT_SIZE_AXIS, color="#333333",
                )

    # Axis formatting
    ax.set_xticks(positions)
    ax.set_xticklabels(group_labels, fontsize=FONT_SIZE_TICK)
    ax.set_ylabel(metric_name, fontsize=FONT_SIZE_AXIS)
    if show_title:
        ax.set_title(metric_name, fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=8)
    ax.tick_params(axis="y", labelsize=FONT_SIZE_TICK)
    ax.grid(axis="y", alpha=GRID_ALPHA, linewidth=0.5)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # n label below each group on the x-axis
    for i, (vals, cat) in enumerate(zip(group_data, group_labels)):
        ax.text(
            positions[i], ax.get_ylim()[0],
            f"n={len(vals)}",
            ha="center", va="top",
            fontsize=FONT_SIZE_ANNOT, color="#555555",
        )

    if standalone:
        handles = [
            mpatches.Patch(facecolor=color_dict[cat], label=cat)
            for cat in group_labels
        ]
        ax.legend(
            handles=handles,
            title=group_col,
            fontsize=FONT_SIZE_LEGEND,
            title_fontsize=FONT_SIZE_LEGEND,
            framealpha=0.8,
            edgecolor="#cccccc",
            loc="upper right",
        )
        fig.tight_layout()

    return fig


def cmd_alpha(args: argparse.Namespace) -> None:
    """Generate alpha diversity strip/box plots from one or more alpha vector artifacts."""
    out_dir  = Path(args.output_dir)
    metadata = load_metadata(Path(args.metadata))

    if args.group_by not in metadata.columns:
        raise ValueError(
            f"--group-by column '{args.group_by}' not found in metadata.\n"
            f"  Available columns: {', '.join(metadata.columns.tolist())}"
        )

    artifacts = [Path(a) for a in args.artifact]
    for a in artifacts:
        if not a.exists():
            raise FileNotFoundError(f"Artifact not found: {a}")

    all_data = []
    for artifact_path in artifacts:
        print(f"\n  Reading: {artifact_path.name}")
        series, metric_name = read_alpha_artifact(artifact_path)
        matched_meta = match_sample_ids(series.index, metadata)

        if args.group_by not in matched_meta.columns:
            raise ValueError(
                f"Column '{args.group_by}' not found after metadata matching."
            )
        group_col = matched_meta[args.group_by]
        all_data.append((series, metric_name, group_col, artifact_path.stem))

    palette     = _resolve_palette(args.palette)
    group_order = getattr(args, "group_order", None)

    if args.panel and len(all_data) > 1:
        n_groups = len(all_data[0][2].unique())
        n        = len(all_data)
        ncols    = min(4, n)
        nrows    = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(max(3, n_groups * 1.6) * ncols, 5 * nrows),
        )
        axes_flat = np.array(axes).flatten()

        for i, (series, metric_name, group_col, stem) in enumerate(all_data):
            ax    = axes_flat[i]
            label = "ABCDEFGHIJ"[i]
            # Panel letter in top-left corner
            ax.text(
                0.02, 0.98, label, transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="top", ha="left", zorder=10,
                bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.85),
            )
            plot_alpha(
                series, group_col, metric_name, args.group_by,
                ax=ax, show_stats=not args.no_stats,
                palette=palette,
                show_title=not args.no_title,
                group_order=group_order,
            )

        for j in range(len(all_data), len(axes_flat)):
            axes_flat[j].set_visible(False)

        # Add legend to the last visible axis
        last_series, last_metric, last_group_col, last_stem = all_data[-1]
        if group_order:
            last_group_col = last_group_col[last_group_col.isin(group_order)]
        color_dict_last, group_labels_last = _get_color_map(last_group_col, palette)
        if group_order:
            group_labels_last = [c for c in group_order if c in group_labels_last]
        handles = [
            mpatches.Patch(facecolor=color_dict_last[cat], label=cat)
            for cat in group_labels_last
        ]
        axes_flat[len(all_data) - 1].legend(
            handles=handles,
            title=args.group_by,
            fontsize=FONT_SIZE_LEGEND,
            title_fontsize=FONT_SIZE_LEGEND,
            framealpha=0.8,
            edgecolor="#cccccc",
            loc="upper right",
        )

        if not args.no_title:
            fig.suptitle(
                f"Alpha Diversity — grouped by {args.group_by}",
                fontsize=FONT_SIZE_TITLE + 1, fontweight="bold", y=1.01,
            )
        fig.tight_layout()
        out_stem = args.output_stem or f"alpha_panel_{args.group_by}"
        _save_figure(fig, out_dir, out_stem)
        plt.close(fig)

    else:
        for series, metric_name, group_col, stem in all_data:
            fig = plot_alpha(
                series, group_col, metric_name, args.group_by,
                group_order=group_order,
                show_stats=not args.no_stats,
                palette=palette,
                show_title=not args.no_title,
            )
            out_stem = args.output_stem or f"alpha_{stem}_{args.group_by}"
            _save_figure(fig, out_dir, out_stem)
            plt.close(fig)

    print(f"\n  All alpha figures saved to: {out_dir.resolve()}")


# ===========================================================================
# list-columns helper
# ===========================================================================

def cmd_list_columns(args: argparse.Namespace) -> None:
    """Print available metadata columns and their unique values."""
    metadata = load_metadata(Path(args.metadata))
    print(f"\nMetadata file: {args.metadata}")
    print(f"Samples: {len(metadata)}")
    print(f"\nAvailable columns for --color-by / --group-by:\n")
    for col in metadata.columns:
        unique_vals = metadata[col].dropna().unique()
        if len(unique_vals) <= 12:
            print(f"  {col:<30} {sorted(unique_vals, key=str)}")
        else:
            print(f"  {col:<30} ({len(unique_vals)} unique values — continuous/ID column)")


# ===========================================================================
# Argument parsing
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    """
    Build and return the top-level argument parser for 09_plot_diversity.py.

    Two subcommands are available:
      pcoa    — 2D PCoA scatter plots from beta diversity QZAs, with optional
                95% confidence ellipses and --panel mode for multi-metric figures.
      alpha   — Strip/box plots from alpha diversity vector QZAs, with
                non-parametric significance annotations (Mann-Whitney U or
                Kruskal-Wallis).
      columns — List available metadata columns and unique values.
    """
    p = argparse.ArgumentParser(
        prog="plot_diversity.py",
        description=(
            "Generate publication-ready diversity figures from QIIME2 artifacts.\n"
            "No QIIME2 installation required — reads .qza files directly.\n"
            "Outputs: PNG (300 dpi) + SVG (Illustrator-editable)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n\n"
            "  # Check what columns are available in your metadata\n"
            "  python 09_plot_diversity.py columns \\\n"
            "    --metadata metadata/qiime/metadata_16S.tsv\n\n"
            "  # 2D PCoA panel of all four beta metrics, colored by Group\n"
            "  python 09_plot_diversity.py pcoa \\\n"
            "    --artifact weighted_unifrac_pcoa_results.qza \\\n"
            "               unweighted_unifrac_pcoa_results.qza \\\n"
            "               bray_curtis_pcoa_results.qza \\\n"
            "               jaccard_pcoa_results.qza \\\n"
            "    --metadata metadata/qiime/metadata_16S.tsv \\\n"
            "    --color-by Group --panel \\\n"
            "    --group-order Diseased Trauma Marine \\\n"
            "    --output-dir results/figures/\n\n"
            "  # Alpha diversity panel, all four metrics\n"
            "  python 09_plot_diversity.py alpha \\\n"
            "    --artifact faith_pd_vector.qza shannon_vector.qza \\\n"
            "               observed_features_vector.qza evenness_vector.qza \\\n"
            "    --metadata metadata/qiime/metadata_16S.tsv \\\n"
            "    --group-by Group --panel \\\n"
            "    --group-order Diseased Trauma \\\n"
            "    --output-dir results/figures/\n"
        ),
    )

    sp = p.add_subparsers(dest="command", required=True, metavar="COMMAND")

    # ── pcoa ──────────────────────────────────────────────────────────────────
    pc = sp.add_parser("pcoa", help="2D PCoA scatter plots from PCoA results artifacts")
    pc.add_argument(
        "--artifact", nargs="+", required=True, metavar="QZA",
        help="One or more PCoA results .qza files",
    )
    pc.add_argument(
        "--metadata", required=True, metavar="CSV/TSV",
        help="Metadata file (CSV or TSV). Must contain a TV or SampleID column.",
    )
    pc.add_argument(
        "--color-by", required=True, metavar="COLUMN",
        help="Metadata column to color points by (e.g. Group, 'Cadaver Condition')",
    )
    pc.add_argument(
        "--group-order", nargs="+", default=None,
        help=(
            "Explicit ordering of groups for legend and color assignment "
            "(e.g. Diseased Trauma Marine). Default: alphabetical. "
            "Groups absent from data are silently skipped."
        ),
    )
    pc.add_argument(
        "--pc-axes", default="1,2", metavar="X,Y",
        help="Which PC axes to plot, 1-indexed (default: 1,2 — PC1 vs PC2)",
    )
    pc.add_argument(
        "--panel", action="store_true",
        help="Combine all artifacts into one panel figure with A/B/C/D labels",
    )
    pc.add_argument(
        "--no-ellipses", action="store_true",
        help="Disable 95%% confidence ellipses around groups",
    )
    pc.add_argument(
        "--output-dir", default="figures", metavar="DIR",
        help="Output directory for figures (default: figures/)",
    )
    pc.add_argument(
        "--palette", default=None, metavar="PALETTE",
        help=(
            "Named palette: purple, redblue, or wong (Wong 2011 colorblind-safe). "
            "Or pass raw hex codes: '#B22222,#2E86C1'. Default: wong."
        ),
    )
    pc.add_argument(
        "--stats-dir", default=None, metavar="DIR",
        help=(
            "Directory containing PERMANOVA QZVs (from 08_run_diversity_stats.py). "
            "If provided, F-stat and p-value are annotated on each PCoA panel. "
            "Note: always report PERMDISP alongside PERMANOVA in the manuscript."
        ),
    )
    pc.add_argument(
        "--output-stem", default=None, metavar="STEM",
        help="Base filename stem (no extension). Default: pcoa_panel_{color_by}.",
    )
    pc.add_argument(
        "--no-title", action="store_true",
        help="Suppress subplot titles (use for journal submission figures)",
    )
    pc.set_defaults(func=cmd_pcoa)

    # ── alpha ─────────────────────────────────────────────────────────────────
    al = sp.add_parser(
        "alpha", help="Alpha diversity strip/box plots from alpha vector artifacts"
    )
    al.add_argument(
        "--artifact", nargs="+", required=True, metavar="QZA",
        help="One or more alpha diversity vector .qza files",
    )
    al.add_argument(
        "--metadata", required=True, metavar="CSV/TSV",
        help="Metadata file (CSV or TSV)",
    )
    al.add_argument(
        "--group-by", required=True, metavar="COLUMN",
        help="Metadata column to group samples by (e.g. Group, Season)",
    )
    al.add_argument(
        "--group-order", nargs="+", default=None,
        help=(
            "Explicit ordering of groups for x-axis "
            "(e.g. Diseased Trauma). Default: alphabetical."
        ),
    )
    al.add_argument(
        "--panel", action="store_true",
        help="Combine all metrics into one panel figure with A/B/C/D labels",
    )
    al.add_argument(
        "--no-stats", action="store_true",
        help="Disable statistical annotation (Mann-Whitney U / Kruskal-Wallis)",
    )
    al.add_argument(
        "--output-dir", default="figures", metavar="DIR",
        help="Output directory for figures (default: figures/)",
    )
    al.add_argument(
        "--palette", default=None, metavar="PALETTE",
        help=(
            "Named palette: purple, redblue, or wong (Wong 2011 colorblind-safe). "
            "Or pass raw hex codes: '#B22222,#2E86C1'. Default: wong."
        ),
    )
    al.add_argument(
        "--output-stem", default=None, metavar="STEM",
        help="Base filename stem (no extension). Default: alpha_panel_{group_by}.",
    )
    al.add_argument(
        "--no-title", action="store_true",
        help="Suppress subplot titles (use for journal submission figures)",
    )
    al.set_defaults(func=cmd_alpha)

    # ── columns ───────────────────────────────────────────────────────────────
    co = sp.add_parser(
        "columns",
        help="List available metadata columns and their unique values",
    )
    co.add_argument("--metadata", required=True, metavar="CSV/TSV", help="Metadata file")
    co.set_defaults(func=cmd_list_columns)

    return p


# ===========================================================================
# Entry point
# ===========================================================================

def main(argv=None) -> int:
    """
    Entry point for 09_plot_diversity.py.

    Parses the subcommand (pcoa, alpha, or columns), dispatches to the
    appropriate function, and writes PNG + SVG outputs to --output-dir.
    Returns 0 on success, 1 on any error.
    """
    parser = build_parser()
    args   = parser.parse_args(argv)

    try:
        args.func(args)
        return 0

    except KeyboardInterrupt:
        print("\n[INTERRUPTED]")
        return 130

    except FileNotFoundError as exc:
        print(f"\n[ERROR] File not found:\n  {exc}", file=sys.stderr)
        return 1

    except ValueError as exc:
        print(f"\n[ERROR] Invalid input:\n  {exc}", file=sys.stderr)
        return 1

    except OSError as exc:
        print(f"\n[ERROR] File system error:\n  {exc}", file=sys.stderr)
        return 1

    except Exception as exc:
        print(
            f"\n[UNEXPECTED ERROR] {type(exc).__name__}: {exc}\n"
            f"Please report this with the full traceback below.",
            file=sys.stderr,
        )
        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
