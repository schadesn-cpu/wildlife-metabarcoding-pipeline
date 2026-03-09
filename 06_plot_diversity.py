#!/usr/bin/env python3
"""
plot_diversity.py

Standalone diversity visualization script for loon metabarcoding data.
Reads QIIME2 artifacts directly (.qza files) — no QIIME2 installation needed.

Generates:
  - 2D PCoA plots (beta diversity) — one per distance matrix artifact
  - Alpha diversity strip/box plots — one per alpha vector artifact

Output formats: PNG (slides/SharePoint) + SVG (publication/Illustrator editing)

Usage:
  # Single beta metric, colored by Group
  python plot_diversity.py pcoa \\
    --artifact weighted_unifrac_pcoa_results.qza \\
    --metadata full_hcgs_sample_metadata.csv \\
    --color-by Group \\
    --output-dir figures/

  # All beta metrics at once
  python plot_diversity.py pcoa \\
    --artifact weighted_unifrac_pcoa_results.qza unweighted_unifrac_pcoa_results.qza \\
                bray_curtis_pcoa_results.qza jaccard_pcoa_results.qza \\
    --metadata full_hcgs_sample_metadata.csv \\
    --color-by Group \\
    --output-dir figures/

  # Alpha diversity
  python plot_diversity.py alpha \\
    --artifact faith_pd_vector.qza shannon_vector.qza \\
               observed_features_vector.qza evenness_vector.qza \\
    --metadata full_hcgs_sample_metadata.csv \\
    --group-by Group \\
    --output-dir figures/

  # Panel of all four beta metrics in one figure
  python plot_diversity.py pcoa \\
    --artifact weighted_unifrac_pcoa_results.qza unweighted_unifrac_pcoa_results.qza \\
                bray_curtis_pcoa_results.qza jaccard_pcoa_results.qza \\
    --metadata full_hcgs_sample_metadata.csv \\
    --color-by Group --panel \\
    --output-dir figures/
"""
from __future__ import annotations

import argparse
import csv
import io
import json
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

warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Style constants — edit these to match your publication style
# ---------------------------------------------------------------------------

FIGURE_DPI = 300
POINT_SIZE = 80
POINT_ALPHA = 0.85
FONT_FAMILY = "Arial"
FONT_SIZE_TITLE = 13
FONT_SIZE_AXIS = 11
FONT_SIZE_TICK = 9
FONT_SIZE_LEGEND = 9
FONT_SIZE_ANNOT = 8
ELLIPSE_ALPHA = 0.12        # transparency of 95% confidence ellipses
GRID_ALPHA = 0.25

# Colorblind-friendly palette (Wong 2011)
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
# Named palettes — edit here to change figures globally
# Use --palette purple / --palette redblue / --palette wong on the CLI,
# or pass raw hex codes as before: --palette "#B22222,#2E86C1"
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
        "#B22222",  # red          — recommended for multi-marker or colorblind
        "#2E86C1",  # blue
        "#E74C3C",  # light red    — Group 3
        "#7FB3D3",  # light blue   — Group 4
        "#7B241C", "#1A5276", "#F1948A", "#2980B9",
    ],
    "wong": [
        "#0072B2",  # blue         — Wong 2011 colorblind-safe (citable)
        "#E69F00",  # orange
        "#009E73",  # green
        "#CC79A7",  # pink
        "#56B4E9",  # sky blue
        "#D55E00",  # vermillion
        "#F0E442",  # yellow
        "#000000",  # black
    ],
}


def _resolve_palette(palette_arg: Optional[str]) -> List[str]:
    """
    Resolve --palette argument to a list of hex color strings.
    Accepts either a named palette ('purple', 'redblue', 'wong')
    or a comma-separated list of hex codes ('#B22222,#2E86C1').
    Falls back to PALETTE_DEFAULT (Wong) if None.
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
    """Find a file by name anywhere inside the zip."""
    for name in zf.namelist():
        if name.endswith(f"/data/{filename}") or name.endswith(f"/{filename}"):
            return name
    return None


def read_pcoa_artifact(path: Path) -> Tuple[pd.DataFrame, List[float], str]:
    """
    Read a QIIME2 PCoA results artifact (.qza).
    Returns:
      coords_df   : DataFrame with sample IDs as index, PC columns
      prop_expl   : List of proportion explained per axis
      metric_name : Guessed metric name from filename
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
        raise ValueError(f"Could not open {path.name} as a zip file: {exc}") from exc

    coords, prop_expl = _parse_ordination(content)
    metric_name = _guess_metric_name(path.name)
    return coords, prop_expl, metric_name


def _parse_ordination(content: str) -> Tuple[pd.DataFrame, List[float]]:
    """Parse QIIME2 ordination.txt format."""
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
                raise ValueError(f"Could not parse 'Proportion explained' section: {exc}") from exc

        elif header.startswith("Site"):
            try:
                parts = header.split("\t")
                n_samples = int(parts[1])
                n_axes = int(parts[2])
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Could not parse 'Site' header: {exc}") from exc

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
    Returns:
      series      : Series with sample IDs as index, diversity values
      metric_name : Metric name from the TSV header
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
        raise ValueError(f"Could not open {path.name} as a zip file: {exc}") from exc

    try:
        df = pd.read_csv(io.StringIO(content), sep="\t", index_col=0)
    except Exception as exc:
        raise ValueError(f"Could not parse alpha-diversity.tsv from {path.name}: {exc}") from exc

    if df.empty or df.shape[1] == 0:
        raise ValueError(f"Alpha diversity TSV appears empty in {path.name}")

    metric_col = df.columns[0]
    series = df[metric_col].copy()
    series.index.name = "SampleID"
    metric_name = _pretty_metric_name(metric_col)
    return series, metric_name



def read_permanova_qzv(qzv_path: Path) -> Optional[Dict]:
    """
    Parse a QIIME2 beta-group-significance QZV (PERMANOVA or PERMDISP).
    Returns dict with keys: method, pseudo_f, p_value, group_column
    or None if parsing fails.
    """
    import re
    try:
        with zipfile.ZipFile(qzv_path) as zf:
            names = zf.namelist()

            # Try to find a results TSV/CSV first
            for name in names:
                if name.endswith("permanova_results.csv") or name.endswith("results.tsv"):
                    with zf.open(name) as f:
                        content = f.read().decode("utf-8", errors="replace")
                    # Parse first data row
                    lines = [l for l in content.splitlines() if l.strip() and not l.startswith("#")]
                    if len(lines) >= 2:
                        header = lines[0].split(",")
                        vals = lines[1].split(",")
                        row = dict(zip(header, vals))
                        pseudo_f = float(row.get("pseudo-F", row.get("test statistic", 0)))
                        p_val = float(row.get("p-value", row.get("p_value", 1)))
                        return {"pseudo_f": pseudo_f, "p_value": p_val}

            # Fallback: parse index.html
            html_files = [n for n in names if n.endswith("index.html")]
            if not html_files:
                return None
            with zf.open(html_files[0]) as f:
                html = f.read().decode("utf-8", errors="replace")

            # Look for pseudo-F and p-value in the HTML table
            # QIIME2 renders these as table rows
            pseudo_f = None
            p_value = None

            # Pattern: numbers in table cells near "pseudo-F" or "test statistic"
            f_match = re.search(
                r"pseudo-F.*?<td[^>]*>\s*([0-9]+\.?[0-9]*(?:e[+-]?[0-9]+)?)",
                html, re.IGNORECASE | re.DOTALL
            )
            p_match = re.search(
                r"p-value.*?<td[^>]*>\s*([0-9]+\.?[0-9]*(?:e[+-]?[0-9]+)?)",
                html, re.IGNORECASE | re.DOTALL
            )
            
            # Alternative: look for the data embedded as JSON/JS variable
            json_match = re.search(r'"pseudo-F"\s*:\s*([0-9.eE+\-]+)', html)
            pjson_match = re.search(r'"p-value"\s*:\s*([0-9.eE+\-]+)', html)
            
            if json_match:
                pseudo_f = float(json_match.group(1))
            elif f_match:
                pseudo_f = float(f_match.group(1))
                
            if pjson_match:
                p_value = float(pjson_match.group(1))
            elif p_match:
                p_value = float(p_match.group(1))

            if pseudo_f is not None and p_value is not None:
                return {"pseudo_f": pseudo_f, "p_value": p_value}

    except Exception as e:
        pass
    return None


# ===========================================================================
# Metadata loading
# ===========================================================================

def load_metadata(path: Path) -> pd.DataFrame:
    """
    Load metadata CSV or TSV.
    Handles:
      - BOM characters (Excel exports)
      - QIIME2 TSV format (#SampleID)
      - CSV with TV column for matching to full sample IDs
    """
    if not path.exists():
        raise FileNotFoundError(f"Metadata file not found: {path}")

    suffix = path.suffix.lower()
    sep = "\t" if suffix in (".tsv", ".txt") else ","

    try:
        df = pd.read_csv(path, sep=sep, dtype=str, encoding="utf-8-sig")
    except Exception as exc:
        raise ValueError(f"Could not read metadata file {path}: {exc}") from exc

    # Normalize column names — strip whitespace, handle BOM, handle #SampleID
    df.columns = [c.strip().lstrip("\ufeff").lstrip("#") for c in df.columns]

    # Find the sample ID column
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

    # Strip whitespace from all values
    for col in df.columns:
        df[col] = df[col].str.strip() if df[col].dtype == object else df[col]

    return df


def match_sample_ids(
    artifact_ids: pd.Index,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """
    Match artifact sample IDs to metadata.

    Handles the common QIIME2 pattern where artifact IDs look like
    'TV230007-GI-16S_S1483' but metadata uses short IDs like 'TV230007'.

    Strategy:
      1. Direct match (exact)
      2. Prefix match — artifact ID starts with metadata ID
      3. Suffix match — metadata ID is contained in artifact ID
    """
    # 1. Direct match
    direct = artifact_ids.intersection(metadata.index)
    if len(direct) == len(artifact_ids):
        return metadata.loc[artifact_ids]

    # 2. Build a mapping: artifact_id → metadata_id via prefix
    mapping: Dict[str, str] = {}
    meta_ids = metadata.index.tolist()

    for art_id in artifact_ids:
        # Direct
        if art_id in metadata.index:
            mapping[art_id] = art_id
            continue
        # Prefix: artifact ID starts with metadata ID (e.g. TV230007-GI-16S...)
        matched = [m for m in meta_ids if art_id.startswith(m)]
        if matched:
            # Take longest match
            mapping[art_id] = max(matched, key=len)
            continue
        # Substring: metadata ID is contained anywhere in artifact ID
        matched = [m for m in meta_ids if m in art_id]
        if matched:
            mapping[art_id] = max(matched, key=len)
            continue
        # No match — will be NaN in output
        mapping[art_id] = None

    unmatched = [k for k, v in mapping.items() if v is None]
    if unmatched:
        print(
            f"  ⚠  Warning: {len(unmatched)} sample(s) in artifact not found in metadata:\n"
            + "\n".join(f"       {u}" for u in unmatched[:5])
            + (f"\n       ... and {len(unmatched)-5} more" if len(unmatched) > 5 else "")
        )

    matched_meta_ids = [mapping[a] for a in artifact_ids if mapping[a] is not None]
    art_ids_matched = [a for a in artifact_ids if mapping[a] is not None]

    result = metadata.loc[matched_meta_ids].copy()
    result.index = art_ids_matched
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
    """Convert snake_case metric names to readable labels."""
    mapping = {
        "weighted_unifrac": "Weighted UniFrac",
        "unweighted_unifrac": "Unweighted UniFrac",
        "bray_curtis": "Bray-Curtis",
        "jaccard": "Jaccard",
        "faith_pd": "Faith's PD",
        "faith_pd_vector": "Faith's PD",
        "shannon": "Shannon Diversity",
        "shannon_entropy": "Shannon Diversity",
        "observed_features": "Observed Features",
        "pielou_evenness": "Pielou's Evenness",
        "evenness": "Pielou's Evenness",
    }
    return mapping.get(name, name.replace("_", " ").title())


def _get_color_map(
    values: pd.Series,
    palette: Optional[List[str]] = None,
) -> Tuple[Dict[str, str], List[str]]:
    """
    Build a color mapping for unique values in a series.
    Returns (color_dict, ordered_categories).
    """
    if palette is None:
        palette = PALETTE_DEFAULT

    # Try to treat as categorical
    unique_vals = sorted(values.dropna().unique(), key=str)

    if len(unique_vals) > len(palette):
        # Fall back to matplotlib colormap for many categories
        cmap = plt.cm.get_cmap("tab20", len(unique_vals))
        palette = [matplotlib.colors.to_hex(cmap(i)) for i in range(len(unique_vals))]

    color_dict = {str(v): palette[i % len(palette)] for i, v in enumerate(unique_vals)}
    return color_dict, [str(v) for v in unique_vals]


def _confidence_ellipse(
    x: np.ndarray,
    y: np.ndarray,
    ax: plt.Axes,
    n_std: float = 2.0,
    facecolor: str = "none",
    **kwargs,
) -> None:
    """
    Draw a covariance confidence ellipse around a set of 2D points.
    n_std=2.0 approximates the 95% confidence region.
    """
    if len(x) < 3:
        return  # Can't draw ellipse with fewer than 3 points

    cov = np.cov(x, y)
    if np.any(np.isnan(cov)) or np.linalg.det(cov) == 0:
        return

    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
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
    """Save figure as both PNG and SVG."""
    out_dir.mkdir(parents=True, exist_ok=True)

    png_path = out_dir / f"{stem}.png"
    svg_path = out_dir / f"{stem}.svg"

    fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, format="svg", bbox_inches="tight", facecolor="white")

    print(f"  ✓ Saved: {png_path.name}")
    print(f"  ✓ Saved: {svg_path.name}")


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
) -> plt.Figure:
    """
    Draw a 2D PCoA scatter plot.

    Parameters
    ----------
    coords        : DataFrame of PC coordinates (index=SampleID)
    prop_expl     : Proportion of variance explained per axis
    metadata_col  : Series mapping SampleID → group label
    metric_name   : e.g. 'Weighted UniFrac' for axis/title labels
    color_col     : Metadata column name used for coloring
    pc_axes       : Which PC axes to plot (1-indexed), default (1,2)
    ellipses      : Draw 95% confidence ellipses around groups
    ax            : Existing axes to draw into (for panel figures)
    permanova_result : Optional dict with 'pseudo_f' and 'p_value' to annotate
    """
    pc_x = f"PC{pc_axes[0]}"
    pc_y = f"PC{pc_axes[1]}"

    # Align metadata with coords
    common = coords.index.intersection(metadata_col.index)
    if len(common) == 0:
        raise ValueError(
            f"No overlapping sample IDs between PCoA coords and metadata column '{color_col}'.\n"
            f"  PCoA IDs (first 3): {list(coords.index[:3])}\n"
            f"  Metadata IDs (first 3): {list(metadata_col.index[:3])}"
        )
    coords_plot = coords.loc[common]
    meta_plot = metadata_col.loc[common].fillna("Unknown").astype(str)

    color_dict, categories = _get_color_map(meta_plot, palette)

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.get_figure()

    # Plot each group
    for i, group in enumerate(categories):
        mask = meta_plot == group
        if not mask.any():
            continue
        x = coords_plot.loc[mask, pc_x].values
        y = coords_plot.loc[mask, pc_y].values
        color = color_dict[group]
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
            _confidence_ellipse(
                x, y, ax,
                n_std=2.0,
                facecolor=color,
                alpha=ELLIPSE_ALPHA,
                edgecolor=color,
                linewidth=1.2,
                linestyle="--",
                zorder=2,
            )

    # Axis labels with variance explained
    pct_x = prop_expl[pc_axes[0] - 1] * 100 if len(prop_expl) >= pc_axes[0] else 0
    pct_y = prop_expl[pc_axes[1] - 1] * 100 if len(prop_expl) >= pc_axes[1] else 0
    ax.set_xlabel(f"{pc_x} ({pct_x:.1f}%)", fontsize=FONT_SIZE_AXIS)
    ax.set_ylabel(f"{pc_y} ({pct_y:.1f}%)", fontsize=FONT_SIZE_AXIS)

    # Title — add star if PERMANOVA significant
    sig = permanova_result and permanova_result.get("p_value", 1) <= 0.05
    plot_title = title or f"{metric_name}"
    if sig:
        plot_title = f"{plot_title} *"
    title_color = "#6A0DAD" if sig else "black"
    ax.set_title(plot_title, fontsize=FONT_SIZE_TITLE, fontweight="bold",
                 pad=8, color=title_color)

    # PERMANOVA annotation — match 16S style
    if permanova_result:
        pf = permanova_result.get("pseudo_f")
        pv = permanova_result.get("p_value")
        if pf is not None and pv is not None:
            star = " *" if pv <= 0.05 else ""
            annot = f"PERMANOVA  F={pf:.3f}, p={pv:.3f}{star}"
            annot_color = "#6A0DAD" if pv <= 0.05 else "#444444"
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

    # Grid and spines
    ax.axhline(0, color="grey", linewidth=0.5, alpha=0.4, zorder=1)
    ax.axvline(0, color="grey", linewidth=0.5, alpha=0.4, zorder=1)
    ax.grid(True, alpha=GRID_ALPHA, linewidth=0.5)
    ax.tick_params(labelsize=FONT_SIZE_TICK)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # Legend
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
    out_dir = Path(args.output_dir)
    metadata = load_metadata(Path(args.metadata))
    palette = _resolve_palette(args.palette)
    stats_dir = Path(args.stats_dir) if args.stats_dir else None

    if args.color_by not in metadata.columns:
        available = ", ".join(metadata.columns.tolist())
        raise ValueError(
            f"--color-by column '{args.color_by}' not found in metadata.\n"
            f"  Available columns: {available}"
        )

    pc_axes = tuple(int(x) for x in args.pc_axes.split(","))
    if len(pc_axes) != 2:
        raise ValueError(f"--pc-axes must be two comma-separated integers, e.g. '1,2', got: {args.pc_axes}")

    artifacts = [Path(a) for a in args.artifact]
    for a in artifacts:
        if not a.exists():
            raise FileNotFoundError(f"Artifact not found: {a}")

    # Load all artifacts
    all_data = []
    for artifact_path in artifacts:
        print(f"\n  Reading: {artifact_path.name}")
        coords, prop_expl, metric_name = read_pcoa_artifact(artifact_path)
        matched_meta = match_sample_ids(coords.index, metadata)
        color_col = matched_meta[args.color_by] if args.color_by in matched_meta.columns else None

        if color_col is None:
            raise ValueError(f"Column '{args.color_by}' not found after metadata matching.")

        # Load PERMANOVA result if stats_dir provided
        permanova_result = None
        if stats_dir:
            group_col = args.color_by
            # Try both naming conventions
            candidates = [
                stats_dir / f"{metric_name}_permanova_{group_col}.qzv",
                stats_dir / f"{metric_name}_permanova.qzv",
            ]
            for qzv in candidates:
                if qzv.exists():
                    permanova_result = read_permanova_qzv(qzv)
                    if permanova_result:
                        print(f"  PERMANOVA loaded: F={permanova_result['pseudo_f']:.3f}, p={permanova_result['p_value']:.3f}")
                    break

        all_data.append((coords, prop_expl, metric_name, color_col, artifact_path.stem, permanova_result))

    if args.panel and len(all_data) > 1:
        # Panel figure — all metrics in a grid
        n = len(all_data)
        ncols = min(2, n)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
        axes_flat = np.array(axes).flatten()

        for i, (coords, prop_expl, metric_name, color_col, stem, permanova_result) in enumerate(all_data):
            ax = axes_flat[i]
            plot_pcoa(
                coords, prop_expl, color_col, metric_name,
                args.color_by, pc_axes=pc_axes,
                ellipses=not args.no_ellipses,
                palette=palette,
                permanova_result=permanova_result,
                ax=ax,
            )
            # Only show legend on last subplot
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

        # Hide any unused axes
        for j in range(len(all_data), len(axes_flat)):
            axes_flat[j].set_visible(False)

        fig.suptitle(
            f"Beta Diversity PCoA — colored by {args.color_by}",
            fontsize=FONT_SIZE_TITLE + 1,
            fontweight="bold",
            y=1.01,
        )
        fig.tight_layout()
        stem = f"pcoa_panel_{args.color_by}"
        _save_figure(fig, out_dir, stem)
        plt.close(fig)

    else:
        # Individual figures
        for coords, prop_expl, metric_name, color_col, stem, permanova_result in all_data:
            fig = plot_pcoa(
                coords, prop_expl, color_col, metric_name,
                args.color_by, pc_axes=pc_axes,
                ellipses=not args.no_ellipses,
                palette=palette,
                permanova_result=permanova_result,
            )
            out_stem = f"pcoa_{stem}_{args.color_by}"
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
) -> plt.Figure:
    """
    Draw an alpha diversity strip + box plot.

    Points are jittered over a box, with each sample visible.
    Optional Mann-Whitney U p-value annotated for 2-group comparisons,
    or Kruskal-Wallis for 3+ groups.
    """
    # Align
    common = series.index.intersection(metadata_col.index)
    if len(common) == 0:
        raise ValueError(
            f"No overlapping sample IDs between alpha diversity and metadata column '{group_col}'."
        )
    series_plot = series.loc[common]
    meta_plot = metadata_col.loc[common].fillna("Unknown").astype(str)

    color_dict, categories = _get_color_map(meta_plot, palette)

    # Build groups list
    group_data = []
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

    # Box plots (no outlier markers — individual points shown instead)
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

    # Strip plot — jittered individual points
    rng = np.random.default_rng(seed=42)
    for i, (vals, cat) in enumerate(zip(group_data, group_labels)):
        color = color_dict[cat]
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
    if show_stats and len(group_data) >= 2:
        if len(group_data) == 2:
            stat, p_val = stats.mannwhitneyu(
                group_data[0], group_data[1], alternative="two-sided"
            )
            stat_label = f"Mann-Whitney U\np = {p_val:.3f}"
            if p_val < 0.001:
                stat_label = f"Mann-Whitney U\np < 0.001"
        else:
            stat, p_val = stats.kruskal(*group_data)
            stat_label = f"Kruskal-Wallis\np = {p_val:.3f}"
            if p_val < 0.001:
                stat_label = f"Kruskal-Wallis\np < 0.001"

        y_max = max(np.max(v) for v in group_data)
        y_range = y_max - min(np.min(v) for v in group_data)
        ax.annotate(
            stat_label,
            xy=(0.98, 0.97),
            xycoords="axes fraction",
            fontsize=FONT_SIZE_ANNOT,
            color="#333333",
            style="italic",
            ha="right",
            va="top",
        )

        # Significance bracket for 2-group comparison
        if len(group_data) == 2:
            sig_str = _p_to_stars(p_val)
            if sig_str:
                y_bracket = y_max + y_range * 0.08
                ax.plot(
                    [1, 1, 2, 2],
                    [y_bracket, y_bracket + y_range * 0.03,
                     y_bracket + y_range * 0.03, y_bracket],
                    lw=1.2, color="#333333"
                )
                ax.text(
                    1.5, y_bracket + y_range * 0.04,
                    sig_str, ha="center", va="bottom",
                    fontsize=FONT_SIZE_AXIS, color="#333333"
                )

    # Axes formatting
    ax.set_xticks(positions)
    ax.set_xticklabels(group_labels, fontsize=FONT_SIZE_TICK)
    ax.set_ylabel(metric_name, fontsize=FONT_SIZE_AXIS)
    ax.set_title(metric_name, fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=8)
    ax.tick_params(axis="y", labelsize=FONT_SIZE_TICK)
    ax.grid(axis="y", alpha=GRID_ALPHA, linewidth=0.5)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # n label under each group
    for i, (vals, cat) in enumerate(zip(group_data, group_labels)):
        ax.text(
            positions[i], ax.get_ylim()[0],
            f"n={len(vals)}",
            ha="center", va="top",
            fontsize=FONT_SIZE_ANNOT,
            color="#555555",
        )

    if standalone:
        fig.tight_layout()

    return fig


def _p_to_stars(p: float) -> str:
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return ""


def cmd_alpha(args: argparse.Namespace) -> None:
    """Generate alpha diversity strip/box plots from one or more alpha vector artifacts."""
    out_dir = Path(args.output_dir)
    metadata = load_metadata(Path(args.metadata))

    if args.group_by not in metadata.columns:
        available = ", ".join(metadata.columns.tolist())
        raise ValueError(
            f"--group-by column '{args.group_by}' not found in metadata.\n"
            f"  Available columns: {available}"
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
        group_col = matched_meta[args.group_by] if args.group_by in matched_meta.columns else None

        if group_col is None:
            raise ValueError(f"Column '{args.group_by}' not found after metadata matching.")

        all_data.append((series, metric_name, group_col, artifact_path.stem))

    if args.panel and len(all_data) > 1:
        n = len(all_data)
        ncols = min(4, n)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(max(3, len(all_data[0][2].unique()) * 1.6) * ncols, 5 * nrows))
        axes_flat = np.array(axes).flatten()

        palette = _resolve_palette(args.palette)
        for i, (series, metric_name, group_col, stem) in enumerate(all_data):
            ax = axes_flat[i]
            plot_alpha(
                series, group_col, metric_name, args.group_by,
                ax=ax, show_stats=not args.no_stats,
                palette=palette,
            )

        for j in range(len(all_data), len(axes_flat)):
            axes_flat[j].set_visible(False)

        fig.suptitle(
            f"Alpha Diversity — grouped by {args.group_by}",
            fontsize=FONT_SIZE_TITLE + 1,
            fontweight="bold",
            y=1.01,
        )
        fig.tight_layout()
        out_stem = f"alpha_panel_{args.group_by}"
        _save_figure(fig, out_dir, out_stem)
        plt.close(fig)

    else:
        palette = _resolve_palette(args.palette)
        for series, metric_name, group_col, stem in all_data:
            fig = plot_alpha(
                series, group_col, metric_name, args.group_by,
                show_stats=not args.no_stats,
                palette=palette,
            )
            out_stem = f"alpha_{stem}_{args.group_by}"
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
    p = argparse.ArgumentParser(
        prog="plot_diversity.py",
        description=(
            "Generate publication-ready diversity figures from QIIME2 artifacts.\n"
            "No QIIME2 installation required — reads .qza files directly.\n"
            "Outputs: PNG (slides) + SVG (Illustrator-editable)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n\n"
            "  # Check what columns are available in your metadata\n"
            "  python plot_diversity.py columns \\\n"
            "    --metadata full_hcgs_sample_metadata.csv\n\n"
            "  # 2D PCoA panel of all four beta metrics, colored by Group\n"
            "  python plot_diversity.py pcoa \\\n"
            "    --artifact weighted_unifrac_pcoa_results.qza \\\n"
            "               unweighted_unifrac_pcoa_results.qza \\\n"
            "               bray_curtis_pcoa_results.qza \\\n"
            "               jaccard_pcoa_results.qza \\\n"
            "    --metadata full_hcgs_sample_metadata.csv \\\n"
            "    --color-by Group --panel \\\n"
            "    --output-dir figures/\n\n"
            "  # Same PCoA but colored by Cadaver Condition\n"
            "  python plot_diversity.py pcoa \\\n"
            "    --artifact weighted_unifrac_pcoa_results.qza \\\n"
            "    --metadata full_hcgs_sample_metadata.csv \\\n"
            "    --color-by 'Cadaver Condition' \\\n"
            "    --output-dir figures/\n\n"
            "  # Alpha diversity panel, all four metrics\n"
            "  python plot_diversity.py alpha \\\n"
            "    --artifact faith_pd_vector.qza shannon_vector.qza \\\n"
            "               observed_features_vector.qza evenness_vector.qza \\\n"
            "    --metadata full_hcgs_sample_metadata.csv \\\n"
            "    --group-by Group --panel \\\n"
            "    --output-dir figures/\n"
        ),
    )

    sp = p.add_subparsers(dest="command", required=True, metavar="COMMAND")

    # ---- pcoa ----
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
        help="Metadata column to color points by (e.g. Group, 'Cadaver Condition', Location)",
    )
    pc.add_argument(
        "--pc-axes", default="1,2", metavar="X,Y",
        help="Which PC axes to plot, 1-indexed (default: 1,2 — i.e. PC1 vs PC2)",
    )
    pc.add_argument(
        "--panel", action="store_true",
        help="Combine all artifacts into one panel figure instead of separate files",
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
            "Or pass raw hex codes: '#B22222,#2E86C1'. "
            "Default: wong. "
            "purple = dark purple/lavender (single-marker manuscript). "
            "redblue = red/blue (multi-marker or colorblind). "
            "wong = citable 8-color colorblind-safe palette."
        ),
    )
    pc.add_argument(
        "--stats-dir", default=None, metavar="DIR",
        help="Directory containing PERMANOVA QZVs. If provided, F-stat and p-value "
             "are annotated on each PCoA panel.",
    )
    pc.set_defaults(func=cmd_pcoa)

    # ---- alpha ----
    al = sp.add_parser("alpha", help="Alpha diversity strip/box plots from alpha vector artifacts")
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
        help="Metadata column to group samples by (e.g. Group, Location)",
    )
    al.add_argument(
        "--panel", action="store_true",
        help="Combine all metrics into one panel figure",
    )
    al.add_argument(
        "--no-stats", action="store_true",
        help="Disable statistical annotation (Mann-Whitney U / Kruskal-Wallis p-values)",
    )
    al.add_argument(
        "--output-dir", default="figures", metavar="DIR",
        help="Output directory for figures (default: figures/)",
    )
    al.add_argument(
        "--palette", default=None, metavar="PALETTE",
        help=(
            "Named palette: purple, redblue, or wong (Wong 2011 colorblind-safe). "
            "Or pass raw hex codes: '#B22222,#2E86C1'. "
            "Default: wong. "
            "purple = dark purple/lavender (single-marker manuscript). "
            "redblue = red/blue (multi-marker or colorblind). "
            "wong = citable 8-color colorblind-safe palette."
        ),
    )
    al.set_defaults(func=cmd_alpha)

    # ---- columns ----
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
    parser = build_parser()
    args = parser.parse_args(argv)

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

    except Exception as exc:
        print(f"\n[UNEXPECTED ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
