#!/usr/bin/env python3
"""
09c_visualize_diversity.py
=========================
Generate publication-quality alpha and beta diversity figures directly
from QIIME 2 QZA files. No QIIME 2 installation required at plot time:
QZA files are ZIP archives that this script reads directly.

Beta diversity:
  - PCoA scatter plots for each distance metric
  - 95% confidence ellipses per group
  - PERMANOVA statistic and p-value shown in each panel title
  - Stats sourced from a pre-computed summary TSV, selected by the
    group_column that matches the active --group-column, or computed
    on-the-fly from distance matrix QZAs.

Dependencies:
  pip install matplotlib numpy scipy pandas
"""

from __future__ import annotations

import argparse
import csv
import io
import logging
import sys
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Palettes
# ---------------------------------------------------------------------------

PALETTES = {
    "purple": {
        "colors": [
            "#7B2D8B", "#C19FD8", "#4B1369", "#D09EE0",
            "#2D0A40", "#E0BAEC", "#9870B0", "#F0D6F5",
        ],
        "markers": ["o", "s", "^", "D", "v", "P", "X", "*"],
        "sig_color": "#4A0060",
        "ns_color": "#888888",
    },
    "redblue": {
        "colors": [
            "#B22222", "#2E86C1", "#E74C3C", "#7FB3D3",
            "#7B241C", "#1A5276", "#F1948A", "#2980B9",
        ],
        "markers": ["o", "s", "^", "D", "v", "P", "X", "*"],
        "sig_color": "#1A1A1A",
        "ns_color": "#777777",
    },
    "wong": {
        # Wong 2011 citable colorblind-safe palette.
        "colors": [
            "#0072B2",  # blue
            "#E69F00",  # orange
            "#009E73",  # green
            "#CC79A7",  # pink
            "#56B4E9",  # sky blue
            "#D55E00",  # vermillion
            "#F0E442",  # yellow
            "#000000",  # black
        ],
        "markers": ["o", "s", "^", "D", "v", "P", "X", "*"],
        "sig_color": "#000000",
        "ns_color": "#777777",
    },
}


# ---------------------------------------------------------------------------
# QZA reading helpers
# ---------------------------------------------------------------------------

def _find_in_zip(zf: zipfile.ZipFile, suffix: str) -> Optional[str]:
    """Find the first entry in a ZipFile whose name ends with suffix."""
    return next((n for n in zf.namelist() if n.endswith(suffix)), None)


def read_alpha_vector(qza_path: Path) -> pd.Series:
    """Read an alpha diversity vector QZA. Returns pd.Series indexed by sample-id."""
    with zipfile.ZipFile(qza_path) as zf:
        tsv_name = _find_in_zip(zf, "alpha-diversity.tsv")
        if tsv_name is None:
            raise FileNotFoundError(f"alpha-diversity.tsv not found in {qza_path}")
        with zf.open(tsv_name) as fh:
            df = pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"), sep="\t", index_col=0)
    return df.iloc[:, 0]


def read_pcoa_ordination(qza_path: Path) -> Tuple[List[float], Dict[str, List[float]]]:
    """Read a PCoA ordination QZA. Returns (proportion_explained, {sample_id: [pc1, pc2, ...]})."""
    with zipfile.ZipFile(qza_path) as zf:
        txt_name = _find_in_zip(zf, "ordination.txt")
        if txt_name is None:
            raise FileNotFoundError(f"ordination.txt not found in {qza_path}")
        with zf.open(txt_name) as fh:
            content = io.TextIOWrapper(fh, encoding="utf-8").read()

    lines = content.strip().split("\n")
    prop_exp: List[float] = []
    coords: Dict[str, List[float]] = {}

    for i, line in enumerate(lines):
        if line.startswith("Proportion explained"):
            prop_exp = [float(v) * 100 for v in lines[i + 1].split("\t")]
        if line.startswith("Site\t"):
            n = int(line.split("\t")[1])
            for j in range(i + 1, i + 1 + n):
                parts = lines[j].split("\t")
                coords[parts[0]] = [float(x) for x in parts[1:]]

    return prop_exp, coords


def read_distance_matrix(qza_path: Path) -> pd.DataFrame:
    """Read a distance matrix QZA. Returns a square DataFrame."""
    with zipfile.ZipFile(qza_path) as zf:
        tsv_name = _find_in_zip(zf, "distance-matrix.tsv")
        if tsv_name is None:
            raise FileNotFoundError(f"distance-matrix.tsv not found in {qza_path}")
        with zf.open(tsv_name) as fh:
            dm = pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"), sep="\t", index_col=0)
    return dm


def infer_metric_name(qza_path: Path) -> str:
    """Guess a human-readable metric name from the QZA filename."""
    stem = qza_path.stem.lower()
    MAP = {
        "unweighted_unifrac": "Unweighted UniFrac",
        "weighted_unifrac":   "Weighted UniFrac",
        "bray_curtis":        "Bray-Curtis",
        "jaccard":            "Jaccard",
        "faith_pd":           "Faith's PD",
        "faith":              "Faith's PD",
        "shannon":            "Shannon Entropy",
        "observed_features":  "Observed Features",
        "evenness":           "Pielou's Evenness",
        "pielou":             "Pielou's Evenness",
    }
    for key, label in MAP.items():
        if key in stem:
            return label
    return qza_path.stem.replace("_", " ").title()


# ---------------------------------------------------------------------------
# Metadata loading
# ---------------------------------------------------------------------------

def load_metadata(metadata_path: Path, group_column: str) -> pd.Series:
    """Load a QIIME 2 metadata TSV and return a Series {sample_id: group_value}."""
    with metadata_path.open(encoding="utf-8") as f:
        lines = f.readlines()

    header = lines[0]
    data_lines = [l for l in lines[1:] if not l.startswith("#")]

    content = header + "".join(data_lines)
    df = pd.read_csv(io.StringIO(content), sep="\t", index_col=0, dtype=str)
    df.index.name = "sample-id"

    if group_column not in df.columns:
        available = list(df.columns)
        log.error(
            "Group column '%s' not found in metadata. Available columns: %s",
            group_column, available,
        )
        sys.exit(1)

    return df[group_column].dropna()


# ---------------------------------------------------------------------------
# PERMANOVA / PERMDISP (manual fallbacks; approximate, used only when no
# --stats-tsv is supplied). Prefer the QZV-derived stats TSV for reporting.
# ---------------------------------------------------------------------------

def compute_permanova(
    dm: pd.DataFrame,
    group_series: pd.Series,
    n_perm: int = 999,
    seed: int = 42,
) -> Tuple[float, float]:
    """Compute PERMANOVA F and p from a distance matrix (Anderson 2001)."""
    common = [s for s in dm.index if s in group_series.index
              and group_series[s] not in (None, "", float("nan"))]
    if len(common) < 3:
        return float("nan"), float("nan")

    groups = group_series.loc[common]
    mat = dm.loc[common, common].values.astype(float)
    unique_groups = sorted(groups.unique())
    g = np.array([unique_groups.index(x) for x in groups])

    def f_stat(dist_mat, grp):
        n = len(grp)
        total_ss = np.sum(dist_mat ** 2) / n
        within_ss = sum(
            np.sum(dist_mat[np.ix_(np.where(grp == gi)[0],
                                   np.where(grp == gi)[0])] ** 2)
            / max(len(np.where(grp == gi)[0]), 1)
            for gi in np.unique(grp)
        )
        between_ss = total_ss - within_ss
        a = len(np.unique(grp))
        return (between_ss / (a - 1)) / (within_ss / (n - a))

    obs_f = f_stat(mat, g)
    np.random.seed(seed)
    perm_fs = [f_stat(mat, np.random.permutation(g)) for _ in range(n_perm)]
    p = (np.sum(np.array(perm_fs) >= obs_f) + 1) / (n_perm + 1)
    return obs_f, p


def compute_permdisp(
    dm: pd.DataFrame,
    group_series: pd.Series,
    n_perm: int = 999,
    seed: int = 42,
) -> Tuple[float, float]:
    """Approximate PERMDISP F and p (homogeneity of dispersion)."""
    common = [s for s in dm.index if s in group_series.index
              and group_series[s] not in (None, "", float("nan"))]
    if len(common) < 3:
        return float("nan"), float("nan")

    groups = group_series.loc[common]
    mat = dm.loc[common, common].values.astype(float)
    unique_groups = sorted(groups.unique())

    def group_dispersions(dist_mat, grp_labels):
        dispersions = []
        for g in unique_groups:
            idxs = [i for i, s in enumerate(common) if grp_labels.iloc[i] == g]
            if len(idxs) < 2:
                continue
            sub = dist_mat[np.ix_(idxs, idxs)]
            mean_dists = sub.mean(axis=1)
            dispersions.extend(mean_dists.tolist())
        return dispersions

    def f_from_dispersions(disp, grp_labels):
        vals = np.array(disp)
        overall_mean = vals.mean()
        between = sum(
            len([i for i, s in enumerate(common) if grp_labels.iloc[i] == g])
            * (np.mean([vals[i] for i, s in enumerate(common) if grp_labels.iloc[i] == g]) - overall_mean) ** 2
            for g in unique_groups
        ) / (len(unique_groups) - 1)
        within = np.var(vals, ddof=len(unique_groups))
        return between / within if within > 0 else float("nan")

    obs_disp = group_dispersions(mat, groups)
    obs_f = f_from_dispersions(obs_disp, groups)

    np.random.seed(seed)
    perm_fs = []
    for _ in range(n_perm):
        perm_groups = groups.copy()
        perm_groups.values[:] = np.random.permutation(groups.values)
        pd_disp = group_dispersions(mat, perm_groups)
        perm_fs.append(f_from_dispersions(pd_disp, perm_groups))

    p = (np.sum(np.array(perm_fs) >= obs_f) + 1) / (n_perm + 1)
    return obs_f, p


# ---------------------------------------------------------------------------
# Stats TSV loading
# ---------------------------------------------------------------------------

def load_stats_tsv(tsv_path: Path) -> Dict[str, Dict]:
    """
    Load a beta_stats_summary_*.tsv.
    Returns {metric_label: {group_column_value: {method: {F, p}}}}.

    Keeping the group_column value is essential: a single metric (e.g.
    Bray-Curtis) can have one row per grouping variable (batch, state),
    and the figure must pick the row matching the grouping it is drawing.
    """
    stats: Dict[str, Dict] = {}
    METRIC_NORM = {
        "unweighted unifrac":  "Unweighted UniFrac",
        "unweighted_unifrac":  "Unweighted UniFrac",
        "weighted unifrac":    "Weighted UniFrac",
        "weighted_unifrac":    "Weighted UniFrac",
        "bray-curtis":         "Bray-Curtis",
        "bray_curtis":         "Bray-Curtis",
        "jaccard":             "Jaccard",
    }
    with tsv_path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            metric_raw = row.get("metric", "").strip().lower()
            metric = METRIC_NORM.get(metric_raw, row.get("metric", "").strip())
            method = row.get("method", "").strip().lower()
            # If method column is blank, infer from statistic_name:
            # pseudo-F -> permanova, F-value -> permdisp
            if not method:
                stat_name = row.get("statistic_name", "").strip().lower()
                if "pseudo" in stat_name:
                    method = "permanova"
                elif "f-value" in stat_name or "f_value" in stat_name:
                    method = "permdisp"
            try:
                f_val = float(row.get("statistic", "nan"))
                p_val = float(row.get("p_value", "nan"))
            except ValueError:
                continue
            group_val = row.get("group_column", "").strip()
            stats.setdefault(metric, {}).setdefault(group_val, {})[method] = {"F": f_val, "p": p_val}
    return stats


def _select_group_stats(metric_entry, group_col_name):
    """Pick the stats entry whose group_column matches the active grouping.
    Tolerant of label differences, e.g. 'batch' matches 'batch_short'."""
    if not metric_entry:
        return None
    gcl = (group_col_name or "").strip().lower()
    if gcl:
        for g, v in metric_entry.items():          # exact match first
            if g.strip().lower() == gcl:
                return v
        for g, v in metric_entry.items():          # then substring either way
            gl = g.strip().lower()
            if gl and (gl in gcl or gcl in gl):
                return v
    if len(metric_entry) == 1:                     # single-group TSV: unambiguous
        return next(iter(metric_entry.values()))
    return None


# ---------------------------------------------------------------------------
# Confidence ellipse helper
# ---------------------------------------------------------------------------

def _confidence_ellipse(x, y, ax, n_std=2.448, **kwargs):
    """
    Draw a covariance-based confidence ellipse around a set of 2D points.
    n_std = sqrt(chi2.ppf(0.95, df=2)) ~= 2.448 gives a true 95% ellipse
    for bivariate normal data (n_std=2.0 would give ~86%).
    """
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as transforms
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    if np.any(np.isnan(cov)) or np.linalg.det(cov) == 0:
        return
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell = Ellipse(
        (0, 0),
        width=np.sqrt(1 + pearson) * 2,
        height=np.sqrt(1 - pearson) * 2,
        **kwargs,
    )
    transf = (
        transforms.Affine2D()
        .rotate_deg(45)
        .scale(np.sqrt(cov[0, 0]) * n_std, np.sqrt(cov[1, 1]) * n_std)
        .translate(np.mean(x), np.mean(y))
    )
    ell.set_transform(transf + ax.transData)
    ax.add_patch(ell)


# ---------------------------------------------------------------------------
# Alpha diversity plot
# ---------------------------------------------------------------------------

def plot_alpha(
    alpha_qzas: List[Path],
    group_series: pd.Series,
    group_order: List[str],
    palette: Dict,
    marker: str,
    outpath: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_panels = len(alpha_qzas)
    ncols = min(n_panels, 2)
    nrows = (n_panels + 1) // 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5.5 * nrows))
    if n_panels == 1:
        axes = [[axes]]
    elif nrows == 1:
        axes = [axes]
    axes_flat = [ax for row in axes for ax in (row if hasattr(row, "__iter__") else [row])]

    fig.suptitle(
        f"Alpha Diversity - {marker}\n"
        f"Grouped by: {group_series.name or 'group'}",
        fontsize=13, fontweight="bold", y=0.99,
    )

    np.random.seed(0)
    positions = list(range(1, len(group_order) + 1))
    group_ns = {g: int((group_series == g).sum()) for g in group_order}

    for idx, (qza_path, ax) in enumerate(zip(alpha_qzas, axes_flat)):
        metric_name = infer_metric_name(qza_path)
        try:
            vector = read_alpha_vector(qza_path)
        except Exception as e:
            log.warning("Could not read %s: %s", qza_path.name, e)
            ax.set_visible(False)
            continue

        group_vals = [
            vector.loc[[s for s in vector.index if s in group_series.index and group_series[s] == g]].values
            for g in group_order
        ]

        valid = [(g, v) for g, v in zip(group_order, group_vals) if len(v) >= 2]
        if len(valid) == 2:
            _, p_val = scipy_stats.mannwhitneyu(valid[0][1], valid[1][1], alternative="two-sided")
            p_label = f"Mann-Whitney U  p = {p_val:.3f}" + ("" if p_val >= 0.05 else " *")
        elif len(valid) > 2:
            _, p_val = scipy_stats.kruskal(*[v for _, v in valid])
            p_label = f"Kruskal-Wallis  p = {p_val:.3f}" + ("" if p_val >= 0.05 else " *")
        else:
            p_val, p_label = None, ""

        box_positions = [pos for pos, v in zip(positions, group_vals) if len(v) >= 2]
        box_data = [v for v in group_vals if len(v) >= 2]
        box_groups = [g for g, v in zip(group_order, group_vals) if len(v) >= 2]

        if box_data:
            bp = ax.boxplot(
                box_data, positions=box_positions, widths=0.42,
                patch_artist=True, notch=False,
                medianprops=dict(color="black", linewidth=2.2),
                whiskerprops=dict(color="#444", linewidth=1.1),
                capprops=dict(color="#444", linewidth=1.1),
                flierprops=dict(marker="", linestyle="none"),
                boxprops=dict(linewidth=1.1),
            )
            for patch, g in zip(bp["boxes"], box_groups):
                color = palette["colors"][group_order.index(g) % len(palette["colors"])]
                patch.set_facecolor(color)
                patch.set_alpha(0.5)
                patch.set_edgecolor(color)

        for pos, g, vals in zip(positions, group_order, group_vals):
            if len(vals) == 0:
                continue
            color = palette["colors"][group_order.index(g) % len(palette["colors"])]
            marker_sym = palette["markers"][group_order.index(g) % len(palette["markers"])]
            jitter = np.random.uniform(-0.13, 0.13, size=len(vals))
            ax.scatter(pos + jitter, vals, marker=marker_sym, color=color,
                       s=55, zorder=5, edgecolors="white", linewidths=0.6, alpha=0.95)

        ax.set_xticks(positions)
        ax.set_xticklabels(
            [f"{g}\n(n={group_ns.get(g, 0)})" for g in group_order], fontsize=10
        )
        ax.set_ylabel(metric_name, fontsize=11)
        ax.set_title(metric_name, fontsize=12, fontweight="bold")
        ax.set_xlim(0.4, len(group_order) + 0.6)
        ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if p_label:
            is_sig = p_val is not None and p_val <= 0.05
            y_top = ax.get_ylim()[1]
            col = palette["sig_color"] if is_sig else "#555"
            ax.text(np.mean(positions), y_top * 0.97, p_label,
                    ha="center", va="top", fontsize=9, style="italic", color=col)

    for ax in axes_flat[n_panels:]:
        ax.set_visible(False)

    legend_handles = [
        plt.scatter([], [], marker=palette["markers"][i % len(palette["markers"])],
                    color=palette["colors"][i % len(palette["colors"])], s=60,
                    edgecolors="white", linewidths=0.6, label=g)
        for i, g in enumerate(group_order)
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=len(group_order),
               fontsize=10, frameon=True, bbox_to_anchor=(0.5, 0.005),
               edgecolor="#ccc")

    plt.tight_layout(rect=[0, 0.04, 1, 0.98])
    try:
        fig.savefig(outpath, dpi=180, bbox_inches="tight", facecolor="white")
    except OSError as exc:
        raise OSError(
            f"Could not write alpha figure to {outpath}.\n"
            f"  Original error: {exc}"
        ) from exc
    plt.close()
    log.info("Alpha diversity figure saved: %s", outpath)


# ---------------------------------------------------------------------------
# Beta diversity plot
# ---------------------------------------------------------------------------

def plot_beta(
    beta_pcoa_qzas: List[Path],
    beta_dm_qzas: Optional[List[Path]],
    group_series: pd.Series,
    group_order: List[str],
    palette: Dict,
    marker: str,
    outpath: Path,
    group_column_name: str = "",
    stats_tsv: Optional[Path] = None,
    n_perm: int = 999,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_panels = len(beta_pcoa_qzas)
    ncols = min(n_panels, 2)
    nrows = (n_panels + 1) // 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 6 * nrows))
    if n_panels == 1:
        axes = [[axes]]
    elif nrows == 1:
        axes = [axes]
    axes_flat = [ax for row in axes for ax in (row if hasattr(row, "__iter__") else [row])]

    fig.suptitle(
        f"Beta Diversity - {marker}\n"
        f"Grouped by: {group_series.name or 'group'}",
        fontsize=13, fontweight="bold", y=0.96,
    )

    precomp_stats: Dict = {}
    if stats_tsv is not None and stats_tsv.exists():
        precomp_stats = load_stats_tsv(stats_tsv)
        log.info("Loaded pre-computed stats from: %s", stats_tsv)

    dm_lookup: Dict[str, pd.DataFrame] = {}
    if beta_dm_qzas:
        for dm_qza in beta_dm_qzas:
            try:
                dm = read_distance_matrix(dm_qza)
                dm_lookup[infer_metric_name(dm_qza)] = dm
            except Exception as e:
                log.warning("Could not read DM %s: %s", dm_qza.name, e)

    plotted_ids: set = set()  # sample IDs actually drawn, for accurate legend n's

    for idx, (pcoa_qza, ax) in enumerate(zip(beta_pcoa_qzas, axes_flat)):
        metric_name = infer_metric_name(pcoa_qza)

        try:
            prop_exp, coords = read_pcoa_ordination(pcoa_qza)
        except Exception as e:
            log.warning("Could not read PCoA %s: %s", pcoa_qza.name, e)
            ax.set_visible(False)
            continue

        plotted_ids.update(coords.keys())

        # Ellipses
        for gi, g in enumerate(group_order):
            color = palette["colors"][gi % len(palette["colors"])]
            xs = np.array([coords[s][0] for s in coords if group_series.get(s) == g])
            ys = np.array([coords[s][1] for s in coords if group_series.get(s) == g])
            if len(xs) >= 3:
                _confidence_ellipse(xs, ys, ax, facecolor=color, alpha=0.13, edgecolor="none")
                _confidence_ellipse(xs, ys, ax, facecolor="none", edgecolor=color,
                                    linewidth=1.0, alpha=0.65)

        # Points
        for gi, g in enumerate(group_order):
            color = palette["colors"][gi % len(palette["colors"])]
            mk = palette["markers"][gi % len(palette["markers"])]
            xs = np.array([coords[s][0] for s in coords if group_series.get(s) == g])
            ys = np.array([coords[s][1] for s in coords if group_series.get(s) == g])
            if len(xs) > 0:
                ax.scatter(xs, ys, marker=mk, color=color, s=65,
                           zorder=5, edgecolors="white", linewidths=0.6)

        pc1 = prop_exp[0] if len(prop_exp) > 0 else 0
        pc2 = prop_exp[1] if len(prop_exp) > 1 else 0
        ax.set_xlabel(f"PC1 ({pc1:.1f}%)", fontsize=11)
        ax.set_ylabel(f"PC2 ({pc2:.1f}%)", fontsize=11)
        ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.45)
        ax.axhline(0, color="grey", linewidth=0.3)
        ax.axvline(0, color="grey", linewidth=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Stats: pick the row matching the active grouping, else compute on the fly.
        perm_f, perm_p = float("nan"), float("nan")
        disp_f, disp_p = float("nan"), float("nan")

        ms = _select_group_stats(precomp_stats.get(metric_name), group_column_name)
        if ms is not None:
            if "permanova" in ms:
                perm_f, perm_p = ms["permanova"]["F"], ms["permanova"]["p"]
            if "permdisp" in ms:
                disp_f, disp_p = ms["permdisp"]["F"], ms["permdisp"]["p"]
        elif metric_name in dm_lookup:
            log.info("Computing PERMANOVA for %s...", metric_name)
            perm_f, perm_p = compute_permanova(dm_lookup[metric_name], group_series, n_perm)
            log.info("Computing PERMDISP for %s...", metric_name)
            try:
                disp_f, disp_p = compute_permdisp(dm_lookup[metric_name], group_series, n_perm)
            except Exception as e:
                log.warning("PERMDISP failed for %s (will show NaN): %s", metric_name, e)
                disp_f, disp_p = float("nan"), float("nan")

        perm_sig = not np.isnan(perm_p) and perm_p <= 0.05

        def fmt(val):
            return f"{val:.3f}" if not np.isnan(val) else "n/a"

        # PERMANOVA goes in the panel title so it can never overlap points.
        # PERMDISP is omitted from the figure; report it in the caption if needed.
        title_col = palette["sig_color"] if perm_sig else "black"
        if not np.isnan(perm_f):
            stat_line = f"PERMANOVA pseudo-F={perm_f:.2f}, p={fmt(perm_p)}" + (" *" if perm_sig else "")
            title_str = f"{metric_name}\n{stat_line}"
        else:
            title_str = metric_name
        ax.set_title(title_str, fontsize=12, fontweight="bold", color=title_col, pad=6)

    for ax in axes_flat[n_panels:]:
        ax.set_visible(False)

    # Legend: group colors only (stats now live in each panel title).
    # Count n's from samples actually plotted, not the full metadata, so the
    # per-group n's sum to what is on screen (and to the PERMANOVA n).
    group_ns = {
        g: sum(1 for s in plotted_ids if group_series.get(s) == g)
        for g in group_order
    }
    scatter_handles = [
        plt.scatter([], [], marker=palette["markers"][i % len(palette["markers"])],
                    color=palette["colors"][i % len(palette["colors"])], s=65,
                    edgecolors="white", linewidths=0.6,
                    label=f"{g} (n={group_ns.get(g, 0)})")
        for i, g in enumerate(group_order)
    ]
    fig.legend(handles=scatter_handles, loc="lower center",
               ncol=4, fontsize=9,
               frameon=True, bbox_to_anchor=(0.5, 0.0), edgecolor="#ccc")

    plt.tight_layout(rect=[0, 0.12, 1, 0.99])
    try:
        fig.savefig(outpath, dpi=180, bbox_inches="tight", facecolor="white")
    except OSError as exc:
        raise OSError(
            f"Could not write beta figure to {outpath}.\n"
            f"  Original error: {exc}"
        ) from exc
    plt.close()
    log.info("Beta diversity figure saved: %s", outpath)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Generate alpha and beta diversity figures from QIIME 2 QZA files.\n"
            "No QIIME 2 installation required; reads QZA files directly."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    req = p.add_argument_group("required arguments")
    req.add_argument("--metadata", required=True, type=Path,
                     help="QIIME 2 metadata TSV. First column must be 'sample-id'.")
    req.add_argument("--group-column", required=True,
                     help="Metadata column to color and group samples by.")

    alpha = p.add_argument_group("alpha diversity inputs")
    alpha.add_argument("--alpha-qza", nargs="*", type=Path, default=None,
                       help="One or more alpha diversity vector QZA files.")

    beta = p.add_argument_group("beta diversity inputs")
    beta.add_argument("--beta-pcoa", nargs="*", type=Path, default=None,
                      help="One or more PCoA ordination QZA files (*_pcoa_results.qza).")
    beta.add_argument("--beta-dm", nargs="*", type=Path, default=None,
                      help="Distance matrix QZA files (for on-the-fly stats). "
                           "Optional if --stats-tsv is provided.")
    beta.add_argument("--stats-tsv", type=Path, default=None,
                      help="Pre-computed beta stats summary TSV. If provided, "
                           "stats are read from this file instead of recomputed.")
    beta.add_argument("--n-perm", type=int, default=999,
                      help="Permutations for on-the-fly PERMANOVA/PERMDISP. Default: 999.")

    opt = p.add_argument_group("options")
    opt.add_argument("--group-order", nargs="*", default=None,
                     help="Explicit group order for plots. Default: alphabetical.")
    opt.add_argument("--palette", choices=list(PALETTES.keys()), default="purple",
                     help="Color palette (default: purple). wong = colorblind-safe.")
    opt.add_argument("--group-label", default=None,
                     help="Display label for the grouping variable in figure titles. "
                          "Defaults to --group-column value.")
    opt.add_argument("--marker", default="",
                     help="Marker label for figure titles and filenames (e.g., 16S).")
    opt.add_argument("--outdir", type=Path, default=Path("results/figures"),
                     help="Output directory. Default: results/figures/")
    opt.add_argument("--suffix", default="",
                     help="Optional extra suffix for output filenames.")

    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    if args.alpha_qza is None and args.beta_pcoa is None:
        log.error("Provide at least one of --alpha-qza or --beta-pcoa.")
        return 1

    metadata_path = args.metadata.resolve()
    if not metadata_path.exists():
        log.error("Metadata not found: %s", metadata_path)
        return 1

    group_series = load_metadata(metadata_path, args.group_column)
    group_series.name = args.group_label if args.group_label else args.group_column

    if args.group_order:
        group_order = args.group_order
        missing = [g for g in group_order if g not in group_series.values]
        if missing:
            log.warning("Groups in --group-order not found in metadata: %s", missing)
    else:
        group_order = sorted(group_series.unique())

    log.info("Groups: %s", group_order)
    for g in group_order:
        n = int((group_series == g).sum())
        log.info("  %s: n=%d", g, n)

    palette = PALETTES[args.palette]
    marker_label = args.marker or ""
    suffix = f"_{args.marker}" if args.marker else ""
    suffix += f"_{args.suffix}" if args.suffix else ""
    suffix += f"_{args.palette}"

    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.alpha_qza:
        log.info("Generating alpha diversity figure (%d metrics)...", len(args.alpha_qza))
        alpha_out = args.outdir / f"alpha_diversity{suffix}.png"
        plot_alpha(
            alpha_qzas=[p.resolve() for p in args.alpha_qza],
            group_series=group_series,
            group_order=group_order,
            palette=palette,
            marker=marker_label,
            outpath=alpha_out,
        )

    if args.beta_pcoa:
        log.info("Generating beta diversity figure (%d metrics)...", len(args.beta_pcoa))
        beta_out = args.outdir / f"beta_diversity{suffix}.png"
        plot_beta(
            beta_pcoa_qzas=[p.resolve() for p in args.beta_pcoa],
            beta_dm_qzas=[p.resolve() for p in args.beta_dm] if args.beta_dm else None,
            group_series=group_series,
            group_order=group_order,
            palette=palette,
            marker=marker_label,
            outpath=beta_out,
            group_column_name=args.group_column,
            stats_tsv=args.stats_tsv.resolve() if args.stats_tsv else None,
            n_perm=args.n_perm,
        )

    log.info("Done. Outputs in: %s", args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
