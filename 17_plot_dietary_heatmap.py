#!/usr/bin/env python3
"""
17_plot_dietary_heatmap.py
==========================
Multi-marker dietary detection heatmap — species × season.

Rows = 35 prey species in the seasonal cohort, grouped by habitat:
    Freshwater → Anadromous → Estuarine → Marine →
    Marine/Anadromous → Marine/Estuary

Columns = 3 ecological seasons:
    Breeding | Freshwater_Nonbreeding | Saltwater

Cells encode per-bird detection frequency within season cohort
(0 = white, 1 = full saturation). Cell text shows raw counts as
"n / N_cohort" for unambiguous reading. Right-side symbols mark the
marker that recovered each species.

Inputs (from 16_pooled_dietary_diversity.py outputs):
    taxonomy_by_season.tsv  species × season counts + frequencies

Usage:
    python scripts/17_plot_dietary_heatmap.py \\
        --indir  results/multimarker/pooled/ \\
        --outdir results/multimarker/pooled/figures/
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Wong palette + sequential colormap
# ─────────────────────────────────────────────────────────────────────────────
WONG = {
    "blue":       "#0072B2",
    "orange":     "#E69F00",
    "green":      "#009E73",
    "sky":        "#56B4E9",
    "vermillion": "#D55E00",
    "purple":     "#CC79A7",
    "yellow":     "#F0E442",
    "black":      "#000000",
    "grey":       "#999999",
}

FREQ_CMAP = LinearSegmentedColormap.from_list(
    "wong_blues", ["#FFFFFF", "#C5DCED", WONG["blue"]],
)

# ─────────────────────────────────────────────────────────────────────────────
# Layout constants
# ─────────────────────────────────────────────────────────────────────────────
HABITAT_ORDER = [
    "Freshwater",
    "Anadromous",
    "Estuarine",
    "Marine",
    "Marine/Anadromous",
    "Marine/Estuary",
]

# Light pastel row-band backgrounds. Marine is neutral pale grey so the
# frequency-colored (blue) cells stand out against it.
HABITAT_BG = {
    "Freshwater":         "#E5F4ED",
    "Anadromous":         "#FDF1D9",
    "Estuarine":          "#E8F4FA",
    "Marine":             "#EFEFEF",
    "Marine/Anadromous":  "#FAE3D6",
    "Marine/Estuary":     "#EFE3F0",
}

SEASON_ORDER = ["Breeding", "Freshwater_Nonbreeding", "Saltwater"]
SEASON_HEADERS = {
    "Breeding":               "Breeding\n(May–Aug)",
    "Freshwater_Nonbreeding": "FW Non-Breeding\n(Apr + Sep)",
    "Saltwater":              "Saltwater\n(Oct–Mar)",
}

MARKER_ORDER   = ["MiFish", "cytb", "18S_invert"]
MARKER_SYMBOLS = {"MiFish": "o", "cytb": "^", "18S_invert": "s"}
MARKER_COLORS  = {
    "MiFish":     WONG["blue"],
    "cytb":       WONG["vermillion"],
    "18S_invert": WONG["green"],
}
MARKER_LABELS  = {
    "MiFish":     "MiFish 12S",
    "cytb":       "Cytochrome b",
    "18S_invert": "18S V9 (inverts)",
}


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def italicize(name: str) -> str:
    """Render species name in italics via matplotlib mathtext."""
    if name.endswith(" sp."):
        return r"$\mathit{" + name[:-4] + "}$ sp."
    if " " in name:
        return r"$\mathit{" + name.replace(" ", r"\ ") + "}$"
    return r"$\mathit{" + name + "}$"


def load_taxonomy(indir: Path) -> pd.DataFrame:
    tax = pd.read_csv(indir / "taxonomy_by_season.tsv", sep="\t")
    log.info("Loaded taxonomy_by_season: %d species × seasons", len(tax))
    return tax


def order_species(tax: pd.DataFrame) -> pd.DataFrame:
    """Sort by habitat (defined order), then total_n desc, then name."""
    tax = tax.copy()
    tax["habitat_rank"] = tax["habitat"].map(
        {h: i for i, h in enumerate(HABITAT_ORDER)}
    ).fillna(99).astype(int)
    return tax.sort_values(
        ["habitat_rank", "total_n", "species"],
        ascending=[True, False, True],
    ).reset_index(drop=True)


def _strip_off(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# ─────────────────────────────────────────────────────────────────────────────
# Plot components
# ─────────────────────────────────────────────────────────────────────────────
def plot_main_heatmap(ax, tax: pd.DataFrame):
    species  = tax["species"].tolist()
    habitats = tax["habitat"].tolist()
    n_species = len(species)
    n_seasons = len(SEASON_ORDER)

    freq = np.zeros((n_species, n_seasons))
    n    = np.zeros((n_species, n_seasons), dtype=int)
    N    = np.zeros((n_species, n_seasons), dtype=int)
    for i, sp_row in tax.iterrows():
        for j, season in enumerate(SEASON_ORDER):
            freq[i, j] = sp_row[f"{season}_freq"]
            n[i, j]    = int(sp_row[f"{season}_n"])
            N[i, j]    = int(sp_row[f"{season}_n_cohort"])

    # Habitat row-band backgrounds
    for habitat in HABITAT_ORDER:
        rows = [i for i, h in enumerate(habitats) if h == habitat]
        if not rows:
            continue
        ax.add_patch(Rectangle(
            (0, min(rows)), n_seasons, max(rows) - min(rows) + 1,
            facecolor=HABITAT_BG.get(habitat, "#FFFFFF"),
            edgecolor="none", zorder=0,
        ))

    # Frequency-colored cells with n/N text overlay
    for i in range(n_species):
        for j in range(n_seasons):
            if freq[i, j] > 0:
                color = FREQ_CMAP(freq[i, j])
                ax.add_patch(Rectangle(
                    (j + 0.02, i + 0.02), 0.96, 0.96,
                    facecolor=color, edgecolor="none", zorder=2,
                ))
                text_color = "white" if freq[i, j] >= 0.35 else WONG["black"]
                ax.text(
                    j + 0.5, i + 0.5, f"{n[i, j]}/{N[i, j]}",
                    ha="center", va="center", fontsize=8.5,
                    color=text_color, fontweight="bold", zorder=3,
                )

    # Subtle grid
    for x in range(n_seasons + 1):
        ax.axvline(x, color="#CCCCCC", linewidth=0.4, zorder=1)
    for y in range(n_species + 1):
        ax.axhline(y, color="#DDDDDD", linewidth=0.3, zorder=1)

    # Bold horizontal dividers between habitat groups
    for i in range(1, n_species):
        if habitats[i] != habitats[i - 1]:
            ax.axhline(i, color=WONG["black"], linewidth=1.2, zorder=3)

    for s in ax.spines.values():
        s.set_linewidth(0.7)

    # Column headers (season + cohort size, multi-line)
    cohort_sizes = {s: int(tax[f"{s}_n_cohort"].iloc[0]) for s in SEASON_ORDER}
    headers = [
        f"{SEASON_HEADERS[s]}\nn = {cohort_sizes[s]}"
        for s in SEASON_ORDER
    ]
    ax.set_xticks([j + 0.5 for j in range(n_seasons)])
    ax.set_xticklabels(headers, fontsize=10)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")

    ax.set_yticks([i + 0.5 for i in range(n_species)])
    ax.set_yticklabels([italicize(s) for s in species], fontsize=9)

    ax.set_xlim(0, n_seasons)
    ax.set_ylim(n_species, 0)


def plot_marker_strip(ax, tax: pd.DataFrame):
    """Right-side single-column strip: marker symbol per species."""
    n_species = len(tax)
    for i, sp_row in tax.iterrows():
        marker_name = sp_row["marker"]
        ax.plot(
            0.5, i + 0.5,
            marker=MARKER_SYMBOLS.get(marker_name, "?"),
            color=MARKER_COLORS.get(marker_name, WONG["grey"]),
            markerfacecolor=MARKER_COLORS.get(marker_name, WONG["grey"]),
            markeredgecolor=MARKER_COLORS.get(marker_name, WONG["grey"]),
            markersize=9, linestyle="None",
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(n_species, 0)
    _strip_off(ax)


def plot_marker_legend(ax):
    ax.axis("off")
    handles = [
        mlines.Line2D(
            [], [], marker=MARKER_SYMBOLS[m], linestyle="None",
            markerfacecolor=MARKER_COLORS[m],
            markeredgecolor=MARKER_COLORS[m],
            markersize=9, label=MARKER_LABELS[m],
        )
        for m in MARKER_ORDER
    ]
    ax.legend(
        handles=handles, title="Marker",
        loc="center left", bbox_to_anchor=(0.0, 0.5),
        frameon=False, fontsize=9, title_fontsize=10,
        handletextpad=0.6, labelspacing=0.4,
    )


def plot_habitat_legend(ax):
    ax.axis("off")
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=HABITAT_BG[h],
                   edgecolor=WONG["grey"], linewidth=0.6, label=h)
        for h in HABITAT_ORDER
    ]
    ax.legend(
        handles=handles, title="Habitat (row band)",
        loc="center left", bbox_to_anchor=(0.0, 0.5),
        frameon=False, fontsize=9, title_fontsize=10,
        handletextpad=0.6, labelspacing=0.4,
    )


def plot_colorbar(ax):
    """Compact horizontal colorbar with caption below."""
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    ax.imshow(gradient, aspect="auto", cmap=FREQ_CMAP, extent=[0, 1, 0, 1])
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1.0"], fontsize=8)
    ax.set_yticks([])
    ax.set_title("Detection frequency", fontsize=9, pad=4)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)


# ─────────────────────────────────────────────────────────────────────────────
# Figure assembly
# ─────────────────────────────────────────────────────────────────────────────
def plot_figure(tax: pd.DataFrame, outdir: Path):
    n_species = len(tax)

    # Wider cells so the longest column header ("FW Non-Breeding") fits cleanly
    cell_w = 1.4
    cell_h = 0.28
    fig_w  = cell_w * 3 + 5.4
    fig_h  = cell_h * n_species + 2.4

    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")

    # Two top-level cols: (heatmap + marker strip) | (legends + colorbar stack)
    gs = GridSpec(
        nrows=1, ncols=3,
        figure=fig,
        width_ratios=[3 * cell_w, 0.45, 2.6],
        wspace=0.08,
        left=0.20, right=0.97, top=0.86, bottom=0.06,
    )
    ax_main         = fig.add_subplot(gs[0, 0])
    ax_marker_strip = fig.add_subplot(gs[0, 1])

    # Right column: Marker legend → colorbar → Habitat legend (vertical stack)
    gs_right = GridSpecFromSubplotSpec(
        nrows=3, ncols=1,
        subplot_spec=gs[0, 2],
        height_ratios=[1.0, 0.7, 2.0],
        hspace=0.8,
    )
    ax_marker_legend  = fig.add_subplot(gs_right[0])
    ax_cbar           = fig.add_subplot(gs_right[1])
    ax_habitat_legend = fig.add_subplot(gs_right[2])

    plot_main_heatmap(ax_main, tax)
    plot_marker_strip(ax_marker_strip, tax)
    plot_marker_legend(ax_marker_legend)
    plot_colorbar(ax_cbar)
    plot_habitat_legend(ax_habitat_legend)

    fig.suptitle(
        "Multi-marker dietary detections by ecological season\n"
        f"({n_species} prey species; MiFish 12S + cytochrome b + 18S V9 invertebrates)",
        fontsize=12, y=0.985,
    )

    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        out = outdir / f"dietary_detection_heatmap.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        log.info("Saved: %s", out)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--indir", type=Path, required=True,
                    help="Directory containing taxonomy_by_season.tsv")
    p.add_argument("--outdir", type=Path, required=True,
                    help="Directory to write the figure files")
    args = p.parse_args()

    tax = load_taxonomy(args.indir)
    tax = order_species(tax)
    log.info("Species in heatmap: %d (sorted by habitat → total_n)", len(tax))

    plot_figure(tax, args.outdir)
    log.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
