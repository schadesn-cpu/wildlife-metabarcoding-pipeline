#!/usr/bin/env python3
"""
11_plot_summary_figures.py
==========================
Per-sample summary figures for the tick project, generated from the
BLAST-verified, collapsed sample table.

  fig1 - Blood meal host composition (cytb)
  fig2 - Pathogen / parasite prevalence (16S + 18S)
  fig3 - Pathogen co-detection (pooled) vs co-infection (single-tick)
  fig5 - Pathogen-host co-occurrence heatmap

Outputs both annotated (lab) and manuscript (clean) versions, mirroring
the convention in run_all_figures.py.

All paths, column names, category mappings, and group orders are read
from pipeline_config.yml via config_loader. No hardcoded project-specific
values remain in this script.

Usage (run from project root):
    python scripts/11_plot_summary_figures.py [--dry-run]
    python scripts/11_plot_summary_figures.py --figures fig1 fig2
    python scripts/11_plot_summary_figures.py --palette wong
    python scripts/11_plot_summary_figures.py --config /path/to/pipeline_config.yml

Expected pipeline_config.yml additions (under a 'tick' section):

    tick:
      sample_table:        results/cooccurrence/adrianna_sample_table_collapsed.tsv
      single_tick_table:   results/cooccurrence/single_tick_cooccurrence.tsv
      host_column:         cytb_host
      pathogen_columns:    [Anaplasma, Borreliella, Babesia_odocoilei, Bartonella,
                            Babesia_microti, Borrelia, Francisella, Babesia_sp]
      pathogen_categories:
        human:    [Anaplasma, Borreliella, Babesia_microti, Borrelia]
        wildlife: [Babesia_odocoilei, Bartonella, Babesia_sp]
        endo:     [Francisella]
      host_order:          [Deer, Moose, Bear, Deer_ked, Porcupine, Opossum]
      heatmap_pathogens:   [Anaplasma, Borreliella, Babesia_odocoilei, Bartonella,
                            Babesia_microti, Borrelia, Francisella]

Author: Samantha Schade - MEED Lab, UNH - 2026-05-01
"""

import argparse
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# ---------------------------------------------------------------------------
# Config loader auto-discovery (matches run_all_figures.py)
# ---------------------------------------------------------------------------

_PIPELINE_DIR = Path(__file__).resolve().parent
for _search_dir in [_PIPELINE_DIR, _PIPELINE_DIR.parent]:
    if (_search_dir / "config_loader.py").exists():
        if str(_search_dir) not in sys.path:
            sys.path.insert(0, str(_search_dir))
        break

import config_loader as _cl  # noqa: E402


# ---------------------------------------------------------------------------
# Palettes — colorblind-safe categorical + matched sequential
# ---------------------------------------------------------------------------

# Bang Wong, Nature Methods 8:441 (2011). Same as Okabe-Ito.
WONG = {
    "black":      "#000000",
    "orange":     "#E69F00",
    "sky":        "#56B4E9",
    "green":      "#009E73",
    "yellow":     "#F0E442",
    "blue":       "#0072B2",
    "vermillion": "#D55E00",
    "purple":     "#CC79A7",
    "neutral":    "#9A9A9A",
}


def get_palette(name: str) -> dict:
    """Return a categorical palette dict. Extendable for other named palettes."""
    if name in ("wong", "okabe", "okabe_ito"):
        return WONG
    raise ValueError(
        f"Unknown palette: {name!r}. Known: 'wong'."
    )


def get_sequential_cmap(name: str):
    """Return a sequential colormap that pairs visually with the categorical palette."""
    if name in ("wong", "okabe", "okabe_ito"):
        return LinearSegmentedColormap.from_list(
            "wong_seq",
            ["#FFFFFF", WONG["yellow"], WONG["orange"], WONG["vermillion"]],
            N=256,
        )
    # cividis is also colorblind-safe; fallback for unknown palettes
    return mpl.colormaps["cividis"]


# ---------------------------------------------------------------------------
# Global rcParams
# ---------------------------------------------------------------------------

mpl.rcParams.update({
    "font.size":          11,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.titlesize":     13,
    "axes.labelsize":     11,
})


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def missing_path(path: Path, label: str) -> bool:
    """Return True (and warn) if a required input file is missing."""
    if not path.exists():
        print(f"  Skip {label} — missing: {path}")
        return True
    return False


def missing_columns(df: pd.DataFrame, cols: list, label: str) -> bool:
    """Return True (and warn) if any required columns are missing from df."""
    bad = [c for c in cols if c not in df.columns]
    if bad:
        print(f"  Skip {label} — sample table missing columns: {bad}")
        return True
    return False


def make_dirs(ann: Path, ms: Path, sub: str = "summary"):
    """Create and return the (annotated, manuscript) subdirectories."""
    a = ann / "tick" / sub
    m = ms  / "tick" / sub
    a.mkdir(parents=True, exist_ok=True)
    m.mkdir(parents=True, exist_ok=True)
    return a, m


def save_both(fig, ann_path: Path, ms_path: Path,
              dry_run: bool, label: str) -> None:
    """
    Save the annotated version (with title), then strip the title and
    save the manuscript version. Mirrors the both() pattern in
    run_all_figures.py but for direct matplotlib figures.
    """
    if dry_run:
        print(f"  [DRY] would write:")
        print(f"        {ann_path}")
        print(f"        {ms_path}")
        plt.close(fig)
        return

    fig.savefig(ann_path, dpi=200, bbox_inches="tight")
    fig.savefig(ann_path.with_suffix(".pdf"), bbox_inches="tight")

    # manuscript: strip title from every axis
    for ax in fig.axes:
        ax.set_title("")
    fig.savefig(ms_path, dpi=200, bbox_inches="tight")
    fig.savefig(ms_path.with_suffix(".pdf"), bbox_inches="tight")

    plt.close(fig)
    print(f"  Done {label}")


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def _tick_cfg(cfg) -> dict:
    """Pull the 'tick' subsection from pipeline_config.yml, with defaults."""
    sub = getattr(cfg, "tick", None) or cfg.raw.get("tick", {})
    if not sub:
        raise ValueError(
            "pipeline_config.yml has no 'tick:' section. See module docstring."
        )
    return sub


def load_sample_table(cfg) -> pd.DataFrame:
    """Load the collapsed per-sample table (one row per biological sample)."""
    tick = _tick_cfg(cfg)
    p = Path(cfg.resolve(tick.get(
        "sample_table",
        "results/cooccurrence/adrianna_sample_table_collapsed.tsv",
    )))
    if missing_path(p, "sample table"):
        return None
    return pd.read_csv(p, sep="\t")


def load_single_tick_table(cfg) -> pd.DataFrame:
    """Load the single-tick co-occurrence subset (n=124 in current data)."""
    tick = _tick_cfg(cfg)
    p = Path(cfg.resolve(tick.get(
        "single_tick_table",
        "results/cooccurrence/single_tick_cooccurrence.tsv",
    )))
    if missing_path(p, "single-tick table"):
        return None
    return pd.read_csv(p, sep="\t")


# ---------------------------------------------------------------------------
# Per-figure generators
# ---------------------------------------------------------------------------

def fig1_host_composition(cfg, ann: Path, ms: Path,
                          dry_run: bool, palette: dict) -> None:
    """fig1 — Blood meal host composition from cytb."""
    df = load_sample_table(cfg)
    if df is None:
        return
    tick = _tick_cfg(cfg)
    host_col = tick.get("host_column", "cytb_host")
    if missing_columns(df, [host_col], "fig1 host composition"):
        return

    counts = df[host_col].dropna().value_counts()
    n = int(counts.sum())
    pct = counts / n * 100

    fig, ax = plt.subplots(figsize=(9, 5.5))
    y = np.arange(len(counts))[::-1]
    ax.barh(y, pct.values, color=palette["blue"],
            edgecolor="white", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(counts.index)
    ax.set_xlabel("% of samples")
    ax.set_title(f"Blood meal host composition (cytb, n={n})")
    ax.set_xlim(0, 100)
    for yi, p, c in zip(y, pct.values, counts.values):
        ax.text(p + 1.2, yi, f"{p:.1f}% (n={int(c)})",
                va="center", fontsize=10)
    fig.tight_layout()

    save_both(fig,
              ann / "fig1_host_composition_annotated.png",
              ms  / "fig1_host_composition.png",
              dry_run, "fig1 host composition")


def fig2_pathogen_prevalence(cfg, ann: Path, ms: Path,
                             dry_run: bool, palette: dict) -> None:
    """fig2 — Pathogen / parasite prevalence, colored by category."""
    df = load_sample_table(cfg)
    if df is None:
        return
    tick = _tick_cfg(cfg)
    cols = tick.get("pathogen_columns", [])
    if not cols:
        print("  Skip fig2 — no pathogen_columns in config")
        return
    if missing_columns(df, cols, "fig2 pathogen prevalence"):
        return

    cats = tick.get("pathogen_categories", {})  # human / wildlife / endo
    cat_color = {
        "human":    palette["vermillion"],
        "wildlife": palette["green"],
        "endo":     palette["neutral"],
    }
    cat_label = {
        "human":    "Human-relevant",
        "wildlife": "Wildlife / non-human",
        "endo":     "Endosymbiont",
    }
    name_to_cat = {name: cat for cat, names in cats.items() for name in names}

    n = len(df)
    rows = []
    for c in cols:
        # treat numeric > 0 OR truthy strings as positive
        col_pos = pd.to_numeric(df[c], errors="coerce").fillna(0).gt(0)
        if col_pos.sum() == 0:
            col_pos = df[c].astype(str).str.lower().isin(
                {"true", "1", "yes", "positive"}
            )
        pos = int(col_pos.sum())
        rows.append({
            "pathogen": c,
            "n":        pos,
            "pct":      pos / n * 100 if n else 0.0,
            "cat":      name_to_cat.get(c, "endo"),
        })
    rows.sort(key=lambda r: r["pct"], reverse=True)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    y = np.arange(len(rows))[::-1]
    colors = [cat_color[r["cat"]] for r in rows]
    ax.barh(y, [r["pct"] for r in rows], color=colors,
            edgecolor="white", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([r["pathogen"] for r in rows])
    ax.set_xlabel("% of samples positive")
    ax.set_title(f"Pathogen / parasite prevalence (n={n})")
    xmax = max((r["pct"] for r in rows), default=1) * 1.18
    ax.set_xlim(0, xmax)
    for yi, r in zip(y, rows):
        ax.text(r["pct"] + xmax * 0.01, yi,
                f"{r['pct']:.1f}% (n={r['n']})",
                va="center", fontsize=10)

    used_cats = sorted({r["cat"] for r in rows},
                       key=["human", "wildlife", "endo"].index)
    handles = [plt.Rectangle((0, 0), 1, 1, color=cat_color[c],
                             label=cat_label[c]) for c in used_cats]
    ax.legend(handles=handles, loc="lower right",
              frameon=False, fontsize=9)
    fig.tight_layout()

    save_both(fig,
              ann / "fig2_pathogen_prevalence_annotated.png",
              ms  / "fig2_pathogen_prevalence.png",
              dry_run, "fig2 pathogen prevalence")


def fig3_coinfection_rates(cfg, ann: Path, ms: Path,
                           dry_run: bool, palette: dict) -> None:
    """fig3 — Co-detection (all samples) vs co-infection (single-tick)."""
    df_all    = load_sample_table(cfg)
    df_single = load_single_tick_table(cfg)
    if df_all is None or df_single is None:
        return
    tick = _tick_cfg(cfg)
    cols = tick.get("pathogen_columns", [])
    if missing_columns(df_all,    cols, "fig3 (all)") or \
       missing_columns(df_single, cols, "fig3 (single)"):
        return

    def _bucket(df: pd.DataFrame) -> list:
        pos = pd.DataFrame({
            c: pd.to_numeric(df[c], errors="coerce").fillna(0).gt(0)
            for c in cols
        }).sum(axis=1)
        n = len(df)
        return [
            (pos == 0).sum() / n * 100,
            (pos == 1).sum() / n * 100,
            (pos == 2).sum() / n * 100,
            (pos >= 3).sum() / n * 100,
        ]

    cats = ["0 pathogens", "1 pathogen", "2 pathogens", "3+ pathogens"]
    all_pct    = _bucket(df_all)
    single_pct = _bucket(df_single)
    n_all      = len(df_all)
    n_single   = len(df_single)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(cats))
    w = 0.4
    b1 = ax.bar(x - w/2, all_pct,    w, color=palette["sky"],
                edgecolor="white", label=f"All samples (n={n_all})")
    b2 = ax.bar(x + w/2, single_pct, w, color=palette["blue"],
                edgecolor="white", label=f"Single-tick (n={n_single})")
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylabel("% of samples")
    ax.set_title("Pathogen co-detection (pooled) vs co-infection (single-tick)")
    ymax = max(all_pct + single_pct) * 1.15
    ax.set_ylim(0, ymax)
    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2,
                    h + ymax * 0.012,
                    f"{h:.1f}%", ha="center", fontsize=10)
    ax.legend(loc="upper right", frameon=True)
    fig.tight_layout()

    save_both(fig,
              ann / "fig3_coinfection_rates_annotated.png",
              ms  / "fig3_coinfection_rates.png",
              dry_run, "fig3 coinfection rates")


def fig5_cooccurrence_heatmap(cfg, ann: Path, ms: Path,
                              dry_run: bool, palette_name: str) -> None:
    """fig5 — Pathogen × host co-occurrence heatmap (% of pathogen-positives)."""
    df = load_sample_table(cfg)
    if df is None:
        return
    tick = _tick_cfg(cfg)
    host_col   = tick.get("host_column", "cytb_host")
    pathogens  = tick.get("heatmap_pathogens",
                          tick.get("pathogen_columns", []))
    host_order = tick.get("host_order", [])
    if not pathogens or not host_order:
        print("  Skip fig5 — heatmap_pathogens / host_order missing in config")
        return
    if missing_columns(df, pathogens + [host_col], "fig5 heatmap"):
        return

    # rows: pathogens; cols: hosts; value: % of pathogen-positive samples
    # whose host is this column.
    M = np.zeros((len(pathogens), len(host_order)))
    for i, p in enumerate(pathogens):
        pos_mask = pd.to_numeric(df[p], errors="coerce").fillna(0).gt(0)
        sub = df.loc[pos_mask, host_col].dropna()
        denom = len(sub)
        if denom == 0:
            continue
        for j, h in enumerate(host_order):
            M[i, j] = (sub == h).sum() / denom * 100

    cmap = get_sequential_cmap(palette_name)

    fig, ax = plt.subplots(figsize=(8.5, 6))
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(host_order)))
    ax.set_xticklabels(host_order, rotation=30, ha="right")
    ax.set_yticks(range(len(pathogens)))
    ax.set_yticklabels(pathogens)
    ax.set_title("Pathogen-host co-occurrence (%)")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if v <= 0:
                continue
            txt_color = "white" if v >= 60 else "black"
            ax.text(j, i, f"{v:.1f}",
                    ha="center", va="center",
                    color=txt_color, fontsize=10)
    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("% co-occurrence")
    fig.tight_layout()

    save_both(fig,
              ann / "fig5_cooccurrence_heatmap_annotated.png",
              ms  / "fig5_cooccurrence_heatmap.png",
              dry_run, "fig5 cooccurrence heatmap")


# ---------------------------------------------------------------------------
# Figure registry
# ---------------------------------------------------------------------------

# fig5 takes the palette name (it needs the matched sequential cmap),
# the others take the categorical palette dict.
GENS = {
    "fig1": ("dict", fig1_host_composition),
    "fig2": ("dict", fig2_pathogen_prevalence),
    "fig3": ("dict", fig3_coinfection_rates),
    "fig5": ("name", fig5_cooccurrence_heatmap),
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(
        description="Generate per-sample summary figures (annotated + manuscript) "
                    "for the tick project from the BLAST-verified sample table.",
    )
    p.add_argument("--config", default=None, metavar="YML",
                   help="Path to pipeline_config.yml. Default: auto-discovered.")
    p.add_argument("--dry-run", action="store_true",
                   help="Print what would be written without saving any files.")
    p.add_argument("--figures", nargs="+", default=list(GENS),
                   choices=list(GENS),
                   help="Subset of figures to generate. Default: all.")
    p.add_argument("--palette", default=None,
                   help="Override palette (e.g. 'wong'). "
                        "Default: figures.palette in config, fallback 'wong'.")
    args = p.parse_args()

    try:
        cfg = _cl.load_config(args.config)
    except (FileNotFoundError, ValueError) as e:
        print(f"\n[ERROR] {e}", file=sys.stderr)
        return 1

    pal_name = (args.palette
                or cfg.figures.get("palette", "wong"))
    try:
        palette = get_palette(pal_name)
    except ValueError as e:
        print(f"\n[ERROR] {e}", file=sys.stderr)
        return 1

    ann = cfg.resolve(cfg.figures.get("annotated_dir",
                                      "results/figures_annotated"))
    ms  = cfg.resolve(cfg.figures.get("manuscript_dir",
                                      "results/figures_manuscript"))
    a, m = make_dirs(ann, ms)

    print(f"Root:       {cfg.root}")
    print(f"Annotated:  {a}")
    print(f"Manuscript: {m}")
    print(f"Palette:    {pal_name}")
    print(f"Figures:    {args.figures}  |  Dry: {args.dry_run}")

    n_total = 0
    n_failed = 0
    for fig_id in args.figures:
        kind, fn = GENS[fig_id]
        print(f"\n>> {fig_id}")
        try:
            if kind == "dict":
                fn(cfg, a, m, args.dry_run, palette)
            else:  # 'name' — pass the palette name for cmap selection
                fn(cfg, a, m, args.dry_run, pal_name)
            n_total += 1
        except Exception as exc:
            print(f"\n  [UNEXPECTED ERROR in {fig_id}] "
                  f"{type(exc).__name__}: {exc}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            n_failed += 1
            n_total += 1

    print(f"\nDONE  |  ann: {a}  |  ms: {m}")
    print(f"  Total: {n_total}  |  Failed: {n_failed}")
    return 1 if n_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
