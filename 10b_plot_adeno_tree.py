#!/usr/bin/env python3
"""
10b_plot_adeno_tree.py
======================
Publication-quality phylogenetic tree for loon adenovirus OTUs.
Expanded reference set covering diverse Aviadenovirus species.

Tree topology and bootstrap support values were generated with IQ-TREE 2
using the GTR+G model (--model-test selected). UFBoot=96 for the loon OTU
clade provides strong statistical support for a novel aviadenovirus lineage.

Changes from previous version:
  - New tree: adeno_tree_expanded.nwk (14 taxa, GTR model)
  - No x/y axes — scale bar only
  - Asterisk (*) next to loon OTUs to flag putative novel species
  - Human adenovirus 2 as outgroup (Mastadenovirus — outgroup to all Aviadenovirus)
  - Bootstrap values displayed only for nodes with support >= 0.80
  - Two output versions: full (OTU1 + OTU2) and OTU1-only (--drop-otu2)

Input:
  results/adenovirus/adeno_tree_expanded.nwk

Outputs (written to results/adenovirus/):
  adeno_phylo_tree.png / .svg           — full tree (both OTUs)
  adeno_phylo_tree_otu1only.png / .svg  — OTU1 only (use --drop-otu2)

Usage:
  python scripts/10b_plot_adeno_tree.py
  python scripts/10b_plot_adeno_tree.py --drop-otu2
"""

from __future__ import annotations

import argparse
import logging
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from Bio import Phylo

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TREE_PATH  = "results/adenovirus/adeno_tree_expanded.nwk"
OUT_PNG    = "results/adenovirus/adeno_phylo_tree.png"
OUT_SVG    = "results/adenovirus/adeno_phylo_tree.svg"
OTU1_PNG   = "results/adenovirus/adeno_phylo_tree_otu1only.png"
OTU1_SVG   = "results/adenovirus/adeno_phylo_tree_otu1only.svg"

# Wong 2011 colorblind-safe palette (citable)
LOON_COLOR   = "#0072B2"   # blue  — loon OTUs (this study)
REF_COLOR    = "#444444"   # dark grey — reference aviadenoviruses
OUTGRP_COLOR = "#999999"   # light grey — outgroup (Mastadenovirus)
CONF_COLOR   = "#CC0000"   # red — bootstrap support values

FONT_FAMILY    = "Arial"
FONT_SIZE_TIP  = 9
FONT_SIZE_BOOT = 7
FIGURE_DPI     = 300

# Bootstrap threshold below which support values are suppressed.
# UFBoot values < 0.80 (80%) are considered insufficient to report.
BOOTSTRAP_THRESHOLD = 0.80

# Newick tip ID → display label mapping.
# Accession numbers and identity values are included for manuscript traceability.
LABEL_MAP = {
    # Loon OTUs — asterisk flags putative novel species pending formal description
    "adeno_OTU1_total9123":       "Loon aviadenovirus OTU1  (9,123 reads) *",
    "adeno_OTU6_total128":        "Loon aviadenovirus OTU2  (128 reads) *",

    # Closest BLAST hits — percent identity to loon OTU1
    "Aviadenovirus_YN06_76.8pct": "Aviadenovirus sp. YN06  (PP319115.1)  76.8% identity",
    "Aviadenovirus_YN10_76.6pct": "Aviadenovirus sp. YN10  (PP319121.1)  76.6% identity",

    # Broader Aviadenovirus reference set (NCBI)
    "Goose_adenovirus_4":         "Goose aviadenovirus A  (NC_017979.1)",
    "Goose_adenovirus_5":         "Goose adenovirus 5  (JQ178216.1)",
    "Duck_adenovirus_2":          "Duck aviadenovirus B  (KJ469653.1)",
    "Crane_adenovirus_1":         "Crane-associated adenovirus 1  (LC469780.1)",
    "Turkey_adenovirus_1":        "Turkey aviadenovirus B  (GU936707.2)",
    "Pigeon_adenovirus_1":        "Pigeon aviadenovirus A  (MW286325.1)",
    "Fowl_adenovirus_A":          "Fowl aviadenovirus A  (NC_001720.1)",
    "Psittacine_adenovirus_4":    "Psittacine aviadenovirus B  (KX577802.1)",

    # Outgroup — Mastadenovirus (mammalian), sister to all Aviadenovirus
    "Human_adenovirus_2_OUTGROUP": "Human adenovirus 2  (AC_000007.1)  [outgroup]",
}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Build and parse the CLI argument parser."""
    p = argparse.ArgumentParser(
        prog="10b_plot_adeno_tree.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--drop-otu2", action="store_true",
        help=(
            "Exclude OTU2 (adeno_OTU6_total128, 128 reads) from the tree. "
            "Use for the single-OTU figure variant."
        ),
    )
    p.add_argument(
        "--tree", default=TREE_PATH, metavar="NWK",
        help=f"Path to the Newick tree file. Default: {TREE_PATH}",
    )
    p.add_argument(
        "--outdir", default="results/adenovirus", metavar="DIR",
        help="Output directory for PNG and SVG files. Default: results/adenovirus/",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Core plotting function
# ---------------------------------------------------------------------------

def build_and_save_tree(tree_path: Path, drop_otu2: bool, outdir: Path) -> None:
    """
    Read the Newick tree, root on the outgroup, optionally prune OTU2,
    apply display labels, and save PNG + SVG figures.

    Parameters
    ----------
    tree_path : Path to the IQ-TREE Newick output file.
    drop_otu2 : If True, prune OTU2 before plotting.
    outdir    : Directory to write output files.
    """
    if not tree_path.exists():
        raise FileNotFoundError(
            f"Tree file not found: {tree_path}\n"
            f"  Run IQ-TREE on the adenovirus alignment first, or check --tree path."
        )

    # Parse Newick — BioPython Phylo handles IQ-TREE UFBoot format
    try:
        tree = Phylo.read(str(tree_path), "newick")
    except Exception as exc:
        raise ValueError(
            f"Could not parse Newick tree from {tree_path}: {exc}\n"
            f"  Ensure the file is a valid Newick format produced by IQ-TREE."
        ) from exc

    log.info("Tree loaded: %d terminals", tree.count_terminals())

    # Optionally prune OTU2 before rooting
    if drop_otu2:
        otu2_id = "adeno_OTU6_total128"
        try:
            tree.prune(otu2_id)
            log.info("OTU2 pruned (%s) — single-OTU version", otu2_id)
        except Exception as exc:
            log.warning(
                "Could not prune OTU2 (%s): %s — continuing without pruning.",
                otu2_id, exc,
            )

    # Root on human adenovirus outgroup BEFORE relabelling tips.
    # We search for "OUTGROUP" in the clade name to match the Newick label.
    outgroup = next(
        (c for c in tree.find_clades() if c.name and "OUTGROUP" in c.name),
        None,
    )
    if outgroup is None:
        raise ValueError(
            "Could not find the outgroup clade (expected a tip name containing 'OUTGROUP').\n"
            "  Check that the Newick file contains 'Human_adenovirus_2_OUTGROUP' or similar."
        )
    tree.root_with_outgroup(outgroup)
    log.info("Tree rooted on outgroup: %s", outgroup.name)

    # Apply display labels (Newick IDs → manuscript-ready labels)
    for clade in tree.find_clades():
        if clade.name and clade.name in LABEL_MAP:
            clade.name = LABEL_MAP[clade.name]

    # Suppress bootstrap values below threshold — reduces visual clutter
    # for weakly supported nodes that should not be reported in the manuscript.
    for clade in tree.find_clades():
        if clade.confidence is not None and clade.confidence < BOOTSTRAP_THRESHOLD:
            clade.confidence = None

    # ── Drawing ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 7.5))
    plt.subplots_adjust(top=0.93)
    plt.rcParams["font.family"] = FONT_FAMILY

    loon_labels = {v for v in LABEL_MAP.values() if "Loon" in v}
    all_labels  = set(LABEL_MAP.values())

    def label_colors(label: str) -> str:
        """Assign tip label color by taxon category."""
        if not label:
            return REF_COLOR
        if "Loon" in label:
            return LOON_COLOR
        if "outgroup" in label:
            return OUTGRP_COLOR
        return REF_COLOR

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda x: x.name if x.name else "",
        label_colors=label_colors,
    )

    # Style tip labels and bootstrap support values
    for text in ax.texts:
        t = text.get_text().strip()
        if t in loon_labels:
            text.set_fontsize(FONT_SIZE_TIP)
            text.set_color(LOON_COLOR)
            text.set_fontweight("bold")
        elif t in all_labels:
            is_outgroup = "outgroup" in t
            text.set_fontsize(FONT_SIZE_TIP)
            text.set_color(OUTGRP_COLOR if is_outgroup else REF_COLOR)
            text.set_fontstyle("italic" if is_outgroup else "normal")
        else:
            # Text not in label map — check if it's a bootstrap value
            try:
                val = float(t)
                if val >= BOOTSTRAP_THRESHOLD:
                    text.set_fontsize(FONT_SIZE_BOOT)
                    text.set_color(CONF_COLOR)
                    text.set_fontweight("bold")
                else:
                    text.set_text("")
            except ValueError:
                pass  # Not a number — leave as-is

    # ── Remove axes, add scale bar ────────────────────────────────────────────
    for attr in ("ylabel", "xlabel"):
        getattr(ax, f"set_{attr}")("")
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False)

    # Scale bar: 0.1 substitutions/site in the bottom-left corner.
    # Placed at 2% and 4% of axis range so it doesn't overlap with any tip label.
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    sb_len = 0.1
    sb_x   = xlim[0] + (xlim[1] - xlim[0]) * 0.02
    sb_y   = ylim[0] + (ylim[1] - ylim[0]) * 0.04
    ax.plot([sb_x, sb_x + sb_len], [sb_y, sb_y],
            color="black", linewidth=1.5, solid_capstyle="butt")
    ax.text(
        sb_x + sb_len / 2, sb_y - (ylim[1] - ylim[0]) * 0.025,
        "0.1 substitutions/site",
        ha="center", va="top", fontsize=7.5, color="black",
    )

    # ── Legend ────────────────────────────────────────────────────────────────
    loon_patch  = mpatches.Patch(color=LOON_COLOR,   label="Loon aviadenovirus OTUs (this study)")
    ref_patch   = mpatches.Patch(color=REF_COLOR,    label="Reference aviadenoviruses (NCBI)")
    out_patch   = mpatches.Patch(color=OUTGRP_COLOR, label="Outgroup (Mastadenovirus)")
    star_handle = mlines.Line2D(
        [], [], color="none", marker="*", markerfacecolor=LOON_COLOR,
        markersize=8, label="* Putative novel aviadenovirus species",
    )
    ax.legend(
        handles=[loon_patch, ref_patch, out_patch, star_handle],
        fontsize=8, framealpha=0.9, edgecolor="#CCCCCC",
        loc="upper right",
    )

    fig.tight_layout()

    # ── Save outputs ──────────────────────────────────────────────────────────
    outdir.mkdir(parents=True, exist_ok=True)

    stem = "adeno_phylo_tree_otu1only" if drop_otu2 else "adeno_phylo_tree"
    out_png = outdir / f"{stem}.png"
    out_svg = outdir / f"{stem}.svg"

    for out_path, fmt, kwargs in [
        (out_png, "png", {"dpi": FIGURE_DPI, "facecolor": "white"}),
        (out_svg, "svg", {"facecolor": "white"}),
    ]:
        try:
            fig.savefig(out_path, format=fmt, bbox_inches="tight", **kwargs)
            log.info("✓ Saved: %s", out_path)
        except OSError as exc:
            raise OSError(
                f"Could not write figure to {out_path}.\n"
                f"  Check that {outdir} exists and is writable.\n"
                f"  Original error: {exc}"
            ) from exc

    plt.close(fig)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    """
    Entry point for 10b_plot_adeno_tree.py.

    Parses arguments, reads and roots the Newick tree on the human adenovirus
    outgroup, applies publication labels, and saves PNG + SVG figures.
    Returns 0 on success, 1 on error.
    """
    args = parse_args()

    try:
        build_and_save_tree(
            tree_path=Path(args.tree),
            drop_otu2=args.drop_otu2,
            outdir=Path(args.outdir),
        )
        return 0

    except FileNotFoundError as exc:
        log.error("File not found:\n  %s", exc)
        return 1

    except ValueError as exc:
        log.error("Invalid input:\n  %s", exc)
        return 1

    except OSError as exc:
        log.error("File system error:\n  %s", exc)
        return 1

    except Exception as exc:
        log.error("Unexpected error: %s: %s", type(exc).__name__, exc)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
