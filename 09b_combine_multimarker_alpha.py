#!/usr/bin/env python3
"""
09b_combine_multimarker_alpha.py
=================================
Combine four single-marker Observed Features alpha diversity plots into one
publication-ready 4-panel figure (panels A–D, left to right).

This script runs after 09_plot_diversity.py has generated the per-marker
alpha diversity PNG files. It reads those PNGs and assembles them into a
single figure suitable for a methods comparison figure or supplementary panel.

Expected input files (default location: results/multimarker/figures/):
  multimarker_16S_observed.png
  multimarker_MiFish_observed.png
  multimarker_cytb_observed.png
  multimarker_18S_observed.png

ADAPT FOR YOUR STUDY:
  If your project uses different markers, update the 'panels' list below
  to match your marker names and file names. You can use 2, 3, or 4 panels —
  just adjust the list and the figsize accordingly.

Usage:
  python scripts/06b_combine_multimarker_alpha.py
  python scripts/06b_combine_multimarker_alpha.py --fig-dir results/myproject/figures
  python scripts/06b_combine_multimarker_alpha.py --outdir results/multimarker/figures
"""

import argparse
import sys
import pathlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# ── Panel configuration ───────────────────────────────────────────────────────
# ADAPT FOR YOUR STUDY: Change labels and filenames to match your markers.
# Each tuple is (panel label, filename relative to --fig-dir).
# Remove or add tuples for fewer/more panels. Update figsize in build_figure()
# if you change the number of panels (width = n_panels * 4.5).
DEFAULT_PANELS = [
    ("16S rRNA V4\n(Bacteria)",        "multimarker_16S_observed.png"),
    ("MiFish 12S\n(Fish prey)",         "multimarker_MiFish_observed.png"),
    ("Cytochrome b\n(Vertebrate prey)", "multimarker_cytb_observed.png"),
    ("18S rRNA V9\n(Eukaryotes)",       "multimarker_18S_observed.png"),
]


def build_figure(panels: list, fig_dir: pathlib.Path) -> plt.Figure:
    """
    Load each panel image from disk and assemble them into a multi-panel figure.

    Each panel gets a bold panel label (A, B, C, D) in the top-left corner
    and a subtitle below showing the marker name. The figure background is
    white for clean export to PDF or SVG.

    Args:
        panels:  List of (label, filename) tuples to display left-to-right.
        fig_dir: Directory containing the per-marker PNG files.

    Returns:
        A matplotlib Figure object ready to save.

    Raises:
        FileNotFoundError: If any input image file is missing.
        SystemExit: With a clear message listing the missing file path.
    """
    n = len(panels)

    # Check all input files exist before attempting to build the figure.
    # This gives a clean error message instead of a cryptic matplotlib traceback.
    missing = []
    for label, fname in panels:
        path = fig_dir / fname
        if not path.exists():
            missing.append(str(path))
    if missing:
        print("[ERROR] The following input images were not found:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)
        print("\nRun 09_plot_diversity.py first to generate per-marker figures,",
              file=sys.stderr)
        print("or check --fig-dir points to the correct directory.",
              file=sys.stderr)
        sys.exit(1)

    fig, axes = plt.subplots(1, n, figsize=(n * 4.5, 5))
    if n == 1:
        axes = [axes]  # ensure axes is always a list
    fig.patch.set_facecolor('white')

    for i, (ax, (label, fname)) in enumerate(zip(axes, panels)):
        img = mpimg.imread(fig_dir / fname)
        ax.imshow(img)
        ax.axis('off')
        ax.set_title(label, fontsize=13, fontweight='bold', pad=8,
                     fontfamily='Arial', color='#1a1a1a')
        # Bold panel letter (A, B, C, D) in top-left corner
        ax.text(-0.04, 1.04, chr(65 + i), transform=ax.transAxes,
                fontsize=16, fontweight='bold', va='top', ha='right',
                fontfamily='Arial', color='#1a1a1a')

    plt.tight_layout(w_pad=1.5)
    return fig


def main() -> int:
    """
    Parse arguments, build the 4-panel figure, and save PNG + SVG outputs.

    Returns 0 on success, 1 on error. Errors print a clear message to stderr
    rather than showing a Python traceback.
    """
    p = argparse.ArgumentParser(
        prog="06b_combine_multimarker_alpha.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--fig-dir",
        type=pathlib.Path,
        default=pathlib.Path("results/multimarker/figures"),
        help="Directory containing per-marker PNG files. Default: results/multimarker/figures/",
    )
    p.add_argument(
        "--outdir",
        type=pathlib.Path,
        default=None,
        help="Output directory. Default: same as --fig-dir.",
    )
    p.add_argument(
        "--output-stem",
        default="multimarker_alpha_observed_4panel",
        help="Output filename stem (without extension). Default: multimarker_alpha_observed_4panel",
    )
    args = p.parse_args()

    fig_dir = args.fig_dir
    out_dir = args.outdir or fig_dir

    if not fig_dir.exists():
        print(f"[ERROR] --fig-dir not found: {fig_dir}", file=sys.stderr)
        print("  Run 09_plot_diversity.py first to generate per-marker figures.",
              file=sys.stderr)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        fig = build_figure(DEFAULT_PANELS, fig_dir)
    except Exception as e:
        print(f"[ERROR] Failed to build figure: {e}", file=sys.stderr)
        return 1

    out_png = out_dir / f"{args.output_stem}_wong.png"
    out_svg = out_dir / f"{args.output_stem}_wong.svg"

    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_svg, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    print(f"✓ Saved: {out_png}")
    print(f"✓ Saved: {out_svg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
