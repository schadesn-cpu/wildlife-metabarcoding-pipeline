#!/usr/bin/env python3
"""
04_rarefaction.py
=================
Standalone rarefaction curve generator for QIIME 2 feature tables.

Run this BEFORE committing to a sampling depth in the main pipeline.
It generates an alpha-rarefaction QZV (viewable at view.qiime2.org) and
a static PNG figure with guidance annotations to help you choose a threshold.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
HOW TO CHOOSE A RAREFACTION DEPTH — what to look for
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

A rarefaction (subsampling) curve plots alpha diversity (e.g., observed
features, Faith's PD, Shannon entropy) against sequencing depth. Good
sequencing produces curves that rise steeply then flatten — indicating
that additional reads do not discover new diversity.

STEP 1 — Find the plateau ("knee point")
  Look for the depth at which most samples' curves level off and become
  approximately horizontal. This is where you have captured the true
  diversity of the community and more reads add little information.

STEP 2 — Balance depth vs. sample retention
  Any sample with fewer reads than the chosen depth is DROPPED from all
  downstream analyses. You must decide how many samples you can afford
  to lose:
    • Losing 0–1 samples: generally acceptable if they are low-quality
    • Losing >20% of samples: reconsider lowering the depth
    • Losing samples unevenly across groups: a serious bias risk

  A common rule of thumb is to choose the depth at the 10th percentile
  of sample read counts (keeping ~90% of samples), provided curves have
  plateaued by that point.

STEP 3 — Check that curves have plateaued at your chosen depth
  If curves are still rising steeply at the chosen depth, your diversity
  estimates will be artificially low. Consider deeper sequencing or
  acknowledge this limitation.

STEP 4 — Consistency across groups
  Look at curves colored by your grouping variable. If one group
  consistently has lower sequencing depth than the other, a high
  threshold could disproportionately drop samples from that group.

GOOD DEPTH: curves plateau well before it, ≥80–90% of samples retained,
            both groups represented evenly.
BAD DEPTH:  curves still rising at that depth, or >20% of samples lost.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Usage:
  # Basic (requires rooted tree for Faith's PD)
  python 04_rarefaction.py \\
      --table    qiime2/dada2/table_16S.qza \\
      --tree     qiime2/diversity/rooted-tree_16S.qza \\
      --metadata metadata/qiime/metadata_16S.tsv \\
      --max-depth 20000 \\
      --outdir   results/rarefaction/

  # With group coloring and candidate depth marked
  python 04_rarefaction.py \\
      --table        qiime2/dada2/table_16S.qza \\
      --tree         qiime2/diversity/rooted-tree_16S.qza \\
      --metadata     metadata/qiime/metadata_16S.tsv \\
      --max-depth    20000 \\
      --candidate-depth 8000 \\
      --group-column Group \\
      --marker       16S \\
      --outdir       results/rarefaction/

  # Without tree (skips Faith's PD)
  python 04_rarefaction.py \\
      --table    qiime2/dada2/table_16S.qza \\
      --metadata metadata/qiime/metadata_16S.tsv \\
      --max-depth 20000 \\
      --no-tree \\
      --outdir   results/rarefaction/

Dependencies:
  Python stdlib only for core logic; matplotlib + numpy for static figure.
  QIIME 2 conda environment must be active (provides qiime CLI).

  pip install matplotlib numpy   (or conda install -c conda-forge matplotlib numpy)

Outputs:
  <outdir>/alpha_rarefaction_<marker>.qzv   — Interactive QIIME viz
  <outdir>/rarefaction_summary_<marker>.png — Static annotated figure
  <outdir>/rarefaction_depth_report_<marker>.tsv — Per-sample read counts + drop status
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import shlex
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def run_cmd(cmd: List[str], dry_run: bool = False) -> None:
    cmd_str = " ".join(shlex.quote(c) for c in cmd)
    log.info("$ %s", cmd_str)
    if dry_run:
        log.info("[DRY RUN] Command not executed.")
        return
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        log.error("Command failed (exit %d): %s", result.returncode, cmd_str)
        sys.exit(result.returncode)


def require_exists(path: Path, label: str) -> None:
    if not path.exists():
        log.error("%s not found: %s", label, path)
        sys.exit(1)


def which_or_die(exe: str) -> None:
    import shutil
    if not shutil.which(exe):
        log.error(
            "'%s' not found on PATH. Activate your QIIME 2 conda environment first.\n"
            "  e.g.: conda activate qiime2-amplicon-2024.5", exe
        )
        sys.exit(1)


# ---------------------------------------------------------------------------
# Read count extraction (from feature-table summarize QZV)
# ---------------------------------------------------------------------------

def get_sample_read_counts(table_qza: Path, dry_run: bool) -> Dict[str, int]:
    """
    Run qiime feature-table summarize and parse per-sample counts from the QZV.
    Returns {sample_id: read_count}.
    """
    import tempfile
    if dry_run:
        log.info("[DRY RUN] Skipping read count extraction.")
        return {}

    with tempfile.TemporaryDirectory() as td:
        summary_qzv = Path(td) / "summary.qzv"
        run_cmd([
            "qiime", "feature-table", "summarize",
            "--i-table", str(table_qza),
            "--o-visualization", str(summary_qzv),
        ], dry_run=False)

        counts: Dict[str, int] = {}
        with zipfile.ZipFile(summary_qzv) as zf:
            csv_name = next(
                (n for n in zf.namelist() if n.endswith("sample-frequency-detail.csv")),
                None,
            )
            if csv_name is None:
                log.warning("sample-frequency-detail.csv not found in summary QZV.")
                return {}
            with zf.open(csv_name) as fh:
                reader = csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8"))
                for row in reader:
                    vals = list(row.values())
                    sid = row.get("index") or vals[0]
                    freq = row.get("0") or vals[1]
                    counts[sid] = int(float(freq))
    return counts


# ---------------------------------------------------------------------------
# Depth recommendation
# ---------------------------------------------------------------------------

def recommend_depth(
    counts: Dict[str, int],
    candidate_depth: Optional[int] = None,
) -> Tuple[int, Dict]:
    """
    Suggest a sampling depth and report retention statistics.
    Returns (recommended_depth, stats_dict).
    """
    if not counts:
        return candidate_depth or 0, {}

    sorted_counts = sorted(counts.values())
    n = len(sorted_counts)

    stats = {
        "n_samples": n,
        "min_reads": sorted_counts[0],
        "max_reads": sorted_counts[-1],
        "median_reads": sorted_counts[n // 2],
        "p10_reads": sorted_counts[max(0, int(n * 0.10))],  # 10th percentile
        "p25_reads": sorted_counts[max(0, int(n * 0.25))],
    }

    # Recommended depth = 10th percentile (retains ~90% of samples) as a starting point
    recommended = stats["p10_reads"]

    if candidate_depth is not None:
        retained = sum(1 for v in sorted_counts if v >= candidate_depth)
        dropped  = n - retained
        stats["candidate_depth"]          = candidate_depth
        stats["candidate_retained"]       = retained
        stats["candidate_dropped"]        = dropped
        stats["candidate_pct_retained"]   = 100.0 * retained / n

    stats["recommended_depth"] = recommended
    return recommended, stats


def print_depth_report(counts: Dict[str, int], stats: Dict, candidate_depth: Optional[int]) -> None:
    """Print a formatted depth selection report to the terminal."""
    if not counts:
        return

    sep = "─" * 70
    print(f"\n{sep}")
    print("  RAREFACTION DEPTH SELECTION REPORT")
    print(sep)
    print(f"  Total samples:   {stats['n_samples']}")
    print(f"  Min reads:       {stats['min_reads']:,}")
    print(f"  Median reads:    {stats['median_reads']:,}")
    print(f"  Max reads:       {stats['max_reads']:,}")
    print(f"  10th percentile: {stats['p10_reads']:,}  ← retains ~90% of samples")
    print(f"  25th percentile: {stats['p25_reads']:,}  ← retains ~75% of samples")
    print()

    if candidate_depth is not None:
        pct = stats.get("candidate_pct_retained", 0)
        retained = stats.get("candidate_retained", 0)
        dropped = stats.get("candidate_dropped", 0)
        flag = "✓" if pct >= 80 else ("⚠" if pct >= 50 else "✗")
        print(f"  Candidate depth: {candidate_depth:,}")
        print(f"  Samples retained: {retained}/{stats['n_samples']} ({pct:.0f}%)  {flag}")
        if dropped > 0:
            print(f"  Samples DROPPED ({dropped}):")
            for sid, c in sorted(counts.items(), key=lambda x: x[1]):
                if c < candidate_depth:
                    print(f"    {sid:<45} {c:,} reads")
        print()
        if pct < 50:
            print("  ⚠ WARNING: >50% of samples will be dropped at this depth.")
            print(f"    Consider lowering to ~{stats['p10_reads']:,} (10th percentile).")
        elif pct < 80:
            print("  ⚠ CAUTION: >20% of samples will be dropped at this depth.")
            print(f"    Consider lowering to ~{stats['p10_reads']:,} (10th percentile).")
        else:
            print("  ✓ Retention looks acceptable (≥80% of samples retained).")

    print()
    print("  HOW TO READ THE RAREFACTION CURVES:")
    print("    1. Look for where curves FLATTEN (plateau) — that is your minimum depth.")
    print("    2. Choose a depth AFTER the plateau where you retain ≥80–90% of samples.")
    print("    3. Curves still rising steeply = undersampled; consider deeper sequencing.")
    print("    4. Compare groups: uneven curve heights may indicate a group-depth confound.")
    print(f"{sep}\n")


def write_depth_tsv(counts: Dict[str, int], candidate_depth: Optional[int], out_path: Path) -> None:
    """Write per-sample read count + drop status to TSV."""
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "read_count", "dropped_at_candidate_depth"])
        for sid, cnt in sorted(counts.items(), key=lambda x: x[1]):
            drop = ""
            if candidate_depth is not None:
                drop = "yes" if cnt < candidate_depth else "no"
            writer.writerow([sid, cnt, drop])
    log.info("Depth report written: %s", out_path)


# ---------------------------------------------------------------------------
# Static PNG figure
# ---------------------------------------------------------------------------

def make_rarefaction_figure(
    counts: Dict[str, int],
    stats: Dict,
    candidate_depth: Optional[int],
    marker: str,
    out_path: Path,
) -> None:
    """
    Generate an annotated static figure showing the read depth distribution
    and where the candidate/recommended threshold falls.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        log.warning("matplotlib/numpy not available — skipping static figure. "
                    "Install with: pip install matplotlib numpy")
        return

    if not counts:
        return

    sorted_counts = sorted(counts.values())
    n = len(sorted_counts)
    sample_ids = [k for k, _ in sorted(counts.items(), key=lambda x: x[1])]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        f"Rarefaction Depth Selection — {marker}\n"
        "Use this to choose --sampling-depth for core-metrics-phylogenetic",
        fontsize=13, fontweight="bold", y=0.99
    )

    # ── Left panel: sorted read counts per sample ─────────────────────────
    ax = axes[0]
    colors = ["#7B2D8B"] * n
    if candidate_depth is not None:
        colors = ["#C19FD8" if c < candidate_depth else "#7B2D8B" for c in sorted_counts]

    ax.barh(range(n), sorted_counts, color=colors, height=0.75)

    if candidate_depth is not None:
        ax.axvline(candidate_depth, color="#B22222", linewidth=1.8, linestyle="--",
                   label=f"Candidate depth: {candidate_depth:,}")
    ax.axvline(stats["p10_reads"], color="#E59866", linewidth=1.4, linestyle=":",
               label=f"10th pct: {stats['p10_reads']:,}")
    ax.axvline(stats["median_reads"], color="#2E86C1", linewidth=1.4, linestyle=":",
               label=f"Median: {stats['median_reads']:,}")

    ax.set_yticks(range(n))
    ax.set_yticklabels(sample_ids, fontsize=6.5)
    ax.set_xlabel("Read count", fontsize=11)
    ax.set_title("Per-sample read counts\n(purple = retained, lavender = dropped)", fontsize=11)
    ax.legend(fontsize=8.5, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Right panel: cumulative retention curve ───────────────────────────
    ax2 = axes[1]
    depths = np.linspace(0, max(sorted_counts), 500)
    pct_retained = [100 * sum(1 for c in sorted_counts if c >= d) / n for d in depths]

    ax2.plot(depths, pct_retained, color="#7B2D8B", linewidth=2.2)
    ax2.axhline(90, color="#E59866", linewidth=1.2, linestyle=":", label="90% retention")
    ax2.axhline(80, color="#B8860B", linewidth=1.2, linestyle=":", label="80% retention")

    if candidate_depth is not None:
        pct = stats.get("candidate_pct_retained", 0)
        ax2.axvline(candidate_depth, color="#B22222", linewidth=1.8, linestyle="--",
                    label=f"Candidate: {candidate_depth:,} ({pct:.0f}% retained)")
        ax2.scatter([candidate_depth], [pct], color="#B22222", s=80, zorder=5)

    ax2.axvline(stats["p10_reads"], color="#E59866", linewidth=1.2, linestyle=":")

    ax2.set_xlabel("Sampling depth", fontsize=11)
    ax2.set_ylabel("% samples retained", fontsize=11)
    ax2.set_title("Sample retention vs. rarefaction depth\n"
                  "Choose a depth where curves have plateaued AND retention ≥80%",
                  fontsize=11)
    ax2.set_ylim(0, 105)
    ax2.legend(fontsize=8.5)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.5)

    # Guidance annotation
    fig.text(
        0.5, 0.01,
        "After viewing this plot: open the .qzv in view.qiime2.org to inspect individual curves. "
        "Look for the 'knee' where most samples plateau, then use --sampling-depth in 03_run_full_metabarcoding_pipeline.py.",
        ha="center", fontsize=8.5, color="#555", style="italic"
    )

    plt.tight_layout(rect=[0, 0.04, 1, 0.97])
    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close()
    log.info("Static figure saved: %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Generate QIIME 2 alpha-rarefaction curves and a static depth-selection figure.\n\n"
            "Run this BEFORE the main pipeline to determine an appropriate --sampling-depth.\n"
            "Open the output .qzv at view.qiime2.org to inspect per-sample curves interactively."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    req = p.add_argument_group("required arguments")
    req.add_argument("--table",    required=True, type=Path, help="Feature table QZA.")
    req.add_argument("--metadata", required=True, type=Path, help="QIIME 2 metadata TSV.")

    tree = p.add_argument_group("tree arguments")
    tree.add_argument("--tree",    default=None, type=Path,
                      help="Rooted phylogenetic tree QZA (required for Faith's PD). "
                           "Omit with --no-tree to skip phylogenetic metrics.")
    tree.add_argument("--no-tree", action="store_true",
                      help="Skip Faith's PD (no rooted tree available).")

    opt = p.add_argument_group("options")
    opt.add_argument("--max-depth",       type=int, default=20000,
                     help="Maximum rarefaction depth for the curve. Default: 20000.")
    opt.add_argument("--steps",           type=int, default=20,
                     help="Number of depth steps in rarefaction curve. Default: 20.")
    opt.add_argument("--iterations",      type=int, default=10,
                     help="Rarefaction iterations per step. Default: 10.")
    opt.add_argument("--candidate-depth", type=int, default=None,
                     help="Candidate sampling depth to evaluate and mark on the figure.")
    opt.add_argument("--group-column",    default=None,
                     help="Metadata column to color samples by in the QZV (optional).")
    opt.add_argument("--marker",          default="",
                     help="Marker name for output file naming (e.g., 16S).")
    opt.add_argument("--outdir",          type=Path, default=Path("results/rarefaction"),
                     help="Output directory. Default: results/rarefaction/")
    opt.add_argument("--dry-run",         action="store_true",
                     help="Print commands without executing.")

    return p.parse_args()


def main() -> int:
    args = parse_args()

    which_or_die("qiime")

    table = args.table.resolve()
    metadata = args.metadata.resolve()
    require_exists(table, "Feature table")
    require_exists(metadata, "Metadata file")

    if not args.no_tree:
        if args.tree is None:
            log.error(
                "No --tree provided. Provide a rooted tree QZA for Faith's PD, "
                "or pass --no-tree to skip phylogenetic metrics."
            )
            return 1
        tree = args.tree.resolve()
        require_exists(tree, "Rooted tree")
    else:
        tree = None

    args.outdir.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.marker}" if args.marker else ""

    # ── Step 1: Get per-sample read counts and print report ───────────────
    log.info("Extracting per-sample read counts...")
    counts = get_sample_read_counts(table, args.dry_run)
    _, stats = recommend_depth(counts, args.candidate_depth)
    print_depth_report(counts, stats, args.candidate_depth)

    # Write depth TSV
    depth_tsv = args.outdir / f"rarefaction_depth_report{suffix}.tsv"
    if not args.dry_run and counts:
        write_depth_tsv(counts, args.candidate_depth, depth_tsv)

    # ── Step 2: Run alpha-rarefaction ─────────────────────────────────────
    qzv_out = args.outdir / f"alpha_rarefaction{suffix}.qzv"
    cmd = [
        "qiime", "diversity", "alpha-rarefaction",
        "--i-table", str(table),
        "--m-metadata-file", str(metadata),
        "--p-max-depth", str(args.max_depth),
        "--p-steps", str(args.steps),
        "--p-iterations", str(args.iterations),
        "--o-visualization", str(qzv_out),
    ]
    if tree is not None:
        cmd += ["--i-phylogeny", str(tree)]

    run_cmd(cmd, dry_run=args.dry_run)

    # ── Step 3: Static annotated figure ──────────────────────────────────
    fig_out = args.outdir / f"rarefaction_summary{suffix}.png"
    if not args.dry_run:
        make_rarefaction_figure(counts, stats, args.candidate_depth, args.marker or "Unknown", fig_out)

    # ── Summary ───────────────────────────────────────────────────────────
    print("\n── Outputs ─────────────────────────────────────────────────────────")
    print(f"  Interactive QZV:  {qzv_out}")
    print(f"    → Open at: https://view.qiime2.org")
    print(f"  Static figure:    {fig_out}")
    print(f"  Depth report TSV: {depth_tsv}")
    if counts and "p10_reads" in stats:
        print(f"\n  Suggested starting depth: {stats['p10_reads']:,} reads")
        print(f"  (10th percentile — retains ~90% of samples)")
        print(f"  Confirm by inspecting curves in the QZV before proceeding.\n")

    print(
        "Next step — run diversity analysis with 03_run_full_metabarcoding_pipeline.py:\n"
        "  python 03_run_full_metabarcoding_pipeline.py diversity \\\n"
        "    --project-dir . \\\n"
        "    --marker <MARKER> \\\n"
        "    --rarefaction-depth <CHOSEN_DEPTH>"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
