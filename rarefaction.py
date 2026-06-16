#!/usr/bin/env python3
"""
rarefaction.py
==============
Step: rarefaction (optional `analyses.rarefaction` stage; runs before diversity)

Purpose:
    Generate QIIME 2 alpha-rarefaction curves and a static depth-selection
    figure so you can choose a sampling depth BEFORE committing to it in the
    diversity stage. Produces an interactive .qzv (view.qiime2.org) plus an
    annotated .png and a per-sample read-count TSV.

How to choose a depth (kept from the original — genuinely useful):
    A rarefaction curve plots alpha diversity vs sequencing depth; good data
    rises then plateaus. Pick the depth where most samples' curves flatten,
    balanced against sample retention — any sample below the depth is dropped
    from downstream analyses. A common starting point is the 10th percentile of
    per-sample read counts (keeps ~90% of samples), confirmed by eyeballing the
    curves. Watch for uneven dropout across groups, which biases comparisons.

Inputs:
    --table     FeatureTable[Frequency] .qza   (or derived from --marker)
    --metadata  QIIME 2 metadata .tsv          (or derived from --marker)
    --tree      rooted tree .qza for Faith's PD (optional; --no-tree to skip)
    pipeline_config.yml  analyses.rarefaction (max_depth, steps, iterations)

Outputs (in --outdir, default results/<marker>/all/rarefaction/):
    alpha_rarefaction[_<marker>].qzv     interactive curve
    rarefaction_summary[_<marker>].png   annotated depth-selection figure
    rarefaction_depth_report[_<marker>].tsv  per-sample counts + retention
    logs/run_manifest.jsonl              run appended on completion

Usage:
    python rarefaction.py --marker 16S
    python rarefaction.py --marker 16S --candidate-depth 8000
    python rarefaction.py --table t.qza --metadata m.tsv --no-tree --outdir out/

Requirements:
    QIIME 2 (provides the `qiime` CLI), matplotlib for the figure.
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

# --- make config_loader and the utils package importable regardless of cwd ---
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from config_loader import load_config, get_paths, get_metadata_path  # noqa: E402
from utils import checkpoint, provenance, validate                   # noqa: E402

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
    """Log and execute a shell command; print without running in dry-run mode."""
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
    """Log an error and exit if path does not exist."""
    if not path.exists():
        log.error("%s not found: %s", label, path)
        sys.exit(1)


def which_or_die(exe: str) -> None:
    """Log an error and exit if exe is not found on PATH."""
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

def parse_args(argv=None) -> argparse.Namespace:
    """Build and parse the argument parser for rarefaction.py."""
    p = argparse.ArgumentParser(
        description=(
            "Generate QIIME 2 alpha-rarefaction curves and a static depth-selection figure.\n\n"
            "Run this BEFORE the main pipeline to determine an appropriate --sampling-depth.\n"
            "Open the output .qzv at view.qiime2.org to inspect per-sample curves interactively."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    req = p.add_argument_group("inputs (derived from --marker if omitted)")
    req.add_argument("--marker",    default="", help="Marker name (e.g. 16S). Drives derived paths and config.")
    req.add_argument("--config",    default=None, help="Path to pipeline_config.yml.")
    req.add_argument("--table",    default=None, type=Path, help="Feature table QZA. Derived from --marker if omitted.")
    req.add_argument("--metadata", default=None, type=Path, help="QIIME 2 metadata TSV. Derived from --marker if omitted.")

    tree = p.add_argument_group("tree arguments")
    tree.add_argument("--tree",    default=None, type=Path,
                      help="Rooted phylogenetic tree QZA (required for Faith's PD). "
                           "Omit with --no-tree to skip phylogenetic metrics.")
    tree.add_argument("--no-tree", action="store_true",
                      help="Skip Faith's PD (no rooted tree available).")

    opt = p.add_argument_group("options (default to analyses.rarefaction in config)")
    opt.add_argument("--max-depth",       type=int, default=None,
                     help="Maximum rarefaction depth for the curve.")
    opt.add_argument("--steps",           type=int, default=None,
                     help="Number of depth steps in rarefaction curve.")
    opt.add_argument("--iterations",      type=int, default=None,
                     help="Rarefaction iterations per step.")
    opt.add_argument("--candidate-depth", type=int, default=None,
                     help="Candidate sampling depth to evaluate and mark on the figure.")
    opt.add_argument("--group-column",    default=None,
                     help="Metadata column to color samples by in the QZV (optional).")
    opt.add_argument("--outdir",          type=Path, default=None,
                     help="Output directory. Default: results/<marker>/all/rarefaction/")
    opt.add_argument("--dry-run",         action="store_true",
                     help="Print commands without executing.")

    return p.parse_args(argv)


def main(argv=None) -> int:
    """
    Parse arguments and run the rarefaction curve analysis.

    Resolves the table/metadata/tree/outdir and curve parameters from --marker
    and config when not given explicitly, validates inputs, generates the
    alpha-rarefaction QZV + depth-selection figure + per-sample TSV, then
    records the run. Returns 0 on success, 1 on error.
    """
    args = parse_args(argv)

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    rcfg = cfg.analyses.get("rarefaction", {})
    marker = args.marker

    # Derive inputs/outputs from the marker when not given explicitly.
    if not marker and (args.table is None or args.metadata is None):
        log.error("Provide --marker (to derive paths) or both --table and --metadata.")
        return 1
    table = (args.table or paths.engine_table_qza(marker, "all", nocontrols=True)).resolve()
    metadata = (args.metadata
                or get_metadata_path(cfg, marker, "all")).resolve()
    outdir = args.outdir or paths.engine_rarefaction_results_dir(marker, "all")

    # Curve parameters: CLI overrides config, config overrides built-in defaults.
    max_depth  = args.max_depth  if args.max_depth  is not None else rcfg.get("max_depth", 20000)
    steps      = args.steps      if args.steps      is not None else rcfg.get("steps", 20)
    iterations = args.iterations if args.iterations is not None else rcfg.get("iterations", 10)
    candidate  = (args.candidate_depth if args.candidate_depth is not None
                  else cfg.markers.get(marker, {}).get("rarefaction_depth") or None)

    validate.require_qiime()
    require_exists(table, "Feature table")
    require_exists(metadata, "Metadata file")

    # Tree: explicit --tree wins; else for phylo markers try the engine's rooted
    # tree; otherwise fall back to no-tree (skip Faith's PD) with a loud note.
    if args.no_tree:
        tree = None
    elif args.tree is not None:
        tree = args.tree.resolve()
        require_exists(tree, "Rooted tree")
    else:
        candidate_tree = paths.rooted_tree_qza(marker) if marker else None
        if candidate_tree and candidate_tree.exists():
            tree = candidate_tree
        else:
            log.warning("No rooted tree found for %s — skipping Faith's PD "
                        "(pass --tree to include it).", marker or "(no marker)")
            tree = None

    outdir.mkdir(parents=True, exist_ok=True)
    suffix = f"_{marker}" if marker else ""

    # ── Step 1: per-sample read counts + report ──────────────────────────────
    log.info("Extracting per-sample read counts...")
    counts = get_sample_read_counts(table, args.dry_run)
    _, stats = recommend_depth(counts, candidate)
    print_depth_report(counts, stats, candidate)

    depth_tsv = outdir / f"rarefaction_depth_report{suffix}.tsv"
    if not args.dry_run and counts:
        write_depth_tsv(counts, candidate, depth_tsv)

    # ── Step 2: alpha-rarefaction QZV ────────────────────────────────────────
    qzv_out = outdir / f"alpha_rarefaction{suffix}.qzv"
    cmd = [
        "qiime", "diversity", "alpha-rarefaction",
        "--i-table", str(table),
        "--m-metadata-file", str(metadata),
        "--p-max-depth", str(max_depth),
        "--p-steps", str(steps),
        "--p-iterations", str(iterations),
        "--o-visualization", str(qzv_out),
    ]
    if tree is not None:
        cmd += ["--i-phylogeny", str(tree)]

    run_cmd(cmd, dry_run=args.dry_run)

    # ── Step 3: static annotated figure ──────────────────────────────────────
    fig_out = outdir / f"rarefaction_summary{suffix}.png"
    if not args.dry_run:
        make_rarefaction_figure(counts, stats, candidate, marker or "Unknown", fig_out)

    if not args.dry_run:
        produced = [str(p) for p in (qzv_out, fig_out, depth_tsv) if p.exists()]
        checkpoint.print_checkpoint(
            cfg, "rarefaction",
            marker=marker or None,
            produced=produced,
            provenance={
                "inputs": {"table": table, "metadata": metadata, "tree": tree},
                "outputs": produced,
                "command": "python " + " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
                "extra": {"max_depth": max_depth, "candidate_depth": candidate,
                          "p10_reads": stats.get("p10_reads")},
            },
        )
        if counts and "p10_reads" in stats:
            log.info("Suggested starting depth: %s reads (10th percentile, ~90%% retained). "
                     "Confirm against the curves before setting the diversity depth.",
                     f"{stats['p10_reads']:,}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except validate.ValidationError as exc:
        log.error("%s", exc)
        sys.exit(1)
