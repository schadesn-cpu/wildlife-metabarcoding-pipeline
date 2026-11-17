#!/usr/bin/env python3
"""
08_run_diversity_stats.py

Run alpha and beta diversity group-significance tests on any existing
QIIME2 core-metrics output directory.

Designed to be run AFTER core-metrics or core-metrics-phylogenetic, at
whatever rarefaction depth was chosen for that analysis. Point it at the
directory, specify your grouping column, and it produces all QZVs.

Alpha tests (Kruskal-Wallis):
    - observed_features
    - shannon
    - evenness
    - faith_pd  (if --phylo is set)

Beta tests (PERMANOVA + PERMDISP, pairwise):
    - bray_curtis
    - jaccard
    - weighted_unifrac    (if --phylo is set)
    - unweighted_unifrac  (if --phylo is set)

All outputs are written to:
    results/{marker}/{dataset}/diversity/

Usage examples:

  # MiFish (no tree, non-phylogenetic):
  python scripts/08_run_diversity_stats.py \\
      --marker MiFish \\
      --dataset all \\
      --metadata metadata/qiime/metadata_MiFish.tsv \\
      --metrics-dir qiime2/MiFish/all/diversity/core-metrics-17000 \\
      --group-column Group

  # 16S (phylogenetic, different depth):
  python scripts/08_run_diversity_stats.py \\
      --marker 16S \\
      --dataset all \\
      --metadata metadata/qiime/metadata_16S.tsv \\
      --metrics-dir qiime2/16S/all/diversity/core_metrics_depth5000 \\
      --group-column Group \\
      --phylo

  # Dry run to preview commands without executing:
  python scripts/08_run_diversity_stats.py \\
      --marker MiFish --dataset all \\
      --metadata metadata/qiime/metadata_MiFish.tsv \\
      --metrics-dir qiime2/MiFish/all/diversity/core-metrics-17000 \\
      --group-column Group \\
      --dry-run
"""
from __future__ import annotations

import argparse
import datetime as dt
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def eprint(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def require_exists(path: Path, what: str = "Path") -> None:
    """Raise FileNotFoundError with a descriptive message if path does not exist."""
    if not path.exists():
        raise FileNotFoundError(f"{what} does not exist: {path}")


def safe_mkdir(path: Path) -> None:
    """Create path and any missing parents; no-op if it already exists."""
    path.mkdir(parents=True, exist_ok=True)


def which_or_die(exe: str) -> str:
    """Return the full path to exe, or raise RuntimeError if it is not on PATH."""
    found = shutil.which(exe)
    if not found:
        raise RuntimeError(
            f"Required executable not found on PATH: {exe}\n"
            f"Activate your QIIME2 conda env first."
        )
    return found


def maybe_overwrite(path: Path, force: bool) -> None:
    """Remove existing output if --force, otherwise raise."""
    if path.exists():
        if not force:
            raise FileExistsError(
                f"Output already exists: {path}\n"
                f"Use --force to overwrite."
            )
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def run_cmd(
    cmd: List[str],
    log_path: Path,
    commands_log: Path,
    dry_run: bool = False,
    verbose: bool = False,
) -> None:
    """Run a shell command, streaming output and writing to log file."""
    safe_mkdir(log_path.parent)

    if verbose or dry_run:
        eprint("\n$ " + " \\\n    ".join(cmd))
        eprint(f"  (log: {log_path})")

    if dry_run:
        log_path.write_text("DRY RUN\n" + " ".join(cmd) + "\n")
        safe_mkdir(commands_log.parent)
        with commands_log.open("a") as cf:
            cf.write(" ".join(cmd) + "\n")
        return

    with log_path.open("a", encoding="utf-8") as lf:
        lf.write("\n$ " + " ".join(cmd) + "\n")
        lf.flush()
        proc = subprocess.run(
            list(cmd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        lf.write(proc.stdout)
        lf.flush()

    safe_mkdir(commands_log.parent)
    with commands_log.open("a") as cf:
        cf.write(" ".join(cmd) + "\n")

    sys.stdout.write(proc.stdout)
    sys.stdout.flush()

    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit={proc.returncode}). See log: {log_path}"
        )


# ---------------------------------------------------------------------------
# Core stats logic
# ---------------------------------------------------------------------------

def run_stats(
    metrics_dir: Path,
    metadata: Path,
    results_dir: Path,
    logs_dir: Path,
    commands_log: Path,
    group_column: str,
    phylo: bool,
    force: bool,
    dry_run: bool,
    verbose: bool,
) -> None:
    """Run all alpha and beta group-significance tests."""

    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")

    # ── Alpha diversity ────────────────────────────────────────────────────
    alpha_metrics = ["observed_features", "shannon", "evenness"]
    if phylo:
        alpha_metrics.append("faith_pd")

    log.info("── Alpha diversity (Kruskal-Wallis) ──────────────────────────")
    for metric in alpha_metrics:
        vec = metrics_dir / f"{metric}_vector.qza"
        if not vec.exists():
            log.warning("  Skipping (not found): %s", vec.name)
            continue

        out = results_dir / f"{metric}_group_sig_{group_column}.qzv"
        maybe_overwrite(out, force)
        safe_mkdir(out.parent)

        log_path = logs_dir / f"stats_alpha_{metric}_{stamp}.log"

        cmd = [
            "qiime", "diversity", "alpha-group-significance",
            "--i-alpha-diversity", str(vec),
            "--m-metadata-file", str(metadata),
            "--o-visualization", str(out),
        ]
        run_cmd(cmd, log_path, commands_log, dry_run=dry_run, verbose=verbose)
        log.info("  Saved: %s", out.name)

    # ── Beta diversity ─────────────────────────────────────────────────────
    beta_metrics = ["bray_curtis", "jaccard"]
    if phylo:
        beta_metrics += ["weighted_unifrac", "unweighted_unifrac"]

    log.info("── Beta diversity (PERMANOVA + PERMDISP) ─────────────────────")
    for metric in beta_metrics:
        dm = metrics_dir / f"{metric}_distance_matrix.qza"
        if not dm.exists():
            log.warning("  Skipping (not found): %s", dm.name)
            continue

        for method in ("permanova", "permdisp"):
            out = results_dir / f"{metric}_{method}_{group_column}.qzv"
            maybe_overwrite(out, force)
            safe_mkdir(out.parent)

            log_path = logs_dir / f"stats_beta_{metric}_{method}_{stamp}.log"

            cmd = [
                "qiime", "diversity", "beta-group-significance",
                "--i-distance-matrix", str(dm),
                "--m-metadata-file", str(metadata),
                "--m-metadata-column", group_column,
                "--p-method", method,
                "--p-pairwise",
                "--o-visualization", str(out),
            ]
            run_cmd(cmd, log_path, commands_log, dry_run=dry_run, verbose=verbose)
            log.info("  Saved: %s", out.name)

    log.info("── Done ──────────────────────────────────────────────────────")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser for 08_run_diversity_stats.py."""
    p = argparse.ArgumentParser(
        prog="08_run_diversity_stats.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required
    p.add_argument("--marker",      required=True, help="Marker name (e.g. MiFish, 16S).")
    p.add_argument("--dataset",     required=True, help="Dataset slug (e.g. all, DvT).")
    p.add_argument("--metadata",    required=True, help="Path to QIIME2 metadata TSV.")
    p.add_argument("--metrics-dir", required=True,
                   help="Path to core-metrics output dir (e.g. qiime2/MiFish/all/diversity/core-metrics-17000).")
    p.add_argument("--group-column", required=True,
                   help="Metadata column to test group significance (e.g. Group).")

    # Optional
    p.add_argument("--phylo", action="store_true", default=False,
                   help="Include phylogenetic metrics (faith_pd, UniFrac). Default: off.")
    p.add_argument("--project-root", default=".",
                   help="Project root directory. Default: current directory.")
    p.add_argument("--results-root", default="results",
                   help="Results root relative to project root. Default: results/")
    p.add_argument("--logs-root", default="logs",
                   help="Logs root relative to project root. Default: logs/")
    p.add_argument("--force", action="store_true", default=False,
                   help="Overwrite existing outputs.")
    p.add_argument("--dry-run", action="store_true", default=False,
                   help="Print commands without executing them.")
    p.add_argument("--verbose", action="store_true", default=False,
                   help="Print each command to stderr before running.")

    return p


def main(argv: Optional[list] = None) -> int:
    """
Parse arguments and run all alpha and beta diversity group-significance tests.

    Writes one QZV per metric/method combination to results/{marker}/{dataset}/diversity/.
    Alpha tests use Kruskal-Wallis; beta tests run both PERMANOVA and PERMDISP
    with pairwise comparisons. Returns 0 on success, 1 on error.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    project_root = Path(args.project_root).resolve()
    os.chdir(project_root)

    metadata    = Path(args.metadata)
    if not metadata.is_absolute():
        metadata = (project_root / metadata).resolve()
    require_exists(metadata, "Metadata file (--metadata)")

    metrics_dir = Path(args.metrics_dir)
    if not metrics_dir.is_absolute():
        metrics_dir = (project_root / metrics_dir).resolve()
    require_exists(metrics_dir, "Core-metrics directory (--metrics-dir)")

    results_dir  = project_root / args.results_root / args.marker / args.dataset / "diversity"
    logs_dir     = project_root / args.logs_root    / args.marker / args.dataset
    commands_log = results_dir / "commands.sh"

    which_or_die("qiime")

    log.info("Marker:      %s", args.marker)
    log.info("Dataset:     %s", args.dataset)
    log.info("Metrics dir: %s", metrics_dir)
    log.info("Group col:   %s", args.group_column)
    log.info("Phylo:       %s", args.phylo)
    log.info("Output dir:  %s", results_dir)
    if args.dry_run:
        log.info("DRY RUN — no commands will be executed")

    try:
        run_stats(
            metrics_dir  = metrics_dir,
            metadata     = metadata,
            results_dir  = results_dir,
            logs_dir     = logs_dir,
            commands_log = commands_log,
            group_column = args.group_column,
            phylo        = args.phylo,
            force        = args.force,
            dry_run      = args.dry_run,
            verbose      = args.verbose,
        )
    except Exception as ex:
        eprint(f"\n[ERROR] {ex}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
