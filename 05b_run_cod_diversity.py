#!/usr/bin/env python3
"""
05b_run_cod_diversity.py

Run diversity group-significance tests stratified by cause-of-death (COD)
category and collection source. Designed to run after 05_run_diversity_stats.py
using the same core-metrics directories and updated metadata TSVs.

Produces:
  - PERMANOVA + PERMDISP: COD_broad (Lead / Parasitic_Infectious / Trauma)
  - PERMANOVA + PERMDISP: Collection_source (NHVDL / LPC / CFW) — confound check
  - Alpha group-significance: COD_broad and Collection_source
  - PCoA figures colored by COD_broad (both palettes)
  - Alpha figures grouped by COD_broad (both palettes)

COD_broad categories:
  Lead                  n=5  (confirmed or suspected Pb toxicity)
  Parasitic_Infectious  n=11 (malaria, coccidiosis, aspergillosis, sepsis, etc.)
  Trauma                n=16 (controls)
  Marine                n=5  (excluded from DvT COD analyses)
  Unknown_Other         n=3  (open COD or neoplasia — included in metadata,
                               not the primary comparison group)

Usage (from loon_project root):
  conda activate qiime2-amplicon-2024.5

  python scripts/05b_run_cod_diversity.py \\
      --marker 16S \\
      --metrics-dir qiime2/16S/rarefied_8000/DvT/diversity/core_metrics_depth8000 \\
      --metadata metadata/qiime/metadata_16S_updated.tsv \\
      --phylo

  python scripts/05b_run_cod_diversity.py \\
      --marker MiFish \\
      --metrics-dir qiime2/MiFish/rarefied_17000/DvT/diversity/core-metrics-17000 \\
      --metadata metadata/qiime/metadata_MiFish_updated.tsv

  python scripts/05b_run_cod_diversity.py \\
      --marker cytb \\
      --metrics-dir qiime2/cytb/all/diversity/core-metrics-200 \\
      --metadata metadata/qiime/metadata_cytb_updated.tsv

  python scripts/05b_run_cod_diversity.py \\
      --marker 18S \\
      --metrics-dir qiime2/18S/all/diversity/core-metrics-1000 \\
      --metadata metadata/qiime/metadata_18S_updated.tsv

  # Dry run to verify paths:
  python scripts/05b_run_cod_diversity.py \\
      --marker MiFish \\
      --metrics-dir qiime2/MiFish/rarefied_17000/DvT/diversity/core-metrics-17000 \\
      --metadata metadata/qiime/metadata_MiFish_updated.tsv \\
      --dry-run

Outputs:
  results/{marker}/COD/diversity/   QZV files (view at view.qiime2.org)
  results/{marker}/COD/figures/     PNG + SVG figures

Notes:
  - Requires updated metadata TSVs with COD_broad and Collection_source columns
  - Lead n=5 yields low statistical power — treat results as exploratory
  - After running, use parse_beta_stats.py on COD/diversity/ to extract F and p
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
# Logging — identical setup to 05_run_diversity_stats.py
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Utilities — identical to 05_run_diversity_stats.py
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
            f"Activate your QIIME2 conda env first:\n"
            f"  conda activate qiime2-amplicon-2024.5"
        )
    return found


def maybe_overwrite(path: Path, force: bool) -> None:
    """Remove existing output if --force is set; raise FileExistsError otherwise."""
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
    """Run a shell command, write stdout to log_path, and append the command to commands_log."""
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
# COD stats
# ---------------------------------------------------------------------------

def run_cod_stats(
    metrics_dir: Path,
    metadata: Path,
    results_dir: Path,
    logs_dir: Path,
    commands_log: Path,
    phylo: bool,
    force: bool,
    dry_run: bool,
    verbose: bool,
) -> None:

    """
Run PERMANOVA, PERMDISP, and alpha group-significance tests for COD_broad
    and Collection_source groupings.

    Iterates over all beta metrics present in metrics_dir for each grouping
    column, running both PERMANOVA and PERMDISP. Alpha tests are run with
    Kruskal-Wallis. All QZV outputs are written to results_dir. Skips any
    metric whose distance matrix or vector QZA is not found.
    """
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")

    alpha_metrics = ["observed_features", "shannon", "evenness"]
    if phylo:
        alpha_metrics.append("faith_pd")

    beta_metrics = ["bray_curtis", "jaccard"]
    if phylo:
        beta_metrics += ["weighted_unifrac", "unweighted_unifrac"]

    # Alpha group-significance — run once per metric (no --m-metadata-column;
    # QIIME2 2024.5 includes all metadata columns in the visualization for
    # interactive selection). One QZV covers both COD_broad and Collection_source.
    log.info("── Alpha group-significance (all metadata columns) ──")
    for metric in alpha_metrics:
        vec = metrics_dir / f"{metric}_vector.qza"
        if not vec.exists():
            log.warning("  Skipping (not found): %s", vec.name)
            continue
        out = results_dir / f"{metric}_group_sig_all_columns.qzv"
        maybe_overwrite(out, force)
        safe_mkdir(out.parent)
        cmd = [
            "qiime", "diversity", "alpha-group-significance",
            "--i-alpha-diversity", str(vec),
            "--m-metadata-file",   str(metadata),
            "--o-visualization",   str(out),
        ]
        run_cmd(cmd,
                logs_dir / f"alpha_all_{metric}_{stamp}.log",
                commands_log, dry_run=dry_run, verbose=verbose)
        log.info("  Saved: %s (select COD_broad or Collection_source in viewer)", out.name)

    for group_col in ("COD_broad", "Collection_source"):

        log.info("── Beta group-significance: %s (PERMANOVA + PERMDISP) ──", group_col)
        for metric in beta_metrics:
            dm = metrics_dir / f"{metric}_distance_matrix.qza"
            if not dm.exists():
                log.warning("  Skipping (not found): %s", dm.name)
                continue
            for method in ("permanova", "permdisp"):
                out = results_dir / f"{metric}_{method}_{group_col}.qzv"
                maybe_overwrite(out, force)
                cmd = [
                    "qiime", "diversity", "beta-group-significance",
                    "--i-distance-matrix",   str(dm),
                    "--m-metadata-file",     str(metadata),
                    "--m-metadata-column",   group_col,
                    "--p-method",            method,
                    "--p-pairwise",
                    "--p-permutations",      "999",
                    "--o-visualization",     str(out),
                ]
                run_cmd(cmd,
                        logs_dir / f"beta_{group_col}_{metric}_{method}_{stamp}.log",
                        commands_log, dry_run=dry_run, verbose=verbose)
                log.info("  Saved: %s", out.name)

    log.info("── Done (stats) ─────────────────────────────────────────────")


def run_cod_figures(
    marker: str,
    metrics_dir: Path,
    metadata: Path,
    figures_dir: Path,
    logs_dir: Path,
    commands_log: Path,
    phylo: bool,
    palettes: List[str],
    dry_run: bool,
    verbose: bool,
) -> None:

    """
Generate PCoA and alpha diversity figures coloured by COD_broad.

    Calls 06_plot_diversity.py as a subprocess for both purple and wong
    palettes. Skips if the metrics directory does not exist. Output files
    follow the naming convention {marker}_r{depth}_cod_{type}_{palette}.png/svg.
    """
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_mkdir(figures_dir)

    pcoa_names = ["bray_curtis", "jaccard"]
    if phylo:
        pcoa_names += ["weighted_unifrac", "unweighted_unifrac"]

    pcoa_artifacts = []
    for name in pcoa_names:
        p = metrics_dir / f"{name}_pcoa_results.qza"
        if p.exists():
            pcoa_artifacts.append(str(p))
        else:
            log.warning("  PCoA artifact not found (skipping): %s", p.name)

    alpha_names = ["observed_features", "shannon", "evenness"]
    if phylo:
        alpha_names.append("faith_pd")

    alpha_artifacts = [
        str(metrics_dir / f"{n}_vector.qza")
        for n in alpha_names
        if (metrics_dir / f"{n}_vector.qza").exists()
    ]

    log.info("── Figures: COD_broad PCoA + alpha ─────────────────────────")
    for palette in palettes:

        if pcoa_artifacts:
            stem = f"{marker}_COD_pcoa_{palette}"
            cmd = (
                ["python", "scripts/06_plot_diversity.py", "pcoa"]
                + ["--artifact"] + pcoa_artifacts
                + [
                    "--metadata",    str(metadata),
                    "--color-by",    "COD_broad",
                    "--panel",
                    "--palette",     palette,
                    "--no-title",
                    "--output-stem", stem,
                    "--output-dir",  str(figures_dir),
                ]
            )
            run_cmd(cmd,
                    logs_dir / f"fig_COD_pcoa_{palette}_{stamp}.log",
                    commands_log, dry_run=dry_run, verbose=verbose)
            log.info("  PCoA %s -> %s/%s.png", palette, figures_dir, stem)

        if alpha_artifacts:
            stem = f"{marker}_COD_alpha_{palette}"
            cmd = (
                ["python", "scripts/06_plot_diversity.py", "alpha"]
                + ["--artifact"] + alpha_artifacts
                + [
                    "--metadata",    str(metadata),
                    "--group-by",    "COD_broad",
                    "--panel",
                    "--palette",     palette,
                    "--no-title",
                    "--output-stem", stem,
                    "--output-dir",  str(figures_dir),
                ]
            )
            run_cmd(cmd,
                    logs_dir / f"fig_COD_alpha_{palette}_{stamp}.log",
                    commands_log, dry_run=dry_run, verbose=verbose)
            log.info("  Alpha %s -> %s/%s.png", palette, figures_dir, stem)

    log.info("── Done (figures) ───────────────────────────────────────────")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser for 05b_run_cod_diversity.py."""
    p = argparse.ArgumentParser(
        prog="05b_run_cod_diversity.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--marker",      required=True,
                   help="Marker name (16S, MiFish, cytb, 18S).")
    p.add_argument("--metrics-dir", required=True,
                   help="Path to core-metrics output directory.")
    p.add_argument("--metadata",    required=True,
                   help="Path to updated QIIME2 metadata TSV (with COD_broad column).")
    p.add_argument("--phylo", action="store_true", default=False,
                   help="Include phylogenetic metrics (faith_pd, UniFrac). 16S only.")
    p.add_argument("--palette", nargs="+", default=["purple", "wong"],
                   choices=["purple", "wong", "redblue"],
                   help="Figure palette(s). Default: purple wong.")
    p.add_argument("--project-root", default=".",
                   help="Project root. Default: current directory.")
    p.add_argument("--results-root", default="results",
                   help="Results root. Default: results/")
    p.add_argument("--logs-root", default="logs",
                   help="Logs root. Default: logs/")
    p.add_argument("--no-figures", action="store_true", default=False,
                   help="Skip figure generation (stats QZVs only).")
    p.add_argument("--force", action="store_true", default=False,
                   help="Overwrite existing outputs.")
    p.add_argument("--dry-run", action="store_true", default=False,
                   help="Print commands without executing.")
    p.add_argument("--verbose", action="store_true", default=False,
                   help="Print each command to stderr.")
    return p


def main(argv: Optional[List] = None) -> int:
    """
Parse arguments and run COD-stratified diversity stats and figures.

    Validates inputs, sets up output directories, calls run_cod_stats for
    each marker/metrics-dir pair, and optionally calls run_cod_figures.
    Returns 0 on success, 1 on error.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    project_root = Path(args.project_root).resolve()
    os.chdir(project_root)

    metadata = Path(args.metadata)
    if not metadata.is_absolute():
        metadata = (project_root / metadata).resolve()
    require_exists(metadata, "Metadata file (--metadata)")

    metrics_dir = Path(args.metrics_dir)
    if not metrics_dir.is_absolute():
        metrics_dir = (project_root / metrics_dir).resolve()
    require_exists(metrics_dir, "Core-metrics directory (--metrics-dir)")

    results_dir  = project_root / args.results_root / args.marker / "COD" / "diversity"
    figures_dir  = project_root / args.results_root / args.marker / "COD" / "figures"
    logs_dir     = project_root / args.logs_root    / args.marker / "COD"
    commands_log = results_dir / "commands.sh"

    which_or_die("qiime")

    log.info("Marker      : %s", args.marker)
    log.info("Metrics dir : %s", metrics_dir)
    log.info("Metadata    : %s", metadata)
    log.info("Phylo       : %s", args.phylo)
    log.info("Palettes    : %s", args.palette)
    log.info("Results dir : %s", results_dir)
    log.info("Figures dir : %s", figures_dir)
    if args.dry_run:
        log.info("DRY RUN — no commands will be executed")

    try:
        run_cod_stats(
            metrics_dir=metrics_dir, metadata=metadata,
            results_dir=results_dir, logs_dir=logs_dir,
            commands_log=commands_log, phylo=args.phylo,
            force=args.force, dry_run=args.dry_run, verbose=args.verbose,
        )
        if not args.no_figures:
            run_cod_figures(
                marker=args.marker, metrics_dir=metrics_dir,
                metadata=metadata, figures_dir=figures_dir,
                logs_dir=logs_dir, commands_log=commands_log,
                phylo=args.phylo, palettes=args.palette,
                dry_run=args.dry_run, verbose=args.verbose,
            )
    except Exception as ex:
        eprint(f"\n[ERROR] {ex}")
        return 1

    log.info("")
    log.info("=== Complete: %s COD diversity ===", args.marker)
    log.info("  QZV stats  : %s/", results_dir)
    log.info("  Figures    : %s/", figures_dir)
    log.info("")
    log.info("Next steps:")
    log.info("  1. Open QZV files at view.qiime2.org")
    log.info("  2. Run parse_beta_stats.py on %s/COD/diversity/", args.marker)
    log.info("  3. Check Collection_source PERMANOVA — if significant, note as confound in Methods")
    log.info("  4. Lead n=5: report as exploratory, avoid overinterpreting significance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
