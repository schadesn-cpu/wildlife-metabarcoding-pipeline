#!/usr/bin/env python3
"""
group_diversity.py
==================
Step: diversity (group-significance comparisons on an existing core-metrics dir)

Purpose:
    Compare community diversity across grouping variables on a QIIME 2
    core-metrics output. Built around two roles:

      - PRIMARY   the variable you are testing (e.g. Group, Treatment, Site).
      - CONFOUND  optional extra variables you want to show are NOT driving the
                  signal (e.g. CollectionBatch, SequencingRun, Source). Each is
                  tested the same way; if a confound comes back significant, that
                  is your cue to report it as a potential confounder.

    For each variable it runs beta-group-significance (PERMANOVA + PERMDISP,
    pairwise) over every beta metric present, and alpha-group-significance
    (Kruskal-Wallis). It then draws PCoA + alpha figures coloured by the PRIMARY
    variable. Confound checking is optional: list zero confounds (or comment them
    out in the config) and it simply tests the primary variable.

Inputs (derived from --marker + config when omitted):
    --metadata       QIIME 2 metadata TSV (must contain every grouping column)
    --metrics-dir    core-metrics dir (qiime2/<marker>/<dataset>/diversity/core-metrics-<depth>)
    --primary-column    primary variable (default: groups.primary.column)
    --confounds         confound variables (default: groups.confounds)
    --phylo          include phylogenetic metrics (default: markers.<marker>.phylo)

Outputs (in results/<marker>/<dataset>/diversity/):
    <metric>_<method>_<column>.qzv     PERMANOVA/PERMDISP per metric/variable
    <metric>_group_sig_all_columns.qzv alpha group-significance (all columns)
    figures/<marker>_<primary>_*.png    PCoA + alpha figures by primary
    commands.sh                        exact qiime commands run
    logs/run_manifest.jsonl            run appended on completion

Usage:
    python group_diversity.py --marker MiFish
    python group_diversity.py --marker 16S --phylo
    python group_diversity.py --marker cytb --primary-column Group --confounds Source Batch
    python group_diversity.py --marker MiFish --no-figures   # skip figures

Requirements:
    QIIME 2 on PATH. The core-metrics directory must already exist. Figures are
    drawn by the plot_diversity helper.
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
from typing import List, Optional, Tuple

# --- make config_loader and the utils package importable regardless of cwd ---
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from config_loader import load_config, get_paths, get_metadata_path  # noqa: E402
from utils import checkpoint, provenance, validate, steps as _steps  # noqa: E402

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

def run_group_stats(
    metrics_dir: Path,
    metadata: Path,
    results_dir: Path,
    logs_dir: Path,
    commands_log: Path,
    group_columns: List[Tuple[str, str]],
    phylo: bool,
    force: bool,
    dry_run: bool,
    verbose: bool,
) -> None:

    """
    Run PERMANOVA, PERMDISP, and alpha group-significance for each grouping
    variable in group_columns (a list of (column, role) pairs, role being
    'PRIMARY' or 'CONFOUND').

    Iterates over all beta metrics present in metrics_dir for each grouping
    column, running both PERMANOVA and PERMDISP. Alpha tests are run once with
    Kruskal-Wallis (the QZV exposes every metadata column for interactive
    selection). Skips any metric whose distance matrix or vector QZA is absent.
    """
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    columns = [c for c, _ in group_columns]

    alpha_metrics = ["observed_features", "shannon", "evenness"]
    if phylo:
        alpha_metrics.append("faith_pd")

    beta_metrics = ["bray_curtis", "jaccard"]
    if phylo:
        beta_metrics += ["weighted_unifrac", "unweighted_unifrac"]

    # Alpha group-significance — run once per metric (no --m-metadata-column;
    # QIIME2 2024.5 includes all metadata columns in the visualization for
    # interactive selection). One QZV covers every grouping column.
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
        log.info("  Saved: %s (select %s in viewer)", out.name, " / ".join(columns))

    for group_col, role in group_columns:

        log.info("── Beta group-significance: %s [%s] (PERMANOVA + PERMDISP) ──",
                 group_col, role)
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
        if role == "CONFOUND":
            log.info("  ^ confound check: if %s is significant above, your primary "
                     "signal may be confounded — report it in Methods.", group_col)

    log.info("── Done (stats) ─────────────────────────────────────────────")


def run_group_figures(
    marker: str,
    primary_col: str,
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
    Generate PCoA and alpha diversity figures coloured by the primary variable.

    Calls the plot_diversity helper as a subprocess for each palette. Skips
    artifacts that are not present. Output files follow the naming convention
    {marker}_{primary}_{type}_{palette}.png/svg.
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

    log.info("── Figures: %s PCoA + alpha ─────────────────────────", primary_col)

    # Resolve the plot helper via the registry, falling back to a path next to
    # this script. sys.executable keeps the subprocess in the same environment.
    try:
        plot_script = str(_SCRIPTS_DIR / _steps.resolve_script("plot_diversity"))
    except Exception:
        plot_script = str(_SCRIPTS_DIR / "09_plot_diversity.py")

    for palette in palettes:

        if pcoa_artifacts:
            stem = f"{marker}_{primary_col}_pcoa_{palette}"
            cmd = (
                [sys.executable, plot_script, "pcoa"]
                + ["--artifact"] + pcoa_artifacts
                + [
                    "--metadata",    str(metadata),
                    "--color-by",    primary_col,
                    "--panel",
                    "--palette",     palette,
                    "--no-title",
                    "--output-stem", stem,
                    "--output-dir",  str(figures_dir),
                ]
            )
            run_cmd(cmd,
                    logs_dir / f"fig_{primary_col}_pcoa_{palette}_{stamp}.log",
                    commands_log, dry_run=dry_run, verbose=verbose)
            log.info("  PCoA %s -> %s/%s.png", palette, figures_dir, stem)

        if alpha_artifacts:
            stem = f"{marker}_{primary_col}_alpha_{palette}"
            cmd = (
                [sys.executable, plot_script, "alpha"]
                + ["--artifact"] + alpha_artifacts
                + [
                    "--metadata",    str(metadata),
                    "--group-by",    primary_col,
                    "--panel",
                    "--palette",     palette,
                    "--no-title",
                    "--output-stem", stem,
                    "--output-dir",  str(figures_dir),
                ]
            )
            run_cmd(cmd,
                    logs_dir / f"fig_{primary_col}_alpha_{palette}_{stamp}.log",
                    commands_log, dry_run=dry_run, verbose=verbose)
            log.info("  Alpha %s -> %s/%s.png", palette, figures_dir, stem)

    log.info("── Done (figures) ───────────────────────────────────────────")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser for group_diversity.py."""
    p = argparse.ArgumentParser(
        prog="group_diversity.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--marker",      required=True,
                   help="Marker name (e.g. 16S, MiFish, cytb, 18S).")
    p.add_argument("--config",      default=None, help="Path to pipeline_config.yml.")
    p.add_argument("--dataset",     default="all", help="Dataset slug (e.g. all, DvT). Default: all.")
    p.add_argument("--metadata",    default=None,
                   help="QIIME2 metadata TSV (with every grouping column). Derived if omitted.")
    p.add_argument("--metrics-dir", default=None,
                   help="Core-metrics output directory. Derived from --marker + depth if omitted.")
    p.add_argument("--primary-column", default=None,
                   help="Primary variable of interest. Default: groups.primary.column.")
    p.add_argument("--confounds", nargs="*", default=None,
                   help="Confound-check variables. Default: groups.confounds. "
                        "Pass none (or leave empty in config) to test only the primary.")
    p.add_argument("--phylo", action="store_const", const=True, default=None,
                   help="Include phylogenetic metrics (faith_pd, UniFrac). "
                        "Default: markers.<marker>.phylo from config.")
    p.add_argument("--palette", nargs="+", default=None,
                   choices=["purple", "wong", "redblue"],
                   help="Figure palette(s). Default: figures.palette from config.")
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
    Parse arguments and run group-significance stats (+ figures) for the primary
    variable and any confound variables.

    Resolves metadata, the core-metrics dir, grouping columns, phylo flag, and
    palettes from --marker + config when not given. Returns 0 on success, 2 for
    missing inputs, 1 on error.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    marker = args.marker
    dataset = args.dataset
    marker_cfg = cfg.markers.get(marker, {})

    # Derive inputs from the marker + config when not given explicitly.
    metadata = Path(args.metadata) if args.metadata else None
    if metadata is None:
        try:
            metadata = Path(get_metadata_path(cfg, marker, dataset))
        except (KeyError, ValueError):
            metadata = None

    if args.metrics_dir:
        metrics_dir = Path(args.metrics_dir)
    else:
        depth = marker_cfg.get("rarefaction_depth", 5000)
        metrics_dir = paths.engine_diversity_dir(marker, dataset) / f"core-metrics-{depth}"

    primary = args.primary_column or cfg.groups.get("primary", {}).get("column", "Group")
    confounds = args.confounds if args.confounds is not None else cfg.groups.get("confounds", []) or []
    group_columns: List[Tuple[str, str]] = [(primary, "PRIMARY")] + \
        [(c, "CONFOUND") for c in confounds]

    phylo = args.phylo if args.phylo is not None else marker_cfg.get("phylo", False)

    palette = args.palette
    if palette is None:
        cfg_pal = cfg.figures.get("palette", "wong")
        palette = cfg_pal if isinstance(cfg_pal, list) else [cfg_pal]

    results_dir  = paths.engine_diversity_results_dir(marker, dataset)
    figures_dir  = results_dir / "figures"
    logs_dir     = cfg.root / "logs" / marker / dataset
    commands_log = results_dir / "commands.sh"

    # Validate
    validate.require_qiime()
    if metadata is None:
        log.error("No metadata for %s/%s — pass --metadata or set it in the config.",
                  marker, dataset)
        return 2
    if not metadata.exists():
        log.error("Metadata file not found: %s", metadata)
        return 2
    if not metrics_dir.exists():
        log.error("Core-metrics directory not found: %s — run the diversity stage "
                  "(core-metrics) first.", metrics_dir)
        return 2

    log.info("Marker        : %s", marker)
    log.info("Dataset       : %s", dataset)
    log.info("Metrics dir   : %s", metrics_dir)
    log.info("Primary       : %s", primary)
    log.info("Confounds     : %s", ", ".join(confounds) if confounds else "(none)")
    log.info("Phylo         : %s", phylo)
    log.info("Palettes      : %s", palette)
    log.info("Results dir   : %s", results_dir)
    if args.dry_run:
        log.info("DRY RUN — no commands will be executed")

    try:
        run_group_stats(
            metrics_dir=metrics_dir, metadata=metadata,
            results_dir=results_dir, logs_dir=logs_dir,
            commands_log=commands_log, group_columns=group_columns,
            phylo=phylo, force=args.force, dry_run=args.dry_run, verbose=args.verbose,
        )
        if not args.no_figures:
            run_group_figures(
                marker=marker, primary_col=primary, metrics_dir=metrics_dir,
                metadata=metadata, figures_dir=figures_dir,
                logs_dir=logs_dir, commands_log=commands_log,
                phylo=phylo, palettes=palette,
                dry_run=args.dry_run, verbose=args.verbose,
            )
    except Exception as ex:
        eprint(f"\n[ERROR] {ex}")
        return 1

    log.info("")
    log.info("=== Complete: %s diversity (%s) ===", marker, primary)
    log.info("  QZV stats : %s/", results_dir)
    if not args.no_figures:
        log.info("  Figures   : %s/", figures_dir)
    log.info("Next steps:")
    log.info("  1. Open the QZVs at https://view.qiime2.org")
    log.info("  2. Read the PERMANOVA/PERMDISP p-values for %s", primary)
    if confounds:
        log.info("  3. Check the confound PERMANOVAs (%s) — if significant, report "
                 "as a potential confounder in Methods", ", ".join(confounds))

    if not args.dry_run:
        produced = sorted(str(p) for p in results_dir.glob("*.qzv"))
        checkpoint.print_checkpoint(
            cfg, "diversity",
            marker=marker,
            produced=produced,
            provenance={
                "inputs": {"metadata": str(metadata), "metrics_dir": str(metrics_dir)},
                "outputs": produced,
                "command": "python " + " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
                "extra": {"dataset": dataset, "primary": primary,
                          "confounds": confounds, "phylo": phylo},
            },
        )

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except validate.ValidationError as e:
        log.error("%s", e)
        sys.exit(1)
