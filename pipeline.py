#!/usr/bin/env python3
"""
pipeline.py
===========
Wildlife Metabarcoding Pipeline — Main Entry Point
MEED Lab, University of New Hampshire

A modular amplicon sequencing pipeline for dietary and gut community
metabarcoding in wildlife. Supports 16S, MiFish 12S, cytochrome b,
18S, ITS, and viral detection (adenovirus, herpesvirus).

Usage
-----
  # Create a starter config for a new project
  python pipeline.py init

  # Validate your config and check file paths
  python pipeline.py check

  # Run the full pipeline from scratch
  python pipeline.py run --steps all

  # Run specific steps only
  python pipeline.py run --steps taxonomy diversity figures

  # Regenerate figures from existing diversity results
  python pipeline.py figures

  # Use a non-default config file
  python pipeline.py --config /path/to/my_config.yml run --steps all

  # Preview commands without executing (dry run)
  python pipeline.py run --steps all --dry-run

Pipeline steps (run in order)
------------------------------
  qc          Pre-DADA2 QC: primer detection + demux report
  import      Merge run dirs + make QIIME2 manifests + import FASTQs
  denoise     Cutadapt primer trimming + DADA2 denoising
  taxonomy    Taxonomic classification + export count tables
  diversity   Rarefaction curves + core-metrics + group significance tests
  cod         Cause-of-death stratified diversity (secondary analysis)
  figures     All publication figures (both annotated and manuscript sets)

Each step is idempotent — re-running skips files that already exist.
Use --force to overwrite existing outputs.

Requirements
------------
  conda activate qiime2-amplicon-2024.5
  python pipeline.py --help

See README.md for full installation and configuration instructions.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
import traceback
from pathlib import Path
from typing import List, Optional

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  [%(levelname)s]  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Directory layout
#
#   loon_project/               ← PIPELINE_DIR (where pipeline.py lives)
#   loon_project/pipeline.py
#   loon_project/pipeline_config.yml
#   loon_project/scripts/       ← SCRIPTS_DIR (numbered scripts + config_loader.py)
#   loon_project/scripts/config_loader.py
#   loon_project/scripts/00_build_classifiers.py
#   loon_project/scripts/01_make_manifests.py
#   ...
# ---------------------------------------------------------------------------

PIPELINE_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR  = PIPELINE_DIR / "scripts"

# Fall back to the same directory as pipeline.py if scripts/ doesn't exist
# (supports running from a flat layout during development)
if not SCRIPTS_DIR.exists():
    SCRIPTS_DIR = PIPELINE_DIR

# Add scripts/ to sys.path so config_loader (which lives there) can be
# imported from any working directory without PYTHONPATH manipulation.
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


def _script(name: str) -> Path:
    """
    Resolve a named pipeline script.

    Looks in scripts/ first, then in the same directory as pipeline.py.
    Raises FileNotFoundError with a helpful message if not found.
    """
    for d in [SCRIPTS_DIR, PIPELINE_DIR]:
        p = d / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"Pipeline script not found: {name}\n"
        f"  Searched: {SCRIPTS_DIR}  and  {PIPELINE_DIR}\n"
        f"  Make sure all pipeline scripts are in the scripts/ directory."
    )


# ---------------------------------------------------------------------------
# Subprocess runner
# ---------------------------------------------------------------------------

def run_script(
    script_name: str,
    args: List[str],
    dry_run: bool = False,
    label: Optional[str] = None,
) -> None:
    """
    Run a pipeline script as a subprocess using the current Python interpreter.

    Using sys.executable ensures the same conda environment is used for
    all subprocess calls — critical on shared HPC systems where multiple
    Python environments may be active.

    Parameters
    ----------
    script_name : filename of the script in scripts/ (e.g. '08_taxonomy_table.py')
    args        : list of CLI arguments to pass to the script
    dry_run     : if True, print the command without executing
    label       : human-readable step name for logging

    Raises
    ------
    RuntimeError  if the script exits with a non-zero return code
    """
    script_path = _script(script_name)
    cmd         = [sys.executable, str(script_path)] + args
    display     = label or script_name

    log.info("▶ %s", display)
    if dry_run:
        log.info("  [DRY RUN] %s", " ".join(cmd))
        return

    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(
            f"Step '{display}' failed (exit {result.returncode}).\n"
            f"  Script: {script_path}\n"
            f"  Check the output above for details."
        )
    log.info("  ✓ Done: %s", display)


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def step_qc(cfg, dry_run: bool) -> None:
    """
    Pre-DADA2 QC: primer detection from raw reads.

    Runs primer_advisor.py detect --all to identify primers in each
    marker's reads directory and write a TSV report. This step should
    be run before any trimming to confirm primer sequences and detect
    potential issues (low detection rate, heterogeneity spacers, dimers).
    """
    reads_dir   = cfg.resolve("reads")
    report_path = cfg.resolve("reports/primers_detected.tsv")
    report_path.parent.mkdir(parents=True, exist_ok=True)

    run_script(
        "primer_advisor.py",
        ["detect", "--all",
         "--reads-dir", str(reads_dir),
         "--report",    str(report_path)],
        dry_run=dry_run,
        label="QC: primer detection",
    )
    log.info("  Primer report: %s", report_path)


def step_import(cfg, dry_run: bool) -> None:
    """
    Build QIIME2 manifests and import demultiplexed FASTQs.

    Calls 01_make_manifests.py to build PairedEndFastqManifestPhred33V2
    TSVs from the reads/ directory, then calls 03_run_full_metabarcoding_pipeline.py
    import for each active marker.
    """
    reads_dir  = cfg.resolve("reads")
    import_dir = cfg.resolve("qiime2/imported")

    # 1. Build manifests
    run_script(
        "01_make_manifests.py",
        ["--reads-dir", str(reads_dir),
         "--outdir",    str(import_dir),
         "--markers"] + cfg.active_markers,
        dry_run=dry_run,
        label="Build QIIME2 manifests",
    )

    # 2. Import each marker into QIIME2
    for marker in cfg.active_markers:
        manifest = import_dir / f"manifest_{marker}.tsv"
        if not manifest.exists() and not dry_run:
            log.warning("  Manifest not found for %s — skipping import", marker)
            continue
        run_script(
            "03_run_full_metabarcoding_pipeline.py",
            ["import",
             "--marker",   marker,
             "--manifest", str(manifest),
             "--outdir",   str(cfg.resolve("qiime2"))],
            dry_run=dry_run,
            label=f"Import FASTQs: {marker}",
        )


def step_denoise(cfg, dry_run: bool) -> None:
    """
    Cutadapt primer trimming + DADA2 denoising for all active markers.

    Calls 03_run_full_metabarcoding_pipeline.py cutadapt and dada2 for
    each marker. Truncation parameters must be set per-marker in the
    pipeline call or configured in the relevant sub-script.

    Note: primer lengths come from primer_advisor detect output (reports/
    primers_detected.tsv). Run the qc step first if they are not known.
    """
    for marker in cfg.active_markers:
        marker_cfg = cfg.markers.get(marker, {})
        run_script(
            "03_run_full_metabarcoding_pipeline.py",
            ["cutadapt",
             "--marker", marker,
             "--outdir", str(cfg.resolve("qiime2"))],
            dry_run=dry_run,
            label=f"Cutadapt: {marker}",
        )
        run_script(
            "03_run_full_metabarcoding_pipeline.py",
            ["dada2",
             "--marker", marker,
             "--outdir", str(cfg.resolve("qiime2"))],
            dry_run=dry_run,
            label=f"DADA2 denoise: {marker}",
        )


def step_taxonomy(cfg, dry_run: bool) -> None:
    """
    Taxonomic classification and table export for all active markers.

    Calls 03_run_full_metabarcoding_pipeline.py taxonomy for classification,
    then 08_taxonomy_table.py to export human-readable count TSVs.
    """
    for marker in cfg.active_markers:
        marker_cfg  = cfg.markers.get(marker, {})
        classifier  = cfg.resolve(marker_cfg.get("classifier", ""))
        results_dir = cfg.resolve(f"results/{marker}/all/taxonomy")

        # Classify
        run_script(
            "03_run_full_metabarcoding_pipeline.py",
            ["taxonomy",
             "--marker",     marker,
             "--classifier", str(classifier),
             "--outdir",     str(cfg.resolve("qiime2"))],
            dry_run=dry_run,
            label=f"Classify taxonomy: {marker}",
        )

        # Export to TSV count tables
        run_script(
            "08_taxonomy_table.py",
            ["--taxonomy", str(cfg.resolve(f"qiime2/{marker}/all/taxonomy/taxonomy.qza")),
             "--table",    str(cfg.resolve(f"qiime2/{marker}/all/dada2/table_no_controls.qza")),
             "--marker",   marker,
             "--outdir",   str(results_dir)],
            dry_run=dry_run,
            label=f"Export taxonomy tables: {marker}",
        )


def step_diversity(cfg, dry_run: bool) -> None:
    """
    Rarefaction curves + core-metrics + group significance tests.

    For each active marker:
      1. Generates alpha-rarefaction QZV
      2. Runs core-metrics-[phylogenetic] at the configured depth
      3. Runs PERMANOVA + alpha group-significance via 05_run_diversity_stats.py
    """
    for marker in cfg.active_markers:
        marker_cfg = cfg.markers.get(marker, {})
        depth      = marker_cfg.get("rarefaction_depth", 5000)
        phylo      = marker_cfg.get("phylo", False)
        meta_path  = cfg.resolve(cfg.metadata.get(marker, {}).get("all", ""))

        # Core metrics
        div_dir = cfg.resolve(f"qiime2/{marker}/all/diversity")
        run_script(
            "03_run_full_metabarcoding_pipeline.py",
            ["diversity",
             "--marker",   marker,
             "--depth",    str(depth),
             "--metadata", str(meta_path),
             "--outdir",   str(cfg.resolve("qiime2"))]
            + (["--phylo"] if phylo else []),
            dry_run=dry_run,
            label=f"Core metrics: {marker}",
        )

        # Group significance tests
        from config_loader import get_diversity_dir
        core_dir = get_diversity_dir(cfg, marker, "DvT")
        results_dir = cfg.resolve(f"results/{marker}/DvT/diversity")

        run_script(
            "05_run_diversity_stats.py",
            ["--marker",       marker,
             "--dataset",      "DvT",
             "--metadata",     str(meta_path),
             "--metrics-dir",  str(core_dir),
             "--group-column", cfg.groups.get("primary", {}).get("column", "Group")]
            + (["--phylo"] if phylo else []),
            dry_run=dry_run,
            label=f"Diversity stats (DvT): {marker}",
        )


def step_cod(cfg, dry_run: bool) -> None:
    """
    Cause-of-death stratified diversity analysis (secondary analysis).

    Runs 05b_run_cod_diversity.py for each active marker, using the
    COD-filtered metadata file specified in pipeline_config.yml.
    """
    from config_loader import get_diversity_dir, get_metadata_path
    for marker in cfg.active_markers:
        marker_cfg = cfg.markers.get(marker, {})
        phylo      = marker_cfg.get("phylo", False)
        try:
            meta_path = get_metadata_path(cfg, marker, "cod")
            core_dir  = get_diversity_dir(cfg, marker, "DvT")
        except ValueError as e:
            log.warning("  Skipping COD for %s: %s", marker, e)
            continue

        if not meta_path.exists():
            log.warning("  COD metadata not found for %s: %s — skipping", marker, meta_path)
            continue

        run_script(
            "05b_run_cod_diversity.py",
            ["--marker",      marker,
             "--metrics-dir", str(core_dir),
             "--metadata",    str(meta_path)]
            + (["--phylo"] if phylo else []),
            dry_run=dry_run,
            label=f"COD diversity: {marker}",
        )


def step_figures(cfg, dry_run: bool) -> None:
    """
    Generate all publication figures — both annotated (lab) and manuscript sets.

    Calls run_all_figures_v2.py, which reads this same config file and
    generates PCoA and alpha diversity panels for all markers and analyses.
    """
    try:
        fig_script = _script("run_all_figures_v2.py")
    except FileNotFoundError:
        fig_script = _script("run_all_figures.py")

    cmd = [sys.executable, str(fig_script),
           "--config", str(_find_config_path())]
    if dry_run:
        cmd.append("--dry-run")

    log.info("▶ Generate all figures")
    if dry_run:
        log.info("  [DRY RUN] %s", " ".join(cmd))
        return

    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(
            "Figure generation failed. Check output above for details."
        )
    log.info("  ✓ Done: all figures")


def _find_config_path() -> Path:
    """Return the path to the active config file."""
    from config_loader import find_config
    return find_config()


# ---------------------------------------------------------------------------
# Step registry
# ---------------------------------------------------------------------------

ALL_STEPS = ["qc", "import", "denoise", "taxonomy", "diversity", "cod", "figures"]

STEP_FUNCTIONS = {
    "qc":        step_qc,
    "import":    step_import,
    "denoise":   step_denoise,
    "taxonomy":  step_taxonomy,
    "diversity": step_diversity,
    "cod":       step_cod,
    "figures":   step_figures,
}

STEP_DESCRIPTIONS = {
    "qc":        "Primer detection + pre-DADA2 QC report",
    "import":    "Build manifests + QIIME2 FASTQ import",
    "denoise":   "Cutadapt primer trimming + DADA2 denoising",
    "taxonomy":  "Taxonomic classification + export count tables",
    "diversity": "Rarefaction curves + core-metrics + group significance",
    "cod":       "Cause-of-death stratified diversity (secondary)",
    "figures":   "All publication figures (annotated + manuscript)",
}


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def cmd_init(args: argparse.Namespace) -> int:
    """
    Create a starter pipeline_config.yml in the current directory.

    Copies the template config from the pipeline directory. If a config
    already exists, prompts before overwriting (unless --force is set).
    """
    template = PIPELINE_DIR / "pipeline_config.yml"
    if not template.exists():
        log.error(
            "Template config not found: %s\n"
            "  Make sure pipeline_config.yml is in the same directory as pipeline.py.",
            template,
        )
        return 1

    dest = Path.cwd() / "pipeline_config.yml"
    if dest.exists() and not args.force:
        answer = input(f"pipeline_config.yml already exists in {Path.cwd()}. Overwrite? [y/N] ")
        if answer.strip().lower() != "y":
            print("Aborting. Existing config preserved.")
            return 0

    shutil.copy2(template, dest)
    print(f"✓ Created: {dest}")
    print()
    print("Next steps:")
    print("  1. Edit pipeline_config.yml — set project.root, markers, metadata paths")
    print("  2. python pipeline.py check       — validate the config")
    print("  3. python pipeline.py run --steps all")
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    """
    Validate the config file and check that key files exist.

    Checks:
      - YAML parses cleanly
      - project.root exists
      - active_markers all have config blocks
      - metadata files exist (warns if missing)
      - diversity directories exist (warns if missing)
      - classifier files exist (warns if missing)

    Returns 0 if no errors found (warnings do not cause non-zero exit).
    """
    import config_loader as cl

    try:
        cfg = cl.load_config(args.config)
    except (FileNotFoundError, ValueError) as e:
        log.error("%s", e)
        return 1

    print(f"\n{'='*60}")
    print(f"  Config: {cl.find_config(args.config)}")
    print(f"  Project: {cfg.name}")
    print(f"  Root:    {cfg.root}")
    print(f"  Markers: {', '.join(cfg.active_markers)}")
    print(f"{'='*60}\n")

    issues = cl.validate_config(cfg)

    if not issues:
        print("  ✓ All checks passed. Config looks valid.")
    else:
        errors   = [i for i in issues if i.startswith("[ERROR]")]
        warnings = [i for i in issues if i.startswith("[WARN]")]
        for iss in errors + warnings:
            print(f"  {iss}")
        if errors:
            print(f"\n  {len(errors)} error(s) must be resolved before running the pipeline.")
            return 1
        else:
            print(f"\n  {len(warnings)} warning(s) — most relate to files not yet generated.")
            print("  These will resolve as you run each pipeline step.")

    print()
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    """
    Run one or more pipeline steps in order.

    Steps are run in the canonical order regardless of the order provided
    on the command line, to prevent dependency violations.
    """
    import config_loader as cl

    try:
        cfg = cl.load_config(args.config)
    except (FileNotFoundError, ValueError) as e:
        log.error("%s", e)
        return 1

    # Resolve which steps to run
    requested = args.steps
    if "all" in requested:
        steps_to_run = ALL_STEPS
    else:
        invalid = [s for s in requested if s not in STEP_FUNCTIONS]
        if invalid:
            log.error(
                "Unknown step(s): %s\n  Valid steps: %s",
                ", ".join(invalid), ", ".join(ALL_STEPS),
            )
            return 1
        # Run in canonical order
        steps_to_run = [s for s in ALL_STEPS if s in requested]

    print(f"\n{'='*60}")
    print(f"  Pipeline: {cfg.name}")
    print(f"  Root:     {cfg.root}")
    print(f"  Markers:  {', '.join(cfg.active_markers)}")
    print(f"  Steps:    {', '.join(steps_to_run)}")
    if args.dry_run:
        print(f"  Mode:     DRY RUN (no commands will execute)")
    print(f"{'='*60}\n")

    n_failed = 0
    for step_name in steps_to_run:
        step_fn = STEP_FUNCTIONS[step_name]
        desc    = STEP_DESCRIPTIONS[step_name]
        print(f"\n── Step: {step_name.upper()} — {desc} {'─'*(40 - len(step_name) - len(desc))}")
        try:
            step_fn(cfg, dry_run=args.dry_run)
        except RuntimeError as exc:
            log.error("[STEP FAILED] %s\n  %s", step_name, exc)
            n_failed += 1
            if not args.keep_going:
                log.error("Stopping pipeline. Use --keep-going to continue after failures.")
                break
        except Exception as exc:
            log.error("[UNEXPECTED ERROR in %s] %s: %s", step_name, type(exc).__name__, exc)
            traceback.print_exc(file=sys.stderr)
            n_failed += 1
            if not args.keep_going:
                break

    print(f"\n{'='*60}")
    if n_failed == 0:
        print(f"  ✓ Pipeline complete. All steps succeeded.")
    else:
        print(f"  ✗ Pipeline finished with {n_failed} failed step(s).")
        print(f"    Review errors above and re-run failed steps individually.")
    print(f"{'='*60}\n")

    return 0 if n_failed == 0 else 1


def cmd_figures(args: argparse.Namespace) -> int:
    """
    Regenerate all figures without re-running any QIIME2 steps.

    Shortcut for: python pipeline.py run --steps figures
    """
    import config_loader as cl

    try:
        cfg = cl.load_config(args.config)
    except (FileNotFoundError, ValueError) as e:
        log.error("%s", e)
        return 1

    try:
        step_figures(cfg, dry_run=getattr(args, "dry_run", False))
    except RuntimeError as exc:
        log.error("%s", exc)
        return 1
    return 0


def cmd_list(args: argparse.Namespace) -> int:
    """Print all available pipeline steps with descriptions."""
    print("\nAvailable pipeline steps (run in this order):\n")
    for step in ALL_STEPS:
        print(f"  {step:<12}  {STEP_DESCRIPTIONS[step]}")
    print()
    print("Usage:")
    print("  python pipeline.py run --steps all")
    print("  python pipeline.py run --steps taxonomy diversity figures")
    print()
    return 0


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser with all subcommands."""
    p = argparse.ArgumentParser(
        prog="pipeline.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--config", default=None, metavar="YML",
        help=(
            "Path to pipeline_config.yml. "
            "Default: auto-discovered in current dir or pipeline dir."
        ),
    )

    sub = p.add_subparsers(dest="command", required=True, metavar="COMMAND")

    # ── init ──────────────────────────────────────────────────────────────────
    ini = sub.add_parser("init", help="Create a starter pipeline_config.yml")
    ini.add_argument("--force", action="store_true",
                     help="Overwrite existing config without prompting")
    ini.set_defaults(func=cmd_init)

    # ── check ─────────────────────────────────────────────────────────────────
    chk = sub.add_parser("check", help="Validate config and check file paths")
    chk.set_defaults(func=cmd_check)

    # ── run ───────────────────────────────────────────────────────────────────
    run = sub.add_parser("run", help="Run pipeline steps")
    run.add_argument(
        "--steps", nargs="+", required=True,
        metavar="STEP",
        help=(
            f"Steps to run. Use 'all' for the full pipeline, or name "
            f"specific steps: {', '.join(ALL_STEPS)}"
        ),
    )
    run.add_argument(
        "--dry-run", action="store_true",
        help="Print commands without executing them",
    )
    run.add_argument(
        "--keep-going", action="store_true",
        help="Continue to the next step even if one fails (default: stop on failure)",
    )
    run.set_defaults(func=cmd_run)

    # ── figures ───────────────────────────────────────────────────────────────
    fig = sub.add_parser("figures",
                         help="Regenerate all figures without re-running QIIME2 steps")
    fig.add_argument("--dry-run", action="store_true")
    fig.add_argument(
        "--markers", nargs="+", default=None,
        help="Limit to specific markers. Default: all active markers in config.",
    )
    fig.set_defaults(func=cmd_figures)

    # ── list ──────────────────────────────────────────────────────────────────
    lst = sub.add_parser("list", help="List all available pipeline steps")
    lst.set_defaults(func=cmd_list)

    return p


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    """
    Entry point for the wildlife metabarcoding pipeline.

    Parses the subcommand and dispatches to the appropriate handler.
    All project-specific settings are read from pipeline_config.yml.

    Returns 0 on success, 1 on error.
    """
    parser = build_parser()
    args   = parser.parse_args(argv)

    try:
        return args.func(args)
    except KeyboardInterrupt:
        print("\n[INTERRUPTED]")
        return 130
    except Exception as exc:
        log.error("Unexpected error: %s: %s", type(exc).__name__, exc)
        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
