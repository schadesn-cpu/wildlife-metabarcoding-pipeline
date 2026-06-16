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
  qc                Pre-DADA2 QC: primer detection + demux report
  import            Merge run dirs + make QIIME2 manifests + import FASTQs
  denoise           Cutadapt primer trimming + DADA2 denoising
  metadata          Build per-marker QIIME2 metadata from the source sample sheet
  taxonomy          Taxonomic classification + export count tables
  blast             BLAST refinement / QC of taxonomy        (optional; review loop)
  presence_absence  Detection-frequency analysis              (optional; detection markers)
  rarefaction       Alpha-rarefaction curves to pick a depth  (optional)
  diversity         Core-metrics + group/confound significance (group_diversity)
  figures           All publication figures (both annotated and manuscript sets)

The optional steps (blast, presence_absence, rarefaction) are gated by the
`analyses:` block in pipeline_config.yml and skipped loudly when disabled.
Run `python pipeline.py list` to see the live step list and order.

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
import csv
import json
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
#   loon_project/scripts/       ← SCRIPTS_DIR (pipeline scripts + config_loader.py)
#   loon_project/scripts/config_loader.py
#   loon_project/scripts/make_manifests.py
#   loon_project/scripts/utils/ ← the shared backbone (steps, validate, qc, ...)
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

# Backbone lookups: resolve_script maps a stable logical name to the current
# script filename (so renumbering never strands the orchestrator); require_qiime
# is the environment preflight used by `check` and `run`.
from utils.steps import resolve_script          # noqa: E402
from utils.validate import require_qiime, ValidationError  # noqa: E402
from utils import qc as _qc, provenance as _prov, checkpoint as _ckpt  # noqa: E402
from config_loader import get_paths              # noqa: E402


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
    script_name : filename of the script in scripts/ (e.g. 'make_manifests.py')
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
    Pre-DADA2 QC: primer detection from raw reads, then demux read-count and
    primer-dimer QC.

    Runs primer_advisor.py to identify primers per marker and write a TSV, then
    parse_demux.py for the per-sample read-count and dimer checks. The demux QC
    is treated as non-fatal here (it warns; it does not stop the run), so a dimer
    or low-read flag surfaces without halting a long pipeline.
    """
    paths = get_paths(cfg)
    reads_dir   = paths.reads_dir()
    report_path = paths.primers_detected_tsv()
    report_path.parent.mkdir(parents=True, exist_ok=True)

    run_script(
        resolve_script("primer_advisor"),
        ["--all", "--reads-dir", str(reads_dir), "--report", str(report_path)],
        dry_run=dry_run,
        label="QC: primer detection",
    )
    log.info("  Primer report: %s", report_path)

    # Demux read-count + primer-dimer QC. parse_demux exits non-zero when it
    # flags a FAIL; we surface that as a warning rather than stopping the run.
    try:
        run_script(
            resolve_script("parse_demux"),
            [],
            dry_run=dry_run,
            label="QC: demux read-counts + dimers",
        )
    except RuntimeError as exc:
        log.warning("  demux QC flagged issues (non-fatal) — see the report above: %s", exc)


def step_import(cfg, dry_run: bool) -> None:
    """
    Build QIIME 2 manifests and import demultiplexed FASTQs.

    Calls make_manifests.py to build PairedEndFastqManifestPhred33V2 TSVs from
    the reads/ directory (one per marker, written under qiime2/<marker>/imported/
    by config-derived paths), then calls the run_full engine's import for each
    active marker.
    """
    paths = get_paths(cfg)
    reads_dir = paths.reads_dir()

    # 1. Build manifests (output location is derived from config by the script)
    run_script(
        resolve_script("make_manifests"),
        ["--reads-dir", str(reads_dir),
         "--markers"] + cfg.active_markers,
        dry_run=dry_run,
        label="Build QIIME2 manifests",
    )

    # 2. Import each marker into QIIME2
    for marker in cfg.active_markers:
        manifest = paths.engine_manifest_tsv(marker)
        if not manifest.exists() and not dry_run:
            log.warning("  Manifest not found for %s (%s) — skipping import", marker, manifest)
            continue
        run_script(
            resolve_script("run_full"),
            ["import",
             "--marker",   marker,
             "--manifest", str(manifest),
             "--outdir",   str(paths.qiime2_root())],
            dry_run=dry_run,
            label=f"Import FASTQs: {marker}",
        )


def _detected_primers(cfg, marker: str):
    """
    Return (front_f, front_r) for a marker from the primer_advisor's
    reports/primers_detected.tsv, or None if absent/unreadable.
    """
    path = get_paths(cfg).primers_detected_tsv()
    if not path.exists():
        return None
    try:
        with path.open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("marker") == marker:
                    ff, fr = row.get("fwd_seq"), row.get("rev_seq")
                    if ff and fr:
                        return ff, fr
    except (OSError, csv.Error) as exc:
        log.warning("  could not read primers_detected.tsv: %s", exc)
    return None


def _cutadapt_primer_args(cfg, marker: str) -> List[str]:
    """
    Resolve cutadapt primers for a marker as engine CLI args. Prefers primers
    given explicitly in config (markers.<m>.primers.front_f/front_r); otherwise
    falls back to what the primer advisor detected (reports/primers_detected.tsv).
    Warns and returns [] — no trimming — if neither is available, so the gap is
    loud rather than silent.
    """
    primers = cfg.markers.get(marker, {}).get("primers", {}) or {}
    front_f, front_r = primers.get("front_f"), primers.get("front_r")
    source = "config"
    if not (front_f and front_r):
        detected = _detected_primers(cfg, marker)
        if detected:
            front_f, front_r = detected
            source = "detected (primers_detected.tsv)"
    if not (front_f and front_r):
        log.warning("  no primers for %s (config or detected) — cutadapt will NOT "
                    "trim primers; run the qc step or set markers.%s.primers", marker, marker)
        return []
    log.info("  cutadapt primers for %s from %s: F=%s R=%s", marker, source, front_f, front_r)
    return ["--front-f", front_f, "--front-r", front_r]


def _retention_gate(cfg, marker: str, dataset: str = "all") -> None:
    """
    After DADA2, flag samples that lost too many reads (the merge column is the
    tell for an over-aggressive --trunc-len). Reads the denoising-stats the
    engine wrote for this marker/dataset and runs the retention quality gate.

    Non-fatal by design: the denoise already succeeded, so this warns and logs
    rather than stopping. Skips with a clear note if the stats aren't where
    expected, and a gate error is surfaced (logged), never swallowed.
    """
    stats = get_paths(cfg).engine_denoising_stats_qza(marker, dataset)
    if not stats.exists():
        log.info("  retention gate skipped for %s — no stats at %s", marker, stats)
        return
    try:
        result = _qc.check_retention(
            stats, control_prefixes=cfg.samples.get("control_prefixes", []))
    except (ValueError, OSError) as exc:
        log.warning("  retention gate could not run for %s: %s", marker, exc)
        return
    for line in result.report().splitlines():
        log.info("  %s", line)
    if result.status != "pass":
        log.warning("  [QC] %s retention: %s", marker, result.summary)
    _prov.record_run(
        cfg, "denoise", marker=marker,
        inputs={"denoising_stats": stats},
        extra={"retention_status": result.status, "flagged": len(result.flagged)},
    )


def _dada2_trunc_args(cfg, marker: str, dataset: str = "all") -> List[str]:
    """
    Read the advisor's recommended truncation from dada2_params.json (written by
    the dada2_advisor step) and return them as engine CLI args. Falls back to the
    engine's defaults — with a logged note — when the advisor hasn't been run, so
    denoise still works but you know the trunc-len wasn't tuned.
    """
    params = get_paths(cfg).engine_dada2_params_json(marker, dataset)
    if not params.exists():
        log.info("  no dada2_params.json for %s — using engine default trunc-len "
                 "(run the dada2_advisor step to tune it)", marker)
        return []
    try:
        data = json.loads(params.read_text())
        f = int(data["dada2_args"]["trunc_len_f"])
        r = int(data["dada2_args"]["trunc_len_r"])
    except (ValueError, KeyError, OSError) as exc:
        log.warning("  could not read dada2_params.json for %s: %s — using engine defaults",
                    marker, exc)
        return []
    log.info("  using advisor trunc-len for %s: f=%d r=%d", marker, f, r)
    return ["--trunc-len-f", str(f), "--trunc-len-r", str(r)]


def step_denoise(cfg, dry_run: bool) -> None:
    """
    Cutadapt primer trimming + DADA2 denoising for all active markers.

    Calls the run_full engine's cutadapt and dada2 subcommands per marker, then
    runs the retention quality gate on the DADA2 stats. Truncation comes from the
    dada2_advisor (dada2_params.json) when present, else engine defaults.

    Note: primers come from config (markers.<m>.primers) or the primer advisor's
    report. Run the qc step first if they are not known.
    """
    paths = get_paths(cfg)
    for marker in cfg.active_markers:
        marker_cfg = cfg.markers.get(marker, {})
        run_script(
            resolve_script("run_full"),
            ["cutadapt",
             "--marker", marker,
             "--outdir", str(paths.qiime2_root())]
            + _cutadapt_primer_args(cfg, marker),
            dry_run=dry_run,
            label=f"Cutadapt: {marker}",
        )
        run_script(
            resolve_script("run_full"),
            ["dada2",
             "--marker", marker,
             "--outdir", str(paths.qiime2_root())]
            + _dada2_trunc_args(cfg, marker),
            dry_run=dry_run,
            label=f"DADA2 denoise: {marker}",
        )
        if not dry_run:
            _retention_gate(cfg, marker)


def step_metadata(cfg, dry_run: bool) -> None:
    """
    Build per-marker QIIME 2 metadata after denoising.

    For each active marker, matches the sample-ids in the DADA2 feature table to
    the source sample sheet (metadata.source_sheet, keyed on
    metadata.source_id_column) and writes the per-marker metadata TSV the
    diversity stage reads. If the source sheet isn't configured the stage is
    skipped with a clear warning, so the rest of a run can proceed and you can
    supply the metadata by hand.
    """
    source = cfg.metadata.get("source_sheet")
    id_col = cfg.metadata.get("source_id_column")
    if not source or not id_col:
        log.warning(
            "  metadata stage skipped — set metadata.source_sheet and "
            "metadata.source_id_column in the config to build metadata automatically."
        )
        return

    source_path = cfg.resolve(source)
    paths = get_paths(cfg)
    for marker in cfg.active_markers:
        table = paths.engine_table_qza(marker, "all")
        out   = cfg.metadata.get(marker, {}).get("all", "")
        if not out:
            log.warning("  no metadata.%s.all path configured — skipping %s", marker, marker)
            continue
        run_script(
            resolve_script("make_metadata"),
            ["--marker",           marker,
             "--table",            str(table),
             "--source-metadata",  str(source_path),
             "--source-id-column", id_col,
             "--out",              str(cfg.resolve(out))],
            dry_run=dry_run,
            label=f"Build metadata: {marker}",
        )


def _coverage_gate(cfg, marker: str, dataset: str = "all") -> None:
    """
    After taxonomy assignment, flag a classifier that left too many features
    unassigned — the cue to BLAST-verify (07b/07c/07d) rather than trust it.
    Reads the engine's taxonomy.qza. Non-fatal: warns and logs, never stops.
    """
    tax = get_paths(cfg).engine_taxonomy_qza(marker, dataset)
    if not tax.exists():
        log.info("  coverage gate skipped for %s — no taxonomy at %s", marker, tax)
        return
    try:
        result = _qc.check_classifier_coverage(tax)
    except (ValueError, OSError) as exc:
        log.warning("  coverage gate could not run for %s: %s", marker, exc)
        return
    log.info("  %s", result.report())
    if result.status != "pass":
        log.warning("  [QC] %s classifier coverage: %s", marker, result.summary)
    _prov.record_run(cfg, "taxonomy", marker=marker, inputs={"taxonomy": tax},
                     extra={"coverage_status": result.status})


def step_taxonomy(cfg, dry_run: bool) -> None:
    """
    Taxonomic classification and table export for all active markers.

    Calls the run_full engine's taxonomy subcommand to classify, runs the
    classifier-coverage quality gate, then exports human-readable count TSVs.
    """
    paths = get_paths(cfg)
    for marker in cfg.active_markers:
        marker_cfg  = cfg.markers.get(marker, {})
        classifier  = cfg.resolve(marker_cfg.get("classifier", ""))
        results_dir = paths.engine_taxonomy_results_dir(marker, "all")

        # Classify
        run_script(
            resolve_script("run_full"),
            ["taxonomy",
             "--marker",     marker,
             "--classifier", str(classifier),
             "--outdir",     str(paths.qiime2_root())],
            dry_run=dry_run,
            label=f"Classify taxonomy: {marker}",
        )

        if not dry_run:
            _coverage_gate(cfg, marker)

        # Export to TSV count tables
        run_script(
            resolve_script("taxonomy_table"),
            ["--taxonomy", str(paths.engine_taxonomy_qza(marker, "all")),
             "--table",    str(paths.engine_table_qza(marker, "all", nocontrols=True)),
             "--marker",   marker,
             "--outdir",   str(results_dir)],
            dry_run=dry_run,
            label=f"Export taxonomy tables: {marker}",
        )

    if not dry_run:
        _ckpt.print_checkpoint(cfg, "taxonomy")


def step_blast(cfg, dry_run: bool) -> None:
    """
    BLAST refinement + QC (optional), controlled by analyses.blast.enabled.

    Runs three independent, advisory tools per configured marker: blast_refine
    (push under-resolved calls down a rank), blast_qc (flag classifier conflicts),
    and blast_verify (verify suspect taxa). Each self-derives its inputs from the
    marker + config. A failure in one tool is logged but does not block the
    others — they answer different questions.

    Option (b) flow: blast_refine may write refined_taxonomy_<marker>.tsv; on the
    NEXT taxonomy build the export prefers it (loudly — see TAXONOMY_SOURCE_*.md).
    So after this stage: review the reports in results/<marker>/all/blast/, then
    re-run the taxonomy stage (and downstream) to propagate refined assignments.
    Skipped loudly when disabled.
    """
    bcfg = cfg.analyses.get("blast", {})
    if not bcfg.get("enabled", False):
        log.warning("  blast stage skipped — set analyses.blast.enabled: true "
                    "(and analyses.blast.db) in the config.")
        return
    if not bcfg.get("db") and dry_run is False:
        log.warning("  analyses.blast.db is empty — the BLAST tools will fail loud "
                    "until you point it at a local nt database.")
    markers = bcfg.get("markers") or cfg.active_markers
    for marker in markers:
        if marker not in cfg.active_markers:
            log.warning("  blast: %s not in active_markers — skipping.", marker)
            continue
        for tool, label in [("blast_refine", "refine"),
                            ("blast_qc", "QC"),
                            ("blast_verify", "verify")]:
            try:
                run_script(resolve_script(tool), ["--marker", marker],
                           dry_run=dry_run, label=f"BLAST {label}: {marker}")
            except RuntimeError as e:
                log.warning("  BLAST %s for %s did not complete: %s", label, marker, e)
    if not dry_run:
        log.warning("BLAST stage done. Review reports in results/<marker>/all/blast/. "
                    "To propagate any refined taxonomy, re-run the taxonomy stage "
                    "(and downstream): python pipeline.py run --steps taxonomy,...")


def step_presence_absence(cfg, dry_run: bool) -> None:
    """
    Presence/absence detection analysis from the taxonomy count tables.

    Optional stage, controlled by analyses.presence_absence.enabled. Runs the
    presence_absence script for each configured marker (it derives its count
    table, metadata, output dir, and thresholds from the marker + config).
    Skipped loudly when disabled.
    """
    pcfg = cfg.analyses.get("presence_absence", {})
    if not pcfg.get("enabled", False):
        log.warning("  presence/absence stage skipped — set "
                    "analyses.presence_absence.enabled: true in the config.")
        return
    markers = pcfg.get("markers") or cfg.active_markers
    for marker in markers:
        if marker not in cfg.active_markers:
            log.warning("  presence/absence: %s not in active_markers — skipping.", marker)
            continue
        run_script(
            resolve_script("presence_absence"),
            ["--marker", marker],
            dry_run=dry_run,
            label=f"Presence/absence: {marker}",
        )


def step_rarefaction(cfg, dry_run: bool) -> None:
    """
    Alpha-rarefaction curves to choose a sampling depth before diversity.

    Optional stage, controlled by analyses.rarefaction.enabled in the config.
    Runs the rarefaction script per active marker (it derives its table,
    metadata, tree, depth parameters, and output dir from the marker + config).
    Skipped loudly when disabled so a normal run proceeds untouched.
    """
    rcfg = cfg.analyses.get("rarefaction", {})
    if not rcfg.get("enabled", False):
        log.warning("  rarefaction stage skipped — set analyses.rarefaction.enabled: "
                    "true in the config to generate curves.")
        return
    for marker in cfg.active_markers:
        run_script(
            resolve_script("rarefaction"),
            ["--marker", marker],
            dry_run=dry_run,
            label=f"Rarefaction curves: {marker}",
        )


def step_diversity(cfg, dry_run: bool) -> None:
    """
    Rarefaction curves + core-metrics + group significance tests.

    For each active marker:
      1. Generates alpha-rarefaction QZV
      2. Runs core-metrics-[phylogenetic] at the configured depth
      3. Runs PERMANOVA + alpha group-significance (+ confound checks) via
         the group_diversity script
    """
    paths = get_paths(cfg)
    for marker in cfg.active_markers:
        marker_cfg = cfg.markers.get(marker, {})
        depth      = marker_cfg.get("rarefaction_depth", 5000)
        phylo      = marker_cfg.get("phylo", False)
        meta_path  = cfg.resolve(cfg.metadata.get(marker, {}).get("all", ""))

        # Core metrics
        div_dir = paths.engine_diversity_dir(marker, "all")
        run_script(
            resolve_script("run_full"),
            ["diversity",
             "--marker",   marker,
             "--depth",    str(depth),
             "--metadata", str(meta_path),
             "--outdir",   str(paths.qiime2_root())]
            + (["--phylo"] if phylo else []),
            dry_run=dry_run,
            label=f"Core metrics: {marker}",
        )

        # Group significance tests
        from config_loader import get_diversity_dir
        core_dir = get_diversity_dir(cfg, marker, "DvT")
        results_dir = paths.engine_diversity_results_dir(marker, "DvT")

        run_script(
            resolve_script("group_diversity"),
            ["--marker",      marker,
             "--dataset",     "DvT",
             "--metadata",    str(meta_path),
             "--metrics-dir", str(core_dir),
             "--no-figures"]                      # figures handled by step_figures
            + (["--phylo"] if phylo else []),
            dry_run=dry_run,
            label=f"Group diversity (DvT): {marker}",
        )

        if not dry_run:
            _prov.record_run(
                cfg, "diversity", marker=marker, analysis="DvT", read_depth=depth,
                inputs={"metadata": meta_path, "core_metrics": core_dir},
            )
        # group_diversity.py prints the per-marker "diversity" checkpoint itself,
        # so no separate stage checkpoint is printed here (avoids duplicates).


def step_figures(cfg, dry_run: bool) -> None:
    """
    Generate all publication figures — both annotated (lab) and manuscript sets.

    Calls run_all_figures_v2.py, which reads this same config file and
    generates PCoA and alpha diversity panels for all markers and analyses.
    """
    fig_script = _script(resolve_script("figures"))

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

ALL_STEPS = ["qc", "import", "denoise", "metadata", "taxonomy", "blast", "presence_absence", "rarefaction", "diversity", "figures"]

STEP_FUNCTIONS = {
    "qc":               step_qc,
    "import":           step_import,
    "denoise":          step_denoise,
    "metadata":         step_metadata,
    "taxonomy":         step_taxonomy,
    "blast":            step_blast,
    "presence_absence": step_presence_absence,
    "rarefaction":      step_rarefaction,
    "diversity":        step_diversity,
    "figures":          step_figures,
}

STEP_DESCRIPTIONS = {
    "qc":               "Primer detection + pre-DADA2 QC report",
    "import":           "Build manifests + QIIME2 FASTQ import",
    "denoise":          "Cutadapt primer trimming + DADA2 denoising",
    "metadata":         "Build per-marker QIIME2 metadata from the source sample sheet",
    "taxonomy":         "Taxonomic classification + export count tables",
    "blast":            "BLAST refinement / QC of taxonomy (optional; review loop)",
    "presence_absence": "Detection-frequency analysis for detection markers (optional)",
    "rarefaction":      "Alpha-rarefaction curves to choose a depth (optional)",
    "diversity":        "Core-metrics + group/confound significance (group_diversity)",
    "figures":          "All publication figures (annotated + manuscript)",
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

    # Environment preflight — reported here, enforced in `run`. `check` itself
    # does not need QIIME, so a missing install is a prominent warning, not a
    # failure (you can validate a config from a login node).
    try:
        print(f"  ✓ QIIME 2 detected: {require_qiime()}")
    except ValidationError as e:
        print(f"  [WARN]  {e}")
    print()

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

    # Environment preflight: the steps shell out to QIIME 2, so refuse to start
    # a real run without it. Dry-runs only print commands, so they're exempt.
    if not args.dry_run:
        try:
            qiime_version = require_qiime()
            log.info("QIIME 2 detected: %s", qiime_version)
        except ValidationError as e:
            log.error("%s", e)
            return 1
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
