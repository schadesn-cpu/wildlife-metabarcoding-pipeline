#!/usr/bin/env python3
"""
utils/steps.py
==============
The pipeline step registry — the single source of truth for execution order.

Every script's "what comes next" footer, the README quick-start, and the
orchestrator's run order all read from this one list. Add or reorder a step
HERE and nothing downstream needs renumbering: no script hardcodes its own
position or the name of the step that follows it.

A Step records:
    key       short stable identifier (used on the CLI: --step <key>)
    title     human-readable name shown in checkpoints
    produces  artifacts/files this step writes (for the checkpoint "Produced:" line)
    inspect   QZV (or other) outputs worth looking at, paired with WHAT to look for
    requires  one-line note on the input format this step needs (for the
              next-step "Heads-up:" line and for validate.py)
    status    "ready"   — rebuilt onto the clean backbone
              "planned" — still on the legacy script, not yet ported

Public API:
    all_steps()        -> list[Step]      in execution order
    get_step(key)      -> Step
    next_step(key)     -> Step | None
    step_index(key)    -> (position, total)   1-based, e.g. (2, 12)
    run_command(key)   -> str             canonical CLI to run a step
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple


@dataclass(frozen=True)
class Inspect:
    """A single inspectable output and the thing the user should check in it."""
    artifact: str        # e.g. "demux_{marker}.qzv"
    look_for: str        # e.g. "per-sample read counts and quality drop-off"
    is_qzv: bool = True  # QZVs get the QIIME2 View hint appended automatically


@dataclass(frozen=True)
class Step:
    key: str
    title: str
    produces: List[str] = field(default_factory=list)
    inspect: List[Inspect] = field(default_factory=list)
    requires: str = ""
    status: str = "planned"
    # The coarse `pipeline.py run --steps <stage>` this step runs under. Empty
    # means it has no orchestrator stage (a manual / advisory step run directly).
    stage: str = ""


# ---------------------------------------------------------------------------
# THE PIPELINE, in scientific execution order.
#
# This order is deliberate. Note in particular that primer removal (cutadapt)
# precedes DADA2 — the "trim" the advisor recommends for cutadapt is primer
# stripping, which is a different operation from DADA2's --trunc-len quality
# truncation. Keeping them as separate, correctly-ordered steps is what stops
# new users from conflating the two.
# ---------------------------------------------------------------------------

_STEPS: List[Step] = [
    Step(
        key="merge_runs",
        title="Merge sequencing runs",
        produces=["reads/{marker}/  (consolidated per-marker read directories)"],
        inspect=[Inspect("reads/ tree and collision_report.tsv",
                         "any duplicate filenames across run directories",
                         is_qzv=False)],
        requires="one subdirectory of demultiplexed fastq.gz per sequencing run",
        status="ready",
    ),
    Step(
        key="manifest",
        stage="import",
        title="Build manifests",
        produces=["qiime2/{marker}/imported/manifest_{marker}.tsv"],
        requires="reads/<marker>/ folders of paired *_R1_*/*_R2_* fastq.gz files",
        status="ready",
    ),
    Step(
        key="metadata",
        stage="metadata",
        title="Build QIIME 2 metadata",
        produces=["metadata/qiime/metadata_{marker}.tsv"],
        requires="a sample sheet with a sample-id column and your grouping columns",
        status="ready",
    ),
    Step(
        key="import",
        stage="import",
        title="Import reads to QIIME 2",
        produces=["qiime2/{marker}/imported/demux_{marker}.qza",
                  "qiime2/{marker}/imported/demux_{marker}.qzv"],
        inspect=[Inspect("demux_{marker}.qzv",
                         "per-sample read counts and where quality scores drop off")],
        requires="a valid manifest TSV (run the manifest step first)",
        status="planned",
    ),
    Step(
        key="parse_demux",
        stage="qc",
        title="Parse demux / QC summary",
        produces=["results/qc/demux_read_counts.tsv", "results/qc/demux_qc_report.txt"],
        inspect=[Inspect("demux read-count table",
                         "samples below the read-count threshold, and dimer flags",
                         is_qzv=False)],
        requires="the Illumina reports/ directory (Demultiplex_Stats + MultiQC data)",
        status="ready",
    ),
    Step(
        key="primer_advisor",
        stage="qc",
        title="Primer advisor (detect)",
        produces=["reports/primers_detected.tsv"],
        inspect=[Inspect("the detected forward/reverse primer per marker",
                         "low-confidence matches (<30%) and the recommended cutadapt command",
                         is_qzv=False)],
        requires="raw reads/<marker>/ FASTQ.gz (run before cutadapt)",
        status="ready",
    ),
    Step(
        key="trim",
        stage="denoise",
        title="Trim primers (cutadapt)",
        produces=["qiime2/{marker}/imported/demux_{marker}_trimmed.qza",
                  "qiime2/{marker}/imported/demux_{marker}_trimmed.qzv"],
        inspect=[Inspect("demux_{marker}_trimmed.qzv",
                         "primer-removal summary: fraction of reads kept")],
        requires="the imported demux_{marker}.qza and primer sequences",
        status="planned",
    ),
    Step(
        key="dada2_advisor",
        title="DADA2 truncation advisor",
        produces=["qiime2/{marker}/all/imported/dada2_params.json"],
        inspect=[Inspect("the recommended trunc-len-f / trunc-len-r",
                         "whether forward + reverse still overlap enough to merge",
                         is_qzv=False)],
        requires="the trimmed demux_trimmed.qzv (read quality after primer removal)",
        status="ready",
    ),
    Step(
        key="denoise",
        stage="denoise",
        title="Denoise (DADA2)",
        produces=["qiime2/{marker}/dada2/table_{marker}.qza",
                  "qiime2/{marker}/dada2/rep_seqs_{marker}.qza",
                  "qiime2/{marker}/dada2/denoising_stats_{marker}.qza"],
        inspect=[Inspect("denoising_stats_{marker}.qzv",
                         "the retention column — if <60% of reads survive, "
                         "revisit your --trunc-len")],
        requires="the trimmed demux artifact and the DADA2 params from the advisor",
        status="planned",
    ),
    Step(
        key="taxonomy",
        stage="taxonomy",
        title="Assign taxonomy",
        produces=["qiime2/{marker}/taxonomy/taxonomy_{marker}_{db}.qza",
                  "results/{marker}/taxonomy/.../taxonomy_relabund_L*.tsv"],
        inspect=[Inspect("taxa barplot QZV", "dominant taxa per group; obvious contaminants")],
        requires="rep_seqs_{marker}.qza and the marker's trained classifier",
        status="planned",
    ),
    Step(
        key="blast",
        stage="blast",
        title="BLAST refine / QC",
        produces=["results/{marker}/all/blast/blast_summary_{marker}.tsv",
                  "results/{marker}/all/blast/refined_taxonomy_{marker}.tsv (with --apply)",
                  "results/{marker}/all/blast/blast_refine_report_{marker}.txt"],
        inspect=[Inspect("blast_refine_report",
                         "proposed genus calls — review before --apply")],
        requires="exported taxonomy TSV + rep-seqs FASTA, and a local BLAST nt DB",
        status="planned",
    ),
    Step(
        key="presence_absence",
        stage="presence_absence",
        title="Presence / absence",
        produces=["results/{marker}/all/presence_absence/detection_freq_{marker}.tsv",
                  "results/{marker}/all/presence_absence/detection_barplot_{marker}_*.png"],
        inspect=[Inspect("detection barplot",
                         "per-taxon detection frequency across samples / groups")],
        requires="the taxonomy count TSV from the taxonomy stage",
        status="planned",
    ),
    Step(
        key="rarefaction",
        stage="diversity",
        title="Rarefaction",
        produces=["alpha-rarefaction QZV per marker"],
        inspect=[Inspect("alpha-rarefaction QZV",
                         "where curves plateau — this sets your sampling depth")],
        requires="table_{marker}.qza and metadata",
        status="planned",
    ),
    Step(
        key="diversity",
        stage="diversity",
        title="Diversity (core-metrics)",
        produces=["qiime2/{marker}/diversity/r{depth}_{analysis}/core_metrics/",
                  "results/{marker}/diversity/r{depth}_{analysis}/*_stats_*.tsv"],
        inspect=[Inspect("core-metrics QZVs", "PCoA separation and PERMANOVA group stats")],
        requires="table_{marker}.qza, metadata, and a chosen rarefaction depth",
        status="planned",
    ),
]

# Fast lookup by key, preserving order.
_BY_KEY = {s.key: s for s in _STEPS}


def all_steps() -> List[Step]:
    """Return every step in execution order."""
    return list(_STEPS)


def get_step(key: str) -> Step:
    """Return the Step with this key, or raise a clear error listing valid keys."""
    try:
        return _BY_KEY[key]
    except KeyError:
        valid = ", ".join(s.key for s in _STEPS)
        raise KeyError(f"Unknown step {key!r}. Valid steps, in order: {valid}")


def next_step(key: str) -> Optional[Step]:
    """Return the step that follows `key`, or None if `key` is the last step."""
    keys = [s.key for s in _STEPS]
    i = keys.index(key)  # raises ValueError if key is unknown — intentional
    return _STEPS[i + 1] if i + 1 < len(_STEPS) else None


def step_index(key: str) -> Tuple[int, int]:
    """Return (1-based position, total) for `key`, e.g. (2, 12)."""
    keys = [s.key for s in _STEPS]
    return keys.index(key) + 1, len(keys)


def run_command(key: str) -> str:
    """
    A correct, runnable CLI hint for a step.

    Steps that belong to an orchestrator stage return the real
    `pipeline.py run --steps <stage>` command; manual/advisory steps (no stage)
    point at the script to run directly.
    """
    step = get_step(key)
    stage = getattr(step, "stage", "") if step else ""
    if stage:
        return f"python pipeline.py run --steps {stage}"
    script = SCRIPTS.get(key)
    if script:
        return f"python scripts/{script}   (advisory/manual — run directly)"
    return "(no orchestrator command for this step)"


# ---------------------------------------------------------------------------
# Script-name resolution
#
# The single place that maps a stable logical name to the script file that
# currently implements it. pipeline.py looks names up here instead of
# hardcoding filenames, so renaming or renumbering a script is a one-line edit
# in this map — it can never again silently strand the orchestrator.
#
# Filenames only (no paths): the orchestrator's _script() resolves them against
# scripts/ then the pipeline root, so this map stays machine-independent.
#
# As each legacy numbered script is ported onto the clean backbone, update its
# value here to the new descriptive filename. 'make_manifests' already points
# at the ported version.
# ---------------------------------------------------------------------------

SCRIPTS = {
    "make_manifests":  "make_manifests.py",                      # ported (clean backbone)
    "merge_runs":      "merge_runs.py",                          # ported (clean backbone)
    "make_metadata":   "make_metadata.py",                       # ported (clean backbone)
    "primer_advisor":  "primer_advisor.py",                      # ported (detect only; clean backbone)
    "parse_demux":     "parse_demux.py",                         # ported (clean backbone)
    "dada2_advisor":   "dada2_advisor.py",                       # ported (clean backbone)
    "run_full":        "05_run_full_metabarcoding_pipeline.py",
    # rarefaction runs inside the engine's diversity subcommand; no standalone script.
    "taxonomy_table":  "taxonomy_table.py",                      # ported (clean backbone)
    "rarefaction":     "rarefaction.py",                         # ported (clean backbone)
    "presence_absence": "presence_absence.py",                   # ported (clean backbone)
    "blast_refine":    "blast_refine.py",                        # ported (clean backbone)
    "blast_qc":        "blast_qc.py",                            # ported (clean backbone)
    "blast_verify":    "blast_verify.py",                        # ported (clean backbone)
    "group_diversity": "group_diversity.py",                    # ported (unifies diversity_stats + cod_diversity)
    "plot_diversity":  "09_plot_diversity.py",                  # figure helper (called by group_diversity)
    "figures":         "run_all_figures.py",
}


def resolve_script(key: str) -> str:
    """
    Return the current script filename for a logical name.

    Raises KeyError (with the list of valid names) for an unknown key — a typo
    must fail loudly here, not produce a "file not found" deep in a run.
    """
    try:
        return SCRIPTS[key]
    except KeyError:
        valid = ", ".join(sorted(SCRIPTS))
        raise KeyError(f"Unknown script key {key!r}. Valid keys: {valid}")
