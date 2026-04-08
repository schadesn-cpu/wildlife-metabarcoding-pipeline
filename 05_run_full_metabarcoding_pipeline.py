#!/usr/bin/env python3
"""
run_full_metabarcoding_pipeline.py

A reproducible, marker-aware QIIME2 wrapper for metabarcoding workflows.

Designed for:
- Multi-marker projects (16S, 18S, MiFish, etc.)
- Multiple dataset subsets (all, DvT, marine_only, etc.)
- Clean separation of:
    qiime2/      → primary artifacts (qza/qzv)
    results/     → derived tables, figures, stats
    logs/        → execution logs

Core capabilities:
- Import reads
- DADA2 denoising
- Taxonomy assignment
- Feature filtering
- Taxonomic collapse
- Diversity analysis (alpha + beta)
- Differential abundance (ANCOM)
- Taxa barplots (ASV-level canonical plot)
- Export and bundling

All outputs are scoped to:
    {marker}/{dataset}

This ensures reproducibility and prevents cross-contamination
between analysis subsets.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, Tuple
from typing import List, Optional, Sequence

import logging

# ---------------------------
# Logging
# ---------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------
# Utilities
# ---------------------------

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
            f"Activate your QIIME2 conda env first (e.g., `conda activate qiime2-amplicon-2024.5`)."
        )
    return found

def deep_update(base: Dict[str, Any], upd: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge upd into base (returns base)."""
    for k, v in upd.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            deep_update(base[k], v)
        else:
            base[k] = v
    return base


def load_config_file(path: Path) -> Dict[str, Any]:
    """Load a config file. Supports YAML (if PyYAML installed) or JSON."""
    require_exists(path, "Config file")
    suffix = path.suffix.lower()

    if suffix in [".yaml", ".yml"]:
        try:
            import yaml  # type: ignore
        except ImportError as e:
            raise RuntimeError(
                f"PyYAML is required to read YAML configs but is not installed.\n"
                f"Install it with: pip install pyyaml\n"
                f"Or use a JSON config instead.\n"
                f"Config path: {path}"
            ) from e
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
        if not isinstance(data, dict):
            raise ValueError(f"YAML config must be a mapping/dict: {path}")
        return data

    if suffix == ".json":
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        if not isinstance(data, dict):
            raise ValueError(f"JSON config must be an object/dict: {path}")
        return data

    raise ValueError(f"Unsupported config file type: {path} (use .yml/.yaml or .json)")


def load_combined_config(project_root: Path, marker: str, config_path: Optional[str]) -> Dict[str, Any]:
    """
    Load defaults + marker defaults + optional user config.
    Precedence (later overrides earlier):
      defaults.yml -> markers/{marker}.yml -> user --config
    """
    cfg: Dict[str, Any] = {}

    defaults = project_root / "config" / "defaults.yml"
    marker_cfg = project_root / "config" / "markers" / f"{marker}.yml"

    if defaults.exists():
        deep_update(cfg, load_config_file(defaults))
    if marker_cfg.exists():
        deep_update(cfg, load_config_file(marker_cfg))
    if config_path:
        user_cfg = Path(config_path)
        if not user_cfg.is_absolute():
            user_cfg = (project_root / user_cfg).resolve()
        deep_update(cfg, load_config_file(user_cfg))

    return cfg


@dataclass(frozen=True)
class Context:
    project_root: Path
    marker: str
    dataset: str
    metadata: Path
    config: Dict[str, Any]

    qiime2_root: Path
    results_root: Path
    logs_root: Path

    threads: int
    verbose: bool
    dry_run: bool
    force: bool

    def scope(self) -> str:
        """Return 'marker/dataset' scope string used in log messages."""
        return f"{self.marker}/{self.dataset}"

    def qdir(self, *parts: str) -> Path:
        """Return a path under the QIIME2 artifact root for this marker/dataset."""
        return self.qiime2_root / self.marker / self.dataset / Path(*parts)

    def rdir(self, *parts: str) -> Path:
        """Return a path under the results root for this marker/dataset."""
        return self.results_root / self.marker / self.dataset / Path(*parts)

    def ldir(self, *parts: str) -> Path:
        """Return a path under the logs root for this marker/dataset."""
        return self.logs_root / self.marker / self.dataset / Path(*parts)

    def log_file(self, subcmd: str) -> Path:
        """Return a timestamped log file path for a subcommand."""
        stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        return self.ldir(f"{subcmd}_{stamp}.log")
    def commands_file(self) -> Path:
        """Return the path to the running commands shell log for this scope."""
        return self.rdir("commands.sh")

    def get_cfg(self, *keys: str, default=None):
        """
Retrieve a nested value from the YAML config by a chain of keys.

        Returns default if any key in the chain is missing or if the value at
        that level is not a dict. Example: ctx.get_cfg('dada2', 'trunc_len_f').
        """
        cur: Any = self.config
        for k in keys:
            if not isinstance(cur, dict) or k not in cur:
                return default
            cur = cur[k]
        return cur


def cmd_smoke_test(args: argparse.Namespace, ctx: Context) -> None:
    """
Validate the environment and expected input artifacts without running any analysis.

    Checks that the QIIME2 conda environment is active, that the metadata file
    exists and is readable, and that key expected artifacts are present for the
    given marker/dataset scope. Use this after setup to confirm the pipeline is
    ready to run.
    """
    init_layout(ctx)

    # Check QIIME
    qiime = which_or_die("qiime")
    print(f"[smoke-test] qiime: {qiime}")

    # Check metadata exists
    require_exists(ctx.metadata, "Metadata")
    print(f"[smoke-test] metadata: OK ({ctx.metadata})")

    # If manifest provided, check it
    if args.manifest:
        m = Path(args.manifest)
        if not m.is_absolute():
            m = (ctx.project_root / m).resolve()
        require_exists(m, "Manifest")
        print(f"[smoke-test] manifest: OK ({m})")

    # Check expected artifacts if requested
    if args.check_artifacts:
        for p, label in [
            (ctx.qdir("imported", "demux.qza"), "demux.qza"),
            (ctx.qdir("dada2", "table.qza"), "table.qza"),
            (ctx.qdir("dada2", "rep-seqs.qza"), "rep-seqs.qza"),
            (ctx.qdir("taxonomy", "taxonomy.qza"), "taxonomy.qza"),
        ]:
            if p.exists():
                print(f"[smoke-test] {label}: OK ({p})")
            else:
                print(f"[smoke-test] {label}: MISSING ({p})")

    print("[smoke-test] done")

def run_cmd(cmd: Sequence[str], ctx: Context, subcmd: str, cwd: Optional[Path] = None) -> None:
    """
    Run a shell command safely, streaming stdout/stderr to terminal and also to a log file.
    """
    log_path = ctx.log_file(subcmd)
    safe_mkdir(log_path.parent)

    # Pretty print
    if ctx.verbose or ctx.dry_run:
        eprint("\n$ " + " ".join(cmd))
        eprint(f"  (log: {log_path})")

    if ctx.dry_run:
        # Still write a log stub
        log_path.write_text("DRY RUN\n" + " ".join(cmd) + "\n")
        return

    # Open log file for append
    with log_path.open("a", encoding="utf-8") as lf:
        lf.write("\n$ " + " ".join(cmd) + "\n")
        lf.flush()

        proc = subprocess.run(
            list(cmd),
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        lf.write(proc.stdout)
        lf.flush()

    # Append exact command to a reproducible script log
    cmd_sh = ctx.commands_file()
    safe_mkdir(cmd_sh.parent)
    with cmd_sh.open("a", encoding="utf-8") as cf:
        cf.write(" ".join(cmd) + "\n")

        # Stream to terminal as well
        sys.stdout.write(proc.stdout)
        sys.stdout.flush()

        if proc.returncode != 0:
            raise RuntimeError(f"Command failed (exit={proc.returncode}). See log: {log_path}")


def maybe_overwrite(path: Path, ctx: Context) -> None:
    """
    QIIME2 will often refuse to overwrite outputs. If --force is set and output exists, remove it.
    """
    if path.exists():
        if not ctx.force:
            raise FileExistsError(
                f"Output already exists: {path}\n"
                f"Use --force to overwrite, or choose a different dataset slug."
            )
        # QIIME2 outputs are usually files; sometimes directories (core-metrics output)
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def write_json(path: Path, obj) -> None:
    """Serialise obj to JSON and write to path, creating parent directories as needed."""
    safe_mkdir(path.parent)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


# ---------------------------
# Path conventions
# ---------------------------

def init_layout(ctx: Context) -> None:
    """
    Create the directory structure for this marker/dataset scope.
    """
    # qiime2 core
    for p in [
        ctx.qdir("imported"),
        ctx.qdir("dada2"),
        ctx.qdir("taxonomy"),
        ctx.qdir("diversity"),
    ]:
        safe_mkdir(p)

    # results
    for p in [
        ctx.rdir("tables"),
        ctx.rdir("taxonomy"),
        ctx.rdir("diversity"),
        ctx.rdir("figures"),
        ctx.rdir("differential"),
        ctx.rdir("exports"),
        ctx.rdir("bundles"),
    ]:
        safe_mkdir(p)

    # logs
    safe_mkdir(ctx.ldir())

    # record a small manifest of this run scope
    write_json(
        ctx.rdir("run_scope.json"),
        {
            "marker": ctx.marker,
            "dataset": ctx.dataset,
            "metadata": str(ctx.metadata),
            "created": dt.datetime.now().isoformat(),
        },
    )


# ---------------------------
# Subcommand implementations
# ---------------------------

def cmd_init(args: argparse.Namespace, ctx: Context) -> None:
    """Create the directory layout for this marker/dataset scope."""
    init_layout(ctx)
    print(f"[init] Created layout for {ctx.scope()}")

def cmd_import(args: argparse.Namespace, ctx: Context) -> None:
    """
Import paired-end reads from a manifest TSV into a QIIME2 demux artifact.

    Reads the manifest at --manifest (or the default imported/manifest.tsv),
    runs qiime tools import, and writes demux.qza + demux.qzv to the imported/
    subdirectory. The QZV can be inspected at view.qiime2.org to check
    read-quality profiles before choosing DADA2 truncation lengths.
    """
    init_layout(ctx)
    which_or_die("qiime")

    manifest = Path(args.manifest)
    require_exists(manifest, "Manifest file")

    out_demux_qza = ctx.qdir("imported", "demux.qza")
    out_demux_qzv = ctx.qdir("imported", "demux.qzv")
    maybe_overwrite(out_demux_qza, ctx)
    maybe_overwrite(out_demux_qzv, ctx)

    # You can adjust these import parameters to match your actual data type.
    # This is the common paired-end EMP format pattern; you may need different type/format.
    import_type = args.type
    import_format = args.format

    cmd = [
        "qiime", "tools", "import",
        "--type", import_type,
        "--input-path", str(manifest),
        "--output-path", str(out_demux_qza),
        "--input-format", import_format,
    ]
    run_cmd(cmd, ctx, "import")

    # Summarize to qzv (for demux viz)
    summarize_cmd = [
        "qiime", "demux", "summarize",
        "--i-data", str(out_demux_qza),
        "--o-visualization", str(out_demux_qzv),
    ]
    run_cmd(summarize_cmd, ctx, "import")


def cmd_cutadapt(args: argparse.Namespace, ctx: Context) -> None:
    """
    Trim primers from a demux artifact using qiime cutadapt trim-paired
    (or trim-single for single-end data).

    Input:  qiime2/{marker}/{dataset}/imported/demux.qza
    Output: qiime2/{marker}/{dataset}/imported/demux_trimmed.qza
            qiime2/{marker}/{dataset}/imported/demux_trimmed.qzv
    """
    init_layout(ctx)
    which_or_die("qiime")

    demux_in = Path(args.demux) if args.demux else ctx.qdir("imported", "demux.qza")
    require_exists(demux_in, "Demux artifact (--demux or default imported/demux.qza)")

    demux_trimmed = ctx.qdir("imported", "demux_trimmed.qza")
    demux_trimmed_qzv = ctx.qdir("imported", "demux_trimmed.qzv")
    for p in [demux_trimmed, demux_trimmed_qzv]:
        maybe_overwrite(p, ctx)

    single_end = getattr(args, "single_end", False)
    subcommand = "trim-single" if single_end else "trim-paired"

    cmd = [
        "qiime", "cutadapt", subcommand,
        "--i-demultiplexed-sequences", str(demux_in),
        "--o-trimmed-sequences", str(demux_trimmed),
        "--p-cores", str(ctx.threads),
        "--verbose",
    ]

    if single_end:
        if args.front:
            cmd += ["--p-front", args.front]
        if args.adapter:
            cmd += ["--p-adapter", args.adapter]
    else:
        if args.front_f:
            cmd += ["--p-front-f", args.front_f]
        if args.front_r:
            cmd += ["--p-front-r", args.front_r]
        if args.adapter_f:
            cmd += ["--p-adapter-f", args.adapter_f]
        if args.adapter_r:
            cmd += ["--p-adapter-r", args.adapter_r]

    if args.error_rate is not None:
        cmd += ["--p-error-rate", str(args.error_rate)]
    if args.minimum_length is not None:
        cmd += ["--p-minimum-length", str(args.minimum_length)]
    if args.discard_untrimmed:
        cmd += ["--p-discard-untrimmed"]

    run_cmd(cmd, ctx, "cutadapt")

    # Summarize trimmed demux
    viz_cmd = [
        "qiime", "demux", "summarize",
        "--i-data", str(demux_trimmed),
        "--o-visualization", str(demux_trimmed_qzv),
    ]
    run_cmd(viz_cmd, ctx, "cutadapt")
    log.info("Trimmed demux: %s", demux_trimmed)
    log.info("Trimmed summary: %s", demux_trimmed_qzv)


def cmd_dada2(args: argparse.Namespace, ctx: Context) -> None:
    """
Denoise paired-end reads with DADA2.

    Reads demux.qza (or --demux override), applies primer trimming
    (--trim-left-f / --trim-left-r) and quality truncation (--trunc-len-f /
    --trunc-len-r), and writes table.qza, rep-seqs.qza, and
    denoising-stats.qzv to the dada2/ subdirectory. Use primer_advisor.py
    to determine appropriate trim/trunc values from the demux QZV.
    """
    init_layout(ctx)
    which_or_die("qiime")

    demux_qza = Path(args.demux) if args.demux else ctx.qdir("imported", "demux.qza")
    require_exists(demux_qza, "Demux artifact (--demux or default imported/demux.qza)")

    out_table = ctx.qdir("dada2", "table.qza")
    out_repseqs = ctx.qdir("dada2", "rep-seqs.qza")
    out_stats = ctx.qdir("dada2", "denoising-stats.qza")
    out_stats_qzv = ctx.qdir("dada2", "denoising-stats.qzv")
    for p in [out_table, out_repseqs, out_stats, out_stats_qzv]:
        maybe_overwrite(p, ctx)

    # Paired-end DADA2 (adjust if single-end)
    cmd = [
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs", str(demux_qza),
        "--p-trim-left-f", str(args.trim_left_f),
        "--p-trim-left-r", str(args.trim_left_r),
        "--p-trunc-len-f", str(args.trunc_len_f),
        "--p-trunc-len-r", str(args.trunc_len_r),
        "--p-max-ee-f", str(args.max_ee_f),
        "--p-max-ee-r", str(args.max_ee_r),
        "--p-n-threads", str(ctx.threads),
        "--o-table", str(out_table),
        "--o-representative-sequences", str(out_repseqs),
        "--o-denoising-stats", str(out_stats),
    ]
    run_cmd(cmd, ctx, "dada2")

    viz_cmd = [
        "qiime", "metadata", "tabulate",
        "--m-input-file", str(out_stats),
        "--o-visualization", str(out_stats_qzv),
    ]
    run_cmd(viz_cmd, ctx, "dada2")


def cmd_taxonomy(args: argparse.Namespace, ctx: Context) -> None:
    """
Assign taxonomy to representative sequences using a trained classifier.

    Runs qiime feature-classifier classify-sklearn on rep-seqs.qza (or
    --repseqs override) using the classifier at --classifier. Writes
    taxonomy.qza and taxonomy.qzv to the taxonomy/ subdirectory. Marker-aware
    filter defaults are applied automatically; override with --include /
    --exclude. Use 00_build_classifiers.py to build the classifier QZA.
    """
    init_layout(ctx)
    which_or_die("qiime")

    classifier = Path(args.classifier)
    require_exists(classifier, "Classifier artifact (--classifier)")

    repseqs = Path(args.repseqs) if args.repseqs else ctx.qdir("dada2", "rep-seqs.qza")
    require_exists(repseqs, "Rep seqs artifact (--repseqs or default dada2/rep-seqs.qza)")

    if getattr(args, "exclude_controls", False):
        # Filter controls from the table so barplots & downstream collapse exclude them
        table_path = ctx.qdir("dada2", "table.qza")
        if table_path.exists():
            cmd_filter_controls(argparse.Namespace(table=str(table_path), sampletype_col="SampleType"), ctx)
            log.info("[taxonomy] Control-excluded table written to dada2/table_no_controls.qza")

    out_tax_qza = ctx.qdir("taxonomy", "taxonomy.qza")
    out_tax_qzv = ctx.qdir("taxonomy", "taxonomy.qzv")
    for p in [out_tax_qza, out_tax_qzv]:
        maybe_overwrite(p, ctx)

    cmd = [
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-classifier", str(classifier),
        "--i-reads", str(repseqs),
        "--o-classification", str(out_tax_qza),
    ]
    run_cmd(cmd, ctx, "taxonomy")

    viz_cmd = [
        "qiime", "metadata", "tabulate",
        "--m-input-file", str(out_tax_qza),
        "--o-visualization", str(out_tax_qzv),
    ]
    run_cmd(viz_cmd, ctx, "taxonomy")


def cmd_filter(args: argparse.Namespace, ctx: Context) -> None:
    """
    Filter table to exclude certain taxa (mitochondria/chloroplast/Euk/Archaea etc).
    """
    init_layout(ctx)
    which_or_die("qiime")

    table = Path(args.table) if args.table else ctx.qdir("dada2", "table.qza")
    require_exists(table, "Feature table (--table or default dada2/table.qza)")

    taxonomy = Path(args.taxonomy) if args.taxonomy else ctx.qdir("taxonomy", "taxonomy.qza")
    require_exists(taxonomy, "Taxonomy (--taxonomy or default taxonomy/taxonomy.qza)")

    exclude = args.exclude
    out_table = ctx.rdir("tables", f"table_filtered_excl_{exclude.replace(',', '-')}.qza")
    maybe_overwrite(out_table, ctx)

    cmd = [
        "qiime", "taxa", "filter-table",
        "--i-table", str(table),
        "--i-taxonomy", str(taxonomy),
        "--p-exclude", exclude,
        "--o-filtered-table", str(out_table),
    ]
    run_cmd(cmd, ctx, "filter")


def cmd_filter_controls(args: argparse.Namespace, ctx: Context) -> None:
    """
    Remove negative controls (NTC-*), positive controls (PAC-*), blanks (XB-*),
    and any other samples whose SampleType column equals 'Control' from the
    feature table.  Saves the result to dada2/table_no_controls.qza.

    Uses qiime feature-table filter-samples with a SQLite WHERE clause so the
    filtering is driven entirely by the QIIME2 metadata — no hard-coded sample
    IDs needed.

    The SampleType column must be present in the metadata TSV and set to
    'Control' for every control sample (02_make_qiime_metadata.py does this
    automatically for samples whose IDs match --control-prefixes).
    """
    init_layout(ctx)
    which_or_die("qiime")

    table = Path(args.table) if args.table else ctx.qdir("dada2", "table.qza")
    require_exists(table, "Feature table (--table or default dada2/table.qza)")

    out_table = ctx.qdir("dada2", "table_no_controls.qza")
    maybe_overwrite(out_table, ctx)

    sampletype_col = args.sampletype_col if hasattr(args, "sampletype_col") else "SampleType"
    where_clause = f'[{sampletype_col}] NOT IN (\'Control\')'

    cmd = [
        "qiime", "feature-table", "filter-samples",
        "--i-table", str(table),
        "--m-metadata-file", str(ctx.metadata),
        "--p-where", where_clause,
        "--o-filtered-table", str(out_table),
    ]
    run_cmd(cmd, ctx, "filter_controls")
    log.info("Controls excluded. Filtered table: %s", out_table)


def cmd_collapse(args: argparse.Namespace, ctx: Context) -> None:
    """
Collapse a feature table to a specified taxonomic level.

    Runs qiime taxa collapse on the filtered feature table and taxonomy,
    grouping ASVs at the requested level (e.g. 6 = genus for SILVA, 7 for
    species). The collapsed table is written to tables/ and used as input
    for differential abundance and taxonomy plotting steps.
    """
    init_layout(ctx)
    which_or_die("qiime")

    level = int(args.level)

    table = Path(args.table) if args.table else ctx.qdir("dada2", "table.qza")
    require_exists(table, "Feature table (--table or default dada2/table.qza)")

    taxonomy = Path(args.taxonomy) if args.taxonomy else ctx.qdir("taxonomy", "taxonomy.qza")
    require_exists(taxonomy, "Taxonomy (--taxonomy or default taxonomy/taxonomy.qza)")

    out_table = ctx.rdir("tables", f"table_L{level}.qza")
    maybe_overwrite(out_table, ctx)

    cmd = [
        "qiime", "taxa", "collapse",
        "--i-table", str(table),
        "--i-taxonomy", str(taxonomy),
        "--p-level", str(level),
        "--o-collapsed-table", str(out_table),
    ]
    run_cmd(cmd, ctx, "collapse")


def cmd_diversity(args: argparse.Namespace, ctx: Context) -> None:
    """
    Core metrics + UniFrac require a rooted tree.
    """
    init_layout(ctx)
    which_or_die("qiime")

    table = Path(args.table) if args.table else ctx.qdir("dada2", "table.qza")
    require_exists(table, "Feature table (--table or default dada2/table.qza)")

    if getattr(args, "exclude_controls", False):
        cmd_filter_controls(argparse.Namespace(table=str(table), sampletype_col="SampleType"), ctx)
        table = ctx.qdir("dada2", "table_no_controls.qza")
        log.info("[diversity] Using control-excluded table: %s", table)

    rooted_tree = Path(args.tree)
    require_exists(rooted_tree, "Rooted tree (--tree)")

    sampling_depth = int(args.sampling_depth)

    out_dir = ctx.qdir("diversity", f"core_metrics_depth{sampling_depth}")
    maybe_overwrite(out_dir, ctx)

    cmd = [
        "qiime", "diversity", "core-metrics-phylogenetic",
        "--i-phylogeny", str(rooted_tree),
        "--i-table", str(table),
        "--p-sampling-depth", str(sampling_depth),
        "--m-metadata-file", str(ctx.metadata),
        "--output-dir", str(out_dir),
    ]
    run_cmd(cmd, ctx, "diversity")

    # Optional group significance: PERMANOVA/PERMDISP on distance matrices
    # You can run these later via separate commands if you prefer.
    if args.group_column:
        group = args.group_column
        for metric in ["weighted_unifrac", "unweighted_unifrac", "bray_curtis", "jaccard"]:
            dm = out_dir / f"{metric}_distance_matrix.qza"
            if not dm.exists():
                continue

            perma_out = ctx.rdir("diversity", f"{metric}_permanova_{group}.qzv")
            disp_out = ctx.rdir("diversity", f"{metric}_permdisp_{group}.qzv")
            maybe_overwrite(perma_out, ctx)
            maybe_overwrite(disp_out, ctx)

            perma_cmd = [
                "qiime", "diversity", "beta-group-significance",
                "--i-distance-matrix", str(dm),
                "--m-metadata-file", str(ctx.metadata),
                "--m-metadata-column", group,
                "--p-method", "permanova",
                "--p-pairwise",
                "--o-visualization", str(perma_out),
            ]
            run_cmd(perma_cmd, ctx, "diversity")

            disp_cmd = [
                "qiime", "diversity", "beta-group-significance",
                "--i-distance-matrix", str(dm),
                "--m-metadata-file", str(ctx.metadata),
                "--m-metadata-column", group,
                "--p-method", "permdisp",
                "--p-pairwise",
                "--o-visualization", str(disp_out),
            ]
            run_cmd(disp_cmd, ctx, "diversity")



def cmd_diversity_nonphylo(args, ctx):
    """
    Run core-metrics (non-phylogenetic) for markers without a tree.

    Produces: rarefied table, observed_features, shannon, evenness vectors,
    Bray-Curtis and Jaccard distance matrices + Emperor PCoA plots.

    Use this for MiFish (12S), cytb, and any other marker where a
    reliable phylogenetic tree is not available.

    Example:
      diversity-nonphylo --sampling-depth 17000 --group-column Group
    """
    init_layout(ctx)
    which_or_die("qiime")

    table = Path(args.table) if args.table else ctx.qdir("dada2", "table_filtered.qza")
    if not table.exists():
        table = ctx.qdir("dada2", "table.qza")
    require_exists(table, "Feature table (--table or default dada2/table_filtered.qza)")

    if getattr(args, "exclude_controls", False):
        cmd_filter_controls(argparse.Namespace(table=str(table), sampletype_col="SampleType"), ctx)
        table = ctx.qdir("dada2", "table_no_controls.qza")
        log.info("[diversity-nonphylo] Using control-excluded table: %s", table)

    sampling_depth = int(args.sampling_depth)
    out_dir = ctx.qdir("diversity", f"core-metrics-{sampling_depth}")
    maybe_overwrite(out_dir, ctx)

    cmd = [
        "qiime", "diversity", "core-metrics",
        "--i-table", str(table),
        "--p-sampling-depth", str(sampling_depth),
        "--m-metadata-file", str(ctx.metadata),
        "--output-dir", str(out_dir),
    ]
    run_cmd(cmd, ctx, "diversity-nonphylo")

    if args.group_column:
        _run_stats(ctx, out_dir, args.group_column,
                   phylo=False, subcmd="diversity-nonphylo")


def cmd_stats(args, ctx):
    """
    Run alpha-group-significance (Kruskal-Wallis) and beta-group-significance
    (PERMANOVA + PERMDISP) on an existing core-metrics output directory.

    Works with both phylogenetic and non-phylogenetic core-metrics outputs.
    Outputs go to results/{marker}/{dataset}/diversity/.

    Example:
      stats --metrics-dir qiime2/MiFish/all/diversity/core-metrics-17000
            --group-column Group --phylo false
    """
    init_layout(ctx)
    which_or_die("qiime")

    metrics_dir = Path(args.metrics_dir)
    require_exists(metrics_dir, "Core-metrics output directory (--metrics-dir)")

    phylo = args.phylo.lower() not in ("false", "no", "0", "f")
    _run_stats(ctx, metrics_dir, args.group_column, phylo=phylo, subcmd="stats")


def _run_stats(ctx, metrics_dir, group_column, phylo=True, subcmd="stats"):
    """Internal: run alpha + beta group significance on a core-metrics dir."""

    # Alpha diversity
    alpha_metrics = ["observed_features", "shannon", "evenness"]
    if phylo:
        alpha_metrics.append("faith_pd")

    for metric in alpha_metrics:
        vec = metrics_dir / f"{metric}_vector.qza"
        if not vec.exists():
            log.warning("[%s] Skipping alpha metric (not found): %s", subcmd, vec)
            continue
        out = ctx.rdir("diversity", f"{metric}_group_sig_{group_column}.qzv")
        maybe_overwrite(out, ctx)
        cmd = [
            "qiime", "diversity", "alpha-group-significance",
            "--i-alpha-diversity", str(vec),
            "--m-metadata-file", str(ctx.metadata),
            "--o-visualization", str(out),
        ]
        run_cmd(cmd, ctx, subcmd)
        log.info("[%s] Alpha significance saved: %s", subcmd, out)

    # Beta diversity
    beta_metrics = ["bray_curtis", "jaccard"]
    if phylo:
        beta_metrics += ["weighted_unifrac", "unweighted_unifrac"]

    for metric in beta_metrics:
        dm = metrics_dir / f"{metric}_distance_matrix.qza"
        if not dm.exists():
            log.warning("[%s] Skipping beta metric (not found): %s", subcmd, dm)
            continue
        for method in ("permanova", "permdisp"):
            out = ctx.rdir("diversity", f"{metric}_{method}_{group_column}.qzv")
            maybe_overwrite(out, ctx)
            cmd = [
                "qiime", "diversity", "beta-group-significance",
                "--i-distance-matrix", str(dm),
                "--m-metadata-file", str(ctx.metadata),
                "--m-metadata-column", group_column,
                "--p-method", method,
                "--p-pairwise",
                "--o-visualization", str(out),
            ]
            run_cmd(cmd, ctx, subcmd)
            log.info("[%s] Beta %s saved: %s", subcmd, method, out)


def cmd_barplot(args: argparse.Namespace, ctx: Context) -> None:
    """
    QIIME taxa barplot visualization.
    IMPORTANT: This expects that table Feature IDs match taxonomy Feature IDs.
    If you pass a collapsed table created via taxa collapse, you should NOT pass taxonomy.qza that is ASV-level;
    instead, pass taxonomy consistent with the table feature ids (or just barplot the original table).
    """
    init_layout(ctx)
    which_or_die("qiime")

    table = Path(args.table)
    require_exists(table, "--table")

    taxonomy = Path(args.taxonomy)
    require_exists(taxonomy, "--taxonomy")

    out_qzv = ctx.rdir("taxonomy", args.out_name)
    maybe_overwrite(out_qzv, ctx)

    cmd = [
        "qiime", "taxa", "barplot",
        "--i-table", str(table),
        "--i-taxonomy", str(taxonomy),
        "--m-metadata-file", str(ctx.metadata),
        "--o-visualization", str(out_qzv),
    ]
    run_cmd(cmd, ctx, "barplot")

def cmd_barplot_asv(args: argparse.Namespace, ctx: Context) -> None:
    """
    Convenience taxa barplot using the ASV table + matching taxonomy.

    Purpose:
    - Canonical interactive taxa barplot that lets you toggle taxonomic levels (phylum→genus)
    - Avoids the feature-ID mismatch that occurs after taxa collapse (where feature IDs become taxonomy strings)
    """
    init_layout(ctx)
    which_or_die("qiime")

    table = ctx.qdir("dada2", "table.qza")
    taxonomy = ctx.qdir("taxonomy", "taxonomy.qza")
    require_exists(table, "ASV table (expected at qiime2/{marker}/{dataset}/dada2/table.qza)")
    require_exists(taxonomy, "Taxonomy (expected at qiime2/{marker}/{dataset}/taxonomy/taxonomy.qza)")

    out_qzv = ctx.rdir("taxonomy", args.out_name)
    maybe_overwrite(out_qzv, ctx)

    cmd = [
        "qiime", "taxa", "barplot",
        "--i-table", str(table),
        "--i-taxonomy", str(taxonomy),
        "--m-metadata-file", str(ctx.metadata),
        "--o-visualization", str(out_qzv),
    ]
    run_cmd(cmd, ctx, "barplot_asv")

def cmd_taxa_plots(args: argparse.Namespace, ctx: Context) -> None:
    """
    Make common taxonomy visualizations.

    Profiles:
      - asv:      ASV table + taxonomy barplot (canonical)
      - filtered: filtered ASV table (exclude taxa) + barplot
      - both:     asv + filtered
    """
    init_layout(ctx)
    which_or_die("qiime")

    profile = args.profile

    # 1) ASV plot
    if profile in ("asv", "both"):
        cmd_barplot_asv(argparse.Namespace(out_name=args.asv_out_name), ctx)

    # 2) Filtered plot (still ASV IDs, so taxonomy matches!)
    if profile in ("filtered", "both"):
        in_table = ctx.qdir("dada2", "table.qza")
        in_tax = ctx.qdir("taxonomy", "taxonomy.qza")
        require_exists(in_table, "ASV table")
        require_exists(in_tax, "Taxonomy")

        exclude = args.exclude
        filtered_table = ctx.rdir("tables", f"table_filtered_excl_{exclude.replace(',', '-')}.qza")
        maybe_overwrite(filtered_table, ctx)

        filter_cmd = [
            "qiime", "taxa", "filter-table",
            "--i-table", str(in_table),
            "--i-taxonomy", str(in_tax),
            "--p-exclude", exclude,
            "--o-filtered-table", str(filtered_table),
        ]
        run_cmd(filter_cmd, ctx, "taxa_plots")

        out_qzv = ctx.rdir("taxonomy", args.filtered_out_name)
        maybe_overwrite(out_qzv, ctx)

        bar_cmd = [
            "qiime", "taxa", "barplot",
            "--i-table", str(filtered_table),
            "--i-taxonomy", str(in_tax),
            "--m-metadata-file", str(ctx.metadata),
            "--o-visualization", str(out_qzv),
        ]
        run_cmd(bar_cmd, ctx, "taxa_plots")

        # Optional: collapse to a chosen level for downstream ANCOM (produces qza table)
        if args.collapse_level is not None:
            level = int(args.collapse_level)
            collapsed = ctx.rdir("tables", f"table_filtered_L{level}.qza")
            maybe_overwrite(collapsed, ctx)
            col_cmd = [
                "qiime", "taxa", "collapse",
                "--i-table", str(filtered_table),
                "--i-taxonomy", str(in_tax),
                "--p-level", str(level),
                "--o-collapsed-table", str(collapsed),
            ]
            run_cmd(col_cmd, ctx, "taxa_plots")

def cmd_diff_ancom(args: argparse.Namespace, ctx: Context) -> None:
    """
    ANCOM expects FeatureTable[Composition] input.
    """
    init_layout(ctx)
    which_or_die("qiime")

    table = Path(args.table)
    require_exists(table, "--table")

    group_col = args.group_column
    min_frequency = int(args.min_frequency)

    # Filter rare features (optional but common)
    filtered = ctx.rdir("tables", f"ancom_minfreq{min_frequency}.qza")
    pseudo = ctx.rdir("tables", f"ancom_minfreq{min_frequency}_pseudo.qza")
    out_qzv = ctx.rdir("differential", args.out_name)

    for p in [filtered, pseudo, out_qzv]:
        maybe_overwrite(p, ctx)

    filter_cmd = [
        "qiime", "feature-table", "filter-features",
        "--i-table", str(table),
        "--p-min-frequency", str(min_frequency),
        "--o-filtered-table", str(filtered),
    ]
    run_cmd(filter_cmd, ctx, "diff")

    pseudo_cmd = [
        "qiime", "composition", "add-pseudocount",
        "--i-table", str(filtered),
        "--o-composition-table", str(pseudo),
    ]
    run_cmd(pseudo_cmd, ctx, "diff")

    ancom_cmd = [
        "qiime", "composition", "ancom",
        "--i-table", str(pseudo),
        "--m-metadata-file", str(ctx.metadata),
        "--m-metadata-column", group_col,
        "--o-visualization", str(out_qzv),
    ]
    run_cmd(ancom_cmd, ctx, "diff")


def cmd_export(args: argparse.Namespace, ctx: Context) -> None:
    """
Export any QIIME2 artifact to a plain-file format using qiime tools export.

    Writes output to results/{marker}/{dataset}/exports/{out_dir_name}/. For
    feature tables this produces a BIOM file; for most other artifacts it
    produces the native file format (TSV, FASTA, newick, etc.).
    """
    init_layout(ctx)
    which_or_die("qiime")

    artifact = Path(args.artifact)
    require_exists(artifact, "--artifact")

    out_dir = ctx.rdir("exports", args.out_dir_name)
    maybe_overwrite(out_dir, ctx)

    cmd = [
        "qiime", "tools", "export",
        "--input-path", str(artifact),
        "--output-path", str(out_dir),
    ]
    run_cmd(cmd, ctx, "export")


def cmd_bundle(args: argparse.Namespace, ctx: Context) -> None:
    """
Bundle metadata and key results artifacts into a dated tar.gz archive.

    Collects the metadata TSV, taxonomy QZAs, feature tables, and diversity
    outputs for this marker/dataset and writes a compressed archive to
    results/{marker}/{dataset}/bundles/. Useful for sharing results or
    creating reproducible snapshots before a major analysis step.
    """
    init_layout(ctx)

    stamp = dt.date.today().isoformat()
    out_tar = ctx.rdir("bundles", f"{ctx.marker}_{ctx.dataset}_bundle_{stamp}.tar.gz")
    maybe_overwrite(out_tar, ctx)

    # Build include list
    include: List[Path] = []
    if args.include_metadata:
        include.append(ctx.metadata)
    if args.include_qiime2:
        include.append(ctx.qiime2_root / ctx.marker / ctx.dataset)
    if args.include_results:
        include.append(ctx.results_root / ctx.marker / ctx.dataset)

    # Convert to project-root-relative paths (tar wants that)
    rel_paths: List[str] = []
    for p in include:
        require_exists(p, "Bundle input")
        rel_paths.append(str(p.relative_to(ctx.project_root)))

    cmd = ["tar", "-czf", str(out_tar.relative_to(ctx.project_root)), *rel_paths]
    run_cmd(cmd, ctx, "bundle", cwd=ctx.project_root)

    print(f"[bundle] Wrote: {out_tar}")

def cmd_run(args: argparse.Namespace, ctx: Context) -> None:
    """
    Convenience: run a common end-to-end workflow.
    init -> import -> dada2 -> taxonomy -> diversity -> barplot-asv
    """
    init_layout(ctx)

    # Import
    cmd_import(argparse.Namespace(
        manifest=args.manifest,
        type=args.type,
        format=args.format,
    ), ctx)

    # DADA2
    cmd_dada2(argparse.Namespace(
        demux=None,
        trim_left_f=args.trim_left_f,
        trim_left_r=args.trim_left_r,
        trunc_len_f=args.trunc_len_f,
        trunc_len_r=args.trunc_len_r,
        max_ee_f=args.max_ee_f,
        max_ee_r=args.max_ee_r,
    ), ctx)

    # Taxonomy
    cmd_taxonomy(argparse.Namespace(
        classifier=args.classifier,
        repseqs=None,
    ), ctx)

    # Diversity (requires a rooted tree you provide)
    cmd_diversity(argparse.Namespace(
        sampling_depth=args.sampling_depth,
        tree=args.tree,
        table=None,
        group_column=args.group_column,
    ), ctx)

    # Canonical ASV barplot
    cmd_barplot_asv(argparse.Namespace(
        out_name=args.barplot_out_name,
    ), ctx)

# ---------------------------
# Argument parsing
# ---------------------------

def build_parser() -> argparse.ArgumentParser:
    """
Build and return the top-level argument parser for 03_run_full_metabarcoding_pipeline.py.

    The parser requires --marker, --dataset, and --metadata as global flags,
    then dispatches to a subcommand (init, import, dada2, taxonomy, collapse,
    diversity, diff-ancom, export, bundle, smoke-test, run, etc.). Each
    subcommand has its own set of optional arguments documented in --help.
    """
    p = argparse.ArgumentParser(prog="run_full_metabarcoding_pipeline.py", description="QIIME2 pipeline wrapper (marker + dataset scoped).", formatter_class=argparse.ArgumentDefaultsHelpFormatter,)

    # Global flags
    p.add_argument("--project-root", default=".", help="Project root directory.")
    p.add_argument("--marker", required=True, help="Marker name, e.g. 16S, 18S.")
    p.add_argument("--dataset", required=True, help="Dataset slug, e.g. all, DvT.")
    p.add_argument("--metadata", required=True, help="QIIME metadata TSV path.")

    p.add_argument("--qiime2-root", default="qiime2", help="Root for core QIIME artifacts.")
    p.add_argument("--results-root", default="results", help="Root for analysis outputs.")
    p.add_argument("--logs-root", default="logs", help="Root for logs.")

    p.add_argument("--threads", type=int, default=8, help="Threads/jobs for tools.")
    p.add_argument("--dry-run", action="store_true", help="Print commands without running.")
    p.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    p.add_argument("--verbose", action="store_true", help="Print commands as they run.")
    p.add_argument("--config", default=None, help="Optional YAML/JSON config to override defaults.")
    p.add_argument("--print-paths", action="store_true", help="Print computed output paths for this marker/dataset and exit.")

    sp = p.add_subparsers(dest="command", required=True)

    # init
    sp_init = sp.add_parser("init", help="Create directory layout for this marker/dataset.")
    sp_init.set_defaults(func=cmd_init)


    # import
    sp_imp = sp.add_parser("import", help="Import reads from a manifest into a demux artifact.")
    sp_imp.add_argument("--manifest", required=True, help="Manifest TSV path.")
    sp_imp.add_argument("--type", default="SampleData[PairedEndSequencesWithQuality]", help="QIIME2 import type.",)
    sp_imp.add_argument("--format", default="PairedEndFastqManifestPhred33V2", help="QIIME2 import format.",)

    sp_imp.set_defaults(func=cmd_import)

    # cutadapt
    sp_cut = sp.add_parser(
        "cutadapt",
        help="Trim primers from a demux artifact using qiime cutadapt trim-paired (or --single-end).",
    )
    sp_cut.add_argument("--demux", default=None,
                        help="Demux QZA to trim (defaults to imported/demux.qza).")
    sp_cut.add_argument("--single-end", action="store_true", default=False,
                        help="Use trim-single instead of trim-paired.")
    # Paired-end primer args
    sp_cut.add_argument("--front-f", default=None,
                        help="Forward primer sequence (5′→3′, linked or plain). Paired-end only.")
    sp_cut.add_argument("--front-r", default=None,
                        help="Reverse primer sequence (5′→3′). Paired-end only.")
    sp_cut.add_argument("--adapter-f", default=None,
                        help="3′ adapter on forward reads. Paired-end only.")
    sp_cut.add_argument("--adapter-r", default=None,
                        help="3′ adapter on reverse reads. Paired-end only.")
    # Single-end primer args
    sp_cut.add_argument("--front", default=None,
                        help="5′ adapter/primer sequence. Single-end only.")
    sp_cut.add_argument("--adapter", default=None,
                        help="3′ adapter sequence. Single-end only.")
    # Shared args
    sp_cut.add_argument("--error-rate", type=float, default=0.1,
                        help="Maximum error rate for primer matching. Default: 0.1.")
    sp_cut.add_argument("--minimum-length", type=int, default=50,
                        help="Discard reads shorter than this after trimming. Default: 50.")
    sp_cut.add_argument("--discard-untrimmed", action="store_true", default=False,
                        help="Discard reads where no primer was found.")
    sp_cut.set_defaults(func=cmd_cutadapt)
    sp_run = sp.add_parser("run", help="Run init→import→dada2→taxonomy→diversity→barplot-asv")
    sp_run.add_argument("--manifest", required=True, help="Manifest TSV path.")
    sp_run.add_argument("--classifier", required=True, help="Classifier QZA path.")
    sp_run.add_argument("--sampling-depth", required=True, help="Rarefaction depth integer.")
    sp_run.add_argument("--tree", required=True, help="Rooted tree QZA path for UniFrac.")
    sp_run.add_argument("--group-column", default="Group", help="Metadata column for group stats (PERMANOVA/PERMDISP).")

    sp_run.add_argument(
        "--type",
        default="SampleData[PairedEndSequencesWithQuality]",
        help="QIIME2 import type.",
    )
    sp_run.add_argument(
        "--format",
        default="PairedEndFastqManifestPhred33V2",
        help="QIIME2 import format.",
    )

    # DADA2 params (defaults may be overridden by config later if you want)
    sp_run.add_argument("--trim-left-f", type=int, default=0)
    sp_run.add_argument("--trim-left-r", type=int, default=0)
    sp_run.add_argument("--trunc-len-f", type=int, default=0)
    sp_run.add_argument("--trunc-len-r", type=int, default=0)
    sp_run.add_argument("--max-ee-f", type=float, default=2.0)
    sp_run.add_argument("--max-ee-r", type=float, default=2.0)

    sp_run.add_argument("--barplot-out-name", default="taxa_barplot_ASV.qzv")
    sp_run.set_defaults(func=cmd_run)

    # dada2
    sp_dada2 = sp.add_parser("dada2", help="Run DADA2 denoise-paired.")
    sp_dada2.add_argument("--demux", default=None, help="Demux QZA (defaults to qiime2/.../imported/demux.qza).")
    sp_dada2.add_argument("--trim-left-f", type=int, default=0)
    sp_dada2.add_argument("--trim-left-r", type=int, default=0)
    sp_dada2.add_argument("--trunc-len-f", type=int, default=0)
    sp_dada2.add_argument("--trunc-len-r", type=int, default=0)
    sp_dada2.add_argument("--max-ee-f", type=float, default=2.0)
    sp_dada2.add_argument("--max-ee-r", type=float, default=2.0)
    sp_dada2.set_defaults(func=cmd_dada2)

    # taxonomy
    sp_tax = sp.add_parser("taxonomy", help="Assign taxonomy with a sklearn classifier.")
    sp_tax.add_argument("--classifier", required=True, help="Classifier QZA.")
    sp_tax.add_argument("--repseqs", default=None, help="Rep seqs QZA (defaults to qiime2/.../dada2/rep-seqs.qza).")
    sp_tax.add_argument("--exclude-controls", action="store_true", default=False,
                        help="Also write a control-excluded table (dada2/table_no_controls.qza) for downstream use.")
    sp_tax.set_defaults(func=cmd_taxonomy)

    # filter
    sp_filt = sp.add_parser("filter", help="Filter table by excluding taxa using taxonomy.")
    sp_filt.add_argument("--table", default=None, help="Table QZA (defaults to qiime2/.../dada2/table.qza).")
    sp_filt.add_argument("--taxonomy", default=None, help="Taxonomy QZA (defaults to qiime2/.../taxonomy/taxonomy.qza).")
    sp_filt.add_argument("--exclude", default="mitochondria,chloroplast,Eukaryota,Archaea", help="Comma-separated taxa strings for --p-exclude.",)
    sp_filt.set_defaults(func=cmd_filter)

    # filter-controls
    sp_fc = sp.add_parser(
        "filter-controls",
        help="Remove NTC/PAC/XB controls from the feature table (saves dada2/table_no_controls.qza)."
    )
    sp_fc.add_argument("--table", default=None, help="Table QZA (defaults to dada2/table.qza).")
    sp_fc.add_argument("--sampletype-col", default="SampleType",
                       help="Metadata column whose value is 'Control' for control samples. Default: SampleType.")
    sp_fc.set_defaults(func=cmd_filter_controls)

    # collapse
    sp_col = sp.add_parser("collapse", help="Collapse a table to a taxonomic level.")
    sp_col.add_argument("--level", required=True, help="Taxonomic level integer (e.g., 6=genus, 5=family).")
    sp_col.add_argument("--table", default=None, help="Table QZA (defaults to qiime2/.../dada2/table.qza).")
    sp_col.add_argument("--taxonomy", default=None, help="Taxonomy QZA (defaults to qiime2/.../taxonomy/taxonomy.qza).")
    sp_col.set_defaults(func=cmd_collapse)

    # diversity
    sp_div = sp.add_parser("diversity", help="Run core-metrics-phylogenetic.")
    sp_div.add_argument("--sampling-depth", required=True, help="Rarefaction depth integer.")
    sp_div.add_argument("--tree", required=True, help="Rooted tree QZA.")
    sp_div.add_argument("--table", default=None, help="Table QZA (defaults to qiime2/.../dada2/table.qza).")
    sp_div.add_argument("--group-column", default=None, help="If set, also run beta-group-significance PERMANOVA/PERMDISP.")
    sp_div.add_argument("--exclude-controls", action="store_true", default=False,
                        help="Filter NTC/PAC/XB controls from table before diversity analysis.")
    sp_div.set_defaults(func=cmd_diversity)

    # diversity-nonphylo (for markers without a tree: MiFish, cytb, etc.)
    sp_divnp = sp.add_parser(
        "diversity-nonphylo",
        help="Run core-metrics (non-phylogenetic) for markers without a tree."
    )
    sp_divnp.add_argument("--sampling-depth", required=True, help="Rarefaction depth integer.")
    sp_divnp.add_argument("--table", default=None, help="Table QZA (defaults to dada2/table_filtered.qza, falls back to table.qza).")
    sp_divnp.add_argument("--group-column", default=None, help="If set, also run alpha/beta group significance tests.")
    sp_divnp.add_argument("--exclude-controls", action="store_true", default=False,
                          help="Filter NTC/PAC/XB controls from table before diversity analysis.")
    sp_divnp.set_defaults(func=cmd_diversity_nonphylo)

    # stats — alpha + beta significance on an existing core-metrics directory
    sp_stats = sp.add_parser(
        "stats",
        help="Run alpha-group-significance (KW) and beta-group-significance (PERMANOVA+PERMDISP)."
    )
    sp_stats.add_argument("--metrics-dir", required=True, help="Path to existing core-metrics output directory.")
    sp_stats.add_argument("--group-column", required=True, help="Metadata column to test (e.g. Group).")
    sp_stats.add_argument("--phylo", default="false", help="Include phylogenetic metrics (faith_pd, UniFrac). Default: false.")
    sp_stats.set_defaults(func=cmd_stats)


    # barplot
    sp_bar = sp.add_parser("barplot", help="Make a QIIME2 taxa barplot qzv.")
    sp_bar.add_argument("--table", required=True, help="FeatureTable QZA to plot.")
    sp_bar.add_argument("--taxonomy", required=True, help="Taxonomy QZA matching the table feature IDs.")
    sp_bar.add_argument("--out-name", default="taxa_barplot.qzv", help="Output qzv filename (in results/.../taxonomy/).")
    sp_bar.set_defaults(func=cmd_barplot)

    # taxa plots
    sp_tp = sp.add_parser("taxa-plots", help="Generate canonical ASV and/or filtered taxonomy barplots.")
    sp_tp.add_argument("--profile", choices=["asv", "filtered", "both"], default="both")
    sp_tp.add_argument("--exclude", default="mitochondria,chloroplast,Eukaryota,Archaea", help="Comma-separated taxa to exclude for filtered plot (16S defaults).",)
    sp_tp.add_argument("--asv-out-name", default="taxa_barplot_ASV.qzv")
    sp_tp.add_argument("--filtered-out-name", default="taxa_barplot_filtered_ASV.qzv")
    sp_tp.add_argument("--collapse-level", default=None, help="Optional collapse level (e.g., 5=family, 6=genus) for downstream tests.")
    sp_tp.set_defaults(func=cmd_taxa_plots)

    # diff (ANCOM for now)
    sp_diff = sp.add_parser("diff-ancom", help="Run ANCOM differential abundance.")
    sp_diff.add_argument("--table", required=True, help="Input table QZA (typically collapsed to genus/family).")
    sp_diff.add_argument("--group-column", default="Group", help="Metadata column for grouping.")
    sp_diff.add_argument("--min-frequency", type=int, default=10, help="Min feature frequency before ANCOM.")
    sp_diff.add_argument("--out-name", default="ancom.qzv", help="Output qzv filename (in results/.../differential/).")
    sp_diff.set_defaults(func=cmd_diff_ancom)

    # export
    sp_exp = sp.add_parser("export", help="qiime tools export for any artifact.")
    sp_exp.add_argument("--artifact", required=True, help="Artifact path (qza/qzv).")
    sp_exp.add_argument("--out-dir-name", default="exported_artifact", help="Directory name under results/.../exports/")
    sp_exp.set_defaults(func=cmd_export)

    # bundle
    sp_bun = sp.add_parser("bundle", help="Bundle metadata + qiime2 + results into a tar.gz.")
    sp_bun.add_argument("--include-metadata", action="store_true", default=True)
    sp_bun.add_argument("--include-qiime2", action="store_true", default=True)
    sp_bun.add_argument("--include-results", action="store_true", default=True)
    sp_bun.set_defaults(func=cmd_bundle)

    # barplot-asv (convenience)
    sp_ba = sp.add_parser("barplot-asv", help="Make a canonical taxa barplot from the ASV table + taxonomy (toggle levels in QIIME view).")
    sp_ba.add_argument("--out-name", default="taxa_barplot_ASV.qzv", help="Output qzv filename (in results/.../taxonomy/).")
    sp_ba.set_defaults(func=cmd_barplot_asv)

    # smoke test
    sp_st = sp.add_parser("smoke-test", help="Validate environment, metadata, and (optionally) expected artifacts.")
    sp_st.add_argument("--manifest", default=None, help="Optional manifest to validate existence.")
    sp_st.add_argument("--check-artifacts", action="store_true", help="Check for common expected QIIME artifacts.")
    sp_st.set_defaults(func=cmd_smoke_test)

    return p


def make_context(args: argparse.Namespace) -> Context:
    """
Build a Context dataclass from parsed command-line arguments.

    Resolves all paths relative to --project-root, loads an optional YAML
    config file, and applies any per-marker config overrides. The returned
    Context is passed to every subcommand handler.
    """
    project_root = Path(args.project_root).resolve()
    marker = args.marker
    dataset = args.dataset
    metadata = Path(args.metadata)

    # Interpret metadata relative to project root if not absolute
    if not metadata.is_absolute():
        metadata = (project_root / metadata).resolve()
    require_exists(metadata, "Metadata file (--metadata)")

    qiime2_root = (project_root / args.qiime2_root).resolve()
    results_root = (project_root / args.results_root).resolve()
    logs_root = (project_root / args.logs_root).resolve()

    cfg = load_combined_config(project_root, marker, args.config)

    return Context(
        project_root=project_root,
        marker=marker,
        dataset=dataset,
        metadata=metadata,
        qiime2_root=qiime2_root,
        results_root=results_root,
        logs_root=logs_root,
        threads=args.threads,
        verbose=args.verbose,
        dry_run=args.dry_run,
        force=args.force,
        config=cfg,
    )

def main(argv: Optional[Sequence[str]] = None) -> int:
    """
Entry point for 03_run_full_metabarcoding_pipeline.py.

    Parses global flags and the requested subcommand, builds a Context, and
    dispatches to the appropriate handler. Returns 0 on success, 1 on error.
    Add --dry-run to any invocation to print commands without executing them.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    ctx = make_context(args)

    if getattr(args, "print_paths", False):
        init_layout(ctx)
        print(f"Scope: {ctx.scope()}")
        print(f"Project root: {ctx.project_root}")
        print(f"QIIME2 root: {ctx.qiime2_root}")
        print(f"Results root: {ctx.results_root}")
        print(f"Logs root: {ctx.logs_root}")
        print(f"Commands file: {ctx.commands_file()}")
        print("\nKey dirs:")
        print(f"  qiime2 imported:  {ctx.qdir('imported')}")
        print(f"  qiime2 dada2:     {ctx.qdir('dada2')}")
        print(f"  qiime2 taxonomy:  {ctx.qdir('taxonomy')}")
        print(f"  qiime2 diversity: {ctx.qdir('diversity')}")
        print(f"  results figures:  {ctx.rdir('figures')}")
        print(f"  results tables:   {ctx.rdir('tables')}")
        print(f"  results taxonomy: {ctx.rdir('taxonomy')}")
        print(f"  results diversity:{ctx.rdir('diversity')}")
        return 0

    try:
        # Ensure we run from project_root for consistent relative paths/logs
        os.chdir(ctx.project_root)
        args.func(args, ctx)

        # ── Post-run guidance ────────────────────────────────────────────────
        subcommand = getattr(args, "subparser_name", None) or args.__dict__.get("subparser_name", "")
        if "import" in str(subcommand).lower():
            log.info(
                "\nNext step — denoise with DADA2:\n"
                "  python 03_run_full_metabarcoding_pipeline.py dada2 \\\n"
                "    --project-dir . --marker <MARKER> \\\n"
                "    --trunc-len-f <F> --trunc-len-r <R>"
            )
        elif "dada2" in str(subcommand).lower():
            log.info(
                "\nNext step — check rarefaction depth:\n"
                "  python 04_rarefaction.py \\\n"
                "    --table qiime2/<MARKER>/dada2/table.qza \\\n"
                "    --outdir results/<MARKER>/rarefaction/"
            )
        elif "taxonomy" in str(subcommand).lower():
            log.info(
                "\nNext step — generate taxonomy tables and barplots:\n"
                "  python 07_taxonomy_table.py \\\n"
                "    --taxonomy qiime2/<MARKER>/taxonomy/taxonomy.qza \\\n"
                "    --table    qiime2/<MARKER>/dada2/table.qza \\\n"
                "    --marker   <MARKER> \\\n"
                "    --outdir   results/<MARKER>/taxonomy/"
            )
        elif "diversity" in str(subcommand).lower():
            log.info(
                "\nNext step — taxonomy barplots:\n"
                "  python 07_taxonomy_table.py \\\n"
                "    --taxonomy qiime2/<MARKER>/taxonomy/taxonomy.qza \\\n"
                "    --table    qiime2/<MARKER>/dada2/table.qza \\\n"
                "    --marker   <MARKER> \\\n"
                "    --outdir   results/<MARKER>/taxonomy/"
            )
        return 0
    except Exception as ex:
        eprint(f"\n[ERROR] {ex}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
