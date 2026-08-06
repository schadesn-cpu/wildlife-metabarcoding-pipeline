#!/usr/bin/env python3
"""
config_loader.py
================
Shared configuration loader for the wildlife metabarcoding pipeline.

All pipeline scripts import this module to get their parameters from
pipeline_config.yml rather than hardcoding project-specific values.

Usage (in any pipeline script):
    from config_loader import load_config, get_metadata_path, get_diversity_dir

    cfg = load_config()                              # reads pipeline_config.yml
    root = cfg.root                                  # project root Path
    meta = get_metadata_path(cfg, "16S", "all")      # metadata/qiime/metadata_16S.tsv
    div  = get_diversity_dir(cfg, "16S", "DvT")      # qiime2/16S/.../core_metrics_...
    groups = cfg.groups["primary"]["order"]          # ["Diseased", "Trauma", "Marine"]
    regex  = cfg.samples["id_regex"]                 # "(TV\\d+)"

Config file discovery (in priority order):
    1. Path passed explicitly to load_config(path=...)
    2. $PIPELINE_CONFIG environment variable
    3. pipeline_config.yml in the current working directory
    4. pipeline_config.yml in the directory containing this file
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default config file name
# ---------------------------------------------------------------------------
CONFIG_FILENAME = "pipeline_config.yml"


# ---------------------------------------------------------------------------
# Simple config object — avoids requiring attribute-style access boilerplate
# ---------------------------------------------------------------------------

class PipelineConfig:
    """
    Lightweight wrapper around the raw YAML dict.

    Provides:
      .root          — resolved project root Path
      .name          — project name string
      .email         — contact email string
      .active_markers — list of active marker names
      .markers       — dict of per-marker config blocks
      .metadata      — dict of metadata path blocks per marker
      .groups        — dict of group definitions
      .diversity_dirs — dict of diversity directory paths per marker
      .figures       — dict of figure settings
      .ncbi          — dict of NCBI credentials
      .slurm         — dict of SLURM settings
      .samples       — dict of sample naming config

    All path-valued fields are returned as Path objects resolved
    relative to .root.
    """

    def __init__(self, data: Dict[str, Any], config_path: Path) -> None:
        self._data        = data
        self._config_path = config_path

        # Resolve project root
        raw_root = data.get("project", {}).get("root", ".")
        if Path(raw_root).is_absolute():
            self.root = Path(raw_root).resolve()
        else:
            # Relative root is resolved relative to the config file's directory
            self.root = (config_path.parent / raw_root).resolve()

        project = data.get("project", {})
        self.name           = project.get("name", "unnamed_project")
        self.email          = project.get("email", "")
        self.active_markers = data.get("active_markers", [])
        self.markers        = data.get("markers", {})
        self.metadata       = data.get("metadata", {})
        self.groups         = data.get("groups", {})
        self.diversity_dirs = data.get("diversity_dirs", {})
        self.figures        = data.get("figures", {})
        self.analyses       = data.get("analyses", {})
        self.ncbi           = data.get("ncbi", {})
        self.slurm          = data.get("slurm", {})
        self.samples        = data.get("samples", {})
        self.qc             = data.get("qc", {})

    def resolve(self, rel_path: str) -> Path:
        """Resolve a config-relative path against the project root."""
        p = Path(rel_path)
        if p.is_absolute():
            return p
        return (self.root / p).resolve()

    def __repr__(self) -> str:
        return (
            f"PipelineConfig(name={self.name!r}, "
            f"root={self.root}, "
            f"markers={self.active_markers})"
        )


# ---------------------------------------------------------------------------
# YAML loading — PyYAML is available in the qiime2-amplicon-2024.5 environment
# ---------------------------------------------------------------------------

def _load_yaml(path: Path) -> Dict[str, Any]:
    """
    Load a YAML file using PyYAML.

    PyYAML ships with QIIME2's conda environment and is always available
    when the pipeline is run in the correct conda env. If PyYAML is somehow
    missing, we fall back to a minimal hand-written parser that handles the
    simple key: value structure used by pipeline_config.yml.
    """
    try:
        import yaml
        with path.open(encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    except ImportError:
        log.warning(
            "PyYAML not found — falling back to minimal YAML parser.\n"
            "  Install with: pip install pyyaml --break-system-packages\n"
            "  or: conda install -c conda-forge pyyaml"
        )
        return _minimal_yaml_parser(path)
    except Exception as exc:
        raise ValueError(
            f"Could not parse config file {path}:\n  {exc}\n"
            f"  Check that the file is valid YAML (no tabs, correct indentation)."
        ) from exc


def _minimal_yaml_parser(path: Path) -> Dict[str, Any]:
    """
    Extremely minimal YAML parser for simple key: value configs.

    Handles:
      - top-level and nested key: value pairs (string, int, float, bool)
      - list items (lines starting with '  - ')
      - '#' comments and blank lines

    This fallback is intentionally limited — it handles pipeline_config.yml's
    structure but is not a general YAML parser. PyYAML is strongly preferred.
    """
    result: Dict[str, Any] = {}
    stack: List[tuple] = [(result, -1)]  # (dict, indent_level)

    with path.open(encoding="utf-8") as f:
        for line in f:
            # Strip comment and trailing whitespace
            stripped = line.split("#")[0].rstrip()
            if not stripped.strip():
                continue

            indent = len(stripped) - len(stripped.lstrip())

            # List item
            if stripped.lstrip().startswith("- "):
                val = stripped.lstrip()[2:].strip()
                # Find parent dict and current key
                parent, _ = stack[-1]
                if isinstance(parent, list):
                    parent.append(_coerce(val))
                continue

            # Key: value or Key: (nested block)
            if ":" in stripped:
                key, _, val = stripped.lstrip().partition(":")
                key = key.strip()
                val = val.strip()

                # Pop stack to correct indent level
                while len(stack) > 1 and stack[-1][1] >= indent:
                    stack.pop()

                parent, _ = stack[-1]
                if not isinstance(parent, dict):
                    continue

                if val == "":
                    # Nested block — will be populated by subsequent lines
                    new_dict: Dict[str, Any] = {}
                    parent[key] = new_dict
                    stack.append((new_dict, indent))
                elif val.startswith("["):
                    # Inline list: [a, b, c]
                    items = [x.strip().strip('"\'') for x in
                             val.strip("[]").split(",") if x.strip()]
                    parent[key] = items
                else:
                    parent[key] = _coerce(val)

    return result


def _coerce(val: str) -> Any:
    """Convert a YAML scalar string to an appropriate Python type."""
    if val.lower() in ("true", "yes"):
        return True
    if val.lower() in ("false", "no"):
        return False
    if val.lower() in ("null", "~", ""):
        return None
    # Probe int then float; a ValueError just means "not that type", so fall
    # through to the next. Whatever is left is returned as a trimmed string.
    for converter in (int, float):
        try:
            return converter(val)
        except ValueError:
            continue
    return val.strip('"\'')


# ---------------------------------------------------------------------------
# Config discovery and loading
# ---------------------------------------------------------------------------

def find_config(explicit_path: Optional[str] = None) -> Path:
    """
    Locate the pipeline config file.

    Search order:
      1. explicit_path argument (if provided)
      2. $PIPELINE_CONFIG environment variable
      3. pipeline_config.yml in the current working directory
      4. pipeline_config.yml in the directory containing config_loader.py

    Raises FileNotFoundError with a helpful message if no config is found.
    """
    candidates = []

    if explicit_path:
        candidates.append(Path(explicit_path))

    env_path = os.environ.get("PIPELINE_CONFIG")
    if env_path:
        candidates.append(Path(env_path))

    # Current working directory — first choice when running from project root
    candidates.append(Path.cwd() / CONFIG_FILENAME)

    # Directory containing config_loader.py (scripts/) — unlikely to hold
    # the config, but checked for flat/development layouts
    candidates.append(Path(__file__).resolve().parent / CONFIG_FILENAME)

    # Parent of the directory containing config_loader.py — this is the
    # normal case: config_loader.py is in scripts/, pipeline_config.yml
    # is in the project root one level up.
    candidates.append(Path(__file__).resolve().parent.parent / CONFIG_FILENAME)

    for c in candidates:
        if c.exists():
            log.debug("Config found: %s", c)
            return c.resolve()

    searched = "\n  ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        f"Could not find {CONFIG_FILENAME}.\n"
        f"  Searched:\n  {searched}\n\n"
        f"  To create a starter config:\n"
        f"    python pipeline.py init\n\n"
        f"  Or set the path explicitly:\n"
        f"    python pipeline.py --config /path/to/pipeline_config.yml run"
    )


def load_config(path: Optional[str] = None) -> PipelineConfig:
    """
    Load and return the pipeline configuration.

    Parameters
    ----------
    path : optional explicit path to the config file.
           If None, uses find_config() to locate it automatically.

    Returns
    -------
    PipelineConfig object with all settings resolved relative to project root.

    Raises
    ------
    FileNotFoundError  if no config file can be found
    ValueError         if the config file contains invalid YAML
    """
    config_path = find_config(path)
    log.info("Loading config: %s", config_path)

    data = _load_yaml(config_path)

    cfg = PipelineConfig(data, config_path)
    log.info(
        "Config loaded: project=%r  root=%s  markers=%s",
        cfg.name, cfg.root, cfg.active_markers,
    )
    return cfg


# ---------------------------------------------------------------------------
# Convenience accessors used throughout the pipeline
# ---------------------------------------------------------------------------

def get_metadata_path(cfg: PipelineConfig, marker: str,
                      analysis: str = "all") -> Path:
    """
    Return the resolved metadata file path for a marker + analysis combination.

    Parameters
    ----------
    marker   : e.g. "16S", "MiFish"
    analysis : "all" (full sample set) or "cod" (cause-of-death subset)

    Raises ValueError if the marker or analysis key is not in config.
    """
    try:
        raw = cfg.metadata[marker][analysis]
    except KeyError:
        raise ValueError(
            f"No metadata path for marker={marker!r}, analysis={analysis!r}.\n"
            f"  Check the 'metadata' section of {CONFIG_FILENAME}."
        )
    return cfg.resolve(raw)


def get_diversity_dir(cfg: PipelineConfig, marker: str,
                      analysis: str = "DvT") -> Path:
    """
    Return the core-metrics diversity directory for a marker + analysis.

    Paths are DERIVED from the PathBuilder using the marker's rarefaction depth
    and the analysis name, so the config only needs the depth (a parameter),
    not a hand-written directory tree. If an explicit override happens to be
    present under the optional 'diversity_dirs' config block it is honored, but
    that block is no longer required and is omitted by default.

    Parameters
    ----------
    marker   : e.g. "16S", "MiFish"
    analysis : "DvT", "season", etc. — becomes the {analysis} tag in the path
    """
    override = cfg.diversity_dirs.get(marker, {}).get(analysis)
    if override:
        return cfg.resolve(override)
    return get_paths(cfg).core_metrics_dir(marker, analysis)


def get_classifier_path(cfg: PipelineConfig, marker: str) -> Path:
    """Return the resolved classifier QZA path for a marker."""
    try:
        raw = cfg.markers[marker]["classifier"]
    except KeyError:
        raise ValueError(
            f"No classifier path for marker={marker!r}.\n"
            f"  Check the 'markers' section of {CONFIG_FILENAME}."
        )
    return cfg.resolve(raw)


def get_group_order(cfg: PipelineConfig,
                    group_type: str = "primary") -> List[str]:
    """
    Return the group order list for a given group type.

    Parameters
    ----------
    group_type : "primary", "secondary", or "seasonal"
    """
    try:
        return cfg.groups[group_type]["order"]
    except KeyError:
        return []


def get_dvt_order(cfg: PipelineConfig) -> List[str]:
    """Return the two-group (DvT) comparison order from the primary group config."""
    try:
        return cfg.groups["primary"]["dvt_order"]
    except KeyError:
        # Fall back to first two groups from primary order
        order = get_group_order(cfg, "primary")
        return order[:2] if len(order) >= 2 else order


def validate_config(cfg: PipelineConfig) -> List[str]:
    """
    Run basic validation checks on a loaded config.

    Returns a list of warning/error strings.
    An empty list means the config looks valid.

    Checks:
      - project root exists
      - active_markers are all present in markers block
      - metadata files exist (warns, does not error — may not be generated yet)
      - diversity dirs exist (warns, does not error)
      - classifier files exist (warns if not)
    """
    issues = []

    if not cfg.root.exists():
        issues.append(f"[ERROR] project.root does not exist: {cfg.root}")

    for marker in cfg.active_markers:
        if marker not in cfg.markers:
            issues.append(
                f"[ERROR] Active marker '{marker}' has no config block under 'markers'."
            )

    for marker in cfg.active_markers:
        for analysis in ("all", "cod"):
            try:
                p = get_metadata_path(cfg, marker, analysis)
                if not p.exists():
                    issues.append(
                        f"[WARN]  Metadata file not found (may not be generated yet): {p}"
                    )
            except ValueError as e:
                issues.append(f"[WARN]  {e}")

        try:
            d = get_diversity_dir(cfg, marker, "DvT")
            if not d.exists():
                issues.append(
                    f"[WARN]  Diversity dir not found (run diversity step first): {d}"
                )
        except ValueError as e:
            issues.append(f"[WARN]  {e}")

        try:
            c = get_classifier_path(cfg, marker)
            if not c.exists():
                issues.append(
                    f"[WARN]  Classifier not found (run build_classifiers step first): {c}"
                )
        except ValueError as e:
            issues.append(f"[WARN]  {e}")

    if not cfg.email or "@" not in cfg.email:
        issues.append(
            "[WARN]  project.email is not set — required for NCBI Entrez access."
        )

    return issues


# ---------------------------------------------------------------------------
# PathBuilder — single source of truth for all output paths
#
# Naming convention (enforced here, nowhere else):
#
#   Rarefaction:  r{depth}  (e.g. r8000, r17000, r200)  or  unrarefied
#   Separator:    underscore only — no hyphens in any pipeline-created path
#   Marker:       always first element in filenames
#   File order:   {marker}_{rarefaction}_{analysis}_{type}_{palette}.{ext}
#   QZA/QZV:      always in separate subdirs (core_metrics/ vs group_significance/)
#   Stats TSVs:   in results/  not in qiime2/
#   Taxonomy level: always explicit in filename (_L6_ or _L7_)
#
# Usage:
#   paths = get_paths(cfg)
#
#   # QIIME2 artifact paths
#   paths.dada2_table(marker)                    → qiime2/16S/dada2/table_16S_nocontrols.qza
#   paths.taxonomy_qza(marker)                   → qiime2/16S/taxonomy/taxonomy_16S_silva138_v4.qza
#   paths.core_metrics_dir(marker, "DvT")        → qiime2/16S/diversity/r8000_DvT/core_metrics/
#   paths.group_sig_dir(marker, "DvT")           → qiime2/16S/diversity/r8000_DvT/group_significance/
#
#   # Results paths
#   paths.taxonomy_results_dir(marker)           → results/16S/taxonomy/unrarefied/
#   paths.taxonomy_tsv(marker, "relabund", 6)    → results/16S/taxonomy/unrarefied/16S_taxonomy_relabund_L6_unrarefied.tsv
#   paths.diversity_stats_dir(marker, "DvT")     → results/16S/diversity/r8000_DvT/
#   paths.figure_dir(marker, "DvT")              → results/16S/figures/r8000_DvT/
#   paths.figure_stem(marker, "DvT", "pcoa")     → 16S_r8000_DvT_pcoa_wong
#   paths.manuscript_dir()                       → manuscript_figures/
# ---------------------------------------------------------------------------

class PathBuilder:
    """
    Generates all pipeline output paths from config following the
    project naming convention. Every path in the pipeline should be
    obtained from this class rather than constructed in scripts.

    Get an instance with:  paths = get_paths(cfg)
    """

    def __init__(self, cfg: "PipelineConfig") -> None:
        self._cfg = cfg

    # ── helpers ──────────────────────────────────────────────────────────────

    def _r(self, marker: str) -> str:
        """
        Return the rarefaction tag for a marker: 'r{depth}' or 'unrarefied'.

        Uses the rarefaction_depth from the marker config block. If depth is
        0 or not set, returns 'unrarefied'.
        """
        depth = self._cfg.markers.get(marker, {}).get("rarefaction_depth", 0)
        return f"r{depth}" if depth else "unrarefied"

    def _palette(self) -> str:
        return self._cfg.figures.get("palette", "wong")

    # ── QIIME2 artifact paths (qiime2/) ──────────────────────────────────────

    def qiime2_dir(self, marker: str) -> Path:
        """qiime2/{marker}/"""
        return self._cfg.resolve(f"qiime2/{marker}")

    def imported_dir(self, marker: str) -> Path:
        """qiime2/{marker}/imported/"""
        return self._cfg.resolve(f"qiime2/{marker}/imported")

    def demux_qza(self, marker: str, trimmed: bool = False) -> Path:
        """qiime2/{marker}/imported/demux_{marker}[_trimmed].qza"""
        suffix = "_trimmed" if trimmed else ""
        return self._cfg.resolve(f"qiime2/{marker}/imported/demux_{marker}{suffix}.qza")

    def demux_qzv(self, marker: str, trimmed: bool = False) -> Path:
        """qiime2/{marker}/imported/demux_{marker}[_trimmed].qzv"""
        suffix = "_trimmed" if trimmed else ""
        return self._cfg.resolve(f"qiime2/{marker}/imported/demux_{marker}{suffix}.qzv")

    def manifest_tsv(self, marker: str) -> Path:
        """qiime2/{marker}/imported/manifest_{marker}.tsv"""
        return self._cfg.resolve(f"qiime2/{marker}/imported/manifest_{marker}.tsv")

    def dada2_dir(self, marker: str) -> Path:
        """qiime2/{marker}/dada2/"""
        return self._cfg.resolve(f"qiime2/{marker}/dada2")

    def dada2_table(self, marker: str, nocontrols: bool = True) -> Path:
        """qiime2/{marker}/dada2/table_{marker}[_nocontrols].qza"""
        suffix = "_nocontrols" if nocontrols else ""
        return self._cfg.resolve(f"qiime2/{marker}/dada2/table_{marker}{suffix}.qza")

    def rep_seqs_qza(self, marker: str) -> Path:
        """qiime2/{marker}/dada2/rep_seqs_{marker}.qza"""
        return self._cfg.resolve(f"qiime2/{marker}/dada2/rep_seqs_{marker}.qza")

    def denoising_stats_qza(self, marker: str) -> Path:
        """qiime2/{marker}/dada2/denoising_stats_{marker}.qza"""
        return self._cfg.resolve(f"qiime2/{marker}/dada2/denoising_stats_{marker}.qza")

    def denoising_stats_qzv(self, marker: str) -> Path:
        """qiime2/{marker}/dada2/denoising_stats_{marker}.qzv"""
        return self._cfg.resolve(f"qiime2/{marker}/dada2/denoising_stats_{marker}.qzv")

    def taxonomy_dir(self, marker: str) -> Path:
        """qiime2/{marker}/taxonomy/"""
        return self._cfg.resolve(f"qiime2/{marker}/taxonomy")

    def taxonomy_qza(self, marker: str) -> Path:
        """
        qiime2/{marker}/taxonomy/taxonomy_{marker}_{db}.qza

        The database tag is derived from the classifier filename stem, e.g.
        'silva-138-99-nb-classifier-515-806' → 'silva138_v4'.
        Falls back to marker name only if classifier not configured.
        """
        classifier = self._cfg.markers.get(marker, {}).get("classifier", "")
        db_tag = _classifier_to_db_tag(classifier) if classifier else marker
        return self._cfg.resolve(
            f"qiime2/{marker}/taxonomy/taxonomy_{marker}_{db_tag}.qza"
        )

    def taxonomy_qzv(self, marker: str) -> Path:
        """qiime2/{marker}/taxonomy/taxonomy_{marker}_{db}.qzv"""
        classifier = self._cfg.markers.get(marker, {}).get("classifier", "")
        db_tag = _classifier_to_db_tag(classifier) if classifier else marker
        return self._cfg.resolve(
            f"qiime2/{marker}/taxonomy/taxonomy_{marker}_{db_tag}.qzv"
        )

    def rooted_tree_qza(self, marker: str) -> Path:
        """qiime2/{marker}/tree/rooted_tree_{marker}.qza"""
        return self._cfg.resolve(f"qiime2/{marker}/tree/rooted_tree_{marker}.qza")

    def diversity_dir(self, marker: str) -> Path:
        """qiime2/{marker}/diversity/"""
        return self._cfg.resolve(f"qiime2/{marker}/diversity")

    # ── run_full engine I/O contract ─────────────────────────────────────────
    # The engine (05_run_full_metabarcoding_pipeline.py) writes artifacts scoped
    # by marker AND dataset: qiime2/{marker}/{dataset}/{stage}/{generic-name}
    # (dataset defaults to 'all'; subsets like 'DvT' for cause-of-death). These
    # methods describe that on-disk contract exactly — generic filenames, dataset
    # in the path — which is distinct from the marker-suffixed results methods
    # elsewhere in this class. The orchestrator reads engine outputs through
    # these so the layout lives in one place, not in f-strings across pipeline.py.

    def qiime2_root(self) -> Path:
        """qiime2/  (the engine's --outdir)"""
        return self._cfg.resolve("qiime2")

    def reads_dir(self) -> Path:
        """reads/"""
        return self._cfg.resolve("reads")

    def primers_detected_tsv(self) -> Path:
        """reports/primers_detected.tsv  (written by the primer advisor)"""
        return self._cfg.resolve("reports/primers_detected.tsv")

    def engine_scope_dir(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/"""
        return self._cfg.resolve(f"qiime2/{marker}/{dataset}")

    def engine_imported_dir(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/imported/"""
        return self.engine_scope_dir(marker, dataset) / "imported"

    def engine_demux_qzv(self, marker: str, dataset: str = "all",
                         trimmed: bool = False) -> Path:
        """qiime2/{marker}/{dataset}/imported/demux[_trimmed].qzv"""
        name = "demux_trimmed.qzv" if trimmed else "demux.qzv"
        return self.engine_imported_dir(marker, dataset) / name

    def engine_dada2_params_json(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/imported/dada2_params.json (from dada2_advisor)"""
        return self.engine_imported_dir(marker, dataset) / "dada2_params.json"

    def engine_dada2_dir(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/dada2/"""
        return self.engine_scope_dir(marker, dataset) / "dada2"

    def engine_denoising_stats_qza(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/dada2/denoising-stats.qza"""
        return self.engine_dada2_dir(marker, dataset) / "denoising-stats.qza"

    def engine_table_qza(self, marker: str, dataset: str = "all",
                         nocontrols: bool = False) -> Path:
        """qiime2/{marker}/{dataset}/dada2/table[_no_controls].qza"""
        name = "table_no_controls.qza" if nocontrols else "table.qza"
        return self.engine_dada2_dir(marker, dataset) / name

    def engine_rep_seqs_qza(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/dada2/rep-seqs.qza"""
        return self.engine_dada2_dir(marker, dataset) / "rep-seqs.qza"

    def engine_taxonomy_exported_tsv(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/taxonomy-exported/taxonomy.tsv (qiime tools export)"""
        return self.engine_scope_dir(marker, dataset) / "taxonomy-exported" / "taxonomy.tsv"

    def engine_rep_seqs_exported_fasta(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/rep-seqs-exported/dna-sequences.fasta"""
        return self.engine_scope_dir(marker, dataset) / "rep-seqs-exported" / "dna-sequences.fasta"

    def engine_taxonomy_dir(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/taxonomy/"""
        return self.engine_scope_dir(marker, dataset) / "taxonomy"

    def engine_taxonomy_qza(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/taxonomy/taxonomy.qza"""
        return self.engine_taxonomy_dir(marker, dataset) / "taxonomy.qza"

    def engine_diversity_dir(self, marker: str, dataset: str = "all") -> Path:
        """qiime2/{marker}/{dataset}/diversity/"""
        return self.engine_scope_dir(marker, dataset) / "diversity"

    def engine_results_dir(self, marker: str, dataset: str = "all") -> Path:
        """results/{marker}/{dataset}/"""
        return self._cfg.resolve(f"results/{marker}/{dataset}")

    def engine_taxonomy_results_dir(self, marker: str, dataset: str = "all") -> Path:
        """results/{marker}/{dataset}/taxonomy/"""
        return self.engine_results_dir(marker, dataset) / "taxonomy"

    def engine_diversity_results_dir(self, marker: str, dataset: str = "all") -> Path:
        """results/{marker}/{dataset}/diversity/"""
        return self.engine_results_dir(marker, dataset) / "diversity"

    def engine_rarefaction_results_dir(self, marker: str, dataset: str = "all") -> Path:
        """results/{marker}/{dataset}/rarefaction/"""
        return self.engine_results_dir(marker, dataset) / "rarefaction"

    def engine_presence_absence_results_dir(self, marker: str, dataset: str = "all") -> Path:
        """results/{marker}/{dataset}/presence_absence/"""
        return self.engine_results_dir(marker, dataset) / "presence_absence"

    def engine_blast_results_dir(self, marker: str, dataset: str = "all") -> Path:
        """results/{marker}/{dataset}/blast/"""
        return self.engine_results_dir(marker, dataset) / "blast"

    def engine_manifest_tsv(self, marker: str) -> Path:
        """qiime2/{marker}/imported/manifest_{marker}.tsv (no dataset; matches make_manifests)"""
        return self._cfg.resolve(f"qiime2/{marker}/imported/manifest_{marker}.tsv")

    def rarefied_table_qza(self, marker: str, analysis: str) -> Path:
        """qiime2/{marker}/diversity/{r}_{ analysis}/core_metrics/rarefied_table.qza"""
        r = self._r(marker)
        return self._cfg.resolve(
            f"qiime2/{marker}/diversity/{r}_{analysis}/core_metrics/rarefied_table.qza"
        )

    def core_metrics_dir(self, marker: str, analysis: str) -> Path:
        """
        qiime2/{marker}/diversity/r{depth}_{analysis}/core_metrics/

        Contains QZA files: distance matrices, PCoA results, alpha vectors.
        Maps to what QIIME2 calls 'core-metrics[-phylogenetic]' output.
        """
        r = self._r(marker)
        return self._cfg.resolve(
            f"qiime2/{marker}/diversity/{r}_{analysis}/core_metrics"
        )

    def group_sig_dir(self, marker: str, analysis: str) -> Path:
        """
        qiime2/{marker}/diversity/r{depth}_{analysis}/group_significance/

        Contains QZV files: PERMANOVA, PERMDISP, alpha-group-significance.
        Separated from core_metrics/ so artifacts and visualizations are distinct.
        """
        r = self._r(marker)
        return self._cfg.resolve(
            f"qiime2/{marker}/diversity/{r}_{analysis}/group_significance"
        )

    def permanova_qzv(self, marker: str, analysis: str, metric: str) -> Path:
        """qiime2/{marker}/diversity/r{depth}_{analysis}/group_significance/permanova_{group}_{metric}.qzv"""
        r          = self._r(marker)
        group_col  = self._cfg.groups.get("primary", {}).get("column", "Group")
        return self._cfg.resolve(
            f"qiime2/{marker}/diversity/{r}_{analysis}/"
            f"group_significance/permanova_{group_col}_{metric}.qzv"
        )

    # ── Results paths (results/) ──────────────────────────────────────────────

    def results_dir(self, marker: str) -> Path:
        """results/{marker}/"""
        return self._cfg.resolve(f"results/{marker}")

    def qc_dir(self) -> Path:
        """results/qc/"""
        return self._cfg.resolve("results/qc")

    def taxonomy_results_dir(self, marker: str) -> Path:
        """results/{marker}/taxonomy/unrarefied/"""
        return self._cfg.resolve(f"results/{marker}/taxonomy/unrarefied")

    def taxonomy_tsv(self, marker: str,
                     table_type: str = "relabund",
                     level: int = 6) -> Path:
        """
        results/{marker}/taxonomy/unrarefied/{marker}_taxonomy_{type}_L{level}_unrarefied.tsv

        table_type: 'relabund' | 'counts'
        level:      taxonomy level (6 = family, 7 = genus)
        """
        return self._cfg.resolve(
            f"results/{marker}/taxonomy/unrarefied/"
            f"{marker}_taxonomy_{table_type}_L{level}_unrarefied.tsv"
        )

    def diversity_stats_dir(self, marker: str, analysis: str) -> Path:
        """results/{marker}/diversity/r{depth}_{analysis}/"""
        r = self._r(marker)
        return self._cfg.resolve(f"results/{marker}/diversity/{r}_{analysis}")

    def beta_stats_tsv(self, marker: str, analysis: str) -> Path:
        """results/{marker}/diversity/r{depth}_{analysis}/beta_stats_{group}.tsv"""
        r         = self._r(marker)
        group_col = self._cfg.groups.get("primary", {}).get("column", "Group")
        return self._cfg.resolve(
            f"results/{marker}/diversity/{r}_{analysis}/beta_stats_{group_col}.tsv"
        )

    def alpha_stats_tsv(self, marker: str, analysis: str) -> Path:
        """results/{marker}/diversity/r{depth}_{analysis}/alpha_stats_{group}.tsv"""
        r         = self._r(marker)
        group_col = self._cfg.groups.get("primary", {}).get("column", "Group")
        return self._cfg.resolve(
            f"results/{marker}/diversity/{r}_{analysis}/alpha_stats_{group_col}.tsv"
        )

    def figure_dir(self, marker: str, analysis: str,
                   annotated: bool = False) -> Path:
        """
        Annotated:  results/figures_annotated/{marker}/r{depth}_{analysis}/
        Manuscript: results/figures_manuscript/{marker}/r{depth}_{analysis}/

        Taxonomy (unrarefied) figures use 'unrarefied' instead of r{depth}_{analysis}.
        """
        r        = self._r(marker)
        base_key = "annotated_dir" if annotated else "manuscript_dir"
        base     = self._cfg.figures.get(
            base_key,
            "results/figures_annotated" if annotated else "results/figures_manuscript",
        )
        return self._cfg.resolve(f"{base}/{marker}/{r}_{analysis}")

    def taxonomy_figure_dir(self, marker: str,
                            annotated: bool = False) -> Path:
        """
        results/figures_{annotated|manuscript}/{marker}/unrarefied/

        Taxonomy barplots are always unrarefied so get their own subdir.
        """
        base_key = "annotated_dir" if annotated else "manuscript_dir"
        base     = self._cfg.figures.get(
            base_key,
            "results/figures_annotated" if annotated else "results/figures_manuscript",
        )
        return self._cfg.resolve(f"{base}/{marker}/unrarefied")

    def figure_stem(self, marker: str, analysis: str,
                    figure_type: str, annotated: bool = False) -> str:
        """
        Return the output filename stem (no extension, no directory) for a figure.

        Pattern:  {marker}_{rarefaction}_{analysis}_{type}_{palette}[_annotated]
        Example:  16S_r8000_DvT_pcoa_wong
                  16S_r8000_DvT_pcoa_wong_annotated
        """
        r       = self._r(marker)
        palette = self._palette()
        stem    = f"{marker}_{r}_{analysis}_{figure_type}_{palette}"
        return f"{stem}_annotated" if annotated else stem

    def taxonomy_figure_stem(self, marker: str, analysis: str,
                             level: int = 6,
                             annotated: bool = False) -> str:
        """
        {marker}_unrarefied_{analysis}_taxonomy_L{level}_{palette}[_annotated]
        Example:  16S_unrarefied_DvT_taxonomy_L6_wong
        """
        palette = self._palette()
        stem    = f"{marker}_unrarefied_{analysis}_taxonomy_L{level}_{palette}"
        return f"{stem}_annotated" if annotated else stem

    def manuscript_figures_dir(self) -> Path:
        """manuscript_figures/  — flat directory with Fig01_… naming"""
        return self._cfg.resolve("manuscript_figures")


def _classifier_to_db_tag(classifier_path: str) -> str:
    """
    Derive a short human-readable database tag from a classifier filename.

    Examples:
      silva-138-99-nb-classifier-515-806.qza  → silva138_v4
      mitofish-12S-mifish-gavia-classifier.qza → mitofish
      ncbi-cytb-vertebrata-classifier.qza      → ncbi_cytb
      silva-138-V9-classifier.qza              → silva138_v9
      unite-ver10-99-nb-classifier.qza         → unite_v10
    """
    stem = Path(classifier_path).stem.lower()
    # Remove common noise words
    for word in ("nb", "classifier", "99", "97"):
        stem = stem.replace(f"-{word}-", "-").replace(f"-{word}", "")

    if "silva" in stem and "138" in stem:
        if "515" in stem or "v4" in stem or "806" in stem:
            return "silva138_v4"
        elif "v9" in stem or "1391" in stem:
            return "silva138_v9"
        return "silva138"
    if "mitofish" in stem or "mifish" in stem:
        return "mitofish"
    if "ncbi" in stem and "cytb" in stem:
        return "ncbi_cytb"
    if "unite" in stem:
        ver = next((p for p in stem.split("-") if p.startswith("ver")), "")
        num = ver.replace("ver", "") if ver else ""
        return f"unite_v{num}" if num else "unite"
    if "pr2" in stem:
        return "pr2"

    # Generic fallback: take first two non-trivial hyphen-segments
    parts = [p for p in stem.replace("_", "-").split("-")
             if len(p) > 2 and not p.isdigit()]
    return "_".join(parts[:2]) if parts else Path(classifier_path).stem


def get_paths(cfg: "PipelineConfig") -> PathBuilder:
    """
    Return a PathBuilder instance for the given config.

    This is the recommended way to get output paths in pipeline scripts:

        from config_loader import load_config, get_paths
        cfg   = load_config()
        paths = get_paths(cfg)

        core_dir = paths.core_metrics_dir("16S", "DvT")
        fig_stem = paths.figure_stem("16S", "DvT", "pcoa")
    """
    return PathBuilder(cfg)


# ---------------------------------------------------------------------------
# Standalone: print config summary
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.WARNING)

    path_arg = sys.argv[1] if len(sys.argv) > 1 else None
    try:
        cfg = load_config(path_arg)
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        sys.exit(1)

    print(f"\nConfig loaded successfully")
    print(f"  Project:  {cfg.name}")
    print(f"  Root:     {cfg.root}")
    print(f"  Email:    {cfg.email or '(not set)'}")
    print(f"  Markers:  {', '.join(cfg.active_markers)}")
    print(f"  Palette:  {cfg.figures.get('palette', 'wong')}")
    print(f"\nRunning validation checks...")
    issues = validate_config(cfg)
    if not issues:
        print("  All checks passed.")
    else:
        for iss in issues:
            print(f"  {iss}")
    print()
