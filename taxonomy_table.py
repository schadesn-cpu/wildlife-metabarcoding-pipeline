#!/usr/bin/env python3
"""
taxonomy_table.py
=================
Step: taxonomy export (part of the taxonomy stage)

Purpose:
    Export QIIME 2 taxonomy classifications to human-readable TSV tables for any
    marker, using the *correct* order of operations for composition barplots:
    start from the UNRAREFIED table, apply marker-aware filters (remove host,
    off-target kingdoms), collapse to genus/species, then compute relative
    abundance among the remaining classified reads.

    Why unrarefied: rarefying a host-dominated table (e.g. gut 16S) before
    filtering leaves almost only Unassigned reads, producing near-empty barplots.
    Diversity metrics still use the rarefied table; composition does not.

Inputs:
    --taxonomy  FeatureData[Taxonomy] .qza   (or derived from --marker via config)
    --table     FeatureTable[Frequency] .qza (unrarefied; or derived from --marker)
    --marker    marker gene (controls depth + auto-filter defaults)
    pipeline_config.yml  active_markers, output locations

Outputs (in --outdir, default results/<marker>/all/taxonomy/):
    taxonomy_counts_L{N}_{marker}.tsv      collapsed count table
    taxonomy_relabund_L{N}_{marker}.tsv    relative abundance (barplot input)
    taxonomy_top{N}_L{N}_{marker}.tsv      top-N summary
    logs/run_manifest.jsonl                run appended on completion

Marker filter defaults (auto-applied unless --include/--exclude/--no-auto-filter):
    16S    bacteria only      exclude mitochondria, chloroplast, Eukaryota, Archaea
    18S    all eukaryotes     exclude Bacteria, Archaea
    ITS    Fungi              exclude Bacteria, Archaea
    MiFish Actinopteri        exclude Bacteria, Archaea
    cytb   Vertebrata         exclude Bacteria, Viruses, Archaea
    COI    Metazoa            exclude Bacteria, Viruses, Archaea, Viridiplantae, Fungi

Usage:
    python taxonomy_table.py --marker 16S
    python taxonomy_table.py --marker 16S --taxonomy t.qza --table tbl.qza --outdir out/

Requirements:
    Python >= 3.8, pandas, numpy. No QIIME 2 needed (reads .qza zips directly).
"""

from __future__ import annotations

import argparse
import logging
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# --- make config_loader and the utils package importable regardless of cwd ---
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from config_loader import load_config, get_paths            # noqa: E402
from utils import checkpoint, provenance                    # noqa: E402

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Marker configuration
# ---------------------------------------------------------------------------

MARKER_CONFIG: Dict[str, dict] = {
    "16S": {
        "description": "SILVA 138 (16S rRNA V4)",
        "levels": 7,
        "default_level": 6,  # genus
        "level_names": ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
        "separator": ";",
        "prefix_strip": True,   # strip d__, p__, etc.
    },
    "18S": {
        "description": "PR2 v5.0.0 (18S rRNA)",
        "levels": 8,
        "default_level": 7,  # genus
        "level_names": ["Domain", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"],
        "separator": ";",
        "prefix_strip": True,
    },
    "ITS": {
        "description": "UNITE v10 (ITS)",
        "levels": 7,
        "default_level": 6,  # genus
        "level_names": ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
        "separator": ";",
        "prefix_strip": True,
    },
    "MiFish": {
        "description": "MARES v2 / MitoFish (12S MiFish)",
        "levels": 7,
        "default_level": 7,  # species — fish ID is the goal
        "level_names": ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
        "separator": ";",
        "prefix_strip": True,   # MitoFish QIIME DB uses d__, p__, c__, etc.
    },
    "cytb": {
        "description": "NCBI vertebrate cytochrome b",
        "levels": 7,
        "default_level": 7,  # species
        "level_names": ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
        "separator": ";",
        "prefix_strip": False,
    },
    "COI": {
        "description": "MIDORI2 (COI, cytochrome c oxidase I)",
        "levels": 7,
        "default_level": 7,  # species
        "level_names": ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"],
        "separator": ";",
        "prefix_strip": True,
    },
}

# Taxonomy filters applied automatically per marker when --include/--exclude
# are not set on the command line.  Mirrors MARKER_FILTER_DEFAULTS in
# 03_run_full_metabarcoding_pipeline.py — keep these in sync.
MARKER_FILTER_DEFAULTS: Dict[str, Dict[str, Optional[str]]] = {
    "16S": {
        "include": None,
        "exclude": "mitochondria,chloroplast,Eukaryota,Archaea",
    },
    "18S": {
        "include": None,
        "exclude": "Bacteria,Archaea",
    },
    "ITS": {
        "include": "Fungi",
        "exclude": "Bacteria,Archaea",
    },
    "MiFish": {
        # MitoFish QIIME DB uses QIIME-style prefixes (d__/p__/c__) and
        # classifies fish as c__Actinopteri, not Vertebrata.
        # Gavia (loon host) falls back to c__Actinopteri with no species hit,
        # so it is excluded post-collapse via the biogeographic filter in 09_.
        "include": "Actinopteri",
        "exclude": "Bacteria,Archaea",
    },
    "cytb": {
        "include": "Vertebrata",
        "exclude": "Bacteria,Viruses,Archaea",
    },
    "COI": {
        "include": "Metazoa",
        "exclude": "Bacteria,Viruses,Archaea,Viridiplantae,Fungi",
    },
}


# ---------------------------------------------------------------------------
# QZA loading helpers
# ---------------------------------------------------------------------------

def _normalize_taxonomy_columns(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    """Standardize columns to Feature ID / Taxon / Confidence."""
    df.columns = df.columns.str.strip()
    col_map: Dict[str, str] = {}
    for c in df.columns:
        cl = c.lower()
        if "feature" in cl or cl == "id":
            col_map[c] = "Feature ID"
        elif any(x in cl for x in ("taxon", "taxonomy", "classification")):
            col_map[c] = "Taxon"
        elif "confidence" in cl:
            col_map[c] = "Confidence"
    df = df.rename(columns=col_map)
    if "Feature ID" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Feature ID"})
    if "Taxon" not in df.columns and len(df.columns) > 1:
        df = df.rename(columns={df.columns[1]: "Taxon"})
    log.info("  Taxonomy: %d features from %s", len(df), source_name)
    return df


def load_taxonomy_qza(qza_path: Path) -> pd.DataFrame:
    """
    Read a FeatureData[Taxonomy] QZA directly (no qiime CLI needed).
    Returns DataFrame with columns: Feature ID, Taxon, Confidence.
    """
    with zipfile.ZipFile(qza_path, "r") as zf:
        tax_name = next(
            (n for n in zf.namelist() if n.endswith("taxonomy.tsv")), None
        )
        if tax_name is None:
            log.error("No taxonomy.tsv found inside %s", qza_path)
            sys.exit(1)
        with zf.open(tax_name) as fh:
            df = pd.read_csv(fh, sep="\t", dtype=str)
    return _normalize_taxonomy_columns(df, qza_path.name)


def load_taxonomy_tsv(tsv_path: Path) -> pd.DataFrame:
    """
    Read an exported / BLAST-refined taxonomy TSV (Feature ID, Taxon, Confidence).
    Returns DataFrame with the same normalized columns as load_taxonomy_qza.
    """
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    return _normalize_taxonomy_columns(df, tsv_path.name)


def load_taxonomy_source(path: Path) -> pd.DataFrame:
    """Dispatch to the QZA or TSV loader based on file suffix."""
    if str(path).endswith(".qza"):
        return load_taxonomy_qza(path)
    return load_taxonomy_tsv(path)


def load_feature_table_qza(qza_path: Path) -> pd.DataFrame:
    """
    Export a FeatureTable[Frequency] QZA to a temp dir via qiime tools export
    and return it as a DataFrame (features × samples).
    """
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        biom_dir = tmpdir / "biom_export"
        biom_dir.mkdir()

        result = subprocess.run(
            ["qiime", "tools", "export",
             "--input-path", str(qza_path),
             "--output-path", str(biom_dir)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            log.error("qiime tools export failed:\n%s", result.stderr)
            sys.exit(1)

        biom_files = list(biom_dir.glob("*.biom"))
        if not biom_files:
            log.error("No .biom found after exporting %s", qza_path)
            sys.exit(1)

        tsv_path = tmpdir / "feature_table.tsv"
        result = subprocess.run(
            ["biom", "convert", "-i", str(biom_files[0]),
             "-o", str(tsv_path), "--to-tsv"],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            log.error("biom convert failed:\n%s", result.stderr)
            log.error(
                "Install biom-format:  conda install -c conda-forge biom-format"
            )
            sys.exit(1)

        df = pd.read_csv(tsv_path, sep="\t", skiprows=1, index_col=0, dtype={0: str})
        df.index.name = "Feature ID"
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
        log.info(
            "  Feature table: %d features × %d samples from %s",
            len(df), len(df.columns), qza_path.name,
        )
        return df


# ---------------------------------------------------------------------------
# Taxonomy filtering
# ---------------------------------------------------------------------------

def apply_taxonomy_filter(
    table: pd.DataFrame,
    tax_df: pd.DataFrame,
    include: List[str],
    exclude: List[str],
    min_freq: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter feature table and taxonomy by include/exclude taxon strings and
    minimum total read count.

    All matching is case-insensitive against the full Taxon string.
    """
    tax_idx = tax_df.set_index("Feature ID")["Taxon"]
    fids = table.index.tolist()
    n_start = len(fids)

    # Frequency filter
    if min_freq > 1:
        totals = table.sum(axis=1)
        fids = totals[totals >= min_freq].index.tolist()
        log.info("  min-freq filter (≥%d reads): %d → %d features",
                 min_freq, n_start, len(fids))

    # Include filter
    if include:
        fids = [f for f in fids
                if any(s.lower() in str(tax_idx.get(f, "")).lower() for s in include)]
        log.info("  include filter (%s): %d features retained", include, len(fids))

    # Exclude filter
    if exclude:
        fids = [f for f in fids
                if not any(s.lower() in str(tax_idx.get(f, "")).lower() for s in exclude)]
        log.info("  exclude filter (%s): %d features retained", exclude, len(fids))

    if not fids:
        log.error("No features remain after filtering. Check --include/--exclude.")
        sys.exit(1)

    table = table.loc[fids]
    tax_df = tax_df[tax_df["Feature ID"].isin(fids)].copy()
    return table, tax_df


# ---------------------------------------------------------------------------
# Genus label parsing
# ---------------------------------------------------------------------------

def _clean_taxon_part(part: str, do_strip: bool) -> str:
    """Strip QIIME2 rank prefixes (e.g. 'g__') and surrounding whitespace/underscores."""
    if do_strip:
        part = re.sub(r"^[a-z]__", "", part)
    return part.strip().strip("_")


def _truncate_taxon(taxon: str, sep: str, do_strip: bool, level: int) -> str:
    """Truncate a full taxonomy string to the requested level, padding with 'Unclassified'."""
    if pd.isna(taxon) or not str(taxon).strip():
        return "Unclassified"
    taxon_str = str(taxon).strip()
    if do_strip:
        taxon_str = re.sub(r"[a-z]__", "", taxon_str)
    parts = [part.strip() for part in taxon_str.split(sep)][:level]
    while len(parts) < level:
        parts.append("Unclassified")
    return sep.join(parts)


def clean_taxon_label(taxon: str, level: int, cfg: dict) -> Optional[str]:
    """
    Parse a full taxonomy string and return a clean display label at `level`.

    For genus level (16S/18S SILVA-style), produces labels like:
      - "Cetobacterium"
      - "uncl. Fusobacteriaceae"
      - "uncl. Clostridiales"
    Returns None for rows that should be dropped (Unassigned, pure host reads).
    """
    sep = cfg["separator"]
    do_strip = cfg.get("prefix_strip", False)

    if pd.isna(taxon) or str(taxon).strip().lower() in ("", "unassigned", "no blast hit"):
        return None

    taxon_str = str(taxon).strip()
    parts = [part.strip() for part in taxon_str.split(sep)]
    cleaned = [_clean_taxon_part(part, do_strip) for part in parts]

    # At the requested level, try to find the most specific non-empty name
    target = cleaned[level - 1] if len(cleaned) >= level else ""

    if target and target not in ("", "__"):
        # Fix common SILVA underscore artifacts
        target = (target
                  .replace("_sensu_stricto_1", " sensu stricto 1")
                  .replace("_sensu_stricto_4", " sensu stricto 4")
                  .replace("_", " "))
        # Treat "uncultured" genus as unclassified — go up to family
        if target.lower() == "uncultured":
            family_name = cleaned[level - 2] if len(cleaned) >= level - 1 else ""
            if family_name and family_name not in ("", "__"):
                family_name = family_name.replace("_", " ")
                return f"uncl. {family_name}"
            return None
        return target

    # Walk back up to find the deepest named level (including Domain/Kingdom)
    for i in range(level - 2, -1, -1):
        name = cleaned[i] if len(cleaned) > i else ""
        if name and name not in ("", "__"):
            name = name.replace("_", " ")
            return f"uncl. {name}"

    return None


# ---------------------------------------------------------------------------
# Table collapse and relative abundance
# ---------------------------------------------------------------------------

def collapse_to_level(
    table: pd.DataFrame,
    tax_df: pd.DataFrame,
    marker: str,
    level: int,
    use_clean_labels: bool = True,
) -> pd.DataFrame:
    """
    Collapse feature table to the given taxonomic level.

    If use_clean_labels=True, rows get short display names (e.g. "Cetobacterium").
    If False, rows get the full truncated taxonomy string (for programmatic use).
    Returns a DataFrame: taxa × samples.
    """
    cfg = MARKER_CONFIG[marker]
    level_names = cfg["level_names"]

    if level < 1 or level > len(level_names):
        log.warning("Level %d out of range for %s — using default %d",
                    level, marker, cfg["default_level"])
        level = cfg["default_level"]

    tax_idx = tax_df.set_index("Feature ID")["Taxon"]

    if use_clean_labels:
        labels = {
            fid: clean_taxon_label(tax_idx.get(fid, ""), level, cfg)
            for fid in table.index
        }
        # Drop features where label is None (unclassified at all levels)
        keep = [f for f, lbl in labels.items() if lbl is not None]
        if len(keep) < len(table):
            log.info(
                "  %d features dropped (unclassifiable at level %d); %d retained",
                len(table) - len(keep), level, len(keep),
            )
        table = table.loc[keep]
        label_series = pd.Series({f: labels[f] for f in keep})
    else:
        # Full taxonomy string truncated to level
        sep = cfg["separator"]
        do_strip = cfg.get("prefix_strip", False)

        label_series = tax_idx.reindex(table.index).apply(
            lambda taxon: _truncate_taxon(taxon, sep, do_strip, level)
        )

    collapsed = table.copy()
    collapsed.index = label_series.values
    collapsed = collapsed.groupby(level=0).sum()
    collapsed.index.name = level_names[level - 1]
    return collapsed


def to_relative_abundance(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Convert counts to relative abundance (column-wise, 0–1).
    Samples with zero total reads are left as zero rather than NaN.
    """
    col_sums = counts.sum(axis=0)
    col_sums[col_sums == 0] = 1
    return counts.div(col_sums, axis=1)


def drop_empty_samples(
    counts: pd.DataFrame,
    min_reads: int = 1,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Remove samples with fewer than min_reads classified reads.
    Returns (filtered_table, list_of_dropped_sample_ids).
    """
    totals = counts.sum(axis=0)
    keep = totals[totals >= min_reads].index.tolist()
    dropped = totals[totals < min_reads].index.tolist()
    if dropped:
        log.warning(
            "  Dropping %d sample(s) with <%d classified reads: %s",
            len(dropped), min_reads, dropped,
        )
    return counts[keep], dropped


def top_n_table(relabund: pd.DataFrame, n: int) -> pd.DataFrame:
    """Top-N taxa by mean relative abundance, with mean_relabund as first column."""
    means = relabund.mean(axis=1).sort_values(ascending=False)
    top = means.head(n)
    result = relabund.loc[top.index].copy()
    result.insert(0, "mean_relabund", top.values)
    return result


# ---------------------------------------------------------------------------
# Per-ASV summary
# ---------------------------------------------------------------------------

def _parse_taxon_row(taxon, sep: str, do_strip: bool, level_names: List[str]) -> List[str]:
    """
    Split a taxonomy string into a fixed-length list aligned to level_names.

    Strips rank prefixes if do_strip is set, splits on the configured separator,
    pads short strings with 'Unclassified', and replaces any empty parts with
    'Unclassified'.
    """
    if pd.isna(taxon) or str(taxon).strip().lower() in ("", "unassigned"):
        return ["Unclassified"] * len(level_names)
    taxon_str = str(taxon).strip()
    if do_strip:
        taxon_str = re.sub(r"[a-z]__", "", taxon_str)
    parts = [part.strip() for part in taxon_str.split(sep)][:len(level_names)]
    while len(parts) < len(level_names):
        parts.append("Unclassified")
    return [part or "Unclassified" for part in parts]


def build_asv_summary(
    table: pd.DataFrame,
    tax_df: pd.DataFrame,
    marker: str,
) -> pd.DataFrame:
    """
    One row per ASV: Feature ID, Taxon, Confidence, parsed level columns,
    total_reads, samples_present.
    """
    cfg = MARKER_CONFIG[marker]
    sep = cfg["separator"]
    do_strip = cfg.get("prefix_strip", False)
    level_names = cfg["level_names"]

    totals   = table.sum(axis=1).rename("total_reads")
    presence = (table > 0).sum(axis=1).rename("samples_present")

    level_df = pd.DataFrame(
        [_parse_taxon_row(tax_df.set_index("Feature ID")["Taxon"].get(feature_id, ""),
                          sep, do_strip, level_names)
         for feature_id in table.index],
        index=table.index,
        columns=level_names,
    )
    summary = tax_df.set_index("Feature ID").loc[table.index].copy()
    summary = pd.concat([summary, level_df, totals, presence], axis=1)
    summary = summary.sort_values("total_reads", ascending=False)
    return summary


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def save_tsv(df: pd.DataFrame, path: Path, dry_run: bool = False) -> None:
    """
Write a DataFrame to a tab-separated file, creating parent dirs as needed.

    Logs the output path and dimensions. In dry-run mode the file is not
    written but the log message is still printed so the intended output can
    be audited before committing to a full run.
    """
    log.info("  Writing: %s  (%d rows × %d cols)", path, len(df), len(df.columns))
    if dry_run:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=True)


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def build_taxonomy_tables(
    taxonomy_qza: Path,
    table_qza: Path,
    marker: str,
    outdir: Path,
    level: Optional[int],
    top_n: int,
    include: List[str],
    exclude: List[str],
    min_freq: int,
    min_sample_reads: int,
    dry_run: bool,
) -> None:
    """
Build and write taxonomy summary TSVs from QIIME2 taxonomy and feature table QZAs.

    Reads both QZAs directly via zipfile (no QIIME2 installation required),
    applies marker-aware include/exclude taxonomy filters, collapses features
    to the target taxonomic level, drops samples below min_sample_reads, and
    writes four TSVs to outdir:

      taxonomy_summary_{marker}.tsv      — per-feature taxonomy with read counts
      taxonomy_counts_L{N}_{marker}.tsv  — collapsed count table (features x samples)
      taxonomy_relabund_L{N}_{marker}.tsv — relative abundance (same shape)
      taxonomy_top{N}_L{N}_{marker}.tsv  — top-N taxa with 'Other' remainder

    Args:
        taxonomy_qza:      QIIME2 FeatureData[Taxonomy] QZA.
        table_qza:         QIIME2 FeatureTable[Frequency] QZA.
        marker:            Marker key (e.g. '16S', 'MiFish') — controls defaults.
        outdir:            Directory to write output TSVs.
        level:             Taxonomic collapse level (None = use marker default).
        top_n:             Number of top taxa to retain; remainder to 'Other'.
        include:           Taxonomy substrings required to keep a feature.
        exclude:           Taxonomy substrings that remove a feature.
        min_freq:          Minimum total read count to retain a feature.
        min_sample_reads:  Drop samples with fewer reads after filtering.
        dry_run:           If True, log intended actions without writing files.
    """
    cfg = MARKER_CONFIG[marker]
    effective_level = level if level is not None else cfg["default_level"]
    level_name = cfg["level_names"][effective_level - 1]

    log.info("=== 07_taxonomy_table: %s — %s ===", marker, cfg["description"])
    log.info("Collapse level  : %d (%s)", effective_level, level_name)
    log.info("Top-N taxa      : %d", top_n)
    log.info("Min sample reads: %d (drop samples below this after filtering)", min_sample_reads)
    log.info("Output dir      : %s", outdir.resolve())
    if dry_run:
        log.info("DRY RUN — no files will be written")

    # ── 1. Auto-apply marker filter defaults if not overridden ────────────
    defaults = MARKER_FILTER_DEFAULTS.get(marker, {"include": None, "exclude": None})
    if not include and defaults.get("include"):
        include = [s.strip() for s in defaults["include"].split(",") if s.strip()]
        log.info("Auto include filter (marker default): %s", include)
    if not exclude and defaults.get("exclude"):
        exclude = [s.strip() for s in defaults["exclude"].split(",") if s.strip()]
        log.info("Auto exclude filter (marker default): %s", exclude)

    # ── 2. Load ───────────────────────────────────────────────────────────
    log.info("Loading data...")
    tax_df = load_taxonomy_source(taxonomy_qza)
    table  = load_feature_table_qza(table_qza)

    # Align — drop features not in both
    shared = table.index.intersection(tax_df["Feature ID"])
    if len(shared) < len(table):
        log.warning("  %d features in table not in taxonomy — dropped",
                    len(table) - len(shared))
    table  = table.loc[shared]
    tax_df = tax_df[tax_df["Feature ID"].isin(shared)].copy()

    log.info("Before filter: %d features × %d samples", len(table), len(table.columns))
    log.info("Total reads  : %s", f"{int(table.values.sum()):,}")

    # ── 3. Filter ─────────────────────────────────────────────────────────
    if include or exclude or min_freq > 1:
        log.info("Applying filters...")
        table, tax_df = apply_taxonomy_filter(table, tax_df, include, exclude, min_freq)

    log.info("After filter : %d features × %d samples", len(table), len(table.columns))
    log.info("Filtered reads: %s", f"{int(table.values.sum()):,}")

    # ── 4. Per-ASV summary ────────────────────────────────────────────────
    summary = build_asv_summary(table, tax_df, marker)
    save_tsv(summary, outdir / f"taxonomy_summary_{marker}.tsv", dry_run)

    # ── 5. Collapse to level — full taxonomy strings (for programmatic use)
    collapsed_full = collapse_to_level(table, tax_df, marker, effective_level,
                                        use_clean_labels=False)
    save_tsv(collapsed_full,
             outdir / f"taxonomy_counts_L{effective_level}_{marker}.tsv",
             dry_run)

    # ── 6. Collapse to level — clean display labels (for barplots) ────────
    collapsed = collapse_to_level(table, tax_df, marker, effective_level,
                                   use_clean_labels=True)

    # ── 7. Drop samples with insufficient classified reads ─────────────────
    collapsed, dropped = drop_empty_samples(collapsed, min_reads=min_sample_reads)
    if dropped:
        log.warning(
            "  %d sample(s) dropped from barplot table (< %d classified reads): %s",
            len(dropped), min_sample_reads, dropped,
        )

    # ── 8. Relative abundance ─────────────────────────────────────────────
    relabund = to_relative_abundance(collapsed)
    save_tsv(relabund,
             outdir / f"taxonomy_relabund_L{effective_level}_{marker}.tsv",
             dry_run)

    # ── 9. Top-N summary ──────────────────────────────────────────────────
    top_table = top_n_table(relabund, top_n)
    save_tsv(top_table,
             outdir / f"taxonomy_top{top_n}_L{effective_level}_{marker}.tsv",
             dry_run)

    # ── 10. Console summary ───────────────────────────────────────────────
    log.info("")
    log.info("=== Summary ===")
    log.info("Features classified at L%d (%s): %d", effective_level, level_name, len(collapsed))
    log.info("Samples retained (≥%d reads)  : %d / %d",
             min_sample_reads, len(collapsed.columns),
             len(collapsed.columns) + len(dropped))
    log.info("Total classified reads        : %s", f"{int(collapsed.values.sum()):,}")

    top5 = relabund.mean(axis=1).sort_values(ascending=False).head(5)
    log.info("Top 5 taxa at L%d by mean rel. abundance:", effective_level)
    for taxon, ra in top5.items():
        log.info("  %-55s  %.2f%%", taxon, ra * 100)

    log.info("")
    log.info("Output files:")
    for fname in [
        f"taxonomy_summary_{marker}.tsv",
        f"taxonomy_counts_L{effective_level}_{marker}.tsv",
        f"taxonomy_relabund_L{effective_level}_{marker}.tsv",
        f"taxonomy_top{top_n}_L{effective_level}_{marker}.tsv",
    ]:
        log.info("  %s", outdir / fname)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """
Build and return the argument parser for 07_taxonomy_table.py.

    Key arguments: --taxonomy (QZA), --table (QZA), --marker, --outdir.
    Marker-aware defaults for taxonomic level, include/exclude filters, and
    prefix stripping are applied automatically; all can be overridden.
    See module docstring for full rationale and usage examples.
    """
    parser = argparse.ArgumentParser(
        prog="taxonomy_table.py",
        description=(
            "Export QIIME 2 taxonomy to TSV tables for barplotting.\n"
            "Uses unrarefied counts + marker-aware filters + relative abundance.\n"
            "See module docstring for full rationale and usage examples."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--taxonomy", default=None, type=Path,
                   help="FeatureData[Taxonomy] QZA. Derived from --marker if omitted.")
    parser.add_argument("--table", default=None, type=Path,
                   help="FeatureTable[Frequency] QZA (unrarefied). Derived from --marker if omitted.")
    parser.add_argument("--marker", required=True, choices=list(MARKER_CONFIG.keys()),
                   help="Marker gene — controls taxonomy depth and auto-filter defaults.")
    parser.add_argument("--config", default=None, help="Path to pipeline_config.yml.")
    parser.add_argument("--outdir", default=None, type=Path,
                   help="Output directory. Default: results/<marker>/all/taxonomy/ from config.")
    parser.add_argument("--level", type=int, default=None,
                   help=(
                       "Taxonomic level to collapse to (1=Kingdom, N=Species). "
                       "Defaults: 16S/18S/ITS→genus, MiFish/cytb/COI→species."
                   ))
    parser.add_argument("--top-n", type=int, default=30,
                   help="Top-N taxa to include in the summary table. Default: 30")
    parser.add_argument("--include", default="",
                   help=(
                       "Comma-separated taxon strings to keep. "
                       "Overrides the marker default. "
                       "Example: --include Metazoa,Chordata"
                   ))
    parser.add_argument("--exclude", default="",
                   help=(
                       "Comma-separated taxon strings to remove. "
                       "Overrides the marker default. "
                       "Example: --exclude mitochondria,chloroplast"
                   ))
    parser.add_argument("--no-auto-filter", action="store_true",
                   help=(
                       "Disable automatic marker-specific filters entirely. "
                       "Use if you want raw unfiltered output or are applying "
                       "custom --include/--exclude."
                   ))
    parser.add_argument("--min-freq", type=int, default=1,
                   help="Min total reads to keep a feature (ASV-level). Default: 1")
    parser.add_argument("--min-sample-reads", type=int, default=100,
                   help=(
                       "Min classified reads to retain a sample in the "
                       "relabund/barplot table. Samples below this threshold are "
                       "dropped from the output (but noted in log). Default: 100"
                   ))
    parser.add_argument("--dry-run", action="store_true",
                   help="Print planned actions without writing files.")
    return parser


def _resolve_taxonomy_source(paths, cfg, marker: str, explicit: Optional[Path]):
    """
    Decide which taxonomy feeds the count tables (pipeline option b).

    An explicit --taxonomy always wins. Otherwise, if BLAST is enabled in the
    config and a refined taxonomy exists for this marker, prefer it; else use the
    classifier output. Returns (taxonomy_path, kind, classifier_path, refined_path)
    where kind is 'explicit' | 'refined' | 'classifier'.
    """
    classifier = paths.engine_taxonomy_qza(marker, "all")
    refined = paths.engine_blast_results_dir(marker, "all") / f"refined_taxonomy_{marker}.tsv"
    refined_exists = refined.exists()
    if explicit is not None:
        return explicit, "explicit", classifier, (refined if refined_exists else None)
    blast_enabled = cfg.analyses.get("blast", {}).get("enabled", False)
    if blast_enabled and refined_exists:
        return refined, "refined", classifier, refined
    if refined_exists and not blast_enabled:
        log.warning("A BLAST-refined taxonomy exists (%s) but analyses.blast.enabled "
                    "is false — using the classifier taxonomy. Enable BLAST to use it.",
                    refined)
    return classifier, "classifier", classifier, (refined if refined_exists else None)


def _write_taxonomy_provenance(outdir: Path, marker: str, kind: str,
                               classifier_path: Path, refined_path: Optional[Path],
                               cfg, dry_run: bool) -> None:
    """
    Write TAXONOMY_SOURCE_<marker>.md next to the count tables so it is always
    obvious — to anyone, months later — whether the downstream data and stats
    rest on raw classifier output or BLAST-refined taxonomy.
    """
    from datetime import datetime, timezone
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    commit = provenance.git_commit(cfg) or "(not a git checkout)"
    out = outdir / f"TAXONOMY_SOURCE_{marker}.md"

    lines = [f"# Taxonomy provenance — {marker}", "",
             f"- Written: {now}", f"- Pipeline commit: `{commit}`", ""]

    if kind == "refined":
        # Count how many ASVs the refinement actually changed.
        n_total = n_changed = None
        try:
            clf = load_taxonomy_source(classifier_path).set_index("Feature ID")["Taxon"]
            ref = load_taxonomy_source(refined_path).set_index("Feature ID")["Taxon"]
            joined = clf.to_frame("clf").join(ref.to_frame("ref"), how="outer")
            n_total = len(joined)
            n_changed = int((joined["clf"].fillna("") != joined["ref"].fillna("")).sum())
        except Exception as e:  # documentation must never crash the run
            log.warning("Could not diff classifier vs refined taxonomy: %s", e)
        bcfg = cfg.analyses.get("blast", {})
        try:
            mtime = datetime.fromtimestamp(
                Path(refined_path).stat().st_mtime, timezone.utc
            ).strftime("%Y-%m-%d %H:%M UTC")
        except OSError:
            mtime = "(unknown)"
        lines += [
            "## Source: BLAST-REFINED taxonomy", "",
            "> **These count tables — and everything downstream that reads them",
            "> (presence/absence, diversity, figures, stats) — are built on",
            "> BLAST-refined taxonomy, not the raw classifier output.**", "",
            f"- Refined taxonomy used: `{refined_path}`",
            f"- Refined file generated: {mtime}",
            f"- Classifier taxonomy (superseded): `{classifier_path}`",
        ]
        if n_total is not None:
            lines.append(f"- ASVs changed by refinement: **{n_changed} of {n_total}**")
        lines += [
            "", "### BLAST settings (from `analyses.blast`)",
            f"- db: `{bcfg.get('db', '')}`",
            f"- target_rank: {bcfg.get('target_rank', 'genus')}",
            f"- min_pident: {bcfg.get('min_pident', '')}",
            f"- max_target_seqs: {bcfg.get('max_target_seqs', '')}",
            f"- apply: {bcfg.get('apply', False)}",
            "", "### To fall back to the classifier taxonomy",
            "Set `analyses.blast.enabled: false` (or remove the refined file) and",
            "re-run the taxonomy stage. The count tables will be rebuilt from the",
            "classifier output and this note will update accordingly.",
        ]
    elif kind == "explicit":
        lines += ["## Source: explicit --taxonomy override", "",
                  f"- Taxonomy used: `{classifier_path}` (passed on the command line)"]
    else:
        lines += ["## Source: classifier taxonomy", "",
                  f"- Taxonomy used: `{classifier_path}`",
                  "- No BLAST refinement was applied to these tables."]
        if refined_path is not None:
            lines.append(f"- Note: a refined taxonomy exists at `{refined_path}` but "
                         "BLAST is disabled, so it was not used.")

    text = "\n".join(lines) + "\n"
    if dry_run:
        log.info("[dry-run] would write taxonomy provenance: %s", out)
        return
    out.write_text(text)
    log.info("Taxonomy provenance written: %s", out)


def main(argv: Optional[List[str]] = None) -> int:
    """
    Parse arguments and run build_taxonomy_tables for the requested marker.

    Resolves the taxonomy/table/outdir from --marker via config when not given,
    validates the input QZAs exist, merges any include/exclude overrides with the
    marker defaults, builds the tables, then records the run. Returns 0 on
    success, 2 if required inputs are missing.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    paths = get_paths(cfg)

    # Derive inputs/outputs from the marker when not given explicitly.
    # Taxonomy source follows pipeline option (b): prefer a BLAST-refined
    # taxonomy when BLAST is enabled and one exists, else the classifier output.
    taxonomy, tax_source, classifier_path, refined_path = _resolve_taxonomy_source(
        paths, cfg, args.marker, args.taxonomy)
    table    = args.table or paths.engine_table_qza(args.marker, "all", nocontrols=True)
    outdir   = args.outdir or paths.engine_taxonomy_results_dir(args.marker, "all")

    if tax_source == "refined":
        log.warning("=" * 64)
        log.warning("Using BLAST-REFINED taxonomy for %s:", args.marker)
        log.warning("  %s", taxonomy)
        log.warning("Downstream count tables, diversity, and figures will be built")
        log.warning("on refined assignments. See TAXONOMY_SOURCE_%s.md in the output.",
                    args.marker)
        log.warning("=" * 64)

    if not taxonomy.exists():
        log.error("Taxonomy not found: %s", taxonomy)
        return 2
    if not table.exists():
        log.error("Feature table QZA not found: %s", table)
        return 2

    include = [s.strip() for s in args.include.split(",") if s.strip()] \
        if args.include else []
    exclude = [s.strip() for s in args.exclude.split(",") if s.strip()] \
        if args.exclude else []

    if args.no_auto_filter:
        log.info("--no-auto-filter: marker-specific defaults disabled.")
        MARKER_FILTER_DEFAULTS[args.marker] = {"include": None, "exclude": None}

    outdir.mkdir(parents=True, exist_ok=True)

    _write_taxonomy_provenance(outdir, args.marker, tax_source,
                               classifier_path, refined_path, cfg, args.dry_run)

    build_taxonomy_tables(
        taxonomy_qza      = taxonomy,
        table_qza         = table,
        marker            = args.marker,
        outdir            = outdir,
        level             = args.level,
        top_n             = args.top_n,
        include           = include,
        exclude           = exclude,
        min_freq          = args.min_freq,
        min_sample_reads  = args.min_sample_reads,
        dry_run           = args.dry_run,
    )

    if not args.dry_run:
        produced = sorted(str(p) for p in outdir.glob(f"taxonomy_*_{args.marker}.tsv"))
        checkpoint.print_checkpoint(
            cfg, "taxonomy",
            marker=args.marker,
            produced=produced,
            provenance={
                "inputs": {"taxonomy": taxonomy, "table": table},
                "outputs": produced,
                "command": "python " + " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
                "extra": {"level": args.level, "top_n": args.top_n,
                          "taxonomy_source": tax_source},
            },
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
