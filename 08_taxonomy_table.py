#!/usr/bin/env python3
"""
08_taxonomy_table.py
====================
Export taxonomy classifications to human-readable TSV tables for any marker
in the metabarcoding pipeline.

Philosophy
----------
Diversity statistics (alpha/beta) are computed on a rarefied feature table
that includes all reads — host DNA and all. This is correct for those metrics.

Taxonomy barplots are different: rarefying a table dominated by host reads
(e.g. loon gut 16S) before filtering means almost all subsampled reads are
Unassigned, leaving only a handful of bacterial reads per sample. The result
is barplots that are nearly empty — not a bug, but the wrong order of operations.

This script takes the correct approach:
  1. Start from the UNRAREFIED feature table.
  2. Apply marker-aware filters (remove host, off-target kingdoms, etc.).
  3. Collapse to genus (or specified level).
  4. Compute relative abundance among the remaining classified reads.
  5. Export: counts TSV, relative-abundance TSV, and top-N summary.

The relative-abundance TSV is the direct input to barplot scripts.
Methods note for manuscript: "Taxonomic composition was visualized using
unrarefied relative abundance among [target]-assigned reads."

Marker filter defaults (applied automatically unless --include/--exclude set):
  Marker   Include          Exclude
  -------  ---------------  -----------------------------------------
  16S      (bacteria only)  mitochondria,chloroplast,Eukaryota,Archaea
  18S      (all eukaryotes) Bacteria,Archaea
  ITS      Fungi            Bacteria,Archaea
  MiFish   Actinopteri      Bacteria,Archaea  (MitoFish QIIME DB uses QIIME prefixes)
  cytb     Vertebrata       Bacteria,Viruses,Archaea
  COI      Metazoa          Bacteria,Viruses,Archaea,Viridiplantae,Fungi

Outputs (all in --outdir):
  taxonomy_counts_L{N}_{marker}.tsv     — collapsed count table
  taxonomy_relabund_L{N}_{marker}.tsv   — relative abundance table (0–1)
  taxonomy_top{N}_L{N}_{marker}.tsv     — top-N taxa by mean relative abundance
  taxonomy_summary_{marker}.tsv         — per-ASV summary with parsed levels

Usage examples:

  # 16S — genus level, auto bacteria filter
  python 08_taxonomy_table.py \\
      --taxonomy  qiime2/taxonomy/taxonomy_16S_silva138.qza \\
      --table     qiime2/dada2/table_16S_DvT.qza \\
      --marker    16S \\
      --outdir    results/16S/DvT/taxonomy/

  # MiFish — species level, auto Vertebrata filter
  python 08_taxonomy_table.py \\
      --taxonomy  qiime2/MiFish/all/taxonomy/taxonomy.qza \\
      --table     qiime2/MiFish/all/dada2/table_filtered.qza \\
      --marker    MiFish \\
      --outdir    results/MiFish/all/taxonomy/

  # ITS — family level (override default genus), Fungi only
  python 08_taxonomy_table.py \\
      --taxonomy  qiime2/taxonomy/taxonomy_ITS.qza \\
      --table     qiime2/dada2/table_ITS.qza \\
      --marker    ITS \\
      --level     5 \\
      --outdir    results/ITS/all/taxonomy/

  # Override auto-filters completely
  python 08_taxonomy_table.py \\
      --taxonomy  qiime2/taxonomy/taxonomy_COI.qza \\
      --table     qiime2/dada2/table_COI.qza \\
      --marker    COI \\
      --include   Arthropoda \\
      --exclude   Insecta \\
      --outdir    results/COI/all/taxonomy/

  # Dry run — see what would happen without writing files
  python 08_taxonomy_table.py \\
      --taxonomy  qiime2/taxonomy/taxonomy_16S_silva138.qza \\
      --table     qiime2/dada2/table_16S_DvT.qza \\
      --marker    16S \\
      --outdir    results/16S/DvT/taxonomy/ \\
      --dry-run
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
    log.info("  Taxonomy: %d features from %s", len(df), qza_path.name)
    return df


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

    t = str(taxon).strip()
    parts = [p.strip() for p in t.split(sep)]

    def clean(s: str) -> str:
        if do_strip:
            s = re.sub(r"^[a-z]__", "", s)
        return s.strip().strip("_")

    cleaned = [clean(p) for p in parts]

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
            f = cleaned[level - 2] if len(cleaned) >= level - 1 else ""
            if f and f not in ("", "__"):
                f = f.replace("_", " ")
                return f"uncl. {f}"
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

        def truncate(taxon: str) -> str:
            if pd.isna(taxon) or not str(taxon).strip():
                return "Unclassified"
            t = str(taxon).strip()
            if do_strip:
                t = re.sub(r"[a-z]__", "", t)
            parts = [p.strip() for p in t.split(sep)][:level]
            while len(parts) < level:
                parts.append("Unclassified")
            return sep.join(parts)

        label_series = tax_idx.reindex(table.index).apply(truncate)

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

    def parse_row(taxon):
        if pd.isna(taxon) or str(taxon).strip().lower() in ("", "unassigned"):
            return ["Unclassified"] * len(level_names)
        t = str(taxon).strip()
        if do_strip:
            t = re.sub(r"[a-z]__", "", t)
        parts = [p.strip() for p in t.split(sep)][:len(level_names)]
        while len(parts) < len(level_names):
            parts.append("Unclassified")
        return [p or "Unclassified" for p in parts]

    totals   = table.sum(axis=1).rename("total_reads")
    presence = (table > 0).sum(axis=1).rename("samples_present")

    level_df = pd.DataFrame(
        [parse_row(tax_df.set_index("Feature ID")["Taxon"].get(f, ""))
         for f in table.index],
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
    cfg = MARKER_CONFIG[marker]
    effective_level = level if level is not None else cfg["default_level"]
    level_name = cfg["level_names"][effective_level - 1]

    log.info("=== 08_taxonomy_table: %s — %s ===", marker, cfg["description"])
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
    tax_df = load_taxonomy_qza(taxonomy_qza)
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
    p = argparse.ArgumentParser(
        prog="08_taxonomy_table.py",
        description=(
            "Export QIIME 2 taxonomy to TSV tables for barplotting.\n"
            "Uses unrarefied counts + marker-aware filters + relative abundance.\n"
            "See module docstring for full rationale and usage examples."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--taxonomy", required=True, type=Path,
                   help="FeatureData[Taxonomy] QZA.")
    p.add_argument("--table", required=True, type=Path,
                   help="FeatureTable[Frequency] QZA (unrarefied recommended).")
    p.add_argument("--marker", required=True, choices=list(MARKER_CONFIG.keys()),
                   help="Marker gene — controls taxonomy depth and auto-filter defaults.")
    p.add_argument("--outdir", default="taxonomy_tables", type=Path,
                   help="Output directory. Default: ./taxonomy_tables")
    p.add_argument("--level", type=int, default=None,
                   help=(
                       "Taxonomic level to collapse to (1=Kingdom, N=Species). "
                       "Defaults: 16S/18S/ITS→genus, MiFish/cytb/COI→species."
                   ))
    p.add_argument("--top-n", type=int, default=30,
                   help="Top-N taxa to include in the summary table. Default: 30")
    p.add_argument("--include", default="",
                   help=(
                       "Comma-separated taxon strings to keep. "
                       "Overrides the marker default. "
                       "Example: --include Metazoa,Chordata"
                   ))
    p.add_argument("--exclude", default="",
                   help=(
                       "Comma-separated taxon strings to remove. "
                       "Overrides the marker default. "
                       "Example: --exclude mitochondria,chloroplast"
                   ))
    p.add_argument("--no-auto-filter", action="store_true",
                   help=(
                       "Disable automatic marker-specific filters entirely. "
                       "Use if you want raw unfiltered output or are applying "
                       "custom --include/--exclude."
                   ))
    p.add_argument("--min-freq", type=int, default=1,
                   help="Min total reads to keep a feature (ASV-level). Default: 1")
    p.add_argument("--min-sample-reads", type=int, default=100,
                   help=(
                       "Min classified reads to retain a sample in the "
                       "relabund/barplot table. Samples below this threshold are "
                       "dropped from the output (but noted in log). Default: 100"
                   ))
    p.add_argument("--dry-run", action="store_true",
                   help="Print planned actions without writing files.")
    return p


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if not args.taxonomy.exists():
        log.error("Taxonomy QZA not found: %s", args.taxonomy)
        return 2
    if not args.table.exists():
        log.error("Feature table QZA not found: %s", args.table)
        return 2

    include = [s.strip() for s in args.include.split(",") if s.strip()] \
        if args.include else []
    exclude = [s.strip() for s in args.exclude.split(",") if s.strip()] \
        if args.exclude else []

    # --no-auto-filter: pass sentinel values so auto-defaults are skipped
    if args.no_auto_filter:
        # We signal "user explicitly wants no filter" by setting include/exclude
        # to a non-empty sentinel that won't match anything meaningful, then
        # just skip the auto-apply block.  Simplest: set a flag on the marker.
        log.info("--no-auto-filter: marker-specific defaults disabled.")
        # Monkey-patch: temporarily empty the defaults for this marker
        MARKER_FILTER_DEFAULTS[args.marker] = {"include": None, "exclude": None}

    args.outdir.mkdir(parents=True, exist_ok=True)

    build_taxonomy_tables(
        taxonomy_qza      = args.taxonomy,
        table_qza         = args.table,
        marker            = args.marker,
        outdir            = args.outdir,
        level             = args.level,
        top_n             = args.top_n,
        include           = include,
        exclude           = exclude,
        min_freq          = args.min_freq,
        min_sample_reads  = args.min_sample_reads,
        dry_run           = args.dry_run,
    )

    log.info("=== Done ===")
    log.info(
        "\nNext step — generate barplots with 09_plot_taxonomy.py:\n"
        "  python 09_plot_taxonomy.py \\\n"
        "    --relabund  %s/taxonomy_relabund_L*_%s.tsv \\\n"
        "    --metadata  metadata/qiime/metadata_%s.tsv \\\n"
        "    --group-by  Group \\\n"
        "    --marker    %s \\\n"
        "    --palette   purple \\\n"
        "    --outdir    %s/../figures/taxonomy/",
        str(args.outdir), args.marker, args.marker, args.marker, str(args.outdir),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
