#!/usr/bin/env python3
"""
09b_clean_diet_table.py
=======================
Clean raw QIIME2 taxonomy count tables from dietary metabarcoding markers
(MiFish 12S, cytb) before running 08b_presence_absence.py.

Performs the following steps in order, logging every decision:

  1. Remove host reads (loon Gaviidae for MiFish; loon Aves for cytb)
  2. Remove negative controls and extraction blanks from sample columns
  3. Flag and optionally remove contamination-suspect samples (e.g., tropical
     taxa in New England samples, extreme outlier read counts)
  4. Remove artefact/misclassification taxa (extinct species, non-target
     mammals, humans, domestic animals) with full logging of what was removed
  5. Collapse multiple OTUs to species level (or lowest resolved level)
     so each species appears once per sample in the output
  6. Remove Unassigned / no-taxonomy rows
  7. Write cleaned count table ready for 08b_presence_absence.py
  8. Write a cleaning report documenting every decision

Marker-specific rules
---------------------
MiFish 12S:
  - Host filter: remove rows where Taxon contains 'Gaviidae'
  - Artefact taxa: extinct species (Pinguinus impennis / great auk),
    non-fish vertebrate families (Canidae, Bovidae, Felidae, Hominidae,
    Cervidae, Phocidae, Gruidae, Laridae non-seabird), Salamandroidea
  - Contamination flag: Siganidae (Indo-Pacific rabbitfish) and Lutjanidae
    (tropical snapper) are not expected in New England loon gut contents --
    flag samples where these are the dominant signal for manual review

cytb:
  - Host filter: remove rows where Taxon contains 'Gaviidae' OR broader
    Aves (cytb amplifies birds broadly -- loon host reads are common)
  - Artefact taxa: domestic mammals (Canidae, Bovidae, Felidae, Hominidae,
    Suidae, Equidae), likely environmental DNA rather than prey
  - Note: non-fish vertebrates (Cervidae, Phocidae, Mustelidae) ARE
    potentially real cytb prey signals -- do NOT auto-remove, flag only

Usage
-----
  # Clean MiFish table
  python 09b_clean_diet_table.py \\
      --counts  results/MiFish/all/taxonomy/taxonomy_counts_L7_MiFish.tsv \\
      --marker  MiFish \\
      --outdir  results/MiFish/all/taxonomy_cleaned/

  # Clean cytb table
  python 09b_clean_diet_table.py \\
      --counts  results/cytb/all/taxonomy/notrim/taxonomy_counts_L7_cytb.tsv \\
      --marker  cytb \\
      --outdir  results/cytb/all/taxonomy_cleaned/

  # Override contamination-suspect samples (comma-separated)
  python 09b_clean_diet_table.py \\
      --counts   results/MiFish/all/taxonomy/taxonomy_counts_L7_MiFish.tsv \\
      --marker   MiFish \\
      --exclude-samples TV250064 \\
      --outdir   results/MiFish/all/taxonomy_cleaned/

  # Dry run -- show what would be removed without writing files
  python 09b_clean_diet_table.py \\
      --counts  results/MiFish/all/taxonomy/taxonomy_counts_L7_MiFish.tsv \\
      --marker  MiFish \\
      --outdir  results/MiFish/all/taxonomy_cleaned/ \\
      --dry-run

Dependencies: pip install pandas numpy
"""
from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Marker-specific filter rules
# ---------------------------------------------------------------------------

# Families always removed as host reads (loon itself)
HOST_FAMILIES: Dict[str, List[str]] = {
    "MiFish": ["Gaviidae"],
    "cytb":   ["Gaviidae"],
}

# Taxon substrings that flag host reads (catches variants above family level)
HOST_TAXON_STRINGS: Dict[str, List[str]] = {
    "MiFish": ["Gaviidae", "Gaviiformes", "Esocidae"],  # Esocidae added: BLAST confirmed Esox ASV 7d61a0be = 100% Gavia immer (loon host DNA miscalled as pike)
    "cytb":   ["Gaviidae", "Gaviiformes", "Aves; Neognathae"],
}

# Families always removed as artefacts / misclassifications
# These are non-prey taxa that should not appear in gut content analysis
ARTEFACT_FAMILIES: Dict[str, List[str]] = {
    "MiFish": [
        "Canidae",       # dog/wolf/fox — domestic or environmental DNA
        "Bovidae",       # goat/cattle — domestic animal contamination
        "Felidae",       # cat — domestic contamination
        "Hominidae",     # human — lab or handler contamination
        "Cervidae",      # deer — not a fish, not MiFish prey
        "Phocidae",      # seal — not expected MiFish prey in gut contents
        "Gruidae",       # crane — bird, not fish prey
        "Lutjanidae",    # tropical/deepwater snapper
        "Apogonidae",    # cardinalfish — Indo-Pacific tropical, BLAST-confirmed artefact (Pristiapogon kallopterus) — BLAST confirmed Pacific species (Etelis coruscans), not NE Atlantic
        # Carangidae family NOT removed — Trachurus (horse mackerel) is a real NE Atlantic fish
        # Only specific subtropical genera are artefacts (confirmed by BLAST):
        "Suidae",        # pig — domestic contamination
        "Equidae",       # horse — domestic contamination
    ],
    "cytb": [
        "Canidae",       # domestic dog/wolf — likely contamination
        "Bovidae",       # domestic cattle — likely contamination
        "Felidae",       # cat — domestic contamination
        "Hominidae",     # human — lab contamination
        "Suidae",        # pig — domestic contamination
        "Equidae",       # horse — domestic contamination
        "Cottidae",      # cytb: BLAST confirmed avian off-target (fd316c6b = raptor genomic DNA, not sculpin)
        # NOTE: Cervidae, Phocidae, Mustelidae NOT removed for cytb --
        # these are real wild vertebrates that could be ingested or
        # detected as secondary prey. Flag them but do not auto-remove.
    ],
}

# Specific species that are known misclassifications
# Key: species name substring; Value: reason
ARTEFACT_SPECIES: Dict[str, str] = {
    "Pinguinus impennis": "Great auk — extinct since 1844, classifier artefact",
    "Pinguinus":          "Great auk genus — extinct since 1844, classifier artefact",
    "Trachinotus":        "Pompano — subtropical Carangidae, BLAST confirmed not NE Atlantic (T. paitensis/carolinus)",
    # NOTE: Esox ASV (7d61a0be) = 100% Gavia immer by BLAST — loon host DNA not filtered by Gaviidae string match
    # NOTE: Brevoortia patronus in classifier = B. tyrannus by BLAST (99.5%) — Atlantic menhaden, relabeled
    # NOTE: Sprattus sprattus in classifier = Clupea harengus by BLAST (99.5%) — Atlantic herring, relabeled
    "Petromyzon_polyA": "Lamprey ASV db1004c6 = poly-A artefact sequence, BLAST no-hit, exclude",
    "Petromyzontidae": "Lamprey — BLAST confirmed raptor/bird genomic DNA (Buteo, Haliaeetus, Aquila); cytb avian off-target amplification",
    "Cottidae_cytb_bird": "Sculpin (cytb) — same ASV as Petromyzontidae misclassification; BLAST confirmed raptor genomic DNA (fd316c6b)",
    "Pekania":   "Fisher — terrestrial mustelid, lab contamination from concurrent fisher study (TV230067, Aug 2023)",
}

# BLAST-verified taxon relabeling
# Classifier assigned wrong species due to reference database gaps.
# Key: substring of classifier taxon string; Value: (replacement string, reason)
# ADAPT: add entries here when BLAST identifies misclassifications in your study.
TAXON_RELABELING: Dict[str, tuple] = {
    "Brevoortia patronus": (
        "k__Metazoa;p__Chordata;c__Actinopteri;o__Clupeiformes;f__Clupeidae;g__Brevoortia;s__tyrannus",
        "BLAST: 99.5% B. tyrannus (Atlantic menhaden); B. patronus is Gulf of Mexico species"
    ),
    "Sprattus sprattus": (
        "k__Metazoa;p__Chordata;c__Actinopteri;o__Clupeiformes;f__Clupeidae;g__Clupea;s__harengus",
        "BLAST: 99.5% Clupea harengus (Atlantic herring); Sprattus is Eastern Atlantic only"
    ),
    # Esox ASV 7d61a0be = 100% Gavia immer by BLAST — handled in host filter below
}

# Families to flag as potentially unexpected but not auto-remove
# These generate a WARNING in the log and are noted in the report
FLAG_FAMILIES: Dict[str, List[str]] = {
    "MiFish": [
        "Siganidae",       # rabbitfish — Indo-Pacific, not NE Atlantic
        "Lutjanidae",      # tropical snapper — not expected in NE loon diet
        "Apogonidae",      # cardinalfish — tropical
        "Eleginopidae",    # Patagonian toothfish relatives — southern hemisphere
    ],
    "cytb": [
        "Cervidae",        # deer — real but worth flagging for review
        "Phocidae",        # seal — real but worth flagging for review
    ],
}

# Control sample prefixes -- these columns are dropped from output
CONTROL_PREFIXES = ("NTC", "PAC", "XB", "ntc", "pac", "xb")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def extract_sample_cols(df: pd.DataFrame, marker: str) -> List[str]:
    """
    Return sample columns -- TV-prefixed columns only, dropping controls.
    Also drops any column that contains only the marker suffix without a TV ID.
    """
    return [
        c for c in df.columns
        if c.startswith("TV") and not any(
            c.startswith(p) for p in CONTROL_PREFIXES
        )
    ]


def extract_species(taxon_str: str) -> Optional[str]:
    """
    Extract the lowest resolved taxonomic name from a full NCBI lineage string.
    Returns the last semicolon-delimited token that looks like a species binomial
    or genus name, not a higher-rank descriptor.
    """
    if not isinstance(taxon_str, str) or taxon_str in ("Unassigned", ""):
        return None
    parts = [p.strip() for p in taxon_str.split(";")]
    # Walk from end -- find last part that isn't a rank descriptor
    skip_patterns = re.compile(
        r"(incertae sedis|cellular organisms|Opisthokonta|Eumetazoa|"
        r"Bilateria|Deuterostomia|Vertebrata|Gnathostomata|Teleostomi|"
        r"Euteleostomi|Actinopterygii|Teleostei|Neopterygii|unclassified)",
        re.IGNORECASE
    )
    for part in reversed(parts):
        if part and not skip_patterns.search(part):
            return part
    return parts[-1] if parts else None


def get_taxon_col(df: pd.DataFrame) -> str:
    """
    Auto-detect the taxonomy string column.
    08_taxonomy_table.py outputs 'Species' as the first column.
    Raw CSV files may use 'Taxon'. Fall back to first column.
    """
    for candidate in ("Species", "Taxon"):
        if candidate in df.columns:
            return candidate
    return df.columns[0]


def strip_rank_prefix(part: str) -> str:
    """
    Strip Greengenes-style rank prefixes (k__, p__, c__, o__, f__, g__, s__).
    'f__Hominidae' → 'Hominidae', 's__sapiens' → 'sapiens'
    """
    return re.sub(r'^[kpcofgsr]__', '', part).strip()


def get_family(row: pd.Series, taxon_col: str = "Species") -> Optional[str]:
    """
    Get family from the Family column, or parse from the taxonomy string.
    Handles both NCBI style ('Hominidae') and Greengenes style ('f__Hominidae').
    """
    if pd.notna(row.get("Family")) and str(row.get("Family", "")).strip():
        return str(row["Family"]).strip()
    taxon = str(row.get(taxon_col, ""))
    if not taxon:
        return None
    # Search for family name — with or without rank prefix
    m = re.search(r'\b(?:f__)?(\w+idae|\w+inae)\b', taxon)
    if m:
        return m.group(1)
    return None


# ---------------------------------------------------------------------------
# Cleaning steps
# ---------------------------------------------------------------------------

def step_remove_unassigned(df: pd.DataFrame, taxon_col: str) -> Tuple[pd.DataFrame, List[str]]:
    """Remove rows with no taxonomy assignment."""
    mask = (df[taxon_col] == "Unassigned") | df[taxon_col].isna()
    n = mask.sum()
    log.info("  Unassigned rows removed: %d", n)
    return df[~mask].copy(), [f"Unassigned ({n} rows)"]


def step_remove_host(
    df: pd.DataFrame,
    marker: str,
    taxon_col: str,
) -> Tuple[pd.DataFrame, List[str]]:
    """Remove host (loon) reads. Handles both NCBI and Greengenes k__/c__ prefix formats."""
    strings = HOST_TAXON_STRINGS.get(marker, [])
    # Also add prefix variants for Greengenes-style lineages
    prefix_variants = []
    for s in strings:
        prefix_variants.append(s)
        # e.g. "Gaviidae" also matches "f__Gaviidae"; "Aves" also matches "c__Aves"
        prefix_variants.append(f"__{s}")
    mask = pd.Series(False, index=df.index)
    for s in prefix_variants:
        mask |= df[taxon_col].str.contains(s, na=False)
    n_rows = mask.sum()
    n_reads = df[mask]["Total"].sum() if "Total" in df.columns else 0
    log.info("  Host reads removed: %d rows, %s reads", n_rows, f"{n_reads:,}")
    return df[~mask].copy(), [f"Host (Gaviidae/Aves) — {n_rows} rows, {n_reads:,} reads"]


def step_remove_artefacts(
    df: pd.DataFrame,
    marker: str,
    taxon_col: str,
) -> Tuple[pd.DataFrame, List[str]]:
    """Remove known artefact and misclassification taxa."""
    removed_log: List[str] = []
    mask = pd.Series(False, index=df.index)

    # Family-level removals — handles both plain and k__/f__ prefix formats
    artefact_fams = set(ARTEFACT_FAMILIES.get(marker, []))
    for idx, row in df.iterrows():
        fam = get_family(row, taxon_col)
        if fam in artefact_fams:
            mask[idx] = True

    # Species-level removals
    for species_str, reason in ARTEFACT_SPECIES.items():
        species_mask = df[taxon_col].str.contains(species_str, na=False)
        if species_mask.any():
            log.warning("  Removing extinct/artefact species '%s': %s",
                        species_str, reason)
            removed_log.append(f"{species_str}: {reason}")
        mask |= species_mask

    n = mask.sum()
    if n > 0:
        for idx in df[mask].index:
            fam = get_family(df.loc[idx], taxon_col)
            sp  = extract_species(str(df.loc[idx][taxon_col]))
            reads = df.loc[idx]["Total"] if "Total" in df.columns else "?"
            log.info("    Removing artefact: %s / %s (%s reads)", fam, sp, reads)
            removed_log.append(f"{fam} / {sp} — {reads} reads")

    log.info("  Artefact rows removed: %d", n)
    return df[~mask].copy(), removed_log


def step_flag_unexpected(
    df: pd.DataFrame,
    marker: str,
    sample_cols: List[str],
    taxon_col: str = "Species",
) -> List[str]:
    """Flag families that are unexpected but not auto-removed."""
    warnings: List[str] = []
    flag_fams = set(FLAG_FAMILIES.get(marker, []))
    for idx, row in df.iterrows():
        fam = get_family(row, taxon_col)
        if fam in flag_fams:
            detections = [c for c in sample_cols if c in df.columns and row.get(c, 0) > 0]
            sp = extract_species(str(row.get(taxon_col, "")))
            reads = row.get("Total", "?")
            msg = (f"UNEXPECTED FAMILY {fam} / {sp}: "
                   f"{reads} total reads in {len(detections)} samples — "
                   f"review before including in analysis")
            log.warning("  ⚠  %s", msg)
            warnings.append(msg)
    return warnings


def step_exclude_samples(
    df: pd.DataFrame,
    exclude: List[str],
    sample_cols: List[str],
) -> Tuple[pd.DataFrame, List[str]]:
    """Remove specified samples from the count table."""
    removed: List[str] = []
    for sample in exclude:
        # Match by exact name or by partial TV ID
        matches = [c for c in sample_cols if sample in c]
        for col in matches:
            if col in df.columns:
                total_reads = df[col].sum()
                df = df.drop(columns=[col])
                log.warning("  Excluded sample: %s (%s reads)", col, f"{total_reads:,}")
                removed.append(f"{col} ({total_reads:,} reads) — user-excluded")
        if not matches:
            log.warning("  --exclude-samples: '%s' not found in columns", sample)
    return df, removed


def step_collapse_to_species(
    df: pd.DataFrame,
    sample_cols: List[str],
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Collapse multiple OTUs that resolve to the same species (or lowest
    taxonomic level) into a single row by summing read counts.

    Multiple ASVs from the same species are real haplotype variation but
    should count as a single detection for presence/absence analysis.
    The collapsed row retains the most abundant OTU's sequence and OTU ID.
    """
    log.info("  Collapsing OTUs to species level...")
    notes: List[str] = []

    # Assign a collapse key: Species column if populated, else Genus, else Family,
    # else lowest resolvable name from Taxon string
    def collapse_key(row: pd.Series) -> str:
        if pd.notna(row.get("Species")) and str(row["Species"]).strip():
            return str(row["Species"]).strip()
        if pd.notna(row.get("Genus")) and str(row["Genus"]).strip():
            return str(row["Genus"]).strip()
        if pd.notna(row.get("Family")) and str(row["Family"]).strip():
            return str(row["Family"]).strip()
        return extract_species(str(row.get(taxon_col, ""))) or "Unknown"

    df = df.copy()
    df["_collapse_key"] = df.apply(collapse_key, axis=1)

    groups = df.groupby("_collapse_key")
    collapsed_rows = []
    n_collapsed = 0

    for key, group in groups:
        if len(group) == 1:
            collapsed_rows.append(group.iloc[0].drop("_collapse_key"))
            continue

        # Multiple OTUs for same taxon — sum reads, keep dominant OTU metadata
        n_collapsed += len(group) - 1
        dominant = group.loc[group["Total"].idxmax()]

        new_row = dominant.copy().drop("_collapse_key")
        for col in sample_cols:
            if col in group.columns:
                new_row[col] = group[col].sum()
        if "Total" in group.columns:
            new_row["Total"] = group["Total"].sum()

        collapsed_rows.append(new_row)
        msg = (f"Collapsed {len(group)} OTUs → 1 for '{key}' "
               f"(total {int(new_row.get('Total', 0)):,} reads)")
        log.info("    %s", msg)
        notes.append(msg)

    df_out = pd.DataFrame(collapsed_rows)
    log.info("  OTU collapse: %d rows removed, %d unique taxa remain",
             n_collapsed, len(df_out))
    return df_out, notes


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_cleaning(
    counts_path: Path,
    marker: str,
    outdir: Path,
    exclude_samples: List[str],
    dry_run: bool,
) -> None:
    """Full cleaning pipeline."""
    safe_mkdir(outdir)

    log.info("=== 09b_clean_diet_table: %s ===", marker)
    log.info("Input : %s", counts_path)
    log.info("Outdir: %s", outdir)
    if dry_run:
        log.info("DRY RUN — no files will be written")

    # ── Load ──────────────────────────────────────────────────────────────
    # Handle both TSV and CSV inputs
    sep = "\t" if counts_path.suffix in (".tsv", ".txt") else ","
    df = pd.read_csv(counts_path, sep=sep, dtype=str)
    df.columns = df.columns.str.strip()

    # Standardise numeric columns
    numeric_skip = {"Sequence", "Taxon", "OTU", "Species", "Genus", "Family"}
    for col in df.columns:
        if col not in numeric_skip:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    log.info("Loaded: %d rows", len(df))

    # Identify sample columns (TV-prefixed, non-control)
    all_sample_cols = extract_sample_cols(df, marker)
    log.info("Sample columns: %d TV birds", len(all_sample_cols))

    # ── Step 0: drop control columns ──────────────────────────────────────
    control_cols = [c for c in df.columns if any(
        c.startswith(p) for p in CONTROL_PREFIXES)]
    if control_cols:
        log.info("  Dropping %d control columns: %s", len(control_cols), control_cols)
        df = df.drop(columns=control_cols, errors="ignore")

    # Auto-detect taxonomy string column ('Species' from 08_taxonomy_table.py, 'Taxon' from raw CSV)
    taxon_col = get_taxon_col(df)
    log.info("Taxonomy column: '%s'", taxon_col)

    report_lines: List[str] = [
        f"=== 09b_clean_diet_table cleaning report: {marker} ===",
        f"Input file  : {counts_path}",
        f"Input rows  : {len(df)}",
        f"Marker      : {marker}",
        f"Taxon column: {taxon_col}",
        "",
    ]

    # ── Step 1: remove unassigned ─────────────────────────────────────────
    log.info("── Step 1: remove unassigned ─────────────────────────────────")
    df, removed = step_remove_unassigned(df, taxon_col)
    report_lines.append("Step 1 — Unassigned removed:")
    report_lines += [f"  {r}" for r in removed]
    report_lines.append(f"  Rows remaining: {len(df)}")
    report_lines.append("")

    # ── Step 2: remove host reads ──────────────────────────────────────────
    log.info("── Step 2: remove host reads ──────────────────────────────────")
    df, removed = step_remove_host(df, marker, taxon_col)
    report_lines.append("Step 2 — Host reads removed:")
    report_lines += [f"  {r}" for r in removed]
    report_lines.append(f"  Rows remaining: {len(df)}")
    report_lines.append("")

    # ── Step 3: remove artefact taxa ──────────────────────────────────────
    log.info("── Step 3: remove artefact taxa ───────────────────────────────")
    df, removed = step_remove_artefacts(df, marker, taxon_col)
    report_lines.append("Step 3 — Artefact/misclassification taxa removed:")
    if removed:
        report_lines += [f"  {r}" for r in removed]
    else:
        report_lines.append("  (none)")
    report_lines.append(f"  Rows remaining: {len(df)}")
    report_lines.append("")

    # ── Step 4: flag unexpected families ──────────────────────────────────
    log.info("── Step 4: flag unexpected families ──────────────────────────")
    current_sample_cols = [c for c in all_sample_cols if c in df.columns]
    warnings = step_flag_unexpected(df, marker, current_sample_cols, taxon_col)
    report_lines.append("Step 4 — Unexpected families flagged (NOT removed — review):")
    if warnings:
        report_lines += [f"  ⚠  {w}" for w in warnings]
    else:
        report_lines.append("  (none)")
    report_lines.append("")

    # ── Step 5: exclude user-specified samples ─────────────────────────────
    if exclude_samples:
        log.info("── Step 5: exclude specified samples ─────────────────────")
        df, removed = step_exclude_samples(df, exclude_samples, current_sample_cols)
        report_lines.append("Step 5 — User-excluded samples:")
        report_lines += [f"  {r}" for r in removed]
    else:
        report_lines.append("Step 5 — No samples excluded by user.")
    report_lines.append("")

    # ── Step 6: collapse to species ────────────────────────────────────────
    log.info("── Step 6: collapse OTUs to species level ─────────────────────")
    current_sample_cols = [c for c in all_sample_cols if c in df.columns]
    df, notes = step_collapse_to_species(df, current_sample_cols)
    report_lines.append("Step 6 — OTU collapse to species level:")
    if notes:
        report_lines += [f"  {n}" for n in notes]
    else:
        report_lines.append("  No collapsing needed (all taxa already unique).")
    report_lines.append(f"  Final unique taxa: {len(df)}")
    report_lines.append("")

    # ── Summary ────────────────────────────────────────────────────────────
    current_sample_cols = [c for c in all_sample_cols if c in df.columns]
    n_samples_with_prey = int((df[current_sample_cols].sum(axis=0) > 0).sum())
    total_prey_reads = int(df[current_sample_cols].values.sum())

    report_lines += [
        "=== Final Summary ===",
        f"  Final taxa rows    : {len(df)}",
        f"  Samples with prey  : {n_samples_with_prey} / {len(current_sample_cols)}",
        f"  Total prey reads   : {total_prey_reads:,}",
        "",
        "Next step — run presence/absence analysis:",
        f"  python 08b_presence_absence.py \\",
        f"    --counts   {outdir / f'taxonomy_counts_cleaned_{marker}.tsv'} \\",
        f"    --metadata metadata/qiime/metadata_{marker}.tsv \\",
        f"    --marker   {marker} \\",
        f"    --group-by Group \\",
        f"    --sample-label loon \\",
        f"    --outdir   results/{marker}/all/presence_absence/",
    ]

    report_text = "\n".join(report_lines) + "\n"
    log.info("")
    log.info("=== Summary ===")
    log.info("Final taxa: %d | Samples with prey: %d/%d | Total reads: %s",
             len(df), n_samples_with_prey, len(current_sample_cols),
             f"{total_prey_reads:,}")

    # ── Write outputs ──────────────────────────────────────────────────────
    if not dry_run:
        out_counts = outdir / f"taxonomy_counts_cleaned_{marker}.tsv"
        out_report = outdir / f"cleaning_report_{marker}.txt"

        # Write count table (taxa as rows, samples as columns)
        # Index on the collapse key / species name for readability
        if "_collapse_key" in df.columns:
            df = df.drop(columns=["_collapse_key"], errors="ignore")

        df.to_csv(out_counts, sep="\t", index=False)
        log.info("Written: %s", out_counts)

        out_report.write_text(report_text, encoding="utf-8")
        log.info("Written: %s", out_report)
    else:
        print("\n" + report_text)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="09b_clean_diet_table.py",
        description=(
            "Clean MiFish or cytb taxonomy count tables before presence/absence "
            "analysis. Removes host reads, artefact taxa, and collapses OTUs to "
            "species level. Run before 08b_presence_absence.py."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--counts", required=True, type=Path,
        help="Taxonomy count TSV or CSV from 08_taxonomy_table.py.",
    )
    p.add_argument(
        "--marker", required=True, choices=["MiFish", "cytb"],
        help="Marker name — determines which cleaning rules are applied.",
    )
    p.add_argument(
        "--outdir", required=True, type=Path,
        help="Directory for cleaned outputs.",
    )
    p.add_argument(
        "--exclude-samples", default=None,
        help=(
            "Comma-separated TV IDs to exclude (e.g. TV250064). "
            "Use for contamination-suspect samples identified during QC. "
            "Partial match: 'TV250064' matches 'TV250064-GI-Mifish'."
        ),
    )
    p.add_argument(
        "--dry-run", action="store_true", default=False,
        help="Print planned actions and report without writing files.",
    )
    return p


def main() -> int:
    args = build_parser().parse_args()

    if not args.counts.exists():
        log.error("Input file not found: %s", args.counts)
        return 2

    exclude = []
    if args.exclude_samples:
        exclude = [s.strip() for s in args.exclude_samples.split(",")]

    try:
        run_cleaning(
            counts_path    = args.counts,
            marker         = args.marker,
            outdir         = args.outdir,
            exclude_samples = exclude,
            dry_run        = args.dry_run,
        )
    except Exception as e:
        log.error("Cleaning failed: %s", e)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
