#!/usr/bin/env python3
"""
09c_blast_qc_unclassified.py
=============================
BLAST-based QC for poorly-classified ASVs from QIIME2 metabarcoding pipelines.

Identifies ASVs that the classifier could not resolve below a specified
taxonomic rank (default: class), extracts their sequences, runs BLASTn
against a local or remote NCBI nt database, resolves taxids to species
names via NCBI E-utilities, and generates a flagged report for manual
review before cleaning (11_clean_diet_table.py).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PIPELINE POSITION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  dada2/rep-seqs.qza  ──┐
  taxonomy/taxonomy.qza ─┼──► 09c_blast_qc_unclassified.py
                         │         ↓
                    blast_qc_report_{marker}.txt   ← human review
                    confirmed_artefacts_{marker}.txt ← pass to 09b
                         ↓
                  11_clean_diet_table.py  (--exclude-asvs flag)
                         ↓
                  downstream figures / diversity

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WHY THIS EXISTS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Short-amplicon classifiers (MiFish 12S ~170bp, cytb ~300bp) sometimes
  assign host DNA or environmental bacteria to broad fish classes (e.g.
  'Actinopteri Unclassified') when the reference database lacks a close
  match. String-based host filters in 09b_ cannot catch these because the
  taxonomy string gives no indication of the true identity.

  This script catches those cases by directly BLASTing residual
  unclassified ASVs against nt, so misclassified host reads (e.g. Gavia
  immer miscalled as Actinopteri) and bacterial off-targets are identified
  before they inflate diversity metrics or appear in barplots.

  ADAPT FOR YOUR STUDY:
  - Change --target-rank to match the level you expect classification
    to stop at for your marker (default: class for MiFish/cytb)
  - Change --min-confidence to match your classifier threshold
  - Change --conflict-taxa to the host taxa expected in your study

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
USAGE EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # MiFish — local nt database
  python scripts/09c_blast_qc_unclassified.py \\
      --taxonomy  qiime2/MiFish/all/taxonomy-exported/taxonomy.tsv \\
      --rep-seqs  qiime2/MiFish/all/rep-seqs-exported/dna-sequences.fasta \\
      --marker    MiFish \\
      --blast-db  /home/share/databases/ncbi_nt/nt \\
      --outdir    qiime2/MiFish/all/blast_qc/

  # cytb — local nt database, stricter confidence threshold
  python scripts/09c_blast_qc_unclassified.py \\
      --taxonomy  qiime2/cytb/all/taxonomy-exported/taxonomy.tsv \\
      --rep-seqs  qiime2/cytb/all/rep-seqs-exported/dna-sequences.fasta \\
      --marker    cytb \\
      --blast-db  /home/share/databases/ncbi_nt/nt \\
      --min-confidence 0.80 \\
      --outdir    qiime2/cytb/all/blast_qc/

  # Remote BLAST (no local nt database — slower)
  python scripts/09c_blast_qc_unclassified.py \\
      --taxonomy  qiime2/MiFish/all/taxonomy-exported/taxonomy.tsv \\
      --rep-seqs  qiime2/MiFish/all/rep-seqs-exported/dna-sequences.fasta \\
      --marker    MiFish \\
      --remote \\
      --outdir    qiime2/MiFish/all/blast_qc/

  # Dry run — just show which ASVs would be BLASTed
  python scripts/09c_blast_qc_unclassified.py \\
      --taxonomy  qiime2/MiFish/all/taxonomy-exported/taxonomy.tsv \\
      --rep-seqs  qiime2/MiFish/all/rep-seqs-exported/dna-sequences.fasta \\
      --marker    MiFish \\
      --blast-db  /home/share/databases/ncbi_nt/nt \\
      --outdir    qiime2/MiFish/all/blast_qc/ \\
      --dry-run

Dependencies:
  - BLAST+ (blastn must be on PATH)
  - NCBI E-utilities (efetch, xtract must be on PATH for taxid resolution)
  - pip install pandas
"""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
import urllib.request
import urllib.parse
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
# Rank order — used to determine whether an ASV is "unresolved"
# ADAPT: extend if your marker resolves to a different level
# ---------------------------------------------------------------------------
RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

RANK_PREFIXES = {
    "k__": "kingdom",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}

# ---------------------------------------------------------------------------
# Taxa that signal a conflict with expected marker targets
# ADAPT: change these to match your study organism and marker
# ---------------------------------------------------------------------------
CONFLICT_TAXA_DEFAULT = [
    # Loon host (Gaviiformes)
    "Gavia", "Gaviidae", "Gaviiformes", "Gaviimorphae",
    # Bacteria — should not appear in dietary metabarcoding
    "Bacteria", "Proteobacteria", "Firmicutes", "Enterobacteriaceae",
    "Rahnella", "Enterococcus", "Escherichia", "uncultured bacterium",
    "Cetobacterium",
    # Raptors / off-target birds (cytb known issue)
    "Buteo", "Haliaeetus", "Aquila", "Accipitridae",
]


# ---------------------------------------------------------------------------
# Taxonomy parsing
# ---------------------------------------------------------------------------

def parse_deepest_rank(taxonomy_str: str) -> Tuple[str, str]:
    """
    Return (deepest_rank_name, deepest_taxon_value) for a QIIME2 taxonomy
    string like 'k__Metazoa;p__Chordata;c__Actinopteri'.

    Returns ('unassigned', '') if the string is empty or Unassigned.
    """
    if not taxonomy_str or taxonomy_str.strip().lower() in ("unassigned", ""):
        return "unassigned", ""

    parts = [p.strip() for p in taxonomy_str.split(";")]
    deepest_rank = "kingdom"
    deepest_val  = ""

    for part in parts:
        for prefix, rank in RANK_PREFIXES.items():
            if part.startswith(prefix):
                val = part[len(prefix):].strip()
                if val and val.lower() not in ("", "unclassified", "uncultured",
                                               "uncl.", "x", "nan"):
                    deepest_rank = rank
                    deepest_val  = val
                break

    return deepest_rank, deepest_val


def find_unresolved_asvs(
    taxonomy_df: pd.DataFrame,
    target_rank: str,
    min_confidence: float,
) -> pd.DataFrame:
    """
    Return rows where the deepest classified rank does not reach target_rank
    AND confidence is below min_confidence.

    taxonomy_df columns: Feature ID, Taxon, Confidence
    """
    target_idx = RANK_ORDER.index(target_rank)
    unresolved = []

    for _, row in taxonomy_df.iterrows():
        asv_id   = row["Feature ID"]
        taxon    = str(row.get("Taxon", ""))
        try:
            conf = float(row.get("Confidence", 0.0))
        except (ValueError, TypeError):
            conf = 0.0

        deepest_rank, deepest_val = parse_deepest_rank(taxon)

        if deepest_rank == "unassigned":
            unresolved.append({
                "asv_id":       asv_id,
                "taxon":        taxon,
                "confidence":   conf,
                "deepest_rank": deepest_rank,
                "deepest_val":  deepest_val,
                "reason":       "Unassigned",
            })
            continue

        depth_idx = RANK_ORDER.index(deepest_rank) if deepest_rank in RANK_ORDER else -1
        if depth_idx < target_idx or conf < min_confidence:
            unresolved.append({
                "asv_id":       asv_id,
                "taxon":        taxon,
                "confidence":   conf,
                "deepest_rank": deepest_rank,
                "deepest_val":  deepest_val,
                "reason":       (
                    f"Only classified to {deepest_rank} (need {target_rank})"
                    if depth_idx < target_idx
                    else f"Low confidence ({conf:.3f} < {min_confidence})"
                ),
            })

    return pd.DataFrame(unresolved)


# ---------------------------------------------------------------------------
# FASTA extraction
# ---------------------------------------------------------------------------

def extract_sequences(
    fasta_path: Path,
    asv_ids: Set[str],
) -> Dict[str, str]:
    """Extract sequences for given ASV IDs from a FASTA file."""
    seqs: Dict[str, str] = {}
    current_id: Optional[str] = None
    current_seq: List[str] = []

    with fasta_path.open() as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and current_id in asv_ids:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id and current_id in asv_ids:
        seqs[current_id] = "".join(current_seq)

    return seqs


def write_fasta(seqs: Dict[str, str], out_path: Path) -> None:
    """Write a dict of {id: sequence} to a FASTA file."""
    with out_path.open("w") as f:
        for seq_id, seq in seqs.items():
            f.write(f">{seq_id}\n{seq}\n")


# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------

def run_blast(
    query_fasta: Path,
    blast_db: Optional[str],
    out_path: Path,
    num_threads: int,
    max_target_seqs: int,
    remote: bool,
) -> bool:
    """
    Run BLASTn and write results to out_path.
    Returns True on success, False on failure.
    """
    cmd = [
        "blastn",
        "-query",   str(query_fasta),
        "-out",     str(out_path),
        "-max_target_seqs", str(max_target_seqs),
        "-outfmt",  "6 qseqid sseqid staxid pident length evalue bitscore",
    ]

    if remote:
        cmd += ["-remote"]
        log.info("Using remote BLAST (this may take several minutes per sequence)")
    else:
        if not blast_db:
            log.error("--blast-db is required unless --remote is specified")
            return False
        cmd += ["-db", blast_db, "-num_threads", str(num_threads)]

    log.info("Running BLAST: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log.error("BLAST failed:\n%s", result.stderr)
        return False

    if result.stderr:
        for line in result.stderr.strip().split("\n"):
            if line.strip():
                log.warning("BLAST stderr: %s", line)

    log.info("BLAST complete. Results: %s", out_path)
    return True


# ---------------------------------------------------------------------------
# Taxid resolution via NCBI E-utilities
# ---------------------------------------------------------------------------

def resolve_taxids(taxids: List[str]) -> Dict[str, str]:
    """
    Resolve a list of NCBI taxids to scientific names using E-utilities.
    Returns {taxid: scientific_name}.
    Falls back gracefully if efetch is unavailable or network fails.
    """
    if not taxids:
        return {}

    resolved: Dict[str, str] = {}
    unique_ids = list(set(t for t in taxids if t and t != "N/A"))

    if not unique_ids:
        return {}

    # Check efetch availability
    check = subprocess.run(["which", "efetch"], capture_output=True, text=True)
    if check.returncode != 0:
        log.warning("efetch not found on PATH — taxids will not be resolved to names. "
                    "Install NCBI E-utilities for full species name resolution.")
        return {tid: f"taxid:{tid}" for tid in unique_ids}

    log.info("Resolving %d unique taxids via NCBI E-utilities...", len(unique_ids))

    # Process in batches to avoid overloading NCBI
    batch_size = 50
    for i in range(0, len(unique_ids), batch_size):
        batch = unique_ids[i:i + batch_size]
        try:
            efetch_cmd = ["efetch", "-db", "taxonomy", "-id",
                          ",".join(batch), "-format", "xml"]
            xtract_cmd = ["xtract", "-pattern", "Taxon",
                          "-element", "TaxId,ScientificName,Rank"]

            efetch_proc = subprocess.run(efetch_cmd, capture_output=True, text=True,
                                         timeout=60)
            xtract_proc = subprocess.run(xtract_cmd, input=efetch_proc.stdout,
                                         capture_output=True, text=True, timeout=30)

            for line in xtract_proc.stdout.strip().split("\n"):
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    tid, name = parts[0], parts[1]
                    rank = parts[2] if len(parts) > 2 else ""
                    resolved[tid] = f"{name} ({rank})" if rank else name

            # Be polite to NCBI
            time.sleep(0.4)

        except subprocess.TimeoutExpired:
            log.warning("E-utilities timeout for batch starting at index %d", i)
        except Exception as e:
            log.warning("E-utilities error: %s", e)

    # Fill in any that failed
    for tid in unique_ids:
        if tid not in resolved:
            resolved[tid] = f"taxid:{tid} (unresolved)"

    log.info("Resolved %d / %d taxids", len(resolved), len(unique_ids))
    return resolved


# ---------------------------------------------------------------------------
# Result parsing and conflict detection
# ---------------------------------------------------------------------------

def parse_blast_results(
    blast_path: Path,
    taxid_names: Dict[str, str],
    conflict_taxa: List[str],
    min_pident: float,
) -> pd.DataFrame:
    """
    Parse BLAST tabular output and flag conflicts.

    Returns DataFrame with one row per ASV (top hit only) with columns:
    asv_id, top_hit_accession, taxid, species_name, pident, evalue,
    bitscore, conflict, conflict_reason
    """
    if not blast_path.exists() or blast_path.stat().st_size == 0:
        log.warning("BLAST output file is empty or missing: %s", blast_path)
        return pd.DataFrame()

    cols = ["asv_id", "accession", "taxid", "pident", "length",
            "evalue", "bitscore"]
    df = pd.read_csv(blast_path, sep="\t", header=None, names=cols,
                     dtype=str)
    df["pident"]   = pd.to_numeric(df["pident"],   errors="coerce")
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")

    # Keep only top hit per ASV (highest bitscore)
    df = df.sort_values("bitscore", ascending=False)
    df = df.groupby("asv_id").first().reset_index()

    # Resolve taxids to names
    df["species_name"] = df["taxid"].map(
        lambda t: taxid_names.get(str(t), f"taxid:{t}")
    )

    # Detect conflicts
    conflict_lower = [c.lower() for c in conflict_taxa]

    def check_conflict(row):
        name = str(row["species_name"]).lower()
        if row["pident"] < min_pident:
            return False, f"Low identity ({row['pident']:.1f}% < {min_pident}%)"
        for ct in conflict_lower:
            if ct in name:
                return True, f"BLAST hit '{row['species_name']}' conflicts with expected target taxa"
        return False, ""

    conflicts = df.apply(check_conflict, axis=1)
    df["conflict"]        = [c[0] for c in conflicts]
    df["conflict_reason"] = [c[1] for c in conflicts]

    return df


# ---------------------------------------------------------------------------
# Report writing
# ---------------------------------------------------------------------------

def write_report(
    unresolved_df: pd.DataFrame,
    blast_df: pd.DataFrame,
    marker: str,
    target_rank: str,
    min_confidence: float,
    outdir: Path,
) -> Tuple[Path, Path]:
    """
    Write human-readable QC report and machine-readable confirmed artefacts list.
    Returns (report_path, artefacts_path).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    report_path    = outdir / f"blast_qc_report_{marker}.txt"
    artefacts_path = outdir / f"confirmed_artefacts_{marker}.txt"

    # Merge unresolved info with BLAST results
    if not blast_df.empty:
        merged = unresolved_df.merge(
            blast_df[["asv_id", "accession", "taxid", "species_name",
                       "pident", "evalue", "bitscore", "conflict",
                       "conflict_reason"]],
            on="asv_id", how="left"
        )
    else:
        merged = unresolved_df.copy()
        for col in ["accession", "taxid", "species_name", "pident",
                    "evalue", "bitscore", "conflict", "conflict_reason"]:
            merged[col] = "N/A"

    # Confirmed artefacts = conflicting BLAST hits with high identity
    confirmed = merged[merged.get("conflict", False) == True]  # noqa: E712

    lines = [
        f"BLAST QC Report — {marker}",
        f"{'=' * 60}",
        f"Generated by: 09c_blast_qc_unclassified.py",
        f"Marker          : {marker}",
        f"Target rank     : {target_rank}",
        f"Min confidence  : {min_confidence}",
        f"",
        f"Summary",
        f"-------",
        f"  ASVs below target rank / low confidence : {len(unresolved_df)}",
        f"  ASVs with BLAST results                 : {len(blast_df)}",
        f"  Confirmed conflicts (recommend removal)  : {len(confirmed)}",
        f"",
        f"{'=' * 60}",
        f"CONFIRMED CONFLICTS — recommended for removal",
        f"{'=' * 60}",
    ]

    if confirmed.empty:
        lines.append("  None found.")
    else:
        for _, row in confirmed.iterrows():
            lines += [
                f"",
                f"  ASV ID     : {row['asv_id']}",
                f"  Classifier : {row['taxon']} (conf={row['confidence']:.3f})",
                f"  BLAST hit  : {row.get('species_name', 'N/A')} "
                f"({row.get('pident', 'N/A')}% identity)",
                f"  Accession  : {row.get('accession', 'N/A')}",
                f"  Reason     : {row.get('conflict_reason', 'N/A')}",
                f"  Action     : ADD to ARTEFACT_SPECIES in 11_clean_diet_table.py",
            ]

    lines += [
        f"",
        f"{'=' * 60}",
        f"ALL UNRESOLVED ASVs (for reference)",
        f"{'=' * 60}",
    ]

    for _, row in merged.iterrows():
        lines += [
            f"",
            f"  ASV ID     : {row['asv_id']}",
            f"  Classifier : {row['taxon']} (conf={row['confidence']:.3f})",
            f"  Deepest    : {row['deepest_rank']} = {row['deepest_val']}",
            f"  BLAST hit  : {row.get('species_name', 'N/A')} "
            f"({row.get('pident', 'N/A')}% identity)",
            f"  Conflict   : {'YES ⚠' if row.get('conflict') else 'no'}",
            f"  Reason     : {row.get('conflict_reason', '') or row['reason']}",
        ]

    lines += [
        f"",
        f"{'=' * 60}",
        f"NEXT STEPS",
        f"{'=' * 60}",
        f"",
        f"  1. Review confirmed conflicts above",
        f"  2. Add confirmed artefact ASV IDs to ARTEFACT_SPECIES in",
        f"     scripts/11_clean_diet_table.py",
        f"  3. Rerun 11_clean_diet_table.py to regenerate cleaned counts",
        f"  4. If ASVs are in QIIME2 feature table (diversity analyses),",
        f"     also filter with:",
        f"       qiime feature-table filter-features \\",
        f"         --i-table qiime2/{marker}/all/dada2/table.qza \\",
        f"         --m-metadata-file {artefacts_path} \\",
        f"         --p-exclude-ids \\",
        f"         --o-filtered-table qiime2/{marker}/all/dada2/table_filtered.qza",
        f"  5. Rerun core-metrics and regenerate all diversity figures",
        f"",
    ]

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log.info("Written report: %s", report_path)

    # Write machine-readable artefacts list (compatible with qiime feature-table)
    if not confirmed.empty:
        with artefacts_path.open("w") as f:
            f.write("#ASV-ID\tBlast-species\tpident\tReason\n")
            for _, row in confirmed.iterrows():
                f.write(
                    f"{row['asv_id']}\t"
                    f"{row.get('species_name', 'N/A')}\t"
                    f"{row.get('pident', 'N/A')}\t"
                    f"{row.get('conflict_reason', 'N/A')}\n"
                )
        log.info("Written artefacts list: %s", artefacts_path)
    else:
        log.info("No confirmed artefacts — artefacts file not written")
        artefacts_path = None

    return report_path, artefacts_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="09c_blast_qc_unclassified.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    req = p.add_argument_group("required")
    req.add_argument("--taxonomy",  required=True, type=Path,
                     help="Exported QIIME2 taxonomy TSV (Feature ID, Taxon, Confidence)")
    req.add_argument("--rep-seqs",  required=True, type=Path,
                     help="Exported QIIME2 rep-seqs FASTA (dna-sequences.fasta)")
    req.add_argument("--marker",    required=True,
                     help="Marker name for output labeling (e.g. MiFish, cytb)")
    req.add_argument("--outdir",    required=True, type=Path,
                     help="Output directory for FASTA, BLAST results, and report")

    blast = p.add_argument_group("BLAST options (one of --blast-db or --remote required)")
    blast.add_argument("--blast-db", default=None,
                       help="Path to local BLAST nt database (e.g. /share/databases/ncbi_nt/nt)")
    blast.add_argument("--remote",   action="store_true", default=False,
                       help="Use NCBI remote BLAST instead of local database (slower)")
    blast.add_argument("--num-threads", type=int, default=4,
                       help="BLAST threads (local only). Default: 4")
    blast.add_argument("--max-target-seqs", type=int, default=5,
                       help="Max BLAST hits per ASV. Default: 5")
    blast.add_argument("--min-pident", type=float, default=95.0,
                       help="Minimum BLAST identity (%%) to call a conflict. Default: 95.0")

    filt = p.add_argument_group("filter options")
    filt.add_argument("--target-rank",    default="order",
                      choices=RANK_ORDER,
                      help="ASVs not classified to this rank are BLASTed. Default: order")
    filt.add_argument("--min-confidence", type=float, default=0.80,
                      help="ASVs below this classifier confidence are BLASTed. Default: 0.80")
    filt.add_argument("--conflict-taxa",  nargs="+", default=None,
                      help="Species/genus names that signal a conflict with target taxa. "
                           "Default: loon (Gavia), bacteria, raptors. "
                           "ADAPT for your study organism.")
    filt.add_argument("--max-asvs",       type=int, default=200,
                      help="Maximum ASVs to BLAST (safety limit). Default: 200")

    util = p.add_argument_group("utility")
    util.add_argument("--skip-blast", action="store_true", default=False,
                      help="Skip BLAST step — just identify and report unresolved ASVs")
    util.add_argument("--dry-run",    action="store_true", default=False,
                      help="Print planned actions without running BLAST or writing files")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    # ── Validate inputs ───────────────────────────────────────────────────
    for path, name in [(args.taxonomy, "--taxonomy"),
                       (args.rep_seqs, "--rep_seqs")]:
        if not path.exists():
            log.error("%s not found: %s", name, path)
            return 2

    if not args.remote and not args.blast_db and not args.skip_blast:
        log.error("Provide --blast-db (local) or --remote, or use --skip-blast")
        return 2

    if args.blast_db and not Path(args.blast_db).parent.exists():
        log.error("BLAST database path not found: %s", args.blast_db)
        return 2

    conflict_taxa = args.conflict_taxa or CONFLICT_TAXA_DEFAULT

    # ── Load taxonomy ─────────────────────────────────────────────────────
    log.info("Loading taxonomy: %s", args.taxonomy)
    tax_df = pd.read_csv(args.taxonomy, sep="\t", dtype=str)
    tax_df.columns = [c.strip() for c in tax_df.columns]

    # Handle QIIME2 #q2:types header row
    tax_df = tax_df[~tax_df.iloc[:, 0].str.startswith("#", na=False)]

    # Normalise column names
    col_map = {}
    for col in tax_df.columns:
        cl = col.lower().strip()
        if cl in ("feature id", "featureid", "#featureid", "asv"):
            col_map[col] = "Feature ID"
        elif cl == "taxon":
            col_map[col] = "Taxon"
        elif cl == "confidence":
            col_map[col] = "Confidence"
    tax_df = tax_df.rename(columns=col_map)

    if "Feature ID" not in tax_df.columns:
        log.error("Could not find Feature ID column. Columns: %s",
                  list(tax_df.columns))
        return 1

    log.info("Loaded %d ASVs", len(tax_df))

    # ── Find unresolved ASVs ──────────────────────────────────────────────
    log.info("Finding ASVs not classified to '%s' or below confidence %.2f ...",
             args.target_rank, args.min_confidence)
    unresolved = find_unresolved_asvs(tax_df, args.target_rank,
                                      args.min_confidence)
    log.info("Found %d unresolved ASVs", len(unresolved))

    if unresolved.empty:
        log.info("No unresolved ASVs found. Pipeline looks clean for this marker.")
        return 0

    if len(unresolved) > args.max_asvs:
        log.warning(
            "%d unresolved ASVs exceeds --max-asvs=%d. "
            "BLASTing first %d only. Consider raising --min-confidence or "
            "--target-rank to narrow the set.",
            len(unresolved), args.max_asvs, args.max_asvs,
        )
        unresolved = unresolved.head(args.max_asvs)

    if args.dry_run:
        log.info("DRY RUN — would BLAST %d ASVs:", len(unresolved))
        for _, row in unresolved.iterrows():
            log.info("  %s  %s  (conf=%.3f)",
                     row["asv_id"], row["taxon"], row["confidence"])
        return 0

    # ── Extract sequences ─────────────────────────────────────────────────
    args.outdir.mkdir(parents=True, exist_ok=True)
    asv_ids = set(unresolved["asv_id"].tolist())

    log.info("Extracting %d sequences from %s ...", len(asv_ids), args.rep_seqs)
    seqs = extract_sequences(args.rep_seqs, asv_ids)
    log.info("Extracted %d sequences", len(seqs))

    missing = asv_ids - set(seqs.keys())
    if missing:
        log.warning("%d ASV IDs not found in rep-seqs FASTA: %s",
                    len(missing), list(missing)[:5])

    query_fasta = args.outdir / f"unresolved_{args.marker}.fasta"
    write_fasta(seqs, query_fasta)
    log.info("Written query FASTA: %s", query_fasta)

    # ── BLAST ─────────────────────────────────────────────────────────────
    blast_out = args.outdir / f"blast_results_{args.marker}.txt"
    blast_df  = pd.DataFrame()

    if args.skip_blast:
        log.info("--skip-blast specified — skipping BLAST step")
    else:
        success = run_blast(
            query_fasta    = query_fasta,
            blast_db       = args.blast_db,
            out_path       = blast_out,
            num_threads    = args.num_threads,
            max_target_seqs= args.max_target_seqs,
            remote         = args.remote,
        )

        if success and blast_out.exists():
            # Resolve taxids
            raw_df   = pd.read_csv(blast_out, sep="\t", header=None,
                                   names=["asv_id","accession","taxid","pident",
                                          "length","evalue","bitscore"],
                                   dtype=str)
            taxids   = raw_df["taxid"].dropna().unique().tolist()
            taxid_names = resolve_taxids(taxids)

            blast_df = parse_blast_results(
                blast_path   = blast_out,
                taxid_names  = taxid_names,
                conflict_taxa= conflict_taxa,
                min_pident   = args.min_pident,
            )
            log.info("Parsed BLAST results: %d ASVs with hits, %d conflicts",
                     len(blast_df),
                     blast_df["conflict"].sum() if not blast_df.empty else 0)
        else:
            log.warning("BLAST failed or produced no output — report will lack BLAST data")

    # ── Write report ──────────────────────────────────────────────────────
    report_path, artefacts_path = write_report(
        unresolved_df  = unresolved,
        blast_df       = blast_df,
        marker         = args.marker,
        target_rank    = args.target_rank,
        min_confidence = args.min_confidence,
        outdir         = args.outdir,
    )

    log.info("")
    log.info("=== Done ===")
    log.info("Report       : %s", report_path)
    if artefacts_path:
        log.info("Artefacts    : %s", artefacts_path)
        log.info("")
        log.info("Next steps:")
        log.info("  1. Review %s", report_path)
        log.info("  2. Add confirmed ASV IDs to ARTEFACT_SPECIES in "
                 "scripts/11_clean_diet_table.py")
        log.info("  3. Filter QIIME2 feature table using %s", artefacts_path)
        log.info("  4. Rerun core-metrics and regenerate diversity figures")
    else:
        log.info("No confirmed artefacts found — no action needed")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
