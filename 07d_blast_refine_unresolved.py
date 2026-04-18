#!/usr/bin/env python3
"""
07d_blast_refine_unresolved.py
===============================
BLAST-based refinement of ASVs that classifiers could not resolve to a
target taxonomic rank (typically genus).

Complementary to 07c_blast_qc_unclassified.py:
  - 07c_ flags ASVs where classifier output CONFLICTS with expected taxa
    (e.g. host DNA miscalled as fish in MiFish/cytb).
  - 07d_ refines ASVs where the classifier STOPPED TOO HIGH (e.g.
    "f__Enterobacteriaceae" with no genus) and BLAST against NCBI nt can
    push them down using more comprehensive references.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PIPELINE POSITION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  dada2/rep-seqs.qza ──┐
  taxonomy/taxonomy.qza ┼──► 07d_blast_refine_unresolved.py
                        │         ↓
                  blast_summary_{marker}.tsv         ← per-ASV genus calls
                  refined_taxonomy_{marker}.tsv      ← patched taxonomy (w/ --apply)
                  blast_refine_report_{marker}.txt   ← human review
                        ↓
                  07_taxonomy_table.py (re-run with refined TSV)
                        ↓
                  10_plot_taxonomy.py

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WHY THIS EXISTS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SILVA and similar classifiers stop at family or order level for common
gut/environmental bacteria when V4 16S resolution hits species limits.
Examples:

  - "f__Enterobacteriaceae" with no genus: V4 cannot fully separate
    E. coli / Shigella / Salmonella / Klebsiella, but CAN usually pick
    Escherichia vs Citrobacter vs Klebsiella as the dominant call.
  - "o__Lactobacillales" with no family: V4 typically resolves
    Carnobacterium, Vagococcus, Enterococcus against NCBI nt.

BLAST against nt uses a far larger reference set than SILVA's
representatives. For well-represented genera this yields clean genus
calls. For poorly-represented taxa (common in wildlife gut microbiomes),
BLAST hits are dominated by "uncultured bacterium" deposits and no
refinement is possible — these are correctly flagged as unresolved
rather than forcing a call.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WHAT THIS SCRIPT DOES NOT DO
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  - Does NOT rescue ASVs whose top hits are "uncultured bacterium" at
    high identity — these are genuinely unresolved in public databases.
  - Does NOT separate genera that share identical V4 (e.g. classic
    Escherichia/Shigella/Salmonella ambiguity). These remain ambiguous.
  - Does NOT replace manual review. Review blast_summary_{marker}.tsv
    and blast_refine_report_{marker}.txt before using --apply output
    downstream.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
USAGE EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # 16S — refine all ASVs SILVA stopped above genus, write patched TSV
  python scripts/07d_blast_refine_unresolved.py \\
      --taxonomy  qiime2/16S/all/taxonomy-exported/taxonomy.tsv \\
      --rep-seqs  qiime2/16S/all/rep-seqs-exported/dna-sequences.fasta \\
      --marker    16S \\
      --blast-db  /home/share/databases/ncbi_nt/nt \\
      --outdir    qiime2/16S/all/blast_refine/ \\
      --apply

  # Target specific uncl. labels only
  python scripts/07d_blast_refine_unresolved.py \\
      --taxonomy     qiime2/16S/all/taxonomy-exported/taxonomy.tsv \\
      --rep-seqs     qiime2/16S/all/rep-seqs-exported/dna-sequences.fasta \\
      --marker       16S \\
      --blast-db     /home/share/databases/ncbi_nt/nt \\
      --target-labels f__Enterobacteriaceae o__Lactobacillales \\
      --outdir       qiime2/16S/all/blast_refine/

  # Dry run — see candidate breakdown before running BLAST
  python scripts/07d_blast_refine_unresolved.py \\
      --taxonomy  qiime2/16S/all/taxonomy-exported/taxonomy.tsv \\
      --rep-seqs  qiime2/16S/all/rep-seqs-exported/dna-sequences.fasta \\
      --marker    16S \\
      --blast-db  /home/share/databases/ncbi_nt/nt \\
      --outdir    qiime2/16S/all/blast_refine/ \\
      --dry-run

  # Re-score with different dominance threshold (BLAST already done)
  python scripts/07d_blast_refine_unresolved.py \\
      --taxonomy        qiime2/16S/all/taxonomy-exported/taxonomy.tsv \\
      --rep-seqs        qiime2/16S/all/rep-seqs-exported/dna-sequences.fasta \\
      --marker          16S \\
      --existing-blast  qiime2/16S/all/blast_refine/blast_results_16S.tsv \\
      --dominance-threshold 0.60 \\
      --outdir          qiime2/16S/all/blast_refine/

Dependencies:
  - BLAST+ (blastn on PATH)
  - Local NCBI nt database with taxdb files alongside it
    (taxdb.btd, taxdb.bti, taxonomy4blast.sqlite3) so sscinames resolve.
    Set BLASTDB=/path/to/nt_dir or rely on the DB path.
  - pip install pandas
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RANK_PREFIXES: Dict[str, str] = {
    "kingdom": "d__",
    "phylum":  "p__",
    "class":   "c__",
    "order":   "o__",
    "family":  "f__",
    "genus":   "g__",
    "species": "s__",
}

RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

# Hits whose first token looks like a genus but is actually a higher rank
# or environmental descriptor. Keep expanding as new cases surface.
NON_GENUS_LABELS = {
    "Bacteria", "Archaea", "Eukaryota",
    "Enterobacteriaceae", "Lactobacillales", "Bacillales", "Clostridiales",
    "Actinobacteria", "Firmicutes", "Proteobacteria", "Bacteroidota",
    "Gamma", "Alpha", "Beta", "Delta", "Epsilon",
    "uncultured", "unclassified", "unknown", "unidentified",
}

NON_GENUS_PREFIXES = ("uncultured", "bacterium", "unidentified", "unknown")


# ---------------------------------------------------------------------------
# Taxonomy parsing and filtering
# ---------------------------------------------------------------------------

def load_taxonomy(path: Path) -> pd.DataFrame:
    """Load QIIME2 taxonomy TSV (Feature ID, Taxon, Confidence)."""
    df = pd.read_csv(path, sep="\t")
    required = {"Feature ID", "Taxon", "Confidence"}
    if not required.issubset(df.columns):
        log.error(
            "Taxonomy file missing expected columns. "
            "Found: %s, required: %s",
            list(df.columns), sorted(required),
        )
        sys.exit(2)
    df["Confidence"] = pd.to_numeric(df["Confidence"], errors="coerce")
    log.info("Loaded %d taxonomy entries from %s", len(df), path)
    return df


def taxon_reaches_rank(taxon: str, target_rank: str) -> bool:
    """
    True if the QIIME2 taxonomy string has a non-empty field at target_rank.
    Handles 'Unassigned' and missing fields.
    """
    if not isinstance(taxon, str) or taxon.strip() == "Unassigned":
        return False
    target_prefix = RANK_PREFIXES[target_rank]
    # "; g__Something" — at least one alphanumeric char after the prefix
    pattern = rf"(^|;\s*){re.escape(target_prefix)}[A-Za-z0-9]"
    return bool(re.search(pattern, taxon))


def taxon_matches_label(taxon: str, label: str) -> bool:
    """
    True if `label` is the deepest classified rank in the taxon string
    (i.e., the rank after `label` is empty or absent).

    Examples with label='f__Enterobacteriaceae':
      'd__Bacteria; ...; f__Enterobacteriaceae'                 → True
      'd__Bacteria; ...; f__Enterobacteriaceae; g__Escherichia' → False
    """
    if not isinstance(taxon, str):
        return False
    if label not in taxon:
        return False
    idx = taxon.index(label)
    tail = taxon[idx + len(label):].strip()
    if not tail:
        return True
    if tail.startswith(";"):
        rest = tail.lstrip(";").strip()
        if not rest:
            return True
        # Next rank present but empty (e.g. "g__" alone)
        for prefix in ("g__", "s__", "f__", "o__", "c__"):
            first_field = rest.split(";")[0].strip()
            if first_field == prefix:  # empty rank marker
                return True
    return False


def select_unresolved(
    tax_df: pd.DataFrame,
    target_rank: str,
    target_labels: Optional[List[str]] = None,
    min_confidence: float = 0.0,
    require_parent_rank: Optional[str] = None,
) -> pd.DataFrame:
    """
    Select ASVs that are candidates for BLAST refinement.

    An ASV is a candidate if it:
      1. Does NOT reach `target_rank`, AND
      2. DOES reach `require_parent_rank` (defaults to rank-1 above target,
         excluding 'Unassigned' ASVs that likely aren't bacteria), AND
      3. (If target_labels given) its taxon string matches one of them
         at its deepest classified rank, AND
      4. Has classifier confidence >= min_confidence.
    """
    if require_parent_rank is None:
        idx = RANK_ORDER.index(target_rank)
        require_parent_rank = RANK_ORDER[max(idx - 1, 0)]

    # If user provided explicit target_labels, relax the parent-rank
    # requirement to the SHALLOWEST rank among the labels — otherwise an
    # order-level label (e.g. 'o__Lactobacillales') would be excluded by
    # a family-level parent requirement.
    effective_parent = require_parent_rank
    if target_labels:
        label_ranks = []
        for lb in target_labels:
            for rname, rprefix in RANK_PREFIXES.items():
                if lb.startswith(rprefix):
                    label_ranks.append(rname)
                    break
        if label_ranks:
            shallowest_idx = min(RANK_ORDER.index(r) for r in label_ranks)
            if shallowest_idx < RANK_ORDER.index(effective_parent):
                effective_parent = RANK_ORDER[shallowest_idx]

    mask = tax_df["Taxon"].apply(
        lambda t: (
            not taxon_reaches_rank(t, target_rank)
            and taxon_reaches_rank(t, effective_parent)
        )
    )
    if min_confidence > 0:
        mask &= tax_df["Confidence"] >= min_confidence
    if target_labels:
        label_mask = tax_df["Taxon"].apply(
            lambda t: any(taxon_matches_label(t, lb) for lb in target_labels)
        )
        mask &= label_mask

    selected = tax_df[mask].copy()
    log.info(
        "Selected %d candidate ASVs (target_rank=%s, require_parent=%s, "
        "min_confidence=%.2f, target_labels=%s)",
        len(selected), target_rank, effective_parent,
        min_confidence, target_labels,
    )
    return selected


# ---------------------------------------------------------------------------
# FASTA extraction
# ---------------------------------------------------------------------------

def extract_target_sequences(
    ids: List[str],
    fasta_path: Path,
    out_path: Path,
) -> int:
    """Extract sequences with matching IDs from a FASTA file."""
    id_set = set(ids)
    found = set()
    printing = False
    with fasta_path.open() as fin, out_path.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                name = line[1:].split()[0].strip()
                printing = name in id_set
                if printing:
                    found.add(name)
                    fout.write(line)
            elif printing:
                fout.write(line)

    missing = id_set - found
    if missing:
        log.warning(
            "%d ASV(s) in taxonomy not found in rep-seqs FASTA — "
            "taxonomy and rep-seqs may be out of sync. First few missing: %s",
            len(missing), list(missing)[:3],
        )
    log.info("Wrote %d sequences to %s", len(found), out_path)
    return len(found)


# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------

def run_blast(
    query_fasta: Path,
    blast_db: str,
    out_path: Path,
    max_target_seqs: int = 50,
    num_threads: int = 4,
    min_pident: Optional[float] = None,
) -> None:
    """Run blastn with a standardized output format."""
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", blast_db,
        "-outfmt",
        "6 qseqid sacc sscinames staxids pident length evalue bitscore stitle",
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(num_threads),
        "-out", str(out_path),
    ]
    if min_pident is not None:
        cmd += ["-perc_identity", str(min_pident)]

    log.info("Running BLAST against %s with %d threads, max_target_seqs=%d",
             blast_db, num_threads, max_target_seqs)
    log.info("  %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        log.error("blastn not found on PATH. Activate BLAST+ env / module.")
        sys.exit(127)
    except subprocess.CalledProcessError as e:
        log.error("blastn exited with code %d", e.returncode)
        sys.exit(e.returncode)

    n_hits = sum(1 for _ in out_path.open())
    log.info("BLAST complete. %d hits written to %s", n_hits, out_path)


# ---------------------------------------------------------------------------
# Result scoring
# ---------------------------------------------------------------------------

def extract_genus(sciname: str) -> Optional[str]:
    """Pull the genus from an NCBI sscinames string. None if not a named genus."""
    if not isinstance(sciname, str):
        return None
    low = sciname.lower()
    for prefix in NON_GENUS_PREFIXES:
        if low.startswith(prefix):
            return None
    if "uncultured" in low:
        return None
    parts = sciname.split()
    if not parts:
        return None
    first = parts[0]
    if first in NON_GENUS_LABELS:
        return None
    if not first[0].isupper():
        return None
    return first


def score_blast_results(
    blast_tsv: Path,
    dominance_threshold: float = 0.70,
) -> pd.DataFrame:
    """Score BLAST results into per-ASV genus calls with dominance promotion."""
    cols = ["asv", "acc", "sciname", "taxid", "pident",
            "length", "evalue", "bitscore", "stitle"]
    df = pd.read_csv(blast_tsv, sep="\t", header=None, names=cols)
    df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
    df["genus"] = df["sciname"].apply(extract_genus)

    rows = []
    for asv, g in df.groupby("asv"):
        top_pident = g["pident"].max()
        top_hits = g[g["pident"] == top_pident]
        top_genera = [x for x in top_hits["genus"] if x]
        genus_counts = Counter(top_genera)
        total_named = sum(genus_counts.values())

        if not top_genera:
            call = "UNRESOLVED (only uncultured hits)"
            category, called = "unresolved", None
        elif len(genus_counts) == 1:
            called = list(genus_counts.keys())[0]
            call, category = f"GENUS: {called}", "resolved"
        else:
            top_genus, top_count = genus_counts.most_common(1)[0]
            frac = top_count / total_named
            if frac >= dominance_threshold:
                minor = ", ".join(f"{x}({n})" for x, n in
                                  genus_counts.most_common()[1:])
                call = (f"GENUS (dominant {top_count}/{total_named}): "
                        f"{top_genus}  [minor: {minor}]")
                category, called = "dominant", top_genus
            else:
                top_list = ", ".join(f"{x}({n})" for x, n in
                                     genus_counts.most_common())
                call = f"AMBIGUOUS: {top_list}"
                category, called = "ambiguous", None

        rows.append({
            "asv": asv,
            "top_pident": top_pident,
            "n_top_hits": len(top_hits),
            "n_genera_top": len(genus_counts),
            "called_genus": called,
            "category": category,
            "call": call,
        })

    summary = pd.DataFrame(rows)
    cat_order = {"resolved": 0, "dominant": 1, "ambiguous": 2, "unresolved": 3}
    summary["_sort"] = summary["category"].map(cat_order)
    summary = summary.sort_values(["_sort", "top_pident"],
                                  ascending=[True, False])
    summary = summary.drop(columns=["_sort"]).reset_index(drop=True)
    return summary


# ---------------------------------------------------------------------------
# Taxonomy patching
# ---------------------------------------------------------------------------

# Genus → family lookup for common bacterial refinements.
# Used to fill in the family rank when appending a BLAST-derived genus to a
# taxonomy string that stopped at order (e.g. "o__Lactobacillales"). Keeping
# the SILVA-style full path (d__;p__;c__;o__;f__;g__) ensures the refined
# rows collapse onto the same L6 labels as SILVA-native rows in 07.
# Expand as new genera appear in BLAST output.
GENUS_TO_FAMILY: Dict[str, str] = {
    # Enterobacteriaceae (Enterobacterales)
    "Escherichia":   "Enterobacteriaceae",
    "Shigella":      "Enterobacteriaceae",
    "Salmonella":    "Enterobacteriaceae",
    "Citrobacter":   "Enterobacteriaceae",
    "Klebsiella":    "Enterobacteriaceae",
    "Enterobacter":  "Enterobacteriaceae",
    "Atlantibacter": "Enterobacteriaceae",
    "Kluyvera":      "Enterobacteriaceae",
    "Raoultella":    "Enterobacteriaceae",
    "Cedecea":       "Enterobacteriaceae",
    "Edwardsiella":  "Hafniaceae",
    "Hafnia":        "Hafniaceae",
    # Lactobacillales families
    "Carnobacterium":   "Carnobacteriaceae",
    "Vagococcus":       "Enterococcaceae",
    "Enterococcus":     "Enterococcaceae",
    "Catellicoccus":    "Carnobacteriaceae",
    "Lactobacillus":    "Lactobacillaceae",
    "Lactococcus":      "Streptococcaceae",
    "Streptococcus":    "Streptococcaceae",
    "Leuconostoc":      "Leuconostocaceae",
}


def build_refined_taxon(original: str, called_genus: str) -> str:
    """
    Append g__<genus> (and f__<family> if missing) to the taxonomy string.

    When the refined genus has a known family in GENUS_TO_FAMILY, the family
    rank is filled in too, so downstream QIIME2 taxa-collapse treats the
    refined row the same as SILVA-native rows with that genus (preventing
    split 'Carnobacterium' / 'uncl. Carnobacterium' rows in L6 output).
    """
    # Overwrite any existing populated g__ (shouldn't happen for targeted
    # unresolved ASVs, but handle defensively).
    if re.search(r"g__[A-Za-z0-9]", original):
        return re.sub(r"g__[^;]*", f"g__{called_genus}", original)

    trimmed = original.rstrip().rstrip(";").rstrip()
    family = GENUS_TO_FAMILY.get(called_genus)

    # If family is known and not already present in the string, insert it.
    if family and not re.search(rf"f__{re.escape(family)}(\b|;|$)", trimmed):
        if re.search(r"f__[A-Za-z0-9]", trimmed):
            # Some other family is already populated — don't overwrite, just
            # append the genus (preserves weird upstream cases where an
            # unusual family was assigned).
            return f"{trimmed}; g__{called_genus}"
        # No family yet — insert f__<family> before appending g__<genus>.
        return f"{trimmed}; f__{family}; g__{called_genus}"

    return f"{trimmed}; g__{called_genus}"


def apply_refinements(
    tax_df: pd.DataFrame,
    summary: pd.DataFrame,
    out_path: Path,
) -> Tuple[int, int]:
    """Write a new taxonomy TSV with genus calls applied."""
    calls = summary[summary["called_genus"].notna()].set_index("asv")["called_genus"]
    out_df = tax_df.copy()

    n_refined = 0
    for asv, genus in calls.items():
        mask = out_df["Feature ID"] == asv
        if not mask.any():
            log.warning("ASV %s not found in taxonomy — skipping.", asv[:20])
            continue
        original = out_df.loc[mask, "Taxon"].values[0]
        out_df.loc[mask, "Taxon"] = build_refined_taxon(original, genus)
        n_refined += 1

    out_df.to_csv(out_path, sep="\t", index=False)
    log.info("Refined taxonomy: %d / %d ASVs updated → %s",
             n_refined, len(out_df), out_path)
    return n_refined, len(out_df)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report(
    summary: pd.DataFrame,
    marker: str,
    out_path: Path,
    target_rank: str,
    target_labels: Optional[List[str]],
    blast_db: str,
    max_target_seqs: int,
    dominance_threshold: float,
    refined_out: Optional[Path] = None,
) -> None:
    """Human-readable report."""
    L = []
    L.append("=" * 72)
    L.append(f"BLAST refinement report — {marker}")
    L.append("=" * 72)
    L.append("")
    L.append(f"Target rank         : {target_rank}")
    L.append(f"Target labels       : {target_labels or '(all unresolved at target rank)'}")
    L.append(f"BLAST database      : {blast_db}")
    L.append(f"Max target seqs     : {max_target_seqs}")
    L.append(f"Dominance threshold : {dominance_threshold:.2f}")
    L.append("")

    L.append("-" * 72)
    L.append("Per-ASV outcome")
    L.append("-" * 72)
    total = len(summary)
    for cat in ["resolved", "dominant", "ambiguous", "unresolved"]:
        n = (summary["category"] == cat).sum()
        pct = 100 * n / total if total else 0.0
        L.append(f"  {cat:12s}: {n:4d}  ({pct:5.1f}%)")
    called = summary["called_genus"].notna().sum()
    L.append("")
    if total:
        L.append(f"  Total with genus call: {called} / {total} "
                 f"({100*called/total:.1f}%)")
    else:
        L.append("  (no candidates)")
    L.append("")

    if called > 0:
        L.append("-" * 72)
        L.append("Genus calls (resolved + dominant combined)")
        L.append("-" * 72)
        counts = summary[summary["called_genus"].notna()]["called_genus"].value_counts()
        for genus, n in counts.items():
            L.append(f"  {n:4d}  {genus}")
        L.append("")

    ambig = summary[summary["category"] == "ambiguous"]
    if len(ambig) > 0:
        L.append("-" * 72)
        L.append("Ambiguous cases (manual review)")
        L.append("-" * 72)
        L.append("These ASVs have multiple genera tied at top identity and no")
        L.append("single genus reaches the dominance threshold. Typical causes:")
        L.append("V4 cannot separate closely-related genera (e.g. Escherichia/")
        L.append("Shigella/Salmonella). Left unchanged in --apply output.")
        L.append("")
        for _, row in ambig.iterrows():
            L.append(f"  {row['asv'][:20]}  {row['call'][:120]}")
        L.append("")

    L.append("-" * 72)
    L.append("Downstream usage")
    L.append("-" * 72)
    L.append("")
    if refined_out:
        L.append(f"Refined taxonomy TSV written: {refined_out}")
        L.append("")
        L.append("To use in 07_taxonomy_table.py:")
        L.append("")
        L.append("  # Option A — import back into QIIME2 as a .qza:")
        L.append("  qiime tools import \\")
        L.append(f"    --input-path {refined_out} \\")
        L.append("    --type 'FeatureData[Taxonomy]' \\")
        L.append("    --input-format HeaderlessTSVTaxonomyFormat \\")
        L.append(f"    --output-path {refined_out.with_suffix('.qza')}")
        L.append("")
        L.append("  # Option B — point 07_taxonomy_table.py at the TSV directly")
        L.append("  # (if it supports TSV input).")
        L.append("")
    else:
        L.append("Run with --apply to write a patched taxonomy TSV.")
        L.append("")

    L.append("Then:")
    L.append("  1. Re-run 07_taxonomy_table.py with the refined taxonomy.")
    L.append("  2. Re-run 10_plot_taxonomy.py — uncl. bars should split into")
    L.append("     named genera where BLAST resolved them.")
    L.append("")

    out_path.write_text("\n".join(L))
    log.info("Report written: %s", out_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="07d_blast_refine_unresolved.py",
        description=(
            "BLAST-based refinement of ASVs poorly-classified at target rank. "
            "Runs blastn against a local NCBI nt database, scores top hits with "
            "a dominance-threshold rule, and optionally writes a patched "
            "taxonomy TSV for downstream re-plotting."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    req = p.add_argument_group("required arguments")
    req.add_argument("--taxonomy", required=True, type=Path,
                     help="Exported QIIME2 taxonomy TSV (Feature ID, Taxon, Confidence).")
    req.add_argument("--rep-seqs", required=True, type=Path, dest="rep_seqs",
                     help="Exported QIIME2 rep-seqs FASTA (dna-sequences.fasta).")
    req.add_argument("--marker", required=True,
                     help="Marker name for output labeling (e.g. 16S, cytb).")
    req.add_argument("--outdir", required=True, type=Path,
                     help="Output directory for FASTA, BLAST results, summary, report.")

    bl = p.add_argument_group("BLAST options (one of --blast-db or --existing-blast required)")
    bl.add_argument("--blast-db", default=None, dest="blast_db",
                    help="Path to local BLAST nt database "
                         "(e.g. /home/share/databases/ncbi_nt/nt).")
    bl.add_argument("--existing-blast", default=None, type=Path, dest="existing_blast",
                    help="Use a pre-existing BLAST results TSV (skip BLAST step).")
    bl.add_argument("--num-threads", type=int, default=None, dest="num_threads",
                    help="BLAST threads. Default: os.cpu_count().")
    bl.add_argument("--max-target-seqs", type=int, default=50, dest="max_target_seqs",
                    help="Max BLAST hits per ASV. Default: 50. "
                         "(10 is too few for common taxa; 50 balances coverage "
                         "and runtime.)")
    bl.add_argument("--min-pident", type=float, default=None, dest="min_pident",
                    help="Minimum BLAST percent identity. Default: unset.")

    fi = p.add_argument_group("filter options")
    fi.add_argument("--target-rank", choices=list(RANK_PREFIXES.keys()),
                    default="genus", dest="target_rank",
                    help="ASVs not classified to this rank are refined. Default: genus.")
    fi.add_argument("--target-labels", nargs="*", default=None, dest="target_labels",
                    help="Optional: restrict to specific taxonomy labels "
                         "(e.g. 'f__Enterobacteriaceae o__Lactobacillales'). "
                         "Applied as 'deepest classified rank equals this label'.")
    fi.add_argument("--min-confidence", type=float, default=0.0, dest="min_confidence",
                    help="Min classifier confidence for inclusion. Default: 0.0. "
                         "Set to ~0.7 to skip noisy low-confidence classifications.")
    fi.add_argument("--max-asvs", type=int, default=500, dest="max_asvs",
                    help="Safety cap on candidate ASVs. Default: 500.")

    sc = p.add_argument_group("scoring options")
    sc.add_argument("--dominance-threshold", type=float, default=0.70,
                    dest="dominance_threshold",
                    help="Fraction of top-identity hits a single genus must hold "
                         "to be promoted from 'ambiguous' to 'dominant'. Default: 0.70.")

    ou = p.add_argument_group("output options")
    ou.add_argument("--apply", action="store_true",
                    help="Write a patched taxonomy TSV with refined genus calls. "
                         "Unresolved/ambiguous ASVs are preserved unchanged.")
    ou.add_argument("--dry-run", action="store_true", dest="dry_run",
                    help="Report candidate counts without running BLAST.")

    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    # Validate
    if not args.taxonomy.exists():
        log.error("Taxonomy file not found: %s", args.taxonomy); return 2
    if not args.rep_seqs.exists():
        log.error("Rep-seqs FASTA not found: %s", args.rep_seqs); return 2
    if not args.blast_db and not args.existing_blast:
        log.error("Must provide either --blast-db or --existing-blast."); return 2
    if args.existing_blast and not args.existing_blast.exists():
        log.error("Existing BLAST TSV not found: %s", args.existing_blast); return 2

    args.outdir.mkdir(parents=True, exist_ok=True)
    threads = args.num_threads or (os.cpu_count() or 4)

    # Load + filter
    tax_df = load_taxonomy(args.taxonomy)
    candidates = select_unresolved(
        tax_df,
        target_rank=args.target_rank,
        target_labels=args.target_labels,
        min_confidence=args.min_confidence,
    )

    if len(candidates) == 0:
        log.warning("No candidate ASVs match the filters. Nothing to do.")
        return 0
    if len(candidates) > args.max_asvs:
        log.error(
            "Candidate count %d exceeds --max-asvs %d. "
            "Tighten filters or raise the cap.",
            len(candidates), args.max_asvs,
        )
        return 1

    # Candidate breakdown by deepest label
    log.info("Candidate breakdown (top 20 deepest-label groups):")
    label_counts: Counter = Counter()
    for t in candidates["Taxon"]:
        fields = [f.strip() for f in str(t).split(";") if f.strip()]
        label_counts[fields[-1] if fields else "(empty)"] += 1
    for lab, n in label_counts.most_common(20):
        log.info("  %4d  %s", n, lab[:100])

    if args.dry_run:
        log.info("Dry run — would BLAST %d ASVs. Exiting.", len(candidates))
        return 0

    # BLAST (or load existing)
    blast_tsv = args.outdir / f"blast_results_{args.marker}.tsv"
    if args.existing_blast:
        blast_tsv = args.existing_blast
        log.info("Using existing BLAST results: %s", blast_tsv)
    else:
        target_fasta = args.outdir / f"blast_targets_{args.marker}.fasta"
        extract_target_sequences(
            candidates["Feature ID"].tolist(),
            args.rep_seqs, target_fasta,
        )
        run_blast(
            query_fasta=target_fasta, blast_db=args.blast_db, out_path=blast_tsv,
            max_target_seqs=args.max_target_seqs, num_threads=threads,
            min_pident=args.min_pident,
        )

    # Score
    summary = score_blast_results(blast_tsv, args.dominance_threshold)
    summary_path = args.outdir / f"blast_summary_{args.marker}.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    log.info("Summary written: %s", summary_path)

    # Apply
    refined_out: Optional[Path] = None
    if args.apply:
        refined_out = args.outdir / f"refined_taxonomy_{args.marker}.tsv"
        apply_refinements(tax_df, summary, refined_out)

    # Report
    report_path = args.outdir / f"blast_refine_report_{args.marker}.txt"
    write_report(
        summary=summary, marker=args.marker, out_path=report_path,
        target_rank=args.target_rank, target_labels=args.target_labels,
        blast_db=args.blast_db or f"(existing: {args.existing_blast})",
        max_target_seqs=args.max_target_seqs,
        dominance_threshold=args.dominance_threshold,
        refined_out=refined_out,
    )

    # Stdout summary
    log.info("=" * 60)
    for cat in ["resolved", "dominant", "ambiguous", "unresolved"]:
        n = (summary["category"] == cat).sum()
        log.info("  %-12s: %4d", cat, n)
    called = summary["called_genus"].notna().sum()
    log.info("  %-12s: %4d", "total called", called)
    log.info("=" * 60)
    if not args.apply:
        log.info("Run with --apply to write a patched taxonomy TSV.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
