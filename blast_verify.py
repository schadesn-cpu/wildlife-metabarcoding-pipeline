#!/usr/bin/env python3
"""
blast_verify.py
===============
Step: blast (optional `analyses.blast` stage; targeted verification after taxonomy)

Purpose:
    BLAST-verify specific suspect or ambiguous taxa from the taxonomy count table
    against a local NCBI nt (or remote) database, then compare the BLAST top hit
    to the classifier's assignment and flag agreements vs disagreements. Use it to
    confirm ecologically implausible IDs, pin down unresolved taxa, or check
    candidate artefacts before trusting them. Complements blast_refine (fills in
    under-resolved calls) and blast_qc (flags conflicts): this one answers
    "is THIS specific call real?" for taxa you name or that clear a read threshold.

Target selection (at least one; from config or CLI):
    --taxa       taxon-name substrings to verify (e.g. Brevoortia Lutjanidae)
    --min-reads  verify every taxon with >= this many total reads

Inputs:
    --taxonomy  FeatureData[Taxonomy] .qza   (exported internally; from --marker)
    --rep-seqs  FeatureData[Sequence] .qza   (exported internally; from --marker)
    --counts    taxonomy count TSV (from taxonomy_table.py; derived from --marker)
    --db        path to a local BLAST nt DB, no extension (from analyses.blast.db)
    pipeline_config.yml  analyses.blast.verify (min_reads, identity, hits, evalue)

Outputs (in --outdir, default results/<marker>/all/blast/):
    blast_verify_<marker>.tsv             per-target BLAST top hits vs classifier
    blast_verify_report_<marker>.txt      human-review report
    logs/run_manifest.jsonl               run appended on completion

Usage:
    python blast_verify.py --marker MiFish                    # min_reads from config
    python blast_verify.py --marker cytb --taxa Brevoortia    # verify named taxa
    python blast_verify.py --marker MiFish --remote           # NCBI remote BLAST

Requirements:
    QIIME 2 (exports the .qza inputs) and BLAST+ (`blastn`). A local nt DB is
    strongly preferred over --remote.
"""
from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# --- make config_loader and the utils package importable regardless of cwd ---
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from config_loader import load_config, get_paths  # noqa: E402
from utils import checkpoint, provenance, validate  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Geographic plausibility flags
# These taxa are ecologically implausible in the northeastern US and should
# be flagged as likely artefacts regardless of BLAST result.
# ADAPT: update for your study region.
# ---------------------------------------------------------------------------
IMPLAUSIBLE_TAXA = [
    # Tropical/Indo-Pacific fish families — not in Gulf of Maine
    "Lutjanidae",      # tropical snappers
    "Siganidae",       # rabbitfish — Indo-Pacific
    "Apogonidae",      # cardinalfish — tropical
    "Carangidae",      # jacks/pompano — subtropical
    "Haemulidae",      # grunts — tropical Atlantic
    "Serranidae",      # sea bass — mostly tropical
    # European/Southern Hemisphere fish
    "Sprattus sprattus",  # European sprat — Eastern Atlantic only
    "Eleginopidae",       # Patagonian toothfish relatives
    # Terrestrial mammals (lab contamination)
    "Mustelidae",      # weasels/fishers/otters
    "Bovidae",         # cattle
    "Canidae",         # dogs
    "Felidae",         # cats
    "Hominidae",       # human
    "Suidae",          # pigs
    "Equidae",         # horses
]


# ---------------------------------------------------------------------------
# QIIME2 export helpers
# ---------------------------------------------------------------------------

def export_qza(qza_path: Path, out_dir: Path) -> Path:
    """Export a QIIME2 artifact to a directory and return the output path."""
    out_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        ["qiime", "tools", "export",
         "--input-path", str(qza_path),
         "--output-path", str(out_dir)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to export {qza_path}:\n{result.stderr}"
        )
    return out_dir


def load_taxonomy_tsv(taxonomy_dir: Path) -> pd.DataFrame:
    """Load exported taxonomy TSV into a DataFrame indexed by feature ID."""
    tsv = taxonomy_dir / "taxonomy.tsv"
    if not tsv.exists():
        raise FileNotFoundError(f"taxonomy.tsv not found in {taxonomy_dir}")
    df = pd.read_csv(tsv, sep="\t", index_col=0)
    return df


def load_rep_seqs_fasta(rep_seqs_dir: Path) -> Dict[str, str]:
    """Load rep-seqs FASTA into a dict of {feature_id: sequence}."""
    fasta_path = rep_seqs_dir / "dna-sequences.fasta"
    if not fasta_path.exists():
        raise FileNotFoundError(f"dna-sequences.fasta not found in {rep_seqs_dir}")
    seqs = {}
    current_id = None
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            current_id = line[1:].strip()
        elif current_id:
            seqs[current_id] = line.strip()
    return seqs


# ---------------------------------------------------------------------------
# Target selection
# ---------------------------------------------------------------------------

def find_target_asvs(
    counts_tsv: Path,
    taxonomy_df: pd.DataFrame,
    taxa_patterns: Optional[List[str]],
    min_reads: Optional[int],
) -> Dict[str, Tuple[str, int]]:
    """
    Find ASV feature IDs matching the requested taxa or read threshold.

    Returns a dict of {feature_id: (taxon_string, total_reads)}.

    Args:
        counts_tsv:      Path to the taxonomy counts TSV from 07_taxonomy_table.py
        taxonomy_df:     DataFrame from exported taxonomy.tsv (feature_id -> Taxon)
        taxa_patterns:   List of substrings to match against taxon strings
        min_reads:       Include any taxon with total reads >= this threshold

    The function builds a reverse map from taxon string to feature IDs,
    then filters based on the requested criteria.
    """
    counts = pd.read_csv(counts_tsv, sep="\t", index_col=0)
    # Total reads per taxon row
    sample_cols = [c for c in counts.columns if not c.startswith("#")]
    totals = counts[sample_cols].sum(axis=1)

    # Build taxon string -> list of feature IDs map
    # taxonomy_df is indexed by feature_id; Taxon column has the assignment
    taxon_col = "Taxon" if "Taxon" in taxonomy_df.columns else taxonomy_df.columns[0]

    # The counts table uses taxon strings as index, not feature IDs.
    # We need to map taxon strings back to feature IDs via the taxonomy TSV.
    taxon_to_features: Dict[str, List[str]] = {}
    for fid, row in taxonomy_df.iterrows():
        taxon = str(row[taxon_col])
        if taxon not in taxon_to_features:
            taxon_to_features[taxon] = []
        taxon_to_features[taxon].append(fid)

    targets: Dict[str, Tuple[str, int]] = {}

    for taxon_str, total in totals.items():
        total = int(total)
        taxon_str = str(taxon_str)

        # Check if this taxon matches any of our criteria
        matches = False

        if taxa_patterns:
            for pat in taxa_patterns:
                if pat.lower() in taxon_str.lower():
                    matches = True
                    break

        if min_reads and total >= min_reads:
            matches = True

        if matches:
            # Find the feature IDs for this taxon
            # Try exact match first, then substring
            feature_ids = taxon_to_features.get(taxon_str, [])
            if not feature_ids:
                # Try matching by substring of taxon string
                for tax, fids in taxon_to_features.items():
                    if taxon_str in tax or tax in taxon_str:
                        feature_ids = fids
                        break

            for fid in feature_ids[:1]:  # Take first ASV per taxon
                targets[fid] = (taxon_str, total)

    return targets


# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------

def write_query_fasta(
    targets: Dict[str, Tuple[str, int]],
    seqs: Dict[str, str],
    fasta_path: Path,
) -> List[str]:
    """
    Write target ASV sequences to a FASTA file for BLASTn.

    Returns list of feature IDs that were successfully written.
    Labels each sequence with a short ID for easy cross-referencing.
    """
    written = []
    with fasta_path.open("w") as f:
        for fid, (taxon, total_reads) in targets.items():
            seq = seqs.get(fid)
            if not seq:
                log.warning("Sequence not found for feature %s — skipping", fid[:12])
                continue
            # Short label: first 8 chars of hash + taxon family/genus
            parts = taxon.replace("d__", "").replace("p__", "").replace(
                "c__", "").replace("o__", "").replace("f__", "").replace(
                "g__", "").replace("s__", "").split(";")
            short_taxon = "_".join(p.strip() for p in parts[-2:] if p.strip())[:30]
            label = f"{fid[:8]}_{short_taxon}_{total_reads}reads"
            f.write(f">{label}\n{seq}\n")
            written.append(fid)
    log.info("Written %d query sequences to %s", len(written), fasta_path)
    return written


def run_blast(
    query_fasta: Path,
    db: Optional[str],
    remote: bool,
    outfile: Path,
    threads: int = 8,
    max_hits: int = 5,
    perc_identity: float = 85.0,
    evalue: float = 0.001,
) -> bool:
    """
    Run BLASTn against local database or remote NCBI.

    Returns True on success, False on failure.

    Uses output format 6 (tabular) with custom fields for easy parsing.
    The perc_identity threshold is intentionally permissive (85%) to catch
    cases where the correct species isn't in the database but a close relative
    is — useful for novel species detection.
    """
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-max_target_seqs", str(max_hits),
        "-outfmt", "6 qseqid sseqid stitle pident length evalue bitscore",
        "-out", str(outfile),
        "-perc_identity", str(perc_identity),
        "-evalue", str(evalue),
    ]

    if remote:
        cmd += ["-db", "nt", "-remote"]
        log.info("Running remote BLASTn against NCBI nt (this may take several minutes)...")
    elif db:
        cmd += ["-db", db, "-num_threads", str(threads)]
        log.info("Running local BLASTn against %s with %d threads...", db, threads)
    else:
        log.error("Either --db or --remote must be specified")
        return False

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error("BLASTn failed: %s", result.stderr[:500])
        return False

    n_hits = sum(1 for line in outfile.read_text().splitlines() if line.strip())
    log.info("BLAST complete: %d hits written to %s", n_hits, outfile)
    return True


# ---------------------------------------------------------------------------
# Result parsing and report
# ---------------------------------------------------------------------------

def parse_blast_results(
    blast_tsv: Path,
    targets: Dict[str, Tuple[str, int]],
) -> pd.DataFrame:
    """
    Parse BLASTn tabular output and join with classifier assignments.

    For each query, takes the top hit by bitscore and compares it to the
    original classifier assignment. Determines AGREE/DISAGREE/NO_HIT status
    and flags implausible taxa.
    """
    cols = ["qseqid", "sseqid", "stitle", "pident", "length", "evalue", "bitscore"]
    blast_text = blast_tsv.read_text().strip()

    if not blast_text:
        log.warning("BLAST output is empty — no hits found above thresholds")
        blast_df = pd.DataFrame(columns=cols)
    else:
        blast_df = pd.read_csv(blast_tsv, sep="\t", header=None, names=cols)

    # Keep only top hit per query
    if not blast_df.empty:
        blast_df = blast_df.sort_values("bitscore", ascending=False)
        blast_df = blast_df.drop_duplicates(subset=["qseqid"], keep="first")

    rows = []
    for fid, (classifier_taxon, total_reads) in targets.items():
        # Find corresponding BLAST row by matching feature ID prefix in qseqid
        fid_short = fid[:8]
        blast_row = blast_df[blast_df["qseqid"].str.startswith(fid_short)]

        if blast_row.empty:
            status = "NO_HIT"
            blast_hit = ""
            pident = ""
            evalue_val = ""
        else:
            r = blast_row.iloc[0]
            blast_hit = str(r["stitle"])[:80]
            pident = r["pident"]
            evalue_val = r["evalue"]

            # Compare classifier to BLAST at genus level
            classifier_genus = _extract_genus(classifier_taxon)
            blast_genus = _extract_genus_from_title(blast_hit)

            # Check for implausible taxa
            is_implausible = any(
                imp.lower() in classifier_taxon.lower()
                for imp in IMPLAUSIBLE_TAXA
            )

            if is_implausible:
                status = "ARTEFACT_FLAG"
            elif classifier_genus and blast_genus:
                if classifier_genus.lower() == blast_genus.lower():
                    status = "AGREE"
                elif _same_family(classifier_taxon, blast_hit):
                    status = "AGREE_FAMILY"
                else:
                    status = "DISAGREE"
            else:
                status = "UNRESOLVED"

        rows.append({
            "feature_id":         fid,
            "classifier_taxon":   classifier_taxon,
            "total_reads":        total_reads,
            "blast_top_hit":      blast_hit,
            "pident":             pident,
            "evalue":             evalue_val,
            "status":             status,
            "action_recommended": _recommend_action(status, classifier_taxon, blast_hit),
        })

    return pd.DataFrame(rows)


def _extract_genus(taxon_str: str) -> str:
    """Extract genus from a Greengenes-style or NCBI-style taxonomy string."""
    # Greengenes: ...;g__Brevoortia;s__tyrannus
    for part in taxon_str.split(";"):
        part = part.strip()
        if part.startswith("g__"):
            return part[3:]
    # Plain: Brevoortia tyrannus
    parts = taxon_str.split()
    if len(parts) >= 1:
        return parts[-2] if len(parts) >= 2 else parts[0]
    return ""


def _extract_genus_from_title(blast_title: str) -> str:
    """Extract likely genus from a BLAST hit title string."""
    # BLAST titles look like: "Brevoortia tyrannus voucher ... 12S ..."
    # First word after removing accession is usually genus
    parts = blast_title.strip().split()
    # Skip accession if present (all caps + numbers)
    for i, p in enumerate(parts):
        if p[0].isupper() and not p.isupper():
            return p
    return parts[0] if parts else ""


def _same_family(classifier_taxon: str, blast_title: str) -> bool:
    """Check if classifier family appears in BLAST title (loose match)."""
    for part in classifier_taxon.split(";"):
        part = part.strip()
        if part.startswith("f__"):
            family = part[3:].lower()
            if family and family in blast_title.lower():
                return True
    return False


def _recommend_action(status: str, classifier_taxon: str, blast_hit: str) -> str:
    """Generate a plain-English action recommendation based on BLAST status."""
    if status == "AGREE":
        return "Classifier confirmed — keep in table"
    elif status == "AGREE_FAMILY":
        return "Same family, different genus — check species-level reference coverage"
    elif status == "DISAGREE":
        return f"BLAST disagrees — review manually; consider adding to artefact list"
    elif status == "NO_HIT":
        return "No BLAST hit — sequence may be novel or database gap; keep but note uncertainty"
    elif status == "ARTEFACT_FLAG":
        return "Geographically/ecologically implausible — add to artefact exclusion list"
    elif status == "UNRESOLVED":
        return "Could not compare — check sequence quality"
    return ""


def write_report(results: pd.DataFrame, report_path: Path, marker: str) -> None:
    """Write a human-readable verification report."""
    with report_path.open("w") as f:
        f.write(f"BLAST Verification Report — {marker}\n")
        f.write("=" * 70 + "\n\n")

        status_counts = results["status"].value_counts()
        f.write("Summary:\n")
        for status, count in status_counts.items():
            f.write(f"  {status:<20} {count}\n")
        f.write("\n")

        # Flag artefacts and disagreements prominently
        problems = results[results["status"].isin(
            ["ARTEFACT_FLAG", "DISAGREE", "NO_HIT"]
        )]
        if not problems.empty:
            f.write("ACTION REQUIRED:\n")
            f.write("-" * 70 + "\n")
            for _, row in problems.iterrows():
                f.write(f"\n  Taxon:      {row['classifier_taxon']}\n")
                f.write(f"  Reads:      {row['total_reads']:,}\n")
                f.write(f"  BLAST hit:  {row['blast_top_hit']}\n")
                f.write(f"  Identity:   {row['pident']}%\n")
                f.write(f"  Status:     {row['status']}\n")
                f.write(f"  Action:     {row['action_recommended']}\n")
            f.write("\n")

        # Confirmed taxa
        confirmed = results[results["status"].isin(["AGREE", "AGREE_FAMILY"])]
        if not confirmed.empty:
            f.write("CONFIRMED TAXA:\n")
            f.write("-" * 70 + "\n")
            for _, row in confirmed.iterrows():
                f.write(f"  {row['classifier_taxon'].split(';')[-1].strip():<35} "
                        f"{row['pident']}% identity — {row['blast_top_hit'][:50]}\n")

    log.info("Report written to %s", report_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="blast_verify.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    req = p.add_argument_group("inputs (derived from --marker + config if omitted)")
    req.add_argument("--marker",    required=True,
                     help="Marker name (e.g. MiFish, cytb). Drives derived paths and config.")
    req.add_argument("--config",    default=None,
                     help="Path to pipeline_config.yml.")
    req.add_argument("--taxonomy",  default=None, type=Path,
                     help="QIIME2 taxonomy .qza (exported internally). Derived if omitted.")
    req.add_argument("--rep-seqs",  default=None, type=Path,
                     help="QIIME2 rep-seqs .qza (DADA2 output). Derived if omitted.")
    req.add_argument("--counts",    default=None, type=Path,
                     help="Taxonomy count TSV from taxonomy_table.py. Derived if omitted.")
    req.add_argument("--outdir",    default=None, type=Path,
                     help="Output directory. Default: results/<marker>/all/blast/")

    sel = p.add_argument_group("target selection (one or both; default from config)")
    sel.add_argument("--taxa", nargs="+", default=None,
                     help="Taxon name substrings to BLAST (e.g. Brevoortia Lutjanidae).")
    sel.add_argument("--min-reads", type=int, default=None,
                     help="BLAST all taxa with total reads >= this threshold "
                          "(default: analyses.blast.verify.min_reads).")

    db_grp = p.add_argument_group("BLAST database (--db from analyses.blast.db, or --remote)")
    db_grp.add_argument("--db", default=None,
                         help="Path to local BLAST nt database (no extension). "
                              "Default: analyses.blast.db from config.")
    db_grp.add_argument("--remote", action="store_true", default=False,
                         help="Use remote NCBI BLAST instead of a local database (slower).")

    opt = p.add_argument_group("BLAST parameters (default from analyses.blast.verify)")
    opt.add_argument("--threads",       type=int,   default=None,
                     help="CPU threads for local BLAST. Default: blast.num_threads.")
    opt.add_argument("--max-hits",      type=int,   default=None,
                     help="Maximum BLAST hits per query. Default: 5 (config: verify.max_hits).")
    opt.add_argument("--perc-identity", type=float, default=None,
                     help="Minimum percent identity for BLAST hits. Default: 85 (config: verify.perc_identity).")
    opt.add_argument("--evalue",        type=float, default=None,
                     help="Maximum e-value for BLAST hits. Default: 0.001 (config: verify.evalue).")
    opt.add_argument("--dry-run", action="store_true",
                     help="Export sequences and write FASTA but do not run BLAST.")
    return p


def _derive_counts_tsv(paths, marker: str) -> Optional[Path]:
    """
    Find the taxonomy count TSV for a marker, picking the deepest taxonomic
    level available (e.g. L7 over L6). Returns None if none exist yet.
    """
    import re
    tax_dir = paths.engine_taxonomy_results_dir(marker, "all")
    hits = sorted(tax_dir.glob(f"taxonomy_counts_L*_{marker}.tsv"))
    if not hits:
        return None
    def _level(p: Path) -> int:
        m = re.search(r"_L(\d+)_", p.name)
        return int(m.group(1)) if m else 0
    return max(hits, key=_level)


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    bcfg = cfg.analyses.get("blast", {})
    vcfg = bcfg.get("verify", {})       # tool-specific knobs for blast_verify
    marker = args.marker

    # Derive inputs/outputs from the marker + config when not given explicitly.
    if args.taxonomy is None:
        args.taxonomy = paths.engine_taxonomy_qza(marker, "all")
    if args.rep_seqs is None:
        args.rep_seqs = paths.engine_rep_seqs_qza(marker, "all")
    if args.counts is None:
        args.counts = _derive_counts_tsv(paths, marker)
    if args.outdir is None:
        args.outdir = paths.engine_blast_results_dir(marker, "all")
    if args.db is None and not args.remote:
        args.db = bcfg.get("db") or None
    # Config fills unset knobs (CLI still wins).
    if args.threads is None:
        args.threads = bcfg.get("num_threads", 8)
    if args.max_hits is None:
        args.max_hits = vcfg.get("max_hits", 5)
    if args.perc_identity is None:
        args.perc_identity = vcfg.get("perc_identity", 85.0)
    if args.evalue is None:
        args.evalue = vcfg.get("evalue", 0.001)
    if args.taxa is None and args.min_reads is None:
        args.min_reads = vcfg.get("min_reads")        # config default target selection
        args.taxa = vcfg.get("taxa")

    # Validate inputs
    validate.require_qiime()
    if args.counts is None:
        log.error("No taxonomy count TSV found for %s — run the taxonomy stage "
                  "first, or pass --counts.", marker)
        return 2
    for path, name in [
        (args.taxonomy,  "--taxonomy"),
        (args.rep_seqs,  "--rep-seqs"),
        (args.counts,    "--counts"),
    ]:
        if not Path(path).exists():
            log.error("%s not found: %s", name, path)
            return 2

    if not args.taxa and not args.min_reads:
        log.error("No targets selected: set analyses.blast.verify.min_reads (or "
                  ".taxa) in the config, or pass --taxa / --min-reads.")
        return 2

    if not args.remote and not args.db:
        log.error("No BLAST database: set analyses.blast.db, pass --db, or use --remote.")
        return 2

    if not args.remote and args.db and \
       not Path(str(args.db) + ".nal").exists() and \
       not Path(str(args.db) + ".nsi").exists():
        log.warning("Local BLAST database not found at %s", args.db)
        log.warning("Use --remote for NCBI BLAST, or check --db path")
        log.error("Cannot proceed without a valid database. Add --remote to use NCBI instead.")
        return 2

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Export QIIME2 artifacts
    log.info("Exporting taxonomy artifact...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        tax_dir = export_qza(args.taxonomy, tmp / "taxonomy")
        taxonomy_df = load_taxonomy_tsv(tax_dir)
        log.info("Loaded %d feature taxonomy assignments", len(taxonomy_df))

        seq_dir = export_qza(args.rep_seqs, tmp / "rep_seqs")
        seqs = load_rep_seqs_fasta(seq_dir)
        log.info("Loaded %d representative sequences", len(seqs))

        # Find target ASVs
        log.info("Selecting target ASVs...")
        targets = find_target_asvs(
            args.counts, taxonomy_df, args.taxa, args.min_reads
        )
        log.info("Found %d target ASVs to BLAST", len(targets))

        if not targets:
            log.warning("No matching taxa found. Check --taxa patterns or --min-reads threshold.")
            return 0

        # Write query FASTA
        fasta_out = args.outdir / f"blast_query_{args.marker}.fasta"
        written_ids = write_query_fasta(targets, seqs, fasta_out)

        if args.dry_run:
            log.info("DRY RUN — sequences written to %s, BLAST not run", fasta_out)
            log.info("To BLAST manually, upload %s to https://blast.ncbi.nlm.nih.gov/", fasta_out)
            return 0

        # Run BLAST
        blast_out = args.outdir / f"blast_raw_{args.marker}.tsv"
        success = run_blast(
            fasta_out,
            db=None if args.remote else args.db,
            remote=args.remote,
            outfile=blast_out,
            threads=args.threads,
            max_hits=args.max_hits,
            perc_identity=args.perc_identity,
            evalue=args.evalue,
        )

        if not success:
            log.error("BLAST failed — check logs above")
            return 1

        # Parse results
        results = parse_blast_results(blast_out, targets)

        # Write outputs
        results_path = args.outdir / f"blast_verify_{args.marker}.tsv"
        results.to_csv(results_path, sep="\t", index=False)
        log.info("Results table: %s", results_path)

        report_path = args.outdir / f"blast_verify_{args.marker}_report.txt"
        write_report(results, report_path, args.marker)

        # Print summary to console
        print("\n" + "=" * 60)
        print(f"  BLAST Verification — {args.marker}")
        print("=" * 60)
        for _, row in results.iterrows():
            taxon_short = row["classifier_taxon"].split(";")[-1].strip()[:35]
            print(f"  {row['status']:<18} {taxon_short:<35} "
                  f"{row['total_reads']:>8,} reads")
            if row["status"] in ("ARTEFACT_FLAG", "DISAGREE"):
                print(f"    → {row['action_recommended']}")
        print("=" * 60)
        print(f"\nFull report: {report_path}")

        produced = [str(results_path), str(report_path)]

    checkpoint.print_checkpoint(
        cfg, "blast",
        marker=marker,
        produced=produced,
        provenance={
            "inputs": {"taxonomy": str(args.taxonomy), "rep_seqs": str(args.rep_seqs),
                       "counts": str(args.counts),
                       "blast_db": str(args.db) if args.db else None},
            "outputs": produced,
            "command": "python " + " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
            "extra": {"perc_identity": args.perc_identity, "max_hits": args.max_hits,
                      "remote": bool(args.remote),
                      "min_reads": args.min_reads,
                      "taxa": args.taxa},
        },
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except validate.ValidationError as e:
        log.error("%s", e)
        sys.exit(1)
