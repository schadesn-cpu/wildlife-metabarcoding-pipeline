#!/usr/bin/env python3
"""
07e_entrez_verify_accessions.py

Verify BLAST-assigned taxonomy by looking up ALL hit accessions via NCBI Entrez.
Produces a verified remap table and clean taxonomy for downstream analysis.

Designed to work with outputs from:
  - 07b_blast_verify.py (blast_verify_{marker}.tsv)
  - 07c_blast_qc_unclassified.py (blast_qc_{marker}.tsv)
  - 07d_blast_refine_unresolved.py (blast_refine_{marker}.tsv)
  - NCBI web BLAST CSV downloads (Alignment-HitTable.csv)

Workflow
--------
  1. Parse BLAST results (supports multiple input formats)
  2. Extract top-hit accession per query
  3. Batch-lookup all accessions via Bio.Entrez
  4. Compare classifier label vs BLAST title vs Entrez-verified species
  5. Flag mismatches and generate remap table
  6. Optionally apply remapping to produce clean taxonomy

Usage
-----
  # Verify BLAST results from 07b
  python 07e_entrez_verify_accessions.py \\
      --blast-results results/cytb/all/blast_verify/blast_verify_cytb.tsv \\
      --format pipeline \\
      --taxonomy qiime2/cytb/all/taxonomy-exported-temp/taxonomy.tsv \\
      --marker cytb \\
      --email your@email.edu \\
      --outdir results/cytb/all/blast_verify/

  # Verify NCBI web BLAST CSV
  python 07e_entrez_verify_accessions.py \\
      --blast-results YN4CKJ3D014-Alignment-HitTable.csv \\
      --format ncbi-csv \\
      --taxonomy qiime2/cytb/all/taxonomy-exported-temp/taxonomy.tsv \\
      --marker cytb \\
      --email your@email.edu \\
      --outdir results/cytb/all/blast_verify/

  # Apply remap table to produce clean taxonomy
  python 07e_entrez_verify_accessions.py \\
      --apply-remap results/cytb/all/blast_verify/verified_remap.tsv \\
      --taxonomy qiime2/cytb/all/taxonomy-exported-temp/taxonomy.tsv \\
      --outdir results/cytb/all/taxonomy_verified/
"""

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Entrez accession lookup
# ---------------------------------------------------------------------------

def verify_accessions(
    accessions: List[str],
    email: str,
    batch_size: int = 20,
    delay: float = 0.4,
) -> pd.DataFrame:
    """
    Look up accession numbers via NCBI Entrez and return verified species info.
    
    Returns DataFrame with columns:
        accession, organism, title, mol_type, length
    """
    try:
        from Bio import Entrez
    except ImportError:
        log.error("Biopython not installed. Run: pip install biopython")
        sys.exit(1)

    Entrez.email = email
    unique_accs = list(set(a for a in accessions if a and a != "N/A"))
    log.info("Verifying %d unique accessions via NCBI Entrez...", len(unique_accs))

    records = []
    failed = []

    for i in range(0, len(unique_accs), batch_size):
        batch = unique_accs[i:i + batch_size]
        try:
            handle = Entrez.efetch(
                db="nucleotide", id=batch,
                rettype="docsum", retmode="xml"
            )
            result = Entrez.read(handle)
            handle.close()

            for rec in result:
                org = rec.get("Organism", "")
                title = rec.get("Title", "")
                # Organism field is often empty in nucleotide docsums
                # Extract from Title: "Ursus americanus cytochrome b..." → "Ursus americanus"
                if not org and title:
                    import re
                    # Remove leading "PREDICTED: " if present
                    clean_title = re.sub(r'^PREDICTED:\s*', '', title)
                    # Remove leading "Uncultured " prefix
                    clean_title = re.sub(r'^[Uu]ncultured\s+', '', clean_title)
                    # Split on common delimiters and take first 2 words as organism
                    # Stop before: gene, isolate, voucher, mitochondrion, genome, DNA,
                    # chromosome, strain, haplotype, clone, small, large, partial, complete
                    stop_words = r'\b(gene|isolate|voucher|mitochondrion|mitochondrial|genome|DNA|dna|'
                    stop_words += r'chromosome|strain|haplotype|clone|small|large|partial|complete|'
                    stop_words += r'cytochrome|ribosomal|16S|18S|28S|cytb|COI|NADH)\b'
                    m = re.split(stop_words, clean_title, maxsplit=1)
                    org = m[0].strip().rstrip(',').strip() if m else clean_title[:50]
                    # Sanity: if it looks like just 1 word, keep up to 3
                    if ' ' not in org:
                        words = clean_title.split()[:3]
                        org = ' '.join(words)

                records.append({
                    "accession": rec.get("AccessionVersion", ""),
                    "organism": org,
                    "title": title,
                    "length": rec.get("Length", ""),
                    "mol_type": rec.get("MolType", ""),
                })

            time.sleep(delay)  # Be polite to NCBI

        except Exception as e:
            log.warning("Entrez batch %d failed: %s", i, e)
            failed.extend(batch)

    if failed:
        log.warning("Failed to resolve %d accessions: %s",
                    len(failed), failed[:5])
        for acc in failed:
            records.append({
                "accession": acc,
                "organism": "LOOKUP_FAILED",
                "title": "LOOKUP_FAILED",
                "length": "",
                "mol_type": "",
            })

    df = pd.DataFrame(records)
    log.info("Verified %d / %d accessions", len(df) - len(failed), len(unique_accs))
    return df


# ---------------------------------------------------------------------------
# BLAST result parsers
# ---------------------------------------------------------------------------

def parse_pipeline_blast(path: Path) -> pd.DataFrame:
    """
    Parse TSV from 07b/07c/07d pipeline scripts.
    Auto-detects column format from header or column count.
    """
    # Try reading with header first
    try:
        df = pd.read_csv(path, sep="\t")
        if "qseqid" in df.columns or "query_id" in df.columns:
            return df
    except Exception:
        pass

    # Try standard BLAST outfmt 6 columns
    col_options = {
        7: ["qseqid", "sseqid", "stitle", "pident", "length", "evalue", "bitscore"],
        8: ["qseqid", "sseqid", "staxid", "pident", "length", "evalue", "bitscore", "stitle"],
        9: ["qseqid", "sacc", "sscinames", "staxids", "pident", "length", "evalue", "bitscore", "stitle"],
    }

    df = pd.read_csv(path, sep="\t", header=None)
    ncols = len(df.columns)

    if ncols in col_options:
        df.columns = col_options[ncols]
    else:
        # Best guess
        log.warning("Unexpected column count (%d). Assuming first 7 standard BLAST cols.", ncols)
        base_cols = ["qseqid", "sseqid", "pident", "length", "mismatch",
                     "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        df.columns = base_cols[:ncols] + [f"col_{i}" for i in range(ncols - len(base_cols))] \
            if ncols > len(base_cols) else base_cols[:ncols]

    return df


def parse_ncbi_csv(path: Path) -> pd.DataFrame:
    """Parse NCBI web BLAST CSV download (Alignment-HitTable.csv)."""
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    df = pd.read_csv(path, names=cols, comment="#")
    return df


# ---------------------------------------------------------------------------
# Comparison and flagging logic
# ---------------------------------------------------------------------------

def extract_label_from_query(query_id: str) -> str:
    """Extract the classifier label from a query ID like 'abc123__Rickettsia__s__sp.'."""
    parts = query_id.split("__")
    if len(parts) >= 2:
        return "__".join(parts[1:])
    return query_id


def classify_match(
    classifier_label: str,
    blast_title: str,
    entrez_organism: str,
    entrez_title: str,
    pident: float,
    length: int,
) -> Tuple[str, str, str]:
    """
    Compare classifier, BLAST, and Entrez results.
    Returns (status, verified_species, action).
    
    status: confirmed, mismatch, artifact, contamination, review
    action: keep, remap, exclude, review
    """
    label_lower = classifier_label.lower()
    entrez_org_lower = entrez_organism.lower()
    entrez_title_lower = entrez_title.lower()

    # Check for known artifacts
    if length < 50:
        return "artifact", f"Short alignment ({length}bp)", "exclude"

    if "adapter" in entrez_title_lower or "GATCGGAAGAGC" in entrez_title:
        return "artifact", "Adapter contamination", "exclude"

    # Check for contamination
    if "homo sapiens" in entrez_org_lower and "sapiens" not in label_lower:
        return "contamination", "Homo sapiens", "exclude"

    if any(x in entrez_org_lower for x in ["sus scrofa", "gallus gallus"]) \
            and not any(x in label_lower for x in ["sus", "gallus", "scrofa"]):
        return "contamination", entrez_organism, "exclude"

    # Check for pseudogenes / NUMTs
    if "pseudogene" in entrez_title_lower:
        return "pseudogene", f"{entrez_organism} (pseudogene)", "review"

    if "genome assembly, chromosome" in entrez_title_lower:
        return "numt", f"{entrez_organism} (NUMT)", "review"

    # Check for no significant hit
    if pident < 80:
        return "weak_hit", f"{entrez_organism} ({pident:.1f}%)", "review"

    if entrez_organism == "LOOKUP_FAILED":
        return "lookup_failed", "LOOKUP_FAILED", "review"

    # Check if classifier genus matches Entrez organism
    # Extract genus from both
    # Handle labels like "s__porcinus annamiticus__n=1078" or "g__Hydromantes__n=49"
    import re
    # Strip n=xxx suffix and rank prefixes
    clean_label = re.sub(r'__n=\d+$', '', classifier_label)
    clean_label = clean_label.split("__")[-1]  # Get last meaningful part
    for prefix in ["g__", "s__", "k__", "p__", "c__", "o__", "f__"]:
        clean_label = clean_label.replace(prefix, "")
    classifier_genus = clean_label.strip().split()[0].split("_")[0].lower()
    
    entrez_genus = entrez_org_lower.split()[0] if entrez_org_lower else ""

    if classifier_genus and entrez_genus and classifier_genus == entrez_genus:
        return "confirmed", entrez_organism, "keep"
    elif pident >= 97:
        return "mismatch_high_id", entrez_organism, "remap"
    elif pident >= 90:
        return "mismatch_moderate_id", entrez_organism, "remap"
    else:
        return "mismatch_low_id", entrez_organism, "review"


# ---------------------------------------------------------------------------
# Remap table application
# ---------------------------------------------------------------------------

def apply_remap(
    taxonomy_path: Path,
    remap_path: Path,
    outdir: Path,
) -> None:
    """Apply a verified remap table to a taxonomy TSV."""
    tax = pd.read_csv(taxonomy_path, sep="\t", comment="#",
                      names=["Feature_ID", "Taxon", "Confidence"], header=0)
    remap = pd.read_csv(remap_path, sep="\t")

    log.info("Loaded %d ASVs and %d remap rules", len(tax), len(remap))

    tax["Taxon_verified"] = None
    tax["Action"] = None

    remapped = excluded = reviewed = 0

    for _, rule in remap.iterrows():
        search = rule["search_string"]
        action = rule["action"]
        verified = rule["verified_species"]

        mask = tax["Taxon"].str.contains(search, case=False, na=False)
        n = mask.sum()
        if n == 0:
            continue

        if action == "remap":
            tax.loc[mask, "Taxon_verified"] = verified
            tax.loc[mask, "Action"] = "remapped"
            remapped += n
            log.info("  REMAP: %s -> %s (%d ASVs)", search, verified, n)
        elif action == "exclude":
            tax.loc[mask, "Taxon_verified"] = "EXCLUDED"
            tax.loc[mask, "Action"] = "excluded"
            excluded += n
            log.info("  EXCLUDE: %s (%d ASVs)", search, n)
        elif action == "review":
            tax.loc[mask, "Taxon_verified"] = f"REVIEW:{verified}"
            tax.loc[mask, "Action"] = "review"
            reviewed += n
            log.info("  REVIEW: %s (%d ASVs)", search, n)

    # Fill unmatched
    unmatched = tax["Taxon_verified"].isna()
    unresolved = unmatched & (tax["Taxon"].str.count(";") < 4)
    confirmed = unmatched & ~unresolved

    tax.loc[unresolved, "Taxon_verified"] = "EXCLUDED"
    tax.loc[unresolved, "Action"] = "unresolved_excluded"
    tax.loc[confirmed, "Taxon_verified"] = tax.loc[confirmed, "Taxon"]
    tax.loc[confirmed, "Action"] = "kept_as_is"

    log.info("Remapped: %d | Excluded: %d | Review: %d | Unresolved: %d | Kept: %d",
             remapped, excluded, reviewed, unresolved.sum(), confirmed.sum())

    # Save outputs
    outdir.mkdir(parents=True, exist_ok=True)

    tax.to_csv(outdir / "taxonomy_blast_verified.tsv", sep="\t", index=False)

    clean = tax[tax["Taxon_verified"] != "EXCLUDED"].copy()
    clean = clean[~clean["Taxon_verified"].str.startswith("REVIEW", na=False)]
    clean = clean[["Feature_ID", "Taxon_verified", "Confidence"]].rename(
        columns={"Taxon_verified": "Taxon"})
    clean.to_csv(outdir / "taxonomy_clean.tsv", sep="\t", index=False)

    summary = clean["Taxon"].value_counts()
    summary.to_csv(outdir / "host_summary.tsv", sep="\t", header=["ASV_count"])

    log.info("Saved: %s/taxonomy_blast_verified.tsv", outdir)
    log.info("Saved: %s/taxonomy_clean.tsv (%d ASVs)", outdir, len(clean))
    log.info("Saved: %s/host_summary.tsv (%d taxa)", outdir, len(summary))


# ---------------------------------------------------------------------------
# Main verify workflow
# ---------------------------------------------------------------------------

def verify_workflow(args) -> None:
    """Full verification workflow: parse BLAST → Entrez lookup → flag → remap table."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Parse BLAST results
    log.info("Parsing BLAST results: %s (format: %s)", args.blast_results, args.format)
    if args.format == "ncbi-csv":
        blast = parse_ncbi_csv(Path(args.blast_results))
    else:
        blast = parse_pipeline_blast(Path(args.blast_results))

    log.info("Loaded %d BLAST hits for %d queries",
             len(blast), blast["qseqid"].nunique())

    # 2. Get top hit per query
    if "bitscore" in blast.columns:
        top = blast.sort_values("bitscore", ascending=False).groupby("qseqid").first().reset_index()
    else:
        top = blast.groupby("qseqid").first().reset_index()

    log.info("Top hits: %d queries", len(top))

    # 3. Get accession column
    acc_col = "sacc" if "sacc" in top.columns else "sseqid"
    accessions = top[acc_col].unique().tolist()

    # 4. Entrez verification
    entrez_df = verify_accessions(accessions, args.email)

    # Save raw Entrez results
    entrez_df.to_csv(outdir / f"entrez_verified_{args.marker}.tsv", sep="\t", index=False)

    # 5. Merge BLAST + Entrez (rename Entrez cols to avoid collision)
    entrez_df = entrez_df.rename(columns={
        "length": "entrez_seq_length",
        "title": "entrez_title",
        "organism": "entrez_organism",
        "mol_type": "entrez_mol_type",
    })
    top = top.merge(entrez_df, left_on=acc_col, right_on="accession", how="left")

    # 6. Classify each hit
    results = []
    for _, row in top.iterrows():
        query = row["qseqid"]
        classifier_label = extract_label_from_query(query)
        blast_title = row.get("stitle", "")
        entrez_org = row.get("entrez_organism", "")
        entrez_title = row.get("entrez_title", "")
        pident = float(row.get("pident", 0))
        length = int(row.get("length", 0))

        status, verified, action = classify_match(
            classifier_label, str(blast_title),
            str(entrez_org), str(entrez_title),
            pident, length
        )

        results.append({
            "query_id": query,
            "classifier_label": classifier_label,
            "top_accession": row[acc_col],
            "entrez_organism": entrez_org,
            "entrez_title": str(entrez_title)[:80],
            "pident": pident,
            "length": length,
            "status": status,
            "verified_species": verified,
            "action": action,
        })

    results_df = pd.DataFrame(results)

    # 7. Save verification table
    results_df.to_csv(outdir / f"verification_table_{args.marker}.tsv",
                      sep="\t", index=False)

    # 8. Generate remap table
    remap_rows = []
    for _, row in results_df.iterrows():
        # Use the classifier label as search string
        search = row["classifier_label"]
        # Strip n=NNN suffix and rank prefixes to get the taxon name
        import re
        search_clean = re.sub(r'__n=\d+$', '', search)
        # Remove rank prefixes like s__, g__, f__
        search_clean = re.sub(r'^[a-z]__', '', search_clean)
        # If still has __ parts, take the last meaningful one
        parts = [p for p in search_clean.split("__") if p and not p.startswith("n=")]
        search_simple = parts[-1] if parts else search_clean

        if row["action"] in ("remap", "exclude", "review"):
            remap_rows.append({
                "search_string": search_simple,
                "verified_species": row["verified_species"],
                "action": row["action"],
                "notes": f"{row['top_accession']} = {row['entrez_organism']} "
                         f"({row['pident']:.1f}%, {row['length']}bp) [{row['status']}]",
            })

    remap_df = pd.DataFrame(remap_rows)
    remap_path = outdir / f"verified_remap_{args.marker}.tsv"
    remap_df.to_csv(remap_path, sep="\t", index=False)

    # 9. Print summary
    print(f"\n{'='*70}")
    print(f"ENTREZ VERIFICATION SUMMARY — {args.marker}")
    print(f"{'='*70}")
    print(f"\nTotal queries:  {len(results_df)}")
    for status, count in results_df["status"].value_counts().items():
        print(f"  {status:<25} {count:>4}")

    print(f"\nAction breakdown:")
    for action, count in results_df["action"].value_counts().items():
        print(f"  {action:<10} {count:>4}")

    print(f"\nOutputs:")
    print(f"  Verification table: {outdir}/verification_table_{args.marker}.tsv")
    print(f"  Entrez lookups:     {outdir}/entrez_verified_{args.marker}.tsv")
    print(f"  Remap table:        {remap_path}")

    # Show flagged items
    flagged = results_df[results_df["action"] != "keep"]
    if len(flagged) > 0:
        print(f"\n{'='*70}")
        print("FLAGGED ITEMS (review before applying remap):")
        print(f"{'='*70}")
        for _, row in flagged.iterrows():
            icon = {"remap": "🔄", "exclude": "❌", "review": "⚠️"}.get(row["action"], "?")
            print(f"  {icon} [{row['action']:>7}] {row['classifier_label'][:30]:<30} "
                  f"→ {row['entrez_organism'][:35]:<35} "
                  f"({row['pident']:.1f}%, {row['status']})")

    print(f"\n{'='*70}")
    print("NEXT STEPS:")
    print(f"  1. Review flagged items above")
    print(f"  2. Edit {remap_path} if needed")
    print(f"  3. Apply with:")
    print(f"     python 07e_entrez_verify_accessions.py \\")
    print(f"       --apply-remap {remap_path} \\")
    print(f"       --taxonomy <your_taxonomy.tsv> \\")
    print(f"       --outdir <output_dir>/")
    print(f"{'='*70}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="07e_entrez_verify_accessions.py",
        description="Verify BLAST accessions via Entrez and generate remap tables.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    sub = p.add_subparsers(dest="command")

    # --- verify subcommand ---
    v = sub.add_parser("verify", help="Verify BLAST results via Entrez lookup")
    v.add_argument("--blast-results", required=True, type=Path,
                   help="BLAST results file (TSV from pipeline or CSV from NCBI web)")
    v.add_argument("--format", choices=["pipeline", "ncbi-csv", "auto"],
                   default="auto",
                   help="Input format: pipeline (07b/c/d output), ncbi-csv, or auto-detect")
    v.add_argument("--taxonomy", type=Path, default=None,
                   help="Exported taxonomy.tsv (optional, for context)")
    v.add_argument("--marker", required=True,
                   help="Marker gene name (16S, 18S, cytb, etc.)")
    v.add_argument("--email", required=True,
                   help="Email for NCBI Entrez (required by NCBI)")
    v.add_argument("--outdir", required=True, type=Path,
                   help="Output directory")

    # --- apply subcommand ---
    a = sub.add_parser("apply", help="Apply a verified remap table to taxonomy")
    a.add_argument("--remap", required=True, type=Path,
                   help="Verified remap TSV (from verify step)")
    a.add_argument("--taxonomy", required=True, type=Path,
                   help="Exported taxonomy.tsv")
    a.add_argument("--outdir", required=True, type=Path,
                   help="Output directory")

    # --- common ---
    p.add_argument("-v", "--verbose", action="store_true")

    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if getattr(args, "verbose", False) else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.command == "verify":
        # Auto-detect format
        if args.format == "auto":
            path = str(args.blast_results)
            if path.endswith(".csv") or "HitTable" in path:
                args.format = "ncbi-csv"
            else:
                args.format = "pipeline"
            log.info("Auto-detected format: %s", args.format)

        verify_workflow(args)

    elif args.command == "apply":
        apply_remap(args.taxonomy, args.remap, args.outdir)

    else:
        parser.print_help()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
