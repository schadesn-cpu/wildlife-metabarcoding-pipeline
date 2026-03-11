#!/usr/bin/env python3
"""
00_build_classifiers.py
=======================
QIIME 2 classifier training wrapper for metabarcoding workflows.

Supports:
  - 16S   : SILVA 138 (downloads pre-trained classifier)
  - 18S   : PR2 v5.0.0
  - ITS   : UNITE v10 (developer release)
  - MiFish: MitoFish 12S rRNA via Mitohelper (Mar 2025)
  - cytb  : NCBI vertebrate cytochrome b (via RESCRIPt + Entrez)

Usage:
  python 00_build_classifiers.py --markers 16S 18S ITS MiFish cytb --outdir classifiers/
  python 00_build_classifiers.py --markers cytb --outdir classifiers/ --ncbi-api-key YOUR_KEY
  python 00_build_classifiers.py --markers 16S --outdir classifiers/ --threads 8

Options:
  --markers       One or more of: 16S, 18S, ITS, MiFish, cytb   [required]
  --outdir        Output directory for classifiers                [default: ./classifiers]
  --f-primer      Forward primer (optional, for trimming)
  --r-primer      Reverse primer (optional, for trimming)
  --threads       Number of threads for training                  [default: 4]
  --ncbi-api-key  NCBI Entrez API key (required for cytb)
  --ncbi-email    Email for NCBI Entrez (required for cytb)
  --cytb-taxon    NCBI taxon for cytb query                      [default: Vertebrata]
  --cytb-min-len  Minimum cytb sequence length to retain         [default: 200]
  --cytb-max-len  Maximum cytb sequence length to retain         [default: 1200]
  --skip-download Skip download if files already exist
  --dry-run       Print commands without executing them
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

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
# Database / URL config
# ---------------------------------------------------------------------------
DB_CONFIG = {
    "16S": {
        "description": "SILVA 138 99% OTUs (16S rRNA)",
        # ── Full-length (pre-trained, QIIME2 2024.5) ──────────────────────────
        # 209 MB. Classifies full-length 16S accurately but performs POORLY on
        # V4 amplicons (~0.1% genus resolution vs ~40-60% for V4-trained).
        # Only use this if you sequenced full-length 16S.
        "pretrained_url": "https://data.qiime2.org/2024.5/common/silva-138-99-nb-classifier.qza",
        "pretrained_filename": "silva-138-99-nb-classifier.qza",
        # ── Full-length seqs + taxonomy (for custom primer training) ──────────
        "ref_seqs_url": "https://data.qiime2.org/2024.5/common/silva-138-99-seqs.qza",
        "ref_seqs_filename": "silva-138-99-seqs.qza",
        "ref_tax_url": "https://data.qiime2.org/2024.5/common/silva-138-99-tax.qza",
        "ref_tax_filename": "silva-138-99-tax.qza",
        # ── V4 amplicon (515F/806R) pre-trimmed seqs — RECOMMENDED for V4 ────
        # ~50 MB. Pre-trimmed to the 515F/806R region by the QIIME2 team.
        # Training from these instead of full seqs takes ~2-4 hours and produces
        # a classifier with 40-60% genus resolution on V4 reads vs <1% for full.
        # Use --16s-v4 to select this path.
        "v4_seqs_url": "https://data.qiime2.org/2024.5/common/silva-138-99-seqs-515-806.qza",
        "v4_seqs_filename": "silva-138-99-seqs-515-806.qza",
        "v4_classifier_filename": "silva-138-99-nb-classifier-515-806.qza",
        "pretrained": True,  # Can use pre-trained unless primers are specified
    },
    "18S": {
        "description": "PR2 v5.0.0 (18S rRNA)",
        "fasta_url": "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_QIIME.fasta.gz",
        "tax_url": "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_QIIME.tax.gz",
        "fasta_filename": "pr2_version_5.0.0_SSU_QIIME.fasta",
        "tax_filename": "pr2_version_5.0.0_SSU_QIIME.tax",
        "seqs_qza": "pr2-18S-seqs.qza",
        "tax_qza": "pr2-18S-tax.qza",
        "classifier_qza": "pr2-18S-classifier.qza",
        "tax_format": "HeaderlessTSVTaxonomyFormat",
        "pretrained": False,
    },
    "ITS": {
        "description": "UNITE v10 developer release (ITS1-2, all eukaryotes)",
        "fasta_url": "https://files.plutof.ut.ee/public/orig/EB/E9/EBE9A788A03E4E3B8C53B4E83C36B36E0CC0DCC5BDA2B9CA19BB86A4AF8DF53.gz",
        "tax_url": None,  # bundled in same archive
        "fasta_filename": "sh_refs_qiime_ver10_dynamic_dev.fasta",
        "tax_filename": "sh_taxonomy_qiime_ver10_dynamic_dev.txt",
        "seqs_qza": "unite-ITS-seqs.qza",
        "tax_qza": "unite-ITS-tax.qza",
        "classifier_qza": "unite-ITS-classifier.qza",
        "tax_format": "HeaderlessTSVTaxonomyFormat",
        "pretrained": False,
        "manual_download": True,  # UNITE requires manual download (see note)
    },
    "MiFish": {
        # Mitohelper pre-formatted MitoFish 12S reference database (Mar 2025).
        # Pre-built QIIME 2 QZA files hosted on Zenodo — no import step needed.
        # Replaces the previous MARES v2 entry — MARES is a COI database and was
        # never appropriate for MiFish 12S rRNA. Mitohelper extracts the 12S rRNA
        # region from the MitoFish mitogenome database and formats it for QIIME 2.
        # Source:   https://github.com/aomlomics/mitohelper
        # Zenodo:   https://zenodo.org/records/15028392
        # Citation: Lim et al. (2021) Environmental DNA. doi:10.1002/edn3.187
        "description": "MitoFish 12S rRNA via Mitohelper (Mar 2025, QIIME2-compatible)",
        "seqs_url":    "https://zenodo.org/records/15028392/files/12S-seqs-derep-uniq.qza",
        "tax_url":     "https://zenodo.org/records/15028392/files/12S-tax-derep-uniq.qza",
        "seqs_qza":    "mitofish-12S-seqs-derep-uniq.qza",
        "tax_qza":     "mitofish-12S-tax-derep-uniq.qza",
        "classifier_qza": "mitofish-12S-mifish-classifier.qza",
        "pretrained": False,
        # ── Gavia (loon host) augmentation ────────────────────────────────────
        # Without Gavia sequences in the reference, host reads from Common Loon
        # (Gavia immer) and related species fall through classification to
        # "uncl. Actinopteri", swamping the barplots (~81% of reads in loon data).
        # Use --add-gavia to fetch these from NCBI and add them before training.
        # Once in the reference, host reads are labeled as Aves/Gavia and can be
        # excluded with --exclude Gavia in 08_taxonomy_table.py.
        "gavia_classifier_qza": "mitofish-12S-mifish-gavia-classifier.qza",
        "gavia_merged_seqs_qza": "mitofish-gavia-merged-seqs.qza",
        "gavia_merged_tax_qza":  "mitofish-gavia-merged-tax.qza",
        "gavia_trimmed_qza":     "mitofish-gavia-seqs-mifish-trimmed.qza",
        "gavia_species": [
            # (NCBI organism name,  QIIME2 taxonomy string)
            ("Gavia immer",    "d__Eukaryota;p__Chordata;c__Aves;o__Gaviiformes;f__Gaviidae;g__Gavia;s__Gavia immer"),
            ("Gavia stellata", "d__Eukaryota;p__Chordata;c__Aves;o__Gaviiformes;f__Gaviidae;g__Gavia;s__Gavia stellata"),
            ("Gavia arctica",  "d__Eukaryota;p__Chordata;c__Aves;o__Gaviiformes;f__Gaviidae;g__Gavia;s__Gavia arctica"),
            ("Gavia pacifica", "d__Eukaryota;p__Chordata;c__Aves;o__Gaviiformes;f__Gaviidae;g__Gavia;s__Gavia pacifica"),
            ("Gavia adamsii",  "d__Eukaryota;p__Chordata;c__Aves;o__Gaviiformes;f__Gaviidae;g__Gavia;s__Gavia adamsii"),
        ],
    },
    "cytb": {
        "description": "NCBI vertebrate cytochrome b (fetched via RESCRIPt + Entrez)",
        # Output filenames — sequences are fetched at runtime via RESCRIPt, not wget
        "seqs_qza": "ncbi-cytb-vertebrata-seqs.qza",
        "tax_qza": "ncbi-cytb-vertebrata-tax.qza",
        "seqs_filtered_qza": "ncbi-cytb-vertebrata-seqs-filtered.qza",
        "tax_filtered_qza": "ncbi-cytb-vertebrata-tax-filtered.qza",
        "classifier_qza": "ncbi-cytb-vertebrata-classifier.qza",
        "pretrained": False,
        "requires_rescript": True,
        # Default NCBI Entrez query — targets cytochrome b gene in all vertebrates.
        # Uses the [Gene] field (not just [Title]) to capture sequences annotated
        # as "cytochrome b", "cytb", or "COB" regardless of record title wording.
        # Sequence length bounds (200–1200 bp) match the ~324 bp loon amplicon while
        # retaining full-length cytb (~1140 bp) references for good training coverage.
        # Whole mitogenomes are excluded by length (>16 kb) rather than by title
        # keyword, which avoids silently dropping short legitimate cytb partials.
        "entrez_query": (
            '"{taxon}"[Organism] AND '
            '("cytochrome b"[Gene] OR "cytb"[Gene] OR "COB"[Gene]) '
            'AND {min_len}:{max_len}[SLEN]'
        ),
        "entrez_db": "nuccore",
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def run(cmd: list[str], dry_run: bool = False) -> None:
    """Run a subprocess command with logging."""
    log.info("Running: %s", " ".join(cmd))
    if dry_run:
        log.info("[DRY RUN] — command not executed.")
        return
    result = subprocess.run(cmd, capture_output=False, text=True)
    if result.returncode != 0:
        log.error("Command failed with return code %d", result.returncode)
        sys.exit(result.returncode)


def wget(url: str, outpath: Path, dry_run: bool = False) -> None:
    """Download a file with wget."""
    if outpath.exists():
        log.info("File already exists, skipping download: %s", outpath)
        return
    run(["wget", "-O", str(outpath), url], dry_run=dry_run)


def gunzip(gz_path: Path, dry_run: bool = False) -> None:
    """Decompress a .gz file in place."""
    out_path = gz_path.with_suffix("")  # strip .gz
    if out_path.exists():
        log.info("Already decompressed: %s", out_path)
        return
    run(["gunzip", "-k", str(gz_path)], dry_run=dry_run)


def qiime_import_seqs(fasta: Path, qza: Path, dry_run: bool = False) -> None:
    if qza.exists():
        log.info("Already exists, skipping import: %s", qza)
        return
    run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Sequence]",
        "--input-path", str(fasta),
        "--output-path", str(qza),
    ], dry_run=dry_run)


def qiime_import_tax(tax: Path, qza: Path, tax_format: str, dry_run: bool = False) -> None:
    if qza.exists():
        log.info("Already exists, skipping import: %s", qza)
        return
    run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Taxonomy]",
        "--input-format", tax_format,
        "--input-path", str(tax),
        "--output-path", str(qza),
    ], dry_run=dry_run)


def qiime_extract_reads(
    seqs_qza: Path,
    trimmed_qza: Path,
    f_primer: str,
    r_primer: str,
    dry_run: bool = False,
) -> None:
    if trimmed_qza.exists():
        log.info("Already exists, skipping extraction: %s", trimmed_qza)
        return
    run([
        "qiime", "feature-classifier", "extract-reads",
        "--i-sequences", str(seqs_qza),
        "--p-f-primer", f_primer,
        "--p-r-primer", r_primer,
        "--p-min-length", "50",
        "--p-max-length", "1000",
        "--o-reads", str(trimmed_qza),
    ], dry_run=dry_run)


def qiime_train_classifier(
    seqs_qza: Path,
    tax_qza: Path,
    classifier_qza: Path,
    threads: int = 4,
    dry_run: bool = False,
) -> None:
    if classifier_qza.exists():
        log.info("Classifier already exists, skipping training: %s", classifier_qza)
        return
    run([
        "qiime", "feature-classifier", "fit-classifier-naive-bayes",
        "--i-reference-reads", str(seqs_qza),
        "--i-reference-taxonomy", str(tax_qza),
        "--o-classifier", str(classifier_qza),
    ], dry_run=dry_run)


def qiime_merge_seqs(qzas: list, merged_qza: Path, dry_run: bool = False) -> None:
    """Merge multiple FeatureData[Sequence] QZAs into one."""
    if merged_qza.exists():
        log.info("Already exists, skipping merge: %s", merged_qza)
        return
    cmd = ["qiime", "feature-table", "merge-seqs"]
    for q in qzas:
        cmd += ["--i-data", str(q)]
    cmd += ["--o-merged-data", str(merged_qza)]
    run(cmd, dry_run=dry_run)


def qiime_merge_taxa(qzas: list, merged_qza: Path, dry_run: bool = False) -> None:
    """Merge multiple FeatureData[Taxonomy] QZAs into one."""
    if merged_qza.exists():
        log.info("Already exists, skipping merge: %s", merged_qza)
        return
    cmd = ["qiime", "feature-table", "merge-taxa"]
    for q in qzas:
        cmd += ["--i-data", str(q)]
    cmd += ["--o-merged-data", str(merged_qza)]
    run(cmd, dry_run=dry_run)


def fetch_gavia_seqs_ncbi(
    outdir: Path,
    email: str,
    species_list: list,
    f_primer: str,
    r_primer_rc: str,
    max_per_species: int = 5,
    dry_run: bool = False,
) -> tuple:
    """
    Fetch Gavia (loon) 12S amplicon sequences from NCBI complete mitogenomes
    and write QIIME2-formatted FASTA + taxonomy files.

    Strategy
    --------
    MiFish-U primers do NOT bind to bird 12S with standard identity thresholds
    (too many mismatches between fish and bird primer sites).  QIIME2
    extract-reads and cutadapt at default error rates both return 0 sequences.

    The correct approach is:
    1. Fetch COMPLETE MITOGENOMES (16-17 kb) — these reliably exist for all
       Gavia species and contain the full 12S gene.
    2. Run cutadapt with HIGH error tolerance (-e 0.3, ~6 mismatches on 21bp)
       to extract the amplicon window even with bird-specific primer divergence.
    3. If cutadapt still fails (e.g. unusual primer orientation), fall back to
       HARDCODED reference amplicons confirmed by BLAST for each species.

    The output FASTA contains pre-trimmed amplicon sequences ready to be
    appended directly to the already-trimmed MitoFish QZA — no further
    primer extraction needed.

    Returns (fasta_path, taxonomy_path).
    """
    import tempfile as _tempfile

    # Hardcoded fallback amplicons confirmed by NCBI BLAST
    # These are MiFish-region sequences extracted from known Gavia mitogenomes
    HARDCODED_AMPLICONS = {
        # Gavia immer — confirmed BLAST hit (Gavia immer voucher 1B-105, etc.)
        # Sequence from rep-seqs of this study, BLAST identity >99% to G. immer
        "Gavia immer": "CACCGCGGTCACACAAGAGGCCCAAATTAACCGTATACACGGCGTAAAGAGTGGTACCATGCTATCCCATCAACTAGGATCAAAGTGCAACTGAGCTGTCGTAAGCCCAAGATGCATTAAAAGCCACCCTCAAGACGATCTTAGCACCCCCGATCAATTGAACCCCACGAAAGCTGGGACACAAACTGGGATTAGATAC",
    }

    fasta_path = outdir / "gavia-12S-seqs.fasta"
    tax_path   = outdir / "gavia-12S-taxonomy.tsv"

    if fasta_path.exists() and tax_path.exists():
        # Validate that existing file is non-empty
        n = sum(1 for l in open(fasta_path) if l.startswith(">"))
        if n > 0:
            log.info("Gavia sequences already present (%d seqs), skipping NCBI fetch.", n)
            log.info("  Delete %s to re-fetch.", fasta_path)
            return fasta_path, tax_path
        log.warning("Existing Gavia FASTA is empty — re-fetching.")
        fasta_path.unlink()
        tax_path.unlink()

    if dry_run:
        log.info("[DRY RUN] Would fetch Gavia mitogenomes from NCBI for: %s",
                 [sp for sp, _ in species_list])
        fasta_path.touch()
        tax_path.write_text("Feature ID\tTaxon\n")
        return fasta_path, tax_path

    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        log.error("Biopython is required for Gavia NCBI fetch: pip install biopython")
        sys.exit(1)

    if not email:
        log.error("--ncbi-email is required for Gavia NCBI fetch (NCBI policy).")
        sys.exit(1)

    Entrez.email = email
    total_written = 0
    tax_rows = ["Feature ID\tTaxon"]

    with open(fasta_path, "w") as fasta_out:
        for species_name, taxon_string in species_list:
            log.info("  Processing: %s", species_name)
            extracted_this_species = 0

            # ── Step 1: fetch complete mitogenomes ────────────────────────────
            handle = Entrez.esearch(
                db="nuccore",
                term=f'"{species_name}"[Organism] AND mitochondrion[Title] AND complete[Title]',
                retmax=max_per_species,
            )
            ids = Entrez.read(handle)["IdList"]
            handle.close()
            log.info("    Found %d complete mitogenomes: %s", len(ids), ids)

            for ncbi_id in ids:
                try:
                    fh = Entrez.efetch(db="nuccore", id=ncbi_id, rettype="fasta", retmode="text")
                    rec = SeqIO.read(fh, "fasta")
                    fh.close()
                    log.info("    %s: %d bp", ncbi_id, len(rec.seq))
                except Exception as exc:
                    log.warning("    Failed to fetch %s: %s", ncbi_id, exc)
                    continue

                # ── Step 2: cutadapt with high error tolerance ─────────────
                # -e 0.3 = 30% error rate ≈ 6 mismatches on a 21bp primer.
                # This is necessary because MiFish-U primers were designed for
                # fish and have 4-6 mismatches in the Gavia primer binding site.
                with _tempfile.NamedTemporaryFile(mode="w", suffix=".fasta",
                                                  delete=False) as tmp_in:
                    tmp_in.write(f">{ncbi_id}\n{str(rec.seq)}\n")
                    tmp_in_path = tmp_in.name
                tmp_out_path = tmp_in_path + "_amplicon.fasta"

                result = subprocess.run([
                    "cutadapt",
                    "-g", f_primer,
                    "-a", r_primer_rc,
                    "-e", "0.3",
                    "--discard-untrimmed",
                    "--minimum-length", "50",
                    "--maximum-length", "300",
                    "-o", tmp_out_path,
                    tmp_in_path,
                ], capture_output=True, text=True)
                os.unlink(tmp_in_path)

                amplicon_seq = None
                if os.path.exists(tmp_out_path):
                    content = open(tmp_out_path).read().strip()
                    if content:
                        lines = content.split("\n")
                        amplicon_seq = "".join(l for l in lines if not l.startswith(">"))
                    os.unlink(tmp_out_path)

                if amplicon_seq:
                    feat_id = f"Gavia_{ncbi_id}"
                    fasta_out.write(f">{feat_id}\n{amplicon_seq}\n")
                    tax_rows.append(f"{feat_id}\t{taxon_string}")
                    extracted_this_species += 1
                    total_written += 1
                    log.info("    ✓ amplicon extracted: %d bp", len(amplicon_seq))
                else:
                    log.debug("    ✗ cutadapt found no amplicon in %s (expected for some species)", ncbi_id)

            # ── Step 3: hardcoded fallback if cutadapt got nothing ────────────
            if extracted_this_species == 0:
                fallback_seq = HARDCODED_AMPLICONS.get(species_name)
                if fallback_seq:
                    feat_id = f"Gavia_hardcoded_{species_name.replace(' ', '_')}"
                    fasta_out.write(f">{feat_id}\n{fallback_seq}\n")
                    tax_rows.append(f"{feat_id}\t{taxon_string}")
                    total_written += 1
                    log.info(
                        "    ✓ used hardcoded fallback amplicon for %s "
                        "(BLAST-confirmed, cutadapt primer mismatch too high)", species_name
                    )
                else:
                    log.warning(
                        "    No amplicon extracted and no hardcoded fallback for %s. "
                        "Host reads from this species will remain as uncl. Actinopteri.",
                        species_name,
                    )

    with open(tax_path, "w") as tax_out:
        tax_out.write("\n".join(tax_rows) + "\n")

    log.info("Gavia fetch complete: %d amplicon sequences written to %s",
             total_written, fasta_path)

    if total_written == 0:
        log.error(
            "FATAL: No Gavia amplicon sequences obtained from NCBI or hardcoded fallbacks.\n"
            "You can manually place a pre-trimmed FASTA at:\n"
            "  %s\n"
            "with a matching taxonomy TSV at:\n"
            "  %s\n"
            "then re-run.", fasta_path, tax_path
        )
        sys.exit(1)

    return fasta_path, tax_path


# ---------------------------------------------------------------------------
# Per-marker build functions
# ---------------------------------------------------------------------------
def build_16S(outdir: Path, f_primer: str, r_primer: str, threads: int,
              dry_run: bool, use_v4: bool = False) -> None:
    cfg = DB_CONFIG["16S"]
    log.info("=== 16S: %s ===", cfg["description"])

    if use_v4:
        # ── Fast V4 path (recommended for V4 amplicons) ───────────────────────
        # Downloads pre-trimmed 515F/806R seqs from QIIME2 resources (~50 MB)
        # and trains a classifier on them.  This produces 40-60% genus
        # resolution on V4 reads vs <1% for the full-length classifier.
        # Training still takes ~2-4 hours but there is no in-silico PCR step.
        log.info("V4 mode — downloading pre-trimmed 515F/806R seqs and training.")
        log.info("  This will take 2-4 hours. Use tmux/screen or submit as a job.")
        v4_seqs_qza   = outdir / cfg["v4_seqs_filename"]
        tax_qza        = outdir / cfg["ref_tax_filename"]
        classifier_qza = outdir / cfg["v4_classifier_filename"]

        wget(cfg["v4_seqs_url"],  v4_seqs_qza, dry_run)
        wget(cfg["ref_tax_url"],  tax_qza,      dry_run)
        qiime_train_classifier(v4_seqs_qza, tax_qza, classifier_qza, threads, dry_run)
        log.info("V4 16S classifier saved to: %s", classifier_qza)
        log.info(
            "\nUpdate your classify step to use:\n"
            "  --i-classifier classifiers/%s", cfg["v4_classifier_filename"]
        )

    elif f_primer and r_primer:
        # ── Custom primer path ────────────────────────────────────────────────
        # Downloads full SILVA seqs, runs in-silico PCR (slow, ~6-8 hours),
        # then trains.  Use this only if your primers differ from 515F/806R.
        log.info("Custom primers provided — training region-specific 16S classifier.")
        log.info("  WARNING: Full SILVA extraction is slow (6-8+ hours).")
        log.info("  If you used 515F/806R, use --16s-v4 instead (2-4 hours).")
        seqs_qza       = outdir / cfg["ref_seqs_filename"]
        tax_qza        = outdir / cfg["ref_tax_filename"]
        trimmed_qza    = outdir / "silva-138-99-seqs-trimmed.qza"
        classifier_qza = outdir / "silva-138-99-classifier-trimmed.qza"

        wget(cfg["ref_seqs_url"], seqs_qza, dry_run)
        wget(cfg["ref_tax_url"],  tax_qza,  dry_run)
        qiime_extract_reads(seqs_qza, trimmed_qza, f_primer, r_primer, dry_run)
        qiime_train_classifier(trimmed_qza, tax_qza, classifier_qza, threads, dry_run)
        log.info("16S classifier saved to: %s", classifier_qza)

    else:
        # ── Pre-trained full-length fallback ──────────────────────────────────
        # Downloads the pre-trained full-length SILVA 138 classifier (209 MB).
        # DO NOT USE FOR V4 AMPLICONS — use --16s-v4 instead.
        log.warning(
            "Downloading full-length SILVA 138 classifier.\n"
            "  This will NOT give good genus resolution on V4 amplicons.\n"
            "  Re-run with --16s-v4 for V4 (515F/806R) data."
        )
        classifier_qza = outdir / cfg["pretrained_filename"]
        wget(cfg["pretrained_url"], classifier_qza, dry_run)
        log.info("16S classifier saved to: %s", classifier_qza)


def build_18S(outdir: Path, f_primer: str, r_primer: str, threads: int, dry_run: bool) -> None:
    cfg = DB_CONFIG["18S"]
    log.info("=== 18S: %s ===", cfg["description"])

    fasta_gz = outdir / (cfg["fasta_filename"] + ".gz")
    tax_gz = outdir / (cfg["tax_filename"] + ".gz")
    fasta = outdir / cfg["fasta_filename"]
    tax = outdir / cfg["tax_filename"]
    seqs_qza = outdir / cfg["seqs_qza"]
    tax_qza = outdir / cfg["tax_qza"]
    classifier_qza = outdir / cfg["classifier_qza"]

    wget(cfg["fasta_url"], fasta_gz, dry_run)
    wget(cfg["tax_url"], tax_gz, dry_run)
    gunzip(fasta_gz, dry_run)
    gunzip(tax_gz, dry_run)
    qiime_import_seqs(fasta, seqs_qza, dry_run)
    qiime_import_tax(tax, tax_qza, cfg["tax_format"], dry_run)

    if f_primer and r_primer:
        trimmed_qza = outdir / "pr2-18S-seqs-trimmed.qza"
        qiime_extract_reads(seqs_qza, trimmed_qza, f_primer, r_primer, dry_run)
        seqs_qza = trimmed_qza
        classifier_qza = outdir / "pr2-18S-classifier-trimmed.qza"

    qiime_train_classifier(seqs_qza, tax_qza, classifier_qza, threads, dry_run)
    log.info("18S classifier saved to: %s", classifier_qza)


def build_ITS(outdir: Path, f_primer: str, r_primer: str, threads: int, dry_run: bool) -> None:
    cfg = DB_CONFIG["ITS"]
    log.info("=== ITS: %s ===", cfg["description"])

    fasta = outdir / cfg["fasta_filename"]
    tax = outdir / cfg["tax_filename"]
    seqs_qza = outdir / cfg["seqs_qza"]
    tax_qza = outdir / cfg["tax_qza"]
    classifier_qza = outdir / cfg["classifier_qza"]

    # Check if the user has manually downloaded UNITE files
    if not fasta.exists() or not tax.exists():
        log.warning(
            "\n"
            "  *** UNITE requires a manual download ***\n"
            "  1. Go to: https://unite.ut.ee/repository.php\n"
            "  2. Download the QIIME release (developer version, all eukaryotes, dynamic)\n"
            "  3. Extract the archive and place these files in: %s\n"
            "       - %s\n"
            "       - %s\n"
            "  Then re-run this script.\n",
            outdir, cfg["fasta_filename"], cfg["tax_filename"],
        )
        return

    qiime_import_seqs(fasta, seqs_qza, dry_run)
    qiime_import_tax(tax, tax_qza, cfg["tax_format"], dry_run)

    if f_primer and r_primer:
        trimmed_qza = outdir / "unite-ITS-seqs-trimmed.qza"
        qiime_extract_reads(seqs_qza, trimmed_qza, f_primer, r_primer, dry_run)
        seqs_qza = trimmed_qza
        classifier_qza = outdir / "unite-ITS-classifier-trimmed.qza"

    qiime_train_classifier(seqs_qza, tax_qza, classifier_qza, threads, dry_run)
    log.info("ITS classifier saved to: %s", classifier_qza)


def build_MiFish(outdir: Path, f_primer: str, r_primer: str, threads: int,
                 dry_run: bool, add_gavia: bool = False,
                 ncbi_email: str = "") -> None:
    cfg = DB_CONFIG["MiFish"]
    log.info("=== MiFish 12S: %s ===", cfg["description"])

    # Default to MiFish-U primers (Miya et al. 2015) if none provided
    f_primer = f_primer or "GTCGGTAAAACTCGTGCCAGC"
    r_primer = r_primer or "CATAGTGGGGTATCTAATCCCAGTTTG"

    seqs_qza = outdir / cfg["seqs_qza"]
    tax_qza  = outdir / cfg["tax_qza"]

    # Download pre-built QIIME 2 QZAs directly from Zenodo — no import step needed
    wget(cfg["seqs_url"], seqs_qza, dry_run)
    wget(cfg["tax_url"],  tax_qza,  dry_run)

    if add_gavia:
        # ── Gavia-augmented classifier ────────────────────────────────────────
        # Adds Gavia (loon) 12S amplicon sequences to the already-trimmed
        # MitoFish reference so host reads are labeled as Aves/Gavia instead
        # of falling through to "uncl. Actinopteri".
        #
        # IMPORTANT: We do NOT run extract-reads on Gavia sequences.
        # MiFish-U primers have too many mismatches with bird 12S for
        # extract-reads to work at any reasonable identity threshold.
        # Instead, fetch_gavia_seqs_ncbi() extracts amplicons directly
        # from complete mitogenomes using cutadapt -e 0.3, and falls back
        # to hardcoded BLAST-confirmed amplicons if cutadapt also fails.
        # The resulting FASTA is pre-trimmed and appended directly to the
        # already-trimmed MitoFish sequences.
        log.info("Gavia augmentation enabled — fetching loon 12S amplicons...")

        # Reverse complement of reverse primer (for cutadapt -a)
        r_primer_rc = str(
            __import__("Bio.Seq", fromlist=["Seq"]).Seq(r_primer).reverse_complement()
        )

        gavia_fasta, gavia_tax = fetch_gavia_seqs_ncbi(
            outdir=outdir,
            email=ncbi_email,
            species_list=cfg["gavia_species"],
            f_primer=f_primer,
            r_primer_rc=r_primer_rc,
            dry_run=dry_run,
        )

        # Combine already-trimmed MitoFish FASTA with Gavia amplicons
        # We use the pre-trimmed QZA (not the full seqs) as the base.
        import zipfile as _zf, tempfile as _tf
        trimmed_mitofish_qza = outdir / "mitofish-12S-seqs-mifish-trimmed.qza"
        combined_fasta = outdir / "mitofish-gavia-combined-trimmed.fasta"

        if not trimmed_mitofish_qza.exists():
            log.error(
                "Cannot find already-trimmed MitoFish QZA at %s.\n"
                "Run without --add-gavia first to generate it, then re-run with --add-gavia.",
                trimmed_mitofish_qza,
            )
            sys.exit(1)

        log.info("Combining trimmed MitoFish seqs + Gavia amplicons...")
        with _zf.ZipFile(trimmed_mitofish_qza) as zf:
            fasta_member = next(n for n in zf.namelist() if n.endswith("dna-sequences.fasta"))
            mitofish_seqs = zf.read(fasta_member).decode()

        with open(combined_fasta, "w") as out:
            out.write(mitofish_seqs)
            if not mitofish_seqs.endswith("\n"):
                out.write("\n")
            out.write(open(gavia_fasta).read())

        n_mitofish = mitofish_seqs.count(">")
        n_gavia    = open(gavia_fasta).read().count(">")
        log.info("  MitoFish (trimmed): %d seqs", n_mitofish)
        log.info("  Gavia amplicons:    %d seqs", n_gavia)
        log.info("  Combined total:     %d seqs", n_mitofish + n_gavia)

        # Build combined taxonomy TSV
        combined_tax = outdir / "mitofish-gavia-combined-trimmed-tax.tsv"
        with _zf.ZipFile(outdir / "mitofish-12S-tax-derep-uniq.qza") as zf:
            tax_member = next(n for n in zf.namelist() if n.endswith("taxonomy.tsv"))
            mitofish_tax = zf.read(tax_member).decode().rstrip("\n")

        gavia_tax_lines = open(gavia_tax).read().strip().split("\n")[1:]  # skip header
        with open(combined_tax, "w") as out:
            out.write(mitofish_tax + "\n")
            out.write("\n".join(gavia_tax_lines) + "\n")

        # Import combined FASTA + taxonomy
        combined_seqs_qza = outdir / "mitofish-gavia-merged-trimmed-seqs.qza"
        combined_tax_qza  = outdir / "mitofish-gavia-merged-trimmed-tax.qza"
        qiime_import_seqs(combined_fasta, combined_seqs_qza, dry_run)
        qiime_import_tax(combined_tax, combined_tax_qza, "TSVTaxonomyFormat", dry_run)

        # Train directly on combined trimmed seqs — NO extract-reads needed
        classifier_qza = outdir / cfg["gavia_classifier_qza"]
        qiime_train_classifier(combined_seqs_qza, combined_tax_qza, classifier_qza, threads, dry_run)
        log.info("MiFish+Gavia classifier saved to: %s", classifier_qza)
        log.info(
            "\nUpdate your classify step to use:\n"
            "  --i-classifier classifiers/%s\n"
            "\nThe default MiFish exclude filter in 08_taxonomy_table.py\n"
            "will catch Gavia reads (Aves/Gaviidae in taxonomy string).",
            cfg["gavia_classifier_qza"],
        )

    else:
        # ── Standard MitoFish classifier (no host augmentation) ──────────────
        trimmed_qza    = outdir / "mitofish-12S-seqs-mifish-trimmed.qza"
        classifier_qza = outdir / cfg["classifier_qza"]

        log.info("Extracting MiFish amplicon region (primers: %s / %s)", f_primer, r_primer)
        if not trimmed_qza.exists():
            run([
                "qiime", "feature-classifier", "extract-reads",
                "--i-sequences",  str(seqs_qza),
                "--p-f-primer",   f_primer,
                "--p-r-primer",   r_primer,
                "--p-min-length", "100",
                "--p-max-length", "300",
                "--o-reads",      str(trimmed_qza),
            ], dry_run=dry_run)
        else:
            log.info("Already exists, skipping extraction: %s", trimmed_qza)

        qiime_train_classifier(trimmed_qza, tax_qza, classifier_qza, threads, dry_run)
        log.info("MiFish 12S classifier saved to: %s", classifier_qza)
        log.info(
            "NOTE: This classifier has NO Gavia sequences. Host reads will be\n"
            "classified as 'uncl. Actinopteri' and cannot be excluded cleanly.\n"
            "Re-run with --add-gavia --ncbi-email you@email.edu to fix this."
        )


# ---------------------------------------------------------------------------
# RESCRIPt helpers
# ---------------------------------------------------------------------------

def check_or_install_rescript(dry_run: bool = False) -> bool:
    """
    Check if the rescript QIIME 2 plugin is available.
    If not, attempt to install it via conda into the active environment.
    Returns True if rescript is available after the check/install attempt.
    """
    import shutil

    # Check if qiime is on PATH at all first
    if not shutil.which("qiime"):
        log.error("qiime not found on PATH. Activate your QIIME 2 conda environment first.")
        return False

    # Check if rescript plugin is registered
    result = subprocess.run(
        ["qiime", "rescript", "--help"],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        log.info("RESCRIPt plugin is available.")
        return True

    log.warning("RESCRIPt plugin not found. Attempting to install via conda...")
    if dry_run:
        log.info("[DRY RUN] Would run: conda install -c conda-forge -c bioconda -c qiime2 -c defaults xmltodict rescript")
        return True  # assume success in dry run

    install_cmd = [
        "conda", "install", "-y",
        "-c", "conda-forge",
        "-c", "bioconda",
        "-c", "qiime2",
        "-c", "defaults",
        "xmltodict", "rescript",
    ]
    log.info("Running: %s", " ".join(install_cmd))
    result = subprocess.run(install_cmd, text=True)
    if result.returncode != 0:
        log.error(
            "RESCRIPt installation failed.\n"
            "Install manually with:\n"
            "  conda install -c conda-forge -c bioconda -c qiime2 -c defaults xmltodict rescript\n"
            "Then re-run this script."
        )
        return False

    # Verify again
    result = subprocess.run(["qiime", "rescript", "--help"], capture_output=True, text=True)
    if result.returncode == 0:
        log.info("RESCRIPt installed successfully.")
        return True

    log.error("RESCRIPt installed but still not detected. Try reactivating your conda environment.")
    return False


def rescript_get_ncbi_data(
    query: str,
    outdir: Path,
    seqs_qza: Path,
    tax_qza: Path,
    ncbi_api_key: str,
    ncbi_email: str,
    entrez_db: str = "nuccore",
    dry_run: bool = False,
) -> None:
    """
    Use qiime rescript get-ncbi-data to fetch sequences and taxonomy from NCBI.
    Skips if both output QZAs already exist.
    """
    if seqs_qza.exists() and tax_qza.exists():
        log.info("NCBI data already downloaded — skipping fetch.")
        log.info("  Seqs: %s", seqs_qza)
        log.info("  Tax:  %s", tax_qza)
        return

    cmd = [
        "qiime", "rescript", "get-ncbi-data",
        "--p-query", query,
        "--p-n-jobs", "1",
        "--o-sequences", str(seqs_qza),
        "--o-taxonomy", str(tax_qza),
    ]

    # Set NCBI credentials as environment variables (safer than CLI args)
    env = os.environ.copy()
    if ncbi_api_key:
        env["NCBI_API_KEY"] = ncbi_api_key
    if ncbi_email:
        env["ENTREZ_EMAIL"] = ncbi_email

    log.info("Fetching NCBI cytb sequences — this may take 10–30 minutes...")
    log.info("Query: %s", query)

    if dry_run:
        log.info("[DRY RUN] Would run: %s", " ".join(cmd))
        return

    log.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, text=True, env=env)
    if result.returncode != 0:
        log.error("NCBI data fetch failed with return code %d", result.returncode)
        sys.exit(result.returncode)


def rescript_filter_seqs(
    seqs_qza: Path,
    tax_qza: Path,
    seqs_filtered_qza: Path,
    tax_filtered_qza: Path,
    min_len: int = 200,
    max_len: int = 1200,
    dry_run: bool = False,
) -> None:
    """
    Filter sequences by length, remove sequences with excess ambiguous bases,
    then dereplicate to remove redundant sequences.
    Skips if both final filtered outputs already exist.

    Intermediate files:
      ncbi-cytb-len-filtered.qza  — after length filter
      ncbi-cytb-culled.qza        — after cull-seqs (N-content + homopolymer)
    Final outputs (passed in as args):
      seqs_filtered_qza           — dereplicated sequences
      tax_filtered_qza            — dereplicated taxonomy
    """
    if seqs_filtered_qza.exists() and tax_filtered_qza.exists():
        log.info("Filtered sequences already exist — skipping filter step.")
        return

    outdir = seqs_filtered_qza.parent

    # ── Step 1: Length filter ──────────────────────────────────────────────
    len_filtered_qza     = outdir / "ncbi-cytb-len-filtered.qza"
    len_discarded_qza    = outdir / "ncbi-cytb-len-discarded.qza"

    if not len_filtered_qza.exists():
        run([
            "qiime", "rescript", "filter-seqs-length-by-taxon",
            "--i-sequences",      str(seqs_qza),
            "--i-taxonomy",       str(tax_qza),
            "--p-labels",         "Vertebrata",
            "--p-min-lens",       str(min_len),
            "--p-max-lens",       str(max_len),
            "--o-filtered-seqs",  str(len_filtered_qza),
            "--o-discarded-seqs", str(len_discarded_qza),
        ], dry_run=dry_run)
    else:
        log.info("Length-filtered sequences already exist — skipping.")

    # ── Step 2: Cull sequences with excess ambiguous bases / homopolymers ─
    # Output must be a DIFFERENT path from the final seqs_filtered_qza
    # because dereplicate will read it as input and write seqs_filtered_qza
    # as output.  Writing cull-seqs output directly to seqs_filtered_qza and
    # then passing it as both --i-sequences and --o-dereplicated-sequences to
    # dereplicate causes QIIME2 to refuse the overwrite (file already exists).
    culled_qza = outdir / "ncbi-cytb-culled.qza"

    if not culled_qza.exists():
        run([
            "qiime", "rescript", "cull-seqs",
            "--i-sequences",      str(len_filtered_qza),
            "--p-num-degenerates",   "5",
            "--p-homopolymer-length","8",
            "--o-clean-sequences",   str(culled_qza),
        ], dry_run=dry_run)
    else:
        log.info("Culled sequences already exist — skipping cull step.")

    # ── Step 3: Dereplicate ───────────────────────────────────────────────
    # Input:  culled_qza (distinct from seqs_filtered_qza)
    # Output: seqs_filtered_qza + tax_filtered_qza  (final products)
    run([
        "qiime", "rescript", "dereplicate",
        "--i-sequences",              str(culled_qza),
        "--i-taxa",                   str(tax_qza),
        "--p-mode",                   "uniq",
        "--o-dereplicated-sequences", str(seqs_filtered_qza),
        "--o-dereplicated-taxa",      str(tax_filtered_qza),
    ], dry_run=dry_run)


# ---------------------------------------------------------------------------
# cytb build function
# ---------------------------------------------------------------------------

def build_cytb(
    outdir: Path,
    f_primer: str,
    r_primer: str,
    threads: int,
    dry_run: bool,
    ncbi_api_key: str = "",
    ncbi_email: str = "",
    cytb_taxon: str = "Vertebrata",
    cytb_min_len: int = 200,
    cytb_max_len: int = 1200,
) -> None:
    """
    Build a QIIME 2 naive Bayes classifier for cytochrome b (cytb) using
    vertebrate sequences fetched from NCBI via the RESCRIPt plugin.

    Default primers are L14841 / H15149 (Kocher et al. universal vertebrate
    cytb), which match the ~324 bp amplicon detected in the loon GI data.
    These are used for read extraction before training if no custom primers
    are provided via --f-primer / --r-primer.

    IMPORTANT — Environment requirements
    -------------------------------------
    This function requires the RESCRIPt QIIME 2 plugin. The shared
    qiime2-amplicon-2024.5 environment on the cluster cannot be modified,
    so RESCRIPt must be installed in a SEPARATE conda environment:

        conda create -n qiime2-rescript --clone qiime2-amplicon-2024.5
        conda activate qiime2-rescript
        conda install -c conda-forge -c bioconda -c qiime2 -c defaults \\
            xmltodict rescript

    Or use a standalone QIIME2 + RESCRIPt env:

        conda create -n qiime2-rescript -c conda-forge -c bioconda \\
            -c qiime2 -c defaults qiime2 rescript

    Run this script from that env. The resulting classifier QZA is
    env-agnostic and can be used with qiime2-amplicon-2024.5 for all
    downstream steps (classify-sklearn, filter-table, etc.).

    NCBI fetch takes 10–30 minutes depending on connection speed and the
    number of vertebrate cytb sequences in GenBank (~500k records with
    default query). All intermediate files are cached so re-runs skip
    completed steps.
    """
    cfg = DB_CONFIG["cytb"]
    log.info("=== cytb: %s ===", cfg["description"])

    # L14841 / H15149 — Kocher et al. universal vertebrate cytb
    # These match the amplicon sequence data from the loon genome center output
    DEFAULT_F_PRIMER = "CGAAGCTTGATATGAAAAACCATCGTTG"   # L14841
    DEFAULT_R_PRIMER = "GGAAACAGCTATGACATTGATGGYGGTTTCG"  # H15149

    f = f_primer or DEFAULT_F_PRIMER
    r = r_primer or DEFAULT_R_PRIMER

    if not f_primer:
        log.info("No primers provided — using default L14841 / H15149 (Kocher universal vertebrate cytb).")
        log.info("  Forward (L14841): %s", f)
        log.info("  Reverse (H15149): %s", r)

    # ------------------------------------------------------------------
    # Step 1: Ensure RESCRIPt is installed
    # ------------------------------------------------------------------
    if not check_or_install_rescript(dry_run=dry_run):
        log.error("Cannot proceed without RESCRIPt. Exiting.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 2: Validate NCBI credentials
    # ------------------------------------------------------------------
    if not ncbi_api_key and not dry_run:
        log.warning(
            "No NCBI API key provided (--ncbi-api-key). Downloads will be rate-limited "
            "to 3 requests/second and may fail for large queries.\n"
            "  Get a free API key at: https://www.ncbi.nlm.nih.gov/account/\n"
            "  Then re-run with: --ncbi-api-key YOUR_KEY"
        )
    if not ncbi_email and not dry_run:
        log.warning(
            "No NCBI email provided (--ncbi-email). NCBI requires an email for Entrez access.\n"
            "  Add: --ncbi-email your@email.edu"
        )

    # ------------------------------------------------------------------
    # Step 3: Fetch sequences from NCBI via RESCRIPt
    # ------------------------------------------------------------------
    seqs_qza = outdir / cfg["seqs_qza"]
    tax_qza = outdir / cfg["tax_qza"]
    seqs_filtered_qza = outdir / cfg["seqs_filtered_qza"]
    tax_filtered_qza = outdir / cfg["tax_filtered_qza"]
    classifier_qza = outdir / cfg["classifier_qza"]

    # Format the Entrez query with runtime parameters
    query = cfg["entrez_query"].format(
        taxon=cytb_taxon,
        min_len=cytb_min_len,
        max_len=cytb_max_len,
    )

    rescript_get_ncbi_data(
        query=query,
        outdir=outdir,
        seqs_qza=seqs_qza,
        tax_qza=tax_qza,
        ncbi_api_key=ncbi_api_key,
        ncbi_email=ncbi_email,
        dry_run=dry_run,
    )

    # ------------------------------------------------------------------
    # Step 4: Filter sequences (length + ambiguous bases + dereplicate)
    # ------------------------------------------------------------------
    rescript_filter_seqs(
        seqs_qza=seqs_qza,
        tax_qza=tax_qza,
        seqs_filtered_qza=seqs_filtered_qza,
        tax_filtered_qza=tax_filtered_qza,
        min_len=cytb_min_len,
        max_len=cytb_max_len,
        dry_run=dry_run,
    )

    # ------------------------------------------------------------------
    # Step 5: Skip extract-reads for cytb — train on filtered seqs directly
    # ------------------------------------------------------------------
    # NCBI cytb records fetched by the [Gene] query are already gene fragments
    # (200-1200 bp), not full mitogenomes, so primer extraction adds no value.
    # More importantly, extract-reads returns "No matches found" against a broad
    # Vertebrata reference because Kocher L14841/H15149 primers vary too much
    # across vertebrate classes. Skipping is the correct approach — naive Bayes
    # classifies amplicon reads accurately from k-mer frequencies regardless.
    train_seqs_qza = seqs_filtered_qza
    log.info("Skipping extract-reads for cytb — training on filtered sequences directly.")

    # ------------------------------------------------------------------
    # Step 6: Train the naive Bayes classifier
    # ------------------------------------------------------------------
    qiime_train_classifier(
        seqs_qza=train_seqs_qza,
        tax_qza=tax_filtered_qza,
        classifier_qza=classifier_qza,
        threads=threads,
        dry_run=dry_run,
    )

    log.info("cytb classifier saved to: %s", classifier_qza)
    log.info(
        "\nPost-classification tip: cytb primers amplify some bacterial DNA in GI samples.\n"
        "After classifying, filter your feature table to retain only Vertebrata:\n"
        "  qiime taxa filter-table \\\n"
        "    --i-table qiime2/dada2/table_cytb.qza \\\n"
        "    --i-taxonomy qiime2/taxonomy/taxonomy_cytb.qza \\\n"
        "    --p-include Vertebrata \\\n"
        "    --p-exclude Bacteria,Viruses,Archaea \\\n"
        "    --o-filtered-table qiime2/dada2/table_cytb_vertebrata.qza"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
BUILDERS = {
    "16S": build_16S,
    "18S": build_18S,
    "ITS": build_ITS,
    "MiFish": build_MiFish,
    "cytb": build_cytb,
}

MARKER_ALIASES = {
    "its": "ITS",
    "its1": "ITS",
    "its2": "ITS",
    "its1-2": "ITS",
    "its1/2": "ITS",
    "mifish": "MiFish",
    "12s": "MiFish",
    "16s": "16S",
    "18s": "18S",
    "cytb": "cytb",
    "cytochrome b": "cytb",
    "cob": "cytb",
}


def resolve_marker(marker: str) -> str:
    """Normalize marker names to canonical keys."""
    resolved = MARKER_ALIASES.get(marker.lower(), marker)
    if resolved not in BUILDERS:
        log.error("Unknown marker: '%s'. Valid options: %s", marker, list(BUILDERS.keys()))
        sys.exit(1)
    return resolved


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build QIIME 2 classifiers for metabarcoding markers.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--markers", nargs="+", required=True,
        metavar="MARKER",
        help="Markers to build: 16S, 18S, ITS, MiFish (case-insensitive)",
    )
    parser.add_argument(
        "--outdir", default="classifiers",
        help="Output directory for classifiers and intermediate files (default: ./classifiers)",
    )
    parser.add_argument(
        "--f-primer", default=None,
        help="Forward primer sequence (optional). If provided with --r-primer, trains a region-specific classifier.",
    )
    parser.add_argument(
        "--r-primer", default=None,
        help="Reverse primer sequence (optional).",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Number of threads for classifier training (default: 4)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print commands without executing them.",
    )
    # 16S-specific
    s16 = parser.add_argument_group("16S options")
    s16.add_argument(
        "--16s-v4", action="store_true", dest="use_v4",
        help=(
            "Train a V4-region (515F/806R) SILVA 138 classifier instead of using "
            "the pre-trained full-length classifier. Strongly recommended if your "
            "16S data is V4 amplicons — improves genus resolution from <1%% to 40-60%%."
            " Takes 2-4 hours. Downloads silva-138-99-seqs-515-806.qza (~50 MB)."
        ),
    )
    # MiFish-specific
    mf = parser.add_argument_group("MiFish options")
    mf.add_argument(
        "--add-gavia", action="store_true",
        help=(
            "Add Gavia (loon) 12S sequences to the MitoFish reference before training. "
            "Required for loon diet studies — without Gavia in the DB, host reads fall "
            "through to 'uncl. Actinopteri' and swamp the barplots. "
            "Requires --ncbi-email."
        ),
    )
    # Shared NCBI args (used by both cytb and MiFish --add-gavia)
    ncbi = parser.add_argument_group("NCBI options (cytb and MiFish --add-gavia)")
    ncbi.add_argument(
        "--ncbi-api-key", default="",
        metavar="KEY",
        help=(
            "NCBI Entrez API key. Required for reliable cytb downloads. "
            "Get one free at https://www.ncbi.nlm.nih.gov/account/"
        ),
    )
    ncbi.add_argument(
        "--ncbi-email", default="",
        metavar="EMAIL",
        help="Email address for NCBI Entrez access (required for cytb and --add-gavia).",
    )
    # cytb-specific
    cytb = parser.add_argument_group("cytb options")
    cytb.add_argument(
        "--cytb-taxon", default="Vertebrata",
        metavar="TAXON",
        help="NCBI organism filter for cytb query. Default: Vertebrata",
    )
    cytb.add_argument(
        "--cytb-min-len", type=int, default=200,
        metavar="BP",
        help="Minimum cytb sequence length to retain. Default: 200",
    )
    cytb.add_argument(
        "--cytb-max-len", type=int, default=1200,
        metavar="BP",
        help="Maximum cytb sequence length to retain. Default: 1200",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Validate primer args
    if bool(args.f_primer) != bool(args.r_primer):
        log.error("Please provide both --f-primer and --r-primer, or neither.")
        sys.exit(1)

    # Set up output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log.info("Output directory: %s", outdir.resolve())

    if args.dry_run:
        log.info("DRY RUN mode — no commands will be executed.")

    if args.f_primer:
        log.info("Primers provided — will train region-specific classifiers.")
        log.info("  Forward: %s", args.f_primer)
        log.info("  Reverse: %s", args.r_primer)
    else:
        log.info("No primers provided — using full-length references (or pre-trained where available).")

    # Resolve and build each marker
    markers = [resolve_marker(m) for m in args.markers]
    log.info("Markers to process: %s", markers)

    for marker in markers:
        try:
            if marker == "cytb":
                BUILDERS[marker](
                    outdir=outdir,
                    f_primer=args.f_primer,
                    r_primer=args.r_primer,
                    threads=args.threads,
                    dry_run=args.dry_run,
                    ncbi_api_key=args.ncbi_api_key,
                    ncbi_email=args.ncbi_email,
                    cytb_taxon=args.cytb_taxon,
                    cytb_min_len=args.cytb_min_len,
                    cytb_max_len=args.cytb_max_len,
                )
            elif marker == "16S":
                BUILDERS[marker](
                    outdir=outdir,
                    f_primer=args.f_primer,
                    r_primer=args.r_primer,
                    threads=args.threads,
                    dry_run=args.dry_run,
                    use_v4=args.use_v4,
                )
            elif marker == "MiFish":
                if args.add_gavia and not args.ncbi_email:
                    log.error("--add-gavia requires --ncbi-email your@institution.edu")
                    sys.exit(1)
                BUILDERS[marker](
                    outdir=outdir,
                    f_primer=args.f_primer,
                    r_primer=args.r_primer,
                    threads=args.threads,
                    dry_run=args.dry_run,
                    add_gavia=args.add_gavia,
                    ncbi_email=args.ncbi_email,
                )
            else:
                BUILDERS[marker](
                    outdir=outdir,
                    f_primer=args.f_primer,
                    r_primer=args.r_primer,
                    threads=args.threads,
                    dry_run=args.dry_run,
                )
        except Exception as e:
            log.error("Failed to build classifier for %s: %s", marker, e)
            sys.exit(1)

    log.info("=== All done! ===")
    log.info("Classifiers are in: %s", outdir.resolve())
    log.info(
        "\nTo classify your rep-seqs, run:\n"
        "  qiime feature-classifier classify-sklearn \\\n"
        "    --i-classifier classifiers/<classifier>.qza \\\n"
        "    --i-reads your-rep-seqs.qza \\\n"
        "    --o-classification your-taxonomy.qza\n"
    )


if __name__ == "__main__":
    main()
