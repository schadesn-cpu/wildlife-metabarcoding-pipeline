#!/usr/bin/env python3
"""
00_build_classifiers.py
=======================
QIIME 2 classifier training and reference database builder for metabarcoding
and viral detection workflows.

Supported markers
-----------------
  16S        SILVA 138 (V4-specific recommended; --16s-v4)
  18S        PR2 v5.0.0
  ITS        UNITE v10 (developer release; manual download required)
  MiFish     MitoFish 12S rRNA via Mitohelper (Mar 2025)
             Use --add-gavia for loon diet studies (augments with Gavia seqs)
  cytb       NCBI vertebrate cytochrome b (via RESCRIPt + Entrez)
             Requires a separate qiime2-rescript conda environment — see below
  COI        MIDORI2 UNIQ NUC COI (Leray amplicon by default)
             Large download (~2-4 GB); extract-reads takes 2-4 hr
  adenovirus NCBI avian adenovirus (hexon gene) — builds BLAST DB only
             Presence/absence detection via BLAST (not NB classifier)
  herpesvirus NCBI avian herpesvirus (DNA polymerase) — builds BLAST DB only
             Presence/absence detection via BLAST (not NB classifier)

Usage examples
--------------
  # Train all metabarcoding classifiers
  python 00_build_classifiers.py --markers 16S 18S MiFish cytb COI \\
      --outdir classifiers/ --ncbi-email you@unh.edu --ncbi-api-key YOUR_KEY

  # 16S V4 (recommended for 515F/806R data)
  python 00_build_classifiers.py --markers 16S --16s-v4 --outdir classifiers/

  # MiFish with Gavia augmentation (required for loon studies)
  python 00_build_classifiers.py --markers MiFish --add-gavia \\
      --ncbi-email you@unh.edu --outdir classifiers/

  # COI Leray amplicon
  python 00_build_classifiers.py --markers COI --outdir classifiers/

  # cytb (requires qiime2-rescript env — see cytb section)
  conda activate qiime2-rescript
  python 00_build_classifiers.py --markers cytb --outdir classifiers/ \\
      --ncbi-email you@unh.edu --ncbi-api-key YOUR_KEY

  # Viral BLAST reference databases
  python 00_build_classifiers.py --markers adenovirus herpesvirus \\
      --outdir references/ --ncbi-email you@unh.edu --ncbi-api-key YOUR_KEY

  # Dry run to preview all commands
  python 00_build_classifiers.py --markers 16S MiFish COI --dry-run

Viral detection note
--------------------
adenovirus and herpesvirus build BLAST nucleotide databases, NOT QIIME2
naive Bayes classifiers. Presence/absence is determined by running BLAST
in a downstream detection script (e.g. 10_detect_viruses.py):

  blastn -query rep-seqs.fasta -db references/adenovirus/adenovirus_blast_db \\
         -perc_identity 80 -qcov_hsp_perc 60 -out results.txt -outfmt 6

cytb environment note
---------------------
cytb requires the RESCRIPt QIIME 2 plugin to fetch sequences from NCBI.
The shared qiime2-amplicon-2024.5 environment on the cluster cannot be
modified. Create a separate environment:

  conda create -n qiime2-rescript --clone qiime2-amplicon-2024.5
  conda activate qiime2-rescript
  conda install -c conda-forge -c bioconda -c qiime2 -c defaults xmltodict rescript

The resulting classifier QZA is env-agnostic and works in qiime2-amplicon-2024.5.
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

    # ── 16S ──────────────────────────────────────────────────────────────────
    "16S": {
        "description": "SILVA 138 99% OTUs (16S rRNA)",
        # Full-length pre-trained (209 MB). Poor on V4 amplicons (<1% genus).
        "pretrained_url":      "https://data.qiime2.org/2024.5/common/silva-138-99-nb-classifier.qza",
        "pretrained_filename": "silva-138-99-nb-classifier.qza",
        # Full-length seqs + taxonomy (for custom primer training)
        "ref_seqs_url":      "https://data.qiime2.org/2024.5/common/silva-138-99-seqs.qza",
        "ref_seqs_filename": "silva-138-99-seqs.qza",
        "ref_tax_url":       "https://data.qiime2.org/2024.5/common/silva-138-99-tax.qza",
        "ref_tax_filename":  "silva-138-99-tax.qza",
        # V4 amplicon (515F/806R) pre-trimmed — RECOMMENDED for V4 data
        # ~50 MB. Gives 40-60% genus resolution vs <1% for full-length.
        "v4_seqs_url":          "https://data.qiime2.org/2024.5/common/silva-138-99-seqs-515-806.qza",
        "v4_seqs_filename":     "silva-138-99-seqs-515-806.qza",
        "v4_classifier_filename": "silva-138-99-nb-classifier-515-806.qza",
        "pretrained": True,
    },

    # ── 18S ──────────────────────────────────────────────────────────────────
    "18S": {
        "description": "PR2 v5.0.0 (18S rRNA)",
        "fasta_url": "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_QIIME.fasta.gz",
        "tax_url":   "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_QIIME.tax.gz",
        "fasta_filename": "pr2_version_5.0.0_SSU_QIIME.fasta",
        "tax_filename":   "pr2_version_5.0.0_SSU_QIIME.tax",
        "seqs_qza":       "pr2-18S-seqs.qza",
        "tax_qza":        "pr2-18S-tax.qza",
        "classifier_qza": "pr2-18S-classifier.qza",
        "tax_format":     "HeaderlessTSVTaxonomyFormat",
        "pretrained": False,
    },

    # ── ITS ───────────────────────────────────────────────────────────────────
    "ITS": {
        "description": "UNITE v10 developer release (ITS1-2, all eukaryotes)",
        # UNITE requires manual download — see build_ITS() warning
        "fasta_filename": "sh_refs_qiime_ver10_dynamic_dev.fasta",
        "tax_filename":   "sh_taxonomy_qiime_ver10_dynamic_dev.txt",
        "seqs_qza":       "unite-ITS-seqs.qza",
        "tax_qza":        "unite-ITS-tax.qza",
        "classifier_qza": "unite-ITS-classifier.qza",
        "tax_format":     "HeaderlessTSVTaxonomyFormat",
        "manual_download": True,
    },

    # ── MiFish ───────────────────────────────────────────────────────────────
    "MiFish": {
        "description": "MitoFish 12S rRNA via Mitohelper (Mar 2025, QIIME2-compatible)",
        # Pre-built QIIME2 QZAs from Zenodo — no import step needed.
        # Lim et al. (2021) Environmental DNA. doi:10.1002/edn3.187
        "seqs_url":    "https://zenodo.org/records/15028392/files/12S-seqs-derep-uniq.qza",
        "tax_url":     "https://zenodo.org/records/15028392/files/12S-tax-derep-uniq.qza",
        "seqs_qza":    "mitofish-12S-seqs-derep-uniq.qza",
        "tax_qza":     "mitofish-12S-tax-derep-uniq.qza",
        "classifier_qza":       "mitofish-12S-mifish-classifier.qza",
        "gavia_classifier_qza": "mitofish-12S-mifish-gavia-classifier.qza",
        "gavia_species": [
            # (NCBI organism name, QIIME2 taxonomy string)
            ("Gavia immer",    "d__Eukaryota;p__Chordata;c__Aves;o__Gaviiformes;f__Gaviidae;g__Gavia;s__Gavia immer"),
            ("Gavia stellata", "d__Eukaryota;p__Chordata;c__Aves;o__Gaviiformes;f__Gaviidae;g__Gavia;s__Gavia stellata"),
            ("Gavia arctica",  "d__Eukaryota;p__Chordata;c__Aves;o__Gaviiformes;f__Gaviidae;g__Gavia;s__Gavia arctica"),
            ("Gavia pacifica", "d__Eukaryota;p__Chordata;c__Aves;o__Gaviiformes;f__Gaviidae;g__Gavia;s__Gavia pacifica"),
            ("Gavia adamsii",  "d__Eukaryota;p__Chordata;c__Aves;o__Gaviiformes;f__Gaviidae;g__Gavia;s__Gavia adamsii"),
        ],
    },

    # ── cytb ──────────────────────────────────────────────────────────────────
    "cytb": {
        "description": "NCBI vertebrate cytochrome b (fetched via RESCRIPt + Entrez)",
        "seqs_qza":          "ncbi-cytb-vertebrata-seqs.qza",
        "tax_qza":           "ncbi-cytb-vertebrata-tax.qza",
        "seqs_filtered_qza": "ncbi-cytb-vertebrata-seqs-filtered.qza",
        "tax_filtered_qza":  "ncbi-cytb-vertebrata-tax-filtered.qza",
        "classifier_qza":    "ncbi-cytb-vertebrata-classifier.qza",
        "requires_rescript": True,
        # Entrez query: [Gene] field captures "cytochrome b", "cytb", "COB".
        # Length bounds (200-1200 bp) capture both fragments and full-length cytb
        # while excluding whole mitogenomes (>16 kb). extract-reads is skipped
        # because NCBI [Gene]-queried records are already gene fragments, not
        # complete mitogenomes, so primer extraction adds no value and fails
        # across broad vertebrate diversity.
        "entrez_query": (
            '"{taxon}"[Organism] AND '
            '("cytochrome b"[Gene] OR "cytb"[Gene] OR "COB"[Gene]) '
            'AND {min_len}:{max_len}[SLEN]'
        ),
        "entrez_db": "nuccore",
    },

    # ── COI ───────────────────────────────────────────────────────────────────
    "COI": {
        "description": "MIDORI2 UNIQ NUC COI (cytochrome c oxidase I)",
        # MIDORI2 (Machida et al. 2017, Sci. Data) — standard COI reference.
        # Full 7-rank Eukaryota taxonomy; updated with each GenBank release.
        # QIIME_sp/uniq = one sequence per unique taxon × species combination.
        # NOTE: MIDORI2 filenames use "CO1" (older notation); we call it COI.
        # Source: https://www.reference-midori.info/
        "url_base":     "https://www.reference-midori.info/download/Databases/GenBank{version}.0/QIIME_sp/uniq/",
        "fasta_gz_tpl": "MIDORI2_UNIQ_NUC_GB{version}_CO1_QIIME.fasta.gz",
        "tax_gz_tpl":   "MIDORI2_UNIQ_NUC_GB{version}_CO1_QIIME_tax.txt.gz",
        "tax_format":   "HeaderlessTSVTaxonomyFormat",
        # Default primers: Leray mlCOIintF / jgHCO2198 (~313 bp, Leray et al. 2013)
        # Most widely used COI primers for marine/estuarine Metazoa.
        # For freshwater invertebrates consider BF3/BR2 or fwhF2/fwhR2n.
        # Provide --f-primer / --r-primer to override.
        "default_f_primer":     "GGWACWGGWTGAACWGTWTAYCCYCC",   # mlCOIintF
        "default_r_primer":     "TAIACYTCIGGRTGICCRAARAAYCA",   # jgHCO2198
        "default_primer_label": "leray",
        "default_primer_name":  "Leray mlCOIintF/jgHCO2198 (Leray et al. 2013)",
        # In-silico PCR is STRONGLY RECOMMENDED for COI:
        # Full gene = ~1.5 kb; Illumina reads = 150-300 bp.
        # Training on the amplicon window is far more accurate than full-gene.
        "amplicon_min_len": 250,
        "amplicon_max_len": 380,
    },

    # ── Adenovirus ────────────────────────────────────────────────────────────
    "adenovirus": {
        "description": "Avian adenovirus (hexon gene) — BLAST database for presence/absence",
        # Aviadenovirus and Siadenovirus are the two genera found in birds.
        # Hexon is the most conserved and sequenced gene across Adenoviridae.
        # NCBI query fetches avian adenovirus hexon sequences (200-2000 bp).
        # Presence/absence: BLAST rep-seqs against this DB.
        "entrez_query": (
            '("Aviadenovirus"[Organism] OR "Siadenovirus"[Organism] OR '
            '"fowl adenovirus"[Organism] OR "duck adenovirus"[Organism] OR '
            '"turkey adenovirus"[Organism]) AND '
            '("hexon"[Gene] OR "hexon protein"[Title] OR "penton"[Gene]) '
            'AND 200:2000[SLEN]'
        ),
        "entrez_db":      "nuccore",
        "fasta_filename": "avian-adenovirus-hexon.fasta",
        "blastdb_name":   "adenovirus_blast_db",
        # BLAST detection thresholds (use in downstream 10_detect_viruses.py)
        "blast_perc_identity": 80,
        "blast_qcov":          60,
        "blast_evalue":        "1e-10",
        "note": (
            "Builds a BLAST nucleotide database only — no QIIME2 NB classifier.\n"
            "Presence/absence: run blastn against this DB in 10_detect_viruses.py.\n"
            "Recommend: perc_identity >= 80, qcov_hsp_perc >= 60, evalue <= 1e-10.\n"
            "Positive = any rep-seq hit above these thresholds in a sample.\n"
            "Avian adenoviruses documented in Gaviiformes (loons), Anseriformes,\n"
            "Pelecaniformes; hexon gene is the diagnostic target of choice."
        ),
    },

    # ── Herpesvirus ───────────────────────────────────────────────────────────
    "herpesvirus": {
        "description": "Avian herpesvirus (DNA polymerase) — BLAST database for presence/absence",
        # Aviherpesviridae includes Marek's disease virus (MDV/GaHV-2),
        # Duck plague virus / Anatid herpesvirus 1 (DPV/AHV-1),
        # Infectious laryngotracheitis virus (ILTV/GaHV-1),
        # Columbid herpesvirus 1, and Loon herpesvirus (documented but rare DB).
        # DNA polymerase (UL30) is the most conserved gene across Herpesviridae
        # and is the standard phylogenetic/diagnostic target.
        "entrez_query": (
            '("Aviherpesviridae"[Organism] OR "Mardivirus"[Organism] OR '
            '"Iltovirus"[Organism] OR "Anatid herpesvirus"[Organism] OR '
            '"Gallid herpesvirus"[Organism] OR "Columbid herpesvirus"[Organism] OR '
            '"Psittacid herpesvirus"[Organism]) AND '
            '("DNA polymerase"[Gene] OR "UL30"[Gene] OR '
            '"glycoprotein B"[Gene] OR "gB"[Gene]) '
            'AND 300:3500[SLEN]'
        ),
        "entrez_db":      "nuccore",
        "fasta_filename": "avian-herpesvirus-dnap.fasta",
        "blastdb_name":   "herpesvirus_blast_db",
        "blast_perc_identity": 75,
        "blast_qcov":          60,
        "blast_evalue":        "1e-10",
        "note": (
            "Builds a BLAST nucleotide database only — no QIIME2 NB classifier.\n"
            "Presence/absence: run blastn against this DB in 10_detect_viruses.py.\n"
            "Recommend: perc_identity >= 75 (herpesviruses diverge more than adenoviruses),\n"
            "qcov_hsp_perc >= 60, evalue <= 1e-10.\n"
            "Avian herpesviruses documented in Gaviiformes, Anseriformes, Galliformes.\n"
            "DNA polymerase (UL30) is the preferred target; glycoprotein B is an alt.\n"
            "NOTE: Gavia-specific herpesvirus sequences are sparse in GenBank.\n"
            "Any positive hit at relaxed thresholds should be BLAST-confirmed manually."
        ),
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def run(cmd: list[str], dry_run: bool = False) -> None:
    """Run a subprocess command with logging."""
    log.info("Running: %s", " ".join(str(c) for c in cmd))
    if dry_run:
        log.info("[DRY RUN] — command not executed.")
        return
    result = subprocess.run(cmd, capture_output=False, text=True)
    if result.returncode != 0:
        log.error("Command failed with return code %d", result.returncode)
        sys.exit(result.returncode)


def wget(url: str, outpath: Path, dry_run: bool = False) -> None:
    """Download url to outpath using wget. Skip if outpath already exists."""
    if outpath.exists():
        log.info("File already exists, skipping download: %s", outpath)
        return
    run(["wget", "-O", str(outpath), url], dry_run=dry_run)


def gunzip(gz_path: Path, dry_run: bool = False) -> None:
    """Decompress gz_path in-place using gunzip -k. Skip if output already exists."""
    out_path = gz_path.with_suffix("")
    if out_path.exists():
        log.info("Already decompressed: %s", out_path)
        return
    run(["gunzip", "-k", str(gz_path)], dry_run=dry_run)


def qiime_import_seqs(fasta: Path, qza: Path, dry_run: bool = False) -> None:
    """Import a FASTA file as a QIIME2 FeatureData[Sequence] artifact. Skip if qza exists."""
    if qza.exists():
        log.info("Already exists, skipping import: %s", qza)
        return
    run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Sequence]",
        "--input-path", str(fasta),
        "--output-path", str(qza),
    ], dry_run=dry_run)


def qiime_import_tax(tax: Path, qza: Path, tax_format: str,
                     dry_run: bool = False) -> None:
    """Import a taxonomy file as a QIIME2 FeatureData[Taxonomy] artifact. Skip if qza exists."""
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


def qiime_extract_reads(seqs_qza: Path, trimmed_qza: Path,
                        f_primer: str, r_primer: str,
                        min_len: int = 50, max_len: int = 1000,
                        threads: int = 1,
                        dry_run: bool = False) -> None:
    """
Extract amplicon reads from a reference sequence QZA using primer sequences.

    Wraps qiime feature-classifier extract-reads. Produces a trimmed QZA
    containing only sequences that span the primer pair within min/max length
    bounds. Skip if trimmed_qza already exists.
    """
    if trimmed_qza.exists():
        log.info("Already exists, skipping extraction: %s", trimmed_qza)
        return
    run([
        "qiime", "feature-classifier", "extract-reads",
        "--i-sequences", str(seqs_qza),
        "--p-f-primer", f_primer,
        "--p-r-primer", r_primer,
        "--p-min-length", str(min_len),
        "--p-max-length", str(max_len),
        "--p-n-jobs", str(threads),
        "--o-reads", str(trimmed_qza),
    ], dry_run=dry_run)


def qiime_train_classifier(seqs_qza: Path, tax_qza: Path,
                            classifier_qza: Path, threads: int = 4,
                            dry_run: bool = False) -> None:
    """
Train a QIIME2 naive Bayes taxonomic classifier.

    Wraps qiime feature-classifier fit-classifier-naive-bayes. The resulting
    classifier QZA can be used directly with the taxonomy subcommand in
    03_run_full_metabarcoding_pipeline.py. Skip if classifier_qza exists.
    """
    if classifier_qza.exists():
        log.info("Classifier already exists, skipping training: %s", classifier_qza)
        return
    run([
        "qiime", "feature-classifier", "fit-classifier-naive-bayes",
        "--i-reference-reads", str(seqs_qza),
        "--i-reference-taxonomy", str(tax_qza),
        "--o-classifier", str(classifier_qza),
    ], dry_run=dry_run)


def makeblastdb(fasta: Path, db_name: Path, dry_run: bool = False) -> None:
    """Build a BLAST nucleotide database from a FASTA file."""
    nhr = db_name.with_suffix(".nhr")
    if nhr.exists():
        log.info("BLAST DB already exists, skipping: %s", db_name)
        return
    run([
        "makeblastdb",
        "-in", str(fasta),
        "-dbtype", "nucl",
        "-out", str(db_name),
        "-parse_seqids",
    ], dry_run=dry_run)


def fetch_ncbi_seqs(query: str, outpath: Path, entrez_db: str,
                    ncbi_email: str, ncbi_api_key: str = "",
                    max_records: int = 10000,
                    dry_run: bool = False) -> None:
    """
    Fetch sequences from NCBI Entrez via Biopython and write to a FASTA file.
    Used for adenovirus and herpesvirus reference building where the record
    counts are small enough to avoid the need for RESCRIPt batch fetch.
    Skips if outpath already exists and is non-empty.
    """
    if outpath.exists() and outpath.stat().st_size > 0:
        n = sum(1 for l in open(outpath) if l.startswith(">"))
        log.info("FASTA already exists (%d sequences) — skipping: %s", n, outpath)
        return

    if dry_run:
        log.info("[DRY RUN] Would fetch NCBI sequences to %s", outpath)
        log.info("  Query: %s", query)
        outpath.touch()
        return

    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        log.error("Biopython is required: pip install biopython --break-system-packages")
        sys.exit(1)

    if not ncbi_email:
        log.error("--ncbi-email is required for NCBI Entrez access (NCBI policy).")
        sys.exit(1)

    Entrez.email = ncbi_email
    if ncbi_api_key:
        Entrez.api_key = ncbi_api_key

    log.info("Searching NCBI %s — query: %s", entrez_db, query)
    handle = Entrez.esearch(db=entrez_db, term=query, retmax=max_records)
    result = Entrez.read(handle)
    handle.close()
    ids = result["IdList"]
    log.info("Found %d records. Fetching sequences...", len(ids))

    if not ids:
        log.error(
            "No sequences found for query: %s\n"
            "  Check the query in DB_CONFIG or broaden the search terms.", query
        )
        sys.exit(1)

    # Fetch in batches of 500 to avoid NCBI timeouts
    batch_size = 500
    total = 0
    with open(outpath, "w") as out_fh:
        for i in range(0, len(ids), batch_size):
            batch = ids[i:i + batch_size]
            log.info("  Fetching records %d–%d of %d...", i + 1,
                     min(i + batch_size, len(ids)), len(ids))
            try:
                fh = Entrez.efetch(
                    db=entrez_db,
                    id=",".join(batch),
                    rettype="fasta",
                    retmode="text",
                )
                for rec in SeqIO.parse(fh, "fasta"):
                    out_fh.write(f">{rec.id} {rec.description}\n{str(rec.seq)}\n")
                    total += 1
                fh.close()
            except Exception as exc:
                log.warning("  Batch %d failed: %s — continuing.", i, exc)

    log.info("Fetched %d sequences to %s", total, outpath)
    if total == 0:
        log.error("No sequences were written. Check NCBI query and credentials.")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Gavia 12S amplicon fetcher (MiFish augmentation)
# ---------------------------------------------------------------------------
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

    MiFish-U primers do NOT bind to bird 12S at standard identity thresholds
    (4-6 mismatches between fish and bird primer sites). The correct approach:
    1. Fetch COMPLETE MITOGENOMES (16-17 kb) — reliable for all Gavia species.
    2. Run cutadapt with -e 0.3 (~6 mismatches) to extract the amplicon window.
    3. Fall back to HARDCODED BLAST-confirmed amplicons if cutadapt also fails.

    Returns (fasta_path, taxonomy_path).
    """
    import tempfile as _tempfile

    # Hardcoded fallback — BLAST-confirmed amplicon from G. immer
    HARDCODED_AMPLICONS = {
        "Gavia immer": (
            "CACCGCGGTCACACAAGAGGCCCAAATTAACCGTATACACGGCGTAAAGAGTGGTACCATGCTATCCC"
            "ATCAACTAGGATCAAAGTGCAACTGAGCTGTCGTAAGCCCAAGATGCATTAAAAGCCACCCTCAAGAC"
            "GATCTTAGCACCCCCGATCAATTGAACCCCACGAAAGCTGGGACACAAACTGGGATTAGATAC"
        ),
    }

    fasta_path = outdir / "gavia-12S-seqs.fasta"
    tax_path   = outdir / "gavia-12S-taxonomy.tsv"

    if fasta_path.exists() and tax_path.exists():
        n = sum(1 for l in open(fasta_path) if l.startswith(">"))
        if n > 0:
            log.info("Gavia sequences already present (%d seqs), skipping fetch.", n)
            return fasta_path, tax_path
        log.warning("Existing Gavia FASTA is empty — re-fetching.")
        fasta_path.unlink()
        tax_path.unlink()

    if dry_run:
        log.info("[DRY RUN] Would fetch Gavia mitogenomes for: %s",
                 [sp for sp, _ in species_list])
        fasta_path.touch()
        tax_path.write_text("Feature ID\tTaxon\n")
        return fasta_path, tax_path

    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        log.error("Biopython is required: pip install biopython --break-system-packages")
        sys.exit(1)

    if not email:
        log.error("--ncbi-email is required for Gavia NCBI fetch.")
        sys.exit(1)

    Entrez.email = email
    total_written = 0
    tax_rows = ["Feature ID\tTaxon"]

    with open(fasta_path, "w") as fasta_out:
        for species_name, taxon_string in species_list:
            log.info("  Processing: %s", species_name)
            extracted_this_species = 0

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
                    fh = Entrez.efetch(db="nuccore", id=ncbi_id,
                                       rettype="fasta", retmode="text")
                    rec = SeqIO.read(fh, "fasta")
                    fh.close()
                    log.info("    %s: %d bp", ncbi_id, len(rec.seq))
                except Exception as exc:
                    log.warning("    Failed to fetch %s: %s", ncbi_id, exc)
                    continue

                with _tempfile.NamedTemporaryFile(
                        mode="w", suffix=".fasta", delete=False) as tmp_in:
                    tmp_in.write(f">{ncbi_id}\n{str(rec.seq)}\n")
                    tmp_in_path = tmp_in.name
                tmp_out_path = tmp_in_path + "_amplicon.fasta"

                subprocess.run([
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

            if extracted_this_species == 0:
                fallback_seq = HARDCODED_AMPLICONS.get(species_name)
                if fallback_seq:
                    feat_id = f"Gavia_hardcoded_{species_name.replace(' ', '_')}"
                    fasta_out.write(f">{feat_id}\n{fallback_seq}\n")
                    tax_rows.append(f"{feat_id}\t{taxon_string}")
                    total_written += 1
                    log.info("    ✓ used hardcoded fallback amplicon for %s", species_name)
                else:
                    log.warning(
                        "    No amplicon for %s. "
                        "Host reads will remain as uncl. Actinopteri.", species_name)

    with open(tax_path, "w") as tax_out:
        tax_out.write("\n".join(tax_rows) + "\n")

    log.info("Gavia fetch complete: %d amplicon sequences written.", total_written)

    if total_written == 0:
        log.error(
            "FATAL: No Gavia amplicon sequences obtained.\n"
            "Manually place a pre-trimmed FASTA at: %s\n"
            "with matching taxonomy TSV at: %s\nthen re-run.", fasta_path, tax_path
        )
        sys.exit(1)

    return fasta_path, tax_path


# ---------------------------------------------------------------------------
# RESCRIPt helpers (cytb only)
# ---------------------------------------------------------------------------
def check_rescript_available() -> bool:
    """
    Check if the RESCRIPt QIIME 2 plugin is available.
    Does NOT attempt auto-install (unsafe on shared cluster envs).
    Returns True if available, exits with instructions if not.
    """
    import shutil
    if not shutil.which("qiime"):
        log.error("qiime not found on PATH. Activate your QIIME 2 conda environment.")
        return False

    result = subprocess.run(
        ["qiime", "rescript", "--help"],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        log.info("RESCRIPt plugin is available.")
        return True

    log.error(
        "RESCRIPt plugin not found in the current conda environment.\n"
        "\n"
        "cytb requires RESCRIPt, which cannot be installed into the shared\n"
        "qiime2-amplicon-2024.5 environment on the cluster.\n"
        "\n"
        "Create a separate environment and run from it:\n"
        "\n"
        "  conda create -n qiime2-rescript --clone qiime2-amplicon-2024.5\n"
        "  conda activate qiime2-rescript\n"
        "  conda install -c conda-forge -c bioconda -c qiime2 -c defaults \\\n"
        "      xmltodict rescript\n"
        "\n"
        "Then re-run this script from qiime2-rescript.\n"
        "The resulting classifier QZA is env-agnostic and works in\n"
        "qiime2-amplicon-2024.5 for all downstream steps."
    )
    return False


def rescript_get_ncbi_data(
    query: str, outdir: Path, seqs_qza: Path, tax_qza: Path,
    ncbi_api_key: str, ncbi_email: str, dry_run: bool = False,
) -> None:
    """
Fetch sequences and taxonomy from NCBI via the RESCRIPt QIIME2 plugin.

    Runs qiime rescript get-ncbi-data with the provided Entrez query string.
    Requires the RESCRIPt plugin (conda install rescript) and NCBI credentials
    (--ncbi-email, --ncbi-api-key). Skip if both output QZAs already exist.
    """
    if seqs_qza.exists() and tax_qza.exists():
        log.info("NCBI data already downloaded — skipping fetch.")
        return

    cmd = [
        "qiime", "rescript", "get-ncbi-data",
        "--p-query", query,
        "--p-n-jobs", "1",
        "--o-sequences", str(seqs_qza),
        "--o-taxonomy", str(tax_qza),
    ]
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
    result = subprocess.run(cmd, text=True, env=env)
    if result.returncode != 0:
        log.error("NCBI data fetch failed with return code %d", result.returncode)
        sys.exit(result.returncode)


def rescript_filter_seqs(
    seqs_qza: Path, tax_qza: Path,
    seqs_filtered_qza: Path, tax_filtered_qza: Path,
    min_len: int = 200, max_len: int = 1200,
    dry_run: bool = False,
) -> None:
    """
    Filter cytb sequences: length filter → cull (N/homopolymer) → dereplicate.
    Uses three distinct intermediate paths to avoid QIIME2 overwrite errors.
    """
    if seqs_filtered_qza.exists() and tax_filtered_qza.exists():
        log.info("Filtered sequences already exist — skipping filter step.")
        return

    outdir = seqs_filtered_qza.parent
    len_filtered_qza  = outdir / "ncbi-cytb-len-filtered.qza"
    len_discarded_qza = outdir / "ncbi-cytb-len-discarded.qza"
    culled_qza        = outdir / "ncbi-cytb-culled.qza"

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

    if not culled_qza.exists():
        run([
            "qiime", "rescript", "cull-seqs",
            "--i-sequences",         str(len_filtered_qza),
            "--p-num-degenerates",   "5",
            "--p-homopolymer-length","8",
            "--o-clean-sequences",   str(culled_qza),
        ], dry_run=dry_run)
    else:
        log.info("Culled sequences already exist — skipping.")

    # Input: culled_qza (DISTINCT from seqs_filtered_qza — avoids overwrite error)
    run([
        "qiime", "rescript", "dereplicate",
        "--i-sequences",              str(culled_qza),
        "--i-taxa",                   str(tax_qza),
        "--p-mode",                   "uniq",
        "--o-dereplicated-sequences", str(seqs_filtered_qza),
        "--o-dereplicated-taxa",      str(tax_filtered_qza),
    ], dry_run=dry_run)


# ---------------------------------------------------------------------------
# Per-marker build functions
# ---------------------------------------------------------------------------

# ── 16S ─────────────────────────────────────────────────────────────────────
def build_16S(outdir: Path, f_primer: str, r_primer: str, threads: int,
              dry_run: bool, use_v4: bool = False) -> None:
    """
    Build a SILVA 138 16S classifier.

    --16s-v4 (recommended): Downloads pre-trimmed 515F/806R sequences and
    trains a V4-specific classifier. Improves genus resolution from <1% to
    40-60% on V4 amplicon data. Takes 2-4 hours.

    Custom primers: Downloads full SILVA seqs, runs in-silico PCR, then
    trains. Use only if primers differ from 515F/806R (6-8+ hours).

    No flags: Downloads the full-length pre-trained classifier (fast but
    NOT suitable for V4 amplicons).
    """
    cfg = DB_CONFIG["16S"]
    log.info("=== 16S: %s ===", cfg["description"])

    if use_v4:
        log.info("V4 mode (515F/806R) — downloading pre-trimmed seqs and training.")
        log.info("  Estimated time: 2-4 hours. Run in tmux/screen or as a cluster job.")
        v4_seqs_qza   = outdir / cfg["v4_seqs_filename"]
        tax_qza       = outdir / cfg["ref_tax_filename"]
        classifier_qza = outdir / cfg["v4_classifier_filename"]
        wget(cfg["v4_seqs_url"], v4_seqs_qza, dry_run)
        wget(cfg["ref_tax_url"], tax_qza,      dry_run)
        qiime_train_classifier(v4_seqs_qza, tax_qza, classifier_qza, threads, dry_run)
        log.info("V4 16S classifier saved to: %s", classifier_qza)

    elif f_primer and r_primer:
        log.info("Custom primers — training region-specific classifier (6-8+ hours).")
        log.info("  If you used 515F/806R, use --16s-v4 instead (2-4 hours).")
        seqs_qza      = outdir / cfg["ref_seqs_filename"]
        tax_qza       = outdir / cfg["ref_tax_filename"]
        trimmed_qza   = outdir / "silva-138-99-seqs-trimmed.qza"
        classifier_qza = outdir / "silva-138-99-classifier-trimmed.qza"
        wget(cfg["ref_seqs_url"], seqs_qza, dry_run)
        wget(cfg["ref_tax_url"],  tax_qza,  dry_run)
        qiime_extract_reads(seqs_qza, trimmed_qza, f_primer, r_primer,
                            threads=threads, dry_run=dry_run)
        qiime_train_classifier(trimmed_qza, tax_qza, classifier_qza, threads, dry_run)
        log.info("16S classifier saved to: %s", classifier_qza)

    else:
        log.warning(
            "Downloading full-length SILVA 138 pre-trained classifier.\n"
            "  This gives <1%% genus resolution on V4 amplicons.\n"
            "  Re-run with --16s-v4 for 515F/806R data."
        )
        classifier_qza = outdir / cfg["pretrained_filename"]
        wget(cfg["pretrained_url"], classifier_qza, dry_run)
        log.info("16S classifier saved to: %s", classifier_qza)


# ── 18S ─────────────────────────────────────────────────────────────────────
def build_18S(outdir: Path, f_primer: str, r_primer: str, threads: int,
              dry_run: bool) -> None:
    """
Build a QIIME2 naive Bayes classifier for 18S rRNA from the PR2 database.

    Downloads the PR2 FASTA and taxonomy from the configured URL, imports them
    into QIIME2, extracts the amplicon region using the supplied primer pair,
    and trains a classifier. All intermediate and final files are written to
    outdir. Steps that have already completed are skipped automatically.
    """
    cfg = DB_CONFIG["18S"]
    log.info("=== 18S: %s ===", cfg["description"])

    fasta_gz = outdir / (cfg["fasta_filename"] + ".gz")
    tax_gz   = outdir / (cfg["tax_filename"] + ".gz")
    fasta    = outdir / cfg["fasta_filename"]
    tax      = outdir / cfg["tax_filename"]
    seqs_qza = outdir / cfg["seqs_qza"]
    tax_qza  = outdir / cfg["tax_qza"]
    classifier_qza = outdir / cfg["classifier_qza"]

    wget(cfg["fasta_url"], fasta_gz, dry_run)
    wget(cfg["tax_url"],   tax_gz,   dry_run)
    gunzip(fasta_gz, dry_run)
    gunzip(tax_gz,   dry_run)
    qiime_import_seqs(fasta, seqs_qza, dry_run)
    qiime_import_tax(tax, tax_qza, cfg["tax_format"], dry_run)

    if f_primer and r_primer:
        trimmed_qza    = outdir / "pr2-18S-seqs-trimmed.qza"
        classifier_qza = outdir / "pr2-18S-classifier-trimmed.qza"
        qiime_extract_reads(seqs_qza, trimmed_qza, f_primer, r_primer,
                            threads=threads, dry_run=dry_run)
        seqs_qza = trimmed_qza

    qiime_train_classifier(seqs_qza, tax_qza, classifier_qza, threads, dry_run)
    log.info("18S classifier saved to: %s", classifier_qza)


# ── ITS ──────────────────────────────────────────────────────────────────────
def build_ITS(outdir: Path, f_primer: str, r_primer: str, threads: int,
              dry_run: bool) -> None:
    """
Build a QIIME2 naive Bayes classifier for ITS from the UNITE database.

    Expects the UNITE FASTA and taxonomy files to already be present in outdir
    (UNITE requires a manual download from unite.ut.ee — see README). Imports,
    optionally extracts the amplicon region, and trains a classifier. ITS
    amplicons are variable-length so primer extraction is skipped when primers
    are absent; pass --f-primer / --r-primer to enable region-specific extraction.
    """
    cfg = DB_CONFIG["ITS"]
    log.info("=== ITS: %s ===", cfg["description"])

    fasta = outdir / cfg["fasta_filename"]
    tax   = outdir / cfg["tax_filename"]
    seqs_qza = outdir / cfg["seqs_qza"]
    tax_qza  = outdir / cfg["tax_qza"]
    classifier_qza = outdir / cfg["classifier_qza"]

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
        trimmed_qza    = outdir / "unite-ITS-seqs-trimmed.qza"
        classifier_qza = outdir / "unite-ITS-classifier-trimmed.qza"
        qiime_extract_reads(seqs_qza, trimmed_qza, f_primer, r_primer,
                            threads=threads, dry_run=dry_run)
        seqs_qza = trimmed_qza

    qiime_train_classifier(seqs_qza, tax_qza, classifier_qza, threads, dry_run)
    log.info("ITS classifier saved to: %s", classifier_qza)


# ── MiFish ───────────────────────────────────────────────────────────────────
def build_MiFish(outdir: Path, f_primer: str, r_primer: str, threads: int,
                 dry_run: bool, add_gavia: bool = False,
                 ncbi_email: str = "") -> None:
    """
    Build a MitoFish 12S MiFish classifier.

    --add-gavia (required for loon studies): Augments the MitoFish reference
    with Gavia (loon) 12S amplicon sequences fetched from NCBI. Without this,
    host loon reads fall through to "uncl. Actinopteri" (~81% of reads) and
    swamp the diet barplots. Requires --ncbi-email.

    Note: extract-reads is NOT run on Gavia sequences — MiFish-U primers
    have 4-6 mismatches with bird 12S. Pre-trimmed amplicons are appended
    directly to the already-trimmed MitoFish QZA.
    """
    cfg = DB_CONFIG["MiFish"]
    log.info("=== MiFish 12S: %s ===", cfg["description"])

    # Default to MiFish-U primers (Miya et al. 2015)
    f_primer = f_primer or "GTCGGTAAAACTCGTGCCAGC"
    r_primer = r_primer or "CATAGTGGGGTATCTAATCCCAGTTTG"

    seqs_qza = outdir / cfg["seqs_qza"]
    tax_qza  = outdir / cfg["tax_qza"]
    wget(cfg["seqs_url"], seqs_qza, dry_run)
    wget(cfg["tax_url"],  tax_qza,  dry_run)

    if add_gavia:
        log.info("Gavia augmentation enabled — fetching loon 12S amplicons...")

        r_primer_rc = str(
            __import__("Bio.Seq", fromlist=["Seq"]).Seq(r_primer).reverse_complement()
        )
        gavia_fasta, gavia_tax = fetch_gavia_seqs_ncbi(
            outdir=outdir, email=ncbi_email,
            species_list=cfg["gavia_species"],
            f_primer=f_primer, r_primer_rc=r_primer_rc,
            dry_run=dry_run,
        )

        import zipfile as _zf
        trimmed_mitofish_qza = outdir / "mitofish-12S-seqs-mifish-trimmed.qza"
        if not trimmed_mitofish_qza.exists():
            log.error(
                "Cannot find already-trimmed MitoFish QZA at %s.\n"
                "Run without --add-gavia first to generate it, then re-run.",
                trimmed_mitofish_qza,
            )
            sys.exit(1)

        combined_fasta = outdir / "mitofish-gavia-combined-trimmed.fasta"
        log.info("Combining trimmed MitoFish seqs + Gavia amplicons...")
        with _zf.ZipFile(trimmed_mitofish_qza) as zf:
            fasta_member = next(n for n in zf.namelist()
                                if n.endswith("dna-sequences.fasta"))
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

        combined_tax = outdir / "mitofish-gavia-combined-trimmed-tax.tsv"
        with _zf.ZipFile(outdir / cfg["tax_qza"]) as zf:
            tax_member = next(n for n in zf.namelist() if n.endswith("taxonomy.tsv"))
            mitofish_tax = zf.read(tax_member).decode().rstrip("\n")

        gavia_tax_lines = open(gavia_tax).read().strip().split("\n")[1:]
        with open(combined_tax, "w") as out:
            out.write(mitofish_tax + "\n")
            out.write("\n".join(gavia_tax_lines) + "\n")

        combined_seqs_qza = outdir / "mitofish-gavia-merged-trimmed-seqs.qza"
        combined_tax_qza  = outdir / "mitofish-gavia-merged-trimmed-tax.qza"
        qiime_import_seqs(combined_fasta, combined_seqs_qza, dry_run)
        qiime_import_tax(combined_tax, combined_tax_qza, "TSVTaxonomyFormat", dry_run)

        classifier_qza = outdir / cfg["gavia_classifier_qza"]
        qiime_train_classifier(combined_seqs_qza, combined_tax_qza,
                               classifier_qza, threads, dry_run)
        log.info("MiFish+Gavia classifier saved to: %s", classifier_qza)

    else:
        trimmed_qza    = outdir / "mitofish-12S-seqs-mifish-trimmed.qza"
        classifier_qza = outdir / cfg["classifier_qza"]

        log.info("Extracting MiFish amplicon region (primers: %s / %s)",
                 f_primer, r_primer)
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
        log.warning(
            "No Gavia sequences in this classifier. Host reads will be\n"
            "classified as 'uncl. Actinopteri' and cannot be excluded cleanly.\n"
            "Re-run with --add-gavia --ncbi-email you@email.edu to fix this."
        )


# ── cytb ─────────────────────────────────────────────────────────────────────
def build_cytb(
    outdir: Path, f_primer: str, r_primer: str, threads: int, dry_run: bool,
    ncbi_api_key: str = "", ncbi_email: str = "",
    cytb_taxon: str = "Vertebrata",
    cytb_min_len: int = 200, cytb_max_len: int = 1200,
) -> None:
    """
    Build a vertebrate cytb classifier via RESCRIPt + NCBI Entrez.

    Default primers: L14841 / H15149 (Kocher et al. universal vertebrate cytb).
    extract-reads is SKIPPED for cytb — NCBI [Gene]-queried records are already
    gene fragments (not full mitogenomes), so primer extraction adds no value
    and fails unreliably across broad vertebrate diversity.

    REQUIRES: qiime2-rescript conda environment (see module docstring).
    """
    cfg = DB_CONFIG["cytb"]
    log.info("=== cytb: %s ===", cfg["description"])

    # L14841 / H15149 — Kocher et al. universal vertebrate cytb
    DEFAULT_F = "CGAAGCTTGATATGAAAAACCATCGTTG"    # L14841
    DEFAULT_R = "GGAAACAGCTATGACATTGATGGYGGTTTCG"  # H15149
    f = f_primer or DEFAULT_F
    r = r_primer or DEFAULT_R

    if not f_primer:
        log.info("No primers provided — using default L14841/H15149 (Kocher cytb).")
        log.info("  Forward (L14841): %s", f)
        log.info("  Reverse (H15149): %s", r)

    if not check_rescript_available():
        sys.exit(1)

    if not ncbi_api_key:
        log.warning(
            "No --ncbi-api-key provided. Downloads rate-limited to 3 req/sec.\n"
            "  Get a free key at: https://www.ncbi.nlm.nih.gov/account/"
        )

    seqs_qza          = outdir / cfg["seqs_qza"]
    tax_qza           = outdir / cfg["tax_qza"]
    seqs_filtered_qza = outdir / cfg["seqs_filtered_qza"]
    tax_filtered_qza  = outdir / cfg["tax_filtered_qza"]
    classifier_qza    = outdir / cfg["classifier_qza"]

    query = cfg["entrez_query"].format(
        taxon=cytb_taxon, min_len=cytb_min_len, max_len=cytb_max_len
    )

    rescript_get_ncbi_data(
        query=query, outdir=outdir, seqs_qza=seqs_qza, tax_qza=tax_qza,
        ncbi_api_key=ncbi_api_key, ncbi_email=ncbi_email, dry_run=dry_run,
    )
    rescript_filter_seqs(
        seqs_qza=seqs_qza, tax_qza=tax_qza,
        seqs_filtered_qza=seqs_filtered_qza, tax_filtered_qza=tax_filtered_qza,
        min_len=cytb_min_len, max_len=cytb_max_len, dry_run=dry_run,
    )

    # Skipping extract-reads — see function docstring for rationale
    log.info("Skipping extract-reads for cytb — training on filtered sequences directly.")
    qiime_train_classifier(seqs_filtered_qza, tax_filtered_qza,
                           classifier_qza, threads, dry_run)
    log.info("cytb classifier saved to: %s", classifier_qza)
    log.info(
        "\nPost-classification tip: cytb primers amplify some bacterial DNA in GI samples.\n"
        "Filter your feature table to retain only Vertebrata after classifying:\n"
        "  qiime taxa filter-table \\\n"
        "    --i-table qiime2/dada2/table_cytb.qza \\\n"
        "    --i-taxonomy qiime2/taxonomy/taxonomy_cytb.qza \\\n"
        "    --p-include Vertebrata \\\n"
        "    --p-exclude Bacteria,Viruses,Archaea \\\n"
        "    --o-filtered-table qiime2/dada2/table_cytb_vertebrata.qza"
    )


# ── COI ──────────────────────────────────────────────────────────────────────
def build_COI(
    outdir: Path, f_primer: str, r_primer: str, threads: int, dry_run: bool,
    coi_version: str = "261", skip_extract: bool = False,
) -> None:
    """
    Build a MIDORI2 COI (cytochrome c oxidase I) classifier.

    Reference: MIDORI2 UNIQ NUC (Machida et al. 2017, Sci. Data).
    Full 7-rank Eukaryota taxonomy; updated with each GenBank release.

    Default primers: Leray mlCOIintF / jgHCO2198 (~313 bp, Leray et al. 2013).
    Most widely adopted for marine/estuarine Metazoa. For freshwater
    invertebrates consider BF3/BR2 or fwhF2/fwhR2n (use --f-primer/--r-primer).

    in-silico PCR is STRONGLY RECOMMENDED (do not use --coi-skip-extract unless
    you have a specific reason). Full COI gene is ~1.5 kb; training on the
    amplicon window dramatically improves accuracy on Illumina short reads.

    Timing (16-core cluster node):
      Download (~2-4 GB compressed) : 10-30 min
      gunzip                         :  5-10 min
      QIIME2 import                  : 15-30 min
      extract-reads (Leray)          :  2-4 hr
      classifier training            :  2-4 hr
    """
    cfg = DB_CONFIG["COI"]
    log.info("=== COI: %s (MIDORI2 GB%s) ===", cfg["description"], coi_version)

    url_base    = cfg["url_base"].format(version=coi_version)
    fasta_gz_fn = cfg["fasta_gz_tpl"].format(version=coi_version)
    tax_gz_fn   = cfg["tax_gz_tpl"].format(version=coi_version)

    fasta_gz = outdir / fasta_gz_fn
    tax_gz   = outdir / tax_gz_fn
    fasta    = outdir / fasta_gz_fn.replace(".gz", "")
    tax      = outdir / tax_gz_fn.replace(".gz", "")
    seqs_qza = outdir / f"midori2-COI-seqs-GB{coi_version}.qza"
    tax_qza  = outdir / f"midori2-COI-tax-GB{coi_version}.qza"

    log.info(
        "Downloading MIDORI2 COI reference (~2-4 GB compressed).\n"
        "  Use a stable connection or submit as a cluster job.\n"
        "  Download is skipped if files already exist."
    )
    wget(url_base + fasta_gz_fn, fasta_gz, dry_run)
    wget(url_base + tax_gz_fn,   tax_gz,   dry_run)
    gunzip(fasta_gz, dry_run)
    gunzip(tax_gz,   dry_run)
    qiime_import_seqs(fasta, seqs_qza, dry_run)
    qiime_import_tax(tax, tax_qza, cfg["tax_format"], dry_run)

    f = f_primer or cfg["default_f_primer"]
    r = r_primer or cfg["default_r_primer"]
    primer_label = "custom" if f_primer else cfg["default_primer_label"]

    if not f_primer:
        log.info("No primers provided — using default %s.", cfg["default_primer_name"])
        log.info("  Forward (mlCOIintF): %s", f)
        log.info("  Reverse (jgHCO2198): %s", r)

    if skip_extract:
        log.warning(
            "--coi-skip-extract: training on full COI gene (~1.5 kb).\n"
            "  NOT recommended for Illumina short-read data.\n"
            "  Amplicon-trained classifiers substantially outperform full-length."
        )
        train_seqs_qza = seqs_qza
        classifier_qza = outdir / f"midori2-COI-fulllength-classifier-GB{coi_version}.qza"
    else:
        trimmed_qza = outdir / f"midori2-COI-seqs-GB{coi_version}-{primer_label}-trimmed.qza"
        if not trimmed_qza.exists():
            log.info(
                "Running in-silico PCR with %s primers.\n"
                "  Takes 2-4 hours. Submit as a cluster job or use tmux/screen.",
                primer_label,
            )
            run([
                "qiime", "feature-classifier", "extract-reads",
                "--i-sequences",  str(seqs_qza),
                "--p-f-primer",   f,
                "--p-r-primer",   r,
                "--p-min-length", str(cfg["amplicon_min_len"]),
                "--p-max-length", str(cfg["amplicon_max_len"]),
                "--p-n-jobs",     str(threads),
                "--o-reads",      str(trimmed_qza),
            ], dry_run=dry_run)
        else:
            log.info("Already exists, skipping extraction: %s", trimmed_qza)

        train_seqs_qza = trimmed_qza
        classifier_qza = outdir / f"midori2-COI-{primer_label}-classifier-GB{coi_version}.qza"

    log.info("Training COI classifier — 2-4 hours. Use tmux/screen or cluster job.")
    qiime_train_classifier(train_seqs_qza, tax_qza, classifier_qza, threads, dry_run)
    log.info("COI classifier saved to: %s", classifier_qza)
    log.info(
        "\nPost-classification tip: filter to Metazoa before analysis:\n"
        "  qiime taxa filter-table \\\n"
        "    --i-table qiime2/COI/dada2/table.qza \\\n"
        "    --i-taxonomy qiime2/COI/taxonomy/taxonomy.qza \\\n"
        "    --p-include Metazoa \\\n"
        "    --p-exclude Bacteria,Viruses,Archaea,Viridiplantae,Fungi \\\n"
        "    --o-filtered-table qiime2/COI/dada2/table_metazoa.qza\n"
        "\nMIDORI2 taxonomy uses no rank prefixes;\n"
        "08_taxonomy_table.py handles this automatically."
    )


# ── Adenovirus ────────────────────────────────────────────────────────────────
def build_adenovirus(
    outdir: Path, dry_run: bool,
    ncbi_email: str = "", ncbi_api_key: str = "",
) -> None:
    """
    Build a BLAST nucleotide database for avian adenovirus detection.

    Target: hexon gene of Aviadenovirus and Siadenovirus (the two genera
    found in birds). Hexon is the most conserved and sequenced gene across
    Adenoviridae and is the standard diagnostic target.

    Output: BLAST database files in outdir/adenovirus/
    Detection: Run blastn in 10_detect_viruses.py:
      blastn -query rep-seqs.fasta -db adenovirus/adenovirus_blast_db \\
             -perc_identity 80 -qcov_hsp_perc 60 -evalue 1e-10 \\
             -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore"

    Presence threshold: any rep-seq with perc_identity >= 80 AND qcov >= 60
    is considered a positive. All positives should be manually BLASTed at
    NCBI to confirm species identity.
    """
    cfg = DB_CONFIG["adenovirus"]
    log.info("=== adenovirus: %s ===", cfg["description"])
    log.info(cfg["note"])

    virus_dir = outdir / "adenovirus"
    virus_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = virus_dir / cfg["fasta_filename"]
    db_path    = virus_dir / cfg["blastdb_name"]

    fetch_ncbi_seqs(
        query=cfg["entrez_query"],
        outpath=fasta_path,
        entrez_db=cfg["entrez_db"],
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        max_records=5000,
        dry_run=dry_run,
    )
    makeblastdb(fasta_path, db_path, dry_run)

    log.info("Adenovirus BLAST database saved to: %s", db_path)
    log.info(
        "\nDetection command (run in 10_detect_viruses.py or manually):\n"
        "  blastn -query your-rep-seqs.fasta \\\n"
        "         -db %s \\\n"
        "         -perc_identity %d -qcov_hsp_perc %d -evalue %s \\\n"
        "         -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' \\\n"
        "         -out adenovirus_hits.tsv\n"
        "\n"
        "A sample is POSITIVE if any rep-seq has a hit above the thresholds.\n"
        "Confirm all positives with manual NCBI BLAST.",
        db_path,
        cfg["blast_perc_identity"],
        cfg["blast_qcov"],
        cfg["blast_evalue"],
    )


# ── Herpesvirus ───────────────────────────────────────────────────────────────
def build_herpesvirus(
    outdir: Path, dry_run: bool,
    ncbi_email: str = "", ncbi_api_key: str = "",
) -> None:
    """
    Build a BLAST nucleotide database for avian herpesvirus detection.

    Targets: DNA polymerase (UL30) and glycoprotein B (gB) of Aviherpesviridae.
    UL30 is the most conserved gene across Herpesviridae and is the standard
    phylogenetic and diagnostic target. gB is an alternative target with
    more sequences available for some genera.

    Includes: Marek's disease virus (GaHV-2), Infectious laryngotracheitis
    virus (GaHV-1), Duck plague virus (Anatid herpesvirus 1), Columbid
    herpesvirus 1, and Psittacid herpesviruses.

    Note: Gavia-specific herpesvirus sequences are sparse in GenBank.
    Any positive hit should be BLAST-confirmed manually at NCBI.

    Output: BLAST database files in outdir/herpesvirus/
    Detection: Run blastn in 10_detect_viruses.py:
      blastn -query rep-seqs.fasta -db herpesvirus/herpesvirus_blast_db \\
             -perc_identity 75 -qcov_hsp_perc 60 -evalue 1e-10 \\
             -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore"

    Use perc_identity >= 75 (not 80) — herpesviruses diverge more than
    adenoviruses across host taxa and thresholds should be slightly relaxed.
    """
    cfg = DB_CONFIG["herpesvirus"]
    log.info("=== herpesvirus: %s ===", cfg["description"])
    log.info(cfg["note"])

    virus_dir = outdir / "herpesvirus"
    virus_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = virus_dir / cfg["fasta_filename"]
    db_path    = virus_dir / cfg["blastdb_name"]

    fetch_ncbi_seqs(
        query=cfg["entrez_query"],
        outpath=fasta_path,
        entrez_db=cfg["entrez_db"],
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        max_records=5000,
        dry_run=dry_run,
    )
    makeblastdb(fasta_path, db_path, dry_run)

    log.info("Herpesvirus BLAST database saved to: %s", db_path)
    log.info(
        "\nDetection command (run in 10_detect_viruses.py or manually):\n"
        "  blastn -query your-rep-seqs.fasta \\\n"
        "         -db %s \\\n"
        "         -perc_identity %d -qcov_hsp_perc %d -evalue %s \\\n"
        "         -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' \\\n"
        "         -out herpesvirus_hits.tsv\n"
        "\n"
        "A sample is POSITIVE if any rep-seq has a hit above the thresholds.\n"
        "Confirm all positives with manual NCBI BLAST.\n"
        "NOTE: Gavia-specific herpesvirus seqs are sparse in GenBank;\n"
        "relax thresholds to 70%% identity if initial screen finds nothing.",
        db_path,
        cfg["blast_perc_identity"],
        cfg["blast_qcov"],
        cfg["blast_evalue"],
    )


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
BUILDERS = {
    "16S":        build_16S,
    "18S":        build_18S,
    "ITS":        build_ITS,
    "MiFish":     build_MiFish,
    "cytb":       build_cytb,
    "COI":        build_COI,
    "adenovirus": build_adenovirus,
    "herpesvirus": build_herpesvirus,
}

MARKER_ALIASES = {
    # 16S
    "16s": "16S",
    # 18S
    "18s": "18S",
    # ITS
    "its": "ITS", "its1": "ITS", "its2": "ITS",
    "its1-2": "ITS", "its1/2": "ITS",
    # MiFish
    "mifish": "MiFish", "12s": "MiFish",
    # cytb
    "cytb": "cytb", "cytochrome b": "cytb", "cob": "cytb",
    # COI
    "coi": "COI", "co1": "COI", "cox1": "COI",
    "cox-1": "COI", "cytochrome oxidase i": "COI",
    "cytochrome c oxidase i": "COI",
    # Viruses
    "adeno": "adenovirus", "adv": "adenovirus",
    "herpes": "herpesvirus", "hv": "herpesvirus",
}


def resolve_marker(marker: str) -> str:
    """
Normalise a marker name via MARKER_ALIASES and validate it against BUILDERS.

    Allows common aliases (e.g. '12S' -> 'MiFish', 'cytochrome_b' -> 'cytb') so
    users are not required to know the exact internal key. Exits with an error
    message listing valid options if the marker is unrecognised after aliasing.
    """
    resolved = MARKER_ALIASES.get(marker.lower(), marker)
    if resolved not in BUILDERS:
        log.error(
            "Unknown marker: '%s'. Valid options: %s",
            marker, list(BUILDERS.keys())
        )
        sys.exit(1)
    return resolved


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Build and parse the command-line argument parser for 00_build_classifiers.py."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--markers", nargs="+", required=True, metavar="MARKER",
        help=(
            "Markers to build. Choose from: "
            "16S, 18S, ITS, MiFish, cytb, COI, adenovirus, herpesvirus "
            "(case-insensitive; aliases: adeno, herpes, co1, cox1, etc.)"
        ),
    )
    parser.add_argument(
        "--outdir", default="classifiers/", type=Path,
        help="Output directory for classifiers and reference databases (default: classifiers/)",
    )
    parser.add_argument(
        "--f-primer", default=None,
        help="Forward primer sequence (optional; uses marker-specific default if omitted).",
    )
    parser.add_argument(
        "--r-primer", default=None,
        help="Reverse primer sequence (optional; uses marker-specific default if omitted).",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Threads for classifier training (default: 4).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print commands without executing them.",
    )

    # 16S
    s16 = parser.add_argument_group("16S options")
    s16.add_argument(
        "--16s-v4", action="store_true", dest="use_v4",
        help=(
            "Train a V4-region (515F/806R) SILVA 138 classifier. "
            "Strongly recommended for V4 amplicon data — improves genus "
            "resolution from <1%% to 40-60%%. Takes 2-4 hours."
        ),
    )

    # MiFish
    mf = parser.add_argument_group("MiFish options")
    mf.add_argument(
        "--add-gavia", action="store_true",
        help=(
            "Add Gavia (loon) 12S sequences to the MitoFish reference. "
            "Required for loon diet studies — without Gavia in the DB, "
            "host reads fall through to 'uncl. Actinopteri'. Requires --ncbi-email."
        ),
    )

    # COI
    coi = parser.add_argument_group("COI options")
    coi.add_argument(
        "--coi-version", default="261", metavar="VERSION",
        help=(
            "MIDORI2 GenBank version number (default: 261). "
            "Update to the latest release number as new DBs are published."
        ),
    )
    coi.add_argument(
        "--coi-skip-extract", action="store_true", dest="coi_skip_extract",
        help=(
            "Skip in-silico PCR for COI (NOT recommended). "
            "Training on the full ~1.5 kb COI gene degrades accuracy on "
            "Illumina short reads. Use only if providing custom amplicon-length seqs."
        ),
    )

    # cytb
    cytb_g = parser.add_argument_group("cytb options")
    cytb_g.add_argument(
        "--cytb-taxon", default="Vertebrata",
        help="NCBI taxon for cytb Entrez query (default: Vertebrata).",
    )
    cytb_g.add_argument(
        "--cytb-min-len", type=int, default=200,
        help="Minimum cytb sequence length to retain (default: 200).",
    )
    cytb_g.add_argument(
        "--cytb-max-len", type=int, default=1200,
        help="Maximum cytb sequence length to retain (default: 1200).",
    )

    # NCBI shared
    ncbi = parser.add_argument_group(
        "NCBI options (required for cytb, adenovirus, herpesvirus; "
        "and MiFish --add-gavia)"
    )
    ncbi.add_argument(
        "--ncbi-api-key", default="", metavar="KEY",
        help=(
            "NCBI Entrez API key. Required for reliable cytb downloads. "
            "Strongly recommended for virus reference building. "
            "Get one free at https://www.ncbi.nlm.nih.gov/account/"
        ),
    )
    ncbi.add_argument(
        "--ncbi-email", default="", metavar="EMAIL",
        help=(
            "Email for NCBI Entrez access (NCBI policy requirement). "
            "Required for cytb, adenovirus, herpesvirus, and MiFish --add-gavia."
        ),
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """
Entry point for 00_build_classifiers.py.

    Parses arguments, iterates over the requested markers, resolves each to a
    canonical name, and calls the appropriate builder function. Dry-run mode
    prints all commands without executing them.
    """
    args = parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    if args.f_primer or args.r_primer:
        log.info("Custom primers provided:")
        if args.f_primer:
            log.info("  Forward: %s", args.f_primer)
        if args.r_primer:
            log.info("  Reverse: %s", args.r_primer)
    else:
        log.info("No primers provided — using marker-specific defaults.")

    markers = [resolve_marker(m) for m in args.markers]
    log.info("Markers to process: %s", markers)

    for marker in markers:
        try:
            if marker == "16S":
                BUILDERS[marker](
                    outdir=outdir, f_primer=args.f_primer, r_primer=args.r_primer,
                    threads=args.threads, dry_run=args.dry_run, use_v4=args.use_v4,
                )
            elif marker == "MiFish":
                if args.add_gavia and not args.ncbi_email:
                    log.error("--add-gavia requires --ncbi-email your@institution.edu")
                    sys.exit(1)
                BUILDERS[marker](
                    outdir=outdir, f_primer=args.f_primer, r_primer=args.r_primer,
                    threads=args.threads, dry_run=args.dry_run,
                    add_gavia=args.add_gavia, ncbi_email=args.ncbi_email,
                )
            elif marker == "cytb":
                BUILDERS[marker](
                    outdir=outdir, f_primer=args.f_primer, r_primer=args.r_primer,
                    threads=args.threads, dry_run=args.dry_run,
                    ncbi_api_key=args.ncbi_api_key, ncbi_email=args.ncbi_email,
                    cytb_taxon=args.cytb_taxon,
                    cytb_min_len=args.cytb_min_len, cytb_max_len=args.cytb_max_len,
                )
            elif marker == "COI":
                BUILDERS[marker](
                    outdir=outdir, f_primer=args.f_primer, r_primer=args.r_primer,
                    threads=args.threads, dry_run=args.dry_run,
                    coi_version=args.coi_version,
                    skip_extract=args.coi_skip_extract,
                )
            elif marker in ("adenovirus", "herpesvirus"):
                if not args.ncbi_email:
                    log.error(
                        "--ncbi-email is required for viral reference building.\n"
                        "  Add: --ncbi-email your@unh.edu"
                    )
                    sys.exit(1)
                BUILDERS[marker](
                    outdir=outdir, dry_run=args.dry_run,
                    ncbi_email=args.ncbi_email, ncbi_api_key=args.ncbi_api_key,
                )
            else:
                # 18S, ITS — generic path
                BUILDERS[marker](
                    outdir=outdir, f_primer=args.f_primer, r_primer=args.r_primer,
                    threads=args.threads, dry_run=args.dry_run,
                )
        except Exception as e:
            log.error("Failed to build classifier for %s: %s", marker, e, exc_info=True)
            sys.exit(1)

    log.info("=== All done! ===")
    log.info("Outputs are in: %s", outdir.resolve())
    log.info(
        "\nTo classify rep-seqs (NB classifiers):\n"
        "  qiime feature-classifier classify-sklearn \\\n"
        "    --i-classifier classifiers/<classifier>.qza \\\n"
        "    --i-reads your-rep-seqs.qza \\\n"
        "    --o-classification your-taxonomy.qza\n"
        "\n"
        "For viral presence/absence (BLAST databases):\n"
        "  blastn -query your-rep-seqs.fasta \\\n"
        "         -db references/adenovirus/adenovirus_blast_db \\\n"
        "         -perc_identity 80 -qcov_hsp_perc 60 -evalue 1e-10 \\\n"
        "         -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' \\\n"
        "         -out adenovirus_hits.tsv"
    )


if __name__ == "__main__":
    main()
