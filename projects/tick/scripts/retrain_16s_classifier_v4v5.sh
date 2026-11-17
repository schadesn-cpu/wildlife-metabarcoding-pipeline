#!/usr/bin/env bash
# retrain_16S_classifier_v4v5.sh
# Region-match the 16S classifier to the V4-V5 amplicon (515F/926R) and reclassify.
# Provenance (2026-06-18) confirmed the current classifier is full-length SILVA 138,
# RESCRIPt-built but never extracted to the amplicon region. This builds the matched one.
# The old taxonomy is NOT overwritten; it stays as the before/after baseline.

set -euo pipefail

# ============================ CONFIG (confirm these) ============================
# Full-length SILVA 138 reference behind the current classifier.
#   SILVA_SEQS should be the dereplicated full-length seqs (provenance reference_reads 3d95f079).
#   SILVA_TAX  should be the matching taxonomy (provenance taxa 2c4f2fe8).
#   If you have the standard RESCRIPt outputs, these are usually silva-138-99-seqs.qza /
#   silva-138-99-tax.qza. Locate with: find ~/meedlab -name 'silva-138*seqs*.qza'
SILVA_SEQS="/home/users/sne24/meedlab/tick_project/classifiers/silva-138-SSU-seqs-dna.qza"
SILVA_TAX="/home/users/sne24/meedlab/tick_project/classifiers/silva-138-SSU-tax.qza"

# The 16S representative sequences to reclassify (the ASVs, not the table).
REP_SEQS="qiime2/16S/all/dada2/rep-seqs.qza"        # CONFIRM exact filename

# Output directory and thread count.
OUTDIR="qiime2/16S/all/taxonomy"
NJOBS=8                                              # set to your allocation

# Primers (verified). Reverse primer is the canonical full 926R; swap to the 14 bp core
# CCGYCAATTYMTTT if you want to mirror the exact cutadapt trim (more permissive matching).
F_PRIMER="GTGYCAGCMGCCGCGGTAA"      # 515F
R_PRIMER="CCGYCAATTYMTTTRAGTTT"     # 926R (Parada 2015)
# ===============================================================================

mkdir -p "${OUTDIR}"

echo "[1/4] Extracting V4-V5 region from full-length SILVA..."
qiime feature-classifier extract-reads \
  --i-sequences "${SILVA_SEQS}" \
  --p-f-primer "${F_PRIMER}" \
  --p-r-primer "${R_PRIMER}" \
  --p-min-length 250 \
  --p-max-length 500 \
  --p-n-jobs "${NJOBS}" \
  --o-reads "${OUTDIR}/silva-138-99-seqs-515f-926r.qza"

echo "[2/4] Dereplicating the extracted region (uniq mode)..."
qiime rescript dereplicate \
  --i-sequences "${OUTDIR}/silva-138-99-seqs-515f-926r.qza" \
  --i-taxa "${SILVA_TAX}" \
  --p-mode 'uniq' \
  --o-dereplicated-sequences "${OUTDIR}/silva-138-99-seqs-515f-926r-derep.qza" \
  --o-dereplicated-taxa "${OUTDIR}/silva-138-99-tax-515f-926r-derep.qza"

echo "[3/4] Training the region-matched naive Bayes classifier..."
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "${OUTDIR}/silva-138-99-seqs-515f-926r-derep.qza" \
  --i-reference-taxonomy "${OUTDIR}/silva-138-99-tax-515f-926r-derep.qza" \
  --o-classifier "${OUTDIR}/silva-138-99-nb-515f-926r-classifier.qza"

echo "[4/4] Reclassifying the 16S rep-seqs (baseline kept untouched)..."
qiime feature-classifier classify-sklearn \
  --i-classifier "${OUTDIR}/silva-138-99-nb-515f-926r-classifier.qza" \
  --i-reads "${REP_SEQS}" \
  --p-reads-per-batch 10000 \
  --p-n-jobs "${NJOBS}" \
  --o-classification "${OUTDIR}/taxonomy_v4v5_regionmatched.qza"

echo "Done. New taxonomy: ${OUTDIR}/taxonomy_v4v5_regionmatched.qza"
echo "Compare unclassified fraction against the baseline taxonomy.qza."
