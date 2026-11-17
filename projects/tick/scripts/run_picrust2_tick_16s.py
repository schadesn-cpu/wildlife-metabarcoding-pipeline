#!/usr/bin/env bash
# run_picrust2_tick_16S.sh
# PiCRUSt2 on the tick 16S community, NSTI-FIRST.
# Phase 1 (placement + NSTI) runs, then STOPS so you can check the go/no-go
# with picrust2_nsti_check.py before predicting any functions.
# Phase 2 (KO/EC/pathway) is guarded below — flip RUN_PHASE2=1 only after the
# NSTI check says GO or MARGINAL.
#
# Run from the tick_project root.
set -euo pipefail

THREADS=4
RUN_PHASE2=0                       # <-- set to 1 only after the NSTI check passes
OUT=results/picrust2
mkdir -p "$OUT" "$OUT/inputs"

# ---------------------------------------------------------------------------
# 0. locate picrust2 (it is usually its OWN conda env, not inside qiime2)
# ---------------------------------------------------------------------------
if ! command -v place_seqs.py >/dev/null 2>&1; then
  echo "picrust2 is not on PATH. Install once (separate env), then re-run:"
  echo "    conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.5.2"
  echo "    conda activate picrust2"
  echo "Note: you'll need qiime2 active only for the two 'qiime tools export'"
  echo "lines below — run those first in the qiime2 env, then switch to picrust2."
  exit 1
fi
echo "picrust2: $(place_seqs.py --version 2>&1 | head -1)"

# ---------------------------------------------------------------------------
# 1. inputs — rep-seqs FASTA + feature-table BIOM
#    (the two 'qiime tools export' commands need the qiime2 env; if you're in
#     the picrust2 env and the files already exist, they're skipped.)
# ---------------------------------------------------------------------------
SEQS=qiime2/16S/all/rep-seqs-exported-temp/dna-sequences.fasta
if [ ! -f "$SEQS" ]; then
  qiime tools export --input-path qiime2/16S/all/dada2/rep-seqs.qza \
        --output-path "$OUT/inputs/repseqs"
  SEQS="$OUT/inputs/repseqs/dna-sequences.fasta"
fi
echo "rep-seqs: $SEQS ($(grep -c '^>' "$SEQS") ASVs)"

# canonical decontam table -> biom. (Swap to table_bacteria.qza or an
# all-samples table here if you want a different scope; decontam = canonical.)
if [ ! -f "$OUT/inputs/table/feature-table.biom" ]; then
  qiime tools export --input-path results/decontam/16S/table_decontam_16S.qza \
        --output-path "$OUT/inputs/table"
fi
BIOM="$OUT/inputs/table/feature-table.biom"
echo "table: $BIOM"

# ---------------------------------------------------------------------------
# 2. PHASE 1 — placement + 16S marker NSTI  (THE GO/NO-GO)
# ---------------------------------------------------------------------------
place_seqs.py -s "$SEQS" -o "$OUT/placed_seqs.tre" -p "$THREADS" \
      --intermediate "$OUT/intermediate/place_seqs"

hsp.py -i 16S -t "$OUT/placed_seqs.tre" \
      -o "$OUT/marker_nsti_predicted.tsv.gz" -p "$THREADS" -n

echo
echo "==================================================================="
echo " PHASE 1 done. STOP and check NSTI before predicting any function:"
echo "     python scripts/picrust2_nsti_check.py"
echo " If it says GO/MARGINAL, set RUN_PHASE2=1 at the top and re-run."
echo "==================================================================="

# ---------------------------------------------------------------------------
# 3. PHASE 2 — KO / EC / pathway  (only after NSTI check passes)
# ---------------------------------------------------------------------------
if [ "$RUN_PHASE2" -eq 1 ]; then
  echo "PHASE 2: predicting KO + pathways ..."
  hsp.py -i KO -t "$OUT/placed_seqs.tre" -o "$OUT/KO_predicted.tsv.gz" -p "$THREADS"

  metagenome_pipeline.py -i "$BIOM" -m "$OUT/marker_nsti_predicted.tsv.gz" \
        -f "$OUT/KO_predicted.tsv.gz" -o "$OUT/KO_metagenome_out" \
        --strat_out --max_nsti 2.0

  pathway_pipeline.py -i "$OUT/KO_metagenome_out/pred_metagenome_contrib.tsv.gz" \
        -o "$OUT/pathways_out" -p "$THREADS"

  add_descriptions.py -i "$OUT/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
        -m KO -o "$OUT/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz"
  add_descriptions.py -i "$OUT/pathways_out/path_abun_unstrat.tsv.gz" \
        -m METACYC -o "$OUT/pathways_out/path_abun_unstrat_descrip.tsv.gz"

  echo "PHASE 2 done. Human-readable tables:"
  echo "  $OUT/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz"
  echo "  $OUT/pathways_out/path_abun_unstrat_descrip.tsv.gz"
fi
