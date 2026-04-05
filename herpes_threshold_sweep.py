#!/bin/bash
# herpes_threshold_sweep.sh
# Run from: ~/meedlab/loon_project/
# Usage: bash scripts/herpes_threshold_sweep.sh

XLSX="metadata/loon_amplicon_analysis.xlsx"
META="metadata/qiime/metadata_16S_updated.tsv"

echo "=============================================="
echo "  HERPESVIRUS DETECTION THRESHOLD SWEEP"
echo "  $(date)"
echo "=============================================="

for THRESH in 0.01 0.05 0.10 0.25 0.50 0.75; do
    PCT=$(echo "$THRESH * 100" | bc)
    echo ""
    echo "----------------------------------------------"
    echo "  Threshold: ${PCT}%"
    echo "----------------------------------------------"
    python scripts/10_plot_viral.py \
        --xlsx "$XLSX" \
        --sheet herpes \
        --metadata "$META" \
        --group-by Group \
        --group-order Diseased Trauma \
        --threshold "$THRESH" \
        --min-sample-reads 500 \
        --palette wong \
        --dry-run 2>&1 | grep -E "INFO|WARNING" | grep -E "n=|Dropping|Samples after|Matched"
done

echo ""
echo "=============================================="
echo "  SWEEP COMPLETE"
echo "=============================================="
