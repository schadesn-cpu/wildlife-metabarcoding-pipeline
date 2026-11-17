#!/usr/bin/env bash
# ALL-SAMPLES (all biological ticks, host-agnostic) 16S diversity on the
# DECONTAMINATED table. Companion to the WTD-only rebuild: same pipeline, no
# host filter. Batch-colored, for Figure 1 panel A (context) alongside the
# deer-only panel B (controlled).
#
# Fresh timestamped folder; touches nothing old.
# Run: cd ~/meedlab/tick_project && bash scripts/rebuild_allsamples_diversity.sh
set -euo pipefail
ROOT=~/meedlab/tick_project; cd "$ROOT"
STAMP=$(date +%Y%m%d_%H%M)
OUT="results/diversity_allsamples_${STAMP}"
Q="qiime2/16S/all"
mkdir -p "$OUT"; echo "output dir: $OUT"

DECON_QZA="results/decontam/16S/table_decontam_16S.qza"
TREE="$Q/diversity/wtd_depth1000/rooted-tree.qza"
META_RAW="reads_staged/16S_metadata.tsv"   # has ALL samples w/ host + batch, not just WTD
DEPTH=1000

# --- 1. build all-biological-sample metadata (drop controls, keep every host) ---
python3 - <<PYEOF
import pandas as pd, re
m = pd.read_csv("$META_RAW", sep="\t")
# keep non-controls, ALL hosts
m = m[m['is_control'].astype(str).str.lower().isin(['no','false','nan',''])]
host = m['sheet__Host Spp'].fillna('unknown').replace({'0':'unknown','not recorded':'unknown','Deer':'WTD'})
out = pd.DataFrame({
    'sample-id': m['sample-id'],
    'host': host,
    'batch': m['batch_short'],
    'state': m['sheet__State'].fillna('unknown'),
})
out.to_csv("$OUT/meta_all.tsv", sep="\t", index=False)
lines=open("$OUT/meta_all.tsv").readlines()
lines.insert(1, "#q2:types\tcategorical\tcategorical\tcategorical\tcategorical\n")
open("$OUT/meta_all.tsv","w").writelines(lines)
print("all-sample metadata rows:", len(out))
print(host.value_counts().head(12).to_string())
PYEOF

# --- 2. filter decontam table to those samples, then tree-match ---
qiime feature-table filter-samples \
  --i-table "$DECON_QZA" --m-metadata-file "$OUT/meta_all.tsv" \
  --o-filtered-table "$OUT/table_all_decontam.qza"
qiime phylogeny filter-table \
  --i-table "$OUT/table_all_decontam.qza" --i-tree "$TREE" \
  --o-filtered-table "$OUT/table_all_treematched.qza"
qiime feature-table summarize \
  --i-table "$OUT/table_all_treematched.qza" \
  --o-visualization "$OUT/table_all_treematched.qzv"

# --- 3. core metrics phylogenetic ---
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$TREE" --i-table "$OUT/table_all_treematched.qza" \
  --p-sampling-depth $DEPTH --m-metadata-file "$OUT/meta_all.tsv" \
  --output-dir "$OUT/core_metrics" --verbose 2>&1 | tail -8

# --- 4. adonis: batch, host, state on the full set ---
qiime diversity adonis \
  --i-distance-matrix "$OUT/core_metrics/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "$OUT/meta_all.tsv" --p-formula "batch+host" \
  --o-visualization "$OUT/adonis_bray_batch_host.qzv"
qiime diversity adonis \
  --i-distance-matrix "$OUT/core_metrics/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "$OUT/meta_all.tsv" --p-formula "batch+state" \
  --o-visualization "$OUT/adonis_bray_batch_state.qzv"

echo ""
echo "=== DONE. All-samples diversity in $OUT ==="
echo "  table_all_treematched.qzv        (n surviving rarefaction)"
echo "  adonis_bray_batch_host.qzv       (does HOST structure community across all ticks?)"
echo "  core_metrics/bray_curtis_pcoa_results.qza  (Figure 1A ordination)"
