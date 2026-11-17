#!/usr/bin/env bash
# Host-controlled (WTD-only) 16S diversity rebuild on the DECONTAMINATED table.
# Writes everything to a fresh timestamped folder; touches NOTHING in the old
# wtd_depth1000 results. Compare new vs old afterward.
#
# Chain: WTD metadata -> filter decontam table to WTD -> rarefy 1000 ->
#        core-metrics-phylogenetic (Faith PD, UniFrac via existing rooted tree) ->
#        adonis host-controlled models (batch / state / town) -> export stats.
#
# Run: cd ~/meedlab/tick_project && bash scripts/rebuild_wtd_diversity.sh
set -euo pipefail

ROOT=~/meedlab/tick_project
cd "$ROOT"
STAMP=$(date +%Y%m%d_%H%M)
OUT="results/diversity_rebuild_${STAMP}"
Q="qiime2/16S/all"
mkdir -p "$OUT"
echo "output dir: $OUT"

DECON_QZA="results/decontam/16S/table_decontam_16S.qza"
TREE="$Q/diversity/wtd_depth1000/rooted-tree.qza"   # existing rooted tree (superset)
META_RAW="reads_staged/16S_metadata_wtd.tsv"
DEPTH=1000

# --- 1. build a clean WTD, non-control QIIME metadata (State + Town grouping) ---
python3 - <<PYEOF
import pandas as pd
m = pd.read_csv("$META_RAW", sep="\t")
# 'WTD' and 'Deer' are the SAME host (white-tailed deer), labeled two ways in
# the metadata. Include BOTH and normalize to WTD, else ~160 deer are dropped.
m = m[m['sheet__Host Spp'].isin(['WTD','Deer'])]
m = m[m['is_control'].astype(str).str.lower().isin(['no','false','nan',''])]
out = pd.DataFrame({
    'sample-id': m['sample-id'],
    'host': 'WTD',   # normalized (was WTD or Deer)
    'batch': m['batch_short'],
    'state': m['sheet__State'],
    'town':  m['sheet__Town'].fillna('not_recorded').replace('not recorded','not_recorded'),
})
out.to_csv("$OUT/meta_wtd.tsv", sep="\t", index=False)
lines=open("$OUT/meta_wtd.tsv").readlines()
lines.insert(1, "#q2:types\tcategorical\tcategorical\tcategorical\tcategorical\n")
open("$OUT/meta_wtd.tsv","w").writelines(lines)
print("WTD+Deer (normalized) metadata rows:", len(out))
PYEOF

# --- 2. filter decontam table to WTD samples ---
qiime feature-table filter-samples \
  --i-table "$DECON_QZA" \
  --m-metadata-file "$OUT/meta_wtd.tsv" \
  --o-filtered-table "$OUT/table_wtd_decontam.qza"

echo "=== WTD decontam table summary ==="
qiime feature-table summarize \
  --i-table "$OUT/table_wtd_decontam.qza" \
  --o-visualization "$OUT/table_wtd_decontam.qzv"

# --- 2b. filter table to features present in the tree (phylo metrics require this) ---
# The rooted tree was built from a different sample/ASV set, so the decontam WTD
# table can contain ASVs absent from the tree, which crashes core-metrics-phylogenetic.
# Dropping tree-absent features (typically low-abundance) lets phylo metrics run.
qiime phylogeny filter-table \
  --i-table "$OUT/table_wtd_decontam.qza" \
  --i-tree "$TREE" \
  --o-filtered-table "$OUT/table_wtd_treematched.qza"

echo "=== tree-matched table summary (compare feature count to above) ==="
qiime feature-table summarize \
  --i-table "$OUT/table_wtd_treematched.qza" \
  --o-visualization "$OUT/table_wtd_treematched.qzv"

# --- 3. core metrics phylogenetic at depth 1000 (on tree-matched table) ---
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$TREE" \
  --i-table "$OUT/table_wtd_treematched.qza" \
  --p-sampling-depth $DEPTH \
  --m-metadata-file "$OUT/meta_wtd.tsv" \
  --output-dir "$OUT/core_metrics" \
  --verbose 2>&1 | tail -8

# --- 4. adonis: host-controlled batch / state / town ---
# batch + state jointly (the key model: does state matter after batch?)
qiime diversity adonis \
  --i-distance-matrix "$OUT/core_metrics/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "$OUT/meta_wtd.tsv" \
  --p-formula "batch+state" \
  --o-visualization "$OUT/adonis_bray_batch_state.qzv"

qiime diversity adonis \
  --i-distance-matrix "$OUT/core_metrics/weighted_unifrac_distance_matrix.qza" \
  --m-metadata-file "$OUT/meta_wtd.tsv" \
  --p-formula "batch+state" \
  --o-visualization "$OUT/adonis_wunifrac_batch_state.qzv"

# batch + town (Lawrence's finer geographic grouping)
qiime diversity adonis \
  --i-distance-matrix "$OUT/core_metrics/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "$OUT/meta_wtd.tsv" \
  --p-formula "batch+town" \
  --o-visualization "$OUT/adonis_bray_batch_town.qzv"

echo ""
echo "=== DONE. New results in $OUT ==="
echo "Key outputs to compare against OLD wtd_depth1000:"
echo "  $OUT/table_wtd_decontam.qzv        (sample count surviving rarefaction)"
echo "  $OUT/adonis_bray_batch_state.qzv   (batch% / state% R2 -- compare to old 10.9% / 0.3%)"
echo "  $OUT/core_metrics/*_vector.qza     (alpha: Faith PD, Shannon, evenness)"
echo "  $OUT/core_metrics/*_pcoa_results.qza (for PCoA plots)"
