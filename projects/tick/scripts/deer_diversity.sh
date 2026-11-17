#!/usr/bin/env bash
# =============================================================================
# Deer-set 16S diversity, regenerated on the CORRECTED host collapse
# Run after the host-scope test, which showed collection host structures the
# community (host R2 ~8.75%, p=0.001), so the main microbiome figures are built
# on the deer collection-host set only.
#
# This collapses host correctly (WTD + Deer + White-tailed deer -> deer, which
# the old Odocoileus-virginianus filter missed), filters the table to deer,
# runs the core-metrics suite, redoes the state effect WITHIN deer, and pulls
# the core microbiota. Nothing existing is overwritten; output lands in $OUTDIR.
# =============================================================================
set -euo pipefail

# ===== CONFIG =================================================================
TABLE="qiime2/16S/all/dada2/table.qza"
METADATA="reads_staged/16S_metadata.tsv"
DEPTH=1000
OUTDIR="qiime2/16S/all/deer_diversity"
# Optional: a rooted tree enables UniFrac + Faith PD. Leave empty to skip.
TREE=""    # e.g. qiime2/16S/all/phylogeny/rooted-tree.qza
# =============================================================================

mkdir -p "$OUTDIR"

# ---- Step 1: minimal deer-tagged metadata (host collapsed, state+batch kept, free text dropped)
python3 - "$METADATA" "$OUTDIR/meta_deer.tsv" <<'PY'
import sys, csv
src, out = sys.argv[1], sys.argv[2]
rows = list(csv.reader(open(src, newline=''), delimiter='\t'))
hdr = rows[0]
def idx(*subs):
    for i, h in enumerate(hdr):
        hl = str(h).strip().lower()
        if all(s in hl for s in subs):
            return i
    return None
hi, bi, si = idx('host', 'sp'), idx('batch'), idx('state')
body = rows[1:]
if body and str(body[0][0]).startswith('#q2:types'):
    body = body[1:]
def collapse(v):
    s = ('' if v is None else str(v)).strip().lower()
    if s in ('wtd', 'deer', 'white-tailed deer', 'white tailed deer',
             'whitetailed deer', 'odocoileus virginianus'):
        return 'deer'
    for k in ('moose', 'bear', 'porcupine', 'opossum'):
        if k in s:
            return k
    if s in ('', 'na', 'none', 'blank', 'control', 'ntc', 'nan'):
        return 'control_or_blank'
    return 'other'
with open(out, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['sample-id', 'batch_short', 'state', 'host_collapsed'])
    w.writerow(['#q2:types', 'categorical', 'categorical', 'categorical'])
    for r in body:
        if not r or all(x == '' for x in r):
            continue
        host = collapse(r[hi])
        batch = (r[bi] if (bi is not None and bi < len(r)) else '')
        state = (r[si] if (si is not None and si < len(r)) else '')
        w.writerow([r[0], str(batch).strip(), str(state).strip(), host])
print('wrote', out)
PY
echo ">> deer-tagged metadata: $OUTDIR/meta_deer.tsv"

# ---- Step 2: filter table to deer ------------------------------------------
qiime feature-table filter-samples \
  --i-table "$TABLE" \
  --m-metadata-file "$OUTDIR/meta_deer.tsv" \
  --p-where "[host_collapsed]='deer'" \
  --o-filtered-table "$OUTDIR/table_deer.qza"

qiime feature-table summarize \
  --i-table "$OUTDIR/table_deer.qza" \
  --o-visualization "$OUTDIR/table_deer_summary.qzv"

# ---- Step 3: core-metrics suite at DEPTH -----------------------------------
# core-metrics refuses to write into an existing directory, so clear a prior run.
rm -rf "$OUTDIR/core_metrics"
qiime diversity core-metrics \
  --i-table "$OUTDIR/table_deer.qza" \
  --p-sampling-depth "$DEPTH" \
  --m-metadata-file "$OUTDIR/meta_deer.tsv" \
  --output-dir "$OUTDIR/core_metrics"

BC="$OUTDIR/core_metrics/bray_curtis_distance_matrix.qza"

# ---- Step 4: state effect WITHIN deer (state after batch) ------------------
# restrict to the three real states so empty/NA state cannot break adonis
qiime diversity filter-distance-matrix \
  --i-distance-matrix "$BC" \
  --m-metadata-file "$OUTDIR/meta_deer.tsv" \
  --p-where "[state] IN ('VT','NH','ME')" \
  --o-filtered-distance-matrix "$OUTDIR/dm_deer_state.qza"

qiime diversity adonis \
  --i-distance-matrix "$OUTDIR/dm_deer_state.qza" \
  --m-metadata-file "$OUTDIR/meta_deer.tsv" \
  --p-formula "batch_short+state" \
  --p-n-jobs 4 \
  --o-visualization "$OUTDIR/adonis_deer_batch_then_state.qzv"

# ---- Step 5: alpha diversity by state --------------------------------------
qiime diversity alpha-group-significance \
  --i-alpha-diversity "$OUTDIR/core_metrics/shannon_vector.qza" \
  --m-metadata-file "$OUTDIR/meta_deer.tsv" \
  --o-visualization "$OUTDIR/alpha_shannon_by_state.qzv"

# ---- Step 6: core microbiota of the deer set -------------------------------
qiime feature-table core-features \
  --i-table "$OUTDIR/table_deer.qza" \
  --o-visualization "$OUTDIR/core_features_deer.qzv"

# ---- Optional: phylogenetic metrics if a rooted tree is set ----------------
if [ -n "$TREE" ]; then
  rm -rf "$OUTDIR/core_metrics_phylo"
  qiime diversity core-metrics-phylogenetic \
    --i-table "$OUTDIR/table_deer.qza" \
    --i-phylogeny "$TREE" \
    --p-sampling-depth "$DEPTH" \
    --m-metadata-file "$OUTDIR/meta_deer.tsv" \
    --output-dir "$OUTDIR/core_metrics_phylo"
fi

cat <<EOF

DONE. The deer-set diversity is regenerated on the corrected 748-deer collapse.

Read out at https://view.qiime2.org :
  $OUTDIR/table_deer_summary.qzv          confirm the deer sample count (~748 before rarefaction)
  $OUTDIR/core_metrics/                   PCoA (bray_curtis_emperor.qzv), alpha vectors
  $OUTDIR/adonis_deer_batch_then_state.qzv  state effect WITHIN deer (replaces the old all-host 0.35%)
  $OUTDIR/alpha_shannon_by_state.qzv      alpha diversity by state
  $OUTDIR/core_features_deer.qzv          core microbiota of deer ticks

Pull the adonis numbers straight out with:
  unzip -p $OUTDIR/adonis_deer_batch_then_state.qzv '*adonis.tsv'
EOF
