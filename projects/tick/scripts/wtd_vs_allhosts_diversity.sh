#!/usr/bin/env bash
# =============================================================================
# Host-scope test for the tick 16S microbiome
# Question: does COLLECTION host structure the community, and should the main
#           figures/stats be deer-only or all-hosts?
#
# It does four things:
#   1. Collapses the host labels (WTD + Deer + White-tailed deer -> "deer"),
#      lumps the singleton hosts, and flags controls/blanks.
#   2. Rarefies the 16S table and builds a Bray-Curtis distance matrix.
#   3. Runs three adonis (PERMANOVA) tests:
#        (a) host overall
#        (b) host AFTER accounting for batch  (mirrors the state analysis)
#        (c) deer vs non-deer                 (the scope call, directly)
#   4. Makes a PCoA colored by host so you can see whether hosts separate.
#
# You run this on ron in the qiime2-amplicon-2024.5 env. It does not change
# any existing files; everything lands in $OUTDIR.
# =============================================================================
set -euo pipefail

# ===== CONFIG: edit to match your paths/columns, then run =====================
TABLE="qiime2/16S/all/dada2/table.qza"      # 16S feature table (907 libraries)
METADATA="reads_staged/16S_metadata.tsv"    # per-library 16S metadata (carries host+batch+state)
HOST_COL="sheet__Host Spp"                   # host-species column in METADATA
BATCH_COL="batch_short"                      # batch column in METADATA
DEPTH=1000                                  # rarefaction depth (matches the diversity work)
OUTDIR="qiime2/16S/all/host_scope"
# =============================================================================

mkdir -p "$OUTDIR"

# ---- Step 1: add collapsed host columns, print group counts -----------------
python3 - "$METADATA" "$HOST_COL" "$OUTDIR/meta_host.tsv" <<'PY'
import sys, csv, collections
meta_in, host_col, meta_out = sys.argv[1], sys.argv[2], sys.argv[3]
rows = list(csv.reader(open(meta_in, newline=''), delimiter='\t'))
if not rows:
    sys.exit('ERROR: metadata file is empty: %s' % meta_in)
header = rows[0]

# optional QIIME "#q2:types" directive row
body_start = 1
types_row = None
if len(rows) > 1 and rows[1] and str(rows[1][0]).startswith('#q2:types'):
    types_row = rows[1]; body_start = 2

# locate host column: exact match first, then fuzzy (contains "host" + "sp")
def find_col(name):
    for i, h in enumerate(header):
        if str(h).strip().lower() == name.strip().lower():
            return i
    for i, h in enumerate(header):
        hl = str(h).lower()
        if 'host' in hl and 'sp' in hl:
            return i
    return None

hi = find_col(host_col)
if hi is None:
    sys.exit('ERROR: could not find a host column matching "%s".\n'
             '       Columns present: %s\n'
             '       Your table-level metadata may not carry host. If so, host needs\n'
             '       to be joined from All Run Data / Tick Rosetta Stone by the 16S\n'
             '       library id first.' % (host_col, header))

def collapse(v):
    s = ('' if v is None else str(v)).strip().lower()
    if s in ('wtd', 'deer', 'white-tailed deer', 'white tailed deer',
             'whitetailed deer', 'odocoileus virginianus'):
        return 'deer'
    if 'moose' in s:     return 'moose'
    if 'bear' in s:      return 'bear'
    if 'porcupine' in s: return 'porcupine'
    if 'opossum' in s:   return 'opossum'
    if s in ('', 'na', 'none', 'blank', 'control', 'ntc', 'nan'):
        return 'control_or_blank'
    return 'other'   # singletons: squirrel, cottontail, birds, etc.

# well-represented biological hosts that can support a group-level test
KEEP = ('deer', 'bear', 'moose', 'porcupine', 'opossum')

collapsed = []
cnt = collections.Counter()
for r in rows[body_start:]:
    if not r or all(x == '' for x in r):
        collapsed.append(None); continue
    c = collapse(r[hi]); collapsed.append(c)
    if c != 'control_or_blank':
        cnt[c] += 1

sys.stderr.write('\n  Host group counts (collapsed):\n')
for h, n in cnt.most_common():
    flag = '' if h in KEEP else '   (excluded from multi-host test: too few or mixed)'
    sys.stderr.write('    {:>5}  {}{}\n'.format(n, h, flag))
sys.stderr.write('\n')

new_header = header + ['host_collapsed', 'deer_vs_other']
out = [new_header]
if types_row is not None:
    out.append(types_row + ['categorical', 'categorical'])

i = 0
for r in rows[body_start:]:
    c = collapsed[i]; i += 1
    if c is None:
        continue
    host_collapsed = c                       # full label kept; WHERE clauses select
    if c == 'control_or_blank':
        deer_vs_other = 'control_or_blank'
    elif c == 'deer':
        deer_vs_other = 'deer'
    else:
        deer_vs_other = 'nondeer'
    out.append(r + [host_collapsed, deer_vs_other])

with open(meta_out, 'w', newline='') as f:
    csv.writer(f, delimiter='\t').writerows(out)
sys.stderr.write('  Wrote %s (added host_collapsed, deer_vs_other)\n\n' % meta_out)
PY

echo ">> cleaned metadata: $OUTDIR/meta_host.tsv"

# ---- Step 2: rarefy + Bray-Curtis -------------------------------------------
# (If you already have a Bray-Curtis matrix from core-metrics, point Step 3 at it
#  and skip these two lines.)
qiime feature-table rarefy \
  --i-table "$TABLE" \
  --p-sampling-depth "$DEPTH" \
  --o-rarefied-table "$OUTDIR/rarefied.qza"

qiime diversity beta \
  --i-table "$OUTDIR/rarefied.qza" \
  --p-metric braycurtis \
  --o-distance-matrix "$OUTDIR/braycurtis_dm.qza"

DM="$OUTDIR/braycurtis_dm.qza"

# ---- Step 3: three adonis tests ---------------------------------------------
# (a) host overall, restricted to the well-represented hosts
qiime diversity filter-distance-matrix \
  --i-distance-matrix "$DM" \
  --m-metadata-file "$OUTDIR/meta_host.tsv" \
  --p-where "[host_collapsed] IN ('deer','bear','moose','porcupine','opossum')" \
  --o-filtered-distance-matrix "$OUTDIR/dm_hosts.qza"

qiime diversity adonis \
  --i-distance-matrix "$OUTDIR/dm_hosts.qza" \
  --m-metadata-file "$OUTDIR/meta_host.tsv" \
  --p-formula "host_collapsed" \
  --p-n-jobs 4 \
  --o-visualization "$OUTDIR/adonis_host_only.qzv"

# (b) host AFTER batch (sequential SS: host's marginal R2 once batch is accounted for)
qiime diversity adonis \
  --i-distance-matrix "$OUTDIR/dm_hosts.qza" \
  --m-metadata-file "$OUTDIR/meta_host.tsv" \
  --p-formula "${BATCH_COL}+host_collapsed" \
  --p-n-jobs 4 \
  --o-visualization "$OUTDIR/adonis_batch_then_host.qzv"

# (c) deer vs non-deer: the direct scope question
qiime diversity filter-distance-matrix \
  --i-distance-matrix "$DM" \
  --m-metadata-file "$OUTDIR/meta_host.tsv" \
  --p-where "[deer_vs_other] IN ('deer','nondeer')" \
  --o-filtered-distance-matrix "$OUTDIR/dm_deer_vs_other.qza"

qiime diversity adonis \
  --i-distance-matrix "$OUTDIR/dm_deer_vs_other.qza" \
  --m-metadata-file "$OUTDIR/meta_host.tsv" \
  --p-formula "deer_vs_other" \
  --p-n-jobs 4 \
  --o-visualization "$OUTDIR/adonis_deer_vs_other.qzv"

# ---- Step 4: PCoA colored by host (visual check) ----------------------------
qiime diversity pcoa \
  --i-distance-matrix "$DM" \
  --o-pcoa "$OUTDIR/braycurtis_pcoa.qza"

qiime emperor plot \
  --i-pcoa "$OUTDIR/braycurtis_pcoa.qza" \
  --m-metadata-file "$OUTDIR/meta_host.tsv" \
  --o-visualization "$OUTDIR/pcoa_host.qzv"

# ---- Done -------------------------------------------------------------------
cat <<EOF

DONE. Open these at https://view.qiime2.org (the adonis qzv has the R2 and p):

  $OUTDIR/adonis_host_only.qzv        host overall (deer/bear/moose/porcupine/opossum)
  $OUTDIR/adonis_batch_then_host.qzv  host R2/p AFTER accounting for batch
  $OUTDIR/adonis_deer_vs_other.qzv    deer vs non-deer  <- the scope call
  $OUTDIR/pcoa_host.qzv               visual: do the hosts separate?

HOW TO READ IT:
  If host R2 is small (in the ballpark of state's 0.35%) and p is not
  significant -> host is not structuring the data; keep all hosts and report
  the full set, and the deer-vs-other test will agree.

  If host R2 is sizable and significant -> restrict the main microbiome and
  pathogen figures to the deer collection-host set, or stratify by host.

Either way: the cytb feeding-ecology (blood-meal) table stays multi-host.
That is a separate axis from collection host and does not change here.
EOF
