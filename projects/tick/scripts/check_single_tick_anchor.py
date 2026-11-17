#!/usr/bin/env python3
"""
check_single_tick_anchor.py
Purpose: join per-sample 16S read depth onto the single-tick table, then test whether
single-tick pathogen prevalence is stable as you raise the read-depth floor.
If prevalence stays flat, the single-tick "anchor" is not depth-fragile (good).

Design: DISCOVERS column names instead of hardcoding them, and prints a CHECKPOINT after
each step. If a step can't proceed it stops with an explICIT message telling you what to fix,
so nothing silently joins wrong.

Run from the project root:  python scripts/check_single_tick_anchor.py
"""
import sys, re
import pandas as pd, numpy as np

# ---- paths (edit here if your layout differs) ----
P_SINGLE = 'results/single_tick_v2/single_tick_v2_full.tsv'
P_XWALK  = 'results/metadata/crosswalk_VERIFIED.tsv'
P_BIOM   = 'results/decontam/16S/table_decontam_16S.biom'   # already on ron
PATHS    = ['Anaplasma', 'Borreliella', 'Borrelia', 'Bartonella']
FLOORS   = [0, 10000, 20000, 40000]

def ck(msg): print(f'[CHECKPOINT] {msg}')
def die(msg):
    print(f'\n[STOP] {msg}')
    sys.exit(1)

def pick_col(cols, patterns, label, sample_row=None):
    """Return the first column whose name matches any regex in patterns; else report options."""
    for pat in patterns:
        for c in cols:
            if re.search(pat, str(c), re.I):
                return c
    print(f'\n[STOP] could not auto-detect the {label} column.')
    print('       columns available:', list(cols))
    if sample_row is not None:
        print('       first row:', sample_row.to_dict())
    print(f'       -> edit pick_col patterns for {label}, or set it explicitly.')
    sys.exit(1)

# ---------------------------------------------------------------- 1. load single-tick
try:
    st = pd.read_csv(P_SINGLE, sep='\t')
except FileNotFoundError:
    die(f'single-tick file not found at {P_SINGLE}')
ck(f'loaded single-tick table: {st.shape[0]} rows, {st.shape[1]} cols')

# confirm pathogen columns exist (they are prefixed 16S_ in the v2 file)
patho_cols = {}
for p in PATHS:
    hit = [c for c in st.columns if re.fullmatch(rf'(16S_)?{p}(_pos)?', str(c), re.I)]
    if not hit:
        die(f'no column for pathogen {p} in single-tick file. cols: {list(st.columns)}')
    patho_cols[p] = hit[0]
ck(f'pathogen columns resolved: {patho_cols}')

st_bc = pick_col(st.columns, [r'^barcode$'], 'single-tick barcode', st.iloc[0])
st['_bc'] = st[st_bc].astype(str).str.strip()

# ---------------------------------------------------------------- 2. load crosswalk, find keys
try:
    xw = pd.read_csv(P_XWALK, sep='\t')
except FileNotFoundError:
    die(f'crosswalk not found at {P_XWALK}')
ck(f'loaded crosswalk: {xw.shape[0]} rows, {xw.shape[1]} cols')
print('       crosswalk columns:', list(xw.columns))

xw_bc  = pick_col(xw.columns, [r'^barcode$', r'barcode'], 'crosswalk barcode', xw.iloc[0])
# the 16S sample-id column: try common spellings
xw_16s = pick_col(xw.columns,
                  [r'^s?_?16s.*sample', r'sample.*16s', r'^16s.*sample', r'^s16s$', r'16s.*name', r'sample_?16s'],
                  'crosswalk 16S-sample-id', xw.iloc[0])
xw['_bc'] = xw[xw_bc].astype(str).str.strip()
ck(f'crosswalk keys: barcode="{xw_bc}", 16S-sample-id="{xw_16s}"')

# barcode overlap sanity
overlap = len(set(st['_bc']) & set(xw['_bc']))
ck(f'barcode overlap single-tick vs crosswalk: {overlap}/{len(st)}')
if overlap < 0.5 * len(st):
    die('barcode overlap is low. keys probably need cleaning (whitespace/case/suffix). '
        'Inspect the two _bc columns before trusting any join.')

# ---------------------------------------------------------------- 3. per-sample depth from biom
try:
    import biom
    t = biom.load_table(P_BIOM)
    depth = pd.Series(np.asarray(t.sum(axis='sample')).ravel(), index=t.ids('sample'), name='total_reads')
    depth = depth.reset_index().rename(columns={'index': '_16s_id'})
except Exception as e:
    die(f'could not load/parse biom at {P_BIOM}: {e}')
ck(f'loaded depth for {len(depth)} samples from feature table '
   f'(median={depth.total_reads.median():.0f})')

# ---------------------------------------------------------------- 4. join
st2 = (st.merge(xw[['_bc', xw_16s]], on='_bc', how='left')
         .merge(depth, left_on=xw_16s, right_on='_16s_id', how='left'))
missing = st2['total_reads'].isna().sum()
ck(f'depth joined; {len(st2)-missing}/{len(st2)} single ticks have a depth value '
   f'({missing} missing)')
if missing > 0:
    ex = st2.loc[st2.total_reads.isna(), ['_bc', xw_16s]].head(8)
    print('       examples with no depth (barcode -> 16S id):')
    print(ex.to_string(index=False))
    print('       (missing usually = barcode present but its 16S sample id not in the feature table)')


# after the join, before the stability table:
# collapse replicates to one row per specimen (Barcode), pathogen-positive if ANY replicate positive,
# and take the MAX depth across replicates (best-sequenced replicate represents the specimen)
agg = {c: 'max' for c in [patho_cols[p] for p in PATHS]}
agg['total_reads'] = 'max'
st2 = st2.groupby('_bc', as_index=False).agg(agg)
print(f'[CHECKPOINT] collapsed to {len(st2)} specimens')  # should print 103

# ---------------------------------------------------------------- 5. anchor stability
print('\n=== ANCHOR STABILITY: single-tick prevalence at rising depth floors ===')
print('(if these stay flat across floors, the single-tick anchor is NOT depth-fragile)')
hdr = 'pathogen'.ljust(12) + ' | ' + ' | '.join(f'>= {f}'.rjust(11) for f in FLOORS)
print(hdr); print('-'*len(hdr))
for p in PATHS:
    col = patho_cols[p]
    cells = []
    for f in FLOORS:
        g = st2[st2['total_reads'].fillna(0) >= f]
        k = int((g[col] == 1).sum()); n = len(g)
        cells.append(f'{100*k/n:4.1f}% (n={n})'.rjust(11) if n else 'NA'.rjust(11))
    print(p.ljust(12) + ' | ' + ' | '.join(cells))

print('\n[DONE] if the row values hold roughly constant left-to-right, the anchor is solid.')
print('       write these numbers into the validation memo / master log.')
