#!/usr/bin/env python3
"""
Join single-tick pathogen calls to CLEAN collection host via the positional
crosswalk (handles ALHT WTD_x <-> ALHT-Deer-x and USGSTICKS-xxxx naming).
Outputs a pathogen-by-host table for the single ticks (per-tick, no pooling issue).

Usage:
  python scripts/pathogen_by_host.py \
     --single results/single_tick_v2/single_tick_v2_full.tsv \
     --crosswalk metadata/usgs_tick_metadata.xlsx \
     --out results/single_tick_v2/single_tick_host_pathogen.tsv
"""
import argparse, re, sys
import pandas as pd
import openpyxl

def norm(s): return re.sub(r'[^a-z0-9]','',str(s).lower())

DEER_ALIASES = {'deer','wtd','white-tailed deer','whitetailed deer','white tailed deer','odocoileus virginianus'}
def norm_host(h):
    if h is None: return None
    s=str(h).strip()
    if s in ('','0','nan','not recorded','unknown'): return None
    if s.lower() in DEER_ALIASES or s=='Deer': return 'White-tailed deer'
    return s

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--single', required=True)
    ap.add_argument('--crosswalk', required=True)
    ap.add_argument('--out', required=True)
    a=ap.parse_args()

    # crosswalk: every naming string on a row -> that row's host
    wb=openpyxl.load_workbook(a.crosswalk, read_only=True, data_only=True)
    sheet='All Run Data - w Metadata'
    ws=wb[sheet] if sheet in wb.sheetnames else wb[wb.sheetnames[0]]
    rows=list(ws.iter_rows(values_only=True)); hdr=list(rows[0])
    host_i=hdr.index('Host Spp')
    name_idxs=[i for i,x in enumerate(hdr) if x in ('Barcode','Barcode2','Sample')]
    name2host={}
    for r in rows[1:]:
        h=norm_host(r[host_i])
        if h is None: continue
        for i in name_idxs:
            v=r[i]
            if v is None: continue
            k=norm(v)
            if k and k not in name2host: name2host[k]=h

    st=pd.read_csv(a.single, sep='\t')
    # I read host_bin straight from the single-tick file. It is written by
    # build_single_tick_set.py from crosswalk_4marker.tsv and is the normalised
    # host group, so 'WTD' and 'Deer' are one row rather than two. Re-deriving
    # host from the workbook here gave Deer 78 where the bloodmeal matrix gave
    # 79, and left 6 of 93 unmatched. The old prefix fallback also took
    # cands[0] from an unordered dict, which could assign a confidently wrong
    # host rather than leaving a blank.
    if 'host_bin' not in st.columns:
        sys.exit("%s has no host_bin column. Re-run build_single_tick_set.py."
                 % a.single)
    st['host']=st['host_bin'].fillna('unknown').replace('', 'unknown')

    print('single ticks:', len(st))
    print('matched to a host:', (st['host']!='unknown').sum())
    print()
    print('=== collection host of single ticks (via crosswalk) ===')
    print(st['host'].value_counts(dropna=False).to_string())

    # pathogen columns (exclude endosymbionts)
    path_cols=[c for c in st.columns if (c.startswith('16S_') or c.startswith('18S_'))
               and 'Rickettsia' not in c and 'Francisella' not in c]
    # binarize
    for p in path_cols:
        st[p+'_pos']=(pd.to_numeric(st[p],errors='coerce').fillna(0)>0).astype(int)

    print()
    print('=== pathogen positives by host (counts) ===')
    grp=st.groupby('host')[[p+'_pos' for p in path_cols]].sum()
    grp['n_ticks']=st.groupby('host').size()
    grp.columns=[c.replace('_pos','') for c in grp.columns]
    print(grp.to_string())

    st.to_csv(a.out, sep='\t', index=False)
    print('\nwrote', a.out)

if __name__=='__main__':
    main()
