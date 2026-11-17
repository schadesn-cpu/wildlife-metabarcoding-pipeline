#!/usr/bin/env python3
"""
Build ONE authoritative host assignment per sample, for all figures.
Precedence (most reliable first), recorded in host_source:
  1. collection : crosswalk 'Host Spp' (folding Deer->WTD)
  2. name       : unambiguous animal in the sample-id (Deer/WTD->WTD, Bear, Moose,
                  Porcupine, Opossum). Conservative: ambiguous VT-12 site codes
                  (TL/TW/LC/WC) are NOT guessed.
  3. unresolved : nothing reliable -> 'unknown' (bloodmeal cytb can infer separately,
                  but that is FEEDING host, not collection host; kept out of here.)

Usage:
  python scripts/build_master_host_table.py \
      --crosswalk metadata/usgs_tick_metadata.xlsx \
      --staged reads_staged/16S_metadata.tsv \
      --out results/metadata/master_host_table.tsv
"""
import argparse, re
import pandas as pd

def norm(s): return re.sub(r'[^a-z0-9]','',str(s).lower())

# conservative name-regex -> host. Only unambiguous animal words.
NAME_RULES = [
    (re.compile(r'\b(deer|wtd)\b', re.I), 'White-tailed deer'),
    (re.compile(r'deer', re.I),           'White-tailed deer'),   # ALHT-Deer-*
    (re.compile(r'bear', re.I),           'Bear'),
    (re.compile(r'moose', re.I),          'Moose'),
    (re.compile(r'porcupine', re.I),      'Porcupine'),
    (re.compile(r'opossum', re.I),        'Opossum'),
]
# host-string normalization for collection records
DEER_ALIASES = {'deer', 'wtd', 'white-tailed deer', 'whitetailed deer', 'white tailed deer', 'odocoileus virginianus'}
def norm_host(h):
    if h is None: return None
    s = str(h).strip()
    if s in ('', '0', 'nan', 'not recorded', 'unknown'): return None
    if s.lower() in DEER_ALIASES: return 'White-tailed deer'
    if s == 'Deer': return 'White-tailed deer'
    return s

CONTROL_RX = re.compile(r'XB|NTC|POS|PC-WELL|BLANK|EXTRACTION|PAC', re.I)
def is_control(sid):
    return bool(CONTROL_RX.search(str(sid)))

def name_host(sid):
    s = str(sid)
    # guard: do NOT match the VT-12 site codes as hosts
    for rx, host in NAME_RULES:
        if rx.search(s):
            return host
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--crosswalk', required=True)
    ap.add_argument('--staged', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    # --- crosswalk: build normalized-name -> collection host across all naming cols ---
    import openpyxl
    wb = openpyxl.load_workbook(a.crosswalk, read_only=True, data_only=True)
    # prefer the compiled sheet if present
    sheet = 'All Run Data - w Metadata'
    ws = wb[sheet] if sheet in wb.sheetnames else wb[wb.sheetnames[0]]
    rows = list(ws.iter_rows(values_only=True))
    hdr = list(rows[0])
    host_i = hdr.index('Host Spp')
    name_idxs = [i for i,h in enumerate(hdr) if h in ('Barcode','Barcode2','Sample')]
    xwalk = {}
    for r in rows[1:]:
        h = norm_host(r[host_i])
        if h is None: continue
        for i in name_idxs:
            v = r[i]
            if v is None: continue
            k = norm(v)
            if k and k not in xwalk:
                xwalk[k] = h

    # --- staged sample list ---
    m = pd.read_csv(a.staged, sep='\t', dtype=str)
    out = []
    for sid in m['sample-id']:
        if is_control(sid):
            out.append({'sample-id': sid, 'host': 'CONTROL', 'host_source': 'control'})
            continue
        k = norm(re.sub(r'-16s.*|-18s.*|-its.*|\..*','',str(sid), flags=re.I))
        host, src = None, 'unresolved'
        # 1. collection via crosswalk (exact then prefix)
        h = xwalk.get(k)
        if h is None:
            cands = [v for kk,v in xwalk.items() if kk.startswith(k)]
            h = cands[0] if cands else None
        if h is not None:
            host, src = h, 'collection'
        else:
            # 2. name-regex
            nh = name_host(sid)
            if nh is not None:
                host, src = nh, 'name'
            else:
                host, src = 'unknown', 'unresolved'
        out.append({'sample-id': sid, 'host': host, 'host_source': src})

    res = pd.DataFrame(out)
    import os
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    res.to_csv(a.out, sep='\t', index=False)
    print('wrote', a.out, 'rows:', len(res))
    print('\nhost:')
    print(res['host'].value_counts(dropna=False).to_string())
    print('\nhost_source:')
    print(res['host_source'].value_counts().to_string())
    print('\nstill unresolved BIOLOGICAL (controls excluded; sample of up to 30):')
    bio_unres = res[(res['host']=='unknown') & (res['host_source']=='unresolved')]
    for s in bio_unres['sample-id'].head(30):
        print('   ', s)
    print('\n  biological unresolved total:', len(bio_unres))
    print('  controls flagged:', (res['host']=='CONTROL').sum())

if __name__ == '__main__':
    main()
