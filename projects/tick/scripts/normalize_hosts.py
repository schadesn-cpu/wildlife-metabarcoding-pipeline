#!/usr/bin/env python3
"""
normalize_hosts.py — add a clean Host_clean column to crosswalk_VERIFIED.tsv.
Does NOT overwrite the raw Host column (reversible). Merges the WTD/Deer variants and standardizes.
CONFIRM the 'Deer' -> White-tailed deer call: safe here because 'Moose' is a separate label,
so 'Deer' means WTD specifically (collectors distinguished moose).
"""
import pandas as pd

HOST_MAP = {
    'WTD':'White-tailed deer', 'WTD ':'White-tailed deer', 'Deer':'White-tailed deer',
    'Moose':'Moose',
    'Bear':'Bear',
    'NA Porcupine':'Porcupine',
    'Virginia opossum':'Opossum',
    'Eastern Cottontail':'Eastern cottontail',
    'Eastern gray squirrel':'Eastern gray squirrel',
    'American Robin':'American robin',
    'American Goldfinch':'American goldfinch',
    'Common Yellowthroat':'Common yellowthroat',
}

def main():
    cw=pd.read_csv('results/metadata/crosswalk_VERIFIED.tsv',sep='\t')
    def clean(h):
        if pd.isna(h): return None
        return HOST_MAP.get(str(h).strip(), HOST_MAP.get(str(h), str(h).strip()))
    cw['Host_clean']=cw['Host'].apply(clean)
    # a coarse grouping for stratified stats (where n>=~10)
    def grp(h):
        if h is None: return None
        if h=='White-tailed deer': return 'White-tailed deer'
        if h in ('Moose','Bear'): return h
        return 'Other wildlife'
    cw['Host_group']=cw['Host_clean'].apply(grp)
    cw.to_csv('results/metadata/crosswalk_VERIFIED.tsv',sep='\t',index=False)
    print("added Host_clean + Host_group to crosswalk_VERIFIED.tsv\n")
    print("=== Host_clean counts ===")
    print(cw['Host_clean'].value_counts(dropna=False).to_string())
    print("\n=== Host_group counts (for stratified stats) ===")
    print(cw['Host_group'].value_counts(dropna=False).to_string())

if __name__=='__main__':
    main()
