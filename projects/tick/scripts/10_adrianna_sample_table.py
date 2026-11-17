#!/usr/bin/env python3
"""Build Adrianna's per-sample table with host, pathogen, parasite, and location."""

import pandas as pd

pa = pd.read_csv('results/cooccurrence/cooccurrence_presence_absence.tsv', sep='\t', index_col=0)
pa = pa.groupby(pa.index).max()
counts = pd.read_csv('results/cooccurrence/cooccurrence_counts.tsv', sep='\t', index_col=0)
counts = counts.groupby(counts.index).sum()

meta = pd.read_csv('metadata/metadata_ticks.tsv', sep='\t')
meta['Barcode_clean'] = meta['Barcode '].str.strip()

sample_meta = {}
for _, row in meta.iterrows():
    bc = row['Barcode_clean']
    if pd.notna(bc):
        for s in pa.index:
            if bc in s:
                sample_meta[s] = {
                    'Coordinates': row.get('Coordinates', ''),
                    'Town': row.get('Town ', ''),
                    'State': str(row.get('State', '')).strip(),
                    'Host_Spp': row.get('Host Spp', ''),
                    'N_Ticks': row.get('# of Ticks', ''),
                }

cytb_cols = [c for c in counts.columns if c.startswith('cytb_')]
s16_cols = [c for c in counts.columns if c.startswith('16S_') and 'Rickettsia' not in c]
s18_cols = [c for c in counts.columns if c.startswith('18S_')]
rick_cols = [c for c in counts.columns if 'Rickettsia' in c and '16S' in c]

rows = []
for sample in pa.index:
    row = {'Sample_ID': sample}

    # cytb hosts (top 2)
    cv = counts.loc[sample, cytb_cols].squeeze()
    cp = cv[cv > 0].sort_values(ascending=False)
    row['cytb_host_1'] = cp.index[0].replace('cytb_', '') if len(cp) > 0 else ''
    row['cytb_host_2'] = cp.index[1].replace('cytb_', '') if len(cp) > 1 else ''

    # 16S pathogens (top 2, exclude Rickettsia endosymbiont)
    sv = counts.loc[sample, s16_cols].squeeze()
    sp = sv[sv > 0].sort_values(ascending=False)
    row['pathogen_16S_1'] = sp.index[0].replace('16S_', '') if len(sp) > 0 else ''
    row['pathogen_16S_2'] = sp.index[1].replace('16S_', '') if len(sp) > 1 else ''

    # 18S parasites (top 2)
    ev = counts.loc[sample, s18_cols].squeeze()
    ep = ev[ev > 0].sort_values(ascending=False)
    row['parasite_18S_1'] = ep.index[0].replace('18S_', '') if len(ep) > 0 else ''
    row['parasite_18S_2'] = ep.index[1].replace('18S_', '') if len(ep) > 1 else ''

    # Rickettsia
    row['Rickettsia_present'] = 1 if any(counts.loc[sample, rick_cols] > 0) else 0

    # Metadata
    if sample in sample_meta:
        m = sample_meta[sample]
        row['Town'] = m.get('Town', '')
        row['State'] = m.get('State', '')
        row['Coordinates'] = m.get('Coordinates', '')
        row['Host_Spp_metadata'] = m.get('Host_Spp', '')
        row['N_Ticks'] = m.get('N_Ticks', '')

    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv('results/cooccurrence/adrianna_sample_table.tsv', sep='\t', index=False)

print(f'Saved: results/cooccurrence/adrianna_sample_table.tsv')
print(f'Samples: {len(df)}')
print(f'Columns: {list(df.columns)}')
print()
has_data = df[(df['cytb_host_1'] != '') | (df['pathogen_16S_1'] != '')]
print(f'Samples with any data: {len(has_data)}')
print()
print('Preview:')
print(has_data[['Sample_ID', 'cytb_host_1', 'pathogen_16S_1', 'parasite_18S_1', 'State']].head(15).to_string())

