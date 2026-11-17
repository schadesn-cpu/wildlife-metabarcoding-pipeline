#!/usr/bin/env python3
"""
Apply BLAST-verified remappings to cytb taxonomy table.
Reads cytb_remap.tsv and applies remap/exclude/review actions
to the exported taxonomy, then rebuilds the taxonomy table.
"""

import argparse
import pandas as pd
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Apply cytb BLAST remap')
    parser.add_argument('--taxonomy', required=True, help='Exported taxonomy.tsv')
    parser.add_argument('--remap', required=True, help='cytb_remap.tsv')
    parser.add_argument('--table', required=True, help='Feature table (counts TSV or QZA)')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--min-sample-reads', type=int, default=100)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load taxonomy
    tax = pd.read_csv(args.taxonomy, sep='\t', comment='#',
                      names=['Feature_ID','Taxon','Confidence'], header=0)
    print(f'Loaded {len(tax)} ASVs')

    # Load remap table
    remap = pd.read_csv(args.remap, sep='\t')
    print(f'Loaded {len(remap)} remap rules')

    # Track changes
    remapped = 0
    excluded = 0
    reviewed = 0

    # Apply remaps
    for _, rule in remap.iterrows():
        search = rule['search_string']
        action = rule['action']
        verified = rule['verified_species']

        mask = tax['Taxon'].str.contains(search, case=False, na=False)
        n_match = mask.sum()

        if n_match == 0:
            continue

        if action == 'remap':
            tax.loc[mask, 'Taxon_verified'] = verified
            tax.loc[mask, 'Action'] = 'remapped'
            remapped += n_match
            print(f'  REMAP: {search} -> {verified} ({n_match} ASVs)')
        elif action == 'exclude':
            tax.loc[mask, 'Taxon_verified'] = 'EXCLUDED'
            tax.loc[mask, 'Action'] = 'excluded'
            excluded += n_match
            print(f'  EXCLUDE: {search} ({n_match} ASVs) — {verified}')
        elif action == 'review':
            tax.loc[mask, 'Taxon_verified'] = f'REVIEW:{verified}'
            tax.loc[mask, 'Action'] = 'review'
            reviewed += n_match
            print(f'  REVIEW: {search} ({n_match} ASVs)')

    # Fill unmatched rows
    unmatched = tax['Taxon_verified'].isna()
    # Check if unresolved (stuck at Chordata)
    unresolved = unmatched & (tax['Taxon'].str.count(';') < 4)
    confirmed = unmatched & ~unresolved

    tax.loc[unresolved, 'Taxon_verified'] = 'EXCLUDED'
    tax.loc[unresolved, 'Action'] = 'unresolved_excluded'
    tax.loc[confirmed, 'Taxon_verified'] = tax.loc[confirmed, 'Taxon']
    tax.loc[confirmed, 'Action'] = 'kept_as_is'

    print(f'\n=== SUMMARY ===')
    print(f'Remapped:   {remapped:>6}')
    print(f'Excluded:   {excluded:>6}')
    print(f'Review:     {reviewed:>6}')
    print(f'Unresolved: {unresolved.sum():>6} (excluded)')
    print(f'Kept as-is: {confirmed.sum():>6}')
    print(f'Total:      {len(tax):>6}')

    if args.dry_run:
        print('\nDRY RUN — no files written')
        return

    # Save full annotated taxonomy
    tax.to_csv(outdir / 'taxonomy_blast_verified.tsv', sep='\t', index=False)
    print(f'\nSaved: {outdir}/taxonomy_blast_verified.tsv')

    # Save clean taxonomy (exclude removed, remap applied)
    clean = tax[tax['Taxon_verified'] != 'EXCLUDED'].copy()
    clean = clean[~clean['Taxon_verified'].str.startswith('REVIEW')]
    clean = clean[['Feature_ID', 'Taxon_verified', 'Confidence']].rename(
        columns={'Taxon_verified': 'Taxon'})
    clean.to_csv(outdir / 'taxonomy_clean.tsv', sep='\t', index=False)
    print(f'Saved: {outdir}/taxonomy_clean.tsv ({len(clean)} ASVs)')

    # Build summary table (collapse to verified species)
    species_counts = clean['Taxon'].value_counts()
    species_counts.to_csv(outdir / 'cytb_host_summary.tsv', sep='\t',
                          header=['ASV_count'])
    print(f'Saved: {outdir}/cytb_host_summary.tsv ({len(species_counts)} taxa)')

if __name__ == '__main__':
    main()
