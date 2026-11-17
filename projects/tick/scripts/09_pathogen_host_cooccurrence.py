#!/usr/bin/env python3
"""
Cross-marker pathogen-host co-occurrence analysis.
Matches 16S/18S pathogens with cytb blood meal hosts per sample.
"""

import argparse
import re
import pandas as pd
import subprocess
import os
from pathlib import Path


# Helper for genus-token matching (16S only when --genus-exact is passed)
_GENUS_RE = re.compile(r"g__([^;]*)")


def _genus_token(taxon):
    """Return the g__<token> from a SILVA-style Taxon string, or '' if absent."""
    if not isinstance(taxon, str):
        return ""
    m = _GENUS_RE.search(taxon)
    return m.group(1) if m else ""


def export_feature_table(qza_path, outdir):
    """Export QZA feature table to TSV."""
    biom_dir = Path(outdir) / 'biom_export'
    biom_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(['qiime', 'tools', 'export', '--input-path', qza_path,
                    '--output-path', str(biom_dir)], check=True, capture_output=True)
    tsv_path = biom_dir / 'feature-table.tsv'
    subprocess.run(['biom', 'convert', '-i', str(biom_dir / 'feature-table.biom'),
                    '-o', str(tsv_path), '--to-tsv'], check=True, capture_output=True)
    return pd.read_csv(tsv_path, sep='\t', skiprows=1, index_col=0)


def detect_in_sample(counts, taxonomy, targets, min_reads=1, match_mode='substring'):
    """For each sample, detect which target taxa are present.

    match_mode:
      'substring'   - case-insensitive substring on the Taxon string (original
                      behavior; appropriate for 18S species-level patterns and
                      cytb species/family names).
      'genus_exact' - match terms exactly against the g__<token> field of Taxon.
                      Use for 16S to avoid family/order-name substring traps
                      (e.g. 'Anaplasma' substring catching f__Anaplasmataceae,
                      thereby grouping Wolbachia/Ehrlichia under Anaplasma).
                      See results/16S/all/blast_verify/audit_substring_traps_summary.tsv
    """
    if match_mode not in ('substring', 'genus_exact'):
        raise ValueError(f"match_mode must be 'substring' or 'genus_exact', got {match_mode!r}")

    # Map feature IDs to target taxa
    target_features = {}
    for fid, row in taxonomy.iterrows():
        taxon = row['Taxon']
        if not isinstance(taxon, str):
            continue
        if match_mode == 'genus_exact':
            gen = _genus_token(taxon)
            for target_name, search_terms in targets.items():
                if any(gen == term for term in search_terms):
                    target_features[fid] = target_name
                    break
        else:  # 'substring'
            lo = taxon.lower()
            for target_name, search_terms in targets.items():
                if any(term.lower() in lo for term in search_terms):
                    target_features[fid] = target_name
                    break

    # Build per-sample detection matrix
    results = {}
    for sample in counts.columns:
        sample_counts = counts[sample]
        detected = {}
        for fid, target_name in target_features.items():
            if fid in sample_counts.index and sample_counts[fid] >= min_reads:
                if target_name not in detected:
                    detected[target_name] = 0
                detected[target_name] += sample_counts[fid]
        results[sample] = detected

    return pd.DataFrame(results).T.fillna(0)


def main():
    parser = argparse.ArgumentParser(description='Pathogen-host co-occurrence')
    parser.add_argument('--s16-table', required=True, help='16S feature table QZA')
    parser.add_argument('--s16-taxonomy', required=True, help='16S taxonomy TSV')
    parser.add_argument('--s18-table', required=True, help='18S feature table QZA')
    parser.add_argument('--s18-taxonomy', required=True, help='18S taxonomy TSV')
    parser.add_argument('--cytb-table', required=True, help='cytb feature table QZA')
    parser.add_argument('--cytb-taxonomy', required=True, help='cytb verified taxonomy TSV')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--min-reads', type=int, default=5,
                        help='Min reads to call a taxon present in a sample')
    parser.add_argument('--genus-exact', action='store_true',
                        help="For 16S only: match target names against the "
                             "g__<token> field of the SILVA Taxon string, "
                             "instead of as case-insensitive substrings. "
                             "Recommended. Default OFF for reproducibility "
                             "with prior runs. See "
                             "results/16S/all/blast_verify/audit_substring_traps_summary.tsv "
                             "for what changes when this flag is on.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Define targets
    pathogens_16s = {
        'Anaplasma': ['Anaplasma'],
        'Borreliella': ['Borreliella'],
        'Borrelia': ['Borrelia'],
        'Francisella': ['Francisella'],
        'Bartonella': ['Bartonella'],
        'Rickettsia': ['Rickettsia'],
        'Arsenophonus': ['Arsenophonus'],
        'Wolbachia': ['Wolbachia'],
    }

    pathogens_18s = {
        'Babesia_microti': ['Babesia_microti', 'Babesia;Babesia_microti'],
        'Babesia_odocoilei': ['Babesia10;Babesia10_odocoilei', 'odocoilei'],
        'Babesia_sp': ['Babesia4', 'Babesiidae'],
    }

    hosts_cytb = {
        'Deer': ['Odocoileus virginianus'],
        'Moose': ['Alces alces'],
        'Cervidae_unresolved': ['Cervidae'],
        'Porcupine': ['Erethizon', 'Erethizontidae'],
        'Bear': ['Ursus'],
        'Opossum': ['Didelphis', 'Didelphidae'],
        'Fisher': ['Pekania', 'pennanti'],
        'Deer_ked': ['Lipoptena'],
        'Rabbit': ['Sylvilagus'],
        'Raccoon': ['Procyon'],
        'Bird': ['Passeriformes', 'Geothlypis', 'Turdus', 'Spinus',
                 'Melospiza', 'Junco', 'Corvus', 'Meleagris', 'Vermivora'],
    }

    # Export and load feature tables
    print('Exporting 16S feature table...')
    counts_16s = export_feature_table(args.s16_table, str(outdir / '16S'))
    print('Exporting 18S feature table...')
    counts_18s = export_feature_table(args.s18_table, str(outdir / '18S'))
    print('Exporting cytb feature table...')
    counts_cytb = export_feature_table(args.cytb_table, str(outdir / 'cytb'))

    # Load taxonomies
    tax_16s = pd.read_csv(args.s16_taxonomy, sep='\t', comment='#',
                          names=['Feature_ID', 'Taxon', 'Confidence'], header=0).set_index('Feature_ID')
    tax_18s = pd.read_csv(args.s18_taxonomy, sep='\t', comment='#',
                          names=['Feature_ID', 'Taxon', 'Confidence'], header=0).set_index('Feature_ID')
    tax_cytb = pd.read_csv(args.cytb_taxonomy, sep='\t').set_index('Feature_ID')

    # 16S uses genus-exact matching when the flag is on; 18S and cytb always
    # use substring (their targets are species/family names, not g__ tokens).
    s16_mode = 'genus_exact' if args.genus_exact else 'substring'
    print(f'\n16S match mode: {s16_mode}')

    # Detect per sample
    print('\nDetecting 16S pathogens per sample...')
    det_16s = detect_in_sample(counts_16s, tax_16s, pathogens_16s,
                               args.min_reads, match_mode=s16_mode)
    det_16s.columns = ['16S_' + c for c in det_16s.columns]

    print('Detecting 18S parasites per sample...')
    det_18s = detect_in_sample(counts_18s, tax_18s, pathogens_18s,
                               args.min_reads, match_mode='substring')
    det_18s.columns = ['18S_' + c for c in det_18s.columns]

    print('Detecting cytb hosts per sample...')
    det_cytb = detect_in_sample(counts_cytb, tax_cytb, hosts_cytb,
                                args.min_reads, match_mode='substring')
    det_cytb.columns = ['cytb_' + c for c in det_cytb.columns]

    # Normalize sample names (strip marker suffixes for matching)
    def clean_sample(name):
        for suffix in ['-16s', '-16S', '-18s', '-18S', '-cytb', '-CYTB', '.16S', '.18S', '.cytb']:
            name = name.replace(suffix, '')
        return name.split('.')[0]  # Take base sample ID before run info

    det_16s.index = det_16s.index.map(clean_sample)
    det_18s.index = det_18s.index.map(clean_sample)
    det_cytb.index = det_cytb.index.map(clean_sample)

    # Merge on sample ID
    print('\nMerging across markers...')
    merged = det_16s.join(det_18s, how='outer').join(det_cytb, how='outer').fillna(0)

    # Convert to presence/absence
    pa = (merged > 0).astype(int)

    print(f'\nSamples with data across markers:')
    print(f'  16S: {(det_16s.sum(axis=1) > 0).sum()} samples with detections')
    print(f'  18S: {(det_18s.sum(axis=1) > 0).sum()} samples with detections')
    print(f'  cytb: {(det_cytb.sum(axis=1) > 0).sum()} samples with detections')
    print(f'  Merged: {len(merged)} total unique samples')

    # Co-occurrence summary
    print(f'\n{"=" * 60}')
    print('PATHOGEN-HOST CO-OCCURRENCE SUMMARY')
    print(f'{"=" * 60}')

    pathogen_cols = [c for c in pa.columns if c.startswith('16S_') or c.startswith('18S_')]
    host_cols = [c for c in pa.columns if c.startswith('cytb_')]

    for pathogen in pathogen_cols:
        p_name = pathogen.replace('16S_', '').replace('18S_', '')
        p_positive = pa[pa[pathogen] == 1]
        if len(p_positive) == 0:
            continue
        print(f'\n{p_name} ({len(p_positive)} positive samples):')
        for host in host_cols:
            h_name = host.replace('cytb_', '')
            cooccur = (p_positive[host] == 1).sum()
            if cooccur > 0:
                print(f'  + {h_name}: {cooccur} samples ({cooccur / len(p_positive) * 100:.1f}%)')

    # Save outputs
    merged.to_csv(outdir / 'cooccurrence_counts.tsv', sep='\t')
    pa.to_csv(outdir / 'cooccurrence_presence_absence.tsv', sep='\t')
    print(f'\nSaved: {outdir}/cooccurrence_counts.tsv')
    print(f'Saved: {outdir}/cooccurrence_presence_absence.tsv')


if __name__ == '__main__':
    main()
