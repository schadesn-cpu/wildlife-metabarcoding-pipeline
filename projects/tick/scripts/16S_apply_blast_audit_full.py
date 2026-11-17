#!/usr/bin/env python3
"""
16S_apply_blast_audit_FULL.py
Apply the 16S BLAST audit corrections to the feature-level taxonomy, then
re-collapse the decontam feature table to genus using the CORRECTED taxonomy.
Mirrors the cytB apply pattern (08_apply_cytb_remap.py -> taxonomy_verified/).

Run on ron:
  python scripts/16S_apply_blast_audit_FULL.py

Inputs (all verified to align, 2026-07-08):
  qiime2/16S/all/taxonomy-exported-temp/taxonomy.tsv         (89,970 features: Feature ID, Taxon, Confidence)
  results/16S/all/blast_verify/audit_fp_Rickettsia.tsv       (741 audited ASVs)
  results/16S/all/blast_verify/audit_fp_Anaplasma.tsv        (10 audited ASVs)
  results/decontam/16S/table_decontam_16S.qza               (89,919 x 907; all 741 audit IDs present)

Outputs:
  results/16S/all/taxonomy_verified/taxonomy_16S_blast_verified.tsv   (corrected feature-level taxonomy)
  results/16S/all/taxonomy_verified/blast_corrections_applied.tsv     (audit trail: what changed)
  results/16S/all/taxonomy_verified/relabund_genus_BLASTcorrected.tsv (corrected genus x sample relabund)
  results/16S/all/taxonomy_verified/genus_before_after.tsv            (top genera: raw vs corrected)

CONSERVATIVE: only reassigns audited ASVs with BLAST confidence >= MIN_CONF and a
resolved actual_taxon (not '(none)'/blank). Others are left as their original SILVA
label and listed with applied=False.
"""
import os, glob, tempfile, zipfile
import pandas as pd
import numpy as np
import biom

MIN_CONF = 0.90
TAX_IN   = 'qiime2/16S/all/taxonomy-exported-temp/taxonomy.tsv'
AUDITS   = ['results/16S/all/blast_verify/audit_fp_Rickettsia.tsv',
            'results/16S/all/blast_verify/audit_fp_Anaplasma.tsv']
TABLE_QZA= 'results/decontam/16S/table_decontam_16S.qza'
OUTDIR   = 'results/16S/all/taxonomy_verified'


def load_audits():
    frames = []
    for f in AUDITS:
        a = pd.read_csv(f, sep='\t')
        a['source_audit'] = os.path.basename(f)
        frames.append(a)
    a = pd.concat(frames, ignore_index=True)
    is_none = a['actual_taxon'].astype(str).str.contains('(none)', regex=False)
    a['applicable'] = (a['confidence'] >= MIN_CONF) & a['actual_taxon'].notna() & (~is_none)
    return a


def genus_of(taxon):
    """last populated level of a SILVA-style lineage string."""
    t = str(taxon)
    for lvl in ('g__', 'f__', 'o__', 'c__', 'p__'):
        segs = [x.strip() for x in t.split(';') if x.strip().startswith(lvl) and not x.strip().endswith('__')]
        if segs:
            return segs[-1].split('__', 1)[1]
    return 'Unassigned'


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # 1. load feature-level taxonomy
    tax = pd.read_csv(TAX_IN, sep='\t')
    tax.columns = ['Feature ID', 'Taxon', 'Confidence'][:len(tax.columns)]
    tax = tax.set_index('Feature ID')
    print('loaded taxonomy: %d features' % len(tax))

    # 2. load audits, build the remap for applicable ASVs
    aud = load_audits()
    appl = aud[aud['applicable']].copy()
    print('audited ASVs: %d | applicable (conf>=%.2f): %d' % (len(aud), MIN_CONF, len(appl)))

    # 3. apply: overwrite Taxon for applicable audited ASVs with actual_taxon
    tax_corrected = tax.copy()
    trail = []
    for _, row in appl.iterrows():
        fid = row['feature_id']
        if fid in tax_corrected.index:
            old = tax_corrected.at[fid, 'Taxon']
            new = row['actual_taxon']
            tax_corrected.at[fid, 'Taxon'] = new
            trail.append({'feature_id': fid, 'old_taxon': old, 'new_taxon': new,
                          'current_target': row['current_target'],
                          'confidence': row['confidence'], 'total_reads': row['total_reads']})
    trail_df = pd.DataFrame(trail)
    print('reassigned %d features' % len(trail_df))

    # write corrected taxonomy + audit trail
    tax_corrected.reset_index().to_csv(f'{OUTDIR}/taxonomy_16S_blast_verified.tsv', sep='\t', index=False)
    trail_df.to_csv(f'{OUTDIR}/blast_corrections_applied.tsv', sep='\t', index=False)

    # 4. collapse decontam feature table to genus using CORRECTED taxonomy
    td = tempfile.mkdtemp()
    with zipfile.ZipFile(TABLE_QZA) as z:
        z.extractall(td)
    t = biom.load_table(glob.glob(f'{td}/*/data/feature-table.biom')[0])
    tdf = t.to_dataframe(dense=True)  # features x samples
    print('decontam table: %d features x %d samples' % tdf.shape)

    # map each feature to CORRECTED genus
    fid_to_genus = tax_corrected['Taxon'].map(genus_of)
    tdf['genus'] = tdf.index.map(fid_to_genus).fillna('Unassigned')
    genus_counts = tdf.groupby('genus').sum()

    # relative abundance per sample, then mean across samples
    relab = genus_counts.div(genus_counts.sum(axis=0), axis=1)
    relab.to_csv(f'{OUTDIR}/relabund_genus_BLASTcorrected.tsv', sep='\t')
    corrected_mean = (relab.mean(axis=1) * 100).sort_values(ascending=False)

    # 5. also compute the RAW (uncorrected) version for before/after comparison
    fid_to_genus_raw = tax['Taxon'].map(genus_of)
    tdf['genus_raw'] = tdf.index.map(fid_to_genus_raw).fillna('Unassigned')
    genus_raw = tdf.drop(columns=['genus']).groupby('genus_raw').sum()
    relab_raw = genus_raw.div(genus_raw.sum(axis=0), axis=1)
    raw_mean = (relab_raw.mean(axis=1) * 100)

    ba = pd.DataFrame({'raw_%': raw_mean, 'corrected_%': corrected_mean}).fillna(0)
    ba = ba.sort_values('corrected_%', ascending=False).head(25)
    ba.round(2).to_csv(f'{OUTDIR}/genus_before_after.tsv', sep='\t')

    print()
    print('=== TOP GENERA: raw (uncorrected) vs BLAST-corrected (mean relabund %) ===')
    print(ba.round(2).to_string())
    print()
    print('Wrote:')
    print(' ', f'{OUTDIR}/taxonomy_16S_blast_verified.tsv')
    print(' ', f'{OUTDIR}/blast_corrections_applied.tsv')
    print(' ', f'{OUTDIR}/relabund_genus_BLASTcorrected.tsv')
    print(' ', f'{OUTDIR}/genus_before_after.tsv')


if __name__ == '__main__':
    main()
