#!/usr/bin/env python3
"""
Apply 16S BLAST audit corrections to the taxonomy/abundance table.
Analogous to 08_apply_cytb_remap.py, but reassigns by feature_id (each audited
ASV -> its BLAST-verified actual_taxon).

Reads:  results/16S/all/blast_verify/audit_fp_Rickettsia.tsv
        results/16S/all/blast_verify/audit_fp_Anaplasma.tsv
Writes: results/16S/all/taxonomy_verified/blast_corrections_applied.tsv (audit trail)

CONSERVATIVE: only counts reassignments with BLAST confidence >= MIN_CONF and a
resolved actual_taxon (not '(none)'/blank). Low-confidence audited ASVs are listed
but flagged applicable=False.

MODES:
  report : show what WOULD change and the corrected genus split, using the audit's
           own read counts. No per-ASV feature table needed. Safe/fast/honest.
  apply  : (stub) rewrite per-ASV taxonomy by feature_id then re-aggregate. Needs
           the ASV x sample feature table; wire in once its path is confirmed.
"""
import argparse
import os
import sys
import pandas as pd

MIN_CONF = 0.90
AUDITS = [
    'results/16S/all/blast_verify/audit_fp_Rickettsia.tsv',
    'results/16S/all/blast_verify/audit_fp_Anaplasma.tsv',
]


def genus_of(taxon):
    t = str(taxon)
    for lvl in ('g__', 'f__', 'o__'):
        segs = [x for x in t.split(';') if x.startswith(lvl) and not x.endswith('__')]
        if segs:
            return segs[-1].split('__')[1]
    return None


def load_audits():
    frames = []
    for f in AUDITS:
        if os.path.exists(f):
            a = pd.read_csv(f, sep='\t')
            a['source_audit'] = os.path.basename(f)
            frames.append(a)
    if not frames:
        sys.exit('no audit files found')
    a = pd.concat(frames, ignore_index=True)
    a['actual_genus_parsed'] = a['actual_taxon'].map(genus_of)
    is_none = a['actual_taxon'].astype(str).str.contains('(none)', regex=False)
    a['applicable'] = (a['confidence'] >= MIN_CONF) & a['actual_taxon'].notna() & (~is_none)
    return a


def report(a):
    print('audited ASVs: %d | applicable (conf>=%.2f, resolved): %d'
          % (len(a), MIN_CONF, a['applicable'].sum()))
    print()
    ap = a[a['applicable']].copy()
    print('=== what the audited calls ACTUALLY are (applicable only) ===')
    grp = ap.groupby(['current_target', 'actual_genus_parsed'])['total_reads'].agg(['sum', 'count'])
    print(grp.to_string())
    print()
    print('=== corrected picture per target ===')
    for tgt in ap['current_target'].unique():
        sub = ap[ap['current_target'] == tgt]
        tot = sub['total_reads'].sum()
        print('%s: %s audited reads reassign to ->' % (tgt, format(int(tot), ',')))
        for g, v in sub.groupby('actual_genus_parsed')['total_reads'].sum().sort_values(ascending=False).items():
            print('    %s: %s reads (%.1f%%)' % (g, format(int(v), ','), 100 * v / tot))
    print()
    print('To correct the taxonomy table: subtract these reassigned reads from the')
    print('raw genus totals (Rickettsia, Anaplasma) and add them to the corrected')
    print('genus (Diplorickettsia, Wolbachia, etc.), then re-normalize.')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=['report', 'apply'], default='report')
    p.add_argument('--feature-table', help='ASV x sample table (required for apply)')
    args = p.parse_args()

    a = load_audits()
    outdir = 'results/16S/all/taxonomy_verified'
    os.makedirs(outdir, exist_ok=True)
    a.to_csv('%s/blast_corrections_applied.tsv' % outdir, sep='\t', index=False)

    if args.mode == 'report':
        report(a)
        print('\naudit trail written: %s/blast_corrections_applied.tsv' % outdir)
    else:
        if not args.feature_table:
            sys.exit('--mode apply requires --feature-table (ASV x sample). Run report first.')
        print('apply mode stub: reassign per-ASV taxonomy by feature_id, then re-aggregate.')
        print('Implement once the ASV-level feature-table path is confirmed.')


if __name__ == '__main__':
    main()
