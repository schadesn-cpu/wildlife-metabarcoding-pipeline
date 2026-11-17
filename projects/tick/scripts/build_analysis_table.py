#!/usr/bin/env python3
"""
build_analysis_table.py
Build the verified analysis-ready table that all pooled analyses (PoolTestR prevalence,
co-occurrence, host x pathogen) read from. ONE input, built once, verified before use.

For each qualifying tick sample, derive presence/absence of each REAL pathogen from the
feature table + verified taxonomy, join host + pool size from the crosswalk, apply exclusions.

EXCLUSIONS (things that would corrupt prevalence):
  - controls (XB/NTC/POS/PC-WELL/EXTRACTIONBLANK) — not ticks
  - VT-12 fungal-trial ticks — experimental, belong in fungal analysis
  - samples below the read threshold (pooled: >=1000; single ticks handled separately with low floor)
  - samples with no crosswalk host/pool (can't stratify or pool-model)

PATHOGENS (real, reportable) — Anaplasma, Borreliella, Borrelia, Ehrlichia (16S)
  NOT counted as pathogens (verified ENDOSYMBIONTS): Rickettsia, Francisella. Reported separately.
  NOT counted as a pathogen (CO-DETECTION, 2026-09-18): Bartonella. Section 7 reports these as
  ruminant-associated and ked-associated organisms, outside the pathogen results; 6 of the 10
  Bartonella-positive single ticks also carry deer ked reads (OR 23.4, p = 0.0001). Still counted
  and still written out, under its own heading. See ground truth 1.16.
  (Babesia comes from 18S — added in a later step / separate file; this script is 16S pathogens.)

PRESENCE RULE: a sample is positive for pathogen P if it has >= MIN_READS reads assigned to any
ASV whose verified taxon contains P. MIN_READS default 1 (any read) — can raise to be conservative.

Outputs results/analysis/analysis_table_16S.tsv:
  sample, Barcode, Host, Host_group, NumTicks, total_reads, <one 0/1 col per pathogen>, endosym flags
"""
import sys
import pandas as pd, numpy as np, zipfile, glob, tempfile, os, re
import biom

TABLE = 'results/decontam/16S/table_decontam_16S.qza'
TAX   = 'results/16S/all/taxonomy_verified/taxonomy_16S_blast_verified.tsv'
CW    = 'results/metadata/crosswalk_4marker.tsv'   # CANONICAL 2026-08-07
OUT   = 'results/analysis/analysis_table_16S.tsv'

THRESH = 1000          # pooled read-depth floor
MIN_READS_CALL = 1     # reads needed to call a pathogen present (>=1 = any read)

PATHOGENS   = ['Anaplasma','Borreliella','Borrelia','Ehrlichia']
ENDOSYMBIONT= ['Rickettsia','Francisella']   # verified endosymbionts — NOT pathogens
CO_DETECTION= ['Bartonella']                 # ked/blood-meal associated — NOT a pathogen call
CONTROL_RE  = re.compile(r'(XB|NTC|POS|PC-WELL|EXTRACTIONBL|XB-OR|BLANK)', re.I)

MARKER_TAG = re.compile(r'-(16s|18s|cytob|cb)([A-Za-z]*\d*)(-\d+)?$', re.I)


def load_table(qza):
    td=tempfile.mkdtemp()
    with zipfile.ZipFile(qza) as z: z.extractall(td)
    return biom.load_table(glob.glob(td+'/*/data/feature-table.biom')[0])

def main():
    os.makedirs('results/analysis', exist_ok=True)
    t=load_table(TABLE)

    # taxonomy: map feature_id -> taxon
    tax=pd.read_csv(TAX, sep='\t')
    idc=tax.columns[0]; txc=[c for c in tax.columns if 'axon' in c.lower()][0]
    tax_map=dict(zip(tax[idc].astype(str), tax[txc].astype(str)))

    # for each pathogen/endosymbiont, the set of feature-ids whose taxon matches
    def ids_for(name):
        # I anchor on the GENUS rank. A bare substring matches the ranks above
        # it: 'Anaplasma' is inside 'Anaplasmataceae' and picks up 9 Wolbachia
        # ASVs (91 vs 76 here), and 'Rickettsia' is inside 'Rickettsiales'.
        # Verified 2026-09-10: only Anaplasma changes; Borreliella, Borrelia,
        # Bartonella and Ehrlichia are identical either way.
        key = 'g__' + name.lower() + ';'
        return {fid for fid,tx in tax_map.items() if key in tx.lower()}
    path_ids={p: ids_for(p) for p in PATHOGENS}
    endo_ids={e: ids_for(e) for e in ENDOSYMBIONT}
    codet_ids={c: ids_for(c) for c in CO_DETECTION}

    # Name resolution is NOT done here. build_4marker_crosswalk.py resolves every
    # name in the workbook and writes results/metadata/library_lookup.tsv; this
    # script reads it through scripts/tick_lookup.py.
    #
    # The private cwmap that used to sit here was the most careful of the local
    # resolvers, and it still lost 44 samples across 10 specimens. MARKER_TAG
    # strips '-16s-2' but leaves a bare '-2' alone, so the replicate-suffix shape
    # (EQ02441291-2 through -10) matched only when the crosswalk carried that
    # exact string. Measured 2026-09-18: 623 of 672 qualifying samples against
    # the lookup's 667. None of the 44 were single ticks.
    #
    # Do not reintroduce a local resolver. An unresolved name is a crosswalk
    # problem and is reported below by name.
    sys.path.insert(0, 'scripts')
    import tick_lookup
    lookup = tick_lookup.load()

    cw = pd.read_csv(CW, sep='\t')
    specimen_meta = {}
    for _, r in cw.iterrows():
        specimen_meta[str(r['sample_id'])] = {
            'Barcode': r['sample_id'], 'Host': r['host_spp'],
            'Host_clean': r['host_spp'], 'Host_group': r['host_bin'],
            'NumTicks': r['num_ticks'],
            'is_single_tick': bool(r.get('is_single_tick', False))}
    print(f"crosswalk: {len(cw)} specimens; resolving through library_lookup.tsv")

    batches = tick_lookup.load_batches()

    def resolve(sample_id, base_id):
        """I turn a feature-table sample id into the crosswalk record, or None."""
        for marker, value in (('16S', sample_id), ('16S', base_id),
                              ('specimen', base_id)):
            hit = tick_lookup.resolve(lookup, marker, value)
            if hit is not None:
                return specimen_meta.get(str(hit))
        return None

    def resolve_batch(sample_id, base_id):
        """I return the 16S run batch for a sample, or None.

        Batch is a LIBRARY property, so I only ask the 16S namespace. A
        specimen key carries no batch by construction, and falling back to it
        would silently attach one library's plate to another's reads.
        """
        for value in (sample_id, base_id):
            hit = tick_lookup.resolve_batch(batches, '16S', value)
            if hit is not None:
                return hit
        return None

    features=set(t.ids('observation'))
    # pre-filter each pathogen's ids to those actually in the table
    path_ids={p:(v & features) for p,v in path_ids.items()}
    endo_ids={e:(v & features) for e,v in endo_ids.items()}
    codet_ids={c:(v & features) for c,v in codet_ids.items()}

    rows=[]
    excluded={'control':0,'fungal_trial':0,'below_thresh':0,
              'unresolved_name':0,'no_host_or_pool':0}
    unresolved=[]
    for s in t.ids('sample'):
        base=s.split('.')[0]
        if CONTROL_RE.search(base): excluded['control']+=1; continue
        if 'VT-12' in base: excluded['fungal_trial']+=1; continue
        col=t.data(s, 'sample')                     # reads per feature for this sample
        obs=t.ids('observation')
        total=col.sum()
        if total < THRESH: excluded['below_thresh']+=1; continue
        meta = resolve(s, base)
        # These two exclusions used to share one counter, which is how a join
        # bug hid inside a count that also held legitimate exclusions. A name
        # that resolves to nothing is a BUG; a specimen with no recorded host or
        # tick count is a DATA GAP and a correct exclusion.
        if meta is None:
            excluded['unresolved_name']+=1
            unresolved.append(base)
            continue
        if pd.isna(meta['Host']) or pd.isna(meta['NumTicks']):
            excluded['no_host_or_pool']+=1; continue
        # reads per feature as a dict only for the ids we care about (fast)
        # build index once
        rec={'sample':s,'Barcode':meta['Barcode'],'Host':meta['Host'],
             'Host_clean':meta.get('Host_clean'),'Host_group':meta['Host_group'],
             'NumTicks':meta['NumTicks'],'total_reads':int(total),
             'RunBatch':resolve_batch(s, base)}
        # map obs->index once
        # (biom returns dense array aligned to t.ids('observation'))
        idx={fid:i for i,fid in enumerate(obs)}
        for p,ids in path_ids.items():
            reads=sum(col[idx[f]] for f in ids if f in idx)
            rec[p]=int(reads); rec[p+'_pos']=int(reads>=MIN_READS_CALL)
        for e,ids in endo_ids.items():
            reads=sum(col[idx[f]] for f in ids if f in idx)
            rec[e+'_reads']=int(reads)   # endosymbiont reads (context, not a pathogen call)
        for c,ids in codet_ids.items():
            reads=sum(col[idx[f]] for f in ids if f in idx)
            rec[c]=int(reads); rec[c+'_pos']=int(reads>=MIN_READS_CALL)
        rows.append(rec)

    df=pd.DataFrame(rows)
    df.to_csv(OUT, sep='\t', index=False)

    # ---- report ----
    print("=== exclusions ===")
    for k,v in excluded.items(): print(f"  {k:16s}: {v}")
    print(f"\nanalysis table: {len(df)} samples -> {OUT}")
    print("\n=== crude pathogen POSITIVITY (any read) across all qualifying samples ===")
    for p in PATHOGENS:
        pos=df[p+'_pos'].sum()
        print(f"  {p:14s}: {pos}/{len(df)} = {100*pos/len(df):.1f}%")
    print("\n(Rickettsia/Francisella are endosymbionts — reads recorded as context, NOT pathogen calls)")
    if unresolved:
        print(f"\n!! {len(unresolved)} sample(s) did not resolve to any crosswalk specimen.")
        print("   These are NOT controls and NOT data gaps; they are unresolved names.")
        for u in sorted(set(unresolved))[:20]:
            print("      " + str(u))
        if len(set(unresolved)) > 20:
            print(f"      ... and {len(set(unresolved))-20} more")
        print("   Fix the crosswalk; do not add a fallback here.")
    print("\n=== 16S run batch, for conditioning (NOT yet stratified on) ===")
    if 'RunBatch' in df.columns:
        counts = df['RunBatch'].fillna('(none)').value_counts()
        print(f"  {len(counts)} distinct value(s) across {len(df)} samples")
        for name, n in counts.items():
            print(f"    {str(name):<44} {n:>4}")
        print("  CAUTION: these are PLATE-AND-MARKER strings from the run sheet,")
        print("  not the ~10 analysis batches used in the ked table (ground truth")
        print("  1.12) or in the R2 10.9% diversity result. Look at the grouping")
        print("  before stratifying anything on it; strata of five samples produce")
        print("  noise that looks like an answer.")

    print("\n=== co-detections (NOT pathogens; see ground truth 1.16) ===")
    for c in CO_DETECTION:
        pos=df[c+'_pos'].sum()
        print(f"  {c:14s}: {pos}/{len(df)} = {100*pos/len(df):.1f}%")
    print("  Bartonella is ked-associated and blood-meal associated. Reported in")
    print("  Results 3.7 as a co-detection, not in 3.2 and not in any co-infection")
    print("  count.")
    print("\n=== single ticks (NumTicks==1) in this table (crude prevalence sanity vs locked file) ===")
    st=df[df['NumTicks']==1]
    print("single ticks:", len(st), "| unique Barcodes:", st['Barcode'].nunique())
    for p in PATHOGENS:
        if len(st): print(f"  {p:14s}: {100*st[p+'_pos'].mean():.1f}%")
    print("\nNOTE: single-tick prevalences here use the >=1000 filter; the LOCKED single-tick")
    print("figures use n=93 (results/single_tick_v3/single_tick_v3_full.tsv, low floor), not 91.")
    print("The cohort has been 93 since the 2026-08-07 exclusion rule; see ground truth 1.6.")
    print("Expect small differences. single_tick_v3_full.tsv is canonical for single-tick")
    print("descriptives. This table is for the POOLED analyses.")

if __name__=='__main__':
    main()
