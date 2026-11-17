#!/usr/bin/env python3
"""
cytb_human_checkpoint.py
Reusable QC + epidemiology checkpoint: detect and quantify HUMAN cytb signal in tick data.

Serves two purposes:
  (1) QC — how much human handler contamination is in this run.
  (2) Epidemiology — which ticks may have taken a HUMAN blood meal (with contamination caveat).

KEY DESIGN (learned from the 2024-2025 USGS tick data):
  - Human DNA is ubiquitous in labs; trace human across many samples = handler contamination.
  - A real human blood meal shows human as a LARGE fraction of a WELL-POPULATED sample.
  - So a candidate requires BOTH: human_pct >= FRACTION_MIN AND human_reads >= READS_MIN.
    (A 27% human read in a 47-read sample is noise, not a meal.)
  - Contamination baseline must come from CLEAN blanks only. Blanks that contain a stray
    tick/ked have host DNA and are COMPROMISED — exclude them. Positive controls may contain
    human by design — exclude them too.
  - Report single (pool=1) vs pooled separately.
  - ALWAYS output the caveat: handler contamination cannot be fully excluded.

Usage:
  python cytb_human_checkpoint.py \
      --table qiime2/cytb/all/dada2/table.qza \
      --taxonomy qiime2/cytb/all/taxonomy-exported-temp/taxonomy.tsv \
      --metadata metadata/tick_specimen_master.csv \
      [--fraction-min 10] [--reads-min 1000] [--host-flag-reads 500] \
      --out results/human_checkpoint

Metadata needs (best effort): a sample-name column, pool size (Num_Ticks), and enough info
to join; controls are auto-detected by name (XB/NTC/BLANK/POS/NEG/EXTRACTION).
"""
import argparse, os, glob, zipfile, tempfile
import pandas as pd
import biom

HOST_TERMS = 'Cervidae|Odocoileus|Alces|Cervus|Ursus|Mammalia|Chordata|Rodentia|Aves|Sciurus|Didelphis'
CONTROL_KEYS = ['XB','NTC','BLANK','POS','NEG','EXTRACTION']

def load_table(qza):
    td=tempfile.mkdtemp()
    with zipfile.ZipFile(qza) as z: z.extractall(td)
    return biom.load_table(glob.glob(td+'/*/data/feature-table.biom')[0]).to_dataframe(dense=True)

def is_control(s):
    return any(k in s.upper() for k in CONTROL_KEYS)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--table', required=True)
    ap.add_argument('--taxonomy', required=True)
    ap.add_argument('--metadata', default=None)
    ap.add_argument('--fraction-min', type=float, default=10.0, help='min %% of sample reads to call a human-bloodmeal candidate')
    ap.add_argument('--reads-min', type=int, default=1000, help='min absolute human reads to call a candidate (filters low-count noise)')
    ap.add_argument('--host-flag-reads', type=int, default=500, help='non-human host reads above which a CONTROL is deemed compromised')
    ap.add_argument('--out', default='results/human_checkpoint')
    a=ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    tax=pd.read_csv(a.taxonomy, sep='\t')
    tax.columns=['Feature ID','Taxon','Confidence'][:len(tax.columns)]; tax=tax.set_index('Feature ID')
    t=load_table(a.table)

    human=[f for f in tax[tax['Taxon'].str.contains('Homo|sapiens',na=False,case=False)].index if f in t.index]
    hostf=[f for f in tax[tax['Taxon'].str.contains(HOST_TERMS,na=False,case=False)
                          & ~tax['Taxon'].str.contains('Homo|sapiens',na=False,case=False)].index if f in t.index]

    htot=t.loc[human].sum(axis=0) if human else pd.Series(0,index=t.columns)
    hostt=t.loc[hostf].sum(axis=0) if hostf else pd.Series(0,index=t.columns)
    stot=t.sum(axis=0)
    hfrac=100*htot/stot.replace(0,1)

    ctrl=[s for s in t.columns if is_control(s)]
    bio=[s for s in t.columns if s not in ctrl]

    # classify controls
    clean_blanks=[]; compromised=[]; posctrl=[]
    for s in ctrl:
        if 'POS' in s.upper(): posctrl.append(s)
        elif hostt[s] > a.host_flag_reads: compromised.append(s)
        else: clean_blanks.append(s)

    # contamination baseline from CLEAN blanks only.
    # Ignore near-empty blanks (<100 total reads) where fraction is meaningless noise.
    cb_ok = [s for s in clean_blanks if stot[s] >= 100]
    base_reads = htot[cb_ok].max() if cb_ok else 0
    base_frac  = hfrac[cb_ok].max() if cb_ok else 0

    # pool sizes from metadata (best effort)
    poolsize={}
    if a.metadata and os.path.exists(a.metadata):
        md=pd.read_csv(a.metadata)
        namecol=next((c for c in md.columns if 'sample' in c.lower() or c.lower() in ('accession','barcode')), md.columns[0])
        sizecol=next((c for c in md.columns if 'tick' in c.lower() and ('num' in c.lower() or 'count' in c.lower())), None)
        if sizecol:
            for _,r in md.iterrows():
                poolsize[str(r[namecol]).strip()]=r[sizecol]
    def psize(s):
        pref=s.split('.')[0]
        for key in (s, pref):
            if key in poolsize:
                try: return int(poolsize[key])
                except: return None
        return None

    # per biological sample
    rows=[]
    for s in bio:
        rows.append({'sample':s,'human_reads':int(htot[s]),'sample_total':int(stot[s]),
                     'human_pct':round(hfrac[s],1),'pool_size':psize(s)})
    df=pd.DataFrame(rows).sort_values('human_pct',ascending=False)

    # candidates: BOTH thresholds
    df['candidate']=(df['human_pct']>=a.fraction_min)&(df['human_reads']>=a.reads_min)
    df.to_csv(f'{a.out}/human_per_sample.tsv',sep='\t',index=False)

    cand=df[df['candidate']]
    singles=df[df['pool_size']==1]; pooled=df[(df['pool_size'].notna())&(df['pool_size']>1)]
    cand_single=cand[cand['pool_size']==1]; cand_pool=cand[(cand['pool_size'].notna())&(cand['pool_size']>1)]

    # ---- report ----
    L=[]
    L.append("# cytb HUMAN checkpoint\n")
    L.append("## QC — handler contamination")
    L.append(f"- Human reads total: {int(htot.sum())} ({100*htot.sum()/stot.sum():.3f}% of all cytb reads)")
    L.append(f"- Biological samples with ANY human: {(df['human_reads']>0).sum()} / {len(bio)}")
    L.append(f"- Controls: {len(clean_blanks)} clean blanks, {len(compromised)} COMPROMISED (host DNA), {len(posctrl)} positive")
    L.append(f"- Clean-blank human baseline: max {int(base_reads)} reads / {base_frac:.1f}% of sample")
    if len(clean_blanks) < 3:
        L.append(f"  ** WARNING: only {len(clean_blanks)} clean blanks — contamination baseline is weak. **")
    L.append("")
    L.append("## Human blood-meal CANDIDATES")
    L.append(f"Threshold: human >= {a.fraction_min}% of reads AND >= {a.reads_min} human reads (both required).")
    L.append(f"- Candidates: {len(cand)} total | {len(cand_single)} single ticks | {len(cand_pool)} pooled")
    if len(singles):
        L.append(f"- Single-tick human-feeding rate: {len(cand_single)}/{len(singles)} = {100*len(cand_single)/len(singles):.1f}%")
    if len(pooled):
        L.append(f"- Pooled pool-positivity (human): {len(cand_pool)}/{len(pooled)} = {100*len(cand_pool)/len(pooled):.1f}%")
    L.append("")
    L.append("### Candidate samples")
    for _,r in cand.iterrows():
        L.append(f"- {r['sample']}: {r['human_reads']} human reads = {r['human_pct']}% (pool={r['pool_size']})")
    L.append("")
    L.append("## MANDATORY CAVEAT")
    L.append("Human cytb cannot distinguish a genuine human blood meal from handler/lab contamination")
    L.append("(identical sequence). Candidates are samples where human is a large fraction of a")
    L.append("well-populated sample (less consistent with trace contamination), but contamination")
    L.append("cannot be fully excluded — especially with few clean blanks. Report cautiously.")
    open(f'{a.out}/human_checkpoint_report.md','w').write("\n".join(L))
    print("\n".join(L))
    print(f"\nWrote {a.out}/human_checkpoint_report.md and human_per_sample.tsv")

if __name__=='__main__':
    main()
