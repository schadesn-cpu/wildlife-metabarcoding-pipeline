#!/usr/bin/env python3
"""
build_crosswalk_from_compiled.py
Build the VERIFIED specimen crosswalk from the Compiled metadata file, which (unlike the Rosetta
Stone) HAS the USGSTICKS host mappings filled in and already tracks replicates.

Structure of the Compiled file ('All Run Data - w Metadata' sheet, one row per specimen-replicate):
  col 0  Barcode      (specimen id — COLLAPSES replicates: LHT002-1 & -2 share one Barcode)
  col 1  Num
  col 2  Barcode2
  col 4  Num Ticks    (pool size)
  col 6  Host Spp     (host)
  col 10 Town, 11 State
  Then 4 marker blocks (Run Batch / Sample / Plate / Num / Amp):
  col 14 s16S | col 20 s18S | col 26 scytb | col 32 sITS   (each = that marker's exact sample name)

Join logic: feature-table sample id  --strip run-suffix '.{run}'-->  match the marker's Sample col
(which already includes the '-16s'/'-18s' suffix and replicate '-N').  Read Host/NumTicks/Barcode.
Collapse by Barcode for per-specimen prevalence (a specimen is + if ANY replicate is +).

Outputs results/metadata/crosswalk_VERIFIED.tsv with columns:
  table_sample_16S, table_sample_18S, table_sample_cytb, Barcode, Num, NumTicks, Host, Town, State
and an explicit unmatched list.
"""
import pandas as pd, openpyxl, zipfile, glob, tempfile, os
try:
    import biom
    HAVE_BIOM=True
except ImportError:
    HAVE_BIOM=False

COMPILED = "metadata/tick_full_metadata.xlsx"  # put on ron
SHEET = "All Run Data - w Metadata"
OUT = "results/metadata/crosswalk_VERIFIED.tsv"

def load_compiled():
    wb=openpyxl.load_workbook(COMPILED,data_only=True,read_only=True)
    rows=list(wb[SHEET].iter_rows(values_only=True))
    recs=[]
    for r in rows[1:]:
        recs.append({'Barcode':r[0],'Num':r[1],'NumTicks':r[4],'Host':r[6],
                     'Town':r[10],'State':r[11],
                     's16S':r[14],'s18S':r[20],'scytb':r[26],'sITS':r[32]})
    return pd.DataFrame(recs)

def table_ids(qza):
    td=tempfile.mkdtemp()
    with zipfile.ZipFile(qza) as z: z.extractall(td)
    return list(biom.load_table(glob.glob(td+'/*/data/feature-table.biom')[0]).ids('sample'))

def norm(s):
    return str(s).strip().lower() if s is not None else None

def main():
    d=load_compiled()
    os.makedirs('results/metadata',exist_ok=True)
    # build lookup: normalized marker-sample -> row index, for each marker
    lut={}
    for col in ['s16S','s18S','scytb','sITS']:
        for i,v in d[col].items():
            if v is not None and str(v).strip():
                lut[norm(v)]=i

    report={}
    for marker,qza in [('16S','results/decontam/16S/table_decontam_16S.qza'),
                       ('18S',[x for x in glob.glob('qiime2/18S/all/**/table*.qza',recursive=True) if 'raref' not in x][0] if HAVE_BIOM else None),
                       ('cytb','qiime2/cytb/all/dada2/table.qza')]:
        if not HAVE_BIOM or qza is None or not os.path.exists(qza): 
            report[marker]='(table not found / biom missing)'; continue
        tab=table_ids(qza)
        matched=[]; unmatched=[]
        for t in tab:
            pref=norm(t.split('.')[0])
            # try exact match, then with common marker suffixes (18S names are mostly bare)
            suffs={'16S':['-16s',''],'18S':['','-18s'],'cytb':['-cytob','-cytoB','']}.get(marker,[''])
            hit=lut.get(pref)
            if hit is None:
                for sfx in suffs:
                    hit=lut.get(pref+sfx)
                    if hit is not None: break
            if hit is not None: matched.append((t,hit))
            else: unmatched.append(t.split('.')[0])
        report[marker]=(len(tab),len(matched),len(unmatched))
        # write unmatched
        with open(f'results/metadata/unmatched_{marker}.txt','w') as f:
            for u in sorted(unmatched): f.write(u+'\n')
        # for 16S, write the full crosswalk
        if marker=='16S':
            out=[]
            for t,idx in matched:
                row=d.loc[idx]
                out.append({'table_sample_16S':t,'Barcode':row['Barcode'],'Num':row['Num'],
                            'NumTicks':row['NumTicks'],'Host':row['Host'],
                            'Town':row['Town'],'State':row['State']})
            cw=pd.DataFrame(out)
            cw.to_csv(OUT,sep='\t',index=False)

    print("=== crosswalk build report ===")
    for m,r in report.items():
        if isinstance(r,tuple):
            print(f"  {m}: {r[1]}/{r[0]} matched, {r[2]} unmatched")
        else:
            print(f"  {m}: {r}")
    if os.path.exists(OUT):
        cw=pd.read_csv(OUT,sep='\t')
        print(f"\nwrote {OUT}: {len(cw)} rows, {cw['Host'].notna().sum()} with host")
        print("unique specimens (by Barcode):", cw['Barcode'].nunique())
        print("single ticks (NumTicks==1):", (cw['NumTicks']==1).sum())

if __name__=='__main__':
    main()
