#!/usr/bin/env python3
"""
make_validation_numbers.py
Reads the actual analysis files on ron and writes a NUMBERS-ONLY memo (markdown) for Joe & Lawrence
to review. Four tables, each with a one-line plain header saying what it shows. No recommendation,
no interpretation paragraphs -- the numbers, laid out, so they can decide.

Also prints everything to screen so you can eyeball before sending.
Run from project root:  python scripts/make_validation_numbers.py
Output: results/validation/pooltestr_validation_numbers.md
"""
import os, sys, re
import pandas as pd, numpy as np
from scipy import stats

P_AT     = 'results/analysis/analysis_table_16S.tsv'
P_SINGLE = 'results/single_tick_v2/single_tick_v2_full.tsv'
P_XWALK  = 'results/metadata/crosswalk_VERIFIED.tsv'
P_MLE    = 'results/analysis/prevalence_overall_16S.tsv'
P_BIOM   = 'results/decontam/16S/table_decontam_16S.biom'
OUT      = 'results/validation/pooltestr_validation_numbers.md'
PATHS    = ['Anaplasma', 'Borreliella', 'Borrelia', 'Bartonella']
FLOORS   = [0, 10000, 20000, 40000]

def die(m): print(f'[STOP] {m}'); sys.exit(1)
for pth in [P_AT, P_SINGLE, P_XWALK, P_MLE, P_BIOM]:
    if not os.path.exists(pth): die(f'missing input: {pth}')

# ---------- load ----------
at = pd.read_csv(P_AT, sep='\t')
for c in ['total_reads','NumTicks','Anaplasma_pos','Borreliella_pos','Borrelia_pos','Bartonella_pos',
          'Anaplasma','Borreliella','Borrelia','Bartonella','Rickettsia_reads','Francisella_reads']:
    if c in at.columns: at[c] = pd.to_numeric(at[c], errors='coerce')
mle = pd.read_csv(P_MLE, sep='\t').set_index('pathogen')
st  = pd.read_csv(P_SINGLE, sep='\t')
xw  = pd.read_csv(P_XWALK, sep='\t')

# ---------- single-tick anchor: join depth, collapse replicates to 103 specimens ----------
import biom
t = biom.load_table(P_BIOM)
depth = pd.Series(np.asarray(t.sum(axis='sample')).ravel(), index=t.ids('sample'),
                  name='total_reads').reset_index().rename(columns={'index':'_16s'})
st['_bc'] = st['Barcode'].astype(str).str.strip()
xw['_bc'] = xw['Barcode'].astype(str).str.strip()
st2 = (st.merge(xw[['_bc','table_sample_16S']], on='_bc', how='left')
         .merge(depth, left_on='table_sample_16S', right_on='_16s', how='left'))
pc = {p: ('16S_'+p) for p in PATHS}
agg = {pc[p]:'max' for p in PATHS}; agg['total_reads']='max'
sp = st2.groupby('_bc', as_index=False).agg(agg)          # 103 specimens, any-replicate-positive, best depth
assert len(sp) == 103, f'expected 103 specimens, got {len(sp)}'

def wilson(k,n):
    lo,hi = stats.binomtest(k,n).proportion_ci(method='wilson'); return 100*lo,100*hi

# ================= TABLE 1: single-tick vs pooled prevalence =================
t1 = []
for p in PATHS:
    k = int((sp[pc[p]]==1).sum()); n = len(sp); pct = 100*k/n; lo,hi = wilson(k,n)
    m = mle.loc[p]
    t1.append((p, f'{pct:.1f}% ({k}/{n})', f'{lo:.1f}-{hi:.1f}',
               f'{m["MLE"]:.2f}%', f'{m["lo"]:.2f}-{m["hi"]:.2f}', f'{m["crude"]:.1f}%'))

# ================= TABLE 2: pool-size positivity =================
bins = [(1,'1 (single)'),(2,'2-3'),(4,'4-6'),(7,'7-12'),(13,'13+')]
rng  = {'1 (single)':(1,1),'2-3':(2,3),'4-6':(4,6),'7-12':(7,12),'13+':(13,9999)}
t2 = []
for p in PATHS:
    row=[p]
    for _,lab in bins:
        lo,hi = rng[lab]; s = at[(at.NumTicks>=lo)&(at.NumTicks<=hi)]
        row.append(f'{100*(s[p+"_pos"]==1).mean():.0f}%')
    t2.append(row)
poolsize_n = [len(at[(at.NumTicks>=rng[lab][0])&(at.NumTicks<=rng[lab][1])]) for _,lab in bins]

# ================= TABLE 3: detection vs read depth =================
at['dq'] = pd.qcut(at.total_reads, 4, labels=['Q1 low','Q2','Q3','Q4 high'])
qmed = [f'{g.total_reads.median():.0f}' for _,g in at.groupby('dq', observed=True)]
t3 = []
for p in PATHS:
    row=[p]
    for _,g in at.groupby('dq', observed=True):
        row.append(f'{100*(g[p+"_pos"]==1).mean():.0f}%')
    t3.append(row)

# ================= TABLE 4: co-occurrence robustness to depth =================
def OR(df,a,b):
    A=((df[a+'_pos']==1)&(df[b+'_pos']==1)).sum(); B=((df[a+'_pos']==1)&(df[b+'_pos']==0)).sum()
    C=((df[a+'_pos']==0)&(df[b+'_pos']==1)).sum(); D=((df[a+'_pos']==0)&(df[b+'_pos']==0)).sum()
    orr,_=stats.fisher_exact([[A,B],[C,D]]); return orr
def MH(df,a,b,strat):
    num=den=0
    for _,g in df.groupby(strat, observed=True):
        A=((g[a+'_pos']==1)&(g[b+'_pos']==1)).sum(); B=((g[a+'_pos']==1)&(g[b+'_pos']==0)).sum()
        C=((g[a+'_pos']==0)&(g[b+'_pos']==1)).sum(); D=((g[a+'_pos']==0)&(g[b+'_pos']==0)).sum()
        n=len(g)
        if n: num+=A*D/n; den+=B*C/n
    return num/den if den else float('nan')
copairs=[('Anaplasma','Borreliella'),('Anaplasma','Borrelia'),('Borreliella','Bartonella')]
t4=[(f'{a} + {b}', f'{OR(at,a,b):.2f}', f'{MH(at,a,b,"dq"):.2f}') for a,b in copairs]

# ---------- anchor stability (goes under table 1) ----------
anc=[]
for p in PATHS:
    row=[p]
    for f in FLOORS:
        g=sp[sp.total_reads.fillna(0)>=f]; k=int((g[pc[p]]==1).sum()); n=len(g)
        row.append(f'{100*k/n:.1f}% (n={n})')
    anc.append(row)

# ================= WRITE =================
os.makedirs(os.path.dirname(OUT), exist_ok=True)
L=[]
L.append('# Tick pooled-prevalence: the numbers\n')
L.append('Generated directly from the analysis files on ron. No interpretation added; '
         'each table has a one-line note on what it shows.\n')
L.append(f'Inputs: `{P_AT}`, `{P_SINGLE}` (collapsed to 103 specimens), `{P_MLE}`, `{P_BIOM}`.\n')

L.append('\n## 1. Per-tick prevalence: single ticks vs pooled estimate\n')
L.append('What it shows: how common each pathogen is, measured two ways. Single ticks are counted directly; '
         'the pooled column is the PoolTestR per-tick estimate back-calculated from the pools.\n')
L.append('| Pathogen | Single-tick % (n=103) | 95% CI | Pooled per-tick est. | 95% CI | Pooled crude pool-positivity |')
L.append('|---|---|---|---|---|---|')
for r in t1: L.append('| '+' | '.join(r)+' |')

L.append('\nSingle-tick prevalence at rising read-depth cutoffs (does the single-tick number stay stable?):\n')
L.append('| Pathogen | '+' | '.join(f'>= {f} reads' for f in FLOORS)+' |')
L.append('|---|'+'---|'*len(FLOORS))
for r in anc: L.append('| '+' | '.join(r)+' |')

L.append('\n## 2. Pooled positivity by pool size\n')
L.append('What it shows: for each pathogen, the % of pools testing positive as pools get bigger. '
         'More ticks per pool means more chances a pool contains a positive, so a real signal should rise with size.\n')
L.append('| Pathogen | '+' | '.join(lab for _,lab in bins)+' |')
L.append('|---|'+'---|'*len(bins))
for r in t2: L.append('| '+' | '.join(r)+' |')
L.append(f'\n(pools per size bin: '+', '.join(f'{lab} n={nn}' for (_,lab),nn in zip(bins,poolsize_n))+')\n')

L.append('\n## 3. Detection by sequencing depth\n')
L.append('What it shows: % of samples positive, split into four groups from lowest to highest sequencing depth.\n')
L.append('| Pathogen | '+' | '.join(f'{q}' for q in ["Q1 low","Q2","Q3","Q4 high"])+' |')
L.append('|---|---|---|---|---|')
for r in t3: L.append('| '+' | '.join(r)+' |')
L.append(f'\n(median reads/sample per group: '+', '.join(f'{q}={m}' for q,m in zip(["Q1","Q2","Q3","Q4"],qmed))+')\n')

L.append('\n## 4. Co-occurrence, before and after adjusting for sequencing depth\n')
L.append('What it shows: how strongly two pathogens occur together (odds ratio), and whether that holds '
         'once sequencing depth is accounted for. Values near or above 2 = tend to co-occur; ~1 = independent.\n')
L.append('| Pathogen pair | Odds ratio (unadjusted) | Odds ratio (depth-adjusted) |')
L.append('|---|---|---|')
for r in t4: L.append('| '+' | '.join(r)+' |')

txt='\n'.join(L)+'\n'
open(OUT,'w').write(txt)
print(txt)
print(f'\n[WROTE] {OUT}')
