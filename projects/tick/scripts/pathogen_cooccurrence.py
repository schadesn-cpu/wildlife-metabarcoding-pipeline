#!/usr/bin/env python3
"""
pathogen_cooccurrence.py
Pool-level pathogen co-occurrence from the verified 16S analysis table.
For each pair of real pathogens: does co-occurrence differ from what independence predicts?
Fisher exact test (2x2), odds ratio, BH-corrected q-values. Writes an edge list for the network fig.

Input : results/analysis/analysis_table_16S.tsv  (has <pathogen>_pos 0/1 columns)
Output: results/analysis/cooccurrence_pairs_16S.tsv  (pair, counts, OR, p, q)

HONEST NOTES:
 - Uses presence/absence per SAMPLE (single tick or pool). A pool co-occurrence means both pathogens
   were detected in the same pool — NOT necessarily the same tick. That's a real caveat (co-detection
   in a pool != co-infection of one tick), documented in the output + for the manuscript.
 - Rickettsia/Francisella excluded (endosymbionts). Ehrlichia dropped if all-zero.
 - Small cell counts -> Fisher (exact), not chi-square. BH correction across all pairs.
"""
import pandas as pd, numpy as np
from itertools import combinations
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

d = pd.read_csv('results/analysis/analysis_table_16S.tsv', sep='\t')

# This list MUST match PATHOGENS in build_analysis_table.py. Nothing enforced
# that until 2026-09-18, which is how Bartonella stayed a pathogen here after
# it had been reclassified in build_single_tick_set.py.
pathogens = ['Anaplasma','Borreliella','Borrelia','Ehrlichia']

# Bartonella is a CO-DETECTION, not a pathogen. Ground truth 1.16: these are
# ruminant-associated and ked-associated organisms reported outside the
# pathogen results, and 6 of the 10 Bartonella-positive single ticks carry deer
# ked reads (OR 23.4, Fisher p = 0.0001). The column is still written by
# build_analysis_table.py; it is simply not tested here.
co_detections = ['Bartonella']

missing = [p for p in pathogens if p + '_pos' not in d.columns]
if missing:
    raise SystemExit(
        "analysis_table_16S.tsv has no _pos column for: " + ", ".join(missing)
        + "\nThis list must match PATHOGENS in build_analysis_table.py. "
          "Rebuild the table or fix the list; do not silently test a "
          "different set.")
for c in co_detections:
    if c + '_pos' not in d.columns:
        print("note: " + c + "_pos is not in the table. It is excluded here by "
              "classification, not by absence; check build_analysis_table.py.")
    else:
        print("excluded as a co-detection, not a pathogen: " + c
              + " (" + str(int(d[c + '_pos'].sum())) + "/" + str(len(d)) + ")")

poscols = {p: p+'_pos' for p in pathogens if p+'_pos' in d.columns}
# drop all-zero pathogens (e.g., Ehrlichia)
poscols = {p:c for p,c in poscols.items() if d[c].sum() > 0}
paths = list(poscols)
print("pathogens with detections:", paths)
print("samples:", len(d), "\n")

rows=[]
for a,b in combinations(paths,2):
    ca, cb = poscols[a], poscols[b]
    both   = int(((d[ca]==1)&(d[cb]==1)).sum())
    a_only = int(((d[ca]==1)&(d[cb]==0)).sum())
    b_only = int(((d[ca]==0)&(d[cb]==1)).sum())
    neither= int(((d[ca]==0)&(d[cb]==0)).sum())
    table = [[both, a_only],[b_only, neither]]
    OR, p = fisher_exact(table, alternative='two-sided')
    # expected co-occurrence under independence
    exp_both = (d[ca].mean()*d[cb].mean())*len(d)
    rows.append({'pathogen_A':a,'pathogen_B':b,'both':both,'A_only':a_only,'B_only':b_only,
                 'neither':neither,'expected_both':round(exp_both,1),
                 'odds_ratio':round(OR,3) if np.isfinite(OR) else np.inf,'p_value':p})

res = pd.DataFrame(rows)
res['q_value'] = multipletests(res['p_value'], method='fdr_bh')[1]
res = res.sort_values('p_value')
res.to_csv('results/analysis/cooccurrence_pairs_16S.tsv', sep='\t', index=False)

print("=== pairwise co-occurrence (Fisher exact, BH-corrected) ===")
print(res[['pathogen_A','pathogen_B','both','expected_both','odds_ratio','p_value','q_value']]
      .to_string(index=False))
print("\nSignificant pairs (q<0.05):")
sig=res[res['q_value']<0.05]
if len(sig):
    for _,r in sig.iterrows():
        direction = "MORE than expected" if r['both']>r['expected_both'] else "LESS than expected"
        print(f"  {r['pathogen_A']} & {r['pathogen_B']}: OR={r['odds_ratio']}, q={r['q_value']:.4f} "
              f"({r['both']} obs vs {r['expected_both']} exp — co-occur {direction})")
else:
    print("  none after BH correction")
print("\nwrote results/analysis/cooccurrence_pairs_16S.tsv")
print("\nCAVEAT: co-occurrence in a POOL = both detected in the same pool, NOT proof of co-infection")
print("in one tick. Single-tick co-infection is the direct measure: n = 93, from")
print("results/single_tick_v3/single_tick_v3_full.tsv, computed by")
print("scripts/recompute_cooccurrence.py. NOT n=91; the cohort has been 93 since the")
print("2026-08-07 exclusion rule (ground truth 1.6).")
print("Bartonella is excluded from BOTH analyses, here and there, for the same reason")
print("(ground truth 1.16). Neither is a pathogen co-occurrence result if it is in.")
