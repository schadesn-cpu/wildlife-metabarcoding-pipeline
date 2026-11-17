#!/usr/bin/env python3
"""
fig_cooccurrence_matrix.py
Pool-level pathogen co-occurrence as a MATRIX heatmap.
Retrofitted from the Listeria ST-1 SNP heatmap styling (annotated imshow, adaptive text colour,
white gridlines, colorbar, caption). Swaps to a DIVERGING colormap on log2(odds ratio):
  red = co-occur MORE than expected (positive), blue = LESS (negative), white = independence.
Significant pairs (q<0.05) marked with * and a bold cell border. Diagonal greyed.
Reads results/analysis/cooccurrence_pairs_16S.tsv + analysis_table_16S.tsv.
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm

pairs = pd.read_csv('results/analysis/cooccurrence_pairs_16S.tsv', sep='\t')
at    = pd.read_csv('results/analysis/analysis_table_16S.tsv', sep='\t')

# pathogen order by prevalence (most prevalent first) — nicer to read
paths = ['Anaplasma','Borreliella','Bartonella','Borrelia']
paths = [p for p in paths if p+'_pos' in at.columns]
n = len(paths)
idx = {p:i for i,p in enumerate(paths)}
npos = {p:int(at[p+'_pos'].sum()) for p in paths}

# build symmetric matrices: log2(OR), OR (for labels), q
L = np.full((n,n), np.nan)   # log2 OR
ORm = np.full((n,n), np.nan)
Q = np.full((n,n), np.nan)
for _,r in pairs.iterrows():
    a,b = r['pathogen_A'], r['pathogen_B']
    if a not in idx or b not in idx: continue
    i,j = idx[a], idx[b]
    orr = r['odds_ratio']
    l = np.log2(orr) if (np.isfinite(orr) and orr>0) else 0.0
    for (x,y) in [(i,j),(j,i)]:
        L[x,y]=l; ORm[x,y]=orr; Q[x,y]=r['q_value']

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
fig = plt.figure(figsize=(6.6,6.4))
gs = GridSpec(1,2, width_ratios=[6.2,0.32], wspace=0.05,
              left=0.22, right=0.88, top=0.86, bottom=0.20)

axh = fig.add_subplot(gs[0,0])
vmax = np.nanmax(np.abs(L)) if np.isfinite(np.nanmax(np.abs(L))) else 2
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
cmap = plt.cm.RdBu_r   # colorblind-safe diverging: blue=neg (less), red=pos (more), white=~independence
im = axh.imshow(np.where(np.isnan(L),0,L), cmap=cmap, norm=norm, aspect='equal')

for i in range(n):
    for j in range(n):
        if i==j:
            axh.add_patch(Rectangle((j-.5,i-.5),1,1, facecolor='#eeeeee', edgecolor='white', lw=1.5, zorder=2))
            continue
        orr = ORm[i,j]; q = Q[i,j]
        if not np.isfinite(orr): continue
        sig = q < 0.05
        star = '*' if sig else ''
        # adaptive text colour (like the listeria script): white on dark cells, dark on light.
        # needed because RdBu_r drives strong-OR cells to deep red where #111 text is low-contrast.
        r_, g_, b_, _ = cmap(norm(L[i,j]))
        lum = 0.299*r_ + 0.587*g_ + 0.114*b_   # perceived luminance
        tcol = '#f7f7f7' if lum < 0.5 else '#111111'
        txt = f"OR {orr:.2f}{star}\nq={q:.1e}"
        axh.text(j, i, txt, ha='center', va='center', fontsize=8.2,
                 color=tcol, zorder=4,
                 fontweight=('bold' if sig else 'normal'))
        if sig:  # bold border on significant cells
            axh.add_patch(Rectangle((j-.5,i-.5),1,1, fill=False, edgecolor='#111111', lw=2.2, zorder=3))

axh.set_xticks(range(n)); axh.set_yticks(range(n))
lab = [f"$\\it{{{p}}}$\nn={npos[p]}" for p in paths]
axh.set_xticklabels(lab, rotation=45, ha='right', fontsize=8.8)
axh.set_yticklabels(lab, fontsize=8.8)
axh.set_xticks(np.arange(-.5,n,1), minor=True)
axh.set_yticks(np.arange(-.5,n,1), minor=True)
axh.grid(which='minor', color='white', linewidth=1.5)
axh.tick_params(which='both', length=0)

axc = fig.add_subplot(gs[0,1])
cb = fig.colorbar(im, cax=axc)
cb.set_label('log\u2082(odds ratio)\n\u2190 less   |   more \u2192 co-occurrence', fontsize=8.2)
cb.ax.tick_params(labelsize=8)

fig.suptitle('Pool-level pathogen co-occurrence (16S)', fontsize=12, y=0.95, x=0.55)
fig.text(0.55, 0.10,
         'Pairwise Fisher exact tests on pathogen presence/absence across 603 pooled samples; '
         'odds ratio > 1 = co-occur more than\nexpected under independence. * = significant after '
         'Benjamini-Hochberg correction (q<0.05). Co-occurrence within a pool reflects shared\n'
         'ecology/vector, not confirmed co-infection of a single tick.',
         ha='center', va='top', fontsize=7.0, color='#444444')

import os; os.makedirs('results/figures', exist_ok=True)
OUT='results/figures/fig_cooccurrence_matrix'
fig.savefig(OUT+'.png', dpi=300, bbox_inches='tight')
fig.savefig(OUT+'.pdf', bbox_inches='tight')
fig.savefig(OUT+'.svg', bbox_inches='tight')
print('wrote', OUT+'.png/.pdf/.svg')
print('\nmatrix (OR, * = q<0.05):')
for i in range(n):
    for j in range(i+1,n):
        if np.isfinite(ORm[i,j]):
            s='*' if Q[i,j]<0.05 else ''
            print(f"  {paths[i]:12s} x {paths[j]:12s}: OR={ORm[i,j]:.2f}{s}  q={Q[i,j]:.1e}")
