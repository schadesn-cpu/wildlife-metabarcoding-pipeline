#!/usr/bin/env python3
"""
Figure 1A: Shannon diversity by host with Dunn post-hoc significance bars.
Run on ron. Reads shannon_vector.qza + meta_all_recovered.tsv, draws Wong-colored
boxplot + jitter, annotates Kruskal-Wallis and significance brackets for the pairs
found significant by Dunn (Bonferroni).
  python scripts/fig1a_shannon_sig.py
Writes results/figures/fig1a_host_shannon_sig.{png,pdf,svg}
"""
import pandas as pd, numpy as np
import zipfile, glob, tempfile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ALL = 'results/diversity_allsamples_20260706_2052'
qza = f'{ALL}/core_metrics/shannon_vector.qza'
td = tempfile.mkdtemp()
with zipfile.ZipFile(qza) as z:
    z.extractall(td)
s = pd.read_csv(glob.glob(f'{td}/*/data/alpha-diversity.tsv')[0], sep='\t', index_col=0)
s.columns = ['shannon']
meta = pd.read_csv(f'{ALL}/meta_all_recovered.tsv', sep='\t')
meta = meta[meta['sample-id'] != '#q2:types']
m = s.merge(meta.set_index('sample-id'), left_index=True, right_index=True)
m = m[m['host_grp'].isin(['WTD', 'Bear', 'Moose', 'Other wildlife'])].copy()

order = ['WTD', 'Bear', 'Moose', 'Other wildlife']
WONG = {'WTD': '#CC79A7', 'Bear': '#0072B2', 'Moose': '#E69F00', 'Other wildlife': '#009E73'}
data = [m.loc[m['host_grp'] == g, 'shannon'].values for g in order]
ns = {g: len(m[m['host_grp'] == g]) for g in order}

fig, ax = plt.subplots(figsize=(8.4, 6.2))
bp = ax.boxplot(data, positions=range(len(order)), widths=0.6, patch_artist=True,
                showfliers=False, medianprops=dict(color='black', linewidth=2))
for patch, g in zip(bp['boxes'], order):
    patch.set_facecolor(WONG[g]); patch.set_alpha(0.45); patch.set_edgecolor(WONG[g])
for i, g in enumerate(order):
    y = data[i]; x = np.random.normal(i, 0.07, len(y))
    ax.scatter(x, y, s=18, color=WONG[g], edgecolor='white', linewidth=0.4, alpha=0.8, zorder=3)
    ax.text(i, ax.get_ylim()[0], f'n={ns[g]}', ha='center', va='top', fontsize=9, color='#555')

ax.set_xticks(range(len(order))); ax.set_xticklabels(order, fontsize=11)
ax.set_ylabel('Shannon diversity', fontsize=12)
ax.set_title('Shannon diversity by collection host', fontsize=13, fontweight='bold')
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)

# significance brackets (Dunn, Bonferroni): pairs of (i,j,stars)
idx = {g: i for i, g in enumerate(order)}
sig = [('WTD', 'Moose', '***'), ('WTD', 'Other wildlife', '***'),
       ('Bear', 'Moose', '***'), ('Bear', 'Other wildlife', '**')]
ymax = m['shannon'].max()
h = ymax * 0.06
lvl = ymax + h
for k, (a, b, star) in enumerate(sig):
    i, j = sorted([idx[a], idx[b]])
    yy = lvl + k * h * 1.6
    ax.plot([i, i, j, j], [yy, yy + h*0.4, yy + h*0.4, yy], lw=1.1, c='#333')
    ax.text((i + j) / 2, yy + h*0.4, star, ha='center', va='bottom', fontsize=11)
ax.set_ylim(top=lvl + len(sig) * h * 1.6 + h)

ax.text(0.02, 0.97, 'Kruskal-Wallis H=44.1, p=1.4e-09\nDunn post-hoc, Bonferroni',
        transform=ax.transAxes, fontsize=9, style='italic', va='top', color='#444')

fig.tight_layout()
for e in ('png', 'pdf', 'svg'):
    fig.savefig(f'results/figures/fig1a_host_shannon_sig.{e}', dpi=300, bbox_inches='tight', facecolor='white')
print('saved results/figures/fig1a_host_shannon_sig.{png,pdf,svg}')
print('significant pairs marked:', sig)
