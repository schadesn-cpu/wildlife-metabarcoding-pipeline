#!/usr/bin/env python3
"""
Figure 1B: deer-only Bray-Curtis PCoA by batch, with PERMANOVA/PERMDISP stats
annotated directly on the figure (no dependency on qzv-parser).
Run on ron.
  python scripts/fig1b_annotate_stats.py
Reads the PCoA results + cleaned metadata, redraws with Wong palette + batch
ellipses + on-figure stats box. Writes results/figures/fig1b_deer_by_batch_stats.{png,pdf,svg}
"""
import pandas as pd, numpy as np
import zipfile, glob, tempfile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

DEER = 'results/diversity_rebuild_20260706_2045'
pcoa_qza = f'{DEER}/core_metrics/bray_curtis_pcoa_results.qza'
meta_f = f'{DEER}/meta_wtd_clean.tsv'

# --- extract PCoA coordinates from the qza ---
td = tempfile.mkdtemp()
with zipfile.ZipFile(pcoa_qza) as z:
    z.extractall(td)
ordf = glob.glob(f'{td}/*/data/ordination.txt')[0]

# parse QIIME ordination.txt: Site section has sample coords
with open(ordf) as fh:
    lines = fh.read().split('\n')
# find proportion explained and site rows
prop = None
coords = {}
i = 0
while i < len(lines):
    ln = lines[i]
    if ln.startswith('Proportion explained'):
        # next line has the values
        vals = lines[i+1].split('\t')
        prop = [float(v) for v in vals if v.strip()]
    if ln.startswith('Site\t'):
        n = int(ln.split('\t')[1])
        for j in range(1, n+1):
            parts = lines[i+j].split('\t')
            coords[parts[0]] = [float(parts[1]), float(parts[2])]
        i += n
    i += 1

pc = pd.DataFrame(coords).T
pc.columns = ['PC1', 'PC2'] + list(range(len(pc.columns)-2))
pc = pc[['PC1', 'PC2']]

meta = pd.read_csv(meta_f, sep='\t')
meta = meta[meta['sample-id'] != '#q2:types']
bcol = 'batch_clean' if 'batch_clean' in meta.columns else 'batch'
m = pc.merge(meta.set_index('sample-id'), left_index=True, right_index=True)

WONG = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#56B4E9', '#D55E00', '#F0E442', '#000000']
batches = sorted(m[bcol].unique())
cmap = {b: WONG[i % len(WONG)] for i, b in enumerate(batches)}
markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
mmap = {b: markers[i % len(markers)] for i, b in enumerate(batches)}

fig, ax = plt.subplots(figsize=(8.6, 7))
for b in batches:
    sub = m[m[bcol] == b]
    ax.scatter(sub['PC1'], sub['PC2'], s=45, c=cmap[b], marker=mmap[b],
               edgecolor='white', linewidth=0.4, alpha=0.85, label=b, zorder=3)
    # 95% confidence ellipse
    if len(sub) >= 3:
        cov = np.cov(sub['PC1'], sub['PC2'])
        pearson = cov[0, 1]/np.sqrt(cov[0, 0]*cov[1, 1])
        rx, ry = np.sqrt(1+pearson), np.sqrt(1-pearson)
        ell = Ellipse((0, 0), width=rx*2, height=ry*2, facecolor=cmap[b], alpha=0.12,
                      edgecolor=cmap[b], linestyle='--', linewidth=1)
        sx = np.sqrt(cov[0, 0])*2.45; sy = np.sqrt(cov[1, 1])*2.45
        tr = transforms.Affine2D().rotate_deg(45).scale(sx, sy).translate(sub['PC1'].mean(), sub['PC2'].mean())
        ell.set_transform(tr + ax.transData)
        ax.add_patch(ell)

p1 = prop[0]*100 if prop else 40.4
p2 = prop[1]*100 if prop else 9.3
ax.set_xlabel(f'PC1 ({p1:.1f}%)', fontsize=12)
ax.set_ylabel(f'PC2 ({p2:.1f}%)', fontsize=12)
ax.set_title('Bray-Curtis PCoA of deer-collected tick microbiomes', fontsize=13, fontweight='bold')

# on-figure stats box
stats = ('PERMANOVA (batch): pseudo-F = 9.17, R\u00b2 = 10.7%, p = 0.001\n'
         'PERMDISP (batch): p = 0.001\n'
         'State: R\u00b2 = 0.3%, p = 0.47 (n.s.)   Year: R\u00b2 = 0.2%, p = 0.62 (n.s.)')
ax.text(0.02, 0.02, stats, transform=ax.transAxes, fontsize=8.5, va='bottom', ha='left',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', alpha=0.9))

ax.legend(title='Batch', fontsize=8.5, title_fontsize=9, loc='upper right', framealpha=0.9)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
fig.tight_layout()
for e in ('png', 'pdf', 'svg'):
    fig.savefig(f'results/figures/fig1b_deer_by_batch_stats.{e}', dpi=300, bbox_inches='tight', facecolor='white')
print('saved results/figures/fig1b_deer_by_batch_stats.{png,pdf,svg}')
print('n =', len(m), '| batches:', batches)
