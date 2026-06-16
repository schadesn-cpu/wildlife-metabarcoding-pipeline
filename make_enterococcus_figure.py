"""
Enterococcus by season figure for AVIC grant preliminary data.
Run from project root: python scripts/make_enterococcus_figure.py
Styled to match 09c_visualize_diversity.py (Wong palette).
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats as scipy_stats

# === Paths ===
TAXA_TSV = 'results/16S/all/taxonomy_refined/taxonomy_relabund_L6_16S.tsv'
META_TSV = 'metadata/qiime/metadata_16S_cod.tsv'
OUT_DIR  = 'results/16S/figures'
os.makedirs(OUT_DIR, exist_ok=True)

# === Wong 2011 colorblind-safe palette (matches lab visualization standard) ===
WONG_COLORS  = ["#0072B2", "#E69F00", "#009E73", "#CC79A7",
                "#56B4E9", "#D55E00", "#F0E442", "#000000"]
WONG_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

# === Load taxonomy ===
taxa = pd.read_csv(TAXA_TSV, sep='\t')
print(f"Taxa table: {taxa.shape[0]} genera x {taxa.shape[1]-1} samples")

ent_row = taxa[taxa['Genus'] == 'Enterococcus']
assert len(ent_row) == 1, f"Expected 1 Enterococcus row, found {len(ent_row)}"
print("Found Enterococcus row.")

sample_cols = [c for c in taxa.columns if c != 'Genus']
df = pd.DataFrame({
    'sample_id_taxa': sample_cols,
    'enterococcus': ent_row[sample_cols].iloc[0].values
})

# Detect proportion vs percentage
if df['enterococcus'].max() <= 1.5:
    df['enterococcus_pct'] = df['enterococcus'] * 100
    print("Values were 0-1; converted to %.")
else:
    df['enterococcus_pct'] = df['enterococcus']
    print("Values already in %.")

# === Load metadata ===
meta = pd.read_csv(META_TSV, sep='\t', comment='#')
if 'sample-id' in meta.columns:
    meta = meta.rename(columns={'sample-id': 'sample_id'})
print(f"Metadata: {meta.shape[0]} rows")

# === Sample ID matching (direct first, then fallbacks) ===
df['sample_id'] = df['sample_id_taxa']
n_match = df['sample_id'].isin(meta['sample_id']).sum()
print(f"Sample ID direct match: {n_match} of {len(df)}")

if n_match == 0:
    df['sample_id'] = df['sample_id_taxa'].str.replace(r'_S\d+$', '', regex=True)
    n_match = df['sample_id'].isin(meta['sample_id']).sum()
    print(f"Sample ID match after _S strip: {n_match} of {len(df)}")

if n_match == 0:
    df['sample_id'] = df['sample_id_taxa'].str.split('-').str[0]
    n_match = df['sample_id'].isin(meta['sample_id']).sum()
    print(f"Sample ID match after bird-ID-only strip: {n_match} of {len(df)}")

# === Merge ===
df = df.merge(meta[['sample_id', 'Season']], on='sample_id', how='inner')
df = df.dropna(subset=['Season'])
print(f"\nMerged: n={len(df)} samples with both Enterococcus and Season")
print("\nPer-season Enterococcus relative abundance:")
print(df.groupby('Season')['enterococcus_pct'].agg(['count', 'median', 'min', 'max']).round(3))

# === Figure ===
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 11,
    'ytick.labelsize': 10,
})

seasons = ['Breeding', 'Freshwater_Nonbreeding', 'Saltwater']
labels  = ['Breeding\n(freshwater)', 'Non-breeding\n(freshwater)', 'Non-breeding\n(saltwater)']

# Wong palette index-based assignment (matches lab convention)
group_colors  = {g: WONG_COLORS[i]  for i, g in enumerate(seasons)}
group_markers = {g: WONG_MARKERS[i] for i, g in enumerate(seasons)}

# Build per-group data with log-floor for visualization of true zeros
FLOOR = 0.001
group_data, group_n = {}, {}
for g in seasons:
    vals = df[df['Season'] == g]['enterococcus_pct'].values
    group_data[g] = np.maximum(vals, FLOOR)
    group_n[g] = len(vals)

# Kruskal-Wallis omnibus test
valid = [group_data[g] for g in seasons if len(group_data[g]) >= 2]
kw_stat, kw_p = scipy_stats.kruskal(*valid)
if kw_p < 0.001:
    p_label = f"Kruskal-Wallis  H = {kw_stat:.2f},  p < 0.001"
else:
    p_label = f"Kruskal-Wallis  H = {kw_stat:.2f},  p = {kw_p:.3f}"
if kw_p <= 0.05:
    p_label += "  *"

# === Plot ===
fig, ax = plt.subplots(figsize=(8.5, 6.5))
positions = list(range(1, len(seasons) + 1))

# Boxplots (color-keyed, light alpha)
box_data = [group_data[g] for g in seasons]
bp = ax.boxplot(
    box_data, positions=positions, widths=0.45,
    patch_artist=True, notch=False,
    medianprops=dict(color="black", linewidth=2.2),
    whiskerprops=dict(color="#444", linewidth=1.1),
    capprops=dict(color="#444", linewidth=1.1),
    flierprops=dict(marker="", linestyle="none"),
    boxprops=dict(linewidth=1.1),
)
for patch, g in zip(bp["boxes"], seasons):
    patch.set_facecolor(group_colors[g])
    patch.set_alpha(0.35)
    patch.set_edgecolor(group_colors[g])

# Jittered points with wider spread to separate floor values
np.random.seed(42)
for pos, g in zip(positions, seasons):
    vals = group_data[g]
    if len(vals) == 0:
        continue
    jitter = np.random.uniform(-0.20, 0.20, size=len(vals))
    ax.scatter(pos + jitter, vals,
               marker=group_markers[g], color=group_colors[g],
               s=85, zorder=5, edgecolors="white", linewidths=0.8, alpha=0.95)

# X-axis labels with n=
ax.set_xticks(positions)
ax.set_xticklabels([f"{lab}\n(n={group_n[g]})" for lab, g in zip(labels, seasons)])

# Y-axis with extended range for breathing room
ax.set_yscale('log')
ax.set_ylim(FLOOR * 0.5, 300)
ax.set_ylabel('$\\it{Enterococcus}$ relative abundance (%)')

# Title
ax.set_title('Gut $\\it{Enterococcus}$ enrichment in saltwater-resident Common Loons',
             pad=14, fontweight='bold')

# Kruskal-Wallis p-value annotation
ax.text(0.5, 0.97, p_label,
        transform=ax.transAxes, ha='center', va='top',
        fontsize=10, style='italic',
        color="#000000" if kw_p <= 0.05 else "#555")

# Detection floor note
ax.text(0.99, 0.015, 'values <0.001% plotted at floor',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=8, color='gray', style='italic')

# Spines and grid
ax.grid(True, alpha=0.45, axis='y', linestyle='--', linewidth=0.4, zorder=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(0.4, len(seasons) + 0.6)

plt.tight_layout()
out_png = f'{OUT_DIR}/enterococcus_by_season.png'
out_pdf = f'{OUT_DIR}/enterococcus_by_season.pdf'
plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
print(f"\nSaved:\n  {out_png}\n  {out_pdf}")
