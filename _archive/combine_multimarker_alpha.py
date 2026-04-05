#!/usr/bin/env python3
"""
Combine four single-marker Observed Features plots into one 4-panel figure.
Usage: python scripts/combine_multimarker_alpha.py
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pathlib

FIG_DIR = pathlib.Path("results/multimarker/figures")
OUT_DIR = FIG_DIR

panels = [
    ("16S rRNA V4\n(Bacteria)",        FIG_DIR / "multimarker_16S_observed.png"),
    ("MiFish 12S\n(Fish prey)",         FIG_DIR / "multimarker_MiFish_observed.png"),
    ("Cytochrome b\n(Vertebrate prey)", FIG_DIR / "multimarker_cytb_observed.png"),
    ("18S rRNA V9\n(Eukaryotes)",       FIG_DIR / "multimarker_18S_observed.png"),
]

fig, axes = plt.subplots(1, 4, figsize=(18, 5))
fig.patch.set_facecolor('white')

for ax, (label, imgpath) in zip(axes, panels):
    img = mpimg.imread(imgpath)
    ax.imshow(img)
    ax.axis('off')
    ax.set_title(label, fontsize=13, fontweight='bold', pad=8,
                 fontfamily='Arial', color='#1a1a1a')

# Panel labels A B C D
for i, ax in enumerate(axes):
    ax.text(-0.04, 1.04, chr(65+i), transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='right',
            fontfamily='Arial', color='#1a1a1a')

plt.tight_layout(w_pad=1.5)

out_png = OUT_DIR / "multimarker_alpha_observed_4panel_wong.png"
out_svg = OUT_DIR / "multimarker_alpha_observed_4panel_wong.svg"
fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out_svg, bbox_inches='tight', facecolor='white')
print(f"✓ Saved: {out_png}")
print(f"✓ Saved: {out_svg}")
plt.close()
