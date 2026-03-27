#!/usr/bin/env python3
"""
plot_adeno_tree.py
Polished publication-quality phylogenetic tree for loon adenovirus OTUs.
"""
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# ---------------------------------------------------------------------------
# Load and relabel tree
# ---------------------------------------------------------------------------
tree = Phylo.read("results/adenovirus/adeno_tree.nwk", "newick")

label_map = {
    "adeno_OTU1_total9123": "Loon aviadenovirus OTU1  (9,123 reads)",
    "adeno_OTU6_total128":  "Loon aviadenovirus OTU2  (128 reads)",
    "PP319115.1":           "Aviadenovirus sp. YN06  (PP319115.1)  76.8% identity",
    "PP319121.1":           "Aviadenovirus sp. YN10  (PP319121.1)  76.6% identity",
    "PP319055.1":           "Goose adenovirus 4 HeB08  (PP319055.1)  76.8% identity",
    "NC_017979.1":          "Goose adenovirus 4  (NC_017979.1)  — complete genome",
}

for clade in tree.find_clades():
    if clade.name and clade.name in label_map:
        clade.name = label_map[clade.name]

# Hide low-confidence internal node labels (< 0.8)
for clade in tree.find_clades():
    if clade.confidence is not None and clade.confidence < 0.8:
        clade.confidence = None

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
LOON_COLOR   = "#0072B2"   # Wong blue — loon OTUs
REF_COLOR    = "#555555"   # dark grey — reference sequences
HIGHLIGHT_BG = "#DDEEFF"   # light blue highlight box for loon clade
CONF_COLOR   = "#CC0000"   # red for high-confidence bootstrap values

FONT_FAMILY  = "Arial"
FONT_SIZE_TIP    = 10
FONT_SIZE_BOOT   = 8
FIGURE_DPI       = 300

# ---------------------------------------------------------------------------
# Draw tree
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 5))
plt.rcParams["font.family"] = FONT_FAMILY

def label_colors(clade):
    if clade.name and "Loon" in clade.name:
        return LOON_COLOR
    return REF_COLOR

Phylo.draw(
    tree,
    axes=ax,
    do_show=False,
    label_func=lambda x: x.name if x.name else "",
    label_colors=label_colors,
)

# ---------------------------------------------------------------------------
# Style tip labels and bootstrap values
# ---------------------------------------------------------------------------
for text in ax.texts:
    t = text.get_text().strip()
    # Tip labels
    if t in [v for v in label_map.values()]:
        is_loon = "Loon" in t
        text.set_fontsize(FONT_SIZE_TIP)
        text.set_color(LOON_COLOR if is_loon else REF_COLOR)
        if is_loon:
            text.set_fontweight("bold")
    # Bootstrap values — only show if >= 0.8
    else:
        try:
            val = float(t)
            if val >= 0.8:
                text.set_fontsize(FONT_SIZE_BOOT)
                text.set_color(CONF_COLOR)
                text.set_fontweight("bold")
            else:
                text.set_text("")
        except ValueError:
            pass

# ---------------------------------------------------------------------------
# Highlight loon OTU clade with a background box
# ---------------------------------------------------------------------------
# Find y positions of loon tip labels
loon_ys = []
loon_xs = []
for text in ax.texts:
    if "Loon" in text.get_text():
        loon_ys.append(text.get_position()[1])
        loon_xs.append(text.get_position()[0])

if loon_ys:
    y_min = min(loon_ys) - 0.35
    y_max = max(loon_ys) + 0.35
    x_min = 0.0
    x_max = ax.get_xlim()[1] * 1.15
    rect = FancyBboxPatch(
        (x_min, y_min), x_max - x_min, y_max - y_min,
        boxstyle="round,pad=0.02",
        linewidth=1.2,
        edgecolor=LOON_COLOR,
        facecolor=HIGHLIGHT_BG,
        alpha=0.4,
        zorder=0,
        transform=ax.transData,
    )
    ax.add_patch(rect)
    # Label the box
    ax.text(
        x_min + 0.002, y_max - 0.05,
        "Loon aviadenovirus (putative novel species)",
        fontsize=8, color=LOON_COLOR, fontstyle="italic",
        va="top", ha="left", zorder=5,
    )

# ---------------------------------------------------------------------------
# Axes formatting
# ---------------------------------------------------------------------------
ax.set_ylabel("")
ax.set_xlabel("Branch length (substitutions per site)", fontsize=10)
ax.yaxis.set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", left=False)
ax.tick_params(axis="x", labelsize=9)

# Legend
loon_patch = mpatches.Patch(color=LOON_COLOR, label="Loon aviadenovirus OTUs (this study)")
ref_patch   = mpatches.Patch(color=REF_COLOR,  label="Reference aviadenoviruses (NCBI)")
ax.legend(
    handles=[loon_patch, ref_patch],
    fontsize=9, framealpha=0.9, edgecolor="#CCCCCC",
    loc="lower right",
)

fig.tight_layout()
fig.savefig("results/adenovirus/adeno_phylo_tree.png", dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
fig.savefig("results/adenovirus/adeno_phylo_tree.svg", bbox_inches="tight", facecolor="white")
print("Saved: adeno_phylo_tree.png + .svg")
