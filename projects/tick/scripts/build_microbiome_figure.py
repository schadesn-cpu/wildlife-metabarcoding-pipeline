#!/usr/bin/env python3
"""
build_microbiome_figure.py  —  run ON RON, from the project root:
    cd ~/meedlab/tick_project
    python3 scripts/build_microbiome_figure.py

16S bacterial community by collection host, two panels:
    A  full community
    B  with the two locked endosymbionts (Rickettsia, Francisella) removed

WHY NOT 10_plot_taxonomy.py
---------------------------
That script draws one bar per library. With 745 libraries the result is a
diagnostic strip, not a figure. Keep it for QC; this draws the four-group
comparison the figure is actually making.

TWO BUGS THIS FIXES (both live, 2026-08-11)
-------------------------------------------
1. RANKING. Selecting the top taxa by their mean across the four HOST GROUPS
   dropped Rickettsia out of the top 11 entirely — because it is ~0% in moose,
   which dragged its four-group mean below taxa that are moderate everywhere.
   Rickettsia is the single most abundant organism in the dataset. Taxa are now
   ranked by OVERALL abundance and the group means are plotted.
2. ITALICS. Genus names are italicised; Unassigned, Other, uncultured, family
   names (-aceae) and hyphenated group names (Escherichia-Shigella) are not.

THE HEADLINE, and it is worth its own sentence in the results
-------------------------------------------------------------
Rickettsia per-library relative abundance, mean / median:
    Deer            64.7% / 86.9%   (n=690)
    Bear            69.9% / 81.6%   (n=24)
    Other wildlife   8.4% /  0.0%   (n=17)
    Moose            0.0% /  0.0%   (n=14)   <-- ZERO in all 14
482 of 745 libraries are >50% Rickettsia.

Rickettsia buchneri is an Ixodes scapularis endosymbiont. The moose ticks are
Dermacentor. So the 16S community CROSS-VALIDATES the 18S species call from a
completely independent direction: 18S says the moose carry Dermacentor, and the
endosymbiont profile says the same thing. Pasteurella dominates the moose bar
instead.

NOTE ON DENOMINATORS
--------------------
Groups here are LIBRARIES (Deer 690, Bear 24, Other wildlife 17, Moose 14), not
specimens (506 / 10 / 16 / 12). Composition is a per-library property so that is
correct, but the caption must say so or a reader will try to reconcile these with
the species figure and fail.

Mean-of-proportions weights every library equally regardless of depth. That is
standard and is NOT the same quantity as pooling reads: pooled, Rickettsia is
49.5% of all reads and ~73% after removing host and unclassified. Both are true;
say which one any given number is.

Reads:
    results/16S/all/microbiome/relabund_A_full.tsv
    results/16S/all/microbiome/relabund_B_noendo.tsv
    results/16S/all/microbiome/meta_host.tsv
Writes:
    results/figures/fig_microbiome_by_host.png/.pdf
    results/analysis/microbiome_by_host.tsv
"""
import os
import sys
import re
import pandas as pd

M    = "results/16S/all/microbiome/"
FDIR = "results/figures"
ODIR = "results/analysis"
ORDER = ["Deer", "Moose", "Bear", "Other wildlife"]
TOP_N = 11

# Okabe-Ito, colourblind-safe
WONG = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00",
        "#F0E442", "#7570B3", "#1B9E77", "#A6761D", "#666666"]
GREY = "#DDDDDD"

# Not genera: leave roman.
NO_ITALIC = {"Unassigned", "Other", "uncultured", "Escherichia-Shigella"}

def label(t):
    """Genus names italicise; ranks above genus, placeholders and hyphenated
    group names stay roman."""
    t = str(t)
    if t in NO_ITALIC or re.search(r"(aceae|ales|ceae|bacteria)$", t):
        return t
    # NB: build the escaped string OUTSIDE the f-string. Python < 3.12 rejects a
    # backslash inside an f-string expression, and ron runs 3.9.
    esc = t.replace(" ", "\\ ")
    return "$\\it{" + esc + "}$"

rep = []
def say(m=""):
    print(m); rep.append(str(m))

os.makedirs(FDIR, exist_ok=True)
os.makedirs(ODIR, exist_ok=True)

meta = pd.read_csv(M + "meta_host.tsv", sep="\t").set_index("sample-id")["Host"]

panels = [("relabund_A_full.tsv",   "A. Full community"),
          ("relabund_B_noendo.tsv", "B. Rickettsia and Francisella removed")]
for f, _ in panels:
    if not os.path.exists(M + f):
        sys.exit(f"ERROR: {M+f} not found")

# Not wrapped: an ImportError here means the environment is wrong and the
# traceback names the missing package better than I can.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.5))
tables = {}

for ax, (f, ttl) in zip(axes, panels):
    d = pd.read_csv(M + f, sep="\t", index_col=0)
    d = d[[c for c in d.columns if c in meta.index]]
    grp = d.T.groupby(meta).mean().T.reindex(columns=ORDER)

    # RANK BY OVERALL ABUNDANCE, not by the mean of the group means. Ranking on
    # group means is what dropped Rickettsia: ~0% in moose pulled its four-group
    # average below taxa that are moderate in every group.
    overall = d.mean(axis=1).sort_values(ascending=False)
    top = [t for t in overall.head(TOP_N).index]

    plot = grp.loc[top].copy()
    plot.loc["Other"] = 1 - plot.sum(axis=0)
    tables[ttl] = plot

    bottom = [0.0] * len(ORDER)
    for i, t in enumerate(plot.index):
        v = (100 * plot.loc[t]).values
        ax.bar(range(len(ORDER)), v, bottom=bottom, width=0.72,
               label=label(t), color=GREY if t == "Other" else WONG[i % len(WONG)],
               edgecolor="white", linewidth=0.4)
        bottom = [b + x for b, x in zip(bottom, v)]

    n = {h: int((meta == h).sum()) for h in ORDER}
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([f"{h}\n(n = {n[h]})" for h in ORDER])
    ax.set_ylabel("Mean relative abundance (%)")
    ax.set_ylim(0, 100)
    ax.set_title(ttl, loc="left", fontweight="bold")
    ax.legend(bbox_to_anchor=(0.5, -0.13), loc="upper center", ncol=3,
              frameon=False, fontsize=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

fig.suptitle("16S bacterial community by collection host", fontweight="bold", y=0.98)
fig.tight_layout(rect=[0, 0.02, 1, 0.97])
for ext in ("png", "pdf"):
    p = f"{FDIR}/fig_microbiome_by_host.{ext}"
    fig.savefig(p, dpi=300, bbox_inches="tight")
    say(f"wrote {p}")
plt.close(fig)

out = pd.concat({k: v for k, v in tables.items()}, names=["panel", "taxon"])
out.to_csv(f"{ODIR}/microbiome_by_host.tsv", sep="\t")
say(f"wrote {ODIR}/microbiome_by_host.tsv")

for ttl, t in tables.items():
    say(f"\n=== {ttl} — mean relative abundance (%) ===")
    say((100 * t).round(1).to_string())

say("\nCAPTION MUST STATE:")
say("  - groups are LIBRARIES (690/14/24/17), not specimens (506/10/16/12);")
say("    composition is a per-library property")
say("  - mean-of-proportions, weighting every library equally regardless of depth;")
say("    this is NOT the pooled-read figure (Rickettsia 49.5% of all reads,")
say("    ~73% after removing host and unclassified)")
say("  - DESCRIPTIVE ONLY. No statistic is attached and the group sizes are far")
say("    too unbalanced to support one informally (690 deer vs 14 moose).")
say("    Per the analysis plan, community composition is presence/description;")
say("    statistics live in the WTD-only controlled set.")
say("  - the moose bar is a mean over 14 libraries and should not be read as")
say("    precisely as the deer bar")

with open(f"{ODIR}/microbiome_report.txt", "w") as fh:
    fh.write("\n".join(rep))
say(f"\nreport -> {ODIR}/microbiome_report.txt")
