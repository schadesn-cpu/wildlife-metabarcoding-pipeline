#!/usr/bin/env python3
"""
Anchor figure: single-tick pathogen prevalence + co-infection (n=91).

Source: results/single_tick_v3/single_tick_v3_full.tsv  (CANONICAL, rebuilt 2026-08-10)
        Built by scripts/build_single_tick_set.py from the crosswalk's is_single_tick.
        v2 (n=103) is SUPERSEDED: it included 8 Batch-12 fungal-treatment ticks, which
        are experimental and not wild-caught, and 4 specimens whose extraction record
        contradicts the collection record on tick count.
Output: results/figures/single_tick_prevalence_coinfection.{png,pdf,svg}

KEY DESIGN DECISIONS (interpretive, not cosmetic):
  1. Rickettsia and Francisella are ENDOSYMBIONTS of Ixodes, NOT pathogens.
     They are plotted in a separate, de-emphasized category so no reader
     mistakes 93% Rickettsia for a pathogen signal.
  2. Co-infection is computed on TRUE PATHOGENS ONLY (endosymbionts excluded).
     Computed two ways (per-detection and per-genus, collapsing the two
     Babesia species) so we can see how the Babesia co-detection affects it.
  3. Prevalence shown with Wilson 95% CIs (correct for proportions at small n).

Prints every number to stdout for interpretation and cross-check against
the locked values before the figure is even opened.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import sqrt

SRC = "results/single_tick_v3/single_tick_v3_full.tsv"
OUTBASE = "results/figures/single_tick_prevalence_coinfection"

# Wong colorblind-safe palette
C_PATH = "#0072B2"   # blue: pathogens
C_ENDO = "#999999"   # gray: endosymbionts (de-emphasized)
C_CO   = "#009E73"   # green: co-infection

PATHOGENS = {
    "16S_Anaplasma":          "Anaplasma",
    "16S_Bartonella":         "Bartonella",
    "16S_Borreliella":        "Borreliella",
    "16S_Borrelia":           "Borrelia (miyamotoi)",
    "18S_Babesia_odocoilei":  "Babesia odocoilei",
    "18S_Babesia_microti":    "Babesia microti",
}
ENDOSYMBIONTS = {
    "16S_Rickettsia":  "Rickettsia (endosym.)",
    "16S_Francisella": "Francisella (endosym.)",
}
# Italic display labels (mathtext): genus/species italic, roman qualifiers upright.
# Reviewers flag roman binomials, so every genus/species name is italicized here.
ITALIC_LABEL = {
    "Anaplasma":              r"$\mathit{Anaplasma}$",
    "Bartonella":             r"$\mathit{Bartonella}$",
    "Borreliella":            r"$\mathit{Borreliella}$",
    "Borrelia (miyamotoi)":   r"$\mathit{Borrelia}$ ($\mathit{miyamotoi}$)",
    "Babesia odocoilei":      r"$\mathit{Babesia\ odocoilei}$",
    "Babesia microti":        r"$\mathit{Babesia\ microti}$",
    "Rickettsia (endosym.)":  r"$\mathit{Rickettsia}$ (endosym.)",
    "Francisella (endosym.)": r"$\mathit{Francisella}$ (endosym.)",
}
# genus grouping for co-infection (collapse the two Babesia)
GENUS_GROUPS = {
    "Anaplasma":  ["16S_Anaplasma"],
    "Bartonella": ["16S_Bartonella"],
    "Borreliella":["16S_Borreliella"],
    "Borrelia":   ["16S_Borrelia"],
    "Babesia":    ["18S_Babesia_odocoilei", "18S_Babesia_microti"],
}


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z*z/n
    c = p + z*z/(2*n)
    h = z*sqrt(p*(1-p)/n + z*z/(4*n*n))
    return ((c-h)/d, (c+h)/d)


def main():
    d = pd.read_csv(SRC, sep="\t")
    n = len(d)
    print("n =", n)

    # --- prevalence ---
    def prev_row(col, label):
        k = int((d[col] > 0).sum())
        lo, hi = wilson(k, n)
        return dict(col=col, label=label, k=k, pct=100*k/n,
                    lo=100*lo, hi=100*hi)

    path = [prev_row(c, l) for c, l in PATHOGENS.items()]
    endo = [prev_row(c, l) for c, l in ENDOSYMBIONTS.items()]
    path.sort(key=lambda r: r["k"], reverse=True)
    endo.sort(key=lambda r: r["k"], reverse=True)

    print("\n=== PATHOGEN PREVALENCE (n=%d) ===" % n)
    for r in path:
        print("  %-22s %3d  %5.1f%%  [%4.1f, %4.1f]" %
              (r["label"], r["k"], r["pct"], r["lo"], r["hi"]))
    print("=== ENDOSYMBIONTS (context, not pathogens) ===")
    for r in endo:
        print("  %-22s %3d  %5.1f%%  [%4.1f, %4.1f]" %
              (r["label"], r["k"], r["pct"], r["lo"], r["hi"]))

    # --- co-infection, pathogens only ---
    pcols = list(PATHOGENS.keys())
    per_detection = (d[pcols] > 0).sum(axis=1)   # counts both Babesia sp separately
    genus_mat = pd.DataFrame({g: (d[cols].gt(0).any(axis=1)).astype(int)
                              for g, cols in GENUS_GROUPS.items()})
    per_genus = genus_mat.sum(axis=1)

    def cume(series, label):
        print("\n=== CO-INFECTION (%s) ===" % label)
        for k in range(0, int(series.max())+1):
            c = int((series == k).sum())
            print("  exactly %d: %3d  (%4.1f%%)" % (k, c, 100*c/n))
        for t in (1, 2, 3):
            c = int((series >= t).sum())
            print("  >= %d:      %3d  (%4.1f%%)" % (t, c, 100*c/n))

    cume(per_detection, "per detection, Babesia split")
    cume(per_genus, "per genus, Babesia collapsed")

    # --- pairwise associations (Fisher) on genus matrix ---
    try:
        from scipy.stats import fisher_exact
        cols = list(genus_mat.columns)
        res = []
        for i in range(len(cols)):
            for j in range(i+1, len(cols)):
                a = genus_mat[cols[i]]; b = genus_mat[cols[j]]
                t11=int(((a==1)&(b==1)).sum()); t10=int(((a==1)&(b==0)).sum())
                t01=int(((a==0)&(b==1)).sum()); t00=int(((a==0)&(b==0)).sum())
                if t11+t10==0 or t11+t01==0:
                    continue
                orr, p = fisher_exact([[t11,t10],[t01,t00]])
                res.append((cols[i], cols[j], t11, orr, p))
        res.sort(key=lambda x: x[4])
        m = len(res)
        print("\n=== PAIRWISE CO-OCCURRENCE (Fisher, genus-level) ===")
        print("  pair                         both  OR     p_raw   p_BH")
        for rank,(x,y,both,orr,p) in enumerate(res, 1):
            p_bh = min(p * m / rank, 1.0)
            print("  %-12s x %-12s %3d  %5.2f  %6.3f  %6.3f" %
                  (x, y, both, orr, p, p_bh))
    except Exception as e:
        print("\n(pairwise skipped:", e, ")")

    # --- figure ---
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6),
                                   gridspec_kw={"width_ratios": [1.7, 1]})

    # Panel A: prevalence, pathogens then a gap then endosymbionts
    rows = path + [None] + endo
    ylabels, yvals, colors, errs, ks, his = [], [], [], [], [], []
    y = 0; ypos = []
    for r in rows:
        if r is None:
            y += 1.0          # wider gap between pathogen and endosymbiont blocks
            continue
        ypos.append(y)
        ylabels.append(ITALIC_LABEL.get(r["label"], r["label"]))
        yvals.append(r["pct"])
        colors.append(C_ENDO if "endosym" in r["label"] else C_PATH)
        errs.append([r["pct"]-r["lo"], r["hi"]-r["pct"]])
        ks.append(r["k"])
        his.append(r["hi"])   # upper CI, so labels clear the whisker
        y += 1
    errs = np.array(errs).T
    axA.barh(ypos, yvals, color=colors, height=0.72, zorder=3)
    axA.errorbar(yvals, ypos, xerr=errs, fmt="none", ecolor="#333333",
                 elinewidth=1, capsize=2.5, zorder=4)
    for yy, hi, k in zip(ypos, his, ks):
        axA.text(hi + 2.5, yy, "n=%d" % k, va="center", fontsize=8, color="#555")
    axA.set_yticks(ypos); axA.set_yticklabels(ylabels, fontsize=9)
    axA.invert_yaxis()
    axA.set_xlabel("Prevalence (%%), n=%d" % n, fontsize=9.5)
    axA.set_xlim(0, 108)
    # divider sits in the widened gap between the two blocks
    axA.axhline(len(path)-0.5+0.5, color="#bbb", lw=0.9, ls="--")
    axA.text(106, len(path)-0.5+0.5, "endosymbionts\n(not pathogens)",
             fontsize=6.5, color="#999", va="center", ha="right", style="italic")
    axA.set_title("A. Pathogen prevalence in single ticks", fontsize=10, loc="left")
    for s in ("top","right"): axA.spines[s].set_visible(False)

    # Panel B: co-infection distribution (per-genus)
    maxk = int(per_genus.max())
    counts = [int((per_genus==k).sum()) for k in range(maxk+1)]
    pcts = [100*c/n for c in counts]
    xs = list(range(maxk+1))
    axB.bar(xs, pcts, color=C_CO, width=0.7, zorder=3)
    for x,c,p in zip(xs, counts, pcts):
        axB.text(x, p+1, "%d" % c, ha="center", fontsize=8.5, color="#333")
    axB.set_xticks(xs)
    axB.set_xlabel("Co-infecting pathogen genera per tick", fontsize=9.5)
    axB.set_ylabel("Ticks (%)", fontsize=9.5)
    axB.set_title("B. Co-infection (distinct pathogen genera)", fontsize=10, loc="left")
    atleast1 = 100*int((per_genus>=1).sum())/n
    atleast2 = 100*int((per_genus>=2).sum())/n
    axB.text(0.97, 0.95, "%.1f%% carry \u22651\n%.1f%% carry \u22652" %
             (atleast1, atleast2),
             transform=axB.transAxes, ha="right", va="top", fontsize=8.5,
             bbox=dict(boxstyle="round,pad=0.3", fc="#f2f2f2", ec="none"))
    for s in ("top","right"): axB.spines[s].set_visible(False)

    fig.tight_layout()
    for ext in ("png","pdf","svg"):
        fig.savefig("%s.%s" % (OUTBASE, ext), dpi=200, bbox_inches="tight")
    print("\nwrote %s.{png,pdf,svg}" % OUTBASE)


if __name__ == "__main__":
    main()
