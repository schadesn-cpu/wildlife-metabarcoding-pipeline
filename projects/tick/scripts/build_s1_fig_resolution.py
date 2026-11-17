#!/usr/bin/env python3
"""build_s1_fig_resolution.py

I build S1 Fig: where each marker's resolution stops, against the rank the
epidemiological question needs.

WHY THIS FORM
-------------
Each row is a before-and-after on one ordinal scale, which is a dumbbell. A
filled mark is the rank the assay reaches, a hollow mark of the same colour is
the rank the question sits at, and the rule between them is the gap. One hue
in two states rather than two hues, because this is not two categories; it is
one quantity in two positions.

Five of the six rows have a gap of exactly one rank. That uniformity is the
point, and it is also what makes equal tick spacing honest here: an ordinal
scale cannot support a claim that one gap is larger than another, and nothing
in this figure makes one. If a row ever needed a two-step gap the spacing
would start to lie and this form would have to change.

The sixth row has no gap. 18S resolves the piroplasms to species, which is why
Babesia odocoilei is a finding in this study rather than a null, and showing
the one assay that reaches is what keeps the figure from reading as a
complaint about markers.

WHY IT IS SUPPLEMENTARY
-----------------------
It carries no new data. Every fact in it is stated in the Results or the
Discussion, and unlike Figures 1 to 3 it is not generated from a table, so it
encodes claims rather than measurements. That makes it a schematic, and a
schematic supporting the future-directions argument belongs in the supplement.

COLOUR
------
The blue and the greys of the main figures, so the supplement looks like the
paper. One hue in two states adds no adjacent pair to the palette already
validated for build_fig_core_microbiome.py at worst-pair deltaE 24.7, so
nothing here needs revalidating. Every row is labelled on the axis and the
two states are distinguished by fill as well as by position, so nothing rests
on colour.

USAGE
-----
    python3 build_s1_fig_resolution.py
    python3 build_s1_fig_resolution.py --write
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "S1_Fig_resolution"

BLUE = "#2a78d6"
INK = "#0b0b0b"
MUTED = "#52514e"
RULE = "#d8d7d2"
SURFACE = "#fcfcfb"

RANKS = ["group or\nfamily", "genus", "species", "below\nspecies"]

# marker, what the assay returns, rank reached, rank needed, what the question is
ROWS = [
    ("16S V4-V5\nbacterial pathogens",
     "$\\it{Borreliella}$", 1, 2, "$\\it{B.\\ burgdorferi}$"),
    ("qPCR panel\n$\\it{Anaplasma}$",
     "$\\it{A.\\ phagocytophilum}$", 2, 3, "Ap-ha / Ap-V1"),
    ("NS5 amplicon\nPowassan virus",
     "Powassan virus", 2, 3, "lineage"),
    ("18S\npiroplasms",
     "$\\it{B.\\ odocoilei}$, $\\it{B.\\ microti}$", 2, 2, None),
    ("18S\nticks",
     "$\\it{Dermacentor}$", 1, 2, "species"),
    ("cytb\nblood meal host",
     "cervid group", 0, 1, "$\\it{Alces}$ / $\\it{Odocoileus}$"),
]


def draw(axis):
    """I draw one row per assay and return the legend handles."""
    handles = {}
    for position, row in enumerate(reversed(ROWS)):
        _, reached_label, reached, needed, needed_label = row
        axis.hlines(position, -0.30, len(RANKS) - 1 + 0.30, color=RULE,
                    linewidth=0.8, zorder=1)
        if needed != reached:
            axis.hlines(position, reached, needed, color=MUTED,
                        linewidth=1.6, zorder=2)
            handles.setdefault(
                "needed",
                axis.scatter([needed], [position], s=64, facecolors=SURFACE,
                             edgecolors=BLUE, linewidths=1.6, zorder=3))
            axis.text(needed + 0.1, position - 0.27, needed_label,
                      fontsize=7.0, color=MUTED, va="center", ha="left")
        handles.setdefault(
            "reached",
            axis.scatter([reached], [position], s=58, c=BLUE, linewidths=0.6,
                         edgecolors="white", zorder=4))
        axis.text(reached - 0.1, position - 0.27, reached_label, fontsize=7.0,
                  color=MUTED, va="center", ha="right")

    closed = len(ROWS) - 1 - [r[3] == r[2] for r in ROWS].index(True)
    axis.text(2.14, closed + 0.05, "no gap", fontsize=7.2, color=BLUE,
              style="italic", va="center", ha="left")
    return handles


def main():
    write = "--write" in sys.argv[1:]
    figure, axis = plt.subplots(figsize=(7.4, 4.0))
    figure.patch.set_facecolor(SURFACE)
    axis.set_facecolor(SURFACE)
    handles = draw(axis)

    axis.set_xticks(range(len(RANKS)))
    axis.set_xticklabels(RANKS, fontsize=7.8)
    axis.set_yticks(range(len(ROWS)))
    axis.set_yticklabels([r[0] for r in reversed(ROWS)], fontsize=7.8)
    axis.set_xlim(-0.95, len(RANKS) - 1 + 0.95)
    axis.set_ylim(-0.70, len(ROWS) - 0.35)
    axis.set_xlabel("Taxonomic resolution", fontsize=8.4, color=INK)
    axis.grid(axis="x", color="#ecebe7", linewidth=0.6, zorder=0)
    axis.set_axisbelow(True)
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_color("#c9c8c3")
    axis.tick_params(axis="x", labelsize=7.8, color="#c9c8c3")
    axis.tick_params(axis="y", length=0)

    order = [("reached", "Resolution this marker reaches"),
             ("needed", "Resolution the epidemiological question needs")]
    figure.legend([handles[k] for k, _ in order if k in handles],
                  [t for k, t in order if k in handles],
                  loc="lower center", ncol=2, frameon=False, fontsize=8,
                  bbox_to_anchor=(0.5, -0.01))
    figure.subplots_adjust(left=0.235, right=0.975, bottom=0.235, top=0.975)

    gaps = sum(1 for r in ROWS if r[3] != r[2])
    print("%d rows, %d with a gap, %d closed" % (len(ROWS), gaps,
                                                 len(ROWS) - gaps))
    sizes = sorted(set(r[3] - r[2] for r in ROWS if r[3] != r[2]))
    print("gap sizes present: %s" % sizes)
    if sizes != [1]:
        print("WARNING: a gap is wider than one rank. Equal tick spacing on an")
        print("ordinal axis can no longer be defended; change the form.")
    if not write:
        print()
        print("Nothing written. Rerun with --write to save %s.png and .pdf"
              % OUT)
        return
    figure.savefig(OUT + ".png", dpi=300, facecolor=figure.get_facecolor())
    figure.savefig(OUT + ".pdf", facecolor=figure.get_facecolor())
    print("wrote %s.png and %s.pdf" % (OUT, OUT))


if __name__ == "__main__":
    main()
