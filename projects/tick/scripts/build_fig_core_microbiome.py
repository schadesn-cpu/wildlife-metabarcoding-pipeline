#!/usr/bin/env python3
"""build_fig_core_microbiome.py

I build the core microbiome figure: the most abundant genera in each tick and
host pair, ranked, with prevalence beside them.

WHY THIS FORM
-------------
Abundance runs from 84% down to 0.2%, which rules out bars. A linear bar makes
Rickettsia the only visible bar and hides everything the section is about, and
a log bar chart is worse, because a bar encodes length from zero and a log axis
has no zero. A dot encodes position, so a log axis is honest for it.

I tried a prevalence-against-abundance scatter of the whole community first. It
failed on real data for reasons worth recording: with 11 moose, prevalence can
only take eleven values, so the moose panel became a comb that showed the
denominator rather than the biology, and a shared x axis implied the two panels
were comparable when one percent is three deer and nine percent is one moose.

Prevalence is printed rather than plotted. It is a second measure on a
different scale, and two scales on one chart is the mistake this avoids.

COLOUR AND SHAPE
----------------
Two hues and a neutral, validated for colour-vision deficiency before use
(worst adjacent pair deltaE 24.7, both above 3:1 contrast on the surface):

    orange   the endosymbiont of that panel
    blue     a genus containing recognised tick-borne human pathogens
    grey     every other genus

A hollow mark means the genus is reagent background: enriched in the negative
controls relative to its share of the sequencing pool, by 54-fold for
Escherichia-Shigella and 14-fold for Bacillus. They stay in the figure because
the rule is the twelve most abundant genera and quietly dropping the
inconvenient ones would be selection. Shape carries this rather than a fourth
hue, so the palette stays at three and the distinction survives greyscale.

Every row is labelled, so identity never rests on colour.

WHAT I READ
-----------
    results/analysis/core_microbiome_ixodes_deer.tsv
    results/analysis/core_microbiome_dermacentor_moose.tsv

WHAT I WRITE
------------
    results/figures/fig3_core_microbiome.png
    results/figures/fig3_core_microbiome.pdf

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/build_fig_core_microbiome.py
    python3 scripts/build_fig_core_microbiome.py --write
"""
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PANELS = [
    ("results/analysis/core_microbiome_ixodes_deer.tsv",
     "A. $\\it{Ixodes}$ from white-tailed deer", "Rickettsia"),
    ("results/analysis/core_microbiome_dermacentor_moose.tsv",
     "B. $\\it{Dermacentor}$ from moose", "Francisella"),
]
OUT = "results/figures/fig3_core_microbiome"
# Bartonella and Rickettsiella were removed 2026-10-13: the taxon decision
# table calls the first ked-associated and the second an endosymbiont, and
# Section 3.3 reports Bartonella outside the pathogen results. The three
# endosymbiont genera that remain each contain a human pathogen elsewhere
# in the genus, which is the point Section 4 makes about genus-level calls.
PATHOGEN = ("Anaplasma", "Borreliella", "Borrelia", "Ehrlichia", "Babesia",
            "Coxiella", "Francisella", "Rickettsia")
BACKGROUND = ("Escherichia-Shigella", "Bacillus")
TOP = 12
ORANGE = "#eb6834"
BLUE = "#2a78d6"
GREY = "#9a9994"
INK = "#0b0b0b"
MUTED = "#52514e"
RULE = "#d8d7d2"


def read_panel(path):
    """I return the genus table, or exit saying what is wrong."""
    if not os.path.isfile(path):
        sys.exit("missing %s. Run build_core_microbiome.py --write first."
                 % path)
    frame = pd.read_csv(path, sep="\t", index_col=0)
    for column in ("animals", "animals_with", "prevalence_pct", "mean_pct"):
        if column not in frame.columns:
            sys.exit("%s has no '%s' column; it has %s"
                     % (path, column, ", ".join(map(str, frame.columns))))
    frame = frame[frame["mean_pct"] > 0]
    if frame.empty:
        sys.exit("%s holds no genus with abundance above zero." % path)
    return frame


def label_for(name):
    """I italicise a genus but leave a family or a hyphenated bin roman."""
    text = str(name)
    if "(" in text or "-" in text or " " in text:
        return text
    return "$\\it{%s}$" % text


def draw(axis, frame, title, endosymbiont, limits):
    """I draw one ranked panel and return the legend handles it used."""
    animals = int(frame["animals"].iloc[0])
    shown = frame.nlargest(TOP, "mean_pct").iloc[::-1]
    handles = {}

    for position, (name, row) in enumerate(shown.iterrows()):
        if name in BACKGROUND:
            kind, colour, size = "background", GREY, 38
        elif name == endosymbiont:
            kind, colour, size = "endosymbiont", ORANGE, 62
        elif name in PATHOGEN:
            kind, colour, size = "pathogen", BLUE, 44
        else:
            kind, colour, size = "other", GREY, 34
        axis.hlines(position, limits[0], row["mean_pct"], color=RULE,
                    linewidth=0.8, zorder=1)
        if kind == "background":
            handle = axis.scatter([row["mean_pct"]], [position], s=size,
                                  facecolors="none", edgecolors=colour,
                                  linewidths=1.3, zorder=3)
        else:
            handle = axis.scatter([row["mean_pct"]], [position], s=size,
                                  c=colour, linewidths=0.6,
                                  edgecolors="white", zorder=3)
        handles.setdefault(kind, handle)
        axis.text(limits[1] * 0.78, position, "%.0f%%" % row["prevalence_pct"],
                  ha="right", va="center", fontsize=7.4, color=MUTED)

    axis.set_yticks(range(len(shown)))
    axis.set_yticklabels([label_for(n) for n in shown.index], fontsize=7.8)
    axis.set_xscale("log")
    axis.set_xlim(*limits)
    axis.set_ylim(-0.8, len(shown) - 0.2)
    axis.set_title("%s (n = %d)" % (title, animals), fontsize=9.5,
                   loc="left", x=-0.42, pad=14)
    axis.set_xlabel("Mean relative abundance (%)", fontsize=8)
    axis.set_xticks([0.1, 1, 10, 100])
    axis.set_xticklabels(["0.1", "1", "10", "100"])
    axis.text(limits[1] * 0.78, len(shown) - 0.45, "detected in", ha="right",
              va="center", fontsize=7.4, color=MUTED, style="italic")
    axis.axvline(limits[1] * 0.28, color=RULE, linewidth=0.7, zorder=0)
    axis.grid(axis="x", color="#ecebe7", linewidth=0.6, zorder=0)
    axis.set_axisbelow(True)
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_color("#c9c8c3")
    axis.tick_params(axis="x", labelsize=7.6, color="#c9c8c3")
    axis.tick_params(axis="y", length=0)
    return handles


def main():
    write = "--write" in sys.argv[1:]
    frames = [(read_panel(path), title, endo) for path, title, endo in PANELS]
    smallest = min(float(f.nlargest(TOP, "mean_pct")["mean_pct"].min())
                   for f, _, _ in frames)
    limits = (smallest * 0.45, 700)

    figure, axes = plt.subplots(1, 2, figsize=(7.8, 4.2))
    figure.patch.set_facecolor("#fcfcfb")
    handles = {}
    for axis, (frame, title, endo) in zip(axes, frames):
        axis.set_facecolor("#fcfcfb")
        handles.update(draw(axis, frame, title, endo, limits))

    order = [("endosymbiont", "Endosymbiont of this tick"),
             ("pathogen", "Genus containing tick-borne human pathogens"),
             ("other", "Other genera"),
             ("background", "Reagent background, enriched in the negative "
                            "controls")]
    figure.legend([handles[k] for k, _ in order if k in handles],
                  [t for k, t in order if k in handles],
                  loc="lower center", ncol=2, frameon=False, fontsize=8,
                  bbox_to_anchor=(0.5, 0.005))
    figure.subplots_adjust(left=0.215, right=0.985, bottom=0.255,
                           top=0.885, wspace=0.80)

    for frame, title, _ in frames:
        shown = frame.nlargest(TOP, "mean_pct")
        print("%s: %d genera, showing top %d from %.2f%% down to %.3f%%"
              % (title.split(".")[1].strip(), len(frame), len(shown),
                 shown["mean_pct"].max(), shown["mean_pct"].min()))
    if not write:
        print()
        print("Nothing written. Rerun with --write to save %s.png and .pdf"
              % OUT)
        return
    for extension in ("png", "pdf"):
        path = "%s.%s" % (OUT, extension)
        if os.path.exists(path):
            sys.exit("%s already exists. Move it aside first." % path)
    figure.savefig(OUT + ".png", dpi=300, facecolor=figure.get_facecolor())
    figure.savefig(OUT + ".pdf", facecolor=figure.get_facecolor())
    print("wrote %s.png and %s.pdf" % (OUT, OUT))


if __name__ == "__main__":
    main()
