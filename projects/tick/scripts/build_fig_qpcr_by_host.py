#!/usr/bin/env python3
"""build_fig_qpcr_by_host.py

I build the qPCR detection figure: each target's detection rate in white-tailed
deer and in moose, with single ticks and pools shown separately.

WHY THIS FORM
-------------
Two marks joined by a line, one panel per host. The paired mark is the point:
a single tick and a pool are different units, one tick against a mean of nine,
and a reader needs to see both without me asserting which is the better
estimate. Joining them makes the pair the visual unit instead of inviting a
target-against-target comparison that the denominators do not support.

Bars would be wrong here for a different reason than in the core microbiome
figure. These are proportions of small denominators, and a bar's area implies a
precision that 1 of 75 does not have.

The moose panel is a column of zeros. That is the result, it agrees with the
metabarcoding in Section 3.4, and showing it is more honest than dropping moose
for having nothing in it. Both panels share one x scale, because they measure
the same quantity in the same units and a reader should be able to compare them
directly.

WHAT IS NOT HERE
----------------
Black bear, opossum, porcupine, other wildlife and the unresolved rows. They
are 89 animals between them across both cohorts, and giving them panels beside
335 deer would misrepresent the study. The one result among them worth stating,
8 of 9 bears positive for Babesia microti, belongs in the caption.

COLOUR
------
The validated pair from the core microbiome figure, blue and orange, both above
3:1 contrast on the surface and separated well beyond the CVD threshold. The
meaning is figure-specific and the legend is always present, so nothing rests
on colour alone.

WHAT I READ
-----------
    results/analysis/qpcr_by_host.tsv

WHAT I WRITE
------------
    results/figures/fig2_qpcr_by_host.png
    results/figures/fig2_qpcr_by_host.pdf

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/build_fig_qpcr_by_host.py
    python3 scripts/build_fig_qpcr_by_host.py --write
"""
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SOURCE = "results/analysis/qpcr_by_host.tsv"
OUT = "results/figures/fig2_qpcr_by_host"
SINGLE = "single ticks, one tick per assay"
POOLS = "pools, one slurry per assay"
PANELS = [("White-tailed deer", "A. White-tailed deer"),
          ("Moose", "B. Moose")]
TARGETS = [
    ("a-phagocytophilum", "$\\it{Anaplasma\\ phagocytophilum}$"),
    ("b-microti", "$\\it{Babesia\\ microti}$"),
    ("b-burgdorferi", "$\\it{Borreliella\\ burgdorferi}$"),
    ("b-miyamotoi", "$\\it{Borrelia\\ miyamotoi}$"),
    ("powassan-virus", "Powassan virus"),
    ("b-mayonii", "$\\it{Borreliella\\ mayonii}$"),
    ("e-muris", "$\\it{Ehrlichia\\ muris}$"),
]
ORANGE = "#eb6834"
BLUE = "#2a78d6"
INK = "#0b0b0b"
MUTED = "#52514e"
RULE = "#d8d7d2"
SURFACE = "#fcfcfb"


def read_source():
    """I return the table, or exit saying what is wrong with it."""
    if not os.path.isfile(SOURCE):
        sys.exit("missing %s. Run check_qpcr_by_host.py --write first."
                 % SOURCE)
    frame = pd.read_csv(SOURCE, sep="\t", dtype=str)
    for column in ["cohort", "host", "n", "ticks"]:
        if column not in frame.columns:
            sys.exit("%s has no '%s' column; it has %s"
                     % (SOURCE, column, ", ".join(map(str, frame.columns))))
    for key, _ in TARGETS:
        if key not in frame.columns or key + "_pct" not in frame.columns:
            sys.exit("%s has no column for target '%s'. Its columns are %s"
                     % (SOURCE, key, ", ".join(map(str, frame.columns))))
    return frame


def row_for(frame, cohort, host):
    """I return the one row for this cohort and host, or exit."""
    part = frame[(frame["cohort"] == cohort) & (frame["host"] == host)]
    if len(part) != 1:
        sys.exit("expected one row for %s / %s, found %d. The cohort labels "
                 "in %s are %s" % (cohort, host, len(part), SOURCE,
                                   "; ".join(sorted(set(frame["cohort"])))))
    return part.iloc[0]


def number(value):
    """I turn a table cell into a float, or exit rather than guess."""
    try:
        return float(value)
    except (TypeError, ValueError):
        sys.exit("cannot read '%s' as a number. Check %s for empty cells."
                 % (value, SOURCE))


def plural(count, word):
    """I return '1 tick' and '75 ticks' without writing 'tick(s)'."""
    try:
        many = float(count) != 1
    except (TypeError, ValueError):
        many = True
    return "%s %s%s" % ("{:,}".format(int(float(count))), word,
                        "s" if many else "")


def draw(axis, frame, host, title, limits, labelled):
    """I draw one host panel and return the legend handles it used.

    Only the leftmost panel carries the target names. Both panels show the
    same seven targets in the same order, so repeating them costs half the
    width of the figure and tells the reader nothing.
    """
    single = row_for(frame, SINGLE, host)
    pools = row_for(frame, POOLS, host)
    handles = {}

    for position, (key, _) in enumerate(reversed(TARGETS)):
        one = number(single[key + "_pct"])
        many = number(pools[key + "_pct"])
        axis.hlines(position, limits[0], max(one, many), color=RULE,
                    linewidth=0.8, zorder=1)
        if one != many:
            axis.hlines(position, min(one, many), max(one, many), color=MUTED,
                        linewidth=1.4, zorder=2)
        handles.setdefault("pools", axis.scatter([many], [position], s=58,
                                                 c=ORANGE, linewidths=0.6,
                                                 edgecolors="white", zorder=3))
        handles.setdefault("single", axis.scatter([one], [position], s=32,
                                                  c=BLUE, linewidths=0.6,
                                                  edgecolors="white", zorder=4))

    axis.set_yticks(range(len(TARGETS)))
    axis.set_yticklabels([label for _, label in reversed(TARGETS)],
                         fontsize=7.8)
    axis.set_xlim(*limits)
    axis.set_ylim(-0.8, len(TARGETS) - 0.2)
    axis.set_title("%s\n%s, %s of %s"
                   % (title, plural(single["n"], "single tick"),
                      plural(pools["n"], "pool"),
                      plural(pools["ticks"], "tick")),
                   fontsize=9.5, loc="left", x=0.0, pad=12, color=INK)
    axis.set_xlabel("Animals with a detection (%)", fontsize=8)
    axis.grid(axis="x", color="#ecebe7", linewidth=0.6, zorder=0)
    axis.set_axisbelow(True)
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_color("#c9c8c3")
    axis.tick_params(axis="x", labelsize=7.6, color="#c9c8c3")
    axis.tick_params(axis="y", length=0)
    if not labelled:
        # Shared axes share a tick formatter, so clearing the labels here
        # would clear them on the labelled panel too. Hide them instead.
        axis.tick_params(axis="y", labelleft=False)
    return handles


def main():
    write = "--write" in sys.argv[1:]
    frame = read_source()

    highest = 0.0
    for host, _ in PANELS:
        for cohort in (SINGLE, POOLS):
            row = row_for(frame, cohort, host)
            for key, _ in TARGETS:
                highest = max(highest, number(row[key + "_pct"]))
    limits = (-1.5, highest * 1.12 + 2)

    figure, axes = plt.subplots(1, 2, figsize=(7.8, 3.4), sharey=True)
    figure.patch.set_facecolor(SURFACE)
    handles = {}
    for position, (axis, (host, title)) in enumerate(zip(axes, PANELS)):
        axis.set_facecolor(SURFACE)
        handles.update(draw(axis, frame, host, title, limits,
                            labelled=position == 0))

    order = [("single", "Single ticks, one tick per assay"),
             ("pools", "Pools, one slurry per assay")]
    figure.legend([handles[k] for k, _ in order if k in handles],
                  [t for k, t in order if k in handles],
                  loc="lower center", ncol=2, frameon=False, fontsize=8,
                  bbox_to_anchor=(0.5, -0.02))
    figure.subplots_adjust(left=0.275, right=0.985, bottom=0.245,
                           top=0.800, wspace=0.10)

    for host, title in PANELS:
        single = row_for(frame, SINGLE, host)
        pools = row_for(frame, POOLS, host)
        found = [name for key, name in TARGETS
                 if number(single[key]) or number(pools[key])]
        print("%s: %s, %s of %s. Detected: %s"
              % (host, plural(single["n"], "single tick"),
                 plural(pools["n"], "pool"), plural(pools["ticks"], "tick"),
                 ", ".join(n.replace("$\\it{", "").replace("}$", "")
                           .replace("\\ ", " ") for n in found) or "nothing"))
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
