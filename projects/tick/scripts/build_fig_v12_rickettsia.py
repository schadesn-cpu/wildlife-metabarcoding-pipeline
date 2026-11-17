#!/usr/bin/env python3
"""build_fig_vt12_rickettsia.py

One figure for the VT-12 update to C. Sullivan: the endosymbiont's share of
classified bacterial reads in each specimen that carries enough 16S to support
a composition.

WHY THIS FORM AND NOT A STACKED BAR
-----------------------------------
A stacked composition bar was the obvious choice and is the wrong one here.

  * It draws proportions, which hide their denominator. 59TW2 has 73 total 16S
    reads and would render as 53% Borreliella off seventeen of them, at the
    same visual weight as BearC's 486,337. vt12_report.txt says so in terms:
    "Do not read them as community profiles."
  * Two of the four categories it would need are Unassigned and Other, which
    are absences of identity rather than taxa, and they would occupy most of
    the legend.
  * The question the figure answers is a single comparison: Rickettsia
    dominates every adequately sequenced specimen except one. That is one
    series, so no legend is needed and no categorical palette has to be
    validated.

The full composition travels with the tables (vt12_genus_matrix.tsv), which is
where per-genus detail belongs. The figure carries the pattern.

CHOICES THAT ARE NOT COSMETIC
-----------------------------
  * Specimens under 1,000 bacterial reads are excluded, not drawn faint.
    51TL (223 reads) and 59TW2 (73) are named in the caption instead. Drawing
    them would invite exactly the reading the report warns against.
  * Read depth is printed beside each specimen label, so the denominator
    travels with the proportion rather than being a caption footnote.
  * BearA and BearC are marked as bear pools. They are not part of the coded
    set and should not be read as treatment specimens, but they give the
    coded ticks something to sit against.
  * Palette, ink and surface are taken from build_fig_core_microbiome.py so
    this is visibly the same family of figure as the manuscript's.

INPUTS, both regenerated 2026-10-20
  results/analysis/vt12/vt12_genus_matrix.tsv   specimen x genus, relative abundance
  results/analysis/vt12/vt12_report.txt         the 16S reads-per-specimen block
Cross-checked against each other; the script stops if they disagree on which
specimens exist.
"""
import os
import re
import sys

MATRIX = "results/analysis/vt12/vt12_genus_matrix.tsv"
REPORT = "results/analysis/vt12/vt12_report.txt"
OUT = "results/figures_vt12/fig_vt12_rickettsia"
FLOOR = 1000
BEARS = ("BearA", "BearC")

BLUE = "#2a78d6"
GREY = "#9a9994"
INK = "#0b0b0b"
MUTED = "#52514e"
RULE = "#d8d7d2"
SURFACE = "#fcfcfb"


def read_matrix(path):
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        specimens = header[1:]
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if row and row[0].strip() == "Rickettsia":
                vals = {}
                for name, v in zip(specimens, row[1:]):
                    try:
                        vals[name.strip()] = float(v)
                    except ValueError:
                        pass
                return specimens, vals
    sys.exit("no Rickettsia row in %s" % path)


def read_depths(path):
    """The '16S reads per specimen' block, which the report prints before its
    composition table precisely so the table is not read without it."""
    depths, inside = {}, False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "16S reads per specimen" in line:
                inside = True
                continue
            if inside:
                # the thin specimens carry a trailing TOO FEW READS note,
                # so anything after the number is allowed and ignored
                m = re.match(r'^\s{2,}(\S+)\s+([\d,]+)(?:\s+\S.*)?\s*$', line)
                if m:
                    depths[m.group(1)] = int(m.group(2).replace(",", ""))
                elif depths and not line.strip():
                    break
    if not depths:
        sys.exit("could not parse the reads-per-specimen block in %s" % path)
    return depths


def main():
    root = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1
                              else "~/meedlab/tick_project")
    os.chdir(root)
    for need in (MATRIX, REPORT):
        if not os.path.isfile(need):
            sys.exit("missing %s" % need)

    specimens, rick = read_matrix(MATRIX)
    depths = read_depths(REPORT)

    only_matrix = sorted(set(rick) - set(depths))
    only_report = sorted(set(depths) - set(rick))
    if only_matrix or only_report:
        print("the two inputs disagree on which specimens exist:")
        print("  only in the matrix: %s" % (only_matrix or "none"))
        print("  only in the report: %s" % (only_report or "none"))
        sys.exit("stopping rather than drawing a partial set")

    kept = [s for s in rick if depths[s] >= FLOOR]
    dropped = [s for s in rick if depths[s] < FLOOR]
    kept.sort(key=lambda s: rick[s])

    print("drawing %d specimens, %d excluded below %d reads: %s"
          % (len(kept), len(dropped), FLOOR,
             ", ".join("%s (%d)" % (s, depths[s]) for s in dropped) or "none"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    height = 0.60 + 0.46 * len(kept)
    fig, ax = plt.subplots(figsize=(7.6, height))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    y = range(len(kept))
    colours = [GREY if s in BEARS else BLUE for s in kept]
    ax.barh(list(y), [rick[s] for s in kept], height=0.52,
            color=colours, linewidth=0)

    # the read count rides in the tick label itself; drawing it as separate
    # text at negative x collided with the specimen names
    labels = []
    for s in kept:
        name = "%s (bear pool)" % s if s in BEARS else s
        labels.append("%s\n%s reads" % (name, "{:,}".format(depths[s])))
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, color=INK, fontsize=8.5, linespacing=1.5)
    ax.tick_params(axis="y", length=0)

    for i, s in enumerate(kept):
        ax.text(rick[s] + 1.5, i, "%.1f%%" % rick[s], va="center",
                ha="left", fontsize=8.5, color=MUTED)

    low = min(kept, key=lambda s: rick[s])
    ax.annotate("95% Mycobacterium", xy=(rick[low] + 9.0, 0),
                xytext=(64, 0), textcoords="offset points",
                va="center", fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle="-", color=RULE, linewidth=1.2))

    ax.set_xlim(0, 108)
    ax.set_xlabel("Rickettsia, % of classified bacterial reads",
                  fontsize=9, color=MUTED)
    ax.xaxis.set_ticks([0, 25, 50, 75, 100])
    ax.tick_params(axis="x", colors=MUTED, labelsize=8.5, length=3)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(RULE)
    ax.xaxis.grid(True, color=RULE, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.margins(y=0.06)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT + ".png", dpi=300, facecolor=fig.get_facecolor())
    fig.savefig(OUT + ".pdf", facecolor=fig.get_facecolor())
    print("wrote %s.png and .pdf" % OUT)


if __name__ == "__main__":
    main()
