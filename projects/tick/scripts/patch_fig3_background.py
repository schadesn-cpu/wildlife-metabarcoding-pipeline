#!/usr/bin/env python3
"""patch_fig3_background.py

I edit build_fig_core_microbiome.py in place so the two genera the negative
controls identified as reagent background are drawn as hollow marks, with a
legend entry saying so. I change nothing else and I refuse rather than guess.

WHY A PATCH
-----------
The figure script is 190 lines and retyping it to change four places invites a
transcription error in the 186 lines that are already correct. Every
replacement below is an exact string match and I stop if one is missing, so a
file that is not the version I expect fails loudly instead of being half
edited.

WHAT CHANGES
------------
    a BACKGROUND tuple naming the two genera
    hollow marks for them in draw()
    a fourth legend entry, on two rows instead of one
    a taller figure and a deeper bottom margin to fit that legend
    the docstring, to say what a hollow mark means

Escherichia-Shigella and Bacillus stay in the figure. The rule is the twelve
most abundant genera, and dropping the ones that turned out to be contaminants
would be selection.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/patch_fig3_background.py
    python3 scripts/patch_fig3_background.py --write
"""
import os
import sys

TARGET = "scripts/build_fig_core_microbiome.py"
EXPECT = "0e94df015ea905dfb2d2c4b92571679f"

EDITS = [
("""COLOUR
------
Two hues and a neutral, validated for colour-vision deficiency before use""",
 """COLOUR AND SHAPE
----------------
Two hues and a neutral, validated for colour-vision deficiency before use"""),

("""    grey     every other genus

Every row is labelled, so identity never rests on colour.""",
 """    grey     every other genus

A hollow mark means the genus is reagent background: enriched in the negative
controls relative to its share of the sequencing pool, by 54-fold for
Escherichia-Shigella and 14-fold for Bacillus. They stay in the figure because
the rule is the twelve most abundant genera and quietly dropping the
inconvenient ones would be selection. Shape carries this rather than a fourth
hue, so the palette stays at three and the distinction survives greyscale.

Every row is labelled, so identity never rests on colour."""),

("""            "Rickettsiella")
TOP = 12""",
 """            "Rickettsiella")
BACKGROUND = ("Escherichia-Shigella", "Bacillus")
TOP = 12"""),

("""        if name == endosymbiont:
            kind, colour, size = "endosymbiont", ORANGE, 62
        elif name in PATHOGEN:
            kind, colour, size = "pathogen", BLUE, 44
        else:
            kind, colour, size = "other", GREY, 34
        axis.hlines(position, limits[0], row["mean_pct"], color=RULE,
                    linewidth=0.8, zorder=1)
        handle = axis.scatter([row["mean_pct"]], [position], s=size, c=colour,
                              linewidths=0.6, edgecolors="white", zorder=3)
        handles.setdefault(kind, handle)""",
 """        if name in BACKGROUND:
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
        handles.setdefault(kind, handle)"""),

("""    figure, axes = plt.subplots(1, 2, figsize=(7.8, 3.9))""",
 """    figure, axes = plt.subplots(1, 2, figsize=(7.8, 4.2))"""),

("""             ("other", "Other genera")]""",
 """             ("other", "Other genera"),
             ("background", "Reagent background, enriched in the negative "
                            "controls")]"""),

("""                  loc="lower center", ncol=3, frameon=False, fontsize=8,
                  bbox_to_anchor=(0.5, -0.015))
    figure.subplots_adjust(left=0.215, right=0.985, bottom=0.21,
                           top=0.87, wspace=0.80)""",
 """                  loc="lower center", ncol=2, frameon=False, fontsize=8,
                  bbox_to_anchor=(0.5, 0.005))
    figure.subplots_adjust(left=0.215, right=0.985, bottom=0.255,
                           top=0.885, wspace=0.80)"""),
]


def main():
    write = "--write" in sys.argv[1:]
    if not os.path.isfile(TARGET):
        sys.exit("missing %s. Run me from the project root." % TARGET)
    text = open(TARGET, encoding="utf-8").read()

    for number, (old, new) in enumerate(EDITS, start=1):
        found = text.count(old)
        if found != 1:
            sys.exit("edit %d matches %d places in %s, not 1. Either the file "
                     "has already been patched or it is not the version I "
                     "expect. Nothing written." % (number, found, TARGET))
        text = text.replace(old, new)
        print("  edit %d: applied" % number)

    print()
    print("all %d edits apply cleanly to %s" % (len(EDITS), TARGET))
    if not write:
        print("Nothing written. Rerun with --write to save the change.")
        return
    open(TARGET, "w", encoding="utf-8").write(text)
    print("wrote %s" % TARGET)
    print("Check it: sed '/^[[:space:]]*$/d' %s | md5sum" % TARGET)
    print("should give %s" % EXPECT)


if __name__ == "__main__":
    main()
