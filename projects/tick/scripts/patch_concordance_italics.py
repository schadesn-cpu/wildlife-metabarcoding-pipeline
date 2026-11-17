#!/usr/bin/env python3
"""patch_concordance_italics.py

I italicise the organism names on the x-axis of fig3_qpcr_concordance.

WHY
---
WRITING_STYLE_GUIDE.md requires genus and species names in italic. fig2 does
this already (its y-tick labels are drawn with fontstyle='italic') and
fig_qpcr_prevalence_coinfection does it with mathtext. This figure draws
"Anaplasma / A. phagocytophilum", "Babesia microti" and "Borreliella /
B. burgdorferi" in roman. Three figures in one manuscript, two conventions.

I use mathtext rather than fontstyle='italic' because these labels are part
organism and part punctuation: the slash between a genus and a species should
stay upright, and fontstyle applies to the whole label.

The TSV is deliberately left alone. results/concordance/concordance_2x2.tsv
keeps plain text; mathtext markup belongs in a figure, not in a data file that
something else may read.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_concordance_italics.py
    python3 scripts/build_concordance.py
"""
import os
import sys

TARGET = "scripts/build_concordance.py"

OLD_PAIRS = '''rep = []
def say(m=""):'''

NEW_PAIRS = '''# Figure-only display labels. WRITING_STYLE_GUIDE.md requires italic genus and
# species names; the slash and the abbreviating period stay upright, which is
# why this is mathtext rather than fontstyle='italic' on the whole label.
# Keyed on the plain label so the TSV keeps plain text.
DISPLAY_LABELS = {
    "Anaplasma / A. phagocytophilum":
        "$\\\\mathit{Anaplasma}$ /\\n$\\\\mathit{A.\\\\ phagocytophilum}$",
    "Babesia microti":
        "$\\\\mathit{Babesia\\\\ microti}$",
    "Borreliella / B. burgdorferi":
        "$\\\\mathit{Borreliella}$ /\\n$\\\\mathit{B.\\\\ burgdorferi}$",
}

rep = []
def say(m=""):'''

OLD_LABELS = '    labels = [r["pathogen"] for r in rows]'
NEW_LABELS = ('    labels = [DISPLAY_LABELS.get(r["pathogen"], r["pathogen"])\n'
              '              for r in rows]\n'
              '    _unstyled = [r["pathogen"] for r in rows\n'
              '                 if r["pathogen"] not in DISPLAY_LABELS]\n'
              '    if _unstyled:\n'
              '        say("  WARN no italic display label for: "\n'
              '            + ", ".join(_unstyled)\n'
              '            + ". Add it to DISPLAY_LABELS; it will plot in roman.")')

REPLACEMENTS = [
    (OLD_PAIRS, NEW_PAIRS, "the say() definition, where DISPLAY_LABELS goes"),
    (OLD_LABELS, NEW_LABELS, "the figure label list"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "DISPLAY_LABELS" in source:
        sys.exit("already patched: " + TARGET + " defines DISPLAY_LABELS. "
                 "Nothing changed.")

    for old, new, label in REPLACEMENTS:
        count = source.count(old)
        if count != 1:
            sys.exit("expected exactly one occurrence of " + label
                     + ", found " + str(count)
                     + ". Patch by hand rather than trusting this.")
        source = source.replace(old, new, 1)

    compile(source, TARGET, "exec")
    with open(TARGET, "w", encoding="utf-8") as handle:
        handle.write(source)
    print("patched " + TARGET)
    print("now run: python3 " + TARGET)
    print("the 2x2 numbers must not move: 16/15/3/55, 4/8/0/77, 17/3/12/57 on n=89")


if __name__ == "__main__":
    main()
