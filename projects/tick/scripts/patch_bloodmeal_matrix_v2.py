#!/usr/bin/env python3
"""Patch tick_bloodmeal_matrix.py: presence/absence per specimen, and stop hiding human.

Two changes, both about what the figure is actually measuring.

1. PRESENCE/ABSENCE, not read counts. The colour is currently row-normalised read
   counts, which is the wrong unit for host detection: 69,485 human reads and 200
   human reads both mean "host DNA present", and the difference between them is
   mostly how much host tissue was in the tube. The cell number is now the count of
   SPECIMENS in which that meal was detected, and the colour is the percentage of
   that host's specimens.

2. HUMAN IS NO LONGER HIDDEN. Human cytb is in 162 specimens and the figure did not
   show it, which was a choice nobody had written down. It becomes a column, with
   the caveat in the footnote: human cytb CANNOT distinguish a blood meal from
   handler or laboratory contamination, human is 0.41% of all cytb reads, and 252 of
   837 biological samples carry trace human. Showing it with the caveat is more
   honest than omitting it.

Run on ron from the project root:
    python3 patch_bloodmeal_matrix_v2.py
"""
import sys

PATH = "scripts/tick_bloodmeal_matrix.py"
source = open(PATH).read()
original_length = len(source)
applied = []


def sub(old, new, label):
    global source
    if source.count(old) != 1:
        sys.exit(f"{label}: found {source.count(old)} matches, expected 1. "
                 f"The file is not what I expected; nothing has been changed.")
    source = source.replace(old, new)
    applied.append(label)


# --- human joins the blood-meal columns -------------------------------------
sub('BM_ORDER = ["Cervidae", "Bird", "Bear", "Opossum", "Raccoon", "Porcupine", "Rabbit"]',
    'BM_ORDER = ["Cervidae", "Bird", "Bear", "Opossum", "Raccoon", "Porcupine",\n'
    '            "Rabbit", "Human"]',
    "Human in BM_ORDER")

sub('BM_LABEL = {"Cervidae":"Cervidae\\n(deer/moose)", "Bird":"Bird",',
    'BM_LABEL = {"Cervidae":"Cervidae\\n(deer/moose)", "Bird":"Bird",\n'
    '            "Human":"Human\\n(see note)",',
    "Human display label")

# --- count specimens, not detections ----------------------------------------
sub("""        hosts=d[host_col].value_counts().index.tolist()
        rows=[]; labels=[]
        for h in hosts:
            idx=d.index[d[host_col]==h]
            sub=bm.loc[idx]
            counts=[int(sub[c].sum()) for c in bm_cols]""",
    """        hosts=d[host_col].value_counts().index.tolist()
        rows=[]; labels=[]
        for h in hosts:
            idx=d.index[d[host_col]==h]
            sub=bm.loc[idx]
            # Presence/absence per specimen. sub[c] is already 0/1 per specimen, so
            # the sum is a specimen count and not a read count; keeping it explicit
            # because the distinction is the whole point of this change.
            counts=[int((sub[c]>0).sum()) for c in bm_cols]""",
    "count specimens with the meal")

# --- colour by percentage of that host's specimens --------------------------
sub("""    rowsum=M.sum(axis=1,keepdims=True); rowsum[rowsum==0]=1
    pct=100*M/rowsum""",
    """    # Percentage of THAT HOST'S specimens in which the meal was detected. I divide
    # by the specimen count parsed out of the row label rather than by the row sum:
    # a specimen with two meals contributes to two cells, so a row can sum to more
    # than 100% and normalising by the row total would silently rescale it away.
    import re as _re
    n_specimens=[]
    for label in rowlab:
        found=_re.search(r"n=(\\d+)", label)
        n_specimens.append(int(found.group(1)) if found else max(int(M[len(n_specimens)].sum()), 1))
    n_specimens=np.array(n_specimens, dtype=float).reshape(-1, 1)
    n_specimens[n_specimens == 0] = 1
    pct=100*M/n_specimens""",
    "colour = % of that host's specimens")

sub('    cbar=fig.colorbar(im,ax=ax,fraction=0.035,pad=0.02); cbar.set_label("% of row",fontsize=8.5)',
    '    cbar=fig.colorbar(im,ax=ax,fraction=0.035,pad=0.02)\n'
    '    cbar.set_label("% of that host\'s specimens",fontsize=8.5)',
    "colourbar label")

# --- footnote ---------------------------------------------------------------
sub("""    note=("Cervidae merges Deer/Moose/Cervidae (cytb can't resolve cervid species). "
          "Diagonal = meal matches collection host; off-diagonal = prior-host meals.")""",
    """    note=("Presence/absence per specimen: the cell is the number of specimens in "
          "which that blood meal was detected and the colour is the percentage of that "
          "host's specimens. Cervidae merges Deer/Moose/Cervidae (cytb cannot resolve "
          "cervid species). Diagonal = meal matches collection host; off-diagonal = "
          "prior-host meals.\\nHuman cytb CANNOT distinguish a blood meal from handler or "
          "laboratory contamination. Human is 0.41% of all cytb reads and 252 of 837 "
          "biological samples carry trace human, so these are human-positive above "
          "baseline, not confirmed human feeding.")""",
    "footnote: presence/absence + human caveat")

open(PATH, "w").write(source)
print(f"{original_length} -> {len(source)}")
for entry in applied:
    print("  +", entry)
print()
print("Now rerun:")
print("  python3 scripts/tick_bloodmeal_matrix.py \\")
print("    --input results/single_tick_v3/single_tick_v3_full.tsv \\")
print("    --host-col host_bin \\")
print("    --out results/figures/fig_bloodmeal_matrix")
print()
print("The single-tick table needs a cytb_Human column for the Human column to")
print("populate. If it does not have one, the column will render empty and that is")
print("the next thing to fix, in build_single_tick_set.py.")
