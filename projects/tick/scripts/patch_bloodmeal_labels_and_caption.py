#!/usr/bin/env python3
"""patch_bloodmeal_labels_and_caption.py

I fix three things in tick_bloodmeal_matrix.py, in one pass so the figure is
rebuilt once rather than three times.

1. HOST LABELS. The row labels come straight from the raw 'Host Spp' string,
   so the figure reads "NA Porcupine" where the species figure reads "North
   American porcupine". HOST_LABEL_ALIASES already exists for exactly this
   reason (it folds wtd / deer / Odocoileus virginianus into one row); it just
   never grew past the deer case. I extend it, and I add a report of every
   label that was NOT normalised, so a new spelling is visible on the next run
   instead of appearing as a new row in the figure.

   I deliberately do NOT fold 'Bear' to 'Black bear'. The collection record
   says 'Bear'. The species figure says 'Black bear', which asserts a species
   the record does not give. One of the two is adding information and I am not
   going to guess which is intended; see the note printed at the end of the run.

2. WHICH 89. The title says n=89, and so does the concordance figure, and they
   are not the same 89 ticks. This one is the 93 minus 4 with no recorded
   collection host (EQ02440606, EQ02440618, EQ02440767, EQ02442642). The
   concordance 89 is the 93 minus 4 whose only qPCR record is a pool
   (USGSTICKS-4103, 4108, 4118, 7309). Zero overlap between the two sets of
   four. A reader seeing 89 on two figures will assume one cohort, so the
   caption now says which.

3. CELLS EXCEED THE ROW n. Moose reads (n=1) with a 1 in Cervidae and a 1 in
   Bird; Virginia opossum reads (n=2) with 2 and 1. That is correct, because a
   tick can return more than one blood meal, and the code already accounts for
   it (plot_matrix normalises colour by the specimen count rather than the row
   total, with a comment saying why). The caption never said so, so the
   arithmetic looks broken.

I also point the USAGE example at the v3 table. It currently names
single_tick_v2_full.tsv, which is the n=103-era input.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_bloodmeal_labels_and_caption.py
    python3 scripts/tick_bloodmeal_matrix.py \\
        --input results/single_tick_v3/single_tick_v3_full.tsv \\
        --out results/figures/fig_bloodmeal_matrix
"""
import os
import sys

TARGET = "scripts/tick_bloodmeal_matrix.py"

OLD_ALIASES = '''HOST_LABEL_ALIASES = {
    "wtd": "White-tailed deer",
    "deer": "White-tailed deer",
    "white-tailed deer": "White-tailed deer",
    "odocoileus virginianus": "White-tailed deer",
}'''

NEW_ALIASES = '''HOST_LABEL_ALIASES = {
    "wtd": "White-tailed deer",
    "deer": "White-tailed deer",
    "white-tailed deer": "White-tailed deer",
    "odocoileus virginianus": "White-tailed deer",
    # Abbreviation expansions only. These add no information the collection
    # record does not already give; they make this figure agree with
    # fig_tick_species_by_host, which spells the names out.
    "na porcupine": "North American porcupine",
    "n.a. porcupine": "North American porcupine",
    "porcupine": "North American porcupine",
    "virginia opossum": "Virginia opossum",
    "eastern gray squirrel": "Eastern gray squirrel",
    "moose": "Moose",
    # NOT folded: 'Bear'. The collection record says Bear;
    # fig_tick_species_by_host says Black bear, which asserts a species the
    # record does not give. Resolve that deliberately rather than here.
}'''

OLD_REPORT = '''    _merged = int((_raw.str.lower() == "deer").sum())
    if _merged:
        print("normalised %d 'Deer' label(s) into 'White-tailed deer'" % _merged)'''

NEW_REPORT = '''    _changed = {}
    for _before, _after in zip(_raw, d[a.host_col]):
        if _after is not None and str(_before) != str(_after):
            _changed[str(_before)] = str(_after)
    if _changed:
        print("normalised host labels:")
        for _before in sorted(_changed):
            _n = int((_raw == _before).sum())
            print("    %-28s -> %-28s (%d specimen(s))"
                  % (_before, _changed[_before], _n))
    # Any label NOT normalised is printed too. A new spelling should be visible
    # here on the run that introduces it, not as a surprise extra row in the
    # figure six weeks later.
    _untouched = sorted(set(str(v) for v in _raw
                            if str(v).lower() not in HOST_LABEL_ALIASES
                            and str(v).lower() not in ("nan", "", "none")))
    if _untouched:
        print("host labels used as-is (add to HOST_LABEL_ALIASES if any of "
              "these are the same animal under two names):")
        for _u in _untouched:
            print("    %-28s (%d specimen(s))" % (_u, int((_raw == _u).sum())))'''

OLD_NOTE = '''    note=("Cervidae merges Deer/Moose/Cervidae (cytb can't resolve cervid species). "
          "Diagonal = meal matches collection host; off-diagonal = prior-host meals.")'''

NEW_NOTE = '''    note=("Cervidae merges Deer/Moose/Cervidae (cytb can't resolve cervid species). "
          "Diagonal = meal matches collection host; off-diagonal = prior-host meals. "
          "A tick can return more than one blood meal, so the counts in a row can "
          "sum to more than that row's specimen count; colour is the percentage of "
          "that host's specimens, not of the row.")'''

OLD_NHOST = '''    if n_nohost:
        KED_NOTE += (" %d specimen(s) with no recorded collection host are not shown; "
                     "the matrix plots %d of %d." % (n_nohost, n_plotted, len(d)))'''

NEW_NHOST = '''    if n_nohost:
        KED_NOTE += (" %d specimen(s) with no recorded collection host are not shown; "
                     "the matrix plots %d of %d." % (n_nohost, n_plotted, len(d)))
        # Two figures in this manuscript carry n=89 and they are NOT the same
        # ticks. This one is the sequenced set minus the specimens with no
        # recorded collection host. The qPCR concordance figure is the same
        # sequenced set minus the specimens whose only qPCR record is a pool.
        # Zero overlap between the two excluded sets. Say which, or a reader
        # reasonably assumes one cohort.
        KED_NOTE += (" This is not the same subset as the qPCR concordance figure, "
                     "which also plots %d: that one excludes specimens whose only "
                     "qPCR record is a pooled subsample, while this one excludes "
                     "specimens with no recorded collection host. No specimen is "
                     "in both excluded groups." % n_plotted)'''

OLD_USAGE = '      --input results/single_tick_v2/single_tick_v2_full.tsv \\'
NEW_USAGE = '      --input results/single_tick_v3/single_tick_v3_full.tsv \\'

REPLACEMENTS = [
    (OLD_USAGE, NEW_USAGE, "the USAGE input path"),
    (OLD_ALIASES, NEW_ALIASES, "the HOST_LABEL_ALIASES dict"),
    (OLD_REPORT, NEW_REPORT, "the normalisation report"),
    (OLD_NOTE, NEW_NOTE, "the footnote opening"),
    (OLD_NHOST, NEW_NHOST, "the no-host footnote clause"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "na porcupine" in source:
        sys.exit("already patched: " + TARGET + " has the expanded aliases. "
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
    print()
    print("run:  python3 scripts/tick_bloodmeal_matrix.py \\")
    print("        --input results/single_tick_v3/single_tick_v3_full.tsv \\")
    print("        --out results/figures/fig_bloodmeal_matrix")
    print()
    print("The 79 on the diagonal and the nine rows must not change; only the")
    print("porcupine label and the footnote should move.")
    print()
    print("STILL TO DECIDE, not patched: the collection record says 'Bear' and")
    print("fig_tick_species_by_host says 'Black bear'. One of the two adds a")
    print("species the record does not give. Pick one and I will apply it to")
    print("both figures.")


if __name__ == "__main__":
    main()
