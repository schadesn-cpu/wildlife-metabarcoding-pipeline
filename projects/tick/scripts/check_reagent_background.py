#!/usr/bin/env python3
"""check_reagent_background.py

I ask, for every genus in the core microbiome figure, whether it is as abundant
in the negative controls as it is in the ticks. I write nothing.

WHY NOT BLAST
-------------
BLAST answers whether the classifier named the organism correctly. The claim in
Sections 3.3 and 3.4 is that certain genera come from the kit rather than from
the tick, which is a claim about provenance, and no BLAST result settles it. The
negative controls do: a genus as abundant in an extraction blank as in a tick
came from the reagents, whatever its name turns out to be.

WHICH LIBRARIES COUNT AS WHAT
-----------------------------
    negative   EXTRACTIONBLANK, NTC, XB-      the background estimate
    positive   POS, PC-WELL                   the mock community, held out of
                                              the background and reported
                                              separately as a sanity check
    sample     everything else

Pooling the positives into the background would import the mock community and
make its genera look like contamination.

HOW TO READ IT
--------------
A genus in most negative controls at an abundance at or above its abundance in
the ticks is background. A genus absent from every control is real. Anything
between is printed as unresolved rather than guessed at.

Control depth is printed because a no-template control with 200 reads gives a
noisy relative abundance and should not outvote an extraction blank with
thousands.

WHAT I READ
-----------
    qiime2/16S/all/dada2/table.qza     all libraries, controls included
    the 16S taxonomy, through check_batch_arms
    the two core microbiome tables, for the genera to test

I use the pre-decontam table on purpose. The decontam table has the controls
removed, which is the whole quantity I am after.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_reagent_background.py
"""
import os
import re
import shutil
import sys
import tempfile
import zipfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

TABLE = "qiime2/16S/all/dada2/table.qza"
PANELS = [("results/analysis/core_microbiome_ixodes_deer.tsv", "deer"),
          ("results/analysis/core_microbiome_dermacentor_moose.tsv", "moose")]
NEGATIVE = re.compile(r"EXTRACTIONBLANK|NTC|XB-", re.I)
POSITIVE = re.compile(r"POS|PC-WELL", re.I)
TOP = 12


def suppression():
    """I return check_endosymbiont_suppression, for its shared helpers."""
    try:
        import check_endosymbiont_suppression as module
    except ImportError as problem:
        sys.exit("cannot import check_endosymbiont_suppression (%s). It must "
                 "sit in scripts/ beside me." % problem)
    return module


def open_table():
    """I return the raw feature table as a biom Table, kept sparse."""
    if not os.path.isfile(TABLE):
        sys.exit("missing %s. Run me from the project root." % TABLE)
    try:
        import biom
    except ImportError:
        sys.exit("biom is not importable. Run me inside the qiime2 "
                 "environment.")
    workdir = tempfile.mkdtemp(prefix="raw16s_")
    try:
        with zipfile.ZipFile(TABLE) as archive:
            wanted = [n for n in archive.namelist()
                      if n.endswith("data/feature-table.biom")]
            if not wanted:
                sys.exit("no feature-table.biom inside %s. Is it a table?"
                         % TABLE)
            archive.extract(wanted[0], workdir)
            return biom.load_table(os.path.join(workdir, wanted[0]))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def collapse_to_genus(table, genus_of):
    """I sum the sparse table to genus and only then make it a frame.

    The raw table is about 90,000 features across 900 libraries. Densifying
    that first costs several gigabytes, so I collapse while it is still
    sparse and densify the few thousand genera that come out.
    """
    collapsed = table.collapse(lambda feature, metadata: genus_of[feature],
                               axis="observation", norm=False,
                               include_collapsed_metadata=False)
    return pd.DataFrame(collapsed.matrix_data.toarray().T,
                        index=[str(s) for s in collapsed.ids(axis="sample")],
                        columns=[str(g) for g
                                 in collapsed.ids(axis="observation")])


def wanted_genera():
    """I return {genus: [panels it appears in]} for the figure's genera."""
    found = {}
    for path, panel in PANELS:
        if not os.path.isfile(path):
            sys.exit("missing %s. Run build_core_microbiome.py --write first."
                     % path)
        frame = pd.read_csv(path, sep="\t", index_col=0)
        if "mean_pct" not in frame.columns:
            sys.exit("%s has no 'mean_pct' column; it has %s"
                     % (path, ", ".join(map(str, frame.columns))))
        for name in frame.nlargest(TOP, "mean_pct").index:
            found.setdefault(str(name), []).append(panel)
    return found


def moose_animals(module, table, lookup_module):
    """I return the moose animals carrying Dermacentor."""
    if not os.path.isfile(module.SPECIES):
        sys.exit("missing %s" % module.SPECIES)
    species = pd.read_csv(module.SPECIES, sep="\t", dtype=str)
    for column in ("specimen", "host_bin", "Dermacentor_pos"):
        if column not in species.columns:
            sys.exit("%s has no column '%s'" % (module.SPECIES, column))
    keep = set()
    for specimen, host, derm in zip(species["specimen"], species["host_bin"],
                                    species["Dermacentor_pos"]):
        if "moose" in str(host).lower() and module.truthy(derm):
            animal = module.animal_of(table, lookup_module, specimen)
            if animal:
                keep.add(animal)
    return keep


def summarise(shares, rows, genus):
    """I return (share of libraries holding it, mean relative abundance)."""
    if not rows or genus not in shares.columns:
        return 0.0, 0.0
    part = shares.loc[rows, genus]
    return 100.0 * float((part > 0).mean()), 100.0 * float(part.mean())


def main():
    module = suppression()
    arms = module.arms()
    lookup_table, lookup_module = module.lookup()

    raw = open_table()
    features = [str(f) for f in raw.ids(axis="observation")]
    libraries = [str(s) for s in raw.ids(axis="sample")]
    print("raw 16S table: %d libraries x %d features"
          % (len(libraries), len(features)))

    path, covered, mapping = arms.read_taxonomy(features)
    print("taxonomy: %s (%d of %d features named)"
          % (path, covered, len(features)))
    keep = [f for f in features
            if not any(word in str(mapping.get(f, "")).lower()
                       for word in module.HOST_WORDS)]
    print("dropped %d host, chloroplast and eukaryotic features"
          % (len(features) - len(keep)))
    raw = raw.filter(set(keep), axis="observation", inplace=False)

    genus_of = {}
    for feature in keep:
        genus_of[feature] = arms.genus_label(mapping.get(feature, ""))
    by_genus = collapse_to_genus(raw, genus_of)
    depth = by_genus.sum(axis=1)
    shares = by_genus.div(depth.replace(0, float("nan")), axis=0).fillna(0)
    print("collapsed to %d genera" % by_genus.shape[1])

    negatives = [i for i in by_genus.index
                 if NEGATIVE.search(i) and not POSITIVE.search(i)]
    positives = [i for i in by_genus.index if POSITIVE.search(i)]
    held = set(negatives) | set(positives)
    others = [i for i in by_genus.index if i not in held]
    animals = {i: module.animal_of(lookup_table, lookup_module, i)
               for i in others}
    deer = module.deer_with_ixodes(lookup_table, lookup_module)
    moose = moose_animals(module, lookup_table, lookup_module)
    deer_rows = [i for i in others if animals[i] in deer]
    moose_rows = [i for i in others if animals[i] in moose]

    if not negatives:
        sys.exit("no library name matched %s. Check the control naming."
                 % NEGATIVE.pattern)
    print()
    print("negative controls %d, median depth %.0f reads, range %.0f to %.0f"
          % (len(negatives), depth[negatives].median(),
             depth[negatives].min(), depth[negatives].max()))
    print("positive controls %d. Deer libraries %d, moose libraries %d."
          % (len(positives), len(deer_rows), len(moose_rows)))
    thin = [i for i in negatives if depth[i] < 100]
    if thin:
        print("%d negative controls hold under 100 reads; their abundances "
              "are noise and I count them as present or absent only."
              % len(thin))
    print()

    print("=" * 94)
    print("EACH FIGURE GENUS, IN THE NEGATIVE CONTROLS AND IN THE TICKS")
    print("=" * 94)
    print("  %-24s %-11s %17s %15s %15s"
          % ("genus", "panel", "neg controls", "deer", "moose"))
    print("  %-24s %-11s %17s %15s %15s"
          % ("", "", "in %  mean %", "in %  mean %", "in %  mean %"))
    verdicts = []
    for name, panels in sorted(wanted_genera().items()):
        neg_in, neg_mean = summarise(shares, negatives, name)
        deer_in, deer_mean = summarise(shares, deer_rows, name)
        moose_in, moose_mean = summarise(shares, moose_rows, name)
        print("  %-24s %-11s %8.0f %8.3f %6.0f %8.3f %6.0f %8.3f"
              % (name[:24], "/".join(panels), neg_in, neg_mean,
                 deer_in, deer_mean, moose_in, moose_mean))
        here = moose_mean if panels == ["moose"] else deer_mean
        if neg_in == 0:
            verdicts.append((name, "absent from controls"))
        elif neg_in >= 50 and neg_mean >= here:
            verdicts.append((name, "background"))
        else:
            verdicts.append((name, "unresolved"))
    print()

    for label in ("background", "unresolved", "absent from controls"):
        names = [n for n, v in verdicts if v == label]
        if names:
            print("  %-22s %s" % (label + ":", ", ".join(names)))
    print()
    print("  background    in at least half the negative controls, and no less")
    print("                abundant there than in the ticks")
    print("  unresolved    in some controls but less abundant than in the")
    print("                ticks; read it against the control depth above")
    print("  absent        not one read in any negative control")
    print()

    if positives:
        top = shares.loc[positives].mean().nlargest(8)
        print("  Sanity check. The positive controls are dominated by: %s"
              % ", ".join("%s %.1f%%" % (g, 100 * v) for g, v in top.items()))
        print("  If that does not look like the mock community, my control")
        print("  naming is wrong and nothing above should be trusted.")
        print()


if __name__ == "__main__":
    main()
