#!/usr/bin/env python3
"""check_cohort_funnel.py

I walk an animal from the collection record down to the analysis cohort and
show where each one drops out. I write nothing.

WHY
---
Table 1 counts 371 deer carrying Ixodes and not Dermacentor. Section 3.3
reports 335. Either those are two definitions of one cohort, which is a
problem, or they are two different quantities, which is fine but has to be
said. This tells us which, and it gives the Methods a sentence that currently
does not exist: how 563 animals became 335.

HOW I FIND AN ANIMAL'S LIBRARIES
--------------------------------
By resolving every library name in the decontaminated table to an animal,
which is what the analysis itself does. My first attempt went the other way,
matching the crosswalk's name_16S against the table's library names, and
every animal dropped, because those two are not written the same way. I print
a few of each below so that difference stays visible rather than being a
thing I worked around.

An animal with several libraries is represented by its deepest, so the read
floor is applied to the deepest library rather than to all of them.

THE STEPS
---------
    the cohort rule      deer with Ixodes and not Dermacentor, or moose with
                         Dermacentor
    in the table         at least one library resolves to that animal
    deep enough          its deepest library is at or above the read floor

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_cohort_funnel.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
SHOW = 3


def suppression():
    """I return check_endosymbiont_suppression, which defines the cohorts."""
    try:
        import check_endosymbiont_suppression as module
    except ImportError as problem:
        sys.exit("cannot import check_endosymbiont_suppression (%s)." % problem)
    return module


def main():
    module = suppression()
    table, finder = module.lookup()

    counts = module.read_counts(table, finder)
    depth = counts.sum(axis=1)
    print("decontaminated table: %d libraries, %d at or above %d reads"
          % (len(counts.index), int((depth >= module.MIN_READS).sum()),
             module.MIN_READS))

    if os.path.isfile(CROSSWALK):
        crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
        theirs = [str(v) for v in crosswalk["name_16S"].dropna()][:SHOW]
        print("  crosswalk name_16S looks like: %s" % ", ".join(theirs))
    print("  table library names look like: %s"
          % ", ".join(str(i) for i in list(counts.index)[:SHOW]))
    print()

    libraries = {}
    for library in counts.index:
        animal = module.animal_of(table, finder, library)
        if animal:
            libraries.setdefault(str(animal), []).append(str(library))
    print("those libraries resolve to %d animals" % len(libraries))

    # The core microbiome works in classified bacterial reads, so the floor
    # is applied after the host, chloroplast and eukaryotic lineages come
    # out, not to the library's total depth.
    arms = module.arms()
    path, covered, mapping = arms.read_taxonomy(list(counts.columns))
    host_like = [f for f in counts.columns
                 if any(word in str(mapping.get(f, "")).lower()
                        for word in module.HOST_WORDS)]
    bacterial = counts.drop(columns=host_like).sum(axis=1)
    print("taxonomy: %s (%d of %d features named)"
          % (path, covered, counts.shape[1]))
    print("dropped %d host, chloroplast and eukaryotic features"
          % len(host_like))
    print()

    deer = module.deer_with_ixodes(table, finder)
    species = pd.read_csv(module.SPECIES, sep="\t", dtype=str)
    moose = set()
    for specimen, host, derm in zip(species["specimen"], species["host_bin"],
                                    species["Dermacentor_pos"]):
        if "moose" in str(host).lower() and module.truthy(derm):
            animal = module.animal_of(table, finder, specimen)
            if animal:
                moose.add(str(animal))

    for label, cohort, section in (("deer with Ixodes, no Dermacentor",
                                    deer, "3.3"),
                                   ("moose with Dermacentor", moose, "3.4")):
        step = {str(a) for a in cohort}
        print("=" * 68)
        print("%s (Section %s)" % (label.upper(), section))
        print("=" * 68)
        print("  %-36s %5d" % ("meet the cohort rule", len(step)))
        present = {a for a in step if a in libraries}
        print("  %-36s %5d   (%d lost)"
              % ("have a library in the table", len(present),
                 len(step) - len(present)))
        chosen = {a: max(libraries[a], key=lambda l: float(depth[l]))
                  for a in present}
        deep = {a for a in present
                if float(depth[chosen[a]]) >= module.MIN_READS}
        print("  %-36s %5d   (%d lost)"
              % ("deepest library at or above %d" % module.MIN_READS,
                 len(deep), len(present) - len(deep)))
        classified = {a for a in deep
                      if float(bacterial[chosen[a]]) >= module.MIN_READS}
        print("  %-36s %5d   (%d lost)"
              % ("and %d bacterial reads" % module.MIN_READS,
                 len(classified), len(deep) - len(classified)))
        several = sum(1 for a in classified if len(libraries[a]) > 1)
        print("  %-36s %5d" % ("of those, with several libraries", several))
        print()

    print("  335 on the deer line and 11 on the moose line means Table 1's")
    print("  371 is the sample set and 335 is the analysis cohort. They are")
    print("  different quantities and Table 1 needs a footnote, not a fix.")
    print()


if __name__ == "__main__":
    main()
