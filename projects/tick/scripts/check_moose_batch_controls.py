#!/usr/bin/env python3
"""check_moose_batch_controls.py

I ask whether the batches holding the moose ticks also held a negative control,
because a genus absent from every other plate's blank is not cleared. I write
nothing.

WHY
---
check_control_enrichment.py showed Enhydrobacter, Deinococcus, Chryseobacterium
and the other moose genera are absent from the negative controls at expected
counts that make the absence real. That clears them of being in the reagents
we sampled. It does not clear them of a contamination event on the moose plate
itself, because reagent contamination happens per plate, and a blank from
another plate cannot see it.

There are 13 moose libraries. If they sit on one or two plates with no blank of
their own, the whole result rests on nothing and I would rather find that out
here than in review.

WHAT WOULD SETTLE IT
--------------------
A blank on the same batch as the moose libraries, with enough reads to have
seen a genus at the abundance the moose carry. I print the batch layout, then
the reads each moose genus has in the blanks of the moose batches specifically,
and say which of the two situations we are in.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_moose_batch_controls.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

FLOOR = 5.0


def background():
    """I return check_reagent_background, whose loading I reuse whole."""
    try:
        import check_reagent_background as module
    except ImportError as problem:
        sys.exit("cannot import check_reagent_background (%s). It must sit "
                 "in scripts/ beside me." % problem)
    return module


def main():
    module = background()
    shared = module.suppression()
    arms = shared.arms()
    lookup_table, lookup_module = shared.lookup()

    raw = module.open_table()
    features = [str(f) for f in raw.ids(axis="observation")]
    _, _, mapping = arms.read_taxonomy(features)
    keep = [f for f in features
            if not any(word in str(mapping.get(f, "")).lower()
                       for word in shared.HOST_WORDS)]
    raw = raw.filter(set(keep), axis="observation", inplace=False)
    genus_of = {}
    for feature in keep:
        genus_of[feature] = arms.genus_label(mapping.get(feature, ""))
    counts = module.collapse_to_genus(raw, genus_of)
    depth = counts.sum(axis=1)

    negatives = [i for i in counts.index
                 if module.NEGATIVE.search(i) and not module.POSITIVE.search(i)]
    held = set(negatives) | set(i for i in counts.index
                                if module.POSITIVE.search(i))
    samples = [i for i in counts.index if i not in held]
    animals = {i: shared.animal_of(lookup_table, lookup_module, i)
               for i in samples}
    deer = shared.deer_with_ixodes(lookup_table, lookup_module)
    moose = module.moose_animals(shared, lookup_table, lookup_module)
    moose_rows = [i for i in samples if animals[i] in moose]
    deer_rows = [i for i in samples if animals[i] in deer]
    if not moose_rows:
        sys.exit("no moose libraries found. Check moose_animals.")

    batch_of = {i: shared.run_of(i) for i in counts.index}
    batches = sorted(set(batch_of.values()))

    print("=" * 80)
    print("WHERE THE MOOSE SAT, AND WHICH PLATES HAD A BLANK")
    print("=" * 80)
    print("  %-26s %7s %7s %7s %12s" % ("batch", "moose", "deer", "blanks",
                                        "blank reads"))
    moose_batches = set(batch_of[i] for i in moose_rows)
    for batch in batches:
        here_moose = sum(1 for i in moose_rows if batch_of[i] == batch)
        here_deer = sum(1 for i in deer_rows if batch_of[i] == batch)
        here_blank = [i for i in negatives if batch_of[i] == batch]
        if not (here_moose or here_deer or here_blank):
            continue
        mark = " <-- moose" if batch in moose_batches else ""
        print("  %-26s %7d %7d %7d %12.0f%s"
              % (batch[:26], here_moose, here_deer, len(here_blank),
                 float(depth[here_blank].sum()) if here_blank else 0, mark))
    print()

    covered = sorted(b for b in moose_batches
                     if any(batch_of[i] == b for i in negatives))
    uncovered = sorted(moose_batches - set(covered))
    print("  moose sit on %d batch(es): %s"
          % (len(moose_batches), ", ".join(sorted(moose_batches))))
    print("  with a blank: %s" % (", ".join(covered) if covered else "none"))
    print("  without:      %s" % (", ".join(uncovered) if uncovered else "none"))
    print()

    same_plate = [i for i in negatives if batch_of[i] in moose_batches]
    if not same_plate:
        print("  NO BLANK SHARES A PLATE WITH THE MOOSE. The genera cleared by")
        print("  check_control_enrichment.py are cleared of the reagents we")
        print("  sampled and not of a contamination event on the moose plate.")
        print("  That is a limitation to state, not a result to publish.")
        return

    plate_reads = float(depth[same_plate].sum())
    print("=" * 80)
    print("THE MOOSE GENERA IN THE BLANKS THAT SHARED THEIR PLATES")
    print("=" * 80)
    print("  %d blank(s) on the moose plates, holding %.0f reads: %s"
          % (len(same_plate), plate_reads,
             ", ".join(str(s)[:28] for s in same_plate)))
    print()
    moose_depth = float(depth[moose_rows].sum())
    print("  %-26s %9s %10s %9s  %s"
          % ("genus", "in moose", "expected", "in blank", "reading"))
    for name, panels in sorted(module.wanted_genera().items()):
        if "moose" not in panels or name not in counts.columns:
            continue
        rate = float(counts.loc[moose_rows, name].sum() / moose_depth)
        expected = rate * plate_reads
        observed = float(counts.loc[same_plate, name].sum())
        if expected < FLOOR:
            reading = "plate blank too shallow to say"
        elif observed == 0:
            reading = "not on the plate, so not from it"
        elif observed >= expected:
            reading = "on the plate: suspect"
        else:
            reading = "present but far below the ticks"
        print("  %-26s %8.3f%% %10.0f %9.0f  %s"
              % (name[:26], 100 * rate, expected, observed, reading))
    print()
    print("  expected  reads the plate's own blank should hold if the plate")
    print("            were contaminated at the rate the moose ticks show")
    print()


if __name__ == "__main__":
    main()
