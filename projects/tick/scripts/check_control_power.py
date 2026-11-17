#!/usr/bin/env python3
"""check_control_power.py

I say whether the negative controls were deep enough for a genus's absence
from them to mean anything. I write nothing.

WHY
---
check_reagent_background.py reports relative abundance in the controls, and
that was the wrong currency. The blanks have a median depth of 246 reads, so a
handful of hopped reads reads as a large percentage and every genus looks
contaminated. Absence is the part of that output that does not depend on the
denominator, but absence is only evidence when detection was possible.

So I ask the question the other way round. If the blanks carried the same
community as the ticks, how many reads of each genus should they hold? I sum
that expectation across the blanks, using each blank's own depth, and set it
beside what is actually there.

    expected many, observed none     absent, and the blanks could have seen it
    expected few, observed none      underpowered, no conclusion
    observed at or above expected    background, or cross-talk

I print the roster of blanks with their depths too, because one of them holds
58,695 reads and a blank that deep is either doing all the work in this test
or is not a blank.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_control_power.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DETECTABLE = 20.0
UNDERPOWERED = 5.0


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
    shares = counts.div(depth.replace(0, float("nan")), axis=0).fillna(0)

    negatives = [i for i in counts.index
                 if module.NEGATIVE.search(i) and not module.POSITIVE.search(i)]
    held = set(negatives) | set(i for i in counts.index
                                if module.POSITIVE.search(i))
    others = [i for i in counts.index if i not in held]
    animals = {i: shared.animal_of(lookup_table, lookup_module, i)
               for i in others}
    deer = shared.deer_with_ixodes(lookup_table, lookup_module)
    moose = module.moose_animals(shared, lookup_table, lookup_module)
    cohorts = {"deer": [i for i in others if animals[i] in deer],
               "moose": [i for i in others if animals[i] in moose]}
    if not negatives:
        sys.exit("no negative controls found. Check the naming.")

    print("=" * 78)
    print("THE NEGATIVE CONTROLS, DEEPEST FIRST")
    print("=" * 78)
    print("  %-46s %10s %12s" % ("library", "reads", "top genus"))
    for name in depth[negatives].sort_values(ascending=False).index:
        row = counts.loc[name]
        top = row.idxmax() if row.sum() > 0 else "none"
        print("  %-46s %10.0f %12s" % (str(name)[:46], depth[name], str(top)[:12]))
    total_depth = float(depth[negatives].sum())
    print()
    print("  %d blanks, %.0f reads in total" % (len(negatives), total_depth))
    print()

    print("=" * 78)
    print("COULD THE BLANKS HAVE SEEN IT?")
    print("=" * 78)
    print("  %-24s %-6s %9s %9s %9s  %s"
          % ("genus", "panel", "in ticks", "expected", "observed", "verdict"))
    for name, panels in sorted(module.wanted_genera().items()):
        panel = panels[-1]
        rows = cohorts[panel]
        if not rows or name not in counts.columns:
            continue
        rate = float(shares.loc[rows, name].mean())
        expected = rate * total_depth
        observed = float(counts.loc[negatives, name].sum())
        if observed >= expected:
            verdict = "background or cross-talk"
        elif observed == 0 and expected >= DETECTABLE:
            verdict = "absent, and detectable"
        elif expected < UNDERPOWERED:
            verdict = "underpowered, no call"
        else:
            verdict = "far below expectation"
        print("  %-24s %-6s %8.3f%% %9.0f %9.0f  %s"
              % (name[:24], panel, 100 * rate, expected, observed, verdict))
    print()
    print("  expected  reads the blanks should hold if they carried the same")
    print("            community as that panel's ticks, at each blank's depth")
    print("  detectable threshold %.0f reads; below %.0f I make no call"
          % (DETECTABLE, UNDERPOWERED))
    print()


if __name__ == "__main__":
    main()
