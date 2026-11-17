#!/usr/bin/env python3
"""check_moose_plate_confounding.py

I ask whether the moose genera are a property of the moose or of the plates the
moose sat on, by comparing them against the deer libraries on those same
plates. I write nothing.

WHY NOT THE BLANKS
------------------
check_moose_batch_controls.py found that 8 of the 13 moose libraries sit on
241209_P3-5, which carries no negative control, and that the blanks on the two
other moose plates hold 260 reads between them. The blanks cannot answer this.

The deer can. A reagent contamination event on a plate reaches every well on
that plate, so if Enhydrobacter and Deinococcus came from the P3-5 reagents
they must also be in the 123 deer libraries on P3-5. If they are in the moose
and not in the deer beside them, the plate is exonerated by a comparison with
123 libraries rather than none.

    high in moose, high in deer on the same plate     the plate
    high in moose, absent from deer on the same plate the moose

I also print whether each genus appears in the moose on all three plates. One
plate could be unlucky. Three agreeing is replication.

WHAT THIS CANNOT DO
-------------------
Every moose tick is Dermacentor and every deer tick in this comparison is
Ixodes, so host and tick genus are confounded by the study design and nothing
here separates them. I am testing the plate, not the host.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_moose_plate_confounding.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def background():
    """I return check_reagent_background, whose loading I reuse whole."""
    try:
        import check_reagent_background as module
    except ImportError as problem:
        sys.exit("cannot import check_reagent_background (%s). It must sit "
                 "in scripts/ beside me." % problem)
    return module


def profile(shares, rows, genus):
    """I return (libraries holding it, of how many, mean percent)."""
    if not rows:
        return 0, 0, 0.0
    part = shares.loc[rows, genus]
    return int((part > 0).sum()), len(rows), 100.0 * float(part.mean())


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

    held = set(i for i in counts.index
               if module.NEGATIVE.search(i) or module.POSITIVE.search(i))
    samples = [i for i in counts.index if i not in held]
    animals = {i: shared.animal_of(lookup_table, lookup_module, i)
               for i in samples}
    deer = shared.deer_with_ixodes(lookup_table, lookup_module)
    moose = module.moose_animals(shared, lookup_table, lookup_module)
    moose_rows = [i for i in samples if animals[i] in moose]
    deer_rows = [i for i in samples if animals[i] in deer]
    if not moose_rows:
        sys.exit("no moose libraries found. Check moose_animals.")

    batch_of = {i: shared.run_of(i) for i in samples}
    moose_batches = sorted(set(batch_of[i] for i in moose_rows),
                           key=lambda b: -sum(1 for i in moose_rows
                                              if batch_of[i] == b))
    genera = [n for n, panels in sorted(module.wanted_genera().items())
              if "moose" in panels and n in counts.columns]
    if not genera:
        sys.exit("no moose-panel genera survive in the table.")

    for batch in moose_batches:
        here_moose = [i for i in moose_rows if batch_of[i] == batch]
        here_deer = [i for i in deer_rows if batch_of[i] == batch]
        print("=" * 84)
        print("%s: %d moose libraries against %d deer libraries on the same "
              "plate" % (batch, len(here_moose), len(here_deer)))
        print("=" * 84)
        if not here_deer:
            print("  No deer on this plate, so it has no internal comparison.")
            print()
            continue
        print("  %-26s %16s %16s  %s"
              % ("genus", "moose", "deer here", "reading"))
        print("  %-26s %16s %16s" % ("", "in / n   mean %", "in / n   mean %"))
        for name in genera:
            m_in, m_n, m_mean = profile(shares, here_moose, name)
            d_in, d_n, d_mean = profile(shares, here_deer, name)
            m_rate = m_in / float(m_n) if m_n else 0.0
            d_rate = d_in / float(d_n) if d_n else 0.0
            if m_rate < 0.5:
                reading = "not general in these moose"
            elif d_rate >= 0.5 and d_mean >= 0.1 * m_mean:
                reading = "the plate"
            elif d_rate <= 0.25:
                reading = "the moose, not the plate"
            else:
                reading = "mixed"
            print("  %-26s %5d /%-4d %7.3f %5d /%-4d %7.3f  %s"
                  % (name[:26], m_in, m_n, m_mean, d_in, d_n, d_mean, reading))
        print()

    print("=" * 84)
    print("DOES EACH GENUS APPEAR IN THE MOOSE ON EVERY PLATE?")
    print("=" * 84)
    print("  %-26s %s" % ("genus", "  ".join(b[:14] for b in moose_batches)))
    for name in genera:
        cells = []
        for batch in moose_batches:
            here = [i for i in moose_rows if batch_of[i] == batch]
            hits, total, _ = profile(shares, here, name)
            cells.append("%d/%d" % (hits, total))
        every = all(int(c.split("/")[0]) > 0 for c in cells)
        print("  %-26s %s   %s"
              % (name[:26], "  ".join("%-14s" % c for c in cells),
                 "on every plate" if every else ""))
    print()
    print("  A genus in the moose on all three plates cannot be one plate's")
    print("  contamination. A genus on P3-5 only needs the deer column above")
    print("  to clear it, and if the deer there carry it too, it is the plate.")
    print()


if __name__ == "__main__":
    main()
