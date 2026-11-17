#!/usr/bin/env python3
"""check_specimen_bloodmeal.py

I take one specimen whose blood meal does not match its collection host and
decide between the three explanations. I write nothing.

WHY
---
USGSTICKS-4109 was collected from a deer and 29,473 of its 30,482 cytb host
reads are porcupine. That is 96.7%, and it is not a shallow library, so the
read floor that explains the other four mismatches cannot reach it. Three
things could produce it:

    a specimen mix-up      the tick came off a porcupine
    index hopping          porcupine reads leaked from another library
    a prior host           the tick fed on a porcupine before this meal

They have very different consequences. A mix-up is a data-integrity item. A
prior host in a three-host tick would be the one observation in this dataset
that bears on the multi-host life cycle, which Section 3.1 currently says
nothing does.

HOW I SEPARATE THEM
-------------------
Hopping sprays. If porcupine reads reached this library from elsewhere on its
run, other libraries on that run should carry them too, and this library's
porcupine share should sit near the run's own share. So I compute what this
library would hold if it were a copy of its run, with itself excluded from
that calculation, and set it beside what it actually holds. A ratio near 1 is
hopping. A ratio far above 1 means the reads are specific to this specimen,
which leaves a mix-up or a real prior meal, and those two are told apart by
the collection record rather than by sequence.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_specimen_bloodmeal.py
    python3 scripts/check_specimen_bloodmeal.py SOME-OTHER-SPECIMEN
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SPECIMEN = "USGSTICKS-4109"
UNNAMED = "Unassigned or other"
SHOW = 15


def deer_module():
    try:
        import check_deer_bloodmeal as module
    except ImportError as problem:
        sys.exit("cannot import check_deer_bloodmeal (%s). It must sit in "
                 "scripts/ beside me." % problem)
    return module


def main():
    wanted = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") \
        else SPECIMEN
    deer_side = deer_module()
    base = deer_side.helpers()
    module = base.bloodmeal()
    table, finder = module.lookup()
    counts = module.read_table()
    path, covered, mapping = module.read_taxonomy(list(counts.columns))
    print("cytb table: %d libraries x %d ASVs" % counts.shape)

    groups = pd.Series({f: module.group_of(mapping.get(f, ""))
                        for f in counts.columns})
    by_group = counts.T.groupby(groups).sum().T
    named = [c for c in by_group.columns if c != UNNAMED]
    animals = pd.Series([module.animal_of(table, finder, s)
                         for s in counts.index], index=counts.index)
    runs = pd.Series([base.run_of(s) for s in counts.index], index=counts.index)
    hosts = base.read_hosts()

    mine = list(animals.index[animals == wanted])
    if not mine:
        sys.exit("no cytb library resolves to %r. Check the spelling, or pass "
                 "a different specimen." % wanted)

    print("specimen: %s, collected from %r"
          % (wanted, hosts.get(wanted, "(not recorded)")))
    print("libraries: %s" % ", ".join(mine))
    print("runs:      %s" % ", ".join(sorted(set(runs[mine]))))
    print()

    totals = by_group.loc[mine, named].sum(axis=0).sort_values(ascending=False)
    total = float(totals.sum())
    print("=" * 70)
    print("ITS HOST READS")
    print("=" * 70)
    for group in totals.index:
        if totals[group] > 0:
            print("  %-28s %10s %7.1f%%"
                  % (group[:28], "{:,}".format(int(totals[group])),
                     100.0 * totals[group] / total))
    print("  %-28s %10s" % ("total", "{:,}".format(int(total))))
    print()
    subject = totals.index[0]

    for run in sorted(set(runs[mine])):
        fellows = [l for l in runs.index[runs == run] if l not in mine]
        print("=" * 70)
        print("RUN %s: %d OTHER LIBRARIES" % (run, len(fellows)))
        print("=" * 70)
        if not fellows:
            print("  none. This specimen is alone on its run, so nothing")
            print("  could have hopped into it from a neighbour.")
            print()
            continue

        collected = {}
        for library in fellows:
            group = base.collected_group(hosts.get(animals[library], ""))
            collected[group or "(not recorded)"] = \
                collected.get(group or "(not recorded)", 0) + 1
        print("  collection hosts of the neighbours: %s"
              % ", ".join("%s %d" % (k, v) for k, v in
                          sorted(collected.items(), key=lambda kv: -kv[1])))

        elsewhere = by_group.loc[fellows, subject]
        carrying = elsewhere[elsewhere > 0].sort_values(ascending=False)
        print("  neighbours carrying %s reads: %d of %d"
              % (subject, len(carrying), len(fellows)))
        for library in list(carrying.index)[:SHOW]:
            print("    %-34s %10s  (%s)"
                  % (str(library)[:34], "{:,}".format(int(carrying[library])),
                     base.collected_group(hosts.get(animals[library], ""))
                     or "host not recorded"))
        print()

        # Every group, not just the dominant one. A group sitting at run
        # rate arrived by hopping; a group far above it belongs to this
        # library. Doing this for all of them says which parts of the
        # specimen's profile are its own.
        pool_named = by_group.loc[fellows, named]
        pool_total = float(pool_named.values.sum())
        print("  if this library were a copy of its run:")
        print("    %-26s %10s %10s %10s"
              % ("host group", "observed", "run rate", "ratio"))
        for group in totals.index:
            observed = float(totals[group])
            if observed <= 0:
                continue
            share = (float(pool_named[group].sum()) / pool_total
                     if pool_total else 0.0)
            expected = share * total
            ratio = ("%.0f" % (observed / expected) if expected >= 1
                     else "no source on run")
            print("    %-26s %10s %10s %10s"
                  % (group[:26], "{:,}".format(int(observed)),
                     "{:,}".format(int(expected)), ratio))
        print()

    print("=" * 70)
    print("READING IT")
    print("=" * 70)
    print("  A ratio near 1 means hopping and the specimen is unremarkable.")
    print("  A ratio in the hundreds or thousands means the %s reads belong"
          % subject)
    print("  to this library, which leaves a specimen mix-up or a genuine")
    print("  prior host. Sequence cannot tell those apart; the collection")
    print("  record can.")
    print()


if __name__ == "__main__":
    main()
