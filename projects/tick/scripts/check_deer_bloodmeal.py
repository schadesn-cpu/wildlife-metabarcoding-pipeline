#!/usr/bin/env python3
"""check_deer_bloodmeal.py

I run the second-blood-meal threshold on the Ixodes specimens from
white-tailed deer, so the claim in Section 3.1 covers the whole dataset rather
than the 20 Dermacentor. I write nothing.

WHY THIS ONE IS DIFFERENT
-------------------------
Ixodes scapularis is a three-host tick. Its larvae and nymphs feed on small
mammals and birds before the adult feeds on a deer, so a prior host is
expected here in a way it is not in the one-host winter tick. If bird or rodent
hosts clear the threshold in deer ticks and nothing clears it in the winter
ticks, the contrast is in the data rather than only in the literature.

I use the same rule as check_second_bloodmeal.py: a secondary host group must
reach 100 reads and 1% of that specimen's host reads, and I say whether the
same sequencing run carried specimens collected from that host, which is what
separates a prior meal from cross-contamination.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_deer_bloodmeal.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SPECIES = "results/analysis/tick_species_by_host.tsv"
THRESHOLD = 0.01
MIN_READS = 100
# A specimen needs this many host reads before any blood meal is called.
# Below it the dominant group is unstable: the four specimens whose dominant
# group disagreed with the collection record all sat between 74 and 2,017
# host reads, against a median of 80,187.
MIN_HOST_READS = 2500
UNNAMED = "Unassigned or other"
SHOW = 20


def helpers():
    try:
        import check_second_bloodmeal as base
    except ImportError as problem:
        sys.exit("cannot import check_second_bloodmeal (%s). It must be in "
                 "scripts/ beside me." % problem)
    return base


def deer_specimens(module, table, lookup_module):
    """I return the deer animals carrying Ixodes and not Dermacentor."""
    if not os.path.isfile(SPECIES):
        sys.exit("missing %s" % SPECIES)
    species = pd.read_csv(SPECIES, sep="\t", dtype=str)
    positive = ("1", "1.0", "true", "TRUE", "True")
    found = set()
    for specimen, host, ixodes, derm in zip(species["specimen"],
                                            species["host_bin"],
                                            species["Ixodes_pos"],
                                            species["Dermacentor_pos"]):
        lowered = str(host).lower()
        if "deer" not in lowered and "wtd" not in lowered:
            continue
        if str(ixodes).strip() not in positive:
            continue
        if str(derm).strip() in positive:
            continue
        animal = module.animal_of(table, lookup_module, specimen)
        if animal:
            found.add(animal)
    return found


def main():
    base = helpers()
    module = base.bloodmeal()
    table, lookup_module = module.lookup()
    counts = module.read_table()
    path, covered, mapping = module.read_taxonomy(list(counts.columns))
    print("cytb table: %d libraries x %d ASVs" % counts.shape)
    print("taxonomy:   %s (%d of %d ASVs)" % (path, covered, counts.shape[1]))

    groups = pd.Series({f: module.group_of(mapping.get(f, ""))
                        for f in counts.columns})
    by_group = counts.T.groupby(groups).sum().T
    named = [c for c in by_group.columns if c != UNNAMED]
    animals = pd.Series([module.animal_of(table, lookup_module, s)
                         for s in counts.index], index=counts.index)
    runs = pd.Series([base.run_of(s) for s in counts.index], index=counts.index)
    hosts = base.read_hosts()

    on_run = {}
    for library, animal in animals.items():
        group = base.collected_group(hosts.get(animal, ""))
        if group:
            on_run.setdefault(runs[library], set()).add(group)

    deer = deer_specimens(module, table, lookup_module)
    print("Ixodes specimens from deer: %d" % len(deer))

    rows, with_reads, too_thin, called = [], 0, 0, 0
    for animal in sorted(deer):
        libraries = animals.index[animals == animal]
        if not len(libraries):
            continue
        totals = by_group.loc[libraries, named].sum(axis=0)
        total = float(totals.sum())
        if total <= 0:
            continue
        with_reads += 1
        if total < MIN_HOST_READS:
            too_thin += 1
            continue
        called += 1
        ordered = totals.sort_values(ascending=False)
        explained_here = set()
        for one in set(runs[libraries]):
            explained_here |= on_run.get(one, set())
        for group in ordered.index[1:]:
            reads = float(ordered[group])
            share = reads / total
            if reads < MIN_READS or share < THRESHOLD:
                continue
            rows.append({"specimen": animal, "dominant": ordered.index[0],
                         "second": group, "reads": int(reads), "share": share,
                         "on_run": group in explained_here})
    print("with cytb host reads: %d" % with_reads)
    print("below the %d host-read floor, not called: %d"
          % (MIN_HOST_READS, too_thin))
    print("blood meal called: %d" % called)
    print()

    print("=" * 78)
    print("SECONDARY HOSTS ABOVE %d READS AND %.0f%% OF HOST READS"
          % (MIN_READS, 100 * THRESHOLD))
    print("=" * 78)
    if not rows:
        print("  none. Every deer tick has a single host above the threshold.")
        print()
        return
    frame = pd.DataFrame(rows)
    specimens = len(set(frame["specimen"]))
    print("  %d of %d specimens carry at least one, %d rows in total"
          % (specimens, called, len(frame)))
    print()
    print("  %-22s %-24s %10s %8s %8s"
          % ("specimen", "second host", "reads", "share", "on run"))
    for _, row in frame.sort_values("share", ascending=False).head(SHOW).iterrows():
        print("  %-22s %-24s %10d %7.1f%% %8s"
              % (str(row["specimen"])[:22], str(row["second"])[:24],
                 row["reads"], 100 * row["share"],
                 "yes" if row["on_run"] else "NO"))
    if len(frame) > SHOW:
        print("  ... and %d more" % (len(frame) - SHOW))
    print()

    print("  by host group:")
    print("  %-24s %8s %16s %14s"
          % ("second host", "rows", "explained by run", "unexplained"))
    for group, part in frame.groupby("second"):
        print("  %-24s %8d %16d %14d"
              % (str(group)[:24], len(part), int(part["on_run"].sum()),
                 int((~part["on_run"]).sum())))
    print()
    print("  Human cannot be explained by run, because no specimen was")
    print("  collected from a human. Read those rows as handling, not as a")
    print("  meal. Bird and small mammal rows that no run explains are the")
    print("  ones consistent with a prior host in a three-host tick.")
    print()


if __name__ == "__main__":
    main()
