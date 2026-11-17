#!/usr/bin/env python3
"""check_second_bloodmeal.py

I ask whether any Dermacentor specimen carries a real second blood meal, or
whether the second hosts in the blood-meal matrix are contamination. I write
nothing.

WHY
---
The existing matrix counts a host group as a meal on one read, so a specimen
with 74,235 cervid reads and 11 bird reads appears twice and reads as two
meals. Two things make that worth resolving rather than caveating.

    Dermacentor albipictus, the winter tick on moose, is a one-host tick: it
    feeds on one animal for its whole life cycle, so a second host is not
    expected.
    Dermacentor variabilis, on opossum and porcupine, is a three-host tick, so
    a host from an earlier life stage is entirely plausible.

If genuine second hosts appear on the three-host ticks and not on the one-host
ticks, that is independent support for assigning species by collection host.
If they appear on neither above a sensible threshold, there are no second
meals and the figure was reading contamination.

HOW I SEPARATE THE TWO
----------------------
For every secondary host group I ask whether the same sequencing run carried
specimens actually collected from that host. A bear signal in a tick from a
moose, on a run full of bear ticks, is cross-contamination. The same signal on
a run with no bear is not explained that way and deserves a look.

I report at two thresholds so the answer does not rest on one cutoff.

WHAT I READ
-----------
    qiime2/cytb/all/dada2/table.qza and the verified cytb taxonomy, through
    check_winter_tick_bloodmeal
    results/metadata/crosswalk_4marker.tsv        collection host
    results/analysis/tick_species_by_host.tsv     which specimens are which

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_second_bloodmeal.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
SPECIES = "results/analysis/tick_species_by_host.tsv"
THRESHOLDS = (0.01, 0.05)
MIN_READS = 100
CERVID = "Cervidae (deer or moose)"
UNNAMED = "Unassigned or other"
# I map a collection host to the blood-meal group it should produce.
COLLECTED_AS = [
    (("moose", "alces"), CERVID),
    (("deer", "wtd", "odocoileus"), CERVID),
    (("bear", "ursus"), "Bear"),
    (("opossum", "didelph"), "Opossum"),
    (("porcupine", "erethizon"), "Porcupine"),
    (("raccoon", "procyon"), "Raccoon"),
    (("rabbit", "hare", "cottontail", "sylvilagus"), "Rabbit"),
]


def bloodmeal():
    try:
        import check_winter_tick_bloodmeal as module
    except ImportError as problem:
        sys.exit("cannot import check_winter_tick_bloodmeal (%s). Run me from "
                 "the project root." % problem)
    return module


def run_of(library):
    text = str(library)
    return text.rsplit(".", 1)[1] if "." in text else "unknown"


def collected_group(text):
    """I return the blood-meal group a collection host should produce."""
    lowered = str(text).lower()
    for words, group in COLLECTED_AS:
        if any(word in lowered for word in words):
            return group
    return ""


def read_hosts():
    """I return {specimen: collection host text} from the crosswalk."""
    if not os.path.isfile(CROSSWALK):
        sys.exit("missing %s. Run me from the project root." % CROSSWALK)
    frame = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    for column in ("sample_id", "host_spp", "host_bin"):
        if column not in frame.columns:
            sys.exit("column '%s' missing from %s" % (column, CROSSWALK))
    hosts = {}
    for specimen, spp, binned in zip(frame["sample_id"], frame["host_spp"],
                                     frame["host_bin"]):
        text = str(spp).strip()
        if text.lower() in ("", "nan", "none"):
            text = str(binned).strip()
        hosts[specimen] = text
    return hosts


def dermacentor_specimens(module, table, lookup_module):
    """I return {animal: collection host} for the Dermacentor specimens."""
    if not os.path.isfile(SPECIES):
        sys.exit("missing %s" % SPECIES)
    species = pd.read_csv(SPECIES, sep="\t", dtype=str)
    positive = ("1", "1.0", "true", "TRUE", "True")
    found = {}
    for specimen, derm in zip(species["specimen"], species["Dermacentor_pos"]):
        if str(derm).strip() not in positive:
            continue
        animal = module.animal_of(table, lookup_module, specimen)
        if animal:
            found[animal] = None
    return found


def main():
    module = bloodmeal()
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
    runs = pd.Series([run_of(s) for s in counts.index], index=counts.index)
    hosts = read_hosts()

    # What each run actually carried, so contamination is separable.
    on_run = {}
    for library, animal in animals.items():
        group = collected_group(hosts.get(animal, ""))
        if group:
            on_run.setdefault(runs[library], set()).add(group)

    derm = dermacentor_specimens(module, table, lookup_module)
    for animal in derm:
        derm[animal] = hosts.get(animal, "")
    print("Dermacentor specimens: %d" % len(derm))
    print()

    rows = []
    for animal, host in sorted(derm.items(), key=lambda kv: str(kv[1])):
        libraries = animals.index[animals == animal]
        if not len(libraries):
            continue
        totals = by_group.loc[libraries, named].sum(axis=0)
        total = float(totals.sum())
        if total <= 0:
            continue
        ordered = totals.sort_values(ascending=False)
        dominant = ordered.index[0]
        expected = collected_group(host)
        species = "albipictus" if "moose" in str(host).lower() else "variabilis"
        seen_runs = set(runs[libraries])
        explained_here = set()
        for one in seen_runs:
            explained_here |= on_run.get(one, set())
        for group in ordered.index[1:]:
            reads = float(ordered[group])
            share = reads / total
            if reads < MIN_READS or share < min(THRESHOLDS):
                continue
            rows.append({
                "specimen": animal, "host": host, "species": species,
                "dominant": dominant, "dominant_share": ordered[dominant] / total,
                "second": group, "second_reads": int(reads),
                "second_share": share,
                "on_run": group in explained_here,
                "matches_host": dominant == expected,
            })

    print("=" * 78)
    print("SECONDARY HOST GROUPS IN DERMACENTOR SPECIMENS")
    print("=" * 78)
    if not rows:
        print("  none reaches %d reads and %.0f%% of host reads."
              % (MIN_READS, 100 * min(THRESHOLDS)))
        print("  Every Dermacentor specimen has one blood meal at this "
              "threshold.")
        print()
        return
    frame = pd.DataFrame(rows)
    print("  %-20s %-12s %-22s %9s %7s %8s"
          % ("specimen", "species", "second host", "reads", "share", "on run"))
    for _, row in frame.sort_values("second_share", ascending=False).iterrows():
        print("  %-20s %-12s %-22s %9d %6.1f%% %8s"
              % (str(row["specimen"])[:20], row["species"],
                 str(row["second"])[:22], row["second_reads"],
                 100 * row["second_share"], "yes" if row["on_run"] else "NO"))
    print()

    print("=" * 78)
    print("HOW MANY SPECIMENS HAVE A SECOND HOST, BY THRESHOLD")
    print("=" * 78)
    print("  %-12s %-10s %8s %14s %16s"
          % ("threshold", "species", "any", "explained by run", "unexplained"))
    for cut in THRESHOLDS:
        for species in ("albipictus", "variabilis"):
            part = frame[(frame["species"] == species)
                         & (frame["second_share"] >= cut)]
            total = len(set(part["specimen"]))
            explained = len(set(part.loc[part["on_run"], "specimen"]))
            unexplained = len(set(part.loc[~part["on_run"], "specimen"]))
            print("  %-12s %-10s %8d %14d %16d"
                  % ("%.0f%%" % (100 * cut), species, total, explained,
                     unexplained))
    print()
    print("  albipictus is a one-host tick, so an unexplained second host "
          "there")
    print("  is the surprising result. variabilis is a three-host tick, where "
          "a")
    print("  second host is expected rather than remarkable.")
    print()


if __name__ == "__main__":
    main()
