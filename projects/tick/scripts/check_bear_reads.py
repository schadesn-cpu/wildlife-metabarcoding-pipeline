#!/usr/bin/env python3
"""check_bear_reads.py

I ask why 42 cytb libraries from white-tailed deer carry black bear reads. I
write nothing.

THE QUESTION
------------
check_winter_tick_bloodmeal.py reports 2,051,050 bear cytb reads across 42
libraries whose specimens are deer. That is far above the trace level of the
human reads, so it is either a host mislabel or cross-contamination between
samples on a plate.

The two explanations look different:

  index hopping or well bleed   bear reads are a small share of the library and
                                the library sits on a sequencing run that also
                                carries a genuine bear sample with many bear
                                reads.
  host mislabel                 bear reads dominate the library, and the run
                                need not carry a bear sample at all.

So I print, for every deer library with bear reads, the bear share of its host
reads and its run, and then per run I put the bear-read deer libraries next to
the genuine bear libraries on the same run.

I reuse check_winter_tick_bloodmeal for the table, the taxonomy and the host
groups, so the calls here are the same ones the blood meal result rests on.

WHAT I READ
-----------
    qiime2/cytb/all/dada2/table.qza
    the verified cytb taxonomy
    results/metadata/crosswalk_4marker.tsv   host per specimen
    scripts/tick_lookup.py

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_bear_reads.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
BEAR = "Bear"
CERVID = "Cervidae (deer or moose)"
SHOW = 45


def bloodmeal():
    try:
        import check_winter_tick_bloodmeal as module
    except ImportError as problem:
        sys.exit("cannot import check_winter_tick_bloodmeal (%s). Run me from "
                 "the project root." % problem)
    return module


def run_of(library):
    """I take the run label off a library name, or return 'unknown'.

    Library names look like ALHT-Deer-09.241209_P3-5, so the run is whatever
    follows the last dot.
    """
    text = str(library)
    return text.rsplit(".", 1)[1] if "." in text else "unknown"


def hosts_by_specimen():
    """I return {specimen: host text} from the crosswalk."""
    if not os.path.exists(CROSSWALK):
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


def is_bear(text):
    return "bear" in str(text).lower() or "ursus" in str(text).lower()


def is_deer(text):
    lowered = str(text).lower()
    return "deer" in lowered or "wtd" in lowered or "odocoileus" in lowered


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
    if BEAR not in by_group.columns:
        sys.exit("no ASV is assigned to the '%s' host group. Nothing to "
                 "explain." % BEAR)

    hosts = hosts_by_specimen()
    specimens = [module.animal_of(table, lookup_module, s) for s in counts.index]
    work = pd.DataFrame(index=counts.index)
    work["specimen"] = specimens
    work["host"] = [hosts.get(s, "") for s in specimens]
    work["run"] = [run_of(s) for s in counts.index]
    work["bear"] = by_group[BEAR].reindex(counts.index).fillna(0).values
    work["cervid"] = (by_group[CERVID].reindex(counts.index).fillna(0).values
                      if CERVID in by_group.columns else 0)
    named = [c for c in by_group.columns if c != "Unassigned or other"]
    work["host_reads"] = by_group[named].reindex(counts.index).fillna(0).sum(axis=1).values
    work["bear_share"] = work["bear"] / work["host_reads"].replace(0, float("nan"))

    bear_libraries = work[work["host"].map(is_bear)]
    deer_with_bear = work[work["host"].map(is_deer) & (work["bear"] > 0)]
    print("libraries resolved to a specimen: %d of %d"
          % (int((work["specimen"] != "").sum()), len(work)))
    print("libraries from bear:              %d (%d bear reads)"
          % (len(bear_libraries), int(bear_libraries["bear"].sum())))
    print("libraries from deer with any bear read: %d (%d bear reads)"
          % (len(deer_with_bear), int(deer_with_bear["bear"].sum())))
    print()

    print("=" * 78)
    print("DEER LIBRARIES CARRYING BEAR READS")
    print("=" * 78)
    ordered = deer_with_bear.sort_values("bear", ascending=False)
    print("  %-34s %-18s %9s %9s %7s" % ("library", "run", "bear", "cervid",
                                         "bear %"))
    for name, row in ordered.head(SHOW).iterrows():
        share = row["bear_share"]
        print("  %-34s %-18s %9d %9d %6s"
              % (str(name)[:34], str(row["run"])[:18], int(row["bear"]),
                 int(row["cervid"]),
                 "%.1f%%" % (100 * share) if pd.notna(share) else "?"))
    if len(ordered) > SHOW:
        print("  ... and %d more" % (len(ordered) - SHOW))
    print()
    dominant = int((ordered["bear_share"] > 0.5).sum())
    print("  bear reads are more than half the host reads in %d of %d libraries"
          % (dominant, len(ordered)))
    print("  median bear share %.2f%%, maximum %.1f%%"
          % (100 * ordered["bear_share"].median(),
             100 * ordered["bear_share"].max()))
    print()

    print("=" * 78)
    print("BY SEQUENCING RUN: ARE THE BEAR READS NEXT TO A REAL BEAR?")
    print("=" * 78)
    print("  %-18s %6s %8s %11s %9s %9s"
          % ("run", "libs", "bear libs", "bear reads", "deer+bear", "median %"))
    rows = []
    for run, part in work.groupby("run"):
        on_run_bear = part[part["host"].map(is_bear)]
        on_run_deer = part[part["host"].map(is_deer) & (part["bear"] > 0)]
        if on_run_bear.empty and on_run_deer.empty:
            continue
        median = on_run_deer["bear_share"].median()
        print("  %-18s %6d %8d %11d %9d %9s"
              % (str(run)[:18], len(part), len(on_run_bear),
                 int(on_run_bear["bear"].sum()), len(on_run_deer),
                 "%.2f%%" % (100 * median) if pd.notna(median) else "-"))
        rows.append((run, len(on_run_bear), len(on_run_deer)))
    print()

    with_bear_sample = sum(n for _, bears, n in rows if bears > 0)
    without = sum(n for _, bears, n in rows if bears == 0)
    print("  deer libraries with bear reads on a run that also holds a bear "
          "sample: %d" % with_bear_sample)
    print("  on a run that holds no bear sample:                              "
          "   %d" % without)
    print()
    print("  If most sit on runs with a real bear sample and carry a small")
    print("  share, this is cross-contamination and the hosts are fine. If")
    print("  they sit on runs with no bear, or bear dominates the library,")
    print("  the collection records need checking.")
    print()


if __name__ == "__main__":
    main()
