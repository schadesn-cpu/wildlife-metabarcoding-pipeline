#!/usr/bin/env python3
"""I test whether the Bartonella signal tracks deer ked co-detection.

WHY
---
The taxon decision table classifies Bartonella as ked_associated_bacterium:
ruminant-associated and hippoboscid-linked, with 34 of 45 ASVs tying
Bartonella schoenbuchensis and two more hitting B. melophagi, reclassified
on 2026-08-04 as not a tick-borne human pathogen and to be reported as a
co-detection outside the pathogen results. The same table records the deer
ked, Lipoptena cervi, as an ectoparasite co-detection noted in the
accessioning records for 29 specimens and detected on three markers.

That attribution comes from the literature. It has not been tested against
this collection, which can test it: if these Bartonella reads came from keds
in the tube, Bartonella should be found in the specimens where a ked was
also detected.

WHO CAN BE SCORED
-----------------
The specimen taxon table lists detections only, so a specimen with no
Lipoptena row either carried no ked or was never sequenced on a marker that
could see one. Those are different and I do not merge them. A specimen is
scored only if the crosswalk gives it a 16S library, without which
Bartonella cannot be found, and an 18S library, without which the ked
cannot. Specimens missing either are counted as unscoreable, never as
negative.

NAMES
-----
Every identifier goes through tick_lookup. The specimen taxon table and the
crosswalk do not always name the same animal the same way, so a lookup miss
is reported as a miss and never repaired by string surgery.

WHAT I DO NOT CLAIM
-------------------
Co-detection in one tube is not a transmission route, and this cannot
separate Bartonella carried by a ked from Bartonella in the deer's blood
that the ked and the tick both acquired. The Methods also record that keds
were removed from tubes partway through collection, so ked detection depends
on when a sample was handled. I report the association; I do not adjust for
that.

OUTPUT
------
A summary to stdout and one table at
results/analysis/ked_bartonella_association.tsv. I create no directories.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import tick_lookup

SPECIMEN_TAXA = "results/tables/specimen_taxon_table.tsv"
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
OUT_PATH = "results/analysis/ked_bartonella_association.tsv"

KED_PATTERN = "lipoptena|hippobosc"
BARTONELLA = "bartonella"
# A ked is visible on these markers; Bartonella only on 16S.
KED_MARKERS = ("18S", "cytb", "ITS")


def blank(value):
    """I report whether a cell is empty."""
    if value is None or pd.isna(value):
        return True
    return str(value).strip() == ""


def fisher(a, b, c, d):
    """I return the two-sided Fisher exact p, or None if scipy is absent."""
    try:
        from scipy.stats import fisher_exact
    except ImportError:
        return None
    return float(fisher_exact([[a, b], [c, d]])[1])


def main():
    if not os.path.isdir("results/analysis"):
        sys.exit("run me from the project root; there is no results/analysis here.")
    for path in (SPECIMEN_TAXA, CROSSWALK):
        if not os.path.isfile(path):
            sys.exit("missing input: " + path)

    table = tick_lookup.load()

    taxa = pd.read_csv(SPECIMEN_TAXA, sep="\t", dtype=str)
    for column in ("specimen", "marker", "taxon", "reads"):
        if column not in taxa.columns:
            sys.exit("%s has no %s column; it has %s"
                     % (SPECIMEN_TAXA, column, ", ".join(taxa.columns)))
    taxa["reads"] = pd.to_numeric(taxa["reads"], errors="coerce").fillna(0)
    taxa["animal"] = [tick_lookup.resolve(table, "specimen", s)
                      for s in taxa["specimen"]]
    missed = sorted(set(taxa.loc[taxa["animal"].isna(), "specimen"].astype(str)))
    print("specimen-taxon rows %d, of which identifiers that did not resolve %d"
          % (len(taxa), int(taxa["animal"].isna().sum())))
    if missed:
        print("   unresolved: %s" % ", ".join(missed[:10]))
        print("   I leave these out rather than guess which animal they are.")

    ked_rows = taxa[taxa["taxon"].str.lower().str.contains(KED_PATTERN, na=False)]
    bart_rows = taxa[(taxa["marker"] == "16S")
                     & (taxa["taxon"].str.lower() == BARTONELLA)]
    print()
    print("ked rows %d on markers %s"
          % (len(ked_rows), ", ".join(sorted(set(ked_rows["marker"].astype(str))))))
    print("Bartonella rows %d on 16S" % len(bart_rows))

    ked_reads, bart_reads, bart_asvs = {}, {}, {}
    for animal, group in ked_rows.dropna(subset=["animal"]).groupby("animal"):
        ked_reads[animal] = float(group["reads"].max())
    for animal, group in bart_rows.dropna(subset=["animal"]).groupby("animal"):
        bart_reads[animal] = float(group["reads"].sum())
        if "n_asvs" in group.columns:
            bart_asvs[animal] = pd.to_numeric(group["n_asvs"],
                                              errors="coerce").fillna(0).sum()

    cross = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    for column in ("sample_id", "name_16S", "name_18S"):
        if column not in cross.columns:
            sys.exit("%s has no %s column" % (CROSSWALK, column))

    rows, unscoreable = [], {"no 16S library": [], "no ked marker": []}
    for _, record in cross.iterrows():
        animal = tick_lookup.resolve(table, "specimen", record["sample_id"])
        if animal is None:
            continue
        has_16S = not blank(record["name_16S"])
        ked_cols = [c for c in ("name_18S", "name_cytb", "name_ITS")
                    if c in cross.columns]
        has_ked_marker = any(not blank(record[c]) for c in ked_cols)
        if not has_16S:
            unscoreable["no 16S library"].append(animal)
            continue
        if not has_ked_marker:
            unscoreable["no ked marker"].append(animal)
            continue
        rows.append({
            "animal": animal,
            "host": "" if blank(record.get("host_bin")) else str(record["host_bin"]),
            "state": "" if blank(record.get("state")) else str(record["state"]),
            "n_ticks": record.get("num_ticks"),
            "ked_detected": animal in ked_reads,
            "ked_reads": ked_reads.get(animal, 0.0),
            "bartonella_detected": animal in bart_reads,
            "bartonella_reads": bart_reads.get(animal, 0.0),
            "bartonella_asvs": bart_asvs.get(animal, 0),
        })
    frame = pd.DataFrame(rows)
    if frame.empty:
        sys.exit("no specimen carries both a 16S library and a ked-capable "
                 "marker; there is nothing to test.")
    frame.to_csv(OUT_PATH, sep="\t", index=False)

    print()
    print("scoreable specimens %d, needing both a 16S library and one of %s"
          % (len(frame), ", ".join(KED_MARKERS)))
    for reason, names in unscoreable.items():
        print("   not scored, %-16s %d" % (reason + ":", len(names)))

    both = int((frame["ked_detected"] & frame["bartonella_detected"]).sum())
    ked_only = int((frame["ked_detected"] & ~frame["bartonella_detected"]).sum())
    bart_only = int((~frame["ked_detected"] & frame["bartonella_detected"]).sum())
    neither = int((~frame["ked_detected"] & ~frame["bartonella_detected"]).sum())

    print()
    print("=" * 70)
    print("DOES BARTONELLA TRACK THE KED?")
    print("=" * 70)
    print("                     ked detected   no ked")
    print("  Bartonella           %6d       %6d" % (both, bart_only))
    print("  no Bartonella        %6d       %6d" % (ked_only, neither))
    ked_total = both + ked_only
    no_ked_total = bart_only + neither
    if ked_total:
        print()
        print("  Bartonella in %d of %d with a ked (%.0f%%)"
              % (both, ked_total, 100.0 * both / ked_total))
    if no_ked_total:
        print("  Bartonella in %d of %d without one (%.0f%%)"
              % (bart_only, no_ked_total, 100.0 * bart_only / no_ked_total))
    p = fisher(both, bart_only, ked_only, neither)
    if p is None:
        print()
        print("  scipy is not importable here, so no exact test was run.")
    else:
        print("  Fisher exact, two-sided: p = %.4g" % p)

    if both or bart_only:
        print()
        print("  Bartonella reads where a ked was detected, and where not:")
        for label, mask in (("with ked", frame["ked_detected"]),
                            ("without", ~frame["ked_detected"])):
            part = frame.loc[mask & frame["bartonella_detected"], "bartonella_reads"]
            if len(part):
                print("    %-9s n %-4d median %-12.0f max %.0f"
                      % (label, len(part), part.median(), part.max()))
            else:
                print("    %-9s none" % label)

    odd = frame[frame["bartonella_detected"] & ~frame["ked_detected"]]
    print()
    print("Bartonella with no ked detected: %d specimens" % len(odd))
    if len(odd):
        print("  These need an explanation other than a ked in the tube.")
        for _, r in odd.sort_values("bartonella_reads", ascending=False).head(12).iterrows():
            print("    %-22s %-10s %-4s reads %-10.0f asvs %s"
                  % (r["animal"], r["host"], r["state"],
                     r["bartonella_reads"], r["bartonella_asvs"]))
    print()
    print("wrote " + OUT_PATH)


if __name__ == "__main__":
    main()
