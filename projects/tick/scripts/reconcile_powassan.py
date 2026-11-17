#!/usr/bin/env python3
"""reconcile_powassan.py

I reconcile the Powassan qPCR screening record against Lawrence's
sequencing run, one row per animal.

WHY
---
resolve_powassan_sequencing.py established that 14 animals carry Powassan
virus NS5 reads and that 6 of them have no positive qPCR well. That is one
direction of the disagreement. I build the other direction too: which
animals the qPCR called positive, whether each was sequenced at all, and
what sequencing said. Without that, "14 confirmed" could be hiding qPCR
positives that sequencing failed to reproduce.

I also attach, for every animal in the table, the state and the 18S tick
genus calls, because the manuscript needs to report where each detection
was and what tick carried it, and neither can be read off an identifier.

WHAT I WILL NOT DO
------------------
I do not parse names. Every qPCR host_id and every genus-table specimen
goes through tick_lookup. The sequencing side is read from the table the
previous script wrote, which already resolved its names the same way.

The qPCR target columns are COUNTS of positive subsamples, not 0/1 flags,
so I parse them as numbers and never test a string against "1". The
genus-table positive flags are read against an explicit set of accepted
values and I exit if I meet one I do not recognise, rather than treating
the unknown as false.

OUTPUT
------
A summary to stdout and one table at
results/analysis/powassan_reconciliation.tsv. I create no directories.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import tick_lookup

SEQ_PATH = "results/analysis/powassan_sequencing_resolved.tsv"
GENUS_PATH = "results/analysis/tick_species_by_host.tsv"
QPCR_PATH = "metadata/tick_full_metadata.xlsx"
QPCR_SHEET = "qPCR Results"
OUT_PATH = "results/analysis/powassan_reconciliation.tsv"
READ_FLOOR = 50

TRUE_VALUES = ("1", "1.0", "true", "yes", "y", "t")
FALSE_VALUES = ("0", "0.0", "false", "no", "n", "f", "", "nan", "none")


def count_of(value):
    """I read a positive-subsample count as a number, or 0 if it is blank.

    The screening record stores counts, so '2' and '2.0' are both two. A
    blank means the animal was not screened for that target and I return
    None so the caller can tell blank from zero.
    """
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if text.lower() in ("", "nan", "none"):
        return None
    try:
        return int(float(text))
    except ValueError:
        sys.exit("cannot read '" + text + "' as a positive-subsample count in "
                 + QPCR_PATH + " sheet '" + QPCR_SHEET + "'")


def flag_of(value, column, specimen):
    """I read a genus positive flag against an explicit vocabulary."""
    text = "" if value is None or pd.isna(value) else str(value).strip().lower()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    sys.exit("unrecognised value '" + str(value) + "' in column '" + column
             + "' for specimen '" + str(specimen) + "' of " + GENUS_PATH
             + "\nI accept " + ", ".join(TRUE_VALUES) + " as true and "
             + ", ".join(v for v in FALSE_VALUES if v) + " as false."
             + "\nI will not guess which one this is.")


def main():
    if not os.path.isdir("results"):
        sys.exit("run me from the project root; there is no results/ here.")
    for path in (SEQ_PATH, GENUS_PATH, QPCR_PATH):
        if not os.path.isfile(path):
            sys.exit("missing input: " + path
                     + ("\nRun scripts/resolve_powassan_sequencing.py first."
                        if path == SEQ_PATH else ""))

    table = tick_lookup.load()

    # Sequencing side, per animal, from the table the previous script wrote.
    seq = pd.read_csv(SEQ_PATH, sep="\t", dtype=str)
    for column in ("library", "powassan_reads", "specimen", "parent_specimen"):
        if column not in seq.columns:
            sys.exit("column '" + column + "' missing from " + SEQ_PATH
                     + "\nRerun scripts/resolve_powassan_sequencing.py.")
    seq["reads"] = seq["powassan_reads"].astype(float).astype(int)
    seq["animal"] = seq["specimen"].fillna("")
    blanks = []
    sequenced = {}
    for _, row in seq.iterrows():
        animal = str(row["specimen"]) if not pd.isna(row["specimen"]) else ""
        if not animal:
            animal = str(row["parent_specimen"]) if not pd.isna(row["parent_specimen"]) else ""
        if not animal:
            blanks.append(str(row["library"]))
            continue
        record = sequenced.setdefault(animal, {"libraries": 0, "positive": 0,
                                               "cleared": 0, "reads": 0,
                                               "names": []})
        record["libraries"] += 1
        record["reads"] += int(row["reads"])
        record["names"].append(str(row["library"]))
        if int(row["reads"]) > 0:
            record["positive"] += 1
        if int(row["reads"]) >= READ_FLOOR:
            record["cleared"] += 1

    # qPCR side, per animal.
    qpcr = pd.read_excel(QPCR_PATH, sheet_name=QPCR_SHEET, dtype=str)
    for column in ("host_id", "powassan-virus", "subsample_count"):
        if column not in qpcr.columns:
            sys.exit("column '" + column + "' missing from " + QPCR_PATH
                     + " sheet '" + QPCR_SHEET + "'"
                     + "\ncolumns present: " + ", ".join(map(str, qpcr.columns)))
    screened = {}
    unresolved = []
    for host_id, powassan, subsamples in zip(qpcr["host_id"],
                                             qpcr["powassan-virus"],
                                             qpcr["subsample_count"]):
        wells = count_of(powassan)
        count = count_of(subsamples)
        specimen = tick_lookup.resolve(table, "specimen", host_id)
        if specimen is None:
            if not pd.isna(host_id):
                unresolved.append((str(host_id), wells, count))
            continue
        screened[specimen] = (wells, count)

    # 18S tick genus and state.
    genus = pd.read_csv(GENUS_PATH, sep="\t", dtype=str)
    needed = ("specimen", "state", "informative", "Ixodes_pos",
              "Dermacentor_pos", "Amblyomma_pos")
    missing = [c for c in needed if c not in genus.columns]
    if missing:
        sys.exit("column(s) " + ", ".join(missing) + " missing from " + GENUS_PATH
                 + "\ncolumns present: " + ", ".join(map(str, genus.columns)))
    genus_of = {}
    unresolved_genus = 0
    for _, row in genus.iterrows():
        specimen = tick_lookup.resolve(table, "specimen", row["specimen"])
        if specimen is None:
            unresolved_genus += 1
            continue
        found = [name for name, column in (("Ixodes", "Ixodes_pos"),
                                           ("Dermacentor", "Dermacentor_pos"),
                                           ("Amblyomma", "Amblyomma_pos"))
                 if flag_of(row[column], column, row["specimen"])]
        genus_of[specimen] = (
            "" if pd.isna(row["state"]) else str(row["state"]).strip(),
            ", ".join(found),
            flag_of(row["informative"], "informative", row["specimen"]),
        )

    animals = sorted(set(sequenced) | set(k for k, v in screened.items()
                                          if v[0] not in (None, 0)))
    rows = []
    for animal in animals:
        seq_record = sequenced.get(animal)
        wells, subsamples = screened.get(animal, (None, None))
        state, genera, informative = genus_of.get(animal, ("", "", False))
        if seq_record is None:
            call = "not sequenced"
        elif seq_record["positive"]:
            call = "positive"
        else:
            call = "negative"
        rows.append({
            "animal": animal,
            "sequencing_call": call,
            "libraries_sequenced": seq_record["libraries"] if seq_record else 0,
            "libraries_positive": seq_record["positive"] if seq_record else 0,
            "libraries_clearing_floor": seq_record["cleared"] if seq_record else 0,
            "powassan_reads": seq_record["reads"] if seq_record else 0,
            "libraries": ";".join(seq_record["names"]) if seq_record else "",
            "qpcr_powassan_wells": "" if wells is None else wells,
            "qpcr_subsample_count": "" if subsamples is None else subsamples,
            "state": state,
            "tick_genus_18S": genera,
            "genus_informative": "yes" if informative else "no",
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT_PATH, sep="\t", index=False)

    def subset(call, positive_qpcr):
        out = []
        for row in rows:
            wells = row["qpcr_powassan_wells"]
            has = isinstance(wells, int) and wells > 0
            if row["sequencing_call"] == call and has == positive_qpcr:
                out.append(row)
        return out

    qpcr_positive_animals = [k for k, v in screened.items()
                             if v[0] not in (None, 0)]

    print("libraries with no animal (expected: the extraction blank): %d  %s"
          % (len(blanks), ", ".join(blanks)))
    print("genus-table specimens that did not resolve: %d of %d"
          % (unresolved_genus, len(genus)))
    print()
    print("animals with at least one positive Powassan qPCR well: %d"
          % len(qpcr_positive_animals))
    print("animals with Powassan sequencing reads: %d"
          % len([r for r in rows if r["sequencing_call"] == "positive"]))
    print()
    print("2x2 over animals that were sequenced")
    print("                        seq positive   seq negative")
    print("  qPCR well(s) >0       %12d   %12d"
          % (len(subset("positive", True)), len(subset("negative", True))))
    print("  qPCR wells == 0       %12d   %12d"
          % (len(subset("positive", False)), len(subset("negative", False))))
    print("  qPCR positive animals never sequenced: %d"
          % len(subset("not sequenced", True)))
    print()
    print("qPCR-positive animals and what sequencing said")
    for row in sorted(rows, key=lambda r: str(r["animal"])):
        wells = row["qpcr_powassan_wells"]
        if isinstance(wells, int) and wells > 0:
            print("  %-24s wells %-4s of %-4s  sequencing %-14s reads %-9d state %-4s genus %s"
                  % (row["animal"], wells, row["qpcr_subsample_count"],
                     row["sequencing_call"], row["powassan_reads"],
                     row["state"] or "-", row["tick_genus_18S"] or "none"))
    print()
    print("sequence-confirmed animals with no positive qPCR well")
    for row in sorted(rows, key=lambda r: -int(r["powassan_reads"])):
        wells = row["qpcr_powassan_wells"]
        if row["sequencing_call"] == "positive" and not (isinstance(wells, int) and wells > 0):
            print("  %-24s reads %-9d wells %-4s of %-4s  state %-4s genus %s"
                  % (row["animal"], row["powassan_reads"],
                     "blank" if wells == "" else wells,
                     row["qpcr_subsample_count"] or "-",
                     row["state"] or "-", row["tick_genus_18S"] or "none"))
    print()
    print("unresolved qPCR host_id rows (host_id, powassan wells, subsamples): %d"
          % len(unresolved))
    for host_id, wells, count in unresolved:
        print("  %-20s wells %-6s subsamples %s"
              % (host_id, "blank" if wells is None else wells,
                 "blank" if count is None else count))
    print()
    print("wrote " + OUT_PATH)


if __name__ == "__main__":
    main()
