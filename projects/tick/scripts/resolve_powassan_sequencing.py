#!/usr/bin/env python3
"""resolve_powassan_sequencing.py

I turn Lawrence's Powassan sequencing library names into animals.

WHY
---
Copy_of_Powassan_Virus_eDNA.xlsx reports 66 libraries: 19 with Powassan
virus reads (one of which is the extraction blank) and 47 with none. Those
names arrive in four surface forms, 'ALHT-Deer-02', 'EQ02442188-1',
'USGSTICKS-LHT-Deer-121-01' and 'Extraction-BLank-1', and the manuscript
counts animals, not libraries. I will not read a host or an animal off an
identifier string, so every name goes through tick_lookup, which is the
only sanctioned resolver in this project.

WHAT I DO
---------
For each sequencing library I try every marker in tick_lookup.MARKERS and
record which one matched and what specimen it gave. A name that matches
nothing is reported as absent, not quietly repaired: that is the contract
in tick_lookup's own docstring.

Many of these names carry a subsample suffix ('EQ02442188-1' is one of
seven subsamples of EQ02442188). If the full name is absent I try the
parent in a SEPARATE column, clearly labelled, so nothing silently folds a
subsample into its animal. I never strip anything else.

I then attach, per resolved specimen, the qPCR screening record's Powassan
positive-well count and subsample count, and the host from the crosswalk.
The qPCR sheet's host_id goes through the resolver too; it is stored
lowercased and in the name_18S spelling, which is not the crosswalk's
sample_id spelling.

READ FLOOR
----------
The extraction blank carries 21 reads. The lowest specimen carries 33,
which is 1.6x the blank, so I report the count both ways: every library
with any reads, and only those clearing 50 reads, the floor this project
already uses for 18S genus calls (about 2x the highest blank). I do not
pick one for you.

OUTPUT
------
A short summary to stdout and one table at
results/analysis/powassan_sequencing_resolved.tsv. I create no directories.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import tick_lookup

CROSSWALK_PATH = "results/metadata/crosswalk_4marker.tsv"
QPCR_PATH = "metadata/tick_full_metadata.xlsx"
QPCR_SHEET = "qPCR Results"
OUT_PATH = "results/analysis/powassan_sequencing_resolved.tsv"
READ_FLOOR = 50

# name, total Powassan reads (0 means no Powassan reads at all), top reference.
# Taken from the workbook's Amplicon feature table, which I recomputed and
# which agrees read-for-read with its Summary sheet.
SEQUENCING = [
    ("USGSTICKS-LHT-Deer-121-01", 1717260, "isolate DTVWiB08"),
    ("EQ02442188-6", 868890, "isolate P0375"),
    ("EQ02440618", 868389, "isolate 8"),
    ("EQ02442194", 689646, "isolate 18103-1419"),
    ("EQ02442188-1", 399144, "isolate P0760"),
    ("ALHT-Deer-15", 354792, "strain NFS001 isolate DTVMa96"),
    ("EQ02442227-4", 233048, "isolate 8"),
    ("ALHT-Deer-84", 228369, "isolate 8"),
    ("USGSTICKS-LHT-Deer-131-01", 184013, "isolate DTVWiB08"),
    ("USGSTICKS-LHT123-3", 137254, "isolate 19103-37-002"),
    ("USGSTICKS-2920-01", 37443, "isolate 8"),
    ("EQ02442188-7", 25734, "isolate P0375"),
    ("EQ02442220", 19202, "isolate LI327"),
    ("USGSTICKS-LHT123-2", 18409, "isolate 19103-37-002"),
    ("ALHT-Deer-02", 12524, "isolate 8"),
    ("eq02440635", 2897, "isolate 8"),
    ("ALHT-Deer-148-1", 150, "isolate 8"),
    ("EQ02442188-3", 33, "isolate P0375"),
    ("Extraction-BLank-1", 21, "isolate LI327"),
    ("ALHT-Deer-121-2", 0, ""),
    ("ALHT-Deer-126", 0, ""),
    ("ALHT-Deer-134", 0, ""),
    ("ALHT-Deer-140-2", 0, ""),
    ("ALHT-Deer-146-2", 0, ""),
    ("ALHT-Deer-32", 0, ""),
    ("ALHT-Deer-37", 0, ""),
    ("ALHT-Deer-53", 0, ""),
    ("ALHT-Deer-55", 0, ""),
    ("ALHT-Deer-80", 0, ""),
    ("eq02440612", 0, ""),
    ("eq02440616-2", 0, ""),
    ("eq02440624", 0, ""),
    ("eq02440626", 0, ""),
    ("eq02440649", 0, ""),
    ("eq02440775-1", 0, ""),
    ("EQ02441291-2", 0, ""),
    ("EQ02441291-3", 0, ""),
    ("EQ02441291-4", 0, ""),
    ("EQ02441291-5", 0, ""),
    ("EQ02441291-6", 0, ""),
    ("EQ02441305", 0, ""),
    ("EQ02441306-1", 0, ""),
    ("EQ02441478-1", 0, ""),
    ("EQ02441506", 0, ""),
    ("EQ02442188-2", 0, ""),
    ("EQ02442188-4", 0, ""),
    ("EQ02442188-5", 0, ""),
    ("EQ02442197-2", 0, ""),
    ("EQ02442197-3", 0, ""),
    ("EQ02442218-1", 0, ""),
    ("EQ02442227-5", 0, ""),
    ("USGS-tick-0601", 0, ""),
    ("USGS-tick-0603", 0, ""),
    ("USGSTICKS-4101-1", 0, ""),
    ("USGSTICKS-8404-03", 0, ""),
    ("USGSTICKS-LHT094-1", 0, ""),
    ("USGSTICKS-LHT105-1", 0, ""),
    ("USGSTICKS-LHT105-2", 0, ""),
    ("USGSTICKS-LHT109-1", 0, ""),
    ("USGSTICKS-LHT113-1", 0, ""),
    ("USGSTICKS-LHT113-2", 0, ""),
    ("USGSTICKS-LHT123-1", 0, ""),
    ("USGSTICKS-LHT-Deer-007-01", 0, ""),
    ("USGSTICKS-LHT-Deer-010-01", 0, ""),
    ("USGSTICKS-LHT-Deer-168-01", 0, ""),
    ("USGSTICKS-LHT-Deer-170-01", 0, ""),
]


def parent_name(name):
    """I drop one trailing subsample suffix, and only that.

    'EQ02442188-1' becomes 'EQ02442188'. A name whose last hyphen-delimited
    piece is not all digits comes back unchanged, so I never chew into a
    name looking for a match.
    """
    pieces = name.rsplit("-", 1)
    if len(pieces) == 2 and pieces[1].isdigit():
        return pieces[0]
    return name


def resolve_every_marker(table, name):
    """I try one name under every marker and return (specimen, markers)."""
    hits = {}
    for marker in tick_lookup.MARKERS:
        specimen = tick_lookup.resolve(table, marker, name)
        if specimen is not None:
            hits[marker] = specimen
    if not hits:
        return None, []
    found = sorted(set(hits.values()))
    if len(found) > 1:
        sys.exit("'" + name + "' resolves to more than one specimen: "
                 + ", ".join(found)
                 + "\nThe lookup is inconsistent. Rebuild the crosswalk:"
                 + "\n    python3 scripts/build_4marker_crosswalk.py")
    return found[0], sorted(hits)


def host_columns(frame):
    """I require both host columns the crosswalk carries and report both.

    host_spp is the finer label, matching the manuscript's porcupine and
    opossum rows. host_bin is the five-way binning build_species_figure.py
    uses. They disagree for the small host rows, so I carry both rather
    than choose one and have the choice disappear into a single column.
    """
    missing = [c for c in ("host_spp", "host_bin") if c not in frame.columns]
    if missing:
        sys.exit("column(s) " + ", ".join(missing) + " missing from "
                 + CROSSWALK_PATH
                 + "\ncolumns present: " + ", ".join(map(str, frame.columns)))
    return "host_spp", "host_bin"


def main():
    if not os.path.isdir("results"):
        sys.exit("run me from the project root; there is no results/ here.")
    for path in (CROSSWALK_PATH, QPCR_PATH):
        if not os.path.isfile(path):
            sys.exit("missing input: " + path)

    table = tick_lookup.load()

    crosswalk = pd.read_csv(CROSSWALK_PATH, sep="\t", dtype=str)
    if "sample_id" not in crosswalk.columns:
        sys.exit("no sample_id column in " + CROSSWALK_PATH
                 + "\ncolumns present: " + ", ".join(map(str, crosswalk.columns)))
    spp_col, bin_col = host_columns(crosswalk)
    spp_of, bin_of = {}, {}
    unresolved_crosswalk = 0
    for sample_id, spp, binned in zip(crosswalk["sample_id"],
                                      crosswalk[spp_col],
                                      crosswalk[bin_col]):
        specimen = tick_lookup.resolve(table, "specimen", sample_id)
        if specimen is None:
            unresolved_crosswalk += 1
            continue
        spp_of[specimen] = "" if pd.isna(spp) else str(spp).strip()
        bin_of[specimen] = "" if pd.isna(binned) else str(binned).strip()

    qpcr = pd.read_excel(QPCR_PATH, sheet_name=QPCR_SHEET, dtype=str)
    for column in ("host_id", "powassan-virus", "subsample_count"):
        if column not in qpcr.columns:
            sys.exit("column '" + column + "' missing from " + QPCR_PATH
                     + " sheet '" + QPCR_SHEET + "'"
                     + "\ncolumns present: " + ", ".join(map(str, qpcr.columns)))
    qpcr_of = {}
    unresolved_qpcr = []
    for host_id, powassan, subsamples in zip(qpcr["host_id"],
                                             qpcr["powassan-virus"],
                                             qpcr["subsample_count"]):
        specimen = tick_lookup.resolve(table, "specimen", host_id)
        if specimen is None:
            if not pd.isna(host_id):
                unresolved_qpcr.append(str(host_id))
            continue
        qpcr_of[specimen] = (
            "" if pd.isna(powassan) else str(powassan).strip(),
            "" if pd.isna(subsamples) else str(subsamples).strip(),
        )

    rows = []
    for name, reads, reference in SEQUENCING:
        specimen, markers = resolve_every_marker(table, name)
        parent, parent_markers = None, []
        if specimen is None:
            candidate = parent_name(name)
            if candidate != name:
                parent, parent_markers = resolve_every_marker(table, candidate)
        animal = specimen if specimen is not None else parent
        powassan, subsamples = qpcr_of.get(animal, ("", "")) if animal else ("", "")
        rows.append({
            "library": name,
            "powassan_reads": reads,
            "top_reference": reference,
            "sequencing_call": "positive" if reads > 0 else "negative",
            "clears_floor": "yes" if reads >= READ_FLOOR else "no",
            "specimen": specimen if specimen is not None else "",
            "matched_markers": ",".join(markers),
            "parent_tried": parent_name(name) if specimen is None else "",
            "parent_specimen": parent if parent is not None else "",
            "parent_markers": ",".join(parent_markers),
            "host_spp": spp_of.get(animal, "") if animal else "",
            "host_bin": bin_of.get(animal, "") if animal else "",
            "qpcr_powassan_wells": powassan,
            "qpcr_subsample_count": subsamples,
        })

    frame = pd.DataFrame(rows)
    frame.to_csv(OUT_PATH, sep="\t", index=False)

    positives = frame[frame["sequencing_call"] == "positive"]
    blank = positives[positives["library"].str.lower().str.startswith("extraction")]
    specimens = positives[~positives["library"].isin(blank["library"])]
    cleared = specimens[specimens["clears_floor"] == "yes"]

    def animals(sub):
        found = set(sub["specimen"]) | set(sub["parent_specimen"])
        found.discard("")
        return found

    print("crosswalk rows whose sample_id did not resolve: %d of %d"
          % (unresolved_crosswalk, len(crosswalk)))
    print("qPCR host_id values that did not resolve: %d%s"
          % (len(unresolved_qpcr),
             (" (" + ", ".join(unresolved_qpcr[:5]) + ")") if unresolved_qpcr else ""))
    print()
    print("libraries: %d  sequencing-positive: %d  negative: %d  blank in positives: %d"
          % (len(frame), len(positives), len(frame) - len(positives), len(blank)))
    print("positive specimens (blank excluded): %d, of which clearing %d reads: %d"
          % (len(specimens), READ_FLOOR, len(cleared)))
    print()
    print("resolved on the full library name: %d of %d libraries"
          % (int((frame["specimen"] != "").sum()), len(frame)))
    print("resolved only via the parent name:  %d"
          % int(((frame["specimen"] == "") & (frame["parent_specimen"] != "")).sum()))
    print("not resolved at all:                %d  %s"
          % (int(((frame["specimen"] == "") & (frame["parent_specimen"] == "")).sum()),
             ", ".join(frame.loc[(frame["specimen"] == "")
                                 & (frame["parent_specimen"] == ""), "library"])))
    print()
    print("distinct animals, sequencing-positive specimens: %d" % len(animals(specimens)))
    print("distinct animals, clearing the read floor:       %d" % len(animals(cleared)))
    print("distinct animals, sequencing-negative:           %d"
          % len(animals(frame[frame["sequencing_call"] == "negative"])))
    print()
    print("positive specimens against the qPCR record:")
    for _, row in specimens.sort_values("powassan_reads", ascending=False).iterrows():
        print("  %-26s %9d reads  animal %-16s %-22s [%s]  qPCR powassan wells %-4s of %s"
              % (row["library"], row["powassan_reads"],
                 row["specimen"] or row["parent_specimen"] or "UNRESOLVED",
                 row["host_spp"] or "-",
                 row["host_bin"] or "-",
                 row["qpcr_powassan_wells"] or "-",
                 row["qpcr_subsample_count"] or "-"))
    print()
    print("wrote " + OUT_PATH)


if __name__ == "__main__":
    main()
