#!/usr/bin/env python3
"""build_table1_sample_set.py

I build Table 1, the sample set by host: how many animals, how many ticks, what
tick genera they carried, and which assays they went through.

WHY THIS TABLE
--------------
Every cohort in the paper is a subset of these rows. A reader should be able to
get from 563 animals to the 330 deer carrying Ixodes and the 11 moose carrying
Dermacentor without being asked to take either number on faith. This table is
also what replaced the tick-species-by-host figure in the restructure.

WHAT THE COLUMNS MEAN
---------------------
    Animals        source animals, the unit of analysis throughout
    Ticks          ticks collected, from the crosswalk rather than the run
                   sheet, whose Num Ticks is corrupted by the barcode clobber
    Ixodes         animals with at least one Ixodes identification
    Dermacentor    animals with at least one Dermacentor identification
    qPCR           animals with a qPCR result
    16S            animals with a 16S library

The Ixodes and Dermacentor columns overlap, because an animal can carry both,
and they do not sum to the animal count, because some animals have no species
call. That is stated in the caption rather than hidden by a made-up category.

WHAT I READ
-----------
    results/metadata/crosswalk_4marker.tsv
    results/analysis/tick_species_by_host.tsv
    results/analysis/qpcr_screening_record_corrected.tsv
      (the 'qPCR Results' sheet, recounted under Section 2.5 on 2026-10-08)

WHAT I WRITE, WITH --write
--------------------------
    results/analysis/table1_sample_set.tsv

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/build_table1_sample_set.py
    python3 scripts/build_table1_sample_set.py --write
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
SPECIES = "results/analysis/tick_species_by_host.tsv"
BOOK = "metadata/tick_full_metadata.xlsx"
SHEET = "qPCR Results"   # supplied record; read via QPCR below
QPCR = "results/analysis/qpcr_screening_record_corrected.tsv"
OUT = "results/analysis/table1_sample_set.tsv"
TRUE_WORDS = ("1", "1.0", "true", "yes", "y")


def truthy(value):
    return str(value).strip().lower() in TRUE_WORDS


def grouper():
    """I return check_qpcr_by_host, for the host grouping already agreed."""
    try:
        import check_qpcr_by_host as module
    except ImportError as problem:
        sys.exit("cannot import check_qpcr_by_host (%s). It must sit in "
                 "scripts/ beside me." % problem)
    for wanted in ("host_group", "HOST_GROUPS"):
        if not hasattr(module, wanted):
            sys.exit("check_qpcr_by_host has no %s. Tell me what it calls the "
                     "host grouping and I will use that instead." % wanted)
    return module


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s)." % problem)
    return tick_lookup.load(), tick_lookup


def resolve(table, module, value):
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    hit = module.resolve(table, "specimen", text)
    return "" if hit is None else str(hit)


def group_of(module, spp, binned):
    """I return the host group, taking the richest text available.

    host_bin is already binned, and host_group matches free text against word
    lists that do not include its own group names, so an animal whose
    host_bin reads 'Other wildlife' falls through to 'Other or unrecorded'
    unless I check for a group name first. host_spp carries the original
    collection text and is tried before host_bin for the same reason.
    """
    names = {}
    for entry in module.HOST_GROUPS:
        names[str(entry[0]).strip().lower()] = entry[0]
    for text in (spp, binned):
        value = str(text).strip()
        if not value or value.lower() == "nan":
            continue
        if value.lower() in names:
            return names[value.lower()]
        hit = module.host_group(value)
        if hit and hit != module.UNMAPPED:
            return hit
    return "Unknown"


def order_of(module):
    """I return the host group names in the order they should be printed."""
    names = []
    for entry in module.HOST_GROUPS:
        names.append(entry[0] if isinstance(entry, (list, tuple)) else entry)
    return names


def main():
    write = "--write" in sys.argv[1:]
    for path in (CROSSWALK, SPECIES, BOOK):
        if not os.path.isfile(path):
            sys.exit("missing %s. Run me from the project root." % path)
    module = grouper()
    table, finder = lookup()

    crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    for column in ("sample_id", "barcode", "host_bin", "num_ticks",
                   "name_16S"):
        if column not in crosswalk.columns:
            sys.exit("%s has no '%s' column" % (CROSSWALK, column))
    crosswalk["animal"] = [resolve(table, finder, b) or resolve(table, finder, s)
                           for b, s in zip(crosswalk["barcode"],
                                           crosswalk["sample_id"])]
    if "host_spp" not in crosswalk.columns:
        sys.exit("%s has no 'host_spp' column" % CROSSWALK)
    crosswalk["group"] = [group_of(module, spp, binned) for spp, binned
                          in zip(crosswalk["host_spp"], crosswalk["host_bin"])]

    species = pd.read_csv(SPECIES, sep="\t", dtype=str)
    for column in ("specimen", "Ixodes_pos", "Dermacentor_pos"):
        if column not in species.columns:
            sys.exit("%s has no '%s' column" % (SPECIES, column))
    ixodes, derm = set(), set()
    for specimen, has_ixodes, has_derm in zip(species["specimen"],
                                              species["Ixodes_pos"],
                                              species["Dermacentor_pos"]):
        animal = resolve(table, finder, specimen)
        if not animal:
            continue
        if truthy(has_ixodes):
            ixodes.add(animal)
        if truthy(has_derm):
            derm.add(animal)

    book = pd.read_csv(QPCR, sep="\t", dtype=str)
    columns = {str(c).strip().lower(): c for c in book.columns}
    screened = set()
    for _, row in book.iterrows():
        for key in ("barcode", "host_id"):
            if key in columns:
                animal = resolve(table, finder, row[columns[key]])
                if animal:
                    screened.add(animal)
                    break

    rows = []
    order = order_of(module) + ["Unknown", getattr(module, "UNMAPPED", ""), ""]
    present = [g for g in order if g in set(crosswalk["group"])]
    present += [g for g in sorted(set(crosswalk["group"])) if g not in present]
    for group in present:
        part = crosswalk[crosswalk["group"] == group]
        animals = set(a for a in part["animal"] if a)
        ticks = pd.to_numeric(part["num_ticks"], errors="coerce").sum()
        sequenced = set(a for a, n in zip(part["animal"], part["name_16S"])
                        if a and pd.notna(n) and str(n).strip()
                        and str(n).strip().lower() != "nan")
        rows.append({
            "Host": group or "Not recorded",
            "Animals": len(animals),
            "Ticks": int(ticks) if pd.notna(ticks) else 0,
            "Ixodes": len(animals & ixodes),
            "Dermacentor": len(animals & derm),
            "qPCR": len(animals & screened),
            "16S": len(sequenced),
        })

    frame = pd.DataFrame(rows).sort_values("Animals", ascending=False)
    total = {"Host": "All hosts"}
    for column in ("Animals", "Ticks", "Ixodes", "Dermacentor", "qPCR", "16S"):
        total[column] = int(frame[column].sum())
    frame = pd.concat([frame, pd.DataFrame([total])], ignore_index=True)

    widths = {c: max(len(str(c)), *(len(str(v)) for v in frame[c]))
              for c in frame.columns}
    header = "  ".join(str(c).ljust(widths[c]) if c == "Host"
                       else str(c).rjust(widths[c]) for c in frame.columns)
    print(header)
    print("-" * len(header))
    for _, row in frame.iterrows():
        print("  ".join(str(row[c]).ljust(widths[c]) if c == "Host"
                        else str(row[c]).rjust(widths[c])
                        for c in frame.columns))
    print()

    deer_only = len({a for a in ixodes if a not in derm}
                    & set(crosswalk[crosswalk["group"].astype(str)
                                    .str.contains("deer", case=False,
                                                  na=False)]["animal"]))
    print("  deer carrying Ixodes and not Dermacentor: %d" % deer_only)
    print("  This counts specimens carrying a genus call. Section 3.3's")
    print("  cohort sits three filters downstream of it -- nameable by the")
    print("  identifier resolver, then two read floors -- and is computed")
    print("  by the core microbiome builder, not here, so the two are not")
    print("  expected to match. The funnel is Table 2's row for the core")
    print("  microbiome in Ixodes from deer.")
    if deer_only != 371:
        print("  BUT 371 was expected here. The genus calls have moved,")
        print("  so every count downstream of them needs recomputing.")
    print()

    if not write:
        print("Nothing written. Rerun with --write to save %s" % OUT)
        return
    if os.path.exists(OUT):
        sys.exit("%s already exists. Move it aside first." % OUT)
    frame.to_csv(OUT, sep="\t", index=False)
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
