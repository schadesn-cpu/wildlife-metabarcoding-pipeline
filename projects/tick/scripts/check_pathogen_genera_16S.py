#!/usr/bin/env python3
"""check_pathogen_genera_16S.py

I ask whether genera the qPCR panel covers, or that sit on the
endosymbiont-pathogen boundary, appear in the 16S data at all.

WHY
---
Three questions, each answerable either way:

  Ehrlichia    The qPCR screened for E. muris and never detected it. If 16S
               also finds nothing, that is a concordance worth one sentence.
               If 16S finds reads, the null in Section 3.2 needs qualifying.
  Coxiella     Named in the Introduction as a genus where the endosymbiont
               and pathogen boundary has been crossed. If it is here, the
               Introduction's example is also a result.
  Diplorickettsia, Rickettsiella
               Both in family Diplorickettsiaceae. The deer table reports
               the family rather than either genus, which is the paper's own
               argument about resolution happening in its own output. I want
               to know whether the genus-level rows carry anything.

Nocardia is included only because it was raised. It is a soil actinomycete
and I expect any reads to be surface or reagent background, which is what
the control columns are for.

Anaplasma and Francisella are included as method controls. Both are known to
be abundant here, so if my sample-column resolution or my read counting is
broken, those two will show it rather than a silent zero on the targets.

WHAT I DO NOT DO
----------------
I do not parse sample names beyond the one split this project already
sanctions. The taxonomy table's columns are '<library>.<batch>', and
build_core_microbiome.py line 51 resolves that form by taking the part
before the first dot and passing it to tick_lookup. I do the same, in
resolve_library below. Everything after that split is the resolver's job.

A column that does not resolve is reported by name rather than dropped. On
the run of 2026-10-06 that was 54 of 907: 48 controls, four specimens the
manuscript excludes by name, and two, ALHT-Deer-10 and EQ02442982, that have
a 16S library and a qPCR row but no crosswalk entry. The second group is
ground truth section 8's watch condition and is not mine to resolve here.

These are classifier labels, not BLAST-verified calls. Section 2.4 verified
223 pathogen-labelled ASVs; anything I report here that matters would need
the same treatment before it goes in the manuscript.

OUTPUT
------
A summary to stdout and one table at
results/analysis/pathogen_genera_16S_check.tsv. I create no directories.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import tick_lookup

COUNTS_PATH = "results/16S/all/taxonomy/taxonomy_counts_L6_16S.tsv"
CROSSWALK_PATH = "results/metadata/crosswalk_4marker.tsv"
QPCR_PATH = "metadata/tick_full_metadata.xlsx"
QPCR_SHEET = "qPCR Results"
OUT_PATH = "results/analysis/pathogen_genera_16S_check.tsv"

READ_FLOOR = 50
FRACTION_FLOOR = 0.01

TARGETS = ("Ehrlichia", "Coxiella", "Diplorickettsia", "Rickettsiella",
           "Nocardia", "Anaplasma", "Francisella")


def genus_of(lineage):
    """I take the last rank of a SILVA lineage string, whatever its depth."""
    return str(lineage).split(";")[-1].strip()


def resolve_library(table, column):
    """I resolve one taxonomy column to an animal, as build_core_microbiome does.

    The taxonomy table's columns are '<library>.<batch>' and the lookup does
    not hold that composite form, so resolving the column verbatim returns
    None for all 907 of them. build_core_microbiome.py line 51 takes the part
    before the first dot and resolves that through tick_lookup, and this is
    the same file, so I follow it rather than give one table two answers.

    The split is on the batch separator the pipeline writes, not a search for
    something that happens to match. If the stem is absent from the lookup I
    return None and the caller reports it, which is the contract in
    tick_lookup's own docstring.
    """
    stem = str(column).split(".")[0]
    for marker in ("16S", "specimen"):
        hit = tick_lookup.resolve(table, marker, stem)
        if hit is not None:
            return hit
    return None


def main():
    if not os.path.isdir("results"):
        sys.exit("run me from the project root; there is no results/ here.")
    for path in (COUNTS_PATH, CROSSWALK_PATH, QPCR_PATH):
        if not os.path.isfile(path):
            sys.exit("missing input: " + path)

    table = tick_lookup.load()

    counts = pd.read_csv(COUNTS_PATH, sep="\t", dtype=str)
    first = counts.columns[0]
    samples = [c for c in counts.columns if c != first]
    if not samples:
        sys.exit(COUNTS_PATH + " has no sample columns besides '" + first + "'")
    for column in samples:
        counts[column] = pd.to_numeric(counts[column], errors="coerce").fillna(0)
    counts["genus"] = counts[first].map(genus_of)

    resolved, unresolved = {}, []
    for column in samples:
        specimen = resolve_library(table, column)
        if specimen is None:
            unresolved.append(column)
        else:
            resolved[column] = specimen

    crosswalk = pd.read_csv(CROSSWALK_PATH, sep="\t", dtype=str)
    host_of = {}
    for sample_id, host in zip(crosswalk["sample_id"], crosswalk["host_bin"]):
        specimen = tick_lookup.resolve(table, "specimen", sample_id)
        if specimen is not None:
            host_of[specimen] = "" if pd.isna(host) else str(host).strip()

    qpcr = pd.read_excel(QPCR_PATH, sheet_name=QPCR_SHEET, dtype=str)
    muris_of = {}
    for host_id, muris in zip(qpcr["host_id"], qpcr["e-muris"]):
        specimen = tick_lookup.resolve(table, "specimen", host_id)
        if specimen is None or pd.isna(muris):
            continue
        try:
            muris_of[specimen] = int(float(str(muris).strip()))
        except ValueError:
            sys.exit("cannot read e-muris value '" + str(muris) + "' as a count")

    specimen_columns = sorted(resolved)
    totals = counts[specimen_columns].sum(axis=0)

    print("taxonomy rows: %d   sample columns: %d" % (len(counts), len(samples)))
    print("columns resolving to a specimen: %d   not resolving: %d"
          % (len(resolved), len(unresolved)))
    if unresolved:
        print("  unresolved (treated as blanks and controls):")
        for column in unresolved:
            print("    %-46s %d reads total"
                  % (column, int(counts[column].sum())))
    print()

    rows = []
    for target in TARGETS:
        hits = counts[counts["genus"] == target]
        if hits.empty:
            print("%-16s NO ROW in the taxonomy at all" % target)
            continue
        per_specimen = hits[specimen_columns].sum(axis=0)
        control = int(hits[unresolved].sum(axis=0).sum()) if unresolved else 0
        any_read = per_specimen[per_specimen > 0]
        cleared = per_specimen[per_specimen >= READ_FLOOR]
        fraction = per_specimen / totals.replace(0, pd.NA)
        both = per_specimen[(per_specimen >= READ_FLOOR)
                            & (fraction >= FRACTION_FLOOR)]
        animals = set(resolved[c] for c in any_read.index)
        print("%-16s rows %-3d reads %-12d specimens >0 %-4d >=%d reads %-4d "
              "and >=%.0f%% %-4d controls %d"
              % (target, len(hits), int(per_specimen.sum()), len(any_read),
                 READ_FLOOR, len(cleared), FRACTION_FLOOR * 100, len(both),
                 control))
        for column, value in cleared.sort_values(ascending=False).items():
            specimen = resolved[column]
            rows.append({
                "genus": target,
                "library": column,
                "specimen": specimen,
                "host_bin": host_of.get(specimen, ""),
                "reads": int(value),
                "fraction_of_specimen": round(float(value) / float(totals[column]), 5)
                                        if totals[column] else "",
                "qpcr_e_muris_wells": muris_of.get(specimen, ""),
            })
        del animals

    if rows:
        pd.DataFrame(rows).to_csv(OUT_PATH, sep="\t", index=False)
        print()
        print("specimens clearing %d reads, by genus and host:" % READ_FLOOR)
        frame = pd.DataFrame(rows)
        print(pd.crosstab(frame["genus"], frame["host_bin"]).to_string())
        print()
        print("wrote " + OUT_PATH)
    else:
        print()
        print("no specimen cleared %d reads for any target genus; nothing written"
              % READ_FLOOR)


if __name__ == "__main__":
    main()
