#!/usr/bin/env python3
"""check_anaplasma_batch.py

I ask whether the Anaplasma signal is a property of the animals or of the
sequencing batch.

WHY
---
Two numbers do not sit well together. Paragraph 73 establishes the paper's
control test on 29 negative controls holding 139,681 reads. Counting every
control column in the 16S genus table gives 48 columns carrying 53,441
Anaplasma reads, and NTC.250609_Batch-11 alone holds 58,695 reads of
everything. Batch-11 is the 89-animal Vermont-only batch, and Vermont is
where A. phagocytophilum prevalence is highest in the qPCR.

That could be three different things: a narrower control set in paragraph
73 than the one I counted, a genuinely contaminated NTC in one batch, or an
Anaplasma signal that tracks batch rather than animal. The third would
matter, because A. phagocytophilum at 61% of deer pools is the paper's
largest pathogen number.

I do not try to adjudicate. I lay the three out so the question can be
answered by looking.

WHAT I DO
---------
Controls: every taxonomy column that does not resolve to an animal, with
its total reads, its Anaplasma reads and the fraction. Named individually,
because which control is dirty is the whole question.

Specimens: Anaplasma reads and fraction per library, grouped by sequencing
batch, then restricted to deer to hold host constant the way the published
PCoA does.

And the qPCR A. phagocytophilum result per sequencing batch. The qPCR plate
is a different batch variable from the sequencing run, so agreement between
them would be evidence the signal is real and disagreement would not settle
anything. I report it as context, not as a test.

NAMES
-----
Column identity goes through tick_lookup on the stem before the first dot,
the pattern build_core_microbiome.py line 51 uses for this file. The part
after the dot is used only as a display label for grouping controls, never
to identify an animal.

OUTPUT
------
A summary to stdout and one table at
results/analysis/anaplasma_batch_check.tsv. I create no directories.
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
OUT_PATH = "results/analysis/anaplasma_batch_check.tsv"

GENUS = "Anaplasma"
READ_FLOOR = 50
FRACTION_FLOOR = 0.01


def genus_of(lineage):
    """I take the last rank of a SILVA lineage string, whatever its depth."""
    return str(lineage).split(";")[-1].strip()


def resolve_library(table, column):
    """I resolve one taxonomy column to an animal, as build_core_microbiome does."""
    stem = str(column).split(".")[0]
    for marker in ("16S", "specimen"):
        hit = tick_lookup.resolve(table, marker, stem)
        if hit is not None:
            return hit
    return None


def batch_label(column):
    """I take the display label after the first dot. This never identifies anything."""
    pieces = str(column).split(".", 1)
    return pieces[1] if len(pieces) == 2 else "(none)"


def is_deer(host):
    """I match both host_bin spellings.

    The crosswalk writes this field as "Deer" for 333 animals and "WTD Deer"
    for 58, so an equality test against "Deer" silently drops a seventh of
    the cohort. The cohort builders in build_core_microbiome.py match on a
    substring for the same reason.
    """
    lowered = str(host).strip().lower()
    return "deer" in lowered or "wtd" in lowered


def count_of(value):
    """I read a positive-subsample count as a number, or None if it is blank."""
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if text.lower() in ("", "nan", "none"):
        return None
    try:
        return int(float(text))
    except ValueError:
        sys.exit("cannot read '" + text + "' as a count in " + QPCR_PATH)


def main():
    if not os.path.isdir("results"):
        sys.exit("run me from the project root; there is no results/ here.")
    for path in (COUNTS_PATH, CROSSWALK_PATH, QPCR_PATH):
        if not os.path.isfile(path):
            sys.exit("missing input: " + path)

    table = tick_lookup.load()
    batches = tick_lookup.load_batches()

    counts = pd.read_csv(COUNTS_PATH, sep="\t", dtype=str)
    first = counts.columns[0]
    samples = [c for c in counts.columns if c != first]
    for column in samples:
        counts[column] = pd.to_numeric(counts[column], errors="coerce").fillna(0)
    counts["genus"] = counts[first].map(genus_of)

    hits = counts[counts["genus"] == GENUS]
    if hits.empty:
        sys.exit("no '" + GENUS + "' row in " + COUNTS_PATH)
    per_column = hits[samples].sum(axis=0)
    totals = counts[samples].sum(axis=0)

    resolved, controls = {}, []
    for column in samples:
        specimen = resolve_library(table, column)
        if specimen is None:
            controls.append(column)
        else:
            resolved[column] = specimen

    crosswalk = pd.read_csv(CROSSWALK_PATH, sep="\t", dtype=str)
    host_of, state_of = {}, {}
    for sample_id, host, state in zip(crosswalk["sample_id"],
                                      crosswalk["host_bin"],
                                      crosswalk["state"]):
        specimen = tick_lookup.resolve(table, "specimen", sample_id)
        if specimen is None:
            continue
        host_of[specimen] = "" if pd.isna(host) else str(host).strip()
        state_of[specimen] = "" if pd.isna(state) else str(state).strip()

    qpcr = pd.read_excel(QPCR_PATH, sheet_name=QPCR_SHEET, dtype=str)
    if "a-phagocytophilum" not in qpcr.columns:
        sys.exit("column 'a-phagocytophilum' missing from " + QPCR_PATH
                 + "\ncolumns present: " + ", ".join(map(str, qpcr.columns)))
    qpcr_of = {}
    for host_id, wells in zip(qpcr["host_id"], qpcr["a-phagocytophilum"]):
        specimen = tick_lookup.resolve(table, "specimen", host_id)
        if specimen is not None:
            qpcr_of[specimen] = count_of(wells)

    print("=== control columns, %d, every one named ===" % len(controls))
    print("%-46s %10s %10s %8s" % ("column", "all reads", GENUS, "fraction"))
    control_rows = []
    for column in sorted(controls, key=lambda c: -per_column[c]):
        total = float(totals[column])
        anap = float(per_column[column])
        share = (anap / total) if total else 0.0
        print("%-46s %10d %10d %7.2f%%" % (column, total, anap, 100.0 * share))
        control_rows.append({"column": column, "kind": "control",
                             "batch": batch_label(column), "specimen": "",
                             "host_bin": "", "state": "",
                             "all_reads": int(total), "anaplasma_reads": int(anap),
                             "anaplasma_fraction": round(share, 6),
                             "qpcr_a_phago_wells": ""})
    print()
    print("control totals: %d reads, of which %s %d (%.2f%%)"
          % (int(totals[controls].sum()), GENUS, int(per_column[controls].sum()),
             100.0 * per_column[controls].sum() / max(totals[controls].sum(), 1)))

    rows = list(control_rows)
    for column, specimen in resolved.items():
        total = float(totals[column])
        anap = float(per_column[column])
        stem = str(column).split(".")[0]
        batch = (tick_lookup.resolve_batch(batches, "16S", stem)
                 or tick_lookup.resolve_batch(batches, "specimen", stem)
                 or batch_label(column))
        wells = qpcr_of.get(specimen)
        rows.append({"column": column, "kind": "specimen", "batch": batch,
                     "specimen": specimen, "host_bin": host_of.get(specimen, ""),
                     "state": state_of.get(specimen, ""),
                     "all_reads": int(total), "anaplasma_reads": int(anap),
                     "anaplasma_fraction": round((anap / total) if total else 0.0, 6),
                     "qpcr_a_phago_wells": "" if wells is None else wells})

    frame = pd.DataFrame(rows)
    frame.to_csv(OUT_PATH, sep="\t", index=False)

    spec = frame[frame["kind"] == "specimen"].copy()
    spec["called"] = ((spec["anaplasma_reads"] >= READ_FLOOR)
                      & (spec["anaplasma_fraction"] >= FRACTION_FLOOR))
    spec["qpcr_pos"] = pd.to_numeric(spec["qpcr_a_phago_wells"],
                                     errors="coerce").fillna(0) > 0

    def summarise(sub, title):
        print()
        print("=== %s: %d libraries ===" % (title, len(sub)))
        print("%-48s %5s %8s %9s %9s %9s"
              % ("batch", "libs", "16S +", "median %", "mean %", "qPCR +"))
        grouped = sub.groupby("batch")
        for batch, part in grouped:
            print("%-48s %5d %8d %8.2f%% %8.2f%% %8d"
                  % (str(batch)[:48], len(part), int(part["called"].sum()),
                     100.0 * part["anaplasma_fraction"].median(),
                     100.0 * part["anaplasma_fraction"].mean(),
                     int(part["qpcr_pos"].sum())))

    summarise(spec, "all specimen libraries by sequencing batch")
    summarise(spec[spec["host_bin"].map(is_deer)],
              "deer only, host held constant")

    deer = spec[spec["host_bin"].map(is_deer)]
    print()
    print("=== deer, 16S Anaplasma calls by state ===")
    print(pd.crosstab(deer["state"].replace("", "(none)"), deer["called"]).to_string())
    print()
    print("=== deer, qPCR A. phagocytophilum by state ===")
    print(pd.crosstab(deer["state"].replace("", "(none)"), deer["qpcr_pos"]).to_string())
    print()
    print("wrote " + OUT_PATH)


if __name__ == "__main__":
    main()
