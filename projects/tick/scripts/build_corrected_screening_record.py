#!/usr/bin/env python3
"""build_corrected_screening_record.py

I write the screening record again, with the detection counts recomputed
from the well-level run files under the criteria of Section 2.5.

WHY
---
The qPCR Results sheet carries one row per animal and, in each target
column, the number of subsamples the testing laboratory recorded as
positive. Those recorded calls do not consistently implement the rule the
laboratory states, so Section 2.5 counts detections from the run files
instead: the laboratory's two quality criteria as stated, amp_score above
1.0 and cq_conf above 0.6, with the cycle-threshold boundary left as the
supplied calls implement it.

Rather than teach that rule to each of the eight scripts that count
detections, I write one corrected record in the sheet's own shape.
check_qpcr_by_host.py reads this instead of the sheet, and Table 2,
Figure 2 and everything downstream inherit the counts with no change to any
cohort or host logic. The sheet itself is never modified: it stays as the
record of what the laboratory supplied, and the correction is a visible
step in the chain.

WHAT CHANGES
------------
Twenty-one wells of 6,664. Twenty recorded inconclusive clear both quality
criteria and are counted as detections; one recorded positive clears
neither and is not, which is the study's only B. mayonii call.

WHAT I WILL NOT DO
------------------
I will not alter a row I cannot match. Every host_id goes through
tick_lookup; a row whose identifier does not resolve, or that has no wells
in the run file, keeps its original counts and is named in the summary. A
silent repair is worse than a visible gap.

If the sheet carries a totals row, meaning one row whose subsample count
equals the sum of all the others, I recompute it from the corrected columns
and say so. Leaving a stale total in place would be its own quiet error.

OUTPUT
------
results/analysis/qpcr_screening_record_corrected.tsv, same columns and same
row order as the sheet. I create no directories.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import tick_lookup

WELLS_PATH = "metadata/USGS-qpcr-combined-results.xlsx"
WELLS_SHEET = "in"
WORKBOOK = "metadata/tick_full_metadata.xlsx"
QPCR_SHEET = "qPCR Results"
OUT_PATH = "results/analysis/qpcr_screening_record_corrected.tsv"

AMP_CUTOFF = 1.0
CQ_CUTOFF = 0.6

# Well-file target name -> screening-record column name.
TARGETS = (
    ("A. phagocytophilum", "a-phagocytophilum"),
    ("B. microti", "b-microti"),
    ("B. burgdorferi", "b-burgdorferi"),
    ("B. mayonii", "b-mayonii"),
    ("B. miyamotoi", "b-miyamotoi"),
    ("E. muris", "e-muris"),
    ("Powassan virus", "powassan-virus"),
)


def narrow_rule(recorded, amp, cq):
    """I apply the two quality criteria to the positive boundary only.

    A cycle-threshold cutoff of 35 taken literally would promote 25 wells
    between cycles 34 and 35 that the supplied calls record as negative, and
    40 of the 61 places where the stated rule and the supplied calls differ
    sit in that one-cycle band, so the implemented boundary was stricter
    than the round number. I leave every well recorded negative alone.
    """
    clears = ((not pd.isna(amp) and amp > AMP_CUTOFF)
              and (not pd.isna(cq) and cq > CQ_CUTOFF))
    if recorded == "Inconclusive" and clears:
        return "Positive"
    if recorded == "Positive" and not clears:
        return "Inconclusive"
    return recorded


def count_of(value):
    """I read a count, or None if the cell is blank."""
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if text.lower() in ("", "nan", "none"):
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def main():
    if not os.path.isdir("results"):
        sys.exit("run me from the project root; there is no results/ here.")
    for path in (WELLS_PATH, WORKBOOK):
        if not os.path.isfile(path):
            sys.exit("missing input: " + path
                     + ("\nCopy USGS-qpcr-combined-results.xlsx into metadata/."
                        if path == WELLS_PATH else ""))

    table = tick_lookup.load()

    wells = pd.read_excel(WELLS_PATH, sheet_name=WELLS_SHEET, dtype=str)
    for column in ("amp_score", "cq_conf"):
        wells[column] = pd.to_numeric(wells[column], errors="coerce")
    wells["call"] = [narrow_rule(p, a, c) for p, a, c
                     in zip(wells["pass_check"], wells["amp_score"],
                            wells["cq_conf"])]
    wells["animal"] = [tick_lookup.resolve(table, "specimen", h)
                       for h in wells["host_id"]]
    moved = int((wells["call"] != wells["pass_check"]).sum())
    print("wells %d, of which the criteria of Section 2.5 move %d"
          % (len(wells), moved))

    positives = {}
    for (well_name, column) in TARGETS:
        part = wells[(wells["target_name"] == well_name)
                     & (wells["call"] == "Positive")
                     & wells["animal"].notna()]
        for animal, n in part.groupby("animal").size().items():
            positives[(animal, column)] = int(n)
    have_wells = set(wells.loc[wells["animal"].notna(), "animal"])

    record = pd.read_excel(WORKBOOK, sheet_name=QPCR_SHEET, dtype=str)
    missing = [c for _, c in TARGETS if c not in record.columns]
    if missing or "host_id" not in record.columns:
        sys.exit("columns missing from " + WORKBOOK + " sheet '" + QPCR_SHEET
                 + "': " + ", ".join(missing + (["host_id"] if "host_id"
                                                not in record.columns else []))
                 + "\ncolumns present: " + ", ".join(map(str, record.columns)))

    # find a totals row: one whose subsample count equals the sum of the rest
    totals_index = None
    if "subsample_count" in record.columns:
        counts = [count_of(v) for v in record["subsample_count"]]
        numbers = [c for c in counts if c is not None]
        for i, c in enumerate(counts):
            if c is not None and c == sum(numbers) - c:
                totals_index = i
                break

    changed, untouched = [], []
    out = record.copy()
    for i, row in record.iterrows():
        if i == totals_index:
            continue
        animal = tick_lookup.resolve(table, "specimen", row["host_id"])
        if animal is None or animal not in have_wells:
            untouched.append(str(row["host_id"]))
            continue
        for _, column in TARGETS:
            before = count_of(row[column])
            after = positives.get((animal, column), 0)
            if before != after:
                changed.append((str(row["host_id"]), column, before, after))
            out.at[i, column] = str(after)

    if totals_index is not None:
        for _, column in TARGETS:
            total = sum(int(out.at[i, column]) for i in out.index
                        if i != totals_index and str(out.at[i, column]).isdigit())
            out.at[totals_index, column] = str(total)
        print("recomputed the totals row at index %d from the corrected columns"
              % totals_index)

    out.to_csv(OUT_PATH, sep="\t", index=False)

    print()
    print("rows %d, of which left untouched %d" % (len(record), len(untouched)))
    if untouched:
        print("   no wells or no resolved identifier: %s"
              % ", ".join(sorted(set(untouched))))
    print()
    print("cells changed: %d" % len(changed))
    for host_id, column, before, after in sorted(changed, key=lambda r: (r[1], r[0])):
        print("   %-24s %-20s %s -> %d"
              % (host_id, column, "blank" if before is None else before, after))
    print()
    print("wrote " + OUT_PATH)
    print("Next: point check_qpcr_by_host.py at this file, then rerun it,")
    print("build_tables_2_3_4.py and build_fig_qpcr_by_host.py in that order.")


if __name__ == "__main__":
    main()
