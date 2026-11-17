#!/usr/bin/env python3
"""recount_qpcr_by_rule.py

I recompute the qPCR detections twice, once from the recorded pass_check
column and once from the rule the testing laboratory states, and show what
differs.

WHY
---
Seth Goodnight ran these plates and on 2026-10-07 gave the rule:

    Positive      Ct < 35 and amp_score > 1.0 and cq_conf > 0.6
    Inconclusive  Ct < 35 but one or both of the other two failed
    Negative      Ct >= 35

Applied to all 6,664 wells that reproduces 6,603 of the recorded calls. Of
the 61 that differ, 40 sit within one cycle of the Ct 35 cutoff and are a
rounding question. Twenty are wells the rule calls positive that were
recorded inconclusive. One is recorded positive and the rule calls
inconclusive: the study's only B. mayonii detection.

Section 3.2 and Table 2 count the recorded column. Before anyone decides
whether to keep doing that, the numbers under both rules have to be visible
side by side, resolved to animals and cohorts rather than wells.

I CHECK MYSELF FIRST
--------------------
Table 2 was built from the screening record, not from this well file. So
before I report anything about the alternative rule, I recompute Table 2
from the wells using pass_check as recorded and compare it to the screening
record's own counts. If those disagree, the two sources do not describe the
same thing and nothing below can be trusted; I say so and stop.

NAMES
-----
Every host_id goes through tick_lookup. The map from the well file's target
names to the screening record's column names is a fixed table below, not
derived by string manipulation.

OUTPUT
------
A summary to stdout and one table at
results/analysis/qpcr_recount_by_rule.tsv. I create no directories.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import tick_lookup

WELLS_PATH = "metadata/USGS-qpcr-combined-results.xlsx"
WELLS_SHEET = "in"
QPCR_PATH = "metadata/tick_full_metadata.xlsx"
QPCR_SHEET = "qPCR Results"
CROSSWALK_PATH = "results/metadata/crosswalk_4marker.tsv"
OUT_PATH = "results/analysis/qpcr_recount_by_rule.tsv"

CT_CUTOFF = 35.0
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


def stated_rule(ct, amp, cq):
    """I apply the laboratory's stated rule to one well."""
    if pd.isna(ct) or not (ct < CT_CUTOFF):
        return "Negative"
    if (not pd.isna(amp) and amp > AMP_CUTOFF) and (not pd.isna(cq) and cq > CQ_CUTOFF):
        return "Positive"
    return "Inconclusive"


def narrow_rule(recorded, amp, cq):
    """I apply amp_score and cq_conf to the positive boundary only.

    The stated rule's Ct < 35 promotes 25 B. burgdorferi wells recorded
    Negative at Ct 34 to 35, which is the band where 40 of the 61 recorded
    disagreements sit, so the implemented cutoff was evidently stricter than
    the round number. I leave every well recorded Negative alone and apply
    the two quality metrics where they are unambiguous: a well recorded
    inconclusive that clears both becomes positive, and a well recorded
    positive that clears neither stops being one.
    """
    clears = ((amp is not None and not pd.isna(amp) and amp > AMP_CUTOFF)
              and (cq is not None and not pd.isna(cq) and cq > CQ_CUTOFF))
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
        sys.exit("cannot read '" + text + "' as a count")


def main():
    if not os.path.isdir("results"):
        sys.exit("run me from the project root; there is no results/ here.")
    for path in (WELLS_PATH, QPCR_PATH, CROSSWALK_PATH):
        if not os.path.isfile(path):
            sys.exit("missing input: " + path
                     + ("\nCopy USGS-qpcr-combined-results.xlsx into metadata/."
                        if path == WELLS_PATH else ""))

    table = tick_lookup.load()

    wells = pd.read_excel(WELLS_PATH, sheet_name=WELLS_SHEET, dtype=str)
    for column in ("ct", "amp_score", "cq_conf"):
        wells[column] = pd.to_numeric(wells[column], errors="coerce")
    wells["stated"] = [stated_rule(a, b, c) for a, b, c
                       in zip(wells["ct"], wells["amp_score"], wells["cq_conf"])]
    wells["narrow"] = [narrow_rule(p, b, c) for p, b, c
                       in zip(wells["pass_check"], wells["amp_score"],
                              wells["cq_conf"])]
    wells["animal"] = [tick_lookup.resolve(table, "specimen", h)
                       for h in wells["host_id"]]
    unresolved = wells[wells["animal"].isna()]
    print("wells %d   resolving to an animal %d   not resolving %d"
          % (len(wells), int(wells["animal"].notna().sum()), len(unresolved)))
    if len(unresolved):
        print("   host_id values that did not resolve: %s"
              % ", ".join(sorted(set(unresolved["host_id"].astype(str)))[:8]))

    # screening record: cohort membership and the counts Table 2 was built on
    record = pd.read_excel(QPCR_PATH, sheet_name=QPCR_SHEET, dtype=str)
    subsamples, recorded_counts = {}, {}
    for _, row in record.iterrows():
        animal = tick_lookup.resolve(table, "specimen", row["host_id"])
        if animal is None:
            continue
        subsamples[animal] = count_of(row["subsample_count"])
        for _, column in TARGETS:
            if column in record.columns:
                recorded_counts[(animal, column)] = count_of(row[column]) or 0

    crosswalk = pd.read_csv(CROSSWALK_PATH, sep="\t", dtype=str)
    host_of, spp_of, state_of, ticks_of, single_of = {}, {}, {}, {}, {}
    for _, row in crosswalk.iterrows():
        animal = tick_lookup.resolve(table, "specimen", row["sample_id"])
        if animal is None:
            continue
        host_of[animal] = "" if pd.isna(row["host_bin"]) else str(row["host_bin"]).strip()
        spp_of[animal] = "" if pd.isna(row["host_spp"]) else str(row["host_spp"]).strip()
        state_of[animal] = "" if pd.isna(row["state"]) else str(row["state"]).strip()
        ticks_of[animal] = count_of(row["num_ticks"])
        single_of[animal] = str(row["is_single_tick"]).strip().lower() in ("true", "1", "yes")

    def cohort(animal):
        """I apply the cohort rules of Section 2.5 to one animal."""
        n = subsamples.get(animal)
        if n != 1:
            return "outside both cohorts"
        if single_of.get(animal):
            return "Single ticks"
        if (ticks_of.get(animal) or 0) > 1:
            return "Pools"
        return "outside both cohorts"

    # --- self-check: does the well file reproduce the screening record? -----
    print()
    print("=== self-check: wells under the recorded call against the screening record ===")
    mismatch = 0
    checked = 0
    for (well_name, column) in TARGETS:
        part = wells[(wells["target_name"] == well_name) & wells["animal"].notna()]
        tally = part[part["pass_check"] == "Positive"].groupby("animal").size().to_dict()
        animals = set(tally) | set(a for (a, c) in recorded_counts if c == column)
        for animal in animals:
            checked += 1
            if tally.get(animal, 0) != recorded_counts.get((animal, column), 0):
                mismatch += 1
                if mismatch <= 5:
                    print("   %-22s %-20s wells say %d, record says %d"
                          % (animal, column, tally.get(animal, 0),
                             recorded_counts.get((animal, column), 0)))
    print("   animal-by-target comparisons %d, disagreements %d" % (checked, mismatch))
    if mismatch:
        print()
        print("   The well file and the screening record do not agree. Everything")
        print("   below compares two readings of the well file only, and the")
        print("   relationship to the published Table 2 is not established.")

    # --- the detection tables, both ways -----------------------------------
    rows = []
    for well_name, column in TARGETS:
        part = wells[(wells["target_name"] == well_name) & wells["animal"].notna()]
        for animal, group in part.groupby("animal"):
            rows.append({
                "animal": animal,
                "target": column,
                "cohort": cohort(animal),
                "host_bin": host_of.get(animal, ""),
                "host_spp": spp_of.get(animal, ""),
                "state": state_of.get(animal, ""),
                "wells": len(group),
                "recorded_positive": int((group["pass_check"] == "Positive").sum()),
                "stated_positive": int((group["stated"] == "Positive").sum()),
                "narrow_positive": int((group["narrow"] == "Positive").sum()),
            })
    frame = pd.DataFrame(rows)
    frame["detected_recorded"] = frame["recorded_positive"] > 0
    frame["detected_stated"] = frame["stated_positive"] > 0
    frame["detected_narrow"] = frame["narrow_positive"] > 0
    frame.to_csv(OUT_PATH, sep="\t", index=False)

    print()
    print("=== how many animal-target pairs each rule moves, against the record ===")
    for label, flag in (("stated rule, verbatim", "detected_stated"),
                        ("narrow rule, positive boundary only", "detected_narrow")):
        moved = frame[frame["detected_recorded"] != frame[flag]]
        print("   %-38s gained %3d  lost %3d"
              % (label, int(moved[flag].sum()), int((~moved[flag]).sum())))

    changed = frame[frame["detected_recorded"] != frame["detected_stated"]]
    print()
    print("=== animals whose detection status changes under the stated rule ===")
    print("%d animal-target pairs change, across %d animals"
          % (len(changed), changed["animal"].nunique()))
    if len(changed):
        print()
        gained = changed[changed["detected_stated"]]
        lost = changed[~changed["detected_stated"]]
        print("gained %d, lost %d" % (len(gained), len(lost)))
        print()
        print("%-22s %-20s %-22s %-16s %-6s" % ("animal", "target", "cohort", "host", "state"))
        for _, r in changed.sort_values(["target", "cohort", "animal"]).iterrows():
            print("%-22s %-20s %-22s %-16s %-6s  %s"
                  % (r["animal"], r["target"], r["cohort"],
                     (r["host_spp"] or r["host_bin"] or "-")[:16], r["state"] or "-",
                     "GAINED" if r["detected_stated"] else "LOST"))

    print()
    sub = frame[frame["cohort"].isin(("Single ticks", "Pools"))]
    totals = sub.drop_duplicates(["cohort", "host_bin", "animal"]) \
                .groupby(["cohort", "host_bin"]).size().rename("animals")
    columns = [c for _, c in TARGETS]
    for label, flag in (("recorded pass_check", "detected_recorded"),
                        ("laboratory's stated rule, verbatim", "detected_stated"),
                        ("narrow rule, positive boundary only", "detected_narrow")):
        print("=== Table 2 under %s: animals with a detection ===" % label)
        pivot = sub[sub[flag]].pivot_table(index=["cohort", "host_bin"],
                                           columns="target", values="animal",
                                           aggfunc="count", fill_value=0)
        # every cohort and host row, including those with no detection at all
        pivot = pivot.reindex(index=totals.index, columns=columns, fill_value=0)
        print(pivot.fillna(0).astype(int).join(totals).to_string())
        print()
    print("wrote " + OUT_PATH)


if __name__ == "__main__":
    main()
