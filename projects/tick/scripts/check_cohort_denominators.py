"""check_cohort_denominators.py

I report every animal-level denominator the manuscript states, each beside the
rule that sets it, so the cohort table can be filled from one run.

WHAT THE FIRST VERSION LEARNED
------------------------------
The crosswalk holds sample_id, barcode, nums, host_spp, host_bin, state,
num_ticks, is_single_tick and the four name_<marker> columns. It does NOT hold
subsample_count. So the qPCR cohort rule stated in Section 2.6, subsample_count
equal to 1, cannot be evaluated against the crosswalk alone: the rule spans two
files. Rather than guess which one, I search every table under results/ for a
column whose name contains "subsample" and report where it lives, so the next
run can finish the job.

WHAT I COMPUTE FROM THE CROSSWALK
---------------------------------
Everything that does not need subsample_count: the animal total, the tubes
holding exactly one tick, the single-tick flag, the set sequenced as single
ticks, the host breakdown behind Table 1, and the two blank-field counts
Section 2.1 quotes (35 with no tick count, 19 with no host).

THE TWO QUESTIONS THIS IS FOR
-----------------------------
First, Section 2.1 says 93 samples were sequenced as single ticks, 89 passed the
qPCR single-tick rule, and 11 one-tick tubes failed it. 93 minus 11 is 82, not
89, so the 11 come from a set the paragraph never names. I expect that set to be
the tubes holding exactly one tick, which should print as 100.

Second, Table 1 counts 559 animals with a qPCR result while the two cohorts hold
386. The remaining 173 appear nowhere in the manuscript. Closing that needs
subsample_count, which is why I go looking for it.

Run from the project root:
    python3 scripts/check_cohort_denominators.py
"""
import csv
import os
import sys
from collections import Counter

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CROSSWALK = os.path.join(PROJECT_ROOT, "results", "metadata", "crosswalk_4marker.tsv")
RESULTS = os.path.join(PROJECT_ROOT, "results")
MARKERS = ("name_16S", "name_18S", "name_cytb", "name_ITS")


def load_crosswalk(path):
    """I return (rows, fieldnames) for the crosswalk, exiting with a message if it is unusable."""
    if not os.path.exists(path):
        sys.exit(
            "I cannot find the crosswalk at %s.\n"
            "Run this from the project root, or correct CROSSWALK at the top of this file."
            % path
        )
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [row for row in reader]
        fieldnames = list(reader.fieldnames or [])
    if not rows:
        sys.exit("The crosswalk at %s has a header but no rows." % path)
    required = ["sample_id", "host_bin", "num_ticks", "is_single_tick"]
    absent = [name for name in required if name not in fieldnames]
    if absent:
        sys.exit(
            "The crosswalk is missing these columns, which every rule below needs: %s\n"
            "Columns present: %s" % (", ".join(absent), ", ".join(fieldnames))
        )
    return rows, fieldnames


def truthy(value):
    """I read the several ways this project has written a boolean into a TSV."""
    return str(value).strip().lower() in ("1", "true", "t", "yes", "y")


def integer(value):
    """I return an int, or None where the cell is blank or not a number."""
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def filled(row, key):
    """I say whether a cell carries anything, treating the usual null spellings as empty."""
    return str(row.get(key, "")).strip().lower() not in ("", "na", "nan", "none", "null")


def find_subsample_column(root):
    """I return [(path, column, delimiter)] for every table under root with a subsample column."""
    hits = []
    for base, _dirs, files in os.walk(root):
        if os.sep + "_archive" in base:
            continue
        for name in sorted(files):
            if not name.endswith((".tsv", ".csv", ".txt")):
                continue
            path = os.path.join(base, name)
            try:
                with open(path, newline="", encoding="utf-8", errors="replace") as handle:
                    header = handle.readline()
            except OSError as problem:
                print("  could not read %s: %s" % (path, problem))
                continue
            delimiter = "\t" if "\t" in header else ","
            for column in header.rstrip("\n").split(delimiter):
                if "subsample" in column.strip().lower():
                    hits.append((os.path.relpath(path, root), column.strip(), delimiter))
    return hits


def main():
    rows, fieldnames = load_crosswalk(CROSSWALK)
    print("crosswalk: %s" % CROSSWALK)
    print("rows: %d" % len(rows))
    print()

    print("FROM THE CROSSWALK")
    animals = set(row["sample_id"].strip() for row in rows if row["sample_id"].strip())
    print("  %-5d animals, distinct sample_id" % len(animals))

    one_tick = [r for r in rows if integer(r["num_ticks"]) == 1]
    flagged = [r for r in rows if truthy(r["is_single_tick"])]
    both = [r for r in one_tick if truthy(r["is_single_tick"])]
    print("  %-5d tubes with num_ticks equal to 1" % len(one_tick))
    print("  %-5d rows with is_single_tick set" % len(flagged))
    print("  %-5d with both, which should be the 89 of Section 3.2" % len(both))
    print("  %-5d one-tick tubes without the flag, which should be the 11" % (len(one_tick) - len(both)))

    print()
    print("  sequenced as a single tick, by marker:")
    for marker in MARKERS:
        if marker in fieldnames:
            n = sum(1 for r in flagged if filled(r, marker))
            print("    %-12s %-4d  (Section 2.1 says 93 for metabarcoding)" % (marker, n))
        else:
            print("    %-12s absent from the crosswalk" % marker)
    any_marker = [r for r in flagged
                  if any(filled(r, m) for m in MARKERS if m in fieldnames)]
    print("    %-12s %-4d  any marker" % ("union", len(any_marker)))

    print()
    print("  blank fields Section 2.1 quotes:")
    print("    %-5d rows with no num_ticks, stated as 35" % sum(1 for r in rows if integer(r["num_ticks"]) is None))
    print("    %-5d rows with no host_spp, stated as 19" % sum(1 for r in rows if not filled(r, "host_spp")))

    print()
    print("  host_bin, which should reproduce Table 1:")
    for host, count in Counter(r["host_bin"].strip() or "(blank)" for r in rows).most_common():
        print("    %-24s %d" % (host, count))

    print()
    print("WHERE subsample_count LIVES")
    hits = find_subsample_column(RESULTS)
    if not hits:
        print("  No table under results/ has a column whose name contains 'subsample'.")
        print("  The qPCR cohort rule in Section 2.6 therefore cites a field that is not")
        print("  in the results tree. It may be in the metadata workbook instead, which")
        print("  would mean the rule cannot be reproduced from results/ alone.")
    else:
        for path, column, delimiter in hits:
            print("  %-58s column %s (%s)" % (path, column, "tsv" if delimiter == "\t" else "csv"))
        print()
        print("  Point the next run at whichever of these is the qPCR record, and the")
        print("  173 residual can be broken down by subsample_count.")

    print()
    print("STILL NOT ANSWERED HERE, AND WHERE EACH ONE LIVES")
    print("  343, 335, 52   scripts/check_endosymbiont_suppression.py")
    print("  360, 333       the blood-meal step, results/cytb/")
    print("  413            Section 2.6 already carries a VERIFY on this pair")


if __name__ == "__main__":
    main()
