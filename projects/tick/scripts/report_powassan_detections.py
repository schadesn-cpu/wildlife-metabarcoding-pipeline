"""report_powassan_detections.py

I report every Powassan virus detection in every table that carries a Powassan
column, with the collection host, state, tick count and single-tick flag read
from the crosswalk, so a rare-agent detection can be reported specifically.

WHY VERSION TWO
---------------
Version one picked the first candidate table and crashed on a successful
crosswalk match, because "index.get(a) or index.get(b)" asks pandas for the
truth value of a Series. Fixed with an explicit None check.

It also found 9 detections in results/analysis/feeding_ecology_table.tsv while
the manuscript reports 5, one single tick and four deer pools. I do not know
yet whether the tables disagree or whether four positives sit outside both
reported cohorts. So I now report every table that carries a Powassan column
rather than only the first, and I count the single-tick and pooled subsets of
each so the comparison with the manuscript is on the screen.

WHAT I DO NOT DO
----------------
I do not infer host or identity from an identifier string. Host, state and tick
count come from the crosswalk.

I do not decide which table is authoritative, and I do not decide whether a
detection is confirmed. I print each table's own value so the encoding and any
disagreement are visible.
"""
import argparse
import os
import sys

try:
    import pandas
except ImportError:
    sys.exit("I need pandas to read these tables. Run this with the project "
             "environment active, or install pandas first.")

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, ".."))
CROSSWALK = os.path.join(PROJECT_ROOT, "results", "metadata", "crosswalk_4marker.tsv")

SEARCH_DIRS = ["metadata", "data", "docs", "results/metadata", "results/analysis", "."]
TABLE_EXTENSIONS = (".tsv", ".csv", ".txt", ".xlsx", ".xls")
POWASSAN_HINTS = ("powassan", "powv", "pow_v", "dtv", "deer_tick_virus")
SPECIMEN_HINTS = ("sample_id", "specimen", "sample", "barcode", "tube")
COHORT_HINTS = ("subsample_count", "cohort", "is_single_tick", "num_ticks")
NEGATIVE_LABELS = {"", "neg", "negative", "nd", "n/a", "na", "none", "0", "0.0",
                   "undetermined", "undet", "false", "no"}


def load_table(path, header_only=False):
    """I read a tsv, csv or xlsx without guessing the separator."""
    rows = 0 if header_only else None
    extension = os.path.splitext(path)[1].lower()
    if extension in (".xlsx", ".xls"):
        return pandas.read_excel(path, nrows=rows, dtype=str)
    separator = "," if extension == ".csv" else "\t"
    return pandas.read_csv(path, sep=separator, nrows=rows, dtype=str)


def column_matching(columns, hints):
    """I return the first column whose lowercased name contains one of the hints."""
    for hint in hints:
        for column in columns:
            if hint in str(column).strip().lower():
                return column
    return None


def columns_matching(columns, hints):
    """I return every column whose lowercased name contains one of the hints."""
    found = []
    for column in columns:
        lowered = str(column).strip().lower()
        for hint in hints:
            if hint in lowered and column not in found:
                found.append(column)
    return found


def candidate_tables():
    """I collect every table one level deep under the search directories."""
    found = []
    for relative in SEARCH_DIRS:
        directory = os.path.join(PROJECT_ROOT, relative)
        if not os.path.isdir(directory):
            continue
        for entry in sorted(os.listdir(directory)):
            full = os.path.join(directory, entry)
            if os.path.isfile(full) and entry.lower().endswith(TABLE_EXTENSIONS):
                if full not in found:
                    found.append(full)
    return found


def is_positive(value):
    """I treat anything that is not an explicit negative or blank as a detection."""
    if value is None:
        return False
    text = str(value).strip().lower()
    if text in NEGATIVE_LABELS or text.startswith("nan"):
        return False
    return True


def build_index(crosswalk):
    """I key the crosswalk by every identifier spelling I can derive from it."""
    keys = [c for c in ("sample_id", "barcode", "nums") if c in crosswalk.columns]
    if not keys:
        sys.exit("The crosswalk has none of sample_id, barcode or nums. Columns "
                 "are: %s" % ", ".join(str(c) for c in crosswalk.columns))
    index = {}
    for _, row in crosswalk.iterrows():
        for key in keys:
            value = str(row.get(key, "")).strip()
            if not value or value.lower() == "nan":
                continue
            for spelling in (value, value.replace(" ", "-"), value.replace("-", " "),
                             value.replace(" ", "_"), value.replace("_", " ")):
                if spelling not in index:
                    index[spelling] = row
    return index


def resolve(index, specimen):
    """I look up one specimen without ever testing a Series for truthiness."""
    for spelling in (specimen, specimen.replace(" ", "-"), specimen.replace("-", " "),
                     specimen.replace(" ", "_"), specimen.replace("_", " ")):
        row = index.get(spelling)
        if row is not None:
            return row
    return None


def report_table(path, index):
    """I print every Powassan detection in one table and return its count."""
    relative = os.path.relpath(path, PROJECT_ROOT)
    try:
        table = load_table(path)
    except (ValueError, OSError, pandas.errors.ParserError) as problem:
        print("=== %s ===" % relative)
        print("I could not read it: %s" % problem)
        return None

    powassan_column = column_matching(table.columns, POWASSAN_HINTS)
    if powassan_column is None:
        return None

    print("=== %s ===" % relative)
    print("rows: %d   columns (%d): %s"
          % (len(table), len(table.columns),
             ", ".join(str(c) for c in table.columns)))

    specimen_column = column_matching(table.columns, SPECIMEN_HINTS)
    if specimen_column is None:
        print("I cannot tell which column holds the specimen identifier, so I "
              "cannot resolve these rows.")
        return None

    values = table[powassan_column].fillna("")
    distinct = sorted(set(str(v).strip() for v in values))
    print("powassan column: %s   specimen column: %s" % (powassan_column, specimen_column))
    print("distinct values: %s" % ", ".join(distinct[:12]))

    own_cohort = columns_matching(table.columns, COHORT_HINTS)
    positives = table[[is_positive(v) for v in values]]
    print("detections: %d" % len(positives))
    if positives.empty:
        print()
        return 0

    header = ["specimen", "value", "host_spp", "host_bin", "state",
              "num_ticks", "is_single_tick", "name_18S"]
    header.extend(str(c) for c in own_cohort)
    lines = ["\t".join(header)]
    unmatched = []
    singles = pools = neither = 0
    for _, row in positives.iterrows():
        specimen = str(row.get(specimen_column, "")).strip()
        matched = resolve(index, specimen)
        if matched is None:
            unmatched.append(specimen)
            continue
        flag = str(matched.get("is_single_tick", "")).strip()
        count = str(matched.get("num_ticks", "")).strip()
        if flag.lower() in ("true", "1", "yes"):
            singles = singles + 1
        elif count.isdigit() and int(count) > 1:
            pools = pools + 1
        else:
            neither = neither + 1
        fields = [specimen,
                  str(row.get(powassan_column, "")).strip(),
                  str(matched.get("host_spp", "")).strip(),
                  str(matched.get("host_bin", "")).strip(),
                  str(matched.get("state", "")).strip(),
                  count,
                  flag,
                  str(matched.get("name_18S", "")).strip()]
        fields.extend(str(row.get(c, "")).strip() for c in own_cohort)
        lines.append("\t".join(fields))

    print()
    print("\n".join(lines))
    print()
    print("crosswalk single-tick flag true: %d    num_ticks > 1: %d    neither: %d"
          % (singles, pools, neither))
    if unmatched:
        print("not found in the crosswalk: %s" % ", ".join(unmatched))
    print()
    return len(positives)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", default=None,
                        help="report only this one table instead of searching")
    arguments = parser.parse_args()

    if not os.path.isfile(CROSSWALK):
        sys.exit("I cannot find the crosswalk at %s. Run me from the project "
                 "root, or fix the path." % CROSSWALK)
    index = build_index(load_table(CROSSWALK))

    if arguments.table:
        if not os.path.isfile(arguments.table):
            sys.exit("That table does not exist: %s" % arguments.table)
        paths = [os.path.abspath(arguments.table)]
    else:
        paths = candidate_tables()

    counts = {}
    for path in paths:
        result = report_table(path, index)
        if result is not None:
            counts[os.path.relpath(path, PROJECT_ROOT)] = result

    if not counts:
        sys.exit("No table under %s carries a Powassan column. Pass --table "
                 "with the path to the screening record." % ", ".join(SEARCH_DIRS))

    print("SUMMARY, detections per table")
    for relative in sorted(counts):
        print("  %-55s %d" % (relative, counts[relative]))
    print()
    print("The manuscript reports 5: one single tick and four deer pools. Any "
          "table above that disagrees needs explaining before the number goes "
          "in the abstract. Host, state and tick count came from the crosswalk, "
          "not from the identifier string.")


if __name__ == "__main__":
    main()
