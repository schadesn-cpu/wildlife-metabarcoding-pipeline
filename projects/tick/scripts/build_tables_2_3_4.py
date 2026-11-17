#!/usr/bin/env python3
"""build_tables_2_3_4.py

I build Tables 2, 3 and 4 from analysis files that already exist. No new
numbers are computed here; I reshape and round.

TABLE 2, qPCR DETECTION BY HOST
-------------------------------
Single ticks and pools, not the "all qPCR rows" block. That block counts rows
rather than animals and includes records excluded from both cohorts, so its
columns do not add up, which in a table reads as an arithmetic error rather
than as the deliberate thing it is.

TABLES 3 AND 4, CORE MICROBIOME
-------------------------------
The twenty most abundant genera, against the figure's twelve, because showing
more than the figure is what a table is for. Each genus carries its enrichment
in the negative controls from the control comparison, so a reader can see why
two rows are drawn hollow in Figure 3 without going to Table S3. A genus the
controls could not speak to, because fewer than five reads were expected,
shows a dash rather than a ratio.

WHAT I READ
-----------
    results/analysis/qpcr_by_host.tsv
    results/analysis/core_microbiome_ixodes_deer.tsv
    results/analysis/core_microbiome_dermacentor_moose.tsv
    results/analysis/control_enrichment.tsv

WHAT I WRITE, WITH --write
--------------------------
    results/analysis/table2_qpcr_by_host.tsv
    results/analysis/table3_core_microbiome_deer.tsv
    results/analysis/table4_core_microbiome_moose.tsv

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/build_tables_2_3_4.py
    python3 scripts/build_tables_2_3_4.py --write
"""
import os
import sys

import pandas as pd

QPCR = "results/analysis/qpcr_by_host.tsv"
DEER = "results/analysis/core_microbiome_ixodes_deer.tsv"
MOOSE = "results/analysis/core_microbiome_dermacentor_moose.tsv"
CONTROLS = "results/analysis/control_enrichment.tsv"

SINGLE = "single ticks, one tick per assay"
POOLS = "pools, one slurry per assay"
TARGETS = [
    ("a-phagocytophilum", "A. phagocytophilum"),
    ("b-microti", "B. microti"),
    ("b-burgdorferi", "B. burgdorferi"),
    ("b-miyamotoi", "B. miyamotoi"),
    ("powassan-virus", "Powassan virus"),
    ("b-mayonii", "B. mayonii"),
    ("e-muris", "E. muris"),
]
TOP = 20
FLOOR = 5.0


def read(path, index=False):
    if not os.path.isfile(path):
        sys.exit("missing %s. Run me from the project root." % path)
    return pd.read_csv(path, sep="\t", index_col=0 if index else None,
                       dtype=None if index else str)


def abundance(value):
    """I round an abundance so 84%, 3.31% and 0.235% all read cleanly."""
    number = float(value)
    if number >= 10:
        return "%.1f" % number
    if number >= 1:
        return "%.2f" % number
    return "%.3f" % number


def show(frame, gap=2):
    """I print a frame as aligned columns, first column left, rest right."""
    columns = list(frame.columns)
    widths = {c: max(len(str(c)), *(len(str(v)) for v in frame[c]))
              for c in columns}
    pad = " " * gap
    header = pad.join(str(c).ljust(widths[c]) if i == 0
                      else str(c).rjust(widths[c])
                      for i, c in enumerate(columns))
    print(header)
    print("-" * len(header))
    for _, row in frame.iterrows():
        print(pad.join(str(row[c]).ljust(widths[c]) if i == 0
                       else str(row[c]).rjust(widths[c])
                       for i, c in enumerate(columns)))
    print()


def table2():
    """I return Table 2, qPCR detection by host in each cohort."""
    source = read(QPCR)
    for column in ("cohort", "host", "n", "ticks"):
        if column not in source.columns:
            sys.exit("%s has no '%s' column" % (QPCR, column))
    rows = []
    for cohort, label in ((SINGLE, "Single ticks"), (POOLS, "Pools")):
        part = source[source["cohort"] == cohort]
        if part.empty:
            sys.exit("%s holds no rows for cohort %r. Its cohorts are %s"
                     % (QPCR, cohort, "; ".join(sorted(set(source["cohort"])))))
        for _, row in part.iterrows():
            entry = {"Cohort": label,
                     # Table 1 calls this group Unknown; matching it here
                     # keeps two adjacent tables from naming one thing twice.
                     "Host": ("Unknown" if str(row["host"]).strip().lower()
                              in ("unresolved", "", "nan")
                              else row["host"]),
                     "Animals": row["n"],
                     "Ticks": row["ticks"]}
            for key, name in TARGETS:
                if key not in source.columns:
                    sys.exit("%s has no column %r" % (QPCR, key))
                count = int(float(row[key]))
                total = int(float(row["n"]))
                entry[name] = ("%d (%.0f%%)" % (count, 100.0 * count / total)
                               if total else "0")
            rows.append(entry)
    return pd.DataFrame(rows)


def enrichment():
    """I return {genus: printable ratio}, or an empty map if unavailable."""
    if not os.path.isfile(CONTROLS):
        print("no %s, so the control column will be blank. Run "
              "check_control_enrichment.py --write first." % CONTROLS)
        return {}
    frame = pd.read_csv(CONTROLS, sep="\t", index_col=0)
    for column in ("expected", "ratio"):
        if column not in frame.columns:
            sys.exit("%s has no '%s' column" % (CONTROLS, column))
    found = {}
    for genus, row in frame.iterrows():
        if float(row["expected"]) < FLOOR:
            found[str(genus)] = "-"
        elif float(row["ratio"]) == 0:
            found[str(genus)] = "0"
        elif float(row["ratio"]) >= 1:
            found[str(genus)] = "%.1f" % float(row["ratio"])
        else:
            found[str(genus)] = "%.2f" % float(row["ratio"])
    return found


def core_table(path, ratios):
    """I return one core microbiome table and the genus count behind it."""
    frame = read(path, index=True)
    for column in ("animals", "animals_with", "prevalence_pct", "mean_pct"):
        if column not in frame.columns:
            sys.exit("%s has no '%s' column; it has %s"
                     % (path, column, ", ".join(map(str, frame.columns))))
    frame = frame[frame["mean_pct"] > 0]
    total = int(frame["animals"].iloc[0])
    shown = frame.nlargest(TOP, "mean_pct")
    rows = []
    for genus, row in shown.iterrows():
        rows.append({
            "Genus": str(genus),
            "Detected in": "%d of %d (%.0f%%)"
                           % (int(row["animals_with"]), total,
                              float(row["prevalence_pct"])),
            "Mean abundance (%)": abundance(row["mean_pct"]),
            "Controls": ratios.get(str(genus), "-"),
        })
    return pd.DataFrame(rows), len(frame), total


def main():
    write = "--write" in sys.argv[1:]
    ratios = enrichment()

    two = table2()
    print("=" * 78)
    print("TABLE 2. qPCR DETECTION BY HOST, SINGLE TICKS AND POOLS")
    print("=" * 78)
    show(two)

    three, deer_genera, deer_animals = core_table(DEER, ratios)
    print("=" * 78)
    print("TABLE 3. CORE MICROBIOME OF IXODES FROM WHITE-TAILED DEER")
    print("=" * 78)
    show(three)
    print("  %d genera detected across %d animals; the %d most abundant are "
          "shown." % (deer_genera, deer_animals, len(three)))
    print()

    four, moose_genera, moose_animals = core_table(MOOSE, ratios)
    print("=" * 78)
    print("TABLE 4. CORE MICROBIOME OF DERMACENTOR FROM MOOSE")
    print("=" * 78)
    show(four)
    print("  %d genera detected across %d animals; the %d most abundant are "
          "shown." % (moose_genera, moose_animals, len(four)))
    print()

    print("  Controls column: the genus's share of the negative controls")
    print("  divided by its share of the whole sequencing pool. Near 1 means")
    print("  the reads in the blanks arrived by index hopping. Well above 1")
    print("  means reagent contamination. A dash means fewer than %d reads"
          % FLOOR)
    print("  were expected in the blanks, so they cannot speak to it.")
    print()

    if not write:
        print("Nothing written. Rerun with --write to save the three tables.")
        return
    for frame, path in ((two, "results/analysis/table2_qpcr_by_host.tsv"),
                        (three,
                         "results/analysis/table3_core_microbiome_deer.tsv"),
                        (four,
                         "results/analysis/table4_core_microbiome_moose.tsv")):
        if os.path.exists(path):
            sys.exit("%s already exists. Move it aside first." % path)
        frame.to_csv(path, sep="\t", index=False)
        print("wrote %s" % path)


if __name__ == "__main__":
    main()
