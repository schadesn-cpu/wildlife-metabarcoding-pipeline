#!/usr/bin/env python3
"""check_threshold_independence.py

I rerun the threshold-independence check for the tick-genus result on the full
library join, replacing the 2026-08-10 check that was run before
build_species_figure.py was put on the shared lookup.

WHY THIS NEEDS RERUNNING
------------------------
The caption claims: "across five inclusion rules (>=1, >=10, >=100 reads;
>=1%, >=10% relative) moose stayed 11/11 and deer 0-1/222-273". That was true
of the join as it stood on 2026-08-10, which resolved 498 of 907 libraries and
covered 399 of 563 specimens. The join now resolves 853 and covers all 563, and
the deer denominator went from 273 to 424. A claim about how a result behaves
across thresholds has to be made on the data the result is computed from.

HOW I AVOID BEING A SEVENTEENTH RESOLVER
----------------------------------------
I do not join anything. build_species_figure.py already writes
results/analysis/tick_species_by_host.tsv with one row per informative
specimen, carrying per-genus read counts AND per-genus fractions. Every
inclusion rule in the claim is a threshold on one of those two columns, so I
read that file and apply the rules to it. The join, the summing across a
specimen's libraries, and the host assignment all happen exactly once, in
build_species_figure.py, and I inherit them.

That also means this script is only as current as the last run of
build_species_figure.py. I check the file's timestamp against the script's and
say so if it is older.

WHAT COUNTS AS THE RESULT SURVIVING
-----------------------------------
Two claims, and I report them separately rather than as one verdict:

    1. every informative moose specimen is Dermacentor-positive
    2. deer are at or near zero Dermacentor

Claim 2 is the one the diversity restriction rests on. At the loosest absolute
rule a single deer has historically appeared (ALHT-Deer-64, 15 Dermacentor
reads against 53,102 Ixodes = 0.03%, inside the blank range), so "0 or 1" is
the honest form of the claim, not "0".

USAGE
-----
    cd ~/meedlab/tick_project
    python3 scripts/build_species_figure.py     # if it has not been run since the join fix
    python3 check_threshold_independence.py
"""
import os
import sys

import pandas as pd

TABLE = "results/analysis/tick_species_by_host.tsv"
SOURCE_SCRIPT = "scripts/build_species_figure.py"
OUT = "results/analysis/threshold_independence.tsv"
REPORT = "results/analysis/threshold_independence_report.txt"

BLANK_MAX = 28          # highest Dermacentor count in a blank or control
GENERA = ["Ixodes", "Dermacentor"]

# The five rules named in the caption, plus the production rule for reference.
RULES = [
    ("absolute", ">= 1 read", "reads", 1),
    ("absolute", ">= 10 reads", "reads", 10),
    ("absolute", ">= 100 reads", "reads", 100),
    ("relative", ">= 1% of tick reads", "frac", 0.01),
    ("relative", ">= 10% of tick reads", "frac", 0.10),
]

report_lines = []


def say(message=""):
    print(message)
    report_lines.append(str(message))


def main():
    if not os.path.isfile(TABLE):
        sys.exit("not found: " + TABLE
                 + "\nRun build_species_figure.py first; it writes this file.")
    if os.path.isfile(SOURCE_SCRIPT):
        if os.path.getmtime(TABLE) < os.path.getmtime(SOURCE_SCRIPT):
            say("WARNING: " + TABLE + " is older than " + SOURCE_SCRIPT + ".")
            say("         Rerun build_species_figure.py before trusting this.")
            say("")

    frame = pd.read_csv(TABLE, sep="\t")
    needed = ["host_bin", "total_tick_reads"]
    for genus in GENERA:
        needed += [genus + "_reads", genus + "_frac"]
    missing = [c for c in needed if c not in frame.columns]
    if missing:
        sys.exit("column(s) missing from " + TABLE + ": " + ", ".join(missing)
                 + "\nThis script expects the per-genus reads and frac columns "
                   "that build_species_figure.py writes.")

    say("threshold independence, rerun " + pd.Timestamp.today().strftime("%Y-%m-%d"))
    say("source: " + TABLE + "  (" + str(len(frame))
        + " informative specimens)")
    say("")
    say("Informative means the specimen produced ANY tick DNA. That is not a")
    say("threshold and does not move between rules, so the denominators below")
    say("are constant and only the numerators respond to the rule.")
    say("")

    moose = frame[frame["host_bin"] == "Moose"]
    deer = frame[frame["host_bin"] == "Deer"]
    say("denominators: moose " + str(len(moose)) + ", deer " + str(len(deer))
        + ", all hosts " + str(len(frame)))
    say("")

    rows = []
    header = ("%-22s %-12s %10s %10s %14s %14s"
              % ("rule", "kind", "Ixodes", "Dermac.", "moose Derm", "deer Derm"))
    say(header)
    say("-" * len(header))

    for kind, label, column, cutoff in RULES:
        def positive(sub, genus):
            return sub[genus + "_" + column].fillna(0) >= cutoff

        ix = int(positive(frame, "Ixodes").sum())
        de = int(positive(frame, "Dermacentor").sum())
        moose_de = int(positive(moose, "Dermacentor").sum())
        deer_de = int(positive(deer, "Dermacentor").sum())
        say("%-22s %-12s %10d %10d %14s %14s"
            % (label, kind, ix, de,
               "%d/%d" % (moose_de, len(moose)),
               "%d/%d" % (deer_de, len(deer))))
        rows.append({"rule": label, "kind": kind, "column": column,
                     "cutoff": cutoff, "Ixodes_pos": ix, "Dermacentor_pos": de,
                     "moose_Dermacentor": moose_de, "moose_n": len(moose),
                     "deer_Dermacentor": deer_de, "deer_n": len(deer)})

    # the production rule, both floors at once
    prod = ((frame["Dermacentor_reads"].fillna(0) >= 50)
            & (frame["Dermacentor_frac"].fillna(0) >= 0.01))
    prod_ix = ((frame["Ixodes_reads"].fillna(0) >= 50)
               & (frame["Ixodes_frac"].fillna(0) >= 0.01))
    moose_prod = prod.loc[moose.index]
    deer_prod = prod.loc[deer.index]
    say("%-22s %-12s %10d %10d %14s %14s"
        % (">= 50 AND >= 1%", "production", int(prod_ix.sum()), int(prod.sum()),
           "%d/%d" % (int(moose_prod.sum()), len(moose)),
           "%d/%d" % (int(deer_prod.sum()), len(deer))))
    rows.append({"rule": ">= 50 AND >= 1%", "kind": "production",
                 "column": "both", "cutoff": "50 / 0.01",
                 "Ixodes_pos": int(prod_ix.sum()),
                 "Dermacentor_pos": int(prod.sum()),
                 "moose_Dermacentor": int(moose_prod.sum()), "moose_n": len(moose),
                 "deer_Dermacentor": int(deer_prod.sum()), "deer_n": len(deer)})

    say("")
    results = pd.DataFrame(rows)

    moose_all = bool((results["moose_Dermacentor"] == results["moose_n"]).all())
    deer_max = int(results["deer_Dermacentor"].max())

    say("=== the two claims, separately ===")
    if moose_all:
        say("  MOOSE: every informative moose is Dermacentor-positive under every")
        say("         rule, including >= 10%% relative. %d/%d throughout."
            % (len(moose), len(moose)))
    else:
        worst = results.loc[results["moose_Dermacentor"].idxmin()]
        say("  MOOSE: NOT constant. Lowest is %d/%d under '%s'."
            % (worst["moose_Dermacentor"], worst["moose_n"], worst["rule"]))
        say("         The 11/11 claim does not survive on the full join.")

    if deer_max == 0:
        say("  DEER:  ZERO Dermacentor under every rule, including >= 1 read,")
        say("         across %d deer specimens." % len(deer))
    elif deer_max <= 1:
        loose = results[results["deer_Dermacentor"] > 0]["rule"].tolist()
        say("  DEER:  at most 1 of %d, appearing only under: %s"
            % (len(deer), ", ".join(loose)))
        say("         Consistent with index bleed; blanks carried up to %d "
            "Dermacentor reads." % BLANK_MAX)
    else:
        say("  DEER:  up to %d of %d Dermacentor-positive. This is MORE than the"
            % (deer_max, len(deer)))
        say("         0-or-1 the caption claims. Investigate before quoting the")
        say("         zero-deer result or the diversity restriction that rests")
        say("         on it.")

    say("")
    if deer_max > 0:
        say("the deer specimen(s) involved, at the loosest rule:")
        loose_deer = deer[deer["Dermacentor_reads"].fillna(0) >= 1]
        cols = [c for c in ("specimen", "state", "total_tick_reads",
                            "Ixodes_reads", "Dermacentor_reads",
                            "Dermacentor_frac") if c in loose_deer.columns]
        for _, row in loose_deer[cols].iterrows():
            say("    " + "  ".join(str(row[c]) for c in cols))
        say("")

    results.to_csv(OUT, sep="\t", index=False)
    say("wrote " + OUT)
    with open(REPORT, "w") as handle:
        handle.write("\n".join(report_lines))
    print("report -> " + REPORT)


if __name__ == "__main__":
    main()
