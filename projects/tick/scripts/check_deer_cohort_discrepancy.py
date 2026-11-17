#!/usr/bin/env python3
"""check_deer_cohort_discrepancy.py

I find out why two paragraphs of the manuscript disagree about how many deer
carried Ixodes and not Dermacentor, and I count the decontam negative set.
I write nothing.

THE DISAGREEMENT
----------------
Section 3.1 says 424 white-tailed deer had an 18S result, that 371 of them
carried Ixodes, and that none carried Dermacentor. That makes 371 deer with
Ixodes and without Dermacentor, and it reconciles against the 53 uninformative
deer, because 371 plus 53 is 424.

Section 2.6 says 366 deer carried Ixodes and not Dermacentor, and that 330 of
those cleared the read floor. The 330 cohort of Section 3.3 descends from the
smaller number, so if 366 is wrong then every count in Section 3.3 rests on
the wrong denominator.

371 minus 366 is 5. My first guess is that five deer rows fail to resolve to
an animal through tick_lookup and are dropped by the cohort builder but
counted in Section 3.1. I do not assume it. I report the count at each step so
the five are visible wherever they are lost, and if they are lost somewhere
else the step-by-step will say so.

WHY I DO NOT GUESS AT HOSTS
---------------------------
Host comes from the host_bin field, matched on a substring because the
crosswalk writes it as both "Deer" and "WTD Deer". I never infer a host from
an identifier string. Where the two host records disagree for one animal I
report the disagreement rather than choosing between them.

THE SECOND QUESTION
-------------------
The cohorts table needs the size of the decontam negative set, which the
manuscript describes by rule and never counts. decontam_16S.R prints it as
"vetted negatives" when it runs, but it needs /tmp/16s_feature_table.tsv,
which may be gone. I therefore reconstruct its selection rule here and say so
plainly: this is a reconstruction, and the authoritative number is whatever
the R script prints. If the two disagree, trust the R script and tell me.

WHAT I READ
-----------
    results/analysis/tick_species_by_host.tsv
    results/metadata/crosswalk_4marker.tsv
    results/16S/all/taxonomy/taxonomy_counts_L6_16S.tsv   (column names only)
    scripts/tick_lookup.py

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_deer_cohort_discrepancy.py
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import tick_lookup

SPECIES_PATH = "results/analysis/tick_species_by_host.tsv"
CROSSWALK_PATH = "results/metadata/crosswalk_4marker.tsv"
COUNTS_PATH = "results/16S/all/taxonomy/taxonomy_counts_L6_16S.tsv"

TRUE_WORDS = ("1", "1.0", "true", "yes")

# decontam_16S.R selects its negatives by these name patterns, then removes
# four blanks found to hold real template and three Rickettsia-laden NTCs.
CONTROL_WORDS = ("xb", "ntc", "pos", "pac", "pc-well", "blank", "extraction")
POSITIVE_WORDS = ("pos", "pc-well", "pac")
MISLABELED = ("XB1-2", "XBG12", "XB-A1", "XB-D8")
HOT_CONTROLS = ("NTC.250609_Batch-11", "NTC-2.250609_Batch-11",
                "NTC-H11.250327_Batch-08")


def is_deer(host):
    """I match both host_bin spellings, "Deer" and "WTD Deer"."""
    lowered = str(host).strip().lower()
    return "deer" in lowered or "wtd" in lowered


def truthy(value):
    """I read a positive-call flag, however the file spells it."""
    return str(value).strip().lower() in TRUE_WORDS


def animal_of(table, name):
    """I resolve one specimen name to an animal, or return ''.

    The order matches the cohort builders: the specimen marker first, then the
    16S marker on the stem before the first dot.
    """
    text = str(name).strip()
    stem = text.split(".")[0]
    for marker, value in (("specimen", text), ("specimen", stem),
                          ("16S", stem)):
        hit = tick_lookup.resolve(table, marker, value)
        if hit is not None:
            return str(hit)
    return ""


def require(path):
    """I refuse to carry on without an input, and say which one."""
    if not os.path.isfile(path):
        sys.exit("missing input: " + path
                 + "\nRun me from the project root, /home/users/sne24/"
                   "meedlab/tick_project.")


def column_check(frame, path, wanted):
    """I check the columns I need are present and name what is there if not."""
    for column in wanted:
        if column not in frame.columns:
            sys.exit(path + " has no column '" + column + "'\ncolumns "
                     "present: " + ", ".join(map(str, frame.columns)))


def deer_steps(table):
    """I count deer with Ixodes at each step, so a loss of five is visible."""
    species = pd.read_csv(SPECIES_PATH, sep="\t", dtype=str)
    column_check(species, SPECIES_PATH,
                 ("specimen", "host_bin", "Ixodes_pos", "Dermacentor_pos"))

    rows = len(species)
    deer_rows = species[species["host_bin"].map(is_deer)]
    ixodes = deer_rows[deer_rows["Ixodes_pos"].map(truthy)]
    both = ixodes[ixodes["Dermacentor_pos"].map(truthy)]
    ixodes_only = ixodes[~ixodes["Dermacentor_pos"].map(truthy)]

    resolved, unresolved = {}, []
    for specimen in ixodes_only["specimen"]:
        animal = animal_of(table, specimen)
        if animal:
            resolved.setdefault(animal, []).append(specimen)
        else:
            unresolved.append(specimen)

    print("=== deer carrying Ixodes, step by step ===")
    print("%-56s %6s" % ("step", "count"))
    print("%-56s %6d" % ("rows in tick_species_by_host.tsv", rows))
    print("%-56s %6d" % ("  host_bin matches deer or wtd", len(deer_rows)))
    print("%-56s %6d" % ("    of those, Ixodes_pos is true", len(ixodes)))
    print("%-56s %6d" % ("      and Dermacentor_pos also true", len(both)))
    print("%-56s %6d" % ("      and Dermacentor_pos not true", len(ixodes_only)))
    print("%-56s %6d" % ("        of those, tick_lookup resolves", len(ixodes_only) - len(unresolved)))
    print("%-56s %6d" % ("        distinct animals after resolution", len(resolved)))
    print()
    print("Section 3.1 reports 371 and Section 2.6 reports 366.")
    print("The step where the count falls to 366, if it does, is the answer.")
    print()

    if both.empty:
        print("No deer row carries both genera, so 'Ixodes and not")
        print("Dermacentor' and 'Ixodes' are the same set for deer, as")
        print("Section 3.1 states.")
    else:
        print("%d deer rows carry BOTH genera, which Section 3.1 says does"
              % len(both))
        print("not happen. Those rows are:")
        for specimen in both["specimen"]:
            print("    " + str(specimen))
    print()

    if unresolved:
        print("%d deer specimens with Ixodes do not resolve to an animal and"
              % len(unresolved))
        print("are therefore dropped by the cohort builders while still being")
        print("counted in Section 3.1:")
        for specimen in unresolved:
            print("    " + str(specimen))
    else:
        print("Every deer specimen with Ixodes resolves to an animal, so the")
        print("five are lost somewhere other than resolution.")
    print()

    repeats = {a: s for a, s in resolved.items() if len(s) > 1}
    if repeats:
        print("%d animals contribute more than one specimen row, which is a"
              % len(repeats))
        print("second way a row count and an animal count can differ:")
        for animal in sorted(repeats)[:10]:
            print("    %s: %s" % (animal, ", ".join(repeats[animal])))
        if len(repeats) > 10:
            print("    and %d more" % (len(repeats) - 10))
    print()
    return set(resolved)


def host_agreement(table, animals):
    """I compare host_bin in the crosswalk against the one used above."""
    crosswalk = pd.read_csv(CROSSWALK_PATH, sep="\t", dtype=str)
    column_check(crosswalk, CROSSWALK_PATH, ("sample_id", "host_bin"))

    crosswalk_host = {}
    for sample_id, host in zip(crosswalk["sample_id"], crosswalk["host_bin"]):
        animal = animal_of(table, sample_id)
        if animal:
            crosswalk_host[animal] = "" if pd.isna(host) else str(host).strip()

    missing = sorted(a for a in animals if a not in crosswalk_host)
    disagree = sorted(a for a in animals
                      if a in crosswalk_host and not is_deer(crosswalk_host[a]))

    print("=== host record agreement for those animals ===")
    print("animals checked                              %6d" % len(animals))
    print("  absent from crosswalk_4marker.tsv          %6d" % len(missing))
    print("  present but host_bin is not a deer         %6d" % len(disagree))
    for animal in missing[:10]:
        print("    absent: " + animal)
    for animal in disagree[:10]:
        print("    not deer in crosswalk: %s -> %s"
              % (animal, crosswalk_host[animal]))
    print()


def negative_set():
    """I reconstruct decontam_16S.R's negative set and count it."""
    frame = pd.read_csv(COUNTS_PATH, sep="\t", dtype=str, nrows=1)
    columns = [c for c in frame.columns[1:]]

    chosen = []
    for name in columns:
        lowered = str(name).lower()
        if not any(word in lowered for word in CONTROL_WORDS):
            continue
        if any(word in lowered for word in POSITIVE_WORDS):
            continue
        if str(name).split(".")[0] in MISLABELED:
            continue
        if name in HOT_CONTROLS:
            continue
        chosen.append(name)

    print("=== decontam negative set, reconstructed from decontam_16S.R ===")
    print("This is a reconstruction of the R script's name rule applied to")
    print("the columns of %s." % COUNTS_PATH)
    print("The authoritative number is the 'vetted negatives' line that")
    print("decontam_16S.R prints. If they disagree, trust the R script.")
    print()
    print("columns in the table                         %6d" % len(columns))
    print("selected as vetted negatives                 %6d" % len(chosen))
    for name in sorted(chosen):
        print("    " + name)
    print()


def main():
    if not os.path.isdir("results"):
        sys.exit("run me from the project root; there is no results/ here.")
    for path in (SPECIES_PATH, CROSSWALK_PATH, COUNTS_PATH):
        require(path)

    table = tick_lookup.load()
    animals = deer_steps(table)
    host_agreement(table, animals)
    negative_set()


if __name__ == "__main__":
    main()
