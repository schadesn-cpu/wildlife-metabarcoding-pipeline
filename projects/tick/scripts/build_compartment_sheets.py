#!/usr/bin/env python3
"""I build one workbook of the sample compartments, side by side.

WHY
---
The collection holds more than one kind of thing in a tube. Ixodes from deer
and Dermacentor from moose are the two tick and host pairs the paper reports.
Twelve tubes held deer keds and no ticks at all, recorded as KEDS ONLY in the
extraction bench sheets for batch 07, and six tubes held spiders alongside
their ticks. Those four groups have never been put next to each other, and
the keds group is the one that matters: if Bartonella dominates libraries
built from keds with no ticks in them, the ked attribution in the taxon
decision table is settled by observation rather than by association.

HOW I STAY COMPARABLE
---------------------
I do not reimplement the genus calling. I import build_core_microbiome and
use its read_counts, one_per_animal and summarise, so every group here is
built on the same MIN_READS floor, the same dropped lineages and the same
renormalisation as Table 3. The taxonomy and host-lineage steps sit inside
that script's main() rather than in a function, so those lines are copied
verbatim and are marked below.

Because they are copied, they can drift. So I rebuild the deer and moose
tables and compare them against the saved core_microbiome files. If they do
not reproduce, the copy has drifted from its source and nothing I say about
the other groups can be trusted; I report the difference and stop.

MEMBERSHIP
----------
The ked and spider groups come from the accessioning and bench records, not
from sequence. They are listed below as identifiers and every one goes
through tick_lookup. A member that does not resolve, or that has no library
above the read floor, is named and dropped rather than quietly missing.

WHAT I DO NOT CLAIM
-------------------
A tube recorded KEDS ONLY is a tube in which no ticks were found when it was
opened. The accessioning record for those same twelve animals gives tick
counts of one to three, so the two records disagree and I report both. The
spiders are not a compartment: every spider tube also held ticks, so that
sheet is a co-detection check rather than a group.

OUTPUT
------
results/analysis/compartment_sheets.xlsx, one sheet per group plus a README
sheet and a membership sheet. I create no directories and I refuse to
overwrite. Without --write I print everything and save nothing.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import tick_lookup
import build_core_microbiome as core
import check_batch_arms as arms

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
CORRECTED = "results/analysis/qpcr_screening_record_corrected.tsv"
SAVED_DEER = "results/analysis/core_microbiome_ixodes_deer.tsv"
SAVED_MOOSE = "results/analysis/core_microbiome_dermacentor_moose.tsv"
OUT_PATH = "results/analysis/compartment_sheets.xlsx"

# Recorded KEDS ONLY in the extraction bench sheet for batch 07, meaning keds
# were present and no ticks were found when the tube was opened.
KEDS_ONLY = ["USGSTICKS-4101", "USGSTICKS-4103", "USGSTICKS-4104",
             "USGSTICKS-4107", "USGSTICKS-4109", "USGSTICKS-4111",
             "USGSTICKS-4114", "USGSTICKS-4115", "USGSTICKS-4117",
             "USGSTICKS-6303", "USGSTICKS-6306", "USGSTICKS-6312"]
# Spiders recorded in the accessioning notes. Every one of these tubes also
# held ticks, so this is a co-detection check and not a compartment.
SPIDERS = ["EQ02441292", "EQ02441293", "EQ02441301",
           "EQ02442187", "EQ02442201", "EQ02442208"]

TARGETS = ["a-phagocytophilum", "b-burgdorferi", "b-microti", "b-miyamotoi",
           "powassan-virus", "b-mayonii", "e-muris"]


def resolve_all(table, names, label):
    """I turn a list of recorded identifiers into animals, naming any miss."""
    out, missed = [], []
    for name in names:
        hit = tick_lookup.resolve(table, "specimen", name)
        if hit is None:
            missed.append(name)
        else:
            out.append(str(hit))
    print("%s: %d identifiers, resolved %d" % (label, len(names), len(out)))
    if missed:
        print("   did NOT resolve: %s" % ", ".join(missed))
        print("   I drop these rather than guess which animal they are.")
    return out


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


def qpcr_sheets(table):
    """I return per-animal qPCR calls split into the two cohorts.

    I apply the cohort rule of Section 2.5 directly rather than importing it,
    because I read the crosswalk myself: a single tick is one tested
    subsample with the is_single_tick flag set, and a pool is one tested
    subsample of more than one tick.
    """
    for path in (CROSSWALK, CORRECTED):
        if not os.path.isfile(path):
            sys.exit("missing %s" % path)
    cross = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    single, ticks, host, state = {}, {}, {}, {}
    for _, row in cross.iterrows():
        animal = tick_lookup.resolve(table, "specimen", row["sample_id"])
        if animal is None:
            continue
        animal = str(animal)
        single[animal] = str(row.get("is_single_tick")).strip().lower() in (
            "true", "1", "yes")
        ticks[animal] = count_of(row.get("num_ticks"))
        host[animal] = "" if pd.isna(row.get("host_bin")) else str(row["host_bin"])
        state[animal] = "" if pd.isna(row.get("state")) else str(row["state"])

    record = pd.read_csv(CORRECTED, sep="\t", dtype=str)
    missing = [c for c in TARGETS if c not in record.columns]
    if missing:
        sys.exit("%s is missing %s" % (CORRECTED, missing))
    rows, seen, controls, repeats = [], set(), [], []
    for _, r in record.iterrows():
        animal = tick_lookup.resolve(table, "specimen", r.get("host_id"))
        if animal is None:
            continue
        animal = str(animal)
        if animal in seen:
            repeats.append(animal)
            continue
        seen.add(animal)
        detected = sum(1 for t in TARGETS if (count_of(r.get(t)) or 0) > 0)
        if detected >= 5:
            controls.append(animal)
            continue
        subs = count_of(r.get("subsample_count"))
        if subs != 1:
            cohort = "outside both cohorts"
        elif single.get(animal):
            cohort = "single tick"
        elif (ticks.get(animal) or 0) > 1:
            cohort = "pool"
        else:
            cohort = "outside both cohorts"
        entry = {"animal": animal, "cohort": cohort, "host": host.get(animal, ""),
                 "state": state.get(animal, ""), "num_ticks": ticks.get(animal),
                 "subsample_count": subs}
        total = 0
        for t in TARGETS:
            n = count_of(r.get(t)) or 0
            entry[t] = n
            total += 1 if n > 0 else 0
        entry["agents_detected"] = total
        rows.append(entry)
    if controls:
        print("   excluded %d row(s) positive for five or more targets, which "
              "is a positive control and not a tick: %s"
              % (len(controls), ", ".join(sorted(controls))))
    if repeats:
        print("   %d repeated animal row(s), keeping the first of each: %s"
              % (len(repeats), ", ".join(sorted(set(repeats)))))
    frame = pd.DataFrame(rows)
    singles = frame[frame["cohort"] == "single tick"].copy()
    pools = frame[frame["cohort"] == "pool"].copy()
    print("\nqPCR cohorts from the corrected record: single ticks %d, pools %d"
          % (len(singles), len(pools)))
    if len(singles) != 89 or len(pools) != 297:
        sys.exit("\n   Section 2.5 gives 89 single ticks and 297 pools and I "
                 "got %d and %d.\n   These sheets are meant to sit beside "
                 "Table 2, so a membership that\n   differs from it is not "
                 "something to warn about and continue past.\n   I have "
                 "written nothing." % (len(singles), len(pools)))
    print("   matches the 89 and 297 of Section 2.5.")
    return singles, pools


def main():
    write = "--write" in sys.argv[1:]
    if not os.path.isdir("results/analysis"):
        sys.exit("run me from the project root; there is no results/analysis here.")

    table, module = core.lookup()
    counts = core.read_counts(table, module)
    animals = [core.animal_of(table, module, s) for s in counts.index]
    print("counts: %d libraries x %d features; resolved to an animal %d"
          % (counts.shape[0], counts.shape[1], sum(1 for a in animals if a)))

    # --- copied verbatim from build_core_microbiome.main(), see the docstring
    tax_path, tax_covered, mapping = arms.read_taxonomy(list(counts.columns))
    print("taxonomy: %s (%d of %d features)" % (tax_path, tax_covered, counts.shape[1]))
    host_like = [f for f in counts.columns
                 if any(word in str(mapping.get(f, "")).lower()
                        for word in ("mitochondri", "chloroplast", "eukaryota"))]
    print("dropping %d features on a mitochondrial, chloroplast or eukaryotic "
          "lineage" % len(host_like))
    counts = counts.drop(columns=host_like)
    genus = pd.Series({f: arms.genus_label(mapping.get(f, ""))
                       for f in counts.columns})
    # --- end of copied block

    deer, moose, species = core.group_members(table, module)
    chosen = core.one_per_animal(counts, animals)
    print("libraries above the read floor, one per animal: %d\n" % len(chosen))

    deer_table = core.summarise(counts, genus,
                                chosen[chosen["animal"].isin(deer)],
                                "Ixodes from white-tailed deer")
    moose_table = core.summarise(counts, genus,
                                 chosen[chosen["animal"].isin(moose)],
                                 "Dermacentor from moose")

    # --- the guard: do I reproduce what is already on disk? ----------------
    print("=" * 72)
    print("GUARD: do I reproduce the saved core microbiome tables?")
    print("=" * 72)
    ok = True
    for label, built, path in (("deer", deer_table, SAVED_DEER),
                               ("moose", moose_table, SAVED_MOOSE)):
        if built is None or not os.path.isfile(path):
            print("  %-6s cannot compare (missing %s)" % (label, path)); ok = False; continue
        saved = pd.read_csv(path, sep="\t", index_col=0)
        shared = [c for c in saved.columns if c in built.columns]
        common = saved.index.intersection(built.index)
        same_rows = len(common) == len(saved.index) == len(built.index)
        diff = (saved.loc[common, shared].astype(float)
                - built.loc[common, shared].astype(float)).abs().max().max()
        print("  %-6s saved %d genera, built %d, shared %d, largest difference %.3g"
              % (label, len(saved), len(built), len(common), diff))
        if not same_rows or diff > 1e-6:
            ok = False
    if not ok:
        sys.exit("\nI do not reproduce the saved tables. The block I copied from\n"
                 "build_core_microbiome.main() has drifted from its source, so the\n"
                 "ked and spider sheets would not be comparable to Table 3.\n"
                 "I have written nothing.")
    print("  reproduced exactly; the other groups are built on the same machinery.\n")

    # --- the groups that are new ------------------------------------------
    ked_animals = resolve_all(table, KEDS_ONLY, "KEDS ONLY tubes")
    spider_animals = resolve_all(table, SPIDERS, "tubes with spiders recorded")
    ked_chosen = chosen[chosen["animal"].isin(ked_animals)]
    spider_chosen = chosen[chosen["animal"].isin(spider_animals)]
    for label, names, got in (("ked", ked_animals, ked_chosen),
                              ("spider", spider_animals, spider_chosen)):
        lost = sorted(set(names) - set(got["animal"]))
        if lost:
            print("%s group: %d of %d have no library above the %d-read floor: %s"
                  % (label, len(lost), len(names), core.MIN_READS, ", ".join(lost)))
    print()
    ked_table = core.summarise(counts, genus, ked_chosen,
                               "keds only, no ticks found in the tube")
    spider_table = core.summarise(counts, genus, spider_chosen,
                                  "tubes that also held spiders")

    # --- does the ked group sit inside the deer cohort of Table 3? ---------
    overlap = sorted(set(ked_animals) & set(deer))
    print("=" * 72)
    print("OVERLAP: keds-only animals inside the Table 3 deer cohort")
    print("=" * 72)
    print("  %d of %d" % (len(overlap), len(ked_animals)))
    if overlap:
        print("  %s" % ", ".join(overlap))
        print("  Table 3 therefore includes libraries from tubes in which no")
        print("  ticks were found. That is worth a decision, not a silent fix.")
    print()

    # --- the twelve, animal by animal ----------------------------------
    print("=" * 72)
    print("THE KEDS-ONLY TUBES, ANIMAL BY ANIMAL")
    print("=" * 72)
    part = counts.loc[ked_chosen.index]
    by_genus = part.T.groupby(genus).sum().T
    drop_these = [c for c in by_genus.columns
                  if c in ("Mitochondria", "Chloroplast", "Unassigned")
                  or c.endswith("(domain)")]
    kept = by_genus.drop(columns=drop_these)
    frac = 100 * kept.div(kept.sum(axis=1), axis=0)
    # I repeated what summarise() does internally to get per-animal values,
    # so I check that my column means reproduce the table it returned.
    if ked_table is not None:
        worst = 0.0
        for name in ("Rickettsia", "Bartonella", "Arsenophonus"):
            if name in frac.columns and name in ked_table.index:
                worst = max(worst, abs(frac[name].mean()
                                       - ked_table.loc[name, "mean_pct"]))
        print("  per-animal means reproduce the group table to %.3g points"
              % worst)
        if worst > 1e-6:
            sys.exit("  they do not. My per-animal recomputation is wrong "
                     "and nothing below it can be trusted.")
    called = {}
    for a, ip, dp in zip(species["animal"], species["Ixodes_pos"],
                         species["Dermacentor_pos"]):
        called[str(a)] = (str(ip), str(dp))
    print()
    print("  %-22s %8s %11s %11s %13s %7s %s"
          % ("animal", "depth", "Rickettsia", "Bartonella", "Arsenophonus",
             "in T3?", "Ixodes_pos"))
    for lib, row in ked_chosen.iterrows():
        animal = str(row["animal"])
        ip, dp = called.get(animal, ("", ""))
        print("  %-22s %8.0f %10.1f%% %10.1f%% %12.1f%% %7s %s"
              % (animal, row["depth"],
                 frac.loc[lib].get("Rickettsia", 0.0),
                 frac.loc[lib].get("Bartonella", 0.0),
                 frac.loc[lib].get("Arsenophonus", 0.0),
                 "yes" if animal in deer else "no", ip))
    print()
    print("  A library that is mostly Rickettsia plausibly held a tick the")
    print("  bench note missed. One that is mostly Bartonella did not.")
    print()
    singles, pools = qpcr_sheets(table)

    membership = []
    for animal in sorted(set(list(deer) + list(moose) + ked_animals + spider_animals)):
        membership.append({
            "animal": animal,
            "ixodes_deer_cohort": animal in deer,
            "dermacentor_moose_cohort": animal in moose,
            "keds_only_tube": animal in ked_animals,
            "spider_in_tube": animal in spider_animals,
            "library_above_floor": animal in set(chosen["animal"]),
        })
    membership = pd.DataFrame(membership)

    readme = pd.DataFrame({"sheet": [
        "ixodes_deer", "dermacentor_moose", "keds_only", "spiders",
        "qpcr_single_tick", "qpcr_pooled", "membership"], "what it holds": [
        "Genus table for Ixodes from white-tailed deer. Same rows as "
        "core_microbiome_ixodes_deer.tsv, verified cell by cell.",
        "Genus table for Dermacentor from moose, verified the same way.",
        "Genus table for the 12 tubes recorded KEDS ONLY in batch 07. No "
        "ticks were found in these when opened; the accessioning record "
        "gives them 1 to 3 ticks, and the two records disagree.",
        "Genus table for the 6 tubes with spiders recorded. Every one also "
        "held ticks, so this is a co-detection check, not a compartment.",
        "Per-animal qPCR calls for the single-tick cohort, counted under the "
        "criteria of Section 2.5 from the corrected screening record.",
        "Per-animal qPCR calls for the pooled cohort, same source.",
        "Which animals belong to which group, and whether each has a library "
        "above the %d-read floor." % core.MIN_READS]})

    if not write:
        print("Nothing written. Rerun with --write to save " + OUT_PATH)
        return
    if os.path.exists(OUT_PATH):
        sys.exit("%s already exists. Move it aside first." % OUT_PATH)
    try:
        writer = pd.ExcelWriter(OUT_PATH, engine="openpyxl")
    except (ImportError, ValueError) as problem:
        sys.exit("cannot open a workbook writer (%s). Install openpyxl or ask "
                 "me for TSVs instead." % problem)
    with writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        deer_table.to_excel(writer, sheet_name="ixodes_deer")
        moose_table.to_excel(writer, sheet_name="dermacentor_moose")
        if ked_table is not None:
            ked_table.to_excel(writer, sheet_name="keds_only")
        if spider_table is not None:
            spider_table.to_excel(writer, sheet_name="spiders")
        singles.to_excel(writer, sheet_name="qpcr_single_tick", index=False)
        pools.to_excel(writer, sheet_name="qpcr_pooled", index=False)
        membership.to_excel(writer, sheet_name="membership", index=False)
    print("wrote " + OUT_PATH)


if __name__ == "__main__":
    main()
