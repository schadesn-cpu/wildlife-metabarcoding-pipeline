#!/usr/bin/env python3
"""build_core_microbiome.py

I build the two core microbiome tables the restructured paper needs: Ixodes
from white-tailed deer, and Dermacentor from moose, at genus level, one library
per animal. With --write I save them under results/analysis/.

Genus level, because at ASV level the sequencing run explains a third of the
variation among ordinary deer ticks (Rickettsia reads split between two
near-identical ASVs in run-dependent proportions); collapsed to genus that
effect is 2.9% and not significant.

I read results/decontam/16S/table_decontam_16S.tsv (counts, all hosts), the
verified 16S taxonomy, results/analysis/tick_species_by_host.tsv (host and tick
genus) and the library lookup. A library needs 1,000 reads to enter.

USAGE
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/build_core_microbiome.py [--write]
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_batch_arms as arms  # noqa: E402

COUNTS = "results/decontam/16S/table_decontam_16S.tsv"
SPECIES = "results/analysis/tick_species_by_host.tsv"
OUT_DEER = "results/analysis/core_microbiome_ixodes_deer.tsv"
OUT_MOOSE = "results/analysis/core_microbiome_dermacentor_moose.tsv"
MIN_READS = 1000
SHOW = 15
ALWAYS = ("Rickettsia", "Francisella", "Coxiella", "Anaplasma", "Borreliella",
          "Borrelia", "Bartonella", "Arsenophonus", "Rickettsiella")

# Tubes the extraction bench recorded as holding keds and no ticks, and whose
# composition agrees: Rickettsia under 3% of classified bacterial reads with
# Bartonella and Arsenophonus over 95%. These are deer tubes with no tick
# material, so they are not Ixodes from white-tailed deer. Four further
# KEDS ONLY tubes in this cohort carry Rickettsia at 50.5% or above and are
# deliberately kept. The per-tube figures are in build_compartment_sheets.py.
EXCLUDED_NO_TICK_MATERIAL = ("USGSTICKS-4101", "USGSTICKS-4103",
                             "USGSTICKS-4109", "USGSTICKS-4111",
                             "USGSTICKS-6303")


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s). Run from the project root."
                 % problem)
    return tick_lookup.load(), tick_lookup


def animal_of(table, module, name):
    """I turn a library or specimen name into its animal, or ''."""
    text = str(name).strip()
    stem = text.split(".")[0]
    for marker, value in (("16S", text), ("16S", stem), ("specimen", stem)):
        hit = module.resolve(table, marker, value)
        if hit is not None:
            return str(hit)
    return ""


def read_counts(table, module):
    """I return counts as libraries x features, dropping annotation columns.

    I decide which axis holds the libraries by asking the lookup: the axis
    whose names resolve to animals is the library axis.
    """
    if not os.path.isfile(COUNTS):
        sys.exit("missing %s" % COUNTS)
    frame = pd.read_csv(COUNTS, sep="\t", index_col=0, low_memory=False)
    frame = frame.loc[:, [c for c in frame.columns
                          if pd.to_numeric(frame[c], errors="coerce").notna().all()]]
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0)
    rows = sum(1 for name in frame.index if animal_of(table, module, name))
    columns = sum(1 for name in frame.columns if animal_of(table, module, name))
    if columns > rows:
        frame = frame.T
    return frame


def group_members(table, module):
    """I return the deer-with-Ixodes and moose-with-Dermacentor animal sets."""
    if not os.path.isfile(SPECIES):
        sys.exit("missing %s" % SPECIES)
    species = pd.read_csv(SPECIES, sep="\t", dtype=str)
    needed = ["specimen", "host_bin", "Ixodes_pos", "Dermacentor_pos"]
    missing = [c for c in needed if c not in species.columns]
    if missing:
        sys.exit("%s has no column(s) %s; it has: %s"
                 % (SPECIES, missing, ", ".join(species.columns[:12])))
    species["animal"] = [animal_of(table, module, s) for s in species["specimen"]]
    species = species[species["animal"] != ""]
    host = species["host_bin"].astype(str).str.lower()
    ixodes = species["Ixodes_pos"].astype(str).isin(["1", "1.0", "True", "TRUE", "true"])
    derm = species["Dermacentor_pos"].astype(str).isin(["1", "1.0", "True", "TRUE", "true"])
    deer = set(species.loc[host.str.contains("deer|wtd") & ixodes & ~derm, "animal"])
    moose = set(species.loc[host.str.contains("moose") & derm, "animal"])
    dropped = set()
    for listed in EXCLUDED_NO_TICK_MATERIAL:
        animal = animal_of(table, module, listed)
        if not animal:
            sys.exit("the exclusion list names %s and tick_lookup does not "
                     "resolve it. I will not drop an animal I cannot name."
                     % listed)
        if animal in deer:
            deer.discard(animal)
            dropped.add(animal)
    print("excluded %d of the %d listed libraries holding no tick material%s"
          % (len(dropped), len(EXCLUDED_NO_TICK_MATERIAL),
             (": " + ", ".join(sorted(dropped))) if dropped else ""))
    return deer, moose, species


def one_per_animal(counts, animals):
    """I keep the deepest library per animal, above the read floor."""
    depth = counts.sum(axis=1)
    frame = pd.DataFrame({"animal": animals, "depth": depth.values},
                         index=counts.index)
    frame = frame[(frame["animal"] != "") & (frame["depth"] >= MIN_READS)]
    frame = frame.sort_values("depth", ascending=False)
    return frame.drop_duplicates("animal", keep="first")


def summarise(counts, genus, chosen, label):
    """I print and return the genus table for one group."""
    print("=" * 72)
    print(label.upper())
    print("=" * 72)
    if chosen.empty:
        print("  no libraries qualify.\n")
        return None
    part = counts.loc[chosen.index]
    by_genus = part.T.groupby(genus).sum().T
    # Mitochondria and chloroplast are tick and host reads, not microbiome, and
    # the domain-level and unassigned bins have no genus identity. I drop them
    # and renormalize to reads classified to a genus.
    dropped = [c for c in by_genus.columns
               if c in ("Mitochondria", "Chloroplast", "Unassigned")
               or c.endswith("(domain)")]
    kept = by_genus.drop(columns=dropped)
    print("  dropped %s: %.1f%% of reads" % (", ".join(sorted(dropped)),
          100 * by_genus[dropped].sum().sum() / by_genus.sum().sum()))
    by_genus = kept.div(kept.sum(axis=1), axis=0)
    table = pd.DataFrame({
        "animals_with": (by_genus > 0).sum(axis=0),
        "prevalence_pct": 100 * (by_genus > 0).mean(axis=0),
        "mean_pct": 100 * by_genus.mean(axis=0),
        "median_pct": 100 * by_genus.median(axis=0),
    }).sort_values("mean_pct", ascending=False)
    table.insert(0, "animals", len(chosen))
    print("  %d animals, one library each, median depth %.0f reads"
          % (len(chosen), chosen["depth"].median()))
    print()
    print("  %-26s %8s %10s %9s %9s" % ("genus", "animals", "prevalence",
                                        "mean", "median"))
    shown = list(table.index[:SHOW])
    for name in ALWAYS:
        if name in table.index and name not in shown:
            shown.append(name)
    for name in shown:
        row = table.loc[name]
        print("  %-26s %8d %9.0f%% %8.1f%% %8.1f%%"
              % (name[:26], row["animals_with"], row["prevalence_pct"],
                 row["mean_pct"], row["median_pct"]))
    print()
    return table


def main():
    write = "--write" in sys.argv[1:]
    table, module = lookup()
    counts = read_counts(table, module)
    print("counts: %d libraries x %d features (%s)"
          % (counts.shape[0], counts.shape[1], COUNTS))
    animals = [animal_of(table, module, s) for s in counts.index]
    print("libraries resolved to an animal: %d of %d"
          % (sum(1 for a in animals if a), len(animals)))
    tax_path, tax_covered, mapping = arms.read_taxonomy(list(counts.columns))
    print("taxonomy: %s (%d of %d features)" % (tax_path, tax_covered, counts.shape[1]))
    # Mitochondrial and chloroplast reads are tick, host or plant, not
    # microbiome, and SILVA classifies some of them down to Acari or Arachnida.
    # I drop the whole lineage rather than the genus label alone.
    host_like = [f for f in counts.columns
                 if any(word in str(mapping.get(f, "")).lower()
                        for word in ("mitochondri", "chloroplast", "eukaryota"))]
    print("dropping %d features on a mitochondrial, chloroplast or eukaryotic "
          "lineage" % len(host_like))
    counts = counts.drop(columns=host_like)
    genus = pd.Series({f: arms.genus_label(mapping.get(f, ""))
                       for f in counts.columns})
    deer, moose, species = group_members(table, module)
    print("specimens with host and tick calls: %d; deer with Ixodes: %d; "
          "moose with Dermacentor: %d\n" % (len(species), len(deer), len(moose)))

    chosen = one_per_animal(counts, animals)
    deer_table = summarise(counts, genus, chosen[chosen["animal"].isin(deer)],
                           "Ixodes from white-tailed deer")
    moose_table = summarise(counts, genus, chosen[chosen["animal"].isin(moose)],
                            "Dermacentor from moose")

    if not write:
        print("Nothing written. Rerun with --write to save both tables.")
        return
    for path, result in ((OUT_DEER, deer_table), (OUT_MOOSE, moose_table)):
        if result is None:
            continue
        if os.path.exists(path):
            sys.exit("%s already exists. Move it aside first." % path)
        result.to_csv(path, sep="\t")
        print("wrote %s (%d genera)" % (path, len(result)))


if __name__ == "__main__":
    main()
