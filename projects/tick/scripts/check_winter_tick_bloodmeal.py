#!/usr/bin/env python3
"""check_winter_tick_bloodmeal.py

I report the cytb blood-meal reads of every Dermacentor (winter tick) specimen,
and the same totals for Ixodes specimens from deer. I write nothing.

Cervids are reported at family level. The classifier emits Alces and Odocoileus
labels, but BLAST verification does not support them: ASVs labelled alces
(n = 868) top-match Rangifer tarandus at 96%, and a label of sp. matches
Odocoileus at 100%. build_single_tick.py says the same. The moose attribution
rests on the collection record, not on cytb.

I read qiime2/cytb/all/dada2/table.qza, the verified cytb taxonomy under
results/cytb/all/taxonomy_verified/, results/analysis/tick_species_by_host.tsv
and the library lookup.

USAGE
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_winter_tick_bloodmeal.py
"""
import glob
import os
import shutil
import sys
import tempfile
import zipfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

TABLE = "qiime2/cytb/all/dada2/table.qza"
TAXONOMY_GLOB = "results/cytb/all/taxonomy_verified/*.tsv"
SPECIES = "results/analysis/tick_species_by_host.tsv"
MIN_READS = 1
GROUPS = [
    ("Cervidae (deer or moose)", ("alces", "odocoileus", "cervidae", "rangifer",
                                  "cervus", "capreolus")),
    ("Deer ked", ("lipoptena", "hippobosc")),
    ("Bird", ("aves", "geothlypis", "turdus", "carduelis", "spinus", "junco",
              "melospiza", "meleagris", "corvus", "setophaga", "vermivora",
              "parulidae", "passer", "sturnus")),
    ("Bear", ("ursus",)),
    ("Porcupine", ("erethizon",)),
    ("Rabbit", ("sylvilagus", "lepus")),
    ("Opossum", ("didelph",)),
    ("Raccoon", ("procyon",)),
    ("Human", ("homo",)),
    ("Other mammal", ("canis", "felis", "bos", "ovis", "sus", "mus", "rattus",
                      "peromyscus", "tamias", "sciurus", "vulpes", "mephitis",
                      "marmota", "castor", "mustela", "martes")),
]


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
    for marker, value in (("cytb", text), ("cytb", stem), ("16S", stem),
                          ("specimen", stem)):
        hit = module.resolve(table, marker, value)
        if hit is not None:
            return str(hit)
    return ""


def read_table():
    """I return the cytb feature table as libraries x ASVs."""
    if not os.path.isfile(TABLE):
        sys.exit("missing %s" % TABLE)
    try:
        import biom
    except ImportError:
        sys.exit("biom is not importable. Run inside the qiime2 environment.")
    workdir = tempfile.mkdtemp(prefix="cytb_")
    try:
        with zipfile.ZipFile(TABLE) as archive:
            wanted = [n for n in archive.namelist()
                      if n.endswith("data/feature-table.biom")]
            if not wanted:
                sys.exit("no feature-table.biom inside %s" % TABLE)
            archive.extract(wanted[0], workdir)
            table = biom.load_table(os.path.join(workdir, wanted[0]))
            return pd.DataFrame(
                table.matrix_data.toarray().T,
                index=[str(s) for s in table.ids(axis="sample")],
                columns=[str(f) for f in table.ids(axis="observation")])
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def read_taxonomy(features):
    """I return the verified cytb taxonomy that covers the most ASVs."""
    best = None
    for path in sorted(glob.glob(TAXONOMY_GLOB)):
        try:
            frame = pd.read_csv(path, sep="\t", dtype=str)
        except (ValueError, OSError, pd.errors.ParserError):
            continue
        tax_col = next((c for c in frame.columns
                        if c.lower() in ("taxon", "taxonomy", "verified_taxon",
                                         "taxon_lineage")), None)
        if tax_col is None:
            continue
        mapping = dict(zip(frame[frame.columns[0]], frame[tax_col]))
        covered = sum(1 for f in features if f in mapping)
        if best is None or covered > best[1]:
            best = (path, covered, mapping)
    if best is None or best[1] == 0:
        sys.exit("no cytb taxonomy in %s matches the feature table."
                 % TAXONOMY_GLOB)
    return best


def group_of(taxon):
    """I return the host group for one taxon string."""
    text = str(taxon).lower()
    for name, words in GROUPS:
        if any(word in text for word in words):
            return name
    return "Unassigned or other"


def main():
    table, module = lookup()
    counts = read_table()
    path, covered, mapping = read_taxonomy(list(counts.columns))
    print("cytb table: %d libraries x %d ASVs" % counts.shape)
    print("taxonomy:   %s (%d of %d ASVs)" % (path, covered, counts.shape[1]))
    groups = pd.Series({f: group_of(mapping.get(f, "")) for f in counts.columns})
    by_group = counts.T.groupby(groups).sum().T

    species = pd.read_csv(SPECIES, sep="\t", dtype=str)
    species["animal"] = [animal_of(table, module, s) for s in species["specimen"]]
    positive = ("1", "1.0", "True", "TRUE", "true")
    derm = set(species.loc[species["Dermacentor_pos"].isin(positive), "animal"])
    ixodes = set(species.loc[species["Ixodes_pos"].isin(positive), "animal"])
    derm.discard("")
    ixodes.discard("")

    animals = pd.Series([animal_of(table, module, s) for s in counts.index],
                        index=counts.index)
    print("cytb libraries matched to an animal: %d of %d"
          % (int((animals != "").sum()), len(animals)))
    print()

    print("=" * 72)
    print("WINTER TICK (Dermacentor) SPECIMENS: cytb READS BY HOST")
    print("=" * 72)
    columns = [c for c in by_group.columns if c != "Unassigned or other"]
    rows = []
    for animal in sorted(derm):
        libraries = animals.index[animals == animal]
        if not len(libraries):
            rows.append((animal, 0, "no cytb library", ""))
            continue
        totals = by_group.loc[libraries].sum(axis=0)
        named = totals[columns]
        named = named[named >= MIN_READS].sort_values(ascending=False)
        rows.append((animal, int(totals.sum()),
                     "; ".join("%s %d" % (k, v) for k, v in named.items())
                     or "no host reads",
                     "%d" % int(totals.get("Unassigned or other", 0))))
    print("  %-22s %9s %-46s %s" % ("specimen", "cytb reads", "host reads",
                                    "other"))
    for animal, total, named, other in rows:
        print("  %-22s %9d %-46s %s" % (animal[:22], total, named[:46], other))

    with_reads = [r for r in rows if r[1] > 0]
    cervid_only = [r for r in with_reads
                   if r[2].startswith("Cervidae") and ";" not in r[2]]
    print()
    dominant = 0
    for animal in sorted(derm):
        libraries = animals.index[animals == animal]
        if not len(libraries):
            continue
        named = by_group.loc[libraries].sum(axis=0)[columns]
        if named.sum() > 0 and named.get("Cervidae (deer or moose)", 0) >= 0.99 * named.sum():
            dominant += 1
    print("  %d Dermacentor specimens, %d with cytb reads, %d with cervid reads"
          " and nothing else," % (len(rows), len(with_reads), len(cervid_only)))
    print("  and %d where cervid reads are at least 99%% of the host reads."
          % dominant)
    print()

    print("=" * 72)
    print("IXODES SPECIMENS FROM DEER, FOR COMPARISON")
    print("=" * 72)
    # `ixodes - derm` is every animal carrying Ixodes and not Dermacentor,
    # which includes the ten black bears and is not a deer comparison. Bear
    # reads once looked alarming here for exactly that reason.
    deer = set(species.loc[[("deer" in str(h).lower()
                             or "wtd" in str(h).lower())
                            for h in species["host_bin"]], "animal"])
    deer.discard("")
    deer_libraries = animals.index[animals.isin((ixodes & deer) - derm)]
    if len(deer_libraries):
        totals = by_group.loc[deer_libraries].sum(axis=0).sort_values(ascending=False)
        carrying = (by_group.loc[deer_libraries] > 0).sum(axis=0)
        print("  %d libraries from %d specimens" % (len(deer_libraries),
              animals[deer_libraries].nunique()))
        print("  %-26s %12s %10s" % ("host group", "libraries", "reads"))
        for name in totals.index:
            if totals[name] > 0:
                print("  %-26s %12d %10d" % (name[:26], carrying[name],
                                             totals[name]))
    print()


if __name__ == "__main__":
    main()
