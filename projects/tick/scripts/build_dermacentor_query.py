#!/usr/bin/env python3
"""build_dermacentor_query.py

I collect the sequences that can tell which Dermacentor species each specimen
is, and write them as one FASTA for BLAST. Read-only unless given --write.

The 18S call is genus level. Of the 20 Dermacentor specimens, 11 are from moose
(expected D. albipictus) and 9 from opossum or porcupine (expected
D. variabilis). Two sequence sets can settle it: the abundant cytb ASVs the
vertebrate classifier left unassigned (likely tick mitochondrial cytb), and the
Dermacentor-assigned 18S ASVs, included to show what 18S can resolve.

Identical ASVs are written once, with the hosts and specimen count in the
header, so the BLAST result reads back without another lookup.

USAGE
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/build_dermacentor_query.py [--write]
"""
import glob
import os
import shutil
import sys
import tempfile
import zipfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

MARKERS = {
    "cytb": {
        "table": "qiime2/cytb/all/dada2/table.qza",
        "seqs": "qiime2/cytb/all/dada2/rep-seqs.qza",
        "taxonomy": "results/cytb/all/taxonomy_verified/*.tsv",
    },
    "18S": {
        "table": "qiime2/18S/all/dada2/table.qza",
        "seqs": "qiime2/18S/all/dada2/rep-seqs.qza",
        "taxonomy": "results/18S/all/taxonomy_verified/*.tsv",
    },
}
SPECIES = "results/analysis/tick_species_by_host.tsv"
OUT = "results/analysis/dermacentor_species_query.fasta"
PER_SPECIMEN = 2
MIN_READS = 100
VERTEBRATE_WORDS = ("chordata", "mammalia", "aves", "cervidae", "alces",
                    "odocoileus", "didelph", "erethizon", "ursus", "procyon",
                    "homo", "sylvilagus", "rangifer", "vertebrata")


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s). Run from the project root."
                 % problem)
    return tick_lookup.load(), tick_lookup


def animal_of(table, module, name, marker):
    """I turn a library or specimen name into its animal, or ''."""
    text = str(name).strip()
    stem = text.split(".")[0]
    for key, value in ((marker, text), (marker, stem), ("specimen", stem)):
        hit = module.resolve(table, key, value)
        if hit is not None:
            return str(hit)
    return ""


def extract(qza, suffix):
    """I return the path of one extracted file and its temp dir."""
    if not os.path.isfile(qza):
        sys.exit("missing %s" % qza)
    workdir = tempfile.mkdtemp(prefix="query_")
    with zipfile.ZipFile(qza) as archive:
        wanted = [n for n in archive.namelist() if n.endswith(suffix)]
        if not wanted:
            shutil.rmtree(workdir, ignore_errors=True)
            sys.exit("no %s inside %s" % (suffix, qza))
        archive.extract(wanted[0], workdir)
    return os.path.join(workdir, wanted[0]), workdir


def read_counts(qza):
    """I return the feature table as libraries x ASVs."""
    try:
        import biom
    except ImportError:
        sys.exit("biom is not importable. Run inside the qiime2 environment.")
    path, workdir = extract(qza, "data/feature-table.biom")
    try:
        table = biom.load_table(path)
        return pd.DataFrame(
            table.matrix_data.toarray().T,
            index=[str(s) for s in table.ids(axis="sample")],
            columns=[str(f) for f in table.ids(axis="observation")])
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def read_sequences(qza):
    """I return {asv: sequence} from a rep-seqs artifact."""
    path, workdir = extract(qza, "data/dna-sequences.fasta")
    sequences, name = {}, None
    try:
        with open(path) as handle:
            for line in handle:
                line = line.strip()
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    sequences[name] = ""
                elif name:
                    sequences[name] += line
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return sequences


def read_taxonomy(pattern, features):
    """I return the taxonomy covering the most features, and its path."""
    best = None
    for path in sorted(glob.glob(pattern)):
        try:
            frame = pd.read_csv(path, sep="\t", dtype=str)
        except (ValueError, OSError, pd.errors.ParserError):
            continue
        column = next((c for c in frame.columns
                       if c.lower() in ("taxon", "taxonomy", "verified_taxon",
                                        "taxon_lineage")), None)
        if column is None:
            continue
        mapping = dict(zip(frame[frame.columns[0]], frame[column]))
        covered = sum(1 for f in features if f in mapping)
        if best is None or covered > best[1]:
            best = (path, covered, mapping)
    if best is None:
        sys.exit("no taxonomy matched %s" % pattern)
    return best


def wanted_asvs(marker, counts, taxonomy, libraries):
    """I pick the ASVs to BLAST for one specimen in one marker."""
    totals = counts.loc[libraries].sum(axis=0)
    totals = totals[totals >= MIN_READS].sort_values(ascending=False)
    picked = []
    for asv in totals.index:
        text = str(taxonomy.get(asv, "")).lower()
        if marker == "cytb":
            if any(word in text for word in VERTEBRATE_WORDS):
                continue
        else:
            if "dermacentor" not in text:
                continue
        picked.append((asv, int(totals[asv])))
        if len(picked) >= PER_SPECIMEN:
            break
    return picked


def main():
    write = "--write" in sys.argv[1:]
    table, module = lookup()
    species = pd.read_csv(SPECIES, sep="\t", dtype=str)
    positive = ("1", "1.0", "True", "TRUE", "true")
    derm = species[species["Dermacentor_pos"].isin(positive)].copy()
    derm["animal"] = [animal_of(table, module, s, "16S") for s in derm["specimen"]]
    derm = derm[derm["animal"] != ""]
    hosts = dict(zip(derm["animal"], derm["host_bin"]))
    print("Dermacentor specimens: %d (%s)"
          % (len(derm), ", ".join("%s %d" % (h, n) for h, n
                                  in derm["host_bin"].value_counts().items())))

    records, seen = [], {}
    for marker, paths in MARKERS.items():
        counts = read_counts(paths["table"])
        sequences = read_sequences(paths["seqs"])
        tax_path, covered, taxonomy = read_taxonomy(paths["taxonomy"],
                                                    list(counts.columns))
        animals = pd.Series([animal_of(table, module, s, marker)
                             for s in counts.index], index=counts.index)
        print("%-4s table %d x %d, taxonomy %s (%d matched), libraries with an "
              "animal %d" % (marker, counts.shape[0], counts.shape[1],
                             os.path.basename(tax_path), covered,
                             int((animals != "").sum())))
        for animal, host in sorted(hosts.items()):
            libraries = animals.index[animals == animal]
            if not len(libraries):
                continue
            for asv, reads in wanted_asvs(marker, counts, taxonomy, libraries):
                if asv not in sequences:
                    continue
                entry = seen.setdefault((marker, asv),
                                        {"reads": 0, "hosts": set(),
                                         "specimens": set(),
                                         "sequence": sequences[asv]})
                entry["reads"] += reads
                entry["hosts"].add(str(host).replace(" ", "_"))
                entry["specimens"].add(animal)

    for (marker, asv), entry in sorted(seen.items(),
                                       key=lambda kv: -kv[1]["reads"]):
        records.append(("%s|%s|%s|specimens=%d|reads=%d"
                        % (marker, asv[:8], "+".join(sorted(entry["hosts"])),
                           len(entry["specimens"]), entry["reads"]),
                        entry["sequence"]))

    per_marker = {}
    for header, _ in records:
        key = header.split("|")[0]
        per_marker[key] = per_marker.get(key, 0) + 1
    print()
    print("queries: %s" % ", ".join("%s %d" % (m, n) for m, n in per_marker.items()))
    covered_specimens = len(set().union(*[e["specimens"] for e in seen.values()])) if seen else 0
    print("specimens represented: %d of %d" % (covered_specimens, len(hosts)))

    if not write:
        print()
        print("Nothing written. Rerun with --write to save %s" % OUT)
        return
    if os.path.exists(OUT):
        sys.exit("%s already exists. Move it aside first." % OUT)
    with open(OUT, "w") as handle:
        for header, sequence in records:
            handle.write(">%s\n%s\n" % (header, sequence))
    print("wrote %s (%d sequences)" % (OUT, len(records)))


if __name__ == "__main__":
    main()
