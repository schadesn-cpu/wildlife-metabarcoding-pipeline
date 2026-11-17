#!/usr/bin/env python3
"""build_endosymbiont_query.py

I pull the Rickettsia sequences from the deer ticks and the Francisella
sequences from the moose ticks and write them as one FASTA for BLAST.

WHY
---
The paper says these are endosymbionts rather than the pathogenic members of
their genera, and at the moment that rests on what is known about ticks rather
than on anything measured here. The BLAST verification already run recorded the
genus each ASV was checked against and the classifier's confidence, not percent
identity to named species, so it cannot settle it.

What I want is not a match to Rickettsia buchneri. A high hit to the
endosymbiont means nothing on its own. I want the identity of the SAME sequence
against the pathogenic species, because if they score alike then 16S cannot
separate them, and that is the paper's own argument stated with a number.

I take the ASVs carrying at least 1% of their genus's reads, at most six per
genus. The rest are error variants of those.

COHORTS
-------
    Rickettsia    deer with Ixodes and not Dermacentor, one library per animal
    Francisella   moose with Dermacentor, one library per animal

Both match build_core_microbiome.py, so the sequences come from the same ticks
the core microbiome tables describe.

WHAT I WRITE
------------
    results/analysis/endosymbiont_query.fasta

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/build_endosymbiont_query.py
    python3 scripts/build_endosymbiont_query.py --write
"""
import os
import shutil
import sys
import tempfile
import zipfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REP_SEQS = "qiime2/16S/all/dada2/rep-seqs.qza"
OUT = "results/analysis/endosymbiont_query.fasta"
TARGETS = ("Rickettsia", "Francisella")
MIN_SHARE = 0.01
PER_GENUS = 6


def suppression():
    try:
        import check_endosymbiont_suppression as module
    except ImportError as problem:
        sys.exit("cannot import check_endosymbiont_suppression (%s). It must be "
                 "in scripts/ beside me." % problem)
    return module


def moose_with_dermacentor(module, table, lookup_module):
    """I return the moose animals carrying Dermacentor."""
    if not os.path.isfile(module.SPECIES):
        sys.exit("missing %s" % module.SPECIES)
    species = pd.read_csv(module.SPECIES, sep="\t", dtype=str)
    keep = set()
    for specimen, host, derm in zip(species["specimen"], species["host_bin"],
                                    species["Dermacentor_pos"]):
        animal = module.animal_of(table, lookup_module, specimen)
        if animal and "moose" in str(host).lower() and module.truthy(derm):
            keep.add(animal)
    return keep


def read_sequences():
    """I return {asv: sequence} from the 16S rep-seqs artifact."""
    if not os.path.isfile(REP_SEQS):
        sys.exit("missing %s" % REP_SEQS)
    workdir = tempfile.mkdtemp(prefix="repseq_")
    sequences, name = {}, None
    try:
        with zipfile.ZipFile(REP_SEQS) as archive:
            wanted = [n for n in archive.namelist()
                      if n.endswith("data/dna-sequences.fasta")]
            if not wanted:
                sys.exit("no dna-sequences.fasta inside %s" % REP_SEQS)
            archive.extract(wanted[0], workdir)
            with open(os.path.join(workdir, wanted[0])) as handle:
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


def one_per_animal(module, counts, table, lookup_module, animals_wanted):
    """I return the deepest library per wanted animal, above the read floor."""
    animals = pd.Series([module.animal_of(table, lookup_module, s)
                         for s in counts.index], index=counts.index)
    depth = counts.sum(axis=1)
    frame = pd.DataFrame({"animal": animals, "depth": depth})
    frame = frame[(frame["animal"] != "") & (frame["depth"] >= module.MIN_READS)]
    frame = frame.sort_values("depth", ascending=False).drop_duplicates("animal")
    return frame[frame["animal"].isin(animals_wanted)]


def main():
    write = "--write" in sys.argv[1:]
    module = suppression()
    arms = module.arms()
    table, lookup_module = module.lookup()
    counts = module.read_counts(table, lookup_module)
    path, covered, mapping = arms.read_taxonomy(list(counts.columns))
    print("counts: %d libraries x %d features" % counts.shape)
    print("taxonomy: %s (%d of %d)" % (path, covered, counts.shape[1]))
    genus = pd.Series({f: arms.genus_label(mapping.get(f, ""))
                       for f in counts.columns})

    deer = module.deer_with_ixodes(table, lookup_module)
    moose = moose_with_dermacentor(module, table, lookup_module)
    cohorts = {"Rickettsia": one_per_animal(module, counts, table,
                                            lookup_module, deer),
               "Francisella": one_per_animal(module, counts, table,
                                             lookup_module, moose)}
    print("deer with Ixodes: %d; moose with Dermacentor: %d"
          % (len(cohorts["Rickettsia"]), len(cohorts["Francisella"])))
    print()

    sequences = read_sequences()
    records = []
    for target in TARGETS:
        chosen = cohorts[target]
        if chosen.empty:
            print("%s: no libraries in this cohort." % target)
            continue
        features = [f for f in counts.columns if genus[f] == target]
        if not features:
            print("%s: no ASV carries this genus label." % target)
            continue
        part = counts.loc[chosen.index, features]
        totals = part.sum(axis=0).sort_values(ascending=False)
        overall = float(totals.sum())
        print("%s: %d ASVs, %d reads across %d animals"
              % (target, len(features), int(overall), len(chosen)))
        kept = 0
        for asv in totals.index:
            share = totals[asv] / overall if overall else 0
            if share < MIN_SHARE or kept >= PER_GENUS:
                break
            if asv not in sequences:
                print("  %s has no representative sequence, skipping" % asv[:8])
                continue
            animals = int((part[asv] > 0).sum())
            print("  %s  %9d reads (%4.1f%%)  in %d animals"
                  % (asv[:8], int(totals[asv]), 100 * share, animals))
            records.append(("16S|%s|%s|animals=%d|reads=%d"
                            % (target, asv[:8], animals, int(totals[asv])),
                            sequences[asv]))
            kept += 1
        print()

    if not records:
        sys.exit("nothing to write. Check the genus labels in the taxonomy.")
    print("queries: %d" % len(records))
    if not write:
        print()
        print("Nothing written. Rerun with --write to save %s" % OUT)
        return
    if os.path.exists(OUT):
        sys.exit("%s already exists. Move it aside first." % OUT)
    with open(OUT, "w") as handle:
        for header, sequence in records:
            handle.write(">%s\n%s\n" % (header, sequence))
    print("wrote %s" % OUT)
    print()
    print("BLAST it at https://blast.ncbi.nlm.nih.gov, nucleotide blast:")
    print("  database core_nt, program megablast, organism box EMPTY")
    print("  then send the Descriptions tab, not the hit table, so the species")
    print("  names come through.")


if __name__ == "__main__":
    main()
