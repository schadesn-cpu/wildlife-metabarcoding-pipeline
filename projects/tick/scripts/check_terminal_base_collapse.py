#!/usr/bin/env python3
"""check_terminal_base_collapse.py

I count how many 16S ASVs are duplicates of each other once the last base or
two of the amplicon is removed. I write nothing.

WHY
---
The two dominant Rickettsia ASVs differ at exactly one position, the terminal
base of the 377 bp amplicon, and the two dominant Francisella ASVs do the same
at the terminal base of their 382 bp amplicon. The terminal base is the
lowest-quality position in a merged read, so those are inference artifacts
rather than organisms. The question this answers is how widespread that is.

HOW TO READ THE RESULT
----------------------
I trim 0, 1, 2 and 3 bases and count distinct sequences each time. Trimming 0
is the control. If trimming ONE base collapses many ASVs and trimming a second
adds little, the problem sits at the terminal base specifically. If each extra
base collapses about as many again, sequences simply resemble each other and
nothing special is happening at the end.

I weight everything by reads as well as by ASV count, because thousands of
one-read ASVs collapsing matters far less than two ASVs carrying six million
reads each.

WHAT I READ
-----------
    qiime2/16S/all/dada2/rep-seqs.qza
    results/decontam/16S/table_decontam_16S.tsv, for read weights
    the 16S taxonomy, through check_batch_arms

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_terminal_base_collapse.py
"""
import os
import shutil
import sys
import tempfile
import zipfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REP_SEQS = "qiime2/16S/all/dada2/rep-seqs.qza"
TRIMS = (0, 1, 2, 3)
SHOW = 12


def suppression():
    try:
        import check_endosymbiont_suppression as module
    except ImportError as problem:
        sys.exit("cannot import check_endosymbiont_suppression (%s). It must be "
                 "in scripts/ beside me." % problem)
    return module


def read_sequences():
    """I return {asv: sequence} from the 16S rep-seqs artifact."""
    if not os.path.isfile(REP_SEQS):
        sys.exit("missing %s. Run me from the project root." % REP_SEQS)
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


def collapse(sequences, reads, trim):
    """I group ASVs by sequence with `trim` bases removed from the end."""
    groups = {}
    for asv, sequence in sequences.items():
        key = sequence[:-trim] if trim else sequence
        groups.setdefault(key, []).append(asv)
    merged = {k: v for k, v in groups.items() if len(v) > 1}
    absorbed = sum(len(v) - 1 for v in merged.values())
    reads_in = sum(float(reads.get(a, 0)) for v in merged.values() for a in v)
    return groups, merged, absorbed, reads_in


def main():
    module = suppression()
    arms = module.arms()
    table, lookup_module = module.lookup()
    sequences = read_sequences()
    print("representative sequences: %d" % len(sequences))
    lengths = pd.Series([len(s) for s in sequences.values()])
    print("amplicon length: median %d, range %d to %d"
          % (lengths.median(), lengths.min(), lengths.max()))

    counts = module.read_counts(table, lookup_module)
    reads = counts.sum(axis=0).to_dict()
    total_reads = float(sum(reads.values()))
    covered = sum(1 for a in sequences if a in reads)
    print("features in the decontaminated table: %d (%d of them have a "
          "representative sequence)" % (len(reads), covered))
    print("total reads: %.0f" % total_reads)
    print()

    _, _, mapping = arms.read_taxonomy(list(counts.columns))
    genus = {f: arms.genus_label(mapping.get(f, "")) for f in counts.columns}

    print("=" * 74)
    print("HOW MANY ASVs ARE DUPLICATES ONCE THE END IS TRIMMED")
    print("=" * 74)
    print("  %5s %12s %12s %14s %12s"
          % ("trim", "distinct", "ASVs merged", "reads merged", "% of reads"))
    results = {}
    for trim in TRIMS:
        groups, merged, absorbed, reads_in = collapse(sequences, reads, trim)
        results[trim] = (groups, merged, absorbed, reads_in)
        print("  %5d %12d %12d %14.0f %11.1f%%"
              % (trim, len(groups), absorbed, reads_in,
                 100 * reads_in / total_reads if total_reads else 0))
    print()
    one = results[1][2]
    two = results[2][2]
    print("  Trimming the first base merges %d ASVs. Trimming a second merges"
          % one)
    print("  %d more. If the first number is much the larger, the terminal"
          % max(two - one, 0))
    print("  base is where the artifact sits.")
    print()

    print("=" * 74)
    print("THE BIGGEST GROUPS THAT MERGE ON THE TERMINAL BASE ALONE")
    print("=" * 74)
    merged = results[1][1]
    ranked = sorted(merged.items(),
                    key=lambda kv: -sum(float(reads.get(a, 0)) for a in kv[1]))
    print("  %-22s %6s %14s  %s" % ("genus", "ASVs", "reads", "members"))
    for _, members in ranked[:SHOW]:
        labels = set(genus.get(a, "") for a in members)
        label = sorted(labels)[0] if len(labels) == 1 else "mixed: " + \
            ", ".join(sorted(x for x in labels if x)[:2])
        total = sum(float(reads.get(a, 0)) for a in members)
        print("  %-22s %6d %14.0f  %s"
              % (label[:22], len(members), total,
                 " ".join(a[:8] for a in members[:4])))
    print()

    families = {}
    for members in merged.values():
        for asv in members:
            name = genus.get(asv, "")
            if not name:
                continue
            families.setdefault(name, [0, 0.0])
            families[name][0] += 1
            families[name][1] += float(reads.get(asv, 0))
    print("  genera most affected, by reads in merging ASVs:")
    for name, (n, total) in sorted(families.items(), key=lambda kv: -kv[1][1])[:SHOW]:
        print("    %-24s %5d ASVs %14.0f reads" % (name[:24], n, total))
    print()


if __name__ == "__main__":
    main()
