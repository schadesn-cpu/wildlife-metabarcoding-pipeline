#!/usr/bin/env python3
"""check_terminal_base_all_markers.py

I run the terminal-base check on every marker, not just 16S. I write nothing.

WHY
---
In the 16S data, trimming one base from the end of each amplicon made 5,204 of
89,970 variants identical to another variant, and those variants carried 76.9%
of all reads. Trimming a second base merged only 244 more, so the artifact sits
at the terminal position. The same thing would split the tick species calls and
the blood meal calls, so I check 18S, cytb and ITS the same way.

HOW TO READ IT
--------------
Trim 0 is the control. A large jump from trim 0 to trim 1 and a small one from
trim 1 to trim 2 means the last base is the problem. A steady decline across
all four means sequences simply resemble each other near the end and nothing
marker-specific is happening.

I weight by reads as well as by variant count, because a thousand one-read
variants merging matters less than two variants holding millions.

WHAT I READ, PER MARKER
-----------------------
    qiime2/<marker>/all/dada2/rep-seqs.qza
    qiime2/<marker>/all/dada2/table.qza
    results/<marker>/all/taxonomy_verified/*.tsv   for labels, if present

A marker whose files are missing is reported and skipped.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_terminal_base_all_markers.py
    python3 scripts/check_terminal_base_all_markers.py 18S cytb
"""
import glob
import os
import shutil
import sys
import tempfile
import zipfile

import pandas as pd

MARKERS = ("16S", "18S", "cytb", "ITS")
TRIMS = (0, 1, 2, 3)
SHOW = 8


def extract(qza, suffix):
    """I return the path of one extracted file and its temporary directory."""
    workdir = tempfile.mkdtemp(prefix="marker_")
    with zipfile.ZipFile(qza) as archive:
        wanted = [n for n in archive.namelist() if n.endswith(suffix)]
        if not wanted:
            shutil.rmtree(workdir, ignore_errors=True)
            return None, None
        archive.extract(wanted[0], workdir)
    return os.path.join(workdir, wanted[0]), workdir


def read_sequences(qza):
    """I return {variant: sequence} from a rep-seqs artifact."""
    path, workdir = extract(qza, "data/dna-sequences.fasta")
    if path is None:
        sys.exit("no dna-sequences.fasta inside %s" % qza)
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


def read_totals(qza):
    """I return {variant: total reads} from a feature table artifact."""
    try:
        import biom
    except ImportError:
        sys.exit("biom is not importable. Run me inside the qiime2 environment.")
    path, workdir = extract(qza, "data/feature-table.biom")
    if path is None:
        sys.exit("no feature-table.biom inside %s" % qza)
    try:
        table = biom.load_table(path)
        return dict(zip([str(f) for f in table.ids(axis="observation")],
                        table.sum(axis="observation")))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def read_taxonomy(marker, features):
    """I return {variant: label} from the marker's verified taxonomy, or {}."""
    pattern = "results/%s/all/taxonomy_verified/*.tsv" % marker
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
    return best[2] if best else {}


def label_of(lineage):
    """I return the lowest named rank of a lineage, without its prefix."""
    text = str(lineage)
    if text in ("", "nan", "None"):
        return "unassigned"
    for piece in reversed(text.split(";")):
        piece = piece.strip()
        if not piece:
            continue
        if "__" in piece:
            piece = piece.split("__", 1)[1].strip()
        if piece:
            return piece
    return "unassigned"


def collapse(sequences, totals, trim):
    """I group variants by sequence with `trim` bases cut from the end."""
    groups = {}
    for name, sequence in sequences.items():
        key = sequence[:-trim] if trim else sequence
        groups.setdefault(key, []).append(name)
    merged = {k: v for k, v in groups.items() if len(v) > 1}
    absorbed = sum(len(v) - 1 for v in merged.values())
    reads_in = sum(float(totals.get(a, 0)) for v in merged.values() for a in v)
    return len(groups), merged, absorbed, reads_in


def run(marker):
    """I report one marker."""
    seqs_qza = "qiime2/%s/all/dada2/rep-seqs.qza" % marker
    table_qza = "qiime2/%s/all/dada2/table.qza" % marker
    print("=" * 74)
    print(marker.upper())
    print("=" * 74)
    missing = [p for p in (seqs_qza, table_qza) if not os.path.isfile(p)]
    if missing:
        print("  skipped, missing: %s" % ", ".join(missing))
        print()
        return
    sequences = read_sequences(seqs_qza)
    totals = read_totals(table_qza)
    lengths = pd.Series([len(s) for s in sequences.values()])
    everything = float(sum(totals.values()))
    print("  variants %d, amplicon length median %d (range %d to %d), reads %.0f"
          % (len(sequences), lengths.median(), lengths.min(), lengths.max(),
             everything))
    if not everything:
        print("  the feature table holds no reads. Nothing to weigh.")
        print()
        return

    results = {}
    print("  %5s %12s %13s %14s %11s"
          % ("trim", "distinct", "variants", "reads merged", "% of reads"))
    for trim in TRIMS:
        distinct, merged, absorbed, reads_in = collapse(sequences, totals, trim)
        results[trim] = (merged, absorbed, reads_in)
        print("  %5d %12d %13d %14.0f %10.1f%%"
              % (trim, distinct, absorbed, reads_in,
                 100 * reads_in / everything))
    first = results[1][1]
    second = results[2][1] - first
    print("  first base merges %d variants, the second merges %d more"
          % (first, second))
    print()

    merged = results[1][0]
    if not merged:
        print("  nothing merges on the terminal base.")
        print()
        return
    taxonomy = read_taxonomy(marker, list(sequences))
    ranked = sorted(merged.items(),
                    key=lambda kv: -sum(float(totals.get(a, 0)) for a in kv[1]))
    print("  biggest groups merging on the terminal base alone:")
    print("    %-28s %6s %13s  %s" % ("label", "variants", "reads", "members"))
    for _, members in ranked[:SHOW]:
        labels = set(label_of(taxonomy.get(a, "")) for a in members)
        label = sorted(labels)[0] if len(labels) == 1 \
            else "mixed: " + ", ".join(sorted(labels)[:2])
        total = sum(float(totals.get(a, 0)) for a in members)
        print("    %-28s %6d %13.0f  %s"
              % (label[:28], len(members), total,
                 " ".join(a[:8] for a in members[:3])))
    print()


def main():
    wanted = [a for a in sys.argv[1:] if not a.startswith("-")] or list(MARKERS)
    unknown = [m for m in wanted if m not in MARKERS]
    if unknown:
        sys.exit("I do not know the marker(s) %s. I know: %s"
                 % (", ".join(unknown), ", ".join(MARKERS)))
    if not os.path.isdir("qiime2"):
        sys.exit("no qiime2/ directory here. Run me from the project root.")
    for marker in wanted:
        run(marker)


if __name__ == "__main__":
    main()
