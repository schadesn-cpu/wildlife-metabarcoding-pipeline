#!/usr/bin/env python3
"""check_16s_genera.py

The two 16S taxonomies disagree at genus for 41,015 of 89,970 ASVs. I ask the
only question that decides what to do about it: does the disagreement touch
the genera this paper reports? I write nothing.

WHY THIS IS THE QUESTION
------------------------
46% disagreement sounds fatal and may be irrelevant. The examples seen so far
are mostly low-information churn -- Unassigned against d__Bacteria,
g__Mitochondria against Unassigned -- on ASVs the paper never counts. What the
paper counts is a short list of genera, in a cohort of 330 deer and 11 moose,
each animal represented by its deepest library.

So the test is not "how different are the tables" but "would any reported
number move". Three things decide that, and I report all three:

  1. For each reported genus, how many ASVs each taxonomy assigns to it, and
     how many they agree on. A genus both tables call the same ASVs is safe
     whichever table was used.
  2. Where the disagreements go. An ASV that one table calls Rickettsia and
     the other calls nothing is a different problem from one the two tables
     assign to different real genera.
  3. The abundance at stake: what share of reads sits on the ASVs that
     disagree. One disputed ASV holding 6 million reads matters; four hundred
     holding twelve do not.

I take read counts from the feature table when I can find one, because an ASV
count without abundance cannot tell you whether a result would move.

WHAT I DO NOT DO
----------------
Choose. If the reported genera are untouched, the fix is a sentence in Section
2.3 naming which artifact was used. If they are touched, it is a rerun, and
that is the author's call with the numbers in front of her.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_16s_genera.py
"""
import os
import sys
import zipfile

PAIR = [("full-length", "qiime2/16S/all/taxonomy/taxonomy.qza"),
        ("region-matched",
         "qiime2/16S/all/taxonomy/taxonomy_v4v5_regionmatched.qza")]

# every genus or family the manuscript reports a number for
REPORTED = [
    "Rickettsia", "Anaplasma", "Borreliella", "Borrelia", "Francisella",
    "Bartonella", "Arsenophonus", "Rickettsiella", "Diplorickettsia",
    "Coxiella", "Ehrlichia", "Escherichia-Shigella", "Cutibacterium",
    "Bacillus", "Pseudomonas", "Mycobacterium", "Luteibacter",
    "Enhydrobacter", "Alkanindiges", "Acinetobacter", "Deinococcus",
    "Chryseobacterium", "Sphingobacterium", "Nocardia",
]
TABLE_GUESSES = [
    "qiime2/16S/all/table/table.qza",
    "qiime2/16S/all/dada2/table.qza",
    "qiime2/16S/all/table.qza",
]


def taxonomy_rows(path):
    with zipfile.ZipFile(path) as z:
        hit = [n for n in z.namelist() if n.endswith("data/taxonomy.tsv")]
        if not hit:
            return {}
        with z.open(sorted(hit, key=len)[0]) as fh:
            lines = fh.read().decode("utf-8", "replace").splitlines()
    out = {}
    for line in lines[1:]:
        p = line.split("\t")
        if len(p) >= 2:
            out[p[0]] = p[1]
    return out


def genus_of(taxon):
    for rank in (r.strip() for r in taxon.split(";")):
        if rank.startswith("g__") and len(rank) > 3:
            return rank[3:]
    return None


def deepest(taxon):
    ranks = [r.strip() for r in taxon.split(";") if r.strip()]
    return ranks[-1] if ranks else "(empty)"


def read_counts():
    """I return {feature_id: total reads} from a biom table, or {}."""
    path = next((p for p in TABLE_GUESSES if os.path.exists(p)), None)
    if not path:
        return {}, None
    try:
        import h5py                                       # noqa: PLC0415
    except ImportError:
        return {}, path + "  (h5py not available, abundance skipped)"
    try:
        with zipfile.ZipFile(path) as z:
            hit = [n for n in z.namelist() if n.endswith(".biom")]
            if not hit:
                return {}, path + "  (no biom inside)"
            z.extract(hit[0], "/tmp")
            local = os.path.join("/tmp", hit[0])
        with h5py.File(local, "r") as f:
            ids = [i.decode() for i in f["observation/ids"][:]]
            import numpy as np                            # noqa: PLC0415
            data = f["observation/matrix/data"][:]
            indptr = f["observation/matrix/indptr"][:]
            sums = [float(np.sum(data[indptr[i]:indptr[i + 1]]))
                    for i in range(len(ids))]
        return dict(zip(ids, sums)), path
    except Exception as exc:                              # noqa: BLE001
        return {}, "%s  (unreadable: %s)" % (path, exc)


def main():
    if not os.path.isdir("qiime2"):
        sys.exit("Run me from the project root; I see no qiime2/ directory.")

    tables = {}
    for label, path in PAIR:
        if not os.path.exists(path):
            sys.exit("missing %s" % path)
        tables[label] = taxonomy_rows(path)
        print("%-15s %6d ASVs   %s" % (label, len(tables[label]), path))

    counts, where = read_counts()
    print("read counts    %s" % (where or "no feature table found"))
    total_reads = sum(counts.values()) if counts else 0

    a, b = tables["full-length"], tables["region-matched"]
    shared = set(a) & set(b)

    print()
    print("=" * 78)
    print("DOES THE DISAGREEMENT TOUCH THE GENERA THE PAPER REPORTS?")
    print("=" * 78)
    print("%-22s %7s %7s %7s %9s" %
          ("genus", "full", "region", "agree", "reads@risk"))
    trouble = []
    for g in REPORTED:
        ia = {k for k in shared if genus_of(a[k]) == g}
        ib = {k for k in shared if genus_of(b[k]) == g}
        both = ia & ib
        risk = (ia ^ ib)
        reads = sum(counts.get(k, 0) for k in risk) if counts else 0
        share = (100.0 * reads / total_reads) if total_reads else 0
        print("%-22s %7d %7d %7d %8.3f%%"
              % (g, len(ia), len(ib), len(both), share))
        if ia != ib:
            trouble.append((g, ia, ib, reads, share))

    print()
    print("=" * 78)
    print("WHERE THE DISPUTED ASVs GO")
    print("=" * 78)
    if not trouble:
        print("  None. Every reported genus is assigned the same ASVs by both")
        print("  taxonomies, so the choice cannot move a reported number.")
    for g, ia, ib, reads, share in sorted(trouble, key=lambda t: -t[4])[:8]:
        print()
        print("  %s   %.3f%% of reads in dispute" % (g, share))
        only_a = sorted(ia - ib, key=lambda k: -counts.get(k, 0))[:4]
        only_b = sorted(ib - ia, key=lambda k: -counts.get(k, 0))[:4]
        for k in only_a:
            print("    full-length says %s, region-matched says %s   (%d reads)"
                  % (g, deepest(b[k]), counts.get(k, 0)))
        for k in only_b:
            print("    region-matched says %s, full-length says %s   (%d reads)"
                  % (g, deepest(a[k]), counts.get(k, 0)))

    print()
    print("Send this whole output back.")


if __name__ == "__main__":
    main()
