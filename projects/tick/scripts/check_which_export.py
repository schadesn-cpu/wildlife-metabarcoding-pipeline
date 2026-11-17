#!/usr/bin/env python3
"""check_which_export.py

I find out which 16S taxonomy the analysis actually read, by matching every
exported taxonomy table on disk against the two artifacts. I write nothing.

THE QUESTION
------------
Two 16S taxonomies exist:

    taxonomy.qza                     full-length SILVA 138
    taxonomy_v4v5_regionmatched.qza  SILVA 138 extracted to 515F/926R

The arithmetic already says the manuscript is using the full-length one:
Section 3.4 reports Rickettsiella at 316,134 reads, and the region-matched
taxonomy pulls at least 532,586 further reads into that genus from
f__Diplorickettsiaceae. But that is inference from one number. This checks it.

It also answers the question that decides the cost of switching: was the
region-matched taxonomy ever EXPORTED? If the pipeline only ever exported
taxonomy.qza, then switching means rerunning the export, the BLAST
verification and the genus tables, not relabelling a column.

HOW
---
Every exported taxonomy is a TSV of feature id and taxon string. I read both
artifacts, then every candidate TSV under results/, and report what fraction
of each TSV's rows match each artifact exactly.

A BLAST-verified table will not match either at 100%, because verification
rewrote some assignments on purpose. That is fine and is the point of
reporting both fractions side by side: the gap between them is what names the
parent, not the absolute number. I report the gap.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_which_export.py
"""
import os
import sys
import time
import zipfile

ARTIFACTS = [
    ("full-length", "qiime2/16S/all/taxonomy/taxonomy.qza"),
    ("region-matched",
     "qiime2/16S/all/taxonomy/taxonomy_v4v5_regionmatched.qza"),
]
ROOTS = ["results", "qiime2"]
MAXBYTES = 60_000_000


def from_qza(path):
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
            out[p[0]] = p[1].strip()
    return out


def from_tsv(path):
    try:
        if os.path.getsize(path) > MAXBYTES:
            return None, "too large"
        with open(path, "r", errors="replace") as fh:
            head = fh.readline()
            cols = [c.strip().lower() for c in head.split("\t")]
            # which column holds the taxon string?
            ti = None
            for cand in ("taxon", "taxonomy", "taxon_verified",
                         "taxon_final", "consensus_taxonomy"):
                if cand in cols:
                    ti = cols.index(cand)
                    break
            if ti is None:
                return None, "no taxon column (%s)" % ",".join(cols[:6])
            out = {}
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) > ti and p[0]:
                    out[p[0]] = p[ti].strip()
            return out, None
    except Exception as exc:                              # noqa: BLE001
        return None, str(exc)


def main():
    if not os.path.isdir("results"):
        sys.exit("Run me from the project root; I see no results/ directory.")

    ref = {}
    for label, path in ARTIFACTS:
        if not os.path.exists(path):
            sys.exit("missing %s" % path)
        ref[label] = from_qza(path)
        print("%-15s %6d rows   %s" % (label, len(ref[label]), path))

    cands = []
    for root in ROOTS:
        for here, dirs, names in os.walk(root):
            dirs[:] = [d for d in dirs if d != "__pycache__"]
            for n in names:
                if n.endswith((".tsv", ".txt")) and "tax" in n.lower():
                    cands.append(os.path.join(here, n))
    cands.sort()
    print("%d candidate exported taxonomy files" % len(cands))

    print()
    print("=" * 78)
    print("WHICH ARTIFACT DOES EACH EXPORT COME FROM?")
    print("=" * 78)
    verdicts = []
    for path in cands:
        rows, why = from_tsv(path)
        if rows is None:
            continue
        if not rows:
            continue
        scores = {}
        for label in ref:
            shared = set(rows) & set(ref[label])
            if not shared:
                scores[label] = 0.0
                continue
            same = sum(1 for k in shared if rows[k] == ref[label][k])
            scores[label] = 100.0 * same / len(shared)
        full, region = scores["full-length"], scores["region-matched"]
        gap = full - region
        if abs(gap) < 1.0:
            call = "ambiguous"
        elif gap > 0:
            call = "FULL-LENGTH"
        else:
            call = "REGION-MATCHED"
        verdicts.append((path, call, full, region, len(rows)))
        print()
        print("  %s" % path)
        print("    %s rows, modified %s"
              % (format(len(rows), ","),
                 time.strftime("%Y-%m-%d", time.localtime(os.path.getmtime(path)))))
        print("    matches full-length %5.1f%%   region-matched %5.1f%%   -> %s"
              % (full, region, call))

    print()
    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    fl = [v for v in verdicts if v[1] == "FULL-LENGTH"]
    rm = [v for v in verdicts if v[1] == "REGION-MATCHED"]
    print("  exports traceable to full-length    : %d" % len(fl))
    print("  exports traceable to region-matched : %d" % len(rm))
    print("  ambiguous                           : %d"
          % (len(verdicts) - len(fl) - len(rm)))
    print()
    if rm:
        print("  The region-matched taxonomy HAS been exported:")
        for p, _, _, _, n in rm:
            print("      %s  (%s rows)" % (p, format(n, ",")))
        print("  Switching would be a change of input file, not a rerun.")
    else:
        print("  The region-matched taxonomy has NEVER been exported.")
        print("  Switching means rerunning the export, the BLAST verification")
        print("  and the genus tables that depend on them -- not a relabel.")
    print()
    print("Send this whole output back.")


if __name__ == "__main__":
    main()
