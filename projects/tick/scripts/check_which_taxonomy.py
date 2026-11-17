#!/usr/bin/env python3
"""check_which_taxonomy.py

Two questions, both raised by check_classifier_provenance.py. I write nothing.

QUESTION ONE: which 16S taxonomy fed the results?
-------------------------------------------------
Provenance found TWO 16S taxonomy artifacts built from DIFFERENT classifiers:

    taxonomy.qza                     <- silva-138-99-nb-classifier.qza
                                        (SILVA 138, full length)
    taxonomy_v4v5_regionmatched.qza  <- silva-138-99-nb-515f-926r-classifier.qza
                                        (SILVA 138 extracted to the amplicon)

Section 2.3 says only "SILVA 138". It has to name one, and the manuscript's
every genus-level result rests on whichever it was. So I compare the two
assignment for assignment and report how far apart they are: if they agree on
nearly every ASV the choice is cosmetic, and if they do not it is not.

QUESTION TWO: what region is each classifier trained on?
--------------------------------------------------------
A classifier built by extract_reads records the primers it was extracted with,
in its own provenance. So the artifact says what region it covers, rather than
anyone inferring it from the file name. This matters most for 18S: the taxonomy
was made by pr2-18S-classifier-v4v8.qza while PrimerTable.xlsx gives the 18S
amplicon as V9, 1391F/1510R. V9 sits downstream of V8. If both are right the
reads were classified against a reference that does not span them, and that
would reach the tick genus calls and the piroplasm calls alike.

I do not assume the file names are honest. extract_reads parameters are.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_which_taxonomy.py
"""
import os
import re
import sys
import zipfile

PAIR = ["qiime2/16S/all/taxonomy/taxonomy.qza",
        "qiime2/16S/all/taxonomy/taxonomy_v4v5_regionmatched.qza"]
CLASSIFIERS = [
    "classifiers/silva-138-99-nb-classifier.qza",
    "qiime2/16S/all/taxonomy/silva-138-99-nb-515f-926r-classifier.qza",
    "qiime2/18S/all/classifier/pr2-18S-classifier-v4v8.qza",
    "classifiers/pr2-18S-classifier-trimmed.qza",
    "classifiers/silva-138-V9-classifier.qza",
]
PARAMS = ("f_primer", "r_primer", "trunc_len", "trim_left", "min_length",
          "max_length", "identity", "region")


def taxonomy_rows(path):
    """I return {feature_id: taxon} from a FeatureData[Taxonomy] artifact."""
    with zipfile.ZipFile(path) as z:
        hit = [n for n in z.namelist() if n.endswith("data/taxonomy.tsv")]
        if not hit:
            return None
        with z.open(sorted(hit, key=len)[0]) as fh:
            lines = fh.read().decode("utf-8", "replace").splitlines()
    out = {}
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) >= 2:
            out[parts[0]] = parts[1]
    return out


def genus(taxon):
    """I return the genus field, or the deepest rank present."""
    ranks = [r.strip() for r in taxon.split(";")]
    for r in ranks:
        if r.startswith("g__") and len(r) > 3:
            return r
    return ranks[-1] if ranks else ""


def extraction_params(path):
    """I return every extract_reads-like action in this artifact's history."""
    found = []
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.endswith("action/action.yaml")]
        for n in names:
            with z.open(n) as fh:
                text = fh.read().decode("utf-8", "replace")
            got = {}
            for key in PARAMS:
                # QIIME 2 writes parameters as YAML list items, so the line
                # reads "    -   f_primer: GTGY..." and a whitespace-only
                # prefix will not reach the key. The dash is optional because
                # top-level keys carry none.
                m = re.search(
                    r'^[^\S\n]*(?:-[^\S\n]*)?%s:[^\S\n]*(\S.*?)[^\S\n]*$' % key,
                    text, flags=re.M)
                if m:
                    got[key] = m.group(1).strip().strip("'\"")
            if got:
                act = re.search(r'^[^\S\n]+action:[^\S\n]*(\S+)', text,
                                flags=re.M)
                found.append((act.group(1) if act else "?", got))
    return found


def main():
    if not os.path.isdir("qiime2"):
        sys.exit("Run me from the project root; I see no qiime2/ directory.")

    print("=" * 74)
    print("ONE: how far apart are the two 16S taxonomies?")
    print("=" * 74)
    tables = {}
    for p in PAIR:
        if not os.path.exists(p):
            print("  MISSING %s" % p)
            continue
        tables[p] = taxonomy_rows(p)
        print("  %-52s %6d ASVs" % (p.split("/")[-1], len(tables[p])))

    if len(tables) == 2:
        a, b = [tables[p] for p in PAIR if p in tables]
        shared = set(a) & set(b)
        same_full = sum(1 for k in shared if a[k] == b[k])
        same_genus = sum(1 for k in shared if genus(a[k]) == genus(b[k]))
        print()
        print("  ASVs in both                %6d" % len(shared))
        print("  identical full lineage      %6d  (%.1f%%)"
              % (same_full, 100.0 * same_full / max(len(shared), 1)))
        print("  identical genus             %6d  (%.1f%%)"
              % (same_genus, 100.0 * same_genus / max(len(shared), 1)))
        print("  DIFFER at genus             %6d" % (len(shared) - same_genus))
        print()
        print("  first 12 genus-level disagreements:")
        shown = 0
        for k in shared:
            if genus(a[k]) != genus(b[k]):
                print("    %s" % k[:34])
                print("      full-length  %s" % (genus(a[k]) or "(none)"))
                print("      region-match %s" % (genus(b[k]) or "(none)"))
                shown += 1
                if shown == 12:
                    break

    print()
    print("=" * 74)
    print("TWO: what region is each classifier trained on?")
    print("=" * 74)
    for p in CLASSIFIERS:
        print()
        print("  %s" % p)
        if not os.path.exists(p):
            print("    not on disk")
            continue
        hits = extraction_params(p)
        if not hits:
            print("    no extract_reads parameters in its provenance")
            print("    (consistent with a classifier trained on full-length "
                  "references)")
            continue
        for action, got in hits[:4]:
            print("    action %s" % action)
            for key in PARAMS:
                if key in got:
                    print("      %-11s %s" % (key, got[key]))

    print()
    print("Send this whole output back.")


if __name__ == "__main__":
    main()
