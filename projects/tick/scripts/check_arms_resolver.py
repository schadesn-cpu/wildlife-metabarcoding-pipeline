#!/usr/bin/env python3
"""check_arms_resolver.py  --  READ ONLY. I write nothing and I import nothing.

THE QUESTION
------------
Section 3.3's numbers come from build_core_microbiome.py, which does not name
a taxonomy. It calls check_batch_arms.read_taxonomy(), and that function picks
a file at run time:

    for path in TAXONOMY_CANDIDATES + sorted(glob.glob(
            "results/**/taxonomy*16S*.tsv", recursive=True)):
        ...
        tax_col = next((c for c in frame.columns if c.lower() == "taxon"), None)
        if tax_col is None:
            continue                      # <- silent skip
        covered = sum(1 for f in features if f in mapping)
        if best is None or covered > best[1]:
            best = (path, covered, mapping)
        if covered == len(features):
            break                         # <- first full-coverage file wins

So which taxonomy the paper's core-microbiome numbers rest on is a property of
the filesystem on the day the script ran, not of anything written down.

WHY I DO NOT NEED THE FEATURE LIST TO ANSWER IT
-----------------------------------------------
The comparison is `covered > best[1]` -- STRICTLY greater. So a later file
never displaces an earlier one on a tie, whatever the feature list is. And
the three verified tables all carry the same 89,970 feature ids. If their id
sets are identical, every candidate scores the same coverage against any
feature list at all, the first one wins, the break fires, and the answer is
positional: candidates[0], provided it survives the taxon-column test.

That makes this decidable from the files alone, with no feature table and no
assumption about what build_core_microbiome loads. I check:

  1. the real candidate order, parsed out of check_batch_arms.py rather than
     retyped here, so this cannot drift from the script it is auditing
  2. which files the glob adds, in the order sorted() gives them
  3. for each candidate: does it have an exactly-'taxon' column, or is it
     silently skipped
  4. whether the id sets are identical, which is what makes coverage a tie
  5. the winner, and whether anything could displace it

I also report the two hazards even when today's answer is benign, because
both change the answer without changing any script:
  - a renamed taxon column silently drops a candidate
  - a new results/**/taxonomy*16S*.tsv file joins the pool on its own
"""
import ast
import glob
import hashlib
import os
import sys

ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1
                          else "~/meedlab/tick_project")
ARMS = os.path.join(ROOT, "scripts", "check_batch_arms.py")
PATTERN = "results/**/taxonomy*16S*.tsv"


def parse_candidates(path):
    """I read TAXONOMY_CANDIDATES out of the module's AST, not by retyping it."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        tree = ast.parse(fh.read(), filename=path)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == "TAXONOMY_CANDIDATES":
                try:
                    value = ast.literal_eval(node.value)
                except ValueError:
                    return None, "TAXONOMY_CANDIDATES is not a literal list"
                return list(value), None
    return None, "no TAXONOMY_CANDIDATES assignment found"


def header_and_ids(path):
    """First line for columns, first column for ids. I never load the frame."""
    ids = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        first = fh.readline().rstrip("\n").rstrip("\r")
        cols = first.split("\t")
        for line in fh:
            if not line.strip():
                continue
            ids.add(line.split("\t", 1)[0])
    return cols, ids


def main():
    if not os.path.isfile(ARMS):
        sys.exit("cannot find %s" % ARMS)
    os.chdir(ROOT)                      # the glob in the real code is relative

    candidates, err = parse_candidates(ARMS)
    if err:
        sys.exit("parse failed: %s\nRead the list by hand and rerun." % err)

    globbed = sorted(glob.glob(PATTERN, recursive=True))
    pool = list(candidates) + globbed

    print("=" * 78)
    print("THE CANDIDATE POOL, IN THE ORDER read_taxonomy() WALKS IT")
    print("=" * 78)
    print("  TAXONOMY_CANDIDATES, parsed from check_batch_arms.py:")
    for i, p in enumerate(candidates):
        print("    [%d] %s" % (i, p))
    print("  then sorted(glob('%s')):" % PATTERN)
    for p in globbed:
        print("        %s" % p)
    print()

    print("=" * 78)
    print("EACH CANDIDATE: IS IT ELIGIBLE, AND WHAT IDS DOES IT CARRY?")
    print("=" * 78)
    rows, seen = [], []
    for p in pool:
        if any(q[0] == p for q in rows):
            continue                      # the same path appears twice in the pool
        if not os.path.isfile(p):
            print("  %-62s MISSING" % p)
            rows.append((p, None, None, "missing"))
            continue
        if ".bak" in p or "_archive" in p:
            print("  %-62s skipped by the .bak/_archive filter" % p)
            rows.append((p, None, None, "filtered"))
            continue
        cols, ids = header_and_ids(p)
        taxcol = next((c for c in cols if c.lower() == "taxon"), None)
        digest = hashlib.sha1(
            "\n".join(sorted(ids)).encode("utf-8")).hexdigest()[:12]
        status = "ELIGIBLE" if taxcol else "SILENTLY SKIPPED (no 'taxon' column)"
        print("  %s" % p)
        print("      %-36s %s" % (status, "taxon col: %r" % taxcol))
        print("      %d ids, id-set %s" % (len(ids), digest))
        print("      columns: %s" % ", ".join(cols[:8]))
        rows.append((p, ids, digest, "eligible" if taxcol else "no-taxon-col"))
        if taxcol:
            seen.append((p, ids, digest))
    print()

    print("=" * 78)
    print("DO THE ELIGIBLE FILES TIE ON COVERAGE?")
    print("=" * 78)
    if not seen:
        print("  Nothing is eligible. read_taxonomy() would exit with")
        print("  'no 16S taxonomy table matches the rarefied table's features.'")
        return
    groups = {}
    for p, ids, digest in seen:
        groups.setdefault(digest, []).append((p, len(ids)))
    for digest, members in groups.items():
        print("  id-set %s  (%d ids)" % (digest, members[0][1]))
        for p, _ in members:
            print("      %s" % p)
    print()
    first_path, first_ids, first_digest = seen[0]
    bigger = [(p, len(ids)) for p, ids, d in seen[1:]
              if not ids.issubset(first_ids)]

    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print("  First eligible candidate: %s" % first_path)
    if not bigger:
        print()
        print("  No later candidate carries an id the first one lacks, so no")
        print("  later candidate can ever score higher coverage. The comparison")
        print("  is strictly-greater, so ties never displace. THIS FILE ALWAYS")
        print("  WINS, for any feature list.")
        print()
        print("  => Section 3.3's core-microbiome numbers rest on:")
        print("     %s" % first_path)
    else:
        print()
        print("  *** The following carry ids the first candidate lacks, so the")
        print("      winner depends on which features are passed in: ***")
        for p, n in bigger:
            print("        %-58s %d ids" % (p, n))
        print("      This cannot be settled statically. Capture the line")
        print("      build_core_microbiome.py prints at its line 186.")

    skipped = [p for p, ids, d, st in rows if st == "no-taxon-col"]
    print()
    print("  HAZARDS, true regardless of today's answer:")
    if skipped:
        print("    These sit in the pool and are dropped in silence, because")
        print("    the column test demands exactly 'taxon':")
        for p in skipped:
            print("      %s" % p)
    print("    Renaming a taxon column, or dropping a new")
    print("    results/**/taxonomy*16S*.tsv into the tree, changes which")
    print("    taxonomy the paper's numbers come from without changing one")
    print("    line of any script. Nothing records the choice: line 186")
    print("    prints it and no caller captures stdout.")
    print()
    print("    Fix worth making before submission: have read_taxonomy() write")
    print("    the path it chose next to its outputs, or pass the taxonomy in")
    print("    explicitly so the paper's provenance is a declaration rather")
    print("    than a search.")
    print()
    print("Send this whole output back.")


if __name__ == "__main__":
    main()
