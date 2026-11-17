#!/usr/bin/env python3
"""check_funnel_floors.py  --  READ ONLY. I write nothing.

THE DISCREPANCY
---------------
Section 2.6 says the decontaminated 16S table holds 907 libraries, of which
853 resolve to an animal and 602 clear the read-depth floor, covering 413
animals. check_cohort_funnel.py reports 712 of 907 libraries at or above
1,000 reads, and a different animal count. Both numbers have sat in the
manuscript behind a marker since before 16 October.

WHY THE OBVIOUS READING DOES NOT WORK
-------------------------------------
If both used the same floor on the same table, then at most 54 of the 712
could be libraries that do not resolve to an animal, because only 54 of the
907 fail to resolve. That would leave at least 658 resolved libraries above
the floor, not 602. So the two numbers are not the same measurement and the
question is which measurement each one is.

THE LIKELY ANSWER, WHICH THIS SCRIPT TESTS RATHER THAN ASSERTS
--------------------------------------------------------------
Section 2.6 uses TWO different thousand-read floors in two different
sentences:

    "1,000 reads once host, chloroplast and eukaryotic reads were removed"
    "1,000 classified bacterial reads once the domain-only and unassigned
     reads were removed as well"

and check_cohort_funnel.py says only "1,000 reads", which most likely means
1,000 TOTAL reads with nothing removed. 515F/926R co-amplifies the tick's own
18S, and in this collection the eukaryotic fraction is frequently well over
90% of a library, so a library can hold tens of thousands of total reads and
far fewer than 1,000 bacterial ones. That gap is the right size and the right
shape to explain 712 against 602.

If that is the explanation then NEITHER number is wrong and nothing needs
recomputing: the manuscript needs each sentence to name which floor it means.
If it is not the explanation, this prints what the real one is.

WHAT I DO
---------
1. Print how check_cohort_funnel.py actually defines its floor, by reading the
   script rather than by inferring from its output. Six times this month a
   document's account of the work has turned out not to match the work.
2. Compute, for every library, three depths: total reads; reads after host,
   chloroplast and mitochondrial reads are removed; and classified bacterial
   reads.
3. Count libraries and distinct animals clearing 1,000 at each depth, for all
   907 and for the subset that resolves to an animal.
4. Say which cell of that table reproduces each published number.

I assert nothing about which number belongs in the manuscript. I print the
nine counts and name the ones that match.
"""
import os
import re
import sys
from collections import defaultdict

ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1
                          else "~/meedlab/tick_project")
TAX = "results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv"
TABLE = "results/decontam/16S/table_decontam_16S.tsv"
CROSSWALK = "results/metadata/crosswalk_VERIFIED.tsv"
FUNNEL = "scripts/check_cohort_funnel.py"
FLOOR = 1000
NONBACT = ("Mitochondria", "Chloroplast")

PUBLISHED = {"907 libraries in the table": 907,
             "853 resolve to an animal": 853,
             "602 clear the floor": 602,
             "413 animals": 413,
             "712 at or above 1,000 reads": 712}


def parse(taxon):
    d = re.search(r'd__([^;]+)', taxon)
    domain = d.group(1).strip() if d and d.group(1).strip() else "Unassigned"
    named = None
    for rank in ("g__", "f__"):
        m = re.search(re.escape(rank) + r'([^;]+)', taxon)
        if m and m.group(1).strip():
            named = m.group(1).strip()
            break
    return domain, named


def main():
    os.chdir(ROOT)
    bar = "=" * 78
    for need in (TAX, TABLE):
        if not os.path.isfile(need):
            sys.exit("missing %s" % need)

    # ---- 1. what does the other script actually measure? ----------------
    print(bar); print("1. HOW check_cohort_funnel.py DEFINES ITS FLOOR"); print(bar)
    if not os.path.isfile(FUNNEL):
        print("  %s not found. Find whatever reports 712 and read its floor." % FUNNEL)
    else:
        src = open(FUNNEL, encoding="utf-8", errors="replace").read().splitlines()
        hits = [i for i, L in enumerate(src)
                if re.search(r'1000|1,000|FLOOR|floor|sum\(|total|depth', L)]
        if not hits:
            print("  no floor-like line found; read it by hand.")
        for i in hits[:30]:
            print("  %4d | %s" % (i + 1, src[i].rstrip()[:110]))
        print()
        print("  Read those lines before trusting anything below. The question")
        print("  is whether its 1,000 is on TOTAL reads or on reads after host")
        print("  removal.")
    print()

    # ---- 2. three depths per library ------------------------------------
    tax = {}
    with open(TAX, encoding="utf-8", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        tcol = next(i for i, c in enumerate(head) if c.lower() == "taxon")
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if len(r) > tcol:
                tax[r[0]] = parse(r[tcol])

    with open(TABLE, encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
        if first.startswith("# Constructed from biom"):
            first = fh.readline()
        samples = first.rstrip("\n").lstrip("#").split("\t")[1:]
        total = [0] * len(samples)
        nonhost = [0] * len(samples)
        classbact = [0] * len(samples)
        unmapped = seen = 0
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if len(r) < 2:
                continue
            seen += 1
            hit = tax.get(r[0])
            if hit is None:
                unmapped += 1
                hit = ("Unassigned", None)
            domain, named = hit
            is_host = (domain == "Eukaryota") or (named in NONBACT)
            is_cb = (domain == "Bacteria") and named is not None and named not in NONBACT
            for k, v in enumerate(r[1:]):
                if not v or v in ("0", "0.0"):
                    continue
                try:
                    n = int(float(v))
                except ValueError:
                    continue
                if not n:
                    continue
                total[k] += n
                if not is_host:
                    nonhost[k] += n
                if is_cb:
                    classbact[k] += n

    print(bar); print("2. TABLE AND TAXONOMY AGREE?"); print(bar)
    print("  %d features in the table, %d absent from the taxonomy (%.1f%%)"
          % (seen, unmapped, 100.0 * unmapped / max(seen, 1)))
    if unmapped > seen * 0.05:
        print("  *** over 5%% unmapped: these may not be the same feature space ***")
    print("  %d libraries in the table" % len(samples))
    print()

    # ---- 3. which libraries resolve to an animal ------------------------
    animal = {}
    if os.path.isfile(CROSSWALK):
        with open(CROSSWALK, encoding="utf-8", errors="replace") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            lib_i = 0
            an_i = next((i for i, c in enumerate(head)
                         if re.search(r'specimen|sample_id|animal', c, re.I)
                         and i != lib_i), 1)
            for line in fh:
                r = line.rstrip("\n").split("\t")
                if len(r) > max(lib_i, an_i) and r[lib_i].strip():
                    animal[r[lib_i].strip()] = r[an_i].strip()
    else:
        print("  %s not found; animal counts will be blank." % CROSSWALK)

    def resolve(name):
        if name in animal:
            return animal[name]
        base = re.sub(r'\.\d{6}_.*$', '', name)
        base = re.sub(r'-(16s|16S)$', '', base)
        for k, v in animal.items():
            if k.startswith(base):
                return v
        return None

    resolved = [resolve(s) for s in samples]
    n_res = sum(1 for a in resolved if a)
    print("  %d of %d libraries resolve to an animal" % (n_res, len(samples)))
    print()

    # ---- 4. the nine counts ---------------------------------------------
    print(bar); print("3. THE FUNNEL AT EACH FLOOR DEFINITION"); print(bar)
    print("  %-34s %10s %10s %10s" % ("", "libraries", "of resolved", "animals"))
    depths = (("total reads, nothing removed", total),
              ("host/chloro/mito removed", nonhost),
              ("classified bacterial reads", classbact))
    results = {}
    for label, vec in depths:
        libs = sum(1 for v in vec if v >= FLOOR)
        libs_res = sum(1 for k, v in enumerate(vec) if v >= FLOOR and resolved[k])
        animals = len({resolved[k] for k, v in enumerate(vec)
                       if v >= FLOOR and resolved[k]})
        results[label] = (libs, libs_res, animals)
        print("  >= %-31s %10d %10d %10d" % (FLOOR, libs, libs_res, animals))
        print("     %s" % label)
    print()

    print(bar); print("4. WHICH PUBLISHED NUMBER IS WHICH"); print(bar)
    allnums = {len(samples): "libraries in the table", n_res: "resolve to an animal"}
    for label, (libs, libs_res, animals) in results.items():
        allnums.setdefault(libs, "libraries >= floor on %s" % label)
        allnums.setdefault(libs_res, "resolved libraries >= floor on %s" % label)
        allnums.setdefault(animals, "animals >= floor on %s" % label)
    for claim, value in PUBLISHED.items():
        match = allnums.get(value)
        print("  %-34s %5d   %s" % (claim, value, match or "*** NOTHING MATCHES ***"))
    print()
    print("  If every published number matches a cell above, nothing is wrong")
    print("  with the analysis and the manuscript only has to name which floor")
    print("  each sentence means. If one matches nothing, that is the one to")
    print("  chase.")
    print()
    print("Send this whole output back.")


if __name__ == "__main__":
    main()
