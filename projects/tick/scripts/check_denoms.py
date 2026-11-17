#!/usr/bin/env python3
"""
check_denoms.py -- READ ONLY. Writes nothing, changes nothing.

The three denominators that are counted two ways. For each one it reports the
count under every plausible rule rather than picking one, because the problem
is not the arithmetic, it is that nobody recorded which rule produced which
number.

  1. TICK GENUS. 404 (Table 1 column sum, and Section 3.1's own host-by-host
     prose) against 406 (Section 3.1's opening sentence). And 420 against
     Table 2's 422, and 56 against Section 3.1's 54 uninformative.

  2. THE 541. Section 3.2 says 88 of 541 animals whose 18S library yielded a
     resolved taxon, against Section 3.1's 476 animals with an 18S result.

  3. THE 16S DATASET. Section 2.6 and Table 2 say 853 of 907 libraries
     resolve to an animal and 712 of 907 carry 1,000 reads.
     cohorts_and_denominators.tsv says 413 animals from 602 of 853 and
     flags the two as unreconciled.

Usage:  python3 -I check_denoms.py [project_root]
"""

import csv
import os
import re
import sys
from collections import Counter, defaultdict

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())


def say(s=""):
    print(s)


def head(t):
    say()
    say("=" * 76)
    say(t)
    say("=" * 76)


def rel(p):
    try:
        return os.path.relpath(p, ROOT)
    except ValueError:
        return p


def find(pattern, exts=None, maxhits=200):
    rx = re.compile(pattern, re.I)
    out = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        if "/_archive" in dirpath or "/reads_" in dirpath:
            continue
        for fn in filenames:
            if exts and not fn.lower().endswith(tuple(exts)):
                continue
            p = os.path.join(dirpath, fn)
            if rx.search(p):
                out.append(p)
                if len(out) >= maxhits:
                    return sorted(out)
    return sorted(out)


def load(path, cap=500000):
    with open(path, errors="replace") as fh:
        first = fh.readline()
    d = "\t" if first.count("\t") >= first.count(",") else ","
    with open(path, errors="replace") as fh:
        rows = []
        for i, r in enumerate(csv.DictReader(fh, delimiter=d)):
            if i >= cap:
                break
            rows.append(r)
    return rows


def num(v):
    try:
        return float(str(v).strip() or 0)
    except ValueError:
        return 0.0


# ---------------------------------------------------- library -> animal map
head("0.  LIBRARY TO ANIMAL MAP")
say("Everything below needs this. Looking for a table carrying both.")
say()
amap = {}
src = None
for p in find(r"analysis_table_16S|specimen_taxon_table|master_host|crosswalk",
              exts=[".tsv", ".csv"]):
    try:
        rows = load(p, 5000)
    except OSError:
        continue
    if not rows:
        continue
    cols = list(rows[0].keys())
    lib = next((c for c in cols if re.search(r"sample_id|library|specimen_id", c, re.I)), None)
    ani = next((c for c in cols if re.search(r"^animal$|animal_id|barcode", c, re.I)), None)
    if lib and ani:
        say("   %s" % rel(p))
        say("      library col %r   animal col %r   %d rows" % (lib, ani, len(rows)))
        if not amap:
            rows = load(p)
            for r in rows:
                k = (r.get(lib) or "").strip()
                v = (r.get(ani) or "").strip()
                if k and v:
                    amap[k] = v
            src = p
if amap:
    say()
    say("   using %s : %d libraries -> %d animals"
        % (rel(src), len(amap), len(set(amap.values()))))
else:
    say("   NONE FOUND. Counts below will be at library level only.")


def to_animal(lib):
    if lib in amap:
        return amap[lib]
    # fall back: strip the .batch suffix and any -NN subsample tail
    base = lib.split(".")[0]
    return re.sub(r"-\d+$", "", base)


# ------------------------------------------------------------------ 1
head("1.  TICK GENUS: 404 / 406, 420 / 422, 54 / 56")

tg = find(r"tick_genus_per_sample_corrected", exts=[".tsv"])
tg += [p for p in find(r"tick_genus_per_sample", exts=[".tsv"]) if p not in tg]
for p in tg[:2]:
    rows = load(p)
    cols = list(rows[0].keys())
    say()
    say("   --- %s  (%d libraries) ---" % (rel(p), len(rows)))
    say("   columns: %s" % ", ".join(cols))
    gcols = [c for c in cols
             if c not in ("sample_id", "tick_reads", "dominant")
             and not re.search(r"sample|read|dominant", c, re.I)]
    say("   genus columns: %s" % gcols)

    RULES = {
        "reads >= 50 AND >= 1% of tick reads  (Section 2.1 rule)":
            lambda r, g: num(r.get(g)) >= 50 and num(r.get("tick_reads")) > 0
            and num(r.get(g)) / num(r.get("tick_reads")) >= 0.01,
        "reads >= 50 only":
            lambda r, g: num(r.get(g)) >= 50,
        "reads > 0 (any)":
            lambda r, g: num(r.get(g)) > 0,
        "dominant column equals the genus":
            lambda r, g: (r.get("dominant") or "").strip().lower() == g.lower(),
    }
    for label, test in RULES.items():
        say()
        say("   RULE: %s" % label)
        peranimal = defaultdict(set)
        for r in rows:
            a = to_animal((r.get("sample_id") or "").strip())
            for g in gcols:
                if test(r, g):
                    peranimal[a].add(g)
        allan = {to_animal((r.get("sample_id") or "").strip()) for r in rows}
        for g in gcols:
            n = sum(1 for a, s in peranimal.items() if g in s)
            say("      animals carrying %-16s %5d" % (g, n))
        both = sum(1 for a, s in peranimal.items()
                   if "Ixodes" in s and "Dermacentor" in s)
        withcall = sum(1 for a, s in peranimal.items() if s)
        say("      animals carrying BOTH Ixodes and Dermacentor  %5d" % both)
        say("      animals with ANY genus call                   %5d" % withcall)
        say("      animals with NO genus call                    %5d"
            % (len(allan) - withcall))
        say("      animals in this table                         %5d" % len(allan))

say()
say("   Target: Table 1's column sums to 404 Ixodes and 20 Dermacentor,")
say("   overlapping on 4. Section 3.1 opens with 406 and 54 uninformative;")
say("   Table 2 says 422 of 476 clear the floors. Whichever rule reproduces")
say("   one of those sets is the rule the manuscript should state.")

# ------------------------------------------------------------------ 2
head("2.  THE 541: ANIMALS WHOSE 18S YIELDED A RESOLVED TAXON")

comm = find(r"community_18S_corrected|counts_taxon_18S_BLASTcorrected",
            exts=[".tsv"])
for p in comm[:1]:
    say("   %s" % rel(p))
    with open(p, errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        samples = hdr[1:]
        say("   %d taxa rows x %d sample columns" % (0, len(samples)))
        nonzero = [0] * len(samples)
        resolved = [0] * len(samples)
        pir = defaultdict(lambda: [0] * len(samples))
        ntax = 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(hdr):
                continue
            ntax += 1
            taxon = parts[0]
            deep = len([x for x in taxon.split(";")
                        if x.strip() and not x.strip().endswith("__")
                        and "unassigned" not in x.lower()])
            vals = [num(v) for v in parts[1:]]
            for i, v in enumerate(vals):
                if v > 0:
                    nonzero[i] += 1
                    if deep >= 2:
                        resolved[i] += 1
            if re.search(r"Babesia", taxon, re.I):
                key = "B. odocoilei" if re.search(r"odocoilei", taxon, re.I) else (
                    "B. microti" if re.search(r"microti", taxon, re.I) else "Babesia other")
                for i, v in enumerate(vals):
                    if v > 0:
                        pir[key][i] += 1
        say("   %d taxa" % ntax)

        def animals(mask):
            return len({to_animal(samples[i]) for i, m in enumerate(mask) if m})

        say()
        say("   libraries with any non-zero taxon        %5d of %d"
            % (sum(1 for x in nonzero if x), len(samples)))
        say("   ANIMALS  with any non-zero taxon         %5d" % animals(nonzero))
        say("   libraries with a taxon below domain      %5d"
            % sum(1 for x in resolved if x))
        say("   ANIMALS  with a taxon below domain       %5d" % animals(resolved))
        say()
        for k in sorted(pir):
            say("   ANIMALS carrying %-16s %5d" % (k, animals(pir[k])))
        say()
        say("   Section 3.2 says 88 of 541 carry B. odocoilei and 18 of the")
        say("   same 541 carry B. microti. Section 3.1 says 476 have an 18S")
        say("   result. Whichever count above is 541 names the rule.")

# ------------------------------------------------------------------ 3
head("3.  THE 16S DATASET: 712 OF 907 AGAINST 602 OF 853")

at = find(r"analysis_table_16S", exts=[".tsv"])
for p in at[:1]:
    rows = load(p)
    cols = list(rows[0].keys())
    say("   %s  (%d rows)" % (rel(p), len(rows)))
    say("   columns: %s" % ", ".join(cols[:16]))
    lib = next((c for c in cols if re.search(r"sample_id|library", c, re.I)), cols[0])
    ani = next((c for c in cols if re.search(r"^animal$|animal_id|barcode", c, re.I)), None)
    rc = [c for c in cols if re.search(r"total_reads|reads_total|^reads$|n_reads|depth", c, re.I)]
    say("   library col %r   animal col %r   read cols %s" % (lib, ani, rc))
    say()
    say("   libraries in table                      %5d" % len(rows))
    if ani:
        res = [r for r in rows if (r.get(ani) or "").strip()]
        say("   libraries resolving to an animal        %5d" % len(res))
        say("   distinct animals                        %5d"
            % len({(r.get(ani) or "").strip() for r in res}))
    for c in rc:
        for floor in (1000,):
            keep = [r for r in rows if num(r.get(c)) >= floor]
            say("   libraries with %s >= %d           %5d" % (c, floor, len(keep)))
            if ani:
                say("      distinct animals among them          %5d"
                    % len({(r.get(ani) or "").strip() for r in keep
                           if (r.get(ani) or "").strip()}))

say()
say("=" * 76)
say("done. nothing was written.")
say("=" * 76)
