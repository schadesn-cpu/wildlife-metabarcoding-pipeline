#!/usr/bin/env python3
"""
check_18S_denoms.py -- READ ONLY. Writes nothing, changes nothing.

Three questions, one pass.

  1. The 18S denominators. Section 3.1 says 476 animals have an 18S result.
     Section 3.2 says 541 animals' 18S libraries yielded a resolved taxon.
     Both cannot be "animals with an 18S result". This counts the funnel so
     the manuscript can name each step.

  2. Ixodes: 406 in Section 3.1 against 404 in Table 1.

  3. cutadapt error rate per log file, with the batch each log belongs to.
     Section 2.3 says 10%; the logs carry both 0.1 and 0.2. If 0.2 belongs
     to the batch whose reverse primer differed at 3 nucleotides, the fix is
     one clause.

Usage:  python3 -I check_18S_denoms.py [project_root]
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
    say("=" * 74)
    say(t)
    say("=" * 74)


def find(pattern, exts=None, maxhits=400):
    rx = re.compile(pattern, re.I)
    out = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        for fn in filenames:
            if exts and not fn.lower().endswith(tuple(exts)):
                continue
            p = os.path.join(dirpath, fn)
            if rx.search(p):
                out.append(p)
                if len(out) >= maxhits:
                    return out
    return sorted(out)


def rel(p):
    try:
        return os.path.relpath(p, ROOT)
    except ValueError:
        return p


def sniff(path):
    with open(path, errors="replace") as fh:
        first = fh.readline()
    return "\t" if first.count("\t") >= first.count(",") else ","


def load(path, cap=400000):
    d = sniff(path)
    with open(path, errors="replace") as fh:
        r = csv.DictReader(fh, delimiter=d)
        rows = []
        for i, row in enumerate(r):
            if i >= cap:
                break
            rows.append(row)
    return rows


# ------------------------------------------------------------------ 1 & 2
head("1 & 2.  THE 18S FUNNEL AND THE IXODES COUNT")

cands = [p for p in find(r"tick_genus_per_sample|tick_genus", exts=[".tsv", ".csv"])]
cands += [p for p in find(r"/18S/.*per_sample", exts=[".tsv", ".csv"])
          if p not in cands]
say("candidate tick-genus tables:")
for p in cands:
    say("   %10d  %s" % (os.path.getsize(p), rel(p)))

for p in cands[:3]:
    rows = load(p)
    if not rows:
        continue
    cols = list(rows[0].keys())
    say()
    say("   --- %s  (%d rows) ---" % (rel(p), len(rows)))
    say("   columns: %s" % ", ".join(cols[:14]))

    def col(*names):
        for n in names:
            for c in cols:
                if n in c.lower():
                    return c
        return None

    idc = col("animal", "sample", "barcode", "specimen", "host_id")
    genc = col("genus", "taxon", "tick_genus", "assign")
    if not idc:
        say("   (no id column found)")
        continue
    say("   id column: %s    genus column: %s" % (idc, genc))

    ids = {(r.get(idc) or "").strip() for r in rows if (r.get(idc) or "").strip()}
    say("   distinct ids: %d" % len(ids))

    if genc:
        gen = Counter()
        byid = defaultdict(set)
        for r in rows:
            v = (r.get(genc) or "").strip()
            gen[v or "(blank)"] += 1
            if (r.get(idc) or "").strip():
                byid[(r.get(idc)).strip()].add(v)
        say()
        say("   %-34s %7s rows" % ("genus value", "n"))
        for v, n in gen.most_common(14):
            say("   %-34s %7d" % (v[:34], n))
        say()
        for want in ("Ixodes", "Dermacentor", "Amblyomma"):
            n = len({k for k, vs in byid.items()
                     if any(want.lower() in (x or "").lower() for x in vs)})
            say("   animals with any %-12s assignment: %d" % (want, n))

# the resolved-taxon side
head("1b.  HOW MANY ANIMALS HAVE *ANY* RESOLVED 18S TAXON")
say("This is the 541 in Section 3.2. It should be countable from the per-sample")
say("18S community table rather than the tick-genus table.")
say()
comm = find(r"/18S/.*(community|counts_taxon|relabund|per_sample)",
            exts=[".tsv", ".csv"], maxhits=30)
for p in comm[:6]:
    say("   %10d  %s" % (os.path.getsize(p), rel(p)))
    try:
        with open(p, errors="replace") as fh:
            say("        header: %s" % fh.readline().strip()[:190])
    except OSError:
        pass

# the Babesia odocoilei side, which is where 541 is used
head("1c.  WHERE THE 541 COMES FROM: B. odocoilei AND B. microti BY ANIMAL")
say("Section 3.2: 88 of 541 animals carry B. odocoilei, 18 of the same 541")
say("carry B. microti. Looking for the table those were computed from.")
say()
pir = find(r"piroplasm|babesia|odocoilei", exts=[".tsv", ".csv"], maxhits=30)
for p in pir[:8]:
    say("   %10d  %s" % (os.path.getsize(p), rel(p)))
    try:
        with open(p, errors="replace") as fh:
            say("        header: %s" % fh.readline().strip()[:190])
    except OSError:
        pass

# ------------------------------------------------------------------ 3
head("3.  CUTADAPT ERROR RATE, PER LOG, WITH ITS BATCH")
say("Section 2.3 says 10%. If 0.2 belongs only to the batch whose reverse")
say("primer differed at 3 nucleotides, one clause fixes the sentence.")
say()
logs = find(r"cutadapt.*\.log$", maxhits=200)
say("   %-52s %-10s %s" % ("log", "err rates", "primers / batch hints"))
say("   " + "-" * 104)
for p in logs:
    rates, prim, batches = Counter(), Counter(), Counter()
    try:
        with open(p, errors="replace") as fh:
            for i, L in enumerate(fh):
                if i > 300000:
                    break
                m = re.search(r"--error-rate\s+([0-9.]+)", L)
                if m:
                    rates[m.group(1)] += 1
                m = re.search(r"Sequence:\s*([ACGTRYSWKMBDHVN]{10,30})\s*;", L)
                if m:
                    prim[m.group(1)] += 1
                for m in re.finditer(r"(batch[-_ ]?\d+|Batch[-_ ]?\d+|P\d+-\d+|VDL[-\w]*)", L):
                    batches[m.group(1)] += 1
    except OSError:
        continue
    say("   %-52s %-10s %s" % (rel(p)[-52:], ",".join(sorted(rates)),
                               list(prim)[:2]))
    if batches:
        say("   %-52s %-10s batches: %s" % ("", "", [b for b, _ in batches.most_common(4)]))

say()
say("   If every 0.2 log is the same batch that shows the 14 bp core")
say("   CCGYCAATTYMTTT and no full-length 926R, the manuscript sentence")
say("   becomes: 'at 10% error tolerance, raised to 20% in the batch whose")
say("   reverse primer differed at 3 nucleotides, where only 14 bp of that")
say("   primer were matched.'")

say()
say("=" * 74)
say("done. nothing was written.")
say("=" * 74)
