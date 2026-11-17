#!/usr/bin/env python3
"""
check_three.py -- READ ONLY. Writes nothing, changes nothing.

  A. The cross-contamination rule. Section 2.6 now says only that non-target
     host reads "were assessed". Find the script that did it and print the
     actual threshold, so the sentence can state a rule.

  B. The 330 funnel. Section 2.6 gives 366 -> 338 -> 330 and then says five
     ked tubes were excluded BEFORE those floors, also ending at 330. Both
     cannot be right. Also the 98 against 114 animals with several libraries.

  C. The terminal base. If the 14 bp core CCGYCAATTYMTTT was the reverse
     primer, every R2 read keeps RAGTTT, which becomes AAACTY at the 3' end
     of the merged amplicon, and Y is C or T. That would make the
     "denoising artifact" pairs primer degeneracy instead. Also checks the
     16S merged length range, which Section 2.6 gives as 220 to 408 bp.

Usage:  python3 -I check_three.py [project_root]
"""

import os
import re
import sys
import zipfile
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


def find(pattern, exts=None, maxhits=600):
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
                    return sorted(out)
    return sorted(out)


def lines(path, limit=200000):
    try:
        with open(path, "r", errors="replace") as fh:
            for i, L in enumerate(fh):
                if i >= limit:
                    return
                yield i + 1, L
    except OSError:
        return


# ------------------------------------------------------------------ A
head("A.  THE CROSS-CONTAMINATION RULE")
say("Looking for where non-target host reads were judged, and on what number.")
say()

scripts = find(r"\.(py|sh|R|r)$", maxhits=900)
scripts = [p for p in scripts if "/reads_" not in p]

say("  -- scripts whose NAME suggests this step --")
named = [p for p in scripts
         if re.search(r"contam|crosstalk|cross_|hop|nontarget|non_target|host", os.path.basename(p), re.I)]
for p in named:
    say("     %s" % rel(p))
if not named:
    say("     (none by name)")

say()
say("  -- lines that set a threshold on non-target host reads --")
PAT = (r"cross[-_ ]?contam|non[-_ ]?target|hopp?ed|hopping|"
       r"other_host|offtarget|off[-_ ]target")
hits = 0
for p in scripts:
    for n, L in lines(p, 40000):
        if re.search(PAT, L, re.I):
            # only show lines that also carry a number or comparison
            if re.search(r"[<>]=?|\b0\.\d+|\b\d+\s*%|threshold|cutoff|min_|frac", L):
                say("     %s:%d: %s" % (rel(p), n, L.strip()[:150]))
                hits += 1
                if hits >= 22:
                    break
    if hits >= 22:
        say("     (stopping at 22)")
        break
if not hits:
    say("     no thresholded line found; showing every mention instead")
    h2 = 0
    for p in scripts:
        for n, L in lines(p, 40000):
            if re.search(PAT, L, re.I):
                say("     %s:%d: %s" % (rel(p), n, L.strip()[:150]))
                h2 += 1
                if h2 >= 18:
                    break
        if h2 >= 18:
            break

say()
say("  -- any results file recording the outcome --")
for p in find(r"cross[-_]?contam|nontarget|non_target", exts=[".tsv", ".csv", ".txt", ".md"], maxhits=20):
    say("     %10d  %s" % (os.path.getsize(p), rel(p)))

# ------------------------------------------------------------------ B
head("B.  THE 330 FUNNEL AND THE 98 AGAINST 114")

fun = find(r"check_cohort_funnel", maxhits=10)
say("  funnel script: %s" % ([rel(p) for p in fun] or "not found"))
for p in fun[:1]:
    say()
    say("  -- lines that count or filter --")
    for n, L in lines(p, 40000):
        if re.search(r"\b(366|338|330|114|98|1000|1,000)\b|len\(|nunique|drop|exclude|ked|floor", L):
            s = L.rstrip()
            if s.strip().startswith("#") and len(s) > 120:
                continue
            say("     %4d: %s" % (n, s[:150]))

say()
say("  -- where 366 / 338 / 330 / 114 / 98 appear anywhere in the project --")
for num in ("366", "338", "330", "114", "98"):
    found = []
    for p in find(r"\.(py|md|txt|log|tsv|csv|yml)$", maxhits=900):
        if "/reads_" in p or "/_archive/" in p:
            continue
        for n, L in lines(p, 20000):
            if re.search(r"\b%s\b" % num, L) and re.search(r"animal|deer|cohort|librar|deer|floor", L, re.I):
                found.append("%s:%d: %s" % (rel(p), n, L.strip()[:110]))
                break
        if len(found) >= 4:
            break
    say("     %-5s %s" % (num, found[0] if found else "(no contextual hit)"))
    for f in found[1:4]:
        say("           %s" % f)

say()
say("  -- the ked tubes: are the five excluded before or after the floors? --")
for p in find(r"\.(py|md|txt)$", maxhits=900):
    if "/reads_" in p:
        continue
    for n, L in lines(p, 40000):
        if re.search(r"\bked", L, re.I) and re.search(r"exclud|drop|remov|filter|before|after", L, re.I):
            say("     %s:%d: %s" % (rel(p), n, L.strip()[:150]))

# ------------------------------------------------------------------ C
head("C.  THE TERMINAL BASE, AND THE 16S LENGTH RANGE")

cand = [p for p in find(r"rep[-_]seqs.*\.qza$", maxhits=60)
        if "/16S/" in p and "dada2" in p and "aligned" not in p]
if not cand:
    cand = [p for p in find(r"rep[-_]seqs.*\.qza$", maxhits=60) if "/16S/" in p]
say("  16S rep-seqs candidates:")
for p in cand:
    say("     %10d  %s" % (os.path.getsize(p), rel(p)))

if cand:
    path = sorted(cand, key=os.path.getsize, reverse=True)[0]
    say()
    say("  using %s" % rel(path))
    fa = None
    try:
        z = zipfile.ZipFile(path)
        for n in z.namelist():
            if re.search(r"/data/.*\.fasta$", n):
                fa = z.read(n).decode("utf8", "replace")
                break
        z.close()
    except (OSError, zipfile.BadZipFile) as e:
        say("  could not open: %s" % e)

    if fa:
        seqs = {}
        sid, buf = None, []
        for L in fa.splitlines():
            if L.startswith(">"):
                if sid:
                    seqs[sid] = "".join(buf)
                sid = L[1:].split()[0]
                buf = []
            else:
                buf.append(L.strip())
        if sid:
            seqs[sid] = "".join(buf)
        say("  %d ASVs" % len(seqs))

        lens = sorted(len(s) for s in seqs.values())
        say()
        say("  LENGTH: min %d  p05 %d  median %d  p95 %d  max %d"
            % (lens[0], lens[int(.05 * (len(lens) - 1))], lens[len(lens) // 2],
               lens[int(.95 * (len(lens) - 1))], lens[-1]))
        say("  Section 2.6 says the 16S merged length varied from 220 to 408 bp.")

        say()
        say("  LAST SIX BASES  (expect AAACTC / AAACTT if the 14 bp core was used)")
        tail = Counter(s[-6:] for s in seqs.values() if len(s) >= 6)
        tot = sum(tail.values())
        for t, n in tail.most_common(10):
            say("     %-10s %7d  %5.1f%%" % (t, n, 100.0 * n / tot))

        say()
        say("  LAST BASE ALONE")
        last = Counter(s[-1:] for s in seqs.values() if s)
        for t, n in last.most_common(6):
            say("     %-4s %7d  %5.1f%%" % (t, n, 100.0 * n / len(seqs)))

        say()
        say("  PAIRS DIFFERING ONLY AT THE TERMINAL BASE")
        stem = defaultdict(list)
        for sid, s in seqs.items():
            if s:
                stem[s[:-1]].append((sid, s[-1]))
        multi = {k: v for k, v in stem.items() if len(v) > 1}
        say("     %d stems carry more than one terminal base" % len(multi))
        say("     covering %d ASVs of %d (%.1f%%)"
            % (sum(len(v) for v in multi.values()), len(seqs),
               100.0 * sum(len(v) for v in multi.values()) / len(seqs)))
        combos = Counter(tuple(sorted(b for _, b in v)) for v in multi.values())
        say()
        say("     which terminal bases pair up:")
        for c, n in combos.most_common(8):
            say("        %-14s %6d stems" % ("/".join(c), n))
        say()
        say("     If C/T dominates, these are the primer's Y position and not")
        say("     a denoising artifact. Section 2.6's explanation then changes.")

say()
say("=" * 76)
say("done. nothing was written.")
say("=" * 76)
