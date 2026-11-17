#!/usr/bin/env python3
"""
read_primers.py -- READ ONLY. Writes nothing, changes nothing.

Looks at the actual read sequences instead of arguing from manifests.

  PART A. The manifests, by CONTENT not by size. Same files or not, with the
          filenames printed so you can see them.
  PART B. The most common first 25 bases of real reads, tabulated. This assumes
          nothing about which primers to expect: whatever is actually on the
          front of the reads shows up at the top of the list.
  PART C. Then, and only then, the candidate primers matched against those
          same reads, IUPAC-aware, allowing mismatches.
  PART D. cutadapt's own "Trimmed: N times" totals per primer per marker.
          This is cutadapt telling you how many reads carried each primer.

Usage:  python3 -I read_primers.py [project_root]
        python3 -I read_primers.py . 50000     # deeper read sample
"""

import gzip
import os
import re
import sys
from collections import Counter, defaultdict

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
NREADS = int(sys.argv[2]) if len(sys.argv) > 2 else 20000

CANDIDATES = {
    "515F  (manuscript, fwd)":      "GTGYCAGCMGCCGCGGTAA",
    "926R  (manuscript, rev)":      "CCGYCAATTYMTTTRAGTTT",
    "926R 14bp core":               "CCGYCAATTYMTTT",
    "18S-V4 fwd (Lawrence)":        "CCAGCASCYGCGGTAATTC",
    "18S-V4 rev (Lawrence)":        "ACTTTCGTTCTTGATYRA",
    "1391F (the broken V9 ref)":    "GTACACACCGCCCGTC",
    "EukBr (the broken V9 ref)":    "TGATCCTTCTGCAGGTTCACCTAC",
    "L14841 cytb fwd":              "CCATCCAACATCTCAGCATGATGAAA",
    "H15149 cytb rev":              "CCCCTCAGAATGATATTTGTCCTCA",
}

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def expand(p):
    return [set(IUPAC.get(c, c)) for c in p.upper()]


def mismatches(seq, pat):
    """How many positions of seq fail the IUPAC pattern. Stops early."""
    n = 0
    for base, allowed in zip(seq, pat):
        if base not in allowed:
            n += 1
            if n > 4:
                return n
    return n


def say(s=""):
    print(s)


def head(t):
    say()
    say("=" * 76)
    say(t)
    say("=" * 76)


def walk(pat, maxhits=4000):
    rx = re.compile(pat, re.I)
    out = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            if rx.search(p):
                out.append(p)
                if len(out) >= maxhits:
                    return out
    return out


def rel(p):
    try:
        return os.path.relpath(p, ROOT)
    except ValueError:
        return p


# ------------------------------------------------------------------- PART A
head("A.  THE MANIFESTS, COMPARED BY CONTENT")

m16 = [p for p in walk(r"manifest_16S\.tsv$")]
m18 = [p for p in walk(r"manifest_18S\.tsv$")]


def read_manifest(p):
    rows = []
    with open(p, errors="replace") as fh:
        for L in fh:
            rows.append(L.rstrip("\n"))
    return rows


if m16 and m18:
    a, b = read_manifest(m16[0]), read_manifest(m18[0])
    say("   16S: %s  (%d lines)" % (rel(m16[0]), len(a)))
    say("   18S: %s  (%d lines)" % (rel(m18[0]), len(b)))
    say()
    say("   identical line-for-line: %s" % ("YES" if a == b else "NO"))
    if a != b:
        diffs = [(i, x, y) for i, (x, y) in enumerate(zip(a, b), 1) if x != y]
        say("   %d lines differ. First few:" % len(diffs))
        for i, x, y in diffs[:6]:
            say("      line %d" % i)
            say("        16S: %s" % x[:150])
            say("        18S: %s" % y[:150])
    say()
    say("   first 4 lines of the 16S manifest:")
    for L in a[:4]:
        say("      %s" % L[:170])
    say()
    say("   first 4 lines of the 18S manifest:")
    for L in b[:4]:
        say("      %s" % L[:170])
else:
    say("   manifests not found at those names. What exists:")
    for p in walk(r"manifest.*\.(tsv|csv)$", 40):
        say("      %s" % rel(p))

# ------------------------------------------------------------------- PART B
head("B.  WHAT IS ACTUALLY ON THE FRONT OF THE READS")
say("No assumption about which primers to expect. The commonest prefix at")
say("three widths, over %d reads per file. A primer is a CONSTANT prefix:" % NREADS)
say("it dominates at its own length, then dissolves into hundreds of distinct")
say("prefixes once you read past it into the variable region.")

fq = walk(r"\.f(ast)?q\.gz$", 4000)
raw = [p for p in fq if "/reads_original/" in p]
staged = [p for p in fq if "/reads_staged/" in p]
other = [p for p in fq if p not in raw and p not in staged]

say()
say("   fastq.gz found: %d total  (%d under reads_original, %d under"
    " reads_staged, %d elsewhere)" % (len(fq), len(raw), len(staged), len(other)))


def pick(pool, n=3):
    r1 = sorted(p for p in pool if re.search(r"_R1[_.]", p))
    return r1[:n] if r1 else sorted(pool)[:n]


groups = []
if raw:
    groups.append(("reads_original  (RAW, primers should still be present)",
                   pick(raw)))
if staged:
    groups.append(("reads_staged    (may already be trimmed)", pick(staged)))
if other:
    groups.append(("elsewhere", pick(other)))


def heads(path, n, width=40):
    """First `width` bases of up to n reads."""
    out = []
    op = gzip.open if path.endswith(".gz") else open
    try:
        with op(path, "rt", errors="replace") as fh:
            for i, L in enumerate(fh):
                if i % 4 == 1:
                    out.append(L.strip()[:width].upper())
                    if len(out) >= n:
                        break
    except OSError as e:
        say("      (could not read: %s)" % e)
    return out


samples = {}
for label, paths in groups:
    say()
    say("   --- %s ---" % label)
    for p in paths:
        hs = heads(p, NREADS)
        samples[p] = hs
        if not hs:
            continue
        say("   %s   (%d reads sampled)" % (rel(p), len(hs)))
        # Tabulate at several widths. A primer is a constant prefix, so it
        # dominates at its own length and dissolves once you read past it
        # into the variable region. That contrast is the signal.
        for width in (13, 19, 24):
            c = Counter(h[:width] for h in hs if len(h) >= width)
            if not c:
                continue
            top, n = c.most_common(1)[0]
            say("        first %2d bases: %-26s %6d  %5.1f%%  (%d distinct)"
                % (width, top, n, 100.0 * n / len(hs), len(c)))
            for prefix, k in c.most_common(4)[1:]:
                say("        %15s  %-26s %6d  %5.1f%%"
                    % ("", prefix, k, 100.0 * k / len(hs)))
        say()

# ------------------------------------------------------------------- PART C
head("C.  THE CANDIDATE PRIMERS, MATCHED AGAINST THOSE SAME READS")
say("IUPAC-aware, at read start, allowing up to 2 mismatches.")

pats = {k: (v, expand(v)) for k, v in CANDIDATES.items()}
for p, hs in samples.items():
    if not hs:
        continue
    say()
    say("   %s" % rel(p))
    anyhit = False
    for name, (seq, pat) in pats.items():
        hit = sum(1 for h in hs if len(h) >= len(seq)
                  and mismatches(h[:len(seq)], pat) <= 2)
        if hit:
            anyhit = True
            say("      %-30s %-26s %6d  %5.1f%%"
                % (name, seq, hit, 100.0 * hit / len(hs)))
    if not anyhit:
        say("      no candidate primer matches at read start.")
        say("      (reads are probably already primer-trimmed; read Part B"
            " and Part D instead)")

# ------------------------------------------------------------------- PART D
head("D.  CUTADAPT'S OWN TRIM COUNTS, PER PRIMER PER MARKER")
say("cutadapt prints 'Sequence: X; ...; Trimmed: N times' for each primer it")
say("was given. This is the tool reporting how many reads carried each one.")

logs = walk(r"cutadapt.*\.log$", 400)
by = defaultdict(lambda: defaultdict(int))
seen = defaultdict(set)
rx = re.compile(
    r"Sequence:\s*([ACGTRYSWKMBDHVNacgtryswkmbdhvn]{10,45})\s*;"
    r".*?Trimmed:\s*([0-9]+)\s*times", re.S)

for p in logs:
    marker = "other"
    for m in ("16S", "18S", "cytb", "ITS"):
        if "/%s/" % m in p or m.lower() in os.path.basename(p).lower():
            marker = m
            break
    try:
        txt = open(p, errors="replace").read(8_000_000)
    except OSError:
        continue
    seen[marker].add(p)
    for m in rx.finditer(txt):
        by[marker][m.group(1).upper()] += int(m.group(2))

known = {v.upper(): k for k, v in CANDIDATES.items()}
for marker in ("16S", "18S", "cytb", "ITS", "other"):
    if marker not in by:
        continue
    say()
    say("   --- %s   (%d log files) ---" % (marker, len(seen[marker])))
    tot = sum(by[marker].values())
    for seq, n in sorted(by[marker].items(), key=lambda kv: -kv[1])[:8]:
        tag = known.get(seq, "")
        if not tag:
            for k, v in CANDIDATES.items():
                if seq.startswith(v.upper()[:14]) or v.upper().startswith(seq[:14]):
                    tag = k + " (partial)"
                    break
        say("      %-28s trimmed %12d times  %5.1f%%   %s"
            % (seq, n, 100.0 * n / tot if tot else 0, tag))

say()
say("=" * 76)
say("HOW TO READ THIS")
say("  Part A: if the two manifests are identical line for line, 16S and 18S")
say("          are the same libraries and there is no separate 18S amplicon.")
say("  Part B: the truth, with nothing assumed. Whatever primer is on the")
say("          reads is the commonest prefix in a RAW file.")
say("  Part D: whichever primer cutadapt trimmed from millions of reads is")
say("          the primer that was used. One trimmed from zero was not.")
say("=" * 76)
say("done. nothing was written.")
