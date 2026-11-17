#!/usr/bin/env python3
"""
check_merge.py -- READ ONLY. Writes nothing, changes nothing.

Diagnoses the 18S resolution problem.

Hypothesis: 515F/926R gives ~411 bp in bacteria but a longer product on
eukaryotic 18S. Truncating at 220 forward and 200 reverse leaves 420 bp of
sequence. If the 18S amplicon is longer than that, the reads cannot overlap
and DADA2's merge step throws them away without erroring. The survivors
would be whichever eukaryotes happen to have a short enough amplicon, which
is a biased remnant, not a sample of the community.

Four things, read straight out of the QIIME 2 artifacts (a .qza is a zip):

  1. Denoising stats per marker: input -> filtered -> denoised -> merged ->
     non-chimeric. If 16S merges at 80% and 18S at 5%, that is the answer.
  2. ASV length distribution per marker, from rep-seqs.
  3. How deep 18S assignments actually go: species, genus, family, higher.
  4. The Ixodes ASVs specifically: length, assignment, read count.

Usage:  python3 -I check_merge.py [project_root]
"""

import os
import re
import sys
import zipfile
from collections import Counter, defaultdict

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
MARKERS = ("16S", "18S", "cytb")


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


def qza_member(path, pattern):
    """Read the first file inside a .qza whose name matches. A qza is a zip."""
    try:
        z = zipfile.ZipFile(path)
    except (OSError, zipfile.BadZipFile):
        return None
    rx = re.compile(pattern)
    for n in z.namelist():
        if rx.search(n):
            try:
                return z.read(n).decode("utf8", "replace")
            except OSError:
                return None
    return None


def find(pattern):
    rx = re.compile(pattern)
    out = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            if rx.search(p):
                out.append(p)
    return sorted(out)


# ------------------------------------------------------------------ 1
head("1.  DENOISING STATS: WHERE DO THE READS GO?")
say("The merged column is the one that matters. A marker whose amplicon is")
say("too long to overlap loses nearly everything at that step.")

for marker in MARKERS:
    cands = [p for p in find(r"denoising[-_]stats.*\.qza$")
             if "/%s/" % marker in p]
    if not cands:
        say()
        say("   %-5s no denoising-stats.qza found" % marker)
        continue
    path = sorted(cands, key=os.path.getmtime)[-1]
    txt = qza_member(path, r"/data/.*\.tsv$")
    if not txt:
        say()
        say("   %-5s could not read stats from %s" % (marker, rel(path)))
        continue

    lines = [L for L in txt.splitlines() if L.strip()]
    hdr = lines[0].split("\t")
    # QIIME writes a second header row of column types; skip it
    body = [L.split("\t") for L in lines[1:] if not L.startswith("#")]
    body = [r for r in body if len(r) == len(hdr)
            and not r[1].strip().lower().startswith("numeric")]

    cols = {h.strip().lower(): i for i, h in enumerate(hdr)}

    def total(*names):
        for n in names:
            for key, i in cols.items():
                if n in key and "percentage" not in key:
                    s = 0.0
                    ok = False
                    for r in body:
                        try:
                            s += float(r[i])
                            ok = True
                        except (ValueError, IndexError):
                            pass
                    if ok:
                        return s
        return None

    inp = total("input")
    filt = total("filtered")
    den = total("denoised")
    mer = total("merged")
    nonchim = total("non-chimeric", "nonchimeric")

    say()
    say("   --- %s   (%d samples)   %s" % (marker, len(body), rel(path)))
    prev = inp
    for name, val in (("input", inp), ("filtered", filt), ("denoised", den),
                      ("merged", mer), ("non-chimeric", nonchim)):
        if val is None:
            say("      %-14s (column not present)" % name)
            continue
        pct_in = 100.0 * val / inp if inp else 0.0
        step = (" step %5.1f%%" % (100.0 * val / prev)) if prev else ""
        say("      %-14s %14.0f   %5.1f%% of input%s" % (name, val, pct_in, step))
        prev = val if val else prev
    if inp and mer is not None and den:
        say()
        say("      >>> MERGE SURVIVAL: %.1f%% of denoised reads merged"
            % (100.0 * mer / den))
        if 100.0 * mer / den < 50:
            say("      >>> Most reads died at the merge. The amplicon is too")
            say("      >>> long to overlap at these truncation lengths.")

# ------------------------------------------------------------------ 2
head("2.  ASV LENGTH DISTRIBUTION PER MARKER")
say("A merged amplicon has a characteristic length. A marker that barely")
say("merged will show a narrow spike at the shortest lengths that could")
say("overlap, instead of the real amplicon length.")

for marker in MARKERS:
    cands = [p for p in find(r"rep[-_]seqs.*\.qza$") if "/%s/" % marker in p]
    if not cands:
        say()
        say("   %-5s no rep-seqs.qza found" % marker)
        continue
    path = sorted(cands, key=os.path.getmtime)[-1]
    fa = qza_member(path, r"/data/.*\.fasta$")
    if not fa:
        say()
        say("   %-5s could not read sequences from %s" % (marker, rel(path)))
        continue
    lens, cur = [], 0
    for L in fa.splitlines():
        if L.startswith(">"):
            if cur:
                lens.append(cur)
            cur = 0
        else:
            cur += len(L.strip())
    if cur:
        lens.append(cur)
    if not lens:
        continue
    lens.sort()
    say()
    say("   --- %s   %d ASVs   %s" % (marker, len(lens), rel(path)))
    say("      min %d   p05 %d   median %d   p95 %d   max %d"
        % (lens[0], lens[int(.05 * (len(lens) - 1))], lens[len(lens) // 2],
           lens[int(.95 * (len(lens) - 1))], lens[-1]))
    buckets = Counter((L // 20) * 20 for L in lens)
    say("      length histogram (20 bp bins, top 8):")
    for b, n in buckets.most_common(8):
        bar = "#" * min(44, int(44.0 * n / max(buckets.values())))
        say("        %4d-%4d  %7d  %s" % (b, b + 19, n, bar))

# ------------------------------------------------------------------ 3
head("3.  HOW DEEP DO 18S ASSIGNMENTS GO?")

tax = None
for cand in (os.path.join(ROOT, "qiime2/18S/all/taxonomy-exported-temp/taxonomy.tsv"),):
    if os.path.exists(cand):
        tax = cand
if tax is None:
    hits = [p for p in find(r"/18S/.*taxonomy\.tsv$")]
    tax = hits[0] if hits else None

if tax is None:
    say("   no 18S taxonomy.tsv found")
else:
    say("   %s" % rel(tax))
    depth = Counter()
    example = {}
    rows = 0
    with open(tax, errors="replace") as fh:
        for i, L in enumerate(fh):
            if i == 0:
                continue
            parts = L.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rows += 1
            lin = [x for x in parts[1].split(";") if x.strip()
                   and not re.match(r"^[a-z]?_*$", x.strip())
                   and "unassigned" not in x.lower()
                   and not x.strip().endswith("__")]
            depth[len(lin)] += 1
            if lin and len(lin) not in example:
                example[len(lin)] = lin[-1].strip()
    say("   %d ASVs" % rows)
    say()
    say("   PR2 lineages run domain;supergroup;division;subdivision;class;")
    say("   order;family;genus;species, so 8 levels is a genus call and 9 is")
    say("   a species call. The example shows what the deepest name looks like.")
    say()
    say("   %-8s %9s %8s   %s" % ("levels", "ASVs", "pct", "example deepest name"))
    for d in sorted(depth):
        say("   %-8d %9d %7.1f%%   %s"
            % (d, depth[d], 100.0 * depth[d] / rows, example.get(d, "")[:40]))

# ------------------------------------------------------------------ 4
head("4.  THE IXODES ASVs SPECIFICALLY")
say("Length, assignment and whether they are the short survivors of a")
say("failed merge.")

if tax:
    seqlen = {}
    cands = [p for p in find(r"rep[-_]seqs.*\.qza$") if "/18S/" in p]
    if cands:
        fa = qza_member(sorted(cands, key=os.path.getmtime)[-1],
                        r"/data/.*\.fasta$")
        if fa:
            fid, cur = None, []
            for L in fa.splitlines():
                if L.startswith(">"):
                    if fid:
                        seqlen[fid] = len("".join(cur))
                    fid = L[1:].split()[0]
                    cur = []
                else:
                    cur.append(L.strip())
            if fid:
                seqlen[fid] = len("".join(cur))

    hits = []
    with open(tax, errors="replace") as fh:
        for i, L in enumerate(fh):
            if i == 0:
                continue
            parts = L.rstrip("\n").split("\t")
            if len(parts) >= 2 and re.search(r"Ixodes|Dermacentor|Amblyomma",
                                             parts[1], re.I):
                hits.append(parts)
    say("   %d tick-genus ASVs" % len(hits))
    if hits:
        byname = Counter()
        for p in hits:
            leaf = [x for x in p[1].split(";") if x.strip()][-1].strip()
            byname[leaf] += 1
        say()
        say("   %-46s %6s" % ("deepest assignment", "ASVs"))
        for name, n in byname.most_common(14):
            say("   %-46s %6d" % (name[:46], n))
        if seqlen:
            ls = sorted(seqlen.get(p[0], 0) for p in hits
                        if seqlen.get(p[0]))
            if ls:
                say()
                say("   their lengths: min %d  median %d  max %d  (n=%d)"
                    % (ls[0], ls[len(ls) // 2], ls[-1], len(ls)))

say()
say("=" * 76)
say("done. nothing was written.")
say("=" * 76)
