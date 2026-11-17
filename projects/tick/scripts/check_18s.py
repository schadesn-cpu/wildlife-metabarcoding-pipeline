#!/usr/bin/env python3
"""
check_18S.py -- READ ONLY. Writes nothing, changes nothing.

Settles whether 18S was its own amplification or the eukaryotic product of
515F/926R, and what the 18S classifier was actually trained on.

Five tests, cheapest first:

  1. Do manifest_16S.tsv and manifest_18S.tsv point at the SAME fastq files?
     If yes, 16S and 18S came out of one amplification and one set of
     libraries, and the manuscript is right. If they point at different
     files, 18S was amplified separately and Lawrence is right.
  2. What sequencing projects exist? A separate 18S amplification would have
     its own project names, the way cytb has "-cytoB" ones.
  3. What primer sequences did cutadapt actually trim, per marker? cutadapt
     prints "Sequence: ...; Type: ...; Trimmed: N times" for each one.
  4. What error rate, per marker? The manuscript says 10%; the 18S log said 0.2.
  5. What region was the PR2 reference extracted to before training?

Usage:  python3 -I check_18S.py [project_root]
"""

import os
import re
import subprocess
import sys
from collections import Counter

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())

LAWRENCE_F = "CCAGCASCYGCGGTAATTC"
LAWRENCE_R = "ACTTTCGTTCTTGATYRA"
F515 = "GTGYCAGCMGCCGCGGTAA"
R926 = "CCGYCAATTYMTTTRAGTTT"


def say(s=""):
    print(s)


def head(t):
    say()
    say("=" * 74)
    say(t)
    say("=" * 74)


def find(pattern, subdir="", maxhits=200):
    hits = []
    base = os.path.join(ROOT, subdir) if subdir else ROOT
    rx = re.compile(pattern, re.I)
    for dirpath, dirnames, filenames in os.walk(base):
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        for fn in filenames:
            if rx.search(fn):
                hits.append(os.path.join(dirpath, fn))
                if len(hits) >= maxhits:
                    return hits
    return hits


def rel(p):
    try:
        return os.path.relpath(p, ROOT)
    except ValueError:
        return p


def lines(path, limit=400000):
    try:
        with open(path, "r", errors="replace") as fh:
            for i, L in enumerate(fh):
                if i >= limit:
                    return
                yield L
    except OSError:
        return


# ---------------------------------------------------------------- TEST 1
head("TEST 1.  DO 16S AND 18S READ THE SAME FASTQ FILES?")
say("This is the whole question. One amplification means one set of libraries.")
say()

m16 = find(r"manifest.*16S|16S.*manifest")
m18 = find(r"manifest.*18S|18S.*manifest")
say("16S manifests found:")
for p in m16:
    say("   %10d  %s" % (os.path.getsize(p), rel(p)))
say("18S manifests found:")
for p in m18:
    say("   %10d  %s" % (os.path.getsize(p), rel(p)))


def manifest_paths(path):
    """Every fastq path mentioned, and every sample id."""
    paths, samples = set(), set()
    for i, L in enumerate(lines(path)):
        parts = L.rstrip("\n").split("\t")
        if len(parts) == 1:
            parts = L.rstrip("\n").split(",")
        if i == 0 and not re.search(r"fastq|fq", L):
            continue
        for f in parts:
            f = f.strip()
            if re.search(r"\.f(ast)?q(\.gz)?$", f):
                paths.add(os.path.basename(f))
            elif f and not re.match(r"^(forward|reverse|absolute)", f, re.I):
                samples.add(f)
    return paths, samples


if m16 and m18:
    a, sa = manifest_paths(m16[0])
    b, sb = manifest_paths(m18[0])
    say()
    say("   16S manifest: %d fastq basenames, %d sample ids  (%s)"
        % (len(a), len(sa), rel(m16[0])))
    say("   18S manifest: %d fastq basenames, %d sample ids  (%s)"
        % (len(b), len(sb), rel(m18[0])))
    say()
    if a and a == b:
        say("   >>> IDENTICAL fastq sets.")
        say("   >>> 16S and 18S are ONE amplification, read twice.")
        say("   >>> The MANUSCRIPT is right. Lawrence is misremembering.")
        say("   >>> PR2 extracted to 515F/926R is then CORRECT.")
    elif a and b and not (a & b):
        say("   >>> COMPLETELY DIFFERENT fastq sets, no overlap.")
        say("   >>> 18S was amplified SEPARATELY.")
        say("   >>> LAWRENCE is right, and sections 2.2 and 2.3 are wrong,")
        say("   >>> and the PR2 extraction region is wrong with them.")
    elif a and b:
        say("   >>> PARTIAL overlap: %d shared, %d only in 16S, %d only in 18S."
            % (len(a & b), len(a - b), len(b - a)))
        say("   >>> Mixed. Look at the examples below before concluding.")
        for x in sorted(a - b)[:4]:
            say("        only 16S: %s" % x)
        for x in sorted(b - a)[:4]:
            say("        only 18S: %s" % x)
    else:
        say("   >>> Could not parse fastq paths. First lines of each:")
        for p in (m16[0], m18[0]):
            say("      %s" % rel(p))
            for L in list(lines(p))[:3]:
                say("         %s" % L.rstrip()[:150])
else:
    say("   One or both manifests not found. Showing what exists:")
    for p in find(r"manifest"):
        say("   %10d  %s" % (os.path.getsize(p), rel(p)))

# ---------------------------------------------------------------- TEST 2
head("TEST 2.  WHAT SEQUENCING PROJECTS EXIST?")
say("cytb has its own '-cytoB' projects. A separate 18S amplification would")
say("have '-18s' projects of its own.")
say()
projects = Counter()
for p in find(r"metadata_template|manifest|samplesheet"):
    for L in lines(p, limit=5000):
        for m in re.finditer(r"AMP-[0-9]{6}-[A-Za-z0-9\-]+", L):
            projects[m.group(0)] += 1
if projects:
    suff = Counter()
    for name in projects:
        m = re.search(r"-(16s|18s|cytob|its|cytb)$", name, re.I)
        suff[m.group(1).lower() if m else "other"] += 1
    say("   project-name endings:")
    for k, v in suff.most_common():
        say("      %-8s %d distinct projects" % (k, v))
    say()
    say("   examples:")
    for name, n in projects.most_common(12):
        say("      %s" % name)
else:
    say("   none found")

# ---------------------------------------------------------------- TEST 3
head("TEST 3.  WHAT PRIMERS DID CUTADAPT ACTUALLY TRIM?")
say("Lawrence's 18S forward : %s" % LAWRENCE_F)
say("Lawrence's 18S reverse : %s" % LAWRENCE_R)
say("Manuscript 515F        : %s" % F515)
say("Manuscript 926R        : %s" % R926)
say()
for marker in ("16S", "18S", "cytb"):
    logs = [p for p in find(r"cutadapt.*\.log$") if "/%s/" % marker in p]
    if not logs:
        logs = [p for p in find(r"cutadapt.*\.log$") if marker.lower() in p.lower()]
    say("   --- %s  (%d log files) ---" % (marker, len(logs)))
    seqs, rates = Counter(), Counter()
    for p in logs[:12]:
        for L in lines(p, limit=200000):
            m = re.search(r"Sequence:\s*([ACGTRYSWKMBDHVNacgtryswkmbdhvn]{10,40})\s*;", L)
            if m:
                seqs[m.group(1).upper()] += 1
            m = re.search(r"--error-rate\s+([0-9.]+)", L)
            if m:
                rates[m.group(1)] += 1
            for m in re.finditer(r"--p-front-[fr]\s+([ACGTRYSWKMBDHVN]{10,40})", L):
                seqs[m.group(1).upper()] += 1
            for m in re.finditer(r"\s-[gGaA]\s+([ACGTRYSWKMBDHVN]{10,40})", L):
                seqs[m.group(1).upper()] += 1
    if not seqs:
        say("      no primer sequence found in these logs")
    for s, n in seqs.most_common(8):
        tag = ""
        if s == LAWRENCE_F:
            tag = "  <== LAWRENCE's 18S FORWARD"
        elif s == LAWRENCE_R:
            tag = "  <== LAWRENCE's 18S REVERSE"
        elif s == F515:
            tag = "  <== 515F"
        elif s == R926:
            tag = "  <== 926R"
        elif s.startswith("CCGYCAATTYMTTT"):
            tag = "  <== 926R 14bp core"
        say("      %-26s seen in %d log sections%s" % (s, n, tag))
    if rates:
        say("      error rates: %s" % dict(rates))
    say()

# ---------------------------------------------------------------- TEST 4
head("TEST 4.  CUTADAPT ERROR RATE, PER MARKER")
say("Manuscript 2.3 says 10 percent. The 18S log said --error-rate 0.2.")
say()
for marker in ("16S", "18S", "cytb"):
    logs = [p for p in find(r"cutadapt.*\.log$")
            if ("/%s/" % marker) in p or marker.lower() in os.path.basename(p).lower()]
    rr = Counter()
    for p in logs[:12]:
        for L in lines(p, limit=200000):
            m = re.search(r"--error-rate\s+([0-9.]+)", L)
            if m:
                rr[m.group(1)] += 1
    say("   %-6s %s" % (marker, dict(rr) if rr else "(no error-rate line found)"))

# ---------------------------------------------------------------- TEST 5
head("TEST 5.  WHAT REGION WAS THE PR2 REFERENCE EXTRACTED TO?")
say("If 18S really is the 515F/926R eukaryotic product, extracting PR2 to")
say("515F/926R is correct. If not, the 18S classifier is trained on the")
say("wrong region and every 18S assignment is suspect.")
say()
pat = (r"extract-reads|extract_reads|p-f-primer|p-r-primer|"
       r"pr2|PR2|silva.*extract|feature-classifier")
hits = 0
for p in find(r"\.(py|sh|log|md|txt|yaml|yml)$", maxhits=1200):
    if "/reads_original/" in p or "/reads_staged/" in p:
        continue
    for i, L in enumerate(lines(p, limit=60000), 1):
        if re.search(pat, L) and re.search(r"primer|extract|pr2", L, re.I):
            say("   %s:%d: %s" % (rel(p), i, L.strip()[:165]))
            hits += 1
            if hits >= 26:
                say("   (stopping at 26)")
                break
    if hits >= 26:
        break
if not hits:
    say("   nothing found. Look for the classifier build script by name:")
    for p in find(r"classifier|rescript|train"):
        say("      %s" % rel(p))

# ---------------------------------------------------------------- TEST 6
head("TEST 6.  DENOISE COMMANDS WITH TIMESTAMPS (which run is final)")
rows = []
for p in find(r"dada2.*\.log$|denoise.*\.log$", maxhits=120):
    txt = "".join(lines(p, limit=4000))
    m = re.search(r"denoise-(paired|single)", txt)
    tf = re.search(r"--p-trunc-len-f\s+(\d+)", txt)
    tr = re.search(r"--p-trunc-len-r\s+(\d+)", txt)
    rows.append((os.path.getmtime(p), rel(p),
                 m.group(1) if m else "?",
                 tf.group(1) if tf else "-",
                 tr.group(1) if tr else "-",
                 os.path.getsize(p)))
rows.sort()
say("   %-46s %-7s %-6s %-6s %s" % ("log (oldest first)", "mode", "tr-f", "tr-r", "bytes"))
for mt, p, mode, f, r, sz in rows:
    say("   %-46s %-7s %-6s %-6s %d" % (p[-46:], mode, f, r, sz))

say()
say("=" * 74)
say("done. nothing was written.")
say("=" * 74)
