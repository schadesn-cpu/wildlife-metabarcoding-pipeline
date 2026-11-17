#!/usr/bin/env python3
"""
find_answers.py -- READ ONLY. Writes nothing, changes nothing.

Purpose: find out which of the open Methods questions can be answered from the
project itself, so we stop asking people for things the data already holds.

This pass only LOCATES and SAMPLES. It does not compute final numbers. Run it,
paste the output back, and the targeted queries come next.

Usage:   python3 -I find_answers.py
         python3 -I find_answers.py /path/to/project/root
"""

import os
import re
import sys

CAP = 14          # most lines printed per section
PEEK = 3          # most files sampled per section


def say(s=""):
    print(s)


def head(n, title):
    say()
    say("=" * 72)
    say("%s. %s" % (n, title))
    say("=" * 72)


# ---------------------------------------------------------------- locate root
def find_root(argv):
    if len(argv) > 1:
        return os.path.abspath(argv[1])
    here = os.path.abspath(os.getcwd())
    # walk up looking for something that smells like the project
    probe = here
    for _ in range(6):
        names = set(os.listdir(probe)) if os.path.isdir(probe) else set()
        if {"results", "data"} & names or {"scripts", "data"} & names:
            return probe
        up = os.path.dirname(probe)
        if up == probe:
            break
        probe = up
    return here


ROOT = find_root(sys.argv)

SKIP_DIR = re.compile(r"(^\.|__pycache__|node_modules|\.git$|site-packages)")


def walk(exts=None, pat=None, maxhits=400):
    """Every file under ROOT matching an extension list and/or a name regex."""
    hits = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if not SKIP_DIR.search(d)]
        for fn in filenames:
            if exts and not fn.lower().endswith(tuple(exts)):
                continue
            if pat and not pat.search(fn):
                continue
            p = os.path.join(dirpath, fn)
            try:
                hits.append((os.path.getsize(p), p))
            except OSError:
                pass
            if len(hits) >= maxhits:
                return hits
    return hits


def rel(p):
    try:
        return os.path.relpath(p, ROOT)
    except ValueError:
        return p


def show(hits, cap=CAP, label="files"):
    hits.sort(reverse=True)
    if not hits:
        say("   (none found)")
        return []
    for size, p in hits[:cap]:
        say("   %9d  %s" % (size, rel(p)))
    if len(hits) > cap:
        say("   ... and %d more %s" % (len(hits) - cap, label))
    return [p for _, p in hits]


def headerline(path):
    """First line of a text file, or None."""
    try:
        with open(path, "r", errors="replace") as fh:
            return fh.readline().rstrip("\n")
    except OSError:
        return None


def grep(paths, pattern, cap=CAP, ctx=0):
    """Print matching lines across paths. Read-only, line-bounded."""
    rx = re.compile(pattern, re.I)
    shown = 0
    for p in paths:
        try:
            with open(p, "r", errors="replace") as fh:
                for i, line in enumerate(fh, 1):
                    if i > 200000:
                        break
                    if rx.search(line):
                        say("   %s:%d: %s" % (rel(p), i, line.strip()[:160]))
                        shown += 1
                        if shown >= cap:
                            say("   (stopping at %d matches)" % cap)
                            return shown
        except OSError:
            pass
    if not shown:
        say("   (no matches)")
    return shown


say("ROOT = %s" % ROOT)
say("read-only: this script opens files and prints; it writes nothing")

# ------------------------------------------------- 1. qPCR well-level data
head(1, "qPCR WELL-LEVEL DATA  (items 1, 3, 4, 5, 6, 7)")
say("Anything holding one row per WELL rather than per animal. We need")
say("amplification score and Ct confidence, which the README cites, plus")
say("whatever the sample/well names are.")
say()
qpcr = walk(exts=[".csv", ".tsv", ".txt", ".xlsx", ".xls"],
            pat=re.compile(r"qpcr|well|plate|ct_|amp|adjud|call", re.I))
qpcr_paths = show(qpcr)
say()
say("   -- column headers of the largest few --")
for p in qpcr_paths[:PEEK]:
    if p.lower().endswith((".xlsx", ".xls")):
        say("   %s: (spreadsheet, needs a reader -- noting it exists)" % rel(p))
        continue
    h = headerline(p)
    if h:
        say("   %s:" % rel(p))
        say("      %s" % h[:300])

# ------------------------------------------------- 2. controls on the plates
head(2, "NO-TEMPLATE / RT-MINUS / BLANK WELLS  (items 4, 5)")
say("If an NTC or an RT-minus ever ran, it is a named sample somewhere.")
say()
grep(qpcr_paths[:6], r"\bNTC\b|no[-_ ]?template|RT[-_ ]?minus|RT[-_ ]?neg|"
                     r"\bNRT\b|water|\bblank\b|\bneg[-_ ]?ctrl|negative control")

# ------------------------------------------------- 3. positive controls
head(3, "POSITIVE CONTROLS  (items 1, 6)")
say("The three rows positive for six of seven targets. We want to see WHICH")
say("target they were negative for: that is the assay that failed.")
say()
grep(qpcr_paths[:6], r"pos[-_ ]?ctrl|positive control|\bPC\b|plasmid|"
                     r"gblock|g-block|synthetic|\bstandard\b|\bSTD\b")
say()
say("   -- any file that looks like a plate MAP rather than results --")
show(walk(pat=re.compile(r"plate.*(map|layout|design)|layout|manifest", re.I)), cap=8)

# ------------------------------------------------- 4. re-adjudication
head(4, "RE-ADJUDICATION  (item 7)")
say("Twenty wells moved into detections and one out. Two versions of the")
say("calls would let us diff them and name the targets.")
say()
show(walk(pat=re.compile(r"adjud|re[-_]?call|revis|v1|v2|_old|_orig|final", re.I),
          exts=[".csv", ".tsv", ".txt"]), cap=10)

# ------------------------------------------------- 5. sequencing platform
head(5, "SEQUENCING PLATFORM  (item 8)")
say("Illumina read IDs encode the instrument. One header answers this.")
say()
fq = walk(pat=re.compile(r"\.(fastq|fq)(\.gz)?$", re.I), maxhits=60)
fq_paths = show(fq, cap=6)
say()
for p in fq_paths[:2]:
    if p.endswith(".gz"):
        say("   %s: gzipped. Run:" % rel(p))
        say("      zcat '%s' | head -1" % p)
    else:
        say("   %s" % rel(p))
        say("      %s" % (headerline(p) or "")[:200])
say()
say("   -- run-info / sample-sheet files, which also carry the platform --")
show(walk(pat=re.compile(r"runinfo|run_info|samplesheet|sample_sheet|"
                         r"runparam|RunParameters|stats\.json", re.I)), cap=8)

# ------------------------------------------------- 6. pipeline parameters
head(6, "TRIMMING / DENOISING / ASSIGNMENT PARAMETERS  (item 9)")
say("The settings the 666 variants and the 182-328 bp range rest on.")
say()
code = walk(exts=[".sh", ".py", ".r", ".rmd", ".qmd", ".nf", ".yaml", ".yml",
                  ".json", ".log", ".md"], maxhits=600)
code_paths = [p for _, p in sorted(code, reverse=True)]
say("   -- denoising calls --")
grep(code_paths, r"filterAndTrim|learnErrors|dada2|dada\(|DADA2|deblur|"
                 r"truncLen|maxEE|trimLeft|minOverlap|unoise|usearch", cap=16)
say()
say("   -- primer trimming --")
grep(code_paths, r"cutadapt|trimmomatic|fastp|bbduk|--discard-untrimmed", cap=10)
say()
say("   -- identity / coverage thresholds for assignment --")
grep(code_paths, r"perc_identity|pident|--id\b|minBoot|identity.*9[0-9]|"
                 r"qcov|coverage.*9[0-9]|98\.7|94\.5", cap=12)

# ------------------------------------------------- 7. BLAST database
head(7, "BLAST DATABASE  (item 10)")
say("The paper describes it three ways. One command line settles it.")
say()
grep(code_paths, r"blastn|blastx|makeblastdb|core_nt|\bnt\b.*-db|-db\s|"
                 r"BLASTDB|update_blastdb", cap=16)
say()
say("   -- the database directory itself, whose timestamps do not lie --")
for cand in ("/db", "/data/db", os.path.expanduser("~/db"),
             os.path.join(ROOT, "db"), "/usr/local/share/blastdb"):
    if os.path.isdir(cand):
        say("   %s exists:" % cand)
        try:
            for fn in sorted(os.listdir(cand))[:8]:
                fp = os.path.join(cand, fn)
                say("      %s" % fn)
        except OSError:
            pass

# ------------------------------------------------- 8. Bacillus read fraction
head(8, "BACILLUS SPIKE READ FRACTION  (item 12)")
say("We spiked every extraction. The 16S ASV table knows what it took.")
say()
asv = walk(pat=re.compile(r"asv|otu|feature[-_]?table|taxonom|counts|abundan",
                          re.I),
           exts=[".tsv", ".csv", ".txt", ".biom"])
asv_paths = show(asv, cap=10)
say()
for p in asv_paths[:PEEK]:
    h = headerline(p)
    if h:
        say("   %s:" % rel(p))
        say("      %s" % h[:240])
say()
say("   -- where Bacillus is named --")
grep(asv_paths[:6], r"Bacillus|atrophaeus", cap=8)

# ------------------------------------------------- 9. anything else
head(9, "LAB NOTES / RUN RECORDS  (items 2, 11)")
say("Platform, cycling conditions and the Powassan control are the two that")
say("probably are not in the data. Worth one look for a notebook anyway.")
say()
show(walk(exts=[".md", ".txt", ".docx", ".pdf"],
          pat=re.compile(r"protocol|notebook|lab|method|run|readme|sop|"
                         r"cycl|thermo|master.?mix", re.I)), cap=12)

say()
say("=" * 72)
say("done. nothing was written.")
say("=" * 72)
