#!/usr/bin/env python3
"""
read_qpcr2.py -- READ ONLY. Writes nothing, changes nothing.

Second pass on metadata/USGS-qpcr-combined-results.xlsx, now that we know its
shape: one sheet 'in', 6665 rows, 42 columns, with sample_name in column 5,
target_name in 6 and task in 7. The first pass searched columns 1-4 and so
searched filename/well/omit. That was my error; this looks in the right place.

Answers, or fails loudly trying:

  * every one of the 42 column names, so nothing is searched blind again
  * what values 'task' takes  -> whether STANDARD / NTC rows exist at all
  * what values 'target_name' takes, with row counts -> whether there is an
    internal positive control target alongside the seven pathogens
  * for any control rows: their result on each target -> WHICH ASSAY FAILED
  * positives per target under the file's own call column -> tests whether
    991 pathogen positives + 825 control wells = the paper's 1,816
  * the numeric spread of ct / amplification score / Ct confidence by call,
    which is the detection rule the Methods cannot currently state

Usage:  python3 -I read_qpcr2.py [project_root]
"""

import os
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
XLSX = os.path.join(ROOT, "metadata", "USGS-qpcr-combined-results.xlsx")

CONTROL = re.compile(
    r"\bNTC\b|no[-_ ]?template|RT[-_ ]?minus|RT[-_ ]?neg|\bNRT\b|"
    r"pos[-_ ]?c|\bPC\b|plasmid|gblock|g[-_ ]block|synthetic|"
    r"\bstd\b|standard|\bblank\b|\bwater\b|\bH2O\b|neg[-_ ]?c|control|"
    r"\bcal\b|\bIPC\b|spike",
    re.I)

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def say(s=""):
    print(s)


def head(t):
    say()
    say("=" * 74)
    say(t)
    say("=" * 74)


def colnum(ref):
    n = 0
    for ch in ref:
        if ch.isalpha():
            n = n * 26 + (ord(ch.upper()) - 64)
        else:
            break
    return n


def load(path):
    """Every row of every sheet, as lists of strings. stdlib only."""
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(NS + "si"):
            shared.append("".join(t.text or "" for t in si.iter(NS + "t")))
    rels = {}
    if "xl/_rels/workbook.xml.rels" in z.namelist():
        for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels")):
            rels[r.get("Id")] = r.get("Target").lstrip("/")
    out = []
    for sh in ET.fromstring(z.read("xl/workbook.xml")).iter(NS + "sheet"):
        tgt = rels.get(sh.get(RNS + "id"), "")
        if tgt and not tgt.startswith("xl/"):
            tgt = "xl/" + tgt
        if tgt not in z.namelist():
            continue
        rows = []
        for row in ET.fromstring(z.read(tgt)).iter(NS + "row"):
            cells = {}
            for c in row.findall(NS + "c"):
                ci = colnum(c.get("r") or "A")
                if c.get("t") == "s":
                    v = c.find(NS + "v")
                    i = int(v.text) if v is not None and v.text else -1
                    val = shared[i] if 0 <= i < len(shared) else ""
                elif c.get("t") == "inlineStr":
                    val = "".join(t.text or "" for t in c.iter(NS + "t"))
                else:
                    v = c.find(NS + "v")
                    val = v.text if v is not None and v.text else ""
                cells[ci] = val
            width = max(cells) if cells else 0
            rows.append([cells.get(i + 1, "") for i in range(width)])
        out.append((sh.get("name"), rows))
    z.close()
    return out


if not os.path.exists(XLSX):
    say("NOT FOUND: %s" % XLSX)
    sys.exit(1)

sheets = load(XLSX)
name, rows = sheets[0]
hdr = [h.strip() for h in rows[0]]
data = rows[1:]

# index by column name, tolerant of case and spacing
IX = {}
for i, h in enumerate(hdr):
    IX[h.lower().replace(" ", "_")] = i


def col(row, key, default=""):
    i = IX.get(key)
    if i is None or i >= len(row):
        return default
    return (row[i] or "").strip()


def find_col(*cands):
    """First column whose name contains any candidate substring."""
    for c in cands:
        for k in IX:
            if c in k:
                return k
    return None


head("0.  ALL %d COLUMN NAMES  (sheet %r, %d data rows)" % (len(hdr), name, len(data)))
for i, h in enumerate(hdr, 1):
    sample = ""
    for r in data[:400]:
        v = (r[i - 1] if i - 1 < len(r) else "").strip()
        if v and v.upper() not in ("NA", "NAN"):
            sample = v
            break
    say("   %2d  %-34s  e.g. %s" % (i, h, sample[:28]))

# ---------------------------------------------------------------- task values
head("1.  'task' VALUES  --  are STANDARD / NTC rows declared?")
tk = find_col("task")
if tk is None:
    say("   no task-like column")
else:
    for v, n in Counter(col(r, tk) for r in data).most_common():
        say("   %-24s %6d rows" % (repr(v), n))

# -------------------------------------------------------------- omit values
head("2.  'omit' VALUES  --  rows the instrument or operator excluded")
om = find_col("omit")
if om:
    for v, n in Counter(col(r, om) for r in data).most_common(10):
        say("   %-24s %6d rows" % (repr(v), n))

# ------------------------------------------------------------ target values
head("3.  'target_name' VALUES  --  7 pathogens, or more?")
tg = find_col("target_name", "target")
tgcounts = Counter(col(r, tg) for r in data) if tg else Counter()
for v, n in tgcounts.most_common():
    say("   %-40s %6d rows" % (repr(v), n))
say()
say("   %d distinct targets, %d rows total" % (len(tgcounts), sum(tgcounts.values())))
say("   The recount table had 825 wells x 7 pathogens = 5,775.")
say("   This file has %d rows. The difference is %d."
    % (len(data), len(data) - 5775))

# ----------------------------------------------------------- control rows
head("4.  CONTROL ROWS  --  searching sample_name AND task, properly this time")
sn = find_col("sample_name", "sample")
hits = []
for r in data:
    blob = " ".join([col(r, sn), col(r, tk) if tk else ""])
    if CONTROL.search(blob):
        hits.append(r)
say("   %d rows look like a control" % len(hits))
if hits:
    names = Counter(col(r, sn) for r in hits)
    say()
    say("   distinct control sample names:")
    for v, n in names.most_common(30):
        say("      %-34s %4d rows" % (repr(v), n))

    # the money question: per control sample, the result on each target
    ctcol = find_col("ct")
    callcol = find_col("call", "detect", "result", "verdict", "interp")
    say()
    say("   PER CONTROL SAMPLE, RESULT ON EACH TARGET")
    say("   A positive control negative for exactly one target names the")
    say("   assay that failed.")
    say()
    bysample = defaultdict(dict)
    for r in hits:
        bysample[col(r, sn)][col(r, tg)] = (
            col(r, ctcol) if ctcol else "", col(r, callcol) if callcol else "")
    alltg = sorted({t for d in bysample.values() for t in d})
    for s in sorted(bysample):
        say("   %s" % s)
        for t in alltg:
            ct, cl = bysample[s].get(t, ("--", "--"))
            say("      %-34s ct=%-10s call=%s" % (t, ct, cl))
        say()
else:
    say()
    say("   Still nothing. Then the controls are not named as samples in this")
    say("   file. The three rows the manuscript calls positive controls were")
    say("   identified by their RESULT, which is the lab reviewer's point, and")
    say("   it means no plate map backs them. Check section 5 below: any")
    say("   sample positive on 6 of 7 targets is one of them.")

# ------------------------------------------- samples positive on many targets
head("5.  SAMPLES POSITIVE ON SIX OR SEVEN TARGETS")
say("   These are the de facto positive controls. The target each is NEGATIVE")
say("   for is the assay that failed.")
say()
ctcol = find_col("ct")
callcol = find_col("call", "detect", "result", "verdict", "interp")


def is_pos(r):
    """Positive by the file's own call column if it has one, else by a Ct
    below the no-amplification sentinel (the first rows show ct=41 for
    nothing detected, so 41 and blank both mean negative)."""
    if callcol:
        v = col(r, callcol).lower()
        if v in ("1", "true", "yes", "y", "pos", "positive", "detected"):
            return True
        if v in ("0", "false", "no", "n", "neg", "negative", "not detected", "undetermined"):
            return False
    v = col(r, ctcol)
    try:
        f = float(v)
    except ValueError:
        return False
    return 0 < f < 40.0


pos_by_sample = defaultdict(set)
tg_by_sample = defaultdict(set)
for r in data:
    s = col(r, sn)
    tg_by_sample[s].add(col(r, tg))
    if is_pos(r):
        pos_by_sample[s].add(col(r, tg))

pathogens = [t for t in tgcounts if t]
wide = [(len(v), s) for s, v in pos_by_sample.items() if len(v) >= 6]
wide.sort(reverse=True)
if not wide:
    say("   none at 6+. Top samples by targets positive:")
    wide = sorted(((len(v), s) for s, v in pos_by_sample.items()), reverse=True)[:6]
for n, s in wide[:12]:
    missing = sorted(set(pathogens) - pos_by_sample[s])
    say("   %-28s positive on %d targets" % (s, n))
    say("      NEGATIVE for: %s" % (", ".join(missing) if missing else "(nothing)"))

# --------------------------------------------------- positives per target
head("6.  POSITIVES PER TARGET ACROSS ALL %d ROWS" % len(data))
say("   Tests the arithmetic: does narrow-rule 991 + 825 control wells")
say("   = the paper's 1,816?")
say()
say("   using call column: %s" % (callcol or "(none -- fell back to Ct < 40)"))
say()
tot = 0
say("   %-40s %7s %7s" % ("target", "rows", "pos"))
say("   " + "-" * 58)
for t in sorted(tgcounts):
    rowsN = tgcounts[t]
    posN = sum(1 for r in data if col(r, tg) == t and is_pos(r))
    tot += posN
    say("   %-40s %7d %7d" % (t, rowsN, posN))
say("   " + "-" * 58)
say("   %-40s %7d %7d" % ("TOTAL", len(data), tot))
say()
say("   paper says 1,816 of 6,664.  this file: %d of %d" % (tot, len(data)))

# ---------------------------------------------------------- the actual rule
head("7.  THE DETECTION RULE  --  numeric spread by call")
say("   Whatever column separates a detection from a non-detection is the")
say("   rule the Methods has to state. Looking for the boundary.")
say()
numcols = [k for k in IX
           if any(w in k for w in ("ct", "amp", "conf", "score", "quantity", "cq"))]
for k in sorted(numcols):
    vals_pos, vals_neg = [], []
    for r in data:
        v = col(r, k)
        try:
            f = float(v)
        except ValueError:
            continue
        (vals_pos if is_pos(r) else vals_neg).append(f)
    if not vals_pos and not vals_neg:
        continue

    def spread(v):
        if not v:
            return "n=0"
        v = sorted(v)
        return ("n=%-6d min=%-10.4g p05=%-10.4g med=%-10.4g p95=%-10.4g max=%-10.4g"
                % (len(v), v[0], v[int(.05 * (len(v) - 1))], v[len(v) // 2],
                   v[int(.95 * (len(v) - 1))], v[-1]))
    say("   %s" % hdr[IX[k]])
    say("      detections    %s" % spread(vals_pos))
    say("      non-detections %s" % spread(vals_neg))
    say()

say("=" * 74)
say("done. nothing was written.")
say("=" * 74)
