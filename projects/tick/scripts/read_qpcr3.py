#!/usr/bin/env python3
"""
read_qpcr3.py -- READ ONLY. Writes nothing, changes nothing.

Last pass. Column 42 is 'pass_check', the file's own verdict, which the two
earlier passes both missed. With that plus the QuantStudio quality flags
(noamp / expfail / tholdfail) the detection rule should be statable in one
sentence instead of described as proprietary.

  1. pass_check value counts, overall and per target
  2. pass_check against the quality flags, to see what drives the verdict
  3. the numeric boundary: amp_score and cq_conf by verdict
  4. the 8 subsamples in this file (833) but not in the recount (825)
  5. the four samples positive on six of seven pathogens, every column
  6. every E. muris and B. mayonii well that amplified, with its verdict

Usage:  python3 -I read_qpcr3.py [project_root]
"""

import csv
import os
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
XLSX = os.path.join(ROOT, "metadata", "USGS-qpcr-combined-results.xlsx")
RECOUNT = os.path.join(ROOT, "results", "analysis", "qpcr_recount_by_rule.tsv")

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

WIDE = ["alht-deer-84", "alht-deer-81", "alht-deer-67", "alht-deer-21"]


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
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(NS + "si"):
            shared.append("".join(t.text or "" for t in si.iter(NS + "t")))
    rels = {}
    for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels")):
        rels[r.get("Id")] = r.get("Target").lstrip("/")
    sh = next(ET.fromstring(z.read("xl/workbook.xml")).iter(NS + "sheet"))
    tgt = rels.get(sh.get(RNS + "id"), "")
    if tgt and not tgt.startswith("xl/"):
        tgt = "xl/" + tgt
    out = []
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
        w = max(cells) if cells else 0
        out.append([cells.get(i + 1, "") for i in range(w)])
    z.close()
    return out


if not os.path.exists(XLSX):
    say("NOT FOUND: %s" % XLSX)
    sys.exit(1)

rows = load(XLSX)
hdr = [h.strip() for h in rows[0]]
data = rows[1:]
IX = {h.lower().replace(" ", "_"): i for i, h in enumerate(hdr)}


def g(r, key):
    i = IX.get(key)
    return (r[i] or "").strip() if i is not None and i < len(r) else ""


def fnum(r, key):
    try:
        return float(g(r, key))
    except ValueError:
        return None


def amplified(r):
    c = fnum(r, "ct")
    return c is not None and c < 40.5


# ---------------------------------------------------------------- 1
head("1.  pass_check VALUES")
say("   overall:")
for v, n in Counter(g(r, "pass_check") for r in data).most_common():
    say("      %-26s %6d" % (repr(v), n))

say()
say("   per target:")
pt = defaultdict(Counter)
for r in data:
    pt[g(r, "target_name")][g(r, "pass_check")] += 1
verdicts = sorted({v for c in pt.values() for v in c})
say("      %-24s %s" % ("target", "".join("%-16s" % v[:15] for v in verdicts)))
say("      " + "-" * (24 + 16 * len(verdicts)))
tot = Counter()
for t in sorted(pt):
    say("      %-24s %s" % (t[:24], "".join("%-16d" % pt[t][v] for v in verdicts)))
    for v in verdicts:
        tot[v] += pt[t][v]
say("      " + "-" * (24 + 16 * len(verdicts)))
say("      %-24s %s" % ("TOTAL", "".join("%-16d" % tot[v] for v in verdicts)))
say()
say("   paper's headline numerator is 1,816 of 6,664 wells.")
say("   Pathogen-only totals (excluding Extraction Control):")
for v in verdicts:
    n = sum(pt[t][v] for t in pt if "xtraction" not in t)
    say("      %-26s %6d" % (repr(v), n))

# ---------------------------------------------------------------- 2
head("2.  pass_check AGAINST THE QUALITY FLAGS")
flags = [k for k in ("noamp", "expfail", "tholdfail", "ctfail", "prflow",
                     "ampnc", "prfdrop", "highsd", "badrox", "noise",
                     "spike", "nosignal") if k in IX]
say("   flags present: %s" % ", ".join(flags))
say()
combo = Counter()
for r in data:
    sig = tuple((f, g(r, f)) for f in flags if g(r, f) not in ("", "N"))
    combo[(g(r, "pass_check"), sig)] += 1
for (verdict, sig), n in combo.most_common(22):
    desc = ", ".join("%s=%s" % kv for kv in sig) or "(no flags raised)"
    say("   %-14s %-52s %6d" % (verdict[:14], desc[:52], n))

# ---------------------------------------------------------------- 3
head("3.  THE NUMERIC BOUNDARY, BY VERDICT")
for key in ("ct", "amp_score", "cq_conf"):
    if key not in IX:
        continue
    say("   %s" % key)
    byv = defaultdict(list)
    for r in data:
        v = fnum(r, key)
        if v is not None:
            byv[g(r, "pass_check")].append(v)
    for v in sorted(byv):
        s = sorted(byv[v])
        say("      %-14s n=%-6d min=%-10.5g med=%-10.5g max=%-10.5g"
            % (v[:14], len(s), s[0], s[len(s) // 2], s[-1]))
    say()
say("   If min(amp_score) among positives is above max(amp_score) among")
say("   negatives, the rule is a clean threshold and can simply be stated.")
say("   If they overlap, something else is deciding and we need to say what.")

# ---------------------------------------------------------------- 4
head("4.  THE 8 SUBSAMPLES DROPPED BETWEEN THIS FILE AND THE RECOUNT")
xl_samples = {g(r, "sample_name") for r in data}
xl_hosts = {g(r, "host_id") for r in data}
say("   distinct sample_name in xlsx : %d" % len(xl_samples))
say("   distinct host_id   in xlsx   : %d" % len(xl_hosts))
if os.path.exists(RECOUNT):
    rc = list(csv.DictReader(open(RECOUNT), delimiter="\t"))
    rc_animals = {(r.get("animal") or "").strip() for r in rc}
    say("   distinct animal    in recount: %d" % len(rc_animals))
    say()
    lo = {a.lower() for a in rc_animals}
    missing = sorted(s for s in xl_samples if s and s.lower() not in lo)
    say("   in xlsx, NOT in recount (%d):" % len(missing))
    for s in missing[:40]:
        npos = sum(1 for r in data
                   if g(r, "sample_name") == s and amplified(r)
                   and "xtraction" not in g(r, "target_name"))
        say("      %-30s  pathogen targets amplified: %d" % (s, npos))
    if len(missing) > 40:
        say("      ... and %d more" % (len(missing) - 40))
    say()
    say("   NOTE: names may differ in form between the two files, so a long")
    say("   list here means a naming mismatch rather than 8 real exclusions.")
else:
    say("   recount not found at %s" % RECOUNT)

# ---------------------------------------------------------------- 5
head("5.  THE FOUR SAMPLES POSITIVE ON SIX OF SEVEN PATHOGENS")
say("   Every column, so we can see plate, well and flags together.")
for s in WIDE:
    sub = [r for r in data if g(r, "sample_name").lower() == s]
    if not sub:
        say()
        say("   %s -- no rows under that exact name" % s)
        continue
    say()
    say("   %s   (%d rows, plate %s)"
        % (s, len(sub), ", ".join(sorted({g(r, "qpcrplate") for r in sub}))))
    say("      %-22s %-6s %-8s %-9s %-10s %s"
        % ("target", "well", "ct", "amp_score", "cq_conf", "pass_check"))
    for r in sorted(sub, key=lambda x: g(x, "target_name")):
        say("      %-22s %-6s %-8s %-9s %-10s %s"
            % (g(r, "target_name")[:22], g(r, "well_position"),
               g(r, "ct")[:8], g(r, "amp_score")[:9],
               g(r, "cq_conf")[:10], g(r, "pass_check")))

# ---------------------------------------------------------------- 6
head("6.  EVERY AMPLIFIED WELL FOR THE TWO 'ABSENT' TARGETS")
say("   The paper reports both as not detected. These wells amplified.")
for target in ("E. muris", "B. mayonii"):
    sub = [r for r in data if g(r, "target_name") == target and amplified(r)]
    say()
    say("   %s -- %d wells with a real Ct" % (target, len(sub)))
    say("      %-26s %-6s %-8s %-9s %-10s %s"
        % ("sample", "well", "ct", "amp_score", "cq_conf", "pass_check"))
    for r in sorted(sub, key=lambda x: fnum(x, "ct") or 99)[:45]:
        say("      %-26s %-6s %-8s %-9s %-10s %s"
            % (g(r, "sample_name")[:26], g(r, "well_position"),
               g(r, "ct")[:8], g(r, "amp_score")[:9],
               g(r, "cq_conf")[:10], g(r, "pass_check")))
    if len(sub) > 45:
        say("      ... and %d more" % (len(sub) - 45))

say()
say("=" * 74)
say("done. nothing was written.")
say("=" * 74)
