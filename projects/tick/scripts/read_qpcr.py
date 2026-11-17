#!/usr/bin/env python3
"""
read_qpcr.py -- READ ONLY. Writes nothing, changes nothing.

Two jobs:

  A. Open metadata/USGS-qpcr-combined-results.xlsx and show its real shape:
     sheet names, headers, sample rows, and every row whose sample name looks
     like a control. That is where the answer to "which assay failed in every
     positive control" lives, and whether a no-template control ever ran on a
     qPCR plate.

  B. Summarise results/analysis/qpcr_recount_by_rule.tsv, which already holds
     detections under three rules side by side. That is the detection-rule
     question and probably the re-adjudication question too.

Uses openpyxl if present, falls back to parsing the xlsx as a zip of XML, so
it runs in any environment.

Usage:  python3 -I read_qpcr.py
        python3 -I read_qpcr.py /path/to/project/root
"""

import csv
import os
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import defaultdict

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())

XLSX = os.path.join(ROOT, "metadata", "USGS-qpcr-combined-results.xlsx")
RECOUNT = os.path.join(ROOT, "results", "analysis", "qpcr_recount_by_rule.tsv")

# a sample name that looks like a control rather than an animal
CONTROL = re.compile(
    r"\bNTC\b|no[-_ ]?template|RT[-_ ]?minus|RT[-_ ]?neg|\bNRT\b|"
    r"pos[-_ ]?c|\bPC\b|plasmid|gblock|g[-_ ]block|synthetic|"
    r"\bstd\b|standard|\bblank\b|\bwater\b|\bH2O\b|neg[-_ ]?c|control",
    re.I)


def say(s=""):
    print(s)


def head(t):
    say()
    say("=" * 72)
    say(t)
    say("=" * 72)


# --------------------------------------------------------------- xlsx reading
def rows_openpyxl(path, cap):
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    out = []
    for ws in wb.worksheets:
        rows = []
        for i, row in enumerate(ws.iter_rows(values_only=True)):
            if i >= cap:
                break
            rows.append(["" if c is None else str(c) for c in row])
        out.append((ws.title, ws.max_row, ws.max_column, rows))
    wb.close()
    return out


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def colnum(ref):
    """'BC12' -> 54 (1-based column)."""
    n = 0
    for ch in ref:
        if ch.isalpha():
            n = n * 26 + (ord(ch.upper()) - 64)
        else:
            break
    return n


def rows_stdlib(path, cap):
    """Parse the xlsx as what it is: a zip of XML."""
    z = zipfile.ZipFile(path)

    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall(NS + "si"):
            shared.append("".join(t.text or "" for t in si.iter(NS + "t")))

    rels = {}
    if "xl/_rels/workbook.xml.rels" in z.namelist():
        for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels")):
            rels[r.get("Id")] = r.get("Target").lstrip("/")

    sheets = []
    for sh in ET.fromstring(z.read("xl/workbook.xml")).iter(NS + "sheet"):
        tgt = rels.get(sh.get(RNS + "id"), "")
        if tgt and not tgt.startswith("xl/"):
            tgt = "xl/" + tgt
        sheets.append((sh.get("name"), tgt))

    out = []
    for name, tgt in sheets:
        if tgt not in z.namelist():
            out.append((name, 0, 0, []))
            continue
        rows, nrow, ncol = [], 0, 0
        root = ET.fromstring(z.read(tgt))
        for row in root.iter(NS + "row"):
            nrow += 1
            if len(rows) < cap:
                cells = {}
                for c in row.findall(NS + "c"):
                    ci = colnum(c.get("r") or "A")
                    if c.get("t") == "s":
                        v = c.find(NS + "v")
                        idx = int(v.text) if v is not None and v.text else -1
                        val = shared[idx] if 0 <= idx < len(shared) else ""
                    elif c.get("t") == "inlineStr":
                        val = "".join(t.text or "" for t in c.iter(NS + "t"))
                    else:
                        v = c.find(NS + "v")
                        val = v.text if v is not None and v.text else ""
                    cells[ci] = val
                    ncol = max(ncol, ci)
                rows.append([cells.get(i + 1, "") for i in range(max(cells) if cells else 0)])
        out.append((name, nrow, ncol, rows))
    z.close()
    return out


def read_xlsx(path, cap=40):
    try:
        return rows_openpyxl(path, cap), "openpyxl"
    except ImportError:
        return rows_stdlib(path, cap), "stdlib zip/XML"
    except Exception as e:
        say("   openpyxl failed (%s); falling back" % e)
        return rows_stdlib(path, cap), "stdlib zip/XML"


# ------------------------------------------------------------------ job A
head("A.  metadata/USGS-qpcr-combined-results.xlsx")

if not os.path.exists(XLSX):
    say("   NOT FOUND at %s" % XLSX)
    say("   Pass the project root as an argument, or tell me the real path.")
else:
    sheets, how = read_xlsx(XLSX)
    say("   read with: %s" % how)
    say()
    for name, nrow, ncol, rows in sheets:
        say("   SHEET %-28s %6s rows x %3s cols" % (repr(name), nrow, ncol))
    say()

    for name, nrow, ncol, rows in sheets:
        if not rows:
            continue
        say("   " + "-" * 66)
        say("   SHEET %r -- header and first rows" % name)
        say("   " + "-" * 66)
        for i, r in enumerate(rows[:6]):
            say("   %2d | %s" % (i + 1, " | ".join(x[:22] for x in r[:14])))
        say()

    # now hunt for control rows across the whole of every sheet
    say("   " + "-" * 66)
    say("   ROWS WHOSE SAMPLE NAME LOOKS LIKE A CONTROL")
    say("   (positive controls, NTC, RT-minus, plasmid, standards)")
    say("   " + "-" * 66)
    sheets_full, _ = read_xlsx(XLSX, cap=200000)
    found = 0
    for name, nrow, ncol, rows in sheets_full:
        hdr = rows[0] if rows else []
        for r in rows[1:]:
            joined = " ".join(r[:4])
            if CONTROL.search(joined):
                if found == 0 and hdr:
                    say("   HEADER: %s" % " | ".join(x[:20] for x in hdr[:16]))
                say("   [%s] %s" % (name, " | ".join(x[:20] for x in r[:16])))
                found += 1
                if found >= 60:
                    say("   (stopping at 60)")
                    break
        if found >= 60:
            break
    if not found:
        say("   NONE. No row on any sheet carries a control-looking name.")
        say("   If the controls are not in this file they are somewhere else,")
        say("   or they were never recorded as named samples.")


# ------------------------------------------------------------------ job B
head("B.  results/analysis/qpcr_recount_by_rule.tsv")

if not os.path.exists(RECOUNT):
    say("   NOT FOUND at %s" % RECOUNT)
else:
    rows = list(csv.DictReader(open(RECOUNT), delimiter="\t"))
    say("   %d rows" % len(rows))
    say("   columns: %s" % ", ".join(rows[0].keys()) if rows else "")
    say()

    def num(x):
        x = (x or "").strip()
        try:
            return float(x)
        except ValueError:
            return 0.0

    def truthy(x):
        return (x or "").strip().lower() in ("1", "true", "yes", "y", "t")

    RULES = ["recorded", "stated", "narrow"]

    say("   PER TARGET: wells tested, positive wells under each rule,")
    say("   and animals called detected under each rule")
    say()
    say("   %-34s %7s | %7s %7s %7s | %6s %6s %6s"
        % ("target", "wells", "w:rec", "w:stat", "w:narr",
           "a:rec", "a:stat", "a:narr"))
    say("   " + "-" * 92)

    by = defaultdict(lambda: defaultdict(float))
    animals = defaultdict(lambda: defaultdict(set))
    for r in rows:
        t = (r.get("target") or "?").strip()
        by[t]["wells"] += num(r.get("wells"))
        for k in RULES:
            by[t]["w_" + k] += num(r.get(k + "_positive"))
            if truthy(r.get("detected_" + k)):
                animals[t][k].add((r.get("animal") or "").strip())
        animals[t]["all"].add((r.get("animal") or "").strip())

    tot = defaultdict(float)
    for t in sorted(by):
        say("   %-34s %7.0f | %7.0f %7.0f %7.0f | %6d %6d %6d"
            % (t, by[t]["wells"],
               by[t]["w_recorded"], by[t]["w_stated"], by[t]["w_narrow"],
               len(animals[t]["recorded"]), len(animals[t]["stated"]),
               len(animals[t]["narrow"])))
        for k in ("wells", "w_recorded", "w_stated", "w_narrow"):
            tot[k] += by[t][k]
    say("   " + "-" * 92)
    say("   %-34s %7.0f | %7.0f %7.0f %7.0f"
        % ("TOTAL", tot["wells"], tot["w_recorded"],
           tot["w_stated"], tot["w_narrow"]))
    say()
    say("   The paper's headline is 1,816 detections of 6,664 wells.")
    say("   Whichever column above equals 1,816 is the rule that was used.")
    say()

    # where the rules disagree, per target
    say("   WHERE THE RULES DISAGREE, by target and animal")
    say("   (this is the re-adjudication: which calls moved, and on what)")
    say()
    n = 0
    for r in rows:
        d = {k: truthy(r.get("detected_" + k)) for k in RULES}
        if len(set(d.values())) > 1:
            n += 1
            if n <= 40:
                say("   %-28s %-30s rec=%-5s stat=%-5s narr=%-5s  wells=%s"
                    % ((r.get("animal") or "")[:28],
                       (r.get("target") or "")[:30],
                       d["recorded"], d["stated"], d["narrow"],
                       r.get("wells")))
    if n > 40:
        say("   ... and %d more disagreeing animal/target pairs" % (n - 40))
    say()
    say("   %d animal/target pairs change status between rules." % n)
    if n == 0:
        say("   The three rules agree everywhere, so the rule choice does not")
        say("   move any animal's status. That is worth one sentence in Methods.")

say()
say("=" * 72)
say("done. nothing was written.")
say("=" * 72)
