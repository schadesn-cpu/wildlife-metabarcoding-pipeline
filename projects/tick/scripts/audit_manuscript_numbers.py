#!/usr/bin/env python3
"""audit_manuscript_numbers.py  --  READ ONLY. I open files for reading and
nothing else. No writes, no renames, no temp files.

WHY THIS EXISTS RATHER THAN THE EARLIER ONE
-------------------------------------------
The earlier version read the USGS package, which is a copy frozen on 23
October and since edited by hand. Numbers for the manuscript have to come
from the live project, and the live project is here.

WHAT IT SETTLES
---------------
    Section 2.2   "tick size or engorgement for 310 of 561 specimens and are
                   silent for the remaining 253"              310 + 253 = 563
    Section 3.4   "the detection appeared in 13 of their 18 tested
                   subsamples against 263 of 807 elsewhere (33%)"
    Section 3.4   "All 7 black bear with a qPCR result"   against Table 1's 10
    Section 3.4   "214 of 559" and "207 of 552"           already resolved,
                                                          recomputed here as
                                                          a cross-check

IT NAMES THE FILE IT READ
-------------------------
This project holds more than one generation of several tables, and
check_which_verified_table.py exists because of it. So part 0 lists EVERY
copy of each input it can find, with size and date, and refuses to continue
when there is more than one unless you name which. It will not quietly pick
the newest and hand you a number built on the wrong generation.

    cd ~/meedlab/tick_project
    python3 scripts/audit_manuscript_numbers.py
    python3 scripts/audit_manuscript_numbers.py --recount results/analysis/qpcr_recount_by_rule.tsv
"""
import csv
import os
import re
import sys
from collections import Counter, defaultdict

ROOT = os.path.expanduser(os.environ.get("AUDIT_ROOT","~/meedlab/tick_project"))
SKIP = {".git", "work", "node_modules", "__pycache__", ".snakemake", "tmp"}
WANTED = {
    "recount": "qpcr_recount_by_rule.tsv",
    "covar": "accessioning_covariates.tsv",
    "book": "tick_full_metadata.xlsx",
}
bar = "=" * 78


def argval(flag):
    return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else None


def locate():
    hits = defaultdict(list)
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in SKIP
                       and not d.startswith(".")]
        for name in filenames:
            for tag, want in WANTED.items():
                if name == want:
                    hits[tag].append(os.path.join(dirpath, name))
    return hits


print(bar)
print("0. WHICH FILE AM I ACTUALLY READING")
print(bar)
if not os.path.isdir(ROOT):
    sys.exit("no such directory: %s" % ROOT)
hits = locate()
chosen, stop = {}, False
for tag, want in WANTED.items():
    override = argval("--" + tag)
    found = hits.get(tag, [])
    print("\n  %s" % want)
    if override:
        print("     named on the command line: %s" % override)
        if not os.path.isfile(override):
            sys.exit("     ...which does not exist.")
        chosen[tag] = override
        continue
    if not found:
        print("     NOT FOUND anywhere under %s" % ROOT)
        continue
    for p in sorted(found):
        st = os.stat(p)
        import datetime
        when = datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M")
        print("     %10d bytes  %s  %s" % (st.st_size, when,
                                           os.path.relpath(p, ROOT)))
    if len(found) == 1:
        chosen[tag] = found[0]
    else:
        print("     *** %d copies. Name one with --%s <path> and run again. ***"
              % (len(found), tag))
        stop = True
if stop:
    print("\nStopping rather than guessing which generation is canonical.")
    sys.exit(1)
print()

found_values = {}


def note(label, value):
    found_values.setdefault(value, label)


def key(animal):
    """Collapse the two duplicate bear accessions, leave everything else alone."""
    a = (animal or "").strip()
    m = re.search(r'bear\s*([a-e])(?![a-z])', a, re.I)
    return "BEAR_" + m.group(1).upper() if m else a


def num(v):
    try:
        return float((v or "").strip() or 0)
    except ValueError:
        return 0.0


def tsv(path):
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------- 1. animals and bears
print(bar)
print("1. ANIMALS, BEARS, AND WHAT THE MERGE CHANGES")
print(bar)
rec = tsv(chosen["recount"])
animals = sorted({r["animal"].strip() for r in rec if r.get("animal")})
groups = defaultdict(list)
for a in animals:
    groups[key(a)].append(a)
merged = {k: v for k, v in groups.items() if len(v) > 1}
print("  accessions with a qPCR result : %d" % len(animals))
print("  animals after the merge       : %d" % len(groups))
note("animals with a qPCR result, pre-merge", len(animals))
note("animals with a qPCR result, post-merge", len(groups))
for k, v in sorted(merged.items()):
    print("     %-10s <- %s" % (k, "  +  ".join(v)))

hostof = {}
for r in rec:
    hostof.setdefault(key(r["animal"]),
                      "%s %s" % (r.get("host_bin", ""), r.get("host_spp", "")))
isbear = {k: bool(re.search(r'bear', v, re.I)) for k, v in hostof.items()}
nb_pre = sum(1 for a in animals if isbear.get(key(a)))
nb_post = sum(1 for k in groups if isbear.get(k))
print("  black bear accessions         : %d" % nb_pre)
print("  black bear animals            : %d" % nb_post)
note("black bear accessions", nb_pre)
note("black bear animals, post-merge", nb_post)
print()

# --------------------------------------- 2. wells, which are subsamples
print(bar)
print("2. TESTED SUBSAMPLES FROM THE WELLS COLUMN")
print(bar)
print("  Section 2.5: each tested subsample occupied eight wells, the seven")
print("  targets and a control. So on one target's rows, wells == subsamples.")
print()
targets = sorted({r["target"].strip() for r in rec if r.get("target")})
RULECOLS = [c for c in rec[0] if c.startswith(("stated_", "recorded_", "narrow_"))]
print("  %-22s %8s %10s %10s %10s" % ("target", "animals", "wells",
                                      "stated+", "recorded+"))
for t in targets:
    rows = [r for r in rec if r.get("target", "").strip() == t]
    w = sum(num(r.get("wells")) for r in rows)
    s = sum(num(r.get("stated_positive")) for r in rows)
    rr = sum(num(r.get("recorded_positive")) for r in rows)
    print("  %-22s %8d %10d %10d %10d" % (t, len(rows), w, s, rr))
    note("total tested subsamples (wells on %s)" % t, int(w))
print()
mic = [r for r in rec if re.search(r'microti', r.get("target", ""), re.I)]
if mic:
    bw = sum(num(r.get("wells")) for r in mic if isbear.get(key(r["animal"])))
    ow = sum(num(r.get("wells")) for r in mic if not isbear.get(key(r["animal"])))
    bs = sum(num(r.get("stated_positive")) for r in mic if isbear.get(key(r["animal"])))
    os_ = sum(num(r.get("stated_positive")) for r in mic if not isbear.get(key(r["animal"])))
    br = sum(num(r.get("recorded_positive")) for r in mic if isbear.get(key(r["animal"])))
    orr = sum(num(r.get("recorded_positive")) for r in mic if not isbear.get(key(r["animal"])))
    print("  B. MICROTI, BEAR AGAINST EVERYTHING ELSE")
    print("  %-30s %10s %10s" % ("", "bear", "other"))
    print("  %-30s %10d %10d" % ("tested subsamples", bw, ow))
    print("  %-30s %10d %10d" % ("positive, stated rule", bs, os_))
    print("  %-30s %10d %10d" % ("positive, recorded rule", br, orr))
    for lbl, v in (("bear subsamples", bw), ("other subsamples", ow),
                   ("bear positive subsamples, stated", bs),
                   ("other positive subsamples, stated", os_)):
        note(lbl, int(v))
    # animal-level, to cross-check the two sentences already replaced
    for rule in ("stated_positive", "recorded_positive", "narrow_positive"):
        if rule not in rec[0]:
            continue
        pos = {key(r["animal"]) for r in mic if num(r.get(rule)) > 0}
        nonbear = {k for k in pos if not isbear.get(k)}
        nonbear_n = sum(1 for k in groups if not isbear.get(k))
        note("B. microti positive animals post-merge (%s)" % rule, len(pos))
        note("B. microti positive non-bear (%s)" % rule, len(nonbear))
        note("non-bear animals post-merge", nonbear_n)
print()

# ------------------------------------- 3. the accessioning notes, looked at
print(bar)
print("3. WHAT THE ACCESSIONING NOTES ACTUALLY SAY")
print(bar)
print("  310 is not the engorged_count column, which covers far fewer animals.")
print("  Section 2.2 says 'accessioning notes describe tick size or")
print("  engorgement', so this prints what is in the note fields instead of")
print("  guessing a keyword list. Read them, then the rule can be written.")
print()
book = chosen.get("book")
if not book:
    print("  tick_full_metadata.xlsx not found; skipping.")
else:
    try:
        import openpyxl
    except ImportError:
        print("  openpyxl is not importable in this environment.")
        print("  Try:  conda activate qiime2-amplicon-2024.5 && pip install openpyxl")
        openpyxl = None
    if openpyxl:
        wb = openpyxl.load_workbook(book, read_only=True, data_only=True)
        for sheet in wb.sheetnames:
            ws = wb[sheet]
            rows_iter = ws.iter_rows(values_only=True)
            try:
                head = [str(c).strip() if c is not None else "" for c in next(rows_iter)]
            except StopIteration:
                continue
            notecols = [i for i, h in enumerate(head)
                        if re.search(r'note|comment|descr|remark|condition|'
                                     r'engorg|size|stage|life', h, re.I)]
            if not notecols:
                continue
            print("  sheet %r, %d candidate columns" % (sheet, len(notecols)))
            vals = defaultdict(list)
            total = 0
            for row in rows_iter:
                total += 1
                for i in notecols:
                    v = row[i] if i < len(row) else None
                    if v is not None and str(v).strip():
                        vals[head[i]].append(str(v).strip())
            for h in [head[i] for i in notecols]:
                got = vals.get(h, [])
                print("\n    %-28s %d of %d rows non-empty" % (h, len(got), total))
                for text, n in Counter(got).most_common(12):
                    print("        %4d  %s" % (n, text[:88]))
                note("non-empty %r" % h, len(got))
        wb.close()
print()

# ----------------------------------------------------------- 4. the verdict
print(bar)
print("4. WHICH PUBLISHED NUMBER IS WHICH")
print(bar)
PUBLISHED = [
    ("Sec 2.2  size or engorgement recorded", 310),
    ("Sec 2.2  silent", 253),
    ("Sec 3.4  B. microti positive animals", 214),
    ("Sec 3.4  animals with a qPCR result", 559),
    ("Sec 3.4  B. microti positive, non-bear", 207),
    ("Sec 3.4  non-bear animals", 552),
    ("Sec 3.4  bear subsamples with the detection", 13),
    ("Sec 3.4  bear subsamples tested", 18),
    ("Sec 3.4  other subsamples with the detection", 263),
    ("Sec 3.4  other subsamples tested", 807),
    ("Sec 3.4  black bear with a qPCR result", 7),
    ("Sec 2.5  subsamples with a positive control", 833),
    ("Table 1  black bear accessions", 10),
]
for label, value in PUBLISHED:
    print("  %-46s %5d   %s"
          % (label, value, found_values.get(value) or "*** reproduces nothing ***"))
print()
print("Send this whole output back.")
