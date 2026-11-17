#!/usr/bin/env python3
"""check_methods_gaps.py

I close the Methods gaps that the files on ron can settle, and print a short
answer for each. I write nothing and change nothing.

WHAT I CHECK
------------
1. READ-THRESHOLD SENSITIVITY. Section 2.6 calls a pathogen present at 1 read
   and promises a sensitivity check at 2, 5 and 10 reads. I recompute the
   animal-level prevalence (one library per animal, 2,000 draws, seed 20260918,
   the same scheme the paper uses) and the Anaplasma x Borreliella odds ratio
   at each threshold. This needs per-pathogen READ COUNTS, not 0/1 calls. If
   the analysis table carries only calls I say so and stop that section.

2. PERMDISP. Section 2.6 says dispersion was checked. I look for any PERMDISP
   visualization under results/ and print its result. If there is none, the
   sentence has to come out of the manuscript or the test has to be run.

3. KEDS. The notes flag 29 samples with keds; the manuscript accounts for 28
   (21 with zero ked reads, 7 with reads). I list every noted sample and say
   which one is missing from the 16S data, and why.

4. MISSING TICK COUNT AND HOST. The manuscript says 35 of 563 samples have no
   tick count and 19 no host. Those were counted before the join fixes. I
   recount.

5. COLLECTION DATES. Section 2.1 has the date range PENDING. If the metadata
   has a date column I print the range per state.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_methods_gaps.py
"""
import glob
import html
import os
import re
import shutil
import sys
import tempfile
import zipfile
from html.parser import HTMLParser

import numpy as np
import pandas as pd

TABLE = "results/analysis/analysis_table_16S.tsv"
PATHOGENS = ["Anaplasma", "Borreliella", "Borrelia", "Ehrlichia", "Bartonella"]
THRESHOLDS = [1, 2, 5, 10]
DRAWS = 2000
SEED = 20260918
EXPECTED_SPECIMENS = 563
EXPECTED_NO_TICKS = 35
EXPECTED_NO_HOST = 19
EXPECTED_KED_NOTES = 29
KED_WORDS = re.compile(r"\bked|keds\b|hippobosc|lipoptena", re.IGNORECASE)
KED_READ_HINTS = ["ked", "lipoptena", "hippobosc"]
RUN_SUFFIX = re.compile(r"\..*$")


def blank(value):
    """I return '' for anything meaning 'no value', and the text otherwise."""
    text = str(value).strip()
    return "" if text.lower() in ("nan", "none", "<na>", "na", "") else text


def heading(text):
    print("=" * 72)
    print(text)
    print("=" * 72)


def load_lookup():
    """I return (table, module) for the library lookup, or (None, None)."""
    try:
        import tick_lookup
    except ImportError as problem:
        print("  could not import tick_lookup (%s)." % problem)
        return None, None
    if not hasattr(tick_lookup, "load") or not hasattr(tick_lookup, "resolve"):
        print("  tick_lookup has no load()/resolve().")
        return None, None
    return tick_lookup.load(), tick_lookup


def resolve(table, module, value):
    """I resolve one library or specimen name to its animal, or ''."""
    stem = RUN_SUFFIX.sub("", str(value).strip())
    for marker, candidate in (("16S", value), ("16S", stem),
                              ("specimen", stem), ("specimen", value)):
        hit = module.resolve(table, marker, candidate)
        if hit is not None:
            return blank(hit)
    return ""


# ---------------------------------------------------------------- 1 sensitivity

def read_count_columns(frame):
    """I map each pathogen to a column of read counts, where one exists."""
    found = {}
    for pathogen in PATHOGENS:
        for column in frame.columns:
            name = str(column)
            if not name.lower().startswith(pathogen.lower()):
                continue
            if name != pathogen and "read" not in name.lower() and "count" not in name.lower():
                continue
            values = pd.to_numeric(frame[column], errors="coerce").dropna()
            if len(values) and values.max() > 1:
                found[pathogen] = column
                break
    return found


def odds_ratio(x, y):
    """I return the Haldane-corrected odds ratio of two boolean arrays."""
    a = np.sum(x & y) + 0.5
    b = np.sum(x & ~y) + 0.5
    c = np.sum(~x & y) + 0.5
    d = np.sum(~x & ~y) + 0.5
    return (a * d) / (b * c)


def sensitivity(frame, animals):
    heading("1. READ-THRESHOLD SENSITIVITY (one library per animal, %d draws)"
            % DRAWS)
    counts = read_count_columns(frame)
    if not counts:
        print("  %s has no per-pathogen read-count columns, only calls." % TABLE)
        print("  Columns that name a pathogen:")
        named = [c for c in frame.columns
                 if any(p.lower() in str(c).lower() for p in PATHOGENS)]
        print("    " + (", ".join(named[:12]) or "none"))
        print("  I need the genus-level count table this table was built from.")
        print("  Paste the output of:")
        print("    grep -n 'read_csv\\|_pos\\|MIN_READS' "
              "scripts/build_analysis_table.py | head -30")
        print()
        return
    print("  read-count columns: %s" % ", ".join(
        "%s=%s" % (p, c) for p, c in counts.items()))
    work = frame.assign(_animal=animals)
    work = work[work["_animal"] != ""].reset_index(drop=True)
    groups = list(work.groupby("_animal").indices.values())
    generator = np.random.default_rng(SEED)
    picks = np.array([[generator.choice(g) for g in groups]
                      for _ in range(DRAWS)])
    reads = {p: pd.to_numeric(work[c], errors="coerce").fillna(0).values
             for p, c in counts.items()}
    print("  %d libraries, %d animals" % (len(work), len(groups)))
    print()
    header = "  %-8s" % "reads>="
    for pathogen in counts:
        header += " %12s" % pathogen[:12]
    pair = [p for p in ("Anaplasma", "Borreliella") if p in counts]
    if len(pair) == 2:
        header += "   A x B OR"
    print(header)
    for threshold in THRESHOLDS:
        line = "  %-8d" % threshold
        for pathogen in counts:
            called = reads[pathogen] >= threshold
            shares = called[picks].mean(axis=1) * 100
            line += " %11.1f%%" % np.median(shares)
        if len(pair) == 2:
            first = reads[pair[0]] >= threshold
            second = reads[pair[1]] >= threshold
            ratios = [odds_ratio(first[row], second[row]) for row in picks]
            low, mid, high = np.percentile(ratios, [2.5, 50, 97.5])
            line += "   %.2f (%.2f-%.2f)" % (mid, low, high)
        print(line)
    print()
    print("  Threshold 1 should reproduce the paper (Anaplasma 38.5%,")
    print("  Borreliella 35.6%, Borrelia 3.6%, Bartonella 13.1%; A x B 1.96).")
    print("  If it does not, the columns above are not the ones the calls")
    print("  came from, and the rest of this table should not be used.")
    print()


# ---------------------------------------------------------------- 2 permdisp

class TableReader(HTMLParser):
    """I collect the cells of every table in an HTML document."""

    def __init__(self):
        HTMLParser.__init__(self)
        self.tables, self._rows, self._cells, self._text = [], None, None, None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self._rows = []
        elif tag == "tr" and self._rows is not None:
            self._cells = []
        elif tag in ("td", "th") and self._cells is not None:
            self._text = []

    def handle_data(self, data):
        if self._text is not None:
            self._text.append(data)

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self._text is not None:
            self._cells.append(html.unescape("".join(self._text)).strip())
            self._text = None
        elif tag == "tr" and self._cells is not None:
            if any(self._cells):
                self._rows.append(self._cells)
            self._cells = None
        elif tag == "table" and self._rows is not None:
            if self._rows:
                self.tables.append(self._rows)
            self._rows = None


def first_table(path):
    """I return the first short table in a .qzv's index.html, or None."""
    workdir = tempfile.mkdtemp(prefix="permdisp_")
    try:
        with zipfile.ZipFile(path) as archive:
            wanted = [n for n in archive.namelist() if n.endswith("index.html")]
            if not wanted:
                return None
            archive.extract(wanted[0], workdir)
            with open(os.path.join(workdir, wanted[0]), encoding="utf-8",
                      errors="replace") as handle:
                document = handle.read()
    except (zipfile.BadZipFile, OSError):
        return None
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    reader = TableReader()
    reader.feed(document)
    for rows in reader.tables:
        if len(rows) <= 12:
            return rows
    return None


def permdisp():
    heading("2. WAS PERMDISP RUN?")
    paths = sorted(p for p in glob.glob("results/**/*.qzv", recursive=True)
                   if "permdisp" in os.path.basename(p).lower()
                   or "disp" in os.path.basename(p).lower())
    logs = sorted(p for p in glob.glob("results/**/*", recursive=True)
                  if os.path.isfile(p) and p.endswith((".sh", ".log", ".txt"))
                  and os.path.getsize(p) < 2000000)
    mentions = []
    for path in logs + glob.glob("scripts/*.sh") + glob.glob("scripts/*.py"):
        if os.path.basename(path) == os.path.basename(__file__):
            continue
        try:
            with open(path, encoding="utf-8", errors="replace") as handle:
                if "permdisp" in handle.read().lower():
                    mentions.append(path)
        except OSError:
            continue
    if not paths:
        print("  no PERMDISP visualization under results/.")
    for path in paths:
        rows = first_table(path)
        print("  %s" % path)
        if rows is None:
            print("    (no readable result table)")
            continue
        for row in rows:
            print("    " + " | ".join(row))
    print("  files that mention permdisp: %s"
          % (", ".join(sorted(set(mentions))[:6]) or "none"))
    if not paths:
        print("  VERDICT: no evidence it was run. Either run it or take the")
        print("  sentence out of section 2.6.")
    print()


# ---------------------------------------------------------------- 3,4,5 metadata

def metadata_frames():
    """I return (path, frame) for every readable table under metadata/."""
    found = []
    patterns = ["metadata/**/*.tsv", "metadata/**/*.csv", "metadata/**/*.txt",
                "metadata/**/*.xlsx"]
    for pattern in patterns:
        for path in sorted(glob.glob(pattern, recursive=True)):
            if ".bak" in path or "/_archive" in path:
                continue
            try:
                if path.endswith(".xlsx"):
                    sheets = pd.read_excel(path, sheet_name=None, dtype=str)
                    for name, frame in sheets.items():
                        found.append(("%s[%s]" % (path, name), frame))
                else:
                    sep = "," if path.endswith(".csv") else "\t"
                    found.append((path, pd.read_csv(path, sep=sep, dtype=str)))
            except (ValueError, OSError, ImportError,
                    pd.errors.ParserError) as problem:
                print("  skipped %s (%s)" % (path, str(problem)[:50]))
    return found


def keds(frames, frame, animals, table, module):
    heading("3. THE 29 KED NOTES: WHICH ONE IS MISSING?")
    noted = {}
    for path, meta in frames:
        text_columns = [c for c in meta.columns
                        if any(w in str(c).lower()
                               for w in ("note", "comment", "remark"))]
        id_columns = [c for c in meta.columns
                      if any(w in str(c).lower()
                             for w in ("barcode", "sample", "specimen", "id"))]
        if not text_columns or not id_columns:
            continue
        for _, row in meta.iterrows():
            text = " ".join(blank(row[c]) for c in text_columns)
            if not KED_WORDS.search(text):
                continue
            name = blank(row[id_columns[0]])
            if name and name not in noted:
                noted[name] = (path, text[:40])
    print("  samples with a ked note: %d (manuscript says %d)"
          % (len(noted), EXPECTED_KED_NOTES))
    if not noted:
        print("  no notes column mentions keds. Where did 29 come from?")
        print()
        return
    ked_columns = [c for c in frame.columns
                   if any(h in str(c).lower() for h in KED_READ_HINTS)]
    ked_reads = None
    if ked_columns:
        ked_reads = pd.to_numeric(frame[ked_columns[0]],
                                  errors="coerce").fillna(0)
        print("  ked read column in the 16S table: %s" % ked_columns[0])
    else:
        print("  no ked read column in %s" % TABLE)
    by_animal = {}
    for index, animal in enumerate(animals):
        if animal:
            value = float(ked_reads.iloc[index]) if ked_reads is not None else 0
            by_animal[animal] = max(by_animal.get(animal, 0), value)
    zero, some, absent = [], [], []
    for name in sorted(noted):
        animal = resolve(table, module, name) or name
        if animal not in by_animal:
            absent.append((name, animal))
        elif by_animal[animal] > 0:
            some.append(name)
        else:
            zero.append(name)
    print("  zero ked reads: %d    ked reads: %d    not in 16S data: %d"
          % (len(zero), len(some), len(absent)))
    for name, animal in absent:
        print("    NOT IN 16S: %s (resolves to %s)  note: %s"
              % (name, animal, noted[name][1]))
    print()


def missing_fields(frames):
    heading("4. SAMPLES WITH NO TICK COUNT OR NO HOST")
    best = None
    for path, meta in frames:
        ticks = [c for c in meta.columns
                 if str(c).lower().replace("_", "") in ("numticks", "ticks",
                                                        "tickcount")]
        host = [c for c in meta.columns
                if str(c).lower() in ("host", "host_spp", "hostspecies",
                                      "host_species", "species")]
        ids = [c for c in meta.columns if "barcode" in str(c).lower()]
        if ticks and host and ids:
            if best is None or len(meta) > len(best[1]):
                best = (path, meta, ids[0], ticks[0], host[0])
    if best is None:
        print("  no metadata table with barcode, tick-count and host columns.")
        print()
        return
    path, meta, id_column, tick_column, host_column = best
    specimens = meta[[id_column, tick_column, host_column]].copy()
    specimens[id_column] = [blank(v) for v in specimens[id_column]]
    specimens = specimens[specimens[id_column] != ""]
    specimens = specimens.drop_duplicates(subset=id_column)
    no_ticks = sum(1 for v in specimens[tick_column]
                   if pd.isna(pd.to_numeric(blank(v) or "x", errors="coerce")))
    no_host = sum(1 for v in specimens[host_column] if not blank(v))
    print("  table: %s" % path)
    print("  specimens: %d (manuscript %d)" % (len(specimens),
                                               EXPECTED_SPECIMENS))
    print("  no tick count: %d (manuscript %d)" % (no_ticks, EXPECTED_NO_TICKS))
    print("  no host:       %d (manuscript %d)" % (no_host, EXPECTED_NO_HOST))
    if len(specimens) != EXPECTED_SPECIMENS:
        print("  The specimen count differs, so this may not be the table the")
        print("  crosswalk used. Treat the two counts as provisional.")
    print()


def dates(frames):
    heading("5. COLLECTION DATE RANGE")
    for path, meta in frames:
        date_columns = [c for c in meta.columns
                        if "date" in str(c).lower()
                        and "extract" not in str(c).lower()
                        and "run" not in str(c).lower()]
        if not date_columns:
            continue
        column = date_columns[0]
        parsed = pd.to_datetime(meta[column], errors="coerce")
        if parsed.notna().sum() < 20:
            continue
        state = [c for c in meta.columns if str(c).lower() == "state"]
        print("  %s  column '%s'" % (path, column))
        print("    all: %s to %s (%d dated, %d not)"
              % (parsed.min().date(), parsed.max().date(),
                 parsed.notna().sum(), parsed.isna().sum()))
        if state:
            for value, rows in meta.groupby(meta[state[0]].fillna("?")):
                block = parsed[rows.index]
                if block.notna().any():
                    print("    %-4s %s to %s" % (str(value)[:4],
                                                 block.min().date(),
                                                 block.max().date()))
        print()
        return
    print("  no metadata table has a usable collection date column.")
    print()


def main():
    if not os.path.isfile(TABLE):
        sys.exit("cannot find %s. Run this from "
                 "/home/users/sne24/meedlab/tick_project." % TABLE)
    frame = pd.read_csv(TABLE, sep="\t", dtype=str)
    if "sample" not in frame.columns:
        sys.exit("%s has no 'sample' column." % TABLE)
    table, module = load_lookup()
    if module is None:
        sys.exit("the library lookup is needed. Build it with "
                 "python3 scripts/build_4marker_crosswalk.py")
    animals = [resolve(table, module, s) for s in frame["sample"]]
    sensitivity(frame, animals)
    permdisp()
    frames = metadata_frames()
    if not frames:
        print("no readable tables under metadata/, so sections 3 to 5 are "
              "skipped.")
        return
    keds(frames, frame, animals, table, module)
    missing_fields(frames)
    dates(frames)


if __name__ == "__main__":
    main()
