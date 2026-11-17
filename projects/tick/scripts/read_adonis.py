#!/usr/bin/env python3
"""read_adonis.py

I print the result table from a QIIME 2 adonis visualization, in about ten
lines, so nobody has to unzip a .qzv and page through the metadata that ships
inside it.

WHY
---
A .qzv is a zip. `qiime diversity adonis` puts its result in data/index.html
and also bundles the full sample metadata, which is thousands of lines. An
ordinary unzip-and-cat gives you the metadata and buries the eight rows that
matter.

WHAT I PRINT
------------
One block per visualization: the term, degrees of freedom, sum of squares, R2
and p, exactly as the visualization states them, plus the file's date so a
stale result is visible as stale.

I write nothing and change nothing.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project

    # one file
    python3 scripts/read_adonis.py results/diversity_rebuild_20260706_2045/adonis_bray_batch_state.qzv

    # or every adonis visualization under results/, newest first
    python3 scripts/read_adonis.py
"""
import glob
import html
import os
import shutil
import sys
import tempfile
import time
import zipfile
from html.parser import HTMLParser


class TableReader(HTMLParser):
    """I collect the cells of every table in an HTML document, row by row."""

    def __init__(self):
        HTMLParser.__init__(self)
        self.tables = []
        self._rows = None
        self._cells = None
        self._text = None

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
            if any(cell for cell in self._cells):
                self._rows.append(self._cells)
            self._cells = None
        elif tag == "table" and self._rows is not None:
            if self._rows:
                self.tables.append(self._rows)
            self._rows = None


def tables_in(path):
    """I return the tables inside one .qzv, without unpacking it where you live."""
    workdir = tempfile.mkdtemp(prefix="adonis_")
    try:
        with zipfile.ZipFile(path) as archive:
            wanted = [name for name in archive.namelist()
                      if name.endswith("index.html")]
            if not wanted:
                return None, "no index.html inside"
            archive.extract(wanted[0], workdir)
            with open(os.path.join(workdir, wanted[0]),
                      encoding="utf-8", errors="replace") as handle:
                document = handle.read()
    except (zipfile.BadZipFile, OSError) as problem:
        return None, str(problem)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    reader = TableReader()
    reader.feed(document)
    return reader.tables, None


def looks_like_adonis(rows):
    """I recognise the adonis term table by its column names."""
    header = " ".join(rows[0]).lower()
    return "r2" in header or "sumofsqs" in header or "df" in header


def show(path):
    """I print one visualization's result table."""
    stamp = time.strftime("%Y-%m-%d", time.localtime(os.path.getmtime(path)))
    print("=" * 72)
    print("%s    (written %s)" % (path, stamp))
    print("=" * 72)
    tables, problem = tables_in(path)
    if tables is None:
        print("  could not read: %s" % problem)
        print()
        return
    chosen = [rows for rows in tables if looks_like_adonis(rows)]
    if not chosen:
        print("  no table with Df / SumOfSqs / R2 columns. Tables found: %d"
              % len(tables))
        for rows in tables[:2]:
            print("    first row: %s" % " | ".join(rows[0])[:70])
        print()
        return
    for rows in chosen:
        widths = []
        for index in range(max(len(row) for row in rows)):
            widths.append(max(len(row[index]) if index < len(row) else 0
                              for row in rows))
        for row in rows:
            cells = []
            for index, cell in enumerate(row):
                cells.append(cell.ljust(widths[index]) if index < len(row) - 1
                             else cell)
            print(("  " + "  ".join(cells)).rstrip())
        print()


def main():
    paths = sys.argv[1:]
    if not paths:
        paths = sorted(glob.glob("results/**/adonis*.qzv", recursive=True),
                       key=os.path.getmtime, reverse=True)
        if not paths:
            sys.exit("no adonis*.qzv found under results/. Give me a path.")
        print("found %d adonis visualization(s), newest first\n" % len(paths))
    for path in paths:
        if not os.path.isfile(path):
            print("not found: %s\n" % path)
            continue
        show(path)
    print("R2 is the share of variance the term explains. In a sequential")
    print("model the order of terms matters: a term is fitted after the ones")
    print("above it, so its R2 is what it adds, not what it would explain")
    print("alone.")


if __name__ == "__main__":
    main()
