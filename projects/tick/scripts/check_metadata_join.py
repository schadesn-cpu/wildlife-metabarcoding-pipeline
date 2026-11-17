#!/usr/bin/env python3
"""check_metadata_join.py

I try several ways of joining the crosswalk to the Metadata sheet and report
how many animals each one resolves. I write nothing.

WHY
---
check_town_to_county.py matched 501 of 563 animals by comparing the barcode
strings directly. Every other script in this pipeline goes through
tick_lookup.resolve instead, which exists because the barcodes are not
consistent between sources. 62 animals with no collection record at all is
more likely to be my join than a gap in the data, so I test the join before
anyone goes looking for missing paperwork.

WHAT I TRY
----------
    raw            barcode against barcode, as written
    stripped       leading and trailing whitespace removed
    upper          case folded as well
    alphanumeric   punctuation and spaces removed as well
    no trailing .0 for barcodes Excel turned into floats
    tick_lookup    both sides resolved to an animal first

The last one is the one the rest of the pipeline uses. If it wins, the fix is
to use it and the 62 disappear.

I also print unmatched examples from both sides, because the shape of the
difference usually says what happened.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_metadata_join.py
"""
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
BOOK = "metadata/tick_full_metadata.xlsx"
SHEETS = ("Metadata", "All Run Data - w Metadata")
SHOW = 12


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s)." % problem)
    return tick_lookup.load(), tick_lookup


def raw(value):
    return str(value)


def stripped(value):
    return str(value).strip()


def upper(value):
    return str(value).strip().upper()


def alphanumeric(value):
    return re.sub(r"[^A-Z0-9]", "", str(value).strip().upper())


def no_float(value):
    text = str(value).strip()
    if re.match(r"^\d+\.0$", text):
        text = text[:-2]
    return re.sub(r"[^A-Z0-9]", "", text.upper())


RULES = [("raw", raw), ("stripped", stripped), ("upper", upper),
         ("alphanumeric", alphanumeric), ("no trailing .0", no_float)]


def barcode_column(frame):
    """I find the barcode column whatever spaces or case it was saved with."""
    for column in frame.columns:
        if str(column).strip().lower() == "barcode":
            return column
    return None


def main():
    if not os.path.isfile(CROSSWALK):
        sys.exit("missing %s. Run me from the project root." % CROSSWALK)
    if not os.path.isfile(BOOK):
        sys.exit("missing %s." % BOOK)
    crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    print("crosswalk animals: %d" % len(crosswalk))
    print()

    book = pd.ExcelFile(BOOK)
    for sheet in SHEETS:
        if sheet not in book.sheet_names:
            print("no sheet named %r, skipping" % sheet)
            continue
        frame = pd.read_excel(book, sheet_name=sheet, dtype=str)
        column = barcode_column(frame)
        if column is None:
            print("%r has no barcode column, skipping" % sheet)
            continue
        print("=" * 72)
        print("%s: %d rows, barcode column %r" % (sheet, len(frame), column))
        print("=" * 72)

        for name, rule in RULES:
            theirs = set(rule(v) for v in frame[column].dropna())
            hits = sum(1 for v in crosswalk["barcode"]
                       if rule(v) in theirs)
            print("  %-16s %4d of %d" % (name, hits, len(crosswalk)))

        theirs = set(no_float(v) for v in frame[column].dropna())
        missing = [str(v) for v in crosswalk["barcode"]
                   if no_float(v) not in theirs]
        if missing:
            print()
            print("  crosswalk barcodes with no match (%d), first %d:"
                  % (len(missing), SHOW))
            print("    %s" % ", ".join(missing[:SHOW]))
            spare = [str(v) for v in frame[column].dropna()
                     if no_float(v) not in set(no_float(b) for b
                                               in crosswalk["barcode"])]
            print("  sheet barcodes not in the crosswalk (%d), first %d:"
                  % (len(spare), SHOW))
            print("    %s" % ", ".join(spare[:SHOW]))
        print()

    table, module = lookup()
    print("=" * 72)
    print("THROUGH tick_lookup, WHICH IS WHAT THE PIPELINE USES")
    print("=" * 72)
    animals = set()
    for value in crosswalk["barcode"]:
        hit = module.resolve(table, "specimen", str(value).strip())
        if hit is not None:
            animals.add(str(hit))
    by_id = set()
    for value in crosswalk["sample_id"]:
        hit = module.resolve(table, "specimen", str(value).strip())
        if hit is not None:
            by_id.add(str(hit))
    print("  crosswalk barcode  resolves to an animal: %d of %d"
          % (len(animals), len(crosswalk)))
    print("  crosswalk sample_id resolves to an animal: %d of %d"
          % (len(by_id), len(crosswalk)))

    for sheet in SHEETS:
        if sheet not in book.sheet_names:
            continue
        frame = pd.read_excel(book, sheet_name=sheet, dtype=str)
        column = barcode_column(frame)
        if column is None:
            continue
        theirs = set()
        for value in frame[column].dropna():
            hit = module.resolve(table, "specimen", str(value).strip())
            if hit is not None:
                theirs.add(str(hit))
        shared = len(animals & theirs) if animals else 0
        print("  %-28s resolves %d animals, %d shared with the crosswalk"
              % (sheet, len(theirs), shared))
    print()
    print("  If the shared count reaches 563, the join is tick_lookup and the")
    print("  62 were never missing. If it does not, the gap is real and the")
    print("  unmatched barcodes above say what kind of gap it is.")
    print()


if __name__ == "__main__":
    main()
