#!/usr/bin/env python3
"""I report what a copy of the tick metadata workbook actually holds.

Two copies of this workbook exist with the same 22 sheets and different
checksums, so one is a later save. Rather than guess which, I count the
fields that matter for the ked and accessioning work and print them, so the
same script run against both copies says whether they are equivalent for
this purpose.

I read only. Pass the workbook path as the one argument.
"""
import os
import re
import sys

import pandas as pd

BATCHES = [("2025-01-24 Batch 07", "KED"), ("2025-02-12 Batch 08", '"ked" count'),
           ("2025-02-13 Batch 09", '"ked" count'), ("2025-02-19 Batch 10", '"ked" count')]


def tidy(frame):
    frame.columns = [str(c).strip() for c in frame.columns]
    return frame


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: check_workbook_vintage.py <path to workbook>")
    path = sys.argv[1]
    if not os.path.isfile(path):
        sys.exit("not found: " + path)
    try:
        book = pd.ExcelFile(path)
    except (ValueError, OSError, ImportError) as problem:
        sys.exit("cannot open %s (%s)" % (path, problem))
    print("file   %s" % path)
    print("sheets %d" % len(book.sheet_names))

    meta = tidy(book.parse("Metadata", dtype=str))
    notes = meta["Accessioning Notes"].dropna().astype(str)
    notes = notes[notes.str.strip().str.lower() != "nan"]
    print()
    print("Metadata sheet: rows %d" % len(meta))
    print("   Accessioning Notes filled      %4d" % len(notes))
    for label, pattern in (("mentioning a ked", r"\bkeds?\b"),
                           ("mentioning engorged", r"engorged"),
                           ("mentioning males", r"\bmales?\b"),
                           ("mentioning a spider", r"spiders?")):
        print("   %-30s %4d" % (label,
              sum(1 for t in notes if re.search(pattern, t, re.I))))
    pairs = 0
    for t in notes:
        if re.search(r"\b\d+\s+(?:fully |partially )?engorged\b[^0-9]{0,12}\d+\s+males?\b",
                     t, re.I):
            pairs += 1
    print("   with an engorged and male count %4d" % pairs)

    print()
    print("extraction bench sheets that carry a ked column")
    for sheet, column in BATCHES:
        if sheet not in book.sheet_names:
            print("   %-10s SHEET MISSING" % sheet[-8:])
            continue
        frame = tidy(book.parse(sheet, dtype=str))
        name = column.strip()
        if name not in frame.columns:
            print("   %-10s no column %r; columns: %s"
                  % (sheet[-8:], name, ", ".join(frame.columns[:6])))
            continue
        values = pd.to_numeric(frame[name], errors="coerce").dropna()
        print("   %-10s ked column filled %3d, of which above zero %3d"
              % (sheet[-8:], len(values), int((values > 0).sum())))
        notes_col = "Additional Notes"
        if notes_col in frame.columns:
            text = frame[notes_col].dropna().astype(str)
            print("   %-10s notes %3d, of which KEDS ONLY %3d"
                  % ("", len(text),
                     sum(1 for t in text if re.search(r"keds?\s*only", t, re.I))))


if __name__ == "__main__":
    main()
