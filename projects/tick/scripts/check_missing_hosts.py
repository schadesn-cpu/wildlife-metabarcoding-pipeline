#!/usr/bin/env python3
"""check_missing_hosts.py

I list every animal with no usable host and show everything in the project
that might still say what it was. I write nothing.

WHY
---
Table 1 puts 35 of 563 animals in "Not recorded" or "Other or unrecorded",
which is 6% of the sample and the paper's main stratification. Before that
becomes a limitation in the manuscript it is worth asking whether the host is
recorded somewhere I did not look. Table 1 read host_bin only, and there are
at least four other places the answer could be.

WHERE I LOOK
------------
    crosswalk host_spp        a second host column Table 1 ignored
    crosswalk host_bin        what Table 1 used
    Metadata 'Host Spp'       the collection sheet
    Metadata 'Host Category'
    the identifiers           ALHT-Deer-102 and ALHT WTD_23 name the host in
                              the barcode itself

I scan the identifiers with the same word lists check_qpcr_by_host uses for
grouping, so a hit here means the existing grouper would have recognised it
had it been given that text.

WHAT I DO NOT DO
----------------
Change anything. If the host turns out to be recoverable, that is a decision
about the crosswalk, not something to bolt onto a table script.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_missing_hosts.py
"""
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
BOOK = "metadata/tick_full_metadata.xlsx"
SHEETS = ("Metadata", "All Run Data - w Metadata")


def grouper():
    try:
        import check_qpcr_by_host as module
    except ImportError as problem:
        sys.exit("cannot import check_qpcr_by_host (%s)." % problem)
    return module


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s)." % problem)
    return tick_lookup.load(), tick_lookup


def resolve(table, module, value):
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    hit = module.resolve(table, "specimen", text)
    return "" if hit is None else str(hit)


def shown(value):
    text = str(value).strip()
    return "" if not text or text.lower() == "nan" else text


def hint(module, text):
    """I return the host group a piece of free text names, or ''."""
    lowered = str(text).lower()
    for entry in module.HOST_GROUPS:
        name, words = entry[0], entry[1]
        for word in (words if isinstance(words, (list, tuple)) else [words]):
            if re.search(r"(?<![a-z])%s" % re.escape(str(word).lower()),
                         lowered):
                return name
    return ""


def main():
    module = grouper()
    table, finder = lookup()
    for path in (CROSSWALK, BOOK):
        if not os.path.isfile(path):
            sys.exit("missing %s. Run me from the project root." % path)

    crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    crosswalk["animal"] = [resolve(table, finder, b) or resolve(table, finder, s)
                           for b, s in zip(crosswalk["barcode"],
                                           crosswalk["sample_id"])]
    crosswalk["group"] = [module.host_group(v) for v in crosswalk["host_bin"]]
    unusable = crosswalk[[g in ("", module.UNMAPPED)
                          for g in crosswalk["group"]]]
    print("animals with no usable host in host_bin: %d of %d"
          % (len(unusable), len(crosswalk)))
    print()

    book = pd.ExcelFile(BOOK)
    sheet_host = {}
    for sheet in SHEETS:
        if sheet not in book.sheet_names:
            continue
        frame = pd.read_excel(book, sheet_name=sheet, dtype=str)
        columns = {str(c).strip().lower(): c for c in frame.columns}
        keys = [columns[k] for k in ("host spp", "host category", "host_spp")
                if k in columns]
        if not keys or "barcode" not in columns:
            continue
        for _, row in frame.iterrows():
            animal = resolve(table, finder, row[columns["barcode"]])
            if not animal:
                continue
            for key in keys:
                value = shown(row[key])
                if value:
                    sheet_host.setdefault(animal, []).append(
                        "%s=%s" % (str(key).strip(), value))

    print("%-22s %-16s %-16s %-28s %s"
          % ("animal", "host_spp", "host_bin", "collection sheet",
             "named in the id"))
    print("-" * 110)
    recoverable = {}
    for _, row in unusable.iterrows():
        animal = row["animal"]
        spp = shown(row.get("host_spp"))
        binned = shown(row.get("host_bin"))
        sheet = "; ".join(sheet_host.get(animal, []))[:28]
        guess = (hint(module, spp) or hint(module, sheet)
                 or hint(module, row.get("sample_id"))
                 or hint(module, row.get("barcode")))
        if guess:
            recoverable[animal] = guess
        print("%-22s %-16s %-16s %-28s %s"
              % (str(animal)[:22], spp[:16], binned[:16], sheet, guess))
    print()

    print("=" * 70)
    print("WHAT THIS MEANS")
    print("=" * 70)
    print("  animals with no usable host_bin:      %d" % len(unusable))
    print("  of those, a host is named elsewhere:  %d" % len(recoverable))
    print("  genuinely unknown:                    %d"
          % (len(unusable) - len(recoverable)))
    if recoverable:
        counted = {}
        for name in recoverable.values():
            counted[name] = counted.get(name, 0) + 1
        print()
        print("  recoverable as: %s"
              % ", ".join("%s %d" % (k, v)
                          for k, v in sorted(counted.items(),
                                             key=lambda kv: -kv[1])))
        print()
        print("  A host named in a barcode is weaker evidence than one on the")
        print("  collection sheet. Both are printed above so you can judge")
        print("  each row rather than take a count from me.")
    print()


if __name__ == "__main__":
    main()
