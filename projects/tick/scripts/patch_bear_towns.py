#!/usr/bin/env python3
"""patch_bear_towns.py

Writes the collection towns for two bears into the metadata workbook,
from the collector's own shipment sheet. The town and county values, and
the coordinates supplied with them, are withheld; they are read from
metadata/bear_town_fixes.tsv, which is not in this repository. See
projects/tick/README.md.

The collector's record disagrees with PROJECT_LOG.md for one of the two
animals, by about fifty miles. She collected it herself and supplied
coordinates, so her record is the one used here and the log needs
correcting separately, by hand.

These two animals otherwise have no town, so they sit among the animals
the collection map leaves out. With a town they resolve to a county
through the ordinary town-to-county path and appear on the map.

BEFORE REBUILDING THE FIGURE, read the count this prints. If the batch-10 and
batch-12 accessions for the same bear are being counted as two animals, this
places FOUR rather than TWO, and the caption number has to match whichever is
true. That duplicate question is still open.

Backs the workbook up first. Refuses to overwrite a town that is already set
to something else, and says so rather than silently winning.

    cd ~/meedlab/tick_project
    python3 scripts/patch_bear_towns.py            # show what it would do
    python3 scripts/patch_bear_towns.py --apply    # do it
"""
import os
import shutil
import sys
import datetime

BOOK = "metadata/tick_full_metadata.xlsx"
SHEET = "Metadata"
TOWNS = {}
_FIX = "metadata/bear_town_fixes.tsv"
if os.path.isfile(_FIX):
    with open(_FIX) as _fh:
        next(_fh, None)
        for _ln in _fh:
            _p = _ln.rstrip("\n").split("\t")
            if len(_p) >= 3 and _p[0].strip():
                TOWNS[_p[0].strip()] = (_p[1].strip(), _p[2].strip())
else:
    sys.exit("missing %s. The collection town and county for these specimens "
             "are withheld; see projects/tick/README.md." % _FIX)


def main():
    apply = "--apply" in sys.argv
    if not os.path.isfile(BOOK):
        sys.exit("missing %s. Run me from the project root." % BOOK)
    import openpyxl
    wb = openpyxl.load_workbook(BOOK)
    if SHEET not in wb.sheetnames:
        sys.exit("%s has no sheet %r" % (BOOK, SHEET))
    ws = wb[SHEET]

    head = [str(c.value).strip().lower() if c.value else ""
            for c in next(ws.iter_rows(min_row=1, max_row=1))]
    try:
        bc = head.index("barcode") + 1
        tc = head.index("town") + 1
    except ValueError:
        sys.exit("need 'barcode' and 'town' columns; found %s" % head)

    changes, conflicts, missing = [], [], dict(TOWNS)
    for r in range(2, ws.max_row + 1):
        raw = ws.cell(row=r, column=bc).value
        key = str(raw).strip().lower() if raw else ""
        if key not in TOWNS:
            continue
        missing.pop(key, None)
        town, county = TOWNS[key]
        now = ws.cell(row=r, column=tc).value
        now = str(now).strip() if now else ""
        if not now:
            changes.append((r, raw, town, county))
        elif now.lower() != town.lower():
            conflicts.append((r, raw, now, town))

    for r, bar, town, county in changes:
        print("  row %-5d %-32s town -> %-12s (%s County)" % (r, bar, town, county))
    for r, bar, now, town in conflicts:
        print("  row %-5d %-32s ALREADY %r, not overwriting with %r"
              % (r, bar, now, town))
    for key in missing:
        print("  %s: no row with this barcode" % key)

    animals = len({str(b).strip().lower() for _, b, _, _ in changes})
    print()
    print("  %d row(s) would gain a town, covering %d distinct accession(s)."
          % (len(changes), animals))
    print("  The map caption's 'N of 563 shown' moves by that second number.")
    if conflicts:
        print("  %d conflict(s) above need a decision before this is trusted."
              % len(conflicts))

    if not apply:
        print("\n  Dry run. Re-run with --apply to write.")
        return
    if not changes:
        print("\n  Nothing to write.")
        return

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    backup = "%s.%s.bak" % (BOOK, stamp)
    shutil.copy2(BOOK, backup)
    for r, _, town, _ in changes:
        ws.cell(row=r, column=tc).value = town
    wb.save(BOOK)
    print("\n  backup: %s" % backup)
    print("  wrote %d town(s) into %s" % (len(changes), BOOK))
    print("\n  Now re-run scripts/check_town_to_county.py to confirm they")
    print("  resolve to their county, then rebuild the collection map and")
    print("  update the caption number to match.")


if __name__ == "__main__":
    main()
