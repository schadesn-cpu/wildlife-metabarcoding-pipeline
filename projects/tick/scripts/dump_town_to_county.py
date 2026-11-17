#!/usr/bin/env python3
"""dump_town_to_county.py  --  READ ONLY except for one output file.

WHY THIS EXISTS
---------------
check_town_to_county.py resolves towns to counties and prints ANIMALS PER
COUNTY, which is an aggregate. The spreadsheet C. Sullivan asked for needs the
county on each animal's row, so it needs the lookup itself rather than the
totals built from it.

Rather than invent a Vermont town-to-county table, this reuses exactly what
that script already does: the Census TIGER county subdivision files, matched on
state and normalised town name. tidy() and the shapefile reading are copied
from it unchanged, so this cannot disagree with the collection map. If the map
places an animal in a named county, so does this.

WHAT IT READS
    metadata/tick_full_metadata.xlsx, sheet 'Metadata'   barcode, town, state
    metadata/tl_2023_{50,33,23}_cousub.zip               town to county
    results/metadata/town_aliases.tsv                    optional, if present

WHAT IT WRITES
    results/metadata/town_to_county.tsv
        state, town_as_recorded, town_matched, county_name, county_fips

Towns it cannot place are written with empty county fields and listed at the
end, so a gap is visible in the file rather than silently absent. It will not
overwrite an existing output.

USAGE
    cd ~/meedlab/tick_project
    python3 scripts/dump_town_to_county.py
"""
import io
import os
import re
import sys
import zipfile

BOOK = "metadata/tick_full_metadata.xlsx"
SHEET = "Metadata"
ALIASES = "results/metadata/town_aliases.tsv"
STATES = {"50": "VT", "33": "NH", "23": "ME"}
COUSUB = "metadata/tl_2023_%s_cousub.zip"
OUT = "results/metadata/town_to_county.tsv"


def tidy(value):
    """Copied verbatim from check_town_to_county.py. Do not 'improve' it:
    the whole point is that this agrees with the map."""
    text = re.sub(r"[^a-z0-9 ]", "", str(value).strip().lower())
    return re.sub(r"\s+", " ", text)


def read_towns(shapefile, code):
    """{normalised town name: (county FIPS, county subdivision name)}"""
    path = COUSUB % code
    if not os.path.isfile(path):
        sys.exit("missing %s" % path)
    with zipfile.ZipFile(path) as archive:
        names = {os.path.splitext(n)[1].lower(): n for n in archive.namelist()}
        for needed in (".shp", ".dbf"):
            if needed not in names:
                sys.exit("%s has no %s inside it." % (path, needed))
        reader = shapefile.Reader(shp=io.BytesIO(archive.read(names[".shp"])),
                                  dbf=io.BytesIO(archive.read(names[".dbf"])))
        fields = [f[0] for f in reader.fields[1:]]
        for wanted in ("NAME", "COUNTYFP"):
            if wanted not in fields:
                sys.exit("%s has no %s field" % (path, wanted))
        found = {}
        for record in reader.records():
            row = dict(zip(fields, record))
            found[tidy(row["NAME"])] = (str(row["COUNTYFP"]), str(row["NAME"]))
        return found


def county_names(shapefile, code):
    """County FIPS -> county name, from the county file if it is here.
    A FIPS code is not a thing to hand a collaborator, so this turns it into
    a name. If the county shapefile is absent the FIPS still goes in the file
    and the name column is left empty rather than guessed."""
    path = "metadata/tl_2023_%s_county.zip" % code
    if not os.path.isfile(path):
        return {}
    try:
        with zipfile.ZipFile(path) as archive:
            names = {os.path.splitext(n)[1].lower(): n
                     for n in archive.namelist()}
            reader = shapefile.Reader(
                shp=io.BytesIO(archive.read(names[".shp"])),
                dbf=io.BytesIO(archive.read(names[".dbf"])))
            fields = [f[0] for f in reader.fields[1:]]
            out = {}
            for record in reader.records():
                row = dict(zip(fields, record))
                out[str(row.get("COUNTYFP"))] = str(row.get("NAME"))
            return out
    except (zipfile.BadZipFile, OSError, KeyError):
        return {}


def main():
    try:
        import shapefile
    except ImportError:
        sys.exit("pyshp is not importable. pip install --user pyshp")
    try:
        import pandas as pd
    except ImportError:
        sys.exit("pandas is not importable.")

    if not os.path.isfile(BOOK):
        sys.exit("missing %s. Run me from the project root." % BOOK)
    if os.path.exists(OUT):
        sys.exit("%s already exists. Move it aside first." % OUT)

    alias = {}
    if os.path.isfile(ALIASES):
        with open(ALIASES, encoding="utf-8", errors="replace") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                row = line.rstrip("\n").split("\t")
                if len(row) >= 2 and row[0].strip() and row[1].strip():
                    alias[tidy(row[0])] = row[1].strip()
        print("read %d town aliases from %s" % (len(alias), ALIASES))

    book = pd.read_excel(BOOK, sheet_name=SHEET, dtype=str)
    columns = {str(c).strip().lower(): c for c in book.columns}
    for wanted in ("town", "state"):
        if wanted not in columns:
            sys.exit("the '%s' sheet has no '%s' column" % (SHEET, wanted))

    pairs = set()
    for town, state in zip(book[columns["town"]], book[columns["state"]]):
        if pd.notna(town) and str(town).strip() and pd.notna(state):
            pairs.add((str(state).strip().upper(), str(town).strip()))
    print("distinct (state, town) pairs in the metadata: %d" % len(pairs))

    lookup, cnames = {}, {}
    for code, abbreviation in STATES.items():
        lookup[abbreviation] = read_towns(shapefile, code)
        cnames[abbreviation] = county_names(shapefile, code)
        print("%s: %d county subdivisions, %d county names"
              % (abbreviation, len(lookup[abbreviation]),
                 len(cnames[abbreviation])))

    rows, unplaced = [], []
    for state, town in sorted(pairs):
        key = tidy(town)
        if key in alias:
            key = tidy(alias[key])
        hit = lookup.get(state, {}).get(key)
        if hit is None:
            unplaced.append((state, town))
            rows.append((state, town, "", "", ""))
        else:
            fips, matched = hit
            rows.append((state, town, matched,
                         cnames.get(state, {}).get(fips, ""), fips))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("state\ttown_as_recorded\ttown_matched\tcounty_name\t"
                 "county_fips\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    print()
    print("wrote %s  (%d rows, %d placed, %d unplaced)"
          % (OUT, len(rows), len(rows) - len(unplaced), len(unplaced)))
    if unplaced:
        print()
        print("could not place, and these need a real answer before the")
        print("county column is trusted:")
        for state, town in unplaced:
            print("    %-4s %s" % (state, town))


if __name__ == "__main__":
    main()
