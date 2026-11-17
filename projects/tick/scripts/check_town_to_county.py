#!/usr/bin/env python3
"""check_town_to_county.py

I resolve each animal's collection town to a county and report how many fail.
I write nothing.

WHY FIRST
---------
The collection map stands or falls on town name matching, so I check that
before anything is drawn. If 40 towns do not resolve, the map is wrong in a way
that looks fine, and I would rather see the list now.

WHY TOWNS AND NOT COORDINATES
-----------------------------
The Coordinates column holds 3 values out of 573. The Town column holds 541.

WHY COUNTY AND NOT TOWN
-----------------------
These are hunter-harvest records. A town-level map publishes where individual
animals were taken, which is not ours to publish. County is coarse enough to be
safe and fine enough to show the sampling pattern across the three states.

HOW
---
In New England a town is a legal county subdivision, so the Census county
subdivision file already carries the county for each town and I do not have to
invent a lookup. I match on state and name, so the several towns that share a
name across states cannot collide.

I read the shapefiles straight out of their zips. Nothing is extracted and no
directory is created.

WHAT I READ
-----------
    results/metadata/crosswalk_4marker.tsv        the animal set
    metadata/tick_full_metadata.xlsx, 'Metadata'  barcode, town, state
    metadata/tl_2023_{50,33,23}_cousub.zip        town to county

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_town_to_county.py
"""
import io
import os
import re
import sys
import zipfile

import pandas as pd

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
BOOK = "metadata/tick_full_metadata.xlsx"
SHEET = "Metadata"
STATES = {"50": "VT", "33": "NH", "23": "ME"}
COUSUB = "metadata/tl_2023_%s_cousub.zip"
SHOW = 25


def shapefile_module():
    try:
        import shapefile
    except ImportError:
        sys.exit("pyshp is not importable. Install it with "
                 "'pip install --user pyshp'.")
    return shapefile


def tidy(value):
    """I normalise a place name so 'Milan ' and 'MILAN' are one thing."""
    text = re.sub(r"[^a-z0-9 ]", "", str(value).strip().lower())
    return re.sub(r"\s+", " ", text)


def read_towns(shapefile, code):
    """I return {normalised town name: (county FIPS, county-free name)}."""
    path = COUSUB % code
    if not os.path.isfile(path):
        sys.exit("missing %s. Download the Census county subdivision files "
                 "first; the command is in the session notes." % path)
    with zipfile.ZipFile(path) as archive:
        names = {os.path.splitext(n)[1].lower(): n for n in archive.namelist()}
        for needed in (".shp", ".dbf"):
            if needed not in names:
                sys.exit("%s has no %s inside it." % (path, needed))
        reader = shapefile.Reader(
            shp=io.BytesIO(archive.read(names[".shp"])),
            dbf=io.BytesIO(archive.read(names[".dbf"])))
        fields = [f[0] for f in reader.fields[1:]]
        for wanted in ("NAME", "COUNTYFP"):
            if wanted not in fields:
                sys.exit("%s has no %s field; it has %s"
                         % (path, wanted, ", ".join(fields)))
        found = {}
        for record in reader.records():
            row = dict(zip(fields, record))
            found[tidy(row["NAME"])] = (str(row["COUNTYFP"]),
                                        str(row["NAME"]))
        return found


def animal_lookup():
    """I return (table, module) for tick_lookup, the project's resolver."""
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s)." % problem)
    return tick_lookup.load(), tick_lookup


def animal_of(table, module, value):
    """I turn any barcode or sample id into its animal, or ''."""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    hit = module.resolve(table, "specimen", text)
    return "" if hit is None else str(hit)


def read_animals():
    """I return the analysis animals with their town and state.

    The crosswalk writes an animal as 'ALHT-Deer-102' and the Metadata sheet
    writes the same animal as 'ALHT WTD_23', so comparing the two barcode
    strings directly loses 62 animals. tick_lookup is the project's resolver
    for exactly this and it recovers 562 of the 563.

    I return the same frame twice so that callers written against the older
    two-value signature keep working.
    """
    if not os.path.isfile(CROSSWALK):
        sys.exit("missing %s. Run me from the project root." % CROSSWALK)
    if not os.path.isfile(BOOK):
        sys.exit("missing %s." % BOOK)
    table, module = animal_lookup()

    crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    for column in ("sample_id", "barcode", "state"):
        if column not in crosswalk.columns:
            sys.exit("%s has no '%s' column" % (CROSSWALK, column))
    crosswalk["animal"] = [animal_of(table, module, barcode)
                           or animal_of(table, module, sample)
                           for barcode, sample in zip(crosswalk["barcode"],
                                                      crosswalk["sample_id"])]
    stranded = int((crosswalk["animal"] == "").sum())
    if stranded:
        print("%d crosswalk rows resolve to no animal at all." % stranded)

    book = pd.read_excel(BOOK, sheet_name=SHEET, dtype=str)
    columns = {str(c).strip().lower(): c for c in book.columns}
    for wanted in ("barcode", "town", "state"):
        if wanted not in columns:
            sys.exit("the '%s' sheet has no '%s' column; it has %s"
                     % (SHEET, wanted, ", ".join(map(str, book.columns))))
    towns, states, clashes = {}, {}, {}
    for barcode, town, state in zip(book[columns["barcode"]],
                                    book[columns["town"]],
                                    book[columns["state"]]):
        animal = animal_of(table, module, barcode)
        if not animal:
            continue
        if pd.notna(town) and str(town).strip():
            here = str(town).strip()
            if animal in towns and tidy(towns[animal]) != tidy(here):
                clashes[animal] = (towns[animal], here)
            towns.setdefault(animal, here)
        if pd.notna(state) and str(state).strip():
            states.setdefault(animal, str(state).strip())
    if clashes:
        print("%d animals carry two different towns; I keep the first. "
              "For example %s"
              % (len(clashes), "; ".join("%s %s vs %s" % (a, x, y) for a,
                                         (x, y) in list(clashes.items())[:3])))

    crosswalk["town"] = [towns.get(a) for a in crosswalk["animal"]]
    crosswalk["book_state"] = [states.get(a) for a in crosswalk["animal"]]
    return crosswalk, crosswalk


def main():
    shapefile = shapefile_module()
    crosswalk, merged = read_animals()
    print("animals in the crosswalk:  %d" % len(crosswalk))
    print("matched to a metadata row: %d"
          % int(merged["book_state"].notna().sum()))
    print("with a town recorded:      %d" % int(merged["town"].notna().sum()))
    print()

    lookup = {}
    for code, abbreviation in STATES.items():
        towns = read_towns(shapefile, code)
        lookup[abbreviation] = towns
        print("%s: %d county subdivisions" % (abbreviation, len(towns)))
    print()

    resolved, failures, no_state, no_town = {}, {}, 0, 0
    for sample, state, town in zip(merged["sample_id"], merged["state"],
                                   merged["town"]):
        if pd.isna(town) or not str(town).strip():
            no_town += 1
            continue
        abbreviation = str(state).strip().upper()
        if abbreviation not in lookup:
            no_state += 1
            failures.setdefault(("no state", str(town).strip()), 0)
            failures[("no state", str(town).strip())] += 1
            continue
        hit = lookup[abbreviation].get(tidy(town))
        if hit is None:
            key = (abbreviation, str(town).strip())
            failures[key] = failures.get(key, 0) + 1
            continue
        county = "%s-%s" % (abbreviation, hit[0])
        resolved[county] = resolved.get(county, 0) + 1

    placed = sum(resolved.values())
    print("=" * 70)
    print("RESOLUTION")
    print("=" * 70)
    print("  animals placed in a county: %d of %d" % (placed, len(merged)))
    print("  no town recorded:           %d" % no_town)
    print("  state not one of VT/NH/ME:  %d" % no_state)
    print("  town not found:             %d" % sum(v for (s, _), v
                                                   in failures.items()
                                                   if s != "no state"))
    print()

    if failures:
        print("  towns I could not place, worst first:")
        for (state, town), count in sorted(failures.items(),
                                           key=lambda kv: -kv[1])[:SHOW]:
            print("    %-6s %-28s %4d animals" % (state, town[:28], count))
        print()
        print("  Each of these needs a real answer before the map is drawn.")
        print("  An unorganised township or a gore has a Census name that")
        print("  differs from the one on the collection sheet.")
        print()

    print("=" * 70)
    print("ANIMALS PER COUNTY")
    print("=" * 70)
    for county, count in sorted(resolved.items(), key=lambda kv: -kv[1]):
        print("  %-10s %4d" % (county, count))
    print()
    print("  %d counties hold animals." % len(resolved))
    print()


if __name__ == "__main__":
    main()
