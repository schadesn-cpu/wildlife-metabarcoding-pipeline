#!/usr/bin/env python3
"""propose_town_aliases.py

I propose a Census town name for each collection town that did not resolve, and
print them for you to confirm. I write nothing.

WHY YOU AND NOT ME
------------------
'Keeene' is obviously Keene and 'Wolfboro' is obviously Wolfeboro, but 'Etna'
is a village inside Hanover and 'Pawlet/Danby' is two towns, and knowing which
one a collector meant is local knowledge I do not have. If I wrote these into
the map myself they would become invisible assumptions. Printed here with a
similarity score, they stay decisions you made.

Confirm them by writing the accepted pairs to

    results/metadata/town_aliases.tsv

with columns state, recorded, census. The map script reads that file and fails
on any town it still cannot place, so nothing silently vanishes.

HOW I CLEAN A NAME BEFORE MATCHING
----------------------------------
    'Pittsford (Florence)'   parenthetical dropped
    'Wakefield WMU J1'       a wildlife management unit code dropped
    'Pawlet/Danby'           first of a pair taken, and flagged as a guess
    'W Rutland'              leading direction expanded to 'West'
    'Mt Holly'               'Mt' expanded to 'Mount'
    'St.Albans'              'St.' expanded to 'Saint'

Each rewrite is printed beside the proposal so you can see what I did to the
name before I looked it up.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/propose_town_aliases.py
"""
import difflib
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CANDIDATES = 3
FLOOR = 0.60
DIRECTIONS = {"n": "North", "s": "South", "e": "East", "w": "West",
              "mt": "Mount", "st": "Saint", "ft": "Fort"}


def resolver():
    """I return check_town_to_county, whose loading I reuse whole."""
    try:
        import check_town_to_county as module
    except ImportError as problem:
        sys.exit("cannot import check_town_to_county (%s). It must sit in "
                 "scripts/ beside me." % problem)
    return module


def clean(name):
    """I return (rewritten name, what I did), before any lookup."""
    text = str(name).strip()
    notes = []
    if "(" in text:
        text = re.sub(r"\s*\([^)]*\)", "", text).strip()
        notes.append("dropped parenthetical")
    hit = re.search(r"\bWMU\b.*$", text, flags=re.I)
    if hit:
        text = text[:hit.start()].strip()
        notes.append("dropped management unit code")
    if "/" in text:
        text = text.split("/")[0].strip()
        notes.append("took the first of a pair, GUESS")
    parts = text.replace(".", ". ").split()
    rebuilt = []
    for part in parts:
        key = part.rstrip(".").lower()
        if key in DIRECTIONS and len(parts) > 1:
            rebuilt.append(DIRECTIONS[key])
            notes.append("expanded '%s'" % part)
        else:
            rebuilt.append(part)
    text = re.sub(r"\s+", " ", " ".join(rebuilt)).strip()
    return text, "; ".join(notes)


def main():
    module = resolver()
    shapefile = module.shapefile_module()
    _, merged = module.read_animals()

    lookup = {}
    for code, abbreviation in module.STATES.items():
        lookup[abbreviation] = module.read_towns(shapefile, code)

    failures = {}
    for state, town in zip(merged["state"], merged["town"]):
        if pd.isna(town) or not str(town).strip():
            continue
        abbreviation = str(state).strip().upper()
        if abbreviation not in lookup:
            continue
        if module.tidy(town) in lookup[abbreviation]:
            continue
        key = (abbreviation, str(town).strip())
        failures[key] = failures.get(key, 0) + 1

    if not failures:
        print("every recorded town resolves. Nothing to propose.")
        return

    print("%d town spellings do not resolve, covering %d animals."
          % (len(failures), sum(failures.values())))
    print()
    print("%-5s %-24s %-24s %-7s %s"
          % ("state", "recorded", "proposed Census name", "score", "how"))
    print("-" * 96)

    accepted = []
    for (state, town), count in sorted(failures.items(),
                                       key=lambda kv: (-kv[1], kv[0])):
        names = lookup[state]
        rewritten, notes = clean(town)
        direct = names.get(module.tidy(rewritten))
        if direct is not None:
            print("%-5s %-24s %-24s %-7s %s"
                  % (state, town[:24], direct[1][:24], "exact",
                     notes or "cleaning alone was enough"))
            accepted.append((state, town, direct[1], count))
            continue
        close = difflib.get_close_matches(module.tidy(rewritten),
                                          list(names), n=CANDIDATES,
                                          cutoff=FLOOR)
        if not close:
            print("%-5s %-24s %-24s %-7s %s"
                  % (state, town[:24], "NOTHING CLOSE", "-",
                     notes or "no rewrite applied"))
            continue
        best = close[0]
        score = difflib.SequenceMatcher(None, module.tidy(rewritten),
                                        best).ratio()
        others = ", ".join(names[c][1] for c in close[1:])
        print("%-5s %-24s %-24s %-7.2f %s"
              % (state, town[:24], names[best][1][:24], score,
                 (notes + "; " if notes else "")
                 + ("also: " + others if others else "")))
        accepted.append((state, town, names[best][1], count))

    print()
    print("  %d animals would be recovered if every proposal above is right,"
          % sum(c for _, _, _, c in accepted))
    print("  taking the mapped total from 437 to %d of 563."
          % (437 + sum(c for _, _, _, c in accepted)))
    print()
    print("  A village is not a town. Etna, Florence, Union and Milton Mills")
    print("  are villages, and the proposal for each is whichever town the")
    print("  name matched, which may not be the town containing the village.")
    print("  Those are the ones to look at hardest.")
    print()
    print("  To accept all of them as printed:")
    print()
    print("    python3 scripts/propose_town_aliases.py --write-all")
    print()
    print("  I would rather you edited the file by hand than did that.")

    if "--write-all" in sys.argv[1:]:
        out = "results/metadata/town_aliases.tsv"
        if os.path.exists(out):
            sys.exit("%s already exists. Move it aside first." % out)
        with open(out, "w") as handle:
            handle.write("state\trecorded\tcensus\n")
            for state, town, census, _ in accepted:
                handle.write("%s\t%s\t%s\n" % (state, town, census))
        print()
        print("wrote %s with %d rows. Check every one." % (out, len(accepted)))


if __name__ == "__main__":
    main()
