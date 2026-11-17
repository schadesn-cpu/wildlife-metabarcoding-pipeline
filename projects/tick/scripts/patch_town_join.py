#!/usr/bin/env python3
"""patch_town_join.py

I change check_town_to_county.py so it joins the crosswalk to the Metadata
sheet through tick_lookup instead of comparing barcode strings. I change
nothing else and I refuse rather than guess.

WHY
---
The crosswalk writes an animal as 'ALHT-Deer-102' and the Metadata sheet
writes the same animal as 'ALHT WTD_23'. Comparing those strings matched 501
of 563 animals and left 62 looking as though they had no collection record at
all. They have one. tick_lookup is the project's resolver for exactly this
disagreement and it reaches 562 of 563.

Every script that reads collection towns goes through read_animals, so fixing
it here fixes the alias proposal and the map as well.

WHAT ELSE CHANGES
-----------------
One animal can now carry two metadata rows. I keep the first town and print a
warning naming the animals whose two rows disagree, rather than picking one
silently.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/patch_town_join.py
    python3 scripts/patch_town_join.py --write
"""
import os
import sys

TARGET = "scripts/check_town_to_county.py"
BEFORE = "e1a9573b934670af192ce77bd39e795e"
EXPECT = "7fa14cf5b1223ed82a4295f43c2496d1"

EDITS = [
("""def read_animals():
    \"\"\"I return the analysis animals with their town and state.\"\"\"
    if not os.path.isfile(CROSSWALK):""",
 """def animal_lookup():
    \"\"\"I return (table, module) for tick_lookup, the project's resolver.\"\"\"
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s)." % problem)
    return tick_lookup.load(), tick_lookup


def animal_of(table, module, value):
    \"\"\"I turn any barcode or sample id into its animal, or ''.\"\"\"
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    hit = module.resolve(table, "specimen", text)
    return "" if hit is None else str(hit)


def read_animals():
    \"\"\"I return the analysis animals with their town and state.

    The crosswalk writes an animal as 'ALHT-Deer-102' and the Metadata sheet
    writes the same animal as 'ALHT WTD_23', so comparing the two barcode
    strings directly loses 62 animals. tick_lookup is the project's resolver
    for exactly this and it recovers 562 of the 563.

    I return the same frame twice so that callers written against the older
    two-value signature keep working.
    \"\"\"
    if not os.path.isfile(CROSSWALK):"""),

("""        sys.exit("missing %s." % BOOK)
    crosswalk = pd.read_csv(""",
 """        sys.exit("missing %s." % BOOK)
    table, module = animal_lookup()

    crosswalk = pd.read_csv("""),

("""            sys.exit("%s has no '%s' column" % (CROSSWALK, column))

    book = pd.read_excel(BOOK, sheet_name=SHEET, dtype=str)""",
 """            sys.exit("%s has no '%s' column" % (CROSSWALK, column))
    crosswalk["animal"] = [animal_of(table, module, barcode)
                           or animal_of(table, module, sample)
                           for barcode, sample in zip(crosswalk["barcode"],
                                                      crosswalk["sample_id"])]
    stranded = int((crosswalk["animal"] == "").sum())
    if stranded:
        print("%d crosswalk rows resolve to no animal at all." % stranded)

    book = pd.read_excel(BOOK, sheet_name=SHEET, dtype=str)"""),

("""    book = book.rename(columns={columns["barcode"]: "barcode",
                                columns["town"]: "town",
                                columns["state"]: "book_state"})
    book["barcode"] = book["barcode"].astype(str).str.strip()
    crosswalk["barcode"] = crosswalk["barcode"].astype(str).str.strip()
    merged = crosswalk.merge(book[["barcode", "town", "book_state"]],
                             on="barcode", how="left")
    return crosswalk, merged""",
 """    towns, states, clashes = {}, {}, {}
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
    return crosswalk, crosswalk"""),
]


def main():
    write = "--write" in sys.argv[1:]
    if not os.path.isfile(TARGET):
        sys.exit("missing %s. Run me from the project root." % TARGET)
    text = open(TARGET, encoding="utf-8").read()

    for number, (old, new) in enumerate(EDITS, start=1):
        found = text.count(old)
        if found != 1:
            sys.exit("edit %d matches %d places in %s, not 1. Either it has "
                     "been patched already or it is not the version I expect "
                     "(%s). Nothing written." % (number, found, TARGET, BEFORE))
        text = text.replace(old, new)
        print("  edit %d: applied" % number)

    print()
    print("all %d edits apply cleanly to %s" % (len(EDITS), TARGET))
    if not write:
        print("Nothing written. Rerun with --write to save the change.")
        return
    open(TARGET, "w", encoding="utf-8").write(text)
    print("wrote %s" % TARGET)
    print("Check it: sed '/^[[:space:]]*$/d' %s | md5sum" % TARGET)
    print("should give %s" % EXPECT)


if __name__ == "__main__":
    main()
