#!/usr/bin/env python3
"""patch_dermacentor_to_lookup.py

I put dermacentor_description.py on the shared library lookup and on the
species figure's own genus calls, so the project has one name resolver and one
tick-genus caller instead of three.

WHY
---
Checked 2026-09-19. The script reports 715 Ixodes plus 23 Dermacentor, which
is 738 keys from a project holding 563 specimens, against the species figure's
476 informative specimens with 406 Ixodes and 20 Dermacentor. Two causes, and
both make its numbers wrong:

  NO READ FLOOR. It takes the 'dominant' column of
  tick_genus_per_sample_corrected.tsv with no minimum. build_species_figure.py
  requires 50 reads AND 1% of that specimen's tick reads, because the blanks
  carried up to 28 Dermacentor reads. The three extra Dermacentor specimens
  sit in exactly the range the blanks reach.

  A BROKEN SPECIMEN KEY. normalise() strips hyphens, spaces and underscores
  alike, which is the bug that lost six single ticks in ground truth 1.14.
  to_barcode() strips a trailing '-NN' of exactly TWO digits, so
  EQ02441291-10 collapses to its specimen and EQ02441291-2 does not, which is
  the replicate-suffix bug of ground truth 1.20 in a new variant.

The outputs in results/analysis/dermacentor/ are dated 21 August and predate
the genus anchoring fix, the join fix and both crosswalk rebuilds.

WHAT I CHANGE
-------------
1. The tick-genus call is read from results/analysis/tick_species_by_host.tsv,
   which is what build_species_figure.py writes after applying the floor.
   There should be ONE genus caller in the project, and that is it.

2. Library names resolve to specimens through scripts/tick_lookup.py, the same
   way the other six consumers do. normalise() and to_barcode() are removed
   rather than left in place to be reused.

3. The collection host comes from the same species table, which already
   carries host_bin, instead of a second join into the crosswalk.

4. A guard: more specimens than the crosswalk holds means the key is wrong,
   and the script stops instead of reporting it.

5. The hard-coded 20 and ~265 in the prose are computed. A number typed into a
   sentence cannot go stale silently; one that is computed can only be wrong
   if the data is.

6. back_up_existing() is removed. It wrote a .bak beside the output, and a
   .bak from a previous run has been shared by accident once already. The
   report carries its own provenance, so the output is simply overwritten.

WHAT I DO NOT CHANGE
--------------------
The community analysis, the core-microbiota thresholds, the Rickettsia and
Francisella reciprocity, and every caveat in the report. Those are the content
and they are fine; only the join and the counting were wrong.

EXPECT ON THE NEXT RUN
----------------------
Dermacentor 20, not 23, and a specimen count at or below 476 rather than 738.
If Dermacentor comes back at anything other than 20 the species figure and
this script still disagree and that needs looking at before either is used.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/patch_dermacentor_to_lookup.py
    python3 scripts/dermacentor_description.py
"""
import os
import sys

TARGET = "scripts/dermacentor_description.py"

OLD_DOC = '''**There is no carriage contrast to test.** All 20 Dermacentor specimens are qPCR-negative on
every panel target, so there is no positive group. Nothing here is a hypothesis test and no
p-value is reported. n = 20 against ~265 Ixodes would not support one even if a contrast
existed.'''

NEW_DOC = '''**There is no carriage contrast to test.** The Dermacentor specimens are qPCR-negative on
every panel target, so there is no positive group. Nothing here is a hypothesis test and no
p-value is reported, and the group sizes would not support one even if a contrast existed.

Tick genus is read from results/analysis/tick_species_by_host.tsv, which build_species_figure.py
writes after applying the read floor, and library names resolve through scripts/tick_lookup.py.
This script used to do both itself and got both wrong; see patch_dermacentor_to_lookup.py.'''

OLD_CONST = '''TICK_GENUS = "results/18S/all/taxonomy_verified/tick_genus_per_sample_corrected.tsv"
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"'''

NEW_CONST = '''TICK_GENUS = "results/18S/all/taxonomy_verified/tick_genus_per_sample_corrected.tsv"
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
# The species figure's per-specimen genus calls, with the read floor already
# applied. One genus caller in the project, not two.
SPECIES_TABLE = "results/analysis/tick_species_by_host.tsv"
# The crosswalk holds 563 specimens. More keys than that means the key is not a
# specimen, which is how this script came to report 738 of them.
MAX_SPECIMENS = 563'''

OLD_IMPORT = '''import numpy as np
import pandas as pd'''

NEW_IMPORT = '''import numpy as np
import pandas as pd

import tick_lookup'''

OLD_HELPERS = '''def normalise(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def to_barcode(library):
    """Library name -> specimen barcode: strips the run suffix and the -NN replicate index."""
    return re.sub(r"-\\d{2}$", "", str(library).split(".")[0])


def back_up_existing(out_path):
    if not os.path.exists(out_path):
        return
    backup_path = out_path + ".bak"
    if os.path.exists(backup_path):
        sys.exit(f"{out_path} and {backup_path} both exist. "
                 f"Move or delete {backup_path} before rerunning.")
    os.replace(out_path, backup_path)'''

NEW_HELPERS = '''def resolve_specimen(lookup, library):
    """I turn a 16S library name into its specimen, through the shared lookup.

    The two helpers that used to live here did this themselves and did it
    wrong: one stripped hyphens so ALHT WTD_09 and alht-deer-09 became
    different keys, the other stripped a replicate index only when it was
    exactly two digits. Between them they made 738 keys out of 563 specimens.
    """
    text = str(library).strip()
    stem = text.split(".")[0]
    for marker, value in (("16S", text), ("16S", stem), ("specimen", stem)):
        hit = tick_lookup.resolve(lookup, marker, value)
        if hit is not None:
            return str(hit)
    return None'''

OLD_PATHS = '''    for path in (TABLE_QZA, TAXONOMY, TICK_GENUS):
        if not os.path.exists(path):
            sys.exit(f"no such file: {path}")
    os.makedirs(args.outdir, exist_ok=True)'''

NEW_PATHS = '''    for path in (TABLE_QZA, TAXONOMY, SPECIES_TABLE):
        if not os.path.exists(path):
            sys.exit(f"no such file: {path}")
    os.makedirs(args.outdir, exist_ok=True)
    lookup = tick_lookup.load()'''

OLD_GENUS = '''    # ---- tick genus per specimen -------------------------------------------
    calls = {}
    with open(TICK_GENUS, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\\t"):
            raw = (row.get("sample_id") or "").strip()
            value = (row.get("dominant") or "").strip()
            if raw and value:
                calls.setdefault(normalise(raw.split(".")[0]), set()).add(value)

    genus_by_barcode = {}
    for key, values in calls.items():
        ticks = {v for v in values if v in ("Ixodes", "Dermacentor", "Amblyomma")}
        if len(ticks) == 1:
            genus_by_barcode[key] = ticks.pop()
        elif len(ticks) > 1:
            genus_by_barcode[key] = "both/mixed"

    counts = pd.Series(genus_by_barcode).value_counts()
    say("specimens with a tick-genus call:")
    say(counts.to_string())

    if counts.get("Dermacentor", 0) < 3:
        sys.exit("fewer than 3 Dermacentor specimens; nothing to describe")

    # ---- host of each, for context -----------------------------------------
    host_by_barcode = {}
    if os.path.exists(CROSSWALK):
        with open(CROSSWALK, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\\t"):
                for key in (row.get("barcode"), row.get("sample_id")):
                    if key:
                        host_by_barcode[normalise(key)] = row.get("host_bin", "")'''

NEW_GENUS = '''    # ---- tick genus per specimen, from the species figure -------------------
    species = pd.read_csv(SPECIES_TABLE, sep="\\t", dtype=str)
    for column in ("specimen", "Ixodes_pos", "Dermacentor_pos"):
        if column not in species.columns:
            sys.exit(f"{SPECIES_TABLE} has no '{column}' column. Rebuild it:\\n"
                     f"    python3 scripts/build_species_figure.py")

    def is_true(value):
        """I read the species table's booleans, whatever casing they carry."""
        return str(value).strip().lower() in ("true", "1", "yes", "t")

    genus_by_barcode = {}
    host_by_barcode = {}
    for _, row in species.iterrows():
        specimen = str(row["specimen"]).strip()
        if not specimen or specimen.lower() in ("nan", "none"):
            continue
        positives = [name for name, column in (("Ixodes", "Ixodes_pos"),
                                               ("Dermacentor", "Dermacentor_pos"),
                                               ("Amblyomma", "Amblyomma_pos"))
                     if column in species.columns and is_true(row.get(column))]
        if len(positives) == 1:
            genus_by_barcode[specimen] = positives[0]
        elif len(positives) > 1:
            genus_by_barcode[specimen] = "both/mixed"
        host_by_barcode[specimen] = str(row.get("host_bin") or "").strip()

    counts = pd.Series(genus_by_barcode).value_counts()
    say(f"tick genus from {SPECIES_TABLE}, read floor already applied")
    say(f"specimens with a tick-genus call: {len(genus_by_barcode)}")
    say(counts.to_string())

    if len(genus_by_barcode) > MAX_SPECIMENS:
        sys.exit(f"{len(genus_by_barcode)} specimens carry a genus call, which is "
                 f"more than the {MAX_SPECIMENS} in the crosswalk. The key is not "
                 f"a specimen. Do not trust anything below this line.")

    if counts.get("Dermacentor", 0) < 3:
        sys.exit("fewer than 3 Dermacentor specimens; nothing to describe")

    # ---- host of each, for context, from the same table ---------------------'''

OLD_KEY = '''        key = normalise(to_barcode(column))
        if key not in genus_by_barcode:
            continue'''

NEW_KEY = '''        key = resolve_specimen(lookup, column)
        if key is None or key not in genus_by_barcode:
            continue'''

OLD_PROSE = '''    say("Descriptive ranking, not a test. n = 20 Dermacentor against ~265 Ixodes would not")
    say("support one, and the two groups differ in host, collection and sequencing depth.")

    ixodes = labels[labels == "Ixodes"].index
    dermacentor = labels[labels == "Dermacentor"].index'''

NEW_PROSE = '''    ixodes = labels[labels == "Ixodes"].index
    dermacentor = labels[labels == "Dermacentor"].index
    say(f"Descriptive ranking, not a test. n = {len(dermacentor)} Dermacentor against "
        f"{len(ixodes)} Ixodes would not")
    say("support one, and the two groups differ in host, collection and sequencing depth.")'''

OLD_BACKUP = '''        out_path = f"{args.outdir}/core_microbiota_by_genus.tsv"
        back_up_existing(out_path)
        pd.DataFrame(core_rows).to_csv(out_path, sep="\\t", index=False)'''

NEW_BACKUP = '''        out_path = f"{args.outdir}/core_microbiota_by_genus.tsv"
        # No .bak. The old backup step left one beside the output, and a .bak
        # from a previous run has been shared by accident once already. The
        # report below records what this run did.
        pd.DataFrame(core_rows).to_csv(out_path, sep="\\t", index=False)'''

OLD_CAVEAT = '''    say("  - NO hypothesis test is run here and no p-value is reported. All 20 Dermacentor")
    say("    specimens are qPCR-negative on every panel target, so no carriage contrast")
    say("    exists within Dermacentor to test.")'''

NEW_CAVEAT = '''    say(f"  - NO hypothesis test is run here and no p-value is reported. The "
        f"{len(dermacentor)} Dermacentor")
    say("    specimens were qPCR-negative on every panel target when last checked, so no")
    say("    carriage contrast exists within Dermacentor to test. [VERIFY that against the")
    say("    current qPCR join before quoting it; it was established on an older cohort.]")'''

REPLACEMENTS = [
    (OLD_DOC, NEW_DOC, "the module docstring"),
    (OLD_CONST, NEW_CONST, "the constants block"),
    (OLD_IMPORT, NEW_IMPORT, "the imports"),
    (OLD_HELPERS, NEW_HELPERS, "normalise / to_barcode / back_up_existing"),
    (OLD_PATHS, NEW_PATHS, "the input path check"),
    (OLD_GENUS, NEW_GENUS, "the tick-genus and host blocks"),
    (OLD_KEY, NEW_KEY, "the library-to-specimen key"),
    (OLD_PROSE, NEW_PROSE, "the hard-coded 20 and 265"),
    (OLD_BACKUP, NEW_BACKUP, "the .bak step"),
    (OLD_CAVEAT, NEW_CAVEAT, "the closing caveat"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from the project root.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "SPECIES_TABLE" in source:
        sys.exit("already patched: " + TARGET + " reads SPECIES_TABLE. "
                 "Nothing changed.")
    if not os.path.isfile("scripts/tick_lookup.py"):
        sys.exit("scripts/tick_lookup.py not found. Build the lookup first:\n"
                 "    python3 scripts/build_4marker_crosswalk.py")

    for old, new, label in REPLACEMENTS:
        count = source.count(old)
        if count != 1:
            sys.exit("in " + TARGET + ": expected exactly one occurrence of "
                     + label + ", found " + str(count)
                     + ". Patch by hand rather than trusting this.")
        source = source.replace(old, new, 1)

    compile(source, TARGET, "exec")
    with open(TARGET, "w", encoding="utf-8") as handle:
        handle.write(source)
    print("patched " + TARGET)
    print()
    print("Expect Dermacentor 20, not 23, and a specimen count at or below 476")
    print("rather than 738. Anything else means the species figure and this")
    print("script still disagree.")
    print()
    print("run: python3 " + TARGET)


if __name__ == "__main__":
    main()
