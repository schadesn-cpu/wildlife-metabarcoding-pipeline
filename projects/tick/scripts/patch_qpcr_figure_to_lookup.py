#!/usr/bin/env python3
"""patch_qpcr_figure_to_lookup.py

I put build_qpcr_figure.py on the shared library lookup and on the crosswalk's
tick counts, replacing a local normaliser and a run-sheet tick-count map.

TWO BUGS BEING REMOVED, AND THE SECOND IS THE SERIOUS ONE
---------------------------------------------------------
1. _norm() stripped every non-alphanumeric character, so 'ALHT WTD_09' became
   'alhtwtd09' and the qPCR sheet's 'alht-deer-09' became 'alhtdeer09'. Two
   keys, one specimen. This is the same defect that cost build_concordance.py
   six ticks.

2. The tick-count map was built from the RAW run sheet's 'Num Ticks' column,
   with the crosswalk consulted only where the run sheet was blank. That is
   backwards. 'Num Ticks' is the column the ALHT/LHT barcode clobber
   corrupted, and build_4marker_crosswalk.py exists to correct it. Checked in
   the workbook:

       ALHT WTD_09   run sheet: 18 ticks   registry and crosswalk: 1 tick

   The run-sheet row carries Barcode 'USGSTICKS-LHT-Deer-091', which is the
   clobbered row Lawrence Gordon identified on 2026-08-07. So a genuine single
   tick is read as an 18-tick pool and drops out of the per-tick cohort, and
   the crosswalk fallback never fires because the run sheet is not blank, it
   is wrong.

   The crosswalk is the authority on tick counts. After this patch the run
   sheet is not consulted for them at all.

COHORT DEFINITION
-----------------
I select on the crosswalk's is_single_tick rather than on num_ticks == 1.
That flag is the project's single stated rule (build_4marker_crosswalk.py
docstring note 4): one tick by the collection record, nothing contradicting
it, wild-caught in the study region. It additionally excludes the 8 Batch-12
fungal-treatment ticks and the 4 specimens whose extraction record contradicts
the collection record. The sequencing cohort already obeys that rule, so
applying it here makes the two cohorts comparable instead of nearly
comparable. Both counts are printed so the difference is visible.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_qpcr_figure_to_lookup.py
    python3 scripts/build_qpcr_figure.py
"""
import os
import sys

TARGET = "scripts/build_qpcr_figure.py"
NEEDS = "scripts/tick_lookup.py"

OLD_RESOLVE = '''import re as _re
def _norm(value):
    return _re.sub(r"[^a-z0-9]", "", str(value).lower())

# The qPCR sheet keys on host_id ('alht-deer-02') while the crosswalk keys on
# the accession ('ALHT WTD_02'). Matching the crosswalk directly leaves 64 of
# the 423 unresolved. The compiled run sheet carries BOTH forms on one row, so
# I build the tick-count lookup from every alias there instead.
_run = pd.read_excel(WORKBOOK, sheet_name="All Run Data - w Metadata")
_run.columns = [str(c) for c in _run.columns]
_sample_cols = [c for c in _run.columns if c.startswith("Sample")]
_tick_count = {}
for _, _row in _run.iterrows():
    _n = pd.to_numeric(pd.Series([_row["Num Ticks"]]), errors="coerce")[0]
    for _alias in [_row["Barcode"]] + [_row[c] for c in _sample_cols]:
        if pd.notna(_alias):
            _tick_count.setdefault(_norm(_alias), _n)

# crosswalk fallback: build_4marker_crosswalk.py fills a tick count for a few
# specimens the run sheet leaves blank (24-3889 among them, from the registry
# fill). Without this they drop out of the per-tick cohort even though the
# crosswalk knows they are single ticks.
_CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
if os.path.exists(_CROSSWALK):
    _cw = pd.read_csv(_CROSSWALK, sep="\\t", dtype=str)
    _cw_ticks = pd.to_numeric(_cw["num_ticks"], errors="coerce")
    for _sid, _n in zip(_cw["sample_id"], _cw_ticks):
        _key = _norm(_sid)
        if pd.notna(_n) and (_key not in _tick_count or pd.isna(_tick_count[_key])):
            _tick_count[_key] = _n

_key_column = "Barcode" if "Barcode" in qpcr_all.columns else "host_id"
qpcr_all["_key"] = [
    _norm(row["Barcode"]) if pd.notna(row.get("Barcode")) else _norm(row.get("host_id"))
    for _, row in qpcr_all.iterrows()]
qpcr_all["num_ticks"] = qpcr_all["_key"].map(_tick_count)

'''

NEW_RESOLVE = '''# Name resolution is NOT done here, and neither is tick counting. Both come
# from build_4marker_crosswalk.py, via results/metadata/library_lookup.tsv and
# results/metadata/crosswalk_4marker.tsv.
#
# What used to sit in this slot did two harmful things. It normalised names by
# stripping every non-alphanumeric character, so 'ALHT WTD_09' and the qPCR
# sheet's 'alht-deer-09' became two different keys. And it took tick counts
# from the RAW run sheet, falling back to the crosswalk only where the run
# sheet was blank. The run sheet's 'Num Ticks' is the column the ALHT/LHT
# barcode clobber corrupted: ALHT WTD_09 reads 18 ticks there against 1 in the
# registry and the crosswalk, because its run row carries the clobbered
# 'USGSTICKS-LHT-Deer-091' barcode. A genuine single tick was therefore read
# as an 18-tick pool, and the fallback never fired because the value was not
# blank, only wrong.
#
# The crosswalk is the authority. Do not reintroduce a run-sheet lookup.
import tick_lookup

LOOKUP = tick_lookup.load()
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
if not os.path.exists(CROSSWALK):
    sys.exit("ERROR: " + CROSSWALK + " not found. Build it first:\\n"
             "    python3 scripts/build_4marker_crosswalk.py")

crosswalk = pd.read_csv(CROSSWALK, sep="\\t", dtype=str)
for _needed in ("sample_id", "num_ticks", "is_single_tick"):
    if _needed not in crosswalk.columns:
        sys.exit("ERROR: column '" + _needed + "' missing from " + CROSSWALK)

tick_count_of = {}
single_tick_flag = {}
for _specimen, _n, _flag in zip(crosswalk["sample_id"], crosswalk["num_ticks"],
                                crosswalk["is_single_tick"]):
    tick_count_of[_specimen] = pd.to_numeric(pd.Series([_n]), errors="coerce")[0]
    single_tick_flag[_specimen] = str(_flag).strip().lower() == "true"


def specimen_of(row):
    """I resolve one qPCR row to a crosswalk specimen, or None if it is absent.

    None is a real answer. Nine rows in this sheet have no crosswalk specimen,
    two of them pathogen-positive, and they should stay visible rather than
    being forced into a denominator.
    """
    for column in ("Barcode", "host_id"):
        hit = tick_lookup.resolve(LOOKUP, "specimen", row.get(column))
        if hit is not None:
            return hit
    return None


qpcr_all["specimen"] = [specimen_of(row) for _, row in qpcr_all.iterrows()]
qpcr_all["num_ticks"] = qpcr_all["specimen"].map(
    lambda s: tick_count_of.get(s) if s is not None else None)
qpcr_all["is_single_tick"] = qpcr_all["specimen"].map(
    lambda s: bool(single_tick_flag.get(s, False)) if s is not None else False)

_no_specimen = qpcr_all[qpcr_all["specimen"].isna()]
say("qPCR rows with no crosswalk specimen: " + str(len(_no_specimen)))
for _, _row in _no_specimen.iterrows():
    say("    " + str(_row.get("host_id")) + "  (Barcode "
        + str(_row.get("Barcode")) + ")")

'''

OLD_SELECT = '''pooled = qpcr_all[qpcr_all["subsamples"] == 1].copy()
single_ticks = pooled[pooled["num_ticks"] == 1].copy()
'''

NEW_SELECT = '''pooled = qpcr_all[qpcr_all["subsamples"] == 1].copy()
# is_single_tick, not num_ticks == 1. The flag is the project's one stated rule
# (build_4marker_crosswalk.py docstring note 4) and additionally excludes the
# Batch-12 fungal-treatment ticks and the specimens whose extraction record
# contradicts the collection record. The sequencing cohort already obeys it.
_by_count = pooled[pooled["num_ticks"] == 1]
single_ticks = pooled[pooled["is_single_tick"]].copy()
say("one-subsample rows from a one-tick tube:          " + str(len(_by_count)))
say("of those, eligible by the is_single_tick rule:    " + str(len(single_ticks)))
if len(_by_count) != len(single_ticks):
    _dropped = sorted(set(_by_count["specimen"].dropna())
                      - set(single_ticks["specimen"].dropna()))
    say("  excluded by the rule: " + ", ".join(_dropped))
'''

OLD_WARN = '''⚠️ THIS IS A DIFFERENT POPULATION FROM THE SEQUENCING FIGURE. 99 qPCR single
ticks vs 93 sequenced single ticks. Panel B here and
panel B of the metabarcoding figure are NOT the same ticks, and the qPCR targets
are species while the metabarcoding co-infection counts genera. I expect qPCR
co-infection to run higher, because qPCR finds A. phagocytophilum at 50.1% and
B. microti at 30.0% where sequencing finds 18.7% and 4.4%. That difference is a
result worth writing about, not an inconsistency to smooth over.'''

NEW_WARN = '''⚠️ THIS IS A DIFFERENT POPULATION FROM THE SEQUENCING FIGURE. The qPCR
per-tick cohort and the 93 sequenced single ticks overlap in 89, which is the
concordance cohort, but they are not the same set. Panel B here and panel B of
the metabarcoding figure are NOT the same ticks, and the qPCR targets are
species while the metabarcoding co-infection counts genera. I expect qPCR
co-infection to run higher, because qPCR finds A. phagocytophilum and B. microti
at pool level at 50.1% and 30.0% where sequencing finds 22.6% and 4.3% per tick.
That difference is a result worth writing about, not an inconsistency to smooth
over.'''

OLD_CAP = '''say("  - a DIFFERENT population from the metabarcoding figure (93 sequenced single")
say("    ticks, overlapping in 83). The two panel-B co-infection figures are not")
say("    the same ticks and are not directly comparable")'''

NEW_CAP = '''say("  - a DIFFERENT population from the metabarcoding figure (93 sequenced single")
say("    ticks, overlapping in 89). The two panel-B co-infection figures are not")
say("    the same ticks and are not directly comparable")'''

OLD_DEF = '''say("  - this panel counts SPECIES across a 7-target qPCR panel; the")
say("    metabarcoding panel counts GENERA with both Babesia collapsed.")
say("    The difference between 14.0% here and 20.4% there is definitional.")'''

NEW_DEF = '''say("  - this panel counts SPECIES across a 7-target qPCR panel; the")
say("    metabarcoding panel counts GENERA with both Babesia collapsed and")
say("    Bartonella excluded as a ked-associated co-detection (15.1% there).")
say("    Part of the difference between the two is definitional, not biological.")'''

OLD_NAMES = '''    "b-burgdorferi":     r"$\\mathit{Borrelia\\ burgdorferi}$",'''
NEW_NAMES = '''    "b-burgdorferi":     r"$\\mathit{Borreliella\\ burgdorferi}$",'''

OLD_MAY = '''    "b-mayonii":         r"$\\mathit{Borrelia\\ mayonii}$",'''
NEW_MAY = '''    "b-mayonii":         r"$\\mathit{Borreliella\\ mayonii}$",'''

REPLACEMENTS = [
    (OLD_WARN, NEW_WARN, "the different-population warning"),
    (OLD_RESOLVE, NEW_RESOLVE, "the name-resolution and tick-count block"),
    (OLD_SELECT, NEW_SELECT, "the cohort selection"),
    (OLD_CAP, NEW_CAP, "the overlap line in the caption block"),
    (OLD_DEF, NEW_DEF, "the definitional-difference caption line"),
    (OLD_NAMES, NEW_NAMES, "the B. burgdorferi label"),
    (OLD_MAY, NEW_MAY, "the B. mayonii label"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    if not os.path.isfile(NEEDS):
        sys.exit("not found: " + NEEDS + "\nWrite tick_lookup.py first.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "import tick_lookup" in source:
        sys.exit("already patched: " + TARGET + " imports tick_lookup. "
                 "Nothing changed.")

    for old, new, label in REPLACEMENTS:
        count = source.count(old)
        if count != 1:
            sys.exit("expected exactly one occurrence of " + label
                     + ", found " + str(count)
                     + ". Patch by hand rather than trusting this.")
        source = source.replace(old, new, 1)

    compile(source, TARGET, "exec")
    with open(TARGET, "w", encoding="utf-8") as handle:
        handle.write(source)
    print("patched " + TARGET)
    print("now run: python3 " + TARGET)
    print("ALHT WTD_09 should now be IN the per-tick cohort; it was read as an")
    print("18-tick pool before. Expect the cohort to be larger than 100.")


if __name__ == "__main__":
    main()
