#!/usr/bin/env python3
"""
build_4marker_crosswalk.py  —  run ON RON, from project root:
    cd ~/meedlab/tick_project
    python3 scripts/build_4marker_crosswalk.py

Ties all four marker library names to each specimen, so no future join has to
regex a name again. Anchors on the ALREADY-VERIFIED 16S crosswalk and adds the
other three markers from the compiled workbook.

═══════════════════════════════════════════════════════════════════════════════
REV 2026-08-07 — three changes. Read these before trusting any downstream number.

1. TWO OF KELLEN'S FIXES WERE WRONG. Corrected per Lawrence Gordon 2026-08-07:
     - Num 247 is ALHT WTD_09, NOT a merge into USGSTICKS-LHT-Deer-091
     - Nums 617/618/619 are USGSTICKS-2213, NOT EQ02442213
   Both had the wrong specimen's host, town, date and TICK COUNT attached.

2. NEW: specimen metadata now comes from the 'Metadata' sheet (the SPECIMEN
   REGISTRY, 573 rows, one per barcode) when a barcode has no run-data row.
   'All Run Data - w Metadata' is a LIBRARY registry with specimen fields copied
   in; where a barcode appears in both they agree exactly (verified 2026-08-07,
   0 tick-count mismatches), but two corrected specimens appear only in the
   registry.
   NOTE: the registry's barcode column is named 'Barcode ' WITH A TRAILING SPACE.

3. NEW: match-rate assertions. Every join reports its coverage and the build FAILS
   LOUDLY on an unexpected drop. Had this existed in April, the 2913 typo, the
   LHT-091 merge and the Rosetta Stone offset would all have surfaced then.

4. SINGLE-TICK INCLUSION RULE (one rule, stated once):
   A specimen is a single tick if, and only if, (a) its collection record gives a
   tick count of exactly 1, (b) no other record contradicts that, and (c) it is a
   wild-caught animal from the study region. Everything failing any clause is kept
   in the crosswalk and in pooled/community analyses but is NOT eligible to be
   called a single tick. Excluded on that basis: 9 Batch-12 fungal-treatment ticks
   (experimental, not wild-caught) and 4 specimens whose extraction record
   contradicts the collection record.

5. TICK COUNTS: 35 specimens have a blank 'Num Ticks' (the USGSTICKS-VT / 23xx
   series and the five VT bears). The batch extraction sheets DO carry a per-tube
   'Tick count', but for 105 of 156 specimens where both exist the two columns
   DISAGREE (e.g. USGSTICKS-6503: Num Ticks 1, Batch-10 Tick count 2, one tube,
   no duplication). We could not establish what the two columns each mean, and
   read depth cannot arbitrate (depth reflects library pooling, not tick number).
   DECISION 2026-08-07: do not guess. Specimens with a blank or contradicted tick
   count stay out of the single-tick set. Reported as a WARNING, not a failure, so
   everything not downstream of tick counts still builds.

Reads:
    results/metadata/crosswalk_VERIFIED.tsv           (trusted 16S anchor, July 13)
    metadata/tick_full_metadata.xlsx  sheet 'All Run Data - w Metadata'  (libraries)
    metadata/tick_full_metadata.xlsx  sheet 'Metadata'                   (specimens)
Writes:
    results/metadata/crosswalk_4marker.tsv            (the full crosswalk — CANONICAL)
    results/metadata/crosswalk_label_mismatches.tsv   (should be empty)
    results/metadata/crosswalk_build_report.txt       (the assertions, for the record)
"""
import sys
import pandas as pd, re

XW    = "results/metadata/crosswalk_VERIFIED.tsv"
WB    = "metadata/tick_full_metadata.xlsx"
SHEET = "All Run Data - w Metadata"
REG   = "Metadata"          # the specimen registry

_report = []
def say(msg=""):
    print(msg)
    _report.append(str(msg))

def check(condition, msg):
    """Assertion that explains itself. Fails the build rather than warning."""
    if condition:
        say(f"  OK   {msg}")
    else:
        say(f"  FAIL {msg}")
        say("\nBUILD FAILED. Nothing written. Fix the above before re-running.")
        open("results/metadata/crosswalk_build_report.txt", "w").write("\n".join(_report))
        sys.exit(1)

# ═══════════════════════════════════════════════════════════════════════════
# LABEL FIXES
# Kellen Sawyer 2026-07-31 (David: "Kellen's answers are good"), as corrected by
# Lawrence Gordon 2026-08-07. Applied explicitly, never fuzzy-matched.
# ═══════════════════════════════════════════════════════════════════════════

# (a) Same specimen under two naming schemes -> collapse to one specimen.
#     VERIFIED 2026-08-07: for all five pairs BOTH names appear as separate
#     run-data rows and the 'ALHT WTD_*' form is the registry barcode, so
#     collapsing to it is correct.
merge_map = {
    'ALHT-Deer-102':'ALHT WTD_102', 'ALHT-Deer-64':'ALHT WTD_64',
    'ALHT-Deer-92':'ALHT WTD_92',   'ALHT-Deer-93':'ALHT WTD_93',
    'ALHT-Deer-98':'ALHT WTD_98',
    # Kellen: "4065 = TB-90024065: agree, it should be the moose."
    # (The workbook has this specimen twice, Num 95 hyphen + Num 97 underscore.)
    'TB_90024065':'TB-90024065',

    # ─── REMOVED 2026-08-07 ────────────────────────────────────────────────
    # 'USGSTICKS-LHT-Deer-091':'ALHT-Deer-09',
    # Kellen read these as the same specimen. They are NOT. Lawrence:
    #   ALHT-Deer-09           = 1 tick,  collected 2023-10-21
    #   USGSTICKS-LHT-Deer-091 = 18 ticks, collected 2024-11-16, different town
    # Num 247's library is named ALHT-Deer-09-16s and was sequenced 12/19/24,
    # BEFORE the LHT batch was accessioned 2/18/2025. Num 815 is the real LHT-091.
    # The merge made ALHT-Deer-09 inherit num_ticks=18, i.e. it was recorded as a
    # POOL when it is a SINGLE TICK -> it was wrongly excluded from the n=103 set.
    # Now handled as a re-attribution in LABEL_FIXES below.
}

# (b) Mislabeled library to IGNORE completely (drop the row).
#     Num 4: Barcode 'ALHT WTD_102' with 'ALHT-Deer-10' names.
#     Kellen: "ALHT WTD_102: Ignore completely."
DROP_NUMS = {4}

# (c) Library typed under the wrong specimen -> re-attribute (Num -> correct Barcode).
LABEL_FIXES = {
    # Num 537: Barcode 'USGSTICKS-7304' but library named '7303-03'.
    #   Kellen: "USGSTICKS-7304: OK to return to 7303, seems right."
    537:'USGSTICKS-7303',

    # Num 247: Barcode column says USGSTICKS-LHT-Deer-091, but the library is
    #   named ALHT-Deer-09-16s and was sequenced before the LHT batch existed.
    #   Lawrence 2026-08-07: "Num 247 is most certainly ALHT-Deer-09 from the
    #   metadata." ALHT WTD_09 has ZERO run-data rows -> metadata comes from the
    #   registry (1 tick, WTD, accessioned 2024-04-23).
    247:'ALHT WTD_09',

    # Nums 617/618/619: Barcode says USGSTICKS-2913; every library on those rows
    #   is named USGSTICKS-2213-0*. Single-cell typo at barcoding.
    #   Lawrence 2026-08-07, corroborated by: barcode run order 2210/2211/2212/
    #   [2913]/2214/2215 (2213 is exactly the gap); Batch 09 extraction sheet
    #   wells B7/C7/D7 = USGSTICKS-2213-01/02/03, tick counts 2+1+3 = 6; and a
    #   real USGSTICKS-2913 existing separately at Num 500 (Effingham NH, 1 tick).
    #   SUPERSEDES Kellen's EQ02442213 (5 ticks, fawn, 2023-11-09) —
    #   wrong specimen. USGSTICKS-2213 has ZERO run-data rows -> registry supplies
    #   6 ticks, Brookfield NH, WTD M, collected 2024-11-14.
    617:'USGSTICKS-2213', 618:'USGSTICKS-2213', 619:'USGSTICKS-2213',
}

# (c2) SYSTEMATIC RULE — the ALHT / LHT barcode clobber.
#      Discovered 2026-08-07 by checking all 108 ALHT-Deer library numbers against
#      the specimen registry. 90 of them agree with the registry EXACTLY. Every one
#      of the 18 that disagree is a row whose Barcode column reads
#      'USGSTICKS-LHT-Deer-*' while the library is named 'ALHT-Deer-N'.
#
#      Cause: the LHT batch was accessioned 2025-02-18 and the Barcode column was
#      backfilled from that accession list, overwriting rows whose libraries had
#      already been sequenced 2024-12-19. Lawrence Gordon found this at Num 247 and
#      reasoned from exactly those dates; it is not a one-off, it is 17 specimens
#      across 26 library rows.
#
#      SAFETY (verified 2026-08-07):
#        - all 17 target 'ALHT WTD_N' specimens exist in the registry
#        - the rule keys on the LIBRARY NAME, not the barcode, so the 92 rows that
#          carry an LHT barcode WITH an LHT-named library (89 distinct genuine LHT
#          specimens) are untouched
#        - Num 4 is a different error (ALHT WTD_102 barcode, ALHT-Deer-10 library)
#          and is handled by DROP_NUMS
# 2026-09-18: hyphen-tolerant. 138 library names in the workbook are
# written ALHTDeer-NN or ALHT-DeerNN rather than ALHT-Deer-NN, and the
# strict form missed row Num 263 (ALHT-Deer143-2-16s), leaving one
# specimen's two libraries on two different specimens. Exactly one row
# changes; verified against the workbook 2026-09-18.
ALHT_LIB_RE = re.compile(r'^ALHT-?Deer-?(\d+)')
LHT_BC_PREFIX = 'USGSTICKS-LHT-Deer'

# (d) Specimens whose metadata exists in NEITHER sheet. Hand-entered from a
#     named source. Keyed on Num, applied after the barcode fixes.
MANUAL_SPECIMENS = {
    # Num 245: VDL submission, blank specimen block in the workbook.
    #   David Needle 2026-08-07: "24-3889 - Moose adult (head) cow, from Wilton, NH."
    #   In-region, so it counts. qPCR negative on all seven targets.
    245: {'bc':'24-3889', 'Host Spp':'Moose', 'Town':'Wilton',
          'State':'NH', 'Num Ticks':1},
}

# (e2) Specimens EXCLUDED FROM THE SINGLE-TICK SET only (they stay in the crosswalk
#      and in pooled/community analyses; they are simply not eligible to be called a
#      single tick). Two reasons, both decided 2026-08-07.
NOT_SINGLE_TICK = {
    # (i) Batch 12 = Cheryl Sullivan's UVM fungal-treatment trial. Every Batch-12
    #     single tick carries a treatment code (TL/TW/LC/WC/VTW); there are no
    #     untreated controls in that batch. Experimental animals do not belong in a
    #     wild-caught surveillance prevalence denominator. Confirms the 2026-07-13
    #     session-lock decision to route the fungal trial to the ITS arm.
    'USGSTICKS-VT-12-51TL','USGSTICKS-VT-12-53TL','USGSTICKS-VT-12-59TL1',
    'USGSTICKS-VT-12-59TL2','USGSTICKS-VT-12-57TW','USGSTICKS-VT-12-59TW2',
    'USGSTICKS-VT-12-51LC','USGSTICKS-VT-12-54WC','USGSTICKS-VT-12-VTW',
    # (ii) Contradicted tick count: run-data 'Num Ticks' says 1 but the batch
    #      extraction sheet's per-tube 'Tick count' says more. One tube each, no
    #      duplication, so this is a genuine source disagreement we could not
    #      resolve. Rule: if the collection record and the extraction record do not
    #      agree that the tube held one tick, it is not a single tick.
    'USGSTICKS-6503','USGSTICKS-6508','USGSTICKS-6510','USGSTICKS-9216',
}

# (e) Rows to EXCLUDE outright, with the reason on the record.
#     Both already fall out because their Barcode cell is blank; listed here so
#     the decision is visible rather than incidental.
EXCLUDED_NUMS = {
    244: "24-3775 — WTD from Martha's Vineyard MA (ascending-paralysis project, "
         "David Needle 2026-08-07). OUT OF STUDY REGION and not a randomly "
         "sampled animal. Excluded from all denominators. Its Amblyomma sp. "
         "detection is reported separately as a bounded species observation.",
    246: "4067 — David: not a VDL number. Kellen: 'shorthand for some other "
         "sample... too murky, I would exclude it.' Known discard.",
}

# ---- 1. trusted 16S anchor ------------------------------------------------
anchor = pd.read_csv(XW, sep="\t")
say(f"anchor (verified 16S crosswalk): {len(anchor)} rows")

# ---- 2. specimen REGISTRY --------------------------------------------------
reg = pd.read_excel(WB, sheet_name=REG, dtype=str)
reg.columns = [str(c).strip() for c in reg.columns]      # 'Barcode ' has a trailing space
reg['bc'] = reg['Barcode'].astype(str).str.strip()
reg = reg[~reg['bc'].str.lower().isin(['nan',''])]
REG_MAP = reg.set_index('bc')
say(f"specimen registry ('{REG}'): {len(REG_MAP)} specimens")
check(REG_MAP.index.is_unique, "registry barcodes are unique")

# ---- 3. workbook with all four marker Sample columns ------------------------
m = pd.read_excel(WB, sheet_name=SHEET)
m = m.rename(columns={'Sample':'name_16S','Sample.1':'name_18S',
                      'Sample.2':'name_cytb','Sample.3':'name_ITS'})
n_raw = len(m)

# apply manual specimens BEFORE the blank-barcode filter would drop them
for _num, _spec in MANUAL_SPECIMENS.items():
    sel = m['Num'] == _num
    check(sel.sum() == 1, f"manual specimen Num {_num} matches exactly one row")
    for _col, _val in _spec.items():
        if _col == 'bc':
            m.loc[sel, 'Barcode'] = _val
        else:
            m.loc[sel, _col] = _val
say(f"hand-entered {len(MANUAL_SPECIMENS)} specimen record(s): "
    f"Num {sorted(MANUAL_SPECIMENS)}")

for _num, _why in EXCLUDED_NUMS.items():
    say(f"excluded Num {_num}: {_why.splitlines()[0]}")
m = m[~m['Num'].isin(EXCLUDED_NUMS)]

m['bc'] = m['Barcode'].astype(str).str.strip()
n_blank = m['bc'].str.lower().isin(['nan','']).sum()
m = m[~m['bc'].str.lower().isin(['nan',''])]
m = m[~m['bc'].str.upper().str.startswith('AMP')]
ctrl = m['name_16S'].astype(str).str.upper().str.contains('BLANK|NTC|POS|CONTROL', na=False)
n_ctrl = int(ctrl.sum())
m = m[~ctrl]
say(f"workbook: {n_raw} rows -> {len(m)} specimen-library rows "
    f"({n_blank} blank-barcode, {n_ctrl} control)")

# ---- 4. APPLY the label fixes ----------------------------------------------
before = len(m)
m = m[~m['Num'].isin(DROP_NUMS)].copy()
say(f"dropped {before-len(m)} ignored row(s): Num {sorted(DROP_NUMS)}")
for _num, _correct_bc in LABEL_FIXES.items():
    sel = m['Num'] == _num
    check(sel.sum() == 1, f"LABEL_FIXES Num {_num} matches exactly one row")
    m.loc[sel, 'bc'] = _correct_bc
say(f"re-attributed {len(LABEL_FIXES)} library row(s) to correct specimens")

# apply the systematic ALHT/LHT rule
_alht_fixed = {}
for _i, _row in m.iterrows():
    _mt = ALHT_LIB_RE.match(str(_row['name_16S']))
    if _mt and str(_row['bc']).startswith(LHT_BC_PREFIX):
        _target = f"ALHT WTD_{_mt.group(1)}"
        m.at[_i, 'bc'] = _target
        _alht_fixed.setdefault(_target, 0)
        _alht_fixed[_target] += 1
say(f"ALHT/LHT barcode rule: {sum(_alht_fixed.values())} library row(s) across "
    f"{len(_alht_fixed)} specimen(s) re-attributed from an LHT barcode to ALHT WTD_N")
for _t in sorted(_alht_fixed): say(f"    {_t}  ({_alht_fixed[_t]} row(s))")

# ---- 5. pull specimen metadata from the registry where the sheet lacks it ---
# A re-attributed barcode may have no run-data row of its own, in which case the
# specimen fields still hold the OLD specimen's values. The registry is authoritative.
SPEC_COLS = {'Host Spp':'Host Spp', 'Town':'Town', 'State':'State', 'Num Ticks':'# of Ticks'}
n_filled = 0
for bc in sorted(set(LABEL_FIXES.values()) | set(_alht_fixed)):
    if bc not in REG_MAP.index:
        say(f"  note: {bc} not in the registry; keeping run-data values")
        continue
    sel = m['bc'] == bc
    for _sheet_col, _reg_col in SPEC_COLS.items():
        if _reg_col in REG_MAP.columns:
            m.loc[sel, _sheet_col] = REG_MAP.loc[bc, _reg_col]
    n_filled += 1
    say(f"  registry -> {bc}: {REG_MAP.loc[bc,'Host Spp']}, "
        f"{REG_MAP.loc[bc,'Town']} {REG_MAP.loc[bc,'State']}, "
        f"{REG_MAP.loc[bc,'# of Ticks']} tick(s)")
say(f"filled specimen metadata from the registry for {n_filled} barcode(s)")
say("  (the registry is the specimen authority: where both sources carry a tick")
say("   count they agree on 100% of specimens, checked 2026-08-07)")

# General pass: the registry is the specimen authority. Where a run-data row has a
# blank specimen field and the registry has a value, take the registry's. Safe:
# where BOTH sources carry a tick count they agree on 100% of specimens (checked
# 2026-08-07), so the registry is not in conflict with run-data, only fuller.
_gen = 0
for _col, _regcol in SPEC_COLS.items():
    if _regcol not in REG_MAP.columns:
        continue
    _blank = m[_col].isna() | m[_col].astype(str).str.strip().str.lower().isin(['nan','','none'])
    for _bc in sorted(set(m.loc[_blank, 'bc']) & set(REG_MAP.index)):
        _val = REG_MAP.loc[_bc, _regcol]
        if str(_val).strip().lower() not in ('nan','','none'):
            m.loc[_blank & (m['bc'] == _bc), _col] = _val
            _gen += 1
say(f"general registry backfill: {_gen} blank specimen field(s) filled")

# ---- 6. residual mismatch check (should be empty) --------------------------
def num_of(s):
    s=str(s); mt=re.search(r'(\d+)', s.replace('16s','').replace('18s','').replace('cyto','').replace('ITS',''))
    return mt.group(1) if mt else ''
m['bc_num']=m['bc'].map(num_of); m['s16_num']=m['name_16S'].map(num_of)
m['mismatch']=(m['bc_num']!='') & (m['s16_num']!='') & (m['bc_num']!=m['s16_num'])
mism = m[m['mismatch']][['Num','bc','name_16S','name_18S','name_cytb','name_ITS','Num Ticks','Host Spp']]
mism.to_csv("results/metadata/crosswalk_label_mismatches.tsv", sep="\t", index=False)
say(f"\nresidual label mismatches: {len(mism)}")
say("  (these are the naive number-check's FALSE ALARMS on intentional merges "
    "and re-attributions, e.g. EQ024'42213' contains '2213'. Inspect, don't panic.)")
if len(mism): say(mism.to_string(index=False))

# ---- 7. per-specimen aggregation on Barcode --------------------------------
m['specimen'] = m['bc'].replace(merge_map)

# ---- 7b. LIBRARY LOOKUP: emit the resolution instead of discarding it -------
# Every surface form of every specimen name, already split, casefolded and
# stripped, one row per form. No consumer should ever parse a name again.
#
# This adds NO new naming knowledge. Each key comes from a name already in the
# workbook, normalised by rules this project already states:
#   - pipeline_config.yml  samples.id_regex "^([^.]+)"  strips the run-id suffix
#   - TICK_GROUND_TRUTH section 3 requires split-on-semicolon, casefold, and a
#     REGEX marker-tag strip rather than a fixed suffix list
#   - merge_map above collapses ALHT-Deer-N into ALHT WTD_N, and the ALHT/LHT
#     rule re-attributes ALHT-Deer-N libraries to ALHT WTD_N specimens, so this
#     script has already asserted that the two forms name one specimen
LOOKUP_MARKER_TAG = re.compile(r'-(16s|18s|cytob|cb|its)([A-Za-z]*\d*)(-\d+)?$',
                               re.I)
LOOKUP_RUN_ID = re.compile(r'\..*$')
LOOKUP_ALHT_ALIAS = re.compile(r'^alht-(?:wtd|deer)-(\d+)$')
LOOKUP_ALHT_SHAPE = re.compile(r'^alht-?deer-?(\d.*)$')


def lookup_key_forms(raw_name):
    """I return every normalised key one library or barcode name can present.

    I apply the transformations cumulatively because the workbook uses all of
    these shapes for the same specimen, and a join that handles only some of
    them loses specimens silently rather than raising.
    """
    text = str(raw_name).strip()
    if text.lower() in ('nan', '', 'none'):
        return set()
    forms = set()
    for variant in (text, LOOKUP_RUN_ID.sub('', text)):
        flat = re.sub(r'[\s_]+', '-', variant).lower().strip('-')
        if not flat:
            continue
        forms.add(flat)
        tag = LOOKUP_MARKER_TAG.search(flat)
        if tag and flat[:tag.start()]:
            forms.add(flat[:tag.start()])
    # Missing-hyphen shape first: alhtdeer-02 and alht-deer02 become
    # alht-deer-02, so all three spellings are one key. 138 library names carry
    # it. Then the WTD/Deer alias, which needs the hyphenated form to match.
    expanded = set(forms)
    for form in list(expanded):
        shape = LOOKUP_ALHT_SHAPE.match(form)
        if shape:
            expanded.add('alht-deer-' + shape.group(1))
    for form in list(expanded):
        alias = LOOKUP_ALHT_ALIAS.match(form)
        if alias:
            expanded.add('alht-wtd-' + alias.group(1))
            expanded.add('alht-deer-' + alias.group(1))
    return expanded


LOOKUP_NAME_COLS = {'16S': 'name_16S', '18S': 'name_18S',
                    'cytb': 'name_cytb', 'ITS': 'name_ITS'}
# The run sheet carries one Run Batch column per marker, in marker order, so
# pandas suffixes the repeats. Batch is a property of a LIBRARY, not of a
# specimen: two libraries of one specimen can run on different plates, which is
# why this sits beside the library key and not in the crosswalk row.
LOOKUP_BATCH_COLS = {'16S': 'Run Batch', '18S': 'Run Batch.1',
                     'cytb': 'Run Batch.2', 'ITS': 'Run Batch.3'}
_lookup = {}          # (marker, key) -> specimen
_lookup_batch = {}    # (marker, key) -> run batch
_lookup_clashes = []

def _lookup_add(marker, key, specimen, batch=None):
    previous = _lookup.get((marker, key))
    if previous is not None and previous != specimen:
        _lookup_clashes.append((marker, key, previous, specimen))
    else:
        _lookup[(marker, key)] = specimen
    if batch is not None and str(batch).strip().lower() not in ('nan', '',
                                                               'none'):
        _lookup_batch.setdefault((marker, key), str(batch).strip())

for _i, _row in m.iterrows():
    _specimen = _row['specimen']
    # 'specimen' is the barcode namespace: barcodes, and the merge_map
    # pre-images that resolve to them. This is the namespace the qPCR sheet
    # and the accessioning records key on.
    for _form in lookup_key_forms(_row['bc']):
        _lookup_add('specimen', _form, _specimen)
    for _form in lookup_key_forms(_specimen):
        _lookup_add('specimen', _form, _specimen)
    for _marker, _col in LOOKUP_NAME_COLS.items():
        _batch = _row.get(LOOKUP_BATCH_COLS.get(_marker))
        for _piece in str(_row.get(_col, '')).split(';'):
            for _form in lookup_key_forms(_piece):
                _lookup_add(_marker, _form, _specimen, _batch)
                # a library name is also a valid specimen-namespace alias
                _lookup_add('specimen', _form, _specimen)

for _pre, _post in merge_map.items():
    for _form in lookup_key_forms(_pre):
        _lookup_add('specimen', _form, _post)

say("")
say("=== library lookup ===")
say("  {} (marker, key) -> specimen entries".format(len(_lookup)))
check(not _lookup_clashes,
      "no lookup key maps to two specimens ({} clash(es))"
      .format(len(_lookup_clashes)))
if _lookup_clashes:
    for _c in _lookup_clashes[:20]:
        say("    CLASH {} '{}' -> {} and {}".format(*_c))
def agg(s):
    v=[str(x).strip() for x in s if str(x).strip().lower() not in ('nan','','none')]
    return ';'.join(sorted(set(v)))
g=m.groupby('specimen')
x=pd.DataFrame({
  'barcode':g['bc'].first(),
  'nums':g['Num'].apply(lambda s:';'.join(sorted(set(str(int(v)) for v in s if pd.notna(v))))),
  'host_spp':g['Host Spp'].first(),'state':g['State'].first(),'num_ticks':g['Num Ticks'].first(),
  'name_16S':g['name_16S'].apply(agg),'name_18S':g['name_18S'].apply(agg),
  'name_cytb':g['name_cytb'].apply(agg),'name_ITS':g['name_ITS'].apply(agg),
}).reset_index().rename(columns={'specimen':'sample_id'})

def hb(v):
    v=str(v).lower()
    if any(k in v for k in ['deer','wtd','odocoile']): return 'Deer'
    if 'moose' in v or 'alces' in v: return 'Moose'
    if 'bear' in v or 'ursus' in v: return 'Bear'
    if v in ('nan','','none'): return 'Unknown'
    return 'Other wildlife'
x['host_bin']=x['host_spp'].map(hb)
x['state']=x['state'].astype(str).str.strip().str.upper().replace(
    {'NEW HAMPSHIRE':'NH','VERMONT':'VT','MAINE':'ME','NAN':''})
x['num_ticks']=pd.to_numeric(x['num_ticks'],errors='coerce')
x['is_single_tick']=(x['num_ticks']==1) & (~x['sample_id'].isin(NOT_SINGLE_TICK))
_blocked = sorted(set(x.loc[x['num_ticks']==1,'sample_id']) & NOT_SINGLE_TICK)
say(f"\nsingle-tick exclusions applied: {len(_blocked)} specimen(s) had "
    f"num_ticks==1 but are not eligible")
for _b in _blocked: say(f"    {_b}")

# ---- 8. ASSERTIONS ---------------------------------------------------------
say("\n=== build checks ===")
check(x['sample_id'].is_unique, "one row per specimen (sample_id unique)")
_nan = x[x['num_ticks'].isna()]
say(f"  WARN {len(_nan)} specimen(s) have NO tick count -> excluded from the "
    f"single-tick set by rule (see docstring note 4)")
if len(_nan):
    for _b in sorted(_nan['sample_id'])[:40]:
        say(f"         {_b}")
    _nan[['sample_id','barcode','host_spp','state']].to_csv(
        "results/metadata/crosswalk_no_tick_count.tsv", sep="\t", index=False)
    say("       -> results/metadata/crosswalk_no_tick_count.tsv")
_nohost = x[x['host_bin']=='Unknown']
if len(_nohost):
    say(f"  WARN {len(_nohost)} specimen(s) have NO collection host. They stay in the "
        f"crosswalk but are dropped by build_analysis_table.py (line 82) and cannot "
        f"enter any host-stratified analysis. Same rule as the tick counts: do not guess.")
    for _b in sorted(_nohost['sample_id'])[:40]:
        say(f"         {_b}")
    _nohost[['sample_id','barcode','state','num_ticks']].to_csv(
        "results/metadata/crosswalk_no_host.tsv", sep="\t", index=False)
    say("       -> results/metadata/crosswalk_no_host.tsv")
else:
    say("  OK   every specimen has a host")

# 2026-09-18: one specimen's libraries must not split across two specimens.
# Rows Num 262 and 263 are the same animal and differ only by a missing hyphen
# in two library names; before the hyphen-tolerant ALHT_LIB_RE, row 263 stayed
# on USGSTICKS-LHT-Deer-143 while row 262 moved to ALHT WTD_143.
_143 = x[x['sample_id'] == 'ALHT WTD_143']
if len(_143):
    _names = str(_143['name_16S'].iloc[0])
    check('143-1' in _names and '143-2' in _names,
          "ALHT WTD_143 holds BOTH 16S libraries (-143-1 and -143-2); got: "
          + _names)
else:
    say("  WARN ALHT WTD_143 is not in the crosswalk; skipping the "
        "split-specimen check")

for _bc in ('ALHT WTD_09','USGSTICKS-2213'):
    row = x[x['sample_id']==_bc]
    check(len(row)==1, f"{_bc} is present as its own specimen")
check(int(x.loc[x['sample_id']=='ALHT WTD_09','num_ticks'].iloc[0])==1,
      "ALHT WTD_09 is a SINGLE TICK (was wrongly 18 before this fix)")
check(int(x.loc[x['sample_id']=='USGSTICKS-2213','num_ticks'].iloc[0])==6,
      "USGSTICKS-2213 has 6 ticks (was wrongly EQ02442213's 5 before this fix)")
check('USGSTICKS-LHT-Deer-091' in set(x['sample_id']),
      "USGSTICKS-LHT-Deer-091 survives as its own specimen (unmerged)")

if 'Barcode' in anchor.columns:
    a=set(anchor['Barcode'].astype(str).str.strip()); b=set(x['barcode'].astype(str).str.strip())
    say(f"  cross-check vs verified 16S: {len(a&b)} shared; {len(a-b)} anchor-only; {len(b-a)} here-only")
    check(len(a&b)/max(len(a),1) > 0.95,
          f"anchor coverage {100*len(a&b)/max(len(a),1):.1f}% (>95%)")

x=x[['sample_id','barcode','nums','host_spp','host_bin','state','num_ticks',
     'is_single_tick','name_16S','name_18S','name_cytb','name_ITS']]
x.to_csv("results/metadata/crosswalk_4marker.tsv", sep="\t", index=False)

# ---- 8b. write the lookup and verify it against the qPCR sheet -------------
_lookup_rows = pd.DataFrame(
    [{'marker': _k[0], 'key': _k[1], 'specimen': _v,
      'batch': _lookup_batch.get(_k, '')}
     for _k, _v in sorted(_lookup.items())])
_lookup_rows.to_csv("results/metadata/library_lookup.tsv", sep="\t",
                    index=False)
say("")
say("wrote results/metadata/library_lookup.tsv  ({} rows)"
    .format(len(_lookup_rows)))
say("  per marker: {}".format(
    _lookup_rows['marker'].value_counts().to_dict()))
_with_batch = int((_lookup_rows['batch'].astype(str).str.strip() != '').sum())
say("  keys carrying a run batch: {} of {}".format(_with_batch,
                                                   len(_lookup_rows)))
_batches = sorted(set(b for b in _lookup_rows['batch'].astype(str)
                      if b.strip()))
say("  distinct run batches: {}".format(len(_batches)))
for _b in _batches:
    say("      " + _b)
say("  Batch is a LIBRARY property. A specimen with libraries on two plates")
say("  has two batches, which is the point: it is what lets a co-occurrence")
say("  be conditioned on batch.")
say("  consumers: join on (marker, key) after lowercasing the incoming name.")
say("  A key that is absent returns nothing, which is visible. Do NOT add a")
say("  fallback that strips characters until something matches.")

# The qPCR sheet is the join this crosswalk never covered, and it is where
# build_concordance.py lost six single ticks. Verify it here, at build time.
_specimen_keys = {_k[1]: _v for _k, _v in _lookup.items()
                  if _k[0] == 'specimen'}
try:
    _q = pd.read_excel(WB, sheet_name="qPCR Results", dtype=str)
except ValueError as _err:
    say("  WARN could not read the 'qPCR Results' sheet: " + str(_err))
    _q = None

if _q is not None:
    def _norm_one(value):
        forms = lookup_key_forms(value)
        for form in forms:
            if form in _specimen_keys:
                return _specimen_keys[form]
        return None

    _resolved, _unresolved = 0, []
    for _j, _qrow in _q.iterrows():
        _hit = _norm_one(_qrow.get('host_id')) or _norm_one(_qrow.get('Barcode'))
        if _hit:
            _resolved += 1
        else:
            _unresolved.append(str(_qrow.get('host_id')))
    say("  qPCR sheet: {} of {} rows resolve to a crosswalk specimen"
        .format(_resolved, len(_q)))
    if _unresolved:
        say("  WARN {} qPCR row(s) do not resolve; first 15:"
            .format(len(_unresolved)))
        for _u in _unresolved[:15]:
            say("         " + _u)

    _single = set(x.loc[x['is_single_tick'], 'sample_id'])
    _qspecimens = set()
    for _j, _qrow in _q.iterrows():
        _hit = _norm_one(_qrow.get('host_id')) or _norm_one(_qrow.get('Barcode'))
        if _hit:
            _qspecimens.add(_hit)
    _paired = _single & _qspecimens
    say("  single ticks present in the qPCR sheet: {} of {}"
        .format(len(_paired), len(_single)))
    # I compare against the cohort size rather than a hardcoded 93. GOLDEN in
    # build_single_tick_v2.py is what a hardcoded cohort constant becomes: it
    # still holds the n=103 values and has been failing silently for weeks.
    check(len(_paired) == len(_single),
          "every single tick resolves against the qPCR sheet "
          "(build_concordance.py resolved only 83 before this lookup existed)")
    _subsampled = set()
    for _j, _qrow in _q.iterrows():
        if str(_qrow.get('subsample_count')).strip() != '1':
            continue
        _hit = _norm_one(_qrow.get('host_id')) or _norm_one(_qrow.get('Barcode'))
        if _hit:
            _subsampled.add(_hit)
    say("  CONCORDANCE COHORT: {} single ticks have a one-subsample qPCR row."
        .format(len(_single & _subsampled)))
    say("  The remainder resolve but their only qPCR row is a pool, so they")
    say("  cannot be paired per tick. This number, not 83, is the n for")
    say("  build_concordance.py and for Results 3.4.")

say(f"\nwrote results/metadata/crosswalk_4marker.tsv  ({len(x)} specimens)")
say(f"  total ticks:  {int(x['num_ticks'].sum())}")
say(f"  single ticks: {int(x['is_single_tick'].sum())}  "
    f"(was 103 before the 2026-08-07 exclusions)")
say(f"  host_bin: {x['host_bin'].value_counts().to_dict()}")
say(f"  state:    {x['state'].value_counts().to_dict()}")
say("\nNOTE: specimen / tick / single-tick counts have CHANGED from the pre-2026-08-07")
say("crosswalk (was 544 / 3703 / 103). Propagate the new numbers everywhere; do not")
say("assume the old ones survive. See TICK_GROUND_TRUTH sections 1.6 and 3.4.")
open("results/metadata/crosswalk_build_report.txt","w").write("\n".join(_report))
say("\nbuild report -> results/metadata/crosswalk_build_report.txt")
