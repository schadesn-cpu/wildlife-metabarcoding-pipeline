#!/usr/bin/env python3
"""patch_add_library_lookup.py

I add a library-lookup output to build_4marker_crosswalk.py so the crosswalk
stops discarding the name resolution it already performs.

WHY
---
build_4marker_crosswalk.py works out, with assertions and a documented
decision trail, that ALHT-Deer-09 and ALHT WTD_09 are one specimen, that
Nums 617/618/619 belong to USGSTICKS-2213, and which of 108 ALHT library
names were clobbered by the LHT accession backfill. It then writes only the
specimen-level summary, with the library names aggregated back into
semicolon-delimited raw strings. Every consumer re-derives the resolution
from those strings, fifteen times, differently. build_concordance.py does
not know the ALHT WTD alias and loses six single ticks because of it.

This patch writes the resolution out. It adds no new naming rules: every key
it emits comes from a name already present in the workbook, normalised by
the transformations this script and pipeline_config.yml already specify.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_add_library_lookup.py

I patch scripts/build_4marker_crosswalk.py in place, refuse to run twice,
and refuse to guess if either anchor line has moved.
"""
import os
import sys

TARGET = "scripts/build_4marker_crosswalk.py"

ANCHOR_BUILD = "m['specimen'] = m['bc'].replace(merge_map)\n"
ANCHOR_WRITE = ('x.to_csv("results/metadata/crosswalk_4marker.tsv", '
                'sep="\\t", index=False)\n')

BLOCK_BUILD = '''
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
LOOKUP_MARKER_TAG = re.compile(r'-(16s|18s|cytob|cb|its)([A-Za-z]*\\d*)(-\\d+)?$',
                               re.I)
LOOKUP_RUN_ID = re.compile(r'\\..*$')
LOOKUP_ALHT_ALIAS = re.compile(r'^alht-(?:wtd|deer)-(\\d+)$')


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
        flat = re.sub(r'[\\s_]+', '-', variant).lower().strip('-')
        if not flat:
            continue
        forms.add(flat)
        tag = LOOKUP_MARKER_TAG.search(flat)
        if tag and flat[:tag.start()]:
            forms.add(flat[:tag.start()])
    expanded = set(forms)
    for form in forms:
        alias = LOOKUP_ALHT_ALIAS.match(form)
        if alias:
            expanded.add('alht-wtd-' + alias.group(1))
            expanded.add('alht-deer-' + alias.group(1))
    return expanded


LOOKUP_NAME_COLS = {'16S': 'name_16S', '18S': 'name_18S',
                    'cytb': 'name_cytb', 'ITS': 'name_ITS'}
_lookup = {}          # (marker, key) -> specimen
_lookup_clashes = []

def _lookup_add(marker, key, specimen):
    previous = _lookup.get((marker, key))
    if previous is not None and previous != specimen:
        _lookup_clashes.append((marker, key, previous, specimen))
    else:
        _lookup[(marker, key)] = specimen

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
        for _piece in str(_row.get(_col, '')).split(';'):
            for _form in lookup_key_forms(_piece):
                _lookup_add(_marker, _form, _specimen)
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
'''

BLOCK_WRITE = '''
# ---- 8b. write the lookup and verify it against the qPCR sheet -------------
_lookup_rows = pd.DataFrame(
    [{'marker': _k[0], 'key': _k[1], 'specimen': _v}
     for _k, _v in sorted(_lookup.items())])
_lookup_rows.to_csv("results/metadata/library_lookup.tsv", sep="\\t",
                    index=False)
say("")
say("wrote results/metadata/library_lookup.tsv  ({} rows)"
    .format(len(_lookup_rows)))
say("  per marker: {}".format(
    _lookup_rows['marker'].value_counts().to_dict()))
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
'''


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "library_lookup.tsv" in source:
        sys.exit("already patched: " + TARGET + " mentions library_lookup.tsv. "
                 "Nothing changed.")
    for anchor, label in ((ANCHOR_BUILD, "the specimen-assignment line"),
                          (ANCHOR_WRITE, "the crosswalk write line")):
        if source.count(anchor) != 1:
            sys.exit("expected exactly one " + label + ", found "
                     + str(source.count(anchor))
                     + ". The script has changed since this patch was written; "
                     "patch by hand rather than trusting this.")

    patched = source.replace(ANCHOR_BUILD, ANCHOR_BUILD + BLOCK_BUILD, 1)
    patched = patched.replace(ANCHOR_WRITE, ANCHOR_WRITE + BLOCK_WRITE, 1)

    compile(patched, TARGET, "exec")        # syntax-check before writing
    with open(TARGET, "w", encoding="utf-8") as handle:
        handle.write(patched)
    print("patched " + TARGET)
    print("added " + str(len(BLOCK_BUILD.splitlines())
                         + len(BLOCK_WRITE.splitlines())) + " lines")
    print("now run: python3 " + TARGET)


if __name__ == "__main__":
    main()
