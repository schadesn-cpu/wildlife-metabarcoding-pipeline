#!/usr/bin/env python3
"""patch_add_batch_to_lookup.py

I add the sequencing run batch to the library lookup, and give tick_lookup a
way to read it, so the pooled co-occurrence can be conditioned on batch.

WHY
---
The pooled Anaplasma x Borreliella association (OR 2.03, q = 0.0002) survives
conditioning on pool size (Mantel-Haenszel 2.05) and on read depth (1.89), and
has been stable across the anchoring fix, the join fix and two crosswalk
rebuilds: 2.48 crude and 2.27 depth-adjusted in August on 558 samples, against
2.03 and 1.89 today on 602.

Batch is the confounder it has NOT been conditioned on, and it is the largest
effect in this dataset: R2 10.9%, p = 0.001, with 530 of 551 animals sitting
entirely within one batch (ground truth 2b, lens 2 constraints). Until the
association is tested against batch it cannot be attributed to anything.

Batch is a property of a LIBRARY, not of a specimen, so it belongs beside the
library key rather than in the crosswalk's per-specimen row. The run sheet
carries one per marker: 'Run Batch' is 16S, '.1' is 18S, '.2' is cytb, '.3' is
ITS.

WHAT I CHANGE
-------------
1. build_4marker_crosswalk.py emits a fourth column, batch, in
   library_lookup.tsv. Existing readers are unaffected: tick_lookup reads
   columns by name.

2. scripts/tick_lookup.py gains load_batches() and resolve_batch(), written
   the same way as load() and resolve(): a miss returns None and nothing
   guesses.

I do NOT touch build_analysis_table.py here. Adding the column and consuming
it are separate changes, so if the batch values turn out to be messy the
lookup can be inspected before anything downstream moves.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_add_batch_to_lookup.py
    python3 scripts/build_4marker_crosswalk.py
"""
import os
import sys

TARGET = "scripts/build_4marker_crosswalk.py"
LOOKUP_MODULE = "scripts/tick_lookup.py"

OLD_HEAD = """LOOKUP_NAME_COLS = {'16S': 'name_16S', '18S': 'name_18S',
                    'cytb': 'name_cytb', 'ITS': 'name_ITS'}
_lookup = {}          # (marker, key) -> specimen
_lookup_clashes = []

def _lookup_add(marker, key, specimen):
    previous = _lookup.get((marker, key))
    if previous is not None and previous != specimen:
        _lookup_clashes.append((marker, key, previous, specimen))
    else:
        _lookup[(marker, key)] = specimen
"""

NEW_HEAD = """LOOKUP_NAME_COLS = {'16S': 'name_16S', '18S': 'name_18S',
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
"""

OLD_LOOP = """    for _marker, _col in LOOKUP_NAME_COLS.items():
        for _piece in str(_row.get(_col, '')).split(';'):
            for _form in lookup_key_forms(_piece):
                _lookup_add(_marker, _form, _specimen)
                # a library name is also a valid specimen-namespace alias
                _lookup_add('specimen', _form, _specimen)"""

NEW_LOOP = """    for _marker, _col in LOOKUP_NAME_COLS.items():
        _batch = _row.get(LOOKUP_BATCH_COLS.get(_marker))
        for _piece in str(_row.get(_col, '')).split(';'):
            for _form in lookup_key_forms(_piece):
                _lookup_add(_marker, _form, _specimen, _batch)
                # a library name is also a valid specimen-namespace alias
                _lookup_add('specimen', _form, _specimen)"""

OLD_WRITE = """_lookup_rows = pd.DataFrame(
    [{'marker': _k[0], 'key': _k[1], 'specimen': _v}
     for _k, _v in sorted(_lookup.items())])"""

NEW_WRITE = """_lookup_rows = pd.DataFrame(
    [{'marker': _k[0], 'key': _k[1], 'specimen': _v,
      'batch': _lookup_batch.get(_k, '')}
     for _k, _v in sorted(_lookup.items())])"""

OLD_REPORT = """say("  per marker: {}".format(
    _lookup_rows['marker'].value_counts().to_dict()))"""

NEW_REPORT = """say("  per marker: {}".format(
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
say("  be conditioned on batch.")"""

REPLACEMENTS = [
    (OLD_HEAD, NEW_HEAD, "the lookup header and _lookup_add"),
    (OLD_LOOP, NEW_LOOP, "the per-marker key loop"),
    (OLD_WRITE, NEW_WRITE, "the lookup write"),
    (OLD_REPORT, NEW_REPORT, "the per-marker report line"),
]

TICK_LOOKUP_ADDITION = '''

def load_batches(path=LOOKUP_PATH):
    """I load the (marker, key) -> run batch map from the same lookup file.

    Batch was added on 2026-09-18 so the pooled co-occurrence could be
    conditioned on it. A lookup built before that has no batch column, and I
    say so rather than returning an empty map that looks like "no batches".
    """
    if not os.path.isfile(path):
        sys.exit("library lookup not found: " + path
                 + "\\nBuild it first, from the project root:"
                 + "\\n    python3 scripts/build_4marker_crosswalk.py")
    frame = pd.read_csv(path, sep="\\t", dtype=str)
    if "batch" not in frame.columns:
        sys.exit(path + " has no 'batch' column. It was written before "
                 "2026-09-18. Rebuild the crosswalk:\\n"
                 "    python3 scripts/build_4marker_crosswalk.py")
    table = {}
    for marker, key, batch in zip(frame["marker"], frame["key"],
                                  frame["batch"]):
        if isinstance(batch, str) and batch.strip():
            table[(str(marker), str(key))] = batch.strip()
    return table


def resolve_batch(batches, marker, value):
    """I return the run batch for one library name, or None if it is absent."""
    if marker not in MARKERS:
        sys.exit("unknown marker '" + str(marker) + "'. Use one of: "
                 + ", ".join(MARKERS))
    key = normalise(value)
    if not key:
        return None
    return batches.get((marker, key))
'''


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    if not os.path.isfile(LOOKUP_MODULE):
        sys.exit("not found: " + LOOKUP_MODULE)

    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()
    if "LOOKUP_BATCH_COLS" in source:
        sys.exit("already patched: " + TARGET + " defines LOOKUP_BATCH_COLS. "
                 "Nothing changed.")
    if "library_lookup.tsv" not in source:
        sys.exit(TARGET + " has no library lookup. Run "
                 "patch_add_library_lookup.py first.")
    for old, new, label in REPLACEMENTS:
        count = source.count(old)
        if count != 1:
            sys.exit("expected exactly one occurrence of " + label
                     + ", found " + str(count)
                     + ". Patch by hand rather than trusting this.")
        source = source.replace(old, new, 1)
    compile(source, TARGET, "exec")

    with open(LOOKUP_MODULE, encoding="utf-8") as handle:
        module = handle.read()
    if "def load_batches(" in module:
        print("note: " + LOOKUP_MODULE + " already has load_batches; "
              "leaving it alone")
    else:
        module = module.rstrip("\n") + TICK_LOOKUP_ADDITION
        compile(module, LOOKUP_MODULE, "exec")
        with open(LOOKUP_MODULE, "w", encoding="utf-8") as handle:
            handle.write(module)
        print("extended " + LOOKUP_MODULE + " with load_batches / resolve_batch")

    with open(TARGET, "w", encoding="utf-8") as handle:
        handle.write(source)
    print("patched " + TARGET)
    print()
    print("Expect on the next crosswalk build: a 'batch' column in")
    print("library_lookup.tsv, a count of keys carrying one, and the list of")
    print("distinct run batches. 563 / 3,939 / 93 must not move.")
    print()
    print("run: python3 scripts/build_4marker_crosswalk.py")


if __name__ == "__main__":
    main()
