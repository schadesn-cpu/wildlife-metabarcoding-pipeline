#!/usr/bin/env python3
"""patch_record_taxonomy.py

Makes check_batch_arms.read_taxonomy() record which taxonomy file it chose.

WHY
---
read_taxonomy() picks a taxonomy at run time by whichever file covers the most
features, prints the choice at build_core_microbiome.py line 186, and nothing
captures that line. So the provenance of Sections 3.3 and 3.4 is a filesystem
search whose outcome exists only in a terminal that has since closed.

Section 2.4 says each call rests on a BLAST-verified ASV. That is true today --
check_arms_resolver.py proved taxonomy_16S_droplist_applied.tsv always wins --
but it is true by accident of file ordering, and a renamed column or a new
results/**/taxonomy*16S*.tsv would change it silently. This makes the choice a
record instead of a reconstruction.

WHAT IT CHANGES
---------------
One append-only log line per call, written to
results/16S/all/taxonomy_verified/TAXONOMY_USED.txt:

    2026-10-20T09:14:22    results/.../taxonomy_16S_droplist_applied.tsv    89970 of 89970 features

Wrapped in try/except OSError, so a read-only filesystem or a missing directory
can never break a pipeline run. Nothing else in the function changes: same
return value, same selection logic, same behaviour.

SAFETY
------
I refuse rather than guess:
  * the anchor must appear exactly once, or I stop
  * I write a .bak first
  * I parse the result with ast before saving, and discard it if it does not
    compile
  * I print a diff and change nothing unless you pass --write

Usage:
    python3 scripts/patch_record_taxonomy.py            # show the diff
    python3 scripts/patch_record_taxonomy.py --write    # apply it
"""
import ast
import difflib
import os
import shutil
import sys

TARGET = "scripts/check_batch_arms.py"
ANCHOR = "    if best is None or best[1] == 0:"
LOGLINE = '''    try:
        os.makedirs("results/16S/all/taxonomy_verified", exist_ok=True)
        with open("results/16S/all/taxonomy_verified/TAXONOMY_USED.txt", "a") as _fh:
            _fh.write("%s\\t%s\\t%d of %d features\\n" % (
                _dt.datetime.now().isoformat(timespec="seconds"),
                best[0] if best else "NONE",
                best[1] if best else 0, len(features)))
    except OSError:
        pass
'''


def main():
    root = os.path.expanduser("~/meedlab/tick_project")
    os.chdir(root)
    if not os.path.isfile(TARGET):
        sys.exit("cannot find %s under %s" % (TARGET, root))

    with open(TARGET, encoding="utf-8") as fh:
        original = fh.read()

    if "TAXONOMY_USED.txt" in original:
        print("Already patched -- TAXONOMY_USED.txt appears in the file.")
        print("Nothing to do.")
        return

    if original.count(ANCHOR) != 1:
        sys.exit("anchor appears %d times, expected 1. Stopping rather than "
                 "guessing:\n  %s" % (original.count(ANCHOR), ANCHOR))

    patched = original.replace(ANCHOR, LOGLINE + ANCHOR)

    # datetime under an alias, so a plain `import datetime` elsewhere is untouched
    if "import datetime as _dt" not in patched:
        lines = patched.split("\n")
        last_import = max((i for i, l in enumerate(lines)
                           if l.startswith(("import ", "from "))), default=0)
        lines.insert(last_import + 1, "import datetime as _dt")
        patched = "\n".join(lines)

    if "\nimport os" not in patched and not patched.startswith("import os"):
        sys.exit("this module does not import os, which the patch needs. "
                 "Add `import os` by hand and rerun.")

    try:
        ast.parse(patched)
    except SyntaxError as exc:
        sys.exit("the patched file does not parse (%s). Nothing written." % exc)

    diff = list(difflib.unified_diff(
        original.split("\n"), patched.split("\n"),
        fromfile=TARGET, tofile=TARGET + " (patched)", lineterm="", n=3))
    print("\n".join(diff))
    print()

    if "--write" not in sys.argv:
        print("Dry run. Re-run with --write to apply.")
        return

    shutil.copy2(TARGET, TARGET + ".bak")
    with open(TARGET, "w", encoding="utf-8") as fh:
        fh.write(patched)
    print("written. backup at %s.bak" % TARGET)
    print()
    print("Next time any script calls read_taxonomy(), the choice lands in")
    print("results/16S/all/taxonomy_verified/TAXONOMY_USED.txt")


if __name__ == "__main__":
    main()
