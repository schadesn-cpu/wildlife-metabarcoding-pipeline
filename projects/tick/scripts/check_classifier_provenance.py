#!/usr/bin/env python3
"""check_classifier_provenance.py

I read QIIME 2 provenance to say which classifier actually produced each
taxonomy artifact. I write nothing.

WHY
---
Section 2.3 carries three markers that rest on the same question, and the
classifiers/ listing makes them sharper rather than settling them:

  16S  The text says SILVA 138. classifiers/ holds no V4-V5 classifier at all,
       only silva-138-99-nb-classifier-515-806.qza (V4), the full-length
       silva-138-99-nb-classifier.qza, and silva-138-V9-classifier.qza. The
       reads are V4-V5 at 373 to 375 bp, so whichever was used was trained on
       a different fragment and the text has to say which.

  18S  The text says PR2 v5.0.0 was "trained on full-length SSU sequences
       rather than extracted to the amplicon". classifiers/ holds
       pr2-18S-classifier-trimmed.qza beside pr2-18S-seqs-trimmed.qza, and a
       trimmed classifier is the opposite of that claim.

A file listing cannot answer either. Provenance can: every QIIME 2 artifact
records the action that made it and the UUIDs of that action's inputs, so the
taxonomy artifact names its own classifier and nobody has to remember.

HOW
---
Each .qza is a zip. Its metadata.yaml carries the artifact's own uuid and
type; its provenance/action/action.yaml records the action and its inputs. I
read every artifact's uuid once, then every taxonomy artifact's recorded
inputs, and match the two. A uuid I cannot match is reported as unmatched
rather than guessed at.

FIXED SINCE THE FIRST VERSION
-----------------------------
The first version opened each zip inside a `with` block and returned the
handle after the block had closed it, so the first read raised "Attempt to use
ZIP archive that was already closed". Every read now happens inside the block
that opened it, and each artifact is identified once rather than three times.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_classifier_provenance.py
"""
import os
import re
import sys
import zipfile

ROOTS = ["results", "classifiers", "qiime2"]
UUID = re.compile(r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-'
                  r'[0-9a-f]{4}-[0-9a-f]{12}')


def slurp(path, wanted):
    """I open the zip, read the first member `wanted` accepts, and close it."""
    try:
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist() if wanted(n)]
            if not names:
                return None
            with z.open(sorted(names, key=len)[0]) as fh:
                return fh.read().decode("utf-8", "replace")
    except Exception as exc:                              # noqa: BLE001
        return "<<unreadable: %s>>" % exc


def field(text, key):
    # [^\S\n] rather than \s so the value cannot be pulled off the NEXT line:
    # "action:\n    type: method" otherwise reports the action as "type: method".
    m = re.search(r'^[^\S\n]*%s:[^\S\n]*(\S.*?)[^\S\n]*$' % key,
                  text or "", flags=re.M)
    return m.group(1).strip().strip("'\"") if m else None


def action_name(text):
    """I return the nested action name, e.g. classify_sklearn."""
    m = re.search(r'^[^\S\n]+action:[^\S\n]*(\S+)[^\S\n]*$', text or "",
                  flags=re.M)
    if m:
        return m.group(1).strip("'\"")
    plugin = field(text, "plugin") or "?"
    return plugin.replace("!ref ", "").strip("'\"")


def identify(path):
    text = slurp(path, lambda n: n.endswith("/metadata.yaml")
                 and n.count("/") == 1)
    if not text or text.startswith("<<"):
        return None, None
    return field(text, "uuid"), field(text, "type")


def inputs_of(path):
    text = slurp(path, lambda n: n.endswith("provenance/action/action.yaml"))
    if not text or text.startswith("<<"):
        return "?", []
    return action_name(text), UUID.findall(text)


def walk():
    for root in ROOTS:
        if not os.path.isdir(root):
            continue
        for here, dirs, names in os.walk(root):
            dirs[:] = [d for d in dirs if d != "__pycache__"]
            for name in sorted(names):
                if name.endswith(".qza"):
                    yield os.path.join(here, name)


def main():
    if not os.path.isdir("results"):
        sys.exit("Run me from the project root; I see no results/ directory.")

    every = sorted(walk())
    print("%d QIIME 2 artifacts found" % len(every))

    seen = {}          # uuid -> (path, type)
    kinds = {}         # path -> type
    for path in every:
        uuid, kind = identify(path)
        kinds[path] = kind
        if uuid:
            seen[uuid] = (path, kind)
    print("%d of them carry a readable uuid" % len(seen))

    taxonomies = [p for p in every if "Taxonomy" in str(kinds.get(p))]
    print("%d are FeatureData[Taxonomy]" % len(taxonomies))

    for path in taxonomies:
        action, uuids = inputs_of(path)
        print()
        print("=" * 72)
        print(path)
        print("  made by    %s" % action)
        named = False
        for u in uuids:
            if u in seen and "Classifier" in str(seen[u][1]):
                print("  CLASSIFIER %s" % seen[u][0])
                print("             uuid %s" % u)
                named = True
        if not named:
            print("  no classifier among its recorded inputs matched a .qza")
            print("  on disk. Its first input uuids:")
            for u in uuids[:8]:
                print("      %s  %s" % (u, seen.get(u, ("not on disk",))[0]))

    print()
    print("=" * 72)
    print("classifiers on disk, for reference")
    for path in every:
        if "Classifier" in str(kinds.get(path)):
            uuid = next((u for u, v in seen.items() if v[0] == path), "?")
            print("  %-58s %s" % (path, uuid))

    print()
    print("Send this whole output back.")


if __name__ == "__main__":
    main()
