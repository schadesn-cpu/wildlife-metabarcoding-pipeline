#!/usr/bin/env python3
"""tick_lookup.py

I am the one place in this project where a name becomes a specimen.

WHY I EXIST
-----------
Fifteen scripts each carried their own name-parsing rule, and every one of
them was a slightly different reimplementation of the same three steps:
split on semicolon, casefold, strip the marker tag. Each reimplementation
lost specimens quietly rather than raising. build_concordance.py stripped
hyphens along with spaces, so 'ALHT WTD_09' and the qPCR sheet's
'alht-deer-09' became two different keys and six single ticks fell out of
the concordance cohort without a warning.

I do not parse names. build_4marker_crosswalk.py already resolved every
name in the workbook and now writes that resolution out as
results/metadata/library_lookup.tsv, one row per surface form. All I do is
load that table and look a name up in it, applying the one normalisation
the table's keys are stored in: lowercase, and whitespace or underscores
collapsed to a single hyphen. Every other variant, including the marker
tags, the run-id suffixes and the ALHT WTD / ALHT-Deer alias, is already a
row.

If a name is not in the table, I return None. That is the point. A lookup
miss is a visible absence, not a partial match. Do NOT add a fallback that
strips characters until something matches; that is the bug this replaces.

USE
---
    import tick_lookup
    table = tick_lookup.load()
    specimen = tick_lookup.resolve(table, "specimen", "alht-deer-09")
    specimen = tick_lookup.resolve(table, "16S", "ALHT-Deer-09-16s")

Markers: '16S', '18S', 'cytb', 'ITS' for library names, and 'specimen' for
barcodes, accessioning IDs and qPCR keys.
"""
import os
import re
import sys

import pandas as pd

LOOKUP_PATH = "results/metadata/library_lookup.tsv"
MARKERS = ("16S", "18S", "cytb", "ITS", "specimen")


def load(path=LOOKUP_PATH):
    """I load the lookup and refuse to run without it."""
    if not os.path.isfile(path):
        sys.exit("library lookup not found: " + path
                 + "\nBuild it first, from the project root:"
                 + "\n    python3 scripts/build_4marker_crosswalk.py")
    frame = pd.read_csv(path, sep="\t", dtype=str)
    for column in ("marker", "key", "specimen"):
        if column not in frame.columns:
            sys.exit("column '" + column + "' missing from " + path
                     + ". Rebuild the crosswalk rather than editing the "
                       "lookup by hand.")
    table = {}
    for marker, key, specimen in zip(frame["marker"], frame["key"],
                                     frame["specimen"]):
        table[(str(marker), str(key))] = str(specimen)
    if not table:
        sys.exit(path + " is empty. Rebuild the crosswalk.")
    return table


def normalise(value):
    """I apply the one normalisation the lookup's keys are stored in.

    Lowercase, and runs of whitespace or underscore collapsed to a single
    hyphen. Nothing else. The table holds every other variant, so a name
    that still does not match is genuinely absent.
    """
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in ("nan", "", "none"):
        return ""
    return re.sub(r"[\s_]+", "-", text).lower().strip("-")


def resolve(table, marker, value):
    """I return the specimen for one name, or None if it is not in the table."""
    if marker not in MARKERS:
        sys.exit("unknown marker '" + str(marker) + "'. Use one of: "
                 + ", ".join(MARKERS))
    key = normalise(value)
    if not key:
        return None
    return table.get((marker, key))


def resolve_field(table, marker, value):
    """I resolve a semicolon-delimited name field to the set of specimens in it.

    The crosswalk's name_* columns hold several library names per specimen,
    semicolon-delimited. Callers should never split these themselves.
    """
    found = set()
    for piece in str(value).split(";"):
        hit = resolve(table, marker, piece)
        if hit is not None:
            found.add(hit)
    return found

def load_batches(path=LOOKUP_PATH):
    """I load the (marker, key) -> run batch map from the same lookup file.

    Batch was added on 2026-09-18 so the pooled co-occurrence could be
    conditioned on it. A lookup built before that has no batch column, and I
    say so rather than returning an empty map that looks like "no batches".
    """
    if not os.path.isfile(path):
        sys.exit("library lookup not found: " + path
                 + "\nBuild it first, from the project root:"
                 + "\n    python3 scripts/build_4marker_crosswalk.py")
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if "batch" not in frame.columns:
        sys.exit(path + " has no 'batch' column. It was written before "
                 "2026-09-18. Rebuild the crosswalk:\n"
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
