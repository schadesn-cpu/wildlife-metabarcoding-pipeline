#!/usr/bin/env python3
"""Add Homo to the cytb host map in build_single_tick_v2.py.

Human was never in the map. Lines 101-112 list ten wildlife hosts and no comment
marks the omission as deliberate, so human cytb has been dropped from the
single-tick table since it was written, and the EXPECTED dict on line 135 has ten
keys too, which meant nothing ever flagged it as missing.

Human is in 162 specimens across 167,318 cytb reads, so this is not a trace signal.
I add it with the caveat attached rather than as just another host, because human
cytb CANNOT distinguish a blood meal from handler or laboratory contamination.

I do NOT add an EXPECTED count for it: the July checkpoint gives 8 candidates at
fraction >= 10% and reads >= 1000, but the single-tick subset has changed since and
I would be asserting a number I have not recomputed. A missing EXPECTED entry is a
blank, not a failure.

Run on ron from the project root:
    python3 patch_add_human_host.py
"""
import sys

PATH = "scripts/build_single_tick_v2.py"
source = open(PATH).read()
original_length = len(source)
applied = []


def sub(old, new, label):
    global source
    if source.count(old) != 1:
        sys.exit(f"{label}: found {source.count(old)} matches, expected 1. "
                 f"Nothing has been changed.")
    source = source.replace(old, new)
    applied.append(label)


sub('    "cytb_Raccoon":            ["Procyon"],',
    '    "cytb_Raccoon":            ["Procyon"],\n'
    '    # Human. NOT a confirmed blood meal: human cytb cannot distinguish feeding\n'
    '    # from handler or laboratory contamination. Human is 0.41% of all cytb reads\n'
    '    # and 252 of 837 biological samples carry trace human, so a detection here\n'
    '    # means "human-positive above baseline". Reported, never called a meal.\n'
    '    # See Human_Checkpoint_Notes, 2026-07-09.\n'
    '    "cytb_Human":              ["Homo"],',
    "Homo in the host map")

sub('              "cytb_Opossum","cytb_Raccoon"]',
    '              "cytb_Opossum","cytb_Raccoon","cytb_Human"]',
    "cytb_Human in OUT_ORDER")

open(PATH, "w").write(source)
print(f"{original_length} -> {len(source)}")
for entry in applied:
    print("  +", entry)
print()
print("Rerun, in this order:")
print("  python3 scripts/build_single_tick_v2.py     # regenerates v2 with cytb_Human")
print("  python3 scripts/build_single_tick_set.py    # v3 inherits the new column")
print()
print("Then check it landed:")
print("  head -1 results/single_tick_v3/single_tick_v3_full.tsv | tr '\\\\t' '\\\\n' | grep -i human")
print("  awk -F'\\\\t' 'NR>1 && $NF>0' results/single_tick_v3/single_tick_v3_full.tsv | wc -l")
