#!/usr/bin/env python3
"""patch_bartonella_out_of_coinfection.py

I take Bartonella out of the pathogen set in build_single_tick_set.py, so the
co-infection count and fig2 agree with how the manuscript already classifies
it.

WHY
---
Section 3.7 of the manuscript says these are ruminant-associated Bartonella,
not tick-borne human pathogens, and reports them "as a co-detection, outside
the pathogen results in Section 3.2". The co-infection count did the opposite
and counted Bartonella as a pathogen genus, and fig2 panel A carried a
Bartonella bar under the title "Pathogen prevalence in single ticks". The
paper contradicted itself in two places.

The evidence behind the classification, measured on the 93 single ticks on
2026-09-18: Bartonella is 10 of 93, deer ked is 11 of 93, and 6 of the 10
Bartonella-positive ticks also carry ked reads (odds ratio 23.4, Fisher
p = 0.0001). A seventh carries Wolbachia, the other ked endosymbiont. Eight of
the ten are USGSTICKS specimens from the NH 41xx/63xx/65xx series, the batch
cluster where keds were still being run before they were removed at
accessioning (Kellen Sawyer, 2026-09-12). The three with no ked signal are
explained in 3.7 by blood-meal acquisition from bacteremic deer, which is also
not a tick-borne human infection.

WHAT CHANGES
------------
Bartonella moves from PATHOGEN_GENUS to a new CO_DETECTIONS list, reported
alongside the endosymbionts. It therefore leaves the prevalence table, leaves
fig2 panel A, and leaves the co-infection denominatorless count.

    co-infection, distinct pathogen genera, n = 93
        before   0:38 (40.9%)  1:36 (38.7%)  2:14 (15.1%)  3:5 (5.4%)
                 >=1 55 (59.1%)   >=2 19 (20.4%)
        after    0:42 (45.2%)  1:37 (39.8%)  2:10 (10.8%)  3:4 (4.3%)
                 >=1 51 (54.8%)   >=2 14 (15.1%)

Those after-values are already confirmed twice: by
scripts/recompute_cooccurrence.py on ron on 2026-09-18, and independently in
the assistant's sandbox from the same table.

Bartonella's own prevalence, 10 of 93 (10.8%), is still printed. It is not
deleted, it is reclassified.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_bartonella_out_of_coinfection.py
    python3 scripts/build_single_tick_set.py
"""
import os
import sys

TARGET = "scripts/build_single_tick_set.py"

OLD_MAP = """    '16S_Anaplasma':          'Anaplasma',
    '16S_Bartonella':         'Bartonella',
    '16S_Borreliella':        'Borreliella',"""

NEW_MAP = """    '16S_Anaplasma':          'Anaplasma',
    '16S_Borreliella':        'Borreliella',"""

OLD_ENDO = """# Wolbachia is an endosymbiont here, not a pathogen: these ASVs are deer ked
# (Lipoptena cervi) symbionts that SILVA had filed under Anaplasma.
ENDOSYMBIONTS = ['16S_Rickettsia', '16S_Francisella', '16S_Wolbachia']"""

NEW_ENDO = """# Wolbachia is an endosymbiont here, not a pathogen: these ASVs are deer ked
# (Lipoptena cervi) symbionts that SILVA had filed under Anaplasma.
ENDOSYMBIONTS = ['16S_Rickettsia', '16S_Francisella', '16S_Wolbachia']

# Bartonella is a CO-DETECTION, not a pathogen, and is therefore not in
# PATHOGEN_GENUS above. Section 3.7 of the manuscript already says so: these
# are ruminant-associated Bartonella, not tick-borne human pathogens, reported
# outside the pathogen results. The co-infection count used to contradict that
# by counting them as a pathogen genus, which put 'this tick had a ked on it'
# into a co-infection statistic.
#
# Measured on these 93 single ticks, 2026-09-18: Bartonella 10/93, deer ked
# 11/93, and 6 of the 10 Bartonella-positive ticks carry ked reads (odds ratio
# 23.4, Fisher p = 0.0001). A seventh carries Wolbachia, the other ked
# endosymbiont. Keds were removed from tubes after the first processing batch
# (Kellen Sawyer, 2026-09-12), so the rate is partly a handling artifact.
#
# Reclassified, NOT deleted: the prevalence is still reported below.
CO_DETECTIONS = ['16S_Bartonella']"""

OLD_REPORT = '''say(f"\\n=== endosymbionts (context, NOT pathogens) ===")
for c in ENDOSYMBIONTS:
    if c in new.columns:
        pos = int((pd.to_numeric(new[c], errors='coerce').fillna(0) > 0).sum())
        say(f"  {c:<24} {pos:>3}/{n} = {100*pos/n:5.1f}%")'''

NEW_REPORT = '''say(f"\\n=== endosymbionts (context, NOT pathogens) ===")
for c in ENDOSYMBIONTS:
    if c in new.columns:
        pos = int((pd.to_numeric(new[c], errors='coerce').fillna(0) > 0).sum())
        say(f"  {c:<24} {pos:>3}/{n} = {100*pos/n:5.1f}%")

say(f"\\n=== co-detections (context, NOT pathogens, NOT in co-infection) ===")
for c in CO_DETECTIONS:
    if c in new.columns:
        pos = int((pd.to_numeric(new[c], errors='coerce').fillna(0) > 0).sum())
        say(f"  {c:<24} {pos:>3}/{n} = {100*pos/n:5.1f}%")
        if '16S_Bartonella' == c and 'cytb_Deer_ked' in new.columns:
            bart = pd.to_numeric(new[c], errors='coerce').fillna(0) > 0
            ked = pd.to_numeric(new['cytb_Deer_ked'], errors='coerce').fillna(0) > 0
            both = int((bart & ked).sum())
            say(f"      of which {both} also carry deer ked reads "
                f"({100*both/max(pos,1):.0f}% of Bartonella positives); "
                f"ked is {int(ked.sum())}/{n}")
say("  Reported in Results 3.7 as a co-detection, not in 3.2 or 3.3.")'''

OLD_COINF = 'say(f"\\n=== co-infection, distinct pathogen GENERA (n = {n}) ===")'
NEW_COINF = ('say(f"\\n=== co-infection, distinct pathogen GENERA (n = {n}) ===")\n'
             'say("  Bartonella EXCLUDED as a co-detection (see CO_DETECTIONS above).")')

OLD_NOTE = ("NOTE ON CO-INFECTION: counted at GENUS level. Babesia odocoilei + B. microti in one\n"
            "tick is ONE genus, not two. This is what the published fig2 does, and it is why the\n"
            "old set gave 22 ticks with >=2 genera rather than 23 with >=2 species.")

NEW_NOTE = ("NOTE ON CO-INFECTION: counted at GENUS level. Babesia odocoilei + B. microti in one\n"
            "tick is ONE genus, not two. This is what the published fig2 does, and it is why the\n"
            "old set gave 22 ticks with >=2 genera rather than 23 with >=2 species.\n"
            "\n"
            "Bartonella is NOT counted, as of 2026-09-18. See CO_DETECTIONS below for the\n"
            "reasoning and the numbers. Counting it gave >=1 55/93 (59.1%) and >=2 19/93\n"
            "(20.4%); excluding it gives >=1 51/93 (54.8%) and >=2 14/93 (15.1%). If you are\n"
            "reading a document that quotes 59.1% or 20.4%, it predates this change.")

REPLACEMENTS = [
    (OLD_NOTE, NEW_NOTE, "the co-infection note in the docstring"),
    (OLD_MAP, NEW_MAP, "the Bartonella entry in PATHOGEN_GENUS"),
    (OLD_ENDO, NEW_ENDO, "the ENDOSYMBIONTS block"),
    (OLD_REPORT, NEW_REPORT, "the endosymbiont report block"),
    (OLD_COINF, NEW_COINF, "the co-infection heading"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "CO_DETECTIONS" in source:
        sys.exit("already patched: " + TARGET + " defines CO_DETECTIONS. "
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
    print("expect co-infection 42 / 37 / 10 / 4, with 51 (54.8%) at >=1 "
          "and 14 (15.1%) at >=2")
    print("fig2 panel A loses its Bartonella bar and keeps six")


if __name__ == "__main__":
    main()
