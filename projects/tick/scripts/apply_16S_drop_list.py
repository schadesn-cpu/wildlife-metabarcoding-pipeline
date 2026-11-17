#!/usr/bin/env python3
"""
apply_16S_drop_list.py  —  run ON RON, from the project root:
    cd ~/meedlab/tick_project
    python3 scripts/apply_16S_drop_list.py --report      # writes nothing
    python3 scripts/apply_16S_drop_list.py --apply

WHAT THIS IS
------------
The 16S pathogen-ASV drop list was decided 2026-08-04 from the per-ASV BLAST of all
223 pathogen-labelled ASVs, and never applied. It sits UPSTREAM of every number
locked since: single_tick_v3 inherits its calls from single_tick_v2, which was built
on the pre-drop-list taxonomy. So Anaplasma 18.7%, Borreliella 34.1% and the
2026-08-10 concordance are all still computed with these ASVs included.

NOTHING IS DELETED. Every change is a RELABEL. Wolbachia is a real organism and a
real finding (the deer ked, TICK_GROUND_TRUTH section 6); it is simply not a
pathogen, so those samples move out of the pathogen count and into an endosymbiont
row. Sub-threshold ASVs are relabelled to the rank they actually support.

REDUNDANCY, tested 2026-08-10 (this is what actually moves):
  3 sub-threshold Borreliella  ->   3 samples,   0 lost (fully redundant)
  9 Wolbachia-as-Anaplasma     -> 118 samples,  55 lost   <-- the real change
  2 sub-threshold Anaplasma    ->   2 samples,   1 lost
  1 Ca. Neoehrlichia           ->   1 sample,    1 lost

Borreliella is free. Anaplasma is not: 55 samples are Anaplasma-positive ONLY via a
Wolbachia ASV, against 267 positive via a clean ASV (~17% of positives).

Reads:
    results/16S/all/taxonomy_verified/taxonomy_16S_blast_verified.tsv
Writes:
    results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv
    results/16S/all/taxonomy_verified/droplist_changes.tsv
    results/16S/all/taxonomy_verified/droplist_report.txt
    (the previous taxonomy is copied to *_PRE_DROPLIST.tsv before anything is written)

AFTER APPLYING, RE-RUN IN ORDER:
    build_analysis_table.py -> build_single_tick_set.py -> build_concordance.py
    -> pathogen_cooccurrence.py
and point build_analysis_table.py's TAX at the new file.
"""
import argparse
import os
import shutil
import sys
import pandas as pd

TAXDIR = "results/16S/all/taxonomy_verified"
TAX    = f"{TAXDIR}/taxonomy_16S_blast_verified.tsv"
OUT    = f"{TAXDIR}/taxonomy_16S_droplist_applied.tsv"
CHG    = f"{TAXDIR}/droplist_changes.tsv"
BAK    = f"{TAXDIR}/taxonomy_16S_blast_verified_PRE_DROPLIST.tsv"

WOLBACHIA = ("d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rickettsiales;"
             "f__Anaplasmataceae;g__Wolbachia;s__")
ANAPLASMATACEAE = ("d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rickettsiales;"
                   "f__Anaplasmataceae;g__unassigned_Anaplasmataceae;s__")
BORRELIACEAE = ("d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Spirochaetales;"
                "f__Borreliaceae;g__unassigned_Borreliaceae;s__")
NEOEHRLICHIA = ("d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rickettsiales;"
                "f__Anaplasmataceae;g__Candidatus_Neoehrlichia;s__")

# prefix -> (new lineage, why). Prefixes are the 8-char forms used throughout the notes.
RELABEL = {
    # --- 9 Wolbachia mislabelled as Anaplasma. BLAST 98.9-100% to Wolbachia
    #     endosymbionts of hippoboscid flies (L. cervi, M. ovinus, Icosta sp.).
    #     SILVA collapsed them because Wolbachia and Anaplasma are both
    #     Anaplasmataceae. The deer ked is confirmed on two markers (section 6),
    #     which is what explains their presence. THE REAL CHANGE: 55 samples.
    "d89c872f": (WOLBACHIA, "Wolbachia 100%"),
    "0115599b": (WOLBACHIA, "Wolbachia 99.734%"),
    "f75f006b": (WOLBACHIA, "Wolbachia 99.733%"),
    "9656c248": (WOLBACHIA, "Wolbachia 99.733%"),
    "f95671c1": (WOLBACHIA, "Wolbachia 99.732%"),
    "017311f0": (WOLBACHIA, "Wolbachia 99.468%"),
    "511087b7": (WOLBACHIA, "Wolbachia 99.467%"),
    "91a688e0": (WOLBACHIA, "Wolbachia 98.936%"),
    "40d30de6": (WOLBACHIA, "Wolbachia 99.653% but 288bp - short, flagged"),

    # --- below the 94.5% genus threshold (Yarza et al. 2014): best match in ALL of
    #     nt is under the genus floor, so no genus is defensible. Family only.
    "bbd4b04e": (ANAPLASMATACEAE, "best nt match 91.6% - below genus floor"),
    "3175edc8": (ANAPLASMATACEAE, "95.8% uncultured alphaproteobacterium - not genus-callable"),
    "9c33bb21": (BORRELIACEAE,    "91.4% - below genus floor"),
    "ac679d98": (BORRELIACEAE,    "92.4% - below genus floor"),
    "278dfe41": (BORRELIACEAE,    "94.0% - below genus floor"),

    # --- NOT a drop. Off-panel pathogen genus, its own label, a Goal 3 observation.
    #     95.5% so genus-consistent but species is NOT callable. Do not fold into
    #     Anaplasma and do not delete.
    "865afed6": (NEOEHRLICHIA, "Ca. Neoehrlichia 95.5% - off-panel genus, genus level only"),
}

# NOT applied here. These need a RULE, not a list, and the rule is not settled.
DEFERRED = """
~17 Anaplasma ASVs sit at 91.6-96.3% identity, many best-matching "Uncultured
Anaplasma clone IS25": cdfea937 91.6 - 91434118 92.3 - ffae4480 92.6 - 26e9e5dd 93.0
9c85780a 93.1 - a4c98ffb 93.7 - 394abae9 93.9 - a1f614c2 93.9 - b7fea086 93.9
bb3bfc0a 94.2 - db26d1ba 94.3 - d46fedd7 94.5 - 1e4b9479 94.6 - 2c6babdf 95.0
6098e67a 95.1 - 69be8933 96.0 - a25b19f2 96.3
Decide ONE threshold rule for the whole set rather than adjudicating each. Settle the
confident removals above first and see where the numbers land.
"""

rep = []
def say(m=""):
    print(m); rep.append(str(m))

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
g = ap.add_mutually_exclusive_group(required=True)
g.add_argument("--report", action="store_true", help="show what would change, write nothing")
g.add_argument("--apply", action="store_true", help="write the corrected taxonomy")
args = ap.parse_args()

if not os.path.exists(TAX):
    sys.exit(f"ERROR: {TAX} not found")

tax = pd.read_csv(TAX, sep="\t", dtype=str)
idcol, taxcol = tax.columns[0], tax.columns[1]
say(f"taxonomy: {len(tax):,} features ({TAX})")

# resolve 8-char prefixes to full feature ids
say("\n=== resolving prefixes ===")
resolved, missing, ambiguous = {}, [], []
for pre, (new, why) in RELABEL.items():
    hits = tax[tax[idcol].astype(str).str.startswith(pre)]
    if len(hits) == 0:
        missing.append(pre)
    elif len(hits) > 1:
        ambiguous.append((pre, len(hits)))
    else:
        resolved[hits.iloc[0][idcol]] = (new, why, pre)
say(f"  resolved {len(resolved)} of {len(RELABEL)}")
if missing:
    say(f"  NOT FOUND: {missing}")
if ambiguous:
    say(f"  AMBIGUOUS PREFIX (would rewrite the wrong feature): {ambiguous}")
    sys.exit("ERROR: ambiguous prefixes. Nothing written.")
if missing:
    sys.exit("ERROR: some ASVs are not in this taxonomy file. Check you are applying "
             "to the right file. Nothing written.")

# what changes
rows = []
for fid, (new, why, pre) in resolved.items():
    old = tax.loc[tax[idcol] == fid, taxcol].iloc[0]
    og = old.split(";")[-2] if ";" in old else old
    ng = new.split(";")[-2]
    rows.append({"feature_id": fid, "prefix": pre, "old_genus": og,
                 "new_genus": ng, "reason": why, "old_taxon": old, "new_taxon": new})
chg = pd.DataFrame(rows)

say("\n=== changes ===")
for _, r in chg.iterrows():
    say(f"  {r['prefix']}  {r['old_genus']:<28} -> {r['new_genus']:<34} {r['reason']}")

say("\n=== summary by destination ===")
for ng, sub_ in chg.groupby("new_genus"):
    say(f"  {ng:<38} {len(sub_)} ASV(s)")

say("\nEXPECTED EFFECT (redundancy tested 2026-08-10):")
say("  Borreliella: 3 ASVs, 3 samples, 0 LOST  -> 34.1% unchanged")
say("  Anaplasma:  12 ASVs, ~119 samples, 56 LOST (55 of them from the Wolbachia group)")
say("              267 samples remain positive via a clean Anaplasma ASV")
say("  Wolbachia becomes a new ENDOSYMBIONT row carrying ~118 samples")
say("\n  The Anaplasma concordance will likely get WORSE, not better: sequencing")
say("  already calls 18.7% vs qPCR 50.1%, with 16 qPCR-only vs 2 metabarcoding-only.")
say("  Removing false positives widens that gap. Compute it; do not predict it.")

say("\nDEFERRED, not applied here:")
say(DEFERRED)

if args.report:
    say("\n" + "=" * 70)
    say("REPORT MODE. Nothing written. Re-run with --apply.")
    say("=" * 70)
    sys.exit(0)

# apply
if not os.path.exists(BAK):
    shutil.copy2(TAX, BAK)
    say(f"\nbacked up -> {BAK}")
else:
    say(f"\nbackup already exists, left alone -> {BAK}")

new_tax = tax.copy()
for fid, (new, why, pre) in resolved.items():
    new_tax.loc[new_tax[idcol] == fid, taxcol] = new
new_tax.to_csv(OUT, sep="\t", index=False)
chg.to_csv(CHG, sep="\t", index=False)
say(f"wrote {OUT}")
say(f"wrote {CHG}")

say("\n" + "=" * 70)
say("NEXT, IN THIS ORDER:")
say("  1. point build_analysis_table.py's TAX at:")
say(f"     {OUT}")
say("  2. python3 scripts/build_analysis_table.py")
say("  3. python3 scripts/build_single_tick_set.py")
say("  4. python3 scripts/build_concordance.py")
say("  5. python3 scripts/pathogen_cooccurrence.py")
say("  6. add a Wolbachia row to the endosymbiont block in build_single_tick_set.py")
say("     (ENDOSYMBIONTS list) so the 118 samples are reported, not silently lost")
say("=" * 70)

with open(f"{TAXDIR}/droplist_report.txt", "w") as fh:
    fh.write("\n".join(rep))
say(f"\nreport -> {TAXDIR}/droplist_report.txt")
