#!/usr/bin/env python3
"""
picrust2_nsti_check.py  — the PiCRUSt2 go/no-go for the tick 16S community.

Why abundance-weighted, not median:
Your 16S is ~93% Rickettsia endosymbiont. The per-ASV MEDIAN NSTI (the number
the loon run reported as 0.078) can look fine while the community is actually
untrustworthy, because it treats a rare well-referenced ASV and the dominant
endosymbiont equally. The number that matters here is the ABUNDANCE-WEIGHTED
NSTI, and the NSTI of the dominant taxon specifically. If the endosymbiont that
carries most of your reads has no close sequenced relative, the whole prediction
is that endosymbiont's guessed genome, and PiCRUSt2 is not defensible.

Run from tick_project root, AFTER phase 1 of run_picrust2_tick_16S.sh.
"""
import sys
import pandas as pd

OUT = "results/picrust2"
NSTI_FILE = f"{OUT}/marker_nsti_predicted.tsv.gz"
BIOM_FILE = f"{OUT}/inputs/table/feature-table.biom"
TAX_FILE  = "qiime2/16S/all/taxonomy-exported-temp/taxonomy.tsv"  # for naming top ASVs

# --- load per-ASV NSTI ------------------------------------------------------
nsti = pd.read_csv(NSTI_FILE, sep="\t")
asv_col = nsti.columns[0]
nsti_col = next(c for c in nsti.columns if "NSTI" in c.upper())
nsti = nsti.rename(columns={asv_col: "ASV"}).set_index("ASV")
nsti.index = nsti.index.astype(str)

# --- per-ASV total abundance from the feature table ------------------------
import biom
t = biom.load_table(BIOM_FILE).to_dataframe(dense=True)   # ASV x sample
abund = t.sum(axis=1)
abund.index = abund.index.astype(str)

df = nsti.copy()
df["abund"] = abund.reindex(df.index).fillna(0)
df = df[df["abund"] > 0]                                    # ASVs actually in the table
if df.empty:
    sys.exit("No overlap between NSTI ASVs and table ASVs — check inputs/paths.")

# --- the numbers that decide it --------------------------------------------
total = df["abund"].sum()
wmean = (df[nsti_col] * df["abund"]).sum() / total          # abundance-weighted NSTI
median = df[nsti_col].median()                              # per-ASV median (loon-style)
kept = df[df[nsti_col] <= 2.0]                              # PiCRUSt default cutoff
pct_reads_kept = 100 * kept["abund"].sum() / total
pct_asvs_kept = 100 * len(kept) / len(df)

print(f"ASVs in table            : {len(df):,}   total reads: {int(total):,}")
print(f"per-ASV MEDIAN NSTI       : {median:.3f}    (loon reference was 0.078)")
print(f"ABUNDANCE-WEIGHTED NSTI   : {wmean:.3f}    <-- the number that matters")
print(f"ASVs retained (NSTI<=2)   : {pct_asvs_kept:.1f}%")
print(f"READS retained (NSTI<=2)  : {pct_reads_kept:.1f}%")

# --- who drives the prediction, and can we trust them? ---------------------
try:
    tax = pd.read_csv(TAX_FILE, sep="\t").set_index("Feature ID")["Taxon"]
    df["taxon"] = tax.reindex(df.index)
except Exception:
    df["taxon"] = "(taxonomy not loaded)"

print("\nTop 8 ASVs by abundance (these dominate any function prediction):")
for a, r in df.sort_values("abund", ascending=False).head(8).iterrows():
    frac = 100 * r["abund"] / total
    flag = "  <-- DROPPED (NSTI>2)" if r[nsti_col] > 2 else ""
    print(f"  {frac:5.1f}% reads   NSTI={r[nsti_col]:.3f}{flag}   {str(r['taxon'])[:58]}")

# --- verdict ----------------------------------------------------------------
print("\n----------------------------- verdict -----------------------------")
if wmean > 1.0 or pct_reads_kept < 80:
    print("NO-GO. Abundance-weighted NSTI is high and/or a large read fraction is")
    print("dropped — the taxa carrying your reads have no close sequenced relatives.")
    print("Function predictions would be unreliable. Report ONE sentence: PiCRUSt2")
    print("was not applied because the endosymbiont-dominated community yields high")
    print("NSTI. That is a legitimate, defensible result — not a gap.")
elif wmean > 0.5:
    print("MARGINAL. Usable only as heavily-caveated SUPPLEMENTARY functional")
    print("potential. If you report it, state the abundance-weighted NSTI in the")
    print("caption and frame everything as coarse potential, not activity.")
    print("Consider re-running on a table with Rickettsia removed to describe the")
    print("REST of the community, and report both NSTI values.")
else:
    print("GO. NSTI is acceptable. Set RUN_PHASE2=1 in run_picrust2_tick_16S.sh")
    print("and re-run to produce KO / pathway tables. Still report the")
    print("abundance-weighted NSTI in the methods.")
