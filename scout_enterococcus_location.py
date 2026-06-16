"""Quick scout: per-bird Enterococcus % + location + COD, sorted."""
import pandas as pd
import sys
import re

RELABUND = "results/16S/all/taxonomy_refined/taxonomy_relabund_L6_16S.tsv"
METADATA = "metadata/qiime/metadata_16S_codrelabel.tsv"

relabund = pd.read_csv(RELABUND, sep="\t", index_col=0)
ent_rows = [idx for idx in relabund.index if "Enterococcus" in idx]
print(f"Found {len(ent_rows)} Enterococcus row(s): {ent_rows}", file=sys.stderr)

ent_per_sample = relabund.loc[ent_rows].sum(axis=0) * 100  # to percent

# Strip the _S#### batch suffix so IDs match metadata
def strip_suffix(s):
    return re.sub(r"_S\d+$", "", s)

df = pd.DataFrame({
    "sample_id_raw": ent_per_sample.index,
    "sample_id": [strip_suffix(s) for s in ent_per_sample.index],
    "enterococcus_pct": ent_per_sample.values,
})

meta = pd.read_csv(METADATA, sep="\t", dtype=str)
meta = meta[meta["sample-id"] != "#q2:types"]

df = df.merge(meta[["sample-id", "TV", "COD_broad", "Season", "State Found", "Location", "Date Found"]],
              left_on="sample_id", right_on="sample-id", how="left")
df = df.drop(columns=["sample-id", "sample_id_raw", "sample_id"])
df = df.sort_values("enterococcus_pct", ascending=False).reset_index(drop=True)

pd.set_option("display.max_rows", None)
pd.set_option("display.max_colwidth", 40)
pd.set_option("display.width", 200)
print(df[["TV", "enterococcus_pct", "COD_broad", "Season", "State Found", "Location", "Date Found"]]
      .to_string(index=False, float_format=lambda x: f"{x:6.1f}"))

# Quick sanity check — how many merged successfully?
matched = df["TV"].notna().sum()
total = len(df)
print(f"\n[Merged {matched}/{total} samples to metadata]", file=sys.stderr)

