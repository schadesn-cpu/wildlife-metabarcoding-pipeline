#!/usr/bin/env python3
"""
Build the per-sample tick-genus read table from the 18S ASV-level data.

CORRECTED SOURCE (2026-06-30): the L7 relabund table collapses ticks into
Arachnida/Chelicerata and cannot be used. This builds from the per-ASV
taxonomy joined to the exported PLAIN 18S feature table (read counts).

18S resolves Ixodes vs Dermacentor at GENUS only. Within-Ixodes species
labels (affinis/ricinus/scapularis) are classifier artifacts and are NOT
used here beyond folding into genus Ixodes. cookei is tracked in its own
column for the Powassan follow-up but the real cookei call is cytB.

Inputs:
  results/18S/all/taxonomy/taxonomy_summary_18S.tsv   (Feature ID | Taxon | Confidence)
  /tmp/18s_feature_table.tsv                          (biom convert --to-tsv of table.qza)

  If /tmp cleared, re-export first:
    qiime tools export --input-path qiime2/18S/all/dada2/table.qza --output-path /tmp/18s_table_export
    biom convert -i /tmp/18s_table_export/feature-table.biom -o /tmp/18s_feature_table.tsv --to-tsv

Output:
  results/18S/all/tick_genus_per_sample.tsv   (durable)

Usage:
  python scripts/build_tick_genus_per_sample.py
"""
import sys
import pandas as pd

TAX = "results/18S/all/taxonomy/taxonomy_summary_18S.tsv"
TABLE = "/tmp/18s_feature_table.tsv"
OUT = "results/18S/all/tick_genus_per_sample.tsv"


def tick_genus(taxon):
    """Return Ixodes / Dermacentor / cookei / None for a lineage string."""
    t = str(taxon)
    if "Dermacentor" in t:
        return "Dermacentor"
    if "Ixodes_cookei" in t:
        return "cookei"
    if "Ixodes" in t:
        return "Ixodes"
    return None


def main():
    # --- taxonomy: map each tick ASV to a genus bucket ---
    tax = pd.read_csv(TAX, sep="\t")
    fid_col, taxon_col = tax.columns[0], tax.columns[1]
    tax["_g"] = tax[taxon_col].map(tick_genus)
    tick_map = tax.loc[tax["_g"].notna(), [fid_col, "_g"]].set_index(fid_col)["_g"]
    print("tick ASVs in taxonomy: %d" % len(tick_map))

    # --- feature table: skip the '# Constructed from biom file' comment line ---
    tbl = pd.read_csv(TABLE, sep="\t", skiprows=1)
    tbl = tbl.rename(columns={tbl.columns[0]: "ASV"}).set_index("ASV")

    # keep only tick ASVs present in the table
    shared = [a for a in tick_map.index if a in tbl.index]
    print("tick ASVs found in feature table: %d of %d" % (len(shared), len(tick_map)))
    if not shared:
        sys.exit("No tick ASVs overlap the feature table; check ID formats.")

    sub = tbl.loc[shared].copy()
    sub["_g"] = tick_map.loc[shared].values

    # sum reads per sample per genus bucket -> samples x {Ixodes,Dermacentor,cookei}
    per = sub.groupby("_g").sum().T
    for col in ("Ixodes", "Dermacentor", "cookei"):
        if col not in per.columns:
            per[col] = 0.0

    per["tick_reads"] = per[["Ixodes", "Dermacentor", "cookei"]].sum(axis=1)

    def dominant(r):
        ix = r["Ixodes"] + r["cookei"]   # cookei is within-Ixodes for the deer/moose split
        dm = r["Dermacentor"]
        if ix == 0 and dm == 0:
            return "none"
        return "Ixodes" if ix >= dm else "Dermacentor"

    per["dominant"] = per.apply(dominant, axis=1)
    out = (per.reset_index()
              .rename(columns={"index": "sample_id"})
              [["sample_id", "Ixodes", "Dermacentor", "cookei", "tick_reads", "dominant"]])
    out.to_csv(OUT, sep="\t", index=False)

    print("\nwrote %s" % OUT)
    print("%d samples, %d with tick reads" % (len(out), (out["tick_reads"] > 0).sum()))
    print(out["dominant"].value_counts().to_string())
    print("\ntop XB blanks by tick_reads (control vetting):")
    xb = out[out["sample_id"].str.contains("XB", case=False)].sort_values("tick_reads", ascending=False)
    print(xb.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
