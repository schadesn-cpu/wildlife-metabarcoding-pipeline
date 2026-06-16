#!/usr/bin/env python3
"""verify_mifish_pa_cohort.py — reproduce the n=14 PA cohort derivation."""
import subprocess
import tempfile
import os
import pandas as pd
from biom import load_table

QZA = "qiime2/MiFish/all/dada2/table_filtered_final.qza"
META = "metadata/qiime/metadata_MiFish.tsv"   # was: ..._cod_filtered_clean.tsv
READ_THRESHOLD = 1000   # was 500
TAXON_THRESHOLD = 0.01  # ≥1%
PRIMARY_GROUPS = ["Diseased", "Trauma", "Marine"]

# 1. Extract biom from qza
with tempfile.TemporaryDirectory() as tmp:
    biom_path = os.path.join(tmp, "table.biom")
    out = subprocess.run(
        ["unzip", "-p", QZA, "*/data/feature-table.biom"],
        capture_output=True, check=True
    )
    with open(biom_path, "wb") as f:
        f.write(out.stdout)
    table = load_table(biom_path)

df = table.to_dataframe(dense=True)  # rows=ASVs, cols=samples
print(f"\n=== STEP 1: Raw cleaned table ===")
print(f"  Total entries: {df.shape[1]} samples × {df.shape[0]} ASVs")
print(f"  Total reads: {int(df.values.sum()):,}")

# 2. Drop controls
biological = [s for s in df.columns if s.startswith("TV")]
controls = [s for s in df.columns if not s.startswith("TV")]
df_bio = df[biological]
print(f"\n=== STEP 2: Drop controls ===")
print(f"  Biological samples: {len(biological)}")
print(f"  Controls dropped: {controls}")

# 3. Apply ≥500 reads/sample
totals = df_bio.sum(axis=0)
pass_reads = totals[totals >= READ_THRESHOLD].index.tolist()
df_reads = df_bio[pass_reads]
print(f"\n=== STEP 3: ≥{READ_THRESHOLD} reads/sample ===")
print(f"  Samples passing: {len(pass_reads)}")

# 4. Apply ≥1% per-taxon-per-sample → PA conversion
relabund = df_reads.div(df_reads.sum(axis=0), axis=1)
pa = (relabund >= TAXON_THRESHOLD).astype(int)
print(f"\n=== STEP 4: ≥{TAXON_THRESHOLD*100:.0f}% per-taxon threshold ===")
print(f"  Total taxon-sample detections: {int(pa.values.sum())}")
print(f"  Distinct taxa with ≥1 detection: {(pa.sum(axis=1) > 0).sum()}")

# 5. Drop samples with zero detections after threshold
det_per_sample = pa.sum(axis=0)
with_detections = det_per_sample[det_per_sample > 0].index.tolist()
print(f"\n=== STEP 5: Drop samples with zero detections ===")
print(f"  Samples retained: {len(with_detections)}")
zero_det = [s for s in pass_reads if s not in with_detections]
if zero_det:
    print(f"  Dropped (zero detections): {zero_det}")

# 6. Load metadata and restrict to primary groups
meta = pd.read_csv(META, sep="\t", comment="#", dtype=str)
meta.columns = [c.strip() for c in meta.columns]
print(f"\n=== STEP 6: Group restriction ({', '.join(PRIMARY_GROUPS)}) ===")
print(f"  Metadata entries: {len(meta)}")
print(f"  Group breakdown in metadata: {dict(meta['Group'].value_counts())}")

# Match PA samples (with -GI suffix) to metadata
pa_in_meta = meta[meta['sample-id'].isin(with_detections)]
final_cohort = pa_in_meta[pa_in_meta['Group'].isin(PRIMARY_GROUPS)]
print(f"\n=== FINAL COHORT ===")
print(f"  Total: {len(final_cohort)}")
print(f"  Breakdown: {dict(final_cohort['Group'].value_counts())}")

# Show which samples
print(f"\n  Final cohort samples by group:")
for grp in PRIMARY_GROUPS:
    samples_in_grp = final_cohort[final_cohort['Group'] == grp]['sample-id'].tolist()
    print(f"    {grp} (n={len(samples_in_grp)}): {samples_in_grp}")

# Show what was filtered out at the group step
in_pa_not_in_final = [s for s in with_detections if s not in final_cohort['sample-id'].tolist()]
if in_pa_not_in_final:
    print(f"\n  Passed PA thresholds but excluded by group:")
    for s in in_pa_not_in_final:
        row = meta[meta['sample-id'] == s]
        grp = row['Group'].values[0] if len(row) else "NOT IN METADATA"
        print(f"    {s}: {grp}")

# Show samples in pa_in_meta vs meta total
in_meta_not_pa = [s for s in meta['sample-id'].tolist() if s not in with_detections]
if in_meta_not_pa and len(in_meta_not_pa) <= 10:
    print(f"\n  In metadata but didn't pass PA thresholds:")
    for s in in_meta_not_pa:
        row = meta[meta['sample-id'] == s]
        grp = row['Group'].values[0] if len(row) else "?"
        print(f"    {s}: {grp}")
