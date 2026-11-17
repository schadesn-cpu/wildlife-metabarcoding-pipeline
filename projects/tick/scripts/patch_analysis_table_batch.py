#!/usr/bin/env python3
"""patch_analysis_table_batch.py

I add the 16S sequencing run batch to analysis_table_16S.tsv, and print its
distribution, so the pooled co-occurrence can be conditioned on it.

WHY
---
Anaplasma x Borreliella at pool level is OR 2.03, q = 0.0002. It survives pool
size (Mantel-Haenszel 2.05) and read depth (1.89), and it has been stable
across the anchoring fix, the join fix and two crosswalk rebuilds: 2.48 crude
and 2.27 depth-adjusted in August on 558 samples against 2.03 and 1.89 today
on 602.

Batch is the confounder it has not been tested against, and it is the largest
effect in this dataset: R2 10.9%, p = 0.001, with 530 of 551 animals entirely
within one batch.

A CAUTION ABOUT WHAT 'BATCH' MEANS HERE
---------------------------------------
The run sheet's Run Batch strings are PLATE-AND-MARKER, not the analysis
batch. There are 40 distinct values, and AMP-021925-VDL-Tick-P7-16s, -18s and
-cytoB are one plate under three names. Some carry no marker suffix at all
(AMP-120924-USGS-TICKS-P3-5).

Filtering to 16S narrows it, but the ~10 batches used in the ked table
(section 1.12) and in the R2 10.9% diversity result are a coarser grouping
than these strings. So this patch DOES NOT stratify anything. It adds the
column and prints how many samples fall in each value, so the grouping can be
looked at before any odds ratio is computed on it. Stratifying on 20 strings
with five samples each would produce noise that looks like an answer.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_analysis_table_batch.py
    python3 scripts/build_analysis_table.py
"""
import os
import sys

TARGET = "scripts/build_analysis_table.py"

OLD_RESOLVE = '''    def resolve(sample_id, base_id):
        """I turn a feature-table sample id into the crosswalk record, or None."""
        for marker, value in (('16S', sample_id), ('16S', base_id),
                              ('specimen', base_id)):
            hit = tick_lookup.resolve(lookup, marker, value)
            if hit is not None:
                return specimen_meta.get(str(hit))
        return None
'''

NEW_RESOLVE = '''    batches = tick_lookup.load_batches()

    def resolve(sample_id, base_id):
        """I turn a feature-table sample id into the crosswalk record, or None."""
        for marker, value in (('16S', sample_id), ('16S', base_id),
                              ('specimen', base_id)):
            hit = tick_lookup.resolve(lookup, marker, value)
            if hit is not None:
                return specimen_meta.get(str(hit))
        return None

    def resolve_batch(sample_id, base_id):
        """I return the 16S run batch for a sample, or None.

        Batch is a LIBRARY property, so I only ask the 16S namespace. A
        specimen key carries no batch by construction, and falling back to it
        would silently attach one library's plate to another's reads.
        """
        for value in (sample_id, base_id):
            hit = tick_lookup.resolve_batch(batches, '16S', value)
            if hit is not None:
                return hit
        return None
'''

OLD_REC = """        rec={'sample':s,'Barcode':meta['Barcode'],'Host':meta['Host'],
             'Host_clean':meta.get('Host_clean'),'Host_group':meta['Host_group'],
             'NumTicks':meta['NumTicks'],'total_reads':int(total)}"""

NEW_REC = """        rec={'sample':s,'Barcode':meta['Barcode'],'Host':meta['Host'],
             'Host_clean':meta.get('Host_clean'),'Host_group':meta['Host_group'],
             'NumTicks':meta['NumTicks'],'total_reads':int(total),
             'RunBatch':resolve_batch(s, base)}"""

OLD_TAIL = '''    print("\\n=== co-detections (NOT pathogens; see ground truth 1.16) ==="'''

NEW_TAIL = '''    print("\\n=== 16S run batch, for conditioning (NOT yet stratified on) ===")
    if 'RunBatch' in df.columns:
        counts = df['RunBatch'].fillna('(none)').value_counts()
        print(f"  {len(counts)} distinct value(s) across {len(df)} samples")
        for name, n in counts.items():
            print(f"    {str(name):<44} {n:>4}")
        print("  CAUTION: these are PLATE-AND-MARKER strings from the run sheet,")
        print("  not the ~10 analysis batches used in the ked table (ground truth")
        print("  1.12) or in the R2 10.9% diversity result. Look at the grouping")
        print("  before stratifying anything on it; strata of five samples produce")
        print("  noise that looks like an answer.")

    print("\\n=== co-detections (NOT pathogens; see ground truth 1.16) ==="'''

REPLACEMENTS = [
    (OLD_RESOLVE, NEW_RESOLVE, "the resolve() definition"),
    (OLD_REC, NEW_REC, "the per-sample record"),
    (OLD_TAIL, NEW_TAIL, "the co-detections report heading"),
]


def main():
    if not os.path.isfile(TARGET):
        sys.exit("not found: " + TARGET
                 + "\nRun this from ~/meedlab/tick_project.")
    with open(TARGET, encoding="utf-8") as handle:
        source = handle.read()

    if "RunBatch" in source:
        sys.exit("already patched: " + TARGET + " writes RunBatch. "
                 "Nothing changed.")
    if "import tick_lookup" not in source:
        sys.exit(TARGET + " is not on the lookup yet. Run "
                 "patch_analysis_table_to_lookup.py first.")

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
    print()
    print("602 samples and every pathogen count must not move. This adds a")
    print("column and a listing, nothing else.")
    print()
    print("run: python3 " + TARGET)


if __name__ == "__main__":
    main()
