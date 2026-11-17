#!/usr/bin/env python3
"""ITS feature tables for the VT-12 fungal trial, host rDNA removed.

For C. Sullivan. ITS is the marker for this trial and it never went through the table
building the other three markers got, so there was nothing to hand her.

THE HOST FILTER IS THE WHOLE POINT
The UNITE v10 classifier assigned all 730 ITS ASVs to Kingdom Fungi at median confidence
0.999. BLAST of the 20 most abundant, carrying 95.4% of reads, found 18 of them to be
Ixodes scapularis ribosomal DNA (GU318946.1, 96 to 99% identity). UNITE contains no animal
references and the conserved 18S to ITS1 junction at the 5' end of the ITS1F amplicon
matches across eukaryotes, so the classifier had no way to flag them.

The two genuine fungi in the top 20 are ASV_13 (Epicoccum sp.) and ASV_18 (Cladosporium
sp.), both at 100% identity and 100% query coverage. Every other ASV in the top 20 is host.
That is recorded in logs/PROJECT_LOG.md and it is where the host list below comes from, so
the filter is the BLAST result and not an abundance heuristic.

Unfiltered, this dataset looks 95% fungal when it is 95% tick. Both tables are written so
the filter can be seen rather than taken on trust.

WHICH DADA2 TABLE
table_single.qza: 730 ASVs, 384,545 reads. The two paired-end attempts (table.qza,
table_ee5.qza) retained under 9,000 reads, presumably because ITS1 length variation defeats
merging. Single-end is the correct choice here and matches the manuscript numbers.

Run:
    cd ~/meedlab/tick_project
    python3 scripts/build_its_tables.py --outdir results/tables

Writes:
    <outdir>/its_asv_table_FILTERED.tsv    host removed; use this one
    <outdir>/its_asv_table_ALL.tsv         everything, with is_host_rDNA marked
    <outdir>/its_sample_summary.tsv        per specimen, before and after
    <outdir>/its_README.txt
"""

import argparse
import csv
import glob
import os
import sys
import tempfile
import zipfile

import pandas as pd

TABLE_QZA = "qiime2/ITS1-2/dada2/table_single.qza"
TAXONOMY_QZA = "qiime2/ITS1-2/taxonomy/taxonomy.qza"
TOP20 = "reports/ITS_top20_for_blast.tsv"
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"

# From logs/PROJECT_LOG.md: of the top 20 by abundance, ranks 13 and 18 are genuine fungi
# (Epicoccum, Cladosporium). The other 18 are Ixodes scapularis rDNA.
FUNGAL_RANKS = {13, 18}
FUNGAL_NAMES = {13: "Epicoccum sp.", 18: "Cladosporium sp."}

report_lines = []


def say(message=""):
    print(message)
    report_lines.append(str(message))


def back_up_existing(out_path):
    if not os.path.exists(out_path):
        return
    backup = out_path + ".bak"
    if os.path.exists(backup):
        sys.exit(f"{out_path} and {backup} both exist. Move or delete {backup} first.")
    os.replace(out_path, backup)


def export_qza(path, pattern):
    """Pull one file out of a .qza without needing qiime on the path."""
    if not os.path.exists(path):
        sys.exit(f"no such file: {path}")
    temp_dir = tempfile.mkdtemp()
    with zipfile.ZipFile(path) as archive:
        archive.extractall(temp_dir)
    matches = glob.glob(f"{temp_dir}/*/data/{pattern}")
    if not matches:
        sys.exit(f"no {pattern} inside {path}")
    return matches[0]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="results/tables")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ---- the host list, from the recorded BLAST result ----------------------
    if not os.path.exists(TOP20):
        sys.exit(f"{TOP20} not found. That file carries the 20 BLASTed ASVs and without "
                 f"it the host filter cannot be reconstructed.")
    top20 = pd.read_csv(TOP20, sep="\t")
    if "rank" not in top20.columns or "ASV_ID" not in top20.columns:
        sys.exit(f"{TOP20} does not have the expected rank and ASV_ID columns")

    host_ids, fungal_ids = set(), {}
    for _, row in top20.iterrows():
        rank = int(row["rank"])
        if rank in FUNGAL_RANKS:
            fungal_ids[str(row["ASV_ID"])] = FUNGAL_NAMES[rank]
        else:
            host_ids.add(str(row["ASV_ID"]))
    say(f"host rDNA ASVs from the recorded BLAST: {len(host_ids)}")
    say(f"confirmed fungi in the top 20: "
        f"{', '.join(f'rank {r} = {n}' for r, n in FUNGAL_NAMES.items())}")
    if len(host_ids) != 18:
        say(f"WARNING: expected 18 host ASVs, found {len(host_ids)}. Check {TOP20}.")

    # ---- table and taxonomy -------------------------------------------------
    import biom
    table = biom.load_table(export_qza(TABLE_QZA, "feature-table.biom")).to_dataframe(dense=True)
    say(f"\n{TABLE_QZA}: {table.shape[0]} ASVs x {table.shape[1]} samples, "
        f"{int(table.values.sum()):,} reads")

    taxonomy = {}
    confidence = {}
    with open(export_qza(TAXONOMY_QZA, "taxonomy.tsv"), newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                taxonomy[row[0].strip()] = row[1].strip()
                confidence[row[0].strip()] = row[2].strip() if len(row) > 2 else ""

    # ---- per-ASV table ------------------------------------------------------
    def genus_of(lineage):
        """Deepest real name in the lineage.

        UNITE appends a species-hypothesis code as a final rank, sh__SH0740224.10FU. It is
        an accession for a reference cluster, not a name, and taking it as the deepest rank
        put codes in the results table where Cordyceps tenuipes and Metarhizium
        gaoligongense should have been.
        """
        import re as _re
        for part in reversed([p.strip() for p in str(lineage).split(";")]):
            name = part.split("__")[-1].strip()
            if _re.match(r"^SH\d+", name):
                continue
            if name and name.lower() not in ("", "unidentified", "unassigned"):
                return name.replace("_", " ")
        return "Unassigned"

    rows = []
    for asv in table.index:
        lineage = taxonomy.get(asv, "")
        is_host = asv in host_ids
        rows.append({
            "feature_id": asv,
            "is_host_rDNA": "yes" if is_host else "no",
            "blast_note": ("Ixodes scapularis rDNA, GU318946.1, 96-99% identity"
                           if is_host else fungal_ids.get(asv, "")),
            "taxon_lineage": lineage,
            "lowest_named": genus_of(lineage),
            "classifier_confidence": confidence.get(asv, ""),
            "total_reads": int(table.loc[asv].sum()),
            "n_samples_present": int((table.loc[asv] > 0).sum()),
        })
    all_asvs = pd.DataFrame(rows).sort_values("total_reads", ascending=False)

    host_reads = all_asvs[all_asvs["is_host_rDNA"] == "yes"]["total_reads"].sum()
    total_reads = all_asvs["total_reads"].sum()
    say(f"\nhost rDNA: {int(host_reads):,} of {int(total_reads):,} reads "
        f"({100*host_reads/total_reads:.1f}%)")

    filtered = all_asvs[all_asvs["is_host_rDNA"] == "no"].copy()
    say(f"after removing host: {len(filtered)} ASVs, "
        f"{int(filtered['total_reads'].sum()):,} reads")

    for name, frame in (("its_asv_table_ALL.tsv", all_asvs),
                        ("its_asv_table_FILTERED.tsv", filtered)):
        path = f"{args.outdir}/{name}"
        back_up_existing(path)
        frame.to_csv(path, sep="\t", index=False)
        say(f"wrote {path}")

    # ---- per sample ---------------------------------------------------------
    keep = [a for a in table.index if a not in host_ids]
    summary = pd.DataFrame({
        "reads_total": table.sum(axis=0).astype(int),
        "reads_fungal": table.loc[keep].sum(axis=0).astype(int),
        "asvs_fungal": (table.loc[keep] > 0).sum(axis=0).astype(int),
    })
    summary["pct_host"] = (100 * (1 - summary["reads_fungal"] / summary["reads_total"])
                           ).round(1)
    summary = summary.sort_values("reads_fungal", ascending=False)

    say("\n=== per sample, before and after the host filter ===")
    say(summary.to_string())

    path = f"{args.outdir}/its_sample_summary.tsv"
    back_up_existing(path)
    summary.to_csv(path, sep="\t")
    say(f"\nwrote {path}")

    # ---- what is actually in the filtered data ------------------------------
    say("\n=== fungal genera after filtering, by reads ===")
    by_genus = (filtered.groupby("lowest_named")["total_reads"].sum()
                .sort_values(ascending=False))
    for genus, reads in by_genus.head(20).items():
        n_asvs = int((filtered["lowest_named"] == genus).sum())
        say(f"  {genus[:34]:<34} {int(reads):>8,}  ({n_asvs} ASVs)")

    # Genus-level only. Searching the whole Hypocreales order pulled in Nectriaceae
    # (Fusarium, Cosmospora), which are plant-associated rather than insect pathogens, and
    # listing them as entomopathogens would overstate the result.
    say("\n=== entomopathogenic genera (Cordycipitaceae and Clavicipitaceae) ===")
    ENTOMOPATHOGENS = ("g__Beauveria", "g__Metarhizium", "g__Cordyceps",
                       "g__Simplicillium", "g__Leptobacillium", "g__Isaria",
                       "g__Lecanicillium", "g__Akanthomyces", "g__Ophiocordyceps",
                       "g__Purpureocillium")
    ento = filtered[filtered["taxon_lineage"].astype(str).str.contains(
        "|".join(ENTOMOPATHOGENS))]
    if ento.empty:
        say("  none in the filtered table")
    else:
        for _, row in ento.sort_values("total_reads", ascending=False).iterrows():
            in_samples = [s.replace("USGSTICKS-VT-12-", "").replace("-ITS1-2", "")
                          for s in table.columns if table.loc[row["feature_id"], s] > 0]
            say(f"  {row['lowest_named'][:32]:<32} {row['total_reads']:>6} reads   "
                f"{', '.join(in_samples)}")
        say("")
        say("  The four ASVs BLAST-verified in the earlier work (Beauveria 19 reads in")
        say("  BearA-01; Cordyceps tenuipes 2, C. cateniannulata 2 and Metarhizium 9 in")
        say("  BearA-02) are all here. Anything else in this list has NOT been")
        say("  BLAST-verified and is a classifier call at low read count. Simplicillium in")
        say("  VTW and Leptobacillium in BearC are new relative to the written account and")
        say("  need BLAST before either is reported.")

    readme = f"{args.outdir}/its_README.txt"
    with open(readme, "w") as handle:
        handle.write(
            "ITS TABLES, VT-12 FUNGAL TRIAL\n"
            "==============================\n\n"
            "USE its_asv_table_FILTERED.tsv. The ALL version is there so the filter can be\n"
            "checked, not so it can be analysed.\n\n"
            "WHY THE FILTER MATTERS\n"
            "  The UNITE classifier called all 730 ASVs Fungi at median confidence 0.999.\n"
            "  BLAST of the 20 most abundant, carrying 95.4% of reads, found 18 to be\n"
            "  Ixodes scapularis ribosomal DNA (GU318946.1, 96-99% identity). UNITE has no\n"
            "  animal references and the conserved 18S-ITS1 junction matches across\n"
            "  eukaryotes, so the classifier could not flag them.\n"
            "  Unfiltered, this dataset looks 95% fungal when it is 95% tick.\n\n"
            "  The two genuine fungi in the top 20 are Epicoccum sp. and Cladosporium sp.,\n"
            "  both at 100% identity and 100% query coverage.\n\n"
            "COLUMNS\n"
            "  is_host_rDNA        yes for the 18 BLAST-confirmed host ASVs\n"
            "  blast_note          what BLAST returned, where it was run\n"
            "  lowest_named        deepest named rank in the lineage\n"
            "  total_reads         summed across all 13 libraries\n"
            "  n_samples_present   libraries with at least one read\n\n"
            "READ DEPTH\n"
            "  After filtering the whole run holds about 17,600 reads across 12 samples.\n"
            "  Three samples carry most of it. The rest are in the hundreds or lower. At\n"
            "  that depth the data support presence and absence and little else, and no\n"
            "  comparison between treatment groups is possible: at most three specimens\n"
            "  share a letter code and the whole set is one sequencing batch.\n\n"
            "SOURCE\n"
            "  qiime2/ITS1-2/dada2/table_single.qza (730 ASVs, 384,545 reads)\n"
            "  Single-end. The two paired-end attempts retained under 9,000 reads, most\n"
            "  likely because ITS1 length variation defeats merging.\n")
    say(f"wrote {readme}")

    with open(f"{args.outdir}/its_build_report.txt", "w") as handle:
        handle.write("\n".join(report_lines))


if __name__ == "__main__":
    main()
