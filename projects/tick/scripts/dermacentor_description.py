#!/usr/bin/env python3
"""Dermacentor, described separately. No hypothesis test.

The carriage analysis was restricted to Ixodes because Dermacentor is a different vector
with a different endosymbiont and no panel pathogens at all. This is the other half: what
Dermacentor looks like, on its own terms.

**There is no carriage contrast to test.** The Dermacentor specimens are qPCR-negative on
every panel target, so there is no positive group. Nothing here is a hypothesis test and no
p-value is reported, and the group sizes would not support one even if a contrast existed.

Tick genus is read from results/analysis/tick_species_by_host.tsv, which build_species_figure.py
writes after applying the read floor, and library names resolve through scripts/tick_lookup.py.
This script used to do both itself and got both wrong; see patch_dermacentor_to_lookup.py.

What this does produce:
  - core microbiota of each genus, at three prevalence thresholds
  - the taxa that separate them, ranked by prevalence difference
  - the Francisella / Rickettsia reciprocity, quantified per specimen rather than as a
    group mean
  - the mixed-genus pools, which are neither and should not be silently averaged into one

Run:
    cd ~/meedlab/tick_project
    python3 scripts/dermacentor_description.py \
        --outdir results/analysis/dermacentor

Writes:
    <outdir>/core_microbiota_by_genus.tsv
    <outdir>/dermacentor_report.txt
"""

import argparse
import csv
import glob
import os
import re
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd

import tick_lookup

TABLE_QZA = "results/decontam/16S/table_decontam_16S.qza"
TAXONOMY = "results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv"
TICK_GENUS = "results/18S/all/taxonomy_verified/tick_genus_per_sample_corrected.tsv"
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
# The species figure's per-specimen genus calls, with the read floor already
# applied. One genus caller in the project, not two.
SPECIES_TABLE = "results/analysis/tick_species_by_host.tsv"
# The crosswalk holds 563 specimens. More keys than that means the key is not a
# specimen, which is how this script came to report 738 of them.
MAX_SPECIMENS = 563

# Non-bacterial sequence. Organellar, plus the 9,333 eukaryotic ASVs the 16S primers
# pick up: the tick (o__Acari), the deer ked (o__Diptera), the host blood meal
# (c__Mammalia) and 6,648 resolved no deeper than d__Eukaryota. All correctly
# labelled, none part of a bacterial community.
ORGANELLAR = ("f__Mitochondria", "f__Chloroplast", ";Mitochondria", ";Chloroplast",
              "d__Eukaryota", "d__Archaea")

# "Core" has no agreed definition, so I report three thresholds rather than picking one and
# letting the choice do the work invisibly.
CORE_THRESHOLDS = [0.50, 0.75, 0.90]

report_lines = []


def say(message=""):
    print(message)
    report_lines.append(str(message))


def resolve_specimen(lookup, library):
    """I turn a 16S library name into its specimen, through the shared lookup.

    The two helpers that used to live here did this themselves and did it
    wrong: one stripped hyphens so ALHT WTD_09 and alht-deer-09 became
    different keys, the other stripped a replicate index only when it was
    exactly two digits. Between them they made 738 keys out of 563 specimens.
    """
    text = str(library).strip()
    stem = text.split(".")[0]
    for marker, value in (("16S", text), ("16S", stem), ("specimen", stem)):
        hit = tick_lookup.resolve(lookup, marker, value)
        if hit is not None:
            return str(hit)
    return None


def genus_of(lineage):
    for part in reversed([p.strip() for p in str(lineage).split(";")]):
        name = part.split("__")[-1].strip()
        if name and name.lower() not in ("", "uncultured", "unassigned", "unclassified"):
            return name
    return "Unassigned"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="results/analysis/dermacentor")
    args = parser.parse_args()

    for path in (TABLE_QZA, TAXONOMY, SPECIES_TABLE):
        if not os.path.exists(path):
            sys.exit(f"no such file: {path}")
    os.makedirs(args.outdir, exist_ok=True)
    lookup = tick_lookup.load()

    # ---- tick genus per specimen, from the species figure -------------------
    species = pd.read_csv(SPECIES_TABLE, sep="\t", dtype=str)
    for column in ("specimen", "Ixodes_pos", "Dermacentor_pos"):
        if column not in species.columns:
            sys.exit(f"{SPECIES_TABLE} has no '{column}' column. Rebuild it:\n"
                     f"    python3 scripts/build_species_figure.py")

    def is_true(value):
        """I read the species table's booleans, whatever casing they carry."""
        return str(value).strip().lower() in ("true", "1", "yes", "t")

    genus_by_barcode = {}
    host_by_barcode = {}
    for _, row in species.iterrows():
        specimen = str(row["specimen"]).strip()
        if not specimen or specimen.lower() in ("nan", "none"):
            continue
        positives = [name for name, column in (("Ixodes", "Ixodes_pos"),
                                               ("Dermacentor", "Dermacentor_pos"),
                                               ("Amblyomma", "Amblyomma_pos"))
                     if column in species.columns and is_true(row.get(column))]
        if len(positives) == 1:
            genus_by_barcode[specimen] = positives[0]
        elif len(positives) > 1:
            genus_by_barcode[specimen] = "both/mixed"
        host_by_barcode[specimen] = str(row.get("host_bin") or "").strip()

    counts = pd.Series(genus_by_barcode).value_counts()
    say(f"tick genus from {SPECIES_TABLE}, read floor already applied")
    say(f"specimens with a tick-genus call: {len(genus_by_barcode)}")
    say(counts.to_string())

    if len(genus_by_barcode) > MAX_SPECIMENS:
        sys.exit(f"{len(genus_by_barcode)} specimens carry a genus call, which is "
                 f"more than the {MAX_SPECIMENS} in the crosswalk. The key is not "
                 f"a specimen. Do not trust anything below this line.")

    if counts.get("Dermacentor", 0) < 3:
        sys.exit("fewer than 3 Dermacentor specimens; nothing to describe")

    # ---- host of each, for context, from the same table ---------------------
    say("\nDermacentor specimens by collection host:")
    derm_hosts = pd.Series({b: host_by_barcode.get(b, "")
                            for b, g in genus_by_barcode.items() if g == "Dermacentor"})
    say(derm_hosts.replace("", "(not recorded)").value_counts().to_string())

    # ---- community ---------------------------------------------------------
    organellar = set()
    lineage_by_id = {}
    with open(TAXONOMY, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            if any(marker in row[1] for marker in ORGANELLAR):
                organellar.add(row[0].strip())
            else:
                lineage_by_id[row[0].strip()] = row[1].strip()

    import biom
    with tempfile.TemporaryDirectory() as temp_dir:
        with zipfile.ZipFile(TABLE_QZA) as archive:
            archive.extractall(temp_dir)
        table = biom.load_table(
            glob.glob(f"{temp_dir}/*/data/feature-table.biom")[0]).to_dataframe(dense=True)
    table = table.loc[[i for i in table.index if i not in organellar]]
    say(f"\ncommunity table: {table.shape[0]:,} ASVs x {table.shape[1]:,} libraries "
        f"({len(organellar):,} organellar ASVs removed)")

    summed = {}
    genus_of_specimen = {}
    for column in table.columns:
        key = resolve_specimen(lookup, column)
        if key is None or key not in genus_by_barcode:
            continue
        summed[key] = summed.get(key, 0) + table[column]
        genus_of_specimen[key] = genus_by_barcode[key]
    if not summed:
        sys.exit("no library matched a specimen with a tick-genus call")

    community = pd.DataFrame(summed).T
    depth = community.sum(axis=1)
    community = community[depth > 0]
    labels = pd.Series({k: genus_of_specimen[k] for k in community.index})

    say("\nspecimens with community data, by genus:")
    say(labels.value_counts().to_string())
    say("\nread depth by genus (organellar removed):")
    for genus in labels.unique():
        values = depth[labels[labels == genus].index]
        say(f"  {genus:<14} n={len(values):>3}  median={int(values.median()):>9,}")

    # ---- collapse to genus level -------------------------------------------
    relative = community.div(community.sum(axis=1), axis=0).fillna(0)
    genus_names = [genus_of(lineage_by_id.get(a, "")) for a in relative.columns]
    by_genus = relative.T.groupby(genus_names).sum().T

    # ---- core microbiota ----------------------------------------------------
    say("\n" + "=" * 70)
    say("CORE MICROBIOTA BY TICK GENUS")
    say("=" * 70)
    say("'Core' has no agreed definition, so three prevalence thresholds are reported")
    say("rather than one chosen threshold doing the work invisibly. A taxon is present in")
    say("a specimen if it has any reads at all.")

    core_rows = []
    for tick_genus in ("Ixodes", "Dermacentor"):
        members = labels[labels == tick_genus].index
        if len(members) < 3:
            continue
        subset = by_genus.loc[members]
        prevalence = (subset > 0).mean()
        say(f"\n  {tick_genus} (n = {len(members)}):")
        for threshold in CORE_THRESHOLDS:
            core = prevalence[prevalence >= threshold]
            say(f"    >={int(100*threshold)}% of specimens: {len(core)} taxa")
            if len(core) and threshold == 0.75:
                top = core.sort_values(ascending=False).head(12)
                for taxon, value in top.items():
                    mean_abundance = 100 * subset[taxon].mean()
                    say(f"        {taxon[:34]:<34} {100*value:5.0f}%  "
                        f"mean {mean_abundance:5.2f}%")
        for taxon, value in prevalence.items():
            core_rows.append({"tick_genus": tick_genus, "taxon": taxon,
                              "prevalence": round(float(value), 4),
                              "mean_rel_abundance_pct": round(100 * float(subset[taxon].mean()), 4),
                              "n_specimens": len(members)})

    # ---- what separates them ------------------------------------------------
    say("\n" + "=" * 70)
    say("TAXA THAT SEPARATE THE TWO GENERA (by prevalence difference)")
    say("=" * 70)
    ixodes = labels[labels == "Ixodes"].index
    dermacentor = labels[labels == "Dermacentor"].index
    say(f"Descriptive ranking, not a test. n = {len(dermacentor)} Dermacentor against "
        f"{len(ixodes)} Ixodes would not")
    say("support one, and the two groups differ in host, collection and sequencing depth.")
    if len(ixodes) >= 3 and len(dermacentor) >= 3:
        comparison = pd.DataFrame({
            "prev_Ixodes": (by_genus.loc[ixodes] > 0).mean(),
            "prev_Dermacentor": (by_genus.loc[dermacentor] > 0).mean(),
            "mean_pct_Ixodes": 100 * by_genus.loc[ixodes].mean(),
            "mean_pct_Dermacentor": 100 * by_genus.loc[dermacentor].mean(),
        })
        comparison["prev_diff"] = (comparison["prev_Dermacentor"]
                                   - comparison["prev_Ixodes"])
        top = comparison.reindex(comparison["prev_diff"].abs()
                                 .sort_values(ascending=False).index).head(15)
        say("")
        say(top.round(3).to_string())

    # ---- the reciprocity, per specimen -------------------------------------
    say("\n" + "=" * 70)
    say("RICKETTSIA / FRANCISELLA, PER SPECIMEN")
    say("=" * 70)
    say("Reported per specimen rather than as a group mean: a mean can be produced by a")
    say("few extreme libraries, and the claim is that this separates EVERY specimen.")
    for tick_genus in ("Ixodes", "Dermacentor", "both/mixed"):
        members = labels[labels == tick_genus].index
        if not len(members):
            continue
        rickettsia = (by_genus.loc[members, "Rickettsia"]
                      if "Rickettsia" in by_genus.columns else pd.Series(0.0, index=members))
        francisella = (by_genus.loc[members, "Francisella"]
                       if "Francisella" in by_genus.columns else pd.Series(0.0, index=members))
        say(f"\n  {tick_genus} (n = {len(members)}):")
        say(f"    Rickettsia  >0 in {int((rickettsia > 0).sum()):>3} of {len(members)}   "
            f"median {100*rickettsia.median():6.2f}%   max {100*rickettsia.max():6.2f}%")
        say(f"    Francisella >0 in {int((francisella > 0).sum()):>3} of {len(members)}   "
            f"median {100*francisella.median():6.2f}%   max {100*francisella.max():6.2f}%")

    mixed = labels[labels == "both/mixed"].index
    if len(mixed):
        say(f"\n  NOTE: {len(mixed)} specimen(s) have libraries calling BOTH genera. These")
        say("  are pooled samples containing both tick species and are neither group; they")
        say("  are shown separately above and excluded from the two-genus comparison.")

    if core_rows:
        out_path = f"{args.outdir}/core_microbiota_by_genus.tsv"
        # No .bak. The old backup step left one beside the output, and a .bak
        # from a previous run has been shared by accident once already. The
        # report below records what this run did.
        pd.DataFrame(core_rows).to_csv(out_path, sep="\t", index=False)
        say(f"\nwrote {out_path}")

    say("\nREAD BEFORE USING:")
    say(f"  - NO hypothesis test is run here and no p-value is reported. The "
        f"{len(dermacentor)} Dermacentor")
    say("    specimens were qPCR-negative on every panel target when last checked, so no")
    say("    carriage contrast exists within Dermacentor to test. [VERIFY that against the")
    say("    current qPCR join before quoting it; it was established on an older cohort.]")
    say("  - Dermacentor and Ixodes differ in collection host as well as species, so the")
    say("    two are confounded and the comparison cannot separate them.")
    say("  - depth differs between the groups; prevalence at a given threshold is not")
    say("    depth-corrected and a shallow library will show fewer taxa as 'present'.")

    report_path = f"{args.outdir}/dermacentor_report.txt"
    with open(report_path, "w") as handle:
        handle.write("\n".join(report_lines))
    print(f"report -> {report_path}")


if __name__ == "__main__":
    main()
