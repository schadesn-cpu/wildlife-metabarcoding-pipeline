#!/usr/bin/env python3
"""One ASV table per marker, for sharing.

Everything anyone needs to check a call, in one place per marker: the sequence, the
classifier's assignment and confidence, what BLAST said where we checked it, the
category from the decision table, how many reads and how many samples. Someone who
disagrees with a call can see the evidence and disagree with it.

I write one file per marker rather than one combined table. 16S alone has ~90,000
ASVs, and the read floors and verification status differ enough between markers that
a combined file would need a column explaining which rules applied to which rows.

I do NOT filter. Low-abundance ASVs, unassigned ASVs and cross-amplified ASVs are all
in the file with a column saying what they are. A filtered table is a table where I
already decided what someone gets to see.

Run:
  cd ~/meedlab/tick_project
  python3 scripts/make_asv_tables.py --outdir results/tables

Writes:
  results/tables/asv_table_16S.tsv
  results/tables/asv_table_18S.tsv
  results/tables/asv_table_cytb.tsv
  results/tables/asv_table_README.txt
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile

TAXONOMY_PATHS = {
    "16S": ["results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv",
            "results/16S/all/taxonomy_verified/taxonomy_16S_blast_verified.tsv"],
    "18S": ["results/18S/all/taxonomy_verified/taxonomy_18S_blast_verified.tsv"],
    # cytb: the FULL verified taxonomy, not taxonomy_clean.tsv. clean.tsv is the
    # 3,709 ASVs that survived filtering; the verified file has all 52,671 with an
    # Action column saying why each was excluded. Sharing only the survivors would
    # hide 93% of the marker and the reasons for every exclusion.
    "cytb": ["results/cytb/all/taxonomy_verified/taxonomy_blast_verified.tsv",
             "results/cytb/all/taxonomy_verified/taxonomy_clean.tsv"],
}

REPSEQ_PATHS = {
    "16S": ["qiime2/16S/all/rep-seqs-exported-temp/dna-sequences.fasta",
            "qiime2/16S/all/dada2/rep-seqs.qza"],
    "18S": ["qiime2/18S/all/rep-seqs-exported-temp/dna-sequences.fasta",
            "qiime2/18S/all/dada2/rep-seqs.qza"],
    "cytb": ["qiime2/cytb/all/rep-seqs-exported-temp/dna-sequences.fasta",
             "qiime2/cytb/all/dada2/rep-seqs.qza"],
}

FEATURE_TABLE_PATHS = {
    "16S": ["results/cooccurrence/16S/biom_export/feature-table.tsv"],
    "18S": ["results/cooccurrence/18S/biom_export/feature-table.tsv"],
    "cytb": ["results/cooccurrence/cytb/biom_export/feature-table.tsv"],
}

# Per-ASV BLAST results, where we ran them. Only a subset of ASVs were checked, and
# the file says which, so a blank verification column means "not checked" and not
# "failed".
BLAST_PATHS = {
    "16S": [("results/16S/all/blast_verify/pathogen_verification_report.tsv",
             "pathogen ASVs, 2026-08-03"),
            ("results/16S/all/blast_verify/audit_fp_Rickettsia.tsv",
             "Rickettsia false-positive audit, 2026-05-01"),
            ("results/16S/all/blast_verify/audit_fp_Anaplasma.tsv",
             "Anaplasma false-positive audit, 2026-05-01")],
    "18S": [("results/18S/all/blast_verify/piroplasm_asvs_blast.tsv",
             "piroplasm ASVs, 2026-08-10")],
    "cytb": [],
}

DECISION_TABLE = "results/tables/taxon_decision_table.tsv"


def first_existing(paths):
    return next((p for p in paths if os.path.exists(p)), None)


def back_up_existing(out_path):
    """Moves an existing output aside to .bak. Refuses to clobber a backup."""
    if not os.path.exists(out_path):
        return
    backup_path = out_path + ".bak"
    if os.path.exists(backup_path):
        sys.exit(f"{out_path} and {backup_path} both exist. "
                 f"Move or delete {backup_path} before rerunning.")
    os.replace(out_path, backup_path)


def load_sequences(paths, marker):
    """feature_id -> sequence. Exports the .qza if only that is available."""
    path = first_existing(paths)
    if path is None:
        print(f"warning: no rep-seqs for {marker}; sequence column will be blank",
              file=sys.stderr)
        return {}
    if path.endswith(".qza"):
        temp_dir = tempfile.mkdtemp()
        subprocess.run(["qiime", "tools", "export", "--input-path", path,
                        "--output-path", temp_dir], check=True,
                       stdout=subprocess.DEVNULL)
        path = os.path.join(temp_dir, "dna-sequences.fasta")
    sequences = {}
    feature_id = None
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                feature_id = line[1:].split()[0]
                sequences[feature_id] = ""
            elif feature_id:
                sequences[feature_id] += line.strip()
    return sequences


def load_abundance(paths, marker):
    """feature_id -> (total_reads, n_libraries_present).

    Read straight off the biom export rather than through pandas: the 16S table is
    330 MB and only two summary numbers per row are needed.
    """
    path = first_existing(paths)
    if path is None:
        print(f"warning: no feature table for {marker}; read counts will be blank",
              file=sys.stderr)
        return {}
    abundance = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            values = [float(value) for value in fields[1:] if value]
            abundance[fields[0]] = (int(sum(values)),
                                    sum(1 for value in values if value > 0))
    return abundance


def load_blast(entries, marker):
    """feature_id -> (best_hit, percent_identity, which_search).

    Different BLAST outputs have different columns, so I look for the fields by name
    where there is a header and fall back to position for raw blastn output.
    """
    results = {}
    for path, label in entries:
        if not os.path.exists(path):
            continue
        with open(path, newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header is None:
                continue
            index = {name: position for position, name in enumerate(header)}
            hit_column = next((index[name] for name in
                               ("BLAST_top_hit", "actual_taxon") if name in index),
                              None)
            identity_column = index.get("pident")
            for row in reader:
                if not row or not row[0].strip():
                    continue
                feature_id = row[0].strip()
                if feature_id in results:
                    continue          # first search wins; they are ordered
                if hit_column is not None and len(row) > hit_column:
                    hit = row[hit_column]
                    identity = (row[identity_column]
                                if identity_column is not None
                                and len(row) > identity_column else "")
                elif len(row) >= 7:   # raw blastn -outfmt 6 with stitle last
                    hit, identity = row[-1], row[1]
                else:
                    continue
                results[feature_id] = (hit.strip(), identity.strip(), label)
    return results


def lookup_by_taxon_name(decisions_by_name, name):
    """Find a decision by taxon NAME rather than by full lineage.

    Only used for remapped cytb ASVs, where the lineage is stale but the verified
    taxon name is right. I match the whole name first, then the genus alone, so
    'Lipoptena cervi' finds a decision filed under 'Lipoptena_cervi' or 'Lipoptena'.
    I do NOT fall back to a substring search: that is what put 67 'uncultured' taxa
    into endosymbiont in the decision table, and I am not reintroducing it here.
    """
    key = name.strip().lower().replace(" ", "_")
    if key in decisions_by_name:
        return decisions_by_name[key]
    genus = key.split("_")[0]
    return decisions_by_name.get(genus)


def normalise_lineage(lineage):
    """A lineage key that ignores trailing empty ranks.

    'f__Rickettsiaceae;g__;s__' and 'f__Rickettsiaceae' are the same assignment
    written two ways, and an exact string join drops the category for one of them
    without saying so. I strip empty trailing ranks from both sides instead.
    """
    parts = [part.strip() for part in str(lineage).split(";")]
    while parts and (not parts[-1] or parts[-1].endswith("__")):
        parts.pop()
    return ";".join(parts)


def load_decisions():
    """full_lineage -> (category, reason, resolution_limit, rank_resolved)."""
    if not os.path.exists(DECISION_TABLE):
        print(f"warning: {DECISION_TABLE} not found. Run "
              f"make_taxon_decision_table.py first so the category column is "
              f"populated.", file=sys.stderr)
        return {}
    decisions = {}
    with open(DECISION_TABLE, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            value = (row.get("category", ""), row.get("reason", ""),
                     row.get("resolution_limit", ""), row.get("rank_resolved", ""))
            decisions[row["full_lineage"]] = value
            decisions.setdefault(normalise_lineage(row["full_lineage"]), value)
    return decisions


def load_taxonomy(paths, marker):
    """Returns (path, [(feature_id, lineage, confidence, verified, action)]).

    The cytb verified taxonomy carries two extra columns, Taxon_verified and
    Action, recording what the Entrez remap decided about each ASV and why it was
    kept or excluded. I carry those through rather than dropping them: 'EXCLUDED,
    unresolved_excluded' is information someone reading this table wants.
    """
    path = first_existing(paths)
    if path is None:
        return None, []
    rows = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None) or []
        index = {name.strip(): position for position, name in enumerate(header)}
        verified_column = index.get("Taxon_verified")
        action_column = index.get("Action")
        for row in reader:
            if len(row) < 2 or not row[0].strip():
                continue
            rows.append((
                row[0].strip(),
                row[1].strip(),
                row[2].strip() if len(row) > 2 else "",
                row[verified_column].strip()
                if verified_column is not None and len(row) > verified_column else "",
                row[action_column].strip()
                if action_column is not None and len(row) > action_column else ""))
    return path, rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="results/tables")
    parser.add_argument("--markers", nargs="*", default=list(TAXONOMY_PATHS))
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    decisions = load_decisions()
    # A second index keyed on the taxon NAME, for the remapped-cytb case above.
    decisions_by_name = {}
    for lineage, value in decisions.items():
        leaf = lineage.split(";")[-1].strip() or lineage.split(";")[-2].strip()
        leaf = leaf.split("__")[-1].strip().lower().replace(" ", "_")
        if leaf and leaf not in decisions_by_name:
            decisions_by_name[leaf] = value
    summary = []

    for marker in args.markers:
        print(f"\n=== {marker} ===")
        taxonomy_path, taxonomy_rows = load_taxonomy(TAXONOMY_PATHS[marker], marker)
        if taxonomy_path is None:
            print(f"  no taxonomy found; skipping")
            continue
        print(f"  taxonomy   {taxonomy_path}  ({len(taxonomy_rows):,} ASVs)")

        sequences = load_sequences(REPSEQ_PATHS[marker], marker)
        abundance = load_abundance(FEATURE_TABLE_PATHS[marker], marker)
        blast = load_blast(BLAST_PATHS[marker], marker)
        print(f"  sequences  {len(sequences):,}")
        print(f"  abundance  {len(abundance):,}")
        print(f"  BLAST      {len(blast):,} ASVs checked")

        out_rows = []
        for feature_id, lineage, confidence, verified, action in taxonomy_rows:
            # 0 not blank: an ASV absent from the feature table genuinely has zero
            # reads. A blank would read as "not measured", which is different.
            total_reads, n_libraries = abundance.get(feature_id, (0, 0))
            best_hit, identity, search = blast.get(feature_id, ("", "", ""))
            # A remapped cytb ASV still carries its ORIGINAL lineage in the taxon
            # column, because the remap corrects identity without rewriting the
            # lineage string. The 121 deer-ked ASVs are still filed under
            # Amphibia;Caudata;Plethodontidae. Joining the decision table on that
            # lineage matched only the 691 kept_as_is ASVs and missed all 3,018
            # remapped ones, so the ked was invisible to anyone filtering the shared
            # table by category. Where a verified taxon exists, I categorise on that.
            category, reason, limit, rank = decisions.get(
                lineage, decisions.get(normalise_lineage(lineage), ("", "", "", "")))
            if verified and verified.upper() != "EXCLUDED":
                remapped = lookup_by_taxon_name(decisions_by_name, verified)
                if remapped:
                    category, reason, limit, rank = remapped
            out_rows.append({
                "feature_id": feature_id,
                "marker": marker,
                "taxon_lineage": lineage,
                "rank_resolved": rank,
                "classifier_confidence": confidence,
                "category": category,
                "category_reason": reason,
                "resolution_limit": limit,
                "blast_best_hit": best_hit,
                "blast_percent_identity": identity,
                "blast_search": search,
                # cytb only: what the Entrez accession-level remap decided.
                "remap_verified_taxon": verified,
                "remap_action": action,
                "total_reads": total_reads,
                "n_libraries_present": n_libraries,
                "sequence": sequences.get(feature_id, ""),
            })

        # sort on reads, treating a missing count as zero rather than as an error:
        # an ASV in the taxonomy but not in the feature table has no reads by
        # definition, and I want it at the bottom rather than crashing the sort.
        out_rows.sort(key=lambda row: -(row["total_reads"]
                                        if isinstance(row["total_reads"], int) else 0))
        out_path = os.path.join(args.outdir, f"asv_table_{marker}.tsv")
        back_up_existing(out_path)
        with open(out_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t",
                                    fieldnames=list(out_rows[0]),
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(out_rows)
        print(f"  wrote {out_path}  ({len(out_rows):,} rows)")

        n_categorised = sum(1 for row in out_rows if row["category"])
        n_blasted = sum(1 for row in out_rows if row["blast_best_hit"])
        summary.append((marker, len(out_rows), n_categorised, n_blasted))

    readme = os.path.join(args.outdir, "asv_table_README.txt")
    with open(readme, "w") as handle:
        handle.write(
            "ASV TABLES — one file per marker\n"
            "================================\n\n"
            "One row per ASV. Nothing is filtered out: low-abundance, unassigned and\n"
            "cross-amplified ASVs are all here, with a column saying what they are.\n\n"
            "COLUMNS\n"
            "  feature_id              the ASV hash, stable across all project files\n"
            "  taxon_lineage           full classifier lineage, post-correction\n"
            "  rank_resolved           deepest rank actually assigned\n"
            "  classifier_confidence   the classifier's own confidence\n"
            "  category                from taxon_decision_table.tsv\n"
            "  category_reason         why that call was made\n"
            "  resolution_limit        where we have decided a species call is NOT\n"
            "                          possible, and why\n"
            "  blast_best_hit          BLAST result where we checked this ASV;\n"
            "                          BLANK MEANS NOT CHECKED, not failed\n"
            "  blast_percent_identity  identity for that hit\n"
            "  blast_search            which search it came from, with the date\n"
            "  total_reads             summed across all libraries\n"
            "  n_libraries_present     libraries with at least one read\n"
            "  sequence                the representative sequence\n\n"
            "IMPORTANT\n"
            "  - A blank category is not an error. Most of the community is neither a\n"
            "    pathogen nor an endosymbiont and needs no call.\n"
            "  - A blank blast_best_hit means that ASV was not part of a BLAST search.\n"
            "    Only pathogen-labelled and audited ASVs were checked.\n"
            "  - cytb shows 0 BLAST-checked, but cytb WAS verified: by accession-level\n"
            "    Entrez remap rather than by a blastn run, so there is no per-ASV BLAST\n"
            "    file to join. See results/cytb/all/blast_verify/cytb_remap.tsv, and the\n"
            "    remap_verified_taxon and remap_action columns in the cytb table.\n"
            "  - The cytb table includes ASVs the remap EXCLUDED, with the reason in\n"
            "    remap_action. They are here on purpose: 3,709 of 52,671 cytb ASVs were\n"
            "    retained, and the reasons for the rest are part of the record.\n"
            "  - Identity floors, Yarza et al. 2014: 94.5% for a genus call, 98.7% for\n"
            "    a species call.\n"
            "  - Read counts come from the co-occurrence biom export, which is\n"
            "    PRE-DECONTAM. 51 ASVs flagged by decontam are still present here.\n"
            "    They are marked reagent_contaminant_flagged in the category column.\n\n"
            "SOURCE\n")
        for marker, rows, categorised, blasted in summary:
            handle.write(f"  {marker}: {rows:,} ASVs, {categorised:,} categorised, "
                         f"{blasted:,} BLAST-checked\n")
    print(f"\nwrote {readme}")


if __name__ == "__main__":
    main()
