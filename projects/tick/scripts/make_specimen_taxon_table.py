#!/usr/bin/env python3
"""What was found in each animal: one row per specimen per taxon detected.

This is the table someone opens first. The ASV tables answer "what is this
sequence"; this one answers "what was in this deer", which is the question a
collaborator actually has. Filter on host, on state, on category, on a taxon name,
and read across.

One row per specimen x marker x taxon, with the collection metadata attached and
the category from the decision table. Detections only: a row exists because
something was found, so this is long and sparse rather than a wide matrix of mostly
zeros.

I do NOT apply a read threshold here. Every detection is a row, with its read count,
and the reader decides. A table where I have already dropped everything below some
cutoff is a table where I have made that decision on someone else's behalf and not
told them. The exception is that I mark rather than remove: is_single_tick and
category are columns, so the pathogen-prevalence subset is a filter away.

Run:
  cd ~/meedlab/tick_project
  python3 scripts/make_specimen_taxon_table.py --out results/tables/specimen_taxon_table.tsv

Writes:
  results/tables/specimen_taxon_table.tsv
  results/tables/specimen_taxon_README.txt
"""

import argparse
import csv
import re
import os
import sys
from collections import defaultdict

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
DECISION_TABLE = "results/tables/taxon_decision_table.tsv"

ASV_TABLES = {
    "16S": "results/tables/asv_table_16S.tsv",
    "18S": "results/tables/asv_table_18S.tsv",
    "cytb": "results/tables/asv_table_cytb.tsv",
}

FEATURE_TABLES = {
    "16S": "results/cooccurrence/16S/biom_export/feature-table.tsv",
    "18S": "results/cooccurrence/18S/biom_export/feature-table.tsv",
    "cytb": "results/cooccurrence/cytb/biom_export/feature-table.tsv",
}

# Which crosswalk column holds the library names for each marker.
LIBRARY_COLUMNS = {"16S": "name_16S", "18S": "name_18S", "cytb": "name_cytb"}

# Suffixes a library name may or may not carry in the crosswalk.
# Every spelling a marker tag takes in the crosswalk. cytb is written both
# "-cytoB" (696 names) and "-cb" (90 names); knowing only the first dropped
# 49 libraries, 9 of them single ticks, out of the blood-meal analysis.
# The FIRST entry is the canonical form appended to a bare name.
MARKER_SUFFIXES = {"16S": ("-16s",), "18S": ("-18s",), "cytb": ("-cytoB", "-cb")}

# A marker tag is written several ways: "-18s", "-18sV4" through "-18sV51"
# (variable region), "-cytoB" or "-cb", and sometimes with a replicate number
# after the tag ("-CytoB-10"). Matching only the plain spelling dropped 49 cytb
# and 314 18S libraries, including the six single ticks behind the Anaplasma
# and Babesia discrepancies. I strip whatever the pattern matches and register
# the stem under every canonical spelling.
MARKER_TAG_PATTERN = re.compile(r"-(16s|18s|cytob|cb)([A-Za-z]*\d*)(-\d+)?$", re.I)


def back_up_existing(out_path):
    """Moves an existing output aside to .bak. Refuses to clobber a backup."""
    if not os.path.exists(out_path):
        return
    backup_path = out_path + ".bak"
    if os.path.exists(backup_path):
        sys.exit(f"{out_path} and {backup_path} both exist. "
                 f"Move or delete {backup_path} before rerunning.")
    os.replace(out_path, backup_path)


def normalise_library_name(name):
    """I casefold a library name before it is used as a join key.

    16 specimens carry lowercase library names in the crosswalk while the
    sequencing outputs are uppercase. A case-sensitive lookup dropped 15 of
    them from this table, including 2 single ticks, which is why the table
    reported 91 single ticks against the crosswalk's 93.
    """
    return name.strip().lower()


def load_crosswalk():
    """Returns (specimen_metadata, library_to_specimen_per_marker).

    name_16S and friends are SEMICOLON-DELIMITED lists, because one specimen can
    have several libraries, and the names sometimes carry a marker suffix and
    sometimes do not. Splitting on ';' and registering both forms is not optional:
    joining on the raw field returns about 13% of the true matches, which cost us
    two wrong answers before it was noticed.
    """
    if not os.path.exists(CROSSWALK):
        sys.exit(f"{CROSSWALK} not found. Run build_4marker_crosswalk.py first.")

    specimens = {}
    library_to_specimen = {marker: {} for marker in LIBRARY_COLUMNS}
    with open(CROSSWALK, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            specimen_id = row["sample_id"]
            specimens[specimen_id] = row
            for marker, column in LIBRARY_COLUMNS.items():
                for name in str(row.get(column, "")).split(";"):
                    name = name.strip()
                    if not name or name.lower() == "nan":
                        continue
                    suffixes = MARKER_SUFFIXES[marker]
                    forms = {name}
                    tag = MARKER_TAG_PATTERN.search(name)
                    stem = name[: tag.start()] if tag else name
                    forms.add(stem)
                    for suffix in suffixes:
                        forms.add(stem + suffix)
                    for form in forms:
                        key = normalise_library_name(form)
                        previous = library_to_specimen[marker].get(key)
                        if previous is not None and previous != specimen_id:
                            sys.exit("library '" + form + "' maps to two "
                                     "specimens once casefolded: " + previous
                                     + " and " + specimen_id)
                        library_to_specimen[marker][key] = specimen_id
    return specimens, library_to_specimen


def load_asv_metadata(path):
    """feature_id -> its row from the ASV table, for taxon and category."""
    if not os.path.exists(path):
        return {}
    metadata = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metadata[row["feature_id"]] = row
    return metadata


def taxon_label(asv_row):
    """The name to show for an ASV.

    For a remapped cytb ASV the lineage is the ORIGINAL bad call and the corrected
    identity is in remap_verified_taxon, so I show the verified name where there is
    one. Showing the lineage would file 121 deer-ked detections under
    Plethodontidae, which is a salamander family.
    """
    verified = (asv_row.get("remap_verified_taxon") or "").strip()
    if verified and verified.upper() != "EXCLUDED":
        return verified

    # A cytb ASV the remap REJECTED must not be shown under its original lineage. The
    # remap looked at these and threw them out; displaying the discarded label makes a
    # reject look like a result. 1,635 fish ASVs (Actinopteri, Gobiidae) sat in this
    # table as tick blood meals for exactly this reason, and I read them as
    # contamination before checking remap_action.
    action = (asv_row.get("remap_action") or "").strip()
    if action in ("unresolved_excluded", "excluded"):
        lineage_leaf = asv_row.get("taxon_lineage", "").split(";")[-1].strip()
        lineage_leaf = lineage_leaf.split("__")[-1].strip() or "unresolved"
        return f"[REJECTED by remap] {lineage_leaf}"
    lineage = asv_row.get("taxon_lineage", "")
    parts = [part.strip() for part in lineage.split(";") if part.strip()]
    for part in reversed(parts):
        name = part.split("__")[-1].strip()
        if name and name.lower() not in ("", "unassigned", "unclassified", "uncultured"):
            return name
    return "Unassigned"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="results/tables/specimen_taxon_table.tsv")
    parser.add_argument("--markers", nargs="*", default=list(ASV_TABLES))
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    specimens, library_to_specimen = load_crosswalk()
    print(f"crosswalk: {len(specimens):,} specimens")

    out_rows = []
    unmatched_libraries = defaultdict(set)

    for marker in args.markers:
        asv_path = ASV_TABLES.get(marker)
        feature_path = FEATURE_TABLES.get(marker)
        if not asv_path or not os.path.exists(asv_path):
            print(f"  {marker}: no ASV table, skipping", file=sys.stderr)
            continue
        if not feature_path or not os.path.exists(feature_path):
            print(f"  {marker}: no feature table, skipping", file=sys.stderr)
            continue

        asv_metadata = load_asv_metadata(asv_path)
        marker_map = library_to_specimen[marker]

        # (specimen, taxon) -> [reads, n_libraries, category, ...]
        totals = {}
        header = None
        with open(feature_path) as handle:
            for line in handle:
                if line.startswith("# "):
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields[0] in ("#OTU ID", "OTU ID"):
                    header = fields
                    continue
                if header is None or fields[0] not in asv_metadata:
                    continue
                asv_row = asv_metadata[fields[0]]
                label = taxon_label(asv_row)
                for position in range(1, len(fields)):
                    if not fields[position]:
                        continue
                    reads = float(fields[position])
                    if reads <= 0:
                        continue
                    library = header[position].split(".")[0]
                    specimen_id = marker_map.get(normalise_library_name(library))
                    if specimen_id is None:
                        unmatched_libraries[marker].add(library)
                        continue
                    # I group on (specimen, label, CATEGORY), not on the label
                    # alone. Two different lineages can produce the same display
                    # label with different categories: 'Rickettsiaceae' is both a
                    # family-level Rickettsiaceae assignment and the leaf of any
                    # lineage whose genus and species ranks are empty. Grouping on
                    # the label let whichever ASV arrived first set the category for
                    # all of them, which blanked a category that was correct in the
                    # ASV table. Same failure as every other label-keyed join in
                    # this project.
                    category = asv_row.get("category", "")
                    key = (specimen_id, label, category)
                    if key not in totals:
                        totals[key] = {
                            "reads": 0, "libraries": set(), "asvs": set(),
                            "category": category,
                            "resolution_limit": asv_row.get("resolution_limit", ""),
                            "rank": asv_row.get("rank_resolved", ""),
                        }
                    totals[key]["reads"] += reads
                    totals[key]["libraries"].add(library)
                    totals[key]["asvs"].add(fields[0])

        for (specimen_id, label, _category), value in totals.items():
            specimen = specimens[specimen_id]
            out_rows.append({
                "specimen": specimen_id,
                "barcode": specimen.get("barcode", ""),
                "host": specimen.get("host_bin", ""),
                "host_species": specimen.get("host_spp", ""),
                "state": specimen.get("state", ""),
                "n_ticks_in_sample": specimen.get("num_ticks", ""),
                "is_single_tick": specimen.get("is_single_tick", ""),
                "marker": marker,
                "taxon": label,
                "category": value["category"],
                "rank_resolved": value["rank"],
                "resolution_limit": value["resolution_limit"],
                "reads": int(value["reads"]),
                "n_asvs": len(value["asvs"]),
                "n_libraries": len(value["libraries"]),
            })
        print(f"  {marker}: {len(totals):,} specimen x taxon detections")

    if not out_rows:
        sys.exit("no detections assembled; check the ASV and feature table paths")

    out_rows.sort(key=lambda row: (row["specimen"], row["marker"], -row["reads"]))
    back_up_existing(args.out)
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t",
                                fieldnames=list(out_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"\nwrote {args.out}  ({len(out_rows):,} rows)")

    # I list the unmatched libraries rather than only counting them. The old
    # message asserted they were controls and excluded specimens, which was a
    # hardcoded string and was never checked. It hid 15 real specimens that a
    # case-sensitive join was dropping, so the count now comes with names.
    # guessed_specimen is a GUESS, produced by stripping the run suffix and any
    # trailing marker tag. It is for triage only. Nothing reads this file.
    unmatched_path = os.path.join(os.path.dirname(args.out) or ".",
                                  "unmatched_libraries.tsv")
    with open(unmatched_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["marker", "library", "guessed_specimen",
                         "guess_in_crosswalk"])
        for marker, libraries in sorted(unmatched_libraries.items()):
            for library in sorted(libraries):
                guess = str(library).split(".")[0]
                lowered = guess.lower()
                for tag in ("-16s", "-18s", "-cytob", "-cb"):
                    position = lowered.rfind(tag)
                    if position > 0:
                        guess = guess[:position]
                        break
                in_cw = guess in specimens
                writer.writerow([marker, library, guess, in_cw])

    for marker, libraries in sorted(unmatched_libraries.items()):
        print(f"  note: {len(libraries)} {marker} libraries had no crosswalk row",
              file=sys.stderr)
    print(f"wrote {unmatched_path}", file=sys.stderr)

    readme = os.path.join(os.path.dirname(args.out) or ".",
                          "specimen_taxon_README.txt")
    with open(readme, "w") as handle:
        handle.write(
            "SPECIMEN x TAXON TABLE — what was found in each animal\n"
            "=====================================================\n\n"
            "One row per specimen per marker per taxon detected. Detections only:\n"
            "a row exists because something was found, so absence is the absence of\n"
            "a row rather than a zero.\n\n"
            "COLUMNS\n"
            "  specimen            the animal, as in crosswalk_4marker.tsv\n"
            "  host / host_species collection host; host is the grouped form\n"
            "  n_ticks_in_sample   how many ticks went into the tube\n"
            "  is_single_tick      TRUE only where the collection record gives one\n"
            "                      tick, nothing contradicts it, and the animal was\n"
            "                      wild-caught in the study region\n"
            "  taxon               the organism. For remapped cytb ASVs this is the\n"
            "                      VERIFIED name, not the original classifier call\n"
            "  category            from taxon_decision_table.tsv\n"
            "  resolution_limit    where a species-level call is NOT possible, why\n"
            "  reads / n_asvs      summed across that specimen's libraries\n\n"
            "READ THIS BEFORE USING IT\n"
            "  - NO read threshold has been applied. Every detection is here,\n"
            "    including single-read ones. Extraction blanks carried up to 28\n"
            "    Dermacentor reads, so treat anything in that range as noise.\n"
            "  - Pooled samples cannot give per-tick prevalence: a pathogen found in\n"
            "    a tube of twelve ticks cannot be attributed to any one of them.\n"
            "    Filter is_single_tick == True for per-tick questions.\n"
            "  - Rickettsia and Francisella are ENDOSYMBIONTS, not pathogens, and\n"
            "    Bartonella here is a ked-associated ruminant organism. The category\n"
            "    column says which is which; do not count rows by taxon name alone.\n"
            "  - A category of 'crossamp' means bacterial DNA amplified by the 18S\n"
            "    primers. It is not a detection on that marker.\n")
    print(f"wrote {readme}")


if __name__ == "__main__":
    main()
