#!/usr/bin/env python3
"""The single source of truth for what each taxon IS in this project.

Every call we have made about a taxon currently lives somewhere different: pathogen
vs endosymbiont is hardcoded in build_analysis_table.py, the Yarza identity floors
are in a comment in apply_16S_drop_list.py, the Bartonella and ked reclassifications
are in the ground truth as prose, and the Borreliella species limit is in the
manuscript. That is four places, and they have already drifted from each other more
than once.

This writes one table. Scripts and figures read it instead of re-deciding, so they
cannot disagree.

I key on FULL LINEAGE, not on the display label. The label repeats: 'uncultured'
appears under several families and means something different each time, and
'Anaplasma' is a substring of 'Anaplasmataceae', which is exactly how the
SILVA_genus column in pathogen_verification_report.tsv came to say Anaplasma for
eleven Wolbachia ASVs. A label key would reproduce that bug.

Nothing is filtered. Every taxon in the corrected taxonomies gets a row, with its
category and the reason. Someone reading this should be able to see everything and
disagree with a call, rather than getting a view where I already decided.

Run:
  cd ~/meedlab/tick_project
  python3 scripts/make_taxon_decision_table.py \
    --out results/tables/taxon_decision_table.tsv

Writes:
  results/tables/taxon_decision_table.tsv
"""

import argparse
import csv
import os
import re
import sys

# ── the corrected taxonomies, one per marker ──────────────────────────────────
# I prefer the drop-list-applied 16S taxonomy where it exists and fall back to the
# BLAST-verified one otherwise, so this table cannot silently read a taxonomy that
# predates the corrections. If it falls back, it says so.
TAXONOMY_PATHS = {
    "16S": ["results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv",
            "results/16S/all/taxonomy_verified/taxonomy_16S_blast_verified.tsv"],
    "18S": ["results/18S/all/taxonomy_verified/taxonomy_18S_blast_verified.tsv"],
    # The FULL verified cytb taxonomy, matching make_asv_tables.py. clean.tsv is the
    # 3,709 that survived filtering and does not contain the human lineages, so
    # reading it here left 162 human cytb detections uncategorised in a table I was
    # about to share.
    "cytb": ["results/cytb/all/taxonomy_verified/taxonomy_blast_verified.tsv",
             "results/cytb/all/taxonomy_verified/taxonomy_clean.tsv"],
}

# ── categories ────────────────────────────────────────────────────────────────
# I keep endosymbionts and pathogens apart because conflating them is the single
# most consequential error available here: Rickettsia is 64.7% of the deer bar and
# is NOT a pathogen call, and Francisella is 45.8% of the moose bar and is a
# Francisella-like endosymbiont, not F. tularensis.
CATEGORY_RULES = [
    # (category, reason, [genus or lineage substrings])
    ("pathogen_bacterial",
     "Bacterial pathogen on the study panel. Reported as a detection frequency, "
     "not as confirmed infection of the tick.",
     ["g__Anaplasma", "g__Borreliella", "g__Borrelia", "g__Ehrlichia"]),

    ("pathogen_protozoan",
     "Protozoan pathogen. Babesia species separable at 220 bp V4 (BLAST-verified "
     "2026-08-10): B. microti and B. odocoilei do not tie.",
     ["g__Babesia", "Babesia"]),

    ("piroplasm_nonhuman",
     "Piroplasm detected but not on the qPCR panel and not a human pathogen in "
     "this system. Theileria cervi is a deer piroplasm; both ASVs came from "
     "24-3775, the out-of-region Martha's Vineyard deer excluded from the study. "
     "Hepatozoon canis was 7 reads in one library and is vectored by "
     "Rhipicephalus sanguineus, which is not present in this dataset.",
     ["Theileria", "Hepatozoon"]),

    ("offpanel_pathogen_genus",
     "Pathogen genus NOT on the qPCR panel. Genus-level only; the identity is "
     "below the 98.7% species threshold. A Goal 3 observation: the agnostic method "
     "flagged an organism the targeted assay does not cover.",
     ["Neoehrlichia", "Candidatus_Neoehrlichia"]),

    ("endosymbiont",
     "Tick endosymbiont. NOT a pathogen and never counted as one. Rickettsia "
     "buchneri is the I. scapularis endosymbiont; Francisella-like endosymbionts "
     "are characteristic of Dermacentor. All 26 Francisella ASVs BLAST to "
     "endosymbiont references; F. tularensis is the best hit for none of them.",
     ["g__Rickettsia", "g__Francisella", "g__Rickettsiella", "g__Arsenophonus",
      "g__Wolbachia", "Cardinium", "Diplorickettsia", "g__Coxiella",
      "Midichloria", "Spiroplasma"]),

    ("ectoparasite_codetection",
     "Co-collected ectoparasite, not a blood meal and not a pathogen. The deer ked "
     "was recorded in the accessioning notes for 29 specimens and detected on "
     "three markers.",
     ["Lipoptena", "Hippobosc", "Melophagus", "Icosta"]),

    ("ked_associated_bacterium",
     "Ruminant-associated Bartonella, hippoboscid-linked. RECLASSIFIED 2026-08-04: "
     "NOT a tick-borne human pathogen. 34 of 45 ASVs tie B. schoenbuchensis; two "
     "more hit B. melophagi. Reported as a co-detection, outside the pathogen "
     "results.",
     ["g__Bartonella"]),

    ("nontarget_tick_dna_cytb",
     "The tick's OWN genome amplified by the cytb primers and assigned to the nearest "
     "available vertebrate, because the classifier is vertebrate-only and has no tick "
     "reference. BLAST-verified April 2026: the 89 Gobiidae ASVs are Ixodes ricinus "
     "genomic DNA (OZ285889.1). Actinopteri and the Rangifer tarandus genome hits (NUMTs) "
     "are the same artefact. NOT a blood meal and NOT contamination from another project.",
     ["Gobiidae", "Actinopteri", "Actinopterygii", "Teleostei", "Petromyzont",
      "Rangifer tarandus"]),

    ("classifier_misassignment_cytb",
     "Assigned to a taxon with no local reference, corrected by accession-level Entrez "
     "verification in April 2026. Plethodontidae and Hydromantes (about 153 ASVs) are the "
     "deer ked Lipoptena cervi at 94.8 to 99.1% identity (PV090979.1, PV090977.1); the "
     "classifier had no Hippoboscidae references. Axis porcinus annamiticus (1,078 ASVs) "
     "is North American cervid. Use remap_verified_taxon, never the original lineage.",
     ["Plethodontidae", "Hydromantes", "Axis porcinus", "Dama mesopotamica"]),

    ("tick_host",
     "The tick itself. Used for species assignment, excluded from community "
     "composition.",
     ["Ixodes", "Dermacentor", "Amblyomma", "Haemaphysalis", "Rhipicephalus",
      "g__Acari", "Arachnida"]),

    ("human_dna",
     "Human mitochondrial cytb. This CANNOT distinguish a human blood meal from "
     "handler or laboratory contamination, and it must not be reported as a blood "
     "meal. Human is 0.41% of all cytb reads and 252 of 837 biological samples "
     "carry trace human, which is a contamination background rather than a "
     "finding. Treat an individual detection as 'human-positive above baseline'. "
     "See Human_Checkpoint_Notes, 2026-07-09.",
     ["Homo", "Homo_sapiens", "Homo sapiens"]),

    ("vertebrate_host",
     "Vertebrate DNA: the blood meal or the collection host. Excluded from the "
     "eukaryotic community table. Abundant in 18S from engorged ticks (7.74% of "
     "the 18S run).",
     ["Odocoileus", "Craniata", "Chordata", "Vertebrata", "Mammalia", "Aves",
      "Archosauria", "Ursus", "Alces", "Didelphis", "Erethizon",
      "Tarsius", "Cervidae"]),

    ("crossamp",
     "Bacterial sequence amplified by the 18S primers. 40.4% of the 18S run. "
     "Excluded from the eukaryotic community.",
     ["CROSSAMP"]),

    ("reagent_contaminant_flagged",
     "Flagged by decontam as more abundant in negatives than in samples. NOTE: "
     "decontam works PER ASV, not per genus, so a genus can have both flagged and "
     "genuine ASVs. 51 ASVs were removed from the decontam table; Lactococcus was "
     "NOT flagged (p = 0.189) despite an earlier caption saying so.",
     ["Escherichia-Shigella", "g__Pseudomonas"]),

    ("mitochondrial",
     "Mitochondrial reference. Off-target for an SSU survey whatever organism it "
     "came from.",
     [":mito"]),
]

# Species-level calls we have decided we CANNOT make, with the reason. These ride
# alongside the category rather than replacing it.
RESOLUTION_LIMITS = {
    "g__Borreliella":
        "Genus only. All 29 Borreliella ASVs tie B. burgdorferi with B. bissettiae "
        "at identical percent identity (confirmed at -max_target_seqs 50 on "
        "full-length queries, 2026-08-04). V4-V5 carries no species-level "
        "information for the Lyme group. Species-level Lyme calls rest on qPCR.",
    "g__Anaplasma":
        "Genus only from 16S. Ap-Variant 1 and the human-active Ap-ha strain "
        "differ by 2 bp in 16S rRNA; V4-V5 cannot be assumed to resolve them. "
        "Ap-V1 is deer-associated and not linked to human disease.",
    "Cervidae":
        "Family only. cytb cannot separate white-tailed deer from moose.",
    "Amblyomma":
        "Genus only. 18S does not give species; COI was proposed but never run.",
}

# Identity floors, Yarza et al. 2014. Recorded here so every script uses the same
# numbers rather than each carrying its own copy.
GENUS_IDENTITY_FLOOR = 94.5
SPECIES_IDENTITY_FLOOR = 98.7

UNASSIGNED_TOKENS = {"", "unclassified", "unassigned", "na", "none", "__"}


def back_up_existing(out_path):
    """Moves an existing output aside to .bak before I overwrite it.

    I don't clobber an existing backup either, so a second rerun stops rather than
    losing the original.
    """
    if not os.path.exists(out_path):
        return
    backup_path = out_path + ".bak"
    if os.path.exists(backup_path):
        sys.exit(f"{out_path} and {backup_path} both exist. "
                 f"Move or delete {backup_path} before rerunning.")
    os.replace(out_path, backup_path)
    print(f"backed up previous output to {backup_path}", file=sys.stderr)


def clean_rank(rank):
    """One lineage rank, or empty if it was never assigned."""
    rank = rank.strip()
    rank = re.sub(r"^[a-z]__", "", rank)
    rank = rank.strip("_")
    if rank.lower() in UNASSIGNED_TOKENS:
        return ""
    return rank


def leaf_label(lineage):
    """Deepest NAMED rank, bare.

    A bare species epithet is meaningless on its own: 's__microti' could be Babesia
    or Cryptosporidium. Where the deepest rank is a species I prepend the genus.
    """
    parts = [part.strip() for part in str(lineage).split(";")]
    named = [part for part in parts if clean_rank(part)]
    if not named:
        return "Unassigned"
    last = named[-1]
    # 'Homo;sapiens environmental sample' leafs to 'sapiens environmental sample',
    # which is both unreadable and misleading. Where the deepest rank is a bare
    # epithet or an NCBI 'environmental sample' string, I prepend the parent rank.
    if len(named) >= 2 and ("environmental sample" in last.lower()
                            or clean_rank(last).islower()):
        parent = clean_rank(named[-2])
        epithet = clean_rank(last).replace(" environmental sample", "")
        if parent and epithet and not epithet.lower().startswith(parent.lower()[:4]):
            return f"{parent} {epithet}"
    if last.startswith("s__") and len(named) >= 2:
        genus = clean_rank(named[-2])
        epithet = clean_rank(last)
        if genus and not epithet.lower().startswith(genus.lower()[:4]):
            return f"{genus} {epithet}"
        return epithet
    return clean_rank(last)


def lineage_ranks(lineage):
    """Every named rank of a lineage as a list, bare and lowercased.

    I match against these rather than against the raw string. A substring test on
    the whole lineage put 67 'uncultured' taxa into endosymbiont, because their
    lineage passes through an endosymbiont family on the way down. It is the same
    failure as the SILVA_genus column, which called eleven Wolbachia ASVs Anaplasma
    because 'Anaplasma' is a substring of 'Anaplasmataceae'.
    """
    return [clean_rank(part).lower()
            for part in str(lineage).split(";") if clean_rank(part)]


RANK_MARKERS = [("species", ("s__",)), ("genus", ("g__",)), ("family", ("f__",)),
                ("order", ("o__",)), ("class", ("c__",)), ("phylum", ("p__",)),
                ("domain", ("d__", "k__"))]


def rank_reached(lineage):
    """The deepest rank this lineage actually reached.

    Where the lineage carries no rank prefixes (cytb writes bare binomials) I say
    'unprefixed' rather than guessing from position.
    """
    parts = [part.strip() for part in str(lineage).split(";")]
    if not any(re.match(r"^[a-z]__", part) for part in parts):
        return "unprefixed"
    for name, prefixes in RANK_MARKERS:
        for part in parts:
            if part.startswith(prefixes) and clean_rank(part):
                return name
    return "unassigned"


def categorise(lineage):
    """Returns (category, reason) for a lineage, or ('', '') if no rule matches.

    Rules are checked IN ORDER and the first match wins, so the ordering in
    CATEGORY_RULES is load-bearing: 'ked_associated_bacterium' must be checked
    before anything that would catch Bartonella as a pathogen.

    A pattern matches only if it equals a whole rank, or if it is a prefix marker
    like 'CROSSAMP' or ':mito' that is not a rank at all. Substring matching against
    the whole lineage is what miscategorised 'uncultured'.
    """
    ranks = lineage_ranks(lineage)
    raw = str(lineage)

    # CROSSAMP first, unconditionally. It is a statement about the whole lineage,
    # not a taxonomic rank, and it outranks every taxonomic call: a CROSSAMP
    # Anaplasma in the 18S run is bacterial DNA amplified by the eukaryotic
    # primers, NOT a pathogen detection on that marker. Checking it inside the
    # pattern loop let 'Anaplasma' match a rank first and put 60 cross-amplified
    # ASVs into pathogen_bacterial.
    # Eukaryotic sequence in the 16S run. This is the mirror image of CROSSAMP and it was
    # missing: SILVA contains 18S, the 515F/926R primers co-amplify it, and 9,333 ASVs in
    # the 16S taxonomy are eukaryotic. 1,217 are the tick itself (o__Acari), 26 the deer ked
    # (o__Diptera), 73 the host blood meal (c__Mammalia), and 6,648 resolve no deeper than
    # d__Eukaryota. Every one is CORRECTLY labelled, which is why no BLAST verification
    # flagged it, but none is a member of a bacterial community. Left in, they sit in the
    # denominator of every relative abundance computed from the 16S table.
    if raw.startswith("d__Eukaryota") or ";d__Eukaryota" in raw:
        return ("eukaryotic_in_16S",
                "Eukaryotic sequence recovered by the 16S primers, which co-amplify 18S. "
                "The tick itself, the deer ked, the host blood meal, and eukaryotic ASVs "
                "resolved no further. Correctly identified, but NOT part of the bacterial "
                "community: exclude from any 16S relative abundance or diversity measure. "
                "Mirror image of the crossamp category, which covers bacteria appearing in "
                "the 18S run.")

    if raw.startswith("CROSSAMP") or ";CROSSAMP" in raw:
        return ("crossamp",
                "Bacterial sequence amplified by the 18S primers. 40.4% of the 18S "
                "run. NOT a detection on this marker: the same organism is called "
                "from 16S, where it belongs. Excluded from the eukaryotic community.")

    # A genus field that says 'unassigned_<family>' is the classifier telling us it
    # could not resolve a genus. The rank prefix reads g__, but the content means
    # family-level, and it must not be counted as a genus-level call.
    if any(rank.startswith("unassigned_") or rank.startswith("unclassified_")
           for rank in ranks):
        return ("rank_above_genus",
                "The genus field records no resolved genus (unassigned_<family>), "
                "so no genus-level call is possible. Produced by the identity-floor "
                "corrections: an ASV whose best match in nt is below the 94.5% "
                "genus threshold (Yarza et al. 2014) is relabelled to family.")

    # A lineage that stops at family or above cannot carry a genus-level call. I
    # test this BEFORE the pattern rules but only for the pathogen-adjacent
    # families, because a family-level 'Anaplasmataceae' is not a pathogen call and
    # must not be counted as one. Testing it as a normal pattern rule swallowed
    # Wolbachia and Rickettsia, whose lineages pass through those same families on
    # the way down to a real genus.
    AMBIGUOUS_FAMILIES = {"anaplasmataceae", "rickettsiaceae", "francisellaceae",
                          "babesiidae", "borreliaceae", "rickettsiales",
                          "piroplasmida", "piroplasmorida", "theileriidae"}
    if rank_reached(lineage) in ("family", "order", "class", "phylum", "domain") \
            or (rank_reached(lineage) == "unprefixed" and ranks and ranks[-1] in AMBIGUOUS_FAMILIES):
        if ranks and ranks[-1] in AMBIGUOUS_FAMILIES:
            return ("rank_above_genus",
                    "Assigned no deeper than family or order, so no genus-level "
                    "call is possible. Kept visible rather than dropped: a "
                    "family-level Anaplasmataceae or Rickettsiaceae assignment is "
                    "NOT a pathogen call and must not be counted as one.")

    for category, reason, patterns in CATEGORY_RULES:
        for pattern in patterns:
            bare = pattern.replace("g__", "").replace("f__", "").lower()
            if not pattern.startswith(("g__", "f__")) and not pattern[0].isalpha():
                if pattern in raw:                       # CROSSAMP, :mito
                    return category, reason
                continue
            for rank in ranks:
                # exact rank, or the rank with a trailing _sp. / _X_sp. / _cervi
                if rank == bare or rank.startswith(bare + "_") or rank.startswith(bare + " "):
                    return category, reason
            if pattern in ("CROSSAMP", ":mito") and pattern in raw:
                return category, reason
    return "", ""


def resolution_limit(lineage):
    """Returns the recorded resolution limit for this lineage, or ''."""
    for pattern, note in RESOLUTION_LIMITS.items():
        if pattern in lineage:
            return note
    return ""


def load_taxonomy(paths, marker):
    """Returns [(feature_id, lineage, confidence)] from the first path that exists.

    Takes a list so the caller can express a preference order. I print which one I
    used, because reading a pre-correction taxonomy without noticing is exactly the
    failure this table is meant to prevent.
    """
    path = next((candidate for candidate in paths if os.path.isfile(candidate)), None)
    if path is None:
        print(f"warning: no {marker} taxonomy found; tried {paths}", file=sys.stderr)
        return []
    if path != paths[0]:
        print(f"NOTE: {marker} fell back to {path} because {paths[0]} does not "
              f"exist. That file predates the drop-list corrections.", file=sys.stderr)
    else:
        print(f"  {marker}: {path}", file=sys.stderr)
    rows = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            print(f"warning: {path} is empty", file=sys.stderr)
            return []
        for row in reader:
            if len(row) < 2 or not row[0].strip():
                continue
            confidence = row[2].strip() if len(row) > 2 else ""
            rows.append((row[0].strip(), row[1].strip(), confidence))
    if not rows:
        print(f"warning: {path} has a header but no data rows", file=sys.stderr)
    return rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="results/tables/taxon_decision_table.tsv")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    # One row per distinct lineage per marker. I collapse feature IDs here because
    # this table is about what a TAXON is; the per-ASV detail belongs in the ASV
    # table, which is a separate file.
    rows_by_key = {}
    for marker, path in TAXONOMY_PATHS.items():
        for feature_id, lineage, confidence in load_taxonomy(path, marker):
            key = (marker, lineage)
            if key not in rows_by_key:
                category, reason = categorise(lineage)
                rows_by_key[key] = {
                    "marker": marker,
                    "taxon": leaf_label(lineage),
                    "category": category,
                    "reason": reason,
                    "resolution_limit": resolution_limit(lineage),
                    "is_pathogen": "yes" if category.startswith("pathogen") else (
                        "no" if category else ""),
                    "in_community_table": "no" if category in (
                        "tick_host", "vertebrate_host", "crossamp", "mitochondrial",
                        "ectoparasite_codetection") else "yes",
                    "rank_resolved": rank_reached(lineage),
                    "n_asvs": 0,
                    "full_lineage": lineage,
                }
            rows_by_key[key]["n_asvs"] += 1

    if not rows_by_key:
        sys.exit("no taxonomies could be read; check the paths in TAXONOMY_PATHS")

    out_rows = sorted(
        rows_by_key.values(),
        key=lambda row: (row["marker"], row["category"] or "zz_uncategorised",
                         -row["n_asvs"], row["taxon"]))

    back_up_existing(args.out)
    with open(args.out, "w", newline="") as handle:
        # lineterminator="\n": the csv default is \r\n, which leaves a stray \r on
        # the last field and quietly breaks awk and grep on ron.
        writer = csv.DictWriter(handle, delimiter="\t",
                                fieldnames=list(out_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"wrote {args.out}")
    print(f"  {len(out_rows)} distinct taxa across {len(TAXONOMY_PATHS)} markers")
    print(f"  identity floors recorded: genus >= {GENUS_IDENTITY_FLOOR}%, "
          f"species >= {SPECIES_IDENTITY_FLOOR}% (Yarza et al. 2014)")
    print()
    print("  by category:")
    counts = {}
    for row in out_rows:
        label = row["category"] or "(no rule matched)"
        counts[label] = counts.get(label, 0) + 1
    for label, n in sorted(counts.items(), key=lambda item: -item[1]):
        print(f"    {label:32s} {n:>5}")

    uncategorised = [row for row in out_rows if not row["category"]]
    if uncategorised:
        print(f"\nnote: {len(uncategorised)} taxa matched no rule, so their category "
              f"and is_pathogen are blank. That is expected: most of the community "
              f"is neither a pathogen nor an endosymbiont and does not need a call. "
              f"They are in the file with everything else.", file=sys.stderr)


if __name__ == "__main__":
    main()
