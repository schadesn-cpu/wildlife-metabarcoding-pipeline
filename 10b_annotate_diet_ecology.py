#!/usr/bin/env python3
"""
10b_annotate_diet_ecology.py
============================
Add ecological and lay-friendly annotations to a cleaned dietary
metabarcoding taxonomy table produced by 09b_clean_diet_table.py.

Purpose
-------
Taxonomy strings from reference databases (NCBI, MitoFish, PR2) are
accurate but opaque to non-specialist audiences. This script adds three
annotation columns to the cleaned count table:

  habitat       — where the prey lives (e.g. Marine, Freshwater, Anadromous)
  trophic_role  — functional role in the food web (e.g. Forage/schooling,
                  Benthic/demersal, Apex predator)
  common_group  — a single lay-friendly label suitable for figure legends
                  (e.g. "Marine forage fish", "Freshwater prey fish")

The annotated table can then be used to:
  - Generate barplots grouped by common_group rather than species name
  - Produce summary tables showing diet composition by ecological category
  - Support figures for general audiences (conservation reports, press releases)

─────────────────────────────────────────────────────────────────────────────
ADAPTING THIS SCRIPT FOR OTHER STUDY SYSTEMS
─────────────────────────────────────────────────────────────────────────────
This script was written for Common Loon (Gavia immer) gut content
metabarcoding using the MiFish 12S and cytb markers. The annotation
logic is completely separated from the lookup tables, so retrofitting
it to a different study system requires only four steps:

STEP 1 — Replace the species/family lookup table (PREY_LOOKUP below).
  Each entry maps a species name or family name (as it appears in your
  cleaned taxonomy table's first column) to three annotation values.
  You do not need to match every taxon — the script fills unmapped
  taxa with configurable defaults and logs a warning so you can review
  them.

  Example entries for a CANID DIET study (coyote scat metabarcoding):
    "Sylvilagus floridanus":   ("Mammal", "Small mammal prey", "Rabbits and hares"),
    "Lepus americanus":        ("Mammal", "Small mammal prey", "Rabbits and hares"),
    "Peromyscus leucopus":     ("Mammal", "Small mammal prey", "Rodents"),
    "Microtus pennsylvanicus": ("Mammal", "Small mammal prey", "Rodents"),
    "Castor canadensis":       ("Mammal", "Large mammal prey", "Large mammals"),
    "Odocoileus virginianus":  ("Mammal", "Large mammal prey", "Large mammals"),
    "Cervidae":                ("Mammal", "Large mammal prey", "Large mammals"),
    "Meleagris gallopavo":     ("Bird",   "Ground bird prey",  "Ground birds"),
    "Anatinae":                ("Bird",   "Waterfowl",         "Waterfowl"),
    "Vaccinium":               ("Plant",  "Fruit/berry",       "Fruit and berries"),
    "Rosaceae":                ("Plant",  "Fruit/berry",       "Fruit and berries"),
    "Poaceae":                 ("Plant",  "Grass/forb",        "Grasses and forbs"),
    "Homo sapiens":            ("Human",  "Anthropogenic",     "Human food/garbage"),

  The three values are: (habitat_or_origin, trophic_role, common_group)
  Rename habitat to "prey_origin" or "prey_type" via --habitat-label if
  "habitat" doesn't make sense for your system.

STEP 2 — Update the FAMILY_FALLBACK table.
  If a taxon isn't matched at species level, the script tries its family.
  Add your study system's common families here as a fallback.

STEP 3 — Update DEFAULT_HABITAT, DEFAULT_TROPHIC, DEFAULT_GROUP to
  sensible defaults for unmapped taxa in your system. These are used
  as fallback when neither species nor family match.

STEP 4 — Run with --study-system to label outputs and reports clearly
  (e.g. --study-system "Coyote scat MiFish 12S").

That is all that needs to change. The matching logic, output format,
and report generation work identically regardless of the study system.

─────────────────────────────────────────────────────────────────────────────
MATCHING LOGIC
─────────────────────────────────────────────────────────────────────────────
For each row in the cleaned count table, the script tries matches in
priority order:

  1. Exact species name match (case-insensitive) against PREY_LOOKUP
  2. Partial species name match — if the taxon string contains a key
     from PREY_LOOKUP (handles "Urophycis tenuis" matching "tenuis")
  3. Family name match against FAMILY_FALLBACK
  4. Family name parsed from full taxonomy string against FAMILY_FALLBACK
  5. DEFAULT_* values if nothing else matches

All unmatched taxa are logged as warnings and listed in the report so
you can decide whether to add them to the lookup table.

─────────────────────────────────────────────────────────────────────────────
USAGE
─────────────────────────────────────────────────────────────────────────────
  # Loon MiFish (default lookup table)
  python 10b_annotate_diet_ecology.py \\
      --counts   results/MiFish/all/taxonomy_cleaned/taxonomy_counts_cleaned_MiFish.tsv \\
      --marker   MiFish \\
      --outdir   results/MiFish/all/taxonomy_annotated/

  # Loon cytb
  python 10b_annotate_diet_ecology.py \\
      --counts   results/cytb/all/taxonomy_cleaned/taxonomy_counts_cleaned_cytb.tsv \\
      --marker   cytb \\
      --outdir   results/cytb/all/taxonomy_annotated/

  # Different study system with custom lookup table
  python 10b_annotate_diet_ecology.py \\
      --counts        results/coyote/taxonomy_counts_cleaned_MiFish.tsv \\
      --marker        MiFish \\
      --lookup        config/coyote_prey_lookup.tsv \\
      --study-system  "Coyote scat MiFish" \\
      --habitat-label "prey_origin" \\
      --outdir        results/coyote/taxonomy_annotated/

  # Dry run — shows matched/unmatched taxa without writing files
  python 10b_annotate_diet_ecology.py \\
      --counts   results/MiFish/all/taxonomy_cleaned/taxonomy_counts_cleaned_MiFish.tsv \\
      --marker   MiFish \\
      --outdir   results/MiFish/all/taxonomy_annotated/ \\
      --dry-run

  # Export the built-in lookup table as TSV so you can edit it
  python 10b_annotate_diet_ecology.py --export-lookup lookup_template.tsv

Dependencies: pip install pandas
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# LOOKUP TABLE — LOON PREY (MiFish 12S + cytb markers)
# ---------------------------------------------------------------------------
# Format: "species or family name": (habitat, trophic_role, common_group)
#
# To adapt for another study system, replace the contents of this dict or
# pass --lookup pointing to a TSV file with columns:
#   taxon | habitat | trophic_role | common_group
#
# Notes on this loon-specific table:
#   - "habitat" = where the prey lives / migrates through
#   - "trophic_role" = functional feeding guild
#   - "common_group" = single label for lay-audience figures
#
# Anadromous species (alewife, menhaden, shad) are categorised as
# "Anadromous" because their seasonality in loon gut contents directly
# reflects migration phenology — this distinction would be lost if they
# were grouped with either "Marine" or "Freshwater".
#
# "Unresolved" family-level rows (e.g. "Clupeidae (family)") are treated
# as catch-alls for their most ecologically likely constituent species.

PREY_LOOKUP: Dict[str, Tuple[str, str, str]] = {

    # ── Marine forage / schooling ──────────────────────────────────────────
    "Atlantic menhaden":        ("Marine",       "Forage/schooling",   "Marine forage fish"),
    "Brevoortia tyrannus":      ("Marine",       "Forage/schooling",   "Marine forage fish"),
    "Brevoortia patronus":      ("Marine",       "Forage/schooling",   "Marine forage fish"),
    "Clupea harengus":          ("Marine",       "Forage/schooling",   "Marine forage fish"),
    "Atlantic herring":         ("Marine",       "Forage/schooling",   "Marine forage fish"),
    "Pacific herring":          ("Marine",       "Forage/schooling",   "Marine forage fish"),
    "Sprattus sprattus":        ("Marine",       "Forage/schooling",   "Marine forage fish"),
    "Engraulidae":              ("Marine",       "Forage/schooling",   "Marine forage fish"),
    "Anchoa mitchilli":         ("Estuarine",    "Forage/schooling",   "Estuarine forage fish"),
    "Atlantic silverside":      ("Estuarine",    "Forage/schooling",   "Estuarine forage fish"),
    "Menidia menidia":          ("Estuarine",    "Forage/schooling",   "Estuarine forage fish"),
    "Atherinopsidae":           ("Estuarine",    "Forage/schooling",   "Estuarine forage fish"),

    # ── Marine benthic / demersal ─────────────────────────────────────────
    "Urophycis tenuis":         ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "white hake":               ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Phycidae":                 ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Merluccius bilinearis":    ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Merlucciidae":             ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Pseudopleuronectes":       ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Pleuronectidae":           ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Paralichthys dentatus":    ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Paralichthyidae":          ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Cyclopsettidae":           ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Myoxocephalus scorpius":   ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Cottidae":                 ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Centropristis striata":    ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Prionotus carolinus":      ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Sciaenidae":               ("Marine/Estuary","Benthic/demersal",  "Marine bottom fish"),
    "Leiostomus xanthurus":     ("Estuarine",    "Benthic/demersal",   "Estuarine fish"),
    "Pholidae":                 ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Tautogolabrus adspersus":  ("Marine",       "Benthic/demersal",   "Marine bottom fish"),
    "Ctenolabrus rupestris":    ("Marine",       "Benthic/demersal",   "Marine bottom fish"),

    # ── Marine sand lance (ecologically distinct — key loon prey) ─────────
    "Ammodytes americanus":     ("Marine",       "Forage/schooling",   "Sand lance"),
    "Ammodytes dubius":         ("Marine",       "Forage/schooling",   "Sand lance"),
    "Ammodytidae":              ("Marine",       "Forage/schooling",   "Sand lance"),

    # ── Anadromous (key for seasonal ecology narrative) ───────────────────
    "Alosa pseudoharengus":     ("Anadromous",   "Forage/schooling",   "Anadromous fish"),
    "Alosa aestivalis":         ("Anadromous",   "Forage/schooling",   "Anadromous fish"),
    "Alosa sapidissima":        ("Anadromous",   "Forage/schooling",   "Anadromous fish"),
    "Alosa":                    ("Anadromous",   "Forage/schooling",   "Anadromous fish"),
    "alewife":                  ("Anadromous",   "Forage/schooling",   "Anadromous fish"),
    "Salvelinus fontinalis":    ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Salvelinus":               ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Salmonidae":               ("Anadromous",   "Prey fish",          "Anadromous fish"),
    "Morone americana":         ("Anadromous",   "Prey fish",          "Anadromous fish"),
    "Moronidae":                ("Anadromous",   "Prey fish",          "Anadromous fish"),

    # ── Freshwater prey fish ─────────────────────────────────────────────
    "Notemigonus crysoleucas":  ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "golden shiner":            ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Rhinichthys atratulus":    ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Semotilus atromaculatus":  ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Semotilus":                ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Leuciscidae":              ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Fundulus heteroclitus":    ("Estuarine",    "Prey fish",          "Estuarine fish"),
    "Fundulidae":               ("Estuarine",    "Prey fish",          "Estuarine fish"),
    "mummichog":                ("Estuarine",    "Prey fish",          "Estuarine fish"),
    "Lepomis gibbosus":         ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "pumpkinseed":              ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Perca flavescens":         ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "yellow perch":             ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Percidae":                 ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Ameiurus natalis":         ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Ictaluridae":              ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Gasterosteus aculeatus":   ("Estuarine",    "Prey fish",          "Estuarine fish"),
    "threespine stickleback":   ("Estuarine",    "Prey fish",          "Estuarine fish"),
    "Gasterosteidae":           ("Estuarine",    "Prey fish",          "Estuarine fish"),

    # ── Freshwater predators (loons sometimes eat these, notable detections)
    "Micropterus salmoides":    ("Freshwater",   "Apex predator",      "Freshwater predators"),
    "largemouth bass":          ("Freshwater",   "Apex predator",      "Freshwater predators"),
    "Micropterus dolomieu":     ("Freshwater",   "Apex predator",      "Freshwater predators"),
    "Micropterus":              ("Freshwater",   "Apex predator",      "Freshwater predators"),
    "Centrarchidae":            ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Esox americanus":          ("Freshwater",   "Apex predator",      "Freshwater predators"),
    "Esocidae":                 ("Freshwater",   "Apex predator",      "Freshwater predators"),

    # ── Clupeidae catch-all (family-level rows that could be marine/anadromous)
    "Clupeidae":                ("Marine/Anadromous", "Forage/schooling", "Marine or anadromous fish"),

    # ── Unresolved / broad groups ─────────────────────────────────────────
    "Eupercaria":               ("Marine",       "Unresolved",         "Marine fish (unresolved)"),
    "Carangidae":               ("Marine",       "Forage/schooling",   "Marine forage fish"),
    "Apogonidae":               ("Marine",       "Prey fish",          "Marine fish (tropical flag)"),

    # ── cytb-specific non-fish prey ───────────────────────────────────────
    "Pekania pennanti":         ("Terrestrial",  "Secondary prey",     "Terrestrial vertebrates"),
    "fisher":                   ("Terrestrial",  "Secondary prey",     "Terrestrial vertebrates"),
    "Mustelidae":               ("Terrestrial",  "Secondary prey",     "Terrestrial vertebrates"),
    "Cervidae":                 ("Terrestrial",  "Secondary prey",     "Terrestrial vertebrates"),
    "Hyperoartia":              ("Freshwater",   "Prey fish",          "Freshwater prey fish"),
    "Petromyzontidae":          ("Anadromous",   "Prey fish",          "Anadromous fish"),
}


# ---------------------------------------------------------------------------
# FAMILY FALLBACK TABLE
# ---------------------------------------------------------------------------
# Used when species-level match fails. Map family name → annotation.
# This handles "Cottidae (family)" rows from 09b_clean_diet_table.py
# where the species could not be resolved.

FAMILY_FALLBACK: Dict[str, Tuple[str, str, str]] = {
    "Clupeidae":       ("Marine/Anadromous", "Forage/schooling",  "Marine or anadromous fish"),
    "Phycidae":        ("Marine",            "Benthic/demersal",  "Marine bottom fish"),
    "Merlucciidae":    ("Marine",            "Benthic/demersal",  "Marine bottom fish"),
    "Pleuronectidae":  ("Marine",            "Benthic/demersal",  "Marine bottom fish"),
    "Paralichthyidae": ("Marine",            "Benthic/demersal",  "Marine bottom fish"),
    "Cyclopsettidae":  ("Marine",            "Benthic/demersal",  "Marine bottom fish"),
    "Cottidae":        ("Marine",            "Benthic/demersal",  "Marine bottom fish"),
    "Sciaenidae":      ("Marine/Estuary",    "Benthic/demersal",  "Marine bottom fish"),
    "Pholidae":        ("Marine",            "Benthic/demersal",  "Marine bottom fish"),
    "Ammodytidae":     ("Marine",            "Forage/schooling",  "Sand lance"),
    "Atherinopsidae":  ("Estuarine",         "Forage/schooling",  "Estuarine forage fish"),
    "Atherinidae":     ("Marine",            "Forage/schooling",  "Marine forage fish"),
    "Engraulidae":     ("Marine",            "Forage/schooling",  "Marine forage fish"),
    "Leuciscidae":     ("Freshwater",        "Prey fish",         "Freshwater prey fish"),
    "Centrarchidae":   ("Freshwater",        "Prey fish",         "Freshwater prey fish"),
    "Esocidae":        ("Freshwater",        "Apex predator",     "Freshwater predators"),
    "Salmonidae":      ("Anadromous",        "Prey fish",         "Anadromous fish"),
    "Moronidae":       ("Anadromous",        "Prey fish",         "Anadromous fish"),
    "Percidae":        ("Freshwater",        "Prey fish",         "Freshwater prey fish"),
    "Ictaluridae":     ("Freshwater",        "Prey fish",         "Freshwater prey fish"),
    "Gasterosteidae":  ("Estuarine",         "Prey fish",         "Estuarine fish"),
    "Fundulidae":      ("Estuarine",         "Prey fish",         "Estuarine fish"),
    "Carangidae":      ("Marine",            "Forage/schooling",  "Marine forage fish"),
    "Siganidae":       ("Marine (tropical)", "Forage/schooling",  "Tropical fish (flag)"),
    "Lutjanidae":      ("Marine (tropical)", "Benthic/demersal",  "Tropical fish (flag)"),
    "Mustelidae":      ("Terrestrial",       "Secondary prey",    "Terrestrial vertebrates"),
    "Cervidae":        ("Terrestrial",       "Secondary prey",    "Terrestrial vertebrates"),
    "Petromyzontidae": ("Anadromous",        "Prey fish",         "Anadromous fish"),
}

# Defaults for completely unmatched taxa
DEFAULT_HABITAT  = "Unclassified"
DEFAULT_TROPHIC  = "Unclassified"
DEFAULT_GROUP    = "Unclassified prey"


# ---------------------------------------------------------------------------
# Lookup logic
# ---------------------------------------------------------------------------

def _normalise(s: str) -> str:
    """Lowercase, strip whitespace and trailing punctuation for matching."""
    return str(s).strip().lower().rstrip(".")


def annotate_row(
    taxon_str: str,
    lookup: Dict[str, Tuple[str, str, str]],
    family_fallback: Dict[str, Tuple[str, str, str]],
) -> Tuple[str, str, str, str]:
    """
    Return (habitat, trophic_role, common_group, match_method) for a taxon.

    match_method values:
      'exact'   — species name found directly in lookup
      'partial' — lookup key found as substring of taxon_str
      'family'  — family name matched in family_fallback
      'default' — no match, defaults used
    """
    norm = _normalise(taxon_str)

    # 1. Exact match (case-insensitive)
    for key, vals in lookup.items():
        if _normalise(key) == norm:
            return (*vals, "exact")

    # 2. Partial match — key appears anywhere in taxon string
    # Useful for "Urophycis tenuis" matching "tenuis", or taxonomy path
    # strings matching a genus or family keyword
    for key, vals in lookup.items():
        if _normalise(key) in norm:
            return (*vals, "partial")

    # 3. Family-level fallback — extract family name from taxon string
    # Handles "Cottidae (family)" and "k__Metazoa;...;f__Cottidae;..."
    import re
    fam_match = re.search(r'\b(\w+idae|\w+inae)\b', taxon_str)
    if fam_match:
        fam = fam_match.group(1)
        if fam in family_fallback:
            return (*family_fallback[fam], "family")
        # Also try normalised lookup keys
        for key, vals in lookup.items():
            if _normalise(key) == _normalise(fam):
                return (*vals, "family")

    return (DEFAULT_HABITAT, DEFAULT_TROPHIC, DEFAULT_GROUP, "default")


# ---------------------------------------------------------------------------
# Load external lookup TSV (optional --lookup flag)
# ---------------------------------------------------------------------------

def load_lookup_tsv(path: Path) -> Dict[str, Tuple[str, str, str]]:
    """
    Load a custom lookup table from a TSV file.
    Required columns: taxon, habitat, trophic_role, common_group
    Lines starting with # are treated as comments.
    """
    lookup: Dict[str, Tuple[str, str, str]] = {}
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(
            (row for row in f if not row.startswith("#")),
            delimiter="\t"
        )
        required = {"taxon", "habitat", "trophic_role", "common_group"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                f"Lookup TSV must have columns: {required}. "
                f"Found: {reader.fieldnames}"
            )
        for row in reader:
            lookup[row["taxon"].strip()] = (
                row["habitat"].strip(),
                row["trophic_role"].strip(),
                row["common_group"].strip(),
            )
    log.info("Loaded %d entries from custom lookup: %s", len(lookup), path)
    return lookup


def export_lookup_tsv(path: Path) -> None:
    """Export the built-in lookup table as an editable TSV template."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "# Ecological lookup table for 10b_annotate_diet_ecology.py",
            "", "", ""
        ])
        writer.writerow([
            "# Columns: taxon, habitat, trophic_role, common_group",
            "", "", ""
        ])
        writer.writerow([
            "# Edit habitat/trophic_role/common_group for your study system.",
            "", "", ""
        ])
        writer.writerow(["taxon", "habitat", "trophic_role", "common_group"])
        for taxon, (habitat, trophic, group) in sorted(PREY_LOOKUP.items()):
            writer.writerow([taxon, habitat, trophic, group])
    log.info("Exported lookup template to: %s", path)


# ---------------------------------------------------------------------------
# Main annotation pipeline
# ---------------------------------------------------------------------------

def run_annotation(
    counts_path: Path,
    outdir: Path,
    lookup: Dict[str, Tuple[str, str, str]],
    family_fallback: Dict[str, Tuple[str, str, str]],
    habitat_label: str,
    study_system: str,
    dry_run: bool,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("=== 10b_annotate_diet_ecology ===")
    log.info("Input       : %s", counts_path)
    log.info("Study system: %s", study_system)
    log.info("Habitat col : %s", habitat_label)
    if dry_run:
        log.info("DRY RUN — no files will be written")

    # ── Load count table ───────────────────────────────────────────────────
    sep = "\t" if counts_path.suffix in (".tsv", ".txt") else ","
    df = pd.read_csv(counts_path, sep=sep, dtype=str)
    df.columns = df.columns.str.strip()

    # Auto-detect taxon column
    taxon_col = next(
        (c for c in ("Species", "Taxon") if c in df.columns),
        df.columns[0]
    )
    log.info("Taxon column: %s (%d rows)", taxon_col, len(df))

    # ── Annotate ───────────────────────────────────────────────────────────
    habitats, trophics, groups, methods = [], [], [], []
    unmatched: List[str] = []

    for taxon in df[taxon_col]:
        h, t, g, m = annotate_row(str(taxon), lookup, family_fallback)
        habitats.append(h)
        trophics.append(t)
        groups.append(g)
        methods.append(m)
        if m == "default":
            unmatched.append(str(taxon))

    df[habitat_label] = habitats
    df["trophic_role"] = trophics
    df["common_group"] = groups
    df["_match_method"] = methods

    # ── Log match summary ──────────────────────────────────────────────────
    from collections import Counter
    method_counts = Counter(methods)
    log.info("Match summary: %d exact, %d partial, %d family, %d default",
             method_counts["exact"], method_counts["partial"],
             method_counts["family"], method_counts["default"])

    if unmatched:
        log.warning("%d taxa unmatched — using defaults:", len(unmatched))
        for u in unmatched:
            log.warning("  UNMATCHED: %s", u)

    # ── Common group summary ───────────────────────────────────────────────
    log.info("Common group distribution:")
    sample_cols = [c for c in df.columns
                   if c.startswith("TV") and "-GI" in c]
    for group in sorted(df["common_group"].unique()):
        n_taxa = (df["common_group"] == group).sum()
        log.info("  %-35s %d taxa", group, n_taxa)

    # ── Build report ───────────────────────────────────────────────────────
    report_lines = [
        f"=== 10b_annotate_diet_ecology annotation report ===",
        f"Study system : {study_system}",
        f"Input file   : {counts_path}",
        f"Input rows   : {len(df)}",
        f"Taxon column : {taxon_col}",
        f"Habitat col  : {habitat_label}",
        "",
        "Match summary:",
        f"  Exact match  : {method_counts['exact']}",
        f"  Partial match: {method_counts['partial']}",
        f"  Family match : {method_counts['family']}",
        f"  No match     : {method_counts['default']}",
        "",
        "Unmatched taxa (added to lookup to improve coverage):",
    ]
    if unmatched:
        for u in unmatched:
            report_lines.append(f"  {u}")
    else:
        report_lines.append("  (none — all taxa matched)")

    report_lines += [
        "",
        "Annotation table (taxon → common_group):",
    ]
    for _, row in df[[taxon_col, "common_group", "_match_method"]].iterrows():
        report_lines.append(
            f"  {str(row[taxon_col]):<45} → {row['common_group']}"
            f"  [{row['_match_method']}]"
        )

    report_lines += [
        "",
        "Common group summary (for figure legends):",
    ]
    for group in sorted(df["common_group"].unique()):
        taxa = df[df["common_group"] == group][taxon_col].tolist()
        report_lines.append(f"  {group}:")
        for t in taxa:
            report_lines.append(f"    - {t}")

    report_lines += [
        "",
        "Next step — generate annotated barplot:",
        f"  python 09_plot_taxonomy.py \\",
        f"    --relabund  <relabund_file> \\",
        f"    --annotation {outdir / f'annotation_table.tsv'} \\",
        f"    --group-col common_group \\",
        f"    --metadata  metadata/qiime/metadata_MiFish.tsv \\",
        f"    --outdir    {outdir}",
    ]

    report_text = "\n".join(report_lines) + "\n"
    log.info("\n" + report_text) if dry_run else None

    # ── Write outputs ──────────────────────────────────────────────────────
    if not dry_run:
        # Full annotated count table
        out_counts = outdir / f"taxonomy_counts_annotated_{counts_path.stem.split('_')[-1]}.tsv"
        df_out = df.drop(columns=["_match_method"])
        df_out.to_csv(out_counts, sep="\t", index=False)
        log.info("Written: %s", out_counts)

        # Compact annotation lookup (taxon → annotations only, no sample counts)
        out_annot = outdir / "annotation_table.tsv"
        annot_cols = [taxon_col, habitat_label, "trophic_role", "common_group"]
        df_annot = df[annot_cols].drop_duplicates()
        df_annot.to_csv(out_annot, sep="\t", index=False)
        log.info("Written: %s", out_annot)

        # Report
        out_report = outdir / "annotation_report.txt"
        out_report.write_text(report_text, encoding="utf-8")
        log.info("Written: %s", out_report)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="10b_annotate_diet_ecology.py",
        description=(
            "Add ecological annotations to a cleaned dietary metabarcoding "
            "taxonomy table. Produces habitat, trophic_role, and common_group "
            "columns for lay-friendly figure labels. Fully adaptable to any "
            "prey-diet study system via --lookup."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--counts", type=Path,
                   help="Cleaned taxonomy count TSV from 09b_clean_diet_table.py")
    p.add_argument("--marker", default="MiFish",
                   help="Marker name for output file naming (default: MiFish)")
    p.add_argument("--outdir", type=Path,
                   help="Output directory")
    p.add_argument("--lookup", type=Path, default=None,
                   help=(
                       "Optional TSV lookup table to override the built-in "
                       "loon prey table. Must have columns: taxon, habitat, "
                       "trophic_role, common_group. Lines starting with # "
                       "are treated as comments. See --export-lookup to "
                       "generate a starting template."
                   ))
    p.add_argument("--study-system", default="Common Loon MiFish/cytb",
                   help="Study system description for report headers")
    p.add_argument("--habitat-label", default="habitat",
                   help=(
                       "Column name for the habitat/origin annotation. "
                       "Override for non-aquatic study systems, e.g. "
                       "'prey_origin' for canid scat or 'foraging_zone' "
                       "for bat diet studies. (default: habitat)"
                   ))
    p.add_argument("--export-lookup", type=Path, default=None,
                   help=(
                       "Export the built-in lookup table as an editable TSV "
                       "template to this path, then exit. Use this as a "
                       "starting point when adapting for a new study system."
                   ))
    p.add_argument("--dry-run", action="store_true",
                   help="Show annotation results without writing files")
    return p


def main() -> int:
    args = build_parser().parse_args()

    # Export lookup template and exit
    if args.export_lookup:
        export_lookup_tsv(args.export_lookup)
        return 0

    if not args.counts:
        print("Error: --counts is required unless using --export-lookup",
              file=sys.stderr)
        return 2

    if not args.counts.exists():
        log.error("Input file not found: %s", args.counts)
        return 2

    # Load lookup table — custom file or built-in
    if args.lookup:
        lookup = load_lookup_tsv(args.lookup)
        family_fallback: Dict[str, Tuple[str, str, str]] = {}
    else:
        lookup = PREY_LOOKUP
        family_fallback = FAMILY_FALLBACK

    try:
        run_annotation(
            counts_path   = args.counts,
            outdir        = args.outdir,
            lookup        = lookup,
            family_fallback = family_fallback,
            habitat_label = args.habitat_label,
            study_system  = args.study_system,
            dry_run       = args.dry_run,
        )
    except Exception as e:
        log.error("Annotation failed: %s", e)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
