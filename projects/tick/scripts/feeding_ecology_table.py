#!/usr/bin/env python3
"""The feeding-ecology cross-table: tick species x host x blood meal x pathogens.

D. Needle's task list asks for this layout directly: for each tick species, on each
collection host, what did it feed on and which of the seven panel pathogens did it carry,
with qPCR and metabarcoding side by side. Every component already exists; nothing here is a
new measurement, it is an assembly.

One row per specimen. Wide rather than long, because this table is meant to be read across:
someone should be able to find a tick and see everything known about it without joining
anything. The long form for filtering is already `specimen_taxon_table.tsv`.

WHAT IS AND IS NOT IN A CELL
  - qPCR columns are per-species assays. A blank means the assay was not run on that
    specimen, NOT that it was negative; the qpcr_screened column says which.
  - metabarcoding columns are per-genus for 16S and per-species for 18S. Borreliella cannot
    be resolved to species by the V4-V5 fragment, so a Borreliella-positive tick is not
    necessarily B. burgdorferi-positive.
  - blood meal is cytb. Cervidae merges deer, moose and unresolved cervid, because cytb
    cannot separate them.
  - deer ked is a co-detection, not a blood meal, and has its own column.
  - human cytb is NOT included: the highest human read count in the study is in a positive
    control, so it is contamination rather than feeding.

Run:
    cd ~/meedlab/tick_project
    python3 scripts/feeding_ecology_table.py \
        --out results/analysis/feeding_ecology_table.tsv

Writes:
    results/analysis/feeding_ecology_table.tsv    one row per specimen
    results/analysis/feeding_ecology_summary.txt  the cross-tabulations
"""

import argparse
import csv
import os
import re
import sys

import pandas as pd

CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
SINGLE_TICK = "results/single_tick_v3/single_tick_v3_full.tsv"
TICK_GENUS = "results/18S/all/taxonomy_verified/tick_genus_per_sample_corrected.tsv"
WORKBOOK = "metadata/tick_full_metadata.xlsx"
QPCR_SHEET = "qPCR Results"

QPCR_TARGETS = {
    "a-phagocytophilum": "qpcr_A_phagocytophilum",
    "b-burgdorferi": "qpcr_B_burgdorferi",
    "b-microti": "qpcr_B_microti",
    "b-miyamotoi": "qpcr_B_miyamotoi",
    "powassan-virus": "qpcr_Powassan",
    "b-mayonii": "qpcr_B_mayonii",
    "e-muris": "qpcr_E_muris",
}

SEQ_PATHOGENS = {
    "16S_Anaplasma": "seq_Anaplasma",
    "16S_Borreliella": "seq_Borreliella",
    "16S_Borrelia": "seq_Borrelia",
    "18S_Babesia_odocoilei": "seq_Babesia_odocoilei",
    "18S_Babesia_microti": "seq_Babesia_microti",
}

# Bartonella is kept apart from the pathogen block: it is B. schoenbuchensis, a
# ruminant-associated organism carried by deer keds, not a tick-borne human pathogen.
SEQ_CODETECTIONS = {"16S_Bartonella": "codetect_Bartonella"}

BLOODMEAL_COLUMNS = {
    "cytb_Deer": "Cervidae", "cytb_Moose": "Cervidae",
    "cytb_Cervidae_unresolved": "Cervidae",
    "cytb_Bird": "Bird", "cytb_Bear": "Bear", "cytb_Opossum": "Opossum",
    "cytb_Raccoon": "Raccoon", "cytb_Porcupine": "Porcupine", "cytb_Rabbit": "Rabbit",
}
KED_COLUMN = "cytb_Deer_ked"

report_lines = []


def say(message=""):
    print(message)
    report_lines.append(str(message))


def normalise(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def is_positive(value):
    return str(value).strip() in ("1", "1.0")


def back_up_existing(out_path):
    if not os.path.exists(out_path):
        return
    backup_path = out_path + ".bak"
    if os.path.exists(backup_path):
        sys.exit(f"{out_path} and {backup_path} both exist. "
                 f"Move or delete {backup_path} before rerunning.")
    os.replace(out_path, backup_path)


def load_tick_genus():
    """barcode -> tick genus, reconciling a specimen's libraries.

    sample_id is a LIBRARY name with a run suffix; the barcode is the part before the first
    dot. A specimen with several libraries can disagree with itself, and in a pooled sample
    that is a real thing rather than an error, so disagreement becomes 'both/mixed'.
    """
    if not os.path.exists(TICK_GENUS):
        print(f"warning: {TICK_GENUS} not found; tick_genus will be blank", file=sys.stderr)
        return {}
    calls = {}
    with open(TICK_GENUS, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            raw = (row.get("sample_id") or "").strip()
            value = (row.get("dominant") or "").strip()
            if raw and value:
                calls.setdefault(normalise(raw.split(".")[0]), set()).add(value)
    genus = {}
    for key, values in calls.items():
        ticks = {v for v in values if v in ("Ixodes", "Dermacentor", "Amblyomma")}
        genus[key] = (ticks.pop() if len(ticks) == 1
                      else "both/mixed" if len(ticks) > 1 else "")
    return genus


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="results/analysis/feeding_ecology_table.tsv")
    args = parser.parse_args()

    for path in (CROSSWALK, SINGLE_TICK, WORKBOOK):
        if not os.path.isfile(path):
            sys.exit(f"no such file: {path}")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    crosswalk = list(csv.DictReader(open(CROSSWALK, newline=""), delimiter="\t"))
    say(f"crosswalk: {len(crosswalk)} specimens")

    single_ticks = {normalise(r["Barcode"]): r for r in
                    csv.DictReader(open(SINGLE_TICK, newline=""), delimiter="\t")}
    say(f"single ticks with calls: {len(single_ticks)}")

    genus_by_barcode = load_tick_genus()
    say(f"specimens with a tick-genus call: {sum(1 for v in genus_by_barcode.values() if v)}")

    qpcr_raw = pd.read_excel(WORKBOOK, sheet_name=QPCR_SHEET, dtype=str)
    missing = [t for t in QPCR_TARGETS if t not in qpcr_raw.columns]
    if missing:
        sys.exit(f"{QPCR_SHEET} is missing target columns: {missing}")
    qpcr = {normalise(r["host_id"]): r for _, r in qpcr_raw.iterrows()}
    say(f"qPCR rows: {len(qpcr)}")

    rows = []
    for specimen in crosswalk:
        barcode = normalise(specimen.get("barcode") or specimen.get("sample_id") or "")
        if not barcode:
            continue
        tick = single_ticks.get(barcode)
        qpcr_row = qpcr.get(barcode)

        row = {
            "specimen": specimen.get("sample_id", ""),
            "barcode": specimen.get("barcode", ""),
            "tick_genus": genus_by_barcode.get(barcode, ""),
            "collection_host": specimen.get("host_bin", ""),
            "host_species": specimen.get("host_spp", ""),
            "state": specimen.get("state", ""),
            "town": specimen.get("town", ""),
            "n_ticks_in_sample": specimen.get("num_ticks", ""),
            "is_single_tick": specimen.get("is_single_tick", ""),
            "sequenced": "yes" if tick else "no",
            "qpcr_screened": "yes" if qpcr_row is not None else "no",
        }

        # blood meal, collapsed
        meals = set()
        if tick:
            for column, label in BLOODMEAL_COLUMNS.items():
                if float(tick.get(column) or 0) > 0:
                    meals.add(label)
        row["blood_meal"] = ";".join(sorted(meals)) if meals else (
            "none detected" if tick else "")
        row["deer_ked_codetected"] = (
            "yes" if tick and float(tick.get(KED_COLUMN) or 0) > 0 else
            ("no" if tick else ""))

        # metabarcoding pathogens
        n_seq = 0
        for column, label in SEQ_PATHOGENS.items():
            if tick is None:
                row[label] = ""
            else:
                hit = float(tick.get(column) or 0) > 0
                row[label] = "1" if hit else "0"
                n_seq += int(hit)
        row["n_pathogens_seq"] = n_seq if tick else ""
        for column, label in SEQ_CODETECTIONS.items():
            row[label] = ("" if tick is None else
                          ("1" if float(tick.get(column) or 0) > 0 else "0"))

        # qPCR
        n_qpcr = 0
        for target, label in QPCR_TARGETS.items():
            if qpcr_row is None:
                row[label] = ""
            else:
                hit = is_positive(qpcr_row.get(target))
                row[label] = "1" if hit else "0"
                n_qpcr += int(hit)
        row["n_targets_qpcr"] = n_qpcr if qpcr_row is not None else ""

        rows.append(row)

    table = pd.DataFrame(rows)
    back_up_existing(args.out)
    table.to_csv(args.out, sep="\t", index=False)
    say(f"\nwrote {args.out}  ({len(table)} specimens)")

    # ---- the cross-tabulations D. Needle asked for -------------------------
    say("\n" + "=" * 70)
    say("TICK SPECIES x COLLECTION HOST x BLOOD MEAL")
    say("=" * 70)
    described = table[(table["tick_genus"] != "") &
                      (table["blood_meal"] != "")].copy()
    if described.empty:
        say("  no specimen has both a tick-genus call and blood-meal data")
    else:
        say(f"  specimens with both a genus call and blood-meal data: {len(described)}")
        say("")
        say(pd.crosstab([described["tick_genus"], described["collection_host"]],
                        described["blood_meal"]).to_string())

    say("\n" + "=" * 70)
    say("PATHOGEN CARRIAGE BY TICK SPECIES AND HOST")
    say("=" * 70)
    say("(metabarcoding; qPCR shown separately below because the assays are")
    say(" species-level and the sequencing calls are genus-level)")
    with_calls = table[(table["tick_genus"] != "") & (table["sequenced"] == "yes")].copy()
    if with_calls.empty:
        say("  no specimen has both a genus call and sequencing pathogen calls")
    else:
        for column in SEQ_PATHOGENS.values():
            with_calls[column] = pd.to_numeric(with_calls[column], errors="coerce")
        grouped = with_calls.groupby(["tick_genus", "collection_host"])
        summary = grouped[list(SEQ_PATHOGENS.values())].sum()
        summary["n_specimens"] = grouped.size()
        say(summary.to_string())

    say("\n" + "=" * 70)
    say("qPCR CARRIAGE BY TICK SPECIES AND HOST")
    say("=" * 70)
    q_calls = table[(table["tick_genus"] != "") & (table["qpcr_screened"] == "yes")].copy()
    if q_calls.empty:
        say("  no specimen has both a genus call and qPCR results")
    else:
        for column in QPCR_TARGETS.values():
            q_calls[column] = pd.to_numeric(q_calls[column], errors="coerce")
        grouped = q_calls.groupby(["tick_genus", "collection_host"])
        summary = grouped[list(QPCR_TARGETS.values())].sum()
        summary["n_specimens"] = grouped.size()
        say(summary.to_string())

    say("\n" + "=" * 70)
    say("COVERAGE — what is actually populated")
    say("=" * 70)
    say(f"  specimens in the crosswalk          {len(table)}")
    say(f"  with a tick-genus call              {(table['tick_genus'] != '').sum()}")
    say(f"  sequenced (pathogen calls)          {(table['sequenced'] == 'yes').sum()}")
    say(f"  qPCR screened                       {(table['qpcr_screened'] == 'yes').sum()}")
    say(f"  BOTH methods                        "
        f"{((table['sequenced'] == 'yes') & (table['qpcr_screened'] == 'yes')).sum()}")
    say(f"  genus + both methods                "
        f"{((table['tick_genus'] != '') & (table['sequenced'] == 'yes') & (table['qpcr_screened'] == 'yes')).sum()}")
    say("")
    say("  The last figure is the honest denominator for a full four-way cross-tabulation.")
    say("  Most specimens are missing at least one component, and the table above shows")
    say("  each combination on whatever specimens have the parts it needs, which is why")
    say("  the totals differ between blocks.")

    say("\nREAD BEFORE USING:")
    say("  - a BLANK qPCR cell means the assay was not run, not that it was negative")
    say("  - Borreliella cannot be resolved to species by 16S; a Borreliella-positive")
    say("    tick is not necessarily B. burgdorferi-positive")
    say("  - Cervidae merges deer, moose and unresolved cervid: cytb cannot separate them")
    say("  - the deer ked is a co-detection, not a blood meal, and has its own column")
    say("  - Bartonella here is B. schoenbuchensis, ked-associated, NOT a tick-borne")
    say("    human pathogen; it is in a codetect_ column, outside the pathogen block")
    say("  - human cytb is excluded: the highest human read count in the study is in a")
    say("    positive control, so it is contamination rather than feeding")

    summary_path = os.path.join(os.path.dirname(args.out) or ".",
                                "feeding_ecology_summary.txt")
    with open(summary_path, "w") as handle:
        handle.write("\n".join(report_lines))
    print(f"\nsummary -> {summary_path}")


if __name__ == "__main__":
    main()
