#!/usr/bin/env python3
"""Step 4a: describe the design before testing anything.

This is the gate. It builds the three carriage-confidence groups and cross-tabulates them
against batch, host, state and tick genus. Nothing here is a hypothesis test, and it does
not write a community table — the point is to find out whether the comparison in
ANALYSIS_PLAN_pathogen_carriage_vs_microbiome.md is answerable before I spend effort on it.

The specific thing I am looking for: batch explains 10.7% of community variation (p = 0.001)
while collection state explains 0.3% and is not significant. If carriage status is confounded
with batch, then any community difference between carriage groups is partly a batch
difference and no model separates them after the fact. Better to find that out now.

THE THREE GROUPS
  confirmed_carrier    positive by BOTH qPCR and sequencing
  single_method        positive by one only, and the two halves stay labelled because
                       sequencing-only and qPCR-only mean different things
  confirmed_negative   negative by both

  Out-of-scope qPCR targets (Powassan, B. mayonii, B. miyamotoi, E. muris) are NOT
  discordance: 16S and 18S cannot detect a virus or resolve those species, so a qPCR-only
  positive for them is a detection by the only capable method. Those ticks are confirmed
  carriers.

Ixodes only. Dermacentor is a different vector with a different endosymbiont and no panel
pathogens at all, so it is described separately, not pooled in.

Run:
    cd ~/meedlab/tick_project
    python3 scripts/carriage_design_check.py \
        --out results/analysis/carriage_design_check.tsv

Writes:
    results/analysis/carriage_design_check.tsv   one row per tick, with its group
    results/analysis/carriage_design_report.txt  the cross-tabulations
"""

import argparse
import csv
import os
import re
import sys

import pandas as pd

SINGLE_TICK = "results/single_tick_v3/single_tick_v3_full.tsv"
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
WORKBOOK = "metadata/tick_full_metadata.xlsx"
QPCR_SHEET = "qPCR Results"
TICK_GENUS = "results/18S/all/taxonomy_verified/tick_genus_per_sample_corrected.tsv"

# qPCR targets the sequencing markers CAN also see. Only these define discordance.
SHARED_TARGETS = ["a-phagocytophilum", "b-burgdorferi", "b-microti"]

# qPCR targets they CANNOT. Powassan is a virus; the rest are species-level calls below
# what the V4-V5 fragment resolves. A qPCR-only positive here is not a disagreement.
OUT_OF_SCOPE_TARGETS = ["powassan-virus", "b-mayonii", "b-miyamotoi", "e-muris"]

# Sequencing pathogen columns in the single-tick table.
SEQ_PATHOGEN_COLUMNS = [
    "16S_Anaplasma", "16S_Borreliella", "16S_Borrelia", "16S_Bartonella",
    "18S_Babesia_odocoilei", "18S_Babesia_microti",
]

report_lines = []


def say(message=""):
    print(message)
    report_lines.append(str(message))


def normalise(value):
    """Barcodes are written several ways across files; compare on letters and digits only."""
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def is_positive(value):
    return str(value).strip() in ("1", "1.0")


def back_up_existing(out_path):
    """Moves an existing output aside to .bak. Refuses to clobber a backup."""
    if not os.path.exists(out_path):
        return
    backup_path = out_path + ".bak"
    if os.path.exists(backup_path):
        sys.exit(f"{out_path} and {backup_path} both exist. "
                 f"Move or delete {backup_path} before rerunning.")
    os.replace(out_path, backup_path)


def load_tick_genus():
    """barcode -> 'Ixodes' / 'Dermacentor' / 'both' / '' .

    Returns an empty dict with a warning rather than exiting if the file is absent: the
    cross-tabulation is still worth seeing without the genus split, and I would rather
    print a partial answer than nothing.
    """
    if not os.path.exists(TICK_GENUS):
        print(f"warning: {TICK_GENUS} not found; the Ixodes/Dermacentor split will be "
              f"blank and every tick will be treated as unassigned", file=sys.stderr)
        return {}
    # sample_id here is a LIBRARY name carrying a run suffix, e.g.
    # '24-3775.241209_P3-5'. The barcode is everything before the first dot.
    # One specimen can have several libraries, so I collect all their calls and
    # reconcile: disagreement between libraries of the same specimen becomes
    # 'both', which is a real thing in a pooled sample and not an error.
    calls_by_barcode = {}
    with open(TICK_GENUS, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            raw = (row.get("sample_id") or row.get("Barcode") or "").strip()
            if not raw:
                continue
            key = normalise(raw.split(".")[0])
            value = (row.get("dominant") or row.get("tick_genus") or "").strip()
            if value:
                calls_by_barcode.setdefault(key, set()).add(value)

    genus_by_barcode = {}
    for key, values in calls_by_barcode.items():
        tick_calls = {v for v in values if v in ("Ixodes", "Dermacentor", "Amblyomma")}
        if len(tick_calls) == 1:
            genus_by_barcode[key] = tick_calls.pop()
        elif len(tick_calls) > 1:
            genus_by_barcode[key] = "both/mixed"
        else:
            genus_by_barcode[key] = ""
    return genus_by_barcode


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="results/analysis/carriage_design_check.tsv")
    args = parser.parse_args()

    for path in (SINGLE_TICK, CROSSWALK, WORKBOOK):
        if not os.path.isfile(path):
            sys.exit(f"no such file: {path}")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    # ---- the three inputs ---------------------------------------------------
    single_ticks = {normalise(row["Barcode"]): row for row in
                    csv.DictReader(open(SINGLE_TICK, newline=""), delimiter="\t")}
    say(f"single-tick set: {len(single_ticks)} specimens")

    crosswalk = {}
    with open(CROSSWALK, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for key in (row.get("barcode"), row.get("sample_id")):
                if key:
                    crosswalk[normalise(key)] = row

    qpcr = pd.read_excel(WORKBOOK, sheet_name=QPCR_SHEET, dtype=str)
    qpcr["subsamples"] = pd.to_numeric(qpcr["subsample_count"], errors="coerce")
    qpcr = qpcr[qpcr["subsamples"] == 1]
    missing = [t for t in SHARED_TARGETS + OUT_OF_SCOPE_TARGETS if t not in qpcr.columns]
    if missing:
        sys.exit(f"{QPCR_SHEET} is missing expected target columns: {missing}")
    say(f"qPCR single ticks: {len(qpcr)}")

    genus_by_barcode = load_tick_genus()

    # Batch is the factor that explains 10.7% of community variation, so it is worth
    # looking for in more than one place. If the crosswalk does not carry it, the
    # diversity metadata does.
    batch_by_barcode = {}
    for candidate in ("results/diversity_rebuild_20260706_2045/meta_wtd_clean.tsv",
                      "results/diversity_allsamples_20260706_2052/meta_all_recovered.tsv"):
        if not os.path.exists(candidate):
            continue
        with open(candidate, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = normalise(str(list(row.values())[0]).split(".")[0])
                value = (row.get("batch_clean") or row.get("batch") or "").strip()
                if key and value and key not in batch_by_barcode:
                    batch_by_barcode[key] = value
        if batch_by_barcode:
            say(f"batch labels from {candidate}: {len(batch_by_barcode)}")
            break
    if not batch_by_barcode:
        print("warning: no batch labels found in any candidate file", file=sys.stderr)

    # ---- build the groups ---------------------------------------------------
    rows = []
    for _, qpcr_row in qpcr.iterrows():
        barcode = normalise(qpcr_row["host_id"])
        tick = single_ticks.get(barcode)
        if tick is None:
            continue                      # qPCR-only tick, not sequenced

        shared_positive = any(is_positive(qpcr_row.get(t)) for t in SHARED_TARGETS)
        out_of_scope_positive = any(is_positive(qpcr_row.get(t))
                                    for t in OUT_OF_SCOPE_TARGETS)
        qpcr_positive = shared_positive or out_of_scope_positive
        seq_positive = any(float(tick.get(c) or 0) > 0 for c in SEQ_PATHOGEN_COLUMNS)

        if seq_positive and qpcr_positive:
            group, detail = "confirmed_carrier", "both methods"
        elif out_of_scope_positive and not shared_positive and not seq_positive:
            # Only detectable by qPCR. Not a disagreement.
            group, detail = "confirmed_carrier", "qPCR only, out-of-scope target"
        elif seq_positive:
            group, detail = "single_method", "sequencing only"
        elif qpcr_positive:
            group, detail = "single_method", "qPCR only"
        else:
            group, detail = "confirmed_negative", "both negative"

        n_seq = sum(1 for c in SEQ_PATHOGEN_COLUMNS if float(tick.get(c) or 0) > 0)
        n_qpcr = sum(1 for t in SHARED_TARGETS + OUT_OF_SCOPE_TARGETS
                     if is_positive(qpcr_row.get(t)))

        meta = crosswalk.get(barcode, {})
        rows.append({
            "barcode": tick["Barcode"],
            "group": group,
            "detail": detail,
            "tick_genus": genus_by_barcode.get(barcode, ""),
            "host": meta.get("host_bin", ""),
            "state": meta.get("state", ""),
            "batch": (meta.get("batch", "") or meta.get("Batch", "")
                      or batch_by_barcode.get(barcode, "")),
            "n_pathogens_seq": n_seq,
            "n_targets_qpcr": n_qpcr,
        })

    if not rows:
        sys.exit("no ticks matched between the single-tick set and the qPCR sheet")

    table = pd.DataFrame(rows)
    say(f"\nticks screened by BOTH methods: {len(table)}")

    # ---- the groups ---------------------------------------------------------
    say("\n=== carriage confidence groups ===")
    for group in ("confirmed_carrier", "single_method", "confirmed_negative"):
        subset = table[table["group"] == group]
        say(f"  {group:<20} {len(subset):>3}")
        for detail, n in subset["detail"].value_counts().items():
            say(f"      {detail:<32} {n:>3}")

    # ---- Ixodes only, which is the analysis set -----------------------------
    say("\n=== tick genus ===")
    if table["tick_genus"].astype(str).str.strip().eq("").all():
        say("  no genus assignments joined; cannot restrict to Ixodes.")
        say("  CHECK the key column in " + TICK_GENUS)
    else:
        say(table["tick_genus"].replace("", "unassigned").value_counts().to_string())

    # ---- the confound check, which is the point of this script --------------
    for column, label in (("batch", "BATCH — the one that matters"),
                          ("host", "host"),
                          ("state", "state")):
        say(f"\n=== group x {label} ===")
        if table[column].astype(str).str.strip().eq("").all():
            say(f"  '{column}' is empty for every tick; cannot check this confound.")
            continue
        crosstab = pd.crosstab(table[column].replace("", "(blank)"), table["group"])
        say(crosstab.to_string())
        # A batch that is entirely one group cannot be separated from that group.
        for value, counts in crosstab.iterrows():
            nonzero = counts[counts > 0]
            if counts.sum() >= 3 and len(nonzero) == 1:
                say(f"  ** {column}={value} is ENTIRELY {nonzero.index[0]} "
                    f"(n={counts.sum()}) — not separable from group **")

    # ---- is the gradient also a pathogen-load gradient? ---------------------
    say("\n=== pathogens per tick, by group ===")
    say("(if confirmed carriers are mostly the co-infected ticks, the gradient is partly")
    say(" a load gradient and should be described that way)")
    say(table.groupby("group")[["n_pathogens_seq", "n_targets_qpcr"]]
             .agg(["mean", "median", "max"]).round(2).to_string())

    back_up_existing(args.out)
    table.to_csv(args.out, sep="\t", index=False)
    say(f"\nwrote {args.out}  ({len(table)} rows)")

    report_path = os.path.join(os.path.dirname(args.out) or ".",
                               "carriage_design_report.txt")
    with open(report_path, "w") as handle:
        handle.write("\n".join(report_lines))
    print(f"report -> {report_path}")

    say("\nGATE: if any batch is entirely one carriage group, or if the groups are")
    say("badly unbalanced across batches, STOP. Report the confound instead of")
    say("modelling through it.")


if __name__ == "__main__":
    main()
