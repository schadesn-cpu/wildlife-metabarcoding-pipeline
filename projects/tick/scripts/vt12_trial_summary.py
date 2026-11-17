#!/usr/bin/env python3
"""VT-12 trial ticks: what is in each one, and how they differ from each other.

For C. Sullivan. Eleven VT-12 specimens: nine carrying letter codes with no recorded host
(one tick each), and two bear-collected pools.

WHAT THIS IS NOT
No hypothesis test and no p-value. At most three specimens share a letter code, 59TL1 and
59TL2 look like replicates of one animal, and the whole set sits in one sequencing batch,
which explains 10.7% of community variation elsewhere in this project. Treatment cannot be
separated from batch, host or tick count in this design. Everything here is description.

The letter codes (LC, TL, TW, WC, VTW) are treated as opaque labels. Nothing in the project
metadata decodes them, and I would rather print them as they appear than guess.

This is BACTERIAL (16S) and eukaryotic (18S) and blood meal (cytb). It is NOT the fungal
data: ITS is a separate marker and is not in these tables.

Run:
    cd ~/meedlab/tick_project
    python3 scripts/vt12_trial_summary.py --outdir results/analysis/vt12

Writes:
    <outdir>/vt12_taxonomy.tsv        every detection, one row each
    <outdir>/vt12_genus_matrix.tsv    specimen x genus relative abundance, 16S
    <outdir>/vt12_report.txt
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

import pandas as pd

SPECIMEN_TAXON = "results/tables/specimen_taxon_table.tsv"
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"

# The eleven. Written out rather than pattern-matched on "VT-12" so that a specimen added
# later does not silently join the set.
VT12 = [
    "USGSTICKS-VT-12-51LC", "USGSTICKS-VT-12-51TL", "USGSTICKS-VT-12-53TL",
    "USGSTICKS-VT-12-54WC", "USGSTICKS-VT-12-57TW", "USGSTICKS-VT-12-59TL1",
    "USGSTICKS-VT-12-59TL2", "USGSTICKS-VT-12-59TW2", "USGSTICKS-VT-12-VTW",
    "USGSTICKS-VT-12-BearA", "USGSTICKS-VT-12-BearC",
]

report_lines = []


def say(message=""):
    print(message)
    report_lines.append(str(message))


def normalise(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def parse_code(specimen):
    """Animal or site number and letter code from the specimen name.

    'USGSTICKS-VT-12-59TL2' -> ('59', 'TL', '2'). I do not interpret the letters: nothing in
    the project metadata says what LC, TL, TW, WC or VTW mean, and guessing would put a
    wrong label in front of the person who set the trial up.
    """
    tail = specimen.split("VT-12-")[-1]
    match = re.match(r"^(\d+)([A-Za-z]+)(\d*)$", tail)
    if match:
        return match.group(1), match.group(2).upper(), match.group(3)
    return "", tail.upper(), ""


def back_up_existing(out_path):
    if not os.path.exists(out_path):
        return
    backup = out_path + ".bak"
    if os.path.exists(backup):
        sys.exit(f"{out_path} and {backup} both exist. Move or delete {backup} first.")
    os.replace(out_path, backup)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="results/analysis/vt12")
    parser.add_argument("--min-reads", type=int, default=10,
                        help="reads below this are listed but flagged as trace (default 10)")
    args = parser.parse_args()

    for path in (SPECIMEN_TAXON, CROSSWALK):
        if not os.path.isfile(path):
            sys.exit(f"no such file: {path}. Run make_specimen_taxon_table.py first.")
    os.makedirs(args.outdir, exist_ok=True)

    wanted = {normalise(s): s for s in VT12}

    metadata = {}
    with open(CROSSWALK, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for key in (row.get("barcode"), row.get("sample_id")):
                if key and normalise(key) in wanted:
                    metadata[normalise(key)] = row

    say("=== the VT-12 set ===")
    say(f"{'specimen':<26} {'code':<6} {'host':<14} {'ticks':>5}  markers")
    for name in VT12:
        meta = metadata.get(normalise(name), {})
        _, code, replicate = parse_code(name)
        markers = [m for m in ("16S", "18S", "cytb")
                   if str(meta.get(f"name_{m}", "")).strip() not in ("", "nan")]
        say(f"  {name.replace('USGSTICKS-VT-12-',''):<24} {code+replicate:<6} "
            f"{meta.get('host_bin','(none)'):<14} {meta.get('num_ticks','?'):>5}  "
            f"{','.join(markers) if markers else 'NONE'}")

    say("")
    no_data = [n.replace("USGSTICKS-VT-12-", "") for n in VT12
               if not any(str(metadata.get(normalise(n), {}).get(f"name_{m}", "")).strip()
                          not in ("", "nan") for m in ("16S", "18S", "cytb"))]
    if no_data:
        say("")
        say(f"NO DATA: {', '.join(no_data)} has no library on any marker. It is in the")
        say("workbook but never made it through, so the coded set is "
            f"{9 - len(no_data)} specimens with data, not 9.")
    say("")
    say("Nine specimens carry a letter code, have no recorded host and hold one tick each.")
    say("BearA and BearC are bear-collected pools of 16 and 4 ticks and are not part of the")
    say("coded set. 59TL1 and 59TL2 share an animal number and are likely replicates.")
    say("")
    say("Letter codes as they appear, with no interpretation:")
    codes = defaultdict(list)
    for name in VT12:
        number, code, replicate = parse_code(name)
        if number:
            codes[code].append(number + replicate)
    for code in sorted(codes):
        say(f"  {code:<5} n = {len(codes[code])}   animals {', '.join(sorted(codes[code]))}")
    say("")
    say("Nothing in the project metadata decodes these. C. Sullivan set the trial up and is")
    say("the person who knows what they mean.")

    # ---- detections ---------------------------------------------------------
    # The specimen-taxon table carries EVERY cytb ASV, including the 48,650 the Entrez
    # remap rejected, and it displays them under their ORIGINAL lineage because
    # remap_verified_taxon reads EXCLUDED and the label falls back to the lineage leaf.
    # That is how 1,635 fish ASVs (Actinopteri, Gobiidae) appear as blood meals in a tick.
    # They are not contamination we discovered; they are calls the pipeline already threw
    # out. Reporting them would have been reporting the pipeline's rejects as results.
    rejected_cytb = set()
    asv_cytb = "results/tables/asv_table_cytb.tsv"
    if os.path.exists(asv_cytb):
        with open(asv_cytb, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("remap_action", "").strip() in ("unresolved_excluded", "excluded"):
                    leaf = row.get("taxon_lineage", "").split(";")[-1].strip()
                    leaf = leaf.split("__")[-1].strip()
                    if leaf:
                        rejected_cytb.add(leaf)
    else:
        print(f"warning: {asv_cytb} not found; cannot filter rejected cytb calls",
              file=sys.stderr)

    rows, dropped = [], 0
    with open(SPECIMEN_TAXON, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not (normalise(row["specimen"]) in wanted
                    or normalise(row["barcode"]) in wanted):
                continue
            if row["marker"] == "cytb" and row["taxon"].strip() in rejected_cytb:
                dropped += 1
                continue
            rows.append(row)
    if dropped:
        say(f"\nfiltered out {dropped} cytb detection(s) the Entrez remap had already "
            f"rejected")
    if not rows:
        sys.exit("no VT-12 specimen appears in the specimen-taxon table")

    table = pd.DataFrame(rows)
    table["reads"] = pd.to_numeric(table["reads"], errors="coerce").fillna(0).astype(int)
    table["code"] = table["specimen"].map(lambda s: parse_code(s)[1] + parse_code(s)[2])

    out_path = f"{args.outdir}/vt12_taxonomy.tsv"
    back_up_existing(out_path)
    table.to_csv(out_path, sep="\t", index=False)
    say(f"\nwrote {out_path}  ({len(table)} detections, "
        f"{table['specimen'].nunique()} specimens)")

    # ---- what is in each ----------------------------------------------------
    say("\n" + "=" * 74)
    say("WHAT IS IN EACH SPECIMEN")
    say("=" * 74)
    say(f"(top detections per specimen; anything under {args.min_reads} reads marked trace)")
    for specimen in sorted(table["specimen"].unique()):
        subset = table[table["specimen"] == specimen].sort_values("reads", ascending=False)
        short = specimen.replace("USGSTICKS-VT-12-", "")
        say(f"\n--- {short}   ({subset['reads'].sum():,} reads across "
            f"{subset['marker'].nunique()} marker(s)) ---")
        for _, row in subset.head(10).iterrows():
            trace = "  trace" if row["reads"] < args.min_reads else ""
            category = row["category"] or "-"
            say(f"    {row['marker']:5} {row['taxon'][:32]:32} {category[:24]:24} "
                f"{row['reads']:>8,}{trace}")

    # ---- pathogens, which is the first thing anyone will look for -----------
    say("\n" + "=" * 74)
    say("PATHOGEN AND ENDOSYMBIONT DETECTIONS")
    say("=" * 74)
    flagged = table[table["category"].astype(str).str.startswith(
        ("pathogen", "endosymbiont", "ked_associated", "offpanel", "piroplasm"))]
    if flagged.empty:
        say("  none in any VT-12 specimen")
    else:
        pivot = flagged.pivot_table(index="taxon", columns="specimen", values="reads",
                                    aggfunc="sum", fill_value=0)
        pivot.columns = [c.replace("USGSTICKS-VT-12-", "") for c in pivot.columns]
        say(pivot.to_string())

    # ---- 16S composition, side by side --------------------------------------
    say("\n" + "=" * 74)
    say("16S COMPOSITION, SIDE BY SIDE")
    say("=" * 74)
    say("Relative abundance within each specimen's 16S detections. Descriptive only: these")
    say("are single ticks in one batch, and nothing here separates a treatment effect from")
    say("a batch, host or tick-count difference.")
    # 16S depth first. Two specimens have almost no bacterial reads, and a percentage on
    # 73 reads is noise wearing the clothes of a composition. Printing the depths before the
    # table is the only way that caveat travels with the numbers.
    depth_16s = (table[table["marker"] == "16S"].groupby("specimen")["reads"].sum()
                 .sort_values(ascending=False))
    say("")
    say("16S reads per specimen. READ THIS BEFORE THE TABLE BELOW:")
    thin = []
    for specimen, reads in depth_16s.items():
        short = specimen.replace("USGSTICKS-VT-12-", "")
        if reads < 1000:
            thin.append(short)
            say(f"    {short:<8} {int(reads):>9,}   TOO FEW READS TO INTERPRET")
        else:
            say(f"    {short:<8} {int(reads):>9,}")
    if thin:
        say("")
        say(f"    {' and '.join(thin)} have too little 16S to support a composition. Their")
        say("    percentages below are computed from tens or low hundreds of reads and a")
        say("    single ASV can dominate a column. Do not read them as community profiles.")

    bacterial = table[(table["marker"] == "16S") &
                      (~table["category"].isin(["crossamp", "eukaryotic_in_16S",
                                                "mitochondrial", "tick_host",
                                                "vertebrate_host", "human_dna"]))]
    if bacterial.empty:
        say("  no 16S detections")
    else:
        matrix = bacterial.pivot_table(index="taxon", columns="specimen", values="reads",
                                       aggfunc="sum", fill_value=0)
        matrix = 100 * matrix / matrix.sum()
        matrix.columns = [c.replace("USGSTICKS-VT-12-", "") for c in matrix.columns]
        top = matrix.loc[matrix.max(axis=1).sort_values(ascending=False).index].head(20)
        say(top.round(2).to_string())
        matrix_path = f"{args.outdir}/vt12_genus_matrix.tsv"
        back_up_existing(matrix_path)
        matrix.round(4).to_csv(matrix_path, sep="\t")
        say(f"\nwrote {matrix_path}")

    say("\n" + "=" * 74)
    say("READ BEFORE USING")
    say("=" * 74)
    say("  - No statistical test is run and no p-value is reported. At most three specimens")
    say("    share a letter code, 59TL1 and 59TL2 are likely replicates of one animal, and")
    say("    the whole set sits in one sequencing batch. Treatment cannot be separated from")
    say("    batch, host or tick count.")
    say("  - The letter codes are printed as they appear. Nothing in the project metadata")
    say("    decodes them.")
    say("  - BearA and BearC are pools of 16 and 4 ticks with a recorded host. They are not")
    say("    part of the coded set and should not be averaged with it.")
    say("  - This is bacterial, eukaryotic and blood meal. It is NOT the fungal data: ITS is")
    say("    a separate marker and is not in these tables.")
    say("  - These specimens are EXCLUDED from the main analyses as fungal-trial samples,")
    say("    because they are experimental rather than wild-caught.")

    report_path = f"{args.outdir}/vt12_report.txt"
    with open(report_path, "w") as handle:
        handle.write("\n".join(report_lines))
    print(f"\nreport -> {report_path}")


if __name__ == "__main__":
    main()
