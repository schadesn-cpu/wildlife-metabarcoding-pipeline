#!/usr/bin/env python3
"""
Metadata build + join script
=============================

Joins Adrianna's sample-tracking spreadsheet against the two metadata
templates produced by build_symlink_farm.py, applying comprehensive
canonical-ID normalization so the many naming variants on both sides
resolve correctly.

Inputs (defaults point at the standard project layout on ron):

  - Spreadsheet (TSV or CSV) with columns including at least
    'Barcode' and ideally the biological columns too
  - reads_staged/16S_metadata_template.tsv
  - reads_staged/cytb_metadata_template.tsv

Outputs:

  - reads_staged/16S_metadata.tsv   (joined, ready for QIIME2)
  - reads_staged/cytb_metadata.tsv  (joined, ready for QIIME2)
  - reads_staged/metadata_join_report.csv   audit of matched/unmatched

Normalization rules applied to BOTH sides before joining:

  - Lowercase eq -> uppercase EQ
  - 'ALHT WTD_##'  and  'ALHTDeer-##'  and  'ALHT-DEER-##'   ->  'ALHT-Deer-##'
  - 'TB_#########'                                            ->  'TB-#########'
  - Strip trailing replicate (-1, -2, -3, -4) so ALHT-Deer-04-1 matches
    spreadsheet ALHT-Deer-04. Preserved in the 'subsample_tag' column
    so we don't lose the distinction.
  - Leading-zero normalization (ALHT-Deer-04 ~ ALHT-Deer-4)

Usage:

    python3 build_metadata.py --spreadsheet /path/to/adrianna_sheet.tsv

    # explicit paths
    python3 build_metadata.py \\
        --spreadsheet /path/to/adrianna_sheet.tsv \\
        --stage-dir /home/users/sne24/meedlab/tick_project/reads_staged

No non-standard dependencies.
"""

import re
import csv
import sys
import argparse
from pathlib import Path
from collections import defaultdict, Counter


# ============================================================
# Defaults
# ============================================================
DEFAULT_STAGE_DIR = Path(
    "/home/users/sne24/meedlab/tick_project/reads_staged"
)


# ============================================================
# Extra control-pattern detection
# ============================================================
# The staging-farm control detector catches the canonical control names
# (NTC, POSH, XB-A1, XBG12, etc.) but the genome center delivered some
# additional variants that slipped through. We catch those here at the
# metadata-join step so they're still correctly classified downstream.

EXTRA_CONTROL_PATTERNS = [
    (re.compile(r"^EXTRACTIONBLANK",           re.I), "extraction_blank"),
    (re.compile(r"^POS(-\d+)?$",               re.I), "positive_control"),
    (re.compile(r"^PC-WELL",                   re.I), "positive_control"),
    (re.compile(r"^H\d+PC-WELL",               re.I), "positive_control"),
    (re.compile(r"^XB-OR-NTC",                 re.I), "extraction_blank_or_ntc"),
    (re.compile(r"^H\d+XB-OR-NTC",             re.I), "extraction_blank_or_ntc"),
    (re.compile(r"^XB\d*-(18sV|ITS)",          re.I), "extraction_blank"),
]

# Patterns that indicate samples belonging to a different marker study
# (e.g. 18S/ITS reads that were delivered in a combined 16S demux).
# These get set_aside=yes so they can be excluded from the 16S analysis
# without deleting anything on disk.
SET_ASIDE_PATTERNS = [
    (re.compile(r"-18sV\d+($|-)", re.I),  "18S_marker_wrong_tree"),
    (re.compile(r"-ITS\d?-\d+($|-)", re.I), "ITS_marker_wrong_tree"),
    (re.compile(r"-ITS1-2$",              re.I), "ITS_marker_wrong_tree"),
]


def detect_extra_control(biological_id):
    """Return (is_control, kind) if this ID matches an additional control pattern."""
    for pat, kind in EXTRA_CONTROL_PATTERNS:
        if pat.match(biological_id):
            return (True, kind)
    return (False, "")


def detect_set_aside(biological_id):
    """Return (should_set_aside, reason) for samples from other marker studies."""
    for pat, reason in SET_ASIDE_PATTERNS:
        if pat.search(biological_id):
            return (True, reason)
    return (False, "")


# ============================================================
# Canonical ID normalization
# ============================================================

def canonicalize(raw):
    """Produce a canonical barcode form for matching across naming variants.

    Returns (canonical_id, subsample_tag).
      canonical_id  : normalized ID suitable for joining
      subsample_tag : '' or the stripped replicate suffix ('1', '2', ...)
    """
    if raw is None:
        return ("", "")
    s = str(raw).strip()
    if not s:
        return ("", "")

    # --- Global character-class normalizations ---

    # ALHT WTD_##  ->  ALHT-Deer-##   (spreadsheet form)
    m = re.match(r"^ALHT\s+WTD[_\s]*(\d+)$", s, re.I)
    if m:
        return (f"ALHT-Deer-{int(m.group(1))}", "")

    # ALHTDeer-##  ->  ALHT-Deer-##   (some FASTQ filenames, no hyphen)
    m = re.match(r"^ALHTDeer[-_](\d.*)$", s, re.I)
    if m:
        s = f"ALHT-Deer-{m.group(1)}"

    # ALHT-DEER-##  ->  ALHT-Deer-##  (all-caps variant)
    if "ALHT-DEER" in s.upper() and s.startswith("ALHT-DEER"):
        s = "ALHT-Deer-" + s[10:]

    # ALHT-Deer### (missing hyphen) -> ALHT-Deer-###
    m = re.match(r"^(ALHT-Deer)(\d.*)$", s)
    if m:
        s = f"{m.group(1)}-{m.group(2)}"

    # TB_#########  ->  TB-#########
    if re.match(r"^TB_\d", s):
        s = "TB-" + s[3:]

    # lowercase eq##### -> EQ#####
    if re.match(r"^eq\d", s):
        s = "EQ" + s[2:]

    # USGS-tick-0603 / USGSTICKS-#### / USGSTICKS-LHT* / USGSTICKS-VT*
    # preserve as-is (they align across spreadsheet and FASTQ)

    # --- Replicate-suffix stripping ---
    # Strip a trailing "-N" if:
    #   (a) the root (everything before the final -N) ends in a digit
    #       — e.g. ALHT-Deer-04-1, EQ02441291-2, USGSTICKS-6301-1,
    #              USGSTICKS-LHT002-1 — or
    #   (b) the tail is EXACTLY 2 zero-padded digits (01..09)
    #       — e.g. USGSTICKS-VTBearA-01, USGSTICKS-VT-12-VTW-01
    #
    # The "exactly 2 digits" constraint is important: IDs like
    # USGS-tick-0604 and USGSTICKS-LHT-Deer-091 have zero-padded tails
    # longer than 2 digits, but those are full sample numbers, not
    # subsample tags, and must NOT be stripped. Real subsample tags in
    # this dataset are always 1-2 digits (1, 2, 01, 02, 03, ...).
    subsample = ""
    m = re.match(r"^(.*)-(\d+)$", s)
    if m:
        root = m.group(1)
        tail = m.group(2)
        root_ends_digit = bool(re.search(r"\d$", root))
        tail_is_padded_2digit = len(tail) == 2 and tail.startswith("0")
        if root_ends_digit or tail_is_padded_2digit:
            s = root
            subsample = tail

    # --- Leading-zero normalization on the final numeric segment ---
    m = re.match(r"^(.*-)(\d+)$", s)
    if m:
        prefix, num = m.group(1), m.group(2)
        s = prefix + str(int(num))

    return (s, subsample)


# ============================================================
# IO helpers
# ============================================================

def read_spreadsheet(path):
    """Read the lab spreadsheet. Auto-detects TSV vs CSV.

    Handles the common export quirks:
      - trailing empty rows (Excel sum footer)
      - columns with trailing whitespace in header
      - completely blank rows
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    # Sniff tab vs comma
    first_line = text.split("\n", 1)[0]
    delim = "\t" if first_line.count("\t") > first_line.count(",") else ","

    rows = []
    reader = csv.DictReader(text.splitlines(), delimiter=delim)
    # Normalize header (strip whitespace)
    orig_fieldnames = reader.fieldnames or []
    clean_fieldnames = [h.strip() if h else h for h in orig_fieldnames]
    field_map = dict(zip(orig_fieldnames, clean_fieldnames))

    for r in reader:
        # Re-key with stripped names
        clean_row = {field_map[k]: (v.strip() if isinstance(v, str) else v)
                     for k, v in r.items() if k in field_map}
        # Skip rows with no Barcode
        if not clean_row.get("Barcode"):
            continue
        rows.append(clean_row)

    return rows, clean_fieldnames


def read_metadata_template(path):
    """Read a metadata template TSV from reads_staged/."""
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)
    return rows, fields


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spreadsheet", type=Path, required=True,
                    help="Path to Adrianna's spreadsheet (TSV or CSV)")
    ap.add_argument("--stage-dir", type=Path, default=DEFAULT_STAGE_DIR)
    ap.add_argument("--barcode-col", default="Barcode",
                    help="Column name in spreadsheet containing barcodes")
    args = ap.parse_args()

    if not args.spreadsheet.is_file():
        sys.exit(f"ERROR: spreadsheet not found: {args.spreadsheet}")

    # --- Load spreadsheet ---
    ss_rows, ss_fields = read_spreadsheet(args.spreadsheet)
    print(f"Spreadsheet: {args.spreadsheet}")
    print(f"  {len(ss_rows)} non-blank rows")
    print(f"  Columns: {len(ss_fields)}")
    if args.barcode_col not in ss_fields:
        sys.exit(f"ERROR: barcode column '{args.barcode_col}' not in "
                 f"spreadsheet headers: {ss_fields}")

    # Build canonical lookup: canonical_id -> first matching row
    # (if multiple rows share a canonical id, later ones are flagged)
    ss_canon = {}
    ss_duplicates = defaultdict(list)
    for r in ss_rows:
        canon, _ = canonicalize(r.get(args.barcode_col, ""))
        if not canon:
            continue
        if canon in ss_canon:
            ss_duplicates[canon].append(r[args.barcode_col])
        else:
            ss_canon[canon] = r

    if ss_duplicates:
        print(f"  WARNING: {len(ss_duplicates)} canonical IDs map from >1 "
              f"spreadsheet row; keeping the first. Affected IDs:")
        for canon, extras in list(ss_duplicates.items())[:5]:
            print(f"    {canon}  (also: {', '.join(extras)})")
    print()

    # --- Build which spreadsheet columns we want to merge in ---
    # Drop the join key and keep everything else
    ss_cols_to_merge = [c for c in ss_fields if c != args.barcode_col]

    # --- Process each marker's metadata template ---
    report_rows = []
    for marker in ("16S", "cytb"):
        tpl = args.stage_dir / f"{marker}_metadata_template.tsv"
        out = args.stage_dir / f"{marker}_metadata.tsv"
        if not tpl.is_file():
            print(f"  [skip] template not found: {tpl}")
            continue

        tpl_rows, tpl_fields = read_metadata_template(tpl)

        # Extend the fieldnames: add canonical columns up front, then spreadsheet cols
        new_fields = list(tpl_fields) + [
            "spreadsheet_barcode",      # the raw barcode as typed in the sheet
            "canonical_barcode",        # normalized match key
            "subsample_tag",            # replicate suffix that was stripped (if any)
            "metadata_matched",         # yes/no
            "reclassified_as_control",  # if we caught a control that slipped the staging detector
            "set_aside",                # yes if belongs to a different marker study
            "set_aside_reason",
        ] + [f"sheet__{c}" for c in ss_cols_to_merge]

        matched = 0
        reclassified = 0
        set_aside_count = 0
        for r in tpl_rows:
            bio_id = r.get("biological_sample_id", "")
            canon, subsample = canonicalize(bio_id)
            r["canonical_barcode"] = canon
            r["subsample_tag"] = subsample

            # Catch additional control patterns the staging-step detector missed
            extra_is_ctrl, extra_kind = detect_extra_control(bio_id)
            if extra_is_ctrl and r.get("is_control") != "yes":
                r["is_control"] = "yes"
                r["control_kind"] = extra_kind
                r["reclassified_as_control"] = "yes"
                reclassified += 1
            else:
                r["reclassified_as_control"] = ""

            # Set-aside flag for wrong-marker samples (18S / ITS in 16S tree)
            set_aside, reason = detect_set_aside(bio_id)
            r["set_aside"] = "yes" if set_aside else ""
            r["set_aside_reason"] = reason
            if set_aside:
                set_aside_count += 1

            # Controls don't need spreadsheet metadata; treat them as matched=controls
            if r.get("is_control") == "yes":
                r["spreadsheet_barcode"] = ""
                r["metadata_matched"] = "control"
                for c in ss_cols_to_merge:
                    r[f"sheet__{c}"] = ""
                continue

            # Set-aside samples don't need spreadsheet metadata either
            if set_aside:
                r["spreadsheet_barcode"] = ""
                r["metadata_matched"] = "set_aside"
                for c in ss_cols_to_merge:
                    r[f"sheet__{c}"] = ""
                continue

            hit = ss_canon.get(canon)
            if hit:
                matched += 1
                r["spreadsheet_barcode"] = hit.get(args.barcode_col, "")
                r["metadata_matched"] = "yes"
                for c in ss_cols_to_merge:
                    r[f"sheet__{c}"] = hit.get(c, "")
            else:
                r["spreadsheet_barcode"] = ""
                r["metadata_matched"] = "no"
                for c in ss_cols_to_merge:
                    r[f"sheet__{c}"] = ""

            report_rows.append({
                "marker": marker,
                "qiime_sample_id": r.get("sample-id", ""),
                "biological_sample_id": bio_id,
                "canonical_barcode": canon,
                "matched": r["metadata_matched"],
                "spreadsheet_barcode": r["spreadsheet_barcode"],
            })

        # Write joined output
        with out.open("w") as fh:
            writer = csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(tpl_rows)
        extras = []
        if reclassified:
            extras.append(f"+{reclassified} reclassified as controls")
        if set_aside_count:
            extras.append(f"{set_aside_count} set aside (wrong marker)")
        extras_str = f"  ({', '.join(extras)})" if extras else ""
        print(f"  {marker}:  {matched}/{len(tpl_rows)} samples matched "
              f"({100*matched/max(len(tpl_rows),1):.1f}%)   ->  {out.name}{extras_str}")

    # --- Coverage report ---
    # Which spreadsheet barcodes were never used in a join?
    used_canons = set(r["canonical_barcode"] for r in report_rows
                      if r["matched"] == "yes")
    orphan_sheet_rows = [(canon, r) for canon, r in ss_canon.items()
                         if canon not in used_canons]

    report_path = args.stage_dir / "metadata_join_report.csv"
    with report_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["section", "marker", "qiime_sample_id",
                         "biological_sample_id", "canonical_barcode",
                         "matched", "spreadsheet_barcode"])
        for r in report_rows:
            writer.writerow(["fastq", r["marker"], r["qiime_sample_id"],
                             r["biological_sample_id"], r["canonical_barcode"],
                             r["matched"], r["spreadsheet_barcode"]])
        for canon, r in orphan_sheet_rows:
            writer.writerow(["spreadsheet_orphan", "", "", "", canon, "",
                             r.get(args.barcode_col, "")])

    print()
    print(f"Spreadsheet orphans (rows with metadata but no FASTQ match): "
          f"{len(orphan_sheet_rows)}")
    if orphan_sheet_rows:
        print("  First 10:")
        for canon, r in orphan_sheet_rows[:10]:
            town = r.get("Town", "") or r.get("Town ", "")
            state = r.get("State", "")
            print(f"    {r.get(args.barcode_col, ''):30}  "
                  f"[{town.strip()} {state}]")
    print()
    print(f"Join report written to: {report_path}")


if __name__ == "__main__":
    main()
