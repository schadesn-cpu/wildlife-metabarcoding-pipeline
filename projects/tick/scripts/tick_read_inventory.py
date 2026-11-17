#!/usr/bin/env python3
"""
Tick project READ-ONLY inventory script
=========================================

Walks the Hubbard Center project export directories for 16S and cytb and
produces a single CSV listing every paired-end FASTQ pair, with detected
naming anomalies, control flags, and duplicate flags (within-batch and
cross-batch).

This script does NOT rename, move, or modify any files. It only reads
file names. Nothing in your reads/ directories is touched.

Usage:
    # Scan both trees at their default locations:
    python3 tick_read_inventory.py

    # Override paths:
    python3 tick_read_inventory.py \\
        --projects-16s  /path/to/16S/cobb-data.sr.unh.edu/projects \\
        --projects-cytb /path/to/cytb/cobb-data.sr.unh.edu/projects

    # Only one tree:
    python3 tick_read_inventory.py --skip-cytb
    python3 tick_read_inventory.py --skip-16s

Output:
    A CSV at <output-dir>/tick_read_inventory.csv. Open in Excel and
    review the flagged rows.

No non-standard Python dependencies — runs fine in your qiime2 env or
any Python 3.6+.
"""

import re
import csv
import sys
import argparse
from pathlib import Path
from collections import defaultdict


# ============================================================
# DEFAULTS - can be overridden on the command line
# ============================================================
DEFAULT_16S_PROJECTS = Path(
    "/home/users/sne24/meedlab/tick_project/reads/16S/"
    "cobb-data.sr.unh.edu/projects"
)
DEFAULT_CYTB_PROJECTS = Path(
    "/home/users/sne24/meedlab/tick_project/reads/cytb/"
    "cobb-data.sr.unh.edu/projects"
)
DEFAULT_OUTPUT_DIR = Path(
    "/home/users/sne24/meedlab/tick_project"
)

# Directory-name substrings that should flag a batch as needs_review = yes.
# Empty by default now that Fungal-Batch12 is confirmed to be a combined
# 16S+cytb demux (per Lawrence's wget lists), not an ITS/fungal run.
NEEDS_REVIEW_SUBSTRINGS = []


# ============================================================
# Patterns
# ============================================================

# Control/blank prefixes. First match wins. Applied to cleaned_sample_id
# (i.e. AFTER stripping -16s / -CB suffix).
CONTROL_PATTERNS = [
    (re.compile(r"^NTCH?\d*$", re.I),          "negative_control"),  # NTC, NTCH5, NTCH11
    (re.compile(r"^NTC[-_]", re.I),            "negative_control"),
    (re.compile(r"^POSH?\d", re.I),            "positive_control"),  # POSH6, POSH12
    (re.compile(r"^POS[-_\d]", re.I),          "positive_control"),
    (re.compile(r"^XBG\d", re.I),              "extraction_blank"),  # XBG12
    (re.compile(r"^XB[-_]?[A-Z0-9]*$", re.I),  "extraction_blank"),  # XB, XB1, XB-A9, XB-G6
    (re.compile(r"^BLK", re.I),                "blank"),
    (re.compile(r"^NEG", re.I),                "negative_control"),
]

# Illumina FASTQ filename: SAMPLEID_S##_L00#_R#_001.fastq.gz
FASTQ_RE = re.compile(
    r"^(?P<sample>.+?)_S(?P<snum>\d+)_L(?P<lane>\d+)_R(?P<read>[12])_001\.fastq\.gz$"
)


# ============================================================
# Helpers
# ============================================================

def classify_marker(sample_raw, expected_marker=None):
    """
    Parse the marker suffix off the raw sample name.

    expected_marker : 'cytb' or '16S' or None
        A hint based on which tree the file was found in. Used only to
        resolve ambiguity (e.g. a bare sample name in a cytb-labeled run
        is presumed to be cytb rather than unknown).

    Returns (marker_tag, marker_guess, cleaned_id):
      marker_tag   : literal suffix found, e.g. '-16s', '-16s-N', '-CB',
                     '-cytb', '-cytob', '-cytopb', or ''
      marker_guess : '16S' | 'cytb' | 'unknown_needs_verification'
                     | 'unknown_no_suffix'
      cleaned_id   : sample_raw with the marker suffix stripped
    """
    lower = sample_raw.lower()

    # Pattern like EQ02441291-16s-2 (marker, then replicate number)
    m = re.search(r"-16s-(\d+)$", lower)
    if m:
        # Keep a hyphen before the replicate digits so we don't collapse
        # USGSTICKS-6301-1-cytopB-2 into USGSTICKS-6301-12.
        cleaned = sample_raw[:m.start()] + "-" + m.group(1)
        return ("-16s-N", "16S", cleaned)

    # Similar pattern for cytb/cytoB/cytopB with a trailing replicate number
    m = re.search(r"-(cyt[oa-z]*b)-(\d+)$", lower)
    if m:
        cleaned = sample_raw[:m.start()] + "-" + m.group(2)
        return (f"-{m.group(1)}-N", "cytb", cleaned)

    if lower.endswith("-16s"):
        return ("-16s", "16S", sample_raw[:-4])

    # Accept the range of cytb spellings seen in directory names:
    # -cytb, -cytob, -cytopb  (case-insensitive)
    m = re.search(r"-(cyt[oa-z]*b)$", lower)
    if m:
        return (f"-{m.group(1)}", "cytb", sample_raw[:m.start()])

    if lower.endswith("-cb"):
        # '-CB' is the genome center's older abbreviation for cytoB, confirmed
        # via Lawrence's wget lists: the combined-demux folders (P3-5 241209,
        # Fungal-Batch12) appear on BOTH the 16S and cytb download lists, and
        # the cytoB-suffixed files from those lists are the -CB files in the
        # 16S tree. Treat as cytb regardless of which tree it's in.
        return ("-CB", "cytb", sample_raw[:-3])

    # No marker suffix at all. If we were told the tree, trust it.
    if expected_marker in ("16S", "cytb"):
        return ("", expected_marker, sample_raw)
    return ("", "unknown_no_suffix", sample_raw)


def detect_control(cleaned_id):
    for pattern, kind in CONTROL_PATTERNS:
        if pattern.match(cleaned_id):
            return (True, kind)
    return (False, "")


def detect_anomalies(sample_raw):
    issues = []

    # "ALHT-Deer" glued directly to digits with no hyphen (e.g. ALHT-Deer143-2)
    if re.search(r"ALHT-Deer\d", sample_raw):
        issues.append("missing_hyphen_after_Deer")

    # All-caps DEER variant
    if "ALHT-DEER" in sample_raw:
        issues.append("all_caps_DEER")

    # Lowercase eq prefix (vs standard EQ)
    if re.match(r"^eq\d", sample_raw):
        issues.append("lowercase_eq_prefix")

    return issues


def parse_batch_short(batch_dir_name):
    """Short readable batch label from the long Hubbard directory name.

    Always prepended with the six-digit run date (YYMMDD) to guarantee
    uniqueness — two runs can share a short name (e.g. P3-5 appears as
    both 241209_... and 241219_...) and collapsing them would hide
    resequencing duplicates.
    """
    # Extract the date prefix (first token before the first underscore)
    date_m = re.match(r"(\d{6})_", batch_dir_name)
    date_prefix = date_m.group(1) if date_m else "??????"

    # Strip the trailing marker suffix so plate-pattern matching doesn't
    # accidentally capture '16s' / 'cytoB' / 'cytopB' as part of a plate
    # designation (e.g. 'P7-16s' -> 'P7' not 'P7-16').
    stripped = re.sub(r"-(?:16[sS]|cyt[oa-z]*B)$", "", batch_dir_name)

    # Figure out the semantic label (most specific pattern wins)
    m = re.search(r"[Bb]atch[-_]?(\d+)", stripped)
    if m:
        label = f"Batch-{int(m.group(1)):02d}"
    elif (m := re.search(r"-P(\d+)-(\d+)", stripped)):
        label = f"P{m.group(1)}-{m.group(2)}"
    elif (m := re.search(r"-P(\d+)", stripped)):
        label = f"P{m.group(1)}"
    else:
        # Fallback: strip the date prefix + Illumina machine/run tokens
        # (A#####_####_XXXXXXDRX#_) and any leading AMP-###### prefix,
        # keeping the human-descriptive tail. E.g.
        # '240430_A01346_0132_BHNWCWDRX3_AMP-043024-DN-USGSTICKS-2023-2024'
        # -> 'USGSTICKS-2023-2024'
        tail = batch_dir_name
        tail = re.sub(r"^\d{6}_", "", tail)                     # strip date
        tail = re.sub(r"^A\d+_\d+_[A-Z0-9]+_", "", tail)        # strip machine+run+flowcell
        tail = re.sub(r"^AMP-\d+-", "", tail)                   # strip AMP-YYMMDD-
        tail = re.sub(r"^(DN|VDL)-", "", tail)                  # strip submitter prefix
        label = tail[:30] if tail else batch_dir_name[:20]

    return f"{date_prefix}_{label}"


# ============================================================
# Main
# ============================================================

def scan_tree(projects_dir, expected_marker, sample_locations, rows_all):
    """
    Scan one projects/ tree. Appends rows to rows_all in place and
    updates cross-tree sample_locations in place.
    """
    if not projects_dir.is_dir():
        print(f"  [skip tree] {projects_dir} does not exist")
        return

    batch_dirs = sorted(d for d in projects_dir.iterdir() if d.is_dir())
    if not batch_dirs:
        print(f"  [skip tree] no batch dirs under {projects_dir}")
        return

    print(f"\nTree: {projects_dir}")
    print(f"  expected marker: {expected_marker}")
    print(f"  batches found:   {len(batch_dirs)}")

    for batch_dir in batch_dirs:
        reads_dir = batch_dir / "reads"
        if not reads_dir.is_dir():
            print(f"    [skip]   {batch_dir.name}: no reads/ subdirectory")
            continue

        batch_short = parse_batch_short(batch_dir.name)
        needs_review = any(s in batch_dir.name for s in NEEDS_REVIEW_SUBSTRINGS)

        pairs = defaultdict(dict)
        for f in sorted(reads_dir.iterdir()):
            if not f.is_file():
                continue
            m = FASTQ_RE.match(f.name)
            if not m:
                continue
            key = (m.group("sample"), m.group("snum"), m.group("lane"))
            pairs[key][m.group("read")] = f

        fastq_count = sum(
            1 for f in reads_dir.iterdir()
            if f.is_file() and f.name.endswith(".fastq.gz")
        )
        print(f"    [{batch_short:12}] {fastq_count:>4} FASTQ files, "
              f"{len(pairs):>4} pairs  ({batch_dir.name})")

        within_batch_ids = defaultdict(list)
        rows_this_batch = []

        for (sample_raw, snum, lane), reads in pairs.items():
            r1 = reads.get("1")
            r2 = reads.get("2")

            marker_tag, marker_guess, cleaned_id = classify_marker(
                sample_raw, expected_marker=expected_marker
            )
            is_control, control_kind = detect_control(cleaned_id)
            anomalies = detect_anomalies(sample_raw)

            # Key cross-tree dupes by (cleaned_id, marker_guess) so that a
            # sample appearing in both the 16S tree and cytb tree is NOT
            # treated as a duplicate — different marker = different library.
            if not is_control:
                dup_key = (cleaned_id.lower(), marker_guess)
                sample_locations[dup_key].append((batch_short, expected_marker))
                within_batch_ids[(cleaned_id.lower(), marker_tag)].append(sample_raw)

            row = {
                "tree": expected_marker or "unknown",
                "batch_short": batch_short,
                "batch_dir": batch_dir.name,
                "original_sample_raw": sample_raw,
                "cleaned_sample_id": cleaned_id,
                "marker_tag": marker_tag,
                "marker_guess": marker_guess,
                "s_number": snum,
                "lane": lane,
                "r1_path": str(r1) if r1 else "",
                "r2_path": str(r2) if r2 else "",
                "paired_complete": "yes" if (r1 and r2) else "NO",
                "is_control": "yes" if is_control else "",
                "control_kind": control_kind,
                "needs_review": "yes" if needs_review else "",
                "anomalies": ";".join(anomalies),
                "within_batch_duplicate": "",
                "cross_batch_duplicate": "",
                "notes": "",
            }
            rows_all.append(row)
            rows_this_batch.append(row)

        for (cid_lower, m_tag), raws in within_batch_ids.items():
            if len(raws) > 1:
                for row in rows_this_batch:
                    if (row["cleaned_sample_id"].lower() == cid_lower
                            and row["marker_tag"] == m_tag
                            and not row["is_control"]):
                        row["within_batch_duplicate"] = f"yes ({len(raws)} entries)"


def main():
    parser = argparse.ArgumentParser(
        description="Read-only inventory of tick project FASTQ files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--projects-16s", type=Path, default=DEFAULT_16S_PROJECTS,
        help=f"Path to 16S projects/ directory (default: {DEFAULT_16S_PROJECTS})",
    )
    parser.add_argument(
        "--projects-cytb", type=Path, default=DEFAULT_CYTB_PROJECTS,
        help=f"Path to cytb projects/ directory (default: {DEFAULT_CYTB_PROJECTS})",
    )
    parser.add_argument(
        "--skip-16s", action="store_true", help="Skip the 16S tree")
    parser.add_argument(
        "--skip-cytb", action="store_true", help="Skip the cytb tree")
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help=f"Where to write the CSV (default: {DEFAULT_OUTPUT_DIR})",
    )
    args = parser.parse_args()

    if args.skip_16s and args.skip_cytb:
        sys.exit("ERROR: both --skip-16s and --skip-cytb given, nothing to do")

    rows = []
    # Key: (cleaned_id_lower, marker_guess) -> list of (batch_short, tree_marker)
    sample_locations = defaultdict(list)

    if not args.skip_16s:
        scan_tree(args.projects_16s, "16S", sample_locations, rows)
    if not args.skip_cytb:
        scan_tree(args.projects_cytb, "cytb", sample_locations, rows)

    if not rows:
        sys.exit("ERROR: no FASTQ pairs found in any tree")

    # Cross-batch duplicates: same (cleaned_id, marker_guess) in >1 batch
    cross_batch = {
        key: sorted(set(b for b, _ in batches))
        for key, batches in sample_locations.items()
        if len(set(b for b, _ in batches)) > 1
    }
    for row in rows:
        key = (row["cleaned_sample_id"].lower(), row["marker_guess"])
        if key in cross_batch:
            row["cross_batch_duplicate"] = "yes (" + ",".join(cross_batch[key]) + ")"

    # Write CSV
    output_csv = args.output_dir / "tick_read_inventory.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "tree", "batch_short", "batch_dir",
        "original_sample_raw", "cleaned_sample_id",
        "marker_tag", "marker_guess",
        "s_number", "lane",
        "r1_path", "r2_path", "paired_complete",
        "is_control", "control_kind",
        "needs_review",
        "anomalies",
        "within_batch_duplicate", "cross_batch_duplicate",
        "notes",
    ]
    with output_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # ============================================================
    # Summary
    # ============================================================
    print("\n" + "=" * 64)
    print("SUMMARY")
    print("=" * 64)
    print(f"Total FASTQ pairs:       {len(rows)}")
    print(f"Output CSV:              {output_csv}")
    print()

    by_marker = defaultdict(int)
    for r in rows:
        by_marker[r["marker_guess"]] += 1
    print("By marker guess:")
    for k in sorted(by_marker):
        print(f"  {k:35} {by_marker[k]}")
    print()

    by_tree_batch_marker = defaultdict(lambda: defaultdict(int))
    for r in rows:
        key = f"{r['tree']} / {r['batch_short']}"
        by_tree_batch_marker[key][r["marker_guess"]] += 1
    print("By tree / batch / marker:")
    for tb in sorted(by_tree_batch_marker):
        parts = ", ".join(f"{k}={v}" for k, v in sorted(by_tree_batch_marker[tb].items()))
        print(f"  {tb:30} {parts}")
    print()

    print(f"Controls:                {sum(1 for r in rows if r['is_control'])}")
    print(f"Needs-review pairs:      {sum(1 for r in rows if r['needs_review'])}")
    print(f"Rows with anomalies:     {sum(1 for r in rows if r['anomalies'])}")
    print(f"Within-batch duplicates: {sum(1 for r in rows if r['within_batch_duplicate'])}")
    print(f"Cross-batch dup keys:    {len(cross_batch)} "
          f"(affecting {sum(1 for r in rows if r['cross_batch_duplicate'])} rows)")
    print(f"Unpaired (R1 or R2 missing): {sum(1 for r in rows if r['paired_complete'] == 'NO')}")

    print()
    print("Next: open the CSV, sort/filter by the flag columns, and eyeball")
    print("anything unexpected. No files were modified.")


if __name__ == "__main__":
    main()
