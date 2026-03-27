#!/usr/bin/env python3
"""
add_season_to_metadata.py

Adds a 'Season' column to a QIIME2 metadata TSV based on the 'Date Found' column.

Season assignments (meteorological):
  Winter : Dec, Jan, Feb  (months 12, 1, 2)
  Spring : Mar, Apr, May  (months 3, 4, 5)
  Summer : Jun, Jul, Aug  (months 6, 7, 8)
  Fall   : Sep, Oct, Nov  (months 9, 10, 11)

Samples with missing/unknown dates get Season = '' (empty),
which QIIME2 will skip during group-significance tests.

Usage:
    python add_season_to_metadata.py \
        --input  metadata/qiime/metadata_MiFish.tsv \
        --output metadata/qiime/metadata_MiFish.tsv \
        --date-col "Date Found" \
        --season-col Season
"""

import argparse
import sys
from pathlib import Path


def month_to_season(month: int) -> str:
    if month in (12, 1, 2):
        return "Winter"
    elif month in (3, 4, 5):
        return "Spring"
    elif month in (6, 7, 8):
        return "Summer"
    elif month in (9, 10, 11):
        return "Fall"
    return ""


def parse_date(date_str: str):
    """Parse M/D/YY or M/D/YYYY. Returns month int or None."""
    date_str = date_str.strip()
    if not date_str:
        return None
    try:
        parts = date_str.split("/")
        if len(parts) >= 1:
            return int(parts[0])
    except (ValueError, IndexError):
        pass
    return None


def main():
    parser = argparse.ArgumentParser(description="Add Season column to QIIME2 metadata TSV.")
    parser.add_argument("--input",      required=True, help="Input metadata TSV path.")
    parser.add_argument("--output",     required=True, help="Output metadata TSV path (can be same as input).")
    parser.add_argument("--date-col",   default="Date Found", help="Column name containing dates.")
    parser.add_argument("--season-col", default="Season", help="Name for the new season column.")
    args = parser.parse_args()

    in_path  = Path(args.input)
    out_path = Path(args.output)

    lines = in_path.read_text(encoding="utf-8").splitlines()
    if not lines:
        print("[ERROR] Empty file.", file=sys.stderr)
        sys.exit(1)

    # Parse header — handle #q2:types row if present
    header_line = lines[0]
    headers = header_line.split("\t")

    # Find date column
    try:
        date_idx = headers.index(args.date_col)
    except ValueError:
        print(f"[ERROR] Column '{args.date_col}' not found in header.", file=sys.stderr)
        print(f"  Available columns: {headers}", file=sys.stderr)
        sys.exit(1)

    # Check if season column already exists
    if args.season_col in headers:
        print(f"[WARN] Column '{args.season_col}' already exists — it will be overwritten.")
        season_idx = headers.index(args.season_col)
        replace_mode = True
    else:
        season_idx = None
        replace_mode = False

    counts = {"Winter": 0, "Spring": 0, "Summer": 0, "Fall": 0, "unknown": 0}
    out_lines = []

    for i, line in enumerate(lines):
        cols = line.split("\t")

        # Header row
        if i == 0:
            if replace_mode:
                out_lines.append(line)
            else:
                cols.append(args.season_col)
                out_lines.append("\t".join(cols))
            continue

        # #q2:types row
        if cols[0].startswith("#q2:types"):
            if replace_mode:
                out_lines.append(line)
            else:
                cols.append("categorical")
                out_lines.append("\t".join(cols))
            continue

        # Data rows
        date_val = cols[date_idx].strip() if date_idx < len(cols) else ""
        month = parse_date(date_val)

        if month is not None:
            season = month_to_season(month)
        else:
            season = ""
            counts["unknown"] += 1

        if season:
            counts[season] += 1

        if replace_mode:
            while len(cols) <= season_idx:
                cols.append("")
            cols[season_idx] = season
            out_lines.append("\t".join(cols))
        else:
            cols.append(season)
            out_lines.append("\t".join(cols))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")

    print(f"✅ Written to: {out_path}")
    print(f"\nSeason counts:")
    for season in ["Summer", "Fall", "Winter", "Spring"]:
        print(f"  {season:8s}: {counts[season]}")
    print(f"  {'unknown':8s}: {counts['unknown']} (empty date — will be skipped in QIIME2)")


if __name__ == "__main__":
    main()
