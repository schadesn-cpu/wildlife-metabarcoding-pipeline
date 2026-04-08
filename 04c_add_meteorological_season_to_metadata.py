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
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


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
    """
    Parse M/D/YY or M/D/YYYY date string and return the month as an integer.

    Returns None for empty strings (expected — sample has no collection date).
    Raises ValueError for non-empty strings that don't match the expected
    format, so the caller can log which sample had the bad date and why.

    Args:
        date_str: Raw date string from the metadata TSV cell.

    Returns:
        Integer month (1–12), or None if date_str is empty.

    Raises:
        ValueError: If date_str is non-empty but cannot be parsed as M/D/YY[YY].
    """
    date_str = date_str.strip()

    # Empty string is the expected case for samples with no collection date.
    # Return None quietly — the caller will assign an empty Season value.
    if not date_str:
        return None

    # Expected format is M/D/YY or M/D/YYYY. Split on "/" and parse the month
    # (first field). Raise ValueError if the format doesn't match so the caller
    # can log which sample caused the problem.
    parts = date_str.split("/")
    if len(parts) < 2:
        raise ValueError(
            f"Expected M/D/YY format but got '{date_str}' "
            f"(no '/' separator found)"
        )
    try:
        return int(parts[0])
    except ValueError:
        raise ValueError(
            f"Could not parse month from '{date_str}' "
            f"(first field '{parts[0]}' is not an integer)"
        )


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
        log.error("Empty file: %s", in_path)
        sys.exit(1)

    # Parse header — handle #q2:types row if present
    header_line = lines[0]
    headers = header_line.split("\t")

    # Find date column
    try:
        date_idx = headers.index(args.date_col)
    except ValueError:
        log.error(
            "Column '%s' not found in header. Available columns: %s",
            args.date_col, headers,
        )
        sys.exit(1)

    # Check if season column already exists
    if args.season_col in headers:
        log.warning("Column '%s' already exists — it will be overwritten.", args.season_col)
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
        # cols[0] is the sample-id — include it in any warnings so it's clear
        # which sample had a date that could not be parsed.
        date_val = cols[date_idx].strip() if date_idx < len(cols) else ""
        try:
            month = parse_date(date_val)
        except ValueError as e:
            log.warning(
                "Row %d (sample '%s'): %s — assigning empty Season.",
                i, cols[0], e,
            )
            month = None

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
