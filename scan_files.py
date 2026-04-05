#!/usr/bin/env python3
"""
scan_files.py
=============
Scan a directory tree for output/results files and save a metadata inventory
to CSV, sorted by most-recently-modified first.

Each row records: directory, filename, extension, size (KB), metadata-change
time, last-modified time, and full path.

By default only common analysis output extensions are included (see
OUTPUT_EXTENSIONS). Pass --all to include every file type.

Usage
-----
    # Scan current directory, write file_inventory.csv
    python scan_files.py

    # Scan a specific results directory
    python scan_files.py results/

    # Scan all file types and write to a custom CSV
    python scan_files.py results/ --all -o full_inventory.csv

Notes
-----
- Hidden files and directories (names starting with '.') are always skipped.
- On Linux, true file creation time is not stored by most filesystems.
  The 'created_or_metadata_change' column contains st_ctime, which records
  the last time the file's metadata (permissions, ownership, etc.) changed —
  not necessarily when the file was first created.
- Symlinks are followed by os.walk; the reported size is the target file size.
"""

import os
import csv
import argparse
import logging
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# Extensions considered "output/results" files for the default scan mode.
# Edit this set to add or remove types for your project. Pass --all at the
# command line to bypass this filter entirely.
OUTPUT_EXTENSIONS = {
    ".csv", ".tsv", ".xlsx", ".xls",
    ".json", ".jsonl",
    ".txt", ".log",
    ".png", ".jpg", ".jpeg", ".pdf", ".svg",
    ".html", ".htm",
    ".parquet", ".feather", ".hdf5", ".h5",
    ".npy", ".npz",
    ".pkl", ".pickle",
    ".out", ".results", ".summary",
}

def get_file_times(filepath: str):
    """
    Return (ctime_str, mtime_str, size_kb) for a file.

    ctime_str: st_ctime formatted as 'YYYY-MM-DD HH:MM:SS'. On Linux this is
               the metadata-change time, not true creation time.
    mtime_str: st_mtime formatted as 'YYYY-MM-DD HH:MM:SS'. This is the last
               time the file's content was modified.
    size_kb:   File size in kilobytes, rounded to two decimal places.

    Raises OSError if the file cannot be stat'd (caller should handle).
    """
    stat = os.stat(filepath)
    modified = datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S")
    # On Linux, st_ctime is the metadata-change time, not true creation time.
    # It's the closest available proxy; stored as-is.
    ctime = datetime.fromtimestamp(stat.st_ctime).strftime("%Y-%m-%d %H:%M:%S")
    size_kb = round(stat.st_size / 1024, 2)
    return ctime, modified, size_kb

def scan_directory(root_dir: str, all_files: bool = False) -> list:
    """
    Recursively walk root_dir and collect metadata for matching files.

    Hidden directories and files (names starting with '.') are always skipped.
    If all_files is False, only files whose extension is in OUTPUT_EXTENSIONS
    are included. Results are sorted by last_modified descending (newest first).

    Args:
        root_dir:  Absolute path to the directory to scan.
        all_files: If True, include every file regardless of extension.

    Returns:
        List of dicts with keys: directory, filename, extension, size_kb,
        created_or_metadata_change, last_modified, full_path.
        Files that cannot be stat'd (permission errors, race conditions) are
        silently skipped.
    """
    rows = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # Skip hidden directories
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        for fname in filenames:
            if fname.startswith("."):
                continue
            ext = os.path.splitext(fname)[1].lower()
            if all_files or ext in OUTPUT_EXTENSIONS:
                full_path = os.path.join(dirpath, fname)
                try:
                    ctime, mtime, size_kb = get_file_times(full_path)
                    rows.append({
                        "directory": dirpath,
                        "filename": fname,
                        "extension": ext,
                        "size_kb": size_kb,
                        "created_or_metadata_change": ctime,
                        "last_modified": mtime,
                        "full_path": full_path,
                    })
                except PermissionError:
                    # File exists but cannot be read — common on shared HPC
                    # filesystems where some directories have restricted access.
                    # Log at WARNING so it's visible without cluttering output.
                    log.warning("Permission denied — skipping: %s", full_path)
                except FileNotFoundError:
                    # File disappeared between os.walk listing and os.stat call.
                    # This is a harmless race condition (e.g. a temp file deleted
                    # by another process). Log at DEBUG only — not actionable.
                    log.debug("File vanished before stat — skipping: %s", full_path)
    # Sort by last_modified descending (newest first)
    rows.sort(key=lambda x: x["last_modified"], reverse=True)
    return rows

def main():
    """
    Parse command-line arguments, run the directory scan, and write the CSV.

    Prints a summary to stdout including total file count, output path, and
    the five most recently modified files found.
    """
    parser = argparse.ArgumentParser(
        description="Scan a directory for output/results files and save to CSV."
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="Root directory to scan (default: current directory)",
    )
    parser.add_argument(
        "-o", "--output",
        default="file_inventory.csv",
        help="Output CSV filename (default: file_inventory.csv)",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Include ALL file types, not just output/results extensions",
    )
    args = parser.parse_args()

    root = os.path.abspath(args.directory)
    print(f"Scanning: {root}")
    print(f"Mode: {'All files' if args.all else 'Output/results files only'}")

    rows = scan_directory(root, all_files=args.all)

    if not rows:
        print("No matching files found.")
        return

    output_path = os.path.abspath(args.output)
    fieldnames = ["directory", "filename", "extension", "size_kb",
                  "created_or_metadata_change", "last_modified", "full_path"]

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n✅ Done! Found {len(rows)} files.")
    print(f"📄 Saved to: {output_path}")
    print(f"\nTop 5 most recently modified:")
    for r in rows[:5]:
        print(f"  {r['last_modified']}  {r['filename']}  ({r['directory']})")

if __name__ == "__main__":
    main()
