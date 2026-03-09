#!/usr/bin/env python3
"""
Scan a directory for output/results files and save metadata to CSV.
Lists: directory, filename, size, created time, modified time.
"""

import os
import csv
import argparse
from datetime import datetime

# Common output/results file extensions
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

def get_file_times(filepath):
    stat = os.stat(filepath)
    modified = datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S")
    # Linux doesn't reliably store true "created" time; use the earlier of mtime/ctime
    ctime = datetime.fromtimestamp(stat.st_ctime).strftime("%Y-%m-%d %H:%M:%S")
    size_kb = round(stat.st_size / 1024, 2)
    return ctime, modified, size_kb

def scan_directory(root_dir, all_files=False):
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
                except (PermissionError, FileNotFoundError):
                    pass
    # Sort by last_modified descending (newest first)
    rows.sort(key=lambda x: x["last_modified"], reverse=True)
    return rows

def main():
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
