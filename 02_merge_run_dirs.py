#!/usr/bin/env python3
"""
00_merge_run_dirs.py
====================
Merge multiple sequencing run directories (each containing flat FASTQ files)
into a single reads/ directory tree ready for 01_make_manifests.py.

Expected input layout (one or more run dirs, FASTQs flat inside):
    run1/
        TV230084-GI-16S_S1492_L002_R1_001.fastq.gz
        TV230084-GI-16S_S1492_L002_R2_001.fastq.gz
        TV230084-GI-cytb_S10_L002_R1_001.fastq.gz
        TV230084-GI-cytb_S10_L002_R2_001.fastq.gz
    run2/
        TV240010-GI-16S_S5_R1_001.fastq.gz
        ...

Output layout (ready for 01_make_manifests.py):
    reads/
    ├── 16S/
    │   ├── TV230084-GI-16S_S1492_L002_R1_001.fastq.gz  (symlink → run1/...)
    │   └── ...
    └── cytb/
        ├── TV230084-GI-cytb_S10_L002_R1_001.fastq.gz   (symlink → run1/...)
        └── ...

Marker detection:
    Files are assigned to a marker directory based on a case-insensitive token
    search in the filename.  Default tokens: 16S, cytb.
    Override with --markers.

Collision handling:
    If the same filename appears in two run dirs, the script aborts by default.
    Use --on-collision warn  to keep the first copy and log a warning.
    Use --on-collision skip  to silently skip duplicates.
    (Collisions are also written to <outdir>/collision_report.tsv)

Usage examples:

    # Dry run — see what would be linked, nothing written:
    python 00_merge_run_dirs.py \\
        --run-dirs run1/ run2/ run3/ \\
        --out-reads reads/ \\
        --markers 16S cytb \\
        --dry-run

    # Real run with symlinks (default, saves disk):
    python 00_merge_run_dirs.py \\
        --run-dirs /data/run1 /data/run2 \\
        --out-reads reads/ \\
        --markers 16S cytb

    # Copy files instead of symlinking (safer for long-term archiving):
    python 00_merge_run_dirs.py \\
        --run-dirs run1/ run2/ \\
        --out-reads reads/ \\
        --copy

    # Then build manifests:
    python 01_make_manifests.py --reads-dir reads/ --outdir qiime2/imported/

Requirements:
    Python >= 3.8  (standard library only)
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
FASTQ_SUFFIXES = {".fastq", ".fastq.gz", ".fq", ".fq.gz"}

DEFAULT_MARKERS = ["16S", "cytb"]


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
class FileEntry(NamedTuple):
    filename: str       # original filename (basename)
    src: Path           # absolute resolved source path
    run_dir: Path       # which run directory it came from
    marker: str         # assigned marker token


class Collision(NamedTuple):
    filename: str
    first_run: Path
    second_run: Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def is_fastq(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(ext) for ext in FASTQ_SUFFIXES)


def detect_marker(filename: str, markers: List[str]) -> Optional[str]:
    """
    Return the first marker token found in *filename* (case-insensitive).
    Matches whole tokens delimited by [-_.] or start/end of string so that
    e.g. '16S' does not accidentally match inside a longer token.
    Returns None if no marker matches.
    """
    name_lower = filename.lower()
    for marker in markers:
        pattern = rf"(?<![a-zA-Z0-9]){re.escape(marker.lower())}(?![a-zA-Z0-9])"
        if re.search(pattern, name_lower):
            return marker
    return None


def collect_fastqs(run_dir: Path, markers: List[str]) -> Tuple[List[FileEntry], List[str]]:
    """
    Scan a single flat run directory.
    Returns (entries, warnings).
    """
    if not run_dir.is_dir():
        return [], [f"Run directory does not exist or is not a directory: {run_dir}"]

    entries: List[FileEntry] = []
    warnings: List[str] = []
    skipped_no_marker: List[str] = []
    skipped_not_fastq: int = 0

    for f in sorted(run_dir.iterdir()):
        if not (f.is_file() or f.is_symlink()):
            continue
        if not is_fastq(f):
            skipped_not_fastq += 1
            continue

        marker = detect_marker(f.name, markers)
        if marker is None:
            skipped_no_marker.append(f.name)
            continue

        try:
            src = f.resolve(strict=True)
        except (OSError, FileNotFoundError):
            warnings.append(f"Broken symlink or unresolvable path — skipping: {f}")
            continue

        entries.append(FileEntry(
            filename=f.name,
            src=src,
            run_dir=run_dir,
            marker=marker,
        ))

    if skipped_no_marker:
        warnings.append(
            f"{run_dir.name}: {len(skipped_no_marker)} FASTQ(s) skipped — "
            f"no marker token found. First few: {skipped_no_marker[:5]}"
        )
    if skipped_not_fastq:
        log.debug("%s: %d non-FASTQ files ignored", run_dir.name, skipped_not_fastq)

    return entries, warnings


def detect_collisions(
    all_entries: List[FileEntry],
) -> Tuple[Dict[str, FileEntry], List[Collision]]:
    """
    Build a filename→entry map; collect collisions where the same filename
    appears in more than one run directory.
    Returns (unique_map, collisions).
    """
    unique: Dict[str, FileEntry] = {}
    collisions: List[Collision] = []

    for entry in all_entries:
        if entry.filename in unique:
            existing = unique[entry.filename]
            if existing.run_dir != entry.run_dir:
                collisions.append(Collision(
                    filename=entry.filename,
                    first_run=existing.run_dir,
                    second_run=entry.run_dir,
                ))
                # keep the first occurrence
        else:
            unique[entry.filename] = entry

    return unique, collisions


def write_collision_report(collisions: List[Collision], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["filename", "first_run_dir", "second_run_dir"])
        for c in collisions:
            writer.writerow([c.filename, str(c.first_run), str(c.second_run)])
    log.info("Collision report: %s", out_path)


# ---------------------------------------------------------------------------
# Main merge logic
# ---------------------------------------------------------------------------

def merge(
    run_dirs: List[Path],
    out_reads: Path,
    markers: List[str],
    on_collision: str,      # "abort" | "warn" | "skip"
    copy: bool,
    dry_run: bool,
) -> int:
    """
    Returns exit code (0 = success, 1 = error).
    """
    # ── Collect all FASTQs from all run dirs ──────────────────────────────
    all_entries: List[FileEntry] = []
    had_warnings = False

    for run_dir in run_dirs:
        entries, warnings = collect_fastqs(run_dir, markers)
        for w in warnings:
            log.warning(w)
            had_warnings = True
        log.info("  %s: %d FASTQ files found for markers %s",
                 run_dir.name, len(entries), markers)
        all_entries.extend(entries)

    if not all_entries:
        log.error("No FASTQ files matched any marker token across all run dirs.")
        log.error("Check --markers and ensure filenames contain the marker token.")
        return 1

    # ── Collision detection ───────────────────────────────────────────────
    unique_map, collisions = detect_collisions(all_entries)

    if collisions:
        collision_report = out_reads / "collision_report.tsv"
        log.warning("─" * 60)
        log.warning("COLLISIONS DETECTED: %d filename(s) appear in >1 run dir", len(collisions))
        for c in collisions[:10]:
            log.warning("  %s", c.filename)
            log.warning("    first : %s", c.first_run)
            log.warning("    second: %s", c.second_run)
        if len(collisions) > 10:
            log.warning("  ... and %d more (see collision report)", len(collisions) - 10)
        log.warning("─" * 60)

        if not dry_run:
            write_collision_report(collisions, collision_report)

        if on_collision == "abort":
            log.error(
                "Aborting due to collisions. Use --on-collision warn|skip to proceed, "
                "or resolve the duplicate files manually.\n"
                "Collision report: %s", collision_report
            )
            return 1
        elif on_collision == "warn":
            log.warning("Proceeding with first occurrence of each duplicate (--on-collision warn).")
        else:
            log.info("Skipping duplicates silently (--on-collision skip).")

    # ── Group by marker ───────────────────────────────────────────────────
    by_marker: Dict[str, List[FileEntry]] = {m: [] for m in markers}
    unassigned = 0
    for entry in unique_map.values():
        if entry.marker in by_marker:
            by_marker[entry.marker].append(entry)
        else:
            unassigned += 1

    # ── Summary before writing ────────────────────────────────────────────
    log.info("=" * 60)
    log.info("Merge plan:")
    total = 0
    for marker in markers:
        n = len(by_marker[marker])
        log.info("  %-10s  %d files → %s/%s/", marker, n, out_reads, marker)
        total += n
    log.info("  Total: %d files", total)
    if collisions:
        log.info("  Collisions: %d (first copy kept)", len(collisions))
    if dry_run:
        log.info("DRY RUN — no files or directories will be created")
    log.info("=" * 60)

    if dry_run:
        # Print full plan
        for marker in markers:
            entries = sorted(by_marker[marker], key=lambda e: e.filename)
            print(f"\n{'─'*60}")
            print(f"  Marker: {marker}  ({len(entries)} files)")
            print(f"{'─'*60}")
            print(f"  {'Filename':<60}  Source run dir")
            for e in entries:
                print(f"  {e.filename:<60}  {e.run_dir.name}")
        print()
        return 0

    # ── Create output dirs and link/copy ──────────────────────────────────
    action = "copy" if copy else "symlink"
    linked = 0
    skipped_existing = 0
    errors = 0

    for marker in markers:
        marker_dir = out_reads / marker
        marker_dir.mkdir(parents=True, exist_ok=True)

        for entry in sorted(by_marker[marker], key=lambda e: e.filename):
            dst = marker_dir / entry.filename

            if dst.exists() or dst.is_symlink():
                log.debug("Already exists, skipping: %s", dst.name)
                skipped_existing += 1
                continue

            try:
                if copy:
                    shutil.copy2(entry.src, dst)
                else:
                    dst.symlink_to(entry.src)
                linked += 1
            except Exception as exc:
                log.error("Failed to %s %s → %s: %s", action, entry.src, dst, exc)
                errors += 1

    # ── Done ──────────────────────────────────────────────────────────────
    log.info("Done.")
    log.info("  %s: %d files", action.capitalize() + "d", linked)
    if skipped_existing:
        log.info("  Already existed (skipped): %d", skipped_existing)
    if errors:
        log.error("  Errors: %d — check messages above", errors)
        return 1

    if had_warnings:
        log.warning("Completed with warnings — review messages above.")

    log.info("")
    log.info("Next step:")
    log.info("  python 01_make_manifests.py \\")
    log.info("      --reads-dir %s \\", out_reads)
    log.info("      --outdir qiime2/imported/ \\")
    log.info("      --markers %s", " ".join(markers))

    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="00_merge_run_dirs.py",
        description=(
            "Merge multiple flat sequencing run directories into a single "
            "reads/ tree ready for 01_make_manifests.py."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Dry run — preview only:
  python 00_merge_run_dirs.py --run-dirs run1/ run2/ --out-reads reads/ --dry-run

  # Real merge with symlinks (default):
  python 00_merge_run_dirs.py --run-dirs /data/run1 /data/run2 --out-reads reads/

  # Copy files (for archiving or portability):
  python 00_merge_run_dirs.py --run-dirs run1/ run2/ --out-reads reads/ --copy

  # Warn instead of aborting on duplicate filenames:
  python 00_merge_run_dirs.py --run-dirs run1/ run2/ --out-reads reads/ --on-collision warn

  # Then build manifests:
  python 01_make_manifests.py --reads-dir reads/ --outdir qiime2/imported/ --markers 16S cytb
""",
    )

    p.add_argument(
        "--run-dirs",
        nargs="+",
        required=True,
        type=Path,
        metavar="DIR",
        help="One or more run directories, each containing flat FASTQ files.",
    )
    p.add_argument(
        "--out-reads",
        type=Path,
        default=Path("reads"),
        metavar="DIR",
        help=(
            "Output reads directory. Per-marker subdirectories are created inside. "
            "Default: ./reads"
        ),
    )
    p.add_argument(
        "--markers",
        nargs="+",
        default=DEFAULT_MARKERS,
        metavar="MARKER",
        help=(
            f"Marker tokens to detect in filenames. "
            f"Default: {' '.join(DEFAULT_MARKERS)}"
        ),
    )
    p.add_argument(
        "--on-collision",
        choices=["abort", "warn", "skip"],
        default="abort",
        help=(
            "What to do when the same filename appears in >1 run dir. "
            "'abort' (default) stops immediately; 'warn' keeps the first copy "
            "and logs a warning; 'skip' keeps the first copy silently. "
            "A collision_report.tsv is always written."
        ),
    )
    p.add_argument(
        "--copy",
        action="store_true",
        default=False,
        help=(
            "Copy files instead of creating symlinks. "
            "Slower and uses more disk, but more portable. "
            "Default: symlink."
        ),
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Preview the merge plan without creating any files or directories.",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Enable debug logging (shows every file processed).",
    )

    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Resolve and validate run dirs
    run_dirs: List[Path] = []
    for rd in args.run_dirs:
        resolved = rd.resolve()
        if not resolved.is_dir():
            log.error("Run directory not found: %s", rd)
            return 1
        run_dirs.append(resolved)

    out_reads = args.out_reads.resolve()
    markers = args.markers

    log.info("Run directories (%d):", len(run_dirs))
    for rd in run_dirs:
        log.info("  %s", rd)
    log.info("Output reads dir : %s", out_reads)
    log.info("Markers          : %s", ", ".join(markers))
    log.info("On collision     : %s", args.on_collision)
    log.info("File action      : %s", "copy" if args.copy else "symlink")
    if args.dry_run:
        log.info("Mode             : DRY RUN")
    log.info("")

    return merge(
        run_dirs=run_dirs,
        out_reads=out_reads,
        markers=markers,
        on_collision=args.on_collision,
        copy=args.copy,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    sys.exit(main())
