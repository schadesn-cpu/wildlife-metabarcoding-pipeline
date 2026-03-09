#!/usr/bin/env python3
"""
01_make_manifests.py
====================
Auto-build QIIME 2 PairedEndFastqManifestPhred33V2 manifest TSVs from a
reads directory organized into per-marker subdirectories.

Expected reads directory layout:
    reads/
    ├── 16S/
    │   ├── TV230084-GI-16S_S1492_L002_R1_001.fastq.gz
    │   ├── TV230084-GI-16S_S1492_L002_R2_001.fastq.gz
    │   └── ...
    ├── 18S/
    ├── ITS1-2/
    ├── MiFish/
    └── cytb/

Output manifest format (PairedEndFastqManifestPhred33V2):
    sample-id    forward-absolute-filepath    reverse-absolute-filepath

One manifest TSV is written per marker to the output directory:
    qiime2/imported/manifest_16S.tsv
    qiime2/imported/manifest_18S.tsv
    ...

Features:
  - Auto-detects R1/R2 pairs using common Illumina naming patterns
  - Resolves symlinks to real absolute paths (prevents QIIME import errors)
  - Extracts clean sample IDs by stripping marker/lane/index suffixes
  - Validates every sample has both R1 and R2 before writing
  - Reports unpaired files and skipped non-FASTQ files
  - Dry-run mode to preview without writing

Usage:
    # Auto-detect all markers
    python scripts/01_make_manifests.py \\
        --reads-dir reads/ \\
        --outdir qiime2/imported/

    # Specific markers only
    python scripts/01_make_manifests.py \\
        --reads-dir reads/ \\
        --outdir qiime2/imported/ \\
        --markers 16S 18S MiFish

    # Preview without writing
    python scripts/01_make_manifests.py \\
        --reads-dir reads/ \\
        --outdir qiime2/imported/ \\
        --dry-run

Requirements:
    Python >= 3.8 (standard library only)
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
class ReadPair(NamedTuple):
    sample_id: str
    forward: Path   # resolved absolute path
    reverse: Path   # resolved absolute path


# ---------------------------------------------------------------------------
# R1/R2 detection
# ---------------------------------------------------------------------------

# Patterns that identify the read direction in a filename.
# Each tuple is (compiled_regex, r1_value, r2_value).
# The regex must capture a group that matches one of r1_value or r2_value.
_DIRECTION_PATTERNS: List[Tuple[re.Pattern, str, str]] = [
    # Illumina standard:  _R1_001 / _R2_001
    (re.compile(r"_(R[12])_\d+"), "R1", "R2"),
    # Short form:         _R1 / _R2  (before extension or end of stem)
    (re.compile(r"_(R[12])(?=[._]|$)"), "R1", "R2"),
    # Numeric only:       _1 / _2   (less common but seen on some platforms)
    (re.compile(r"_([12])(?=[._]|$)"), "1", "2"),
]


def _get_direction(stem: str) -> Optional[str]:
    """
    Return 'R1', 'R2', or None if direction cannot be determined from the
    file stem (filename without extensions).
    """
    for pattern, r1_val, r2_val in _DIRECTION_PATTERNS:
        m = pattern.search(stem)
        if m:
            val = m.group(1)
            if val == r1_val:
                return "R1"
            if val == r2_val:
                return "R2"
    return None


def _strip_direction(stem: str) -> str:
    """
    Remove the read direction token and everything after it
    (lane, index, suffix tokens) to produce a base sample key used for
    pairing R1 and R2 files that belong to the same sample.

    Example:
        TV230084-GI-16S_S1492_L002_R1_001  →  TV230084-GI-16S_S1492
    """
    # Try each direction pattern; strip from the first match onward
    for pattern, _, _ in _DIRECTION_PATTERNS:
        m = pattern.search(stem)
        if m:
            return stem[: m.start()]
    return stem


# ---------------------------------------------------------------------------
# Sample ID extraction
# ---------------------------------------------------------------------------

def _extract_sample_id(base_key: str, marker: str) -> str:
    """
    Derive a clean sample ID from a base pairing key by removing common
    Illumina/sequencer suffixes and the marker token.

    Example inputs → outputs:
        TV230084-GI-16S_S1492_L002   →  TV230084-GI
        TV230084-GI-16S_S1492        →  TV230084-GI
        TV230084_S1492               →  TV230084
        NTC-16S_S99_L002             →  NTC
        TV230084-GI-16S              →  TV230084-GI

    Strategy (applied in order):
        1. Remove _L00N lane suffix
        2. Remove _S<number> Illumina sample index suffix
        3. Remove -<MARKER> or _<MARKER> token (case-insensitive)
        4. Strip trailing separators
    """
    sid = base_key

    # Remove _L00N lane identifier (e.g., _L002)
    sid = re.sub(r"_L\d{3}$", "", sid)

    # Remove _S<digits> Illumina sample index
    sid = re.sub(r"_S\d+$", "", sid)

    # Remove marker token: -16S, _16S, -ITS1-2, _MiFish, etc.
    marker_escaped = re.escape(marker)
    sid = re.sub(rf"[-_]{marker_escaped}$", "", sid, flags=re.IGNORECASE)

    # Strip trailing separators
    sid = sid.rstrip("-_")

    return sid if sid else base_key  # fallback: keep original if stripping left nothing


# ---------------------------------------------------------------------------
# Core pairing logic
# ---------------------------------------------------------------------------

FASTQ_EXTENSIONS = {".fastq", ".fastq.gz", ".fq", ".fq.gz"}


def _is_fastq(path: Path) -> bool:
    """Return True if path looks like a FASTQ file."""
    name = path.name.lower()
    return any(name.endswith(ext) for ext in FASTQ_EXTENSIONS)


def _strip_fastq_ext(name: str) -> str:
    """Remove FASTQ extension(s) to get the file stem."""
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.lower().endswith(ext):
            return name[: -len(ext)]
    return name


def _resolve(path: Path) -> Path:
    """
    Resolve a path to its real absolute path, following symlinks.
    Raises FileNotFoundError if the resolved target does not exist.
    """
    resolved = path.resolve()
    if not resolved.exists():
        raise FileNotFoundError(
            f"Symlink target does not exist: {path} -> {resolved}"
        )
    return resolved


def find_pairs(
    marker_dir: Path,
    marker: str,
    resolve_symlinks: bool = True,
) -> Tuple[List[ReadPair], List[str]]:
    """
    Scan marker_dir for R1/R2 FASTQ pairs.

    Returns:
        pairs   — list of ReadPair (sample_id, forward, reverse)
        warnings — list of warning strings for unpaired / unrecognized files
    """
    if not marker_dir.is_dir():
        raise NotADirectoryError(f"Reads directory does not exist: {marker_dir}")

    # Collect all FASTQ files (non-recursive — reads should be flat per marker)
    fastq_files = [f for f in marker_dir.iterdir() if f.is_file() or f.is_symlink()]
    fastq_files = [f for f in fastq_files if _is_fastq(f)]

    if not fastq_files:
        return [], [f"No FASTQ files found in {marker_dir}"]

    # Deduplicate: if a symlink and an original resolve to the same real path,
    # keep only the original (non-symlink). This prevents double-counting when
    # both TV230084_R1.fastq.gz (symlink) and TV230084_L002_R1_001.fastq.gz
    # (original) exist in the same directory.
    seen_real: Dict[Path, Path] = {}  # real_path -> chosen path
    deduped: List[Path] = []
    for f in fastq_files:
        try:
            real = f.resolve()
        except OSError:
            real = f  # broken symlink — keep it so we can warn later
        if real in seen_real:
            existing = seen_real[real]
            # Prefer the non-symlink (original)
            if f.is_symlink() and not existing.is_symlink():
                continue  # skip this symlink
            elif not f.is_symlink() and existing.is_symlink():
                seen_real[real] = f  # replace symlink with real file
                deduped = [f if x == existing else x for x in deduped]
                continue
            else:
                log.debug("Duplicate real path %s — keeping %s, skipping %s", real, existing.name, f.name)
                continue
        seen_real[real] = f
        deduped.append(f)
    fastq_files = deduped

    # Bucket by base_key → {R1: path, R2: path}
    buckets: Dict[str, Dict[str, Path]] = {}
    unrecognized: List[str] = []

    for f in sorted(fastq_files):
        stem = _strip_fastq_ext(f.name)
        direction = _get_direction(stem)

        if direction is None:
            unrecognized.append(f.name)
            continue

        base_key = _strip_direction(stem)

        if base_key not in buckets:
            buckets[base_key] = {}

        if direction in buckets[base_key]:
            # Duplicate — prefer non-symlink (original file)
            existing = buckets[base_key][direction]
            if f.is_symlink() and not existing.is_symlink():
                log.debug("Skipping symlink duplicate: %s (keeping %s)", f.name, existing.name)
                continue
            elif not f.is_symlink() and existing.is_symlink():
                log.debug("Replacing symlink with real file: %s -> %s", existing.name, f.name)
                buckets[base_key][direction] = f
            else:
                log.warning(
                    "Duplicate %s for key '%s': %s and %s — keeping first",
                    direction, base_key, existing.name, f.name,
                )
        else:
            buckets[base_key][direction] = f

    # Build pairs and collect warnings for unpaired files
    pairs: List[ReadPair] = []
    warnings: List[str] = []

    if unrecognized:
        warnings.append(
            f"Could not determine read direction for {len(unrecognized)} file(s): "
            + ", ".join(unrecognized)
        )

    for base_key in sorted(buckets.keys()):
        bucket = buckets[base_key]
        r1 = bucket.get("R1")
        r2 = bucket.get("R2")

        if r1 is None and r2 is not None:
            warnings.append(f"No R1 found for key '{base_key}' (R2: {r2.name})")
            continue
        if r2 is None and r1 is not None:
            warnings.append(f"No R2 found for key '{base_key}' (R1: {r1.name})")
            continue

        # Resolve symlinks if requested
        if resolve_symlinks:
            try:
                r1 = _resolve(r1)
                r2 = _resolve(r2)
            except FileNotFoundError as e:
                warnings.append(str(e))
                continue
        else:
            r1 = r1.resolve()
            r2 = r2.resolve()

        sample_id = _extract_sample_id(base_key, marker)

        pairs.append(ReadPair(
            sample_id=sample_id,
            forward=r1,
            reverse=r2,
        ))

    return pairs, warnings


# ---------------------------------------------------------------------------
# Manifest writing
# ---------------------------------------------------------------------------

MANIFEST_HEADER = "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"


def write_manifest(pairs: List[ReadPair], out_path: Path) -> None:
    """Write a QIIME 2 PairedEndFastqManifestPhred33V2 TSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(MANIFEST_HEADER + "\n")
        for pair in pairs:
            fh.write(f"{pair.sample_id}\t{pair.forward}\t{pair.reverse}\n")
    log.info("Wrote manifest: %s (%d samples)", out_path, len(pairs))


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_manifest(pairs: List[ReadPair]) -> List[str]:
    """
    Run basic sanity checks on a list of pairs before writing.
    Returns a list of error strings (empty = all good).
    """
    errors: List[str] = []
    seen_ids: Dict[str, str] = {}  # sample_id -> base_key (for duplicate reporting)

    for pair in pairs:
        # Check for duplicate sample IDs (can happen if ID extraction is too aggressive)
        if pair.sample_id in seen_ids:
            errors.append(
                f"Duplicate sample-id '{pair.sample_id}' — "
                f"check that --markers matches the marker token in your filenames"
            )
        else:
            seen_ids[pair.sample_id] = str(pair.forward)

        # Check files are readable
        for filepath in (pair.forward, pair.reverse):
            if not filepath.exists():
                errors.append(f"File does not exist: {filepath}")
            elif not os.access(filepath, os.R_OK):
                errors.append(f"File is not readable: {filepath}")

    return errors


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="01_make_manifests.py",
        description=(
            "Auto-build QIIME 2 PairedEndFastqManifestPhred33V2 manifest TSVs "
            "from a reads directory organized into per-marker subdirectories."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect all markers
  python scripts/01_make_manifests.py --reads-dir reads/ --outdir qiime2/imported/

  # Specific markers only
  python scripts/01_make_manifests.py --reads-dir reads/ --outdir qiime2/imported/ --markers 16S 18S

  # Preview without writing
  python scripts/01_make_manifests.py --reads-dir reads/ --outdir qiime2/imported/ --dry-run
""",
    )

    p.add_argument(
        "--reads-dir",
        type=Path,
        default=Path("reads"),
        metavar="DIR",
        help="Root reads directory containing per-marker subdirectories. Default: ./reads",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("qiime2/imported"),
        metavar="DIR",
        help="Directory to write manifest TSVs. Default: ./qiime2/imported",
    )
    p.add_argument(
        "--markers",
        nargs="+",
        default=None,
        metavar="MARKER",
        help=(
            "Markers to process. Must match subdirectory names under --reads-dir "
            "(e.g., 16S 18S ITS1-2 MiFish cytb). "
            "If omitted, all subdirectories under --reads-dir are processed."
        ),
    )
    p.add_argument(
        "--no-resolve-symlinks",
        action="store_true",
        default=False,
        help=(
            "Do NOT resolve symlinks to their real paths. "
            "By default symlinks are resolved so QIIME import does not fail on HPC systems."
        ),
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Print detected pairs without writing any files.",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Print each detected pair during processing.",
    )

    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    reads_dir = args.reads_dir.resolve()
    outdir = args.outdir.resolve()
    resolve_symlinks = not args.no_resolve_symlinks

    if not reads_dir.is_dir():
        log.error("Reads directory does not exist: %s", reads_dir)
        return 1

    # Determine which markers to process
    if args.markers:
        markers = args.markers
        # Validate that the subdirectories exist
        missing = [m for m in markers if not (reads_dir / m).is_dir()]
        if missing:
            log.error(
                "The following marker directories were not found under %s: %s",
                reads_dir,
                ", ".join(missing),
            )
            log.error(
                "Available subdirectories: %s",
                ", ".join(d.name for d in sorted(reads_dir.iterdir()) if d.is_dir()),
            )
            return 1
    else:
        markers = sorted(d.name for d in reads_dir.iterdir() if d.is_dir())
        if not markers:
            log.error("No subdirectories found under %s", reads_dir)
            return 1
        log.info("Auto-detected markers: %s", ", ".join(markers))

    # -----------------------------------------------------------------------
    # Process each marker
    # -----------------------------------------------------------------------
    overall_ok = True

    for marker in markers:
        marker_dir = reads_dir / marker
        log.info("--- Processing marker: %s (%s) ---", marker, marker_dir)

        try:
            pairs, warnings = find_pairs(
                marker_dir,
                marker=marker,
                resolve_symlinks=resolve_symlinks,
            )
        except NotADirectoryError as e:
            log.error("%s", e)
            overall_ok = False
            continue

        # Report warnings (unpaired files, unrecognized names, broken symlinks)
        for w in warnings:
            log.warning("[%s] %s", marker, w)

        if not pairs:
            log.error("[%s] No valid pairs found — skipping manifest", marker)
            overall_ok = False
            continue

        # Print pairs if verbose or dry-run
        if args.verbose or args.dry_run:
            print(f"\n{'sample-id':<30} {'R1 filename':<50} {'R2 filename'}")
            print("-" * 120)
            for pair in pairs:
                print(
                    f"{pair.sample_id:<30} "
                    f"{pair.forward.name:<50} "
                    f"{pair.reverse.name}"
                )

        # Validate
        errors = validate_manifest(pairs)
        if errors:
            log.error("[%s] Validation failed:", marker)
            for e in errors:
                log.error("  %s", e)
            overall_ok = False
            continue

        log.info("[%s] Found %d valid sample pairs", marker, len(pairs))

        # Write (or skip in dry-run)
        out_path = outdir / f"manifest_{marker}.tsv"

        if args.dry_run:
            log.info("[%s] DRY RUN — would write: %s", marker, out_path)
        else:
            write_manifest(pairs, out_path)

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print()
    if args.dry_run:
        log.info("Dry run complete — no files written.")
    elif overall_ok:
        log.info("All manifests written successfully to: %s", outdir)
    else:
        log.warning(
            "Completed with errors — check warnings above. "
            "Manifests for failed markers were not written."
        )

    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
