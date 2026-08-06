#!/usr/bin/env python3
"""
make_manifests.py
=================
Step: manifest

Purpose:
    Build one QIIME 2 PairedEndFastqManifestPhred33V2 manifest TSV per marker
    from a reads directory organized into per-marker subfolders. First step of
    an end-to-end run: point it at the raw reads and it produces the manifests
    the import step consumes.

Inputs:
    reads/<marker>/        paired *_R1_*/*_R2_* fastq.gz files (one folder per marker)
    pipeline_config.yml    active_markers, samples.control_prefixes, project.root

Outputs:
    qiime2/<marker>/imported/manifest_<marker>.tsv    one per marker
    logs/run_manifest.jsonl                           run appended on success

Read pairing handles R1/R2 detection across common Illumina naming schemes,
de-duplicates symlink/original collisions, and derives clean sample ids by
stripping lane/index/marker tokens (no hardcoded sample prefix).

Usage:
    python pipeline.py run --step manifest
    python make_manifests.py [--reads-dir reads/] [--markers 16S MiFish] [--dry-run]

Requirements:
    Python >= 3.8; project config_loader and the utils package on the path.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

# --- make config_loader and the utils package importable regardless of cwd ---
# This script lives in scripts/, alongside config_loader.py and the utils/ package.
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from config_loader import load_config, get_paths              # noqa: E402
from utils import validate, checkpoint, provenance            # noqa: E402

log = logging.getLogger("make_manifests")


# ===========================================================================
# Read pairing
# ===========================================================================

class ReadPair(NamedTuple):
    sample_id: str
    forward: Path
    reverse: Path


_DIRECTION_PATTERNS: List[Tuple[re.Pattern, str, str]] = [
    (re.compile(r"_(R[12])_\d+"), "R1", "R2"),       # Illumina:  _R1_001 / _R2_001
    (re.compile(r"_(R[12])(?=[._]|$)"), "R1", "R2"),  # short:     _R1 / _R2
    (re.compile(r"_([12])(?=[._]|$)"), "1", "2"),     # numeric:   _1 / _2
]

_FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _is_fastq(path: Path) -> bool:
    return any(path.name.lower().endswith(ext) for ext in _FASTQ_EXTS)


def _strip_fastq_ext(name: str) -> str:
    for ext in _FASTQ_EXTS:
        if name.lower().endswith(ext):
            return name[: -len(ext)]
    return name


def _get_direction(stem: str) -> Optional[str]:
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
    for pattern, _, _ in _DIRECTION_PATTERNS:
        m = pattern.search(stem)
        if m:
            return stem[: m.start()]
    return stem


def _extract_sample_id(base_key: str, marker: str) -> str:
    """
    Derive a clean sample id from a pairing key by stripping common
    Illumina/sequencer suffixes and the marker token. General by design — it
    keys off the marker name and standard _L###/_S### suffixes, not a fixed
    sample prefix, so it works for any naming scheme.

        TV230084-GI-16S_S1492_L002  ->  TV230084-GI
        NTC-16S_S99_L002            ->  NTC
    """
    sid = base_key
    sid = re.sub(r"_L\d{3}$", "", sid)                 # lane
    sid = re.sub(r"_S\d+$", "", sid)                   # Illumina sample index
    sid = re.sub(rf"[-_]{re.escape(marker)}$", "", sid, flags=re.IGNORECASE)  # marker
    sid = sid.rstrip("-_")
    return sid if sid else base_key


def _resolve(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"Symlink target does not exist: {path} -> {resolved}")
    return resolved


def find_pairs(marker_dir: Path, marker: str) -> Tuple[List[ReadPair], List[str]]:
    """Scan marker_dir for R1/R2 fastq pairs. Returns (pairs, warnings)."""
    fastq_files = [f for f in marker_dir.iterdir()
                   if (f.is_file() or f.is_symlink()) and _is_fastq(f)]
    if not fastq_files:
        return [], [f"No FASTQ files found in {marker_dir}"]

    # De-duplicate symlink/original pairs that resolve to the same real file.
    seen_real: Dict[Path, Path] = {}
    deduped: List[Path] = []
    for f in fastq_files:
        try:
            real = f.resolve()
        except OSError as e:
            log.debug("could not resolve %s (%s); kept for pairing, validated later", f.name, e)
            real = f
        if real in seen_real:
            existing = seen_real[real]
            if f.is_symlink() and not existing.is_symlink():
                continue
            if not f.is_symlink() and existing.is_symlink():
                seen_real[real] = f
                deduped = [f if x == existing else x for x in deduped]
            continue
        seen_real[real] = f
        deduped.append(f)

    buckets: Dict[str, Dict[str, Path]] = {}
    unrecognized: List[str] = []
    for f in sorted(deduped):
        stem = _strip_fastq_ext(f.name)
        direction = _get_direction(stem)
        if direction is None:
            unrecognized.append(f.name)
            continue
        base_key = _strip_direction(stem)
        slot = buckets.setdefault(base_key, {})
        if direction in slot:
            existing = slot[direction]
            if not f.is_symlink() and existing.is_symlink():
                slot[direction] = f
        else:
            slot[direction] = f

    pairs: List[ReadPair] = []
    warnings: List[str] = []
    if unrecognized:
        warnings.append(
            f"Could not determine read direction for {len(unrecognized)} file(s): "
            + ", ".join(unrecognized)
        )
    for base_key in sorted(buckets):
        r1, r2 = buckets[base_key].get("R1"), buckets[base_key].get("R2")
        if r1 is None:
            warnings.append(f"No R1 found for key '{base_key}' (R2: {r2.name})")
            continue
        if r2 is None:
            warnings.append(f"No R2 found for key '{base_key}' (R1: {r1.name})")
            continue
        try:
            r1, r2 = _resolve(r1), _resolve(r2)
        except FileNotFoundError as e:
            warnings.append(str(e))
            continue
        pairs.append(ReadPair(_extract_sample_id(base_key, marker), r1, r2))

    return pairs, warnings


_MANIFEST_HEADER = "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"


def write_manifest(pairs: List[ReadPair], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(_MANIFEST_HEADER + "\n")
        for p in pairs:
            fh.write(f"{p.sample_id}\t{p.forward}\t{p.reverse}\n")


# ===========================================================================
# Orchestration
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="make_manifests.py",
        description="Build QIIME 2 paired-end manifests, one per marker.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--reads-dir", default=None,
                   help="Reads directory (default: 'reads/' under project root).")
    p.add_argument("--markers", nargs="+", default=None,
                   help="Markers to build (default: active_markers from config).")
    p.add_argument("--dry-run", action="store_true",
                   help="Report pairs without writing manifests.")
    p.add_argument("--config", default=None, help="Path to pipeline_config.yml.")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s  %(message)s",
    )

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    markers = args.markers or cfg.active_markers
    reads_dir = Path(args.reads_dir) if args.reads_dir else cfg.resolve("reads")

    # --- validate on entry: fail loud with a clear message ------------------
    try:
        counts = validate.validate_reads_dir(reads_dir, markers)
    except validate.ValidationError as e:
        log.error("Input validation failed:\n%s", e)
        return 2
    log.info("Reads directory OK: %s",
             ", ".join(f"{m}={counts[m]} files" for m in markers))

    control_prefixes = tuple(cfg.samples.get("control_prefixes", []))
    written: List[Path] = []

    for marker in markers:
        marker_dir = reads_dir / marker
        pairs, warnings = find_pairs(marker_dir, marker)
        for w in warnings:
            log.warning("[%s] %s", marker, w)
        if not pairs:
            log.error("[%s] no read pairs found — skipping.", marker)
            continue

        n_control = sum(1 for p in pairs
                        if any(p.sample_id.startswith(c) for c in control_prefixes))
        log.info("[%s] %d sample pairs (%d controls)", marker, len(pairs), n_control)

        out_path = paths.engine_manifest_tsv(marker)
        if args.dry_run:
            log.info("[%s] DRY RUN — would write %s", marker, out_path)
            continue

        write_manifest(pairs, out_path)
        # validate what we just wrote — catches bad/relative/missing paths now,
        # not at qiime import time.
        try:
            n = validate.validate_manifest(out_path)
        except validate.ValidationError as e:
            log.error("[%s] manifest failed validation after writing:\n%s", marker, e)
            return 3
        log.info("[%s] wrote and validated manifest: %s (%d samples)", marker, out_path, n)
        written.append(out_path)

    if args.dry_run:
        log.info("Dry run complete — no files written.")
        return 0
    if not written:
        log.error("No manifests were written.")
        return 1

    # --- checkpoint + provenance on exit ------------------------------------
    checkpoint.print_checkpoint(
        cfg,
        "manifest",
        produced=written,
        provenance={
            "outputs": written,
            "command": "python " + " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
            "extra": {"markers": markers},
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
