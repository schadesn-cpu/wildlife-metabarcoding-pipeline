#!/usr/bin/env python3
"""
utils/validate.py
=================
Fail-loud input validation. Every script calls the relevant check on the way
in, so a malformed input stops the run *here* with a message that says exactly
what was expected and what was found — instead of surfacing as a cryptic
QIIME 2 traceback 200 lines later, or (worse) silently producing wrong output.

All checks raise ValidationError on failure. The message is written for the
person running the pipeline, not for a stack trace.

Public API:
    validate_reads_dir(reads_dir, markers)   -> dict[marker -> fastq count]
    validate_manifest(path)                  -> int (sample count)
    validate_metadata(path)                  -> list[str] (column names)
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Union

PathLike = Union[str, Path]


class ValidationError(Exception):
    """Raised when an input file or directory does not meet a step's contract."""


# ---------------------------------------------------------------------------
# Environment preflight
# ---------------------------------------------------------------------------

def require_qiime(min_version: Optional[str] = None) -> str:
    """
    Confirm the `qiime` CLI is available before any step that shells out to it.

    The whole pipeline drives QIIME 2 as a subprocess, so a missing or
    un-activated QIIME 2 environment means nothing past the pure-Python steps
    can run. This fails loud, up front, with a message that points at the fix
    rather than letting a 'command not found' surface mid-run.

    Returns the detected version string (or "unknown" if it could not be
    parsed — qiime being present is what matters; version is informational).
    Raises ValidationError if `qiime` is not on PATH.
    """
    exe = shutil.which("qiime")
    if exe is None:
        raise ValidationError(
            "QIIME 2 not found: the `qiime` command is not on your PATH.\n"
            "  This pipeline runs QIIME 2 as a subprocess, so it must be "
            "installed and activated first.\n"
            "  - conda:   conda activate qiime2-amplicon-<version>\n"
            "  - cluster: module load qiime2   (or your site's equivalent)\n"
            "  Then re-run."
        )

    try:
        out = subprocess.run(
            ["qiime", "--version"],
            capture_output=True, text=True, timeout=30,
        )
    except (OSError, subprocess.SubprocessError) as e:
        # qiime exists on PATH but would not report a version — surface it,
        # but don't block: presence is the hard requirement, version is a nicety.
        import logging
        logging.getLogger("validate").warning(
            "found qiime at %s but could not read its version: %s", exe, e
        )
        return "unknown"

    version = out.stdout.strip().splitlines()[0] if out.stdout.strip() else "unknown"

    if min_version and version != "unknown" and version < min_version:
        raise ValidationError(
            f"QIIME 2 version {version!r} is older than the required "
            f"{min_version!r}. Activate a newer environment and re-run."
        )
    return version


# ---------------------------------------------------------------------------
# Reads directory
# ---------------------------------------------------------------------------

_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _is_fastq(name: str) -> bool:
    n = name.lower()
    return any(n.endswith(s) for s in _FASTQ_SUFFIXES)


def validate_reads_dir(reads_dir: PathLike, markers: List[str]) -> Dict[str, int]:
    """
    Confirm the reads directory has a subfolder per marker, each containing
    fastq files. Returns {marker: n_fastq_files}.

    Expected layout:
        reads/
        ├── 16S/      *_R1_*.fastq.gz  *_R2_*.fastq.gz
        ├── MiFish/   ...
        └── ...
    """
    root = Path(reads_dir)
    if not root.is_dir():
        raise ValidationError(
            f"Reads directory not found: {root}\n"
            f"  Expected a directory with one subfolder per marker "
            f"({', '.join(markers)})."
        )

    counts: Dict[str, int] = {}
    missing: List[str] = []
    empty: List[str] = []
    for marker in markers:
        sub = root / marker
        if not sub.is_dir():
            missing.append(marker)
            continue
        n = sum(1 for p in sub.iterdir() if p.is_file() and _is_fastq(p.name))
        counts[marker] = n
        if n == 0:
            empty.append(marker)

    problems = []
    if missing:
        present = sorted(p.name for p in root.iterdir() if p.is_dir())
        problems.append(
            f"  missing marker subfolders: {', '.join(missing)}\n"
            f"  subfolders present in {root}: {', '.join(present) or '(none)'}"
        )
    if empty:
        problems.append(f"  marker subfolders contain no fastq files: {', '.join(empty)}")
    if problems:
        raise ValidationError("Reads directory layout problem:\n" + "\n".join(problems))

    return counts


# ---------------------------------------------------------------------------
# QIIME 2 paired-end manifest (PairedEndFastqManifestPhred33V2)
# ---------------------------------------------------------------------------

_MANIFEST_COLUMNS = [
    "sample-id",
    "forward-absolute-filepath",
    "reverse-absolute-filepath",
]


def validate_manifest(path: PathLike) -> int:
    """
    Confirm a manifest TSV has exactly the required header columns and that
    every read path is absolute and exists. Returns the sample count.

    This catches the two failures that account for most manifest import errors:
    a relative path, and a path that points at a symlink whose target is gone.
    """
    p = Path(path)
    if not p.exists():
        raise ValidationError(f"Manifest not found: {p}")

    with p.open("r", encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        raise ValidationError(f"Manifest is empty: {p}")

    header = lines[0].split("\t")
    if header != _MANIFEST_COLUMNS:
        raise ValidationError(
            f"Manifest header is wrong: {p}\n"
            f"  expected (tab-separated): {_MANIFEST_COLUMNS}\n"
            f"  found:                    {header}"
        )

    bad_paths: List[str] = []
    for ln in lines[1:]:
        cells = ln.split("\t")
        if len(cells) != 3:
            raise ValidationError(
                f"Manifest row does not have 3 columns: {p}\n  offending row: {ln!r}"
            )
        sid, fwd, rev = cells
        for direction, fp in (("forward", fwd), ("reverse", rev)):
            fpp = Path(fp)
            if not fpp.is_absolute():
                bad_paths.append(f"  {sid} ({direction}): not an absolute path -> {fp}")
            elif not fpp.exists():
                bad_paths.append(f"  {sid} ({direction}): file does not exist -> {fp}")

    if bad_paths:
        raise ValidationError(
            f"Manifest references {len(bad_paths)} unusable read path(s): {p}\n"
            + "\n".join(bad_paths[:20])
            + ("\n  ... (more)" if len(bad_paths) > 20 else "")
        )

    return len(lines) - 1


# ---------------------------------------------------------------------------
# QIIME 2 sample metadata
# ---------------------------------------------------------------------------

# QIIME 2 accepts any one of these (case-insensitive) as the first-column id header.
_VALID_ID_HEADERS = {
    "id", "sampleid", "sample id", "sample-id", "sample_name",
    "featureid", "feature id", "feature-id", "#sampleid", "#otuid", "#otu id",
}


def validate_metadata(path: PathLike) -> List[str]:
    """
    Confirm a QIIME 2 metadata TSV is usable: the FIRST column header must be a
    recognized id column, and if a '#q2:types' directive row is present it must
    be the SECOND line. Returns the list of column names.
    """
    p = Path(path)
    if not p.exists():
        raise ValidationError(f"Metadata not found: {p}")

    with p.open("r", encoding="utf-8") as fh:
        rows = [ln.rstrip("\n") for ln in fh]
    rows = [r for r in rows if r.strip() != ""]
    if not rows:
        raise ValidationError(f"Metadata is empty: {p}")

    columns = rows[0].split("\t")
    first = columns[0].strip().lower()
    if first not in _VALID_ID_HEADERS:
        raise ValidationError(
            f"Metadata first column is not a recognized id header: {p}\n"
            f"  found first column:  {columns[0]!r}\n"
            f"  must be one of:      sample-id, id, #SampleID, sampleid, "
            f"featureid (case-insensitive)\n"
            f"  (QIIME 2 keys every sample on the first column — it cannot be a "
            f"data column.)"
        )

    # If a types directive exists, it must be row 2.
    for i, row in enumerate(rows[1:], start=2):
        if row.split("\t")[0].strip().lower() == "#q2:types":
            if i != 2:
                raise ValidationError(
                    f"Metadata '#q2:types' directive must be the second line "
                    f"(found on line {i}): {p}"
                )
            break

    return columns
