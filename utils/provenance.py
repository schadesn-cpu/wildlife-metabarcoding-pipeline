#!/usr/bin/env python3
"""
utils/provenance.py
===================
Append-only run ledger for the pipeline.

Every time a step finishes it appends ONE JSON line to logs/run_manifest.jsonl
recording exactly what it ran with: the marker, the analysis, the read depth,
the *identity* of every input artifact (QIIME 2 UUID, not just a path that
could later be overwritten), the outputs produced, the pipeline git commit,
and the command line. This is what turns "which rep-seqs produced this PERMANOVA
at depth 17000?" into a one-line grep.

This complements — does not replace — QIIME 2's own provenance. Every .qza/.qzv
already embeds its full internal command graph (open one in QIIME 2 View →
Provenance tab). What QIIME 2 does NOT track is the figure/stats half of the
pipeline (matplotlib, pandas, R). The ledger spans the whole thing.

Public API:
    record_run(cfg, step, ...)   append one run record; returns the dict written
    ledger_path(cfg)             Path to logs/run_manifest.jsonl
    qza_uuid(path)               read a QIIME 2 artifact UUID (no QIIME install needed)
    read_runs(cfg, **filters)    load records back, optionally filtered
"""

from __future__ import annotations

import json
import subprocess
import uuid as _uuid
import zipfile
from datetime import datetime, timezone
from pathlib import Path
import logging
from typing import Any, Dict, List, Optional, Sequence, Union

log = logging.getLogger("provenance")

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Where the ledger lives
# ---------------------------------------------------------------------------

def ledger_path(cfg) -> Path:
    """logs/run_manifest.jsonl under the project root."""
    return cfg.resolve("logs/run_manifest.jsonl")


# ---------------------------------------------------------------------------
# Reading a QIIME 2 artifact's identity without importing QIIME 2
#
# A .qza/.qzv is a zip archive whose every entry lives under a top-level
# directory named with the artifact's UUID. So the UUID is simply the first
# path component of the first archive member — readable with the stdlib alone.
# ---------------------------------------------------------------------------

def qza_uuid(path: PathLike) -> Optional[str]:
    """
    Return the QIIME 2 artifact UUID for a .qza/.qzv file, or None if the file
    is missing or not a readable QIIME 2 archive. Never raises — provenance
    capture must not be able to crash a successful run.
    """
    p = Path(path)
    if not p.exists():
        log.warning("provenance: input artifact not found, UUID not captured: %s", p)
        return None
    try:
        with zipfile.ZipFile(p) as zf:
            names = zf.namelist()
    except (zipfile.BadZipFile, OSError) as e:
        log.warning("provenance: could not read artifact UUID from %s: %s", p, e)
        return None
    if not names:
        log.warning("provenance: artifact archive is empty: %s", p)
        return None
    top = names[0].split("/", 1)[0]
    if len(top) == 36 and top.count("-") == 4:
        return top
    log.warning("provenance: first archive entry is not a UUID directory in %s: %r", p, top)
    return top or None


# ---------------------------------------------------------------------------
# Pipeline version stamp
# ---------------------------------------------------------------------------

def git_commit(cfg) -> Optional[str]:
    """Short git commit of the pipeline at project root, or None if not a repo."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(cfg.root),
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError) as e:
        log.warning("provenance: git unavailable, commit not captured: %s", e)
        return None
    if out.returncode != 0:
        log.warning("provenance: not a git repository at %s, commit not captured", cfg.root)
        return None
    return out.stdout.strip() or None


# ---------------------------------------------------------------------------
# Describing an input artifact (path + identity)
# ---------------------------------------------------------------------------

def _describe_input(path: PathLike, cfg) -> Dict[str, Any]:
    """
    Build an {path, uuid?} record for one input. For .qza/.qzv inputs the UUID
    is read so the run is tied to the artifact's identity, not just its location.
    """
    p = Path(path)
    rec: Dict[str, Any] = {"path": _relativize(p, cfg)}
    if p.suffix in (".qza", ".qzv"):
        u = qza_uuid(p)
        if u:
            rec["uuid"] = u
    return rec


def _relativize(p: PathLike, cfg) -> str:
    """Store paths relative to project root when possible — keeps the ledger portable."""
    p = Path(p)
    try:
        return str(p.resolve().relative_to(cfg.root.resolve()))
    except (ValueError, OSError):
        return str(p)


# ---------------------------------------------------------------------------
# The main entry point
# ---------------------------------------------------------------------------

def record_run(
    cfg,
    step: str,
    *,
    marker: Optional[str] = None,
    analysis: Optional[str] = None,
    read_depth: Optional[int] = None,
    inputs: Optional[Dict[str, PathLike]] = None,
    outputs: Optional[Sequence[PathLike]] = None,
    command: Optional[str] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Append one run record to the ledger and return the dict written.

    Parameters
    ----------
    step        step key, e.g. "denoise"
    marker      e.g. "MiFish"
    analysis    e.g. "DvT", "season", "all"
    read_depth  rarefaction depth used, if any
    inputs      mapping of role -> path, e.g.
                {"table": paths.dada2_table("MiFish"),
                 "rep_seqs": paths.rep_seqs_qza("MiFish")}
                .qza/.qzv inputs get their UUID captured automatically.
    outputs     list of output paths produced
    command     the exact command line that produced this run
    extra       any extra key/values to fold into the record

    Failure to write the ledger is logged to stderr but never raised — a
    bookkeeping failure must not sink a successful analysis step.
    """
    now = datetime.now(timezone.utc)
    record: Dict[str, Any] = {
        "run_id": now.strftime("%Y%m%d-%H%M%S") + "-" + _uuid.uuid4().hex[:4],
        "timestamp": now.isoformat(timespec="seconds"),
        "step": step,
    }
    if marker is not None:
        record["marker"] = marker
    if analysis is not None:
        record["analysis"] = analysis
    if read_depth is not None:
        record["read_depth"] = read_depth
    if inputs:
        record["inputs"] = {role: _describe_input(p, cfg)
                            for role, p in inputs.items() if p is not None}
    if outputs:
        record["outputs"] = [_relativize(p, cfg) for p in outputs]
    if command is not None:
        record["command"] = command

    commit = git_commit(cfg)
    if commit:
        record["pipeline_commit"] = commit
    if extra:
        record.update(extra)

    try:
        path = ledger_path(cfg)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(record) + "\n")
    except OSError as e:
        log.error("provenance: could not write run ledger %s: %s", ledger_path(cfg), e)

    return record


def read_runs(cfg, **filters: Any) -> List[Dict[str, Any]]:
    """
    Load all run records, optionally filtered by exact-match top-level fields.

    Example:
        read_runs(cfg, step="denoise", marker="MiFish")
    """
    path = ledger_path(cfg)
    if not path.exists():
        return []
    records: List[Dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError as e:
                log.warning("provenance: skipping malformed ledger line %d: %s", lineno, e)
                continue
            if all(rec.get(k) == v for k, v in filters.items()):
                records.append(rec)
    return records
