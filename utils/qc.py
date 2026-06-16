#!/usr/bin/env python3
"""
utils/qc.py
===========
Output quality gates. Where validate.py asks "is this input well-formed?",
qc.py asks "is this result scientifically acceptable?" — the checks that catch
a run which completed without erroring but produced bad data.

Each function returns a QCResult describing what it found. It does NOT decide
whether to stop the pipeline: that policy belongs to the calling step, which
can warn-and-continue or escalate to fatal based on config. Keeping assessment
(here) separate from policy (the step) means thresholds and severity stay
configurable without touching the check logic.

Quality gates, by the step that owns them:
    denoise   -> check_retention          (this module)
    trim      -> check_primer_dimers       (lands with the primer-advisor port)
    taxonomy  -> check_classifier_coverage (lands with the taxonomy port; trips
                                            the BLAST verification fallback)

Public API:
    check_retention(stats, ...) -> QCResult
"""

from __future__ import annotations

import csv
import io
import logging
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Sequence, Tuple, Union

log = logging.getLogger("qc")

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

@dataclass
class QCResult:
    """
    The outcome of one quality gate.

    status   "pass"  — nothing flagged
             "warn"  — something worth a human look, but not necessarily fatal
             "fail"  — a hard problem (a check sets this only when the result is
                       unusable regardless of policy)
    flagged  list of (sample_id, value) the check singled out
    controls list of (sample_id, value) for negative controls, reported
             separately so expected low-retention NTCs are not counted as
             problems
    """
    check: str
    status: str
    summary: str
    flagged: List[Tuple[str, float]] = field(default_factory=list)
    controls: List[Tuple[str, float]] = field(default_factory=list)

    def ok(self) -> bool:
        return self.status == "pass"

    def report(self) -> str:
        lines = [f"[QC: {self.check}] {self.status.upper()} — {self.summary}"]
        for sid, val in self.flagged:
            lines.append(f"    flagged: {sid}  ({val:.1f}%)")
        for sid, val in self.controls:
            lines.append(f"    control: {sid}  ({val:.1f}%)  [expected low]")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Loading a QIIME 2 denoising-stats table from a .qza or an exported .tsv
# ---------------------------------------------------------------------------

def _load_stats_rows(stats: PathLike) -> Tuple[List[str], List[List[str]]]:
    """
    Return (header, data_rows) from a DADA2 denoising-stats table.

    Accepts either an exported stats .tsv or the denoising-stats .qza directly
    (the .qza is a zip containing <uuid>/data/stats.tsv). The QIIME '#q2:types'
    directive row is dropped. Raises on a missing file or unreadable table —
    a QC check given no data should not pass by default.
    """
    p = Path(stats)
    if not p.exists():
        raise FileNotFoundError(f"denoising-stats not found: {p}")

    if p.suffix == ".qza":
        with zipfile.ZipFile(p) as zf:
            members = [n for n in zf.namelist() if n.endswith("data/stats.tsv")]
            if not members:
                raise ValueError(f"no data/stats.tsv inside artifact: {p}")
            text = zf.read(members[0]).decode("utf-8")
    else:
        text = p.read_text(encoding="utf-8")

    reader = csv.reader(io.StringIO(text), delimiter="\t")
    rows = [r for r in reader if r and r[0].strip() != ""]
    if not rows:
        raise ValueError(f"denoising-stats table is empty: {p}")

    header = rows[0]
    data = [r for r in rows[1:] if not r[0].lower().startswith("#q2:types")]
    return header, data


def _find_col(header: List[str], *needles: str) -> int:
    """Return the index of the first column whose name contains all needles."""
    for i, name in enumerate(header):
        low = name.lower()
        if all(n in low for n in needles):
            return i
    raise ValueError(
        f"could not find a column matching {needles} in: {header}"
    )


# ---------------------------------------------------------------------------
# The retention gate
# ---------------------------------------------------------------------------

def check_retention(
    stats: PathLike,
    *,
    min_merged_frac: float = 0.60,
    min_nonchimeric_frac: float = 0.50,
    control_prefixes: Sequence[str] = (),
) -> QCResult:
    """
    Flag samples that lose too many reads through DADA2.

    The merge percentage is the most diagnostic: low merge retention almost
    always means --trunc-len was set so short that forward/reverse reads no
    longer overlap. Non-chimeric percentage catches excessive chimera removal.

    Negative controls (matched by control_prefixes, e.g. 'NTC-') are reported
    separately and never counted as failures — you *want* a control to lose its
    reads.

    Returns a QCResult with status:
        "pass" if no real sample falls below thresholds
        "warn" if one or more real samples do (the denoise step decides whether
               that warning is fatal, per config)
    """
    header, data = _load_stats_rows(stats)
    merged_i = _find_col(header, "percentage", "merged")
    nonchim_i = _find_col(header, "percentage", "non-chimeric")

    flagged: List[Tuple[str, float]] = []
    controls: List[Tuple[str, float]] = []
    n_samples = 0

    for row in data:
        sid = row[0]
        is_control = any(sid.startswith(c) for c in control_prefixes)
        try:
            merged = float(row[merged_i])
            nonchim = float(row[nonchim_i])
        except (ValueError, IndexError):
            log.warning("qc: could not parse retention for sample %r — skipping row", sid)
            continue

        if is_control:
            controls.append((sid, merged))
            continue

        n_samples += 1
        worst = min(merged, nonchim)
        if merged < min_merged_frac * 100 or nonchim < min_nonchimeric_frac * 100:
            flagged.append((sid, worst))

    flagged.sort(key=lambda t: t[1])

    if flagged:
        summary = (
            f"{len(flagged)}/{n_samples} samples below retention thresholds "
            f"(merge ≥ {min_merged_frac:.0%}, non-chimeric ≥ {min_nonchimeric_frac:.0%}); "
            f"worst {flagged[0][0]} at {flagged[0][1]:.1f}% — check --trunc-len"
        )
        status = "warn"
    else:
        summary = (
            f"all {n_samples} samples meet retention thresholds "
            f"(merge ≥ {min_merged_frac:.0%}, non-chimeric ≥ {min_nonchimeric_frac:.0%})"
        )
        status = "pass"

    return QCResult(
        check="retention",
        status=status,
        summary=summary,
        flagged=flagged,
        controls=controls,
    )


# ---------------------------------------------------------------------------
# Loading a tab-separated table from inside a .qza (or a plain .tsv)
# ---------------------------------------------------------------------------

def _load_qza_table(path: PathLike, member_suffix: str) -> Tuple[List[str], List[List[str]]]:
    """
    Return (header, data_rows) for a tab-separated table, read either from a
    plain .tsv or from the named member inside a .qza zip (e.g.
    'data/taxonomy.tsv'). The QIIME '#q2:types' directive row, if present, is
    dropped. Raises on a missing file or unreadable/empty table.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"table not found: {p}")
    if p.suffix == ".qza":
        with zipfile.ZipFile(p) as zf:
            members = [n for n in zf.namelist() if n.endswith(member_suffix)]
            if not members:
                raise ValueError(f"no {member_suffix} inside artifact: {p}")
            text = zf.read(members[0]).decode("utf-8")
    else:
        text = p.read_text(encoding="utf-8")
    rows = [r for r in csv.reader(io.StringIO(text), delimiter="\t") if r and r[0].strip() != ""]
    if not rows:
        raise ValueError(f"table is empty: {p}")
    header = rows[0]
    data = [r for r in rows[1:] if not r[0].lower().startswith("#q2:types")]
    return header, data


# ---------------------------------------------------------------------------
# The classifier-coverage gate
# ---------------------------------------------------------------------------

def check_classifier_coverage(
    taxonomy: PathLike,
    *,
    max_unassigned_frac: float = 0.20,
    unassigned_label: str = "Unassigned",
) -> QCResult:
    """
    Flag a classifier that left too many features unassigned.

    Reads the QIIME 2 taxonomy table (the Taxon column of data/taxonomy.tsv,
    from a .qza or an exported .tsv) and computes the fraction of features whose
    taxon is blank or begins with the unassigned label. A high fraction means
    the reference/classifier doesn't cover these sequences — the cue to BLAST
    the unassigned features (the 07b/07c/07d verification path) rather than
    trust the assignment.

    Returns a QCResult: "warn" when the unassigned fraction exceeds the
    threshold, else "pass". Non-fatal by design — the caller decides policy.
    """
    header, data = _load_qza_table(taxonomy, "data/taxonomy.tsv")
    taxon_i = _find_col(header, "taxon")

    total = 0
    unassigned = 0
    for row in data:
        total += 1
        taxon = row[taxon_i].strip() if taxon_i < len(row) else ""
        if (not taxon) or taxon.lower().startswith(unassigned_label.lower()):
            unassigned += 1

    frac = (unassigned / total) if total else 0.0
    status = "warn" if frac > max_unassigned_frac else "pass"
    summary = (
        f"{unassigned}/{total} features unassigned ({frac:.0%}); "
        f"threshold {max_unassigned_frac:.0%}"
        + ("  — consider BLAST verification of the unassigned features"
           if status == "warn" else "")
    )
    return QCResult(check="classifier_coverage", status=status, summary=summary)
