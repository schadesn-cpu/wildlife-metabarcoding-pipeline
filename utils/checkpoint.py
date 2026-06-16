#!/usr/bin/env python3
"""
utils/checkpoint.py
===================
The standard footer every pipeline script prints when it finishes.

It answers the three questions a user always has at a checkpoint:
    - what did this step produce?
    - is there a QZV I should look at, and what am I looking for?
    - what do I run next?

The "next step" is read from utils/steps.py, NOT hardcoded here or in the
calling script — so reordering the pipeline can never leave a stale "now run
script 07" footer behind.

If a `provenance` dict is passed, the run is also appended to the ledger
(logs/run_manifest.jsonl) in the same call, so logging the run and telling the
user what's next are one action, not two things to remember.

Public API:
    print_checkpoint(cfg, step_key, produced=None, provenance=None)
    qzv_view_hint()      -> str   the standard "how to open a QZV" line
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from utils import steps as _steps
from utils import provenance as _prov

PathLike = Union[str, Path]

_BAR = "─" * 64


def qzv_view_hint() -> str:
    """The standard instruction for opening a QIIME 2 visualization."""
    return ("qiime tools view <file>.qzv   "
            "or drag onto https://view.qiime2.org  "
            "(renders in your browser, nothing uploaded)")


def _fmt(template: str, marker: Optional[str]) -> str:
    """Fill the {marker} placeholder used in registry artifact templates."""
    return template.replace("{marker}", marker) if marker else template


def print_checkpoint(
    cfg,
    step_key: str,
    *,
    marker: Optional[str] = None,
    produced: Optional[List[PathLike]] = None,
    provenance: Optional[Dict[str, Any]] = None,
    stream=None,
) -> None:
    """
    Print the checkpoint footer for a completed step.

    Parameters
    ----------
    step_key    key in utils/steps.py for the step that just finished
    marker      fills {marker} in the registry's artifact templates
    produced    actual output paths written this run; falls back to the
                registry's generic `produces` list if not given
    provenance  if provided, a dict of kwargs forwarded to
                provenance.record_run(cfg, step_key, **provenance) so the run
                is logged to the ledger as part of printing the checkpoint
    """
    out = stream or sys.stdout
    step = _steps.get_step(step_key)
    pos, total = _steps.step_index(step_key)

    # Log to the ledger first, so a printed checkpoint always corresponds to a
    # recorded run. The display marker is forwarded so every run is attributable
    # without each caller having to repeat it in the provenance dict.
    if provenance is not None:
        prov = dict(provenance)
        if marker is not None:
            prov.setdefault("marker", marker)
        _prov.record_run(cfg, step_key, **prov)

    lines: List[str] = []
    lines.append("")
    lines.append(_BAR)
    lines.append(f"  CHECKPOINT {pos}/{total} — {step.title}  ✓")
    lines.append(_BAR)

    # Produced
    items = [str(p) for p in produced] if produced else [
        _fmt(t, marker) for t in step.produces
    ]
    if items:
        lines.append("  Produced:")
        for it in items:
            lines.append(f"    - {it}")

    # Inspect
    if step.inspect:
        lines.append("  Inspect:")
        for ins in step.inspect:
            lines.append(f"    - {_fmt(ins.artifact, marker)}")
            lines.append(f"        look for: {ins.look_for}")
            if ins.is_qzv:
                lines.append(f"        open with: {qzv_view_hint()}")

    # Next
    nxt = _steps.next_step(step_key)
    if nxt is None:
        lines.append("  Next:")
        lines.append("    - pipeline complete — nothing further in the registry.")
    else:
        lines.append("  Next:")
        lines.append(f"    - {nxt.title}:  {_steps.run_command(nxt.key)}")
        if nxt.requires:
            lines.append(f"        heads-up: {nxt.requires}")
        if nxt.status == "planned":
            lines.append("        (this step is not yet ported to the clean "
                         "backbone — see legacy script)")

    lines.append(_BAR)
    lines.append("")

    out.write("\n".join(lines) + "\n")
    out.flush()
