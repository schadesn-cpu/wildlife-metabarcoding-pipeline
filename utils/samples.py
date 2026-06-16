#!/usr/bin/env python3
"""
utils/samples.py
================
Sample-id extraction, done once and made to fail loud.

QIIME 2 feature tables carry long ids like "TV230084-GI-16S_S1492", while a
source sample sheet keys on the short biological id "TV230084". Every step that
needed to bridge the two used to re-implement `str.extract(r"(TV\\d+)")` inline,
which silently returns nothing for an id that doesn't match — so a mistyped
regex or an off-scheme id quietly dropped a sample into a blank, unlabeled row.

This is the one config-driven extractor, and it RAISES on a non-matching id
instead of losing it. Controls are expected not to match and are handled
explicitly.

Public API:
    extract_sample_ids(ids, cfg, *, regex=None, control_prefixes=None)
        -> dict[str, str | None]
"""

from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional, Sequence


class SampleIDError(Exception):
    """Raised when a non-control sample id cannot be reduced to a key."""


def extract_sample_ids(
    ids: Iterable[str],
    cfg,
    *,
    regex: Optional[str] = None,
    control_prefixes: Optional[Sequence[str]] = None,
) -> Dict[str, Optional[str]]:
    """
    Map each sample id to its extracted biological key.

    Parameters
    ----------
    ids               the (long) sample ids to reduce
    cfg               loaded PipelineConfig — supplies samples.id_regex and
                      samples.control_prefixes by default
    regex             optional override for samples.id_regex (one capture group)
    control_prefixes  optional override for samples.control_prefixes

    Returns
    -------
    dict mapping each id to its key. Controls map to None (no biological key —
    they're expected not to match). Every non-control id MUST match; if any do
    not, SampleIDError is raised listing them, rather than returning a blank
    that would become an unlabeled sample downstream.
    """
    pattern_str = regex if regex is not None else cfg.samples.get("id_regex")
    if not pattern_str:
        raise SampleIDError(
            "No samples.id_regex set in config — cannot extract sample keys.\n"
            "  Add e.g.  id_regex: \"(TV\\\\d+)\"  under the 'samples' block, or "
            "pass regex=..."
        )
    pattern = re.compile(pattern_str)
    controls = tuple(
        control_prefixes if control_prefixes is not None
        else cfg.samples.get("control_prefixes", [])
    )

    keys: Dict[str, Optional[str]] = {}
    failures: List[str] = []
    for sid in ids:
        if controls and any(sid.startswith(c) for c in controls):
            keys[sid] = None
            continue
        m = pattern.search(sid)
        if not m:
            failures.append(sid)
            continue
        keys[sid] = m.group(1) if m.groups() else m.group(0)

    if failures:
        shown = "\n".join("  " + f for f in failures[:20])
        more = f"\n  ... and {len(failures) - 20} more" if len(failures) > 20 else ""
        raise SampleIDError(
            f"{len(failures)} sample id(s) did not match id_regex {pattern_str!r}:\n"
            f"{shown}{more}\n"
            "  Fix samples.id_regex in your config to match your naming, or check "
            "whether these ids belong in this run."
        )
    return keys
