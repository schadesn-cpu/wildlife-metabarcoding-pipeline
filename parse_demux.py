#!/usr/bin/env python3
"""
parse_demux.py
==============
Step: parse_demux

Purpose:
    Parse Illumina demultiplex reports and MultiQC FastQC data into a pre-DADA2
    QC report — without cutadapt logs or a QIIME 2 environment. Answers, per
    sample:
      1. Read count vs a minimum threshold (low-depth flag)
      2. Adapter detection rate (did cutadapt have something to trim?)
      3. Effective amplicon length (from where the adapter onsets in the read)
      4. Dimer signal (adapter appearing very early — probable primer-dimer)
      5. Quality drop-off position (informs DADA2 --trunc-len)
      6. Primer identity cross-check

Inputs:
    reports/
    ├── primers_detected.tsv                 (from the primer advisor)
    └── demultiplex/
        ├── Demultiplex_Stats.csv            (reads per sample)
        └── additional-reports/
            ├── Adapter_Metrics.csv          (% adapter bases per sample)
            ├── Adapter_Cycle_Metrics.csv    (per-cycle counts — dimer signal)
            └── multiqc_data/
                ├── fastqc_adapter_content_plot.txt
                └── fastqc_per_base_sequence_quality_plot.txt
    pipeline_config.yml   active_markers, samples.control_prefixes, and optional
                          markers.<m>.amplicon_length / qc.min_reads

Outputs:
    results/qc/demux_qc_report.txt       (full formatted report)
    results/qc/demux_read_counts.tsv     (per-sample summary)
    logs/run_manifest.jsonl              (run appended on completion)

Markers are inferred from sample names using active_markers from config.
Exit code is non-zero if any FAIL-level (❌) issues are detected.

Usage:
    python parse_demux.py [--reports-dir reports/] [--min-reads 10000]
    python parse_demux.py --amplicon-lens 16S=253,MiFish=180,cytb=307

Requirements:
    Python >= 3.8, pandas. No QIIME 2 needed.
"""
from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# --- make config_loader and the utils package importable regardless of cwd ---
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from config_loader import load_config, get_paths            # noqa: E402
from utils import checkpoint, provenance                    # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Marker inference ──────────────────────────────────────────────────────────
# Default token→marker map. Overridden at runtime by tokens built from
# active_markers in config (see build_marker_tokens), so this is only a
# fallback for standalone use without a config.
MARKER_TOKENS = [
    ("MiFish",  "MiFish"),
    ("Mifish",  "MiFish"),
    ("mifish",  "MiFish"),
    ("cytb",    "cytb"),
    ("CYTB",    "cytb"),
    ("ITS1-2",  "ITS"),
    ("ITS",     "ITS"),
    ("18S",     "18S"),
    ("16S",     "16S"),
    ("Virus",   "Virus"),
    ("virus",   "Virus"),
]

# ---------------------------------------------------------------------------
# Amplicon length defaults (bp, excluding primers)
# ---------------------------------------------------------------------------
# Standard expected amplicon sizes for common markers, used in section 3 to
# flag deviations between where the adapter actually appears in the read and
# where it should given the amplicon length. These are generic defaults — set
# markers.<m>.amplicon_length in config, or pass --amplicon-lens, to override.
#   16S 253 bp (V4 515F/806R) · MiFish 180 bp (12S MiFish-U) · cytb ~307 bp
#   18S ~130 bp (V9 1391F/EukBr) · ITS variable (rough midpoint, unreliable)
DEFAULT_AMPLICON_LENS: Dict[str, int] = {
    "16S":    253,
    "18S":    130,
    "ITS":    250,
    "MiFish": 180,
    "cytb":   307,
    "Virus":  110,
}

# Control-sample prefixes; overridden at runtime from samples.control_prefixes.
CONTROL_PREFIXES = ("NTC-", "PAC-", "XB-", "NTC_", "PAC_", "XB_")


def build_marker_tokens(active_markers: List[str]) -> List[Tuple[str, str]]:
    """
    Build a token→marker map from the project's active markers, covering common
    case variants, so marker inference works for any project without a hardcoded
    list. Longer marker names are checked first so e.g. 'MiFish' wins before a
    substring could match something shorter.
    """
    tokens: List[Tuple[str, str]] = []
    for m in sorted(active_markers, key=len, reverse=True):
        for variant in (m, m.lower(), m.upper(), m.capitalize()):
            if (variant, m) not in tokens:
                tokens.append((variant, m))
    return tokens


def build_control_prefixes(cfg_prefixes: List[str]) -> Tuple[str, ...]:
    """Expand configured control prefixes to cover both '-' and '_' separators."""
    out: List[str] = []
    for p in cfg_prefixes:
        base = p.rstrip("-_")
        for sep in ("-", "_"):
            if (base + sep) not in out:
                out.append(base + sep)
    return tuple(out)


def infer_marker(sample_name: str) -> str:
    """
    Infer the marker from a sample name using MARKER_TOKENS.
    Returns 'unknown' if no token matches.
    """
    for token, marker in MARKER_TOKENS:
        if token in sample_name:
            return marker
    return "unknown"


def is_control(sample_name: str) -> bool:
    """Return True if the sample name starts with a known control prefix."""
    return any(sample_name.startswith(p) for p in CONTROL_PREFIXES)


def is_r1(sample_name: str) -> bool:
    """Return True if the sample name refers to the R1 (forward) read."""
    return "_R1_" in sample_name or sample_name.endswith("_R1")


# ── File parsers ──────────────────────────────────────────────────────────────

def parse_demultiplex_stats(path: Path) -> Dict[str, int]:
    """
    Parse Demultiplex_Stats.csv and return {Sample_ID: total_reads}.

    The CSV has columns including Sample_ID and # Reads (or similar).
    We try several known column name variants across Illumina pipeline
    versions rather than hardcoding a single name.
    """
    df = pd.read_csv(path, dtype=str)
    df.columns = df.columns.str.strip()

    # Try known column name variants for sample ID
    id_col = next((c for c in df.columns
                   if c.lower() in ("sample_id", "sampleid", "sample id",
                                    "sample", "sample name")), None)
    if id_col is None:
        id_col = df.columns[0]
        log.warning("Could not identify sample ID column in %s; using '%s'",
                    path.name, id_col)

    # Try known column name variants for read count
    count_col = next((c for c in df.columns
                      if any(x in c.lower() for x in
                             ("# reads", "reads", "clusters", "count"))),
                     None)
    if count_col is None:
        log.warning("Could not identify read count column in %s", path.name)
        return {}

    counts: Dict[str, int] = {}
    for _, row in df.iterrows():
        sid = str(row[id_col]).strip()
        try:
            counts[sid] = int(str(row[count_col]).replace(",", "").strip())
        except ValueError:
            # Non-numeric cells (e.g. header repetitions, summary rows) are
            # expected in some Illumina CSV formats — skip silently at debug.
            log.debug("Skipping non-numeric read count for '%s'", sid)

    log.info("  Demultiplex_Stats: %d samples", len(counts))
    return counts


def parse_adapter_metrics(path: Path) -> Dict[str, float]:
    """
    Parse Adapter_Metrics.csv and return {Sample_ID_R: pct_adapter_bases}.

    The '% Adapter Bases' column gives the fraction of bases in that
    sample/read direction that were identified as adapter sequence by
    the Illumina DRAGEN/BCL2FASTQ demultiplexer. A high value confirms
    the library has adapter bleed-through (expected for short amplicons).

    Returns a dict keyed by "{Sample_ID}_R{ReadNumber}" so R1 and R2
    are kept separate.
    """
    df = pd.read_csv(path, dtype=str)
    df.columns = df.columns.str.strip()

    result: Dict[str, float] = {}
    for _, row in df.iterrows():
        sid    = str(row.get("Sample_ID", "")).strip()
        readn  = str(row.get("ReadNumber", "1")).strip()
        pct    = str(row.get("% Adapter Bases", "0")).strip()
        try:
            result[f"{sid}_R{readn}"] = float(pct)
        except ValueError:
            # Non-numeric pct cells (header rows, empty cells) — skip at debug.
            log.debug("Skipping non-numeric adapter pct for '%s'", sid)

    log.info("  Adapter_Metrics: %d sample×read entries", len(result))
    return result


def parse_adapter_cycle_metrics(path: Path,
                                dimer_threshold_bp: int = 40
                                ) -> Dict[str, float]:
    """
    Parse Adapter_Cycle_Metrics.csv and estimate the dimer rate per sample.

    The file gives the number of clusters where the adapter was detected
    at each cycle position. Reads with adapter detected at positions ≤
    dimer_threshold_bp are classified as probable dimers (primer–primer
    or adapter–primer), since a genuine amplicon insert would be longer.

    Returns {Sample_ID_R: dimer_rate} where dimer_rate is the fraction
    of adapter-containing reads whose adapter appeared at a very short
    position. This is a proxy — the true dimer rate requires the total
    read count for normalisation.

    Args:
        path:                 Path to Adapter_Cycle_Metrics.csv.
        dimer_threshold_bp:   Cycle positions ≤ this are flagged as dimers.
                              Default 40 bp — safe margin above the longest
                              primer in the dataset (cytb at 35 bp).
    """
    df = pd.read_csv(path, dtype=str)
    df.columns = df.columns.str.strip()

    # Accumulate per sample: total adapter clusters and dimer-position clusters
    total_adapter: Dict[str, int] = {}
    dimer_clusters: Dict[str, int] = {}

    for _, row in df.iterrows():
        sid   = str(row.get("Sample_ID", "")).strip()
        readn = str(row.get("ReadNumber", "1")).strip()
        key   = f"{sid}_R{readn}"

        try:
            cycle = int(float(str(row.get("Cycle", "999")).strip()))
            n     = int(float(str(row.get("NumClustersWithAdapterAtCycle",
                                          "0")).strip()))
        except ValueError:
            continue

        total_adapter[key] = total_adapter.get(key, 0) + n
        if cycle <= dimer_threshold_bp:
            dimer_clusters[key] = dimer_clusters.get(key, 0) + n

    result: Dict[str, float] = {}
    for key, total in total_adapter.items():
        if total > 0:
            result[key] = dimer_clusters.get(key, 0) / total
        else:
            result[key] = 0.0

    log.info("  Adapter_Cycle_Metrics: %d sample×read entries", len(result))
    return result


def _parse_tuple_row(row_str: str) -> List[Tuple[float, float]]:
    """
    Parse a MultiQC FastQC data row of the form:
        (pos1, val1)  (pos2, val2)  ...
    Returns a list of (position, value) tuples.

    MultiQC sometimes has truncated closing parens due to line wrapping
    in the terminal paste — we handle partial tuples gracefully.
    """
    # Match all (number, number) patterns, tolerating minor formatting issues
    pairs = re.findall(r'\((\d+(?:\.\d+)?),\s*([\d.]+(?:e[+-]?\d+)?)\)',
                       row_str)
    return [(float(p), float(v)) for p, v in pairs]


def parse_adapter_content(path: Path,
                          adapter_type: str = "illumina_universal_adapter",
                          onset_threshold_pct: float = 5.0,
                          ) -> Dict[str, dict]:
    """
    Parse fastqc_adapter_content_plot.txt and extract, per sample + direction:
      - onset_position:  first cycle where adapter % exceeds onset_threshold_pct
                         This estimates the effective insert (amplicon) length.
      - plateau_pct:     the maximum adapter % reached (proxy for detection rate)
      - is_r1:           True if R1 read direction

    Only rows matching adapter_type are parsed. The "illumina_universal_adapter"
    rows represent the Illumina TruSeq adapter that is ligated to the end of the
    insert — its onset position directly encodes the insert length.

    Args:
        path:               Path to fastqc_adapter_content_plot.txt.
        adapter_type:       Adapter type to extract (default: illumina_universal_adapter).
        onset_threshold_pct: % at which we call the adapter as "appeared". Default 5%.

    Returns:
        Dict keyed by clean sample name (without adapter type suffix), with
        sub-keys: onset_position, plateau_pct, is_r1.
    """
    result: Dict[str, dict] = {}

    with path.open(encoding="utf-8", errors="replace") as fh:
        # First line is the header (position column labels) — skip it
        header = fh.readline()

        for line in fh:
            line = line.strip()
            if not line:
                continue

            # Tab-split: first field is "SampleName - adapter_type"
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue

            label = parts[0].strip()

            # Only process the adapter type we care about
            if f" - {adapter_type}" not in label:
                continue

            sample_name = label.replace(f" - {adapter_type}", "").strip()
            r1 = is_r1(sample_name)

            tuples = _parse_tuple_row(parts[1])
            if not tuples:
                continue

            # Find onset position: first position where adapter % > threshold
            onset_pos = None
            for pos, pct in tuples:
                if pct >= onset_threshold_pct:
                    onset_pos = int(pos)
                    break

            plateau = max(v for _, v in tuples) if tuples else 0.0

            result[sample_name] = {
                "onset_position": onset_pos,    # None if adapter never reaches threshold
                "plateau_pct":    plateau,       # max % adapter in this read
                "is_r1":          r1,
            }

    log.info("  Adapter content: %d samples (%s)", len(result), adapter_type)
    return result


def parse_quality_profiles(path: Path,
                            q_threshold: float = 25.0,
                            window: int = 5,
                            ) -> Dict[str, Optional[int]]:
    """
    Parse fastqc_per_base_sequence_quality_plot.txt and return the quality
    dropoff position for each sample.

    The dropoff position is the first cycle where the median Phred score
    stays below q_threshold for window consecutive positions. This is the
    same sliding-window approach used by primer_advisor suggest, and it
    informs the --trunc-len-f / --trunc-len-r values for DADA2.

    Returns {sample_name: dropoff_position} where dropoff_position is
    None if quality never drops below the threshold.
    """
    result: Dict[str, Optional[int]] = {}

    with path.open(encoding="utf-8", errors="replace") as fh:
        fh.readline()  # skip header

        for line in fh:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue

            sample_name = parts[0].strip()
            tuples = _parse_tuple_row(parts[1])
            if not tuples:
                continue

            positions = [int(p) for p, _ in tuples]
            qualities = [q for _, q in tuples]

            # Sliding window: find first position where q < threshold for
            # `window` consecutive cycles.
            dropoff = None
            n = len(qualities)
            for i in range(n - window + 1):
                if all(q < q_threshold for q in qualities[i:i + window]):
                    dropoff = positions[i]
                    break

            result[sample_name] = dropoff

    log.info("  Quality profiles: %d samples", len(result))
    return result


def load_primers(path: Path) -> Dict[str, Tuple[int, int, float, float]]:
    """
    Load primers_detected.tsv and return per-marker primer info.

    Returns {marker: (fwd_len, rev_len, fwd_score_pct, rev_score_pct)}.
    The score is already a percentage string like '44%' — we strip the
    percent sign and convert to float.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    result: Dict[str, Tuple[int, int, float, float]] = {}

    for _, row in df.iterrows():
        marker = str(row.get("marker", "")).strip()
        if not marker:
            continue

        # Normalise marker name to match the rest of the pipeline
        for token, canonical in MARKER_TOKENS:
            if token.lower() == marker.lower():
                marker = canonical
                break

        try:
            fwd_len = int(str(row.get("fwd_len", "0")).strip())
            rev_len = int(str(row.get("rev_len", "0")).strip())
            fwd_score = float(str(row.get("fwd_score", "0%")).strip()
                              .replace("%", ""))
            rev_score = float(str(row.get("rev_score", "0%")).strip()
                              .replace("%", ""))
        except ValueError:
            continue

        result[marker] = (fwd_len, rev_len, fwd_score, rev_score)

    log.info("  Primers loaded: %s", list(result.keys()))
    return result


# ── Report assembly ───────────────────────────────────────────────────────────

PASS  = "PASS"
WARN  = "WARN"
FAIL  = "FAIL"
INFO  = "INFO"

SYMBOLS = {PASS: "✅", WARN: "⚠️ ", FAIL: "❌", INFO: "ℹ️ "}


def flag(value: float, pass_t: float, warn_t: float,
         higher_is_better: bool = True) -> str:
    if higher_is_better:
        return PASS if value >= pass_t else (WARN if value >= warn_t else FAIL)
    else:
        return PASS if value <= pass_t else (WARN if value <= warn_t else FAIL)


def build_report(
    demux_counts:      Dict[str, int],
    adapter_metrics:   Dict[str, float],
    adapter_cycle:     Dict[str, float],
    adapter_content:   Dict[str, dict],
    quality_profiles:  Dict[str, Optional[int]],
    primers:           Dict[str, Tuple[int, int, float, float]],
    amplicon_lens:     Dict[str, int],
    min_reads:         int,
) -> Tuple[str, pd.DataFrame]:
    """
    Assemble the full QC report.

    Returns (report_text, summary_dataframe).
    The dataframe has one row per sample (R1 only for per-read metrics)
    with all QC fields for export to TSV.
    """
    SEP  = "=" * 72
    DASH = "─" * 72

    lines: List[str] = []
    rows:  List[dict] = []
    issues: List[str] = []

    lines += [
        SEP,
        "  PRE-DADA2 QC REPORT  —  MultiQC / Illumina Demux",
        SEP,
        "",
    ]

    # ── Section 1: Primer detection summary ───────────────────────────────
    lines.append(f"{DASH}")
    lines.append("  1. Primer Detection Summary  (from primers_detected.tsv)")
    lines.append(DASH)

    if not primers:
        lines.append(f"  {SYMBOLS[INFO]} No primers_detected.tsv provided.")
    else:
        lines.append(f"  {'Marker':<12} {'Fwd primer':<12} {'Fwd len':>8} "
                     f"{'Fwd score':>10} {'Rev primer':<12} {'Rev len':>8} "
                     f"{'Rev score':>10}  Status")
        lines.append(f"  {'-'*12} {'-'*12} {'-'*8} {'-'*10} {'-'*12} "
                     f"{'-'*8} {'-'*10}  ------")

        for marker, (fl, rl, fs, rs) in sorted(primers.items()):
            # Flag: both forward AND reverse should be >30% for confidence.
            # 18S at 0%/0% is a clear FAIL. cytb at 100%/98% is a clear PASS.
            min_score = min(fs, rs)
            status = flag(min_score, 30, 10)

            # Special warning: if fwd is good but rev is 0, likely wrong rev primer
            if fs > 30 and rs < 5:
                status = WARN
                issues.append(f"{marker}: forward primer detected but reverse is 0%")

            sym = SYMBOLS[status]
            lines.append(
                f"  {marker:<12} {'':<12} {fl:>8} {fs:>9.0f}% "
                f"{'':<12} {rl:>8} {rs:>9.0f}%  {sym} {status}"
            )

            if status == FAIL:
                issues.append(f"{marker}: low primer detection score "
                               f"(fwd={fs:.0f}%, rev={rs:.0f}%)")
    lines.append("")

    # ── Section 2: Per-sample read counts ────────────────────────────────
    lines.append(DASH)
    lines.append(f"  2. Per-Sample Read Counts  (threshold: {min_reads:,})")
    lines.append(DASH)

    all_samples = sorted(demux_counts.keys())
    n_pass = n_warn = n_fail = n_ctrl = 0

    lines.append(f"  {'Sample':<42} {'Reads':>10}  {'Control':>7}  Status")
    lines.append(f"  {'-'*42} {'-'*10}  {'-'*7}  ------")

    for sid in all_samples:
        count   = demux_counts[sid]
        ctrl    = is_control(sid)
        marker  = infer_marker(sid)

        if ctrl:
            status = INFO
            n_ctrl += 1
        elif count >= min_reads:
            status = PASS
            n_pass += 1
        elif count >= min_reads * 0.5:
            status = WARN
            n_warn += 1
        else:
            status = FAIL
            n_fail += 1

        sym = SYMBOLS[status]
        lines.append(
            f"  {sid:<42} {count:>10,}  {'yes' if ctrl else 'no':>7}  "
            f"{sym} {status}"
        )

        rows.append({
            "sample_id":     sid,
            "marker":        marker,
            "is_control":    ctrl,
            "read_count":    count,
            "read_count_ok": count >= min_reads,
        })

    lines.append(
        f"\n  {n_pass} PASS  {n_warn} WARN  {n_fail} FAIL  "
        f"{n_ctrl} controls (not counted)"
    )
    if n_fail > 0:
        issues.append(f"{n_fail} sample(s) below {min_reads:,} read threshold")
    lines.append("")

    # ── Section 3: Adapter content → amplicon length ─────────────────────
    lines.append(DASH)
    lines.append("  3. Effective Amplicon Length  (adapter onset position in R1)")
    lines.append(DASH)
    lines.append(
        "  The position where the Illumina adapter first appears in a read\n"
        "  equals the insert length. For a 110 bp virus fragment at 250 bp\n"
        "  sequencing, expect onset at ~110. For 16S V4 (~253 bp amplicon +\n"
        "  39 bp primers = 292 bp total), expect NO adapter in 250 bp reads.\n"
        "  Very short onset (<40 bp) = probable primer/adapter dimer.\n"
    )

    lines.append(f"  {'Sample':<42} {'Direction':>9} {'Onset (bp)':>11} "
                 f"{'Plateau%':>9} {'Expected':>9}  Status")
    lines.append(f"  {'-'*42} {'-'*9} {'-'*11} {'-'*9} {'-'*9}  ------")

    for sample_name, info in sorted(adapter_content.items()):
        if not info["is_r1"]:
            continue   # show R1 only in this table; R2 in quality section

        marker      = infer_marker(sample_name)
        ctrl        = is_control(sample_name)
        fwd_primer_len = primers.get(marker, (0, 0, 0, 0))[0]
        amp_len     = amplicon_lens.get(marker)
        expected_onset = amp_len if amp_len else None

        onset   = info["onset_position"]
        plateau = info["plateau_pct"]
        onset_str = str(onset) if onset is not None else "none"

        # Determine status
        if onset is None:
            # No adapter detected — fine for long amplicons (cytb, 16S),
            # unexpected for short ones (Virus, MiFish)
            if amp_len and amp_len < 200:
                status = WARN
                issues.append(f"{sample_name}: no adapter detected for "
                               f"short amplicon ({amp_len} bp)")
            else:
                status = PASS  # expected for 16S / cytb
        elif onset < 40:
            # Very short onset = dimer
            status = FAIL
            issues.append(f"{sample_name}: adapter at position {onset} bp "
                           "(probable dimer)")
        elif expected_onset and abs(onset - expected_onset) > 20:
            status = WARN
            issues.append(f"{sample_name}: onset {onset} bp vs expected "
                           f"{expected_onset} bp (diff {abs(onset-expected_onset)} bp)")
        else:
            status = PASS

        # Controls with adapter at very short positions are expected
        if ctrl and onset and onset < 40:
            status = INFO  # expected for NTC dimers

        sym = SYMBOLS[status]
        expected_str = str(expected_onset) if expected_onset else "N/A"
        lines.append(
            f"  {sample_name:<42} {'R1':>9} {onset_str:>11} "
            f"{plateau:>8.1f}% {expected_str:>9}  {sym} {status}"
        )

        # Update the row dict
        for row in rows:
            if row["sample_id"] == sample_name:
                row["adapter_onset_r1"]   = onset
                row["adapter_plateau_r1"] = round(plateau, 1)
                row["expected_amplicon"]  = expected_onset
                break
    lines.append("")

    # ── Section 4: Dimer rate from Adapter_Cycle_Metrics ─────────────────
    lines.append(DASH)
    lines.append("  4. Adapter/Primer Dimer Rate  (Adapter_Cycle_Metrics.csv)")
    lines.append(DASH)
    lines.append(
        "  Fraction of adapter-containing reads where adapter appeared at\n"
        "  position ≤40 bp. Reads this short are primer–primer or adapter–\n"
        "  primer dimers. >15% among adapter-containing reads is a concern.\n"
    )

    if not adapter_cycle:
        lines.append(f"  {SYMBOLS[INFO]} Adapter_Cycle_Metrics.csv not found or empty.")
    else:
        lines.append(f"  {'Sample (R1)':<42} {'Dimer rate':>11}  Status")
        lines.append(f"  {'-'*42} {'-'*11}  ------")

        for key in sorted(adapter_cycle.keys()):
            if "_R1" not in key:
                continue
            sample_name = key.replace("_R1", "")
            rate = adapter_cycle[key]
            ctrl = is_control(sample_name)

            if ctrl:
                status = INFO
            else:
                status = flag(rate, 0.05, 0.15, higher_is_better=False)
                if status == FAIL:
                    issues.append(f"{sample_name}: dimer rate {rate:.1%} > 15%")

            sym = SYMBOLS[status]
            lines.append(
                f"  {sample_name:<42} {rate:>10.1%}  {sym} {status}"
            )
    lines.append("")

    # ── Section 5: Quality dropoff → trunc-len guidance ──────────────────
    lines.append(DASH)
    lines.append("  5. Quality Dropoff → DADA2 Truncation Guidance")
    lines.append(DASH)
    lines.append(
        "  Position where median Phred quality drops below Q25 and stays\n"
        "  there for ≥5 consecutive cycles. Use as starting point for\n"
        "  --trunc-len-f (R1) and --trunc-len-r (R2) in DADA2. The final\n"
        "  values also need to guarantee sufficient paired-end overlap.\n"
    )

    if not quality_profiles:
        lines.append(f"  {SYMBOLS[INFO]} Quality profile data not found.")
    else:
        # Group R1 and R2 by clean sample name
        r1_samples: Dict[str, int] = {}
        r2_samples: Dict[str, int] = {}
        for sname, dropoff in quality_profiles.items():
            if dropoff is None:
                continue
            if is_r1(sname):
                clean = re.sub(r"_R1_\d+", "", sname).rstrip("_")
                r1_samples[clean] = dropoff
            else:
                clean = re.sub(r"_R2_\d+", "", sname).rstrip("_")
                r2_samples[clean] = dropoff

        # Summarise by marker
        marker_r1: Dict[str, List[int]] = {}
        marker_r2: Dict[str, List[int]] = {}
        for sname, dp in r1_samples.items():
            m = infer_marker(sname)
            marker_r1.setdefault(m, []).append(dp)
        for sname, dp in r2_samples.items():
            m = infer_marker(sname)
            marker_r2.setdefault(m, []).append(dp)

        all_markers = sorted(set(list(marker_r1.keys()) +
                                 list(marker_r2.keys())))
        lines.append(f"  {'Marker':<12} {'R1 Q25 dropoff':>16} "
                     f"{'R2 Q25 dropoff':>16}  Suggested trunc-len")
        lines.append(f"  {'-'*12} {'-'*16} {'-'*16}  -------------------")

        for m in all_markers:
            r1_vals = marker_r1.get(m, [])
            r2_vals = marker_r2.get(m, [])
            r1_med  = int(sorted(r1_vals)[len(r1_vals) // 2]) if r1_vals else None
            r2_med  = int(sorted(r2_vals)[len(r2_vals) // 2]) if r2_vals else None

            r1_str = str(r1_med) if r1_med else "N/A"
            r2_str = str(r2_med) if r2_med else "N/A"

            # Overlap check using amplicon length and primer lengths
            amp  = amplicon_lens.get(m)
            pf   = primers.get(m, (0, 0, 0, 0))
            pfl, prl = pf[0], pf[1]

            suggestion = ""
            if r1_med and r2_med and amp:
                overlap = (r1_med - pfl) + (r2_med - prl) - amp
                if overlap < 12:
                    suggestion = (f"⚠️  overlap={overlap} bp — too short, "
                                  f"consider longer trunc")
                elif overlap < 20:
                    suggestion = f"⚠️  overlap={overlap} bp — marginal"
                else:
                    suggestion = f"✅ overlap={overlap} bp — adequate"

            lines.append(
                f"  {m:<12} {r1_str:>16} {r2_str:>16}  {suggestion}"
            )
    lines.append("")

    # ── Section 6: Issues summary ─────────────────────────────────────────
    lines.append(DASH)
    lines.append("  6. Issues Summary")
    lines.append(DASH)
    if issues:
        lines.append(f"  {len(issues)} issue(s) detected:")
        for iss in issues:
            lines.append(f"    • {iss}")
    else:
        lines.append(f"  {SYMBOLS[PASS]} No issues detected.")
    lines.append("")
    lines.append(SEP)

    report_text = "\n".join(lines)
    df = pd.DataFrame(rows) if rows else pd.DataFrame()
    return report_text, df


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_amplicon_lens(s: str) -> Dict[str, int]:
    """
    Parse a comma-separated 'Marker=length' string into a dict.
    Example: '16S=253,MiFish=180,cytb=307'
    """
    result: Dict[str, int] = {}
    for pair in s.split(","):
        pair = pair.strip()
        if "=" not in pair:
            continue
        marker, length = pair.split("=", 1)
        try:
            result[marker.strip()] = int(length.strip())
        except ValueError:
            log.warning("Could not parse amplicon length: '%s'", pair)
    return result


def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser for parse_multiqc_demux.py."""
    p = argparse.ArgumentParser(
        prog="parse_multiqc_demux.py",
        description=(
            "Parse Illumina demultiplex reports and MultiQC FastQC data\n"
            "to produce a pre-DADA2 QC report. Requires the reports/ directory\n"
            "structure from the pipeline. No QIIME2 environment needed."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--reports-dir", default=None, type=Path,
        help="Path to reports/ directory (default: 'reports/' under project root).",
    )
    p.add_argument(
        "--config", default=None,
        help="Path to pipeline_config.yml.",
    )
    p.add_argument(
        "--primers", default=None, type=Path,
        help="Path to primers_detected.tsv (default: reports/primers_detected.tsv).",
    )
    p.add_argument(
        "--amplicon-lens", default=None,
        help=(
            "Override default amplicon lengths. Format: Marker=bp,Marker=bp ...\n"
            "Example: --amplicon-lens 16S=253,MiFish=180,cytb=307,Virus=110\n"
            f"Defaults: {', '.join(f'{k}={v}' for k,v in DEFAULT_AMPLICON_LENS.items())}"
        ),
    )
    p.add_argument(
        "--min-reads", type=int, default=None,
        help="Minimum acceptable read count per sample "
             "(default: qc.min_reads in config, else 10000).",
    )
    p.add_argument(
        "--out-txt", default=None, type=Path,
        help="Write report text to this path.",
    )
    p.add_argument(
        "--out-tsv", default=None, type=Path,
        help="Write per-sample summary TSV to this path.",
    )
    return p


def main(argv: Optional[List[str]] = None) -> int:
    """
    Load demultiplex + MultiQC data, run all QC checks, print and write the
    report, and append the run to the ledger. Returns non-zero if any
    FAIL-level (❌) issue is detected.
    """
    global MARKER_TOKENS, CONTROL_PREFIXES
    args = build_parser().parse_args(argv)

    cfg = load_config(args.config)
    paths = get_paths(cfg)

    # Marker inference and control detection come from config, not a hardcoded
    # list, so this works for any project's markers and naming scheme.
    MARKER_TOKENS = build_marker_tokens(cfg.active_markers)
    control_cfg = cfg.samples.get("control_prefixes", [])
    if control_cfg:
        CONTROL_PREFIXES = build_control_prefixes(control_cfg)

    reports_dir = (args.reports_dir or cfg.resolve("reports")).resolve()
    if not reports_dir.is_dir():
        log.error("--reports-dir not found: %s", reports_dir)
        return 1

    demux_dir = reports_dir / "demultiplex"
    addl_dir  = demux_dir / "additional-reports"
    mqc_dir   = addl_dir / "multiqc_data"

    demux_stats_path     = demux_dir / "Demultiplex_Stats.csv"
    adapter_metrics_path = addl_dir / "Adapter_Metrics.csv"
    adapter_cycle_path   = addl_dir / "Adapter_Cycle_Metrics.csv"
    adapter_content_path = mqc_dir / "fastqc_adapter_content_plot.txt"
    quality_path         = mqc_dir / "fastqc_per_base_sequence_quality_plot.txt"
    primers_path = args.primers or (reports_dir / "primers_detected.tsv")

    # Amplicon lengths: generic defaults <- per-marker config <- CLI override
    amplicon_lens = dict(DEFAULT_AMPLICON_LENS)
    for m in cfg.active_markers:
        cfg_len = cfg.markers.get(m, {}).get("amplicon_length")
        if cfg_len:
            amplicon_lens[m] = int(cfg_len)
    if args.amplicon_lens:
        amplicon_lens.update(parse_amplicon_lens(args.amplicon_lens))

    # min reads: CLI <- qc.min_reads in config <- 10000
    qc_cfg = cfg.qc if isinstance(getattr(cfg, "qc", {}), dict) else {}
    min_reads = args.min_reads if args.min_reads is not None else int(qc_cfg.get("min_reads", 10000))

    log.info("Loading data files from: %s", reports_dir)
    demux_counts = (parse_demultiplex_stats(demux_stats_path) if demux_stats_path.exists()
                    else (log.warning("Not found: %s", demux_stats_path) or {}))
    adapter_metrics = (parse_adapter_metrics(adapter_metrics_path) if adapter_metrics_path.exists()
                       else (log.warning("Not found: %s", adapter_metrics_path) or {}))
    adapter_cycle = (parse_adapter_cycle_metrics(adapter_cycle_path) if adapter_cycle_path.exists()
                     else (log.warning("Not found: %s", adapter_cycle_path) or {}))
    adapter_content = (parse_adapter_content(adapter_content_path) if adapter_content_path.exists()
                       else (log.warning("Not found: %s", adapter_content_path) or {}))
    quality_profiles = (parse_quality_profiles(quality_path) if quality_path.exists()
                        else (log.warning("Not found: %s", quality_path) or {}))
    primers = (load_primers(primers_path) if primers_path.exists()
               else (log.warning("Not found: %s", primers_path) or {}))

    report_text, summary_df = build_report(
        demux_counts=demux_counts, adapter_metrics=adapter_metrics,
        adapter_cycle=adapter_cycle, adapter_content=adapter_content,
        quality_profiles=quality_profiles, primers=primers,
        amplicon_lens=amplicon_lens, min_reads=min_reads,
    )
    print(report_text)

    # Default output locations via the PathBuilder (results/qc/), overridable.
    qc_dir  = paths.qc_dir()
    out_txt = args.out_txt or (qc_dir / "demux_qc_report.txt")
    out_tsv = args.out_tsv or (qc_dir / "demux_read_counts.tsv")

    out_txt.parent.mkdir(parents=True, exist_ok=True)
    out_txt.write_text(report_text, encoding="utf-8")
    log.info("Report written: %s", out_txt)
    written = [out_txt]
    if not summary_df.empty:
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(out_tsv, sep="\t", index=False)
        log.info("TSV written: %s", out_tsv)
        written.append(out_tsv)

    has_fail = any("❌" in line for line in report_text.splitlines())

    checkpoint.print_checkpoint(
        cfg, "parse_demux",
        produced=written,
        provenance={
            "outputs": written,
            "command": "python " + " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
            "extra": {"min_reads": min_reads, "qc_fail": has_fail},
        },
    )

    if has_fail:
        log.warning("QC flagged FAIL-level (❌) issues — review the report above.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
