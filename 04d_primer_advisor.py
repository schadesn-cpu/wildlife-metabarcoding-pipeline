#!/usr/bin/env python3
"""
primer_advisor.py
=================
Combined primer detection, DADA2 parameter suggestion, and post-cutadapt
QC tool. Replaces detect_primers.py and suggest_dada2_params.py.

Subcommands
-----------
detect  — Identify primers from raw R1/R2 reads and recommend a cutadapt command.
suggest — Suggest DADA2 truncation parameters from a QIIME2 demux QZV.
run     — Chain both: detect primers automatically, then suggest DADA2 params.
check   — Verify cutadapt worked: read counts per sample, adapter detection rate,
          dimer rate, and adapter position histogram (confirms amplicon length).

Workflow position — run these IN ORDER before touching DADA2
------------------------------------------------------------
Step 1  BEFORE cutadapt, directly on raw FASTQ reads:
    python primer_advisor.py detect \\
        --marker MiFish --reads-dir reads/

Step 2  AFTER qiime tools import, on the demux QZV:
    python primer_advisor.py suggest \\
        --demux qiime2/MiFish/all/imported/demux.qzv \\
        --primer-f 21 --primer-r 27 \\
        --amplicon-length 180 --marker MiFish

Step 3  AFTER cutadapt, on the trimmed QZV + pipeline log:
    python primer_advisor.py check \\
        --marker        MiFish \\
        --raw-qzv       qiime2/MiFish/all/imported/demux.qzv \\
        --trimmed-qzv   qiime2/MiFish/all/imported/demux_trimmed.qzv \\
        --cutadapt-log  logs/MiFish/all/cutadapt_*.log \\
        --amplicon-len  180 \\
        --min-reads     10000

    (run chains detect + suggest — it does NOT include the check step)

Full automatic detect + suggest in one call:
    python primer_advisor.py run \\
        --marker MiFish \\
        --reads-dir reads/ \\
        --demux qiime2/MiFish/all/imported/demux.qzv \\
        --amplicon-length 180

Detect all markers and write a TSV report:
    python primer_advisor.py detect \\
        --all --reads-dir reads/ --report primers_detected.tsv

Dependencies
------------
  - Python >= 3.8 (standard library only — no third-party packages required)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import sys
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  [%(levelname)s]  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ===========================================================================
# SECTION 1 — Primer database and IUPAC matching
# ===========================================================================

def rc(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


PRIMER_DB = {
    # ── 16S ──────────────────────────────────────────────────────────────────
    "515F":       {"seq": "GTGYCAGCMGCCGCGGTAA",              "marker": "16S",    "region": "V4",    "pair": "806R"},
    "806R":       {"seq": "GGACTACNVGGGTWTCTAAT",             "marker": "16S",    "region": "V4",    "pair": "515F"},
    "515F_EMP":   {"seq": "GTGCCAGCMGCCGCGGTAA",              "marker": "16S",    "region": "V4",    "pair": "806R"},
    "27F":        {"seq": "AGAGTTTGATCMTGGCTCAG",             "marker": "16S",    "region": "V1-V3", "pair": "519R"},
    "519R":       {"seq": "GWATTACCGCGGCKGCTG",               "marker": "16S",    "region": "V1-V3", "pair": "27F"},
    "341F":       {"seq": "CCTACGGGNGGCWGCAG",                "marker": "16S",    "region": "V3-V4", "pair": "806R"},

    # ── 18S ──────────────────────────────────────────────────────────────────
    "1391F":      {"seq": "GTACACACCGCCCGTC",                 "marker": "18S",    "region": "V9",    "pair": "EukBr"},
    "EukBr":      {"seq": "TGATCCTTCTGCAGGTTCACCTAC",         "marker": "18S",    "region": "V9",    "pair": "1391F"},
    "565F":       {"seq": "CCAGCASCYGCGGTAATTCC",             "marker": "18S",    "region": "V4",    "pair": "981R"},
    "981R":       {"seq": "ACTTTCGTTCTTGATYRA",               "marker": "18S",    "region": "V4",    "pair": "565F"},
    "EukA":       {"seq": "AACCTGGTTGATCCTGCCAGT",            "marker": "18S",    "region": "V1-V2", "pair": "EukB"},
    "EukB":       {"seq": "TGATCCTTCTGCAGGTTCACCTAC",         "marker": "18S",    "region": "V1-V2", "pair": "EukA"},

    # ── ITS ──────────────────────────────────────────────────────────────────
    "ITS1F":      {"seq": "CTTGGTCATTTAGAGGAAGTAA",           "marker": "ITS",    "region": "ITS1",  "pair": "ITS2"},
    "ITS2":       {"seq": "GCTGCGTTCTTCATCGATGC",             "marker": "ITS",    "region": "ITS1",  "pair": "ITS1F"},
    "ITS3":       {"seq": "GCATCGATGAAGAACGCAGC",             "marker": "ITS",    "region": "ITS2",  "pair": "ITS4"},
    "ITS4":       {"seq": "TCCTCCGCTTATTGATATGC",             "marker": "ITS",    "region": "ITS2",  "pair": "ITS3"},
    "ITS1":       {"seq": "TCCGTAGGTGAACCTGCGG",              "marker": "ITS",    "region": "ITS1",  "pair": "ITS4"},
    "fITS7":      {"seq": "GTGARTCATCGAATCTTTG",              "marker": "ITS",    "region": "ITS2",  "pair": "ITS4"},
    "ITS86F":     {"seq": "GTGAATCATCGAATCTTTGAA",            "marker": "ITS",    "region": "ITS2",  "pair": "ITS4"},

    # ── MiFish (12S) ─────────────────────────────────────────────────────────
    "MiFish-U-F": {"seq": "GTCGGTAAAACTCGTGCCAGC",            "marker": "MiFish", "region": "12S",   "pair": "MiFish-U-R"},
    "MiFish-U-R": {"seq": "CATAGTGGGGTATCTAATCCCAGTTTG",      "marker": "MiFish", "region": "12S",   "pair": "MiFish-U-F"},
    "MiFish-E-F": {"seq": "GTCGGTAAAACTCGTGCCAGC",            "marker": "MiFish", "region": "12S",   "pair": "MiFish-E-R"},
    "MiFish-E-R": {"seq": "CATAGTGGGGTATCTAATCCCAGTTTG",      "marker": "MiFish", "region": "12S",   "pair": "MiFish-E-F"},

    # ── cytb ─────────────────────────────────────────────────────────────────
    "L14841":     {"seq": "AAAAAGCTTCCATCCAACATCTCAGCATGATGAAA", "marker": "cytb", "region": "cytb", "pair": "H15149"},
    "H15149":     {"seq": "AAACTGCAGCCCCTCAGAATGATATTTGTCCTCA",  "marker": "cytb", "region": "cytb", "pair": "L14841"},
    "cytb-F":     {"seq": "CGAAACTTGATCAACGAACC",              "marker": "cytb",   "region": "cytb",  "pair": "cytb-R"},
    "cytb-R":     {"seq": "GGGTTGTTTGATCCTGTTTCG",             "marker": "cytb",   "region": "cytb",  "pair": "cytb-F"},
}

# IUPAC degenerate base expansion — maps each ambiguity code to the set of
# bases it represents.  iupac_match() uses this dict directly; there is no
# parallel if/elif chain.
IUPAC_BASES: Dict[str, set] = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'},
}


def iupac_match(primer: str, read: str, mismatches: int = 2) -> bool:
    """
    Return True if primer matches the 5' start of read, allowing IUPAC
    degenerate bases and up to `mismatches` mismatches.

    Uses IUPAC_BASES for all ambiguity resolution — no hard-coded per-code
    conditionals.
    """
    primer = primer.upper()
    read = read.upper()
    check_len = min(len(primer), len(read))
    if check_len < 8:
        return False
    errors = 0
    for p, r in zip(primer[:check_len], read[:check_len]):
        allowed = IUPAC_BASES.get(p, {p})   # unknown codes match only themselves
        if r not in allowed:
            errors += 1
            if errors > mismatches:
                return False
    return True


def score_primer(primer_seq: str, reads: List[str], mismatches: int = 2,
                 prefix_len: Optional[int] = None) -> float:
    """
    Return the fraction of reads that match primer_seq at their 5' end.

    prefix_len: if set, only check the first N bases — useful for
    high-degeneracy primers like ITS1F where the 3' end is variable.
    """
    if not reads:
        return 0.0
    check_seq = primer_seq[:prefix_len] if prefix_len else primer_seq
    matches = sum(1 for r in reads if iupac_match(check_seq, r, mismatches))
    return matches / len(reads)


# ===========================================================================
# SECTION 2 — Read sampling
# ===========================================================================

def get_sample_reads(reads_dir: Path, marker: str,
                     n: int = 200) -> Tuple[List[str], List[str]]:
    """
    Return up to n sequence strings from R1 and R2 FASTQ.gz files for a marker.

    Reads are sampled from the first file(s) found — sufficient for primer
    detection. Control samples (prefixes NTC-, PAC-, XB-) are skipped because
    they typically have low or no target DNA and would depress match scores.

    R1/R2 files are identified by '_R1_', '_R1.', or '_r1.fastq.gz' patterns
    in the filename (case-insensitive). The marker subdirectory is matched
    case-insensitively so 'mifish' and 'MiFish' both work.

    FASTQ format: every 4th line starting from line 2 (0-indexed) is the
    sequence. This function reads only those lines to avoid loading quality
    scores or headers into memory.

    Args:
        reads_dir: Root directory containing per-marker subdirectories.
        marker:    Marker name (e.g. 'MiFish', '16S'). Matched case-insensitively.
        n:         Maximum number of reads to return per direction.

    Returns:
        Tuple of (r1_sequences, r2_sequences) as lists of strings.

    Raises:
        FileNotFoundError: If no subdirectory matching marker is found.
    """
    marker_dir = reads_dir / marker
    if not marker_dir.exists():
        for d in reads_dir.iterdir():
            if d.name.lower() == marker.lower():
                marker_dir = d
                break
        else:
            raise FileNotFoundError(
                f"Reads directory not found for marker: {marker}\n"
                f"Expected: {reads_dir / marker}"
            )

    control_prefixes = ("NTC-", "PAC-", "XB-", "ntc-", "pac-", "xb-")
    r1_seqs: List[str] = []
    r2_seqs: List[str] = []

    for f in sorted(marker_dir.iterdir()):
        if f.name.startswith(".") or any(f.name.startswith(p) for p in control_prefixes):
            continue
        fname = f.name.lower()
        if "_r1_" in fname or "_r1." in fname or fname.endswith("_r1.fastq.gz"):
            try:
                with gzip.open(f, "rt") as fh:
                    for i, line in enumerate(fh):
                        if i % 4 == 1:
                            r1_seqs.append(line.strip())
                        if len(r1_seqs) >= n:
                            break
            except Exception:
                pass
            if len(r1_seqs) >= n:
                break

    for f in sorted(marker_dir.iterdir()):
        if f.name.startswith(".") or any(f.name.startswith(p) for p in control_prefixes):
            continue
        fname = f.name.lower()
        if "_r2_" in fname or "_r2." in fname or fname.endswith("_r2.fastq.gz"):
            try:
                with gzip.open(f, "rt") as fh:
                    for i, line in enumerate(fh):
                        if i % 4 == 1:
                            r2_seqs.append(line.strip())
                        if len(r2_seqs) >= n:
                            break
            except Exception:
                pass
            if len(r2_seqs) >= n:
                break

    return r1_seqs[:n], r2_seqs[:n]


# ===========================================================================
# SECTION 3 — Primer detection
# ===========================================================================

def detect_primers_for_marker(
    marker: str,
    reads_dir: Path,
    n_reads: int = 200,
    mismatch: int = 2,
    verbose: bool = False,
) -> dict:
    """
    Detect primers for a single marker by scoring all database entries
    against sampled R1 and R2 reads.

    Returns a dict with keys: marker, fwd_name, fwd_seq, fwd_len, fwd_score,
    rev_name, rev_seq, rev_len, rev_score, region, variable_length,
    spacer_detected, cutadapt_cmd, notes.
    """
    print(f"\n{'='*60}")
    print(f"  Marker: {marker}")
    print(f"{'='*60}")

    try:
        r1_seqs, r2_seqs = get_sample_reads(reads_dir, marker, n=n_reads)
    except FileNotFoundError as e:
        print(f"  [ERROR] {e}")
        return {"marker": marker, "error": str(e)}

    print(f"  Loaded {len(r1_seqs)} R1 reads, {len(r2_seqs)} R2 reads")

    r1_scores: Dict[str, float] = {}
    r2_scores: Dict[str, float] = {}

    for name, info in PRIMER_DB.items():
        seq = info["seq"]
        r1_scores[name] = max(
            score_primer(seq, r1_seqs, mismatch),
            score_primer(seq, r1_seqs, mismatch, prefix_len=15),
        )
        r2_scores[name] = max(
            score_primer(seq, r2_seqs, mismatch),
            score_primer(seq, r2_seqs, mismatch, prefix_len=15),
        )

    best_fwd = max(r1_scores, key=r1_scores.get)
    best_fwd_score = r1_scores[best_fwd]
    best_rev = max(r2_scores, key=r2_scores.get)
    best_rev_score = r2_scores[best_rev]

    if verbose:
        print("\n  Top R1 matches (forward primer candidates):")
        for name, score in sorted(r1_scores.items(), key=lambda x: -x[1])[:5]:
            if score > 0.01:
                print(f"    {name:<20} {score:.1%}  [{PRIMER_DB[name]['seq']}]")
        print("\n  Top R2 matches (reverse primer candidates):")
        for name, score in sorted(r2_scores.items(), key=lambda x: -x[1])[:5]:
            if score > 0.01:
                print(f"    {name:<20} {score:.1%}  [{PRIMER_DB[name]['seq']}]")

    # Prefer the known biological pair when the forward primer is identified
    # with reasonable confidence — handles cases where R2 primer score is low
    # due to short read coverage but the biology is unambiguous.
    fwd_info = PRIMER_DB[best_fwd]
    expected_rev = fwd_info.get("pair")
    if expected_rev and expected_rev in PRIMER_DB and best_fwd_score > 0.2:
        rev_name  = expected_rev
        rev_seq   = PRIMER_DB[expected_rev]["seq"]
        rev_score = r2_scores.get(expected_rev, 0)
        if verbose:
            print(f"\n  Using expected pair for {best_fwd}: {rev_name}")
    else:
        rev_name  = best_rev
        rev_seq   = PRIMER_DB[best_rev]["seq"]
        rev_score = best_rev_score

    fwd_seq = fwd_info["seq"]
    region  = fwd_info.get("region", "unknown")

    is_its = "ITS" in marker.upper()
    variable_length = is_its

    print(f"\n  Forward primer: {best_fwd} ({best_fwd_score:.0%} match)")
    print(f"    Sequence: {fwd_seq} ({len(fwd_seq)} bp)")
    print(f"  Reverse primer: {rev_name} ({rev_score:.0%} match)")
    print(f"    Sequence: {rev_seq} ({len(rev_seq)} bp)")
    print(f"  Region:   {region}")

    if best_fwd_score < 0.3:
        print(f"\n  [WARN] Low confidence match ({best_fwd_score:.0%}) — verify primers manually")

    notes = []
    if is_its:
        notes.append("ITS: variable-length amplicon — use --p-trunc-len-f 0 --p-trunc-len-r 0 in DADA2")
        notes.append("ITS: read-through may occur — cutadapt will also trim RC of reverse primer from 3' end")
        cutadapt_cmd = (
            f"qiime cutadapt trim-paired \\\n"
            f"    --i-demultiplexed-sequences qiime2/imported/demux_{marker}.qza \\\n"
            f"    --p-front-f {fwd_seq} \\\n"
            f"    --p-front-r {rev_seq} \\\n"
            f"    --p-adapter-f {rc(rev_seq)} \\\n"
            f"    --p-adapter-r {rc(fwd_seq)} \\\n"
            f"    --p-discard-untrimmed \\\n"
            f"    --p-error-rate 0.2 \\\n"
            f"    --o-trimmed-sequences qiime2/imported/demux_{marker}_trimmed.qza \\\n"
            f"    --verbose"
        )
    else:
        cutadapt_cmd = (
            f"qiime cutadapt trim-paired \\\n"
            f"    --i-demultiplexed-sequences qiime2/imported/demux_{marker}.qza \\\n"
            f"    --p-front-f {fwd_seq} \\\n"
            f"    --p-front-r {rev_seq} \\\n"
            f"    --p-discard-untrimmed \\\n"
            f"    --p-error-rate 0.1 \\\n"
            f"    --o-trimmed-sequences qiime2/imported/demux_{marker}_trimmed.qza \\\n"
            f"    --verbose"
        )

    spacer_detected = False
    if r1_seqs:
        unique_starts = len({r[:4] for r in r1_seqs[:20]})
        if unique_starts > 3:
            spacer_detected = True
            notes.append("Heterogeneity spacers detected — primer position varies; cutadapt handles this correctly")

    print(f"\n  Recommended cutadapt command:")
    for line in cutadapt_cmd.splitlines():
        print(f"    {line}")

    if notes:
        print(f"\n  Notes:")
        for note in notes:
            print(f"    • {note}")

    return {
        "marker":          marker,
        "fwd_name":        best_fwd,
        "fwd_seq":         fwd_seq,
        "fwd_len":         len(fwd_seq),
        "fwd_score":       f"{best_fwd_score:.0%}",
        "rev_name":        rev_name,
        "rev_seq":         rev_seq,
        "rev_len":         len(rev_seq),
        "rev_score":       f"{rev_score:.0%}",
        "region":          region,
        "variable_length": str(variable_length),
        "spacer_detected": str(spacer_detected),
        "cutadapt_cmd":    cutadapt_cmd,
        "notes":           "; ".join(notes),
    }


def write_detect_report(results: List[dict], out_path: Path) -> None:
    """
    Write primer detection results to a TSV file.

    Columns: marker, fwd_name, fwd_seq, fwd_len, fwd_score, rev_name, rev_seq,
    rev_len, rev_score, region, variable_length, spacer_detected, notes.

    Args:
        results:  List of result dicts from detect_primers_for_marker().
        out_path: Destination path for the TSV file.
    """
    cols = ["marker", "fwd_name", "fwd_seq", "fwd_len", "fwd_score",
            "rev_name", "rev_seq", "rev_len", "rev_score",
            "region", "variable_length", "spacer_detected", "notes"]
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    print(f"\nReport written to: {out_path}")


def discover_markers(reads_dir: Path) -> List[str]:
    """List all marker subdirectories in the reads directory."""
    return [d.name for d in sorted(reads_dir.iterdir())
            if d.is_dir() and not d.name.startswith(".")]


# ===========================================================================
# SECTION 4 — DADA2 parameter suggestion
# ===========================================================================

MEDIAN_COL        = "50%"
DEFAULT_MIN_QUAL  = 25
DEFAULT_MIN_OVER  = 20
DROPOFF_WINDOW    = 5


def find_quality_tsvs(zf: zipfile.ZipFile) -> Tuple[Optional[str], Optional[str]]:
    """
    Locate forward and reverse seven-number-summary TSVs inside a QIIME2 QZV.

    QZV archives have the structure:
        <uuid>/data/forward-seven-number-summaries.tsv
        <uuid>/data/reverse-seven-number-summaries.tsv   (paired only)
    """
    names = zf.namelist()
    fwd = next((n for n in names if "forward-seven-number-summaries" in n), None)
    rev = next((n for n in names if "reverse-seven-number-summaries" in n), None)
    return fwd, rev


def parse_quality_tsv(zf: zipfile.ZipFile,
                      tsv_path: str) -> Tuple[List[int], List[float]]:
    """
    Parse a QIIME2 seven-number-summary TSV and return (positions, medians).

    The TSV columns are: position  2%  9%  25%  50%  75%  91%  98%
    We use the 50th percentile (median) as the quality indicator.
    """
    content = zf.read(tsv_path).decode("utf-8")
    lines = content.strip().splitlines()
    header = lines[0].split("\t")

    try:
        median_idx = header.index(MEDIAN_COL)
    except ValueError:
        log.warning(
            "Could not find '%s' column in quality TSV; falling back to index 4.",
            MEDIAN_COL,
        )
        median_idx = 4

    positions: List[int] = []
    medians:   List[float] = []

    for line in lines[1:]:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        try:
            positions.append(int(float(parts[0])))
            medians.append(float(parts[median_idx]))
        except (ValueError, IndexError) as e:
            log.debug("Skipping malformed TSV line: %s (%s)", line[:60], e)

    if not positions:
        raise ValueError(f"No quality data parsed from {tsv_path}")

    return positions, medians


def find_quality_dropoff(positions: List[int], medians: List[float],
                         min_quality: float,
                         window: int = DROPOFF_WINDOW) -> int:
    """
    Find the cycle position where quality drops below min_quality and stays
    there for at least `window` consecutive positions.

    Uses a sliding window to avoid truncating aggressively on a single dip.
    Returns the last position if quality never drops below the threshold.
    """
    n = len(positions)
    for i in range(n - window + 1):
        if all(q < min_quality for q in medians[i: i + window]):
            return positions[i]
    return positions[-1]


def check_overlap(trunc_f: int, trunc_r: int,
                  primer_f: int, primer_r: int,
                  amplicon_length: int,
                  min_overlap: int) -> Tuple[int, bool, str]:
    """
    Calculate expected overlap between F and R reads after DADA2 truncation.

    Formula:  overlap = (trunc_f − primer_f) + (trunc_r − primer_r) − amplicon_length

    Returns (overlap_bp, can_merge, human_readable_message).
    """
    overlap = (trunc_f - primer_f) + (trunc_r - primer_r) - amplicon_length

    if overlap >= min_overlap:
        return overlap, True, (
            f"Overlap = {overlap} bp (≥ {min_overlap} bp minimum) — "
            f"paired-end merging should succeed."
        )
    elif overlap >= 0:
        return overlap, False, (
            f"Overlap = {overlap} bp (< {min_overlap} bp minimum) — "
            f"merging is marginal; consider single-end mode."
        )
    else:
        return overlap, False, (
            f"Overlap = {overlap} bp (NEGATIVE) — "
            f"reads cannot bridge the amplicon; single-end mode required."
        )


def _round_down(value: int, base: int = 5) -> int:
    """Round value down to the nearest multiple of base (e.g. 148 → 145)."""
    return (value // base) * base


def suggest_params(
    demux_path: Path,
    primer_f: int,
    primer_r: int,
    amplicon_length: Optional[int],
    min_quality: float,
    min_overlap: int,
    force_single_end: bool,
) -> dict:
    """
    Core suggestion logic: parse quality scores, find dropoffs, check overlap,
    and return a recommendation dict.

    Returns keys: mode, trim_left_f, trim_left_r, trunc_len_f, trunc_len_r,
    overlap_bp, can_merge, dropoff_f, dropoff_r, warnings, quality_profile.
    """
    result: dict = {
        "mode": None, "trim_left_f": primer_f, "trim_left_r": primer_r,
        "trunc_len_f": None, "trunc_len_r": None,
        "overlap_bp": None, "can_merge": None,
        "dropoff_f": None,  "dropoff_r": None,
        "warnings": [], "quality_profile": {},
    }

    if not demux_path.exists():
        raise FileNotFoundError(f"Demux file not found: {demux_path}")

    log.info("Opening demux archive: %s", demux_path)

    try:
        zf = zipfile.ZipFile(demux_path, "r")
    except zipfile.BadZipFile:
        raise ValueError(f"File does not appear to be a valid ZIP/QZV: {demux_path}")

    with zf:
        fwd_tsv, rev_tsv = find_quality_tsvs(zf)
        if fwd_tsv is None:
            raise ValueError(
                "Could not find forward quality summary TSV inside the QZV. "
                "Provide a demux.qzv (visualization), not a demux.qza (artifact)."
            )

        log.info("Parsing forward read quality profile...")
        pos_f, med_f = parse_quality_tsv(zf, fwd_tsv)
        result["quality_profile"]["forward"] = {"positions": pos_f, "medians": med_f}

        is_paired = (rev_tsv is not None) and not force_single_end
        if is_paired:
            log.info("Parsing reverse read quality profile...")
            pos_r, med_r = parse_quality_tsv(zf, rev_tsv)
            result["quality_profile"]["reverse"] = {"positions": pos_r, "medians": med_r}
        else:
            reason = "--single-end forced" if force_single_end else "no reverse data found"
            log.info("Single-end mode (%s).", reason)

    # Find dropoff positions
    log.info("Finding quality dropoff (Q%.0f, window %d bp)...", min_quality, DROPOFF_WINDOW)

    dropoff_f = find_quality_dropoff(pos_f, med_f, min_quality)
    trunc_f   = _round_down(dropoff_f)
    result["dropoff_f"]   = dropoff_f
    result["trunc_len_f"] = trunc_f
    log.info("Forward: Q%.0f dropoff at %d → trunc_len_f = %d", min_quality, dropoff_f, trunc_f)

    if is_paired:
        dropoff_r = find_quality_dropoff(pos_r, med_r, min_quality)
        trunc_r   = _round_down(dropoff_r)
        result["dropoff_r"]   = dropoff_r
        result["trunc_len_r"] = trunc_r
        log.info("Reverse: Q%.0f dropoff at %d → trunc_len_r = %d", min_quality, dropoff_r, trunc_r)
    else:
        trunc_r = None

    # Overlap check
    if is_paired and amplicon_length is not None:
        overlap, can_merge, msg = check_overlap(
            trunc_f, trunc_r, primer_f, primer_r, amplicon_length, min_overlap
        )
        result["overlap_bp"] = overlap
        result["can_merge"]  = can_merge
        result["mode"] = "paired" if can_merge else "single"
        (log.info if can_merge else log.warning)("Overlap check: %s", msg)
        if not can_merge:
            result["warnings"].append(msg)
    elif is_paired:
        result["mode"] = "paired"
        result["warnings"].append(
            "No --amplicon-length provided; cannot verify paired-end overlap. "
            "Check that trunc_len_f + trunc_len_r > amplicon_length + min_overlap."
        )
    else:
        result["mode"] = "single"

    # Sanity: warn on very short effective reads
    MIN_USABLE = 80
    for direction, tlen, plen in [("forward", trunc_f, primer_f),
                                   ("reverse", trunc_r, primer_r)]:
        if tlen is not None and (tlen - plen) < MIN_USABLE:
            w = (
                f"Effective {direction} read after primer trim is only {tlen - plen} bp "
                f"(trunc_len={tlen} − primer={plen}). May be too short for reliable taxonomy."
            )
            result["warnings"].append(w)
            log.warning(w)

    return result


def print_suggestion(result: dict, marker: str) -> None:
    """Print a human-readable DADA2 parameter recommendation."""
    mode, trunc_f, trunc_r = result["mode"], result["trunc_len_f"], result["trunc_len_r"]
    trim_f, trim_r = result["trim_left_f"], result["trim_left_r"]
    sep = "=" * 65

    print(f"\n{sep}")
    print(f"  DADA2 Parameter Recommendations — {marker}")
    print(sep)
    print(f"  Mode:              {mode.upper()}-END")
    print(f"  --trim-left-f:     {trim_f}  (forward primer length)")
    if mode == "paired":
        print(f"  --trim-left-r:     {trim_r}  (reverse primer length)")
    print(f"  --trunc-len-f:     {trunc_f}")
    if mode == "paired":
        print(f"  --trunc-len-r:     {trunc_r}")
    if result["overlap_bp"] is not None:
        print(f"  Expected overlap:  {result['overlap_bp']} bp")
    print()

    print("  ── Suggested DADA2 command fragment ──────────────────────")
    if mode == "paired":
        print(f"  qiime dada2 denoise-paired \\")
        print(f"      --p-trim-left-f {trim_f} --p-trim-left-r {trim_r} \\")
        print(f"      --p-trunc-len-f {trunc_f} --p-trunc-len-r {trunc_r}")
    else:
        print(f"  qiime dada2 denoise-single \\")
        print(f"      --p-trim-left {trim_f} \\")
        print(f"      --p-trunc-len {trunc_f}")

    print()
    print("  ── For 03_run_full_metabarcoding_pipeline.py ─────────────")
    if mode == "paired":
        print(f"  dada2 \\")
        print(f"      --trim-left-f {trim_f} --trim-left-r {trim_r} \\")
        print(f"      --trunc-len-f {trunc_f} --trunc-len-r {trunc_r}")
    else:
        print(f"  dada2 \\")
        print(f"      --trim-left-f {trim_f} \\")
        print(f"      --trunc-len-f {trunc_f}")

    if result["warnings"]:
        print()
        print("  ── Warnings ──────────────────────────────────────────────")
        for w in result["warnings"]:
            print(f"  ⚠  {w}")

    print(sep)
    print()


def write_suggest_json(result: dict, marker: str, out_path: Path) -> None:
    """Write parameter recommendations to JSON for downstream pipeline use."""
    out = {k: result[k] for k in (
        "mode", "trim_left_f", "trim_left_r", "trunc_len_f", "trunc_len_r",
        "overlap_bp", "can_merge", "dropoff_f", "dropoff_r", "warnings",
    )}
    out["marker"] = marker
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(out, f, indent=2)
    log.info("Parameter JSON written to: %s", out_path)



# ===========================================================================
# SECTION 5 — Post-cutadapt sanity check
# ===========================================================================
# These functions answer the bioinformatician's key QC questions AFTER cutadapt
# has run but BEFORE you commit to DADA2 parameters.  The three things we need
# to verify are:
#
#   1. READ COUNTS — did every sample keep ≥ min_reads after trimming?
#      (Your bioinformatician said ≥10 k; low-count samples should be flagged
#       before denoising, not discovered as empty tables afterward.)
#
#   2. ADAPTER DETECTION RATE — did cutadapt actually find our primers?
#      If <80% of reads had an adapter detected, either the primer sequences
#      are wrong or reads are genuinely shorter than the amplicon (no bleed-
#      through), so trimming was a no-op.
#
#   3. DIMER RATE — how many reads were too short after trimming?
#      These are primer–primer or adapter–primer dimers: the two primers
#      ligated to each other without a real template in between.  They
#      appear as a peak at very short positions in the adapter histogram.
#      A high dimer rate means real data was lost and gel cleanup may be needed.
#
#   4. ADAPTER POSITION HISTOGRAM — where in the read did the adapter appear?
#      The modal position equals the effective insert length (amplicon bp).
#      For a 110 bp virus fragment sequenced at 250 bp, you expect the peak
#      near position 110.  A peak nowhere near the expected length means
#      something went wrong in primer design or library prep.
#
# Note: quality profile / trunc-len guidance deliberately lives in 'suggest'
# (Section 4), not here.  This section is only about whether cutadapt itself
# did what we asked it to do.
# ---------------------------------------------------------------------------

import glob  # needed for expanding wildcard log paths (e.g. cutadapt_*.log)
import re    # needed for parsing cutadapt's plain-text summary lines


def _find_in_zip(zf: zipfile.ZipFile, filename: str) -> Optional[str]:
    """
    Return the internal ZIP path of the first entry whose basename matches
    `filename`, or None if not found.

    QIIME2 QZVs are standard ZIP files but their internal paths are prefixed
    by a UUID directory, e.g.:
        a1b2c3d4-e5f6-7890-abcd-ef1234567890/data/per-sample-fastq-counts.tsv

    We cannot know the UUID ahead of time, so we scan all entries and compare
    only the final path component.  There should never be more than one file
    with a given basename inside a single QZV.
    """
    for name in zf.namelist():
        if Path(name).name == filename:
            return name
    return None


def parse_per_sample_counts(qzv_path: Path) -> Dict[str, int]:
    """
    Extract {sample_id: forward_read_count} from a QIIME2 demux QZV.

    The file 'per-sample-fastq-counts.tsv' is written by `qiime demux
    summarize` and contains one row per sample with forward and reverse
    read counts.  We use the forward count because that is the one that
    drives the ≥10 k threshold check — paired reads are always equal or
    lower after quality filtering.

    The file is TAB-separated (not comma).  QIIME2 has used slightly
    different column names across versions; we try a ranked list so the
    parser stays robust to minor schema changes.

    Args:
        qzv_path: Path to a demux.qzv or demux_trimmed.qzv file.

    Returns:
        Dict mapping sample-id strings to integer forward read counts.

    Raises:
        FileNotFoundError: If the TSV is not present inside the QZV.
        ValueError:        If no counts can be parsed from the TSV.
    """
    counts: Dict[str, int] = {}

    with zipfile.ZipFile(qzv_path) as zf:
        entry = _find_in_zip(zf, "per-sample-fastq-counts.tsv")
        if entry is None:
            raise FileNotFoundError(
                f"'per-sample-fastq-counts.tsv' not found inside {qzv_path}.\n"
                f"Make sure this is a demux summarize QZV, not a raw demux QZA."
            )

        # Read the file and decode as UTF-8; QZV files are always UTF-8.
        content = zf.read(entry).decode("utf-8")

    lines = [l for l in content.splitlines() if l.strip()]
    if not lines:
        raise ValueError(f"Empty per-sample-fastq-counts.tsv in {qzv_path}")

    # Parse the tab-separated header to find the column indices we need.
    # QIIME2 has used both "forward sequence count" and "Sequence Count"
    # (and variants) across releases, so we try multiple names in order.
    header = lines[0].split("\t")

    # Ranked list of known column name variants for sample ID and fwd count.
    id_candidates  = ["sample-id", "Sample ID", "SampleID"]
    cnt_candidates = [
        "forward sequence count",
        "Forward Sequence Count",
        "sequence count",
        "Sequence Count",
    ]

    id_col  = next((h for h in id_candidates  if h in header), None)
    cnt_col = next((h for h in cnt_candidates if h in header), None)

    if id_col is None or cnt_col is None:
        raise ValueError(
            f"Could not locate required columns in per-sample-fastq-counts.tsv.\n"
            f"Found headers: {header}\n"
            f"Expected one of {id_candidates} and one of {cnt_candidates}."
        )

    id_idx  = header.index(id_col)
    cnt_idx = header.index(cnt_col)

    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) <= max(id_idx, cnt_idx):
            continue  # skip malformed rows rather than crashing
        sid = parts[id_idx].strip()
        try:
            # Counts may be formatted with commas (e.g. "12,345") in some versions.
            counts[sid] = int(parts[cnt_idx].replace(",", "").strip())
        except ValueError:
            log.debug("Could not parse count for sample '%s': %r", sid, parts[cnt_idx])

    if not counts:
        raise ValueError(f"No sample counts parsed from {qzv_path}")

    return counts


# ---------------------------------------------------------------------------
# Cutadapt log parsing
# ---------------------------------------------------------------------------
# cutadapt's --verbose output (which the pipeline already passes) writes a
# plain-text summary at the end of each run.  We extract four numbers from it.
#
# Paired-end summary format (from qiime cutadapt trim-paired):
#   Total read pairs processed:      1,234,567
#     Read 1 with adapter:             987,654 (80.0%)
#     Read 2 with adapter:             876,543 (71.0%)
#   Pairs that were too short:          12,345 (1.0%)
#   Pairs written (passing filters):   987,654 (80.0%)
#
# Single-end summary format:
#   Total reads processed:           1,234,567
#   Reads with adapters:               987,654 (80.0%)
#   Reads that were too short:          12,345 (1.0%)
#   Reads written (passing filters):   987,654 (80.0%)
#
# Additionally, when --verbose is set, cutadapt prints one
# "Overview of removed sequences" block per adapter sequence.  Each block
# is a histogram: rows are (position, count, ...) where 'position' is
# where in the read the adapter was found.  The modal position equals the
# effective insert length — this is the most direct way to confirm that
# your amplicon is the length you think it is.
# ---------------------------------------------------------------------------

# Pre-compiled regex patterns for the cutadapt summary block.
# These cover both paired-end and single-end output formats.
_RE_PAIRS_TOTAL   = re.compile(r"Total read pairs processed:\s+([\d,]+)")
_RE_SINGLE_TOTAL  = re.compile(r"Total reads processed:\s+([\d,]+)")
_RE_R1_ADAPTER    = re.compile(r"Read 1 with adapter:\s+([\d,]+)")
_RE_R2_ADAPTER    = re.compile(r"Read 2 with adapter:\s+([\d,]+)")
_RE_SE_ADAPTER    = re.compile(r"Reads with adapters:\s+([\d,]+)")
_RE_PAIRS_SHORT   = re.compile(r"Pairs that were too short:\s+([\d,]+)")
_RE_SINGLE_SHORT  = re.compile(r"Reads that were too short:\s+([\d,]+)")
_RE_PAIRS_WRITTEN = re.compile(r"Pairs written \(passing filters\):\s+([\d,]+)")
_RE_SINGLE_WRITTEN= re.compile(r"Reads written \(passing filters\):\s+([\d,]+)")

# The histogram block starts with this header line (case-insensitive).
_RE_HIST_HEADER   = re.compile(r"^length\s+count\s+expect", re.IGNORECASE)
# Each data row: "   5   12345   ..."  — we only need the first two columns.
_RE_HIST_ROW      = re.compile(r"^\s*(\d+)\s+(\d+)")


def _parse_int(match: Optional[re.Match]) -> int:
    """Extract a comma-formatted integer from a regex match, or return 0."""
    if match is None:
        return 0
    return int(match.group(1).replace(",", ""))


def parse_cutadapt_log(log_path: Path) -> dict:
    """
    Parse a cutadapt verbose log and return a structured summary dict.

    The pipeline's cmd_cutadapt already passes --verbose, so this histogram
    data is always present in the log files written to logs/{marker}/{dataset}/.

    Handles both paired-end and single-end log formats automatically by
    trying the paired-end regex first (more informative) and falling back
    to single-end patterns.

    Returns a dict with keys:
        total          (int)   — total read pairs/reads processed
        r1_adapter     (int)   — R1 reads where forward primer was found
        r2_adapter     (int)   — R2 reads where reverse primer was found (0 if SE)
        too_short      (int)   — pairs/reads discarded as too short after trim
        written        (int)   — pairs/reads that passed all filters
        detect_rate_r1 (float) — r1_adapter / total as a fraction 0–1
        detect_rate_r2 (float) — r2_adapter / total (0.0 if single-end)
        dimer_rate     (float) — too_short / total as a fraction 0–1
        histograms     (list)  — one dict per adapter block, each with:
                                   {positions: [int, ...], counts: [int, ...],
                                    modal_pos: int, label: str}
    """
    text = log_path.read_text(encoding="utf-8", errors="replace")

    # Try paired-end patterns first; fall back to single-end if no match.
    paired = _RE_PAIRS_TOTAL.search(text)
    if paired:
        total     = _parse_int(paired)
        r1_adapt  = _parse_int(_RE_R1_ADAPTER.search(text))
        r2_adapt  = _parse_int(_RE_R2_ADAPTER.search(text))
        too_short = _parse_int(_RE_PAIRS_SHORT.search(text))
        written   = _parse_int(_RE_PAIRS_WRITTEN.search(text))
    else:
        total     = _parse_int(_RE_SINGLE_TOTAL.search(text))
        r1_adapt  = _parse_int(_RE_SE_ADAPTER.search(text))
        r2_adapt  = 0   # single-end has no R2
        too_short = _parse_int(_RE_SINGLE_SHORT.search(text))
        written   = _parse_int(_RE_SINGLE_WRITTEN.search(text))

    # Avoid division-by-zero if the log is empty or parsing failed.
    denom = total if total > 0 else 1

    # ── Adapter position histograms ───────────────────────────────────────
    # cutadapt prints one "Overview of removed sequences" block per adapter
    # sequence, in the order they were passed on the command line:
    #   block 0 → forward primer (--front-f)
    #   block 1 → reverse primer (--front-r)
    #   block 2 → RC of reverse primer (--adapter-f, read-through trimming)
    #   block 3 → RC of forward primer (--adapter-r, read-through trimming)
    # For a short amplicon like a 110 bp virus fragment, the read-through
    # adapters (blocks 2/3) are the most important: their modal position
    # directly tells you the insert length.
    histograms = []
    current_positions: List[int] = []
    current_counts:    List[int] = []
    in_histogram = False
    block_index  = 0

    for line in text.splitlines():
        # A new adapter section always starts with "=== Adapter" or "=== Read"
        if line.startswith("=== "):
            # Save the previous histogram if it had any data.
            if current_positions:
                modal = current_positions[
                    current_counts.index(max(current_counts))
                ]
                histograms.append({
                    "label":      f"Adapter block {block_index}",
                    "positions":  current_positions,
                    "counts":     current_counts,
                    "modal_pos":  modal,
                })
                block_index += 1
            current_positions = []
            current_counts    = []
            in_histogram      = False
            continue

        # The histogram block starts after its own header line.
        if _RE_HIST_HEADER.match(line):
            in_histogram = True
            continue

        if in_histogram:
            m = _RE_HIST_ROW.match(line)
            if m:
                pos   = int(m.group(1))
                count = int(m.group(2))
                if count > 0:   # skip zero-count rows to keep histograms compact
                    current_positions.append(pos)
                    current_counts.append(count)
            elif line.strip() == "":
                # A blank line ends the histogram block.
                in_histogram = False

    # Flush the last histogram block if the file ended without a blank line.
    if current_positions:
        modal = current_positions[current_counts.index(max(current_counts))]
        histograms.append({
            "label":     f"Adapter block {block_index}",
            "positions": current_positions,
            "counts":    current_counts,
            "modal_pos": modal,
        })

    return {
        "total":          total,
        "r1_adapter":     r1_adapt,
        "r2_adapter":     r2_adapt,
        "too_short":      too_short,
        "written":        written,
        "detect_rate_r1": r1_adapt  / denom,
        "detect_rate_r2": r2_adapt  / denom,
        "dimer_rate":     too_short / denom,
        "histograms":     histograms,
    }


# ---------------------------------------------------------------------------
# Check report printer
# ---------------------------------------------------------------------------

# Display symbols — short constants keep the report code readable.
_PASS = "✅ PASS"
_WARN = "⚠️  WARN"
_FAIL = "❌ FAIL"
_INFO = "ℹ️  INFO"


def _flag(value: float, pass_threshold: float, warn_threshold: float,
          higher_is_better: bool = True) -> str:
    """
    Return a PASS / WARN / FAIL symbol based on a numeric value and two
    thresholds.

    higher_is_better=True  (default): used for detection rate — higher is good.
    higher_is_better=False:           used for dimer rate — lower is good.
    """
    if higher_is_better:
        if value >= pass_threshold:
            return _PASS
        if value >= warn_threshold:
            return _WARN
        return _FAIL
    else:
        if value <= pass_threshold:
            return _PASS
        if value <= warn_threshold:
            return _WARN
        return _FAIL


def _bar(count: int, max_count: int, width: int = 30) -> str:
    """
    Build a proportional ASCII bar for histogram rows.
    Using block characters (█) rather than dashes so peaks are visually obvious.
    """
    if max_count == 0:
        return ""
    filled = int(width * count / max_count)
    return "█" * filled


def print_check_report(
    marker:         str,
    raw_counts:     Dict[str, int],
    trimmed_counts: Dict[str, int],
    cutadapt:       Optional[dict],
    amplicon_len:   Optional[int],
    primer_f_len:   int,
    min_reads:      int,
) -> List[str]:
    """
    Print a formatted post-cutadapt QC report and return a list of issue keys.

    The report has four sections:
        1. Per-sample read counts (raw vs trimmed, pass/warn/fail per sample)
        2. Adapter detection rate (did cutadapt find the primers?)
        3. Dimer rate (how many reads were discarded as too short?)
        4. Adapter position histograms (where did adapters appear in the read?)

    Returns a list of string issue keys so that cmd_check can exit non-zero
    when a FAIL-level problem is detected.  Current issue keys:
        "samples_below_threshold" — one or more samples < min_reads
        "low_detection_rate"      — R1 adapter found in < 50% of reads
        "high_dimer_rate"         — > 15% of reads discarded as too short
        "amplicon_length_mismatch"— modal adapter position differs > 15 bp
                                    from the expected amplicon length

    Args:
        marker:         Marker name string (for display only).
        raw_counts:     {sample_id: count} from pre-trim demux QZV.
        trimmed_counts: {sample_id: count} from post-trim demux QZV.
        cutadapt:       Parsed cutadapt log dict from parse_cutadapt_log(),
                        or None if no log was provided.
        amplicon_len:   Expected amplicon length in bp, or None if unknown.
        primer_f_len:   Forward primer length in bp (used to flag dimer peaks).
        min_reads:      Minimum acceptable read count per sample.
    """
    issues: List[str] = []
    SEP  = "=" * 70
    DASH = "─" * 70

    print(SEP)
    print(f"  POST-CUTADAPT QC REPORT   marker = {marker}")
    print(SEP)

    # ── Section 1: Per-sample read counts ────────────────────────────────
    # We compare raw vs trimmed counts for every sample.  A sample that
    # drops below min_reads after trimming will produce sparse or empty
    # DADA2 output — better to know now.
    print(f"\n{DASH}")
    print(f"  1. Per-Sample Read Counts   (threshold: {min_reads:,} reads)")
    print(DASH)

    all_ids  = sorted(set(list(raw_counts) + list(trimmed_counts)))
    n_pass = n_warn = n_fail = 0
    low_samples: List[Tuple[str, int]] = []

    # Print a fixed-width table with one row per sample.
    print(f"  {'Sample':<42} {'Raw':>10}  {'Trimmed':>10}  {'Kept':>7}  Status")
    print(f"  {'-'*42} {'-'*10}  {'-'*10}  {'-'*7}  ------")

    for sid in all_ids:
        raw  = raw_counts.get(sid)
        trim = trimmed_counts.get(sid)
        # Use trimmed count for the threshold check if available; otherwise raw.
        check = trim if trim is not None else raw

        # Calculate retention percentage for display.
        kept = f"{100 * trim / raw:.0f}%" if (raw and trim) else "N/A"

        if check is None:
            status = _INFO
        elif check >= min_reads:
            status = _PASS
            n_pass += 1
        elif check >= min_reads * 0.5:
            # 50–100% of threshold: warn rather than hard-fail because the
            # sample may still produce useful ASVs at lower depth.
            status = _WARN
            n_warn += 1
            low_samples.append((sid, check))
        else:
            status = _FAIL
            n_fail += 1
            low_samples.append((sid, check))

        raw_str  = f"{raw:>10,}"  if raw  is not None else f"{'N/A':>10}"
        trim_str = f"{trim:>10,}" if trim is not None else f"{'N/A':>10}"
        print(f"  {sid:<42} {raw_str}  {trim_str}  {kept:>7}  {status}")

    print(f"\n  Summary: {n_pass} PASS  {n_warn} WARN  {n_fail} FAIL")
    if low_samples:
        issues.append("samples_below_threshold")
        print(f"\n  Low-read samples:")
        for sid, cnt in low_samples:
            print(f"    {sid}: {cnt:,} reads (threshold: {min_reads:,})")

    # ── Section 2: Adapter detection rate ────────────────────────────────
    # If cutadapt found adapters in < 80% of reads, something is wrong:
    # either the primer sequences passed to cutadapt are incorrect, or the
    # reads are genuinely shorter than the amplicon (no bleed-through) and
    # trimming was unnecessary.  In the latter case the trimmed QZV read
    # counts will match the raw counts exactly.
    print(f"\n{DASH}")
    print(f"  2. Adapter Detection Rate")
    print(DASH)

    if cutadapt is None:
        print(f"  {_INFO} No cutadapt log provided.")
        print(f"       Add --cutadapt-log logs/{marker}/*/cutadapt_*.log to enable this check.")
    else:
        total = cutadapt["total"]
        det   = cutadapt["detect_rate_r1"]
        det_r2= cutadapt["detect_rate_r2"]

        print(f"  Total read pairs processed : {total:>12,}")
        print(f"  R1 with adapter detected   : {cutadapt['r1_adapter']:>12,}  "
              f"({det:.1%})  {_flag(det, 0.8, 0.5)}")
        if cutadapt["r2_adapter"] > 0:
            print(f"  R2 with adapter detected   : {cutadapt['r2_adapter']:>12,}  "
                  f"({det_r2:.1%})  {_flag(det_r2, 0.8, 0.5)}")
        print(f"  Pairs written              : {cutadapt['written']:>12,}  "
              f"({cutadapt['written'] / (total or 1):.1%})")

        if det < 0.5:
            issues.append("low_detection_rate")
            print(f"\n  {_FAIL} R1 adapter detection < 50%.")
            print(f"       Possible causes:")
            print(f"         (a) Wrong primer sequences passed to cutadapt.")
            print(f"         (b) Reads shorter than amplicon — no adapter bleed-through.")
            print(f"         (c) Very high sequencing error rate in the primer region.")
        elif det < 0.8:
            print(f"\n  {_WARN} Detection 50–80%. Check primer sequences match your library prep.")
        else:
            print(f"\n  {_PASS} Good adapter detection (>80%).")

    # ── Section 3: Dimer rate ─────────────────────────────────────────────
    # Reads flagged "too short" by cutadapt are primer dimers or adapter
    # dimers — the primers ligated to each other with no real template
    # in between, producing extremely short inserts.  They typically appear
    # as a peak at very small positions (< 2× primer length) in the adapter
    # histogram.  A high dimer rate indicates a library prep problem.
    print(f"\n{DASH}")
    print(f"  3. Adapter / Primer Dimer Rate")
    print(DASH)

    if cutadapt is None:
        print(f"  {_INFO} No cutadapt log provided — skipping dimer check.")
    else:
        dr = cutadapt["dimer_rate"]
        # Thresholds: < 5% is fine, 5–15% warrants monitoring, > 15% is a problem.
        status = _flag(dr, 0.05, 0.15, higher_is_better=False)
        print(f"  Pairs too short after trim : {cutadapt['too_short']:>12,}  "
              f"({dr:.1%})  {status}")

        if dr > 0.15:
            issues.append("high_dimer_rate")
            dimer_len = primer_f_len * 2 if primer_f_len else "?"
            print(f"\n  {_FAIL} > 15% reads discarded as too short.")
            print(f"       These are likely primer–primer dimers (~{dimer_len} bp).")
            print(f"       Options: gel-purify the library, reduce primer concentrations,")
            print(f"                or add --p-minimum-length to cutadapt to set a hard floor.")
        elif dr > 0.05:
            print(f"\n  {_WARN} 5–15% short reads. Moderate dimer contamination — monitor.")
        else:
            print(f"\n  {_PASS} Dimer rate < 5% — acceptable.")

    # ── Section 4: Adapter position histograms ───────────────────────────
    # Each adapter block from cutadapt's verbose output is printed as a
    # compact ASCII bar chart.  The position axis is "where in the read was
    # the adapter found," so the modal position = effective insert length.
    # For your 110 bp virus fragment you expect the peak near position 110.
    # A peak at very low positions confirms dimer contamination.
    print(f"\n{DASH}")
    print(f"  4. Adapter Position Histograms  (modal position = insert length)")
    print(DASH)

    if cutadapt is None:
        print(f"  {_INFO} No cutadapt log provided — skipping histogram check.")
    elif not cutadapt["histograms"]:
        print(f"  {_INFO} No histogram data found in the log.")
        print(f"       The pipeline's cmd_cutadapt already passes --verbose,")
        print(f"       so check that the log file path is correct.")
    else:
        for hist in cutadapt["histograms"]:
            modal  = hist["modal_pos"]
            label  = hist["label"]
            positions = hist["positions"]
            counts    = hist["counts"]
            max_count = max(counts) if counts else 1

            print(f"\n  {label}   modal position = {modal} bp")

            # Flag a mismatch between the modal position and the expected
            # amplicon length — a discrepancy > 15 bp suggests something
            # unexpected happened (wrong amplicon, contamination, etc.).
            if amplicon_len is not None:
                diff = abs(modal - amplicon_len)
                if diff <= 5:
                    print(f"  {_PASS} Modal position matches expected amplicon "
                          f"({amplicon_len} bp, diff = {diff} bp)")
                elif diff <= 15:
                    print(f"  {_WARN} Modal {modal} bp vs expected {amplicon_len} bp "
                          f"(diff = {diff} bp) — within tolerance but worth checking.")
                else:
                    issues.append("amplicon_length_mismatch")
                    print(f"  {_FAIL} Modal {modal} bp differs from expected {amplicon_len} bp "
                          f"by {diff} bp — verify primer sequences and amplicon size.")

            # Print the histogram rows, limiting to top positions by count
            # so the table stays compact even for very diverse distributions.
            top_n  = 15
            sorted_by_pos = sorted(zip(positions, counts), key=lambda x: x[0])
            # Show up to top_n positions; if the histogram is dense we show
            # the ones with the highest counts (most informative).
            if len(sorted_by_pos) > top_n:
                top_by_count = sorted(sorted_by_pos, key=lambda x: -x[1])[:top_n]
                display_rows = sorted(top_by_count, key=lambda x: x[0])
                print(f"  (showing {top_n} of {len(sorted_by_pos)} positions by count)")
            else:
                display_rows = sorted_by_pos

            print(f"  {'Position':>10}  {'Count':>10}  Distribution")
            for pos, cnt in display_rows:
                bar  = _bar(cnt, max_count)
                note = ""
                # Flag positions that are suspiciously short — likely dimers.
                if primer_f_len and pos < primer_f_len * 2:
                    note = "  ← probable primer/adapter dimer"
                print(f"  {pos:>10}  {cnt:>10,}  {bar}{note}")

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\n{SEP}")
    if issues:
        print(f"  {len(issues)} issue(s) detected:")
        for iss in issues:
            print(f"    • {iss.replace('_', ' ')}")
    else:
        print(f"  {_PASS} All checks passed.  Proceed to DADA2.")
    print(SEP)
    print()

    return issues


# ===========================================================================
# SECTION 6 — CLI
# ===========================================================================

def _add_common_detect_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--reads-dir", default="reads",
                   help="Root reads directory with per-marker subdirectories. Default: reads/")
    p.add_argument("--n-reads", type=int, default=200,
                   help="Reads to sample per marker. Default: 200")
    p.add_argument("--mismatches", type=int, default=2,
                   help="Allowed mismatches when matching primers. Default: 2")
    p.add_argument("--verbose", action="store_true",
                   help="Show all primer scores, not just the best match.")


def _add_common_suggest_args(p: argparse.ArgumentParser) -> None:
    # Note: --demux is NOT added here so that suggest and run can set
    # different required= values without touching argparse internals.
    p.add_argument("--amplicon-length", type=int, default=None,
                   help="Expected amplicon length in bp (excluding primers). "
                        "Required for overlap check. Examples: cytb=324, 16S V4=253, MiFish=180.")
    p.add_argument("--min-quality", type=float, default=DEFAULT_MIN_QUAL,
                   help=f"Minimum acceptable median Phred Q score. Default: {DEFAULT_MIN_QUAL}")
    p.add_argument("--min-overlap", type=int, default=DEFAULT_MIN_OVER,
                   help=f"Minimum bp overlap for paired-end merging. Default: {DEFAULT_MIN_OVER}")
    p.add_argument("--single-end", action="store_true",
                   help="Force single-end mode recommendation.")
    p.add_argument("--out-json", type=Path, default=None,
                   help="Optional path to write DADA2 params as JSON.")


def build_parser() -> argparse.ArgumentParser:
    """
    Build and return the top-level argument parser with all four subcommands.

    Subcommand structure:
      detect  --marker / --all  --reads-dir  [--n-reads] [--mismatches] [--report] [--verbose]
      suggest --demux --primer-f --primer-r  [--marker] [--amplicon-length] ...
      run     --marker --reads-dir           [--demux]  [--amplicon-length] ...
      check   --marker --trimmed-qzv         [--raw-qzv] [--cutadapt-log] [--amplicon-len] ...
    """
    parser = argparse.ArgumentParser(
        prog="primer_advisor.py",
        description=(
            "Primer detection and DADA2 parameter suggestion for metabarcoding.\n\n"
            "Subcommands:\n"
            "  detect  — identify primers from raw reads\n"
            "  suggest — recommend DADA2 params from a demux QZV\n"
            "  run     — detect primers then suggest DADA2 params in one step\n"
            "  check   — verify cutadapt results: read counts, detection rate, dimers"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="subcommand", required=True)

    # ── detect ────────────────────────────────────────────────────────────────
    det = sub.add_parser("detect",
                         help="Identify primers from raw R1/R2 reads.")
    det_group = det.add_mutually_exclusive_group(required=True)
    det_group.add_argument("--marker", default=None,
                           help="Single marker to check (e.g. 18S, ITS1-2, MiFish).")
    det_group.add_argument("--all", action="store_true",
                           help="Check all markers found in --reads-dir.")
    _add_common_detect_args(det)
    det.add_argument("--report", default=None,
                     help="Optional path to write a TSV summary report.")

    # ── suggest ───────────────────────────────────────────────────────────────
    sug = sub.add_parser("suggest",
                         help="Suggest DADA2 truncation params from a demux QZV.")
    sug.add_argument("--demux", required=True, type=Path,
                     help="Path to QIIME2 demux.qzv (preferred) or demux.qza.")
    sug.add_argument("--marker", default="unknown",
                     help="Marker name for display and output naming. Default: unknown")
    sug.add_argument("--primer-f", required=True, type=int,
                     help="Forward primer length in bp (used as --trim-left-f in DADA2).")
    sug.add_argument("--primer-r", required=True, type=int,
                     help="Reverse primer length in bp (used as --trim-left-r in DADA2).")
    _add_common_suggest_args(sug)
    sug.add_argument("--verbose", action="store_true", help="Enable debug logging.")

    # ── run ───────────────────────────────────────────────────────────────────
    run = sub.add_parser("run",
                         help="Detect primers then suggest DADA2 params automatically.")
    run.add_argument("--marker", required=True,
                     help="Marker to process (e.g. MiFish, 16S, cytb).")
    run.add_argument("--demux", required=False, type=Path, default=None,
                     help="Path to QIIME2 demux.qzv. If omitted, only primer detection is run.")
    _add_common_detect_args(run)
    _add_common_suggest_args(run)

    # ── check ─────────────────────────────────────────────────────────────────
    # This subcommand runs AFTER cutadapt and answers: did trimming actually
    # work?  It requires at least the trimmed QZV; the raw QZV and cutadapt
    # log are optional but unlock the most useful checks.
    chk = sub.add_parser(
        "check",
        help="Verify cutadapt results: read counts, adapter detection rate, dimer rate.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Post-cutadapt QC check.\n\n"
            "Run AFTER 'qiime cutadapt trim-paired' and 'qiime demux summarize'\n"
            "on the trimmed artifact.  Checks:\n"
            "  1. Per-sample read counts vs --min-reads threshold\n"
            "  2. Adapter detection rate (were primers found?)\n"
            "  3. Dimer rate (reads discarded as too short)\n"
            "  4. Adapter position histograms (confirms amplicon length)\n"
        ),
    )
    chk.add_argument("--marker", required=True,
                     help="Marker name (for display and log path pattern matching).")
    chk.add_argument("--trimmed-qzv", required=True, type=Path,
                     help="demux_trimmed.qzv produced by qiime demux summarize on the "
                          "cutadapt output.  Required.")
    chk.add_argument("--raw-qzv", default=None, type=Path,
                     help="demux.qzv BEFORE cutadapt.  Optional but enables before/after "
                          "read count comparison.")
    chk.add_argument("--cutadapt-log", default=None,
                     help="Path or glob to the cutadapt log written by the pipeline, e.g. "
                          "logs/MiFish/all/cutadapt_*.log.  Required for detection rate, "
                          "dimer rate, and histogram checks.  Quote globs to prevent shell "
                          "expansion before Python sees them.")
    chk.add_argument("--amplicon-len", type=int, default=None,
                     help="Expected amplicon length in bp (e.g. 110 for the virus fragment, "
                          "180 for MiFish).  Used to validate the modal adapter position.")
    chk.add_argument("--primer-f-len", type=int, default=0,
                     help="Forward primer length in bp.  Used to flag dimer peaks "
                          "(peaks at < 2× primer length are likely dimers).  Default: 0")
    chk.add_argument("--min-reads", type=int, default=10000,
                     help="Minimum acceptable read count per sample after trimming. "
                          "Default: 10000")
    chk.add_argument("--out-tsv", type=Path, default=None,
                     help="Optional path to write per-sample read counts as TSV.")

    return parser


# ===========================================================================
# SECTION 7 — Subcommand handlers
# ===========================================================================

def cmd_detect(args: argparse.Namespace) -> int:
    """
Handler for the 'detect' subcommand.

    Samples reads for one or all markers, scores every primer in PRIMER_DB,
    and prints the best-matching forward/reverse pair with a recommended
    cutadapt command. Optionally writes a TSV summary report.

    Returns 0 on success, 1 if the reads directory is not found.
    """
    reads_dir = Path(args.reads_dir)
    if not reads_dir.exists():
        print(f"[ERROR] Reads directory not found: {reads_dir}", file=sys.stderr)
        return 1

    markers = discover_markers(reads_dir) if args.all else [args.marker]
    if args.all:
        print(f"Found markers: {markers}")

    results = []
    for marker in markers:
        result = detect_primers_for_marker(
            marker=marker, reads_dir=reads_dir,
            n_reads=args.n_reads, mismatch=args.mismatches,
            verbose=args.verbose,
        )
        results.append(result)

    if args.report:
        write_detect_report(results, Path(args.report))

    return 0


def cmd_suggest(args: argparse.Namespace) -> int:
    """
Handler for the 'suggest' subcommand.

    Opens a QIIME2 demux QZV, extracts per-position quality profiles, finds
    the quality dropoff positions, checks paired-end overlap against the
    expected amplicon length, and prints recommended DADA2 parameters.
    Optionally writes a JSON params file for pipeline use.

    Returns 0 on success, 1 on any error parsing the demux file.
    """
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    demux_path = Path(args.demux).resolve()
    log.info("Marker: %s | demux: %s", args.marker, demux_path)
    log.info("primer_f=%d  primer_r=%d  amplicon=%s  min_q=%.0f  min_overlap=%d",
             args.primer_f, args.primer_r, args.amplicon_length,
             args.min_quality, args.min_overlap)

    try:
        result = suggest_params(
            demux_path=demux_path,
            primer_f=args.primer_f, primer_r=args.primer_r,
            amplicon_length=args.amplicon_length,
            min_quality=args.min_quality, min_overlap=args.min_overlap,
            force_single_end=args.single_end,
        )
    except Exception as e:
        log.error("Failed to analyse demux file: %s", e)
        if args.verbose:
            import traceback; traceback.print_exc()
        return 1

    print_suggestion(result, marker=args.marker)

    if args.out_json:
        write_suggest_json(result, marker=args.marker, out_path=Path(args.out_json))

    if result["warnings"]:
        log.warning("%d warning(s). Review before running DADA2.", len(result["warnings"]))

    return 0


def cmd_run(args: argparse.Namespace) -> int:
    """Detect primers, then use detected lengths to suggest DADA2 params."""
    reads_dir = Path(args.reads_dir)
    if not reads_dir.exists():
        print(f"[ERROR] Reads directory not found: {reads_dir}", file=sys.stderr)
        return 1

    # Step 1 — detect
    detection = detect_primers_for_marker(
        marker=args.marker, reads_dir=reads_dir,
        n_reads=args.n_reads, mismatch=args.mismatches,
        verbose=args.verbose,
    )

    if "error" in detection:
        return 1

    primer_f_len = detection["fwd_len"]
    primer_r_len = detection["rev_len"]
    print(f"\n  Detected primer lengths: fwd={primer_f_len} bp, rev={primer_r_len} bp")
    print(f"  These will be used as --primer-f / --primer-r for DADA2 suggestion.\n")

    # Step 2 — suggest (only if --demux was provided)
    if not args.demux:
        print("No --demux provided; skipping DADA2 parameter suggestion.")
        print("Re-run with --demux path/to/demux.qzv to get truncation recommendations.")
        return 0

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    demux_path = Path(args.demux).resolve()
    log.info("Proceeding to DADA2 suggestion with detected primer lengths...")

    try:
        result = suggest_params(
            demux_path=demux_path,
            primer_f=primer_f_len, primer_r=primer_r_len,
            amplicon_length=args.amplicon_length,
            min_quality=args.min_quality, min_overlap=args.min_overlap,
            force_single_end=args.single_end,
        )
    except Exception as e:
        log.error("Failed to analyse demux file: %s", e)
        if args.verbose:
            import traceback; traceback.print_exc()
        return 1

    print_suggestion(result, marker=args.marker)

    if args.out_json:
        write_suggest_json(result, marker=args.marker, out_path=Path(args.out_json))

    if result["warnings"]:
        log.warning("%d warning(s). Review before running DADA2.", len(result["warnings"]))

    return 0


def cmd_check(args: argparse.Namespace) -> int:
    """
    Handler for the 'check' subcommand.

    Loads per-sample read counts from the post-trim QZV (and optionally the
    pre-trim QZV for comparison), parses the cutadapt verbose log if provided,
    then delegates all reporting to print_check_report().

    Exits with code 1 if any FAIL-level issue is detected so the pipeline
    can catch failures in a shell script with 'set -e' or an explicit
    exit-code check.  WARN-level issues print but still exit 0.

    Returns 0 (success / warn only) or 1 (one or more FAIL-level issues).
    """
    # ── Load trimmed read counts (required) ───────────────────────────────
    trimmed_path = Path(args.trimmed_qzv)
    if not trimmed_path.exists():
        print(f"[ERROR] --trimmed-qzv not found: {trimmed_path}", file=sys.stderr)
        return 1

    try:
        trimmed_counts = parse_per_sample_counts(trimmed_path)
    except Exception as e:
        print(f"[ERROR] Could not read trimmed QZV: {e}", file=sys.stderr)
        return 1

    # ── Load raw read counts (optional — enables before/after comparison) ─
    raw_counts: Dict[str, int] = {}
    if args.raw_qzv is not None:
        raw_path = Path(args.raw_qzv)
        if not raw_path.exists():
            log.warning("--raw-qzv not found (%s); skipping before/after comparison.", raw_path)
        else:
            try:
                raw_counts = parse_per_sample_counts(raw_path)
            except Exception as e:
                log.warning("Could not read raw QZV (%s); skipping: %s", raw_path, e)

    # ── Load cutadapt log (optional — enables detection/dimer/histogram) ──
    # The pipeline writes one log per cutadapt run with a timestamp suffix,
    # so we support glob patterns here (e.g. "logs/MiFish/all/cutadapt_*.log").
    # We always use the most recent match to avoid double-counting.
    cutadapt: Optional[dict] = None
    if args.cutadapt_log is not None:
        matches = sorted(glob.glob(args.cutadapt_log))
        if not matches:
            log.warning(
                "No files matched --cutadapt-log pattern '%s'.  "
                "Quote the pattern to prevent shell expansion, e.g.: "
                "--cutadapt-log 'logs/MiFish/all/cutadapt_*.log'",
                args.cutadapt_log,
            )
        else:
            log_path = Path(matches[-1])   # most recent log file
            log.info("Reading cutadapt log: %s", log_path)
            try:
                cutadapt = parse_cutadapt_log(log_path)
            except Exception as e:
                log.warning("Could not parse cutadapt log (%s): %s", log_path, e)

    # ── Run the report ─────────────────────────────────────────────────────
    issues = print_check_report(
        marker        = args.marker,
        raw_counts    = raw_counts,
        trimmed_counts= trimmed_counts,
        cutadapt      = cutadapt,
        amplicon_len  = args.amplicon_len,
        primer_f_len  = args.primer_f_len,
        min_reads     = args.min_reads,
    )

    # ── Optional TSV output ────────────────────────────────────────────────
    # Write a simple per-sample TSV useful for pasting into a lab notebook
    # or loading into R/Excel for a quick sanity check figure.
    if args.out_tsv:
        out_tsv = Path(args.out_tsv)
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        all_ids = sorted(set(list(raw_counts) + list(trimmed_counts)))
        with out_tsv.open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["sample_id", "raw_reads", "trimmed_reads",
                        "retention_pct", "passes_threshold"])
            for sid in all_ids:
                raw  = raw_counts.get(sid, "")
                trim = trimmed_counts.get(sid, "")
                kept = f"{100 * trim / raw:.1f}" if (raw and trim) else ""
                passes = "yes" if (trim or 0) >= args.min_reads else "no"
                w.writerow([sid, raw, trim, kept, passes])
        log.info("Per-sample TSV written: %s", out_tsv)

    # FAIL-level issues cause a non-zero exit so shell pipelines can catch them.
    fail_issues = {"low_detection_rate", "high_dimer_rate", "amplicon_length_mismatch"}
    if any(iss in fail_issues for iss in issues):
        return 1
    return 0


# ===========================================================================
# SECTION 8 — Entry point
# ===========================================================================

def main(argv: Optional[List[str]] = None) -> int:
    """
    Entry point. Parses argv (or sys.argv if None) and dispatches to the
    appropriate subcommand handler.

    Subcommand → handler mapping:
        detect  → cmd_detect   (primer detection from raw reads)
        suggest → cmd_suggest  (DADA2 trunc params from demux QZV)
        run     → cmd_run      (detect + suggest chained)
        check   → cmd_check    (post-cutadapt QC: counts, detection, dimers)

    Returns the integer exit code (0 = success, 1 = error or FAIL-level issue).
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.subcommand == "detect":
        return cmd_detect(args)
    elif args.subcommand == "suggest":
        return cmd_suggest(args)
    elif args.subcommand == "run":
        return cmd_run(args)
    elif args.subcommand == "check":
        return cmd_check(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
