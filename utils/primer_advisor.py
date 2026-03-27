#!/usr/bin/env python3
"""
primer_advisor.py
=================
Combined primer detection and DADA2 parameter suggestion tool.
Replaces detect_primers.py and suggest_dada2_params.py.

Subcommands
-----------
detect  — Identify primers from raw R1/R2 reads and recommend a cutadapt command.
suggest — Suggest DADA2 truncation parameters from a QIIME2 demux QZV.
run     — Chain both: detect primers automatically, then suggest DADA2 params.

Usage
-----
# Detect primers for one marker:
python primer_advisor.py detect --marker MiFish --reads-dir reads/

# Detect all markers and write a TSV report:
python primer_advisor.py detect --all --reads-dir reads/ --report primers_detected.tsv

# Suggest DADA2 params (primer lengths already known):
python primer_advisor.py suggest \
    --demux qiime2/MiFish/all/imported/demux.qzv \
    --primer-f 21 --primer-r 27 \
    --amplicon-length 180 --marker MiFish

# Full automatic run — detect primers then suggest params:
python primer_advisor.py run \
    --marker MiFish \
    --reads-dir reads/ \
    --demux qiime2/MiFish/all/imported/demux.qzv \
    --amplicon-length 180

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
# SECTION 5 — CLI
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
Build and return the top-level argument parser with all three subcommands.

    Subcommand structure:
      detect  --marker / --all  --reads-dir  [--n-reads] [--mismatches] [--report] [--verbose]
      suggest --demux --primer-f --primer-r  [--marker] [--amplicon-length] ...
      run     --marker --reads-dir           [--demux]  [--amplicon-length] ...
    """
    parser = argparse.ArgumentParser(
        prog="primer_advisor.py",
        description=(
            "Primer detection and DADA2 parameter suggestion for metabarcoding.\n\n"
            "Subcommands:\n"
            "  detect  — identify primers from raw reads\n"
            "  suggest — recommend DADA2 params from a demux QZV\n"
            "  run     — detect primers then suggest DADA2 params in one step"
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

    return parser


# ===========================================================================
# SECTION 6 — Subcommand handlers
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


# ===========================================================================
# SECTION 7 — Entry point
# ===========================================================================

def main(argv: Optional[List[str]] = None) -> int:
    """
Entry point. Parses argv (or sys.argv if None) and dispatches to the
    appropriate subcommand handler.

    Returns the integer exit code (0 = success, 1 = error).
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.subcommand == "detect":
        return cmd_detect(args)
    elif args.subcommand == "suggest":
        return cmd_suggest(args)
    elif args.subcommand == "run":
        return cmd_run(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
