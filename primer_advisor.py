#!/usr/bin/env python3
"""
primer_advisor.py
=================
Step: primer_advisor

Purpose:
    Identify the PCR primers in your raw R1/R2 reads, per marker, by scoring a
    reference primer database against sampled reads (IUPAC-aware, mismatch
    tolerant). Reports the best-matching forward/reverse pair, the amplicon
    region, and a ready-to-run cutadapt command, and writes a TSV the rest of
    the pipeline reads.

    This is the "detect" stage. Two related jobs that an older combined script
    also did now live in their own ported steps: DADA2 truncation suggestions
    are in dada2_advisor.py, and post-cutadapt read-count / dimer QC is in
    parse_demux.py. This script does not duplicate them.

Inputs:
    reads/<marker>/   raw paired FASTQ.gz (the merge_runs / manifest layout)
    pipeline_config.yml   active_markers, samples.control_prefixes, reads/ location

Outputs:
    reports/primers_detected.tsv   per-marker forward/reverse primer calls
    stdout                          per-marker detail + recommended cutadapt command
    logs/run_manifest.jsonl         run appended on completion

Runs before cutadapt, directly on raw reads. Control samples (config
control_prefixes) are skipped so empty negatives don't depress match scores.

Usage:
    python primer_advisor.py                 # all active markers from config
    python primer_advisor.py --marker MiFish # a single marker
    python primer_advisor.py --all --reads-dir reads/

Requirements:
    Python >= 3.8 (standard library only for detection; pyyaml via config_loader)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import logging
import sys
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# --- make config_loader and the utils package importable regardless of cwd ---
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from config_loader import load_config, get_paths            # noqa: E402
from utils import checkpoint, provenance                    # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ===========================================================================
# Primer database, IUPAC matching, read sampling, and detection
# (carried verbatim from the v1 primer advisor — the reference primer set and
#  detection logic are unchanged; only control-sample handling is now wired to
#  config, see get_sample_reads / detect_primers_for_marker.)
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
                     n: int = 200,
                     control_prefixes: Tuple[str, ...] = ()) -> Tuple[List[str], List[str]]:
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

    _ctrl = tuple(pfx.lower() for pfx in control_prefixes)
    r1_seqs: List[str] = []
    r2_seqs: List[str] = []

    for f in sorted(marker_dir.iterdir()):
        if f.name.startswith(".") or any(f.name.lower().startswith(p) for p in _ctrl):
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
        if f.name.startswith(".") or any(f.name.lower().startswith(p) for p in _ctrl):
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
    control_prefixes: Tuple[str, ...] = (),
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
        r1_seqs, r2_seqs = get_sample_reads(reads_dir, marker, n=n_reads,
                                            control_prefixes=control_prefixes)
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
            f"    --i-demultiplexed-sequences qiime2/{marker}/all/imported/demux.qza \\\n"
            f"    --p-front-f {fwd_seq} \\\n"
            f"    --p-front-r {rev_seq} \\\n"
            f"    --p-adapter-f {rc(rev_seq)} \\\n"
            f"    --p-adapter-r {rc(fwd_seq)} \\\n"
            f"    --p-discard-untrimmed \\\n"
            f"    --p-error-rate 0.2 \\\n"
            f"    --o-trimmed-sequences qiime2/{marker}/all/imported/demux_trimmed.qza \\\n"
            f"    --verbose"
        )
    else:
        cutadapt_cmd = (
            f"qiime cutadapt trim-paired \\\n"
            f"    --i-demultiplexed-sequences qiime2/{marker}/all/imported/demux.qza \\\n"
            f"    --p-front-f {fwd_seq} \\\n"
            f"    --p-front-r {rev_seq} \\\n"
            f"    --p-discard-untrimmed \\\n"
            f"    --p-error-rate 0.1 \\\n"
            f"    --o-trimmed-sequences qiime2/{marker}/all/imported/demux_trimmed.qza \\\n"
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
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore", delimiter="\t")
        w.writeheader()
        w.writerows(results)
    print(f"\nReport written to: {out_path}")


def discover_markers(reads_dir: Path) -> List[str]:
    """List all marker subdirectories in the reads directory."""
    return [d.name for d in sorted(reads_dir.iterdir())
            if d.is_dir() and not d.name.startswith(".")]


# ===========================================================================
# CLI
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="primer_advisor.py",
        description="Identify PCR primers in raw reads and recommend a cutadapt command.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    group = p.add_mutually_exclusive_group()
    group.add_argument("--marker", default=None,
                       help="Single marker to check (e.g. MiFish, 16S).")
    group.add_argument("--all", action="store_true",
                       help="Check every marker subdirectory found in --reads-dir.")
    p.add_argument("--reads-dir", default=None, type=Path,
                   help="Root reads directory with per-marker subdirectories. "
                        "Default: 'reads/' under the project root.")
    p.add_argument("--config", default=None, help="Path to pipeline_config.yml.")
    p.add_argument("--report", default=None, type=Path,
                   help="Path to write the TSV summary. "
                        "Default: reports/primers_detected.tsv under the project root.")
    p.add_argument("--n-reads", type=int, default=200,
                   help="Reads to sample per marker. Default: 200")
    p.add_argument("--mismatches", type=int, default=2,
                   help="Allowed mismatches when matching primers. Default: 2")
    p.add_argument("--verbose", action="store_true",
                   help="Show all primer scores, not just the best match.")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    cfg = load_config(args.config)
    paths = get_paths(cfg)

    reads_dir = (args.reads_dir or cfg.resolve("reads")).resolve()
    if not reads_dir.is_dir():
        log.error("Reads directory not found: %s", reads_dir)
        return 1

    control_prefixes = tuple(cfg.samples.get("control_prefixes", []))

    # Which markers: --all discovers subdirectories; --marker is one; default is
    # the project's active_markers from config.
    if args.all:
        markers = discover_markers(reads_dir)
        log.info("Discovered markers: %s", ", ".join(markers) if markers else "(none)")
    elif args.marker:
        markers = [args.marker]
    else:
        markers = list(cfg.active_markers)
        log.info("No --marker/--all given; using active_markers: %s", ", ".join(markers))

    if not markers:
        log.error("No markers to check (empty --reads-dir and no active_markers).")
        return 1

    results = []
    for marker in markers:
        results.append(detect_primers_for_marker(
            marker=marker, reads_dir=reads_dir,
            n_reads=args.n_reads, mismatch=args.mismatches,
            control_prefixes=control_prefixes, verbose=args.verbose,
        ))

    report_path = args.report or (cfg.resolve("reports") / "primers_detected.tsv")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    write_detect_report(results, report_path)

    detected = [r for r in results if "error" not in r]
    checkpoint.print_checkpoint(
        cfg, "primer_advisor",
        produced=[report_path],
        provenance={
            "inputs": {"reads_dir": reads_dir},
            "outputs": [report_path],
            "command": "python " + " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
            "extra": {"markers": markers, "detected": len(detected)},
        },
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
