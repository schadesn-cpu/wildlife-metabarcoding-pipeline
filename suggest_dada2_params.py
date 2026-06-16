#!/usr/bin/env python3
"""
suggest_dada2_params.py
=======================

Read a QIIME2 demux QZV (ideally a demux_trimmed.qzv — i.e., post-cutadapt)
and recommend DADA2 truncation parameters that won't silently filter your
data to zero.

Combines TWO separate checks:

    1. QUALITY CHECK: find the last position where median quality >= Qmin
    2. LENGTH CHECK: find the last position where >= length_pct of reads
       still exist (i.e., reads haven't ended)

Returns the STRICTER (shorter) of the two checks per direction, with safety
margin. Also validates that forward + reverse truncations together cover the
expected amplicon with enough overlap for DADA2 merging.

Why this script exists
----------------------
DADA2 is known for silently filtering every read to zero if the truncation
length exceeds actual read length. Looking only at quality scores (the
conventional approach) misses this entirely: a 150-bp read with quality Q37
throughout will still be rejected if you pass --trunc-len-f 200. This script
checks BOTH conditions so you catch the failure mode before spending hours
on DADA2.

Pipeline context
----------------
Designed for the wildlife-metabarcoding-pipeline. See:
    https://github.com/schadesn-cpu/wildlife-metabarcoding-pipeline

Usage
-----
    # Basic: just get suggestions
    python suggest_dada2_params.py \
        --qzv qiime2/16S/all/imported/demux_trimmed.qzv

    # Specify amplicon size for overlap checking
    python suggest_dada2_params.py \
        --qzv qiime2/16S/all/imported/demux_trimmed.qzv \
        --amplicon-length 253

    # Tighter quality threshold
    python suggest_dada2_params.py \
        --qzv qiime2/cytb/all/imported/demux_trimmed.qzv \
        --amplicon-length 307 \
        --min-quality 30

    # JSON output for pipeline integration
    python suggest_dada2_params.py \
        --qzv qiime2/16S/all/imported/demux_trimmed.qzv \
        --out-json qiime2/16S/all/imported/dada2_params.json

Exit codes
----------
    0: Suggestions produced successfully (may or may not be safe; see warnings)
    1: Fatal error — QZV could not be parsed or read
    2: Suggestions produced but overlap is negative (merging will fail);
       the caller may want to treat this as a failure
"""

import argparse
import json
import math
import sys
import zipfile
from pathlib import Path


# ---------- Argument parsing ----------

def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Recommend DADA2 truncation parameters based on read quality and length.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--qzv", required=True, type=Path,
                    help="Path to demux QZV (preferably demux_trimmed.qzv, post-cutadapt).")
    ap.add_argument("--amplicon-length", type=int, default=None,
                    help="Expected amplicon length in bp (after primer removal). "
                         "Used to check forward+reverse overlap is sufficient. "
                         "Common values: 16S V4 = 253, 18S V4 = 253, "
                         "cytb Kocher = ~307.")
    ap.add_argument("--min-overlap", type=int, default=20,
                    help="Minimum base-pair overlap required between merged F/R reads (default: 20).")
    ap.add_argument("--min-quality", type=int, default=25,
                    help="Minimum median quality score for quality-based truncation (default: 25). "
                         "Use 30 for stricter, 20 for more permissive.")
    ap.add_argument("--length-pct", type=float, default=95.0,
                    help="Required percentage of reads that must still exist at the truncation "
                         "position (default: 95). Reads shorter than trunc-len are discarded by DADA2.")
    ap.add_argument("--safety-margin", type=int, default=3,
                    help="Subtract this many bp from the raw quality/length cutoff as a safety "
                         "buffer (default: 3). Prevents edge effects near the cliff.")
    ap.add_argument("--out-json", type=Path, default=None,
                    help="Optional: write recommendations as JSON to this path.")
    ap.add_argument("--quiet", action="store_true",
                    help="Suppress the human-readable report (JSON still written if requested).")
    return ap.parse_args(argv)


# ---------- QZV parsing ----------

def load_quality_data(qzv_path, direction):
    """
    Extract the seven-number-summary table from a demux QZV for one read direction.

    Returns: (positions, summaries) where:
        positions = list of int read positions (1-indexed)
        summaries = dict with keys 'count', '2%', '9%', '25%', '50%', '75%', '91%', '98%'
                    each mapping to a list of floats aligned to positions
    """
    target_suffix = f"{direction}-seven-number-summaries.tsv"
    target_suffix_csv = f"{direction}-seven-number-summaries.csv"

    with zipfile.ZipFile(qzv_path) as zf:
        target_name = None
        for name in zf.namelist():
            if name.endswith(target_suffix) or name.endswith(target_suffix_csv):
                target_name = name
                break
        if target_name is None:
            raise ValueError(
                f"Could not find {target_suffix} in {qzv_path}. "
                f"Is this a valid demux QZV?"
            )
        with zf.open(target_name) as fh:
            content = fh.read().decode("utf-8")

    sep = "\t" if target_name.endswith(".tsv") else ","
    lines = content.strip().split("\n")
    header = lines[0].split(sep)

    # First column is the label ("count", "2%", etc); the rest are position indices
    positions = []
    for h in header[1:]:
        h = h.strip().strip('"')
        if h.isdigit():
            positions.append(int(h))

    n_pos = len(positions)
    summaries = {}
    for line in lines[1:]:
        parts = line.split(sep)
        label = parts[0].strip().strip('"')
        if label in ("count", "2%", "9%", "25%", "50%", "75%", "91%", "98%"):
            values = []
            for v in parts[1:1 + n_pos]:
                v = v.strip().strip('"')
                try:
                    values.append(float(v))
                except ValueError:
                    values.append(float("nan"))
            summaries[label] = values

    return positions, summaries


# ---------- Analysis ----------

def quality_cutoff(positions, medians, min_q):
    """
    Find the last position where median quality stays >= min_q.
    If quality never drops, returns the max position.
    """
    last_good = 0
    for pos, q in zip(positions, medians):
        if q >= min_q:
            last_good = pos
        else:
            # Stop at first sustained drop
            break
    return last_good


def length_cutoff(positions, counts, length_pct):
    """
    Find the last position where at least length_pct% of reads are present.
    Required because DADA2 silently drops reads shorter than trunc_len.
    """
    if not counts:
        return 0
    max_count = max(counts)
    threshold = (length_pct / 100.0) * max_count
    last_good = 0
    for pos, c in zip(positions, counts):
        if c >= threshold:
            last_good = pos
        else:
            break
    return last_good


def analyze_direction(positions, summaries, min_q, length_pct, safety_margin):
    """
    Run both quality and length checks and return analysis dict.
    """
    medians = summaries.get("50%", [])
    counts = summaries.get("count", [])

    q_cutoff = quality_cutoff(positions, medians, min_q)
    l_cutoff = length_cutoff(positions, counts, length_pct)

    # The constraining cutoff is the stricter (shorter) of the two
    raw_cutoff = min(q_cutoff, l_cutoff)
    recommended = max(0, raw_cutoff - safety_margin)

    # Also compute the max position for context
    max_pos = max(positions) if positions else 0
    max_count = max(counts) if counts else 0

    # Quality at the cutoff for reporting
    q_at_cutoff = None
    c_at_cutoff = None
    for pos, q, c in zip(positions, medians, counts):
        if pos == recommended:
            q_at_cutoff = q
            c_at_cutoff = c
            break

    return {
        "max_position": max_pos,
        "max_read_count": int(max_count),
        "quality_cutoff": q_cutoff,
        "length_cutoff": l_cutoff,
        "raw_cutoff": raw_cutoff,
        "recommended_truncation": recommended,
        "constraining_factor": "length" if l_cutoff < q_cutoff else "quality",
        "median_quality_at_recommendation": q_at_cutoff,
        "read_count_at_recommendation": int(c_at_cutoff) if c_at_cutoff is not None else None,
    }


def check_overlap(fwd_trunc, rev_trunc, amplicon_length, min_overlap):
    """
    Check that forward + reverse truncations will leave enough overlap for DADA2 merging.

    Returns dict with overlap info.
    """
    total = fwd_trunc + rev_trunc
    overlap = total - amplicon_length if amplicon_length else None

    status = "ok"
    if overlap is not None:
        if overlap < 0:
            status = "fatal_no_overlap"
        elif overlap < min_overlap:
            status = "warning_insufficient_overlap"

    return {
        "forward_truncation": fwd_trunc,
        "reverse_truncation": rev_trunc,
        "combined_length": total,
        "amplicon_length": amplicon_length,
        "overlap_bp": overlap,
        "min_required": min_overlap,
        "status": status,
    }


# ---------- Reporting ----------

def print_human_report(fwd, rev, overlap_info, qzv_path, min_q, length_pct):
    print(f"\n{'=' * 70}")
    print(f"  DADA2 TRUNCATION RECOMMENDATIONS")
    print(f"{'=' * 70}")
    print(f"  Source QZV:      {qzv_path}")
    print(f"  Min median Q:    {min_q}")
    print(f"  Min read %:      {length_pct}%")
    print()

    for direction, info in [("FORWARD (R1)", fwd), ("REVERSE (R2)", rev)]:
        print(f"  {direction}:")
        print(f"    Read length:              up to position {info['max_position']}")
        print(f"    Quality-based cutoff:     {info['quality_cutoff']} (last pos where Q >= {min_q})")
        print(f"    Length-based cutoff:      {info['length_cutoff']} (last pos where >= {length_pct}% reads present)")
        print(f"    Constraining factor:      {info['constraining_factor']}")
        print(f"    Recommended truncation:   {info['recommended_truncation']}")
        if info['median_quality_at_recommendation'] is not None:
            print(f"    At recommended position:  Q={info['median_quality_at_recommendation']:.1f}, "
                  f"{info['read_count_at_recommendation']} reads")
        print()

    if overlap_info["amplicon_length"]:
        print(f"  OVERLAP CHECK:")
        print(f"    Amplicon length:      {overlap_info['amplicon_length']} bp")
        print(f"    Forward + Reverse:    {overlap_info['combined_length']} bp "
              f"({fwd['recommended_truncation']} + {rev['recommended_truncation']})")
        print(f"    Available overlap:    {overlap_info['overlap_bp']} bp "
              f"(need >= {overlap_info['min_required']})")

        if overlap_info["status"] == "fatal_no_overlap":
            print(f"    *** FATAL: Forward + reverse truncations are SHORTER than the amplicon.")
            print(f"    ***        DADA2 merging will fail. Increase --trunc-len-f or --trunc-len-r")
            print(f"    ***        or use single-end denoising.")
        elif overlap_info["status"] == "warning_insufficient_overlap":
            print(f"    *** WARNING: Overlap is below minimum ({overlap_info['overlap_bp']} < {overlap_info['min_required']} bp).")
            print(f"    ***          Merging may fail for reads with indels or variants.")
        else:
            print(f"    Status:               OK")
    else:
        print(f"  OVERLAP CHECK: skipped (--amplicon-length not provided)")

    print()
    print(f"  --- SUGGESTED DADA2 COMMAND FRAGMENT ---")
    print(f"  --trunc-len-f {fwd['recommended_truncation']} --trunc-len-r {rev['recommended_truncation']}")
    print()


# ---------- Main ----------

def main(argv=None):
    args = parse_args(argv)

    if not args.qzv.exists():
        print(f"Error: QZV file not found: {args.qzv}", file=sys.stderr)
        return 1

    try:
        fwd_positions, fwd_summaries = load_quality_data(args.qzv, "forward")
        rev_positions, rev_summaries = load_quality_data(args.qzv, "reverse")
    except (ValueError, zipfile.BadZipFile) as e:
        print(f"Error reading QZV: {e}", file=sys.stderr)
        return 1

    fwd = analyze_direction(fwd_positions, fwd_summaries, args.min_quality, args.length_pct, args.safety_margin)
    rev = analyze_direction(rev_positions, rev_summaries, args.min_quality, args.length_pct, args.safety_margin)

    overlap_info = check_overlap(
        fwd["recommended_truncation"],
        rev["recommended_truncation"],
        args.amplicon_length,
        args.min_overlap,
    )

    if not args.quiet:
        print_human_report(fwd, rev, overlap_info, args.qzv, args.min_quality, args.length_pct)

    if args.out_json is not None:
        result = {
            "qzv_path": str(args.qzv),
            "parameters": {
                "min_quality": args.min_quality,
                "length_pct": args.length_pct,
                "safety_margin": args.safety_margin,
                "min_overlap": args.min_overlap,
                "amplicon_length": args.amplicon_length,
            },
            "forward": fwd,
            "reverse": rev,
            "overlap": overlap_info,
            "dada2_args": {
                "trunc_len_f": fwd["recommended_truncation"],
                "trunc_len_r": rev["recommended_truncation"],
            },
        }
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.out_json, "w") as fh:
            json.dump(result, fh, indent=2)
        if not args.quiet:
            print(f"  Wrote JSON: {args.out_json}")

    # Exit code: 2 if overlap is fatal, 0 otherwise
    if overlap_info["status"] == "fatal_no_overlap":
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
