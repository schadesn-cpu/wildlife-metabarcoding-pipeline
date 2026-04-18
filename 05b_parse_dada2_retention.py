#!/usr/bin/env python3
"""
05b_parse_dada2_retention.py
========================
Diagnostic utility for DADA2 read retention across all markers in a
QIIME2 metabarcoding pipeline. Exports and parses denoising-stats.qza
files, computes per-step retention rates, flags problematic samples,
and identifies which step is losing reads.

This is especially useful for:
  - Diagnosing ITS/18S amplicons with variable-length inserts
    (poor merging = truncation parameters too strict)
  - Identifying samples with overall low retention (bad extraction,
    failed PCR, adapter contamination)
  - Comparing retention across markers in a multi-marker study

Usage
-----
  # Check all markers using default paths
  python scripts/utils/parse_dada2_retention.py

  # Check specific markers
  python scripts/utils/parse_dada2_retention.py --markers 16S MiFish cytb ITS

  # Point at specific stats files directly
  python scripts/utils/parse_dada2_retention.py \\
      --stats qiime2/dada2/stats_16S.qza \\
              qiime2/MiFish/all/dada2/denoising-stats.qza \\
              qiime2/cytb/all/dada2/denoising-stats.qza

  # Set custom thresholds
  python scripts/utils/parse_dada2_retention.py \\
      --min-retention 50 --min-merged 30

  # Write TSV report
  python scripts/utils/parse_dada2_retention.py \\
      --out results/qc/dada2_retention_report.tsv

Interpreting the output
-----------------------
The script identifies the BOTTLENECK step for each marker — the step
where the most reads are lost. Common patterns and their causes:

  BOTTLENECK: filtering
    - Quality scores too low → raise --p-max-ee
    - Reads too short for truncation → lower --p-trunc-len
    - Primer contamination remaining → recheck cutadapt step

  BOTTLENECK: merging
    - Most common cause: --p-trunc-len too aggressive for variable-
      length amplicons (ITS, 18S). Fix: set --p-trunc-len-f 0 and
      --p-trunc-len-r 0 to disable truncation.
    - Reads don't overlap: amplicon longer than 2×read length minus
      minimum overlap. Fix: increase read length or reduce amplicon.
    - Wrong orientation: try --p-read-orientation

  BOTTLENECK: chimeras
    - High chimera rate (>20%) suggests PCR over-amplification or
      very low-input samples. Usually acceptable at <10%.

  BOTTLENECK: none / good retention
    - All steps retained >80% → no action needed.

Adapting for other study systems
---------------------------------
The default marker search paths follow the pipeline directory
structure. For a different project, either:
  1. Use --stats to point directly at .qza files, or
  2. Use --markers with --project-root to auto-search, or
  3. Edit MARKER_SEARCH_PATHS below.

The thresholds (--min-retention, --min-merged) can be adjusted for
your system. Low-biomass samples (tick blood meal, museum specimens)
may have legitimately low input and lower acceptable thresholds.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Default marker search paths (relative to project root)
# Edit these for your project structure
# ---------------------------------------------------------------------------
MARKER_SEARCH_PATHS = {
    "16S":    ["qiime2/dada2/stats_16S.qza",
                "qiime2/16S/all/dada2/denoising-stats.qza"],
    "MiFish": ["qiime2/MiFish/all/dada2/denoising-stats.qza",
                "qiime2/dada2/stats_Mifish.qza"],
    "cytb":   ["qiime2/cytb/all/dada2/denoising-stats.qza",
                "qiime2/dada2/stats_cytb.qza"],
    "18S":    ["qiime2/18S/all/dada2/denoising-stats.qza",
                "qiime2/dada2/stats_18S.qza"],
    "ITS":    ["qiime2/ITS/all/dada2/denoising-stats.qza",
                "qiime2/dada2/stats_ITS1-2.qza"],
}

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

class SampleStats:
    def __init__(self, sample_id: str, row: dict):
        self.sample_id = sample_id
        self.input       = int(row.get("input", 0) or 0)
        self.filtered    = int(row.get("filtered", 0) or 0)
        self.denoised    = int(row.get("denoised", 0) or 0)
        self.merged      = int(row.get("merged", 0) or 0)
        self.non_chimeric = int(row.get("non-chimeric", 0) or 0)

    @property
    def pct_filtered(self) -> float:
        return 100 * self.filtered / self.input if self.input else 0

    @property
    def pct_merged(self) -> float:
        return 100 * self.merged / self.input if self.input else 0

    @property
    def pct_final(self) -> float:
        return 100 * self.non_chimeric / self.input if self.input else 0

    @property
    def bottleneck(self) -> str:
        """Identify which step lost the most reads."""
        if self.input == 0:
            return "no_input"
        lost_filter   = self.input - self.filtered
        lost_merge    = self.filtered - self.merged
        lost_chimera  = self.merged - self.non_chimeric
        steps = {
            "filtering": lost_filter,
            "merging":   lost_merge,
            "chimeras":  lost_chimera,
        }
        worst = max(steps, key=steps.get)
        # Only flag as bottleneck if >20% of input lost at that step
        if steps[worst] / self.input > 0.20:
            return worst
        return "none"

    def is_control(self) -> bool:
        sid = self.sample_id.upper()
        return any(sid.startswith(p) for p in ("NTC", "PAC", "XB-", "BLANK"))


# ---------------------------------------------------------------------------
# QIIME2 export
# ---------------------------------------------------------------------------

def export_stats_qza(qza_path: Path, tmpdir: Path) -> Path:
    """Export a denoising-stats.qza to a TSV file."""
    outdir = tmpdir / qza_path.stem
    outdir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        ["qiime", "tools", "export",
         "--input-path", str(qza_path),
         "--output-path", str(outdir)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  [warn] Export failed for {qza_path}: {result.stderr[:200]}",
              file=sys.stderr)
        return None
    tsv = outdir / "stats.tsv"
    if not tsv.exists():
        print(f"  [warn] stats.tsv not found in {outdir}", file=sys.stderr)
        return None
    return tsv


# ---------------------------------------------------------------------------
# TSV parsing
# ---------------------------------------------------------------------------

def parse_stats_tsv(tsv_path: Path) -> List[SampleStats]:
    """Parse a DADA2 denoising stats TSV into SampleStats objects."""
    samples = []
    with open(tsv_path, encoding="utf-8") as f:
        lines = f.readlines()

    # Find header line (skip #q2:types line)
    header_line = lines[0].strip()
    headers = [h.strip() for h in header_line.split("\t")]

    # Normalize header names — QIIME2 uses spaces, we want underscores
    def norm(h: str) -> str:
        return h.lower().replace(" ", "-").replace("_", "-")

    headers_norm = [norm(h) for h in headers]

    for line in lines[1:]:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        row = dict(zip(headers_norm, parts))
        sample_id = parts[0]
        try:
            s = SampleStats(sample_id, row)
            samples.append(s)
        except (ValueError, KeyError):
            continue

    return samples


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

COLORS = {
    "red":    "\033[91m",
    "yellow": "\033[93m",
    "green":  "\033[92m",
    "reset":  "\033[0m",
    "bold":   "\033[1m",
}

def colorize(text: str, color: str) -> str:
    if not sys.stdout.isatty():
        return text
    return f"{COLORS.get(color, '')}{text}{COLORS['reset']}"

def flag(pct: float, min_val: float) -> str:
    if pct < min_val * 0.5:
        return colorize("✗ FAIL", "red")
    elif pct < min_val:
        return colorize("⚠ WARN", "yellow")
    return colorize("✓", "green")


def print_marker_report(
    marker: str,
    samples: List[SampleStats],
    min_retention: float,
    min_merged: float,
    show_controls: bool,
) -> Tuple[int, int]:
    """Print report for one marker. Returns (n_warn, n_fail)."""

    bio_samples = [s for s in samples if not s.is_control()]
    controls    = [s for s in samples if s.is_control()]

    if not bio_samples:
        print(f"\n{colorize(marker, 'bold')}: no biological samples found")
        return 0, 0

    # Summary stats
    inputs   = [s.input for s in bio_samples]
    finals   = [s.pct_final for s in bio_samples]
    mergeds  = [s.pct_merged for s in bio_samples]

    avg_input  = sum(inputs) / len(inputs)
    min_input  = min(inputs)
    avg_final  = sum(finals) / len(finals)
    avg_merged = sum(mergeds) / len(mergeds)

    # Identify overall bottleneck
    bottlenecks = [s.bottleneck for s in bio_samples
                   if s.bottleneck != "none" and s.bottleneck != "no_input"]
    if bottlenecks:
        from collections import Counter
        most_common_bottleneck = Counter(bottlenecks).most_common(1)[0][0]
    else:
        most_common_bottleneck = "none"

    print(f"\n{'='*70}")
    print(f"  {colorize(marker, 'bold')}  ({len(bio_samples)} biological samples)")
    print(f"{'='*70}")
    print(f"  Input reads:    avg={avg_input:,.0f}   min={min_input:,}")
    print(f"  Avg merged:     {avg_merged:.1f}%  {flag(avg_merged, min_merged)}")
    print(f"  Avg final:      {avg_final:.1f}%  {flag(avg_final, min_retention)}")

    if most_common_bottleneck != "none":
        bottleneck_advice = {
            "filtering": "Check --p-max-ee and --p-trunc-len settings",
            "merging":   "Try --p-trunc-len-f 0 --p-trunc-len-r 0 (especially for ITS/18S)",
            "chimeras":  "High chimera rate — check PCR cycles and input DNA quality",
        }
        advice = bottleneck_advice.get(most_common_bottleneck, "")
        print(f"  Bottleneck:     {colorize(most_common_bottleneck.upper(), 'yellow')}"
              f"  →  {advice}")
    else:
        print(f"  Bottleneck:     {colorize('none — retention looks good', 'green')}")

    # Per-sample table
    n_warn = n_fail = 0
    problem_samples = [s for s in bio_samples
                       if s.pct_final < min_retention or s.pct_merged < min_merged]

    if problem_samples:
        print(f"\n  Samples below threshold (final<{min_retention}% or merged<{min_merged}%):")
        print(f"  {'Sample':<40} {'Input':>8} {'Filtered%':>10} "
              f"{'Merged%':>10} {'Final%':>8} {'Bottleneck':<12}")
        print(f"  {'-'*40} {'-'*8} {'-'*10} {'-'*10} {'-'*8} {'-'*12}")
        for s in sorted(problem_samples, key=lambda x: x.pct_final):
            final_f = flag(s.pct_final, min_retention)
            merged_f = flag(s.pct_merged, min_merged)
            print(f"  {s.sample_id:<40} {s.input:>8,} "
                  f"{s.pct_filtered:>9.1f}% {s.pct_merged:>9.1f}% "
                  f"{s.pct_final:>7.1f}% {s.bottleneck:<12}")
            if s.pct_final < min_retention * 0.5:
                n_fail += 1
            else:
                n_warn += 1
    else:
        print(f"\n  {colorize('All biological samples pass thresholds.', 'green')}")

    # Controls summary
    if show_controls and controls:
        print(f"\n  Controls ({len(controls)}):")
        for c in controls:
            print(f"    {c.sample_id:<40} input={c.input:>6,}  "
                  f"final={c.non_chimeric:>6,} ({c.pct_final:.1f}%)")
        ntc_reads = [c.non_chimeric for c in controls
                     if c.sample_id.upper().startswith("NTC")]
        if ntc_reads and max(ntc_reads) > 100:
            print(f"  {colorize('⚠ NTC has >100 reads after denoising — check for contamination', 'yellow')}")

    return n_warn, n_fail


def write_tsv_report(
    marker_results: Dict[str, List[SampleStats]],
    out_path: Path,
) -> None:
    """Write a flat TSV with all markers and samples."""
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("marker\tsample_id\tis_control\tinput\tfiltered\t"
                "pct_filtered\tmerged\tpct_merged\tnon_chimeric\t"
                "pct_final\tbottleneck\n")
        for marker, samples in sorted(marker_results.items()):
            for s in samples:
                f.write(
                    f"{marker}\t{s.sample_id}\t{s.is_control()}\t"
                    f"{s.input}\t{s.filtered}\t{s.pct_filtered:.1f}\t"
                    f"{s.merged}\t{s.pct_merged:.1f}\t"
                    f"{s.non_chimeric}\t{s.pct_final:.1f}\t"
                    f"{s.bottleneck}\n"
                )
    print(f"\nReport written to: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="parse_dada2_retention.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--markers", nargs="+",
        default=list(MARKER_SEARCH_PATHS.keys()),
        help="Markers to check. Default: all known markers."
    )
    p.add_argument(
        "--stats", nargs="+", type=Path, default=None,
        help="Direct paths to denoising-stats.qza files. "
             "Overrides --markers auto-search."
    )
    p.add_argument(
        "--project-root", type=Path, default=Path("."),
        help="Project root directory. Default: current directory."
    )
    p.add_argument(
        "--min-retention", type=float, default=60.0,
        help="Minimum acceptable %% of input reads in final output. "
             "Default: 60%%."
    )
    p.add_argument(
        "--min-merged", type=float, default=50.0,
        help="Minimum acceptable %% of input reads passing merging. "
             "Default: 50%%."
    )
    p.add_argument(
        "--show-controls", action="store_true", default=False,
        help="Show control samples (NTC, PAC, XB) in output."
    )
    p.add_argument(
        "--out", type=Path, default=None,
        help="Write TSV report to this path."
    )
    return p


def main() -> int:
    args = build_parser().parse_args()
    project_root = args.project_root.resolve()
    os.chdir(project_root)

    # Collect QZA paths
    qza_paths: Dict[str, Path] = {}

    if args.stats:
        # Direct paths provided
        for p in args.stats:
            qza_paths[p.stem] = p
    else:
        # Auto-search by marker name
        for marker in args.markers:
            candidates = MARKER_SEARCH_PATHS.get(marker, [])
            found = None
            for c in candidates:
                cp = project_root / c
                if cp.exists():
                    found = cp
                    break
            if found:
                qza_paths[marker] = found
            else:
                print(f"  [skip] {marker}: no stats QZA found "
                      f"(searched {len(candidates)} locations)",
                      file=sys.stderr)

    if not qza_paths:
        print("No DADA2 stats files found. Use --stats to specify paths directly.",
              file=sys.stderr)
        return 2

    print(f"\nDADA2 Retention Diagnostic")
    print(f"Project: {project_root}")
    print(f"Markers: {', '.join(qza_paths.keys())}")
    print(f"Thresholds: final≥{args.min_retention}%  merged≥{args.min_merged}%")

    marker_results: Dict[str, List[SampleStats]] = {}
    total_warn = total_fail = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        for marker, qza_path in sorted(qza_paths.items()):
            print(f"\n  Exporting {marker} stats from {qza_path}...",
                  file=sys.stderr, end=" ")
            tsv = export_stats_qza(qza_path, Path(tmpdir))
            if tsv is None:
                print("FAILED", file=sys.stderr)
                continue
            print("done", file=sys.stderr)

            samples = parse_stats_tsv(tsv)
            if not samples:
                print(f"  [warn] No samples parsed from {tsv}", file=sys.stderr)
                continue

            marker_results[marker] = samples
            n_warn, n_fail = print_marker_report(
                marker, samples,
                args.min_retention, args.min_merged,
                args.show_controls,
            )
            total_warn += n_warn
            total_fail += n_fail

    # Overall summary
    print(f"\n{'='*70}")
    print(f"  SUMMARY across {len(marker_results)} markers")
    print(f"{'='*70}")
    if total_fail > 0:
        print(colorize(f"  ✗ {total_fail} samples FAIL (< {args.min_retention*0.5:.0f}% retention)", "red"))
    if total_warn > 0:
        print(colorize(f"  ⚠ {total_warn} samples WARN (< {args.min_retention:.0f}% retention)", "yellow"))
    if total_fail == 0 and total_warn == 0:
        print(colorize("  ✓ All samples pass retention thresholds", "green"))

    print(f"\n  Bottleneck guide:")
    print(f"    filtering → raise --p-max-ee or check truncation lengths")
    print(f"    merging   → set --p-trunc-len-f 0 --p-trunc-len-r 0 (ITS/18S)")
    print(f"    chimeras  → check PCR cycle count and input DNA quality")

    if args.out:
        write_tsv_report(marker_results, args.out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
