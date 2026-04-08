#!/usr/bin/env python3
"""
14_viral_stats.py
=================
Compute detection statistics for viral amplicon data produced by the
loon metabarcoding pipeline. Designed to be run alongside 10_plot_viral.py
using identical parameters so figures and stats are always consistent.

Supports two signal modes (same as 10_plot_viral.py):
  - No-hit mode (default): signal = reads with no BLAST match
    Use for herpesvirus (TGF-IYG pan-herpesvirus primers)
  - Taxon-filter mode: signal = reads matching a specific BLAST taxon
    Use for adenovirus (--taxon-filter Aviadenovirus)

Outputs
-------
  results/{marker}/stats/{stem}_stats.tsv   — detection rate table
  results/{marker}/stats/{stem}_fisher.tsv  — pairwise Fisher's exact results
  Console: formatted summary table

Usage examples
--------------
  # Herpesvirus (no-hit signal, relative abundance threshold)
  python scripts/14_viral_stats.py \\
      --xlsx loon_amplicon_analysis.xlsx \\
      --sheet herpes \\
      --metadata metadata/qiime/metadata_MiFish.tsv \\
      --threshold 0.01 \\
      --outdir results/herpesvirus/stats/

  # Adenovirus (taxon-filtered signal, absolute read count threshold)
  python scripts/14_viral_stats.py \\
      --xlsx loon_amplicon_analysis.xlsx \\
      --sheet adeno \\
      --metadata metadata/qiime/metadata_MiFish.tsv \\
      --taxon-filter Aviadenovirus \\
      --min-reads 10 \\
      --group-order Diseased Trauma Marine \\
      --outdir results/adenovirus/stats/
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from itertools import combinations
from typing import List, Optional

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

BLAST_COL  = "BLAST nt"
NO_HIT_TAG = "no-hit"
# ADAPT FOR YOUR STUDY: This regex extracts the short sample ID from the
# full sample name in your metadata. For loons, IDs look like TV230084.
# Change this pattern to match your own sample naming convention.
TV_REGEX   = r"(TV\d+)"


# ── Data loading ──────────────────────────────────────────────────────────────

def load_sheet(xlsx: Path, sheet: str) -> pd.DataFrame:
    """Load BLAST amplicon sheet from Excel workbook."""
    df = pd.read_excel(xlsx, sheet_name=sheet, header=0)
    if BLAST_COL not in df.columns:
        raise ValueError(f"Column '{BLAST_COL}' not found in sheet '{sheet}'.")
    return df


def load_metadata(path: Path, group_col: str) -> pd.DataFrame:
    """Load QIIME2 metadata TSV and extract TV ID + group column."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[~df.iloc[:, 0].str.startswith("#", na=False)].reset_index(drop=True)
    df["_TV"] = df.iloc[:, 0].str.extract(TV_REGEX, expand=False)
    if group_col not in df.columns:
        raise ValueError(f"Group column '{group_col}' not found in metadata.")
    return df[["_TV", group_col]].dropna(subset=["_TV"])


def compute_per_sample(
    df: pd.DataFrame,
    taxon_filter: Optional[str] = None,
) -> pd.DataFrame:
    """
    Compute per-sample signal read counts and relative abundance.
    Signal = no-hit reads (herpesvirus) or taxon-filtered reads (adenovirus).
    """
    sample_cols = [c for c in df.columns if str(c).startswith("TV")]

    if taxon_filter:
        signal_mask = df[BLAST_COL].str.contains(taxon_filter, case=False, na=False)
        n_match = int(signal_mask.sum())
        log.info("Taxon filter '%s': %d / %d OTUs match", taxon_filter, n_match, len(df))
        if n_match == 0:
            log.warning("No OTUs matched taxon filter '%s' — check spelling.", taxon_filter)
    else:
        signal_mask = df[BLAST_COL].str.startswith(NO_HIT_TAG, na=False)
        log.info("No-hit mode: %d / %d OTUs are no-hit", signal_mask.sum(), len(df))

    total  = df[sample_cols].sum()
    signal = df.loc[signal_mask, sample_cols].sum()

    result = pd.DataFrame({
        "total_reads":   total,
        "signal_reads":  signal,
        "signal_relabund": signal / total.replace(0, np.nan),
    })
    result.index = pd.Series(result.index).str.extract(TV_REGEX, expand=False).values
    result.index.name = "TV"
    return result


def call_detection(
    per_sample: pd.DataFrame,
    min_reads: Optional[int],
    threshold: float,
    min_sample_reads: int,
) -> pd.DataFrame:
    """
    Add 'detected' boolean column to per-sample DataFrame.
    Drops samples below min_sample_reads total depth first.
    """
    df = per_sample.copy()

    # Drop low-depth samples
    if min_sample_reads > 0:
        before = len(df)
        df = df[df["total_reads"] >= min_sample_reads]
        dropped = before - len(df)
        if dropped:
            log.warning("Dropped %d samples with < %d total reads", dropped, min_sample_reads)

    # Detection calling
    if min_reads is not None:
        df["detected"] = df["signal_reads"] >= min_reads
        log.info("Detection threshold: >= %d reads (absolute)", min_reads)
    else:
        df["detected"] = df["signal_relabund"] >= threshold
        log.info("Detection threshold: >= %.1f%% relative abundance", threshold * 100)

    return df


# ── Statistics ────────────────────────────────────────────────────────────────

def detection_rates(
    per_sample: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    group_order: Optional[List[str]],
) -> pd.DataFrame:
    """Build per-group detection rate table."""
    merged = per_sample.reset_index().merge(
        meta.rename(columns={"_TV": "TV"}),
        on="TV", how="inner"
    )
    if merged.empty:
        raise RuntimeError("No samples matched between Excel and metadata. Check TV ID format.")

    groups = group_order or sorted(merged[group_col].dropna().unique().tolist())
    rows = []
    for g in groups:
        grp = merged[merged[group_col] == g]
        n_total = len(grp)
        n_det   = int(grp["detected"].sum())
        pct     = (n_det / n_total * 100) if n_total > 0 else 0.0
        rows.append({
            "group":       g,
            "n_total":     n_total,
            "n_detected":  n_det,
            "n_not_detected": n_total - n_det,
            "pct_detected": round(pct, 1),
        })
    return pd.DataFrame(rows)


def run_fisher(rates: pd.DataFrame) -> pd.DataFrame:
    """
    Run pairwise Fisher's exact tests (two-sided) for all group pairs.
    Also runs chi-square for the overall 3+ group comparison.
    """
    rows = []
    groups = rates["group"].tolist()

    for g1, g2 in combinations(groups, 2):
        r1 = rates[rates["group"] == g1].iloc[0]
        r2 = rates[rates["group"] == g2].iloc[0]

        # 2x2 contingency: [[det1, notdet1], [det2, notdet2]]
        table = [
            [int(r1["n_detected"]), int(r1["n_not_detected"])],
            [int(r2["n_detected"]), int(r2["n_not_detected"])],
        ]
        _, p = fisher_exact(table, alternative="two-sided")

        sig = "***" if p <= 0.001 else "**" if p <= 0.01 else "*" if p <= 0.05 else "ns"
        rows.append({
            "comparison":  f"{g1} vs {g2}",
            "group_1":     g1,
            "group_2":     g2,
            "detected_1":  int(r1["n_detected"]),
            "total_1":     int(r1["n_total"]),
            "detected_2":  int(r2["n_detected"]),
            "total_2":     int(r2["n_total"]),
            "p_value":     round(p, 4),
            "sig":         sig,
            "test":        "Fisher exact (two-sided)",
        })

    # Overall chi-square (if 3+ groups)
    if len(groups) >= 3:
        contingency = rates[["n_detected", "n_not_detected"]].values
        # Check minimum expected cell counts
        chi2, p_overall, dof, expected = chi2_contingency(contingency)
        min_expected = expected.min()
        caveat = "" if min_expected >= 5 else f" [CAUTION: min expected cell = {min_expected:.1f} < 5]"
        rows.append({
            "comparison":  "Overall (all groups)",
            "group_1":     "",
            "group_2":     "",
            "detected_1":  "",
            "total_1":     "",
            "detected_2":  "",
            "total_2":     "",
            "p_value":     round(p_overall, 4),
            "sig":         ("***" if p_overall <= 0.001 else "**" if p_overall <= 0.01
                            else "*" if p_overall <= 0.05 else "ns") + caveat,
            "test":        f"Chi-square (df={dof}){caveat}",
        })

    return pd.DataFrame(rows)


# ── Output formatting ─────────────────────────────────────────────────────────

def print_summary(
    rates: pd.DataFrame,
    fisher: pd.DataFrame,
    sheet: str,
    taxon_filter: Optional[str],
    min_reads: Optional[int],
    threshold: float,
) -> None:
    """Print a clean summary table to stdout."""
    signal_desc = (
        f"taxon = '{taxon_filter}', >= {min_reads} reads"
        if taxon_filter and min_reads
        else f"taxon = '{taxon_filter}'"
        if taxon_filter
        else f">= {threshold * 100:.0f}% relative abundance (no-hit)"
    )

    print()
    print("=" * 60)
    print(f"  Viral detection stats — sheet: {sheet}")
    print(f"  Signal: {signal_desc}")
    print("=" * 60)
    print()
    print("Detection rates:")
    print(f"  {'Group':<20} {'Detected':>10} {'Total':>8} {'%':>8}")
    print(f"  {'-'*20} {'-'*10} {'-'*8} {'-'*8}")
    for _, row in rates.iterrows():
        print(f"  {row['group']:<20} {row['n_detected']:>5}/{row['n_total']:<4} "
              f"{row['n_total']:>8} {row['pct_detected']:>7.1f}%")
    print()
    print("Pairwise Fisher's exact tests (two-sided):")
    print(f"  {'Comparison':<30} {'p-value':>10} {'Sig':>6}")
    print(f"  {'-'*30} {'-'*10} {'-'*6}")
    for _, row in fisher.iterrows():
        print(f"  {row['comparison']:<30} {row['p_value']:>10.4f} {row['sig']:>6}")
    print()


# ── CLI ───────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="14_viral_stats.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    req = p.add_argument_group("required arguments")
    req.add_argument("--xlsx",     required=True, type=Path,
                     help="Excel workbook (loon_amplicon_analysis.xlsx).")
    req.add_argument("--sheet",    required=True,
                     help="Sheet name (e.g. herpes, adeno).")
    req.add_argument("--metadata", required=True, type=Path,
                     help="QIIME2 metadata TSV.")

    opt = p.add_argument_group("optional arguments")
    opt.add_argument("--group-by",    default="Group",
                     help="Metadata column for grouping. Default: Group.")
    opt.add_argument("--group-order", nargs="+", default=None,
                     help="Explicit group order. Default: alphabetical.")
    opt.add_argument("--taxon-filter", default=None, metavar="TAXON",
                     help="Filter signal to OTUs whose BLAST nt contains TAXON "
                          "(e.g. Aviadenovirus). Default: None (uses no-hit reads).")
    opt.add_argument("--threshold",   type=float, default=0.01,
                     help="Relative abundance threshold for detection (0-1). "
                          "Default: 0.01 (1%%). Ignored if --min-reads is set.")
    opt.add_argument("--min-reads",   type=int, default=None,
                     help="Absolute read count threshold. Overrides --threshold.")
    opt.add_argument("--min-sample-reads", type=int, default=500, metavar="N",
                     help="Drop samples with fewer than N total reads. Default: 500.")
    opt.add_argument("--outdir",      type=Path, default=Path("."),
                     help="Output directory. Default: current directory.")
    opt.add_argument("--output-stem", default=None,
                     help="Output filename stem. Default: {sheet}_stats.")
    opt.add_argument("--dry-run", action="store_true",
                     help="Print stats to console only; do not write files.")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    # Check input files exist before doing any work
    if not args.xlsx.exists():
        log.error("Excel file not found: %s", args.xlsx)
        log.error("  Check --xlsx points to your analysis workbook.")
        return 1
    if not args.metadata.exists():
        log.error("Metadata file not found: %s", args.metadata)
        log.error("  Check --metadata points to your QIIME2 metadata TSV.")
        return 1

    # Load data
    log.info("Loading sheet '%s' from %s", args.sheet, args.xlsx)
    try:
        df = load_sheet(args.xlsx, args.sheet)
    except ValueError as e:
        log.error("%s", e)
        log.error("  Available sheets can be checked by opening the workbook in Excel.")
        return 1
    except Exception as e:
        log.error("Failed to read Excel file: %s", e)
        return 1

    log.info("Loading metadata from %s", args.metadata)
    try:
        meta = load_metadata(args.metadata, args.group_by)
    except ValueError as e:
        log.error("%s", e)
        return 1

    # Compute signal
    per_sample = compute_per_sample(df, taxon_filter=args.taxon_filter)

    # Call detection
    per_sample = call_detection(
        per_sample,
        min_reads=args.min_reads,
        threshold=args.threshold,
        min_sample_reads=args.min_sample_reads,
    )

    # Detection rates per group
    rates = detection_rates(per_sample, meta, args.group_by, args.group_order)

    # Fisher's exact + chi-square
    fisher = run_fisher(rates)

    # Print summary
    print_summary(
        rates, fisher,
        sheet=args.sheet,
        taxon_filter=args.taxon_filter,
        min_reads=args.min_reads,
        threshold=args.threshold,
    )

    if args.dry_run:
        log.info("DRY RUN — no files written.")
        return 0

    # Write outputs
    stem = args.output_stem or f"{args.sheet}_stats"
    args.outdir.mkdir(parents=True, exist_ok=True)

    rates_path  = args.outdir / f"{stem}_detection_rates.tsv"
    fisher_path = args.outdir / f"{stem}_fisher.tsv"

    rates.to_csv(rates_path, sep="\t", index=False)
    fisher.to_csv(fisher_path, sep="\t", index=False)

    log.info("Saved: %s", rates_path)
    log.info("Saved: %s", fisher_path)
    log.info("=== Done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
