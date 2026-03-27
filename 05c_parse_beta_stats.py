#!/usr/bin/env python3
"""
parse_beta_stats.py
───────────────────
Parse QIIME 2 beta-group-significance QZV files (PERMANOVA and PERMDISP)
and write a clean summary TSV.

Usage
-----
# Point at a directory containing .qzv files:
python parse_beta_stats.py --stats-dir results/16S/diversity/core_metrics_depth8000 \
                           --output beta_stats_summary.tsv

# Or point at the raw extracted QZV folders (if already unzipped):
python parse_beta_stats.py --stats-dir stats/ --output beta_stats_summary.tsv

# Dry-run: just print table to stdout without writing file:
python parse_beta_stats.py --stats-dir results/ --print-only

Output columns
--------------
metric | method | statistic_name | statistic | p_value | n_permutations | sample_size | group_column
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import zipfile
from pathlib import Path
from typing import Optional, List


# ── HTML parser ───────────────────────────────────────────────────────────────

def parse_qzv_html(qzv_path: Path) -> "Optional[dict]":
    """
    Open a QZV (ZIP archive), find index.html, and extract the overview table.
    Returns a dict with keys: method, statistic_name, statistic, p_value,
    n_permutations, sample_size.  Returns None if not parseable.
    """
    try:
        with zipfile.ZipFile(qzv_path, "r") as zf:
            html_name = next(
                (n for n in zf.namelist() if n.endswith("/data/index.html")),
                None,
            )
            if html_name is None:
                return None
            html = zf.read(html_name).decode("utf-8", errors="ignore")
    except Exception as e:
        print(f"  [warn] Could not open {qzv_path.name}: {e}", file=sys.stderr)
        return None

    return _parse_html(html)


def parse_extracted_html(html_path: Path) -> "Optional[dict]":
    """Parse an already-extracted index.html."""
    try:
        html = html_path.read_text(encoding="utf-8", errors="ignore")
    except Exception as e:
        print(f"  [warn] Could not read {html_path}: {e}", file=sys.stderr)
        return None
    return _parse_html(html)


def _parse_html(html: str) -> "Optional[dict]":
    """
    Extract key/value pairs from the QIIME2 beta-group-significance HTML table.
    Handles both PERMANOVA (pseudo-F) and PERMDISP (F-value).
    """
    rows = re.findall(
        r'<th[^>]*>(.*?)</th>\s*<td[^>]*>(.*?)</td>',
        html, re.DOTALL
    )
    data = {}
    for k, v in rows:
        key = re.sub(r'<[^>]+>', '', k).strip().lower().replace(' ', '_')
        val = re.sub(r'<[^>]+>', '', v).strip()
        data[key] = val

    if not data:
        return None

    return {
        "method":           data.get("method_name", ""),
        "statistic_name":   data.get("test_statistic_name", ""),
        "statistic":        data.get("test_statistic", ""),
        "p_value":          data.get("p-value", data.get("p_value", "")),
        "n_permutations":   data.get("number_of_permutations", "999"),
        "sample_size":      data.get("sample_size", ""),
    }


# ── File discovery ─────────────────────────────────────────────────────────────

# Recognise metric and group from typical QIIME2 output filenames:
#   weighted_unifrac_permanova_Group.qzv
#   unweighted_permanova/  (extracted folder name)
METRIC_PATTERNS = {
    "unweighted_unifrac": ["unweighted_unifrac", "unweighted"],
    "weighted_unifrac":   ["weighted_unifrac", "weighted"],
    "bray_curtis":        ["bray_curtis", "braycurtis"],
    "jaccard":            ["jaccard"],
}


def _infer_metric(name: str) -> str:
    """
    Infer the beta diversity metric from a QZV filename or folder name.

    Checks the name against METRIC_PATTERNS in order; returns the first match,
    or 'unknown' if none match. Order matters: 'weighted_unifrac' is checked
    before 'weighted' to avoid the shorter pattern matching first.

    Args:
        name: QZV stem or directory name (e.g. 'weighted_unifrac_permanova_Group').

    Returns:
        One of 'unweighted_unifrac', 'weighted_unifrac', 'bray_curtis',
        'jaccard', or 'unknown'.
    """
    name_l = name.lower()
    for metric, patterns in METRIC_PATTERNS.items():
        if any(p in name_l for p in patterns):
            return metric
    return "unknown"


def _infer_group(name: str) -> str:
    """
    Extract the metadata group column name from a QZV filename or folder name.

    Strips known metric and method tokens, then joins the remaining parts with
    underscores. For example:
        'weighted_unifrac_permanova_Group' -> 'Group'
        'bray_curtis_permdisp_DiseaseStatus' -> 'DiseaseStatus'

    Returns an empty string if all parts are known tokens.
    """
    # e.g. weighted_unifrac_permanova_Group.qzv -> Group
    parts = Path(name).stem.split("_")
    skip = {"weighted", "unifrac", "unweighted", "bray", "curtis", "jaccard",
            "permanova", "permdisp"}
    group_parts = [p for p in parts if p.lower() not in skip]
    return "_".join(group_parts) if group_parts else ""


def collect_results(stats_dir: Path) -> List[dict]:
    """
    Walk stats_dir looking for either:
      - *.qzv files (packed)
      - subdirectories containing index.html (already extracted by QIIME view)
    """
    results = []

    # Packed QZVs
    for qzv in sorted(stats_dir.rglob("*.qzv")):
        stem = qzv.stem
        if "permanova" not in stem.lower() and "permdisp" not in stem.lower():
            continue
        parsed = parse_qzv_html(qzv)
        if parsed:
            parsed["metric"]       = _infer_metric(stem)
            parsed["group_column"] = _infer_group(stem)
            parsed["source"]       = str(qzv.relative_to(stats_dir))
            results.append(parsed)

    # Extracted folders (index.html inside data/ or directly)
    for html in sorted(stats_dir.rglob("index.html")):
        parent = html.parent.name  # e.g. "unweighted_permanova"
        if "permanova" not in parent.lower() and "permdisp" not in parent.lower():
            continue
        # Skip if we already captured this via QZV
        already = any(r["source"].endswith(parent) for r in results)
        if already:
            continue
        parsed = parse_extracted_html(html)
        if parsed:
            parsed["metric"]       = _infer_metric(parent)
            parsed["group_column"] = _infer_group(parent)
            parsed["source"]       = str(html.relative_to(stats_dir))
            results.append(parsed)

    return results


# ── Formatting ────────────────────────────────────────────────────────────────

METRIC_LABELS = {
    "weighted_unifrac":   "Weighted UniFrac",
    "unweighted_unifrac": "Unweighted UniFrac",
    "bray_curtis":        "Bray-Curtis",
    "jaccard":            "Jaccard",
    "unknown":            "Unknown",
}

COLUMNS = [
    "metric", "group_column", "method", "statistic_name",
    "statistic", "p_value", "n_permutations", "sample_size", "source"
]


def _sig(p: str) -> str:
    """
    Convert a p-value string to a significance star notation.

    Returns '***' (p<=0.001), '**' (p<=0.01), '*' (p<=0.05), 'ns' (p>0.05),
    or '' if the string cannot be parsed as a float.
    """
    try:
        v = float(p)
        if v <= 0.001: return "***"
        if v <= 0.01:  return "**"
        if v <= 0.05:  return "*"
        return "ns"
    except ValueError:
        return ""


def print_table(results: list[dict]) -> None:
    """
    Print a formatted summary table of PERMANOVA/PERMDISP results to stdout.

    Results are sorted by metric (Weighted UniFrac -> Unweighted -> Bray-Curtis
    -> Jaccard) and then by method (PERMANOVA before PERMDISP). Each row shows
    the metric, group column, method, test statistic name and value, p-value,
    number of permutations, sample size, and significance stars.

    Args:
        results: List of result dicts from collect_results().
    """
    if not results:
        print("No PERMANOVA/PERMDISP results found.")
        return

    # Sort: metric order, then permanova before permdisp
    order = list(METRIC_LABELS.keys())
    results.sort(key=lambda r: (
        order.index(r["metric"]) if r["metric"] in order else 99,
        0 if r["method"] == "permanova" else 1
    ))

    header = (
        f"{'Metric':<22}  {'Group':<12}  {'Method':<10}  "
        f"{'Statistic':<10}  {'Value':>8}  {'p-value':>8}  "
        f"{'Permutations':>14}  {'n':>6}  {'Sig':>4}"
    )
    rule = "─" * len(header)
    print(rule)
    print(header)
    print(rule)
    for r in results:
        label = METRIC_LABELS.get(r["metric"], r["metric"])
        stat  = r["statistic_name"].replace("-value", "")
        try:
            val = f"{float(r['statistic']):.4f}"
        except ValueError:
            val = r["statistic"]
        try:
            pval = f"{float(r['p_value']):.3f}"
        except ValueError:
            pval = r["p_value"]
        sig = _sig(r["p_value"])
        print(
            f"{label:<22}  {r['group_column']:<12}  {r['method'].upper():<10}  "
            f"{stat:<10}  {val:>8}  {pval:>8}  "
            f"{r['n_permutations']:>14}  {r['sample_size']:>6}  {sig:>4}"
        )
    print(rule)
    print("Significance: *** p≤0.001  ** p≤0.01  * p≤0.05  ns = not significant")


def write_tsv(results: list[dict], output: Path) -> None:
    """
    Write PERMANOVA/PERMDISP results to a tab-separated file.

    Columns: metric, group_column, method, statistic_name, statistic,
    p_value, n_permutations, sample_size, source. Extra keys are ignored.
    The file can be passed to 07_visualize_diversity.py via --stats-tsv.

    Args:
        results: List of result dicts from collect_results().
        output:  Destination path for the TSV file.
    """
    with output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"\nWrote {len(results)} rows to {output}")


# ── CLI ────────────────────────────────────────────────────────────────────────

def main() -> None:
    """
    Parse command-line arguments, scan for PERMANOVA/PERMDISP results, print
    a formatted table to stdout, and optionally write a summary TSV.

    Exit codes: sys.exit(1) if --stats-dir does not exist; otherwise exits 0.
    """
    ap = argparse.ArgumentParser(
        description="Parse QIIME 2 beta-group-significance QZVs into a summary table."
    )
    ap.add_argument(
        "--stats-dir", required=True, type=Path,
        help="Directory to search for PERMANOVA/PERMDISP QZVs or extracted folders."
    )
    ap.add_argument(
        "--output", type=Path, default=None,
        help="Output TSV path. Defaults to <stats-dir>/beta_stats_summary.tsv"
    )
    ap.add_argument(
        "--print-only", action="store_true",
        help="Print table to stdout only, do not write file."
    )
    args = ap.parse_args()

    if not args.stats_dir.exists():
        print(f"Error: --stats-dir not found: {args.stats_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning: {args.stats_dir}")
    results = collect_results(args.stats_dir)
    print(f"Found {len(results)} result(s).\n")

    print_table(results)

    if not args.print_only:
        out = args.output or args.stats_dir / "beta_stats_summary.tsv"
        write_tsv(results, out)


if __name__ == "__main__":
    main()
