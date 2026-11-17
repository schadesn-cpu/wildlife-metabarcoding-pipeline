#!/usr/bin/env python3
"""
diagnose_substring_traps.py
===========================
Audit how 09_pathogen_host_cooccurrence.py's substring matching compares
to true genus-level (g__X) classifications in the SILVA taxonomy export.

For each pathogen target, reports:
  - ASVs matched by the current substring rule  (what 09 counts now)
  - ASVs whose genus is actually g__<target>    (what 09 should count)
  - The delta — false positives, with their actual genus and read counts
  - Sample-level impact:
      n_samples currently scored positive,
      n_samples that would still be positive after correction,
      n_samples whose positivity is driven *only* by false-positive ASVs

Read-only. Writes a single audit TSV; produces no qza / cooccurrence outputs.

Usage (from project root):
    # Most common — relies on conventional paths
    python scripts/diagnose_substring_traps.py

    # Override anything
    python scripts/diagnose_substring_traps.py \\
        --taxonomy qiime2/16S/all/taxonomy-exported-temp/taxonomy.tsv \\
        --feature-table /tmp/16s_export/feature-table.tsv \\
        --min-reads 5 \\
        --outdir results/16S/all/blast_verify

The targets and substring rules are read from 09_pathogen_host_cooccurrence.py
by default, so the audit always reflects the script's actual behavior. Pass
--targets-from-config to override with values in pipeline_config.yml.

Author: Samantha Schade - MEED Lab, UNH - 2026-05-01
"""

import argparse
import re
import sys
import traceback
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Defaults — match the dict in 09_pathogen_host_cooccurrence.py
# Kept here as a fallback; the script tries to import from 09 first.
# ---------------------------------------------------------------------------

DEFAULT_TARGETS_16S = {
    "Anaplasma":   ["Anaplasma"],
    "Borreliella": ["Borreliella"],
    "Borrelia":    ["Borrelia"],
    "Francisella": ["Francisella"],
    "Bartonella":  ["Bartonella"],
    "Rickettsia":  ["Rickettsia"],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

GENUS_RE = re.compile(r"g__([^;]*)")


def extract_genus(taxon: str) -> str:
    """Return the genus token from a SILVA-style Taxon string, or '' if absent."""
    if not isinstance(taxon, str):
        return ""
    m = GENUS_RE.search(taxon)
    return m.group(1) if m else ""


def substring_match(taxon: str, terms: list) -> bool:
    """Replicate 09's matching logic: case-insensitive substring on Taxon."""
    if not isinstance(taxon, str):
        return False
    lo = taxon.lower()
    return any(term.lower() in lo for term in terms)


def genus_match(taxon: str, genus_name: str) -> bool:
    """Strict genus-level match: g__<genus_name> (case-sensitive)."""
    return extract_genus(taxon) == genus_name


def load_taxonomy(path: Path) -> pd.DataFrame:
    """Load the exported SILVA taxonomy.tsv (Feature ID, Taxon, Confidence)."""
    if not path.exists():
        raise FileNotFoundError(f"taxonomy.tsv not found: {path}")
    df = pd.read_csv(path, sep="\t")
    # SILVA exports use 'Feature ID' (with space). Normalize.
    df.columns = [c.strip().replace(" ", "_") for c in df.columns]
    if "Feature_ID" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Feature_ID"})
    if "Taxon" not in df.columns:
        df = df.rename(columns={df.columns[1]: "Taxon"})
    df["genus"] = df["Taxon"].map(extract_genus)
    return df.set_index("Feature_ID")


def load_feature_table(path: Path) -> pd.DataFrame:
    """
    Load the BIOM-converted feature table TSV.
    Format: first row '# Constructed from biom file' (skip), header, then
    Feature ID + per-sample read counts.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"feature-table.tsv not found: {path}\n"
            "  Export it first with:\n"
            "    qiime tools export --input-path <table.qza> "
            "--output-path /tmp/16s_export\n"
            "    biom convert -i /tmp/16s_export/feature-table.biom "
            "-o /tmp/16s_export/feature-table.tsv --to-tsv"
        )
    return pd.read_csv(path, sep="\t", skiprows=1, index_col=0)


def load_targets_from_09(scripts_dir: Path) -> dict:
    """
    Try to read pathogens_16s out of 09_pathogen_host_cooccurrence.py so the
    audit reflects whatever's actually in the script. Falls back to defaults.
    """
    nine = scripts_dir / "09_pathogen_host_cooccurrence.py"
    if not nine.exists():
        return DEFAULT_TARGETS_16S
    src = nine.read_text()
    # crude but adequate — pull the literal dict assigned to pathogens_16s
    m = re.search(r"pathogens_16s\s*=\s*(\{[^}]*\})", src, re.DOTALL)
    if not m:
        return DEFAULT_TARGETS_16S
    try:
        # safe-ish eval: only literal dicts of strings/lists
        ns = {}
        exec("d = " + m.group(1), {"__builtins__": {}}, ns)
        d = ns["d"]
        if isinstance(d, dict) and all(isinstance(v, list) for v in d.values()):
            return d
    except Exception:
        pass
    return DEFAULT_TARGETS_16S


# ---------------------------------------------------------------------------
# The audit itself
# ---------------------------------------------------------------------------

def audit_target(target: str, terms: list,
                 tax: pd.DataFrame, table: pd.DataFrame,
                 min_reads: int) -> dict:
    """
    Return a dict with per-target counts and a per-ASV breakdown for the
    false-positive ASVs (the ones substring catches but genus-level wouldn't).
    """
    sub_mask = tax["Taxon"].map(lambda t: substring_match(t, terms))
    gen_mask = tax["Taxon"].map(lambda t: genus_match(t, target))

    sub_ids = tax.index[sub_mask]
    gen_ids = tax.index[gen_mask]
    fp_ids  = sub_ids.difference(gen_ids)   # caught by substring, not by genus
    fn_ids  = gen_ids.difference(sub_ids)   # rare but possible: genus exists but rule misses

    # restrict to ASVs that are actually in the feature table
    in_table = table.index
    sub_present = sub_ids.intersection(in_table)
    gen_present = gen_ids.intersection(in_table)
    fp_present  = fp_ids.intersection(in_table)

    # per-sample positivity using the same threshold 09 uses
    def n_samples_with(asvs) -> int:
        if len(asvs) == 0:
            return 0
        return int((table.loc[asvs] >= min_reads).any(axis=0).sum())

    n_sub  = n_samples_with(sub_present)
    n_gen  = n_samples_with(gen_present)
    # samples that go positive ONLY because of false positives — i.e.,
    # currently scored positive, but would not be after the fix
    fp_only_samples = 0
    if len(sub_present) > 0:
        any_sub = (table.loc[sub_present] >= min_reads).any(axis=0)
        any_gen = ((table.loc[gen_present] >= min_reads).any(axis=0)
                   if len(gen_present) > 0
                   else pd.Series(False, index=table.columns))
        fp_only_samples = int((any_sub & ~any_gen).sum())

    # per-ASV breakdown of the false positives
    fp_rows = []
    for fid in fp_present:
        reads_total = int(table.loc[fid].sum())
        n_samples_pos = int((table.loc[fid] >= min_reads).sum())
        fp_rows.append({
            "feature_id":          fid,
            "current_target":      target,
            "actual_genus":        tax.at[fid, "genus"] or "(none)",
            "actual_taxon":        tax.at[fid, "Taxon"],
            "confidence":          tax.at[fid, "Confidence"]
                                    if "Confidence" in tax.columns else "",
            "total_reads":         reads_total,
            f"n_samples>={min_reads}": n_samples_pos,
        })
    fp_df = (pd.DataFrame(fp_rows)
             .sort_values(["actual_genus", "total_reads"],
                          ascending=[True, False])
             if fp_rows else pd.DataFrame())

    return {
        "target":           target,
        "rule":             ", ".join(terms),
        "n_asvs_substring": len(sub_present),
        "n_asvs_genus":     len(gen_present),
        "n_asvs_fp":        len(fp_present),
        "n_asvs_fn":        len(fn_ids.intersection(in_table)),
        "n_samples_now":    n_sub,
        "n_samples_fixed":  n_gen,
        "n_samples_fp_only": fp_only_samples,
        "fp_breakdown":     fp_df,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(
        description="Audit substring-matching false positives in "
                    "09_pathogen_host_cooccurrence.py against true "
                    "genus-level classifications.",
    )
    p.add_argument("--taxonomy",
                   default="qiime2/16S/all/taxonomy-exported-temp/taxonomy.tsv",
                   help="Path to exported SILVA taxonomy.tsv")
    p.add_argument("--feature-table",
                   default="/tmp/16s_export/feature-table.tsv",
                   help="Path to BIOM-converted feature-table.tsv")
    p.add_argument("--scripts-dir", default="scripts",
                   help="Directory containing 09_pathogen_host_cooccurrence.py "
                        "(used to read the actual targets dict)")
    p.add_argument("--targets",
                   choices=["from-09", "defaults"], default="from-09",
                   help="Source of the pathogens_16s dict")
    p.add_argument("--min-reads", type=int, default=5,
                   help="Read threshold for sample positivity (matches 09's default)")
    p.add_argument("--outdir", default="results/16S/all/blast_verify",
                   help="Where to write the audit TSV(s)")
    args = p.parse_args()

    try:
        tax_path   = Path(args.taxonomy)
        table_path = Path(args.feature_table)
        outdir     = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"Taxonomy:      {tax_path}")
        print(f"Feature table: {table_path}")
        print(f"Outdir:        {outdir}")

        tax   = load_taxonomy(tax_path)
        table = load_feature_table(table_path)
        print(f"  Loaded {len(tax):>6} taxonomy rows, "
              f"{len(table):>6} ASVs in feature table, "
              f"{table.shape[1]:>4} samples")

        if args.targets == "from-09":
            targets = load_targets_from_09(Path(args.scripts_dir))
        else:
            targets = DEFAULT_TARGETS_16S
        print(f"\nAuditing {len(targets)} targets: {list(targets)}")

    except FileNotFoundError as e:
        print(f"\n[ERROR] {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"\n[UNEXPECTED ERROR loading inputs] "
              f"{type(e).__name__}: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return 1

    summary_rows = []
    n_failed = 0
    for tgt, terms in targets.items():
        print(f"\n>> {tgt}  (rule: {terms})")
        try:
            r = audit_target(tgt, terms, tax, table, args.min_reads)
        except Exception as e:
            print(f"  [FAILED] {type(e).__name__}: {e}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            n_failed += 1
            continue

        print(f"  ASVs substring-matched: {r['n_asvs_substring']:>5}")
        print(f"  ASVs genus-matched:     {r['n_asvs_genus']:>5}")
        print(f"  False positives:        {r['n_asvs_fp']:>5}  "
              f"(missed by genus rule: {r['n_asvs_fn']})")
        print(f"  Samples positive now:   {r['n_samples_now']:>5}")
        print(f"  Samples positive fixed: {r['n_samples_fixed']:>5}")
        print(f"  Samples lost on fix:    {r['n_samples_fp_only']:>5}  "
              f"(ONLY positive due to FP ASVs)")

        if not r["fp_breakdown"].empty:
            fp_path = outdir / f"audit_fp_{tgt}.tsv"
            r["fp_breakdown"].to_csv(fp_path, sep="\t", index=False)
            print(f"  Wrote per-ASV breakdown -> {fp_path.name}")

        summary_rows.append({k: v for k, v in r.items() if k != "fp_breakdown"})

    summary_path = outdir / "audit_substring_traps_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    print(f"\nSummary -> {summary_path}")
    print(f"DONE  |  targets: {len(targets)}  |  failed: {n_failed}")
    return 1 if n_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
