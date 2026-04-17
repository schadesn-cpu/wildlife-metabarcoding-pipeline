#!/usr/bin/env python3
"""
run_all_figures.py
==================
Generates two complete figure sets for all markers.

  SET 1 — results/figures_annotated/  : Full labels, PERMANOVA F+p stats,
           ecological season ordering. For lab use and verification.

  SET 2 — results/figures_manuscript/ : Clean, no titles. For submission.

All paths and group names are read from pipeline_config.yml via config_loader.
No hardcoded project-specific values remain in this script.

Usage (run from project root):
    python scripts/run_all_figures.py [--dry-run]
    python scripts/run_all_figures.py --markers 16S MiFish
    python scripts/run_all_figures.py --analyses DvT COD
    python scripts/run_all_figures.py --config /path/to/pipeline_config.yml

Author: Samantha Schade - MEED Lab, UNH - 2026-04-03
"""

import argparse
import subprocess
import sys
import traceback
from pathlib import Path

# ---------------------------------------------------------------------------
# Config loader — find config_loader.py in scripts/ or parent directory
# ---------------------------------------------------------------------------

_PIPELINE_DIR = Path(__file__).resolve().parent
for _search_dir in [_PIPELINE_DIR, _PIPELINE_DIR.parent]:
    if (_search_dir / "config_loader.py").exists():
        if str(_search_dir) not in sys.path:
            sys.path.insert(0, str(_search_dir))
        break

import config_loader as _cl


def _find_plot_script() -> Path:
    """Locate 09_plot_diversity.py in scripts/ or the same directory."""
    for d in [_PIPELINE_DIR, _PIPELINE_DIR.parent]:
        p = d / "09_plot_diversity.py"
        if p.exists():
            return p
    raise FileNotFoundError(
        "09_plot_diversity.py not found.\n"
        f"  Searched: {_PIPELINE_DIR}  and  {_PIPELINE_DIR.parent}"
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _pcoa_arts(cfg, marker: str, core: Path) -> list:
    """Return PCoA QZA paths for a marker based on its config pcoa_metrics list."""
    metrics = cfg.markers.get(marker, {}).get(
        "pcoa_metrics", ["bray_curtis", "jaccard"])
    return [str(core / f"{m}_pcoa_results.qza") for m in metrics]


def _alpha_vecs(cfg, marker: str, core: Path) -> list:
    """Return alpha vector QZA paths for a marker based on its config alpha_metrics list."""
    metrics = cfg.markers.get(marker, {}).get(
        "alpha_metrics", ["observed_features", "shannon", "evenness"])
    return [str(core / f"{m}_vector.qza") for m in metrics]


def missing(files: list, label: str) -> bool:
    """Return True and print a warning if any file in the list does not exist."""
    m = [Path(f).name for f in files if not Path(f).exists()]
    if m:
        print(f"  Skip {label} - missing: {m}")
        return True
    return False


def no_meta(path, label: str) -> bool:
    """Return True and print a warning if the metadata file does not exist."""
    if not Path(path).exists():
        print(f"  Skip {label} - no metadata: {Path(path).name}")
        return True
    return False


def run(args: list, dry_run: bool, label: str) -> None:
    """
    Run 09_plot_diversity.py as a subprocess.

    Uses sys.executable so the same conda environment is used for the
    subprocess. Raises RuntimeError on non-zero exit so failures propagate
    to main() rather than being silently swallowed.
    """
    script = _find_plot_script()
    cmd = [sys.executable, str(script)] + args
    print(f"\n{'[DRY] ' if dry_run else ''}>> {label}")
    if dry_run:
        print(f"  {' '.join(str(x) for x in args[:6])} ...")
        return
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise RuntimeError(
            f"Figure generation failed for: {label}\n"
            f"  Exit code: {r.returncode}\n"
            f"  Check the output above for error details."
        )
    print(f"  Done")


def both(base: list, stem: str, stats_dir: Path,
         out_ann: Path, out_ms: Path, dry_run: bool, label: str) -> None:
    """Generate both the annotated and manuscript versions of a figure."""
    ann_extras = (["--stats-dir", str(stats_dir)] if base[0] == "pcoa"
                  else ["--stats-qzv-dir", str(stats_dir)] if base[0] == "alpha"
                  else [])
    run(
        base + ann_extras
             + ["--output-stem", f"{stem}_annotated",
                "--output-dir", str(out_ann)],
        dry_run, f"{label} - annotated",
    )
    run(
        base + ["--no-title",
                "--output-stem", stem,
                "--output-dir", str(out_ms)],
        dry_run, f"{label} - manuscript",
    )

def make(marker: str, analysis: str, ann: Path, ms: Path):
    """Create output directories and return (annotated_dir, manuscript_dir)."""
    a = ann / marker / analysis
    m = ms  / marker / analysis
    a.mkdir(parents=True, exist_ok=True)
    m.mkdir(parents=True, exist_ok=True)
    return a, m


# ---------------------------------------------------------------------------
# Per-analysis figure generators
# All paths and group names come from config - no hardcodes
# ---------------------------------------------------------------------------

def dvt_pcoa(cfg, ann: Path, ms: Path, dr: bool, marker: str) -> None:
    """Generate DvT PCoA panel figures (annotated + manuscript)."""
    try:
        core = _cl.get_diversity_dir(cfg, marker, "DvT")
    except ValueError as e:
        print(f"  Skip {marker} DvT PCoA: {e}")
        return
    try:
        meta = _cl.get_metadata_path(cfg, marker, "dvt")
    except (ValueError, KeyError):
        meta = _cl.get_metadata_path(cfg, marker, "all")        
    arts      = _pcoa_arts(cfg, marker, core)
    stats_dir = cfg.resolve(f"results/{marker}/DvT/diversity")
    group_col = cfg.groups.get("primary", {}).get("column", "Group")
    dvt_order = _cl.get_dvt_order(cfg)
    palette   = cfg.figures.get("palette", "wong")
    oa, om    = make(marker, "DvT", ann, ms)

    if missing(arts, f"{marker} DvT PCoA") or no_meta(meta, f"{marker} DvT PCoA"):
        return
    base = (["pcoa", "--artifact"] + arts
            + ["--metadata", str(meta),
               "--color-by", group_col,
               "--group-order"] + dvt_order
            + ["--panel", "--palette", palette])
    both(base, f"{marker}_DvT_pcoa_{palette}", stats_dir, oa, om, dr,
         f"{marker} DvT PCoA")


def dvt_alpha(cfg, ann: Path, ms: Path, dr: bool, marker: str) -> None:
    """Generate DvT alpha diversity panel figures (annotated + manuscript)."""
    try:
        core = _cl.get_diversity_dir(cfg, marker, "DvT")
    except ValueError as e:
        print(f"  Skip {marker} DvT PCoA: {e}")
        return
    try:
        meta = _cl.get_metadata_path(cfg, marker, "dvt")
    except (ValueError, KeyError):
        meta = _cl.get_metadata_path(cfg, marker, "all")
    vecs      = _alpha_vecs(cfg, marker, core)
    stats_dir = cfg.resolve(f"results/{marker}/DvT/diversity")
    group_col = cfg.groups.get("primary", {}).get("column", "Group")
    dvt_order = _cl.get_dvt_order(cfg)
    palette   = cfg.figures.get("palette", "wong")
    oa, om    = make(marker, "DvT", ann, ms)

    if missing(vecs, f"{marker} DvT alpha") or no_meta(meta, f"{marker} DvT alpha"):
        return
    base = (["alpha", "--artifact"] + vecs
            + ["--metadata", str(meta),
               "--group-by", group_col,
               "--group-order"] + dvt_order
            + ["--panel", "--palette", palette])
    both(base, f"{marker}_DvT_alpha_{palette}", stats_dir, oa, om, dr,
         f"{marker} DvT alpha")


def cod_pcoa(cfg, ann: Path, ms: Path, dr: bool, marker: str) -> None:
    """Generate COD PCoA panel figures (annotated + manuscript)."""
    try:
        core = _cl.get_diversity_dir(cfg, marker, "DvT")
        meta = _cl.get_metadata_path(cfg, marker, "cod")
    except ValueError as e:
        print(f"  Skip {marker} COD PCoA: {e}")
        return
    arts      = _pcoa_arts(cfg, marker, core)
    stats_dir = cfg.resolve(f"results/{marker}/COD/diversity")
    group_col = cfg.groups.get("secondary", {}).get("column", "COD_broad")
    cod_order = _cl.get_group_order(cfg, "secondary")
    palette   = cfg.figures.get("palette", "wong")
    oa, om    = make(marker, "COD", ann, ms)

    if missing(arts, f"{marker} COD PCoA") or no_meta(meta, f"{marker} COD PCoA"):
        return
    base = (["pcoa", "--artifact"] + arts
            + ["--metadata", str(meta),
               "--color-by", group_col]
            + (["--group-order"] + cod_order if cod_order else [])
            + ["--panel", "--palette", palette])
    both(base, f"{marker}_COD_pcoa_{palette}_filtered", stats_dir, oa, om, dr,
         f"{marker} COD PCoA")


def cod_alpha(cfg, ann: Path, ms: Path, dr: bool, marker: str) -> None:
    """Generate COD alpha diversity panel figures (annotated + manuscript)."""
    try:
        core = _cl.get_diversity_dir(cfg, marker, "DvT")
        meta = _cl.get_metadata_path(cfg, marker, "cod")
    except ValueError as e:
        print(f"  Skip {marker} COD alpha: {e}")
        return
    vecs      = _alpha_vecs(cfg, marker, core)
    stats_dir = cfg.resolve(f"results/{marker}/COD/diversity")
    group_col = cfg.groups.get("secondary", {}).get("column", "COD_broad")
    cod_order = _cl.get_group_order(cfg, "secondary")
    palette   = cfg.figures.get("palette", "wong")
    oa, om    = make(marker, "COD", ann, ms)

    if missing(vecs, f"{marker} COD alpha") or no_meta(meta, f"{marker} COD alpha"):
        return
    base = (["alpha", "--artifact"] + vecs
            + ["--metadata", str(meta),
               "--group-by", group_col]
            + (["--group-order"] + cod_order if cod_order else [])
            + ["--panel", "--palette", palette])
    both(base, f"{marker}_COD_alpha_{palette}_filtered", stats_dir, oa, om, dr,
         f"{marker} COD alpha")


def season_pcoa(cfg, ann: Path, ms: Path, dr: bool, marker: str) -> None:
    """Generate seasonal PCoA panel figures (annotated + manuscript)."""
    try:
        core = _cl.get_diversity_dir(cfg, marker, "season")
        meta = _cl.get_metadata_path(cfg, marker, "all")
    except ValueError as e:
        print(f"  Skip {marker} season PCoA: {e}")
        return
    arts         = _pcoa_arts(cfg, marker, core)
    stats_dir    = cfg.resolve(f"results/{marker}/season/diversity")
    group_col    = cfg.groups.get("seasonal", {}).get("column", "Season")
    season_order = _cl.get_group_order(cfg, "seasonal")
    palette      = cfg.figures.get("palette", "wong")
    oa, om       = make(marker, "season", ann, ms)

    if missing(arts, f"{marker} season PCoA") or no_meta(meta, f"{marker} season PCoA"):
        return
    base = (["pcoa", "--artifact"] + arts
            + ["--metadata", str(meta),
               "--color-by", group_col]
            + (["--group-order"] + season_order if season_order else [])
            + ["--panel", "--palette", palette])
    both(base, f"{marker}_season_pcoa_{palette}", stats_dir, oa, om, dr,
         f"{marker} season PCoA")


def season_alpha(cfg, ann: Path, ms: Path, dr: bool, marker: str) -> None:
    """Generate seasonal alpha diversity panel figures (annotated + manuscript)."""
    try:
        core = _cl.get_diversity_dir(cfg, marker, "season")
        meta = _cl.get_metadata_path(cfg, marker, "all")
    except ValueError as e:
        print(f"  Skip {marker} season alpha: {e}")
        return
    vecs         = _alpha_vecs(cfg, marker, core)
    stats_dir    = cfg.resolve(f"results/{marker}/season/diversity")
    group_col    = cfg.groups.get("seasonal", {}).get("column", "Season")
    season_order = _cl.get_group_order(cfg, "seasonal")
    palette      = cfg.figures.get("palette", "wong")
    oa, om       = make(marker, "season", ann, ms)

    if missing(vecs, f"{marker} season alpha") or no_meta(meta, f"{marker} season alpha"):
        return
    base = (["alpha", "--artifact"] + vecs
            + ["--metadata", str(meta),
               "--group-by", group_col]
            + (["--group-order"] + season_order if season_order else [])
            + ["--panel", "--palette", palette])
    both(base, f"{marker}_season_alpha_{palette}", stats_dir, oa, om, dr,
         f"{marker} season alpha")

# ---------------------------------------------------------------------------
# Taxonomy figure generators (unrarefied RRA — all groupings)
# ---------------------------------------------------------------------------

_TAXONOMY_RELABUND = {
    "16S":    "results/16S/all/taxonomy_refined/taxonomy_relabund_L6_16S.tsv",
    "MiFish": "results/MiFish/all/taxonomy_cleaned/taxonomy_relabund_L7_MiFish_cleaned.tsv",
    "cytb":   "results/cytb/all/taxonomy_cleaned/taxonomy_relabund_L7_cytb_cleaned.tsv",
    "18S":    "results/18S/all/taxonomy/taxonomy_relabund_L6_18S.tsv",
}


def _find_taxonomy_script() -> Path:
    for d in [_PIPELINE_DIR, _PIPELINE_DIR.parent]:
        p = d / "10_plot_taxonomy.py"
        if p.exists():
            return p
    raise FileNotFoundError("10_plot_taxonomy.py not found.")


def run_tax(args: list, dry_run: bool, label: str) -> None:
    script = _find_taxonomy_script()
    cmd = [sys.executable, str(script)] + args
    print(f"\n{'[DRY] ' if dry_run else ''}>> {label}")
    if dry_run:
        print(f"  {' '.join(str(x) for x in args[:6])} ...")
        return
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise RuntimeError(
            f"Taxonomy figure failed: {label}\n"
            f"  Exit code: {r.returncode}"
        )
    print(f"  Done")


def tax_both(base: list, stem: str,
             out_ann: Path, out_ms: Path,
             dry_run: bool, label: str) -> None:
    run_tax(base + ["--output-stem", f"{stem}_annotated",
                    "--outdir", str(out_ann)],
            dry_run, f"{label} - annotated")
    run_tax(base + ["--no-title",
                    "--output-stem", stem,
                    "--outdir", str(out_ms)],
            dry_run, f"{label} - manuscript")


def _tax_relabund(cfg, marker: str):
    rel = _TAXONOMY_RELABUND.get(marker)
    if rel is None:
        print(f"  Skip {marker} taxonomy — no relabund path configured")
        return None
    p = Path(cfg.resolve(rel))
    if not p.exists():
        print(f"  Skip {marker} taxonomy — relabund not found: {p.name}")
        return None
    return p


def tax_all(cfg, ann: Path, ms: Path, dr: bool, marker: str) -> None:
    """All samples grouped by Group (Diseased / Trauma / Marine)."""
    relabund = _tax_relabund(cfg, marker)
    if relabund is None:
        return
    try:
        meta = _cl.get_metadata_path(cfg, marker, "all")
    except ValueError as e:
        print(f"  Skip {marker} tax_all: {e}")
        return
    if no_meta(meta, f"{marker} tax_all"):
        return
    palette = cfg.figures.get("palette", "wong")
    oa, om = make(marker, "all", ann, ms)
    stem = f"{marker}_unrarefied_all_barplot_{palette}"
    base = ["--relabund", str(relabund),
            "--metadata", str(meta),
            "--group-by", "Group",
            "--group-order", "Diseased", "Trauma", "Marine",
            "--marker", marker,
            "--palette", palette]
    tax_both(base, stem, oa, om, dr, f"{marker} taxonomy all")


def tax_dvt(cfg, ann: Path, ms: Path, dr: bool, marker: str) -> None:
    """Diseased vs Trauma (Marine excluded via group-order)."""
    relabund = _tax_relabund(cfg, marker)
    if relabund is None:
        return
    try:
        meta = _cl.get_metadata_path(cfg, marker, "all")
    except ValueError as e:
        print(f"  Skip {marker} tax_dvt: {e}")
        return
    if no_meta(meta, f"{marker} tax_dvt"):
        return
    dvt_order = _cl.get_dvt_order(cfg)
    palette = cfg.figures.get("palette", "wong")
    oa, om = make(marker, "DvT", ann, ms)
    stem = f"{marker}_unrarefied_DvT_barplot_{palette}"
    base = (["--relabund", str(relabund),
             "--metadata", str(meta),
             "--group-by", "Group",
             "--group-order"] + dvt_order +
            ["--marker", marker,
             "--palette", palette])
    tax_both(base, stem, oa, om, dr, f"{marker} taxonomy DvT")


def tax_cod(cfg, ann: Path, ms: Path, dr: bool, marker: str) -> None:
    """COD subgroups: Lead / Parasitic_Infectious / Trauma."""
    relabund = _tax_relabund(cfg, marker)
    if relabund is None:
        return
    try:
        meta = _cl.get_metadata_path(cfg, marker, "cod")
    except ValueError as e:
        print(f"  Skip {marker} tax_cod: {e}")
        return
    if no_meta(meta, f"{marker} tax_cod"):
        return
    cod_order = _cl.get_group_order(cfg, "secondary") or \
                ["Lead", "Parasitic_Infectious", "Trauma"]
    palette = cfg.figures.get("palette", "wong")
    oa, om = make(marker, "COD", ann, ms)
    stem = f"{marker}_unrarefied_COD_barplot_{palette}"
    base = (["--relabund", str(relabund),
             "--metadata", str(meta),
             "--group-by", "COD_broad",
             "--group-order"] + cod_order +
            ["--marker", marker,
             "--palette", palette])
    tax_both(base, stem, oa, om, dr, f"{marker} taxonomy COD")


def tax_season(cfg, ann: Path, ms: Path, dr: bool, marker: str) -> None:
    """Ecological season grouping (Breeding / Freshwater_Nonbreeding / Saltwater)."""
    relabund = _tax_relabund(cfg, marker)
    if relabund is None:
        return
    try:
        meta = _cl.get_metadata_path(cfg, marker, "all")
    except ValueError as e:
        print(f"  Skip {marker} tax_season: {e}")
        return
    if no_meta(meta, f"{marker} tax_season"):
        return
    season_order = _cl.get_group_order(cfg, "seasonal") or \
                   ["Breeding", "Freshwater_Nonbreeding", "Saltwater"]
    palette = cfg.figures.get("palette", "wong")
    oa, om = make(marker, "season", ann, ms)
    stem = f"{marker}_unrarefied_season_barplot_{palette}"
    base = (["--relabund", str(relabund),
             "--metadata", str(meta),
             "--group-by", "Season",
             "--group-order"] + season_order +
            ["--marker", marker,
             "--palette", palette])
    tax_both(base, stem, oa, om, dr, f"{marker} taxonomy season")


# ---------------------------------------------------------------------------
# Analysis registry
# ---------------------------------------------------------------------------

GENS = {
    "DvT":    [dvt_pcoa,    dvt_alpha, tax_dvt],
    "COD":    [cod_pcoa,    cod_alpha, tax_cod],
    "season": [season_pcoa, season_alpha, tax_season],
    "all":    [tax_all],
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    """
    Entry point for run_all_figures.py.

    Reads all paths, group names, and markers from pipeline_config.yml
    via config_loader. Generates both annotated (lab) and manuscript
    (clean) figure sets for each marker and analysis combination.
    Returns 0 on success, 1 if any figures failed.
    """
    p = argparse.ArgumentParser(
        description="Generate annotated + manuscript figure sets from pipeline_config.yml"
    )
    p.add_argument("--config",   default=None, metavar="YML",
                   help="Path to pipeline_config.yml. Default: auto-discovered.")
    p.add_argument("--dry-run",  action="store_true")
    p.add_argument("--markers",  nargs="+", default=None,
                   help="Limit to specific markers. Default: active_markers in config.")
    p.add_argument("--analyses", nargs="+", default=["DvT", "COD", "season", "all"],
                    choices=["DvT", "COD", "season", "all"])
    args = p.parse_args()

    try:
        cfg = _cl.load_config(args.config)
    except (FileNotFoundError, ValueError) as e:
        print(f"\n[ERROR] {e}", file=sys.stderr)
        return 1

    markers = args.markers or cfg.active_markers
    ann     = cfg.resolve(cfg.figures.get("annotated_dir", "results/figures_annotated"))
    ms      = cfg.resolve(cfg.figures.get("manuscript_dir", "results/figures_manuscript"))

    print(f"Root:       {cfg.root}")
    print(f"Annotated:  {ann}")
    print(f"Manuscript: {ms}")
    print(f"Markers:    {markers}  |  Analyses: {args.analyses}  |  Dry: {args.dry_run}")

    try:
        _find_plot_script()
    except FileNotFoundError as e:
        print(f"\n[ERROR] {e}", file=sys.stderr)
        return 1

    n_total  = 0
    n_failed = 0

    for marker in markers:
        print(f"\n{'='*55}\nMARKER: {marker}\n{'='*55}")
        for analysis in args.analyses:
            for fn in GENS[analysis]:
                try:
                    fn(cfg, ann, ms, args.dry_run, marker)
                    n_total += 1
                except RuntimeError as exc:
                    print(f"\n  [FAILED] {exc}", file=sys.stderr)
                    n_failed += 1
                    n_total  += 1
                except Exception as exc:
                    print(f"\n  [UNEXPECTED ERROR] {type(exc).__name__}: {exc}",
                          file=sys.stderr)
                    traceback.print_exc(file=sys.stderr)
                    n_failed += 1
                    n_total  += 1

    print(f"\nDONE  |  annotated: {ann}  |  manuscript: {ms}")
    print(f"  Total: {n_total}  |  Failed: {n_failed}")

    if n_failed > 0:
        print(f"\n  {n_failed} figure(s) failed - check output above.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
