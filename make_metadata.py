#!/usr/bin/env python3
"""
make_metadata.py
================
Step: metadata

Purpose:
    Build a QIIME 2-ready sample metadata TSV whose `sample-id` column matches
    the sample ids actually present in a feature table exactly. Bridges the long
    table ids (e.g. TV230084-GI-16S_S1492) to a source sample sheet keyed on the
    short biological id (e.g. TV230084), using the config's id_regex.

Inputs:
    --table             feature table .qza (or exported feature-table.biom)
    --source-metadata   source-of-truth sample sheet (TSV or CSV)
    --source-id-column  the column in the sheet holding the short ids
    pipeline_config.yml samples.id_regex, samples.control_prefixes, metadata paths

Outputs:
    metadata/qiime/metadata_<marker>.tsv   (or --out)
    a mapping report TSV when any table samples are missing from the sheet
    logs/run_manifest.jsonl                (run appended on success)

Unmatched handling (no silent drops):
    - a non-control id that does not match id_regex stops the run (format/regex
      problem) — see utils/samples.py
    - a non-control id that matches but is absent from the source sheet is
      reported loudly and written to the mapping report; the row is kept blank
      so you can see exactly which samples your sheet is missing

Usage:
    python make_metadata.py --marker 16S \\
        --table qiime2/16S/dada2/table_16S.qza \\
        --source-metadata metadata/source_metadata.tsv \\
        --source-id-column TV
    # --out defaults to the metadata path configured for the marker

Requirements:
    Python >= 3.8, pandas; biom-format and QIIME 2 only when --table is a .qza.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# --- make config_loader and the utils package importable regardless of cwd ---
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from config_loader import load_config, get_metadata_path          # noqa: E402
from utils import validate, checkpoint, provenance, samples       # noqa: E402

log = logging.getLogger("make_metadata")


# ===========================================================================
# Feature-table sample-id reading  (logic preserved; QIIME/biom at runtime)
# ===========================================================================

def export_qza_to_biom(table_qza: Path, tmpdir: Path) -> Path:
    """Export a FeatureTable[Frequency] .qza and return the feature-table.biom."""
    outdir = tmpdir / "exported_table"
    outdir.mkdir(parents=True, exist_ok=True)
    log.info("exporting %s", table_qza)
    subprocess.run(
        ["qiime", "tools", "export",
         "--input-path", str(table_qza), "--output-path", str(outdir)],
        check=True,
    )
    biom_fp = outdir / "feature-table.biom"
    if not biom_fp.exists():
        raise FileNotFoundError(
            f"Expected BIOM at {biom_fp} after exporting {table_qza}, but it was not found."
        )
    return biom_fp


def read_sample_ids_from_biom(biom_fp: Path) -> List[str]:
    """Read the sample ids from a BIOM feature table."""
    try:
        from biom import load_table  # type: ignore
    except ImportError as ex:
        raise RuntimeError(
            "biom-format is required to read BIOM files. Install with:\n"
            "  conda install -c conda-forge biom-format"
        ) from ex
    return list(load_table(str(biom_fp)).ids(axis="sample"))


def detect_delimiter(fp: Path) -> str:
    """Prefer tab, fall back to comma, for the source sample sheet."""
    head = "\n".join(fp.read_text(encoding="utf-8", errors="replace").splitlines()[:3])
    if "\t" in head:
        return "\t"
    if "," in head:
        return ","
    return "\t"


# ===========================================================================
# Row building
# ===========================================================================

def build_rows(
    table_sample_ids: List[str],
    source_df: pd.DataFrame,
    source_id_col: str,
    keys: Dict[str, Optional[str]],
    marker: Optional[str],
    control_sampletype_value: str,
    sampletype_column: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[Tuple[str, str]]]:
    """
    Build the output metadata rows and a mapping report.

    `keys` maps each table id to its extracted biological key (None for controls)
    and is produced by utils.samples.extract_sample_ids, which has already failed
    loud on any non-control id that didn't match the regex.

    Returns (metadata_df, mapping_df, missing_from_source) where missing_from_source
    lists (table_id, key) for non-control samples whose key is absent from the sheet.
    """
    if source_id_col not in source_df.columns:
        raise validate.ValidationError(
            f"Source metadata has no column {source_id_col!r}.\n"
            f"  Available columns: {list(source_df.columns)}"
        )

    source_df = source_df.copy()
    source_df[source_id_col] = source_df[source_id_col].astype(str)

    dupes = source_df[source_id_col][source_df[source_id_col].duplicated()].unique().tolist()
    if dupes:
        log.warning("duplicate ids in source column %r (keeping first): %s%s",
                    source_id_col, dupes[:10], " (truncated)" if len(dupes) > 10 else "")
        source_df = source_df.drop_duplicates(subset=[source_id_col], keep="first")

    source_df = source_df.set_index(source_id_col)
    data_columns = [c for c in source_df.reset_index().columns if c != source_id_col]

    rows: List[dict] = []
    map_rows: List[dict] = []
    missing_from_source: List[Tuple[str, str]] = []

    for sid in table_sample_ids:
        key = keys.get(sid)
        is_control = key is None
        matched = False

        if key is not None and key in source_df.index:
            src = source_df.loc[key]
            r = src.to_dict() if hasattr(src, "to_dict") else dict(src)
            matched = True
        else:
            r = {c: "" for c in data_columns}
            if is_control:
                r[sampletype_column] = control_sampletype_value
            else:
                missing_from_source.append((sid, key))

        r["sample-id"] = sid
        if marker:
            r["Marker"] = marker
        rows.append(r)
        map_rows.append({
            "table_sample_id": sid,
            "extracted_key": key if key is not None else "",
            "is_control": str(is_control),
            "matched": str(matched),
        })

    out_df = pd.DataFrame(rows)
    out_df = out_df[["sample-id"] + [c for c in out_df.columns if c != "sample-id"]]
    return out_df, pd.DataFrame(map_rows), missing_from_source


# ===========================================================================
# Orchestration
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="make_metadata.py",
        description="Build QIIME 2 metadata matching feature-table sample ids.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--table", required=True,
                   help="Feature table .qza or exported feature-table.biom.")
    p.add_argument("--source-metadata", required=True,
                   help="Source sample sheet (TSV or CSV).")
    p.add_argument("--source-id-column", required=True,
                   help="Column in the sheet holding the short join ids.")
    p.add_argument("--marker", default=None,
                   help="Marker label; sets the default --out and adds a Marker column.")
    p.add_argument("--out", default=None,
                   help="Output metadata TSV (default: configured path for --marker).")
    p.add_argument("--mapping-report", default=None,
                   help="Where to write the mapping report (default: alongside --out).")
    p.add_argument("--key-regex", default=None,
                   help="Override samples.id_regex (one capture group).")
    p.add_argument("--sampletype-column", default="SampleType")
    p.add_argument("--control-sampletype-value", default="Control")
    p.add_argument("--config", default=None, help="Path to pipeline_config.yml.")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s  %(message)s",
    )

    cfg = load_config(args.config)
    table_path = Path(args.table)
    source_path = Path(args.source_metadata)

    # --- resolve output path: explicit --out, else the configured metadata path
    if args.out:
        out_path = Path(args.out)
    elif args.marker:
        out_path = get_metadata_path(cfg, args.marker, "all")
    else:
        log.error("Provide --out, or --marker so the output path can be read from config.")
        return 2

    # --- validate inputs on entry ------------------------------------------
    if not source_path.exists():
        log.error("Source metadata not found: %s", source_path)
        return 2
    if not table_path.exists():
        log.error("Feature table not found: %s", table_path)
        return 2
    if table_path.suffix == ".qza":
        try:
            validate.require_qiime()
        except validate.ValidationError as e:
            log.error("%s", e)
            return 2

    # --- read the source sheet and the table's sample ids ------------------
    source_df = pd.read_csv(source_path, sep=detect_delimiter(source_path), dtype=str)
    if table_path.suffix == ".qza":
        with tempfile.TemporaryDirectory() as td:
            table_sample_ids = read_sample_ids_from_biom(
                export_qza_to_biom(table_path, Path(td)))
    else:
        table_sample_ids = read_sample_ids_from_biom(table_path)

    if not table_sample_ids:
        log.error("No sample ids found in %s — is it a valid feature table?", table_path)
        return 2

    # --- extract keys (fails loud on a non-control id that doesn't match) ---
    try:
        keys = samples.extract_sample_ids(table_sample_ids, cfg, regex=args.key_regex)
    except samples.SampleIDError as e:
        log.error("%s", e)
        return 3

    out_df, map_df, missing = build_rows(
        table_sample_ids, source_df, args.source_id_column, keys,
        args.marker, args.control_sampletype_value, args.sampletype_column,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)

    # --- validate what we wrote --------------------------------------------
    try:
        validate.validate_metadata(out_path)
    except validate.ValidationError as e:
        log.error("metadata failed validation after writing:\n%s", e)
        return 3

    n_total = len(table_sample_ids)
    n_control = sum(1 for v in keys.values() if v is None)
    n_matched = int((map_df["matched"] == "True").sum())
    log.info("%d table samples: %d matched, %d controls, %d missing from sheet",
             n_total, n_matched, n_control, len(missing))

    # --- surface samples present in the table but absent from the sheet -----
    if missing:
        report = Path(args.mapping_report) if args.mapping_report else \
            out_path.with_name(out_path.stem + "_mapping_report.tsv")
        report.parent.mkdir(parents=True, exist_ok=True)
        map_df.to_csv(report, sep="\t", index=False)
        log.warning(
            "%d non-control sample(s) are in the table but missing from the "
            "source sheet — they have BLANK metadata. See %s\n  e.g. %s",
            len(missing), report,
            ", ".join(f"{sid}->{key}" for sid, key in missing[:8]),
        )

    checkpoint.print_checkpoint(
        cfg,
        "metadata",
        marker=args.marker,
        produced=[out_path],
        provenance={
            "marker": args.marker,
            "inputs": {"table": table_path, "source_metadata": source_path},
            "outputs": [out_path],
            "command": "python " + " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
            "extra": {"matched": n_matched, "controls": n_control,
                      "missing_from_sheet": len(missing)},
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
