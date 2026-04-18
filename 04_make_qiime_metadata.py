#!/usr/bin/env python3
"""
04_make_qiime_metadata.py

Create a QIIME2-ready sample metadata TSV whose `sample-id` matches the sample IDs
present in a QIIME2 feature table (table.qza) or exported BIOM (feature-table.biom).

Typical problem solved:
- Feature table sample IDs look like:  TV230084-GI-16S_S1492
- Source metadata sample IDs look like: TV230084
QIIME2 requires exact matches, so we deterministically build a derived metadata file.

Outputs:
- A TSV with first column named `sample-id` that matches the table sample IDs exactly.
- Optional mapping report TSV with (table_sample_id -> extracted_key -> matched_source_id).

Requirements:
- python >= 3.8
- pandas
- biom-format  (only needed if input is BIOM or if exporting QZA produces BIOM)

Install (conda):
  conda install -c conda-forge pandas biom-format

Usage examples:
  python scripts/make_qiime_metadata.py \
    --table qiime2/dada2/table_16S.qza \
    --source-metadata metadata/source_metadata.tsv \
    --source-id-column TV \
    --out metadata/qiime/metadata_16S.tsv \
    --marker 16S
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd


def eprint(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def run_cmd(cmd: List[str]) -> None:
    """Log and execute a shell command, raising on non-zero exit."""
    eprint("[cmd]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def export_qza_to_biom(table_qza: Path, tmpdir: Path) -> Path:
    """
    Export QIIME2 FeatureTable[Frequency] .qza to a directory that includes feature-table.biom
    """
    outdir = tmpdir / "exported_table"
    outdir.mkdir(parents=True, exist_ok=True)
    run_cmd(["qiime", "tools", "export", "--input-path", str(table_qza), "--output-path", str(outdir)])
    biom_fp = outdir / "feature-table.biom"
    if not biom_fp.exists():
        raise FileNotFoundError(f"Expected BIOM at {biom_fp} after exporting {table_qza}, but it was not found.")
    return biom_fp


def read_sample_ids_from_biom(biom_fp: Path) -> List[str]:
    """
    Read sample IDs from BIOM file.
    """
    try:
        from biom import load_table  # type: ignore
    except Exception as ex:
        raise RuntimeError(
            "biom-format is required to read BIOM files. Install with:\n"
            "  conda install -c conda-forge biom-format\n"
        ) from ex

    table = load_table(str(biom_fp))
    return list(table.ids(axis="sample"))


def detect_delimiter(fp: Path) -> str:
    """
    Basic delimiter detection for source metadata.
    Prefers tab, falls back to comma.
    """
    head = fp.read_text(encoding="utf-8", errors="replace").splitlines()[:3]
    joined = "\n".join(head)
    if "\t" in joined:
        return "\t"
    if "," in joined:
        return ","
    # Default to tab; pandas will still read single-column files but we will validate later
    return "\t"


def normalize_column_names(cols: List[str]) -> List[str]:
    # Keep as-is; QIIME is picky only about output first column.
    """
Return column names unchanged.

    QIIME2 is strict only about the first column being named 'sample-id';
    all other column names are passed through as-is to preserve the source
    metadata schema.
    """
    return cols


def build_metadata(
    table_sample_ids: List[str],
    source_df: pd.DataFrame,
    source_id_col: str,
    out_fp: Path,
    marker: Optional[str],
    key_regex: str,
    control_prefixes: List[str],
    control_sampletype_value: str,
    sampletype_column: str,
    write_mapping_report: Optional[Path],
) -> Tuple[int, int, int]:
    """
    Returns: (n_table_samples, n_matched, n_unmatched)
    """
    if source_id_col not in source_df.columns:
        raise ValueError(
            f"Source metadata does not contain column '{source_id_col}'. "
            f"Available columns: {list(source_df.columns)}"
        )

    # make sure source ids are strings
    source_df = source_df.copy()
    source_df[source_id_col] = source_df[source_id_col].astype(str)

    # de-duplicate source IDs (warn if duplicates)
    dupes = source_df[source_id_col][source_df[source_id_col].duplicated()].unique().tolist()
    if dupes:
        eprint(f"[warn] Duplicate IDs found in source metadata column '{source_id_col}': {dupes[:10]} "
               f"{'(truncated)' if len(dupes) > 10 else ''}")
        # Keep first occurrence
        source_df = source_df.drop_duplicates(subset=[source_id_col], keep="first")

    source_df = source_df.set_index(source_id_col)

    key_re = re.compile(key_regex)

    rows = []
    map_rows = []

    for sid in table_sample_ids:
        # Determine if it's a control
        is_control = any(sid.startswith(pfx) for pfx in control_prefixes)

        # Extract join key (e.g., TV230084) from table sample ID
        m = key_re.search(sid)
        join_key = m.group(1) if m else ""

        matched = False
        src_row = None
        if join_key and join_key in source_df.index:
            src_row = source_df.loc[join_key]
            matched = True

        if matched:
            # src_row can be Series
            r = src_row.to_dict() if hasattr(src_row, "to_dict") else dict(src_row)
        else:
            # Fill with blanks (same columns as source)
            r = {c: "" for c in source_df.reset_index().columns if c != source_id_col}
            # Add/override SampleType for controls if requested
            if is_control:
                r[sampletype_column] = control_sampletype_value

        # Ensure sample-id matches table
        r["sample-id"] = sid
        if marker:
            r["Marker"] = marker

        rows.append(r)
        map_rows.append(
            {
                "table_sample_id": sid,
                "extracted_key": join_key,
                "matched_source_id": join_key if matched else "",
                "is_control": str(is_control),
                "matched": str(matched),
            }
        )

    out_df = pd.DataFrame(rows)

    # Put sample-id first
    cols = ["sample-id"] + [c for c in out_df.columns if c != "sample-id"]
    out_df = out_df[cols]

    out_fp.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_fp, sep="\t", index=False)

    if write_mapping_report:
        mr = pd.DataFrame(map_rows)
        write_mapping_report.parent.mkdir(parents=True, exist_ok=True)
        mr.to_csv(write_mapping_report, sep="\t", index=False)

    n_total = len(table_sample_ids)
    n_matched = sum(1 for r in map_rows if r["matched"] == "True")
    n_unmatched = n_total - n_matched
    return n_total, n_matched, n_unmatched


def parse_args() -> argparse.Namespace:
    """Build and parse the argument parser for 02_make_qiime_metadata.py."""
    p = argparse.ArgumentParser(
        description="Build QIIME2-ready metadata that matches feature table sample IDs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--table",
        required=True,
        help="QIIME2 feature table .qza OR exported BIOM feature-table.biom",
    )
    p.add_argument(
        "--source-metadata",
        required=True,
        help="Source-of-truth metadata file (TSV preferred; CSV ok).",
    )
    p.add_argument(
        "--source-id-column",
        required=True,
        help="Column in source metadata containing join keys (e.g., TV).",
    )
    p.add_argument(
        "--out",
        required=True,
        help="Output QIIME2 metadata TSV path.",
    )
    p.add_argument(
        "--marker",
        default=None,
        help="Optional marker label to add as a column (e.g., 16S).",
    )
    p.add_argument(
        "--key-regex",
        default=r"(TV\d+)",
        help="Regex with ONE capturing group used to extract join key from table sample-id.",
    )
    p.add_argument(
        "--control-prefixes",
        default="NTC-,PAC-,XB-",
        help="Comma-separated list of table sample-id prefixes considered controls.",
    )
    p.add_argument(
        "--sampletype-column",
        default="SampleType",
        help="Column name to use for control labeling (only filled for unmatched controls unless present in source).",
    )
    p.add_argument(
        "--control-sampletype-value",
        default="Control",
        help="Value assigned to sampletype-column for controls.",
    )
    p.add_argument(
        "--mapping-report",
        default=None,
        help="Optional path to write a mapping report TSV.",
    )
    return p.parse_args()


def main() -> int:
    """
Parse arguments, load the source metadata and feature table sample IDs,
    join them on the extracted key, and write a QIIME2-ready metadata TSV.

    Returns 0 on success, 2 if required input files are missing or the table
    contains no sample IDs. Prints a summary of matched vs unmatched samples
    to stderr. An optional mapping report TSV can be written for debugging.
    """
    args = parse_args()

    table_path = Path(args.table)
    source_meta_path = Path(args.source_metadata)
    out_path = Path(args.out)
    mapping_report_path = Path(args.mapping_report) if args.mapping_report else None

    if not source_meta_path.exists():
        eprint(f"[error] Source metadata not found: {source_meta_path}")
        return 2

    # Load source metadata
    delim = detect_delimiter(source_meta_path)
    source_df = pd.read_csv(source_meta_path, sep=delim, dtype=str)
    source_df.columns = normalize_column_names(list(source_df.columns))

    # Get table sample IDs
    biom_fp: Optional[Path] = None
    if table_path.suffix == ".qza":
        with tempfile.TemporaryDirectory() as td:
            tmpdir = Path(td)
            biom_fp = export_qza_to_biom(table_path, tmpdir)
            table_sample_ids = read_sample_ids_from_biom(biom_fp)
    else:
        # assume BIOM
        if not table_path.exists():
            eprint(f"[error] Table not found: {table_path}")
            return 2
        biom_fp = table_path
        table_sample_ids = read_sample_ids_from_biom(biom_fp)

    if not table_sample_ids:
        eprint("[error] No sample IDs found in table. Check that the input is a valid feature table.")
        return 2

    control_prefixes = [pfx.strip() for pfx in args.control_prefixes.split(",") if pfx.strip()]

    n_total, n_matched, n_unmatched = build_metadata(
        table_sample_ids=table_sample_ids,
        source_df=source_df,
        source_id_col=args.source_id_column,
        out_fp=out_path,
        marker=args.marker,
        key_regex=args.key_regex,
        control_prefixes=control_prefixes,
        control_sampletype_value=args.control_sampletype_value,
        sampletype_column=args.sampletype_column,
        write_mapping_report=mapping_report_path,
    )

    eprint(f"[done] Wrote: {out_path}")
    eprint(f"[summary] table_samples={n_total} matched_to_source={n_matched} unmatched={n_unmatched}")
    if mapping_report_path:
        eprint(f"[done] Mapping report: {mapping_report_path}")

    # Helpful warning if many unmatched
    if n_unmatched > 0:
        eprint(
            "[warn] Some table sample IDs did not match source metadata. "
            "This is expected for controls and any samples missing from source metadata.\n"
            "       Inspect the mapping report if provided."
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
