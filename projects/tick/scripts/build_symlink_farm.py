#!/usr/bin/env python3
"""
Tick project SYMLINK FARM builder
====================================

Reads tick_read_inventory.csv and plans a clean marker-organized symlink
tree plus QIIME2 manifest files. Never modifies the original FASTQ files.

By default runs in DRY-RUN mode — no symlinks are created, only a plan
is printed. Pass --execute to actually create the symlinks.

Design:
  - Single source of truth: the inventory CSV. Edit it before running
    (e.g. add 'skip' to the notes column) to override decisions.
  - Normalizes known cosmetic anomalies so biological samples aren't
    fragmented: eq -> EQ, Deer143-2 -> Deer-143-2, DEER -> Deer.
  - Resolves the 241209 P3-5 physical mirror (same files present in both
    /16S/ and /cytb/ trees) by preferring the /16S/ tree path as source.
  - KEEPS ALL RESEQUENCING VERSIONS. Each physical sequencing gets its
    own QIIME2 sample ID tagged with the run: biological_id.batch_short
    (e.g. ALHT-Deer-64.241209_P3-5). This lets you decide which version
    to use based on actual DADA2 retention rather than pre-committing.

Outputs (default layout):
    reads_staged/
      16S/                      <symlinks to 16S FASTQs, one pair per sample>
      cytb/                     <symlinks to cytb FASTQs, one pair per sample>
      16S_manifest.tsv          QIIME2 PairedEndFastqManifestPhred33V2
      cytb_manifest.tsv         QIIME2 PairedEndFastqManifestPhred33V2
      16S_metadata_template.tsv starter QIIME2 sample metadata
      cytb_metadata_template.tsv starter QIIME2 sample metadata
      staging_audit.csv         full audit: symlink -> source, with decision column

Usage:
    # Dry run: print the plan, create nothing
    python3 build_symlink_farm.py

    # Execute: create the symlink farm + manifests + metadata templates
    python3 build_symlink_farm.py --execute

    # Override paths
    python3 build_symlink_farm.py --inventory /path/to.csv --stage-dir /path

No non-standard Python dependencies.
"""

import re
import csv
import sys
import argparse
from pathlib import Path
from collections import defaultdict, Counter


# ============================================================
# Defaults
# ============================================================
DEFAULT_INVENTORY = Path(
    "/home/users/sne24/meedlab/tick_project/tick_read_inventory.csv"
)
DEFAULT_STAGE_DIR = Path(
    "/home/users/sne24/meedlab/tick_project/reads_staged"
)

# For the 241209 P3-5 combined demux, the same file is present in both
# the 16S and cytb filesystem trees. Prefer this tree as canonical source.
CANONICAL_MIRROR_TREE = "16S"


# ============================================================
# Normalization
# ============================================================

def normalize_sample_id(cleaned_id):
    """Fix cosmetic anomalies so biological samples aren't fragmented.

    Returns (normalized_id, list_of_fixes_applied).
    """
    fixes = []
    out = cleaned_id

    # eq###### -> EQ######  (Batch-02 lowercased the VDL prefix)
    if re.match(r"^eq\d", out):
        out = "EQ" + out[2:]
        fixes.append("eq->EQ")

    # ALHT-DEER-## -> ALHT-Deer-##  (one batch used all-caps DEER)
    if "ALHT-DEER" in out:
        out = out.replace("ALHT-DEER", "ALHT-Deer")
        fixes.append("DEER->Deer")

    # ALHT-Deer### -> ALHT-Deer-###  (missing hyphen at one submission)
    m = re.match(r"^(ALHT-Deer)(\d.*)$", out)
    if m:
        out = f"{m.group(1)}-{m.group(2)}"
        fixes.append("insert_hyphen")

    return out, fixes


# ============================================================
# Plan building
# ============================================================

def build_plan(rows, canonical_mirror_tree):
    """Return list of planned symlink operations.

    Each planned op is a dict with:
        marker                : '16S' or 'cytb'
        qiime_sample_id       : biological_id.batch_short (after normalization)
        normalized_id         : biological_id after cosmetic normalization
        original_sample_raw   : untouched name from disk
        batch_short, batch_dir
        tree, s_number, lane
        r1_src, r2_src        : Path objects to the original FASTQs
        r1_dst, r2_dst        : Path objects where symlinks will land
        is_control            : bool
        control_kind
        normalization_fixes   : str (semicolon-joined)
        skipped               : bool
        skip_reason           : str
    """
    ops = []

    # First pass: figure out which rows correspond to the 241209 P3-5 mirror
    # (same file physically present in both trees) so we can deduplicate.
    # Key on (batch_dir, original_sample_raw, s_number) — if this tuple
    # occurs in BOTH trees, it's a mirror and we keep only the canonical one.
    seen_by_phys_key = defaultdict(list)
    for r in rows:
        phys_key = (r["batch_dir"], r["original_sample_raw"], r["s_number"])
        seen_by_phys_key[phys_key].append(r)

    mirror_keys = {k for k, rs in seen_by_phys_key.items() if len(rs) > 1}

    for r in rows:
        phys_key = (r["batch_dir"], r["original_sample_raw"], r["s_number"])
        skipped = False
        skip_reason = ""

        # Mirror dedup: if this row is part of a mirror pair, keep only the
        # copy whose tree matches CANONICAL_MIRROR_TREE.
        if phys_key in mirror_keys and r["tree"] != canonical_mirror_tree:
            skipped = True
            skip_reason = f"mirror_of_{canonical_mirror_tree}_copy"

        # Respect user-added 'skip' hint in the notes column
        if (r.get("notes") or "").strip().lower().startswith("skip"):
            skipped = True
            skip_reason = "user_marked_skip"

        # Incomplete pair - should be 0 everywhere but defend anyway
        if r["paired_complete"] != "yes":
            skipped = True
            skip_reason = "incomplete_pair"

        if not r["r1_path"] or not r["r2_path"]:
            skipped = True
            skip_reason = skip_reason or "missing_path"

        marker = r["marker_guess"]
        if marker not in ("16S", "cytb"):
            skipped = True
            skip_reason = skip_reason or f"unknown_marker_{marker}"

        normalized_id, fixes = normalize_sample_id(r["cleaned_sample_id"])
        qiime_sample_id = f"{normalized_id}.{r['batch_short']}"

        op = {
            "marker": marker,
            "qiime_sample_id": qiime_sample_id,
            "normalized_id": normalized_id,
            "original_sample_raw": r["original_sample_raw"],
            "batch_short": r["batch_short"],
            "batch_dir": r["batch_dir"],
            "tree": r["tree"],
            "s_number": r["s_number"],
            "lane": r["lane"],
            "r1_src": Path(r["r1_path"]) if r["r1_path"] else None,
            "r2_src": Path(r["r2_path"]) if r["r2_path"] else None,
            "r1_dst": None,
            "r2_dst": None,
            "is_control": bool(r["is_control"]),
            "control_kind": r["control_kind"],
            "normalization_fixes": ";".join(fixes),
            "skipped": skipped,
            "skip_reason": skip_reason,
        }
        ops.append(op)

    return ops


def detect_sample_id_collisions(ops):
    """Find cases where multiple kept ops share the same qiime_sample_id.

    This indicates a bug in batch labeling (likely an outdated inventory CSV
    from before the date-prefix fix) which would cause silent data loss in
    the manifest writer. Mark all but the first collision as skipped with
    a clear reason.
    """
    by_qid = defaultdict(list)
    for op in ops:
        if not op["skipped"]:
            by_qid[(op["marker"], op["qiime_sample_id"])].append(op)

    collisions = {k: v for k, v in by_qid.items() if len(v) > 1}
    for (marker, qid), ops_group in collisions.items():
        for op in ops_group[1:]:  # keep first, skip rest
            op["skipped"] = True
            op["skip_reason"] = "duplicate_qiime_sample_id_rerun_inventory"

    return collisions


def assign_destinations(ops, stage_dir):
    """Fill in r1_dst/r2_dst for each non-skipped op.

    Also detects destination collisions (would two ops want the same
    symlink name?) and marks the duplicates as skipped.
    """
    claimed = {}  # dst path -> op that got it first
    collisions = 0

    for op in ops:
        if op["skipped"]:
            continue
        base = f"{op['qiime_sample_id']}_S{op['s_number']}_L{int(op['lane']):03d}"
        r1_dst = stage_dir / op["marker"] / f"{base}_R1_001.fastq.gz"
        r2_dst = stage_dir / op["marker"] / f"{base}_R2_001.fastq.gz"

        if r1_dst in claimed or r2_dst in claimed:
            op["skipped"] = True
            op["skip_reason"] = "destination_collision"
            collisions += 1
            continue

        claimed[r1_dst] = op
        claimed[r2_dst] = op
        op["r1_dst"] = r1_dst
        op["r2_dst"] = r2_dst

    return collisions


# ============================================================
# Execution
# ============================================================

def create_symlinks(ops, execute):
    """If execute=True, actually create the symlinks. Otherwise no-op."""
    created = 0
    already_ok = 0
    errors = []

    for op in ops:
        if op["skipped"]:
            continue

        for src, dst in [(op["r1_src"], op["r1_dst"]),
                         (op["r2_src"], op["r2_dst"])]:
            if not execute:
                continue

            try:
                dst.parent.mkdir(parents=True, exist_ok=True)
                if dst.is_symlink():
                    # Idempotent: if it already points to the right place, skip
                    if dst.resolve() == src.resolve():
                        already_ok += 1
                        continue
                    else:
                        errors.append(f"symlink exists with different target: {dst} -> {dst.resolve()} (wanted {src})")
                        continue
                if dst.exists():
                    errors.append(f"regular file already at destination: {dst}")
                    continue

                dst.symlink_to(src.resolve())
                created += 1
            except OSError as e:
                errors.append(f"failed {src} -> {dst}: {e}")

    return created, already_ok, errors


def write_manifests(ops, stage_dir, execute):
    """Write QIIME2 PairedEndFastqManifestPhred33V2 per marker."""
    by_marker = defaultdict(list)
    for op in ops:
        if op["skipped"]:
            continue
        by_marker[op["marker"]].append(op)

    written = []
    for marker, marker_ops in sorted(by_marker.items()):
        manifest_path = stage_dir / f"{marker}_manifest.tsv"
        if execute:
            manifest_path.parent.mkdir(parents=True, exist_ok=True)
            with manifest_path.open("w") as fh:
                fh.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
                # Dedup: if a qiime_sample_id somehow appears twice (shouldn't
                # after our uniqueness tag, but defend) we keep only first.
                seen = set()
                for op in sorted(marker_ops, key=lambda o: o["qiime_sample_id"]):
                    if op["qiime_sample_id"] in seen:
                        continue
                    seen.add(op["qiime_sample_id"])
                    fh.write(f"{op['qiime_sample_id']}\t{op['r1_dst']}\t{op['r2_dst']}\n")
        written.append((marker, manifest_path, len(marker_ops)))

    return written


def write_metadata_templates(ops, stage_dir, execute):
    """Write a starter QIIME2 sample metadata TSV per marker.

    Includes every column we can derive from the inventory. The researcher
    will extend this with experimental condition columns.
    """
    by_marker = defaultdict(list)
    for op in ops:
        if op["skipped"]:
            continue
        by_marker[op["marker"]].append(op)

    written = []
    header = [
        "sample-id",              # must be first for QIIME2
        "biological_sample_id",   # unnormalized-but-cleaned version for grouping reseqs
        "normalized_id",          # normalized version
        "original_sample_raw",
        "batch_short",
        "batch_dir",
        "run_date",               # YYMMDD extracted from batch_short
        "tree_source",
        "marker",
        "is_control",
        "control_kind",
        "s_number",
        "lane",
        "normalization_fixes",
    ]
    for marker, marker_ops in sorted(by_marker.items()):
        meta_path = stage_dir / f"{marker}_metadata_template.tsv"
        if execute:
            meta_path.parent.mkdir(parents=True, exist_ok=True)
            with meta_path.open("w") as fh:
                fh.write("\t".join(header) + "\n")
                # QIIME2 convention: a "categorical"/"numeric" types row can go
                # second, but it's optional. Leave it out to keep simple.
                for op in sorted(marker_ops, key=lambda o: o["qiime_sample_id"]):
                    run_date = op["batch_short"].split("_", 1)[0]
                    row = [
                        op["qiime_sample_id"],
                        op["normalized_id"],
                        op["normalized_id"],
                        op["original_sample_raw"],
                        op["batch_short"],
                        op["batch_dir"],
                        run_date,
                        op["tree"],
                        op["marker"],
                        "yes" if op["is_control"] else "no",
                        op["control_kind"],
                        op["s_number"],
                        op["lane"],
                        op["normalization_fixes"],
                    ]
                    fh.write("\t".join(row) + "\n")
        written.append((marker, meta_path, len(marker_ops)))

    return written


def write_audit(ops, stage_dir, execute):
    """Write the full audit CSV: every inventory row, whether kept or skipped."""
    audit_path = stage_dir / "staging_audit.csv"
    if execute:
        audit_path.parent.mkdir(parents=True, exist_ok=True)
        with audit_path.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "marker", "qiime_sample_id", "normalized_id",
                "original_sample_raw", "batch_short", "batch_dir", "tree",
                "s_number", "lane",
                "r1_src", "r2_src", "r1_dst", "r2_dst",
                "is_control", "control_kind", "normalization_fixes",
                "skipped", "skip_reason",
            ])
            for op in ops:
                writer.writerow([
                    op["marker"], op["qiime_sample_id"], op["normalized_id"],
                    op["original_sample_raw"], op["batch_short"],
                    op["batch_dir"], op["tree"],
                    op["s_number"], op["lane"],
                    op["r1_src"] or "", op["r2_src"] or "",
                    op["r1_dst"] or "", op["r2_dst"] or "",
                    "yes" if op["is_control"] else "",
                    op["control_kind"],
                    op["normalization_fixes"],
                    "yes" if op["skipped"] else "",
                    op["skip_reason"],
                ])
    return audit_path


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    ap.add_argument("--stage-dir", type=Path, default=DEFAULT_STAGE_DIR)
    ap.add_argument("--canonical-mirror-tree", choices=["16S", "cytb"],
                    default=CANONICAL_MIRROR_TREE,
                    help="For the 241209 P3-5 demux that is physically present in "
                         "both filesystem trees, use this tree as the symlink source")
    ap.add_argument("--execute", action="store_true",
                    help="Actually create symlinks, manifests, and metadata files. "
                         "Without this flag, the script runs in dry-run mode.")
    args = ap.parse_args()

    if not args.inventory.is_file():
        sys.exit(f"ERROR: inventory CSV not found: {args.inventory}")

    with args.inventory.open() as fh:
        rows = list(csv.DictReader(fh))

    print(f"Loaded {len(rows)} rows from {args.inventory}")
    print(f"Stage dir:              {args.stage_dir}")
    print(f"Canonical mirror tree:  {args.canonical_mirror_tree}")
    print(f"Mode:                   {'EXECUTE' if args.execute else 'DRY-RUN'}")
    print()

    # 1) Build plan
    ops = build_plan(rows, args.canonical_mirror_tree)

    # 1.5) Detect qiime_sample_id collisions (usually signals stale inventory)
    qid_collisions = detect_sample_id_collisions(ops)

    # 2) Assign destination paths and detect filesystem-level collisions
    collisions = assign_destinations(ops, args.stage_dir)

    # 3) Execute (or no-op in dry-run)
    created, already_ok, errors = create_symlinks(ops, args.execute)

    # 4) Write manifests + metadata templates + audit
    manifests = write_manifests(ops, args.stage_dir, args.execute)
    metas = write_metadata_templates(ops, args.stage_dir, args.execute)
    audit_path = write_audit(ops, args.stage_dir, args.execute)

    # ============================================================
    # Summary
    # ============================================================
    kept = [o for o in ops if not o["skipped"]]
    skipped = [o for o in ops if o["skipped"]]

    print("=" * 64)
    print("PLAN SUMMARY")
    print("=" * 64)
    print(f"Inventory rows:              {len(ops)}")
    print(f"  to be symlinked:           {len(kept)} pairs ({len(kept)*2} files)")
    print(f"  skipped:                   {len(skipped)}")
    print()

    skip_counts = Counter(o["skip_reason"] for o in skipped)
    if skip_counts:
        print("Skip reasons:")
        for reason, n in skip_counts.most_common():
            print(f"  {reason:35} {n}")
        print()

    by_marker = Counter(o["marker"] for o in kept)
    print("Will be staged by marker:")
    for m, n in sorted(by_marker.items()):
        print(f"  {m:6} {n} pairs")
    print()

    fixes_counts = Counter()
    for o in kept:
        for f in o["normalization_fixes"].split(";"):
            if f:
                fixes_counts[f] += 1
    if fixes_counts:
        print("Normalizations applied:")
        for f, n in fixes_counts.most_common():
            print(f"  {f:20} {n}")
        print()

    controls = sum(1 for o in kept if o["is_control"])
    print(f"Controls in staging plan: {controls}")

    # How many biological samples have multiple kept runs?
    bio_to_runs = defaultdict(set)
    for o in kept:
        if not o["is_control"]:
            bio_to_runs[(o["normalized_id"], o["marker"])].add(o["batch_short"])
    multi_run = {k: v for k, v in bio_to_runs.items() if len(v) > 1}
    print(f"Biological samples with multiple kept runs (resequencings preserved): "
          f"{len(multi_run)}")
    if multi_run:
        print("  (first 10 shown — see staging_audit.csv for full list)")
        for (bio, marker), runs in sorted(multi_run.items())[:10]:
            print(f"    {bio:30} [{marker}]  runs: {','.join(sorted(runs))}")
    print()

    if qid_collisions:
        print("!" * 64)
        print(f"STOP: {len(qid_collisions)} QIIME2 sample-id collisions detected.")
        print("This means two different physical sequencings would have been")
        print("assigned the same QIIME2 sample ID — silent data loss.")
        print()
        print("Most likely cause: your tick_read_inventory.csv is from the")
        print("OLD version of tick_read_inventory.py (before the date-prefix")
        print("fix that disambiguates 241209_P3-5 from 241219_P3-5).")
        print()
        print("Please re-run tick_read_inventory.py with the latest version,")
        print("then re-run this script.")
        print()
        print("First few colliding sample IDs:")
        for (marker, qid), ops_group in list(qid_collisions.items())[:5]:
            print(f"  {marker} :: {qid}  ({len(ops_group)} rows)")
            for op in ops_group:
                print(f"      {op['batch_dir']}  S{op['s_number']}")
        print("!" * 64)
        print()

    if collisions:
        print(f"WARNING: {collisions} destination collisions. Review staging_audit.csv.")
        print()

    if errors:
        print(f"ERRORS during symlink creation: {len(errors)}")
        for e in errors[:10]:
            print(f"  {e}")
        print()

    if args.execute:
        print("EXECUTED:")
        print(f"  symlinks created:    {created}")
        print(f"  symlinks already OK: {already_ok}")
        print(f"  errors:              {len(errors)}")
        for m, p, n in manifests:
            print(f"  {m} manifest:  {p} ({n} samples)")
        for m, p, n in metas:
            print(f"  {m} metadata:  {p}")
        print(f"  audit log:     {audit_path}")
    else:
        print("DRY-RUN — no files were created.")
        print(f"Would write:")
        print(f"  {args.stage_dir}/16S/         ({by_marker.get('16S',0)} symlinked pairs)")
        print(f"  {args.stage_dir}/cytb/        ({by_marker.get('cytb',0)} symlinked pairs)")
        print(f"  {args.stage_dir}/16S_manifest.tsv")
        print(f"  {args.stage_dir}/cytb_manifest.tsv")
        print(f"  {args.stage_dir}/16S_metadata_template.tsv")
        print(f"  {args.stage_dir}/cytb_metadata_template.tsv")
        print(f"  {args.stage_dir}/staging_audit.csv")
        print()
        print("Re-run with --execute to create the symlink farm.")


if __name__ == "__main__":
    main()
