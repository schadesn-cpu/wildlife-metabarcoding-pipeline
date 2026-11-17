#!/usr/bin/env python3
"""
apply_18S_blast_verified_full.py

Re-collapse the 18S feature table using the BLAST-VERIFIED taxonomy.

WHY THIS EXISTS
---------------
scripts/apply_18S_blast_verified.py writes a corrected *taxonomy* file and stops.
Its own docstring says: "Nothing is deleted; the feature table is untouched."
Every downstream number is computed from TABLES, so corrections that live only in
a taxonomy file change nothing. That is why Lipoptena cervi (327,934 reads, 100%
BLAST) is absent from results/18S/all/taxonomy/taxonomy_relabund_L7_18S.tsv, and
why 23.2M reads of bacterial cross-amplification are still counted as eukaryotic
community.

This script is the missing second half. It mirrors the 16S pattern
(16S_apply_blast_audit_full.py): corrected feature-level taxonomy + feature table
-> corrected collapsed tables.

RULE ESTABLISHED 2026-08-03: a marker is not "applied" until a table NEWER than
its verified taxonomy exists. Check with:
    find results/<marker> -type f -newer <verified_taxonomy.tsv>

USAGE
-----
    # look first, write nothing
    python3 scripts/apply_18S_blast_verified_full.py --report

    # write outputs
    python3 scripts/apply_18S_blast_verified_full.py --apply

    # explicit paths if autodetect picks wrong
    python3 scripts/apply_18S_blast_verified_full.py --apply \
        --table results/cooccurrence/18S/biom_export/feature-table.tsv

OUTPUTS (results/18S/all/taxonomy_verified/)
--------------------------------------------
  counts_taxon_18S_BLASTcorrected.tsv    taxon x sample read counts, EVERYTHING
  relabund_taxon_18S_BLASTcorrected.tsv  same, relative abundance within sample
  community_18S_corrected.tsv            eukaryotic community: CROSSAMP, tick
                                         host, VERTEBRATE blood-meal host,
                                         mitochondrial refs and non-target
                                         arthropods all REMOVED
  tick_genus_per_sample_corrected.tsv    rebuild of tick_genus_per_sample.tsv
                                         (Ixodes / Dermacentor / other), plus
                                         Lipoptena, from CORRECTED taxonomy
  taxon_before_after_18S.tsv             what the corrections moved
  apply_18S_provenance.json              inputs, hashes, mtimes, counts

NOTHING IS DELETED. The exclusions are applied only to community_18S_corrected.tsv.
The counts/relabund tables retain every feature so host and cross-amplification
reads stay auditable.
"""

import argparse
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime

import pandas as pd

# --------------------------------------------------------------------------
# config
# --------------------------------------------------------------------------

VERIFIED_TAX = "results/18S/all/taxonomy_verified/taxonomy_18S_blast_verified.tsv"
CORRECTIONS = "results/18S/all/taxonomy_verified/blast_corrections_applied.tsv"
OUTDIR = "results/18S/all/taxonomy_verified"

# candidate feature tables, in preference order
TABLE_CANDIDATES = [
    "qiime2/18S/all/**/table*.qza",
    "results/18S/all/tables/*.qza",
    "results/cooccurrence/18S/biom_export/feature-table.biom",
    "results/cooccurrence/18S/biom_export/feature-table.tsv",
]

# Rows whose corrected taxon starts with this are bacterial reads picked up by
# the 18S primers. Not eukaryotic community. Excluded from community table only.
CROSSAMP_PREFIX = "CROSSAMP"

# Tick host genera. The specimen itself, not community. Excluded from community
# table, but broken out separately in tick_genus_per_sample_corrected.tsv.
TICK_GENERA = ["Ixodes", "Dermacentor", "Amblyomma", "Haemaphysalis", "Rhipicephalus"]

# Non-target arthropods that are real detections but not eukaryotic microbial
# community. Kept visible, excluded from community table.
OTHER_ARTHROPOD_MARKERS = ["Lipoptena", "Hippoboscidae", "Melophagus", "Icosta"]

# Vertebrate blood-meal DNA. Abundant in 18S from engorged ticks, and NOT
# eukaryotic microbial community. Added 2026-08-06 after a --report run showed
# Odocoileus (2.1M reads), Craniata, Tarsius and Archosauria sitting at the top
# of the "community" table.
VERTEBRATE_MARKERS = ["Craniata", "Chordata", "Vertebrata", "Mammalia", "Aves",
                      "Archosauria", "Odocoileus", "Tarsius", "Ursus", "Alces",
                      "Didelphis", "Erethizon", "Homo"]

# PR2 appends ":mito" to mitochondrial references. Off-target for an SSU survey
# whatever organism they came from.
MITO_SUFFIX = ":mito"


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def die(msg):
    print(f"\nERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def sha256(path, blocks=64):
    """Hash first N MB. Full hash on a 400MB tsv is not worth the wait."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for _ in range(blocks):
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()[:16]


def mtime(path):
    return datetime.fromtimestamp(os.path.getmtime(path)).isoformat(timespec="seconds")


def find_table(explicit=None):
    if explicit:
        if not os.path.exists(explicit):
            die(f"--table given but not found: {explicit}")
        return explicit
    for pattern in TABLE_CANDIDATES:
        hits = [h for h in glob.glob(pattern, recursive=True) if "raref" not in h]
        if hits:
            hits.sort(key=os.path.getmtime, reverse=True)
            return hits[0]
    die("no 18S feature table found. Pass one with --table.")


def load_table(path):
    """Return features x samples DataFrame of counts."""
    print(f"  loading table: {path}")

    if path.endswith(".qza"):
        tmp = tempfile.mkdtemp()
        try:
            subprocess.run(["unzip", "-qq", "-o", path, "-d", tmp], check=True)
            biom_hits = glob.glob(f"{tmp}/*/data/feature-table.biom")
            if not biom_hits:
                die(f"no feature-table.biom inside {path}")
            df = _load_biom(biom_hits[0])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        return df

    if path.endswith(".biom"):
        return _load_biom(path)

    if path.endswith((".tsv", ".txt")):
        # biom-exported tsv: line 1 is a '# Constructed from biom file' comment
        with open(path) as fh:
            first = fh.readline()
        skip = 1 if first.startswith("# Constructed") else 0
        df = pd.read_csv(path, sep="\t", skiprows=skip, index_col=0, low_memory=False)
        df.index.name = "Feature ID"
        return df

    die(f"unrecognized table format: {path}")


def _load_biom(path):
    try:
        import biom
    except ImportError:
        die("the biom package is not importable. Run inside the qiime2 env, "
            "or pass the exported .tsv with --table.")
    t = biom.load_table(path)
    df = t.to_dataframe(dense=True)
    df.index.name = "Feature ID"
    return df


def terminal_label(taxon):
    """Deepest non-empty rank in a semicolon lineage.

    18S lineages here are ragged (some skip ranks, some are bare binomials),
    so a fixed level like L6 or L7 silently mislabels. Deepest-available is
    honest about what each feature actually resolved to.
    """
    if not isinstance(taxon, str) or not taxon.strip():
        return "Unassigned"
    parts = [p.strip() for p in taxon.split(";")]
    parts = [p for p in parts if p and p not in ("__", "NA", "unassigned")]
    if not parts:
        return "Unassigned"

    def strip_rank(p):
        return p[3:] if len(p) > 3 and p[1:3] == "__" else p

    last = parts[-1]
    # A bare species epithet ("s__microti") is meaningless alone: Babesia
    # microti and Cryptosporidium microti would collide. Prepend the genus.
    if len(last) > 3 and last[:3] == "s__" and len(parts) >= 2:
        genus = strip_rank(parts[-2])
        epithet = strip_rank(last)
        if not epithet.lower().startswith(genus.lower()[:4]):
            last = f"{genus} {epithet}"
        else:
            last = epithet
    else:
        last = strip_rank(last)

    return last.replace("_", " ").strip() or "Unassigned"


def classify_row(taxon):
    """Bucket a lineage.

    Order matters: the most specific claim wins. A vertebrate mitochondrial
    reference is reported as vertebrate_host, not lumped into mitochondrial.
    """
    if not isinstance(taxon, str):
        return "community"
    if taxon.startswith(CROSSAMP_PREFIX):
        return "crossamp"
    for g in TICK_GENERA:
        if g in taxon:
            return "tick_host"
    for m in OTHER_ARTHROPOD_MARKERS:
        if m in taxon:
            return "other_arthropod"
    for v in VERTEBRATE_MARKERS:
        if v in taxon:
            return "vertebrate_host"
    if MITO_SUFFIX in taxon:
        return "mitochondrial"
    return "community"


def tick_genus_of(taxon):
    for g in TICK_GENERA:
        if isinstance(taxon, str) and g in taxon:
            return g
    return None


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--report", action="store_true",
                      help="show what would change, write nothing")
    mode.add_argument("--apply", action="store_true",
                      help="write the corrected tables")
    ap.add_argument("--table", default=None, help="explicit feature table path")
    ap.add_argument("--taxonomy", default=VERIFIED_TAX)
    ap.add_argument("--outdir", default=OUTDIR)
    args = ap.parse_args()

    print("=" * 74)
    print("18S apply: re-collapse feature table with BLAST-verified taxonomy")
    print("=" * 74)

    # ---- inputs -----------------------------------------------------------
    if not os.path.exists(args.taxonomy):
        die(f"verified taxonomy not found: {args.taxonomy}")

    print(f"\n[1] verified taxonomy: {args.taxonomy}")
    print(f"    mtime {mtime(args.taxonomy)}")
    tax = pd.read_csv(args.taxonomy, sep="\t", low_memory=False)

    idcol = tax.columns[0]
    taxcol = [c for c in tax.columns if c.lower() in ("taxon", "taxonomy")]
    if not taxcol:
        die(f"no Taxon column in {args.taxonomy}; saw {list(tax.columns)}")
    taxcol = taxcol[0]
    tax = tax.rename(columns={idcol: "feature_id", taxcol: "taxon"})
    tax = tax[["feature_id", "taxon"]].drop_duplicates("feature_id")
    print(f"    {len(tax):,} features")

    table_path = find_table(args.table)
    print(f"\n[2] feature table")
    tab = load_table(table_path)
    print(f"    {tab.shape[0]:,} features x {tab.shape[1]:,} samples")
    print(f"    mtime {mtime(table_path)}")

    # ---- alignment check --------------------------------------------------
    tab_ids = set(tab.index.astype(str))
    tax_ids = set(tax["feature_id"].astype(str))
    shared = tab_ids & tax_ids
    only_tab = tab_ids - tax_ids
    only_tax = tax_ids - tab_ids

    print(f"\n[3] alignment")
    print(f"    in both:          {len(shared):,}")
    print(f"    table only:       {len(only_tab):,}  (will label Unassigned)")
    print(f"    taxonomy only:    {len(only_tax):,}  (not in this table)")
    if len(shared) == 0:
        die("zero shared feature IDs. Wrong table for this taxonomy.")
    if len(shared) / max(len(tab_ids), 1) < 0.5:
        print("    WARNING: under half of table features have taxonomy. "
              "Check that the table matches the marker.")

    # ---- join -------------------------------------------------------------
    tab.index = tab.index.astype(str)
    lut = dict(zip(tax["feature_id"].astype(str), tax["taxon"]))
    lineage = pd.Series({fid: lut.get(fid, "Unassigned") for fid in tab.index})

    label = lineage.map(terminal_label)
    bucket = lineage.map(classify_row)

    # ---- what moved -------------------------------------------------------
    print(f"\n[4] corrections")
    if os.path.exists(CORRECTIONS):
        corr = pd.read_csv(CORRECTIONS, sep="\t", low_memory=False)
        ccol = corr.columns[0]
        corr_ids = set(corr[ccol].astype(str)) & set(tab.index)
        moved_reads = int(tab.loc[sorted(corr_ids)].values.sum()) if corr_ids else 0
        print(f"    {len(corr):,} corrections logged, "
              f"{len(corr_ids):,} present in this table")
        print(f"    reads under corrected features: {moved_reads:,}")
    else:
        corr = None
        print(f"    (no corrections log at {CORRECTIONS})")

    # ---- bucket summary ---------------------------------------------------
    read_by_bucket = tab.groupby(bucket.values).sum().sum(axis=1)
    total_reads = int(tab.values.sum())
    print(f"\n[5] read partition (total {total_reads:,})")
    for b in ["community", "tick_host", "vertebrate_host", "crossamp",
              "mitochondrial", "other_arthropod"]:
        r = int(read_by_bucket.get(b, 0))
        print(f"    {b:<18} {r:>15,}  ({100*r/max(total_reads,1):5.2f}%)")
    print("\n    community table excludes tick_host, vertebrate_host, crossamp,")
    print("    mitochondrial and other_arthropod. counts/relabund keep everything.")

    # ---- collapse ---------------------------------------------------------
    counts = tab.groupby(label.values).sum()
    counts.index.name = "Taxon"
    counts = counts.loc[counts.sum(axis=1).sort_values(ascending=False).index]

    comm_mask = (bucket == "community").reindex(tab.index).fillna(False).values
    comm = tab.loc[comm_mask].groupby(label[comm_mask].values).sum()
    comm.index.name = "Taxon"
    comm = comm.loc[comm.sum(axis=1).sort_values(ascending=False).index]

    denom = comm.sum(axis=0).replace(0, pd.NA)
    relabund = (comm / denom).fillna(0.0).astype(float)

    print(f"\n[6] collapsed")
    print(f"    all taxa:       {counts.shape[0]:,} x {counts.shape[1]:,}")
    print(f"    community taxa: {comm.shape[0]:,} x {comm.shape[1]:,}")
    print("\n    top 12 community taxa by reads:")
    for t, v in comm.sum(axis=1).head(12).items():
        print(f"      {t:<42} {int(v):>12,}")

    # ---- tick genus per sample -------------------------------------------
    tg = pd.Series({fid: tick_genus_of(lineage[fid]) for fid in tab.index})
    tick_rows = {}
    for g in TICK_GENERA:
        m = (tg == g).reindex(tab.index).fillna(False).values
        if m.sum():
            tick_rows[g] = tab.loc[m].sum(axis=0)
    ked_m = label.str.contains("Lipoptena", case=False, na=False).values
    if ked_m.sum():
        tick_rows["Lipoptena_cervi"] = tab.loc[ked_m].sum(axis=0)

    tickdf = pd.DataFrame(tick_rows)
    genus_cols = [c for c in tickdf.columns if c in TICK_GENERA]
    tickdf["tick_reads"] = tickdf[genus_cols].sum(axis=1) if genus_cols else 0
    tickdf["dominant"] = (
        tickdf[genus_cols].idxmax(axis=1).where(tickdf["tick_reads"] > 0, "none")
        if genus_cols else "none"
    )
    tickdf.index.name = "sample_id"

    print(f"\n[6b] tick-genus read totals (for reconciling library counts)")
    for g in TICK_GENERA:
        m = (tg == g).reindex(tab.index).fillna(False).values
        if m.sum():
            sub_ = tab.loc[m]
            print(f"    {g:<14} {int(m.sum()):>4} ASVs  {int(sub_.values.sum()):>12,} reads"
                  f"   max in one library: {int(sub_.sum(axis=0).max()):>10,}")

    print(f"\n[7] tick genus per sample (rebuilt from corrected taxonomy)")
    if genus_cols:
        for g in genus_cols:
            n = int((tickdf[g] > 0).sum())
            n100 = int((tickdf[g] >= 100).sum())
            print(f"    {g:<16} {n:>4} samples >0 reads, {n100:>4} at >=100")
    if "Lipoptena_cervi" in tickdf:
        n = int((tickdf["Lipoptena_cervi"] > 0).sum())
        n100 = int((tickdf["Lipoptena_cervi"] >= 100).sum())
        print(f"    {'Lipoptena cervi':<16} {n:>4} samples >0 reads, {n100:>4} at >=100")
    print("\n    NOTE: these are LIBRARIES, not specimens, and include controls.")
    print("    Join name_18S in crosswalk_4marker.tsv before quoting any count.")
    print("    2026-08-03: 25 Dermacentor libraries = 20 specimens + 5 controls.")

    # ---- before / after ---------------------------------------------------
    ba = None
    if corr is not None and len(corr):
        oldcol = [c for c in corr.columns if "old" in c.lower()]
        newcol = [c for c in corr.columns if "new" in c.lower()]
        if oldcol and newcol:
            ccol = corr.columns[0]
            sub = corr[corr[ccol].astype(str).isin(tab.index)].copy()
            if len(sub):
                sub["reads_in_table"] = [
                    int(tab.loc[str(f)].sum()) for f in sub[ccol]
                ]
                sub["old_label"] = sub[oldcol[0]].map(terminal_label)
                sub["new_label"] = sub[newcol[0]].map(terminal_label)
                ba = (sub.groupby(["old_label", "new_label"])["reads_in_table"]
                      .agg(["sum", "count"])
                      .rename(columns={"sum": "reads", "count": "features"})
                      .sort_values("reads", ascending=False)
                      .reset_index())
                print(f"\n[8] top label moves")
                for _, r in ba.head(10).iterrows():
                    print(f"    {r['old_label'][:26]:<26} -> {r['new_label'][:26]:<26} "
                          f"{int(r['reads']):>12,} reads ({int(r['features'])} feat)")

    # ---- write ------------------------------------------------------------
    if args.report:
        print("\n" + "=" * 74)
        print("REPORT MODE. Nothing written. Re-run with --apply to write.")
        print("=" * 74)
        return

    os.makedirs(args.outdir, exist_ok=True)
    written = []

    def w(df, name, **kw):
        p = os.path.join(args.outdir, name)
        df.to_csv(p, sep="\t", **kw)
        written.append(p)
        print(f"    wrote {p}")

    print(f"\n[9] writing to {args.outdir}")
    w(counts, "counts_taxon_18S_BLASTcorrected.tsv")
    w(relabund, "relabund_taxon_18S_BLASTcorrected.tsv")
    w(comm, "community_18S_corrected.tsv")
    w(tickdf, "tick_genus_per_sample_corrected.tsv")
    if ba is not None:
        w(ba, "taxon_before_after_18S.tsv", index=False)

    prov = {
        "script": os.path.basename(__file__),
        "run_at": datetime.now().isoformat(timespec="seconds"),
        "inputs": {
            "taxonomy": {"path": args.taxonomy, "mtime": mtime(args.taxonomy),
                         "sha256_head": sha256(args.taxonomy),
                         "n_features": int(len(tax))},
            "table": {"path": table_path, "mtime": mtime(table_path),
                      "sha256_head": sha256(table_path),
                      "shape": [int(tab.shape[0]), int(tab.shape[1])]},
            "corrections": ({"path": CORRECTIONS, "mtime": mtime(CORRECTIONS),
                             "n": int(len(corr))} if corr is not None else None),
        },
        "alignment": {"shared": len(shared), "table_only": len(only_tab),
                      "taxonomy_only": len(only_tax)},
        "reads": {"total": total_reads,
                  **{b: int(read_by_bucket.get(b, 0))
                     for b in ["community", "tick_host", "crossamp", "other_arthropod"]}},
        "exclusions": {"crossamp_prefix": CROSSAMP_PREFIX,
                       "tick_genera": TICK_GENERA,
                       "other_arthropod_markers": OTHER_ARTHROPOD_MARKERS,
                       "vertebrate_markers": VERTEBRATE_MARKERS,
                       "mito_suffix": MITO_SUFFIX},
        "outputs": written,
    }
    pp = os.path.join(args.outdir, "apply_18S_provenance.json")
    with open(pp, "w") as fh:
        json.dump(prov, fh, indent=2)
    print(f"    wrote {pp}")

    print("\n" + "=" * 74)
    print("DONE. Verify the apply actually landed:")
    print(f"  find results/18S -type f -newer {args.taxonomy} | head")
    print("Outputs above should appear. If nothing does, the write failed.")
    print("=" * 74)


if __name__ == "__main__":
    main()
