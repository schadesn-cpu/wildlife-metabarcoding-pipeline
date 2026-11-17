#!/usr/bin/env python3
"""
build_single_tick_v2.py
=======================
Rebuild results/single_tick_v2/single_tick_v2_full.tsv (n=103) FROM SOURCE.

WHY THIS EXISTS
---------------
The single-tick file is the per-tick prevalence anchor for the manuscript, but
its builder script was lost (the reproducibility gap noted in
Tick_Ground_Truth_7_8_2026.md). This script regenerates it from the canonical
QIIME2 feature tables + BLAST-verified taxonomies + the crosswalk, so n=103 and
every 0/1 call become reproducible from raw instead of a frozen mystery file.

HOW IT PROVES ITSELF
--------------------
It self-checks two ways, so you never have to trust it blind:
  1. GOLDEN COUNTS  — per-column positive counts must match the known-good file
     (embedded below, taken from the current single_tick_v2_full.tsv).
  2. CELL DIFF      — if the old file is present, every cell of the rebuild is
     compared to it and any mismatch is printed (barcode + column).
If both pass, the rebuild reproduces the file exactly and the gap is closed.
If they don't, the script tells you WHICH specimens/columns differ — that is a
map to the remaining mapping question, not a failure.

WHAT IS CONFIDENT vs WHAT TO CONFIRM ON FIRST RUN
-------------------------------------------------
Confident (verified together this session):
  - output schema + column order
  - 16S / 18S pathogen taxon matches (genus / species strings)
  - the golden counts
Confirm on first run (marked CONFIRM below):
  - exact crosswalk column names (Num, Num Ticks, per-marker Sample cols)
  - the cytb host taxon strings (host mapping is marker-specific)
  - the specimen<->biom-sample join (16S Sample col is authoritative per ground
    truth; 18S/cytb are rebuilt here by barcode match and CHECKED by the diff)
  - MIN_READS: the per-call read threshold was never documented. Run with
    --discover to have the script find the threshold that reproduces golden.

USAGE
-----
  # dry structure check, no writing:
  python build_single_tick_v2.py --discover
  # write once it reproduces:
  python build_single_tick_v2.py --out results/single_tick_v2/single_tick_v2_full.REBUILD.tsv

Run from the tick_project root (paths below are relative to it).
"""
from __future__ import annotations
import argparse, glob, re, sys, tempfile, zipfile
from pathlib import Path
import pandas as pd

# ---------------------------------------------------------------------------
# CONFIG  — paths are the canonical ones from Tick_Ground_Truth_7_8_2026.md
# ---------------------------------------------------------------------------
CROSSWALK_XLSX  = "metadata/tick_full_metadata.xlsx"
CROSSWALK_SHEET = "All Run Data - w Metadata"          # CONFIRM sheet name

# feature tables (.qza) and BLAST-verified taxonomies (.tsv) per marker
MARKERS = {
    "16S": dict(
        # canonical per ground truth: calls come off the DECONTAMINATED table,
        # not raw dada2/table.qza (contaminant reads stripped before calling)
        table="results/decontam/16S/table_decontam_16S.qza",
        taxonomy="results/16S/all/taxonomy_verified/taxonomy_16S_blast_verified.tsv",
        tax_col="Taxon", fallback_col=None,
        sample_col="Sample",        # crosswalk block 1 = 16S (confirmed)
    ),
    "18S": dict(
        table="qiime2/18S/all/dada2/table.qza",
        taxonomy="results/18S/all/taxonomy_verified/taxonomy_18S_blast_verified.tsv",
        tax_col="Taxon", fallback_col=None,
        sample_col="Sample.1",      # crosswalk block 2 = 18S (confirmed)
    ),
    "cytb": dict(
        table="qiime2/cytb/all/dada2/table.qza",
        taxonomy="results/cytb/all/taxonomy_verified/taxonomy_blast_verified.tsv",
        # cytb verified id is in Taxon_verified; Taxon holds the OLD wrong label
        tax_col="Taxon_verified", fallback_col="Taxon",
        sample_col="Sample.2",      # crosswalk block 3 = cytB (confirmed)
    ),
}

# Taxon match rules: output_column -> list of case-insensitive substrings.
# A specimen is positive for a column if any matching ASV has reads >= MIN_READS.
TARGETS = {
    # ---- 16S pathogens/endosymbionts (CONFIDENT) ----
    # Anchored on the GENUS rank, matching the 16S_Borrelia line below.
    # A bare substring matches the rank ABOVE the genus: "Anaplasma" is inside
    # "Anaplasmataceae" and caught 9 Wolbachia ASVs (58 libraries); "Rickettsia"
    # is inside "Rickettsiales" and caught 1,758 mitochondrial ASVs (84
    # libraries). Verified against the decontam table 2026-09-09: Bartonella,
    # Borreliella and Francisella are unchanged by anchoring.
    "16S_Rickettsia":  ["g__Rickettsia;"],
    "16S_Anaplasma":   ["g__Anaplasma;"],
    "16S_Francisella": ["g__Francisella;"],
    "16S_Bartonella":  ["g__Bartonella;"],
    "16S_Borreliella": ["g__Borreliella;"],
    # Wolbachia as its own endosymbiont row: item 6 of the drop-list follow-up.
    # These are deer ked (Lipoptena cervi) endosymbionts, 98.9-100% identity,
    # and are one of the five evidence lines for the ked (ground truth section 6).
    "16S_Wolbachia":   ["g__Wolbachia;"],
    "16S_Borrelia":    [";Borrelia", "g__Borrelia", "Borrelia "],  # avoid matching Borreliella
    # ---- 18S piroplasms (CONFIDENT) ----
    "18S_Babesia_odocoilei": ["Babesia_odocoilei", "odocoilei"],
    "18S_Babesia_microti":   ["Babesia_microti", "microti"],
    # generic Babesia = Babesiidae/Babesia NOT resolved to odocoilei/microti:
    "18S_Babesia_sp":        ["__BABESIA_OTHER__"],   # handled specially below
    # ---- cytb hosts (CONFIRM strings against your cytb verified taxonomy) ----
    "cytb_Moose":              ["Alces"],
    "cytb_Deer":               ["Odocoileus"],
    "cytb_Cervidae_unresolved":["Cervidae"],           # family-level, unresolved
    "cytb_Deer_ked":           ["Lipoptena", "Hippoboscidae"],
    "cytb_Bird":               ["Aves", "Geothlypis", "Turdus", "Carduelis", "Spinus",
                                "Junco", "Melospiza", "Meleagris", "Corvus", "Setophaga",
                                "Vermivora", "Parulidae"],
    "cytb_Bear":               ["Ursus"],
    "cytb_Porcupine":          ["Erethizon"],
    "cytb_Rabbit":             ["Sylvilagus"],
    "cytb_Opossum":            ["Didelph"],
    "cytb_Raccoon":            ["Procyon"],
    # Human. NOT a confirmed blood meal: human cytb cannot distinguish feeding
    # from handler or laboratory contamination. Human is 0.41% of all cytb reads
    # and 252 of 837 biological samples carry trace human, so a detection here
    # means "human-positive above baseline". Reported, never called a meal.
    # See Human_Checkpoint_Notes, 2026-07-09.
    "cytb_Human":              ["Homo"],
}
# NOTE on cytb Moose/Deer: ground truth says cytb cannot reliably separate deer
# from moose. These columns reproduce the classifier-level calls as-is; report
# cervids at family level in the manuscript regardless of what these split into.

MIN_READS_DEFAULT = 1

META_COLS  = ["Barcode", "State", "Town", "Host Spp", "State_clean", "barcode"]
CALL_ORDER = list(TARGETS.keys())   # emitted in this order after META_COLS
OUT_ORDER  = ["Barcode", "State", "Town", "Host Spp", "State_clean", "barcode",
              "16S_Rickettsia","16S_Anaplasma","16S_Francisella","16S_Bartonella",
              "16S_Borreliella","16S_Borrelia","16S_Wolbachia","18S_Babesia_sp",
              "18S_Babesia_odocoilei","18S_Babesia_microti",
              "cytb_Moose","cytb_Cervidae_unresolved","cytb_Deer","cytb_Bird",
              "cytb_Deer_ked","cytb_Bear","cytb_Porcupine","cytb_Rabbit",
              "cytb_Opossum","cytb_Raccoon","cytb_Human"]

# GOLDEN per-column positive counts from the current single_tick_v2_full.tsv
GOLDEN = {"__n__": 103,
    "16S_Rickettsia": 96, "16S_Anaplasma": 18, "16S_Francisella": 2,
    "16S_Bartonella": 13, "16S_Borreliella": 32, "16S_Borrelia": 2,
    "18S_Babesia_sp": 0, "18S_Babesia_odocoilei": 12, "18S_Babesia_microti": 4,
    "cytb_Moose": 75, "cytb_Cervidae_unresolved": 76, "cytb_Deer": 85,
    "cytb_Bird": 4, "cytb_Deer_ked": 9, "cytb_Bear": 3, "cytb_Porcupine": 0,
    "cytb_Rabbit": 0, "cytb_Opossum": 2, "cytb_Raccoon": 1}

EXISTING_FILE = "results/single_tick_v2/single_tick_v2_full.tsv"  # for cell-diff


# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------
def load_qza_biom(qza_path: str) -> pd.DataFrame:
    """Extract a QIIME2 .qza feature table to an ASV x sample DataFrame."""
    import biom
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(qza_path) as z:
            z.extractall(td)
        biom_path = glob.glob(f"{td}/*/data/feature-table.biom")[0]
        return biom.load_table(biom_path).to_dataframe(dense=True)  # rows=ASV, cols=sample


def load_taxonomy(tsv_path: str, tax_col: str = "Taxon",
                  fallback_col: str | None = None) -> pd.Series:
    """Load a BLAST-verified taxonomy TSV -> Series(feature_id -> identity string).

    Reads by HEADER NAME (files vary: 16S/18S use 'Taxon' with the correction
    written in; cytb keeps the wrong label in 'Taxon' and the verified species in
    'Taxon_verified'). Uses `tax_col`, falling back to `fallback_col` where empty.
    """
    df = pd.read_csv(tsv_path, sep="\t")
    df.columns = [str(c).strip() for c in df.columns]
    id_col = df.columns[0]  # 'Feature ID' or 'Feature_ID'
    if tax_col not in df.columns:
        sys.exit(f"CONFIRM: '{tax_col}' not in {tsv_path}. Columns: {list(df.columns)}")
    ident = df[tax_col].astype(str).replace({"nan": "", "None": ""})
    if fallback_col and fallback_col in df.columns:
        fb = df[fallback_col].astype(str)
        ident = ident.where(ident.str.strip() != "", fb)
    return pd.Series(ident.values, index=df[id_col].astype(str))


def biom_prefix(sample_id: str) -> str:
    """The crosswalk Sample value == biom sample id BEFORE the first dot."""
    return str(sample_id).split(".")[0]


# ---------------------------------------------------------------------------
# presence
# ---------------------------------------------------------------------------
def asvs_for(taxa: pd.Series, patterns: list[str], all_targets: dict) -> list[str]:
    """Feature IDs whose taxon matches any pattern (special-case generic Babesia)."""
    if patterns == ["__BABESIA_OTHER__"]:
        is_bab = taxa.str.contains("Babesia", case=False, na=False)
        named  = taxa.str.contains("odocoilei|microti", case=False, na=False)
        return list(taxa[is_bab & ~named].index)
    m = pd.Series(False, index=taxa.index)
    for p in patterns:
        m |= taxa.str.contains(p, case=False, na=False, regex=False)
    return list(taxa[m].index)


def presence_by_sample(table: pd.DataFrame, taxa: pd.Series, patterns: list[str],
                       min_reads: int, all_targets: dict) -> pd.Series:
    """Series(sample -> 0/1) = any matching ASV has >= min_reads in that sample."""
    ids = [f for f in asvs_for(taxa, patterns, all_targets) if f in table.index]
    if not ids:
        return pd.Series(0, index=table.columns, dtype=int)
    reads = table.loc[ids].sum(axis=0)          # reads per sample across matching ASVs
    return (reads >= min_reads).astype(int)


# ---------------------------------------------------------------------------
# crosswalk / specimen selection + join
# ---------------------------------------------------------------------------
def load_single_ticks(xlsx: str, sheet: str) -> pd.DataFrame:
    """Single-tick specimens, taken from crosswalk_4marker.tsv.

    I read the CROSSWALK, not the raw workbook. Reading the workbook bypassed every
    correction the crosswalk holds: the ALHT/LHT barcode clobber (17 specimens), the
    LHT-091 unmerge, USGSTICKS-2213, the registry fills, and the VDL exclusions. It
    also missed 24-3889 (blank Barcode cell in the workbook) and ALHT WTD_09 (which
    exists as a specimen only because LABEL_FIXES reassigns Num 247 to it), which is
    why the set landed at 91 when the rule says 93.

    is_single_tick is the ONE definition: tick count of exactly 1, nothing
    contradicting it, wild-caught in the study region. The xlsx and sheet arguments
    are kept so callers do not break, and ignored.
    """
    cw = pd.read_csv("results/metadata/crosswalk_4marker.tsv", sep="\t")
    if "is_single_tick" not in cw.columns:
        sys.exit("crosswalk_4marker.tsv has no is_single_tick column. "
                 "Run build_4marker_crosswalk.py first.")
    flag = cw["is_single_tick"].astype(str).str.lower().isin(["true", "1"])
    single = cw[flag].copy()
    # downstream code expects a Barcode column
    if "Barcode" not in single.columns:
        single["Barcode"] = single["sample_id"]
    print(f"[crosswalk] single-tick specimens (is_single_tick): {len(single)}")
    return single




# --- specimen identity via scheme + number (per ID_Reconciliation sheet) -------
# Reduces any name (crosswalk Barcode OR biom sample prefix) to a canonical
# (scheme, id) key. Survives: marker tags (-16s/-18s/-cytoB), suffix-order
# variation, case, USGSTICKS library expansion (-NN), and EQ sub-samples. This
# matches the biom prefix straight to the Barcode, so it does NOT depend on the
# 18S/cytb Sample columns (which the ground-truth doc flags as incomplete).
_MARKER_TAG = re.compile(r"-(16S|18S|CYTB|CYTOB|ITS)(?=-|$|_)")
_SCHEME_RULES = [
    ("ALHT",     r"ALHT[\s_-]*DEER[\s_-]*0*(\d+)"),
    ("USGSVT12", r"USGSTICKS-VT-12-([A-Z0-9]+)"),
    ("LHTDEER",  r"USGSTICKS-LHT-DEER-0*(\d+)"),
    ("LHT",      r"USGSTICKS-LHT0*(\d+)"),
    ("USGSTICK", r"USGS-TICK-0*(\d+)"),
    ("USGS",     r"USGSTICKS-0*(\d+)"),
    ("EQ",       r"(EQ\d+)"),
    ("TB",       r"TB[_-](\d+)"),
]


def core(name: str):
    """Canonical (scheme, id) key for a specimen name or biom prefix."""
    s = str(name).upper().strip()
    s = _MARKER_TAG.sub("", s)          # drop -16s/-18s/-cytoB wherever it sits
    s = s.replace("WTD", "DEER")
    for tag, pat in _SCHEME_RULES:
        m = re.search(pat, s)
        if m:
            g = m.group(1)
            return (tag, int(g) if g.isdigit() else g)
    m = re.search(r"^(\d[\d-]*)$", s)  # bare numeric (4065, 24-3775)
    if m:
        return ("NUM", m.group(1).replace("-", ""))
    return ("RAW", s)


def specimen_barcode_map(table_cols, single_barcodes) -> dict:
    """Map each biom sample-id -> specimen Barcode by matching core() keys.
    Multi-library biom ids (ALHT-Deer-04-1/-2) share a core and fold to one
    barcode (OR'd downstream). Warns on any barcode-core collision."""
    core_to_bc = {}
    for b in single_barcodes:
        k = core(b)
        if k in core_to_bc and core_to_bc[k] != b:
            print(f"[warn] core collision: {k} <- {core_to_bc[k]} and {b}")
        core_to_bc[k] = b
    out = {}
    for col in table_cols:
        bc = core_to_bc.get(core(biom_prefix(col)))
        if bc is not None:
            out[col] = bc
    return out


# ---------------------------------------------------------------------------
# assemble + validate
# ---------------------------------------------------------------------------
def build(min_reads: int) -> pd.DataFrame:
    single_all = load_single_ticks(CROSSWALK_XLSX, CROSSWALK_SHEET)   # all library rows
    single = single_all.drop_duplicates(subset="Barcode")            # 1 row/barcode
    barcodes = list(dict.fromkeys(single["Barcode"].astype(str)))    # unique, ordered
    print(f"[crosswalk] unique single-tick barcodes: {len(barcodes)}")

    # start the output frame from crosswalk metadata, keyed by Barcode
    meta = (single.set_index(single["Barcode"].astype(str))
                  .reindex(columns=[c for c in META_COLS if c in single.columns]))
    out = pd.DataFrame(index=barcodes)
    for c in META_COLS:
        out[c] = meta[c] if c in meta.columns else 0

    for marker, cfg in MARKERS.items():
        table = load_qza_biom(cfg["table"])
        taxa  = load_taxonomy(cfg["taxonomy"], cfg.get("tax_col", "Taxon"),
                              cfg.get("fallback_col"))
        samp2bc = specimen_barcode_map(table.columns, barcodes)
        matched = len(set(samp2bc.values()))
        print(f"[{marker}] biom samples matched to single ticks: {matched}/{len(barcodes)}")

        for col in [c for c in CALL_ORDER if c.startswith(marker)]:
            pres = presence_by_sample(table, taxa, TARGETS[col], min_reads, TARGETS)
            # fold sample-level presence up to specimen (Barcode), OR across libraries
            by_bc = {}
            for samp, v in pres.items():
                bc = samp2bc.get(samp)
                if bc is not None:
                    by_bc[bc] = max(by_bc.get(bc, 0), int(v))
            out[col] = [by_bc.get(bc, 0) for bc in out.index]

    return out.reindex(columns=OUT_ORDER).fillna(0)


def check_golden(df: pd.DataFrame) -> bool:
    ok = True
    if len(df) != GOLDEN["__n__"]:
        print(f"  [n] {len(df)} != {GOLDEN['__n__']}  <-- FAIL"); ok = False
    else:
        print(f"  [n] {len(df)} == {GOLDEN['__n__']}  OK")
    for col, want in GOLDEN.items():
        if col == "__n__":
            continue
        got = int(df[col].sum()) if col in df.columns else -1
        flag = "OK" if got == want else "<-- FAIL"
        if got != want:
            ok = False
        print(f"  [{col:26}] got {got:>3}  want {want:>3}  {flag}")
    return ok


def cell_diff(df: pd.DataFrame, existing_path: str) -> None:
    if not Path(existing_path).exists():
        print(f"[diff] {existing_path} not found — skipping cell diff")
        return
    old = pd.read_csv(existing_path, sep="\t").set_index("Barcode")
    new = df.set_index("Barcode")
    call_cols = [c for c in OUT_ORDER if c.startswith(("16S_", "18S_", "cytb_"))]
    diffs = 0
    CAP = 40  # don't flood the terminal; report the rest as a count
    def show(msg):
        nonlocal diffs
        if diffs < CAP:
            print(msg)
        elif diffs == CAP:
            print("  ... (further mismatches suppressed; see total below)")
        diffs += 1
    extra = [bc for bc in new.index if bc not in old.index]
    for bc in extra:
        show(f"  [extra row — in rebuild, not in existing 103] {bc}")
    for bc in new.index:
        if bc not in old.index:
            continue
        for c in call_cols:
            # A column the older file does not have is a NEW column, not a
            # difference. cytb_Human was added 2026-08-14.
            if c not in old.columns:
                continue
            o = 1 if str(old.loc[bc, c]).strip() in ("1", "1.0") else 0
            n = int(new.loc[bc, c])
            if o != n:
                show(f"  [cell] {bc:24} {c:24} old={o} new={n}")
    for bc in (set(old.index) - set(new.index)):
        show(f"  [missing row — in existing 103, not in rebuild] {bc}")
    print(f"[diff] total cell/row mismatches: {diffs}"
          + ("  <-- REPRODUCES EXACTLY" if diffs == 0 else ""))
    if extra:
        print(f"[diff] {len(extra)} 'extra' rows = the specimens the real 116->103 "
              f"filter dropped. These are the filter we still need to encode.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-reads", type=int, default=MIN_READS_DEFAULT,
                    help=f"reads to call a taxon present (default {MIN_READS_DEFAULT})")
    ap.add_argument("--discover", action="store_true",
                    help="try several thresholds and report which reproduces golden")
    ap.add_argument("--out", default=None, help="write rebuilt TSV to this path")
    args = ap.parse_args()

    if args.discover:
        for thr in (1, 2, 3, 5, 10):
            print(f"\n===== MIN_READS = {thr} =====")
            df = build(thr)
            check_golden(df)
        return 0

    print(f"===== building single_tick_v2 (min_reads={args.min_reads}) =====")
    df = build(args.min_reads)
    print("\n----- golden count check -----")
    golden_ok = check_golden(df)
    print("\n----- cell diff vs existing file -----")
    cell_diff(df, EXISTING_FILE)

    if args.out:
        df.to_csv(args.out, sep="\t", index=False)
        print(f"\n[write] {args.out}")
        if not golden_ok:
            print("[warn] golden counts did NOT all match — inspect diffs above "
                  "before trusting this file")
    else:
        print("\n[dry run] no --out given; nothing written")
    # A golden-count mismatch is a WARNING, not a failure. The goldens were
    # copied from the output file itself (see line ~136), so they assert only
    # that the file still matches itself and cannot tell a correction from a
    # regression. Seven of them differ from June because the taxonomy, the
    # crosswalk and the cytb remap have all been corrected since. Failing the
    # build on that stopped a caller from ever regenerating this file, which
    # is how it came to be two months stale.
    # Update the goldens once the read threshold is settled with L. Gordon,
    # and record which run they came from.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
