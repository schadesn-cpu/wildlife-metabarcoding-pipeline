#!/usr/bin/env python3
"""
Build the canonical Num-keyed marker crosswalk from the positional compiled sheet.

WHY THIS EXISTS: the 'Tick Rosetta Stone' sheet is a stub (only 16S-name is
populated; 18S/ITS/cytoB columns are 0/880). The real crosswalk is the
positional block layout in 'All Run Data - w Metadata':
  metadata cols 0-11, then four marker blocks each [Run Batch, Sample, Plate, Num, Amp]:
    16S  Sample = col 14
    18S  Sample = col 20
    cytb Sample = col 26
    ITS  Sample = col 32
  The universal key 'Num' (col 1) is stamped into every block (cols 16,22,28,34 repeat it).

RULE: join everything on Num. NEVER parse a number out of a Sample name:
  e.g. barcode 'ALHT WTD_102' (Num 4) has 16S Sample 'ALHT-Deer-10' -- the
  name-number (10) is NOT the barcode-number (102). Names are unreliable; Num is not.

Each marker Sample column is 879/879 filled = the intended name for every specimen,
whether or not that specimen was actually sequenced on that marker. So a name here
does not guarantee a feature-table sample exists; that is checked separately against
each marker's table.

OUTPUT: results/metadata/num_crosswalk.tsv  (durable), columns:
  Num, Barcode, Num_Ticks, Host, Town, State,
  s16S, s18S, scytb, sITS   (raw Sample names as written in the sheet)

Also provides normalize_* helpers to map a raw Sample name to the format used in
each marker's exported feature table (strip marker suffix, lowercase where needed).
Extend the normalizers as feature-table mismatches are found.

Usage:
    python scripts/build_crosswalk.py
"""
import os
import pandas as pd
from openpyxl import load_workbook

SRC = "metadata/tick_full_metadata.xlsx"
SHEET = "All Run Data - w Metadata"
OUTDIR = "results/metadata"
OUT = os.path.join(OUTDIR, "num_crosswalk.tsv")

# column indices in the positional sheet (0-based)
C = dict(num=1, barcode=0, nticks=4, host=6, town=10, state=11,
         s16S=14, s18S=20, scytb=26, sITS=32)


def load_rows():
    wb = load_workbook(SRC, read_only=True, data_only=True)
    ws = wb[SHEET]
    rows = list(ws.iter_rows(values_only=True))
    return rows[2:]  # skip the two header rows (block header + first field header row)


def clean(v):
    return None if v is None else str(v).strip()


def build():
    data = load_rows()
    recs = []
    for r in data:
        def g(k):
            i = C[k]
            return clean(r[i]) if i < len(r) else None
        num = g("num")
        if num in (None, ""):
            continue
        host = g("host")
        if host == "Deer":          # normalize the Deer/WTD split flagged earlier
            host = "WTD"
        state = g("state")
        if state:
            state = state.strip()   # 'NH ' whitespace seen in raw data
        recs.append(dict(
            Num=num, Barcode=g("barcode"), Num_Ticks=g("nticks"),
            Host=host, Town=g("town"), State=state,
            s16S=g("s16S"), s18S=g("s18S"), scytb=g("scytb"), sITS=g("sITS"),
        ))
    df = pd.DataFrame(recs)
    os.makedirs(OUTDIR, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    print("wrote %s  (%d specimens)" % (OUT, len(df)))
    print("host counts:\n" + df["Host"].value_counts(dropna=False).to_string())
    print("\nstate counts:\n" + df["State"].value_counts(dropna=False).to_string())
    return df


# --- normalizers: raw sheet Sample name -> feature-table ID format ---
# These are best-effort and MUST be validated against each marker's actual table.
# The sheet names often lack the run-suffix that the feature-table IDs carry
# (e.g. sheet '18S' name 'ALHT-Deer-02' vs table 'ALHT-Deer-02.240430_USGSTICKS-2023-2024'),
# so matching is usually done by prefix, not exact string. See match_to_table().

def strip_marker_suffix(name):
    if name is None:
        return None
    s = str(name)
    for suf in ("-16s", "-16S", "-cytoB", "-cytopB", "-cytob", "-ITS", "-ITS1-2"):
        if s.endswith(suf):
            s = s[: -len(suf)]
    return s


def match_to_table(sheet_name, table_ids):
    """Find the feature-table ID whose prefix matches this sheet Sample name.
    table_ids: iterable of real column IDs from the marker's exported table.
    Returns the matching ID or None. Prefix match because table IDs append
    '.<run>' to the sheet name."""
    if sheet_name is None:
        return None
    base = strip_marker_suffix(sheet_name)
    # exact first
    for t in table_ids:
        if t == base or t == sheet_name:
            return t
    # then prefix: table id starts with base + '.'
    for t in table_ids:
        if str(t).startswith(base + "."):
            return t
    return None


if __name__ == "__main__":
    build()
