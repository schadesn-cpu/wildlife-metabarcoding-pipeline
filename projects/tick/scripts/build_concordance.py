#!/usr/bin/env python3
"""
build_concordance.py  —  run ON RON, from the project root:
    cd ~/meedlab/tick_project
    python3 scripts/build_concordance.py

WHY THIS EXISTS
---------------
The qPCR-vs-metabarcoding concordance existed ONLY as a PNG (78 / 90 / 86%, n=93
as then reported). The values were later recomputed as 79.5 / 90.4 / 85.5% on
n=83, and that n was wrong: a local key normaliser stripped hyphens and lost the
six 'ALHT WTD_nn' specimens. On the shared lookup the cohort is n=89, being the
93 sequenced single ticks minus the 4 whose only qPCR rows are pools.
No script, no output table, nothing on disk that could regenerate it. It is the
entire evidence base for Goal 3, so it needed rebuilding from source rather than
re-running.

WHAT IT COMPARES
----------------
Per specimen screened by BOTH methods, a 2x2 of qPCR call vs metabarcoding call:

    qPCR target                metabarcoding column        note
    a-phagocytophilum   <->    16S_Anaplasma               GENUS vs species
    b-burgdorferi       <->    16S_Borreliella             GENUS vs species
    b-microti           <->    18S_Babesia_microti         species vs species

⚠️ TWO OF THE THREE PAIRS ARE ASYMMETRIC. 16S V4-V5 resolves Borreliella to genus
only (it cannot separate B. burgdorferi from B. bissettiae — see TICK_GROUND_TRUTH
section 1.2), and Anaplasma likewise. So a "metabarcoding-positive, qPCR-negative"
call may be a non-target congener rather than a disagreement. This MUST be stated
in the caption; it is a real limit on what concordance means here.

Reports raw agreement AND Cohen's kappa. Raw agreement is inflated by the large
both-negative cell; kappa corrects for agreement expected by chance and is the more
honest headline number.

Reads:
    results/metadata/library_lookup.tsv   (via tick_lookup; the ONLY name resolver)
    results/single_tick_v3/single_tick_v3_full.tsv
    results/analysis/qpcr_screening_record_corrected.tsv
Writes:
    results/concordance/concordance_2x2.tsv
    results/concordance/concordance_report.txt
    results/figures/fig4_qpcr_concordance.png/.pdf
"""
import os
import re
import sys
import pandas as pd

ST   = "results/single_tick_v3/single_tick_v3_full.tsv"
WB   = "metadata/tick_full_metadata.xlsx"
QPCR = "results/analysis/qpcr_screening_record_corrected.tsv"
ODIR = "results/concordance"
FDIR = "results/figures"

PAIRS = [
    ("Anaplasma /\nA. phagocytophilum", "a-phagocytophilum", "16S_Anaplasma",       "genus vs species"),
    ("Babesia microti",                 "b-microti",         "18S_Babesia_microti", "species vs species"),
    ("Borreliella /\nB. burgdorferi",   "b-burgdorferi",     "16S_Borreliella",     "genus vs species"),
]

# Figure-only display labels. WRITING_STYLE_GUIDE.md requires italic genus and
# species names; the slash and the abbreviating period stay upright, which is
# why this is mathtext rather than fontstyle='italic' on the whole label.
# Keyed on the plain label so the TSV keeps plain text.
DISPLAY_LABELS = {
    "Anaplasma / A. phagocytophilum":
        "$\\mathit{Anaplasma}$ /\n$\\mathit{A.\\ phagocytophilum}$",
    "Babesia microti":
        "$\\mathit{Babesia\\ microti}$",
    "Borreliella / B. burgdorferi":
        "$\\mathit{Borreliella}$ /\n$\\mathit{B.\\ burgdorferi}$",
}

rep = []
def say(m=""):
    print(m); rep.append(str(m))

# Name resolution is NOT done here. It is done once, by
# build_4marker_crosswalk.py, and read from results/metadata/library_lookup.tsv.
# The previous norm() in this slot stripped hyphens as well as spaces, so
# 'ALHT WTD_09' and the qPCR sheet's 'alht-deer-09' became two different keys
# and six single ticks silently failed to pair. Do not reintroduce a local
# normaliser; if a name does not resolve, fix the crosswalk.
import tick_lookup

LOOKUP = tick_lookup.load()


def specimen_of(value):
    """I resolve one accessioning ID, barcode or qPCR key to a specimen."""
    return tick_lookup.resolve(LOOKUP, "specimen", value)

os.makedirs(ODIR, exist_ok=True)
os.makedirs(FDIR, exist_ok=True)

# ---- inputs ----------------------------------------------------------------
st = pd.read_csv(ST, sep="\t")
say(f"single-tick set: {len(st)} specimens ({ST})")

q = pd.read_csv(QPCR, sep="\t", dtype=str)
q["n"] = pd.to_numeric(q["subsample_count"], errors="coerce")
q1 = q[q["n"] == 1].copy()
say(f"qPCR single-tick rows: {len(q1)} of {len(q)}")

# qPCR has BOTH a Barcode column (sometimes blank) and host_id. Either one
# resolves through the lookup, so both are tried and the result is keyed on the
# SPECIMEN rather than on whichever string that row happened to carry.
qmap = {}
qclash = []
for _, r in q1.iterrows():
    hit = specimen_of(r["Barcode"]) or specimen_of(r["host_id"])
    if hit is None:
        continue
    if hit in qmap:
        qclash.append(hit)
    else:
        qmap[hit] = r

say(f"qPCR single-subsample rows resolving to a specimen: {len(qmap)}")
if qclash:
    say(f"  WARN {len(qclash)} specimen(s) have more than one single-subsample "
        f"qPCR row; the first was kept: {sorted(set(qclash))}")

st["_specimen"] = st["Barcode"].map(specimen_of)
_unresolved = st.loc[st["_specimen"].isna(), "Barcode"].astype(str).tolist()
if _unresolved:
    say(f"  WARN {len(_unresolved)} sequencing specimen(s) are not in the "
        f"lookup at all: {_unresolved}")
    say("       That is a crosswalk problem, not a join problem. Rebuild it.")
st["_q"] = st["_specimen"].map(lambda s: qmap.get(s) if s is not None else None)

paired = st[st["_q"].notna()].copy()
say(f"\nspecimens screened by BOTH methods: {len(paired)}")
unpaired = sorted(st.loc[st["_q"].isna(), "Barcode"].astype(str))
say(f"in the sequencing set but with NO single-subsample qPCR row: {len(unpaired)}")
say("  These specimens ARE in the qPCR sheet. Their only rows are pools, and a")
say("  pooled qPCR result cannot be compared against a single tick. This is the")
say("  reason the cohort is smaller than the sequencing set, and it is the only")
say("  reason: every single tick resolves against the qPCR sheet.")
for b in unpaired[:25]:
    say(f"      {b}")
if len(unpaired) > 25:
    say(f"      ... and {len(unpaired)-25} more")

if len(paired) == 0:
    sys.exit("ERROR: nothing joined. Check the key normalisation before trusting anything.")

# ---- 2x2 per pathogen ------------------------------------------------------
def kappa(a, b, c, d):
    """Cohen's kappa from a 2x2: a=both+, b=qPCR only, c=seq only, d=both-."""
    n = a + b + c + d
    if n == 0:
        return float("nan")
    po = (a + d) / n
    pe = ((a + b) * (a + c) + (c + d) * (b + d)) / (n * n)
    return (po - pe) / (1 - pe) if pe != 1 else float("nan")

rows = []
say("\n=== concordance, per pathogen ===")
for label, qcol, scol, note in PAIRS:
    if qcol not in paired.columns and qcol not in q1.columns:
        say(f"  {label}: qPCR column '{qcol}' missing, skipped"); continue
    if scol not in paired.columns:
        say(f"  {label}: metabarcoding column '{scol}' missing, skipped"); continue

    qv = paired["_q"].map(lambda r: pd.to_numeric(r[qcol], errors="coerce")).fillna(0) > 0
    sv = pd.to_numeric(paired[scol], errors="coerce").fillna(0) > 0

    a = int((qv & sv).sum())        # both positive
    b = int((qv & ~sv).sum())       # qPCR only
    c = int((~qv & sv).sum())       # metabarcoding only
    d = int((~qv & ~sv).sum())      # both negative
    n = a + b + c + d
    agree = 100 * (a + d) / n
    k = kappa(a, b, c, d)

    say(f"\n  {label.replace(chr(10),' ')}   [{note}]")
    say(f"    both positive        {a:>3}")
    say(f"    qPCR only            {b:>3}")
    say(f"    metabarcoding only   {c:>3}")
    say(f"    both negative        {d:>3}")
    say(f"    raw agreement        {agree:.1f}%   (n = {n})")
    say(f"    Cohen's kappa        {k:.2f}")

    rows.append({"pathogen": label.replace(chr(10), " "), "comparison": note,
                 "both_pos": a, "qpcr_only": b, "seq_only": c, "both_neg": d,
                 "n": n, "agreement_pct": round(agree, 1), "kappa": round(k, 3)})

res = pd.DataFrame(rows)
out = os.path.join(ODIR, "concordance_2x2.tsv")
res.to_csv(out, sep="\t", index=False)
say(f"\nwrote {out}")

say("\nDIRECTION OF DISAGREEMENT — read this before writing any prose:")
for r in rows:
    if r["qpcr_only"] > r["seq_only"]:
        say(f"  {r['pathogen']}: qPCR finds more ({r['qpcr_only']} vs {r['seq_only']}) "
            f"-> sequencing is missing detections, consistent with read-depth limits")
    elif r["seq_only"] > r["qpcr_only"]:
        say(f"  {r['pathogen']}: sequencing finds more ({r['seq_only']} vs {r['qpcr_only']}) "
            f"-> may be non-target congeners, since 16S resolves genus only")
    else:
        say(f"  {r['pathogen']}: symmetric ({r['qpcr_only']} each way)")

# ---- figure ----------------------------------------------------------------
# Top-level import on purpose: if plotting fails the run has failed, and I would
# rather crash than write a report that reads like a successful run.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if True:   # keeps the plotting block indented
    GREEN, BLUE, ORANGE, GREY = "#009E73", "#0072B2", "#E69F00", "#CCCCCC"
    fig, ax = plt.subplots(figsize=(10, 6))
    labels = [DISPLAY_LABELS.get(r["pathogen"], r["pathogen"])
              for r in rows]
    _unstyled = [r["pathogen"] for r in rows
                 if r["pathogen"] not in DISPLAY_LABELS]
    if _unstyled:
        say("  WARN no italic display label for: "
            + ", ".join(_unstyled)
            + ". Add it to DISPLAY_LABELS; it will plot in roman.")
    xs = range(len(rows))
    bp  = [r["both_pos"] for r in rows]
    qo  = [r["qpcr_only"] for r in rows]
    so  = [r["seq_only"] for r in rows]
    bn  = [r["both_neg"] for r in rows]

    ax.bar(xs, bp, color=GREEN, label="Both positive (concordant)")
    ax.bar(xs, qo, bottom=bp, color=BLUE, label="qPCR positive only")
    ax.bar(xs, so, bottom=[i+j for i, j in zip(bp, qo)], color=ORANGE,
           label="Metabarcoding positive only")
    ax.bar(xs, bn, bottom=[i+j+k for i, j, k in zip(bp, qo, so)], color=GREY,
           label="Both negative (concordant)")

    for i, r in enumerate(rows):
        ax.text(i, r["n"] * 1.03, f"{r['agreement_pct']:.1f}% agree\n\u03ba = {r['kappa']:.2f}",
                ha="center", fontweight="bold", fontsize=11)
        y = 0
        for v, col in ((r["both_pos"], "white"), (r["qpcr_only"], "white"),
                       (r["seq_only"], "white"), (r["both_neg"], "#333333")):
            if v > max(rows[0]["n"] * 0.03, 1):
                ax.text(i, y + v / 2, str(v), ha="center", va="center",
                        color=col, fontweight="bold", fontsize=10)
            y += v

    ax.set_xticks(list(xs)); ax.set_xticklabels(labels)
    ax.set_ylabel(f"Number of ticks (n = {rows[0]['n']} screened by both methods)")
    ax.set_title("qPCR vs. metabarcoding pathogen detection, paired per-tick concordance",
                 fontweight="bold")
    ax.set_ylim(0, rows[0]["n"] * 1.18)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2, frameon=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = os.path.join(FDIR, f"fig4_qpcr_concordance.{ext}")
        fig.savefig(p, dpi=300, bbox_inches="tight")
        say(f"wrote {p}")
    plt.close(fig)

say("\nCAPTION MUST STATE: two of the three comparisons are genus (16S) against")
say("species (qPCR). A metabarcoding-positive, qPCR-negative call may be a")
say("non-target congener rather than a disagreement between methods.")

with open(os.path.join(ODIR, "concordance_report.txt"), "w") as fh:
    fh.write("\n".join(rep))
say(f"\nreport -> {ODIR}/concordance_report.txt")
