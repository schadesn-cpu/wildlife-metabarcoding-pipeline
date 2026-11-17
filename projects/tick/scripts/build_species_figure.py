#!/usr/bin/env python3
"""
build_species_figure.py  —  run ON RON, from the project root:
    cd ~/meedlab/tick_project
    python3 scripts/build_species_figure.py

Tick species by collection host, at SPECIMEN level, from the corrected 18S taxonomy.

WHY REBUILT
-----------
The old figure was at TICK level (n=495, 18 Dermacentor) and predates every
correction since. This reads `tick_genus_per_sample_corrected.tsv` (written by the
18S apply, 2026-08-06) and joins to `crosswalk_4marker.tsv` (rebuilt 2026-08-07),
so it agrees with the locked numbers by construction rather than by luck.

THE RULE (revised 2026-08-10)
-----------------------------
1. READS ARE SUMMED ACROSS A SPECIMEN'S LIBRARIES before thresholding. A specimen
   split into two tubes must not fail a threshold a single tube would pass.
   TB_90024117 is why: 65 + 36 Dermacentor reads across two libraries, so neither
   library cleared 100 on its own and a genuinely positive moose was being dropped.
2. PRESENCE, NOT DOMINANCE. A specimen is positive for a genus if it carries
   >= MIN_READS of it, whether or not another genus is more abundant. These are
   POOLED samples by design; a pool off one animal can hold both species, and
   "this animal carried winter ticks" is the biologically meaningful statement,
   not "winter ticks outnumbered deer ticks in this tube". Categories therefore
   OVERLAP and the panels are not stacked.
3. SPECIMENS WITH ZERO TICK READS are reported as uninformative, not as negatives.
   TB_90024100 is why: its 18S library returned no tick reads at all, which is a
   failed library, not a moose without winter ticks.

WHAT TO EXPECT
--------------
The counts are NOT restated here. A hardcoded expectation in a docstring goes
stale silently, which is what happened to this one: it said 18 Dermacentor
specimens with 9 moose while the script had been printing 20 and 11 for weeks,
and nothing noticed because prose cannot fail. The same failure is recorded for
GOLDEN in build_single_tick_v2.py.

The current counts are in results/analysis/tick_species_report.txt, written by
this script on every run.

The result that MATTERS, and that this script does assert rather than describe,
is ZERO Dermacentor from deer. That is what justifies restricting the diversity
statistics to the deer subset: it controls host and tick species at once. The
script prints a WARNING instead if it ever stops being true.

2026-09-18: the library join moved from a private resolver to the shared
library lookup, which added 185 specimens the old join was silently dropping.
Any count quoted from a run before that date is on two-thirds of the data.

Reads:
    results/18S/all/taxonomy_verified/tick_genus_per_sample_corrected.tsv
    results/metadata/crosswalk_4marker.tsv
Writes:
    results/figures/fig_tick_species_by_host.png/.pdf
    results/analysis/tick_species_by_host.tsv
    results/analysis/tick_species_report.txt
"""
import os
import sys
import pandas as pd

TG   = "results/18S/all/taxonomy_verified/tick_genus_per_sample_corrected.tsv"
CW   = "results/metadata/crosswalk_4marker.tsv"
FDIR = "results/figures"
ODIR = "results/analysis"

# ── THE FLOOR IS DERIVED FROM THE BLANKS, NOT CHOSEN ────────────────────────
# Extraction blanks and controls carry a measurable Dermacentor signal from index
# bleed: XB-B11 = 15 reads, XB-B2 = 16, XB-B1 = 28, POS-H6 = 25 (2026-08-10).
# Every genuine Dermacentor specimen is 3,489-197,719 reads. Three orders of
# magnitude separate them, so the cutoff does not need to be argued for.
#
# MIN_READS is set an order of magnitude above the highest blank (28 -> 280 -> 300).
# A relative floor is applied as well, because absolute counts mean different
# things at different sequencing depths and the metabarcoding literature works in
# relative abundance rather than raw counts.
#
# Verified threshold-independent 2026-08-10. Across five inclusion rules
# (>=1 read, >=10, >=100 absolute; >=1%, >=10% relative) the result did not move:
#   moose 11/11 under EVERY rule
#   deer  0/222-273 (the single deer that appears at the loosest absolute cutoffs
#         is ALHT-Deer-64: 53,102 Ixodes vs 15 Dermacentor = 0.03%, squarely
#         inside the blank range, i.e. bleed)
#   other wildlife 9/14 under four of five
BLANK_MAX = 28           # highest Dermacentor count seen in a blank/control
MIN_READS = 50           # ~2x the highest blank
MIN_FRAC  = 0.01         # AND >= 1% of that specimen's tick reads
#
# Why 50 and not 10x the blank: the RELATIVE floor is what separates bleed from
# signal, and the two failure modes look completely different.
#   ALHT-Deer-64   53,102 Ixodes / 15 Dermacentor = 0.03%  -> bleed into a large
#                  library. Fails both floors.
#   TB_90024117         0 Ixodes / 101 Dermacentor = 100%  -> a real winter-tick
#                  moose sequenced shallowly, with NO Ixodes at all. Bleed would
#                  have carried the dominant genus along with it; this did not.
# So the rule is: exceed the highest blank, and be a real fraction of the tick DNA
# in that specimen. A 10x absolute floor alone would discard TB_90024117.
GENERA = ["Ixodes", "Dermacentor", "Amblyomma"]

# Host display order and grouping. Opossum and porcupine are collapsed for the
# figure only; the underlying table keeps them separate.
# Panel B groups by the crosswalk's host_bin (opossum + porcupine are already
# "Other wildlife" there). The written table keeps host_spp separate.

rep = []
def say(m=""):
    print(m); rep.append(str(m))

os.makedirs(FDIR, exist_ok=True)
os.makedirs(ODIR, exist_ok=True)

for f in (TG, CW):
    if not os.path.exists(f):
        sys.exit(f"ERROR: {f} not found")

tg = pd.read_csv(TG, sep="\t")
cw = pd.read_csv(CW, sep="\t")
say(f"tick-genus table: {len(tg)} libraries")
say(f"crosswalk:        {len(cw)} specimens")

# ---- library -> specimen, via name_18S -------------------------------------
# Name resolution is NOT done here. It is done once, by
# build_4marker_crosswalk.py, and read from results/metadata/library_lookup.tsv.
#
# What used to sit here built a private lib2spec from name_18S by splitting on
# semicolons and swapping a literal "-18s" suffix, then discarded every miss as
# "(controls and excluded specimens — expected)". Measured 2026-09-18: it
# resolved 498 of 907 libraries where the lookup resolves 853, missing 185
# distinct specimens, 33 of them single ticks. The misses were bare accessions
# (EQ02440612), replicate suffixes (EQ02441291-2 .. -10) and the ALHT alias
# (ALHT-Deer-103 -> ALHT WTD_103): the three shapes in TICK_GROUND_TRUTH
# section 3.
#
# Do not reintroduce a local resolver. If a library does not resolve, that is a
# crosswalk problem and it is now named below rather than assumed away.
sys.path.insert(0, "scripts")
import tick_lookup

LOOKUP = tick_lookup.load()

tg["lib"] = tg["sample_id"].astype(str).str.replace(r"\..*$", "", regex=True)
tg["specimen"] = [
    tick_lookup.resolve(LOOKUP, "18S", lib)
    or tick_lookup.resolve(LOOKUP, "specimen", lib)
    for lib in tg["lib"]]

unmatched = tg[tg["specimen"].isna()]
say(f"\nlibraries with no crosswalk specimen: {len(unmatched)} of {len(tg)}")

# No claim about what these are. The list is written out so it can be checked.
_ctrl_markers = ("NTC", "NEG", "BLK", "POS", "XB", "EXTRACTIONBLANK",
                 "PC-WELL", "BLANK", "CONTROL")
if len(unmatched):
    _u = unmatched[["sample_id", "lib"]].copy()
    _u["looks_like_control"] = _u["lib"].astype(str).str.upper().apply(
        lambda s: any(m in s for m in _ctrl_markers))
    _u.to_csv(f"{ODIR}/species_figure_unmatched.tsv", sep="\t", index=False)
    _n_ctrl = int(_u["looks_like_control"].sum())
    say(f"  {_n_ctrl} carry a control prefix; {len(_u) - _n_ctrl} do not")
    say(f"  named in {ODIR}/species_figure_unmatched.tsv")
    if len(_u) - _n_ctrl:
        say("  The non-control ones are libraries with no crosswalk row. That is")
        say("  a crosswalk gap, not an expected loss. Check them before quoting")
        say("  any number from this figure.")

tg = tg[tg["specimen"].notna()].copy()
say(f"libraries joined to a specimen:  {len(tg)}")

# ---- specimen-level genus calls (summed across libraries) -------------------
present = [g for g in GENERA if g in tg.columns]
for g in present:
    tg[g] = pd.to_numeric(tg[g], errors="coerce").fillna(0)
agg = tg.groupby("specimen")[present].sum()
say(f"\nsummed to {len(agg)} specimens")

meta = cw.set_index("sample_id")
rows = []
for spec, r in agg.iterrows():
    tot = int(r[present].sum())
    host = meta.loc[spec, "host_bin"] if spec in meta.index else "Unknown"
    hsp  = meta.loc[spec, "host_spp"] if spec in meta.index else "Unknown"
    st   = meta.loc[spec, "state"] if spec in meta.index else ""
    # "informative" = the specimen produced ANY tick DNA. A deer with 40 reads that
    # are all Ixodes IS informative about that deer; discarding it would shrink the
    # denominator for no reason. Only specimens with zero tick reads are excluded.
    rec = {"specimen": spec, "host_bin": host, "host_spp": hsp, "state": st,
           "total_tick_reads": tot, "informative": tot > 0}
    for g in present:
        rec[g + "_reads"] = int(r[g])
        rec[g + "_frac"] = round(r[g] / tot, 4) if tot else 0.0
        # BOTH conditions: clears the blank-derived absolute floor AND is a real
        # fraction of this specimen's tick DNA. Either alone lets bleed through.
        rec[g + "_pos"] = bool(r[g] >= MIN_READS and (r[g] / tot if tot else 0) >= MIN_FRAC)
    rows.append(rec)

df = pd.DataFrame(rows)

say(f"\nfloor: >= {MIN_READS} reads AND >= {100*MIN_FRAC:.0f}% of the specimen's tick "
    f"reads (blanks carried up to {BLANK_MAX} Dermacentor reads)")
uninf = df[~df["informative"]]
say(f"uninformative (ZERO tick reads, nothing to call): {len(uninf)}")
for h, n in uninf["host_bin"].value_counts().items():
    say(f"      {h:<18} {n}")
df = df[df["informative"]].copy()
say(f"informative specimens: {len(df)}")

for g in present:
    say(f"  {g}-positive: {int(df[g+'_pos'].sum())}")
mixed = df[df.get("Ixodes_pos", False) & df.get("Dermacentor_pos", False)]
if len(mixed):
    say(f"\n  MIXED (both genera >= {MIN_READS}): {len(mixed)} specimen(s) — these are")
    say(f"  pools containing both species. Counted in BOTH categories, not forced to one.")
    for _, m in mixed.iterrows():
        say(f"      {m['specimen']:<24} Ix={m['Ixodes_reads']:>8}  Derm={m['Dermacentor_reads']:>8}")
df.to_csv(f"{ODIR}/tick_species_by_host.tsv", sep="\t", index=False)
say(f"\nwrote {ODIR}/tick_species_by_host.tsv  ({len(df)} specimens with a genus call)")

say("\n=== overall (categories overlap: a pool can hold both) ===")
for g in present:
    n = int(df[g + "_pos"].sum())
    say(f"  {g:<14} {n:>4}  ({100*n/len(df):5.1f}% of {len(df)} informative specimens)")

say("\n=== by collection host: n positive / n informative (%) ===")
say(f"  {'host':<20} {'n':>4}  " + "  ".join(f"{g:>22}" for g in present))
for h, sub_ in df.groupby("host_bin"):
    cells = []
    for g in present:
        k = int(sub_[g + "_pos"].sum())
        cells.append(f"{k:>3}/{len(sub_):<3} ({100*k/len(sub_):5.1f}%)".rjust(22))
    say(f"  {h:<20} {len(sub_):>4}  " + "  ".join(cells))

say("\n=== Dermacentor detail (the finding) ===")
d = df[df["Dermacentor_pos"]]
say(f"  {len(d)} Dermacentor specimens")
for (h, s), n in d.groupby(["host_spp", "state"]).size().items():
    say(f"    {str(h):<22} {s}  {n}")
if (d["host_bin"] == "Deer").sum() == 0:
    say("  ZERO from deer -> the deer subset is confirmed single-species Ixodes,")
    say("  which is what justifies the WTD-only restriction on the diversity stats.")
else:
    say(f"  WARNING: {(d['host_bin']=='Deer').sum()} Dermacentor from DEER. "
        f"This contradicts the locked result. Investigate before using.")

# ---- figure ----------------------------------------------------------------
# I import matplotlib at the top level rather than inside a try. If plotting is
# unavailable the run has failed and I want to know, not get a report file and a
# "figure skipped" line that looks like success.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if True:   # keeps the plotting block indented; see note above
    BLUE, ORANGE, GREEN = "#0072B2", "#E69F00", "#009E73"
    COLOR = {"Ixodes": BLUE, "Dermacentor": ORANGE, "Amblyomma": GREEN}

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.5),
                                   gridspec_kw={"width_ratios": [1, 1.5]})

    # A: overall (overlapping categories)
    vc = pd.Series({g: int(df[g + "_pos"].sum()) for g in present})
    vc = vc[vc > 0].sort_values(ascending=False)
    # Mathtext, not fontstyle='italic': the abbreviating period and the
    # spacing are not part of the organism name. Matches fig2 and fig3.
    GENUS_LABEL = {
        "Ixodes": "$\\mathit{I.\\ scapularis}$",
        "Dermacentor": "$\\mathit{Dermacentor}$",
        "Amblyomma": "$\\mathit{Amblyomma}$",
    }
    axA.bar([GENUS_LABEL.get(g, g) for g in vc.index], vc.values,
            color=[COLOR.get(g, "#999999") for g in vc.index])
    for i, (g, n) in enumerate(vc.items()):
        axA.text(i, n + max(vc.values) * 0.02, f"{n}\n{100*n/len(df):.1f}%",
                 ha="center", fontweight="bold")
    axA.set_ylabel(f"Specimens (n = {len(df)} informative)")
    # The bars do not sum to n, for two reasons, and a reader who subtracts
    # deserves to be told both rather than left with unexplained specimens.
    _pos_cols = [g + "_pos" for g in present]
    _n_nocall = int((~df[_pos_cols].any(axis=1)).sum())
    _n_both = int((df[_pos_cols].sum(axis=1) > 1).sum())
    axA.text(0.5, -0.20,
             "Categories overlap: %d specimen(s) positive for two genera.\n"
             "%d of %d informative specimens cleared no genus floor."
             % (_n_both, _n_nocall, len(df)),
             transform=axA.transAxes, ha="center", va="top",
             fontsize=8, style="italic", color="#666")
    say(f"\nspecimens positive for BOTH genera: {_n_both}")
    say(f"informative specimens clearing NO genus floor: {_n_nocall} "
        f"of {len(df)} (they have tick reads but below the floor)")
    axA.set_title("A. Tick genus overall", loc="left", fontweight="bold")
    axA.set_ylim(0, max(vc.values) * 1.18)
    for sp in ("top", "right"):
        axA.spines[sp].set_visible(False)

    # B: by host, collapsed
    # Use the crosswalk's own host_bin rather than string-matching host_spp:
    # host_spp carries variants like "NA Porcupine" that a literal map misses.
    LABEL = {"Deer": "White-tailed deer", "Moose": "Moose", "Bear": "Black bear",
             "Other wildlife": "Other wildlife", "Unknown": "Host not recorded"}
    df["host_fig"] = df["host_bin"].map(lambda h: LABEL.get(h, h))
    # B: PROPORTION positive within each host. Counts are unreadable because deer
    # outnumber every other host by two orders of magnitude; proportion is the
    # biology. n per host is printed on the axis so the reader can weight it.
    gh = df.groupby("host_fig")
    hosts = sorted(gh.groups, key=lambda h: -len(gh.get_group(h)))
    ypos = range(len(hosts))
    w = 0.38
    for j, g in enumerate([c for c in ["Dermacentor", "Ixodes"] if c in present]):
        pct = [100 * gh.get_group(h)[g + "_pos"].mean() for h in hosts]
        off = (j - 0.5) * w
        axB.barh([y + off for y in ypos], pct, height=w,
                 color=COLOR.get(g, "#999999"),
                 label=GENUS_LABEL.get(g, g))
        for y, p in zip(ypos, pct):
            if p > 0:
                axB.text(p + 1.5, y + off, f"{p:.0f}%", va="center", fontsize=9)
    axB.set_yticks(list(ypos))
    axB.set_yticklabels([f"{h}\n(n = {len(gh.get_group(h))})" for h in hosts])
    axB.invert_yaxis()
    axB.set_xlim(0, 118)
    axB.set_xlabel("% of informative specimens positive")
    axB.set_title("B. Tick genus by collection host", loc="left", fontweight="bold")
    # legend outside the axes: at 100% the bars run the full width and any
    # in-axes placement collides with a bar or its label
    axB.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    for sp in ("top", "right"):
        axB.spines[sp].set_visible(False)

    fig.suptitle("Tick genus confirmed from 18S metabarcoding, not recorded at collection",
                 fontweight="bold", y=1.00)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = f"{FDIR}/fig_tick_species_by_host.{ext}"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        say(f"wrote {p}")
    plt.close(fig)

say("\nCAPTION SHOULD STATE:")
say("  - tick genus was NOT recorded at collection; it was recovered from 18S")
say(f"  - reads summed across a specimen's libraries; positive requires >= {MIN_READS}")
say(f"    reads AND >= {100*MIN_FRAC:.0f}% of that specimen's tick reads. The floor is")
say(f"    derived from the blanks, which carried up to {BLANK_MAX} Dermacentor reads,")
say("    against 3,489-197,719 in genuine positives")
_n_moose = int((df["host_bin"] == "Moose").sum())
_n_deer = int((df["host_bin"] == "Deer").sum())
say("  - the result is threshold-independent, rerun 2026-09-18 on the full")
say("    library join: across five rules (>=1, >=10, >=100 reads; >=1%, >=10%")
say(f"    relative) moose stayed {_n_moose}/{_n_moose} and deer 0-1 of {_n_deer}.")
say("    The single deer detection is ALHT WTD_64, 15 Dermacentor reads against")
say("    53,119 Ixodes (0.03%), inside the blank ceiling of 28 reads.")
say("    Source: results/analysis/threshold_independence.tsv")
say("  - the bars do not sum to n. Categories OVERLAP (a pool off one animal")
say("    can hold both genera, so a specimen may be counted twice) and some")
say("    informative specimens clear neither genus floor. Both counts are")
say("    printed above and drawn under panel A.")
say("  - panel B is PROPORTION within host, with n per host on the axis, because")
say("    raw counts are unreadable when deer outnumber every other host 40-fold")
say("  - Amblyomma is absent: its single detection came from an out-of-region")
say("    sample excluded from the study")

with open(f"{ODIR}/tick_species_report.txt", "w") as fh:
    fh.write("\n".join(rep))
say(f"\nreport -> {ODIR}/tick_species_report.txt")
