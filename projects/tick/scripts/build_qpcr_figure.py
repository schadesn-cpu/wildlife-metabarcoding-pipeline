#!/usr/bin/env python3
"""
build_qpcr_figure.py

qPCR pathogen prevalence and co-infection in individually screened ticks.

Run:
    cd ~/meedlab/tick_project
    python3 scripts/build_qpcr_figure.py

Why this exists: the qPCR prevalences are the Goal 1 headline numbers and nothing
in scripts/ produced them. I reproduced them by hand on 2026-08-07 and they were
correct, but a number nobody can regenerate is a number nobody can check.

WHICH TICKS
-----------
CORRECTED 2026-09-14. I previously used the 423 rows with subsample_count == 1
and described them as one tick per assay. That was wrong. Kellen Sawyer confirmed
in person that a subsample is roughly 3 ticks taken from the tube and slurried
together, so subsample_count == 1 means one slurry was tested, not one tick. 265
of those 423 came from tubes holding more than one tick.

I now require BOTH subsample_count == 1 AND the tube to have held exactly one
tick, which gives 100 genuine per-tick assays. The 423-row figures are pool-level
positivity and are printed for comparison but not plotted as per-tick prevalence.
That is a per-tick prevalence and it is the quantity comparable in kind to the
metabarcoding figure. I deliberately do NOT plot the all-sample figures (569
samples, A. phagocytophilum 55.9%): those are pool-level positivity, a different
question, and mixing the two in one figure is how denominators get confused.

⚠️ THIS IS A DIFFERENT POPULATION FROM THE SEQUENCING FIGURE. The qPCR
per-tick cohort and the 93 sequenced single ticks overlap in 89, which is the
concordance cohort, but they are not the same set. Panel B here and panel B of
the metabarcoding figure are NOT the same ticks, and the qPCR targets are
species while the metabarcoding co-infection counts genera. I expect qPCR
co-infection to run higher, because qPCR finds A. phagocytophilum and B. microti
at pool level at 50.1% and 30.0% where sequencing finds 22.6% and 4.3% per tick.
That difference is a result worth writing about, not an inconsistency to smooth
over.

I show all seven panel targets including the zeros. E. muris at 0/423 and
B. mayonii at 1/423 document what the panel covers, which is exactly the point
when the argument is that qPCR catches things metabarcoding cannot.

Reads:
    metadata/tick_full_metadata.xlsx, sheet 'qPCR Results'
Writes:
    results/figures/fig_qpcr_prevalence_coinfection.{png,pdf}
    results/analysis/qpcr_prevalence.tsv
    results/analysis/qpcr_report.txt
"""
import os
import sys
from math import sqrt

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORKBOOK = "metadata/tick_full_metadata.xlsx"
SHEET = "qPCR Results"
FIGURE_DIR = "results/figures"
ANALYSIS_DIR = "results/analysis"

# qPCR column -> display label. Order is the panel as run, not by prevalence;
# I re-sort for the plot but keep this as the record of what was screened.
TARGET_LABELS = {
    "a-phagocytophilum": r"$\mathit{Anaplasma\ phagocytophilum}$",
    "b-burgdorferi":     r"$\mathit{Borreliella\ burgdorferi}$",
    "b-microti":         r"$\mathit{Babesia\ microti}$",
    "b-miyamotoi":       r"$\mathit{Borrelia\ miyamotoi}$",
    "powassan-virus":    "Powassan virus",
    "b-mayonii":         r"$\mathit{Borreliella\ mayonii}$",
    "e-muris":           r"$\mathit{Ehrlichia\ muris}$",
}

BLUE = "#0072B2"
GREEN = "#009E73"

report_lines = []
def say(message=""):
    print(message)
    report_lines.append(str(message))

def wilson_interval(n_positive, n_total):
    """95% Wilson interval, in percent. Correct for proportions at small n, which
    matters here: Powassan is 1 of 100 and a normal-approximation interval on
    that is nonsense."""
    if n_total == 0:
        return 0.0, 0.0
    z = 1.96
    p = n_positive / n_total
    denom = 1 + z * z / n_total
    centre = (p + z * z / (2 * n_total)) / denom
    half = z * sqrt(p * (1 - p) / n_total + z * z / (4 * n_total * n_total)) / denom
    return 100 * max(0.0, centre - half), 100 * min(1.0, centre + half)

if not os.path.exists(WORKBOOK):
    sys.exit(f"ERROR: {WORKBOOK} not found. Run me from the project root.")

os.makedirs(FIGURE_DIR, exist_ok=True)
os.makedirs(ANALYSIS_DIR, exist_ok=True)

qpcr_all = pd.read_excel(WORKBOOK, sheet_name=SHEET, dtype=str)
missing_columns = [c for c in TARGET_LABELS if c not in qpcr_all.columns]
if missing_columns:
    sys.exit(f"ERROR: {SHEET} is missing expected columns: {missing_columns}")

qpcr_all["subsamples"] = pd.to_numeric(qpcr_all["subsample_count"], errors="coerce")
n_unparseable = int(qpcr_all["subsamples"].isna().sum())
if n_unparseable:
    print(f"warning: {n_unparseable} rows have a non-numeric subsample_count "
          f"and are excluded", file=sys.stderr)

# A subsample is a slurry of roughly 3 ticks, so subsample_count == 1 is NOT a
# single tick. I additionally require the tube to have held exactly one tick,
# taken from the crosswalk, which is the only place num_ticks lives.
# Name resolution is NOT done here, and neither is tick counting. Both come
# from build_4marker_crosswalk.py, via results/metadata/library_lookup.tsv and
# results/metadata/crosswalk_4marker.tsv.
#
# What used to sit in this slot did two harmful things. It normalised names by
# stripping every non-alphanumeric character, so 'ALHT WTD_09' and the qPCR
# sheet's 'alht-deer-09' became two different keys. And it took tick counts
# from the RAW run sheet, falling back to the crosswalk only where the run
# sheet was blank. The run sheet's 'Num Ticks' is the column the ALHT/LHT
# barcode clobber corrupted: ALHT WTD_09 reads 18 ticks there against 1 in the
# registry and the crosswalk, because its run row carries the clobbered
# 'USGSTICKS-LHT-Deer-091' barcode. A genuine single tick was therefore read
# as an 18-tick pool, and the fallback never fired because the value was not
# blank, only wrong.
#
# The crosswalk is the authority. Do not reintroduce a run-sheet lookup.
import tick_lookup

LOOKUP = tick_lookup.load()
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
if not os.path.exists(CROSSWALK):
    sys.exit("ERROR: " + CROSSWALK + " not found. Build it first:\n"
             "    python3 scripts/build_4marker_crosswalk.py")

crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
for _needed in ("sample_id", "num_ticks", "is_single_tick"):
    if _needed not in crosswalk.columns:
        sys.exit("ERROR: column '" + _needed + "' missing from " + CROSSWALK)

tick_count_of = {}
single_tick_flag = {}
for _specimen, _n, _flag in zip(crosswalk["sample_id"], crosswalk["num_ticks"],
                                crosswalk["is_single_tick"]):
    tick_count_of[_specimen] = pd.to_numeric(pd.Series([_n]), errors="coerce")[0]
    single_tick_flag[_specimen] = str(_flag).strip().lower() == "true"


def specimen_of(row):
    """I resolve one qPCR row to a crosswalk specimen, or None if it is absent.

    None is a real answer. Nine rows in this sheet have no crosswalk specimen,
    two of them pathogen-positive, and they should stay visible rather than
    being forced into a denominator.
    """
    for column in ("Barcode", "host_id"):
        hit = tick_lookup.resolve(LOOKUP, "specimen", row.get(column))
        if hit is not None:
            return hit
    return None


qpcr_all["specimen"] = [specimen_of(row) for _, row in qpcr_all.iterrows()]
qpcr_all["num_ticks"] = qpcr_all["specimen"].map(
    lambda s: tick_count_of.get(s) if s is not None else None)
qpcr_all["is_single_tick"] = qpcr_all["specimen"].map(
    lambda s: bool(single_tick_flag.get(s, False)) if s is not None else False)

_no_specimen = qpcr_all[qpcr_all["specimen"].isna()]
say("qPCR rows with no crosswalk specimen: " + str(len(_no_specimen)))
for _, _row in _no_specimen.iterrows():
    say("    " + str(_row.get("host_id")) + "  (Barcode "
        + str(_row.get("Barcode")) + ")")

pooled = qpcr_all[qpcr_all["subsamples"] == 1].copy()
# is_single_tick, not num_ticks == 1. The flag is the project's one stated rule
# (build_4marker_crosswalk.py docstring note 4) and additionally excludes the
# Batch-12 fungal-treatment ticks and the specimens whose extraction record
# contradicts the collection record. The sequencing cohort already obeys it.
_by_count = pooled[pooled["num_ticks"] == 1]
single_ticks = pooled[pooled["is_single_tick"]].copy()
say("one-subsample rows from a one-tick tube:          " + str(len(_by_count)))
say("of those, eligible by the is_single_tick rule:    " + str(len(single_ticks)))
if len(_by_count) != len(single_ticks):
    _dropped = sorted(set(_by_count["specimen"].dropna())
                      - set(single_ticks["specimen"].dropna()))
    say("  excluded by the rule: " + ", ".join(_dropped))
if len(single_ticks) == 0:
    sys.exit("ERROR: no rows with one subsample from a one-tick tube. Nothing to plot.")

n_unmatched = int(pooled["num_ticks"].isna().sum())
say(f"qPCR rows: {len(qpcr_all)}")
say(f"  subsample_count == 1 (POOLED, roughly 3 ticks each): {len(pooled)}")
say(f"  of those, tube held exactly 1 tick (TRUE per tick):  {len(single_ticks)}")
if n_unmatched:
    say(f"  {n_unmatched} of the {len(pooled)} could not be matched to a tube tick count")

say("\n=== pool-level positivity, for comparison, NOT plotted ===")
for target in TARGET_LABELS:
    k = int((pd.to_numeric(pooled[target], errors="coerce").fillna(0) > 0).sum())
    say(f"  {target:22} {k:>4} / {len(pooled)} = {100.0*k/len(pooled):5.1f}%")

# For a single tick every target is 0 or 1 (checked 2026-08-11). I assert it
# rather than trust it: a stray Ct value or a count would silently inflate
# every prevalence in this figure.
for target in TARGET_LABELS:
    values = set(single_ticks[target].dropna().unique())
    unexpected = values - {"0", "1", "0.0", "1.0"}
    if unexpected:
        sys.exit(f"ERROR: {target} holds values other than 0/1 in single ticks: "
                 f"{sorted(unexpected)}. I assumed presence/absence; stopping "
                 f"rather than counting these as positives.")

is_positive = pd.DataFrame(
    {t: pd.to_numeric(single_ticks[t], errors="coerce").fillna(0) > 0
     for t in TARGET_LABELS},
    index=single_ticks.index)

n_ticks = len(single_ticks)

say(f"\n=== qPCR prevalence, {n_ticks} individually screened ticks ===")
prevalence_rows = []
for target in TARGET_LABELS:
    n_positive = int(is_positive[target].sum())
    low, high = wilson_interval(n_positive, n_ticks)
    prevalence_rows.append({"target": target, "n_positive": n_positive,
                            "n_screened": n_ticks,
                            "prevalence_pct": round(100 * n_positive / n_ticks, 1),
                            "wilson_low": round(low, 1), "wilson_high": round(high, 1)})
    say(f"  {target:<20} {n_positive:>4}/{n_ticks}  "
        f"{100*n_positive/n_ticks:5.1f}%  [{low:4.1f}, {high:4.1f}]")

prevalence = pd.DataFrame(prevalence_rows)
prevalence.to_csv(f"{ANALYSIS_DIR}/qpcr_prevalence.tsv", sep="\t", index=False)
say(f"\nwrote {ANALYSIS_DIR}/qpcr_prevalence.tsv")

n_targets_per_tick = is_positive.sum(axis=1)
say(f"\n=== co-infection, distinct qPCR targets per tick (n = {n_ticks}) ===")
for k in sorted(n_targets_per_tick.unique()):
    count = int((n_targets_per_tick == k).sum())
    say(f"  {k} targets: {count:>4}  ({100*count/n_ticks:5.1f}%)")
say(f"  >=1 target : {int((n_targets_per_tick>=1).sum()):>4}  "
    f"({100*(n_targets_per_tick>=1).mean():5.1f}%)")
say(f"  >=2 targets: {int((n_targets_per_tick>=2).sum()):>4}  "
    f"({100*(n_targets_per_tick>=2).mean():5.1f}%)")

figure, (axis_prevalence, axis_coinfection) = plt.subplots(1, 2, figsize=(13, 5.5))

plot_order = prevalence.sort_values("prevalence_pct", ascending=True)
y_positions = range(len(plot_order))
axis_prevalence.barh(list(y_positions), plot_order["prevalence_pct"], color=BLUE)
axis_prevalence.errorbar(
    plot_order["prevalence_pct"], list(y_positions),
    xerr=[plot_order["prevalence_pct"] - plot_order["wilson_low"],
          plot_order["wilson_high"] - plot_order["prevalence_pct"]],
    fmt="none", ecolor="#444444", capsize=3, linewidth=1)
axis_prevalence.set_yticks(list(y_positions))
axis_prevalence.set_yticklabels([TARGET_LABELS[t] for t in plot_order["target"]])
for y, (_, row) in zip(y_positions, plot_order.iterrows()):
    axis_prevalence.text(row["wilson_high"] + 1.5, y,
                         f"n={int(row['n_positive'])}", va="center",
                         fontsize=9, color="#444444")
axis_prevalence.set_xlabel(f"Prevalence (%), n = {n_ticks} individually screened ticks")
axis_prevalence.set_title("A. qPCR pathogen prevalence", loc="left", fontweight="bold")
axis_prevalence.set_xlim(0, max(60, plot_order["wilson_high"].max() * 1.2))
for spine in ("top", "right"):
    axis_prevalence.spines[spine].set_visible(False)

target_counts = [int((n_targets_per_tick == k).sum())
                 for k in range(int(n_targets_per_tick.max()) + 1)]
percentages = [100 * c / n_ticks for c in target_counts]
axis_coinfection.bar([str(k) for k in range(len(percentages))], percentages, color=GREEN)
for k, pct in enumerate(percentages):
    axis_coinfection.text(str(k), pct + 0.8, f"{pct:.1f}%", ha="center", fontsize=9)
axis_coinfection.set_xlabel("qPCR targets detected per tick")
axis_coinfection.set_ylabel("Ticks (%)")
axis_coinfection.set_title("B. Co-infection (qPCR panel targets)", loc="left",
                           fontweight="bold")
axis_coinfection.set_ylim(0, max(percentages) * 1.2)
axis_coinfection.text(
    0.97, 0.95,
    f"{100*(n_targets_per_tick>=1).mean():.1f}% carry \u22651\n"
    f"{100*(n_targets_per_tick>=2).mean():.1f}% carry \u22652",
    transform=axis_coinfection.transAxes, ha="right", va="top", fontsize=10,
    bbox=dict(boxstyle="round", fc="#EEEEEE", ec="none"))
for spine in ("top", "right"):
    axis_coinfection.spines[spine].set_visible(False)

figure.tight_layout()
for extension in ("png", "pdf"):
    path = f"{FIGURE_DIR}/fig_qpcr_prevalence_coinfection.{extension}"
    figure.savefig(path, dpi=300, bbox_inches="tight")
    say(f"wrote {path}")
plt.close(figure)

say("\nCAPTION MUST STATE:")
say(f"  - {n_ticks} individually screened ticks; this is a per-tick prevalence")
say("  - these ticks are a SUBSET of the sequenced single ticks, not a")
say("    different population. The cohort selects on is_single_tick, the same")
say("    rule the sequencing cohort uses, so these are exactly the same ticks")
say("    as the qPCR/metabarcoding concordance cohort. The sequenced specimens")
say("    absent here have only pooled qPCR rows. See ground truth 1.15.")
say("  - the two panel-B co-infection figures ARE the same ticks, but they")
say("    count different things; see the next point.")
say("  - this panel counts SPECIES across a 7-target qPCR panel; the")
say("    metabarcoding panel counts GENERA with both Babesia collapsed and")
say("    Bartonella excluded as a ked-associated co-detection (15.1% there).")
say("    Part of the difference between the two is definitional, not biological.")
say("  - all seven panel targets shown, including the zeros, so panel coverage is")
say("    visible")
say(f"  - Wilson 95% intervals; correct for proportions at small n (Powassan")
say(f"    is 1 of {n_ticks})")

with open(f"{ANALYSIS_DIR}/qpcr_report.txt", "w") as handle:
    handle.write("\n".join(report_lines))
say(f"\nreport -> {ANALYSIS_DIR}/qpcr_report.txt")
