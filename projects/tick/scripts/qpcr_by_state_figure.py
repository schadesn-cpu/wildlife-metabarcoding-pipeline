#!/usr/bin/env python3
"""qPCR pathogen detection by state, single ticks.

Built from crosswalk_4marker.tsv, which is the corrected state source. An earlier version of
this figure showed NH 120 / VT 193 / ME 38; those counts are not reproducible from any file
in the project and do not match the crosswalk (NH 125 / VT 212 / ME 37). This script exists
so the numbers have a source someone else can re-derive.

Two things the earlier figure did not show and this one does:
  - 47 of the 423 single ticks have NO recorded state, and 2 are recorded 'NONE'. That is
    11% of the set. Leaving them out silently makes the bars look like full coverage.
  - Maine is 37 ticks against Vermont's 212. A percentage on 37 moves 2.7 points per tick,
    so the ME bars are marked and the caption says not to read them against the others.

Percentages are within-state: the denominator is that state's single ticks, not the whole
set. That is the only way three unequal groups can be compared at all, and it is still a
comparison of three very different sample sizes.

Run:
    cd ~/meedlab/tick_project
    python3 scripts/qpcr_by_state_figure.py \
        --out results/figures/fig_qpcr_by_state

Writes:
    results/figures/fig_qpcr_by_state.{png,pdf}
    results/analysis/qpcr_by_state.tsv
"""

import argparse
import os
import re
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORKBOOK = "metadata/tick_full_metadata.xlsx"
QPCR_SHEET = "qPCR Results"
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"

TARGETS = [
    ("a-phagocytophilum", r"$\mathit{A.\ phagocytophilum}$"),
    ("b-burgdorferi",     r"$\mathit{B.\ burgdorferi}$"),
    ("b-microti",         r"$\mathit{B.\ microti}$"),
    ("b-miyamotoi",       r"$\mathit{B.\ miyamotoi}$"),
    ("powassan-virus",    "Powassan virus"),
    ("b-mayonii",         r"$\mathit{B.\ mayonii}$"),
    ("e-muris",           r"$\mathit{E.\ muris}$"),
]

STATE_ORDER = ["NH", "VT", "ME"]
# Wong colourblind-safe
STATE_COLOUR = {"NH": "#0072B2", "VT": "#009E73", "ME": "#E69F00"}

# Below this a within-state percentage moves more than 2 points per tick, so the bar is
# marked in the legend and the caption says not to compare it against the others.
SMALL_GROUP = 50

report_lines = []


def say(message=""):
    print(message)
    report_lines.append(str(message))


def normalise(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def back_up_existing(out_path):
    if not os.path.exists(out_path):
        return
    backup_path = out_path + ".bak"
    if os.path.exists(backup_path):
        sys.exit(f"{out_path} and {backup_path} both exist. "
                 f"Move or delete {backup_path} before rerunning.")
    os.replace(out_path, backup_path)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="results/figures/fig_qpcr_by_state")
    args = parser.parse_args()

    for path in (WORKBOOK, CROSSWALK):
        if not os.path.isfile(path):
            sys.exit(f"no such file: {path}")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    os.makedirs("results/analysis", exist_ok=True)

    qpcr = pd.read_excel(WORKBOOK, sheet_name=QPCR_SHEET, dtype=str)
    qpcr["subsamples"] = pd.to_numeric(qpcr["subsample_count"], errors="coerce")
    single = qpcr[qpcr["subsamples"] == 1].copy()
    say(f"qPCR rows: {len(qpcr)}   single ticks: {len(single)}")

    missing = [t for t, _ in TARGETS if t not in single.columns]
    if missing:
        sys.exit(f"{QPCR_SHEET} is missing target columns: {missing}")

    crosswalk = pd.read_csv(CROSSWALK, sep="\t")
    state_by_barcode = {}
    for _, row in crosswalk.iterrows():
        for key in (row.get("barcode"), row.get("sample_id")):
            if isinstance(key, str) and key.strip():
                state_by_barcode[normalise(key)] = str(row.get("state", "")).strip()

    single["state"] = single["host_id"].map(lambda v: state_by_barcode.get(normalise(v), ""))
    single["state"] = single["state"].replace({"nan": "", "NONE": ""})

    say("\nsingle ticks by state (from the crosswalk):")
    counts = single["state"].replace("", "(no state recorded)").value_counts()
    say(counts.to_string())

    n_no_state = int((single["state"] == "").sum())
    say(f"\n{n_no_state} of {len(single)} single ticks have no recorded state "
        f"({100*n_no_state/len(single):.0f}%). They are excluded from the figure and the "
        f"caption says so.")

    rows = []
    for state in STATE_ORDER:
        subset = single[single["state"] == state]
        if subset.empty:
            say(f"warning: no single ticks for {state}")
            continue
        for target, label in TARGETS:
            positive = (subset[target].astype(str).str.strip().isin(["1", "1.0"])).sum()
            rows.append({"state": state, "target": target, "label": label,
                         "n_positive": int(positive), "n_screened": len(subset),
                         "pct": round(100 * positive / len(subset), 1)})
    if not rows:
        sys.exit("no state had any single ticks")

    table = pd.DataFrame(rows)
    out_tsv = "results/analysis/qpcr_by_state.tsv"
    back_up_existing(out_tsv)
    table.to_csv(out_tsv, sep="\t", index=False)
    say(f"\nwrote {out_tsv}")

    say("\n=== positives per state ===")
    say(table.pivot(index="target", columns="state", values="n_positive")
             .reindex([t for t, _ in TARGETS]).to_string())
    say("\n=== percent within state ===")
    say(table.pivot(index="target", columns="state", values="pct")
             .reindex([t for t, _ in TARGETS]).to_string())

    # ---- figure -------------------------------------------------------------
    present = [s for s in STATE_ORDER if (table["state"] == s).any()]
    n_by_state = {s: int(table[table["state"] == s]["n_screened"].iloc[0]) for s in present}

    figure, axis = plt.subplots(figsize=(11, 6))
    height = 0.8 / len(present)
    positions = range(len(TARGETS))

    for offset, state in enumerate(present):
        subset = table[table["state"] == state].set_index("target")
        values = [subset.loc[t, "pct"] for t, _ in TARGETS]
        counts_ = [subset.loc[t, "n_positive"] for t, _ in TARGETS]
        y = [p + (len(present) / 2 - offset - 0.5) * height for p in positions]
        small = n_by_state[state] < SMALL_GROUP
        axis.barh(y, values, height=height * 0.9, color=STATE_COLOUR[state],
                  label=f"{state} (n = {n_by_state[state]})" + ("*" if small else ""))
        for yi, value, count in zip(y, values, counts_):
            axis.text(value + 0.6, yi, str(int(count)), va="center", fontsize=8,
                      color="#333333")

    axis.set_yticks(list(positions))
    axis.set_yticklabels([label for _, label in TARGETS], fontsize=10)
    axis.invert_yaxis()
    axis.set_xlabel("Ticks positive within that state (%)", fontsize=11)
    axis.set_title("qPCR pathogen detection by state, individually screened ticks",
                   fontsize=12, fontweight="bold", loc="left")
    axis.legend(loc="lower right", frameon=True, fontsize=9)
    axis.set_xlim(0, max(60, table["pct"].max() * 1.18))
    axis.grid(axis="x", color="#DDDDDD", lw=0.6)
    axis.set_axisbelow(True)
    for spine in ("top", "right"):
        axis.spines[spine].set_visible(False)

    small_states = [s for s in present if n_by_state[s] < SMALL_GROUP]
    note = ("Bars are the percentage of that state's individually screened ticks; numbers "
            "beside each bar are tick counts. Denominators differ, so the bars compare "
            "rates within a state and not sample sizes.\n"
            f"{n_no_state} of {len(single)} single ticks have no recorded state and are not "
            "shown. State assignments come from crosswalk_4marker.tsv.")
    if small_states:
        note += ("\n* " + ", ".join(small_states) +
                 f" has fewer than {SMALL_GROUP} ticks: one tick moves the percentage by "
                 "more than two points, so read those bars as counts rather than rates.")
    figure.text(0.01, -0.02, note, fontsize=7.5, style="italic", color="#666666", va="top")

    figure.tight_layout()
    for extension in ("png", "pdf"):
        path = f"{args.out}.{extension}"
        figure.savefig(path, dpi=300, bbox_inches="tight")
        say(f"wrote {path}")
    plt.close(figure)

    report_path = "results/analysis/qpcr_by_state_report.txt"
    with open(report_path, "w") as handle:
        handle.write("\n".join(report_lines))
    print(f"report -> {report_path}")


if __name__ == "__main__":
    main()
