#!/usr/bin/env python3
"""recompute_cooccurrence.py

I recompute the single-tick co-infection and pairwise co-occurrence statistics
from the post-anchoring single-tick table, so the manuscript's section 3.3 and
discussion numbers can be replaced with values that match the locked n = 93.

WHY
---
The published pairwise numbers (Borreliella x Anaplasma 11 of 91, OR 2.20,
Fisher p = 0.130) were computed before two things changed: the cohort went
from 91 to 93, and the Anaplasma pattern was anchored on the genus rank,
which removed 3 Wolbachia-driven false positives. Both move these cells, so
the old values cannot be carried into the manuscript.

INPUT
-----
single_tick_host_pathogen.tsv, 93 rows, the *_pos columns, which are the
anchored calls regenerated 2026-09-10.

WHAT I COUNT AS A PATHOGEN GENUS
--------------------------------
Anaplasma, Borreliella, Borrelia and Babesia, with the three Babesia columns
collapsed to one genus as the locked numbers do. Bartonella is EXCLUDED: it
was reclassified as a ked-associated co-detection and not a pathogen finding
(ground truth section 7, restated 2026-09-17). Rickettsia, Francisella and
Wolbachia are endosymbionts and are excluded.

OUTPUT
------
Prints the co-infection distribution, every pairwise 2x2 with odds ratio,
two-sided Fisher exact p and Benjamini-Hochberg q, and writes a TSV.
"""
import itertools
import os
import sys
from math import comb

import pandas as pd

TABLE = "results/single_tick_v3/single_tick_host_pathogen.tsv"
OUT = "results/single_tick_v3/cooccurrence_n93.tsv"

GENERA = {
    "Anaplasma": ["16S_Anaplasma_pos"],
    "Borreliella": ["16S_Borreliella_pos"],
    "Borrelia": ["16S_Borrelia_pos"],
    "Babesia": ["18S_Babesia_odocoilei_pos", "18S_Babesia_microti_pos",
                "18S_Babesia_sp_pos"],
}


def fisher_two_sided(a, b, c, d):
    """I compute the two-sided Fisher exact p for a 2x2 by total enumeration.

    scipy is available here, but I want this script to run on ron inside the
    qiime2 environment without depending on it, and a 2x2 with n = 93 is
    small enough that summing the hypergeometric tail directly is exact and
    instant.
    """
    row1, row2 = a + b, c + d
    col1 = a + c
    total = row1 + row2

    def probability(x):
        return (comb(row1, x) * comb(row2, col1 - x)) / comb(total, col1)

    observed = probability(a)
    low = max(0, col1 - row2)
    high = min(row1, col1)
    # I add a relative tolerance so tables equally extreme as the observed one
    # are not dropped by floating-point noise.
    return sum(probability(x) for x in range(low, high + 1)
               if probability(x) <= observed * (1 + 1e-9))


def odds_ratio(a, b, c, d):
    """I return the Haldane-corrected odds ratio, and None if it is undefined."""
    if 0 in (a, b, c, d):
        if a + b == 0 or c + d == 0 or a + c == 0 or b + d == 0:
            return None
        return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return (a * d) / (b * c)


def benjamini_hochberg(pvalues):
    """I return BH-adjusted q values in the original order."""
    indexed = sorted(enumerate(pvalues), key=lambda pair: pair[1])
    total = len(pvalues)
    adjusted = [0.0] * total
    running = 1.0
    for rank in range(total, 0, -1):
        position, pvalue = indexed[rank - 1]
        running = min(running, pvalue * total / rank)
        adjusted[position] = min(running, 1.0)
    return adjusted


def main():
    if not os.path.isfile(TABLE):
        sys.exit("input not found: " + TABLE)
    frame = pd.read_csv(TABLE, sep="\t")

    calls = pd.DataFrame(index=frame.index)
    for genus, columns in GENERA.items():
        missing = [c for c in columns if c not in frame.columns]
        if missing:
            sys.exit("column(s) missing from " + TABLE + ": "
                     + ", ".join(missing))
        calls[genus] = frame[columns].fillna(0).astype(int).max(axis=1)

    n = len(calls)
    if n != 93:
        sys.exit("cohort is " + str(n) + " rows, not the locked 93. "
                 "TABLE points at " + TABLE + ". The copy under "
                 "results/single_tick_v2/ is the n=103 pre-exclusion "
                 "table; use the v3 one.")
    print("cohort n =", n)
    print()
    print("=== per-genus prevalence, pathogens only ===")
    for genus in GENERA:
        positives = int(calls[genus].sum())
        print("  {:<12} {:>2} / {}  = {:>5.1f}%".format(
            genus, positives, n, 100 * positives / n))
    print()

    burden = calls.sum(axis=1)
    print("=== co-infection, distinct pathogen genera ===")
    for level in range(0, int(burden.max()) + 1):
        count = int((burden == level).sum())
        print("  {} genera: {:>3}  ({:>5.1f}%)".format(
            level, count, 100 * count / n))
    for threshold in (1, 2, 3):
        count = int((burden >= threshold).sum())
        print("  >= {}    : {:>3}  ({:>5.1f}%)".format(
            threshold, count, 100 * count / n))
    infected = int((burden >= 1).sum())
    print("  >= 2 as a share of INFECTED ticks: {} / {} = {:.1f}%".format(
        int((burden >= 2).sum()), infected,
        100 * int((burden >= 2).sum()) / infected))
    print()

    rows = []
    for first, second in itertools.combinations(GENERA, 2):
        x, y = calls[first] == 1, calls[second] == 1
        a = int((x & y).sum())
        b = int((x & ~y).sum())
        c = int((~x & y).sum())
        d = int((~x & ~y).sum())
        rows.append({
            "pair": first + " x " + second,
            "both": a, first + "_only": b, second + "_only": c, "neither": d,
            "both_pct": round(100 * a / n, 1),
            "odds_ratio": odds_ratio(a, b, c, d),
            "fisher_p": fisher_two_sided(a, b, c, d),
        })

    qvalues = benjamini_hochberg([r["fisher_p"] for r in rows])
    for row, q in zip(rows, qvalues):
        row["bh_q"] = q

    print("=== pairwise co-occurrence ===")
    for row in sorted(rows, key=lambda r: -r["both"]):
        ratio = row["odds_ratio"]
        print("  {:<26} both {:>2} of {} ({:>4.1f}%)  OR {}  p = {:.3f}  q = {:.3f}"
              .format(row["pair"], row["both"], n, row["both_pct"],
                      "undef" if ratio is None else "{:.2f}".format(ratio),
                      row["fisher_p"], row["bh_q"]))

    frame_out = pd.DataFrame(rows)
    frame_out.to_csv(OUT, sep="\t", index=False, float_format="%.4f")
    print()
    print("wrote " + OUT)


if __name__ == "__main__":
    main()
