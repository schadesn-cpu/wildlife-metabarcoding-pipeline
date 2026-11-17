#!/usr/bin/env python3
"""check_control_enrichment.py

I compare each genus's share of the negative controls against its share of the
whole sequencing pool, which separates reagent contamination from index
hopping. I write only with --write.

WHY THIS AND NOT THE LAST ONE
-----------------------------
check_control_power.py asked how many reads a blank should hold if it carried
the tick community of one panel. That was the wrong model. A blank holds no
template, so nearly everything in it arrived from other wells, and its
composition follows the pool rather than any one panel. With 604 deer
libraries against 13 moose, expecting a blank to look like a moose is why that
script predicted 14,808 Enhydrobacter reads. The real expectation is smaller by
about two orders of magnitude, and every verdict it printed for the moose panel
was scored against an inflated number.

The ratio below needs no such assumption:

    blank share / pool share    about 1   arrived by hopping
                                well above 1   enriched in blanks: reagents
                                well below 1   present in ticks, not in blanks

Rickettsia is 72% of the deer reads and should land near 1. If it does, the
method is calibrated on this dataset and the other rows can be read. If it
does not, nothing here is trustworthy and I say so rather than printing
verdicts.

POWER
-----
A ratio is meaningless when the expected count is a read or two. I print the
expected count beside every ratio and decline to call anything below five.

WHAT I READ
-----------
    qiime2/16S/all/dada2/table.qza, through check_reagent_background

WHAT I WRITE, WITH --write
--------------------------
    results/analysis/control_enrichment.tsv

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_control_enrichment.py
    python3 scripts/check_control_enrichment.py --write
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = "results/analysis/control_enrichment.tsv"
CALIBRATOR = "Rickettsia"
FLOOR = 5.0
ENRICHED = 3.0
DEPLETED = 0.2
SHOW = 20


def background():
    """I return check_reagent_background, whose loading I reuse whole."""
    try:
        import check_reagent_background as module
    except ImportError as problem:
        sys.exit("cannot import check_reagent_background (%s). It must sit "
                 "in scripts/ beside me." % problem)
    return module


def verdict(expected, ratio):
    """I name what a genus's blank-to-pool ratio means, or decline to."""
    if expected < FLOOR:
        return "too few expected to call"
    if ratio >= ENRICHED:
        return "enriched in blanks: reagent background"
    if ratio <= DEPLETED:
        return "in ticks, not in blanks"
    return "at pool rate: index hopping"


def main():
    write = "--write" in sys.argv[1:]
    module = background()
    shared = module.suppression()
    arms = shared.arms()

    raw = module.open_table()
    features = [str(f) for f in raw.ids(axis="observation")]
    _, _, mapping = arms.read_taxonomy(features)
    keep = [f for f in features
            if not any(word in str(mapping.get(f, "")).lower()
                       for word in shared.HOST_WORDS)]
    raw = raw.filter(set(keep), axis="observation", inplace=False)
    genus_of = {}
    for feature in keep:
        genus_of[feature] = arms.genus_label(mapping.get(feature, ""))
    counts = module.collapse_to_genus(raw, genus_of)

    negatives = [i for i in counts.index
                 if module.NEGATIVE.search(i) and not module.POSITIVE.search(i)]
    positives = [i for i in counts.index if module.POSITIVE.search(i)]
    held = set(negatives) | set(positives)
    samples = [i for i in counts.index if i not in held]
    if not negatives or not samples:
        sys.exit("need both blanks and samples; found %d and %d."
                 % (len(negatives), len(samples)))

    blank = counts.loc[negatives].sum(axis=0)
    pool = counts.loc[samples].sum(axis=0)
    blank_total, pool_total = float(blank.sum()), float(pool.sum())
    pool_share = pool / pool_total
    expected = pool_share * blank_total
    ratio = (blank / expected.replace(0, float("nan")))
    prevalence = (counts.loc[negatives] > 0).sum(axis=0)

    print("%d blanks holding %.0f reads, against %d sample libraries holding "
          "%.0f reads" % (len(negatives), blank_total, len(samples), pool_total))
    print()

    if CALIBRATOR not in counts.columns or expected[CALIBRATOR] < FLOOR:
        sys.exit("%s is not in the table often enough to calibrate on. Stop "
                 "here." % CALIBRATOR)
    calibration = float(ratio[CALIBRATOR])
    print("CALIBRATION. %s is %.1f%% of the pool and %.1f%% of the blanks, a "
          "ratio of %.2f."
          % (CALIBRATOR, 100 * pool_share[CALIBRATOR],
             100 * blank[CALIBRATOR] / blank_total, calibration))
    if 0.5 <= calibration <= 2.0:
        print("That is near one, so the blanks track the pool and hopping is "
              "the dominant source in them. The rows below can be read.")
    else:
        print("That is not near one. The blanks do not track the pool, so the "
              "ratios below mean something other than what I claim, and I "
              "would not build an argument on them.")
    print()

    table = pd.DataFrame({"blank_reads": blank, "pool_reads": pool,
                          "pool_pct": 100 * pool_share, "expected": expected,
                          "ratio": ratio, "blanks_with": prevalence})
    table = table[table["pool_reads"] > 0]
    table["verdict"] = [verdict(e, r) for e, r
                        in zip(table["expected"], table["ratio"])]

    figure = module.wanted_genera()
    print("=" * 92)
    print("THE GENERA IN THE FIGURE")
    print("=" * 92)
    print("  %-24s %-6s %9s %9s %7s %6s  %s"
          % ("genus", "panel", "expected", "in blanks", "ratio", "of 29",
             "verdict"))
    for name, panels in sorted(figure.items()):
        if name not in table.index:
            print("  %-24s %-6s %s" % (name[:24], panels[-1], "not in the pool"))
            continue
        row = table.loc[name]
        print("  %-24s %-6s %9.0f %9.0f %7.2f %6d  %s"
              % (name[:24], panels[-1], row["expected"], row["blank_reads"],
                 row["ratio"], row["blanks_with"], row["verdict"]))
    print()

    callable_rows = table[table["expected"] >= FLOOR]
    worst = callable_rows.nlargest(SHOW, "ratio")
    print("=" * 92)
    print("THE %d GENERA MOST ENRICHED IN THE BLANKS, WHETHER OR NOT WE ASKED"
          % SHOW)
    print("=" * 92)
    print("  %-30s %9s %9s %7s %6s" % ("genus", "expected", "in blanks",
                                       "ratio", "of 29"))
    for name, row in worst.iterrows():
        print("  %-30s %9.0f %9.0f %7.2f %6d"
              % (str(name)[:30], row["expected"], row["blank_reads"],
                 row["ratio"], row["blanks_with"]))
    print()
    print("  Anything here at a ratio above %.0f is a reagent organism, and "
          "any of them" % ENRICHED)
    print("  that also appears in the manuscript needs saying so.")
    print()

    if not write:
        print("Nothing written. Rerun with --write to save %s" % OUT)
        return
    if os.path.exists(OUT):
        sys.exit("%s already exists. Move it aside first." % OUT)
    if not os.path.isdir(os.path.dirname(OUT)):
        sys.exit("%s does not exist and I will not create it."
                 % os.path.dirname(OUT))
    table.sort_values("ratio", ascending=False).to_csv(OUT, sep="\t")
    print("wrote %s (%d genera)" % (OUT, len(table)))


if __name__ == "__main__":
    main()
