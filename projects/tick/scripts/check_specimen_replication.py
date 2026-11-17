#!/usr/bin/env python3
"""check_specimen_replication.py

I ask how many of the 602 rows in the analysis table are the same animal
counted more than once, and I recompute prevalence and co-occurrence one row
per animal. I write nothing and change nothing.

WHY THIS EXISTS
---------------
check_cooccurrence_batch_heterogeneity.py, 2026-09-18, listed the 43 samples
carrying no run batch. Their names:

    EQ02441291-2.241209_P3-5   ...   EQ02441291-10.241209_P3-5     9 rows
    EQ02442207-2 ... -9                                            8 rows
    EQ02442188-2 ... -7                                            6 rows
    EQ02442189-2 ... -6                                            5 rows
    EQ02442227-2 ... -5                                            4 rows
    EQ02442197, EQ02442213, EQ02442218                             3 rows each
    EQ02441306-2, -3                                               2 rows

Forty-three rows from nine EQ numbers. Two things follow.

First, those 43 are not unbatched. '241209' is 2024-12-09 and the batch string
is AMP-120924-USGS-TICKS-P3-5, the same date the other way round, and both say
P3-5. So the plate that carries the co-occurrence is 96 + 43 = 139 rows, not
96, and my earlier script held 43 of them aside as unknown.

Second, and worse: if several rows are the same animal, they are not
independent. One deer whose ticks carry both Anaplasma and Borreliella
contributes nine co-positive rows, and that manufactures co-occurrence with no
biology in it. It is what an odds ratio of 8.23 in that group looks like.

That reaches past co-occurrence. If the 602 holds repeated animals, every
prevalence on 602 is weighted by how often each animal was subsampled rather
than by how many animals there are.

WHAT I DO
---------
1. Identify each row's specimen THREE ways and check they agree: the Barcode
   column, the library lookup, and the bare name prefix. If they disagree I
   say so and do not pick a winner, because the whole point is not to guess.

2. Count rows per specimen, and show every specimen holding more than one.

3. Report whether the repeated rows are SUBSAMPLES (the animal's ticks split
   between tubes) or REPLICATES (the same material sequenced twice), by
   looking at whether NumTicks differs between the rows. I do not assume. The
   answer changes what collapsing means and I cannot get it from this file
   alone.

4. Recompute, one row per specimen, taking a specimen positive if ANY of its
   rows is positive: per-pathogen prevalence, and the co-occurrence odds
   ratios with Fisher exact tests. Then print row-level and specimen-level
   side by side.

WHAT I DO NOT DO
----------------
I do not write a corrected analysis table, and I do not decide what NumTicks
should be for a collapsed specimen. Summing is right for subsamples and wrong
for replicates.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_specimen_replication.py
"""
import itertools
import os
import re
import sys

import pandas as pd
from scipy import stats

TABLE = "results/analysis/analysis_table_16S.tsv"
PATHOGENS = ["Anaplasma", "Borreliella", "Borrelia", "Ehrlichia"]
CO_DETECTIONS = ["Bartonella"]
BATCH_COLUMN = "RunBatch"

# A sample id like 'EQ02441291-10.241209_P3-5' is specimen EQ02441291,
# subsample 10, run 241209_P3-5. I strip the run first, then a trailing
# subsample index, and I only strip an index of one or two digits so that a
# specimen whose name genuinely ends in a long number is left alone.
RUN_SUFFIX = re.compile(r"\..*$")
SUBSAMPLE_SUFFIX = re.compile(r"-(\d{1,2})$")


def blank(value):
    """I return '' for anything meaning 'no value', and the text otherwise."""
    text = str(value).strip()
    return "" if text.lower() in ("nan", "none", "<na>", "") else text


def prefix_of(sample):
    """I return the specimen prefix of a sample id, by name alone."""
    stem = RUN_SUFFIX.sub("", str(sample).strip())
    return SUBSAMPLE_SUFFIX.sub("", stem)


def subsample_index_of(sample):
    """I return the trailing subsample index of a sample id, or '' if none."""
    stem = RUN_SUFFIX.sub("", str(sample).strip())
    found = SUBSAMPLE_SUFFIX.search(stem)
    return found.group(1) if found else ""


def load_lookup_specimens(frame):
    """I resolve each row to a specimen through the library lookup, or None.

    I import tick_lookup the way the build scripts do. If it is not importable
    I say so and carry on with the other two identifications rather than
    stopping, because two out of three still answers the question.
    """
    try:
        import tick_lookup
    except ImportError as problem:
        print("note: could not import tick_lookup (%s). Continuing without it."
              % problem)
        return None
    loader = getattr(tick_lookup, "load", None)
    if loader is None or not hasattr(tick_lookup, "resolve"):
        print("note: tick_lookup has no load()/resolve(). Continuing without"
              " it.")
        return None
    table = loader()
    resolved = []
    for sample in frame["sample"]:
        stem = RUN_SUFFIX.sub("", str(sample).strip())
        hit = None
        for marker, value in (("16S", sample), ("16S", stem),
                              ("specimen", stem),
                              ("specimen", prefix_of(sample))):
            hit = tick_lookup.resolve(table, marker, value)
            if hit is not None:
                break
        resolved.append(blank(hit))
    return resolved


def agreement_report(frame):
    """I check the three identifications against each other and report."""
    print("=" * 72)
    print("HOW MANY WAYS DO WE GET THE SAME ANSWER?")
    print("=" * 72)
    columns = [name for name in ("_barcode", "_lookup", "_prefix")
               if frame[name].str.len().sum() > 0]
    for name in columns:
        filled = int((frame[name] != "").sum())
        print("  %-10s fills %d of %d rows, %d distinct values"
              % (name.lstrip("_"), filled, len(frame),
                 frame.loc[frame[name] != "", name].nunique()))
    for first, second in itertools.combinations(columns, 2):
        both = frame[(frame[first] != "") & (frame[second] != "")]
        if not len(both):
            continue
        # I compare the GROUPINGS, not the strings: two schemes agree if they
        # cut the rows into the same sets, whatever they call them. So I count
        # groups that one scheme splits and the other does not, in both
        # directions. Comparing group NUMBERS would be nonsense, because the
        # numbering follows each scheme's own sort order.
        one, two = first.lstrip("_"), second.lstrip("_")
        split_by_second = int((both.groupby(first)[second].nunique() > 1).sum())
        split_by_first = int((both.groupby(second)[first].nunique() > 1).sum())
        print("  %s vs %s, on the %d rows carrying both:" % (one, two, len(both)))
        if not split_by_first and not split_by_second:
            print("      identical groupings. %d groups either way."
                  % both[first].nunique())
        else:
            print("      %d %s groups are split by %s; %d %s groups are split"
                  " by %s" % (split_by_second, one, two,
                              split_by_first, two, one))
        difference = both[first].nunique() - both[second].nunique()
        if difference > 0:
            print("      %s makes %d FEWER groups than %s, so it is MERGING"
                  % (two, difference, one))
            print("      rows. If %s is the name prefix, suspect animals whose"
                  % two)
            print("      names end in a number: ALHT-Deer-10 and ALHT-Deer-11")
            print("      collapse into ALHT-Deer.")
        elif difference < 0:
            print("      %s makes %d MORE groups than %s, so it is SPLITTING"
                  % (two, -difference, one))
            print("      rows that %s treats as one animal." % one)
    print()


def replication_report(frame, key):
    """I count rows per specimen and show every specimen with more than one."""
    print("=" * 72)
    print("ROWS PER SPECIMEN, by %s" % key.lstrip("_"))
    print("=" * 72)
    usable = frame[frame[key] != ""]
    missing = len(frame) - len(usable)
    if missing:
        print("  %d rows carry no %s and are counted on their own."
              % (missing, key.lstrip("_")))
    counts = usable.groupby(key).size().sort_values(ascending=False)
    distribution = counts.value_counts().sort_index()
    print("  %d rows resolve to %d specimens." % (len(usable), len(counts)))
    print()
    print("  rows per specimen   specimens   rows")
    for rows_each, specimens in distribution.items():
        print("      %-15d %-11d %d" % (rows_each, specimens,
                                        rows_each * specimens))
    repeated = counts[counts > 1]
    print()
    print("  %d specimens appear more than once, accounting for %d of %d rows"
          " (%.1f%%)."
          % (len(repeated), int(repeated.sum()), len(frame),
             100.0 * repeated.sum() / len(frame)))
    if not len(repeated):
        return counts
    print()
    print("  SUBSAMPLES OR REPLICATES?")
    print("  If every row of a specimen carries the same NumTicks, the tick")
    print("  count is the ANIMAL's, copied onto each row, and these rows are")
    print("  not independent draws of different ticks. If NumTicks differs")
    print("  between rows, the animal's ticks were split between tubes.")
    print()
    header = "  %-22s %5s %8s %-20s %s"
    print(header % ("specimen", "rows", "batches", "NumTicks per row",
                    "subsample ids"))
    identical = differing = 0
    for name in repeated.index[:40]:
        sub = usable[usable[key] == name]
        ticks = sorted(set(pd.to_numeric(sub["NumTicks"],
                                         errors="coerce").dropna()))
        if len(ticks) == 1:
            identical += 1
        else:
            differing += 1
        batches = set(b for b in sub["_batch"] if b)
        batch_text = str(len(batches))
        indices = ",".join(sorted(
            (subsample_index_of(s) or "-") for s in sub["sample"]))
        ticks_text = ",".join("%g" % t for t in ticks)
        print(header % (str(name)[:22], len(sub), batch_text,
                        ticks_text[:20], indices[:40]))
    if len(repeated) > 40:
        print("  ... and %d more" % (len(repeated) - 40))
    print()
    print("  %d repeated specimens carry ONE NumTicks across their rows"
          " (replicate-shaped)." % identical)
    print("  %d carry DIFFERENT NumTicks across their rows (subsample-shaped)."
          % differing)
    print()
    return counts


def collapse(frame, key, columns):
    """I return one row per specimen, positive if ANY of its rows is positive."""
    work = frame.copy()
    work.loc[work[key] == "", key] = ["_row_" + str(i) for i in
                                      work.index[work[key] == ""]]
    aggregated = work.groupby(key)[columns].max()
    return aggregated


def prevalence_report(frame, key, present):
    """I print row-level and specimen-level prevalence side by side."""
    columns = [p + "_pos" for p in present]
    collapsed = collapse(frame, key, columns)
    print("=" * 72)
    print("PREVALENCE: one row per SAMPLE against one row per SPECIMEN")
    print("=" * 72)
    print("  A specimen counts as positive if ANY of its rows is positive.")
    print()
    print("  %-14s %-22s %-22s %s"
          % ("", "per sample (n=%d)" % len(frame),
             "per specimen (n=%d)" % len(collapsed), "shift"))
    for pathogen in present:
        column = pathogen + "_pos"
        rows_k = int((frame[column] == 1).sum())
        spec_k = int((collapsed[column] == 1).sum())
        rows_pct = 100.0 * rows_k / len(frame)
        spec_pct = 100.0 * spec_k / len(collapsed)
        print("  %-14s %-22s %-22s %+.1f pts"
              % (pathogen, "%.1f%% (%d/%d)" % (rows_pct, rows_k, len(frame)),
                 "%.1f%% (%d/%d)" % (spec_pct, spec_k, len(collapsed)),
                 spec_pct - rows_pct))
    print()
    return collapsed


def cooccurrence_report(frame, collapsed, present):
    """I print the co-occurrence at sample level and at specimen level."""
    print("=" * 72)
    print("CO-OCCURRENCE: one row per SAMPLE against one row per SPECIMEN")
    print("=" * 72)
    print("  %-28s %-30s %s" % ("pair", "per sample", "per specimen"))
    for first, second in itertools.combinations(present, 2):
        line = []
        for source, label in ((frame, "sample"), (collapsed, "specimen")):
            x = source[first + "_pos"] == 1
            y = source[second + "_pos"] == 1
            a = int((x & y).sum())
            b = int((x & ~y).sum())
            c = int((~x & y).sum())
            d = int((~x & ~y).sum())
            ratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
            line.append("OR %.2f p=%.3g (both %d)" % (ratio, pvalue, a))
        print("  %-28s %-30s %s"
              % ("%s x %s" % (first, second), line[0], line[1]))
    print()
    print("  Fisher exact, two-sided, uncorrected. These are the same")
    print("  comparisons pathogen_cooccurrence.py runs; I print them here only")
    print("  to show what collapsing does to them, not to replace that script.")
    print()


def plate_report(frame, key):
    """I report the P3-5 plate, which is where the association was living."""
    print("=" * 72)
    print("THE P3-5 PLATE")
    print("=" * 72)
    is_p35 = (frame["_batch"].str.contains("P3-5", case=False, regex=False)
              | frame["sample"].astype(str).str.contains("P3-5", case=False,
                                                         regex=False))
    plate = frame[is_p35]
    if not len(plate):
        print("  no rows match P3-5. Nothing to report.")
        print()
        return
    named = plate[plate["_batch"] != ""]
    unnamed = plate[plate["_batch"] == ""]
    specimens = plate.loc[plate[key] != "", key].nunique()
    print("  %d of %d rows are P3-5 by batch string or by sample name."
          % (len(plate), len(frame)))
    print("    %d carry the batch string, %d carry it only in the sample name"
          % (len(named), len(unnamed)))
    print("  those rows resolve to %d specimens" % specimens)
    if specimens:
        print("  which is %.1f rows per animal, against %.1f for the rest"
              % (len(plate) / float(specimens),
                 _rows_per_animal(frame[~is_p35], key)))
    print()


def _rows_per_animal(frame, key):
    """I return rows per specimen for a subset, or 0 if it cannot be measured."""
    usable = frame[frame[key] != ""]
    if not len(usable):
        return 0.0
    return len(usable) / float(usable[key].nunique())


def main():
    if not os.path.isfile(TABLE):
        sys.exit("not found: " + TABLE
                 + "\nRun this from the project root, after "
                   "python3 scripts/build_analysis_table.py")
    frame = pd.read_csv(TABLE, sep="\t")
    for column in ("sample",):
        if column not in frame.columns:
            sys.exit(TABLE + " has no '" + column + "' column.")
    frame = frame.copy()
    frame["_batch"] = ([blank(v) for v in frame[BATCH_COLUMN]]
                       if BATCH_COLUMN in frame.columns else "")
    frame["_barcode"] = ([blank(v) for v in frame["Barcode"]]
                         if "Barcode" in frame.columns else "")
    frame["_prefix"] = [prefix_of(s) for s in frame["sample"]]
    resolved = load_lookup_specimens(frame)
    frame["_lookup"] = resolved if resolved is not None else ""

    present = [p for p in PATHOGENS + CO_DETECTIONS
               if p + "_pos" in frame.columns
               and frame[p + "_pos"].sum() > 0]
    if len(present) < 2:
        sys.exit("fewer than two organisms with detections in " + TABLE)

    print("samples: %d" % len(frame))
    print("organisms with detections: %s" % ", ".join(present))
    print()

    agreement_report(frame)

    # Order matters. The name prefix is a CROSS-CHECK, not an identity: it
    # strips a trailing one- or two-digit index, which is right for
    # 'EQ02441291-10' and catastrophic for 'ALHT-Deer-10', where it would merge
    # two unrelated animals. So it is used only if nothing else is available,
    # and it says so when it is.
    key = ""
    for candidate in ("_lookup", "_barcode", "_prefix"):
        if frame[candidate].str.len().sum() > 0:
            key = candidate
            break
    if not key:
        sys.exit("no way to identify a specimen in " + TABLE
                 + ": no Barcode column, no library lookup, no usable names.")
    if key == "_prefix":
        print("CAUTION: falling back to the NAME PREFIX for specimen identity.")
        print("It strips a trailing 1-2 digit index, which merges names like")
        print("ALHT-Deer-10 and ALHT-Deer-11 into one animal. Read the counts")
        print("below as an upper bound on replication, not as the answer.")
        print()
    print("Using '%s' as the specimen identity below. Check it against the"
          % key.lstrip("_"))
    print("agreement report above before trusting any number that follows.")
    print()

    replication_report(frame, key)
    plate_report(frame, key)
    collapsed = prevalence_report(frame, key, present)
    cooccurrence_report(frame, collapsed, present)

    print("=" * 72)
    print("WHAT THIS DOES AND DOES NOT SETTLE")
    print("=" * 72)
    print("  If the specimen-level co-occurrence falls to about 1, the pooled")
    print("  association was repeated animals, not pathogens, and it does not")
    print("  go in the manuscript at all.")
    print()
    print("  If prevalence moves between the two columns, then every")
    print("  prevalence reported on 602 samples is weighted by how often each")
    print("  animal was subsampled. That reaches Table 1, not just the")
    print("  co-occurrence, and the analysis table needs a decision about its")
    print("  unit before anything else is written on it.")
    print()
    print("  Collapsing with ANY-positive is not obviously right either. It")
    print("  raises prevalence whenever a specimen has several rows, because")
    print("  more tubes means more chances to detect. The honest comparison")
    print("  needs the unit fixed first, which is a question for the run sheet")
    print("  and for whoever set the subsampling up, not for this script.")


if __name__ == "__main__":
    main()
