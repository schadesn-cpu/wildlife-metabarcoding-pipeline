#!/usr/bin/env python3
"""patch_confounding_add_batch.py

I teach check_pooled_cooccurrence_confounding.py to stratify on sequencing run
batch, which is the last confounder standing between the pooled
Anaplasma x Borreliella association and a reportable result.

WHY BATCH
---------
The association is OR 2.03 crude, q = 0.0002, and it survives:

    pool size    Mantel-Haenszel 2.05   (flat across all four quartiles)
    read depth   Mantel-Haenszel 1.89   (and the gradient runs the wrong way
                                         for a confounding story)

and it has been stable across the genus anchoring fix, the library-join fix
and two crosswalk rebuilds (2.48 / 2.27 in August on 558 samples against
2.03 / 1.89 today on 602).

Batch has not been tested and is the largest effect in this dataset: R2 10.9%,
p = 0.001, with 530 of 551 animals sitting entirely within one batch (ground
truth 2b, lens 2). If the association is a batch effect, it will collapse here.

WHY A SEPARATE CODE PATH
------------------------
The existing report() cuts a numeric column into quartiles with pd.qcut. Batch
is categorical, so it needs grouping by value rather than binning. Forcing a
categorical through qcut would either fail or silently order strata by string,
which is the kind of quiet wrongness this whole exercise has been about.

I skip strata below 10 samples and say which ones I skipped, and I report how
many strata are INFORMATIVE, meaning they contain at least one discordant pair
so they can contribute to the Mantel-Haenszel ratio at all. A summary built on
two informative strata is not the same claim as one built on nine.

USAGE
-----
    cd ~/meedlab/tick_project
    python3 patch_confounding_add_batch.py
    python3 check_pooled_cooccurrence_confounding.py
"""
import os
import sys

TARGET = "scripts/check_pooled_cooccurrence_confounding.py"
FALLBACK = "check_pooled_cooccurrence_confounding.py"

OLD_REPORT_END = '''    summary = mantel_haenszel(tables)
    if summary is not None:
        print("      Mantel-Haenszel summary OR, holding %s constant: %.2f"
              % (label, summary))
'''

NEW_REPORT_END = '''    summary = mantel_haenszel(tables)
    if summary is not None:
        print("      Mantel-Haenszel summary OR, holding %s constant: %.2f"
              % (label, summary))


def report_categorical(frame, pair, column, label, min_n=10):
    """I stratify on a categorical column, grouping by value rather than binning.

    pd.qcut is for numbers. Pushing a categorical through it would either fail
    or order the strata by string, which is exactly the sort of quiet wrongness
    this script exists to catch.

    Samples with no value get their own line but do NOT enter the summary. They
    are not one stratum, they are an unknown number of strata pooled together,
    so including them would put the very confounding I am trying to remove back
    into the answer.
    """
    first, second = pair
    if column not in frame.columns:
        print("    %-14s column not present; rebuild the analysis table"
              % label)
        return
    work = frame.copy()
    # Not .astype(str).replace(): pandas turns the replacement of the string
    # "nan" back into a real NaN, so the blanks vanish and the missing group
    # silently rejoins the strata. This project has been bitten by the
    # string-versus-NA distinction twice already. A comprehension has no
    # opinions.
    work["_stratum"] = ["" if str(v).strip().lower() in ("nan", "none",
                                                         "<na>", "")
                        else str(v).strip() for v in work[column]]
    missing = work[work["_stratum"] == ""]
    work = work[work["_stratum"] != ""]

    tables, skipped, informative = [], [], 0
    print("    by %s:" % label)
    for level, sub in work.groupby("_stratum", observed=True):
        if len(sub) < min_n:
            skipped.append((str(level), len(sub)))
            continue
        a, b, c, d = two_by_two(sub, first, second)
        tables.append((a, b, c, d))
        if b > 0 and c > 0:
            informative += 1
        print("      %-46s n=%-4d both=%-4d OR=%.2f"
              % (str(level)[:46], len(sub), a, odds_ratio(a, b, c, d)))
    for level, n in skipped:
        print("      %-46s n=%-4d SKIPPED, below %d"
              % (str(level)[:46], n, min_n))
    if len(missing):
        a, b, c, d = two_by_two(missing, first, second)
        print("      %-46s n=%-4d both=%-4d OR=%.2f  NOT in the summary"
              % ("(no " + label + " recorded)", len(missing), a,
                 odds_ratio(a, b, c, d)))
    summary = mantel_haenszel(tables)
    if summary is None:
        print("      no stratum large enough to summarise")
        return
    print("      Mantel-Haenszel summary OR, holding %s constant: %.2f"
          % (label, summary))
    print("      strata in the summary: %d, covering %d of %d samples"
          % (len(tables), sum(sum(t) for t in tables), len(frame)))
    print("      informative strata (both discordant cells non-zero): "
          "%d of %d" % (informative, len(tables)))
    if informative < 3:
        print("      CAUTION: fewer than three informative strata. This summary"
              " rests on very little.")
'''

OLD_CALLS = '''        report(frame, (first, second), "NumTicks", "pool size")
        report(frame, (first, second), "total_reads", "read depth")
        print()'''

NEW_CALLS = '''        report(frame, (first, second), "NumTicks", "pool size")
        report(frame, (first, second), "total_reads", "read depth")
        report_categorical(frame, (first, second), "RunBatch", "run batch")
        print()'''

OLD_TAIL = '''    print("  Either way, animal and batch are almost perfectly confounded in this")
    print("  design (ground truth 2b, lens 2 constraints: 530 of 551 animals sit")
    print("  entirely within one batch), so a surviving association still cannot")
    print("  be attributed to biology without addressing that.")'''

NEW_TAIL = '''    print("  BATCH is now tested too. It is the largest effect in this dataset")
    print("  (R2 10.9%, p = 0.001) and 530 of 551 animals sit entirely within one")
    print("  batch (ground truth 2b, lens 2), so animal and batch cannot be")
    print("  separated here. If the association survives batch as well, the three")
    print("  obvious artifacts are ruled out and it is a real pool-level")
    print("  co-detection. It is STILL not co-infection of one tick, because a")
    print("  pool is a tube, and the animal-batch confound means it still cannot")
    print("  be attributed to tick biology rather than to something shared by")
    print("  ticks off the same animal.")'''

REPLACEMENTS = [
    (OLD_REPORT_END, NEW_REPORT_END, "the end of report()"),
    (OLD_CALLS, NEW_CALLS, "the per-pair report calls"),
    (OLD_TAIL, NEW_TAIL, "the closing note"),
]


def main():
    target = TARGET if os.path.isfile(TARGET) else FALLBACK
    if not os.path.isfile(target):
        sys.exit("not found: " + TARGET + " or " + FALLBACK
                 + "\nRun this from ~/meedlab/tick_project.")
    with open(target, encoding="utf-8") as handle:
        source = handle.read()

    if "report_categorical" in source:
        sys.exit("already patched: " + target + " defines report_categorical. "
                 "Nothing changed.")

    for old, new, label in REPLACEMENTS:
        count = source.count(old)
        if count != 1:
            sys.exit("in " + target + ": expected exactly one occurrence of "
                     + label + ", found " + str(count)
                     + ". Patch by hand rather than trusting this.")
        source = source.replace(old, new, 1)

    compile(source, target, "exec")
    with open(target, "w", encoding="utf-8") as handle:
        handle.write(source)
    print("patched " + target)
    print()
    print("Needs analysis_table_16S.tsv to carry RunBatch, which it does as of")
    print("the 2026-09-18 rebuild. Ten values, n 14 to 96, plus 43 with none.")
    print()
    print("run: python3 " + target)


if __name__ == "__main__":
    main()
