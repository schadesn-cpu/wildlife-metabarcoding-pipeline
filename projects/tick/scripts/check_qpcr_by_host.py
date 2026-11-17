#!/usr/bin/env python3
"""check_qpcr_by_host.py

I tabulate qPCR detection by collection host, separately for true single ticks
and for pools, so the paper can say which pathogen was found on which animal.
I print a summary and write nothing unless asked.

WHERE EVERYTHING COMES FROM
---------------------------
The crosswalk, and only the crosswalk. build_qpcr_figure.py records why:

    the run sheet's 'Num Ticks' is the column the ALHT/LHT barcode clobber
    corrupted, and normalising names by stripping punctuation turns
    'ALHT WTD_09' and 'alht-deer-09' into two different keys.

So I resolve each qPCR row to a crosswalk specimen with tick_lookup, then read
num_ticks, is_single_tick, host_spp and host_bin from the crosswalk. I do not
read the run sheet at all.

COHORTS
-------
    single tick   subsample_count == 1 AND is_single_tick
    pool          subsample_count == 1 AND num_ticks > 1

is_single_tick, not num_ticks == 1: the flag is the project's stated rule and
additionally drops the Batch-12 fungal-treatment ticks and the specimens whose
extraction record contradicts the collection record. Using num_ticks == 1 here
gave 100 ticks against the 89 in results/analysis/qpcr_prevalence.tsv, and the
sequencing cohort obeys the flag, so the flag it is. A one-tick tube that fails
the flag belongs to neither cohort and I say how many there are.

WHAT I READ
-----------
    metadata/tick_full_metadata.xlsx, sheet 'qPCR Results'
    results/metadata/crosswalk_4marker.tsv
    scripts/tick_lookup.py

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_qpcr_by_host.py
    python3 scripts/check_qpcr_by_host.py --write
"""
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

WORKBOOK = "metadata/tick_full_metadata.xlsx"
QPCR_SHEET = "qPCR Results"
CORRECTED = "results/analysis/qpcr_screening_record_corrected.tsv"
CROSSWALK = "results/metadata/crosswalk_4marker.tsv"
OUT = "results/analysis/qpcr_by_host.tsv"

TARGETS = ["a-phagocytophilum", "b-burgdorferi", "b-microti", "b-miyamotoi",
           "powassan-virus", "b-mayonii", "e-muris"]

# I keep the host groups coarse because that is what the sample sizes support.
HOST_GROUPS = [
    ("White-tailed deer", ("white-tailed deer", "whitetail", "white tail",
                           "wtd", "odocoileus", "deer")),
    ("Moose", ("moose", "alces")),
    ("Black bear", ("bear", "ursus")),
    ("Virginia opossum", ("opossum", "didelph")),
    ("Porcupine", ("porcupine", "erethizon")),
    ("Raccoon", ("raccoon", "procyon")),
    ("Domestic animal", ("dog", "canis", "cat", "felis", "horse", "equus",
                         "cow", "bos", "sheep", "goat")),
    ("Human", ("human", "homo")),
    ("Dragged or flagged", ("drag", "flag", "vegetation", "environment")),
    # Everything else that turned up on a check station sheet. These sit last
    # so the named hosts above win when a label carries both.
    ("Other wildlife", ("bobcat", "lynx", "fisher", "coyote", "fox", "vulpes",
                        "otter", "mink", "weasel", "marten", "beaver", "skunk",
                        "woodchuck", "marmot", "squirrel", "chipmunk", "hare",
                        "rabbit", "cottontail", "sylvilagus", "lepus", "turkey",
                        "grouse", "bird", "robin", "goldfinch", "yellowthroat",
                        "warbler", "sparrow", "thrush", "finch", "junco",
                        "chickadee", "mouse", "mice", "vole", "shrew",
                        "muskrat", "coywolf", "wolf", "bat")),
]
UNMAPPED = "Other or unrecorded"
# A row positive for this many targets is a positive control, not a tick.
CONTROL_POSITIVES = 5
SHOW = 12


def host_group(text):
    """I map free host text to one of the groups above, or ''.

    I require a word start before each keyword, so 'cat' does not match
    'Bobcat', but I allow anything after it, so 'didelph' still matches
    'Didelphis virginiana'.
    """
    lowered = str(text).strip().lower()
    if lowered in ("", "nan", "none", "na", "unknown"):
        return ""
    for name, words in HOST_GROUPS:
        for word in words:
            if re.search(r"(?<![a-z])" + re.escape(word), lowered):
                return name
    return UNMAPPED


def named(group):
    """I say whether a group is a real host rather than a fallback."""
    return group not in ("", UNMAPPED)


def lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s). Run me from the project root."
                 % problem)
    return tick_lookup.load(), tick_lookup


def read_qpcr():
    """I return the screening record, counted under the Section 2.5 criteria.

    This read the qPCR Results sheet directly until 2026-10-08. It now reads
    the corrected record written by build_corrected_screening_record.py,
    which is that sheet with the per-target counts recomputed from the
    well-level run files under the testing laboratory's own quality
    criteria, an amplification score above 1.0 and a cycle-threshold
    confidence above 0.6. Twenty-one wells of 6,664 move. The sheet itself
    is untouched and remains the record of what the laboratory supplied.
    """
    if not os.path.exists(CORRECTED):
        sys.exit("missing %s.\nBuild it first, from the project root:"
                 "\n    python3 scripts/build_corrected_screening_record.py"
                 % CORRECTED)
    try:
        frame = pd.read_csv(CORRECTED, sep="\t", dtype=str)
    except (ValueError, OSError) as problem:
        sys.exit("cannot read %s (%s)" % (CORRECTED, problem))
    missing = [c for c in TARGETS if c not in frame.columns]
    if missing:
        sys.exit("%s is missing target column(s): %s" % (CORRECTED, missing))
    if "subsample_count" not in frame.columns:
        sys.exit("%s has no subsample_count column" % QPCR_SHEET)
    return frame


def read_crosswalk():
    """I return per-specimen tick count, single-tick flag and host text."""
    if not os.path.exists(CROSSWALK):
        sys.exit("missing %s. Build it first:\n"
                 "    python3 scripts/build_4marker_crosswalk.py" % CROSSWALK)
    frame = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    for column in ("sample_id", "num_ticks", "is_single_tick", "host_spp",
                   "host_bin"):
        if column not in frame.columns:
            sys.exit("column '%s' missing from %s" % (column, CROSSWALK))
    counts, flags, hosts = {}, {}, {}
    for specimen, number, flag, spp, binned in zip(
            frame["sample_id"], frame["num_ticks"], frame["is_single_tick"],
            frame["host_spp"], frame["host_bin"]):
        counts[specimen] = pd.to_numeric(pd.Series([number]),
                                         errors="coerce")[0]
        flags[specimen] = str(flag).strip().lower() == "true"
        text = str(spp).strip()
        if text.lower() in ("", "nan", "none"):
            text = str(binned).strip()
        hosts[specimen] = text
    return counts, flags, hosts


def specimen_of(table, module, row):
    """I resolve one qPCR row to a crosswalk specimen, or '' if it is absent."""
    for column in ("Barcode", "host_id"):
        hit = module.resolve(table, "specimen", row.get(column))
        if hit is not None:
            return str(hit)
    return ""


def positives(frame, target):
    """I count rows where the target is a positive call."""
    values = pd.to_numeric(frame[target], errors="coerce").fillna(0)
    return int((values > 0).sum())


def report(frame, label):
    """I print one host by target table and return it."""
    print("=" * 78)
    print(label.upper())
    print("=" * 78)
    if frame.empty:
        print("  no rows in this cohort.")
        print()
        return None
    rows = []
    order = [name for name, _ in HOST_GROUPS] + [UNMAPPED, ""]
    groups = [g for g in order if g in set(frame["host_group"])]
    header = "  %-20s %5s %6s" % ("host", "n", "ticks")
    for target in TARGETS:
        header += " %13s" % target[:13]
    print(header)
    for group in groups:
        part = frame[frame["host_group"] == group]
        ticks = pd.to_numeric(part["num_ticks"], errors="coerce").sum()
        line = "  %-20s %5d %6s" % ((group or "unresolved")[:20], len(part),
                                    "%.0f" % ticks if pd.notna(ticks) else "?")
        record = {"cohort": label, "host": group or "unresolved",
                  "n": len(part),
                  "ticks": None if pd.isna(ticks) else int(ticks)}
        for target in TARGETS:
            count = positives(part, target)
            record[target] = count
            record[target + "_pct"] = round(100.0 * count / len(part), 1)
            line += " %13s" % ("%d (%.0f%%)"
                               % (count, 100.0 * count / len(part)))
        rows.append(record)
        print(line)
    all_ticks = pd.to_numeric(frame["num_ticks"], errors="coerce").sum()
    total = "  %-20s %5d %6s" % ("ALL", len(frame),
                                 "%.0f" % all_ticks if pd.notna(all_ticks)
                                 else "?")
    for target in TARGETS:
        count = positives(frame, target)
        total += " %13s" % ("%d (%.0f%%)"
                            % (count, 100.0 * count / len(frame)))
    print(total)
    print()
    return pd.DataFrame(rows)


def main():
    write = "--write" in sys.argv[1:]
    table, module = lookup()
    qpcr = read_qpcr()
    counts, flags, hosts = read_crosswalk()

    specimens = [specimen_of(table, module, row) for _, row in qpcr.iterrows()]
    qpcr = qpcr.copy()
    qpcr["specimen"] = specimens
    qpcr["num_ticks"] = [counts.get(s, float("nan")) for s in specimens]
    qpcr["single"] = [bool(flags.get(s, False)) for s in specimens]
    qpcr["host_text"] = [hosts.get(s, "") for s in specimens]
    qpcr["host_group"] = [host_group(t) for t in qpcr["host_text"]]
    qpcr["subsamples"] = pd.to_numeric(qpcr["subsample_count"],
                                       errors="coerce")

    missing = [s for s in specimens if not s]
    print("qPCR rows: %d; resolved to a crosswalk specimen: %d; not in the "
          "crosswalk: %d" % (len(qpcr), len(qpcr) - len(missing), len(missing)))
    print("rows with a recognised host: %d"
          % sum(1 for g in qpcr["host_group"] if named(g)))

    unmapped = qpcr.loc[~qpcr["host_group"].map(named), "host_text"]
    if len(unmapped):
        print()
        print("host values I could not map (%d rows). Tell me which of these "
              "are real hosts and I will add them:" % len(unmapped))
        for value, count in unmapped.astype(str).value_counts().head(SHOW).items():
            print("  %-46s %4d" % (value[:46] or "(blank)", count))
    print()

    scored = sum((pd.to_numeric(qpcr[t], errors="coerce").fillna(0) > 0)
                 .astype(int) for t in TARGETS)
    controls = qpcr[scored >= CONTROL_POSITIVES]
    if len(controls):
        print("rows positive for %d or more of the %d targets, which is a "
              "positive control and not a tick:"
              % (CONTROL_POSITIVES, len(TARGETS)))
        for specimen, text, single in zip(controls["specimen"],
                                          controls["host_text"],
                                          controls["single"]):
            print("  %-28s %-20s single tick: %s"
                  % (str(specimen)[:28] or "(no specimen)", str(text)[:20],
                     "yes" if single else "no"))
        print("  I exclude them below.")
        qpcr = qpcr[scored < CONTROL_POSITIVES]
        print()

    pooled = qpcr[qpcr["subsamples"] == 1]
    # One animal, one result. The qPCR sheet holds TB-90024065 twice, once as
    # '4065' and once as 'tb-90024065', and tick_lookup resolves both to the
    # same moose. Rows with no specimen are left alone: they are unresolved,
    # not known to be one animal.
    names = pooled["specimen"].astype(str).str.strip().str.lower()
    has_name = ~names.isin(("", "nan", "none"))
    repeated = has_name & pooled.duplicated(subset="specimen", keep="first")
    if repeated.any():
        print("animals appearing in more than one row: %s"
              % ", ".join(sorted(set(pooled.loc[repeated, "specimen"]))))
        print("  keeping the first row of each; %d rows dropped"
              % int(repeated.sum()))
        print()
        pooled = pooled[~repeated]
    single = pooled[pooled["single"]]
    pools = pooled[(~pooled["single"]) & (pooled["num_ticks"] > 1)]
    one_tick_not_flagged = pooled[(~pooled["single"])
                                  & (pooled["num_ticks"] == 1)]
    unknown = int(pooled["num_ticks"].isna().sum())
    print("one subsample tested: %d rows" % len(pooled))
    print("  single ticks (is_single_tick):                %d" % len(single))
    print("  pools (more than one tick):                   %d" % len(pools))
    print("  one-tick tubes the flag excludes:             %d"
          % len(one_tick_not_flagged))
    print("  no tick count in the crosswalk:               %d" % unknown)
    if len(one_tick_not_flagged):
        excluded = sorted(set(one_tick_not_flagged["specimen"].dropna()))
        print("  excluded by the flag: %s" % ", ".join(excluded)[:200])
    print()

    parts = []
    for frame, label in ((single, "single ticks, one tick per assay"),
                         (pools, "pools, one slurry per assay"),
                         (qpcr, "all qPCR rows")):
        result = report(frame, label)
        if result is not None:
            parts.append(result)

    print("The single-tick row above should read 31, 20, 13, 1, 2, 0, 0 of 89")
    print("under the criteria of Section 2.5. It read 31, 20, 12, 1, 1, 0, 0")
    print("from the supplied pass_check column. The difference is one")
    print("B. microti animal and one Powassan virus animal, both promoted on")
    print("2026-10-08 by build_corrected_screening_record.py.")
    print()
    print("Do not check against results/analysis/qpcr_prevalence.tsv. That is")
    print("retired output of build_qpcr_figure.py, whose figure was archived")
    print("to results/_archive/2026-09-24_restructure and is not in the")
    print("paper. Its denominator is 89, though paragraph 95 describes the")
    print("retired figure as using 93, which is unresolved and does not")
    print("affect anything here.")
    print()

    if not write:
        print("Nothing written. Rerun with --write to save %s" % OUT)
        return
    if os.path.exists(OUT):
        sys.exit("%s already exists. Move it aside first." % OUT)
    pd.concat(parts, ignore_index=True).to_csv(OUT, sep="\t", index=False)
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
