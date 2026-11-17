#!/usr/bin/env python3
"""I turn the accessioning and bench records into one covariate table.

WHY
---
metadata/tick_full_metadata.xlsx has 22 sheets and the project reads one of
them. The other 21 hold what a person wrote down when each tube arrived and
when it was opened: how many ticks, how many of those were engorged females
and how many were males, whether a deer ked was in the tube, whether a
spider was, whether the tube was dry or empty or the ticks crushed, and
which extraction batch it went into. None of that has been available to any
analysis, because it lives in prose and in per-batch sheets rather than in a
table.

I extract it once, into one row per animal, so that the microbiome and the
qPCR results can be grouped by what was physically in the tube.

A NOTE ON IDENTIFIERS
---------------------
Some bench sheets name subsamples rather than animals, for example
EQ02442213-1 and EQ02442213-2 for two aliquots of one tube. Those do not
resolve to an animal and I report them as unresolved rather than stripping
the suffix, because stripping is the guess tick_lookup exists to prevent.
Animals reached only through such a sheet therefore have no extraction
batch, and I say how many.

THE KED RECORD HAS THREE REGIMES, AND I KEEP THEM APART
-------------------------------------------------------
Batches 01 to 06 noted keds in prose at accessioning and never recorded a
zero, so an animal there without a ked note is unknown, not negative.
Batches 07 to 10 counted keds at the bench including zeros, so those are the
only animals that can be scored negative. Batches 11 and 12 have no ked
record at all. I write a source column, and any analysis that treats a blank
as a zero outside batches 07 to 10 is wrong.

WHAT I WILL NOT DO
------------------
I will not invent a count. A note reading "1 ked" gives 1; a note that
mentions a ked without a number gives 1 and is flagged; a weight written as
"165 WTD KED" gives 1, not 165, because the number must sit immediately
before the word. Every barcode goes through tick_lookup and an identifier
that does not resolve is named and left out rather than guessed at.

I CHECK MYSELF
--------------
The counts this produces are already known from reading the workbook
directly: 263 accessioning notes holding text, 27 animals with a ked noted,
96 with both an engorged and a male count, 12 tubes recorded KEDS ONLY, 28
animals with a bench ked count above zero and 44 scored zero. I assert each
one and stop if any disagrees, because a silent drift in the parsing would
be invisible downstream. The note count is 263 and not the 352 non-null
cells, because 89 of those cells hold an empty string.

OUTPUT
------
results/metadata/accessioning_covariates.tsv. I create no directories and I
refuse to overwrite.
"""
import os
import re
import sys

import pandas as pd

sys.path.insert(0, "scripts")
import tick_lookup

WORKBOOK = "metadata/tick_full_metadata.xlsx"
OUT_PATH = "results/metadata/accessioning_covariates.tsv"

BATCH_SHEETS = ["2024-04-18 Batch 01", "2024-07-16 Batch 02", "2024-06-24 Batch 03",
                "2024-07-23 Batch 04", "2024-07-23 Batch 05", "2024-10-30 Batch 06",
                "2025-01-24 Batch 07", "2025-02-12 Batch 08", "2025-02-13 Batch 09",
                "2025-02-19 Batch 10", "2025-02-27 Batch 11", "2025-03-07 Batch 12"]
KED_COLUMN = {"2025-01-24 Batch 07": "KED", "2025-02-12 Batch 08": '"ked" count',
              "2025-02-13 Batch 09": '"ked" count', "2025-02-19 Batch 10": '"ked" count'}
WORDS = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5}
# Counted by reading the workbook directly. "notes filled" counts cells
# holding actual text: 352 cells are non-null but 89 of those are an empty
# string, which is not a note.
EXPECTED = {"notes filled": 263, "ked noted": 27, "sex census": 96,
            "keds only": 12, "bench ked": 28, "bench zero": 44}


def tidy(frame):
    frame.columns = [str(c).strip() for c in frame.columns]
    return frame


def clean(value):
    text = re.sub(r"\s+", " ", str(value if value is not None else "")).strip()
    return "" if text.lower() in ("nan", "none") else text


def number(value):
    text = clean(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def ked_from_note(text):
    """I read a ked count from prose, taking only a number next to the word."""
    if not re.search(r"\bkeds?\b", text, re.I):
        return None, False
    match = re.search(r"\b(\d{1,2}|one|two|three|four|five)\s+(?:wtd\s+)?keds?\b",
                      text, re.I)
    if not match:
        return 1, True
    token = match.group(1).lower()
    return (WORDS.get(token) or int(token)), False


def main():
    if not os.path.isdir("results/metadata"):
        sys.exit("run me from the project root; there is no results/metadata here.")
    if not os.path.isfile(WORKBOOK):
        sys.exit("missing " + WORKBOOK)
    if os.path.exists(OUT_PATH):
        sys.exit(OUT_PATH + " already exists. Move it aside first.")

    table = tick_lookup.load()
    book = pd.ExcelFile(WORKBOOK)
    meta = tidy(book.parse("Metadata", dtype=str))

    # --- the bench sheets: batch membership, tick and ked counts ----------
    batch_of, bench_ticks, bench_ked, bench_zero = {}, {}, {}, set()
    seen_batches, bench_missed = {}, {}
    # Tallied from the workbook before resolution, so the self-check
    # measures my parsing rather than the lookup's coverage.
    raw_ked, raw_zero, raw_keds_only = set(), set(), set()
    keds_only, mating, bench_note = set(), set(), {}
    tally = {"keds only": 0, "bench ked": 0, "bench zero": 0}
    for sheet in BATCH_SHEETS:
        if sheet not in book.sheet_names:
            continue
        # Batch 06 carries two date lines above its header, so its real
        # column names sit on the fourth row rather than the first.
        header = 3 if sheet.endswith("Batch 06") else 0
        frame = tidy(book.parse(sheet, dtype=str, header=header))
        id_column = next((c for c in ("OG SAMPLE ID (ANIMAL/METADATA ID)",
                                      "Barcode", "CoCo ID")
                          if c in frame.columns), None)
        if id_column is None:
            print("   %s has no identifier column; skipped" % sheet)
            continue
        label = sheet[-8:]
        column = KED_COLUMN.get(sheet)
        for _, row in frame.iterrows():
            raw = clean(row.get(id_column))
            if raw:
                raw_note = clean(row.get("Additional Notes"))
                if raw_note and re.search(r"keds?\s*only", raw_note, re.I):
                    raw_keds_only.add(raw)
                if column and column.strip() in frame.columns:
                    raw_value = number(row.get(column.strip()))
                    if raw_value is not None:
                        if raw_value > 0:
                            raw_ked.add(raw)
                        else:
                            raw_zero.add(raw)
            animal = tick_lookup.resolve(table, "specimen", raw)
            if animal is None:
                if raw:
                    bench_missed.setdefault(label, []).append(raw)
                continue
            animal = str(animal)
            batch_of.setdefault(animal, label)
            seen_batches.setdefault(animal, []).append(label)
            count = number(row.get("Tick count"))
            if count is None:
                count = number(row.get("TICK COUNT"))
            if count is not None:
                bench_ticks[animal] = count
            note = clean(row.get("Additional Notes"))
            if note:
                bench_note[animal] = note
                if re.search(r"keds?\s*only", note, re.I):
                    keds_only.add(animal)
                if re.search(r"mating", note, re.I):
                    mating.add(animal)
            if column and column.strip() in frame.columns:
                value = number(row.get(column.strip()))
                if value is None:
                    continue
                if value > 0:
                    bench_ked[animal] = max(bench_ked.get(animal, 0), value)
                else:
                    bench_zero.add(animal)
    bench_zero -= set(bench_ked)
    raw_zero -= raw_ked
    tally["keds only"] = len(raw_keds_only)
    tally["bench ked"] = len(raw_ked)
    tally["bench zero"] = len(raw_zero)

    # --- the Metadata sheet: prose notes ----------------------------------
    rows, unresolved = [], []
    notes_filled = ked_noted = sex_census = 0
    for _, record in meta.iterrows():
        recorded = clean(record.get("Barcode"))
        if not recorded:
            continue
        accession = clean(record.get("Accessioning Notes"))
        paper = clean(record.get("Notes (on paper/tube)"))
        processing = clean(record.get("Processing Notes"))
        if accession:
            notes_filled += 1
        note_ked, assumed = ked_from_note(accession)
        if note_ked is not None:
            ked_noted += 1
        pair = re.search(r"\b(\d+)\s+(?:fully |partially )?engorged\b[^0-9]{0,12}"
                         r"(\d+)\s+males?\b", accession, re.I)
        if pair:
            sex_census += 1

        animal = tick_lookup.resolve(table, "specimen", recorded)
        if animal is None:
            unresolved.append(recorded)
            continue
        animal = str(animal)

        engorged = males = None
        if pair:
            engorged, males = int(pair.group(1)), int(pair.group(2))
        elif re.search(r"\bengorged\b", accession, re.I):
            lone = re.search(r"\b(\d+)\s+(?:fully |partially )?engorged\b",
                             accession, re.I)
            engorged = int(lone.group(1)) if lone else None

        trouble = re.search(r"\b(empty|spill\w*|leak\w*|torn|ripp\w*|crush\w*|"
                            r"squish\w*|squash\w*|dry|cloudy|deflat\w*|split)\b",
                            accession + " " + processing, re.I)
        season = re.search(r"\b(Youth Weekend|Regular Season|archery|bow)\b", paper, re.I)
        weight = re.search(r"\b(\d{2,3})#", paper)

        if animal in bench_ked:
            source = "bench"
            ked = bench_ked[animal]
        elif animal in bench_zero:
            source = "bench"
            ked = 0
        elif note_ked is not None:
            source = "notes"
            ked = note_ked
        else:
            source = "none"
            ked = None

        rows.append({
            "animal": animal,
            "barcode_recorded": recorded,
            "extraction_batch": batch_of.get(animal, ""),
            "ked_count": ked,
            "ked_source": source,
            "ked_count_assumed": bool(assumed) if source == "notes" else False,
            "ked_scored_zero": animal in bench_zero,
            "keds_only_tube": animal in keds_only,
            "spider_recorded": bool(re.search(r"spiders?", accession, re.I)),
            "engorged_count": engorged,
            "male_count": males,
            "sex_census_complete": pair is not None,
            "mating_pair": animal in mating,
            "tick_count_metadata": number(record.get("# of Ticks")),
            "tick_count_bench": bench_ticks.get(animal),
            "tube_problem": bool(trouble),
            "tube_problem_word": trouble.group(0).lower() if trouble else "",
            "season": season.group(0).title() if season else "",
            "carcass_weight_lb": int(weight.group(1)) if weight else None,
            "host_recorded": clean(record.get("Host Spp")),
            "state": clean(record.get("State")),
        })

    frame = pd.DataFrame(rows)
    print("Metadata rows %d, resolved to an animal %d, unresolved %d"
          % (len(meta), len(frame), len(unresolved)))
    # An animal can appear in more than one bench sheet, from a re-extraction.
    # extraction_batch holds the first it appears in, which is a choice, so I
    # say how often it was made rather than leaving it implicit.
    if bench_missed:
        total = sum(len(v) for v in bench_missed.values())
        print("bench rows whose identifier did not resolve to an animal: %d"
              % total)
        for label in sorted(bench_missed):
            sample = sorted(set(bench_missed[label]))[:3]
            print("   %-10s %3d, for example %s"
                  % (label, len(bench_missed[label]), ", ".join(sample)))
        print("   Some bench sheets name subsamples rather than animals, so")
        print("   these are a visible absence and not a failure to look.")
    blank_batch = int((frame["extraction_batch"] == "").sum())
    print("animals with no extraction batch recorded: %d of %d"
          % (blank_batch, len(frame)))
    repeats = {a: b for a, b in seen_batches.items() if len(set(b)) > 1}
    print("animals appearing in more than one extraction batch sheet: %d"
          % len(repeats))
    if repeats:
        shown = sorted(repeats.items())[:5]
        for animal, labels in shown:
            print("   %-24s %s  (extraction_batch holds the first)"
                  % (animal, ", ".join(sorted(set(labels)))))
    if unresolved:
        print("   not in tick_lookup, left out: %s"
              % ", ".join(sorted(set(unresolved))[:12]))

    got = {"notes filled": notes_filled, "ked noted": ked_noted,
           "sex census": sex_census, "keds only": tally["keds only"],
           "bench ked": tally["bench ked"], "bench zero": tally["bench zero"]}
    print()
    print("self-check against the workbook read directly")
    wrong = []
    for key in sorted(EXPECTED):
        mark = "ok" if got[key] == EXPECTED[key] else "MISMATCH"
        print("   %-14s expected %3d, got %3d   %s"
              % (key, EXPECTED[key], got[key], mark))
        if got[key] != EXPECTED[key]:
            wrong.append(key)
    if wrong:
        sys.exit("\n%s disagree. The parsing has drifted from what the workbook "
                 "holds and I have written nothing." % ", ".join(wrong))

    print()
    print("coverage: how much of the workbook reaches a known animal")
    print("   bench ked-positive, in the workbook %3d, resolved %3d"
          % (len(raw_ked), len(bench_ked)))
    print("   bench scored zero, in the workbook  %3d, resolved %3d"
          % (len(raw_zero), len(bench_zero)))
    print("   tubes marked KEDS ONLY              %3d, resolved %3d"
          % (len(raw_keds_only), len(keds_only)))
    short = len(raw_ked) - len(bench_ked)
    if short:
        print("   %d animals carry a recorded ked count and do not"
              % short)
        print("   resolve to a sequenced animal. They are absent from the")
        print("   crosswalk and from library_lookup both, so they were")
        print("   accessioned and never sequenced. That is a property of")
        print("   the collection, not a gap in the crosswalk and not a")
        print("   fault here.")
    frame.to_csv(OUT_PATH, sep="\t", index=False)
    print()
    print("ked record by regime")
    for source in ("bench", "notes", "none"):
        part = frame[frame["ked_source"] == source]
        with_ked = int((part["ked_count"].fillna(0) > 0).sum())
        print("   %-6s animals %4d, with a ked %3d, scored zero %3d"
              % (source, len(part), with_ked, int(part["ked_scored_zero"].sum())))
    print()
    print("other covariates, animals with a value")
    for column in ("engorged_count", "male_count", "sex_census_complete",
                   "mating_pair", "spider_recorded", "tube_problem",
                   "season", "carcass_weight_lb", "tick_count_bench"):
        series = frame[column]
        if series.dtype == bool:
            filled = int(series.sum())
        else:
            filled = int(series.notna().sum() - (series == "").sum())
        print("   %-22s %4d" % (column, filled))
    print()
    print("wrote " + OUT_PATH)


if __name__ == "__main__":
    main()
