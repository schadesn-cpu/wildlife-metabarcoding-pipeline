#!/usr/bin/env python3
"""prepare_wtd_diversity_subsets.py

I build the two sample sets for redoing the deer diversity analysis, and report
how they spread across batches BEFORE anything is run, so we know whether each
set can estimate a batch effect at all.

A  ONE LIBRARY PER DEER (deepest library kept; pool size carried for the model).
B  SINGLE TICKS from deer (project single-tick table if found, else tick count 1).

I write nothing unless given --write, and never over an existing file.

USAGE
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/prepare_wtd_diversity_subsets.py           # report only
    python3 scripts/prepare_wtd_diversity_subsets.py --write   # also write sets
"""
import glob
import os
import re
import shutil
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd

RUN = "results/diversity_rebuild_20260706_2045"
META = RUN + "/meta_wtd_clean.tsv"
BRAY = RUN + "/core_metrics/bray_curtis_distance_matrix.qza"
SHANNON = RUN + "/core_metrics/shannon_vector.qza"
TABLE = RUN + "/table_wtd_decontam.qza"
ANALYSIS = "results/analysis/analysis_table_16S.tsv"
OUT_A = RUN + "/meta_wtd_one_per_animal.tsv"
OUT_B = RUN + "/meta_wtd_single_tick.tsv"
SINGLE_PATTERNS = ["results/**/*single_tick*.tsv", "results/**/*single*tick*.tsv"]
RUN_SUFFIX = re.compile(r"\..*$")


def blank(value):
    """I return '' for anything meaning 'no value', and the text otherwise."""
    text = str(value).strip()
    return "" if text.lower() in ("nan", "none", "<na>", "na", "") else text


def extract(qza, pattern):
    """I pull one file out of a QIIME artifact into a temp dir; caller cleans."""
    workdir = tempfile.mkdtemp(prefix="qza_")
    with zipfile.ZipFile(qza) as archive:
        wanted = [n for n in archive.namelist() if re.search(pattern, n)]
        if not wanted:
            shutil.rmtree(workdir, ignore_errors=True)
            return None, None
        archive.extract(wanted[0], workdir)
    return os.path.join(workdir, wanted[0]), workdir


def matrix_ids(qza):
    """I return the sample IDs in a distance matrix, reading only its header."""
    path, workdir = extract(qza, r"data/distance-matrix\.tsv$")
    if path is None:
        sys.exit("no distance-matrix.tsv inside %s" % qza)
    try:
        with open(path) as handle:
            header = handle.readline().rstrip("\n").split("\t")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return [h for h in header[1:] if h]


def shannon_values():
    """I return the Shannon vector as a Series, or None."""
    if not os.path.isfile(SHANNON):
        return None
    path, workdir = extract(SHANNON, r"data/alpha-diversity\.tsv$")
    if path is None:
        return None
    try:
        frame = pd.read_csv(path, sep="\t", index_col=0)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return frame.iloc[:, 0]


def library_depths():
    """I return reads per library before rarefaction, or {} if I cannot."""
    try:
        import biom
    except ImportError:
        print("  (biom not importable; I will break ties by name, not depth)")
        return {}
    path, workdir = extract(TABLE, r"data/feature-table\.biom$")
    if path is None:
        return {}
    try:
        table = biom.load_table(path)
        return dict(zip(table.ids(axis="sample"), table.sum(axis="sample")))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def column(frame, *names):
    """I return the first column matching any name, ignoring case."""
    lowered = {str(c).lower(): c for c in frame.columns}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return None


def load_lookup():
    try:
        import tick_lookup
    except ImportError as problem:
        sys.exit("cannot import tick_lookup (%s). Run from the project root."
                 % problem)
    return tick_lookup.load(), tick_lookup


def resolve(table, module, value):
    stem = RUN_SUFFIX.sub("", str(value).strip())
    for marker, candidate in (("16S", value), ("16S", stem),
                              ("specimen", stem), ("specimen", value)):
        hit = module.resolve(table, marker, candidate)
        if hit is not None:
            return blank(hit)
    return ""


def pool_sizes(table, module):
    """I return the tick count per animal from the analysis table."""
    if not os.path.isfile(ANALYSIS):
        return {}
    frame = pd.read_csv(ANALYSIS, sep="\t", dtype=str)
    if "NumTicks" not in frame.columns or "sample" not in frame.columns:
        return {}
    sizes = {}
    for sample, ticks in zip(frame["sample"], frame["NumTicks"]):
        animal = resolve(table, module, sample)
        value = pd.to_numeric(ticks, errors="coerce")
        if animal and not pd.isna(value):
            sizes[animal] = max(sizes.get(animal, 0), float(value))
    return sizes


def official_single_ticks(table, module):
    """I look for the project's single-tick table and return its animals."""
    candidates = []
    for pattern in SINGLE_PATTERNS:
        candidates.extend(glob.glob(pattern, recursive=True))
    best = None
    for path in sorted(set(candidates)):
        if ".bak" in path or "_archive" in path or "stale" in path.lower():
            continue
        try:
            frame = pd.read_csv(path, sep="\t", dtype=str)
        except (ValueError, OSError, pd.errors.ParserError):
            continue
        if not 80 <= len(frame) <= 110:
            continue
        for name in frame.columns:
            animals = set(resolve(table, module, v) for v in frame[name])
            animals.discard("")
            if len(animals) >= 0.9 * len(frame):
                if best is None or len(animals) > len(best[1]):
                    best = (path, animals)
                break
    return best


def main():
    write = "--write" in sys.argv[1:]
    for path in (META, BRAY):
        if not os.path.isfile(path):
            sys.exit("missing %s. Run from the project root." % path)
    meta = pd.read_csv(META, sep="\t", dtype=str)
    id_column = meta.columns[0]
    meta = meta[~meta[id_column].astype(str).str.startswith("#q2")]
    batch = column(meta, "batch_clean", "batch")
    state = column(meta, "state", "State")
    if batch is None or state is None:
        sys.exit("need batch and state columns in %s; it has: %s"
                 % (META, ", ".join(meta.columns[:15])))
    meta = meta.set_index(id_column)

    ids = matrix_ids(BRAY)
    missing = [i for i in ids if i not in meta.index]
    if missing:
        sys.exit("%d matrix IDs are not in %s, e.g. %s"
                 % (len(missing), META, missing[:3]))
    work = meta.loc[ids, [batch, state]].copy()
    work.columns = ["batch", "state"]

    table, module = load_lookup()
    work["animal"] = [resolve(table, module, i) for i in work.index]
    unresolved = int((work["animal"] == "").sum())
    sizes = pool_sizes(table, module)
    work["pool_size"] = [sizes.get(a, np.nan) for a in work["animal"]]
    depths = library_depths()
    work["depth"] = [float(depths.get(i, 0)) for i in work.index]

    print("=" * 72)
    print("1. THE JULY DEER MATRIX")
    print("=" * 72)
    print("  %d libraries, %d deer, %d libraries not resolved to an animal"
          % (len(work), work.loc[work["animal"] != "", "animal"].nunique(),
             unresolved))
    print("  pool size known for %d of %d libraries"
          % (work["pool_size"].notna().sum(), len(work)))
    print()

    usable = work[work["animal"] != ""].copy()
    usable = usable.sort_values(["animal", "depth"], ascending=[True, False])
    every_deer = usable.drop_duplicates("animal", keep="first")
    no_pool = int((~(every_deer["pool_size"] > 0)).sum())
    set_a = every_deer[every_deer["pool_size"] > 0]
    conflicts = int((usable.groupby("animal")["batch"].nunique() > 1).sum())

    official = official_single_ticks(table, module)
    if official is not None:
        path, singles = official
        set_b = every_deer[every_deer["animal"].isin(singles)]
        source = "project table %s (%d animals)" % (path, len(singles))
    else:
        set_b = every_deer[every_deer["pool_size"] == 1]
        source = "tick count = 1 (no single-tick table found)"

    print("=" * 72)
    print("2. HOW THE SETS SPREAD ACROSS BATCHES")
    print("=" * 72)
    print("  set A: one library per deer (deepest)      n = %d" % len(set_a))
    print("  set B: single ticks, %s" % source)
    print("         deer single ticks in the matrix     n = %d" % len(set_b))
    print("  deer whose libraries sit in more than one batch: %d" % conflicts)
    print("  deer left out of set A for having no pool size: %d" % no_pool)
    print()
    print("  %-26s %6s %5s %6s %6s %9s" % ("batch", "libs", "deer", "set A",
                                          "set B", "med pool"))
    for name, rows in work.groupby("batch"):
        in_a = set_a[set_a["batch"] == name]
        median = ("%.0f" % in_a["pool_size"].median()
                  if in_a["pool_size"].notna().any() else "-")
        print("  %-26s %6d %5d %6d %6d %9s"
              % (str(name)[:26], len(rows), rows["animal"].nunique(),
                 len(in_a), int((set_b["batch"] == name).sum()), median))
    print()
    small = [n for n, c in set_b["batch"].value_counts().items() if c < 5]
    if small or set_b["batch"].nunique() < 3:
        print("  WARNING: set B has %d batches, %d with fewer than 5 ticks."
              % (set_b["batch"].nunique(), len(small)))
        print("  Batch is weakly estimable in set B; read it as a check on A.")
        print()

    print("=" * 72)
    print("3. IS THE SHANNON PILE-UP AT 1.0 IN THE DATA?")
    print("=" * 72)
    shannon = shannon_values()
    if shannon is None:
        print("  could not read %s" % SHANNON)
    else:
        deer = shannon.reindex(work.index).dropna()
        print("  %d deer values: min %.3f, 5%% %.3f, median %.3f, max %.3f"
              % (len(deer), deer.min(), deer.quantile(0.05), deer.median(),
                 deer.max()))
        print("  below 1.0: %d    between 0.99 and 1.01: %d"
              % ((deer < 1).sum(), ((deer > 0.99) & (deer < 1.01)).sum()))
    print()

    if not write:
        print("Nothing written. Rerun with --write to write the two sets.")
        return
    counts = set_b["batch"].value_counts()
    kept = counts[counts >= 5].index
    print("set B keeps batches with >=5 ticks: %d ticks in %d batches"
          % (int(set_b["batch"].isin(kept).sum()), len(kept)))
    set_b = set_b[set_b["batch"].isin(kept)]
    for target, frame in ((OUT_A, set_a), (OUT_B, set_b)):
        if os.path.exists(target):
            sys.exit("%s already exists. Move it aside first." % target)
        out = pd.DataFrame({
            "sample-id": frame.index,
            "batch_clean": frame["batch"].values,
            "state": frame["state"].values,
            "pool_size": frame["pool_size"].values,
            "log2_pool_size": np.log2(frame["pool_size"].values),
        })
        types = pd.DataFrame([["#q2:types", "categorical", "categorical",
                               "numeric", "numeric"]], columns=out.columns)
        pd.concat([types, out]).to_csv(target, sep="\t", index=False,
                                       na_rep="")
        print("wrote %s (%d samples)" % (target, len(out)))


if __name__ == "__main__":
    main()
