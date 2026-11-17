#!/usr/bin/env python3
"""check_batch_arms.py

I ask what the two "arms" in Figure 1B are made of, and whether Batch 3 (which
supplies most of them) is the P3-5 plate. Read-only, except that --write adds
results/diversity_rebuild_20260706_2045/meta_wtd_one_per_animal_no_batch3.tsv
so the batch model can be rerun without Batch 3.

Groups use the same ordination as fig1b_deer_batch_pcoa.py:
    main PCoA1 < 0.1;  upper PCoA1 >= 0.1, PCoA2 >= 0.1;
    lower PCoA1 >= 0.1, PCoA2 <= -0.15;  other: the rest.

USAGE
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/check_batch_arms.py [--write]
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
import datetime as _dt

RUN = "results/diversity_rebuild_20260706_2045"
MATRIX = RUN + "/dm_bray_curtis_one_per_animal.qza"
META = RUN + "/meta_wtd_one_per_animal.tsv"
META_FULL = RUN + "/meta_wtd_clean.tsv"
RAREFIED = RUN + "/core_metrics/rarefied_table.qza"
ANALYSIS = "results/analysis/analysis_table_16S.tsv"
TAXONOMY_CANDIDATES = [
    "results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv",
    "results/16S/all/taxonomy_verified/taxonomy_16S_blast_verified.tsv",
]
OUT_NO_B3 = RUN + "/meta_wtd_one_per_animal_no_batch3.tsv"
PLATE = re.compile(r"P3[-_ ]?5", re.IGNORECASE)
TOP = 6


def extract(qza, suffix):
    """I pull one file out of a QIIME artifact into a temp dir; caller cleans."""
    if not os.path.isfile(qza):
        sys.exit("missing %s" % qza)
    workdir = tempfile.mkdtemp(prefix="arms_")
    with zipfile.ZipFile(qza) as archive:
        wanted = [n for n in archive.namelist() if n.endswith(suffix)]
        if not wanted:
            shutil.rmtree(workdir, ignore_errors=True)
            sys.exit("no %s inside %s" % (suffix, qza))
        archive.extract(wanted[0], workdir)
    return os.path.join(workdir, wanted[0]), workdir


def read_matrix():
    path, workdir = extract(MATRIX, "data/distance-matrix.tsv")
    try:
        frame = pd.read_csv(path, sep="\t", index_col=0)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.loc[frame.index, frame.index]


def principal_coordinates(distances):
    """The same classical scaling as fig1b_deer_batch_pcoa.py."""
    d = distances.values.astype(float)
    n = d.shape[0]
    centre = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * centre @ (d ** 2) @ centre
    values, vectors = np.linalg.eigh(b)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    return vectors[:, :2] * np.sqrt(np.maximum(values[:2], 0))


def group_of(x, y):
    if x < 0.1:
        return "main"
    if y >= 0.1:
        return "upper"
    if y <= -0.15:
        return "lower"
    return "other"


def read_rarefied():
    """I return the rarefied table as samples x features (proportions)."""
    try:
        import biom
    except ImportError:
        sys.exit("biom is not importable. Run inside the qiime2 environment.")
    path, workdir = extract(RAREFIED, "data/feature-table.biom")
    try:
        table = biom.load_table(path)
        frame = pd.DataFrame(table.matrix_data.toarray().T,
                             index=[str(s) for s in table.ids(axis="sample")],
                             columns=[str(f) for f in table.ids(axis="observation")])
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    totals = frame.sum(axis=1).replace(0, np.nan)
    return frame.div(totals, axis=0).fillna(0)


def genus_label(taxon):
    """I return the genus from a SILVA string, or the deepest named rank."""
    parts = [p.strip() for p in str(taxon).split(";") if p.strip()]
    named = []
    for part in parts:
        rank, _, name = part.partition("__")
        if name and name.lower() not in ("uncultured", "unassigned", ""):
            named.append((rank, name))
    for rank, name in named:
        if rank == "g":
            return name
    if named:
        rank, name = named[-1]
        return "%s (%s)" % (name, {"f": "family", "o": "order", "c": "class",
                                   "p": "phylum", "d": "domain"}.get(rank, rank))
    return "Unassigned"


def read_taxonomy(features):
    """I pick the taxonomy file that covers the most features, and say so."""
    best = None
    for path in TAXONOMY_CANDIDATES + sorted(glob.glob("results/**/taxonomy*16S*.tsv",
                                                        recursive=True)):
        if not os.path.isfile(path) or ".bak" in path or "_archive" in path:
            continue
        try:
            frame = pd.read_csv(path, sep="\t", dtype=str)
        except (ValueError, OSError, pd.errors.ParserError):
            continue
        id_col = frame.columns[0]
        tax_col = next((c for c in frame.columns if c.lower() == "taxon"), None)
        if tax_col is None:
            continue
        mapping = dict(zip(frame[id_col], frame[tax_col]))
        covered = sum(1 for f in features if f in mapping)
        if best is None or covered > best[1]:
            best = (path, covered, mapping)
        if covered == len(features):
            break
    try:
        os.makedirs("results/16S/all/taxonomy_verified", exist_ok=True)
        with open("results/16S/all/taxonomy_verified/TAXONOMY_USED.txt", "a") as _fh:
            _fh.write("%s\t%s\t%d of %d features\n" % (
                _dt.datetime.now().isoformat(timespec="seconds"),
                best[0] if best else "NONE",
                best[1] if best else 0, len(features)))
    except OSError:
        pass
    if best is None or best[1] == 0:
        sys.exit("no 16S taxonomy table matches the rarefied table's features.")
    return best


def plate_hits(samples):
    """I return {sample: where P3-5 was found} from names and metadata."""
    hits = {}
    for sample in samples:
        if PLATE.search(sample):
            hits[sample] = "sample name"
    for path, key in ((META_FULL, None), (ANALYSIS, "sample")):
        if not os.path.isfile(path):
            continue
        frame = pd.read_csv(path, sep="\t", dtype=str)
        frame = frame[~frame[frame.columns[0]].astype(str).str.startswith("#q2")]
        ids = frame[key] if key else frame[frame.columns[0]]
        for column in frame.columns:
            matches = frame[column].astype(str).str.contains(PLATE, na=False)
            for sample in ids[matches]:
                sample = str(sample)
                for candidate in samples:
                    if candidate == sample or candidate.split(".")[0] == sample.split(".")[0]:
                        hits.setdefault(candidate, "%s:%s" % (os.path.basename(path), column))
    return hits


def main():
    write = "--write" in sys.argv[1:]
    distances = read_matrix()
    meta = pd.read_csv(META, sep="\t", dtype=str)
    meta = meta[meta["sample-id"] != "#q2:types"].set_index("sample-id")
    samples = list(distances.index)
    coords = principal_coordinates(distances)
    work = pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1]}, index=samples)
    work["group"] = [group_of(x, y) for x, y in zip(work["x"], work["y"])]
    work["batch"] = meta.loc[samples, "batch_clean"].values
    work["pool"] = pd.to_numeric(meta.loc[samples, "pool_size"], errors="coerce").values
    short = {b: re.sub(r"\s*\(.*\)", "", b).replace("Batch ", "B")
             for b in work["batch"].unique()}

    shares = read_rarefied()
    missing = [s for s in samples if s not in shares.index]
    if missing:
        sys.exit("%d deer are not in the rarefied table, e.g. %s"
                 % (len(missing), missing[:3]))
    shares = shares.loc[samples]
    shares = shares.loc[:, shares.sum(axis=0) > 0]
    path, covered, mapping = read_taxonomy(list(shares.columns))
    genera = shares.T.groupby([genus_label(mapping.get(f, "")) for f in shares.columns]).sum().T

    print("=" * 72)
    print("1. WHAT THE GROUPS ARE MADE OF (%d deer, rarefied reads)" % len(samples))
    print("=" * 72)
    print("  taxonomy: %s (%d of %d features)" % (path, covered, shares.shape[1]))
    for name in ("main", "upper", "lower", "other"):
        rows = work[work["group"] == name]
        if rows.empty:
            continue
        counts = rows["batch"].value_counts()
        mix = ", ".join("%s %d" % (short[b], counts[b])
                        for b in sorted(counts.index, key=lambda b: short[b]))
        mean = genera.loc[rows.index].mean().sort_values(ascending=False)
        rick = mean.get("Rickettsia", 0.0)
        print()
        print("  %-5s n = %3d   median pool %s" % (name.upper(), len(rows),
              "%.0f" % rows["pool"].median() if rows["pool"].notna().any() else "-"))
        print("        batches: %s" % mix)
        print("        Rickettsia: %.1f%% of reads" % (100 * rick))
        top = [(g, v) for g, v in mean.items() if g != "Rickettsia"][:TOP]
        print("        next: " + "; ".join("%s %.1f%%" % (g, 100 * v) for g, v in top))
    print()

    print("=" * 72)
    print("2. HOW MUCH OF EACH BATCH LEAVES THE MAIN CLUSTER")
    print("=" * 72)
    for batch in sorted(work["batch"].unique(), key=lambda b: short[b]):
        rows = work[work["batch"] == batch]
        off = int((rows["group"] != "main").sum())
        print("  %-22s %3d of %3d off the main cluster (%.0f%%)"
              % (batch[:22], off, len(rows), 100.0 * off / len(rows)))
    print()

    print("=" * 72)
    print("3. WHERE IS THE P3-5 PLATE?")
    print("=" * 72)
    hits = plate_hits(samples)
    if not hits:
        print("  no deer library carries 'P3-5' in its name, the diversity")
        print("  metadata or the analysis table. Plate membership cannot be read")
        print("  from these files; the run sheet would settle it.")
    else:
        found = work.loc[list(hits)]
        print("  %d of the %d deer libraries are marked P3-5 (first found in %s)"
              % (len(hits), len(samples), sorted(set(hits.values()))[0]))
        for batch, rows in found.groupby("batch"):
            groups = rows["group"].value_counts().to_dict()
            print("    %-22s %3d   %s" % (batch[:22], len(rows), groups))
        b3 = work[work["batch"].str.startswith("Batch 3")]
        if len(b3):
            in_plate = b3.index.isin(list(hits))
            off = b3["group"] != "main"
            print("  within Batch 3: off-cluster %d of %d on P3-5, %d of %d elsewhere"
                  % (int((off & in_plate).sum()), int(in_plate.sum()),
                     int((off & ~in_plate).sum()), int((~in_plate).sum())))
    print()

    if not write:
        print("Nothing written. Rerun with --write to write the set without Batch 3.")
        return
    if os.path.exists(OUT_NO_B3):
        sys.exit("%s already exists. Move it aside first." % OUT_NO_B3)
    full = pd.read_csv(META, sep="\t", dtype=str)
    keep = (full["sample-id"] == "#q2:types") | ~full["batch_clean"].astype(str).str.startswith("Batch 3")
    full[keep].to_csv(OUT_NO_B3, sep="\t", index=False)
    print("wrote %s (%d deer)" % (OUT_NO_B3, int(keep.sum()) - 1))


if __name__ == "__main__":
    main()
