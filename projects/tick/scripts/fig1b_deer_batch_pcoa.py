#!/usr/bin/env python3
"""fig1b_deer_batch_pcoa.py

I draw Figure 1B: the deer tick microbiome ordination by sequencing batch, one
library per deer, as small multiples (every deer in grey, one batch in blue per
panel, cross at the batch centroid). Eight colours in one scatter cannot stay
distinguishable to colour-blind readers, so I facet instead.

I read results/diversity_rebuild_20260706_2045/dm_bray_curtis_one_per_animal.qza
and meta_wtd_one_per_animal.tsv, compute principal coordinates by classical
scaling (as qiime diversity pcoa does), and write
results/figures/fig1b_deer_batch_pcoa.png and .pdf. The old fig1b is untouched.

The statistics on the figure are the 2026-09-21 results. If the matrix does not
hold 343 deer I stop, because they would describe a different set.

USAGE
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/fig1b_deer_batch_pcoa.py
"""
import os
import re
import shutil
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN = "results/diversity_rebuild_20260706_2045"
MATRIX = RUN + "/dm_bray_curtis_one_per_animal.qza"
META = RUN + "/meta_wtd_one_per_animal.tsv"
OUT = "results/figures/fig1b_deer_batch_pcoa"
EXPECTED_N = 343
STATS = ("Bray-Curtis, one library per deer (n = 343).  "
         "Batch R\u00b2 = 13.8%, p = 0.001;  state R\u00b2 = 0.6%, p = 0.31;  "
         "dispersion differs between batches (PERMDISP F = 12.20, p = 0.001)")
HIGHLIGHT = "#2a78d6"
CONTEXT = "#d4d3cd"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#ebeae6"


def read_matrix(qza):
    """I return the distance matrix inside a QIIME artifact as a DataFrame."""
    if not os.path.isfile(qza):
        sys.exit("missing %s. Run the filter-distance-matrix step first." % qza)
    workdir = tempfile.mkdtemp(prefix="dm_")
    try:
        with zipfile.ZipFile(qza) as archive:
            wanted = [n for n in archive.namelist()
                      if n.endswith("data/distance-matrix.tsv")]
            if not wanted:
                sys.exit("no distance-matrix.tsv inside %s" % qza)
            archive.extract(wanted[0], workdir)
            frame = pd.read_csv(os.path.join(workdir, wanted[0]), sep="\t",
                                index_col=0)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.loc[frame.index, frame.index]


def principal_coordinates(distances):
    """I return the first two axes and each axis's share of variation."""
    d = distances.values.astype(float)
    n = d.shape[0]
    centre = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * centre @ (d ** 2) @ centre
    values, vectors = np.linalg.eigh(b)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    coords = vectors[:, :2] * np.sqrt(np.maximum(values[:2], 0))
    share = values[:2] / values[values > 0].sum()
    return coords, share


def batch_key(label):
    """I sort batches by their number."""
    found = re.search(r"(\d+)", str(label))
    return (int(found.group(1)) if found else 999, str(label))


def main():
    distances = read_matrix(MATRIX)
    meta = pd.read_csv(META, sep="\t", dtype=str)
    meta = meta[meta["sample-id"] != "#q2:types"].set_index("sample-id")
    missing = [i for i in distances.index if i not in meta.index]
    if missing:
        sys.exit("%d matrix samples are not in %s, e.g. %s"
                 % (len(missing), META, missing[:3]))
    if len(distances) != EXPECTED_N:
        sys.exit("the matrix has %d deer, not %d; update STATS first."
                 % (len(distances), EXPECTED_N))
    coords, share = principal_coordinates(distances)
    batch = meta.loc[distances.index, "batch_clean"].values
    batches = sorted(set(batch), key=batch_key)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
    columns = 4
    rows = int(np.ceil(len(batches) / float(columns)))
    fig, axes = plt.subplots(rows, columns, figsize=(10.5, 2.75 * rows + 0.9),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for index, name in enumerate(batches):
        ax = axes[index]
        inside = batch == name
        ax.scatter(coords[~inside, 0], coords[~inside, 1], s=9, color=CONTEXT,
                   linewidths=0, zorder=1)
        ax.scatter(coords[inside, 0], coords[inside, 1], s=16, color=HIGHLIGHT,
                   edgecolors="white", linewidths=0.5, zorder=2)
        ax.scatter([coords[inside, 0].mean()], [coords[inside, 1].mean()],
                   marker="+", s=110, color=INK, linewidths=1.8, zorder=3)
        ax.set_title("%s   n = %d" % (name, int(inside.sum())), fontsize=9,
                     color=INK, loc="left")
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(INK_SECONDARY)
            ax.spines[side].set_linewidth(0.6)
        ax.tick_params(colors=INK_SECONDARY, labelsize=8, length=2)
    for ax in axes[len(batches):]:
        ax.set_visible(False)
    for ax in axes[(rows - 1) * columns:]:
        ax.set_xlabel("PCoA 1 (%.1f%%)" % (100 * share[0]), color=INK_SECONDARY)
    for ax in axes[::columns]:
        ax.set_ylabel("PCoA 2 (%.1f%%)" % (100 * share[1]), color=INK_SECONDARY)
    fig.suptitle("Deer tick microbiome by sequencing batch", x=0.01, y=0.995,
                 ha="left", fontsize=12, fontweight="bold", color=INK)
    fig.text(0.01, 0.955, STATS, ha="left", va="top", fontsize=8.5,
             color=INK_SECONDARY)
    fig.text(0.01, 0.005, "Blue: deer in the named batch. Grey: all other deer. "
             "Cross: batch centroid.", ha="left", va="bottom", fontsize=8,
             color=INK_SECONDARY)
    fig.tight_layout(rect=(0, 0.03, 1, 0.945))
    for extension in ("png", "pdf"):
        fig.savefig("%s.%s" % (OUT, extension), dpi=300)
    print("wrote %s.png and .pdf (%d deer, %d batches; PCoA 1 %.1f%%, "
          "PCoA 2 %.1f%%)" % (OUT, len(distances), len(batches),
                              100 * share[0], 100 * share[1]))


if __name__ == "__main__":
    main()
