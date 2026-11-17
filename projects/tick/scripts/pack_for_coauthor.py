#!/usr/bin/env python3
"""pack_for_coauthor.py

I pack the tables and figures a coauthor needs to check the tick analysis into
ONE file, results/for_lawrence_<date>.tar.gz, with a MANIFEST.txt inside that
says what each file is. I create no directory and change no file.

WHAT GOES IN
------------
ASV tables for all four markers, the taxonomy each marker was classified
against, the BLAST verification tables, the per-sample tables the paper is
built from, and the figures the manuscript uses (results/figures_manuscript/,
with the matching PDF from results/figures/ where one exists).

WHAT NEVER GOES IN
------------------
Anything under metadata/. The metadata carries collection towns, coordinates,
collector names and hunter tag numbers. If the coauthor needs host or tick
counts, send them separately and deliberately, not by accident inside this
archive. I also refuse .bak files.

If a listed file is missing I stop and name it, rather than ship a partial
package that looks complete.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/pack_for_coauthor.py            # show what would go in
    python3 scripts/pack_for_coauthor.py --write    # write the archive
"""
import glob
import io
import os
import sys
import tarfile
import time

FILES = [
    ("results/tables/asv_table_16S.tsv",
     "16S ASV table: reads per ASV per library."),
    ("results/tables/asv_table_18S.tsv",
     "18S ASV table (tick species, piroplasms)."),
    ("results/tables/asv_table_cytb.tsv",
     "cytb ASV table (bloodmeal host)."),
    ("results/tables/its_asv_table_FILTERED.tsv",
     "ITS ASV table after filtering."),
    ("results/16S/all/taxonomy_verified/taxonomy_16S_droplist_applied.tsv",
     "16S taxonomy after BLAST verification and the drop list."),
    ("results/18S/all/taxonomy_verified/taxonomy_18S_blast_verified.tsv",
     "18S taxonomy after BLAST verification."),
    ("results/16S/all/blast_verify/pathogen_asvs_blast.tsv",
     "BLAST hits for every 16S ASV assigned to a pathogen genus."),
    ("results/18S/all/blast_verify/piroplasm_asvs_blast.tsv",
     "BLAST hits for 18S piroplasm ASVs."),
    ("results/18S/all/taxonomy_verified/tick_genus_per_sample_corrected.tsv",
     "Tick genus per sample from 18S."),
    ("results/analysis/tick_species_by_host.tsv",
     "Tick species by host, one row per specimen (476)."),
    ("results/analysis/analysis_table_16S.tsv",
     "Per-library pathogen reads and calls (602 libraries, 413 animals). "
     "Calls use >=1 read; libraries need >=1,000 reads to enter."),
]
# Only the figures the manuscript uses. results/figures/ also holds exploratory
# plots from June to August that a reader would mistake for current results.
FIGURE_GLOB = "results/figures_manuscript/*.png"
# Add a figure here only once it is checked against the manuscript numbers.
EXTRA_FIGURES = []
FORBIDDEN = ("metadata/",)


def human(size):
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024:
            return "%.0f %s" % (size, unit)
        size /= 1024.0
    return "%.1f TB" % size


def main():
    if not os.path.isdir("results"):
        sys.exit("no results/ here. Run from "
                 "/home/users/sne24/meedlab/tick_project.")
    write = "--write" in sys.argv[1:]
    entries = list(FILES)
    figures = sorted(glob.glob(FIGURE_GLOB))
    if not figures:
        sys.exit("no figures matched %s." % FIGURE_GLOB)
    for path in figures:
        entries.append((path, "Manuscript figure (PNG)."))
        stem = os.path.splitext(os.path.basename(path))[0]
        vector = os.path.join("results/figures", stem + ".pdf")
        if os.path.isfile(vector):
            entries.append((vector, "Same figure as PDF."))
    for path in EXTRA_FIGURES:
        entries.append((path, "Manuscript figure, added by hand."))

    missing = [path for path, _ in entries if not os.path.isfile(path)]
    if missing:
        sys.exit("missing, so I am not packing anything:\n  "
                 + "\n  ".join(missing))
    for path, _ in entries:
        if path.startswith(FORBIDDEN) or ".bak" in path:
            sys.exit("refusing %s: metadata and .bak files never go in." % path)

    total = 0
    lines = ["Tick project files for coauthor review, packed %s."
             % time.strftime("%Y-%m-%d"),
             "Sample metadata is NOT included: it carries collection "
             "locations and hunter records.",
             "Library names map to animals through the project lookup; "
             "several libraries can come from one animal.",
             ""]
    for path, note in entries:
        size = os.path.getsize(path)
        total += size
        stamp = time.strftime("%Y-%m-%d",
                              time.localtime(os.path.getmtime(path)))
        print("  %s  %8s  %s" % (stamp, human(size), path))
        lines.append("%s  (%s, %s)\n    %s" % (path, stamp, human(size), note))
    print("  %d files, %s before compression" % (len(entries), human(total)))

    target = "results/for_lawrence_%s.tar.gz" % time.strftime("%Y%m%d")
    if not write:
        print()
        print("Nothing written. Rerun with --write to make %s" % target)
        return
    if os.path.exists(target):
        sys.exit("%s already exists. Move it aside first." % target)
    manifest = ("\n".join(lines) + "\n").encode("utf-8")
    with tarfile.open(target, "w:gz") as archive:
        info = tarfile.TarInfo("MANIFEST.txt")
        info.size = len(manifest)
        info.mtime = int(time.time())
        archive.addfile(info, io.BytesIO(manifest))
        for path, _ in entries:
            archive.add(path)
    print()
    print("wrote %s (%s)" % (target, human(os.path.getsize(target))))


if __name__ == "__main__":
    main()
