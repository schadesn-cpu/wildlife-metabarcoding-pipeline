#!/usr/bin/env python3
"""archive_retired_figures.py

I move the figures the restructure retired into results/_archive, with a
manifest saying what each one was and why it went. I move, never delete, and I
print what I would do unless you pass --write.

WHERE
-----
    results/_archive/2026-09-24_restructure

That follows the pattern already in results/_archive (2026-09-18_superseded,
diversity_july2026). I create that one directory and nothing else.

WHAT I LEAVE ALONE
------------------
fig3_qpcr_concordance, because Section 3.5 is still undecided and archiving a
figure we may need back is worse than leaving it where it is.

Every analysis table. The co-occurrence and diversity TSVs are the evidence
behind decisions the manuscript now takes in one sentence, and they cost
nothing to keep.

USAGE
-----
    cd /home/users/sne24/meedlab/tick_project
    python3 scripts/archive_retired_figures.py
    python3 scripts/archive_retired_figures.py --write
"""
import glob
import os
import shutil
import sys
from datetime import date

SOURCES = ("results/figures", "results/figures_manuscript")
DEST = "results/_archive/2026-09-24_restructure"
KEEP = ("fig3_qpcr_concordance",)

RETIRE = {
    "16S_baseline_dominant_taxon_wong":
        "exploratory baseline, per library, superseded by the animal-level "
        "core microbiome",
    "fig1a_host_shannon_sig":
        "diversity section cut in the restructure",
    "fig1b_deer_batch_pcoa":
        "batch effect section cut; the finding survives as two Methods "
        "sentences on terminal-base variants",
    "fig1b_deer_by_batch_stats":
        "same as fig1b_deer_batch_pcoa",
    "fig2_single_tick_prevalence_coinfection":
        "co-infection framing dropped from the paper",
    "fig_qpcr_prevalence_coinfection":
        "replaced by qPCR detection by host, single ticks and pools",
    "fig_tick_species_by_host":
        "replaced by Table 1, sample set by host",
    "fig_microbiome_A_full":
        "built 2026-08-11 on libraries rather than animals, before the "
        "host-lineage filter and the classified-reads denominator",
    "fig_microbiome_B_noendo":
        "same as fig_microbiome_A_full",
    "fig_microbiome_by_host":
        "same as fig_microbiome_A_full",
    "fig_bloodmeal_matrix":
        "blood meal is now a paragraph; the matrix counted any read as a "
        "meal, so trace contamination appeared as a second host",
}


def collect():
    """I return [(path, stem, reason)] for every file I would move."""
    found = []
    for folder in SOURCES:
        if not os.path.isdir(folder):
            continue
        for path in sorted(glob.glob(os.path.join(folder, "*"))):
            if not os.path.isfile(path):
                continue
            stem = os.path.splitext(os.path.basename(path))[0]
            if stem in KEEP:
                continue
            if stem in RETIRE:
                found.append((path, stem, RETIRE[stem]))
    return found


def main():
    write = "--write" in sys.argv[1:]
    if not os.path.isdir("results"):
        sys.exit("no results/ directory here. Run me from the project root.")
    if not os.path.isdir(os.path.dirname(DEST)):
        sys.exit("%s does not exist. I will not create it."
                 % os.path.dirname(DEST))
    found = collect()
    if not found:
        sys.exit("nothing to archive. Either it has been done already or the "
                 "file names have changed.")

    total = sum(os.path.getsize(p) for p, _, _ in found)
    print("moving %d files (%.1f MB) to %s" % (len(found), total / 1e6, DEST))
    print()
    print("  %-52s %8s" % ("file", "size"))
    for path, _, _ in found:
        print("  %-52s %7.1fM" % (path[:52], os.path.getsize(path) / 1e6))
    print()
    kept = []
    for folder in SOURCES:
        for stem in KEEP:
            kept += sorted(glob.glob(os.path.join(folder, stem + ".*")))
    print("leaving in place: %s" % (", ".join(kept) if kept else "nothing"))
    print()

    if not write:
        print("Nothing moved. Rerun with --write to do it.")
        return

    clashes = [p for p, _, _ in found
               if os.path.exists(os.path.join(DEST, os.path.basename(p)))]
    if clashes:
        sys.exit("these already exist in %s: %s" % (DEST, ", ".join(clashes)))
    os.makedirs(DEST, exist_ok=True)

    lines = ["RETIRED FIGURES, %s" % date.today().isoformat(),
             "=" * 60,
             "",
             "Retired when the manuscript was restructured onto the agreed",
             "outline: pathogen detection by host, then the core microbiome of",
             "each tick and host pair. Nothing here was wrong when it was",
             "built. Each one either belongs to a section that was cut or was",
             "built on a superseded cohort.",
             "",
             "The surviving figures are the collection map, qPCR detection by",
             "host, and the core microbiome, all built after 2026-09-24.",
             ""]
    for path, stem, reason in found:
        shutil.move(path, os.path.join(DEST, os.path.basename(path)))
        lines.append("%-52s %s" % (os.path.basename(path), reason))
    lines.append("")
    lines.append("fig3_qpcr_concordance was left in place: Section 3.5 is "
                 "still undecided.")
    with open(os.path.join(DEST, "MANIFEST.txt"), "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("moved %d files and wrote %s/MANIFEST.txt" % (len(found), DEST))


if __name__ == "__main__":
    main()
