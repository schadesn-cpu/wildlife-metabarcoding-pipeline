#!/usr/bin/env python3
"""
Archive confirmed-superseded files into a mirrored _archive/ tree, reversibly.

SAFE BY DESIGN:
  * Dry-run by default. Nothing moves unless you pass --commit.
  * Reference guard: any file referenced by an active script in scripts/
    is SKIPPED (never moved), so the pipeline can't be orphaned.
  * Every move is logged to _archive/ARCHIVE_MANIFEST.tsv (orig -> archived),
    so it is fully reversible.

SCOPE (deliberately conservative):
  CONFIRMED  -> moved (after reference check). Only files we have positively
                established as superseded this session.
  REVIEW     -> printed for manual decision, NEVER moved automatically.

Run from project root:
  python scripts/archive_superseded.py            # dry run, shows plan
  python scripts/archive_superseded.py --commit    # actually move
"""
import os, sys, shutil, glob, datetime

ROOT = os.path.expanduser("~/meedlab/tick_project")
ARCH = os.path.join(ROOT, "_archive")
MANIFEST = os.path.join(ARCH, "ARCHIVE_MANIFEST.tsv")

# --- CONFIRMED superseded (this session). Paths relative to ROOT. ---
CONFIRMED = [
    "qiime2/18S/all/dada2/table_paired.qza",
    "qiime2/18S/all/dada2/rep-seqs_paired.qza",
    "qiime2/18S/all/dada2/denoising-stats_paired.qza",
    "qiime2/18S/all/taxonomy/taxonomy_OLD_BAD.qza",
    # superseded 2026-07-02 by results/figures/single_tick_prevalence_coinfection.*
    # (verified, italic binomials, endosymbionts separated, per-genus co-infection)
    "results/figures/fig2_pathogen_prevalence.png",
    "results/figures/fig2_pathogen_prevalence.pdf",
    "results/figures/fig3_coinfection_rates.png",
    "results/figures/fig3_coinfection_rates.pdf",
]

# --- REVIEW: likely superseded but NOT auto-moved. Printed for a human call. ---
REVIEW = [
    # old 16S taxonomy (515-806 era) vs region-matched; confirm which before moving
    "results/16S/all/taxonomy/taxonomy_relabund_L6_16S.tsv",
    "results/16S/all/taxonomy/taxonomy_counts_L6_16S.tsv",
    "results/16S/all/taxonomy/relabund_taxonomy.tsv",  # non-regionmatched
    # fig1/fig5 stale May figures: archive only AFTER regenerating replacements
    "results/figures/fig1_host_composition.png",
    "results/figures/fig5_cooccurrence_heatmap.png",
    # cytb entrez_test iterations: trace which fed cytb_host_summary.tsv first
    "results/cytb/all/blast_verify/entrez_test",
    "results/cytb/all/blast_verify/entrez_test2",
    "results/cytb/all/blast_verify/entrez_test3",
]

def basename(p): return os.path.basename(p.rstrip("/"))

def referenced_by_scripts(fname):
    """Return list of active scripts referencing this basename.
    Skips this archive script itself (its CONFIRMED list names the files,
    which would otherwise be a false self-reference)."""
    hits = []
    self_name = os.path.basename(__file__)
    for pat in ("scripts/*.py", "scripts/*.sh"):
        for f in glob.glob(os.path.join(ROOT, pat)):
            if os.path.basename(f) == self_name:
                continue
            try:
                with open(f, "r", errors="ignore") as fh:
                    if fname in fh.read():
                        hits.append(os.path.relpath(f, ROOT))
            except Exception:
                pass
    return hits

def main():
    commit = "--commit" in sys.argv
    force = "--force" in sys.argv
    print("MODE:", "COMMIT (files will move)" if commit else "DRY RUN (no changes)",
          "+ FORCE" if force else "")
    print("archive root:", ARCH)

    moves, skipped = [], []
    for rel in CONFIRMED:
        src = os.path.join(ROOT, rel)
        if not os.path.exists(src):
            print("  [missing] %s (already moved or never existed)" % rel); continue
        refs = referenced_by_scripts(basename(rel))
        if refs and not force:
            skipped.append((rel, refs))
            print("  [SKIP: referenced] %s  <- %s" % (rel, ", ".join(refs)))
            continue
        if refs and force:
            print("  [FORCE move: ref overridden] %s  <- %s" % (rel, ", ".join(refs)))
            moves.append((rel, refs))
        else:
            print("  [move] %s" % rel)
            moves.append((rel, []))

    print("\n--- REVIEW (never auto-moved; decide by hand) ---")
    for rel in REVIEW:
        src = os.path.join(ROOT, rel)
        exists = "exists" if os.path.exists(src) else "absent"
        print("  [review:%s] %s" % (exists, rel))

    if not commit:
        print("\nDry run only. Re-run with --commit to move the [move] items above.")
        print("%d to move, %d skipped for references." % (len(moves), len(skipped)))
        return

    os.makedirs(ARCH, exist_ok=True)
    stamp = datetime.datetime.now().isoformat(timespec="seconds")
    newlog = not os.path.exists(MANIFEST)
    with open(MANIFEST, "a") as log:
        if newlog:
            log.write("timestamp\toriginal_path\tarchived_path\toverridden_refs\n")
        for rel, refs in moves:
            src = os.path.join(ROOT, rel)
            dst = os.path.join(ARCH, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.move(src, dst)
            refnote = ";".join(refs) if refs else ""
            log.write("%s\t%s\t%s\t%s\n" % (stamp, rel, os.path.relpath(dst, ROOT), refnote))
            print("  moved %s -> %s%s" % (rel, os.path.relpath(dst, ROOT),
                  ("  (ref overridden: %s)" % refnote) if refnote else ""))
    print("\ndone. %d moved. manifest: %s" % (len(moves), os.path.relpath(MANIFEST, ROOT)))
    print("to reverse: read ARCHIVE_MANIFEST.tsv and shutil.move back.")

if __name__ == "__main__":
    main()
