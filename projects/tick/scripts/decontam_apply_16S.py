#!/usr/bin/env python3
"""
Categorize decontam-flagged 16S ASVs into REMOVE / KEEP / BORDERLINE and
(only with --commit) write a filtered feature table. DRY RUN by default.

Rules (loon-informed, taxon-aware):
  KEEP (never remove):
    - host DNA: taxon contains Mammalia / Vertebrata / Craniata / Chloroplast /
      Mitochondria  (handled by the host-filter step, not decontam)
    - pathogens/endosymbionts of interest: Rickettsia, Anaplasma, Borrelia,
      Borreliella, Bartonella, Francisella, Babesia (safety net; none were
      flagged, but we guarantee it)
  REMOVE (clear reagent/kitome contaminant):
    - a known kitome genus AND blank-enriched enough:
      reads_neg >= 0.5 * reads_samp  (present in blanks at >= half the sample reads)
  BORDERLINE (report, do NOT auto-remove; decide by hand):
    - flagged but sample-enriched (reads_samp > 2 * reads_neg) and not clearly host
      -> e.g. Escherichia here: abundant in samples, also in blanks. Kept by default,
         listed so it gets an explicit call + a methods sentence.

Run:
  python scripts/decontam_apply_16S.py            # dry run: prints the 3 buckets
  python scripts/decontam_apply_16S.py --commit    # writes filtered table + removal list
"""
import sys, re
import pandas as pd

INSP  = "results/decontam/16S/inspection.tsv"
TABLE = "/tmp/16s_feature_table.tsv"
OUT_REMOVE = "results/decontam/16S/remove_asvs.txt"
OUT_TABLE  = "results/decontam/16S/table_decontam_16S.tsv"

KEEP_TAXA = re.compile(
    r"Mammalia|Vertebrata|Craniata|Chloroplast|Mitochondria|"
    r"Rickettsia|Anaplasma|Borreliella|Borrelia|Bartonella|Francisella|Babesia",
    re.I)

# well-known reagent / skin (kitome) genera
KITOME = re.compile(
    r"Cutibacterium|Propionibacter|Corynebacterium|Micrococcus|Staphylococcus|"
    r"Streptococc|Lactococcus|Lactobacillus|Pseudomonas|Acinetobacter|Enhydrobacter|"
    r"Sphingomonas|Sphingobium|Bacillus|Comamonadaceae|Massilia|Methylobacterium|"
    r"Bradyrhizobium|Ralstonia|Paracoccus|Halomonas|Mycobacterium|Chryseobacterium|"
    r"Flavobacterium|Amaricoccus|Tepidimonas|Pseudorhodoferax|Finegoldia|Delftia|"
    r"Rhodococcus|Brevundimonas|Stenotrophomonas|Oxalobacteraceae|Rhodobacteraceae|"
    r"Sphingomonadaceae|Neisseriaceae|Oscillospiraceae|Halomonadaceae|Spirosomaceae",
    re.I)

BLANK_ENRICH = 0.5   # reads_neg >= this * reads_samp  -> contaminant side
SAMPLE_ENRICH = 2.0  # reads_samp > this * reads_neg   -> sample side (keep/borderline)


# genera that could plausibly be REAL tick-associated / gut taxa, not just kitome.
# These stay BORDERLINE (hand decision) even when flagged, rather than auto-remove.
PLAUSIBLE_REAL = re.compile(
    r"Escherichia|Shigella|Lactococcus|Lactobacillus|Enterococcus|"
    r"Enterobacter|Klebsiella|Serratia|Rickettsiella|Spiroplasma|Wolbachia|"
    r"Coxiella|Midichloria|Arsenophonus",
    re.I)


def categorize(row):
    tax = str(row.get("Taxon", ""))
    rn, rs = row["reads_neg"], row["reads_samp"]
    if KEEP_TAXA.search(tax):
        return "KEEP", "host/pathogen taxon"
    # a flagged, named kitome genus is a contaminant regardless of sample/blank
    # ratio: kitome organisms amplify in real samples too, so sample-enrichment
    # does NOT exonerate them. Only spare plausibly-real genera.
    if KITOME.search(tax) and not PLAUSIBLE_REAL.search(tax):
        return "REMOVE", "flagged kitome genus"
    if PLAUSIBLE_REAL.search(tax):
        return "BORDERLINE", "plausibly real (e.g. Escherichia); keep, decide by hand"
    if rn >= BLANK_ENRICH * max(rs, 1):
        return "BORDERLINE", "blank-enriched but untyped; likely contaminant, review"
    # untypable (bare d__Bacteria / Unassigned), sample-enriched: keep as borderline
    return "BORDERLINE", "untyped; keep, review"


def main():
    commit = "--commit" in sys.argv
    insp = pd.read_csv(INSP, sep="\t")
    cats = insp.apply(categorize, axis=1, result_type="expand")
    insp["decision"], insp["reason"] = cats[0], cats[1]

    for d in ("REMOVE", "BORDERLINE", "KEEP"):
        sub = insp[insp["decision"] == d]
        print("\n=== %s  (%d ASVs) ===" % (d, len(sub)))
        for _, r in sub.iterrows():
            tax = str(r.get("Taxon", ""))
            short = tax.split(";")[-2] if ";" in tax else tax
            short = short[:40]
            print("  %s  samp=%-6d neg=%-5d  %-34s  [%s]" %
                  (str(r["ASV"])[:8], int(r["reads_samp"]), int(r["reads_neg"]),
                   short, r["reason"]))

    remove = insp.loc[insp["decision"] == "REMOVE", "ASV"].tolist()
    print("\n----------------------------------------")
    print("REMOVE: %d   BORDERLINE(keep): %d   KEEP: %d   of %d flagged"
          % ((insp.decision=="REMOVE").sum(), (insp.decision=="BORDERLINE").sum(),
             (insp.decision=="KEEP").sum(), len(insp)))

    if not commit:
        print("\nDRY RUN. Nothing written. Re-run with --commit to:")
        print("  - write %s (the %d REMOVE ASVs)" % (OUT_REMOVE, len(remove)))
        print("  - write %s (feature table minus those ASVs)" % OUT_TABLE)
        print("BORDERLINE ASVs are KEPT (not removed); review them by hand.")
        return

    # commit: write removal list and filtered table
    with open(OUT_REMOVE, "w") as fh:
        fh.write("\n".join(remove) + "\n")
    tbl = pd.read_csv(TABLE, sep="\t", skiprows=1)
    tbl = tbl.rename(columns={tbl.columns[0]: "ASV"}).set_index("ASV")
    before = tbl.shape[0]
    tbl = tbl.drop(index=[a for a in remove if a in tbl.index])
    tbl.to_csv(OUT_TABLE, sep="\t")
    print("\nremoved %d ASVs (%d -> %d). wrote:\n  %s\n  %s"
          % (before - tbl.shape[0], before, tbl.shape[0], OUT_REMOVE, OUT_TABLE))
    print("Note: filtered table is TSV. Re-import to .qza for QIIME diversity if needed.")


if __name__ == "__main__":
    main()
