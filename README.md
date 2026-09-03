# README Corrections — wildlife-metabarcoding-pipeline (loon deposit)
Paste these over the stale sections before making the repo public / minting the Zenodo DOI.
All numbers reconciled to the final manuscript + GROUND_TRUTH.

============================================================
## FIX 1 — Header date
------------------------------------------------------------
REPLACE:  **Last updated:** 2026-03-27
WITH:     **Last updated:** 2026-09-03

============================================================
## FIX 2 — Overview paragraph (18S is INCLUDED)
------------------------------------------------------------
The current Overview says 18S is under a "see note below" caveat, and the note
below says 18S was excluded. That is stale — 18S is a core marker in the final
paper (parasite/eukaryote detection, Fig 5). Replace the Overview with:

Amplicon sequencing study of the common loon (*Gavia immer*) gut microbiome and
dietary composition. Four markers were sequenced: 16S rRNA gene (bacteria),
MiFish 12S (fish diet), cytochrome b (vertebrate diet), and 18S rRNA V9
(eukaryotes and parasites). Lung tissue was screened for herpesviruses and
adenoviruses by consensus PCR, and swabs for influenza A by RT-PCR.

Dietary markers (MiFish, cytb) are reported as presence/absence and detection
frequency, not relative abundance. Relative read abundance is not a defensible
proxy for dietary biomass because PCR efficiency differs across prey taxa,
mitochondrial copy number varies by tissue and species, and digestion state
affects DNA yield independently of consumption (Deagle et al. 2019, Mol Ecol).

============================================================
## FIX 3 — Delete / replace the "18S status: Excluded" note
------------------------------------------------------------
DELETE the entire "18S status: Excluded from the main analysis... failed
amplification..." paragraph. It is from an early project phase and contradicts
the final paper. Replace with:

**18S marker:** Included as the eukaryote/parasite marker. 18S rRNA V9 detections
(>=10 reads per sample, or independently identified at gross necropsy) were used
for parasite detection and necropsy concordance (Fig 5). 18S cohort n = 35.
Host reads were identified and removed by BLAST-based sequence assignment.

============================================================
## FIX 4 — Sample Groups table (numbers reconciled to final paper)
------------------------------------------------------------
The current table shows Diseased 13 / Trauma 13 / Marine 3-5, which is stale.
The final manuscript uses Disease 13 / Trauma 12 for the 16S DvT contrast, and
the cohorts differ by analysis. Replace the Sample Groups section with:

Cohorts differ by analysis (rarefaction and season-exclusion rules apply):

| Analysis                        | n   | Composition                                   |
|---------------------------------|-----|-----------------------------------------------|
| 16S season PERMANOVA (all)      | 25  | Breeding / FW-nonbreeding / Saltwater         |
| 16S season PERMANOVA (adults)   | 15  | Breeding 10 / Saltwater 5 (Br-vs-Sw contrast) |
| 16S Disease-vs-Trauma (DvT)     | 25  | Disease 13 / Trauma 12                        |
| 16S taxonomy barplot (adults)   | 17  | unrarefied adult cohort (Fig 3)               |
| 18S parasite detection          | 35  | TV240106 RETAINED (non-season analysis)       |
| Diet (MiFish/cytb pooled)       | 20  | season diet transition (Fig 4)                |
| Viral screen (lung / swab)      | 40  | full cohort                                   |

Key exclusion rule: **TV240106 is excluded from all SEASON analyses**
(EcoSeason = Unknown; the bird froze in a lake, out of season) but is RETAINED
in non-season analyses (parasites/18S, diet, viral). This is why season cohorts
are one lower than the raw marker counts.

Canonical season stats: EcoSeason PERMANOVA n=25, pseudo-F=1.89, p=0.005;
adult Br-vs-Sw n=15, F=2.71, p=0.006; PERMDISP p=0.205 / 0.108 (ns).
DvT PERMANOVA n=25, F=1.37, p=0.131 (ns). Diet KW p=0.013.

============================================================
## FIX 5 — Remove the stale "verify 18S n / expected n=30" note
------------------------------------------------------------
DELETE: "18S sample counts (rarefaction depth 1,000): n=13 Diseased + 12 Trauma
+ 5 Marine = 30 total. Note: Verify final 18S n..."
The final 18S parasite cohort is n=35 (see Fix 3/4). The n=30 was an interim
DvT-subset count, not the parasite-detection cohort.

============================================================
## OPTIONAL — add a short "Correct taxonomy table" caveat
------------------------------------------------------------
If any taxonomy_refined output is present in the deposit, add:

**Note on taxonomy tables:** Figure 3 and reported genus relative abundances use
the plain taxonomy table (results/16S/all/taxonomy/taxonomy_relabund_L6_16S.tsv).
A taxonomy_refined/ variant exists from an exploratory reassignment and is
superseded; do not use it to reproduce the figures.
