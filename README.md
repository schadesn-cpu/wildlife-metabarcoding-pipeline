# Common Loon Gut Microbiome & Dietary Metabarcoding — Project README

**Last updated:** 2026-09-03

---

## Overview

Amplicon sequencing study of the common loon (*Gavia immer*) gut microbiome and dietary composition. Four markers were sequenced: **16S rRNA gene** (bacteria), **MiFish 12S** (fish diet), **cytochrome b** (vertebrate diet), and **18S rRNA V9** (eukaryotes and parasites). Lung tissue was screened for herpesviruses and adenoviruses by consensus PCR, and swabs for influenza A by RT-PCR.

This pipeline is a reproducible, marker-agnostic QIIME2-based workflow for wildlife amplicon sequencing, including multi-marker support, diversity statistics (PERMANOVA, PERMDISP), automated summary reporting, and scalable HPC deployment.

**Dietary markers (MiFish, cytb) are reported as presence/absence and detection frequency, not relative abundance.** Relative read abundance is not a defensible proxy for dietary biomass because PCR efficiency differs across prey taxa, mitochondrial copy number varies by tissue and species, and digestion state affects DNA yield independently of consumption (Deagle et al. 2019, *Molecular Ecology*). cytb additionally rarefies to a shallow depth at which relative-abundance estimates have variance too large to interpret.

**18S marker:** Included as the eukaryote/parasite marker. 18S rRNA V9 detections (≥10 reads per sample, or independently identified at gross necropsy) were used for parasite detection and necropsy concordance (Fig 5). Host reads were identified and removed by BLAST-based sequence assignment.

---

## Environment Setup

Two conda environments — one for QIIME2 pipeline steps, one for Python analysis and plotting scripts.

**QIIME2 environment (scripts 00–05, run_all_figures.sh):**
```
conda env create -f environment_qiime2.yml
conda activate metabarcoding-qiime2
```

**Analysis environment (scripts 06–14 and utilities):**
```
conda env create -f environment_analysis.yml
conda activate metabarcoding-analysis
```

---

## Sample Groups & Cohorts

Cohorts differ by analysis, because rarefaction and season-exclusion rules apply differently to each.

| Analysis                          | n  | Composition                                         |
|-----------------------------------|----|-----------------------------------------------------|
| 16S season PERMANOVA (all)        | 25 | Breeding / Freshwater non-breeding / Saltwater      |
| 16S season PERMANOVA (adults)     | 15 | Breeding 10 / Saltwater 5 (breeding-vs-wintering)   |
| 16S Disease-vs-Trauma (DvT)       | 25 | Disease 13 / Trauma 12                              |
| 16S taxonomy barplot (adults)     | 17 | unrarefied adult cohort (Fig 3)                     |
| 18S parasite detection            | 35 | TV240106 retained (non-season analysis)             |
| Diet (MiFish/cytb pooled)         | 20 | seasonal diet transition (Fig 4)                    |
| Viral screen (lung / swab)        | 40 | full cohort                                         |

**Key exclusion rule — TV240106:** Excluded from all **season** analyses (its EcoSeason is "Unknown"; the bird froze in a lake, out of season) but **retained** in non-season analyses (parasites/18S, diet, viral screen). This is why season cohorts are one lower than the raw marker counts.

**Canonical statistics:**
- Season PERMANOVA (16S): n = 25, pseudo-F = 1.89, p = 0.005; adult breeding-vs-wintering n = 15, F = 2.71, p = 0.006; PERMDISP p = 0.205 and 0.108 (ns).
- Disease-vs-Trauma PERMANOVA (16S): n = 25, F = 1.37, p = 0.131 (ns).
- Diet (marine-ratio Kruskal–Wallis): p = 0.013.

---

## Pipeline Structure

Scripts are numbered in execution order.

**QIIME2 pipeline (00–05):**
- `00_build_classifiers.py` — build marker-specific taxonomy classifiers
- `00_merge_run_dirs.py` — merge multi-run sequencing directories
- `01_make_manifests.py` — build QIIME2 import manifests
- `02_make_qiime_metadata.py` — assemble per-marker metadata
- `02b_add_season_to_metadata.py` — add ecological season assignments
- `03_run_full_metabarcoding_pipeline.py` — import, denoise (DADA2), taxonomy
- `04_rarefaction.py` — rarefaction curves and depth selection
- `05_run_diversity_stats.py`, `05b_run_cod_diversity.py`, `05c_parse_beta_stats.py` — diversity metrics and PERMANOVA/PERMDISP

**Analysis & figures (06–14):**
- `06_plot_diversity.py`, `06b_combine_multimarker_alpha.py`, `07_visualize_diversity.py` — diversity visualization (Fig 2 PCoA)
- `08_taxonomy_table.py`, `08b_presence_absence.py`, `08c_blast_verify.py` — taxonomy tables and BLAST-based host/taxon verification
- `09_plot_taxonomy.py` — taxonomy barplots (Fig 3, Fig S1)
- `09b_clean_diet_table.py`, `09b_plot_mifish_season_ecology.py` — diet tables and plots
- `10_plot_viral.py`, `10b_plot_adeno_tree.py`, `10b_plot_herpes_cutadapt.py`, `10c_plot_adeno_cutadapt.py` — viral screening and adenovirus phylogeny (Fig 6)
- `10b_annotate_diet_ecology.py`, `11_plot_mifish_season_ecology.py`, `12_plot_habitat_season.py` — dietary ecology (Fig 4)
- `14_viral_stats.py` — viral detection summaries
- `run_all_figures.py` — config-driven figure regeneration
- `utils/` — shared helpers; `scan_files.py` — inventory utility

---

## Note on taxonomy tables

Figure 3 and all reported genus relative abundances use the **plain** taxonomy table (`results/16S/all/taxonomy/taxonomy_relabund_L6_16S.tsv`). A `taxonomy_refined/` variant exists from an exploratory reassignment and is **superseded** — do not use it to reproduce the figures.

---

## Data Availability

- Raw amplicon sequence reads: NCBI Sequence Read Archive (BioProject accession pending deposit).
- Adenovirus DNA polymerase (DPOL) sequence: GenBank (accession pending deposit).
- Sample metadata are formatted to the MIxS/MIMARKS standard.
- This repository is archived on Zenodo (DOI on release).

---

## Citation

If you use this pipeline, please cite the associated manuscript (in preparation) and this repository.
