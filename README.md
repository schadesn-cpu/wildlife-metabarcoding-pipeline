# Common Loon Gut Microbiome & Dietary Metabarcoding — Project README

**Last updated:** 2026-03-19

---

## Overview

Amplicon sequencing study of Common Loon (*Gavia immer*) gut microbiome and dietary
composition. Four markers: **16S rRNA** (bacteria), **MiFish 12S** (fish diet),
**cytochrome b** (vertebrate diet), **18S rRNA** (eukaryotes/parasites).
Viral detection via pan-herpesvirus TGF-IYG primers on lung tissue.

---

## Environment Setup

Two conda environments are used — one for QIIME2 pipeline steps, one for
Python analysis and plotting scripts (no QIIME2 required).

### QIIME2 environment (scripts 00–05, run_all_figures.sh)
```bash
# Follow official QIIME2 installation for your platform:
# https://docs.qiime2.org/2024.5/install/
conda env create -f environment_qiime2.yml
conda activate metabarcoding-qiime2
```

### Analysis environment (scripts 06–10 and utilities)
```bash
conda env create -f environment_analysis.yml
conda activate metabarcoding-analysis
```

---

## Sample Groups

| Group | n (rarefied) | Notes |
|---|---|---|
| Diseased | 13 (16S/MiFish/cytb) | Lead toxicosis or infectious disease COD |
| Trauma | 13 (16S/MiFish/cytb) | Trauma COD |
| Marine | 3–5 | Excluded from DvT diversity analyses |

Three samples dropped at rarefaction (8,000 reads): TV240046 (Diseased), TV230067 + TV240036 (Trauma).
Marine samples (TV220031, TV230063, TV240057) excluded from 16S/MiFish/cytb DvT diversity.

**18S sample counts (rarefaction depth 1,000):** n=13 Diseased + 12 Trauma + 5 Marine = 30 total.
⚠ Note: Verify final 18S n — two samples may lack Collection_source values; expected n=30.

---

## Directory Structure

```
project_root/
├── metadata/
│   └── qiime/
│       ├── metadata_{marker}.tsv              ← QIIME2 metadata (DvT analyses)
│       ├── metadata_{marker}_cod.tsv          ← with COD_broad + Collection_source
│       └── metadata_cytb_cod_filtered.tsv     ← Lead + Parasitic_Infectious + Trauma only
├── qiime2/                                    ← QIIME2 artifacts (.qza/.qzv)
│   ├── 16S/rarefied_8000/DvT/diversity/core_metrics_depth8000/
│   ├── MiFish/rarefied_17000/DvT/diversity/core-metrics-17000/
│   ├── cytb/all/diversity/core-metrics-200/
│   ├── 18S/all/diversity/core-metrics-1000/
│   └── cytb/COD/                              ← filtered distance matrices for cytb COD stats
├── results/                                   ← all outputs (flat, clean structure)
│   ├── 16S/
│   │   ├── DvT/{diversity,figures,taxonomy}/
│   │   ├── season/{diversity,figures}/
│   │   └── COD/{diversity,figures}/           ← COD_broad + Collection_source stats
│   ├── MiFish/  (same layout)
│   ├── cytb/    (same layout)
│   ├── 18S/all/{diversity,figures,taxonomy}/
│   ├── herpesvirus/figures/                   ← herpesvirus detection figures
│   ├── adenovirus/figures/                    ← aviadenovirus detection figures
│   └── _archive_old_structure/                ← old files, safe to delete
├── scripts/
│   ├── 00_build_classifiers.py
│   ├── 01_make_manifests.py
│   ├── 02_make_qiime_metadata.py
│   ├── 03_run_full_metabarcoding_pipeline.py
│   ├── 04_rarefaction.py
│   ├── 05_run_diversity_stats.py              ← PERMANOVA/PERMDISP/alpha-sig QZVs
│   ├── 05b_run_cod_diversity.py               ← COD_broad + Collection_source stats + figures
│   ├── 06_plot_diversity.py                   ← PCoA and alpha figures (patched)
│   ├── 08_taxonomy_table.py
│   ├── 09_plot_taxonomy.py                    ← stacked barplots (patched)
│   ├── 10_plot_viral.py                       ← herpesvirus + adenovirus figures
│   ├── primer_advisor.py                      ← primer detection + DADA2 param suggestion
│   ├── parse_beta_stats.py                    ← extract PERMANOVA/PERMDISP from QZVs
│   ├── add_season_to_metadata.py
│   ├── scan_files.py
│   ├── run_all_figures.sh                     ← batch generate all DvT + season figures
│   ├── reorganize_results.sh                  ← one-time restructure (already run)
│   └── _archive/
│       └── 07_visualize_diversity.py          ← retired; superseded by 06_plot_diversity.py
└── herpesvirus/                               ← herpesvirus analysis raw outputs
```

---

## Figure Naming Convention

```
{marker}_r{rarefaction}_{dataset}_{type}_{palette}.png/svg
```

Examples:
- `16S_r8000_DvT_pcoa_purple.png`
- `MiFish_r17000_season_alpha_wong.svg`
- `cytb_r200_DvT_pcoa_wong.png`
- `18S_r1000_all_alpha_purple.png`

Barplots: `{marker}_unrarefied_{dataset}_barplot_{palette}.png/svg`

COD figures: `{marker}_r{depth}_cod_{type}_{palette}.png/svg`

Viral figures: `herpes_{type}_{palette}.png/svg`, `adeno_{type}_{palette}.png/svg`

---

## Rarefaction Depths

| Marker | Depth | DvT n | Notes |
|---|---|---|---|
| 16S | 8,000 reads | 13+13 | 3 samples dropped |
| MiFish | 17,000 reads | 13+13 | Unrarefied: 14+16 |
| cytb | 200 reads | 13+15 | Single-end mode; Marine excluded |
| 18S | 1,000 reads | 13+12+5 | All groups; verify n=30 |

---

## Confirmed Statistics

### 16S Alpha Diversity (Mann-Whitney U, DvT n=13+13)
| Metric | p | Sig |
|---|---|---|
| Faith's PD | 0.720 | ns |
| Shannon | 0.505 | ns |
| Observed Features | 0.682 | ns |
| Pielou's Evenness | 0.473 | ns |

### 16S Beta Diversity PERMANOVA + PERMDISP
| Metric | F | p | PERMDISP p |
|---|---|---|---|
| Unweighted UniFrac | 1.154 | 0.059 ns | 0.166 ns |
| Weighted UniFrac | 1.176 | 0.308 ns | 0.091 ns |
| Bray-Curtis | 1.427 | 0.088 ns | **0.013*** (dispersion confound — caveat) |
| Jaccard | 1.139 | **0.012*** | 0.495 ns |

### MiFish PERMANOVA (DvT n=13+13)
| Metric | F | p |
|---|---|---|
| Bray-Curtis | 2.268 | **0.017*** |
| Jaccard | 1.577 | **0.002*** |

### MiFish Alpha (DvT)
| Metric | p |
|---|---|
| Observed Features | **0.019*** |
| Shannon | 0.059 (trend) |
| Pielou's Evenness | 0.065 |

### cytb PERMANOVA (DvT n=13+15)
| Metric | F | p |
|---|---|---|
| Bray-Curtis | 1.200 | 0.085 ns |
| Jaccard | 1.185 | **0.025*** |

### cytb Alpha (DvT)
| Metric | p |
|---|---|
| Observed Features | **0.019*** |
| Shannon | **0.043*** |
| Pielou's Evenness | 0.167 ns |

### 18S PERMANOVA (all groups n=13+12+5)
| Metric | F | p |
|---|---|---|
| Bray-Curtis | 1.093 | 0.146 ns |
| Jaccard | 1.053 | **0.029*** |

### 18S Alpha (DvT)
| Metric | p |
|---|---|
| Observed Features | 0.080 ns (trend) |
| Shannon | 0.068 ns (trend) |
| Pielou's Evenness | 0.204 ns |

### MiFish COD_broad Analysis (Lead / Parasitic_Infectious / Trauma)
| Test | Metric | F | p | Notes |
|---|---|---|---|---|
| PERMANOVA overall | Jaccard | 1.255 | **0.009*** | 3-group |
| Lead vs Trauma pairwise | Jaccard | 1.555 | **0.002*** | q=0.012 — survives FDR |
| Parasitic vs Trauma | Jaccard | — | 0.029* | q=0.087 ns — does not survive FDR |
| Lead vs Parasitic | Jaccard | — | 0.207 ns | disease subgroups do not differ |
| Alpha Observed Features | — | — | **0.029*** | after removing Unknown_Other |

⚠ cytb and 18S COD_broad p-values pending (see Pending Tasks).

### Collection_source Confound Check (all markers confirmed ✓)
| Marker | Metric | PERMANOVA p | PERMDISP p | Notes |
|---|---|---|---|---|
| 16S | Bray-Curtis | 0.163 ns | **0.035*** | Dispersion differs — note limitation |
| 16S | Jaccard | 0.148 ns | ns | |
| 18S | Bray-Curtis | 0.188 ns | ns | |
| 18S | Jaccard | 0.056 ns (trend) | ns | CFW vs NHVDL pairwise p=0.040, q=0.120 |
| MiFish | Bray-Curtis | 0.196 ns | ns | |
| MiFish | Jaccard | 0.076 ns | ns | **CFW vs NHVDL pairwise q=0.048*** — note limitation |
| cytb | Bray-Curtis | 0.321 ns | ns | Cleanest result |
| cytb | Jaccard | 0.107 ns | ns | |

### Herpesvirus Detection (TGF-IYG, lung n=40)
- No-hit reads (putative novel herpesvirus): Diseased 19/19, Trauma 16/16
- Fisher's exact p=1.000 — ubiquitous detection, not disease-associated
- Figures: `results/herpesvirus/figures/`

### Aviadenovirus Detection (TGF-IYG amplicon, lung n=40) — NEW
- Classified reads matching Aviadenovirus: 2 OTUs detected at low prevalence
- Diseased: 2/19, Trauma: 1/16 — Fisher's exact p=1.000
- Top BLAST hits: ~77% identity to uncharacterized Aviadenovirus spp. (below 85% species threshold)
- Interpretation: putative novel aviadenovirus species
- Figures: `results/adenovirus/figures/`

---

## Pending Tasks

**Blocking manuscript:**
- [ ] Upload cytb COD_broad PERMANOVA QZVs and extract p-values (`results/cytb/COD/diversity/`)
- [ ] Upload 18S COD_broad PERMANOVA QZVs and extract p-values (`results/18S/COD/diversity/`)
- [ ] Verify 18S n=28 vs expected n=30 — check `metadata_18S_cod.tsv` for missing Collection_source
- [ ] Update `results.docx` — add COD findings, herpesvirus ubiquity, adenovirus novelty, collection_source caveats
- [ ] Update `abstracts_three_versions.docx` — add COD p-values and viral findings

**Near-term:**
- [ ] Adenovirus phylogenetic tree — MAFFT alignment + FastTree NJ; supports novelty claim
- [ ] `10_plot_viral.py` — add `--signal-mode`, `--taxon-filter`, `--min-reads` flags for adenovirus support (currently only supports no-hit signal mode)
- [ ] Partial PERMANOVA (adonis) controlling for Season — MiFish and cytb
- [ ] Fill manuscript placeholders — institution, extraction protocol, PCR params, citations
- [ ] Install `openpyxl` (conda clone env or pip --user + PYTHONPATH fix)
- [ ] Graduate Research Conference poster — incorporate COD + viral findings

**Completed this session (2026-03-19):**
- [x] All DvT + season diversity figures regenerated (purple + wong, standardized filenames)
- [x] COD_broad + Collection_source stats run for all markers (MiFish complete; cytb/18S QZVs need upload)
- [x] Collection_source confound check — all markers confirmed ns except noted caveats
- [x] Herpesvirus figures generated (both palettes)
- [x] Adenovirus figures generated (both palettes)
- [x] BLAST results for aviadenovirus OTUs — novelty claim supported
- [x] cytb barplot verified
- [x] MiFish barplot verified

---

## Running the Pipeline

### Generate all diversity figures (both palettes)
```bash
conda activate metabarcoding-analysis
bash scripts/run_all_figures.sh
```

### Run COD_broad + Collection_source diversity stats and figures
```bash
conda activate metabarcoding-qiime2
python scripts/05b_run_cod_diversity.py \
    --marker MiFish \
    --metrics-dir qiime2/MiFish/rarefied_17000/DvT/diversity/core-metrics-17000 \
    --metadata metadata/qiime/metadata_MiFish_cod_filtered.tsv \
    --group-column COD_broad
```

### Generate viral detection figures
```bash
conda activate metabarcoding-analysis
# Herpesvirus (no-hit signal):
python scripts/10_plot_viral.py \
    --xlsx loon_amplicon_analysis.xlsx --sheet herpes \
    --metadata metadata/qiime/metadata_16S.tsv \
    --palette purple --outdir results/herpesvirus/figures

# Adenovirus (classified signal — requires --signal-mode update to 10_plot_viral.py):
# See Pending Tasks
```

### Regenerate a single marker/dataset
```bash
python scripts/06_plot_diversity.py pcoa \
    --artifact qiime2/MiFish/rarefied_17000/DvT/diversity/core-metrics-17000/bray_curtis_pcoa_results.qza \
               qiime2/MiFish/rarefied_17000/DvT/diversity/core-metrics-17000/jaccard_pcoa_results.qza \
    --metadata metadata/qiime/metadata_MiFish.tsv \
    --color-by Group --panel --palette wong --no-title \
    --stats-dir results/MiFish/DvT/diversity \
    --output-stem MiFish_r17000_DvT_pcoa_wong \
    --output-dir results/MiFish/DvT/figures
```
