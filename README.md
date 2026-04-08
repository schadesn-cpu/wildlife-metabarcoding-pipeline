# Common Loon Gut Microbiome & Dietary Metabarcoding — Project README

**Last updated:** 2026-03-27

---

## Overview

Amplicon sequencing study of Common Loon (*Gavia immer*) gut microbiome and dietary
composition. Four markers sequenced: **16S rRNA** (bacteria), **MiFish 12S** (fish diet),
**cytochrome b** (vertebrate diet), **18S rRNA** (eukaryotes/parasites — see note below).
Viral detection via pan-herpesvirus TGF-IYG primers on lung tissue.

**Analytical framework decision (2026-03-27):** MiFish and cytb results are reported as
**presence/absence and detection frequency**, not relative abundance. Relative read
abundance is not a defensible proxy for dietary biomass because PCR efficiency differs
across prey taxa, mitochondrial copy number varies by tissue and species, and digestion
state affects DNA yield independently of consumption. See Deagle et al. (2019, *Molecular
Ecology*) for the standard treatment of this issue. cytb has the additional constraint of
rarefaction depth 200 reads, at which relative abundance estimates have variance too large
to be interpretable.

**18S status:** Excluded from the main analysis. Primer detection returned 0% for both
forward and reverse primers in primers_detected.tsv, indicating failed amplification or
wrong primers. If included in any revision it is treated as explicitly exploratory.

---

## Environment Setup

Two conda environments — one for QIIME2 pipeline steps, one for Python analysis and
plotting scripts.

### QIIME2 environment (scripts 00–05, run_all_figures.sh)
```bash
conda env create -f environment_qiime2.yml
conda activate metabarcoding-qiime2
```

### Analysis environment (scripts 06–11 and utilities)
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

### Important: Two-phase artifact layout

The project was run in two phases with different directory conventions. Know which
phase each marker belongs to before pointing scripts at paths.

**Phase 1 — flat structure (16S, 18S, ITS, original MiFish/cytb runs):**
```
qiime2/dada2/        ← table_16S.qza, rep-seqs_16S.qza, table_Mifish.qza, table_cytb.qza
qiime2/imported/     ← demux_{marker}.qza, manifest_{marker}.tsv
qiime2/taxonomy/     ← taxonomy_16S_silva138_v4.qza, taxonomy_18S.qza (flat, no subdirs)
```

**Phase 2 — nested structure (MiFish and cytb reruns — use these):**
```
qiime2/MiFish/all/dada2/table.qza
qiime2/MiFish/all/taxonomy/taxonomy.qza
qiime2/cytb/all/dada2/table.qza
qiime2/cytb/all/taxonomy/taxonomy.qza
```

**16S:** All artifacts in Phase 1 flat structure. Always use
`qiime2/taxonomy/taxonomy_16S_silva138_v4.qza` (V4-specific, 515F/806R region).
**Never** use the full-length `taxonomy_16S_silva138.qza` for this data — genus
resolution drops from ~50% to <1%.

**cytb notrim:** The cytb taxonomy results are under `results/cytb/all/taxonomy/notrim/`.
This is correct — cytb was not run through cutadapt by design. L14841/H15149 primers
(35 + 34 bp) plus the ~307 bp amplicon = ~376 bp total insert, which exceeds the 250 bp
read length. No adapter bleed-through occurs. This is documented, not a bug.

### Full project tree
```
project_root/
├── metadata/
│   └── qiime/
│       ├── metadata_{marker}.tsv              ← QIIME2 metadata (DvT analyses)
│       ├── metadata_{marker}_cod.tsv          ← with COD_broad + Collection_source
│       └── metadata_cytb_cod_filtered.tsv     ← Lead + Parasitic_Infectious + Trauma only
├── qiime2/                                    ← QIIME2 artifacts (.qza/.qzv)
│   ├── dada2/                                 ← Phase 1: 16S, 18S, ITS flat artifacts
│   ├── imported/                              ← Phase 1: demux QZAs, manifests
│   ├── taxonomy/                              ← Phase 1: taxonomy QZAs
│   ├── 16S/rarefied_8000/DvT/diversity/core_metrics_depth8000/
│   ├── MiFish/all/dada2/ + taxonomy/          ← Phase 2: use these for MiFish
│   ├── MiFish/rarefied_17000/DvT/diversity/core-metrics-17000/
│   ├── cytb/all/dada2/ + taxonomy/            ← Phase 2: use these for cytb
│   ├── cytb/all/diversity/core-metrics-200/
│   ├── 18S/all/dada2/ + taxonomy/             ← exists but 18S excluded from main analysis
│   ├── 18S/all/diversity/core-metrics-1000/
│   └── cytb/COD/                              ← filtered distance matrices for cytb COD stats
├── results/
│   ├── 16S/
│   │   ├── DvT/{diversity,figures,taxonomy}/
│   │   ├── season/{diversity,figures}/
│   │   └── COD/{diversity,figures}/
│   ├── MiFish/
│   │   ├── all/taxonomy/                      ← taxonomy_counts_L7_MiFish.tsv (source for PA)
│   │   ├── all/presence_absence/              ← 08b output (detection freq, binary table)
│   │   ├── DvT/{diversity,figures}/
│   │   ├── season/{diversity,figures}/
│   │   └── COD/{diversity,figures}/
│   ├── cytb/
│   │   ├── all/taxonomy/notrim/               ← taxonomy_counts_L7_cytb.tsv (source for PA)
│   │   ├── all/presence_absence/              ← 08b output (detection freq, binary table)
│   │   ├── DvT/{diversity,figures}/
│   │   ├── season/{diversity,figures}/
│   │   └── COD/{diversity,figures}/
│   ├── 18S/all/{diversity,figures,taxonomy}/  ← exists; excluded from main analysis
│   ├── herpesvirus/figures/
│   ├── adenovirus/figures/
│   ├── qc/                                    ← demux_qc_report.txt + .tsv (new)
│   └── _archive_old_structure/
├── reports/
│   ├── primers_detected.tsv                   ← from primer_advisor detect --all
│   └── demultiplex/
│       ├── Demultiplex_Stats.csv
│       └── additional-reports/
│           ├── Adapter_Metrics.csv
│           ├── Adapter_Cycle_Metrics.csv
│           └── multiqc_data/
│               ├── fastqc_adapter_content_plot.txt
│               └── fastqc_per_base_sequence_quality_plot.txt
└── scripts/
    ├── 00_build_classifiers.py
    ├── 00a_primer_advisor.py                  ← primer detect + suggest + check subcommands
    ├── 01_make_manifests.py
    ├── 02_make_qiime_metadata.py
    ├── 03_run_full_metabarcoding_pipeline.py
    ├── 04_rarefaction.py
    ├── 08_run_diversity_stats.py
    ├── 05b_run_cod_diversity.py               ← (fixed: sys.executable + relative path)
    ├── 09_plot_diversity.py                   ← (fixed: PERMANOVA parser no longer silent)
    ├── 07_taxonomy_table.py
    ├── 11b_presence_absence.py                ← universal PA; use for MiFish + cytb
    ├── 10_plot_taxonomy.py
    ├── 10_plot_viral.py
    ├── 11_plot_mifish_season_ecology.py
    ├── add_season_to_metadata.py              ← (fixed: bad dates now logged with sample ID)
    ├── parse_beta_stats.py
    ├── parse_multiqc_demux.py                 ← NEW: QC report from Illumina/MultiQC data
    ├── plot_adeno_tree.py                     ← (fixed: bootstrap label loop)
    ├── run_all_figures.sh
    ├── scan_files.py                          ← (fixed: PermissionError now logged)
    └── _archive/
        ├── 09c_visualize_diversity.py          ← retired; superseded by 09_plot_diversity.py
        ├── combine_multimarker_alpha.py       ← retired; use 06_ --panel instead
        ├── parse_demux_report.py              ← retired; absorbed into primer_advisor check
        ├── reorganize_loon.sh                 ← one-time migration, complete
        └── reorganize_results.sh             ← one-time migration, complete
```

---

## Pre-DADA2 QC Workflow

Three steps run before DADA2. All require the `metabarcoding-analysis` environment
(no QIIME2 needed).

### Step 1 — Detect primers from raw reads
```bash
python scripts/00a_primer_advisor.py detect \
    --all \
    --reads-dir reads/ \
    --report reports/primers_detected.tsv
```

### Step 2 — Suggest DADA2 parameters from demux QZV
```bash
python scripts/00a_primer_advisor.py suggest \
    --demux    qiime2/imported/demux_MiFish.qza \
    --marker   MiFish \
    --primer-f 21 --primer-r 27 \
    --amplicon-length 180
```

### Step 3 — QC report from Illumina demultiplex + MultiQC data
Run after cutadapt (or in this project's case, after demultiplexing when no cutadapt
logs exist). Reads the reports/ directory directly.
```bash
mkdir -p logs results/qc

python scripts/parse_multiqc_demux.py \
    --reports-dir  reports/ \
    --primers      reports/primers_detected.tsv \
    --amplicon-lens 16S=253,MiFish=180,cytb=307,Virus=110 \
    --min-reads    10000 \
    --out-txt      results/qc/demux_qc_report.txt \
    --out-tsv      results/qc/demux_qc_report.tsv
```

Note: amplicon lengths in parse_multiqc_demux.py are loon-project-specific defaults
(`LOON_PROJECT_AMPLICON_LENS`). Verify before use on any other project.

---

## Presence/Absence Analysis (MiFish and cytb)

`11b_presence_absence.py` converts taxonomy count tables to binary presence/absence
and computes detection frequencies per group. Use this — not relative abundance — for
MiFish and cytb results.

```bash
# MiFish
nohup python scripts/11b_presence_absence.py \
    --counts            results/MiFish/all/taxonomy/taxonomy_counts_L7_MiFish.tsv \
    --metadata          metadata/qiime/metadata_MiFish.tsv \
    --marker            MiFish \
    --group-by          Group \
    --min-sample-reads  10000 \
    --min-taxon-reads   10 \
    --min-relabund      0.01 \
    --sample-label      loon \
    --outdir            results/MiFish/all/presence_absence/ \
    > logs/mifish_presence_absence.log 2>&1 &

# cytb — thresholds set to match rarefaction depth of 200 reads
# --min-relabund 0.01 requires >=2 reads for a detection call at this depth
nohup python scripts/11b_presence_absence.py \
    --counts            results/cytb/all/taxonomy/notrim/taxonomy_counts_L7_cytb.tsv \
    --metadata          metadata/qiime/metadata_cytb.tsv \
    --marker            cytb \
    --group-by          Group \
    --min-sample-reads  50 \
    --min-taxon-reads   5 \
    --min-relabund      0.01 \
    --sample-label      loon \
    --outdir            results/cytb/all/presence_absence/ \
    > logs/cytb_presence_absence.log 2>&1 &
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
| 18S | 1,000 reads | 13+12+5 | Excluded from main analysis |

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

⚠ MiFish diversity stats above are from rarefied relative abundance analysis.
Presence/absence (Jaccard on binary table) will be rerun via 08b and 08_run_diversity_stats.py.
These p-values may change slightly — update this table when complete.

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

⚠ cytb diversity stats above are from rarefied relative abundance analysis.
Presence/absence rerun pending — update when complete.

### cytb Alpha (DvT)
| Metric | p |
|---|---|
| Observed Features | **0.019*** |
| Shannon | **0.043*** |
| Pielou's Evenness | 0.167 ns |

### 18S PERMANOVA (all groups n=13+12+5) — exploratory only
| Metric | F | p |
|---|---|---|
| Bray-Curtis | 1.093 | 0.146 ns |
| Jaccard | 1.053 | **0.029*** |

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

### Aviadenovirus Detection (TGF-IYG amplicon, lung n=40)
- Classified reads matching Aviadenovirus: 2 OTUs at low prevalence
- Diseased: 2/19, Trauma: 1/16 — Fisher's exact p=1.000
- Top BLAST hits: ~77% identity to uncharacterized Aviadenovirus spp. (below 85% species threshold)
- Interpretation: putative novel aviadenovirus species

---

## Pending Tasks

**Blocking manuscript:**
- [ ] Run 11b_presence_absence.py for MiFish and cytb — generate detection_freq TSVs and
      binary presence/absence tables (commands above)
- [ ] Rerun Jaccard PERMANOVA on binary presence/absence table for MiFish and cytb —
      update stats tables above with new p-values
- [ ] Run taxonomy test run for 16S (table_16S.qza + taxonomy_16S_silva138_v4.qza) —
      first time 16S taxonomy table has been generated; check output before using
- [ ] Upload cytb COD_broad PERMANOVA QZVs and extract p-values
- [ ] Upload 18S COD_broad PERMANOVA QZVs and extract p-values
- [ ] Verify 18S n=28 vs expected n=30 — check metadata_18S_cod.tsv for missing Collection_source
- [ ] Run parse_multiqc_demux.py QC report — document demux QC in results/qc/
- [ ] Update results.docx — add COD findings, herpesvirus ubiquity, adenovirus novelty,
      collection_source caveats, presence/absence framework for MiFish/cytb
- [ ] Update abstracts_three_versions.docx — add COD p-values and viral findings

**Near-term:**
- [ ] Adenovirus phylogenetic tree — MAFFT alignment + FastTree NJ; supports novelty claim
- [ ] 10_plot_viral.py — add --signal-mode, --taxon-filter, --min-reads flags for
      adenovirus support (currently only supports no-hit signal mode)
- [ ] Partial PERMANOVA (adonis) controlling for Season — MiFish and cytb
- [ ] Fill manuscript placeholders — institution, extraction protocol, PCR params, citations
- [ ] Install openpyxl (conda clone env or pip --user + PYTHONPATH fix)
- [ ] Graduate Research Conference poster — incorporate COD + viral findings
- [ ] Decide 18S fate: explicitly exploratory (state limitations) or excluded entirely

**Completed this session (2026-03-27):**
- [x] Analytical decision: MiFish and cytb → presence/absence (Deagle et al. 2019)
- [x] 18S exclusion confirmed — 0%/0% primer detection in primers_detected.tsv
- [x] parse_multiqc_demux.py written — reads Illumina demux + MultiQC data directly
- [x] 11b_presence_absence.py made universal — added --sample-label, removed dead code,
      renamed project-specific defaults
- [x] primer_advisor.py updated — added check subcommand for post-cutadapt QC
- [x] Bug fixes: 06_plot_diversity (silent except), 07_visualize_diversity (dead branch),
      05b_run_cod_diversity (hard-coded paths), add_season_to_metadata (silent date parse),
      scan_files (undifferentiated exception handler), plot_adeno_tree (exception as type probe)
- [x] Scripts archived: reorganize_loon.sh, reorganize_results.sh,
      combine_multimarker_alpha.py, parse_demux_report.py
- [x] Project structure fully audited — two-phase layout documented
- [x] Lab notebook sprint entry written (sprint_notebook.docx)

**Completed previous session (2026-03-19):**
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

# Adenovirus: requires --signal-mode update to 10_plot_viral.py (see Pending Tasks)
```

### Regenerate a single marker/dataset figure
```bash
python scripts/09_plot_diversity.py pcoa \
    --artifact qiime2/MiFish/rarefied_17000/DvT/diversity/core-metrics-17000/bray_curtis_pcoa_results.qza \
               qiime2/MiFish/rarefied_17000/DvT/diversity/core-metrics-17000/jaccard_pcoa_results.qza \
    --metadata metadata/qiime/metadata_MiFish.tsv \
    --color-by Group --panel --palette wong --no-title \
    --stats-dir results/MiFish/DvT/diversity \
    --output-stem MiFish_r17000_DvT_pcoa_wong \
    --output-dir results/MiFish/DvT/figures
```

### Taxonomy table (16S — flat structure, V4 classifier)
```bash
conda activate metabarcoding-analysis
python scripts/07_taxonomy_table.py \
    --taxonomy  qiime2/taxonomy/taxonomy_16S_silva138_v4.qza \
    --table     qiime2/dada2/table_16S.qza \
    --marker    16S \
    --outdir    results/16S/all/taxonomy/
```
