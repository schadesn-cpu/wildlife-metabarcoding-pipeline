# Changelog

All notable analysis decisions, parameter changes, and script updates are recorded here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

---

## [1.2.0] — 2026-03-08 — Classifiers trained + barplots validated end-to-end

### Classifiers built

**16S — SILVA 138 V4-specific classifier**
- Downloaded pre-trimmed 515F/806R SILVA 138 seqs and trained `silva-138-99-nb-classifier-515-806.qza` (~57 min on ron.sr.unh.edu)
- Result: genus resolution jumped from <1% (full-length classifier) to ~21%; 335 genera, 28/30 samples retained
- Top genera: Cetobacterium 17.8%, uncl. Enterobacteriaceae 13.7%, Enterococcus 11.8% — consistent with fish-eating waterbird GI microbiome

**MiFish — Gavia-augmented classifier**
- Initial NCBI fetch returned 3 partial 12S seqs from wrong region; cutadapt primer extraction returned 0 sequences
- Root cause confirmed: MiFish-U primers have 4-6 mismatches with bird 12S primer binding sites — no extraction method works at standard thresholds
- BLAST confirmed dominant ASV (47dfe1c3, 4.3M reads, 34 samples) is Gavia immer (>99% identity, placed within Gavia clade in BLAST tree)
- Fix: appended that ASV sequence directly to already-trimmed MitoFish FASTA, bypassing extract-reads entirely
- Trained `mitofish-12S-mifish-gavia-classifier.qza`; Gavia reads now correctly excluded, remaining uncl. Actinopteri is genuine unclassified fish prey

### Pipeline runs completed

- **16S `08_`**: 335 genera × 28 samples → `results/16S/DvT/taxonomy/`
- **MiFish `08_`**: 48 species × 34 samples → `results/MiFish/all/taxonomy/` (TV240036-GI, TV240100-GI dropped, <100 reads after host removal)
- **`09_`** both markers: `results/16S/DvT/figures/taxonomy/barplot_16S_Group_purple.png` and `results/MiFish/all/figures/taxonomy/barplot_MiFish_Group_purple.png`
- Both figures validated against Excel barplots — genera and proportions consistent

### Bug fixes

**`08_taxonomy_table.py` — `load_feature_table_qza` (portability)**
- Was: called `qiime tools export` + `biom convert` via subprocess — required QIIME2 env for pure analysis tasks
- Fixed: reads QZA directly via `zipfile` + `biom` Python API; now works in `wildlife-metabarcoding-analysis` without QIIME2

### Environment

- Created `envs/wildlife-metabarcoding-analysis.yml` — lightweight env for analysis scripts, no QIIME2 required
- `qiime2-amplicon-2024.5` is a shared system install at `/home/share/anaconda/envs/`; cannot be modified, export only
- Use `wildlife-metabarcoding-analysis` for `08_`, `09_`, diversity plot scripts; use `qiime2-amplicon-2024.5` for all `qiime` commands

### Updates — `00_build_classifiers.py` (Gavia fetch rewrite)

- `fetch_gavia_seqs_ncbi()` now fetches complete mitogenomes; removed `>20000 bp` filter that silently skipped all mitogenomes
- cutadapt now runs with `-e 0.3` (30% error, ~6 mismatches) to handle fish→bird primer site divergence
- Added hardcoded fallback for *Gavia immer*: BLAST-confirmed 199bp amplicon from this study; used if cutadapt fails
- Empty output FASTA now causes immediate hard abort with actionable error message
- `build_MiFish --add-gavia`: removed `extract-reads` step entirely; reads already-trimmed MitoFish QZA via zipfile, appends pre-trimmed Gavia amplicons, trains directly

### Decision log additions

**Why uncl. Actinopteri persists after Gavia augmentation**
The remaining 71.6% "uncl. Actinopteri" is genuine unclassified fish prey — MitoFish has incomplete species coverage for many coastal fish. Species that ARE resolved (Brevoortia, Fundulus, Alosa, Notemigonus, Urophycis, Prionotus) are consistent with loon coastal/estuarine diet.

**Why MiFish primers fail to extract Gavia sequences**
MiFish-U was designed for Actinopterygii; 4-6 mismatches with Aves in primer binding sites. cutadapt -e 0.1 and QIIME2 extract-reads at identity=0.8 both return 0 sequences from complete Gavia mitogenomes. Pre-trimmed BLAST-confirmed amplicons are the correct reference source.

**Why the Gavia reference uses this study's rep-seqs**
The dominant uncl. Actinopteri ASV (47dfe1c3) is pre-trimmed to the exact amplicon region, matches the sequencing platform and primers used here, and was confirmed >99% identity to G. immer by BLAST phylogenetic placement. Using an observed ASV as reference is more appropriate than a database partial sequence for this application.

---

## [1.1.0] — 2026-03-07 — Taxonomy barplots + classifier fixes

### Scripts added
- `08_taxonomy_table.py` — collapses QIIME2 feature table to relative abundance TSV
  at a given taxonomic level, with include/exclude filtering and per-marker defaults
- `09_plot_taxonomy.py` — publication-quality stacked barplots from relabund TSV;
  supports purple/redblue/wong palettes, group ordering, top-N with Other collapse

### Bug fixes — `08_taxonomy_table.py`

**1. `clean_taxon_label` walk-back range (critical)**
- Was: `range(level-2, 0, -1)` — stopped before index 0 (Domain), so reads
  classified only to Bacteria/Eukaryota returned `None` and were silently dropped
- Fixed: `range(level-2, -1, -1)` — now returns e.g. "uncl. Bacteria" instead
  of dropping. Affected 16S samples with few genus-level reads most severely.

**2. MiFish `prefix_strip` default (critical)**
- Was: `False` — QIIME2 rank prefixes (`d__`, `p__`, `c__`, `s__`) were left in
  labels, producing display strings like "s  Brevoortia patronus"
- Fixed: `True` — prefixes stripped at label-parse time; labels now clean

**3. MiFish `MARKER_FILTER_DEFAULTS` include term (critical)**
- Was: `"Vertebrata"` — matches 0 features in the MitoFish QIIME2 DB, which uses
  `d__Eukaryota;p__Chordata;c__Actinopteri` (not Vertebrata/Actinopterygii)
- Fixed: `"Actinopteri"` — matches all 1024 features

**4. f-string in `drop_empty_samples` warning (cosmetic)**
- Was: f-string used `<{min_reads}>` literal instead of the variable value
- Fixed: uses `%d` formatting to show the actual threshold number

### Diagnosis — classifier root causes

**MiFish: Gavia not in reference DB**
- 81% of MiFish reads fall through to "uncl. Actinopteri" because the MitoFish
  QIIME2 DB (Mitohelper Mar 2025, 1053 features) contains zero Gavia sequences
- Host reads cannot be labeled or excluded without Gavia in the reference
- Fix: `00_build_classifiers.py --markers MiFish --add-gavia` — fetches all 5
  Gavia species from NCBI, merges with MitoFish, retrains classifier
- Output: `mitofish-12S-mifish-gavia-classifier.qza`
- After retraining, `08_` default exclude filter ("Aves") cleanly removes host

**16S: full-length SILVA classifier used on V4 amplicons**
- `silva-138-99-nb-classifier.qza` (209 MB) is pre-trained on full-length 16S
- Applied to V4 amplicon reads (~253 bp), genus resolution is <1%
  (~99% of bacterial reads fall to "d__Bacteria;__;__;__;__;__")
- The genome center's classifier was trained on the V4 region, giving 40-60%
  genus resolution — that is why their Excel barplots had informative genera
- Fix: `00_build_classifiers.py --markers 16S --16s-v4` — downloads pre-trimmed
  515F/806R SILVA 138 seqs (~14 MB) and trains a V4-specific classifier
- Output: `silva-138-99-nb-classifier-515-806.qza`
- Training time: ~2-4 hours (submit as cluster job)

### Updates — `00_build_classifiers.py`

- **`--16s-v4` flag**: downloads `silva-138-99-seqs-515-806.qza` (pre-trimmed to
  515F/806R region) and trains a V4 classifier. Replaces broken default behaviour
  of downloading the full-length pre-trained classifier for V4 data.
- **`--add-gavia` flag**: fetches Gavia immer, G. stellata, G. arctica,
  G. pacifica, G. adamsii from NCBI nuccore; imports, merges with MitoFish,
  re-extracts MiFish amplicon region, retrains. Requires `--ncbi-email`.
- **Removed `--p-classify--n-jobs`**: this parameter was removed from
  `qiime feature-classifier fit-classifier-naive-bayes` in QIIME2 2024.5.
  Training now uses all available cores automatically. The old flag caused
  an immediate error and prevented any classifier from being trained.
- Fixed QIIME2 data URLs from `2024.10` to `2024.5` to match server environment.

### Decision log additions

**Why current 16S barplots were generated from Excel, not QZA files**
The SILVA 138 classifier on the server is the full-length pre-trained version,
which gives <1% genus resolution on V4 amplicons. The genome center used a
V4-trimmed SILVA classifier. Until `silva-138-99-nb-classifier-515-806.qza` is
trained and taxonomy is re-run, the Excel-derived barplots are used for figures.

**Why MiFish barplots were generated from Excel, not QZA files**
The MitoFish QIIME2 DB has no Gavia sequences. Host loon reads classify to
"uncl. Actinopteri" and dominate all samples (~81% mean). The Gavia-augmented
classifier must be built and taxonomy re-run before QZA-derived MiFish barplots
are usable. Excel barplots used in the interim.

---

## [1.0.0] — Manuscript submission

### Analysis — 16S rRNA (V4)

**Samples**
- 26 samples total: Diseased (n=13), Trauma (n=13)
- 4 seasonal groups: Summer (n=14), Fall (n=6), Winter (n=4), Spring (n=1)
- 1 sample excluded from seasonal analysis: TV240106 (collection date unknown)

**DADA2 denoising**
- Paired-end denoising via `qiime dada2 denoise-paired`
- Classifier: SILVA 138 99% OTUs, pre-trained full-length
- Taxonomy filter: excluded mitochondria, chloroplast, Eukaryota, Archaea

**Rarefaction**
- Chosen depth: 8,000 reads/sample
- All 26 samples retained at this depth

**Beta diversity**
- Metrics: Unweighted UniFrac, Weighted UniFrac, Bray-Curtis, Jaccard
- PERMANOVA: 999 permutations, seed=42
- PERMDISP: 999 permutations, seed=42
- Key finding: Bray-Curtis PERMDISP significant (F=6.089, p=0.012) — groups
  differ in compositional variance, not centroid position
- Unweighted UniFrac (F=1.177, p=0.050) and Jaccard (F=1.134, p=0.012)
  significant by PERMANOVA

**Alpha diversity**
- Metrics: Faith's PD, Shannon entropy, Observed Features, Pielou's Evenness
- Statistical test: Mann-Whitney U (two-sided)
- All metrics non-significant (Diseased vs. Trauma): p > 0.05

**Differential abundance**
- Method: ANCOM-BC (qiime composition ancombc), Trauma vs. Diseased reference
- Taxonomic level: genus (SILVA level 6)
- Result: No genera significant after FDR correction (all q ≥ 0.05)
- Largest effect sizes (non-significant): Escherichia-Shigella (LFC=−0.582),
  Cetobacterium (LFC=−0.437), Bacteroides (LFC=+0.289)

### Scripts added
- `00_build_classifiers.py` — classifier training for 16S, 18S, ITS, MiFish, cytb
- `01_make_manifests.py` — auto-build QIIME 2 manifests from reads directory
- `02_make_qiime_metadata.py` — generate QIIME metadata TSV matching table sample IDs
- `03_run_full_metabarcoding_pipeline.py` — full marker-aware QIIME 2 pipeline wrapper
- `04_rarefaction.py` — standalone rarefaction curve generator and depth advisor
- `05_plot_diversity.py` — publication figures (PNG + SVG, colorblind palette)
- `06_visualize_diversity.py` — publication figures with PERMANOVA/PERMDISP overlays
- `parse_beta_stats.py` — standalone utility: parse QZV files into stats summary TSV
- `run_qiime_marker_pipeline.py` — post-DADA2 diversity pipeline with auto stats summary

---

## Decision log

Decisions recorded here explain *why* key parameters were chosen, for methods
transparency and reviewer response.

**Why 8,000 reads rarefaction depth?**
8,000 reads represents the 10th percentile of per-sample read counts, retaining
all 26 samples. Rarefaction curves plateaued well before this depth for all
samples, indicating adequate diversity capture.

**Why exclude TV240106 from seasonal analysis?**
Collection date was unknown, preventing seasonal group assignment. Sample was
retained in the Diseased vs. Trauma analysis where date is not a grouping variable.

**Why ANCOM-BC instead of ANCOM?**
ANCOM-BC accounts for unequal sampling fractions, is more robust with small
sample sizes (n=13 per group), and produces interpretable log fold-change
estimates with FDR-corrected q-values. The legacy `qiime composition ancom`
command was not used for final analyses.

**Why report PERMDISP alongside PERMANOVA?**
Bray-Curtis PERMDISP was significant (p=0.012), indicating the groups differ
in compositional variance rather than centroid position. Reporting both tests
allows readers to distinguish true compositional separation from
heterogeneity of dispersion.

**Why filter mitochondria, chloroplast, Eukaryota, and Archaea from 16S table?**
Standard practice for 16S bacterial community analysis. These sequences
represent host/organelle contamination and off-target amplification that
would inflate diversity estimates and confound bacterial community comparisons.
