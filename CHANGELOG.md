# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.0.0] — 2026-04-14

First public release. Pipeline developed and validated on Common Loon
(*Gavia immer*) gut microbiome and dietary metabarcoding study
(Schade et al., in prep), MEED Lab, University of New Hampshire.

### Added

**Pipeline orchestration**
- `pipeline.py` — main entry point with `init`, `check`, `run`, and `figures`
  subcommands; idempotent step execution; `--dry-run` throughout
- `pipeline_config.yml` — single YAML config file; the only file that needs
  editing for a new project
- `config_loader.py` — shared config reader imported by all downstream scripts;
  no project-specific values hardcoded in any script
- `run_all_figures.py` — regenerates both annotated and manuscript figure sets
  from a single command
- `build_manuscript_figures.sh` — assembles flat `manuscript_figures/` with
  numbered figure names (Fig01–Fig13, FigS1–FigS11); checks for missing files
  and flags PENDING/REGENERATE outputs

**Reference database building** (`00_build_classifiers.py`)
- SILVA 138 (16S V4, 515F/806R), PR2 v5 (18S), UNITE v10 (ITS)
- MitoFish 12S with optional Gavia augmentation (`--add-gavia`) for loon
  studies or any species absent from the standard database
- NCBI Vertebrata cytochrome b via RESCRIPt + Entrez (separate conda env
  documented)
- MIDORI2 UNIQ NUC COI (Leray amplicon)
- BLAST nucleotide databases for adenovirus (avian adenovirus hexon) and
  herpesvirus (DPOL gene) — correct tool choice documented in code comments

**Pre-DADA2 QC**
- `04d_primer_advisor.py` — primer detection, suggestion, and post-cutadapt
  check subcommands
- `04e_parse_multiqc_demux.py` — six-question QC report from Illumina demux
  and MultiQC output; outputs inform DADA2 truncation parameters directly
- `suggest_dada2_params.py` — reads per-position quality scores from a demux
  QZV, finds Q-score dropoff via sliding window, calculates paired-end overlap,
  warns if overlap is negative (single-end required), and outputs a
  ready-to-paste DADA2 command fragment; `--out-json` for pipeline integration

**DADA2 QC**
- `05b_parse_dada2_retention.py` — per-marker bottleneck detection (filtering
  / merging / chimeras) with causal interpretation

**Metadata**
- `02_make_qiime_metadata.py` / `04_make_qiime_metadata.py` — QIIME2 metadata
  builder; handles sample ID mismatch between feature table and source metadata
  via configurable regex; supports control prefix detection
- `02b_add_season_to_metadata.py` / `04b_add_season_to_metadata.py` —
  loon ecological seasons (Breeding / Freshwater_Nonbreeding / Saltwater)
- `04c_add_meteorological_season_to_metadata.py` — standard meteorological
  seasons for studies where ecological seasons are not defined

**Taxonomy**
- `07_taxonomy_table.py` / `08_taxonomy_table.py` — marker-aware filtering;
  unrarefied filtered table for barplots (documented rationale), rarefied
  unfiltered for diversity
- `09b_clean_diet_table.py` / `11_clean_diet_table.py` — 8-step dietary table
  cleaning with full per-decision logging; marker-specific artefact lists
  including geographic plausibility flags
- `10b_annotate_diet_ecology.py` / `11c_annotate_diet_ecology.py` — adds
  habitat, trophic role, and lay-friendly common group annotations; fully
  adaptable to non-loon study systems via external lookup table

**BLAST verification** (three-stage pipeline)
- `07c_blast_qc_unclassified.py` / `09a_blast_qc_unclassified.py` — BLAST QC
  for ASVs unresolved below class level; catches host DNA miscalled as target
- `07b_blast_verify.py` / `08c_blast_verify.py` — verification of suspect or
  high-abundance taxa by name substring or read-count threshold; AGREE /
  DISAGREE / NO HIT / ARTEFACT report
- Outputs feed exclusion list in `11_clean_diet_table.py`

**Presence/absence framework**
- `08b_presence_absence.py` / `11b_presence_absence.py` — universal PA
  framework with methodological citations (Deagle et al. 2019; Borland &
  Kading 2021); `--sample-label` for study-specific terminology

**Diversity statistics**
- `05_run_diversity_stats.py` / `08_run_diversity_stats.py` — Kruskal-Wallis
  alpha + PERMANOVA + PERMDISP beta, pairwise, phylo-aware
- `05b_run_cod_diversity.py` / `08b_run_cod_diversity.py` — cause-of-death
  stratified diversity analysis
- `05c_parse_beta_stats.py` / `08c_parse_beta_stats.py` — parses QIIME2
  PERMANOVA/PERMDISP HTML from QZV ZIP archives into clean TSV with
  significance stars

**Viral statistics**
- `14_viral_stats.py` / `15_viral_stats.py` — Fisher's exact (pairwise) and
  chi-square (3+ groups) for viral detection rates; includes chi-square
  validity check (min expected cell warning)

**Visualization**
- `06_plot_diversity.py` / `09_plot_diversity.py` — PCoA and alpha diversity
  plots; reads QZA files directly without QIIME2 installation; PNG (300 dpi)
  + SVG (Illustrator-editable)
- `09_plot_taxonomy.py` / `10_plot_taxonomy.py` — stacked barplots with
  italic taxonomy labels, group dividers, and n= annotations
- `06b_combine_multimarker_alpha.py` / `09b_combine_multimarker_alpha.py` —
  4-panel multi-marker alpha diversity comparison figure
- `12_plot_habitat_season.py` / `13_plot_habitat_season.py` — prey habitat
  composition by ecological season; serves as pipeline validation figure
- `plot_mifish_ecology_pa.py` — 2-panel detection frequency by ecological
  season and cause-of-death
- `10_plot_viral.py` / `14_plot_viral.py` — relative abundance and
  presence/absence barplots for viral markers; Fisher's exact p annotation
- `10b_plot_herpes_cutadapt.py` / `14c_plot_herpes_cutadapt.py` — herpesvirus
  detection from cutadapt primer-confirmed reads (replaces invalid DADA2
  no-hit approach; rationale documented in script header)
- `10b_plot_adeno_tree.py` / `14b_plot_adeno_tree.py` — publication-quality
  IQ-TREE phylogenetic tree via BioPython; outgroup-rooted, bootstrap
  suppressed below 80%, Wong 2011 colorblind-safe palette throughout
- `11_plot_mifish_season_ecology.py` — 3-panel seasonal ecology figure with
  chi-square confound test

**Environment**
- `environment.yml` — single conda environment covering QIIME2 2024.5 +
  all analysis and plotting dependencies (MEED Lab, UNH)

### Key analytical decisions documented

- MiFish 12S and cytb results reported as presence/absence and detection
  frequency (Deagle et al. 2019 rationale; validated by Clucas et al. 2024
  for Gulf of Maine piscivores)
- Herpesvirus detection method corrected from DADA2 no-hit reads (invalid) to
  cutadapt primer-confirmed read pairs; BLAST-validated at 91–97% identity to
  *Gallid alphaherpesvirus 1* DPOL
- Rarefied vs. unrarefied table choice is explicit and documented: unrarefied
  filtered for taxonomy barplots; rarefied unfiltered for diversity statistics
- Extinct species (*Pinguinus impennis*) explicitly listed in artefact filter
  with documented rationale

### Bug fixes (from internal development history)

- `06_plot_diversity.py`: silent except clause replaced with logged warning
- `05b_run_cod_diversity.py`: hardcoded HPC paths replaced with config-driven
  relative paths
- `add_season_to_metadata.py`: unparseable dates now logged with sample ID
  instead of silently assigned empty season
- `plot_adeno_tree.py`: bootstrap label detection changed from `float()` as
  type probe with bare except to explicit numeric string check
- `scan_files.py`: undifferentiated `Exception` handler replaced with specific
  `PermissionError` logging
- `09_plot_diversity.py`: alpha diversity annotations now read Kruskal-Wallis
  H-statistic and p/q-values from QZV pairwise CSV; adds `--stats-qzv-dir`
  argument (replaces scipy Mann-Whitney U recomputation)
- `11b_presence_absence.py`: zero-detection taxa filter applied consistently
  to overall barplot path; group barplot labels now parsed from full taxonomy
  strings into clean binomials via `shorten_taxon_label()`; Wong palette
  applied to `GROUP_COLORS`
- `13_plot_habitat_season.py`: `GROUP_COLORS` updated to Wong 2011; MiFish
  figure now uses nobact-filtered counts table (691K reads) instead of
  annotated table (5.6M reads, 87% bacterial contamination)
- `plot_mifish_ecology_pa.py`: `Unknown_Other` COD category removed from
  Panel B; Wong-colored group identity strip added below x-axis in both panels

---

## Notes on version history

This pipeline was developed iteratively during the loon study (2024–2026).
Script numbering reflects the order of execution in the analysis workflow,
not a version number. Scripts with the same number and different letter
suffixes (e.g. `08b`, `08c`) are companions to the main numbered script,
not replacements.

The `_archive/` directory in `scripts/` contains retired scripts superseded
by the versions documented above.

[1.0.0]: https://github.com/YOUR_USERNAME/YOUR_REPO/releases/tag/v1.0.0
