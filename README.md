# Metabarcoding Pipeline

A reproducible, marker-aware QIIME 2 pipeline for multi-marker metabarcoding projects. Designed for projects combining 16S rRNA (bacteria), 18S rRNA (eukaryotes), ITS (fungi), MiFish 12S (fish), and cytochrome b (vertebrates).

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Pipeline Sequence](#pipeline-sequence)
- [Script Reference](#script-reference)
  - [00_build_classifiers.py](#00_build_classifierspy)
  - [01_make_manifests.py](#01_make_manifestspy)
  - [02_make_qiime_metadata.py](#02_make_qiime_metadatapy)
  - [03_run_full_metabarcoding_pipeline.py](#03_run_full_metabarcoding_pipelinepy)
  - [04_rarefaction.py](#04_rarefactionpy)
  - [06_plot_diversity.py](#05_plot_diversitypy)
  - [07_visualize_diversity.py](#06_visualize_diversitypy)
  - [run_qiime_marker_pipeline.py](#run_qiime_marker_pipelinepy)
- [QIIME 2 Metadata Format](#qiime-2-metadata-format)
- [Choosing a Rarefaction Depth](#choosing-a-rarefaction-depth)
- [Directory Structure](#directory-structure)
- [Color Palettes](#color-palettes)
- [Notes and Known Limitations](#notes-and-known-limitations)

---

## Overview

This pipeline wraps QIIME 2 commands into modular Python scripts with consistent argparse interfaces, logging, dry-run support, and marker-aware defaults. Each script is numbered to reflect its position in the workflow. Scripts can be run independently or as part of an end-to-end pipeline.

All QIIME artifacts are scoped by `{marker}/{dataset}` to prevent cross-contamination between markers and analysis subsets (e.g., `16S/all`, `16S/DvT`, `ITS/all`).

---

## Requirements

### QIIME 2

Activate a QIIME 2 conda environment before running any script that calls `qiime`:

```bash
conda activate qiime2-amplicon-2024.5
```

Verify:
```bash
qiime info
```

### Python dependencies

All scripts require Python ≥ 3.8. Additional packages for visualization:

```bash
pip install matplotlib numpy scipy pandas
# or
conda install -c conda-forge matplotlib numpy scipy pandas
```

For YAML config file support in `03_run_full_metabarcoding_pipeline.py`:

```bash
pip install pyyaml
```

### Optional (cytb only)

The RESCRIPt QIIME 2 plugin is required to build a cytochrome b classifier. `00_build_classifiers.py` will attempt to install it automatically, or install manually:

```bash
conda install -c conda-forge -c bioconda -c qiime2 -c defaults xmltodict rescript
```

---

## Pipeline Sequence

```
reads/
├── 16S/
├── ITS/
└── ...
         │
         ▼
00_build_classifiers.py      ← Download / train taxonomic classifiers
         │
         ▼
01_make_manifests.py         ← Build QIIME manifests from reads directory
         │
         ▼
02_make_qiime_metadata.py    ← Create QIIME metadata TSV matching table sample IDs
         │
         ▼
03_run_full_metabarcoding_pipeline.py
    ├── import                ← Import reads via manifest
    ├── dada2                 ← Denoise (DADA2 paired-end)
    ├── taxonomy              ← Assign taxonomy (sklearn classifier)
    ├── filter                ← Remove off-target taxa
    ├── collapse              ← Collapse to taxonomic level
    ├── taxa-plots            ← Barplot visualizations
    ├── diversity             ← Core metrics + PERMANOVA/PERMDISP
    └── diff-ancom            ← ANCOM differential abundance
         │
         ▼
04_rarefaction.py            ← Inspect rarefaction curves; choose sampling depth
         │
         ▼
06_plot_diversity.py         ← Publication figures: PNG + SVG, colorblind palette
                                (independent of 06 — run either or both)
07_visualize_diversity.py    ← Publication figures: purple/red-blue palettes,
                                PERMANOVA + PERMDISP stats overlaid on PCoA
```

> **Important:** Run `04_rarefaction.py` before committing to a `--sampling-depth` in `03_`. The rarefaction script helps you find an appropriate threshold that balances diversity capture and sample retention.

---

## Script Reference

### `00_build_classifiers.py`

Downloads reference databases and trains QIIME 2 naive Bayes classifiers for one or more markers.

| Marker | Database |
|--------|----------|
| 16S | SILVA 138 (pre-trained or region-specific) |
| 18S | PR2 v5.0.0 |
| ITS | UNITE v10 (requires manual download — see note below) |
| MiFish | MARES v2 nobar (12S rRNA) |
| cytb | NCBI vertebrate cytochrome b via RESCRIPt |

**UNITE note:** UNITE requires a manual download from https://unite.ut.ee/repository.php. Download the QIIME release (developer version, all eukaryotes, dynamic), extract, and place the `.fasta` and `.txt` files in your `--outdir` before running.

```bash
# Download pre-trained 16S + train ITS (after manual UNITE download)
python 00_build_classifiers.py \
    --markers 16S ITS \
    --outdir classifiers/

# Train region-specific 16S classifier with primers
python 00_build_classifiers.py \
    --markers 16S \
    --f-primer GTGYCAGCMGCCGCGGTAA \
    --r-primer GGACTACNVGGGTWTCTAAT \
    --outdir classifiers/

# cytb with NCBI credentials
python 00_build_classifiers.py \
    --markers cytb \
    --ncbi-api-key YOUR_KEY \
    --ncbi-email your@email.edu \
    --outdir classifiers/
```

---

### `01_make_manifests.py`

Auto-builds QIIME 2 `PairedEndFastqManifestPhred33V2` manifest TSVs from a reads directory organized into per-marker subdirectories.

**Expected reads directory structure:**
```
reads/
├── 16S/
│   ├── TV230084-GI-16S_S1492_L002_R1_001.fastq.gz
│   ├── TV230084-GI-16S_S1492_L002_R2_001.fastq.gz
│   └── ...
├── 18S/
├── ITS/
└── MiFish/
```

The script auto-detects R1/R2 pairs, resolves symlinks to absolute paths (important on HPC systems), strips lane/index suffixes to derive clean sample IDs, validates every sample has both R1 and R2, and reports unpaired files before writing.

```bash
# Auto-detect all markers
python 01_make_manifests.py \
    --reads-dir reads/ \
    --outdir qiime2/imported/

# Specific markers only
python 01_make_manifests.py \
    --reads-dir reads/ \
    --outdir qiime2/imported/ \
    --markers 16S ITS

# Preview without writing
python 01_make_manifests.py \
    --reads-dir reads/ \
    --outdir qiime2/imported/ \
    --dry-run --verbose
```

**Output:** One TSV per marker at `qiime2/imported/manifest_<marker>.tsv`.

---

### `02_make_qiime_metadata.py`

Builds a QIIME 2-ready metadata TSV whose `sample-id` column exactly matches the sample IDs present in a QIIME 2 feature table.

**Problem it solves:** Feature table sample IDs look like `TV230084-GI-16S_S1492` but your source metadata uses short IDs like `TV230084`. QIIME 2 requires exact matches; this script extracts the join key from the table ID and maps it to the source metadata.

```bash
python 02_make_qiime_metadata.py \
    --table qiime2/dada2/table_16S.qza \
    --source-metadata metadata/source_metadata.tsv \
    --source-id-column TV \
    --out metadata/qiime/metadata_16S.tsv \
    --marker 16S \
    --mapping-report metadata/qiime/mapping_report_16S.tsv
```

Key options:
- `--key-regex` — regex with one capture group to extract the join key from table sample IDs (default: `(TV\d+)`)
- `--control-prefixes` — prefixes that identify negative/positive controls (default: `NTC-,PAC-,XB-`)
- `--mapping-report` — optional TSV showing which sample IDs matched or failed to match

See the [QIIME 2 Metadata Format](#qiime-2-metadata-format) section for exact formatting requirements.

---

### `03_run_full_metabarcoding_pipeline.py`

The main pipeline runner. All outputs are scoped to `{marker}/{dataset}` to keep analyses organized.

**Global flags (required for all subcommands):**

```bash
python 03_run_full_metabarcoding_pipeline.py \
    --marker 16S \
    --dataset DvT \
    --metadata metadata/qiime/metadata_16S.tsv \
    [subcommand] [subcommand options]
```

**Subcommands:**

| Subcommand | Description |
|------------|-------------|
| `init` | Create directory layout for this marker/dataset |
| `import` | Import reads from manifest into a demux artifact |
| `dada2` | Run DADA2 paired-end denoising |
| `taxonomy` | Assign taxonomy with sklearn classifier |
| `filter` | Remove off-target taxa (marker-aware defaults) |
| `collapse` | Collapse feature table to a taxonomic level |
| `taxa-plots` | Generate ASV and filtered barplot QZVs |
| `diversity` | Core metrics (alpha + beta) + optional PERMANOVA/PERMDISP |
| `diff-ancom` | ANCOM differential abundance |
| `export` | Export any artifact with `qiime tools export` |
| `bundle` | Bundle metadata + artifacts into a tar.gz |
| `barplot-asv` | Canonical ASV-level taxa barplot |
| `smoke-test` | Validate environment and expected artifacts |
| `run` | One-liner: init → import → dada2 → taxonomy → diversity → barplot |

**Marker-aware taxonomy filter defaults** (applied automatically unless overridden):

| Marker | Exclude | Include |
|--------|---------|---------|
| 16S | mitochondria, chloroplast, Eukaryota, Archaea | — |
| 18S | Bacteria, Archaea, Viruses | — |
| ITS / ITS1-2 | Bacteria, Archaea, Viruses | Fungi |
| MiFish | Bacteria, Archaea, Viruses | Actinopterygii |
| cytb | Bacteria, Archaea, Viruses | Vertebrata |

**Example — full run:**

```bash
python 03_run_full_metabarcoding_pipeline.py \
    --marker 16S --dataset DvT \
    --metadata metadata/qiime/metadata_16S.tsv \
    run \
    --manifest    qiime2/imported/manifest_16S.tsv \
    --classifier  classifiers/silva-138-99-nb-classifier.qza \
    --tree        qiime2/diversity/rooted-tree_16S.qza \
    --sampling-depth 8000 \
    --group-column   Group \
    --trunc-len-f 240 --trunc-len-r 200
```

**Example — diversity only (after DADA2 is complete):**

```bash
python 03_run_full_metabarcoding_pipeline.py \
    --marker 16S --dataset DvT \
    --metadata metadata/qiime/metadata_16S.tsv \
    diversity \
    --sampling-depth 8000 \
    --tree qiime2/diversity/rooted-tree_16S.qza \
    --group-column Group
```

**Dry-run mode:** Add `--dry-run` to any command to print what would be executed without running it.

**Note on `diff-ancom`:** This subcommand uses the legacy `qiime composition ancom`. For new analyses, ANCOM-BC (`qiime composition ancombc`) is preferred as it accounts for sampling fraction and produces log fold-change estimates with FDR-corrected q-values. ANCOM-BC support will be added in a future version.

---

### `04_rarefaction.py`

**Run this before the main pipeline** to explore rarefaction curves and choose an appropriate `--sampling-depth`.

```bash
# With rooted tree (recommended — includes Faith's PD)
python 04_rarefaction.py \
    --table    qiime2/dada2/table_16S.qza \
    --tree     qiime2/diversity/rooted-tree_16S.qza \
    --metadata metadata/qiime/metadata_16S.tsv \
    --max-depth       20000 \
    --candidate-depth 8000 \
    --group-column    Group \
    --marker          16S \
    --outdir          results/rarefaction/

# Without tree (skips Faith's PD)
python 04_rarefaction.py \
    --table    qiime2/dada2/table_ITS.qza \
    --metadata metadata/qiime/metadata_ITS.tsv \
    --no-tree \
    --max-depth 10000 \
    --marker    ITS \
    --outdir    results/rarefaction/
```

**Outputs:**
- `alpha_rarefaction_<marker>.qzv` — interactive QIIME visualization (open at [view.qiime2.org](https://view.qiime2.org))
- `rarefaction_summary_<marker>.png` — annotated static figure showing read count distribution and retention curve
- `rarefaction_depth_report_<marker>.tsv` — per-sample read counts and drop status at candidate depth

**Terminal output** includes a full depth selection report listing samples that would be dropped, percentages retained, and the 10th and 25th percentile depths.

See [Choosing a Rarefaction Depth](#choosing-a-rarefaction-depth) for guidance on interpreting the output.

---

### `06_plot_diversity.py`

Generate publication-ready PCoA and alpha diversity figures from QIIME 2 QZA files. Outputs both **PNG** (for slides/SharePoint) and **SVG** (for Illustrator editing) for every figure. No QIIME 2 installation required.

This script uses a **colorblind-friendly palette** (Wong 2011) and automatically handles the common case where artifact sample IDs are long (e.g., `TV230084-GI-16S_S1492`) but metadata uses short IDs (e.g., `TV230084`) — no pre-processing needed.

**Subcommands:**

| Subcommand | Description |
|------------|-------------|
| `pcoa` | 2D PCoA scatter plots with 95% confidence ellipses |
| `alpha` | Alpha diversity strip + box plots with significance brackets |
| `columns` | List available metadata columns and their unique values |

```bash
# Check what columns are available in your metadata
python 06_plot_diversity.py columns \
    --metadata metadata/qiime/metadata_16S.tsv

# PCoA panel — all four beta metrics in one figure
python 06_plot_diversity.py pcoa \
    --artifact \
        qiime2/diversity/core_metrics_16S/weighted_unifrac_pcoa_results.qza \
        qiime2/diversity/core_metrics_16S/unweighted_unifrac_pcoa_results.qza \
        qiime2/diversity/core_metrics_16S/bray_curtis_pcoa_results.qza \
        qiime2/diversity/core_metrics_16S/jaccard_pcoa_results.qza \
    --metadata    metadata/qiime/metadata_16S.tsv \
    --color-by    Group \
    --panel \
    --output-dir  results/figures/

# Individual PCoA figures (one file per metric)
python 06_plot_diversity.py pcoa \
    --artifact qiime2/diversity/core_metrics_16S/weighted_unifrac_pcoa_results.qza \
    --metadata metadata/qiime/metadata_16S.tsv \
    --color-by Group \
    --output-dir results/figures/

# Alpha diversity panel — all four metrics in one figure
python 06_plot_diversity.py alpha \
    --artifact \
        qiime2/diversity/core_metrics_16S/faith_pd_vector.qza \
        qiime2/diversity/core_metrics_16S/shannon_vector.qza \
        qiime2/diversity/core_metrics_16S/observed_features_vector.qza \
        qiime2/diversity/core_metrics_16S/evenness_vector.qza \
    --metadata   metadata/qiime/metadata_16S.tsv \
    --group-by   Group \
    --panel \
    --output-dir results/figures/
```

**Statistical annotations:**
- 2 groups: Mann-Whitney U p-value + significance bracket and stars (`*`, `**`, `***`)
- 3+ groups: Kruskal-Wallis p-value
- Disable with `--no-stats`

**Relationship to `07_visualize_diversity.py`:** These two scripts are **independent** — neither depends on output from the other. Use `05_` for colorblind-friendly figures with PNG+SVG output and automatic short-ID matching; use `06_` when you need the project-specific purple/red-blue palettes and PERMANOVA/PERMDISP statistics overlaid directly on PCoA panels.

---

### `07_visualize_diversity.py`

Generate publication-quality alpha and beta diversity figures directly from QIIME 2 QZA files. **No QIIME 2 installation required** — QZA files are ZIP archives that this script reads directly using only Python standard library + matplotlib/numpy/scipy.

```bash
python 07_visualize_diversity.py \
    --metadata     metadata/qiime/metadata_16S.tsv \
    --group-column Group \
    --marker       16S \
    --alpha-qza \
        qiime2/diversity/core_metrics_16S/faith_pd_vector.qza \
        qiime2/diversity/core_metrics_16S/shannon_vector.qza \
        qiime2/diversity/core_metrics_16S/observed_features_vector.qza \
        qiime2/diversity/core_metrics_16S/evenness_vector.qza \
    --beta-pcoa \
        qiime2/diversity/core_metrics_16S/unweighted_unifrac_pcoa_results.qza \
        qiime2/diversity/core_metrics_16S/weighted_unifrac_pcoa_results.qza \
        qiime2/diversity/core_metrics_16S/bray_curtis_pcoa_results.qza \
        qiime2/diversity/core_metrics_16S/jaccard_pcoa_results.qza \
    --stats-tsv \
        qiime2/diversity/core_metrics_16S/beta_stats_summary_Group.tsv \
    --palette purple \
    --outdir  results/figures/
```

If `--stats-tsv` is not provided but `--beta-dm` distance matrix QZAs are, PERMANOVA and PERMDISP are computed on-the-fly. If neither is provided, the stats box shows `n/a`.

**Statistical tests:**
- Alpha diversity: Mann-Whitney U (2 groups) or Kruskal-Wallis (3+ groups)
- Beta diversity: PERMANOVA + PERMDISP from pre-computed TSV or computed from distance matrices

**Palettes:** `purple` (default) or `redblue`. See [Color Palettes](#color-palettes).

**Group ordering:** Groups appear alphabetically by default; override with `--group-order Diseased Trauma`.

---

### `run_qiime_marker_pipeline.py`

An alternative post-DADA2 entry point focused on the diversity analysis steps: control filtering, alpha rarefaction, core metrics, and optional PERMANOVA/PERMDISP stats with an auto-generated summary TSV.

Use this script if you already have a denoised feature table and want to run only the diversity steps, or if you prefer a simpler interface than `03_`'s subcommand system.

```bash
python run_qiime_marker_pipeline.py \
    --project-root /path/to/project \
    --marker       16S \
    --table        qiime2/dada2/table_16S.qza \
    --phylogeny    qiime2/diversity/rooted-tree_16S.qza \
    --metadata     metadata/qiime/metadata_16S.tsv \
    --controls-ids metadata/qiime/controls_16S_ids.txt \
    --sampling-depth 8000 \
    --max-depth    20000 \
    --do-beta-stats \
    --do-alpha-group-significance \
    --group-column Group
```

The `--do-beta-stats` flag runs both PERMANOVA and PERMDISP for all four distance metrics and automatically writes a `beta_stats_summary_<group>.tsv` file to the core metrics directory. This TSV can be passed directly to `07_visualize_diversity.py` via `--stats-tsv`.

---

## QIIME 2 Metadata Format

QIIME 2 metadata files must follow a precise format. Using the wrong format is one of the most common causes of import failures.

### Required format

```
sample-id<TAB>Group<TAB>Season<TAB>Date_Collected
#q2:types<TAB>categorical<TAB>categorical<TAB>categorical
TV230007-GI-16S_S1483<TAB>Diseased<TAB>Fall<TAB>9/30/22
TV230018-GI-16S_S1484<TAB>Diseased<TAB>Winter<TAB>1/20/23
TV230045-GI-16S_S1485<TAB>Trauma<TAB>Summer<TAB>7/14/22
```

### Rules

- **First column must be named exactly `sample-id`** (case-sensitive, hyphenated). Any other name will cause QIIME 2 to reject the file.
- **Tab-separated only.** CSV format will not import correctly.
- **The `#q2:types` row is optional** but strongly recommended. It tells QIIME 2 which columns are `categorical` (text, grouping variables) and which are `numeric`. Without it, QIIME 2 may misinterpret numeric-looking columns.
- **Sample IDs must exactly match** the IDs in your feature table — including suffixes like `_S1483` or `-GI-16S`. Use `02_make_qiime_metadata.py` to auto-generate a matching metadata file from an existing table.
- **Missing values should be left empty** — do not use `NA`, `nan`, `NULL`, or `N/A`. QIIME 2 interprets these differently across versions and they can cause silent errors in group-significance tests.
- **No leading or trailing whitespace** in column names or values.
- **Do not use special characters** (`/`, `\`, `#`, spaces) in the `sample-id` column.

### Generating metadata with `02_make_qiime_metadata.py`

If your source metadata uses short IDs (e.g., `TV230084`) but your feature table uses long IDs (e.g., `TV230084-GI-16S_S1492`), run:

```bash
python 02_make_qiime_metadata.py \
    --table              qiime2/dada2/table_16S.qza \
    --source-metadata    metadata/source_metadata.tsv \
    --source-id-column   TV \
    --out                metadata/qiime/metadata_16S.tsv \
    --marker             16S
```

The script extracts the `TV\d+` key from each table sample ID (configurable via `--key-regex`), looks it up in your source metadata, and writes a new TSV with `sample-id` values that exactly match the table.

---

## Choosing a Rarefaction Depth

Rarefaction subsamples all samples to the same sequencing depth so that diversity comparisons are not biased by unequal library sizes. Choosing the right depth requires balancing two competing goals: capturing true community diversity (favors higher depth) and retaining as many samples as possible (favors lower depth).

### Step 1 — Generate rarefaction curves

```bash
python 04_rarefaction.py \
    --table    qiime2/dada2/table_16S.qza \
    --tree     qiime2/diversity/rooted-tree_16S.qza \
    --metadata metadata/qiime/metadata_16S.tsv \
    --max-depth 20000 \
    --marker 16S \
    --outdir results/rarefaction/
```

Open the output `.qzv` at [view.qiime2.org](https://view.qiime2.org) to inspect individual sample curves interactively.

### Step 2 — Find the plateau

Look for the depth at which most samples' curves level off and become approximately horizontal. This is the "knee point" — where additional reads no longer discover new diversity. Your chosen depth should be at or beyond this plateau for the majority of samples.

**Curves still rising steeply** at a given depth means diversity estimates will be artificially low at that threshold. Consider deeper sequencing or acknowledge this as a limitation.

### Step 3 — Check sample retention

Any sample with fewer reads than the chosen depth is **permanently dropped** from all downstream diversity analyses. Evaluate the tradeoff:

| % Samples Retained | Interpretation |
|---|---|
| ≥ 90% | Generally acceptable |
| 80–90% | Acceptable; note which samples are dropped |
| 50–80% | Caution — check whether dropped samples are biased toward a particular group |
| < 50% | Reconsider; lower the depth or investigate low-yield samples |

A useful starting point is the **10th percentile** of sample read counts, which retains approximately 90% of samples. The rarefaction script reports this value automatically.

### Step 4 — Check for group-depth confounding

If samples from one group consistently have lower sequencing depth than samples from another, a high threshold will disproportionately drop samples from the lower-depth group. Inspect the per-group breakdown in the rarefaction QZV before finalizing a threshold.

### Step 5 — Evaluate your candidate depth

```bash
python 04_rarefaction.py \
    ... \
    --candidate-depth 8000
```

The script will report exactly which samples would be dropped at your candidate depth and flag whether retention is acceptable.

---

## Directory Structure

After running the pipeline, outputs are organized as:

```
project_root/
├── classifiers/
│   ├── silva-138-99-nb-classifier.qza
│   └── unite-ITS-classifier.qza
├── reads/
│   ├── 16S/
│   └── ITS/
├── metadata/
│   └── qiime/
│       ├── metadata_16S.tsv
│       └── metadata_ITS.tsv
├── qiime2/
│   └── {marker}/
│       └── {dataset}/
│           ├── imported/
│           │   ├── demux.qza
│           │   └── demux.qzv
│           ├── dada2/
│           │   ├── table.qza
│           │   ├── rep-seqs.qza
│           │   └── denoising-stats.qzv
│           ├── taxonomy/
│           │   └── taxonomy.qza
│           └── diversity/
│               └── core_metrics_depth{N}/
├── results/
│   └── {marker}/
│       └── {dataset}/
│           ├── tables/
│           ├── taxonomy/
│           ├── diversity/
│           ├── figures/
│           ├── differential/
│           └── exports/
├── logs/
│   └── {marker}/
│       └── {dataset}/
└── scripts/
```

---

## Color Palettes

Two palettes are available across `07_visualize_diversity.py`:

**`purple`** (default) — recommended for single-marker papers
- Group 1: dark purple `#7B2D8B`
- Group 2: lavender `#C19FD8`
- Significant stats: deep purple `#4A0060`
- PERMDISP: amber `#E59866`

**`redblue`** — recommended for multi-marker comparisons or colorblind accessibility
- Group 1: red `#B22222`
- Group 2: blue `#2E86C1`
- Significant stats: near-black `#1A1A1A`
- PERMDISP: amber `#E59866`

Both palettes support up to 8 groups with distinct colors and marker shapes.

---

## Notes and Known Limitations


**`04_run_diversity_pipeline.py` is retired.** This file was a partially-completed
evolution of `03_` that was never finished and contained a broken function reference.
All of its intended improvements (marker-aware filter defaults, `--include` positive
filtering, rarefaction depth preflight, non-phylogenetic COI/ITS diversity) have been
properly implemented in `03_run_full_metabarcoding_pipeline.py`.

**Tree building is not included.** Scripts `03_` and `run_qiime_marker_pipeline.py` both require a rooted phylogenetic tree QZA (`--tree`) for UniFrac metrics and Faith's PD. This tree must be built separately using QIIME 2's alignment and phylogeny pipeline:

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences  qiime2/{marker}/{dataset}/dada2/rep-seqs.qza \
    --o-alignment  qiime2/diversity/aligned-rep-seqs_{marker}.qza \
    --o-masked-alignment qiime2/diversity/masked-aligned-rep-seqs_{marker}.qza \
    --o-tree       qiime2/diversity/unrooted-tree_{marker}.qza \
    --o-rooted-tree qiime2/diversity/rooted-tree_{marker}.qza
```

**ANCOM vs. ANCOM-BC.** The `diff-ancom` subcommand in `03_` uses the legacy `qiime composition ancom`. For new analyses, ANCOM-BC is preferred: it accounts for unequal sampling fractions, is more robust with small sample sizes, and produces interpretable log fold-change estimates with FDR-corrected q-values. Run ANCOM-BC separately:

```bash
qiime composition ancombc \
    --i-table       results/{marker}/{dataset}/tables/table_filtered_L6.qza \
    --m-metadata-file metadata/qiime/metadata_{marker}.tsv \
    --p-formula     Group \
    --o-differentials qiime2/{marker}/{dataset}/ancombc_genus_{marker}_Group.qza
```

The output QZA can be parsed and plotted using the analysis scripts in this repository.

**`run_full_metabarcoding_pipeline.py` (unnumbered) is a superseded version** of `03_run_full_metabarcoding_pipeline.py`. It lacks marker-aware filter defaults, the rarefaction preflight check, and the `--include` positive filter argument. Use `03_` for all new analyses.


**COI diversity note:** COI amplicons are variable in length, making phylogenetic
alignment unreliable. `03_run_full_metabarcoding_pipeline.py` automatically
uses `core-metrics` (non-phylogenetic) when `--marker COI` is set, skipping
Faith's PD and UniFrac. Use Bray-Curtis and Jaccard for beta diversity only.
This behaviour can be overridden with `--no-phylo` (force non-phylo on any marker)
or by providing `--tree` manually.

**ITS rarefaction and phylogenetic metrics.** ITS amplicons are variable in length, making multiple sequence alignment unreliable. For ITS data, omit `--tree` / use `--no-tree` in `04_rarefaction.py`, and use non-phylogenetic metrics only (Bray-Curtis, Jaccard, Shannon, Observed Features) for diversity analyses.

**PERMDISP interpretation.** A significant PERMDISP result means groups differ in compositional *variance* (dispersion), not just in community composition (centroid position). When PERMDISP is significant and PERMANOVA is not — or vice versa — interpret results carefully. A significant PERMANOVA result accompanied by a significant PERMDISP should be noted in the text, as the PERMANOVA may be partially driven by unequal dispersion rather than true centroid separation. `05_visualize_diversity.py` displays PERMANOVA and PERMDISP stats separately per panel using distinct colors to make this distinction clear.
