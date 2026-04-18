# Wildlife Metabarcoding Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![QIIME2 2024.5](https://img.shields.io/badge/QIIME2-2024.5-green)](https://qiime2.org)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

A modular amplicon sequencing pipeline for dietary and gut community
metabarcoding in wildlife. Built around QIIME2 (2024.5) with Python analysis
and visualization scripts that are study-system agnostic.

**Developed by:** MEED Lab, University of New Hampshire  
**Contact:** Samantha Schade, Department of Natural Resources & the Environment  
**Reference study:** Common Loon (*Gavia immer*) gut microbiome and dietary
metabarcoding (Schade et al., in prep)

---

## Quick start

```bash
# 1. Install the conda environment
conda env create -f environment.yml
conda activate metabarcoding-pipeline

# 2. Copy the config template and edit for your project
cp pipeline_config.yml my_project_config.yml
# (edit: project root, markers, metadata paths, rarefaction depths, groups)

# 3. Validate your config
python pipeline.py --config my_project_config.yml check

# 4. Run
python pipeline.py --config my_project_config.yml run --steps all

# 5. Regenerate all figures
python pipeline.py --config my_project_config.yml figures
```

`pipeline_config.yml` is the only file you need to edit for a new project.
All scripts read their parameters from it — no hardcoded project-specific
values remain in any script.

---

## What this pipeline does

Starting from demultiplexed FASTQ reads, this pipeline produces:

- Species-level or genus-level taxonomy count tables per marker
- Presence/absence detection tables for dietary markers with ecological
  annotations (habitat, trophic role, lay-friendly common group labels)
- Alpha and beta diversity analyses with group significance testing
  (Kruskal-Wallis, PERMANOVA, PERMDISP — pairwise, BH-corrected)
- Publication-quality figures in SVG and PNG (300 dpi) with consistent
  Wong 2011 colorblind-safe palette throughout
- Pre-DADA2 QC report documenting primer detection, adapter onset, dimer
  rates, and quality drop-off per sample
- BLAST verification reports for artefact detection and novel taxon discovery
- Viral detection statistics (Fisher's exact + chi-square) from cutadapt
  primer-confirmed amplicons

---

## Markers supported

| Marker | Target | Reference database | Framework |
|---|---|---|---|
| 16S rRNA V4 | Bacteria (gut microbiome) | SILVA 138 (515F/806R) | Diversity (alpha/beta) |
| MiFish 12S | Fish prey | MitoFish + optional host spike-in | Presence/absence |
| Cytochrome b | Vertebrate prey | NCBI custom (RESCRIPt) | Presence/absence |
| 18S rRNA V9 | Eukaryotes / parasites | PR2 v5 | Presence/absence (exploratory) |
| ITS1-2 | Fungi | UNITE v10 | See note below |
| COI | Invertebrates | MIDORI2 UNIQ NUC | Presence/absence |
| Adenovirus | Avian adenovirus (hexon) | NCBI custom BLAST DB | Presence/absence |
| Herpesvirus | Herpesvirus DPOL | Cutadapt primer-confirmed | Presence/absence |

**Note on ITS1-2:** Fungal metabarcoding from gut contents of piscivorous
species typically yields insufficient read depth for diversity analysis. In
the loon study, mean DADA2 retention was 16.6% with NTC contamination
comparable to biological samples. ITS is appropriate for herbivores or
omnivores with dietary fungal components but not for fish-eating birds. If
you attempt ITS, run `00_build_classifiers.py --markers ITS` and check NTC
read counts via `05b_parse_dada2_retention.py` before proceeding.

---

## The key analytical decision: presence/absence vs relative abundance

For dietary markers (MiFish, cytb), both frameworks have merit and should be
reported when they disagree.

Clucas et al. (2024, bioRxiv 10.1101/2024.03.22.586275) validated that MiFish
12S relative read abundance correlates r=0.94 with prey biomass proportions in
wild piscivorous seabirds in the Gulf of Maine — the same ecosystem and primer
as the loon study. Deagle et al. (2019, Mol Ecol 28:391-406) showed through
simulation that frequency of occurrence systematically overestimates rare prey
items.

**Recommended approach:** Report rarefied relative abundance as the primary
analysis for MiFish (supported by Clucas et al. 2024) and binary
presence/absence as a secondary comparison. When the two frameworks disagree,
report both and interpret the discrepancy.

**For cytb:** Binary presence/absence is preferred at low rarefaction depths
(200 reads) where relative abundance estimates have variance too large to be
meaningfully interpreted.

**For microbiome markers (16S, 18S):** Use standard diversity analysis.

---

## Pipeline steps

Steps are run in order by `pipeline.py run --steps all`, or individually:

```
qc        Pre-DADA2 QC: primer detection + demux report
import    Merge sequencing runs + QIIME2 manifests + import FASTQs
denoise   Cutadapt primer trimming + DADA2 denoising
taxonomy  Taxonomic classification + export count tables
diversity Rarefaction curves + core-metrics + group significance tests
cod       Cause-of-death stratified diversity (secondary analysis)
figures   All publication figures (annotated and manuscript sets)
```

Each step is idempotent — re-running skips files that already exist.
Use `--force` to overwrite. Use `--dry-run` to preview commands.

---

## Three worked examples

### Example 1: Common Loon gut contents (multi-marker, aquatic piscivore)

**Study design:** n=40 loons (NH/ME 2022–2025); n=35 with gastrointestinal
tissue (16S, MiFish, cytb, 18S); n=40 with lung tissue (herpesvirus,
adenovirus). 5 loons had lung tissue only — the pipeline handles uneven
marker coverage across samples gracefully.

**Key decisions:**
- MiFish and cytb → presence/absence (Deagle et al. 2019)
- cytb rarefaction depth = 200 reads (constrained by sample depth)
- MiFish rarefaction depth = 17,000 reads; 16S = 8,000 reads
- Host spike-in classifier: loon 12S sequences added to MitoFish via
  `00_build_classifiers.py --markers MiFish --add-gavia`
- TV250064 excluded from MiFish (contamination confirmed: 226k reads
  dominated by Indo-Pacific taxa not present in New England loon diet)

**Typical output:** Jaccard PERMANOVA Diseased vs Trauma p=0.002 (MiFish);
Lead vs Trauma dietary composition p=0.002 surviving BH correction.

---

### Example 2: Tick blood meal identification (single marker, single host)

**Study design:** DNA from individual ticks (*Ixodes scapularis*), cytb
amplification to identify vertebrate host.

**Key decisions:**
- Single marker only — you are identifying one host, not a diet community
- Do NOT use the presence/absence framework — report the dominant signal
- Many samples will fail the minimum read threshold; plan for high dropout
- Exclude tick mitochondrial sequences: `--exclude Acari` in
  `07_taxonomy_table.py`
- Raise classifier confidence threshold: `--p-confidence 0.8`

```bash
python scripts/07_taxonomy_table.py \
    --taxonomy  qiime2/cytb/tick/taxonomy/taxonomy.qza \
    --table     qiime2/cytb/tick/dada2/table.qza \
    --marker    cytb \
    --include   Chordata \
    --exclude   Acari,Arachnida,Bacteria,Viruses,Archaea \
    --outdir    results/cytb/tick/taxonomy/
```

**Note:** Tick samples often contain human handler DNA. Add `Hominidae`
to your exclude list and log the removal.

---

### Example 3: Fisher scat (terrestrial predator, mixed prey)

**Study design:** Fisher (*Pekania pennanti*) scat, MiFish + cytb, generalist
carnivore diet.

**Key decisions:**
- cytb is more informative than MiFish for terrestrial prey — many prey
  items (squirrels, voles, rabbits, grouse) will not be amplified by MiFish
- Presence/absence framework (same rationale as loon dietary markers)
- Update `HOST_TAXON_STRINGS` in `11_clean_diet_table.py` for Mustelidae
- Cervidae and Leporidae are real prey — do NOT add to artefact list
- Update `PREY_LOOKUP` in `10b_annotate_diet_ecology.py` with terrestrial
  prey categories (the script header has a worked coyote scat example)

---

## Adapting for a new study system

### 1. Edit `pipeline_config.yml`

```yaml
project:
  name: "my_study"
  root: "."
  email: "your@institution.edu"

samples:
  id_regex: "(MYID\\d+)"   # regex to extract short ID from QIIME2 sample names
  control_prefixes: ["NTC-", "BLK-"]

active_markers: ["16S", "MiFish"]

markers:
  16S:
    classifier: "classifiers/silva-138-99-nb-classifier-515-806.qza"
    rarefaction_depth: 8000
    phylo: true
  MiFish:
    classifier: "classifiers/mitofish-12S-mifish-classifier.qza"
    rarefaction_depth: 5000
    phylo: false

groups:
  primary:
    column: "Group"
    order: ["Treatment", "Control"]
```

### 2. Adapt the host filter

Edit `HOST_TAXON_STRINGS` in `11_clean_diet_table.py`:

```python
# Loon default:
HOST_TAXON_STRINGS = {"MiFish": ["Gaviidae"], "cytb": ["Gaviidae", "Aves"]}

# Fisher scat:
HOST_TAXON_STRINGS = {"cytb": ["Mustelidae", "Pekania"]}

# Tick blood meal:
HOST_TAXON_STRINGS = {"cytb": ["Acari", "Ixodidae"]}
```

### 3. Adapt the ecological annotation lookup

Update `PREY_LOOKUP` in `10b_annotate_diet_ecology.py` with your prey
categories. The script header has step-by-step instructions and a worked
coyote scat example. Each entry maps to three values:
`(habitat_or_origin, trophic_role, common_group)`.

### 4. Adding a new marker

1. Add marker block to `pipeline_config.yml`
2. Train classifier: `python 00_build_classifiers.py --markers YOURMARKER`
3. Run pipeline: `python pipeline.py run --steps all`
4. Choose rarefaction depth from `06_rarefaction.py` output
5. For dietary markers: run `11_clean_diet_table.py` →
   `11b_presence_absence.py` → `10b_annotate_diet_ecology.py`
6. Optional BLAST verification: `07c_blast_qc_unclassified.py` for poorly
   classified ASVs; `08c_blast_verify.py` for suspect taxa by name

---

## Project directory structure

```
project_root/
├── pipeline_config.yml          ← edit this for your project
├── pipeline.py                  ← main entry point
├── environment.yml              ← conda environment
├── reads/                       ← demultiplexed FASTQs
├── classifiers/                 ← QIIME2 classifier QZAs
├── metadata/
│   ├── source_metadata.csv      ← your source-of-truth metadata
│   └── qiime/                   ← QIIME2-ready TSVs (built by pipeline)
├── qiime2/
│   └── {marker}/all/
│       ├── dada2/               ← table.qza, rep-seqs.qza
│       ├── taxonomy/            ← taxonomy.qza
│       └── diversity/           ← core-metrics output
├── results/
│   └── {marker}/
│       ├── all/taxonomy/        ← taxonomy count TSVs
│       ├── all/taxonomy_cleaned/← cleaned count TSVs
│       ├── all/presence_absence/← binary detection table
│       ├── all/taxonomy_annotated/ ← with common_group labels
│       ├── DvT/diversity/       ← PERMANOVA/PERMDISP QZVs
│       └── DvT/figures/         ← PCoA + alpha plots
├── manuscript_figures/          ← assembled by build_manuscript_figures.sh
├── reports/
│   ├── primers_detected.tsv
│   └── demultiplex/
├── logs/
└── scripts/
    ├── 00_build_classifiers.py
    ├── 04d_primer_advisor.py
    ├── 04e_parse_multiqc_demux.py
    ├── 05_run_full_metabarcoding_pipeline.py
    ├── 05b_parse_dada2_retention.py
    ├── suggest_dada2_params.py
    ├── 06_rarefaction.py
    ├── 07_taxonomy_table.py
    ├── 08_run_diversity_stats.py
    ├── 08b_run_cod_diversity.py
    ├── 08c_parse_beta_stats.py
    ├── 09_plot_diversity.py
    ├── 09_plot_taxonomy.py
    ├── 11_clean_diet_table.py
    ├── 11b_presence_absence.py
    ├── 11c_annotate_diet_ecology.py
    ├── 13_plot_habitat_season.py
    ├── plot_mifish_ecology_pa.py
    ├── 14_viral_stats.py
    ├── 14b_plot_adeno_tree.py
    ├── 14c_plot_herpes_cutadapt.py
    ├── build_manuscript_figures.sh
    ├── run_all_figures.py
    └── _archive/                ← retired scripts
```

---

## Key scripts reference

| Script | Purpose | Input | Output |
|---|---|---|---|
| `pipeline.py` | Main entry point | `pipeline_config.yml` | Orchestrates all steps |
| `00_build_classifiers.py` | Train classifiers + BLAST DBs | Reference FASTA | `classifier.qza` / BLAST DB |
| `04d_primer_advisor.py` | Pre-DADA2 QC | Raw FASTQs | `primers_detected.tsv` |
| `04e_parse_multiqc_demux.py` | Demux QC report | Illumina demux + MultiQC | `demux_qc_report.txt` |
| `04_make_qiime_metadata.py` | QIIME2 metadata builder | Source CSV + feature table | `metadata_{marker}.tsv` |
| `05_run_full_metabarcoding_pipeline.py` | Import → DADA2 → taxonomy | Reads + classifier | QIIME2 artifacts |
| `05b_parse_dada2_retention.py` | Post-DADA2 bottleneck detection | `denoising-stats.qza` | Retention report |
| `suggest_dada2_params.py` | DADA2 parameter advisor | demux QZV + primer/amplicon lengths | Ready-to-paste DADA2 command |
| `06_rarefaction.py` | Rarefaction guidance | Feature table | Curves + depth recommendation |
| `07_taxonomy_table.py` | Taxonomy count table | Taxonomy + table QZAs | Count + relabund TSVs |
| `07c_blast_qc_unclassified.py` | BLAST QC for unresolved ASVs | Rep-seqs + taxonomy | BLAST QC report |
| `08_run_diversity_stats.py` | PERMANOVA + Kruskal-Wallis | Core-metrics + metadata | Stats QZVs |
| `08b_run_cod_diversity.py` | COD-stratified diversity | Core-metrics + COD metadata | Stats QZVs |
| `08c_parse_beta_stats.py` | Parse QZV stats to TSV | Results directory | Summary TSV with sig stars |
| `08c_blast_verify.py` | Verify suspect taxa by name/reads | Rep-seqs + taxonomy | AGREE/DISAGREE report |
| `09_plot_diversity.py` | PCoA + alpha figures | Core-metrics QZAs | PNG + SVG |
| `09_plot_taxonomy.py` | Taxonomy barplots | Relabund TSV + metadata | PNG + SVG |
| `11_clean_diet_table.py` | Remove host/artefacts + log | Count TSV | Cleaned TSV + cleaning report |
| `11b_presence_absence.py` | Detection frequency analysis | Cleaned TSV + metadata | Binary table + barplot |
| `11c_annotate_diet_ecology.py` | Ecological annotations | Cleaned TSV | Annotated TSV |
| `13_plot_habitat_season.py` | Habitat by season figure | Annotated TSV + metadata | PNG + SVG |
| `plot_mifish_ecology_pa.py` | Detection frequency by season and COD | PA table + annotation + metadata | PNG + SVG |
| `14_viral_stats.py` | Viral Fisher's exact + chi-square | Excel workbook + metadata | Stats TSV |
| `14b_plot_adeno_tree.py` | Adenovirus phylogenetic tree | IQ-TREE Newick | PNG + SVG |
| `14c_plot_herpes_cutadapt.py` | Herpesvirus detection figures | Cutadapt summary + metadata | PNG + SVG |
| `run_all_figures.py` | Regenerate all figures | Config | Full figure sets |
| `build_manuscript_figures.sh` | Assemble numbered figure directory | Results | `manuscript_figures/` |

---

## Citation

If you use this pipeline, please cite:

> Schade SN et al. (in prep). Lead toxicosis disrupts foraging ecology and gut
> microbiome diversity in Common Loons (*Gavia immer*): a multi-marker
> metabarcoding study. MEED Lab, University of New Hampshire.

And cite the pipeline itself:

> Schade SN (2026). Wildlife Metabarcoding Pipeline v1.0.0. MEED Lab,
> University of New Hampshire. https://doi.org/10.5281/zenodo.XXXXXXX

Please also cite the key frameworks this pipeline depends on:

> Deagle BE, Thomas AC, McInnes JC, et al. (2019). Counting with DNA in
> metabarcoding studies: How should we convert sequence reads to dietary data?
> *Molecular Ecology*, 28(2), 391–406.

> Callahan BJ, McMurdie PJ, Rosen MJ, et al. (2016). DADA2: High-resolution
> sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583.

> Bolyen E, Rideout JR, Dillon MR, et al. (2019). Reproducible, interactive,
> scalable and extensible microbiome data science using QIIME 2.
> *Nature Biotechnology*, 37, 852–857.

> Wong B (2011). Points of view: Color blindness. *Nature Methods*, 8, 441.

---

## License

MIT License — see [LICENSE](LICENSE).

Reference databases (SILVA, UNITE, MitoFish, PR2, MIDORI2, NCBI) are subject
to their own licensing terms; see LICENSE for details.

---

## Contact

Open a GitHub issue or contact the MEED Lab at the University of New Hampshire
for questions about adapting this pipeline to new study systems.
