# Wildlife Metabarcoding Pipeline

A modular amplicon sequencing pipeline for dietary and gut community
metabarcoding in wildlife. Built for QIIME2 (2024.5) with Python analysis
scripts that are study-system agnostic.

**Developed by:** MEED Lab, University of New Hampshire  
**Contact:** Samantha Schade, Department of Natural Resources & the Environment  
**Reference study:** Common Loon (*Gavia immer*) gut microbiome and dietary
metabarcoding (Schade et al., in prep)

---

## What this pipeline does

Starting from demultiplexed FASTQ reads, this pipeline produces:

- Species-level or genus-level taxonomy count tables per marker
- Presence/absence detection tables for dietary markers with ecological
  annotations (habitat, trophic role, lay-friendly common group labels)
- Alpha and beta diversity analyses with group significance testing
- Publication-quality figures in SVG and PNG (300 dpi)
- Pre-DADA2 QC report documenting primer detection, adapter onset, and
  dimer rates per sample

**Markers supported out of the box:**

| Marker | Target | Reference database | Framework |
|---|---|---|---|
| 16S rRNA V4 | Bacteria (gut microbiome) | SILVA 138 (515F/806R) | Diversity (alpha/beta) |
| MiFish 12S | Fish prey | MitoFish + host spike-in | Presence/absence |
| Cytochrome b | Vertebrate prey | NCBI cytb (custom) | Presence/absence |
| 18S rRNA | Eukaryotes / parasites | PR2 v5 | Presence/absence (exploratory) |
| ITS1-2 | Fungi | UNITE v10 | *Not recommended for gut contents* — see note below |

**Note on ITS1-2:** Fungal metabarcoding from gut contents of piscivorous
species typically yields insufficient read depth for diversity analysis.
In the loon study, mean DADA2 retention was 16.6% with NTC contamination
comparable to biological samples. ITS is appropriate for herbivores or
omnivores with dietary fungal components, but not for fish-eating birds.
If you attempt ITS, run `parse_dada2_retention.py` before proceeding and
check NTC read counts carefully.

Adding a new marker requires one configuration block in `05_run_full_metabarcoding_pipeline.py`
and a reference database in QIIME2 format. See [Adding a new marker](#adding-a-new-marker).

---

## The key analytical decision: presence/absence vs relative abundance

**For dietary markers (MiFish, cytb), both frameworks have merit and should
be reported when they disagree.**

Clucas et al. (2024, bioRxiv doi:10.1101/2024.03.22.586275) validated that
MiFish 12S relative read abundance (RRA) correlates r=0.94 with prey biomass
proportions in wild piscivorous seabirds (Gulf of Maine common terns) — the
same ecosystem and primer as the loon study. Deagle et al. (2019, Mol Ecol
28:391-406) showed through simulation that frequency of occurrence (FOO/
presence-absence) systematically overestimates rare prey items.

**Recommended approach:** Report rarefied RRA as the primary analysis for
MiFish (supported by Clucas et al. 2024 validation) and binary presence-
absence as a secondary comparison. When the two frameworks disagree, report
both and interpret the discrepancy.

**For cytb:** Binary presence-absence is preferred at low rarefaction depths
(200 reads) where RRA estimates are noisy. PA is more defensible when read
depth does not support reliable proportional estimates.

See: Deagle et al. (2019) *Molecular Ecology* — the standard treatment of this
issue for dietary metabarcoding studies.

**For microbiome markers (16S, 18S), use standard diversity analysis.** The
relative abundance framework is appropriate because you are characterizing a
community, not estimating prey consumption.

---

## Three worked examples

### Example 1: Common Loon gut contents (multi-marker, aquatic)

**Study design:** n=35 loons, four markers (16S, MiFish, cytb, 18S), 
gut content from GI tract dissection.

**Key decisions:**
- MiFish and cytb → presence/absence (Deagle et al. 2019 rationale)
- cytb rarefaction depth = 200 reads (constrained by sample depth)
- MiFish rarefaction depth = 17,000 reads
- 16S rarefaction depth = 8,000 reads
- Host spike-in classifier (loon 12S sequences added to MitoFish) to
  distinguish loon host reads from prey
- Ecological annotation using `10b_annotate_diet_ecology.py` with prey
  categorised as marine forage fish, anadromous fish, freshwater prey fish,
  marine bottom fish, estuarine fish, and sand lance
- TV250064 excluded from MiFish analysis (contamination confirmed: 226k reads
  dominated by Indo-Pacific taxa not expected in New England loon diet)

**Reference species/taxa:** Gaviidae excluded as host from MiFish;
broader Aves excluded from cytb (Greengenes k__/p__ prefix format requires
`--include Chordata --exclude Aves`).

**Typical output:** Jaccard PERMANOVA Diseased vs Trauma p=0.002 (MiFish),
Lead vs Trauma dietary composition p=0.002 surviving FDR correction.

---

### Example 2: Tick blood meal identification (single-marker, single-host)

**Study design:** DNA extracted from individual ticks (*Ixodes scapularis*),
single 12S or cytb amplification to identify vertebrate host.

**Key decisions:**
- Single marker only — you are identifying one host, not characterising a diet
- Do NOT use presence/absence framework — you want the dominant signal, which
  is the host. Report the top BLAST hit or highest-confidence classifier call.
- Rarefaction depth is very low (often 200-500 reads per tick)
- Many samples will fail the minimum read threshold — plan for high dropout rate
- Exclude tick mitochondrial sequences as host reads (add tick sequences to
  classifier or use `--exclude Acari` in `07_taxonomy_table.py`)
- Classifier confidence threshold matters more here than in dietary studies —
  use `--p-confidence 0.8` or higher in `classify-sklearn`

**Configuration changes from loon default:**
```bash
# 07_taxonomy_table.py for tick blood meal (cytb):
python scripts/07_taxonomy_table.py \
    --taxonomy  qiime2/cytb/tick/taxonomy/taxonomy.qza \
    --table     qiime2/cytb/tick/dada2/table.qza \
    --marker    cytb \
    --include   Chordata \
    --exclude   Acari,Arachnida,Bacteria,Viruses,Archaea \
    --outdir    results/cytb/tick/taxonomy/

# No 11_clean_diet_table.py or 11b_presence_absence.py needed
# Report top classification per sample directly from taxonomy_counts TSV
```

**Note on host DNA contamination:** Tick samples often contain human handler
DNA. Add `Hominidae` to your exclude list and log the removal. If processing
in the field, also add `Bovidae` and `Canidae` for livestock/dog handler
contamination.

---

### Example 3: Fisher scat and stomach contents (terrestrial predator, mixed prey)

**Study design:** Fisher (*Pekania pennanti*) scat and stomach contents,
MiFish 12S + cytb, aiming to characterise diet of a generalist carnivore.

**Key decisions:**
- Presence/absence framework (same rationale as loon dietary markers)
- cytb is more informative than MiFish for terrestrial prey — many prey items
  (squirrels, voles, rabbits, grouse) are not fish and MiFish will not amplify
  them. Use cytb as primary dietary marker, MiFish as secondary.
- Fisher host DNA must be excluded: add Mustelidae to host filter in
  `11_clean_diet_table.py`. For cytb this requires adding `Mustelidae` to
  `HOST_TAXON_STRINGS` in the script (currently set for Gaviidae/Aves).
- Artefact taxa for terrestrial study: Canidae and Bovidae are still domestic
  contamination. Cervidae and Leporidae are real prey — do NOT remove them.
- Scent marking behaviour means scats may contain prey DNA from multiple meals.
  Presence/absence per scat is the correct unit.
- Environmental fungi (ITS) may be informative for habitat use if
  berries/mushrooms are consumed seasonally.

**Adapting the ecological annotation lookup table:**

The `10b_annotate_diet_ecology.py` script has a built-in lookup table for
loon fish prey. For a fisher study, export the template and replace it:

```bash
# Export the built-in loon lookup as a starting template
python scripts/10b_annotate_diet_ecology.py \
    --export-lookup config/fisher_prey_lookup_template.tsv

# Edit config/fisher_prey_lookup_template.tsv with your prey categories:
# taxon | habitat | trophic_role | common_group
# Sylvilagus floridanus | Terrestrial | Small mammal | Rabbits and hares
# Microtus pennsylvanicus | Terrestrial | Small mammal | Rodents
# Meleagris gallopavo | Terrestrial | Ground bird | Upland birds
# Vaccinium | Plant | Fruit/berry | Berries and fruit

# Run with custom lookup
python scripts/10b_annotate_diet_ecology.py \
    --counts       results/cytb/all/taxonomy_cleaned/taxonomy_counts_cleaned_cytb.tsv \
    --lookup       config/fisher_prey_lookup_template.tsv \
    --study-system "Fisher scat cytb" \
    --habitat-label prey_origin \
    --outdir       results/cytb/all/taxonomy_annotated/
```

Note the `--habitat-label prey_origin` flag — for a terrestrial study "habitat"
is misleading. The column will be named `prey_origin` in the output instead.

---

## Installation

### Requirements

- Conda (Miniconda or Anaconda)
- QIIME2 amplicon distribution 2024.5
- ~50 GB disk space for classifiers

### QIIME2 environment

```bash
conda env create -n qiime2-amplicon-2024.5 \
    --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.5-py39-linux-conda.yml
conda activate qiime2-amplicon-2024.5
```

### Analysis environment

```bash
conda create -n metabarcoding-analysis python=3.10
conda activate metabarcoding-analysis
pip install pandas numpy matplotlib scipy biopython openpyxl
```

### Reference databases

Download pre-built classifiers into `classifiers/`:

```bash
# SILVA 138 V4 (16S) — QIIME2 hosted
wget -O classifiers/silva-138-99-nb-classifier-515-806.qza \
  "https://data.qiime2.org/2024.5/common/silva-138-99-nb-classifier-515-806.qza"

# UNITE v10 (ITS) — hosted by Colin Brislawn
wget -O classifiers/unite-ver10-99-nb-classifier.qza \
  "https://github.com/colinbrislawn/unite-train/releases/download/v10.0-v04.04.2024-qiime2-2024.5/unite_ver10_dynamic_all_04.04.2024-Q2-2024.5.qza"

# MiFish (12S fish) — build from MitoFish database
# See 01_build_classifiers.py for instructions
# Optionally add host sequences to exclude host reads at classification stage
```

---

## Input data requirements

### Metadata file (required)

A CSV or TSV with one row per sample. Required columns:

| Column | Description | Example |
|---|---|---|
| sample ID | Unique identifier matching FASTQ filenames | TV230084 |
| group | Primary biological grouping | Diseased, Trauma, Control |
| collection date | For seasonal analysis (M/D/YYYY) | 8/15/2023 |

**Seasonal analysis note:** Use ecologically meaningful seasons for your
study species rather than meteorological calendar seasons. For loons:
Breeding (May-Aug), Freshwater_Nonbreeding (Apr+Sep), Saltwater (Oct-Mar).
Edit `month_to_season()` in `02b_add_season_to_metadata.py` for your species.
Always verify season is independent of your primary grouping variable using
chi-square before reporting seasonal results.
| collection site | For confound checking | Lake Umbagog, NH |

The script `04_make_qiime_metadata.py` builds a QIIME2-ready TSV from your
source metadata by matching sample IDs between your metadata file and the
QIIME2 feature table.

### FASTQ reads

Demultiplexed, paired-end reads in separate directories per marker:
```
reads/
├── 16S/
│   ├── SAMPLE1-16S_L001_R1_001.fastq.gz
│   └── SAMPLE1-16S_L001_R2_001.fastq.gz
├── MiFish/
│   ├── SAMPLE1-MiFish_L001_R1_001.fastq.gz
│   └── ...
```

If reads come from multiple sequencing runs, use `02_merge_run_dirs.py` first.

---

## Configuring for your study system

### Four things you must change

**1. Marker names and primer sequences** — in `05_run_full_metabarcoding_pipeline.py`:
```python
MARKERS = {
    "MiFish": {
        "primer_f": "GTCGGTAAAACTCGTGCCAGC",
        "primer_r": "CATAGTGGGGTATCTAATCCCAGTTTG",
        "trunc_len_f": 220,
        "trunc_len_r": 180,
        "classifier": "classifiers/mitofish-12S-mifish-gavia-classifier.qza",
    },
    # Add your marker here
}
```

**2. Rarefaction depth** — in `06_rarefaction.py` and diversity commands.
Choose a depth that retains ≥70% of your samples. Check the rarefaction
curve before committing. For low-biomass samples (tick, museum specimens),
acceptable depths may be 200-500 reads.

**3. Taxonomy filter strings** — in `07_taxonomy_table.py`.
The `--include` and `--exclude` flags depend on your reference database format.
SILVA uses plain NCBI-style strings ("Bacteria", "Eukaryota").
UNITE uses plain strings ("Fungi"). Custom cytb databases built from NCBI
may use Greengenes-style k__/p__ prefixes — check with:
```bash
unzip -p your_taxonomy.qza "*/data/taxonomy.tsv" | head -5
```

**4. Detection thresholds** — in `11b_presence_absence.py`.
```
--min-sample-reads: minimum total reads to include a sample (default: 10000)
--min-taxon-reads:  minimum reads for a taxon to count as detected (default: 10)
--min-relabund:     minimum relative abundance threshold (default: 0.01 = 1%)
```
Lower thresholds for sparse markers (cytb: use 50/5/0.01).
Higher thresholds for noisy markers or low-quality extractions.

### Adapting the host filter

The `11_clean_diet_table.py` script filters host reads using family-level
strings. The default is set for loon (Gaviidae). Change `HOST_TAXON_STRINGS`
in the script for your study organism:

```python
HOST_TAXON_STRINGS = {
    "MiFish": ["Gaviidae", "Gaviiformes"],    # Loon default
    "cytb":   ["Gaviidae", "Gaviiformes", "Aves"],
}

# For fisher scat (cytb):
HOST_TAXON_STRINGS = {
    "cytb": ["Mustelidae", "Pekania"],
}

# For tick blood meal (cytb):
HOST_TAXON_STRINGS = {
    "cytb": ["Acari", "Ixodidae", "Ixodes"],
}
```

---

## Adding a new marker

1. Add primer sequences and DADA2 parameters to `05_run_full_metabarcoding_pipeline.py`
2. Add the marker's classifier path to `01_build_classifiers.py` or download pre-built
3. Create a metadata TSV for the marker using `04_make_qiime_metadata.py`
4. Run the full pipeline via `05_run_full_metabarcoding_pipeline.py --marker YOURMARKER`
5. Choose rarefaction depth from `06_rarefaction.py` output
6. Run diversity stats via `08_run_diversity_stats.py`
7. Run taxonomy table via `07_taxonomy_table.py` with appropriate `--include/--exclude`
7b. **Optional BLAST verification** — run `08c_blast_verify.py` on suspect taxa
   (ecologically implausible, low classifier confidence, or high read counts
   with unresolved taxonomy). Requires local NCBI nt database or NCBI remote
   access. Updates artefact exclusion list in `11_clean_diet_table.py`.
8. If dietary marker: run `11_clean_diet_table.py` → `11b_presence_absence.py` →
   `10b_annotate_diet_ecology.py`
9. Add the marker to `run_all_figures.sh` for the diversity figure generation block

---

## Project directory structure

```
project_root/
├── reads/                        ← demultiplexed FASTQs, one subdir per marker
├── classifiers/                  ← QIIME2 classifier QZAs
├── metadata/
│   ├── full_metadata_{study}.csv ← your source-of-truth metadata
│   └── qiime/                   ← QIIME2-ready TSVs (built by 04_make_qiime_metadata.py)
├── qiime2/                       ← QIIME2 artifacts (.qza/.qzv)
│   └── {marker}/all/
│       ├── dada2/               ← table.qza, rep-seqs.qza
│       ├── taxonomy/            ← taxonomy.qza
│       └── diversity/           ← core-metrics output
├── results/
│   └── {marker}/
│       ├── all/taxonomy/        ← taxonomy count TSVs
│       ├── all/taxonomy_cleaned/← cleaned count TSVs (after 09b)
│       ├── all/presence_absence/← detection freq + binary table (after 08b)
│       ├── all/taxonomy_annotated/ ← with common group labels (after 10b)
│       ├── DvT/diversity/       ← PERMANOVA/PERMDISP QZVs
│       └── DvT/figures/         ← PCoA + alpha plots
├── reports/
│   ├── primers_detected.tsv     ← from primer_advisor.py detect
│   └── demultiplex/             ← Illumina demux + MultiQC output
├── logs/                        ← one subdir per marker
├── scripts/
│   ├── PIPELINE.md              ← study-specific execution guide
│   ├── utils/                   ← primer_advisor.py, parse_multiqc_demux.py
│   └── _archive/                ← retired scripts
└── envs/                        ← conda environment YAML files
```

---

## Key scripts

| Script | Purpose | Input | Output |
|---|---|---|---|
| `utils/primer_advisor.py` | Pre-DADA2 QC | Raw FASTQs | primers_detected.tsv |
| `utils/parse_multiqc_demux.py` | Demux QC report | Illumina demux + MultiQC | demux_qc_report.txt |
| `01_build_classifiers.py` | Build/train classifiers | Reference FASTA + taxonomy | classifier.qza |
| `03_make_manifests.py` | FASTQ manifests | reads/ directory | manifest_{marker}.tsv |
| `04_make_qiime_metadata.py` | QIIME2 metadata | Source CSV + feature table | metadata_{marker}.tsv |
| `02b_add_season_to_metadata.py` | Add Season column | Metadata TSV | Updated metadata TSV |
| `05_run_full_metabarcoding_pipeline.py` | Import → DADA2 → taxonomy | Reads + classifier | QIIME2 artifacts |
| `06_rarefaction.py` | Choose rarefaction depth | Feature table | Rarefaction curves |
| `08_run_diversity_stats.py` | PERMANOVA + alpha stats | Core-metrics + metadata | Stats QZVs |
| `05b_run_cod_diversity.py` | COD-filtered diversity | Core-metrics + COD metadata | Stats QZVs |
| `05c_parse_beta_stats.py` | Parse QZV stats | Results directory | Summary TSV |
| `07_taxonomy_table.py` | Taxonomy count table | Taxonomy + table QZAs | Count TSV |
| `11_clean_diet_table.py` | Remove host/artefacts | Count TSV | Cleaned TSV |
| `11b_presence_absence.py` | Detection analysis | Cleaned TSV + metadata | Binary table + barplot |
| `10b_annotate_diet_ecology.py` | Ecological annotations | Cleaned TSV | Annotated TSV |
| `09_plot_diversity.py` | PCoA + alpha figures | Core-metrics QZAs | PNG + SVG figures |
| `10_plot_taxonomy.py` | Taxonomy barplots | Relabund TSV | PNG + SVG figures |
| `11_plot_mifish_season_ecology.py` | Seasonal diet figure | PA table + metadata | PNG + SVG |
| `10_plot_viral.py` | Viral detection figures | Excel workbook | PNG + SVG figures |
| `14_viral_stats.py` | Viral Fisher's exact | Excel + metadata | Stats TSV |
| `run_all_figures.sh` | All diversity figures | Core-metrics QZAs | All PNG + SVG figures |

---

## Citation

If you use this pipeline, please cite:

> Schade SN et al. (in prep). Lead toxicosis disrupts foraging ecology and gut
> microbiome diversity in Common Loons (*Gavia immer*): a multi-marker
> metabarcoding study. MEED Lab, University of New Hampshire.

And the key analytical frameworks:

> Deagle BE, Thomas AC, McInnes JC, et al. (2019). Counting with DNA in
> metabarcoding studies: How should we convert sequence reads to dietary data?
> *Molecular Ecology*, 28(2), 391–406.

> Callahan BJ, McMurdie PJ, Rosen MJ, et al. (2016). DADA2: High-resolution
> sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583.

> Bolyen E, Rideout JR, Dillon MR, et al. (2019). Reproducible, interactive,
> scalable and extensible microbiome data science using QIIME 2.
> *Nature Biotechnology*, 37, 852–857.

---

## License

MIT License — see LICENSE.txt. Reference databases (SILVA, UNITE, MitoFish,
NCBI) are subject to their own licensing terms.

---

## Contact

Open an issue or contact the MEED Lab at UNH for questions about adapting
this pipeline to new study systems.
