#!/usr/bin/env bash
# =============================================================================
#  rename_pipeline_scripts.sh
#  MEED Lab — Wildlife Metabarcoding Pipeline
#
#  Renumbers all pipeline scripts to reflect scientifically correct
#  execution order: BLAST QC and feature table filtering must happen
#  BEFORE diversity analysis, which must happen BEFORE figures.
#
#  Run from the loon_project root:
#    bash scripts/rename_pipeline_scripts.sh
#
#  This script:
#    1. Renames scripts in scripts/ directory
#    2. Removes confirmed duplicate scripts
#    3. Prints a summary of what changed
#    4. Does NOT modify pipeline.py or docstrings — see MANUAL STEPS below
#
#  IMPORTANT: Run on a clean git branch. Review the diff before committing.
# =============================================================================

set -euo pipefail
SCRIPTS="$(dirname "$0")"
cd "$SCRIPTS"

echo "=============================================="
echo "  PIPELINE SCRIPT RENUMBERING"
echo "  Working directory: $(pwd)"
echo "=============================================="
echo ""

# ── Helper ─────────────────────────────────────────────────────────────────
rename_script() {
    local old="$1"
    local new="$2"
    if [[ -f "$old" ]]; then
        mv "$old" "$new"
        echo "  ✅  $old  →  $new"
    else
        echo "  ⚠   NOT FOUND (skipped): $old"
    fi
}

remove_duplicate() {
    local file="$1"
    local reason="$2"
    if [[ -f "$file" ]]; then
        rm "$file"
        echo "  🗑   REMOVED duplicate: $file  ($reason)"
    else
        echo "  —   Already gone: $file"
    fi
}

# =============================================================================
# STAGE 1 — SETUP
# Scripts that prepare the environment, build classifiers, merge raw reads,
# and create metadata before any QIIME2 processing begins.
# =============================================================================
echo "── STAGE 1: Setup ────────────────────────────────────────────────────"
rename_script "00_build_classifiers.py"        "01_build_classifiers.py"
rename_script "00_merge_run_dirs.py"            "02_merge_run_dirs.py"
rename_script "01_make_manifests.py"            "03_make_manifests.py"
rename_script "02_make_qiime_metadata.py"       "04_make_qiime_metadata.py"
rename_script "02b_add_season_to_metadata.py"   "04b_add_season_to_metadata.py"

# add_season_to_metadata.py (no number) uses meteorological seasons (Winter/Spring/Summer/Fall)
# 04b_ uses ecological loon seasons (Breeding/Freshwater_Nonbreeding/Saltwater)
# Both are kept — they serve different purposes.
# Rename the generic meteorological version clearly:
rename_script "add_season_to_metadata.py"       "04c_add_meteorological_season_to_metadata.py"

# Primer detection, DADA2 parameter suggestion, and post-cutadapt QC.
# Three subcommands run in order: detect → suggest → check
# This informs truncation parameters before denoising begins.
rename_script "primer_advisor.py"              "04d_primer_advisor.py"

# Pre-DADA2 QC report from Illumina demultiplex stats and MultiQC FastQC data.
# Run after sequencing, before QIIME2 import, to verify read counts and quality.
rename_script "parse_multiqc_demux.py"         "04e_parse_multiqc_demux.py"

echo ""

# =============================================================================
# STAGE 2 — DENOISING
# Runs cutadapt trimming and DADA2 denoising via QIIME2.
# NOTE: The 'diversity' subcommand of this script should NOT be run here.
# Diversity must be run AFTER BLAST QC (Stage 4). See pipeline.py step_diversity.
# =============================================================================
echo "── STAGE 2: Denoising ────────────────────────────────────────────────"
rename_script "03_run_full_metabarcoding_pipeline.py" "05_run_full_metabarcoding_pipeline.py"

# Diagnostic for DADA2 read retention — runs immediately after DADA2 to verify
# how many reads survived each denoising step. Flags samples with poor retention
# and identifies which step is losing reads before any downstream analysis.
rename_script "parse_dada2_retention.py"       "05b_parse_dada2_retention.py"
echo ""

# =============================================================================
# STAGE 3 — RAREFACTION DECISION
# Generate rarefaction curves to choose sampling depth before core-metrics.
# Must run after denoising, before diversity.
# =============================================================================
echo "── STAGE 3: Rarefaction ──────────────────────────────────────────────"
rename_script "04_rarefaction.py"  "06_rarefaction.py"
echo ""

# =============================================================================
# STAGE 4 — TAXONOMY EXPORT AND BLAST QC
# Export taxonomy count tables, then BLAST-verify assignments before any
# downstream analysis. Filtering the feature table here ensures diversity
# metrics and taxonomy plots are both computed on clean, verified data.
#
# Order within stage:
#   07_  Export QIIME2 taxonomy to human-readable count TSVs
#   07b_ BLAST-verify named/suspect taxa (targeted)
#   07c_ BLAST QC all poorly-classified ASVs (broad safety net)
#   07d_ Filter confirmed artefact ASVs from QIIME2 feature table
# =============================================================================
echo "── STAGE 4: Taxonomy export and BLAST QC ─────────────────────────────"
rename_script "08_taxonomy_table.py"            "07_taxonomy_table.py"
rename_script "08c_blast_verify.py"             "07b_blast_verify.py"
rename_script "09a_blast_qc_unclassified.py"    "07c_blast_qc_unclassified.py"
# 07d_filter_feature_table.py — create this new script (see NOTE below)
echo "  📝  NOTE: 07d_filter_feature_table.py needs to be created."
echo "      This wraps: qiime feature-table filter-features"
echo "      Input: confirmed_artefacts_{marker}.txt from 07c_"
echo "      Output: filtered table.qza ready for diversity"
echo ""

# =============================================================================
# STAGE 5 — DIVERSITY ANALYSIS
# Runs AFTER the feature table has been filtered in Stage 4.
# Core-metrics is called via 05_run_full_metabarcoding_pipeline.py diversity
# subcommand (or directly via qiime), then stats and COD analysis follow.
# =============================================================================
echo "── STAGE 5: Diversity analysis ───────────────────────────────────────"
rename_script "05_run_diversity_stats.py"       "08_run_diversity_stats.py"
rename_script "05b_run_cod_diversity.py"        "08b_run_cod_diversity.py"
rename_script "05c_parse_beta_stats.py"         "08c_parse_beta_stats.py"
echo ""

# =============================================================================
# STAGE 6 — DIVERSITY VISUALIZATION
# Generates publication-quality PCoA and alpha diversity figures.
# Runs after diversity stats are complete.
# =============================================================================
echo "── STAGE 6: Diversity visualization ─────────────────────────────────"
rename_script "06_plot_diversity.py"            "09_plot_diversity.py"
rename_script "06b_combine_multimarker_alpha.py" "09b_combine_multimarker_alpha.py"
rename_script "07_visualize_diversity.py"       "09c_visualize_diversity.py"
echo ""

# =============================================================================
# STAGE 7 — TAXONOMY PLOTS
# Generates stacked barplots from the taxonomy count tables.
# Runs after BLAST QC so labels reflect verified assignments.
# =============================================================================
echo "── STAGE 7: Taxonomy plots ───────────────────────────────────────────"
rename_script "09_plot_taxonomy.py"             "10_plot_taxonomy.py"
echo ""

# =============================================================================
# STAGE 8 — DIET ANALYSIS (dietary markers: MiFish, cytb)
# Clean diet tables, annotate ecology, run presence/absence.
# Must run after taxonomy plots (Stage 7) and before ecological figures (Stage 9).
# Order within stage:
#   11_  Clean diet table (remove host, artefacts, collapse to species)
#   11b_ Presence/absence detection frequency analysis
#   11c_ Annotate taxa with habitat ecology categories
# =============================================================================
echo "── STAGE 8: Diet analysis ────────────────────────────────────────────"
rename_script "09b_clean_diet_table.py"         "11_clean_diet_table.py"
rename_script "08b_presence_absence.py"         "11b_presence_absence.py"
rename_script "10b_annotate_diet_ecology.py"    "11c_annotate_diet_ecology.py"
echo ""

# =============================================================================
# STAGE 9 — ECOLOGICAL FIGURES
# MiFish seasonal ecology and habitat composition figures.
# Require annotated diet tables from Stage 8.
# =============================================================================
echo "── STAGE 9: Ecological figures ───────────────────────────────────────"

# 09b_plot_mifish_season_ecology.py and 11b_plot_mifish_season_ecology.py are
# confirmed identical duplicates of 11_plot_mifish_season_ecology.py.
# Remove duplicates, keep and rename the canonical version.
remove_duplicate "09b_plot_mifish_season_ecology.py" "identical to 11_plot_mifish_season_ecology.py"
remove_duplicate "11b_plot_mifish_season_ecology.py" "identical to 11_plot_mifish_season_ecology.py"
rename_script    "11_plot_mifish_season_ecology.py"  "12_plot_mifish_season_ecology.py"
rename_script    "12_plot_habitat_season.py"         "13_plot_habitat_season.py"
echo ""

# =============================================================================
# STAGE 10 — VIRAL DETECTION AND STATS
# Adenovirus and herpesvirus figures and statistics.
# Independent of dietary analysis — can run in parallel with Stages 7-9.
# =============================================================================
echo "── STAGE 10: Viral ───────────────────────────────────────────────────"
rename_script "10_plot_viral.py"                "14_plot_viral.py"
rename_script "10b_plot_adeno_tree.py"          "14b_plot_adeno_tree.py"
rename_script "10b_plot_herpes_cutadapt.py"     "14c_plot_herpes_cutadapt.py"
rename_script "10c_plot_adeno_cutadapt.py"      "14d_plot_adeno_cutadapt.py"
rename_script "14_viral_stats.py"               "15_viral_stats.py"
echo ""

# =============================================================================
# SUMMARY
# =============================================================================
echo "=============================================="
echo "  RENAMING COMPLETE"
echo "=============================================="
echo ""
echo "MANUAL STEPS REQUIRED AFTER THIS SCRIPT:"
echo ""
echo "  1. UPDATE pipeline.py"
echo "     The following script name strings need updating:"
echo "       '01_make_manifests.py'             → '03_make_manifests.py'"
echo "       '03_run_full_metabarcoding_pipeline.py' → '05_run_full_metabarcoding_pipeline.py'"
echo "       '05_run_diversity_stats.py'         → '08_run_diversity_stats.py'"
echo "       '05b_run_cod_diversity.py'          → '08b_run_cod_diversity.py'"
echo "       '08_taxonomy_table.py'              → '07_taxonomy_table.py'"
echo "     Search for all script name strings:"
echo "       grep -n '\.py\"' scripts/pipeline.py"
echo ""
echo "  2. UPDATE run_all_figures.py"
echo "     Search for any script name references:"
echo "       grep -n '\.py' scripts/run_all_figures.py"
echo ""
echo "  3. UPDATE environment.yml comments"
echo "     Lines referencing script numbers 00-14 need updating."
echo ""
echo "  4. UPDATE docstrings in each renamed script"
echo "     Each script's docstring references its own name and pipeline position."
echo "     Key ones to update:"
echo "       07_taxonomy_table.py     (was 08_)"
echo "       08_run_diversity_stats.py (was 05_)"
echo "       09_plot_diversity.py      (was 06_)"
echo "       09c_visualize_diversity.py (internal docstring says 06_visualize_diversity.py)"
echo "       11_clean_diet_table.py    (was 09b_)"
echo "       11b_presence_absence.py   (was 08b_)"
echo "       11c_annotate_diet_ecology.py (was 10b_)"
echo ""
echo "  5. CREATE 07d_filter_feature_table.py"
echo "     Wraps qiime feature-table filter-features"
echo "     Input: confirmed_artefacts_{marker}.txt from 07c_"
echo "     See 07c_ report output for the exact qiime command."
echo ""
echo "  6. UPDATE README_public.md"
echo "     Pipeline step numbers in documentation need updating."
echo ""
echo "  7. COMMIT on a clean branch, verify pipeline.py check passes."
echo ""
echo "NEW PIPELINE ORDER:"
echo "  Stage 1  Setup:              01-04c"
echo "  Stage 2  Denoising:          05"
echo "  Stage 3  Rarefaction:        06"
echo "  Stage 4  Taxonomy + BLAST QC: 07-07d"
echo "  Stage 5  Diversity analysis: 08-08c"
echo "  Stage 6  Diversity figures:  09-09c"
echo "  Stage 7  Taxonomy plots:     10"
echo "  Stage 8  Diet analysis:      11-11c"
echo "  Stage 9  Ecological figures: 12-13"
echo "  Stage 10 Viral:              14-15"
echo "  Utilities (unnumbered): pipeline.py, run_all_figures.py,"
echo "            config_loader.py, audit_figures.sh,"
echo "            build_manuscript_figures.sh, herpes_threshold_sweep.py,"
echo "            scan_files.py"
