#!/usr/bin/env bash
# =============================================================================
# run_all_figures.sh
#
# Generates all diversity figures (PCoA + alpha) in both purple and wong
# palettes for every marker and dataset combination.
#
# Run AFTER reorganize_results.sh has been applied and stat QZVs are in place.
#
# Output naming convention:
#   results/{marker}/{dataset}/figures/{marker}_r{depth}_{dataset}_{type}_{palette}.png/svg
#
# Examples:
#   16S_r8000_DvT_pcoa_purple.png
#   MiFish_r17000_season_alpha_wong.svg
#   cytb_r200_DvT_pcoa_wong.png
#   18S_r1000_all_alpha_purple.png
#
# Usage:
#   # Run from project root (directory containing qiime2/, results/, scripts/)
#   bash scripts/run_all_figures.sh
#   bash scripts/run_all_figures.sh --dry-run   # preview commands only
#   bash scripts/run_all_figures.sh --project-root /path/to/project
# =============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_ROOT="$(pwd)"
DRY_RUN=false

for arg in "$@"; do
  case $arg in
    --dry-run)         DRY_RUN=true ;;
    --project-root)    shift; PROJECT_ROOT="$1" ;;
    --project-root=*)  PROJECT_ROOT="${arg#*=}" ;;
  esac
done

PROJECT_ROOT="$(cd "$PROJECT_ROOT" && pwd)"
BASE="$PROJECT_ROOT/qiime2"
META="$PROJECT_ROOT/metadata/qiime"
RES="$PROJECT_ROOT/results"
SCRIPT="$PROJECT_ROOT/scripts/06_plot_diversity.py"

# ── Pre-flight checks ─────────────────────────────────────────────────────────
if [[ ! -f "$SCRIPT" ]]; then
  echo "[ERROR] Cannot find 06_plot_diversity.py at: $SCRIPT" >&2
  echo "  Run from project root, or pass --project-root /path/to/project" >&2
  exit 1
fi

if ! python3 -c "import matplotlib, numpy, pandas, scipy" 2>/dev/null; then
  echo "[ERROR] Required Python packages missing (matplotlib, numpy, pandas, scipy)." >&2
  echo "  Activate the analysis environment first:" >&2
  echo "    conda activate metabarcoding-analysis" >&2
  exit 1
fi

echo ""
echo "============================================="
echo "  run_all_figures.sh"
echo "  DRY_RUN      = $DRY_RUN"
echo "  PROJECT_ROOT = $PROJECT_ROOT"
echo "============================================="
echo ""

# ── Helper: generate pcoa + alpha for one marker/dataset/grouping ─────────────
# Usage: plot_marker <marker> <depth> <dataset> <group_col> <metrics_dir> <pcoa_metrics...> -- <alpha_metrics...>
#
# This function handles both the pcoa and alpha subcommands, looping over
# both palettes, and skips gracefully if the metrics directory is missing.
#
# Arguments:
#   $1  marker name (e.g. 16S, MiFish)
#   $2  rarefaction depth (e.g. 8000)
#   $3  dataset slug (e.g. DvT, season, all)
#   $4  grouping column (e.g. Group, Season)
#   $5  metrics directory (absolute path to core-metrics output)
#   $6+ pcoa QZA basenames (e.g. bray_curtis_pcoa_results.qza), then "--", then alpha QZA basenames
#
plot_marker() {
  local MARKER="$1"
  local DEPTH="$2"
  local DATASET="$3"
  local GROUP="$4"
  local METRICS="$5"
  shift 5

  # Split remaining args on "--" into pcoa and alpha lists
  local PCOA_QZAS=()
  local ALPHA_QZAS=()
  local in_alpha=false
  for arg in "$@"; do
    if [[ "$arg" == "--" ]]; then
      in_alpha=true
    elif $in_alpha; then
      ALPHA_QZAS+=("$METRICS/$arg")
    else
      PCOA_QZAS+=("$METRICS/$arg")
    fi
  done

  local OUT_DIR="$RES/$MARKER/$DATASET/figures"
  local META_FILE="$META/metadata_${MARKER}.tsv"
  local STATS_DIR="$RES/$MARKER/$DATASET/diversity"
  local STEM_BASE="${MARKER}_r${DEPTH}_${DATASET}"

  # Check metrics dir exists
  if [[ ! -d "$METRICS" ]]; then
    echo "[SKIP] Metrics directory not found: $METRICS"
    return 0
  fi

  # Check metadata exists
  if [[ ! -f "$META_FILE" ]]; then
    echo "[SKIP] Metadata not found: $META_FILE"
    return 0
  fi

  # Verify all QZA files exist before starting
  local missing=0
  for qza in "${PCOA_QZAS[@]}" "${ALPHA_QZAS[@]}"; do
    if [[ ! -f "$qza" ]]; then
      echo "[WARN] QZA not found: $qza"
      missing=$((missing + 1))
    fi
  done
  if [[ $missing -gt 0 ]]; then
    echo "[SKIP] $missing QZA file(s) missing for $MARKER/$DATASET — skipping this block."
    return 0
  fi

  mkdir -p "$OUT_DIR"

  for PALETTE in purple wong; do
    echo "── $MARKER $DATASET ($GROUP) $PALETTE ──"

    local pcoa_cmd=(python3 "$SCRIPT" pcoa
      --artifact "${PCOA_QZAS[@]}"
      --metadata "$META_FILE"
      --color-by "$GROUP" --panel --palette "$PALETTE" --no-title
      --output-stem "${STEM_BASE}_pcoa_${PALETTE}"
      --output-dir "$OUT_DIR"
    )
    # Add --stats-dir only if it exists (season analyses may not have stat QZVs yet)
    if [[ -d "$STATS_DIR" ]]; then
      pcoa_cmd+=(--stats-dir "$STATS_DIR")
    fi

    local alpha_cmd=(python3 "$SCRIPT" alpha
      --artifact "${ALPHA_QZAS[@]}"
      --metadata "$META_FILE"
      --group-by "$GROUP" --panel --palette "$PALETTE" --no-title
      --output-stem "${STEM_BASE}_alpha_${PALETTE}"
      --output-dir "$OUT_DIR"
    )

    if $DRY_RUN; then
      echo "  [DRY] ${pcoa_cmd[*]}"
      echo "  [DRY] ${alpha_cmd[*]}"
    else
      "${pcoa_cmd[@]}"
      "${alpha_cmd[@]}"
    fi
  done
  echo ""
}

# =============================================================================
# 16S  (phylogenetic; depth 8000)
# =============================================================================
M16S="$BASE/16S/rarefied_8000/DvT/diversity/core_metrics_depth8000"

plot_marker 16S 8000 DvT Group "$M16S" \
  unweighted_unifrac_pcoa_results.qza \
  weighted_unifrac_pcoa_results.qza \
  bray_curtis_pcoa_results.qza \
  jaccard_pcoa_results.qza \
  -- \
  faith_pd_vector.qza \
  shannon_vector.qza \
  observed_features_vector.qza \
  evenness_vector.qza

plot_marker 16S 8000 season Season "$M16S" \
  unweighted_unifrac_pcoa_results.qza \
  weighted_unifrac_pcoa_results.qza \
  bray_curtis_pcoa_results.qza \
  jaccard_pcoa_results.qza \
  -- \
  faith_pd_vector.qza \
  shannon_vector.qza \
  observed_features_vector.qza \
  evenness_vector.qza

# =============================================================================
# MiFish  (non-phylogenetic; depth 17000)
# =============================================================================
MMIFISH="$BASE/MiFish/rarefied_17000/DvT/diversity/core-metrics-17000"

plot_marker MiFish 17000 DvT Group "$MMIFISH" \
  bray_curtis_pcoa_results.qza \
  jaccard_pcoa_results.qza \
  -- \
  observed_features_vector.qza \
  shannon_vector.qza \
  evenness_vector.qza

plot_marker MiFish 17000 season Season "$MMIFISH" \
  bray_curtis_pcoa_results.qza \
  jaccard_pcoa_results.qza \
  -- \
  observed_features_vector.qza \
  shannon_vector.qza \
  evenness_vector.qza

# =============================================================================
# cytb  (non-phylogenetic; depth 200)
# =============================================================================
MCYTB="$BASE/cytb/all/diversity/core-metrics-200"

plot_marker cytb 200 DvT Group "$MCYTB" \
  bray_curtis_pcoa_results.qza \
  jaccard_pcoa_results.qza \
  -- \
  observed_features_vector.qza \
  shannon_vector.qza \
  evenness_vector.qza

plot_marker cytb 200 season Season "$MCYTB" \
  bray_curtis_pcoa_results.qza \
  jaccard_pcoa_results.qza \
  -- \
  observed_features_vector.qza \
  shannon_vector.qza \
  evenness_vector.qza

# =============================================================================
# 18S  (non-phylogenetic; depth 1000)
# =============================================================================
M18S="$BASE/18S/all/diversity/core-metrics-1000"

plot_marker 18S 1000 all Group "$M18S" \
  bray_curtis_pcoa_results.qza \
  jaccard_pcoa_results.qza \
  -- \
  observed_features_vector.qza \
  shannon_vector.qza \
  evenness_vector.qza

plot_marker 18S 1000 season Season "$M18S" \
  bray_curtis_pcoa_results.qza \
  jaccard_pcoa_results.qza \
  -- \
  observed_features_vector.qza \
  shannon_vector.qza \
  evenness_vector.qza

# =============================================================================
# Done
# =============================================================================
echo "============================================="
echo "  All figures complete."
echo ""
echo "  Output locations:"
echo "    results/16S/DvT/figures/"
echo "    results/16S/season/figures/"
echo "    results/MiFish/DvT/figures/"
echo "    results/MiFish/season/figures/"
echo "    results/cytb/DvT/figures/"
echo "    results/cytb/season/figures/"
echo "    results/18S/all/figures/"
echo ""
echo "  Naming: {marker}_r{depth}_{dataset}_{type}_{palette}.png/svg"
if $DRY_RUN; then
  echo ""
  echo "  DRY RUN — no files written. Re-run without --dry-run to apply."
fi
echo "============================================="
