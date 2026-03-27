#!/usr/bin/env bash
# =============================================================================
# reorganize_loon.sh
# Reorganizes loon_project/results/ into a clean rarefied/unrarefied structure.
#
# BEFORE RUNNING: set the base results directory with --base, or edit the
# DEFAULT_BASE variable below to match your project path.
#
# Usage:
#   bash reorganize_loon.sh --base /path/to/loon_project/results --dry-run
#   bash reorganize_loon.sh --base /path/to/loon_project/results
#
#   # If DEFAULT_BASE is set correctly below, --base can be omitted:
#   bash reorganize_loon.sh --dry-run
#   bash reorganize_loon.sh
# =============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
# Edit DEFAULT_BASE to match your project's results directory, or pass
# --base /your/path at the command line to override without editing this file.
DEFAULT_BASE=""   # e.g. "/home/users/you/meedlab/loon_project/results"

DRY_RUN=false
BASE=""

# Parse flags
for arg in "$@"; do
  case $arg in
    --dry-run)    DRY_RUN=true ;;
    --base)       shift; BASE="$1" ;;
    --base=*)     BASE="${arg#*=}" ;;
  esac
done

# Fall back to DEFAULT_BASE if --base was not supplied
if [[ -z "$BASE" ]]; then
  BASE="$DEFAULT_BASE"
fi

# Abort if we still have no path
if [[ -z "$BASE" ]]; then
  echo "[ERROR] No base directory specified." >&2
  echo "  Pass --base /path/to/loon_project/results, or set DEFAULT_BASE in this script." >&2
  exit 1
fi

# Resolve to absolute path and verify it exists
BASE="$(cd "$BASE" && pwd)"
if [[ ! -d "$BASE" ]]; then
  echo "[ERROR] Base directory does not exist: $BASE" >&2
  exit 1
fi

# ── Helpers ───────────────────────────────────────────────────────────────────
log()    { echo "[INFO]  $*"; }
action() { echo "[$(${DRY_RUN} && echo DRY-RUN || echo ACTION)] $*"; }

run_cmd() {
  # Print the command always; only execute if not dry-run
  action "$*"
  if ! $DRY_RUN; then
    eval "$*"
  fi
}

make_dir() {
  action "mkdir -p $1"
  if ! $DRY_RUN; then
    mkdir -p "$1"
  fi
}

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================="
echo "  Loon Project Results Reorganization Script"
echo "  DRY_RUN = $DRY_RUN"
echo "  BASE    = $BASE"
echo "============================================="
echo ""

# =============================================================================
# STEP 1: Create new directory structure
# =============================================================================
log "STEP 1: Creating new folder structure..."

for MARKER in 16S MiFish; do
  make_dir "$BASE/$MARKER/unrarefied/taxonomy"
  make_dir "$BASE/$MARKER/unrarefied/figures/taxonomy"
  make_dir "$BASE/$MARKER/rarefied_17000/DvT/diversity"
  make_dir "$BASE/$MARKER/rarefied_17000/DvT/figures/purple"
  make_dir "$BASE/$MARKER/rarefied_17000/DvT/figures/highcontrast"
  make_dir "$BASE/$MARKER/rarefied_17000/season/diversity"
  make_dir "$BASE/$MARKER/rarefied_17000/season/figures/purple"
  make_dir "$BASE/$MARKER/rarefied_17000/season/figures/highcontrast"
done

echo ""

# =============================================================================
# STEP 2: Move MiFish unrarefied taxonomy files
# =============================================================================
log "STEP 2: Moving MiFish unrarefied taxonomy files..."

MIFISH_TAX_SRC="$BASE/MiFish/all/taxonomy"
MIFISH_TAX_DST="$BASE/MiFish/unrarefied/taxonomy"

for FILE in \
  taxonomy_summary_MiFish.tsv \
  taxonomy_counts_L7_MiFish.tsv \
  taxonomy_relabund_L7_MiFish.tsv \
  taxonomy_top30_L7_MiFish.tsv; do
  run_cmd "mv '$MIFISH_TAX_SRC/$FILE' '$MIFISH_TAX_DST/$FILE'"
done

echo ""

# =============================================================================
# STEP 3: Move MiFish barplots
# =============================================================================
log "STEP 3: Moving MiFish barplots to unrarefied/figures/taxonomy/..."

MIFISH_FIG_SRC="$BASE/MiFish/all/figures/taxonomy"
MIFISH_FIG_DST="$BASE/MiFish/unrarefied/figures/taxonomy"

for FILE in \
  barplot_MiFish_Group_purple.png \
  barplot_MiFish_Group_purple.svg; do
  run_cmd "mv '$MIFISH_FIG_SRC/$FILE' '$MIFISH_FIG_DST/$FILE'"
done

echo ""

# =============================================================================
# STEP 4: Move 16S unrarefied taxonomy files
# =============================================================================
log "STEP 4: Moving 16S unrarefied taxonomy files..."

S16_TAX_SRC="$BASE/16S/DvT/taxonomy"
S16_TAX_DST="$BASE/16S/unrarefied/taxonomy"

for FILE in \
  taxonomy_summary_16S.tsv \
  taxonomy_counts_L6_16S.tsv \
  taxonomy_relabund_L6_16S.tsv \
  taxonomy_top30_L6_16S.tsv; do
  run_cmd "mv '$S16_TAX_SRC/$FILE' '$S16_TAX_DST/$FILE'"
done

echo ""

# =============================================================================
# STEP 5: Move 16S barplots
# =============================================================================
log "STEP 5: Moving 16S barplots to unrarefied/figures/taxonomy/..."

S16_FIG_SRC="$BASE/16S/DvT/figures/taxonomy"
S16_FIG_DST="$BASE/16S/unrarefied/figures/taxonomy"

for FILE in \
  barplot_16S_Group_purple.png \
  barplot_16S_Group_purple.svg; do
  run_cmd "mv '$S16_FIG_SRC/$FILE' '$S16_FIG_DST/$FILE'"
done

echo ""

# =============================================================================
# STEP 6: Move rarefied MiFish diversity figures to proper location
#         (purple/ and highcontrast/ subfolders are the ones to keep)
# =============================================================================
log "STEP 6: Moving MiFish rarefied diversity figures (DvT)..."

CORE_BASE="$BASE/MiFish/all/figures/core_metrics_17000/diversity"
RAR_DST="$BASE/MiFish/rarefied_17000/DvT/figures"

for COLOR in purple highcontrast; do
  for FILE in alpha_panel_Group.png alpha_panel_Group.svg \
              pcoa_panel_Group.png pcoa_panel_Group.svg; do
    SRC="$CORE_BASE/$COLOR/$FILE"
    DST="$RAR_DST/$COLOR/$FILE"
    run_cmd "mv '$SRC' '$DST'"
  done
done

echo ""

# =============================================================================
# STEP 7: Delete stale files
# =============================================================================
log "STEP 7: Deleting stale files..."

STALE_FILES=(
  # Stale genus counts (Mar 7, predates classifier fix)
  "$BASE/16S/DvT/taxonomy/genus_counts_16S_DvT.tsv"

  # Old presence/absence export (Feb 23)
  "$BASE/16S/export/genus_pa_export/genus_pa.tsv"

  # Old run config (Mar 6)
  "$BASE/MiFish/all/run_scope.json"

  # Duplicate parent-level figures (purple/ and highcontrast/ copies are the keepers)
  "$CORE_BASE/pcoa_panel_Group.png"
  "$CORE_BASE/pcoa_panel_Group.svg"
  "$CORE_BASE/alpha_panel_Group.png"
  "$CORE_BASE/alpha_panel_Group.svg"
)

for FILE in "${STALE_FILES[@]}"; do
  run_cmd "rm -f '$FILE'"
done

echo ""

# =============================================================================
# STEP 8: Clean up now-empty old directories
# =============================================================================
log "STEP 8: Removing now-empty old directories..."

OLD_DIRS=(
  "$BASE/MiFish/all/figures/core_metrics_17000/diversity/purple"
  "$BASE/MiFish/all/figures/core_metrics_17000/diversity/highcontrast"
  "$BASE/MiFish/all/figures/core_metrics_17000/diversity"
  "$BASE/MiFish/all/figures/core_metrics_17000"
  "$BASE/MiFish/all/figures/taxonomy"
  "$BASE/MiFish/all/figures"
  "$BASE/MiFish/all/taxonomy"
  "$BASE/MiFish/all"
  "$BASE/16S/DvT/figures/taxonomy"
  "$BASE/16S/DvT/figures"
  "$BASE/16S/DvT/taxonomy"
  "$BASE/16S/DvT"
  "$BASE/16S/export/genus_pa_export"
  "$BASE/16S/export"
)

for DIR in "${OLD_DIRS[@]}"; do
  # Only remove if directory is empty
  run_cmd "rmdir --ignore-fail-on-non-empty '$DIR'"
done

echo ""

# =============================================================================
# Done
# =============================================================================
echo "============================================="
if $DRY_RUN; then
  echo "  DRY RUN COMPLETE — nothing was changed."
  echo "  Re-run without --dry-run to apply changes."
else
  echo "  REORGANIZATION COMPLETE."
  echo ""
  echo "  New structure:"
  echo "  results/"
  echo "  ├── 16S/"
  echo "  │   ├── unrarefied/taxonomy/       ← taxonomy TSVs + barplots"
  echo "  │   └── rarefied_17000/DvT/        ← ready for diversity rerun"
  echo "  │                    /season/"
  echo "  └── MiFish/"
  echo "      ├── unrarefied/taxonomy/"
  echo "      └── rarefied_17000/DvT/"
  echo "                       /season/"
fi
echo "============================================="
