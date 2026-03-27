#!/usr/bin/env bash
# =============================================================================
# reorganize_results.sh
#
# Flattens results/ into a clean, consistent structure:
#
#   results/
#   ├── 16S/
#   │   ├── DvT/
#   │   │   ├── diversity/     ← all stat QZVs (permanova, permdisp, alpha sig)
#   │   │   ├── figures/       ← all PNGs/SVGs: {marker}_r{depth}_*_{palette}
#   │   │   └── taxonomy/      ← TSVs + barplot PNGs/SVGs
#   │   └── season/
#   │       ├── diversity/
#   │       └── figures/
#   ├── MiFish/  (same layout)
#   ├── cytb/    (same layout)
#   └── 18S/
#       └── all/
#           ├── diversity/
#           ├── figures/
#           └── taxonomy/
#
# BEFORE RUNNING: set the base results directory with --base, or edit the
# DEFAULT_BASE variable below to match your project path.
#
# Usage:
#   bash reorganize_results.sh --base /path/to/project/results --dry-run
#   bash reorganize_results.sh --base /path/to/project/results
#   # If DEFAULT_BASE is set below, --base can be omitted.
# =============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
# Edit DEFAULT_BASE to match your project's results directory, or pass
# --base /your/path at the command line to override without editing this file.
DEFAULT_BASE=""   # e.g. "/home/users/you/project/results"

BASE=""
DRY_RUN=false

for arg in "$@"; do
  case $arg in
    --dry-run)  DRY_RUN=true ;;
    --base)     shift; BASE="$1" ;;
    --base=*)   BASE="${arg#*=}" ;;
  esac
done

if [[ -z "$BASE" ]]; then
  BASE="$DEFAULT_BASE"
fi
if [[ -z "$BASE" ]]; then
  echo "[ERROR] No base directory specified." >&2
  echo "  Pass --base /path/to/results, or set DEFAULT_BASE in this script." >&2
  exit 1
fi
BASE="$(cd "$BASE" && pwd)"
if [[ ! -d "$BASE" ]]; then
  echo "[ERROR] Base directory does not exist: $BASE" >&2
  exit 1
fi

ARCHIVE="$BASE/_archive_old_structure"

# ── Helpers ───────────────────────────────────────────────────────────────────
log()  { echo "[INFO]  $*"; }
act()  { echo "[$(${DRY_RUN} && echo DRY || echo MOVE)] $*"; }

do_mkdir() {
  act "mkdir -p $1"
  $DRY_RUN || mkdir -p "$1"
}

do_mv() {
  if [ -e "$1" ]; then
    act "mv $1 -> $2"
    $DRY_RUN || mv "$1" "$2"
  else
    echo "[SKIP]  (not found) $1"
  fi
}

do_mv_glob() {
  local count
  count=$(find "$1" -maxdepth 1 -name "$2" 2>/dev/null | wc -l)
  if [ "$count" -gt 0 ]; then
    act "mv $1/$2 -> $3/"
    $DRY_RUN || find "$1" -maxdepth 1 -name "$2" -exec mv {} "$3/" \;
  else
    echo "[SKIP]  (none found) $1/$2"
  fi
}

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================="
echo "  Results Reorganization"
echo "  DRY_RUN  = $DRY_RUN"
echo "  BASE     = $BASE"
echo "  ARCHIVE  = $ARCHIVE"
echo "============================================="
echo ""

log "Creating new directory structure..."
do_mkdir "$ARCHIVE"

for MARKER in 16S MiFish cytb; do
  do_mkdir "$BASE/$MARKER/DvT/diversity"
  do_mkdir "$BASE/$MARKER/DvT/figures"
  do_mkdir "$BASE/$MARKER/DvT/taxonomy"
  do_mkdir "$BASE/$MARKER/season/diversity"
  do_mkdir "$BASE/$MARKER/season/figures"
done
do_mkdir "$BASE/18S/all/diversity"
do_mkdir "$BASE/18S/all/figures"
do_mkdir "$BASE/18S/all/taxonomy"
echo ""

# =============================================================================
# 16S
# =============================================================================
log "── 16S ──────────────────────────────────────────────────────────────────"
log "  DvT stat QZVs..."
do_mv_glob "$BASE/16S/rarefied_8000/DvT/diversity/16S/DvT/diversity" "*.qzv" "$BASE/16S/DvT/diversity"
do_mv_glob "$BASE/16S/rarefied_8000/DvT/diversity/16S/DvT/diversity" "*.sh"  "$BASE/16S/DvT/diversity"
log "  DvT figures -> archive..."
do_mv "$BASE/16S/rarefied_8000/DvT/figures/purple"         "$ARCHIVE/16S_DvT_purple"
do_mv "$BASE/16S/rarefied_8000/DvT/figures/purple_notitle" "$ARCHIVE/16S_DvT_purple_notitle"
do_mv_glob "$BASE/16S/rarefied_8000/DvT/figures" "*.png" "$ARCHIVE"
do_mv_glob "$BASE/16S/rarefied_8000/DvT/figures" "*.svg" "$ARCHIVE"
log "  Season figures -> archive..."
do_mv "$BASE/16S/rarefied_8000/season/figures_notitle" "$ARCHIVE/16S_season_notitle"
do_mv_glob "$BASE/16S/rarefied_8000/season/figures" "*.png" "$ARCHIVE/16S_season_old"
do_mv_glob "$BASE/16S/rarefied_8000/season/figures" "*.svg" "$ARCHIVE/16S_season_old"
log "  Taxonomy -> flat location..."
do_mv_glob "$BASE/16S/unrarefied/taxonomy"         "*.tsv" "$BASE/16S/DvT/taxonomy"
do_mv_glob "$BASE/16S/unrarefied/figures/taxonomy" "*.png" "$BASE/16S/DvT/taxonomy"
do_mv_glob "$BASE/16S/unrarefied/figures/taxonomy" "*.svg" "$BASE/16S/DvT/taxonomy"
echo ""

# =============================================================================
# MiFish
# =============================================================================
log "── MiFish ───────────────────────────────────────────────────────────────"
log "  DvT stat QZVs..."
do_mv_glob "$BASE/MiFish/rarefied_17000/DvT/diversity/MiFish/DvT/diversity" "*.qzv" "$BASE/MiFish/DvT/diversity"
do_mv_glob "$BASE/MiFish/rarefied_17000/DvT/diversity/MiFish/DvT/diversity" "*.sh"  "$BASE/MiFish/DvT/diversity"
log "  DvT figures -> archive..."
do_mv "$BASE/MiFish/rarefied_17000/DvT/figures/purple"         "$ARCHIVE/MiFish_DvT_purple"
do_mv "$BASE/MiFish/rarefied_17000/DvT/figures/purple_notitle" "$ARCHIVE/MiFish_DvT_purple_notitle"
do_mv_glob "$BASE/MiFish/rarefied_17000/DvT/figures" "*.png" "$ARCHIVE"
do_mv_glob "$BASE/MiFish/rarefied_17000/DvT/figures" "*.svg" "$ARCHIVE"
log "  Season figures -> archive..."
do_mv "$BASE/MiFish/rarefied_17000/season/figures_notitle" "$ARCHIVE/MiFish_season_notitle"
do_mv_glob "$BASE/MiFish/rarefied_17000/season/figures" "*.png" "$ARCHIVE/MiFish_season_old"
do_mv_glob "$BASE/MiFish/rarefied_17000/season/figures" "*.svg" "$ARCHIVE/MiFish_season_old"
log "  Taxonomy..."
do_mv_glob "$BASE/MiFish/unrarefied/taxonomy"         "*.tsv" "$BASE/MiFish/DvT/taxonomy"
do_mv_glob "$BASE/MiFish/unrarefied/figures/taxonomy" "*.png" "$BASE/MiFish/DvT/taxonomy"
do_mv_glob "$BASE/MiFish/unrarefied/figures/taxonomy" "*.svg" "$BASE/MiFish/DvT/taxonomy"
echo ""

# =============================================================================
# cytb
# =============================================================================
log "── cytb ─────────────────────────────────────────────────────────────────"
log "  DvT stat QZVs (merging two source locations)..."
do_mv_glob "$BASE/cytb/all/diversity" "*.qzv" "$BASE/cytb/DvT/diversity"
do_mv_glob "$BASE/cytb/all/diversity" "*.sh"  "$BASE/cytb/DvT/diversity"
do_mv_glob "$BASE/cytb/rarefied_200/DvT/diversity/cytb/DvT/diversity" "*.qzv" "$BASE/cytb/DvT/diversity"
do_mv_glob "$BASE/cytb/rarefied_200/DvT/diversity/cytb/DvT/diversity" "*.sh"  "$BASE/cytb/DvT/diversity"
log "  DvT figures -> archive..."
do_mv "$BASE/cytb/rarefied_200/DvT/figures/purple"         "$ARCHIVE/cytb_DvT_purple"
do_mv "$BASE/cytb/rarefied_200/DvT/figures/purple_notitle" "$ARCHIVE/cytb_DvT_purple_notitle"
do_mv "$BASE/cytb/all/figures/diversity"                   "$ARCHIVE/cytb_all_figures_old"
log "  Taxonomy..."
do_mv_glob "$BASE/cytb/all/taxonomy" "*.tsv" "$BASE/cytb/DvT/taxonomy"
echo ""

# =============================================================================
# 18S
# =============================================================================
log "── 18S ──────────────────────────────────────────────────────────────────"
log "  Stat QZVs (nested fresh run)..."
do_mv_glob "$BASE/18S/all/diversity/18S/all/diversity" "*.qzv" "$BASE/18S/all/diversity"
do_mv_glob "$BASE/18S/all/diversity/18S/all/diversity" "*.sh"  "$BASE/18S/all/diversity"
log "  Figures -> archive..."
do_mv "$BASE/18S/all/figures/purple"         "$ARCHIVE/18S_all_purple"
do_mv "$BASE/18S/all/figures/purple_notitle" "$ARCHIVE/18S_all_purple_notitle"
do_mv "$BASE/18S/all/figures/diversity"      "$ARCHIVE/18S_all_figures_old"
log "  Taxonomy..."
do_mv_glob "$BASE/18S/all/taxonomy"         "*.tsv" "$BASE/18S/all/taxonomy"
do_mv_glob "$BASE/18S/all/figures/taxonomy" "*.png" "$BASE/18S/all/taxonomy"
do_mv_glob "$BASE/18S/all/figures/taxonomy" "*.svg" "$BASE/18S/all/taxonomy"

# =============================================================================
# Done
# =============================================================================
echo ""
log "Done. New flat structure:"
echo ""
echo "  results/"
echo "  ├── 16S/DvT/{diversity,figures,taxonomy}/"
echo "  ├── 16S/season/{diversity,figures}/"
echo "  ├── MiFish/  (same)"
echo "  ├── cytb/    (same)"
echo "  ├── 18S/all/{diversity,figures,taxonomy}/"
echo "  └── _archive_old_structure/  ← old files, safe to delete later"
echo ""
if $DRY_RUN; then
  echo "  DRY RUN — nothing was changed. Re-run without --dry-run to apply."
fi
echo "============================================="
