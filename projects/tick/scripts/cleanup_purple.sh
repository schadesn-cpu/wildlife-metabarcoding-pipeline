#!/usr/bin/env bash
# ============================================================
#  cleanup_purple.sh
#  Removes all _purple figure files from loon_project.
#  Run from: your project root directory
#
#  BEFORE RUNNING:
#    1. Regenerate adenovirus and herpesvirus wong figures
#       (purple versions are 3hrs newer than wong on 2026-03-31)
#    2. Run in DRY_RUN=1 mode first to confirm list
#
#  Usage:
#    DRY_RUN=1 bash cleanup_purple.sh   # preview only
#    DRY_RUN=0 bash cleanup_purple.sh   # actually delete
# ============================================================

set -euo pipefail
DRY_RUN="${DRY_RUN:-1}"
ROOT="${ROOT:-$(pwd)}"

echo "=============================================="
echo "  PURPLE FILE CLEANUP"
echo "  DRY_RUN=${DRY_RUN}  (set to 0 to actually delete)"
echo "  Root: ${ROOT}"
echo "=============================================="

# ── SAFETY: Warn if adeno/herpes wong files are still stale ──
echo ""
echo "── STALE WONG CHECK (adeno + herpes) ─────────"
STALE_WARN=0
for pair in \
    "results/adenovirus/figures/adeno_Group_wong_presence.png:results/adenovirus/figures/adeno_Group_purple_presence.png" \
    "results/adenovirus/figures/adeno_Group_wong_relabund.png:results/adenovirus/figures/adeno_Group_purple_relabund.png" \
    "results/herpesvirus/figures/herpes_Group_wong_presence.png:results/herpesvirus/figures/herpes_Group_purple_presence.png" \
    "results/herpesvirus/figures/herpes_Group_wong_relabund.png:results/herpesvirus/figures/herpes_Group_purple_relabund.png"
do
    WONG="${ROOT}/${pair%%:*}"
    PURPLE="${ROOT}/${pair##*:}"
    if [[ -f "$WONG" && -f "$PURPLE" ]]; then
        WT=$(stat -c %Y "$WONG")
        PT=$(stat -c %Y "$PURPLE")
        if [[ $WT -lt $PT ]]; then
            echo "  ⚠  STALE: $(basename $WONG)"
            echo "     Wong:   $(stat -c %y $WONG | cut -d. -f1)"
            echo "     Purple: $(stat -c %y $PURPLE | cut -d. -f1)"
            STALE_WARN=1
        else
            echo "  ✅  OK: $(basename $WONG) is newer than purple"
        fi
    fi
done

if [[ $STALE_WARN -eq 1 ]]; then
    echo ""
    echo "  ⛔  STALE WONG FILES DETECTED."
    echo "  Regenerate adeno + herpes wong figures first, then rerun."
    echo "  Run 10_plot_viral.py with wong palette for both viral markers."
    echo "  Aborting cleanup."
    exit 1
fi

# ── COLLECT ALL PURPLE FILES ─────────────────────────────────
echo ""
echo "── FILES TO DELETE ────────────────────────────"
mapfile -t PURPLE_FILES < <(find "${ROOT}/results" -type f \( -name "*_purple.png" -o -name "*_purple.svg" \) | sort)

if [[ ${#PURPLE_FILES[@]} -eq 0 ]]; then
    echo "  No purple files found. Nothing to do."
    exit 0
fi

COUNT=0
for f in "${PURPLE_FILES[@]}"; do
    echo "  DEL: ${f#$ROOT/}"
    COUNT=$((COUNT + 1))
done

echo ""
echo "  Total to delete: ${COUNT} files"

# ── DELETE ───────────────────────────────────────────────────
if [[ $DRY_RUN -eq 1 ]]; then
    echo ""
    echo "  DRY RUN — nothing deleted."
    echo "  Rerun with DRY_RUN=0 to delete."
else
    echo ""
    echo "  Deleting..."
    for f in "${PURPLE_FILES[@]}"; do
        rm -f "$f"
        echo "  ✅  Removed: ${f#$ROOT/}"
    done
    echo ""
    echo "  Done. ${COUNT} files deleted."
fi

echo "=============================================="
