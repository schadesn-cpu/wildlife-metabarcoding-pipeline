#!/bin/bash
# Run this from your results/ directory:  bash audit_figures.sh

RESULTS_DIR="${1:-.}"
REPORT="figure_audit_report.txt"

echo "=====================================================" > $REPORT
echo "  FIGURE AUDIT REPORT" >> $REPORT
echo "  Generated: $(date)" >> $REPORT
echo "  Directory: $(realpath $RESULTS_DIR)" >> $REPORT
echo "=====================================================" >> $REPORT

# ── 1. ALL FIGURES WITH TIMESTAMPS ───────────────────────
echo "" >> $REPORT
echo "── ALL FIGURES (sorted by path) ──────────────────────" >> $REPORT
find "$RESULTS_DIR" -path "*/figures/*" \( -name "*.png" -o -name "*.svg" \) \
  | sort \
  | while read f; do
      printf "%-80s  %s\n" "$f" "$(stat -c '%y' "$f" | cut -d'.' -f1)"
    done >> $REPORT

# ── 2. PURPLE FILES (candidates for deletion) ────────────
echo "" >> $REPORT
echo "── PURPLE FILES TO REMOVE ────────────────────────────" >> $REPORT
find "$RESULTS_DIR" -path "*/figures/*" \( -name "*purple*" \) \
  | sort >> $REPORT
PURPLE_COUNT=$(find "$RESULTS_DIR" -path "*/figures/*" -name "*purple*" | wc -l)
echo "(Total: $PURPLE_COUNT files)" >> $REPORT

# ── 3. WONG FILES ─────────────────────────────────────────
echo "" >> $REPORT
echo "── WONG FILES TO KEEP ────────────────────────────────" >> $REPORT
find "$RESULTS_DIR" -path "*/figures/*" \( -name "*wong*" \) \
  | sort >> $REPORT
WONG_COUNT=$(find "$RESULTS_DIR" -path "*/figures/*" -name "*wong*" | wc -l)
echo "(Total: $WONG_COUNT files)" >> $REPORT

# ── 4. TIMESTAMP MISMATCH CHECK ───────────────────────────
# Flag any wong figure that is OLDER than its purple counterpart
echo "" >> $REPORT
echo "── TIMESTAMP MISMATCH (wong older than purple) ───────" >> $REPORT
MISMATCH=0
find "$RESULTS_DIR" -path "*/figures/*" -name "*wong*.png" -o \
                    -path "*/figures/*" -name "*wong*.svg" 2>/dev/null \
  | sort \
  | while read wong_file; do
      # Derive the purple equivalent path
      purple_file="${wong_file/wong/purple}"
      if [ -f "$purple_file" ]; then
          wong_time=$(stat -c '%Y' "$wong_file")
          purple_time=$(stat -c '%Y' "$purple_file")
          if [ "$wong_time" -lt "$purple_time" ]; then
              echo "  ⚠ OLDER WONG: $wong_file" >> $REPORT
              echo "         vs purple: $purple_file" >> $REPORT
              echo "         wong: $(stat -c '%y' "$wong_file" | cut -d'.' -f1)" >> $REPORT
              echo "         purple: $(stat -c '%y' "$purple_file" | cut -d'.' -f1)" >> $REPORT
              MISMATCH=$((MISMATCH+1))
          fi
      fi
    done
echo "" >> $REPORT

# ── 5. FIGURES WITH NO WONG EQUIVALENT ───────────────────
echo "── FIGURES WITH NO WONG VERSION (check these) ────────" >> $REPORT
find "$RESULTS_DIR" -path "*/figures/*" \( -name "*.png" -o -name "*.svg" \) \
  | grep -v "wong\|purple" \
  | sort >> $REPORT

# ── 6. SUMMARY ────────────────────────────────────────────
echo "" >> $REPORT
echo "── SUMMARY ───────────────────────────────────────────" >> $REPORT
TOTAL=$(find "$RESULTS_DIR" -path "*/figures/*" \( -name "*.png" -o -name "*.svg" \) | wc -l)
NO_SCHEME=$(find "$RESULTS_DIR" -path "*/figures/*" \( -name "*.png" -o -name "*.svg" \) | grep -v "wong\|purple" | wc -l)
echo "  Total figure files:           $TOTAL" >> $REPORT
echo "  Purple files (to remove):     $PURPLE_COUNT" >> $REPORT
echo "  Wong files (to keep):         $WONG_COUNT" >> $REPORT
echo "  Files with neither label:     $NO_SCHEME (review manually)" >> $REPORT

echo "" >> $REPORT
echo "=====================================================" >> $REPORT
echo "  NEXT STEP: If no mismatches flagged above," >> $REPORT
echo "  run cleanup_purple.sh to delete purple files." >> $REPORT
echo "=====================================================" >> $REPORT

echo "Audit complete. See: $REPORT"
cat $REPORT
