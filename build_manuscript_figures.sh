#!/usr/bin/env bash
# ============================================================
#  build_manuscript_figures.sh
#  Creates a flat manuscript_figures/ directory with clean,
#  paper-ready names. Copies (not moves) from results/.
#  Run AFTER cleanup_purple.sh.
#
#  Run from: /home/users/sne24/meedlab/loon_project
#  Output:   /home/users/sne24/meedlab/loon_project/manuscript_figures/
#
#  Usage:
#    bash build_manuscript_figures.sh
# ============================================================

set -euo pipefail
ROOT="/home/users/sne24/meedlab/loon_project"
OUT="${ROOT}/manuscript_figures"

mkdir -p "$OUT"

echo "=============================================="
echo "  BUILD MANUSCRIPT FIGURES"
echo "  Output: ${OUT}"
echo "=============================================="

cp_fig() {
    local src="${ROOT}/$1"
    local dst="${OUT}/$2"
    if [[ -f "$src" ]]; then
        cp "$src" "$dst"
        echo "  ✅  $2"
    else
        echo "  ⚠   MISSING SOURCE: $1"
        echo "      → $2  [NOT COPIED — file does not exist]"
        # Write a placeholder text file
        local base="${dst%.png}"
        echo "PLACEHOLDER — source file missing: $1" > "${base}_MISSING.txt"
    fi
}

echo ""
echo "── MAIN FIGURES ───────────────────────────────"

# Fig 1 — Multimarker alpha 4-panel (opening figure)
cp_fig \
    "results/multimarker/figures/multimarker_alpha_observed_4panel_wong.png" \
    "Fig01_multimarker_alpha_4panel.png"
cp_fig \
    "results/multimarker/figures/multimarker_alpha_observed_4panel_wong.svg" \
    "Fig01_multimarker_alpha_4panel.svg"

# Fig 2 — 16S DvT PCoA
cp_fig \
    "results/16S/DvT/figures/16S_r8000_DvT_pcoa_wong.png" \
    "Fig02_16S_DvT_pcoa.png"
cp_fig \
    "results/16S/DvT/figures/16S_r8000_DvT_pcoa_wong.svg" \
    "Fig02_16S_DvT_pcoa.svg"

# Fig 3 — MiFish DvT PCoA (primary dietary result)
cp_fig \
    "results/MiFish/DvT/figures/MiFish_r17000_DvT_pcoa_wong.png" \
    "Fig03_MiFish_DvT_pcoa.png"
cp_fig \
    "results/MiFish/DvT/figures/MiFish_r17000_DvT_pcoa_wong.svg" \
    "Fig03_MiFish_DvT_pcoa.svg"

# Fig 4 — cytb DvT PCoA
cp_fig \
    "results/cytb/DvT/figures/cytb_r200_DvT_pcoa_wong.png" \
    "Fig04_cytb_DvT_pcoa.png"
cp_fig \
    "results/cytb/DvT/figures/cytb_r200_DvT_pcoa_wong.svg" \
    "Fig04_cytb_DvT_pcoa.svg"

# Fig 5 — 18S DvT PCoA  ⚠ PENDING 18S decision
cp_fig \
    "results/18S/all/figures/18S_r1000_all_pcoa_wong.png" \
    "Fig05_18S_DvT_pcoa_PENDING18S.png"
cp_fig \
    "results/18S/all/figures/18S_r1000_all_pcoa_wong.svg" \
    "Fig05_18S_DvT_pcoa_PENDING18S.svg"

# Fig 6 — COD alpha (MiFish, filtered = Unknown_Other removed)
cp_fig \
    "results/MiFish/COD/figures/MiFish_COD_alpha_wong_filtered.png" \
    "Fig06_MiFish_COD_alpha.png"
cp_fig \
    "results/MiFish/COD/figures/MiFish_COD_alpha_wong_filtered.svg" \
    "Fig06_MiFish_COD_alpha.svg"

# Fig 7 — COD PCoA (MiFish, filtered)
cp_fig \
    "results/MiFish/COD/figures/MiFish_COD_pcoa_wong_filtered.png" \
    "Fig07_MiFish_COD_pcoa.png"
cp_fig \
    "results/MiFish/COD/figures/MiFish_COD_pcoa_wong_filtered.svg" \
    "Fig07_MiFish_COD_pcoa.svg"

# Fig 8 — MiFish habitat by season (pipeline validation, PRIMARY version)
cp_fig \
    "results/MiFish/all/figures/habitat_season/MiFish_habitat_season_Diseased_Trauma_identified_only_wong.png" \
    "Fig08_MiFish_habitat_season.png"
cp_fig \
    "results/MiFish/all/figures/habitat_season/MiFish_habitat_season_Diseased_Trauma_identified_only_wong.svg" \
    "Fig08_MiFish_habitat_season.svg"

# Fig 9 — Herpesvirus relative abundance
cp_fig \
    "results/herpesvirus/figures/herpes_Group_wong_relabund.png" \
    "Fig09_herpes_relabund.png"
cp_fig \
    "results/herpesvirus/figures/herpes_Group_wong_relabund.svg" \
    "Fig09_herpes_relabund.svg"

# Fig 10 — Herpesvirus presence/absence
cp_fig \
    "results/herpesvirus/figures/herpes_Group_wong_presence.png" \
    "Fig10_herpes_presence.png"
cp_fig \
    "results/herpesvirus/figures/herpes_Group_wong_presence.svg" \
    "Fig10_herpes_presence.svg"

# Fig 11 — Adenovirus relative abundance
cp_fig \
    "results/adenovirus/figures/adeno_Group_wong_relabund.png" \
    "Fig11_adeno_relabund.png"
cp_fig \
    "results/adenovirus/figures/adeno_Group_wong_relabund.svg" \
    "Fig11_adeno_relabund.svg"

# Fig 12 — Adenovirus presence/absence
cp_fig \
    "results/adenovirus/figures/adeno_Group_wong_presence.png" \
    "Fig12_adeno_presence.png"
cp_fig \
    "results/adenovirus/figures/adeno_Group_wong_presence.svg" \
    "Fig12_adeno_presence.svg"

# Fig 13 — Adenovirus phylogenetic tree (check path — not in audit)
cp_fig \
    "results/adenovirus/figures/adeno_phylotree_wong.png" \
    "Fig13_adeno_phylotree.png"
cp_fig \
    "results/adenovirus/figures/adeno_phylotree_wong.svg" \
    "Fig13_adeno_phylotree.svg"

echo ""
echo "── SUPPLEMENTARY FIGURES ──────────────────────"

# Fig S1 — MiFish COD all-groups PCoA (Lead/Marine/PI/Trauma)
cp_fig \
    "results/MiFish/COD/figures/MiFish_COD_allgroups_pcoa_wong.png" \
    "FigS1_MiFish_COD_allgroups_pcoa.png"
cp_fig \
    "results/MiFish/COD/figures/MiFish_COD_allgroups_pcoa_wong.svg" \
    "FigS1_MiFish_COD_allgroups_pcoa.svg"

# Fig S2 — 16S DvT + Season taxonomy barplot (unrarefied)
# ⚠ The meteorological-season version exists but has wrong labels
# Use this only if regenerated with ecological seasons
cp_fig \
    "results/16S/season/figures/16S_unrarefied_season_barplot_wong.png" \
    "FigS2_16S_DvT_season_barplot_WRONG_SEASON_LABELS.png"

# Fig S3 — MiFish DvT + Season taxonomy barplot
cp_fig \
    "results/MiFish/season/figures/MiFish_unrarefied_season_barplot_wong.png" \
    "FigS3_MiFish_DvT_season_barplot.png"
cp_fig \
    "results/MiFish/season/figures/MiFish_unrarefied_season_barplot_wong.svg" \
    "FigS3_MiFish_DvT_season_barplot.svg"

# Fig S4 — cytb taxonomy barplot
cp_fig \
    "results/cytb/season/figures/cytb_unrarefied_season_barplot_wong.png" \
    "FigS4_cytb_season_barplot.png"
cp_fig \
    "results/cytb/season/figures/cytb_unrarefied_season_barplot_wong.svg" \
    "FigS4_cytb_season_barplot.svg"

# Fig S5 — 18S taxonomy barplot ⚠ PENDING 18S decision
cp_fig \
    "results/18S/all/figures/18S_unrarefied_DvT_barplot_wong.png" \
    "FigS5_18S_DvT_barplot_PENDING18S.png"
cp_fig \
    "results/18S/all/figures/18S_unrarefied_DvT_barplot_wong.svg" \
    "FigS5_18S_DvT_barplot_PENDING18S.svg"

# Fig S6 — 18S season PCoA ⚠ PENDING 18S decision
cp_fig \
    "results/18S/all/figures/18S_r1000_season_pcoa_wong.png" \
    "FigS6_18S_season_pcoa_PENDING18S.png"
cp_fig \
    "results/18S/all/figures/18S_r1000_season_pcoa_wong.svg" \
    "FigS6_18S_season_pcoa_PENDING18S.svg"

# Fig S7 — 16S seasonal PCoA (ecological seasons, corrected)
cp_fig \
    "results/16S/season/figures/pcoa_panel_Season.png" \
    "FigS7_16S_season_pcoa.png"
cp_fig \
    "results/16S/season/figures/pcoa_panel_Season.svg" \
    "FigS7_16S_season_pcoa.svg"

# Fig S8a–d — 16S seasonal alpha (4 individual metrics)
cp_fig \
    "results/16S/season/figures/alpha_faith_pd_vector_Season.png" \
    "FigS8a_16S_season_alpha_faithpd.png"
cp_fig \
    "results/16S/season/figures/alpha_faith_pd_vector_Season.svg" \
    "FigS8a_16S_season_alpha_faithpd.svg"
cp_fig \
    "results/16S/season/figures/alpha_shannon_vector_Season.png" \
    "FigS8b_16S_season_alpha_shannon.png"
cp_fig \
    "results/16S/season/figures/alpha_shannon_vector_Season.svg" \
    "FigS8b_16S_season_alpha_shannon.svg"
cp_fig \
    "results/16S/season/figures/alpha_observed_features_vector_Season.png" \
    "FigS8c_16S_season_alpha_observed.png"
cp_fig \
    "results/16S/season/figures/alpha_observed_features_vector_Season.svg" \
    "FigS8c_16S_season_alpha_observed.svg"
cp_fig \
    "results/16S/season/figures/alpha_evenness_vector_Season.png" \
    "FigS8d_16S_season_alpha_evenness.png"
cp_fig \
    "results/16S/season/figures/alpha_evenness_vector_Season.svg" \
    "FigS8d_16S_season_alpha_evenness.svg"

# Fig S9 — 16S seasonal barplot ⚠ PLACEHOLDER - needs regeneration
echo "  ⚠   FigS9_16S_season_barplot — NEEDS REGENERATION"
echo "      Current file has wrong meteorological season labels." \
    > "${OUT}/FigS9_16S_season_barplot_REGENERATE_NEEDED.txt"
echo "      Run: 10_plot_taxonomy.py with --season-column Season" \
    >> "${OUT}/FigS9_16S_season_barplot_REGENERATE_NEEDED.txt"
echo "      Expected output: Breeding / Freshwater_Nonbreeding / Saltwater labels" \
    >> "${OUT}/FigS9_16S_season_barplot_REGENERATE_NEEDED.txt"
echo "  ✅  FigS9_16S_season_barplot_REGENERATE_NEEDED.txt (placeholder)"

# Fig S10 — MiFish seasonal ecology 3-panel
cp_fig \
    "results/MiFish/figures/mifish_season_ecology_wong.png" \
    "FigS10_MiFish_season_ecology.png"
cp_fig \
    "results/MiFish/figures/mifish_season_ecology_wong.svg" \
    "FigS10_MiFish_season_ecology.svg"

# Fig S11 — MiFish habitat by season ALL groups version (supplementary)
cp_fig \
    "results/MiFish/all/figures/habitat_season/MiFish_habitat_season_all_identified_only_wong.png" \
    "FigS11_MiFish_habitat_season_allgroups.png"
cp_fig \
    "results/MiFish/all/figures/habitat_season/MiFish_habitat_season_all_identified_only_wong.svg" \
    "FigS11_MiFish_habitat_season_allgroups.svg"

echo ""
echo "── UNLABELED FILES (manual check needed) ──────"
echo "  The following files have no _wong or _purple label."
echo "  Check palette manually before including in manuscript."
for f in \
    "results/multimarker/figures/multimarker_16S_observed.png" \
    "results/multimarker/figures/multimarker_18S_observed.png" \
    "results/multimarker/figures/multimarker_cytb_observed.png" \
    "results/multimarker/figures/multimarker_MiFish_observed.png"
do
    if [[ -f "${ROOT}/$f" ]]; then
        echo "  ?  $f"
    else
        echo "  —  $f (not found)"
    fi
done

echo ""
echo "── SUMMARY ─────────────────────────────────────"
TOTAL=$(ls "${OUT}"/*.png 2>/dev/null | wc -l)
MISSING=$(ls "${OUT}"/*_MISSING.txt 2>/dev/null | wc -l || true)
echo "  PNGs copied:   ${TOTAL}"
echo "  Missing files: ${MISSING}"
echo "  Output dir:    ${OUT}"
echo ""
echo "  ⚠  Files with PENDING18S in name: hold until 18S decision"
echo "  ⚠  WRONG_SEASON_LABELS: do not use until regenerated"
echo "  ⚠  REGENERATE_NEEDED: placeholder only"
echo "  ⚠  Adeno + herpes wong files: verify these are fresh before"
echo "     submitting (purple versions were newer on 2026-03-31)"
echo "=============================================="
