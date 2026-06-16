# =============================================================================
# Sensitivity test: does the DvT Jaccard q=0.013 finding survive season adjustment?
#
# Runs on the rarefied 16S Jaccard distance matrix (depth=8000), with optional
# Bray-Curtis, Weighted UniFrac, and Unweighted UniFrac as supporting metrics.
#
# Three flavors:
#   1. Additive PERMANOVA — does DvT retain variance after Season is in the model?
#   2. Strata PERMANOVA   — permutations restricted within season blocks
#   3. Within-season PERMANOVA — Breeding-only and Saltwater-only DvT tests
#
# All tests use 999 permutations, seed=42, FDR within each flavor.
# =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
})

set.seed(42)

# ---- Paths (edit if needed) -------------------------------------------------
DIST_DIR   <- "./qiime2/16S/rarefied_8000/DvT/diversity/core_metrics_depth8000"
METADATA   <- "./metadata/metadata_16S_sensitivity.tsv"   # the file built from the cohort spreadsheet
OUT_PREFIX <- "sensitivity_16S_DvT"

# Distance matrices to test. Primary = Jaccard (matches the headline q=0.013).
DIST_FILES <- list(
  Jaccard            = file.path(DIST_DIR, "jaccard_distance_matrix.qza"),
  Bray_Curtis        = file.path(DIST_DIR, "bray_curtis_distance_matrix.qza"),
  Weighted_UniFrac   = file.path(DIST_DIR, "weighted_unifrac_distance_matrix.qza"),
  Unweighted_UniFrac = file.path(DIST_DIR, "unweighted_unifrac_distance_matrix.qza")
)

# ---- Helper: extract distance matrix TSV from a .qza --------------------------
# QIIME2 .qza files are zip archives. We extract the distance-matrix.tsv from each
# without requiring qiime2-R bindings (which can be painful to install).
load_qza_distance <- function(qza_path) {
  tmp <- tempfile()
  dir.create(tmp)
  utils::unzip(qza_path, exdir = tmp)
  tsv <- list.files(tmp, pattern = "distance-matrix\\.tsv$",
                    recursive = TRUE, full.names = TRUE)
  if (length(tsv) != 1) stop("Could not locate distance-matrix.tsv in ", qza_path)
  dm <- read.table(tsv, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  unlink(tmp, recursive = TRUE)
  as.dist(as.matrix(dm))
}

# ---- Load metadata -----------------------------------------------------------
meta <- read.table(METADATA, header = TRUE, sep = "\t", check.names = FALSE,
                   stringsAsFactors = FALSE)
names(meta)[1] <- "sample_id"
cat(sprintf("Loaded metadata: n=%d samples\n", nrow(meta)))

# ---- Loop over distance matrices --------------------------------------------
results <- list()

for (metric_name in names(DIST_FILES)) {
  qza <- DIST_FILES[[metric_name]]
  if (!file.exists(qza)) {
    cat(sprintf("[SKIP] %s not found at %s\n", metric_name, qza))
    next
  }
  cat(sprintf("\n========================================\n  %s\n========================================\n",
              metric_name))

  dm_full <- load_qza_distance(qza)
  dm_ids  <- labels(dm_full)

  # Normalize sample IDs to TV-prefix form. Handles three conventions:
  #   "TV230007"                   -> "TV230007"
  #   "230007"                     -> "TV230007"
  #   "TV230007-GI-16S_S1483"      -> "TV230007"
  extract_tv <- function(xs) {
    out <- regmatches(xs, regexpr("^TV\\d+", xs))
    # If no TV-prefix match found, try digits-only and prepend TV
    for (i in seq_along(xs)) {
      if (length(out) < i || nchar(out[i]) == 0) {
        m <- regmatches(xs[i], regexpr("^\\d+", xs[i]))
        out[i] <- if (length(m) > 0) paste0("TV", m) else NA_character_
      }
    }
    out
  }

  dm_tv   <- extract_tv(dm_ids)
  meta_tv <- extract_tv(meta$sample_id)

  # Match metadata to the distance matrix by normalized TV prefix
  meta_dm <- meta[match(dm_tv, meta_tv), ]
  # CRITICAL: preserve the original DM label so distance-matrix subsetting works
  meta_dm$sample_id_match <- dm_ids

  if (any(is.na(meta_dm$sample_id))) {
    missing <- dm_ids[is.na(meta_dm$sample_id)]
    cat(sprintf("  WARNING: %d samples in DM have no metadata: %s\n",
                length(missing), paste(missing, collapse = ", ")))
  }

  # ---- Build the DvT analytical subset --------------------------------------
  # Exclude Marine (matches the manuscript's primary DvT analysis).
  # Also exclude samples with Unknown season (can't be modeled as covariate).
  keep <- meta_dm$DvT %in% c("Diseased", "Trauma") &
          meta_dm$Season %in% c("Breeding", "FW_Nonbreeding", "Saltwater") &
          !is.na(meta_dm$sample_id)
  meta_sub <- meta_dm[keep, ]
  dm_sub   <- as.dist(as.matrix(dm_full)[meta_sub$sample_id_match, meta_sub$sample_id_match])

  cat(sprintf("  DvT analytical subset: n=%d (Diseased=%d, Trauma=%d)\n",
              nrow(meta_sub),
              sum(meta_sub$DvT == "Diseased"),
              sum(meta_sub$DvT == "Trauma")))
  cat("  Season distribution:\n")
  print(table(meta_sub$DvT, meta_sub$Season))

  # ---- Flavor 0: DvT alone (recapitulate headline) --------------------------
  cat("\n  [0] DvT alone:\n")
  a0 <- adonis2(dm_sub ~ DvT, data = meta_sub, permutations = 999, by = "margin")
  print(a0)

  # ---- Flavor 1: Additive DvT + Season --------------------------------------
  cat("\n  [1] Additive (DvT + Season), by='margin':\n")
  a1 <- adonis2(dm_sub ~ DvT + Season, data = meta_sub,
                permutations = 999, by = "margin")
  print(a1)

  # ---- Flavor 2: DvT with strata=Season -------------------------------------
  cat("\n  [2] DvT with permutations restricted within Season strata:\n")
  a2 <- adonis2(dm_sub ~ DvT, data = meta_sub,
                strata = meta_sub$Season, permutations = 999)
  print(a2)

  # ---- Flavor 3: Within-season DvT ------------------------------------------
  cat("\n  [3] Within-season DvT:\n")
  within_season <- list()
  for (sn in unique(meta_sub$Season)) {
    sub <- meta_sub[meta_sub$Season == sn, ]
    if (length(unique(sub$DvT)) < 2 || nrow(sub) < 6) {
      cat(sprintf("    [%s] underpowered (n=%d, levels=%d) — skipping\n",
                  sn, nrow(sub), length(unique(sub$DvT))))
      next
    }
    dm_s <- as.dist(as.matrix(dm_full)[sub$sample_id_match, sub$sample_id_match])
    cat(sprintf("    [%s] n=%d (D=%d, T=%d)\n", sn, nrow(sub),
                sum(sub$DvT == "Diseased"), sum(sub$DvT == "Trauma")))
    a_s <- adonis2(dm_s ~ DvT, data = sub, permutations = 999, by = "margin")
    print(a_s)
    within_season[[sn]] <- a_s
  }

  # ---- PERMDISP cross-check -------------------------------------------------
  cat("\n  [PERMDISP] DvT and Season homogeneity of dispersion:\n")
  pd_dvt <- permutest(betadisper(dm_sub, meta_sub$DvT), permutations = 999)
  pd_sea <- permutest(betadisper(dm_sub, meta_sub$Season), permutations = 999)
  cat("  PERMDISP DvT:\n");    print(pd_dvt)
  cat("  PERMDISP Season:\n"); print(pd_sea)

  results[[metric_name]] <- list(
    DvT_alone       = a0,
    additive        = a1,
    strata          = a2,
    within_season   = within_season,
    permdisp_DvT    = pd_dvt,
    permdisp_Season = pd_sea
  )
}

# ---- Save a tidy summary -----------------------------------------------------
saveRDS(results, paste0(OUT_PREFIX, "_results.rds"))

# Build summary table for the headline metric (Jaccard)
if (!is.null(results$Jaccard)) {
  r <- results$Jaccard
  summary_tbl <- data.frame(
    Test = c("DvT alone",
             "DvT (in DvT+Season, marginal)",
             "Season (in DvT+Season, marginal)",
             "DvT with strata=Season"),
    R2   = c(r$DvT_alone$R2[1], r$additive$R2[1], r$additive$R2[2], r$strata$R2[1]),
    F    = c(r$DvT_alone$F[1],  r$additive$F[1],  r$additive$F[2],  r$strata$F[1]),
    p    = c(r$DvT_alone$`Pr(>F)`[1], r$additive$`Pr(>F)`[1],
             r$additive$`Pr(>F)`[2],  r$strata$`Pr(>F)`[1])
  )
  cat("\n\n=========================================================\n")
  cat("HEADLINE SUMMARY — Jaccard, DvT+Trauma subset, Marine excluded\n")
  cat("=========================================================\n")
  print(summary_tbl, row.names = FALSE)
  write.table(summary_tbl, paste0(OUT_PREFIX, "_jaccard_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("\nDone. Saved:\n",
    "  ", paste0(OUT_PREFIX, "_results.rds"), "\n",
    "  ", paste0(OUT_PREFIX, "_jaccard_summary.tsv"), "\n", sep = "")
