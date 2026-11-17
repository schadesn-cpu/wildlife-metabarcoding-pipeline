#!/usr/bin/env Rscript
# pooltestr_prevalence_16S.R
# Pool-aware pathogen prevalence from the VERIFIED 16S analysis table, using PoolTestR MLE.
# Uses the SAME validated call as pooltestr_qpcr_TEST.R: PoolPrev(dd, result, PoolSize=pool_size).
#
# Input: results/analysis/analysis_table_16S.tsv  (built by build_analysis_table.py)
#   one row per qualifying sample; columns include NumTicks (pool size), Host_group,
#   and per-pathogen presence columns (Anaplasma_pos, Borreliella_pos, ... = 0/1, "any read" rule).
# Real pathogens only (Rickettsia/Francisella are endosymbionts, excluded upstream).
#
# NOTE: still pending Joe/Lawrence sign-off on PoolTestR-for-metabarcoding appropriateness.
#       This produces the estimates; framing (primary vs supporting) is a separate decision.

suppressMessages({library(PoolTestR); library(dplyr); library(readr)})

inp <- "results/analysis/analysis_table_16S.tsv"
d <- read_tsv(inp, show_col_types=FALSE)

# pathogen presence columns (the *_pos 0/1 columns from the analysis table)
pathogens <- c("Anaplasma","Borreliella","Borrelia","Bartonella","Ehrlichia")
poscols   <- paste0(pathogens, "_pos")
poscols   <- poscols[poscols %in% names(d)]
pathogens <- sub("_pos$","",poscols)

cat("Analysis table samples:", nrow(d), "\n")
cat("Single ticks (NumTicks==1):", sum(d$NumTicks==1, na.rm=TRUE),
    " | Pools (>1):", sum(d$NumTicks>1, na.rm=TRUE), "\n")
cat("Pool size range:", min(d$NumTicks,na.rm=TRUE), "-", max(d$NumTicks,na.rm=TRUE), "\n\n")

# ---- 1. OVERALL pool-aware prevalence per pathogen (all qualifying samples) ----
cat("=== OVERALL pool-aware prevalence (PoolTestR MLE, all samples) ===\n")
cat(sprintf("%-14s %10s %10s %24s\n","pathogen","n_pos_pools","crude_pos%","MLE_prev%_[95%CI]"))
overall <- list()
for (i in seq_along(poscols)) {
  p <- pathogens[i]; pc <- poscols[i]
  dd <- d[!is.na(d[[pc]]) & !is.na(d$NumTicks), ]
  dd$result <- dd[[pc]]
  crude <- 100*mean(dd$result)
  npos  <- sum(dd$result)
  mle <- tryCatch(PoolPrev(dd, result, poolSize=NumTicks, bayesian=FALSE), error=function(e) NULL)
  if (!is.null(mle)) {
    est <- 100*mle$PrevMLE[1]; lo <- 100*mle$CILow[1]; hi <- 100*mle$CIHigh[1]
    cat(sprintf("%-14s %10d %9.1f   %6.2f [%.2f-%.2f]\n", p, npos, crude, est, lo, hi))
    overall[[p]] <- data.frame(pathogen=p, n_pos=npos, crude=crude, MLE=est, lo=lo, hi=hi)
  } else {
    cat(sprintf("%-14s %10d %9.1f   MLE failed\n", p, npos, crude))
  }
}
overall_df <- bind_rows(overall)
write_tsv(overall_df, "results/analysis/prevalence_overall_16S.tsv")

# ---- 2. HOST-STRATIFIED prevalence (by Host_group) ----
# Only groups with enough samples give stable estimates; PoolTestR handles the grouping.
cat("\n=== HOST-STRATIFIED prevalence (PoolTestR MLE by Host_group) ===\n")
host_results <- list()
for (i in seq_along(poscols)) {
  p <- pathogens[i]; pc <- poscols[i]
  dd <- d[!is.na(d[[pc]]) & !is.na(d$NumTicks) & !is.na(d$Host_group), ]
  dd$result <- dd[[pc]]
  cat("\n--", p, "--\n")
  mle <- tryCatch(PoolPrev(dd, result, poolSize=NumTicks, Host_group), error=function(e) NULL)
  if (!is.null(mle)) {
    out <- mle[, c("Host_group","PrevMLE","CILow","CIHigh","NumberOfPools")]
    out$PrevMLE <- round(100*out$PrevMLE,2); out$CILow <- round(100*out$CILow,2)
    out$CIHigh <- round(100*out$CIHigh,2)
    print(out, row.names=FALSE)
    out$pathogen <- p
    host_results[[p]] <- out
  } else cat("  (host-stratified MLE failed — likely too few pools in a group)\n")
}
if (length(host_results)) {
  write_tsv(bind_rows(host_results), "results/analysis/prevalence_by_host_16S.tsv")
}

cat("\nwrote results/analysis/prevalence_overall_16S.tsv")
cat("\nwrote results/analysis/prevalence_by_host_16S.tsv\n")
cat("\nSANITY CHECK: single-tick-only prevalence (should be near the locked all-host figures\n")
cat("if computed on the same population; here it's the analysis-table single ticks only):\n")
st <- d[d$NumTicks==1, ]
for (i in seq_along(poscols)) {
  cat(sprintf("  %-14s %.1f%%\n", pathogens[i], 100*mean(st[[poscols[i]]], na.rm=TRUE)))
}
