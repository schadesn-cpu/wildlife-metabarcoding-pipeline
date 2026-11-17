#!/usr/bin/env Rscript
# pooltestr_qpcr_TEST.R  — TEST, pending Joe/Lawrence sign-off
suppressMessages({library(PoolTestR); library(dplyr); library(readr)})

d <- read_tsv("results/pooltestr_TEST/qpcr_input.tsv", show_col_types=FALSE)
pathogens <- c("a_phagocytophilum","b_microti","b_burgdorferi","b_mayonii","b_miyamotoi","e_muris","powassan_virus")

cat("Total qPCR samples:", nrow(d), "\n")
cat("Single ticks (pool=1):", sum(d$pool_size==1), " | Pools (>1):", sum(d$pool_size>1), "\n\n")

cat("=== direct single-tick prevalence vs PoolTestR MLE (all data) ===\n")
for (p in pathogens) {
  if (!p %in% names(d)) next
  singles <- d[d$pool_size==1, ]
  sp <- 100*mean(singles[[p]], na.rm=TRUE); npos <- sum(singles[[p]], na.rm=TRUE)
  dd <- d[!is.na(d[[p]]), ]; dd$result <- dd[[p]]
  mle <- tryCatch(PoolPrev(dd, result, poolSize=pool_size, bayesian=FALSE),
                  error=function(e) {cat("ERR",p,conditionMessage(e),"\n"); NULL})
  if (!is.null(mle)) {
    cat(sprintf("%-20s single=%.1f%% (npos=%d)   ", p, sp, npos))
    print(as.data.frame(mle))
  }
}

cat("\n=== by host (MLE) ===\n")
for (p in c("a_phagocytophilum","b_burgdorferi","b_microti","b_miyamotoi")) {
  if (!p %in% names(d)) next
  dd <- d[!is.na(d[[p]]), ]; dd$result <- dd[[p]]
  mle <- tryCatch(PoolPrev(dd, result, poolSize=pool_size, hostgrp, bayesian=FALSE),
                  error=function(e) {cat("ERR",p,conditionMessage(e),"\n"); NULL})
  cat("\n--", p, "--\n")
  if (!is.null(mle)) print(as.data.frame(mle))
}
