#!/usr/bin/env Rscript
# decontam on the 16S table, loon-informed:
#   * negatives = vetted blanks + NTCs + XB-OR-NTC + Batch-12 tagged blanks,
#     EXCLUDING the 4 mislabeled heavy XBs and all positive controls.
#   * runs prevalence method at BOTH threshold 0.1 (default) and 0.5 (stringent).
#   * FLAGS, does not remove. Per-ASV inspection shows read distribution in
#     negatives vs real samples so we apply the loon lesson: an ASV abundant in
#     real ticks and sparse in blanks is a FALSE POSITIVE (keep, e.g. Rickettsia
#     endosymbiont / Ixodes), one abundant in blanks is a TRUE contaminant.
#   * also reports per-negative total reads, to catch any additional heavy
#     (mislabeled) blank that would poison the negative set.
#
# Run:  cd ~/meedlab/tick_project && Rscript scripts/decontam_16S.R
# Needs: /tmp/16s_feature_table.tsv (biom convert of qiime2/16S/all/dada2/table.qza)

suppressMessages(library(decontam))

TABLE <- "/tmp/16s_feature_table.tsv"
# ASV -> Taxon label (for the keep/remove judgement). Swap to the region-matched
# export if you prefer; any ASV->name map works for identifying flagged ASVs.
TAX   <- "/tmp/16s_tax/taxonomy.tsv"
OUT   <- "results/decontam/16S"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MISLABELED <- c("XB1-2", "XBG12", "XB-A1", "XB-D8")   # heavy, held real template

# --- read feature table (features x samples; skip biom comment line) ---
ft <- read.delim(TABLE, skip = 1, row.names = 1, check.names = FALSE)
samples <- colnames(ft)
prefix  <- sub("\\..*$", "", samples)   # token before first "."

is_ctrl <- grepl("XB|NTC|POS|PAC|PC-WELL|BLANK|EXTRACTION", samples, ignore.case = TRUE)
is_pos  <- grepl("POS|PC-WELL|PAC", samples, ignore.case = TRUE)
mislab  <- prefix %in% MISLABELED

is_neg  <- is_ctrl & !is_pos & !mislab   # vetted negatives

# PATCH: drop the 3 confirmed-contaminated NTCs (Rickettsia-laden, ~18k-68k reads)
HOT_CONTROLS <- c("NTC.250609_Batch-11","NTC-2.250609_Batch-11","NTC-H11.250327_Batch-08")
is_neg <- is_neg & !(samples %in% HOT_CONTROLS)
is_true <- !is_ctrl                      # real biological samples
keep    <- is_true | is_neg              # drop positives + mislabeled entirely

cat("=== sample accounting ===\n")
cat("total samples in table :", length(samples), "\n")
cat("vetted negatives       :", sum(is_neg), "\n")
cat("true biological samples:", sum(is_true), "\n")
cat("dropped (pos+mislabeled):", sum(!keep), "\n\n")

cat("negatives used:\n"); print(sort(samples[is_neg]))
cat("\ndropped as positive/mislabeled:\n"); print(sort(samples[!keep]))

# --- guard: per-negative total reads (catch a hidden heavy blank) ---
neg_names <- samples[is_neg]
neg_tot   <- sort(colSums(ft[, neg_names, drop = FALSE]), decreasing = TRUE)
cat("\n=== per-negative total 16S reads (watch for an outlier = possible mislabel) ===\n")
print(neg_tot)
med <- median(neg_tot)
hot <- neg_tot[neg_tot > 20 * pmax(med, 1)]
if (length(hot) > 0) {
  cat("\n** WARNING: these negatives are >20x the median; inspect before trusting **\n")
  print(hot)
}

# --- subset table, build neg vector, run decontam ---
ft2  <- ft[, keep, drop = FALSE]
neg2 <- is_neg[keep]
seqtab <- t(as.matrix(ft2))            # decontam wants samples x features

run_thr <- function(thr) {
  r <- isContaminant(seqtab, neg = neg2, method = "prevalence", threshold = thr)
  r$ASV <- rownames(r); r
}
r01 <- run_thr(0.1)
r05 <- run_thr(0.5)
f01 <- r01$ASV[which(r01$contaminant)]
f05 <- r05$ASV[which(r05$contaminant)]

cat("\n=== decontam results ===\n")
cat("flagged at 0.1 (default) :", length(f01), "\n")
cat("flagged at 0.5 (stringent):", length(f05), "\n")
cat("overlap                  :", length(intersect(f01, f05)), "\n")

# --- inspection: read distribution neg vs samples for every flagged ASV ---
flagged  <- union(f01, f05)
neg_cols <- which(neg2); samp_cols <- which(!neg2)
if (length(flagged) > 0) {
  sub <- ft2[flagged, , drop = FALSE]
  insp <- data.frame(ASV = flagged, stringsAsFactors = FALSE)
  insp$reads_neg  <- rowSums(sub[, neg_cols,  drop = FALSE])
  insp$reads_samp <- rowSums(sub[, samp_cols, drop = FALSE])
  insp$prev_neg   <- rowSums(sub[, neg_cols,  drop = FALSE] > 0)
  insp$prev_samp  <- rowSums(sub[, samp_cols, drop = FALSE] > 0)
  insp$mean_neg   <- insp$reads_neg  / length(neg_cols)
  insp$mean_samp  <- insp$reads_samp / length(samp_cols)
  insp$flag01 <- insp$ASV %in% f01
  insp$flag05 <- insp$ASV %in% f05
  insp$p      <- r01[insp$ASV, "p"]
  insp$hint <- ifelse(insp$mean_samp > 10 * pmax(insp$mean_neg, 1e-9),
                      "KEEP? abundant in real samples",
               ifelse(insp$mean_neg >= insp$mean_samp,
                      "CONTAMINANT? higher in negatives", "review"))
  if (file.exists(TAX)) {
    tax <- read.delim(TAX, row.names = 1, check.names = FALSE)
    insp$Taxon <- tax[insp$ASV, "Taxon"]
  }
  insp <- insp[order(-insp$reads_samp), ]
  write.table(insp, file.path(OUT, "inspection.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\n=== flagged ASV inspection (sorted by reads in real samples) ===\n")
  show <- insp[, c("ASV","reads_neg","reads_samp","prev_neg","prev_samp","flag01","flag05","hint")]
  show$ASV <- substr(show$ASV, 1, 8)
  print(show, row.names = FALSE)
  if ("Taxon" %in% names(insp)) {
    cat("\n=== taxa of flagged ASVs ===\n")
    for (i in seq_len(nrow(insp)))
      cat(sprintf("  %s  reads_samp=%d reads_neg=%d  %s\n",
          substr(insp$ASV[i],1,8), insp$reads_samp[i], insp$reads_neg[i],
          ifelse(is.na(insp$Taxon[i]), "(no taxon)", insp$Taxon[i])))
  }
} else {
  cat("\nNo ASVs flagged at either threshold. Data is clean by decontam prevalence.\n")
}

write.table(r01, file.path(OUT, "decontam_p01_full.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(r05, file.path(OUT, "decontam_p05_full.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nwrote", file.path(OUT, "inspection.tsv"), "and full decontam tables.\n")
cat("NOTHING REMOVED. Decide keep/remove per ASV from the inspection above.\n")
