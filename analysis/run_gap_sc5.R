# run_gap_sc5.R
# Fit K=9 and K=11 for sc5 to resolve whether K=10 is the true perplexity minimum.
# Merges results into existing kselection.rds and re-plots the full sorted curve.

source("D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_v2_functions.R")

qc_dir  <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"
out_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2/sc5"

gpu_state <- init_gpu()

cat("Loading QC objects...\n")
seu_list <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))
sc5      <- seu_list[["sc5"]]
counts5  <- GetAssayData(sc5, layer = "counts")

cat("\n=== SC5 — variance-based gene selection ===\n")
counts5f <- filter_genes_variance(counts5, sc5@meta.data, condition_col = "condition",
                                   n_within = 7000, n_between = 6000)
counts5t <- t(counts5f)
cat("Matrix:", nrow(counts5t), "cells x", ncol(counts5t), "genes\n")

# ---- Fit and score only the gap K values ------------------------------------
gap_ks  <- c(9, 11)
gap_res <- select_k(counts5t,
                    k_range    = gap_ks,
                    out_prefix = out_dir,
                    gpu_state  = gpu_state)

# ---- Merge with existing kselection.rds -------------------------------------
old <- readRDS(file.path(out_dir, "kselection.rds"))

all_k    <- sort(c(old$k_range, gap_ks))
old_idx  <- match(old$k_range, all_k)
gap_idx  <- match(gap_ks,      all_k)

new_perps   <- numeric(length(all_k))
new_logliks <- numeric(length(all_k))
new_perps[old_idx]   <- old$perplexities
new_logliks[old_idx] <- old$logliks
new_perps[gap_idx]   <- gap_res$perplexities
new_logliks[gap_idx] <- gap_res$logliks

ksel_updated <- list(k_range      = all_k,
                     logliks      = new_logliks,
                     perplexities = new_perps)
saveRDS(ksel_updated, file.path(out_dir, "kselection.rds"))

# ---- Print full updated curve -----------------------------------------------
cat("\n=== SC5 — FULL perplexity curve (updated) ===\n")
for (i in seq_along(all_k)) {
  flag <- if (all_k[i] %in% gap_ks) " [NEW]" else ""
  cat(sprintf("  K = %2d : perplexity = %.4f%s\n",
              all_k[i], new_perps[i], flag))
}
best_k <- all_k[which.min(new_perps)]
cat("\nMinimum perplexity at K =", best_k, "\n")

cat("\n=== SC5 gap fills done ===\n")
