source("D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_v2_functions.R")

qc_dir  <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"
out_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2/sc5"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gpu_state <- init_gpu()

cat("Loading QC objects...\n")
seu_list <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))
sc5      <- seu_list[["sc5"]]
counts5  <- GetAssayData(sc5, layer = "counts")

# Print condition distribution for the timecourse sample
cat("\nSC5 condition distribution:\n")
print(table(sc5@meta.data$condition))

cat("\n=== SC5 (timecourse) — variance-based gene selection ===\n")
# n_between is larger here: 7 timepoints means rich between-condition variance
counts5f <- filter_genes_variance(counts5, sc5@meta.data, condition_col = "condition",
                                   n_within = 7000, n_between = 6000)
counts5t <- t(counts5f)
cat("Matrix:", nrow(counts5t), "cells x", ncol(counts5t), "genes\n")

# Timecourse may warrant more topics — extend upper K range
sc5_ks <- select_k(counts5t,
                   k_range    = c(3, 4, 5, 6, 8, 10, 12, 15),
                   out_prefix = out_dir,
                   gpu_state  = gpu_state)

best_k <- c(3, 4, 5, 6, 8, 10, 12, 15)[which.min(sc5_ks$perplexities)]
cat("\nPreliminary best K (min perplexity):", best_k, "\n")
print_top_genes(sc5_ks$fits[[as.character(best_k)]], "sc5")

cat("\n=== SC5 done. Outputs in:", out_dir, "===\n")
