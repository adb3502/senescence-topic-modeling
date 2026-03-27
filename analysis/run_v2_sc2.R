source("D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_v2_functions.R")

qc_dir  <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"
out_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2/sc2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gpu_state <- init_gpu()

cat("Loading QC objects...\n")
seu_list <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))
sc2      <- seu_list[["sc2"]]
counts2  <- GetAssayData(sc2, layer = "counts")

cat("\n=== SC2 — variance-based gene selection ===\n")
counts2f <- filter_genes_variance(counts2, sc2@meta.data, condition_col = "condition",
                                   n_within = 7000, n_between = 5000)
counts2t <- t(counts2f)
cat("Matrix:", nrow(counts2t), "cells x", ncol(counts2t), "genes\n")

sc2_ks <- select_k(counts2t,
                   k_range    = c(3, 4, 5, 6, 8, 10, 12, 15),
                   out_prefix = out_dir,
                   gpu_state  = gpu_state)

best_k <- c(3, 4, 5, 6, 8, 10, 12, 15)[which.min(sc2_ks$perplexities)]
cat("\nPreliminary best K (min perplexity):", best_k, "\n")
print_top_genes(sc2_ks$fits[[as.character(best_k)]], "sc2")

cat("\n=== SC2 done. Outputs in:", out_dir, "===\n")
