source("D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_v2_functions.R")

qc_dir  <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"
out_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2/sc78"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gpu_state <- init_gpu()

cat("Loading QC objects...\n")
seu_list <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))
sc78     <- seu_list[["sc78"]]
counts78 <- GetAssayData(sc78, layer = "counts")

cat("\n=== SC78 — variance-based gene selection ===\n")
counts78f <- filter_genes_variance(counts78, sc78@meta.data, condition_col = "condition",
                                    n_within = 7000, n_between = 5000)
counts78t <- t(counts78f)
cat("Matrix:", nrow(counts78t), "cells x", ncol(counts78t), "genes\n")

sc78_ks <- select_k(counts78t,
                    k_range    = c(3, 4, 5, 6, 8, 10, 12, 15),
                    out_prefix = out_dir,
                    gpu_state  = gpu_state)

best_k <- c(3, 4, 5, 6, 8, 10, 12, 15)[which.min(sc78_ks$perplexities)]
cat("\nPreliminary best K (min perplexity):", best_k, "\n")
print_top_genes(sc78_ks$fits[[as.character(best_k)]], "sc78")

cat("\n=== SC78 done. Outputs in:", out_dir, "===\n")
