# 03_lda_poc.R — LDA topic model POC on sc2 and sc78 independently
# fastTopics (Poisson NMF = multinomial LDA at MLE)
# Raw counts, corpus-style gene filtering, K selection by held-out likelihood
# No timepoint information used — topics emerge unsupervised, conditions checked post hoc

library(Seurat)
library(fastTopics)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)

qc_dir  <- "D:/Projects/Topic Modeling/analysis/qc_output"
out_dir <- "D:/Projects/Topic Modeling/analysis/lda_poc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_publication <- function(base_size = 11) {
  theme_cowplot(base_size) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.text = element_text(colour = "black", size = base_size),
    axis.title = element_text(colour = "black", size = base_size + 1),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = base_size, face = "bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    plot.title = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )
}

# ---- Gene filtering (corpus analogy) ----------------------------------------
# Remove stop words (ubiquitous genes) and rare words (sparse genes)
# Keep genes detected in >= 1% and <= 95% of cells
filter_genes <- function(counts, min_pct = 0.01, max_pct = 0.95) {
  n_cells    <- ncol(counts)
  detection  <- Matrix::rowMeans(counts > 0)
  keep       <- detection >= min_pct & detection <= max_pct
  # Also remove MT and ribosomal genes (stop words)
  keep       <- keep & !grepl("^MT-|^RP[SL]", rownames(counts))
  cat("  Genes before filter:", nrow(counts), "\n")
  cat("  Genes after filter: ", sum(keep), "\n")
  counts[keep, ]
}

# ---- K selection via held-out likelihood ------------------------------------
select_k <- function(counts_t, k_range = c(4, 6, 8, 10, 12), out_prefix) {
  cat("  Fitting models for K =", paste(k_range, collapse = ", "), "...\n")

  fits <- list()
  for (k in k_range) {
    ckpt <- paste0(out_prefix, "_fit_k", k, ".rds")
    if (file.exists(ckpt)) {
      cat("    K =", k, "... loading from checkpoint\n")
      fits[[as.character(k)]] <- readRDS(ckpt)
    } else {
      cat("    K =", k, "... ")
      set.seed(42)
      fit <- fit_topic_model(counts_t, k = k, numiter.main = 100, verbose = "none")
      saveRDS(fit, ckpt)
      fits[[as.character(k)]] <- fit
      cat("done\n")
    }
  }

  # Compare log-likelihoods
  logliks <- sapply(fits, function(f) mean(loglik_multinom_topic_model(counts_t, f)))
  cat("  Log-likelihoods:\n")
  for (i in seq_along(k_range)) {
    cat("    K =", k_range[i], ":", round(logliks[i], 2), "\n")
  }

  # Plot
  ll_df <- data.frame(K = k_range, loglik = logliks)
  p_ll  <- ggplot(ll_df, aes(x = K, y = loglik)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of topics (K)", y = "Log-likelihood") +
    theme_publication()
  ggsave(paste0(out_prefix, "_loglik.svg"), p_ll,
         width = 85, height = 70, units = "mm", bg = "white")
  ggsave(paste0(out_prefix, "_loglik.png"), p_ll,
         width = 85, height = 70, units = "mm", dpi = 300, bg = "white")

  list(fits = fits, logliks = logliks)
}

# ---- Topic visualization helpers --------------------------------------------
plot_topic_structure <- function(fit, meta, condition_col, out_prefix, sample_name) {
  # Structure plot — per-cell topic proportions stacked bar
  # Order cells by dominant topic within each condition
  theta <- poisson2multinom(fit)$L  # cells x topics
  df    <- as.data.frame(theta)
  df$cell      <- rownames(df)
  df$condition <- meta[[condition_col]][match(df$cell, rownames(meta))]

  # Long format
  topic_cols <- grep("^k", names(df), value = TRUE)
  df_long    <- tidyr::pivot_longer(df, cols = all_of(topic_cols),
                                     names_to = "topic", values_to = "weight")

  # Order cells by dominant topic
  dom_topic  <- apply(theta, 1, which.max)
  cell_order <- rownames(theta)[order(df$condition, dom_topic)]
  df_long$cell <- factor(df_long$cell, levels = cell_order)

  p <- ggplot(df_long, aes(x = cell, y = weight, fill = topic)) +
    geom_col(width = 1) +
    scale_fill_viridis_d(option = "H") +
    facet_grid(~condition, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Topic weight", fill = "Topic") +
    theme_publication() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0.5, "mm"))

  ggsave(paste0(out_prefix, "_structure.svg"), p,
         width = 174, height = 70, units = "mm", bg = "white")
  ggsave(paste0(out_prefix, "_structure.png"), p,
         width = 174, height = 70, units = "mm", dpi = 300, bg = "white")
  invisible(theta)
}

plot_top_genes <- function(fit, counts_t, k, out_prefix) {
  # Top genes per topic by largest F (topic-specific gene frequency)
  de  <- de_analysis(fit, counts_t, verbose = "none")
  top <- lapply(seq_len(k), function(i) {
    topic_de <- de$postmean[, i]
    names(topic_de) <- colnames(counts_t)
    head(sort(topic_de, decreasing = TRUE), 10)
  })

  cat("  Top genes per topic:\n")
  for (i in seq_len(k)) {
    cat("    Topic", i, ":", paste(names(top[[i]]), collapse = ", "), "\n")
  }

  # Heatmap of top genes
  top_genes <- unique(unlist(lapply(top, names)))
  F_mat     <- fit$F[top_genes, , drop = FALSE]  # genes x topics

  df_heat <- as.data.frame(as.matrix(F_mat)) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "topic", values_to = "freq")

  p <- ggplot(df_heat, aes(x = topic, y = gene, fill = freq)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma") +
    labs(x = "Topic", y = NULL, fill = "Gene freq.") +
    theme_publication() +
    theme(axis.text.y = element_text(size = 7))

  ggsave(paste0(out_prefix, "_top_genes_heatmap.svg"), p,
         width = 100, height = 120, units = "mm", bg = "white")
  ggsave(paste0(out_prefix, "_top_genes_heatmap.png"), p,
         width = 100, height = 120, units = "mm", dpi = 300, bg = "white")
  invisible(de)
}

# =============================================================================
# Load data
# =============================================================================
cat("Loading QC objects...\n")
seu_list <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))

# =============================================================================
# SC2 — Nutlin, OxStress, Quiescence (3,313 cells)
# =============================================================================
cat("\n=== SC2 ===\n")
sc2      <- seu_list[["sc2"]]
counts2  <- GetAssayData(sc2, layer = "counts")
counts2f <- filter_genes(counts2)
counts2t <- t(counts2f)   # cells x genes (fastTopics convention)

cat("Running K selection for sc2...\n")
sc2_ks <- select_k(counts2t, k_range = c(3, 4, 5, 6, 8),
                   out_prefix = file.path(out_dir, "sc2"))

# Best K — pick the elbow visually; also fit best model with more iterations
best_k_sc2 <- 5
ckpt_sc2 <- file.path(out_dir, "fit_sc2_final.rds")
if (file.exists(ckpt_sc2)) {
  cat("\nLoading sc2 final model from checkpoint...\n")
  fit_sc2 <- readRDS(ckpt_sc2)
} else {
  cat("\nFitting final sc2 model K =", best_k_sc2, "with 200 iterations...\n")
  set.seed(42)
  fit_sc2 <- fit_topic_model(counts2t, k = best_k_sc2,
                              numiter.main = 200, verbose = "none")
  saveRDS(fit_sc2, ckpt_sc2)
}

# Visualize
cat("Plotting sc2 topic structure...\n")
theta_sc2 <- plot_topic_structure(fit_sc2, sc2@meta.data, "condition",
                                   file.path(out_dir, "sc2"), "sc2")
cat("Plotting sc2 top genes...\n")
de_sc2 <- plot_top_genes(fit_sc2, counts2t, best_k_sc2,
                          file.path(out_dir, "sc2"))

# =============================================================================
# SC78 — Control, Nutlin, OxStress, Quiescence x 2 reps (18,429 cells)
# =============================================================================
cat("\n=== SC78 ===\n")
sc78      <- seu_list[["sc78"]]
counts78  <- GetAssayData(sc78, layer = "counts")
counts78f <- filter_genes(counts78)
counts78t <- t(counts78f)

cat("Running K selection for sc78...\n")
sc78_ks <- select_k(counts78t, k_range = c(3, 4, 5, 6, 8),
                    out_prefix = file.path(out_dir, "sc78"))

ckpt_sc78 <- file.path(out_dir, "fit_sc78_final.rds")
if (file.exists(ckpt_sc78)) {
  cat("\nLoading sc78 final model from checkpoint...\n")
  fit_sc78 <- readRDS(ckpt_sc78)
} else {
  cat("\nFitting final sc78 model K =", best_k_sc78, "with 200 iterations...\n")
  set.seed(42)
  fit_sc78 <- fit_topic_model(counts78t, k = best_k_sc78,
                               numiter.main = 200, verbose = "none")
  saveRDS(fit_sc78, ckpt_sc78)
}

cat("Plotting sc78 topic structure...\n")
theta_sc78 <- plot_topic_structure(fit_sc78, sc78@meta.data, "condition",
                                    file.path(out_dir, "sc78"), "sc78")
cat("Plotting sc78 top genes...\n")
de_sc78 <- plot_top_genes(fit_sc78, counts78t, best_k_sc78,
                           file.path(out_dir, "sc78"))

cat("\n=== Done. Outputs in:", out_dir, "===\n")
