# biological_coherence_v2.R
# Bio coherence scoring for v2 models (variance-based gene selection)
# Reactome GSEA on F matrix columns; coherence x specificity to select K
# Run after identifying perplexity plateau from K-selection logs.
#
# Usage: set plateau ranges below per sample, then source or Rscript.

library(fastTopics)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(cowplot)
library(dplyr)

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"

# ---- Reactome gene sets ------------------------------------------------------
h_raw    <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
reactome <- split(h_raw$gene_symbol, h_raw$gs_name)
cat("Loaded", length(reactome), "Reactome gene sets\n")

# ---- Scoring -----------------------------------------------------------------
score_model <- function(fit, gene_sets, nperm = 1000, fdr_thr = 0.1) {
  k <- ncol(fit$F)
  cat("    Scoring K =", k, "...\n")

  nes_list <- lapply(seq_len(k), function(i) {
    rank_vec        <- fit$F[, i]
    names(rank_vec) <- rownames(fit$F)
    rank_vec        <- sort(rank_vec, decreasing = TRUE)
    fgsea(pathways = gene_sets, stats = rank_vec,
          nperm = nperm, minSize = 10, maxSize = 500, nproc = 1)
  })

  # Coherence: mean -log10(best padj per topic)
  best_padj <- sapply(nes_list, function(res) {
    pos <- res[res$NES > 0, ]
    if (nrow(pos) == 0) return(1)
    min(pos$padj, na.rm = TRUE)
  })
  coherence <- mean(-log10(pmax(best_padj, 1e-10)))

  # Specificity: fraction of topics with unique top pathway
  top_paths <- sapply(nes_list, function(res) {
    pos <- res[res$NES > 0, ]
    if (nrow(pos) == 0) return(NA_character_)
    pos$pathway[which.min(pos$padj)]
  })
  n_unique    <- length(unique(na.omit(top_paths)))
  specificity <- n_unique / k
  combined    <- coherence * specificity

  cat("      coherence =", round(coherence, 3),
      " specificity =", round(specificity, 3),
      " combined =", round(combined, 3), "\n")
  cat("      top pathways:", paste(top_paths, collapse = " | "), "\n")

  list(coherence = coherence, specificity = specificity, combined = combined,
       top_pathways = top_paths, nes_list = nes_list)
}

score_k_range <- function(sample_name, plateau, nperm = 1000) {
  model_dir <- file.path(v2_dir, sample_name)
  cat("\n===", sample_name, "— K =", paste(plateau, collapse = ", "), "===\n")

  fits <- lapply(setNames(plateau, plateau), function(k) {
    ckpt <- file.path(model_dir, paste0("fit_k", k, ".rds"))
    if (file.exists(ckpt)) readRDS(ckpt) else { cat("  Missing K=", k, "\n"); NULL }
  })

  scores <- lapply(as.character(plateau), function(ks) {
    if (is.null(fits[[ks]])) return(NULL)
    score_model(fits[[ks]], reactome, nperm = nperm)
  })
  names(scores) <- as.character(plateau)

  df <- data.frame(
    K           = plateau,
    coherence   = sapply(scores, function(s) if (is.null(s)) NA else s$coherence),
    specificity = sapply(scores, function(s) if (is.null(s)) NA else s$specificity),
    combined    = sapply(scores, function(s) if (is.null(s)) NA else s$combined)
  )

  cat("\nSummary:\n"); print(df)
  best_k <- plateau[which.max(df$combined)]
  cat("Best K by coherence x specificity:", best_k, "\n")

  theme_pub <- function(base_size = 11) {
    theme_cowplot(base_size) +
    theme(axis.line  = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text  = element_text(colour = "black", size = base_size),
          axis.title = element_text(colour = "black", size = base_size + 1),
          panel.background = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.background  = element_blank(),
          plot.title       = element_blank(), plot.margin      = margin(8, 8, 8, 8))
  }

  p_comb <- ggplot(df, aes(x = K, y = combined)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
    scale_x_continuous(breaks = plateau) +
    labs(x = "Number of topics (K)", y = "Coherence \u00d7 specificity") + theme_pub()
  p_spec <- ggplot(df, aes(x = K, y = specificity)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
    scale_x_continuous(breaks = plateau) +
    labs(x = "Number of topics (K)", y = "Topic specificity") + theme_pub()
  p_coh <- ggplot(df, aes(x = K, y = coherence)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
    scale_x_continuous(breaks = plateau) +
    labs(x = "Number of topics (K)", y = "Biological coherence") + theme_pub()

  ggsave(file.path(model_dir, "bio_combined.png"),    p_comb, width = 100, height = 75, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(model_dir, "bio_combined.svg"),    p_comb, width = 100, height = 75, units = "mm", bg = "white")
  ggsave(file.path(model_dir, "bio_specificity.png"), p_spec, width = 100, height = 75, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(model_dir, "bio_coherence.png"),   p_coh,  width = 100, height = 75, units = "mm", dpi = 300, bg = "white")

  saveRDS(list(scores = scores, df = df, best_k = best_k),
          file.path(model_dir, "bio_scores.rds"))

  list(scores = scores, df = df, best_k = best_k)
}

# =============================================================================
# SET PLATEAU RANGES HERE after reading perplexity curves
# (fill in once K-selection runs complete)
# =============================================================================

# sc2 is already finalized at K=4; bio_scores.rds intact — do not re-run
plateau_sc2  <- NULL
# sc78 gap fills (K=7, K=9) confirm plateau at K=7–9:
#   K=7=4.487, K=8≈4.481, K=9=4.481, K=10=4.560 (jump)
plateau_sc78 <- c(7, 8, 9)
# sc5 gap fills (K=9, K=11) confirm plateau at K=9–11:
#   K=9=2.363, K=10=2.362, K=11=2.360, K=12=2.440 (jump)
plateau_sc5  <- c(9, 10, 11)

if (!is.null(plateau_sc2))  bio_sc2  <- score_k_range("sc2",  plateau_sc2)
if (!is.null(plateau_sc78)) bio_sc78 <- score_k_range("sc78", plateau_sc78)
if (!is.null(plateau_sc5))  bio_sc5  <- score_k_range("sc5",  plateau_sc5)

cat("\nDone.\n")
