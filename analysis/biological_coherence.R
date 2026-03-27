# biological_coherence.R
# For a fitted set of topic models (one per K), compute biological coherence
# using fgsea on MSigDB Hallmarks gene weights (F matrix columns).
#
# Coherence score per model:
#   - Per topic: mean NES of significantly enriched pathways (FDR < 0.1)
#   - Per model: mean across topics (coherence) * topic specificity
#
# Specificity penalises K values where topics are redundant (similar pathway profiles).
# Run ONLY over the perplexity plateau range — passed in as k_range.

library(fastTopics)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(cowplot)
library(dplyr)

# ---- Gene sets: MSigDB Hallmarks (human) ------------------------------------
h_raw     <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
hallmarks <- split(h_raw$gene_symbol, h_raw$gs_name)

cat("Loaded", length(hallmarks), "Reactome gene sets\n")

# ---- Scoring functions -------------------------------------------------------

# Score a single topic: mean NES of significant pathways (FDR < 0.1)
# Returns 0 if no pathways significant
score_topic <- function(fit, topic_idx, gene_sets, nperm = 1000, fdr_thr = 0.1) {
  rank_vec        <- fit$F[, topic_idx]
  names(rank_vec) <- rownames(fit$F)
  rank_vec        <- sort(rank_vec, decreasing = TRUE)

  res <- fgsea(pathways = gene_sets, stats = rank_vec,
               nperm = nperm, minSize = 10, maxSize = 500,
               nproc = 1)

  sig <- res[res$padj < fdr_thr & res$NES > 0, ]
  if (nrow(sig) == 0) return(0)
  mean(sig$NES)
}

# Specificity: 1 - mean pairwise cosine similarity of NES vectors across topics
# High specificity = topics have distinct pathway profiles
topic_specificity <- function(fit, gene_sets, nperm = 1000) {
  k <- ncol(fit$F)
  if (k == 1) return(1)

  # NES matrix: pathways x topics
  nes_mat <- sapply(seq_len(k), function(i) {
    rank_vec        <- fit$F[, i]
    names(rank_vec) <- rownames(fit$F)
    rank_vec        <- sort(rank_vec, decreasing = TRUE)
    res <- fgsea(pathways = gene_sets, stats = rank_vec,
                 nperm = nperm, minSize = 10, maxSize = 500,
                 nproc = 1)
    # Return NES vector ordered by pathway name
    nes <- setNames(res$NES, res$pathway)
    nes[names(gene_sets)]
  })
  nes_mat[is.na(nes_mat)] <- 0

  # Pairwise cosine similarity
  norms    <- sqrt(colSums(nes_mat^2))
  norms    <- pmax(norms, 1e-10)
  nes_norm <- sweep(nes_mat, 2, norms, "/")
  cos_mat  <- t(nes_norm) %*% nes_norm

  # Mean of off-diagonal entries
  idx      <- lower.tri(cos_mat)
  mean_cos <- mean(cos_mat[idx])
  1 - mean_cos
}

# Score a full model at K
# coherence  : mean -log10(best padj) per topic — how strongly each topic maps to any pathway
# specificity: fraction of topics with a unique top pathway — penalises redundant splits
# combined   : coherence * specificity
score_model <- function(fit, gene_sets, nperm = 1000, fdr_thr = 0.1) {
  k <- ncol(fit$F)
  cat("    Scoring K =", k, "topics...\n")

  nes_list <- lapply(seq_len(k), function(i) {
    rank_vec        <- fit$F[, i]
    names(rank_vec) <- rownames(fit$F)
    rank_vec        <- sort(rank_vec, decreasing = TRUE)
    fgsea(pathways = gene_sets, stats = rank_vec,
          nperm = nperm, minSize = 10, maxSize = 500,
          nproc = 1)
  })

  # Coherence: mean -log10(best padj) across topics
  best_padj <- sapply(nes_list, function(res) {
    pos <- res[res$NES > 0, ]
    if (nrow(pos) == 0) return(1)
    min(pos$padj, na.rm = TRUE)
  })
  coherence <- mean(-log10(pmax(best_padj, 1e-10)))

  # Specificity: fraction of topics whose top pathway is unique across topics
  top_pathways <- sapply(nes_list, function(res) {
    pos <- res[res$NES > 0, ]
    if (nrow(pos) == 0) return(NA_character_)
    pos$pathway[which.min(pos$padj)]
  })
  n_unique    <- length(unique(na.omit(top_pathways)))
  specificity <- n_unique / k

  combined <- coherence * specificity
  cat("      coherence =", round(coherence, 3),
      " specificity =", round(specificity, 3),
      " (", n_unique, "unique top pathways /", k, "topics )",
      " combined =", round(combined, 3), "\n")
  cat("      top pathways:", paste(top_pathways, collapse = " | "), "\n")

  list(coherence    = coherence,
       specificity  = specificity,
       combined     = combined,
       top_pathways = top_pathways,
       nes_list     = nes_list)
}

# ---- Main: score across plateau K range -------------------------------------
score_k_range <- function(fits_list, k_range, out_prefix,
                           gene_sets = hallmarks, nperm = 1000) {
  cat("Scoring biological coherence for K =",
      paste(k_range, collapse = ", "), "...\n")

  scores <- lapply(as.character(k_range), function(ks) {
    fit <- fits_list[[ks]]
    if (is.null(fit)) {
      cat("  K =", ks, "not found in fits_list, skipping\n")
      return(NULL)
    }
    score_model(fit, gene_sets, nperm = nperm)
  })
  names(scores) <- as.character(k_range)

  df <- data.frame(
    K           = k_range,
    coherence   = sapply(scores, `[[`, "coherence"),
    specificity = sapply(scores, `[[`, "specificity"),
    combined    = sapply(scores, `[[`, "combined")
  )

  cat("\nSummary:\n")
  print(df)

  theme_pub <- function(base_size = 11) {
    theme_cowplot(base_size) +
    theme(axis.line        = element_line(colour = "black", linewidth = 0.5),
          axis.ticks       = element_line(colour = "black", linewidth = 0.5),
          axis.text        = element_text(colour = "black", size = base_size),
          axis.title       = element_text(colour = "black", size = base_size + 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background  = element_blank(),
          plot.title       = element_blank(),
          plot.margin      = margin(8, 8, 8, 8))
  }

  p_coh <- ggplot(df, aes(x = K, y = coherence)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of topics (K)", y = "Biological coherence") +
    theme_pub()

  p_spec <- ggplot(df, aes(x = K, y = specificity)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of topics (K)", y = "Topic specificity") +
    theme_pub()

  p_comb <- ggplot(df, aes(x = K, y = combined)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of topics (K)", y = "Coherence × specificity") +
    theme_pub()

  ggsave(paste0(out_prefix, "_bio_coherence.png"),  p_coh,  width = 90, height = 70, units = "mm", dpi = 300, bg = "white")
  ggsave(paste0(out_prefix, "_bio_specificity.png"), p_spec, width = 90, height = 70, units = "mm", dpi = 300, bg = "white")
  ggsave(paste0(out_prefix, "_bio_combined.png"),   p_comb, width = 90, height = 70, units = "mm", dpi = 300, bg = "white")

  best_k <- k_range[which.max(df$combined)]
  cat("\nBest K by coherence x specificity:", best_k, "\n")

  list(scores = scores, df = df, best_k = best_k)
}

# =============================================================================
# SC2 — plateau K=4:10 from perplexity curve
# =============================================================================
lda_dir   <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc"
model_dir <- file.path(lda_dir, "models")
plot_dir  <- file.path(lda_dir, "k_selection")

cat("\n=== SC2 biological coherence (K=4 to K=10) ===\n")
plateau_sc2 <- c(4, 5, 6, 7, 8, 10)

sc2_fits <- lapply(setNames(plateau_sc2, plateau_sc2), function(k) {
  ckpt <- file.path(model_dir, paste0("sc2_fit_k", k, ".rds"))
  if (file.exists(ckpt)) readRDS(ckpt) else NULL
})

sc2_bio <- score_k_range(sc2_fits, k_range = plateau_sc2,
                          out_prefix = file.path(plot_dir, "sc2"))

saveRDS(sc2_bio, file.path(model_dir, "sc2_bio_scores.rds"))
cat("Saved sc2 bio scores\n")
