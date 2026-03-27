# gsea_enrichment_plots.R
# Individual GSEA running-sum enrichment plots for top pathways per topic.
# Re-uses the fgsea nes_list stored in bio_scores.rds — no re-running needed.
# Generates top-N enrichment curves per topic, styled with cowplot/ggplot2.
#
# Usage: set sample_name and best_k, then source or Rscript.
# Outputs: gsea_topic{i}_top{n}.png/svg per topic, plus a combined panel per sample.

library(fgsea)
library(msigdbr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(fastTopics)

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"

# ---- Reactome gene sets (needed for plotEnrichment rank vector) --------------
h_raw    <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
reactome <- split(h_raw$gene_symbol, h_raw$gs_name)
cat("Loaded", length(reactome), "Reactome gene sets\n")

# ---- Theme ------------------------------------------------------------------
theme_gsea <- function(base_size = 9) {
  theme_cowplot(base_size) +
  theme(
    axis.line   = element_line(colour = "black", linewidth = 0.4),
    axis.ticks  = element_line(colour = "black", linewidth = 0.4),
    axis.text   = element_text(colour = "black", size = base_size - 1),
    axis.title  = element_text(colour = "black", size = base_size),
    plot.title  = element_text(size = base_size, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = base_size - 1, colour = "grey40", hjust = 0),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background  = element_blank(),
    plot.margin      = margin(6, 8, 4, 8)
  )
}

# Okabe-Ito topic colours
oi_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
            "#CC79A7", "#000000", "#88CCEE", "#44AA99", "#332288", "#DDCC77")

# ---- Truncate pathway name for plot title -----------------------------------
clean_pathway <- function(x, max_chars = 55) {
  x <- sub("^REACTOME_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 1), "\u2026"), x)
}

# ---- Main function ----------------------------------------------------------
plot_gsea_sample <- function(sample_name, best_k, n_top = 5) {
  model_dir  <- file.path(v2_dir, sample_name)
  out_dir    <- file.path(model_dir, "gsea_plots")
  dir.create(out_dir, showWarnings = FALSE)

  fit  <- readRDS(file.path(model_dir, paste0("fit_k", best_k, ".rds")))
  bio  <- readRDS(file.path(model_dir, "bio_scores.rds"))
  scores_k <- bio$scores[[as.character(best_k)]]

  K           <- best_k
  topic_lbls  <- paste0("Topic ", seq_len(K))
  topic_cols  <- oi_pal[seq_len(K)]

  cat("\n===", sample_name, "K =", K, "— generating GSEA enrichment plots ===\n")

  all_topic_panels <- list()

  for (i in seq_len(K)) {
    topic_lbl <- topic_lbls[i]
    topic_col <- topic_cols[i]

    # Ranked gene vector (F matrix column = gene weight in this topic)
    rank_vec        <- fit$F[, i]
    names(rank_vec) <- rownames(fit$F)
    rank_vec        <- sort(rank_vec, decreasing = TRUE)

    # Top n_top pathways for this topic (NES > 0, sorted by padj)
    res  <- scores_k$nes_list[[i]]
    pos  <- res[!is.na(res$NES) & res$NES > 0, ]
    pos  <- pos[order(pos$padj), ]
    top_paths <- head(pos, n_top)

    if (nrow(top_paths) == 0) {
      cat("  Topic", i, ": no significant pathways (NES>0)\n")
      next
    }

    cat("  Topic", i, ": plotting", nrow(top_paths), "pathways\n")

    # One enrichment plot per pathway
    path_plots <- lapply(seq_len(nrow(top_paths)), function(j) {
      pw   <- top_paths$pathway[j]
      nes  <- round(top_paths$NES[j], 2)
      padj <- signif(top_paths$padj[j], 2)
      sz   <- top_paths$size[j]

      p <- plotEnrichment(reactome[[pw]], rank_vec) +
        labs(
          title    = clean_pathway(pw),
          subtitle = sprintf("NES = %s  |  FDR = %s  |  n = %d genes", nes, padj, sz),
          x        = "Gene rank",
          y        = "Enrichment score"
        ) +
        theme_gsea()

      # plotEnrichment layers: [1] geom_segment (ticks), [2] geom_hline, [3] geom_line (curve)
      # Recolour existing layers rather than adding a new geom_line
      for (li in seq_along(p$layers)) {
        geom_name <- class(p$layers[[li]]$geom)[1]
        if (geom_name == "GeomLine") {
          p$layers[[li]]$aes_params$colour    <- topic_col
          p$layers[[li]]$aes_params$linewidth <- 0.9
        } else if (geom_name == "GeomSegment") {
          p$layers[[li]]$aes_params$colour <- adjustcolor(topic_col, alpha.f = 0.45)
        }
      }
      p
    })

    # Arrange into a column for this topic
    topic_panel <- plot_grid(
      plotlist = path_plots,
      ncol     = 1,
      align    = "v"
    )

    # Add topic label on left
    topic_label <- ggdraw() +
      draw_label(topic_lbl, fontface = "bold", size = 10,
                 colour = topic_col, angle = 90, hjust = 0.5)

    topic_panel_labelled <- plot_grid(topic_label, topic_panel,
                                      ncol = 2, rel_widths = c(0.05, 1))

    # Save individual topic file
    fname_base <- file.path(out_dir, sprintf("gsea_%s_k%d_topic%d_top%d",
                                              sample_name, K, i, n_top))
    ggsave(paste0(fname_base, ".png"), topic_panel_labelled,
           width = 120, height = 55 * nrow(top_paths), units = "mm",
           dpi = 300, bg = "white")
    ggsave(paste0(fname_base, ".svg"), topic_panel_labelled,
           width = 120, height = 55 * nrow(top_paths), units = "mm",
           bg = "white")

    all_topic_panels[[i]] <- topic_panel_labelled
  }

  # ---- Combined panel: all topics side by side -----------------------------
  if (length(all_topic_panels) > 0) {
    combined <- plot_grid(plotlist = all_topic_panels,
                          nrow = 1, align = "h")

    fname_comb <- file.path(out_dir,
                            sprintf("gsea_%s_k%d_all_topics_top%d", sample_name, K, n_top))
    panel_w <- 120 * K
    panel_h <- 55 * n_top
    ggsave(paste0(fname_comb, ".png"), combined,
           width = panel_w, height = panel_h, units = "mm",
           dpi = 300, bg = "white", limitsize = FALSE)
    ggsave(paste0(fname_comb, ".svg"), combined,
           width = panel_w, height = panel_h, units = "mm",
           bg = "white", limitsize = FALSE)
    cat("  Combined panel saved:", basename(paste0(fname_comb, ".png")), "\n")
  }

  invisible(all_topic_panels)
}

# =============================================================================
# Run for sc2 K=4 (finalized)
# =============================================================================
plot_gsea_sample("sc2", best_k = 4, n_top = 5)

plot_gsea_sample("sc78", best_k = 9,  n_top = 5)
plot_gsea_sample("sc5",  best_k = 10, n_top = 5)

cat("\nDone.\n")
