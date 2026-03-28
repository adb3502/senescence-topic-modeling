# ora_hallmark.R
# Over-representation analysis (ORA) using MSigDB Hallmark gene sets
# on top 50 genes per topic for sc2 (K=4), sc78 (K=9), sc5 (K=10).
# Uses Fisher's exact test via clusterProfiler::enricher.
# Outputs: ora_hallmark_k{K}.csv and dot plot per sample.

library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(dplyr)

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"

# ---- Consensus letter maps --------------------------------------------------
sc2_map  <- c("Topic 1"="F","Topic 2"="H","Topic 3"="A","Topic 4"="B")
sc78_map <- c("Topic 1"="G","Topic 2"="D","Topic 3"="E","Topic 4"="F",
              "Topic 5"="C","Topic 6"="H","Topic 7"="B","Topic 8"="I","Topic 9"="A")
sc5_map  <- c("Topic 1"="A","Topic 2"="B","Topic 3"="J","Topic 4"="F",
              "Topic 5"="G","Topic 6"="D","Topic 7"="K","Topic 8"="L",
              "Topic 9"="C","Topic 10"="E")

# Consensus colors
consensus_cols <- c(
  H="#D55E00", A="#0072B2", B="#009E73", F="#E69F00",
  D="#56B4E9", E="#CC79A7", C="#44AA99", G="#882255",
  I="#DDCC77", K="#AA4499", J="#332288", L="#999933"
)

# ---- Hallmark gene sets (human) --------------------------------------------
cat("Loading Hallmark gene sets...\n")
h_df   <- msigdbr(species = "Homo sapiens", collection = "H")
h_t2g  <- h_df[, c("gs_name", "gene_symbol")]
# clean names: strip HALLMARK_ prefix
h_t2g$gs_name <- sub("HALLMARK_", "", h_t2g$gs_name)
h_t2g$gs_name <- gsub("_", " ", h_t2g$gs_name)
universe <- unique(h_df$gene_symbol)
cat("Hallmark sets:", length(unique(h_t2g$gs_name)), "\n")

# ---- ORA function -----------------------------------------------------------
run_ora <- function(fit, topic_map, n_top = 50) {
  K      <- ncol(fit$F)
  t_lbls <- paste0("Topic ", seq_len(K))

  results <- lapply(seq_len(K), function(i) {
    tl     <- t_lbls[i]
    letter <- topic_map[tl]
    genes  <- names(sort(fit$F[, i], decreasing = TRUE))[1:n_top]

    res <- tryCatch(
      enricher(gene      = genes,
               TERM2GENE = h_t2g,
               universe  = universe,
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.20,
               minGSSize = 5,
               maxGSSize = 500),
      error = function(e) NULL
    )

    if (is.null(res) || nrow(as.data.frame(res)) == 0) {
      cat(sprintf("  %s (%s): no significant pathways\n", letter, tl))
      return(NULL)
    }

    df <- as.data.frame(res)
    df$letter <- letter
    df$topic  <- tl
    cat(sprintf("  %s (%s): %d significant pathways\n", letter, tl, nrow(df)))
    df
  })

  bind_rows(results)
}

# ---- Run per sample ---------------------------------------------------------
samples <- list(
  sc2  = list(fit = readRDS(file.path(v2_dir, "sc2",  "fit_k4.rds")),  map = sc2_map,  K = 4),
  sc78 = list(fit = readRDS(file.path(v2_dir, "sc78", "fit_k9.rds")),  map = sc78_map, K = 9),
  sc5  = list(fit = readRDS(file.path(v2_dir, "sc5",  "fit_k10.rds")), map = sc5_map,  K = 10)
)

for (sname in names(samples)) {
  cat("\n===", sname, "===\n")
  s   <- samples[[sname]]
  out <- run_ora(s$fit, s$map, n_top = 50)

  if (is.null(out) || nrow(out) == 0) {
    cat("  No results for", sname, "\n")
    next
  }

  out_dir <- file.path(v2_dir, sname)
  write.csv(out, file.path(out_dir, paste0("ora_hallmark_k", s$K, ".csv")),
            row.names = FALSE)
  cat("  Saved ora_hallmark_k", s$K, ".csv\n", sep = "")

  # ---- Dot plot -------------------------------------------------------------
  # Top 5 pathways per topic, deduplicated for display
  plot_df <- out |>
    group_by(letter) |>
    slice_min(p.adjust, n = 5) |>
    ungroup() |>
    mutate(
      Description = factor(Description),
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]; as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      neg_log10_padj = -log10(p.adjust),
      letter = factor(letter, levels = names(consensus_cols)[names(consensus_cols) %in% unique(out$letter)])
    )

  # Order pathways by letter then padj
  path_order <- plot_df |>
    arrange(letter, p.adjust) |>
    pull(Description) |>
    unique() |>
    as.character()
  plot_df$Description <- factor(plot_df$Description, levels = rev(path_order))

  p <- ggplot(plot_df, aes(x = letter, y = Description,
                            size = GeneRatio_num, color = neg_log10_padj)) +
    geom_point() +
    scale_color_gradient(low = "#DEEBF7", high = "#08306B",
                         name = "-log10(FDR)") +
    scale_size_continuous(name = "Gene ratio", range = c(2, 7)) +
    theme_classic(base_size = 9) +
    theme(axis.text.x  = element_text(colour = sapply(levels(plot_df$letter),
                                                        function(l) consensus_cols[l]),
                                       face = "bold", size = 10),
          axis.text.y  = element_text(size = 7),
          legend.position = "right",
          plot.margin  = margin(6, 6, 6, 6)) +
    labs(x = "Topic", y = NULL,
         title = paste0(sname, " — Hallmark ORA (top 5 per topic, FDR<0.05)"))

  n_paths <- length(levels(plot_df$Description))
  n_topics <- length(unique(plot_df$letter))
  ggsave(file.path(out_dir, paste0("ora_hallmark_dotplot_k", s$K, ".png")),
         p, width = 55 + n_topics * 22, height = 40 + n_paths * 5.5,
         units = "mm", dpi = 300, bg = "white", limitsize = FALSE)
  cat("  Saved ora_hallmark_dotplot_k", s$K, ".png\n", sep = "")
}

cat("\nDone.\n")
