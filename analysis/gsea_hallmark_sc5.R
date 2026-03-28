# gsea_hallmark_sc5.R
# GSEA using MSigDB Hallmark gene sets on sc5 topics (K=10).
# Ranked by F matrix column weights (all genes, not just top N).
# Outputs: gsea_hallmark_k10.csv and dot plot.

library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(fastTopics)

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"

sc5_map <- c("Topic 1"="A","Topic 2"="B","Topic 3"="J","Topic 4"="F",
             "Topic 5"="G","Topic 6"="D","Topic 7"="K","Topic 8"="L",
             "Topic 9"="C","Topic 10"="E")
letter_order <- unname(sc5_map)

consensus_cols <- c(
  H="#D55E00", A="#0072B2", B="#009E73", F="#E69F00",
  D="#56B4E9", E="#CC79A7", C="#44AA99", G="#882255",
  I="#DDCC77", K="#AA4499", J="#332288", L="#999933"
)

# ---- Hallmark gene sets -----------------------------------------------------
cat("Loading Hallmark gene sets...\n")
h_df  <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(h_df$gene_symbol, h_df$gs_name)
hallmark <- lapply(hallmark, unique)
# clean names
names(hallmark) <- sub("HALLMARK_", "", names(hallmark))
names(hallmark) <- gsub("_", " ", names(hallmark))
cat("Hallmark sets:", length(hallmark), "\n")

# ---- Load fit ---------------------------------------------------------------
fit5 <- readRDS(file.path(v2_dir, "sc5", "fit_k10.rds"))
K    <- 10L

# ---- Run fgsea per topic ----------------------------------------------------
cat("Running GSEA...\n")
set.seed(42)

all_res <- lapply(seq_len(K), function(i) {
  tl     <- paste0("Topic ", i)
  letter <- sc5_map[tl]

  rank_vec        <- fit5$F[, i]
  names(rank_vec) <- rownames(fit5$F)
  rank_vec        <- sort(rank_vec, decreasing = TRUE)

  res <- fgsea(pathways   = hallmark,
               stats      = rank_vec,
               minSize    = 10,
               maxSize    = 500,
               scoreType  = "pos",
               nPermSimple = 10000)

  res$letter <- letter
  res$topic  <- tl
  cat(sprintf("  %s (%s): %d significant (padj<0.05)\n",
              letter, tl, sum(res$padj < 0.05, na.rm = TRUE)))
  as.data.frame(res)
})

res_df <- bind_rows(all_res)
# fgsea leadingEdge is a list column — flatten to semicolon-separated string
res_df$leadingEdge <- sapply(res_df$leadingEdge, function(x)
  if (is.null(x)) "" else paste(x, collapse = ";"))

out_dir <- file.path(v2_dir, "sc5")
write.csv(res_df, file.path(out_dir, "gsea_hallmark_k10.csv"), row.names = FALSE)
cat("Saved gsea_hallmark_k10.csv\n")

# ---- Dot plot: NES, significant only ----------------------------------------
sig_df <- res_df |>
  filter(padj < 0.05) |>
  mutate(
    letter = factor(letter, levels = letter_order),
    neg_log10_padj = -log10(pmax(padj, 1e-10))
  )

# Select top 5 per topic by padj, keep unique pathways
top_df <- sig_df |>
  group_by(letter) |>
  slice_min(padj, n = 5) |>
  ungroup()

path_order <- top_df |>
  arrange(letter, NES) |>
  pull(pathway) |>
  unique()
top_df$pathway <- factor(top_df$pathway, levels = rev(path_order))

p <- ggplot(top_df, aes(x = letter, y = pathway, size = abs(NES), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                        midpoint = 0, name = "NES") +
  scale_size_continuous(name = "|NES|", range = c(2, 8)) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x  = element_text(
      colour = sapply(levels(top_df$letter), function(l) consensus_cols[l]),
      face = "bold", size = 10),
    axis.text.y  = element_text(size = 7.5),
    legend.position = "right",
    plot.margin  = margin(6, 6, 6, 6)
  ) +
  labs(x = "Topic", y = NULL,
       title = "sc5 — Hallmark GSEA (top 5 per topic, FDR<0.05)")

n_paths  <- length(levels(top_df$pathway))
n_topics <- length(unique(as.character(top_df$letter)))
ggsave(file.path(out_dir, "gsea_hallmark_dotplot_k10.png"),
       p, width = 55 + n_topics * 22, height = 40 + n_paths * 5.5,
       units = "mm", dpi = 300, bg = "white", limitsize = FALSE)
cat("Saved gsea_hallmark_dotplot_k10.png\n")

# ---- Also print full significant results per topic --------------------------
cat("\n--- Significant Hallmark pathways per topic (padj<0.05) ---\n")
for (lt in letter_order) {
  sub <- sig_df |> filter(letter == lt) |> arrange(padj)
  if (nrow(sub) == 0) { cat(sprintf("\n%s: none\n", lt)); next }
  cat(sprintf("\n%s:\n", lt))
  for (j in seq_len(nrow(sub))) {
    cat(sprintf("  %-45s NES=%+.2f  padj=%.3g\n",
                sub$pathway[j], sub$NES[j], sub$padj[j]))
  }
}

cat("\nDone.\n")
