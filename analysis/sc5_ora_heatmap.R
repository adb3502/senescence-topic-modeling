# sc5_ora_heatmap.R
# Heatmap of Hallmark ORA -log10(FDR) scores across sc5 topics.
# Rows: pathways significant (padj<0.05) in at least one topic.
# Fill: -log10(FDR), gray if not significant.

library(ggplot2)
library(dplyr)
library(tidyr)

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

# ---- Load ORA results -------------------------------------------------------
res <- read.csv(file.path(v2_dir, "sc5", "ora_hallmark_k10.csv"),
                stringsAsFactors = FALSE)

# Keep pathways significant in at least one topic
sig_paths <- res |> filter(p.adjust < 0.05) |> pull(Description) |> unique()
cat("Pathways significant in >= 1 topic:", length(sig_paths), "\n")

# Build full matrix including non-significant entries (all topic x pathway combos)
# First get all topic-pathway combinations from full results, fill missing as NS
all_combos <- expand.grid(
  Description = sig_paths,
  letter      = letter_order,
  stringsAsFactors = FALSE
)

res_sub <- all_combos |>
  left_join(res |> select(Description, letter, p.adjust),
            by = c("Description", "letter")) |>
  mutate(
    neg_log10_fdr = ifelse(!is.na(p.adjust) & p.adjust < 0.05,
                           -log10(pmax(p.adjust, 1e-10)), NA)
  )

# ---- Cluster rows and columns by -log10(FDR), NA=0 -------------------------
score_wide <- res_sub |>
  mutate(score = ifelse(is.na(neg_log10_fdr), 0, neg_log10_fdr)) |>
  select(Description, letter, score) |>
  pivot_wider(names_from = letter, values_from = score, values_fill = 0) |>
  tibble::column_to_rownames("Description")

score_wide <- score_wide[, letter_order[letter_order %in% colnames(score_wide)], drop = FALSE]

row_clust   <- hclust(dist(score_wide),    method = "ward.D2")
col_clust   <- hclust(dist(t(score_wide)), method = "ward.D2")
path_order  <- rownames(score_wide)[row_clust$order]
topic_order <- colnames(score_wide)[col_clust$order]

# ---- Build plot data --------------------------------------------------------
plot_df <- res_sub |>
  mutate(
    Description = factor(Description, levels = path_order),
    letter      = factor(letter,      levels = topic_order)
  )

max_score <- max(plot_df$neg_log10_fdr, na.rm = TRUE)

# ---- Heatmap ----------------------------------------------------------------
p <- ggplot(plot_df, aes(x = letter, y = Description)) +
  geom_tile(fill = "grey88", color = "white", linewidth = 0.3) +
  geom_tile(aes(fill = neg_log10_fdr), color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low      = "#FEE0D2",
    high     = "#99000D",
    limits   = c(0, max_score),
    na.value = "grey88",
    name     = "-log10(FDR)"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x       = element_text(
      colour = sapply(topic_order, function(l) consensus_cols[l]),
      face = "bold", size = 10),
    axis.text.y       = element_text(size = 7.5),
    axis.line         = element_blank(),
    axis.ticks        = element_blank(),
    legend.key.height = unit(12, "mm"),
    plot.margin       = margin(6, 6, 6, 6)
  ) +
  labs(x = "Topic", y = NULL,
       title = "sc5 Hallmark ORA — enrichment heatmap (top 50 genes per topic)",
       subtitle = "Gray = not significant (FDR\u22650.05)")

n_paths <- length(path_order)
ggsave(file.path(v2_dir, "sc5", "ora_hallmark_heatmap_k10.png"),
       p, width = 120, height = 30 + n_paths * 7, units = "mm",
       dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(v2_dir, "sc5", "ora_hallmark_heatmap_k10.svg"),
       p, width = 120, height = 30 + n_paths * 7, units = "mm",
       bg = "white", limitsize = FALSE)
cat("Saved ora_hallmark_heatmap_k10.png/svg\n")
cat("Done.\n")
