# sc5_hallmark_heatmap.R
# Heatmap of Hallmark GSEA NES scores across sc5 topics.
# Rows: pathways significant (padj<0.05) in at least one topic.
# Fill: NES (red=up, blue=down), gray if padj>=0.05.

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

# ---- Load GSEA results ------------------------------------------------------
res <- read.csv(file.path(v2_dir, "sc5", "gsea_hallmark_k10.csv"),
                stringsAsFactors = FALSE)

# Keep pathways significant in at least one topic
sig_paths <- res |>
  filter(padj < 0.05) |>
  pull(pathway) |>
  unique()
cat("Pathways significant in >= 1 topic:", length(sig_paths), "\n")

res_sub <- res |>
  filter(pathway %in% sig_paths) |>
  mutate(
    letter = factor(letter, levels = letter_order),
    # NES_masked: NA (shown as gray) if not significant
    NES_plot = ifelse(padj < 0.05, NES, NA)
  )

# ---- Cluster pathways by NES profile ----------------------------------------
# Build full NES matrix (all topics x sig pathways), NAs as 0 for clustering
nes_wide <- res_sub |>
  select(pathway, letter, NES) |>
  pivot_wider(names_from = letter, values_from = NES, values_fill = 0) |>
  tibble::column_to_rownames("pathway")

# Ensure column order
nes_wide <- nes_wide[, letter_order[letter_order %in% colnames(nes_wide)], drop = FALSE]

row_clust  <- hclust(dist(nes_wide),   method = "ward.D2")
col_clust  <- hclust(dist(t(nes_wide)), method = "ward.D2")
path_order  <- rownames(nes_wide)[row_clust$order]
topic_order <- colnames(nes_wide)[col_clust$order]

# ---- Build plot data --------------------------------------------------------
plot_df <- res_sub |>
  mutate(pathway = factor(pathway, levels = path_order),
         letter  = factor(letter,  levels = topic_order))

# symmetric NES limits
nes_lim <- max(abs(res_sub$NES), na.rm = TRUE)
nes_lim <- ceiling(nes_lim * 10) / 10

# ---- Heatmap ----------------------------------------------------------------
p <- ggplot(plot_df, aes(x = letter, y = pathway)) +
  # gray tile for all cells first
  geom_tile(fill = "grey88", color = "white", linewidth = 0.3) +
  # colored tile only where significant
  geom_tile(aes(fill = NES_plot), color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low      = "#FEE0D2",
    high     = "#99000D",
    limits   = c(1, nes_lim),
    na.value = "grey88",
    name     = "NES"
  ) +
  # topic color bar on x-axis text
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x  = element_text(
      colour = sapply(topic_order, function(l) consensus_cols[l]),
      face   = "bold", size = 10),
    axis.text.y  = element_text(size = 7.5),
    axis.line    = element_blank(),
    axis.ticks   = element_blank(),
    legend.key.height = unit(12, "mm"),
    plot.margin  = margin(6, 6, 6, 6)
  ) +
  labs(x = "Topic", y = NULL,
       title = "sc5 Hallmark GSEA — NES heatmap",
       subtitle = "Gray = not significant (FDR\u22650.05)  |  Color intensity = enrichment strength (NES)")

n_paths <- length(path_order)
ggsave(file.path(v2_dir, "sc5", "gsea_hallmark_heatmap_k10.png"),
       p, width = 120, height = 30 + n_paths * 7, units = "mm",
       dpi = 300, bg = "white", limitsize = FALSE)
cat("Saved gsea_hallmark_heatmap_k10.png\n")

# SVG
ggsave(file.path(v2_dir, "sc5", "gsea_hallmark_heatmap_k10.svg"),
       p, width = 120, height = 30 + n_paths * 7, units = "mm",
       bg = "white", limitsize = FALSE)
cat("Saved gsea_hallmark_heatmap_k10.svg\n")

cat("Done.\n")
