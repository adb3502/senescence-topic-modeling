# topic_correlation.R
# Correlate gene weight vectors (F matrix columns) across sc2 (K=4),
# sc78 (K=9), and sc5 (K=10) to identify conserved and sample-specific topics.
#
# Method:
#   1. Load F matrices, intersect shared genes
#   2. Spearman-correlate all topic pairs (rank-based, robust to scale differences)
#   3. Hierarchical clustering (ward.D2) + pheatmap with sample colour bar
#   4. Save heatmap PNG/SVG and correlation matrix CSV

library(pheatmap)
library(ggplot2)
library(RColorBrewer)

v2_dir  <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"
out_dir <- v2_dir   # save to top-level lda_poc_v2

# ---- Load F matrices --------------------------------------------------------
f2  <- readRDS(file.path(v2_dir, "sc2",  "fit_k4.rds"))$F
f78 <- readRDS(file.path(v2_dir, "sc78", "fit_k9.rds"))$F
f5  <- readRDS(file.path(v2_dir, "sc5",  "fit_k10.rds"))$F

cat("sc2  genes:", nrow(f2),  "  topics:", ncol(f2),  "\n")
cat("sc78 genes:", nrow(f78), "  topics:", ncol(f78), "\n")
cat("sc5  genes:", nrow(f5),  "  topics:", ncol(f5),  "\n")

# ---- Intersect shared genes -------------------------------------------------
shared_genes <- Reduce(intersect, list(rownames(f2), rownames(f78), rownames(f5)))
cat("Shared genes across all three samples:", length(shared_genes), "\n")

f2_s  <- f2[shared_genes, ]
f78_s <- f78[shared_genes, ]
f5_s  <- f5[shared_genes, ]

# ---- Name topics ------------------------------------------------------------
colnames(f2_s)  <- paste0("sc2_t",  seq_len(ncol(f2_s)))
colnames(f78_s) <- paste0("sc78_t", seq_len(ncol(f78_s)))
colnames(f5_s)  <- paste0("sc5_t",  seq_len(ncol(f5_s)))

# ---- Build combined gene x topic matrix -------------------------------------
F_combined <- cbind(f2_s, f78_s, f5_s)
cat("Combined matrix:", nrow(F_combined), "genes x", ncol(F_combined), "topics\n")

# ---- Spearman correlation matrix --------------------------------------------
# cor() on columns = correlation between topics
cor_mat <- cor(F_combined, method = "spearman")

# Save CSV
write.csv(cor_mat, file.path(out_dir, "topic_correlation_spearman.csv"))
cat("Saved topic_correlation_spearman.csv\n")

# ---- Annotation for pheatmap ------------------------------------------------
topic_names <- colnames(cor_mat)
sample_ann  <- data.frame(
  Sample = factor(
    ifelse(grepl("^sc2_",  topic_names), "sc2",
    ifelse(grepl("^sc78_", topic_names), "sc78", "sc5")),
    levels = c("sc2", "sc78", "sc5")
  ),
  row.names = topic_names
)

ann_colors <- list(
  Sample = c(sc2 = "#4477AA", sc78 = "#EE6677", sc5 = "#228833")
)

# ---- Heatmap ----------------------------------------------------------------
# Distance = 1 - r; complete linkage: only merges clusters where ALL pairs
# are strongly correlated, preventing weak cross-correlations from pulling
# distinct topic groups together.
dist_mat <- as.dist(1 - cor_mat)
hclust_obj <- hclust(dist_mat, method = "complete")

# diverging palette centred at 0; Spearman r range ~ -0.3 to 1
breaks_seq <- seq(-0.4, 1, length.out = 101)
col_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

ph <- pheatmap(
  cor_mat,
  color             = col_palette,
  breaks            = breaks_seq,
  cluster_rows      = hclust_obj,
  cluster_cols      = hclust_obj,
  annotation_col    = sample_ann,
  annotation_row    = sample_ann,
  annotation_colors = ann_colors,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize          = 8,
  fontsize_row      = 7,
  fontsize_col      = 7,
  border_color      = NA,
  cellwidth         = 14,
  cellheight        = 14,
  main              = "Topic gene-weight correlation (Spearman, complete linkage)",
  filename          = file.path(out_dir, "topic_correlation_heatmap.png"),
  width             = 9,
  height            = 8
)

# SVG version
svg(file.path(out_dir, "topic_correlation_heatmap.svg"), width = 9, height = 8)
grid::grid.newpage()
grid::grid.draw(ph$gtable)
dev.off()

cat("Saved topic_correlation_heatmap.png/svg\n")

# ---- Print top cross-sample correlations ------------------------------------
cat("\nTop cross-sample topic correlations (r >= 0.6):\n")
samples <- c("sc2", "sc78", "sc5")
for (i in seq_len(nrow(cor_mat))) {
  for (j in seq_len(ncol(cor_mat))) {
    if (j <= i) next
    si <- sub("_t\\d+$", "", rownames(cor_mat)[i])
    sj <- sub("_t\\d+$", "", colnames(cor_mat)[j])
    if (si == sj) next  # skip within-sample
    r <- cor_mat[i, j]
    if (r >= 0.6) {
      cat(sprintf("  %s  <->  %s  r = %.3f\n",
                  rownames(cor_mat)[i], colnames(cor_mat)[j], r))
    }
  }
}

cat("\nDone.\n")
