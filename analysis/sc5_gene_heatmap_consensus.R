# sc5_gene_heatmap_consensus.R
# Top 10 genes per sc5 topic, expression heatmap across cells.
# Topics labeled by consensus letters (A-L), both genes and topics clustered.
# Uses raw normalized log-expression from Seurat object.

library(Seurat)
library(fastTopics)
library(pheatmap)
library(dplyr)
library(RColorBrewer)

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"
qc_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"

# ---- Consensus letter map for sc5 ------------------------------------------
sc5_map <- c("Topic 1"="A", "Topic 2"="B", "Topic 3"="J", "Topic 4"="F",
             "Topic 5"="G", "Topic 6"="D", "Topic 7"="K", "Topic 8"="L",
             "Topic 9"="C", "Topic 10"="E")

# Letter order (original topic order for consistency)
letter_order <- unname(sc5_map)   # A B J F G D K L C E

# Consensus colors
consensus_cols <- c(
  H="#D55E00", A="#0072B2", B="#009E73", F="#E69F00",
  D="#56B4E9", E="#CC79A7", C="#44AA99", G="#882255",
  I="#DDCC77", K="#AA4499", J="#332288", L="#999933"
)
topic_cols <- consensus_cols[letter_order]

# ---- Load data --------------------------------------------------------------
cat("Loading data...\n")
seu  <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))[["sc5"]]
fit5 <- readRDS(file.path(v2_dir, "sc5", "fit_k10.rds"))

K      <- 10L
t_lbls <- paste0("Topic ", seq_len(K))

# ---- Top 10 genes per topic (by F weight) -----------------------------------
top_genes_per_topic <- lapply(seq_len(K), function(i) {
  names(sort(fit5$F[, i], decreasing = TRUE))[1:10]
})
names(top_genes_per_topic) <- letter_order

# Union of all top genes (deduplicated)
all_top_genes <- unique(unlist(top_genes_per_topic))
cat("Unique top-10 genes across all topics:", length(all_top_genes), "\n")

# ---- Gene annotation: which topic is each gene top-weighted in? -------------
gene_topic <- sapply(all_top_genes, function(g) {
  weights <- fit5$F[g, ]
  letter_order[which.max(weights)]
})
gene_ann <- data.frame(
  Topic = factor(gene_topic, levels = letter_order),
  row.names = all_top_genes
)

# Order genes by their primary topic (letter_order), then by weight within topic
gene_order <- unlist(lapply(letter_order, function(lt) {
  i <- which(names(sc5_map) == names(sc5_map)[sc5_map == lt])
  genes_i <- top_genes_per_topic[[lt]]
  # keep only genes whose primary topic is this letter
  genes_i[gene_topic[genes_i] == lt]
}))
# append any genes whose primary is elsewhere but appear in another topic's top10
remaining <- setdiff(all_top_genes, gene_order)
gene_order <- c(gene_order, remaining)

# ---- Pseudobulk: mean log-normalised expression per topic-dominant cell group
# Use L matrix to assign each cell to its dominant topic
L_norm <- fit5$L / rowSums(fit5$L)
colnames(L_norm) <- letter_order
dominant_topic <- letter_order[max.col(L_norm)]
names(dominant_topic) <- rownames(L_norm)

# Normalise Seurat object if not already
if (!"data" %in% names(seu@assays$RNA)) {
  seu <- NormalizeData(seu, verbose = FALSE)
}

# Subset to top genes present in Seurat
genes_present <- intersect(gene_order, rownames(seu))
cat("Genes present in Seurat object:", length(genes_present), "\n")

expr_mat <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "data")[genes_present, ])

# Mean expression per dominant-topic group
topic_expr <- sapply(letter_order, function(lt) {
  cells <- names(dominant_topic)[dominant_topic == lt]
  cells <- intersect(cells, colnames(expr_mat))
  if (length(cells) == 0) return(rep(NA, nrow(expr_mat)))
  rowMeans(expr_mat[, cells, drop = FALSE])
})
colnames(topic_expr) <- letter_order
rownames(topic_expr) <- genes_present

cat("Expression matrix:", nrow(topic_expr), "genes x", ncol(topic_expr), "topics\n")

# ---- Z-score across topics (per gene) for visual contrast -------------------
topic_expr_z <- t(scale(t(topic_expr)))
topic_expr_z[is.nan(topic_expr_z)] <- 0

# ---- Annotations ------------------------------------------------------------
col_ann <- data.frame(
  Topic = factor(letter_order, levels = letter_order),
  row.names = letter_order
)
row_ann <- data.frame(
  Topic = factor(gene_topic[genes_present], levels = letter_order),
  row.names = genes_present
)
ann_colors <- list(
  Topic = setNames(as.character(topic_cols), letter_order)
)

# ---- Plot -------------------------------------------------------------------
breaks_z <- seq(-2.5, 2.5, length.out = 101)
col_pal   <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

out_dir <- file.path(v2_dir, "sc5")

ph <- pheatmap(
  topic_expr_z[genes_present, ],
  color             = col_pal,
  breaks            = breaks_z,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  clustering_method = "ward.D2",
  annotation_col    = col_ann,
  annotation_row    = row_ann,
  annotation_colors = ann_colors,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize          = 8,
  fontsize_row      = 7,
  fontsize_col      = 9,
  border_color      = NA,
  cellwidth         = 18,
  cellheight        = 8,
  main              = "sc5 top-10 genes per topic (mean log-norm expr, z-scored across topics)",
  filename          = file.path(out_dir, "gene_expr_heatmap_consensus.png"),
  width             = 10,
  height            = 14
)

# SVG
svg(file.path(out_dir, "gene_expr_heatmap_consensus.svg"), width = 10, height = 14)
grid::grid.newpage()
grid::grid.draw(ph$gtable)
dev.off()

cat("Saved gene_expr_heatmap_consensus.png/svg\n")
cat("Done.\n")
