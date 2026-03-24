# 02_normalize_integrate.R
# Log-normalize -> HVG selection -> PCA -> Harmony batch correction -> UMAP
# Batch variable: sample (sc2, sc5, sc6, sc78, sc9)
# Timecourse replicates (sc5, sc6, sc9) and controls (sc2, sc78) integrated together

library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)

qc_dir  <- "D:/Projects/Topic Modeling/analysis/qc_output"
out_dir <- "D:/Projects/Topic Modeling/analysis/integration_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Theme ------------------------------------------------------------------
theme_publication <- function(base_size = 11) {
  theme_cowplot(base_size) +
  theme(
    axis.line         = element_line(colour = "black", linewidth = 0.5),
    axis.ticks        = element_line(colour = "black", linewidth = 0.5),
    axis.text         = element_text(colour = "black", size = base_size),
    axis.title        = element_text(colour = "black", size = base_size + 1),
    strip.background  = element_blank(),
    strip.text        = element_text(colour = "black", size = base_size, face = "bold"),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.background   = element_blank(),
    plot.title        = element_blank(),
    plot.margin       = margin(8, 8, 8, 8)
  )
}

sample_colors <- c(sc2 = "#E64B35", sc5 = "#4DBBD5", sc6 = "#00A087",
                   sc78 = "#F39B7F", sc9 = "#3C5488")

tp_order  <- c("Control", "Sen 4hrs", "Sen D1", "Sen D2", "Sen D3", "Sen D4", "Sen D7")
tp_colors <- setNames(viridis::viridis(7, option = "D"), tp_order)

# ---- Load QC objects --------------------------------------------------------
cat("Loading QC objects...\n")
seu_list <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))

# ---- Merge all samples ------------------------------------------------------
cat("Merging samples...\n")
merged <- merge(
  seu_list[[1]],
  y         = seu_list[-1],
  add.cell.ids = names(seu_list),
  project   = "GSE223128"
)
cat("Merged:", ncol(merged), "cells,", nrow(merged), "genes\n")

# Ensure layers are split by sample for per-sample normalization
merged <- JoinLayers(merged)

# ---- Normalize --------------------------------------------------------------
cat("Normalizing...\n")
merged <- NormalizeData(merged, normalization.method = "LogNormalize",
                        scale.factor = 1e4, verbose = FALSE)

# ---- Variable features (timecourse only, to avoid control-driven HVGs) -----
cat("Finding variable features...\n")
# Use all cells but find HVGs — 3000 is robust for SPLIT-seq
merged <- FindVariableFeatures(merged, selection.method = "vst",
                                nfeatures = 3000, verbose = FALSE)

# Remove MT and ribosomal genes from HVG list (technical, not biological)
hvg <- VariableFeatures(merged)
hvg <- hvg[!grepl("^MT-|^RP[SL]", hvg)]
VariableFeatures(merged) <- hvg
cat("HVGs after MT/ribo removal:", length(hvg), "\n")

# ---- Scale and PCA ----------------------------------------------------------
cat("Scaling and running PCA...\n")
merged <- ScaleData(merged, features = hvg, verbose = FALSE)
merged <- RunPCA(merged, features = hvg, npcs = 50, verbose = FALSE)

# Elbow plot to choose PCs
p_elbow <- ElbowPlot(merged, ndims = 50) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "red", alpha = 0.6) +
  labs(x = "PC", y = "Std. dev.") +
  theme_publication()

ggsave(file.path(out_dir, "fig_elbow.svg"), p_elbow,
       width = 85, height = 70, units = "mm", bg = "white")
ggsave(file.path(out_dir, "fig_elbow.png"), p_elbow,
       width = 85, height = 70, units = "mm", dpi = 300, bg = "white")

# ---- Harmony integration (Seurat v5 API) ------------------------------------
cat("Running Harmony (batch = sample)...\n")
# Split layers by sample before integration (required for IntegrateLayers)
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$sample)
merged <- IntegrateLayers(
  object        = merged,
  method        = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = FALSE
)
# Rejoin layers after integration
merged <- JoinLayers(merged)

# ---- UMAP -------------------------------------------------------------------
cat("Running UMAP...\n")
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30,
                  reduction.name = "umap.harmony", verbose = FALSE)

# ---- Clustering (for QC/sanity check) ---------------------------------------
cat("Finding neighbors and clusters...\n")
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.3, verbose = FALSE)

cat("Clusters found:", length(unique(merged$seurat_clusters)), "\n")

# ---- UMAP figures -----------------------------------------------------------
cat("Generating UMAP plots...\n")

umap_df <- as.data.frame(Embeddings(merged, "umap.harmony"))
umap_df$sample    <- merged$sample
umap_df$role      <- merged$role
umap_df$condition <- merged$condition
umap_df$cluster   <- merged$seurat_clusters
umap_df$nCount    <- merged$nCount_RNA
umap_df$pct_mt    <- merged$pct_mt

# By sample
p_sample <- ggplot(umap_df[sample(nrow(umap_df)), ],
                   aes(x = umap_1, y = umap_2, color = sample)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_manual(values = sample_colors) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Sample") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_publication()

# By condition (timecourse only)
tc_umap <- umap_df[umap_df$role == "timecourse" &
                     umap_df$condition %in% tp_order, ]
tc_umap$condition <- factor(tc_umap$condition, levels = tp_order)

p_timepoint <- ggplot(tc_umap[sample(nrow(tc_umap)), ],
                      aes(x = umap_1, y = umap_2, color = condition)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_manual(values = tp_colors) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Timepoint") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_publication()

# By cluster
p_cluster <- ggplot(umap_df[sample(nrow(umap_df)), ],
                    aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_viridis_d(option = "H") +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_publication()

# Combined panel
fig_umap <- (p_sample | p_timepoint | p_cluster) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(file.path(out_dir, "fig_umap_overview.svg"), fig_umap,
       width = 250, height = 85, units = "mm", bg = "white")
ggsave(file.path(out_dir, "fig_umap_overview.png"), fig_umap,
       width = 250, height = 85, units = "mm", dpi = 300, bg = "white")

# Timecourse split by timepoint
p_tp_facet <- ggplot(tc_umap[sample(nrow(tc_umap)), ],
                     aes(x = umap_1, y = umap_2, color = condition)) +
  geom_point(size = 0.2, alpha = 0.5) +
  scale_color_manual(values = tp_colors) +
  facet_wrap(~condition, nrow = 2) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "fig_umap_timepoints_facet.svg"), p_tp_facet,
       width = 174, height = 100, units = "mm", bg = "white")
ggsave(file.path(out_dir, "fig_umap_timepoints_facet.png"), p_tp_facet,
       width = 174, height = 100, units = "mm", dpi = 300, bg = "white")

# ---- Save -------------------------------------------------------------------
cat("\nSaving integrated object...\n")
saveRDS(merged, file.path(out_dir, "seu_integrated.rds"))
cat("Saved to:", file.path(out_dir, "seu_integrated.rds"), "\n")
cat("\nDone.\n")
