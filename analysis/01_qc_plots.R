# 01_qc_plots.R — Publication-quality QC figures for GSE223128
# Follows publication-r-viz skill: cowplot theme, viridis/Nature palettes, SVG export

library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggridges)
library(dplyr)
library(tidyr)

out_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"
seu_list <- readRDS(file.path(out_dir, "seu_list_qc.rds"))

# ---- Theme ------------------------------------------------------------------
theme_publication <- function(base_size = 11) {
  theme_cowplot(base_size) +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.5),
    axis.ticks       = element_line(colour = "black", linewidth = 0.5),
    axis.text        = element_text(colour = "black", size = base_size),
    axis.title       = element_text(colour = "black", size = base_size + 1),
    strip.background = element_blank(),
    strip.text       = element_text(colour = "black", size = base_size, face = "bold"),
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

# ---- Palettes ---------------------------------------------------------------
# Nature journal palette — timecourse blues, controls warm
nature_pal <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                "#F39B7F", "#8491B4", "#91D1C2", "#DC0000")

sample_colors <- c(
  sc2  = "#E64B35",   # control — red
  sc5  = "#4DBBD5",   # timecourse — teal
  sc6  = "#00A087",   # timecourse — green
  sc78 = "#F39B7F",   # control — salmon
  sc9  = "#3C5488"    # timecourse — navy
)

# Timepoint order for timecourse plots
tp_order <- c("Control", "Sen 4hrs", "Sen D1", "Sen D2", "Sen D3", "Sen D4", "Sen D7")
tp_colors <- setNames(
  viridis::viridis(7, option = "D"),
  tp_order
)

# ---- Combine metadata -------------------------------------------------------
meta <- do.call(rbind, lapply(seu_list, function(obj) {
  obj@meta.data[, c("sample", "role", "condition",
                    "nCount_RNA", "nFeature_RNA", "pct_mt", "pct_ribo")]
}))
meta$sample <- factor(meta$sample, levels = c("sc5", "sc6", "sc9", "sc2", "sc78"))

# ---- Fig 1: Per-sample QC violins (3-panel) ---------------------------------
p_counts <- ggplot(meta, aes(x = sample, y = log1p(nCount_RNA), fill = sample)) +
  geom_violin(scale = "width", alpha = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = sample_colors) +
  labs(x = NULL, y = "log\u2081\u208a\u2081(UMI counts)") +
  theme_publication() +
  theme(legend.position = "none")

p_genes <- ggplot(meta, aes(x = sample, y = log1p(nFeature_RNA), fill = sample)) +
  geom_violin(scale = "width", alpha = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = sample_colors) +
  labs(x = NULL, y = "log\u2081\u208a\u2081(genes detected)") +
  theme_publication() +
  theme(legend.position = "none")

p_mt <- ggplot(meta, aes(x = sample, y = pct_mt, fill = sample)) +
  geom_violin(scale = "width", alpha = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = sample_colors) +
  labs(x = NULL, y = "Mitochondrial reads (%)") +
  theme_publication() +
  theme(legend.position = "none")

fig1 <- p_counts | p_genes | p_mt
ggsave(file.path(out_dir, "fig_qc_per_sample.svg"), fig1,
       width = 174, height = 70, units = "mm", bg = "white")
ggsave(file.path(out_dir, "fig_qc_per_sample.png"), fig1,
       width = 174, height = 70, units = "mm", dpi = 300, bg = "white")

# ---- Fig 2: UMI vs genes scatter, colored by MT%, per sample ---------------
p_scatter <- ggplot(meta, aes(x = log1p(nCount_RNA), y = log1p(nFeature_RNA),
                               color = pct_mt)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_viridis_c(option = "magma", name = "MT%") +
  facet_wrap(~sample, nrow = 2) +
  labs(x = "log\u2081\u208a\u2081(UMI counts)", y = "log\u2081\u208a\u2081(genes detected)") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "fig_qc_scatter.svg"), p_scatter,
       width = 174, height = 110, units = "mm", bg = "white")
ggsave(file.path(out_dir, "fig_qc_scatter.png"), p_scatter,
       width = 174, height = 110, units = "mm", dpi = 300, bg = "white")

# ---- Fig 3: Timecourse-only — UMI ridge plot by condition ------------------
tc_meta <- meta[meta$role == "timecourse" & meta$condition %in% tp_order, ]
tc_meta$condition <- factor(tc_meta$condition, levels = rev(tp_order))

p_ridge <- ggplot(tc_meta, aes(x = log1p(nCount_RNA), y = condition,
                                fill = condition)) +
  geom_density_ridges(alpha = 0.8, scale = 0.9, rel_min_height = 0.01,
                      bandwidth = 0.08) +
  scale_fill_manual(values = rev(tp_colors)) +
  facet_wrap(~sample, nrow = 1) +
  labs(x = "log\u2081\u208a\u2081(UMI counts)", y = NULL) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "fig_qc_ridge_timecourse.svg"), p_ridge,
       width = 174, height = 85, units = "mm", bg = "white")
ggsave(file.path(out_dir, "fig_qc_ridge_timecourse.png"), p_ridge,
       width = 174, height = 85, units = "mm", dpi = 300, bg = "white")

# ---- Fig 4: Cell counts per condition (timecourse) -------------------------
tc_counts <- tc_meta %>%
  group_by(sample, condition) %>%
  summarise(n_cells = n(), .groups = "drop")
tc_counts$condition <- factor(tc_counts$condition, levels = tp_order)

p_ncells <- ggplot(tc_counts, aes(x = condition, y = n_cells, fill = condition)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = tp_colors) +
  facet_wrap(~sample, nrow = 1) +
  labs(x = NULL, y = "Cells") +
  theme_publication() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

ggsave(file.path(out_dir, "fig_qc_cell_counts.svg"), p_ncells,
       width = 174, height = 70, units = "mm", bg = "white")
ggsave(file.path(out_dir, "fig_qc_cell_counts.png"), p_ncells,
       width = 174, height = 70, units = "mm", dpi = 300, bg = "white")

cat("QC figures saved to:", out_dir, "\n")
cat("Files: fig_qc_per_sample, fig_qc_scatter, fig_qc_ridge_timecourse, fig_qc_cell_counts\n")
