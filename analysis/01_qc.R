# 01_qc.R — Load GSE223128 into Seurat using rhdf5 + raw MTX
# Strategy: extract author-QC'd cell barcodes from h5ad, filter raw MTX, build Seurat objects

library(Seurat)
library(Matrix)
library(rhdf5)
library(ggplot2)
library(patchwork)

raw_dir <- "D:/Projects/Topic Modeling/data/raw/extracted"
out_dir <- "D:/Projects/Topic Modeling/analysis/qc_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_meta <- list(
  sc2  = list(prefix = "GSM6940120_sc2",
              h5ad   = "GSM6940120_sc2_adata_named_filtered.h5ad",
              bc_csv = "GSM6940120_sc2_bc_samples.csv.gz",
              role   = "control"),
  sc5  = list(prefix = "GSM6940121_sc5",
              h5ad   = "GSM6940121_sc5_anndata_samples_named_filtered.h5ad",
              bc_csv = "GSM6940121_sc5_bc_samples.csv.gz",
              role   = "timecourse"),
  sc6  = list(prefix = "GSM6940122_sc6",
              h5ad   = "GSM6940122_sc6_anndata_named_filtered.h5ad",
              bc_csv = "GSM6940122_sc6_bc_samples.csv.gz",
              role   = "timecourse"),
  sc78 = list(prefix = "GSM6940123_sc7_8",
              h5ad   = "GSM6940123_sc7_8_adata_named_filtered.h5ad",
              bc_csv = "GSM6940123_sc7_8_bc_samples.csv.gz",
              role   = "control"),
  sc9  = list(prefix = "GSM6940124_sc9",
              h5ad   = "GSM6940124_sc9_anndata_samples_named_filtered_FigS3f-j_.h5ad",
              bc_csv = "GSM6940124_sc9_bc_samples.csv.gz",
              role   = "timecourse")
)

# ---- Load one sample --------------------------------------------------------
load_sample <- function(name, info) {
  cat("\nLoading", name, "...\n")
  p <- file.path(raw_dir, info$prefix)

  # 1. Get author-filtered cell barcodes from h5ad
  h5_path      <- file.path(raw_dir, info$h5ad)
  author_cells <- h5read(h5_path, "/obs/_index")
  cat("  Author-filtered cells:", length(author_cells), "\n")

  # 2. Load raw MTX
  mat      <- readMM(gzcon(gzfile(paste0(p, "_matrix.mtx.gz"))))
  barcodes <- read.table(gzfile(paste0(p, "_barcodes.tsv.gz")), header = FALSE)$V1
  features <- read.table(gzfile(paste0(p, "_features.tsv.gz")), header = FALSE,
                         sep = "\t", col.names = c("ensembl", "symbol", "type"))

  rownames(mat) <- make.unique(features$symbol)
  colnames(mat) <- barcodes

  # 3. Filter to author-QC'd cells
  keep <- barcodes %in% author_cells
  mat  <- mat[, keep]
  cat("  Cells retained from MTX:", ncol(mat), "\n")

  # 4. Attach condition labels via bc_samples mapping
  bc_map <- read.csv(file.path(raw_dir, info$bc_csv))
  names(bc_map) <- c("well_bc", "condition")
  bc_map$condition <- trimws(bc_map$condition)
  bc_map$condition <- sub("^control$", "Control", bc_map$condition)

  well_bc   <- sub(".*_", "", colnames(mat))
  condition <- bc_map$condition[match(well_bc, bc_map$well_bc)]
  cat("  Conditions:", paste(sort(unique(na.omit(condition))), collapse = ", "), "\n")

  # 5. Create Seurat object
  obj <- CreateSeuratObject(counts = mat, project = name, min.cells = 0, min.features = 0)
  obj$condition <- condition
  obj$sample    <- name
  obj$role      <- info$role
  obj$pct_mt    <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$pct_ribo  <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

  obj
}

# ---- Load all samples -------------------------------------------------------
cat("=== Loading all samples ===\n")
seu_list <- mapply(load_sample, names(sample_meta), sample_meta, SIMPLIFY = FALSE)

# ---- QC summary -------------------------------------------------------------
cat("\n=== QC Summary ===\n")
for (s in names(seu_list)) {
  obj <- seu_list[[s]]
  cat(sprintf("\n%s (%s): %d cells, %d genes\n", s, unique(obj$role), ncol(obj), nrow(obj)))
  cat(sprintf("  nCount_RNA  : median=%d  [%d, %d]\n",
              as.integer(median(obj$nCount_RNA)), min(obj$nCount_RNA), max(obj$nCount_RNA)))
  cat(sprintf("  nFeature_RNA: median=%d  [%d, %d]\n",
              as.integer(median(obj$nFeature_RNA)), min(obj$nFeature_RNA), max(obj$nFeature_RNA)))
  cat(sprintf("  pct_mt      : median=%.2f%%  max=%.2f%%\n",
              median(obj$pct_mt), max(obj$pct_mt)))
  cat(sprintf("  pct_ribo    : median=%.2f%%  max=%.2f%%\n",
              median(obj$pct_ribo), max(obj$pct_ribo)))
}

# ---- Visualize QC -----------------------------------------------------------
cat("\n=== Generating QC plots ===\n")

meta <- do.call(rbind, lapply(seu_list, function(obj) {
  obj@meta.data[, c("sample", "role", "condition",
                    "nCount_RNA", "nFeature_RNA", "pct_mt", "pct_ribo")]
}))

role_colors <- c("timecourse" = "#4393c3", "control" = "#d6604d")

p_counts <- ggplot(meta, aes(x = sample, y = log1p(nCount_RNA), fill = role)) +
  geom_violin(scale = "width", alpha = 0.85) +
  scale_fill_manual(values = role_colors) +
  labs(title = "UMI counts (author-QC'd cells)", x = NULL, y = "log1p(UMI)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_genes <- ggplot(meta, aes(x = sample, y = log1p(nFeature_RNA), fill = role)) +
  geom_violin(scale = "width", alpha = 0.85) +
  scale_fill_manual(values = role_colors) +
  labs(title = "Genes detected", x = NULL, y = "log1p(genes)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_mt <- ggplot(meta, aes(x = sample, y = pct_mt, fill = role)) +
  geom_violin(scale = "width", alpha = 0.85) +
  scale_fill_manual(values = role_colors) +
  labs(title = "Mitochondrial %", x = NULL, y = "% MT reads") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_scatter <- ggplot(meta, aes(x = log1p(nCount_RNA), y = log1p(nFeature_RNA),
                               color = pct_mt)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_viridis_c(option = "magma") +
  facet_wrap(~sample, nrow = 2) +
  labs(title = "Counts vs genes (colored by MT%)", color = "MT%") +
  theme_bw()

# Per-sample, per-condition cell counts
count_tbl <- as.data.frame(table(meta$sample, meta$condition))
names(count_tbl) <- c("sample", "condition", "n_cells")
p_cellcounts <- ggplot(count_tbl[count_tbl$n_cells > 0, ],
                        aes(x = condition, y = n_cells, fill = sample)) +
  geom_col(position = "dodge") +
  facet_wrap(~sample, scales = "free_x", nrow = 2) +
  labs(title = "Cells per condition", x = NULL, y = "Cell count") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "none")

ggsave(file.path(out_dir, "qc_violin_counts.png"),  p_counts,     width = 8,  height = 5, dpi = 150)
ggsave(file.path(out_dir, "qc_violin_genes.png"),   p_genes,      width = 8,  height = 5, dpi = 150)
ggsave(file.path(out_dir, "qc_violin_mt.png"),      p_mt,         width = 8,  height = 5, dpi = 150)
ggsave(file.path(out_dir, "qc_scatter.png"),        p_scatter,    width = 12, height = 6, dpi = 150)
ggsave(file.path(out_dir, "qc_cells_per_cond.png"), p_cellcounts, width = 12, height = 6, dpi = 150)

# ---- Save -------------------------------------------------------------------
saveRDS(seu_list, file.path(out_dir, "seu_list_qc.rds"))
cat("\nSaved to:", file.path(out_dir, "seu_list_qc.rds"), "\n")
cat("QC plots saved to:", out_dir, "\n")
cat("\nDone.\n")
