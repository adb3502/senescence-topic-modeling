library(Matrix)

base <- "D:/Projects/Topic Modeling/data/raw/extracted"

samples <- list(
  sc5 = "GSM6940121_sc5",
  sc6 = "GSM6940122_sc6",
  sc9 = "GSM6940124_sc9"
)

results <- lapply(names(samples), function(s) {
  p <- file.path(base, samples[[s]])
  cat("Loading", s, "...\n")

  mat      <- readMM(gzcon(gzfile(paste0(p, "_matrix.mtx.gz"))))
  barcodes <- read.table(gzfile(paste0(p, "_barcodes.tsv.gz")), header = FALSE)$V1
  features <- read.table(gzfile(paste0(p, "_features.tsv.gz")), header = FALSE)
  bc_map   <- read.csv(paste0(p, "_bc_samples.csv.gz"))

  rownames(mat) <- features$V2  # gene symbols in column 2
  colnames(mat) <- barcodes

  names(bc_map) <- c("well_bc", "condition")
  bc_map$condition <- sub("^control$", "Control", bc_map$condition)

  # SPLIT-seq: cell barcode = bc1_bc2_well_bc; extract last component
  well_bc_per_cell <- sub(".*_", "", barcodes)

  # Map to condition
  condition <- bc_map$condition[match(well_bc_per_cell, bc_map$well_bc)]

  cat("  Genes:", nrow(mat), "| Total barcodes:", ncol(mat), "\n")
  cat("  Cells with condition assigned:", sum(!is.na(condition)), "\n")
  cat("  Condition breakdown:\n")
  print(table(condition, useNA = "ifany"))

  # Keep only assigned cells
  keep    <- !is.na(condition)
  mat     <- mat[, keep]
  condition <- condition[keep]

  cat("  After filtering:", ncol(mat), "cells retained\n\n")

  list(mat = mat, features = features, condition = condition, sample = s)
})

names(results) <- names(samples)
cat("Done. Summary:\n")
for (s in names(results)) {
  cat(" ", s, ":", ncol(results[[s]]$mat), "cells\n")
}
