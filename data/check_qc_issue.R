library(Matrix)
base <- "D:/Projects/Topic Modeling/data/raw/extracted"

# Load sc9 quickly and check count distribution
p <- file.path(base, "GSM6940124_sc9")
mat      <- readMM(gzcon(gzfile(paste0(p, "_matrix.mtx.gz"))))
barcodes <- read.table(gzfile(paste0(p, "_barcodes.tsv.gz")), header = FALSE)$V1
bc_map   <- read.csv(paste0(p, "_bc_samples.csv.gz"))
names(bc_map) <- c("well_bc", "condition")
bc_map$condition <- sub("^control$", "Control", bc_map$condition)
well_bc <- sub(".*_", "", barcodes)
condition <- bc_map$condition[match(well_bc, bc_map$well_bc)]
keep <- !is.na(condition)
mat <- mat[, keep]
n_counts <- Matrix::colSums(mat)

cat("sc9 count distribution:\n")
print(quantile(n_counts, probs = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)))
cat("\nCells with 0 counts:", sum(n_counts == 0), "\n")
cat("Cells with <10 counts:", sum(n_counts < 10), "\n")
cat("Cells with <50 counts:", sum(n_counts < 50), "\n")
cat("Cells with <100 counts:", sum(n_counts < 100), "\n")
