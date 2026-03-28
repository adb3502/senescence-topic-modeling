library(fastTopics)
library(msigdbr)

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"
fit5 <- readRDS(file.path(v2_dir, "sc5", "fit_k10.rds"))

# Topic L = Topic 8
top50_L <- names(sort(fit5$F[, 8], decreasing = TRUE))[1:50]
cat("Top 50 genes for L (Topic 8):\n")
cat(paste(top50_L, collapse = ", "), "\n\n")

# Load Hallmark Hypoxia
h_df     <- msigdbr(species = "Homo sapiens", collection = "H")
hypoxia  <- h_df$gene_symbol[h_df$gs_name == "HALLMARK_HYPOXIA"]
glycol   <- h_df$gene_symbol[h_df$gs_name == "HALLMARK_GLYCOLYSIS"]
mtorc1   <- h_df$gene_symbol[h_df$gs_name == "HALLMARK_MTORC1_SIGNALING"]

cat("Overlap with Hallmark HYPOXIA (", length(hypoxia), "genes ):\n")
ov_h <- intersect(top50_L, hypoxia)
cat(length(ov_h), "overlapping:", paste(ov_h, collapse = ", "), "\n\n")

cat("Overlap with Hallmark GLYCOLYSIS (", length(glycol), "genes ):\n")
ov_g <- intersect(top50_L, glycol)
cat(length(ov_g), "overlapping:", paste(ov_g, collapse = ", "), "\n\n")

cat("Overlap with Hallmark MTORC1_SIGNALING (", length(mtorc1), "genes ):\n")
ov_m <- intersect(top50_L, mtorc1)
cat(length(ov_m), "overlapping:", paste(ov_m, collapse = ", "), "\n\n")

# How many hypoxia genes are in top 100, 200, 500?
for (n in c(50, 100, 200, 500)) {
  topN <- names(sort(fit5$F[, 8], decreasing = TRUE))[1:n]
  cat(sprintf("Top %d genes: %d / %d Hypoxia genes present\n",
              n, length(intersect(topN, hypoxia)), length(hypoxia)))
}
