library(rhdf5)

h5 <- "D:/Projects/Topic Modeling/data/raw/extracted/GSM6940120_sc2_adata_named_filtered.h5ad"
cat("Top-level keys:\n")
print(h5ls(h5, recursive = FALSE)$name)

cat("\nobs keys (cell metadata):\n")
obs <- h5ls(h5, recursive = FALSE)
print(h5ls(h5)[h5ls(h5)$group == "/obs", c("name", "otype", "dclass")])
