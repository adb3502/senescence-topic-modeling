if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery", update = FALSE, ask = FALSE)
}

library(GEOquery)

# Fetch series metadata (soft file only, no expression matrices yet)
gse <- getGEO("GSE223128", GSEMatrix = FALSE, destdir = "D:/Users/adb/Topic Modeling/senescence-topic-modeling/data/raw")

# Print series-level info
cat("=== Series Title ===\n")
cat(Meta(gse)$title, "\n\n")

cat("=== Summary ===\n")
cat(Meta(gse)$summary, "\n\n")

cat("=== Overall Design ===\n")
cat(Meta(gse)$overall_design, "\n\n")

cat("=== Number of Samples ===\n")
cat(length(GSMList(gse)), "\n\n")

# Print supplementary files at series level
cat("=== Series-level Supplementary Files ===\n")
suppl <- Meta(gse)$supplementary_file
for (f in suppl) cat(" ", f, "\n")

cat("\n=== Sample-level Supplementary Files (first 5 samples) ===\n")
samples <- GSMList(gse)
for (i in seq_len(min(5, length(samples)))) {
  gsm <- samples[[i]]
  cat("\nSample:", names(samples)[i], "-", Meta(gsm)$title, "\n")
  cat("  Characteristics:", paste(Meta(gsm)$characteristics_ch1, collapse = " | "), "\n")
  sf <- Meta(gsm)$supplementary_file
  for (f in sf) cat("  File:", f, "\n")
}
