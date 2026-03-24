for (f in c(
  "GSM6940120_sc2_bc_samples.csv.gz",
  "GSM6940123_sc7_8_bc_samples.csv.gz"
)) {
  path <- paste0("D:/Projects/Topic Modeling/data/raw/extracted/", f)
  df <- read.csv(path)
  names(df) <- c("well_bc", "condition")
  cat("===", f, "===\n")
  print(table(df$condition))
  cat("\n")
}
