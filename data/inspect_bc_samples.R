for (f in c(
  "GSM6940121_sc5_bc_samples.csv.gz",
  "GSM6940122_sc6_bc_samples.csv.gz",
  "GSM6940124_sc9_bc_samples.csv.gz"
)) {
  path <- paste0("D:/Projects/Topic Modeling/data/raw/extracted/", f)
  df <- read.csv(path)
  cat("===", f, "===\n")
  print(table(df[, 2]))
  cat("\n")
}
