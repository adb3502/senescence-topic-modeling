f <- read.table(gzfile("D:/Projects/Topic Modeling/data/raw/extracted/GSM6940120_sc2_features.tsv.gz"),
               header = FALSE, sep = "\t")
cat("Columns:", ncol(f), "\n")
print(head(f, 3))
