v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"
library(fastTopics)

fit5 <- readRDS(file.path(v2_dir, "sc5", "fit_k10.rds"))

sc5_map <- c("Topic 1"="A","Topic 2"="B","Topic 3"="J","Topic 4"="F",
             "Topic 5"="G","Topic 6"="D","Topic 7"="K","Topic 8"="L",
             "Topic 9"="C","Topic 10"="E")

for (tl in names(sc5_map)) {
  letter <- sc5_map[tl]
  i      <- as.integer(sub("Topic ","",tl))
  top    <- names(sort(fit5$F[,i], decreasing=TRUE))[1:15]
  cat(sprintf("\n%s (%s): %s\n", letter, tl, paste(top, collapse=", ")))
}
