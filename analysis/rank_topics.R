library(fastTopics)
suppressMessages(library(Seurat))

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"

fit2  <- readRDS(file.path(v2_dir, "sc2",  "fit_k4.rds"))
fit78 <- readRDS(file.path(v2_dir, "sc78", "fit_k9.rds"))
fit5  <- readRDS(file.path(v2_dir, "sc5",  "fit_k10.rds"))

sc2_map  <- c("Topic 1"="F","Topic 2"="H","Topic 3"="A","Topic 4"="B")
sc78_map <- c("Topic 1"="G","Topic 2"="D","Topic 3"="E","Topic 4"="F",
              "Topic 5"="C","Topic 6"="H","Topic 7"="B","Topic 8"="I","Topic 9"="A")
sc5_map  <- c("Topic 1"="A","Topic 2"="B","Topic 3"="J","Topic 4"="F",
              "Topic 5"="G","Topic 6"="D","Topic 7"="K","Topic 8"="L",
              "Topic 9"="C","Topic 10"="E")

mean_prop <- function(fit, topic_map) {
  L <- fit$L / rowSums(fit$L)
  colnames(L) <- paste0("Topic ", seq_len(ncol(L)))
  v <- sapply(names(topic_map), function(t) mean(L[, t]))
  names(v) <- unname(topic_map)
  v
}

p2  <- mean_prop(fit2,  sc2_map)
p78 <- mean_prop(fit78, sc78_map)
p5  <- mean_prop(fit5,  sc5_map)

all_letters <- LETTERS[1:12]
totals <- sapply(all_letters, function(l) {
  sum(c(p2[names(p2)==l], p78[names(p78)==l], p5[names(p5)==l]))
})
n_samples <- sapply(all_letters, function(l) {
  sum(c(any(sc2_map==l), any(sc78_map==l), any(sc5_map==l)))
})
avg <- totals / n_samples
df <- data.frame(letter=all_letters, avg_proportion=round(avg, 4), n_samples=n_samples)
df <- df[order(-df$avg_proportion), ]
rownames(df) <- NULL
print(df)
