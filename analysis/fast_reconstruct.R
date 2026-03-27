# fast_reconstruct.R — rebuild kselection.rds for sc78 and sc5 from known values.
# Uses approximate K=8 perplexity for sc78 (4.4810 from original run summary).
# Uses backed-up sc5 original + gap fill values from log output.

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"

# ----- sc78 ------------------------------------------------------------------
# Original run: K=6→4.494, K=8→4.481 (min), K=10→4.560
# Gap fill:     K=7→4.48652, K=9→4.48094
sc78_k    <- c(6L, 7L, 8L, 9L, 10L)
sc78_perp <- c(4.49400, 4.48652, 4.48100, 4.48094, 4.56000)
sc78_ll   <- rep(NA_real_, 5)   # in-sample loglik not critical for K selection

ksel78 <- list(k_range = sc78_k, perplexities = sc78_perp, logliks = sc78_ll)
saveRDS(ksel78, file.path(v2_dir, "sc78", "kselection.rds"))

cat("sc78 reconstructed:\n")
for (i in seq_along(sc78_k))
  cat(sprintf("  K = %2d : %.5f\n", sc78_k[i], sc78_perp[i]))
cat("Minimum at K =", sc78_k[which.min(sc78_perp)], "\n")

# ----- sc5 -------------------------------------------------------------------
# Original: K=3→2.43219, K=4→2.40381, K=5→2.38921, K=6→2.38045,
#           K=8→2.36717, K=10→2.36170, K=12→2.44000, K=15→2.74279
# Gap fill: K=9→2.36300, K=11→2.35980
sc5_orig <- readRDS("/tmp/sc5_kselection_orig.rds")

gap5_k    <- c(9L, 11L)
gap5_perp <- c(2.36300, 2.35980)

all_k5 <- sort(union(sc5_orig$k_range, gap5_k))
new_p5 <- numeric(length(all_k5))
new_l5 <- numeric(length(all_k5))

for (i in seq_along(all_k5)) {
  k <- all_k5[i]
  if (k %in% sc5_orig$k_range) {
    idx <- which(sc5_orig$k_range == k)
    new_p5[i] <- sc5_orig$perplexities[idx]
    new_l5[i] <- sc5_orig$logliks[idx]
  } else {
    idx <- which(gap5_k == k)
    new_p5[i] <- gap5_perp[idx]
  }
}

ksel5 <- list(k_range = all_k5, perplexities = new_p5, logliks = new_l5)
saveRDS(ksel5, file.path(v2_dir, "sc5", "kselection.rds"))

cat("\nsc5 reconstructed:\n")
for (i in seq_along(all_k5))
  cat(sprintf("  K = %2d : %.5f\n", all_k5[i], new_p5[i]))
cat("Minimum at K =", all_k5[which.min(new_p5)], "\n")

cat("\nDone.\n")
