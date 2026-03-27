# 03_lda_poc_v2.R — LDA topic model with variance-based gene selection
#
# Gene filtering: within-condition overdispersion UNION between-condition variance
#   - No hard detection cutoffs
#   - MT- and RP[SL] genes included (compete on equal footing)
#   - Library-size normalised before computing variance metrics
#
# Everything downstream identical to v1: GPU MU warm-start -> fastTopics SCD

library(Seurat)
library(fastTopics)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)

qc_dir  <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"
out_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"
dir.create(file.path(out_dir, "models"), recursive = TRUE, showWarnings = FALSE)

repo_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling"
py_script <- file.path(repo_dir, "analysis/gpu_nmf.py")

# ---- GPU detection (unchanged) -----------------------------------------------
detect_r_gpu <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) return(FALSE)
  torch_dir  <- system.file(package = "torch")
  lantern_ok <- any(file.exists(
    file.path(torch_dir, "deps", "lantern.dll"),
    file.path(torch_dir, "deps", "liblantern.so"),
    file.path(torch_dir, "deps", "liblantern.dylib")
  ))
  if (!lantern_ok) return(FALSE)
  rscript <- Sys.which("Rscript")
  result  <- tryCatch(
    system2(rscript,
            c("--no-save", "--no-environ", "-e",
              "library(torch,quietly=TRUE);q(status=if(cuda_is_available())0L else 1L)"),
            stdout = FALSE, stderr = FALSE, timeout = 30),
    error = function(e) 1L, warning = function(e) 1L
  )
  isTRUE(result == 0L)
}

detect_py_gpu <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  py312 <- "C:/Users/adb/AppData/Local/Programs/Python/Python312/python.exe"
  if (file.exists(py312)) reticulate::use_python(py312, required = FALSE)
  py <- tryCatch(reticulate::py_config(), error = function(e) NULL)
  if (is.null(py)) return(FALSE)
  tryCatch({
    reticulate::source_python(py_script, envir = globalenv())
    isTRUE(cuda_available())
  }, error = function(e) FALSE)
}

use_r_gpu  <- detect_r_gpu()
use_py_gpu <- if (!use_r_gpu) detect_py_gpu() else FALSE
use_gpu    <- use_r_gpu || use_py_gpu

if (use_r_gpu) {
  suppressPackageStartupMessages(library(torch))
  cat("GPU backend: R torch |", cuda_get_device_name(), "\n")
} else if (use_py_gpu) {
  cat("GPU backend: Python torch (reticulate) |", cuda_device_name(), "\n")
} else {
  cat("GPU backend: none — running on CPU\n")
}

theme_publication <- function(base_size = 11) {
  theme_cowplot(base_size) +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.5),
    axis.ticks       = element_line(colour = "black", linewidth = 0.5),
    axis.text        = element_text(colour = "black", size = base_size),
    axis.title       = element_text(colour = "black", size = base_size + 1),
    strip.background = element_blank(),
    strip.text       = element_text(colour = "black", size = base_size, face = "bold"),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.background   = element_blank(),
    plot.title        = element_blank(),
    plot.margin       = margin(8, 8, 8, 8)
  )
}

# ---- Variance-based gene selection -------------------------------------------
#
# Two axes:
#   within_disp  : mean overdispersion (var/mean - 1) within each condition
#                  captures cell-to-cell heterogeneity within conditions
#   between_var  : variance of condition-level means
#                  captures systematic condition differences (DE genes, SASP, etc.)
#
# Gene selected if in top n_within by within_disp OR top n_between by between_var.
# No hard detection cutoffs. MT- and RP[SL] genes compete on equal footing.

filter_genes_variance <- function(counts, meta, condition_col = "condition",
                                   n_within = 7000, n_between = 5000) {
  cat("  Total genes:", nrow(counts), "\n")

  # Library-size normalise (preserves sparsity — zeros stay zero)
  lib_size <- Matrix::colSums(counts)
  med_lib  <- median(lib_size)
  norm     <- Matrix::t(Matrix::t(counts) / lib_size * med_lib)

  conditions <- sort(unique(na.omit(meta[[condition_col]])))
  cat("  Conditions:", paste(conditions, collapse = ", "), "\n")

  # ---- Within-condition overdispersion ----------------------------------------
  # For each condition: var/mean - 1 on library-normalised counts
  # sum then divide (faster than Reduce for many conditions)
  disp_sum <- rep(0, nrow(norm))
  for (cond in conditions) {
    cells <- which(meta[[condition_col]] == cond)
    if (length(cells) < 2) next
    sub <- norm[, cells, drop = FALSE]
    m   <- Matrix::rowMeans(sub)
    # Sparse-safe variance: E[X^2] - E[X]^2
    ss  <- Matrix::rowSums(sub * sub)
    v   <- ss / length(cells) - m^2
    disp_sum <- disp_sum + v / pmax(m, 1e-10) - 1
  }
  within_disp <- disp_sum / length(conditions)

  # ---- Between-condition variance ---------------------------------------------
  cond_means  <- vapply(conditions, function(cond) {
    cells <- which(meta[[condition_col]] == cond)
    Matrix::rowMeans(norm[, cells, drop = FALSE])
  }, numeric(nrow(norm)))
  # cond_means is genes x conditions
  between_var <- rowSums((cond_means - rowMeans(cond_means))^2) / (ncol(cond_means) - 1)

  # ---- Union ------------------------------------------------------------------
  keep_within  <- rank(-within_disp,  ties.method = "min") <= n_within
  keep_between <- rank(-between_var,  ties.method = "min") <= n_between
  keep         <- keep_within | keep_between

  mt_kept <- sum(keep & grepl("^MT-",    rownames(counts)))
  rp_kept <- sum(keep & grepl("^RP[SL]", rownames(counts)))

  cat("  Top within-dispersion:  ", sum(keep_within),  "genes\n")
  cat("  Top between-variance:   ", sum(keep_between), "genes\n")
  cat("  Union (kept):           ", sum(keep),         "genes\n")
  cat("  MT- genes kept:         ", mt_kept, "\n")
  cat("  RP[SL] genes kept:      ", rp_kept, "\n")

  counts[keep, ]
}

# ---- GPU warm-start (unchanged) ----------------------------------------------
gpu_warm_start <- function(X, k, n_iter_gpu, seed) {
  if (use_r_gpu) {
    n <- nrow(X); p <- ncol(X)
    set.seed(seed)
    W0 <- matrix(rgamma(n * k, 1, 1), n, k)
    H0 <- matrix(rgamma(k * p, 1, 1), k, p)
    device <- torch_device("cuda")
    V <- torch_tensor(as.matrix(X), dtype = torch_float32(), device = device)
    W <- torch_tensor(W0, dtype = torch_float32(), device = device)
    H <- torch_tensor(H0, dtype = torch_float32(), device = device)
    rm(W0, H0)
    eps <- 1e-10
    for (iter in seq_len(n_iter_gpu)) {
      WH <- torch_mm(W, H) + eps
      H  <- H * torch_mm(W$t(), V / WH) / (W$sum(dim = 1L)$unsqueeze(2L) + eps)
      WH <- torch_mm(W, H) + eps
      W  <- W * torch_mm(V / WH, H$t()) / (H$sum(dim = 2L)$unsqueeze(1L) + eps)
      if (iter %% 50 == 0) cat("    R-GPU iter", iter, "\n")
    }
    L_warm <- as.matrix(W$cpu()); F_warm <- t(as.matrix(H$cpu()))
    rm(V, W, H); torch_cuda_empty_cache()
  } else {
    result <- poisson_nmf_gpu(as.matrix(X), k = as.integer(k),
                              n_iter = as.integer(n_iter_gpu),
                              seed   = as.integer(seed))
    L_warm <- result$L; F_warm <- result$F
  }
  rownames(L_warm) <- rownames(X); rownames(F_warm) <- colnames(X)
  colnames(L_warm) <- paste0("k", seq_len(k))
  colnames(F_warm) <- paste0("k", seq_len(k))
  list(L = L_warm, F = F_warm)
}

fit_topic_model_gpu <- function(X, k, n_iter_gpu = 150, n_iter_final = 100,
                                seed = 42) {
  if (!use_gpu) {
    cat("  [CPU] fit_topic_model k =", k, "\n")
    set.seed(seed)
    return(fit_topic_model(X, k = k,
                           numiter.main = n_iter_final + n_iter_gpu,
                           verbose = "none"))
  }
  cat("  [GPU] MU warm-start k =", k, "for", n_iter_gpu, "iters...\n")
  warm <- gpu_warm_start(X, k, n_iter_gpu, seed)
  cat("  [SCD] Refining with fastTopics k =", k, "for", n_iter_final, "iters...\n")
  set.seed(seed)
  fit0 <- init_poisson_nmf(X, F = warm$F, L = warm$L)
  fit_poisson_nmf(X, fit0 = fit0, numiter = n_iter_final,
                  method = "scd", verbose = "none")
}

# ---- Cross-validated perplexity ----------------------------------------------
cv_perplexity <- function(counts_t, k, n_iter_gpu = 100, n_iter_final = 50,
                          holdout_frac = 0.2, seed = 42) {
  set.seed(seed)
  X    <- as.matrix(counts_t)
  nz   <- which(X > 0, arr.ind = TRUE)
  mask <- nz[sample(nrow(nz), size = floor(nrow(nz) * holdout_frac)), , drop = FALSE]

  X_train       <- X; X_train[mask] <- 0L
  X_test_counts <- X[mask]
  test_total    <- sum(X_test_counts)

  fit    <- fit_topic_model_gpu(Matrix::Matrix(X_train, sparse = TRUE), k = k,
                                n_iter_gpu = n_iter_gpu, n_iter_final = n_iter_final,
                                seed = seed)
  lambda <- pmax((fit$L %*% t(fit$F))[mask], 1e-10)
  ll     <- sum(dpois(X_test_counts, lambda, log = TRUE))
  list(perplexity = exp(-ll / test_total), fit = fit)
}

# ---- K selection -------------------------------------------------------------
select_k <- function(counts_t, k_range = c(3, 4, 5, 6, 8, 10, 12, 15),
                     out_prefix) {
  cat("  Fitting models for K =", paste(k_range, collapse = ", "), "...\n")

  fits    <- list()
  logliks <- numeric(length(k_range))
  perps   <- numeric(length(k_range))

  for (i in seq_along(k_range)) {
    k    <- k_range[i]
    ckpt <- paste0(out_prefix, "_fit_k", k, ".rds")

    if (file.exists(ckpt)) {
      cat("    K =", k, "... loading from checkpoint\n")
      fits[[as.character(k)]] <- readRDS(ckpt)
    } else {
      cat("    K =", k, "...\n")
      fit <- fit_topic_model_gpu(counts_t, k = k,
                                 n_iter_gpu = 100, n_iter_final = 50)
      saveRDS(fit, ckpt)
      fits[[as.character(k)]] <- fit
    }

    cat("      CV perplexity k =", k, "...\n")
    cv       <- cv_perplexity(counts_t, k = k)
    perps[i] <- cv$perplexity
    logliks[i] <- mean(loglik_multinom_topic_model(counts_t, fits[[as.character(k)]]))
  }

  cat("  Log-likelihoods and perplexity:\n")
  for (i in seq_along(k_range)) {
    cat("    K =", k_range[i], ": loglik =", round(logliks[i], 2),
        "  perplexity =", round(perps[i], 3), "\n")
  }

  ll_df  <- data.frame(K = k_range, loglik = logliks)
  p_ll   <- ggplot(ll_df, aes(x = K, y = loglik)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of topics (K)", y = "Log-likelihood") +
    theme_publication()
  ggsave(paste0(out_prefix, "_loglik.png"), p_ll,
         width = 85, height = 70, units = "mm", dpi = 300, bg = "white")
  ggsave(paste0(out_prefix, "_loglik.svg"), p_ll,
         width = 85, height = 70, units = "mm", bg = "white")

  perp_df <- data.frame(K = k_range, perplexity = perps)
  p_perp  <- ggplot(perp_df, aes(x = K, y = perplexity)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of topics (K)", y = "Perplexity") +
    theme_publication()
  ggsave(paste0(out_prefix, "_perplexity.png"), p_perp,
         width = 85, height = 70, units = "mm", dpi = 300, bg = "white")
  ggsave(paste0(out_prefix, "_perplexity.svg"), p_perp,
         width = 85, height = 70, units = "mm", bg = "white")

  list(fits = fits, logliks = logliks, perplexities = perps)
}

# ---- Topic visualisation helpers (unchanged) ---------------------------------
plot_topic_structure <- function(fit, meta, condition_col, out_prefix) {
  theta     <- poisson2multinom(fit)$L
  df        <- as.data.frame(theta)
  df$cell      <- rownames(df)
  df$condition <- meta[[condition_col]][match(df$cell, rownames(meta))]
  topic_cols   <- grep("^k", names(df), value = TRUE)
  df_long      <- tidyr::pivot_longer(df, cols = all_of(topic_cols),
                                      names_to = "topic", values_to = "weight")
  dom_topic    <- apply(theta, 1, which.max)
  cell_order   <- rownames(theta)[order(df$condition, dom_topic)]
  df_long$cell <- factor(df_long$cell, levels = cell_order)

  p <- ggplot(df_long, aes(x = cell, y = weight, fill = topic)) +
    geom_col(width = 1) +
    scale_fill_viridis_d(option = "H") +
    facet_grid(~condition, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Topic weight", fill = "Topic") +
    theme_publication() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.spacing = unit(0.5, "mm"))
  ggsave(paste0(out_prefix, "_structure.png"), p,
         width = 174, height = 70, units = "mm", dpi = 300, bg = "white")
  ggsave(paste0(out_prefix, "_structure.svg"), p,
         width = 174, height = 70, units = "mm", bg = "white")
  invisible(theta)
}

plot_top_genes <- function(fit, k, out_prefix) {
  top <- lapply(seq_len(k), function(i) {
    freq <- fit$F[, i]
    head(sort(freq, decreasing = TRUE), 10)
  })
  cat("  Top genes per topic:\n")
  for (i in seq_len(k)) {
    cat("    Topic", i, ":", paste(names(top[[i]]), collapse = ", "), "\n")
  }
  top_genes <- unique(unlist(lapply(top, names)))
  F_mat     <- fit$F[top_genes, , drop = FALSE]
  df_heat   <- as.data.frame(as.matrix(F_mat)) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "topic", values_to = "freq")
  p <- ggplot(df_heat, aes(x = topic, y = gene, fill = freq)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma") +
    labs(x = "Topic", y = NULL, fill = "Gene freq.") +
    theme_publication() +
    theme(axis.text.y = element_text(size = 7))
  ggsave(paste0(out_prefix, "_top_genes_heatmap.png"), p,
         width = 100, height = 120, units = "mm", dpi = 300, bg = "white")
  ggsave(paste0(out_prefix, "_top_genes_heatmap.svg"), p,
         width = 100, height = 120, units = "mm", bg = "white")
}

# =============================================================================
# Load data
# =============================================================================
cat("Loading QC objects...\n")
seu_list <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))

# =============================================================================
# SC2 — Nutlin, OxStress, Quiescence (3,313 cells)
# =============================================================================
cat("\n=== SC2 ===\n")
sc2     <- seu_list[["sc2"]]
counts2 <- GetAssayData(sc2, layer = "counts")

cat("Variance-based gene selection for sc2...\n")
counts2f <- filter_genes_variance(counts2, sc2@meta.data, condition_col = "condition")
counts2t <- t(counts2f)

cat("Running K selection for sc2...\n")
sc2_ks <- select_k(counts2t,
                   k_range    = c(3, 4, 5, 6, 8, 10, 12, 15),
                   out_prefix = file.path(out_dir, "models", "sc2"))

# Save perplexity data
saveRDS(list(k_range      = c(3, 4, 5, 6, 8, 10, 12, 15),
             logliks      = sc2_ks$logliks,
             perplexities = sc2_ks$perplexities),
        file.path(out_dir, "models", "sc2_kselection.rds"))

# Structure plot for K with min perplexity
best_k_sc2 <- c(3, 4, 5, 6, 8, 10, 12, 15)[which.min(sc2_ks$perplexities)]
cat("Preliminary best K sc2 (min perplexity):", best_k_sc2, "\n")
fit_sc2_best <- sc2_ks$fits[[as.character(best_k_sc2)]]

cat("Plotting sc2 topic structure (K =", best_k_sc2, ")...\n")
plot_topic_structure(fit_sc2_best, sc2@meta.data, "condition",
                     file.path(out_dir, paste0("sc2_k", best_k_sc2)))
plot_top_genes(fit_sc2_best, best_k_sc2,
               file.path(out_dir, paste0("sc2_k", best_k_sc2)))

# =============================================================================
# SC78 — Control, Nutlin, OxStress, Quiescence x 2 reps (18,429 cells)
# =============================================================================
cat("\n=== SC78 ===\n")
sc78     <- seu_list[["sc78"]]
counts78 <- GetAssayData(sc78, layer = "counts")

cat("Variance-based gene selection for sc78...\n")
counts78f <- filter_genes_variance(counts78, sc78@meta.data, condition_col = "condition")
counts78t <- t(counts78f)

cat("Running K selection for sc78...\n")
sc78_ks <- select_k(counts78t,
                    k_range    = c(3, 4, 5, 6, 8, 10, 12, 15),
                    out_prefix = file.path(out_dir, "models", "sc78"))

saveRDS(list(k_range      = c(3, 4, 5, 6, 8, 10, 12, 15),
             logliks      = sc78_ks$logliks,
             perplexities = sc78_ks$perplexities),
        file.path(out_dir, "models", "sc78_kselection.rds"))

best_k_sc78 <- c(3, 4, 5, 6, 8, 10, 12, 15)[which.min(sc78_ks$perplexities)]
cat("Preliminary best K sc78 (min perplexity):", best_k_sc78, "\n")
fit_sc78_best <- sc78_ks$fits[[as.character(best_k_sc78)]]

cat("Plotting sc78 topic structure (K =", best_k_sc78, ")...\n")
plot_topic_structure(fit_sc78_best, sc78@meta.data, "condition",
                     file.path(out_dir, paste0("sc78_k", best_k_sc78)))
plot_top_genes(fit_sc78_best, best_k_sc78,
               file.path(out_dir, paste0("sc78_k", best_k_sc78)))

cat("\n=== Done. Outputs in:", out_dir, "===\n")
