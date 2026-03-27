# 03_lda_poc.R — LDA topic model POC on sc2 and sc78 independently
# fastTopics (Poisson NMF = multinomial LDA at MLE)
# Raw counts, corpus-style gene filtering, K selection by held-out likelihood
# GPU acceleration via torch: MU warm-start → fastTopics SCD refinement

library(Seurat)
library(fastTopics)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)

qc_dir  <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"
out_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Detect GPU --------------------------------------------------------------
# Priority:
#   1. R torch (with Lantern backend installed)
#   2. Python torch via reticulate (requires python3 + torch CUDA in PATH)
#   3. CPU fallback (fastTopics SCD, already fast for our sizes)
#
# R torch 0.16.x segfaults on ANY R operation after library(torch) if the
# Lantern backend DLL is missing. Detect safely via subprocess.

repo_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling"
py_script <- file.path(repo_dir, "analysis/gpu_nmf.py")

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
  # Use Python 3.12 which has torch installed
  py312 <- "C:/Users/adb/AppData/Local/Programs/Python/Python312/python.exe"
  if (file.exists(py312)) {
    reticulate::use_python(py312, required = FALSE)
  }
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
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.text = element_text(colour = "black", size = base_size),
    axis.title = element_text(colour = "black", size = base_size + 1),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = base_size, face = "bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    plot.title = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )
}

# ---- Gene filtering (corpus analogy) ----------------------------------------
filter_genes <- function(counts, min_pct = 0.01, max_pct = 0.95) {
  n_cells   <- ncol(counts)
  detection <- Matrix::rowMeans(counts > 0)
  keep      <- detection >= min_pct & detection <= max_pct
  keep      <- keep & !grepl("^MT-|^RP[SL]", rownames(counts))
  cat("  Genes before filter:", nrow(counts), "\n")
  cat("  Genes after filter: ", sum(keep), "\n")
  counts[keep, ]
}

# ---- GPU-accelerated Poisson NMF (multiplicative updates) -------------------
# Strategy:
#   1. If R torch + CUDA: run MU updates natively
#   2. If Python torch + CUDA (reticulate): call gpu_nmf.py
#   3. CPU fallback: fastTopics SCD for all iterations
# In all cases MU warm-start → fastTopics SCD refinement when GPU is available.

gpu_warm_start <- function(X, k, n_iter_gpu, seed) {
  if (use_r_gpu) {
    # --- R torch path ---
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
    # --- Python torch path via reticulate ---
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

# ---- K selection via held-out likelihood ------------------------------------
# ---- Cross-validated perplexity ---------------------------------------------
cv_perplexity <- function(counts_t, k, n_iter_gpu = 100, n_iter_final = 50,
                          holdout_frac = 0.2, seed = 42) {
  set.seed(seed)
  X <- as.matrix(counts_t)

  # Hold out ~20% of entries (sparse-safe: only mask observed counts)
  nz       <- which(X > 0, arr.ind = TRUE)
  mask_idx <- sample(nrow(nz), size = floor(nrow(nz) * holdout_frac))
  mask     <- nz[mask_idx, , drop = FALSE]

  X_train <- X
  X_train[mask] <- 0L

  X_test_counts <- X[mask]          # held-out observed counts
  test_total    <- sum(X_test_counts)

  # Fit on training data
  fit <- fit_topic_model_gpu(Matrix::Matrix(X_train, sparse = TRUE), k = k,
                             n_iter_gpu = n_iter_gpu, n_iter_final = n_iter_final,
                             seed = seed)

  # Predicted Poisson rates at held-out positions
  lambda <- (fit$L %*% t(fit$F))[mask]
  lambda <- pmax(lambda, 1e-10)

  ll_test <- sum(dpois(X_test_counts, lambda, log = TRUE))
  perp    <- exp(-ll_test / test_total)

  list(perplexity = perp, fit = fit)
}

select_k <- function(counts_t, k_range = c(4, 6, 8, 10, 12), out_prefix) {
  cat("  Fitting models for K =", paste(k_range, collapse = ", "), "...\n")

  fits   <- list()
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

    # Cross-validated perplexity (always recomputed — cheap relative to fitting)
    cat("      CV perplexity k =", k, "...\n")
    cv <- cv_perplexity(counts_t, k = k)
    perps[i] <- cv$perplexity

    logliks[i] <- mean(loglik_multinom_topic_model(counts_t, fits[[as.character(k)]]))
  }

  cat("  Log-likelihoods and perplexity:\n")
  for (i in seq_along(k_range)) {
    cat("    K =", k_range[i], ": loglik =", round(logliks[i], 2),
        "  perplexity =", round(perps[i], 1), "\n")
  }

  # Plot log-likelihood
  ll_df <- data.frame(K = k_range, loglik = logliks)
  p_ll  <- ggplot(ll_df, aes(x = K, y = loglik)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of topics (K)", y = "Log-likelihood") +
    theme_publication()
  ggsave(paste0(out_prefix, "_loglik.svg"), p_ll,
         width = 85, height = 70, units = "mm", bg = "white")
  ggsave(paste0(out_prefix, "_loglik.png"), p_ll,
         width = 85, height = 70, units = "mm", dpi = 300, bg = "white")

  # Plot perplexity
  perp_df <- data.frame(K = k_range, perplexity = perps)
  p_perp  <- ggplot(perp_df, aes(x = K, y = perplexity)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of topics (K)", y = "Perplexity (lower = better)") +
    theme_publication()
  ggsave(paste0(out_prefix, "_perplexity.svg"), p_perp,
         width = 85, height = 70, units = "mm", bg = "white")
  ggsave(paste0(out_prefix, "_perplexity.png"), p_perp,
         width = 85, height = 70, units = "mm", dpi = 300, bg = "white")

  list(fits = fits, logliks = logliks, perplexities = perps)
}

# ---- Topic visualization helpers --------------------------------------------
plot_topic_structure <- function(fit, meta, condition_col, out_prefix) {
  theta <- poisson2multinom(fit)$L   # cells × topics
  df    <- as.data.frame(theta)
  df$cell      <- rownames(df)
  df$condition <- meta[[condition_col]][match(df$cell, rownames(meta))]

  topic_cols <- grep("^k", names(df), value = TRUE)
  df_long    <- tidyr::pivot_longer(df, cols = all_of(topic_cols),
                                    names_to = "topic", values_to = "weight")

  dom_topic  <- apply(theta, 1, which.max)
  cell_order <- rownames(theta)[order(df$condition, dom_topic)]
  df_long$cell <- factor(df_long$cell, levels = cell_order)

  p <- ggplot(df_long, aes(x = cell, y = weight, fill = topic)) +
    geom_col(width = 1) +
    scale_fill_viridis_d(option = "H") +
    facet_grid(~condition, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Topic weight", fill = "Topic") +
    theme_publication() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0.5, "mm"))

  ggsave(paste0(out_prefix, "_structure.svg"), p,
         width = 174, height = 70, units = "mm", bg = "white")
  ggsave(paste0(out_prefix, "_structure.png"), p,
         width = 174, height = 70, units = "mm", dpi = 300, bg = "white")
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

  df_heat <- as.data.frame(as.matrix(F_mat)) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "topic", values_to = "freq")

  p <- ggplot(df_heat, aes(x = topic, y = gene, fill = freq)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma") +
    labs(x = "Topic", y = NULL, fill = "Gene freq.") +
    theme_publication() +
    theme(axis.text.y = element_text(size = 7))

  ggsave(paste0(out_prefix, "_top_genes_heatmap.svg"), p,
         width = 100, height = 120, units = "mm", bg = "white")
  ggsave(paste0(out_prefix, "_top_genes_heatmap.png"), p,
         width = 100, height = 120, units = "mm", dpi = 300, bg = "white")
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
sc2      <- seu_list[["sc2"]]
counts2  <- GetAssayData(sc2, layer = "counts")
counts2f <- filter_genes(counts2)
counts2t <- t(counts2f)   # cells × genes (fastTopics convention)

cat("Running K selection for sc2...\n")
sc2_ks <- select_k(counts2t, k_range = c(3, 4, 5, 6, 8, 10, 12, 15),
                   out_prefix = file.path(out_dir, "sc2"))

best_k_sc2 <- 5
ckpt_sc2   <- file.path(out_dir, "fit_sc2_final.rds")
if (file.exists(ckpt_sc2)) {
  cat("\nLoading sc2 final model from checkpoint...\n")
  fit_sc2 <- readRDS(ckpt_sc2)
} else {
  cat("\nFitting final sc2 model K =", best_k_sc2, "...\n")
  fit_sc2 <- fit_topic_model_gpu(counts2t, k = best_k_sc2,
                                 n_iter_gpu = 150, n_iter_final = 100,
                                 seed = 42)
  saveRDS(fit_sc2, ckpt_sc2)
}

cat("Plotting sc2 topic structure...\n")
theta_sc2 <- plot_topic_structure(fit_sc2, sc2@meta.data, "condition",
                                  file.path(out_dir, "sc2"))
cat("Plotting sc2 top genes...\n")
plot_top_genes(fit_sc2, best_k_sc2, file.path(out_dir, "sc2"))

# =============================================================================
# SC78 — Control, Nutlin, OxStress, Quiescence x 2 reps (18,429 cells)
# =============================================================================
cat("\n=== SC78 ===\n")
sc78      <- seu_list[["sc78"]]
counts78  <- GetAssayData(sc78, layer = "counts")
counts78f <- filter_genes(counts78)
counts78t <- t(counts78f)

cat("Running K selection for sc78...\n")
sc78_ks <- select_k(counts78t, k_range = c(3, 4, 5, 6, 8, 10, 12, 15),
                    out_prefix = file.path(out_dir, "sc78"))

best_k_sc78 <- 5
ckpt_sc78   <- file.path(out_dir, "fit_sc78_final.rds")
if (file.exists(ckpt_sc78)) {
  cat("\nLoading sc78 final model from checkpoint...\n")
  fit_sc78 <- readRDS(ckpt_sc78)
} else {
  cat("\nFitting final sc78 model K =", best_k_sc78, "...\n")
  fit_sc78 <- fit_topic_model_gpu(counts78t, k = best_k_sc78,
                                  n_iter_gpu = 150, n_iter_final = 100,
                                  seed = 42)
  saveRDS(fit_sc78, ckpt_sc78)
}

cat("Plotting sc78 topic structure...\n")
theta_sc78 <- plot_topic_structure(fit_sc78, sc78@meta.data, "condition",
                                   file.path(out_dir, "sc78"))
cat("Plotting sc78 top genes...\n")
plot_top_genes(fit_sc78, best_k_sc78, file.path(out_dir, "sc78"))

cat("\n=== Done. Outputs in:", out_dir, "===\n")
