# lda_v2_functions.R — shared functions for v2 LDA pipeline
# Source this from per-sample run scripts.

library(Seurat)
library(fastTopics)
library(Matrix)
library(ggplot2)
library(cowplot)
library(dplyr)

repo_dir  <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling"
py_script <- file.path(repo_dir, "analysis/gpu_nmf.py")

# ---- GPU detection -----------------------------------------------------------
detect_r_gpu <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) return(FALSE)
  torch_dir  <- system.file(package = "torch")
  lantern_ok <- any(file.exists(
    file.path(torch_dir, "deps", "lantern.dll"),
    file.path(torch_dir, "deps", "liblantern.so"),
    file.path(torch_dir, "deps", "liblantern.dylib")
  ))
  if (!lantern_ok) return(FALSE)
  result <- tryCatch(
    system2(Sys.which("Rscript"),
            c("--no-save", "--no-environ", "-e",
              "library(torch,quietly=TRUE);q(status=if(cuda_is_available())0L else 1L)"),
            stdout = FALSE, stderr = FALSE, timeout = 30),
    error = function(e) 1L, warning = function(e) 1L)
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

init_gpu <- function() {
  use_r_gpu  <- detect_r_gpu()
  use_py_gpu <- if (!use_r_gpu) detect_py_gpu() else FALSE
  use_gpu    <- use_r_gpu || use_py_gpu
  if (use_r_gpu) {
    suppressPackageStartupMessages(library(torch))
    cat("GPU backend: R torch |", cuda_get_device_name(), "\n")
  } else if (use_py_gpu) {
    cat("GPU backend: Python torch |", cuda_device_name(), "\n")
  } else {
    cat("GPU backend: none — CPU\n")
  }
  list(use_r_gpu = use_r_gpu, use_py_gpu = use_py_gpu, use_gpu = use_gpu)
}

# ---- Variance-based gene selection -------------------------------------------
# within_disp : mean overdispersion (var/mean - 1) within each condition
# between_var : variance of condition-level means (after library-size norm)
# Union of top n_within and top n_between genes.
# No hard detection cutoffs; MT- and RP[SL] compete on equal footing.

filter_genes_variance <- function(counts, meta, condition_col = "condition",
                                   n_within = 7000, n_between = 5000) {
  cat("  Total genes:", nrow(counts), "\n")

  # ---- Exclude structural/artefact genes upfront ----------------------------
  # Ribosomal RNA (structural, not protein-coding; dominates libraries)
  is_rrna  <- grepl("^RNA[0-9]|^MT-RNR", rownames(counts))
  # Mitochondrial protein-coding: KEEP (mt dysfunction is senescence signal)
  # RP[SL] ribosomal proteins: KEEP (mTOR/ribosome biogenesis in senescence)
  # Very highly expressed lncRNAs that saturate library fractions
  is_lncrna_sat <- rownames(counts) %in% c("MALAT1", "NEAT1", "XIST", "TSIX",
                                            "KCNQ1OT1", "RMRP", "RPPH1")
  exclude  <- is_rrna | is_lncrna_sat
  counts   <- counts[!exclude, ]
  cat("  Excluded rRNA / saturating lncRNA:", sum(exclude), "genes\n")
  cat("  Remaining:", nrow(counts), "genes\n")

  # ---- Library-size normalise -----------------------------------------------
  lib_size <- Matrix::colSums(counts)
  med_lib  <- median(lib_size)
  norm     <- Matrix::t(Matrix::t(counts) / lib_size * med_lib)

  conditions <- sort(unique(na.omit(meta[[condition_col]])))
  cat("  Conditions:", paste(conditions, collapse = ", "), "\n")

  # ---- Within-condition normalised dispersion --------------------------------
  # Use (var/mean - 1) / sqrt(mean) — normalises for mean-variance relationship.
  # Prevents high-abundance genes from dominating on raw overdispersion alone.
  disp_sum <- rep(0, nrow(norm))
  n_cond   <- 0L
  for (cond in conditions) {
    cells <- which(meta[[condition_col]] == cond)
    if (length(cells) < 2) next
    sub      <- norm[, cells, drop = FALSE]
    m        <- Matrix::rowMeans(sub)
    ss       <- Matrix::rowSums(sub * sub)
    v        <- ss / length(cells) - m^2
    # Normalised excess dispersion: (VMR - 1) / sqrt(mean)
    disp_sum <- disp_sum + (v / pmax(m, 1e-10) - 1) / pmax(sqrt(m), 1)
    n_cond   <- n_cond + 1L
  }
  within_disp <- disp_sum / max(n_cond, 1L)

  # ---- Between-condition variance (mean-normalised) -------------------------
  cond_means  <- vapply(conditions, function(cond) {
    cells <- which(meta[[condition_col]] == cond)
    Matrix::rowMeans(norm[, cells, drop = FALSE])
  }, numeric(nrow(norm)))
  gene_grand_mean <- rowMeans(cond_means)
  raw_between     <- rowSums((cond_means - gene_grand_mean)^2) /
                     max(ncol(cond_means) - 1, 1)
  # Normalise by grand mean so fold-change-like signal dominates, not absolute level
  between_var <- raw_between / pmax(gene_grand_mean, 1e-10)

  keep_within  <- rank(-within_disp, ties.method = "min") <= n_within
  keep_between <- rank(-between_var, ties.method = "min") <= n_between
  keep         <- keep_within | keep_between

  cat("  Top within-dispersion:  ", sum(keep_within),  "genes\n")
  cat("  Top between-variance:   ", sum(keep_between), "genes\n")
  cat("  Union (kept):           ", sum(keep),         "genes\n")
  cat("  MT- protein genes kept: ", sum(keep & grepl("^MT-", rownames(counts))), "\n")
  cat("  RP[SL] genes kept:      ", sum(keep & grepl("^RP[SL]", rownames(counts))), "\n")

  counts[keep, ]
}

# ---- GPU warm-start ----------------------------------------------------------
gpu_warm_start <- function(X, k, n_iter_gpu, seed, gpu_state) {
  if (gpu_state$use_r_gpu) {
    n <- nrow(X); p <- ncol(X)
    set.seed(seed)
    W0 <- matrix(rgamma(n * k, 1, 1), n, k)
    H0 <- matrix(rgamma(k * p, 1, 1), k, p)
    device <- torch_device("cuda")
    V  <- torch_tensor(as.matrix(X), dtype = torch_float32(), device = device)
    W  <- torch_tensor(W0, dtype = torch_float32(), device = device)
    H  <- torch_tensor(H0, dtype = torch_float32(), device = device)
    rm(W0, H0); eps <- 1e-10
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
                              n_iter = as.integer(n_iter_gpu), seed = as.integer(seed))
    L_warm <- result$L; F_warm <- result$F
  }
  rownames(L_warm) <- rownames(X); rownames(F_warm) <- colnames(X)
  colnames(L_warm) <- colnames(F_warm) <- paste0("k", seq_len(k))
  list(L = L_warm, F = F_warm)
}

fit_topic_model_gpu <- function(X, k, n_iter_gpu = 100, n_iter_final = 50,
                                seed = 42, gpu_state) {
  if (!gpu_state$use_gpu) {
    cat("  [CPU] fit_topic_model k =", k, "\n")
    set.seed(seed)
    return(fit_topic_model(X, k = k,
                           numiter.main = n_iter_final + n_iter_gpu,
                           verbose = "none"))
  }
  cat("  [GPU] MU warm-start k =", k, "for", n_iter_gpu, "iters...\n")
  warm <- gpu_warm_start(X, k, n_iter_gpu, seed, gpu_state)
  cat("  [SCD] Refining k =", k, "for", n_iter_final, "iters...\n")
  set.seed(seed)
  fit0 <- init_poisson_nmf(X, F = warm$F, L = warm$L)
  fit_poisson_nmf(X, fit0 = fit0, numiter = n_iter_final,
                  method = "scd", verbose = "none")
}

# ---- Cross-validated perplexity ----------------------------------------------
cv_perplexity <- function(counts_t, k, n_iter_gpu = 100, n_iter_final = 50,
                          holdout_frac = 0.2, seed = 42, gpu_state) {
  set.seed(seed)
  X    <- as.matrix(counts_t)
  nz   <- which(X > 0, arr.ind = TRUE)
  mask <- nz[sample(nrow(nz), floor(nrow(nz) * holdout_frac)), , drop = FALSE]
  X_train       <- X; X_train[mask] <- 0L
  X_test_counts <- X[mask]
  fit    <- fit_topic_model_gpu(Matrix::Matrix(X_train, sparse = TRUE), k = k,
                                n_iter_gpu = n_iter_gpu, n_iter_final = n_iter_final,
                                seed = seed, gpu_state = gpu_state)
  lambda <- pmax((fit$L %*% t(fit$F))[mask], 1e-10)
  ll     <- sum(dpois(X_test_counts, lambda, log = TRUE))
  list(perplexity = exp(-ll / sum(X_test_counts)), fit = fit)
}

# ---- K selection -------------------------------------------------------------
select_k <- function(counts_t, k_range, out_prefix, gpu_state) {
  cat("  K range:", paste(k_range, collapse = ", "), "\n")
  fits    <- list()
  logliks <- numeric(length(k_range))
  perps   <- numeric(length(k_range))

  for (i in seq_along(k_range)) {
    k    <- k_range[i]
    ckpt <- file.path(out_prefix, paste0("fit_k", k, ".rds"))

    if (file.exists(ckpt)) {
      cat("    K =", k, "... loading checkpoint\n")
      fits[[as.character(k)]] <- readRDS(ckpt)
    } else {
      cat("    K =", k, "... fitting\n")
      fit <- fit_topic_model_gpu(counts_t, k = k,
                                 n_iter_gpu = 100, n_iter_final = 50,
                                 gpu_state = gpu_state)
      saveRDS(fit, ckpt)
      fits[[as.character(k)]] <- fit
    }

    cat("      CV perplexity k =", k, "...\n")
    cv       <- cv_perplexity(counts_t, k = k, gpu_state = gpu_state)
    perps[i] <- cv$perplexity
    logliks[i] <- mean(loglik_multinom_topic_model(counts_t, fits[[as.character(k)]]))
  }

  cat("\n  Results:\n")
  for (i in seq_along(k_range)) {
    cat("    K =", k_range[i], ": loglik =", round(logliks[i], 2),
        "  perplexity =", round(perps[i], 3), "\n")
  }

  theme_pub <- function(base_size = 11) {
    theme_cowplot(base_size) +
    theme(axis.line  = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text  = element_text(colour = "black", size = base_size),
          axis.title = element_text(colour = "black", size = base_size + 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background  = element_blank(),
          plot.title       = element_blank(),
          plot.margin      = margin(8, 8, 8, 8))
  }

  ll_df <- data.frame(K = k_range, loglik = logliks)
  p_ll  <- ggplot(ll_df, aes(x = K, y = loglik)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of topics (K)", y = "Log-likelihood") + theme_pub()
  ggsave(file.path(out_prefix, "loglik.png"),  p_ll,
         width = 85, height = 70, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(out_prefix, "loglik.svg"),  p_ll,
         width = 85, height = 70, units = "mm", bg = "white")

  perp_df <- data.frame(K = k_range, perplexity = perps)
  p_perp  <- ggplot(perp_df, aes(x = K, y = perplexity)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of topics (K)", y = "Perplexity") + theme_pub()
  ggsave(file.path(out_prefix, "perplexity.png"), p_perp,
         width = 85, height = 70, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(out_prefix, "perplexity.svg"), p_perp,
         width = 85, height = 70, units = "mm", bg = "white")

  saveRDS(list(k_range = k_range, logliks = logliks, perplexities = perps),
          file.path(out_prefix, "kselection.rds"))

  list(fits = fits, logliks = logliks, perplexities = perps)
}

# ---- Top genes ---------------------------------------------------------------
print_top_genes <- function(fit, sample_name) {
  k <- ncol(fit$F)
  cat("Top genes per topic (", sample_name, "K =", k, "):\n")
  for (i in seq_len(k)) {
    top10 <- names(head(sort(fit$F[, i], decreasing = TRUE), 10))
    cat("  Topic", i, ":", paste(top10, collapse = ", "), "\n")
  }
}
