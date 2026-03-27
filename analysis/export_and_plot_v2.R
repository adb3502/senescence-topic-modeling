# export_and_plot_v2.R
# Export gene weights + pathway tables and generate all tidyplots figures
# for the winning K per sample: sc2/sc78 K=4, sc5 K=6
#
# Outputs per sample directory:
#   gene_weights_k{K}.csv
#   top_pathways_k{K}.csv
#   perplexity_curve.png/svg
#   bio_coherence_panel.png/svg
#   gene_heatmap_k{K}.png/svg
#   structure_plot_k{K}.png/svg
#   pathway_heatmap_k{K}.png/svg
# sc5 only:
#   timecourse_trajectory_k10.png/svg

library(ggplot2)
library(tidyplots)
library(tidyr)
library(dplyr)
library(fastTopics)
library(Seurat)

v2_dir  <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"
qc_dir  <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"

cat("Loading Seurat objects...\n")
seu_list <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))

# ---------------------------------------------------------------------------
# Colour palettes
# ---------------------------------------------------------------------------
# Okabe-Ito 8-colour + 4 extras for K up to 12 (all colour-blind friendly)
oi_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
            "#CC79A7", "#000000", "#88CCEE", "#44AA99", "#332288", "#DDCC77")

metric_pal <- c(combined = "#2C3E7A", coherence = "#2196A0", specificity = "#E07B39")

# ---------------------------------------------------------------------------
# Helper: shorten Reactome pathway names
# ---------------------------------------------------------------------------
clean_pathway <- function(x, max_chars = 45) {
  x <- sub("^REACTOME_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 1), "\u2026"), x)
}

# ---------------------------------------------------------------------------
# Core function: exports + plots for one sample at best_k
# ---------------------------------------------------------------------------
process_sample <- function(sample_name, best_k, conditions_order = NULL) {
  model_dir   <- file.path(v2_dir, sample_name)
  K           <- best_k
  topic_lbls  <- paste0("Topic ", seq_len(K))
  topic_cols  <- oi_pal[seq_len(K)]

  cat("\n===", sample_name, "K =", K, "===\n")

  fit  <- readRDS(file.path(model_dir, paste0("fit_k", K, ".rds")))
  bio  <- readRDS(file.path(model_dir, "bio_scores.rds"))
  ksel <- readRDS(file.path(model_dir, "kselection.rds"))

  colnames(fit$F) <- topic_lbls
  colnames(fit$L) <- topic_lbls

  # ---- 1. CSV EXPORTS -------------------------------------------------------

  # Gene weights (all genes)
  gw_df <- as.data.frame(fit$F)
  gw_df$gene <- rownames(fit$F)
  gw_df <- gw_df[, c("gene", topic_lbls)]
  write.csv(gw_df, file.path(model_dir, paste0("gene_weights_k", K, ".csv")),
            row.names = FALSE)
  cat("  Exported gene_weights_k", K, ".csv (", nrow(gw_df), "genes)\n", sep = "")

  # Top 10 pathways per topic
  scores_k  <- bio$scores[[as.character(K)]]
  path_rows <- lapply(seq_len(K), function(i) {
    res <- scores_k$nes_list[[i]]
    pos <- res[!is.na(res$NES) & res$NES > 0, ]
    pos <- pos[order(pos$padj), ]
    head(pos[, c("pathway", "NES", "padj", "size")], 10) |>
      mutate(topic = topic_lbls[i])
  })
  path_tbl <- bind_rows(path_rows)
  write.csv(path_tbl, file.path(model_dir, paste0("top_pathways_k", K, ".csv")),
            row.names = FALSE)
  cat("  Exported top_pathways_k", K, ".csv\n", sep = "")

  # ---- 2. PERPLEXITY CURVE --------------------------------------------------
  perp_df <- data.frame(K_val = ksel$k_range, perplexity = ksel$perplexities)

  p_perp <- perp_df |>
    tidyplot(x = K_val, y = perplexity) |>
    add_mean_line(linewidth = 0.9, color = "#2C3E7A") |>
    add_mean_dot(size = 2.8, color = "#2C3E7A") |>
    theme_tidyplot(fontsize = 9) +
    ggplot2::scale_x_continuous(breaks = ksel$k_range) +
    ggplot2::labs(x = "Number of topics (K)", y = "CV perplexity")

  save_plot(p_perp, filename = file.path(model_dir, "perplexity_curve.png"),
            width = 90, height = 68, units = "mm")
  save_plot(p_perp, filename = file.path(model_dir, "perplexity_curve.svg"),
            width = 90, height = 68, units = "mm")

  # ---- 3. BIO COHERENCE PANEL -----------------------------------------------
  bio_long <- bio$df |>
    pivot_longer(c(coherence, specificity, combined),
                 names_to = "metric", values_to = "value") |>
    mutate(metric = factor(metric, levels = c("combined", "coherence", "specificity")))

  p_bio <- bio_long |>
    tidyplot(x = K, y = value, color = metric) |>
    add_mean_line(linewidth = 0.9) |>
    add_mean_dot(size = 2.8) |>
    adjust_colors(new_colors = unname(metric_pal)) |>
    theme_tidyplot(fontsize = 9) +
    ggplot2::scale_x_continuous(breaks = bio$df$K) +
    ggplot2::labs(x = "Number of topics (K)", y = "Score", color = NULL)

  save_plot(p_bio, filename = file.path(model_dir, "bio_coherence_panel.png"),
            width = 110, height = 75, units = "mm")
  save_plot(p_bio, filename = file.path(model_dir, "bio_coherence_panel.svg"),
            width = 110, height = 75, units = "mm")

  # ---- 4. GENE WEIGHT HEATMAP -----------------------------------------------
  # Top 10 per topic (unique union) — 10 keeps ~30-35 unique genes, readable at 8mm/row
  top_genes <- unique(unlist(lapply(seq_len(K), function(i) {
    names(sort(fit$F[, i], decreasing = TRUE))[1:10]
  })))
  F_sub    <- fit$F[top_genes, , drop = FALSE]
  hm_long  <- as.data.frame(F_sub) |>
    mutate(gene = rownames(F_sub)) |>
    pivot_longer(-gene, names_to = "topic", values_to = "weight") |>
    mutate(topic = factor(topic, levels = topic_lbls),
           gene  = factor(gene,  levels = rev(top_genes)))

  p_hm <- hm_long |>
    tidyplot(x = topic, y = gene, color = weight) |>
    add_heatmap(scale = "row") |>
    theme_tidyplot(fontsize = 7) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Row z-score")

  hm_w <- 55 + K * 20
  hm_h <- 25 + length(top_genes) * 8   # 8mm per gene row — no overlap at 6pt labels
  save_plot(p_hm, filename = file.path(model_dir, paste0("gene_heatmap_k", K, ".png")),
            width = hm_w, height = hm_h, units = "mm", dpi = 300)
  save_plot(p_hm, filename = file.path(model_dir, paste0("gene_heatmap_k", K, ".svg")),
            width = hm_w, height = hm_h, units = "mm")

  # ---- 5. STRUCTURE PLOT (mean topic proportion per condition) ---------------
  seu_meta  <- seu_list[[sample_name]]@meta.data
  # Normalize L to true topic proportions (fastTopics L rows are unnormalized loadings)
  L_norm    <- fit$L / rowSums(fit$L)
  colnames(L_norm) <- topic_lbls
  L_df      <- as.data.frame(L_norm)
  L_df$cell_id   <- rownames(fit$L)
  L_df$condition <- seu_meta[L_df$cell_id, "condition"]

  if (!is.null(conditions_order)) {
    L_df$condition <- factor(L_df$condition, levels = conditions_order)
  }

  L_cond <- L_df |>
    group_by(condition) |>
    summarise(across(all_of(topic_lbls), mean), .groups = "drop")

  L_cond_long <- L_cond |>
    pivot_longer(-condition, names_to = "topic", values_to = "proportion") |>
    mutate(topic = factor(topic, levels = topic_lbls))

  p_struct <- L_cond_long |>
    tidyplot(x = condition, y = proportion, color = topic) |>
    add_barstack_relative() |>
    adjust_colors(new_colors = topic_cols) |>
    theme_tidyplot(fontsize = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1)) +
    ggplot2::labs(x = NULL, y = "Mean topic proportion", fill = NULL, color = NULL)

  n_conds <- length(unique(L_df$condition))
  save_plot(p_struct,
            filename = file.path(model_dir, paste0("structure_plot_k", K, ".png")),
            width = 55 + n_conds * 18, height = 85, units = "mm", dpi = 300)
  save_plot(p_struct,
            filename = file.path(model_dir, paste0("structure_plot_k", K, ".svg")),
            width = 55 + n_conds * 18, height = 85, units = "mm")

  # ---- 5b. INDIVIDUAL-CELL STRUCTURE PLOT (plain ggplot2) --------------------
  max_per_cond <- 700L
  set.seed(42)
  L_cells <- L_df |>
    group_by(condition) |>
    group_modify(~ slice_sample(.x, n = min(max_per_cond, nrow(.x)))) |>
    ungroup()

  L_cells$dominant <- topic_lbls[max.col(as.matrix(L_cells[, topic_lbls, drop = FALSE]))]
  L_cells <- L_cells |>
    arrange(condition, dominant, desc(.data[[topic_lbls[1]]]))

  # Integer x positions — avoids discrete-factor whitespace issues in ggplot2
  L_cells$x <- seq_len(nrow(L_cells))

  L_cells_long <- L_cells |>
    select(x, condition, all_of(topic_lbls)) |>
    pivot_longer(all_of(topic_lbls), names_to = "topic", values_to = "proportion") |>
    mutate(topic = factor(topic, levels = topic_lbls))

  # Condition dividers and label midpoints
  cond_counts <- table(L_cells$condition)[levels(factor(L_cells$condition,
                                                         levels = if (!is.null(conditions_order))
                                                           conditions_order else unique(L_cells$condition)))]
  cond_counts <- cond_counts[cond_counts > 0]
  cond_ends   <- cumsum(cond_counts)
  cond_starts <- c(1, head(cond_ends, -1) + 1)
  cond_mids   <- (cond_starts + cond_ends) / 2
  dividers    <- cond_ends[-length(cond_ends)] + 0.5

  p_cells <- ggplot(L_cells_long, aes(x = x, y = proportion, fill = topic)) +
    geom_col(width = 1, position = "stack") +
    scale_fill_manual(values = setNames(topic_cols, topic_lbls)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, nrow(L_cells) + 0.5)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    geom_vline(xintercept = dividers, colour = "white", linewidth = 0.8) +
    annotate("text", x = cond_mids, y = 1.08, label = names(cond_counts),
             size = 2.8, fontface = "bold", colour = "grey15", hjust = 0.5, vjust = 0) +
    annotate("text", x = cond_mids, y = 1.03,
             label = paste0("n=", as.integer(cond_counts)),
             size = 2.2, colour = "grey40", hjust = 0.5, vjust = 0) +
    theme_classic(base_size = 8) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x  = element_blank(),
          legend.key.size = unit(3.5, "mm"),
          plot.margin  = margin(14, 4, 4, 4, "mm")) +
    labs(x = NULL, y = "Topic proportion", fill = NULL)

  cell_w <- max(180, nrow(L_cells) * 0.18)
  ggsave(file.path(model_dir, paste0("structure_plot_cells_k", K, ".png")),
         p_cells, width = cell_w, height = 70, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(model_dir, paste0("structure_plot_cells_k", K, ".svg")),
         p_cells, width = cell_w, height = 70, units = "mm", bg = "white")

  # ---- 6. PATHWAY HEATMAP (top 5 per topic, NES fill) -----------------------
  ph_rows <- lapply(seq_len(K), function(i) {
    res <- scores_k$nes_list[[i]]
    pos <- res[!is.na(res$NES) & res$NES > 0, ]
    pos <- pos[order(pos$padj), ]
    top5 <- head(pos, 5)
    data.frame(
      topic          = topic_lbls[i],
      pathway        = clean_pathway(top5$pathway),
      NES            = top5$NES,
      neg_log10_padj = -log10(pmax(top5$padj, 1e-10)),
      stringsAsFactors = FALSE
    )
  })
  ph_df <- bind_rows(ph_rows) |>
    mutate(topic = factor(topic, levels = topic_lbls))

  # Keep pathway order: rank by max NES across topics
  path_order <- ph_df |>
    group_by(pathway) |>
    summarise(max_nes = max(NES)) |>
    arrange(max_nes) |>
    pull(pathway)
  ph_df$pathway <- factor(ph_df$pathway, levels = path_order)

  p_path <- ph_df |>
    tidyplot(x = topic, y = pathway, color = neg_log10_padj) |>
    add_heatmap(scale = "none") |>
    adjust_colors(new_colors = colors_continuous_bluepinkyellow) |>
    theme_tidyplot(fontsize = 7.5) +
    ggplot2::labs(x = NULL, y = NULL, fill = "-log\u2081\u2080(FDR)")

  ph_h <- 40 + nrow(ph_df) * 5
  save_plot(p_path,
            filename = file.path(model_dir, paste0("pathway_heatmap_k", K, ".png")),
            width = 140, height = ph_h, units = "mm", dpi = 300)
  save_plot(p_path,
            filename = file.path(model_dir, paste0("pathway_heatmap_k", K, ".svg")),
            width = 140, height = ph_h, units = "mm")

  cat("  All plots saved to", model_dir, "\n")
  invisible(list(fit = fit, bio = bio, ksel = ksel, L_cond = L_cond, path_tbl = path_tbl))
}

# ---------------------------------------------------------------------------
# Process sc2 and sc78 (K=4)
# ---------------------------------------------------------------------------
sc2_conds  <- c("Nutlin", "OxStress", "Quiescence")
sc78_conds <- c("sc7 Control", "sc7 Nutlin", "sc7 OxStress", "sc7 Quiescence",
                "sc8 Control", "sc8 Nutlin", "sc8 OxStress", "sc8 Quiescence")

res_sc2  <- process_sample("sc2",  best_k = 4, conditions_order = sc2_conds)
res_sc78 <- process_sample("sc78", best_k = 9, conditions_order = sc78_conds)

# ---------------------------------------------------------------------------
# Process sc5 (K=10, timecourse)
# ---------------------------------------------------------------------------
tp_order <- c("Control", "Sen 4hrs", "Sen D1", "Sen D2", "Sen D3", "Sen D4", "Sen D7")
res_sc5  <- process_sample("sc5",  best_k = 10, conditions_order = tp_order)

# ---- sc5 timecourse trajectory (additional plot) --------------------------
cat("\nGenerating sc5 timecourse trajectory...\n")

sc5_dir   <- file.path(v2_dir, "sc5")
fit5      <- res_sc5$fit
K5        <- 10
topic_lbls5 <- paste0("Topic ", seq_len(K5))
topic_cols5 <- oi_pal[seq_len(K5)]

sc5_meta  <- seu_list[["sc5"]]@meta.data
L5_norm   <- fit5$L / rowSums(fit5$L)   # normalize to true proportions
L5_df     <- as.data.frame(L5_norm)
colnames(L5_df) <- topic_lbls5
L5_df$cell_id   <- rownames(fit5$L)
L5_df$timepoint <- factor(sc5_meta[L5_df$cell_id, "condition"], levels = tp_order)
L5_df <- L5_df[!is.na(L5_df$timepoint), ]   # drop cells with unassigned condition

# Map timepoints to numeric days for x-axis
tp_days <- c("Control" = 0, "Sen 4hrs" = 0.17, "Sen D1" = 1,
             "Sen D2" = 2, "Sen D3" = 3, "Sen D4" = 4, "Sen D7" = 7)
L5_df$day <- tp_days[as.character(L5_df$timepoint)]

L5_long <- L5_df |>
  select(cell_id, day, timepoint, all_of(topic_lbls5)) |>
  pivot_longer(all_of(topic_lbls5), names_to = "topic", values_to = "proportion") |>
  mutate(topic = factor(topic, levels = topic_lbls5))

p_tc <- L5_long |>
  tidyplot(x = day, y = proportion, color = topic) |>
  add_mean_line(linewidth = 0.9) |>
  add_sem_ribbon(alpha = 0.18) |>
  add_mean_dot(size = 2.2) |>
  adjust_colors(new_colors = topic_cols5) |>
  theme_tidyplot(fontsize = 9) +
  ggplot2::scale_x_continuous(
    breaks = c(0, 0.17, 1, 2, 3, 4, 7),
    labels = c("Ctrl", "4h", "D1", "D2", "D3", "D4", "D7")
  ) +
  ggplot2::labs(x = "Day post H\u2082O\u2082", y = "Mean topic proportion",
                color = NULL, fill = NULL)

save_plot(p_tc, filename = file.path(sc5_dir, "timecourse_trajectory_k10.png"),
          width = 130, height = 85, units = "mm", dpi = 300)
save_plot(p_tc, filename = file.path(sc5_dir, "timecourse_trajectory_k10.svg"),
          width = 130, height = 85, units = "mm")
cat("  Saved timecourse_trajectory_k10.png/svg\n")

# ---------------------------------------------------------------------------
# Individual-cell structure plot for sc5 (sorted by timepoint, plain ggplot2)
# ---------------------------------------------------------------------------
cat("Generating sc5 individual-cell structure plot...\n")

set.seed(42)
L5_sampled <- L5_df |>
  group_by(timepoint) |>
  group_modify(~ slice_sample(.x, n = min(700L, nrow(.x)))) |>
  ungroup()

L5_sampled$dominant <- topic_lbls5[max.col(as.matrix(L5_sampled[, topic_lbls5, drop = FALSE]))]
L5_sampled <- L5_sampled |>
  arrange(timepoint, dominant, desc(.data[[topic_lbls5[1]]]))
L5_sampled$x <- seq_len(nrow(L5_sampled))

L5_samp_long <- L5_sampled |>
  select(x, timepoint, all_of(topic_lbls5)) |>
  pivot_longer(all_of(topic_lbls5), names_to = "topic", values_to = "proportion") |>
  mutate(topic = factor(topic, levels = topic_lbls5))

tp_counts  <- table(L5_sampled$timepoint)[tp_order]
tp_counts  <- tp_counts[tp_counts > 0]
tp_ends    <- cumsum(tp_counts)
tp_starts  <- c(1, head(tp_ends, -1) + 1)
tp_mids    <- (tp_starts + tp_ends) / 2
tp_divs    <- tp_ends[-length(tp_ends)] + 0.5

p_str5 <- ggplot(L5_samp_long, aes(x = x, y = proportion, fill = topic)) +
  geom_col(width = 1, position = "stack") +
  scale_fill_manual(values = setNames(topic_cols5, topic_lbls5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, nrow(L5_sampled) + 0.5)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  geom_vline(xintercept = tp_divs, colour = "white", linewidth = 0.8) +
  annotate("text", x = tp_mids, y = 1.08, label = names(tp_counts),
           size = 2.8, fontface = "bold", colour = "grey15", hjust = 0.5, vjust = 0) +
  annotate("text", x = tp_mids, y = 1.03,
           label = paste0("n=", as.integer(tp_counts)),
           size = 2.2, colour = "grey40", hjust = 0.5, vjust = 0) +
  theme_classic(base_size = 8) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        legend.key.size = unit(3.5, "mm"),
        plot.margin  = margin(14, 4, 4, 4, "mm")) +
  labs(x = NULL, y = "Topic proportion", fill = NULL)

str5_w <- max(180, nrow(L5_sampled) * 0.18)
ggsave(file.path(sc5_dir, "structure_plot_cells_k10.png"),
       p_str5, width = str5_w, height = 70, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(sc5_dir, "structure_plot_cells_k10.svg"),
       p_str5, width = str5_w, height = 70, units = "mm", bg = "white")

# ---------------------------------------------------------------------------
# Helper: build one row panel of the 2-row structure plot
# conds_subset: character vector of condition labels to include (in order)
# L_cells: full individual-cell data frame (already sampled, sorted, x assigned)
# topic_lbls / topic_cols: labels and colours
# show_legend: include legend on this row?
# ---------------------------------------------------------------------------
make_struct_row <- function(L_cells, conds_subset, topic_lbls, topic_cols,
                             show_legend = FALSE, x_max = NULL) {
  row_cells <- L_cells[L_cells$condition %in% conds_subset, ]
  row_cells$condition <- factor(row_cells$condition, levels = conds_subset)
  row_cells <- row_cells[order(row_cells$condition, row_cells$dominant), ]
  row_cells$x <- seq_len(nrow(row_cells))

  row_long <- row_cells |>
    select(x, condition, all_of(topic_lbls)) |>
    pivot_longer(all_of(topic_lbls), names_to = "topic", values_to = "proportion") |>
    mutate(topic = factor(topic, levels = topic_lbls))

  cond_counts <- table(row_cells$condition)[conds_subset]
  cond_counts <- cond_counts[cond_counts > 0]
  cond_ends   <- cumsum(cond_counts)
  cond_starts <- c(1, head(cond_ends, -1) + 1)
  cond_mids   <- (cond_starts + cond_ends) / 2
  dividers    <- cond_ends[-length(cond_ends)] + 0.5

  xlim_max <- if (!is.null(x_max)) x_max + 0.5 else nrow(row_cells) + 0.5

  p <- ggplot(row_long, aes(x = x, y = proportion, fill = topic)) +
    geom_col(width = 1, position = "stack") +
    scale_fill_manual(values = setNames(topic_cols, topic_lbls)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, xlim_max)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    geom_vline(xintercept = dividers, colour = "white", linewidth = 0.8) +
    annotate("text", x = cond_mids, y = 1.08, label = names(cond_counts),
             size = 2.8, fontface = "bold", colour = "grey15", hjust = 0.5, vjust = 0) +
    annotate("text", x = cond_mids, y = 1.03,
             label = paste0("n=", as.integer(cond_counts)),
             size = 2.2, colour = "grey40", hjust = 0.5, vjust = 0) +
    theme_classic(base_size = 8) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x  = element_blank(),
          legend.key.size = unit(3.5, "mm"),
          plot.margin  = margin(14, 4, 4, 4, "mm")) +
    labs(x = NULL, y = "Topic proportion", fill = NULL)

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

# ---------------------------------------------------------------------------
# 2-row structure plot for sc78 (K=9)
# Row 1: sc7 conditions; Row 2: sc8 conditions
# ---------------------------------------------------------------------------
cat("\nGenerating sc78 2-row individual-cell structure plot...\n")

sc78_dir    <- file.path(v2_dir, "sc78")
fit78       <- res_sc78$fit
K78         <- 9L
topic_lbls78 <- paste0("Topic ", seq_len(K78))
topic_cols78 <- oi_pal[seq_len(K78)]

sc78_meta <- seu_list[["sc78"]]@meta.data
L78_norm  <- fit78$L / rowSums(fit78$L)
colnames(L78_norm) <- topic_lbls78
L78_df    <- as.data.frame(L78_norm)
L78_df$cell_id   <- rownames(fit78$L)
L78_df$condition <- factor(sc78_meta[L78_df$cell_id, "condition"], levels = sc78_conds)
L78_df <- L78_df[!is.na(L78_df$condition), ]

set.seed(42)
L78_cells <- L78_df |>
  group_by(condition) |>
  group_modify(~ slice_sample(.x, n = min(700L, nrow(.x)))) |>
  ungroup()

L78_cells$dominant <- topic_lbls78[max.col(as.matrix(L78_cells[, topic_lbls78, drop = FALSE]))]

sc7_conds <- c("sc7 Control", "sc7 Nutlin", "sc7 OxStress", "sc7 Quiescence")
sc8_conds <- c("sc8 Control", "sc8 Nutlin", "sc8 OxStress", "sc8 Quiescence")

p78_top <- make_struct_row(L78_cells, sc7_conds, topic_lbls78, topic_cols78, show_legend = FALSE)
p78_bot <- make_struct_row(L78_cells, sc8_conds, topic_lbls78, topic_cols78, show_legend = TRUE)

# Shared legend from bottom row; combine 2 rows + legend column
library(cowplot)
legend78 <- cowplot::get_legend(p78_bot)
p78_bot_noleg <- p78_bot + theme(legend.position = "none")

p78_rows <- cowplot::plot_grid(p78_top, p78_bot_noleg, nrow = 2, align = "v", axis = "lr")
p78_2row <- cowplot::plot_grid(p78_rows, legend78, ncol = 2, rel_widths = c(1, 0.12))

n78_cells <- sum(table(L78_cells$condition))
row_w78   <- max(180, n78_cells / 2 * 0.18)

ggsave(file.path(sc78_dir, "structure_plot_cells_k9_2row.png"),
       p78_2row, width = row_w78, height = 130, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(sc78_dir, "structure_plot_cells_k9_2row.svg"),
       p78_2row, width = row_w78, height = 130, units = "mm", bg = "white")
cat("  Saved structure_plot_cells_k9_2row.png/svg\n")

# ---------------------------------------------------------------------------
# 2-row structure plot for sc5 (K=10)
# Row 1: Control, Sen 4hrs, Sen D1; Row 2: Sen D2, Sen D3, Sen D4, Sen D7
# ---------------------------------------------------------------------------
cat("\nGenerating sc5 2-row individual-cell structure plot...\n")

sc5_conds_top <- c("Control", "Sen 4hrs", "Sen D1")
sc5_conds_bot <- c("Sen D2", "Sen D3", "Sen D4", "Sen D7")

# Re-use L5_sampled from the sc5 section above (already normalised + sampled)
# Need to add 'condition' column matching the tp_order labels
L5_sampled$condition <- as.character(L5_sampled$timepoint)
L5_sampled$dominant  <- topic_lbls5[max.col(as.matrix(L5_sampled[, topic_lbls5, drop = FALSE]))]

n5_top  <- sum(L5_sampled$condition %in% sc5_conds_top)
n5_bot  <- sum(L5_sampled$condition %in% sc5_conds_bot)
x5_max  <- max(n5_top, n5_bot)

p5_top <- make_struct_row(L5_sampled, sc5_conds_top, topic_lbls5, topic_cols5, show_legend = FALSE, x_max = x5_max)
p5_bot <- make_struct_row(L5_sampled, sc5_conds_bot, topic_lbls5, topic_cols5, show_legend = TRUE,  x_max = x5_max)

legend5     <- cowplot::get_legend(p5_bot)
p5_bot_noleg <- p5_bot + theme(legend.position = "none")

p5_rows <- cowplot::plot_grid(p5_top, p5_bot_noleg, nrow = 2, align = "v", axis = "lr")
p5_2row <- cowplot::plot_grid(p5_rows, legend5, ncol = 2, rel_widths = c(1, 0.12))

row_w5   <- max(180, x5_max * 0.18)

ggsave(file.path(sc5_dir, "structure_plot_cells_k10_2row.png"),
       p5_2row, width = row_w5, height = 130, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(sc5_dir, "structure_plot_cells_k10_2row.svg"),
       p5_2row, width = row_w5, height = 130, units = "mm", bg = "white")
cat("  Saved structure_plot_cells_k10_2row.png/svg\n")

cat("\nAll exports and plots complete.\n")
