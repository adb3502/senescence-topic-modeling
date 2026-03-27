# structure_plots_consensus.R
# Individual-cell structure plots with consensus topic labels A-L and
# consistent colors across sc2, sc78, sc5.
# Directly adapts the working export_and_plot_v2.R cell plot code —
# only topic labels and colors are changed.

library(ggplot2)
library(dplyr)
library(tidyr)
library(fastTopics)
library(Seurat)
library(cowplot)

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"
qc_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"

cat("Loading Seurat objects...\n")
seu_list <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))

# ---- Consensus color palette (by dominance rank) ----------------------------
consensus_cols <- c(
  H = "#D55E00",   # rank 1
  A = "#0072B2",   # rank 2
  B = "#009E73",   # rank 3
  F = "#E69F00",   # rank 4
  D = "#56B4E9",   # rank 5
  E = "#CC79A7",   # rank 6
  C = "#44AA99",   # rank 7
  G = "#882255",   # rank 8
  I = "#DDCC77",   # rank 9
  K = "#88CCEE",   # rank 10
  J = "#332288",   # rank 11
  L = "#999933"    # rank 12
)

# ---- Topic-to-letter maps (Topic N -> letter) --------------------------------
sc2_map  <- c("Topic 1"="F", "Topic 2"="H", "Topic 3"="A", "Topic 4"="B")
sc78_map <- c("Topic 1"="G", "Topic 2"="D", "Topic 3"="E", "Topic 4"="F",
              "Topic 5"="C", "Topic 6"="H", "Topic 7"="B", "Topic 8"="I", "Topic 9"="A")
sc5_map  <- c("Topic 1"="A", "Topic 2"="B", "Topic 3"="J", "Topic 4"="F",
              "Topic 5"="G", "Topic 6"="D", "Topic 7"="K", "Topic 8"="L",
              "Topic 9"="C", "Topic 10"="E")

# ---- Core plotting function (identical logic to export_and_plot_v2.R 5b) ----
plot_consensus_cells <- function(sample_name, fit, topic_map, conditions_order,
                                  seu_list, x_max = NULL, show_legend = TRUE) {
  K      <- ncol(fit$L)
  t_lbls <- paste0("Topic ", seq_len(K))

  # Normalize
  L_norm <- fit$L / rowSums(fit$L)
  colnames(L_norm) <- t_lbls

  L_df <- as.data.frame(L_norm)
  L_df$cell_id   <- rownames(fit$L)
  L_df$condition <- factor(seu_list[[sample_name]]@meta.data[L_df$cell_id, "condition"],
                            levels = conditions_order)
  L_df <- L_df[!is.na(L_df$condition), ]

  # Rename topic columns to consensus letters
  for (tl in t_lbls) {
    letter <- topic_map[tl]
    L_df[[letter]] <- L_df[[tl]]
    L_df[[tl]]     <- NULL
  }
  letter_lbls <- unname(topic_map)   # letters in original topic-number order

  # Colors and display labels for this sample's letters
  topic_cols  <- consensus_cols[letter_lbls]
  topic_disp  <- paste0(letter_lbls, " (t", seq_len(K), ")")  # e.g. "F (t1)"

  # Downsample
  set.seed(42)
  L_cells <- L_df |>
    group_by(condition) |>
    group_modify(~ slice_sample(.x, n = min(700L, nrow(.x)))) |>
    ungroup()

  # Sort: dominant letter (in original topic order = letter_lbls order),
  # then descending weight of the first topic — exactly as in export_and_plot_v2.R
  L_cells$dominant      <- letter_lbls[max.col(as.matrix(L_cells[, letter_lbls, drop = FALSE]))]
  L_cells$dominant_rank <- match(L_cells$dominant, letter_lbls)  # rank by original topic order
  L_cells <- L_cells |>
    arrange(condition, dominant_rank, desc(.data[[letter_lbls[1]]]))

  L_cells$x <- seq_len(nrow(L_cells))

  L_cells_long <- L_cells |>
    select(x, condition, all_of(letter_lbls)) |>
    pivot_longer(all_of(letter_lbls), names_to = "topic", values_to = "proportion") |>
    mutate(topic = factor(topic, levels = letter_lbls))

  cond_counts <- table(L_cells$condition)[conditions_order]
  cond_counts <- cond_counts[cond_counts > 0]
  cond_ends   <- cumsum(cond_counts)
  cond_starts <- c(1, head(cond_ends, -1) + 1)
  cond_mids   <- (cond_starts + cond_ends) / 2
  dividers    <- cond_ends[-length(cond_ends)] + 0.5

  xlim_max <- if (!is.null(x_max)) x_max + 0.5 else nrow(L_cells) + 0.5

  p <- ggplot(L_cells_long, aes(x = x, y = proportion, fill = topic)) +
    geom_col(width = 1, position = "stack") +
    scale_fill_manual(values = setNames(topic_cols, letter_lbls),
                      labels = setNames(topic_disp, letter_lbls)) +
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
    theme(axis.text.x   = element_blank(),
          axis.ticks.x  = element_blank(),
          axis.line.x   = element_blank(),
          legend.key.size = unit(3.5, "mm"),
          plot.margin   = margin(14, 4, 4, 4, "mm")) +
    labs(x = NULL, y = "Topic proportion", fill = NULL)

  if (!show_legend) p <- p + theme(legend.position = "none")

  list(p = p, L_cells = L_cells, letter_lbls = letter_lbls, n_cells = nrow(L_cells))
}

# ============================================================
# sc2 — single row
# ============================================================
cat("\nsc2...\n")
sc2_conds <- c("Nutlin", "OxStress", "Quiescence")
fit2 <- readRDS(file.path(v2_dir, "sc2", "fit_k4.rds"))

r2 <- plot_consensus_cells("sc2", fit2, sc2_map, sc2_conds, seu_list)
w2 <- max(180, r2$n_cells * 0.18)
ggsave(file.path(v2_dir, "sc2", "structure_consensus.png"),
       r2$p, width = w2, height = 70, units = "mm", dpi = 300, bg = "white")
cat("  Saved sc2/structure_consensus.png\n")

# ============================================================
# sc78 — 2-row (sc7 top / sc8 bottom)
# ============================================================
cat("\nsc78...\n")
sc78_conds <- c("sc7 Control", "sc7 Nutlin", "sc7 OxStress", "sc7 Quiescence",
                "sc8 Control", "sc8 Nutlin", "sc8 OxStress", "sc8 Quiescence")
sc7_conds  <- sc78_conds[1:4]
sc8_conds  <- sc78_conds[5:8]
fit78 <- readRDS(file.path(v2_dir, "sc78", "fit_k9.rds"))

# Generate both rows using full condition order so x_max is consistent
r78_top <- plot_consensus_cells("sc78", fit78, sc78_map, sc78_conds, seu_list,
                                 show_legend = FALSE)
# count cells per row to set shared x_max
n78_top <- sum(r78_top$L_cells$condition %in% sc7_conds)
n78_bot <- sum(r78_top$L_cells$condition %in% sc8_conds)
x78_max <- max(n78_top, n78_bot)

r78_top2 <- plot_consensus_cells("sc78", fit78, sc78_map, sc7_conds, seu_list,
                                  x_max = x78_max, show_legend = FALSE)
r78_bot2 <- plot_consensus_cells("sc78", fit78, sc78_map, sc8_conds, seu_list,
                                  x_max = x78_max, show_legend = TRUE)

leg78    <- cowplot::get_legend(r78_bot2$p)
p78_rows <- cowplot::plot_grid(r78_top2$p, r78_bot2$p + theme(legend.position = "none"),
                                nrow = 2, align = "v", axis = "lr")
p78_2row <- cowplot::plot_grid(p78_rows, leg78, ncol = 2, rel_widths = c(1, 0.12))

w78 <- max(180, x78_max * 0.18)
ggsave(file.path(v2_dir, "sc78", "structure_consensus_2row.png"),
       p78_2row, width = w78, height = 130, units = "mm", dpi = 300, bg = "white")
cat("  Saved sc78/structure_consensus_2row.png\n")

# ============================================================
# sc5 — 2-row (Control-D1 top / D2-D7 bottom)
# ============================================================
cat("\nsc5...\n")
sc5_conds <- c("Control", "Sen 4hrs", "Sen D1", "Sen D2", "Sen D3", "Sen D4", "Sen D7")
sc5_top   <- sc5_conds[1:3]
sc5_bot   <- sc5_conds[4:7]
fit5 <- readRDS(file.path(v2_dir, "sc5", "fit_k10.rds"))

r5_full <- plot_consensus_cells("sc5", fit5, sc5_map, sc5_conds, seu_list, show_legend = FALSE)
n5_top  <- sum(r5_full$L_cells$condition %in% sc5_top)
n5_bot  <- sum(r5_full$L_cells$condition %in% sc5_bot)
x5_max  <- max(n5_top, n5_bot)

r5_top <- plot_consensus_cells("sc5", fit5, sc5_map, sc5_top, seu_list,
                                x_max = x5_max, show_legend = FALSE)
r5_bot <- plot_consensus_cells("sc5", fit5, sc5_map, sc5_bot, seu_list,
                                x_max = x5_max, show_legend = TRUE)

leg5    <- cowplot::get_legend(r5_bot$p)
p5_rows <- cowplot::plot_grid(r5_top$p, r5_bot$p + theme(legend.position = "none"),
                               nrow = 2, align = "v", axis = "lr")
p5_2row <- cowplot::plot_grid(p5_rows, leg5, ncol = 2, rel_widths = c(1, 0.12))

w5 <- max(180, x5_max * 0.18)
ggsave(file.path(v2_dir, "sc5", "structure_consensus_2row.png"),
       p5_2row, width = w5, height = 130, units = "mm", dpi = 300, bg = "white")
cat("  Saved sc5/structure_consensus_2row.png\n")

cat("\nDone.\n")
