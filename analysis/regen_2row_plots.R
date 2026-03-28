# regen_2row_plots.R
# Regenerate 2-row structure plots for sc78 and sc5 with correct cell ordering.
# Uses original topic labels/colors (not consensus), same sort logic as working sc2 plot.

library(ggplot2)
library(dplyr)
library(tidyr)
library(fastTopics)
library(Seurat)
library(cowplot)

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"
qc_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"

oi_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00",
            "#CC79A7","#000000","#88CCEE","#44AA99","#332288","#DDCC77")

cat("Loading Seurat objects...\n")
seu_list <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))

# ---- Helper: one structure row with correct sorting -------------------------
make_row <- function(L_cells, topic_lbls, topic_cols, conds_subset,
                     show_legend = FALSE, x_max = NULL) {

  row_cells <- L_cells[L_cells$condition %in% conds_subset, ]
  row_cells$condition <- factor(row_cells$condition, levels = conds_subset)

  # Sort: by dominant topic (in topic_lbls order = Topic 1, 2, 3...),
  # then descending weight of Topic 1 — identical to working sc2 logic
  row_cells$dominant      <- topic_lbls[max.col(as.matrix(row_cells[, topic_lbls, drop = FALSE]))]
  row_cells$dominant_rank <- match(row_cells$dominant, topic_lbls)
  row_cells <- row_cells[order(row_cells$condition, row_cells$dominant_rank,
                                -row_cells[[topic_lbls[1]]]), ]
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
    theme(axis.text.x   = element_blank(),
          axis.ticks.x  = element_blank(),
          axis.line.x   = element_blank(),
          legend.key.size = unit(3.5, "mm"),
          plot.margin   = margin(14, 4, 4, 4, "mm")) +
    labs(x = NULL, y = "Topic proportion", fill = NULL)

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

# ============================================================
# sc78 K=9
# ============================================================
cat("\nsc78...\n")
K78         <- 9L
topic_lbls78 <- paste0("Topic ", seq_len(K78))
topic_cols78 <- oi_pal[seq_len(K78)]
sc78_conds  <- c("sc7 Control","sc7 Nutlin","sc7 OxStress","sc7 Quiescence",
                 "sc8 Control","sc8 Nutlin","sc8 OxStress","sc8 Quiescence")
sc7_conds   <- sc78_conds[1:4]
sc8_conds   <- sc78_conds[5:8]

fit78    <- readRDS(file.path(v2_dir, "sc78", "fit_k9.rds"))
L78_norm <- fit78$L / rowSums(fit78$L)
colnames(L78_norm) <- topic_lbls78
L78_df   <- as.data.frame(L78_norm)
L78_df$cell_id   <- rownames(fit78$L)
L78_df$condition <- factor(seu_list[["sc78"]]@meta.data[L78_df$cell_id, "condition"],
                            levels = sc78_conds)
L78_df <- L78_df[!is.na(L78_df$condition), ]

set.seed(42)
L78_cells <- L78_df |>
  group_by(condition) |>
  group_modify(~ slice_sample(.x, n = min(700L, nrow(.x)))) |>
  ungroup()

n78_top <- sum(L78_cells$condition %in% sc7_conds)
n78_bot <- sum(L78_cells$condition %in% sc8_conds)
x78_max <- max(n78_top, n78_bot)

p78_top <- make_row(L78_cells, topic_lbls78, topic_cols78, sc7_conds,
                    show_legend = FALSE, x_max = x78_max)
p78_bot <- make_row(L78_cells, topic_lbls78, topic_cols78, sc8_conds,
                    show_legend = TRUE, x_max = x78_max)

leg78    <- cowplot::get_legend(p78_bot)
p78_rows <- cowplot::plot_grid(p78_top, p78_bot + theme(legend.position = "none"),
                                nrow = 2, align = "v", axis = "lr")
p78_2row <- cowplot::plot_grid(p78_rows, leg78, ncol = 2, rel_widths = c(1, 0.12))

w78 <- max(180, x78_max * 0.18)
ggsave(file.path(v2_dir, "sc78", "structure_plot_cells_k9_2row.png"),
       p78_2row, width = w78, height = 130, units = "mm", dpi = 300, bg = "white")
cat("  Saved sc78/structure_plot_cells_k9_2row.png\n")

# ============================================================
# sc5 K=10
# ============================================================
cat("\nsc5...\n")
K5          <- 10L
topic_lbls5 <- paste0("Topic ", seq_len(K5))
topic_cols5 <- oi_pal[seq_len(K5)]
tp_order    <- c("Control","Sen 4hrs","Sen D1","Sen D2","Sen D3","Sen D4","Sen D7")
sc5_top     <- tp_order[1:3]
sc5_bot     <- tp_order[4:7]

fit5     <- readRDS(file.path(v2_dir, "sc5", "fit_k10.rds"))
L5_norm  <- fit5$L / rowSums(fit5$L)
colnames(L5_norm) <- topic_lbls5
L5_df    <- as.data.frame(L5_norm)
L5_df$cell_id   <- rownames(fit5$L)
L5_df$condition <- factor(seu_list[["sc5"]]@meta.data[L5_df$cell_id, "condition"],
                           levels = tp_order)
L5_df <- L5_df[!is.na(L5_df$condition), ]

set.seed(42)
L5_cells <- L5_df |>
  group_by(condition) |>
  group_modify(~ slice_sample(.x, n = min(700L, nrow(.x)))) |>
  ungroup()

n5_top <- sum(L5_cells$condition %in% sc5_top)
n5_bot <- sum(L5_cells$condition %in% sc5_bot)
x5_max <- max(n5_top, n5_bot)

p5_top <- make_row(L5_cells, topic_lbls5, topic_cols5, sc5_top,
                   show_legend = FALSE, x_max = x5_max)
p5_bot <- make_row(L5_cells, topic_lbls5, topic_cols5, sc5_bot,
                   show_legend = TRUE, x_max = x5_max)

leg5    <- cowplot::get_legend(p5_bot)
p5_rows <- cowplot::plot_grid(p5_top, p5_bot + theme(legend.position = "none"),
                               nrow = 2, align = "v", axis = "lr")
p5_2row <- cowplot::plot_grid(p5_rows, leg5, ncol = 2, rel_widths = c(1, 0.12))

w5 <- max(180, x5_max * 0.18)
ggsave(file.path(v2_dir, "sc5", "structure_plot_cells_k10_2row.png"),
       p5_2row, width = w5, height = 130, units = "mm", dpi = 300, bg = "white")
cat("  Saved sc5/structure_plot_cells_k10_2row.png\n")

cat("\nDone.\n")
