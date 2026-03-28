# sc5_timecourse_consensus.R
# Topic proportion trajectories across the H2O2 senescence timecourse (sc5).
# Two plots:
#   1. Overlaid line plot — all topics on one panel
#   2. Faceted plot — one panel per topic, easier to read individual dynamics

library(ggplot2)
library(dplyr)
library(tidyr)
library(fastTopics)
library(Seurat)

v2_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/lda_poc_v2"
qc_dir <- "D:/Users/adb/Topic Modeling/senescence-topic-modeling/analysis/qc_output"

sc5_map <- c("Topic 1"="A","Topic 2"="B","Topic 3"="J","Topic 4"="F",
             "Topic 5"="G","Topic 6"="D","Topic 7"="K","Topic 8"="L",
             "Topic 9"="C","Topic 10"="E")
letter_order <- unname(sc5_map)

consensus_cols <- c(
  H="#D55E00", A="#0072B2", B="#009E73", F="#E69F00",
  D="#56B4E9", E="#CC79A7", C="#44AA99", G="#882255",
  I="#DDCC77", K="#AA4499", J="#332288", L="#999933"
)
topic_cols <- consensus_cols[letter_order]

tp_order <- c("Control","Sen 4hrs","Sen D1","Sen D2","Sen D3","Sen D4","Sen D7")
tp_days  <- c(Control=0, "Sen 4hrs"=0.17, "Sen D1"=1,
              "Sen D2"=2, "Sen D3"=3, "Sen D4"=4, "Sen D7"=7)

# ---- Load and normalise -----------------------------------------------------
cat("Loading data...\n")
seu  <- readRDS(file.path(qc_dir, "seu_list_qc.rds"))[["sc5"]]
fit5 <- readRDS(file.path(v2_dir, "sc5", "fit_k10.rds"))

L_norm <- fit5$L / rowSums(fit5$L)
colnames(L_norm) <- letter_order

L_df <- as.data.frame(L_norm)
L_df$cell_id   <- rownames(fit5$L)
L_df$timepoint <- factor(seu@meta.data[L_df$cell_id, "condition"], levels = tp_order)
L_df <- L_df[!is.na(L_df$timepoint), ]
L_df$day <- tp_days[as.character(L_df$timepoint)]

# ---- Summarise: mean + SEM per topic per timepoint -------------------------
L_long <- L_df |>
  select(day, timepoint, all_of(letter_order)) |>
  pivot_longer(all_of(letter_order), names_to = "topic", values_to = "proportion") |>
  mutate(topic = factor(topic, levels = letter_order))

summary_df <- L_long |>
  group_by(topic, day, timepoint) |>
  summarise(
    mean_prop = mean(proportion),
    sem_prop  = sd(proportion) / sqrt(n()),
    .groups   = "drop"
  )

day_breaks <- c(0, 0.17, 1, 2, 3, 4, 7)
day_labels <- c("Ctrl", "4h", "D1", "D2", "D3", "D4", "D7")

# ---- Plot 1: Overlaid -------------------------------------------------------
p_overlay <- ggplot(summary_df, aes(x = day, y = mean_prop, color = topic, fill = topic)) +
  geom_ribbon(aes(ymin = mean_prop - sem_prop, ymax = mean_prop + sem_prop),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(values = topic_cols) +
  scale_fill_manual(values  = topic_cols) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 9) +
  theme(legend.key.size = unit(4, "mm")) +
  labs(x = "Day post H\u2082O\u2082", y = "Mean topic proportion",
       color = NULL, fill = NULL,
       title = "sc5 topic trajectories across senescence timecourse")

ggsave(file.path(v2_dir, "sc5", "timecourse_consensus_overlay.png"),
       p_overlay, width = 150, height = 90, units = "mm", dpi = 300, bg = "white")
cat("Saved timecourse_consensus_overlay.png\n")

# ---- Plot 2: Faceted --------------------------------------------------------
# Add Hallmark annotation as subtitle per facet
hallmark_ann <- c(
  A = "UV Resp DN / p53-damage",
  B = "EMT / Activated fibroblast",
  J = "EMT (strong)",
  F = "EMT / TGF-\u03b2",
  G = "EMT",
  D = "EMT / Hedgehog",
  K = "Coagulation / SASP",
  L = "Hypoxia / Glycolysis",
  C = "Cycling / G2M",
  E = "EMT (late)"
)

summary_df$label <- paste0(summary_df$topic, ": ", hallmark_ann[as.character(summary_df$topic)])
# keep factor order
lbl_levels <- paste0(letter_order, ": ", hallmark_ann[letter_order])
summary_df$label <- factor(summary_df$label, levels = lbl_levels)

p_facet <- ggplot(summary_df, aes(x = day, y = mean_prop, color = topic, fill = topic)) +
  geom_ribbon(aes(ymin = mean_prop - sem_prop, ymax = mean_prop + sem_prop),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(values = topic_cols) +
  scale_fill_manual(values  = topic_cols) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ label, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold"),
        axis.text.x = element_text(size = 6),
        panel.spacing = unit(3, "mm")) +
  labs(x = "Day post H\u2082O\u2082", y = "Mean proportion",
       title = "sc5 — topic trajectories (mean \u00b1 SEM)")

ggsave(file.path(v2_dir, "sc5", "timecourse_consensus_facet.png"),
       p_facet, width = 240, height = 110, units = "mm", dpi = 300, bg = "white")
cat("Saved timecourse_consensus_facet.png\n")

cat("Done.\n")
