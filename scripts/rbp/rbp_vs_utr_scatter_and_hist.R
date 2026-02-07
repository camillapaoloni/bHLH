#!/usr/bin/env Rscript

# RBP binding vs 3'UTR length (plus per-RBP histogram).
#
# Inputs:
# - data/intermediate/rbp/RBP_binary_matrix.csv
# - data/intermediate/rbp/Violin_plot_3UTR.csv
#
# Outputs:
# - outputs/rbp/figures/rbp_vs_utr_scatter.svg
# - outputs/rbp/figures/rbp_bound_tfs_per_rbp_histogram.svg

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
source(file.path(project_root, "scripts", "lib", "bhlh_utils.R"))

require_pkgs(c("dplyr", "tidyr", "ggplot2", "tibble"))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

in_rbp <- p("data", "intermediate", "rbp", "RBP_binary_matrix.csv")
in_utr <- p("data", "intermediate", "rbp", "Violin_plot_3UTR.csv")
if (!file.exists(in_rbp)) stop("Missing input: ", in_rbp)
if (!file.exists(in_utr)) stop("Missing input: ", in_utr)

out_dir <- p("outputs", "rbp", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_scatter <- file.path(out_dir, "rbp_vs_utr_scatter.svg")
out_hist <- file.path(out_dir, "rbp_bound_tfs_per_rbp_histogram.svg")

pal <- bhlh_class_palette()

rbp_mat <- read.csv(in_rbp, row.names = 1, check.names = FALSE)
if (ncol(rbp_mat) < 1) stop("RBP binary matrix looks empty: ", in_rbp)

rbp_mat[] <- lapply(rbp_mat, function(x) as.numeric(as.character(x)))
stopifnot(all(unlist(rbp_mat) %in% c(0, 1)))

utr <- read.csv(in_utr, stringsAsFactors = FALSE, check.names = FALSE)
names(utr)[names(utr) == ""] <- "row_id"
if ("Ledent2002+Simionato2007" %in% names(utr)) {
  utr <- dplyr::rename(utr, Ledent2002.Simionato2007 = `Ledent2002+Simionato2007`)
}
if ("HGNC symbol" %in% names(utr)) {
  utr <- dplyr::rename(utr, HGNC.symbol = `HGNC symbol`)
}
if (!all(c("HGNC.symbol", "UTR_length_mean", "Ledent2002.Simionato2007") %in% names(utr))) {
  stop("Violin_plot_3UTR.csv missing required columns (need HGNC.symbol, UTR_length_mean, Ledent2002.Simionato2007).")
}
utr <- utr %>%
  mutate(
    Ledent2002.Simionato2007 = ifelse(
      is.na(Ledent2002.Simionato2007) | Ledent2002.Simionato2007 == "",
      "NC",
      Ledent2002.Simionato2007
    ),
    Ledent2002.Simionato2007 = factor(Ledent2002.Simionato2007, levels = c("A", "B", "C", "D", "E", "NC"))
  ) %>%
  distinct(HGNC.symbol, .keep_all = TRUE)

# ---- Plot 1: scatter (per gene) ----
rbp_counts <- rowSums(rbp_mat)
df_scatter <- tibble(
  HGNC.symbol = names(rbp_counts),
  rbp_count = as.numeric(rbp_counts)
) %>%
  left_join(utr %>% select(HGNC.symbol, UTR_length_mean, Ledent2002.Simionato2007), by = "HGNC.symbol") %>%
  filter(!is.na(UTR_length_mean))

p_scatter <- ggplot(df_scatter, aes(x = UTR_length_mean, y = rbp_count, color = Ledent2002.Simionato2007)) +
  geom_point(size = 3.2, alpha = 0.6) +
  scale_color_manual(values = pal, drop = FALSE) +
  labs(
    title = "RBP binding vs 3'UTR length",
    subtitle = "Each point is a bHLH TF (color = class)",
    x = "Mean 3'UTR length",
    y = "Number of bound RBPs",
    color = "bHLH class"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    legend.position = "right"
  )

ggsave(out_scatter, plot = p_scatter, width = 9, height = 6, units = "in", dpi = 300)
message("Wrote: ", out_scatter)

# ---- Plot 2: histogram (per RBP) ----
binding_long <- rbp_mat %>%
  tibble::rownames_to_column("HGNC.symbol") %>%
  tidyr::pivot_longer(-HGNC.symbol, names_to = "RBP", values_to = "bound") %>%
  filter(bound == 1) %>%
  left_join(utr %>% select(HGNC.symbol, Ledent2002.Simionato2007), by = "HGNC.symbol")

plot_df <- binding_long %>%
  group_by(RBP, Ledent2002.Simionato2007) %>%
  summarise(n_tfs = n(), .groups = "drop")

p_hist <- ggplot(plot_df, aes(x = reorder(RBP, -n_tfs), y = n_tfs, fill = Ledent2002.Simionato2007)) +
  geom_col() +
  scale_fill_manual(values = pal, drop = FALSE) +
  labs(
    title = "bHLH TFs bound per RBP",
    subtitle = "Stacked bars show counts by bHLH class",
    x = "RBP",
    y = "Number of bound bHLH TFs",
    fill = "bHLH class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(out_hist, plot = p_hist, width = 14, height = 6, units = "in", dpi = 300)
message("Wrote: ", out_hist)
