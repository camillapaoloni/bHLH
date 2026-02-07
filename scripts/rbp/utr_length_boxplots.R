#!/usr/bin/env Rscript

# 3'UTR length summaries by bHLH class.
#
# Input:
# - data/intermediate/rbp/Violin_plot_3UTR.csv (produced by notebooks/3'UTR.ipynb)
#
# Output:
# - outputs/rbp/figures/utr_length_boxplots.svg

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
source(file.path(project_root, "scripts", "lib", "bhlh_utils.R"))

require_pkgs(c("ggplot2", "patchwork", "dplyr"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

in_path <- p("data", "intermediate", "rbp", "Violin_plot_3UTR.csv")
if (!file.exists(in_path)) stop("Missing input: ", in_path)

out_dir <- p("outputs", "rbp", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_svg <- file.path(out_dir, "utr_length_boxplots.svg")

df <- read.csv(in_path, stringsAsFactors = FALSE, check.names = FALSE)
names(df)[names(df) == ""] <- "row_id"
if ("Ledent2002+Simionato2007" %in% names(df)) {
  df <- dplyr::rename(df, Ledent2002.Simionato2007 = `Ledent2002+Simionato2007`)
}
df$Ledent2002.Simionato2007 <- ifelse(
  is.na(df$Ledent2002.Simionato2007) | df$Ledent2002.Simionato2007 == "",
  "NC",
  df$Ledent2002.Simionato2007
)
df$Ledent2002.Simionato2007 <- factor(df$Ledent2002.Simionato2007, levels = c("A", "B", "C", "D", "E", "NC"))

pal <- bhlh_class_palette()

make_box <- function(ycol, title) {
  ggplot(df, aes(x = Ledent2002.Simionato2007, y = .data[[ycol]], fill = Ledent2002.Simionato2007)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 1.4) +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(title = title, x = "bHLH class", y = "3'UTR length") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title.x = element_text(face = "bold")
    )
}

p_median <- make_box("UTR_length_median", "Median 3'UTR length")
p_mean <- make_box("UTR_length_mean", "Mean 3'UTR length")
p_max <- make_box("UTR_length_max", "Max 3'UTR length")

final_plot <- p_median + p_mean + p_max + plot_layout(ncol = 3)
ggsave(out_svg, plot = final_plot, width = 14, height = 5, units = "in", dpi = 300)
message("Wrote: ", out_svg)
