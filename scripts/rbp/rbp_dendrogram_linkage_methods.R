#!/usr/bin/env Rscript

# Dendrograms of RBP binding profiles across bHLH TFs.
#
# Input:
# - data/intermediate/rbp/RBP_binary_matrix.csv
# - data/raw/LS_classes.csv
#
# Output:
# - outputs/rbp/figures/rbp_dendrograms_linkage_methods.svg
#
# Notes:
# - We color leaf labels by bHLH class (A–E, NC) to keep the dendrogram readable.
# - Branch coloring by class is ambiguous when clusters mix classes, so we avoid it here.

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
source(file.path(project_root, "scripts", "lib", "bhlh_utils.R"))

require_pkgs(c("ggplot2", "dendextend", "patchwork", "dplyr"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dendextend)
  library(patchwork)
  library(dplyr)
})

exclude_genes <- Sys.getenv("BHLH_RBP_EXCLUDE_GENES", unset = "")
exclude_genes <- if (nzchar(exclude_genes)) strsplit(exclude_genes, ",")[[1]] else character(0)
exclude_genes <- trimws(exclude_genes)

in_rbp <- p("data", "intermediate", "rbp", "RBP_binary_matrix.csv")
if (!file.exists(in_rbp)) stop("Missing input: ", in_rbp)

out_dir <- p("outputs", "rbp", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_svg <- file.path(out_dir, "rbp_dendrograms_linkage_methods.svg")

rbp_mat <- read.csv(in_rbp, row.names = 1, check.names = FALSE)
if (length(exclude_genes) > 0) {
  rbp_mat <- rbp_mat[!rownames(rbp_mat) %in% exclude_genes, , drop = FALSE]
}

classes <- read_ls_classes()
class_map <- classes %>% select(HGNC.symbol, Ledent2002.Simionato2007) %>% distinct()
pal <- bhlh_class_palette()

leaf_classes <- class_map$Ledent2002.Simionato2007[match(rownames(rbp_mat), class_map$HGNC.symbol)]
leaf_classes <- ifelse(is.na(leaf_classes) | leaf_classes == "", "NC", leaf_classes)
leaf_cols <- pal[leaf_classes]
leaf_cols[is.na(leaf_cols)] <- pal[["NC"]]

# Jaccard distance on binary profiles.
dist_jaccard <- dist(rbp_mat, method = "binary")

plot_dendro <- function(linkage_method, title) {
  hc <- hclust(dist_jaccard, method = linkage_method)
  dend <- as.dendrogram(hc)
  dend <- dendextend::set(dend, "labels_col", leaf_cols)
  dend <- dendextend::set(dend, "labels_cex", 0.55)

  ggd <- dendextend::as.ggdend(dend)
  ggplot(ggd) +
    labs(title = title, x = NULL, y = "Linkage distance") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
}

p_single <- plot_dendro("single", "Linkage: single")
p_complete <- plot_dendro("complete", "Linkage: complete")
p_average <- plot_dendro("average", "Linkage: average")
p_ward <- plot_dendro("ward.D2", "Linkage: ward.D2")

final_plot <- p_single + p_complete + p_average + p_ward + plot_layout(ncol = 2)
ggsave(out_svg, plot = final_plot, width = 14, height = 10, units = "in", dpi = 300)
message("Wrote: ", out_svg)
