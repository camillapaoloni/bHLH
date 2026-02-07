#!/usr/bin/env Rscript

# t-SNE of RBP binding profiles with K-modes clustering.
#
# Input:
# - data/intermediate/rbp/RBP_binary_matrix.csv
# - data/raw/LS_classes.csv
#
# Output:
# - outputs/rbp/figures/rbp_tsne_kmodes_clusters.svg

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
source(file.path(project_root, "scripts", "lib", "bhlh_utils.R"))

require_pkgs(c("klaR", "Rtsne", "ggplot2", "dplyr", "ggrepel", "tibble", "viridis"))
suppressPackageStartupMessages({
  library(klaR)
  library(Rtsne)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(tibble)
  library(viridis)
})

k <- suppressWarnings(as.integer(Sys.getenv("BHLH_RBP_K", unset = "2")))
if (is.na(k) || k < 2) stop("Invalid BHLH_RBP_K (must be >= 2).")

in_rbp <- p("data", "intermediate", "rbp", "RBP_binary_matrix.csv")
if (!file.exists(in_rbp)) stop("Missing input: ", in_rbp)

out_dir <- p("outputs", "rbp", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_svg <- file.path(out_dir, "rbp_tsne_kmodes_clusters.svg")

rbp_mat <- read.csv(in_rbp, row.names = 1, check.names = FALSE)
rbp_mat[] <- lapply(rbp_mat, function(x) as.numeric(as.character(x)))
stopifnot(all(unlist(rbp_mat) %in% c(0, 1)))

classes <- read_ls_classes()
class_map <- classes %>% dplyr::select(HGNC.symbol, Ledent2002.Simionato2007) %>% dplyr::distinct()
pal <- bhlh_class_palette()

set.seed(42)
kmodes_res <- kmodes(rbp_mat, modes = k, iter.max = 20)
rbp_mat$Cluster <- as.factor(kmodes_res$cluster)

rbp_unique <- rbp_mat[!duplicated(rbp_mat[, -ncol(rbp_mat)]), ]
set.seed(42)
tsne_res <- Rtsne(as.matrix(rbp_unique[, -ncol(rbp_unique)]), perplexity = 10, check_duplicates = FALSE)

tsne_df <- tibble(
  tSNE1 = tsne_res$Y[, 1],
  tSNE2 = tsne_res$Y[, 2],
  Cluster = rbp_unique$Cluster,
  HGNC.symbol = rownames(rbp_unique)
) %>%
  left_join(class_map, by = "HGNC.symbol") %>%
  mutate(
    Ledent2002.Simionato2007 = ifelse(
      is.na(Ledent2002.Simionato2007) | Ledent2002.Simionato2007 == "",
      "NC",
      Ledent2002.Simionato2007
    ),
    Ledent2002.Simionato2007 = factor(Ledent2002.Simionato2007, levels = c("A", "B", "C", "D", "E", "NC"))
  )

centroids <- tsne_df %>% group_by(Cluster) %>% summarise(tSNE1 = mean(tSNE1), tSNE2 = mean(tSNE2), .groups = "drop")
hull_data <- tsne_df %>% group_by(Cluster) %>% slice(chull(tSNE1, tSNE2))

p_tsne <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2)) +
  geom_polygon(data = hull_data, aes(fill = Cluster), alpha = 0.18, show.legend = FALSE) +
  geom_point(aes(color = Ledent2002.Simionato2007, shape = Cluster), size = 2.5, alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = centroids,
    aes(x = tSNE1, y = tSNE2, label = paste0("Cluster ", Cluster)),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey60"
  ) +
  scale_color_manual(values = pal, drop = FALSE) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(
    title = "t-SNE of RBP binding profiles (K-modes clustering)",
    subtitle = paste0("K-modes (k = ", k, ") | point color = bHLH class"),
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "bHLH class",
    shape = "Cluster"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    legend.position = "right"
  )

ggsave(out_svg, plot = p_tsne, width = 9, height = 6, units = "in", dpi = 300)
message("Wrote: ", out_svg)
