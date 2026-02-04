# === Load Libraries ===
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)


# === Load Data ===
binary_matrix <- read.csv(p("data", "intermediate", "rbp", "RBP_binary_matrix.csv"), row.names = 1)
binary_matrix[] <- lapply(binary_matrix, function(x) as.numeric(as.character(x)))
stopifnot(all(unlist(binary_matrix) %in% c(0, 1)))

# === Hierarchical Clustering (Ward.D2) ===
diss <- daisy(binary_matrix, metric = "gower")
hc <- hclust(diss, method = "ward.D2")
k <- 2
cluster_assignments <- cutree(hc, k = k)
binary_matrix$Cluster <- as.factor(cluster_assignments)

# === t-SNE ===
binary_matrix_unique <- binary_matrix[!duplicated(binary_matrix[, -ncol(binary_matrix)]), ]
set.seed(42)
tsne_res <- Rtsne(as.matrix(binary_matrix_unique[, -ncol(binary_matrix_unique)]), perplexity = 10, check_duplicates = FALSE)

tsne_df <- data.frame(
  tSNE1 = tsne_res$Y[, 1],
  tSNE2 = tsne_res$Y[, 2],
  Cluster = binary_matrix_unique$Cluster
)

# === Visualization ===
centroids <- tsne_df %>%
  group_by(Cluster) %>%
  summarise(tSNE1 = mean(tSNE1), tSNE2 = mean(tSNE2))

hull_data <- tsne_df %>%
  group_by(Cluster) %>%
  slice(chull(tSNE1, tSNE2))

p_hierarchical <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
  geom_polygon(data = hull_data, aes(fill = Cluster), alpha = 0.2, show.legend = FALSE) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text_repel(data = centroids, aes(label = paste("Cluster", Cluster)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(
    title = "t-SNE Visualization of Gene Clusters (Hierarchical)",
    x = "t-SNE 1", y = "t-SNE 2"
  ) +
  theme_minimal(base_size = 14)

print(p_hierarchical)
output_dir <- p("outputs", "rbp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "hierarchical_clusters.svg"), p_hierarchical, width = 8, height = 6)