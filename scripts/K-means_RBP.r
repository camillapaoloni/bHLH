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

# Convert all entries to numeric 0/1
binary_matrix[] <- lapply(binary_matrix, function(x) as.numeric(as.character(x)))

# Sanity check: make sure all values are 0 or 1
stopifnot(all(unlist(binary_matrix) %in% c(0, 1)))

# === PAM Clustering using Gower Distance ===
diss <- daisy(binary_matrix, metric = "gower")
k <- 2
pam_res <- pam(diss, k = k, diss = TRUE)
binary_matrix$Cluster <- as.factor(pam_res$clustering)

# === t-SNE for Visualization ===
# Remove duplicated rows for Rtsne
binary_matrix_unique <- binary_matrix[!duplicated(binary_matrix[, -ncol(binary_matrix)]), ]

set.seed(42)
tsne_res <- Rtsne(as.matrix(binary_matrix_unique[, -ncol(binary_matrix_unique)]), perplexity = 10, check_duplicates = FALSE)

# Create data frame for ggplot
tsne_df <- data.frame(
  tSNE1 = tsne_res$Y[, 1],
  tSNE2 = tsne_res$Y[, 2],
  Cluster = binary_matrix_unique$Cluster
)

# === Advanced Visualization with ggplot2 ===

# 1. Calculate cluster centroids for labeling
tsne_centroids <- tsne_df %>%
  group_by(Cluster) %>%
  summarise(tSNE1 = mean(tSNE1), tSNE2 = mean(tSNE2))

# 2. Create convex hulls for each cluster
hull_data <- tsne_df %>%
  group_by(Cluster) %>%
  slice(chull(tSNE1, tSNE2))

# 3. Create the improved plot
p_improved <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
  # Draw convex hulls
  geom_polygon(data = hull_data, aes(fill = Cluster), alpha = 0.2, show.legend = FALSE) +
  # Add points
  geom_point(size = 2.5, alpha = 0.8) +
  # Add cluster labels using ggrepel to avoid overlap
  ggrepel::geom_text_repel(
    data = tsne_centroids,
    aes(label = paste("Cluster", Cluster)),
    fontface = "bold",
    size = 4,
    color = "black", # Use a single color for labels to avoid clutter
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  ) +
  # Use a better color palette (viridis is colorblind-friendly)
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(
    title = "t-SNE Visualization of Gene Clusters",
    subtitle = "Clusters identified using PAM with Gower Distance",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2",
    color = "Gene Cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold")
  )

print(p_improved)
output_dir <- p("outputs", "rbp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "k_means_clusters.svg"), p_improved, width = 8, height = 6)