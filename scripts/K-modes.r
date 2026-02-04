# === Load Libraries ===
library(klaR)       # For K-modes
library(Rtsne)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)
normalize_classes <- function(df) {
  if ("Ledent2002+Simionato2007" %in% names(df)) {
    df <- dplyr::rename(df, Ledent2002.Simionato2007 = `Ledent2002+Simionato2007`)
  }
  if (!"Ledent2002.Simionato2007" %in% names(df)) {
    stop("Class column not found in LS_classes")
  }
  df$Ledent2002.Simionato2007 <- ifelse(is.na(df$Ledent2002.Simionato2007) | df$Ledent2002.Simionato2007 == "",
                                          "NC", df$Ledent2002.Simionato2007)
  df
}



# === Load Data ===
binary_matrix <- read.csv(p("data", "intermediate", "rbp", "RBP_binary_matrix.csv"), row.names = 1)
LS_classes <- read.csv(p("data", "raw", "LS_classes.csv"))
LS_classes <- normalize_classes(LS_classes)

my_colors <- c(
"A" = "#7AC36A",
"B" = "#F7CB45",
"C" = "#9063CD",
"D" = "#D83F3D",
"E" = "#4C86C6",
"NC" = "grey"
)


binary_matrix[] <- lapply(binary_matrix, function(x) as.numeric(as.character(x)))
stopifnot(all(unlist(binary_matrix) %in% c(0, 1)))

# === K-Modes Clustering ===
k <- 2
set.seed(42)
kmodes_res <- kmodes(binary_matrix, modes = k, iter.max = 10)
binary_matrix$Cluster <- as.factor(kmodes_res$cluster)

# === t-SNE ===
binary_matrix_unique <- binary_matrix[!duplicated(binary_matrix[, -ncol(binary_matrix)]), ]
set.seed(42)
tsne_res <- Rtsne(as.matrix(binary_matrix_unique[, -ncol(binary_matrix_unique)]), perplexity = 10, check_duplicates = FALSE)

tsne_df <- data.frame(
  tSNE1 = tsne_res$Y[, 1],
  tSNE2 = tsne_res$Y[, 2],
  Cluster = binary_matrix_unique$Cluster
)
# Add a column with gene names (or key) to tsne_df
tsne_df$Gene <- rownames(binary_matrix_unique)

# Ensure LS_classes has a matching 'Gene' column
# Adjust the key column name if needed
LS_classes$Gene <- as.character(LS_classes$HGNC.symbol) # or the column name with gene names

# Join the data
tsne_df <- left_join(tsne_df, LS_classes[, c("Gene", "Ledent2002.Simionato2007")], by = "Gene")


# === Visualization ===
centroids <- tsne_df %>%
  group_by(Cluster) %>%
  summarise(tSNE1 = mean(tSNE1), tSNE2 = mean(tSNE2))

hull_data <- tsne_df %>%
  group_by(Cluster) %>%
  slice(chull(tSNE1, tSNE2))


# ...existing code...

p_kmodes <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2)) +
  geom_polygon(data = hull_data, aes(fill = Cluster), alpha = 0.2, show.legend = FALSE) +
  geom_point(aes(color = Ledent2002.Simionato2007), size = 2.5, alpha = 0.8) +
  geom_text_repel(data = centroids, aes(label = paste("Cluster", Cluster), x = tSNE1, y = tSNE2), size = 4, inherit.aes = FALSE) +
  scale_color_manual(values = my_colors, na.value = "black") +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(
    title = "t-SNE Visualization of Gene Clusters (K-Modes)",
    x = "t-SNE 1", y = "t-SNE 2",
    color = "Ledent2002.Simionato2007"
  ) +
  theme_minimal(base_size = 14)+

  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(
    title = "t-SNE Visualization of Gene Clusters (K-Modes)",
    x = "t-SNE 1", y = "t-SNE 2",
    color = "Ledent2002.Simionato2007"
  ) +
  theme_minimal(base_size = 14)

print(p_kmodes)
output_dir <- p("outputs", "rbp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "k_modes_clusters.svg"), p_kmodes, width = 8, height = 6)
