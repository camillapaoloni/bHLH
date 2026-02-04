library(dplyr)
library(tidyr)
library(ggplot2)
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



#This plot will be used to create a scatter plot that shows the number of RBPs bound by each gene
#in respect to the length of the 3'UTR of that gene
RBP_matrix <- read.csv(p("data", "intermediate", "rbp", "RBP_binary_matrix.csv"), row.names = 1)
UTR_data <- read.csv(p("data", "intermediate", "rbp", "Violin_plot_3UTR.csv"))

output_dir <- p("outputs", "rbp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Normalize column names from Violin_plot_3UTR.csv
if ("Ledent2002+Simionato2007" %in% names(UTR_data)) {
  names(UTR_data)[names(UTR_data) == "Ledent2002+Simionato2007"] <- "Ledent2002.Simionato2007"
}
if ("HGNC symbol" %in% names(UTR_data)) {
  names(UTR_data)[names(UTR_data) == "HGNC symbol"] <- "HGNC.symbol"
}

# Count the number of RBPs bound per gene
RBP_counts <- rowSums(RBP_matrix)

# Create a data frame with the results
df <- data.frame(
  Gene = rownames(RBP_matrix),
  RBP_count = RBP_counts,
  UTR_length = UTR_data$UTR_length_mean,
  LS_class = UTR_data$Ledent2002.Simionato2007

)
my_colors <- c(
  "A" = "#7AC36A",   # green
  "B" = "#F7CB45",   # yellow
  "C" = "#9063CD",   # purple
  "D" = "#D83F3D",   # red
  "E" = "#4C86C6",   # blue
  "NC" = "grey"      # gray
)
# Create the scatter plot
p <- ggplot(df, aes(x = UTR_length, y = RBP_count, color = LS_class)) +
  geom_point(size = 4, alpha = 0.5) +
  labs(
    x = "Length of 3'UTR",
    y = "Numbers of RBPs bound",
    title = "Number of RBPs bound VS Length of 3'UTR",
    color = "Ledent&Simionato Classification"
  ) +
  scale_color_manual(values = my_colors) +
  theme_minimal()+ 
  theme(
    plot.title = element_text(hjust = 0.6, size = 18, face = "bold"),   # Center and enlarge title
    axis.title = element_text(size = 14, face = "italic"),              # Set axis title size/style
    legend.title = element_text(size = 13, face = "bold"),              # Set legend title style
    legend.text = element_text(size = 11)                               # Set legend text size
  )
#print(p)
ggsave(file.path(output_dir, "RBPvs3UTR_scatter.svg"), plot = p, width = 16, height = 12,dpi = 300, units = "in")

# PLOT 2: For each RBP, show the number of bound genes with a histogram
# whose height is colored based on class.
# Required columns: class, TF names, and bound RBP per TF.


# 1. Convert binary matrix to long format
binding_long <- as.data.frame(RBP_matrix) %>%
  tibble::rownames_to_column("TF") %>%
  pivot_longer(-TF, names_to = "RBP", values_to = "bound") %>%
  filter(bound == 1)

# 2. Add class information
binding_long <- binding_long %>%
  left_join(UTR_data[, c("HGNC.symbol", "Ledent2002.Simionato2007")], by = c("TF" = "HGNC.symbol"))

# 3. Count number of TFs per RBP per class
plot_df <- binding_long %>%
  group_by(RBP, Ledent2002.Simionato2007) %>%
  summarise(n_genes = n(), .groups = "drop")

# 4. Plot
p1 <- ggplot(plot_df, aes(x = reorder(RBP, -n_genes), y = n_genes, fill = Ledent2002.Simionato2007)) +
  geom_col() +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  labs(
    x = "RBPs",
    y = "Number of TFs bound",
    fill = "TF Class",
    title = "TFs Bound per RBP Colored by Class"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.6, size = 18, face = "bold"),   # Center and enlarge title
    axis.title = element_text(size = 14, face = "italic"),              # Set axis title size/style
    legend.title = element_text(size = 14, face = "bold"),              # Set legend title style
    legend.text = element_text(size = 12)                               # Set legend text size
  )
print(p1)
ggsave(file.path(output_dir, "RBPvsTFs_histogram.svg"), plot = p1, width = 30, height = 12,dpi = 300, units = "in")
