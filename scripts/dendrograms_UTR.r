
library(ggplot2)
library(dendextend)
library(patchwork)
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



# Load data
pdata <- read.csv(p("data", "intermediate", "rbp", "RBP_binary_matrix.csv"), row.names = 1)
classes <- read.csv(p("data", "raw", "LS_classes.csv"))
classes <- normalize_classes(classes)
# Remove BMAL1 and BMAL2
pdata <- pdata[!rownames(pdata) %in% c("BMAL1", "BMAL2"), ]

# Assign class to each gene
gene_classes <- setNames(classes$Ledent2002.Simionato2007, classes$HGNC.symbol)
gene_class_vector <- gene_classes[rownames(pdata)]
# Define color palette
my_colors <- c(
  "A" = "#7AC36A",
  "B" = "#F7CB45",
  "C" = "#9063CD",
  "D" = "#D83F3D",
  "E" = "#4C86C6",
  "NC" = "grey"
)
# Function to color branches by class
color_branches_by_class <- function(dend, classes, palette) {
  # Assign a color to each class
  col_vector <- palette[as.character(classes)]
  # Color branches based on leaf class
  dend <- dendextend::color_branches(dend, col = col_vector, groupLabels = FALSE)
  # Color labels as well
  dend <- dendextend::set(dend, "labels_col", col_vector)
  return(dend)
}
# Compute distance matrix (Jaccard)
dist_jaccard <- dist(pdata, method = "binary")

# Function to create a colored dendrogram
plot_dendro <- function(linkage_method, title) {
  hc <- hclust(dist_jaccard, method = linkage_method)
  dend <- as.dendrogram(hc)
  dend <- color_branches_by_class(dend, gene_class_vector, my_colors)
  # Create the plot
  ggd1 <- as.ggdend(dend)
  ggd1$labels$label.size <- 10
  p <- ggplot(ggd1) + 
    labs(title = title, x = "Ledent&Simionato Classification", y= "Linkage Distance") +
    theme_minimal(base_size = 14, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 13, face = "italic"),
      axis.title.y = element_text(size = 13, face = "italic"),
    ) + scale_y_continuous(expand = expansion(mult = c(.2, NA)))
  return(p)
}

# Create dendrograms with different linkage methods
p_single   <- plot_dendro("single",   "Linkage method: single")
p_complete <- plot_dendro("complete", "Linkage method: complete")
p_average  <- plot_dendro("average",  "Linkage method: average")
p_ward     <- plot_dendro("ward.D2",  "Linkage method: ward.D2")

p_single[[2]][[3]]$mapping$size <- 3.5
p_complete[[2]][[3]]$mapping$size <- 3.5
p_average[[2]][[3]]$mapping$size <- 3.5
p_ward[[2]][[3]]$mapping$size <- 3.5



# Show dendrograms side by side
final_plot <- p_single + p_complete + p_average + p_ward + plot_layout(ncol = 2)
print(final_plot)
output_dir <- p("outputs", "rbp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "dendrograms_linkage_methods.svg"), plot = final_plot, width = 20, height = 12)
