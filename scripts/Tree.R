# Load necessary libraries
library(ggtree)       # for phylogenetic tree visualization
library(ape)          # phylogenetic tree tools
library(dplyr)        # data manipulation (optional here)
library(RColorBrewer) # color palettes
library(svglite)      # save plots as SVG
library(ggplot2)      # additional plotting utilities
library(patchwork)    # combine multiple plots
project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)


# 1. Load the phylogenetic tree from Newick file
tree_path <- if (file.exists(p("data", "raw", "Tree23sp.newick"))) p("data", "raw", "Tree23sp.newick") else p("data", "raw", "Tree23sp.nwk")
tree <- read.tree(tree_path)

# 2. Define desired tip order (underscore notation)
species_order <- c(
  "Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae", "Neurospora_crassa",
  "Caenorhabditis_elegans", "Anopheles_gambiae", "Drosophila_melanogaster",
  "Tribolium_castaneum", "Helobdella_robusta",
  "Danio_rerio", "Oryzias_latipes", "Lepisosteus_oculatus",
  "Xenopus_tropicalis", "Gallus_gallus", "Anolis_carolinensis",
  "Monodelphis_domestica", "Bos_taurus", "Canis_lupus_familiaris",
  "Mus_musculus", "Rattus_norvegicus",
  "Macaca_mulatta", "Gorilla_gorilla", "Pan_troglodytes", "Homo_sapiens"
)

# 3. Reverse order so Homo sapiens appears first (top of tree)
species_order_rev <- rev(species_order)

# 4. Reorder tree tips to match desired order
tree <- rotateConstr(tree, species_order_rev)

# 5. Initialize base tree plot with rectangular layout
p <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(size = 3, fontface = "italic", align = TRUE, linesize = 0.5) +  # italic species names, aligned tips
  theme_tree2() +  # cleaner style with axis lines for better readability
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank()
  ) +
  ggtitle("Phylogenetic Tree of Selected Species")

# 6. Define clades and their member species
clades <- list(
  "Primates" = c("Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla", "Macaca_mulatta"),
  "Rodents" = c("Mus_musculus", "Rattus_norvegicus"),
  "Ungulates and Carnivores" = c("Bos_taurus", "Canis_lupus_familiaris"),
  "Marsupials" = "Monodelphis_domestica",
  "Birds and Reptiles" = c("Gallus_gallus", "Anolis_carolinensis"),
  "Amphibians and Fish" = c("Danio_rerio", "Oryzias_latipes", "Xenopus_tropicalis", "Lepisosteus_oculatus"),
  "Arthropods" = c("Drosophila_melanogaster", "Anopheles_gambiae", "Tribolium_castaneum"),
  "Worms" = c("Caenorhabditis_elegans", "Helobdella_robusta"),
  "Fungi" = c("Schizosaccharomyces_pombe", "Saccharomyces_cerevisiae", "Neurospora_crassa")
)

# 7. Select a soft, colorblind-friendly palette for clades
clade_colors <- brewer.pal(n = length(clades), name = "Set3")

# 8. Add highlights and clade labels
i <- 1
for (clade_name in names(clades)) {
  species <- clades[[clade_name]]
  species_in_tree <- species[species %in% tree$tip.label]
  
  if (length(species_in_tree) >= 2) {
  node <- getMRCA(tree, species_in_tree)
  if (!is.na(node)) {
    p <- p +
      geom_hilight(node = node, fill = clade_colors[i], alpha = 0.35) +
      geom_cladelabel(
        node = node, label = clade_name, align = TRUE,
        offset = 0.7, barsize = 1.2, fontsize = 4
      )
  }
} else if (length(species_in_tree) == 1) {
  p <- p + geom_tiplab(
    color = clade_colors[i],
    subset = (label %in% species_in_tree),
    fontface = "italic",
    size = 3.2
  )
}
  i <- i + 1
}

# 9. Create a legend for the clades
legend_df <- data.frame(
  Clade = factor(names(clades), levels = names(clades)),
  Color = clade_colors
)

legend_plot <- ggplot(legend_df, aes(x = Clade, y = 1, fill = Clade)) +
  geom_tile() +
  scale_fill_manual(values = clade_colors) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  ) +
  geom_text(aes(label = Clade), hjust = 0, nudge_x = 0.1, size = 4) +
  coord_fixed(ratio = 0.3) +
  xlab(NULL) + ylab(NULL)

# 10. Display tree and legend side by side using patchwork
combined_plot <- p + legend_plot + plot_layout(widths = c(4, 1))

# Print combined plot
print(combined_plot)

# 11. Save the figure as SVG with high quality for publication
output_dir <- p("outputs", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "phylogenetic_tree.svg"), plot = combined_plot,
       width = 12, height = 9, device = "svg", bg = "white")