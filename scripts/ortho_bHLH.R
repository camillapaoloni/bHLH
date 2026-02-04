library(tidyverse)
library(patchwork)
library(ggplot2)
project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)


# === 1. Load data ===
df <- read_csv(p("data", "intermediate", "orthologs", "annotated_bHLH_merged_data_with_gene_names.csv"))

# Output directory
output_dir <- p("outputs", "orthogroups")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# === 2. Generic plotting function ===
make_domain_plot <- function(data, group_name, group_label, hgnc_symbol = NULL) {
  
  # Create a new column for the Y label
  data <- data %>%
    mutate(
      target_label = ifelse(!is.na(target_gene_name),
         paste0(target_id, " | ", target_gene_name, " (", target_species, ")"),
         paste0(target_id, " (", target_species, ")"))

    )
  
  # Vertical position (Y) for each unique label
  isoform_positions <- data %>%
    distinct(target_label) %>%
    mutate(y_pos = rev(seq_len(n())))
  
  highlight_df <- data %>%
    mutate(
      rel_start = Start_T / Length_target,
      rel_end   = Stop_T / Length_target
    ) %>%
    left_join(isoform_positions, by = "target_label")

  
  highlight_df$y_pos <- as.numeric(highlight_df$target_species)
  plot_title <- paste(group_label, group_name)
  
  p <- ggplot(highlight_df, aes(y = y_pos)) +
    geom_segment(aes(x = 0, xend = 1, yend = y_pos), color = "gray80", size = 0.8) +
    geom_rect(
      aes(xmin = rel_start, xmax = rel_end, ymin = y_pos - 0.3, ymax = y_pos + 0.3),
      fill = "steelblue", color = "black", linewidth = 0.2
    ) +
    scale_y_continuous(
      breaks = isoform_positions$y_pos,
      labels = isoform_positions$target_label,
      expand = expansion(mult = c(0.02, 0.1))
    ) +
    labs(
      title = plot_title,
      x = "Relative Position in Protein",
      y = "Target ID | Gene Name (Species)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 6, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x = element_text(size = 8),
      panel.grid.major.y = element_blank()
    )
  
  return(p)
}


# ====================================================================================
# === 3A. Plot for each SPECIES (without HGNC symbol) ===
# ====================================================================================

df_species_only <- df %>% filter(target_species != "homo_sapiens")

species_list <- unique(df_species_only$target_species)
species_plots <- list()

for (spec in species_list) {
  df_spec <- df_species_only %>% filter(target_species == spec)
  
  if (n_distinct(df_spec$target_id) > 1) {
    p <- make_domain_plot(
      data = df_spec,
      group_name = spec,
      group_label = "Species"
    )
    species_plots[[length(species_plots) + 1]] <- p
  }
}

species_plot_groups <- split(species_plots, ceiling(seq_along(species_plots) / 6))

for (i in seq_along(species_plot_groups)) {
  group_plot <- wrap_plots(species_plot_groups[[i]], ncol = 3) +
    plot_annotation(title = paste("bHLH Domains by Species — Group", i)) &
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
  
  ggsave(
    filename = file.path(output_dir, paste0("species_group_", i, ".png")),
    plot = group_plot,
    width = 16, height = 10, dpi = 300
  )
}

# ====================================================================================
# === 3B. Plot for each ORTHOGROUP (with HGNC symbol) ===
# ====================================================================================
phylo_order <- c(
    'pan_troglodytes',
    'gorilla_gorilla',
    'macaca_mulatta',
    'bos_taurus',
    'canis_lupus_familiaris',
    'rattus_norvegicus',
    'mus_musculus',
    'monodelphis_domestica',
    'gallus_gallus',
    'anolis_carolinensis',
    'lepisosteus_oculatus',
    'oryzias_latipes',
    'danio_rerio',
    'xenopus_tropicalis',
    'drosophila_melanogaster',
    'tribolium_castaneum',
    'anopheles_gambiae',
    'helobdella_robusta',
    'caenorhabditis_elegans',
    'neurospora_crassa',
    'schizosaccharomyces_pombe',
    'saccharomyces_cerevisiae'
)

# Define an orthogroup as all genes sharing the same 'HGNC symbol'
hgnc_list <- df %>%
  filter(!is.na(`HGNC symbol`)) %>%
  distinct(`HGNC symbol`) %>%
  pull()

orthogroup_plots <- list()

df <- df %>% mutate(
  target_species = factor(target_species, levels = phylo_order)
)
for (hgnc in hgnc_list) {
  df_ortho <- df %>% filter(`HGNC symbol` == hgnc)
  
  if (n_distinct(df_ortho$target_id) > 1) {
    p <- make_domain_plot(
      data = df_ortho,
      group_name = hgnc,  # Group name
      group_label = "Orthogroup",
      hgnc_symbol = hgnc
    )
    orthogroup_plots[[length(orthogroup_plots) + 1]] <- p
  }
}

ortho_plot_groups <- split(orthogroup_plots, ceiling(seq_along(orthogroup_plots) / 6))

for (i in seq_along(ortho_plot_groups)) {
  group_plot <- wrap_plots(ortho_plot_groups[[i]], ncol = 3)  &
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))

  ggsave(
    filename = file.path(output_dir, paste0("Orthogroup_", i, ".png")),
    plot = group_plot,
    width = 16, height = 10, dpi = 300
  )
}
