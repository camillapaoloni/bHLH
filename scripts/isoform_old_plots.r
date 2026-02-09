# === Load Libraries ===
library(ggplot2)
library(tidyverse)
library(patchwork)
library(tidyr)
project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)
source(file.path(project_root, "scripts", "lib", "bhlh_utils.R"))
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



# Output directory
output_dir <- p("outputs", "figures", "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read data files
df <- read.csv(intermediate_csv_path("bHLH_StartEnd_withISO.csv")) # nolint
df_classes <- read.csv(p("data", "raw", "LS_classes.csv")) # nolint
df_classes <- normalize_classes(df_classes)
df_dym <- read.csv(p("data", "raw", "Lambert_bHLH.csv")) # nolint
df_dym <- df_dym %>% select(HGNC.symbol, Binding.mode) %>% distinct()

#MY OTHER PLOT 
# Merge class information
df <- df %>% left_join(df_classes, by = "HGNC.symbol")

# Function to generate plots
generate_plots <- function(gene_df) {
  gene <- unique(gene_df$HGNC.symbol)
  if (n_distinct(gene_df$ensembl_transcript_id) == 1) { # nolint
    message("Skipping gene: ", gene, " because it has only one transcript ID")
    return(NULL)
  }
  gene_class <- unique(gene_df$Ledent2002.Simionato2007)
  fill_color <- ifelse(gene_class %in% names(class_colors), class_colors[gene_class], "gray") # nolint
  # Get binding mode info for the gene
  binding_mode <- df_dym %>%
    filter(HGNC.symbol == gene) %>%
    pull(Binding.mode)

  # If not found, use fallback
  binding_mode <- ifelse(length(binding_mode) == 0, "Unknown", binding_mode)
  
  gene_df <- gene_df %>% mutate(protein_length = as.numeric(protein_length))
  
  # Assign consistent y positions
  isoform_positions <- gene_df %>% 
    distinct(ensembl_transcript_id) %>% 
    mutate(y_pos = rev(seq_len(n())))
  
  # Calculate relative positions
  gene_df <- gene_df %>%
    rowwise() %>%
    mutate(relative_position = list(seq(0, 1, length.out = protein_length))) %>%
    unnest(cols = c(relative_position)) %>%
    mutate(
      highlight = (round(relative_position * (protein_length - 1)) + 1) >= interpro_start &
                  (round(relative_position * (protein_length - 1)) + 1) <= interpro_end
    ) %>%
    ungroup()
  
  highlight_df <- gene_df %>% 
    filter(highlight) %>% 
    group_by(ensembl_transcript_id) %>% 
    summarise(xmin = min(relative_position), xmax = max(relative_position), .groups = "drop") %>%
    left_join(isoform_positions, by = "ensembl_transcript_id")
  
  protein_length_df <- gene_df %>% 
    distinct(ensembl_transcript_id, protein_length) %>% 
    left_join(isoform_positions, by = "ensembl_transcript_id")
  
  rectangle_height <- 0.1
  
  # Plot 1: Domain relative positions
  p1 <- ggplot(highlight_df, aes(y = y_pos)) +  
    geom_rect(aes(xmin = xmin, xmax = xmax, 
                  ymin = y_pos - rectangle_height / 2,
                  ymax = y_pos + rectangle_height / 2),
              fill = fill_color, alpha = 0.8) +  
    scale_y_continuous(
      breaks = isoform_positions$y_pos, 
      labels = isoform_positions$ensembl_transcript_id
    ) +  
    labs(x = "Relative Position in Protein", y = "Transcript ID") +  
    theme_classic() +
    theme(axis.text.y = element_text(size = 8), axis.ticks.y = element_blank()) + 
    ggtitle(paste0(gene, " (Binding mode: ", binding_mode, ")")) + 
    xlim(0, 1)
  
  # Plot 2: Protein length as horizontal bars
  p2 <- ggplot(protein_length_df, aes(y = y_pos, x = protein_length)) +  
    geom_bar(stat = "identity", fill = "#008B8B", width = rectangle_height, orientation = "y") +  
    scale_y_continuous(
      breaks = isoform_positions$y_pos, 
      labels = isoform_positions$ensembl_transcript_id
    ) +  
    labs(x = "Protein Length", y = "") +  
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  return(p1 + p2 + plot_layout(ncol = 2, widths = c(1.5, 1)))
}

message("Code updated: plots are now correctly aligned!")
# Order genes by number of isoforms (descending)
gene_isoform_counts <- df %>% 
  group_by(HGNC.symbol) %>% 
  summarise(num_isoforms = n_distinct(ensembl_transcript_id), .groups = "drop") %>%
  arrange(desc(num_isoforms))

# Generate plots and record isoform counts per gene
plot_list <- list()
iso_counts <- c()  # vector to store isoform counts

genes <- gene_isoform_counts$HGNC.symbol
for (gene in genes) {
  message("Generating plot for: ", gene)
  gene_plot <- generate_plots(df %>% filter(HGNC.symbol == gene))
  if (!is.null(gene_plot)) {
    plot_list[[gene]] <- gene_plot
    iso <- gene_isoform_counts$num_isoforms[gene_isoform_counts$HGNC.symbol == gene]
    iso_counts <- c(iso_counts, iso)
  }
}

# Split the plots into 4 parts with dynamic heights
num_plots <- length(plot_list)
quarter <- ceiling(num_plots / 4)

if (num_plots > 0) {
  for (part in 1:4) {
    start_idx <- (part - 1) * quarter + 1
    end_idx <- min(part * quarter, num_plots)
    if (start_idx > end_idx) break
    
    # Extract plots and isoform counts for this chunk
    chunk_plots <- plot_list[start_idx:end_idx]
    chunk_iso_counts <- iso_counts[start_idx:end_idx]
    
    if (length(chunk_plots) == 0) {
      # If truly no plots in this chunk, skip saving an empty file
      next
    }
    
    # Layout: 2 columns per row
    ncol_val <- 2
    nrow_val <- ceiling(length(chunk_plots) / ncol_val)
    
    # Compute row heights: each row's height is proportional to the maximum isoform count in that row
    row_heights <- numeric(nrow_val)
    for (r in seq_len(nrow_val)) {
      row_start <- (r - 1) * ncol_val + 1
      row_end   <- min(r * ncol_val, length(chunk_plots))
      row_heights[r] <- max(chunk_iso_counts[row_start:row_end])
    }
    
    # Multiply by a scale factor to convert isoform count to plotting height
    scale_factor <- 0.3
    row_heights <- row_heights * scale_factor
    
    # Enforce a minimum row height to avoid "invisible" plots:
    row_heights[row_heights < 2] <- 2
    
    total_height <- sum(row_heights)
    
    # Combine the plots using patchwork with dynamic row heights
    fig <- wrap_plots(
      plotlist = chunk_plots,
      ncol = ncol_val,
      heights = row_heights
    )
    
    out_file <- file.path(
      output_dir,
      paste0("ISO_relposbHLH_", part, ".svg")
    )
    
    ggsave(
      filename = out_file,
      plot     = fig,
      width    = 30, 
      height   = total_height,
      device   = "svg"
    )
    
    
    message("Saved: ", out_file)
  }
  
  message("Plots saved in up to four images!")
} else {
  message("No multi-isoform genes to plot. No images saved.")
}
