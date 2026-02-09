#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(viridis)
})

# Human isoform plots with PAS overlay
# ===================================
#
# This script generates class-specific isoform plots (A3 portrait) showing:
# - a 0-1 baseline representing the full protein length for each transcript,
# - the bHLH domain (IPR011598) in the class color,
# - the PAS domain (IPR000014), if present, overlaid in orange.
#
# Inputs:
# - data/intermediate/CSVs/bHLH_StartEnd_withISO.csv (or legacy data/intermediate/...)
# - data/intermediate/Metadata_CSVs/InterPro_Domains_cleaned.csv
# - data/raw/LS_classes.csv
#
# Outputs:
# - outputs/figures/bHLH_human_transcript_class*_A3_with_PAS.svg

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
source(file.path(project_root, "scripts", "lib", "bhlh_utils.R"))

require_pkgs(c("dplyr", "ggplot2", "tidyr", "viridis"))

out_dir <- p("outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

in_iso <- intermediate_csv_path("bHLH_StartEnd_withISO.csv")
in_interpro <- p("data", "intermediate", "Metadata_CSVs", "InterPro_Domains_cleaned.csv")
if (!file.exists(in_interpro)) stop("Missing input: ", in_interpro)

iso <- read.csv(in_iso, stringsAsFactors = FALSE, check.names = FALSE)
names(iso)[names(iso) == ""] <- "row_id"
if ("HGNC symbol" %in% names(iso)) iso <- dplyr::rename(iso, HGNC.symbol = `HGNC symbol`)

req_iso <- c("HGNC.symbol", "ensembl_transcript_id", "protein_length", "interpro_start", "interpro_end")
missing_iso <- setdiff(req_iso, names(iso))
if (length(missing_iso) > 0) stop("bHLH_StartEnd_withISO.csv missing columns: ", paste(missing_iso, collapse = ", "))

iso <- iso %>%
  mutate(
    protein_length = as.numeric(protein_length),
    bhlh_start = as.numeric(interpro_start),
    bhlh_end = as.numeric(interpro_end)
  ) %>%
  select(HGNC.symbol, ensembl_transcript_id, protein_length, bhlh_start, bhlh_end)

classes <- read_ls_classes()
class_map <- classes %>% dplyr::select(HGNC.symbol, Ledent2002.Simionato2007) %>% dplyr::distinct()

iso <- iso %>%
  left_join(class_map, by = "HGNC.symbol") %>%
  mutate(
    class = ifelse(is.na(Ledent2002.Simionato2007) | Ledent2002.Simionato2007 == "", "NC", Ledent2002.Simionato2007),
    class = factor(class, levels = c("A", "B", "C", "D", "E", "NC"))
  ) %>%
  select(-Ledent2002.Simionato2007)

# PAS domains (IPR000014) per transcript: keep min start, max end.
interpro <- read.csv(in_interpro, stringsAsFactors = FALSE, check.names = FALSE)
names(interpro)[names(interpro) == ""] <- "row_id"
if ("HGNC symbol" %in% names(interpro)) interpro <- dplyr::rename(interpro, HGNC.symbol = `HGNC symbol`)

req_ip <- c("HGNC.symbol", "ensembl_transcript_id", "interpro", "interpro_start", "interpro_end")
missing_ip <- setdiff(req_ip, names(interpro))
if (length(missing_ip) > 0) stop("InterPro_Domains_cleaned.csv missing columns: ", paste(missing_ip, collapse = ", "))

pas <- interpro %>%
  filter(interpro == "IPR000014") %>%
  mutate(
    pas_start = as.numeric(interpro_start),
    pas_end = as.numeric(interpro_end)
  ) %>%
  filter(!is.na(ensembl_transcript_id) & is.finite(pas_start) & is.finite(pas_end)) %>%
  group_by(ensembl_transcript_id) %>%
  summarise(pas_start = min(pas_start, na.rm = TRUE), pas_end = max(pas_end, na.rm = TRUE), .groups = "drop")

df <- iso %>%
  left_join(pas, by = "ensembl_transcript_id") %>%
  mutate(
    bhlh_rel_start = bhlh_start / protein_length,
    bhlh_rel_end = bhlh_end / protein_length,
    pas_rel_start = pas_start / protein_length,
    pas_rel_end = pas_end / protein_length
  ) %>%
  filter(is.finite(bhlh_rel_start) & is.finite(bhlh_rel_end) & protein_length > 0) %>%
  filter(bhlh_rel_start >= 0 & bhlh_rel_end <= 1 & bhlh_rel_end >= bhlh_rel_start)

class_colors <- bhlh_class_palette()
pas_color <- "#FF8C00"
fill_values <- c(class_colors, "PAS domain" = pas_color)

make_plot <- function(class_name) {
  d <- df %>% filter(class == class_name)
  if (nrow(d) == 0) return(NULL)

  # Order transcripts within each gene by decreasing length.
  d <- d %>%
    group_by(HGNC.symbol) %>%
    arrange(desc(protein_length), ensembl_transcript_id, .by_group = TRUE) %>%
    ungroup()

  # Use a global factor to support geom_rect y placement with facets.
  d <- d %>% mutate(row_label = factor(ensembl_transcript_id, levels = rev(unique(ensembl_transcript_id))))

  # One row per transcript for length dot.
  d_len <- d %>%
    distinct(HGNC.symbol, ensembl_transcript_id, protein_length, row_label)

  g <- ggplot(d, aes(y = row_label)) +
    geom_segment(aes(x = 0, xend = 1, yend = row_label), color = "gray85", linewidth = 0.7) +
    # bHLH domain in class color
    geom_rect(
      aes(
        xmin = bhlh_rel_start,
        xmax = bhlh_rel_end,
        ymin = as.numeric(row_label) - 0.35,
        ymax = as.numeric(row_label) + 0.35,
        fill = class
      ),
      color = "black", linewidth = 0.2, alpha = 0.95,
      show.legend = TRUE
    ) +
    # PAS overlay (orange), if present
    geom_rect(
      data = d %>% filter(is.finite(pas_rel_start) & is.finite(pas_rel_end)) %>%
        filter(pas_rel_start >= 0 & pas_rel_end <= 1 & pas_rel_end >= pas_rel_start),
      inherit.aes = FALSE,
      aes(
        xmin = pas_rel_start,
        xmax = pas_rel_end,
        ymin = as.numeric(row_label) - 0.20,
        ymax = as.numeric(row_label) + 0.20,
        fill = "PAS domain"
      ),
      color = "black", linewidth = 0.2, alpha = 0.85,
      show.legend = TRUE
    ) +
    # length dot (right side)
    geom_point(
      data = d_len,
      inherit.aes = FALSE,
      aes(x = 1.05, y = row_label, color = protein_length, size = protein_length),
      alpha = 0.8
    ) +
    geom_text(
      data = d_len,
      inherit.aes = FALSE,
      aes(x = 1.10, y = row_label, label = round(protein_length)),
      size = 2.0
    ) +
    scale_fill_manual(values = fill_values, breaks = c(as.character(class_name), "PAS domain"), drop = TRUE) +
    scale_color_viridis_c(option = "turbo", name = "Transcript length (aa)") +
    scale_size_continuous(range = c(1, 3), name = "Transcript length (aa)") +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1, 1.10),
      labels = c("0", "25%", "50%", "75%", "100%", "Transcript\nlength (aa)"),
      limits = c(0, 1.12),
      expand = expansion(mult = c(0, 0.02))
    ) +
    facet_grid(rows = vars(HGNC.symbol), space = "free", scales = "free", switch = "y") +
    labs(
      title = paste0("bHLH domain position across coding isoforms (Class ", class_name, ")"),
      x = "Relative position of bHLH/PAS domains",
      y = "Transcript ID",
      fill = "Feature"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      panel.spacing.y = unit(0.25, "lines"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 6, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )

  g
}

classes_to_plot <- c("A", "B", "C", "D", "E", "NC")
plots <- lapply(classes_to_plot, make_plot)
names(plots) <- classes_to_plot

for (cls in classes_to_plot) {
  g <- plots[[cls]]
  if (is.null(g)) next
  out_svg <- file.path(out_dir, paste0("bHLH_human_transcript_class", cls, "_A3_with_PAS.svg"))
  ggsave(out_svg, plot = g, width = 11.7, height = 16.5, units = "in", dpi = 300)
  message("Wrote: ", out_svg)
}

