#!/usr/bin/env Rscript

# Isoform plots (A3 portrait) with PAS overlay
# ===========================================
#
# This script mirrors the style of `scripts/isoform_plots_A3.r` (stacked "pre / bHLH / post"
# bars scaled to [0, 1] per transcript), and overlays PAS domains when present.
#
# bHLH domain:
# - InterPro: IPR011598 (HLH DNA-binding)
#
# PAS domain:
# - InterPro: IPR000014
# - Displayed as an orange segment on top of the bHLH bar (only when present).
#
# Inputs:
# - data/intermediate/CSVs/bHLH_StartEnd_withISO.csv (or legacy data/intermediate/)
# - data/intermediate/Metadata_CSVs/InterPro_Domains_cleaned.csv
# - data/raw/LS_classes.csv
#
# Outputs:
# - outputs/figures/bHLH_human_transcript_class*_A3_with_PAS.svg

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(tidyr)
})

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)

source(p("scripts", "lib", "bhlh_utils.R"))

# Output directory
output_dir <- p("outputs", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read inputs (isoforms + classes)
df <- read.csv(intermediate_csv_path("bHLH_StartEnd_withISO.csv")) # nolint
df_classes <- read_ls_classes() # nolint

df_classes$Ledent2002.Simionato2007 <- ifelse(
  is.na(df_classes$Ledent2002.Simionato2007) | df_classes$Ledent2002.Simionato2007 == "",
  "NC",
  df_classes$Ledent2002.Simionato2007
)

gene_class_map <- df_classes %>%
  dplyr::select(HGNC.symbol, Ledent2002.Simionato2007) %>%
  distinct()

# Build stacked bins: pre / class / post (sizes are in amino acids).
split(df, paste(df$HGNC.symbol, df$ensembl_transcript_id)) %>%
  lapply(function(tdf) {
    gene_class <- plyr::mapvalues(
      x = tdf$HGNC.symbol,
      from = df_classes$HGNC.symbol,
      to = df_classes$Ledent2002.Simionato2007,
      warn_missing = FALSE
    )

    data.frame(
      HGNC.symbol = tdf$HGNC.symbol,
      ensembl_transcript_id = tdf$ensembl_transcript_id,
      bin = c("pre", gene_class, "post"),
      size = c(
        tdf$interpro_start,
        tdf$interpro_end - tdf$interpro_start,
        tdf$protein_length - tdf$interpro_end
      )
    )
  }) %>%
  do.call("rbind", .) -> bin_sizes

# bHLH class color mapping (used across the repository)
class_colors <- c(
  "A" = "#7AC36A",
  "B" = "#F7CB45",
  "C" = "#9063CD",
  "D" = "#D83F3D",
  "E" = "#4C86C6",
  "NC" = "#cb4289"
)

bin_sizes$bin <- factor(bin_sizes$bin, levels = c("pre", names(class_colors), "post"))
bin_sizes$class <- plyr::mapvalues(
  x = bin_sizes$HGNC.symbol,
  from = df_classes$HGNC.symbol,
  to = df_classes$Ledent2002.Simionato2007,
  warn_missing = FALSE
)

# PAS domain overlay (IPR000014)
in_interpro <- p("data", "intermediate", "Metadata_CSVs", "InterPro_Domains_cleaned.csv")
if (!file.exists(in_interpro)) stop("Missing input: ", in_interpro)

ip <- read.csv(in_interpro, check.names = FALSE, stringsAsFactors = FALSE)
if ("HGNC symbol" %in% names(ip)) ip <- dplyr::rename(ip, HGNC.symbol = `HGNC symbol`)

pas_segments <- ip %>%
  filter(interpro == "IPR000014") %>%
  transmute(
    HGNC.symbol = HGNC.symbol,
    ensembl_transcript_id = ensembl_transcript_id,
    pas_start = as.numeric(interpro_start),
    pas_end = as.numeric(interpro_end)
  ) %>%
  filter(!is.na(ensembl_transcript_id) & is.finite(pas_start) & is.finite(pas_end)) %>%
  group_by(HGNC.symbol, ensembl_transcript_id) %>%
  summarise(pas_start = min(pas_start), pas_end = max(pas_end), .groups = "drop") %>%
  left_join(
    df %>% distinct(HGNC.symbol, ensembl_transcript_id, protein_length),
    by = c("HGNC.symbol", "ensembl_transcript_id")
  ) %>%
  mutate(
    protein_length = as.numeric(protein_length),
    pas_rel_start = pas_start / protein_length,
    pas_rel_end = pas_end / protein_length
  ) %>%
  filter(is.finite(pas_rel_start) & is.finite(pas_rel_end) & protein_length > 0) %>%
  filter(pas_rel_start >= 0 & pas_rel_end <= 1 & pas_rel_end >= pas_rel_start) %>%
  left_join(gene_class_map, by = "HGNC.symbol") %>%
  mutate(
    class = ifelse(is.na(Ledent2002.Simionato2007) | Ledent2002.Simionato2007 == "", "NC", Ledent2002.Simionato2007)
  )

pas_color <- "#FF8C00"

split(bin_sizes, bin_sizes$class) %>%
  lapply(function(pdata) {
    pdata_len <- pdata %>%
      group_by(HGNC.symbol, ensembl_transcript_id) %>%
      summarise(len = sum(size), .groups = "drop")

    cls <- unique(pdata$class)
    if (length(cls) != 1) stop("Expected a single class per split() group, found: ", paste(cls, collapse = ", "))

    pas_cls <- pas_segments %>% filter(class == cls)

    ggplot(pdata, aes(y = ensembl_transcript_id, x = size, fill = bin)) +
      geom_bar(
        color = "black",
        linewidth = 0.25,
        stat = "identity",
        position = position_fill(reverse = TRUE),
        show.legend = FALSE
      ) +
      # PAS overlay: drawn on the same 0-1 scale produced by position_fill().
      geom_segment(
        data = pas_cls,
        inherit.aes = FALSE,
        aes(
          x = pas_rel_start,
          xend = pas_rel_end,
          y = ensembl_transcript_id,
          yend = ensembl_transcript_id
        ),
        color = pas_color,
        linewidth = 2.2,
        alpha = 0.85,
        lineend = "butt"
      ) +
      geom_point(
        data = pdata_len,
        inherit.aes = TRUE,
        mapping = aes(fill = NULL, color = len, size = len, x = 1.05)
      ) +
      geom_text(
        data = pdata_len,
        show.legend = FALSE,
        size = 2,
        inherit.aes = TRUE,
        mapping = aes(fill = NULL, label = len, x = 1.1)
      ) +
      scale_size_continuous(range = c(1, 3)) +
      scale_color_viridis_c(option = "turbo") +
      scale_fill_manual(values = c("pre" = "lightgrey", "post" = "lightgrey", class_colors), guide = guide_none()) +
      facet_grid(
        rows = vars(HGNC.symbol),
        space = "free",
        scales = "free",
        switch = "y"
      ) +
      scale_x_continuous(
        breaks = c(0, 0.25, 0.5, 0.75, 1, 1.1),
        labels = c(0, "25%", "50%", "75%", "100%", "Transcript\nlength (aa)"),
        expand = expansion(mult = c(0, 0.05))
      ) +
      scale_y_discrete(expand = expansion(mult = c(0, 0))) +
      theme(
        strip.placement = "outside",
        panel.spacing.y = unit(0.3, "lines"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5)
      ) +
      labs(
        title = paste0("bHLH isoform plots with PAS overlay (Class ", cls, ")"),
        subtitle = "Orange segment: PAS domain (IPR000014), when present",
        fill = "Gene class",
        x = "Relative position of bHLH/PAS domains",
        y = "Transcript ID",
        color = "Transcript\nlength (aa)",
        size = "Transcript\nlength (aa)"
      )
  }) -> plots_transcripts_human

ggsave(
  file.path(output_dir, "bHLH_human_transcript_classA_A3_with_PAS.svg"),
  plot = plots_transcripts_human$A,
  width = 11.7,
  height = 16.5,
  units = "in"
)
ggsave(
  file.path(output_dir, "bHLH_human_transcript_classB_A3_with_PAS.svg"),
  plot = plots_transcripts_human$B,
  width = 11.7,
  height = 16.5,
  units = "in"
)
ggsave(
  file.path(output_dir, "bHLH_human_transcript_classC_A3_with_PAS.svg"),
  plot = plots_transcripts_human$C,
  width = 11.7,
  height = 16.5,
  units = "in"
)
ggsave(
  file.path(output_dir, "bHLH_human_transcript_classD_A3_with_PAS.svg"),
  plot = plots_transcripts_human$D,
  width = 11.7,
  height = 16.5,
  units = "in"
)
ggsave(
  file.path(output_dir, "bHLH_human_transcript_classE_A3_with_PAS.svg"),
  plot = plots_transcripts_human$E,
  width = 11.7,
  height = 16.5,
  units = "in"
)
ggsave(
  file.path(output_dir, "bHLH_human_transcript_classNC_A3_with_PAS.svg"),
  plot = plots_transcripts_human$NC,
  width = 11.7,
  height = 16.5,
  units = "in"
)

