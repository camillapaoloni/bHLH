#!/usr/bin/env Rscript

# Parallel-coordinates style plot for ortholog bHLH domain positions.
#
# Each vertical "axis" is a species; each polyline corresponds to one human bHLH TF
# (defined by HGNC symbol) and tracks the midpoint of the bHLH domain along the
# protein sequence (relative position in [0, 1]).
#
# Input:
# - data/intermediate/orthologs/annotated_bHLH_merged_data_with_gene_names.csv
# - data/raw/LS_classes.csv
#
# Output:
# - outputs/orthologs/figures/parallel_coordinates_midpoints.svg
#
# Environment variables:
# - BHLH_PROJECT_ROOT: project root (default ".")
# - BHLH_PARCOORD_ONLY: comma-separated HGNC symbols to plot (e.g. "CLOCK,NCOA1")
# - BHLH_PARCOORD_LIMIT: integer, plot only first N HGNC symbols (after filtering)
# - BHLH_PARCOORD_DUPLICATE_STRATEGY:
#     - "select_best" (default): for each (gene, species) with >1 domain hit, keep the hit whose midpoint is
#       closest to the human midpoint for that gene. Ties are broken by best InterPro/Pfam score (lowest),
#       then by longest domain_length.
#     - "collapse_median": collapse duplicates by median midpoint (robust, but can create a biologically
#       meaningless "average" if the multiple hits are far apart).
#     - "drop_gene": drop any gene that has at least one (gene, species) with >1 hit.
#
# Backwards compatibility:
# - BHLH_PARCOORD_DROP_NONUNIQUE (deprecated): if set to "1" -> strategy="drop_gene"; if "0" -> "collapse_median".

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
source(file.path(project_root, "scripts", "lib", "bhlh_utils.R"))

require_pkgs(c("dplyr", "tidyr", "ggplot2"))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

only_hgnc <- Sys.getenv("BHLH_PARCOORD_ONLY", unset = "")
only_hgnc <- if (nzchar(only_hgnc)) strsplit(only_hgnc, ",")[[1]] else character(0)
only_hgnc <- trimws(only_hgnc)
limit_n <- suppressWarnings(as.integer(Sys.getenv("BHLH_PARCOORD_LIMIT", unset = "0")))

dup_strategy <- Sys.getenv("BHLH_PARCOORD_DUPLICATE_STRATEGY", unset = "select_best")
dup_strategy <- tolower(trimws(dup_strategy))

deprecated_drop_nonunique <- Sys.getenv("BHLH_PARCOORD_DROP_NONUNIQUE", unset = "")
if (nzchar(deprecated_drop_nonunique)) {
  if (deprecated_drop_nonunique != "0") {
    dup_strategy <- "drop_gene"
  } else {
    dup_strategy <- "collapse_median"
  }
}

if (!dup_strategy %in% c("select_best", "collapse_median", "drop_gene")) {
  stop("Invalid BHLH_PARCOORD_DUPLICATE_STRATEGY: ", dup_strategy)
}

# Consistent phylogenetic ordering used across the project (plus human).
species_order <- c(
  "schizosaccharomyces_pombe",
  "saccharomyces_cerevisiae",
  "neurospora_crassa",
  "caenorhabditis_elegans",
  "anopheles_gambiae",
  "drosophila_melanogaster",
  "tribolium_castaneum",
  "helobdella_robusta",
  "danio_rerio",
  "oryzias_latipes",
  "lepisosteus_oculatus",
  "xenopus_tropicalis",
  "gallus_gallus",
  "anolis_carolinensis",
  "monodelphis_domestica",
  "bos_taurus",
  "canis_lupus_familiaris",
  "mus_musculus",
  "rattus_norvegicus",
  "macaca_mulatta",
  "gorilla_gorilla",
  "pan_troglodytes",
  "homo_sapiens"
)

in_path <- p("data", "intermediate", "orthologs", "annotated_bHLH_merged_data_with_gene_names.csv")
if (!file.exists(in_path)) stop("Missing input: ", in_path)

out_dir <- p("outputs", "orthologs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_svg <- file.path(out_dir, "parallel_coordinates_midpoints.svg")

df_raw <- read.csv(in_path, stringsAsFactors = FALSE, check.names = FALSE)
names(df_raw)[names(df_raw) == ""] <- "row_id"

if ("HGNC symbol" %in% names(df_raw)) df_raw <- dplyr::rename(df_raw, HGNC.symbol = `HGNC symbol`)

req <- c("HGNC.symbol", "target_species", "Rel_start_T", "Rel_end_T", "Rel_start_Q", "Rel_end_Q", "homology_type")
missing <- setdiff(req, names(df_raw))
if (length(missing) > 0) stop("Ortholog table missing columns: ", paste(missing, collapse = ", "))

classes <- read_ls_classes()
class_map <- classes %>% select(HGNC.symbol, Ledent2002.Simionato2007) %>% distinct()
pal <- bhlh_class_palette()

# ---- Build long table of midpoints ----
df_orth <- df_raw %>%
  filter(homology_type == "ortholog_one2one") %>%
  mutate(
    target_species = tolower(target_species),
    midpoint = (as.numeric(Rel_start_T) + as.numeric(Rel_end_T)) / 2,
    score = if ("Score" %in% names(df_raw)) suppressWarnings(as.numeric(Score)) else NA_real_,
    domain_length = if ("domain_length" %in% names(df_raw)) suppressWarnings(as.numeric(domain_length)) else NA_real_
  ) %>%
  select(HGNC.symbol, target_species, midpoint, score, domain_length)

df_human <- df_raw %>%
  distinct(HGNC.symbol, Rel_start_Q, Rel_end_Q) %>%
  mutate(
    target_species = "homo_sapiens",
    human_midpoint = (as.numeric(Rel_start_Q) + as.numeric(Rel_end_Q)) / 2
  ) %>%
  select(HGNC.symbol, human_midpoint)

df <- df_orth %>%
  left_join(df_human, by = "HGNC.symbol") %>%
  left_join(class_map, by = "HGNC.symbol") %>%
  mutate(
    Ledent2002.Simionato2007 = ifelse(
      is.na(Ledent2002.Simionato2007) | Ledent2002.Simionato2007 == "",
      "NC",
      Ledent2002.Simionato2007
    ),
    Ledent2002.Simionato2007 = factor(Ledent2002.Simionato2007, levels = c("A", "B", "C", "D", "E", "NC"))
  )

write_dup_report <- function(df_all, chosen_df, out_path) {
  dup <- df_all %>% count(HGNC.symbol, target_species) %>% filter(n > 1)
  if (nrow(dup) == 0) return(invisible(FALSE))

  details <- df_all %>%
    semi_join(dup, by = c("HGNC.symbol", "target_species")) %>%
    group_by(HGNC.symbol, target_species) %>%
    summarise(
      n_hits = n(),
      midpoints = paste(sprintf("%.4f", midpoint), collapse = ";"),
      scores = paste(ifelse(is.na(score), "NA", format(score, scientific = TRUE)), collapse = ";"),
      domain_lengths = paste(ifelse(is.na(domain_length), "NA", sprintf("%.1f", domain_length)), collapse = ";"),
      human_midpoint = unique(human_midpoint)[1],
      .groups = "drop"
    ) %>%
    left_join(
      chosen_df %>%
        select(HGNC.symbol, target_species, midpoint) %>%
        rename(chosen_midpoint = midpoint),
      by = c("HGNC.symbol", "target_species")
    ) %>%
    mutate(strategy = dup_strategy)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write.csv(details, out_path, row.names = FALSE, quote = TRUE)
  invisible(TRUE)
}

resolve_duplicates <- function(df_in) {
  dup <- df_in %>% count(HGNC.symbol, target_species) %>% filter(n > 1)
  if (nrow(dup) == 0) return(df_in)

  if (dup_strategy == "drop_gene") {
    bad_genes <- unique(dup$HGNC.symbol)
    message("Duplicate strategy: drop_gene (dropping ", length(bad_genes), " genes).")
    return(df_in %>% filter(!HGNC.symbol %in% bad_genes))
  }

  if (dup_strategy == "collapse_median") {
    message("Duplicate strategy: collapse_median.")
    return(
      df_in %>%
        group_by(HGNC.symbol, target_species, Ledent2002.Simionato2007, human_midpoint) %>%
        summarise(
          midpoint = median(midpoint, na.rm = TRUE),
          score = suppressWarnings(min(score, na.rm = TRUE)),
          domain_length = suppressWarnings(max(domain_length, na.rm = TRUE)),
          .groups = "drop"
        )
    )
  }

  message("Duplicate strategy: select_best (closest to human midpoint).")
  df_in %>%
    group_by(HGNC.symbol, target_species) %>%
    mutate(
      dist_to_human = abs(midpoint - human_midpoint),
      dist_to_human = ifelse(is.na(dist_to_human), Inf, dist_to_human),
      score_rank = ifelse(is.na(score), Inf, score),
      domain_len_rank = ifelse(is.na(domain_length), -Inf, -domain_length)
    ) %>%
    arrange(dist_to_human, score_rank, domain_len_rank, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(-dist_to_human, -score_rank, -domain_len_rank)
}

df_chosen <- resolve_duplicates(df)
write_dup_report(
  df_all = df,
  chosen_df = df_chosen,
  out_path = p("outputs", "orthologs", "reports", "parallel_coordinates_duplicate_hits.csv")
)
df <- df_chosen

# Add the human reference axis after resolving duplicates.
df <- bind_rows(
  df %>% select(HGNC.symbol, target_species, midpoint, Ledent2002.Simionato2007),
  df_human %>%
    left_join(class_map, by = "HGNC.symbol") %>%
    mutate(
      target_species = "homo_sapiens",
      midpoint = human_midpoint,
      Ledent2002.Simionato2007 = ifelse(
        is.na(Ledent2002.Simionato2007) | Ledent2002.Simionato2007 == "",
        "NC",
        Ledent2002.Simionato2007
      ),
      Ledent2002.Simionato2007 = factor(Ledent2002.Simionato2007, levels = c("A", "B", "C", "D", "E", "NC"))
    ) %>%
    select(HGNC.symbol, target_species, midpoint, Ledent2002.Simionato2007)
)

# Optional gene filtering / limiting for readability.
hgnc_list <- df %>% distinct(HGNC.symbol) %>% pull()
if (length(only_hgnc) > 0) hgnc_list <- intersect(hgnc_list, only_hgnc)
if (!is.na(limit_n) && limit_n > 0) hgnc_list <- head(hgnc_list, limit_n)
df <- df %>% filter(HGNC.symbol %in% hgnc_list)

# Add explicit missing species rows to avoid drawing misleading long segments across gaps.
df <- df %>%
  mutate(target_species = factor(target_species, levels = species_order)) %>%
  select(-Ledent2002.Simionato2007) %>%
  tidyr::complete(HGNC.symbol, target_species, fill = list(midpoint = NA_real_)) %>%
  left_join(class_map, by = "HGNC.symbol") %>%
  mutate(
    Ledent2002.Simionato2007 = ifelse(
      is.na(Ledent2002.Simionato2007) | Ledent2002.Simionato2007 == "",
      "NC",
      Ledent2002.Simionato2007
    ),
    Ledent2002.Simionato2007 = factor(Ledent2002.Simionato2007, levels = c("A", "B", "C", "D", "E", "NC"))
  )

df$target_species_label <- species_pretty(as.character(df$target_species))

df <- df %>%
  group_by(HGNC.symbol) %>%
  arrange(target_species, .by_group = TRUE) %>%
  # Split lines at missing species so we do not draw misleading long segments across gaps.
  mutate(segment_id = cumsum(is.na(midpoint))) %>%
  ungroup()

p_plot <- ggplot(
  df,
  aes(
    x = target_species_label,
    y = midpoint,
    group = interaction(HGNC.symbol, segment_id),
    color = Ledent2002.Simionato2007
  )
) +
  geom_line(linewidth = 0.35, alpha = 0.35, na.rm = TRUE) +
  geom_point(size = 0.6, alpha = 0.55, na.rm = TRUE) +
  scale_color_manual(values = pal, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "bHLH domain midpoint across one-to-one orthologs",
    subtitle = paste0(
      "Each line is a human bHLH TF (HGNC symbol); y = relative bHLH midpoint | duplicate strategy: ",
      dup_strategy
    ),
    x = NULL,
    y = "Relative midpoint (0 = N-terminus, 1 = C-terminus)",
    color = "bHLH class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = "bold"),
    legend.position = "right",
    panel.grid.major.x = element_blank()
  )

ggsave(out_svg, plot = p_plot, width = 14, height = 7, units = "in", dpi = 300)
message("Wrote: ", out_svg)
