#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Report: Ensembl vs TOGA (Zoonomia) orthology concordance for mammals
# --------------------------------------------------------------------
# Produces:
# - outputs/reports/ensembl_vs_toga_mammals_mapping.csv
# - outputs/reports/ensembl_vs_toga_mammals_report.md
# - outputs/figures/ensembl_vs_toga_mammals_matrix.svg

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)

in_ensembl <- p("data", "intermediate", "orthologs", "annotated_bHLH_merged_data_with_gene_names.csv")
in_toga <- p("data", "intermediate", "orthologs", "TOGA_orthologs_allSpecies.csv")

out_reports_dir <- p("outputs", "reports")
out_figures_dir <- p("outputs", "figures")
dir.create(out_reports_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figures_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- file.path(out_reports_dir, "ensembl_vs_toga_mammals_mapping.csv")
out_md <- file.path(out_reports_dir, "ensembl_vs_toga_mammals_report.md")
out_svg <- file.path(out_figures_dir, "ensembl_vs_toga_mammals_matrix.svg")

if (!file.exists(in_ensembl)) stop("Missing input: ", in_ensembl)
if (!file.exists(in_toga)) stop("Missing input: ", in_toga)

# Mammal subset (species represented in the Zoonomia alignment set).
zoonomia_mammals <- c(
  "pan_troglodytes",
  "gorilla_gorilla",
  "macaca_mulatta",
  "bos_taurus",
  "canis_lupus_familiaris",
  "rattus_norvegicus",
  "mus_musculus",
  "monodelphis_domestica"
)

species_pretty <- function(x) {
  x <- stringr::str_replace_all(as.character(x), "_", " ")
  parts <- stringr::str_split(x, "\\s+")
  vapply(
    parts,
    function(w) {
      w <- w[w != ""]
      if (length(w) == 0) return("")
      w <- stringr::str_to_lower(w)
      w[1] <- stringr::str_to_title(w[1])
      paste(w, collapse = " ")
    },
    character(1)
  )
}

df_ensembl_raw <- readr::read_csv(in_ensembl, show_col_types = FALSE)
ens_required <- c("HGNC symbol", "target_species", "target_id", "homology_type")
ens_missing <- setdiff(ens_required, names(df_ensembl_raw))
if (length(ens_missing) > 0) stop("Ensembl input missing columns: ", paste(ens_missing, collapse = ", "))

df_ensembl <- df_ensembl_raw %>%
  transmute(
    hgnc = .data[["HGNC symbol"]],
    target_species = as.character(.data[["target_species"]]),
    target_id = as.character(.data[["target_id"]]),
    target_gene_name = if ("target_gene_name" %in% names(df_ensembl_raw)) as.character(.data[["target_gene_name"]]) else NA_character_,
    homology_type = as.character(.data[["homology_type"]])
  ) %>%
  filter(target_species %in% zoonomia_mammals)

df_ensembl <- df_ensembl %>%
  mutate(ensembl_type_simple = str_replace(homology_type, "^ortholog_", ""))

df_ensembl_sum <- df_ensembl %>%
  group_by(hgnc, target_species) %>%
  summarise(
    ensembl_target_ids = paste(sort(unique(na.omit(target_id))), collapse = ";"),
    ensembl_gene_names = paste(sort(unique(na.omit(target_gene_name[target_gene_name != ""]))), collapse = ";"),
    ensembl_homology_types = paste(sort(unique(na.omit(homology_type))), collapse = ";"),
    ensembl_types_simple = paste(sort(unique(na.omit(ensembl_type_simple))), collapse = ";"),
    .groups = "drop"
  )

df_toga_raw <- readr::read_csv(in_toga, show_col_types = FALSE)
toga_required <- c("HGNC symbol", "target_species", "q_gene", "q_transcript", "homology_type")
toga_missing <- setdiff(toga_required, names(df_toga_raw))
if (length(toga_missing) > 0) stop("TOGA input missing columns: ", paste(toga_missing, collapse = ", "))

df_toga <- df_toga_raw %>%
  transmute(
    hgnc = .data[["HGNC symbol"]],
    target_species = as.character(.data[["target_species"]]),
    q_gene = as.character(.data[["q_gene"]]),
    q_transcript = as.character(.data[["q_transcript"]]),
    homology_type = as.character(.data[["homology_type"]])
  ) %>%
  filter(target_species %in% zoonomia_mammals)

df_toga_sum <- df_toga %>%
  group_by(hgnc, target_species) %>%
  summarise(
    toga_q_genes = paste(sort(unique(na.omit(q_gene[q_gene != ""]))), collapse = ";"),
    toga_q_transcripts = paste(sort(unique(na.omit(q_transcript[q_transcript != ""]))), collapse = ";"),
    toga_homology_types = paste(sort(unique(na.omit(homology_type))), collapse = ";"),
    .groups = "drop"
  )

df_map <- full_join(df_ensembl_sum, df_toga_sum, by = c("hgnc", "target_species"))

split_types <- function(x) {
  x <- as.character(x)
  if (is.na(x) || x == "") return(character(0))
  unique(str_split(x, ";", simplify = TRUE))
}

df_map <- df_map %>%
  mutate(
    has_ensembl = !is.na(ensembl_target_ids) & ensembl_target_ids != "",
    has_toga = !is.na(toga_q_genes) & toga_q_genes != "",
    type_match = mapply(
      function(a, b) {
        aa <- split_types(a)
        bb <- split_types(b)
        if (length(aa) == 0 || length(bb) == 0) return(NA)
        length(intersect(aa, bb)) > 0
      },
      ensembl_types_simple,
      toga_homology_types
    ),
    cell_status = case_when(
      has_ensembl & has_toga & type_match == TRUE ~ "Type match",
      has_ensembl & has_toga & type_match == FALSE ~ "Type mismatch",
      has_ensembl & !has_toga ~ "Ensembl only",
      !has_ensembl & has_toga ~ "TOGA only",
      TRUE ~ "None"
    )
  ) %>%
  arrange(hgnc, target_species)

write_csv(df_map, out_csv, na = "")

# ---- Figure: matrix ----
df_plot <- df_map %>%
  mutate(
    target_species = factor(target_species, levels = zoonomia_mammals),
    species_label = species_pretty(target_species),
    cell_status = factor(
      cell_status,
      levels = c("Type match", "Type mismatch", "Ensembl only", "TOGA only", "None")
    )
  )

g <- ggplot(df_plot, aes(x = species_label, y = hgnc, fill = cell_status)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "Type match" = "#4daf4a",
      "Type mismatch" = "#e41a1c",
      "Ensembl only" = "#377eb8",
      "TOGA only" = "#ff7f00",
      "None" = "#f0f0f0"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Ensembl vs TOGA orthology concordance (mammals)",
    subtitle = "Cell shows presence and (when both present) homology type concordance",
    x = NULL,
    y = NULL,
    fill = "Status"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

height_in <- max(6, min(40, 0.18 * length(unique(df_plot$hgnc)) + 2))
ggsave(out_svg, g, width = 12, height = height_in, units = "in", dpi = 300)

# ---- Markdown summary ----
summary_counts <- df_map %>%
  count(cell_status, name = "n") %>%
  arrange(factor(cell_status, levels = c("Type match", "Type mismatch", "Ensembl only", "TOGA only", "None")))

mismatches <- df_map %>%
  filter(cell_status == "Type mismatch") %>%
  transmute(
    hgnc,
    target_species,
    ensembl_target_ids,
    ensembl_types_simple,
    toga_q_genes,
    toga_homology_types
  )

ensembl_only <- df_map %>%
  filter(cell_status == "Ensembl only") %>%
  transmute(
    hgnc,
    target_species,
    ensembl_target_ids,
    ensembl_types_simple
  )

toga_only <- df_map %>%
  filter(cell_status == "TOGA only") %>%
  transmute(
    hgnc,
    target_species,
    toga_q_genes,
    toga_homology_types
  )

md_lines <- c(
  "# Ensembl vs TOGA (Zoonomia) — Mammals mapping report",
  "",
  paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("- Ensembl input: `", in_ensembl, "`"),
  paste0("- TOGA input: `", in_toga, "`"),
  paste0("- Output table: `", out_csv, "`"),
  paste0("- Status matrix: `", out_svg, "`"),
  "",
  "## Summary (cells)",
  ""
)

md_lines <- c(
  md_lines,
  paste0("- ", summary_counts$cell_status, ": ", summary_counts$n)
)

md_lines <- c(md_lines, "", "## Type mismatches (details)", "")
if (nrow(mismatches) == 0) {
  md_lines <- c(md_lines, "- None")
} else {
  for (i in seq_len(nrow(mismatches))) {
    r <- mismatches[i, ]
    md_lines <- c(
      md_lines,
      paste0(
        "- **", r$hgnc, "** | ", r$target_species,
        " — Ensembl: ", ifelse(r$ensembl_target_ids == "", "NA", r$ensembl_target_ids),
        " (", ifelse(r$ensembl_types_simple == "", "NA", r$ensembl_types_simple), ")",
        " vs TOGA q_gene: ", ifelse(r$toga_q_genes == "", "NA", r$toga_q_genes),
        " (", ifelse(r$toga_homology_types == "", "NA", r$toga_homology_types), ")"
      )
    )
  }
}

md_lines <- c(md_lines, "", "## Ensembl-only (no TOGA entry)", "")
if (nrow(ensembl_only) == 0) {
  md_lines <- c(md_lines, "- None")
} else {
  for (i in seq_len(nrow(ensembl_only))) {
    r <- ensembl_only[i, ]
    md_lines <- c(
      md_lines,
      paste0(
        "- **", r$hgnc, "** | ", r$target_species,
        " — Ensembl: ", ifelse(r$ensembl_target_ids == "", "NA", r$ensembl_target_ids),
        " (", ifelse(r$ensembl_types_simple == "", "NA", r$ensembl_types_simple), ")"
      )
    )
  }
}

md_lines <- c(md_lines, "", "## TOGA-only (no Ensembl entry)", "")
if (nrow(toga_only) == 0) {
  md_lines <- c(md_lines, "- None")
} else {
  for (i in seq_len(nrow(toga_only))) {
    r <- toga_only[i, ]
    md_lines <- c(
      md_lines,
      paste0(
        "- **", r$hgnc, "** | ", r$target_species,
        " — TOGA q_gene: ", ifelse(r$toga_q_genes == "", "NA", r$toga_q_genes),
        " (", ifelse(r$toga_homology_types == "", "NA", r$toga_homology_types), ")"
      )
    )
  }
}

md_lines <- c(
  md_lines,
  "",
  "## Related outputs",
  "",
  "- Mammals domain-definition plots (predicted vs projected): `outputs/orthogroups/domain_definitions_mammals/`",
  "- Summary table (per HGNC × species): `outputs/reports/mammals_domain_definitions_summary.csv`"
)

writeLines(md_lines, out_md)

message("Wrote: ", out_csv)
message("Wrote: ", out_md)
message("Wrote: ", out_svg)
