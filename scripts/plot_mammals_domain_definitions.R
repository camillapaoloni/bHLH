#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Plot mammals: two domain definitions
# -----------------------------------
# - "Predicted" = Ensembl/InterProScan-derived bHLH domain coordinates on target protein
# - "Projected" = Human-domain coordinates projected to query protein via Zoonomia alignments
#
# Outputs
# - Per-orthogroup plot: outputs/orthogroups/domain_definitions_mammals/orthogroup_<HGNC>.svg
# - Focus list + summary: outputs/reports/mammals_domain_definitions_summary.csv
#
# Filters/flags (env vars)
# - BHLH_PROJECT_ROOT: project root (default ".")
# - BHLH_SKIP_EXISTING: "1" (default) to skip existing plots
# - BHLH_LIMIT: integer, plot only first N HGNC symbols
# - BHLH_ONLY: comma-separated HGNC symbols
# - BHLH_DOMAINPLOT_ALL: "1" to plot all HGNCs; default plots only HGNCs with any non-Type match status

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)

skip_existing <- Sys.getenv("BHLH_SKIP_EXISTING", unset = "1") != "0"
limit_n <- suppressWarnings(as.integer(Sys.getenv("BHLH_LIMIT", unset = "0")))
only_hgnc <- Sys.getenv("BHLH_ONLY", unset = "")
only_hgnc <- if (nzchar(only_hgnc)) str_split(only_hgnc, ",", simplify = TRUE) else character(0)
plot_all <- Sys.getenv("BHLH_DOMAINPLOT_ALL", unset = "0") != "0"

in_ensembl <- p("data", "intermediate", "orthologs", "annotated_bHLH_merged_data_with_gene_names.csv")
in_zoo <- p("data", "intermediate", "zoonomia", "Zoonomia_Start_End_final_with_relpos.csv")
in_zoo_map <- p("data", "intermediate", "zoonomia", "Zoonomia_filtered_orthologs.csv")

if (!file.exists(in_ensembl)) stop("Missing input: ", in_ensembl)
if (!file.exists(in_zoo)) stop("Missing input: ", in_zoo)

out_dir <- p("outputs", "orthogroups", "domain_definitions_mammals")
out_reports_dir <- p("outputs", "reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_reports_dir, recursive = TRUE, showWarnings = FALSE)

out_summary_csv <- file.path(out_reports_dir, "mammals_domain_definitions_summary.csv")
out_report_md <- file.path(out_reports_dir, "mammals_domain_definitions_report.md")

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

sanitize <- function(df, start_col, end_col) {
  df %>%
    mutate(
      start_col = suppressWarnings(as.numeric(.data[[start_col]])),
      end_col = suppressWarnings(as.numeric(.data[[end_col]]))
    ) %>%
    filter(is.finite(start_col) & is.finite(end_col)) %>%
    filter(start_col >= 0 & end_col <= 1 & end_col > start_col)
}

# ---- Ensembl: predicted domain ----
df_ens_raw <- readr::read_csv(in_ensembl, show_col_types = FALSE)
ens_required <- c("HGNC symbol", "target_species", "homology_type")
ens_missing <- setdiff(ens_required, names(df_ens_raw))
if (length(ens_missing) > 0) stop("Ensembl input missing columns: ", paste(ens_missing, collapse = ", "))

df_ens <- df_ens_raw %>%
  transmute(
    hgnc = .data[["HGNC symbol"]],
    target_species = as.character(.data[["target_species"]]),
    ensembl_homology_type = as.character(.data[["homology_type"]]),
    ensembl_type_simple = str_replace(ensembl_homology_type, "^ortholog_", ""),
    target_id = if ("target_id" %in% names(df_ens_raw)) as.character(.data[["target_id"]]) else NA_character_,
    target_gene_name = if ("target_gene_name" %in% names(df_ens_raw)) as.character(.data[["target_gene_name"]]) else NA_character_,
    prot_len = if ("Length_target" %in% names(df_ens_raw)) suppressWarnings(as.numeric(.data[["Length_target"]])) else NA_real_,
    rel_start = suppressWarnings(as.numeric(.data[["Rel_start_T"]])),
    rel_end = suppressWarnings(as.numeric(.data[["Rel_end_T"]]))
  ) %>%
  filter(target_species %in% zoonomia_mammals)

# If Rel_* columns are missing, fall back to Start_T/Stop_T/Length_target when available.
if (!("Rel_start_T" %in% names(df_ens_raw)) || !("Rel_end_T" %in% names(df_ens_raw))) {
  if (all(c("Start_T", "Stop_T", "Length_target") %in% names(df_ens_raw))) {
    df_ens <- df_ens_raw %>%
      transmute(
        hgnc = .data[["HGNC symbol"]],
        target_species = as.character(.data[["target_species"]]),
        ensembl_homology_type = as.character(.data[["homology_type"]]),
        ensembl_type_simple = str_replace(ensembl_homology_type, "^ortholog_", ""),
        target_id = if ("target_id" %in% names(df_ens_raw)) as.character(.data[["target_id"]]) else NA_character_,
        target_gene_name = if ("target_gene_name" %in% names(df_ens_raw)) as.character(.data[["target_gene_name"]]) else NA_character_,
        prot_len = suppressWarnings(as.numeric(.data[["Length_target"]])),
        rel_start = suppressWarnings(as.numeric(.data[["Start_T"]]) / as.numeric(.data[["Length_target"]])),
        rel_end = suppressWarnings(as.numeric(.data[["Stop_T"]]) / as.numeric(.data[["Length_target"]]))
      ) %>%
      filter(target_species %in% zoonomia_mammals)
  } else {
    stop("Ensembl input missing Rel_* and Start_T/Stop_T/Length_target fallback.")
  }
}

df_ens_clean <- df_ens %>%
  filter(is.finite(rel_start) & is.finite(rel_end)) %>%
  filter(rel_start >= 0 & rel_end <= 1 & rel_end > rel_start)

df_ens_types <- df_ens_clean %>%
  group_by(hgnc, target_species) %>%
  summarise(
    ensembl_types_simple = paste(sort(unique(na.omit(ensembl_type_simple))), collapse = ";"),
    .groups = "drop"
  )

df_ens_rep <- df_ens_clean %>%
  group_by(hgnc, target_species) %>%
  mutate(
    is_one2one = ensembl_type_simple == "one2one",
    prot_len_num = suppressWarnings(as.numeric(prot_len))
  ) %>%
  arrange(desc(is_one2one), desc(prot_len_num), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    hgnc,
    target_species,
    ens_rel_start = rel_start,
    ens_rel_end = rel_end,
    ens_len = prot_len,
    ens_target_id = target_id,
    ens_gene_name = target_gene_name
  )

# ---- Zoonomia: projected domain ----
df_zoo_raw <- readr::read_csv(in_zoo, show_col_types = FALSE)
zoo_required <- c("HGNC", "full_enst_label", "target_species", "rel_query_start", "rel_query_end", "query_length", "homology_type")
zoo_missing <- setdiff(zoo_required, names(df_zoo_raw))
if (length(zoo_missing) > 0) stop("Zoonomia input missing columns: ", paste(zoo_missing, collapse = ", "))

df_zoo_map <- NULL
if (file.exists(in_zoo_map)) {
  df_zoo_map_raw <- readr::read_csv(in_zoo_map, show_col_types = FALSE)
  zoo_map_required <- c("HGNC", "species", "full_enst_label", "q_gene")
  zoo_map_missing <- setdiff(zoo_map_required, names(df_zoo_map_raw))
  if (length(zoo_map_missing) == 0) {
    df_zoo_map <- df_zoo_map_raw %>%
      transmute(
        HGNC = .data[["HGNC"]],
        target_species = as.character(.data[["species"]]),
        full_enst_label = .data[["full_enst_label"]],
        q_gene = .data[["q_gene"]]
      ) %>%
      distinct()
  }
}

df_zoo <- df_zoo_raw %>%
  mutate(target_species = as.character(target_species)) %>%
  filter(target_species %in% zoonomia_mammals) %>%
  mutate(
    rel_start = suppressWarnings(as.numeric(.data[["rel_query_start"]])),
    rel_end = suppressWarnings(as.numeric(.data[["rel_query_end"]])),
    zoo_type = as.character(.data[["homology_type"]]),
    zoo_len = suppressWarnings(as.numeric(.data[["query_length"]]))
  ) %>%
  filter(is.finite(rel_start) & is.finite(rel_end)) %>%
  filter(rel_start >= 0 & rel_end <= 1 & rel_end > rel_start)

if (!is.null(df_zoo_map)) {
  df_zoo <- df_zoo %>%
    left_join(df_zoo_map, by = c("HGNC", "target_species", "full_enst_label"))
}

df_zoo_types <- df_zoo %>%
  group_by(HGNC, target_species) %>%
  summarise(
    toga_types = paste(sort(unique(na.omit(zoo_type))), collapse = ";"),
    n_q_genes = n_distinct(na.omit(q_gene)),
    n_projected_entries = dplyr::n(),
    .groups = "drop"
  ) %>%
  rename(hgnc = HGNC)

df_zoo_rep <- df_zoo %>%
  group_by(HGNC, target_species) %>%
  arrange(desc(zoo_len), rel_start, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    hgnc = HGNC,
    target_species,
    zoo_rel_start = rel_start,
    zoo_rel_end = rel_end,
    zoo_len,
    zoo_id = if ("q_gene" %in% names(df_zoo)) q_gene else NA_character_
  )

# ---- Status join (per HGNC x species) ----
df_status <- full_join(df_ens_types, df_zoo_types, by = c("hgnc", "target_species")) %>%
  mutate(
    has_ensembl = !is.na(ensembl_types_simple) & ensembl_types_simple != "",
    has_toga = !is.na(toga_types) & toga_types != "",
    type_match = mapply(
      function(a, b) {
        if (is.na(a) || a == "" || is.na(b) || b == "") return(NA)
        aa <- unique(str_split(a, ";", simplify = TRUE))
        bb <- unique(str_split(b, ";", simplify = TRUE))
        length(intersect(aa, bb)) > 0
      },
      ensembl_types_simple,
      toga_types
    ),
    cell_status = case_when(
      has_ensembl & has_toga & type_match == TRUE ~ "Type match",
      has_ensembl & has_toga & type_match == FALSE ~ "Type mismatch",
      has_ensembl & !has_toga ~ "Ensembl only",
      !has_ensembl & has_toga ~ "TOGA only",
      TRUE ~ "None"
    )
  )

df_summary <- df_status %>%
  left_join(df_ens_rep, by = c("hgnc", "target_species")) %>%
  left_join(df_zoo_rep, by = c("hgnc", "target_species")) %>%
  mutate(
    delta_start = ifelse(is.finite(ens_rel_start) & is.finite(zoo_rel_start), abs(ens_rel_start - zoo_rel_start), NA_real_),
    delta_end = ifelse(is.finite(ens_rel_end) & is.finite(zoo_rel_end), abs(ens_rel_end - zoo_rel_end), NA_real_)
  )

write_csv(df_summary, out_summary_csv, na = "")

# ---- Write a short markdown report ----
status_counts <- df_summary %>%
  count(cell_status, name = "n") %>%
  arrange(factor(cell_status, levels = c("Type match", "Type mismatch", "Ensembl only", "TOGA only", "None")))

flagged_orthogroups <- df_summary %>%
  group_by(hgnc) %>%
  summarise(any_flag = any(cell_status != "Type match"), .groups = "drop") %>%
  summarise(n_flagged = sum(any_flag), n_total = dplyr::n(), .groups = "drop")

md <- c(
  "# Mammals domain definitions report",
  "",
  paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("- Ensembl input: `", in_ensembl, "`"),
  paste0("- Zoonomia input: `", in_zoo, "`"),
  "",
  "This report contrasts two domain definitions for the mammals subset:",
  "",
  "- **Predicted**: Ensembl/InterProScan-derived domain coordinates on the target protein.",
  "- **Projected**: human-domain coordinates projected onto the query protein via Zoonomia protein alignments.",
  "",
  "## Selection logic (why some entries do not appear)",
  "",
  "- Plots are restricted to the 8 mammals available in the Zoonomia alignment set.",
  "- For each (HGNC, species), a **single representative** predicted entry is selected (preference: `one2one`, then longest protein).",
  "- For each (HGNC, species), a **single representative** projected entry is selected (longest `query_length`).",
  "- Rows with invalid projected/predicted coordinates are dropped (`NA/Inf`, outside `[0,1]`, or `rel_end <= rel_start`).",
  "",
  "## Outputs",
  "",
  paste0("- Plots: `", out_dir, "/orthogroup_<HGNC>.svg`"),
  paste0("- Summary table: `", out_summary_csv, "`"),
  "",
  "## Status summary (cells)",
  ""
)

md <- c(md, paste0("- ", status_counts$cell_status, ": ", status_counts$n))
md <- c(
  md,
  "",
  "## Orthogroups plotted by default",
  "",
  paste0("- Flagged (any non-Type match cell): ", flagged_orthogroups$n_flagged),
  paste0("- Total in table: ", flagged_orthogroups$n_total),
  "",
  "Set `BHLH_DOMAINPLOT_ALL=1` to plot all orthogroups."
)

writeLines(md, out_report_md)

plot_one <- function(h) {
  out_path <- file.path(out_dir, paste0("orthogroup_", h, ".svg"))
  if (skip_existing && file.exists(out_path)) return(invisible(FALSE))

  df_h <- df_summary %>%
    filter(hgnc == h) %>%
    filter(target_species %in% zoonomia_mammals) %>%
    mutate(
      species_label = species_pretty(target_species),
      species_label = factor(species_label, levels = species_pretty(zoonomia_mammals))
    ) %>%
    arrange(species_label)

  if (nrow(df_h) == 0) return(invisible(FALSE))

  df_h <- df_h %>%
    mutate(
      cell_status = factor(
        cell_status,
        levels = c("Type match", "Type mismatch", "Ensembl only", "TOGA only", "None")
      )
    )

  df_h <- df_h %>%
    mutate(
      ens_id = ifelse(!is.na(ens_target_id) & ens_target_id != "", ens_target_id, NA_character_),
      zoo_region = ifelse(!is.na(zoo_id) & zoo_id != "", zoo_id, NA_character_),
      row_label = dplyr::case_when(
        !is.na(ens_id) & !is.na(zoo_region) ~ paste0(species_label, " | ", ens_id, " | ", zoo_region),
        !is.na(ens_id) ~ paste0(species_label, " | ", ens_id),
        !is.na(zoo_region) ~ paste0(species_label, " | ", zoo_region),
        TRUE ~ as.character(species_label)
      ),
      row_label = factor(row_label, levels = row_label)
    )

  status_colors <- c(
    "Type match" = "#4daf4a",
    "Type mismatch" = "#e41a1c",
    "Ensembl only" = "#377eb8",
    "TOGA only" = "#ff7f00",
    "None" = "#bdbdbd"
  )

  definition_colors <- c(
    "Predicted (Ensembl/InterProScan)" = "#2b8cbe",
    "Projected (Human→Query via Zoonomia aln)" = "#f16913"
  )

  df_status <- df_h %>%
    select(row_label, cell_status) %>%
    distinct()

  p_status <- ggplot(df_status, aes(y = row_label, x = 1, color = cell_status)) +
    geom_point(shape = 15, size = 3, show.legend = TRUE) +
    scale_color_manual(values = status_colors, drop = FALSE) +
    labs(color = "Status") +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      plot.margin = margin(5.5, 0, 5.5, 5.5)
    )

  df_rect <- bind_rows(
    df_h %>%
      filter(is.finite(ens_rel_start) & is.finite(ens_rel_end)) %>%
      transmute(
        row_label,
        definition = "Predicted (Ensembl/InterProScan)",
        xmin = ens_rel_start,
        xmax = ens_rel_end,
        y = as.numeric(row_label),
        ymin = y - 0.30,
        ymax = y - 0.06
      ),
    df_h %>%
      filter(is.finite(zoo_rel_start) & is.finite(zoo_rel_end)) %>%
      transmute(
        row_label,
        definition = "Projected (Human→Query via Zoonomia aln)",
        xmin = zoo_rel_start,
        xmax = zoo_rel_end,
        y = as.numeric(row_label),
        ymin = y + 0.06,
        ymax = y + 0.30
      )
  ) %>%
    mutate(definition = factor(definition, levels = names(definition_colors)))

  p_pos <- ggplot(df_h, aes(y = row_label)) +
    geom_segment(aes(x = 0, xend = 1, yend = row_label), color = "gray85", linewidth = 0.7) +
    geom_rect(
      data = df_rect,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = definition),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.2,
      alpha = 0.85
    ) +
    scale_fill_manual(values = definition_colors, drop = FALSE) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(
      title = paste0("bHLH domain position across mammalian orthologs: ", h),
      subtitle = "Predicted (Ensembl/InterProScan) vs projected (human-domain via Zoonomia alignment)",
      x = "Relative position (0 = N-terminus, 1 = C-terminus)",
      y = NULL,
      fill = "Domain"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 6, face = "bold"),
      panel.grid.major.y = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  df_len <- df_h %>%
    select(row_label, ens_len, zoo_len) %>%
    mutate(y = as.numeric(row_label)) %>%
    pivot_longer(c(ens_len, zoo_len), names_to = "definition", values_to = "prot_len") %>%
    mutate(
      definition = recode(
        definition,
        ens_len = "Predicted (Ensembl/InterProScan)",
        zoo_len = "Projected (Human→Query via Zoonomia aln)"
      ),
      definition = factor(definition, levels = names(definition_colors)),
      prot_len = suppressWarnings(as.numeric(prot_len)),
      y0 = ifelse(definition == "Predicted (Ensembl/InterProScan)", y - 0.18, y + 0.18)
    )

  p_len <- ggplot(df_len, aes(y = y0)) +
    geom_rect(
      data = df_len %>% filter(is.finite(prot_len)),
      aes(
        xmin = 0,
        xmax = prot_len,
        ymin = y0 - 0.12,
        ymax = y0 + 0.12,
        fill = definition
      ),
      color = "black",
      linewidth = 0.15,
      alpha = 0.85,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = df_len %>% filter(is.finite(prot_len)),
      aes(x = prot_len, label = round(prot_len)),
      hjust = -0.15,
      size = 2.2,
      color = "gray20",
      inherit.aes = TRUE
    ) +
    scale_fill_manual(values = definition_colors, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    scale_y_continuous(expand = c(0.02, 0.02)) +
    labs(x = "Protein length (aa)", y = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 0)
    )

  g <- (p_status + p_pos + p_len) +
    plot_layout(widths = c(0.18, 3, 1), guides = "collect") &
    theme(legend.position = "right")

  height_in <- max(4.5, 0.55 * length(levels(df_h$species_label)) + 1.5)
  ggsave(out_path, g, width = 14, height = height_in, units = "in", dpi = 300)
  invisible(TRUE)
}

hgnc_list <- df_summary %>%
  distinct(hgnc) %>%
  pull()

if (!plot_all) {
  flagged <- df_summary %>%
    group_by(hgnc) %>%
    summarise(any_flag = any(cell_status != "Type match"), .groups = "drop") %>%
    filter(any_flag) %>%
    pull(hgnc)
  hgnc_list <- intersect(hgnc_list, flagged)
}

if (length(only_hgnc) > 0) hgnc_list <- intersect(hgnc_list, only_hgnc)
if (!is.na(limit_n) && limit_n > 0) hgnc_list <- head(hgnc_list, limit_n)

message("Orthogroups to plot: ", length(hgnc_list))
message("Output dir: ", out_dir)

for (h in hgnc_list) {
  ok <- plot_one(h)
  if (isTRUE(ok)) message("Wrote: ", file.path(out_dir, paste0("orthogroup_", h, ".svg")))
}

message("Wrote summary: ", out_summary_csv)
message("Wrote report: ", out_report_md)
