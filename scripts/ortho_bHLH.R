#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
})

# Orthogroup bHLH domain plots
# ===========================
#
# Goal
# ----
# For each human bHLH TF (defined by HGNC symbol), create plots that show how the bHLH
# domain start/end positions vary across orthologous proteins.
#
# We generate two plot types per orthogroup:
#
# (A) Ensembl-only overview (homogeneous species set)
#   - Ensembl/InterProScan-based orthologs across the full project species list.
#   - Includes a human reference row derived from the query (human) bHLH coordinates.
#   - Output: outputs/orthogroups/domain_positions_ensembl/orthogroup_<HGNC>.svg
#
# (B) Mammals-only integrated zoom (smaller, heterogeneous sources)
#   - Restricts to mammalian species present in the Zoonomia alignments.
#   - Integrates Ensembl orthologs + Zoonomia-mapped orthologs for those mammals.
#   - Includes the same human reference row.
#   - Output: outputs/orthogroups/domain_positions_mammals_integrated/orthogroup_<HGNC>.svg
#
# Biological assumption
# ---------------------
# The bHLH domain is homologous across orthologs. Coordinates are interpreted in amino acids:
# - Ensembl table provides relative coordinates directly (Rel_start_T / Rel_end_T).
# - Zoonomia coordinates are mapped from human onto orthologs via protein alignments and
#   then converted to relative coordinates using query protein lengths.
#
# Class annotation
# ----------------
# Each orthogroup has a single bHLH class label (A/B/C/D/E/NC) from data/raw/LS_classes.csv.
# We display this as:
# - subtitle: "Class: X"
# - a single-item legend keyed by class color (and a dataset legend for integrated plots).
#
# Reproducibility flags (env vars)
# --------------------------------
# - BHLH_PROJECT_ROOT: project root (default ".")
# - BHLH_USE_ZOONOMIA: "1" (default) to include Zoonomia if available, "0" to skip
# - BHLH_SKIP_EXISTING: "1" (default) to skip plots whose output SVG already exists
# - BHLH_LIMIT: integer, plot only first N HGNC symbols (for quick testing)
# - BHLH_ONLY: comma-separated HGNC symbols to plot (e.g. "CLOCK,NCOA1")

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)

use_zoonomia <- Sys.getenv("BHLH_USE_ZOONOMIA", unset = "1") != "0"
skip_existing <- Sys.getenv("BHLH_SKIP_EXISTING", unset = "1") != "0"
limit_n <- suppressWarnings(as.integer(Sys.getenv("BHLH_LIMIT", unset = "0")))
only_hgnc <- Sys.getenv("BHLH_ONLY", unset = "")
only_hgnc <- if (nzchar(only_hgnc)) str_split(only_hgnc, ",", simplify = TRUE) else character(0)

# Consistent phylogenetic ordering used across the project.
# (This is the 22-species set; the 23rd species in plots is homo_sapiens as the human reference row.)
phylo_order <- c(
  "pan_troglodytes",
  "gorilla_gorilla",
  "macaca_mulatta",
  "bos_taurus",
  "canis_lupus_familiaris",
  "rattus_norvegicus",
  "mus_musculus",
  "monodelphis_domestica",
  "gallus_gallus",
  "anolis_carolinensis",
  "lepisosteus_oculatus",
  "oryzias_latipes",
  "danio_rerio",
  "xenopus_tropicalis",
  "drosophila_melanogaster",
  "tribolium_castaneum",
  "anopheles_gambiae",
  "helobdella_robusta",
  "caenorhabditis_elegans",
  "neurospora_crassa",
  "schizosaccharomyces_pombe",
  "saccharomyces_cerevisiae"
)

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

class_colors <- c(
  "A" = "#7AC36A",
  "B" = "#F7CB45",
  "C" = "#9063CD",
  "D" = "#D83F3D",
  "E" = "#4C86C6",
  "NC" = "#cb4289"
)

species_pretty <- function(x) {
  # Convert Ensembl species IDs (snake_case) to correctly formatted Latin names.
  # Examples:
  # - "homo_sapiens" -> "Homo sapiens"
  # - "canis_lupus_familiaris" -> "Canis lupus familiaris"
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

normalize_classes <- function(df) {
  if ("Ledent2002+Simionato2007" %in% names(df)) {
    df <- dplyr::rename(df, Ledent2002.Simionato2007 = `Ledent2002+Simionato2007`)
  }
  if (!"Ledent2002.Simionato2007" %in% names(df)) {
    stop("Class column not found in LS_classes")
  }
  if ("HGNC symbol" %in% names(df)) {
    df <- dplyr::rename(df, HGNC.symbol = `HGNC symbol`)
  } else if (!"HGNC.symbol" %in% names(df)) {
    stop("HGNC column not found in LS_classes")
  }

  df$Ledent2002.Simionato2007 <- ifelse(
    is.na(df$Ledent2002.Simionato2007) | df$Ledent2002.Simionato2007 == "",
    "NC",
    df$Ledent2002.Simionato2007
  )
  df
}

# ---- Inputs ----
in_ensembl <- p("data", "intermediate", "orthologs", "annotated_bHLH_merged_data_with_gene_names.csv")
if (!file.exists(in_ensembl)) stop("Missing input: ", in_ensembl)

ls_classes_path <- p("data", "raw", "LS_classes.csv")
if (!file.exists(ls_classes_path)) stop("Missing input: ", ls_classes_path)

in_zoo <- p("data", "intermediate", "zoonomia", "Zoonomia_Start_End_final_with_relpos.csv")

# ---- Outputs ----
out_dir_ensembl <- p("outputs", "orthogroups", "domain_positions_ensembl")
out_dir_mammals <- p("outputs", "orthogroups", "domain_positions_mammals_integrated")
dir.create(out_dir_ensembl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_mammals, recursive = TRUE, showWarnings = FALSE)

# ---- Load class map ----
df_classes <- readr::read_csv(ls_classes_path, show_col_types = FALSE) %>% normalize_classes()
class_map <- df_classes %>% select(HGNC.symbol, Ledent2002.Simionato2007) %>% distinct()

orthogroup_class <- function(hgnc) {
  v <- class_map %>% filter(HGNC.symbol == hgnc) %>% pull(Ledent2002.Simionato2007)
  if (length(v) == 0 || is.na(v[[1]]) || v[[1]] == "") return("NC")
  v[[1]]
}

# ---- Load Ensembl ortholog table ----
df_ensembl_raw <- readr::read_csv(in_ensembl, show_col_types = FALSE)

required_cols <- c("HGNC symbol", "target_species", "target_id")
missing_cols <- setdiff(required_cols, names(df_ensembl_raw))
if (length(missing_cols) > 0) {
  stop("Ensembl input missing columns: ", paste(missing_cols, collapse = ", "))
}

df_ensembl <- df_ensembl_raw %>%
  mutate(
    source = "Ensembl",
    hgnc = .data[["HGNC symbol"]],
    target_species = as.character(target_species),
    species_label = species_pretty(target_species),
    target_gene_name = if ("target_gene_name" %in% names(df_ensembl_raw)) .data[["target_gene_name"]] else NA_character_,
    rel_start = suppressWarnings(as.numeric(.data[["Rel_start_T"]])),
    rel_end = suppressWarnings(as.numeric(.data[["Rel_end_T"]])),
    row_id = paste0("Ensembl:", target_id),
    row_label = ifelse(
      !is.na(target_gene_name) & target_gene_name != "",
      paste0(species_label, " | ", target_id, " | ", target_gene_name),
      paste0(species_label, " | ", target_id)
    )
  )

# Fallback: compute relative coordinates if explicit columns are missing.
if (!("Rel_start_T" %in% names(df_ensembl_raw)) || !("Rel_end_T" %in% names(df_ensembl_raw))) {
  if (all(c("Start_T", "Stop_T", "Length_target") %in% names(df_ensembl_raw))) {
    df_ensembl <- df_ensembl %>%
      mutate(
        rel_start = suppressWarnings(as.numeric(.data[["Start_T"]]) / as.numeric(.data[["Length_target"]])),
        rel_end = suppressWarnings(as.numeric(.data[["Stop_T"]]) / as.numeric(.data[["Length_target"]]))
      )
  } else {
    stop("Ensembl input has no Rel_* columns and missing Start_T/Stop_T/Length_target fallback.")
  }
}

df_ensembl <- df_ensembl %>%
  mutate(
    gene_class = vapply(hgnc, orthogroup_class, character(1)),
    gene_class = factor(gene_class, levels = c("A", "B", "C", "D", "E", "NC"))
  ) %>%
  select(hgnc, target_species, row_id, row_label, rel_start, rel_end, gene_class, source)

# ---- Human reference row (adds the 23rd species) ----
df_human_ref <- NULL
if (all(c("Rel_start_Q", "Rel_end_Q", "query_gene", "HGNC symbol") %in% names(df_ensembl_raw))) {
  df_human_ref <- df_ensembl_raw %>%
    mutate(
      source = "Human",
      hgnc = .data[["HGNC symbol"]],
      target_species = "homo_sapiens",
      species_label = species_pretty("homo_sapiens"),
      rel_start = suppressWarnings(as.numeric(.data[["Rel_start_Q"]])),
      rel_end = suppressWarnings(as.numeric(.data[["Rel_end_Q"]])),
      gene_class = vapply(.data[["HGNC symbol"]], orthogroup_class, character(1)),
      gene_class = factor(gene_class, levels = c("A", "B", "C", "D", "E", "NC")),
      row_id = paste0("Human:", .data[["query_gene"]]),
      row_label = paste0(species_label, " | ", .data[["query_gene"]], " | Human")
    ) %>%
    select(hgnc, target_species, row_id, row_label, rel_start, rel_end, gene_class, source) %>%
    distinct(hgnc, .keep_all = TRUE)
} else {
  message("Rel_start_Q/Rel_end_Q not found; human reference rows will be skipped.")
}

# ---- Optional Zoonomia dataset ----
df_zoo <- NULL
if (use_zoonomia && file.exists(in_zoo)) {
  df_zoo_raw <- readr::read_csv(in_zoo, show_col_types = FALSE)
  zoo_required <- c("HGNC", "ENST", "target_species", "rel_query_start", "rel_query_end")
  zoo_missing <- setdiff(zoo_required, names(df_zoo_raw))
  if (length(zoo_missing) == 0) {
    df_zoo <- df_zoo_raw %>%
      mutate(
        source = "Zoonomia",
        hgnc = .data[["HGNC"]],
        target_species = as.character(target_species),
        species_label = species_pretty(target_species),
        rel_start = suppressWarnings(as.numeric(.data[["rel_query_start"]])),
        rel_end = suppressWarnings(as.numeric(.data[["rel_query_end"]])),
        gene_class = vapply(hgnc, orthogroup_class, character(1)),
        gene_class = factor(gene_class, levels = c("A", "B", "C", "D", "E", "NC")),
        row_id = paste0("Zoonomia:", .data[["ENST"]], ":", target_species),
        row_label = paste0(species_label, " | ", .data[["ENST"]], " | Zoonomia")
      ) %>%
      select(hgnc, target_species, row_id, row_label, rel_start, rel_end, gene_class, source)
  } else {
    message("Zoonomia file found but missing columns; skipping. Missing: ", paste(zoo_missing, collapse = ", "))
  }
} else if (use_zoonomia) {
  message("Zoonomia relpos file not found (optional): ", in_zoo)
  message("To generate it: python scripts/prepare_zoonomia_relpos.py --project-root .")
}

sanitize_group <- function(df) {
  df %>%
    filter(is.finite(rel_start) & is.finite(rel_end)) %>%
    filter(rel_start >= 0 & rel_end <= 1 & rel_end >= rel_start)
}

plot_ensembl_full <- function(df_group, hgnc, out_path) {
  df_group <- sanitize_group(df_group)
  if (nrow(df_group) == 0) return(invisible(FALSE))

  # Order rows by species then label; show human row first (top).
  df_group <- df_group %>%
    mutate(target_species = factor(target_species, levels = c("homo_sapiens", phylo_order))) %>%
    arrange(target_species, row_label) %>%
    mutate(row_label = factor(row_label, levels = rev(unique(row_label))))

  gene_class <- unique(as.character(na.omit(df_group$gene_class)))
  gene_class <- if (length(gene_class) == 0) "NC" else gene_class[[1]]

  g <- ggplot(df_group, aes(y = row_label)) +
    geom_segment(aes(x = 0, xend = 1, yend = row_label),
                 color = "gray85", linewidth = 0.7) +
    geom_rect(
      aes(
        xmin = rel_start,
        xmax = rel_end,
        ymin = as.numeric(row_label) - 0.35,
        ymax = as.numeric(row_label) + 0.35,
        fill = gene_class
      ),
      color = "black", linewidth = 0.2, alpha = 0.9
    ) +
    scale_fill_manual(values = class_colors, breaks = gene_class, drop = TRUE) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    labs(
      title = paste0("bHLH domain position across orthologs: ", hgnc),
      subtitle = paste0("Ensembl overview (23 species incl. human) | Class: ", gene_class),
      x = "Relative position (0 = N-terminus, 1 = C-terminus)",
      y = NULL,
      fill = "bHLH class"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 6, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )

  n_rows <- length(levels(df_group$row_label))
  height_in <- max(4, min(40, 0.12 * n_rows + 1))
  ggsave(out_path, g, width = 12, height = height_in, units = "in", dpi = 300)
  invisible(TRUE)
}

plot_mammals_integrated <- function(df_group, hgnc, out_path) {
  df_group <- sanitize_group(df_group)
  if (nrow(df_group) == 0) return(invisible(FALSE))

  df_group <- df_group %>%
    mutate(
      target_species = factor(target_species, levels = c("homo_sapiens", zoonomia_mammals)),
      source = factor(source, levels = c("Human", "Ensembl", "Zoonomia"))
    ) %>%
    arrange(target_species, source, row_label) %>%
    mutate(row_label = factor(row_label, levels = rev(unique(row_label))))

  gene_class <- unique(as.character(na.omit(df_group$gene_class)))
  gene_class <- if (length(gene_class) == 0) "NC" else gene_class[[1]]

  # Dummy point to create a single-item class legend in addition to the dataset legend.
  df_class_key <- tibble(
    x = 0.5,
    y = levels(df_group$row_label)[[1]],
    gene_class = factor(gene_class, levels = c("A", "B", "C", "D", "E", "NC"))
  )

  g <- ggplot(df_group, aes(y = row_label)) +
    geom_segment(aes(x = 0, xend = 1, yend = row_label),
                 color = "gray85", linewidth = 0.7) +
    geom_rect(
      aes(
        xmin = rel_start,
        xmax = rel_end,
        ymin = as.numeric(row_label) - 0.35,
        ymax = as.numeric(row_label) + 0.35,
        fill = source
      ),
      color = "black", linewidth = 0.2, alpha = 0.85
    ) +
    geom_point(
      data = df_class_key,
      aes(x = x, y = y, color = gene_class),
      inherit.aes = FALSE,
      size = 0.1,
      show.legend = TRUE
    ) +
    scale_fill_manual(values = c(Human = "#999999", Ensembl = "#2b8cbe", Zoonomia = "#f16913")) +
    scale_color_manual(values = class_colors, breaks = gene_class, drop = TRUE) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    labs(
      title = paste0("bHLH domain position across mammalian orthologs: ", hgnc),
      subtitle = paste0("Mammals-only zoom (Ensembl + Zoonomia) | Class: ", gene_class),
      x = "Relative position (0 = N-terminus, 1 = C-terminus)",
      y = NULL,
      fill = "Dataset",
      color = "bHLH class"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 6, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )

  n_rows <- length(levels(df_group$row_label))
  height_in <- max(4, min(40, 0.12 * n_rows + 1))
  ggsave(out_path, g, width = 12, height = height_in, units = "in", dpi = 300)
  invisible(TRUE)
}

# ---- HGNC list ----
hgnc_list <- df_ensembl %>% distinct(hgnc) %>% pull()
if (length(only_hgnc) > 0) hgnc_list <- intersect(hgnc_list, only_hgnc)
if (!is.na(limit_n) && limit_n > 0) hgnc_list <- head(hgnc_list, limit_n)

message("Orthogroups: ", length(hgnc_list))
message("Output dirs: ", out_dir_ensembl, " | ", out_dir_mammals)

for (h in hgnc_list) {
  # ---- (A) Ensembl overview (23 species incl. human reference) ----
  out_path_ensembl <- file.path(out_dir_ensembl, paste0("orthogroup_", h, ".svg"))
  if (!(skip_existing && file.exists(out_path_ensembl))) {
    df_g1 <- df_ensembl %>%
      filter(hgnc == h) %>%
      filter(target_species %in% phylo_order)
    if (!is.null(df_human_ref)) {
      df_g1 <- bind_rows(df_g1, df_human_ref %>% filter(hgnc == h))
    }
    ok1 <- plot_ensembl_full(df_g1, h, out_path_ensembl)
    if (isTRUE(ok1)) message("Wrote: ", out_path_ensembl)
  }

  # ---- (B) Mammals-only integrated zoom ----
  out_path_mammals <- file.path(out_dir_mammals, paste0("orthogroup_", h, ".svg"))
  if (!(skip_existing && file.exists(out_path_mammals))) {
    df_g2 <- df_ensembl %>%
      filter(hgnc == h) %>%
      filter(target_species %in% zoonomia_mammals)
    if (!is.null(df_human_ref)) {
      df_g2 <- bind_rows(df_g2, df_human_ref %>% filter(hgnc == h))
    }
    if (!is.null(df_zoo)) {
      df_g2 <- bind_rows(df_g2, df_zoo %>% filter(hgnc == h) %>% filter(target_species %in% zoonomia_mammals))
    }
    ok2 <- plot_mammals_integrated(df_g2, h, out_path_mammals)
    if (isTRUE(ok2)) message("Wrote: ", out_path_mammals)
  }
}
