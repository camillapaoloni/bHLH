#!/usr/bin/env Rscript

# Shared utilities for the bHLH project.
# Keep this file dependency-light: it should be safe to source from any script.

suppressPackageStartupMessages({
  library(dplyr)
})

project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)

normalize_ls_classes <- function(df) {
  # Normalize the LS class table column names across the different variants found in the project.
  if ("Ledent2002+Simionato2007" %in% names(df)) {
    df <- dplyr::rename(df, Ledent2002.Simionato2007 = `Ledent2002+Simionato2007`)
  }
  if ("HGNC symbol" %in% names(df)) {
    df <- dplyr::rename(df, HGNC.symbol = `HGNC symbol`)
  }

  if (!("Ledent2002.Simionato2007" %in% names(df))) stop("Class column not found in LS_classes")
  if (!("HGNC.symbol" %in% names(df))) stop("HGNC column not found in LS_classes")

  df$Ledent2002.Simionato2007 <- ifelse(
    is.na(df$Ledent2002.Simionato2007) | df$Ledent2002.Simionato2007 == "",
    "NC",
    df$Ledent2002.Simionato2007
  )
  df
}

read_ls_classes <- function(ls_classes_path = p("data", "raw", "LS_classes.csv")) {
  if (!file.exists(ls_classes_path)) stop("Missing LS classes file: ", ls_classes_path)
  df <- read.csv(ls_classes_path, stringsAsFactors = FALSE, check.names = FALSE)
  normalize_ls_classes(df)
}

bhlh_class_palette <- function() {
  # Consistent class colors across the whole project.
  c(
    "A" = "#7AC36A",
    "B" = "#F7CB45",
    "C" = "#9063CD",
    "D" = "#D83F3D",
    "E" = "#4C86C6",
    "NC" = "#7f7f7f"
  )
}

species_pretty <- function(x) {
  # Convert Ensembl species IDs (snake_case) to correctly formatted Latin names.
  # Examples:
  # - "homo_sapiens" -> "Homo sapiens"
  # - "canis_lupus_familiaris" -> "Canis lupus familiaris"
  x <- gsub("_", " ", as.character(x))
  parts <- strsplit(x, "\\s+")
  vapply(
    parts,
    function(w) {
      w <- w[w != ""]
      if (length(w) == 0) return("")
      w <- tolower(w)
      w[1] <- paste0(toupper(substr(w[1], 1, 1)), substr(w[1], 2, nchar(w[1])))
      paste(w, collapse = " ")
    },
    character(1)
  )
}

require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      "\nInstall with: install.packages(c(", paste(sprintf('\"%s\"', missing), collapse = ", "), "))",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
