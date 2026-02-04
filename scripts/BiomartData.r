library(biomaRt)

# ---- Paths (relative to project root) ----
project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)
input_ids_path <- p("data", "intermediate", "bHLH_Ensembl_IDs.txt")
lambert_csv_path <- p("data", "raw", "Lambert_bHLH.csv")
output_dir <- p("data", "intermediate", "Metadata_CSVs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read Ensembl Gene IDs ----
# Prefer explicit txt file; fallback to data/raw/Lambert_bHLH.csv if missing.
if (file.exists(input_ids_path)) {
  ensembl_ids <- readLines(input_ids_path)
} else if (file.exists(lambert_csv_path)) {
  lambert_df <- read.csv(lambert_csv_path, stringsAsFactors = FALSE)
  if (!"Ensembl ID" %in% colnames(lambert_df)) {
    stop("data/raw/Lambert_bHLH.csv missing 'Ensembl ID' column.")
  }
  ensembl_ids <- lambert_df[["Ensembl ID"]]
} else {
  stop("No Ensembl IDs found. Provide data/intermediate/bHLH_Ensembl_IDs.txt or data/raw/Lambert_bHLH.csv.")
}

ensembl_ids <- unique(na.omit(ensembl_ids))

# Define BioMart dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

### **Query 1: Transcript-related attributes**
transcript_attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length", "cds_length","description")

transcript_data <- getBM(
  attributes = transcript_attributes,
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
write.csv(transcript_data, file = file.path(output_dir, "Transcript_Attributes.csv"), row.names = FALSE) # nolint

### **Query 2: Protein domains - Pfam**
pfam_attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "pfam", "pfam_start", "pfam_end")  #nolint

pfam_data <- getBM(
  attributes = pfam_attributes,
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
write.csv(pfam_data, file = file.path(output_dir, "Pfam_Domains.csv"), row.names = FALSE)

### **Query 3: Protein domains - InterPro**
interpro_attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "interpro", "interpro_start", "interpro_end", "interpro_short_description", "interpro_description")   # nolint

interpro_data <- getBM(
  attributes = interpro_attributes,
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
write.csv(interpro_data, file = file.path(output_dir, "InterPro_Domains.csv"), row.names = FALSE)

### **Query 4: Protein domains - SMART**
smart_attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "smart", "smart_start", "smart_end")  #nolint

smart_data <- getBM(
  attributes = smart_attributes,
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
write.csv(smart_data, file = file.path(output_dir, "SMART_Domains.csv"), row.names = FALSE)

### **Query 5: Cross-references (UniProt, PDB, etc.)**
xref_attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "uniprotswissprot", "pdb")

xref_data <- getBM(
  attributes = xref_attributes,
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
write.csv(xref_data, file = file.path(output_dir, "External_References.csv"), row.names = FALSE)

### **Query 6: Protein domains - Prosite**
prosite_attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "scanprosite", "scanprosite_start", "scanprosite_end") # nolint

prosite_data <- getBM(
  attributes = prosite_attributes,
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
write.csv(prosite_data, file = file.path(output_dir, "Prosite_Domains.csv"), row.names = FALSE)

### **Query 7: Protein domains - CDD**
cdd_attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "cdd", "cdd_start", "cdd_end") # nolint

cdd_data <- getBM(
  attributes = cdd_attributes,
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
write.csv(cdd_data, file = file.path(output_dir, "CDD_Domains.csv"), row.names = FALSE)

cat("Query complete! Separate CSV files have been saved for each category, including all isoforms.\n")
