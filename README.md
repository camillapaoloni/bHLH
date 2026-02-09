# bHLH Project (Pipeline & Reproducibility Guide)

This repository contains the code and documentation used to analyze bHLH transcription factors, with a focus on:

- bHLH domain position along proteins (human longest isoforms and human alternative isoforms).
- evolutionary conservation of bHLH domain position across orthologs.
- integration with TOGA/Zoonomia resources for mammalian comparisons.
- an RBP / 3'UTR section describing post-transcriptional regulation patterns across bHLH classes.

This README is aligned with the thesis *Materials and Methods* section (Chapter 3) and updated to reflect the latest code refactors in this repository.

## Repository Conventions

- Scripts use paths relative to the project root. If running from a different directory, set `BHLH_PROJECT_ROOT=/path/to/repo`.
- `outputs/` and `data/intermediate/` are intentionally **not tracked** by Git (figures and intermediate tables are reproducible).
- Very large raw FASTA files `data/raw/Zoonomia_protaln/*.fa` are excluded from Git (GitHub 100 MB limit). Keep them locally or store them externally (or via Git LFS).

## Quick Re-Runs (Safe / No Re-Download)

```bash
# Orthogroup bHLH domain-position plots (Ensembl-only + mammal zoom with Zoonomia)
BHLH_PROJECT_ROOT=. bash scripts/run_orthogroup_domain_plots.sh

# RBP / 3'UTR figures
Rscript scripts/rbp/utr_length_boxplots.R
Rscript scripts/rbp/rbp_vs_utr_scatter_and_hist.R
Rscript scripts/rbp/rbp_tsne_pam_clusters.R
Rscript scripts/rbp/rbp_tsne_hclust_clusters.R
Rscript scripts/rbp/rbp_tsne_kmodes_clusters.R
Rscript scripts/rbp/rbp_dendrogram_linkage_methods.R

# Parallel-coordinates plot across species (ortholog midpoints)
Rscript scripts/Parallel_coordinates.r
```

Some notebooks and scripts rely on Ensembl REST/BioMart calls. For fully offline re-runs, prefer scripts that consume already-generated tables under `data/intermediate/`.

## Data Sources (High-Level)

- **Human TF reference list**: Lambert et al. (2018) Human TF catalogue; bHLH subset exported as `data/raw/Lambert_bHLH.csv`.
- **bHLH class annotation**: Ledent/Simionato/Atchley-style classification in `data/raw/LS_classes.csv` (A–E + NC).
- **Ensembl metadata**: BioMart tables (Ensembl Genes 114 / GRCh38.p14 used in the thesis) generated locally under `data/intermediate/Metadata_CSVs/`.
- **Domain annotation**: InterPro entry **IPR011598** (HLH DNA-binding; corresponds to Pfam PF00010).
- **Orthology**: Ensembl Compara (Release 114 in the thesis) + TOGA/Zoonomia supplementation for selected mammals.

## Pipeline (Inputs → Code → Outputs)

### 1) Human bHLH set and BioMart metadata (Thesis 3.1.1)

- Starting point: a CSV containing HGNC symbols + Ensembl Gene IDs for **112 human bHLH TFs** (Lambert-derived; see `data/raw/Lambert_bHLH.csv`).
- Metadata retrieval: `scripts/BiomartData.r`
  - Input: `data/intermediate/bHLH_Ensembl_IDs.txt` (or `data/raw/Lambert_bHLH.csv` column `Ensembl ID`)
  - Output directory: `data/intermediate/Metadata_CSVs/` (transcripts + Pfam/InterPro domain tables + additional annotations)

### 2) Cleaning, merging, and selecting the bHLH InterPro entry (Thesis 3.1.2)

- `notebooks/CSVs_parsing.ipynb`
  - Merges BioMart tables with the 112-gene list and removes entries lacking domain annotations.
  - Selects **IPR011598** as the canonical bHLH InterPro identifier for downstream analyses (thesis rationale: broad coverage; IPR011598 is present across the 112 genes and corresponds to PF00010).

### 3) Protein sequences, translation QC, and bHLH domain extraction (Thesis 3.1.3–3.1.4)

- `notebooks/RoadTo_FASTA.ipynb`
  - Retrieves CDS sequences via Ensembl REST for transcripts annotated with IPR011598, translates them to proteins, and performs basic QC (e.g., CDS length multiple-of-3 checks).
  - Extracts the bHLH domain sequence using `interpro_start` / `interpro_end` coordinates (IPR011598).
  - Outputs are written under `data/intermediate/` (FASTA files and filtered tables).

### 4) Core tables for plotting (human longest isoforms) (Thesis 3.1.5)

- `notebooks/bHLH_classes.ipynb`
  - Integrates: longest isoform per gene, protein length, and bHLH coordinates (IPR011598).
  - Produces `data/intermediate/table_input.csv` (used for the main human bHLH domain position figure).

- Human plot scripts (relative coordinates on [0, 1]):
  - `scripts/bHLH_human_position.r` → `outputs/figures/bHLH_human_position.svg`
  - `scripts/bHLH_human_position_A4.r` → `outputs/figures/bHLH_human_position_A4.svg`

### 5) Isoform-level domain position (Thesis 3.2)

- `scripts/make_bHLH_StartEnd_withISO.py`
  - Reconstructs `data/intermediate/bHLH_StartEnd_withISO.csv` (previously an unsaved bash step).

- Isoform plotting (scaled 0–1 per transcript; “pre / domain / post” concept as in the thesis):
  - Inputs: `data/intermediate/bHLH_StartEnd_withISO.csv`, `data/raw/LS_classes.csv`, `data/raw/Lambert_bHLH.csv`
  - Canonical outputs are the `*_A3.svg` versions produced by:
    - `scripts/isoform_plots_A3.r`
    - `scripts/isoform_plot_merged_improved.R` (merged isoform plot)

### 6) Orthologs across the phylogenetic tree (Thesis 3.3)

The thesis describes two ortholog strategies:

- initial attempt: ProteinOrtho on QfO (not used in the final analysis).
- final approach: **Ensembl Compara** gene-tree reconciliation (Release 114 in the thesis).

In this repository, the Ensembl-based workflow is captured in:

- `notebooks/Orthologs_retrival.ipynb` and `notebooks/Ortho_heatmap.ipynb`
  - Outputs a merged ortholog/domain table:
    - `data/intermediate/orthologs/annotated_bHLH_merged_data_with_gene_names.csv`

- Domain annotation step relies on InterProScan output (thesis uses InterProScan 5.74-105.0).
  - Filtering for the bHLH domain is performed using InterPro identifier **IPR011598**.

### 7) Orthogroup domain-position plots (Thesis 3.4; updated implementation)

- Recommended entry point: `scripts/run_orthogroup_domain_plots.sh`
  - Ensembl-only overview SVGs:
    - `outputs/orthogroups/domain_positions_ensembl/orthogroup_<HGNC>.svg`
  - Mammals-only integrated SVGs (Ensembl + Zoonomia):
    - `outputs/orthogroups/domain_positions_mammals_integrated/orthogroup_<HGNC>.svg`
  - Species labels are formatted as Latin binomials (e.g., `Homo sapiens`) for readability.

### 8) Zoonomia / TOGA supplementation (Thesis 3.5; updated implementation)

- Zoonomia mapping notebook:
  - `notebooks/zoonomia_bHLH_mapping.ipynb` → `data/intermediate/zoonomia/Zoonomia_Start_End_final.csv`

- Relative-position conversion for Zoonomia-mapped coordinates:
  - `scripts/prepare_zoonomia_relpos.py` → `data/intermediate/zoonomia/Zoonomia_Start_End_final_with_relpos.csv`

### 9) Parallel coordinates (ortholog midpoints; updated implementation)

- `scripts/Parallel_coordinates.r` → `outputs/orthologs/figures/parallel_coordinates_midpoints.svg`
- Duplicate hits for the same (gene, species) pair are resolved using:
  - `BHLH_PARCOORD_DUPLICATE_STRATEGY=select_best` (default; chooses the hit closest to the human midpoint)
  - Alternatives: `collapse_median`, `drop_gene`
- A duplicate-hit report is written to:
  - `outputs/orthologs/reports/parallel_coordinates_duplicate_hits.csv`

### 10) RBP / 3'UTR section (project extension)

- Data preparation:
  - `notebooks/3'UTR.ipynb` → `data/intermediate/rbp/Violin_plot_3UTR.csv` and `data/intermediate/rbp/RBP_binary_matrix.csv`

- Figures (all color-coded by bHLH class; outputs under `outputs/rbp/figures/`):
  - `scripts/rbp/utr_length_boxplots.R`
  - `scripts/rbp/rbp_vs_utr_scatter_and_hist.R`
  - `scripts/rbp/rbp_tsne_pam_clusters.R`
  - `scripts/rbp/rbp_tsne_hclust_clusters.R`
  - `scripts/rbp/rbp_tsne_kmodes_clusters.R`
  - `scripts/rbp/rbp_dendrogram_linkage_methods.R`
