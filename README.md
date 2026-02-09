# bHLH Project (Pipeline & Reproducibility Guide)

This repository contains the code and documentation used to analyze bHLH transcription factors with an explicit focus on *protein architecture* and the evolutionary constraints that may shape it.

The biological motivation (as described in the thesis) is that co-translational assembly (Co-TA) can impose constraints on where interaction domains are located along a protein: an N-terminal interaction domain may favor Co-TA, whereas a C-terminal one may make it less likely. In this project, the interaction domain of interest is the bHLH DNA-binding/dimerization domain.

## Core Objectives (Thesis Chapter 3)

The project is organized around three main questions:

- **Human architecture (longest isoforms):** where is the bHLH domain located along each of the 112 human bHLH TFs?
- **Isoform plasticity:** does alternative splicing shift the relative position of the bHLH domain across coding isoforms?
- **Evolutionary conservation:** is the domain position conserved across orthologs spanning a broad phylogenetic range?

In addition, this repository includes an **RBP / 3'UTR extension**, which explores post-transcriptional regulation patterns across bHLH classes.

This README is aligned with the thesis *Materials and Methods* section (Chapter 3) and updated to reflect the latest code refactors in this repository.

## Project Layout and Reproducibility

- Scripts use paths relative to the project root. If running from a different directory, set `BHLH_PROJECT_ROOT=/path/to/repo`.
- `outputs/` and `data/intermediate/` are intentionally **not tracked** by Git (figures and intermediate tables are reproducible).
- Very large raw FASTA files `data/raw/Zoonomia_protaln/*.fa` are excluded from Git (GitHub 100 MB limit). Keep them locally or store them externally (or via Git LFS). This data can be downloaded from [Zoonomia](https://zoonomiaproject.org/the-data/).

## Key Data Sources (Methods Summary)

- **Human TF reference list:** Lambert et al. (2018) Human TF catalogue; bHLH subset exported as `data/raw/Lambert_bHLH.csv` (HGNC symbols and Ensembl Gene IDs).
- **Class annotation:** `data/raw/LS_classes.csv` (A-E + NC), used consistently as the color code across human, isoform, ortholog, and RBP plots.
- **Ensembl metadata:** BioMart tables (Ensembl Genes 114 / GRCh38.p14 used in the thesis) generated locally under `data/intermediate/Metadata_CSVs/`.
- **Domain identifier:** InterPro entry **IPR011598** (Helix-loop-helix DNA-binding domain; corresponds to Pfam PF00010).
- **Orthology:** Ensembl Compara (Release 114 in the thesis) + TOGA/Zoonomia supplementation for selected mammals.

## Coordinate Definitions (Important)

Across the whole project, domain coordinates are treated in **amino acids**. For a given protein of length `L`:

- absolute domain coordinates: `start`, `end` (1-based coordinates from BioMart/InterProScan outputs)
- relative coordinates: `rel_start = start / L` and `rel_end = end / L`
- relative midpoint: `midpoint = (rel_start + rel_end) / 2`

Using relative coordinates allows comparisons across proteins with different lengths.

## Quick Re-Runs (Safe / No Re-Download)

```bash
# Orthogroup bHLH domain-position plots (Ensembl-only + mammal zoom with Zoonomia)
BHLH_PROJECT_ROOT=. bash scripts/run_orthogroup_domain_plots.sh

# RBP / 3'UTR figures (extension)
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

## Pipeline (Methods-Aligned; Inputs -> Code -> Outputs)

### 3.1 Relative position of the bHLH domain in human genes

#### 3.1.1 Data retrieval

The starting point is a CSV containing HGNC symbols and Ensembl Gene IDs for **112 human bHLH TFs** (Lambert-derived; see `data/raw/Lambert_bHLH.csv`). Transcript and domain metadata are retrieved from Ensembl BioMart (Ensembl Genes 114; GRCh38.p14) using:

- `scripts/BiomartData.r`
  - Input: `data/intermediate/bHLH_Ensembl_IDs.txt` (or `data/raw/Lambert_bHLH.csv` column `Ensembl ID`)
  - Output directory: `data/intermediate/Metadata_CSVs/`

BioMart attributes used in the thesis include Ensembl gene/transcript IDs, transcript/CDS lengths, and structural domain annotations (Pfam and InterPro) with start/end coordinates.

#### 3.1.2 Filtering and selection of the bHLH InterPro entry

In the thesis, multiple InterPro entries related to the HLH superfamily exist, but **IPR011598** was selected as the canonical identifier for downstream analyses because it captures the bHLH DNA-binding domain consistently across the 112-gene set and provides broad coverage in BioMart outputs.

In this repository the cleaning and normalization step is captured in:

- `notebooks/CSVs_parsing.ipynb`
  - merges BioMart tables with the 112-gene list (ENSG key)
  - removes rows with missing domain annotations
  - keeps the bHLH domain coordinates using InterPro identifier IPR011598

#### 3.1.3 Protein sequence retrieval and translation QC

In the thesis, CDS sequences for all transcripts annotated with IPR011598 are fetched programmatically via the Ensembl REST API and translated to amino acids. A small number of CDS sequences not divisible by 3 (often marked as incomplete) are excluded to avoid translation artifacts.

In this repository these steps are implemented in:

- `notebooks/RoadTo_FASTA.ipynb`
  - outputs FASTA files under `data/intermediate/` (CDS, protein, longest isoforms, and extracted bHLH domain FASTA)

#### 3.1.4 Multiple sequence alignment and tree topology QC

The thesis uses an MSA of extracted bHLH domain sequences as a quality control step, comparing the inferred clustering/topology against the reference bHLH tree from the Human TF database. Clustal Omega (ClustalO) is used to mirror the methodology of the reference study.

In this repository, the alignment and domain extraction logic is documented in:

- `notebooks/RoadTo_FASTA.ipynb`

#### 3.1.5 Human longest-isoform plot (relative bHLH position)

The thesis describes a composite figure built in R (ggplot2 + tidyverse + patchwork) combining:

- a central panel showing the bHLH domain as a rectangle placed on a 0-1 protein scale
- a length bar chart (longest isoform length)
- an isoform count bar chart (transcript diversity)

In this repository:

- `notebooks/bHLH_classes.ipynb` generates the core table `data/intermediate/table_input.csv`
- plotting scripts:
  - `scripts/bHLH_human_position.r` -> `outputs/figures/bHLH_human_position.svg`
  - `scripts/bHLH_human_position_A4.r` -> `outputs/figures/bHLH_human_position_A4.svg`

All panels are colored by class using `data/raw/LS_classes.csv`.

### 3.2 Relative position of the bHLH domain in each human gene's isoforms

The isoform analysis evaluates whether the relative position of IPR011598 is conserved across coding transcripts for a given gene or shifts due to alternative splicing. In the thesis, each transcript is represented as a "pre / domain / post" bar scaled to [0, 1], with an additional dot encoding transcript length.

In this repository:

- table reconstruction (previously an unsaved bash step, now reproducible):
  - `scripts/make_bHLH_StartEnd_withISO.py` -> `data/intermediate/bHLH_StartEnd_withISO.csv`
- canonical isoform plots are the A3 exports:
  - `scripts/isoform_plots_A3.r` (class-specific panels)
  - `scripts/isoform_plot_merged_improved.R` (single merged plot)

### 3.3 Relative position of the bHLH domain across the phylogenetic tree

#### 3.3.2 Revised ortholog strategy: Ensembl Compara

The thesis moves from an initial ProteinOrtho attempt to an Ensembl Compara gene-tree reconciliation strategy (Release 114). This approach clusters homologous families, builds gene trees, and reconciles them to the species tree to label orthologs/paralogs.

#### 3.3.3 Species selection and dataset design

Species selection is inspired by Simionato et al. (2007) and adjusted to what is available in Ensembl Release 114. The working set (22 non-human species) is stored as:

- `data/raw/Species.txt`

and is used across ortholog and orthogroup visualizations (with a human reference axis/row added for plotting).

#### 3.3.4 Ortholog retrieval via Ensembl REST API

The thesis describes a custom Python retrieval pipeline against Ensembl REST `/homology/id/...` endpoints, switching between Compara databases (`vertebrates` vs `pan_homology`) depending on evolutionary distance and collecting homology type, IDs, sequences, and identity metrics.

In this repository, the orthology retrieval and downstream table building are captured in:

- `notebooks/Orthologs_retrival.ipynb` (REST calls + table assembly)
- `notebooks/Ortho_heatmap.ipynb` (table enrichment with display gene names)

The main integrated ortholog table used for plotting is:

- `data/intermediate/orthologs/annotated_bHLH_merged_data_with_gene_names.csv`

#### 3.3.6-3.3.8 Domain annotation via InterProScan and integration

The thesis uses InterProScan (v5.74-105.0) to annotate domains in ortholog proteins and then filters domain hits to the bHLH identifier (IPR011598). Finally, the merged table is enriched with target gene "display_name" via the Ensembl REST `/lookup/id/...` endpoint.

In this repository, the filtered/interpreted results are stored in the merged ortholog table above and used directly for visualization.

### 3.4 Phylogenetic data visualization

The thesis generates per-orthogroup plots where each ortholog is drawn on a 0-1 protein scale (grey baseline) and the bHLH domain is highlighted as a colored rectangle (geom_rect). The same idea is implemented here with two complementary plot sets:

- **Ensembl-only overview** (homogeneous species set):
  - `outputs/orthogroups/domain_positions_ensembl/orthogroup_<HGNC>.svg`
- **Mammals-only integrated zoom** (Ensembl + Zoonomia; heterogeneous sources):
  - `outputs/orthogroups/domain_positions_mammals_integrated/orthogroup_<HGNC>.svg`

Recommended entry point:

- `scripts/run_orthogroup_domain_plots.sh`

Species labels are formatted as Latin binomials (e.g., `Homo sapiens`) for readability.

### 3.5 Orthology refinement using the Zoonomia Project

The thesis motivates supplementing Ensembl orthology calls when key orthologs appear missing (e.g., due to assembly/annotation limitations). Zoonomia provides improved mammalian assemblies via Cactus alignment, and TOGA provides derived annotations for those assemblies.

In this repository:

- Zoonomia mapping notebook:
  - `notebooks/zoonomia_bHLH_mapping.ipynb` -> `data/intermediate/zoonomia/Zoonomia_Start_End_final.csv`
- conversion to relative coordinates for plotting:
  - `scripts/prepare_zoonomia_relpos.py` -> `data/intermediate/zoonomia/Zoonomia_Start_End_final_with_relpos.csv`

### Updated: parallel coordinates (ortholog midpoints across species)

- `scripts/Parallel_coordinates.r` -> `outputs/orthologs/figures/parallel_coordinates_midpoints.svg`

Duplicate hits for the same (gene, species) pair can occur (e.g., multiple Pfam hits on the same target protein). The default strategy is:

- `BHLH_PARCOORD_DUPLICATE_STRATEGY=select_best` (default): select the hit whose midpoint is closest to the human midpoint

Alternatives:

- `collapse_median`: collapse duplicates by median midpoint
- `drop_gene`: drop genes that have any non-unique (gene, species) pairs

A duplicate-hit report is written to:

- `outputs/orthologs/reports/parallel_coordinates_duplicate_hits.csv`

### Project extension: RBP / 3'UTR section

This section is not part of the core thesis Chapter 3 pipeline, but extends the project to explore post-transcriptional regulation via RNA-binding proteins (RBPs) and 3'UTR length.

- Data preparation:
  - `notebooks/3'UTR.ipynb` -> `data/intermediate/rbp/Violin_plot_3UTR.csv` and `data/intermediate/rbp/RBP_binary_matrix.csv`
- Figures (all color-coded by bHLH class; outputs under `outputs/rbp/figures/`):
  - `scripts/rbp/utr_length_boxplots.R`
  - `scripts/rbp/rbp_vs_utr_scatter_and_hist.R`
  - `scripts/rbp/rbp_tsne_pam_clusters.R`
  - `scripts/rbp/rbp_tsne_hclust_clusters.R`
  - `scripts/rbp/rbp_tsne_kmodes_clusters.R`
  - `scripts/rbp/rbp_dendrogram_linkage_methods.R`
