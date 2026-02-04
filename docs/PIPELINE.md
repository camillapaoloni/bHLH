# bHLH Project Pipeline (reconstructed)

This pipeline summarizes the steps used to generate the main outputs in the project, with inputs and outputs.
All paths are relative to the project root.
If you run scripts from a different directory, set the `BHLH_PROJECT_ROOT` environment variable to the project root.
Note: very large files `data/raw/Zoonomia_protaln/*.fa` are excluded from Git (GitHub 100MB limit); keep them locally or store with external storage/Git LFS.

## 1) BioMart metadata
- `scripts/BiomartData.r` → input `data/intermediate/bHLH_Ensembl_IDs.txt` (or `data/raw/Lambert_bHLH.csv` column `Ensembl ID`) → outputs in `data/intermediate/Metadata_CSVs/`:
  - `Transcript_Attributes.csv`
  - `Pfam_Domains.csv`
  - `InterPro_Domains.csv`
  - `SMART_Domains.csv`
  - `External_References.csv`
  - `Prosite_Domains.csv`
  - `CDD_Domains.csv`

## 2) Metadata CSV cleaning
- `notebooks/CSVs_parsing.ipynb` → input `data/raw/Human_TFs.csv`, `data/raw/Lambert_bHLH.csv`, `data/intermediate/Metadata_CSVs/*.csv` → outputs `data/intermediate/Metadata_CSVs/*_cleaned.csv` + CSV with HGNC column.

## 3) FASTA + bHLH domains
- `notebooks/RoadTo_FASTA.ipynb` → input `data/intermediate/Metadata_CSVs/InterPro_Domains_cleaned.csv`, `data/intermediate/Metadata_CSVs/Pfam_Domains_cleaned.csv` → outputs:
  - `data/intermediate/interpro/InterPro_Domains_bHLH_filtered.csv`
  - `data/intermediate/bHLH_transcripts_CDS.fasta` (Ensembl REST)
  - `data/intermediate/bHLH_transcripts_protein.fasta`
  - `data/intermediate/longest_isoform.fasta`
  - `data/intermediate/bHLH_domains.fasta`

## 4) Core bHLH tables
- `notebooks/bHLH_classes.ipynb` → input `data/intermediate/Metadata_CSVs/InterPro_Domains_cleaned.csv`, `data/intermediate/Metadata_CSVs/Transcript_Attributes.csv`, `data/intermediate/Metadata_CSVs/Pfam_Domains_cleaned.csv`, `data/raw/LS_classes.csv`, `data/intermediate/longest_isoform.fasta`, `data/raw/LI_HGNC.fasta` → outputs:
  - `data/intermediate/longest_isoform.csv`
  - `data/intermediate/table_input.csv`
  - `data/intermediate/table_input_withPAS.csv` (later renamed to `data/intermediate/table_input_PAS.csv`)

## 4b) bHLH_StartEnd_withISO.csv (former bash step)
- `scripts/make_bHLH_StartEnd_withISO.py` → input `data/intermediate/Metadata_CSVs/InterPro_Domains_cleaned.csv`, `data/intermediate/Metadata_CSVs/Transcript_Attributes_cleaned.csv` → output `data/intermediate/bHLH_StartEnd_withISO.csv`

## 5) Main domain figures
- `scripts/bHLH_human_position.r` → input `data/intermediate/table_input.csv`, `data/raw/LS_classes.csv` → output `outputs/figures/bHLH_human_position.svg`
- `scripts/bHLH_human_position_A4.r` → input `data/intermediate/table_input.csv`, `data/raw/LS_classes.csv` → output `outputs/figures/bHLH_human_position_A4.svg`
- `scripts/bHLH_human_PAS.r` → input `data/intermediate/table_input_PAS.csv` (fallback: `data/intermediate/table_input_withPAS.csv`) → outputs `outputs/figures/bHLH_human_position_PAS.svg`, `outputs/figures/bHLH_human_position_PAS.pdf`

## 6) Isoforms and domain position
- `scripts/isoforms_plots.R` → input `data/intermediate/bHLH_StartEnd_withISO.csv`, `data/raw/LS_classes.csv`, `data/raw/Lambert_bHLH.csv` → output `outputs/figures/bHLH_human_transcript_class*_A3.svg`
- `scripts/isoform_plots_A3.r` → same input → output `outputs/figures/bHLH_human_transcript_class*_A3.svg`
- `scripts/isoform_plot_merged.R` and `scripts/isoform_plot_merged_improved.R` → same input → output `outputs/figures/bHLH_singleplot_isoforms.svg`
- `scripts/isoform_old_plots.r` → same input → output `outputs/figures/plots/ISO_relposbHLH_*.svg`

## 7) RBP / 3'UTR
- `notebooks/3'UTR.ipynb` → input `data/raw/RBP_3UTR_data/3UTR_merged_clean.csv`, `data/raw/LS_classes.csv` → output `data/intermediate/rbp/Violin_plot_3UTR.csv`, `data/intermediate/rbp/RBP_binary_matrix.csv`
- `scripts/violin_length3utr.r` → input `data/intermediate/rbp/Violin_plot_3UTR.csv` → output `outputs/rbp/UTR_length_boxplots.svg`
- `scripts/RBPvs3UTR_scatter.r` → input `data/intermediate/rbp/RBP_binary_matrix.csv`, `data/intermediate/rbp/Violin_plot_3UTR.csv` → outputs `outputs/rbp/RBPvs3UTR_scatter.svg`, `outputs/rbp/RBPvsTFs_histogram.svg`
- `scripts/Hierarch_gower.r` → input `data/intermediate/rbp/RBP_binary_matrix.csv` → output `outputs/rbp/hierarchical_clusters.svg`
- `scripts/K-means_RBP.r` → input `data/intermediate/rbp/RBP_binary_matrix.csv` → output `outputs/rbp/k_means_clusters.svg`
- `scripts/K-modes.r` → input `data/intermediate/rbp/RBP_binary_matrix.csv`, `data/raw/LS_classes.csv` → output `outputs/rbp/k_modes_clusters.svg`
- `scripts/dendrograms_UTR.r` → input `data/intermediate/rbp/RBP_binary_matrix.csv`, `data/raw/LS_classes.csv` → output `outputs/rbp/dendrograms_linkage_methods.svg`

## 8) Ortho / TOGA / Zoonomia
- `notebooks/Orthologs_retrival.ipynb` → input `data/intermediate/Metadata_CSVs/Transcript_Attributes_cleaned.csv`, `data/intermediate/Metadata_CSVs/InterPro_Domains_cleaned.csv`, `data/intermediate/table_input.csv`, `data/intermediate/interpro/interpro_output.tsv` → outputs:
  - `data/intermediate/orthologs/homology_results.csv`
  - `data/intermediate/interpro/OutputInterProScan.tsv`
  - `data/intermediate/interpro/IPR011598_IPRScan_final.tsv`
  - `data/intermediate/orthologs/in&out_bHLH_data.csv`
  - `data/intermediate/orthologs/annotated_bHLH_merged_data.csv`
- `notebooks/Ortho_heatmap.ipynb` → input `data/intermediate/orthologs/annotated_bHLH_merged_data.csv`, `data/raw/TOGA_orthologs/*.tsv` → outputs:
  - `data/intermediate/orthologs/annotated_bHLH_merged_data_with_gene_names.csv` (Ensembl REST)
  - `data/intermediate/orthologs/TOGA_orthologs_allSpecies.csv`
- `scripts/ortho_bHLH.R` → input `data/intermediate/orthologs/annotated_bHLH_merged_data_with_gene_names.csv` → outputs `outputs/orthogroups/Orthogroup_*.png` and `outputs/orthogroups/species_group_*.png`
- `notebooks/Random_Orthoplots.ipynb` → input `data/intermediate/orthologs/annotated_bHLH_merged_data.csv` → output `outputs/figures/a3_heatmap.svg`
- `scripts/Parallel_coordinates.r` → input `data/intermediate/orthologs/annotated_bHLH_merged_data_with_gene_names.csv` → output interactive plot (not saved)
- `notebooks/zoonomia_bHLH_mapping.ipynb` → input `data/raw/Zoonomia_protaln/*.fa`, `data/intermediate/bHLH_StartEnd_withISO.csv`, `data/intermediate/orthologs/TOGA_orthologs_allSpecies.csv` → outputs `data/intermediate/zoonomia/Zoonomia_Start_End_final.csv`, `outputs/reports/report_comparison.csv`

## 9) Phylogenetic tree (not used in the final project)
- `scripts/Tree.R` → input `data/raw/Tree23sp.newick` → output `outputs/figures/phylo23.svg` (and/or `outputs/figures/phylo23.png`)

## Missing steps (bash not saved)
- `data/intermediate/bHLH_StartEnd_withISO.csv` is now reproducible via `scripts/make_bHLH_StartEnd_withISO.py`.
