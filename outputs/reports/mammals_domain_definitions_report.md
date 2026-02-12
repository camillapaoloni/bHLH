# Mammals domain definitions report

- Generated: 2026-02-11 17:17:20
- Ensembl input: `./data/intermediate/orthologs/annotated_bHLH_merged_data_with_gene_names.csv`
- Zoonomia input: `./data/intermediate/zoonomia/Zoonomia_Start_End_final_with_relpos.csv`

This report contrasts two domain definitions for the mammals subset:

- **Predicted**: Ensembl/InterProScan-derived domain coordinates on the target protein.
- **Projected**: human-domain coordinates projected onto the query protein via Zoonomia protein alignments.

## Selection logic (why some entries do not appear)

- Plots are restricted to the 8 mammals available in the Zoonomia alignment set.
- For each (HGNC, species), a **single representative** predicted entry is selected (preference: `one2one`, then longest protein).
- For each (HGNC, species), a **single representative** projected entry is selected (longest `query_length`).
- Rows with invalid projected/predicted coordinates are dropped (`NA/Inf`, outside `[0,1]`, or `rel_end <= rel_start`).

## Outputs

- Plots: `./outputs/orthogroups/domain_definitions_mammals/orthogroup_<HGNC>.svg`
- Summary table: `./outputs/reports/mammals_domain_definitions_summary.csv`

## Status summary (cells)

- Type match: 695
- Type mismatch: 33
- Ensembl only: 99
- TOGA only: 49

## Orthogroups plotted by default

- Flagged (any non-Type match cell): 83
- Total in table: 111

Set `BHLH_DOMAINPLOT_ALL=1` to plot all orthogroups.
