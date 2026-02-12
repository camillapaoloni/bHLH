# Scripts overview

This folder contains the project’s main executable scripts (R + Python + shell).

## Core orthogroup plots

- `run_orthogroup_domain_plots.sh` — convenience runner for the plotting pipeline.
- `ortho_bHLH.R` — generates per-orthogroup plots:
  - Ensembl-only: `outputs/orthogroups/domain_positions_ensembl/`
  - Mammals integrated (Ensembl + Zoonomia): `outputs/orthogroups/domain_positions_mammals_integrated/`

## Mammals: two domain definitions (QC)

- `plot_mammals_domain_definitions.R` — contrasts:
  - **Predicted** (Ensembl/InterProScan on target proteins)
  - **Projected** (Zoonomia human→query projection)
  - outputs: `outputs/orthogroups/domain_definitions_mammals/` + `outputs/reports/mammals_domain_definitions_*`

## Orthology concordance report

- `report_ensembl_vs_toga_mammals.R` — Ensembl vs TOGA/TOGA summary tables and matrix figure:
  - `outputs/reports/ensembl_vs_toga_mammals_*`
  - `outputs/figures/ensembl_vs_toga_mammals_matrix.svg`

## Zoonomia projection preprocessing

- `prepare_zoonomia_relpos.py` — adds protein lengths and relative coordinates to the Zoonomia mapping table:
  - input: `data/intermediate/zoonomia/Zoonomia_Start_End_final.csv`
  - output: `data/intermediate/zoonomia/Zoonomia_Start_End_final_with_relpos.csv`

## Domain conservation follow-up (local analysis)

- `domain_conservation_analysis.py` — extracts domain sequences (predicted vs projected) and computes identity metrics.
  - outputs are local (git-ignored): `outputs/analysis/domain_conservation_analysis/` (see `README.md`)

## Human / isoform plots

Scripts prefixed with `bHLH_human_*` and `isoform_*` are used for human-focused figures (see `README.md`).
