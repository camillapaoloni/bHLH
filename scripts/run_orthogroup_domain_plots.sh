#!/usr/bin/env bash
set -euo pipefail

# Run only the steps needed to generate orthogroup domain-position plots.
# This script is designed to be re-runnable: it skips outputs that already exist.
#
# It does NOT call Ensembl/BioMart or other network-heavy steps.

PROJECT_ROOT="${BHLH_PROJECT_ROOT:-.}"

ENSEMBL_TABLE="${PROJECT_ROOT}/data/intermediate/orthologs/annotated_bHLH_merged_data_with_gene_names.csv"
ZOO_BASE="${PROJECT_ROOT}/data/intermediate/zoonomia/Zoonomia_Start_End_final.csv"
ZOO_RELPOS="${PROJECT_ROOT}/data/intermediate/zoonomia/Zoonomia_Start_End_final_with_relpos.csv"

if [[ ! -f "${ENSEMBL_TABLE}" ]]; then
  echo "Missing required input: ${ENSEMBL_TABLE}" >&2
  exit 1
fi

if [[ -f "${ZOO_BASE}" && ! -f "${ZOO_RELPOS}" ]]; then
  echo "Preparing Zoonomia relative positions: ${ZOO_RELPOS}"
  python "${PROJECT_ROOT}/scripts/prepare_zoonomia_relpos.py" --project-root "${PROJECT_ROOT}"
fi

echo "Generating orthogroup plots (skipping existing outputs)..."
BHLH_SKIP_EXISTING=1 Rscript "${PROJECT_ROOT}/scripts/ortho_bHLH.R"

echo "Done. Outputs:"
echo "  - ${PROJECT_ROOT}/outputs/orthogroups/domain_positions_ensembl/"
echo "  - ${PROJECT_ROOT}/outputs/orthogroups/domain_positions_mammals_integrated/"
