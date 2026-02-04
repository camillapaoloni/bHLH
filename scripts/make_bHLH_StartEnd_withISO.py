#!/usr/bin/env python3
"""
Build data/intermediate/bHLH_StartEnd_withISO.csv from cleaned InterPro + Transcript_Attributes.

Inputs (relative to project root):
- data/intermediate/Metadata_CSVs/InterPro_Domains_cleaned.csv
- data/intermediate/Metadata_CSVs/Transcript_Attributes_cleaned.csv

Output:
- data/intermediate/bHLH_StartEnd_withISO.csv
"""

import os
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(os.getenv("BHLH_PROJECT_ROOT", ".")).resolve()

def p(*parts: str) -> Path:
    return PROJECT_ROOT.joinpath(*parts)

INTERPRO_PATH = p('data', 'intermediate', 'Metadata_CSVs', 'InterPro_Domains_cleaned.csv')
TRANSCRIPTS_PATH = p('data', 'intermediate', 'Metadata_CSVs', 'Transcript_Attributes_cleaned.csv')
OUTPUT_PATH = p('data', 'intermediate', 'bHLH_StartEnd_withISO.csv')

if not INTERPRO_PATH.exists():
    raise SystemExit(f"Missing {INTERPRO_PATH}")
if not TRANSCRIPTS_PATH.exists():
    raise SystemExit(f"Missing {TRANSCRIPTS_PATH}")
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

# Load
interpro = pd.read_csv(INTERPRO_PATH)
transcripts = pd.read_csv(TRANSCRIPTS_PATH)

# Normalize column names for HGNC symbol
if 'HGNC symbol' in interpro.columns:
    hgnc_col = 'HGNC symbol'
elif 'HGNC.symbol' in interpro.columns:
    hgnc_col = 'HGNC.symbol'
else:
    raise SystemExit('HGNC symbol column not found in InterPro_Domains_cleaned.csv')

# Filter bHLH InterPro domain
bhlh = interpro[interpro['interpro'] == 'IPR011598'].copy()

# Keep only needed columns
bhlh = bhlh[[hgnc_col, 'ensembl_transcript_id', 'interpro_start', 'interpro_end']].copy()

# Ensure numeric
bhlh['interpro_start'] = pd.to_numeric(bhlh['interpro_start'], errors='coerce')
bhlh['interpro_end'] = pd.to_numeric(bhlh['interpro_end'], errors='coerce')

# Aggregate per transcript (min start, max end)
bhlh_agg = (
    bhlh
    .dropna(subset=['ensembl_transcript_id'])
    .groupby(['ensembl_transcript_id', hgnc_col], as_index=False)
    .agg(interpro_start=('interpro_start', 'min'),
         interpro_end=('interpro_end', 'max'))
)

# Protein length from CDS length (floor division)
if 'cds_length' not in transcripts.columns:
    raise SystemExit('cds_length column not found in Transcript_Attributes_cleaned.csv')

transcripts = transcripts[['ensembl_transcript_id', 'cds_length']].copy()
transcripts['protein_length'] = (transcripts['cds_length'] // 3).astype('Int64')

# Merge
out = bhlh_agg.merge(transcripts[['ensembl_transcript_id', 'protein_length']],
                     on='ensembl_transcript_id', how='left')

# Rename HGNC column to match existing CSV
out = out.rename(columns={hgnc_col: 'HGNC symbol'})

# Order columns
out = out[['HGNC symbol', 'ensembl_transcript_id', 'protein_length', 'interpro_start', 'interpro_end']]

# Save
out.to_csv(OUTPUT_PATH, index=False)
print(f"Saved: {OUTPUT_PATH} ({len(out)} rows)")
