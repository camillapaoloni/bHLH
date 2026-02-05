#!/usr/bin/env python3
"""
Prepare Zoonomia bHLH domain coordinates with protein lengths and relative positions.

Why this exists
---------------
`data/intermediate/zoonomia/Zoonomia_Start_End_final.csv` contains mapped absolute
domain coordinates (query_start/query_end) but does not include the ortholog
protein length needed to compute relative positions (0..1).

We recover protein lengths from the Zoonomia filtered alignments in:
  data/raw/Zoonomia_protaln/filtered_alignments/*.fa

Each record header is expected to look like:
  >ENST....GENE.... | PROT | REFERENCE
  >ENST....GENE.... | PROT | QUERY

We compute:
  - human_length: number of non-gap characters in REFERENCE sequence
  - query_length: number of non-gap characters in QUERY sequence
  - rel_*: absolute position divided by corresponding length

Output
------
Writes:
  data/intermediate/zoonomia/Zoonomia_Start_End_final_with_relpos.csv

The output is intended to be consumed by downstream plotting scripts.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re

import pandas as pd


HEADER_RE = re.compile(r"^>?(.+?)\s*\|\s*PROT\s*\|\s*(REFERENCE|QUERY)\s*$", re.IGNORECASE)


SPECIES_MAP = {
    "Mus_protaln": "mus_musculus",
    "Opossum_protaln": "monodelphis_domestica",
    "Macaca_protaln": "macaca_mulatta",
    "Bos_protaln": "bos_taurus",
    "Rat_protaln": "rattus_norvegicus",
    "Canis_protaln": "canis_lupus_familiaris",
    "Chimp_protaln": "pan_troglodytes",
    "Gorilla_protaln": "gorilla_gorilla",
}


@dataclass(frozen=True)
class LengthKey:
    enst: str
    target_species: str


def ungapped_len(seq: str) -> int:
    # Treat '-' as alignment gap; keep everything else (including X).
    return sum(1 for ch in seq if ch != "-")


def parse_fasta_lengths(filtered_alignments_dir: Path) -> dict[LengthKey, dict[str, int]]:
    """
    Returns a mapping:
      (ENST, target_species) -> {'human_length': int, 'query_length': int}
    """
    lengths: dict[LengthKey, dict[str, int]] = {}

    current_header: str | None = None
    current_seq_parts: list[str] = []

    def flush_record(file_species: str) -> None:
        nonlocal current_header, current_seq_parts
        if current_header is None:
            return
        m = HEADER_RE.match(current_header.strip())
        if not m:
            current_header = None
            current_seq_parts = []
            return

        full_enst_label = m.group(1).strip()
        label = m.group(2).upper()
        enst = full_enst_label.split()[0].split(".")[0]
        seq = "".join(current_seq_parts).strip()

        key = LengthKey(enst=enst, target_species=file_species)
        lengths.setdefault(key, {})
        if label == "REFERENCE":
            lengths[key]["human_length"] = ungapped_len(seq)
        elif label == "QUERY":
            lengths[key]["query_length"] = ungapped_len(seq)

        current_header = None
        current_seq_parts = []

    for fp in sorted(filtered_alignments_dir.glob("*.fa")):
        raw_species = fp.stem
        species = SPECIES_MAP.get(raw_species, raw_species)

        current_header = None
        current_seq_parts = []

        with fp.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    flush_record(species)
                    current_header = line[1:]
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line.strip())

        flush_record(species)

    return lengths


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--project-root",
        default=".",
        help="Project root (or set BHLH_PROJECT_ROOT env var).",
    )
    ap.add_argument(
        "--input",
        default="data/intermediate/zoonomia/Zoonomia_Start_End_final.csv",
        help="Input Zoonomia mapping CSV.",
    )
    ap.add_argument(
        "--filtered-alignments",
        default="data/raw/Zoonomia_protaln/filtered_alignments",
        help="Directory containing filtered alignment FASTA files.",
    )
    ap.add_argument(
        "--output",
        default="data/intermediate/zoonomia/Zoonomia_Start_End_final_with_relpos.csv",
        help="Output CSV path.",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output if it already exists.",
    )
    args = ap.parse_args()

    project_root = Path(args.project_root).resolve()
    in_path = project_root / args.input
    fa_dir = project_root / args.filtered_alignments
    out_path = project_root / args.output

    if out_path.exists() and not args.overwrite:
        print(f"Skip (exists): {out_path}")
        return 0

    if not in_path.exists():
        raise FileNotFoundError(in_path)
    if not fa_dir.exists():
        raise FileNotFoundError(fa_dir)

    df = pd.read_csv(in_path)
    required = {"HGNC", "ENST", "target_species", "human_start", "human_end", "query_start", "query_end"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {in_path}: {sorted(missing)}")

    lengths = parse_fasta_lengths(fa_dir)

    def lookup_len(row: pd.Series) -> pd.Series:
        key = LengthKey(enst=str(row["ENST"]).split(".")[0], target_species=str(row["target_species"]))
        info = lengths.get(key, {})
        return pd.Series(
            {
                "human_length": info.get("human_length", pd.NA),
                "query_length": info.get("query_length", pd.NA),
            }
        )

    df[["human_length", "query_length"]] = df.apply(lookup_len, axis=1)

    # Compute relative positions when lengths are available.
    for col in ["human_start", "human_end", "query_start", "query_end", "human_length", "query_length"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["rel_human_start"] = df["human_start"] / df["human_length"]
    df["rel_human_end"] = df["human_end"] / df["human_length"]
    df["rel_query_start"] = df["query_start"] / df["query_length"]
    df["rel_query_end"] = df["query_end"] / df["query_length"]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"Wrote: {out_path} ({df.shape[0]} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

