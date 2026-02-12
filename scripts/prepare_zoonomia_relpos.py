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
import csv
from dataclasses import dataclass
from pathlib import Path
import re


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


def _to_float(value: str) -> float | None:
    try:
        if value is None:
            return None
        s = str(value).strip()
        if s == "":
            return None
        return float(s)
    except Exception:
        return None


def _fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if value != value:  # NaN
            return ""
    return str(value)


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

    with in_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"Empty or invalid CSV: {in_path}")
        fieldnames = list(reader.fieldnames)

        required = {"HGNC", "ENST", "target_species", "human_start", "human_end", "query_start", "query_end"}
        missing = required - set(fieldnames)
        if missing:
            raise ValueError(f"Missing required columns in {in_path}: {sorted(missing)}")

        rows = list(reader)

    lengths = parse_fasta_lengths(fa_dir)

    extra_cols = ["human_length", "query_length", "rel_human_start", "rel_human_end", "rel_query_start", "rel_query_end"]
    out_fieldnames = fieldnames[:]
    for c in extra_cols:
        if c not in out_fieldnames:
            out_fieldnames.append(c)

    out_rows: list[dict[str, str]] = []
    for row in rows:
        enst = str(row.get("ENST", "")).split(".")[0]
        species = str(row.get("target_species", ""))
        info = lengths.get(LengthKey(enst=enst, target_species=species), {})
        human_length = info.get("human_length")
        query_length = info.get("query_length")

        human_start = _to_float(row.get("human_start", ""))
        human_end = _to_float(row.get("human_end", ""))
        query_start = _to_float(row.get("query_start", ""))
        query_end = _to_float(row.get("query_end", ""))

        rel_human_start = (human_start / human_length) if (human_start is not None and human_length) else None
        rel_human_end = (human_end / human_length) if (human_end is not None and human_length) else None
        rel_query_start = (query_start / query_length) if (query_start is not None and query_length) else None
        rel_query_end = (query_end / query_length) if (query_end is not None and query_length) else None

        out_row: dict[str, str] = {k: _fmt(row.get(k, "")) for k in out_fieldnames}
        out_row["human_length"] = _fmt(human_length)
        out_row["query_length"] = _fmt(query_length)
        out_row["rel_human_start"] = _fmt(rel_human_start)
        out_row["rel_human_end"] = _fmt(rel_human_end)
        out_row["rel_query_start"] = _fmt(rel_query_start)
        out_row["rel_query_end"] = _fmt(rel_query_end)
        out_rows.append(out_row)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"Wrote: {out_path} ({len(out_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
