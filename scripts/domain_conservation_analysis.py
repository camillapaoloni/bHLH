#!/usr/bin/env python3
"""
Domain conservation follow-up for mammals (Ensembl vs Zoonomia/TOGA).

This script is designed to support meeting-ready biological validation for
orthogroups flagged as Ensembl-vs-TOGA "Type mismatch" in:
  outputs/reports/mammals_domain_definitions_summary.csv

It extracts bHLH domain sequences from:
  - Ensembl/InterProScan-derived coordinates on target proteins (predicted)
  - Zoonomia projected coordinates on query proteins (projected)

and computes simple pairwise identity metrics via a lightweight global alignment.

Outputs (default):
  outputs/analysis/domain_conservation_analysis/
    overall_summary.csv
    mismatch_cells.csv
    mismatch_cells.md
    meeting_shortlist.csv
    meeting_shortlist.md
    genes/
      <HGNC>/
        mismatch_rows.csv
        predicted_domains.faa
        projected_domains.faa
        all_domains.faa
        predicted_metrics.csv
        projected_metrics.csv
        cross_metrics.csv
        report.md

Notes
-----
- Coordinates are treated as 0-based, end-exclusive (as used throughout the project:
  domain_length ~= end-start and rel_start = start/length).
- The script purposefully keeps projected sequences for all ENST entries present for
  a gene/species (useful for one2many/many2many cases).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional

import pandas as pd


MAMMALS = [
    "pan_troglodytes",
    "gorilla_gorilla",
    "macaca_mulatta",
    "bos_taurus",
    "canis_lupus_familiaris",
    "rattus_norvegicus",
    "mus_musculus",
    "monodelphis_domestica",
]


SPECIES_FILE = {
    "mus_musculus": "Mus_protaln.fa",
    "monodelphis_domestica": "Opossum_protaln.fa",
    "macaca_mulatta": "Macaca_protaln.fa",
    "bos_taurus": "Bos_protaln.fa",
    "rattus_norvegicus": "Rat_protaln.fa",
    "canis_lupus_familiaris": "Canis_protaln.fa",
    "pan_troglodytes": "Chimp_protaln.fa",
    "gorilla_gorilla": "Gorilla_protaln.fa",
}


HEADER_RE = re.compile(r"^>?(.+?)\s*\|\s*PROT\s*\|\s*(REFERENCE|QUERY)\s*$", re.IGNORECASE)


def _is_finite_number(x: object) -> bool:
    try:
        return x is not None and not (isinstance(x, float) and math.isnan(x)) and math.isfinite(float(x))
    except Exception:
        return False


def _as_int(x: object) -> Optional[int]:
    try:
        if x is None:
            return None
        if isinstance(x, float) and math.isnan(x):
            return None
        return int(float(x))
    except Exception:
        return None


def _slice_0based(seq: str, start: int, end: int) -> str:
    if start < 0:
        start = 0
    if end < 0:
        end = 0
    if end < start:
        end = start
    return seq[start:end]


def ungap(seq: str) -> str:
    return "".join(ch for ch in str(seq).strip() if ch != "-")


def read_fasta_records(fp: Path) -> Iterator[tuple[str, str]]:
    header: Optional[str] = None
    parts: list[str] = []
    with fp.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


@dataclass(frozen=True)
class ProtAlnKey:
    species: str
    full_enst_label: str


@dataclass
class ProtAlnRecord:
    reference_seq: Optional[str] = None
    query_seq: Optional[str] = None


def load_alignment_sequences(
    alignments_dir: Path,
    needed: set[ProtAlnKey],
) -> dict[ProtAlnKey, ProtAlnRecord]:
    """
    Loads only the requested (species, full_enst_label) records from the filtered
    Zoonomia alignment FASTA files.
    """
    by_species: dict[str, set[str]] = {}
    for k in needed:
        by_species.setdefault(k.species, set()).add(k.full_enst_label)

    out: dict[ProtAlnKey, ProtAlnRecord] = {}
    for species, labels in by_species.items():
        fn = SPECIES_FILE.get(species)
        if not fn:
            continue
        fp = alignments_dir / fn
        if not fp.exists():
            continue

        for hdr, seq in read_fasta_records(fp):
            m = HEADER_RE.match(hdr)
            if not m:
                continue
            full_label = m.group(1).strip()
            if full_label not in labels:
                continue
            role = m.group(2).upper()
            key = ProtAlnKey(species=species, full_enst_label=full_label)
            rec = out.setdefault(key, ProtAlnRecord())
            seq_u = ungap(seq)
            if role == "REFERENCE":
                rec.reference_seq = seq_u
            elif role == "QUERY":
                rec.query_seq = seq_u

    return out


def nw_align(a: str, b: str, match: int = 1, mismatch: int = 0, gap: int = -1) -> tuple[str, str]:
    """
    Simple Needleman–Wunsch global alignment for short sequences.
    Returns aligned (a_aln, b_aln).
    """
    a = a or ""
    b = b or ""
    n = len(a)
    m = len(b)

    # score matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    bt = [[0] * (m + 1) for _ in range(n + 1)]  # 0 diag, 1 up, 2 left

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap
        bt[i][0] = 1
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap
        bt[0][j] = 2

    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            bj = b[j - 1]
            s_diag = dp[i - 1][j - 1] + (match if ai == bj else mismatch)
            s_up = dp[i - 1][j] + gap
            s_left = dp[i][j - 1] + gap
            best = s_diag
            direction = 0
            if s_up > best:
                best = s_up
                direction = 1
            if s_left > best:
                best = s_left
                direction = 2
            dp[i][j] = best
            bt[i][j] = direction

    # traceback
    i, j = n, m
    a_aln: list[str] = []
    b_aln: list[str] = []
    while i > 0 or j > 0:
        direction = bt[i][j] if i >= 0 and j >= 0 else 0
        if i > 0 and j > 0 and direction == 0:
            a_aln.append(a[i - 1])
            b_aln.append(b[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or direction == 1):
            a_aln.append(a[i - 1])
            b_aln.append("-")
            i -= 1
        else:
            a_aln.append("-")
            b_aln.append(b[j - 1])
            j -= 1

    return "".join(reversed(a_aln)), "".join(reversed(b_aln))


def identity_metrics(ref: str, qry: str) -> dict[str, object]:
    a_aln, b_aln = nw_align(ref, qry)
    matches = 0
    aligned = 0
    nongap = 0
    gaps = 0
    for x, y in zip(a_aln, b_aln):
        aligned += 1
        if x == "-" or y == "-":
            gaps += 1
        else:
            nongap += 1
            if x == y:
                matches += 1
    ident = (matches / nongap) if nongap else float("nan")
    return {
        "ref_len": len(ref),
        "qry_len": len(qry),
        "aligned_len": aligned,
        "nongap_len": nongap,
        "matches": matches,
        "identity": ident,
        "gap_fraction": (gaps / aligned) if aligned else float("nan"),
    }


def write_fasta(fp: Path, records: Iterable[tuple[str, str]]) -> None:
    fp.parent.mkdir(parents=True, exist_ok=True)
    with fp.open("w", encoding="utf-8") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            s = seq.strip()
            for i in range(0, len(s), 80):
                f.write(s[i : i + 80] + "\n")


def md_species(species: str) -> str:
    return species.replace("_", " ")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=".")
    ap.add_argument(
        "--summary-csv",
        default="outputs/reports/mammals_domain_definitions_summary.csv",
    )
    ap.add_argument(
        "--ensembl-csv",
        default="data/intermediate/orthologs/annotated_bHLH_merged_data_with_gene_names.csv",
    )
    ap.add_argument(
        "--zoonomia-csv",
        default="data/intermediate/zoonomia/Zoonomia_Start_End_final_with_relpos.csv",
    )
    ap.add_argument(
        "--zoonomia-map",
        default="data/intermediate/zoonomia/Zoonomia_filtered_orthologs.csv",
        help="Optional mapping file containing q_gene (reg_*) for labels.",
    )
    ap.add_argument(
        "--alignments-dir",
        default="data/raw/Zoonomia_protaln/filtered_alignments",
    )
    ap.add_argument(
        "--outdir",
        default="outputs/analysis/domain_conservation_analysis",
    )
    ap.add_argument("--mode", choices=["all", "mismatch"], default="all", help="Which genes to process.")
    ap.add_argument("--only", default="", help="Comma-separated HGNC symbols.")
    ap.add_argument(
        "--clean",
        action="store_true",
        help="Delete existing outputs under outdir before writing new ones (including legacy outdir/all and outdir/mismatch folders).",
    )
    ap.add_argument(
        "--include-all-mismatches",
        action="store_true",
        help="Meeting shortlist: include all HGNC with any Type mismatch cell.",
    )
    ap.add_argument(
        "--include-pan-mismatches",
        action="store_true",
        default=True,
        help="Meeting shortlist: include HGNC with mismatch in pan_troglodytes (default: true).",
    )
    ap.add_argument(
        "--include-max-rat",
        action="store_true",
        default=True,
        help="Meeting shortlist: always include MAX mismatch in rattus_norvegicus (default: true).",
    )
    ap.add_argument(
        "--top-delta",
        type=int,
        default=5,
        help="Meeting shortlist: also include the top N HGNC by max(delta_start+delta_end) among mismatches.",
    )
    args = ap.parse_args()

    root = Path(args.project_root).resolve()
    p = lambda *xs: (root / Path(*xs)).resolve()

    summary_path = p(args.summary_csv)
    ens_path = p(args.ensembl_csv)
    zoo_path = p(args.zoonomia_csv)
    zoo_map_path = p(args.zoonomia_map)
    aln_dir = p(args.alignments_dir)
    outdir = p(args.outdir)

    if not summary_path.exists():
        raise FileNotFoundError(summary_path)
    if not ens_path.exists():
        raise FileNotFoundError(ens_path)
    if not zoo_path.exists():
        raise FileNotFoundError(zoo_path)

    df_sum = pd.read_csv(summary_path)
    df_mis = df_sum[df_sum["cell_status"].eq("Type mismatch")].copy()

    def mismatch_priority_list() -> list[str]:
        if df_mis.empty:
            return []
        if args.include_all_mismatches:
            return sorted({x for x in df_mis["hgnc"].dropna().astype(str).unique() if x.strip()})
        selected_m: set[str] = set()
        if args.include_pan_mismatches:
            selected_m.update(sorted(df_mis[df_mis["target_species"].eq("pan_troglodytes")]["hgnc"].unique()))

        if args.include_max_rat:
            has_max_rat = not df_mis[
                (df_mis["hgnc"].eq("MAX")) & (df_mis["target_species"].eq("rattus_norvegicus"))
            ].empty
            if has_max_rat:
                selected_m.add("MAX")

        if args.top_delta and args.top_delta > 0:
            df_tmp = df_mis.copy()
            df_tmp["delta_total"] = (
                df_tmp["delta_start"].fillna(0).astype(float) + df_tmp["delta_end"].fillna(0).astype(float)
            )
            top = (
                df_tmp.groupby("hgnc")["delta_total"]
                .max()
                .sort_values(ascending=False)
                .head(args.top_delta)
                .index.tolist()
            )
            selected_m.update(top)

        # Also include genes with mismatch in >=2 species.
        multi = df_mis.groupby("hgnc").size()
        selected_m.update(multi[multi >= 2].index.tolist())
        return sorted({x for x in selected_m if isinstance(x, str) and x.strip()})

    mismatch_selected = mismatch_priority_list()

    selected: set[str] = set()
    only = [x.strip() for x in args.only.split(",") if x.strip()]
    if only:
        selected.update(only)
    elif args.mode == "mismatch":
        selected.update(sorted(df_mis["hgnc"].dropna().astype(str).unique()))
    else:
        selected.update(sorted(df_sum["hgnc"].dropna().astype(str).unique()))

    selected = {x for x in selected if isinstance(x, str) and x.strip()}
    selected_list = sorted(selected)

    outdir.mkdir(parents=True, exist_ok=True)
    genes_dir = outdir / "genes"
    if args.clean:
        import shutil

        shutil.rmtree(genes_dir, ignore_errors=True)
        shutil.rmtree(outdir / "all", ignore_errors=True)
        shutil.rmtree(outdir / "mismatch", ignore_errors=True)
        for fp in [
            outdir / "overall_summary.csv",
            outdir / "mismatch_cells.csv",
            outdir / "mismatch_cells.md",
            outdir / "meeting_shortlist.csv",
            outdir / "meeting_shortlist.md",
        ]:
            try:
                fp.unlink()
            except FileNotFoundError:
                pass
    genes_dir.mkdir(parents=True, exist_ok=True)

    # ---- Mismatch cells (full; restricted to processed gene set) ----
    mismatch_cells = (
        df_mis[df_mis["hgnc"].isin(selected_list)]
        .sort_values(["hgnc", "target_species"])
        .reset_index(drop=True)
    )
    mismatch_cells_csv = outdir / "mismatch_cells.csv"
    mismatch_cells.to_csv(mismatch_cells_csv, index=False)

    mismatch_cells_md = outdir / "mismatch_cells.md"
    with mismatch_cells_md.open("w", encoding="utf-8") as f:
        f.write("# Mammals domain mismatch cells\n\n")
        f.write(f"- Source: `{summary_path.relative_to(root)}`\n")
        f.write(f"- Genes processed: {len(selected_list)}\n")
        f.write(f"- Mismatch rows in processed set: {len(mismatch_cells)}\n\n")
        if mismatch_cells.empty:
            f.write("No mismatches for the processed gene set.\n")
        else:
            for _, r in mismatch_cells.iterrows():
                f.write(
                    f"- **{r['hgnc']}** | {r['target_species']} — Ensembl: {r['ensembl_types_simple']} vs TOGA: {r['toga_types']} "
                    f"(Δstart={float(r['delta_start']):.3f}, Δend={float(r['delta_end']):.3f})\n"
                )

    # ---- Meeting shortlist (prioritized subset; restricted to processed gene set) ----
    meeting_selected = [h for h in mismatch_selected if h in set(selected_list)]
    meeting_rows = (
        df_mis[df_mis["hgnc"].isin(meeting_selected)]
        .sort_values(["hgnc", "target_species"])
        .reset_index(drop=True)
    )
    meeting_csv = outdir / "meeting_shortlist.csv"
    meeting_rows.to_csv(meeting_csv, index=False)

    meeting_md = outdir / "meeting_shortlist.md"
    with meeting_md.open("w", encoding="utf-8") as f:
        f.write("# Mammals domain conservation — meeting shortlist\n\n")
        f.write(f"- Source: `{summary_path.relative_to(root)}`\n")
        f.write(f"- Genes processed: {len(selected_list)}\n")
        f.write(f"- Genes shortlisted: {len(meeting_selected)}\n\n")
        f.write("Heuristics:\n\n")
        f.write("- mismatch in *Pan troglodytes*\n")
        f.write("- mismatch in >=2 species\n")
        f.write("- top coordinate deltas (delta_start+delta_end)\n")
        f.write("- always include MAX mismatch in rat (if present)\n\n")
        for h in meeting_selected:
            pan = not df_mis[(df_mis["hgnc"].eq(h)) & (df_mis["target_species"].eq("pan_troglodytes"))].empty
            nsp = int(df_mis[df_mis["hgnc"].eq(h)]["target_species"].nunique())
            df_h = df_mis[df_mis["hgnc"].eq(h)]
            dmax = float(
                (df_h["delta_start"].fillna(0).astype(float) + df_h["delta_end"].fillna(0).astype(float)).max()
            )
            flags = []
            if pan:
                flags.append("pan")
            if nsp >= 2:
                flags.append(f"{nsp} spp")
            flags.append(f"maxΔ={dmax:.3f}")
            f.write(f"- **{h}** ({', '.join(flags)})\n")

    # ---- Load inputs ----
    df_ens = pd.read_csv(ens_path)
    df_zoo = pd.read_csv(zoo_path)

    df_qgene = None
    if zoo_map_path.exists():
        try:
            df_qgene = pd.read_csv(zoo_map_path)
        except Exception:
            df_qgene = None

    # Subset to relevant genes/species early.
    df_ens = df_ens[df_ens["HGNC symbol"].isin(selected_list)].copy()
    df_zoo = df_zoo[df_zoo["HGNC"].isin(selected_list)].copy()
    df_zoo = df_zoo[df_zoo["target_species"].isin(MAMMALS)].copy()

    # Needed alignment keys
    needed_keys: set[ProtAlnKey] = set()
    for _, r in df_zoo.iterrows():
        sp = str(r.get("target_species", "")).strip()
        lbl = str(r.get("full_enst_label", "")).strip()
        if sp and lbl and sp in SPECIES_FILE:
            needed_keys.add(ProtAlnKey(species=sp, full_enst_label=lbl))

    aln = load_alignment_sequences(aln_dir, needed_keys) if needed_keys else {}

    overall_rows: list[dict[str, object]] = []

    for h in selected_list:
        h_dir = genes_dir / h
        h_dir.mkdir(parents=True, exist_ok=True)

        # Mismatch rows for this HGNC
        mm = df_mis[df_mis["hgnc"].eq(h)].copy().sort_values("target_species")
        mm.to_csv(h_dir / "mismatch_rows.csv", index=False)

        # --- Human predicted reference domain (canonical) ---
        ens_h = df_ens[df_ens["HGNC symbol"].eq(h)].copy()
        human_pred = None
        human_pred_meta = {}
        if not ens_h.empty:
            first = ens_h.iloc[0]
            start_q = _as_int(first.get("Start_Q"))
            stop_q = _as_int(first.get("Stop_Q"))
            src_seq = str(first.get("source_seq") or "")
            if start_q is not None and stop_q is not None and src_seq:
                human_pred = _slice_0based(src_seq, start_q, stop_q)
                human_pred_meta = {
                    "query_gene": str(first.get("query_gene") or ""),
                    "source_protein": str(first.get("source_protein") or ""),
                    "start_q": start_q,
                    "stop_q": stop_q,
                }

        # --- Predicted target domains (Ensembl/InterPro) ---
        pred_records: list[tuple[str, str]] = []
        pred_metrics: list[dict[str, object]] = []

        if human_pred:
            pred_records.append((f"homo_sapiens|{h}|predicted", human_pred))

        if not ens_h.empty:
            ens_h = ens_h[ens_h["target_species"].isin(MAMMALS)].copy()
            # Keep only Pfam PF00010 if present to reduce duplicates.
            if "Signature Accession" in ens_h.columns:
                pf = ens_h[ens_h["Signature Accession"].astype(str).eq("PF00010")]
                if not pf.empty:
                    ens_h = pf

            # Distinct domains per target gene
            ens_h["Start_T_i"] = ens_h["Start_T"].apply(_as_int)
            ens_h["Stop_T_i"] = ens_h["Stop_T"].apply(_as_int)
            ens_h = ens_h.dropna(subset=["target_seq", "target_species", "target_id", "Start_T_i", "Stop_T_i"])
            ens_h = ens_h[ens_h["Stop_T_i"] > ens_h["Start_T_i"]].copy()

            for _, r in ens_h.iterrows():
                sp = str(r["target_species"])
                tid = str(r["target_id"])
                start_t = int(r["Start_T_i"])
                stop_t = int(r["Stop_T_i"])
                tseq = str(r["target_seq"] or "")
                dom = _slice_0based(tseq, start_t, stop_t)
                if not dom:
                    continue
                name = f"{sp}|{tid}|predicted"
                pred_records.append((name, dom))

                row = {
                    "hgnc": h,
                    "species": sp,
                    "target_id": tid,
                    "start": start_t,
                    "end": stop_t,
                    "domain_len": len(dom),
                    "homology_type": str(r.get("homology_type") or ""),
                }
                if human_pred:
                    row.update(identity_metrics(human_pred, dom))
                pred_metrics.append(row)

        write_fasta(h_dir / "predicted_domains.faa", pred_records)
        df_pred_metrics = pd.DataFrame(pred_metrics)
        df_pred_metrics.to_csv(h_dir / "predicted_metrics.csv", index=False)

        # --- Projected query domains (Zoonomia) ---
        zoo_h = df_zoo[df_zoo["HGNC"].eq(h)].copy()

        if df_qgene is not None:
            # Expected columns in this project file:
            # HGNC, species, full_enst_label, q_gene
            cols = set(df_qgene.columns)
            if {"HGNC", "species", "full_enst_label", "q_gene"} <= cols:
                zoo_h = zoo_h.merge(
                    df_qgene[["HGNC", "species", "full_enst_label", "q_gene"]],
                    left_on=["HGNC", "target_species", "full_enst_label"],
                    right_on=["HGNC", "species", "full_enst_label"],
                    how="left",
                )
            # Keep it light: drop the join helper column if present.
            if "species" in zoo_h.columns:
                zoo_h = zoo_h.drop(columns=["species"])

        zoo_h["query_start_i"] = zoo_h["query_start"].apply(_as_int)
        zoo_h["query_end_i"] = zoo_h["query_end"].apply(_as_int)
        zoo_h["human_start_i"] = zoo_h["human_start"].apply(_as_int)
        zoo_h["human_end_i"] = zoo_h["human_end"].apply(_as_int)
        zoo_h = zoo_h.dropna(subset=["target_species", "full_enst_label", "query_start_i", "query_end_i"])
        zoo_h = zoo_h[zoo_h["query_end_i"] > zoo_h["query_start_i"]].copy()

        proj_records: list[tuple[str, str]] = []
        proj_metrics: list[dict[str, object]] = []

        # Add projected human reference per ENST label (useful to validate transcript-specific differences).
        projected_human_refs: dict[str, str] = {}

        for _, r in zoo_h.iterrows():
            sp = str(r["target_species"])
            full_label = str(r["full_enst_label"])
            key = ProtAlnKey(species=sp, full_enst_label=full_label)
            rec = aln.get(key)
            if rec is None or not rec.query_seq:
                continue

            qs = int(r["query_start_i"])
            qe = int(r["query_end_i"])
            dom_q = _slice_0based(rec.query_seq, qs, qe)
            if not dom_q:
                continue

            raw_q_gene = r.get("q_gene")
            if raw_q_gene is None or (isinstance(raw_q_gene, float) and math.isnan(raw_q_gene)):
                q_gene = ""
            else:
                q_gene = str(raw_q_gene).strip()
            tag = q_gene if q_gene else full_label
            name = f"{sp}|{tag}|projected"
            proj_records.append((name, dom_q))

            # ENST-specific projected human reference (if available)
            hs = _as_int(r.get("human_start_i"))
            he = _as_int(r.get("human_end_i"))
            dom_h = None
            if rec.reference_seq and hs is not None and he is not None and he > hs:
                dom_h = _slice_0based(rec.reference_seq, hs, he)
                projected_human_refs.setdefault(full_label, dom_h)

            row = {
                "hgnc": h,
                "species": sp,
                "full_enst_label": full_label,
                "q_gene": q_gene,
                "start": qs,
                "end": qe,
                "domain_len": len(dom_q),
                "homology_type": str(r.get("homology_type") or ""),
                "query_length": _as_int(r.get("query_length")),
            }
            if dom_h:
                row.update({f"ref_{k}": v for k, v in identity_metrics(dom_h, dom_q).items()})
            elif human_pred:
                row.update({f"ref_{k}": v for k, v in identity_metrics(human_pred, dom_q).items()})
            proj_metrics.append(row)

        # Write projected human refs (one per ENST label seen)
        for full_label, dom_h in sorted(projected_human_refs.items()):
            proj_records.append((f"homo_sapiens|{full_label}|projected_ref", dom_h))

        write_fasta(h_dir / "projected_domains.faa", proj_records)
        df_proj_metrics = pd.DataFrame(proj_metrics)
        df_proj_metrics.to_csv(h_dir / "projected_metrics.csv", index=False)

        # Domain-length QC (helps interpret low identities driven by truncated projections).
        pred_dom_lens: list[int] = []
        if not df_pred_metrics.empty and "domain_len" in df_pred_metrics.columns:
            pred_dom_lens = [int(x) for x in df_pred_metrics["domain_len"].dropna().astype(int).tolist() if x >= 0]

        proj_dom_lens: list[int] = []
        if not df_proj_metrics.empty and "domain_len" in df_proj_metrics.columns:
            proj_dom_lens = [int(x) for x in df_proj_metrics["domain_len"].dropna().astype(int).tolist() if x >= 0]

        proj_lt30 = [x for x in proj_dom_lens if x < 30]
        min_pred_dom = min(pred_dom_lens) if pred_dom_lens else float("nan")
        min_proj_dom = min(proj_dom_lens) if proj_dom_lens else float("nan")

        # --- Combined FASTA for optional external alignment tools ---
        combined = []
        # Prefer putting human canonical first, then human projected refs, then others.
        if human_pred:
            combined.append((f"homo_sapiens|{h}|predicted", human_pred))
        for full_label, dom_h in sorted(projected_human_refs.items()):
            combined.append((f"homo_sapiens|{full_label}|projected_ref", dom_h))
        combined.extend([x for x in pred_records if not x[0].startswith("homo_sapiens|")])
        combined.extend([x for x in proj_records if not x[0].startswith("homo_sapiens|")])
        # De-duplicate by name
        seen = set()
        combined_unique = []
        for name, seq in combined:
            if name in seen:
                continue
            seen.add(name)
            combined_unique.append((name, seq))
        write_fasta(h_dir / "all_domains.faa", combined_unique)

        # --- Cross metrics: best predicted vs projected per species ---
        pred_by_species: dict[str, list[tuple[str, str]]] = {}
        for name, seq in pred_records:
            if name.startswith("homo_sapiens|"):
                continue
            sp = name.split("|", 1)[0]
            pred_by_species.setdefault(sp, []).append((name, seq))

        proj_by_species: dict[str, list[tuple[str, str]]] = {}
        for name, seq in proj_records:
            if name.startswith("homo_sapiens|"):
                continue
            sp = name.split("|", 1)[0]
            proj_by_species.setdefault(sp, []).append((name, seq))

        cross_rows: list[dict[str, object]] = []
        for sp in sorted(set(pred_by_species) | set(proj_by_species)):
            for pn, ps in pred_by_species.get(sp, []):
                best = None
                best_row = None
                for qn, qs in proj_by_species.get(sp, []):
                    met = identity_metrics(ps, qs)
                    ident = float(met["identity"]) if _is_finite_number(met["identity"]) else float("nan")
                    if best is None or (ident == ident and ident > best):
                        best = ident
                        best_row = (qn, qs, met)
                if best_row:
                    qn, qs, met = best_row
                    row = {
                        "hgnc": h,
                        "species": sp,
                        "predicted_name": pn,
                        "projected_name": qn,
                        **{f"pred_vs_proj_{k}": v for k, v in met.items()},
                    }
                    cross_rows.append(row)
        pd.DataFrame(cross_rows).to_csv(h_dir / "cross_metrics.csv", index=False)

        # Aggregate cross identities per species (max over predicted entries).
        best_cross_by_species: dict[str, float] = {}
        for row in cross_rows:
            sp = str(row.get("species") or "")
            ident = row.get("pred_vs_proj_identity")
            try:
                v = float(ident)
            except Exception:
                continue
            if sp and math.isfinite(v):
                best_cross_by_species[sp] = max(best_cross_by_species.get(sp, float("-inf")), v)

        # --- Mini report ---
        report_path = h_dir / "report.md"
        with report_path.open("w", encoding="utf-8") as f:
            f.write(f"# Domain conservation: {h}\n\n")
            f.write(f"- Mammals subset: {', '.join(MAMMALS)}\n")
            f.write(f"- Mismatch cells: {len(mm)}\n\n")
            f.write("## Mismatch overview\n\n")
            for _, r in mm.iterrows():
                f.write(
                    f"- {r['target_species']}: Ensembl={r['ensembl_types_simple']} vs TOGA={r['toga_types']} "
                    f"(Δstart={float(r['delta_start']):.3f}, Δend={float(r['delta_end']):.3f})\n"
                )
            f.write("\n## Sequence agreement (predicted vs projected)\n\n")
            if not mm.empty:
                for sp in mm["target_species"].astype(str).tolist():
                    v = best_cross_by_species.get(sp)
                    if v is None or not math.isfinite(v):
                        f.write(f"- {sp}: n/a (missing sequence)\n")
                        continue
                    flag = "high" if v >= 0.95 else ("medium" if v >= 0.85 else "low")
                    f.write(f"- {sp}: best identity={v:.3f} ({flag} agreement)\n")
            f.write("\n## QC flags\n\n")
            f.write(
                f"- min predicted domain length: {('n/a' if not math.isfinite(min_pred_dom) else int(min_pred_dom))}\n"
            )
            f.write(
                f"- min projected domain length: {('n/a' if not math.isfinite(min_proj_dom) else int(min_proj_dom))}\n"
            )
            if proj_lt30:
                f.write(f"- WARNING: projected domain shorter than 30 aa (n={len(proj_lt30)}; min={min(proj_lt30)})\n")
            f.write("\n## Files\n\n")
            f.write("- `predicted_domains.faa` / `predicted_metrics.csv`\n")
            f.write("- `projected_domains.faa` / `projected_metrics.csv`\n")
            f.write("- `cross_metrics.csv` (best predicted-vs-projected identity per species)\n")
            f.write("\n## Quick notes\n\n")
            f.write("- Predicted sequences are extracted from `target_seq` using Start_T/Stop_T (PF00010 only when available).\n")
            f.write("- Projected sequences are extracted from Zoonomia filtered alignments using query_start/query_end.\n")

        mm_species = mm["target_species"].astype(str).tolist()
        mm_best = [best_cross_by_species.get(sp) for sp in mm_species]
        mm_best = [v for v in mm_best if v is not None and math.isfinite(float(v))]

        all_best = [v for v in best_cross_by_species.values() if v is not None and math.isfinite(float(v))]

        # overall row (for top-level scan)
        overall_rows.append(
            {
                "hgnc": h,
                "n_mismatch_cells": int(mm["target_species"].nunique()),
                "has_pan_mismatch": int(
                    not mm[mm["target_species"].eq("pan_troglodytes")].empty
                ),
                "max_delta_total": float(
                    (
                        mm["delta_start"].fillna(0).astype(float)
                        + mm["delta_end"].fillna(0).astype(float)
                    ).max()
                ),
                "min_best_identity_mismatch": (min(mm_best) if mm_best else float("nan")),
                "mean_best_identity_mismatch": (sum(mm_best) / len(mm_best) if mm_best else float("nan")),
                "n_species_compared": int(len(best_cross_by_species)),
                "min_best_identity_any": (min(all_best) if all_best else float("nan")),
                "mean_best_identity_any": (sum(all_best) / len(all_best) if all_best else float("nan")),
                "min_predicted_domain_len": min_pred_dom,
                "min_projected_domain_len": min_proj_dom,
                "n_projected_domains_lt30": int(len(proj_lt30)),
                "n_projected_domains_total": int(len(proj_dom_lens)),
                "out_dir": str(h_dir.relative_to(root)),
            }
        )

    pd.DataFrame(overall_rows).sort_values(
        ["has_pan_mismatch", "n_mismatch_cells", "max_delta_total"],
        ascending=[False, False, False],
    ).to_csv(outdir / "overall_summary.csv", index=False)

    print(f"Wrote: {outdir}")
    print(f"- overall summary: {outdir / 'overall_summary.csv'}")
    print(f"- mismatch cells: {mismatch_cells_csv}")
    print(f"- meeting shortlist: {meeting_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
