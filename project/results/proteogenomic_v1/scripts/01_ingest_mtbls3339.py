#!/usr/bin/env python3
"""Step 1 — ingest MetaboLights MTBLS3339 ISA-Tab into clean TSVs.

Reads all i_*.txt / s_*.txt / a_*.txt / m_*.tsv / m_*.txt / m_*_maf.tsv
files in --raw-dir. Detects HTML-corrupted downloads and skips them with
a clear warning. Produces:
  - mtbls3339_sample_metadata.tsv      (one row per sample)
  - mtbls3339_metabolite_matrix_wide.tsv  (feature × sample)
  - mtbls3339_metabolite_matrix_long.tsv  (feature, sample, value)
  - mtbls3339_feature_annotation.tsv   (metabolite-level annotations)
  - mtbls3339_qc_report.md

Never modifies originals. Logs every parsing decision.
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def is_html(p: Path) -> bool:
    try:
        with p.open("rb") as f:
            head = f.read(1024).lower()
        return b"<!doctype html" in head or b"<html" in head
    except Exception:
        return True


def safe_read_tsv(p: Path) -> Tuple[pd.DataFrame | None, str]:
    if is_html(p):
        return None, f"HTML content detected — likely SPA page, not ISA-Tab data ({p.name})"
    try:
        df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False,
                         na_values=["", "NA", "nan", "NaN"], low_memory=False)
        return df, ""
    except Exception as e:
        return None, f"pd.read_csv tsv failed for {p.name}: {e}"


def parse_investigation(p: Path) -> Tuple[Dict[str, str], List[str]]:
    """ISA-Tab i_* file is a section/key/value tabbed text. Return flat dict."""
    if is_html(p):
        return {}, [f"i_ file is HTML — failed download? ({p.name})"]
    out = {}
    notes = []
    section = "ROOT"
    try:
        with p.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                # Section markers in caps with no tab
                if "\t" not in line and line.strip().isupper() and len(line) < 80:
                    section = line.strip()
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    key = parts[0].strip()
                    val = "\t".join(parts[1:]).strip()
                    if key:
                        out[f"{section}::{key}"] = val
    except Exception as e:
        notes.append(f"investigation parse failed: {e}")
    return out, notes


def detect_sample_columns(maf_df: pd.DataFrame) -> List[str]:
    """In an ISA MAF, sample-abundance columns are everything after the
    last 'standard' MAF column. Heuristic: numeric-coercible columns that are
    not in the known annotation set."""
    annotation_cols_lower = {
        "database_identifier", "chemical_formula", "smiles", "inchi",
        "metabolite_identification", "mass_to_charge", "fragmentation",
        "modifications", "charge", "retention_time", "taxid", "species",
        "database", "database_version", "reliability", "uri",
        "search_engine", "search_engine_score",
        "smallmolecule_abundance_sub", "smallmolecule_abundance_stdev_sub",
        "smallmolecule_abundance_std_error_sub",
    }
    sample_cols = []
    for c in maf_df.columns:
        if c.lower() in annotation_cols_lower:
            continue
        # try numeric coercion on first 100 non-empty
        s = maf_df[c].dropna().head(100)
        if len(s) == 0:
            continue
        n_numeric = sum(1 for v in s if str(v).replace(".", "", 1).replace("-", "", 1)
                        .replace("e", "", 1).replace("+", "", 1).isdigit() or _is_float(v))
        if n_numeric / len(s) >= 0.5:
            sample_cols.append(c)
    return sample_cols


def _is_float(x) -> bool:
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-dir", required=True, type=Path,
                    help="MTBLS3339 raw directory (with i_/s_/a_/m_ files)")
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="where to write mtbls3339_*.tsv outputs")
    args = ap.parse_args()

    if not args.raw_dir.exists():
        sys.exit(f"ERROR: --raw-dir not found: {args.raw_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    qc_lines = ["# MTBLS3339 ingest QC report",
                f"- raw-dir: `{args.raw_dir}`",
                f"- out-dir: `{args.out_dir}`", ""]

    inv_files = sorted(args.raw_dir.glob("i_*.txt"))
    sam_files = sorted(args.raw_dir.glob("s_*.txt"))
    asy_files = sorted(args.raw_dir.glob("a_*.txt"))
    maf_files = (sorted(args.raw_dir.glob("m_*.tsv"))
                 + sorted(args.raw_dir.glob("m_*.txt"))
                 + sorted(args.raw_dir.glob("m_*_maf.tsv"))
                 + sorted(args.raw_dir.glob("m_*.maf")))
    maf_files = sorted(set(maf_files))

    qc_lines.append(f"## File discovery")
    qc_lines.append(f"- investigation files: {[p.name for p in inv_files]}")
    qc_lines.append(f"- sample files: {[p.name for p in sam_files]}")
    qc_lines.append(f"- assay files: {[p.name for p in asy_files]}")
    qc_lines.append(f"- metabolite MAF files: {[p.name for p in maf_files]}")
    qc_lines.append("")

    # ---- Investigation ----
    if inv_files:
        inv, notes = parse_investigation(inv_files[0])
        qc_lines.append("## Investigation")
        if notes:
            qc_lines.append("- WARNINGS: " + "; ".join(notes))
        qc_lines.append(f"- parsed {len(inv)} key/value pairs from {inv_files[0].name}")
        # Save flat
        inv_path = args.out_dir / "mtbls3339_investigation_flat.tsv"
        with inv_path.open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["section_key", "value"])
            for k, v in sorted(inv.items()):
                w.writerow([k, v])
        qc_lines.append(f"- wrote {inv_path.name}")
    else:
        qc_lines.append("## Investigation\n- no i_*.txt found.")

    # ---- Samples ----
    sample_meta = None
    if sam_files:
        df, err = safe_read_tsv(sam_files[0])
        if df is None:
            qc_lines.append(f"## Samples\n- ERROR reading {sam_files[0].name}: {err}")
        else:
            sample_meta = df
            out_path = args.out_dir / "mtbls3339_sample_metadata.tsv"
            df.to_csv(out_path, sep="\t", index=False)
            qc_lines.append(f"## Samples\n- {sam_files[0].name}: {len(df)} rows × {len(df.columns)} cols")
            qc_lines.append(f"- wrote {out_path.name}")
            qc_lines.append(f"- columns: {list(df.columns)[:25]}{'...' if len(df.columns) > 25 else ''}")
    else:
        qc_lines.append("## Samples\n- no s_*.txt found.")

    # ---- Assay ----
    if asy_files:
        for ap_ in asy_files:
            df, err = safe_read_tsv(ap_)
            if df is None:
                qc_lines.append(f"## Assay {ap_.name}\n- ERROR: {err}")
                continue
            out_path = args.out_dir / f"mtbls3339_assay_{ap_.stem}.tsv"
            df.to_csv(out_path, sep="\t", index=False)
            qc_lines.append(f"## Assay {ap_.name}\n- {len(df)} rows × {len(df.columns)} cols → {out_path.name}")
            qc_lines.append(f"- columns: {list(df.columns)[:15]}")
    else:
        qc_lines.append("## Assay\n- no a_*.txt found.")

    # ---- Metabolite abundance ----
    all_long = []
    feat_rows = []
    if not maf_files:
        qc_lines.append("## Metabolite matrices\n- NO m_*.tsv FILES FOUND. Cannot build abundance matrix.")
    for mp in maf_files:
        df, err = safe_read_tsv(mp)
        if df is None:
            qc_lines.append(f"## MAF {mp.name}\n- ERROR: {err}")
            continue
        sample_cols = detect_sample_columns(df)
        anno_cols = [c for c in df.columns if c not in sample_cols]
        qc_lines.append(f"## MAF {mp.name}")
        qc_lines.append(f"- {len(df)} features × {len(sample_cols)} sample-cols (detected); annotation cols: {len(anno_cols)}")
        qc_lines.append(f"- sample col preview: {sample_cols[:5]}")
        # Feature key
        key_cols = [c for c in ("metabolite_identification", "database_identifier",
                                "Metabolite identification", "Metabolite_identification") if c in df.columns]
        feat_key = key_cols[0] if key_cols else df.columns[0]
        df["__feature_id__"] = df[feat_key].astype(str) + "@" + mp.stem
        # Wide: per-MAF wide piece
        wide_piece = df.set_index("__feature_id__")[sample_cols]
        # Coerce numeric
        wide_piece = wide_piece.apply(pd.to_numeric, errors="coerce")
        # Long
        long_piece = wide_piece.stack(dropna=False).rename("value").reset_index()
        long_piece.columns = ["feature_id", "sample_id", "value"]
        long_piece["source_maf"] = mp.name
        all_long.append(long_piece)
        # Feature annotation
        anno = df.set_index("__feature_id__")[anno_cols].copy()
        anno["source_maf"] = mp.name
        feat_rows.append(anno.reset_index())
        # Per-MAF wide TSV (so wide doesn't conflate platforms)
        out_wide = args.out_dir / f"mtbls3339_metabolite_matrix_wide_{mp.stem}.tsv"
        wide_piece.reset_index().to_csv(out_wide, sep="\t", index=False)
        qc_lines.append(f"- wrote {out_wide.name}")

    if all_long:
        long_df = pd.concat(all_long, ignore_index=True)
        long_path = args.out_dir / "mtbls3339_metabolite_matrix_long.tsv"
        long_df.to_csv(long_path, sep="\t", index=False)
        qc_lines.append(f"\n## Aggregated long matrix: {len(long_df)} rows → {long_path.name}")

        # Build a single "wide" by outer join
        wide_concat = (long_df.pivot_table(index="feature_id", columns="sample_id",
                                          values="value", aggfunc="first"))
        wide_path = args.out_dir / "mtbls3339_metabolite_matrix_wide.tsv"
        wide_concat.reset_index().to_csv(wide_path, sep="\t", index=False)
        qc_lines.append(f"## Aggregated wide matrix: {wide_concat.shape} → {wide_path.name}")

    if feat_rows:
        feat_df = pd.concat(feat_rows, ignore_index=True)
        feat_path = args.out_dir / "mtbls3339_feature_annotation.tsv"
        feat_df.to_csv(feat_path, sep="\t", index=False)
        qc_lines.append(f"## Feature annotation: {len(feat_df)} → {feat_path.name}")

    qc_path = args.out_dir / "mtbls3339_qc_report.md"
    qc_path.write_text("\n".join(qc_lines) + "\n")
    print(f"[01] QC → {qc_path}")
    if not maf_files:
        print("[01] WARNING: no MAF files — abundance matrices NOT generated.", file=sys.stderr)


if __name__ == "__main__":
    main()
