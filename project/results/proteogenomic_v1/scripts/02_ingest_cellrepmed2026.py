#!/usr/bin/env python3
"""Step 2 — ingest Cell Rep Med 2026 advanced DTC supplements.

Discovers all .xlsx/.xls/.pdf in --raw-dir, lists every sheet, then
*best-effort* categorizes sheets into clinical / subtype / protein
matrix / feature annotation. **Does not auto-merge.** Each sheet is
emitted as-is with provenance, so the user can verify before joining.

For PDFs: pdftotext extract first 5 pages (caption hints only).

Outputs:
  - cellrepmed2026_sheet_inventory.tsv  (every sheet × file with first-row preview)
  - cellrepmed2026_clinical.tsv         (heuristic clinical sheets, concatenated with provenance)
  - cellrepmed2026_subtypes.tsv         (heuristic subtype/cluster call sheets)
  - cellrepmed2026_protein_matrix_wide.tsv (heuristic large numeric matrix sheets)
  - cellrepmed2026_feature_annotation.tsv (gene/protein annotation sheets)
  - cellrepmed2026_qc_report.md
"""
from __future__ import annotations
import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List

import pandas as pd

CLINICAL_HINTS = {"patient", "patient_id", "case", "caseid", "case_id", "id",
                  "age", "sex", "gender", "stage", "tnm", "histology",
                  "os", "os_status", "survival", "follow", "outcome",
                  "recurrence", "rfs", "dfs", "tert", "braf", "tert_promoter"}
SUBTYPE_HINTS = {"cc1", "cc2", "cc3", "cluster", "subtype", "consensus",
                 "canonical", "stromal", "immunogenic", "ecotype"}
PROTEIN_FEATURE_HINTS = {"gene", "gene_symbol", "gene_name", "protein",
                         "uniprot", "ensembl", "phospho", "phospho_site",
                         "site", "peptide", "sequence", "modification"}


def categorize_sheet(headers: List[str], second_row: List[str], n_rows: int,
                     n_cols: int) -> List[str]:
    cats = []
    headers_l = [str(h).lower() if h is not None else "" for h in headers]
    headers_set = set(headers_l)
    if any(h in CLINICAL_HINTS for h in headers_l):
        cats.append("clinical")
    second_l = [str(v).lower() if v is not None else "" for v in second_row]
    if any(h in SUBTYPE_HINTS for h in headers_l) or \
       any(v in {"cc1", "cc2", "cc3"} for v in second_l):
        cats.append("subtype")
    # Protein matrix heuristic: lots of cols + first col is gene-like
    if n_cols >= 10 and any(h in PROTEIN_FEATURE_HINTS for h in headers_l[:3]):
        cats.append("protein_matrix")
    if any(h in PROTEIN_FEATURE_HINTS for h in headers_l) and n_cols < 15:
        cats.append("feature_annotation")
    if not cats:
        cats.append("uncategorized")
    return cats


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    if not args.raw_dir.exists():
        sys.exit(f"ERROR: --raw-dir not found: {args.raw_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    qc_lines = ["# Cell Rep Med 2026 ingest QC report",
                f"- raw-dir: `{args.raw_dir}`", f"- out-dir: `{args.out_dir}`", ""]

    xlsx_files = sorted(args.raw_dir.glob("*.xlsx")) + sorted(args.raw_dir.glob("*.xls"))
    pdf_files = sorted(args.raw_dir.glob("*.pdf"))
    qc_lines.append(f"## File discovery")
    qc_lines.append(f"- xlsx files: {[p.name for p in xlsx_files]}")
    qc_lines.append(f"- pdf files: {[p.name for p in pdf_files]}")
    qc_lines.append("")

    if not xlsx_files and not pdf_files:
        qc_lines.append("**No supplement files found. Nothing to ingest.**")
        (args.out_dir / "cellrepmed2026_qc_report.md").write_text("\n".join(qc_lines) + "\n")
        print("[02] no xlsx/pdf in raw-dir; QC written.")
        return

    # ---- Sheet inventory ----
    try:
        import openpyxl
    except ImportError:
        sys.exit("ERROR: openpyxl not installed (pip install openpyxl)")

    inv_rows = []
    clin_pieces = []
    sub_pieces = []
    prot_pieces = []
    feat_pieces = []

    for xp in xlsx_files:
        try:
            wb = openpyxl.load_workbook(str(xp), read_only=True, data_only=True)
        except Exception as e:
            qc_lines.append(f"## {xp.name}\n- ERROR loading: {e}")
            inv_rows.append({"file": xp.name, "sheet": "<load-failed>", "error": str(e),
                             "n_rows": -1, "n_cols": -1, "categories": "", "preview": ""})
            continue
        qc_lines.append(f"## {xp.name}  ({len(wb.sheetnames)} sheets)")
        for sname in wb.sheetnames:
            try:
                ws = wb[sname]
                rows_iter = ws.iter_rows(values_only=True)
                first = next(rows_iter, None) or []
                second = next(rows_iter, None) or []
                n_rows = ws.max_row or 0
                n_cols = ws.max_column or 0
                cats = categorize_sheet(list(first), list(second), n_rows, n_cols)
                preview = " | ".join([str(v)[:18] if v is not None else "" for v in first[:8]])
                inv_rows.append({"file": xp.name, "sheet": sname,
                                 "n_rows": n_rows, "n_cols": n_cols,
                                 "categories": ",".join(cats),
                                 "preview": preview, "error": ""})
                qc_lines.append(f"  - [{sname}] {n_rows}×{n_cols} cats={','.join(cats)}")
                # Full read for categorized sheets (only if reasonable size)
                if any(c in {"clinical", "subtype", "protein_matrix", "feature_annotation"}
                       for c in cats) and n_rows * n_cols < 5_000_000:
                    try:
                        df = pd.read_excel(str(xp), sheet_name=sname, dtype=object)
                        df["__source_file__"] = xp.name
                        df["__source_sheet__"] = sname
                        if "clinical" in cats:
                            clin_pieces.append(df)
                        if "subtype" in cats:
                            sub_pieces.append(df)
                        if "protein_matrix" in cats:
                            prot_pieces.append(df)
                        if "feature_annotation" in cats:
                            feat_pieces.append(df)
                    except Exception as e:
                        qc_lines.append(f"    WARN read_excel({sname}) failed: {e}")
            except Exception as e:
                qc_lines.append(f"  - [{sname}] read error: {e}")
                inv_rows.append({"file": xp.name, "sheet": sname,
                                 "n_rows": -1, "n_cols": -1,
                                 "categories": "", "preview": "", "error": str(e)})

    inv_path = args.out_dir / "cellrepmed2026_sheet_inventory.tsv"
    with inv_path.open("w", newline="") as f:
        fields = ["file", "sheet", "n_rows", "n_cols", "categories", "preview", "error"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in inv_rows:
            w.writerow(r)
    qc_lines.append(f"\n## Wrote sheet inventory → {inv_path.name} ({len(inv_rows)} sheets)")

    def _save_pieces(name: str, pieces, fname: str):
        if not pieces:
            qc_lines.append(f"- {name}: 0 sheets matched.")
            return
        cat = pd.concat(pieces, ignore_index=True)
        out = args.out_dir / fname
        cat.to_csv(out, sep="\t", index=False)
        qc_lines.append(f"- {name}: {len(pieces)} sheets, {len(cat)} rows → {out.name}")

    qc_lines.append("\n## Categorized output")
    _save_pieces("clinical", clin_pieces, "cellrepmed2026_clinical.tsv")
    _save_pieces("subtypes", sub_pieces, "cellrepmed2026_subtypes.tsv")
    _save_pieces("protein matrix", prot_pieces, "cellrepmed2026_protein_matrix_wide.tsv")
    _save_pieces("feature annotation", feat_pieces, "cellrepmed2026_feature_annotation.tsv")

    # ---- PDFs ----
    if pdf_files:
        qc_lines.append("\n## PDFs")
        for pp in pdf_files:
            qc_lines.append(f"- {pp.name} ({pp.stat().st_size} B)")
            if shutil.which("pdftotext"):
                try:
                    txt_path = args.out_dir / f"cellrepmed2026_{pp.stem}_first5p.txt"
                    subprocess.check_call(
                        ["pdftotext", "-layout", "-f", "1", "-l", "5",
                         str(pp), str(txt_path)],
                        stderr=subprocess.DEVNULL, timeout=30)
                    qc_lines.append(f"  - first-5-pages text → {txt_path.name}")
                except Exception as e:
                    qc_lines.append(f"  - pdftotext failed: {e}")
            else:
                qc_lines.append("  - pdftotext not installed; skipping text extraction")

    qc_path = args.out_dir / "cellrepmed2026_qc_report.md"
    qc_path.write_text("\n".join(qc_lines) + "\n")
    print(f"[02] QC → {qc_path}")


if __name__ == "__main__":
    main()
