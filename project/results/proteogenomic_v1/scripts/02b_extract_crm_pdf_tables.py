#!/usr/bin/env python3
"""Step 2b — extract Tables S1–S5 from Cell Rep Med 2026 supplementary PDF.

The user's downloaded `mmc1.pdf` (and the article `mmc2.pdf`) are
PDF-only — Cell does not always publish per-table xlsx for Cell Rep
Med supplements. The PDFs DO contain real tables in fixed-width
column layout (see `cellrepmed2026_mmc1_supp_info_full.txt`).

This script reads the previously-extracted full text and parses each
"Table S{N}. <title>" block into a TSV by splitting on runs of 2+
spaces. Output naming:
  - cellrepmed2026_table_S{N}_<slug>.tsv
  - cellrepmed2026_pdf_table_extraction_qc.md

Reads only — never modifies originals.
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
from pathlib import Path
from typing import List, Tuple

# Map table number → conservative output slug + per-table notes.
TABLE_SLUGS = {
    "S1": ("clinical", "Discovery Cohort 113 patients × clinical + CC subtype"),
    "S2": ("ngs", "NGS mutation calls per patient"),
    "S3": ("scrna_spatial_cohort", "External Cohort 1 — scRNA + spatial RNA-seq"),
    "S4": ("realworld_treatment_cohort", "External Cohort 2 — real-world systemic therapy"),
    "S5": ("signature_gene_sets", "TDS / stromal / immune signature gene lists"),
    "S6": ("supp_table_6", "TBD"),
    "S7": ("supp_table_7", "TBD"),
}

TABLE_HEADER_RE = re.compile(
    r"^[\s\f]*Table\s+(S\d+)\.\s*(.+?)\s*$", re.MULTILINE
)


def split_row(line: str) -> List[str]:
    """Split a fixed-width-aligned PDF text line into cells via 2+ spaces."""
    # Strip leading/trailing whitespace and form-feeds first.
    cleaned = line.replace("\f", "").rstrip()
    if not cleaned.strip():
        return []
    return [c.strip() for c in re.split(r"\s{2,}", cleaned.strip()) if c.strip() != ""]


def extract_blocks(text: str) -> List[Tuple[str, str, int, int]]:
    """Find (table_id, title, start_line, end_line_exclusive) blocks."""
    lines = text.splitlines()
    headers = []
    for i, ln in enumerate(lines):
        m = TABLE_HEADER_RE.match(ln)
        if m:
            headers.append((m.group(1), m.group(2), i))
    blocks = []
    for idx, (tid, title, start) in enumerate(headers):
        end = headers[idx + 1][2] if idx + 1 < len(headers) else len(lines)
        blocks.append((tid, title, start, end))
    return blocks


def parse_block(lines: List[str]) -> Tuple[List[str], List[List[str]], List[str]]:
    """Given the lines of a single Table block (header line + data),
    return (column_headers, data_rows, warnings)."""
    warnings: List[str] = []
    # Find the header row: first non-blank line whose split has >=2 cells.
    hdr_idx = None
    hdr_cells: List[str] = []
    for i, ln in enumerate(lines):
        cells = split_row(ln)
        if len(cells) >= 2:
            hdr_idx = i
            hdr_cells = cells
            break
    if hdr_idx is None:
        return [], [], ["No header row found"]
    n_cols = len(hdr_cells)
    rows: List[List[str]] = []
    for ln in lines[hdr_idx + 1:]:
        cells = split_row(ln)
        if not cells:
            continue
        # Skip caption / explanatory continuation lines that look like sentences:
        # detect by: only 1 cell OR cell count is small AND contains lowercase-only words
        if len(cells) == 1:
            warnings.append(f"single-cell line treated as continuation: {cells[0][:60]}")
            continue
        # Pad / truncate
        if len(cells) < n_cols:
            cells = cells + [""] * (n_cols - len(cells))
        elif len(cells) > n_cols:
            # Common when a value contains ≥2 spaces. Merge overflow into last column.
            cells = cells[: n_cols - 1] + [" | ".join(cells[n_cols - 1:])]
        rows.append(cells)
    return hdr_cells, rows, warnings


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--text-file", required=True, type=Path,
                    help="output of pdftotext on mmc1.pdf (e.g. cellrepmed2026_mmc1_supp_info_full.txt)")
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    if not args.text_file.exists():
        sys.exit(f"ERROR: text file not found: {args.text_file}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    text = args.text_file.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    blocks = extract_blocks(text)

    qc_lines = ["# CRM PDF table extraction QC",
                f"- input: `{args.text_file}`", f"- out-dir: `{args.out_dir}`",
                f"- detected {len(blocks)} 'Table S{{N}}.' blocks", ""]

    if not blocks:
        qc_lines.append("**No tables detected. Layout extraction may have failed.**")
        (args.out_dir / "cellrepmed2026_pdf_table_extraction_qc.md").write_text(
            "\n".join(qc_lines) + "\n")
        return

    for tid, title, start, end in blocks:
        block_lines = lines[start + 1:end]   # skip the "Table S{N}." line itself
        hdr, rows, warns = parse_block(block_lines)
        slug, note = TABLE_SLUGS.get(tid, (f"unknown_{tid.lower()}", "unmapped"))
        out_path = args.out_dir / f"cellrepmed2026_table_{tid}_{slug}.tsv"
        with out_path.open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            if hdr:
                w.writerow(hdr)
            for r in rows:
                w.writerow(r)
        qc_lines.append(f"## Table {tid} — {title}")
        qc_lines.append(f"- expected: {note}")
        qc_lines.append(f"- header ({len(hdr)} cols): {hdr}")
        qc_lines.append(f"- data rows: {len(rows)}")
        if warns:
            qc_lines.append(f"- warnings ({len(warns)}): {warns[:3]}")
        qc_lines.append(f"- wrote: `{out_path.name}`\n")
        print(f"[02b] Table {tid}: {len(rows)} rows × {len(hdr)} cols → {out_path.name}")

    qc_path = args.out_dir / "cellrepmed2026_pdf_table_extraction_qc.md"
    qc_path.write_text("\n".join(qc_lines) + "\n")
    print(f"[02b] QC → {qc_path}")


if __name__ == "__main__":
    main()
