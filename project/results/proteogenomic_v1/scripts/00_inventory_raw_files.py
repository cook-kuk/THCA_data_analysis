#!/usr/bin/env python3
"""Step 0 — inventory every file under --raw-root.

Writes processed/raw_file_manifest.tsv with one row per file. For .xlsx,
sheet names + dims are listed. For .tsv/.txt/.csv, line+column counts.
For .pdf, page count via pdfinfo if available. HTML-content text files
(common when a download URL silently returns the SPA) are flagged.

Reads only — never modifies originals.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

TEXT_EXTS = {".txt", ".tsv", ".csv", ".tab", ".maf"}
EXCEL_EXTS = {".xlsx", ".xls", ".xlsm"}


def sha256_head(p: Path, n_bytes: int = 262_144) -> str:
    h = hashlib.sha256()
    try:
        with p.open("rb") as f:
            h.update(f.read(n_bytes))
        return h.hexdigest()
    except Exception as e:
        return f"err:{e.__class__.__name__}"


def is_html(p: Path, n_bytes: int = 1024) -> bool:
    try:
        with p.open("rb") as f:
            head = f.read(n_bytes).lower()
        return b"<!doctype html" in head or b"<html" in head
    except Exception:
        return False


def inspect_excel(p: Path) -> dict:
    try:
        import openpyxl
    except ImportError:
        return {"kind": "xlsx", "warning": "openpyxl_not_installed"}
    try:
        wb = openpyxl.load_workbook(str(p), read_only=True, data_only=True)
        sheets = []
        for name in wb.sheetnames:
            ws = wb[name]
            sheets.append({"name": name, "dim": ws.calculate_dimension(),
                           "max_row": ws.max_row, "max_col": ws.max_column})
        return {"kind": "xlsx", "n_sheets": len(sheets), "sheets": sheets}
    except Exception as e:
        return {"kind": "xlsx", "warning": f"openpyxl_failed:{e}"}


def inspect_text(p: Path) -> dict:
    if is_html(p):
        return {"kind": "text_html_failed_download",
                "warning": "first 1KB looks like HTML — likely a download error, not real data"}
    n_lines = 0
    first_line = ""
    n_cols = None
    delim_guess = "?"
    try:
        with p.open("r", encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f):
                if i == 0:
                    first_line = line.rstrip("\n")
                    if "\t" in first_line:
                        delim_guess = "tab"
                    elif "," in first_line:
                        delim_guess = "comma"
                    else:
                        delim_guess = "ws"
                    cols = first_line.split("\t" if delim_guess == "tab"
                                            else ("," if delim_guess == "comma" else None))
                    n_cols = len(cols)
                n_lines = i + 1
        return {"kind": "text", "n_lines": n_lines, "n_cols": n_cols,
                "delim_guess": delim_guess, "first_line": first_line[:200]}
    except Exception as e:
        return {"kind": "text", "warning": f"read_failed:{e}"}


def inspect_pdf(p: Path) -> dict:
    info = {"kind": "pdf"}
    if shutil.which("pdfinfo"):
        try:
            out = subprocess.check_output(["pdfinfo", str(p)],
                                          stderr=subprocess.DEVNULL,
                                          timeout=15).decode()
            for line in out.splitlines():
                if line.startswith("Pages:"):
                    info["n_pages"] = int(line.split(":")[1].strip())
                    break
        except Exception as e:
            info["warning"] = f"pdfinfo_failed:{e}"
    else:
        info["warning"] = "pdfinfo_not_installed"
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-root", required=True, type=Path,
                    help="root directory containing raw subdirs (MTBLS3339/, cellrepmed2026/, ...)")
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="where to write raw_file_manifest.tsv")
    args = ap.parse_args()

    if not args.raw_root.exists():
        sys.exit(f"ERROR: --raw-root does not exist: {args.raw_root}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for p in sorted(args.raw_root.rglob("*")):
        if p.is_dir():
            continue
        rel = p.relative_to(args.raw_root)
        ext = p.suffix.lower()
        try:
            size = p.stat().st_size
            mtime = p.stat().st_mtime
        except Exception as e:
            rows.append({"relpath": str(rel), "abspath": str(p), "ext": ext,
                         "size_bytes": -1, "mtime": -1, "kind": "stat_failed",
                         "warning": str(e), "details_json": ""})
            continue
        digest = sha256_head(p)
        if ext in EXCEL_EXTS:
            details = inspect_excel(p)
        elif ext in TEXT_EXTS:
            details = inspect_text(p)
        elif ext == ".pdf":
            details = inspect_pdf(p)
        else:
            # Some MTBLS files have no extension or .raw etc.
            details = {"kind": f"other_{ext or 'noext'}"}
            if size < 5_000_000 and ext not in {".mzml", ".raw", ".d", ".gz", ".zip"}:
                # Try treating as text for files that might be ISA-Tab without known ext
                if is_html(p):
                    details["warning"] = "first 1KB HTML — likely failed download"
        rows.append({
            "relpath": str(rel),
            "abspath": str(p),
            "ext": ext,
            "size_bytes": size,
            "mtime": mtime,
            "sha256_head_256k": digest,
            "kind": details.pop("kind", "?"),
            "warning": details.pop("warning", ""),
            "details_json": json.dumps(details, ensure_ascii=False),
        })

    out_path = args.out_dir / "raw_file_manifest.tsv"
    fields = ["relpath", "abspath", "ext", "size_bytes", "mtime",
              "sha256_head_256k", "kind", "warning", "details_json"]
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    n_html = sum(1 for r in rows if "html" in r["kind"] or "html" in r["warning"])
    n_xlsx = sum(1 for r in rows if r["ext"] in EXCEL_EXTS)
    n_text = sum(1 for r in rows if r["ext"] in TEXT_EXTS)
    n_pdf = sum(1 for r in rows if r["ext"] == ".pdf")
    print(f"[00] wrote {out_path}: {len(rows)} files "
          f"(xlsx={n_xlsx}, text={n_text}, pdf={n_pdf}, html_failed_dl={n_html})")
    if n_html:
        print(f"[00] WARNING: {n_html} files look like HTML downloads, not real data:",
              file=sys.stderr)
        for r in rows:
            if "html" in r["kind"] or "html" in r["warning"]:
                print(f"     {r['relpath']}  ({r['size_bytes']} B)", file=sys.stderr)


if __name__ == "__main__":
    main()
