#!/usr/bin/env python3
"""04 — compute 8-gene panel score for a downloaded GEO dataset.

Usage:  python3 scripts/04_score_gene_panel.py <ACCESSION>

Reads expression block from series matrix; maps probes → genes via the platform
annotation (GPL) downloaded alongside; computes within-cohort z-score and mean.

Outputs:
  data/interim/<ACCESSION>_expression.tsv
  data/processed/<ACCESSION>_panel_score.tsv
"""
from __future__ import annotations
import gzip
import re
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"
PANEL_YAML = ROOT / "config" / "eight_gene_panel.yaml"


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_expression_block(path: Path) -> pd.DataFrame:
    in_table = False
    sample_ids: list[str] = []
    rows: list[list[str]] = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            parts = line.split("\t")
            if not sample_ids:
                sample_ids = [p.strip('"') for p in parts[1:]]
                continue
            rows.append(parts)
    if not rows:
        return pd.DataFrame()
    probe_ids = [r[0].strip('"') for r in rows]
    data = []
    for r in rows:
        vals = []
        for v in r[1:]:
            v = v.strip().strip('"')
            try:
                vals.append(float(v))
            except Exception:
                vals.append(np.nan)
        data.append(vals)
    df = pd.DataFrame(data, index=probe_ids, columns=sample_ids)
    df.index.name = "probe_id"
    return df


def load_platform_annotation(src_dir: Path) -> dict[str, str]:
    """Return probe_id → gene_symbol mapping from any GPL annotation file in src_dir."""
    mapping: dict[str, str] = {}
    # Pre-built probe2gene TSV (from _extract_gpl23159_probe2gene.py or similar helpers)
    for p2g in src_dir.glob("GPL*_probe2gene.tsv"):
        try:
            import pandas as pd
            d = pd.read_csv(p2g, sep="\t", dtype=str)
            if "ID" in d.columns and "gene_symbol" in d.columns:
                mapping.update(dict(zip(d["ID"], d["gene_symbol"])))
                print(f"  loaded probe2gene from {p2g.name}: {len(mapping)} entries")
                if mapping:
                    return mapping
        except Exception as e:
            print(f"  probe2gene TSV warning: {e}", file=sys.stderr)
    for gpl in list(src_dir.glob("GPL*.annot.gz")) + list(src_dir.glob("GPL*.annot")) + list(src_dir.glob("GPL*soft.gz")) + list(src_dir.glob("GPL*soft")):
        try:
            with open_text(gpl) as fh:
                in_table = False
                header: list[str] | None = None
                for line in fh:
                    if line.startswith("!platform_table_begin") or line.startswith("!Platform_table_begin"):
                        in_table = True
                        continue
                    if line.startswith("!platform_table_end") or line.startswith("!Platform_table_end"):
                        break
                    if not in_table:
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if header is None:
                        header = [h.lower() for h in parts]
                        continue
                    if len(parts) < len(header):
                        continue
                    row = dict(zip(header, parts))
                    pid = row.get("id") or row.get("probe_id") or row.get("id_ref")
                    sym = (row.get("gene_symbol") or row.get("symbol")
                           or row.get("gene symbol") or row.get("gene_assignment", "").split("//")[1] if "gene_assignment" in row else None)
                    if pid and sym:
                        mapping[pid.strip()] = sym.strip().split(" /// ")[0]
            if mapping:
                return mapping
        except Exception as e:
            print(f"  GPL parse warning ({gpl.name}): {e}", file=sys.stderr)
    return mapping


def collapse_probes_to_genes(expr: pd.DataFrame, probe2gene: dict[str, str]) -> pd.DataFrame:
    if not probe2gene:
        return expr  # caller can decide what to do
    e = expr.copy()
    e["gene"] = e.index.map(probe2gene)
    e = e.dropna(subset=["gene"])
    e = e[e["gene"] != ""]
    # max-variance probe per gene
    e["__var"] = e.iloc[:, :-1].var(axis=1)
    e = e.sort_values("__var", ascending=False).drop_duplicates("gene", keep="first")
    e = e.drop(columns="__var").set_index("gene")
    return e


def panel_score(gene_expr: pd.DataFrame, panel_genes: list[str]) -> pd.DataFrame:
    available = [g for g in panel_genes if g in gene_expr.index]
    missing = [g for g in panel_genes if g not in gene_expr.index]
    if not available:
        raise SystemExit(f"no panel genes found in expression matrix; panel={panel_genes}")
    sub = gene_expr.loc[available]
    # within-cohort z per gene
    zs = sub.subtract(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    panel_z = zs.mean(axis=0)
    out = pd.DataFrame({
        "panel_z": panel_z,
        "n_genes_available": len(available),
    }, index=gene_expr.columns)
    out.attrs["available"] = available
    out.attrs["missing"] = missing
    return out


def main() -> None:
    if len(sys.argv) < 2:
        print("usage: 04_score_gene_panel.py <ACCESSION>", file=sys.stderr); sys.exit(1)
    acc = sys.argv[1].strip()
    src_dir = RAW / acc
    matrices = sorted(src_dir.glob(f"{acc}*series_matrix.txt*"))
    if not matrices:
        print(f"no series_matrix in {src_dir}", file=sys.stderr); sys.exit(2)
    INTERIM.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)

    panel_cfg = yaml.safe_load(PANEL_YAML.read_text())
    panel_genes = [g["symbol"] for g in panel_cfg["genes"]]

    expr_parts = []
    for m in matrices:
        e = parse_expression_block(m)
        if not e.empty:
            print(f"  parsed {e.shape[0]} probes × {e.shape[1]} samples from {m.name}")
            expr_parts.append(e)
    if not expr_parts:
        print("no expression block found in any series_matrix", file=sys.stderr); sys.exit(3)
    # If multiple matrices (e.g. across platforms), keep largest by sample count for v1
    expr = max(expr_parts, key=lambda x: x.shape[1])
    expr.to_csv(INTERIM / f"{acc}_expression.tsv", sep="\t")

    probe2gene = load_platform_annotation(src_dir)
    print(f"# probe2gene mapping size: {len(probe2gene)}")
    if probe2gene:
        gene_expr = collapse_probes_to_genes(expr, probe2gene)
        print(f"# gene-level expression: {gene_expr.shape[0]} genes × {gene_expr.shape[1]} samples")
    else:
        # If no GPL annotation, assume probe IDs == gene symbols (best-effort)
        print("WARN: no GPL annotation; treating probe IDs as gene symbols")
        gene_expr = expr

    score = panel_score(gene_expr, panel_genes)
    available = score.attrs["available"]; missing = score.attrs["missing"]
    print(f"# panel genes available: {len(available)}/{len(panel_genes)} → {available}")
    if missing:
        print(f"# panel genes missing : {missing}")

    out = PROCESSED / f"{acc}_panel_score.tsv"
    score.to_csv(out, sep="\t", index_label="sample_id")
    print(f"# wrote {out}")


if __name__ == "__main__":
    main()
