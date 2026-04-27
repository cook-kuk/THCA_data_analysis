#!/usr/bin/env python3
"""
v17 REAL FIX · R3 — external validation with TRUE ground truth.

R3-A: GSE76039 ATC vs PDTC (histology = ground truth, processed local data exists)
R3-B: GSE151180 RAI-refractory vs RAI-avid (TRUE clinical ground truth)

Outputs:
  - results/v17_realfix/R3A_gse76039_predictions.tsv
  - results/v17_realfix/R3_summary.json
"""
from __future__ import annotations
import io
import json
import re
import sys
import gzip
from pathlib import Path
from urllib.request import urlopen, Request

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

OUT = Path("/opt/thyroid-dash/project/results/v17_realfix")
TPM = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
GSE76039_EXPR = Path("/data/thca/data_processed/microarray/GSE76039_microarray_expression_log2.tsv")
GSE76039_SM = Path("/data/thca/data_raw/geo/GSE76039/GSE76039_series_matrix.txt.gz")
GENE_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]


def fetch_bytes(url: str, timeout: int = 600) -> bytes:
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=timeout) as r:
        return r.read()


def load_r1a_labels() -> pd.DataFrame:
    df = pd.read_csv(OUT / "R1A_cluster_labels.tsv", sep="\t")
    df["y"] = df["cluster"].str.startswith("DM2").astype(int)
    return df


def train_tcga_logreg() -> tuple[LogisticRegression, StandardScaler, list[str]]:
    expr = pd.read_csv(TPM, sep="\t", index_col=0)
    lbl = load_r1a_labels()
    g_have = [g for g in GENE_8 if g in expr.index]
    Xs = expr.loc[g_have, [s for s in lbl["sample_id"] if s in expr.columns]].T
    common = Xs.index.intersection(lbl.set_index("sample_id").index)
    X = Xs.loc[common].values
    y = lbl.set_index("sample_id").loc[common, "y"].values
    sc = StandardScaler().fit(X)
    m = LogisticRegression(C=1.0, max_iter=2000, random_state=42).fit(sc.transform(X), y)
    print(f"[train] LogReg on n={len(common)} TCGA samples; gene order: {g_have}", flush=True)
    return m, sc, g_have


def parse_series_matrix(path: Path) -> tuple[list[str], dict]:
    """Returns (gsm_ids, {char_col_name: [values]}) and the source_name row."""
    with gzip.open(path, "rb") as f:
        txt = f.read().decode("utf-8", errors="replace")
    lines = txt.split("\n")
    info: dict = {}
    gsms = []
    for ln in lines:
        if ln.startswith("!Sample_geo_accession"):
            gsms = [v.strip().strip('"') for v in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_source_name_ch1"):
            info["source_name"] = [v.strip().strip('"') for v in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_title"):
            info["title"] = [v.strip().strip('"') for v in ln.split("\t")[1:]]
    return gsms, info


def r3a_gse76039(model, scaler, gene_order: list[str]) -> dict:
    print("\n=== R3-A · GSE76039 ATC vs PDTC ===", flush=True)
    expr = pd.read_csv(GSE76039_EXPR, sep="\t", index_col=0)
    print(f"[R3A] expression {expr.shape}, sample cols: {list(expr.columns)[:3]}...", flush=True)
    gsms, meta = parse_series_matrix(GSE76039_SM)
    src = meta.get("source_name", [])
    if len(src) != len(gsms):
        return {"error": "source_name length mismatch", "n_gsms": len(gsms), "n_src": len(src)}
    df_meta = pd.DataFrame({"sample_id": gsms, "source": src})
    df_meta["histology"] = df_meta["source"].astype(str).str.lower()
    df_meta["y_ATC"] = df_meta["histology"].str.contains("anaplastic", na=False).astype(int)
    df_meta["is_pdtc"] = df_meta["histology"].str.contains("poorly", na=False).astype(int)
    keep = (df_meta["y_ATC"] == 1) | (df_meta["is_pdtc"] == 1)
    print(f"[R3A] ATC={int(df_meta['y_ATC'].sum())} PDTC={int(df_meta['is_pdtc'].sum())} kept={int(keep.sum())}", flush=True)
    df_keep = df_meta[keep].reset_index(drop=True)

    g_have = [g for g in gene_order if g in expr.index]
    print(f"[R3A] matched 8-gene: {len(g_have)}/{len(gene_order)} ({g_have})", flush=True)
    if len(g_have) < 4:
        return {"error": f"too few matched genes: {g_have}"}

    samples_in_expr = [s for s in df_keep["sample_id"] if s in expr.columns]
    Xext = expr.loc[g_have, samples_in_expr].T  # samples × genes
    Xfull = pd.DataFrame(0.0, index=Xext.index, columns=gene_order)
    for g in g_have:
        Xfull[g] = Xext[g]
    Xt = scaler.transform(Xfull.values)
    p = model.predict_proba(Xt)[:, 1]
    y = df_keep.set_index("sample_id").loc[samples_in_expr, "y_ATC"].astype(int).values
    if len(np.unique(y)) < 2:
        return {"error": "single-class y", "y_unique": list(np.unique(y))}
    auc = roc_auc_score(y, p)
    auc_dc = max(auc, 1 - auc)
    rng = np.random.default_rng(42)
    aucs = []
    for _ in range(1000):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) >= 2:
            aucs.append(roc_auc_score(y[idx], p[idx]))
    ci = (float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))) if aucs else (None, None)
    pd.DataFrame({"sample_id": samples_in_expr, "y_ATC": y, "p_DM2": p}).to_csv(OUT / "R3A_gse76039_predictions.tsv", sep="\t", index=False)
    out = {
        "n": int(len(samples_in_expr)),
        "n_ATC": int(y.sum()),
        "n_PDTC": int(len(y) - y.sum()),
        "n_genes_matched": len(g_have),
        "auc_raw": float(auc),
        "auc_direction_corrected": float(auc_dc),
        "auc_95ci": list(ci),
    }
    print(f"[R3A] AUC raw {auc:.4f}  direction-corrected {auc_dc:.4f}  95% CI [{ci[0]:.3f}, {ci[1]:.3f}]", flush=True)
    return out


def r3b_gse151180(model, scaler, gene_order: list[str]) -> dict:
    print("\n=== R3-B · GSE151180 RAI-refractory vs RAI-avid ===", flush=True)
    # try GEO direct download
    base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151180/matrix/"
    try:
        idx_html = fetch_bytes(base, timeout=30).decode("utf-8", errors="replace")
        files = re.findall(r'href="(GSE151180[^"]*?_series_matrix\.txt\.gz)"', idx_html)
        if not files:
            return {"error": "no series matrix file under " + base, "html_excerpt": idx_html[:500]}
        sm_url = base + files[0]
        print(f"[R3B] downloading {sm_url}", flush=True)
        data = fetch_bytes(sm_url, timeout=300)
        txt = gzip.decompress(data).decode("utf-8", errors="replace")
    except Exception as e:
        return {"error": f"download failed: {e}"}

    lines = txt.split("\n")
    gsms = []
    sample_titles = []
    chars: list[list[str]] = []
    matrix_start = None
    matrix_end = None
    for i, ln in enumerate(lines):
        if ln.startswith("!Sample_geo_accession"):
            gsms = [v.strip().strip('"') for v in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_title"):
            sample_titles = [v.strip().strip('"') for v in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            chars.append([v.strip().strip('"') for v in ln.split("\t")[1:]])
        elif ln.strip() == "!series_matrix_table_begin":
            matrix_start = i + 1
        elif ln.strip() == "!series_matrix_table_end":
            matrix_end = i

    print(f"[R3B] n_samples = {len(gsms)}", flush=True)
    if not chars:
        return {"error": "no characteristics rows"}

    # find RAI-related characteristic
    rai_char = None
    for c in chars:
        if any(re.search(r"refract|avid|rair|radioiod|iodine.*resp", v, re.I) for v in c if v):
            rai_char = c
            break
    if rai_char is None:
        # show samples of chars
        return {"error": "no RAI char row", "chars_examples": [c[:3] for c in chars[:10]]}

    y = []
    for v in rai_char:
        v_low = v.lower()
        if "refract" in v_low or "rair" in v_low or "non-avid" in v_low:
            y.append(1)
        elif "avid" in v_low:
            y.append(0)
        else:
            y.append(-1)
    y = np.array(y)
    keep = y >= 0
    print(f"[R3B] resolved labels — refractory={int((y==1).sum())} avid={int((y==0).sum())} other={int((y==-1).sum())}", flush=True)
    if keep.sum() < 8:
        return {"error": "too few labelled samples", "rai_char_examples": rai_char[:5]}

    # parse matrix
    if matrix_start is None or matrix_end is None:
        return {"error": "no matrix block markers"}
    matrix_text = "\n".join(lines[matrix_start:matrix_end])
    df = pd.read_csv(io.StringIO(matrix_text), sep="\t", index_col=0)
    print(f"[R3B] matrix shape {df.shape}, sample probes? {df.index[:5].tolist()}", flush=True)

    g_in_idx = [g for g in gene_order if g in df.index]
    print(f"[R3B] direct gene-symbol match: {len(g_in_idx)}/{len(gene_order)}", flush=True)
    if len(g_in_idx) < 4:
        # microarray: probe-level → need GPL annotation
        return {
            "error": "probe-level matrix without symbol mapping in this script",
            "n_total_samples": int(len(gsms)),
            "n_labelled": int(keep.sum()),
            "n_refractory": int((y == 1).sum()),
            "n_avid": int((y == 0).sum()),
            "matrix_shape": list(df.shape),
            "matrix_index_examples": list(df.index[:5]),
            "note_for_revision": "GSE151180 is microarray. Need GPL probe→symbol mapping (likely GPL21575 Agilent). Easy 1-day extension for revision round.",
        }
    samples_kept = [g for i, g in enumerate(gsms) if keep[i] and g in df.columns]
    Xext = df.loc[g_in_idx, samples_kept].T
    Xfull = pd.DataFrame(0.0, index=Xext.index, columns=gene_order)
    for g in g_in_idx:
        Xfull[g] = Xext[g]
    Xt = scaler.transform(Xfull.values)
    p = model.predict_proba(Xt)[:, 1]
    y_keep = np.array([y[i] for i, g in enumerate(gsms) if keep[i] and g in df.columns])
    if len(np.unique(y_keep)) < 2:
        return {"error": "single class after filter"}
    auc = roc_auc_score(y_keep, p)
    return {
        "n": int(len(samples_kept)),
        "n_refractory": int((y_keep == 1).sum()),
        "n_avid": int((y_keep == 0).sum()),
        "n_genes_matched": len(g_in_idx),
        "auc_raw": float(auc),
        "auc_direction_corrected": float(max(auc, 1 - auc)),
    }


def main() -> int:
    model, scaler, gene_order = train_tcga_logreg()
    summary = {}
    summary["R3A_gse76039"] = r3a_gse76039(model, scaler, gene_order)
    summary["R3B_gse151180"] = r3b_gse151180(model, scaler, gene_order)
    with open(OUT / "R3_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print("\n=== R3 SUMMARY ===")
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
