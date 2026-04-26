#!/usr/bin/env python3
"""v6 Wave 2 — Geneformer rank-based attention for BRAF-like cells.
Output: $RES/foundation_models/geneformer_ranking.tsv
"""
from __future__ import annotations
import argparse, sys, os, json
from pathlib import Path

DRUGGABLE = ["CYP1B1","TACSTD2","TMPRSS4","PLEKHA6","LDLR","GABRB2","B3GNT3","PTPRE"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", default="/opt/thyroid-dash/project/results/v6_scrna")
    ap.add_argument("--biomarkers_2773_tsv",
                    default="/opt/thyroid-dash/project/results/v8_statgen/v8_biomarker_DE_recomputation.tsv")
    args = ap.parse_args()
    out = Path(args.output_dir); (out/"foundation_models").mkdir(parents=True, exist_ok=True)

    try:
        import torch, scanpy as sc, numpy as np, pandas as pd
    except Exception as e:
        print(f"FATAL core import: {e}", file=sys.stderr); sys.exit(2)
    try:
        from transformers import AutoModel, AutoTokenizer
    except Exception as e:
        print(f"FATAL transformers: {e}", file=sys.stderr); sys.exit(2)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[geneformer] device={device}")

    print("[geneformer] loading model ctheodoris/Geneformer …")
    try:
        model = AutoModel.from_pretrained("ctheodoris/Geneformer", trust_remote_code=True).to(device)
        tok   = AutoTokenizer.from_pretrained("ctheodoris/Geneformer", trust_remote_code=True)
    except Exception as e:
        print(f"FATAL Geneformer load: {e}", file=sys.stderr); sys.exit(2)

    print(f"[geneformer] loading {args.input}")
    adata = sc.read_h5ad(args.input)
    if "feature_name" in adata.var:
        adata.var_names = adata.var["feature_name"].astype(str)

    # Filter to BRAF-like cells if metadata present
    sub_col = "subtype" if "subtype" in adata.obs else None
    if sub_col is not None:
        mask = adata.obs[sub_col].astype(str).str.contains("BRAF", case=False, na=False)
        adata_b = adata[mask].copy()
        print(f"[geneformer] BRAF-like cells: {adata_b.n_obs}")
    else:
        adata_b = adata
        print("[geneformer] no subtype column; using all cells")

    # Geneformer expects rank-tokenized input. The HF model card defines the tokenizer behavior.
    # For each cell: rank genes by expression, take top N (model max_seq_len),
    # tokenize to gene-ID embeddings, forward, capture attention to [CLS].
    n_cells = min(adata_b.n_obs, 5000)  # cap for compute
    sample_idx = np.random.RandomState(0).choice(adata_b.n_obs, n_cells, replace=False)
    expr = adata_b.X[sample_idx]
    if hasattr(expr, "toarray"): expr = expr.toarray()

    gene_names = np.array(adata_b.var_names)
    attn_sum = np.zeros(len(gene_names), dtype="float32")

    print(f"[geneformer] forward pass on {n_cells} cells …")
    import torch.nn.functional as F
    with torch.no_grad():
        for i in range(0, n_cells, 16):
            batch = expr[i:i+16]
            # rank genes, take top 2048 (Geneformer default)
            for c in batch:
                ranked = np.argsort(-c)[:2048]
                # build pseudo-token sequence using gene symbols
                tokens = tok([" ".join(gene_names[ranked].tolist())],
                             return_tensors="pt", truncation=True, padding=True,
                             max_length=2048).to(device)
                out_h = model(**tokens, output_attentions=True)
                # last layer mean attention from CLS to all tokens
                # shape: (batch, heads, q, k)
                last_attn = out_h.attentions[-1][0].mean(0)[0]  # CLS row
                attn_top = last_attn.detach().cpu().numpy()
                m = min(len(ranked), len(attn_top))
                attn_sum[ranked[:m]] += attn_top[:m]
            if (i//16) % 5 == 0:
                print(f"  [geneformer] cells {i+16}/{n_cells}")

    # Rank
    rank = np.argsort(-attn_sum)
    rows = []
    biomarkers_2773 = set()
    try:
        bdf = pd.read_csv(args.biomarkers_2773_tsv, sep="\t")
        if "gene" in bdf.columns and "new_significant" in bdf.columns:
            biomarkers_2773 = set(bdf[bdf["new_significant"]]["gene"].astype(str).tolist()[:2773])
    except Exception:
        pass
    for r, gi in enumerate(rank, 1):
        g = gene_names[gi]
        rows.append({"gene": g, "attention_rank": r,
                     "attention_mean": float(attn_sum[gi] / n_cells),
                     "in_8_druggable": g in DRUGGABLE,
                     "in_2773_biomarker": g in biomarkers_2773})
    pd.DataFrame(rows).to_csv(out/"foundation_models/geneformer_ranking.tsv", sep="\t", index=False)

    # Short summary
    summary = ["# Geneformer attention ranking\n",
               f"- cells used: {n_cells} (capped)",
               f"- top-100 contains druggable: " +
                 ", ".join([r["gene"] for r in rows[:100] if r["in_8_druggable"]]),
               f"- top-100 ∩ 2773 biomarker: " +
                 str(sum(1 for r in rows[:100] if r["in_2773_biomarker"])) ]
    (out/"foundation_models/geneformer_summary.md").write_text("\n".join(summary))
    Path(f"{args.output_dir}/.stamp_20_geneformer").touch()
    print("[geneformer] DONE")

if __name__ == "__main__":
    main()
