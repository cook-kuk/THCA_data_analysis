#!/usr/bin/env python3
"""v6 Wave 2 — scGPT zero-shot cell embedding + 8-target attention + DIAL on embedding.

Produces:
  $RES/foundation_models/scgpt_cell_embeddings.npy
  $RES/foundation_models/scgpt_gene_embeddings.npy
  $RES/foundation_models/scgpt_attention_8targets.tsv
  $RES/foundation_models/scgpt_cell_types.tsv
  $RES/foundation_models/scgpt_dial_analysis.md
  $RES/dial_scrna/scgpt_dial.tsv
"""
from __future__ import annotations
import argparse, sys, os, json, time
from pathlib import Path

DRUGGABLE = ["CYP1B1","TACSTD2","TMPRSS4","PLEKHA6","LDLR","GABRB2","B3GNT3","PTPRE"]
MARKERS = {
    "thyrocyte": ["TG","TPO","TSHR","PAX8","NKX2-1"],
    "immune":    ["CD3D","CD3E","CD8A","CD19","MS4A1","CD68","CD14"],
    "stromal":   ["COL1A1","COL3A1","ACTA2","PECAM1","VWF","FAP"],
}

def stamp(path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).touch()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="AnnData h5ad")
    ap.add_argument("--output_dir", default="/opt/thyroid-dash/project/results/v6_scrna")
    ap.add_argument("--batch_size", type=int, default=32)
    ap.add_argument("--device", default="auto")
    ap.add_argument("--model_repo", default="bowang-lab/scGPT-human")
    ap.add_argument("--model_cache", default="/data/thca/scrna/models/scgpt-human")
    args = ap.parse_args()

    out = Path(args.output_dir)
    (out/"foundation_models").mkdir(parents=True, exist_ok=True)
    (out/"dial_scrna").mkdir(parents=True, exist_ok=True)

    try:
        import torch, scanpy as sc, numpy as np, pandas as pd
        from huggingface_hub import snapshot_download
    except Exception as e:
        print(f"FATAL: core import: {e}", file=sys.stderr); sys.exit(2)
    try:
        import scgpt
        from scgpt.tasks import GeneEmbedding, CellEmbedding
    except Exception as e:
        print(f"FATAL: scGPT not importable. Install with: pip install git+https://github.com/bowang-lab/scGPT.git\n{e}",
              file=sys.stderr); sys.exit(2)

    device = (args.device if args.device != "auto"
              else ("cuda" if torch.cuda.is_available() else "cpu"))
    print(f"[scgpt] device={device}")

    print(f"[scgpt] downloading model {args.model_repo} → {args.model_cache}")
    snapshot_download(repo_id=args.model_repo, local_dir=args.model_cache, local_dir_use_symlinks=False)

    print(f"[scgpt] loading {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[scgpt] adata: {adata.shape}")

    # Ensure HVG and gene symbols present
    if "feature_name" in adata.var:
        adata.var_names = adata.var["feature_name"].astype(str)

    print("[scgpt] embedding cells …")
    t0 = time.time()
    cell_emb = scgpt.tasks.embed_data(adata, model_dir=args.model_cache,
                                      batch_size=args.batch_size, device=device)
    if hasattr(cell_emb, "obsm") and "X_scGPT" in cell_emb.obsm:
        cell_vec = cell_emb.obsm["X_scGPT"]
    else:
        cell_vec = np.asarray(cell_emb)
    np.save(out/"foundation_models/scgpt_cell_embeddings.npy", cell_vec)
    print(f"[scgpt] cell embeddings shape={cell_vec.shape} {time.time()-t0:.1f}s")

    # Gene embedding (best-effort; API may differ across scGPT versions)
    try:
        from scgpt.tasks import GeneEmbedding
        ge = GeneEmbedding(model_dir=args.model_cache, device=device)
        gene_vec = ge.embed(adata.var_names.tolist())
        np.save(out/"foundation_models/scgpt_gene_embeddings.npy", gene_vec)
        print(f"[scgpt] gene embeddings shape={gene_vec.shape}")
    except Exception as e:
        print(f"[scgpt] WARN gene embedding skipped: {e}")
        gene_vec = None

    # Attention for 8 druggable targets — top co-attended partners
    rows = []
    if gene_vec is not None:
        gene_names = list(adata.var_names)
        idx = {g:i for i,g in enumerate(gene_names)}
        for tgt in DRUGGABLE:
            if tgt not in idx:
                rows.append({"target":tgt, "rank":None, "partner":None, "cosine":None, "note":"target_missing_in_dataset"})
                continue
            v = gene_vec[idx[tgt]]
            sims = gene_vec @ v / ((np.linalg.norm(gene_vec, axis=1) * np.linalg.norm(v)) + 1e-9)
            sims[idx[tgt]] = -np.inf
            top = np.argsort(sims)[-10:][::-1]
            for r, j in enumerate(top, 1):
                rows.append({"target":tgt, "rank":r, "partner":gene_names[j], "cosine":float(sims[j]), "note":""})
    pd.DataFrame(rows).to_csv(out/"foundation_models/scgpt_attention_8targets.tsv", sep="\t", index=False)
    print("[scgpt] wrote attention table")

    # Zero-shot cell-type via marker-cosine
    ct_rows = []
    if gene_vec is not None:
        for ct, marks in MARKERS.items():
            mvec = np.mean([gene_vec[idx[m]] for m in marks if m in idx], axis=0)
            scores = cell_vec @ mvec / ((np.linalg.norm(cell_vec, axis=1) * np.linalg.norm(mvec)) + 1e-9)
            for i, s in enumerate(scores):
                ct_rows.append({"cell_id": str(adata.obs_names[i]), "cell_type": ct, "score": float(s)})
        df = pd.DataFrame(ct_rows)
        # call winner per cell
        winner = (df.sort_values("score", ascending=False)
                    .drop_duplicates("cell_id", keep="first")
                    [["cell_id","cell_type","score"]])
        winner.to_csv(out/"foundation_models/scgpt_cell_types.tsv", sep="\t", index=False)
    else:
        Path(out/"foundation_models/scgpt_cell_types.tsv").write_text("# cell-type call skipped (no gene embeddings)\n")

    # DIAL on scGPT embedding (LODO by patient)
    try:
        sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
        from v5p1_common import compute_dial, get_classifier_factories
        from sklearn.model_selection import LeaveOneGroupOut
        Y_col = "subtype" if "subtype" in adata.obs else ("mutation_status" if "mutation_status" in adata.obs else None)
        B_col = "patient_id" if "patient_id" in adata.obs else None
        if Y_col is None or B_col is None:
            print("[scgpt] WARN: cannot compute DIAL — missing Y or B column in obs", file=sys.stderr)
        else:
            Y = adata.obs[Y_col].astype(str).values
            B = adata.obs[B_col].astype(str).values
            logo = LeaveOneGroupOut()
            splits = list(logo.split(np.arange(len(Y)), Y, groups=B))
            facs = get_classifier_factories()
            rows = []
            for clf_name, fac in facs.items():
                res = compute_dial(cell_vec.astype("float32"), Y, B, fac, splits)
                res.update({"layer":"scgpt_cell_embedding","cancer":"THCA","classifier":clf_name})
                rows.append(res)
            pd.DataFrame(rows).to_csv(out/"dial_scrna/scgpt_dial.tsv", sep="\t", index=False)
    except Exception as e:
        print(f"[scgpt] WARN DIAL on embedding skipped: {e}")

    # Markdown summary
    md = ["# scGPT DIAL analysis\n",
          f"- model: `{args.model_repo}`",
          f"- input: `{args.input}` shape {adata.shape}",
          f"- device: {device}",
          f"- cell embeddings: shape {cell_vec.shape}",
          f"- DIAL output: `dial_scrna/scgpt_dial.tsv`"]
    (out/"foundation_models/scgpt_dial_analysis.md").write_text("\n".join(md))

    stamp(f"{args.output_dir}/.stamp_10_scgpt_embedding")
    print("[scgpt] DONE")

if __name__ == "__main__":
    main()
