"""M5 — kNN retrieval over ESM2 peptide embeddings (Khandelwal 2020).

For each ITSNdb test peptide:
  1. Find top-K nearest train peptides by cosine similarity (within same HLA bucket
     when available; back off to global pool otherwise).
  2. Score = sim-weighted vote of neighbor labels.

Sweep K ∈ {5, 10, 20, 50}; pick best on a held-out fold of train pool, then
report on each ITSNdb subset (main / val / combined / no_overlap).
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parent
WAVE1 = ROOT.parent
RNG = np.random.default_rng(0)

EMB_LOCAL = ROOT / "embeddings_local.pt"
BUNDLE = WAVE1 / "bundle.tsv"
OUT_PRED = ROOT / "predictions_knn.tsv"
OUT_LOG = ROOT / "m5_knn.log"


def cosine_topk(query_vec, ref_mat, k):
    """Return top-k indices of ref_mat by cosine sim with query_vec.
    query_vec: [D], ref_mat: [N, D]. Returns idx [k], sims [k]."""
    q = query_vec / (np.linalg.norm(query_vec) + 1e-9)
    r = ref_mat / (np.linalg.norm(ref_mat, axis=1, keepdims=True) + 1e-9)
    sims = r @ q
    idx = np.argpartition(-sims, kth=min(k, len(sims) - 1))[:k]
    order = np.argsort(-sims[idx])
    return idx[order], sims[idx[order]]


def knn_score(test_emb, ref_emb, ref_labels, k):
    """Vectorized cosine kNN. Returns [n_test] sim-weighted vote score in [0,1]."""
    q = test_emb / (np.linalg.norm(test_emb, axis=1, keepdims=True) + 1e-9)
    r = ref_emb / (np.linalg.norm(ref_emb, axis=1, keepdims=True) + 1e-9)
    sims = q @ r.T  # [n_test, n_ref]
    # top-k per row
    topk_idx = np.argpartition(-sims, kth=min(k, sims.shape[1] - 1), axis=1)[:, :k]
    rows = np.arange(sims.shape[0])[:, None]
    topk_sims = sims[rows, topk_idx]  # [n_test, k]
    topk_lbl = ref_labels[topk_idx]   # [n_test, k]
    # softmax-tempered weighting (more numerically stable than raw sim if sims negative)
    w = np.maximum(topk_sims, 0.0)
    denom = w.sum(axis=1) + 1e-9
    score = (w * topk_lbl).sum(axis=1) / denom
    return score


def evaluate_subset(df, score, label_col="label"):
    y = df[label_col].values.astype(int)
    if len(set(y)) < 2 or len(y) < 5:
        return {"AUROC": float("nan"), "n": len(y), "n_pos": int(y.sum())}
    return {"AUROC": float(roc_auc_score(y, score)),
            "n": int(len(y)), "n_pos": int(y.sum())}


def main():
    log = open(OUT_LOG, "w")
    def L(*a):
        s = " ".join(str(x) for x in a)
        print(s); log.write(s + "\n"); log.flush()

    L("loading embeddings + bundle…")
    pkg = torch.load(EMB_LOCAL, map_location="cpu", weights_only=False)
    pep_keys = pkg["pep_keys"]
    pep_emb = pkg["pep_emb"].numpy()
    pep_idx = {k: i for i, k in enumerate(pep_keys)}
    bundle = pd.read_csv(BUNDLE, sep="\t")
    bundle["HLA_norm"] = bundle["HLA_norm"].fillna("")
    L(f"  embeddings: {pep_emb.shape}, bundle: {len(bundle)}")

    # filter rows that have a peptide embedding
    has_emb = bundle["peptide"].map(lambda p: p in pep_idx)
    bundle = bundle[has_emb].copy()
    L(f"  bundle rows with embedding: {len(bundle)}")

    train_df = bundle[bundle["split"] == "train"].reset_index(drop=True).copy()
    itsn_df = bundle[bundle["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].reset_index(drop=True).copy()
    L(f"  train: {len(train_df)} pos_rate={train_df['label'].mean():.3f}")
    L(f"  itsndb: {len(itsn_df)} pos_rate={itsn_df['label'].mean():.3f}")

    # Per-HLA train index
    train_emb = pep_emb[np.array([pep_idx[p] for p in train_df["peptide"].tolist()])]
    train_lbl = train_df["label"].values.astype(int)
    train_hla = train_df["HLA_norm"].values
    by_hla = {}
    for h in np.unique(train_hla):
        by_hla[h] = np.where(train_hla == h)[0]

    # ----- K sweep on train held-out fold (5-fold) -------------------------------
    L("\n=== K sweep on train held-out fold ===")
    Ks = [5, 10, 20, 50]
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    sweep = {k: [] for k in Ks}
    tr_idx, te_idx = next(iter(skf.split(train_df, train_df["label"].values)))
    fold_tr = train_df.iloc[tr_idx].reset_index(drop=True)
    fold_te = train_df.iloc[te_idx].reset_index(drop=True)
    fold_tr_emb = train_emb[tr_idx]
    fold_tr_lbl = train_lbl[tr_idx]
    fold_tr_hla = train_hla[tr_idx]
    fold_te_emb = train_emb[te_idx]
    fold_te_hla = train_hla[te_idx]
    fold_by_hla = {}
    for h in np.unique(fold_tr_hla):
        fold_by_hla[h] = np.where(fold_tr_hla == h)[0]

    for K in Ks:
        scores = np.zeros(len(fold_te))
        for i in range(len(fold_te)):
            h = fold_te_hla[i]
            if h in fold_by_hla and len(fold_by_hla[h]) >= K:
                idxs = fold_by_hla[h]
            else:
                # Fallback to global
                idxs = np.arange(len(fold_tr_emb))
            ref = fold_tr_emb[idxs]
            ref_lbl = fold_tr_lbl[idxs]
            score = knn_score(fold_te_emb[i:i+1], ref, ref_lbl, k=min(K, len(ref)))
            scores[i] = score[0]
        try:
            auc = roc_auc_score(fold_te["label"].values, scores)
        except ValueError:
            auc = float("nan")
        sweep[K] = float(auc)
        L(f"  K={K:>3d}  fold-val AUROC={auc:.4f}")

    best_K = max(sweep, key=lambda k: (sweep[k] if not np.isnan(sweep[k]) else -1))
    L(f"\n  best K (by held-out train AUROC): {best_K}  (auc={sweep[best_K]:.4f})")

    # ----- ITSNdb evaluation with best K -----------------------------------------
    L(f"\n=== ITSNdb scoring with K={best_K} (per-HLA bucket → fallback global) ===")
    itsn_emb = pep_emb[np.array([pep_idx[p] for p in itsn_df["peptide"].tolist()])]
    itsn_hla = itsn_df["HLA_norm"].values

    K_eval = best_K
    scores_eval = np.zeros(len(itsn_df))
    n_per_hla = []
    for i in range(len(itsn_df)):
        h = itsn_hla[i]
        if h in by_hla and len(by_hla[h]) >= K_eval:
            idxs = by_hla[h]
        else:
            idxs = np.arange(len(train_emb))
        ref = train_emb[idxs]
        ref_lbl = train_lbl[idxs]
        s = knn_score(itsn_emb[i:i+1], ref, ref_lbl, k=min(K_eval, len(ref)))
        scores_eval[i] = s[0]
        n_per_hla.append(len(idxs))

    itsn_df["pred_knn"] = scores_eval
    itsn_df["n_ref_pool"] = n_per_hla
    itsn_df["K"] = K_eval
    itsn_df.to_csv(OUT_PRED, sep="\t", index=False)
    L(f"  saved → {OUT_PRED}  ({len(itsn_df)} rows)")

    # subsets
    L("\n--- per-subset AUROC ---")
    in_master = itsn_df["in_master"].astype(bool).values
    for name, mask in [
        ("ITSNdb_main", itsn_df["split"] == "ext_itsndb_main"),
        ("ITSNdb_Val",  itsn_df["split"] == "ext_itsndb_val"),
        ("ITSNdb_combined", np.ones(len(itsn_df), dtype=bool)),
        ("ITSNdb_in_master", in_master),
        ("ITSNdb_no_overlap", ~in_master),
    ]:
        sub = itsn_df[mask]
        m = evaluate_subset(sub, sub["pred_knn"].values)
        L(f"  {name:>22s}: n={m['n']:>4d} pos={m['n_pos']:>3d} AUROC={m['AUROC']:.4f}")

    # save sweep summary
    with open(ROOT / "m5_knn_sweep.json", "w") as f:
        json.dump({"sweep": sweep, "best_K": best_K}, f, indent=2)
    log.close()
    print("\nM5 done.")


if __name__ == "__main__":
    main()
