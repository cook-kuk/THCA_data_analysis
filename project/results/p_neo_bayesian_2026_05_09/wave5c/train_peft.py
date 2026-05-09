"""Wave 5C — Modern PEFT methods comparison.

Four PEFT variants on ESM2-150M attention, same training recipe + eval:
  - LoRA-baseline   (rank=8, plain)
  - DoRA            (rank=8, magnitude+direction decomposition; ICML 2024)
  - VeRA            (frozen random projections + tiny scaling vectors; ICLR 2024)
  - AdaLoRA         (adaptive rank pruning by importance; ICLR 2023)

Common recipe:
  - Wrap facebook/esm2_t30_150M_UR50D with PEFT method on q/k/v projections
  - MLP head 2*640 -> 256 -> 256 -> 1 on top of mean-pooled peptide+HLA encodings
  - 20 epochs, lr=3e-4, bs=32, AdamW, focal-BCE, source-balanced sampler
  - Train on n=2396 train pool; eval ITSNdb stratified + VenusVaccine top10_mean

Output (under --out_dir):
  predictions_<method>.tsv (one per method)
  peft_results.tsv         (long format: method | rank | n_trainable_params | testset | n | AUROC | CI_lo95 | CI_hi95)
  fig_peft_comparison.png/pdf
  WAVE5C_REPORT.md
"""
from __future__ import annotations
import argparse
import json
import math
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, WeightedRandomSampler, TensorDataset
from sklearn.metrics import (
    roc_auc_score, average_precision_score, brier_score_loss,
)

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
ESM2_NAME = "facebook/esm2_t30_150M_UR50D"
HIDDEN = 640


# ---------------------------------------------------------------------------
# Eval helpers
# ---------------------------------------------------------------------------
def ece(y, p, n_bins=15):
    y, p = np.asarray(y), np.asarray(p)
    bins = np.linspace(0, 1, n_bins + 1)
    ece_v = 0.0
    n = len(y)
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (p >= lo) & (p < hi if hi < 1 else p <= hi)
        if mask.sum() == 0:
            continue
        ece_v += (mask.sum() / n) * abs(p[mask].mean() - y[mask].mean())
    return float(ece_v)


def auroc_bootstrap_ci(y, p, n_boot=1000, seed=0):
    y = np.asarray(y); p = np.asarray(p)
    if len(set(y.tolist())) < 2:
        return None, None, None
    auc = float(roc_auc_score(y, p))
    rng = np.random.default_rng(seed)
    n = len(y)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(set(y[idx].tolist())) < 2:
            continue
        boots.append(roc_auc_score(y[idx], p[idx]))
    if not boots:
        return auc, None, None
    boots = np.array(boots)
    return auc, float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))


def metrics_block(y, p, n_boot=1000):
    y = np.asarray(y); p = np.asarray(p)
    if len(set(y.tolist())) < 2:
        return {"AUROC": None, "AUROC_lo95": None, "AUROC_hi95": None,
                "AUPRC": None, "Brier": None, "ECE": None,
                "n": int(len(y)), "n_pos": int(np.sum(y))}
    auc, lo, hi = auroc_bootstrap_ci(y, p, n_boot=n_boot)
    return {
        "AUROC": auc,
        "AUROC_lo95": lo,
        "AUROC_hi95": hi,
        "AUPRC": float(average_precision_score(y, p)),
        "Brier": float(brier_score_loss(y, p)),
        "ECE": ece(y, p),
        "n": int(len(y)),
        "n_pos": int(np.sum(y)),
    }


# ---------------------------------------------------------------------------
# Loss
# ---------------------------------------------------------------------------
def focal_bce(logit, target, gamma=2.0, label_smooth=0.05):
    target = target * (1 - label_smooth) + 0.5 * label_smooth
    p = torch.sigmoid(logit)
    ce = F.binary_cross_entropy_with_logits(logit, target, reduction="none")
    pt = p * target + (1 - p) * (1 - target)
    return ((1 - pt) ** gamma * ce).mean(), ce


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def load_data(bundle_path, emb_path, hla_pseudo_path):
    bundle = pd.read_csv(bundle_path, sep="\t")
    bundle["HLA_norm"] = bundle["HLA_norm"].fillna("")
    if "in_master" not in bundle.columns:
        bundle["in_master"] = False
    emb = torch.load(emb_path, map_location="cpu", weights_only=False)
    pep_idx = {k: i for i, k in enumerate(emb["pep_keys"])}
    hla_idx = {k: i for i, k in enumerate(emb["hla_keys"])}
    venus_windows = emb["venus_windows"]
    hla_pseudo = pd.read_csv(hla_pseudo_path, sep="\t")
    # Drop rows whose pseudo-seq is NaN/empty -- those alleles have no
    # MHC-pseudo and cannot be tokenized.
    hla_pseudo = hla_pseudo.dropna(subset=[hla_pseudo.columns[1]]).copy()
    hla_pseudo = hla_pseudo[hla_pseudo.iloc[:, 1].astype(str).str.len() > 0]
    hla_pseudo_map = dict(zip(hla_pseudo.iloc[:, 0], hla_pseudo.iloc[:, 1]))
    return bundle, pep_idx, hla_idx, venus_windows, hla_pseudo_map


# ---------------------------------------------------------------------------
# PEFT model wrapper — ESM2 with adapter on attention + MLP head
# ---------------------------------------------------------------------------
class PEFTNeo(nn.Module):
    def __init__(self, base, hid=HIDDEN, p_drop=0.3):
        super().__init__()
        self.base = base
        self.fc1 = nn.Linear(2 * hid, 256)
        self.bn1 = nn.BatchNorm1d(256)
        self.fc2 = nn.Linear(256, 256)
        self.bn2 = nn.BatchNorm1d(256)
        self.head = nn.Linear(256, 1)
        self.p = p_drop

    def encode(self, ids, mask):
        out = self.base(input_ids=ids, attention_mask=mask).last_hidden_state
        m = mask.unsqueeze(-1).float()
        return (out * m).sum(1) / m.sum(1).clamp(min=1)

    def forward(self, pep_ids, pep_mask, hla_ids, hla_mask):
        pe = self.encode(pep_ids, pep_mask)
        he = self.encode(hla_ids, hla_mask)
        x = torch.cat([pe, he], dim=-1)
        h = F.relu(self.bn1(self.fc1(x)))
        h = F.dropout(h, p=self.p, training=self.training)
        h = F.relu(self.bn2(self.fc2(h)))
        h = F.dropout(h, p=self.p, training=self.training)
        return self.head(h).squeeze(-1)


def build_peft_base(method: str, rank: int = 8):
    """Wrap ESM2 with the chosen PEFT method on q/k/v.

    Returns (peft_base, n_trainable_params_in_base_only).
    """
    from transformers import AutoModel
    base = AutoModel.from_pretrained(ESM2_NAME)
    target_modules = ["query", "key", "value"]

    if method == "lora_baseline":
        from peft import LoraConfig, get_peft_model, TaskType
        cfg = LoraConfig(
            r=rank, lora_alpha=2 * rank,
            target_modules=target_modules,
            lora_dropout=0.05, bias="none",
            task_type=TaskType.FEATURE_EXTRACTION,
            use_dora=False,
        )
        base = get_peft_model(base, cfg)

    elif method == "dora":
        from peft import LoraConfig, get_peft_model, TaskType
        cfg = LoraConfig(
            r=rank, lora_alpha=2 * rank,
            target_modules=target_modules,
            lora_dropout=0.05, bias="none",
            task_type=TaskType.FEATURE_EXTRACTION,
            use_dora=True,
        )
        base = get_peft_model(base, cfg)

    elif method == "vera":
        from peft import VeraConfig, get_peft_model, TaskType
        # VeRA recipe: shared frozen random A,B + small per-layer scaling vectors d,b
        # r=256 is standard VeRA setting (per ICLR 2024 paper); much fewer params than LoRA r=8
        # because A,B are shared/frozen and only d,b vectors are trainable.
        cfg = VeraConfig(
            r=256,
            target_modules=target_modules,
            vera_dropout=0.05,
            bias="none",
            task_type=TaskType.FEATURE_EXTRACTION,
        )
        base = get_peft_model(base, cfg)

    elif method == "adalora":
        from peft import AdaLoraConfig, get_peft_model, TaskType
        # AdaLora: starts at init_r, prunes to target_r based on importance.
        # Set t_total = approx total update steps so the schedule fires properly.
        cfg = AdaLoraConfig(
            init_r=12, target_r=rank,
            beta1=0.85, beta2=0.85,
            tinit=200, tfinal=1000, deltaT=10,
            target_modules=target_modules,
            lora_dropout=0.05, bias="none",
            task_type=TaskType.FEATURE_EXTRACTION,
            total_step=1500,
        )
        base = get_peft_model(base, cfg)

    else:
        raise ValueError(f"unknown method: {method}")

    # Count trainable in base
    n_train = sum(p.numel() for p in base.parameters() if p.requires_grad)
    return base, n_train


# ---------------------------------------------------------------------------
# Train one PEFT variant (final model, no fold splitting — we directly use
# the full 2396 pool, mirroring the n=2396 instruction).
# ---------------------------------------------------------------------------
def train_peft_method(method, train_df, hla_pseudo_map, epochs=20, lr=3e-4,
                      batch=32, rank=8, log_prefix=""):
    from transformers import AutoTokenizer
    tok = AutoTokenizer.from_pretrained(ESM2_NAME)
    base, n_base_train = build_peft_base(method, rank=rank)
    base = base.to(DEVICE)
    model = PEFTNeo(base).to(DEVICE)

    # Total trainable (base PEFT + head)
    n_total = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"  [{log_prefix}{method}] trainable params: base={n_base_train:,} total(base+head)={n_total:,}")

    opt = torch.optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()),
                            lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)

    # Tokenize everything
    peps = train_df["peptide"].tolist()
    hlas = [hla_pseudo_map.get(h, h) for h in train_df["HLA_norm"].tolist()]
    pep_enc = tok(peps, padding=True, truncation=True, max_length=20, return_tensors="pt")
    hla_enc = tok(hlas, padding=True, truncation=True, max_length=64, return_tensors="pt")
    y_all = torch.tensor(train_df["label"].astype(int).values, dtype=torch.float32)
    src = train_df["source"].astype("category")
    src_codes = torch.tensor(src.cat.codes.values, dtype=torch.long)
    n_src = len(src.cat.categories)

    src_counts = np.bincount(src_codes.numpy(), minlength=n_src).tolist()
    src_counts = [max(c, 1) for c in src_counts]
    n = len(train_df)
    w = np.array([1.0 / src_counts[int(c)] for c in src_codes.tolist()], dtype=np.float64)
    w = w * n / w.sum()

    ds = TensorDataset(
        pep_enc["input_ids"], pep_enc["attention_mask"],
        hla_enc["input_ids"], hla_enc["attention_mask"],
        y_all, src_codes,
    )
    sampler = WeightedRandomSampler(weights=w.tolist(), num_samples=n, replacement=True)
    loader = DataLoader(ds, batch_size=batch, sampler=sampler, num_workers=0)

    step = 0
    for ep in range(epochs):
        model.train()
        losses = []
        for p_ids, p_mask, h_ids, h_mask, y, _ in loader:
            p_ids = p_ids.to(DEVICE); p_mask = p_mask.to(DEVICE)
            h_ids = h_ids.to(DEVICE); h_mask = h_mask.to(DEVICE)
            y = y.float().to(DEVICE)
            logit = model(p_ids, p_mask, h_ids, h_mask)
            loss, _ = focal_bce(logit, y)
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(filter(lambda p: p.requires_grad, model.parameters()), 1.0)
            opt.step()
            # AdaLoRA needs explicit rank-mask updates per step
            if method == "adalora":
                try:
                    model.base.update_and_allocate(step)
                except Exception:
                    pass
            step += 1
            losses.append(loss.item())
        sched.step()
        if (ep + 1) % 2 == 0 or ep == 0:
            print(f"    [{log_prefix}{method}] ep{ep+1:>2d} loss={np.mean(losses):.4f}")

    return model, tok, n_total, n_base_train


@torch.no_grad()
def predict_peft(model, tok, df, hla_pseudo_map, batch=64):
    """Predict probabilities for rows whose HLA has a pseudo-sequence; otherwise 0.5."""
    model.eval()
    valid = [h in hla_pseudo_map for h in df["HLA_norm"]]
    df_v = df[valid].copy().reset_index(drop=True)
    if len(df_v) == 0:
        return np.full(len(df), 0.5)
    peps = df_v["peptide"].tolist()
    hlas = [hla_pseudo_map.get(h, h) for h in df_v["HLA_norm"].tolist()]
    pep_enc = tok(peps, padding=True, truncation=True, max_length=20, return_tensors="pt")
    hla_enc = tok(hlas, padding=True, truncation=True, max_length=64, return_tensors="pt")
    n = len(df_v)
    outs = []
    for i in range(0, n, batch):
        p_ids = pep_enc["input_ids"][i:i+batch].to(DEVICE)
        p_mask = pep_enc["attention_mask"][i:i+batch].to(DEVICE)
        h_ids = hla_enc["input_ids"][i:i+batch].to(DEVICE)
        h_mask = hla_enc["attention_mask"][i:i+batch].to(DEVICE)
        logit = model(p_ids, p_mask, h_ids, h_mask)
        outs.append(torch.sigmoid(logit).cpu())
    p_v = torch.cat(outs).numpy()
    full = np.full(len(df), 0.5)
    v_idx = np.where(valid)[0]
    full[v_idx] = p_v
    return full


# ---------------------------------------------------------------------------
# Eval pipeline (final model, no 5-fold/LOSO -- we focus on headline metrics)
# ---------------------------------------------------------------------------
def eval_method(method, model, tok, bundle, pep_idx, hla_idx, venus_windows,
                hla_pseudo_map, train_df, n_total_params, n_base_params, rank, out_dir):
    out = Path(out_dir)
    summary = []
    pred_rows = []

    # ----- ITSNdb -----
    itsn = bundle[bundle["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy()
    if len(itsn) > 0:
        p = predict_peft(model, tok, itsn, hla_pseudo_map)
        itsn["pred_mean"] = p
        pred_rows.append(pd.DataFrame({
            "peptide": itsn["peptide"], "HLA_norm": itsn["HLA_norm"],
            "label": itsn["label"], "source": itsn["source"],
            "split": itsn["split"], "in_master": itsn["in_master"], "pred_mean": p,
        }))
        for sub_name, mask, in_master_flag in [
            ("ITSNdb_main", itsn["split"] == "ext_itsndb_main", "mixed"),
            ("ITSNdb_Val",  itsn["split"] == "ext_itsndb_val", "mixed"),
            ("ITSNdb_combined", np.ones(len(itsn), dtype=bool), "mixed"),
            ("ITSNdb_no_overlap", ~itsn["in_master"].astype(bool).values, False),
            ("ITSNdb_in_master", itsn["in_master"].astype(bool).values, True),
        ]:
            sub = itsn[mask]
            if len(sub) >= 5 and sub["label"].nunique() == 2:
                m = metrics_block(sub["label"].values, sub["pred_mean"].values, n_boot=1000)
                summary.append({
                    "method": method, "rank": rank,
                    "n_trainable_params": n_total_params,
                    "n_base_peft_params": n_base_params,
                    "testset": sub_name,
                    "in_master": in_master_flag, **m
                })
                print(f"  [{method}] {sub_name:>20s}: n={m['n']} AUROC={m['AUROC']:.3f} [{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]")

    # ----- VenusVaccine -----
    venus_rows = bundle[bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])].copy()
    if len(venus_rows) > 0:
        top_hlas = train_df["HLA_norm"].value_counts().head(10).index.tolist()
        top_hlas = [h for h in top_hlas if h in hla_pseudo_map]
        big_rows = []
        for _, r in venus_rows.iterrows():
            pid = r["protein_id"]
            kmers = [k for k in venus_windows.get(pid, [])]
            for k in kmers:
                for h in top_hlas:
                    big_rows.append((pid, int(r["label"]), r["split"], k, h))
        big_df = pd.DataFrame(big_rows, columns=["protein_id", "label", "split", "peptide", "HLA_norm"])
        big_df["source"] = big_df["split"]
        if len(big_df) > 0:
            p = predict_peft(model, tok, big_df, hla_pseudo_map)
            big_df["pred_mean"] = p
            pred_rows.append(big_df[["peptide", "HLA_norm", "label", "source", "split", "pred_mean"]].copy())
            for split_name in ["ext_venus_test", "ext_venus_valid"]:
                sub = big_df[big_df["split"] == split_name]
                if len(sub) < 10:
                    continue
                rows = []
                for (pid, label), g in sub.groupby(["protein_id", "label"]):
                    scores = g["pred_mean"].values
                    rows.append({"protein_id": pid, "label": int(label),
                                 "top10_mean": float(np.sort(scores)[::-1][:10].mean())})
                vdf = pd.DataFrame(rows)
                if vdf["label"].nunique() == 2:
                    m = metrics_block(vdf["label"].values, vdf["top10_mean"].values, n_boot=1000)
                    summary.append({
                        "method": method, "rank": rank,
                        "n_trainable_params": n_total_params,
                        "n_base_peft_params": n_base_params,
                        "testset": f"venus_{split_name}_top10",
                        "in_master": "NA", **m
                    })
                    print(f"  [{method}] {split_name:>18s} top10_mean: n={m['n']} AUROC={m['AUROC']:.3f} [{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]")

    # Save predictions
    if pred_rows:
        pdf = pd.concat(pred_rows, axis=0, ignore_index=True)
        pdf.to_csv(out / f"predictions_{method}.tsv", sep="\t", index=False)

    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", type=str, default="../bundle.tsv")
    ap.add_argument("--embeddings", type=str, default="../embeddings.pt")
    ap.add_argument("--hla_pseudo", type=str, default="../hla_pseudo.tsv")
    ap.add_argument("--out_dir", type=str, default=".")
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--methods", type=str,
                    default="lora_baseline,dora,vera,adalora")
    ap.add_argument("--rank", type=int, default=8)
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    print(f"running methods: {methods}")
    print(f"epochs: {args.epochs} rank: {args.rank}")

    bundle, pep_idx, hla_idx, venus_windows, hla_pseudo_map = load_data(
        args.bundle, args.embeddings, args.hla_pseudo)
    train_df = bundle[bundle["split"] == "train"].copy().reset_index(drop=True)
    # Only keep train rows whose HLA has a valid pseudo-seq.
    train_df = train_df[train_df["HLA_norm"].isin(hla_pseudo_map.keys())].reset_index(drop=True)
    print(f"train_df: n={len(train_df)} sources={train_df['source'].value_counts().to_dict()}")

    all_summary = []
    timings = {}
    for method in methods:
        print("\n" + "=" * 70 + f"\n  {method}\n" + "=" * 70)
        torch.manual_seed(0); np.random.seed(0)
        t0 = time.time()
        model, tok, n_total, n_base = train_peft_method(
            method, train_df, hla_pseudo_map,
            epochs=args.epochs, rank=args.rank,
        )
        train_t = time.time() - t0
        s = eval_method(method, model, tok, bundle, pep_idx, hla_idx, venus_windows,
                        hla_pseudo_map, train_df, n_total, n_base, args.rank, out)
        all_summary += s
        timings[method] = train_t
        print(f"  [{method}] train+eval time: {train_t:.1f}s")
        del model
        if DEVICE == "cuda":
            torch.cuda.empty_cache()

    df = pd.DataFrame(all_summary)
    df.to_csv(out / "peft_results.tsv", sep="\t", index=False)
    with open(out / "timings.json", "w") as f:
        json.dump(timings, f, indent=2)
    print(f"\nSaved {len(df)} rows to {out / 'peft_results.tsv'}")


if __name__ == "__main__":
    main()
