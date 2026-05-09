"""Wave 4A LoRA fix — re-run after BatchNorm-on-last-batch crash.

Original Wave 4A `train_wave4a.py` LoRA path crashed at fold 1 with:
    ValueError: Expected more than 1 value per channel when training,
                got input size torch.Size([1, 256])

Root cause: a single-sample batch reached `nn.BatchNorm1d` in training
mode. The previous script had `drop_last=True` on the DataLoader and a
`if p_ids.size(0) < 2: continue` guard, but BN can still be flaky on
small training tails when sampler sizes change.

This fix replaces BatchNorm1d with LayerNorm (robust to any batch size,
including bs=1) in the head MLP. We also keep `drop_last=True` and the
batch>=2 guard as belt-and-suspenders.

Same data + same seeds + same eval pipeline as Wave 4A.
LoRA on ESM2-150M attention (q/k/v adapters, rank=8).

Eval (consistent with Wave 1):
  - In-domain 5-fold (within train pool)
  - Cross-source LOSO (4 sources)
  - ITSNdb_no_overlap (n=106, headline)
  - ITSNdb_in_master (n=213)
  - ITSNdb_combined (n=319)
  - VenusVaccine top10_mean

Outputs to --out_dir:
  predictions_lora.tsv
  lora_results.tsv     (long format with bootstrap CI95)
  LORA_FIX_REPORT.md
"""
from __future__ import annotations
import argparse, json, time, os
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, WeightedRandomSampler
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.model_selection import StratifiedKFold

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
HIDDEN = 640  # ESM2-150M
ESM2_NAME = "facebook/esm2_t30_150M_UR50D"


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
    if len(set(y)) < 2:
        return None, None, None
    auc = float(roc_auc_score(y, p))
    rng = np.random.default_rng(seed)
    n = len(y)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(set(y[idx])) < 2:
            continue
        boots.append(roc_auc_score(y[idx], p[idx]))
    if not boots:
        return auc, None, None
    boots = np.array(boots)
    return auc, float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))


def metrics_block(y, p, n_boot=1000):
    y = np.asarray(y); p = np.asarray(p)
    if len(set(y)) < 2:
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
        "Brier": float(brier_score_loss(y, p)),
        "ECE": ece(y, p),
        "n": int(len(y)),
        "n_pos": int(np.sum(y)),
    }


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------
def load_bundle_emb(bundle_path, emb_path):
    bundle = pd.read_csv(bundle_path, sep="\t")
    bundle["HLA_norm"] = bundle["HLA_norm"].fillna("")
    if "in_master" not in bundle.columns:
        bundle["in_master"] = False
    emb = torch.load(emb_path, map_location="cpu", weights_only=False)
    pep_keys = emb["pep_keys"]; hla_keys = emb["hla_keys"]
    pep_idx = {k: i for i, k in enumerate(pep_keys)}
    hla_idx = {k: i for i, k in enumerate(hla_keys)}
    venus_windows = emb["venus_windows"]
    def has_emb(p, h):
        return p in pep_idx and h in hla_idx
    nonvenus = bundle[~bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])].copy()
    nonvenus = nonvenus[[has_emb(p, h) for p, h in zip(nonvenus["peptide"], nonvenus["HLA_norm"])]].copy()
    return bundle, nonvenus, pep_idx, hla_idx, venus_windows


# ---------------------------------------------------------------------------
# LoRA head with LayerNorm (replaces BatchNorm to avoid bs=1 crash)
# ---------------------------------------------------------------------------
def focal_bce(logit, target, gamma=2.0, label_smooth=0.05):
    target = target * (1 - label_smooth) + 0.5 * label_smooth
    p = torch.sigmoid(logit)
    ce = F.binary_cross_entropy_with_logits(logit, target, reduction="none")
    pt = p * target + (1 - p) * (1 - target)
    return ((1 - pt) ** gamma * ce).mean(), ce


class LoraNeo(nn.Module):
    """ESM2 LoRA backbone + MLP head with LayerNorm (NOT BatchNorm)."""
    def __init__(self, base, hid=640, p_drop=0.3):
        super().__init__()
        self.base = base
        self.fc1 = nn.Linear(2 * hid, 256)
        self.ln1 = nn.LayerNorm(256)  # fix: was BatchNorm1d
        self.fc2 = nn.Linear(256, 256)
        self.ln2 = nn.LayerNorm(256)  # fix: was BatchNorm1d
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
        h = F.relu(self.ln1(self.fc1(x)))
        h = F.dropout(h, p=self.p, training=self.training)
        h = F.relu(self.ln2(self.fc2(h)))
        h = F.dropout(h, p=self.p, training=self.training)
        return self.head(h).squeeze(-1)


def train_lora(train_df, hla_pseudo_map, epochs=20, lr_lora=3e-4, lr_head=1e-3,
               batch=32, use_amp=True):
    from transformers import AutoTokenizer, AutoModel
    from peft import LoraConfig, get_peft_model
    tok = AutoTokenizer.from_pretrained(ESM2_NAME)
    base = AutoModel.from_pretrained(ESM2_NAME)
    lora_cfg = LoraConfig(
        r=8, lora_alpha=16,
        target_modules=["query", "key", "value"],
        lora_dropout=0.05, bias="none",
        task_type=None,
    )
    base = get_peft_model(base, lora_cfg)
    base.print_trainable_parameters()
    base = base.to(DEVICE)
    model = LoraNeo(base).to(DEVICE)

    # Two LR groups: LoRA params at 3e-4, head at 1e-3
    lora_params, head_params = [], []
    for n_, p_ in model.named_parameters():
        if not p_.requires_grad:
            continue
        if "base" in n_:
            lora_params.append(p_)
        else:
            head_params.append(p_)
    opt = torch.optim.AdamW(
        [
            {"params": lora_params, "lr": lr_lora},
            {"params": head_params, "lr": lr_head},
        ],
        weight_decay=1e-4,
    )
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    scaler = torch.cuda.amp.GradScaler(enabled=(use_amp and DEVICE == "cuda"))

    # Pre-tokenize
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

    ds = torch.utils.data.TensorDataset(
        pep_enc["input_ids"], pep_enc["attention_mask"],
        hla_enc["input_ids"], hla_enc["attention_mask"],
        y_all, src_codes,
    )
    sampler = WeightedRandomSampler(weights=w.tolist(), num_samples=n, replacement=True)
    # FIX: drop_last=True so we never see batch_size=1 at end of epoch.
    loader = DataLoader(ds, batch_size=batch, sampler=sampler, num_workers=0, drop_last=True)

    for ep in range(epochs):
        model.train()
        losses = []
        for p_ids, p_mask, h_ids, h_mask, y, _ in loader:
            if p_ids.size(0) < 2:
                continue  # extra safety net (LayerNorm doesn't need it but cheap)
            p_ids = p_ids.to(DEVICE); p_mask = p_mask.to(DEVICE)
            h_ids = h_ids.to(DEVICE); h_mask = h_mask.to(DEVICE)
            y = y.float().to(DEVICE)
            opt.zero_grad()
            if use_amp and DEVICE == "cuda":
                with torch.cuda.amp.autocast(dtype=torch.float16):
                    logit = model(p_ids, p_mask, h_ids, h_mask)
                    loss, _ = focal_bce(logit, y)
                scaler.scale(loss).backward()
                scaler.unscale_(opt)
                torch.nn.utils.clip_grad_norm_(
                    [p for p in model.parameters() if p.requires_grad], 1.0
                )
                scaler.step(opt)
                scaler.update()
            else:
                logit = model(p_ids, p_mask, h_ids, h_mask)
                loss, _ = focal_bce(logit, y)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(
                    [p for p in model.parameters() if p.requires_grad], 1.0
                )
                opt.step()
            losses.append(loss.item())
        sched.step()
        if (ep + 1) % 2 == 0 or ep == 0:
            print(f"    [lora-fix] ep{ep+1:>2d} loss={np.mean(losses):.4f}", flush=True)

    return model, tok


@torch.no_grad()
def predict_lora(model, tok, df, hla_pseudo_map, batch=64):
    model.eval()
    peps = df["peptide"].tolist()
    hlas = [hla_pseudo_map.get(h, h) for h in df["HLA_norm"].tolist()]
    pep_enc = tok(peps, padding=True, truncation=True, max_length=20, return_tensors="pt")
    hla_enc = tok(hlas, padding=True, truncation=True, max_length=64, return_tensors="pt")
    n = len(df)
    outs = []
    for i in range(0, n, batch):
        p_ids = pep_enc["input_ids"][i:i+batch].to(DEVICE)
        p_mask = pep_enc["attention_mask"][i:i+batch].to(DEVICE)
        h_ids = hla_enc["input_ids"][i:i+batch].to(DEVICE)
        h_mask = hla_enc["attention_mask"][i:i+batch].to(DEVICE)
        logit = model(p_ids, p_mask, h_ids, h_mask)
        outs.append(torch.sigmoid(logit).cpu())
    return torch.cat(outs).numpy()


# ---------------------------------------------------------------------------
# Eval pipeline
# ---------------------------------------------------------------------------
def run_full_eval(method_name, train_fn, predict_fn,
                  bundle, nonvenus, pep_idx, hla_idx, venus_windows,
                  hla_pseudo_map, out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    train_df = nonvenus[nonvenus["split"] == "train"].copy().reset_index(drop=True)

    summary = []
    all_pred_rows = []

    # In-domain 5-fold
    print(f"\n=== [{method_name}] In-domain 5-fold ===", flush=True)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    fold_aurocs = []
    for fold, (tr, te) in enumerate(skf.split(train_df, train_df["label"].values)):
        tr_df = train_df.iloc[tr].reset_index(drop=True)
        te_df = train_df.iloc[te].reset_index(drop=True)
        torch.manual_seed(0); np.random.seed(0)
        model_pack = train_fn(tr_df)
        p = predict_fn(model_pack, te_df)
        m = metrics_block(te_df["label"].values, p, n_boot=1000)
        m["fold"] = fold
        fold_aurocs.append(m["AUROC"])
        summary.append({"method": method_name, "testset": "in_domain_5fold",
                        "fold": fold, "in_master": "NA", **m})
        all_pred_rows.append(pd.DataFrame({
            "peptide": te_df["peptide"], "HLA_norm": te_df["HLA_norm"],
            "label": te_df["label"], "source": te_df["source"],
            "split": f"in_domain_fold{fold}", "pred_mean": p,
        }))
        print(f"  [{method_name}] fold{fold} AUROC={m['AUROC']:.3f} "
              f"[{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]", flush=True)
        del model_pack
        if DEVICE == "cuda":
            torch.cuda.empty_cache()
    print(f"  [{method_name}] 5-fold mean AUROC={np.mean(fold_aurocs):.3f}", flush=True)

    # LOSO
    print(f"\n=== [{method_name}] Cross-source LOSO ===", flush=True)
    for src in sorted(train_df["source"].unique()):
        tr_df = train_df[train_df["source"] != src].reset_index(drop=True)
        te_df = train_df[train_df["source"] == src].reset_index(drop=True)
        if te_df["label"].nunique() < 2:
            continue
        torch.manual_seed(0); np.random.seed(0)
        model_pack = train_fn(tr_df)
        p = predict_fn(model_pack, te_df)
        m = metrics_block(te_df["label"].values, p, n_boot=1000)
        m["held_out_source"] = src
        summary.append({"method": method_name, "testset": f"loso_{src}",
                        "fold": "NA", "in_master": "NA", **m})
        all_pred_rows.append(pd.DataFrame({
            "peptide": te_df["peptide"], "HLA_norm": te_df["HLA_norm"],
            "label": te_df["label"], "source": te_df["source"],
            "split": f"loso_{src}", "pred_mean": p,
        }))
        print(f"  [{method_name}] LOSO {src:>22s} AUROC={m['AUROC']:.3f} "
              f"[{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]", flush=True)
        del model_pack
        if DEVICE == "cuda":
            torch.cuda.empty_cache()

    # Final on full train pool
    print(f"\n=== [{method_name}] Final model on full train pool ===", flush=True)
    torch.manual_seed(0); np.random.seed(0)
    final_model = train_fn(train_df)

    # ITSNdb
    print(f"\n=== [{method_name}] ITSNdb external ===", flush=True)
    itsn = nonvenus[nonvenus["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy()
    p = predict_fn(final_model, itsn)
    itsn["pred_mean"] = p
    all_pred_rows.append(pd.DataFrame({
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
            summary.append({"method": method_name, "testset": sub_name,
                            "fold": "NA", "in_master": in_master_flag, **m})
            print(f"  [{method_name}] {sub_name:>20s}: n={m['n']} AUROC={m['AUROC']:.3f} "
                  f"[{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]", flush=True)

    # Venus
    print(f"\n=== [{method_name}] VenusVaccine ===", flush=True)
    venus_rows = bundle[bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])].copy()
    if len(venus_rows) > 0:
        train_pool = nonvenus[nonvenus["split"] == "train"]
        top_hlas = train_pool["HLA_norm"].value_counts().head(10).index.tolist()
        top_hlas = [h for h in top_hlas if h in hla_idx]
        big_rows = []
        for _, r in venus_rows.iterrows():
            pid = r["protein_id"]
            kmers = [k for k in venus_windows.get(pid, []) if k in pep_idx]
            for k in kmers:
                for h in top_hlas:
                    big_rows.append((pid, int(r["label"]), r["split"], k, h))
        big_df = pd.DataFrame(big_rows, columns=["protein_id", "label", "split", "peptide", "HLA_norm"])
        big_df["source"] = big_df["split"]
        if len(big_df) > 0:
            p = predict_fn(final_model, big_df)
            big_df["pred_mean"] = p
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
                    summary.append({"method": method_name,
                                    "testset": f"venus_{split_name}_top10",
                                    "fold": "NA", "in_master": "NA", **m})
                    print(f"  [{method_name}] {split_name:>18s} top10_mean: "
                          f"n={m['n']} AUROC={m['AUROC']:.3f} "
                          f"[{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]", flush=True)

    del final_model
    if DEVICE == "cuda":
        torch.cuda.empty_cache()

    pred_df = pd.concat(all_pred_rows, axis=0, ignore_index=True)
    pred_df.to_csv(out / f"predictions_{method_name}.tsv", sep="\t", index=False)

    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", type=str, default="/workspace/neo_bayes/bundle.tsv")
    ap.add_argument("--embeddings", type=str, default="/workspace/neo_bayes/embeddings.pt")
    ap.add_argument("--hla_pseudo", type=str, default="/workspace/neo_bayes/hla_pseudo.tsv")
    ap.add_argument("--out_dir", type=str, default=".")
    ap.add_argument("--epochs", type=int, default=20)
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    print(f"DEVICE={DEVICE}, epochs={args.epochs}", flush=True)
    print(f"out_dir={out.resolve()}", flush=True)

    bundle, nonvenus, pep_idx, hla_idx, venus_windows = load_bundle_emb(
        args.bundle, args.embeddings)
    hla_pseudo = pd.read_csv(args.hla_pseudo, sep="\t")
    hla_pseudo_map = dict(zip(hla_pseudo.iloc[:, 0], hla_pseudo.iloc[:, 1]))
    print(f"bundle={len(bundle)}, nonvenus={len(nonvenus)}, "
          f"pep_emb={len(pep_idx)}, hla_emb={len(hla_idx)}", flush=True)

    def train_fn_lora(train_df):
        return train_lora(train_df, hla_pseudo_map, epochs=args.epochs)

    def predict_fn_lora(model_tok, df):
        model, tok = model_tok
        valid = [h in hla_pseudo_map for h in df["HLA_norm"]]
        df_v = df[valid].copy().reset_index(drop=True)
        if len(df_v) == 0:
            return np.full(len(df), 0.5)
        p_v = predict_lora(model, tok, df_v, hla_pseudo_map)
        full = np.full(len(df), 0.5)
        v_idx = np.where(valid)[0]
        full[v_idx] = p_v
        return full

    t0 = time.time()
    summary = run_full_eval("lora", train_fn_lora, predict_fn_lora,
                            bundle, nonvenus, pep_idx, hla_idx,
                            venus_windows, hla_pseudo_map, out)
    elapsed = time.time() - t0
    print(f"\nTotal LoRA elapsed: {elapsed/60:.1f} min", flush=True)

    df = pd.DataFrame(summary)
    df.to_csv(out / "lora_results.tsv", sep="\t", index=False)
    print(f"Saved {len(df)} rows to {out / 'lora_results.tsv'}", flush=True)


if __name__ == "__main__":
    main()
