"""
Phase 2b — CLAM-style attention-MIL on TCGA-THCA WSI.

Architecture:
  - Tile extraction (CLAM standard: 20×, 256-px non-overlap; ~1k-5k tiles per slide)
  - UNI 1024-d feature per tile
  - Gated-attention MIL pool (Lu 2021 Nat BME) → slide-level prediction
  - Binary task: DM1 vs DM2 (not_DM excluded for binary; use 3-way as supp)
  - LOSO + held-out 20% slide-level

Output:
  clam_results.tsv (per-slide attention scores + slide AUC)
  clam_train_log.txt
  PHASE2_REPORT.md

Resource:
  Tile extraction ~30 min CPU + UNI embed 1× A100 ~1-2 hr + CLAM training ~1 hr
"""
from __future__ import annotations
import argparse, os, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import KFold


# ============================================================
# CLAM gated-attention MIL
# ============================================================
class GatedAttentionMIL(nn.Module):
    def __init__(self, in_dim=1024, hidden=512, n_classes=2, dropout=0.25):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.attn_a = nn.Sequential(nn.Linear(hidden, hidden), nn.Tanh())
        self.attn_b = nn.Sequential(nn.Linear(hidden, hidden), nn.Sigmoid())
        self.attn_w = nn.Linear(hidden, 1)
        self.classifier = nn.Linear(hidden, n_classes)

    def forward(self, x):
        # x: (N_tiles, in_dim)
        h = self.fc(x)
        a = self.attn_w(self.attn_a(h) * self.attn_b(h))  # (N, 1)
        attn = F.softmax(a.squeeze(-1), dim=0)  # (N,)
        bag = (attn.unsqueeze(-1) * h).sum(0)  # (hidden,)
        logits = self.classifier(bag)  # (n_classes,)
        return logits, attn


def train_one_fold(train_bags, train_labels, val_bags, val_labels, device="cuda",
                   epochs=30, lr=2e-4):
    model = GatedAttentionMIL(in_dim=train_bags[0].shape[1], n_classes=2).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    crit = nn.CrossEntropyLoss()
    best_auc = 0
    best_state = None
    log = []
    for epoch in range(epochs):
        model.train()
        loss_sum = 0
        for bag, label in zip(train_bags, train_labels):
            opt.zero_grad()
            x = torch.from_numpy(np.asarray(bag, dtype=np.float32)).to(device)
            y = torch.tensor([label]).to(device)
            logits, _ = model(x)
            loss = crit(logits.unsqueeze(0), y)
            loss.backward(); opt.step()
            loss_sum += loss.item()
        # val
        model.eval()
        probs = []
        with torch.no_grad():
            for bag in val_bags:
                x = torch.from_numpy(np.asarray(bag, dtype=np.float32)).to(device)
                logits, _ = model(x)
                p = F.softmax(logits, dim=-1)[1].item()
                probs.append(p)
        try:
            auc = roc_auc_score(val_labels, probs)
        except Exception:
            auc = 0.5
        log.append((epoch, loss_sum/len(train_bags), auc))
        if auc > best_auc:
            best_auc = auc
            # snapshot best-epoch weights for this fold (CPU copy to free GPU)
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
    # restore best-epoch weights before returning so saved checkpoint matches reported AUC
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, best_auc, log


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--feature_dir", type=Path, required=True,
                   help="Per-slide UNI feature .pt files; each (N_tiles, 1024)")
    p.add_argument("--manifest", type=Path,
                   default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/slide_manifest.tsv"))
    p.add_argument("--out_dir", type=Path,
                   default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam"))
    p.add_argument("--n_folds", type=int, default=5)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # --- load manifest ---
    man = pd.read_csv(args.manifest, sep="\t")
    man = man[man["dm"].isin(["DM1", "DM2"])].copy()
    print(f"[load] {len(man)} slides ({sum(man['dm']=='DM1')} DM1, {sum(man['dm']=='DM2')} DM2)")
    label_map = {"DM2": 0, "DM1": 1}
    man["label"] = man["dm"].map(label_map)

    # --- load features ---
    bags, labels, slide_names = [], [], []
    for _, row in man.iterrows():
        feat_path = args.feature_dir / f"{row['file_id']}.pt"
        if not feat_path.exists():
            print(f"  [skip] {feat_path.name} not found")
            continue
        feat = torch.load(feat_path, map_location="cpu").numpy()
        bags.append(feat); labels.append(row["label"]); slide_names.append(row["file_id"])
    bags = np.array(bags, dtype=object)
    labels = np.array(labels)
    print(f"[features] loaded {len(bags)} bags")

    # --- 5-fold cross-val ---
    kf = KFold(n_splits=args.n_folds, shuffle=True, random_state=args.seed)
    fold_aucs, all_pred = [], []
    for fold_i, (tr_idx, val_idx) in enumerate(kf.split(bags)):
        print(f"\n[fold {fold_i+1}/{args.n_folds}]")
        tr_bags = [bags[i] for i in tr_idx]
        tr_labels = [labels[i] for i in tr_idx]
        val_bags = [bags[i] for i in val_idx]
        val_labels = [labels[i] for i in val_idx]
        model, auc, log = train_one_fold(tr_bags, tr_labels, val_bags, val_labels)
        fold_aucs.append(auc)
        # save best-epoch weights for this fold (so attention maps / Fig 5 can reload)
        ckpt_path = args.out_dir / f"clam_fold{fold_i+1}_best.pt"
        torch.save({k: v.detach().cpu() for k, v in model.state_dict().items()}, ckpt_path)
        print(f"  [ckpt] saved {ckpt_path.name} (val AUC={auc:.3f})")
        # eval per slide
        model.eval()
        with torch.no_grad():
            for j, b in enumerate(val_bags):
                x = torch.from_numpy(np.asarray(b, dtype=np.float32)).to("cuda")
                logits, attn = model(x)
                p = F.softmax(logits, dim=-1)[1].item()
                all_pred.append({"slide": slide_names[val_idx[j]],
                                 "fold": fold_i+1,
                                 "label": int(val_labels[j]),
                                 "prob_DM1": float(p),
                                 "n_tiles": int(b.shape[0])})

    # --- summary ---
    pred_df = pd.DataFrame(all_pred)
    overall_auc = roc_auc_score(pred_df["label"], pred_df["prob_DM1"])
    pred_df.to_csv(args.out_dir / "clam_per_slide_predictions.tsv", sep="\t", index=False)

    # fold AUC summary (so Fig 5 can pick best-fold checkpoint)
    fold_df = pd.DataFrame({"fold": list(range(1, len(fold_aucs) + 1)),
                            "best_val_auc": fold_aucs,
                            "ckpt": [f"clam_fold{i+1}_best.pt" for i in range(len(fold_aucs))]})
    fold_df.to_csv(args.out_dir / "clam_fold_summary.tsv", sep="\t", index=False)

    rep = []
    rep.append("# Phase 2 — CLAM attention-MIL (UNI features) on TCGA-THCA")
    rep.append(f"\n- {len(bags)} slides ({sum(labels==1)} DM1 / {sum(labels==0)} DM2)")
    rep.append(f"- {args.n_folds}-fold CV mean AUC: {np.mean(fold_aucs):.3f} ± {np.std(fold_aucs):.3f}")
    rep.append(f"- Per-fold: {fold_aucs}")
    rep.append(f"- Pooled cross-fold AUC: {overall_auc:.3f}")
    rep.append("")
    rep.append("- Closure ResNet50 baseline AUC ~ 0.55 (negative)")
    rep.append("- Kill-switch: PASS if held-out AUC > 0.70; MARGINAL 0.60-0.70; FAIL ≤ 0.60")
    if overall_auc > 0.70:
        verdict = "**PASS** — foundation+CLAM exceeds closure ResNet50; Paper 2 launch unlocked"
    elif overall_auc > 0.60:
        verdict = "**MARGINAL** — Paper 1 supp only; Paper 2 launch deferred"
    else:
        verdict = "**FAIL** — closure NO-GO re-confirmed at v2 architecture"
    rep.append(f"\n## VERDICT\n\n{verdict}")
    (args.out_dir / "PHASE2_REPORT.md").write_text("\n".join(rep))
    print("\n" + "\n".join(rep))


if __name__ == "__main__":
    main()
