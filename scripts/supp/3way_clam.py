"""
Supp #3 — 3-way CLAM (DM1 vs DM2 vs not_DM), 5-fold stratified CV.

Reads: features/, slide_manifest.tsv
Saves:
  analysis_supp/3way_classification.json
  analysis_supp/figures/figS2_3way_confusion.{png,pdf}
"""
from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import (
    roc_auc_score,
    confusion_matrix,
    accuracy_score,
)
from sklearn.model_selection import StratifiedKFold
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam")
FEATS = ROOT / "features"
MAN = ROOT / "slide_manifest.tsv"
OUT = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp")
FIGS = OUT / "figures"

N_FOLDS = 5
EPOCHS = 30
LR = 2e-4
SEED = 42
DEVICE = "cpu"

LABEL_MAP = {"DM2": 0, "DM1": 1, "not_DM": 2}
INV_LABEL = {v: k for k, v in LABEL_MAP.items()}


class GatedAttentionMIL(nn.Module):
    def __init__(self, in_dim=1024, hidden=512, n_classes=3, dropout=0.25):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.attn_a = nn.Sequential(nn.Linear(hidden, hidden), nn.Tanh())
        self.attn_b = nn.Sequential(nn.Linear(hidden, hidden), nn.Sigmoid())
        self.attn_w = nn.Linear(hidden, 1)
        self.classifier = nn.Linear(hidden, n_classes)

    def forward(self, x):
        h = self.fc(x)
        a = self.attn_w(self.attn_a(h) * self.attn_b(h))
        attn = F.softmax(a.squeeze(-1), dim=0)
        bag = (attn.unsqueeze(-1) * h).sum(0)
        logits = self.classifier(bag)
        return logits, attn


def train_fold(tr_bags, tr_y, val_bags, val_y, n_classes, device=DEVICE,
               epochs=EPOCHS, lr=LR):
    model = GatedAttentionMIL(in_dim=tr_bags[0].shape[1], n_classes=n_classes).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    crit = nn.CrossEntropyLoss()
    best_score = -1.0
    best_state = None
    best_probs = None
    for epoch in range(epochs):
        model.train()
        for bag, label in zip(tr_bags, tr_y):
            opt.zero_grad()
            x = torch.from_numpy(np.asarray(bag, dtype=np.float32)).to(device)
            y_t = torch.tensor([label]).to(device)
            logits, _ = model(x)
            loss = crit(logits.unsqueeze(0), y_t)
            loss.backward()
            opt.step()
        model.eval()
        probs = []
        with torch.no_grad():
            for bag in val_bags:
                x = torch.from_numpy(np.asarray(bag, dtype=np.float32)).to(device)
                logits, _ = model(x)
                p = F.softmax(logits, dim=-1).cpu().numpy()
                probs.append(p)
        probs = np.array(probs)
        try:
            score = roc_auc_score(val_y, probs, multi_class="ovr", average="macro",
                                  labels=list(range(n_classes)))
        except Exception:
            score = -1.0
        if score > best_score:
            best_score = score
            best_probs = probs
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
    return best_state, best_score, best_probs


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    man = pd.read_csv(MAN, sep="\t")
    man = man[man["dm"].isin(LABEL_MAP)].copy()
    man["label"] = man["dm"].map(LABEL_MAP)

    bags, labels, slides = [], [], []
    for _, row in man.iterrows():
        fp = FEATS / f"{row['file_id']}.pt"
        if not fp.exists():
            continue
        feat = torch.load(fp, map_location="cpu").numpy()
        bags.append(feat)
        labels.append(int(row["label"]))
        slides.append(row["file_id"])
    bags = np.array(bags, dtype=object)
    labels = np.array(labels, dtype=int)
    print(f"[3way] N={len(bags)}  per class: {dict(zip(*np.unique(labels, return_counts=True)))}")

    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    n_classes = 3
    fold_scores = []
    all_pred = np.zeros((len(bags), n_classes), dtype=float)
    all_y = labels.copy()

    for fold_i, (tr_idx, val_idx) in enumerate(skf.split(np.arange(len(bags)), labels)):
        tr_bags = [bags[i] for i in tr_idx]
        tr_y = [int(labels[i]) for i in tr_idx]
        val_bags = [bags[i] for i in val_idx]
        val_y = [int(labels[i]) for i in val_idx]
        state, score, probs = train_fold(tr_bags, tr_y, val_bags, val_y, n_classes)
        fold_scores.append(float(score))
        for j, vi in enumerate(val_idx):
            all_pred[vi] = probs[j]
        torch.save(state, OUT / f"3way_clam_fold{fold_i+1}_best.pt")
        print(f"  fold {fold_i+1}: macro-AUC={score:.3f}")

    # Per-class one-vs-rest
    per_class_auc = {}
    for cls_idx in range(n_classes):
        try:
            auc_c = roc_auc_score((all_y == cls_idx).astype(int), all_pred[:, cls_idx])
        except Exception:
            auc_c = float("nan")
        per_class_auc[INV_LABEL[cls_idx]] = float(auc_c)
    macro_auc = float(np.nanmean(list(per_class_auc.values())))
    pred_class = all_pred.argmax(axis=1)
    acc = float(accuracy_score(all_y, pred_class))
    cm = confusion_matrix(all_y, pred_class, labels=list(range(n_classes)))

    # plot confusion matrix
    fig, ax = plt.subplots(figsize=(5.5, 4.8))
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks(range(n_classes))
    ax.set_yticks(range(n_classes))
    ax.set_xticklabels([INV_LABEL[i] for i in range(n_classes)])
    ax.set_yticklabels([INV_LABEL[i] for i in range(n_classes)])
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title(
        f"3-way CLAM confusion (N={len(all_y)})\n"
        f"Acc={acc:.3f}  Macro-AUC={macro_auc:.3f}"
    )
    for i in range(n_classes):
        for j in range(n_classes):
            ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                    color="white" if cm[i, j] > cm.max() / 2 else "black", fontsize=12)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(FIGS / "figS2_3way_confusion.png", dpi=300)
    fig.savefig(FIGS / "figS2_3way_confusion.pdf")
    plt.close(fig)

    res = {
        "n_slides": int(len(all_y)),
        "class_counts": {INV_LABEL[k]: int(v) for k, v in zip(*np.unique(all_y, return_counts=True))},
        "n_folds": N_FOLDS,
        "fold_macro_aucs": fold_scores,
        "per_class_ovr_auc": per_class_auc,
        "macro_auc": macro_auc,
        "accuracy": acc,
        "confusion_matrix": cm.tolist(),
        "label_order": [INV_LABEL[i] for i in range(n_classes)],
        "epochs_per_fold": EPOCHS,
        "lr": LR,
        "seed": SEED,
        "runtime_sec": float(time.time() - t0),
    }
    (OUT / "3way_classification.json").write_text(json.dumps(res, indent=2))
    print(f"[3way] macro-AUC={macro_auc:.3f}  acc={acc:.3f}  "
          f"per-class={per_class_auc}  runtime={res['runtime_sec']:.1f}s")


if __name__ == "__main__":
    main()
