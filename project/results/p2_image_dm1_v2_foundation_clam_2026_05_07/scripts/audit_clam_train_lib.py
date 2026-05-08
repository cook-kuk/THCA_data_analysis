"""
Shared CLAM training library — extracted from phase2_clam_train_eval.py
so audit scripts can call train_one_fold() directly with custom splits.

Identical architecture & hyperparameters to the published model:
  GatedAttentionMIL(in_dim=1024, hidden=512, n_classes=2, dropout=0.25)
  AdamW(lr=2e-4, weight_decay=1e-4), CrossEntropyLoss, 30 epochs
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score


class GatedAttentionMIL(nn.Module):
    def __init__(self, in_dim=1024, hidden=512, n_classes=2, dropout=0.25):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU(),
                                nn.Dropout(dropout))
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


def train_one_fold(train_bags, train_labels, val_bags, val_labels,
                   epochs=30, lr=2e-4, device="cpu", seed=42):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = GatedAttentionMIL(in_dim=train_bags[0].shape[1],
                              n_classes=2).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    crit = nn.CrossEntropyLoss()
    best_auc = 0.0
    best_state = None
    best_probs = None
    for epoch in range(epochs):
        model.train()
        for bag, label in zip(train_bags, train_labels):
            opt.zero_grad()
            x = torch.from_numpy(np.asarray(bag, dtype=np.float32)).to(device)
            y = torch.tensor([label]).to(device)
            logits, _ = model(x)
            loss = crit(logits.unsqueeze(0), y)
            loss.backward(); opt.step()
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
            auc = roc_auc_score(val_labels, probs) if \
                len(set(val_labels)) > 1 else float("nan")
        except Exception:
            auc = float("nan")
        if not np.isnan(auc) and auc > best_auc:
            best_auc = auc
            best_probs = list(probs)
        # always keep latest probs as fallback
        if best_probs is None:
            best_probs = list(probs)
    return best_auc, best_probs


def load_features(feat_dir: Path, file_ids: list[str]):
    bags = []
    for fid in file_ids:
        feat = torch.load(feat_dir / f"{fid}.pt",
                          map_location="cpu", weights_only=False)
        bags.append(feat.numpy() if hasattr(feat, "numpy") else np.asarray(feat))
    return bags
