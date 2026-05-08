"""
Supp #5 — Tile-count sensitivity.

For each DM1/DM2 slide, randomly subsample 50/100/150/200 tiles from the
existing 200-tile feature tensor and re-evaluate using the trained CLAM
fold 3 checkpoint (best held-out AUC=1.00). Pooled AUC at each tile count.

Saves:
  analysis_supp/tile_sensitivity.json
  analysis_supp/figures/figS4_tile_sensitivity.{png,pdf}
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
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

ROOT = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam")
FEATS = ROOT / "features"
MAN = ROOT / "slide_manifest.tsv"
CKPT = ROOT / "clam_fold3_best.pt"
OUT = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp")
FIGS = OUT / "figures"

TILE_COUNTS = [50, 100, 150, 200]
N_REPEAT = 20
SEED = 42
DEVICE = "cpu"


class GatedAttentionMIL(nn.Module):
    def __init__(self, in_dim=1024, hidden=512, n_classes=2, dropout=0.25):
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
        return self.classifier(bag), attn


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    torch.manual_seed(SEED)

    man = pd.read_csv(MAN, sep="\t")
    man = man[man["dm"].isin(["DM1", "DM2"])].copy()
    label_map = {"DM2": 0, "DM1": 1}
    man["label"] = man["dm"].map(label_map)

    # load model with fold 3 weights
    state = torch.load(CKPT, map_location=DEVICE)
    model = GatedAttentionMIL(n_classes=2).to(DEVICE)
    model.load_state_dict(state)
    model.eval()

    # gather features
    bags, labels, slides = [], [], []
    for _, row in man.iterrows():
        fp = FEATS / f"{row['file_id']}.pt"
        if not fp.exists():
            continue
        bags.append(torch.load(fp, map_location="cpu").numpy())
        labels.append(int(row["label"]))
        slides.append(row["file_id"])
    labels = np.array(labels, dtype=int)
    print(f"[tile_sens] {len(bags)} bags loaded; available tiles per slide: "
          f"{np.unique([b.shape[0] for b in bags])}")

    results = []
    for tc in TILE_COUNTS:
        per_repeat_auc = []
        for rep in range(N_REPEAT):
            probs = []
            for b in bags:
                n_avail = b.shape[0]
                k = min(tc, n_avail)
                if k == n_avail:
                    sub = b
                else:
                    idx = rng.choice(n_avail, size=k, replace=False)
                    sub = b[idx]
                with torch.no_grad():
                    x = torch.from_numpy(np.asarray(sub, dtype=np.float32)).to(DEVICE)
                    logits, _ = model(x)
                    p = float(F.softmax(logits, dim=-1)[1].item())
                probs.append(p)
            try:
                auc = float(roc_auc_score(labels, probs))
            except Exception:
                auc = float("nan")
            per_repeat_auc.append(auc)
            if tc == 200:
                # tile-count = 200 = full bag; deterministic so 1 rep enough
                break
        results.append({
            "tile_count": int(tc),
            "n_repeats": len(per_repeat_auc),
            "mean_auc": float(np.mean(per_repeat_auc)),
            "std_auc": float(np.std(per_repeat_auc, ddof=1)) if len(per_repeat_auc) > 1 else 0.0,
            "min_auc": float(np.min(per_repeat_auc)),
            "max_auc": float(np.max(per_repeat_auc)),
            "per_repeat_auc": per_repeat_auc,
        })
        print(f"  tc={tc}: mean_auc={np.mean(per_repeat_auc):.3f}  "
              f"std={np.std(per_repeat_auc, ddof=1) if len(per_repeat_auc)>1 else 0:.3f}")

    # plot
    xs = [r["tile_count"] for r in results]
    means = [r["mean_auc"] for r in results]
    stds = [r["std_auc"] for r in results]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.errorbar(xs, means, yerr=stds, fmt="o-", color="#1f77b4",
                ecolor="#666", capsize=4, markersize=8, linewidth=1.5)
    ax.set_xlabel("Tiles per slide (subsampled)")
    ax.set_ylabel("Pooled AUC (DM1 vs DM2, fold 3 ckpt)")
    ax.set_xticks(xs)
    ax.set_ylim(min(0.5, min(m - s for m, s in zip(means, stds)) - 0.05), 1.02)
    ax.set_title(f"Tile-count sensitivity — N_slides={len(bags)}, {N_REPEAT} reps per count")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "figS4_tile_sensitivity.png", dpi=300)
    fig.savefig(FIGS / "figS4_tile_sensitivity.pdf")
    plt.close(fig)

    out = {
        "n_slides": len(bags),
        "ckpt_used": str(CKPT),
        "tile_counts": TILE_COUNTS,
        "n_repeats_per_count": N_REPEAT,
        "results": results,
        "seed": SEED,
        "runtime_sec": float(time.time() - t0),
    }
    (OUT / "tile_sensitivity.json").write_text(json.dumps(out, indent=2))
    print(f"[tile_sens] runtime={out['runtime_sec']:.1f}s")


if __name__ == "__main__":
    main()
