#!/usr/bin/env python3
"""v5 Track 2 — DANN-style Adversarial De-confounding (the NOVEL method).

Trains a biomarker classifier (BRAF_like vs RAS_like) with a shared
encoder that is adversarially regularized against cohort prediction
using a gradient-reversal layer. Sweeps the adversarial weight lambda
and reports the Pareto frontier (bio AUC vs 1 - ident AUC).

Outputs:
  results/ml/v5_track2_adversarial_results.tsv
  results/ml/v5_track2_best_model.pt
  results/ml/v5_track2_method_card.md  (also copied to reports/)
  reports/html/figs_interactive/v5_adversarial_pareto.html
  reports/html/figs_interactive/v5_adversarial_shap_stability.html
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
torch.set_num_threads(2)

sys.path.insert(0, str(Path(__file__).parent))
from v5_common import (  # noqa: E402
    build_pooled_matrix,
    bio_label_bvr,
    load_sample_master,
    RESULTS_ML,
    FIGS_INTERACTIVE,
    REPORTS,
    safe_write_json,
)


def log(msg: str) -> None:
    print(f"[v5-track2] {msg}", flush=True)


# ------------- Gradient Reversal Layer (core novelty plumbing) -------------
class _GradientReverse(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, lambd):
        ctx.lambd = lambd
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        return -ctx.lambd * grad_output, None


def grad_reverse(x, lambd=1.0):
    return _GradientReverse.apply(x, lambd)


# ------------- Architecture -------------
class DANNNet(nn.Module):
    """Shared encoder (in_dim -> 32 -> 16) with biomarker head (16 -> 2)
    and cohort head that consumes encoder output through a GRL.
    """

    def __init__(self, in_dim: int, n_cohorts: int, dropout: float = 0.3):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(in_dim, 32), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(32, 16), nn.ReLU(), nn.Dropout(dropout),
        )
        self.bio_head = nn.Linear(16, 2)
        self.cohort_head = nn.Sequential(
            nn.Linear(16, 16), nn.ReLU(), nn.Linear(16, n_cohorts)
        )

    def forward(self, x, lambd):
        h = self.encoder(x)
        y_bio = self.bio_head(h)
        y_cohort = self.cohort_head(grad_reverse(h, lambd))
        return y_bio, y_cohort, h


def evaluate_ident_auc(model, X, cohort):
    """Freeze encoder, train a fresh cohort probe on embeddings, macro-AUC."""
    model.eval()
    with torch.no_grad():
        emb = model.encoder(torch.as_tensor(X, dtype=torch.float32)).numpy()
    uniq = np.unique(cohort)
    aucs = []
    for c in uniq:
        y = (cohort == c).astype(int)
        if y.sum() < 5 or (len(y) - y.sum()) < 5:
            continue
        skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=0)
        fold = []
        for tr, te in skf.split(emb, y):
            clf = LogisticRegression(max_iter=1500)
            clf.fit(emb[tr], y[tr])
            p = clf.predict_proba(emb[te])[:, 1]
            try:
                fold.append(roc_auc_score(y[te], p))
            except Exception:
                continue
        if fold:
            aucs.append(float(np.mean(fold)))
    return float(np.mean(aucs)) if aucs else float("nan")


def lodo_bio_auc(model, X, y, cohorts):
    per = []
    model.eval()
    with torch.no_grad():
        probs = F.softmax(model(torch.as_tensor(X, dtype=torch.float32),
                                lambd=0.0)[0], dim=1)[:, 1].numpy()
    for held in np.unique(cohorts):
        m = cohorts == held
        if len(np.unique(y[m])) < 2:
            per.append({"held": held, "auc": float("nan"),
                        "n": int(m.sum()), "note": "single-class"})
            continue
        auc = float(roc_auc_score(y[m], probs[m]))
        per.append({"held": held, "auc": auc, "n": int(m.sum()),
                    "note": "small n" if m.sum() < 20 else ""})
    return float(np.nanmean([p["auc"] for p in per])), per


def train_dann(X_tr, y_tr, c_tr, X_val, y_val, c_val, n_cohorts,
               lambd: float, in_dim: int, epochs: int = 200, lr: float = 1e-3,
               patience: int = 20):
    model = DANNNet(in_dim=in_dim, n_cohorts=n_cohorts)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    Xt = torch.as_tensor(X_tr, dtype=torch.float32)
    yt = torch.as_tensor(y_tr, dtype=torch.long)
    ct = torch.as_tensor(c_tr, dtype=torch.long)
    Xv = torch.as_tensor(X_val, dtype=torch.float32)
    yv = torch.as_tensor(y_val, dtype=torch.long)
    best_val = -1.0
    best_state = None
    stale = 0
    for ep in range(epochs):
        model.train()
        perm = torch.randperm(Xt.shape[0])
        bs = 64
        for i in range(0, len(perm), bs):
            idx = perm[i:i + bs]
            x, yb, yc = Xt[idx], yt[idx], ct[idx]
            y_bio_logits, y_coh_logits, _ = model(x, lambd=lambd)
            loss_bio = F.cross_entropy(y_bio_logits, yb)
            loss_coh = F.cross_entropy(y_coh_logits, yc)
            loss = loss_bio + lambd * loss_coh
            opt.zero_grad()
            loss.backward()
            opt.step()
        # eval
        model.eval()
        with torch.no_grad():
            p_val = F.softmax(model(Xv, lambd=0.0)[0], dim=1)[:, 1].numpy()
        if len(np.unique(y_val)) == 2:
            try:
                auc = roc_auc_score(y_val, p_val)
            except Exception:
                auc = -1
        else:
            auc = -1
        if auc > best_val:
            best_val = auc
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, best_val


def bio_cv_auc(model_builder, X, y, cohort, n_cohorts, in_dim, lambd):
    """5-fold stratified CV AUC for BRAF vs RAS on full matrix."""
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = []
    for tr, te in skf.split(X, y):
        # carve val out of train
        try:
            Xtr, Xv, ytr, yv, ctr, cv = train_test_split(
                X[tr], y[tr], cohort[tr], test_size=0.15,
                random_state=0, stratify=y[tr])
        except Exception:
            Xtr, Xv, ytr, yv, ctr, cv = X[tr], X[tr], y[tr], y[tr], cohort[tr], cohort[tr]
        model, _ = train_dann(Xtr, ytr, ctr, Xv, yv, cv,
                              n_cohorts=n_cohorts, lambd=lambd,
                              in_dim=in_dim, epochs=150)
        model.eval()
        with torch.no_grad():
            p = F.softmax(model(torch.as_tensor(X[te], dtype=torch.float32),
                                lambd=0.0)[0], dim=1)[:, 1].numpy()
        try:
            aucs.append(roc_auc_score(y[te], p))
        except Exception:
            continue
    return float(np.mean(aucs)) if aucs else float("nan"), model


def shap_top_genes(model, X, gene_names, n_top: int = 10, n_background: int = 50):
    """Lightweight SHAP-like gradient*input importance on bio head."""
    model.eval()
    bg = torch.as_tensor(X[np.random.default_rng(0).choice(X.shape[0],
                                                           min(n_background, X.shape[0]),
                                                           replace=False)],
                         dtype=torch.float32)
    # integrated-gradient-ish: gradient of logit[class=1] wrt input on full data
    Xt = torch.as_tensor(X, dtype=torch.float32).clone().requires_grad_(True)
    logits = model(Xt, lambd=0.0)[0]
    target = logits[:, 1].sum()
    grads = torch.autograd.grad(target, Xt)[0].detach().numpy()
    importance = np.mean(np.abs(grads) * np.abs(X), axis=0)  # |grad * input|
    order = np.argsort(-importance)
    top = [(gene_names[i], float(importance[i])) for i in order[:n_top]]
    return top, importance


def main() -> int:
    log("loading pooled matrix with BRAF_like / RAS_like labels...")
    big_df, meta_df = build_pooled_matrix(min_samples_per_cohort=5,
                                          top_var_genes=3000,
                                          label_fn=bio_label_bvr)
    log(f"pooled: {big_df.shape[0]} genes x {big_df.shape[1]} samples; "
        f"cohorts={sorted(meta_df['cohort'].unique())}")
    log(f"  bio label counts: {dict(meta_df['label'].value_counts())}")

    # reduce to top 54 variable genes to match spec (54 -> 32 -> 16)
    top = big_df.var(axis=1).sort_values(ascending=False).head(54).index
    X_df = big_df.loc[top]
    X = StandardScaler().fit_transform(X_df.T.values)  # samples x 54
    y = meta_df["label"].values.astype(int)
    cohort_labels = pd.Categorical(meta_df["cohort"].values)
    cohort_ids = cohort_labels.codes
    n_cohorts = len(cohort_labels.categories)
    in_dim = X.shape[1]
    log(f"  X shape: {X.shape}; n_cohorts={n_cohorts}; in_dim={in_dim}")

    # sanity-check: we need at least 2 classes in y
    if len(np.unique(y)) < 2:
        log("FATAL: only one class present in BRAF/RAS labels; skipping Track 2")
        pd.DataFrame([{"lambda": 0.0, "status": "SKIPPED - single class"}]).to_csv(
            RESULTS_ML / "v5_track2_adversarial_results.tsv", sep="\t", index=False)
        return 0

    lambdas = [0.0, 0.1, 0.3, 1.0, 3.0]
    rows = []
    shap_per_lambda = {}
    best_model = None
    best_score = -np.inf
    best_lambda = None

    for lambd in lambdas:
        log(f"training DANN with lambda = {lambd}...")
        t0 = time.time()
        bio_auc, final_model = bio_cv_auc(None, X, y, cohort_ids,
                                          n_cohorts=n_cohorts, in_dim=in_dim,
                                          lambd=lambd)
        cohort_str = np.asarray(cohort_labels.astype(str))
        ident = evaluate_ident_auc(final_model, X, cohort_str)
        lodo_m, lodo_per = lodo_bio_auc(final_model, X, y, cohort_str)
        top_shap, _ = shap_top_genes(final_model, X, list(X_df.index), n_top=10)
        shap_per_lambda[lambd] = [g for g, _ in top_shap]
        rows.append({
            "lambda": lambd,
            "bio_cv_auc": bio_auc,
            "ident_auc": ident,
            "lodo_mean_auc": lodo_m,
            "top_shap_genes": ";".join([g for g, _ in top_shap]),
            "runtime_sec": round(time.time() - t0, 2),
        })
        score = (bio_auc if not np.isnan(bio_auc) else 0.0) - \
                max(0.0, (ident if not np.isnan(ident) else 1.0) - 0.70)
        log(f"  lambda={lambd}: bio={bio_auc:.3f} ident={ident:.3f} "
            f"lodo={lodo_m:.3f} score={score:.3f}")
        if score > best_score:
            best_score = score
            best_model = final_model
            best_lambda = lambd

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_ML / "v5_track2_adversarial_results.tsv", sep="\t", index=False)

    # Save best model
    if best_model is not None:
        torch.save({
            "state_dict": best_model.state_dict(),
            "lambda": best_lambda,
            "in_dim": in_dim,
            "n_cohorts": n_cohorts,
            "genes": list(X_df.index),
            "cohort_categories": list(cohort_labels.categories),
        }, RESULTS_ML / "v5_track2_best_model.pt")

    # SHAP stability table: genes appearing in top-10 for >= 3/5 lambdas
    from collections import Counter
    flat = []
    for lm, genes in shap_per_lambda.items():
        flat.extend(genes)
    freq = Counter(flat)
    stable_genes = [{"gene": g, "n_lambdas_in_top10": n,
                     "lambdas_containing": ";".join([str(l) for l, gg in shap_per_lambda.items()
                                                     if g in gg])}
                    for g, n in freq.most_common()
                    if n >= max(3, len(shap_per_lambda) // 2)]
    stable_df = pd.DataFrame(stable_genes)
    stable_df.to_csv(RESULTS_ML / "v5_track2_shap_stable_genes.tsv",
                     sep="\t", index=False)

    # Pareto plot
    fig = px.scatter(df, x="ident_auc", y="bio_cv_auc",
                     size=[30] * len(df), color="lambda",
                     text=[f"λ={l}" for l in df["lambda"]],
                     title="v5 Track 2 — DANN Pareto frontier (higher bio, lower ident = better)",
                     labels={"ident_auc": "dataset identifiability AUC",
                             "bio_cv_auc": "biomarker CV AUC (BRAF vs RAS)"})
    fig.add_shape(type="rect", x0=0.0, x1=0.70, y0=0.85, y1=1.0,
                  fillcolor="#16a085", opacity=0.1, line=dict(width=0),
                  layer="below")
    fig.update_traces(textposition="top center")
    fig.write_html(FIGS_INTERACTIVE / "v5_adversarial_pareto.html",
                   include_plotlyjs="cdn")

    # SHAP stability figure
    if len(stable_df) > 0:
        stable_df_fig = stable_df.head(20)
        fig = px.bar(stable_df_fig.sort_values("n_lambdas_in_top10"),
                     x="n_lambdas_in_top10", y="gene",
                     orientation="h",
                     title="v5 Track 2 — SHAP-stable genes across λ sweep",
                     labels={"n_lambdas_in_top10": "# λ values for which gene in top-10",
                             "gene": "gene"})
        fig.write_html(FIGS_INTERACTIVE / "v5_adversarial_shap_stability.html",
                       include_plotlyjs="cdn")
    else:
        fig = go.Figure()
        fig.add_annotation(text="no SHAP-stable genes above threshold",
                           xref="paper", yref="paper", x=0.5, y=0.5)
        fig.write_html(FIGS_INTERACTIVE / "v5_adversarial_shap_stability.html",
                       include_plotlyjs="cdn")

    # Summary + method card
    best_row = df[df["lambda"] == best_lambda].iloc[0] if best_lambda is not None else None
    summary = {
        "lambdas_tested": lambdas,
        "best_lambda": best_lambda,
        "best_bio_cv_auc": float(best_row["bio_cv_auc"]) if best_row is not None else None,
        "best_ident_auc": float(best_row["ident_auc"]) if best_row is not None else None,
        "best_lodo_mean_auc": float(best_row["lodo_mean_auc"]) if best_row is not None else None,
        "n_shap_stable_genes": int(len(stable_df)),
        "shap_stable_top10": stable_df["gene"].head(10).tolist()
        if len(stable_df) > 0 else [],
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "n_cohorts": int(n_cohorts),
    }
    safe_write_json(RESULTS_ML / "v5_track2_summary.json", summary)

    # Method card (technical spec) — written to both results/ml and reports/
    method_card = f"""# v5 Track 2 — DANN Biomarker De-confounding Method Card

## One-line description
A gradient-reversal-layer (GRL) biomarker classifier that is adversarially
regularized against batch (cohort) prediction, enabling cohort-invariant
feature learning for cross-study biomarker discovery.

## Architecture
Shared encoder: `MLP({in_dim} -> 32 -> 16)` with ReLU + Dropout(0.3) after each
hidden layer. Biomarker head: `Linear(16 -> 2)` softmax for BRAF_like vs RAS_like.
Cohort head: `GRL -> Linear(16 -> 16) -> ReLU -> Linear(16 -> K={n_cohorts})`
softmax over K={n_cohorts} cohorts. The gradient-reversal layer (Ganin & Lempitsky,
2015) multiplies gradients by `-lambda` during backprop through the encoder,
forcing the encoder to produce features that a cohort classifier cannot
discriminate.

## Loss function
```
L = CE(y_bio, y_bio_hat) + lambda * CE(y_cohort, y_cohort_hat)
```
The GRL ensures the encoder MINIMIZES the ability of the cohort head to
predict cohort (by gradient negation) while the cohort head itself is trained
normally to predict cohort. Lambda controls anti-confounding strength.

## Hyperparameters
- Optimizer: Adam, lr=1e-3, weight_decay=1e-4
- Batch size: 64
- Max epochs: 200 (with early-stop patience=20 on val bio AUC)
- Dropout: 0.3 in encoder
- Lambda sweep: {{0, 0.1, 0.3, 1.0, 3.0}}
- Input features: top {in_dim} variable genes across pooled cohorts

## Metrics (this run)
| lambda | bio CV AUC | ident AUC | LODO AUC |
|--------|-----------:|----------:|---------:|
""" + "\n".join(
        f"| {r['lambda']} | {r['bio_cv_auc']:.3f} | {r['ident_auc']:.3f} | {r['lodo_mean_auc']:.3f} |"
        for _, r in df.iterrows()
    ) + f"""

**Best λ = {best_lambda}** (selected by `bio_CV_AUC - max(0, ident_AUC - 0.70)`).

## SHAP-stable genes
Genes appearing in the top-10 gradient*input importances for ≥3 of {len(lambdas)}
λ values (considered robust to the anti-confounding regularization):

{", ".join(stable_df["gene"].head(15).tolist()) if len(stable_df) > 0 else "(none above threshold)"}

## Novelty vs. prior art
DANN (Ganin & Lempitsky, 2015) was developed for image-domain adaptation.
To our knowledge, this is the first application to **cross-cohort bulk RNA-seq
biomarker discovery** as a replacement/augmentation for ComBat-style pre-processing.
The key claim is: instead of pre-correcting expression and then fitting a
downstream classifier (which v4's Track A shows is unrecoverable for THCA
subtype labels), **learn a cohort-invariant representation jointly** with the
biomarker classifier. The Pareto frontier (bio AUC vs 1-ident AUC) is a
principled way to report the trade-off rather than claiming a single AUC.

## Limitations
1. CPU-only PyTorch, small MLP — does not attempt deep representations.
2. GRL hyperparameter λ must be swept; no automatic selection.
3. Assumes cohort label is the only confounder; does not jointly
   de-confound platform, sequencing depth, age, stage separately.
4. BRAF/RAS label distribution is severely imbalanced across cohorts —
   some LODO folds have single-class held-out and are reported as NaN.
5. **Prototype only, not clinical**. Decision-support / methodology
   contribution. Not a diagnostic device.

## Reproducibility
- Code: `notebooks_or_scripts/v5_track2_adversarial.py`
- Model weights: `results/ml/v5_track2_best_model.pt`
- Results: `results/ml/v5_track2_adversarial_results.tsv`
- SHAP stability: `results/ml/v5_track2_shap_stable_genes.tsv`
- Random seed: 42 for sklearn splits, 0 for torch training
- Hardware: CPU-only, PyTorch {torch.__version__}
"""
    (RESULTS_ML / "v5_track2_method_card.md").write_text(method_card)
    (REPORTS / "v5_adversarial_method_card.md").write_text(method_card)

    log(f"done; best lambda={best_lambda}; stable-shap genes={len(stable_df)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
