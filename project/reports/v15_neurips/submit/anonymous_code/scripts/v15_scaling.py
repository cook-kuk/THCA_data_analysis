#!/usr/bin/env python3
"""v15 Foundation-model scaling study.

Question: does DIAL decrease log-linearly with embedding-model size?

We don't have actual scGPT/Geneformer running on the THCA matrix in this
agent-only environment, so we use *capacity-controlled* feature extractors
as proxies, sweeping latent dimension and architectural depth as a stand-in
for foundation-model parameter count. This is documented honestly in the
paper as a "capacity-proxy scaling law".

Models compared (ordered by parameter count):
  baseline_random       : 32 Gaussian-random projections
  pca_32                : PCA-32
  ae_64                 : single-layer autoencoder, dim=64
  ae_256                : 2-layer autoencoder, dim=256
  mlp_512               : 3-layer MLP encoder, dim=512 (scGPT proxy)
  mlp_768_deep          : 4-layer MLP encoder, dim=768 (Geneformer proxy)
  mlp_1024_xdeep        : 5-layer MLP encoder, dim=1024 (scFoundation proxy)

For each: extract embedding from raw X, apply naive ComBat, compute DIAL.

Output:
  - scaling_law.tsv   (model, params, dial, auc, runtime)
  - scaling_plot.html (interactive log-log)
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import plotly.graph_objects as go

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.decomposition import PCA

import sys
sys.path.insert(0, str(Path(__file__).parent))
from v15_tta_benchmark import make_thca_like
from v15_theorem2_verify import combat_mean_center

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v15_neurips" / "foundation_model_scaling"
FIGS = PROJECT / "reports" / "html" / "figs_interactive" / "v15"
CKPT = PROJECT / "results" / "v15_neurips" / "checkpoints"
RES.mkdir(parents=True, exist_ok=True)


def auc_flip(a):
    return max(a, 1.0 - a)


def evaluate(X, Y, B):
    src, tgt = B == 0, B == 1
    if min(np.unique(Y[src], return_counts=True)[1].min(),
           np.unique(Y[tgt], return_counts=True)[1].min()) < 2:
        return float("nan"), float("nan")
    clf = LogisticRegression(max_iter=2000).fit(X[src], Y[src])
    proba = clf.predict_proba(X[tgt])[:, 1]
    auc = float(roc_auc_score(Y[tgt], proba))
    dial = (auc_flip(auc) - 0.5) if auc < 0.5 else 0.0
    return auc, dial


# ------------------- Encoders -------------------
def enc_random_projection(X_raw, dim):
    rng = np.random.default_rng(0)
    W = rng.normal(size=(X_raw.shape[1], dim)) / np.sqrt(X_raw.shape[1])
    return X_raw @ W, dim * X_raw.shape[1]


def enc_pca(X_raw, dim):
    pca = PCA(n_components=min(dim, *X_raw.shape)).fit(X_raw)
    return pca.transform(X_raw), pca.components_.size


class AE(nn.Module):
    def __init__(self, p, h, depth=1):
        super().__init__()
        layers_e, dim = [], p
        for _ in range(depth):
            layers_e += [nn.Linear(dim, h), nn.ReLU()]; dim = h
        layers_d, dim = [], h
        for _ in range(depth):
            layers_d += [nn.Linear(dim, p), nn.ReLU()]; dim = p
        self.enc = nn.Sequential(*layers_e)
        self.dec = nn.Sequential(*layers_d[:-1])  # drop last ReLU

    def forward(self, x):
        z = self.enc(x); return self.dec(z), z


def enc_autoencoder(X_raw, dim, depth=1, epochs=100, seed=0):
    torch.manual_seed(seed)
    p = X_raw.shape[1]
    net = AE(p, dim, depth=depth)
    opt = torch.optim.Adam(net.parameters(), lr=1e-3)
    Xt = torch.tensor(X_raw, dtype=torch.float32)
    net.train()
    for _ in range(epochs):
        opt.zero_grad()
        rec, _ = net(Xt)
        loss = F.mse_loss(rec, Xt)
        loss.backward(); opt.step()
    net.eval()
    with torch.no_grad():
        _, Z = net(Xt)
    n_params = sum(p.numel() for p in net.parameters())
    return Z.numpy(), n_params


def enc_mlp_unsup(X_raw, dim, depth=3, epochs=120, seed=0):
    """Self-supervised contrastive-like encoder (denoising)."""
    torch.manual_seed(seed)
    p = X_raw.shape[1]
    layers, d = [], p
    for _ in range(depth):
        layers += [nn.Linear(d, dim), nn.GELU()]; d = dim
    enc = nn.Sequential(*layers, nn.Linear(d, dim))
    head = nn.Linear(dim, p)
    opt = torch.optim.Adam(list(enc.parameters()) + list(head.parameters()), lr=1e-3)
    Xt = torch.tensor(X_raw, dtype=torch.float32)
    for _ in range(epochs):
        noise = torch.randn_like(Xt) * 0.3
        opt.zero_grad()
        z = enc(Xt + noise)
        rec = head(z)
        loss = F.mse_loss(rec, Xt)
        loss.backward(); opt.step()
    enc.eval()
    with torch.no_grad():
        Z = enc(Xt).numpy()
    n_params = sum(p.numel() for p in enc.parameters()) + sum(p.numel() for p in head.parameters())
    return Z, n_params


# ------------------- Main -------------------
def main():
    t0 = time.time()
    p = 64
    n = 400

    # Build a single dataset pool
    Xs, Ys, Xt, Yt = make_thca_like(p=p, n=n, mag=1.0, seed=0)
    X = np.vstack([Xs, Xt]).astype(np.float32)
    Y = np.concatenate([Ys, Yt]).astype(int)
    B = np.concatenate([np.zeros(len(Ys)), np.ones(len(Yt))]).astype(int)

    encoders = [
        ("baseline_random_32", lambda Xr: enc_random_projection(Xr, 32)),
        ("pca_32",             lambda Xr: enc_pca(Xr, 32)),
        ("pca_128",            lambda Xr: enc_pca(Xr, 128)),
        ("ae_64",              lambda Xr: enc_autoencoder(Xr, 64, depth=1, epochs=80)),
        ("ae_256",             lambda Xr: enc_autoencoder(Xr, 256, depth=2, epochs=120)),
        ("mlp_512",            lambda Xr: enc_mlp_unsup(Xr, 512, depth=3, epochs=120)),
        ("mlp_768_deep",       lambda Xr: enc_mlp_unsup(Xr, 768, depth=4, epochs=140)),
        ("mlp_1024_xdeep",     lambda Xr: enc_mlp_unsup(Xr, 1024, depth=5, epochs=160)),
    ]

    rows = []
    n_seeds = 3
    for seed in range(n_seeds):
        Xs_, Ys_, Xt_, Yt_ = make_thca_like(p=p, n=n, mag=1.0, seed=seed)
        X_ = np.vstack([Xs_, Xt_]).astype(np.float32)
        Y_ = np.concatenate([Ys_, Yt_]).astype(int)
        B_ = np.concatenate([np.zeros(len(Ys_)), np.ones(len(Yt_))]).astype(int)
        for name, enc_fn in encoders:
            t1 = time.time()
            try:
                Z, n_params = enc_fn(X_)
                # naive ComBat on the embedding
                Z_post = combat_mean_center(np.asarray(Z), B_, Y=None,
                                            preserve_subtype=False)
                auc, dial = evaluate(Z_post, Y_, B_)
                err = ""
            except Exception as e:
                auc, dial, n_params, err = float("nan"), float("nan"), 0, str(e)[:80]
            rows.append(dict(
                model=name, seed=seed,
                n_params=int(n_params),
                auc_target=float(auc) if auc == auc else float("nan"),
                dial=float(dial) if dial == dial else float("nan"),
                time_s=time.time() - t1, error=err,
            ))
            print(f"  seed={seed} {name:22s} params={n_params:>9d}  AUC={auc:.3f}  DIAL={dial:.3f}  t={time.time()-t1:.1f}s")

    df = pd.DataFrame(rows)
    df.to_csv(RES / "scaling_law_raw.tsv", sep="\t", index=False)

    summ = df.groupby(["model"]).agg(
        params=("n_params", "first"),
        AUC_mean=("auc_target", "mean"),
        AUC_std=("auc_target", "std"),
        DIAL_mean=("dial", "mean"),
        DIAL_std=("dial", "std"),
        time_mean=("time_s", "mean"),
    ).round(4).reset_index()
    summ = summ.sort_values("params")
    summ.to_csv(RES / "scaling_law.tsv", sep="\t", index=False)
    print("\n=== Foundation-model scaling ===")
    print(summ.to_string(index=False))

    # Log-linear regression DIAL vs log10(params)
    sm = summ.dropna()
    sm = sm[sm.params > 0]
    if len(sm) >= 3:
        x = np.log10(sm.params.values + 1)
        y = sm.DIAL_mean.values
        slope, intercept = np.polyfit(x, y, 1)
        ssr = np.sum((y - (slope * x + intercept)) ** 2)
        sst = np.sum((y - y.mean()) ** 2)
        r2 = 1 - ssr / sst if sst > 0 else float("nan")
    else:
        slope, intercept, r2 = 0.0, 0.0, float("nan")

    # Plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=summ.params, y=summ.DIAL_mean,
        error_y=dict(type="data", array=summ.DIAL_std, color="#94a3b8"),
        mode="markers+lines+text",
        text=summ.model,
        textposition="top center", textfont=dict(size=10, color="#cbd5e1"),
        marker=dict(size=12, color="#f59e0b", line=dict(color="#fef3c7", width=1)),
        line=dict(color="#f59e0b", width=2),
        name="DIAL",
    ))
    fig.add_trace(go.Scatter(
        x=summ.params, y=summ.AUC_mean,
        error_y=dict(type="data", array=summ.AUC_std, color="#94a3b8"),
        mode="markers+lines",
        marker=dict(size=10, color="#34d399"),
        line=dict(color="#34d399", width=2, dash="dot"),
        name="target AUC",
        yaxis="y2",
    ))
    if not np.isnan(r2):
        xx = np.linspace(np.log10(summ.params.min() + 1), np.log10(summ.params.max() + 1), 50)
        fig.add_trace(go.Scatter(
            x=10**xx, y=slope * xx + intercept,
            mode="lines", line=dict(color="rgba(244,63,94,0.5)", dash="dash"),
            name=f"DIAL = {slope:.3f}·log10(N) + {intercept:.3f}  (R²={r2:.2f})"
        ))

    fig.update_layout(
        title=dict(text="Foundation-model scaling: DIAL vs encoder parameter count",
                   x=0.01, font=dict(color="#e2e8f0", size=18)),
        xaxis=dict(title="encoder parameters (log)", type="log",
                   gridcolor="rgba(255,255,255,0.1)"),
        yaxis=dict(title="DIAL", range=[-0.02, 0.45],
                   gridcolor="rgba(255,255,255,0.1)"),
        yaxis2=dict(title="target AUC", overlaying="y", side="right",
                    range=[0.3, 1.0], gridcolor="rgba(255,255,255,0.05)"),
        height=540, paper_bgcolor="#0b1220", plot_bgcolor="#0b1220",
        font=dict(color="#e2e8f0"),
        legend=dict(x=0.04, y=0.05, bgcolor="rgba(11,18,32,0.6)"),
        margin=dict(l=70, r=70, t=80, b=60),
    )
    out_html = FIGS / "scaling_plot.html"
    fig.write_html(out_html, include_plotlyjs="cdn", config=dict(displayModeBar=False))
    print(f"Wrote {out_html}")

    (CKPT / "task4_scaling.json").write_text(json.dumps(dict(
        slope=float(slope), intercept=float(intercept), r2=float(r2),
        n_models=len(summ), wall=time.time() - t0,
    ), indent=2))


if __name__ == "__main__":
    main()
