#!/usr/bin/env python3
"""v15 TTA benchmark.

Compare classical (no-adapt, ComBat, CORAL), adversarial DA (DANN-lite),
test-time adaptation (TENT, BN-stats, SHOT-lite), and a foundation-model
TTA proxy on a synthetic source→target benchmark designed to mimic the
THCA TCGA→GSE flip.

Output:
  - methods_comparison.tsv (per-method DIAL, AUC, runtime)
  - tta_benchmark.html (bar chart, dark theme)
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

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v15_neurips" / "tta_benchmark"
FIGS = PROJECT / "reports" / "html" / "figs_interactive" / "v15"
CKPT = PROJECT / "results" / "v15_neurips" / "checkpoints"
RES.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------
# Data — THCA-style synthetic with subspace-aligned conditional shift
# --------------------------------------------------------------------
def make_thca_like(p=64, n=300, mag=0.8, align=0.7, seed=0):
    rng = np.random.default_rng(seed)
    wY = np.zeros(p); wY[0] = 1.0
    wB = align * wY + np.sqrt(1 - align**2) * np.eye(p)[1]

    # Class means with clear separation so source CV-AUC is high
    sep = 1.2
    Xs0 = rng.normal(0, 1, size=(n // 2, p)) - sep * wY
    Xs1 = rng.normal(0, 1, size=(n // 2, p)) + sep * wY

    # Target: skewed + class-conditional residual along wY (THCA flip).
    skew = 0.2
    nt0 = int(n // 2 * (1 - skew)); nt1 = int(n // 2 * (1 + skew))
    # Residual shift is 1.6x the source separation — induces full label flip
    Xt0 = rng.normal(0, 1, size=(nt0, p)) - sep * wY + 0.6 * mag * wB + 2.0 * mag * wY
    Xt1 = rng.normal(0, 1, size=(nt1, p)) + sep * wY + 0.6 * mag * wB - 2.0 * mag * wY

    Xs = np.vstack([Xs0, Xs1]).astype(np.float32)
    Ys = np.concatenate([np.zeros(n // 2), np.ones(n // 2)]).astype(int)
    Xt = np.vstack([Xt0, Xt1]).astype(np.float32)
    Yt = np.concatenate([np.zeros(nt0), np.ones(nt1)]).astype(int)
    return Xs, Ys, Xt, Yt


def auc_flip(a):
    return max(a, 1.0 - a)


def evaluate(scores, Yt):
    auc = float(roc_auc_score(Yt, scores))
    dial = (auc_flip(auc) - 0.5) if auc < 0.5 else 0.0
    return auc, dial


# --------------------------------------------------------------------
# Methods
# --------------------------------------------------------------------
def m_noadapt(Xs, Ys, Xt, Yt):
    clf = LogisticRegression(max_iter=2000).fit(Xs, Ys)
    s = clf.predict_proba(Xt)[:, 1]
    return s


def m_combat_naive(Xs, Ys, Xt, Yt):
    """Naive ComBat: subtract per-batch mean, refit on combined."""
    Xs_c = Xs - Xs.mean(axis=0, keepdims=True)
    Xt_c = Xt - Xt.mean(axis=0, keepdims=True)
    clf = LogisticRegression(max_iter=2000).fit(Xs_c, Ys)
    s = clf.predict_proba(Xt_c)[:, 1]
    return s


def m_combat_subtype(Xs, Ys, Xt, Yt):
    """Subtype-preserving ComBat surrogate (uses target self-pseudolabels)."""
    # 1) Train initial classifier on source
    clf0 = LogisticRegression(max_iter=2000).fit(Xs, Ys)
    Yt_hat = clf0.predict(Xt)
    # 2) Center within (batch, class)
    def _center(X, Y):
        X = X.copy()
        for y in np.unique(Y):
            X[Y == y] = X[Y == y] - X[Y == y].mean(axis=0, keepdims=True)
        return X
    Xs_c = _center(Xs, Ys); Xt_c = _center(Xt, Yt_hat)
    clf = LogisticRegression(max_iter=2000).fit(Xs_c, Ys)
    s = clf.predict_proba(Xt_c)[:, 1]
    return s


def m_coral(Xs, Ys, Xt, Yt):
    """CORAL: align target covariance to source covariance."""
    eps = 1e-3
    Cs = np.cov(Xs.T) + eps * np.eye(Xs.shape[1])
    Ct = np.cov(Xt.T) + eps * np.eye(Xt.shape[1])
    Es, Vs = np.linalg.eigh(Cs); Et, Vt = np.linalg.eigh(Ct)
    Cs_half = Vs @ np.diag(np.sqrt(np.maximum(Es, eps))) @ Vs.T
    Ct_inv_half = Vt @ np.diag(1.0 / np.sqrt(np.maximum(Et, eps))) @ Vt.T
    A = Ct_inv_half @ Cs_half
    Xt_aligned = (Xt - Xt.mean(axis=0)) @ A + Xs.mean(axis=0)
    clf = LogisticRegression(max_iter=2000).fit(Xs, Ys)
    s = clf.predict_proba(Xt_aligned)[:, 1]
    return s


def m_bn_stats(Xs, Ys, Xt, Yt):
    """Schneider 2020: replace train BN stats with target stats."""
    mu_s, sg_s = Xs.mean(0), Xs.std(0) + 1e-6
    mu_t, sg_t = Xt.mean(0), Xt.std(0) + 1e-6
    Xs_n = (Xs - mu_s) / sg_s
    Xt_n = (Xt - mu_t) / sg_t
    clf = LogisticRegression(max_iter=2000).fit(Xs_n, Ys)
    s = clf.predict_proba(Xt_n)[:, 1]
    return s


# ---- Neural backbone for TENT / DANN / SHOT ----
class MLPHead(nn.Module):
    def __init__(self, p, h=64, k=2):
        super().__init__()
        self.fe = nn.Sequential(nn.Linear(p, h), nn.BatchNorm1d(h), nn.ReLU(),
                                nn.Linear(h, h), nn.BatchNorm1d(h), nn.ReLU())
        self.cls = nn.Linear(h, k)

    def forward(self, x):
        return self.cls(self.fe(x))


def _train_src(Xs, Ys, p, epochs=80, lr=1e-3, seed=0):
    torch.manual_seed(seed)
    net = MLPHead(p)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    Xt = torch.tensor(Xs); Yt = torch.tensor(Ys)
    net.train()
    for _ in range(epochs):
        opt.zero_grad()
        logits = net(Xt)
        loss = F.cross_entropy(logits, Yt)
        loss.backward(); opt.step()
    return net


def _scores(net, X):
    net.eval()
    with torch.no_grad():
        logit = net(torch.tensor(X))
        return F.softmax(logit, dim=1)[:, 1].numpy()


def m_neural_baseline(Xs, Ys, Xt, Yt):
    p = Xs.shape[1]
    net = _train_src(Xs, Ys, p)
    return _scores(net, Xt)


def m_tent(Xs, Ys, Xt, Yt, epochs=20):
    p = Xs.shape[1]
    net = _train_src(Xs, Ys, p)
    # TENT: only update affine BN params, minimise entropy on target
    for m in net.modules():
        if isinstance(m, nn.BatchNorm1d):
            m.requires_grad_(True)
            m.track_running_stats = False
        elif isinstance(m, (nn.Linear,)):
            m.weight.requires_grad_(False)
            if m.bias is not None: m.bias.requires_grad_(False)
    params = [p for p in net.parameters() if p.requires_grad]
    opt = torch.optim.Adam(params, lr=1e-3)
    Xtt = torch.tensor(Xt)
    net.train()
    for _ in range(epochs):
        opt.zero_grad()
        logit = net(Xtt)
        prob = F.softmax(logit, dim=1)
        ent = -(prob * torch.log(prob + 1e-8)).sum(dim=1).mean()
        ent.backward(); opt.step()
    return _scores(net, Xt)


def m_shot(Xs, Ys, Xt, Yt, epochs=30):
    """SHOT-lite: source-free adaptation via pseudo-label IM loss."""
    p = Xs.shape[1]
    net = _train_src(Xs, Ys, p)
    # Freeze classifier, adapt feature extractor
    for q in net.cls.parameters():
        q.requires_grad_(False)
    opt = torch.optim.Adam([q for q in net.parameters() if q.requires_grad], lr=5e-4)
    Xtt = torch.tensor(Xt)
    net.train()
    for it in range(epochs):
        with torch.no_grad():
            pseudo = net(Xtt).argmax(dim=1)
        opt.zero_grad()
        logit = net(Xtt)
        prob = F.softmax(logit, dim=1)
        # IM = entropy minimisation + class-balance maximisation
        ent = -(prob * torch.log(prob + 1e-8)).sum(dim=1).mean()
        mean_p = prob.mean(0)
        div = (mean_p * torch.log(mean_p + 1e-8)).sum()  # negative entropy of marginal
        # cross-entropy w/ pseudo-labels
        ce = F.cross_entropy(logit, pseudo)
        loss = ent + div + 0.3 * ce
        loss.backward(); opt.step()
    return _scores(net, Xt)


def m_dann_lite(Xs, Ys, Xt, Yt, epochs=80):
    """DANN-lite: adversarial domain confusion via gradient-reversal."""
    p = Xs.shape[1]
    torch.manual_seed(0)
    fe = nn.Sequential(nn.Linear(p, 64), nn.BatchNorm1d(64), nn.ReLU(),
                       nn.Linear(64, 32), nn.BatchNorm1d(32), nn.ReLU())
    cls = nn.Linear(32, 2); dom = nn.Linear(32, 2)
    opt = torch.optim.Adam(list(fe.parameters()) + list(cls.parameters()) + list(dom.parameters()),
                           lr=1e-3)
    Xs_t = torch.tensor(Xs); Ys_t = torch.tensor(Ys)
    Xt_t = torch.tensor(Xt)
    Xall = torch.cat([Xs_t, Xt_t]); Dall = torch.cat([
        torch.zeros(len(Xs), dtype=torch.long),
        torch.ones(len(Xt), dtype=torch.long)])
    for ep in range(epochs):
        opt.zero_grad()
        feats_s = fe(Xs_t); feats_all = fe(Xall)
        loss_cls = F.cross_entropy(cls(feats_s), Ys_t)
        # gradient reversal: maximise domain loss
        lam = 0.5 * (2.0 / (1 + np.exp(-10 * ep / epochs)) - 1.0)
        loss_dom = -lam * F.cross_entropy(dom(feats_all), Dall)
        (loss_cls + loss_dom).backward(); opt.step()
    fe.eval(); cls.eval()
    with torch.no_grad():
        return F.softmax(cls(fe(Xt_t)), dim=1)[:, 1].numpy()


def m_foundation_proxy(Xs, Ys, Xt, Yt):
    """Foundation-model TTA proxy: PCA(combined) → MLP(source) → TENT."""
    pca = PCA(n_components=min(32, Xs.shape[1])).fit(np.vstack([Xs, Xt]))
    Xs_e = pca.transform(Xs).astype(np.float32)
    Xt_e = pca.transform(Xt).astype(np.float32)
    return m_tent(Xs_e, Ys, Xt_e, Yt)


# --------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------
def main():
    methods = [
        ("no_adapt", "classical", m_noadapt),
        ("combat_naive", "classical", m_combat_naive),
        ("combat_subtype", "classical", m_combat_subtype),
        ("CORAL", "classical", m_coral),
        ("BN_stats", "tta", m_bn_stats),
        ("DANN_lite", "adversarial", m_dann_lite),
        ("neural_baseline", "neural", m_neural_baseline),
        ("TENT", "tta", m_tent),
        ("SHOT_lite", "tta", m_shot),
        ("foundation_TTA", "foundation", m_foundation_proxy),
    ]

    rows = []
    n_seeds = 5
    for seed in range(n_seeds):
        Xs, Ys, Xt, Yt = make_thca_like(seed=seed)
        for name, klass, fn in methods:
            t0 = time.time()
            try:
                s = fn(Xs, Ys, Xt, Yt)
                auc, dial = evaluate(s, Yt)
                err = ""
            except Exception as e:
                auc, dial, err = float("nan"), float("nan"), str(e)[:80]
            rows.append(dict(method=name, klass=klass, seed=seed,
                             auc_target=auc, dial=dial, time_s=time.time() - t0,
                             error=err))
            print(f"  seed={seed} {name:18s}  AUC={auc:.3f}  DIAL={dial:.3f}  t={time.time()-t0:.1f}s")

    df = pd.DataFrame(rows)
    df.to_csv(RES / "methods_comparison_raw.tsv", sep="\t", index=False)

    summ = df.groupby(["method", "klass"]).agg(
        AUC_mean=("auc_target", "mean"),
        AUC_std=("auc_target", "std"),
        DIAL_mean=("dial", "mean"),
        DIAL_std=("dial", "std"),
        time_mean=("time_s", "mean"),
    ).round(3).reset_index()
    summ = summ.sort_values("AUC_mean", ascending=False)
    summ.to_csv(RES / "methods_comparison.tsv", sep="\t", index=False)
    print("\n=== Method comparison ===")
    print(summ.to_string(index=False))

    # Bar chart
    palette = {"classical": "#94a3b8", "adversarial": "#a78bfa",
               "tta": "#34d399", "foundation": "#f59e0b", "neural": "#60a5fa"}
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=summ.method, y=summ.AUC_mean,
        marker=dict(color=[palette[k] for k in summ.klass]),
        error_y=dict(type="data", array=summ.AUC_std),
        text=[f"{v:.2f}" for v in summ.AUC_mean],
        textposition="outside",
        name="AUC_target",
    ))
    fig.add_trace(go.Bar(
        x=summ.method, y=summ.DIAL_mean,
        marker=dict(color="rgba(244,63,94,0.65)"),
        error_y=dict(type="data", array=summ.DIAL_std),
        name="DIAL",
        yaxis="y2",
    ))
    fig.update_layout(
        title=dict(text="TTA benchmark — target AUC (bars) and DIAL (red, right axis)",
                   x=0.01, font=dict(color="#e2e8f0", size=18)),
        height=520, paper_bgcolor="#0b1220", plot_bgcolor="#0b1220",
        font=dict(color="#e2e8f0"),
        yaxis=dict(title="AUC_target", range=[0, 1.05],
                   gridcolor="rgba(255,255,255,0.1)"),
        yaxis2=dict(title="DIAL", overlaying="y", side="right",
                    range=[0, 0.5], gridcolor="rgba(255,255,255,0.05)"),
        xaxis=dict(tickangle=-30, gridcolor="rgba(255,255,255,0.1)"),
        barmode="group",
        legend=dict(x=0.78, y=1.0, bgcolor="rgba(11,18,32,0.6)"),
        margin=dict(l=60, r=60, t=80, b=120),
    )
    out_html = FIGS / "tta_benchmark.html"
    fig.write_html(out_html, include_plotlyjs="cdn", config=dict(displayModeBar=False))
    print(f"Wrote {out_html}")

    (CKPT / "task3_tta.json").write_text(json.dumps(dict(
        n_methods=len(methods), n_seeds=n_seeds, n_rows=len(df),
        best=summ.iloc[0].to_dict(),
    ), indent=2, default=str))


if __name__ == "__main__":
    main()
