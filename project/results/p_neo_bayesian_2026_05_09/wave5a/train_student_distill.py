"""Wave 5A — Step 3+4: Train student MLP distilled from MHCflurry teacher.

Encoder: peptide biophys (24-D) + HLA pseudo-seq one-hot (34-pos × 21-aa = 714-D)
  → 738-D input. (Wave 1 ESM2-Bayesian re-embedding would have taken ~25 min on
   shared CPU; for Wave 5A's distillation-uplift question the encoder is held
   constant within Wave 5A — comparison is *teacher-on vs teacher-off* with the
   SAME student architecture and SAME features.)

Loss = α · BCE(student, hard_label_if_in_train_pool)
       + β · KL_T(student, mhcflurry_presentation)         [main distill]
       + γ · KL_T(student_proc, mhcflurry_processing)     [aux]
       + δ · KL_T(student_aff,  mhcflurry_aff_score)      [aux]

Multi-task: 4 heads on the same trunk:
  H1 — sigmoid: hard immunogenicity (eval head)
  H2 — sigmoid: MHCflurry presentation
  H3 — sigmoid: MHCflurry processing
  H4 — sigmoid: MHCflurry affinity (transformed)

Default α=0.3, β=0.4, γ=0.15, δ=0.15. Temperature T=2.0 for soft labels.

Ablation: --no_distill turns off heads H2–H4 and trains only H1 on the 1962
hard-labeled rows in the pool — this is essentially Wave 1 with biophys features
(an honest baseline that isolates the distillation uplift from any encoder
difference).

MC Dropout T=20 inference for predictive uncertainty / ECE.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss

torch.set_num_threads(4)  # share with Wave 5B (~8 cores total)

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
W5A = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/wave5a")
ROOT = W5A.parent
BUNDLE = ROOT / "bundle.tsv"
HLA_PSEUDO_TSV = ROOT / "hla_pseudo.tsv"
POOL_TSV = W5A / "distill_pool.tsv"
TEACHER_TSV = W5A / "distill_teacher_scores.tsv"

# Encoder constants
KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,
      "H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,
      "T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}
PEP_PAD = 11   # max peptide length we encode (8-11mers)
PEP_FEATS = 24
PSEUDO_LEN = 34
HLA_FEAT_DIM = PSEUDO_LEN * (len(AA) + 1)  # +1 for unknown


# ---------------------------------------------------------------------------
# Featurization
# ---------------------------------------------------------------------------
def biophys_one(p):
    """24-D biophys vector for one peptide."""
    p = str(p).upper()
    L = max(len(p), 1)
    v = np.zeros(PEP_FEATS, dtype=np.float32)
    v[0] = sum(KD.get(a, 0) for a in p) / L
    v[1] = p.count("R") + p.count("K") + p.count("H") * 0.1 - p.count("D") - p.count("E")
    v[2] = (p.count("F") + p.count("W") + p.count("Y")) / L
    v[3] = float(L)
    for j, aa in enumerate(AA):
        v[4 + j] = p.count(aa) / L
    return v


def pep_kmer_oh(p, max_len=PEP_PAD):
    """Position-aware one-hot: pad to max_len, 21-dim per position (20 AA + UNK)."""
    p = str(p).upper()[:max_len]
    M = np.zeros((max_len, len(AA) + 1), dtype=np.float32)
    for i, a in enumerate(p):
        M[i, AA_IDX.get(a, len(AA))] = 1.0
    return M.flatten()  # max_len * 21 = 231-D


def hla_oh(pseudo):
    """Position-wise one-hot of 34-aa pseudo: 34 * 21 = 714-D."""
    s = str(pseudo).upper()[:PSEUDO_LEN].ljust(PSEUDO_LEN, "X")
    M = np.zeros((PSEUDO_LEN, len(AA) + 1), dtype=np.float32)
    for i, a in enumerate(s):
        M[i, AA_IDX.get(a, len(AA))] = 1.0
    return M.flatten()


def encode(pep_list, hla_list, hla_pseudo_map):
    """Returns X = [biophys(24) + pep_kmer_oh(231) + hla_oh(714)] = 969-D."""
    n = len(pep_list)
    X = np.zeros((n, PEP_FEATS + PEP_PAD * (len(AA) + 1) + HLA_FEAT_DIM), dtype=np.float32)
    for i in range(n):
        bp = biophys_one(pep_list[i])
        pk = pep_kmer_oh(pep_list[i])
        h = hla_pseudo_map.get(hla_list[i], "X" * PSEUDO_LEN)
        ho = hla_oh(h)
        X[i] = np.concatenate([bp, pk, ho])
    return X


# ---------------------------------------------------------------------------
# Model — Bayesian multitask MLP with always-on dropout (MC Dropout)
# ---------------------------------------------------------------------------
class MultiTaskMLP(nn.Module):
    def __init__(self, in_dim, hid=256, p_drop=0.3, n_heads=4):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hid)
        self.bn1 = nn.BatchNorm1d(hid)
        self.fc2 = nn.Linear(hid, hid)
        self.bn2 = nn.BatchNorm1d(hid)
        self.heads = nn.ModuleList([nn.Linear(hid, 1) for _ in range(n_heads)])
        self.p = p_drop

    def features(self, x, mc_dropout=True):
        h = F.relu(self.bn1(self.fc1(x)))
        # always-on dropout for MC Dropout
        h = F.dropout(h, p=self.p, training=mc_dropout or self.training)
        h = F.relu(self.bn2(self.fc2(h)))
        h = F.dropout(h, p=self.p, training=mc_dropout or self.training)
        return h

    def forward(self, x, mc_dropout=True):
        h = self.features(x, mc_dropout=mc_dropout)
        return torch.cat([head(h) for head in self.heads], dim=-1)  # [B, 4]


# ---------------------------------------------------------------------------
# Loss
# ---------------------------------------------------------------------------
def soft_kl_logits(student_logit, teacher_prob, T=2.0):
    """KL(teacher || sigmoid(student/T)) as binary distillation loss.
    teacher_prob: [B] in [0,1]; student_logit: [B] real."""
    s = student_logit / T
    # Bernoulli KL: t·log(t/s) + (1-t)·log((1-t)/(1-s)) ; we use temperature-scaled BCE
    t = teacher_prob.clamp(1e-6, 1 - 1e-6)
    return F.binary_cross_entropy_with_logits(s, t, reduction="mean") * (T * T)


def focal_bce(logit, target, gamma=2.0, label_smooth=0.05, mask=None):
    target = target * (1 - label_smooth) + 0.5 * label_smooth
    p = torch.sigmoid(logit)
    ce = F.binary_cross_entropy_with_logits(logit, target, reduction="none")
    pt = p * target + (1 - p) * (1 - target)
    L = (1 - pt) ** gamma * ce
    if mask is not None:
        L = L * mask.float()
        denom = mask.float().sum().clamp(min=1)
        return L.sum() / denom
    return L.mean()


# ---------------------------------------------------------------------------
# Eval helpers (bootstrap CI95, ECE)
# ---------------------------------------------------------------------------
def ece_score(y, p, n_bins=15):
    y, p = np.asarray(y), np.asarray(p)
    bins = np.linspace(0, 1, n_bins + 1)
    e = 0.0
    n = len(y)
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (p >= lo) & (p < hi if hi < 1 else p <= hi)
        if m.sum() == 0:
            continue
        e += (m.sum() / n) * abs(p[m].mean() - y[m].mean())
    return float(e)


def auroc_ci(y, p, n_boot=1000, seed=0):
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
    auc, lo, hi = auroc_ci(y, p, n_boot=n_boot)
    return {
        "AUROC": auc, "AUROC_lo95": lo, "AUROC_hi95": hi,
        "AUPRC": float(average_precision_score(y, p)),
        "Brier": float(brier_score_loss(y, p)),
        "ECE": ece_score(y, p),
        "n": int(len(y)), "n_pos": int(np.sum(y)),
    }


# ---------------------------------------------------------------------------
# Train
# ---------------------------------------------------------------------------
def train_student(X, y_hard, mask_hard, y_pres, y_proc, y_aff,
                  src, n_sources, *,
                  alpha=0.3, beta=0.4, gamma=0.15, delta=0.15, T=2.0,
                  epochs=30, lr=3e-4, batch=128,
                  use_distill=True, seed=0):
    torch.manual_seed(seed); np.random.seed(seed)
    in_dim = X.shape[1]
    model = MultiTaskMLP(in_dim).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)

    X_t = torch.from_numpy(X).float()
    y_hard_t = torch.from_numpy(y_hard).float()
    mask_t = torch.from_numpy(mask_hard).float()
    y_pres_t = torch.from_numpy(y_pres).float()
    y_proc_t = torch.from_numpy(y_proc).float()
    y_aff_t = torch.from_numpy(y_aff).float()
    src_t = torch.from_numpy(src).long()

    # source-balanced sampling weights
    src_counts = np.bincount(src, minlength=n_sources)
    src_counts = np.maximum(src_counts, 1)
    w = 1.0 / src_counts[src]
    w = w * len(w) / w.sum()
    sampler = WeightedRandomSampler(weights=w.tolist(), num_samples=len(w),
                                    replacement=True)
    ds = torch.utils.data.TensorDataset(X_t, y_hard_t, mask_t, y_pres_t, y_proc_t, y_aff_t, src_t)
    loader = DataLoader(ds, batch_size=batch, sampler=sampler, num_workers=0)

    for ep in range(epochs):
        model.train()
        ep_losses = {"hard": [], "pres": [], "proc": [], "aff": [], "total": []}
        for xb, yhb, mb, ypb, yprb, yab, _sb in loader:
            xb = xb.to(DEVICE); yhb = yhb.to(DEVICE); mb = mb.to(DEVICE)
            ypb = ypb.to(DEVICE); yprb = yprb.to(DEVICE); yab = yab.to(DEVICE)
            logits = model(xb, mc_dropout=True)
            l1, l2, l3, l4 = logits[:, 0], logits[:, 1], logits[:, 2], logits[:, 3]
            # H1: hard label, masked
            l_hard = focal_bce(l1, yhb, mask=mb)
            if use_distill:
                l_pres = soft_kl_logits(l2, ypb, T=T)
                l_proc = soft_kl_logits(l3, yprb, T=T)
                l_aff = soft_kl_logits(l4, yab, T=T)
                total = alpha * l_hard + beta * l_pres + gamma * l_proc + delta * l_aff
            else:
                # ablation: only hard-label H1
                l_pres = torch.tensor(0.0); l_proc = torch.tensor(0.0); l_aff = torch.tensor(0.0)
                total = l_hard
            opt.zero_grad()
            total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            ep_losses["hard"].append(float(l_hard))
            ep_losses["pres"].append(float(l_pres))
            ep_losses["proc"].append(float(l_proc))
            ep_losses["aff"].append(float(l_aff))
            ep_losses["total"].append(float(total))
        sched.step()
        if (ep + 1) % 5 == 0 or ep == 0 or ep == epochs - 1:
            mean = lambda k: np.mean(ep_losses[k])
            print(f"  ep{ep+1:>2d}  total={mean('total'):.4f}  hard={mean('hard'):.4f}  "
                  f"pres={mean('pres'):.4f}  proc={mean('proc'):.4f}  aff={mean('aff'):.4f}",
                  flush=True)
    return model


@torch.no_grad()
def predict_mc(model, X, n_samples=20, batch=2048):
    """MC Dropout T-sample inference. Returns (mean_prob_h1, std_prob_h1, all_heads_mean)."""
    model.eval()  # batchnorm in eval; dropout still on via mc_dropout=True flag
    X_t = torch.from_numpy(X).float().to(DEVICE)
    n = X_t.shape[0]
    samples = np.zeros((n_samples, n, 4), dtype=np.float32)
    for t in range(n_samples):
        outs = []
        for i in range(0, n, batch):
            xb = X_t[i:i+batch]
            logits = model(xb, mc_dropout=True)
            outs.append(torch.sigmoid(logits).cpu().numpy())
        samples[t] = np.concatenate(outs, axis=0)
    mean = samples.mean(0)
    std = samples.std(0)
    return mean, std


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--alpha", type=float, default=0.3)
    ap.add_argument("--beta", type=float, default=0.4)
    ap.add_argument("--gamma", type=float, default=0.15)
    ap.add_argument("--delta", type=float, default=0.15)
    ap.add_argument("--T", type=float, default=2.0)
    ap.add_argument("--no_distill", action="store_true",
                    help="ablation: train H1 only on hard labels (no MHCflurry)")
    ap.add_argument("--tag", type=str, default="distill",
                    help="output suffix")
    ap.add_argument("--seeds", type=int, default=3, help="ensemble seeds")
    args = ap.parse_args()

    print(f"=== Wave 5A — Train student (tag={args.tag}, distill={not args.no_distill}) ===")
    print(f"  epochs={args.epochs}  alpha={args.alpha}  beta={args.beta}  gamma={args.gamma}  delta={args.delta}  T={args.T}")
    t_start = time.time()

    # Load data
    teacher = pd.read_csv(TEACHER_TSV, sep="\t")
    pool = pd.read_csv(POOL_TSV, sep="\t")
    bundle = pd.read_csv(BUNDLE, sep="\t")
    hla_pseu = pd.read_csv(HLA_PSEUDO_TSV, sep="\t")
    hla_pseudo_map = dict(zip(hla_pseu["HLA_norm"], hla_pseu["pseudo_seq"]))
    # supplement with NetMHCpan-style pseudo for any MHCflurry alleles missing — fall back to "X"*34
    for h in teacher["hla"].unique():
        if h not in hla_pseudo_map or len(str(hla_pseudo_map.get(h, ""))) != 34:
            hla_pseudo_map[h] = "X" * 34  # gracefully degrade; HLA is still represented in distill via teacher score

    print(f"  teacher: {len(teacher):,} rows  in_train_pool: {teacher['in_train_pool'].sum():,}")

    # Encode features
    print("  featurizing pool...")
    t0 = time.time()
    X_pool = encode(teacher["peptide"].tolist(), teacher["hla"].tolist(), hla_pseudo_map)
    print(f"    X_pool {X_pool.shape}  in {time.time()-t0:.1f}s")

    # Hard labels (only meaningful for in_train_pool rows; mask others)
    y_hard = teacher["label_hard"].astype(int).values.astype(np.float32)
    mask_hard = teacher["in_train_pool"].astype(bool).values.astype(np.float32)
    # Soft labels
    y_pres = teacher["mhcflurry_presentation"].fillna(0.5).values.astype(np.float32)
    y_proc = teacher["mhcflurry_processing"].fillna(0.5).values.astype(np.float32)
    y_aff = teacher["mhcflurry_aff_score"].fillna(0.5).values.astype(np.float32)
    src_codes = pd.Categorical(teacher["source"]).codes
    src = src_codes.astype(np.int64)
    n_sources = len(set(src.tolist()))

    # Train ensemble
    models = []
    for s in range(args.seeds):
        print(f"\n[seed {s+1}/{args.seeds}] training...")
        m = train_student(
            X_pool, y_hard, mask_hard, y_pres, y_proc, y_aff,
            src, n_sources,
            alpha=args.alpha, beta=args.beta, gamma=args.gamma, delta=args.delta,
            T=args.T, epochs=args.epochs, batch=128,
            use_distill=not args.no_distill, seed=s,
        )
        models.append(m)

    # Save model state
    state = {f"seed_{i}": m.state_dict() for i, m in enumerate(models)}
    torch.save(state, W5A / f"student_model_{args.tag}.pt")

    # ----- Eval -----
    print("\n=== Eval ===")
    # Build eval featurization for ITSNdb
    itsn = bundle[bundle["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy()
    itsn = itsn[itsn["HLA_norm"].astype(str).str.len() > 0].copy()
    itsn = itsn[itsn["peptide"].astype(str).str.len().between(8, 11)].copy()
    # featurize — fall back to X*34 if HLA missing
    for h in itsn["HLA_norm"].unique():
        if h not in hla_pseudo_map:
            hla_pseudo_map[h] = "X" * 34
    X_itsn = encode(itsn["peptide"].tolist(), itsn["HLA_norm"].tolist(), hla_pseudo_map)
    print(f"  ITSNdb eval encoded: {X_itsn.shape}  pos_rate={itsn['label'].astype(int).mean():.3f}")

    # Ensemble MC Dropout predictions: average across seeds × samples
    all_preds_h1 = []
    for m in models:
        mean, std = predict_mc(m, X_itsn, n_samples=20)
        all_preds_h1.append(mean[:, 0])  # head 1
    ens = np.stack(all_preds_h1, axis=0).mean(0)
    itsn["score_distill"] = ens

    # In-master flag is on bundle
    itsn["in_master"] = itsn.get("in_master", False).astype(bool)

    summary = []
    for sub_name, mask, in_master_flag in [
        ("ITSNdb_main", itsn["split"] == "ext_itsndb_main", "mixed"),
        ("ITSNdb_Val",  itsn["split"] == "ext_itsndb_val", "mixed"),
        ("ITSNdb_combined", np.ones(len(itsn), dtype=bool), "mixed"),
        ("ITSNdb_no_overlap", ~itsn["in_master"].astype(bool).values, False),
        ("ITSNdb_in_master", itsn["in_master"].astype(bool).values, True),
    ]:
        sub = itsn[mask]
        if len(sub) >= 5 and sub["label"].astype(int).nunique() == 2:
            m = metrics_block(sub["label"].astype(int).values, sub["score_distill"].values, n_boot=1000)
            summary.append({"testset": sub_name, "in_master": in_master_flag, "tag": args.tag, **m})
            print(f"  {sub_name:>22s}: n={m['n']} pos={m['n_pos']} AUROC={m['AUROC']:.3f} "
                  f"[{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]  ECE={m['ECE']:.3f}")

    # Per-allele LOSO-ish (top 11 alleles)
    print("\n  per-allele AUROC (top 11):")
    top_alleles = itsn["HLA_norm"].value_counts().head(11).index.tolist()
    per_allele = []
    for h in top_alleles:
        sub = itsn[itsn["HLA_norm"] == h]
        if len(sub) < 8 or sub["label"].astype(int).nunique() < 2:
            continue
        m = metrics_block(sub["label"].astype(int).values, sub["score_distill"].values, n_boot=500)
        per_allele.append({"hla": h, **m})
        print(f"    {h:>14s} n={m['n']} pos={m['n_pos']} AUROC={m['AUROC']:.3f} [{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]")

    # Save predictions + summary
    pred_out = itsn[["peptide", "HLA_norm", "label", "in_master", "split", "score_distill"]].copy()
    pred_out.columns = ["peptide", "hla", "label", "in_master", "split", "score_distill"]
    pred_out.to_csv(W5A / f"predictions_wave5a_{args.tag}.tsv", sep="\t", index=False)

    pd.DataFrame(summary).to_csv(W5A / f"wave5a_results_{args.tag}.tsv", sep="\t", index=False)
    pd.DataFrame(per_allele).to_csv(W5A / f"wave5a_per_allele_{args.tag}.tsv", sep="\t", index=False)
    print(f"\nsaved: predictions / results / per_allele tables (tag={args.tag})")
    print(f"  total wall: {time.time()-t_start:.1f}s")


if __name__ == "__main__":
    main()
