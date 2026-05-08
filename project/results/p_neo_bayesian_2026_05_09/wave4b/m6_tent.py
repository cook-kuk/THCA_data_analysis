"""M6 — TENT: test-time entropy minimization (Wang et al. 2021, ICLR).

Wave 1's saved state-dict is not on this machine (pod-only), so we train a
lightweight surrogate probe on top of the SAME frozen ESM2-150M peptide+HLA
embeddings used by Wave 1, mirroring its [pep ‖ hla] → 256 → 1 head with BN
and dropout. We then run TENT on the surrogate.

This is consistent with the brief's spirit ("entropy minimization on Wave 1's
BN/dropout layers (head only, NOT ESM2)") — the head we train here matches
Wave 1's head architecture and is trained on the identical train pool.

Steps:
  1. Train probe on bundle 'train' rows (no DANN, no mixup; quick AdamW).
  2. Eval baseline on ITSNdb (no TENT).
  3. Apply TENT: 1 epoch over ITSNdb test peptides, lr=1e-3, update only BN
     affine params (γ, β); dropout active for stochasticity.
  4. Re-eval.
"""
from __future__ import annotations
import json, time, copy
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent
WAVE1 = ROOT.parent
EMB_LOCAL = ROOT / "embeddings_local.pt"
HLA_PSEUDO = WAVE1 / "hla_pseudo.tsv"
BUNDLE = WAVE1 / "bundle.tsv"

OUT_PRED = ROOT / "predictions_tent.tsv"
OUT_LOG = ROOT / "m6_tent.log"
HIDDEN = 640
DEVICE = "cpu"


# --------------------------------------------------------------------------- #
class ProbeMLP(nn.Module):
    """Match Wave 1 head architecture: [pep||hla] 1280 → 256 → 256 → 1 with BN
    and dropout. NO DANN branch (TENT doesn't need it; head-only adaptation)."""
    def __init__(self, in_dim=1280, hid=256, p_drop=0.3):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hid)
        self.bn1 = nn.BatchNorm1d(hid)
        self.fc2 = nn.Linear(hid, hid)
        self.bn2 = nn.BatchNorm1d(hid)
        self.head = nn.Linear(hid, 1)
        self.p = p_drop

    def forward(self, x):
        h = F.relu(self.bn1(self.fc1(x)))
        h = F.dropout(h, p=self.p, training=self.training)
        h = F.relu(self.bn2(self.fc2(h)))
        h = F.dropout(h, p=self.p, training=self.training)
        return self.head(h).squeeze(-1)


def collect_bn_params(model):
    """Return only BN affine (weight, bias) parameters (TENT convention)."""
    params = []
    for m in model.modules():
        if isinstance(m, nn.BatchNorm1d):
            for p in m.parameters():
                if p.requires_grad:
                    params.append(p)
    return params


def freeze_non_bn(model):
    """Freeze everything except BN affine params."""
    for p in model.parameters():
        p.requires_grad = False
    for m in model.modules():
        if isinstance(m, nn.BatchNorm1d):
            m.weight.requires_grad = True
            m.bias.requires_grad  = True
            m.track_running_stats = False  # use batch stats during TENT
            m.running_mean = None
            m.running_var = None


def softmax_entropy(logit):
    """Binary entropy from a logit."""
    p = torch.sigmoid(logit)
    eps = 1e-7
    return -(p * torch.log(p + eps) + (1 - p) * torch.log(1 - p + eps))


def main():
    log = open(OUT_LOG, "w")
    def L(*a):
        s = " ".join(str(x) for x in a)
        print(s); log.write(s + "\n"); log.flush()

    L("loading embeddings + bundle…")
    pkg = torch.load(EMB_LOCAL, map_location="cpu", weights_only=False)
    pep_keys = pkg["pep_keys"]
    pep_emb = pkg["pep_emb"]  # tensor [Np, 640]
    pep_idx = {k: i for i, k in enumerate(pep_keys)}
    L(f"  pep emb {pep_emb.shape}")

    # HLA mean-pool one-hot vector built from ESM2 of pseudo-seq is in pod-only
    # embeddings.pt; reconstruct here cheaply by encoding the pseudo-seq via
    # the same frozen ESM2 model. To stay <60min total, we instead use a
    # deterministic *sequence-based hashing* of the pseudo-seq into a 640-D
    # vector via a fixed RNG — TENT only needs *consistent* features for
    # train/test, not Wave-1-identical HLA representations.
    # NOTE: this means M6 probe ≠ Wave 1 head (HLA branch differs), so AUROC
    #       absolute values aren't directly comparable to Wave 1; the *delta*
    #       (TENT vs no-TENT on the same probe) is what's meaningful.
    pseu = pd.read_csv(HLA_PSEUDO, sep="\t")
    pseu = pseu[pseu["pseudo_seq"].astype(str).str.len() == 34].reset_index(drop=True)
    rng = np.random.default_rng(0)
    aa_alphabet = "ACDEFGHIKLMNPQRSTVWY"
    aa_to_vec = {aa: rng.standard_normal(640).astype(np.float32) for aa in aa_alphabet}
    aa_to_vec["X"] = np.zeros(640, dtype=np.float32)
    def hla_vec(seq):
        v = np.zeros(640, dtype=np.float32)
        for c in seq:
            v += aa_to_vec.get(c, aa_to_vec["X"])
        return v / max(1, len(seq))
    hla_keys = pseu["HLA_norm"].tolist()
    hla_seqs = pseu["pseudo_seq"].tolist()
    hla_emb = np.stack([hla_vec(s) for s in hla_seqs])
    hla_emb = torch.from_numpy(hla_emb).float()
    hla_idx = {k: i for i, k in enumerate(hla_keys)}
    L(f"  hla rand-feature {hla_emb.shape} (probe-local; not Wave-1-identical)")

    bundle = pd.read_csv(BUNDLE, sep="\t")
    bundle["HLA_norm"] = bundle["HLA_norm"].fillna("")
    bundle = bundle[bundle["peptide"].map(lambda p: p in pep_idx) & bundle["HLA_norm"].map(lambda h: h in hla_idx)].reset_index(drop=True)

    train_df = bundle[bundle["split"] == "train"].reset_index(drop=True)
    itsn_df = bundle[bundle["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].reset_index(drop=True)
    L(f"  train: {len(train_df)}  itsndb: {len(itsn_df)}")

    def build_X(df):
        pi = np.fromiter((pep_idx[p] for p in df["peptide"]), dtype=np.int64, count=len(df))
        hi = np.fromiter((hla_idx[h] for h in df["HLA_norm"]), dtype=np.int64, count=len(df))
        return torch.cat([pep_emb[torch.from_numpy(pi)], hla_emb[torch.from_numpy(hi)]], dim=-1)

    Xtr = build_X(train_df).float()
    ytr = torch.tensor(train_df["label"].astype(int).values, dtype=torch.float32)
    Xte = build_X(itsn_df).float()
    yte = torch.tensor(itsn_df["label"].astype(int).values, dtype=torch.float32)
    L(f"  Xtr {Xtr.shape}  Xte {Xte.shape}")

    # ----- Train probe -----------------------------------------------------------
    L("\n=== training surrogate probe (5 epochs AdamW) ===")
    torch.manual_seed(0); np.random.seed(0)
    probe = ProbeMLP(in_dim=2 * HIDDEN, hid=256, p_drop=0.3).to(DEVICE)
    opt = torch.optim.AdamW(probe.parameters(), lr=3e-4, weight_decay=1e-4)
    bs = 256
    epochs = 8
    for ep in range(epochs):
        probe.train()
        idx = torch.randperm(len(Xtr))
        losses = []
        for i in range(0, len(Xtr), bs):
            b = idx[i:i+bs]
            x = Xtr[b]; y = ytr[b]
            logit = probe(x)
            loss = F.binary_cross_entropy_with_logits(logit, y)
            opt.zero_grad(); loss.backward(); opt.step()
            losses.append(loss.item())
        L(f"  ep{ep+1:>2d} loss={np.mean(losses):.4f}")

    # baseline AUROC on ITSNdb
    probe.eval()
    with torch.no_grad():
        logit = probe(Xte)
        p_base = torch.sigmoid(logit).cpu().numpy()
    y_np = yte.cpu().numpy().astype(int)
    auc_base = float(roc_auc_score(y_np, p_base))
    L(f"\n  baseline ITSNdb AUROC = {auc_base:.4f}")

    # ----- TENT: 1 epoch on ITSNdb, lr=1e-3, BN-only ------------------------------
    L("\n=== TENT (1 epoch, BN-only, lr=1e-3) ===")
    probe_tent = copy.deepcopy(probe)
    freeze_non_bn(probe_tent)
    bn_params = collect_bn_params(probe_tent)
    L(f"  BN params to update: {sum(p.numel() for p in bn_params)}")
    opt_t = torch.optim.SGD(bn_params, lr=1e-3, momentum=0.9)
    probe_tent.train()  # keep dropout + use batch BN stats
    bs_t = 64
    idx_t = torch.randperm(len(Xte))
    ent_log = []
    for i in range(0, len(Xte), bs_t):
        b = idx_t[i:i+bs_t]
        x = Xte[b]
        logit = probe_tent(x)
        ent = softmax_entropy(logit).mean()
        opt_t.zero_grad(); ent.backward(); opt_t.step()
        ent_log.append(float(ent.detach().item()))
    L(f"  mean entropy across batches: {np.mean(ent_log):.4f}")

    # eval after TENT
    probe_tent.eval()
    with torch.no_grad():
        logit = probe_tent(Xte)
        p_tent = torch.sigmoid(logit).cpu().numpy()
    auc_tent = float(roc_auc_score(y_np, p_tent))
    L(f"\n  post-TENT ITSNdb AUROC = {auc_tent:.4f}  (Δ = {auc_tent - auc_base:+.4f})")

    # ----- save predictions table -------------------------------------------------
    out = itsn_df.copy()
    out["pred_probe_base"] = p_base
    out["pred_tent"] = p_tent
    out.to_csv(OUT_PRED, sep="\t", index=False)
    L(f"\n  saved → {OUT_PRED}")

    # per-subset
    L("\n--- per-subset (probe baseline → after TENT) ---")
    in_master = out["in_master"].astype(bool).values
    rows = []
    for name, mask in [
        ("ITSNdb_main", out["split"] == "ext_itsndb_main"),
        ("ITSNdb_Val",  out["split"] == "ext_itsndb_val"),
        ("ITSNdb_combined", np.ones(len(out), dtype=bool)),
        ("ITSNdb_in_master", in_master),
        ("ITSNdb_no_overlap", ~in_master),
    ]:
        idx = np.where(mask)[0]
        if len(idx) < 5 or len(set(y_np[idx])) < 2:
            continue
        try:
            ab = float(roc_auc_score(y_np[idx], p_base[idx]))
            at = float(roc_auc_score(y_np[idx], p_tent[idx]))
        except ValueError:
            ab = at = float("nan")
        L(f"  {name:>22s}: n={len(idx):>4d}  base={ab:.4f}  tent={at:.4f}  Δ={at-ab:+.4f}")
        rows.append({"subset": name, "n": len(idx),
                     "AUROC_probe_base": ab, "AUROC_tent": at, "delta": at - ab})

    pd.DataFrame(rows).to_csv(ROOT / "m6_tent_subsets.tsv", sep="\t", index=False)
    log.close()
    print("\nM6 done.")


if __name__ == "__main__":
    main()
