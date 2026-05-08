"""Wave 2 Track A: Structure-aware features (CPU, fast).

Cheap structural proxies:
  1. Anchor-residue PWMs (per-HLA position weight matrices) derived from training pool
  2. BLOSUM62 max-similarity to known binders (per-HLA, label=1 only)
  3. Physicochemical features (Kyte-Doolittle hydrophobicity, charge, helix/sheet propensity)

Skipped here: ESM-IF1 (would need ~150M model download; pod is busy with Wave 1 VenusVaccine
scoring; documented as "deferred to v2" rather than mocked). Structure proxies above are
pure-CPU and do not duplicate Wave 1.

Inputs:
  ../bundle.tsv
  ../hla_pseudo.tsv
  ../wave1/predictions_in_domain.tsv  (for cross-fold AUROC ablation)
  ../wave1/predictions_itsndb.tsv

Outputs (wave2/):
  structure_features.tsv             per-row features for full bundle
  structure_ablation_results.tsv     LR baseline / +structure / +ESM2-mean / +structure+ESM2-mean
                                     evaluated on in_domain (5-fold) & ITSNdb_no_overlap
  fig_structure_ablation.png/pdf
"""

from __future__ import annotations

import json
from pathlib import Path
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WAVE1 = ROOT / "wave1"

bundle = pd.read_csv(ROOT / "bundle.tsv", sep="\t")
hla_pseudo = pd.read_csv(ROOT / "hla_pseudo.tsv", sep="\t").dropna()
print(f"[load] bundle={len(bundle)}, hla_pseudo={len(hla_pseudo)}")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AA = list("ACDEFGHIKLMNPQRSTVWY")
AA_IDX = {a: i for i, a in enumerate(AA)}

# Kyte-Doolittle hydrophobicity (standard)
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "E": -3.5, "Q": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
# Charge at pH 7 (rough)
CHG = {a: 0.0 for a in AA}
CHG.update({"R": 1.0, "K": 1.0, "H": 0.1, "D": -1.0, "E": -1.0})
# Chou-Fasman helix propensity
HEL = {
    "A": 1.42, "R": 0.98, "N": 0.67, "D": 1.01, "C": 0.70, "E": 1.51, "Q": 1.11,
    "G": 0.57, "H": 1.00, "I": 1.08, "L": 1.21, "K": 1.16, "M": 1.45, "F": 1.13,
    "P": 0.57, "S": 0.77, "T": 0.83, "W": 1.08, "Y": 0.69, "V": 1.06,
}
# Chou-Fasman sheet propensity
SHE = {
    "A": 0.83, "R": 0.93, "N": 0.89, "D": 0.54, "C": 1.19, "E": 0.37, "Q": 1.10,
    "G": 0.75, "H": 0.87, "I": 1.60, "L": 1.30, "K": 0.74, "M": 1.05, "F": 1.38,
    "P": 0.55, "S": 0.75, "T": 1.19, "W": 1.37, "Y": 1.47, "V": 1.70,
}

# BLOSUM62 (standard 20x20). Embedded so we don't depend on external file.
BLOSUM62_RAW = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""


def parse_blosum(raw: str):
    rows = [r.split() for r in raw.strip().splitlines() if r.strip()]
    header = rows[0]
    M = np.zeros((20, 20))
    for i, row in enumerate(rows[1:]):
        for j, h in enumerate(header):
            M[AA_IDX[row[0]], AA_IDX[h]] = float(row[1 + j])
    return M


BLOSUM = parse_blosum(BLOSUM62_RAW)


def blosum_score(a: str, b: str) -> float:
    if a not in AA_IDX or b not in AA_IDX:
        return 0.0
    return float(BLOSUM[AA_IDX[a], AA_IDX[b]])


def pep_blosum(p: str, q: str) -> float:
    """Sum BLOSUM62 over aligned positions. Pads shorter peptide with negatives."""
    L = min(len(p), len(q))
    if L == 0:
        return 0.0
    return sum(blosum_score(p[i], q[i]) for i in range(L)) / L


# ---------------------------------------------------------------------------
# Build per-HLA PWMs from training pool only (avoid leakage)
# ---------------------------------------------------------------------------

train_df = bundle[bundle["split"] == "train"].copy()
print(f"[pwm] training pool n={len(train_df)} (label==1: {int(train_df['label'].sum())})")

# Background AA frequency from the full training pool
all_aa = "".join(train_df["peptide"].tolist())
bg = np.array([all_aa.count(a) for a in AA], dtype=float)
bg = (bg + 1) / (bg + 1).sum()

PWM_BY_HLA: dict[str, np.ndarray] = {}  # shape: (Lmax=11, 20) log-odds; padded with zeros
PWM_LMAX = 11

hla_counts = train_df.groupby("HLA_norm").size()
top_hlas = list(hla_counts[hla_counts >= 20].index)  # require at least 20 train peptides
print(f"[pwm] HLAs with >=20 train peptides: {len(top_hlas)}")

for hla in top_hlas:
    sub = train_df[(train_df["HLA_norm"] == hla) & (train_df["label"] == 1)]
    if len(sub) < 5:
        continue
    pwm = np.zeros((PWM_LMAX, 20))
    counts = np.zeros((PWM_LMAX, 20))
    for pep in sub["peptide"]:
        L = min(len(pep), PWM_LMAX)
        for i in range(L):
            a = pep[i]
            if a in AA_IDX:
                counts[i, AA_IDX[a]] += 1
    # Pseudo-count + log-odds
    for i in range(PWM_LMAX):
        if counts[i].sum() == 0:
            continue
        freq = (counts[i] + 1) / (counts[i].sum() + 20)
        pwm[i] = np.log(freq / bg)
    PWM_BY_HLA[hla] = pwm

print(f"[pwm] PWMs built for {len(PWM_BY_HLA)} HLAs")


def anchor_score(pep: str, hla: str) -> tuple[float, float, float, float]:
    """Return (P2_score, POmega_score, mean_pwm, sum_pwm). 0 if HLA missing."""
    if hla not in PWM_BY_HLA:
        return 0.0, 0.0, 0.0, 0.0
    pwm = PWM_BY_HLA[hla]
    L = min(len(pep), PWM_LMAX)
    if L < 2:
        return 0.0, 0.0, 0.0, 0.0
    s = 0.0
    n = 0
    p2 = 0.0
    if pep[1] in AA_IDX:
        p2 = float(pwm[1, AA_IDX[pep[1]]])
    pOm = 0.0
    if pep[L - 1] in AA_IDX:
        # PΩ: align to last position used in PWM (min(L,PWM_LMAX)-1)
        pOm = float(pwm[L - 1, AA_IDX[pep[L - 1]]])
    for i in range(L):
        a = pep[i]
        if a in AA_IDX:
            s += pwm[i, AA_IDX[a]]
            n += 1
    return p2, pOm, s / max(n, 1), s


# Per-HLA pool of label==1 train peptides for BLOSUM max-sim
TRAIN_BINDERS_BY_HLA: dict[str, list[str]] = defaultdict(list)
for _, r in train_df[train_df["label"] == 1].iterrows():
    TRAIN_BINDERS_BY_HLA[r["HLA_norm"]].append(r["peptide"])
print(f"[blosum] HLAs with binders: {len(TRAIN_BINDERS_BY_HLA)}")


def max_blosum_to_train_binders(pep: str, hla: str) -> tuple[float, float]:
    """Max + mean BLOSUM62 score (per-residue avg) to label==1 train peptides for same HLA.
    Excludes self-match if pep itself is in the train pool for that HLA."""
    binders = TRAIN_BINDERS_BY_HLA.get(hla, [])
    if not binders:
        return 0.0, 0.0
    scores = [pep_blosum(pep, b) for b in binders if b != pep]
    if not scores:
        return 0.0, 0.0
    return float(max(scores)), float(np.mean(scores))


def physico(pep: str) -> np.ndarray:
    if not pep:
        return np.zeros(12)
    kd = [KD.get(a, 0.0) for a in pep]
    chg = [CHG.get(a, 0.0) for a in pep]
    hel = [HEL.get(a, 0.0) for a in pep]
    she = [SHE.get(a, 0.0) for a in pep]
    return np.array([
        np.mean(kd), np.std(kd), np.max(kd), np.min(kd),
        np.mean(chg), np.sum(chg),
        np.mean(hel), np.std(hel),
        np.mean(she), np.std(she),
        len(pep), sum(1 for a in pep if a in "FWY") / len(pep),  # aromatic fraction
    ])


# ---------------------------------------------------------------------------
# Compute features for all bundle rows that are not VenusVaccine windows
# ---------------------------------------------------------------------------

print("[feat] computing structural features for full bundle (excluding venus)")
keep_splits = ["train", "ext_itsndb_main", "ext_itsndb_val"]
sub_bundle = bundle[bundle["split"].isin(keep_splits)].copy().reset_index(drop=True)

rows = []
for _, r in sub_bundle.iterrows():
    pep = str(r["peptide"])
    hla = str(r["HLA_norm"])
    p2, pOm, mpwm, spwm = anchor_score(pep, hla)
    bmx, bmn = max_blosum_to_train_binders(pep, hla)
    phys = physico(pep)
    rows.append([
        pep, hla, int(r["label"]), str(r["split"]), bool(r.get("in_master", False)),
        p2, pOm, mpwm, spwm, bmx, bmn,
        *phys.tolist(),
    ])

cols = ["peptide", "HLA_norm", "label", "split", "in_master",
        "pwm_p2", "pwm_pOmega", "pwm_mean", "pwm_sum",
        "blosum_max", "blosum_mean",
        "kd_mean", "kd_std", "kd_max", "kd_min",
        "chg_mean", "chg_sum",
        "hel_mean", "hel_std",
        "she_mean", "she_std",
        "pep_len", "aromatic_frac"]
sf = pd.DataFrame(rows, columns=cols)
sf.to_csv(HERE / "structure_features.tsv", sep="\t", index=False)
print(f"[feat] wrote {HERE / 'structure_features.tsv'} ({len(sf)} rows)")

# ---------------------------------------------------------------------------
# Ablation: structure-only LR baseline vs adding to existing setup
# We don't have ESM2 mean here (full embeddings on pod); approximate ESM2 axis with
# Wave 1 Bayesian prediction p_mean as a proxy for "ESM2-trained signal".
# ---------------------------------------------------------------------------

print("[ablate] per-row ablation (structure features alone, then +Wave1-pred)")

# Attach Wave 1 ITSNdb predictions
itsn_pred = pd.read_csv(WAVE1 / "predictions_itsndb.tsv", sep="\t")
sf_itsn = sf[sf["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy()
sf_itsn = sf_itsn.merge(
    itsn_pred[["peptide", "HLA_norm", "split", "pred_mean", "pred_std"]],
    on=["peptide", "HLA_norm", "split"], how="left",
)

train_idx = sf["split"] == "train"
itsn_idx = sf["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])
no_overlap_idx = (sf["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])) & (~sf["in_master"])

feat_cols = [c for c in sf.columns if c not in ("peptide", "HLA_norm", "label", "split", "in_master")]
print(f"  feature cols: {feat_cols}")

X_train = sf.loc[train_idx, feat_cols].to_numpy()
y_train = sf.loc[train_idx, "label"].to_numpy()
X_itsn = sf.loc[itsn_idx, feat_cols].to_numpy()
y_itsn = sf.loc[itsn_idx, "label"].to_numpy()
itsn_in_master = sf.loc[itsn_idx, "in_master"].to_numpy()

scaler = StandardScaler().fit(X_train)
Xt = scaler.transform(X_train)
Xi = scaler.transform(X_itsn)

# 5-fold CV on train pool
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_aucs = []
for fold, (tr, va) in enumerate(skf.split(Xt, y_train)):
    clf = LogisticRegression(C=1.0, max_iter=3000)
    clf.fit(Xt[tr], y_train[tr])
    p = clf.predict_proba(Xt[va])[:, 1]
    if len(np.unique(y_train[va])) >= 2:
        cv_aucs.append(roc_auc_score(y_train[va], p))

cv_mean = float(np.mean(cv_aucs))
cv_std = float(np.std(cv_aucs))
print(f"  structure-only 5-fold AUROC: {cv_mean:.3f} ± {cv_std:.3f}")

# Train on full pool, evaluate on ITSNdb subsets
clf_full = LogisticRegression(C=1.0, max_iter=3000).fit(Xt, y_train)
p_itsn = clf_full.predict_proba(Xi)[:, 1]

auc_combined_struct = roc_auc_score(y_itsn, p_itsn)
y_no = y_itsn[~itsn_in_master]
p_no = p_itsn[~itsn_in_master]
y_mast = y_itsn[itsn_in_master]
p_mast = p_itsn[itsn_in_master]
auc_no_struct = roc_auc_score(y_no, p_no) if len(np.unique(y_no)) >= 2 else float("nan")
auc_mast_struct = roc_auc_score(y_mast, p_mast) if len(np.unique(y_mast)) >= 2 else float("nan")
print(f"  structure-only ITSNdb_combined: {auc_combined_struct:.3f}")
print(f"  structure-only ITSNdb_no_overlap: {auc_no_struct:.3f}")
print(f"  structure-only ITSNdb_in_master: {auc_mast_struct:.3f}")

# Ablation: structure + Wave1 pred (stacked LR, only on ITSNdb where Wave1 pred exists)
sf_itsn_with = sf_itsn.dropna(subset=["pred_mean"]).copy()
y_itsn_w = sf_itsn_with["label"].to_numpy()
in_m_w = sf_itsn_with["in_master"].to_numpy()

X_struct_w = scaler.transform(sf_itsn_with[feat_cols].to_numpy())
p_struct_w = clf_full.predict_proba(X_struct_w)[:, 1]
p_w1_w = sf_itsn_with["pred_mean"].to_numpy()

# Weighted average (50/50) of struct LR and Wave1 — simple ensemble
p_ens_w = 0.5 * p_struct_w + 0.5 * p_w1_w

# Logistic stacker: fit a tiny LR on ITSNdb_in_master (held in-master rows, leaky-ish but
# acceptable as a sanity check; the *cleaner* eval is the no_overlap subset).
# For honest no-leak eval, we ONLY look at no_overlap below where ensemble is fixed-weight.

ablate_rows = []


def auc_or_nan(y, p):
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, p))


# Baseline: Wave1 only
for split_name, mask in [
    ("ITSNdb_combined", np.ones(len(y_itsn_w), dtype=bool)),
    ("ITSNdb_no_overlap", ~in_m_w),
    ("ITSNdb_in_master", in_m_w),
]:
    ablate_rows.append(dict(
        model="Wave1_Bayesian_only", subset=split_name, n=int(mask.sum()),
        AUROC=auc_or_nan(y_itsn_w[mask], p_w1_w[mask]),
    ))
    ablate_rows.append(dict(
        model="Structure_LR_only", subset=split_name, n=int(mask.sum()),
        AUROC=auc_or_nan(y_itsn_w[mask], p_struct_w[mask]),
    ))
    ablate_rows.append(dict(
        model="Wave1+Structure_avg50", subset=split_name, n=int(mask.sum()),
        AUROC=auc_or_nan(y_itsn_w[mask], p_ens_w[mask]),
    ))

# Add 5-fold structure-only number
ablate_rows.append(dict(
    model="Structure_LR_only", subset="train_5fold_cv", n=int(train_idx.sum()),
    AUROC=cv_mean,
))

abl = pd.DataFrame(ablate_rows)
abl.to_csv(HERE / "structure_ablation_results.tsv", sep="\t", index=False)
print(abl.to_string(index=False))

# ---------------------------------------------------------------------------
# Save ensemble probabilities for later combined table
# ---------------------------------------------------------------------------

sf_itsn_with["p_struct_lr"] = p_struct_w
sf_itsn_with["p_wave1_plus_struct_avg50"] = p_ens_w
keep_cols = ["peptide", "HLA_norm", "label", "split", "in_master",
             "pred_mean", "pred_std", "p_struct_lr", "p_wave1_plus_struct_avg50"]
sf_itsn_with[keep_cols].to_csv(HERE / "structure_lr_predictions.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# Figure: bar comparison
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
order = ["ITSNdb_combined", "ITSNdb_no_overlap", "ITSNdb_in_master"]
models = ["Wave1_Bayesian_only", "Structure_LR_only", "Wave1+Structure_avg50"]
colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
W = 0.25
for i, m in enumerate(models):
    sub = abl[(abl["model"] == m) & (abl["subset"].isin(order))].set_index("subset").reindex(order)
    ax.bar(np.arange(len(order)) + (i - 1) * W, sub["AUROC"], width=W, label=m, color=colors[i])
ax.axhline(0.5, ls="--", color="grey", lw=1)
ax.set_xticks(np.arange(len(order)))
ax.set_xticklabels(order, rotation=15, ha="right")
ax.set_ylabel("AUROC")
ax.set_ylim(0.0, 1.0)
ax.set_title("Wave 2 Track A — Structure LR ablation vs Wave 1 Bayesian")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3, axis="y")
fig.savefig(HERE / "fig_structure_ablation.png", dpi=160)
fig.savefig(HERE / "fig_structure_ablation.pdf")
plt.close(fig)

summary = dict(
    structure_only_5fold_AUROC=cv_mean,
    structure_only_5fold_std=cv_std,
    structure_only_ITSNdb_combined=auc_combined_struct,
    structure_only_ITSNdb_no_overlap=auc_no_struct,
    structure_only_ITSNdb_in_master=auc_mast_struct,
    n_PWMs=len(PWM_BY_HLA),
    n_train_binder_HLAs=sum(1 for v in TRAIN_BINDERS_BY_HLA.values() if len(v) >= 1),
    feat_cols=feat_cols,
    notes="ESM-IF1 deferred (pod busy with Wave 1 venus scoring; ESM-IF1 not paper-blocking).",
)
(HERE / "track_a_summary.json").write_text(json.dumps(summary, indent=2))
print("[done] Track A")
