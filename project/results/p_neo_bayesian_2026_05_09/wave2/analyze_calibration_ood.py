"""Wave 2 Track B: Calibration + OOD detection on Wave 1 Bayesian predictions.

Inputs (from wave1/):
  - predictions_in_domain.tsv  (cols: y, p_bayes, p_ens, p_std)
  - predictions_itsndb.tsv     (cols: peptide, HLA_norm, label, source, split,
                                       protein_id, in_master, pred_mean, pred_std)

Outputs (wave2/):
  - calibration_results.tsv          per-subset ECE / Brier / NLL
  - ood_detection_results.tsv        OOD detection AUROC (entropy & std as scores)
  - selective_prediction_curve.tsv   AUROC vs % retained
  - fig_calibration_panel.png/pdf    reliability diagrams (3 subplots)
  - fig_ood_detection.png/pdf        entropy density + ROC
  - fig_selective_prediction.png/pdf AUROC vs % retained
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


HERE = Path(__file__).resolve().parent
WAVE1 = HERE.parent / "wave1"
WAVE2 = HERE


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

EPS = 1e-9


def predictive_entropy(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, EPS, 1 - EPS)
    return -(p * np.log(p) + (1.0 - p) * np.log(1.0 - p))


def expected_calibration_error(y: np.ndarray, p: np.ndarray, n_bins: int = 15) -> float:
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    inds = np.digitize(p, bins) - 1
    inds = np.clip(inds, 0, n_bins - 1)
    ece = 0.0
    n = len(y)
    for b in range(n_bins):
        mask = inds == b
        if mask.sum() == 0:
            continue
        conf = p[mask].mean()
        acc = y[mask].mean()
        ece += (mask.sum() / n) * abs(acc - conf)
    return float(ece)


def reliability_curve(y: np.ndarray, p: np.ndarray, n_bins: int = 15):
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    centers = 0.5 * (bins[1:] + bins[:-1])
    inds = np.digitize(p, bins) - 1
    inds = np.clip(inds, 0, n_bins - 1)
    accs = np.full(n_bins, np.nan)
    confs = np.full(n_bins, np.nan)
    counts = np.zeros(n_bins, dtype=int)
    for b in range(n_bins):
        mask = inds == b
        if mask.sum() > 0:
            accs[b] = float(y[mask].mean())
            confs[b] = float(p[mask].mean())
            counts[b] = int(mask.sum())
    return centers, confs, accs, counts


def brier(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.mean((p - y) ** 2))


def nll(y: np.ndarray, p: np.ndarray) -> float:
    p = np.clip(p, EPS, 1 - EPS)
    return float(-np.mean(y * np.log(p) + (1 - y) * np.log(1 - p)))


def safe_auroc(y: np.ndarray, p: np.ndarray) -> float:
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, p))


def selective_curve(y: np.ndarray, p: np.ndarray, conf: np.ndarray):
    order = np.argsort(-conf)  # high conf first
    y_sorted = y[order]
    p_sorted = p[order]
    fracs = []
    aurocs = []
    for k in range(20, len(y) + 1):
        yk = y_sorted[:k]
        pk = p_sorted[:k]
        if len(np.unique(yk)) < 2:
            continue
        fracs.append(k / len(y))
        aurocs.append(roc_auc_score(yk, pk))
    return np.array(fracs), np.array(aurocs)


# ---------------------------------------------------------------------------
# Load Wave 1 predictions
# ---------------------------------------------------------------------------

print("[load] reading wave1 predictions")
in_dom = pd.read_csv(WAVE1 / "predictions_in_domain.tsv", sep="\t")
itsn = pd.read_csv(WAVE1 / "predictions_itsndb.tsv", sep="\t")
print(f"  in_domain rows={len(in_dom)} cols={list(in_dom.columns)}")
print(f"  itsndb rows={len(itsn)} cols={list(itsn.columns)}")

# In-domain: 5-fold concatenated. p_ens is the ensemble (better calibrated).
in_dom_y = in_dom["y"].astype(int).to_numpy()
in_dom_p = in_dom["p_ens"].astype(float).to_numpy()
in_dom_pstd = in_dom["p_std"].astype(float).to_numpy()

# ITSNdb subsets
itsn_y = itsn["label"].astype(int).to_numpy()
itsn_p = itsn["pred_mean"].astype(float).to_numpy()
itsn_pstd = itsn["pred_std"].astype(float).to_numpy()
no_overlap_mask = ~itsn["in_master"].astype(bool).to_numpy()

# ---------------------------------------------------------------------------
# Calibration: per-subset ECE / Brier / NLL + reliability data
# ---------------------------------------------------------------------------

print("[calib] computing ECE / Brier / NLL per subset")
calib_rows = []

subsets = [
    ("in_domain_kfold", in_dom_y, in_dom_p, in_dom_pstd),
    ("ITSNdb_combined", itsn_y, itsn_p, itsn_pstd),
    ("ITSNdb_main", itsn_y[itsn["split"] == "ext_itsndb_main"], itsn_p[itsn["split"] == "ext_itsndb_main"], itsn_pstd[itsn["split"] == "ext_itsndb_main"]),
    ("ITSNdb_Val", itsn_y[itsn["split"] == "ext_itsndb_val"], itsn_p[itsn["split"] == "ext_itsndb_val"], itsn_pstd[itsn["split"] == "ext_itsndb_val"]),
    ("ITSNdb_no_overlap", itsn_y[no_overlap_mask], itsn_p[no_overlap_mask], itsn_pstd[no_overlap_mask]),
    ("ITSNdb_in_master", itsn_y[~no_overlap_mask], itsn_p[~no_overlap_mask], itsn_pstd[~no_overlap_mask]),
]

reliability_data = {}
for name, y, p, pstd in subsets:
    if len(y) == 0:
        continue
    ece = expected_calibration_error(y, p, n_bins=15)
    bri = brier(y, p)
    nl = nll(y, p)
    auc = safe_auroc(y, p)
    ent = predictive_entropy(p)
    calib_rows.append(
        dict(
            subset=name,
            n=int(len(y)),
            n_pos=int(y.sum()),
            AUROC=auc,
            ECE=ece,
            Brier=bri,
            NLL=nl,
            mean_pred=float(p.mean()),
            mean_entropy=float(ent.mean()),
            mean_p_std=float(pstd.mean()),
        )
    )
    reliability_data[name] = reliability_curve(y, p, n_bins=15)

calib_df = pd.DataFrame(calib_rows)
calib_df.to_csv(WAVE2 / "calibration_results.tsv", sep="\t", index=False)
print(calib_df.to_string(index=False))

# ---------------------------------------------------------------------------
# OOD detection: in_domain (label 0) vs ITSNdb_no_overlap (label 1)
#   Score 1: predictive entropy (uncertainty in mean p)
#   Score 2: posterior std (epistemic proxy from MC dropout T=30)
# ---------------------------------------------------------------------------

print("[ood] OOD detection AUROC (entropy + std as scores)")

# Build paired ID vs OOD set
ood_p = itsn_p[no_overlap_mask]
ood_pstd = itsn_pstd[no_overlap_mask]
ood_ent = predictive_entropy(ood_p)

id_p = in_dom_p
id_pstd = in_dom_pstd
id_ent = predictive_entropy(id_p)

# OOD detection labels: 0 = ID, 1 = OOD
all_ent = np.concatenate([id_ent, ood_ent])
all_pstd = np.concatenate([id_pstd, ood_pstd])
ood_lab = np.concatenate([np.zeros(len(id_ent)), np.ones(len(ood_ent))]).astype(int)

ood_auroc_ent = safe_auroc(ood_lab, all_ent)
ood_auroc_std = safe_auroc(ood_lab, all_pstd)

# Also compare ID vs ITSNdb_in_master (less OOD) and ID vs ITSNdb_combined
in_master_mask = ~no_overlap_mask
mast_ent = predictive_entropy(itsn_p[in_master_mask])
mast_pstd = itsn_pstd[in_master_mask]
ood_lab_m = np.concatenate([np.zeros(len(id_ent)), np.ones(len(mast_ent))]).astype(int)
ood_auroc_ent_m = safe_auroc(ood_lab_m, np.concatenate([id_ent, mast_ent]))
ood_auroc_std_m = safe_auroc(ood_lab_m, np.concatenate([id_pstd, mast_pstd]))

# Combined ITSNdb as "external"
all_ent_c = np.concatenate([id_ent, predictive_entropy(itsn_p)])
all_pstd_c = np.concatenate([id_pstd, itsn_pstd])
ood_lab_c = np.concatenate([np.zeros(len(id_ent)), np.ones(len(itsn_p))]).astype(int)
ood_auroc_ent_c = safe_auroc(ood_lab_c, all_ent_c)
ood_auroc_std_c = safe_auroc(ood_lab_c, all_pstd_c)

ood_rows = [
    dict(comparison="ID_vs_ITSNdb_no_overlap", n_id=len(id_ent), n_ood=len(ood_ent),
         AUROC_entropy=ood_auroc_ent, AUROC_pstd=ood_auroc_std),
    dict(comparison="ID_vs_ITSNdb_in_master", n_id=len(id_ent), n_ood=len(mast_ent),
         AUROC_entropy=ood_auroc_ent_m, AUROC_pstd=ood_auroc_std_m),
    dict(comparison="ID_vs_ITSNdb_combined", n_id=len(id_ent), n_ood=len(itsn_p),
         AUROC_entropy=ood_auroc_ent_c, AUROC_pstd=ood_auroc_std_c),
]
ood_df = pd.DataFrame(ood_rows)
ood_df.to_csv(WAVE2 / "ood_detection_results.tsv", sep="\t", index=False)
print(ood_df.to_string(index=False))

# ---------------------------------------------------------------------------
# Selective prediction on the headline external test (ITSNdb_no_overlap)
#   Confidence = 1 - entropy_normalized; OR confidence = -p_std (lower std = more confident)
# ---------------------------------------------------------------------------

print("[selective] selective prediction on ITSNdb_no_overlap")
y_no = itsn_y[no_overlap_mask]
p_no = itsn_p[no_overlap_mask]
ps_no = itsn_pstd[no_overlap_mask]
conf_ent_no = -predictive_entropy(p_no)  # high = confident
conf_std_no = -ps_no                      # high = confident

frac_e, auc_e = selective_curve(y_no, p_no, conf_ent_no)
frac_s, auc_s = selective_curve(y_no, p_no, conf_std_no)

# Also do for ITSNdb_combined (more positives, smoother curve)
conf_ent_c = -predictive_entropy(itsn_p)
conf_std_c = -itsn_pstd
frac_ec, auc_ec = selective_curve(itsn_y, itsn_p, conf_ent_c)
frac_sc, auc_sc = selective_curve(itsn_y, itsn_p, conf_std_c)

# Save selective curve table
sel_rows = []
for name, frac, auc in [
    ("ITSNdb_no_overlap_entropy", frac_e, auc_e),
    ("ITSNdb_no_overlap_pstd", frac_s, auc_s),
    ("ITSNdb_combined_entropy", frac_ec, auc_ec),
    ("ITSNdb_combined_pstd", frac_sc, auc_sc),
]:
    for f, a in zip(frac, auc):
        sel_rows.append(dict(curve=name, frac_retained=float(f), AUROC=float(a)))

sel_df = pd.DataFrame(sel_rows)
sel_df.to_csv(WAVE2 / "selective_prediction_curve.tsv", sep="\t", index=False)


def find_threshold(frac, auc, target=0.80):
    """Return min frac such that AUROC >= target. Return NaN if never reached."""
    hits = np.where(auc >= target)[0]
    if len(hits) == 0:
        return float("nan"), float(np.nanmax(auc)) if len(auc) else float("nan")
    idx = hits[-1]  # largest retained frac that still meets threshold
    return float(frac[idx]), float(auc[idx])


thr_e_no, auc_e_no = find_threshold(frac_e, auc_e, 0.80)
thr_s_no, auc_s_no = find_threshold(frac_s, auc_s, 0.80)
thr_e_c, auc_e_c = find_threshold(frac_ec, auc_ec, 0.80)
thr_s_c, auc_s_c = find_threshold(frac_sc, auc_sc, 0.80)

print(f"  no_overlap entropy : K=AUROC>=0.80 retained at frac={thr_e_no} (max AUROC={auc_e_no:.3f})")
print(f"  no_overlap pstd    : K=AUROC>=0.80 retained at frac={thr_s_no} (max AUROC={auc_s_no:.3f})")
print(f"  combined  entropy  : K=AUROC>=0.80 retained at frac={thr_e_c} (max AUROC={auc_e_c:.3f})")
print(f"  combined  pstd     : K=AUROC>=0.80 retained at frac={thr_s_c} (max AUROC={auc_s_c:.3f})")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

print("[fig] reliability diagrams")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
plot_subsets = ["in_domain_kfold", "ITSNdb_combined", "ITSNdb_no_overlap"]
for ax, name in zip(axes, plot_subsets):
    if name not in reliability_data:
        continue
    centers, confs, accs, counts = reliability_data[name]
    valid = counts > 0
    ax.plot([0, 1], [0, 1], color="grey", linestyle="--", lw=1, label="perfect")
    ax.bar(centers[valid], accs[valid], width=1.0 / 15, alpha=0.55, edgecolor="black",
           color="#1f77b4", label="empirical")
    ax.scatter(confs[valid], accs[valid], color="darkred", zorder=5, s=18)
    row = calib_df[calib_df["subset"] == name].iloc[0]
    ax.set_title(f"{name}\nECE={row['ECE']:.3f}  Brier={row['Brier']:.3f}  AUROC={row['AUROC']:.3f}",
                 fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Predicted prob")
    ax.set_ylabel("Empirical prob")
    ax.grid(alpha=0.3)
fig.suptitle("Wave 2 Track B — Reliability diagrams (Wave 1 Bayesian predictions)", fontsize=11)
fig.savefig(WAVE2 / "fig_calibration_panel.png", dpi=160)
fig.savefig(WAVE2 / "fig_calibration_panel.pdf")
plt.close(fig)

print("[fig] OOD entropy density + ROC")
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
ax = axes[0]
ax.hist(id_ent, bins=30, alpha=0.55, color="#1f77b4", label=f"in-domain (n={len(id_ent)})", density=True)
ax.hist(predictive_entropy(itsn_p[in_master_mask]), bins=30, alpha=0.55, color="#ff7f0e",
        label=f"ITSNdb_in_master (n={in_master_mask.sum()})", density=True)
ax.hist(ood_ent, bins=20, alpha=0.55, color="#d62728", label=f"ITSNdb_no_overlap (n={len(ood_ent)})", density=True)
ax.set_xlabel("Predictive entropy H(p)")
ax.set_ylabel("Density")
ax.set_title("Predictive entropy distribution")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

ax = axes[1]
from sklearn.metrics import roc_curve

for lbl, scores, lab, color in [
    ("entropy: ID vs no_overlap", all_ent, ood_lab, "#d62728"),
    ("p_std: ID vs no_overlap", all_pstd, ood_lab, "#1f77b4"),
    ("entropy: ID vs combined", all_ent_c, ood_lab_c, "#ff7f0e"),
    ("p_std: ID vs combined", all_pstd_c, ood_lab_c, "#2ca02c"),
]:
    if len(np.unique(lab)) < 2:
        continue
    fpr, tpr, _ = roc_curve(lab, scores)
    auc = roc_auc_score(lab, scores)
    ax.plot(fpr, tpr, label=f"{lbl}  AUROC={auc:.3f}", color=color)
ax.plot([0, 1], [0, 1], "--", color="grey", lw=1)
ax.set_xlabel("FPR")
ax.set_ylabel("TPR")
ax.set_title("OOD detection ROC")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.suptitle("Wave 2 Track B — OOD detection (predictive entropy & posterior std)", fontsize=11)
fig.savefig(WAVE2 / "fig_ood_detection.png", dpi=160)
fig.savefig(WAVE2 / "fig_ood_detection.pdf")
plt.close(fig)

print("[fig] selective prediction")
fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
ax.plot(frac_e, auc_e, color="#d62728", label="no_overlap, entropy")
ax.plot(frac_s, auc_s, color="#1f77b4", label="no_overlap, p_std")
ax.plot(frac_ec, auc_ec, color="#ff7f0e", label="combined, entropy")
ax.plot(frac_sc, auc_sc, color="#2ca02c", label="combined, p_std")
ax.axhline(0.80, ls="--", color="grey", lw=1, label="AUROC=0.80 (deploy threshold)")
ax.set_xlabel("Fraction retained (most confident)")
ax.set_ylabel("AUROC on retained subset")
ax.set_title("Selective prediction — Wave 1 Bayesian on ITSNdb")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3)
fig.savefig(WAVE2 / "fig_selective_prediction.png", dpi=160)
fig.savefig(WAVE2 / "fig_selective_prediction.pdf")
plt.close(fig)

# ---------------------------------------------------------------------------
# Save selective summary
# ---------------------------------------------------------------------------

summary = dict(
    selective_no_overlap_entropy=dict(retained_at_auc080=thr_e_no, max_auc_at_threshold=auc_e_no),
    selective_no_overlap_pstd=dict(retained_at_auc080=thr_s_no, max_auc_at_threshold=auc_s_no),
    selective_combined_entropy=dict(retained_at_auc080=thr_e_c, max_auc_at_threshold=auc_e_c),
    selective_combined_pstd=dict(retained_at_auc080=thr_s_c, max_auc_at_threshold=auc_s_c),
    ood_AUROC_entropy_no_overlap=float(ood_auroc_ent),
    ood_AUROC_pstd_no_overlap=float(ood_auroc_std),
    ood_AUROC_entropy_combined=float(ood_auroc_ent_c),
    ood_AUROC_pstd_combined=float(ood_auroc_std_c),
)
(WAVE2 / "track_b_summary.json").write_text(json.dumps(summary, indent=2))

print("[done] Track B")
print(json.dumps(summary, indent=2))
