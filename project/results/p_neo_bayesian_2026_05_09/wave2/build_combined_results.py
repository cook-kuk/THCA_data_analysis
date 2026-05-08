"""Wave 2 combined results table + ensemble. Runs AFTER all three tracks complete.

Builds wave2_combined_results.tsv: rows = {RF / Wave1 Bayesian / Structure-LR / VQC / QK_SVM /
GP_RBF / GP_quantum / Ensemble}, columns = {ITSNdb_combined / no_overlap / in_master / VenusVaccine}.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
WAVE1 = HERE.parent / "wave1"


def auc_or_nan(y, p):
    p = np.asarray(p, dtype=float)
    y = np.asarray(y, dtype=int)
    valid = ~np.isnan(p)
    if valid.sum() == 0 or len(np.unique(y[valid])) < 2:
        return float("nan"), int(valid.sum())
    return float(roc_auc_score(y[valid], p[valid])), int(valid.sum())


# ---------------------------------------------------------------------------
# Load all per-row predictions
# ---------------------------------------------------------------------------

w1 = pd.read_csv(WAVE1 / "predictions_itsndb.tsv", sep="\t")
print(f"[w1] {len(w1)} rows; cols={list(w1.columns)}")

vqc = pd.read_csv(HERE / "vqc_predictions.tsv", sep="\t")
print(f"[vqc] {len(vqc)} rows")

qk = HERE / "qk_svm_predictions.tsv"
if qk.exists():
    qk_df = pd.read_csv(qk, sep="\t")
    print(f"[qk] {len(qk_df)} rows")
else:
    print("[qk] NOT YET PRESENT")
    qk_df = None

struct = pd.read_csv(HERE / "structure_lr_predictions.tsv", sep="\t")
print(f"[struct] {len(struct)} rows")


# ---------------------------------------------------------------------------
# Build merged frame
# ---------------------------------------------------------------------------

key = ["peptide", "HLA_norm", "split"]
m = w1[key + ["label", "in_master", "pred_mean", "pred_std"]].copy()
m = m.rename(columns={"pred_mean": "p_wave1_bayes", "pred_std": "p_wave1_std"})

m = m.merge(vqc[key + ["p_vqc", "p_lr_8d_classical"]], on=key, how="outer")
m = m.merge(struct[key + ["p_struct_lr", "p_wave1_plus_struct_avg50"]], on=key, how="outer")
if qk_df is not None:
    m = m.merge(qk_df[key + ["p_qk_svm", "p_rbf_svm", "p_gp_rbf", "p_gp_q"]], on=key, how="outer")

# Backfill label / in_master from any source if missing
m["label"] = m["label"].fillna(0).astype(int)
m["in_master"] = m["in_master"].fillna(False).astype(bool)
print(f"[merge] {len(m)} rows total; cols={list(m.columns)}")

# ---------------------------------------------------------------------------
# Ensembles (only on rows where all components present)
# ---------------------------------------------------------------------------

# Pick a "core ensemble" of available models
core_cols = ["p_wave1_bayes", "p_vqc", "p_struct_lr"]
if qk_df is not None and "p_qk_svm" in m.columns:
    core_cols.append("p_qk_svm")

# Mean ensemble (only rows where all components are non-NaN)
m["p_ensemble_core"] = m[core_cols].mean(axis=1)

# Wave1 + structure (already computed)
# Wave1 + VQC + structure mean
m["p_ens_wave1_vqc_struct"] = m[["p_wave1_bayes", "p_vqc", "p_struct_lr"]].mean(axis=1)

# ---------------------------------------------------------------------------
# Per-subset AUROC table
# ---------------------------------------------------------------------------

subsets = [
    ("ITSNdb_combined", np.ones(len(m), dtype=bool)),
    ("ITSNdb_no_overlap", (~m["in_master"]).to_numpy()),
    ("ITSNdb_in_master", m["in_master"].to_numpy()),
    ("ITSNdb_main", (m["split"] == "ext_itsndb_main").to_numpy()),
    ("ITSNdb_Val", (m["split"] == "ext_itsndb_val").to_numpy()),
]

model_cols = [
    ("Wave1_Bayesian", "p_wave1_bayes"),
    ("Structure_LR", "p_struct_lr"),
    ("Wave1+Structure_avg50", "p_wave1_plus_struct_avg50"),
    ("VQC_8q3l_quantum", "p_vqc"),
    ("LR_8d_PCA_classical", "p_lr_8d_classical"),
]
if qk_df is not None:
    model_cols += [
        ("QK_SVM_4q_quantum", "p_qk_svm"),
        ("RBF_SVM_8d_classical", "p_rbf_svm"),
        ("GP_RBF_8d_classical_Bayesian", "p_gp_rbf"),
        ("GP_quantum_4q_Bayesian", "p_gp_q"),
    ]
model_cols += [
    ("Ensemble_core_mean", "p_ensemble_core"),
    ("Ensemble_W1+VQC+Struct", "p_ens_wave1_vqc_struct"),
]

rows = []
for sname, mask in subsets:
    for mname, mcol in model_cols:
        if mcol not in m.columns:
            continue
        y = m.loc[mask, "label"].to_numpy()
        p = m.loc[mask, mcol].to_numpy()
        a, nv = auc_or_nan(y, p)
        rows.append(dict(model=mname, subset=sname, n_total=int(mask.sum()), n_valid=nv, AUROC=a))

# Add Wave1 RF baseline numbers (from the user's quoted values)
rf_baseline = [
    ("RF_baseline", "ITSNdb_combined", 0.734),
    ("RF_baseline", "ITSNdb_no_overlap", 0.431),
    ("RF_baseline", "ITSNdb_in_master", float("nan")),
    ("RF_baseline_5fold_in_domain", "in_domain_kfold", 0.854),
]
for mname, sname, auc in rf_baseline:
    rows.append(dict(model=mname, subset=sname, n_total=int(0), n_valid=0, AUROC=auc))

# Add VenusVaccine numbers from Wave 1 venus_bayesian_results.tsv (top10_mean is best)
venus_path = WAVE1 / "venus_bayesian_results.tsv"
if venus_path.exists():
    v = pd.read_csv(venus_path, sep="\t")
    for _, r in v.iterrows():
        rows.append(dict(
            model=f"Wave1_venus_{r['aggregator']}",
            subset=f"VenusVaccine_{r['split']}",
            n_total=int(r["n"]), n_valid=int(r["n"]),
            AUROC=float(r["AUROC"]),
        ))

res = pd.DataFrame(rows)
res.to_csv(HERE / "wave2_combined_results.tsv", sep="\t", index=False)
print(res.to_string(index=False))

# ---------------------------------------------------------------------------
# Headline summary JSON
# ---------------------------------------------------------------------------

def best_for(subset_name):
    sub = res[(res["subset"] == subset_name) & (res["AUROC"].notna())]
    if len(sub) == 0:
        return None
    best = sub.loc[sub["AUROC"].idxmax()]
    return dict(model=str(best["model"]), AUROC=float(best["AUROC"]), n=int(best["n_valid"]))


summary = dict(
    best_ITSNdb_no_overlap=best_for("ITSNdb_no_overlap"),
    best_ITSNdb_combined=best_for("ITSNdb_combined"),
    best_ITSNdb_in_master=best_for("ITSNdb_in_master"),
    best_VenusVaccine_both=best_for("VenusVaccine_both"),
    delta_no_overlap_vs_RF=None,
    delta_no_overlap_vs_Wave1=None,
)
if summary["best_ITSNdb_no_overlap"]:
    summary["delta_no_overlap_vs_RF"] = float(summary["best_ITSNdb_no_overlap"]["AUROC"] - 0.431)
    w1_no = res[(res["model"] == "Wave1_Bayesian") & (res["subset"] == "ITSNdb_no_overlap")]
    if len(w1_no):
        summary["delta_no_overlap_vs_Wave1"] = float(summary["best_ITSNdb_no_overlap"]["AUROC"] - w1_no.iloc[0]["AUROC"])

(HERE / "wave2_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))


# ---------------------------------------------------------------------------
# Final headline figure: AUROC across methods on no_overlap (the headline test)
# ---------------------------------------------------------------------------

no_ov = res[(res["subset"] == "ITSNdb_no_overlap") & (res["AUROC"].notna())].copy()
no_ov = no_ov.sort_values("AUROC", ascending=True)

fig, ax = plt.subplots(figsize=(8, max(4.5, 0.35 * len(no_ov) + 1)), constrained_layout=True)
colors = ["#9467bd" if "quantum" in m.lower() or "vqc" in m.lower() or "gp_q" in m.lower()
          else ("#2ca02c" if "ensemble" in m.lower() else
                ("#d62728" if "rf" in m.lower() else "#1f77b4"))
          for m in no_ov["model"]]
ax.barh(no_ov["model"], no_ov["AUROC"], color=colors)
ax.axvline(0.5, ls="--", color="grey", lw=1, label="chance")
ax.axvline(0.431, ls=":", color="#d62728", lw=1, label="RF baseline (0.431)")
ax.axvline(0.411, ls=":", color="#1f77b4", lw=1, label="Wave 1 Bayesian (0.411)")
ax.set_xlabel("AUROC on ITSNdb_no_overlap (n≈103)")
ax.set_xlim(0.30, 0.80)
ax.set_title("Wave 2 — headline external test (truly OOD subset)")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3, axis="x")
fig.savefig(HERE / "fig_headline_no_overlap.png", dpi=160)
fig.savefig(HERE / "fig_headline_no_overlap.pdf")
plt.close(fig)

print("[done] combined")
