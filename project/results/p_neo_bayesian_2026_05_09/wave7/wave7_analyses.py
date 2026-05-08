"""Wave 7 — Paper-grade analyses (shortcut tests, feature importance, OOD, counterfactual,
comparison vs Wave 1-6 baselines, hero figures).

Reads:
  combined_features.tsv
  predictions_w7a.tsv (also has p_w7b)
  wave7_results.tsv  (W7-A/B raw evals)

Pulls in baselines for the comparison table:
  Wave 1: predictions_itsndb.tsv (predictive_mean) + predictions_venus.tsv
  Wave 2: wave2_combined_results.tsv (Structure_LR, GP_quantum, ECE)
  Wave 3 MHCflurry: predictions.tsv
  Wave 3 BigMHC/PRIME/NetMHCpan/DeepImmuno/TransPHLA: per-folder predictions
  Wave 4: results json
  Wave 4C: ensemble file

Writes:
  wave7_results.tsv (extended with Wave1-6 rows)
  wave7_shortcut_tests.tsv
  wave7_feature_importance.tsv
  wave7_counterfactual_edits.tsv
  ood_results.tsv
  fig_wave7_hero_forest.{png,pdf}
  fig_wave7_calibration_panel.{png,pdf}
  fig_wave7_inflation_gap.{png,pdf}
  fig_wave7_feature_importance.{png,pdf}
  fig_wave7_per_allele_forest.{png,pdf}
"""
from __future__ import annotations
import os, sys, time, json, math, glob
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE7 = ROOT / "wave7"
COMB = WAVE7 / "combined_features.tsv"
PRED = WAVE7 / "predictions_w7a.tsv"
RESULTS = WAVE7 / "wave7_results.tsv"

t0 = time.time()
print("[load] combined_features + predictions...", flush=True)
df_feat = pd.read_csv(COMB, sep="\t", dtype={"in_master": "boolean"})
if "mhc_pct_log" in df_feat.columns and df_feat["mhc_pct_log"].isna().all():
    df_feat = df_feat.drop(columns=["mhc_pct_log"])
df_pred = pd.read_csv(PRED, sep="\t", dtype={"in_master": "boolean"})

FEATURE_COLS = [c for c in df_feat.columns if c not in ("peptide","HLA_norm","label","source","split","in_master")]
print(f"[features] {len(FEATURE_COLS)} cols", flush=True)

X = df_feat[FEATURE_COLS].fillna(0).astype(float).to_numpy()
y = df_feat["label"].astype(int).to_numpy()
src = df_feat["source"].to_numpy()
split = df_feat["split"].to_numpy()
in_master = df_feat["in_master"].fillna(False).astype(bool).to_numpy()
hla = df_feat["HLA_norm"].to_numpy()
train_mask = split == "train"

# Standardize once
scaler_full = StandardScaler().fit(X[train_mask])
Xs = scaler_full.transform(X)

# ---------------- Helper: bootstrap AUROC, ECE ----------------
def auroc_with_ci(y_, s_, n_boot=1000, seed=0):
    y_ = np.asarray(y_, int); s_ = np.asarray(s_, float)
    msk = np.isfinite(s_)
    y_ = y_[msk]; s_ = s_[msk]
    if len(np.unique(y_)) < 2: return (np.nan, np.nan, np.nan, len(y_))
    auc = roc_auc_score(y_, s_)
    rng = np.random.default_rng(seed)
    n = len(y_); aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_[idx])) < 2: continue
        aucs.append(roc_auc_score(y_[idx], s_[idx]))
    return (float(auc), float(np.percentile(aucs, 2.5)) if aucs else np.nan,
            float(np.percentile(aucs, 97.5)) if aucs else np.nan, n)

def ece_15bin(y_, p_, n_bins=15):
    y_ = np.asarray(y_, int); p_ = np.clip(np.asarray(p_, float), 1e-9, 1 - 1e-9)
    bins = np.linspace(0, 1, n_bins + 1)
    bin_idx = np.clip(np.digitize(p_, bins) - 1, 0, n_bins - 1)
    n = len(y_); ece = 0.0
    for b in range(n_bins):
        m = bin_idx == b
        if not m.any(): continue
        ece += (m.sum()/n) * abs(p_[m].mean() - y_[m].mean())
    return float(ece)

# ---------------- Subset masks ----------------
ITS_MASK = np.isin(split, ["ext_itsndb_main", "ext_itsndb_val"])
no_overlap_mask = ITS_MASK & (~in_master)
combined_mask = ITS_MASK
in_master_mask = ITS_MASK & in_master

# ---------------- Shortcut tests (DRP §3 analog) ----------------
print("\n[shortcut] running shortcut tests...", flush=True)
shortcut_rows = []

def train_lr_subset(mask_cols, name, train_mask=train_mask):
    Xsub = Xs[:, mask_cols]
    if Xsub.shape[1] == 0:
        return np.zeros(len(Xs))
    sc = StandardScaler().fit(Xsub[train_mask])
    Xss = sc.transform(Xsub)
    clf = LogisticRegression(max_iter=4000, C=1.0, class_weight="balanced", solver="lbfgs", random_state=0)
    clf.fit(Xss[train_mask], y[train_mask])
    return clf.predict_proba(Xss)[:, 1]

# Identify column groups
PEP_COLS = [c for c in FEATURE_COLS if c.startswith(("pwm_","kd_","chg_","hel_","she_","blosum_","aromatic","pep_len"))]
HLA_PCA_COLS = [c for c in FEATURE_COLS if c.startswith("hla_pc")]
MHC_COLS = [c for c in FEATURE_COLS if c.startswith("mhc_")]
STRUCT_ONLY = PEP_COLS  # peptide-derived
print(f"[shortcut] PEP_COLS({len(PEP_COLS)}), HLA_PCA_COLS({len(HLA_PCA_COLS)}), MHC_COLS({len(MHC_COLS)})", flush=True)

def col_idx(cols):
    return [FEATURE_COLS.index(c) for c in cols if c in FEATURE_COLS]

# Full feature set
p_full = train_lr_subset(col_idx(FEATURE_COLS), "Full")
# Drop peptide features
p_no_pep = train_lr_subset(col_idx([c for c in FEATURE_COLS if c not in PEP_COLS]), "DropPep")
# Drop HLA features
p_no_hla = train_lr_subset(col_idx([c for c in FEATURE_COLS if c not in HLA_PCA_COLS]), "DropHLA")
# Drop structure (peptide structural features), keep MHCflurry + HLA
p_no_struct = train_lr_subset(col_idx([c for c in FEATURE_COLS if c not in STRUCT_ONLY]), "DropStruct")
# Drop MHCflurry features
p_no_mhc = train_lr_subset(col_idx([c for c in FEATURE_COLS if c not in MHC_COLS]), "DropMHC")
# MHCflurry only
p_mhc_only = train_lr_subset(col_idx(MHC_COLS), "MHConly")
# Structure only
p_struct_only = train_lr_subset(col_idx(STRUCT_ONLY), "StructOnly")

def eval_all(p, name):
    rows = []
    for label, m in [("ITSNdb_no_overlap", no_overlap_mask),
                     ("ITSNdb_combined", combined_mask),
                     ("ITSNdb_in_master", in_master_mask)]:
        auc, lo, hi, n = auroc_with_ci(y[m], p[m])
        rows.append(dict(model=name, subset=label, n=n, AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi))
    return rows

for nm, p in [("Full", p_full), ("DropPep", p_no_pep), ("DropHLA", p_no_hla),
              ("DropStruct(KeepMHC+HLA)", p_no_struct), ("DropMHC(KeepStruct+HLA)", p_no_mhc),
              ("MHCflurryOnly", p_mhc_only), ("StructOnly", p_struct_only)]:
    shortcut_rows.extend(eval_all(p, nm))

# Per-HLA stratification (top-5 by ITSNdb count)
its_hlas = pd.Series(hla[ITS_MASK]).value_counts().head(5).index.tolist()
for h in its_hlas:
    m = ITS_MASK & (hla == h)
    if m.sum() < 8 or len(np.unique(y[m])) < 2: continue
    auc, lo, hi, n = auroc_with_ci(y[m], df_pred["p_w7b"].values[m])
    shortcut_rows.append(dict(model="W7-B_StackedLR", subset=f"PerHLA_{h}", n=n, AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi))

# Per-length stratification
for L in sorted(np.unique(df_feat["pep_len"])):
    m = ITS_MASK & (df_feat["pep_len"].values == L)
    if m.sum() < 8 or len(np.unique(y[m])) < 2: continue
    auc, lo, hi, n = auroc_with_ci(y[m], df_pred["p_w7b"].values[m])
    shortcut_rows.append(dict(model="W7-B_StackedLR", subset=f"PerLen_{int(L)}", n=n, AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi))

shortcut_df = pd.DataFrame(shortcut_rows)
shortcut_df.to_csv(WAVE7 / "wave7_shortcut_tests.tsv", sep="\t", index=False)
print(f"[shortcut] saved {len(shortcut_df)} rows", flush=True)
print(shortcut_df[shortcut_df["subset"].isin(["ITSNdb_no_overlap","ITSNdb_combined","ITSNdb_in_master"])].to_string(index=False))

# ---------------- Feature importance (LR coef + permutation) ----------------
print("\n[importance] LR coefficients + permutation importance...", flush=True)
sc_full = StandardScaler().fit(X[train_mask])
Xss_full = sc_full.transform(X)
lr_full = LogisticRegression(max_iter=4000, C=1.0, class_weight="balanced", solver="lbfgs", random_state=0)
lr_full.fit(Xss_full[train_mask], y[train_mask])
coefs = lr_full.coef_[0]

# Permutation importance on no_overlap for W7-B
print("[importance] permutation on ITSNdb_no_overlap (W7-B)...", flush=True)
y_test = y[no_overlap_mask]
X_test = Xss_full[no_overlap_mask]
base_p = lr_full.predict_proba(X_test)[:, 1]
base_auc = roc_auc_score(y_test, base_p) if len(np.unique(y_test)) > 1 else np.nan
rng = np.random.default_rng(0)
perm_imp = []
for j, col in enumerate(FEATURE_COLS):
    Xp = X_test.copy()
    drops = []
    for _ in range(10):
        Xp[:, j] = rng.permutation(X_test[:, j])
        p_ = lr_full.predict_proba(Xp)[:, 1]
        if len(np.unique(y_test)) > 1:
            drops.append(base_auc - roc_auc_score(y_test, p_))
    perm_imp.append(np.mean(drops) if drops else 0.0)

imp_df = pd.DataFrame({
    "feature": FEATURE_COLS,
    "lr_coef": coefs,
    "abs_lr_coef": np.abs(coefs),
    "perm_importance_no_overlap": perm_imp,
})
imp_df = imp_df.sort_values("perm_importance_no_overlap", ascending=False)
imp_df.to_csv(WAVE7 / "wave7_feature_importance.tsv", sep="\t", index=False)
print(f"[importance] top-10 by permutation:")
print(imp_df.head(10).to_string(index=False))

# ---------------- OOD detection ----------------
print("\n[ood] OOD detection: in_domain vs no_overlap using W7-A entropy...", flush=True)
# in_domain = ITSNdb_in_master, OOD = ITSNdb_no_overlap
ent = df_pred["p_w7a_entropy"].values
in_dom = in_master_mask
ood_mask = no_overlap_mask
labels_ood = np.concatenate([np.zeros(in_dom.sum()), np.ones(ood_mask.sum())])
scores_ood = np.concatenate([ent[in_dom], ent[ood_mask]])
ood_auc = roc_auc_score(labels_ood, scores_ood)
print(f"[ood] entropy AUROC for OOD (no_overlap vs in_master): {ood_auc:.3f}", flush=True)

# Compare to W7-B prediction-margin entropy
ent_b = -((df_pred["p_w7b"].values * np.log(np.clip(df_pred["p_w7b"].values, 1e-9, 1-1e-9))) +
          ((1-df_pred["p_w7b"].values) * np.log(np.clip(1-df_pred["p_w7b"].values, 1e-9, 1-1e-9))))
scores_ood_b = np.concatenate([ent_b[in_dom], ent_b[ood_mask]])
ood_auc_b = roc_auc_score(labels_ood, scores_ood_b)

ood_df = pd.DataFrame([
    dict(method="W7-A_entropy", AUROC_OOD=ood_auc),
    dict(method="W7-B_entropy", AUROC_OOD=ood_auc_b),
    dict(method="Wave1_predictive_entropy_reference", AUROC_OOD=0.636),
])
ood_df.to_csv(WAVE7 / "ood_results.tsv", sep="\t", index=False)
print(ood_df.to_string(index=False))

# ---------------- Counterfactual analysis (DRP §3.4) ----------------
print("\n[counterfactual] for each negative test peptide, find min L1 edit to flip prediction...", flush=True)
# For each negative test peptide, try every single-position substitution → see if probability flips above 0.5
AA20 = "ACDEFGHIKLMNPQRSTVWY"
its_neg = df_pred[(no_overlap_mask) & (df_pred["label"] == 0)].copy()
print(f"[counterfactual] {len(its_neg)} negatives in no_overlap to perturb", flush=True)

# Need to recompute features per substituted peptide, which requires reusing build_combined_features logic.
# Simpler: estimate via LR coefficients on AA-level features (peptide_len + KD + chg etc).
# For paper we'll approximate counterfactual by flipping single AA in feature space:
# perturb each peptide by simulating an "Ala substitution" at each position (set all peptide-pos features to neutral)
# and report which positions are most informative.

# Practical approach: take the LR linear discriminant in feature space; each peptide-row's
# feature vector is the sum of physico features + position-specific PWM. Flipping AA at position p
# changes pwm_p2/pwm_pOmega/pwm_mean/pwm_sum and physico features.
# We do an APPROXIMATE counterfactual: for each negative, compute Δ to flip = (0.5_logit - current_logit) /
# average per-position-AA-feature gradient norm. That gives a "mutational distance" estimate.

logit_thresh = 0.0  # decision boundary
current_logit = lr_full.decision_function(Xss_full[no_overlap_mask & (df_pred["label"]==0).values])
cf_dist = np.abs(logit_thresh - current_logit)

# Try simple position-level perturbation: zero out each peptide-feature column at a time
# and re-evaluate the prediction
PEP_FEAT_IDX = [FEATURE_COLS.index(c) for c in PEP_COLS]
neg_idx = np.where(no_overlap_mask & (df_pred["label"]==0).values)[0]
top_changes = {col: 0 for col in PEP_COLS}
for i in neg_idx:
    base_p = lr_full.predict_proba(Xss_full[i:i+1])[0, 1]
    for j_idx, col in zip(PEP_FEAT_IDX, PEP_COLS):
        Xtmp = Xss_full[i:i+1].copy()
        Xtmp[0, j_idx] = 0.0  # neutralize
        new_p = lr_full.predict_proba(Xtmp)[0, 1]
        if (base_p < 0.5) and (new_p >= 0.5):
            top_changes[col] += 1

cf_df = pd.DataFrame({
    "peptide_feature": list(top_changes.keys()),
    "n_negatives_flipped_when_neutralized": list(top_changes.values()),
    "n_total_negatives": len(neg_idx),
    "fraction_flipped": [v / max(1, len(neg_idx)) for v in top_changes.values()],
}).sort_values("n_negatives_flipped_when_neutralized", ascending=False)
cf_df.to_csv(WAVE7 / "wave7_counterfactual_edits.tsv", sep="\t", index=False)
print(cf_df.head(10).to_string(index=False))

# ---------------- Comparison table: W7 + Wave 1-6 baselines ----------------
print("\n[compare] building unified comparison table...", flush=True)
# W7 results from train_w7
w7_eval = pd.read_csv(RESULTS, sep="\t")

# Wave 1
W1_FILE = ROOT / "wave1" / "predictions_itsndb.tsv"
def add_baseline_aurocs(predictions_path, model_name, score_col, comp_rows):
    if not predictions_path.exists(): return
    p = pd.read_csv(predictions_path, sep="\t")
    if score_col not in p.columns:
        # Try alt names
        for alt in ["score", "prediction", "predictive_mean", "p_pos"]:
            if alt in p.columns: score_col = alt; break
    if "label" not in p.columns or score_col not in p.columns: return
    if "in_master" not in p.columns: p["in_master"] = False
    p["in_master"] = p["in_master"].fillna(False).astype(bool)
    # Combined
    if len(p) > 5 and len(p["label"].unique()) > 1:
        auc, lo, hi, n = auroc_with_ci(p["label"], p[score_col])
        ece = ece_15bin(p["label"].values, np.asarray(p[score_col].values, float))
        comp_rows.append(dict(model=model_name, subset="ITSNdb_combined", n=n, AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi, ECE=ece))
    nov = p[~p["in_master"]]
    if len(nov) > 5 and len(nov["label"].unique()) > 1:
        auc, lo, hi, n = auroc_with_ci(nov["label"], nov[score_col])
        ece = ece_15bin(nov["label"].values, np.asarray(nov[score_col].values, float))
        comp_rows.append(dict(model=model_name, subset="ITSNdb_no_overlap", n=n, AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi, ECE=ece))
    inm = p[p["in_master"]]
    if len(inm) > 5 and len(inm["label"].unique()) > 1:
        auc, lo, hi, n = auroc_with_ci(inm["label"], inm[score_col])
        ece = ece_15bin(inm["label"].values, np.asarray(inm[score_col].values, float))
        comp_rows.append(dict(model=model_name, subset="ITSNdb_in_master", n=n, AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi, ECE=ece))

baseline_rows = []
# Wave 1 Bayesian
add_baseline_aurocs(ROOT / "wave1" / "predictions_itsndb.tsv", "Wave1_BayesianMLP",
                    "predictive_mean", baseline_rows)
# Wave 3 MHCflurry
add_baseline_aurocs(ROOT / "wave3_mhcflurry" / "predictions.tsv", "Wave3_MHCflurry",
                    "score_mhcflurry", baseline_rows)
# Wave 3 BigMHC, PRIME, NetMHCpan, DeepImmuno, TransPHLA
for tool in ["bigmhc", "prime", "netmhcpan", "deepimmuno", "transphla"]:
    fldr = ROOT / f"wave3_{tool}"
    pred_files = list(fldr.glob("predictions*.tsv")) if fldr.exists() else []
    if not pred_files: continue
    p_path = pred_files[0]
    df_t = pd.read_csv(p_path, sep="\t")
    score_col = None
    for cand in ["score", f"score_{tool}", "prediction", "p_pos", "presentation_score", "predictive_mean", "y_pred", "ImmunoScore", "Immunogenicity"]:
        if cand in df_t.columns: score_col = cand; break
    if score_col is None:
        # try first numeric column not in label
        nums = df_t.select_dtypes("number").columns.tolist()
        nums = [c for c in nums if c not in ("label","in_master")]
        if nums: score_col = nums[0]
    if score_col is None: continue
    add_baseline_aurocs(p_path, f"Wave3_{tool.upper()}", score_col, baseline_rows)

# Add Wave2 from combined results
w2_path = ROOT / "wave2" / "wave2_combined_results.tsv"
if w2_path.exists():
    w2 = pd.read_csv(w2_path, sep="\t")
    for _, r in w2.iterrows():
        sub = r["subset"]
        if sub not in ("ITSNdb_no_overlap","ITSNdb_combined","ITSNdb_in_master"): continue
        if pd.notna(r["AUROC"]) and pd.notna(r.get("n_total")) and r.get("n_total", 0) > 5:
            baseline_rows.append(dict(model=f"Wave2_{r['model']}", subset=sub, n=int(r["n_total"]),
                                      AUROC=float(r["AUROC"]), AUROC_lo95=np.nan, AUROC_hi95=np.nan, ECE=np.nan))

# Wave 4 / 4C: Look for results jsons
for w4 in ["wave4a", "wave4b", "wave4c"]:
    fldr = ROOT / w4
    if not fldr.exists(): continue
    for js in fldr.glob("*.json"):
        try:
            d = json.loads(js.read_text())
        except Exception:
            continue
        # Heuristic: recurse for AUROC keys
        def walk(o, path=""):
            if isinstance(o, dict):
                if "AUROC" in o and "subset" in o:
                    yield (o.get("model", path), o["subset"], o["AUROC"])
                for k, v in o.items():
                    yield from walk(v, f"{path}.{k}" if path else k)
            elif isinstance(o, list):
                for i, v in enumerate(o):
                    yield from walk(v, f"{path}[{i}]")
        for tup in walk(d, w4):
            pass  # not standardized; skip
    # Also check wave4_results.tsv
    rfile = fldr / f"{w4}_results.tsv"
    if not rfile.exists():
        rfile_cands = list(fldr.glob("*results*.tsv"))
        if rfile_cands: rfile = rfile_cands[0]
    if rfile.exists():
        try:
            df_r = pd.read_csv(rfile, sep="\t")
            for _, r in df_r.iterrows():
                if "AUROC" in r and "subset" in r and pd.notna(r["AUROC"]):
                    baseline_rows.append(dict(model=f"{w4.upper()}_{r.get('model','?')}",
                                              subset=str(r["subset"]),
                                              n=int(r.get("n", 0)) if pd.notna(r.get("n")) else 0,
                                              AUROC=float(r["AUROC"]),
                                              AUROC_lo95=float(r.get("AUROC_lo95", np.nan)) if pd.notna(r.get("AUROC_lo95", np.nan)) else np.nan,
                                              AUROC_hi95=float(r.get("AUROC_hi95", np.nan)) if pd.notna(r.get("AUROC_hi95", np.nan)) else np.nan,
                                              ECE=float(r.get("ECE", np.nan)) if pd.notna(r.get("ECE", np.nan)) else np.nan))
        except Exception as e:
            print(f"[compare] skip {rfile}: {e}", flush=True)

# Combine W7 with baselines
w7_eval["ECE"] = w7_eval["ECE"].astype(float) if "ECE" in w7_eval.columns else np.nan
all_rows = w7_eval.to_dict("records") + baseline_rows
all_df = pd.DataFrame(all_rows)
# Reorder columns
cols = ["model","subset","n","AUROC","AUROC_lo95","AUROC_hi95","ECE"]
for c in cols:
    if c not in all_df.columns: all_df[c] = np.nan
all_df = all_df[cols].drop_duplicates(subset=["model","subset"]).reset_index(drop=True)
all_df.to_csv(WAVE7 / "wave7_results.tsv", sep="\t", index=False)
print(f"[compare] saved {len(all_df)} rows to wave7_results.tsv", flush=True)
print(all_df[all_df["subset"]=="ITSNdb_no_overlap"].sort_values("AUROC", ascending=False).to_string(index=False))

# ---------------- Hero figures ----------------
print("\n[fig] generating hero figures...", flush=True)

# Fig 1: hero forest plot of ITSNdb_no_overlap AUROCs across all methods
no_ov = all_df[all_df["subset"] == "ITSNdb_no_overlap"].copy()
no_ov = no_ov.sort_values("AUROC", ascending=True).reset_index(drop=True)
plt.figure(figsize=(10, max(5, len(no_ov)*0.3)))
y_pos = np.arange(len(no_ov))
aucs = no_ov["AUROC"].values
los = np.where(no_ov["AUROC_lo95"].isna(), aucs, no_ov["AUROC_lo95"].values)
his = np.where(no_ov["AUROC_hi95"].isna(), aucs, no_ov["AUROC_hi95"].values)
errs_lo = aucs - los
errs_hi = his - aucs
colors = ["#d62728" if "W7-" in m else ("#1f77b4" if "Wave1" in m else "#7f7f7f") for m in no_ov["model"]]
plt.errorbar(aucs, y_pos, xerr=[errs_lo, errs_hi], fmt='o', ecolor='gray', capsize=3, mec='black')
for i, c in enumerate(colors):
    plt.plot(aucs[i], y_pos[i], 'o', color=c, markersize=8)
plt.axvline(0.5, color="gray", linestyle="--", alpha=0.5, label="chance")
plt.axvline(0.668, color="orange", linestyle="--", alpha=0.5, label="MHCflurry alone")
plt.yticks(y_pos, no_ov["model"], fontsize=8)
plt.xlabel("AUROC (ITSNdb no_overlap; bars = 95% bootstrap CI)")
plt.title("Wave 7 — Hero forest: ITSNdb no_overlap AUROC across methods")
plt.legend(loc="lower right", fontsize=8)
plt.tight_layout()
plt.savefig(WAVE7 / "fig_wave7_hero_forest.png", dpi=150)
plt.savefig(WAVE7 / "fig_wave7_hero_forest.pdf")
plt.close()

# Fig 2: ECE comparison across methods (no_overlap subset where ECE is available)
ece_df = all_df[(all_df["subset"]=="ITSNdb_no_overlap") & all_df["ECE"].notna()].sort_values("ECE")
plt.figure(figsize=(8, max(3, len(ece_df)*0.3)))
plt.barh(np.arange(len(ece_df)), ece_df["ECE"], color=["#d62728" if "W7-" in m else "#7f7f7f" for m in ece_df["model"]])
plt.yticks(np.arange(len(ece_df)), ece_df["model"], fontsize=8)
plt.xlabel("ECE (15-bin)")
plt.title("Wave 7 — ECE on ITSNdb no_overlap (lower = better calibration)")
plt.axvline(0.05, color="green", linestyle="--", alpha=0.5, label="ECE=0.05")
plt.legend(loc="lower right", fontsize=8)
plt.tight_layout()
plt.savefig(WAVE7 / "fig_wave7_calibration_panel.png", dpi=150)
plt.savefig(WAVE7 / "fig_wave7_calibration_panel.pdf")
plt.close()

# Fig 3: Inflation gap (in_master AUROC − no_overlap AUROC)
piv = all_df.pivot_table(index="model", columns="subset", values="AUROC", aggfunc="first")
if "ITSNdb_in_master" in piv.columns and "ITSNdb_no_overlap" in piv.columns:
    piv["inflation_gap"] = piv["ITSNdb_in_master"] - piv["ITSNdb_no_overlap"]
    piv2 = piv.dropna(subset=["inflation_gap"]).sort_values("inflation_gap", ascending=True)
    plt.figure(figsize=(8, max(3, len(piv2)*0.3)))
    plt.barh(np.arange(len(piv2)), piv2["inflation_gap"], color=["#d62728" if "W7-" in m else "#7f7f7f" for m in piv2.index])
    plt.yticks(np.arange(len(piv2)), piv2.index, fontsize=8)
    plt.xlabel("Δ AUROC (in_master − no_overlap)")
    plt.title("Wave 7 — Inflation gap: how much each method's AUROC drops on truly held-out peptides")
    plt.axvline(0, color="black", alpha=0.5)
    plt.tight_layout()
    plt.savefig(WAVE7 / "fig_wave7_inflation_gap.png", dpi=150)
    plt.savefig(WAVE7 / "fig_wave7_inflation_gap.pdf")
    plt.close()

# Fig 4: Feature importance bar plot
top_imp = imp_df.head(20).iloc[::-1]
plt.figure(figsize=(8, 7))
plt.barh(np.arange(len(top_imp)), top_imp["perm_importance_no_overlap"], color="#1f77b4")
plt.yticks(np.arange(len(top_imp)), top_imp["feature"], fontsize=9)
plt.xlabel("Permutation importance Δ AUROC (no_overlap)")
plt.title("Wave 7 — Top 20 features by permutation importance on no_overlap")
plt.tight_layout()
plt.savefig(WAVE7 / "fig_wave7_feature_importance.png", dpi=150)
plt.savefig(WAVE7 / "fig_wave7_feature_importance.pdf")
plt.close()

# Fig 5: Per-allele forest plot (W7-A vs W7-B)
pa = w7_eval[w7_eval["subset"].str.startswith("PerAllele_")].copy()
if len(pa) > 0:
    pa["allele"] = pa["subset"].str.replace("PerAllele_", "")
    fig, ax = plt.subplots(figsize=(9, max(5, pa["allele"].nunique()*0.4)))
    alleles = sorted(pa["allele"].unique())
    for i, h in enumerate(alleles):
        sub = pa[pa["allele"] == h]
        for _, r in sub.iterrows():
            color = "#d62728" if "W7-A" in r["model"] else "#1f77b4"
            offset = -0.15 if "W7-A" in r["model"] else 0.15
            ax.errorbar(r["AUROC"], i + offset,
                        xerr=[[r["AUROC"] - r["AUROC_lo95"]], [r["AUROC_hi95"] - r["AUROC"]]],
                        fmt='o', color=color, capsize=3, mec='black',
                        label=r["model"] if i == 0 else "")
    ax.axvline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_yticks(np.arange(len(alleles)))
    ax.set_yticklabels(alleles, fontsize=8)
    ax.set_xlabel("AUROC (per-allele)")
    ax.set_title("Wave 7 — Per-allele forest plot (W7-A red vs W7-B blue)")
    ax.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(WAVE7 / "fig_wave7_per_allele_forest.png", dpi=150)
    plt.savefig(WAVE7 / "fig_wave7_per_allele_forest.pdf")
    plt.close()

print(f"[done] elapsed {time.time()-t0:.1f}s", flush=True)
