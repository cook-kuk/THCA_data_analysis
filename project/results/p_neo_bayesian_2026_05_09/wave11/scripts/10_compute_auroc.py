#!/usr/bin/env python3
"""Wave 11 — synthesize all per-(algo, bundle) predictions into mega table.

For each predictions TSV at wave11_predictions/<algo>__<bundle>.tsv:
- Compute AUROC + 1000-bootstrap 95% CI on all rows
- Stratify by in_training vs external (using overlap audit table)
- Compute ECE (10-bin)
- Stack into wave11_mega_table_expanded.tsv

Skip AUROC if n_pos < 5 (low-power flag).
For positive-only or negative-only test sets, report n_pos / n only.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE11 = ROOT / "wave11"
PRED = WAVE11 / "wave11_predictions"
BUN = WAVE11 / "test_bundles"

# Load overlap audit
audit = pd.read_csv(WAVE11 / "wave11_overlap_audit.tsv", sep="\t")

# Map algorithm name in filename -> algorithm key in audit
ALGO_FILE_MAP = {
    "MHCflurry": "MHCflurry",
    "BigMHC_IM": "BigMHC_IM",
    "PRIME": "PRIME",
    "DeepImmuno": "DeepImmuno",
    "TransPHLA": "TransPHLA",
    "NetMHCpan_4.1": "NetMHCpan_4.1",
    "NetMHCpan": "NetMHCpan_4.1",
    "RF_biophys": "RF_biophys",
    "LR_biophys": "RF_biophys",  # share overlap rules
    "Structure_LR": "Structure_LR",
    "GP_quantum": "GP_quantum",
    "ESM2_Bayesian": "ESM2_Bayesian",
    "MHCflurry_features_LR": "MHCflurry_features_LR",
    "Wave7_combined": "Wave7_combined",
    "Wave8_self_sim": "Wave8_self_sim",
}

# Score column priority (most-canonical first)
SCORE_COL = {
    "MHCflurry":  "score_mhcflurry",
    "BigMHC_IM":  "score_bigmhc_im",
    "PRIME":      "score_prime",
    "DeepImmuno": "score_deepimmuno",
    "TransPHLA":  "score_transphla",
    "NetMHCpan_4.1": "score_netmhcpan",
    "NetMHCpan":  "score_netmhcpan",
    "RF_biophys": "score_rf_biophys",
    "LR_biophys": "score_lr_biophys",
    "Structure_LR":"score_structure_lr",
    "GP_quantum": "score_gp_quantum",
    "ESM2_Bayesian":"score_esm2_bayes",
    "MHCflurry_features_LR":"score_mhcflurry_features_lr",
    "Wave7_combined":"score_wave7",
    "Wave8_self_sim":"score_wave8",
}

def auroc_with_ci(y, s, n_boot=1000, seed=0):
    y = np.asarray(y, dtype=int); s = np.asarray(s, dtype=float)
    if len(y) == 0 or len(np.unique(y)) < 2:
        return (np.nan, np.nan, np.nan)
    auc = roc_auc_score(y, s)
    rng = np.random.default_rng(seed)
    n = len(y); aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yy = y[idx]; ss = s[idx]
        if len(np.unique(yy)) < 2:
            continue
        aucs.append(roc_auc_score(yy, ss))
    if not aucs:
        return (auc, np.nan, np.nan)
    return (float(auc), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5)))

def ece_score(y, s, n_bins=10):
    y = np.asarray(y, dtype=int); s = np.asarray(s, dtype=float)
    if len(y) == 0:
        return np.nan
    s_min, s_max = s.min(), s.max()
    if s_max - s_min < 1e-9:
        return np.nan
    s_norm = (s - s_min) / (s_max - s_min)
    bins = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for i in range(n_bins):
        mask = (s_norm >= bins[i]) & (s_norm < bins[i+1] if i < n_bins-1 else s_norm <= bins[i+1])
        if mask.sum() == 0: continue
        conf = s_norm[mask].mean()
        acc = y[mask].mean()
        ece += (mask.sum() / len(y)) * abs(conf - acc)
    return float(ece)

def build_in_training_mask(df_pred, algo, bname):
    """Recompute in_training mask using same logic as overlap_audit."""
    sub = audit[(audit["algorithm"] == algo) & (audit["test_bundle"] == bname)]
    if len(sub) == 0:
        return None
    mode = sub.iloc[0]["train_proxy"]
    # We need the master & iedb pool
    master_pool_pep, master_pool_pair, iedb_pep, iedb_pair = _MASTER_CACHE
    pair = df_pred["peptide"] + "|" + df_pred["hla"].fillna("")
    in_master = pair.isin(master_pool_pair) | df_pred["peptide"].isin(master_pool_pep)
    in_iedb = pair.isin(iedb_pair) | df_pred["peptide"].isin(iedb_pep)
    if mode == "master":
        return in_master
    if mode == "iedb":
        return in_iedb
    if mode == "master+iedb":
        return in_master | in_iedb
    return pd.Series([False] * len(df_pred))

print("[load] master + IEDB pools...", flush=True)
bun_full = pd.read_csv(ROOT / "bundle.tsv", sep="\t")
mp = bun_full[bun_full["split"] == "train"]
master_pool_pep = set(mp["peptide"].tolist())
master_pool_pair = set((mp["peptide"] + "|" + mp["HLA_norm"].fillna("")).tolist())
master_full = pd.read_csv("/data/neoantigen_vaccine_hub/experiments/what_matters_neoantigen/cache/benchmark_clean.tsv",
                          sep="\t", low_memory=False)
iedb = master_full[master_full["source"] == "IEDB_tcell_v3"]
iedb_pep = set(iedb["peptide"].astype(str).tolist())
iedb_pair = set((iedb["peptide"].astype(str) + "|" + iedb["HLA"].astype(str)).tolist())
_MASTER_CACHE = (master_pool_pep, master_pool_pair, iedb_pep, iedb_pair)
print(f"  master: {len(master_pool_pep)} peps; iedb: {len(iedb_pep)} peps", flush=True)

# Discover prediction files
pred_files = sorted(PRED.glob("*.tsv"))
print(f"[scan] {len(pred_files)} prediction files", flush=True)

mega = []
for pf in pred_files:
    fn = pf.stem  # e.g. "MHCflurry__cedar_partial"
    parts = fn.split("__")
    if len(parts) != 2: continue
    algo_file, bname = parts
    algo = ALGO_FILE_MAP.get(algo_file, algo_file)
    score_col = SCORE_COL.get(algo, None)
    df = pd.read_csv(pf, sep="\t")
    if score_col is None or score_col not in df.columns:
        # try the most numeric-looking score col
        cands = [c for c in df.columns if c.startswith("score_")]
        if not cands: continue
        score_col = cands[0]

    df = df.dropna(subset=[score_col])
    if len(df) == 0: continue

    # Compute in_training mask
    in_tr = build_in_training_mask(df, algo, bname)
    if in_tr is None:
        in_tr = pd.Series([False] * len(df))

    n = len(df)
    n_pos = int(df["label"].sum())
    y_all = df["label"].astype(int).to_numpy()
    s_all = df[score_col].astype(float).to_numpy()

    # FULL — combined
    if n_pos >= 5 and (n - n_pos) >= 5:
        auc, lo, hi = auroc_with_ci(y_all, s_all)
        ece_all = ece_score(y_all, s_all)
    else:
        auc = lo = hi = ece_all = np.nan

    # in_training stratum
    sub_in = df[in_tr]
    n_in = len(sub_in); npos_in = int(sub_in["label"].sum()) if n_in > 0 else 0
    if n_in >= 10 and npos_in >= 5 and (n_in - npos_in) >= 5:
        auc_in, _, _ = auroc_with_ci(sub_in["label"].astype(int), sub_in[score_col])
    else:
        auc_in = np.nan

    # external stratum
    sub_ext = df[~in_tr]
    n_ext = len(sub_ext); npos_ext = int(sub_ext["label"].sum()) if n_ext > 0 else 0
    if n_ext >= 10 and npos_ext >= 5 and (n_ext - npos_ext) >= 5:
        auc_ext, lo_ext, hi_ext = auroc_with_ci(sub_ext["label"].astype(int), sub_ext[score_col])
    else:
        auc_ext = lo_ext = hi_ext = np.nan

    inflation = auc_in - auc_ext if (not np.isnan(auc_in) and not np.isnan(auc_ext)) else np.nan
    low_power = (n_pos < 5) or (n - n_pos < 5)

    mega.append(dict(
        algorithm=algo, test_bundle=bname,
        n=n, n_pos=n_pos,
        AUROC=round(auc, 4) if not np.isnan(auc) else np.nan,
        CI_lo95=round(lo, 4) if not np.isnan(lo) else np.nan,
        CI_hi95=round(hi, 4) if not np.isnan(hi) else np.nan,
        n_in_training=n_in,
        in_training_AUROC=round(auc_in, 4) if not np.isnan(auc_in) else np.nan,
        n_external=n_ext,
        external_AUROC=round(auc_ext, 4) if not np.isnan(auc_ext) else np.nan,
        external_CI_lo95=round(lo_ext, 4) if not np.isnan(lo_ext) else np.nan,
        external_CI_hi95=round(hi_ext, 4) if not np.isnan(hi_ext) else np.nan,
        inflation_delta=round(inflation, 4) if not np.isnan(inflation) else np.nan,
        ECE=round(ece_all, 4) if not np.isnan(ece_all) else np.nan,
        low_power=int(low_power),
    ))

mega = pd.DataFrame(mega).sort_values(["algorithm", "test_bundle"]).reset_index(drop=True)
mega.to_csv(WAVE11 / "wave11_mega_table_expanded.tsv", sep="\t", index=False)
print(f"\n[done] mega table: {len(mega)} (algo, bundle) cells")
print(mega.to_string(index=False))
