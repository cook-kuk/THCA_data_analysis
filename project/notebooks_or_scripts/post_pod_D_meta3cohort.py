#!/usr/bin/env python3
"""Post-Pod D — Process K2 STAR re-quantification result + 3-cohort meta update.

Trigger: after Pod D's run.log shows "ALL DONE"
Input: /workspace/k2_star/results/k2_counts.tsv (scp'd back to main VM)
Output: project/results/d7p3_k2_starred/{k2_8gene_counts.tsv, 3cohort_meta_updated.json}
"""
from __future__ import annotations
from pathlib import Path
import json, subprocess
import numpy as np
import pandas as pd
from scipy import stats

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/d7p3_k2_starred"
RES.mkdir(parents=True, exist_ok=True)

# Step 1 — Read STAR counts (8-gene focused)
counts_path = RES / "k2_counts.tsv"
if not counts_path.exists():
    print(f"❌ {counts_path} not yet retrieved — run scp from Pod D first")
    print("   ssh -p <port> root@<ip> -- cat /workspace/k2_star/results/k2_counts.tsv > {counts_path}")
    raise SystemExit(1)

raw = pd.read_csv(counts_path, sep="\t", comment="#")
print(f"  raw counts: {raw.shape}")
print(f"  columns: {raw.columns.tolist()[:10]}")

# Step 2 — Convert ENSG to gene_symbol
ensg_map = pd.read_csv(PROJ / "metadata/ensembl_to_symbol.tsv", sep="\t")
ensg_to_sym = dict(zip(ensg_map["ensembl_id"], ensg_map["gene_symbol"]))

# featureCounts output: Geneid (ENSG), Chr, Start, End, Strand, Length, then sample counts
GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']
sym_to_ensg = {v: k for k, v in ensg_to_sym.items()}
target_ensg = [sym_to_ensg.get(g) for g in GENE_8 if g in sym_to_ensg]

raw["ensg_clean"] = raw["Geneid"].astype(str).str.split(".").str[0]
g8_subset = raw[raw["ensg_clean"].isin(target_ensg)].copy()
g8_subset["gene_symbol"] = g8_subset["ensg_clean"].map(ensg_to_sym)
print(f"  8-gene subset: {len(g8_subset)} rows")
print(g8_subset[["gene_symbol", "ensg_clean"]].to_string(index=False))

# Step 3 — Build 8-gene count matrix
sample_cols = [c for c in raw.columns if c not in ["Geneid", "Chr", "Start", "End", "Strand", "Length", "ensg_clean", "gene_symbol"]]
g8_mat = g8_subset.set_index("gene_symbol")[sample_cols]
g8_mat = g8_mat.loc[[g for g in GENE_8 if g in g8_mat.index]]
print(f"\n  8-gene count matrix: {g8_mat.shape}")

# Step 4 — Library-size normalize + log + within-sample center
total_counts = pd.read_csv(counts_path, sep="\t", comment="#", usecols=lambda c: c not in ["Chr", "Start", "End", "Strand", "Length"])
lib_size = total_counts.iloc[:, 6:].sum(axis=0)  # column sums
lib_size_per_sample = lib_size.to_dict()

g8_log_cpm = g8_mat.copy()
for s in g8_mat.columns:
    libN = lib_size_per_sample.get(s, 1)
    g8_log_cpm[s] = np.log2(g8_mat[s] / libN * 1e6 + 1)
print(f"  log2(CPM+1) range: [{g8_log_cpm.values.min():.2f}, {g8_log_cpm.values.max():.2f}]")

# Step 5 — Within-sample-centered profile (subtract sample-level mean of all genes)
# For full-cohort centering need full matrix:
all_log_cpm = total_counts.iloc[:, 6:].copy()
all_log_cpm = all_log_cpm.div(lib_size, axis=1) * 1e6
all_log_cpm = np.log2(all_log_cpm + 1)
sample_mean = all_log_cpm.mean(axis=0)  # per-sample mean over all genes
g8_centered = g8_log_cpm.sub(sample_mean, axis=1)

# Save
g8_log_cpm.to_csv(RES / "k2_starred_8gene_logCPM.tsv", sep="\t")
g8_centered.to_csv(RES / "k2_starred_8gene_centered.tsv", sep="\t")
print(f"  centered profile saved")

# Step 6 — Apply TCGA-trained classifier (centered profile method)
import pickle
classifier_path = PROJ / "results/v17_korean/K2_dm_predictor.pkl"
if classifier_path.exists():
    print(f"\n  Loading TCGA-trained classifier from {classifier_path}")
    with open(classifier_path, "rb") as f:
        clf = pickle.load(f)
    feat = g8_centered.T  # samples × 8 genes
    proba = clf.predict_proba(feat)
    p_dm1 = proba[:, 1] if proba.shape[1] >= 2 else proba[:, 0]
    print(f"  K2 STAR-based P_DM1 distribution: mean={p_dm1.mean():.3f}, sd={p_dm1.std():.3f}")
else:
    print(f"\n  ! Classifier not found, building fresh from TCGA centered profile")
    # Fall-back: train fresh on TCGA
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline

    tcga_expr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                              sep="\t", index_col=0)
    cluster = pd.read_csv(PROJ / "results/v17_realfix/R1A_cluster_labels.tsv", sep="\t")
    cluster["DM"] = cluster["cluster"].str[:3]
    samps_tcga = [s for s in cluster["sample_id"] if s in tcga_expr.columns]
    cluster_t = cluster[cluster["sample_id"].isin(samps_tcga)].set_index("sample_id").loc[samps_tcga]
    tcga_g8 = [g for g in GENE_8 if g in tcga_expr.index]
    tcga_means = tcga_expr.mean(axis=0)
    tcga_centered = tcga_expr.subtract(tcga_means, axis=1)
    Xtr = tcga_centered.loc[tcga_g8, samps_tcga].T.values
    ytr = (cluster_t["DM"] == "DM1").astype(int).values
    pipe = Pipeline([("sc", StandardScaler()), ("lr", LogisticRegression(C=1.0, max_iter=1000))])
    pipe.fit(Xtr, ytr)
    feat = g8_centered.T[tcga_g8].values
    proba = pipe.predict_proba(feat)
    p_dm1 = proba[:, 1]
    print(f"  Fresh classifier — K2 STAR-based P_DM1: mean={p_dm1.mean():.3f}, sd={p_dm1.std():.3f}")

# Step 7 — 3-cohort meta with K2 raw-score d
print("\n=== 3-cohort meta-analysis (TCGA + GSE286332 + K2 STAR-based) ===")
# K2 effect size: high-pDM1 vs low-pDM1 split → mean 8-gene RAI score difference
k2_score = g8_centered.mean(axis=0)  # per-sample 8-gene mean centered
high_p = p_dm1 >= 0.5
low_p = p_dm1 < 0.5
if high_p.sum() > 5 and low_p.sum() > 5:
    a = k2_score[high_p].values
    b = k2_score[low_p].values
    sp = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
    d_k2 = (a.mean() - b.mean()) / max(sp, 1e-9)
    n_high = int(high_p.sum())
    n_low = int(low_p.sum())
    print(f"  K2 STAR DM1-pred (n={n_high}) vs DM2-pred (n={n_low}): Cohen d = {d_k2:+.3f}")
else:
    print(f"  K2 high/low split too unbalanced (high={high_p.sum()}, low={low_p.sum()})")
    d_k2 = None
    n_high = int(high_p.sum())
    n_low = int(low_p.sum())

if d_k2 is not None:
    rows = [
        {"cohort": "TCGA-THCA (DM2 vs DM1, n=500)", "n_case": 360, "n_ctrl": 140, "cohen_d": -1.783},
        {"cohort": "GSE286332 (PTC+HT vs PTC, n=18)", "n_case": 9, "n_ctrl": 9, "cohen_d": -1.602},
        {"cohort": f"K2 STAR (DM1-pred vs DM2-pred, n={n_high}+{n_low})", "n_case": n_low, "n_ctrl": n_high, "cohen_d": -d_k2},
    ]

    def pool_d(rows):
        ds = np.array([r["cohen_d"] for r in rows])
        n1 = np.array([r["n_case"] for r in rows])
        n2 = np.array([r["n_ctrl"] for r in rows])
        v_d = (n1+n2)/(n1*n2) + (ds**2)/(2*(n1+n2))
        w = 1.0 / v_d
        d_fix = (w * ds).sum() / w.sum()
        Q = (w * (ds - d_fix)**2).sum()
        df_h = len(ds) - 1
        tau2 = max(0.0, (Q - df_h) / (w.sum() - (w**2).sum() / w.sum())) if df_h > 0 else 0.0
        w_re = 1.0 / (v_d + tau2)
        d_re = (w_re * ds).sum() / w_re.sum()
        var_re = 1.0 / w_re.sum()
        return dict(d_pooled=round(float(d_re), 3),
                    ci_lo=round(float(d_re - 1.96 * np.sqrt(var_re)), 3),
                    ci_hi=round(float(d_re + 1.96 * np.sqrt(var_re)), 3),
                    Q=round(float(Q), 3), df=df_h, I2_pct=round(max(0, (Q-df_h)/Q*100), 1) if Q > 0 else 0,
                    tau2=round(float(tau2), 4))

    pooled = pool_d(rows)
    print(f"\n  3-cohort pooled d: {pooled['d_pooled']:+.3f} [95% CI {pooled['ci_lo']:+.3f}, {pooled['ci_hi']:+.3f}]")
    print(f"  Cochran's Q: {pooled['Q']}, I²={pooled['I2_pct']}%, τ²={pooled['tau2']}")

    summary = {
        "k2_starred_n": int(g8_log_cpm.shape[1]),
        "k2_8gene_centered_mean": float(k2_score.mean()),
        "k2_p_dm1_distribution": {"mean": float(p_dm1.mean()), "sd": float(p_dm1.std())},
        "k2_d_dm1pred_vs_dm2pred": d_k2,
        "meta_3cohort_rows": rows,
        "pooled": pooled,
    }
    (RES / "3cohort_meta_updated.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\n✓ {RES / '3cohort_meta_updated.json'}")

print("\n✓ Post-Pod D processing complete")
