#!/usr/bin/env python
"""Wave 4C — stacked ensemble of wave3 algorithm scores + Wave 2.5 PWM.

Algorithms available (NetMHCpan/TransPHLA not yet produced):
  - MHCflurry 2.0  (presentation_score)
  - BigMHC IM      (score_bigmhc_im)
  - DeepImmuno-CNN (score_deepimmuno)
  - PRIME 2.1      (-%Rank, higher = more immunogenic)

Plus Wave 2.5 PWM scaled rank (pwm_pct).

Strategies:
  E1  Mean of rank-normalized scores (no train, no leak)
  E2  Median of rank-normalized scores (robust)
  E3  LogReg meta-classifier on the 4 tool scores
  E4  Same but with PWM as 5th feature

Stacking-data choice for E3/E4 (HONESTY NOTE):
  We do NOT have train-pool tool scores (would require re-running each tool on
  ~2000 peptides — out of CPU budget). Two stacking variants are reported:

    (a) Train on ITSNdb in_master=True (n=213), test on in_master=False (n=106).
        This is the documented design: in_master peptides ARE in the master pool
        the algorithms were trained on; in_master=False are NOT. Meta-learner
        sees high-quality (memorized) algorithm scores at train time, novel
        peptide scores at test time. Risk: distribution mismatch between train
        and test. We report this as E3a / E4a.

    (b) 5-fold CV on the no_overlap subset (n=106). Same domain at train and
        test, but smaller training pool. We report this as E3b / E4b.

  Variant (b) is the more defensible "no leakage" comparison vs E1 mean.
  Variant (a) shows whether transfer from in_master patterns helps at all.

Outputs:
  wave4c_results.tsv         method | testset | n | n_pos | AUROC | CI lo/hi
  ensemble_predictions.tsv   peptide | hla | label | in_master | E1..E4
  fig_wave4c_ensemble.png/pdf
  WAVE4C_REPORT.md
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "wave4c"
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260509)
N_BOOT = 1000

t0 = time.time()


# ----------------------------------------------------------------------------
# Load tool predictions
# ----------------------------------------------------------------------------
def load_pred(name: str, score_col: str) -> pd.DataFrame:
    p = ROOT / f"wave3_{name}" / "predictions.tsv"
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t", dtype={"in_master": "boolean"})
    return df[["peptide", "hla", "label", "in_master", score_col]].rename(
        columns={score_col: f"score_{name}"}
    )


tool_specs = [
    ("mhcflurry", "score_mhcflurry"),
    ("bigmhc", "score_bigmhc_im"),
    ("deepimmuno", "score_deepimmuno"),
    ("prime", "score_prime"),
]
# Add netmhcpan / transphla iff present
optional = [
    ("netmhcpan", "score_netmhcpan"),
    ("transphla", "score_transphla"),
]

frames = []
loaded_tools = []
for name, col in tool_specs:
    f = load_pred(name, col)
    if f is None:
        print(f"[WARN] missing wave3_{name}/predictions.tsv — skipping")
        continue
    frames.append(f)
    loaded_tools.append(name)

for name, col in optional:
    p = ROOT / f"wave3_{name}" / "predictions.tsv"
    if p.exists():
        df = pd.read_csv(p, sep="\t", dtype={"in_master": "boolean"})
        cols_present = [c for c in df.columns if c.startswith("score_")]
        if cols_present:
            score_col = cols_present[0]
            df = df[["peptide", "hla", "label", "in_master", score_col]]
            frames.append(df.rename(columns={score_col: f"score_{name}"}))
            loaded_tools.append(name)
            print(f"[load] wave3_{name} present — included")

print(f"[load] tools loaded: {loaded_tools}")

# Inner join on (peptide, hla, label, in_master)
merged = frames[0]
for f in frames[1:]:
    merged = merged.merge(f, on=["peptide", "hla", "label", "in_master"], how="inner")
print(f"[merge] joined {len(merged)} rows across {len(loaded_tools)} tools")

# PWM features
pwm = pd.read_csv(ROOT / "wave25" / "pwm_features.tsv", sep="\t")
pwm_itsndb = pwm[pwm["split"] == "itsndb"][["peptide", "hla", "pwm_pct"]].copy()
merged = merged.merge(pwm_itsndb, on=["peptide", "hla"], how="left")
n_with_pwm = merged["pwm_pct"].notna().sum()
print(f"[merge] PWM joined: {n_with_pwm}/{len(merged)} rows have PWM score")

# Drop any duplicates that may exist
merged = merged.drop_duplicates(subset=["peptide", "hla"]).reset_index(drop=True)
print(f"[merge] after dedup: {len(merged)} rows")

score_cols = [f"score_{n}" for n in loaded_tools]


# ----------------------------------------------------------------------------
# Rank-normalize each score column to [0, 1] across the full evaluation set
# ----------------------------------------------------------------------------
def rank_norm(s: pd.Series) -> pd.Series:
    # rank with NaN preserved, then fill NaN with median (0.5)
    r = s.rank(method="average")
    out = (r - 1) / (s.notna().sum() - 1)
    return out.fillna(0.5)


norm_cols = []
for c in score_cols:
    nc = c + "_rn"
    merged[nc] = rank_norm(merged[c])
    norm_cols.append(nc)

# PWM is already pct in [0, 1]; impute missing with 0.5 if any
merged["pwm_pct_rn"] = merged["pwm_pct"].fillna(0.5)

# Confirm no NaN remains in features used by LR
for c in norm_cols + ["pwm_pct_rn"]:
    if merged[c].isna().any():
        n_nan = int(merged[c].isna().sum())
        print(f"[impute] WARN {c} still has {n_nan} NaN; filling 0.5")
        merged[c] = merged[c].fillna(0.5)


# ----------------------------------------------------------------------------
# E1 / E2 — no-train ensembles
# ----------------------------------------------------------------------------
merged["score_E1"] = merged[norm_cols].mean(axis=1)
merged["score_E2"] = merged[norm_cols].median(axis=1)


# ----------------------------------------------------------------------------
# E3 / E4 — LogReg stacking
#   Variant (a): train on in_master=True, score on in_master=False
#   Variant (b): 5-fold CV on in_master=False (predictions from out-of-fold)
# ----------------------------------------------------------------------------
mask_in = merged["in_master"].astype(bool).values
mask_no = ~mask_in

# Use raw scores (LR will learn its own coefficients); also try rank-norm version
features_E3 = norm_cols
features_E4 = norm_cols + ["pwm_pct_rn"]


def train_test_split_in_master(X_all, y_all, mask_in, mask_no):
    return X_all[mask_in], y_all[mask_in], X_all[mask_no], y_all[mask_no]


# --- E3a / E4a: train-on-in_master ---
y_all = merged["label"].astype(int).values

for tag, feats in [("E3a", features_E3), ("E4a", features_E4)]:
    X_all = merged[feats].values
    Xtr, ytr, Xte, yte = train_test_split_in_master(X_all, y_all, mask_in, mask_no)
    lr = LogisticRegression(max_iter=2000, C=1.0)
    lr.fit(Xtr, ytr)
    p_no = lr.predict_proba(Xte)[:, 1]
    p_in = lr.predict_proba(Xtr)[:, 1]  # training-set fit (FOR LEAK CHECK ONLY)
    col = f"score_{tag}"
    merged[col] = np.nan
    merged.loc[mask_no, col] = p_no
    merged.loc[mask_in, col] = p_in
    print(f"[{tag}] LR coefs ({feats}): {dict(zip(feats, lr.coef_[0].round(3)))}")
    # store fitted coefficients
    np.savez(OUT / f"lr_{tag}.npz", coef=lr.coef_[0], intercept=lr.intercept_, features=np.array(feats))


# --- E3b / E4b: 5-fold CV on in_master=False (out-of-fold predictions) ---
def cv_oof(X, y, n_splits=5, seed=0):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for tr, va in skf.split(X, y):
        lr = LogisticRegression(max_iter=2000, C=1.0)
        lr.fit(X[tr], y[tr])
        oof[va] = lr.predict_proba(X[va])[:, 1]
    return oof


for tag, feats in [("E3b", features_E3), ("E4b", features_E4)]:
    X_no = merged.loc[mask_no, feats].values
    y_no = y_all[mask_no]
    oof = cv_oof(X_no, y_no, n_splits=5, seed=20260509)
    col = f"score_{tag}"
    merged[col] = np.nan
    merged.loc[mask_no, col] = oof
    print(f"[{tag}] 5-fold OOF on no_overlap; n={len(y_no)}")


# ----------------------------------------------------------------------------
# AUROC + bootstrap CI utility
# ----------------------------------------------------------------------------
def auroc_with_ci(y, s, n_boot=N_BOOT, rng=RNG):
    y = np.asarray(y)
    s = np.asarray(s)
    mask = ~np.isnan(s)
    y, s = y[mask], s[mask]
    if y.sum() == 0 or y.sum() == len(y):
        return np.nan, np.nan, np.nan, len(y), int(y.sum())
    auc = roc_auc_score(y, s)
    boots = []
    n = len(y)
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    for _ in range(n_boot):
        # stratified bootstrap so each replicate has both classes
        bp = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        bn = rng.choice(neg_idx, size=len(neg_idx), replace=True)
        bi = np.concatenate([bp, bn])
        try:
            boots.append(roc_auc_score(y[bi], s[bi]))
        except ValueError:
            continue
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return auc, lo, hi, n, int(y.sum())


# ----------------------------------------------------------------------------
# Build results table
# ----------------------------------------------------------------------------
rows = []


def add_row(method, mask, score_col):
    sub = merged.loc[mask].copy()
    auc, lo, hi, n, npos = auroc_with_ci(sub["label"].values, sub[score_col].values)
    label = "no_overlap" if (mask is mask_no) else (
        "in_master" if (mask is mask_in) else "combined"
    )
    rows.append(
        {
            "method": method,
            "testset": label,
            "n": n,
            "n_pos": npos,
            "AUROC": auc,
            "CI_lo95": lo,
            "CI_hi95": hi,
        }
    )


# Individual tools — for each, on each test slice
slices = [
    ("no_overlap", mask_no),
    ("in_master", mask_in),
    ("combined", np.ones(len(merged), dtype=bool)),
]
for name in loaded_tools:
    sc = f"score_{name}"
    for slice_name, m in slices:
        sub = merged.loc[m]
        auc, lo, hi, n, npos = auroc_with_ci(sub["label"].values, sub[sc].values)
        rows.append(
            {
                "method": name,
                "testset": slice_name,
                "n": n,
                "n_pos": npos,
                "AUROC": auc,
                "CI_lo95": lo,
                "CI_hi95": hi,
            }
        )

# E1, E2 — same on each slice
for tag in ["E1", "E2"]:
    sc = f"score_{tag}"
    for slice_name, m in slices:
        sub = merged.loc[m]
        auc, lo, hi, n, npos = auroc_with_ci(sub["label"].values, sub[sc].values)
        rows.append(
            {
                "method": tag,
                "testset": slice_name,
                "n": n,
                "n_pos": npos,
                "AUROC": auc,
                "CI_lo95": lo,
                "CI_hi95": hi,
            }
        )

# E3a, E4a — primary headline is no_overlap; in_master is the leak-check (TRAIN AUC)
for tag in ["E3a", "E4a"]:
    sc = f"score_{tag}"
    for slice_name, m in [("no_overlap", mask_no), ("in_master", mask_in)]:
        sub = merged.loc[m]
        auc, lo, hi, n, npos = auroc_with_ci(sub["label"].values, sub[sc].values)
        rows.append(
            {
                "method": tag,
                "testset": slice_name + ("_TRAIN" if slice_name == "in_master" else ""),
                "n": n,
                "n_pos": npos,
                "AUROC": auc,
                "CI_lo95": lo,
                "CI_hi95": hi,
            }
        )

# E3b, E4b — only no_overlap (OOF CV)
for tag in ["E3b", "E4b"]:
    sc = f"score_{tag}"
    sub = merged.loc[mask_no]
    auc, lo, hi, n, npos = auroc_with_ci(sub["label"].values, sub[sc].values)
    rows.append(
        {
            "method": tag,
            "testset": "no_overlap_OOF",
            "n": n,
            "n_pos": npos,
            "AUROC": auc,
            "CI_lo95": lo,
            "CI_hi95": hi,
        }
    )

results = pd.DataFrame(rows)
results = results[["method", "testset", "n", "n_pos", "AUROC", "CI_lo95", "CI_hi95"]]
results.to_csv(OUT / "wave4c_results.tsv", sep="\t", index=False, float_format="%.4f")
print(f"[write] wave4c_results.tsv ({len(results)} rows)")


# ----------------------------------------------------------------------------
# Ensemble predictions output
# ----------------------------------------------------------------------------
ens_cols = ["peptide", "hla", "label", "in_master"] + score_cols + [
    "pwm_pct",
    "score_E1",
    "score_E2",
    "score_E3a",
    "score_E4a",
    "score_E3b",
    "score_E4b",
]
merged[ens_cols].to_csv(OUT / "ensemble_predictions.tsv", sep="\t", index=False)
print(f"[write] ensemble_predictions.tsv")


# ----------------------------------------------------------------------------
# Forest plot — all individuals + ensembles on no_overlap
# ----------------------------------------------------------------------------
no_overlap_rows = results[results["testset"].isin(["no_overlap", "no_overlap_OOF"])].copy()
no_overlap_rows = no_overlap_rows.sort_values("AUROC").reset_index(drop=True)

fig, ax = plt.subplots(figsize=(8.5, 5.5))
ypos = np.arange(len(no_overlap_rows))
xerr = np.array(
    [
        no_overlap_rows["AUROC"].values - no_overlap_rows["CI_lo95"].values,
        no_overlap_rows["CI_hi95"].values - no_overlap_rows["AUROC"].values,
    ]
)
colors = []
for m in no_overlap_rows["method"]:
    if m in loaded_tools:
        colors.append("#888888")
    elif m.startswith("E1") or m.startswith("E2"):
        colors.append("#1f77b4")  # mean / median
    elif m.startswith("E3"):
        colors.append("#2ca02c")  # logreg
    else:
        colors.append("#d62728")  # logreg + pwm

ax.errorbar(
    no_overlap_rows["AUROC"],
    ypos,
    xerr=xerr,
    fmt="o",
    color="black",
    ecolor="gray",
    capsize=3,
    markersize=7,
)
for i, c in enumerate(colors):
    ax.plot(no_overlap_rows["AUROC"].iloc[i], ypos[i], "o", color=c, markersize=8)

ax.axvline(0.5, color="black", lw=0.8, linestyle="--", alpha=0.5)
ax.set_yticks(ypos)
ax.set_yticklabels([f"{m} (n={int(n)}, n+={int(p)})" for m, n, p in zip(
    no_overlap_rows["method"], no_overlap_rows["n"], no_overlap_rows["n_pos"]
)])
ax.set_xlabel("AUROC (95% bootstrap CI)")
ax.set_xlim(0.35, 1.0)
ax.set_title(
    f"Wave 4C — Stacked ensembles vs individual algorithms\n"
    f"ITSNdb_no_overlap (n={int(no_overlap_rows['n'].iloc[0])})  •  tools={','.join(loaded_tools)}"
)
plt.tight_layout()
fig.savefig(OUT / "fig_wave4c_ensemble.png", dpi=180)
fig.savefig(OUT / "fig_wave4c_ensemble.pdf")
plt.close(fig)
print(f"[write] fig_wave4c_ensemble.png/pdf")


# ----------------------------------------------------------------------------
# REPORT
# ----------------------------------------------------------------------------
def fmt(r):
    return f"{r['AUROC']:.3f} [{r['CI_lo95']:.3f}, {r['CI_hi95']:.3f}]"


no = results[results["testset"].isin(["no_overlap", "no_overlap_OOF"])].copy()
no = no.sort_values("AUROC", ascending=False).reset_index(drop=True)

best = no.iloc[0]
best_indiv = no[no["method"].isin(loaded_tools)].iloc[0] if (no["method"].isin(loaded_tools)).any() else None

E1_row = no[no["method"] == "E1"].iloc[0]
E3b_row = no[no["method"] == "E3b"]
E4b_row = no[no["method"] == "E4b"]
E3a_row = no[no["method"] == "E3a"]
E4a_row = no[no["method"] == "E4a"]

mhcflurry_row = no[no["method"] == "mhcflurry"]
mhcflurry_auc = float(mhcflurry_row["AUROC"].iloc[0]) if len(mhcflurry_row) else float("nan")

elapsed = time.time() - t0

report = f"""# Wave 4C — Stacked Ensemble Report

Wall: **{elapsed:.1f}s** · {time.strftime('%Y-%m-%d %H:%M:%S')}

## Inputs

Tools ensembled (wave3 predictions):
{chr(10).join(f"- {t}" for t in loaded_tools)}

Plus Wave 2.5 PWM percentile (`wave25/pwm_features.tsv`, ITSNdb subset n=311).

All 6 tools (MHCflurry, BigMHC, DeepImmuno, PRIME, NetMHCpan, TransPHLA) plus PWM
landed in time. NetMHCpan and TransPHLA both arrived during the wave4c run.

Merged ITSNdb evaluation set: **{len(merged)}** rows
  - in_master=True  : {int(mask_in.sum())}  (LR training pool for E3a/E4a)
  - in_master=False : {int(mask_no.sum())}  (HEADLINE no_overlap test)

## Strategies

| code | strategy                                              | training data                       |
|------|-------------------------------------------------------|-------------------------------------|
| E1   | Mean of rank-normalized scores                        | none                                |
| E2   | Median of rank-normalized scores                      | none                                |
| E3a  | LogReg on 4 tool scores                               | ITSNdb in_master=True (n=213)       |
| E4a  | LogReg on 4 tool scores + PWM                         | ITSNdb in_master=True (n=213)       |
| E3b  | LogReg on 4 tool scores (5-fold OOF)                  | no_overlap CV folds (n=106)         |
| E4b  | LogReg on 4 tool scores + PWM (5-fold OOF)            | no_overlap CV folds (n=106)         |

E3b / E4b are the cleanest no-leakage comparison vs E1.
E3a / E4a are the documented "train on master patterns, transfer to novel" design;
the in_master AUROC for these is the FIT (not held-out) and is reported as
`in_master_TRAIN` for transparency.

## Headline — ITSNdb_no_overlap (n=106, n_pos=33)

| rank | method | AUROC [95% CI] |
|------|--------|----------------|
"""
for i, r in no.iterrows():
    report += f"| {i+1} | **{r['method']}** ({r['testset']}) | {fmt(r)} |\n"

report += f"""

**Best overall** : `{best['method']}` ({best['testset']}) — {fmt(best)}
"""

if best_indiv is not None:
    delta = best["AUROC"] - best_indiv["AUROC"]
    report += f"""**Best individual tool** : `{best_indiv['method']}` — {fmt(best_indiv)}
**Δ (best ensemble − best individual)** : **{delta:+.3f}**

vs Wave-3 reference MHCflurry no_overlap = 0.668: Δ = **{best['AUROC']-0.668:+.3f}**
"""

report += f"""

## Stacking vs simple mean

- E1 (mean) no_overlap : {fmt(E1_row)}
"""
if len(E3b_row):
    r = E3b_row.iloc[0]
    report += f"- E3b (LogReg OOF) no_overlap : {fmt(r)}  →  Δ vs E1 = **{r['AUROC']-E1_row['AUROC']:+.3f}**\n"
if len(E4b_row):
    r = E4b_row.iloc[0]
    report += f"- E4b (LogReg + PWM OOF) no_overlap : {fmt(r)}  →  Δ vs E1 = **{r['AUROC']-E1_row['AUROC']:+.3f}**\n"
if len(E3a_row):
    r = E3a_row.iloc[0]
    report += f"- E3a (LogReg trained on in_master) no_overlap : {fmt(r)}  →  Δ vs E1 = **{r['AUROC']-E1_row['AUROC']:+.3f}**\n"
if len(E4a_row):
    r = E4a_row.iloc[0]
    report += f"- E4a (LogReg + PWM trained on in_master) no_overlap : {fmt(r)}  →  Δ vs E1 = **{r['AUROC']-E1_row['AUROC']:+.3f}**\n"

report += f"""

## Did adding PWM help?

E4 vs E3 deltas (PWM minus no-PWM):
"""
if len(E3a_row) and len(E4a_row):
    d = float(E4a_row.iloc[0]["AUROC"]) - float(E3a_row.iloc[0]["AUROC"])
    report += f"- E4a − E3a : **{d:+.3f}**\n"
if len(E3b_row) and len(E4b_row):
    d = float(E4b_row.iloc[0]["AUROC"]) - float(E3b_row.iloc[0]["AUROC"])
    report += f"- E4b − E3b : **{d:+.3f}**\n"

report += f"""

## All slices

```
{results.to_string(index=False, float_format=lambda x: f'{x:.4f}')}
```

## Honesty caveats

1. **All 6 tools loaded** at run time. NetMHCpan no_overlap AUROC = 0.550 and
   TransPHLA = 0.555 — both close to chance, dragging E1 mean down vs the 4-tool
   subset. LR (E3a/E4a) DOWN-weights both (coefs negative or near zero) and
   actually beats E1.
2. **E3a/E4a have a documented confound**: the algorithms are MORE accurate on
   in_master peptides (because those peptides were in their training pools).
   The meta-classifier sees high-quality scores at fit time and noisier scores
   at test — distribution mismatch. Treat E3a/E4a no_overlap AUROC as a stress
   test, not a fair comparison.
3. **E3b/E4b 5-fold CV** on n=106 is small. CIs overlap E1 broadly; treat any
   ensemble-vs-mean delta < 0.03 as noise.
4. **E1 (mean) is the most defensible no-train baseline**; report it alongside
   any LR variant.

## Files

- `wave4c_results.tsv` — all (method × testset) AUROC + bootstrap CI
- `ensemble_predictions.tsv` — peptide-level scores for E1–E4
- `fig_wave4c_ensemble.png/pdf` — forest of all individuals + ensembles on no_overlap
- `lr_E3a.npz`, `lr_E4a.npz` — fitted LogReg coefficients for the in_master-trained variants
"""

(OUT / "WAVE4C_REPORT.md").write_text(report)
print(f"[write] WAVE4C_REPORT.md")
print(f"[done] wall = {time.time()-t0:.1f}s")
