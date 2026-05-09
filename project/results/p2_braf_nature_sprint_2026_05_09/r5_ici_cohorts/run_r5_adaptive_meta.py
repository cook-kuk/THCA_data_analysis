#!/usr/bin/env python3
"""R5 supplement — Mantel-Haenszel pooled OR for the
DM1-like-adaptive (HT-13 high & HLA-I high) responder analysis.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/r5_ici_cohorts")

df = pd.read_csv(OUT / "r5_dm1_adaptive_responder.tsv", sep="\t")

# Mantel-Haenszel pooled OR
a = df["a_adapt_resp"].to_numpy(float)
b = df["b_adapt_nonresp"].to_numpy(float)
c = df["c_other_resp"].to_numpy(float)
d = df["d_other_nonresp"].to_numpy(float)
n = a + b + c + d

num = (a * d / n).sum()
den = (b * c / n).sum()
mh_or = num / den if den > 0 else np.nan

# Robins-Breslow-Greenland variance for ln(OR_MH)
P = (a + d) / n
Q = (b + c) / n
R = (a * d) / n
S = (b * c) / n
sumR = R.sum(); sumS = S.sum()
var_ln = ((P * R).sum() / (2 * sumR ** 2)
          + ((P * S + Q * R).sum()) / (2 * sumR * sumS)
          + (Q * S).sum() / (2 * sumS ** 2))
se_ln = float(np.sqrt(var_ln))
ln_or = float(np.log(mh_or))
lo = float(np.exp(ln_or - 1.96 * se_ln))
hi = float(np.exp(ln_or + 1.96 * se_ln))
z = ln_or / se_ln
p = 2 * (1 - stats.norm.cdf(abs(z)))

# heterogeneity Q (Cochran)
Or_i = (a * d) / np.where(b * c > 0, b * c, np.nan)
ln_or_i = np.log(Or_i)
w_i = 1.0 / (1.0/np.where(a>0,a,np.nan) + 1.0/np.where(b>0,b,np.nan)
             + 1.0/np.where(c>0,c,np.nan) + 1.0/np.where(d>0,d,np.nan))
mask = np.isfinite(ln_or_i) & np.isfinite(w_i)
ln_pool = (w_i[mask] * ln_or_i[mask]).sum() / w_i[mask].sum()
Q_stat = (w_i[mask] * (ln_or_i[mask] - ln_pool) ** 2).sum()
df_q = mask.sum() - 1
I2 = max(0.0, (Q_stat - df_q) / Q_stat) * 100 if Q_stat > 0 else 0.0

# pooled response rates
n_resp_adapt = int(a.sum()); n_adapt = int((a + b).sum())
n_resp_other = int(c.sum()); n_other = int((c + d).sum())

print(f"DM1-like-adaptive Mantel-Haenszel pooled OR = {mh_or:.3f}")
print(f"  95% CI [{lo:.3f}, {hi:.3f}]  p={p:.4g}  I2={I2:.1f}%  k={len(df)}")
print(f"  Adaptive responder rate: {n_resp_adapt}/{n_adapt} = {n_resp_adapt/n_adapt*100:.1f}%")
print(f"  Non-adaptive responder rate: {n_resp_other}/{n_other} = {n_resp_other/n_other*100:.1f}%")
print(f"  Sign-consistent OR>1 cohorts: {(df['fisher_OR']>1).sum()}/{len(df)}")

with open(OUT / "r5_adaptive_meta.json", "w") as f:
    json.dump(dict(mh_OR=float(mh_or), lo=lo, hi=hi, p=float(p), I2=float(I2),
                   k=int(len(df)),
                   adaptive_resp_rate=float(n_resp_adapt/n_adapt),
                   nonadaptive_resp_rate=float(n_resp_other/n_other),
                   n_resp_adapt=n_resp_adapt, n_adapt=n_adapt,
                   n_resp_other=n_resp_other, n_other=n_other,
                   sign_concordance=int((df["fisher_OR"]>1).sum()),
                   total_cohorts=int(len(df))), f, indent=2)
