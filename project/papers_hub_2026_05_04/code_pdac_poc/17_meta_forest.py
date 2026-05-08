"""Meta-analysis forest plot with pooled HR diamond."""
import json, math
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/data/pdac_poc"); FIG = ROOT/"figures"
m = json.load(open(ROOT/"results/method_panel/META_ANALYSIS.json"))
rows = [r for r in m["rows"] if r.get("G12D_HR") is not None and r.get("G12D_se_log_HR")]
rows.sort(key=lambda r: -(r.get("n_with_OS") or 0))
meta = m["meta_analysis"]

fig, ax = plt.subplots(figsize=(8.6, 5.0), dpi=140)
y = np.arange(len(rows)+1)
labels = [r["label"] for r in rows] + [f"★ POOLED (n={meta['n_total_OS']}, k={meta['n_studies']})"]
hrs = [r["G12D_HR"] for r in rows] + [meta["pooled_HR"]]
los = [math.exp(math.log(r["G12D_HR"]) - 1.96*r["G12D_se_log_HR"]) for r in rows] + [meta["CI_lo"]]
his = [math.exp(math.log(r["G12D_HR"]) + 1.96*r["G12D_se_log_HR"]) for r in rows] + [meta["CI_hi"]]
ps  = [r["G12D_p"] for r in rows] + [meta["p_value"]]
ns  = [r.get("n_with_OS") for r in rows] + [meta["n_total_OS"]]

for i, (hr, lo, hi, p, n, lbl) in enumerate(zip(hrs, los, his, ps, ns, labels)):
    is_pool = i == len(rows)
    color = "#7B1F2A" if (p is not None and p < 0.05) else "#888"
    if is_pool: color = "#1c8e6d"
    ax.plot([lo, hi], [i, i], color=color, lw=3 if is_pool else 2, alpha=0.92)
    if is_pool:
        ax.plot([lo, hr, hi, hr, lo], [i-.18, i+.22, i-.18, i-.36, i-.18],
                "-", color=color, lw=2.5)
        ax.fill([lo, hr, hi, hr], [i, i+.22, i, i-.22], color=color, alpha=0.5)
    else:
        ax.plot([hr], [i], "s", color=color, markersize=10+min(7, math.log10(max(n,1))*4))
    psstr = f"p={p:.3g}" if p else ""
    ax.text(hi + 0.06, i, f"  HR={hr:.2f}  {psstr}  n={n}",
            va="center", fontsize=9.5 if is_pool else 9,
            color="#0a0a0a", fontweight="bold" if is_pool else "normal")
ax.axvline(1.0, color="#bbb", lw=0.8, linestyle="--")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=10)
ax.invert_yaxis()
ax.set_xlabel("KRAS G12D Cox HR (95% CI)")
ax.set_title("Cross-cohort meta-analysis · KRAS G12D HR in PDAC")
ax.set_xlim(0.6, 4.3)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(FIG/"fig_g12d_meta_forest.png", dpi=160, bbox_inches="tight")
fig.savefig(FIG/"fig_g12d_meta_forest.svg", bbox_inches="tight")
print("[17] meta-forest plot saved.")
print(f"     pooled HR = {meta['pooled_HR']}  CI {meta['CI_lo']}–{meta['CI_hi']}  p = {meta['p_value']}  k={meta['n_studies']}, n={meta['n_total_OS']}")
