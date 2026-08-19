"""
Survival comparison — TSO500-3 (TSHR · PAX8 · NKX2-1) vs Full-8 panel.

Question: 만약 분당 코호트에 TSO500 v2 (3-gene) 만 측정 가능하면,
        그 점수만으로도 환자 OS / PFI / DSS 가 갈리는가?

Endpoints: OS, PFI, DSS (TCGA-THCA pan-cancer clinical, Liu 2018)
Tests:
  - Median split (high vs low) — log-rank + Cox HR
  - Tertiles (T1 vs T3)
  - Continuous score Cox (standardized)
"""
import os, json, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

# Korean font
for f in ["/usr/share/fonts/truetype/nanum/NanumGothic.ttf",
          "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"]:
    if os.path.exists(f):
        font_manager.fontManager.addfont(f)
        break
plt.rcParams["font.family"] = ["NanumGothic", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = "/home/seungho/personal/THCA_data_analysis/project"
OUT  = f"{ROOT}/dm1_story_web/public/figures"

GENES_FULL   = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]
GENES_TSO500 = ["TSHR","PAX8","NKX2-1"]
GENES_C5     = ["PAX8","NKX2-1","TSHR","TPO","DIO1"]

# Load β + clinical
hm = pd.read_csv(f"{ROOT}/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
cl = pd.read_csv(f"{ROOT}/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")

df = hm.merge(cl, left_on="sample_short", right_on="tcga_short", how="inner")
df = df.dropna(subset=GENES_FULL)
print(f"Merged cohort n = {len(df)}")

# Compute panel scores
df["score_full8"]  = df[GENES_FULL].mean(axis=1)
df["score_tso500"] = df[GENES_TSO500].mean(axis=1)
df["score_c5"]     = df[GENES_C5].mean(axis=1)

# Sign convention: higher score → more methylated → DM1-like (bad).
# Round 4 r9_1 master: DM1 = HIGH-β / silenced (cPTC d = +1.50).
# So we'll use "high score = DM1-like" convention.  If HR ends up <1 we'll
# flip the score interpretation in the report.
PANELS = {
    "Full-8":      ("score_full8",   "#0F172A"),
    "TSO500-3":    ("score_tso500",  "#B91C1C"),
    "Compact-5":   ("score_c5",      "#0E7490")
}
ENDPOINTS = {
    "OS":  ("OS.time",  "OS",  "Overall Survival"),
    "PFI": ("PFI.time", "PFI", "Progression-Free Interval"),
    "DSS": ("DSS.time", "DSS", "Disease-Specific Survival")
}

# ─── Statistics ───
def cox_hr(scores, t, e):
    m = pd.DataFrame({"score": scores, "t": t, "e": e}).dropna()
    m = m[m["t"] > 0]
    if m["e"].sum() < 3: return np.nan, np.nan, np.nan, np.nan
    m["z"] = (m["score"] - m["score"].mean()) / m["score"].std()
    cph = CoxPHFitter()
    try:
        cph.fit(m[["t","e","z"]], duration_col="t", event_col="e")
        hr = float(np.exp(cph.params_["z"]))
        p  = float(cph.summary.loc["z","p"])
        ci_lo = float(np.exp(cph.confidence_intervals_.loc["z"].iloc[0]))
        ci_hi = float(np.exp(cph.confidence_intervals_.loc["z"].iloc[1]))
        return hr, ci_lo, ci_hi, p
    except Exception as ex:
        return np.nan, np.nan, np.nan, np.nan

def logrank_split(scores, t, e, split="median"):
    m = pd.DataFrame({"s": scores, "t": t, "e": e}).dropna()
    m = m[m["t"] > 0]
    if split == "median":
        hi = m[m["s"] > m["s"].median()]
        lo = m[m["s"] <= m["s"].median()]
    elif split == "tertile":
        q1, q3 = m["s"].quantile([1/3, 2/3])
        hi = m[m["s"] >= q3]; lo = m[m["s"] <= q1]
    r = logrank_test(hi["t"], lo["t"], event_observed_A=hi["e"], event_observed_B=lo["e"])
    return r.p_value, len(hi), len(lo)

# ─── Compute all metrics ───
results = []
for ep_name, (tcol, ecol, ep_label) in ENDPOINTS.items():
    n_event = int(df[ecol].sum())
    n_total = int(df[ecol].notna().sum())
    for p_name, (score_col, _) in PANELS.items():
        hr, lo, hi, p = cox_hr(df[score_col], df[tcol], df[ecol])
        lr_med_p, n_hi_m, n_lo_m = logrank_split(df[score_col], df[tcol], df[ecol], "median")
        lr_ter_p, n_hi_t, n_lo_t = logrank_split(df[score_col], df[tcol], df[ecol], "tertile")
        results.append({
            "endpoint": ep_name, "ep_label": ep_label, "panel": p_name,
            "n_total": n_total, "n_events": n_event,
            "hr_per_sd": hr, "hr_ci_lo": lo, "hr_ci_hi": hi, "cox_p": p,
            "logrank_median_p": lr_med_p,
            "logrank_tertile_p": lr_ter_p
        })

df_res = pd.DataFrame(results)
df_res.to_csv(f"{OUT}/tso500_survival_compare.tsv", sep="\t", index=False)
with open(f"{OUT}/tso500_survival_compare.json", "w") as f:
    json.dump(df_res.to_dict(orient="records"), f, indent=2)

print("\n=== Survival comparison (HR per +1 SD score) ===")
print(f"{'Endpoint':<6s}  {'Panel':<10s}  {'n':>5s}/{'ev':<4s}  {'HR':>6s}  {'95% CI':>14s}  {'Cox p':>10s}  {'LR-med p':>10s}  {'LR-tert p':>10s}")
for r in results:
    print(f"  {r['endpoint']:<4s}  {r['panel']:<10s}  {r['n_total']:5d}/{r['n_events']:<4d}  {r['hr_per_sd']:6.2f}  [{r['hr_ci_lo']:.2f}-{r['hr_ci_hi']:.2f}]  {r['cox_p']:10.3e}  {r['logrank_median_p']:10.3e}  {r['logrank_tertile_p']:10.3e}")

# ─── KM figure: 3 endpoints × 3 panels (median split) ───
fig, axes = plt.subplots(3, 3, figsize=(15, 13), dpi=160)
for i, (ep_name, (tcol, ecol, ep_label)) in enumerate(ENDPOINTS.items()):
    for j, (p_name, (score_col, color)) in enumerate(PANELS.items()):
        ax = axes[i, j]
        m = df[[score_col, tcol, ecol]].dropna()
        m = m[m[tcol] > 0]
        med = m[score_col].median()
        hi = m[m[score_col] >  med]
        lo = m[m[score_col] <= med]
        kmf_hi = KaplanMeierFitter().fit(hi[tcol]/365.25, hi[ecol], label=f"High score (n={len(hi)})")
        kmf_lo = KaplanMeierFitter().fit(lo[tcol]/365.25, lo[ecol], label=f"Low score  (n={len(lo)})")
        kmf_hi.plot_survival_function(ax=ax, color=color, ci_alpha=0.12, lw=2)
        kmf_lo.plot_survival_function(ax=ax, color="#94A3B8", ci_alpha=0.12, lw=2)
        # Stats annotation
        r = next(x for x in results if x["endpoint"]==ep_name and x["panel"]==p_name)
        ax.text(0.04, 0.10,
                f"HR = {r['hr_per_sd']:.2f}  [{r['hr_ci_lo']:.2f}-{r['hr_ci_hi']:.2f}]\n"
                f"Cox p = {r['cox_p']:.2e}\n"
                f"Log-rank (med) p = {r['logrank_median_p']:.2e}",
                transform=ax.transAxes, fontsize=8.5, family="monospace",
                bbox=dict(facecolor="white", edgecolor="#CBD5E1", boxstyle="round,pad=0.4"))
        ax.set_xlabel("Years"); ax.set_ylabel("Survival probability")
        ax.set_ylim(0.55, 1.02); ax.set_xlim(0, 15)
        ax.set_title(f"{ep_label} · {p_name}", fontsize=11.5, fontweight="bold", color=color, loc="left")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.legend(loc="lower left", fontsize=8.5, frameon=False)
fig.suptitle("TSO500-3 vs Full-8 vs Compact-5 — TCGA-THCA Kaplan-Meier (median split)",
             fontsize=14, fontweight="bold", y=1.00)
plt.tight_layout()
fig.savefig(f"{OUT}/fig_tso500_survival_km.png", bbox_inches="tight")
plt.close(fig)

# ─── Forest comparison ───
fig, ax = plt.subplots(figsize=(11, 5.5), dpi=170)
ys = []; labels = []; hrs = []; los = []; his = []; ps = []; colors = []
y = 0
for ep_name, (tcol, ecol, ep_label) in ENDPOINTS.items():
    for p_name, (score_col, color) in PANELS.items():
        r = next(x for x in results if x["endpoint"]==ep_name and x["panel"]==p_name)
        ys.append(y); labels.append(f"{ep_name}  ·  {p_name}")
        hrs.append(r["hr_per_sd"]); los.append(r["hr_ci_lo"]); his.append(r["hr_ci_hi"])
        ps.append(r["cox_p"]); colors.append(color)
        y -= 1
    y -= 0.4
ys = np.array(ys); hrs = np.array(hrs); los = np.array(los); his = np.array(his)
ax.errorbar(hrs, ys, xerr=[hrs-los, his-hrs], fmt="o", markersize=10,
            ecolor="#94A3B8", elinewidth=1.5, capsize=4, zorder=2)
for i in range(len(ys)):
    ax.plot([hrs[i]], [ys[i]], "o", markersize=10, markerfacecolor=colors[i],
            markeredgecolor="white", markeredgewidth=1.5, zorder=3)
ax.axvline(1.0, color="#0F172A", ls="--", lw=1, alpha=0.5)
ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("Hazard Ratio per +1 SD score  (95 % CI)", fontsize=11)
ax.set_xscale("log"); ax.set_xlim(0.3, 5)
for i, (h, p, lab) in enumerate(zip(hrs, ps, labels)):
    txt = f"HR = {h:.2f}, p = {p:.1e}" if not np.isnan(h) else "n.a."
    ax.text(his[i]*1.18 if not np.isnan(h) else 1.5, ys[i], txt, va="center", fontsize=9, family="monospace")
ax.set_title("Figure SR-1 · TSO500-3 vs Full-8 vs Compact-5 — 생존 분석 forest",
             fontsize=12.5, fontweight="bold", loc="left", pad=12)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
plt.tight_layout()
fig.savefig(f"{OUT}/fig_tso500_survival_forest.png", bbox_inches="tight")
plt.close(fig)

print(f"\nOutputs to {OUT}:")
for f in sorted(os.listdir(OUT)):
    if "tso500_survival" in f: print(f"  {f}")
