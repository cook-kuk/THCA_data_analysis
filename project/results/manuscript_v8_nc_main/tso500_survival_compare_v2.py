"""
Survival comparison v2 — 더 깊은 비교:
  (A) Continuous-score Cox (Full-8 vs TSO500-3 vs Compact-5)
  (B) DM1-derived label by top-quantile score → KM
  (C) Stage III/IV only (event 더 많은 subset)
  (D) TERT-combined stratification (4-way)

Reports both raw and stage-adjusted Cox HR.
"""
import os, json, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

for f in ["/usr/share/fonts/truetype/nanum/NanumGothic.ttf"]:
    if os.path.exists(f):
        font_manager.fontManager.addfont(f); break
plt.rcParams["font.family"] = ["NanumGothic","DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = "/home/seungho/personal/THCA_data_analysis/project"
OUT  = f"{ROOT}/dm1_story_web/public/figures"

GENES_FULL   = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]
GENES_TSO500 = ["TSHR","PAX8","NKX2-1"]
GENES_C5     = ["PAX8","NKX2-1","TSHR","TPO","DIO1"]

hm = pd.read_csv(f"{ROOT}/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
cl = pd.read_csv(f"{ROOT}/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
df = hm.merge(cl, left_on="sample_short", right_on="tcga_short", how="inner").dropna(subset=GENES_FULL)
df["score_full8"]  = df[GENES_FULL].mean(axis=1)
df["score_tso500"] = df[GENES_TSO500].mean(axis=1)
df["score_c5"]     = df[GENES_C5].mean(axis=1)

# Sign check — Round 4 r9_1: DM1 = HIGH β (silenced).  So high score → DM1-like → bad prognosis.
# Verify direction: Cox HR per +1 SD should be ≥ 1 if convention holds.
print(f"Cohort n={len(df)} · OS events={int(df['OS'].sum())} · PFI events={int(df['PFI'].sum())} · DSS events={int(df['DSS'].sum())}")

# Stage info
df["stage_high"] = df["stage"].astype(str).str.contains("III|IV", regex=True, na=False).astype(int)
print(f"  Stage III/IV: n={int(df['stage_high'].sum())}")

# ── Helper: derive DM1-like label by top-quantile score (sign-aware) ──
def derive_dm1_label(scores, frac=0.30):
    """top frac (most DM1-like by score) = 1, rest = 0"""
    thr = scores.quantile(1 - frac)
    return (scores >= thr).astype(int)

# Truth DM1 label (from master_tcga) — need to load
master = pd.read_csv(f"{ROOT}/results/manuscript_v8_nc_main/master_tcga.tsv", sep="\t")[["sample_short","DM"]]
df = df.merge(master, on="sample_short", how="left")
df["DM1_truth"] = (df["DM"]=="DM1").astype(int)
print(f"  DM1 truth: n_DM1={int(df['DM1_truth'].sum())} / total={df['DM1_truth'].notna().sum()}")

# ── Align score direction with truth so HR>1 = bad ──
# If high score correlates with DM1=1, keep; else flip.
for col in ["score_full8","score_tso500","score_c5"]:
    m = df.dropna(subset=[col,"DM1_truth"])
    r = np.corrcoef(m[col], m["DM1_truth"])[0,1]
    if r < 0:
        df[col] = -df[col]
        print(f"  flipped {col} (was r={r:.3f} with DM1)")
    else:
        print(f"  kept    {col} (r={r:.3f} with DM1)")

PANELS = {
    "DM1_truth_label": ("DM1_truth",   "#7C3AED", "category"),
    "Full-8_score":    ("score_full8", "#0F172A", "score"),
    "TSO500-3_score":  ("score_tso500","#B91C1C", "score"),
    "Compact-5_score": ("score_c5",    "#0E7490", "score")
}
ENDPOINTS = {
    "OS":  ("OS.time","OS","Overall Survival"),
    "PFI": ("PFI.time","PFI","Progression-Free Interval"),
    "DSS": ("DSS.time","DSS","Disease-Specific Survival")
}

def cox_hr(scores, t, e, adj=None):
    m = pd.DataFrame({"s":scores,"t":t,"e":e});
    if adj is not None: m["adj"] = adj.values
    m = m.dropna()
    m = m[m["t"]>0]
    if m["e"].sum() < 3: return np.nan,np.nan,np.nan,np.nan,len(m),int(m["e"].sum())
    cph = CoxPHFitter()
    try:
        cph.fit(m, duration_col="t", event_col="e")
        hr = float(np.exp(cph.params_["s"]))
        p  = float(cph.summary.loc["s","p"])
        ci_lo = float(np.exp(cph.confidence_intervals_.loc["s"].iloc[0]))
        ci_hi = float(np.exp(cph.confidence_intervals_.loc["s"].iloc[1]))
        return hr,ci_lo,ci_hi,p,len(m),int(m["e"].sum())
    except Exception: return np.nan,np.nan,np.nan,np.nan,len(m),int(m["e"].sum())

def logrank_by_label(label, t, e):
    m = pd.DataFrame({"l":label,"t":t,"e":e}).dropna()
    m = m[m["t"]>0]
    g1 = m[m["l"]==1]; g0 = m[m["l"]==0]
    if g1["e"].sum() < 2 or g0["e"].sum() < 2: return np.nan
    r = logrank_test(g1["t"], g0["t"], g1["e"], g0["e"])
    return float(r.p_value)

# ─── (A) Continuous-score Cox per endpoint (full cohort + stage-adjusted + high-stage only) ───
rows = []
for ep, (tcol,ecol,lab) in ENDPOINTS.items():
    for p_name, (col, color, ptype) in PANELS.items():
        scores = df[col] if ptype=="score" else df[col].astype(float)
        # Standardize score panels
        s = (scores-scores.mean())/scores.std() if ptype=="score" else scores
        # Raw
        hr,lo,hi,p,n,ev = cox_hr(s, df[tcol], df[ecol])
        # Stage adjusted
        hr_adj,lo_adj,hi_adj,p_adj,_,_ = cox_hr(s, df[tcol], df[ecol], adj=df["stage_high"])
        # High-stage only
        d_hi = df[df["stage_high"]==1]
        s_hi = (d_hi[col]-d_hi[col].mean())/d_hi[col].std() if ptype=="score" else d_hi[col].astype(float)
        hr_h,lo_h,hi_h,p_h,n_h,ev_h = cox_hr(s_hi, d_hi[tcol], d_hi[ecol])
        # Median-split log-rank (use derived label or truth)
        if ptype=="score":
            label = derive_dm1_label(df[col], frac=0.30)
        else:
            label = df[col].astype(int)
        lr_p = logrank_by_label(label, df[tcol], df[ecol])
        rows.append({
            "endpoint":ep,"ep_label":lab,"panel":p_name,
            "n_total":n,"n_events":ev,
            "hr_per_sd":hr,"hr_ci_lo":lo,"hr_ci_hi":hi,"cox_p":p,
            "hr_stageadj":hr_adj,"cox_p_stageadj":p_adj,
            "hr_highstage":hr_h,"cox_p_highstage":p_h,"n_highstage":n_h,"ev_highstage":ev_h,
            "logrank_top30_p":lr_p
        })

df_res = pd.DataFrame(rows)
df_res.to_csv(f"{OUT}/tso500_survival_compare_v2.tsv", sep="\t", index=False)
with open(f"{OUT}/tso500_survival_compare_v2.json","w") as f:
    json.dump(rows, f, indent=2)

print("\n=== (A) Continuous Cox + (B) top-30%-derived log-rank ===")
print(f"{'EP':<4s} {'Panel':<18s} {'n':>4s}/{'ev':<3s}  {'HR':>5s} [{'CI':<13s}] {'CoxP':>9s} {'StgAdjP':>9s} {'HiStgHR':>7s} {'HiStgP':>9s} {'LR-top30':>9s}")
for r in rows:
    print(f"  {r['endpoint']:<3s} {r['panel']:<17s}  {r['n_total']:4d}/{r['n_events']:<3d}  {r['hr_per_sd']:5.2f} [{r['hr_ci_lo']:.2f}-{r['hr_ci_hi']:.2f}]  {r['cox_p']:9.3e}  {r['cox_p_stageadj']:9.3e}  {r['hr_highstage']:7.2f}  {r['cox_p_highstage']:9.3e}  {r['logrank_top30_p']:9.3e}")

# ─── KM by truth-DM1 + top30%-derived (Full-8 vs TSO500-3 vs Compact-5) ───
fig, axes = plt.subplots(3, 4, figsize=(20, 13), dpi=160)
for i, (ep, (tcol,ecol,lab)) in enumerate(ENDPOINTS.items()):
    for j, (p_name, (col, color, ptype)) in enumerate(PANELS.items()):
        ax = axes[i,j]
        if ptype=="score":
            label = derive_dm1_label(df[col], frac=0.30)
        else:
            label = df[col].astype(int)
        m = pd.DataFrame({"l":label,"t":df[tcol],"e":df[ecol]}).dropna()
        m = m[m["t"]>0]
        hi = m[m["l"]==1]; lo = m[m["l"]==0]
        if len(hi)>1 and len(lo)>1:
            kmf_hi = KaplanMeierFitter().fit(hi["t"]/365.25, hi["e"], label=f"DM1-like (n={len(hi)})")
            kmf_lo = KaplanMeierFitter().fit(lo["t"]/365.25, lo["e"], label=f"non-DM1  (n={len(lo)})")
            kmf_hi.plot_survival_function(ax=ax, color=color, ci_alpha=0.12, lw=2.2)
            kmf_lo.plot_survival_function(ax=ax, color="#94A3B8", ci_alpha=0.12, lw=2.2)
        r = next(x for x in rows if x["endpoint"]==ep and x["panel"]==p_name)
        ax.text(0.04, 0.08,
                f"top-30% log-rank p = {r['logrank_top30_p']:.2e}\n"
                f"Cox HR (cont.)   = {r['hr_per_sd']:.2f} [{r['hr_ci_lo']:.2f}-{r['hr_ci_hi']:.2f}]\n"
                f"Cox p            = {r['cox_p']:.2e}",
                transform=ax.transAxes, fontsize=8, family="monospace",
                bbox=dict(facecolor="white", edgecolor="#CBD5E1", boxstyle="round,pad=0.4"))
        ax.set_xlabel("Years"); ax.set_ylabel("Survival probability" if j==0 else "")
        ax.set_ylim(0.55, 1.02); ax.set_xlim(0, 15)
        ax.set_title(f"{lab}\n{p_name}", fontsize=11, fontweight="bold", color=color, loc="left")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.legend(loc="lower left", fontsize=8, frameon=False)
fig.suptitle("DM1-derived label (top 30 % by score) — TCGA-THCA KM",
             fontsize=14, fontweight="bold", y=1.005)
plt.tight_layout()
fig.savefig(f"{OUT}/fig_tso500_survival_km_v2.png", bbox_inches="tight")
plt.close(fig)

# ─── Forest ───
fig, ax = plt.subplots(figsize=(12, 7), dpi=170)
ys = []; labels = []; hrs = []; los = []; his = []; ps = []; colors = []
y = 0
panel_order = ["DM1_truth_label","Full-8_score","Compact-5_score","TSO500-3_score"]
for ep in ENDPOINTS:
    for p_name in panel_order:
        r = next(x for x in rows if x["endpoint"]==ep and x["panel"]==p_name)
        color = PANELS[p_name][1]
        ys.append(y); labels.append(f"{ep:<3s}  ·  {p_name}")
        hrs.append(r["hr_per_sd"]); los.append(r["hr_ci_lo"]); his.append(r["hr_ci_hi"])
        ps.append(r["cox_p"]); colors.append(color)
        y -= 1
    y -= 0.5
ys = np.array(ys); hrs = np.array(hrs); los = np.array(los); his = np.array(his)
for i in range(len(ys)):
    if not np.isnan(hrs[i]):
        ax.plot([los[i], his[i]], [ys[i],ys[i]], color="#94A3B8", lw=1.5, zorder=1)
        ax.plot([hrs[i]], [ys[i]], "o", markersize=11, markerfacecolor=colors[i],
                markeredgecolor="white", markeredgewidth=1.5, zorder=3)
ax.axvline(1.0, color="#0F172A", ls="--", lw=1, alpha=0.5)
ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=10, family="monospace")
ax.set_xlabel("Hazard Ratio per +1 SD score  (95 % CI)", fontsize=11)
ax.set_xscale("log"); ax.set_xlim(0.3, 6)
for i, (h, lo, hi, p) in enumerate(zip(hrs, los, his, ps)):
    if not np.isnan(h):
        sig = "★★★" if p<0.001 else "★★" if p<0.01 else "★" if p<0.05 else ""
        txt = f"HR={h:.2f} [{lo:.2f}-{hi:.2f}] p={p:.1e} {sig}"
        ax.text(7, ys[i], txt, va="center", fontsize=8.5, family="monospace")
ax.set_title("Figure SR-2 · Continuous-score Cox HR forest (TCGA-THCA n=478)",
             fontsize=12.5, fontweight="bold", loc="left", pad=12)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.set_xlim(0.3, 30)
plt.tight_layout()
fig.savefig(f"{OUT}/fig_tso500_survival_forest_v2.png", bbox_inches="tight")
plt.close(fig)

print(f"\nOutputs:")
for f in sorted(os.listdir(OUT)):
    if "tso500_survival" in f and "_v2" in f: print(f"  {f}")
