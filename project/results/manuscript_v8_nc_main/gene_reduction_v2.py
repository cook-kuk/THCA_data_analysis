"""
Gene panel reduction analysis — v2 with story-grouped biomarker combinations.

Groups (4 stories):
  A. Clinical assay tier   — Full-8 / TSO500-3 / IHC-3 / IHC-4 / RT-qPCR-4 / Compact-5
  B. Biology axis          — TF-only / RAI-machinery / Effector-only
  C. Literature convergence — Landa-5-overlap
  D. Statistical parsimony  — Compact-6 / Compact-5 / Compact-4 / Top-5-by-d / Top-4-by-d

Metrics per panel:
  · |Cohen's d| (DM1 vs DM2 on HM450 β)
  · 5-fold CV AUC (LogReg on subset β)
  · Spearman r vs Full-8 mean
  · TSO500-native gene count (of 3: TSHR, PAX8, NKX2-1)
  · IHC-routine gene count (of 3: TG, PAX8, NKX2-1)
  · BRAF+ subset PFI Cox HR + p (clinical KM proxy for RAI-refractory)

Figures:
  · fig_reduction_grouped_bars.png   — grouped horizontal bars by story
  · fig_reduction_perf_matrix.png    — panels × 6 metrics heatmap
"""
import os, json, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

for f in ["/usr/share/fonts/truetype/nanum/NanumGothic.ttf"]:
    if os.path.exists(f): font_manager.fontManager.addfont(f); break
plt.rcParams["font.family"] = ["NanumGothic","DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = "/home/seungho/personal/THCA_data_analysis/project"
OUT  = f"{ROOT}/dm1_story_web/public/figures"

GENES = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]
TSO500  = {"TSHR","PAX8","NKX2-1"}
IHC_ROU = {"TG","PAX8","NKX2-1"}
GENE_D  = {"TPO":2.30, "DIO1":1.24, "TSHR":1.20, "PAX8":0.97, "TG":0.86,
           "FOXE1":0.84, "NKX2-1":0.63, "SLC5A5":0.22}

# ─── Load HM450 + clinical (PFI) ───
hm  = pd.read_csv(f"{ROOT}/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
mst = pd.read_csv(f"{ROOT}/results/manuscript_v8_nc_main/master_tcga.tsv", sep="\t")[["sample_short","DM"]]
cli = pd.read_csv(f"{ROOT}/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
df  = hm.merge(mst, on="sample_short", how="inner").merge(
       cli, left_on="sample_short", right_on="tcga_short", how="left").dropna(subset=GENES)
df["DM1"] = (df["DM"]=="DM1").astype(int)
df["braf_pos"] = df["has_braf_v600e"].astype(str).str.lower().isin(["1","true","yes"])
df["score_full"] = df[GENES].mean(axis=1)

# Sign flip: DM1 label direction in master file — align so higher score = DM1-like
r = np.corrcoef(df["score_full"], df["DM1"])[0,1]
sign = -1 if r < 0 else 1
print(f"Cohort n={len(df)}, DM1={df['DM1'].sum()}, sign flip = {sign}")

# ─── PANEL DEFINITIONS (story-grouped) ───
ranked_by_d = sorted(GENES, key=lambda g: -GENE_D[g])

PANELS = [
    # ─ A. Clinical assay tiers ─
    ("Full 8-gene",                    "A", "★ 기준",            GENES),
    ("TSO500-3 (DNA panel native)",    "A", "표준 NGS",          ["TSHR","PAX8","NKX2-1"]),
    ("IHC routine 3 (병리)",           "A", "IHC 3-plex",       ["TG","PAX8","NKX2-1"]),
    ("IHC 4 (+TPO)",                   "A", "IHC 4-plex",       ["TG","PAX8","NKX2-1","TPO"]),
    ("RT-qPCR 4 (minimal)",            "A", "RT-qPCR",          ["PAX8","NKX2-1","TPO","DIO1"]),
    ("Compact 5 (TSO500 + TPO + DIO1)","A", "★ dual assay",      ["PAX8","NKX2-1","TSHR","TPO","DIO1"]),
    # ─ B. Biology axis subsets ─
    ("TF axis only 3",                 "B", "PAX8+NKX2-1+FOXE1", ["PAX8","NKX2-1","FOXE1"]),
    ("RAI machinery 5",                "B", "TSHR+NIS+TPO+TG+DIO1", ["TSHR","SLC5A5","TPO","TG","DIO1"]),
    ("Effector only 4",                "B", "TG+TPO+NIS+DIO1",    ["TG","TPO","SLC5A5","DIO1"]),
    # ─ C. Literature convergence ─
    ("Landa 5-overlap",                "C", "Landa 2016 silenced ∩ 8", ["TG","TSHR","TPO","DIO1","SLC5A5"]),
    # ─ D. Statistical parsimony ─
    ("Compact 6 (drop FOXE1+NIS)",     "D", "안전 축소",         ["PAX8","NKX2-1","TG","TPO","TSHR","DIO1"]),
    ("Compact 4 (TF-2 + effector-2)",  "D", "최소 viable",       ["PAX8","NKX2-1","TPO","DIO1"]),
    ("Top-5 by per-gene d",            "D", "effector-rich",     ranked_by_d[:5]),
    ("Top-4 by per-gene d",            "D", "aggressive top",    ranked_by_d[:4]),
]

# ─── Compute all metrics ───
def compute(sub):
    score = sign * df[list(sub)].mean(axis=1)
    g1 = score[df["DM1"]==1].values
    g0 = score[df["DM1"]==0].values
    sd_p = np.sqrt((g1.var(ddof=1) + g0.var(ddof=1)) / 2)
    d = abs((g1.mean() - g0.mean()) / sd_p) if sd_p > 0 else 0
    X = StandardScaler().fit_transform(df[list(sub)].values)
    auc = cross_val_score(LogisticRegression(max_iter=400), X, df["DM1"].values, cv=5, scoring="roc_auc").mean()
    r = spearmanr(score, sign * df["score_full"]).statistic
    tso = sum(1 for g in sub if g in TSO500)
    ihc = sum(1 for g in sub if g in IHC_ROU)
    # BRAF+ PFI tertile split log-rank + Cox HR per +1 SD
    bs = df[df["braf_pos"]].dropna(subset=["PFI.time","PFI"])
    bs = bs[bs["PFI.time"] > 0]
    if len(bs) < 20 or bs["PFI"].sum() < 4:
        lrp = coxp = hr = np.nan
    else:
        s = sign * bs[list(sub)].mean(axis=1)
        z = (s - s.mean()) / s.std()
        q1, q2 = s.quantile([1/3, 2/3])
        hi = bs[s >= q2]; lo = bs[s <= q1]
        if len(hi) < 5 or len(lo) < 5:
            lrp = np.nan
        else:
            lrp = float(logrank_test(hi["PFI.time"], lo["PFI.time"], hi["PFI"], lo["PFI"]).p_value)
        try:
            mz = bs[["PFI.time","PFI"]].copy(); mz["z"] = z.values
            cph = CoxPHFitter().fit(mz, duration_col="PFI.time", event_col="PFI")
            hr = float(np.exp(cph.params_["z"]))
            coxp = float(cph.summary.loc["z","p"])
        except Exception:
            hr = coxp = np.nan
    return dict(d=d, auc=auc, r=r, tso=tso, ihc=ihc, n=len(sub),
                lrp=lrp, hr=hr, coxp=coxp)

rows = []
for name, group, sub, sub_genes in PANELS:
    m = compute(sub_genes)
    rows.append({"panel":name, "group":group, "subtitle":sub, "genes":" · ".join(sub_genes), **m})

pd.DataFrame(rows).to_csv(f"{OUT}/reduction_v2_summary.tsv", sep="\t", index=False)
with open(f"{OUT}/reduction_v2_summary.json","w") as f:
    json.dump(rows, f, indent=2)

print(f"\n{'Panel':<38s} {'Grp':<3s} {'n':>2s} {'|d|':>5s} {'AUC':>6s} {'r':>6s} {'TSO':>4s} {'IHC':>4s} {'BRAF+PFI-p':>10s} {'HR':>6s}")
for r in rows:
    print(f"  {r['panel']:<36s} {r['group']:<3s} {r['n']:2d} {r['d']:5.2f} {r['auc']:6.3f} {r['r']:6.3f} {r['tso']}/{r['n']:<2d} {r['ihc']}/{r['n']:<2d} {r['lrp']:10.3e} {r['hr']:6.2f}")

# ═════════════════════════════════════════════════════════════════
# FIG 1 — 2-panel layout: horizontal bars (left) + metrics table (right), no overlap
# ═════════════════════════════════════════════════════════════════
GRP_INFO = {
    "A": ("A. 임상 assay 계층",       "#B91C1C"),
    "B": ("B. 생물학 축 (biology)",    "#0E7490"),
    "C": ("C. 문헌 수렴 (Landa 2016)", "#B45309"),
    "D": ("D. 통계 parsimony",         "#047857")
}

fig = plt.figure(figsize=(22, 15), dpi=170)
gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1], wspace=0.02, top=0.94, bottom=0.06, left=0.14, right=0.98)
ax = fig.add_subplot(gs[0, 0])
ax_tbl = fig.add_subplot(gs[0, 1]); ax_tbl.axis("off")

# Sort within group by |d| descending
rows_sorted = []
for gk in ["A","B","C","D"]:
    grp = [r for r in rows if r["group"]==gk]
    grp.sort(key=lambda x: -x["d"])
    rows_sorted.extend(grp)

y_positions = []
labels = []
values_d = []
values_auc = []
colors = []
group_bounds = []
y = 0
last_g = None
BAR_STEP = 1.4    # space per bar
GRP_GAP  = 1.6    # gap between groups
for r in rows_sorted:
    if r["group"] != last_g:
        if last_g is not None:
            y -= GRP_GAP
            group_bounds.append((last_g, gs, y + GRP_GAP*0.5))
        gs = y
        last_g = r["group"]
    y_positions.append(y)
    labels.append(f"{r['panel']}  ({r['n']}g)")
    values_d.append(r["d"])
    values_auc.append(r["auc"])
    colors.append(GRP_INFO[r["group"]][1])
    y -= BAR_STEP
group_bounds.append((last_g, gs, y + BAR_STEP*0.5))

y_positions = np.array(y_positions)
values_d = np.array(values_d)
values_auc = np.array(values_auc)

# Group background bands (bar chart panel)
for gk, ys, ye in group_bounds:
    ax.axhspan(ye - BAR_STEP*0.5, ys + BAR_STEP*0.5, facecolor=GRP_INFO[gk][1], alpha=0.06)
    ax.text(-0.30, (ys + ye)/2, GRP_INFO[gk][0], rotation=90, va="center", ha="center",
            fontsize=13, fontweight="bold", color=GRP_INFO[gk][1])

# Bars — clean, only |d| value at end (no other overlapping text)
ax.barh(y_positions, values_d, color=colors, edgecolor="white", height=BAR_STEP*0.72, alpha=0.94, linewidth=1.5)
for yi, r, d in zip(y_positions, rows_sorted, values_d):
    ax.text(d + 0.03, yi, f"|d| = {d:.2f}", va="center", fontsize=12, family="monospace",
            color="#0F172A", fontweight="bold")

ax.set_yticks(y_positions)
ax.set_yticklabels(labels, fontsize=12)
ax.set_xlabel("|Cohen's d|   (DM1 vs DM2  ·  TCGA HM450 promoter β  ·  n = 478)",
              fontsize=13, labelpad=12)
ax.set_xlim(-0.6, max(values_d) * 1.35)
ax.axvline(0, color="#64748B", lw=1)
ax.axvline(1.5, color="#94A3B8", ls="--", lw=1, alpha=0.5)
ax.text(1.5, y_positions.min() - BAR_STEP*1.4, "|d| = 1.5\n임상 유용 최소",
        fontsize=10, color="#64748B", style="italic", ha="center")
ax.set_title("Figure GR-1 · Biomarker 조합 성능 비교 — 4 스토리 그룹",
             fontsize=17, fontweight="bold", loc="left", pad=18, color="#0F172A")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)

# ── Right panel: dedicated metrics table (no overlap) ──
ax_tbl.set_xlim(0, 100); ax_tbl.set_ylim(y_positions.min() - BAR_STEP*0.8, y_positions.max() + BAR_STEP*1.4)
# Header row
header_y = y_positions.max() + BAR_STEP*0.75
ax_tbl.text( 2, header_y, "n", fontsize=11.5, fontweight="bold", color="#0F172A")
ax_tbl.text( 8, header_y, "AUC (5-CV)", fontsize=11.5, fontweight="bold", color="#0F172A")
ax_tbl.text(28, header_y, "r vs Full-8", fontsize=11.5, fontweight="bold", color="#0F172A")
ax_tbl.text(48, header_y, "TSO500", fontsize=11.5, fontweight="bold", color="#0F172A")
ax_tbl.text(62, header_y, "IHC-routine", fontsize=11.5, fontweight="bold", color="#0F172A")
ax_tbl.text(78, header_y, "BRAF+ PFI p", fontsize=11.5, fontweight="bold", color="#0F172A")
# Header underline
ax_tbl.plot([0, 100], [header_y - BAR_STEP*0.25, header_y - BAR_STEP*0.25],
            color="#0F172A", lw=1.5, clip_on=False)
# Table rows
for yi, r in zip(y_positions, rows_sorted):
    ax_tbl.text( 2, yi, f"{r['n']}",                fontsize=11.5, family="monospace", color="#0F172A", va="center", fontweight="bold")
    ax_tbl.text( 8, yi, f"{r['auc']:.3f}",          fontsize=11.5, family="monospace", color="#047857", va="center", fontweight="bold")
    ax_tbl.text(28, yi, f"{r['r']:+.2f}",           fontsize=11.5, family="monospace", color="#0F172A", va="center")
    ax_tbl.text(48, yi, f"{r['tso']} / {r['n']}",   fontsize=11.5, family="monospace", color="#B91C1C" if r['tso']==r['n'] else "#0F172A", va="center")
    ax_tbl.text(62, yi, f"{r['ihc']} / {r['n']}",   fontsize=11.5, family="monospace", color="#B45309" if r['ihc']==r['n'] else "#0F172A", va="center")
    if not np.isnan(r["lrp"]):
        sig = "★★★" if r["lrp"]<0.001 else "★★" if r["lrp"]<0.01 else "★" if r["lrp"]<0.05 else "n.s."
        color = "#047857" if r["lrp"]<0.01 else "#B45309" if r["lrp"]<0.05 else "#64748B"
        ax_tbl.text(78, yi, f"{r['lrp']:.2e}  {sig}", fontsize=11.5, family="monospace", color=color, va="center", fontweight="bold")
    else:
        ax_tbl.text(78, yi, "n/a", fontsize=11.5, family="monospace", color="#94A3B8", va="center")
# Group separator lines in table
for gk, ys, ye in group_bounds[:-1]:
    sep_y = ye - BAR_STEP*0.5
    ax_tbl.plot([0, 100], [sep_y, sep_y], color="#E2E8F0", lw=1)
plt.tight_layout()
fig.savefig(f"{OUT}/fig_reduction_grouped_bars.png", bbox_inches="tight")
plt.close(fig)

# ═════════════════════════════════════════════════════════════════
# FIG 2 — Performance matrix heatmap
# ═════════════════════════════════════════════════════════════════
metrics_show = ["d", "auc", "r", "tso_pct", "ihc_pct"]
metric_labels = ["|Cohen's d|", "5-fold CV AUC", "Spearman r vs Full-8", "TSO500 native %", "IHC routine %"]

M = np.zeros((len(rows_sorted), len(metrics_show)))
for i, r in enumerate(rows_sorted):
    M[i,0] = r["d"] / 2.5           # normalize |d| to ~0-1 (2.5 = practical max)
    M[i,1] = (r["auc"] - 0.5) / 0.5 # AUC in [0.5, 1.0] → [0, 1]
    M[i,2] = r["r"]                 # already 0-1
    M[i,3] = r["tso"] / r["n"]
    M[i,4] = r["ihc"] / r["n"]

cmap = LinearSegmentedColormap.from_list("perf", ["#FEF2F2","#FED7AA","#86EFAC","#047857"])

fig, ax = plt.subplots(figsize=(14, 10), dpi=170)
im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")

# Cell annotations with actual values
for i, r in enumerate(rows_sorted):
    vals = [
        f"{r['d']:.2f}",
        f"{r['auc']:.3f}",
        f"{r['r']:.2f}",
        f"{r['tso']} / {r['n']}",
        f"{r['ihc']} / {r['n']}"
    ]
    for j, v in enumerate(vals):
        cell_val = M[i,j]
        txt_color = "#FFFFFF" if cell_val > 0.6 else "#0F172A"
        ax.text(j, i, v, ha="center", va="center", fontsize=11, fontweight="bold", color=txt_color, family="monospace")

# Y labels — panel names with group prefix
y_labels = [f"[{r['group']}]  {r['panel']}" for r in rows_sorted]
ax.set_yticks(range(len(rows_sorted)))
ax.set_yticklabels(y_labels, fontsize=10.5)
# Highlight full-8 and compact-5
for i, r in enumerate(rows_sorted):
    if r["panel"] in ("Full 8-gene", "Compact 5 (TSO500 + TPO + DIO1)"):
        ax.get_yticklabels()[i].set_color("#B91C1C")
        ax.get_yticklabels()[i].set_fontweight("bold")

ax.set_xticks(range(len(metric_labels)))
ax.set_xticklabels(metric_labels, fontsize=11, fontweight="bold")
ax.xaxis.tick_top()
ax.tick_params(axis="both", length=0)

# Group separator lines
prev_g = None
for i, r in enumerate(rows_sorted):
    if prev_g and r["group"] != prev_g:
        ax.axhline(i - 0.5, color="#0F172A", lw=1.5)
    prev_g = r["group"]

# Colorbar
cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.02, shrink=0.85)
cbar.set_label("정규화 성능 (0 = 실패 / 1 = 최적)", fontsize=10)
cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

ax.set_title("Figure GR-2 · Biomarker 조합 × 성능 지표 매트릭스",
             fontsize=15, fontweight="bold", loc="left", pad=18, color="#0F172A")
plt.tight_layout()
fig.savefig(f"{OUT}/fig_reduction_perf_matrix.png", bbox_inches="tight")
plt.close(fig)

print(f"\nOutputs:")
for f in sorted(os.listdir(OUT)):
    if "reduction_v2" in f or "reduction_grouped" in f or "reduction_perf" in f: print(f"  {f}")
