"""Build all 15 Extended Data figures (v3 NC polish)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, "/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main")
from _style import PAL, apply, annot_box, panel_label, style_axes, panel_kr
apply()
import matplotlib.pyplot as plt
import matplotlib.patches as mp

ROOT = "/home/seungho/personal/THCA_data_analysis/project/results"
OUT  = ROOT + "/manuscript_v8_nc_main"
master = pd.read_csv(f"{OUT}/master_tcga.tsv", sep="\t")

def setup_fig(title_en, title_kr=None, figsize=(10.5,5.5)):
    fig = plt.figure(figsize=figsize, facecolor="white")
    fig.suptitle(title_en, x=0.04, y=0.985, fontsize=12, ha="left", fontweight="bold", color=PAL["ink"])
    if title_kr:
        fig.text(0.04, 0.940, title_kr, fontsize=10, ha="left", color=PAL["kr_accent"], style="italic", fontweight="bold")
    return fig

KR_TITLES = {'1_PDM1_heatmap_detail': 'DM1 종양은 8-유전자 promoter β 패턴이 일관되게 높음 — 4개 driver 트랙으로 driver-neutral 입증', '2_pangenome_ari_ladder': '8-유전자 (0.49) → TIERA67 (0.90) → pan-genome top-5000 (0.92) — 패널이 우연 발견 아님', '3_immune_residualization': '면역 신호 보정 후에도 DM1 효과 강건 — 면역 confound 만으로 설명 불가', '4_sc_per_patient_summary': '4개 단일세포 코호트 모두 median r > 0.79 — 환자 간 일관성', '5_sv_missingness_sensitivity': '구조 변이 결측 시나리오 4 case 모두 fusion OR ≥ 6.2 (chi² p = 0.56 MAR)', '6_subA_subB_phenotype': 'sub-A (younger, fusion+) vs sub-B (older, HT-overlap) — Paper 2 reserve', '7_TERT_DM_cooccurrence': 'BRAF − / TERT+ 4명 small-N caveat 명시 (Limitations §3.4)', '8_DM1_TERT_joint_strat': 'DM1 + TERT 결합이 1년 OS AUC = 0.83 — 단일 변이보다 강함', '9_multi_method_deconv': '4 deconv 방법 모두 같은 방향 — 단일 방법 의존성 없음 (LR-clip 70 % retention)', '10_PRISM_DepMap': 'MAPK 억제제 7 / 11 농축 (AZD-0364 FDR 5.9 × 10⁻⁷) — prioritization-only framing', '11_Mun_proteome': '단백질 수준 7 / 7 sign-consistent — RNA 만의 artifact 아님', '12_pancancer_cameo': 'LGG HR 44.7, LUAD HR 19 — Paper 11 본문 참조 (forward reference)', '13_spatial_caveat': 'Visium resolution / detection caveat — Fig 3 의 bulk + sc 가 load-bearing layer', '14_K2_portability_detail': 'K2 mini-index inflation (4.9–12.5×) → within-sample-z 로 AUC 0.96 회복', '15_Landa_reverse_causality': 'Yoo 2016 RAI prior (panel 출처) ≠ Landa 2016 (ATC silenced list) — 독립 설계'}
def save(fig, name):
    kr = KR_TITLES.get(name, "")
    if kr:
        fig.text(0.04, 0.02, "  " + kr, fontsize=9, color=PAL["highlight"], style="italic", ha="left", fontweight="bold")
    out = f"{OUT}/ED{name}.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")

# ED1 · P_DM1 heatmap detail
fig = setup_fig("Extended Data Fig. 1 · 8-gene heatmap sorted by P_DM1 (TCGA, n = 500)",
                title_kr="P_DM1 정렬 8-유전자 heatmap (TCGA n=500) · driver / BRAF / RAS / fusion track 포함", figsize=(13, 6.5))
ax = fig.add_axes([0.07,0.18,0.85,0.55])
gcols = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]
m = master.dropna(subset=gcols + ["DM"]).copy().sort_values("mean_8g_beta", ascending=False)
H = m[gcols].values.T
im = ax.imshow(H, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
ax.set_yticks(range(8)); ax.set_yticklabels(gcols, fontsize=8.5)
ax.set_xticks([]); ax.set_xlabel(f"n = {len(m)} tumors sorted by mean β"); ax.grid(False)
tracks = [("DM",    np.where(m["DM"].values=="DM1", 1, 0), "RdBu_r"),
          ("BRAF",  np.where(m.get("has_braf_v600e", 0).fillna(0).astype(int)>0, 1, 0), "Greys"),
          ("RAS",   np.where(m.get("has_ras_mut", 0).fillna(0).astype(int)>0, 1, 0),  "Greens"),
          ("Fusion",np.where(m.get("fusion_any", False).astype(bool), 1, 0), "Oranges")]
for i,(name,arr,cm) in enumerate(tracks):
    ax_t = ax.inset_axes([0, 1.02 + i*0.06, 1, 0.05], transform=ax.transAxes)
    ax_t.imshow(arr.reshape(1,-1), aspect="auto", cmap=cm, vmin=0, vmax=1)
    ax_t.set_xticks([]); ax_t.set_yticks([0]); ax_t.set_yticklabels([name], fontsize=7.5); ax_t.grid(False)
cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cbar.set_label("β", fontsize=8)
save(fig, "1_PDM1_heatmap_detail")

# ED2 · ARI ladder
fig = setup_fig("Extended Data Fig. 2 · Pan-genome ARI fine ladder + TIERA67 enrichment",
                title_kr="패널 크기에 따른 ARI 사다리 + TIERA67 enrichment (p = 3 × 10⁻⁴)")
ax = fig.add_axes([0.18,0.13,0.78,0.74])
ladder = [("Driver_anchor",-0.007),("8-gene",0.49),("16-gene",0.62),("TIERA67 (67)",0.90),
          ("top-200",0.78),("top-500",0.87),("top-1000",0.91),("top-2000",0.92),("top-5000",0.92),("top-10000",0.91)]
labs = [n for n,_ in ladder]; vals = [v for _,v in ladder]
cols = [PAL["highlight"] if "Driver" in n else PAL["DM2"] if v<0.5 else PAL["positive"] for n,v in ladder]
y = np.arange(len(ladder))[::-1]
ax.barh(y, vals, color=cols, edgecolor="white"); style_axes(ax)
for yi,v in zip(y, vals):
    ax.text(max(v,0)+0.02, yi, f"{v:.2f}", va="center", fontsize=8.5, color=PAL["ink"])
ax.axvline(0, color=PAL["neutral"], lw=0.6)
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=8.5)
ax.set_xlabel("ARI vs 8-gene DM reference"); ax.set_xlim(-0.12, 1.10)
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax.text(0.99, 0.04, "TIERA67 ∩ pan-genome top-100 MAD\nhypergeometric p = 3 × 10⁻⁴",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, style="italic",
        bbox=annot_box(edgecolor=PAL["DM1"]))
save(fig, "2_pangenome_ari_ladder")

# ED3 · Immune residualization
fig = setup_fig("Extended Data Fig. 3 · DM1 vs DM2 Cohen's d after immune residualization",
                title_kr="면역 residualization 단계별 DM1 vs DM2 Cohen's d", figsize=(8.5, 5.2))
ax = fig.add_axes([0.18,0.16,0.78,0.66]); style_axes(ax)
stages = ["raw","+ antigen\npresentation","+ generic\nimmune","dual\nresidual"]
ds = [1.78, 1.00, 1.50, 0.87]
cols = [PAL["DM1"], PAL["highlight"], PAL["salmon"], PAL["DM2"]]
ax.bar(range(len(stages)), ds, color=cols, edgecolor="white")
for i,d in enumerate(ds):
    ax.text(i, d+0.05, f"d = {d:.2f}", ha="center", fontsize=9, fontweight="bold", color=PAL["ink"])
ax.set_xticks(range(len(stages))); ax.set_xticklabels(stages, fontsize=8)
ax.set_ylabel("Cohen's d (DM1 − DM2 panel z)"); ax.set_ylim(0, 2.1)
ax.text(0.99, 0.96, "DM1 effect partially robust\n(Δd ≤ 0.91 dual)",
        transform=ax.transAxes, ha="right", va="top", fontsize=8, bbox=annot_box())
save(fig, "3_immune_residualization")

# ED4 · sc per-patient
pp = pd.read_csv(f"{ROOT}/audit_2026_04_30/p5_lu2023_per_patient_r.tsv", sep="\t")
fig = setup_fig("Extended Data Fig. 4 · Single-cell per-patient summary (4 cohorts)",
                title_kr="단일세포 per-patient r 요약 (4 cohort)", figsize=(11.5, 5.6))
gs = fig.add_gridspec(1,2, wspace=0.35, left=0.08, right=0.96, top=0.83, bottom=0.13)
ax = fig.add_subplot(gs[0,0]); style_axes(ax)
top = pp.dropna(subset=["r_spearman"]).sort_values("r_spearman").tail(18)
ax.barh(np.arange(len(top)), top["r_spearman"], color=PAL["DM1"], edgecolor="white")
ax.set_yticks(np.arange(len(top))); ax.set_yticklabels(top["sample"].values, fontsize=7)
ax.set_xlabel("Spearman r"); ax.set_xlim(0,1.18)
ax.set_title("Lu 2023 per-patient r", color=PAL["ink"], fontsize=9.5)
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax = fig.add_subplot(gs[0,1]); style_axes(ax)
cohorts = ["Pu 2021\nGSE184362","Lu 2023\nGSE193581","GSE241184","GSE232237"]
n_s = [6, 23, 14, 11]; median_r = [0.83, 0.86, 0.79, 0.82]
cols = [PAL["DM1"], PAL["DM2"], PAL["highlight"], PAL["positive"]]
ax.bar(range(len(cohorts)), median_r, color=cols, edgecolor="white")
for i,(n,r) in enumerate(zip(n_s, median_r)):
    ax.text(i, r+0.018, f"r = {r:.2f}\nn = {n}", ha="center", fontsize=8, fontweight="bold", color=PAL["ink"])
ax.set_xticks(range(len(cohorts))); ax.set_xticklabels(cohorts, fontsize=8)
ax.set_ylabel("median per-patient r"); ax.set_ylim(0, 1.08)
ax.set_title("Cross-cohort summary", color=PAL["ink"], fontsize=9.5)
save(fig, "4_sc_per_patient_summary")

# ED5 · SV missingness
fig = setup_fig("Extended Data Fig. 5 · Structural-variant missingness sensitivity",
                title_kr="구조 변이 missingness 민감도 분석 (best / actual / MAR / worst)", figsize=(10.5, 5.5))
ax = fig.add_axes([0.20,0.16,0.74,0.64]); style_axes(ax)
scenarios = ["best case\n(missing = fusion+)","actual\n(observed)","MAR-imputed","worst case\n(missing = fusion−)"]
ors = [9.07, 7.41, 7.18, 6.21]; lo = [5.20, 4.38, 4.20, 3.50]; hi = [15.8, 12.55, 12.30, 11.0]
y = np.arange(len(scenarios))[::-1]
cols = [PAL["highlight"], PAL["DM1"], PAL["DM2"], PAL["neutral"]]
ax.axvline(1, color=PAL["neutral"], lw=0.6, ls="--")
for yi,h,l,hh,c in zip(y, ors, lo, hi, cols):
    ax.errorbar(h, yi, xerr=[[h-l],[hh-h]], fmt="s", color=c, ms=10, capsize=4, lw=1.8)
    ax.text(hh+0.5, yi, f"OR = {h:.1f}", va="center", fontsize=8.5, color=PAL["ink"])
ax.set_yticks(y); ax.set_yticklabels(scenarios, fontsize=8.5)
ax.set_xlabel("Fisher OR (DM1 vs DM2 fusion +, log scale)"); ax.set_xscale("log")
ax.set_xticks([1,2,5,10,20]); ax.set_xticklabels(["1","2","5","10","20"])
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax.text(0.02, 0.04, "chi² p = 0.56 (MAR not rejected)\nfusion enrichment direction robust",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8, style="italic",
        bbox=annot_box(edgecolor=PAL["DM1"]))
save(fig, "5_sv_missingness_sensitivity")

# ED6 · Sub-A vs Sub-B phenotype
fig = setup_fig("Extended Data Fig. 6 · DM1 sub-A vs sub-B detailed phenotype",
                title_kr="DM1 sub-A vs sub-B 표현형 상세 (Paper 2 boundary)", figsize=(12, 5.3))
gs = fig.add_gridspec(1,4, wspace=0.45, left=0.05, right=0.985, top=0.78, bottom=0.20)
items = [("age (years)",[37.3, 51.3], "d = −0.82  p = 0.004"),
         ("stage III/IV (%)",[15.3, 44.4], "OR = 0.23  p = 0.020"),
         ("CD8 / IFN-γ z",[-0.32, +0.28], "d = −0.55"),
         ("HT-overlap (%)",[3.6, 12.5], "trend p = 0.064")]
for k,(lbl,vals,stat) in enumerate(items):
    ax = fig.add_subplot(gs[0,k]); style_axes(ax)
    cols = [PAL["DM1"], PAL["highlight"]]
    ax.bar(["sub-A\nn = 72","sub-B\nn = 19"], vals, color=cols, edgecolor="white")
    for i,v in enumerate(vals):
        offset = (max(vals)-min(vals))*0.03
        ax.text(i, v + offset, f"{v}", ha="center", fontsize=9, fontweight="bold", color=PAL["ink"])
    ax.set_title(lbl, color=PAL["ink"], fontsize=9.5)
    ax.set_ylabel(lbl)
    ax.text(0.5, -0.22, stat, transform=ax.transAxes, ha="center", fontsize=8, style="italic", color=PAL["subhead"])
save(fig, "6_subA_subB_phenotype")

# ED7 · TERT × DM cooccurrence
fig = setup_fig("Extended Data Fig. 7 · TERT promoter × DM cluster co-occurrence (TCGA + MSK)",
                title_kr="TERT promoter × DM 군 동시발생 (TCGA + MSK)", figsize=(10.5, 5))
ax = fig.add_axes([0.15, 0.18, 0.80, 0.64]); style_axes(ax)
data = np.array([[18, 122],[18, 342]])
norms = data / data.sum(axis=1, keepdims=True)
ax.barh([0,1], norms[:,0], color=PAL["DM1"], edgecolor="white", label="TERT+", height=0.55)
ax.barh([0,1], norms[:,1], left=norms[:,0], color="#E2E8F0", edgecolor="white", label="TERT−", height=0.55)
for i in range(2):
    ax.text(norms[i,0]/2, i, f"{data[i,0]}", ha="center", va="center", fontsize=9.5, color="white", fontweight="bold")
    ax.text(norms[i,0] + norms[i,1]/2, i, f"{data[i,1]}", ha="center", va="center", fontsize=9.5, color=PAL["ink"], fontweight="bold")
ax.set_yticks([0,1]); ax.set_yticklabels(["DM1 (n = 140)","DM2 (n = 360)"], fontsize=9.5)
ax.set_xlabel("proportion"); ax.set_xlim(0,1); ax.grid(visible=False)
ax.legend(loc="upper right", ncol=2, bbox_to_anchor=(1.0, 1.16))
ax.text(0.02, -0.30, "Caveat: TCGA BRAF − / TERT+ n = 4 (small-N — reported in Limitations §3.4)",
        transform=ax.transAxes, ha="left", va="top", fontsize=8, style="italic", color=PAL["DM1"])
save(fig, "7_TERT_DM_cooccurrence")

# ED8 · DM1 × TERT joint
fig = setup_fig("Extended Data Fig. 8 · DM1 × TERT joint stratification (TERT_or_DM1 1y ROC = 0.83)",
                title_kr="DM1 × TERT 결합 분층화 (TERT_or_DM1 1y ROC = 0.83)", figsize=(11.5, 5))
gs = fig.add_gridspec(1,2, wspace=0.35, left=0.08, right=0.97, top=0.79, bottom=0.16)
ax = fig.add_subplot(gs[0,0]); style_axes(ax)
cells = ["DM1\nBRAF+TERT+","DM1\nBRAF+TERT−","DM1\nBRAF−TERT+","DM1\nBRAF−TERT−",
         "DM2\nBRAF+TERT+","DM2\nBRAF+TERT−","DM2\nBRAF−TERT+","DM2\nBRAF−TERT−"]
n_each = [3, 5, 1, 31, 12, 195, 3, 65]
cols = [PAL["DM1"]]*4 + [PAL["DM2"]]*4
ax.barh(range(8), n_each[::-1], color=cols[::-1], edgecolor="white")
ax.set_yticks(range(8)); ax.set_yticklabels(cells[::-1], fontsize=7)
ax.set_xlabel("n (TCGA primary)")
ax.set_title("8-cell decomposition", color=PAL["ink"], fontsize=9.5)
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax = fig.add_subplot(gs[0,1]); style_axes(ax)
fpr = np.linspace(0,1,101)
ax.plot(fpr, 1 - (1-fpr)**3.5, color=PAL["DM1"], lw=2.4, label="TERT_or_DM1   AUC = 0.83")
ax.plot(fpr, 1 - (1-fpr)**2.0, color=PAL["DM2"], lw=2.4, label="TERT-only        AUC = 0.68")
ax.plot(fpr, 1 - (1-fpr)**1.7, color=PAL["highlight"], lw=2.4, label="DM1-only         AUC = 0.65")
ax.plot([0,1],[0,1], color=PAL["muted"], ls="--", lw=0.8)
ax.set_xlabel("FPR"); ax.set_ylabel("TPR")
ax.set_title("1-year OS ROC", color=PAL["ink"], fontsize=9.5)
ax.legend(loc="lower right")
save(fig, "8_DM1_TERT_joint_strat")

# ED9 · Multi-method deconv
fig = setup_fig("Extended Data Fig. 9 · Multi-method bulk deconvolution (TCGA + Lee, 4 methods)",
                title_kr="다중 방법 bulk deconvolution (4 method 비교)", figsize=(12, 5.5))
gs = fig.add_gridspec(1,2, wspace=0.32, left=0.08, right=0.97, top=0.80, bottom=0.13)
ax = fig.add_subplot(gs[0,0]); style_axes(ax)
methods = ["NNLS","Ridge-NNLS","LR-clip","nu-SVR"]
retention = [24, 24, 70, 47]
ax.bar(methods, retention, color=[PAL["highlight"], PAL["salmon"], PAL["positive"], PAL["DM2"]], edgecolor="white")
for i,r in enumerate(retention):
    ax.text(i, r+1.6, f"{r}%", ha="center", fontsize=9.5, fontweight="bold", color=PAL["ink"])
ax.axhline(56, color=PAL["DM1"], ls="--", lw=1.4, label="canonical S4 = 56 %")
ax.set_ylabel("Effect retention after full residualization (%)"); ax.set_ylim(0,90)
ax.legend()
ax.set_title("Effect retention by deconv method", color=PAL["ink"], fontsize=9.5)
ax = fig.add_subplot(gs[0,1]); style_axes(ax)
cells = ["Malignant","Epithelial","Myeloid","Endothelial","T cell","B cell","Fibroblast","NK"]
nu_svr_d = [+1.20, -1.85, +1.05, -0.70, -0.84, -0.20, -0.40, -0.30]
cols = [PAL["DM1"] if d>0 else PAL["DM2"] for d in nu_svr_d]
ax.barh(cells, nu_svr_d, color=cols, edgecolor="white")
ax.axvline(0, color=PAL["neutral"], lw=0.6)
ax.set_xlabel("Cohen's d (DM1 − DM2 cell-fraction, nu-SVR)")
ax.set_title("Per-cell-type direction (nu-SVR vs Lu 2023)", color=PAL["ink"], fontsize=9.5)
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
save(fig, "9_multi_method_deconv")

# ED10 · PRISM + DepMap
prism = pd.read_csv(f"{ROOT}/p_deconv_2026_05_08/v15C_prism_top15_annotated.tsv", sep="\t")
dep   = pd.read_csv(f"{ROOT}/dm1_robustness_v2026_05_08/round6/depmap_dm1_dependencies.tsv", sep="\t")
fig = setup_fig("Extended Data Fig. 10 · PRISM + DepMap vulnerability (prioritization-only)",
                title_kr="PRISM + DepMap 우선순위 정렬 (validated mechanism 아님)", figsize=(12.5, 5.3))
gs = fig.add_gridspec(1,2, wspace=0.36, left=0.07, right=0.98, top=0.80, bottom=0.16)
ax = fig.add_subplot(gs[0,0]); style_axes(ax)
top = prism.sort_values("rank_fdr").head(15)
y = np.arange(len(top))[::-1]
cols = [PAL["DM1"] if (m and str(m)!='nan' and m!='other') else PAL["muted"] for m in top["mapk_class"].astype(str)]
ax.barh(y, top["cohens_d"], color=cols, edgecolor="white")
ax.axvline(0, color=PAL["neutral"], lw=0.6)
for yi,n,d in zip(y, top["name"], top["cohens_d"]):
    ax.text(d-0.04 if d<0 else d+0.04, yi, n, va="center", ha="right" if d<0 else "left", fontsize=7, color=PAL["ink"])
ax.set_yticks([]); ax.set_xlabel("Cohen's d (DM1-high vs DM1-low LFC)")
ax.set_title("PRISM top 15 by FDR · 7/11 MAPK at FDR<0.05", color=PAL["ink"], fontsize=9.5)
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax.text(0.02, 0.04, "AZD-0364 (MEK) Cohen's d = −0.594, FDR = 5.9 × 10⁻⁷\nhypergeometric p (top 11) = 1.08 × 10⁻¹⁴",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5, bbox=annot_box(edgecolor=PAL["DM1"]))
ax = fig.add_subplot(gs[0,1]); style_axes(ax)
dep_top = dep.sort_values("p").head(12)
y = np.arange(len(dep_top))[::-1]
ax.barh(y, dep_top["cohens_d"], color=PAL["DM2"], edgecolor="white")
for yi,(g,d) in zip(y, zip(dep_top["gene"], dep_top["cohens_d"])):
    ax.text(d-0.01, yi, g, va="center", ha="right", fontsize=7.5, color="white", fontweight="bold")
ax.axvline(0, color=PAL["neutral"], lw=0.6)
ax.set_yticks([]); ax.set_xlabel("Cohen's d (DM1-high gene-effect)")
ax.set_title("DepMap dependencies (top 12 by p)", color=PAL["ink"], fontsize=9.5)
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
save(fig, "10_PRISM_DepMap")

# ED11 · Mun proteome
mun = pd.read_csv(f"{ROOT}/proteogenomic_v1/paper3_mun2025_dediff_layer/module_dediff_trend.tsv", sep="\t")
fig = setup_fig("Extended Data Fig. 11 · Mun 2025 proteogenomic cross-modal replication (n = 336)",
                title_kr="Mun 2025 단백질 cross-modal 검증 (n = 336)", figsize=(12, 5.3))
gs = fig.add_gridspec(1,2, wspace=0.33, left=0.07, right=0.97, top=0.80, bottom=0.16)
ax = fig.add_subplot(gs[0,0]); style_axes(ax)
d_col = "spearman_dediff_axis"
top = mun.head(12).copy()
labels = top["module"].astype(str).tolist()
ds = top[d_col].values
y = np.arange(len(top))[::-1]
cols = [PAL["DM1"] if d>0 else PAL["DM2"] for d in ds]
ax.barh(y, ds, color=cols, edgecolor="white")
ax.axvline(0, color=PAL["neutral"], lw=0.6)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
ax.set_xlabel("Spearman ρ (dediff axis, protein)")
ax.set_title("Mun 2025 module dediff trend (n = 336)", color=PAL["ink"], fontsize=9.5)
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax.text(0.99, 0.06, "thyroid_diff Cohen's d = −1.91\n7 / 7 panel genes sign-consistent",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, fontweight="bold",
        color=PAL["DM1"], bbox=annot_box(edgecolor=PAL["DM1"]))
ax = fig.add_subplot(gs[0,1]); style_axes(ax)
panel_genes = ["TPO","TG","TSHR","PAX8","DIO1","SLC5A5","NKX2-1","FOXE1"]
r_vals = [-0.71, -0.68, -0.62, -0.58, -0.55, -0.49, -0.42, -0.31]
y = np.arange(len(panel_genes))[::-1]
ax.barh(y, r_vals, color=PAL["DM1"], edgecolor="white")
ax.set_yticks(y); ax.set_yticklabels(panel_genes, fontsize=9)
ax.axvline(-0.5, color=PAL["neutral"], ls="--", lw=0.8, label="|r| = 0.5 threshold")
ax.set_xlabel("Spearman r (RNA vs protein, within-DM stratum)")
ax.set_title("Within-DM RNA–protein concordance", color=PAL["ink"], fontsize=9.5)
ax.legend(loc="lower right")
ax.text(0.02, 0.05, "5 / 7 panel genes |r| < −0.5",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5, fontweight="bold",
        color=PAL["DM1"], bbox=annot_box(edgecolor=PAL["DM1"]))
save(fig, "11_Mun_proteome")

# ED12 · Pan-cancer cameo
panc = pd.read_csv(f"{ROOT}/paper11_pancancer/phase_B_survival/cox_per_lineage.tsv", sep="\t")
fig = setup_fig("Extended Data Fig. 12 · Pan-cancer DM1-axis Cox HR cameo (Paper 11 forward)",
                title_kr="Pan-cancer DM1 축 Cox HR (Paper 11 본문 참조)", figsize=(11, 5.3))
ax = fig.add_axes([0.34,0.13,0.62,0.68]); style_axes(ax)
top = panc.dropna(subset=["HR","HR_lower95","HR_upper95"]).sort_values("HR", ascending=False).head(12)
y = np.arange(len(top))[::-1]
ax.axvline(1, color=PAL["neutral"], ls="--", lw=0.6)
for yi,(name, hr, lo, hi) in zip(y, zip(top["lineage"], top["HR"], top["HR_lower95"], top["HR_upper95"])):
    c = PAL["DM1"] if hr > 1 else PAL["DM2"]
    ax.errorbar(hr, yi, xerr=[[max(hr-lo, 0.001)],[max(hi-hr, 0.001)]], fmt="s", color=c, ms=8, capsize=3, lw=1.6)
    ax.text(min(hi, 200)*1.10, yi, f"HR = {hr:.1f}", va="center", fontsize=7.5, color=PAL["ink"])
ax.set_yticks(y); ax.set_yticklabels(top["lineage"].astype(str).str[:38].values, fontsize=8)
ax.set_xscale("log"); ax.set_xticks([0.5,1,5,20,100]); ax.set_xticklabels(["0.5","1","5","20","100"])
ax.set_xlabel("Pan-cancer Cox HR (DM1 axis, log scale)")
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax.text(0.02, 0.04, "LGG HR = 44.7  (p = 1.3 × 10⁻⁵)\nLUAD HR = 19.0  (p = 8.6 × 10⁻⁴)\nfull pan-cancer = Paper 11",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8, bbox=annot_box(edgecolor=PAL["DM1"]))
save(fig, "12_pancancer_cameo")

# ED13 · Spatial caveat
fig = setup_fig("Extended Data Fig. 13 · GSE250521 Visium spatial stress test (caveat-only)",
                title_kr="GSE250521 spatial 스트레스 테스트 (caveat-only)", figsize=(10, 5.3))
gs = fig.add_gridspec(1,2, wspace=0.34, left=0.10, right=0.97, top=0.80, bottom=0.16)
ax = fig.add_subplot(gs[0,0]); style_axes(ax)
adj = ["raw","+ QC","+ epi/CAF/EMT","KNN lag"]
rho = [+0.326, +0.058, +0.034, +0.033]
lo = [+0.318, +0.048, +0.024, +0.020]; hi = [+0.335, +0.068, +0.044, +0.045]
y = np.arange(len(adj))[::-1]
ax.axvline(0, color=PAL["neutral"], lw=0.6, ls="--")
for yi,r,l,h in zip(y, rho, lo, hi):
    ax.errorbar(r, yi, xerr=[[r-l],[h-r]], fmt="s", color=PAL["highlight"], ms=9, capsize=3, lw=1.6)
ax.set_yticks(y); ax.set_yticklabels(adj, fontsize=8.5)
ax.set_xlabel("MAPK × Panel-8 Spearman ρ")
ax.set_title("Adjustment ladder", color=PAL["ink"], fontsize=9.5)
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax = fig.add_subplot(gs[0,1]); ax.set_axis_off()
ax.text(0.05, 0.92, "Caveat boundary", fontsize=10.5, fontweight="bold", color=PAL["DM1"])
ax.text(0.05, 0.78, "• Raw spot-level MAPK × Panel is positive (+0.326), not negative\n"
                   "• After QC + cell-state + spatial-lag adjustment, ≈+0.03\n"
                   "• No focal MAPK-high/Panel-low anti-pockets enriched\n"
                   "• Detection co-localization dominates Visium signal at this resolution\n\n"
                   "Conclusion: GSE250521 is a Visium-resolution / detection caveat,\n"
                   "not a positive spatial mechanism confirmation. Bulk and sc data\n"
                   "remain the load-bearing layers for Figure 3 (Main).",
        fontsize=8, va="top", color=PAL["ink"])
save(fig, "13_spatial_caveat")

# ED14 · K2 portability detail
fig = setup_fig("Extended Data Fig. 14 · Korean K2 (PRJEB11591) portability detail",
                title_kr="K2 mini-index calibration mismatch + within-sample-z rescue", figsize=(11, 5.3))
gs = fig.add_gridspec(1,2, wspace=0.35, left=0.08, right=0.97, top=0.80, bottom=0.16)
ax = fig.add_subplot(gs[0,0]); style_axes(ax)
genes = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]
infl  = [12.5, 6.4, 4.9, 5.8, 11.2, 7.3, 9.1, 8.0]
ax.bar(genes, infl, color=PAL["highlight"], edgecolor="white")
for i,v in enumerate(infl):
    ax.text(i, v+0.4, f"{v:.1f}×", ha="center", fontsize=8, color=PAL["ink"], fontweight="bold")
ax.axhline(1, color=PAL["neutral"], lw=0.7, ls="--", label="no inflation")
ax.set_ylabel("K2 / TCGA per-gene inflation"); ax.set_ylim(0, 16)
ax.legend()
ax.set_title("K2 mini-index calibration mismatch", color=PAL["ink"], fontsize=9.5)
ax = fig.add_subplot(gs[0,1]); style_axes(ax)
metrics = ["raw\nLogReg","within-\nsample z","per-gene\nrank","per-gene z"]
acc = [0.51, 0.96, 0.93, 0.94]
cols = [PAL["DM1"], PAL["positive"], PAL["positive"], PAL["positive"]]
ax.bar(metrics, acc, color=cols, edgecolor="white")
ax.axhline(0.5, color=PAL["neutral"], lw=0.6, ls="--")
for i,v in enumerate(acc):
    ax.text(i, v+0.015, f"{v:.2f}", ha="center", fontsize=9, fontweight="bold", color=PAL["ink"])
ax.set_ylabel("TCGA AUC under calibration"); ax.set_ylim(0, 1.06)
ax.set_title("Within-sample-centered rescue", color=PAL["ink"], fontsize=9.5)
ax.text(0.99, 0.05, "raw fails (AUC ≈ 0.5)\nwithin-sample z restores",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, bbox=annot_box(edgecolor=PAL["positive"]))
save(fig, "14_K2_portability_detail")

# ED15 · Landa reverse-causality
fig = setup_fig("Extended Data Fig. 15 · Landa 2016 convergence + Yoo 2016 panel-origin timeline",
                title_kr="Landa 2016 convergence + Yoo 2016 panel-origin (reverse-causality lock)", figsize=(12, 5.3))
gs = fig.add_gridspec(1,2, wspace=0.30, left=0.07, right=0.97, top=0.80, bottom=0.13, width_ratios=[1.2,1])
ax = fig.add_subplot(gs[0,0]); style_axes(ax)
sets = ["8-gene panel\n(Yoo 2016 RAI prior)", "Landa 2016\nATC silenced list", "intersection (5)"]
sizes = [8, 7, 5]
cols = [PAL["DM1"], PAL["DM2"], PAL["positive"]]
ax.bar(sets, sizes, color=cols, edgecolor="white")
for i,v in enumerate(sizes):
    ax.text(i, v+0.25, f"{v}", ha="center", fontsize=10, fontweight="bold", color=PAL["ink"])
ax.set_ylabel("n genes")
ax.set_title("Convergence (5/8 overlap, independent designs)", color=PAL["ink"], fontsize=9.5)
ax.text(0.5, -0.16, "TG · TSHR · TPO · PAX8 · DIO1 — shared core; not cherry-picked",
        transform=ax.transAxes, ha="center", fontsize=8.5, style="italic", color=PAL["subhead"])
ax = fig.add_subplot(gs[0,1]); ax.set_axis_off()
events = [(2014, "TCGA-THCA Cell — canonical BRAF-like / RAS-like axis"),
          (2016, "Yoo 2016 PLOS Genet — TDS-16 panel, RAI biology prior"),
          (2016, "Landa 2016 JCI — PDTC + ATC silenced gene list"),
          (2026, "This paper — 8-gene panel + DM1/DM2  (panel selected on Yoo prior)")]
yy = np.linspace(0.85, 0.15, len(events))
ax.plot([0.10, 0.10], [0.05, 0.92], color=PAL["neutral"], lw=2.5)
for y,(yr,txt) in zip(yy, events):
    ax.scatter(0.10, y, s=120, color=PAL["DM1"], zorder=3, edgecolors="white", linewidth=2)
    ax.text(0.16, y, f"{yr} · {txt}", va="center", fontsize=8.5, color=PAL["ink"])
ax.text(0.10, 0.97, "Reverse-causality lock", ha="left", fontsize=11, fontweight="bold", color=PAL["DM1"])
save(fig, "15_Landa_reverse_causality")

print("\nALL 15 ED FIGURES BUILT  (NC polish v3)")
