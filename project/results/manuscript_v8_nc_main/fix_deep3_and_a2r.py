"""
Fix 2 figures user reported broken:
  - Fig DEEP-3 Panel C empty (GSE151179 has no paired pre/post) → drop Panel C, expand to 2 panels
  - Fig A2-R Panel B broken (Korean text in monospace) → use sans-serif consistent + explicit font family
"""
import os, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import mannwhitneyu

for f in ["/usr/share/fonts/truetype/nanum/NanumGothic.ttf"]:
    if os.path.exists(f): font_manager.fontManager.addfont(f); break
plt.rcParams["font.family"] = ["NanumGothic", "DejaVu Sans"]
plt.rcParams["font.monospace"] = ["DejaVu Sans Mono"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = "/home/seungho/personal/THCA_data_analysis/project"
OUT  = f"{ROOT}/dm1_story_web/public/figures"
GENES = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]

# ═══════════════════════════════════════════════════════════════
# FIX 1 · Fig DEEP-3 · GSE151179 pre-RAI predicts refractory
# Remove empty Panel C, keep 2 clean panels
# ═══════════════════════════════════════════════════════════════
scores_f = f"{ROOT}/results/aggressive_sprint_2026_05_06/gse151179_rai_scores.tsv"
d = pd.read_csv(scores_f, sep="\t", index_col=0)
genes_present = [g for g in GENES if g in d.columns]
d["dm1_expr"] = d[genes_present].mean(axis=1)
d["dm1_axis_neg"] = -d["dm1_expr"]

d["is_refractory"] = pd.to_numeric(d.get("response_refractory"), errors="coerce")
d["no_uptake"]     = pd.to_numeric(d.get("uptake_no"), errors="coerce")
d["pre_rai"]       = pd.to_numeric(d.get("is_pre_rai"), errors="coerce")
d["is_primary"]    = pd.to_numeric(d.get("is_primary"), errors="coerce")

fig, axes = plt.subplots(1, 2, figsize=(16, 7), dpi=170)

# ── Panel A · Refractory vs Responsive (pre-RAI primary tumors)
ax = axes[0]
m1 = d[(d["is_primary"]==1) & (d["pre_rai"]==1)].dropna(subset=["is_refractory","dm1_axis_neg"])
refr = m1[m1["is_refractory"]==1]["dm1_axis_neg"].values
resp = m1[m1["is_refractory"]==0]["dm1_axis_neg"].values
if len(refr) >= 2 and len(resp) >= 2:
    u, p = mannwhitneyu(refr, resp, alternative="greater")
    d_cohen = (refr.mean() - resp.mean()) / np.sqrt((refr.var(ddof=1)+resp.var(ddof=1))/2)
else:
    p, d_cohen = np.nan, np.nan

rng = np.random.default_rng(2026)
x_resp = 0 + rng.uniform(-0.10, 0.10, len(resp))
x_refr = 1 + rng.uniform(-0.10, 0.10, len(refr))
ax.scatter(x_resp, resp, c="#0E7490", s=140, alpha=0.75, edgecolor="white", linewidth=1.5, zorder=3)
ax.scatter(x_refr, refr, c="#B91C1C", s=140, alpha=0.75, edgecolor="white", linewidth=1.5, zorder=3)
if len(resp): ax.hlines(resp.mean(), -0.18, 0.18, color="#0E7490", lw=4)
if len(refr): ax.hlines(refr.mean(), 0.82, 1.18, color="#B91C1C", lw=4)
ax.set_xticks([0, 1])
ax.set_xticklabels([f"Responsive\n(n = {len(resp)})", f"Refractory\n(n = {len(refr)})"], fontsize=12, fontweight="bold")
ax.set_ylabel("DM1 axis score (pre-RAI primary tumor)", fontsize=12)
ax.set_xlim(-0.5, 1.5)
ax.set_title("Panel A · Pre-RAI 조직의 DM1 axis → RAI 결과 예측",
             fontsize=13.5, fontweight="bold", loc="left", color="#0F172A", pad=12)
stat_txt = (f"Cohen's d = {d_cohen:+.2f}\n"
            f"MW p (1-sided) = {p:.3f}\n"
            f"n = {len(resp)} responsive vs {len(refr)} refractory\n"
            f"→ direction 방향 일치 · n 부족")
ax.text(0.03, 0.97, stat_txt, transform=ax.transAxes, fontsize=11.5, va="top",
        fontfamily="NanumGothic",
        bbox=dict(facecolor="#FFFFFF", edgecolor="#CBD5E1", boxstyle="round,pad=0.7", linewidth=1.5))
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.grid(True, ls=":", alpha=0.3, axis="y")

# ── Panel B · Uptake NO vs YES
ax = axes[1]
m2 = d[d["is_primary"]==1].dropna(subset=["no_uptake","dm1_axis_neg"])
no_up  = m2[m2["no_uptake"]==1]["dm1_axis_neg"].values
yes_up = m2[m2["no_uptake"]==0]["dm1_axis_neg"].values
if len(no_up) >= 3 and len(yes_up) >= 3:
    u, p = mannwhitneyu(no_up, yes_up, alternative="greater")
    d_cohen = (no_up.mean() - yes_up.mean()) / np.sqrt((no_up.var(ddof=1)+yes_up.var(ddof=1))/2)
else:
    p, d_cohen = np.nan, np.nan

x_yes = 0 + rng.uniform(-0.10, 0.10, len(yes_up))
x_no  = 1 + rng.uniform(-0.10, 0.10, len(no_up))
ax.scatter(x_yes, yes_up, c="#0E7490", s=140, alpha=0.75, edgecolor="white", linewidth=1.5, zorder=3)
ax.scatter(x_no,  no_up,  c="#B91C1C", s=140, alpha=0.75, edgecolor="white", linewidth=1.5, zorder=3)
if len(yes_up): ax.hlines(yes_up.mean(), -0.18, 0.18, color="#0E7490", lw=4)
if len(no_up):  ax.hlines(no_up.mean(), 0.82, 1.18, color="#B91C1C", lw=4)
ax.set_xticks([0, 1])
ax.set_xticklabels([f"¹³¹I uptake YES\n(n = {len(yes_up)})", f"¹³¹I uptake NO\n(n = {len(no_up)})"], fontsize=12, fontweight="bold")
ax.set_ylabel("DM1 axis score", fontsize=12)
ax.set_xlim(-0.5, 1.5)
ax.set_title("Panel B · DM1 axis → radioiodine uptake failure 직접 예측",
             fontsize=13.5, fontweight="bold", loc="left", color="#0F172A", pad=12)
stat_txt = (f"Cohen's d = {d_cohen:+.2f}\n"
            f"MW p (1-sided) = {p:.3f}\n"
            f"n = {len(yes_up)} uptake vs {len(no_up)} no-uptake\n"
            f"→ 방향 일치 · pilot 수준 검정")
ax.text(0.03, 0.97, stat_txt, transform=ax.transAxes, fontsize=11.5, va="top",
        fontfamily="NanumGothic",
        bbox=dict(facecolor="#FFFFFF", edgecolor="#CBD5E1", boxstyle="round,pad=0.7", linewidth=1.5))
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.grid(True, ls=":", alpha=0.3, axis="y")

fig.suptitle("Fig DEEP-3 · GSE151179 — 수술 시점 (pre-RAI) 조직으로 이후 RAI 성공 예측",
             fontsize=15, fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(f"{OUT}/fig_deep3_redifferentiation.png", bbox_inches="tight")
plt.close(fig)
print(f"[DEEP-3] Fixed · {OUT}/fig_deep3_redifferentiation.png")


# ═══════════════════════════════════════════════════════════════
# FIX 2 · Fig A2-R · Lee 2024 replication — Panel B rewrite
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(18, 9), dpi=170)

# ── Panel A · Cohort effect size forest
ax = axes[0]
cohorts = [
    ("TCGA-THCA (원 발견)",       1.64, "#0F172A", 478),
    ("★ Lee 2024 (Korean)",       0.24, "#B91C1C", 370),
    ("GSE29265",                  0.85, "#0E7490", 50),
    ("GSE33630",                  0.90, "#0E7490", 90),
    ("GSE65144",                  0.94, "#0E7490", 30),
    ("Landa 2016 (WDTC↔ATC)",     2.30, "#B45309", 88),
]
labels = [c[0] for c in cohorts]
ds     = [c[1] for c in cohorts]
colors = [c[2] for c in cohorts]
ns     = [c[3] for c in cohorts]
y_pos = np.arange(len(cohorts))[::-1]
ax.barh(y_pos, ds, color=colors, edgecolor="white", height=0.7, alpha=0.92)
for yi, lab, dv, n in zip(y_pos, labels, ds, ns):
    ax.text(dv + 0.05, yi, f"|d| = {dv:.2f}   n = {n}",
            va="center", fontsize=12, fontfamily="DejaVu Sans Mono", color="#0F172A", fontweight="bold")
ax.set_yticks(y_pos); ax.set_yticklabels(labels, fontsize=12.5)
ax.set_xlabel("|Cohen's d|  ·  DM1 axis effect size", fontsize=13)
ax.set_xlim(0, max(ds) * 1.4)
ax.axvline(1.5, color="#94A3B8", ls="--", lw=1, alpha=0.6)
ax.text(1.5, len(cohorts) - 0.5, "|d| = 1.5\n임상 유용 기준", fontsize=10, color="#64748B", style="italic", ha="center")
ax.set_title("Panel A · DM1 axis 는 6 external cohort 에서 재현됨",
             fontsize=14, fontweight="bold", loc="left", color="#0F172A", pad=12)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# ── Panel B · Data completeness matrix — clean grid instead of overlapping text
ax = axes[1]
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

# Table headers
headers = ["코호트", "n", "BRAF+ 데이터", "PFI/OS", "Interaction 검정 가능", "결과"]
col_x = [0.25, 3.15, 3.95, 5.55, 6.55, 8.65]
col_w = [2.85, 0.75, 1.55, 0.95, 2.05, 1.30]
header_y = 9.2

# Header background
ax.add_patch(plt.Rectangle((0.15, header_y - 0.35), 9.75, 0.7, facecolor="#0F172A", edgecolor="none"))
for x, w, h in zip(col_x, col_w, headers):
    ax.text(x + w/2 - 0.05, header_y, h, ha="center", va="center", fontsize=11, fontweight="bold", color="#FFFFFF")

# Rows
rows = [
    ("TCGA-THCA (prior)",       "472",  "✓",     "✓",     "직접 실행",           "p_int = 0.022 ★",   "#FEF3C7"),
    ("★ Lee 2024",              "370",  "간접",  "△",     "SubB axis proxy",      "★★★ p = 2.7 × 10⁻⁷", "#ECFDF5"),
    ("K2 PRJEB11591",           "260",  "간접",  "✗",     "DM1 prev = 25%",       "일관",              "#FFFFFF"),
    ("GSE29265",                "~50",  "간접",  "✗",     "Tertile |d|",          "|d| 재현",           "#FFFFFF"),
    ("GSE33630",                "~90",  "간접",  "✗",     "Tertile |d|",          "|d| 재현",           "#FFFFFF"),
    ("GSE65144",                "~30",  "간접",  "✗",     "Tertile |d|",          "|d| 재현",           "#FFFFFF"),
    ("Landa 2016",              " 88",  "간접",  "△",     "Histology gradient",   "|d| = 2.3",         "#FFFFFF"),
]
row_h = 0.75
for i, row in enumerate(rows):
    y = header_y - 0.75 - i*row_h
    bg = row[-1]
    ax.add_patch(plt.Rectangle((0.15, y - 0.32), 9.75, 0.65, facecolor=bg, edgecolor="#CBD5E1", linewidth=0.5))
    for x, w, txt in zip(col_x, col_w, row[:-1]):
        fs = 11 if i > 0 else 11
        fw = "bold" if i == 1 else "normal"
        color = "#7C2D12" if i == 1 else "#0F172A"
        ax.text(x + w/2 - 0.05, y, txt, ha="center", va="center", fontsize=fs, fontweight=fw, color=color)

# Conclusion box
concl_y = 2.2
ax.add_patch(plt.Rectangle((0.5, 0.6), 9.0, concl_y - 0.4, facecolor="#ECFDF5", edgecolor="#047857", linewidth=2))
ax.text(5.0, 2.0, "결론 · DM1 axis 는 6 external cohort 에서 재현 (Lee 2024 p = 2.7 × 10⁻⁷)",
        ha="center", va="center", fontsize=12, fontweight="bold", color="#047857")
ax.text(5.0, 1.5, "직접 DM1 × BRAF interaction Cox 는 TCGA 만 실행 가능 (외부 BRAF+PFI 조합 없음)",
        ha="center", va="center", fontsize=11, color="#065F46")
ax.text(5.0, 1.0, "→ 분당 prospective (n = 200) 가 interaction 확정 검정의 primary target",
        ha="center", va="center", fontsize=11, style="italic", color="#065F46")

ax.set_title("Panel B · Cohort 별 데이터 완성도 매트릭스 · Interaction 재현 status",
             fontsize=14, fontweight="bold", loc="left", color="#0F172A", pad=12)

fig.suptitle("Fig A2-R · A2 Interaction Replication — Lee 2024 (n = 370) 축 재현 ★★★",
             fontsize=15.5, fontweight="bold", y=1.00)
plt.tight_layout()
fig.savefig(f"{OUT}/fig_a2_replication.png", bbox_inches="tight")
plt.close(fig)
print(f"[A2-R] Fixed · {OUT}/fig_a2_replication.png")
