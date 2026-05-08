"""Add 8 more explanatory figures for the expanded one-page audit."""
from pathlib import Path
import os
import numpy as np
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Circle, Wedge

for _p in ["/home/seungho/.local/share/fonts/NanumGothic-Regular.ttf",
           "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
           "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"]:
    if os.path.exists(_p):
        try: fm.fontManager.addfont(_p)
        except Exception: pass

OUT = Path(__file__).resolve().parent
plt.rcParams.update({
    "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": ["NanumGothic", "Noto Sans CJK JP", "DejaVu Sans"],
    "axes.unicode_minus": False,
})
RED, BLUE, GREY, ORANGE, GREEN, GOLD, NAVY = "#c0392b", "#2c5e9c", "#7f8fa6", "#e67e22", "#3C6B4F", "#d4a017", "#1f3a5f"

def box(ax, x, y, w, h, text, color, fontsize=9, fontweight="normal"):
    r = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04", linewidth=1.2,
                       edgecolor="black", facecolor=color)
    ax.add_patch(r)
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fontsize,
            fontweight=fontweight, wrap=True)
def arrow(ax, x1, y1, x2, y2, color="black", lw=1.5):
    a = FancyArrowPatch((x1,y1),(x2,y2), arrowstyle="-|>", mutation_scale=18,
                        color=color, lw=lw)
    ax.add_patch(a)

# ============================================================
# fig11 — narrative storyline arc
# ============================================================
fig, ax = plt.subplots(figsize=(15, 6))
ax.set_xlim(0,15); ax.set_ylim(0,6); ax.axis("off")
ax.set_title("Paper 1 — narrative storyline arc (Hook → Mechanism → Clinical reframe)",
             fontsize=14, fontweight="bold")
beats = [
    ("HOOK\n임상 문제", "RAI 실패는\n수개월~수년 후\n인지", "#FFE5D9", 0.5),
    ("AXIS\n분자 정체성", "분화도 axis 가\ndriver 와\n직교", "#FFF4D6", 3.0),
    ("READOUT\n측정 도구", "8-gene RT-qPCR\ncompact panel", "#E8F4FD", 5.5),
    ("REPLICATION\n다중 코호트", "GPL570 4× +\nKorean 2× +\nLanda advanced", "#D4F4DD", 8.0),
    ("MECHANISM\n생물학적 근거", "lineage TF 침묵 +\npromoter methylation +\nDNMT/STAT3 활성", GOLD, 10.5),
    ("REFRAME\n임상 의미", "post-surgery\nearly triage scaffold\n(hypothesis only)", "#cfd8dc", 13.0),
]
for h, sub, c, x in beats:
    box(ax, x, 2.5, 2.0, 2.5, f"{h}\n\n{sub}", c, 9, "bold")
for x in [2.5, 5.0, 7.5, 10.0, 12.5]:
    arrow(ax, x, 3.75, x+0.5, 3.75, "#7B1F2A", 2)
ax.text(0.5, 1.3, "이 6 단계가 manuscript 의 narrative spine.", fontsize=11, color=NAVY)
ax.text(0.5, 0.7, "★ 각 단계는 자체로 figure 1개 이상 + 별도 evidence; 빠지면 '무엇이 novel 한가' 가 reviewer 에게 안 보임.", fontsize=10, color=RED, style="italic")
plt.tight_layout(); plt.savefig(OUT/"fig11_narrative_storyline.png", dpi=150, bbox_inches="tight"); plt.close()

# ============================================================
# fig12 — RAI biology pathway (clinically familiar)
# ============================================================
fig, ax = plt.subplots(figsize=(14, 7))
ax.set_xlim(0,14); ax.set_ylim(0,7); ax.axis("off")
ax.set_title("Fig 12 — 갑상선 호르몬 합성 + RAI 흡수 분자 pathway (8-gene + NONOVERLAP 위치)",
             fontsize=12, fontweight="bold")
# follicle schematic
ax.add_patch(Circle((4.5, 3.5), 1.5, facecolor="#FFF8E1", edgecolor="black", lw=1.5))
ax.text(4.5, 4.7, "Thyroid follicle\n(colloid lumen)", ha="center", fontsize=9, fontweight="bold")
# basolateral side (blood)
ax.text(0.5, 5.5, "혈류 측 (basolateral)", fontsize=9, color=BLUE, fontweight="bold")
ax.text(0.5, 5.2, "I⁻ uptake", fontsize=8)
# arrows for I- uptake
arrow(ax, 1.0, 4.0, 2.5, 3.7, BLUE, 2)
ax.text(1.5, 4.2, "I⁻", fontsize=10, color=BLUE, fontweight="bold")
box(ax, 1.5, 3.0, 1.4, 0.6, "NIS\n(SLC5A5)", "#c0392b", 9, "bold")  # RAI_8
# apical side
ax.text(8.5, 5.5, "정점 측 (apical)", fontsize=9, color=ORANGE, fontweight="bold")
arrow(ax, 6.0, 3.7, 7.5, 4.0, ORANGE, 2)
ax.text(7.0, 4.2, "I⁻", fontsize=10, color=ORANGE)
box(ax, 6.5, 2.7, 1.5, 0.6, "Pendrin\n(SLC26A4)", BLUE, 8.5, "bold")  # NONOVERLAP
# H2O2 supply
box(ax, 7.5, 1.5, 1.5, 0.6, "DUOX1/2", BLUE, 8.5, "bold")  # NONOVERLAP
arrow(ax, 7.8, 2.1, 7.0, 2.7, BLUE, 1.5)
ax.text(8.2, 2.0, "H₂O₂", fontsize=8)
# TPO
box(ax, 4.5, 1.0, 1.5, 0.6, "TPO", "#c0392b", 9, "bold")  # RAI_8
arrow(ax, 5.2, 1.6, 5.2, 2.5, "#c0392b", 1.5)
ax.text(3.0, 1.3, "I⁻ 산화 →", fontsize=8)
# TG
box(ax, 9.5, 3.0, 1.5, 0.6, "TG\n(thyroglobulin)", "#c0392b", 8.5, "bold")
arrow(ax, 10.2, 3.0, 10.2, 4.5, "#c0392b", 1.5)
ax.text(10.5, 3.7, "TG-I", fontsize=8)
# T4/T3 release
ax.text(11.0, 5.5, "T4 / T3 분비", fontsize=10, color=GOLD, fontweight="bold")
arrow(ax, 11.5, 5.0, 12.5, 5.5, GOLD, 2)
# IYD recycle
box(ax, 6.0, 0.3, 1.5, 0.6, "IYD\n(요오드 재활용)", BLUE, 8.5, "bold")  # NONOVERLAP
# DIO1/DIO2
box(ax, 11.5, 1.5, 1.5, 0.6, "DIO1 (RAI_8)", "#c0392b", 8.5, "bold")
box(ax, 11.5, 2.5, 1.5, 0.6, "DIO2 (NONOVERLAP)", BLUE, 8, "bold")
ax.text(12.2, 0.9, "T4 → T3", fontsize=9, color=GOLD)

# Lineage TF (regulatory layer above)
box(ax, 0.5, 6.0, 1.6, 0.7, "FOXE1 (TTF-2)", "#94d4a4", 8, "bold")
box(ax, 2.3, 6.0, 1.6, 0.7, "NKX2-1 (TTF-1)", "#94d4a4", 8, "bold")
box(ax, 4.1, 6.0, 1.0, 0.7, "PAX8", "#94d4a4", 9, "bold")
box(ax, 5.3, 6.0, 1.0, 0.7, "HHEX", "#cfd8dc", 9, "bold")  # NONOVERLAP TF
ax.text(0.5, 6.85, "★ Lineage TF backbone (regulatory layer 위에서)", fontsize=10, fontweight="bold", color=GREEN)

# arrows down to enzymes
for x in [1.3, 3.1, 4.6, 5.8]:
    arrow(ax, x, 6.0, x+0.0, 5.5, "#7f8fa6", 1.0)

# TSHR receptor
box(ax, 12.0, 6.0, 1.5, 0.7, "TSHR receptor", "#c0392b", 8.5, "bold")
ax.text(12.7, 6.85, "(RAI_8)", fontsize=8)

# legend
ax.text(0.5, 0.0, "[RED] = RAI_8 panel (8 genes)   |   [BLUE] = THYROID_NONOVERLAP panel (8 genes; zero gene overlap)   |   [GREEN] = TF backbone (TF_collapse)",
        fontsize=9.5)
plt.tight_layout(); plt.savefig(OUT/"fig12_rai_biology_pathway.png", dpi=150, bbox_inches="tight"); plt.close()

# ============================================================
# fig13 — mechanism cascade (lineage collapse → DNMT → epigenetic silencing → RAI fail)
# ============================================================
fig, ax = plt.subplots(figsize=(15, 5.5))
ax.set_xlim(0,15); ax.set_ylim(0,6); ax.axis("off")
ax.set_title("Fig 13 — Mechanism cascade: lineage TF collapse → epigenetic silencing → RAI 흡수 불능",
             fontsize=12, fontweight="bold")

# Layer 1 — TF collapse
box(ax, 0.5, 4.0, 2.5, 1.2, "Lineage TF\nbackbone collapse\n(FOXE1↓ NKX2-1↓\nPAX8↓ HHEX↓)", "#94d4a4", 8.5, "bold")
arrow(ax, 3.0, 4.6, 4.0, 4.6)
# Layer 2 — DNMT/STAT3/AP1 activation
box(ax, 4.0, 4.0, 2.5, 1.2, "DNMT1/3B ↑\nSTAT3 ↑ FOSL1 ↑\nJUNB ↑\n(de-diff TF program)", GOLD, 8.5, "bold")
arrow(ax, 6.5, 4.6, 7.5, 4.6)
# Layer 3 — methylation
box(ax, 7.5, 4.0, 2.5, 1.2, "Promoter\nhypermethylation\n(TPO d=2.30, β +52%)", "#FFE5D9", 8.5, "bold")
arrow(ax, 10.0, 4.6, 11.0, 4.6)
# Layer 4 — RAI fail
box(ax, 11.0, 4.0, 3.5, 1.2, "갑상선 호르몬\n합성 machinery 침묵\n→ RAI 흡수 불능 (RAI refractory)\n→ DM1 임상 표현형", RED, 9, "bold")

# 16-gene readout layer (read out)
box(ax, 1.0, 1.5, 13, 1.5,
    "★ READOUT layer: RAI_8 panel + THYROID_NONOVERLAP panel\n전체 machinery 의 전 layer 가 동시에 침묵하는 것이 측정됨\n(zero-overlap cross-panel ρ = +0.44 in TCGA; ρ ≤ −0.84 with DM1_like in 4 GPL570 cohorts)",
    "#E8F4FD", 10, "bold")
arrow(ax, 12.7, 4.0, 12.7, 3.0, "#7f8fa6", 1)
ax.text(7, 3.3, "↑ 측정", fontsize=9, color="#7f8fa6")

ax.text(0.5, 0.5, "★ 임상 implication: DM1 = 4 가지 동시 발생 → fusion 7×↑ + TPO methylation + RAI 불응 + survival HR 2.5",
        fontsize=10, fontweight="bold", color=RED)
plt.tight_layout(); plt.savefig(OUT/"fig13_mechanism_cascade.png", dpi=150, bbox_inches="tight"); plt.close()

# ============================================================
# fig14 — clinical timeline (post-surgery year 1-3)
# ============================================================
fig, ax = plt.subplots(figsize=(14, 5))
ax.set_xlim(0,14); ax.set_ylim(0,5); ax.axis("off")
ax.set_title("Fig 14 — 수술 후 임상 timeline + DM1 readout 의 가설적 개입 시점",
             fontsize=12, fontweight="bold")
# timeline base
ax.plot([0.5, 13.5], [3.0, 3.0], color="black", lw=2)
for x, label in [(0.5, "Post-op\n0w"), (3.0, "RAI dose\n6-8w"), (5.0, "Tg follow-up\n6mo"), (8.0, "Tg trend\n12mo"), (11.0, "Recurrence\n2-3yr"), (13.0, "Aggressive\n>3yr")]:
    ax.plot([x],[3.0], "ko", markersize=8)
    ax.text(x, 2.5, label, ha="center", fontsize=8.5)
# current (red flags above)
ax.text(0.5, 4.5, "현재 standard 흐름 (RAI 실패 인지 지연)", fontsize=10, color=RED, fontweight="bold")
arrow(ax, 5.0, 3.4, 8.0, 3.8, RED, 1.5)
ax.text(6.5, 3.9, "Tg ↑ 인지\n(6-12mo+)", fontsize=8, color=RED, ha="center")
arrow(ax, 8.0, 3.4, 11.0, 3.7, RED, 1.5)
ax.text(9.5, 3.8, "분자 review", fontsize=8, color=RED, ha="center")
# proposed (blue interventions below)
ax.text(0.5, 1.7, "제안된 가설적 triage scaffold (NOT validated)", fontsize=10, color=BLUE, fontweight="bold")
arrow(ax, 0.5, 2.6, 1.5, 2.0, BLUE, 1.5)
ax.text(1.0, 1.3, "8-gene readout\n(post-op 0w)", fontsize=8, color=BLUE, ha="center")
ax.text(1.0, 0.8, "→ DM1 high-risk flag", fontsize=8, color=NAVY, ha="center")
arrow(ax, 1.5, 1.5, 5.0, 1.5, BLUE, 2)
ax.text(3.0, 1.7, "earlier review", fontsize=8, color=BLUE, ha="center")
ax.text(11.0, 0.5, "기대: 시간 절약 (months)", fontsize=10, color=NAVY, fontweight="bold", ha="center")
ax.text(0.5, 0.0, "★ scaffold 만; prospective trial 후에 임상 implementation 가능. 현재 기준 RAI 표준 절차를 변경 권고하지 않음.",
        fontsize=9, color="#7B1F2A", style="italic")
plt.tight_layout(); plt.savefig(OUT/"fig14_clinical_timeline.png", dpi=150, bbox_inches="tight"); plt.close()

# ============================================================
# fig15 — Korean cohort applicability (East-Asian generalizability)
# ============================================================
fig, ax = plt.subplots(1, 2, figsize=(14, 5))
ax[0].set_title("DM1 prevalence — 대륙별 비교", fontsize=11)
labels = ["TCGA-THCA\n(서양 dominant)", "K2 SNU-GMI\n(Korean)", "Lee 2024\n(Korean)", "GSE286332\n(Korean PTC+HT)"]
prev = [28.4, 37.8, 35.0, 33.0]  # approximate; some descriptive
colors_b = [GREY, RED, RED, RED]
ax[0].bar(labels, prev, color=colors_b, edgecolor="black")
ax[0].set_ylabel("DM1 prevalence (%)"); ax[0].set_ylim(0, 50)
for i, v in enumerate(prev):
    ax[0].text(i, v+1, f"{v}%", ha="center", fontsize=9, fontweight="bold")
ax[0].axhline(33, color="grey", linestyle=":")
ax[0].text(0, 45, "★ Korean ≈ 35%; TCGA ≈ 28% (close)", fontsize=9, color=NAVY)

ax[1].set_title("Hashimoto-like proportion (DM2 axis)", fontsize=11)
ht = [18, 28, 22, 100]
ht_labels = ["TCGA-THCA", "Lee 2024 (632)", "K2 (260)", "GSE286332 (PTC+HT only)"]
ax[1].barh(ht_labels, ht, color=[GREY, BLUE, BLUE, ORANGE], edgecolor="black")
ax[1].set_xlabel("Hashimoto-like %"); ax[1].set_xlim(0, 110)
for i, v in enumerate(ht):
    ax[1].text(v+2, i, f"{v}%", va="center", fontsize=9, fontweight="bold")
ax[1].text(0, -0.7, "★ TCGA 18% ≈ Korean 22-28% — same axis, 같은 distribution", fontsize=9, color=NAVY)

plt.suptitle("Fig 15 — East-Asian generalizability board", fontsize=12, fontweight="bold")
plt.tight_layout(); plt.savefig(OUT/"fig15_korean_cohort_applicability.png", dpi=150, bbox_inches="tight"); plt.close()

# ============================================================
# fig16 — 2-layer architecture diagram (discovery vs deployment)
# ============================================================
fig, ax = plt.subplots(figsize=(13, 6.5))
ax.set_xlim(0,13); ax.set_ylim(0,7); ax.axis("off")
ax.set_title("Fig 16 — 2-layer architecture: discovery axis (위) vs deployment readout (아래)",
             fontsize=12, fontweight="bold")
# discovery layer
ax.add_patch(Rectangle((0.3, 4.0), 12.4, 2.5, facecolor="#fdf5ef", edgecolor=NAVY, lw=2, linestyle="--"))
ax.text(0.5, 6.2, "DISCOVERY LAYER (axis 자체 정의)", fontsize=12, fontweight="bold", color=NAVY)
box(ax, 1.0, 4.5, 3.5, 1.0, "Pan-genome\ntop-5000 (MAD)\nARI = 0.92", "#E8F4FD", 9.5, "bold")
box(ax, 5.0, 4.5, 3.5, 1.0, "TIERA67\n(literature 67-gene)\nARI = 0.90", "#FFF4D6", 9.5, "bold")
box(ax, 9.0, 4.5, 3.5, 1.0, "TDS-16\n(canonical Yoo 16)\nAUC = 0.975", "#cfd8dc", 9.5, "bold")
ax.text(0.5, 3.7, "→ 어느 path 로도 같은 axis 도달 (data ≈ literature ≈ canonical)", fontsize=9, color=NAVY, style="italic")

# deployment layer
ax.add_patch(Rectangle((0.3, 0.5), 12.4, 2.5, facecolor="#fff0f0", edgecolor=RED, lw=2, linestyle="--"))
ax.text(0.5, 2.7, "DEPLOYMENT LAYER (clinical practical readout)", fontsize=12, fontweight="bold", color=RED)
box(ax, 2.0, 1.0, 4.0, 1.2, "RAI_8 (8-gene)\nAUC 0.962  /  ΔAUC vs 16 = 0.013\nRT-qPCR / NanoString viable", RED, 10, "bold")
box(ax, 7.0, 1.0, 4.0, 1.2, "ARI 0.49 vs reference\n(intentional compactness cost;\nnot discovery metric)", "#FFE5D9", 10, "bold")

# vertical separator + arrow
arrow(ax, 6.5, 4.0, 6.5, 2.3, NAVY, 2.5)
ax.text(6.6, 3.2, "compress", fontsize=10, color=NAVY, fontweight="bold")

ax.text(0.5, 0.1, "★ ARI 0.49 = 'unsupervised re-discovery' metric, AUC 0.962 = 'supervised classification' metric. 다른 질문에 답함.",
        fontsize=9, color=RED, style="italic")
plt.tight_layout(); plt.savefig(OUT/"fig16_two_layer_architecture.png", dpi=150, bbox_inches="tight"); plt.close()

# ============================================================
# fig17 — what's lacking / future-validation board
# ============================================================
items = [
    ("Prospective triage trial", "RAI failure flag 능력 prospective 검증", "very_high"),
    ("Korean Bundang FFPE cohort", "Yu 교수 outreach 응답 wait", "high"),
    ("Sub-A/sub-B re-derivation", "master DM1 universe 위에서 재계산 필요", "high"),
    ("Spatial cohort larger n", "GSE250521 16-sample power 한계", "medium"),
    ("EGAS/EGAD verify + apply", "사용처 verify 후 application 결정", "medium"),
    ("Pozdeyev 2018 raw access", "n=779 advanced 데이터 PI request", "medium"),
    ("decitabine + I-131 mechanism trial", "epigenetic re-induction prospective", "high"),
    ("selpercatinib reflex prospective", "81.8% TCGA RET-fusion+ 적용성 검증", "high"),
    ("KEYNOTE-158 raw access", "Merck IIS — months", "low (not Paper 1)"),
    ("Functional perturbation of TF backbone", "causal proof requires CRISPR / lineage rescue", "low (out-of-scope)"),
]
risk_color = {"low (not Paper 1)":"#94d4a4","low (out-of-scope)":"#94d4a4",
              "medium":"#f9c863","high":ORANGE,"very_high":RED}
fig, ax = plt.subplots(figsize=(13, 7))
ax.axis("off")
ax.set_title("Fig 17 — What's lacking / future-validation board (Paper 1 너머)", fontsize=12, fontweight="bold")
for i, (n, desc, r) in enumerate(items):
    y = 0.92 - i*0.085
    ax.add_patch(Rectangle((0.02, y-0.025), 0.15, 0.060, transform=ax.transAxes,
                           facecolor=risk_color[r], edgecolor="black"))
    ax.text(0.095, y+0.005, r.upper().replace("_"," "), fontsize=8, ha="center",
            color="white" if r=="very_high" or r=="high" else "black", transform=ax.transAxes,
            fontweight="bold")
    ax.text(0.20, y+0.015, n, fontsize=11, transform=ax.transAxes, fontweight="bold")
    ax.text(0.20, y-0.012, desc, fontsize=9, transform=ax.transAxes, color="#555")
plt.tight_layout(); plt.savefig(OUT/"fig17_future_validation_board.png", dpi=150, bbox_inches="tight"); plt.close()

# ============================================================
# fig18 — what each main figure means + clinical interpretation summary
# ============================================================
fig, ax = plt.subplots(figsize=(14, 7))
ax.axis("off")
ax.set_title("Fig 18 — 5 main figures 의 'why it matters' summary", fontsize=12, fontweight="bold")
mains = [
    ("F1\nclinical workflow", "현재 RAI 실패 인지 지연 + 가설적 triage scaffold",
     "수술 직후 분자 위험 측정 도구 부재 — 임상의가 즉시 이해 가능", "RAI predictor / 치료선택 도구 아님"),
    ("F2\n8-gene readout +\nARI ladder", "8-gene 은 readout, axis 는 TIERA67/pan-genome 이 정의",
     "8-gene 의 modest ARI 0.49 가 약점이 아닌 의도된 compactness", "panel 이 axis 를 *발견* 했다는 인상 금지"),
    ("F3\nGPL570 4-cohort\nexternal validation", "RAI_8 + NONOVERLAP zero-overlap co-collapse, ρ ≤ -0.84",
     "panel-overlap artifact 가설 직접 차단 — 가장 강한 §4.3 defense", "단일 cohort cherry-pick 으로 보이면 안 됨"),
    ("F4\nPFI + age-adj Cox\n(OS sparse)", "OS event 3.4% 너무 적음 → PFI primary; age + stage adjusted",
     "TCGA-THCA OS 한계 honest disclosure → reviewer trust", "OS 'swap' 무시하면 age confound 의심"),
    ("F5\nHM450 promoter\nmethylation", "DM1 에서 lineage gene promoter β +52% (TPO d=2.30)",
     "분화도 침묵의 mechanistic 근거 + decitabine + I-131 trial 적격성 합리화", "causality 입증 아님 ('consistent with')"),
]
for i, (fig_id, what, why_imp, what_not) in enumerate(mains):
    y = 0.88 - i*0.18
    ax.add_patch(Rectangle((0.02, y-0.085), 0.13, 0.16, transform=ax.transAxes,
                           facecolor="#FFE5D9", edgecolor="black"))
    ax.text(0.085, y, fig_id, fontsize=10, ha="center", transform=ax.transAxes, fontweight="bold")
    ax.text(0.17, y+0.04, "WHAT", fontsize=8, transform=ax.transAxes, color=NAVY, fontweight="bold")
    ax.text(0.23, y+0.04, what, fontsize=10, transform=ax.transAxes)
    ax.text(0.17, y, "WHY", fontsize=8, transform=ax.transAxes, color=GREEN, fontweight="bold")
    ax.text(0.23, y, why_imp, fontsize=9.5, transform=ax.transAxes, color="#222")
    ax.text(0.17, y-0.04, "GUARD", fontsize=8, transform=ax.transAxes, color=RED, fontweight="bold")
    ax.text(0.23, y-0.04, what_not, fontsize=9.5, transform=ax.transAxes, color=RED, style="italic")
plt.tight_layout(); plt.savefig(OUT/"fig18_main_figures_meaning.png", dpi=150, bbox_inches="tight"); plt.close()

print("[done] 8 extra figures")
