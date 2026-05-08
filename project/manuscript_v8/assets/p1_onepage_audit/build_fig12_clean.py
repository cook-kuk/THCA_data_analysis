"""Redesigned Fig 12 — clean, readable RAI biology pathway infographic."""
from pathlib import Path
import os
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Circle, Wedge

for _p in ["/home/seungho/.local/share/fonts/NanumGothic-Regular.ttf",
           "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"]:
    if os.path.exists(_p):
        try: fm.fontManager.addfont(_p)
        except Exception: pass

plt.rcParams.update({
    "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": ["NanumGothic", "Noto Sans CJK JP", "DejaVu Sans"],
    "axes.unicode_minus": False,
})

OUT = Path(__file__).resolve().parent

# Color palette — light pastel for fills, strong for borders/text
RAI_FILL  = "#fff0ed"; RAI_BD  = "#c0392b"
NON_FILL  = "#eef4fb"; NON_BD  = "#2c5e9c"
TF_FILL   = "#eaf6ec"; TF_BD   = "#3C6B4F"
ACCENT    = "#7B1F2A"
INK       = "#1a1a1a"

fig, ax = plt.subplots(figsize=(16, 9))
ax.set_xlim(0, 16); ax.set_ylim(0, 9); ax.axis("off")

# ===== title =====
ax.text(8, 8.55, "갑상선 호르몬 합성 + RAI 흡수 분자 pathway",
        ha="center", fontsize=16, fontweight="bold", color=ACCENT)
ax.text(8, 8.18, "8-gene RAI_8 (red) + THYROID_NONOVERLAP (blue) + TF backbone (green) 의 위치",
        ha="center", fontsize=11, color="#555")

# ===== Layer 1 (top) — TF backbone =====
ax.add_patch(Rectangle((0.4, 6.5), 15.2, 1.2, facecolor=TF_FILL, edgecolor=TF_BD,
                       linewidth=1.5, linestyle="--"))
ax.text(0.7, 7.45, "REGULATORY LAYER · Lineage transcription factors (TF backbone)",
        fontsize=11, fontweight="bold", color=TF_BD)

def tf_box(x, name, sub):
    bx = FancyBboxPatch((x, 6.7), 1.85, 0.55, boxstyle="round,pad=0.04",
                        linewidth=1.5, edgecolor=TF_BD, facecolor="#fff")
    ax.add_patch(bx)
    ax.text(x+0.92, 7.05, name, ha="center", fontsize=10.5, fontweight="bold", color=INK)
    ax.text(x+0.92, 6.83, sub, ha="center", fontsize=8.5, color="#555", style="italic")

tf_box(1.0, "FOXE1", "TTF-2 / FKHL15")
tf_box(3.5, "NKX2-1", "TTF-1")
tf_box(6.0, "PAX8", "PAX-related")
tf_box(8.5, "HHEX", "endoderm TF")
ax.text(11.5, 7.05, "→ NIS / TG / TPO / TSHR promoter binding",
        fontsize=10, color=TF_BD, style="italic")

# down arrow to follicle
ax.annotate("", xy=(8, 6.0), xytext=(8, 6.5),
            arrowprops=dict(arrowstyle="-|>", color=TF_BD, lw=1.5))
ax.text(8.15, 6.25, "drives", fontsize=9, color=TF_BD, style="italic")

# ===== Layer 2 (middle) — Follicle cell =====
# colloid lumen
ax.add_patch(Circle((8, 3.5), 1.7, facecolor="#fff8d6", edgecolor="#b08c1f", lw=2))
ax.text(8, 3.9, "Colloid", ha="center", fontsize=11, fontweight="bold", color="#7a5f0a")
ax.text(8, 3.5, "(lumen)", ha="center", fontsize=9, color="#7a5f0a")
ax.text(8, 3.1, "TG-I", ha="center", fontsize=11, fontweight="bold", color=RAI_BD)
ax.text(8, 2.85, "iodinated", ha="center", fontsize=8, color=RAI_BD, style="italic")

# basolateral side (blood)
ax.add_patch(Rectangle((0.4, 4.6), 4.5, 0.4, facecolor="#ffe1de", edgecolor=RAI_BD, lw=1.0))
ax.text(0.55, 4.8, "BLOOD (basolateral)", fontsize=10, fontweight="bold", color=RAI_BD)
# I- arrow + NIS
ax.annotate("I⁻", xy=(5.5, 4.0), xytext=(2.0, 4.5),
            fontsize=14, color=RAI_BD, fontweight="bold",
            arrowprops=dict(arrowstyle="-|>", color=RAI_BD, lw=2.0))

bx = FancyBboxPatch((4.6, 3.7), 1.2, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=RAI_BD, facecolor=RAI_FILL)
ax.add_patch(bx)
ax.text(5.2, 4.0, "NIS", ha="center", fontsize=12, fontweight="bold", color=RAI_BD)
ax.text(5.2, 3.78, "(SLC5A5)", ha="center", fontsize=8, color="#444")

# apical side (toward colloid)
ax.add_patch(Rectangle((11.1, 4.6), 4.5, 0.4, facecolor="#fcf3e7", edgecolor="#b08c1f", lw=1.0))
ax.text(11.25, 4.8, "APICAL membrane", fontsize=10, fontweight="bold", color="#7a5f0a")
# Pendrin (NONOVERLAP)
bx = FancyBboxPatch((10.0, 3.7), 1.4, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=NON_BD, facecolor=NON_FILL)
ax.add_patch(bx)
ax.text(10.7, 4.0, "Pendrin", ha="center", fontsize=11, fontweight="bold", color=NON_BD)
ax.text(10.7, 3.78, "(SLC26A4)", ha="center", fontsize=8, color="#444")

# I- arrow from follicle interior to colloid
ax.annotate("", xy=(11.3, 3.7), xytext=(9.7, 3.7),
            arrowprops=dict(arrowstyle="-|>", color=NON_BD, lw=2.0))
ax.text(10.3, 3.55, "I⁻", fontsize=11, color=NON_BD, fontweight="bold")

# TPO (RAI_8)
bx = FancyBboxPatch((6.5, 1.6), 1.2, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=RAI_BD, facecolor=RAI_FILL)
ax.add_patch(bx)
ax.text(7.1, 1.9, "TPO", ha="center", fontsize=12, fontweight="bold", color=RAI_BD)
ax.text(7.1, 1.68, "(peroxidase)", ha="center", fontsize=8, color="#444")

# DUOX1/2 (NONOVERLAP) — H2O2 supplier
bx = FancyBboxPatch((8.4, 1.6), 1.5, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=NON_BD, facecolor=NON_FILL)
ax.add_patch(bx)
ax.text(9.15, 1.9, "DUOX1 / 2", ha="center", fontsize=10.5, fontweight="bold", color=NON_BD)
ax.text(9.15, 1.68, "(H₂O₂ supply)", ha="center", fontsize=8, color="#444")

# arrow DUOX → TPO
ax.annotate("", xy=(7.7, 1.9), xytext=(8.4, 1.9),
            arrowprops=dict(arrowstyle="-|>", color=NON_BD, lw=1.5))
ax.text(8.05, 2.05, "H₂O₂", fontsize=8, color=NON_BD)
# arrow TPO upward to colloid
ax.annotate("", xy=(7.5, 3.0), xytext=(7.1, 2.2),
            arrowprops=dict(arrowstyle="-|>", color=RAI_BD, lw=1.5))
ax.text(7.0, 2.6, "I⁻ 산화", fontsize=8, color=RAI_BD)

# TG (RAI_8) — substrate
bx = FancyBboxPatch((4.5, 1.6), 1.4, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=RAI_BD, facecolor=RAI_FILL)
ax.add_patch(bx)
ax.text(5.2, 1.9, "TG", ha="center", fontsize=12, fontweight="bold", color=RAI_BD)
ax.text(5.2, 1.68, "(thyroglobulin)", ha="center", fontsize=8, color="#444")
ax.annotate("", xy=(6.5, 2.0), xytext=(5.9, 2.0),
            arrowprops=dict(arrowstyle="-|>", color="#777", lw=1.0))
ax.annotate("", xy=(7.0, 3.0), xytext=(5.5, 2.2),
            arrowprops=dict(arrowstyle="-|>", color=RAI_BD, lw=1.0, linestyle=":"))

# IYD (NONOVERLAP) — recycle
bx = FancyBboxPatch((10.0, 1.6), 1.4, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=NON_BD, facecolor=NON_FILL)
ax.add_patch(bx)
ax.text(10.7, 1.9, "IYD", ha="center", fontsize=11, fontweight="bold", color=NON_BD)
ax.text(10.7, 1.68, "(I⁻ recycle)", ha="center", fontsize=8, color="#444")

# TSHR (RAI_8) — receptor (basolateral)
bx = FancyBboxPatch((1.0, 3.7), 1.4, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=RAI_BD, facecolor=RAI_FILL)
ax.add_patch(bx)
ax.text(1.7, 4.0, "TSHR", ha="center", fontsize=12, fontweight="bold", color=RAI_BD)
ax.text(1.7, 3.78, "(TSH 수용체)", ha="center", fontsize=8, color="#444")

# Output arrow — T4/T3
ax.annotate("", xy=(15.3, 5.0), xytext=(11.5, 4.5),
            arrowprops=dict(arrowstyle="-|>", color="#7a5f0a", lw=2.5))
ax.text(13.5, 5.15, "T4 / T3 분비", fontsize=12, color="#7a5f0a", fontweight="bold")
ax.text(13.5, 4.85, "(혈류로)", fontsize=9, color="#7a5f0a")

# DIO1 (RAI_8) and DIO2 (NONOVERLAP)
bx = FancyBboxPatch((13.5, 2.8), 1.4, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=RAI_BD, facecolor=RAI_FILL)
ax.add_patch(bx)
ax.text(14.2, 3.1, "DIO1", ha="center", fontsize=11, fontweight="bold", color=RAI_BD)
ax.text(14.2, 2.88, "(deiodinase 1)", ha="center", fontsize=8, color="#444")

bx = FancyBboxPatch((13.5, 1.9), 1.4, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=NON_BD, facecolor=NON_FILL)
ax.add_patch(bx)
ax.text(14.2, 2.2, "DIO2", ha="center", fontsize=11, fontweight="bold", color=NON_BD)
ax.text(14.2, 1.98, "(deiodinase 2)", ha="center", fontsize=8, color="#444")

# TFF3 (NONOVERLAP)
bx = FancyBboxPatch((1.0, 1.6), 1.4, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=NON_BD, facecolor=NON_FILL)
ax.add_patch(bx)
ax.text(1.7, 1.9, "TFF3", ha="center", fontsize=11, fontweight="bold", color=NON_BD)
ax.text(1.7, 1.68, "(secretory)", ha="center", fontsize=8, color="#444")

# GLIS3 (NONOVERLAP TF)
bx = FancyBboxPatch((2.6, 1.6), 1.5, 0.6, boxstyle="round,pad=0.04",
                    linewidth=2, edgecolor=NON_BD, facecolor=NON_FILL)
ax.add_patch(bx)
ax.text(3.35, 1.9, "GLIS3", ha="center", fontsize=11, fontweight="bold", color=NON_BD)
ax.text(3.35, 1.68, "(supportive TF)", ha="center", fontsize=8, color="#444")

# ===== Legend =====
ax.add_patch(Rectangle((0.4, 0.1), 15.2, 0.7, facecolor="#fafafa",
                       edgecolor="#aaa", lw=0.8))
def legend_chip(x, color_fill, color_bd, label):
    bx = FancyBboxPatch((x, 0.30), 0.55, 0.30, boxstyle="round,pad=0.04",
                        linewidth=2, edgecolor=color_bd, facecolor=color_fill)
    ax.add_patch(bx)
    ax.text(x+0.7, 0.45, label, fontsize=10.5, color=color_bd, fontweight="bold", va="center")

legend_chip(0.7, RAI_FILL, RAI_BD,
            "RAI_8 panel (8 genes; main deployable readout)")
legend_chip(7.0, NON_FILL, NON_BD,
            "THYROID_NONOVERLAP (8 genes; zero overlap with RAI_8)")
legend_chip(13.0, TF_FILL, TF_BD, "TF backbone")

plt.tight_layout()
plt.savefig(OUT/"fig12_rai_biology_pathway.png", dpi=160, bbox_inches="tight")
plt.close()
print("Fig 12 redesigned")
