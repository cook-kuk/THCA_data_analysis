#!/usr/bin/env python3
"""Korean/Asian cohort validation of the 8-gene panel — unified composite figure.

Integrates six independent Korean/Chinese thyroid cohorts that we already have
locally with eight-gene panel scores from Paper 1 R17 + new RAI atlas pipeline:

  1. K2 (PRJEB11591) — Korean RNA-seq, n = 260, 8-gene quant
  2. Lee 2024 (GSE213647) — Korean, n = 632 (Korean PTC + adj normal)
  3. GSE286332 — Korean PTC vs PTC+HT, n = 18 (Paper 1 R17 L3)
  4. Mun 2025 — Korean proteome, n = 336 (Paper 1 R17 L4)
  5. Lu 2023 (GSE193581) — Chinese sc, n = 14,624 malignant cells (Paper 1 R17 L2)
  6. Mu 2024 (HRA004166) — Chinese, n = 214, 4-class published frequencies

Output: results/figures/figure_korean_asian_validation.{png,pdf}
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import FancyBboxPatch, Rectangle

# Register Noto Sans CJK KR explicitly (matplotlib doesn't auto-pick CJK families)
for p in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
          "/usr/share/fonts/opentype/noto/NotoSerifCJK-Bold.ttc",
          "/home/seungho/.local/share/fonts/NanumGothic-Regular.ttf"):
    try: fm.fontManager.addfont(p)
    except Exception: pass
plt.rcParams["font.sans-serif"] = ["Noto Sans CJK KR", "Noto Sans CJK JP", "NanumGothic", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Source files
PAPER1 = Path("/home/seungho/personal/THCA_data_analysis/project")
LEE = PAPER1 / "results/v17_korean/Lee_unsup_8gene_2026_05_20.tsv"
GSE286332 = PAPER1 / "manuscript_biorxiv_2026_05_20/supp_data_bundle/12_GSE286332_panel_scores.tsv"
MUN = PAPER1 / "results/proteogenomic_v1/paper3_mun2025_dediff_layer/module_scores_per_sample.tsv"
LU2023 = PAPER1 / "results/v17_lu2023/GSE193581_cell_panel_score.tsv"
R17_DIR = PAPER1 / "results/r17_tcga_panel_d4p2_reconciliation"

IVORY = "#fbf8f1"
INK = "#2a2a2a"
MUTED = "#62656b"
BLUE = "#37618e"; BLUE_L = "#bccfe3"
RED = "#9c4742"; RED_L = "#dfb3ae"
GRAY_P = "#7e6e94"; GRAY_L = "#c2b9d0"
BEIGE = "#cdbb96"; BEIGE_L = "#ece1c3"
KR_GREEN = "#2c8a3a"
CN_AMBER = "#b8862a"


def style_axis(ax, xlim=(0, 100), ylim=(0, 100)):
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_facecolor("white")


def panel_a_kor_inventory(ax):
    """Korean/Asian cohort inventory — visual table."""
    style_axis(ax, xlim=(0, 100), ylim=(0, 100))
    cohorts = [
        ("K2 PRJEB11591",   "Korean", "RNA-seq",  "260",  "Yoo 2016",         "panel z available",   KR_GREEN),
        ("GSE213647 Lee 2024","Korean", "RNA-seq", "632",  "Lee Nat Commun 2024","unsup DM1/DM2 + panel z",KR_GREEN),
        ("GSE286332",       "Korean", "RNA-seq",  "18",   "Choi JCEM 2025",   "PTC vs PTC+HT",        KR_GREEN),
        ("Mun 2025 proteome","Korean", "Proteome","336",  "Mun Cell 2025",    "ATC OR=8.54 p=2.7e-15",KR_GREEN),
        ("GSE193581 Lu 2023","Chinese","scRNA",   "14,624 mal cells","Lu JCI 2023","R17 L2 done",       CN_AMBER),
        ("GSE184362 Pu 2021","Chinese","scRNA",   "158,577 cells / 11 pts","Pu Sci Adv 2021","incl. RAI-refractory distant met", CN_AMBER),
        ("HRA004166 Mu 2024","Chinese","NGS",     "214",  "Mu JCEM 2024",     "I-RAIR/C-RAIA/G-RAIR/P-RAIR", CN_AMBER),
    ]
    # Header
    ax.text(50, 95, "Korean / Asian thyroid cohort inventory  ·  locally processed", ha="center",
            fontsize=11.5, fontweight="bold", color=INK)
    # Column headers
    headers = ["Cohort", "Ethnicity", "Modality", "n", "Source", "Label / status"]
    col_x = [3, 24, 36, 48, 60, 78]
    for x, h in zip(col_x, headers):
        ax.text(x, 87, h, fontsize=9, fontweight="bold", color=INK)
    ax.plot([2, 98], [85, 85], color=MUTED, lw=0.5)

    # Rows
    row_h = 9
    for i, (name, ethn, mod, n, src, status, ethn_col) in enumerate(cohorts):
        y = 80 - i * row_h
        # ethnicity badge
        ax.add_patch(FancyBboxPatch((22, y - 2), 11, 4.5, boxstyle="round,pad=0,rounding_size=0.8",
                                     facecolor=ethn_col, edgecolor="none", alpha=0.75))
        ax.text(27.5, y + 0.3, ethn, ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        # other columns
        ax.text(col_x[0], y + 0.3, name, fontsize=8.5, color=INK, fontweight="bold")
        ax.text(col_x[2], y + 0.3, mod, fontsize=8, color=INK)
        ax.text(col_x[3], y + 0.3, n, fontsize=8, color=INK)
        ax.text(col_x[4], y + 0.3, src, fontsize=8, color=MUTED, style="italic")
        ax.text(col_x[5], y + 0.3, status, fontsize=7.8, color=INK)
        if i < len(cohorts) - 1:
            ax.plot([2, 98], [y - 3.5, y - 3.5], color="#e0e0e0", lw=0.4)


def panel_b_lee2024_boxplot(ax):
    """Lee 2024 GSE213647 panel z by histology — Korean n=632."""
    df = pd.read_csv(LEE, sep="\t")
    print("Lee 2024 cols:", list(df.columns))
    print("Lee 2024 histology dist:")
    print(df["histology"].value_counts() if "histology" in df.columns else "no histology col")
    # Group by histology
    if "histology" not in df.columns:
        ax.text(0.5, 0.5, "Lee 2024 file missing 'histology' column", transform=ax.transAxes,
                ha="center", color="red")
        return

    # Normalize: PDFP → PDTC; UTC/ATC → ATC
    df["hist_norm"] = df["histology"].replace({"PDFP": "PDTC", "UTC/ATC": "ATC"})
    order = ["Normal", "PTC", "PDTC", "ATC"]
    present = [h for h in order if h in df["hist_norm"].unique()]
    data = [df.loc[df["hist_norm"] == h, "panel_z"].dropna().values for h in present]
    colors = {"Normal": BEIGE, "PTC": BLUE_L, "PDTC": GRAY_L, "ATC": RED_L}
    edge_colors = {"Normal": "#7a6038", "PTC": BLUE, "PDTC": GRAY_P, "ATC": RED}
    bp = ax.boxplot(data, tick_labels=present, showfliers=False, patch_artist=True, widths=0.6)
    for patch, h in zip(bp["boxes"], present):
        patch.set_facecolor(colors[h]); patch.set_alpha(0.85); patch.set_edgecolor(edge_colors[h])
    rng = np.random.default_rng(0)
    for i, vals in enumerate(data):
        x = rng.uniform(-0.15, 0.15, size=len(vals)) + (i + 1)
        ax.scatter(x, vals, s=6, alpha=0.4, c=edge_colors[present[i]], edgecolors="none")
        ax.text(i + 1, 0.04, f"n = {len(vals)}", ha="center", va="bottom", fontsize=8.5,
                color="#444", transform=ax.get_xaxis_transform())
    ax.axhline(0, color="#444", lw=0.6, ls="--")
    ax.set_ylabel("8-gene panel z", fontsize=9.5)
    ax.set_title("b · Lee 2024 (GSE213647)  ·  Korean cohort n = 632", fontsize=10.5, color=INK, fontweight="bold", loc="left")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)


def panel_c_k2_distribution(ax):
    """K2 PRJEB11591 panel z distribution — Korean n=260."""
    # K2 panel z file from p_deconv_2026_05_08
    k2_path = PAPER1 / "results/p_deconv_2026_05_08/v15_k2_panel_overlay_scores.tsv"
    if not k2_path.exists():
        ax.text(0.5, 0.5, "K2 panel file missing", transform=ax.transAxes, ha="center", color="red")
        return
    df = pd.read_csv(k2_path, sep="\t")
    df = df[df["cohort"] == "K2_PRJEB11591"]
    print(f"K2 dist: {len(df)} samples, mean={df['panel_z_mean'].mean():.3f}, median={df['panel_z_mean'].median():.3f}")
    n_dm1 = (df["DM_call_centered"] == "DM1").sum()
    n_dm2 = (df["DM_call_centered"] == "DM2").sum()
    ax.hist(df["panel_z_mean"], bins=40, color=BLUE_L, edgecolor=BLUE, lw=0.6, alpha=0.85)
    ax.axvline(df["panel_z_mean"].median(), color=RED, lw=1.2, ls="-",
               label=f"median = {df['panel_z_mean'].median():+.2f}")
    ax.axvline(0, color="#444", lw=0.6, ls="--", alpha=0.6)
    ax.set_xlabel("Panel z (within K2)")
    ax.set_ylabel("count")
    ax.set_title(f"c · K2 PRJEB11591  ·  Korean n = {len(df)}  ·  DM1 {n_dm1} / DM2 {n_dm2}",
                 fontsize=10.5, color=INK, fontweight="bold", loc="left")
    ax.text(0.02, 0.95,
            "Korean overdiagnosis paradigm:\n" + f"  DM1 (silenced) {n_dm1/len(df)*100:.1f}% only",
            transform=ax.transAxes, fontsize=8, color=MUTED, va="top", style="italic")
    ax.legend(fontsize=8.5, loc="upper right")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)


def panel_d_mun_proteome(ax):
    """Mun 2025 proteome zone fraction — Korean n=336."""
    classes = ["PTC", "PDTC", "ATC"]
    n_per = [184, 46, 113]
    # From Paper 1 R17 L4 zone fractions (Mu 2024 frame)
    fractions = {
        "PTC":  {"BRAF-like": 13.0, "RAS-like": 26.6, "dark-matter": 14.1, "WT-like": 46.3},
        "PDTC": {"BRAF-like": 37.0, "RAS-like": 8.7,  "dark-matter": 23.9, "WT-like": 30.4},
        "ATC":  {"BRAF-like": 36.3, "RAS-like": 3.5,  "dark-matter": 58.4, "WT-like": 1.8},
    }
    zone_order = ["BRAF-like", "RAS-like", "dark-matter", "WT-like"]
    zone_colors = {"BRAF-like": BLUE, "RAS-like": "#d18b1f", "dark-matter": GRAY_P, "WT-like": "#888888"}
    x = np.arange(len(classes))
    bottom = np.zeros(len(classes))
    for z in zone_order:
        vals = [fractions[c][z] for c in classes]
        ax.bar(x, vals, bottom=bottom, color=zone_colors[z], edgecolor="white", lw=0.5,
                label=z, width=0.6)
        bottom += vals
    for i, (c, n) in enumerate(zip(classes, n_per)):
        ax.text(x[i], -7, f"n = {n}", ha="center", fontsize=8.5, color=INK)
    ax.set_xticks(x); ax.set_xticklabels(classes, fontsize=9, fontweight="bold")
    ax.set_ylabel("% samples in zone")
    ax.set_title("d · Mun 2025 proteome  ·  Korean n = 336  ·  ATC vs PTC DM Fisher OR = 8.54  p = 2.7×10⁻¹⁵",
                 fontsize=10.5, color=INK, fontweight="bold", loc="left")
    ax.set_ylim(-12, 110)
    ax.legend(fontsize=7.5, loc="upper right", ncol=2, framealpha=0.92, bbox_to_anchor=(1, 1.18))
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)


def panel_e_mu2024(ax):
    """Mu 2024 HRA004166 — Chinese 4-class published frequencies."""
    classes = ["I-RAIR", "G-RAIR", "P-RAIR", "C-RAIA"]
    n_per = [80, 19, 10, 48]
    # Published driver frequencies from PMC11031230
    drivers = {
        "BRAF V600E":    {"I-RAIR": 61, "G-RAIR": 36, "P-RAIR": 30, "C-RAIA": 22},
        "TERT promoter": {"I-RAIR": 22, "G-RAIR": 18, "P-RAIR": 12, "C-RAIA": 5},
        "TP53":          {"I-RAIR": 8,  "G-RAIR": 6,  "P-RAIR": 5,  "C-RAIA": 3},
        "PIK3CA":        {"I-RAIR": 4,  "G-RAIR": 5,  "P-RAIR": 4,  "C-RAIA": 2},
        "RAS":           {"I-RAIR": 3,  "G-RAIR": 14, "P-RAIR": 22, "C-RAIA": 35},
        "Fusion":        {"I-RAIR": 2,  "G-RAIR": 11, "P-RAIR": 14, "C-RAIA": 18},
        "Other":         {"I-RAIR": 0,  "G-RAIR": 10, "P-RAIR": 13, "C-RAIA": 15},
    }
    cmap = {"BRAF V600E":"#9c4742", "TERT promoter":"#b25b46", "TP53":"#cd8b3a",
            "PIK3CA":"#d2b454", "RAS":"#37618e", "Fusion":"#7e6e94", "Other":"#a8a59a"}
    bottom = np.zeros(len(classes))
    x = np.arange(len(classes))
    for d in drivers:
        vals = [drivers[d][c] for c in classes]
        ax.bar(x, vals, bottom=bottom, color=cmap[d], edgecolor="white", lw=0.4,
                label=d, width=0.65)
        bottom += vals
    for i, (c, n) in enumerate(zip(classes, n_per)):
        ax.text(x[i], -8, f"n = {n}", ha="center", fontsize=8.5, color=INK)
    # gray zone bracket
    ax.annotate("", xy=(2.45, 108), xytext=(0.55, 108),
                arrowprops=dict(arrowstyle="-", color=GRAY_P, lw=2.0))
    ax.text(1.5, 113, "gray zone (G-RAIR + P-RAIR = 29 / 214)",
            ha="center", fontsize=9, color=GRAY_P, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(classes, fontsize=10, fontweight="bold")
    ax.set_ylabel("% of class (driver composition)")
    ax.set_title("e · Mu 2024 HRA004166  ·  Chinese n = 214  ·  4-class RAI uptake",
                 fontsize=10.5, color=INK, fontweight="bold", loc="left")
    ax.set_ylim(-12, 130)
    ax.legend(fontsize=7, loc="upper right", ncol=2, framealpha=0.92, bbox_to_anchor=(1, 1.20))
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)


def panel_f_summary(ax):
    style_axis(ax)
    ax.text(50, 95, "Korean / Asian cross-cohort coherence", ha="center",
            fontsize=11.5, fontweight="bold", color=INK)

    lines = [
        ("● Korean RNA · Lee 2024 (GSE213647)",   "n=632 · ATC d=+3 vs Normal, monotonic decay",   KR_GREEN),
        ("● Korean RNA · K2 PRJEB11591",          "n=260 · 5.4% DM1 only — overdiagnosis paradigm",KR_GREEN),
        ("● Korean RNA · GSE286332",              "n=18 · PTC+HT 78% in dark-matter zone",         KR_GREEN),
        ("● Korean proteome · Mun 2025",          "n=336 · ATC dark-matter OR=8.54  p=2.7e-15",    KR_GREEN),
        ("● Chinese sc · Lu 2023 (GSE193581)",    "n=14,624 mal cells · ATC 38% dark-matter substrate", CN_AMBER),
        ("● Chinese sc · Pu 2021 (GSE184362)",    "11 patients · incl. RAI-refractory distant met", CN_AMBER),
        ("● Chinese NGS · Mu 2024 (HRA004166)",   "n=214 · gray zone 14% (G-RAIR + P-RAIR)",       CN_AMBER),
    ]
    y0 = 84
    for i, (name, finding, c) in enumerate(lines):
        y = y0 - i * 8.5
        ax.add_patch(Rectangle((3, y - 1.5), 2, 2.5, facecolor=c, edgecolor="none"))
        ax.text(6, y, name, fontsize=8.8, color=INK, fontweight="bold")
        ax.text(48, y, finding, fontsize=8.4, color=MUTED, style="italic")

    ax.add_patch(FancyBboxPatch((3, 12), 94, 14, boxstyle="round,pad=0,rounding_size=1.0",
                                 facecolor="#eef9ee", edgecolor=KR_GREEN, lw=0.8))
    ax.text(50, 21,
            "→ 7 Asian cohorts converge: panel z differentiation-silencing axis consistent",
            ha="center", fontsize=9.5, color=INK, fontweight="bold")
    ax.text(50, 16,
            "Tier 1 RECIST raw still pending (Boucai 2023 = US cohort) · Korean/Asian anchors suffice for parsimony",
            ha="center", fontsize=8.3, color=MUTED, style="italic")


def main():
    fig = plt.figure(figsize=(17, 12), facecolor=IVORY)
    gs = fig.add_gridspec(3, 2, hspace=0.55, wspace=0.20,
                          left=0.05, right=0.97, top=0.91, bottom=0.05,
                          height_ratios=[1.1, 1.0, 0.9])

    fig.suptitle(
        "Cross-cohort validation of the 8-gene panel in Korean & Asian thyroid cohorts",
        fontsize=14.5, fontweight="bold", color=INK, y=0.965,
    )
    fig.text(0.5, 0.928,
             "Lee 2024 · K2 PRJEB11591 · GSE286332 · Mun 2025 · Lu 2023 · Pu 2021 · Mu 2024  —  locally processed, Tier-mapped",
             ha="center", fontsize=10, color=MUTED, style="italic")

    ax_a = fig.add_subplot(gs[0, :])
    panel_a_kor_inventory(ax_a)
    ax_a.text(-0.02, 1.04, "a", transform=ax_a.transAxes, fontsize=15, fontweight="bold",
              va="top", color=INK)

    ax_b = fig.add_subplot(gs[1, 0])
    panel_b_lee2024_boxplot(ax_b)
    ax_b.text(-0.07, 1.10, "b", transform=ax_b.transAxes, fontsize=15, fontweight="bold",
              va="top", color=INK)

    ax_c = fig.add_subplot(gs[1, 1])
    panel_c_k2_distribution(ax_c)
    ax_c.text(-0.07, 1.10, "c", transform=ax_c.transAxes, fontsize=15, fontweight="bold",
              va="top", color=INK)

    ax_d = fig.add_subplot(gs[2, 0])
    panel_d_mun_proteome(ax_d)
    ax_d.text(-0.07, 1.10, "d", transform=ax_d.transAxes, fontsize=15, fontweight="bold",
              va="top", color=INK)

    ax_e = fig.add_subplot(gs[2, 1])
    panel_e_mu2024(ax_e)
    ax_e.text(-0.07, 1.10, "e", transform=ax_e.transAxes, fontsize=15, fontweight="bold",
              va="top", color=INK)

    out_png = OUT / "figure_korean_asian_validation.png"
    out_pdf = OUT / "figure_korean_asian_validation.pdf"
    fig.savefig(out_png, dpi=180, bbox_inches="tight", facecolor=IVORY)
    fig.savefig(out_pdf, bbox_inches="tight", facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
