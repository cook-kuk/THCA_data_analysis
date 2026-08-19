"""Unified NC-style v4 — refined Nature Cancer 2024-2025 aesthetic.

Color palette (Nature Cancer trim):
  Crimson (DM1)    = #B91C1C  (deep, less orange than #C0392B)
  Royal blue (DM2) = #1E40AF  (deeper, more saturated than #2C5282)
  Amber            = #B45309
  Sage             = #047857
  Slate            = #334155
  Cool gray (muted)= #94A3B8
  Panel bg         = #FFFFFF  (clean white)
  Page bg          = #F8FAFC  (very subtle gray)
  Grid             = #E2E8F0
  Ink              = #0F172A  (rich black)
Typography: Inter-like; bigger panel letters via inset circle badge.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from matplotlib import font_manager

try:
    font_manager.fontManager.addfont("/home/seungho/.local/share/fonts/NanumGothic-Regular.ttf")
except Exception:
    pass
FONT = "DejaVu Sans"
FAMILY = [FONT, "NanumGothic", "DejaVu Sans"]

PAL = {
    "DM1":         "#B91C1C",
    "DM2":         "#1E40AF",
    "highlight":   "#B45309",
    "positive":    "#047857",
    "neutral":     "#334155",
    "muted":       "#94A3B8",
    "panel_bg":    "#FFFFFF",
    "page_bg":     "#F8FAFC",
    "grid":        "#E2E8F0",
    "annot_bg":    "#FFFFFF",
    "ink":         "#0F172A",
    "subhead":     "#64748B",
    "purple":      "#6D28D9",
    "teal":        "#0E7490",
    "salmon":      "#DC2626",
    "rose":        "#BE185D",
    "kr_accent":   "#B45309",
}

def apply():
    plt.rcParams.update({
        "font.family": FAMILY,
        "font.size": 9,
        "axes.titlesize": 10.5,
        "axes.titleweight": "semibold",
        "axes.titlelocation": "left",
        "axes.titlepad": 16,
        "axes.titlecolor": PAL["ink"],
        "axes.labelsize": 9,
        "axes.labelcolor": PAL["ink"],
        "axes.labelpad": 6,
        "axes.edgecolor": PAL["neutral"],
        "axes.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.facecolor": PAL["panel_bg"],
        "axes.axisbelow": True,
        "axes.grid": True,
        "axes.grid.axis": "y",
        "grid.color": PAL["grid"],
        "grid.linewidth": 0.55,
        "grid.alpha": 0.9,
        "grid.linestyle": (0, (1.2, 2.0)),
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.color": PAL["ink"],
        "ytick.color": PAL["ink"],
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.major.pad": 4,
        "ytick.major.pad": 4,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.borderaxespad": 0.6,
        "figure.facecolor": PAL["page_bg"],
        "savefig.facecolor": PAL["page_bg"],
        "savefig.bbox": "tight",
        "savefig.dpi": 180,
        "patch.linewidth": 0.7,
    })

def annot_box(facecolor=None, edgecolor=None, pad=5):
    """Subtle annotation box (Nature 2024 — soft, low-contrast border)."""
    return dict(facecolor=facecolor or "#FFFFFFEE",
                alpha=0.97,
                edgecolor=edgecolor or "#CBD5E1",
                boxstyle=f"round,pad={pad/15:.3f}",
                linewidth=0.7)

def annot_emph(color=None):
    """High-emphasis annotation box (key stat, colored border)."""
    c = color or PAL["DM1"]
    return dict(facecolor="#FFFFFFF2",
                alpha=0.99,
                edgecolor=c,
                boxstyle="round,pad=0.36",
                linewidth=1.2)

def panel_label(ax, letter, x=-0.13, y=1.16, fs=18):
    """Bold panel label (Nature 2024 — slightly larger, more lifted)."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=fs, fontweight="bold", color=PAL["ink"],
            ha="left", va="bottom", family=FAMILY)

def panel_kr(ax, text, y=-0.32, x=0.5, fs=8.8, color=None):
    """Korean caption below panel — italicized accent, sized for readability when zoomed."""
    ax.text(x, y, text, transform=ax.transAxes, ha="center", va="top",
            fontsize=fs, color=color or PAL["kr_accent"], style="italic",
            fontweight="normal", family=FAMILY)

def suptitle(fig, label, kr=None):
    """English title + Korean subtitle with generous breathing room + accent bar."""
    # Accent stripe (thin amber bar) at very top — Nature 2024 hero look
    ax_strip = fig.add_axes([0.04, 0.992, 0.04, 0.004], frameon=False)
    ax_strip.set_xticks([]); ax_strip.set_yticks([])
    ax_strip.set_facecolor(PAL["DM1"])
    fig.suptitle(label, x=0.04, y=0.978, fontsize=15, ha="left", fontweight="bold",
                 color=PAL["ink"], family=FAMILY)
    if kr:
        fig.text(0.04, 0.945, kr, fontsize=12, ha="left", color=PAL["kr_accent"],
                 style="italic", fontweight="bold", family=FAMILY)
    # Soft horizontal rule below title block
    ax_rule = fig.add_axes([0.04, 0.925, 0.92, 0.0005], frameon=False)
    ax_rule.set_xticks([]); ax_rule.set_yticks([])
    ax_rule.axhline(0, color=PAL["grid"], lw=1.0)

def style_axes(ax):
    """Final polish: visible spines bottom/left, subtle grid behind data."""
    for s in ("bottom","left"):
        ax.spines[s].set_color(PAL["neutral"])
        ax.spines[s].set_linewidth(0.9)
    ax.set_axisbelow(True)

def panel_bg(ax, alpha=1.0):
    """Optional: solid white panel background to lift it off the page-bg."""
    ax.set_facecolor(PAL["panel_bg"])
