"""Shared Nature-style helpers — palette, axis styling, typography."""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Register CJK fonts (no Unifont; per weasyprint memory)
for p in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
          "/usr/share/fonts/opentype/noto/NotoSerifCJK-Bold.ttc"):
    try: fm.fontManager.addfont(p)
    except Exception: pass

# Refined Nature Medicine-grade palette
IVORY  = "#fbf8f1"
INK    = "#2a2a2a"
MUTED  = "#62656b"
GRID   = "#d9d4c5"

BLUE   = "#37618e"
BLUE_L = "#bccfe3"
RED    = "#9c4742"
RED_L  = "#dfb3ae"
GRAY_P = "#7e6e94"
GRAY_L = "#c2b9d0"
BEIGE  = "#cdbb96"
BEIGE_L = "#ece1c3"
AMBER  = "#b8862a"
KR_GREEN = "#2c8a3a"

ZONE_COLORS = {
    "BRAF-like":    BLUE,
    "RAS-like":     "#d18b1f",
    "dark-matter":  GRAY_P,
    "WT-like":      "#888888",
}
ZONE_FACES = {
    "BRAF-like":    BLUE_L,
    "RAS-like":     "#f0d4a8",
    "dark-matter":  GRAY_L,
    "WT-like":      "#dddddd",
}


def setup_rc():
    """Set global matplotlib rc params for a consistent Nature-style figure."""
    plt.rcParams.update({
        "font.family": ["Liberation Sans", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9.5,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.titlepad": 6,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK,
        "axes.linewidth": 0.7,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.fontsize": 8.5,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "#d9d4c5",
        "figure.facecolor": IVORY,
        "axes.facecolor": "white",
        "savefig.facecolor": IVORY,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
    })


def panel_letter(ax, text):
    """Nature-style panel letter — lowercase bold, top-left, above axes box."""
    ax.text(-0.10, 1.08, text, transform=ax.transAxes, fontsize=14.5, fontweight="bold",
            color=INK, va="top", ha="left")


def panel_title(ax, text):
    """Nature-style panel header — sans-serif bold."""
    ax.set_title(text, loc="left", fontsize=10.5, color=INK, fontweight="bold", pad=4)


def style_blank(ax, xlim=(0, 100), ylim=(0, 100)):
    """Blank schematic axis."""
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_facecolor(IVORY)
