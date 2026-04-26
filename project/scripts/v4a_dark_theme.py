"""
v4a_dark_theme.py - rewrite exported Plotly HTMLs to dark theme for v4a archive.

For each input .html in a manifest, locate the Plotly.newPlot(...) call, mutate the
layout JSON (template, paper_bgcolor, plot_bgcolor, font, colorway, gridcolor),
and write to reports/html/figs_interactive/v4a_dark/<basename>.
If layout JSON cannot be parsed, copy the file as-is and log a warning.

Visual system inspired by modern biotech platform conventions.
Not affiliated with any drug-discovery AI company.
"""
from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

DARK_COLORWAY = ["#F5A623", "#FFFFFF", "#C8C8C8", "#8A8A8A", "#5A5A5A", "#FFB84D", "#7A5BF0", "#2EB8B0"]
GRID_RGBA = "rgba(255,255,255,0.08)"
AXIS_LINE = "rgba(255,255,255,0.18)"
TICK_RGBA = "rgba(255,255,255,0.6)"
FONT_FAMILY = "Inter, Helvetica Neue, Arial, sans-serif"
FONT_COLOR = "#FFFFFF"


def _apply_dark(layout: dict) -> dict:
    if not isinstance(layout, dict):
        layout = {}
    layout["template"] = "plotly_dark"
    layout["paper_bgcolor"] = "rgba(0,0,0,0)"
    layout["plot_bgcolor"] = "rgba(0,0,0,0)"
    font = layout.get("font") or {}
    if not isinstance(font, dict):
        font = {}
    font["color"] = FONT_COLOR
    font["family"] = FONT_FAMILY
    layout["font"] = font
    layout["colorway"] = DARK_COLORWAY
    # axes
    for ax_key in list(layout.keys()):
        if ax_key.startswith(("xaxis", "yaxis")):
            ax = layout.get(ax_key) or {}
            if isinstance(ax, dict):
                ax["gridcolor"] = GRID_RGBA
                ax["linecolor"] = AXIS_LINE
                ax["zerolinecolor"] = GRID_RGBA
                ax["tickcolor"] = AXIS_LINE
                ax.setdefault("tickfont", {})
                if isinstance(ax["tickfont"], dict):
                    ax["tickfont"]["color"] = TICK_RGBA
                    ax["tickfont"]["family"] = FONT_FAMILY
                layout[ax_key] = ax
    # legend
    legend = layout.get("legend") or {}
    if isinstance(legend, dict):
        legend.setdefault("font", {})
        if isinstance(legend["font"], dict):
            legend["font"]["color"] = FONT_COLOR
            legend["font"]["family"] = FONT_FAMILY
        legend["bgcolor"] = "rgba(0,0,0,0)"
        layout["legend"] = legend
    return layout


NEWPLOT_RE = re.compile(r"Plotly\.newPlot\(\s*(?P<div>\"[^\"]+\"|'[^']+')\s*,\s*(?P<data>\[(?:.|\n)*?\])\s*,\s*(?P<layout>\{(?:.|\n)*?\})\s*(,\s*\{(?:.|\n)*?\}\s*)?\)")


def _balanced_slice(text: str, start: int, open_ch: str, close_ch: str) -> int:
    """Given that text[start] == open_ch, return index one past matching close_ch,
    respecting quotes and escapes. Returns -1 if unbalanced."""
    depth = 0
    i = start
    n = len(text)
    in_str = None
    while i < n:
        c = text[i]
        if in_str:
            if c == "\\":
                i += 2
                continue
            if c == in_str:
                in_str = None
        elif c in ('"', "'"):
            in_str = c
        elif c == open_ch:
            depth += 1
        elif c == close_ch:
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    return -1


def _extract_call_slices(text: str, start_idx: int):
    """Starting at 'Plotly.newPlot(', extract (data_str, layout_str, data_start, data_end, layout_start, layout_end)."""
    # find first '('
    p = text.find("(", start_idx)
    if p < 0:
        return None
    # skip first arg (div id string) — find comma at depth 0
    i = p + 1
    # skip whitespace
    while i < len(text) and text[i] in " \t\n\r":
        i += 1
    # the div id is a quoted string
    if text[i] in ('"', "'"):
        q = text[i]
        i += 1
        while i < len(text):
            if text[i] == "\\":
                i += 2
                continue
            if text[i] == q:
                i += 1
                break
            i += 1
    # skip comma + whitespace
    while i < len(text) and text[i] in " \t\n\r,":
        i += 1
    # data array
    if text[i] != "[":
        return None
    data_start = i
    data_end = _balanced_slice(text, i, "[", "]")
    if data_end < 0:
        return None
    j = data_end
    while j < len(text) and text[j] in " \t\n\r,":
        j += 1
    if j >= len(text) or text[j] != "{":
        return None
    layout_start = j
    layout_end = _balanced_slice(text, j, "{", "}")
    if layout_end < 0:
        return None
    return data_start, data_end, layout_start, layout_end


def rethemed_html(src_html: str) -> tuple[str, bool]:
    """Return (new_html, changed)."""
    # Find first Plotly.newPlot occurrence
    marker = "Plotly.newPlot"
    idx = src_html.find(marker)
    if idx < 0:
        return src_html, False
    slices = _extract_call_slices(src_html, idx)
    if not slices:
        return src_html, False
    data_s, data_e, lay_s, lay_e = slices
    layout_str = src_html[lay_s:lay_e]
    try:
        layout = json.loads(layout_str)
    except Exception:
        return src_html, False
    layout = _apply_dark(layout)
    try:
        new_layout_str = json.dumps(layout)
    except Exception:
        return src_html, False
    new_html = src_html[:lay_s] + new_layout_str + src_html[lay_e:]
    # Rewrite relative asset paths since we live one level deeper (figs_interactive/v4a_dark/)
    new_html = re.sub(r'"\.\./assets/', '"../../assets/', new_html)
    new_html = re.sub(r"'\.\./assets/", "'../../assets/", new_html)
    # also inject dark body background via style tag if present
    if "<body" in new_html and "background:#000" not in new_html:
        new_html = new_html.replace(
            "<body",
            '<style>html,body{background:#000;color:#fff;margin:0;padding:0;}</style><body',
            1,
        )
    else:
        new_html = (
            "<style>html,body{background:#000;color:#fff;margin:0;padding:0;}</style>"
            + new_html
        )
    return new_html, True


def retheme_file(src: Path, dst: Path) -> str:
    """Return 'rethemed' | 'copied' | 'failed'."""
    try:
        text = src.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return "failed"
    new_text, changed = rethemed_html(text)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if changed:
        dst.write_text(new_text, encoding="utf-8")
        return "rethemed"
    # fall-back: copy as-is so links resolve
    try:
        shutil.copy2(src, dst)
        return "copied"
    except Exception:
        return "failed"


def retheme_many(src_dir: Path, out_dir: Path, names: list[str]) -> dict:
    results = {}
    for name in names:
        src = src_dir / name
        if not src.exists():
            results[name] = "missing"
            continue
        dst = out_dir / name
        results[name] = retheme_file(src, dst)
    return results


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("usage: v4a_dark_theme.py <src_dir> <out_dir> [name1 name2 ...]")
        sys.exit(1)
    src = Path(sys.argv[1])
    out = Path(sys.argv[2])
    names = sys.argv[3:]
    if not names:
        names = [p.name for p in src.glob("*.html")]
    res = retheme_many(src, out, names)
    for k, v in res.items():
        print(f"{v:<12} {k}")
