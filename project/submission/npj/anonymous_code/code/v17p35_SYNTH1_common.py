"""v17p35 SYNTH-1 shared helpers: paths, palette, npj layout, save_figure."""
from __future__ import annotations
import os, sys, warnings
from pathlib import Path
warnings.filterwarnings('ignore', category=DeprecationWarning)

# Wire kaleido v1 to playwright's bundled Chromium so PNG/PDF export works without a system Chrome.
_PLAYWRIGHT_CHROME = Path.home() / '.cache' / 'ms-playwright' / 'chromium-1208' / 'chrome-linux64' / 'chrome'
if _PLAYWRIGHT_CHROME.exists():
    os.environ.setdefault('BROWSER_PATH', str(_PLAYWRIGHT_CHROME))

import plotly.io as pio  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PROJ = ROOT / 'project'
RES = PROJ / 'results'
META = PROJ / 'metadata'
SUB = PROJ / 'submission' / 'npj'
FIG_DIR = SUB / 'figures'
TBL_DIR = SUB / 'tables'
REPRO_DIR = SUB / 'reproducibility'
REPORTS = PROJ / 'reports' / 'v17p35'

for d in (FIG_DIR, TBL_DIR, REPRO_DIR, REPORTS):
    d.mkdir(parents=True, exist_ok=True)

DM1_COLOR = '#FF7F0E'
DM2_COLOR = '#1F77B4'
DM_PALETTE = {'DM1': DM1_COLOR, 'DM2': DM2_COLOR, 'DM1_like': DM1_COLOR, 'DM2_like': DM2_COLOR}
DRIVER_PALETTE = {'BRAF': '#D62728', 'RAS': '#2CA02C', 'unknown': '#7F7F7F', 'NTRK': '#9467BD', 'TP53': '#8C564B'}
HIST_PALETTE = {'cPTC': '#1F77B4', 'FVPTC': '#FF7F0E', 'PTC': '#1F77B4', 'PDTC': '#9467BD', 'ATC': '#D62728'}

NPJ_LAYOUT = dict(
    font=dict(family='Arial, Helvetica, sans-serif', size=11, color='#222'),
    title=dict(font=dict(size=13, color='#111')),
    plot_bgcolor='white',
    paper_bgcolor='white',
    margin=dict(l=55, r=20, t=55, b=50),
    legend=dict(font=dict(size=10), bgcolor='rgba(255,255,255,0.9)'),
)

NPJ_AXIS = dict(showline=True, linewidth=1, linecolor='#222', mirror=False, ticks='outside',
                tickfont=dict(size=10), title=dict(font=dict(size=11)),
                gridcolor='#EEE', zerolinecolor='#DDD')

def apply_npj(fig, width=None, height=None, title=None):
    fig.update_layout(**NPJ_LAYOUT)
    if title is not None:
        fig.update_layout(title_text=title, title_x=0.02)
    if width is not None: fig.update_layout(width=width)
    if height is not None: fig.update_layout(height=height)
    fig.update_xaxes(**NPJ_AXIS)
    fig.update_yaxes(**NPJ_AXIS)
    return fig


def save_figure(fig, name: str, width=1200, height=900, scale=3):
    """Save HTML(CDN, small) + PNG(600dpi-equivalent via scale=3) + PDF (vector)."""
    base = FIG_DIR / name
    fig.update_layout(width=width, height=height)
    fig.write_html(str(base.with_suffix('.html')), include_plotlyjs='cdn', full_html=True)
    fig.write_image(str(base.with_suffix('.png')), width=width, height=height, scale=scale)
    fig.write_image(str(base.with_suffix('.pdf')), width=width, height=height)
    sizes = {ext: os.path.getsize(str(base.with_suffix('.' + ext))) for ext in ('html', 'png', 'pdf')}
    return sizes


def log(msg):
    print(f'[SYNTH1] {msg}', flush=True)
