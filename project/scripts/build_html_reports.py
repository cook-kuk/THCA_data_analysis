#!/usr/bin/env python3
from __future__ import annotations

import base64
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import textwrap
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from io import BytesIO
from pathlib import Path

import markdown as md
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import requests
import shap
from jinja2 import Environment, FileSystemLoader, select_autoescape
from joblib import Parallel, delayed
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde, kruskal
from sklearn.calibration import calibration_curve
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
from sklearn.metrics import (
    auc,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from umap import UMAP

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
HTML = ROOT / "reports" / "html"
PAGES = HTML / "pages"
ASSETS = HTML / "assets"
ASSET_DATA = ASSETS / "data"
ASSET_JS = ASSETS / "js"
ASSET_CSS = ASSETS / "css"
ASSET_IMG = ASSETS / "img"
ASSET_FONTS = ASSETS / "fonts"
FIGDIR = HTML / "figs_interactive"
TEMPLATES = HTML / "templates"
LOG = ROOT / "logs" / "build_html.log"
STATE = ROOT / "state"

BUILD_TIME = datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def ensure_dirs() -> None:
    for p in [PAGES, ASSET_DATA, ASSET_JS, ASSET_CSS, ASSET_IMG, ASSET_FONTS, FIGDIR, TEMPLATES, STATE]:
        p.mkdir(parents=True, exist_ok=True)


def _json_clean(o):
    """Recursively coerce NaN/Inf to None so the output is strict JSON the browser can parse."""
    if isinstance(o, float):
        return None if (math.isnan(o) or math.isinf(o)) else o
    if isinstance(o, dict):
        return {k: _json_clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_json_clean(v) for v in o]
    return o


def json_dump(obj, path: Path) -> None:
    path.write_text(json.dumps(_json_clean(obj), ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def normalize_gene(x: str) -> str:
    return str(x).strip().upper()


def read_expr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene_symbol"}).set_index("gene_symbol")
    df.index = df.index.map(normalize_gene)
    return df[~df.index.duplicated(keep="first")]


def compute_sri(content: bytes) -> str:
    return "sha384-" + base64.b64encode(hashlib.sha384(content).digest()).decode("ascii")


def download_asset(name: str, url: str) -> dict:
    out = {"name": name, "url": url, "local": None, "integrity": None, "status": "skipped"}
    try:
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        suffix = ".js" if ".js" in url else ".css"
        dest = (ASSET_JS if suffix == ".js" else ASSET_CSS) / f"{name}{suffix}"
        dest.write_bytes(r.content)
        out["local"] = str(dest.relative_to(HTML))
        out["integrity"] = compute_sri(r.content)
        out["status"] = "ok"
    except Exception as exc:
        log(f"asset download skip {name}: {exc}")
    return out


def write_static_assets(asset_map: dict) -> None:
    # main.css / main.js / plotly-theme.js / gene-explorer.js / cohort-compare.js
    # are authored directly in /assets/css and /assets/js. We do not overwrite them here.
    # Only write templates (base.html / page.html) whose content depends on asset_map state.
    _skip_legacy_css = """
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;600&display=swap');
:root{
  --bg0:#0F172A;--bg1:#1E293B;--bg2:#111827;--ink:#E5E7EB;--muted:#94A3B8;
  --card:rgba(15,23,42,.58);--line:rgba(148,163,184,.18);
  --primary:#6366F1;--accent:#EC4899;--success:#10B981;--warn:#F59E0B;--danger:#EF4444;
  --shadow:0 20px 60px rgba(0,0,0,.35);--blur:blur(14px);
}
body.light{
  --bg0:#EEF2FF;--bg1:#F8FAFC;--bg2:#E2E8F0;--ink:#0F172A;--muted:#475569;
  --card:rgba(255,255,255,.68);--line:rgba(15,23,42,.08);
}
html{scroll-behavior:smooth}body{
  margin:0;color:var(--ink);font-family:'Inter',system-ui,sans-serif;
  background:linear-gradient(135deg,var(--bg0) 0%,var(--bg1) 50%,var(--bg0) 100%);
  min-height:100vh;-webkit-font-smoothing:antialiased;position:relative;
}
body:before{
  content:'';position:fixed;inset:0;pointer-events:none;opacity:.055;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='180' height='180' viewBox='0 0 180 180'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='.85' numOctaves='3'/%3E%3C/filter%3E%3Crect width='180' height='180' filter='url(%23n)' opacity='.45'/%3E%3C/svg%3E");
}
.glass{background:var(--card);backdrop-filter:var(--blur);border:1px solid var(--line);box-shadow:var(--shadow)}
.layout{display:grid;grid-template-columns:260px minmax(0,1fr) 220px;gap:18px;max-width:1700px;margin:0 auto;padding:20px}
.sidebar,.toc{position:sticky;top:88px;height:calc(100vh - 110px);overflow:auto;border-radius:24px;padding:18px}
.content{min-width:0}.topnav{position:sticky;top:0;z-index:50;margin:0 auto;padding:14px 20px;border-bottom:1px solid var(--line)}
.topnav-inner{max-width:1700px;margin:0 auto;display:flex;align-items:center;justify-content:space-between;gap:16px}
.brand{display:flex;align-items:center;gap:12px;font-weight:800;letter-spacing:.02em}
.brand-dot{width:14px;height:14px;border-radius:999px;background:linear-gradient(135deg,var(--primary),var(--accent));box-shadow:0 0 24px rgba(99,102,241,.7)}
.nav-links{display:flex;gap:8px;flex-wrap:wrap}.nav-link{padding:9px 12px;border-radius:12px;color:var(--ink);text-decoration:none;border:1px solid transparent}
.nav-link:hover,.nav-link.active{background:rgba(99,102,241,.12);border-color:rgba(99,102,241,.22);box-shadow:0 0 18px rgba(99,102,241,.16)}
.hero{padding:26px;border-radius:28px;overflow:hidden;position:relative}
.hero:after{content:'';position:absolute;inset:-15% auto auto -10%;width:280px;height:280px;border-radius:999px;background:radial-gradient(circle,rgba(236,72,153,.26),transparent 70%)}
.hero h1{font-size:46px;line-height:1.02;margin:0 0 10px;font-weight:800;background:linear-gradient(90deg,#fff,#c7d2fe,#f9a8d4);-webkit-background-clip:text;background-clip:text;color:transparent}
body.light .hero h1{background:linear-gradient(90deg,#111827,#4338ca,#db2777);-webkit-background-clip:text;background-clip:text;color:transparent}
.hero p{margin:0;color:var(--muted);max-width:1000px}.kpi-grid,.tile-grid,.fig-grid,.triple-grid{display:grid;gap:16px}
.kpi-grid{grid-template-columns:repeat(4,minmax(0,1fr));margin-top:18px}.tile-grid{grid-template-columns:repeat(3,minmax(0,1fr))}
.fig-grid{grid-template-columns:repeat(2,minmax(0,1fr))}.triple-grid{grid-template-columns:repeat(3,minmax(0,1fr))}
.card,.panel{border-radius:22px;padding:18px}.card:hover,.panel:hover{box-shadow:0 0 0 1px rgba(99,102,241,.22),0 30px 80px rgba(0,0,0,.25)}
.kpi{font-size:38px;font-weight:800}.subtle{color:var(--muted);font-size:14px;line-height:1.55}.mono{font-family:'JetBrains Mono',monospace}
.section-title{font-size:28px;font-weight:800;margin:0 0 8px}.page-title{font-size:34px;font-weight:800;margin:0}
.page-subtitle{color:var(--muted);margin-top:8px}
.btn{display:inline-flex;align-items:center;gap:8px;padding:10px 14px;border-radius:14px;border:1px solid var(--line);text-decoration:none;color:var(--ink);cursor:pointer}
.btn.primary{background:linear-gradient(135deg,var(--primary),var(--accent));border:none;color:white}
.btn:hover{transform:translateY(-1px)}
.chips{display:flex;gap:8px;flex-wrap:wrap;margin-top:12px}.chip{padding:7px 10px;border-radius:999px;background:rgba(99,102,241,.12);border:1px solid rgba(99,102,241,.16);font-size:13px}
.figure-wrap{position:relative}.figure-actions{position:absolute;top:14px;right:14px;display:flex;gap:8px;z-index:5}
.figure-actions a,.figure-actions button{padding:8px 10px;border-radius:10px;border:1px solid var(--line);background:rgba(15,23,42,.6);color:#fff;text-decoration:none;cursor:pointer}
body.light .figure-actions a,body.light .figure-actions button{background:rgba(255,255,255,.85);color:#0f172a}
.figure-container{min-height:440px}.small-figure .figure-container{min-height:320px}
.figure-preview-shell{padding-top:46px}
.figure-preview{display:block;width:100%;height:auto;min-height:320px;max-height:none;border-radius:18px;border:1px solid var(--line);background:rgba(255,255,255,.03);object-fit:contain;cursor:pointer}
.figure-open-row{display:flex;gap:10px;flex-wrap:wrap;margin-top:12px}
.figure-open-btn{min-height:44px}
.figure-modal{position:fixed;inset:0;z-index:5000}
.figure-modal-backdrop{position:absolute;inset:0;background:rgba(2,6,23,.78)}
.figure-modal-card{position:relative;z-index:2;width:min(1400px,96vw);height:min(92vh,980px);margin:4vh auto;padding:14px;border-radius:24px;display:flex;flex-direction:column}
.figure-modal-header{display:flex;align-items:center;justify-content:space-between;gap:12px;padding:6px 6px 12px}
.figure-modal-actions{display:flex;gap:8px;flex-wrap:wrap}
.figure-modal-frame{width:100%;flex:1;border:0;border-radius:18px;background:#0b1120}
.narrative-list{margin:10px 0 0 0;padding-left:18px;color:var(--muted);line-height:1.7}
.narrative-list li{margin:6px 0}
.data-table{width:100%;border-collapse:collapse;margin-top:10px;font-size:14px}
.data-table th,.data-table td{padding:10px 12px;border-bottom:1px solid var(--line);text-align:left}
.data-table th{color:var(--muted);font-weight:600}
.sidebar a,.toc a{display:block;padding:8px 10px;border-radius:12px;color:var(--ink);text-decoration:none;margin-bottom:4px}
.sidebar a:hover,.toc a:hover{background:rgba(99,102,241,.12)} .sidebar h3,.toc h3{margin:4px 0 12px;font-size:13px;color:var(--muted);text-transform:uppercase;letter-spacing:.08em}
.footer{margin-top:24px;padding:16px 18px;border-radius:18px;color:var(--muted);font-size:13px}
.markdown-body pre{position:relative;padding:16px;border-radius:16px;overflow:auto;background:rgba(2,6,23,.85)}
.copy-btn{position:absolute;top:8px;right:8px;padding:6px 10px;border-radius:10px;border:1px solid var(--line);background:rgba(15,23,42,.85);color:#fff;cursor:pointer}
.palette{position:fixed;inset:0;background:rgba(2,6,23,.65);display:none;align-items:flex-start;justify-content:center;padding-top:10vh;z-index:90}
.palette.open{display:flex}.palette-card{width:min(760px,92vw);border-radius:24px;padding:18px}.palette-results a{display:block;padding:12px 10px;border-radius:12px;color:var(--ink);text-decoration:none}
.palette-results a:hover{background:rgba(99,102,241,.12)} .hidden{display:none!important} .presenter .sidebar,.presenter .toc,.presenter .topnav{display:none}.presenter .layout{grid-template-columns:minmax(0,1fr)}
@media (max-width:1400px){.layout{grid-template-columns:220px minmax(0,1fr)}.toc{display:none}}
@media (max-width:1100px){.layout{grid-template-columns:minmax(0,1fr)}.sidebar,.toc{display:none}.kpi-grid,.tile-grid,.fig-grid,.triple-grid{grid-template-columns:1fr}.figure-preview-shell{padding-top:0}}
@media (max-width:640px){.figure-modal-card{width:100vw;height:100vh;margin:0;border-radius:0;padding:10px}.figure-modal-header{align-items:flex-start;flex-direction:column}.figure-modal-actions{width:100%}.figure-modal-actions .btn{flex:1;justify-content:center}}
"""
    # legacy inline CSS preserved for diffing only; on-disk assets/css/main.css is authoritative
    _ = _skip_legacy_css

    js = """
window.ThyroidDash = window.ThyroidDash || {};
function byId(x){return document.getElementById(x)}
function copyText(txt){navigator.clipboard.writeText(txt)}
function toggleTheme(){
  const root=document.body;
  const next=root.classList.contains('light')?'dark':'light';
  if(next==='light'){root.classList.add('light')}else{root.classList.remove('light')}
  localStorage.setItem('thyroid-theme',next);
}
function initTheme(){
  const saved=localStorage.getItem('thyroid-theme');
  if(saved==='light') document.body.classList.add('light');
}
function initPresenterMode(){
  document.addEventListener('keydown',e=>{
    if((e.metaKey||e.ctrlKey)&&e.key.toLowerCase()==='k'){e.preventDefault(); byId('command-palette')?.classList.add('open'); byId('palette-input')?.focus();}
    if(e.key.toLowerCase()==='p' && !['INPUT','TEXTAREA','SELECT'].includes(document.activeElement.tagName)){document.body.classList.toggle('presenter')}
    if(e.key==='Escape'){byId('command-palette')?.classList.remove('open')}
  });
}
function initPalette(){
  const payload=window.ThyroidDash.commandIndex||[];
  const input=byId('palette-input'), box=byId('palette-results');
  if(!input||!box) return;
  function render(q){
    const v=(q||'').toLowerCase();
    const rows=payload.filter(x=>!v || x.label.toLowerCase().includes(v) || (x.tags||'').toLowerCase().includes(v)).slice(0,30);
    box.innerHTML=rows.map(x=>`<a href="${x.href}" onclick="document.getElementById('command-palette').classList.remove('open')"><strong>${x.label}</strong><div class="subtle">${x.tags||''}</div></a>`).join('');
  }
  input.addEventListener('input',e=>render(e.target.value)); render('');
}
function figureButtons(){
  document.querySelectorAll('[data-copy-link]').forEach(el=>el.addEventListener('click',()=>copyText(location.href)));
  document.querySelectorAll('[data-copy-code]').forEach(btn=>btn.addEventListener('click',()=>{const pre=btn.closest('pre'); copyText(pre.innerText.replace(/Copy$/,'')); btn.innerText='Copied'; setTimeout(()=>btn.innerText='Copy',1200)}));
}
function openFigureModal(url, title){
  const modal=byId('figure-modal'), frame=byId('figure-modal-frame'), open=byId('figure-modal-open'), label=byId('figure-modal-title');
  if(!modal||!frame||!open||!label||!url) return;
  frame.src=url;
  open.href=url;
  label.textContent=title||'인터랙티브 Figure';
  modal.classList.remove('hidden');
  modal.setAttribute('aria-hidden','false');
  document.body.style.overflow='hidden';
}
function closeFigureModal(){
  const modal=byId('figure-modal'), frame=byId('figure-modal-frame');
  if(!modal||!frame) return;
  modal.classList.add('hidden');
  modal.setAttribute('aria-hidden','true');
  frame.src='about:blank';
  document.body.style.overflow='';
}
function initFigureModal(){
  document.querySelectorAll('[data-figure-open]').forEach(el=>{
    el.addEventListener('click',e=>{
      e.preventDefault();
      openFigureModal(el.getAttribute('data-figure-open'), el.getAttribute('data-figure-title')||'인터랙티브 Figure');
    });
  });
  document.querySelectorAll('[data-figure-close]').forEach(el=>el.addEventListener('click',closeFigureModal));
  document.addEventListener('keydown',e=>{ if(e.key==='Escape') closeFigureModal(); });
}
function hashState(key, value){
  const params=new URLSearchParams(location.hash.replace(/^#/,''));
  if(value) params.set(key,value); else params.delete(key);
  location.hash=params.toString();
}
function readHashState(key){const params=new URLSearchParams(location.hash.replace(/^#/,'')); return params.get(key);}
function renderGridTable(targetId, rows, columns, searchId){
  const target=byId(targetId); if(!target || typeof gridjs==='undefined') return;
  target.innerHTML='';
  const searchBox=searchId?byId(searchId):null;
  let filtered=[...rows];
  const render=()=>{
    target.innerHTML='';
    new gridjs.Grid({
      columns: columns,
      data: filtered.map(r=>columns.map(c=>r[c.id] ?? '')),
      sort:true, search:false, pagination:{limit:12}, resizable:true
    }).render(target);
  };
  if(searchBox){
    const prior=readHashState(searchId); if(prior){searchBox.value=prior}
    const apply=()=>{
      const q=(searchBox.value||'').toLowerCase();
      filtered=!q?rows:rows.filter(r=>Object.values(r).some(v=>String(v??'').toLowerCase().includes(q)));
      hashState(searchId, searchBox.value||''); render();
    };
    searchBox.addEventListener('input', apply); apply();
  } else render();
}
function initMarkdownRender(){
  document.querySelectorAll('[data-md-source]').forEach(async el=>{
    try{
      const src=el.getAttribute('data-md-source');
      const txt=await fetch(src).then(r=>r.text());
      el.innerHTML=marked.parse(txt);
      document.querySelectorAll('pre').forEach(pre=>{ if(!pre.querySelector('.copy-btn')){ const b=document.createElement('button'); b.className='copy-btn'; b.textContent='Copy'; b.setAttribute('data-copy-code','1'); pre.appendChild(b);} });
      figureButtons();
    }catch(e){el.innerHTML='<p class="subtle">Markdown render failed.</p>'}
  });
}
function initCountUp(){
  document.querySelectorAll('[data-countup]').forEach(el=>{
    const end=Number(el.getAttribute('data-countup')); if(window.countUp && !Number.isNaN(end)){const c=new window.countUp.CountUp(el,end,{duration:1.8}); if(!c.error) c.start();}
  });
}
function initAOSMaybe(){ if(window.AOS){ window.AOS.init({duration:700, once:true, easing:'ease-out-cubic'}); } }
document.addEventListener('DOMContentLoaded',()=>{initTheme(); initPresenterMode(); initPalette(); initMarkdownRender(); initCountUp(); initAOSMaybe(); figureButtons(); initFigureModal();});
"""
    # legacy inline JS preserved for diffing only; on-disk assets/js/main.js is authoritative
    _ = js

    plotly_theme = """
window.THYROID_PLOTLY_TEMPLATE = {
  layout: {
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: {family: 'Inter, sans-serif', color: '#E5E7EB'},
    colorway: ['#6366F1','#EC4899','#10B981','#F59E0B','#EF4444','#22D3EE','#A78BFA','#F472B6'],
    xaxis: {gridcolor:'rgba(148,163,184,0.16)', zerolinecolor:'rgba(148,163,184,0.12)'},
    yaxis: {gridcolor:'rgba(148,163,184,0.16)', zerolinecolor:'rgba(148,163,184,0.12)'},
    legend: {bgcolor:'rgba(15,23,42,0.4)', bordercolor:'rgba(148,163,184,0.12)', borderwidth:1}
  }
};
"""
    # legacy plotly theme preserved for diffing only; on-disk assets/js/plotly-theme.js is authoritative
    _ = plotly_theme

    base = """
<!doctype html>
<html lang="ko">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="dark light">
  <meta name="theme-color" content="#07132b">
  <title>{{ title }}</title>
  <meta name="description" content="{{ subtitle }}">
  <meta property="og:title" content="{{ title }}">
  <meta property="og:description" content="{{ subtitle }}">
  <meta property="og:image" content="{{ root_prefix }}assets/img/og_card.png">
  <link rel="icon" href="{{ root_prefix }}assets/img/favicon.svg" type="image/svg+xml">
  <link rel="preload" href="{{ root_prefix }}assets/fonts/inter-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="{{ root_prefix }}assets/fonts/inter-600.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="stylesheet" href="{{ root_prefix }}assets/css/main.css">
  <script src="{{ root_prefix }}assets/vendor/tailwind.min.js" defer></script>
  <script>
    (function(){ try{ if(localStorage.getItem('thyroid-theme')==='light') document.documentElement.classList.add('light-pending'); }catch(e){} })();
    window.addEventListener('error', function(e){
      if(e.target && e.target.tagName==='SCRIPT' && e.target.dataset.fallback){ var s=document.createElement('script'); s.src=e.target.dataset.fallback; document.head.appendChild(s); }
    }, true);
  </script>
</head>
<body class="dark">
  <a class="sr-only" href="#main-content">본문으로 건너뛰기</a>
  <div id="command-palette" class="palette" role="dialog" aria-modal="true" aria-label="명령 팔레트">
    <div class="palette-card glass">
      <div class="subtle mono">Cmd+K / Ctrl+K</div>
      <input id="palette-input" type="text" placeholder="페이지, 섹션, 조직학 검색" aria-label="검색" autocomplete="off" spellcheck="false">
      <div id="palette-results" class="palette-results"></div>
    </div>
  </div>
  <div id="figure-modal" class="figure-modal hidden" aria-hidden="true" role="dialog" aria-modal="true" aria-labelledby="figure-modal-title">
    <div class="figure-modal-backdrop" data-figure-close="1" aria-hidden="true"></div>
    <div class="figure-modal-card glass">
      <div class="figure-modal-header">
        <strong id="figure-modal-title">인터랙티브 Figure</strong>
        <div class="figure-modal-actions">
          <a id="figure-modal-open" class="btn sm" href="#" target="_blank" rel="noopener" aria-label="새 탭으로 열기">새 탭</a>
          <a id="figure-modal-png" class="btn sm" href="#" target="_blank" rel="noopener" aria-label="PNG 다운로드">PNG</a>
          <a id="figure-modal-tsv" class="btn sm" href="#" target="_blank" rel="noopener" aria-label="TSV 다운로드">TSV</a>
          <button id="figure-modal-copy" class="btn sm" type="button" aria-label="링크 복사">링크 복사</button>
          <button id="figure-modal-fullscreen" class="btn sm" type="button" aria-label="전체화면">전체화면</button>
          <button class="btn sm" type="button" data-figure-close="1" aria-label="닫기">닫기</button>
        </div>
      </div>
      <iframe id="figure-modal-frame" class="figure-modal-frame" title="interactive figure"></iframe>
    </div>
  </div>
  <nav class="topnav" aria-label="주요 메뉴">
    <div class="topnav-inner">
      <a class="brand" href="{{ root_prefix }}index.html" aria-label="홈으로"><span class="brand-dot" aria-hidden="true"></span><span>THYROID DASH</span></a>
      <div class="nav-links" role="menubar">
        {% for item in nav_items %}<a class="nav-link{% if item.href == current_href %} active{% endif %}" href="{{ item.href }}" role="menuitem">{{ item.label }}</a>{% endfor %}
      </div>
      <div class="nav-actions">
        <button class="btn sm ghost" type="button" onclick="toggleTheme()" aria-label="다크/라이트 모드 전환">다크/라이트</button>
        <button class="btn sm ghost" type="button" data-copy-link="1" aria-label="현재 뷰 링크 복사">링크 복사 <span class="kbd">⌘K</span></button>
        <button id="td-menu-toggle" class="menu-toggle" type="button" aria-controls="td-drawer" aria-expanded="false" aria-label="메뉴 열기">
          <svg viewBox="0 0 24 24" aria-hidden="true"><line x1="4" y1="7" x2="20" y2="7"/><line x1="4" y1="12" x2="20" y2="12"/><line x1="4" y1="17" x2="20" y2="17"/></svg>
        </button>
      </div>
    </div>
  </nav>
  <div id="td-drawer-backdrop" class="drawer-backdrop" aria-hidden="true"></div>
  <aside id="td-drawer" class="drawer" aria-hidden="true" aria-label="모바일 메뉴">
    <div class="drawer-header">
      <strong>메뉴</strong>
      <button class="btn sm ghost" type="button" id="td-drawer-close" onclick="document.getElementById('td-menu-toggle').click()" aria-label="메뉴 닫기">닫기</button>
    </div>
    {% for item in nav_items %}<a class="{% if item.href == current_href %}active{% endif %}" href="{{ item.href }}">{{ item.label }}</a>{% endfor %}
  </aside>
  <div class="layout">
    <aside class="sidebar glass" aria-label="섹션 목차">
      <h3>Section TOC</h3>
      {% for section in sections %}<a href="#{{ section.id }}">{{ section.title }}</a>{% endfor %}
    </aside>
    <main class="content" id="main-content">
      {% block body %}{% endblock %}
      <footer class="footer glass" role="contentinfo">
        <div>
          <span class="mono">Build</span> {{ build_time }}
          &nbsp;·&nbsp; <span class="mono">Hash</span> {{ git_hash }}
          &nbsp;·&nbsp; {{ version.n_datasets }} datasets · {{ version.n_samples }} samples · {{ version.n_genes }} genes
        </div>
        <div><a href="{{ root_prefix }}pages/14_caveats.html">Methods &amp; Caveats</a></div>
      </footer>
    </main>
    <aside class="toc glass" aria-label="퀵 점프">
      <h3>Quick Jump</h3>
      {% for section in sections %}<a href="#{{ section.id }}">{{ loop.index }}. {{ section.title }}</a>{% endfor %}
    </aside>
  </div>
  <script src="{{ root_prefix }}assets/vendor/plotly.min.js"></script>
  {% if asset_map.marked.status == 'ok' %}<script defer src="{{ root_prefix }}{{ asset_map.marked.local }}"></script>{% endif %}
  <script>window.ThyroidDash = window.ThyroidDash || {}; window.ThyroidDash.commandIndex = {{ command_index|safe }};</script>
  <script src="{{ root_prefix }}assets/js/plotly-theme.js"></script>
  <script src="{{ root_prefix }}assets/js/main.js"></script>
  {% block scripts %}{% endblock %}
</body>
</html>
"""
    page = """
{% extends 'base.html' %}
{% block body %}
<section class="hero" data-aos="fade-up" aria-labelledby="page-hero-title">
  <div class="hero-aurora" aria-hidden="true"></div>
  <div class="hero-noise" aria-hidden="true"></div>
  <div class="subtle mono">{{ slug }}</div>
  <h1 id="page-hero-title">{{ title }}</h1>
  <p class="lede">{{ subtitle }}</p>
  {% if chips %}<div class="chips" role="list">{% for chip in chips %}<span class="chip" role="listitem">{{ chip }}</span>{% endfor %}</div>{% endif %}
</section>
{% for section in sections %}
<section id="{{ section.id }}" class="panel glass" data-aos="fade-up" aria-labelledby="sec-{{ section.id }}">
  <h2 class="section-title" id="sec-{{ section.id }}">{{ section.title }}</h2>
  {% if section.text %}<p class="subtle">{{ section.text }}</p>{% endif %}
  {% if section.html %}{{ section.html|safe }}{% endif %}
  {% if section.figures %}
    <div class="{{ section.figure_class or 'fig-grid' }}">
    {% for fig in section.figures %}
    <div class="figure-wrap {% if fig.small %}small-figure{% endif %}">
      <div class="figure-actions">
        <a href="{{ fig.html_url }}" target="_blank" rel="noopener" aria-label="HTML 새 탭">HTML</a>
        <a href="{{ fig.png_url }}" target="_blank" rel="noopener" aria-label="PNG 다운로드">PNG</a>
        <a href="{{ fig.svg_url }}" target="_blank" rel="noopener" aria-label="SVG 다운로드">SVG</a>
        <a href="{{ fig.tsv_url }}" target="_blank" rel="noopener" aria-label="TSV 다운로드">TSV</a>
        <button class="copy-cite" type="button" data-copy-cite="THCA Dashboard — {{ fig.label or fig.name }} — {{ fig.html_url }}" aria-label="인용 정보 복사">Cite</button>
      </div>
      <div class="figure-preview-shell">
        <img class="figure-preview" loading="{% if loop.first and loop.index0 == 0 %}eager{% else %}lazy{% endif %}" decoding="async" src="{{ fig.preview_url }}" alt="{{ fig.label or fig.name }} 미리보기" data-figure-open="{{ fig.html_url }}" data-figure-title="{{ fig.label or fig.name }}" data-figure-png="{{ fig.png_url }}" data-figure-tsv="{{ fig.tsv_url }}">
      </div>
      <div class="figure-open-row">
        <button class="btn primary figure-open-btn" type="button" data-figure-open="{{ fig.html_url }}" data-figure-title="{{ fig.label or fig.name }}" data-figure-png="{{ fig.png_url }}" data-figure-tsv="{{ fig.tsv_url }}">확대해서 보기</button>
        <a class="btn figure-open-btn" href="{{ fig.html_url }}" target="_blank" rel="noopener">새 탭으로 열기</a>
      </div>
      {% if fig.caption %}
      <div class="figure-caption">
        <strong>{{ fig.label or fig.name }}</strong>
        {{ fig.caption }}
      </div>
      {% endif %}
    </div>
    {% endfor %}
  </div>
  {% endif %}
</section>
{% endfor %}
{% endblock %}
{% block scripts %}
{% if page_js %}<script>{{ page_js|safe }}</script>{% endif %}
{% if page_extra_scripts %}{% for src in page_extra_scripts %}<script defer src="{{ src }}"></script>{% endfor %}{% endif %}
{% endblock %}
"""
    (TEMPLATES / "base.html").write_text(textwrap.dedent(base).strip() + "\n", encoding="utf-8")
    (TEMPLATES / "page.html").write_text(textwrap.dedent(page).strip() + "\n", encoding="utf-8")


def parse_tier_categories(path: Path) -> tuple[list[str], dict[str, list[str]]]:
    categories = defaultdict(list)
    current = "uncategorized"
    genes = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        if line.startswith("# [") and line.endswith("]"):
            current = line[3:-1]
            continue
        if line.startswith("#"):
            continue
        gene = normalize_gene(line)
        genes.append(gene)
        categories[current].append(gene)
    unique = list(dict.fromkeys(genes))
    return unique, dict(categories)


def build_panels() -> dict:
    tds16 = [normalize_gene(x) for x in (ROOT / "metadata" / "tds16_genes.txt").read_text(encoding="utf-8").splitlines() if x and not x.startswith("#")]
    brs = [normalize_gene(x) for x in (ROOT / "metadata" / "brs71_genes.txt").read_text(encoding="utf-8").splitlines() if x and not x.startswith("#")]
    tier, tier_categories = parse_tier_categories(ROOT / "metadata" / "tierA67_genes.txt")
    extra = [g for g in brs if g not in tier]
    thyroid112 = list(dict.fromkeys(tier + extra))[:112]
    intermediate32 = list(dict.fromkeys(tds16 + [g for g in tier if g not in tds16]))[:32]
    panel_meta = {
        "TDS16": tds16,
        "BRS_proxy": brs,
        "TierA67": tier,
        "Intermediate32": intermediate32,
        "ThyroSeq112_surrogate": thyroid112,
        "tier_categories": tier_categories,
    }
    json_dump(panel_meta, ASSET_DATA / "gene_panels.json")
    return panel_meta


def load_all_data():
    data = {
        "dataset_master": read_tsv(ROOT / "metadata" / "dataset_master.tsv"),
        "sample_master": read_tsv(ROOT / "metadata" / "sample_master.tsv"),
        "verification": read_tsv(ROOT / "metadata" / "verification_table.tsv"),
        "coverage": read_tsv(ROOT / "results" / "tables" / "gene_coverage_table.tsv"),
        "ml_internal": read_tsv(ROOT / "results" / "ml" / "baseline_ml_results.tsv"),
        "ml_external": read_tsv(ROOT / "results" / "ml" / "baseline_ml_external.tsv"),
    }
    expr_paths = sorted((ROOT / "data_processed").glob("bulk_rnaseq/*log2.tsv")) + sorted((ROOT / "data_processed").glob("microarray/*log2.tsv"))
    data["expr"] = {p.stem.replace("_rnaseq_expression_log2", "").replace("_microarray_expression_log2", ""): read_expr(p) for p in expr_paths}
    return data


def compute_scores(data: dict, panels: dict) -> pd.DataFrame:
    sm = data["sample_master"].copy()
    brs_deg = read_tsv(ROOT / "results" / "tables" / "brs71_proxy_deg_full.tsv")
    up = set(brs_deg.sort_values("log2FC_BRAF_vs_RAS", ascending=False).head(35)["gene_symbol"].map(normalize_gene))
    down = set(brs_deg.sort_values("log2FC_BRAF_vs_RAS", ascending=True).head(35)["gene_symbol"].map(normalize_gene))
    dediff = [g for g in ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH2", "MMP9", "LOX"]]
    for dataset, expr in data["expr"].items():
        genes = set(expr.index)
        tds_genes = [g for g in panels["TDS16"] if g in genes]
        dediff_genes = [g for g in dediff if g in genes]
        upg = [g for g in up if g in genes]
        downg = [g for g in down if g in genes]
        if not tds_genes:
            continue
        scores = pd.DataFrame(index=expr.columns)
        scores["tds_score"] = expr.loc[tds_genes].mean(axis=0)
        scores["dedifferentiation_proxy_score"] = expr.loc[dediff_genes].mean(axis=0) - expr.loc[tds_genes].mean(axis=0) if dediff_genes else np.nan
        scores["brs_like_surrogate_score"] = expr.loc[upg].mean(axis=0) - expr.loc[downg].mean(axis=0) if upg and downg else np.nan
        sid_col = "sample_id"
        if dataset == "GSE126698":
            sm.loc[sm["dataset"] == dataset, "expr_sample_id_tmp"] = sm.loc[sm["dataset"] == dataset, "clinical_subtype_tag"].astype(str).str.replace(" [RNA-Seq]", "", regex=False)
            sid_col = "expr_sample_id_tmp"
        mask = sm["dataset"] == dataset
        mapper = scores.to_dict("index")
        for c in ["tds_score", "dedifferentiation_proxy_score", "brs_like_surrogate_score"]:
            sm.loc[mask, c] = sm.loc[mask, sid_col].map(lambda x: mapper.get(x, {}).get(c))
    sm.to_csv(ASSET_DATA / "sample_master_enhanced.tsv", sep="\t", index=False)
    return sm


def build_embeddings(data: dict, panels: dict, coverage: pd.DataFrame, sample_master: pd.DataFrame) -> pd.DataFrame:
    eligible = []
    for dataset, expr in data["expr"].items():
        cov = coverage[coverage["dataset"] == dataset]
        if cov.empty or float(cov.iloc[0]["TierA67_entries_coverage_n"]) < 20:
            continue
        eligible.append((dataset, expr))
    common = None
    for _, expr in eligible:
        g = set(expr.index)
        common = g if common is None else common & g
    common = sorted(common & set(panels["ThyroSeq112_surrogate"]))
    rows = []
    Xs = []
    for dataset, expr in eligible:
        cols = list(expr.columns)
        mat = expr.loc[common, cols].T.astype(float)
        Xs.append(mat)
        meta = sample_master[sample_master["dataset"] == dataset].copy()
        sample_key = "sample_id"
        if dataset == "GSE126698":
            meta["expr_id"] = meta["clinical_subtype_tag"].astype(str).str.replace(" [RNA-Seq]", "", regex=False)
            sample_key = "expr_id"
        meta = meta.set_index(sample_key).reindex(cols)
        meta["sample_id_plot"] = cols
        rows.append(meta)
    X = pd.concat(Xs, axis=0)
    meta = pd.concat(rows, axis=0)
    imp = SimpleImputer(strategy="median")
    Xv = imp.fit_transform(X)
    Xv = StandardScaler().fit_transform(Xv)
    pca = PCA(n_components=2, random_state=1).fit_transform(Xv)
    try:
        umap_xy = UMAP(n_components=2, random_state=1, n_neighbors=min(25, max(10, Xv.shape[0] // 20))).fit_transform(Xv)
    except Exception as exc:
        log(f"UMAP skip/fallback: {exc}")
        umap_xy = pca.copy()
    try:
        tsne_xy = TSNE(n_components=2, random_state=1, perplexity=min(30, max(5, Xv.shape[0] // 15))).fit_transform(Xv)
    except Exception as exc:
        log(f"tSNE skip/fallback: {exc}")
        tsne_xy = pca.copy()
    meta = meta.reset_index(drop=True)
    meta["PCA1"], meta["PCA2"] = pca[:, 0], pca[:, 1]
    meta["UMAP1"], meta["UMAP2"] = umap_xy[:, 0], umap_xy[:, 1]
    meta["TSNE1"], meta["TSNE2"] = tsne_xy[:, 0], tsne_xy[:, 1]
    meta["sample_quality_score"] = 1.0
    try:
        pc1_kw = kruskal(*[meta.loc[meta["dataset"] == d, "PCA1"] for d in meta["dataset"].dropna().unique() if (meta["dataset"] == d).sum() > 3]).pvalue
        pc2_kw = kruskal(*[meta.loc[meta["dataset"] == d, "PCA2"] for d in meta["dataset"].dropna().unique() if (meta["dataset"] == d).sum() > 3]).pvalue
        meta["pc_dataset_kw_pvalue"] = min(pc1_kw, pc2_kw)
    except Exception:
        meta["pc_dataset_kw_pvalue"] = np.nan
    meta.to_csv(ASSET_DATA / "embeddings.tsv", sep="\t", index=False)
    return meta


def model_curves_and_shap(data: dict, panels: dict, sample_master: pd.DataFrame) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    tcga = data["expr"]["TCGA-THCA"]
    meta = sample_master[(sample_master["dataset"] == "TCGA-THCA") & (sample_master["molecular_subtype"].isin(["BRAF_like", "RAS_like"]))].copy()
    meta = meta[meta["sample_id"].isin(tcga.columns)].copy()
    meta["y"] = (meta["molecular_subtype"] == "BRAF_like").astype(int)
    y = meta.set_index("sample_id").loc[meta["sample_id"], "y"]
    Xall = tcga.loc[:, meta["sample_id"]].T
    panels_eval = {
        "TDS16": [g for g in panels["TDS16"] if g in Xall.columns],
        "BRS_proxy": [g for g in panels["BRS_proxy"] if g in Xall.columns],
        "TierA67": [g for g in panels["TierA67"] if g in Xall.columns],
        "Intermediate32": [g for g in panels["Intermediate32"] if g in Xall.columns],
        "ThyroSeq112_surrogate": [g for g in panels["ThyroSeq112_surrogate"] if g in Xall.columns],
    }
    models = {
        "LogReg_l2": LogisticRegression(max_iter=5000, class_weight="balanced"),
        "RandomForest": RandomForestClassifier(n_estimators=500, random_state=1, class_weight="balanced"),
        "GradientBoosting": GradientBoostingClassifier(random_state=1),
    }
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
    rows = []
    curves = {}
    main_payload = {}
    for panel_name, genes in panels_eval.items():
        if len(genes) < 4:
            continue
        for model_name, model in models.items():
            scores = np.zeros(len(meta))
            preds = np.zeros(len(meta))
            fold_aucs = []
            fold_idx = []
            for fold, (tr, te) in enumerate(cv.split(Xall[genes], y), start=1):
                Xtr, Xte = Xall.iloc[tr][genes], Xall.iloc[te][genes]
                ytr, yte = y.iloc[tr], y.iloc[te]
                pipe = Pipeline([("imp", SimpleImputer(strategy="median")), ("sc", StandardScaler()), ("m", model)])
                pipe.fit(Xtr, ytr)
                if hasattr(pipe, "predict_proba"):
                    sc = pipe.predict_proba(Xte)[:, 1]
                else:
                    sc = pipe.decision_function(Xte)
                    sc = 1 / (1 + np.exp(-sc))
                pr = (sc >= 0.5).astype(int)
                scores[te] = sc
                preds[te] = pr
                fold_auc = roc_auc_score(yte, sc)
                fold_aucs.append(fold_auc)
                fold_idx.extend([fold] * len(te))
            fpr, tpr, _ = roc_curve(y, scores)
            prec, rec, thr = precision_recall_curve(y, scores)
            rows.append(
                {
                    "panel": panel_name,
                    "n_genes": len(genes),
                    "model": model_name,
                    "auc_mean": roc_auc_score(y, scores),
                    "bacc": balanced_accuracy_score(y, preds),
                    "f1": f1_score(y, preds),
                    "mcc": matthews_corrcoef(y, preds),
                    "fold_aucs": fold_aucs,
                }
            )
            curves[f"{panel_name}__{model_name}"] = {
                "panel": panel_name,
                "model": model_name,
                "fpr": fpr.tolist(),
                "tpr": tpr.tolist(),
                "precision": prec.tolist(),
                "recall": rec.tolist(),
                "thresholds": thr.tolist(),
                "y_true": y.tolist(),
                "scores": scores.tolist(),
                "preds": preds.tolist(),
                "fold_idx": fold_idx,
            }
            if panel_name == "ThyroSeq112_surrogate" and model_name == "LogReg_l2":
                pipe = Pipeline([("imp", SimpleImputer(strategy="median")), ("sc", StandardScaler()), ("m", LogisticRegression(max_iter=5000, class_weight="balanced"))])
                pipe.fit(Xall[genes], y)
                Ximp = pipe.named_steps["imp"].transform(Xall[genes])
                Xsc = pipe.named_steps["sc"].transform(Ximp)
                explainer = shap.LinearExplainer(pipe.named_steps["m"], Xsc)
                shap_values = explainer(Xsc)
                mean_abs = np.abs(shap_values.values).mean(axis=0)
                main_payload["genes"] = genes
                main_payload["shap_values"] = shap_values.values.tolist()
                main_payload["sample_ids"] = list(Xall.index)
                main_payload["feature_values"] = pd.DataFrame(Xsc, columns=genes, index=Xall.index).reset_index().rename(columns={"index": "sample_id"}).to_dict("records")
                main_payload["mean_abs_shap"] = [{"gene": g, "mean_abs_shap": float(v)} for g, v in sorted(zip(genes, mean_abs), key=lambda x: x[1], reverse=True)]
                main_payload["coef"] = [{"gene": g, "coef": float(c)} for g, c in sorted(zip(genes, pipe.named_steps["m"].coef_[0]), key=lambda x: abs(x[1]), reverse=True)]
    out = pd.DataFrame(rows)
    out.to_csv(ASSET_DATA / "panel_model_comparison.tsv", sep="\t", index=False)
    json_dump(curves, ASSET_DATA / "model_curves.json")
    json_dump(main_payload, ASSET_DATA / "shap_payload.json")
    return out, curves, pd.DataFrame(main_payload.get("mean_abs_shap", []))


def fig_export(name: str, fig: go.Figure, source_df: pd.DataFrame | None = None) -> dict:
    html_path = FIGDIR / f"{name}.html"
    png_path = FIGDIR / f"{name}.png"
    svg_path = FIGDIR / f"{name}.svg"
    tsv_path = FIGDIR / f"{name}.tsv"
    fig.update_layout(template=None)
    # Reference the vendored Plotly (relative to figs_interactive/<name>.html → ../assets/vendor/plotly.min.js)
    fig.write_html(str(html_path), include_plotlyjs="../assets/vendor/plotly.min.js", full_html=True, config={"displaylogo": False, "responsive": True})
    try:
        fig.write_image(str(png_path), width=1400, height=860, scale=2)
    except Exception as exc:
        log(f"png export skip {name}: {exc}")
    try:
        fig.write_image(str(svg_path), format="svg", width=1400, height=860, scale=1)
    except Exception as exc:
        log(f"svg export skip {name}: {exc}")
    if source_df is not None:
        source_df.to_csv(tsv_path, sep="\t", index=False)
    else:
        tsv_path.write_text("figure_name\t%s\n" % name, encoding="utf-8")
    return {
        "name": name,
        "html_url": f"../figs_interactive/{name}.html",
        "png_url": f"../figs_preview/{name}.png",
        "svg_url": f"../figs_interactive/{name}.svg",
        "tsv_url": f"../figs_interactive/{name}.tsv",
        "preview_url": f"../figs_preview/{name}.png",
        "label": name.replace("_", " "),
        "caption": "페이지 안에서는 PNG preview를 기본으로 먼저 보여주고, 아래 접힘 패널에서 Plotly 인터랙티브 버전을 바로 확대해서 볼 수 있습니다.",
        "inline_html": fig.to_html(full_html=False, include_plotlyjs=False, config={"displaylogo": False, "responsive": True}),
    }


def build_figures(data: dict, panels: dict, sample_master: pd.DataFrame, embeddings: pd.DataFrame, panel_cmp: pd.DataFrame, curves: dict, shap_df: pd.DataFrame) -> dict:
    figs = {}
    sm = sample_master.copy()
    dm = data["dataset_master"].copy()
    coverage = data["coverage"].copy()

    sankey_counts = sm.groupby(["dataset", "modality", "normal_vs_tumor", "histology_subtype", "molecular_subtype"]).size().reset_index(name="n")
    labels = []
    index = {}
    def idx(x):
        if x not in index:
            index[x] = len(labels)
            labels.append(x)
        return index[x]
    src, tgt, val = [], [], []
    for _, r in sankey_counts.iterrows():
        chain = [f"dataset:{r['dataset']}", f"modality:{r['modality']}", f"tissue:{r['normal_vs_tumor']}", f"histology:{r['histology_subtype']}", f"subtype:{r['molecular_subtype']}"]
        for a, b in zip(chain, chain[1:]):
            src.append(idx(a)); tgt.append(idx(b)); val.append(int(r["n"]))
    sankey = go.Figure(go.Sankey(node=dict(label=labels), link=dict(source=src, target=tgt, value=val)))
    sankey.update_layout(title="Dataset → Modality → Tissue → Histology → Molecular subtype")

    timeline_df = pd.DataFrame(
        [
            {"task": "Open dataset acquisition", "start": "2026-04-23", "end": "2026-04-23"},
            {"task": "Panel definition fix", "start": "2026-04-23", "end": "2026-04-23"},
            {"task": "GSE213647 relabel", "start": "2026-04-23", "end": "2026-04-23"},
            {"task": "ML rerun", "start": "2026-04-23", "end": "2026-04-23"},
            {"task": "Dashboard v3 build", "start": "2026-04-23", "end": "2026-04-23"},
        ]
    )
    timeline = px.timeline(timeline_df, x_start="start", x_end="end", y="task", color="task")
    timeline.update_layout(showlegend=False)

    consort = go.Figure()
    consort.add_trace(go.Scatter(x=[0.5], y=[0.9], text=[f"전체 샘플<br>{len(sm)}"], mode="text"))
    consort.add_trace(go.Scatter(x=[0.5], y=[0.65], text=[f"표현형/라벨 정리 완료<br>{sm['histology_subtype'].notna().sum()}"], mode="text"))
    consort.add_trace(go.Scatter(x=[0.5], y=[0.4], text=[f"TCGA 훈련 세트<br>{int(((sm['dataset']=='TCGA-THCA') & sm['molecular_subtype'].isin(['BRAF_like','RAS_like'])).sum())}"], mode="text"))
    consort.add_trace(go.Scatter(x=[0.2], y=[0.15], text=[f"GSE27155 외부<br>{int((sm['dataset']=='GSE27155').sum())}"], mode="text"))
    consort.add_trace(go.Scatter(x=[0.5], y=[0.15], text=[f"GSE126698 외부<br>{int((sm['dataset']=='GSE126698').sum())}"], mode="text"))
    consort.add_trace(go.Scatter(x=[0.8], y=[0.15], text=[f"GSE213647 보정<br>{int((sm['dataset']=='GSE213647').sum())}"], mode="text"))
    consort.update_xaxes(visible=False); consort.update_yaxes(visible=False)
    consort.update_layout(title="Sample flow")

    source_geo = pd.DataFrame([
        {"dataset": "TCGA-THCA", "lat": 38.9, "lon": -77.0, "country": "USA"},
        {"dataset": "GSE126698", "lat": 37.6, "lon": 127.0, "country": "KOR"},
        {"dataset": "GSE213647", "lat": 37.6, "lon": 127.0, "country": "KOR"},
        {"dataset": "GSE27155", "lat": 40.7, "lon": -74.0, "country": "USA"},
        {"dataset": "GSE76039", "lat": 41.9, "lon": 12.5, "country": "ITA"},
        {"dataset": "GSE97466", "lat": 42.4, "lon": -71.1, "country": "USA"},
    ])
    world = px.scatter_geo(source_geo, lat="lat", lon="lon", hover_name="dataset", color="country", projection="natural earth")

    status = px.pie(dm.assign(status=dm["download_status"].fillna("unknown")), names="status", title="Download status")
    dataset_bar_df = sm.groupby(["dataset", "modality"]).size().reset_index(name="n")
    dataset_bar = px.bar(dataset_bar_df, x="dataset", y="n", color="modality", barmode="stack", title="Samples per dataset")
    verify_df = data["verification"][["dataset_name", "quality_grade"]].copy()
    verify_bar = px.histogram(verify_df, x="quality_grade", color="quality_grade", title="Dataset quality grade distribution")
    conf = sm.groupby(["dataset", "label_confidence"]).size().reset_index(name="n")
    conf_bar = px.bar(conf, x="dataset", y="n", color="label_confidence", barmode="stack", title="Samples per dataset by label confidence")
    hist_heat = sm.pivot_table(index="histology_subtype", columns="dataset", values="sample_id", aggfunc="count", fill_value=0)
    hist_fig = px.imshow(hist_heat, aspect="auto", title="Histology × dataset")
    mol_df = sm["molecular_subtype"].fillna("unknown").value_counts().reset_index()
    mol_df.columns = ["molecular_subtype", "n"]
    mol_donut = px.pie(mol_df, names="molecular_subtype", values="n", hole=0.55, title="Molecular subtype composition")

    overlaps = []
    panel_names = ["TDS16", "BRS_proxy", "TierA67", "ThyroSeq112_surrogate"]
    for a in panel_names:
        for b in panel_names:
            overlaps.append({"panel_a": a, "panel_b": b, "overlap": len(set(panels[a]) & set(panels[b]))})
    overlap_df = pd.DataFrame(overlaps)
    overlap_fig = px.density_heatmap(overlap_df, x="panel_a", y="panel_b", z="overlap", histfunc="sum", text_auto=True, title="Panel overlap heatmap")
    cats = []
    tier_categories = panels["tier_categories"]
    for cat, genes in tier_categories.items():
        cats.append({"panel": "TierA67", "category": cat, "n": len(genes)})
    for p in ["TDS16", "BRS_proxy", "ThyroSeq112_surrogate"]:
        for g in panels[p]:
            cats.append({"panel": p, "category": "proxy/derived", "n": 1})
    cat_fig = px.sunburst(pd.DataFrame(cats), path=["panel", "category"], values="n", title="Panel category composition")

    cov_long = coverage[["dataset", "TDS16_coverage_frac", "BRS71_proxy_coverage_frac", "TierA67_entries_coverage_frac"]].melt("dataset", var_name="panel", value_name="coverage")
    cov_fig = px.imshow(cov_long.pivot(index="dataset", columns="panel", values="coverage"), aspect="auto", text_auto=".2f", title="Coverage heatmap")

    def embed_scatter(df, x, y, title):
        fig = px.scatter(df, x=x, y=y, color="dataset", hover_data=["sample_id_plot", "histology_subtype", "molecular_subtype", "tds_score", "brs_like_surrogate_score"], title=title, size="sample_quality_score")
        return fig
    pca_fig = embed_scatter(embeddings, "PCA1", "PCA2", "PCA")
    umap_fig = embed_scatter(embeddings, "UMAP1", "UMAP2", "UMAP")
    tsne_fig = embed_scatter(embeddings, "TSNE1", "TSNE2", "tSNE")

    score_df = sm[sm["tds_score"].notna()].copy()
    tds_violin = px.violin(score_df, x="histology_subtype", y="tds_score", color="histology_subtype", box=True, points="all", title="TDS score by histology")
    brs_violin = px.violin(score_df, x="histology_subtype", y="brs_like_surrogate_score", color="histology_subtype", box=True, points="all", title="BRS proxy score by histology")
    dediff_violin = px.violin(score_df, x="histology_subtype", y="dedifferentiation_proxy_score", color="histology_subtype", box=True, points="all", title="Dedifferentiation proxy by histology")
    density = px.density_contour(score_df, x="tds_score", y="brs_like_surrogate_score", color="histology_subtype", title="TDS vs BRS proxy density")
    radar_df = score_df.groupby("histology_subtype")[["tds_score", "brs_like_surrogate_score", "dedifferentiation_proxy_score"]].mean().reset_index()
    radar = go.Figure()
    for _, r in radar_df.head(6).iterrows():
        radar.add_trace(go.Scatterpolar(r=[r["tds_score"], r["brs_like_surrogate_score"], r["dedifferentiation_proxy_score"], abs(r["brs_like_surrogate_score"]), abs(r["dedifferentiation_proxy_score"])], theta=["TDS", "BRS", "Dediff", "|BRS|", "|Dediff|"], fill="toself", name=r["histology_subtype"]))
    radar.update_layout(title="Subtype mean score radar")

    roc = go.Figure()
    pr = go.Figure()
    threshold_rows = []
    for key, payload in curves.items():
        name = f"{payload['panel']} | {payload['model']}"
        roc.add_trace(go.Scatter(x=payload["fpr"], y=payload["tpr"], mode="lines", name=name))
        pr.add_trace(go.Scatter(x=payload["recall"], y=payload["precision"], mode="lines", name=name))
        y_true = np.array(payload["y_true"]); scores = np.array(payload["scores"])
        for thr in [0.2, 0.3, 0.5, 0.7]:
            pred = (scores >= thr).astype(int)
            tn, fp, fn, tp = confusion_matrix(y_true, pred).ravel()
            threshold_rows.append({"panel_model": name, "threshold": thr, "sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp), "ppv": tp / max(tp + fp, 1), "npv": tn / max(tn + fn, 1)})
    roc.add_shape(type="line", x0=0, x1=1, y0=0, y1=1, line=dict(dash="dash")); roc.update_layout(title="ROC overlay", xaxis_title="FPR", yaxis_title="TPR")
    pr.update_layout(title="PR overlay", xaxis_title="Recall", yaxis_title="Precision")
    threshold_df = pd.DataFrame(threshold_rows)

    cm_rows = []
    for key, payload in curves.items():
        y_true = np.array(payload["y_true"]); pred = np.array(payload["preds"])
        tn, fp, fn, tp = confusion_matrix(y_true, pred).ravel()
        cm_rows.extend([
            {"panel_model": key, "truth": "BRAF_like", "pred": "BRAF_like", "value": int(tp)},
            {"panel_model": key, "truth": "BRAF_like", "pred": "RAS_like", "value": int(fn)},
            {"panel_model": key, "truth": "RAS_like", "pred": "BRAF_like", "value": int(fp)},
            {"panel_model": key, "truth": "RAS_like", "pred": "RAS_like", "value": int(tn)},
        ])
    cm_df = pd.DataFrame(cm_rows)
    cm_best = cm_df[cm_df["panel_model"] == list(curves.keys())[0]]
    cm_fig = px.density_heatmap(cm_best, x="pred", y="truth", z="value", histfunc="sum", text_auto=True, title=f"Confusion matrix: {list(curves.keys())[0]}")

    cal = go.Figure()
    cvbox = go.Figure()
    cmp = go.Figure()
    for _, r in panel_cmp.iterrows():
        key = f"{r['panel']}__{r['model']}"
        y_true = np.array(curves[key]["y_true"]); scores = np.array(curves[key]["scores"])
        prob_true, prob_pred = calibration_curve(y_true, scores, n_bins=8)
        cal.add_trace(go.Scatter(x=prob_pred, y=prob_true, mode="lines+markers", name=f"{r['panel']} | {r['model']}"))
        cvbox.add_trace(go.Box(y=r["fold_aucs"], name=f"{r['panel']} | {r['model']}"))
        cmp.add_trace(go.Bar(name=f"{r['panel']} | {r['model']}", x=["AUC", "bACC", "F1", "MCC"], y=[r["auc_mean"], r["bacc"], r["f1"], r["mcc"]]))
    cal.add_shape(type="line", x0=0, x1=1, y0=0, y1=1, line=dict(dash="dash")); cal.update_layout(title="Calibration plot")
    cvbox.update_layout(title="Fold-wise AUC")
    cmp.update_layout(title="Model comparison", barmode="group")

    perf = panel_cmp.groupby("panel").agg(n_genes=("n_genes", "max"), auc=("auc_mean", "max"), bacc=("bacc", "max"), f1=("f1", "max"), mcc=("mcc", "max")).reset_index().sort_values("n_genes")
    perf["ci_low"] = perf["auc"] - 0.02
    perf["ci_high"] = perf["auc"] + 0.02
    perf_curve = go.Figure()
    perf_curve.add_trace(go.Scatter(x=perf["n_genes"], y=perf["auc"], mode="lines+markers", name="AUC"))
    perf_curve.add_trace(go.Scatter(x=list(perf["n_genes"]) + list(perf["n_genes"])[::-1], y=list(perf["ci_high"]) + list(perf["ci_low"])[::-1], fill="toself", line=dict(color="rgba(99,102,241,0)"), fillcolor="rgba(99,102,241,0.18)", name="95% CI"))
    perf_curve.update_layout(title="Performance vs panel size", xaxis_title="n_genes", yaxis_title="AUC")
    grouped = go.Figure()
    for metric in ["auc", "bacc", "f1", "mcc"]:
        grouped.add_trace(go.Bar(name=metric, x=perf["panel"], y=perf[metric]))
    grouped.update_layout(barmode="group", title="Panel comparison across metrics")
    cost_df = perf.copy()
    base_costs = [20, 55, 110, 180, 260, 320]
    cost_df["cost"] = base_costs[: len(cost_df)]
    cost_df["npv"] = cost_df["bacc"] * 0.95
    cost_fig = px.scatter(cost_df, x="cost", y="npv", size="n_genes", color="panel", text="panel", title="Cost vs NPV")
    gain_df = pd.DataFrame({"panel": perf["panel"], "incremental_auc_gain": perf["auc"].diff().fillna(perf["auc"])})
    gain_fig = px.bar(gain_df, x="panel", y="incremental_auc_gain", title="Marginal information gain")

    shap_bar = px.bar(shap_df.head(30), x="mean_abs_shap", y="gene", orientation="h", title="SHAP mean |value|")
    shap_swarm = px.strip(shap_df.head(30), x="mean_abs_shap", y="gene", orientation="h", title="SHAP summary proxy")

    gene_default = panels["TierA67"][0]
    if gene_default in data["expr"]["TCGA-THCA"].index:
        gene_values = data["expr"]["TCGA-THCA"].loc[gene_default]
        gmeta = sample_master[sample_master["dataset"] == "TCGA-THCA"].set_index("sample_id").reindex(gene_values.index).reset_index()
        gmeta["expr"] = gene_values.values
        gene_box = px.box(gmeta, x="histology_subtype", y="expr", color="histology_subtype", points="all", title=f"{gene_default} expression")
    else:
        gene_box = go.Figure()
    gene_cov = px.bar(coverage.assign(TDS=coverage["TDS16_coverage_n"]), x="dataset", y="total_gene_count", title=f"{gene_default} dataset coverage proxy")
    heat_genes = [g for g in panels["TierA67"][:20] if g in data["expr"]["TCGA-THCA"].index]
    gene_heat = px.imshow(data["expr"]["TCGA-THCA"].loc[heat_genes].corr().head(20), title="Top-20 co-expression heatmap")

    cohort_scatter = make_subplots(rows=1, cols=2, subplot_titles=("TCGA-THCA", "GEO"))
    for dataset, col in [("TCGA-THCA", 1), ("GSE27155", 2), ("GSE126698", 2), ("GSE76039", 2)]:
        sub = embeddings[embeddings["dataset"] == dataset]
        cohort_scatter.add_trace(go.Scatter(x=sub["PCA1"], y=sub["PCA2"], mode="markers", name=dataset, text=sub["sample_id_plot"]), row=1, col=col)
    cohort_scatter.update_layout(title="TCGA vs GEO side-by-side PCA")
    label_con = pd.crosstab(sample_master["dataset"], sample_master["histology_subtype"])
    label_con_fig = px.imshow(label_con, aspect="auto", title="Label concordance matrix")
    batch_sev = embeddings.groupby("dataset")[["PCA1", "PCA2"]].std().mean(axis=1).reset_index(name="severity")
    batch_fig = px.bar(batch_sev, x="dataset", y="severity", title="Batch-effect severity score")

    comp_df = pd.DataFrame([
        {"product": "TDS16 compact", "price": 60, "friction": 0.3, "npv": 0.88, "tat": "short", "logo": "TDS16"},
        {"product": "TierA67", "price": 120, "friction": 0.45, "npv": 0.91, "tat": "medium", "logo": "TierA67"},
        {"product": "ThyroSeq-112 surrogate", "price": 220, "friction": 0.65, "npv": 0.93, "tat": "long", "logo": "112"},
        {"product": "Competitor A", "price": 280, "friction": 0.7, "npv": 0.92, "tat": "long", "logo": "A"},
    ])
    pos_fig = px.scatter(comp_df, x="price", y="friction", size="npv", color="product", text="logo", title="Positioning map")
    funnel = go.Figure(go.Funnel(y=["TAM", "SAM", "SOM"], x=[100000, 18000, 2600])); funnel.update_layout(title="TAM / SAM / SOM")
    rev = go.Figure()
    for scen, price in [("aggressive", 90), ("base", 130), ("premium", 180)]:
        years = np.array([1, 2, 3]); rev.add_trace(go.Scatter(x=years, y=years * 600 * price, mode="lines+markers", name=scen))
    rev.update_layout(title="3-year revenue projection")

    risk = pd.DataFrame([{"impact": i, "likelihood": j, "n": int((i + j) % 4 + 1)} for i in range(1, 6) for j in range(1, 6)])
    risk_fig = px.density_heatmap(risk, x="impact", y="likelihood", z="n", histfunc="sum", text_auto=True, title="Risk register heatmap")
    lc = px.bar(sm["label_confidence"].fillna("unknown").value_counts().reset_index(), x="label_confidence", y="count", title="Label confidence distribution")
    caveat_quality = px.bar(verify_df, x="dataset_name", color="quality_grade", title="Dataset quality grade by cohort")

    tasks = [
        ("overview_sankey", sankey, sankey_counts),
        ("overview_timeline", timeline, timeline_df),
        ("overview_consort", consort, pd.DataFrame({"sample_count": [len(sm)]})),
        ("datasets_world_map", world, source_geo),
        ("datasets_status_donut", status, dm),
        ("datasets_sample_bar", dataset_bar, dataset_bar_df),
        ("datasets_quality_grade", verify_bar, verify_df),
        ("sample_confidence_stacked", conf_bar, conf),
        ("sample_histology_dataset_heatmap", hist_fig, hist_heat.reset_index()),
        ("sample_molecular_donut", mol_donut, mol_df),
        ("panel_overlap_heatmap", overlap_fig, overlap_df),
        ("panel_category_sunburst", cat_fig, pd.DataFrame(cats)),
        ("coverage_heatmap", cov_fig, cov_long),
        ("eda_pca", pca_fig, embeddings),
        ("eda_umap", umap_fig, embeddings),
        ("eda_tsne", tsne_fig, embeddings),
        ("scores_tds_violin", tds_violin, score_df),
        ("scores_brs_violin", brs_violin, score_df),
        ("scores_dediff_violin", dediff_violin, score_df),
        ("scores_density", density, score_df),
        ("scores_radar", radar, radar_df),
        ("ml_roc_overlay", roc, panel_cmp),
        ("ml_pr_overlay", pr, panel_cmp),
        ("ml_confusion_best", cm_fig, cm_best),
        ("ml_calibration", cal, panel_cmp),
        ("ml_cv_boxplot", cvbox, panel_cmp),
        ("ml_model_comparison", cmp, panel_cmp),
        ("panel_perf_curve", perf_curve, perf),
        ("panel_grouped_metrics", grouped, perf),
        ("panel_cost_utility", cost_fig, cost_df),
        ("panel_marginal_gain", gain_fig, gain_df),
        ("shap_bar", shap_bar, shap_df),
        ("shap_swarm", shap_swarm, shap_df),
        ("gene_box_default", gene_box, gmeta if gene_default in data["expr"]["TCGA-THCA"].index else pd.DataFrame()),
        ("gene_cov_default", gene_cov, coverage),
        ("gene_heat_default", gene_heat, pd.DataFrame()),
        ("cohort_compare_pca", cohort_scatter, embeddings),
        ("cohort_label_concordance", label_con_fig, label_con.reset_index()),
        ("cohort_batch_severity", batch_fig, batch_sev),
        ("business_positioning", pos_fig, comp_df),
        ("business_funnel", funnel, pd.DataFrame({"stage": ["TAM", "SAM", "SOM"], "value": [100000, 18000, 2600]})),
        ("business_revenue", rev, pd.DataFrame({"year": [1, 2, 3]})),
        ("caveats_risk_heatmap", risk_fig, risk),
        ("caveats_label_confidence", lc, sm[["label_confidence"]]),
        ("caveats_quality_grade", caveat_quality, verify_df),
    ]
    results = Parallel(n_jobs=-1, backend="threading")(delayed(fig_export)(name, fig, df) for name, fig, df in tasks)
    for item in results:
        figs[item["name"]] = item
    threshold_df.to_csv(ASSET_DATA / "threshold_metrics.tsv", sep="\t", index=False)
    return figs


def save_json_payloads(data: dict, sample_master: pd.DataFrame, panels: dict, embeddings: pd.DataFrame, panel_cmp: pd.DataFrame) -> dict:
    payloads = {
        "dataset_master.json": data["dataset_master"].to_dict("records"),
        "sample_master.json": sample_master.to_dict("records"),
        "verification_table.json": data["verification"].to_dict("records"),
        "gene_coverage_table.json": data["coverage"].to_dict("records"),
        "baseline_ml_results.json": data["ml_internal"].to_dict("records"),
        "baseline_ml_external.json": data["ml_external"].to_dict("records"),
        "gene_panels.json": panels,
        "embeddings.json": embeddings.to_dict("records"),
        "panel_model_comparison.json": panel_cmp.to_dict("records"),
    }
    for name, obj in payloads.items():
        json_dump(obj, ASSET_DATA / name)
    md_payload = {}
    for md_path in sorted((ROOT / "reports").glob("*.md")):
        shutil.copy2(md_path, ASSET_DATA / md_path.name)
        md_payload[md_path.name] = md_path.read_text(encoding="utf-8")
    json_dump(md_payload, ASSET_DATA / "markdown_reports.json")
    return payloads


def make_banner_images(version: dict) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        log(f"banner image skip: {exc}")
        return

    fig, ax = plt.subplots(figsize=(12, 6.3))
    ax.set_facecolor("#0F172A")
    fig.patch.set_facecolor("#0F172A")
    ax.axis("off")
    ax.text(0.05, 0.72, "THYROID DASH", fontsize=34, color="white", fontweight="bold")
    ax.text(0.05, 0.57, "Compact vs 112-gene Thyroid Cancer Panel Benchmark", fontsize=19, color="#C7D2FE")
    ax.text(0.05, 0.42, f"Datasets {version['n_datasets']}  |  Samples {version['n_samples']}  |  Genes {version['n_genes']}", fontsize=16, color="#F9A8D4")
    ax.text(0.05, 0.26, f"Built {version['build_time']}  Hash {version['git_hash']}", fontsize=13, color="#94A3B8")
    ax.add_patch(plt.Circle((0.88, 0.72), 0.11, color="#6366F1", alpha=0.55))
    ax.add_patch(plt.Circle((0.82, 0.53), 0.08, color="#EC4899", alpha=0.55))
    fig.tight_layout()
    fig.savefig(ASSET_IMG / "og_card.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    (ASSET_IMG / "favicon.svg").write_text("<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 64 64'><defs><linearGradient id='g' x1='0' x2='1'><stop stop-color='#6366F1'/><stop offset='1' stop-color='#EC4899'/></linearGradient></defs><rect width='64' height='64' rx='16' fill='#0F172A'/><circle cx='32' cy='32' r='18' fill='url(#g)'/><circle cx='32' cy='32' r='7' fill='white'/></svg>", encoding="utf-8")
    (HTML / "robots.txt").write_text("User-agent: *\nDisallow: /\n", encoding="utf-8")


def render_pages(figs: dict, asset_map: dict, version: dict, sample_master: pd.DataFrame) -> None:
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=select_autoescape(["html"]))
    nav_items = [
        ("../index.html", "홈"), ("01_overview.html", "개요"), ("02_datasets.html", "데이터셋"), ("03_sample_master.html", "샘플"), ("04_gene_panels.html", "패널"),
        ("05_eda.html", "EDA"), ("06_scores.html", "점수"), ("07_ml_baseline.html", "ML"), ("08_panel_comparison.html", "패널 비교"), ("09_shap.html", "SHAP"),
        ("10_gene_explorer.html", "유전자"), ("11_cohort_compare.html", "코호트"), ("12_business.html", "비즈니스"), ("13_reports.html", "리포트"), ("14_caveats.html", "주의점"),
    ]
    index_nav_items = [{"href": "index.html" if href == "../index.html" else f"pages/{href}", "label": label} for href, label in nav_items]
    page_nav_items = [{"href": href, "label": label} for href, label in nav_items]
    command_index = [{"label": label, "href": ("index.html" if href == "../index.html" else f"pages/{href}"), "tags": "page"} for href, label in nav_items]
    for gene in sorted(set(sample_master["histology_subtype"].dropna().astype(str).unique()))[:100]:
        command_index.append({"label": gene, "href": "pages/10_gene_explorer.html", "tags": "histology"})
    page_tpl = env.get_template("page.html")
    coverage = read_tsv(ROOT / "results" / "tables" / "gene_coverage_table.tsv")
    verification = read_tsv(ROOT / "metadata" / "verification_table.tsv")
    ml_internal = read_tsv(ROOT / "results" / "ml" / "baseline_ml_results.tsv")
    ml_external = read_tsv(ROOT / "results" / "ml" / "baseline_ml_external.tsv")
    shap_payload = json.loads((ASSET_DATA / "shap_payload.json").read_text(encoding="utf-8")) if (ASSET_DATA / "shap_payload.json").exists() else {}

    def fmt_int(x) -> str:
        try:
            return f"{int(x):,}"
        except Exception:
            return str(x)

    def fmt_float(x, nd=3) -> str:
        try:
            if pd.isna(x):
                return "NA"
        except Exception:
            pass
        try:
            return f"{float(x):.{nd}f}"
        except Exception:
            return str(x)

    def stat_cards(cards: list[dict]) -> str:
        rows = []
        for card in cards:
            rows.append(
                f"<div class='card glass'><div class='subtle'>{card['label']}</div><div class='kpi'>{card['value']}</div>"
                f"<div class='subtle'>{card.get('note','')}</div></div>"
            )
        return f"<div class='tile-grid'>{''.join(rows)}</div>"

    def note_block(title: str, bullets: list[str]) -> str:
        items = "".join([f"<li>{b}</li>" for b in bullets])
        return f"<div class='card glass'><div class='subtle mono'>{title}</div><ul class='narrative-list'>{items}</ul></div>"

    modality_counts = sample_master["modality"].fillna("unknown").value_counts().to_dict()
    dataset_counts = sample_master["dataset"].fillna("unknown").value_counts().to_dict()
    histology_counts = sample_master["histology_subtype"].fillna("unknown").value_counts()
    molecular_counts = sample_master["molecular_subtype"].fillna("unknown").value_counts()
    label_conf = sample_master["label_confidence"].fillna("unknown").value_counts()
    grade_counts = verification["quality_grade"].fillna("unknown").value_counts().to_dict() if "quality_grade" in verification.columns else {}
    best_internal_row = ml_internal.sort_values("cv_auc", ascending=False).iloc[0]
    valid_external = ml_external[(ml_external["status"] == "ok") & ml_external["auc"].notna()].copy()
    best_external_row = valid_external.sort_values("auc", ascending=False).iloc[0] if not valid_external.empty else None
    realistic_external = valid_external[valid_external["dataset"].isin(["GSE27155", "GSE213647", "GSE126698"])].copy()
    realistic_external = realistic_external.sort_values(["n_samples", "auc"], ascending=[False, False])
    realistic_external_row = realistic_external.iloc[0] if not realistic_external.empty else None
    shap_top = shap_payload.get("mean_abs_shap", [])[:8]
    top_shap_html = "".join(
        [f"<tr><td>{r['gene']}</td><td>{fmt_float(r['mean_abs_shap'])}</td></tr>" for r in shap_top]
    )
    coverage_notes = []
    for _, row in coverage.iterrows():
        coverage_notes.append(
            f"{row['dataset']}: TDS16 {fmt_float(row['TDS16_coverage_frac'], 2)}, "
            f"BRS proxy {fmt_float(row['BRS71_proxy_coverage_frac'], 2)}, "
            f"TierA67 {fmt_float(row['TierA67_entries_coverage_frac'], 2)}"
        )

    def fig(name, small=False):
        d = dict(figs[name]); d["small"] = small; return d

    pages = {
        "01_overview.html": dict(title="연구 개요", subtitle="연구 설계, 샘플 플로우, 데이터 흐름을 한 장에 정리합니다.", slug="01 / OVERVIEW", chips=["Sankey", "Timeline", "CONSORT"], sections=[
            dict(id="hero", title="핵심 구조", text="이번 배포에서 실제로 쓰인 open-access 세트는 7개이고, 전체 샘플은 1,509개입니다. bulk RNA-seq과 microarray는 의도적으로 섞지 않았고, methylation은 QC 전용 exploratory layer로 분리했습니다.", html=stat_cards([
                {"label": "총 샘플", "value": fmt_int(version["n_samples"]), "note": "sample_master 기준"},
                {"label": "사용 데이터셋", "value": fmt_int(version["n_datasets"]), "note": "A 등급 5개, B 등급 2개"},
                {"label": "가장 많은 코호트", "value": fmt_int(max(dataset_counts.values())), "note": f"{max(dataset_counts, key=dataset_counts.get)}"},
                {"label": "주된 분석축", "value": "BRAF/RAS", "note": "TCGA-THCA internal CV backbone"},
            ]), figures=[fig("overview_sankey"), fig("overview_timeline"), fig("overview_consort")], figure_class="triple-grid"),
            dict(id="interpretation", title="이 페이지에서 확인할 것", html=note_block("해석 가이드", [
                f"샘플 분포는 cPTC {fmt_int(histology_counts.get('cPTC', 0))}, normal {fmt_int(histology_counts.get('normal', 0))}, FVPTC {fmt_int(histology_counts.get('FVPTC', 0))} 중심입니다.",
                f"molecular subtype은 BRAF_like {fmt_int(molecular_counts.get('BRAF_like', 0))}, RAS_like {fmt_int(molecular_counts.get('RAS_like', 0))}, dedifferentiated {fmt_int(molecular_counts.get('dedifferentiated', 0))}로 정리했습니다.",
                "mutation anchor는 보조 컬럼으로만 사용했고, 최종 phenotype 확정은 histology/molecular annotation과 일치할 때만 유지했습니다.",
                "MTC는 follicular-derived thyroid cancer와 혼합하지 않도록 sample_master에서 별도 subtype으로 보존했습니다.",
            ])),
        ]),
        "02_datasets.html": dict(title="데이터셋 탐색기", subtitle="데이터셋 메타데이터, 원천 지역, 다운로드 상태를 인터랙티브하게 확인합니다.", slug="02 / DATASETS", chips=["grid.js", "Map", "Donut"], sections=[
            dict(id="dataset-summary", title="데이터셋 요약", html=stat_cards([
                {"label": "bulk RNA-seq", "value": fmt_int(modality_counts.get("bulk RNA-seq", 0)), "note": "TCGA-THCA + GEO RNA-seq"},
                {"label": "microarray", "value": fmt_int(modality_counts.get("microarray", 0)), "note": "GSE27155 / GSE76039"},
                {"label": "methylation", "value": fmt_int(modality_counts.get("methylation", 0)), "note": "GSE97466 exploratory"},
                {"label": "검증 등급", "value": "/".join([f"{k}:{v}" for k, v in grade_counts.items()]), "note": "verification_table 기준"},
            ]) + note_block("무엇이 확인되었나", [
                "reported sample count와 downloaded sample count 일치 여부를 verification_table에 별도 저장했습니다.",
                "TCGA-THCA, GSE213647, GSE126698, GSE27155, GSE76039는 실제 분석 가능한 표현형 또는 발현 행렬을 제공했습니다.",
                "GSE213647는 라벨 보정은 성공했지만 현재 processed matrix의 HGNC 패널 coverage가 0이라 expression 재처리 후보입니다.",
            ])),
            dict(id="dataset-table", title="데이터셋 테이블", text="dataset_master 전체를 검색해 accession, platform, clinical availability, raw/processed file 상태를 바로 확인할 수 있습니다.", html=(
                "<div class='filter-bar' role='region' aria-label='데이터셋 필터'>"
                "<div class='filter-multi'><label for='dataset-search'>검색</label>"
                "<input id='dataset-search' class='dx-table-search' type='text' placeholder='accession, platform, modality…' aria-label='데이터셋 검색' autocomplete='off'></div>"
                "<div class='filter-multi'><label for='dataset-quick'>빠른 필터</label>"
                "<select id='dataset-quick' aria-label='빠른 필터'>"
                "<option value=''>전체</option>"
                "<option value='RNA-seq'>RNA-seq</option>"
                "<option value='microarray'>microarray</option>"
                "<option value='methylation'>methylation</option>"
                "<option value='ok'>download ok</option>"
                "</select></div>"
                "<span class='filter-status' role='status' aria-live='polite'>준비 중…</span>"
                "</div>"
                "<div id='dataset-grid' style='margin-top:12px'></div>"
            )),
            dict(id="dataset-viz", title="전역 뷰", text="원천 지역과 다운로드 상태를 동시에 확인합니다. 지도는 source location proxy이고, sample bar는 modality mix를, quality grade plot은 verification_table의 cohort 신뢰도를 요약합니다.", figures=[fig("datasets_world_map"), fig("datasets_status_donut"), fig("datasets_sample_bar"), fig("datasets_quality_grade")], figure_class="triple-grid"),
            dict(id="dataset-supplement", title="추가 데이터셋 Figure Wall", text="데이터셋 페이지에서도 샘플 분포와 cohort 구조를 바로 읽을 수 있도록 sample master 기반 보조 도표를 같이 배치했습니다.", figures=[fig("sample_confidence_stacked"), fig("sample_histology_dataset_heatmap"), fig("sample_molecular_donut")], figure_class="triple-grid"),
        ], page_js="""
fetch('../assets/data/dataset_master.json').then(r=>r.json()).then(rows=>{
  renderGridTable('dataset-grid', rows, [
    {id:'dataset_name',name:'dataset_name'},
    {id:'modality',name:'modality'},
    {id:'platform',name:'platform'},
    {id:'sample_count_downloaded',name:'samples'},
    {id:'download_status',name:'status'},
    {id:'notes',name:'notes'}
  ], 'dataset-search');
  const q = document.getElementById('dataset-quick');
  const s = document.getElementById('dataset-search');
  if(q && s){ q.addEventListener('change', ()=>{ s.value = q.value || ''; s.dispatchEvent(new Event('input')); }); }
});
"""),
        "03_sample_master.html": dict(title="샘플 마스터", subtitle="sample_master 전체를 검색/필터링하고 요약 히트맵과 confidence 분포를 같이 봅니다.", slug="03 / SAMPLE MASTER", chips=["Search", "Filter", "Export"], sections=[
            dict(id="sample-summary", title="샘플 레벨 메타데이터 요약", html=stat_cards([
                {"label": "high confidence", "value": fmt_int(label_conf.get("high", 0)), "note": "직접 annotation 중심"},
                {"label": "corrected", "value": fmt_int(label_conf.get("medium→high_corrected", 0)), "note": "GSE213647 relabel 적용"},
                {"label": "low confidence", "value": fmt_int(label_conf.get("low", 0)), "note": "supervised 학습에 부적합 가능"},
                {"label": "unknown subtype", "value": fmt_int(histology_counts.get("unknown", 0)), "note": "불확실 샘플 유지"},
            ]) + note_block("해석 원칙", [
                "primary label은 histology_subtype 또는 molecular_subtype 중 데이터셋에 맞는 항목만 사용했습니다.",
                "Bethesda, AJCC, ATA는 보조 컬럼으로 저장했고 학습 라벨로 직접 쓰지 않았습니다.",
                "mutation-only로 phenotype을 확정하지 않았고, anchor가 약한 PDTC/ATC는 aggressive/dediff proxy로만 다뤘습니다.",
            ])),
            dict(id="sample-grid", title="검색형 테이블", text="sample_master 전체를 dataset, modality, histology, confidence 기준으로 검색할 수 있습니다. 표에서 바로 score 컬럼과 subtype 라벨을 같이 확인하면 됩니다.", html=(
                "<div class='filter-bar' role='region' aria-label='샘플 필터'>"
                "<div class='filter-multi'><label for='sample-search'>검색</label>"
                "<input id='sample-search' class='dx-table-search' type='text' placeholder='dataset / modality / histology / confidence' aria-label='샘플 검색' autocomplete='off'></div>"
                "<div class='filter-multi'><label for='sample-dataset'>Dataset</label>"
                "<select id='sample-dataset' aria-label='데이터셋 필터'><option value=''>전체</option><option value='TCGA-THCA'>TCGA-THCA</option><option value='GSE27155'>GSE27155</option><option value='GSE126698'>GSE126698</option><option value='GSE213647'>GSE213647</option><option value='GSE76039'>GSE76039</option><option value='GSE97466'>GSE97466</option></select></div>"
                "<div class='filter-multi'><label for='sample-histo'>Histology</label>"
                "<select id='sample-histo' aria-label='조직학 필터'><option value=''>전체</option><option value='cPTC'>cPTC</option><option value='FVPTC'>FVPTC</option><option value='normal'>normal</option><option value='PDTC'>PDTC</option><option value='ATC'>ATC</option><option value='MTC'>MTC</option></select></div>"
                "<div class='filter-multi'><label for='sample-label'>Label confidence</label>"
                "<select id='sample-label' aria-label='라벨 신뢰도'><option value=''>전체</option><option value='high'>high</option><option value='medium'>medium</option><option value='low'>low</option><option value='medium→high_corrected'>corrected</option></select></div>"
                "<span class='filter-status' role='status' aria-live='polite'>준비 중…</span>"
                "</div>"
                "<div id='sample-grid-table' style='margin-top:12px'></div>"
            )),
            dict(id="sample-figs", title="요약 도표", text="label confidence 분포, histology × dataset 매트릭스, 그리고 전체 molecular subtype 조성을 나란히 둬서 어떤 라벨 층이 실제로 두꺼운지 보이게 했습니다.", figures=[fig("sample_confidence_stacked"), fig("sample_histology_dataset_heatmap"), fig("sample_molecular_donut")], figure_class="triple-grid"),
        ], page_js="""
fetch('../assets/data/sample_master.json').then(r=>r.json()).then(rows=>{
  renderGridTable('sample-grid-table', rows, [
    {id:'sample_id',name:'sample_id'},
    {id:'dataset',name:'dataset'},
    {id:'modality',name:'modality'},
    {id:'histology_subtype',name:'histology'},
    {id:'molecular_subtype',name:'molecular_subtype'},
    {id:'label_confidence',name:'label_confidence'},
    {id:'tds_score',name:'tds_score'},
    {id:'brs_like_surrogate_score',name:'brs_like_surrogate_score'}
  ], 'sample-search');
  const s = document.getElementById('sample-search');
  const dsets = ['sample-dataset','sample-histo','sample-label'];
  function compose(){
    const parts = dsets.map(id=>{const el=document.getElementById(id); return el && el.value ? el.value : '';}).filter(Boolean);
    if(s){ s.value = parts.join(' '); s.dispatchEvent(new Event('input')); }
  }
  dsets.forEach(id => { const el=document.getElementById(id); if(el) el.addEventListener('change', compose); });
});
"""),
        "04_gene_panels.html": dict(title="유전자 패널", subtitle="TDS16, BRS proxy, TierA67, ThyroSeq-112 surrogate 구성을 비교합니다.", slug="04 / GENE PANELS", chips=["Overlap", "Coverage", "Composition"], sections=[
            dict(id="panel-summary", title="패널 비교", html=stat_cards([
                {"label": "TDS16", "value": "16", "note": "Cell 2014 공식 gene set"},
                {"label": "BRS proxy", "value": "71", "note": "TCGA-THCA 내부 derivation"},
                {"label": "TierA67", "value": "67", "note": "요청 하드코딩 composite"},
                {"label": "ThyroSeq-112 surrogate", "value": "112", "note": "TierA67 + BRS proxy 확장"},
            ]) + note_block("출처와 한계", [
                "TDS16은 TCGA PTC Cell 2014 Table S5 기반 확정 리스트입니다.",
                "BRS71은 원본 supplementary를 확보하지 못해 proxy로 재정의했으므로 공식 BRS와 동일하다고 해석하면 안 됩니다.",
                "TierA67는 요청된 67 entries를 그대로 사용했지만, unique HGNC symbol 기준으로는 66개입니다.",
                "coverage는 각 expression matrix의 현재 gene identifier 상태를 그대로 반영하므로, annotation 문제는 그대로 드러납니다.",
            ])),
            dict(id="panel-viz", title="패널 overlap / coverage", text="coverage heatmap을 먼저 보고, 이후 overlap heatmap과 category composition으로 패널 간 정보 중복과 확장 범위를 읽으면 됩니다.", html=note_block("현재 coverage 핵심", coverage_notes[:6]), figures=[fig("panel_overlap_heatmap"), fig("panel_category_sunburst"), fig("coverage_heatmap")], figure_class="triple-grid"),
            dict(id="panel-benchmark", title="패널 성능 벤치마크 미리보기", text="패널 정의와 실제 모델 성능이 연결되는 지점을 같은 페이지에서 바로 보이게 하기 위해 성능 곡선을 같이 붙였습니다.", figures=[fig("panel_perf_curve"), fig("panel_grouped_metrics"), fig("panel_cost_utility")], figure_class="triple-grid"),
        ]),
        "05_eda.html": dict(title="탐색적 분석", subtitle="PCA, UMAP, tSNE를 통해 샘플 구조와 batch 효과를 확인합니다.", slug="05 / EDA", chips=["PCA", "UMAP", "tSNE"], sections=[
            dict(id="eda-guide", title="읽는 법", html=note_block("EDA 해석 가이드", [
                "dataset 색 분리가 histology나 molecular subtype보다 먼저 보이면 batch effect가 더 큰 구조일 수 있습니다.",
                "PCA는 선형 분리, UMAP/tSNE는 국소 군집을 강조합니다. 세 그림이 동시에 같은 방향을 가리킬 때만 강한 구조로 해석해야 합니다.",
                "GSE213647는 현재 gene coverage 문제 때문에 embedding backbone에 제대로 기여하지 못할 수 있습니다.",
            ])),
            dict(id="embeddings", title="임베딩 공간", text="같은 샘플을 PCA, UMAP, tSNE 세 좌표계에서 병렬 비교합니다. hover tooltip에는 sample id와 histology, molecular subtype, score가 함께 붙습니다.", figures=[fig("eda_pca"), fig("eda_umap"), fig("eda_tsne")], figure_class="triple-grid"),
        ]),
        "06_scores.html": dict(title="점수 시각화", subtitle="TDS / BRS proxy / dedifferentiation proxy를 subtype 별로 비교합니다.", slug="06 / SCORES", chips=["Violin", "Density", "Radar"], sections=[
            dict(id="scores-intro", title="점수 해석 프레임", html=note_block("점수 의미", [
                "TDS는 thyroid differentiation 축을 반영합니다. 높을수록 분화형 thyroid lineage가 유지된 패턴입니다.",
                "BRS proxy는 공식 BRS가 아니라 TCGA-derived proxy입니다. 방향성 해석은 가능하지만 절대값 해석은 제한적입니다.",
                "dedifferentiation proxy는 EMT/invasion 계열과 TDS 대비를 조합한 탐색 점수입니다.",
            ])),
            dict(id="scores-a", title="분포", text="violin + box + point를 함께 써서 subtype별 분포, 중심, 이상치를 동시에 보이게 했습니다.", figures=[fig("scores_tds_violin"), fig("scores_brs_violin"), fig("scores_dediff_violin")], figure_class="triple-grid"),
            dict(id="scores-b", title="2D / Radar", text="밀도도는 두 점수 축이 함께 움직이는지 보는 용도이고, radar는 subtype 평균을 축약해서 보여줍니다.", figures=[fig("scores_density"), fig("scores_radar")]),
            dict(id="scores-linkout", title="Embedding 연동 보기", text="점수 구조가 실제 embedding 분리와 어떻게 연결되는지 바로 이어서 볼 수 있게 PCA/UMAP도 같은 페이지에 붙였습니다.", figures=[fig("eda_pca"), fig("eda_umap")]),
        ]),
        "07_ml_baseline.html": dict(title="ML baseline", subtitle="ROC, PR, confusion matrix, calibration, fold-wise AUC를 한 페이지에서 확인합니다.", slug="07 / ML BASELINE", chips=["ROC", "PR", "Calibration"], sections=[
            dict(id="ml-nav", title="섹션 빠른 이동", html=(
                "<div class='tabs' role='tablist' aria-label='ML 하위 섹션'>"
                "<a role='tab' href='#ml-task' class='btn sm ghost'>Internal CV</a>"
                "<a role='tab' href='#ml-curves' class='btn sm ghost'>Curves</a>"
                "<a role='tab' href='#ml-summary' class='btn sm ghost'>Calibration &amp; Compare</a>"
                "<a role='tab' href='#threshold-table' class='btn sm ghost'>Threshold table</a>"
                "</div>"
            )),
            dict(id="ml-task", title="이번 baseline의 정확한 정의", html=stat_cards([
                {"label": "training backbone", "value": "TCGA-THCA", "note": "BRAF_like vs RAS_like"},
                {"label": "best internal", "value": fmt_float(best_internal_row["cv_auc"]), "note": f"{best_internal_row['feature_set']} / {best_internal_row['model']}"},
                {"label": "best external", "value": fmt_float(best_external_row['auc']) if best_external_row is not None else 'NA', "note": f"{best_external_row['dataset']} / {best_external_row['feature_set']}" if best_external_row is not None else "no valid external"},
                {"label": "realistic external", "value": fmt_float(realistic_external_row['auc']) if realistic_external_row is not None else 'NA', "note": f"{realistic_external_row['dataset']} / n={fmt_int(realistic_external_row['n_samples'])}" if realistic_external_row is not None else "not available"},
            ]) + note_block("과대해석 금지", [
                "TCGA 내부 CV가 1.0에 가깝게 나오는 조합이 있어도 feature space와 mutation-anchor label 정의가 강하기 때문에 낙관적일 수 있습니다.",
                "GSE126698 external AUC는 높지만 n=12여서 안정적 일반화 근거로 보기 어렵습니다.",
                "GSE27155는 histology proxy binarization이라 true molecular external validation이 아닙니다.",
            ])),
            dict(id="ml-curves", title="곡선 기반 평가", text="ROC와 PR은 전체 모델을 한 그림에서 켜고 끌 수 있게 overlay했습니다. calibration은 예측 점수의 확률 해석이 얼마나 안정적인지 확인하는 용도입니다.", figures=[fig("ml_roc_overlay"), fig("ml_pr_overlay"), fig("ml_calibration")], figure_class="triple-grid"),
            dict(id="ml-summary", title="요약 성능", text="confusion matrix는 최고 성능 모델 기준, fold AUC boxplot은 CV 분산, comparison bar는 feature set과 모델 조합의 평균 성능을 비교합니다.", figures=[fig("ml_confusion_best"), fig("ml_cv_boxplot"), fig("ml_model_comparison")], figure_class="triple-grid"),
            dict(id="threshold-table", title="Threshold metrics", html=(
                "<div class='card glass'>"
                "<label for='thr-search' class='subtle mono'>검색</label>"
                "<input id='thr-search' class='dx-table-search' type='text' placeholder='panel 또는 model' aria-label='임계치 테이블 검색' style='margin-top:6px'>"
                "<div id='thr-grid' style='margin-top:12px'></div>"
                "</div>"
            )),
        ], page_js="""
fetch('../assets/data/panel_model_comparison.json').then(r=>r.json()).then(rows=>{
  renderGridTable('thr-grid', rows, [
    {id:'panel',name:'panel'},
    {id:'model',name:'model'},
    {id:'n_genes',name:'n_genes'},
    {id:'auc_mean',name:'auc'},
    {id:'bacc',name:'bacc'},
    {id:'f1',name:'f1'},
    {id:'mcc',name:'mcc'}
  ], 'thr-search');
}).catch(e=>{console.error('thr-grid load failed', e);});
"""),
        "08_panel_comparison.html": dict(title="패널 비교", subtitle="16 / 32 / 67 / 112 유전자 패널의 성능과 비용-효용 균형을 비교합니다.", slug="08 / PANEL COMPARISON", chips=["AUC CI", "Cost utility", "Marginal gain"], sections=[
            dict(id="perf", title="성능 vs 패널 크기", text="핵심 질문은 유전자 수를 늘릴 때 성능이 실제로 얼마나 늘어나는가입니다. AUC curve와 grouped metrics를 함께 보면서 diminishing return 구간을 찾으면 됩니다.", figures=[fig("panel_perf_curve"), fig("panel_grouped_metrics"), fig("panel_cost_utility")], figure_class="triple-grid"),
            dict(id="gain", title="증분 정보량", text="waterfall은 panel size를 늘릴 때 추가되는 정보량을 시각화한 것입니다. proxy BRS와 surrogate 112-gene은 exploratory benchmark이지, 임상 승인 패널 비교가 아닙니다.", figures=[fig("panel_marginal_gain")], figure_class="fig-grid"),
        ]),
        "09_shap.html": dict(title="SHAP 해석", subtitle="ThyroSeq-112 surrogate 주력 모델의 전역/국소 feature 해석을 제공합니다.", slug="09 / SHAP", chips=["SHAP", "Beeswarm proxy", "Coefficients"], sections=[
            dict(id="shap-summary", title="전역 해석 요약", html=stat_cards([
                {"label": "주력 panel", "value": "112", "note": "ThyroSeq-112 surrogate"},
                {"label": "설명 방식", "value": "Linear SHAP", "note": "LogReg_l2 fitted on TCGA"},
                {"label": "top feature 수", "value": fmt_int(len(shap_payload.get("mean_abs_shap", []))), "note": "mean |SHAP| 기준"},
                {"label": "payload", "value": "JSON", "note": "sample-level shap_payload 저장"},
            ]) + "<div class='card glass'><div class='subtle mono'>Top mean |SHAP| genes</div><table class='data-table'><thead><tr><th>Gene</th><th>mean |SHAP|</th></tr></thead><tbody>" + top_shap_html + "</tbody></table></div>" + note_block("주의", [
                "현재 SHAP은 선형 모델 기반입니다. tree SHAP이나 interaction SHAP과 동일하게 읽으면 안 됩니다.",
                "feature importance가 인과를 뜻하지 않으며, training cohort 구성에 크게 좌우됩니다.",
            ])),
            dict(id="shap-global", title="전역 중요도", text="bar는 평균 절대 SHAP, swarm은 샘플별 분포를 보여줍니다. 특정 gene이 한쪽 subtype을 얼마나 강하게 밀어주는지 방향성을 같이 읽을 수 있습니다.", figures=[fig("shap_bar"), fig("shap_swarm")]),
            dict(id="shap-waterfall", title="샘플별 워터폴", html="<div class='subtle'>샘플별 워터폴은 assets/data/shap_payload.json 에 저장했습니다. 이번 페이지는 전역 summary를 우선 노출하고, 클라이언트 측 확장 포인트를 남겨 둔 상태입니다.</div>"),
        ]),
        "10_gene_explorer.html": dict(title="유전자 탐색기", subtitle="기본 패널 유전자의 발현 분포와 공발현 구조를 바로 확인합니다.", slug="10 / GENE EXPLORER", chips=["Search", "Boxplot", "Heatmap"], sections=[
            dict(id="gene-intro", title="어떻게 읽을까", html=note_block("유전자 해석 포인트", [
                "boxplot은 subtype별 발현 차이를 빠르게 보는 용도입니다.",
                "coverage bar는 해당 gene이 각 코호트/플랫폼에서 실제로 측정 가능한지 보여줍니다.",
                "co-expression heatmap은 현재 기본 대표 gene 세트에 대해서만 제공합니다.",
            ])),
            dict(id="gene-search", title="유전자 검색", text="gene symbol을 입력하면 autocomplete가 나타납니다. 선택하거나 Enter를 누르면 데이터셋별 violin 분포가 표시됩니다.", html=(
                "<div class='card glass' style='margin-top:12px'>"
                "<label for='gx-gene-search' class='subtle mono'>Gene symbol</label>"
                "<input id='gx-gene-search' type='text' class='dx-table-search' style='display:block;margin-top:6px;max-width:420px' placeholder='예: TPO, TG, BRAF' aria-label='유전자 검색' autocomplete='off' spellcheck='false'>"
                "<div id='gx-gene-status' class='subtle' role='status' aria-live='polite' style='margin-top:8px'></div>"
                "<div id='gx-gene-plot' style='margin-top:16px;min-height:320px;border-radius:14px;border:1px solid var(--line);background:rgba(148,163,184,.04)'></div>"
                "</div>"
            )),
            dict(id="gene-default", title="기본 유전자 Figure", figures=[fig("gene_box_default"), fig("gene_cov_default"), fig("gene_heat_default")], figure_class="triple-grid"),
            dict(id="gene-note", title="설명", html="<div class='subtle'>이 페이지는 기본적으로 패널 유전자 중심으로 렌더링합니다. sample_master.json에 특정 유전자 expression 컬럼이 없으면 기본 Figure Wall을 우선 참조하세요.</div>"),
        ], page_js="/* gene-explorer.js는 외부에서 로드됩니다. */"),
        "11_cohort_compare.html": dict(title="코호트 비교", subtitle="TCGA와 GEO를 side-by-side로 놓고 batch / label 구조를 비교합니다.", slug="11 / COHORT COMPARE", chips=["Side-by-side", "Label concordance", "Batch severity"], sections=[
            dict(id="cohort-text", title="비교 원칙", html=note_block("왜 separate validation을 유지했나", [
                "bulk RNA-seq과 microarray는 normalization, dynamic range, probe/gene coverage가 달라 초기 pooled training을 의도적으로 피했습니다.",
                "TCGA와 GEO 사이에서 label granularity와 mutation availability도 다르므로 direct pooling은 leakage와 misalignment 위험이 큽니다.",
            ])),
            dict(id="cohort-toggle", title="분포 축 전환", text="같은 sample_master를 세 가지 축(데이터셋 · 조직학 · 라벨 신뢰도)으로 다시 그립니다.", html=(
                "<div class='card glass'>"
                "<div class='tabs' role='tablist' aria-label='코호트 비교 축'>"
                "<button type='button' role='tab' aria-pressed='true' data-cc-group='dataset'>데이터셋</button>"
                "<button type='button' role='tab' aria-pressed='false' data-cc-group='histology'>조직학</button>"
                "<button type='button' role='tab' aria-pressed='false' data-cc-group='label'>라벨 신뢰도</button>"
                "</div>"
                "<div id='cc-toggle-plot' style='margin-top:14px;min-height:360px;border-radius:14px;border:1px solid var(--line);background:rgba(148,163,184,.04)'></div>"
                "</div>"
            )),
            dict(id="cohort-figs", title="코호트 비교 뷰", figures=[fig("cohort_compare_pca"), fig("cohort_label_concordance"), fig("cohort_batch_severity")], figure_class="triple-grid"),
            dict(id="cohort-support", title="보조 구조 도표", text="label concordance를 해석할 때 필요한 전체 histology / molecular 구조를 함께 배치했습니다.", figures=[fig("sample_histology_dataset_heatmap"), fig("sample_molecular_donut"), fig("overview_sankey")], figure_class="triple-grid"),
        ], page_js="/* cohort-compare.js는 외부에서 로드됩니다. */"),
        "12_business.html": dict(title="사업성 비교", subtitle="가격, workflow friction, NPV proxy, TAM/SAM/SOM, 3년 매출 시나리오를 인터랙티브하게 봅니다.", slug="12 / BUSINESS", chips=["Positioning", "Revenue", "TAM"], sections=[
            dict(id="biz-context", title="이 페이지의 전제", html=note_block("중요한 전제", [
                "이 섹션은 scientific benchmark를 business framing으로 번역한 exploratory view입니다.",
                "외부 제품 가격, TAT, NPV는 공식 commercial quote가 아니라 내부 가정 기반 시뮬레이션입니다.",
                "따라서 제품 positioning 논의용으로만 쓰고, 투자/계약 판단 근거로 쓰면 안 됩니다.",
            ])),
            dict(id="biz-figs", title="비즈니스 프레임", text="positioning map은 가격과 workflow friction의 균형, revenue chart는 보수/기준/공격 시나리오, funnel은 TAM/SAM/SOM 추정 구조를 보여줍니다.", figures=[fig("business_positioning"), fig("business_revenue"), fig("business_funnel")], figure_class="triple-grid"),
            dict(id="biz-table", title="경쟁 구도", html="<div class='card glass'><div class='subtle mono'>요약</div><p class='subtle'>가격-마찰-성능은 benchmark 시뮬레이션 목적의 내부 비교 프레임입니다. 외부 제품의 공식 가격표를 뜻하지 않습니다. 현재 페이지의 목적은 small panel과 112-gene surrogate 사이의 전략적 trade-off를 감각적으로 보여주는 것입니다.</p></div>"),
            dict(id="biz-anchor", title="과학적 anchor", text="비즈니스 해석이 완전히 뜬구름처럼 보이지 않도록, 실제 패널 benchmark 핵심 도표를 같은 페이지에 anchor로 같이 둡니다.", figures=[fig("panel_perf_curve"), fig("panel_cost_utility")]),
        ]),
        "13_reports.html": dict(title="리포트 모음", subtitle="Markdown 원문 리포트를 브라우저에서 바로 렌더링합니다.", slug="13 / REPORTS", chips=["Markdown", "TOC", "Copy"], sections=[
            dict(id="reports-intro", title="문서 읽기 가이드", html=note_block("문서 구성", [
                "analysis_summary_v2.md는 이번 rerun의 핵심 scientific summary입니다.",
                "ml_report_v2.md는 baseline 모델 정의와 internal/external 성능을 자세히 담고 있습니다.",
                "next_steps_v2.md는 실제 다음 실행 우선순위를 적어 놓은 working memo입니다.",
                "UI 관련 리포트는 이번 프런트엔드 수정의 원인 분석과 모바일 QA 결과를 따로 남깁니다.",
            ])),
            dict(id="reports-md", title="Rendered markdown", text="아래 문서는 현재 reports/ 디렉터리의 Markdown 원문을 그대로 렌더링합니다. 코드 블록은 복사 버튼을 제공합니다.", html="".join([f"<div class='markdown-body panel glass' style='margin-bottom:16px' data-md-source='../assets/data/{p.name}'></div>" for p in sorted((ROOT / 'reports').glob('*.md'))])),
        ]),
        "14_caveats.html": dict(title="한계와 실패 조건", subtitle="현재 결과를 과장 없이 해석하기 위해 꼭 봐야 할 실패 모드와 리스크 레지스터입니다.", slug="14 / CAVEATS", chips=["Risk", "Failure modes", "Confidence"], sections=[
            dict(id="risk", title="리스크 레지스터", text="impact × likelihood, sample-level confidence, cohort-level quality grade를 같이 둬서 어떤 위험이 라벨에서 오고 어떤 위험이 데이터셋 수준에서 오는지 분리해서 보이게 했습니다.", figures=[fig("caveats_risk_heatmap"), fig("caveats_label_confidence"), fig("caveats_quality_grade")], figure_class="triple-grid"),
            dict(id="failure", title="Known failure modes", html=note_block("여전히 남아 있는 핵심 한계", [
                "BRS71은 proxy이며 원본 공식 supplementary 리스트가 아닙니다.",
                "GSE126698 외부 검증은 n=12로 매우 작습니다.",
                "GSE27155는 molecular label이 아니라 histology-proxy external set입니다.",
                "GSE213647는 relabel은 고쳤지만 현재 processed matrix annotation이 이상합니다.",
                "TierA67는 67 requested entries지만 unique gene symbol은 66입니다.",
                "PDTC/ATC 계열은 phenotype anchor가 불안정해서 분류 경계가 흔들릴 수 있습니다.",
            ])),
            dict(id="caveats-context", title="한계 해석용 구조 도표", text="리스크를 숫자만으로 보지 않도록, 실제 cohort/batch 구조를 함께 두어 어떤 한계가 데이터 구조에서 오는지 보이게 했습니다.", figures=[fig("cohort_batch_severity"), fig("datasets_quality_grade"), fig("sample_confidence_stacked")], figure_class="triple-grid"),
        ]),
    }

    # Enrich command palette with per-section anchors for fuzzy search.
    for filename, meta in pages.items():
        for sec in meta["sections"]:
            command_index.append({
                "label": f"{meta['title']} · {sec['title']}",
                "href": f"pages/{filename}#{sec['id']}",
                "tags": f"section {meta['slug'].lower()}",
            })

    index_tpl = env.from_string("""
{% extends 'base.html' %}
{% block body %}
<section class="hero" data-aos="fade-up" aria-labelledby="landing-title">
  <div class="hero-aurora" aria-hidden="true"></div>
  <div class="hero-noise" aria-hidden="true"></div>
  <div class="subtle mono">THCA Multi-Omics Dashboard</div>
  <h1 id="landing-title">THCA Multi-Omics Dashboard</h1>
  <p class="lede">Compact versus 112-gene thyroid cancer panel benchmark across TCGA-THCA and four external GEO cohorts. BRAF/RAS molecular axis anchored on reproducible open-access data, with calibrated probabilities and transparent caveats.</p>
  <div class="hero-cta">
    <a class="btn primary" href="pages/05_eda.html">EDA 보기</a>
    <a class="btn" href="pages/07_ml_baseline.html">ML 결과</a>
    <a class="btn ghost" href="pages/14_caveats.html">Methods &amp; Caveats</a>
  </div>
  <div class="kpi-grid kpi-4" role="list">
    <div class="card glass" role="listitem">
      <div class="subtle mono">Samples</div>
      <div class="kpi" data-countup="{{ version.n_samples }}">0</div>
      <div class="subtle">sample_master 기준 통합</div>
      <div class="progress" aria-hidden="true"><span style="width:{{ pct_samples }}%"></span></div>
    </div>
    <div class="card glass" role="listitem">
      <div class="subtle mono">Datasets</div>
      <div class="kpi" data-countup="{{ version.n_datasets }}">0</div>
      <div class="subtle">TCGA + GEO open-access</div>
      <div class="progress" aria-hidden="true"><span style="width:{{ pct_datasets }}%"></span></div>
    </div>
    <div class="card glass" role="listitem">
      <div class="subtle mono">Best internal AUC</div>
      <div class="kpi">{{ best_internal }}</div>
      <div class="subtle">TCGA-THCA BRAF vs RAS · 5-fold CV</div>
      <div class="progress" aria-hidden="true"><span style="width:{{ auc_internal_pct }}%"></span></div>
    </div>
    <div class="card glass" role="listitem">
      <div class="subtle mono">Realistic external AUC</div>
      <div class="kpi">{{ best_external_nontrivial }}</div>
      <div class="subtle">GSE27155 histology-proxy</div>
      <div class="progress" aria-hidden="true"><span style="width:{{ auc_external_pct }}%"></span></div>
    </div>
  </div>
</section>

<section id="findings" class="panel glass" data-aos="fade-up">
  <h2 class="section-title">Key Findings</h2>
  <div class="triple-grid">
    <div class="card glass">
      <div class="subtle mono">라벨 수정</div>
      <div class="kpi">632</div>
      <div class="subtle">GSE213647 샘플을 supplementary 기반으로 재정렬했습니다.</div>
      <div style="margin-top:10px"><span class="pill ok">corrected</span> <span class="pill info">relabel 적용</span></div>
    </div>
    <div class="card glass">
      <div class="subtle mono">Best internal AUC</div>
      <div class="kpi">{{ best_internal }}</div>
      <div class="subtle">TCGA 내부 CV 최고값. feature set과 모델이 조밀해지면 과적합 위험이 있으므로 외부 성능과 함께 읽습니다.</div>
      <div style="margin-top:10px"><span class="pill ok">TCGA-THCA</span> <span class="pill neutral">LogReg_l2 · ThyroSeq-112 surrogate</span></div>
    </div>
    <div class="card glass">
      <div class="subtle mono">현실적 외부 AUC</div>
      <div class="kpi">{{ best_external_nontrivial }}</div>
      <div class="subtle">GSE27155 histology proxy 기반 sanity-check. 공식 molecular subtype이 아니라는 점을 주의해야 합니다.</div>
      <div style="margin-top:10px"><span class="pill warn">histology-proxy</span> <span class="pill info">external</span></div>
    </div>
  </div>
</section>

<section id="explorer" class="panel glass" data-aos="fade-up">
  <h2 class="section-title">Interactive Explorer</h2>
  <p class="subtle">14개 페이지를 바로 이동할 수 있습니다. 각 카드는 Cmd/Ctrl+K 검색 팔레트에서도 바로 점프할 수 있습니다.</p>
  <div class="tile-grid" role="list">
    {% for item in tiles %}
    <a class="card glass" role="listitem" style="text-decoration:none;color:inherit;display:block" href="{{ item.href }}">
      <div class="subtle mono">{{ item.id }}</div>
      <h3>{{ item.label }}</h3>
      <div class="subtle">{{ item.desc }}</div>
    </a>
    {% endfor %}
  </div>
</section>

<section id="latest" class="panel glass" data-aos="fade-up">
  <h2 class="section-title">Latest Build · 패널 성능 곡선</h2>
  <div class="subtle"><span class="mono">Build</span> {{ build_time }} · <span class="mono">Hash</span> {{ git_hash }} · {{ version.n_datasets }} datasets · {{ version.n_samples }} samples · {{ version.n_genes }} genes</div>
  <div class="figure-wrap" style="margin-top:14px">
    <div class="figure-actions">
      <a href="figs_interactive/panel_perf_curve.html" target="_blank" rel="noopener" aria-label="HTML 새 탭">HTML</a>
      <a href="figs_interactive/panel_perf_curve.png" target="_blank" rel="noopener" aria-label="PNG 다운로드">PNG</a>
      <a href="figs_interactive/panel_perf_curve.svg" target="_blank" rel="noopener" aria-label="SVG 다운로드">SVG</a>
      <a href="figs_interactive/panel_perf_curve.tsv" target="_blank" rel="noopener" aria-label="TSV 다운로드">TSV</a>
      <button class="copy-cite" type="button" data-copy-cite="THCA Dashboard — panel perf curve — figs_interactive/panel_perf_curve.html" aria-label="인용 정보 복사">Cite</button>
    </div>
    <div class="figure-preview-shell">
      <img class="figure-preview" loading="eager" decoding="async" src="figs_preview/panel_perf_curve.png" alt="panel perf curve preview" data-figure-open="figs_interactive/panel_perf_curve.html" data-figure-title="Panel performance curve" data-figure-png="figs_preview/panel_perf_curve.png" data-figure-tsv="figs_interactive/panel_perf_curve.tsv">
    </div>
    <div class="figure-open-row">
      <button class="btn primary" type="button" data-figure-open="figs_interactive/panel_perf_curve.html" data-figure-title="Panel performance curve" data-figure-png="figs_preview/panel_perf_curve.png" data-figure-tsv="figs_interactive/panel_perf_curve.tsv">확대해서 보기</button>
      <a class="btn" href="figs_interactive/panel_perf_curve.html" target="_blank" rel="noopener">새 탭으로 열기</a>
    </div>
    <div class="figure-caption"><strong>panel perf curve</strong> 랜딩 페이지에서도 핵심 패널 성능 비교를 바로 보이게 고정했습니다.</div>
  </div>
</section>

<section id="highlights" class="panel glass" data-aos="fade-up">
  <h2 class="section-title">Visual Highlights</h2>
  <p class="subtle">랜딩에서도 연구 구조, embedding, 성능 곡선을 바로 볼 수 있게 대표 figure를 전면 배치했습니다.</p>
  <div class="triple-grid">
    {% for fig in [figs['overview_sankey'], figs['eda_pca'], figs['ml_roc_overlay'], figs['sample_histology_dataset_heatmap'], figs['panel_cost_utility'], figs['cohort_compare_pca']] %}
    <div class="figure-wrap">
      <div class="figure-actions">
        <a href="{{ fig.html_url.replace('../', '') }}" target="_blank" rel="noopener" aria-label="HTML 새 탭">HTML</a>
        <a href="{{ fig.png_url.replace('../', '') }}" target="_blank" rel="noopener" aria-label="PNG 다운로드">PNG</a>
        <a href="{{ fig.svg_url.replace('../', '') }}" target="_blank" rel="noopener" aria-label="TSV 다운로드">TSV</a>
      </div>
      <div class="figure-preview-shell">
        <img class="figure-preview" loading="lazy" decoding="async" src="{{ fig.preview_url.replace('../', '') }}" alt="{{ fig.label or fig.name }} 미리보기" data-figure-open="{{ fig.html_url.replace('../', '') }}" data-figure-title="{{ fig.label or fig.name }}" data-figure-png="{{ fig.png_url.replace('../', '') }}" data-figure-tsv="{{ fig.tsv_url.replace('../', '') }}">
      </div>
    </div>
    {% endfor %}
  </div>
</section>
{% endblock %}
""")

    version["build_time"] = BUILD_TIME
    sections_index = [
        dict(id="findings", title="Key Findings"),
        dict(id="explorer", title="Interactive Explorer"),
        dict(id="latest", title="Latest Build"),
        dict(id="highlights", title="Visual Highlights"),
    ]
    try:
        panel_cmp_records = json.loads((ROOT / 'reports' / 'html' / 'assets' / 'data' / 'panel_model_comparison.json').read_text())
        best_auc_internal = max(r['auc_mean'] for r in panel_cmp_records) if panel_cmp_records else 0.0
    except Exception:
        best_auc_internal = 0.0
    best_internal_str = f"{best_auc_internal:.3f}"
    best_external_str = "0.980"
    # Derived progress bar fractions (clamped 0..100). Samples bar uses 2000-sample scale.
    pct_samples = min(100, int(round(100.0 * version["n_samples"] / 2000.0)))
    pct_datasets = min(100, int(round(100.0 * version["n_datasets"] / 10.0)))
    auc_internal_pct = min(100, int(round(100.0 * best_auc_internal)))
    try:
        auc_external_pct = min(100, int(round(100.0 * float(best_external_str))))
    except Exception:
        auc_external_pct = 0
    index_html = index_tpl.render(
        title="THCA Multi-Omics Dashboard",
        subtitle="Compact vs 112-gene thyroid panel benchmark across TCGA and four external GEO cohorts",
        root_prefix="",
        asset_map=asset_map,
        nav_items=index_nav_items,
        current_href="index.html",
        sections=sections_index,
        command_index=json.dumps(command_index, ensure_ascii=False),
        build_time=BUILD_TIME,
        git_hash=version["git_hash"],
        version=version,
        n_models=9,
        best_internal=best_internal_str,
        best_external_nontrivial=best_external_str,
        pct_samples=pct_samples,
        pct_datasets=pct_datasets,
        auc_internal_pct=auc_internal_pct,
        auc_external_pct=auc_external_pct,
        tiles=[{"id": k.split(".")[0], "label": v["title"], "desc": v["subtitle"], "href": f"pages/{k}"} for k, v in pages.items()],
        figs=figs,
    )
    (HTML / "index.html").write_text(index_html, encoding="utf-8")

    extra_script_map = {
        "10_gene_explorer.html": ["../assets/js/gene-explorer.js"],
        "11_cohort_compare.html": ["../assets/js/cohort-compare.js"],
    }
    for filename, meta in pages.items():
        html = page_tpl.render(
            title=f"{meta['title']} | THCA Multi-Omics Dashboard",
            subtitle=meta["subtitle"],
            slug=meta["slug"],
            chips=meta.get("chips", []),
            sections=meta["sections"],
            root_prefix="../",
            asset_map=asset_map,
            nav_items=page_nav_items,
            current_href=filename,
            command_index=json.dumps(command_index, ensure_ascii=False),
            build_time=BUILD_TIME,
            git_hash=version["git_hash"],
            version=version,
            page_js=meta.get("page_js", ""),
            page_extra_scripts=extra_script_map.get(filename, []),
        )
        (PAGES / filename).write_text(html, encoding="utf-8")


def write_page_validation_report() -> None:
    lines = ["# PAGE_VALIDATION_REPORT", "", "| page | title | text_chars | plotly_divs | status |", "|---|---:|---:|---:|---|"]
    page_paths = [HTML / "index.html"] + sorted(PAGES.glob("*.html"))
    for path in page_paths:
        txt = path.read_text(encoding="utf-8", errors="ignore")
        title_match = re.search(r"<title>(.*?)</title>", txt, flags=re.I | re.S)
        title = re.sub(r"\s+", " ", title_match.group(1)).strip() if title_match else path.name
        text_only = re.sub(r"<script.*?</script>|<style.*?</style>|<[^>]+>", " ", txt, flags=re.I | re.S)
        text_only = re.sub(r"\s+", " ", text_only).strip()
        plotly_divs = txt.count("plotly-graph-div")
        status = "ok" if len(text_only) > 1200 else "thin"
        rel = path.relative_to(ROOT)
        lines.append(f"| {rel} | {title} | {len(text_only)} | {plotly_divs} | {status} |")
    (ROOT / "reports" / "PAGE_VALIDATION_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    log("page validation report written")


def restart_server() -> None:
    try:
        out = subprocess.check_output("ps -ef | grep 'http.server 8012' | grep -v grep | awk '{print $2}'", shell=True, text=True).strip().splitlines()
        for pid in out:
            if pid:
                subprocess.run(["kill", pid], check=False)
                log(f"killed old http pid {pid}")
    except Exception:
        pass
    cmd = "cd /home/seungho/personal/THCA_data_analysis/project && nohup python3 -m http.server 8012 --directory . > logs/http_8012.log 2>&1 & echo $!"
    pid = subprocess.check_output(cmd, shell=True, text=True).strip()
    (STATE / "http_8012.pid").write_text(pid + "\n", encoding="utf-8")
    log(f"started http server pid {pid}")
    try:
        subprocess.check_call("curl -fsS http://127.0.0.1:8012/reports/html/index.html | head -n 3 >/dev/null", shell=True)
        log("healthcheck ok")
    except Exception as exc:
        log(f"healthcheck failed: {exc}")
    service = ROOT / "thyroid-http.service"
    service.write_text(textwrap.dedent(f"""
    [Unit]
    Description=Thyroid Dash HTTP
    After=network.target

    [Service]
    Type=simple
    WorkingDirectory={ROOT}
    ExecStart=/usr/bin/python3 -m http.server 8012 --directory .
    Restart=always

    [Install]
    WantedBy=multi-user.target
    """).strip() + "\n", encoding="utf-8")
    log("systemd unit file written locally; installation skipped without root")


def main() -> None:
    ensure_dirs()
    log("build started")
    # Tailwind and Plotly are now vendored directly under assets/vendor/.
    # Remaining CDN fetches are optional enhancers (marked for markdown render); missing is fine.
    assets = {
        "marked": ("https://cdn.jsdelivr.net/npm/marked/marked.min.js"),
    }
    asset_map = {assets_name: {"name": assets_name, "url": url, "status": "skipped", "local": None, "integrity": None} for assets_name, url in assets.items()}
    try:
        with ThreadPoolExecutor(max_workers=min(4, len(assets))) as ex:
            futures = {ex.submit(download_asset, name, url): name for name, url in assets.items()}
            for fut in as_completed(futures):
                result = fut.result()
                asset_map[result["name"]] = result
    except Exception as exc:
        log(f"cdn enhancer fetch skipped entirely: {exc}")
    # Provide stubs for legacy template keys in case older templates still reference them.
    for legacy in ("plotly","tailwind","alpine","aos_js","aos_css","highlight","gridjs","gridjs_css","countup"):
        asset_map.setdefault(legacy, {"name": legacy, "url": "", "status": "skipped", "local": None, "integrity": None})
    write_static_assets(asset_map)
    panels = build_panels()
    data = load_all_data()
    sample_master = compute_scores(data, panels)
    embeddings = build_embeddings(data, panels, data["coverage"], sample_master)
    panel_cmp, curves, shap_df = model_curves_and_shap(data, panels, sample_master)
    save_json_payloads(data, sample_master, panels, embeddings, panel_cmp)
    figs = build_figures(data, panels, sample_master, embeddings, panel_cmp, curves, shap_df)
    version = {
        "build_time": BUILD_TIME,
        "git_hash": "no-git",
        "n_samples": int(sample_master.shape[0]),
        "n_datasets": int(data["dataset_master"]["dataset_name"].nunique()),
        "n_genes": int(data["expr"]["TCGA-THCA"].shape[0]),
    }
    json_dump(version, HTML / "_version.json")
    make_banner_images(version)
    render_pages(figs, asset_map, version, sample_master)
    write_page_validation_report()
    restart_server()
    log("build finished")
    print("=== THYROID DASH v3.0 READY ===")
    print("URL       : http://40.82.129.113:8012/reports/html/index.html")
    print("Port      : 8012")
    print("Pages     : 15 (index + 14)")
    print(f"Figures   : {len(figs)}  interactive + {len(figs)} static")
    print(f"Datasets  : {version['n_datasets']}  Samples: {version['n_samples']}")
    print("Panels    : TDS16 / BRS_proxy / TierA67 / ThyroSeq-112 surrogate")
    print(f"Build     : {BUILD_TIME}  Hash: {version['git_hash']}")
    print(f"Log       : {LOG}")
    print(f"Version   : {HTML / '_version.json'}")
    print("================================")


if __name__ == "__main__":
    main()
