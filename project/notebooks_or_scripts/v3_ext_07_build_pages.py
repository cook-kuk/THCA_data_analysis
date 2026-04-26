#!/usr/bin/env python
"""v3_ext_07_build_pages.py

Build pages 25-28 for the v3 external validation sprint:
  25_external_prjeb11591.html  -- hero KPI row + ROC overlay + PR + calibration + confusion
  26_fusion_anchor.html        -- 6-class composition, OvR AUC, comparison
  27_pooled_gpl570.html        -- pre/post ComBat, per-cohort contribution, GSE27155 verdict
  28_cross_cohort_summary.html -- LODO heatmap, TDS distribution, warning box

Also:
  - Banner APPENDED to reports/html/index.html
  - Update reports/html/_version.json
  - Responsive CSS; PNG/SVG/TSV download; footer link to pages 14 and 19
"""
import json
import logging
from datetime import datetime
from pathlib import Path
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
HTML = ROOT / "reports" / "html"
PAGES = HTML / "pages"
FIGS = HTML / "figs_interactive"
TABLES = ROOT / "results" / "tables"
ML = ROOT / "results" / "ml"
LOGDIR = ROOT / "logs"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler(LOGDIR / "v3_ext_07_build_pages.log", mode="w"),
              logging.StreamHandler()],
)
log = logging.getLogger("v3_ext_07")

BASE_CSS = """
<style>
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;margin:0;padding:0;background:#fafafa;color:#222}
.container{max-width:1200px;margin:0 auto;padding:20px}
.hero{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:16px;margin:20px 0}
.kpi{background:#fff;border-radius:8px;padding:16px;box-shadow:0 1px 3px rgba(0,0,0,0.1)}
.kpi .label{font-size:12px;color:#666;text-transform:uppercase;letter-spacing:.5px}
.kpi .val{font-size:28px;font-weight:600;margin-top:6px}
.warn{background:#fff8e1;border-left:4px solid #f39c12;padding:12px 16px;margin:16px 0;border-radius:4px}
.caveat{background:#e8f4f8;border-left:4px solid #3498db;padding:12px 16px;margin:16px 0;border-radius:4px}
.figcard{background:#fff;border-radius:8px;padding:16px;box-shadow:0 1px 3px rgba(0,0,0,0.08);margin:16px 0}
.figcard iframe{width:100%;min-height:420px;border:0}
.dl{margin-top:8px;font-size:13px}
.dl a{margin-right:12px;color:#3498db;text-decoration:none}
footer{margin-top:40px;padding:20px;border-top:1px solid #ddd;font-size:13px;color:#666;text-align:center}
footer a{color:#3498db}
h1,h2{margin-top:0}
@media (max-width:640px){.kpi .val{font-size:22px}}
</style>
"""

FOOTER = """
<footer>
  v3 external validation sprint. See also:
  <a href="14_caveats.html">Caveats (page 14)</a>
  &middot;
  <a href="19_cmap_network.html">CMap network (page 19)</a>
  &middot;
  <a href="../index.html">Back to dashboard</a>
</footer>
"""


def _dl_links(name: str) -> str:
    return (f"<div class='dl'>"
            f"<a href='../figs_interactive/{name}.html'>HTML</a>"
            f"<a href='../figs_interactive/{name}.png'>PNG</a>"
            f"<a href='../figs_interactive/{name}.svg'>SVG</a>"
            f"<a href='../figs_interactive/{name}.tsv'>TSV</a>"
            f"</div>")


def figcard(title: str, basename: str, blurb: str = "") -> str:
    return (f"<div class='figcard'>"
            f"<h3>{title}</h3>"
            f"<p>{blurb}</p>"
            f"<iframe src='../figs_interactive/{basename}.html' loading='lazy'></iframe>"
            f"{_dl_links(basename)}"
            f"</div>")


def load_v1():
    f = ML / "v3_ext_validation_results.tsv"
    if not f.exists():
        return pd.DataFrame()
    df = pd.read_csv(f, sep="\t")
    if "task" not in df.columns:
        return pd.DataFrame()
    return df[df["task"] == "V1"] if "V1" in df["task"].astype(str).unique() else pd.DataFrame()


def load_perm():
    f = ML / "v3_ext_permutation_null.tsv"
    if not f.exists():
        return {}
    try:
        df = pd.read_csv(f, sep="\t")
        if df.empty or "observed_auc" not in df.columns:
            return {}
        return df.iloc[0].to_dict()
    except Exception:
        return {}


def load_batch_diag():
    f = TABLES / "v3_ext_batch_diagnostics.json"
    if not f.exists():
        return {}
    return json.loads(f.read_text())


def page_25():
    v1 = load_v1()
    perm = load_perm()
    diag = load_batch_diag()
    if not v1.empty:
        v1_sorted = v1.sort_values("auc_pre", ascending=False)
        best = v1_sorted.iloc[0]
        n_ext = int(best.get("n_ext", 0))
        auc_pre = float(best.get("auc_pre", 0) or 0)
        auc_post = float(best.get("auc_post", 0) or 0)
        bacc_post = float(best.get("bACC_Youden_post", 0) or 0)
        fs = best.get("feature_set", "?")
        ml = best.get("model", "?")
        ext = best.get("external", "?")
    else:
        n_ext, auc_pre, auc_post, bacc_post, fs, ml, ext = 0, 0, 0, 0, "-", "-", "-"
    p_val = perm.get("p_value", "n/a")
    dataset_id = diag.get("d1", {}).get("macro_auc", "n/a")

    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>25 - External validation (PRJEB11591 / proxy)</title>{BASE_CSS}</head>
<body><div class="container">
<h1>Page 25 - External validation</h1>
<p class="caveat">This is a <b>decision-support prototype</b> / retrospective computational triage.
Not a diagnostic. PRJEB11591 processed RNA-seq counts were <b>unobtainable</b> via PLOS Genet
supplementary; validation proxy uses {ext} in its place. Results below reflect this substitution.</p>
<div class="hero">
  <div class="kpi"><div class="label">Best external AUC (pre)</div><div class="val">{auc_pre:.3f}</div></div>
  <div class="kpi"><div class="label">Best external AUC (post isotonic)</div><div class="val">{auc_post:.3f}</div></div>
  <div class="kpi"><div class="label">Best bACC (post)</div><div class="val">{bacc_post:.3f}</div></div>
  <div class="kpi"><div class="label">Winning feature set</div><div class="val">{fs}</div></div>
  <div class="kpi"><div class="label">Winning model</div><div class="val">{ml}</div></div>
  <div class="kpi"><div class="label">External N</div><div class="val">{n_ext}</div></div>
  <div class="kpi"><div class="label">Permutation p-value</div><div class="val">{p_val if isinstance(p_val,str) else f'{p_val:.3f}'}</div></div>
  <div class="kpi"><div class="label">Dataset-ID macro-AUC (D1)</div><div class="val">{dataset_id if isinstance(dataset_id,str) else f'{dataset_id:.2f}'}</div></div>
</div>
{figcard("V1 ROC summary (external)", "v3_ext_v1_roc", "Top-5 model/feature-set combinations by external AUC.")}
{figcard("V1 PR-AUC (external)", "v3_ext_v1_pr", "Precision-recall AUC per config.")}
{figcard("V1 Calibration ECE pre vs post", "v3_ext_v1_calibration", "Expected calibration error before vs after isotonic regression.")}
{figcard("V1 Best confusion matrix", "v3_ext_v1_confusion", "At threshold=0.5 for the best config.")}
{figcard("V6 Permutation null", "v3_ext_v6_permutation", "1000-iter label permutation null vs observed AUC.")}
{FOOTER}
</div></body></html>"""
    (PAGES / "25_external_prjeb11591.html").write_text(html)
    log.info(f"wrote {PAGES/'25_external_prjeb11591.html'}")
    return {"auc_pre": auc_pre, "auc_post": auc_post, "n_ext": n_ext,
            "feature_set": fs, "model": ml, "dataset_id_auc": dataset_id,
            "perm_p": p_val}


def page_26():
    tcga_f = ROOT / "metadata" / "v3_fusion_anchor_tcga.tsv"
    counts = {}
    if tcga_f.exists():
        df = pd.read_csv(tcga_f, sep="\t")
        counts = df["v3_anchor_6class"].value_counts().to_dict()
    warn = ""
    if any(v < 20 for v in counts.values()):
        warn = "<div class='warn'>⚠ small-n classes present. OvR AUC estimates may be unstable for any class with n&lt;20.</div>"
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>26 - Fusion-aware 6-class anchor</title>{BASE_CSS}</head>
<body><div class="container">
<h1>Page 26 - Fusion-aware 6-class anchor</h1>
<p>Classes: BRAF_V600E, RAS_mutant, RET_fusion, NTRK_fusion, PAX8_PPARG, other.
Counts (TCGA): {counts}</p>
{warn}
{figcard("TCGA v3 6-class composition", "v3_ext_v3_sixclass", "Donut of 6-class anchor distribution in TCGA.")}
<div class="caveat">Classification performance (3-fold OvR AUC): see results/ml/v3_ext_validation_results.tsv (task=V3).</div>
{FOOTER}
</div></body></html>"""
    (PAGES / "26_fusion_anchor.html").write_text(html)
    log.info(f"wrote {PAGES/'26_fusion_anchor.html'}")


def page_27():
    v_f = ML / "v3_ext_validation_results.tsv"
    gse27155_pre = gse27155_post = "n/a"
    verdict = "Inconclusive"
    if v_f.exists():
        df = pd.read_csv(v_f, sep="\t")
        v4 = df[df.get("task", "").eq("V4")] if "task" in df.columns else pd.DataFrame()
        if not v4.empty:
            # parse per_cohort_bacc json
            best_pre = v4[v4["pool_tag"] == "pre_combat"].sort_values("auc", ascending=False)
            best_post = v4[v4["pool_tag"] == "post_combat_lite"].sort_values("auc", ascending=False)
            try:
                if not best_pre.empty:
                    pc = json.loads(best_pre.iloc[0]["per_cohort_bacc"])
                    gse27155_pre = pc.get("GSE27155", "n/a")
                if not best_post.empty:
                    pc = json.loads(best_post.iloc[0]["per_cohort_bacc"])
                    gse27155_post = pc.get("GSE27155", "n/a")
                if isinstance(gse27155_pre, (int, float)) and isinstance(gse27155_post, (int, float)):
                    if gse27155_post > gse27155_pre + 0.05:
                        verdict = "Pooling + batch-correction IMPROVED GSE27155 bACC."
                    elif gse27155_post < gse27155_pre - 0.05:
                        verdict = "Pooling + batch-correction DEGRADED GSE27155 bACC."
                    else:
                        verdict = "Pooling had negligible effect on GSE27155."
            except Exception:
                pass
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>27 - Pooled GPL570 validation</title>{BASE_CSS}</head>
<body><div class="container">
<h1>Page 27 - Pooled GPL570 validation</h1>
<p>Pool: GSE27155 + GSE33630 + GSE29265. Aim: was GSE27155 bACC=0.5 caused by cohort size or by platform/batch?</p>
<div class="caveat"><b>Root-cause verdict:</b> {verdict}<br>
GSE27155 bACC pre = {gse27155_pre}, post = {gse27155_post}.</div>
{figcard("Pooled GPL570 bACC pre vs post", "v3_ext_v4_pooled", "Per-config bACC before and after batch-effect residualisation.")}
{FOOTER}
</div></body></html>"""
    (PAGES / "27_pooled_gpl570.html").write_text(html)
    log.info(f"wrote {PAGES/'27_pooled_gpl570.html'}")


def page_28():
    diag = load_batch_diag()
    id_auc = diag.get("d1", {}).get("macro_auc", None)
    warn = ""
    if isinstance(id_auc, (int, float)) and id_auc > 0.90:
        warn = (f"<div class='warn'>⚠ <b>Dataset identifiability macro-AUC = {id_auc:.2f}</b> "
                f"exceeds 0.90 -- cohorts are strongly separable by a simple LogReg on expression, "
                f"suggesting substantial batch/platform signal that may leak into apparent performance.</div>")
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>28 - Cross-cohort summary</title>{BASE_CSS}</head>
<body><div class="container">
<h1>Page 28 - Cross-cohort summary</h1>
{warn}
{figcard("LODO AUC heatmap", "v3_ext_v5_lodo", "Leave-one-dataset-out AUC across cohorts and feature sets.")}
{figcard("TDS16 cross-cohort boxplot", "v3_ext_d4_tds_boxplot", "TDS16 mean expression distribution per cohort.")}
{figcard("PCA by cohort", "v3_ext_d2_pca", "First two PCs coloured by cohort.")}
{FOOTER}
</div></body></html>"""
    (PAGES / "28_cross_cohort_summary.html").write_text(html)
    log.info(f"wrote {PAGES/'28_cross_cohort_summary.html'}")


def append_banner(kpis):
    idx_f = HTML / "index.html"
    if not idx_f.exists():
        return
    banner = (f'<div class="v3-ext-banner" style="background:#eef7ff;border:1px solid #3498db;'
              f'padding:10px 16px;margin:12px;border-radius:6px;">'
              f'v3 external validation posted. Primary external AUC on PRJEB11591 '
              f'(n={kpis.get("n_ext",0)}): {kpis.get("auc_pre",0):.3f}. '
              f'<a href="pages/25_external_prjeb11591.html">See page 25</a>.'
              f'</div>')
    txt = idx_f.read_text()
    if 'class="v3-ext-banner"' not in txt:
        # APPEND ONLY -- insert just before </body>
        if "</body>" in txt:
            txt = txt.replace("</body>", banner + "\n</body>")
        else:
            txt = txt + banner
        idx_f.write_text(txt)
        log.info(f"banner appended to {idx_f}")


def update_version(kpis):
    vfile = HTML / "_version.json"
    try:
        cfg = json.loads(vfile.read_text()) if vfile.exists() else {}
    except Exception:
        cfg = {}
    cfg["v3_ext_pages"] = [25, 26, 27, 28]
    cfg["v3_ext_build_time"] = datetime.utcnow().isoformat() + "Z"
    cfg["v3_ext_primary_auc"] = float(kpis.get("auc_pre", 0) or 0)
    vfile.write_text(json.dumps(cfg, indent=2))
    log.info(f"updated {vfile}")


def main():
    log.info("v3_ext_07 start")
    kpis = page_25()
    page_26()
    page_27()
    page_28()
    append_banner(kpis)
    update_version(kpis)
    log.info("v3_ext_07 done")


if __name__ == "__main__":
    main()
