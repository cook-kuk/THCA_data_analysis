#!/usr/bin/env python3
"""Render a self-contained HTML page summarising the GPL570 validation pack.

Reads:
  project/results/p_gpl570_validation/gpl570_score_tests.tsv
  project/results/p_gpl570_validation/gpl570_meta_effect_summary.tsv
  project/results/p_gpl570_validation/sample_metadata.tsv
  project/results/p_gpl570_validation/gpl570_rai_lineage_boxplots.png
  project/results/p_gpl570_validation/gpl570_dm1_nonoverlap_scatter_grid.png
  project/results/p_gpl570_validation/gpl570_direction_consistency_forest.png

Writes:
  project/papers_hub_2026_05_04/gpl570_validation_pack.html

CPU-only; no external network; no manuscript prose.
"""
import base64
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RES  = ROOT / "project/results/p_gpl570_validation"
OUT_HTML = ROOT / "project/papers_hub_2026_05_04/gpl570_validation_pack.html"

def b64(p): return base64.b64encode(p.read_bytes()).decode("ascii")

def fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)): return "—"
    if p < 1e-3: return f"{p:.1e}"
    return f"{p:.3g}"

def fmt_d(d):
    if d is None or (isinstance(d, float) and np.isnan(d)): return "—"
    return f"{d:+.3f}"

# ---------- load ----------
tests = pd.read_csv(RES / "gpl570_score_tests.tsv", sep="\t")
meta  = pd.read_csv(RES / "sample_metadata.tsv", sep="\t")
img_box = b64(RES / "gpl570_rai_lineage_boxplots.png")
img_sca = b64(RES / "gpl570_dm1_nonoverlap_scatter_grid.png")
img_for = b64(RES / "gpl570_direction_consistency_forest.png")

# composition table
comp = (meta.groupby(["dataset","histology"]).size().unstack(fill_value=0)
            .reindex(columns=["Normal","PTC","ATC"], fill_value=0))
comp["Total"] = comp.sum(axis=1)
comp_html = comp.reset_index().to_html(index=False, classes="tbl", border=0)

# headline table
hd = tests[tests["test"].fillna("").str.contains(" :: ATC vs ")].copy()
hd["axis"] = hd["test"].apply(lambda s: s.split(" :: ")[0])
hd["contrast"] = hd["test"].apply(lambda s: s.split(" :: ")[1])
expected = {"RAI_8_score":-1, "DM1_like_score":+1, "THYROID_NONOVERLAP_score":-1,
            "TDS_like_score":-1, "TF_collapse_score":-1,
            "STAT3_AP1_DNMT_score":+1, "TACSTD2_z":+1}
hd["expected_sign"] = hd["axis"].map(expected)
hd["observed_sign"] = np.sign(hd["cohens_d"])
hd["ok"] = (hd["observed_sign"] == hd["expected_sign"]) & hd["cohens_d"].notna()

# Spearman table
sp = tests[tests["test"].fillna("").str.startswith("Spearman")].copy()

# direction consistency totals
n_lineage = hd[hd["axis"] != "TACSTD2_z"].shape[0]
n_lineage_ok = int(hd[hd["axis"] != "TACSTD2_z"]["ok"].sum())
n_trop = hd[hd["axis"] == "TACSTD2_z"].shape[0]
n_trop_ok = int(hd[hd["axis"] == "TACSTD2_z"]["ok"].sum())

# build rows
def headline_rows():
    out = []
    order = ["RAI_8_score","DM1_like_score","THYROID_NONOVERLAP_score","TDS_like_score",
             "TF_collapse_score","STAT3_AP1_DNMT_score","TACSTD2_z"]
    for ds in ["GSE33630","GSE29265","GSE65144"]:
        for axis in order:
            for contrast in ["ATC vs Normal","ATC vs PTC","ATC vs PDTC"]:
                row = hd[(hd["dataset"]==ds) & (hd["axis"]==axis) & (hd["contrast"]==contrast)]
                if row.empty: continue
                r = row.iloc[0]
                cls = "ok" if r["ok"] else "wrong"
                out.append(
                    f"<tr class='{cls}'><td>{ds}</td><td>{axis}</td><td>{contrast}</td>"
                    f"<td>{int(r['n_x'])}</td><td>{int(r['n_y'])}</td>"
                    f"<td>{fmt_d(r['cohens_d'])}</td><td>{fmt_p(r['MW_p_two_sided'])}</td>"
                    f"<td>{int(r['expected_sign']):+d}</td>"
                    f"<td>{'OK' if r['ok'] else 'WRONG'}</td></tr>")
    return "\n".join(out)

def spearman_rows():
    out = []
    for _, r in sp.iterrows():
        out.append(f"<tr><td>{r['dataset']}</td><td>{r['test'].replace('Spearman :: ','')}</td>"
                   f"<td>{r['rho']:.3f}</td><td>{fmt_p(r['p'])}</td><td>{int(r['n'])}</td></tr>")
    # add anchor row reference
    out.append("<tr class='anchor'><td>GSE76039 (anchor)</td><td>DM1_like vs THYROID_NONOVERLAP</td>"
               "<td>-0.925</td><td>3.1e-16</td><td>37</td></tr>")
    out.append("<tr class='anchor'><td>GSE76039 (anchor)</td><td>TF_collapse vs DM1_like</td>"
               "<td>-0.931</td><td>7.3e-17</td><td>37</td></tr>")
    return "\n".join(out)

HTML = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>GPL570 thyroid dedifferentiation validation pack — Paper 1 external</title>
<style>
:root {{
  --bg:#0b0d12; --panel:#141821; --panel2:#1a1f2c; --ink:#e9edf6; --mute:#9aa3b8;
  --accent:#7cf2c0; --accent2:#7cb6ff; --warn:#ffb86b; --bad:#ff6b6b; --ok:#7cf2c0;
  --line:#262c3a;
}}
* {{ box-sizing: border-box; }}
html, body {{ background: var(--bg); color: var(--ink); font: 14.5px/1.55 ui-sans-serif, system-ui, -apple-system, "Inter", "Helvetica Neue", Arial, sans-serif; margin:0; }}
.wrap {{ max-width: 1180px; margin: 0 auto; padding: 28px 22px 80px; }}
h1 {{ font-size: 30px; line-height:1.15; margin: 0 0 6px; letter-spacing:-0.01em; }}
h2 {{ font-size: 21px; margin: 36px 0 10px; padding-bottom: 8px; border-bottom: 1px solid var(--line); }}
h3 {{ font-size: 16px; margin: 20px 0 6px; color: var(--accent2); }}
p, li {{ color: var(--ink); }}
.lede {{ color: var(--mute); margin-bottom: 4px; }}
.dateline {{ color: var(--mute); font-size: 12.5px; margin-top: 0; }}
.kpis {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px; margin: 18px 0 22px; }}
.kpi {{ background: var(--panel); border: 1px solid var(--line); border-radius: 12px; padding: 14px 14px 12px; }}
.kpi .v {{ font-size: 22px; font-weight: 700; color: var(--accent); }}
.kpi .k {{ font-size: 12px; color: var(--mute); text-transform: uppercase; letter-spacing: .06em; }}
.kpi.warn .v {{ color: var(--warn); }} .kpi.bad .v {{ color: var(--bad); }}
.card {{ background: var(--panel); border: 1px solid var(--line); border-radius: 14px; padding: 16px 18px; margin: 14px 0; }}
.card.tight {{ padding: 12px 14px; }}
.card.warn {{ border-color: #5a4516; background: #1c1709; }}
.fig {{ background: var(--panel2); border: 1px solid var(--line); border-radius: 12px; padding: 14px; margin: 16px 0; }}
.fig img {{ width: 100%; height: auto; border-radius: 8px; background: #fff; }}
.fig .cap {{ font-size: 13px; color: var(--mute); margin-top: 10px; line-height:1.6; }}
.fig .cap b {{ color: var(--ink); }}
table.tbl {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
table.tbl th, table.tbl td {{ padding: 6px 8px; border-bottom: 1px solid var(--line); text-align: left; }}
table.tbl th {{ color: var(--mute); font-weight: 600; text-transform: uppercase; letter-spacing: .04em; font-size: 11.5px; }}
tr.ok td {{ color: var(--ok); }} tr.wrong td {{ color: var(--bad); }}
tr.anchor td {{ color: var(--accent2); font-style: italic; }}
.codey {{ font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace; font-size: 12.5px; background: #0f131c; padding: 2px 6px; border-radius: 6px; border: 1px solid var(--line); color: var(--accent); }}
.banner {{ background: linear-gradient(135deg, #0f3b2c, #06231a); border: 1px solid #235a47; border-radius: 14px; padding: 14px 18px; margin: 14px 0 8px; }}
.banner b {{ color: var(--ok); }}
.tag {{ display: inline-block; font-size: 11px; padding: 2px 8px; border-radius: 999px; border: 1px solid var(--line); color: var(--mute); margin-right: 6px; }}
.tag.go {{ color: var(--ok); border-color: #234737; }}
.tag.cau {{ color: var(--warn); border-color: #5a4516; }}
.tag.bad {{ color: var(--bad); border-color: #5a2424; }}
.grid2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 14px; }}
.grid3 {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 14px; }}
ul.fl {{ margin: 6px 0 0 0; padding-left: 18px; }} ul.fl li {{ margin: 2px 0; }}
.small {{ font-size: 12.5px; color: var(--mute); }}
.foot {{ color: var(--mute); font-size: 12px; margin-top: 30px; padding-top: 14px; border-top: 1px solid var(--line); }}
hr {{ border:0; border-top:1px solid var(--line); margin: 22px 0; }}
@media (max-width: 900px) {{ .kpis {{ grid-template-columns: repeat(2,1fr); }} .grid2,.grid3 {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>
<div class="wrap">

<h1>GPL570 thyroid dedifferentiation validation pack</h1>
<p class="lede">External replication of the driver-orthogonal transcriptional differentiation axis (RAI lineage silencing readout) across three independent Affymetrix HG-U133 Plus 2.0 thyroid cancer cohorts, anchored to the GSE76039 / Landa 2016 ATC+PDTC first-pass.</p>
<p class="dateline">Paper 1 — external GPL570 first pass · 2026-05-04 · CPU-only · within-cohort z-score per dataset · no TCGA classifier transfer · no batch pooling · no survival/mutation modelling</p>

<div class="banner">
  <b>Headline.</b> 30/30 lineage-axis × dataset × ATC-contrast cells direction-consistent with the GSE76039 anchor across <b>GSE33630</b>, <b>GSE29265</b>, <b>GSE65144</b>.
  Within-cohort orthogonality between the 8-gene RAI readout and a zero-overlap THYROID_NONOVERLAP module replicates: Spearman ρ in <b>−0.85 to −0.94</b> (anchor −0.93).
  <b>TROP2 (TACSTD2) caveat:</b> only 1/5 ATC contrasts go in the expected direction at GPL570 bulk level — supplement-only / spatial follow-up.
</div>

<div class="kpis">
  <div class="kpi"><div class="v">216</div><div class="k">GPL570 thyroid samples (incl. anchor)</div></div>
  <div class="kpi"><div class="v">{n_lineage_ok}/{n_lineage}</div><div class="k">Lineage axis cells consistent</div></div>
  <div class="kpi"><div class="v">22 / 22</div><div class="k">Panel genes mapped, all 3 datasets</div></div>
  <div class="kpi warn"><div class="v">{n_trop_ok}/{n_trop}</div><div class="k">TACSTD2 cells consistent</div></div>
</div>

<h2>1. Cohort composition</h2>
<div class="card">
  <p class="small">Histology label sources: <span class="codey">GSE33630</span> uses <span class="codey">!Sample_characteristics_ch1: pathological diagnostic</span>;
  <span class="codey">GSE29265</span> uses <span class="codey">!Sample_title</span> (ATC / PTC / Patient-matched non-tumor); <span class="codey">GSE65144</span> uses <span class="codey">tissue type</span>.
  Labels taken at face value from GEO; no central re-review. PDTC absent from the validation pack — GSE76039 remains the only PDTC source.</p>
  {comp_html}
</div>

<h2>2. Boxplot grid — module scores by dataset</h2>
<div class="fig">
  <img alt="GPL570 boxplot grid: 7 axes × 3 datasets, Normal/PTC/ATC" src="data:image/png;base64,{img_box}" />
  <div class="cap">
    <b>Figure 1.</b> Within-cohort z-score distributions for the seven module/single-gene axes (rows) across the three GPL570 cohorts (columns). Each row, top to bottom: <b>RAI_8</b> (8-gene RAI lineage), <b>DM1_like</b> (= −RAI_8), <b>THYROID_NONOVERLAP</b> (zero-overlap 8-gene lineage module), <b>TDS_like</b> (16-gene union), <b>TF collapse</b> (FOXE1 / NKX2-1 / PAX8 / HHEX), <b>STAT3 / AP-1 / DNMT axis</b>, and <b>TACSTD2 (TROP2)</b> z. Box colours: Normal (green), PTC (blue), ATC (red). Sample sizes shown under each box. Score is the per-sample mean z across panel genes; per-gene z is computed within the dataset (no cross-cohort, no cross-platform normalisation). The lineage-silencing pattern in ATC is visible as a downward shift of all five lineage axes in every cohort where ATC is sampled, with parallel upward shift of the inflammation/methylation axis. Note that the TACSTD2 row does <b>not</b> show a consistent ATC-elevated pattern across cohorts at GPL570 bulk level — see Section 4.
  </div>
</div>

<h2>3. Within-cohort orthogonality — DM1_like vs THYROID_NONOVERLAP</h2>
<div class="fig">
  <img alt="DM1_like vs THYROID_NONOVERLAP scatter grid" src="data:image/png;base64,{img_sca}" />
  <div class="cap">
    <b>Figure 2.</b> Per-sample DM1_like score (x) versus THYROID_NONOVERLAP score (y) within each cohort, with Spearman ρ shown in the title. Points coloured by histology (Normal green, PTC blue, ATC red). DM1_like and THYROID_NONOVERLAP share <b>no genes</b> by construction, so a strong negative ρ here means independent zero-overlap evidence that the same biological axis is operative in each cohort. Observed ρ values: GSE33630 −0.941 (n=105, p=2.8e−50), GSE29265 −0.846 (n=49, p=2.1e−14), GSE65144 −0.904 (n=25, p=6.0e−10). The GSE76039 anchor reported ρ=−0.925 (n=37, p=3.1e−16). The orthogonality reproduces in every independent cohort tested.
  </div>
</div>

<h2>4. Direction-consistency forest</h2>
<div class="fig">
  <img alt="Cohen's d forest, color-coded by direction-consistency vs expected sign" src="data:image/png;base64,{img_for}" />
  <div class="cap">
    <b>Figure 3.</b> Cohen's d for headline ATC contrasts (ATC vs PTC and ATC vs Normal where available), one bar per axis × dataset × contrast. Bar colour: <span style="color:#2ca02c">green</span> = observed sign matches the pre-specified expected sign (lineage axes negative; STAT3_AP1_DNMT and TACSTD2 positive); <span style="color:#d62728">red</span> = wrong direction. <b>30/30 lineage and inflammation/methylation cells direction-consistent.</b> The <b>TACSTD2 (TROP2)</b> rows show 4/5 in the wrong direction at GPL570 bulk level — only GSE65144 ATC vs Normal goes in the expected positive direction (d=+0.25, n.s.). The largest lineage effects come from GSE33630 ATC vs Normal (THYROID_NONOVERLAP d=−6.53; TDS_like d=−6.10; RAI_8 d=−5.43); GSE65144 and GSE29265 sit at d ≈ −2 to −3.
  </div>
</div>

<h2>5. Headline contrast table</h2>
<div class="card tight">
<table class="tbl">
<thead><tr><th>Dataset</th><th>Axis</th><th>Contrast</th><th>n_x</th><th>n_y</th><th>Cohen d</th><th>MW p</th><th>Exp. sign</th><th>Status</th></tr></thead>
<tbody>
{headline_rows()}
</tbody>
</table>
<p class="small">Status = OK if observed sign of Cohen d matches the pre-specified expected sign (lineage axes: ATC < comparator → d &lt; 0; STAT3_AP1_DNMT and TACSTD2: ATC &gt; comparator → d &gt; 0). All p-values are two-sided Mann–Whitney U.</p>
</div>

<h2>6. Spearman summary — module-vs-module orthogonality</h2>
<div class="card tight">
<table class="tbl">
<thead><tr><th>Cohort</th><th>Pair</th><th>Spearman ρ</th><th>p</th><th>n</th></tr></thead>
<tbody>
{spearman_rows()}
</tbody>
</table>
<p class="small">All four GPL570 thyroid cohorts (3 validation + GSE76039 anchor) show strong negative ρ for DM1_like vs THYROID_NONOVERLAP — independent zero-overlap evidence of the same axis. TACSTD2 vs DM1_like is positive but moderate in GSE33630 / GSE29265, and ~0 in GSE65144, consistent with the bulk-microarray TROP2 caveat (Figure 3 / §4).</p>
</div>

<h2>7. Verdict per claim</h2>
<div class="grid2">
  <div class="card"><span class="tag go">MAIN-FIGURE EXTENSION</span>
  <h3>Driver-orthogonal differentiation axis replicates externally</h3>
  <p>RAI_8 / DM1_like and the zero-overlap THYROID_NONOVERLAP module both replicate ATC silencing across all 5 ATC contrasts in 3 independent cohorts (10/10 lineage cells consistent). Effect size scales with comparator: ATC vs Normal &gt; ATC vs PTC &gt; ATC vs PDTC.</p></div>

  <div class="card"><span class="tag go">SUPPLEMENT</span>
  <h3>Inflammation / AP-1 / DNMT axis up in ATC</h3>
  <p>STAT3_AP1_DNMT_score: 5/5 ATC contrasts direction-consistent, median d ≈ +1.18. Smaller magnitude than the lineage axes; supportive mechanism-side panel, not a primary claim.</p></div>

  <div class="card warn"><span class="tag cau">INTERNAL ONLY</span>
  <h3>TROP2 elevation in advanced disease</h3>
  <p>TACSTD2 z at GPL570 bulk level does <b>not</b> replicate the ATC-elevated pattern: only 1/5 ATC contrasts in expected direction; GSE33630 ATC vs PTC is significantly negative (d=−1.72, p=8.3e−4). Bulk-mRNA microarray cannot resolve the spatial / clonal TROP2 question — keep as supplement framing or hand off to a dedicated spatial analysis. Do not headline as "TROP2 elevated in advanced disease" from this pack.</p></div>

  <div class="card"><span class="tag go">REPRODUCED</span>
  <h3>Within-cohort module orthogonality</h3>
  <p>DM1_like vs THYROID_NONOVERLAP: ρ in [−0.846, −0.941] across the 3 validation cohorts; ρ=−0.925 in GSE76039 anchor. The orthogonality replicates without overlap of genes between modules.</p></div>
</div>

<h2>8. Risks &amp; framing discipline</h2>
<div class="grid2">
  <div class="card">
    <h3>Platform &amp; calibration</h3>
    <ul class="fl">
      <li>GPL570 is bulk Affymetrix microarray — within-cohort z only, no shared cutoff, no TCGA-trained absolute-form transfer (per <span class="codey">v17_korean_K2_calibration</span>).</li>
      <li>Best-mean-probe collapse per gene; identical methodology to the GSE76039 anchor; NKX2-1 mapped via TITF1 alias.</li>
      <li>22/22 panel genes mapped in all three datasets (no missing-gene caveat).</li>
    </ul>
  </div>
  <div class="card">
    <h3>Histology label quality</h3>
    <ul class="fl">
      <li>Labels taken at face value from GEO; no central pathology re-review.</li>
      <li>No PDTC in any of the three validation cohorts — GSE76039 remains the only PDTC anchor in the GPL570 family.</li>
      <li>GSE65144 lacks PTC; GSE29265 ATC arm is small (n=9); per-arm n must be respected for absolute-magnitude inference.</li>
    </ul>
  </div>
  <div class="card">
    <h3>What we do NOT claim</h3>
    <ul class="fl">
      <li>NOT clinical validation; NOT survival validation; NOT mutation validation; NOT fusion validation — none of these fields exists in the GEO records.</li>
      <li>NOT "monotonic 3-stage progression" — Normal → PTC → ATC effect-size ordering is direction-consistent but not asserted as ordinal.</li>
      <li>NOT "PTC → ATC progression proven" from this pack alone.</li>
      <li>NOT "TROP2 elevated in advanced disease" at GPL570 bulk level — see Figure 3 / §4 / §7.</li>
    </ul>
  </div>
  <div class="card">
    <h3>Safe wording</h3>
    <ul class="fl">
      <li>"External GPL570 expression validation across three independent thyroid cancer cohorts."</li>
      <li>"Direction-consistent lineage-silencing readout across independent cohorts."</li>
      <li>"Independent zero-overlap module replicates the same axis."</li>
      <li>"Within-cohort z-score per dataset; no cross-platform classifier transfer."</li>
    </ul>
  </div>
</div>

<h2>9. Files produced</h2>
<div class="card tight"><pre style="margin:0; white-space: pre-wrap; font-family: ui-monospace, monospace; font-size: 12px; color: var(--mute);">project/results/p_gpl570_validation/
├── sample_metadata.tsv                              (179 rows; dataset / sample_id / histology_raw / histology / disease_group / source_file)
├── probe_to_gene_panel.tsv                          (66 rows; 22 genes × 3 datasets, best-mean-probe collapse)
├── GSE33630_expression_gene_log.tsv.gz              (panel-gene log expression, 22 × 105)
├── GSE29265_expression_gene_log.tsv.gz              (panel-gene log expression, 22 ×  49)
├── GSE65144_expression_gene_log.tsv.gz              (panel-gene log expression, 22 ×  25)
├── GSE33630_scores.tsv                              (per-sample module scores, n=105)
├── GSE29265_scores.tsv                              (per-sample module scores, n= 49)
├── GSE65144_scores.tsv                              (per-sample module scores, n= 25)
├── gpl570_score_tests.tsv                           (all MW + Cohen d + Spearman tests, long format)
├── gpl570_meta_effect_summary.tsv                   (coverage + per-dataset axis summary + overall summary + headline)
├── gpl570_rai_lineage_boxplots.png                  (Figure 1 — 7×3 boxplot grid)
├── gpl570_dm1_nonoverlap_scatter_grid.png           (Figure 2 — orthogonality scatter)
├── gpl570_direction_consistency_forest.png          (Figure 3 — Cohen d forest)
└── raw/  (gitignored: 3 Series Matrix .txt.gz, ~50 MB)
</pre>
</div>

<p class="foot">Generated locally from project/scripts/build_gpl570_validation_html.py · GSE76039 anchor values from project/results/p_landa_2016/score_tests.tsv · No manuscript prose drafted in this pass.</p>
</div>
</body>
</html>
"""

OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
OUT_HTML.write_text(HTML, encoding="utf-8")
print(f"[save] {OUT_HTML}  size={OUT_HTML.stat().st_size//1024} KB")
