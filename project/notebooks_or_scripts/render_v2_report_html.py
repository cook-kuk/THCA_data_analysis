#!/usr/bin/env python3
from __future__ import annotations

import json
from html import escape
from pathlib import Path

import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PROJECT = ROOT / "project"
REPORTS = PROJECT / "reports"
RESULTS = PROJECT / "results"
META = PROJECT / "metadata"
OUTDIR = REPORTS / "html"


REFERENCES = [
    {
        "group": "Paper",
        "title": "Integrated Genomic Characterization of Papillary Thyroid Carcinoma",
        "label": "TCGA Cell 2014 / PMID 25417114",
        "url": "https://pubmed.ncbi.nlm.nih.gov/25417114/",
        "why": "Primary source for the thyroid differentiation score genes and TCGA papillary thyroid carcinoma molecular framework.",
    },
    {
        "group": "Portal",
        "title": "GDC publication landing page for TCGA-THCA",
        "label": "GDC THCA 2014 publication page",
        "url": "https://gdc.cancer.gov/about-data/publications/thca_2014",
        "why": "Public entry point used for TCGA-THCA context and open GDC project linkage.",
    },
    {
        "group": "Paper",
        "title": "Unraveling the role of the mitochondrial one-carbon pathway in undifferentiated thyroid cancer by multi-omics analyses",
        "label": "Nature Communications 2024",
        "url": "https://www.nature.com/articles/s41467-024-45366-0",
        "why": "Source used to correct the GSE213647 cohort interpretation and to locate the supplementary RNA cohort sheet.",
    },
    {
        "group": "Supplementary",
        "title": "Nature Communications supplementary workbook for the SHMT2 cohort",
        "label": "Supplementary RNA cohort sheet workbook",
        "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-45366-0/MediaObjects/41467_2024_45366_MOESM5_ESM.xlsx",
        "why": "Direct source used to relabel GSE213647 sample tissue types.",
    },
    {
        "group": "Dataset",
        "title": "GEO accession GSE213647",
        "label": "GSE213647 GEO page",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213647",
        "why": "Public accession page for the Korean thyroid RNA-seq cohort.",
    },
    {
        "group": "Dataset",
        "title": "GEO accession GSE126698",
        "label": "GSE126698 GEO page",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126698",
        "why": "External RNA-seq cohort used for exploratory restricted validation.",
    },
    {
        "group": "Dataset",
        "title": "GEO accession GSE27155",
        "label": "GSE27155 GEO page",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27155",
        "why": "Microarray cohort used as histology-proxy external validation.",
    },
    {
        "group": "Dataset",
        "title": "GEO accession GSE76039",
        "label": "GSE76039 GEO page",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76039",
        "why": "Dedifferentiated microarray cohort retained as a reference dataset in the rerun outputs.",
    },
]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def as_records(df: pd.DataFrame, float_digits: int = 4) -> list[dict]:
    out = []
    for _, row in df.iterrows():
        rec = {}
        for col, value in row.items():
            if pd.isna(value):
                rec[col] = None
            elif isinstance(value, float):
                rec[col] = round(value, float_digits)
            else:
                rec[col] = value
        out.append(rec)
    return out


def build_dataset_summary(sample_master: pd.DataFrame, coverage: pd.DataFrame) -> list[dict]:
    rows = []
    for dataset in sorted(sample_master["dataset"].dropna().unique()):
        sub = sample_master[sample_master["dataset"] == dataset].copy()
        cov = coverage[coverage["dataset"] == dataset]
        cov_row = cov.iloc[0].to_dict() if not cov.empty else {}
        rows.append(
            {
                "dataset": dataset,
                "samples": int(sub.shape[0]),
                "normal": int(sub["normal_vs_tumor"].astype(str).eq("normal").sum()),
                "tumor": int(sub["normal_vs_tumor"].astype(str).eq("tumor").sum()),
                "unknown_or_other": int(sub.shape[0] - sub["normal_vs_tumor"].astype(str).isin(["normal", "tumor"]).sum()),
                "histology_values": ", ".join(sorted(x for x in sub["histology_subtype"].dropna().astype(str).unique() if x and x != "nan")),
                "label_confidence_values": ", ".join(sorted(x for x in sub["label_confidence"].dropna().astype(str).unique() if x and x != "nan")),
                "tds16_cov": cov_row.get("TDS16_coverage_frac"),
                "brs71_cov": cov_row.get("BRS71_proxy_coverage_frac"),
                "tiera67_cov": cov_row.get("TierA67_entries_coverage_frac"),
            }
        )
    return rows


def build_roc_options(internal: pd.DataFrame, external: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    internal_options = []
    for _, row in internal.sort_values(["feature_set", "model"]).iterrows():
        path = (
            f"../../results/figs/roc_braf_like_vs_ras_like_TCGA-THCA_{row['feature_set']}_{row['model']}.html"
        )
        if not (PROJECT / "results" / "figs" / f"roc_braf_like_vs_ras_like_TCGA-THCA_{row['feature_set']}_{row['model']}.html").exists():
            continue
        internal_options.append(
            {
                "label": f"TCGA internal | {row['feature_set']} | {row['model']}",
                "feature_set": row["feature_set"],
                "model": row["model"],
                "path": path,
                "auc": row["cv_auc"],
                "bacc": row["cv_balanced_accuracy"],
                "f1": row["cv_f1"],
                "mcc": row["cv_mcc"],
            }
        )
    external_options = []
    for dataset in sorted(external["dataset"].dropna().unique()):
        path = f"../../results/figs/roc_braf_like_vs_ras_like_{dataset}.html"
        if not (PROJECT / "results" / "figs" / f"roc_braf_like_vs_ras_like_{dataset}.html").exists():
            continue
        best = external[(external["dataset"] == dataset) & external["auc"].notna()].sort_values(
            ["auc", "balanced_accuracy"], ascending=[False, False]
        )
        if best.empty:
            continue
        row = best.iloc[0]
        external_options.append(
            {
                "dataset": dataset,
                "label": f"{dataset} external | {row['feature_set']} | {row['model']}",
                "feature_set": row["feature_set"],
                "model": row["model"],
                "path": path,
                "auc": row["auc"],
                "bacc": row["balanced_accuracy"],
                "f1": row["f1"],
                "mcc": row["mcc"],
                "n_samples": int(row["n_samples"]),
            }
        )
    return internal_options, external_options


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    internal = read_tsv(RESULTS / "ml" / "baseline_ml_results.tsv")
    external = read_tsv(RESULTS / "ml" / "baseline_ml_external.tsv")
    coverage = read_tsv(RESULTS / "tables" / "gene_coverage_table.tsv")
    sample_master = read_tsv(META / "sample_master.tsv")

    dataset_summary = build_dataset_summary(sample_master, coverage)
    internal_options, external_options = build_roc_options(internal, external)

    best_internal = internal.sort_values(["cv_auc", "cv_balanced_accuracy"], ascending=[False, False]).iloc[0].to_dict()
    best_external = (
        external[external["auc"].notna()].sort_values(["auc", "balanced_accuracy"], ascending=[False, False]).iloc[0].to_dict()
    )
    best_external_nontrivial = (
        external[(external["dataset"] == "GSE27155") & external["auc"].notna()]
        .sort_values(["auc", "balanced_accuracy"], ascending=[False, False])
        .iloc[0]
        .to_dict()
    )

    payload = {
        "best_internal": best_internal,
        "best_external": best_external,
        "best_external_nontrivial": best_external_nontrivial,
        "internal_rows": as_records(internal.sort_values(["cv_auc", "cv_balanced_accuracy"], ascending=[False, False])),
        "external_rows": as_records(external.sort_values(["auc", "balanced_accuracy"], ascending=[False, False], na_position="last")),
        "coverage_rows": as_records(coverage),
        "dataset_rows": dataset_summary,
        "internal_roc_options": internal_options,
        "external_roc_options": external_options,
        "references": REFERENCES,
    }
    payload_json = json.dumps(payload).replace("</", "<\\/")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>THCA v2 Interactive Report</title>
  <style>
    :root {{
      --bg: #f4efe7;
      --paper: #fffdf9;
      --ink: #201a17;
      --muted: #6f6257;
      --line: #dfd1c2;
      --accent: #9e2f2b;
      --accent2: #0d6765;
      --gold: #a36a00;
      --chip: #efe4d6;
      --shadow: 0 16px 40px rgba(49, 34, 22, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      background:
        radial-gradient(circle at 0% 0%, rgba(158,47,43,0.10), transparent 24%),
        radial-gradient(circle at 100% 0%, rgba(13,103,101,0.12), transparent 22%),
        linear-gradient(180deg, #f7f3ed 0%, var(--bg) 100%);
      font-family: "IBM Plex Sans", "Segoe UI", sans-serif;
    }}
    .wrap {{ max-width: 1500px; margin: 0 auto; padding: 24px 22px 40px; }}
    .hero {{
      background: linear-gradient(135deg, #201714 0%, #4a261d 45%, #7b3a24 100%);
      color: #fff8f0;
      border-radius: 24px;
      padding: 26px 28px;
      box-shadow: 0 24px 64px rgba(35, 18, 12, 0.16);
    }}
    .hero h1 {{ margin: 0 0 8px; font-size: 38px; line-height: 1.02; }}
    .hero p {{ margin: 0; max-width: 1050px; color: rgba(255,248,240,0.85); font-size: 15px; }}
    .hero-meta {{
      display: flex; flex-wrap: wrap; gap: 10px; margin-top: 16px;
    }}
    .hero-meta span {{
      display: inline-block; padding: 8px 12px; background: rgba(255,255,255,0.12);
      border: 1px solid rgba(255,255,255,0.14); border-radius: 999px; font-size: 13px;
    }}
    .tabs {{
      display: flex; flex-wrap: wrap; gap: 10px; margin: 18px 0 18px;
    }}
    .tab-btn {{
      border: 1px solid var(--line);
      background: rgba(255,255,255,0.7);
      color: var(--ink);
      border-radius: 999px;
      padding: 10px 14px;
      cursor: pointer;
      font-weight: 600;
    }}
    .tab-btn.active {{ background: var(--accent); color: white; border-color: var(--accent); }}
    .tab-panel {{ display: none; }}
    .tab-panel.active {{ display: block; }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 14px;
      margin-top: 18px;
    }}
    .card, .panel {{
      background: var(--paper);
      border: 1px solid var(--line);
      border-radius: 18px;
      box-shadow: var(--shadow);
    }}
    .card {{ padding: 18px; }}
    .panel {{ padding: 18px; margin-top: 16px; }}
    .kicker {{ color: var(--muted); text-transform: uppercase; font-size: 12px; letter-spacing: 0.10em; }}
    .metric {{ font-size: 30px; font-weight: 700; margin-top: 6px; }}
    .sub {{ color: var(--muted); font-size: 14px; margin-top: 6px; line-height: 1.45; }}
    .two {{ display: grid; grid-template-columns: 1.15fr 0.85fr; gap: 16px; }}
    .three {{ display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 16px; }}
    h2 {{ margin: 0 0 10px; font-size: 24px; }}
    h3 {{ margin: 0 0 8px; font-size: 18px; }}
    p {{ line-height: 1.55; }}
    ul {{ margin: 8px 0 0; padding-left: 18px; }}
    li {{ margin: 6px 0; }}
    .callout {{
      border-left: 4px solid var(--gold);
      background: #fff6e6;
      color: #624800;
      border-radius: 10px;
      padding: 12px 14px;
      margin-top: 12px;
      line-height: 1.5;
    }}
    .chip {{
      display: inline-block;
      padding: 6px 10px;
      border-radius: 999px;
      background: var(--chip);
      margin: 4px 6px 0 0;
      font-size: 13px;
    }}
    .controls {{
      display: flex;
      gap: 12px;
      flex-wrap: wrap;
      align-items: center;
      margin-bottom: 12px;
    }}
    label {{ font-size: 13px; color: var(--muted); font-weight: 600; display: block; margin-bottom: 4px; }}
    select, input {{
      min-width: 220px;
      border: 1px solid var(--line);
      border-radius: 10px;
      padding: 10px 12px;
      background: white;
      font: inherit;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      display: block;
      overflow-x: auto;
      white-space: nowrap;
      font-size: 13px;
    }}
    th {{
      position: sticky;
      top: 0;
      background: #f3eadf;
      border-bottom: 1px solid var(--line);
      text-align: left;
      padding: 10px;
    }}
    td {{
      border-bottom: 1px solid #efe4d7;
      padding: 9px 10px;
      vertical-align: top;
    }}
    tr:nth-child(even) td {{ background: #fffaf4; }}
    .mini {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px;
      margin-top: 12px;
    }}
    .mini .card {{ padding: 14px; }}
    iframe {{
      width: 100%;
      min-height: 520px;
      border: 1px solid var(--line);
      border-radius: 16px;
      background: white;
    }}
    .links a, .ref-url {{
      color: var(--accent2);
      text-decoration: none;
      font-weight: 600;
    }}
    .ref-card {{
      padding: 16px;
      border: 1px solid var(--line);
      border-radius: 16px;
      background: var(--paper);
    }}
    .ref-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 14px;
    }}
    .foot {{ color: var(--muted); font-size: 13px; margin-top: 16px; }}
    @media (max-width: 1100px) {{
      .two, .three {{ grid-template-columns: 1fr; }}
      .mini {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
    }}
    @media (max-width: 700px) {{
      .mini {{ grid-template-columns: 1fr; }}
      .hero h1 {{ font-size: 30px; }}
    }}
  </style>
</head>
<body>
  <div class="wrap">
    <section class="hero">
      <div class="kicker">Reproducible Thyroid Cancer Dashboard</div>
      <h1>THCA rerun v2 interactive report</h1>
      <p>This page is generated from saved rerun outputs in <code>project/results</code>, <code>project/metadata</code>, and <code>project/reports</code>. It separates what is directly supported by the current data from what remains proxy-based or uncertain. The strongest-looking external result is intentionally not treated as definitive because of cohort size and label quality limits.</p>
      <div class="hero-meta">
        <span>TDS16 fixed from Cell 2014</span>
        <span>BRS71 explicitly marked as proxy</span>
        <span>TierA67 hardcoded composite</span>
        <span>GSE213647 relabel corrected from supplementary sheet</span>
        <span>Server path: /reports/html/index.html</span>
      </div>
    </section>

    <div class="tabs">
      <button class="tab-btn active" data-tab="overview">Overview</button>
      <button class="tab-btn" data-tab="methods">Methods</button>
      <button class="tab-btn" data-tab="models">Model Explorer</button>
      <button class="tab-btn" data-tab="datasets">Datasets</button>
      <button class="tab-btn" data-tab="coverage">Coverage</button>
      <button class="tab-btn" data-tab="references">References</button>
    </div>

    <section class="tab-panel active" id="tab-overview">
      <div class="grid">
        <div class="card"><div class="kicker">TDS16</div><div class="metric">16</div><div class="sub">Confirmed from Cell 2014 Table S5.</div></div>
        <div class="card"><div class="kicker">BRS71</div><div class="metric">71</div><div class="sub">Derived internally from TCGA-THCA as a proxy.</div></div>
        <div class="card"><div class="kicker">TierA67</div><div class="metric">67</div><div class="sub">Requested hardcoded entries; 66 unique symbols after HGNC normalization.</div></div>
        <div class="card"><div class="kicker">GSE213647</div><div class="metric">fixed</div><div class="sub">Relabeled using the Nature Communications supplementary RNA cohort sheet.</div></div>
        <div class="card"><div class="kicker">Best Internal CV</div><div class="metric">{best_internal["cv_auc"]:.3f}</div><div class="sub">{escape(str(best_internal["feature_set"]))} + {escape(str(best_internal["model"]))}</div></div>
        <div class="card"><div class="kicker">Best External AUC</div><div class="metric">{best_external["auc"]:.3f}</div><div class="sub">{escape(str(best_external["dataset"]))} + {escape(str(best_external["feature_set"]))} + {escape(str(best_external["model"]))}</div></div>
      </div>

      <div class="two">
        <div class="panel">
          <h2>Headline Readout</h2>
          <ul>
            <li>The rerun fixed the gene panel definitions and replaced the earlier unresolved placeholders.</li>
            <li>TCGA mutation-anchor training labels now separate <code>BRAF_like</code> and <code>RAS_like</code> using open GDC masked MAFs.</li>
            <li>GSE213647 was not kept as a 2-class external validation cohort after relabeling because the corrected cohort lacks a comparable RAS-like differentiated tumor class.</li>
            <li>The best-looking external numeric result is GSE126698 AUC {best_external["auc"]:.3f}, but this uses only {int(best_external["n_samples"])} samples and should be treated as exploratory.</li>
            <li>The more realistic stress-test dataset remains GSE27155, where the best AUC is {best_external_nontrivial["auc"]:.3f} but threshold-dependent metrics still degrade.</li>
          </ul>
          <div class="callout">
            Read AUC together with balanced accuracy and MCC. Several GSE27155 rows show high AUC but poor thresholded performance, which is what you would expect from a label-inferred histology proxy set rather than a clean molecularly labeled validation cohort.
          </div>
        </div>
        <div class="panel">
          <h2>What Is Still Uncertain</h2>
          <ul>
            <li>BRS71 is a proxy list, not the original published supplementary list.</li>
            <li>GSE27155 external labels are phenotype proxies inferred from histology and morphology.</li>
            <li>GSE126698 is too small for strong external claims.</li>
            <li>GSE213647 processed expression currently shows zero TDS16/BRS71/TierA67 coverage in the existing matrix output, which is a processing red flag and should not be overinterpreted biologically.</li>
            <li>TierA67 contains 67 requested entries but one duplicated symbol by design: <code>PAX8</code>.</li>
          </ul>
        </div>
      </div>

      <div class="panel">
        <h2>Direct File Links</h2>
        <div class="links">
          <a href="../analysis_summary_v2.md">analysis_summary_v2.md</a>
          <a href="../ml_report_v2.md">ml_report_v2.md</a>
          <a href="../next_steps_v2.md">next_steps_v2.md</a>
          <a href="../../results/ml/baseline_ml_results.tsv">baseline_ml_results.tsv</a>
          <a href="../../results/ml/baseline_ml_external.tsv">baseline_ml_external.tsv</a>
          <a href="../../results/tables/gene_coverage_table.tsv">gene_coverage_table.tsv</a>
          <a href="../../metadata/tds16_genes.txt">tds16_genes.txt</a>
          <a href="../../metadata/brs71_genes.txt">brs71_genes.txt</a>
          <a href="../../metadata/tierA67_genes.txt">tierA67_genes.txt</a>
          <a href="../../logs/label_correction_GSE213647.log">label_correction_GSE213647.log</a>
        </div>
      </div>
    </section>

    <section class="tab-panel" id="tab-methods">
      <div class="three">
        <div class="panel">
          <h2>Training Data</h2>
          <ul>
            <li>Backbone cohort: <code>TCGA-THCA</code> bulk RNA-seq processed matrix already saved in the project.</li>
            <li>Anchor source: locally available open GDC masked somatic MAF files.</li>
            <li><code>BRAF_like</code>: BRAF hotspot anchor only.</li>
            <li><code>RAS_like</code>: KRAS, HRAS, or NRAS mutation with no BRAF anchor.</li>
            <li>Samples with both anchors or neither anchor were excluded from training.</li>
          </ul>
        </div>
        <div class="panel">
          <h2>Feature Sets</h2>
          <ul>
            <li><code>TDS16</code>: official 16-gene list from Cell 2014.</li>
            <li><code>TierA67</code>: requested hardcoded composite panel. Implemented as 67 entries and 66 unique HGNC symbols.</li>
            <li><code>variance_top50</code>: selected inside each training fold only, not globally before CV.</li>
            <li><code>BRS71</code> was derived and saved as a proxy, but not used as a direct feature set because the rerun task specified TDS16, TierA67, and variance-top50.</li>
          </ul>
        </div>
        <div class="panel">
          <h2>Leakage Control</h2>
          <ul>
            <li>5-fold stratified CV on TCGA only.</li>
            <li>Median imputation, scaling, and model fitting all happen inside the sklearn pipeline.</li>
            <li>Variance-based feature selection is performed inside the training fold.</li>
            <li>No pooled TCGA + GEO training was used.</li>
          </ul>
        </div>
      </div>

      <div class="two">
        <div class="panel">
          <h2>External Validation Logic</h2>
          <ul>
            <li><code>GSE27155</code>: histology and morphology were binarized into BRAF-like and RAS-like proxies. MTC was excluded.</li>
            <li><code>GSE126698</code>: restricted exploratory subset mapping cPTC to BRAF-like proxy and FTC to RAS-like proxy.</li>
            <li><code>GSE213647</code>: relabeling fixed the GEO misinterpretation, but the corrected cohort does not support a comparable BRAF-like vs RAS-like 2-class differentiated tumor test, so it was excluded.</li>
          </ul>
          <div class="callout">
            A high external AUC is not sufficient here. When label construction is proxy-based or sample size is very small, the thresholded metrics and cohort context matter more than a single headline number.
          </div>
        </div>
        <div class="panel">
          <h2>ROC Generation</h2>
          <ul>
            <li>Internal ROC curves were generated from cross-validated out-of-fold scores on TCGA.</li>
            <li>External ROC curves were generated from model scores on the chosen external cohort rows.</li>
            <li>Implementation used sklearn <code>roc_curve</code> and <code>roc_auc_score</code>, then saved HTML wrappers embedding PNG plots.</li>
            <li>Available models in this rerun: Logistic Regression L2, Logistic Regression Elastic Net, Random Forest, Gradient Boosting. <code>xgboost</code> was skipped because the package was not installed.</li>
          </ul>
          <div class="chip">TCGA anchors used in training: 279 BRAF_like, 54 RAS_like</div>
        </div>
      </div>
    </section>

    <section class="tab-panel" id="tab-models">
      <div class="panel">
        <h2>Interactive ROC Viewer</h2>
        <div class="controls">
          <div>
            <label for="rocMode">ROC set</label>
            <select id="rocMode">
              <option value="internal">Internal CV (TCGA)</option>
              <option value="external">External validation</option>
            </select>
          </div>
          <div>
            <label for="rocSelect">Curve</label>
            <select id="rocSelect"></select>
          </div>
        </div>
        <div class="mini">
          <div class="card"><div class="kicker">AUC</div><div class="metric" id="rocAuc">-</div></div>
          <div class="card"><div class="kicker">Balanced Accuracy</div><div class="metric" id="rocBacc">-</div></div>
          <div class="card"><div class="kicker">F1</div><div class="metric" id="rocF1">-</div></div>
          <div class="card"><div class="kicker">MCC</div><div class="metric" id="rocMcc">-</div></div>
        </div>
        <div class="sub" id="rocMeta" style="margin:12px 0 10px;"></div>
        <iframe id="rocFrame" src=""></iframe>
      </div>

      <div class="two">
        <div class="panel">
          <h2>Internal Result Table</h2>
          <div class="controls">
            <div>
              <label for="internalFilter">Filter internal rows</label>
              <input id="internalFilter" type="text" placeholder="feature set or model">
            </div>
          </div>
          <div id="internalTable"></div>
        </div>
        <div class="panel">
          <h2>External Result Table</h2>
          <div class="controls">
            <div>
              <label for="externalFilter">Filter external rows</label>
              <input id="externalFilter" type="text" placeholder="dataset, model, feature set">
            </div>
          </div>
          <div id="externalTable"></div>
        </div>
      </div>
    </section>

    <section class="tab-panel" id="tab-datasets">
      <div class="panel">
        <h2>Dataset Inventory</h2>
        <div class="controls">
          <div>
            <label for="datasetFilter">Filter datasets</label>
            <input id="datasetFilter" type="text" placeholder="dataset name or histology">
          </div>
        </div>
        <div id="datasetTable"></div>
      </div>

      <div class="two">
        <div class="panel">
          <h2>Relabel Outcome: GSE213647</h2>
          <ul>
            <li>Original GEO interpretation created a false normal-vs-follicular-neoplasm issue.</li>
            <li>The rerun used the paper supplementary RNA cohort sheet instead of trusting the ambiguous GEO characteristics alone.</li>
            <li>Corrected tissue-type categories in the supplementary sheet are <code>PTC</code>, <code>Normal</code>, <code>ATC</code>, and <code>PDTC</code>.</li>
            <li>That correction improves metadata integrity, but it also removes GSE213647 from the BRAF-like vs RAS-like 2-class external validation use case.</li>
          </ul>
        </div>
        <div class="panel">
          <h2>Why GSE27155 Still Matters</h2>
          <ul>
            <li>It is not a molecularly anchored validation set.</li>
            <li>It does cover a broader histology spectrum than the tiny restricted GSE126698 subset.</li>
            <li>Its high AUC but poor threshold metrics expose how unstable a proxy-labeled external validation can look.</li>
            <li>That makes it more useful as a sanity-check cohort than as a final benchmark.</li>
          </ul>
        </div>
      </div>
    </section>

    <section class="tab-panel" id="tab-coverage">
      <div class="panel">
        <h2>Gene Coverage Explorer</h2>
        <div class="controls">
          <div>
            <label for="coverageFilter">Filter coverage rows</label>
            <input id="coverageFilter" type="text" placeholder="dataset or gene panel coverage note">
          </div>
        </div>
        <div id="coverageTable"></div>
        <div class="callout">
          Coverage was recalculated against the fixed TDS16 list, the internally derived BRS71 proxy list, and the hardcoded TierA67 composite. The GSE213647 all-zero coverage result should be treated as a matrix annotation problem until the underlying expression identifiers are revalidated.
        </div>
      </div>
    </section>

    <section class="tab-panel" id="tab-references">
      <div class="panel">
        <h2>Reference Tab</h2>
        <p>This tab lists the external papers, portals, and accession pages used or explicitly referenced in the rerun narrative. These are the places to inspect when you want to trace where a panel definition, relabel decision, or dataset context came from.</p>
        <div class="ref-grid" id="referenceGrid"></div>
      </div>
    </section>

    <div class="foot">Generated from project outputs. If you rerun the pipeline, refresh this page after regenerating <code>reports/html/index.html</code>.</div>
  </div>

  <script id="payload" type="application/json">{payload_json}</script>
  <script>
    const data = JSON.parse(document.getElementById("payload").textContent);

    function fmt(val, digits = 4) {{
      if (val === null || val === undefined || Number.isNaN(val)) return "";
      if (typeof val === "number") return val.toFixed(digits);
      return String(val);
    }}

    function tableHtml(rows, columns) {{
      const head = "<thead><tr>" + columns.map(c => `<th>${{c.label}}</th>`).join("") + "</tr></thead>";
      const body = rows.map(row => "<tr>" + columns.map(c => `<td>${{row[c.key] ?? ""}}</td>`).join("") + "</tr>").join("");
      return `<table>${{head}}<tbody>${{body}}</tbody></table>`;
    }}

    function filterRows(rows, query) {{
      if (!query) return rows;
      const q = query.toLowerCase();
      return rows.filter(row => Object.values(row).some(v => String(v ?? "").toLowerCase().includes(q)));
    }}

    function renderTable(targetId, rows, query, columns) {{
      document.getElementById(targetId).innerHTML = tableHtml(filterRows(rows, query), columns);
    }}

    function initTabs() {{
      const btns = document.querySelectorAll(".tab-btn");
      btns.forEach(btn => {{
        btn.addEventListener("click", () => {{
          btns.forEach(x => x.classList.remove("active"));
          document.querySelectorAll(".tab-panel").forEach(x => x.classList.remove("active"));
          btn.classList.add("active");
          document.getElementById(`tab-${{btn.dataset.tab}}`).classList.add("active");
        }});
      }});
    }}

    function initTables() {{
      const internalCols = [
        {{ key: "dataset", label: "dataset" }},
        {{ key: "feature_set", label: "feature_set" }},
        {{ key: "model", label: "model" }},
        {{ key: "n_samples", label: "n" }},
        {{ key: "cv_auc", label: "cv_auc" }},
        {{ key: "cv_balanced_accuracy", label: "cv_bacc" }},
        {{ key: "cv_f1", label: "cv_f1" }},
        {{ key: "cv_mcc", label: "cv_mcc" }},
        {{ key: "median_feature_count", label: "median_feature_count" }},
        {{ key: "status", label: "status" }},
      ];
      const externalCols = [
        {{ key: "dataset", label: "dataset" }},
        {{ key: "feature_set", label: "feature_set" }},
        {{ key: "model", label: "model" }},
        {{ key: "n_samples", label: "n" }},
        {{ key: "auc", label: "auc" }},
        {{ key: "balanced_accuracy", label: "bacc" }},
        {{ key: "f1", label: "f1" }},
        {{ key: "mcc", label: "mcc" }},
        {{ key: "status", label: "status" }},
      ];
      const datasetCols = [
        {{ key: "dataset", label: "dataset" }},
        {{ key: "samples", label: "samples" }},
        {{ key: "normal", label: "normal" }},
        {{ key: "tumor", label: "tumor" }},
        {{ key: "unknown_or_other", label: "unknown_or_other" }},
        {{ key: "tds16_cov", label: "TDS16_cov" }},
        {{ key: "brs71_cov", label: "BRS71_proxy_cov" }},
        {{ key: "tiera67_cov", label: "TierA67_cov" }},
        {{ key: "label_confidence_values", label: "label_confidence_values" }},
        {{ key: "histology_values", label: "histology_values" }},
      ];
      const coverageCols = [
        {{ key: "dataset", label: "dataset" }},
        {{ key: "total_gene_count", label: "total_gene_count" }},
        {{ key: "TDS16_coverage_n", label: "TDS16_n" }},
        {{ key: "TDS16_coverage_frac", label: "TDS16_frac" }},
        {{ key: "BRS71_proxy_coverage_n", label: "BRS71_n" }},
        {{ key: "BRS71_proxy_coverage_frac", label: "BRS71_frac" }},
        {{ key: "TierA67_entries_coverage_n", label: "TierA67_n" }},
        {{ key: "TierA67_entries_coverage_frac", label: "TierA67_frac" }},
      ];

      const normalizeRows = rows => rows.map(r => {{
        const out = {{ ...r }};
        Object.keys(out).forEach(k => {{
          if (typeof out[k] === "number") out[k] = Number.isInteger(out[k]) ? String(out[k]) : out[k].toFixed(4);
        }});
        return out;
      }});

      const internalRows = normalizeRows(data.internal_rows);
      const externalRows = normalizeRows(data.external_rows);
      const datasetRows = normalizeRows(data.dataset_rows);
      const coverageRows = normalizeRows(data.coverage_rows);

      const renderAll = () => {{
        renderTable("internalTable", internalRows, document.getElementById("internalFilter").value, internalCols);
        renderTable("externalTable", externalRows, document.getElementById("externalFilter").value, externalCols);
        renderTable("datasetTable", datasetRows, document.getElementById("datasetFilter").value, datasetCols);
        renderTable("coverageTable", coverageRows, document.getElementById("coverageFilter").value, coverageCols);
      }};

      ["internalFilter", "externalFilter", "datasetFilter", "coverageFilter"].forEach(id => {{
        document.getElementById(id).addEventListener("input", renderAll);
      }});

      renderAll();
    }}

    function initRoc() {{
      const rocMode = document.getElementById("rocMode");
      const rocSelect = document.getElementById("rocSelect");
      const rocFrame = document.getElementById("rocFrame");
      const rocMeta = document.getElementById("rocMeta");
      const rocAuc = document.getElementById("rocAuc");
      const rocBacc = document.getElementById("rocBacc");
      const rocF1 = document.getElementById("rocF1");
      const rocMcc = document.getElementById("rocMcc");

      function options() {{
        return rocMode.value === "internal" ? data.internal_roc_options : data.external_roc_options;
      }}

      function fillSelect() {{
        rocSelect.innerHTML = "";
        options().forEach((opt, idx) => {{
          const el = document.createElement("option");
          el.value = String(idx);
          el.textContent = opt.label;
          rocSelect.appendChild(el);
        }});
        updateRoc();
      }}

      function updateRoc() {{
        const opt = options()[Number(rocSelect.value || 0)];
        if (!opt) return;
        rocFrame.src = opt.path;
        rocAuc.textContent = fmt(opt.auc, 3);
        rocBacc.textContent = fmt(opt.bacc, 3);
        rocF1.textContent = fmt(opt.f1, 3);
        rocMcc.textContent = fmt(opt.mcc, 3);
        rocMeta.textContent = rocMode.value === "internal"
          ? `Internal CV ROC generated from out-of-fold TCGA scores. Feature set: ${{opt.feature_set}}. Model: ${{opt.model}}.`
          : `External ROC shown for the best saved HTML per dataset. Dataset: ${{opt.dataset}}. Feature set: ${{opt.feature_set}}. Model: ${{opt.model}}. n=${{opt.n_samples}}.`;
      }}

      rocMode.addEventListener("change", fillSelect);
      rocSelect.addEventListener("change", updateRoc);
      fillSelect();
    }}

    function initReferences() {{
      const grid = document.getElementById("referenceGrid");
      grid.innerHTML = data.references.map(ref => `
        <div class="ref-card">
          <div class="kicker">${{ref.group}}</div>
          <h3>${{ref.title}}</h3>
          <div class="sub">${{ref.why}}</div>
          <p style="margin-top:10px;"><a class="ref-url" href="${{ref.url}}" target="_blank" rel="noreferrer">${{ref.label}}</a></p>
          <div class="sub">${{ref.url}}</div>
        </div>
      `).join("");
    }}

    initTabs();
    initTables();
    initRoc();
    initReferences();
  </script>
</body>
</html>
"""
    (OUTDIR / "index.html").write_text(html, encoding="utf-8")


if __name__ == "__main__":
    main()
