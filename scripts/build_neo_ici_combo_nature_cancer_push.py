#!/usr/bin/env python3
"""Build the neoantigen-vaccine + ICI public-data push package.

This is a scaffolding/data-block artifact. It does not open Paper 3 Track B
claims or write manuscript prose. The goal is to lock the public combination
dataset registry and show how BioDarwin/CROSS-Neo can be attached to ICI
readiness for a Nature Cancer-style translational package.
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "project/results/p_neo_ici_combo_2026_05_11"
FIG = OUT / "figures"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")

BIODARWIN = ROOT / "project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10"
PHASE_C = ROOT / "project/results/paper11_pancancer/phase_C_ICI"
HLA18 = ROOT / "project/results/hla_deepdive_2026_05_08/track18_ici_hla"


def write_tsv(df: pd.DataFrame, name: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    df.to_csv(path, sep="\t", index=False)
    return path


def build_registry() -> pd.DataFrame:
    rows = [
        {
            "priority": 1,
            "dataset_id": "GSE222011",
            "repository": "GEO/SRA",
            "therapy": "autogene cevumeran + atezolizumab + mFOLFIRINOX",
            "cancer": "PDAC",
            "public_modality": "single-cell GEX + TCR",
            "public_size": "24 GEO samples; public scRNA/TCR subset",
            "clinical_context": "resected PDAC phase 1; sequential atezolizumab, individualized mRNA neoantigen vaccine, chemotherapy",
            "access_level": "open processed files; SRA raw available",
            "algorithm_use": "positive-control longitudinal immune-state and TCR-expansion assay lane",
            "nature_cancer_value": "direct vaccine+ICI+chemo anchor with real neoantigen-specific T cell biology",
            "limitation": "small public scRNA/TCR subset; not a stand-alone response predictor dataset",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222011",
            "direct_combo": 1,
            "open_processed": 1,
            "raw_public_or_controlled": 1,
            "sc_or_tcr": 1,
            "clinical_outcome": 1,
            "wes_rna_for_design": 0,
        },
        {
            "priority": 2,
            "dataset_id": "GSE255830",
            "repository": "GEO/dbGaP",
            "therapy": "GNOS-PV02 + INO-9012 + pembrolizumab",
            "cancer": "HCC",
            "public_modality": "single-cell T-cell GEX + TCR",
            "public_size": "8 GEO samples from 3 patients",
            "clinical_context": "advanced HCC personalized DNA neoantigen vaccine plus IL-12 plasmid and anti-PD-1",
            "access_level": "open processed GEO files; raw data withheld to dbGaP by IRB restriction",
            "algorithm_use": "post-vaccine TCR clonotype and phenotype validation lane",
            "nature_cancer_value": "cleanest public anti-PD-1 plus personalized neoantigen vaccine single-cell/TCR accession",
            "limitation": "week-12 post-treatment only; only 3 patients in GEO",
            "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255830",
            "direct_combo": 1,
            "open_processed": 1,
            "raw_public_or_controlled": 1,
            "sc_or_tcr": 1,
            "clinical_outcome": 1,
            "wes_rna_for_design": 0,
        },
        {
            "priority": 3,
            "dataset_id": "phs003922.v1.p1",
            "repository": "dbGaP",
            "therapy": "PGV001 + atezolizumab",
            "cancer": "urothelial cancer",
            "public_modality": "tumor/normal WES + tumor RNA-seq",
            "public_size": "restricted; dbGaP page reports 9 consented subjects; study text reports 12 enrolled and 10 treated",
            "clinical_context": "adjuvant or metastatic urothelial phase 1 personalized peptide vaccine plus anti-PD-L1",
            "access_level": "controlled access only",
            "algorithm_use": "best WES/RNA design-audit dataset after dbGaP approval",
            "nature_cancer_value": "patient-level manufacture, vaccine content, response, and WES/RNA bridge",
            "limitation": "not open GEO; requires dbGaP authorization",
            "url": "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003922.v1.p1",
            "direct_combo": 1,
            "open_processed": 0,
            "raw_public_or_controlled": 1,
            "sc_or_tcr": 0,
            "clinical_outcome": 1,
            "wes_rna_for_design": 1,
        },
        {
            "priority": 4,
            "dataset_id": "PRJNA1221998/PRJNA1221999",
            "repository": "ENA/BioProject",
            "therapy": "NeoVax-MI + nivolumab + locally administered ipilimumab",
            "cancer": "melanoma",
            "public_modality": "raw sequencing project metadata",
            "public_size": "10-patient Cell 2025 study; ENA project pair located",
            "clinical_context": "multi-adjuvant personalized peptide vaccine on systemic anti-PD-1 plus local CTLA-4 blockade",
            "access_level": "ENA project; practical file inventory still needed",
            "algorithm_use": "CD8 boosting and intratumoral TCR-remodeling comparator",
            "nature_cancer_value": "strong adjuvant-optimized vaccine+ICI biology and tumor TCR repertoire remodeling",
            "limitation": "not a GEO matrix; needs ENA inventory/download pass",
            "url": "https://www.omicsdi.org/dataset/project/PRJNA1221998",
            "direct_combo": 1,
            "open_processed": 0,
            "raw_public_or_controlled": 1,
            "sc_or_tcr": 1,
            "clinical_outcome": 1,
            "wes_rna_for_design": 1,
        },
        {
            "priority": 5,
            "dataset_id": "EVX-01_NCT03715985",
            "repository": "article supplements",
            "therapy": "EVX-01 + anti-PD-1/PD-L1",
            "cancer": "melanoma",
            "public_modality": "clinical table, peptide/T-cell immunomonitoring supplements",
            "public_size": "dose-escalation melanoma cohort; no GEO accession found in quick scan",
            "clinical_context": "AI-selected personalized peptide vaccine with checkpoint inhibitor backbone",
            "access_level": "supplement-level extraction",
            "algorithm_use": "industrial-style immunogenicity response benchmark, not omics reanalysis",
            "nature_cancer_value": "AI-designed vaccine comparator against BioDarwin selection logic",
            "limitation": "no open WES/RNA/scRNA matrix found",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11116868/",
            "direct_combo": 1,
            "open_processed": 0,
            "raw_public_or_controlled": 0,
            "sc_or_tcr": 0,
            "clinical_outcome": 1,
            "wes_rna_for_design": 0,
        },
        {
            "priority": 6,
            "dataset_id": "NEO-PV-01_NT-001/NT-002",
            "repository": "article/patent/supplements",
            "therapy": "NEO-PV-01 + anti-PD-1 +/- chemotherapy",
            "cancer": "melanoma, NSCLC, urothelial",
            "public_modality": "clinical and immune-monitoring tables; accession not confirmed",
            "public_size": "multi-tumor phase 1b plus NSCLC first-line combination",
            "clinical_context": "personalized long peptide vaccine plus nivolumab or pembrolizumab-based regimens",
            "access_level": "supplement-level extraction unless accession is found",
            "algorithm_use": "historical comparator and response/timing logic",
            "nature_cancer_value": "canonical vaccine+ICI translational precedent",
            "limitation": "not an immediate public matrix",
            "url": "https://doi.org/10.1016/j.cell.2020.08.053",
            "direct_combo": 1,
            "open_processed": 0,
            "raw_public_or_controlled": 0,
            "sc_or_tcr": 0,
            "clinical_outcome": 1,
            "wes_rna_for_design": 0,
        },
    ]
    df = pd.DataFrame(rows)
    metric_cols = [
        "direct_combo",
        "open_processed",
        "raw_public_or_controlled",
        "sc_or_tcr",
        "clinical_outcome",
        "wes_rna_for_design",
    ]
    df["reuse_score"] = df[metric_cols].sum(axis=1)
    return df


def build_execution_matrix() -> pd.DataFrame:
    rows = [
        {
            "stage": "A",
            "work_block": "Public combo registry lock",
            "input": "GSE222011, GSE255830, phs003922, PRJNA1221998/1999, EVX-01, NEO-PV-01",
            "output": "auditable dataset registry with access level and claim boundary",
            "decision": "GO now",
        },
        {
            "stage": "B",
            "work_block": "Open GEO single-cell/TCR reanalysis",
            "input": "GSE222011 + GSE255830 processed matrices",
            "output": "T-cell phenotype, clonotype expansion, cytotoxic/exhaustion/memory modules, vaccine timepoint plots",
            "decision": "GO next",
        },
        {
            "stage": "C",
            "work_block": "BioDarwin attach layer",
            "input": "BioDarwin v4 mode bank, CROSS-Neo, BigMHC, public candidate tables",
            "output": "algorithm-to-immune-response concordance where peptide/HLA labels are public enough",
            "decision": "GO after B",
        },
        {
            "stage": "D",
            "work_block": "ICI readiness anchor",
            "input": "Paper 11 Phase C n=421 + HLA Track18",
            "output": "tumor-level HLA-I/IFNG/checkpoint/cytolytic response gate",
            "decision": "already GREEN as external coherence, not thyroid response prediction",
        },
        {
            "stage": "E",
            "work_block": "Controlled WES/RNA expansion",
            "input": "dbGaP phs003922 PGV001 + atezo; ENA NeoVax-MI project inventory",
            "output": "patient-level vaccine design audit and clonality/HLA gate",
            "decision": "submit/access or inventory required",
        },
        {
            "stage": "F",
            "work_block": "Nature Cancer claim package",
            "input": "B-E",
            "output": "vaccine-selection algorithm plus ICI-readiness model; mechanism shown as antigen presentation and T-cell state",
            "decision": "target after open GEO reanalysis + one controlled-data plan",
        },
    ]
    return pd.DataFrame(rows)


def build_claim_boundary() -> pd.DataFrame:
    rows = [
        {
            "claim_class": "Allowed now",
            "claim": "Public neoantigen-vaccine + ICI datasets exist and can anchor a BioDarwin/CROSS-Neo translational extension.",
            "evidence": "GSE222011, GSE255830, phs003922, PRJNA1221998/1999 registry",
        },
        {
            "claim_class": "Allowed after open GEO reanalysis",
            "claim": "BioDarwin-ranked immune-response features align with vaccine-induced T-cell state or clonotype expansion in public combo datasets.",
            "evidence": "requires GSE222011/GSE255830 processed matrix analysis",
        },
        {
            "claim_class": "Allowed as external coherence",
            "claim": "ICI response biology is HLA-I/IFNG/checkpoint/cytolytic-driven across non-vaccine ICI cohorts.",
            "evidence": "Paper 11 Phase C n=421 and Track18 HLA-I OR=1.35, p=0.014",
        },
        {
            "claim_class": "Forbidden now",
            "claim": "Clinical-grade vaccine or ICI treatment selection.",
            "evidence": "small public combo cohorts and post-hoc algorithm development",
        },
        {
            "claim_class": "Forbidden now",
            "claim": "Thyroid cancer ICI response predictor.",
            "evidence": "no thyroid ICI response cohort in this package",
        },
        {
            "claim_class": "Forbidden now",
            "claim": "BioDarwin superiority over industrial platforms in patients.",
            "evidence": "needs frozen external patient-level validation and prospective assay feedback",
        },
    ]
    return pd.DataFrame(rows)


def make_figure(registry: pd.DataFrame) -> tuple[Path, Path]:
    FIG.mkdir(parents=True, exist_ok=True)
    df = registry.sort_values("reuse_score", ascending=True)
    colors = ["#64b5f6" if "GSE" in ds else "#f6c85f" if "phs" in ds else "#7fd1b9" for ds in df["dataset_id"]]
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    ax.barh(df["dataset_id"], df["reuse_score"], color=colors, edgecolor="#1f2933", linewidth=0.8)
    ax.set_xlim(0, 6.3)
    ax.set_xlabel("Reusable evidence axes (0-6)")
    ax.set_title("Neoantigen vaccine + ICI public-data ladder", fontsize=14, weight="bold")
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(
            row["reuse_score"] + 0.08,
            i,
            row["repository"],
            va="center",
            fontsize=8.5,
            color="#263238",
        )
    ax.grid(axis="x", alpha=0.22)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    png = FIG / "Fig_NeoICI_combo_dataset_ladder.png"
    pdf = FIG / "Fig_NeoICI_combo_dataset_ladder.pdf"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def html_escape(text: object) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def table_html(df: pd.DataFrame, cols: list[str]) -> str:
    head = "".join(f"<th>{html_escape(col)}</th>" for col in cols)
    rows = []
    for _, row in df.iterrows():
        cells = "".join(f"<td>{html_escape(row.get(col, ''))}</td>" for col in cols)
        rows.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def build_html(registry: pd.DataFrame, execution: pd.DataFrame, claims: pd.DataFrame, fig_path: Path) -> Path:
    top = registry.sort_values(["reuse_score", "priority"], ascending=[False, True]).head(4)
    cards = []
    for _, row in top.iterrows():
        cards.append(
            f"""
            <article class="card">
              <div class="card-kicker">{html_escape(row['repository'])}</div>
              <h3>{html_escape(row['dataset_id'])}</h3>
              <p>{html_escape(row['therapy'])}</p>
              <b>{html_escape(row['nature_cancer_value'])}</b>
            </article>
            """
        )

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Neoantigen Vaccine + ICI Nature Cancer Push</title>
  <style>
    :root {{
      --bg: #0b1117;
      --panel: #101a24;
      --line: #243447;
      --text: #e7eef7;
      --muted: #9fb1c5;
      --gold: #f6c85f;
      --blue: #64b5f6;
      --green: #7fd1b9;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font-family: Inter, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.55;
    }}
    .wrap {{ width: min(1180px, calc(100vw - 36px)); margin: 0 auto; }}
    header {{
      padding: 48px 0 28px;
      border-bottom: 1px solid var(--line);
      background: linear-gradient(180deg, #111d28 0%, #0b1117 100%);
    }}
    .kicker {{ color: var(--gold); text-transform: uppercase; letter-spacing: .08em; font-size: 12px; font-weight: 800; }}
    h1 {{ margin: 8px 0 10px; font-family: Georgia, serif; font-size: clamp(34px, 5vw, 64px); line-height: 1.02; }}
    .lead {{ max-width: 900px; color: var(--muted); font-size: 18px; }}
    .stats {{ display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 10px; margin-top: 24px; }}
    .stat {{ border: 1px solid var(--line); background: rgba(255,255,255,.03); padding: 14px; border-radius: 8px; }}
    .stat b {{ display:block; font-size: 24px; color: var(--text); }}
    .stat span {{ color: var(--muted); font-size: 12px; }}
    main {{ padding: 26px 0 56px; }}
    section {{ margin: 34px 0; }}
    h2 {{ font-size: 22px; margin: 0 0 14px; }}
    .num {{ color: var(--gold); font-family: "JetBrains Mono", monospace; margin-right: 8px; }}
    .grid {{ display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; }}
    .card {{ border: 1px solid var(--line); background: var(--panel); padding: 16px; border-radius: 8px; min-height: 175px; }}
    .card-kicker {{ color: var(--blue); font-family: "JetBrains Mono", monospace; font-size: 12px; }}
    .card h3 {{ margin: 5px 0; font-size: 20px; }}
    .card p {{ color: var(--muted); min-height: 46px; }}
    .card b {{ color: var(--green); font-weight: 700; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 13px; background: var(--panel); border: 1px solid var(--line); }}
    th, td {{ border-bottom: 1px solid var(--line); padding: 9px 10px; vertical-align: top; }}
    th {{ color: var(--gold); text-align: left; font-size: 12px; }}
    tr:last-child td {{ border-bottom: 0; }}
    .figure {{ border: 1px solid var(--line); background: #f7fafc; padding: 10px; border-radius: 8px; }}
    .figure img {{ display: block; width: 100%; height: auto; }}
    .note {{ border-left: 4px solid var(--gold); background: rgba(246,200,95,.08); padding: 14px 16px; color: #d9e2ec; }}
    code {{ color: var(--green); }}
    a {{ color: var(--blue); }}
    @media (max-width: 900px) {{
      .stats, .grid {{ grid-template-columns: 1fr 1fr; }}
    }}
    @media (max-width: 560px) {{
      .stats, .grid {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  <header>
    <div class="wrap">
      <div class="kicker">Nature Cancer push scaffold · build {datetime.now().strftime('%Y-%m-%d %H:%M')}</div>
      <h1>Neoantigen Vaccine + ICI Public Dataset Lock</h1>
      <p class="lead">BioDarwin/CROSS-Neo can be attached to real combination-therapy datasets now. The safe claim is translational prioritization: epitope selection plus HLA-I/T-cell-inflamed ICI-readiness, not clinical treatment selection.</p>
      <div class="stats">
        <div class="stat"><b>{len(registry)}</b><span>combo anchors registered</span></div>
        <div class="stat"><b>2</b><span>direct open GEO accessions</span></div>
        <div class="stat"><b>421</b><span>external ICI-response samples already scored</span></div>
        <div class="stat"><b>HLA-I</b><span>load-bearing response gate</span></div>
      </div>
    </div>
  </header>
  <main class="wrap">
    <section>
      <h2><span class="num">01</span>Top Usable Anchors</h2>
      <div class="grid">{''.join(cards)}</div>
    </section>
    <section>
      <h2><span class="num">02</span>Evidence Ladder</h2>
      <div class="figure"><img src="assets/neo_ici_combo/Fig_NeoICI_combo_dataset_ladder.png" alt="Neoantigen vaccine plus ICI dataset evidence ladder"></div>
    </section>
    <section>
      <h2><span class="num">03</span>Dataset Registry</h2>
      {table_html(registry, ["priority", "dataset_id", "repository", "therapy", "cancer", "public_modality", "access_level", "algorithm_use", "limitation", "url"])}
    </section>
    <section>
      <h2><span class="num">04</span>Execution Matrix</h2>
      {table_html(execution, ["stage", "work_block", "input", "output", "decision"])}
    </section>
    <section>
      <h2><span class="num">05</span>Claim Boundary</h2>
      {table_html(claims, ["claim_class", "claim", "evidence"])}
      <p class="note">Paper 3 Track A remains frozen. This page is an analysis scaffold for a vaccine+ICI translational package and does not claim a thyroid ICI response predictor.</p>
    </section>
    <section>
      <h2><span class="num">06</span>Local Source Paths</h2>
      <p><code>{OUT.relative_to(ROOT)}/</code></p>
      <p><code>{BIODARWIN.relative_to(ROOT)}/</code></p>
      <p><code>{PHASE_C.relative_to(ROOT)}/</code></p>
      <p><code>{HLA18.relative_to(ROOT)}/</code></p>
    </section>
  </main>
</body>
</html>
"""
    path = HUB / "neo_ici_combo_nature_cancer_push.html"
    path.write_text(html)
    return path


def build_summary(registry: pd.DataFrame, execution: pd.DataFrame, claims: pd.DataFrame) -> Path:
    top = registry.sort_values(["reuse_score", "priority"], ascending=[False, True]).head(3)
    lines = [
        "# Neoantigen vaccine + ICI Nature Cancer push scaffold",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Bottom line",
        "",
        "There are usable public combination-therapy anchors. The immediate open-GEO core is GSE222011 (PDAC autogene cevumeran + atezolizumab + chemotherapy) and GSE255830 (HCC GNOS-PV02 + pembrolizumab + IL-12). dbGaP phs003922 and ENA PRJNA1221998/PRJNA1221999 are the controlled/raw expansion lane.",
        "",
        "## Highest-priority datasets",
        "",
    ]
    for _, row in top.iterrows():
        lines.append(f"- {row['dataset_id']}: {row['therapy']} ({row['repository']}); {row['algorithm_use']}.")
    lines += [
        "",
        "## Safe claim",
        "",
        "BioDarwin/CROSS-Neo can be extended from epitope prioritization into an ICI-combination translational prioritization model by adding HLA-I integrity, IFNG/T-cell-inflamed state, checkpoint/cytolytic modules, and clonotype expansion readouts.",
        "",
        "## Hard boundary",
        "",
        "Do not claim clinical treatment selection, industrial superiority in patients, or thyroid ICI response prediction from this package.",
        "",
        "## Files",
        "",
        "- `neo_ici_combo_dataset_registry.tsv`",
        "- `nature_cancer_execution_matrix.tsv`",
        "- `neo_ici_combo_claim_boundary.tsv`",
        "- `figures/Fig_NeoICI_combo_dataset_ladder.png`",
        "- `project/papers_hub_2026_05_04/neo_ici_combo_nature_cancer_push.html`",
    ]
    path = OUT / "SUMMARY.md"
    path.write_text("\n".join(lines) + "\n")
    return path


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    registry = build_registry()
    execution = build_execution_matrix()
    claims = build_claim_boundary()

    write_tsv(registry, "neo_ici_combo_dataset_registry.tsv")
    write_tsv(execution, "nature_cancer_execution_matrix.tsv")
    write_tsv(claims, "neo_ici_combo_claim_boundary.tsv")
    fig_png, fig_pdf = make_figure(registry)
    summary = build_summary(registry, execution, claims)

    hub_asset = HUB / "assets/neo_ici_combo"
    live_asset = LIVE_HUB / "assets/neo_ici_combo"
    hub_asset.mkdir(parents=True, exist_ok=True)
    live_asset.mkdir(parents=True, exist_ok=True)
    shutil.copy2(fig_png, hub_asset / fig_png.name)
    shutil.copy2(fig_pdf, hub_asset / fig_pdf.name)
    shutil.copy2(fig_png, live_asset / fig_png.name)
    shutil.copy2(fig_pdf, live_asset / fig_pdf.name)

    html = build_html(registry, execution, claims, fig_png)
    if LIVE_HUB.exists():
        shutil.copy2(html, LIVE_HUB / html.name)

    payload = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "registry_rows": int(len(registry)),
        "open_geo": ["GSE222011", "GSE255830"],
        "controlled_or_raw_expansion": ["phs003922.v1.p1", "PRJNA1221998/PRJNA1221999"],
        "outputs": {
            "summary": str(summary),
            "html": str(html),
            "figure_png": str(fig_png),
        },
    }
    (OUT / "neo_ici_combo_summary.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
