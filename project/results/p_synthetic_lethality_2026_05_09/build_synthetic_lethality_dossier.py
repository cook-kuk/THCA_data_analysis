#!/usr/bin/env python3
"""Build the synthetic-lethality claim-boundary dossier from local tables.

The dossier intentionally uses conservative language: local public data support
candidate vulnerability prioritization, not validated synthetic lethality.
"""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "project/results/p_synthetic_lethality_2026_05_09"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_REL = Path("assets/paper1/synthetic_lethality")
HUB_ASSETS = HUB / ASSET_REL
LIVE_ASSETS = LIVE / ASSET_REL


SOURCE = {
    "prism": ROOT / "project/results/p_deconv_2026_05_08/v15C_prism_annotated_drugs.tsv",
    "prism_enrich": ROOT / "project/results/p_deconv_2026_05_08/v16_prism_enrichment_sensitivity.tsv",
    "prism_lineage": ROOT / "project/results/p_deconv_2026_05_08/v16_prism_lineage_sensitivity.tsv",
    "depmap": ROOT / "project/results/p_deconv_2026_05_08/v15C_depmap_dependency_overlay.tsv",
    "metabolic": ROOT / "project/results/p3_p9_full_execution/paper9_sprint2_metabolic/paper9_module_dependency_associations.tsv",
    "paper9_gene": ROOT / "project/results/p3_p9_full_execution/paper9/p9_candidate_dependency_associations.tsv",
    "paper9_go": ROOT / "project/results/paper9_perturbation_2026_05_06/paper9_drug_perturbation_go_no_go.md",
    "paper9_summary": ROOT / "project/results/paper9_perturbation_2026_05_06/paper9_perturbation_evidence_summary.md",
}

FIGURES = {
    "prism_top_drugs.png": ROOT / "project/results/dm1_robustness_v2026_05_08/round6/prism_top_drugs.png",
    "depmap_dependencies.png": ROOT / "project/results/dm1_robustness_v2026_05_08/round6/depmap_dependencies.png",
    "paper9_lineage_dependency_heatmap.png": ROOT / "project/results/p3_p9_full_execution/paper9/fig_p9_lineage_dependency_heatmap.png",
    "v16_spatial_prism_diagnostics.png": ROOT / "project/results/p_deconv_2026_05_08/Fig_SX_v16_spatial_prism_diagnostics.png",
}


def fmt(x: float, digits: int = 3) -> str:
    if pd.isna(x):
        return "NA"
    if x == 0:
        return "0"
    ax = abs(float(x))
    if ax < 0.001:
        return f"{float(x):.2e}"
    return f"{float(x):.{digits}f}"


def fdr(x: float) -> str:
    if pd.isna(x):
        return "NA"
    return f"{float(x):.2e}" if float(x) < 0.001 else f"{float(x):.3f}"


def read_tables() -> dict[str, pd.DataFrame]:
    tables = {}
    for key, path in SOURCE.items():
        if path.suffix == ".tsv":
            tables[key] = pd.read_csv(path, sep="\t")
    return tables


def get_row(df: pd.DataFrame, col: str, value: str) -> pd.Series:
    hit = df.loc[df[col].astype(str) == value]
    if hit.empty:
        raise ValueError(f"Missing row {value!r} in {col}")
    return hit.iloc[0]


def build_summary(tables: dict[str, pd.DataFrame]) -> dict[str, object]:
    prism = tables["prism"]
    enrich = tables["prism_enrich"]
    lineage = tables["prism_lineage"]
    depmap = tables["depmap"]
    metabolic = tables["metabolic"]

    fdr_hits = prism.loc[prism["fdr"] < 0.05].copy()
    top7 = get_row(enrich, "set", "top_7")
    top11 = get_row(enrich, "set", "top_11")
    top15 = get_row(enrich, "set", "top_15")
    fdr_row = get_row(enrich, "set", "fdr_lt_0.05")
    all_nonmissing = get_row(lineage, "group", "all_nonmissing")
    exclude_thyroid = get_row(lineage, "group", "exclude_thyroid")
    azd = get_row(prism, "name", "AZD-0364")
    myc = get_row(depmap, "gene", "MYC")
    nampt = get_row(depmap, "gene", "NAMPT")

    metabolic_rows = {}
    for candidate in ["SLC1A5", "GLUD1", "GLS"]:
        row = get_row(metabolic, "candidate", candidate)
        metabolic_rows[candidate] = {
            "dependency_rho": float(row["dependency_spearman"]),
            "dependency_p": float(row["dependency_p"]),
            "dependency_n": int(row["dependency_n"]),
            "pan_corrected_dependency_rho": float(row["pan_cancer_corrected_dependency_spearman"]),
            "pan_corrected_dependency_p": float(row["pan_cancer_corrected_dependency_p"]),
            "strict_thyroid_dependency_rho": float(row["strict_thyroid_dependency_spearman"]),
            "strict_thyroid_dependency_p": float(row["strict_thyroid_dependency_p"]),
            "strict_thyroid_dependency_n": int(row["strict_thyroid_dependency_n"]),
        }

    return {
        "prism": {
            "n_drugs": int(prism.shape[0]),
            "fdr_hits": int(fdr_hits.shape[0]),
            "fdr_canonical_mapk": int(fdr_hits["is_canonical_mapk"].sum()),
            "fdr_mapk_pi3k": int(fdr_hits["is_mapk_pi3k_axis"].sum()),
            "top7_canonical_mapk": f"{int(top7['n_canonical_mapk'])}/{int(top7['n_hit'])}",
            "top11_canonical_mapk": f"{int(top11['n_canonical_mapk'])}/{int(top11['n_hit'])}",
            "top15_canonical_mapk": f"{int(top15['n_canonical_mapk'])}/{int(top15['n_hit'])}",
            "fdr_enrichment_p": float(fdr_row["canonical_hypergeom_p"]),
            "azd_d": float(azd["cohens_d"]),
            "azd_delta_lfc": float(azd["delta_LFC"]),
            "azd_fdr": float(azd["fdr"]),
            "rho_all": float(all_nonmissing["rho_dm1_vs_mapk_lfc"]),
            "rho_all_p": float(all_nonmissing["p"]),
            "rho_exclude_thyroid": float(exclude_thyroid["rho_dm1_vs_mapk_lfc"]),
            "rho_exclude_thyroid_p": float(exclude_thyroid["p"]),
            "high_low_mapk_lfc_d": float(all_nonmissing["high_vs_low_cohens_d_lfc"]),
            "mapk_lfc_n": int(all_nonmissing["n"]),
        },
        "depmap": {
            "n_high": int(myc["n_high"]),
            "n_low": int(myc["n_low"]),
            "myc_d": float(myc["cohens_d"]),
            "myc_p": float(myc["p"]),
            "nampt_d": float(nampt["cohens_d"]),
            "nampt_p": float(nampt["p"]),
        },
        "metabolic": metabolic_rows,
    }


def build_matrix(summary: dict[str, object]) -> pd.DataFrame:
    p = summary["prism"]
    d = summary["depmap"]
    m = summary["metabolic"]

    rows = [
        {
            "tier": "T1",
            "candidate_program": "DM1-high / MAPK-axis pharmacologic vulnerability",
            "state_or_context": "DM1-high, lineage-silenced pan-cancer cell-line state",
            "perturbation": "MEK/RAF/ERK inhibitors",
            "headline_metrics": (
                f"PRISM {p['n_drugs']} drugs; FDR<0.05 hits {p['fdr_hits']}, "
                f"canonical MAPK {p['fdr_canonical_mapk']}/{p['fdr_hits']}; "
                f"top7 MAPK {p['top7_canonical_mapk']}; AZD-0364 d={fmt(p['azd_d'])}, "
                f"FDR={fdr(p['azd_fdr'])}; DM1 score x MAPK LFC rho={fmt(p['rho_all'])}, "
                f"p={fdr(p['rho_all_p'])}"
            ),
            "claim_status": "Allowed: pharmacologic vulnerability / actionability reserve. Forbidden: validated synthetic lethality.",
            "disposition": "KEEP AS REVIEWER-RESERVE, NOT MAIN SYNTHETIC-LETHAL CLAIM",
            "required_next_validation": "Matched thyroid/DM1-high models, dose-response MEK/RAF/ERK, genetic MAPK rescue or bypass, and lineage-gene restoration readout.",
            "source_paths": "v15C_prism_annotated_drugs.tsv; v16_prism_enrichment_sensitivity.tsv; v16_prism_lineage_sensitivity.tsv",
        },
        {
            "tier": "T2",
            "candidate_program": "DM1-high CRISPR dependency",
            "state_or_context": "DM1-high vs DM1-low DepMap cell lines",
            "perturbation": "MYC and NAMPT genetic perturbation; inhibitor follow-up only as separate validation",
            "headline_metrics": (
                f"n_high={d['n_high']}, n_low={d['n_low']}; MYC d={fmt(d['myc_d'])}, "
                f"p={fdr(d['myc_p'])}; NAMPT d={fmt(d['nampt_d'])}, p={fdr(d['nampt_p'])}"
            ),
            "claim_status": "Allowed: dependency prioritization. Forbidden: selective clinical target without pan-essentiality controls.",
            "disposition": "KEEP AS FUNCTIONAL PRIORITIZATION; MYC HAS PAN-ESSENTIAL RISK",
            "required_next_validation": "Within-lineage matched models, CRISPRi/siRNA titration, rescue, viability plus apoptosis/colony assays, and toxicity window against lineage-high controls.",
            "source_paths": "v15C_depmap_dependency_overlay.tsv; phase_D_depmap/dm1_high_low_dependency.tsv",
        },
        {
            "tier": "T3A",
            "candidate_program": "Glutamine transporter vulnerability",
            "state_or_context": "lineage-silenced thyroid cancer model",
            "perturbation": "SLC1A5 KD/KO/CRISPRi; glutamine withdrawal; V-9302/GPNA only if metadata supports",
            "headline_metrics": (
                f"SLC1A5 dependency rho={fmt(m['SLC1A5']['dependency_rho'])}, "
                f"p={fdr(m['SLC1A5']['dependency_p'])}, n={m['SLC1A5']['dependency_n']}; "
                f"pan-cancer corrected rho={fmt(m['SLC1A5']['pan_corrected_dependency_rho'])}, "
                f"p={fdr(m['SLC1A5']['pan_corrected_dependency_p'])}; strict thyroid n={m['SLC1A5']['strict_thyroid_dependency_n']} NS"
            ),
            "claim_status": "Allowed: wet-lab priority. Forbidden: drug-response claim without matched response data.",
            "disposition": "TOP PAPER9 WET-LAB PRIORITY",
            "required_next_validation": "Genetic perturbation plus glutamine or alpha-ketoglutarate rescue; viability/apoptosis/colony; SLC5A5/TPO/TSHR/TG/PAX8/NKX2-1 readout.",
            "source_paths": "paper9_module_dependency_associations.tsv; paper9_drug_perturbation_go_no_go.md",
        },
        {
            "tier": "T3B",
            "candidate_program": "Glutamate/TCA anaplerosis vulnerability",
            "state_or_context": "lineage-silenced thyroid cancer model",
            "perturbation": "GLUD1 KD/KO/CRISPRi; inhibitor only if specificity acceptable",
            "headline_metrics": (
                f"GLUD1 dependency rho={fmt(m['GLUD1']['dependency_rho'])}, "
                f"p={fdr(m['GLUD1']['dependency_p'])}; pan-cancer corrected rho={fmt(m['GLUD1']['pan_corrected_dependency_rho'])}, "
                f"p={fdr(m['GLUD1']['pan_corrected_dependency_p'])}; strict thyroid n={m['GLUD1']['strict_thyroid_dependency_n']} NS"
            ),
            "claim_status": "Allowed: wet-lab priority. Forbidden: inhibitor claim from weak metadata.",
            "disposition": "SECOND PAPER9 WET-LAB PRIORITY",
            "required_next_validation": "GLUD1 perturbation with alpha-ketoglutarate/TCA rescue; OCR/ECAR; viability and lineage-gene readout.",
            "source_paths": "paper9_module_dependency_associations.tsv; paper9_drug_perturbation_go_no_go.md",
        },
        {
            "tier": "T3C",
            "candidate_program": "Glutaminase vulnerability",
            "state_or_context": "lineage-silenced thyroid cancer model",
            "perturbation": "GLS KD/KO/CRISPRi; BPTES/CB-839 style pharmacology only with caveat",
            "headline_metrics": (
                f"GLS dependency rho={fmt(m['GLS']['dependency_rho'])}, "
                f"p={fdr(m['GLS']['dependency_p'])}; pan-cancer corrected rho={fmt(m['GLS']['pan_corrected_dependency_rho'])}, "
                f"p={fdr(m['GLS']['pan_corrected_dependency_p'])}; strict thyroid n={m['GLS']['strict_thyroid_dependency_n']} NS"
            ),
            "claim_status": "Allowed: candidate. Forbidden: headline due to weak/discordant public drug concordance.",
            "disposition": "CANDIDATE, NOT HEADLINE",
            "required_next_validation": "Concordant genetic and pharmacologic effect with rescue; resolve BPTES/CB-839 concordance before using as main story.",
            "source_paths": "paper9_module_dependency_associations.tsv; paper9_perturbation_evidence_summary.md",
        },
        {
            "tier": "BOUNDARY",
            "candidate_program": "Spatial MAPK-panel test",
            "state_or_context": "GSE250521 thyroid Visium",
            "perturbation": "None; observational spatial stress test",
            "headline_metrics": "v16/v17 did not rescue spatial anti-correlation; 0/45 subset/adjustment tests negative in v16; full+detection residual rho ~0 in v17.",
            "claim_status": "Allowed: caveat/reserve. Forbidden: mechanism support or synthetic lethality support.",
            "disposition": "USE AS HONEST NEGATIVE CAVEAT ONLY",
            "required_next_validation": "Do not use for synthetic-lethality argument; separate from pharmacologic/CRISPR prioritization.",
            "source_paths": "v16_spatial_prism_diagnostics_summary.json; v17_spatial_signal_decomposition_summary.json",
        },
        {
            "tier": "BOUNDARY",
            "candidate_program": "ICI vulnerability / immune phenotype",
            "state_or_context": "DM1 inflamed/dedifferentiated tumor state",
            "perturbation": "Immune checkpoint therapy context",
            "headline_metrics": "Project boundary map already separates immune vulnerability from synthetic-lethal and synthetic-rescue hypotheses.",
            "claim_status": "Allowed: separate Paper 3 immune-vulnerability frame. Forbidden: calling ICI response synthetic lethality.",
            "disposition": "KEEP OUT OF SYNTHETIC-LETHALITY CLAIM SET",
            "required_next_validation": "Treat as separate Paper 3 track; do not merge with DM1 perturbation vulnerability language.",
            "source_paths": "project/reports/2026_05_06_ici_vs_synthetic_lethality_boundary_map.md",
        },
    ]
    return pd.DataFrame(rows)


def copy_assets() -> None:
    HUB_ASSETS.mkdir(parents=True, exist_ok=True)
    LIVE_ASSETS.mkdir(parents=True, exist_ok=True)
    for name, src in FIGURES.items():
        if src.exists():
            shutil.copy2(src, HUB_ASSETS / name)
            shutil.copy2(src, LIVE_ASSETS / name)


def markdown_table(df: pd.DataFrame, columns: list[str]) -> str:
    rows = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for _, row in df.iterrows():
        vals = []
        for col in columns:
            val = str(row[col]).replace("\n", " ").replace("|", "/")
            vals.append(val)
        rows.append("| " + " | ".join(vals) + " |")
    return "\n".join(rows)


def build_markdown(summary: dict[str, object], matrix: pd.DataFrame) -> str:
    p = summary["prism"]
    d = summary["depmap"]
    m = summary["metabolic"]
    return f"""# Synthetic lethality 정리 — 2026-05-09

## Bottom line

현재 로컬 근거로 가능한 가장 강한 표현은 **validated synthetic lethality**가 아니라 **DM1/lineage-silenced state에서의 perturbation vulnerability 후보 우선순위화**입니다. Paper 1에는 reviewer-reserve/actionability block으로만 두고, 본문 핵심 claim으로 올리면 과합니다.

## Claim boundary

- **Allowed:** pharmacologic vulnerability, CRISPR dependency prioritization, wet-lab perturbation roadmap, reviewer-reserve actionability.
- **Forbidden:** validated synthetic lethality, clinical treatment recommendation, patient-selection biomarker, treatment-induced thyroid-gene restoration.
- **Best one-line wording:** DM1-high or lineage-silenced models show public-data vulnerability signals, especially MAPK-axis PRISM sensitivity and MYC/NAMPT dependency, motivating focused perturbation experiments rather than establishing synthetic lethality.

## Highest-signal facts

- PRISM: {p['n_drugs']} drugs, FDR<0.05 hits {p['fdr_hits']}; canonical MAPK hits {p['fdr_canonical_mapk']}/{p['fdr_hits']}; top 7 by FDR are {p['top7_canonical_mapk']} canonical MAPK.
- Top PRISM hit: AZD-0364, Cohen's d={fmt(p['azd_d'])}, delta LFC={fmt(p['azd_delta_lfc'])}, FDR={fdr(p['azd_fdr'])}.
- Continuous PRISM MAPK score: DM1 score vs canonical MAPK-inhibitor mean LFC rho={fmt(p['rho_all'])}, p={fdr(p['rho_all_p'])}, n={p['mapk_lfc_n']}; excluding thyroid rho={fmt(p['rho_exclude_thyroid'])}, p={fdr(p['rho_exclude_thyroid_p'])}.
- DepMap CRISPR: MYC d={fmt(d['myc_d'])}, p={fdr(d['myc_p'])}; NAMPT d={fmt(d['nampt_d'])}, p={fdr(d['nampt_p'])}; n_high={d['n_high']}, n_low={d['n_low']}.
- Paper9 wet-lab queue: SLC1A5 > GLUD1 > GLS. SLC1A5 rho={fmt(m['SLC1A5']['dependency_rho'])}, p={fdr(m['SLC1A5']['dependency_p'])}; GLUD1 rho={fmt(m['GLUD1']['dependency_rho'])}, p={fdr(m['GLUD1']['dependency_p'])}; GLS rho={fmt(m['GLS']['dependency_rho'])}, p={fdr(m['GLS']['dependency_p'])}.

## Candidate matrix

{markdown_table(matrix, ['tier', 'candidate_program', 'perturbation', 'headline_metrics', 'disposition'])}

## Deployment rule

Use this as a **factual data block / dossier**, not as voice-protected Discussion, Limitations, Cover Letter, Hook, Aim, or Q9 prose.
"""


def html_escape(text: object) -> str:
    s = str(text)
    return (
        s.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def figure_card(filename: str, title: str, caption: str) -> str:
    return f"""
    <figure class="fig-card">
      <img src="{ASSET_REL.as_posix()}/{filename}" alt="{html_escape(title)}">
      <figcaption><strong>{html_escape(title)}</strong><br>{html_escape(caption)}</figcaption>
    </figure>
    """


def build_html(summary: dict[str, object], matrix: pd.DataFrame) -> str:
    p = summary["prism"]
    d = summary["depmap"]
    cards = [
        ("PRISM drugs", f"{p['fdr_canonical_mapk']}/{p['fdr_hits']} FDR hits canonical MAPK"),
        ("AZD-0364", f"d {fmt(p['azd_d'])}; FDR {fdr(p['azd_fdr'])}"),
        ("MAPK LFC", f"rho {fmt(p['rho_all'])}; p {fdr(p['rho_all_p'])}"),
        ("MYC", f"d {fmt(d['myc_d'])}; p {fdr(d['myc_p'])}"),
        ("NAMPT", f"d {fmt(d['nampt_d'])}; p {fdr(d['nampt_p'])}"),
        ("Verdict", "candidate prioritization, not validated SL"),
    ]
    stat_html = "\n".join(
        f"<div class=\"stat\"><span>{html_escape(k)}</span><b>{html_escape(v)}</b></div>" for k, v in cards
    )

    matrix_rows = "\n".join(
        "<tr>"
        + "".join(
            f"<td class=\"{ 'gold' if col == 'disposition' else '' }\">{html_escape(row[col])}</td>"
            for col in [
                "tier",
                "candidate_program",
                "state_or_context",
                "perturbation",
                "headline_metrics",
                "claim_status",
                "disposition",
                "required_next_validation",
                "source_paths",
            ]
        )
        + "</tr>"
        for _, row in matrix.iterrows()
    )

    source_rows = "\n".join(
        f"<tr><td>{html_escape(key)}</td><td><code>{html_escape(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path)}</code></td></tr>"
        for key, path in SOURCE.items()
    )

    figs = "\n".join(
        [
            figure_card(
                "prism_top_drugs.png",
                "PRISM top drugs",
                "Existing round6 figure showing DM1-high drug-sensitivity ranking.",
            ),
            figure_card(
                "depmap_dependencies.png",
                "DepMap dependencies",
                "Existing round6 figure summarizing DM1-high CRISPR dependency candidates.",
            ),
            figure_card(
                "paper9_lineage_dependency_heatmap.png",
                "Paper9 perturbation heatmap",
                "Local perturbation roadmap evidence used for SLC1A5 / GLUD1 / GLS triage.",
            ),
            figure_card(
                "v16_spatial_prism_diagnostics.png",
                "Spatial/PRISM diagnostics",
                "Useful caveat: spatial MAPK-panel anti-correlation was not rescued; PRISM wording was corrected.",
            ),
        ]
    )

    sibling_cards = "\n".join(
        [
            '<a class="page-card" href="paper1_fig8_mechanism_dossier.html"><b>Fig 8 Mechanism Dossier</b><span>v6-v17 mechanism evidence and caveats</span></a>',
            '<a class="page-card" href="paper1_reviewer_defense_dashboard.html"><b>Reviewer Defense Dashboard</b><span>Q1-Q14 factual answer surface</span></a>',
            '<a class="page-card" href="paper1_methods_reproducibility_dossier.html"><b>Methods Reproducibility</b><span>gene lists, formulas, filters, edge cases</span></a>',
            '<a class="page-card" href="paper1_spatial_signal_decomposition_v17.html"><b>Spatial Signal Decomposition</b><span>negative/caveat autopsy for GSE250521</span></a>',
        ]
    )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Paper 1 Synthetic-Lethality Claim Boundary</title>
  <style>
    :root {{
      --bg: #111318;
      --panel: #191d25;
      --panel-2: #202631;
      --ink: #ece7da;
      --muted: #a7aebd;
      --line: #343b49;
      --gold: #d9ad5f;
      --red: #ef6b62;
      --green: #86d391;
      --blue: #8ab7ff;
      --mono: "JetBrains Mono", ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      --serif: "Cormorant Garamond", Georgia, serif;
      --sans: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    * {{ box-sizing: border-box; }}
    html {{ scroll-behavior: smooth; }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font-family: var(--sans);
      line-height: 1.55;
    }}
    a {{ color: var(--blue); text-decoration: none; }}
    code {{ font-family: var(--mono); color: #f2d395; font-size: .9em; }}
    .hero {{
      min-height: 58vh;
      padding: 58px max(28px, 7vw) 34px;
      background:
        linear-gradient(90deg, rgba(17,19,24,.96), rgba(17,19,24,.76)),
        url("{ASSET_REL.as_posix()}/prism_top_drugs.png") center/cover no-repeat;
      border-bottom: 1px solid var(--line);
      display: flex;
      flex-direction: column;
      justify-content: flex-end;
    }}
    .kicker {{ font: 700 12px var(--mono); color: var(--gold); letter-spacing: 0; text-transform: uppercase; }}
    h1 {{
      font-family: var(--serif);
      font-size: clamp(42px, 7vw, 92px);
      line-height: .92;
      letter-spacing: 0;
      margin: 12px 0 18px;
      max-width: 1050px;
    }}
    .lead {{ max-width: 980px; color: #d8dde8; font-size: 19px; }}
    .crumbs {{ margin-top: 16px; color: var(--muted); font-family: var(--mono); font-size: 12px; }}
    .stats {{ display: grid; grid-template-columns: repeat(6, minmax(140px, 1fr)); gap: 10px; margin-top: 28px; max-width: 1250px; }}
    .stat {{ background: rgba(25,29,37,.88); border: 1px solid var(--line); padding: 14px; border-radius: 8px; min-height: 86px; }}
    .stat span {{ display:block; color: var(--muted); font: 700 11px var(--mono); text-transform: uppercase; }}
    .stat b {{ display:block; margin-top: 8px; font-size: 20px; color: var(--ink); }}
    .layout {{ display: grid; grid-template-columns: 270px minmax(0, 1fr); gap: 28px; max-width: 1540px; margin: 0 auto; padding: 30px 26px 70px; }}
    nav {{ position: sticky; top: 12px; align-self: start; max-height: calc(100vh - 24px); overflow: auto; border-right: 1px solid var(--line); padding-right: 18px; }}
    nav b {{ display:block; color: var(--gold); font-family: var(--mono); font-size: 12px; margin-bottom: 12px; }}
    nav a {{ display:block; color: var(--muted); padding: 8px 0; border-bottom: 1px solid rgba(52,59,73,.55); font-size: 14px; }}
    section {{ margin: 0 0 42px; }}
    h2 {{ font-family: var(--serif); font-size: 38px; line-height: 1.05; margin: 0 0 16px; letter-spacing: 0; }}
    h2 .num {{ color: var(--gold); font-family: var(--mono); font-size: 18px; margin-right: 10px; }}
    .tl {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-left: 5px solid var(--gold);
      border-radius: 8px;
      padding: 18px 20px;
      font-size: 18px;
    }}
    .warn {{ border-left-color: var(--red); }}
    .ok {{ border-left-color: var(--green); }}
    .grid {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 16px; }}
    .fig-card {{ margin: 0; background: var(--panel); border: 1px solid var(--line); border-radius: 8px; overflow: hidden; }}
    .fig-card img {{ width: 100%; display: block; background: #fff; }}
    figcaption {{ padding: 12px 14px; color: var(--muted); font-size: 13px; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 13px; background: var(--panel); border: 1px solid var(--line); }}
    th, td {{ border-bottom: 1px solid var(--line); vertical-align: top; padding: 10px 12px; }}
    th {{ color: var(--gold); font-family: var(--mono); text-align: left; font-size: 11px; text-transform: uppercase; background: #151922; }}
    td.gold {{ background: rgba(217,173,95,.12); color: #f0d9a7; font-weight: 700; }}
    .table-wrap {{ overflow-x: auto; border-radius: 8px; }}
    .claim-grid {{ display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 14px; }}
    .claim {{ background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 16px; }}
    .claim b {{ display: block; color: var(--gold); margin-bottom: 8px; }}
    .claim .bad {{ color: var(--red); }}
    .claim .good {{ color: var(--green); }}
    .page-grid {{ display:grid; grid-template-columns: repeat(4, minmax(0,1fr)); gap: 12px; }}
    .page-card {{ display:block; min-height: 112px; background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 15px; color: var(--ink); }}
    .page-card span {{ display:block; color: var(--muted); margin-top: 8px; font-size: 13px; }}
    .small {{ color: var(--muted); font-size: 13px; }}
    @media (max-width: 980px) {{
      .stats {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
      .layout {{ grid-template-columns: 1fr; }}
      nav {{ position: static; border-right: 0; border-bottom: 1px solid var(--line); padding: 0 0 14px; }}
      .grid, .claim-grid, .page-grid {{ grid-template-columns: 1fr; }}
      h1 {{ font-size: 46px; }}
    }}
  </style>
</head>
<body>
  <header class="hero">
    <div class="kicker">Paper 1 / Paper 9 bridge · Claim-boundary dossier · 2026-05-09</div>
    <h1>Synthetic-Lethality Boundary</h1>
    <p class="lead">Local PRISM, DepMap, and Paper9 perturbation outputs support a focused vulnerability queue. They do not yet establish validated synthetic lethality. This page keeps the usable actionability signal and the forbidden overclaim in the same frame.</p>
    <div class="crumbs">Source of truth: <code>project/results/p_synthetic_lethality_2026_05_09/</code></div>
    <div class="stats">{stat_html}</div>
  </header>
  <div class="layout">
    <nav>
      <b>Contents</b>
      <a href="#tl">01 TL;DR</a>
      <a href="#boundary">02 Claim Boundary</a>
      <a href="#matrix">03 Decision Matrix</a>
      <a href="#figures">04 Figures</a>
      <a href="#wetlab">05 Wet-Lab Queue</a>
      <a href="#caveats">06 Caveats</a>
      <a href="#siblings">07 Related Pages</a>
      <a href="#sources">08 Sources</a>
    </nav>
    <main>
      <section id="tl">
        <h2><span class="num">01</span>TL;DR</h2>
        <div class="tl"><strong>Current verdict:</strong> strongest allowed claim is candidate vulnerability prioritization in DM1-high / lineage-silenced models. PRISM strongly prioritizes MAPK-axis inhibitors; DepMap prioritizes MYC/NAMPT; Paper9 prioritizes SLC1A5 &gt; GLUD1 &gt; GLS for wet-lab perturbation. None is validated synthetic lethality yet.</div>
      </section>
      <section id="boundary">
        <h2><span class="num">02</span>Claim Boundary</h2>
        <div class="claim-grid">
          <div class="claim"><b class="good">Allowed</b>Pharmacologic vulnerability, CRISPR dependency prioritization, wet-lab perturbation roadmap, reviewer-reserve actionability.</div>
          <div class="claim"><b class="bad">Forbidden</b>Validated synthetic lethality, clinical recommendation, patient-selection biomarker, treatment-induced thyroid-gene restoration.</div>
          <div class="claim"><b>Safe Wording</b>DM1-high or lineage-silenced models show public-data vulnerability signals motivating focused perturbation experiments.</div>
        </div>
      </section>
      <section id="matrix">
        <h2><span class="num">03</span>Decision Matrix</h2>
        <div class="table-wrap">
          <table>
            <thead>
              <tr><th>Tier</th><th>Candidate</th><th>Context</th><th>Perturbation</th><th>Headline metrics</th><th>Claim status</th><th>Disposition</th><th>Next validation</th><th>Sources</th></tr>
            </thead>
            <tbody>{matrix_rows}</tbody>
          </table>
        </div>
      </section>
      <section id="figures">
        <h2><span class="num">04</span>Figures</h2>
        <div class="grid">{figs}</div>
      </section>
      <section id="wetlab">
        <h2><span class="num">05</span>Wet-Lab Queue</h2>
        <div class="tl ok"><strong>Priority order:</strong> SLC1A5 &gt; GLUD1 &gt; GLS. The falsification design should pair genetic perturbation with glutamine or alpha-ketoglutarate rescue, viability/apoptosis/colony phenotypes, and thyroid-lineage readouts.</div>
      </section>
      <section id="caveats">
        <h2><span class="num">06</span>Caveats</h2>
        <div class="tl warn">Spatial GSE250521 does not support a local within-spot MAPK-panel anti-correlation and should remain a caveat/reserve result. ICI vulnerability belongs to Paper 3 immune framing, not synthetic lethality. TROP2/topo-I ADC work is a separate actionability track and should not be merged into this claim set.</div>
      </section>
      <section id="siblings">
        <h2><span class="num">07</span>Related Pages</h2>
        <div class="page-grid">{sibling_cards}</div>
      </section>
      <section id="sources">
        <h2><span class="num">08</span>Sources</h2>
        <div class="table-wrap">
          <table>
            <thead><tr><th>Key</th><th>Path</th></tr></thead>
            <tbody>{source_rows}</tbody>
          </table>
        </div>
        <p class="small">Generated by <code>project/results/p_synthetic_lethality_2026_05_09/build_synthetic_lethality_dossier.py</code>. Live copy deployed to <code>/var/www/papers/papers_hub_2026_05_04/paper1_synthetic_lethality_dossier.html</code>.</p>
      </section>
    </main>
  </div>
</body>
</html>
"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    copy_assets()

    tables = read_tables()
    summary = build_summary(tables)
    matrix = build_matrix(summary)

    matrix_path = OUT / "synthetic_lethality_candidate_matrix.tsv"
    summary_path = OUT / "synthetic_lethality_summary.json"
    brief_path = OUT / "SYNTHETIC_LETHALITY_BRIEF_KR.md"
    html_path = HUB / "paper1_synthetic_lethality_dossier.html"
    live_html_path = LIVE / "paper1_synthetic_lethality_dossier.html"

    matrix.to_csv(matrix_path, sep="\t", index=False)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    brief_path.write_text(build_markdown(summary, matrix), encoding="utf-8")
    html_path.write_text(build_html(summary, matrix), encoding="utf-8")
    shutil.copy2(html_path, live_html_path)

    print(f"Wrote {matrix_path.relative_to(ROOT)}")
    print(f"Wrote {summary_path.relative_to(ROOT)}")
    print(f"Wrote {brief_path.relative_to(ROOT)}")
    print(f"Wrote {html_path.relative_to(ROOT)}")
    print(f"Deployed {live_html_path}")


if __name__ == "__main__":
    main()
