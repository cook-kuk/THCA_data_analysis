#!/usr/bin/env python3
"""Build Paper 9 perturbation-extension registries and web report.

This script intentionally uses cached public/processed local files only. It does
not download protected or controlled datasets.
"""
from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "paper9_perturbation_2026_05_06"
WEB = ROOT / "paper9-perturbation"
RAW = ROOT / "results" / "p3_p9_full_execution" / "paper9" / "raw"
SPRINT2 = ROOT / "results" / "p3_p9_full_execution" / "paper9_sprint2_metabolic"
SPRINT1 = ROOT / "results" / "p3_p9_full_execution" / "paper9_sprint1_gls"
P9 = ROOT / "results" / "p3_p9_full_execution" / "paper9"


TARGET_ORDER = [
    "SLC1A5",
    "GLUD1",
    "GLS",
    "SLC7A5",
    "GLUD2",
    "GOT1",
    "GOT2",
    "GPT2",
    "ASNS",
    "MYC_glutamine_addiction_module",
    "LDHA_glycolysis_comparator",
    "OXPHOS_comparator",
    "NAMPT_NAD_salvage_comparator",
]

DISPLAY = {
    "MYC_glutamine_addiction_module": "MYC glutamine-addiction module",
    "LDHA_glycolysis_comparator": "LDHA/glycolysis comparator",
    "OXPHOS_comparator": "OXPHOS comparator",
    "NAMPT_NAD_salvage_comparator": "NAMPT/NAD salvage comparator",
}


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def safe_str(value) -> str:
    if pd.isna(value):
        return ""
    return str(value)


def write_tsv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_md(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.strip() + "\n")


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")


def evidence_from_local() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ranking = read_tsv(SPRINT2 / "paper9_metabolic_candidate_ranking.tsv")
    assoc = read_tsv(SPRINT2 / "paper9_module_dependency_associations.tsv")
    drug = read_tsv(SPRINT2 / "paper9_drug_concordance_summary.tsv")
    artifact = read_tsv(SPRINT2 / "paper9_artifact_control_summary.tsv")
    return ranking, assoc, drug, artifact


def build_candidate_registry() -> tuple[list[dict], dict]:
    ranking, assoc, drug, artifact = evidence_from_local()
    ranking_by = {r["candidate"]: r for _, r in ranking.iterrows()}
    assoc_by = {r["candidate"]: r for _, r in assoc.iterrows()}
    artifact_by = {r["candidate"]: r for _, r in artifact.iterrows()}
    drug_hits = drug.groupby("candidate").size().to_dict() if not drug.empty else {}

    rows = []
    for cand in TARGET_ORDER:
        rr = ranking_by.get(cand, {})
        ar = assoc_by.get(cand, {})
        af = artifact_by.get(cand, {})
        display = DISPLAY.get(cand, cand)
        if cand == "SLC1A5":
            perturb = "High: genetic KD/KO/CRISPRi, glutamine withdrawal, rescue; pharmacology is metadata-limited."
        elif cand == "GLUD1":
            perturb = "High: genetic KD/KO plus metabolic rescue; pharmacology is metadata-limited."
        elif cand == "GLS":
            perturb = "Moderate-high: genetic KD/KO and BPTES/CB-839-style inhibitor plan; local drug concordance weak."
        elif "MYC" in cand:
            perturb = "Control/secondary only: use as proliferation/artifact guard, not primary target."
        elif "comparator" in cand:
            perturb = "Comparator only: useful as pathway control, not primary Paper 9 target."
        else:
            perturb = "Watchlist: genetic perturbation feasible, lower priority unless wet-lab capacity allows."

        rows.append(
            {
                "candidate": display,
                "source_candidate_id": cand,
                "target_class": safe_str(rr.get("class", ar.get("class", ""))),
                "dependency_association": safe_str(rr.get("dependency_evidence", "")),
                "expression_association": safe_str(rr.get("expression_evidence", "")),
                "pan_cancer_corrected_association": safe_str(rr.get("pan_cancer_corrected_association", "")),
                "thyroid_only_coverage": safe_str(rr.get("thyroid_specific_coverage", "strict thyroid coverage underpowered or unavailable")),
                "drug_response_evidence": safe_str(rr.get("drug_response_evidence", f"matched local drug rows={drug_hits.get(cand, 0)}")),
                "artifact_risk": safe_str(rr.get("artifact_risk", af.get("artifact_risk", ""))),
                "claim_status": safe_str(rr.get("claim_status", "candidate_or_comparator_only")),
                "claim_boundary": safe_str(rr.get("claim_boundary", "hypothesis generation only; no validated synthetic lethality or clinical recommendation")),
                "evidence_tier": safe_str(rr.get("evidence_tier", "LOW")),
                "verdict": safe_str(rr.get("verdict", "WATCHLIST")),
                "perturbation_feasibility": perturb,
                "strict_thyroid_dependency_n": safe_str(ar.get("strict_thyroid_dependency_n", "")),
                "strict_thyroid_expression_n": safe_str(ar.get("strict_thyroid_expression_n", "")),
            }
        )

    payload = {
        "title": "Paper 9 perturbation candidate registry",
        "source_files": [
            str(SPRINT2 / "paper9_metabolic_candidate_ranking.tsv"),
            str(SPRINT2 / "paper9_module_dependency_associations.tsv"),
            str(SPRINT2 / "paper9_drug_concordance_summary.tsv"),
            str(SPRINT2 / "paper9_artifact_control_summary.tsv"),
        ],
        "claim_boundary": "Aggressive prioritization is allowed; validated synthetic lethality, clinical recommendation, and patient-selection biomarker claims are forbidden.",
        "rows": rows,
    }
    return rows, payload


def build_modality_matrix() -> list[dict]:
    specs = {
        "SLC1A5": {
            "target_class": "glutamine transporter / ASCT2",
            "genetic": "siRNA; shRNA; CRISPRi; CRISPR KO",
            "pharmacologic": "V-9302; GPNA if available in local/public metadata; do not claim validated thyroid efficacy",
            "nutrient": "glutamine withdrawal; low-glutamine media titration",
            "rescue": "glutamine rescue; alpha-ketoglutarate rescue; nonessential amino-acid rescue control",
            "readouts": "viability; apoptosis; proliferation; colony formation; TDS/lineage score; SLC5A5/TPO/TSHR/TG/PAX8/NKX2-1 qPCR or RNA-seq",
            "role": "primary wet-lab priority",
        },
        "GLUD1": {
            "target_class": "glutamate dehydrogenase",
            "genetic": "siRNA; shRNA; CRISPRi; CRISPR KO",
            "pharmacologic": "R162; EGCG-like GLUD1 inhibition if source and specificity are confirmed",
            "nutrient": "glutamine withdrawal plus glutamate/TCA stress context",
            "rescue": "alpha-ketoglutarate rescue; TCA intermediate rescue",
            "readouts": "viability; alpha-ketoglutarate; glutamate; OCR/ECAR; TCA/OXPHOS stress; lineage score",
            "role": "primary wet-lab priority",
        },
        "GLS": {
            "target_class": "glutaminase",
            "genetic": "siRNA; shRNA; CRISPRi; CRISPR KO",
            "pharmacologic": "BPTES; CB-839/telaglenastat if available; local BPTES concordance is weak",
            "nutrient": "glutamine withdrawal",
            "rescue": "alpha-ketoglutarate rescue",
            "readouts": "viability; glutamine addiction; apoptosis; lineage/TDS genes",
            "role": "candidate, not central proven target",
        },
        "MYC glutamine-addiction module": {
            "target_class": "proliferation/glutamine-addiction module",
            "genetic": "do not primary-target as Paper 9 mechanism; track MYC and cell-cycle controls",
            "pharmacologic": "not primary; use MYC/proteasome/mTOR-linked compounds only as artifact controls",
            "nutrient": "glutamine withdrawal as shared stress control",
            "rescue": "rescue should separate metabolic dependence from generic proliferation stress",
            "readouts": "MYC; LDHA; MKI67; cell-cycle score; common-essential proxies",
            "role": "secondary/artifact guard",
        },
    }
    rows = []
    for cand, spec in specs.items():
        rows.append({"candidate": cand, **spec, "claim_boundary": "hypothesis-generating perturbation plan only"})
    return rows


def infer_target(term: str, name: str) -> str:
    text = f"{term} {name}".lower()
    if any(x in text for x in ["slc1a5", "asct2", "v-9302", "gpna"]):
        return "SLC1A5/ASCT2 glutamine transport"
    if any(x in text for x in ["glud1", "r162", "egcg"]):
        return "GLUD1/glutamate dehydrogenase"
    if any(x in text for x in ["gls", "glutaminase", "bptes", "cb-839", "telaglenastat"]):
        return "GLS/glutaminase"
    if "ldha" in text:
        return "LDHA/glycolysis comparator"
    if any(x in text for x in ["nampt", "fk866"]):
        return "NAMPT/NAD salvage comparator"
    if any(x in text for x in ["oligomycin", "rotenone", "oxphos"]):
        return "OXPHOS comparator"
    if re.search(r"(^|[^a-z0-9])myc([^a-z0-9]|$)", text) or "bortezomib" in text:
        return "MYC/proliferation artifact control"
    return "other searched perturbagen"


def metadata_term_mask(series: pd.Series, term: str) -> pd.Series:
    if term.upper() == "MYC":
        return series.str.contains(r"(?:^|[^A-Za-z0-9])MYC(?:[^A-Za-z0-9]|$)", case=False, regex=True, na=False)
    return series.str.contains(re.escape(term), case=False, na=False)


def build_drug_registry() -> list[dict]:
    terms = [
        "SLC1A5", "ASCT2", "V-9302", "GPNA", "GLUD1", "R162", "EGCG", "GLS",
        "glutaminase", "BPTES", "CB-839", "telaglenastat", "glutamine",
        "alpha-ketoglutarate", "MYC", "LDHA", "NAMPT", "FK866", "OXPHOS",
        "oligomycin", "rotenone",
    ]
    rows: list[dict] = []

    # Processed association summaries already encode direction-normalized sensitivity.
    processed = read_tsv(SPRINT2 / "paper9_drug_concordance_summary.tsv")
    for _, r in processed.iterrows():
        label = safe_str(r.get("compound_label", ""))
        rows.append(
            {
                "source": safe_str(r.get("source", "")),
                "compound_name": label.split(" DPC-")[0] if label and label != "No metadata match by Sprint 2 regex" else label,
                "compound_id": safe_str(r.get("compound_id", "")),
                "target_annotation": infer_target(safe_str(r.get("candidate", "")), label),
                "dose_time": "not in processed association table",
                "matched_cell_line_count": safe_str(r.get("n_models", "")),
                "thyroid_model_count": safe_str(r.get("strict_thyroid_n", "")),
                "usable_for_paper9": "yes_for_sensitivity_context" if safe_str(r.get("compound_id", "")) != "none_found" else "no_metadata_match",
                "direction_interpretable": "yes; already direction-normalized to sensitivity in Sprint 2 summary" if safe_str(r.get("compound_id", "")) != "none_found" else "no",
                "caveat": "weak/discordant if p not supportive; thyroid subset is underpowered",
                "matched_search_term": safe_str(r.get("candidate", "")),
            }
        )

    # Raw condition/treatment metadata search. This is metadata only, not effect interpretation.
    metadata_files = [
        ("PRISM24Q2_treatment_metadata", RAW / "Repurposing_Public_24Q2_Treatment_Meta_Data.csv", "name"),
        ("GDSC1_conditions", RAW / "GDSC1Log2ViabilityCollapsedConditions.csv", "CompoundName"),
        ("GDSC2_conditions", RAW / "GDSC2Log2ViabilityCollapsedConditions.csv", "CompoundName"),
        ("CTRP_conditions", RAW / "CTRPLog2ViabilityCollapsedConditions.csv", "CompoundName"),
    ]
    seen = {(r["source"], r["compound_id"], r["compound_name"]) for r in rows}
    for source, path, name_col in metadata_files:
        if not path.exists():
            rows.append(
                {
                    "source": source,
                    "compound_name": "metadata file not found",
                    "compound_id": "",
                    "target_annotation": "",
                    "dose_time": "",
                    "matched_cell_line_count": "",
                    "thyroid_model_count": "",
                    "usable_for_paper9": "blocked",
                    "direction_interpretable": "no",
                    "caveat": f"missing local file: {path}",
                    "matched_search_term": "",
                }
            )
            continue
        df = pd.read_csv(path)
        name_series = df[name_col].astype(str)
        for term in terms:
            mask = metadata_term_mask(name_series, term)
            if not mask.any():
                continue
            hit = df.loc[mask].copy()
            id_col = "broad_id" if "broad_id" in hit.columns else "CompoundID"
            sample_col = "SampleID" if "SampleID" in hit.columns else None
            for (cid, cname), sub in hit.groupby([id_col, name_col], dropna=False):
                dose = ""
                if "Dose" in sub.columns:
                    dose = f"{sub['Dose'].min()}-{sub['Dose'].max()} {safe_str(sub.get('DoseUnit', pd.Series([''])).iloc[0])}".strip()
                elif "dose" in sub.columns:
                    dose = f"{sub['dose'].min()}-{sub['dose'].max()}"
                key = (source, safe_str(cid), safe_str(cname))
                if key in seen:
                    continue
                seen.add(key)
                rows.append(
                    {
                        "source": source,
                        "compound_name": safe_str(cname),
                        "compound_id": safe_str(cid),
                        "target_annotation": infer_target(term, safe_str(cname)),
                        "dose_time": dose or "metadata has no dose/time field",
                        "matched_cell_line_count": str(sub[sample_col].nunique()) if sample_col else "not available from treatment metadata alone",
                        "thyroid_model_count": "not available from treatment metadata alone",
                        "usable_for_paper9": "metadata_only; effect direction requires matrix convention audit",
                        "direction_interpretable": "no; metadata match only",
                        "caveat": "raw metadata hit only; do not infer sensitivity without response matrix and convention check",
                        "matched_search_term": term,
                    }
                )
    rows = sorted(rows, key=lambda x: (x["source"], x["target_annotation"], x["compound_name"]))
    return rows


def build_lincs_inputs() -> dict:
    query_dir = OUT / "paper9_lincs_query_inputs"
    query_dir.mkdir(parents=True, exist_ok=True)
    up = [
        "VIM", "MYC", "LDHA", "SLC1A5", "GLUD1", "GLS", "SLC7A5", "ASNS",
        "DNMT1", "FN1", "TACSTD2",
    ]
    down = ["SLC5A5", "TPO", "TSHR", "TG", "PAX8", "NKX2-1", "DIO1", "EPCAM"]
    (query_dir / "paper9_lincs_query_up_genes.txt").write_text("\n".join(up) + "\n")
    (query_dir / "paper9_lincs_query_down_genes.txt").write_text("\n".join(down) + "\n")
    payload = {
        "signature_name": "paper9_lineage_silenced_thyroid_glutamine_axis",
        "up_genes": up,
        "down_genes": down,
        "interpretation": {
            "reversal": "candidate perturbagens may restore thyroid differentiation-like state",
            "mimic": "candidate perturbagens may induce lineage-silenced/metabolic-stress-like state",
            "boundary": "hypothesis generation only; filter toxicity/common-essential/proliferation artifacts",
        },
        "required_gene_check": [
            "SLC5A5", "TPO", "TSHR", "TG", "PAX8", "NKX2-1", "DIO1", "EPCAM",
            "VIM", "MYC", "LDHA", "SLC1A5", "GLUD1", "GLS",
        ],
    }
    write_json(query_dir / "paper9_lincs_query_signature.json", payload)
    write_md(
        query_dir / "paper9_lincs_query_instructions.md",
        """
# Paper 9 LINCS/CMap Query Instructions

Use the up/down gene lists as a query-ready lineage-silenced thyroid cancer signature.

Interpretation:
- Reversal hits are candidates for restoring a thyroid differentiation-like state.
- Mimic hits are candidates for inducing or marking a lineage-silenced/metabolic-stress state.
- All hits must be filtered for toxicity, common-essential biology, MYC/cell-cycle artifacts, and thyroid model relevance.

Claim boundary:
- These inputs do not validate a therapy.
- These inputs do not define patient selection.
- These inputs are hypothesis-generation inputs for perturbation follow-up.
""",
    )
    return payload


def build_perturbseq_plan() -> list[dict]:
    local_files = [
        ROOT / "notebooks_or_scripts" / "v6_gpu" / "20_geneformer.py",
        ROOT / "notebooks_or_scripts" / "v6_gpu" / "30_gears_perturb.py",
        ROOT / "results" / "nature_cancer_feasibility" / "gse232237_target_gene_log1p_cpm10k.tsv",
        ROOT / "results" / "v17_gse241184" / "per_sample_summary.tsv",
    ]
    availability = "; ".join([f"{p.name}:{'present' if p.exists() else 'missing'}" for p in local_files])
    rows = []
    for target in ["SLC1A5", "GLUD1", "GLS"]:
        rows.append(
            {
                "resource_or_workflow": "Replogle-style Perturb-seq metadata/signature search",
                "target": target,
                "local_status": "no local target-level Perturb-seq signature table found",
                "action": "query public metadata for target perturbation signatures; map DE genes to TDS genes and glutamine-axis genes if available",
                "claim_status": "planned only; no perturb-seq validation claim",
            }
        )
        rows.append(
            {
                "resource_or_workflow": "Geneformer/GEARS local workflow feasibility",
                "target": target,
                "local_status": availability,
                "action": "use only if existing workflow can run on processed malignant thyrocyte data; otherwise keep as feasibility plan",
                "claim_status": "simulation feasibility only; no wet-lab substitute",
            }
        )
    return rows


def build_wetlab_matrix() -> list[dict]:
    rows = [
        {
            "model_state": "lineage-high thyroid cancer model",
            "subtype_note": "PTC/DTC-like if available",
            "perturbation": "non-targeting control; glutamine withdrawal control",
            "candidate": "baseline comparator",
            "readouts": "viability; apoptosis; proliferation; TDS genes",
            "go_criterion": "maintains differentiation genes better than lineage-silenced comparator",
            "falsification": "no state separation or severe nonspecific stress",
        },
        {
            "model_state": "lineage-silenced thyroid cancer model",
            "subtype_note": "ATC/PDTC/PTC-derived low-TDS model if available",
            "perturbation": "SLC1A5 KD/KO/CRISPRi; V-9302/GPNA if suitable",
            "candidate": "SLC1A5",
            "readouts": "viability; apoptosis; colony formation; SLC5A5/TPO/TSHR/TG/PAX8/NKX2-1; MYC/LDHA/SLC1A5",
            "go_criterion": "selective vulnerability in lineage-silenced models plus rescue by glutamine or alpha-ketoglutarate",
            "falsification": "equal toxicity in lineage-high models or no rescue",
        },
        {
            "model_state": "lineage-silenced thyroid cancer model",
            "subtype_note": "ATC/PDTC/PTC-derived low-TDS model if available",
            "perturbation": "GLUD1 KD/KO/CRISPRi; R162/EGCG-like inhibitor only if specificity acceptable",
            "candidate": "GLUD1",
            "readouts": "viability; alpha-ketoglutarate; glutamate; OCR/ECAR; TCA rescue; lineage score",
            "go_criterion": "state-selective vulnerability and alpha-ketoglutarate/TCA rescue support metabolic mechanism",
            "falsification": "no viability effect, no metabolic shift, or nonspecific toxicity",
        },
        {
            "model_state": "lineage-silenced thyroid cancer model",
            "subtype_note": "ATC/PDTC/PTC-derived low-TDS model if available",
            "perturbation": "GLS KD/KO/CRISPRi; BPTES; CB-839/telaglenastat if obtainable",
            "candidate": "GLS",
            "readouts": "viability; glutamine-addiction response; apoptosis; alpha-ketoglutarate rescue; TDS genes",
            "go_criterion": "concordant genetic and pharmacologic effect with rescue",
            "falsification": "BPTES/drug effect discordant with GLS dependency or no rescue",
        },
        {
            "model_state": "all tested models",
            "subtype_note": "artifact guard",
            "perturbation": "MYC/proliferation artifact controls; LDHA/OXPHOS/NAMPT comparators",
            "candidate": "artifact/comparator panel",
            "readouts": "MYC; LDHA; MKI67; cell-cycle score; common-essential proxy behavior",
            "go_criterion": "primary candidate effects are separable from generic proliferation collapse",
            "falsification": "signal collapses into common-essential/proliferation biology",
        },
    ]
    return rows


def render_candidate_rows(rows: list[dict]) -> str:
    primary = [r for r in rows if r["candidate"] in ["SLC1A5", "GLUD1", "GLS", "SLC7A5", "GLUD2", "GOT1", "GOT2", "GPT2", "ASNS", "MYC glutamine-addiction module", "LDHA/glycolysis comparator", "OXPHOS comparator", "NAMPT/NAD salvage comparator"]]
    html = []
    for r in primary:
        html.append(
            f"<tr><td><strong>{r['candidate']}</strong></td><td>{r['target_class']}</td><td>{r['dependency_association']}<br>{r['pan_cancer_corrected_association']}</td><td>{r['perturbation_feasibility']}</td><td>{r['expression_association']}</td><td>{r['artifact_risk']}</td><td>{r['verdict']}</td></tr>"
        )
    return "\n".join(html)


def render_wetlab_rows(rows: list[dict]) -> str:
    return "\n".join(
        f"<tr><td>{r['candidate']}</td><td>{r['model_state']}</td><td>{r['perturbation']}</td><td>{r['readouts']}</td><td>{r['go_criterion']}</td><td>{r['falsification']}</td></tr>"
        for r in rows
    )


def build_web_page(candidate_rows: list[dict], wetlab_rows: list[dict], drug_rows: list[dict]) -> None:
    WEB.mkdir(parents=True, exist_ok=True)
    top_drug_hits = [r for r in drug_rows if r["usable_for_paper9"].startswith("yes")][:8]
    if not top_drug_hits:
        top_drug_hits = drug_rows[:8]
    drug_html = "\n".join(
        f"<tr><td>{r['source']}</td><td>{r['compound_name']}</td><td>{r['target_annotation']}</td><td>{r['matched_cell_line_count']}</td><td>{r['thyroid_model_count']}</td><td>{r['usable_for_paper9']}</td><td>{r['caveat']}</td></tr>"
        for r in top_drug_hits
    )
    html = f"""<!doctype html>
<html lang="ko">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Paper 9 Perturbation Extension</title>
<style>
:root{{--ink:#102033;--muted:#546274;--paper:#f6f1e8;--card:#fffdf8;--line:#d8ccbb;--red:#9a2e2e;--redbg:#fff0f0;--blue:#295b84;--bluebg:#eaf3fb;--green:#3d6a50;--greenbg:#edf7ef;--gold:#a86e18;--goldbg:#fff5dc;}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--paper);color:var(--ink);font-family:Arial,'Noto Sans KR',sans-serif;line-height:1.58}} a{{color:var(--blue);text-decoration:none}}
.wrap{{max-width:1220px;margin:auto;padding:24px 24px 72px}} .nav{{display:flex;justify-content:space-between;gap:12px;align-items:center;margin-bottom:18px;border-bottom:1px solid var(--line);padding-bottom:12px;font-size:13px}}
.hero{{background:linear-gradient(135deg,#fffdf8,#ecf4f8);border:1px solid var(--line);border-left:10px solid var(--blue);border-radius:8px;padding:34px 38px;margin-bottom:22px}} .kicker{{font-size:12px;text-transform:uppercase;letter-spacing:.12em;color:var(--blue);font-weight:800}} h1{{font-size:42px;line-height:1.08;margin:10px 0 10px}} .subtitle{{font-size:18px;color:var(--muted);max-width:960px}} .badge{{display:inline-block;margin-top:14px;background:var(--red);color:#fff;border-radius:999px;padding:8px 12px;font-weight:800;font-size:13px}}
section{{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:24px;margin-top:18px}} h2{{font-size:28px;line-height:1.15;margin:0 0 14px}} h3{{font-size:18px;margin:0 0 8px;color:var(--blue)}} .grid{{display:grid;gap:14px}} .g3{{grid-template-columns:repeat(3,minmax(0,1fr))}} .g2{{grid-template-columns:repeat(2,minmax(0,1fr))}}
.card{{background:#fff;border:1px solid var(--line);border-radius:8px;padding:16px}} .green{{border-left:6px solid var(--green);background:var(--greenbg)}} .blue{{border-left:6px solid var(--blue);background:var(--bluebg)}} .red{{border:2px solid var(--red);background:var(--redbg)}} .gold{{border-left:6px solid var(--gold);background:var(--goldbg)}}
.verdict{{font-size:18px;border-left:7px solid var(--green);background:#fff;padding:18px;border-radius:0 8px 8px 0}} .flow{{display:grid;grid-template-columns:repeat(7,minmax(0,1fr));gap:10px}} .step{{background:#fff;border:1px solid var(--line);border-top:4px solid var(--blue);border-radius:8px;padding:12px;font-size:13px;min-height:98px}} .step b{{display:block;color:var(--blue);margin-bottom:5px}}
table{{width:100%;border-collapse:collapse;background:#fff;border:1px solid var(--line);font-size:13px}} th{{background:#142238;color:#fff;text-align:left;padding:9px;font-size:11px;text-transform:uppercase;letter-spacing:.04em}} td{{border-top:1px solid #eadfce;padding:9px;vertical-align:top}} ul{{padding-left:19px}} li{{margin:5px 0}} .kr{{font-size:17px;line-height:1.78;background:#fff;border:2px solid var(--blue);border-radius:8px;padding:20px}} .small{{font-size:13px;color:var(--muted)}} .tablewrap{{overflow:auto}}
@media(max-width:900px){{.wrap{{padding:14px}} h1{{font-size:30px}} .hero{{padding:24px 18px}} .g2,.g3,.flow{{grid-template-columns:1fr}} table{{font-size:12px}}}}
</style>
</head>
<body>
<div class="wrap">
<nav class="nav"><a href="/papers_hub_2026_05_04/index.html">8-Papers Hub</a><span>Paper 9 perturbation roadmap</span><a href="/paper9">Current Paper 9</a></nav>
<header class="hero">
  <div class="kicker">Paper 9 perturbation extension</div>
  <h1>Paper 9 Perturbation Extension: Glutamine-Axis Vulnerability Validation Map</h1>
  <p class="subtitle">Perturbation roadmap for glutamine-axis vulnerabilities in lineage-silenced thyroid cancer.</p>
  <span class="badge">Exploratory translational roadmap — not validated synthetic lethality</span>
</header>

<section>
  <h2>Executive Verdict</h2>
  <div class="verdict"><strong>GO as a perturbation roadmap.</strong> The strongest allowed claim is that public-data dependency and expression evidence prioritize a testable glutamine-axis wet-lab queue in lineage-silenced thyroid cancer. The current order is <strong>SLC1A5 &gt; GLUD1 &gt; GLS</strong>. This is not a validated synthetic-lethality claim, not a clinical treatment recommendation, and not a patient-selection biomarker.</div>
  <div class="grid g3" style="margin-top:14px">
    <div class="card green"><h3>Allowed</h3><p>Candidate prioritization, perturbation roadmap, wet-lab assay design, and falsifiable hypothesis generation.</p></div>
    <div class="card gold"><h3>Still weak</h3><p>Thyroid-only dependency coverage is underpowered, and processed drug response coverage is sparse for SLC1A5/GLUD1.</p></div>
    <div class="card red"><h3>Forbidden</h3><p>Validated synthetic lethality, clinical recommendation, patient selection, and causal therapeutic target discovery.</p></div>
  </div>
</section>

<section>
  <h2>Why Perturbation Matters</h2>
  <div class="grid g3">
    <div class="card"><h3>Correlation to intervention</h3><p>Dependency correlations nominate candidates; perturbation tests whether the lineage-silenced state creates an actionable experimental vulnerability.</p></div>
    <div class="card"><h3>Wet-lab plan</h3><p>The ranking becomes concrete: KD/KO, inhibitor where feasible, glutamine withdrawal, rescue, and lineage/metabolic readouts.</p></div>
    <div class="card"><h3>Falsifiability</h3><p>The hypothesis fails if effects are nonspecific, not stronger in lineage-silenced models, or not rescued by glutamine/TCA pathway rescue.</p></div>
  </div>
</section>

<section>
  <h2>Candidate Perturbation Ranking</h2>
  <div class="tablewrap"><table>
    <thead><tr><th>Candidate</th><th>Target class</th><th>Current evidence</th><th>Perturbation modality</th><th>Expected readout</th><th>Artifact risk</th><th>Verdict</th></tr></thead>
    <tbody>{render_candidate_rows(candidate_rows)}</tbody>
  </table></div>
</section>

<section class="red">
  <h2>Claim Boundary</h2>
  <p><strong>Forbidden:</strong> validated synthetic lethality, clinical recommendation, patient selection, causal therapeutic claim, or drug-response-based treatment proposal. Local evidence supports prioritization only. Negative and weak evidence is retained: GLS BPTES concordance is weak/discordant, SLC1A5 and GLUD1 lack matched processed drug metadata, and strict thyroid CRISPR coverage is small.</p>
</section>

<section>
  <h2>Perturbation Workflow</h2>
  <div class="flow">
    <div class="step"><b>1. Lineage-silenced score</b>Classify lineage-high vs lineage-silenced thyroid models.</div>
    <div class="step"><b>2. Candidate ranking</b>SLC1A5, GLUD1, GLS, then watchlist/comparators.</div>
    <div class="step"><b>3. Perturbagen registry</b>Check cached PRISM/GDSC/CTRP/LINCS metadata.</div>
    <div class="step"><b>4. Genetic/drug perturbation</b>KD/KO/CRISPRi and inhibitor where metadata supports feasibility.</div>
    <div class="step"><b>5. Rescue assay</b>Glutamine or alpha-ketoglutarate rescue.</div>
    <div class="step"><b>6. Readout</b>Viability, apoptosis, OCR/ECAR, TDS genes.</div>
    <div class="step"><b>7. Go/no-go</b>Advance only if state-selective and rescued.</div>
  </div>
</section>

<section>
  <h2>Wet-Lab Validation Matrix</h2>
  <div class="tablewrap"><table>
    <thead><tr><th>Candidate</th><th>Model</th><th>Perturbation</th><th>Readouts</th><th>Go criterion</th><th>Falsification</th></tr></thead>
    <tbody>{render_wetlab_rows(wetlab_rows)}</tbody>
  </table></div>
</section>

<section>
  <h2>Drug/Perturbagen Metadata Snapshot</h2>
  <p class="small">Directionality rule: for AUC data, higher raw AUC generally means resistance; Sprint 2 summaries use direction-normalized sensitivity. For PRISM LFC, lower LFC means stronger depletion/sensitivity, but no SLC1A5/GLUD1/GLS matched PRISM metadata was found in the processed pass.</p>
  <div class="tablewrap"><table>
    <thead><tr><th>Source</th><th>Compound</th><th>Target annotation</th><th>Cell-line count</th><th>Thyroid count</th><th>Use</th><th>Caveat</th></tr></thead>
    <tbody>{drug_html}</tbody>
  </table></div>
</section>

<section>
  <h2>Professor Summary</h2>
  <div class="kr">교수님, Paper 9는 이제 단순 상관분석에서 perturbation 기반 검증 설계로 확장할 수 있습니다. 현재 우선순위는 SLC1A5, GLUD1, GLS 순서이며, 핵심은 lineage-silenced thyroid cancer 모델에서 knockdown, inhibitor, glutamine withdrawal, rescue assay를 통해 실제 취약성이 존재하는지 확인하는 것입니다. 단, 현재 단계에서는 validated synthetic lethality가 아니라 실험 검증 후보 우선순위화로 해석해야 합니다.</div>
</section>

<section>
  <h2>Missing Evidence</h2>
  <div class="card gold">
    <ul>
      <li>Thyroid-specific perturbation data for SLC1A5, GLUD1, and GLS.</li>
      <li>Matched drug-response coverage for validated SLC1A5/GLUD1 inhibitors in thyroid models.</li>
      <li>Genetic and pharmacologic concordance in lineage-high versus lineage-silenced thyroid cancer models.</li>
      <li>Rescue assay with glutamine or alpha-ketoglutarate.</li>
      <li>Single-cell perturbation validation or query-confirmed Perturb-seq signatures.</li>
      <li>Artifact control against MYC, LDHA, common-essential, and generic proliferation effects.</li>
    </ul>
  </div>
</section>

<section>
  <h2>Generated Local Outputs</h2>
  <p class="small">All tables and plans were written under <code>project/results/paper9_perturbation_2026_05_06/</code>. The analysis uses cached local Paper 9/DepMap/PRISM/GDSC/CTRP-derived files and does not use protected data.</p>
</section>
</div>
</body>
</html>
"""
    (WEB / "index.html").write_text(html)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    candidate_rows, candidate_payload = build_candidate_registry()
    write_tsv(candidate_rows, OUT / "paper9_perturbation_candidate_registry.tsv")
    write_json(OUT / "paper9_perturbation_candidate_registry.json", candidate_payload)

    write_md(
        OUT / "paper9_perturbation_evidence_summary.md",
        """
# Paper 9 Perturbation Evidence Summary

Paper 9 is feasible as a perturbation-focused validation roadmap for glutamine-axis vulnerabilities in lineage-silenced thyroid cancer.

Already supported locally:
- SLC1A5 and GLUD1 are the top wet-lab priorities in Sprint 2 metabolic ranking.
- GLS has supportive dependency/expression evidence but weak/discordant BPTES drug concordance, so it remains a candidate rather than the headline.
- MYC glutamine-addiction module is secondary and artifact-guarded.
- LDHA/OXPHOS/NAMPT are comparators only.
- Strict thyroid model dependency coverage is underpowered, so thyroid-specific inference remains weak.

Claim boundary:
- Allowed: metabolic candidate prioritization and perturbation design.
- Forbidden: validated synthetic lethality, therapeutic target discovery, patient selection biomarker, or clinical recommendation.
""",
    )

    modality_rows = build_modality_matrix()
    write_tsv(modality_rows, OUT / "paper9_perturbation_modality_matrix.tsv")
    write_md(
        OUT / "paper9_perturbation_modality_matrix.md",
        "# Paper 9 Perturbation Modality Matrix\n\n"
        + "\n".join(
            f"## {r['candidate']}\n- Genetic: {r['genetic']}\n- Pharmacologic: {r['pharmacologic']}\n- Nutrient: {r['nutrient']}\n- Rescue: {r['rescue']}\n- Readouts: {r['readouts']}\n- Role: {r['role']}\n- Boundary: {r['claim_boundary']}\n"
            for r in modality_rows
        ),
    )

    drug_rows = build_drug_registry()
    write_tsv(drug_rows, OUT / "paper9_drug_perturbagen_registry.tsv")
    write_md(
        OUT / "paper9_drug_directionality_audit.md",
        """
# Paper 9 Drug Directionality Audit

Local directionality rules:
- AUC matrices: higher raw AUC generally indicates less sensitivity/resistance. Paper 9 Sprint 2 uses direction-normalized sensitivity for association summaries.
- LFC matrices: lower LFC indicates stronger depletion/sensitivity. Direction must not be inferred from treatment metadata alone.
- Metadata-only hits are not interpretable as drug response.

Local result:
- GDSC1 BPTES is available for GLS, but the lineage-sensitivity association is weak/discrepant in the processed summary.
- PRISM 24Q2 treatment metadata did not yield processed SLC1A5, GLUD1, or GLS inhibitor matches in the prior Paper 9 pass.
- SLC1A5/GLUD1 pharmacologic validation remains wet-lab feasibility, not public drug-response support.
""",
    )
    write_md(
        OUT / "paper9_drug_perturbation_go_no_go.md",
        """
# Paper 9 Drug/Perturbation Go/No-Go

GO:
- Use SLC1A5 and GLUD1 as genetic perturbation priorities.
- Use GLS as a candidate with BPTES/CB-839-style perturbation only after acknowledging weak public drug concordance.
- Use glutamine withdrawal and rescue assays as core falsification tests.

NO-GO:
- Do not claim validated synthetic lethality.
- Do not claim clinical treatment recommendation.
- Do not treat metadata-only drug matches as sensitivity evidence.
- Do not make SLC1A5/GLUD1 drug claims without matched thyroid model response or wet-lab testing.
""",
    )

    lincs_payload = build_lincs_inputs()
    write_md(
        OUT / "paper9_lincs_cmap_perturbation_plan.md",
        """
# Paper 9 LINCS/CMap Perturbation Plan

Local CMap-style reversal files exist, but the perturbation extension requires a glutamine-axis lineage-silenced query rather than immediate therapeutic claims.

Query construction:
- Up genes represent lineage-silenced/metabolic-stress genes.
- Down genes represent thyroid differentiation genes and lineage-high genes.

Interpretation:
- Reversal perturbagens may restore a thyroid differentiation-like state.
- Mimic perturbagens may induce lineage-silenced/metabolic-stress state.
- Candidate perturbagens must be filtered by toxicity, common-essential and proliferation artifacts.

Boundary:
- Hypothesis generation only.
- No clinical recommendation or patient-selection claim.
""",
    )

    perturbseq_rows = build_perturbseq_plan()
    write_tsv(perturbseq_rows, OUT / "paper9_perturbseq_registry_or_plan.tsv")
    write_md(
        OUT / "paper9_in_silico_perturbation_plan.md",
        """
# Paper 9 In Silico Perturbation Plan

Local workflow signals:
- Geneformer and GEARS-style script names exist locally.
- Processed thyroid single-cell and pseudobulk files exist locally.
- No completed SLC1A5/GLUD1/GLS in silico perturbation result was found.

Plan:
1. Use processed malignant thyrocyte states if available.
2. Compare lineage-high versus lineage-silenced malignant thyrocytes.
3. Simulate perturbation only for SLC1A5, GLUD1, and GLS if the existing workflow can be run reproducibly.
4. Score effects on SLC5A5, TPO, TSHR, TG, PAX8, NKX2-1, DIO1, MYC, LDHA, SLC1A5, GLUD1, and GLS.

Boundary:
- In silico perturbation is a prioritization aid, not validation.
""",
    )
    write_md(
        OUT / "paper9_single_cell_perturbation_claim_boundary.md",
        """
# Paper 9 Single-Cell Perturbation Claim Boundary

Allowed:
- State mapping between lineage-high and lineage-silenced malignant thyrocytes.
- Query plan for Perturb-seq or in silico perturbation.
- Hypothesis that perturbing SLC1A5/GLUD1/GLS may change lineage/metabolic readouts.

Forbidden:
- Claiming perturbation validation without actual perturbation data.
- Treating Geneformer/GEARS simulation as wet-lab evidence.
- Inferring clinical response or patient selection.
""",
    )

    wetlab_rows = build_wetlab_matrix()
    write_tsv(wetlab_rows, OUT / "paper9_wetlab_perturbation_matrix.tsv")
    write_md(
        OUT / "paper9_wetlab_validation_plan.md",
        """
# Paper 9 Wet-Lab Validation Plan

Core design:
- Select lineage-high and lineage-silenced thyroid cancer models where available.
- Prioritize SLC1A5 KD/KO, GLUD1 KD/KO, GLS KD/KO, GLS inhibitor, glutamine withdrawal, and rescue.
- Include MYC/proliferation artifact controls.

Readouts:
- Cell viability, apoptosis, proliferation, colony formation.
- OCR/ECAR if Seahorse is available.
- qPCR or RNA-seq for SLC5A5, TPO, TSHR, TG, PAX8, NKX2-1, MYC, LDHA, SLC1A5, GLUD1, GLS.
- Isotope tracing is optional.

Go criterion:
- Lineage-silenced models show stronger candidate-specific vulnerability and rescue with glutamine or alpha-ketoglutarate.

No-go criterion:
- Effects are nonspecific, not state-selective, not rescued, or collapse into MYC/common-essential proliferation artifacts.
""",
    )
    write_md(
        OUT / "paper9_professor_perturbation_summary_kr.md",
        """
# 교수님 보고용 요약

교수님, Paper 9는 이제 단순 상관분석에서 perturbation 기반 검증 설계로 확장할 수 있습니다. 현재 우선순위는 SLC1A5, GLUD1, GLS 순서이며, 핵심은 lineage-silenced thyroid cancer 모델에서 knockdown, inhibitor, glutamine withdrawal, rescue assay를 통해 실제 취약성이 존재하는지 확인하는 것입니다. 단, 현재 단계에서는 validated synthetic lethality가 아니라 실험 검증 후보 우선순위화로 해석해야 합니다.

현재 가장 강한 claim은 "lineage-silenced thyroid cancer에서 glutamine-axis perturbation 후보를 우선순위화했다"입니다. 아직 말하면 안 되는 claim은 validated synthetic lethality, clinical recommendation, patient selection biomarker입니다.
""",
    )

    build_web_page(candidate_rows, wetlab_rows, drug_rows)

    print(json.dumps({
        "output_dir": str(OUT),
        "web_page": str(WEB / "index.html"),
        "n_candidate_rows": len(candidate_rows),
        "n_drug_rows": len(drug_rows),
        "lincs_signature": lincs_payload["signature_name"],
    }, indent=2))


if __name__ == "__main__":
    main()
