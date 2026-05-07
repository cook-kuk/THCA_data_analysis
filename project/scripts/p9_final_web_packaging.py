#!/usr/bin/env python3
"""Final Paper 9 web packaging and sanity notes.

Uses existing Sprint 1/2 outputs only. No large downloads.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
S2 = ROOT / "project/results/p3_p9_full_execution/paper9_sprint2_metabolic"
S1 = ROOT / "project/results/p3_p9_full_execution/paper9_sprint1_gls"
RAW = ROOT / "project/results/p3_p9_full_execution/paper9/raw"
REPORT = ROOT / "project/reports/p3_p9_full_execution"

KR_SUMMARY = """교수님, Paper 9 Sprint 2 결과를 보면 GLS 단일 synthetic lethality 논문으로 가는 것은 위험하지만, lineage-silenced thyroid cancer의 glutamine-axis metabolic vulnerability를 우선순위화하는 논문 방향은 가능성이 있습니다.

GLS는 dependency와 expression association은 관찰되지만 BPTES drug response concordance가 약해서 단독 validated target으로 주장하기 어렵습니다. 반면 SLC1A5와 GLUD1은 cancer-type correction 후에도 dependency association이 유지되고, 특히 SLC1A5는 glutamine transporter로서 lineage-silenced state와 가장 안정적으로 연결되어 wet-lab 1차 검증 후보로 우선순위가 높습니다.

따라서 Paper 9의 안전한 claim은 “validated synthetic lethality”가 아니라 “lineage-silenced thyroid cancer에서 glutamine-axis vulnerability candidates를 prioritization했다”입니다. 후속 검증은 SLC1A5, GLUD1, GLS를 중심으로 glutamine withdrawal, inhibitor sensitivity, rescue assay를 설계하는 것이 적절합니다.
"""


def metadata_sanity() -> pd.DataFrame:
    terms = ["BPTES", "CB-839", "telaglenastat", "glutaminase", "SLC1A5 inhibitor", "V-9302", "GLUD1 inhibitor", "R162"]
    files = [
        RAW / "Repurposing_Public_24Q2_Treatment_Meta_Data.csv",
        RAW / "GDSC1Log2ViabilityCollapsedConditions.csv",
        RAW / "GDSC2Log2ViabilityCollapsedConditions.csv",
        RAW / "CTRPLog2ViabilityCollapsedConditions.csv",
    ]
    rows = []
    for path in files:
        df = pd.read_csv(path)
        text = df.astype(str).agg(" ".join, axis=1)
        for term in terms:
            hit = df[text.str.contains(re.escape(term), case=False, na=False)]
            example = "" if hit.empty else " ".join(map(str, hit.iloc[0].tolist()))
            rows.append(
                {
                    "metadata_file": path.name,
                    "search_term": term,
                    "n_matches": len(hit),
                    "example_match": example,
                    "testability": "testable_if_matrix_column_present" if len(hit) else "not_testable",
                }
            )
    out = REPORT / "paper9_drug_metadata_sanity_check.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    return pd.DataFrame(rows)


def main() -> None:
    ranking = pd.read_csv(S2 / "paper9_metabolic_candidate_ranking.tsv", sep="\t")
    assoc = pd.read_csv(S2 / "paper9_module_dependency_associations.tsv", sep="\t")
    artifact = pd.read_csv(S2 / "paper9_artifact_control_summary.tsv", sep="\t")
    drug = pd.read_csv(S2 / "paper9_drug_concordance_summary.tsv", sep="\t")
    meta = metadata_sanity()

    web_rank = ranking.copy()
    order = ["SLC1A5", "GLUD1", "GLS", "MYC_glutamine_addiction_module", "SLC7A5", "GOT2"]
    comparator = ranking[ranking["candidate"].str.contains("LDHA|OXPHOS|NAMPT", regex=True)].copy()
    web_rank["Rank"] = web_rank["candidate"].map({c: i + 1 for i, c in enumerate(order)})
    comparator["Rank"] = 7
    web_rank = pd.concat([web_rank[web_rank["candidate"].isin(order)], comparator], ignore_index=True)
    web_rank["Rank"] = web_rank["Rank"].fillna(99).astype(int)
    web_rank = web_rank.sort_values(["Rank", "candidate"])
    web_rank = web_rank.rename(
        columns={
            "candidate": "Candidate",
            "class": "Class",
            "dependency_evidence": "Dependency evidence",
            "expression_evidence": "Expression evidence",
            "drug_response_evidence": "Drug response evidence",
            "artifact_risk": "Artifact risk",
            "thyroid_specific_coverage": "Thyroid-specific coverage",
            "verdict": "Verdict",
            "claim_status": "Claim status",
        }
    )
    web_rank[
        [
            "Rank",
            "Candidate",
            "Class",
            "Dependency evidence",
            "Expression evidence",
            "Drug response evidence",
            "Artifact risk",
            "Thyroid-specific coverage",
            "Verdict",
            "Claim status",
        ]
    ].to_csv(S2 / "paper9_candidate_ranking_web.tsv", sep="\t", index=False)

    boundary = """# Paper 9 Claim Boundary

## Allowed
- Candidate metabolic vulnerability.
- Research-use prioritization.
- Wet-lab validation candidates.

## Weak
- Thyroid-specific inference, because strict thyroid CRISPR overlap is small.

## Forbidden
- Validated synthetic lethality.
- Clinical treatment recommendation.
- Patient selection biomarker.
"""
    (REPORT / "paper9_claim_boundary.md").write_text(boundary)

    wetlab = """# Paper 9 Wet-Lab Validation Plan

## Primary Comparison
- Compare lineage-high vs lineage-silenced thyroid models.

## Primary Candidates
- SLC1A5 inhibition or knockdown.
- GLUD1 knockdown or inhibition.
- GLS inhibitor plus glutamine withdrawal.

## Assay Logic
- Glutamine withdrawal.
- Inhibitor sensitivity.
- Glutamine rescue assay.
- MYC/proliferation artifact control.
- Optional isotope tracing.

## Claim Guard
Successful wet-lab signal would support experimental vulnerability validation, not clinical treatment recommendation.
"""
    (REPORT / "paper9_wetlab_validation_plan.md").write_text(wetlab)
    (REPORT / "paper9_professor_summary_kr.md").write_text("# Paper 9 Professor Summary\n\n" + KR_SUMMARY)

    summary = {
        "title": "Paper 9: Metabolic Vulnerability Prioritization Map of Lineage-Silenced Thyroid Cancer",
        "status": "GO as metabolic vulnerability prioritization map",
        "not_go": "Not GO as validated synthetic lethality",
        "top_wetlab_candidates": ["SLC1A5", "GLUD1", "GLS"],
        "directionality_sanity": {
            "dependency_strength": "positive dependency_strength is -CRISPR_gene_effect, so positive association means stronger dependency with lineage silencing",
            "drug_auc": "AUC matrices are direction-normalized as sensitivity = -AUC; higher raw AUC generally means less sensitivity/resistance, not sensitivity",
            "prism_lfc": "PRISM LFC would require compound match; lower/more negative LFC indicates stronger depletion/sensitivity; no GLS/SLC1A5/GLUD1 processed metadata match found",
            "bptes": "GDSC1 BPTES concordance is weak/non-supportive for GLS: lineage sensitivity r=-0.0153, p=0.694",
        },
        "candidate_snapshot": ranking[ranking["candidate"].isin(["SLC1A5", "GLUD1", "GLS"])].to_dict("records"),
        "drug_metadata_sanity_counts": meta.groupby("search_term")["n_matches"].sum().to_dict(),
        "claim_boundary": {
            "allowed": ["candidate metabolic vulnerability", "research-use prioritization", "wet-lab validation candidates"],
            "weak": ["thyroid-specific inference"],
            "forbidden": ["validated synthetic lethality", "clinical treatment recommendation", "patient selection biomarker"],
        },
    }
    (REPORT / "paper9_web_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"wrote": ["paper9_web_summary.json", "paper9_candidate_ranking_web.tsv", "paper9_claim_boundary.md", "paper9_wetlab_validation_plan.md", "paper9_professor_summary_kr.md"]}, indent=2))


if __name__ == "__main__":
    main()
