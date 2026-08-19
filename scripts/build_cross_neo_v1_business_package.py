#!/usr/bin/env python3
"""Build a practical CROSS-Neo v1 business/operations package.

This is not a clinical claim package. It creates a usable candidate intake
schema, a conservative sequence-only fallback triage scorer, and an internal
demo ranking from existing v1 predictions.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from cross_neo_v1_common import V1, ensure_v1_dirs, load_master, seq_similarity


AA = "ACDEFGHIKLMNPQRSTVWY"
OUT = V1 / "business_package"


def clean_pep(x: object) -> str:
    return "".join(a for a in str(x or "").upper() if a in AA)


def anchor_hydrophobic(pep: str) -> float:
    p = clean_pep(pep)
    if not p:
        return 0.0
    idx = [len(p) - 1]
    if len(p) > 2:
        idx.append(1)
    return float(sum(p[i] in set("AILMFWVY") for i in idx if 0 <= i < len(p)) / max(1, len(idx)))


def manufacturability_penalty(pep: str) -> float:
    p = clean_pep(pep)
    if not p:
        return 0.25
    penalty = 0.0
    if len(p) < 8 or len(p) > 11:
        penalty += 0.18
    if "C" in p:
        penalty += 0.08
    if any(p[i : i + 4].count(p[i]) == 4 for i in range(max(0, len(p) - 3))):
        penalty += 0.10
    hydro_frac = sum(a in "AILMFWVY" for a in p) / len(p)
    if hydro_frac > 0.72:
        penalty += 0.12
    return float(min(0.4, penalty))


def wt_delta_score(mut: str, wt: str) -> float:
    m = clean_pep(mut)
    w = clean_pep(wt)
    if not m or not w:
        return 0.35
    return float(1 - seq_similarity(m, w))


def expression_score(x: object) -> float:
    try:
        v = float(x)
    except Exception:
        return 0.35
    return float(np.clip(np.log1p(max(v, 0.0)) / np.log1p(50.0), 0, 1))


def vaf_score(x: object) -> float:
    try:
        v = float(x)
    except Exception:
        return 0.35
    if v > 1:
        v = v / 100.0
    return float(np.clip(v / 0.5, 0, 1))


def source_window_self_risk(row: pd.Series) -> float:
    pep = clean_pep(row.get("peptide_mut", ""))
    wt = clean_pep(row.get("peptide_wt", ""))
    win15 = clean_pep(row.get("source_window_15aa", ""))
    win30 = clean_pep(row.get("source_window_30aa", ""))
    risks = []
    if wt:
        risks.append(seq_similarity(pep, wt))
    if win15:
        risks.append(seq_similarity(pep, win15))
    if win30:
        risks.append(seq_similarity(pep, win30))
    return float(max(risks) if risks else 0.0)


def business_score(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        pep = clean_pep(r.get("peptide_mut", ""))
        length_ok = float(8 <= len(pep) <= 11)
        anchor = anchor_hydrophobic(pep)
        delta = wt_delta_score(r.get("peptide_mut", ""), r.get("peptide_wt", ""))
        expr = expression_score(r.get("tumor_rna_tpm", np.nan))
        vaf = vaf_score(r.get("tumor_vaf", np.nan))
        self_risk = source_window_self_risk(r)
        manu_pen = manufacturability_penalty(pep)
        missing_wt = int(not clean_pep(r.get("peptide_wt", "")))
        missing_expr = int(pd.isna(r.get("tumor_rna_tpm", np.nan)) or str(r.get("tumor_rna_tpm", "")).strip() == "")
        missing_vaf = int(pd.isna(r.get("tumor_vaf", np.nan)) or str(r.get("tumor_vaf", "")).strip() == "")
        abstention = missing_wt + missing_expr + missing_vaf + int(not str(r.get("hla", "")).strip())
        score = (
            0.22 * length_ok
            + 0.18 * anchor
            + 0.22 * delta
            + 0.16 * expr
            + 0.12 * vaf
            + 0.10 * (1 - self_risk)
            - manu_pen
        )
        score = float(np.clip(score, 0, 1))
        tier = "review"
        if abstention >= 3:
            tier = "abstain_missing_context"
        elif score >= 0.68:
            tier = "priority"
        elif score < 0.42:
            tier = "deprioritize"
        rows.append(
            {
                "candidate_id": r.get("candidate_id", r.get("sample_id", "")),
                "patient_id": r.get("patient_id", ""),
                "peptide_mut": r.get("peptide_mut", ""),
                "peptide_wt": r.get("peptide_wt", ""),
                "hla": r.get("hla", ""),
                "business_triage_score": score,
                "business_tier": tier,
                "length_ok": length_ok,
                "anchor_hydrophobic_score": anchor,
                "mut_wt_delta_score": delta,
                "expression_score": expr,
                "vaf_score": vaf,
                "self_similarity_risk": self_risk,
                "manufacturability_penalty": manu_pen,
                "abstention_missing_context_count": abstention,
                "interpretation": "triage only; requires full CROSS-Neo feature generation and wet-lab review",
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(["business_tier", "business_triage_score"], ascending=[True, False])


def write_templates() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    template = pd.DataFrame(
        [
            {
                "candidate_id": "PAT001_CAND001",
                "patient_id": "PAT001",
                "peptide_mut": "SLLMWITQC",
                "peptide_wt": "SLLTWITQC",
                "hla": "HLA-A*02:01",
                "source_protein": "GENE1",
                "source_window_15aa": "",
                "source_window_30aa": "",
                "tumor_rna_tpm": 12.5,
                "tumor_vaf": 0.32,
                "clonality": "clonal",
                "variant_consequence": "missense",
                "notes": "example only",
            }
        ]
    )
    template.to_csv(OUT / "candidate_input_template.tsv", sep="\t", index=False)
    schema = pd.DataFrame(
        [
            ["candidate_id", "required", "unique row id"],
            ["patient_id", "recommended", "patient or sample id for per-patient top-k"],
            ["peptide_mut", "required", "mutant peptide sequence, class-I length preferred"],
            ["peptide_wt", "strongly_recommended", "matched WT peptide for counterfactual/self-risk"],
            ["hla", "required", "HLA allele such as HLA-A*02:01"],
            ["source_protein", "recommended", "gene or protein id"],
            ["source_window_15aa", "optional", "local source sequence window"],
            ["source_window_30aa", "optional", "larger source sequence window"],
            ["tumor_rna_tpm", "recommended", "tumor expression proxy"],
            ["tumor_vaf", "recommended", "variant allele fraction 0-1 or 0-100"],
            ["clonality", "optional", "clonal/subclonal if available"],
            ["variant_consequence", "optional", "variant annotation"],
        ],
        columns=["column", "requirement", "description"],
    )
    schema.to_csv(OUT / "candidate_input_schema.tsv", sep="\t", index=False)


def write_internal_demo() -> None:
    master = load_master().copy()
    pred_path = V1 / "v1_all_model_predictions.tsv"
    if pred_path.exists():
        pred = pd.read_csv(pred_path, sep="\t")
        use = pred[(pred["split_name"] == "hla_stratified_group_5fold") & (pred["model"] == "prespecified_equal_weight_C_QK_no_anchor")].copy()
        if use.empty:
            use = pred[pred["split_name"] == "hla_stratified_group_5fold"].copy()
        use = use.sort_values("score", ascending=False).head(60)
        demo = use.merge(master, on=["sample_id", "label"], how="left")
        cols = ["sample_id", "label", "score", "model", "peptide_mut", "peptide_wt", "hla", "hla_supertype", "study", "strict_set_flag"]
        demo[[c for c in cols if c in demo.columns]].to_csv(OUT / "internal_locked_demo_priority_queue.tsv", sep="\t", index=False)
    strict = master[master["strict_set_flag"].astype(bool)].copy()
    cand = strict.rename(columns={"sample_id": "candidate_id"})
    cand["tumor_rna_tpm"] = np.nan
    cand["tumor_vaf"] = np.nan
    scored = business_score(cand.head(89))
    scored.to_csv(OUT / "sequence_only_business_triage_demo.tsv", sep="\t", index=False)


def write_docs() -> None:
    files = sorted(str(p.relative_to(V1)) for p in OUT.glob("*"))
    model_card = [
        "# CROSS-Neo Practical Prioritization Package",
        "",
        "Purpose: prioritize neoantigen candidates for internal research review and wet-lab triage.",
        "",
        "This is not a clinical diagnostic, not external validation, and not a public MHC predictor wrapper.",
        "",
        "## What It Uses",
        "- Mutant peptide, WT peptide when available, HLA, source context, expression/VAF if supplied.",
        "- Optional full CROSS-Neo v1 predictions when generated by the locked pipeline.",
        "- No MHCflurry, NetMHCpan, BigMHC, PRIME, MixMHCpred, or NetMHCstabpan scores as model features.",
        "",
        "## Business Workflow",
        "1. Fill `candidate_input_template.tsv` for each patient candidate bag.",
        "2. Run the sequence-only triage as a fast intake sanity check.",
        "3. Run the full CROSS-Neo feature pipeline for candidates selected for deeper review.",
        "4. Report per-patient top-k with evidence flags, self-risk, missing-context abstention, and wet-lab review status.",
        "5. Treat source-heldout and public-overlap warnings as business risk controls.",
        "",
        "## Decision Boundary",
        "- `priority`: good intake score and enough context for review.",
        "- `review`: plausible but missing or mixed context.",
        "- `deprioritize`: weak sequence/context score.",
        "- `abstain_missing_context`: too much missing information for confident prioritization.",
        "",
        "## Files",
        "\n".join(f"- `{f}`" for f in files),
    ]
    (OUT / "BUSINESS_MODEL_CARD.md").write_text("\n".join(model_card) + "\n")
    checklist = pd.DataFrame(
        [
            ["patient_candidate_bag", "required", "rank within patient, not only global score"],
            ["wet_lab_readout", "required_for_product_claim", "ELISpot/TCR/functional assay needed for claims"],
            ["public_overlap_audit", "required_for_benchmarking", "unresolved overlap blocks clean comparator claims"],
            ["full_feature_generation", "required_for_final_rank", "sequence-only score is intake triage only"],
            ["abstention_policy", "required", "missing WT/HLA/expression/VAF must lower confidence"],
            ["regulatory_boundary", "required", "research-use prioritization only unless clinically validated"],
        ],
        columns=["item", "status", "business_reason"],
    )
    checklist.to_csv(OUT / "commercialization_readiness_checklist.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="optional candidate TSV to score with sequence-only triage")
    args = parser.parse_args()
    ensure_v1_dirs()
    write_templates()
    write_internal_demo()
    if args.input:
        df = pd.read_csv(args.input, sep="\t")
        business_score(df).to_csv(OUT / "candidate_sequence_only_triage_scores.tsv", sep="\t", index=False)
    write_docs()
    print(f"[v1-business] wrote {OUT}")


if __name__ == "__main__":
    main()
