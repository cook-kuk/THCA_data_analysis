#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import ensure_run_dir, safe_read_table, write_md, write_tsv


def cv2_or_none():
    try:
        import cv2

        return cv2
    except Exception:
        return None


def put(img, text, x, y, size=0.7, color=(238, 245, 255), thick=2):
    cv2 = cv2_or_none()
    if cv2 is not None:
        cv2.putText(img, text, (x, y), cv2.FONT_HERSHEY_SIMPLEX, size, color, thick, cv2.LINE_AA)


def box(img, x1, y1, x2, y2, fill=(13, 19, 32), line=(38, 51, 72), thick=2):
    cv2 = cv2_or_none()
    if cv2 is not None:
        cv2.rectangle(img, (x1, y1), (x2, y2), fill, -1)
        cv2.rectangle(img, (x1, y1), (x2, y2), line, thick)


def save(img, path: Path):
    cv2 = cv2_or_none()
    path.parent.mkdir(parents=True, exist_ok=True)
    if cv2 is not None:
        cv2.imwrite(str(path), img)


def base(w=1800, h=1100):
    img = np.zeros((h, w, 3), dtype=np.uint8)
    img[:] = (8, 11, 18)
    return img


def load_scores(outdir: Path):
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    local = safe_read_table(outdir / "local_scores" / "all_local_scores.tsv.gz")
    ext = safe_read_table(outdir / "external_scores" / "all_external_scores.tsv.gz")
    prod = safe_read_table(outdir / "production_track" / "production_stack_predictions.tsv")
    clean = safe_read_table(outdir / "clean_track" / "clean_science_predictions.tsv")
    canon["label_immunogenicity"] = pd.to_numeric(canon["label_immunogenicity"], errors="coerce")
    local["normalized_score"] = pd.to_numeric(local["normalized_score"], errors="coerce")
    ext["normalized_score"] = pd.to_numeric(ext["normalized_score"], errors="coerce")
    local_wide = local.pivot_table(index="candidate_id", columns="local_model_name", values="normalized_score", aggfunc="mean").reset_index()
    ext = ext[~ext["candidate_id"].astype(str).str.startswith("__")].copy()
    ext_wide = ext.pivot_table(index="candidate_id", columns="model_name", values="normalized_score", aggfunc="mean").reset_index()
    d = canon.merge(local_wide, on="candidate_id", how="left").merge(ext_wide, on="candidate_id", how="left", suffixes=("", "_external"))
    d = d.merge(prod[["candidate_id", "production_stack_score", "model_disagreement_score", "abstention_coverage"]], on="candidate_id", how="left")
    d = d.merge(clean[["candidate_id", "clean_science_score"]], on="candidate_id", how="left")
    return d, local_wide, ext_wide


def score_cols(df: pd.DataFrame, tokens: list[str]) -> list[str]:
    out = []
    for c in df.columns:
        s = c.lower()
        if any(t.lower() in s for t in tokens):
            out.append(c)
    return out


def rank_pct(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").rank(pct=True, ascending=False)


def build_rescue_tables(outdir: Path, d: pd.DataFrame):
    local_tcr = score_cols(d, ["Wave8", "TCR", "SelfSim"])
    local_quantum = score_cols(d, ["quantum", "W7A_QK", "GP_quantum", "VQC"])
    local_structure = score_cols(d, ["Structure_LR"])
    public = score_cols(d, ["MHCflurry", "BigMHC", "PRIME", "NetMHCpan"])
    for name, cols in {
        "local_tcr_best": local_tcr,
        "local_quantum_best": local_quantum,
        "local_structure_best": local_structure,
        "public_best": public,
    }.items():
        d[name] = d[cols].max(axis=1, skipna=True) if cols else np.nan
    d["public_rank_pct"] = rank_pct(d["public_best"])
    d["tcr_rank_pct"] = rank_pct(d["local_tcr_best"])
    d["quantum_rank_pct"] = rank_pct(d["local_quantum_best"])
    d["structure_rank_pct"] = rank_pct(d["local_structure_best"])
    d["local_rescue_score"] = d[["local_tcr_best", "local_quantum_best", "local_structure_best", "clean_science_score"]].max(axis=1, skipna=True)
    d["local_rescue_delta_vs_public"] = d["local_rescue_score"] - d["public_best"]
    labeled = d[d["label_immunogenicity"].notna()].copy()
    positives = labeled[labeled["label_immunogenicity"] == 1].copy()
    rescued = positives[(positives["local_rescue_score"] >= 0.70) & ((positives["public_best"].isna()) | (positives["public_best"] < 0.50))].copy()
    public_missed = positives[((positives["public_best"].isna()) | (positives["public_best"] < 0.50))].copy()
    binding_only_traps = labeled[(labeled["label_immunogenicity"] == 0) & (labeled["public_best"] >= 0.75) & (labeled["local_rescue_score"] < 0.55)].copy()
    controversial = labeled[(pd.to_numeric(labeled.get("model_disagreement_score"), errors="coerce") >= 0.25)].copy()
    keep = [
        "candidate_id",
        "dataset_source",
        "peptide_mut",
        "peptide_wt",
        "hla_allele",
        "label_immunogenicity",
        "public_best",
        "local_tcr_best",
        "local_quantum_best",
        "local_structure_best",
        "clean_science_score",
        "production_stack_score",
        "local_rescue_delta_vs_public",
        "model_disagreement_score",
        "source_study",
    ]
    for frame, name in [
        (rescued, "local_rescued_public_missed_positives.tsv"),
        (public_missed, "public_missed_positive_pool.tsv"),
        (binding_only_traps, "binding_only_false_positive_traps.tsv"),
        (controversial, "high_disagreement_candidates.tsv"),
    ]:
        for c in keep:
            if c not in frame.columns:
                frame[c] = np.nan
        write_tsv(frame[keep].sort_values("local_rescue_delta_vs_public", ascending=False), outdir / "predictions" / name)
    summary = pd.DataFrame(
        [
            ["labeled_candidates", len(labeled)],
            ["labeled_positives", len(positives)],
            ["public_missed_positives_public_best_lt_0.50", len(public_missed)],
            ["local_rescued_public_missed_positives", len(rescued)],
            ["binding_only_false_positive_traps", len(binding_only_traps)],
            ["high_disagreement_labeled_candidates", len(controversial)],
        ],
        columns=["metric", "value"],
    )
    write_tsv(summary, outdir / "metrics" / "blockbuster_rescue_summary.tsv")
    return summary, rescued, binding_only_traps, d


def draw_rescue(summary: pd.DataFrame, outdir: Path):
    cv2 = cv2_or_none()
    if cv2 is None:
        return
    vals = dict(zip(summary["metric"], summary["value"]))
    img = base()
    put(img, "Blockbuster Claim Engine: local algorithms rescue what binding-only can miss", 60, 80, 0.95, (38, 217, 255), 3)
    put(img, "This figure is a hypothesis-generation audit, not proof of vaccine efficacy.", 60, 130, 0.62, (190, 205, 230), 2)
    cards = [
        ("Labeled candidates", vals.get("labeled_candidates", 0), (38, 217, 255)),
        ("Positive labels", vals.get("labeled_positives", 0), (94, 230, 168)),
        ("Public missed positives", vals.get("public_missed_positives_public_best_lt_0.50", 0), (246, 199, 104)),
        ("Local rescued positives", vals.get("local_rescued_public_missed_positives", 0), (167, 139, 250)),
        ("Binding-only traps", vals.get("binding_only_false_positive_traps", 0), (251, 113, 133)),
    ]
    for i, (lab, val, color) in enumerate(cards):
        x = 80 + i * 340
        box(img, x, 220, x + 290, 440, (13, 19, 32), color, 3)
        put(img, f"{int(val):,}", x + 30, 315, 1.35, color, 4)
        put(img, lab, x + 30, 375, 0.55, (238, 245, 255), 2)
    # Draw conceptual split.
    box(img, 120, 570, 760, 870, (13, 19, 32), (251, 113, 133), 3)
    put(img, "Binding/presentation-only", 160, 640, 0.85, (251, 113, 133), 3)
    put(img, "Can rank presented-looking peptides", 160, 705, 0.62)
    put(img, "But may miss TCR-visible immunogenicity", 160, 765, 0.62, (246, 199, 104), 2)
    box(img, 1040, 570, 1680, 870, (13, 19, 32), (94, 230, 168), 3)
    put(img, "CLEAN-Neo++ local core", 1080, 640, 0.85, (94, 230, 168), 3)
    put(img, "TCR/self-similarity + structure + ESM2 + quantum", 1080, 705, 0.62)
    put(img, "Rescue candidates become wet-lab questions", 1080, 765, 0.62, (246, 199, 104), 2)
    cv2.arrowedLine(img, (790, 720), (1010, 720), (38, 217, 255), 5, tipLength=0.25)
    save(img, outdir / "figures" / "blockbuster_rescue_engine.png")


def draw_wetlab_bridge(outdir: Path):
    cv2 = cv2_or_none()
    if cv2 is None:
        return
    img = base()
    put(img, "From model score to publishable evidence", 60, 80, 1.1, (38, 217, 255), 3)
    put(img, "The next experiment must convert ranked hypotheses into presentation/immunogenicity evidence.", 60, 130, 0.65, (190, 205, 230), 2)
    steps = [
        ("Top-20 per patient", "NeoImmune-Stack ranking\n+ disagreement audit"),
        ("Synthesis / assay queue", "exclude low expression\nHLA LOH/APM loss"),
        ("MS evidence", "presentation claim\nonly if detected"),
        ("T-cell assay", "immunogenicity claim\nonly if reactive"),
        ("Paper-grade endpoint", "patient Recall@20\nhit rate + AUPRC"),
    ]
    for i, (title, body) in enumerate(steps):
        x = 70 + i * 345
        box(img, x, 260, x + 285, 560, (13, 19, 32), (38, 217, 255) if i < 2 else (94, 230, 168), 3)
        put(img, title, x + 20, 330, 0.64, (246, 199, 104), 2)
        for j, line in enumerate(body.split("\n")):
            put(img, line, x + 20, 400 + 45 * j, 0.55)
        if i < len(steps) - 1:
            cv2.arrowedLine(img, (x + 292, 410), (x + 335, 410), (167, 139, 250), 4, tipLength=0.25)
    box(img, 190, 720, 1610, 920, (32, 13, 22), (251, 113, 133), 3)
    put(img, "Hard rule", 240, 790, 0.85, (251, 113, 133), 3)
    put(img, "No MS = no presentation-confirmed claim. No T-cell assay = no immunogenicity-confirmed claim.", 240, 850, 0.72, (238, 245, 255), 2)
    save(img, outdir / "figures" / "wetlab_bridge_to_publishable_evidence.png")


def draw_market_wedge(outdir: Path):
    cv2 = cv2_or_none()
    if cv2 is None:
        return
    img = base()
    put(img, "Why this can be big: not another predictor, an operating layer", 60, 80, 1.0, (38, 217, 255), 3)
    rows = [
        ("Current bottleneck", "too many candidate peptides; binding score alone is insufficient"),
        ("Differentiator", "local immunogenicity branches: TCR-visible, structure, ESM2, quantum"),
        ("Defensibility", "clean track excludes public predictor leakage; production track labels external features"),
        ("Clinical bridge", "hospital data contract converts rankings to per-patient top-N queues"),
        ("Business wedge", "vendor-agnostic layer across NetMHCpan/MHCflurry/BigMHC/PRIME + local algorithms"),
    ]
    for i, (a, b) in enumerate(rows):
        y = 210 + i * 150
        box(img, 110, y, 1690, y + 105, (13, 19, 32), (38, 51, 72), 2)
        put(img, a, 150, y + 65, 0.72, (246, 199, 104), 2)
        put(img, b, 610, y + 65, 0.62, (238, 245, 255), 2)
    save(img, outdir / "figures" / "business_market_wedge.png")


def build_blockbuster_docs(outdir: Path, summary: pd.DataFrame, rescued: pd.DataFrame, traps: pd.DataFrame):
    vals = dict(zip(summary["metric"], summary["value"]))
    memo = f"""# Blockbuster Upgrade Memo KR

## 대박으로 가는 핵심 문장
NeoImmune-Stack / CLEAN-Neo++는 “또 하나의 binding predictor”가 아니라, **patient-specific vaccine 후보를 clean-science track과 production-stack track으로 분리해 ranking, leakage audit, wet-lab handoff까지 끝내는 operating layer**다.

## 지금 새로 만든 임팩트
- Public predictor가 낮게 본 positive pool: {int(vals.get('public_missed_positives_public_best_lt_0.50', 0)):,}
- Local branches가 rescue한 public-missed positives: {int(vals.get('local_rescued_public_missed_positives', 0)):,}
- Binding-only false-positive trap 후보: {int(vals.get('binding_only_false_positive_traps', 0)):,}
- High-disagreement 후보: {int(vals.get('high_disagreement_labeled_candidates', 0)):,}

## 이 숫자의 의미
이것은 “효능 입증”이 아니라, **왜 binding-only가 부족하고 local immunogenicity branch가 필요한지**를 보여주는 analysis story다.

## 논문에서 세게 말할 수 있는 것
- Candidate ranking은 binding/presentation만으로 끝나지 않는다.
- TCR-visible/self-similarity, structure proxy, ESM2, quantum kernel branch를 분리 평가해야 한다.
- Public predictor는 production에서 유용하지만 clean novelty에는 쓰면 안 된다.
- Patient-level top-N은 올바른 endpoint지만 병원 데이터가 필요하다.

## 절대 하면 안 되는 말
- “validated vaccine candidate”
- “clinical efficacy prediction”
- “presentation confirmed” without MS
- “immunogenicity confirmed” without T-cell assay
"""
    write_md(memo, outdir / "reports" / "BLOCKBUSTER_UPGRADE_MEMO_KR.md")

    pitch = """# 10-Slide Blockbuster Pitch Outline KR

## 1. Problem
Cancer vaccine ranking fails when HLA binding is treated as immunogenicity.

## 2. Core Insight
The true bottleneck is patient-level top-N immunogenicity prioritization.

## 3. Our System
CLEAN-Neo++ clean local algorithm track + NeoImmune-Stack production operating layer.

## 4. Local Algorithm Moat
Structure_LR, Wave8_TCR_SelfSim_full, ESM2_Bayesian, W7A/W7B, quantum kernels.

## 5. Public Predictor Integration
NetMHCpan, MHCflurry, BigMHC, PRIME are frozen comparators/features, not clean novelty.

## 6. Leakage Defense
Exact peptide, pHLA, mutant-WT, source/study/patient/HLA, public-tool training risk.

## 7. Rescue Story
Local branches can rescue public-missed positive candidates and expose binding-only traps.

## 8. Wet-lab Bridge
Top-20 per patient, expression/VAF/APM/HLA-LOH gates, MS/T-cell assay validation.

## 9. Product Wedge
Vendor-agnostic AI operating layer for hospital/cohort vaccine candidate prioritization.

## 10. Ask
Give us real patient-level candidate tables and validation labels; we convert this scaffold into a publishable patient-level endpoint.
"""
    write_md(pitch, outdir / "reports" / "blockbuster_10_slide_pitch_outline_kr.md")

    wetlab = """# Minimal Wet-Lab Experiment That Makes This Paper Real KR

## Cohort
5-10 patients first, then 30-50 for paper-grade patient-level metrics.

## Per patient
- HLA class I typing
- somatic mutation table
- mutant and WT peptides
- expression TPM
- VAF/clonality
- HLA LOH
- B2M/TAP/APM context

## Candidate selection
- Top 20 by production stack
- Top 10 clean-local rescue candidates
- 5 binding-only high-score negative-risk controls
- 5 model-disagreement controls

## Assays
- MS immunopeptidomics if feasible for presentation evidence
- T-cell assay for immunogenicity evidence
- Track synthesis success and assay failure separately

## Endpoints
- Primary: patient-level hit rate@20 and Recall@20
- Secondary: strict no-overlap AUPRC and rescued-positive fraction

## Claim after success
Patient-context-aware ranking improves wet-lab candidate prioritization. Still not a clinical efficacy claim.
"""
    write_md(wetlab, outdir / "reports" / "minimal_wetlab_experiment_plan_kr.md")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    d, _, _ = load_scores(outdir)
    summary, rescued, traps, d2 = build_rescue_tables(outdir, d)
    draw_rescue(summary, outdir)
    draw_wetlab_bridge(outdir)
    draw_market_wedge(outdir)
    build_blockbuster_docs(outdir, summary, rescued, traps)
    print(outdir / "reports" / "BLOCKBUSTER_UPGRADE_MEMO_KR.md")


if __name__ == "__main__":
    main()
