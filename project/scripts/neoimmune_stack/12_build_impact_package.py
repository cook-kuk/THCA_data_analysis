#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import REPO_ROOT, ensure_run_dir, safe_read_table, write_md, write_tsv


def try_cv2():
    try:
        import cv2

        return cv2
    except Exception:
        return None


def canvas(w=1800, h=1100, color=(8, 11, 18)):
    img = np.zeros((h, w, 3), dtype=np.uint8)
    img[:] = color
    return img


def put(img, text, x, y, size=0.75, color=(238, 245, 255), thick=2):
    cv2 = try_cv2()
    if cv2 is None:
        return
    cv2.putText(img, text, (x, y), cv2.FONT_HERSHEY_SIMPLEX, size, color, thick, cv2.LINE_AA)


def rect(img, x1, y1, x2, y2, color=(17, 24, 39), border=(38, 51, 72), thick=2):
    cv2 = try_cv2()
    if cv2 is None:
        return
    cv2.rectangle(img, (x1, y1), (x2, y2), color, -1)
    cv2.rectangle(img, (x1, y1), (x2, y2), border, thick)


def save(img, path: Path):
    cv2 = try_cv2()
    path.parent.mkdir(parents=True, exist_ok=True)
    if cv2 is not None:
        cv2.imwrite(str(path), img)


def draw_system_map(outdir: Path):
    cv2 = try_cv2()
    if cv2 is None:
        return
    img = canvas()
    put(img, "NeoImmune-Stack / CLEAN-Neo++", 60, 80, 1.25, (38, 217, 255), 3)
    put(img, "Leakage-aware patient-context ranking layer, not a new foundation model", 60, 130, 0.75, (190, 205, 230), 2)
    nodes = [
        ("Candidate universe", "CLEAN-Neobench\ncross_neo_v0\nneoantigen hub", 70, 240),
        ("Local algorithms", "Structure_LR\nWave8_TCR\nESM2_Bayesian\nQuantum kernels", 410, 240),
        ("Frozen public tools", "NetMHCpan\nMHCflurry\nBigMHC\nPRIME", 750, 240),
        ("Leakage audit", "peptide/pHLA\nstudy/patient/HLA\npublic-train risk", 1090, 240),
        ("Patient top-N", "Recall@20\nhit rate\ndisagreement audit", 1430, 240),
    ]
    for title, body, x, y in nodes:
        rect(img, x, y, x + 280, y + 250)
        put(img, title, x + 18, y + 42, 0.66, (246, 199, 104), 2)
        for i, line in enumerate(body.split("\n")):
            put(img, line, x + 18, y + 95 + i * 38, 0.58, (230, 238, 250), 2)
        if x < 1430:
            cv2.arrowedLine(img, (x + 286, y + 125), (x + 330, y + 125), (38, 217, 255), 4, tipLength=0.25)
    rect(img, 90, 640, 840, 890, (13, 19, 32), (94, 230, 168), 3)
    put(img, "Clean Science Track", 125, 700, 0.9, (94, 230, 168), 3)
    put(img, "Allowed: local peptide/WT/HLA/structure/TCR/ESM2/quantum features", 125, 755, 0.62)
    put(img, "Forbidden: NetMHCpan / MHCflurry / BigMHC / PRIME as training features", 125, 810, 0.62, (251, 113, 133), 2)
    rect(img, 960, 640, 1710, 890, (13, 19, 32), (167, 139, 250), 3)
    put(img, "Production Stack Track", 995, 700, 0.9, (167, 139, 250), 3)
    put(img, "Allowed: frozen public predictors + local branches + patient context", 995, 755, 0.62)
    put(img, "Claim: practical prioritization, not a clean novel model", 995, 810, 0.62, (246, 199, 104), 2)
    save(img, outdir / "figures" / "neoimmune_stack_system_map.png")


def draw_bar(df: pd.DataFrame, path: Path, title: str, metric: str = "AUPRC", n: int = 10):
    cv2 = try_cv2()
    if cv2 is None or df.empty or metric not in df.columns:
        return
    d = df.head(n).copy()
    img = canvas(1800, 1100)
    put(img, title, 60, 80, 1.15, (38, 217, 255), 3)
    put(img, "Strict no-existing-overlap view; use as ranking evidence, not clinical efficacy", 60, 130, 0.65, (190, 205, 230), 2)
    vals = pd.to_numeric(d[metric], errors="coerce").fillna(0).clip(0, 1).tolist()
    maxv = max(vals + [1e-6])
    x0, y0, w, gap = 530, 220, 1040, 62
    for i, (_, r) in enumerate(d.iterrows()):
        y = y0 + i * 82
        name = str(r.get("model", ""))[:42]
        val = float(pd.to_numeric(r.get(metric), errors="coerce") or 0)
        put(img, name, 60, y + 35, 0.55, (238, 245, 255), 2)
        rect(img, x0, y, x0 + w, y + 40, (31, 41, 55), (31, 41, 55), 1)
        fill = int(w * val / max(maxv, 1e-6))
        cv2.rectangle(img, (x0, y), (x0 + fill, y + 40), (38, 217, 255) if i < 3 else (167, 139, 250), -1)
        put(img, f"{metric} {val:.3f}", x0 + w + 25, y + 30, 0.55, (246, 199, 104), 2)
        put(img, f"n={fmt_int(r.get('n'))}, pos={fmt_int(r.get('positives'))}", x0 + w + 25, y + 62, 0.45, (160, 176, 200), 1)
    save(img, path)


def fmt_int(x):
    try:
        return f"{int(float(x)):,}"
    except Exception:
        return "NA"


def draw_claim_boundary(outdir: Path):
    cv2 = try_cv2()
    if cv2 is None:
        return
    img = canvas(1800, 1100)
    put(img, "Claim Ladder: what is strong now vs what needs wet-lab/patient data", 60, 80, 1.0, (38, 217, 255), 3)
    cols = [
        ("Strong now", (94, 230, 168), ["repo integration", "model registry", "canonical schema", "leakage audit", "strict score tables"]),
        ("Promising but bounded", (246, 199, 104), ["Wave8/TCR branch", "quantum branch", "production stack", "top-N workflow", "LLM rationale audit"]),
        ("Not claimable yet", (251, 113, 133), ["clinical efficacy", "true patient endpoint", "presentation without MS", "immunogenicity without T-cell assay", "prospective utility"]),
    ]
    for idx, (head, color, items) in enumerate(cols):
        x = 90 + idx * 570
        rect(img, x, 200, x + 500, 900, (13, 19, 32), color, 3)
        put(img, head, x + 30, 260, 0.85, color, 3)
        for j, item in enumerate(items):
            y = 350 + j * 95
            cv2.circle(img, (x + 44, y - 10), 12, color, -1)
            put(img, item, x + 75, y, 0.67, (238, 245, 255), 2)
    save(img, outdir / "figures" / "claim_ladder.png")


def draw_data_gap(outdir: Path):
    cv2 = try_cv2()
    if cv2 is None:
        return
    img = canvas(1800, 1100)
    put(img, "The single blocker: patient-context data", 60, 80, 1.1, (38, 217, 255), 3)
    put(img, "The system is built; the paper-grade endpoint needs true per-patient candidate sets.", 60, 130, 0.72, (190, 205, 230), 2)
    items = [
        ("Patient ID", "required for top-N hit rate"),
        ("Tumor expression", "exclude silent candidates"),
        ("VAF / clonality", "prioritize clonal targets"),
        ("HLA typing", "patient-specific restriction"),
        ("HLA LOH / APM", "presentation failure risk"),
        ("MS / T-cell assay", "claim presentation/immunogenicity"),
    ]
    for i, (a, b) in enumerate(items):
        x = 120 + (i % 3) * 560
        y = 260 + (i // 3) * 300
        rect(img, x, y, x + 460, y + 210, (13, 19, 32), (246, 199, 104), 3)
        put(img, a, x + 30, y + 70, 0.85, (246, 199, 104), 3)
        put(img, b, x + 30, y + 125, 0.62, (238, 245, 255), 2)
    rect(img, 260, 900, 1540, 1010, (32, 13, 22), (251, 113, 133), 3)
    put(img, "Without these: call it a ranking infrastructure scaffold, not a validated patient vaccine model.", 300, 965, 0.72, (251, 113, 133), 3)
    save(img, outdir / "figures" / "patient_data_gap.png")


def build_reports(outdir: Path, clean: pd.DataFrame, prod: pd.DataFrame, canon: pd.DataFrame):
    best_clean = clean.head(1).to_dict("records")
    best_prod = prod.head(1).to_dict("records")
    if "patient_id" in canon:
        real_mask = canon["patient_id"].fillna("unknown_patient").astype(str).ne("unknown_patient")
        real_patient_rows = int(real_mask.sum())
        real_patient_unique = int(canon.loc[real_mask, "patient_id"].astype(str).nunique())
    else:
        real_patient_rows = 0
        real_patient_unique = 0
    memo = f"""# CLEAN-Neo++ / NeoImmune-Stack Impact Memo KR

## 한 줄 결론
새 foundation model이 아니라, **기존 내 알고리즘 + 공개 predictor + leakage audit + patient top-N handoff**를 묶은 운영 레이어로 가는 것이 가장 임팩트가 큽니다.

## 대박 포인트
- 내 알고리즘이 중심입니다: Structure_LR, Wave8/TCR-self-similarity, ESM2_Bayesian, W7A/W7B, quantum kernel 계열을 clean science track에 배치했습니다.
- 공개 predictor는 보조입니다: NetMHCpan/MHCflurry/BigMHC/PRIME은 production track의 frozen feature/comparator입니다.
- reviewer 방어력이 생겼습니다: clean vs production을 분리했고, strict no-existing-overlap 리더보드를 별도로 만들었습니다.
- 사업화 언어가 생겼습니다: NeoImmune-Stack은 병원/환자별 후보 prioritization operating layer입니다.

## 냉정한 약점
- 현재 public 통합 테이블은 real patient ID가 {real_patient_unique:,}명 / {real_patient_rows:,} rows라서 retrospective public patient-level ranking scaffold는 가능해졌습니다.
- 하지만 hospital-grade prospective endpoint, clinical utility, assay-confirmed presentation/immunogenicity는 아직 주장하면 안 됩니다.
- MS immunopeptidomics 없이 presentation claim을 하면 바로 공격받습니다.
- T-cell assay label 없이 immunogenicity claim을 하면 바로 공격받습니다.
- BAR-Neo 쪽 production 숫자가 강해 보여도 clean novel model claim으로 쓰면 안 됩니다.

## 지금 가장 강한 논문 문장
“CLEAN-Neo++ provides a leakage-aware strategy for integrating peptide-HLA presentation, TCR-visible immunogenicity, structure proxies, and quantum-kernel features into patient-level neoantigen candidate ranking.”

## 바로 필요한 병원 데이터
patient_id, sample_id, tumor type, HLA typing, mutation peptide, WT peptide, expression TPM, VAF, clonality, HLA LOH, B2M/APM, candidate source, assay labels, and top-N validation outcome.
"""
    write_md(memo, outdir / "reports" / "CLEAN_NEO_PLUS_PLUS_IMPACT_MEMO_KR.md")

    paper = """# High-Impact Paper Frame KR

## Title
CLEAN-Neo++: Leakage-aware integration of TCR-visible immunogenicity, structure proxies, and public presentation predictors for patient-level cancer vaccine candidate ranking

## Figure Flow
1. Problem: binding/presentation alone is not immunogenicity.
2. System: clean local algorithm track vs production frozen-predictor stack.
3. Leakage: exact peptide, peptide-HLA, study, patient, HLA allele, public-tool training risk.
4. Clean result: local algorithms under strict no-existing-overlap.
5. Production result: practical stack compared with BigMHC/MHCflurry/PRIME/NetMHCpan artifacts.
6. Failure audit: false positives and rescued positives.
7. Patient handoff: top-20 report format and metadata requirements.

## Reviewer-Proof Main Claim
We do not claim clinical efficacy. We show that candidate ranking needs a leakage-aware, patient-context-aware operating layer and that local immunogenicity branches can be evaluated separately from public presentation predictors.

## Must Not Say
- “validated vaccine candidates”
- “predicts clinical response”
- “presentation confirmed” without MS
- “immunogenicity confirmed” without T-cell assay
- “LLM predicts vaccine efficacy”
"""
    write_md(paper, outdir / "reports" / "paper_frame_high_impact_kr.md")

    llm = """# Latest LLM Copilot Blueprint

## Role
Use a current frontier LLM as a **rationale auditor and report copilot**, not as a scoring feature.

## Allowed LLM jobs
- Summarize candidate evidence into wet-lab-ready rationales.
- Detect contradictions: high score but low expression, HLA LOH, high WT similarity, weak presentation, leakage risk.
- Generate reviewer attack/defense tables.
- Draft patient-specific top-N handoff notes from structured TSV only.
- Convert model disagreement into experimental questions.

## Forbidden LLM jobs
- Directly assign immunogenicity labels.
- Override leakage flags.
- Convert public predictor scores into clean-track features.
- Claim clinical efficacy.

## Runtime contract
Set `NEOIMMUNE_LLM_MODEL` and provider keys externally. Store prompts and outputs with candidate IDs, input hashes, and model name for auditability.
"""
    write_md(llm, outdir / "reports" / "latest_llm_copilot_blueprint.md")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    clean = safe_read_table(outdir / "metrics" / "leaderboard_clean_track_strict_no_overlap.tsv")
    prod = safe_read_table(outdir / "metrics" / "leaderboard_production_track_strict_no_overlap.tsv")
    draw_system_map(outdir)
    draw_bar(clean, outdir / "figures" / "strict_clean_leaderboard.png", "Strict Clean Science Track")
    draw_bar(prod, outdir / "figures" / "strict_production_leaderboard.png", "Strict Production Stack Track")
    draw_claim_boundary(outdir)
    draw_data_gap(outdir)
    build_reports(outdir, clean, prod, canon)
    ladder = pd.DataFrame(
        [
            ["strong_now", "integration architecture", "scripts 00-12 and generated tables", "claim as infrastructure"],
            ["strong_now", "strict no-overlap leaderboard", "leaderboard_*_strict_no_overlap.tsv", "claim as benchmark view, not prospective validation"],
            ["promising", "Wave8/TCR branch", "clean strict leaderboard", "needs source/patient-heldout confirmation"],
            ["promising", "quantum kernel branch", "clean strict leaderboard if present", "must survive source-heldout"],
            ["not_yet", "patient-level endpoint", "0 real patient rows currently", "requires hospital/patient data"],
            ["not_yet", "presentation claim", "no MS evidence in current package", "requires immunopeptidomics"],
            ["not_yet", "immunogenicity claim", "public labels only", "requires T-cell assay"],
        ],
        columns=["tier", "claim", "support", "boundary"],
    )
    write_tsv(ladder, outdir / "reports" / "claim_ladder.tsv")
    print(outdir / "reports" / "CLEAN_NEO_PLUS_PLUS_IMPACT_MEMO_KR.md")


if __name__ == "__main__":
    main()
