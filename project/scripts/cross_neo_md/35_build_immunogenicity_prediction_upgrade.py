#!/usr/bin/env python3
"""Build a CROSS-Neo-I immunogenicity prediction upgrade package.

This package converts the current prioritization/MD/wetlab-control stack into
an explicit immunogenicity-prediction model contract. It does not claim that
the current scores prove immunogenicity; it defines what data, labels, model
layers, and benchmarks are required to make that claim defensible.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
PKG = MD_OUT / "high_impact_decision_package"
LEAD = MD_OUT / "two_lead_impact_package"
PRECLIN = MD_OUT / "preclinical_validation_protocol"
OUT = MD_OUT / "immunogenicity_prediction_upgrade"
FIG = MD_OUT / "figures"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def nfloat(value: object, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default


def md_table(df: pd.DataFrame, cols: list[str], n: int = 30) -> str:
    if df.empty:
        return "_No rows._"
    keep = [c for c in cols if c in df.columns]
    return df[keep].head(n).to_markdown(index=False) if keep else "_No requested columns._"


def load_inputs() -> dict[str, pd.DataFrame]:
    data = {
        "top13": read_tsv(PKG / "no_false_positive_top13_candidates.tsv"),
        "leads": read_tsv(LEAD / "two_lead_external_evidence.tsv"),
        "preclin": read_tsv(PRECLIN / "preclinical_readiness_scorecard.tsv"),
        "evidence": read_tsv(PKG / "evidence_to_claim_matrix.tsv"),
        "claim_ladder": read_tsv(PKG / "claim_ladder.tsv"),
    }
    if data["top13"].empty or data["leads"].empty:
        raise FileNotFoundError("Missing high-impact/top-lead inputs")
    return data


def build_model_layers() -> pd.DataFrame:
    rows = [
        {
            "layer": "presentation",
            "model_or_feature_family": "BigMHC-like MHC-I presentation/immunogenicity transfer",
            "needed_inputs": "mutant peptide, HLA allele, optional expression/processing features",
            "target_signal": "peptide-HLA presentation and presentation-conditioned immunogenicity",
            "current_proxy_in_crossneo": "main_dl_score, bayes_mean, pMHC structural score",
            "what_it_adds": "separates presented candidates from peptides that are unlikely to reach pHLA surface",
            "source_or_precedent": "https://www.nature.com/articles/s42256-023-00694-6",
            "claim_boundary": "presentation is necessary but not sufficient for T-cell activation",
        },
        {
            "layer": "foreignness_and_self_delta",
            "model_or_feature_family": "PRIME/IMPROVE-like immunogenicity and TCR-recognition determinants",
            "needed_inputs": "mutant peptide, WT peptide, anchor positions, hydrophobic/aromatic core features",
            "target_signal": "mutant-over-self recognition potential",
            "current_proxy_in_crossneo": "WT/decoy readiness, anchor-preserved decoy, mutation-position context",
            "what_it_adds": "turns WT/decoy controls into a model feature rather than a post-hoc caveat",
            "source_or_precedent": "https://www.sciencedirect.com/science/article/pii/S2666379121000057; https://pubmed.ncbi.nlm.nih.gov/38633261/",
            "claim_boundary": "foreignness score does not prove a productive T-cell response",
        },
        {
            "layer": "paired_tcr_recognition",
            "model_or_feature_family": "PISTE/EPACT/TCR-HLA-antigen binding expert",
            "needed_inputs": "paired alpha/beta TCR where available, peptide, HLA",
            "target_signal": "TCR-pMHC recognition plausibility",
            "current_proxy_in_crossneo": "paired_tcr_evidence_count, tcr_augmented_score_mean",
            "what_it_adds": "rescues candidates like HMTEVVRHC that main pMHC-only DL may under-rank",
            "source_or_precedent": "https://www.nature.com/articles/s42256-024-00901-y; https://www.nature.com/articles/s42256-024-00913-8",
            "claim_boundary": "TCR expert should abstain when paired TCR evidence is absent",
        },
        {
            "layer": "structure_md_specificity",
            "model_or_feature_family": "OpenMM pMHC/TCR-pMHC stability and mutant-vs-WT/decoy delta",
            "needed_inputs": "mutant, WT, decoy pMHC/TCR-pMHC structures and trajectories",
            "target_signal": "stable pMHC and TCR-facing interface under controls",
            "current_proxy_in_crossneo": "MD_evidence_score, peptide RMSD, TCR-peptide contact tail",
            "what_it_adds": "explains why a candidate is structurally plausible or why a model false positive failed",
            "source_or_precedent": "local OpenMM audit outputs",
            "claim_boundary": "MD supports structural plausibility, not cytokine release",
        },
        {
            "layer": "tumor_context",
            "model_or_feature_family": "expression, clonality, HLA/B2M/TAP, immune-context features",
            "needed_inputs": "RNA expression, mutation clonality/VAF, HLA expression/LOH, antigen-processing context",
            "target_signal": "whether the pHLA target exists in the tumor and can be seen by T cells",
            "current_proxy_in_crossneo": "not yet integrated for this package",
            "what_it_adds": "moves from peptide immunogenicity to vaccine-response plausibility",
            "source_or_precedent": "NeoPrecis-style qualified immunogenicity and clonality-aware landscape concepts",
            "claim_boundary": "tumor context is required before vaccine-response claims",
        },
        {
            "layer": "calibration_and_abstention",
            "model_or_feature_family": "stacked calibrated ensemble with source-heldout and TCR-available abstention gates",
            "needed_inputs": "all layer outputs plus assay labels",
            "target_signal": "well-calibrated probability of mutant-specific T-cell activation in tested context",
            "current_proxy_in_crossneo": "threshold optimizer, uncertainty funnel, evidence-to-claim matrix",
            "what_it_adds": "prevents one impressive retrospective metric from becoming an overclaim",
            "source_or_precedent": "local threshold optimizer and claim ladder",
            "claim_boundary": "calibration must be measured on heldout/external labels",
        },
    ]
    return pd.DataFrame(rows)


def build_label_contract() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "label_name": "presentation_positive",
                "positive_definition": "mutant peptide is detected by MS elution or passes validated HLA binding/stability assay",
                "negative_definition": "matched assay shows no presentation/binding above background",
                "not_allowed_as_positive": "model-predicted binding only",
                "used_for": "presentation layer",
            },
            {
                "label_name": "recognition_positive",
                "positive_definition": "mutant pHLA binds paired TCR or T cells above WT, decoy, and irrelevant controls",
                "negative_definition": "no mutant-over-control binding under matched conditions",
                "not_allowed_as_positive": "peptide-HLA stability without TCR/cell binding",
                "used_for": "TCR recognition layer",
            },
            {
                "label_name": "activation_positive",
                "positive_definition": "ELISpot/ICS/cytokine response to mutant exceeds WT, decoy, and irrelevant controls with replicate support",
                "negative_definition": "mutant response absent, WT-like, decoy-like, or only at nonphysiologic dose",
                "not_allowed_as_positive": "multimer binding alone",
                "used_for": "primary immunogenicity model endpoint",
            },
            {
                "label_name": "functional_killing_positive",
                "positive_definition": "matched HLA/mutation target cells are killed while WT/decoy/mismatched controls are negative",
                "negative_definition": "no specific killing after presentation and activation gates pass",
                "not_allowed_as_positive": "activation without target-cell assay",
                "used_for": "functional validation endpoint",
            },
            {
                "label_name": "unknown_or_weak",
                "positive_definition": "not used as a positive class",
                "negative_definition": "not used as a negative class unless a matched assay proves failure",
                "not_allowed_as_positive": "literature mention, model score, public TCR match without assay context",
                "used_for": "semi-supervised or abstention bucket only",
            },
        ]
    )


def build_candidate_readiness(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    top13 = data["top13"].copy()
    leads = data["leads"].copy()
    preclin = data["preclin"].copy()
    rows = []
    for _, row in top13.sort_values("priority_order").head(13).iterrows():
        peptide = str(row.get("peptide", ""))
        lead = leads[leads["mutant_peptide"].astype(str).eq(peptide)]
        pre = preclin[preclin["mutant_peptide"].astype(str).eq(peptide)]
        has_wt = not lead.empty and bool(str(lead.iloc[0].get("wildtype_peptide", "")).strip())
        has_decoy = not lead.empty and bool(str(lead.iloc[0].get("anchor_preserved_decoy", "")).strip())
        presentation = max(nfloat(row.get("main_dl_score")), nfloat(row.get("bayes_mean")))
        tcr = min(nfloat(row.get("paired_tcr_evidence_count")) / 20.0, 1.0) * 0.45 + min(nfloat(row.get("tcr_augmented_score_mean")), 1.0) * 0.55
        md_score = min(max(nfloat(row.get("md_score")), nfloat(row.get("md_structural_score"))), 1.0)
        if str(row.get("md_label", "")).upper() == "MD_VERY_STRONG":
            md_score = max(md_score, 0.82)
        elif str(row.get("md_label", "")).upper() == "MD_MODERATE":
            md_score = max(md_score, 0.56)
        specificity = 1.0 if has_wt and has_decoy else (0.55 if has_wt or has_decoy else 0.0)
        context = 0.0
        assay_label = 0.0
        readiness = (
            0.20 * presentation
            + 0.20 * tcr
            + 0.15 * md_score
            + 0.20 * specificity
            + 0.10 * context
            + 0.15 * assay_label
        )
        model_tier = "CNI_FLAGSHIP_READY_FOR_LABEL" if readiness >= 0.55 else ("CNI_FEATURE_READY_LABEL_MISSING" if readiness >= 0.35 else "CNI_CURATION_FIRST")
        if not pre.empty:
            model_tier = pre.iloc[0].get("validation_tier", model_tier)
        rows.append(
            {
                "row_id": row.get("row_id", ""),
                "peptide": peptide,
                "hla_4digit": row.get("hla_4digit", ""),
                "source_dataset": row.get("source_dataset", ""),
                "current_label_binary": row.get("label_binary", ""),
                "presentation_proxy": presentation,
                "paired_tcr_proxy": tcr,
                "md_proxy": md_score,
                "wt_decoy_specificity_ready": specificity,
                "tumor_context_available": context,
                "activation_label_available": assay_label,
                "cross_neo_i_readiness_score": readiness,
                "cross_neo_i_tier": model_tier,
                "next_data_needed": "WT/decoy-controlled activation label; tumor expression/HLA context; external heldout labels",
                "claim_boundary": "readiness for immunogenicity prediction, not immunogenicity probability",
            }
        )
    return pd.DataFrame(rows)


def build_feature_contract(layers: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, layer in layers.iterrows():
        rows.append(
            {
                "model_layer": layer["layer"],
                "feature_block": layer["model_or_feature_family"],
                "must_have_before_training": layer["needed_inputs"],
                "primary_label": "activation_positive" if layer["layer"] in {"paired_tcr_recognition", "structure_md_specificity", "calibration_and_abstention"} else "presentation_positive",
                "allowed_missingness_strategy": "abstain or separate expert branch; do not impute as negative",
                "leakage_risk": "public epitope/TCR overlap and source leakage must be audited",
            }
        )
    rows.append(
        {
            "model_layer": "final_endpoint",
            "feature_block": "mutant-over-WT/decoy T-cell activation",
            "must_have_before_training": "matched assay labels with negative controls",
            "primary_label": "activation_positive",
            "allowed_missingness_strategy": "unknown bucket; never train unknown as negative",
            "leakage_risk": "assay batch/source split must be held out",
        }
    )
    return pd.DataFrame(rows)


def build_benchmark_contract() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "benchmark": "presentation_heldout",
                "split": "source-heldout eluted ligand / binding assay",
                "metric": "AUPRC, precision@k, calibration",
                "success_criterion": "improves over presentation-only baseline without public-overlap leakage",
                "claim_if_passed": "presentation prioritization generalizes",
            },
            {
                "benchmark": "activation_heldout",
                "split": "patient/source-heldout T-cell activation labels",
                "metric": "precision@k, recall at fixed FP, enrichment over prevalence",
                "success_criterion": "strict top-k remains enriched after WT/decoy-aware filtering",
                "claim_if_passed": "immunogenicity prioritization, not clinical utility",
            },
            {
                "benchmark": "paired_tcr_expert_ablation",
                "split": "paired-TCR-available rows only; epitope-heldout and TCR-heldout",
                "metric": "delta AUPRC and abstention safety",
                "success_criterion": "TCR expert helps when paired TCR exists and abstains otherwise",
                "claim_if_passed": "TCR-aware expert is useful in supported contexts",
            },
            {
                "benchmark": "structure_md_value_add",
                "split": "matched mutant/WT/decoy structural audit rows",
                "metric": "false-positive reduction and explanation rate",
                "success_criterion": "MD/contact features reduce WT/decoy or unstable-interface false positives",
                "claim_if_passed": "structure layer improves triage/explanation",
            },
            {
                "benchmark": "prospective_locked_panel",
                "split": "pre-registered candidate panel before assay results",
                "metric": "hit rate, FP burden, calibration drift",
                "success_criterion": "pre-registered top candidates exceed baseline hit rate under assay gates",
                "claim_if_passed": "prospective decision-support evidence",
            },
        ]
    )


def make_figures(
    layers: pd.DataFrame,
    label_contract: pd.DataFrame,
    readiness: pd.DataFrame,
    feature_contract: pd.DataFrame,
    benchmark: pd.DataFrame,
) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    fig, ax = plt.subplots(figsize=(12.2, 4.8))
    ax.axis("off")
    blocks = [
        ("Presentation", "BigMHC-like\npHLA"),
        ("Foreignness", "PRIME/IMPROVE\nWT delta"),
        ("TCR", "PISTE/EPACT\nexpert"),
        ("Structure", "OpenMM\nMD delta"),
        ("Context", "RNA/HLA\nclonality"),
        ("Calibration", "heldout\nassay labels"),
    ]
    xs = np.linspace(0.08, 0.92, len(blocks))
    for i, (title, body) in enumerate(blocks):
        ax.add_patch(plt.Rectangle((xs[i] - 0.065, 0.50), 0.13, 0.22, fc="#10243a", ec="#edf5ff", lw=1.2))
        ax.text(xs[i], 0.655, title, ha="center", va="center", color="#f2c46d", weight="bold", fontsize=9.5)
        ax.text(xs[i], 0.555, body, ha="center", va="center", color="white", fontsize=8.5)
        if i < len(blocks) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.075, 0.61), xytext=(xs[i] + 0.075, 0.61), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(0.5, 0.88, "CROSS-Neo-I immunogenicity prediction architecture", ha="center", fontsize=15, weight="bold")
    ax.text(0.5, 0.24, "Endpoint: mutant-specific T-cell activation above WT/decoy controls. Current package is model-readiness, not probability.", ha="center", fontsize=10, color="#9b1c31")
    fig.savefig(FIG / "fig_md51_cross_neo_i_architecture.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md51_cross_neo_i_architecture.pdf", bbox_inches="tight")
    plt.close(fig)

    matrix = np.array(
        [
            [1, 0, 0, 0, 0],
            [1, 1, 0, 0, 0],
            [0, 1, 1, 0, 0],
            [0, 1, 0.8, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 1, 1, 1, 1],
        ]
    )
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    im = ax.imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_yticks(range(len(layers)), layers["layer"].str.replace("_", "\n"))
    ax.set_xticks(range(len(label_contract)), label_contract["label_name"].str.replace("_", "\n"), rotation=25, ha="right")
    ax.set_title("Label-to-feature contract for immunogenicity prediction")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md52_label_feature_contract.png", dpi=220)
    fig.savefig(FIG / "fig_md52_label_feature_contract.pdf")
    plt.close(fig)

    show = readiness.sort_values("cross_neo_i_readiness_score", ascending=False).head(13)
    fig, ax = plt.subplots(figsize=(10.5, 6.0))
    y = np.arange(len(show))
    ax.barh(y, show["cross_neo_i_readiness_score"], color="#4e79a7")
    ax.set_yticks(y, show["peptide"] + "\n" + show["hla_4digit"].astype(str), fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("model-readiness score")
    ax.set_title("CROSS-Neo-I readiness: features present before activation labels")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md53_cross_neo_i_readiness_ranking.png", dpi=220)
    fig.savefig(FIG / "fig_md53_cross_neo_i_readiness_ranking.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    ax.axis("off")
    xs = np.linspace(0.08, 0.92, len(benchmark))
    for i, (_, row) in enumerate(benchmark.iterrows()):
        ax.add_patch(plt.Circle((xs[i], 0.60), 0.058, fc="#17324f", ec="#edf5ff", lw=1.2))
        ax.text(xs[i], 0.60, str(i + 1), ha="center", va="center", color="#f2c46d", weight="bold")
        ax.text(xs[i], 0.34, row["benchmark"].replace("_", "\n"), ha="center", va="center", fontsize=8)
        if i < len(benchmark) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.075, 0.60), xytext=(xs[i] + 0.075, 0.60), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(0.5, 0.86, "Benchmark ladder required before immunogenicity claims", ha="center", fontsize=15, weight="bold")
    fig.savefig(FIG / "fig_md54_immunogenicity_benchmark_ladder.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md54_immunogenicity_benchmark_ladder.pdf", bbox_inches="tight")
    plt.close(fig)

    risk = np.array([[1, 0.8, 0.6], [0.9, 1, 0.7], [0.6, 0.8, 1], [0.5, 0.7, 0.9], [0.4, 0.6, 1]])
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    im = ax.imshow(risk, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_yticks(range(len(benchmark)), benchmark["benchmark"].str.replace("_", "\n"), fontsize=8)
    ax.set_xticks(range(3), ["leakage\nrisk", "source\nshift", "claim\nrisk"])
    ax.set_title("Immunogenicity model risk-control board")
    for i in range(risk.shape[0]):
        for j in range(risk.shape[1]):
            ax.text(j, i, f"{risk[i, j]:.1f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md55_immunogenicity_risk_control_board.png", dpi=220)
    fig.savefig(FIG / "fig_md55_immunogenicity_risk_control_board.pdf")
    plt.close(fig)


def write_reports(
    layers: pd.DataFrame,
    label_contract: pd.DataFrame,
    readiness: pd.DataFrame,
    feature_contract: pd.DataFrame,
    benchmark: pd.DataFrame,
) -> None:
    top = readiness.sort_values("cross_neo_i_readiness_score", ascending=False).head(4)
    lines = [
        "# CROSS-Neo-I Immunogenicity Prediction Upgrade",
        "",
        "## Executive Verdict",
        "",
        "The next impact jump is to make CROSS-Neo-I: a calibrated immunogenicity-prediction layer whose endpoint is mutant-specific T-cell activation above WT/decoy controls. The current package should be framed as model-readiness and experimental design, not as a final immunogenicity probability.",
        "",
        "## Model Layers",
        "",
        md_table(layers, ["layer", "model_or_feature_family", "needed_inputs", "target_signal", "claim_boundary"], 10),
        "",
        "## Label Contract",
        "",
        md_table(label_contract, ["label_name", "positive_definition", "negative_definition", "not_allowed_as_positive", "used_for"], 10),
        "",
        "## Top Candidate Readiness",
        "",
        md_table(top, ["peptide", "hla_4digit", "presentation_proxy", "paired_tcr_proxy", "md_proxy", "wt_decoy_specificity_ready", "cross_neo_i_readiness_score", "cross_neo_i_tier"], 8),
        "",
        "## Feature Contract",
        "",
        md_table(feature_contract, ["model_layer", "feature_block", "must_have_before_training", "primary_label", "allowed_missingness_strategy", "leakage_risk"], 20),
        "",
        "## Benchmark Contract",
        "",
        md_table(benchmark, ["benchmark", "split", "metric", "success_criterion", "claim_if_passed"], 20),
        "",
        "## Claim Boundary",
        "",
        "- Allowed now: CROSS-Neo-I model contract, feature/label schema, benchmark plan, and readiness ranking.",
        "- Not allowed now: final calibrated immunogenicity probability, vaccine efficacy, or clinical response prediction.",
        "- Required next: WT/decoy-controlled activation labels and source-heldout benchmarks.",
    ]
    (OUT / "CROSS_NEO_I_IMMUNOGENICITY_UPGRADE_REPORT.md").write_text("\n".join(lines) + "\n")

    kr = [
        "# CROSS-Neo-I 면역원성 예측 업그레이드 요약",
        "",
        "임팩트를 더 올리려면 지금의 후보선별/MD 패키지를 **면역원성 예측 모델 계약서**로 바꿔야 합니다.",
        "",
        "- endpoint는 `mutant-specific T-cell activation > WT/decoy controls` 입니다.",
        "- BigMHC-like presentation, PRIME/IMPROVE-like foreignness, PISTE/EPACT-like TCR expert, OpenMM MD, tumor context, calibration을 분리합니다.",
        "- unknown을 negative로 학습하면 안 됩니다. assay-confirmed negative만 negative입니다.",
        "- HMTEVVRHC와 GADGVGKSAL은 CROSS-Neo-I의 flagship label-generation cases입니다.",
        "",
        "한 줄 결론: **이제 대박 포인트는 예측 점수 하나가 아니라, 면역원성 label을 만들고 calibrate할 수 있는 모델/실험 contract를 갖췄다는 점입니다.**",
    ]
    (OUT / "CROSS_NEO_I_ONE_PAGE_KR.md").write_text("\n".join(kr) + "\n")

    kakao = (
        "CROSS-Neo-I로 임팩트 올리는 방향은 명확합니다. 지금 후보선별/MD 패키지에 BigMHC-like presentation, "
        "PRIME/IMPROVE-like foreignness, PISTE/EPACT-like TCR expert, OpenMM mutant-vs-WT/decoy MD, tumor context, "
        "source-heldout calibration을 붙여서 endpoint를 'mutant-specific T-cell activation > WT/decoy controls'로 정의합니다. "
        "현재는 면역원성 확률을 확정한 게 아니라, HMTEVVRHC/GADGVGKSAL을 flagship label-generation case로 삼아 면역원성 예측 모델을 만들 수 있는 "
        "label/feature/benchmark contract가 완성된 상태입니다."
    )
    (OUT / "KAKAO_CROSS_NEO_I_IMMUNOGENICITY_KR.md").write_text(kakao + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    layers = build_model_layers()
    label_contract = build_label_contract()
    readiness = build_candidate_readiness(data)
    feature_contract = build_feature_contract(layers, label_contract)
    benchmark = build_benchmark_contract()

    layers.to_csv(OUT / "immunogenicity_model_layers.tsv", sep="\t", index=False)
    label_contract.to_csv(OUT / "immunogenicity_label_contract.tsv", sep="\t", index=False)
    readiness.to_csv(OUT / "cross_neo_i_candidate_readiness.tsv", sep="\t", index=False)
    feature_contract.to_csv(OUT / "immunogenicity_feature_contract.tsv", sep="\t", index=False)
    benchmark.to_csv(OUT / "immunogenicity_benchmark_contract.tsv", sep="\t", index=False)

    make_figures(layers, label_contract, readiness, feature_contract, benchmark)
    write_reports(layers, label_contract, readiness, feature_contract, benchmark)

    summary = {
        "n_model_layers": int(len(layers)),
        "n_label_contracts": int(len(label_contract)),
        "n_candidate_readiness_rows": int(len(readiness)),
        "n_benchmark_axes": int(len(benchmark)),
        "figures": [
            "fig_md51_cross_neo_i_architecture",
            "fig_md52_label_feature_contract",
            "fig_md53_cross_neo_i_readiness_ranking",
            "fig_md54_immunogenicity_benchmark_ladder",
            "fig_md55_immunogenicity_risk_control_board",
        ],
        "boundary": "immunogenicity-prediction model contract, not final immunogenicity probability",
    }
    (OUT / "immunogenicity_prediction_upgrade_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
