#!/usr/bin/env python3
"""Build what-if assay feedback scenarios for CROSS-Neo.

The rows created here are simulated decision scenarios, not real assay labels.
They are kept separate from project/data/cross_neo_assay_feedback/assay_feedback.tsv
so they cannot contaminate the real closed-loop learner.
"""

from __future__ import annotations

import html
import json
import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
IMMUNO = MD_OUT / "immunogenicity_prediction_upgrade"
PRECLIN = MD_OUT / "preclinical_validation_protocol"
OUT = MD_OUT / "assay_feedback_scenarios"
FIG = MD_OUT / "figures"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_assay_feedback_scenarios.html"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
WEB_PAGE = WEB / "cross_neo_assay_feedback_scenarios.html"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def nfloat(value: object, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default


def esc(value: object) -> str:
    return html.escape("" if pd.isna(value) else str(value))


def fmt(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.{digits}f}"
    except Exception:
        return esc(value)


def table_html(df: pd.DataFrame, cols: list[str], n: int = 30) -> str:
    if df.empty:
        return "<p class='muted'>No rows.</p>"
    keep = [c for c in cols if c in df.columns]
    lines = ["<table><thead><tr>"]
    lines.extend(f"<th>{esc(c.replace('_', ' '))}</th>" for c in keep)
    lines.append("</tr></thead><tbody>")
    for _, row in df.head(n).iterrows():
        lines.append("<tr>")
        for c in keep:
            value = row.get(c, "")
            lines.append(f"<td>{fmt(value) if isinstance(value, (float, int)) else esc(value)}</td>")
        lines.append("</tr>")
    lines.append("</tbody></table>")
    return "\n".join(lines)


def md_table(df: pd.DataFrame, cols: list[str], n: int = 30) -> str:
    if df.empty:
        return "_No rows._"
    keep = [c for c in cols if c in df.columns]
    return df[keep].head(n).to_markdown(index=False) if keep else "_No requested columns._"


def safe_copy(src: Path, dst: Path) -> None:
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    except PermissionError:
        print(f"[assay-scenarios] skip permission-denied copy {dst}")


def scenario_defs() -> list[dict[str, object]]:
    return [
        {
            "scenario": "presentation_fail",
            "presentation": 0.0,
            "recognition": math.nan,
            "activation": math.nan,
            "killing": math.nan,
            "specificity": math.nan,
            "wt_risk": 0.0,
            "decoy_risk": 0.0,
            "meaning": "mutant pHLA does not pass binding/stability gate",
        },
        {
            "scenario": "presentation_only",
            "presentation": 1.0,
            "recognition": math.nan,
            "activation": math.nan,
            "killing": math.nan,
            "specificity": math.nan,
            "wt_risk": 0.0,
            "decoy_risk": 0.0,
            "meaning": "presentation supported, recognition still unknown",
        },
        {
            "scenario": "recognition_positive_activation_missing",
            "presentation": 1.0,
            "recognition": 1.0,
            "activation": math.nan,
            "killing": math.nan,
            "specificity": 1.0,
            "wt_risk": 0.0,
            "decoy_risk": 0.0,
            "meaning": "TCR/multimer recognition is mutant-specific; activation not tested",
        },
        {
            "scenario": "activation_specific_positive",
            "presentation": 1.0,
            "recognition": 1.0,
            "activation": 1.0,
            "killing": math.nan,
            "specificity": 1.0,
            "wt_risk": 0.0,
            "decoy_risk": 0.0,
            "meaning": "mutant activation exceeds WT/decoy controls",
        },
        {
            "scenario": "functional_best_case",
            "presentation": 1.0,
            "recognition": 1.0,
            "activation": 1.0,
            "killing": 1.0,
            "specificity": 1.0,
            "wt_risk": 0.0,
            "decoy_risk": 0.0,
            "meaning": "activation and killing are mutant-specific",
        },
        {
            "scenario": "recognition_positive_activation_negative",
            "presentation": 1.0,
            "recognition": 1.0,
            "activation": 0.0,
            "killing": math.nan,
            "specificity": 1.0,
            "wt_risk": 0.0,
            "decoy_risk": 0.0,
            "meaning": "binding/recognition without productive activation",
        },
        {
            "scenario": "wt_cross_reactive",
            "presentation": 1.0,
            "recognition": 1.0,
            "activation": 1.0,
            "killing": math.nan,
            "specificity": 0.0,
            "wt_risk": 1.0,
            "decoy_risk": 0.0,
            "meaning": "WT is positive; mutant-specific claim is blocked",
        },
        {
            "scenario": "decoy_positive",
            "presentation": 1.0,
            "recognition": 1.0,
            "activation": 1.0,
            "killing": math.nan,
            "specificity": 0.0,
            "wt_risk": 0.0,
            "decoy_risk": 1.0,
            "meaning": "decoy is positive; nonspecific or motif artifact risk",
        },
    ]


def posterior(prior: float, sc: dict[str, object]) -> tuple[float, float, str, str, str]:
    alpha = 1.0 + prior * 6.0
    beta = 1.0 + (1.0 - prior) * 6.0
    for key, weight in [
        ("presentation", 1.0),
        ("recognition", 1.5),
        ("activation", 3.0),
        ("killing", 4.0),
        ("specificity", 2.5),
    ]:
        value = nfloat(sc.get(key), math.nan)
        if not math.isnan(value):
            alpha += value * weight
            beta += (1.0 - value) * weight
    wt_risk = nfloat(sc.get("wt_risk"), math.nan)
    decoy_risk = nfloat(sc.get("decoy_risk"), math.nan)
    if not math.isnan(wt_risk):
        beta += wt_risk * 3.0
    if not math.isnan(decoy_risk):
        beta += decoy_risk * 3.0
    post = alpha / (alpha + beta)
    uncertainty = post * (1.0 - post)
    if wt_risk > 0:
        state = "SPECIFICITY_RISK_HOLD_WT_POSITIVE"
    elif decoy_risk > 0:
        state = "SPECIFICITY_RISK_HOLD_DECOY_POSITIVE"
    elif sc.get("killing") == 1.0 and sc.get("specificity") == 1.0:
        state = "FUNCTIONAL_SPECIFICITY_SUPPORTED"
    elif sc.get("activation") == 1.0 and sc.get("specificity") == 1.0:
        state = "ASSAY_SPECIFIC_IMMUNOGENICITY_SUPPORTED"
    elif sc.get("activation") == 0.0:
        state = "RECOGNITION_WITHOUT_ACTIVATION"
    elif sc.get("recognition") == 1.0 and sc.get("specificity") == 1.0:
        state = "RECOGNITION_PLAUSIBLE"
    elif sc.get("presentation") == 1.0:
        state = "PRESENTATION_SUPPORTED"
    elif sc.get("presentation") == 0.0:
        state = "PRESENTATION_FAILED"
    else:
        state = "ASSAY_INCONCLUSIVE"
    allowed = {
        "PRESENTATION_FAILED": "no biological claim; hold or repeat presentation",
        "PRESENTATION_SUPPORTED": "presentation plausibility only",
        "RECOGNITION_PLAUSIBLE": "mutant-specific recognition plausibility in tested context",
        "ASSAY_SPECIFIC_IMMUNOGENICITY_SUPPORTED": "assay-specific immunogenicity in tested context",
        "FUNCTIONAL_SPECIFICITY_SUPPORTED": "functional specificity in tested HLA/mutation system",
        "RECOGNITION_WITHOUT_ACTIVATION": "recognition without productive activation; no immunogenicity claim",
        "SPECIFICITY_RISK_HOLD_WT_POSITIVE": "no mutant-specific claim; WT cross-reactivity risk",
        "SPECIFICITY_RISK_HOLD_DECOY_POSITIVE": "no specificity claim; decoy/nonspecific signal risk",
    }.get(state, "manual review")
    next_step = {
        "PRESENTATION_FAILED": "reagent QC or deprioritize",
        "PRESENTATION_SUPPORTED": "mutant/WT/decoy multimer or TCR reporter",
        "RECOGNITION_PLAUSIBLE": "ELISpot/ICS/CD137 activation assay",
        "ASSAY_SPECIFIC_IMMUNOGENICITY_SUPPORTED": "matched HLA/mutation killing assay",
        "FUNCTIONAL_SPECIFICITY_SUPPORTED": "independent replicate or external validation",
        "RECOGNITION_WITHOUT_ACTIVATION": "dose response, alternate TCR context, or downgrade",
        "SPECIFICITY_RISK_HOLD_WT_POSITIVE": "hold candidate; investigate WT cross-reactivity",
        "SPECIFICITY_RISK_HOLD_DECOY_POSITIVE": "repeat with orthogonal decoy and artifact controls",
    }.get(state, "manual review")
    return post, uncertainty, state, allowed, next_step


def build_scenarios() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    readiness = read_tsv(IMMUNO / "cross_neo_i_candidate_readiness.tsv")
    preclin = read_tsv(PRECLIN / "preclinical_readiness_scorecard.tsv")
    if readiness.empty:
        raise FileNotFoundError(IMMUNO / "cross_neo_i_candidate_readiness.tsv")
    top = readiness.sort_values("cross_neo_i_readiness_score", ascending=False).head(13)
    pre_map = {
        str(r.get("source_row_id", "")): r
        for _, r in preclin.iterrows()
    } if not preclin.empty else {}
    rows = []
    assay_rows = []
    for _, cand in top.iterrows():
        prior = nfloat(cand.get("cross_neo_i_readiness_score"))
        pre = pre_map.get(str(cand.get("row_id", "")), {})
        wt = pre.get("wildtype_peptide", "")
        decoy = pre.get("decoy_peptide", "")
        for sc in scenario_defs():
            post, uncertainty, state, allowed, next_step = posterior(prior, sc)
            delta = post - prior
            rows.append(
                {
                    "simulation_flag": "SIMULATED_NOT_REAL_ASSAY_LABEL",
                    "row_id": cand.get("row_id", ""),
                    "peptide": cand.get("peptide", ""),
                    "hla_4digit": cand.get("hla_4digit", ""),
                    "wildtype_peptide": wt,
                    "decoy_peptide": decoy,
                    "scenario": sc["scenario"],
                    "scenario_meaning": sc["meaning"],
                    "prior_cross_neo_i_readiness": prior,
                    "scenario_posterior": post,
                    "posterior_delta": delta,
                    "scenario_uncertainty": uncertainty,
                    "claim_state": state,
                    "allowed_claim_if_real": allowed,
                    "next_step_if_real": next_step,
                    "presentation_pass": sc["presentation"],
                    "recognition_pass": sc["recognition"],
                    "activation_pass": sc["activation"],
                    "killing_pass": sc["killing"],
                    "specificity_pass": sc["specificity"],
                    "wt_positive": sc["wt_risk"],
                    "decoy_positive": sc["decoy_risk"],
                }
            )
            for assay_type in ["hla_stability", "multimer", "elispot", "killing"]:
                if assay_type == "hla_stability" and math.isnan(nfloat(sc["presentation"], math.nan)):
                    continue
                if assay_type == "multimer" and math.isnan(nfloat(sc["recognition"], math.nan)):
                    continue
                if assay_type == "elispot" and math.isnan(nfloat(sc["activation"], math.nan)):
                    continue
                if assay_type == "killing" and math.isnan(nfloat(sc["killing"], math.nan)):
                    continue
                assay_rows.append(
                    {
                        "simulation_flag": "SIMULATED_NOT_REAL_ASSAY_LABEL",
                        "row_id": cand.get("row_id", ""),
                        "peptide": cand.get("peptide", ""),
                        "hla_4digit": cand.get("hla_4digit", ""),
                        "wildtype_peptide": wt,
                        "decoy_peptide": decoy,
                        "scenario": sc["scenario"],
                        "assay_type": assay_type,
                        "mutant_positive": bool(sc["presentation"] if assay_type == "hla_stability" else sc["recognition"] if assay_type == "multimer" else sc["activation"] if assay_type == "elispot" else sc["killing"]),
                        "wt_positive": bool(sc["wt_risk"]),
                        "decoy_positive": bool(sc["decoy_risk"]),
                        "specificity_pass": bool(sc["specificity"]) if not math.isnan(nfloat(sc["specificity"], math.nan)) else "",
                        "notes": "what-if row; do not copy into real assay_feedback.tsv as observed data",
                    }
                )
    outcomes = pd.DataFrame(rows)
    simulated_assay = pd.DataFrame(assay_rows)
    transition = (
        outcomes.groupby(["scenario", "claim_state", "allowed_claim_if_real"], dropna=False)
        .agg(n_candidates=("row_id", "count"), mean_posterior=("scenario_posterior", "mean"), mean_delta=("posterior_delta", "mean"))
        .reset_index()
        .sort_values(["scenario", "mean_posterior"], ascending=[True, False])
    )
    return outcomes, simulated_assay, transition


def make_figures(outcomes: pd.DataFrame, transition: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    focus = outcomes[outcomes["peptide"].isin(["GADGVGKSAL", "HMTEVVRHC"])].copy()
    if not focus.empty:
        piv = focus.pivot(index="scenario", columns="peptide", values="scenario_posterior")
        fig, ax = plt.subplots(figsize=(10.5, 5.4))
        im = ax.imshow(piv.to_numpy(), aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
        ax.set_yticks(np.arange(len(piv.index)), [x.replace("_", "\n") for x in piv.index], fontsize=8)
        ax.set_xticks(np.arange(len(piv.columns)), piv.columns)
        ax.set_title("What-if posterior by assay scenario")
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                ax.text(j, i, f"{piv.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        fig.tight_layout()
        fig.savefig(FIG / "fig_md60_assay_scenario_posterior_heatmap.png", dpi=220)
        fig.savefig(FIG / "fig_md60_assay_scenario_posterior_heatmap.pdf")
        plt.close(fig)

    claim_counts = transition.pivot_table(index="scenario", columns="claim_state", values="n_candidates", fill_value=0)
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    bottom = np.zeros(len(claim_counts))
    for col in claim_counts.columns:
        ax.bar(np.arange(len(claim_counts)), claim_counts[col], bottom=bottom, label=col)
        bottom += claim_counts[col].to_numpy()
    ax.set_xticks(np.arange(len(claim_counts)), [x.replace("_", "\n") for x in claim_counts.index], rotation=0, fontsize=7)
    ax.set_ylabel("candidate count")
    ax.set_title("Claim-state transition under what-if assay outcomes")
    ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "fig_md61_assay_scenario_claim_transitions.png", dpi=220)
    fig.savefig(FIG / "fig_md61_assay_scenario_claim_transitions.pdf")
    plt.close(fig)

    mean_delta = outcomes.groupby("scenario")["posterior_delta"].mean().sort_values()
    fig, ax = plt.subplots(figsize=(10.5, 5.0))
    colors = ["#e15759" if v < 0 else "#59a14f" for v in mean_delta.values]
    ax.barh(np.arange(len(mean_delta)), mean_delta.values, color=colors)
    ax.set_yticks(np.arange(len(mean_delta)), [x.replace("_", " ") for x in mean_delta.index])
    ax.set_xlabel("mean posterior delta")
    ax.set_title("Scenario impact on posterior readiness")
    ax.axvline(0, color="black", lw=1)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md62_scenario_posterior_delta.png", dpi=220)
    fig.savefig(FIG / "fig_md62_scenario_posterior_delta.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.2, 4.7))
    ax.axis("off")
    boxes = [
        ("simulate", "what-if\nassay outcome"),
        ("posterior", "candidate\nposterior"),
        ("claim", "allowed /\nblocked claim"),
        ("decision", "next wetlab\nmove"),
    ]
    xs = np.linspace(0.12, 0.88, len(boxes))
    for i, (title, body) in enumerate(boxes):
        ax.add_patch(plt.Rectangle((xs[i] - 0.08, 0.50), 0.16, 0.22, fc="#10243a", ec="#edf5ff", lw=1.2))
        ax.text(xs[i], 0.65, title, ha="center", color="#f2c46d", weight="bold", fontsize=10)
        ax.text(xs[i], 0.55, body, ha="center", color="white", fontsize=9)
        if i < len(boxes) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.09, 0.61), xytext=(xs[i] + 0.09, 0.61), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(0.5, 0.85, "Assay feedback scenario simulator", ha="center", fontsize=15, weight="bold")
    ax.text(0.5, 0.25, "All rows are simulated_not_real and cannot be used as training labels.", ha="center", color="#9b1c31", fontsize=10)
    fig.savefig(FIG / "fig_md63_assay_scenario_simulator_flow.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md63_assay_scenario_simulator_flow.pdf", bbox_inches="tight")
    plt.close(fig)


def write_reports(outcomes: pd.DataFrame, simulated_assay: pd.DataFrame, transition: pd.DataFrame) -> None:
    lines = [
        "# CROSS-Neo Assay Feedback Scenario Simulator",
        "",
        "## Purpose",
        "",
        "This package shows how the closed-loop learner will behave once real WT/decoy-controlled assay data arrive. These rows are explicitly simulated and must not be treated as observed labels.",
        "",
        "## Lead Scenario Outcomes",
        "",
        md_table(
            outcomes[outcomes["peptide"].isin(["GADGVGKSAL", "HMTEVVRHC"])],
            [
                "peptide",
                "hla_4digit",
                "scenario",
                "prior_cross_neo_i_readiness",
                "scenario_posterior",
                "posterior_delta",
                "claim_state",
                "allowed_claim_if_real",
                "next_step_if_real",
            ],
            30,
        ),
        "",
        "## Claim Transitions",
        "",
        md_table(transition, ["scenario", "claim_state", "allowed_claim_if_real", "n_candidates", "mean_posterior", "mean_delta"], 40),
        "",
        "## Boundary",
        "",
        "Use this for assay planning and reviewer explanation only. It is not a result and not a substitute for real activation/killing assays.",
    ]
    (OUT / "ASSAY_FEEDBACK_SCENARIO_SIMULATOR_REPORT.md").write_text("\n".join(lines) + "\n")

    kr = [
        "# assay feedback scenario simulator 요약",
        "",
        "이 패키지는 실제 실험값이 들어왔을 때 posterior와 claim이 어떻게 바뀌는지 미리 보여주는 what-if 시뮬레이터입니다.",
        "",
        "- mutant activation positive + WT/decoy negative: assay-specific immunogenicity claim까지 상승 가능",
        "- WT positive: mutant-specific claim 차단",
        "- decoy positive: nonspecific/artifact 위험으로 claim 차단",
        "- presentation positive만 있음: presentation plausibility까지만 가능",
        "- recognition positive but activation negative: TCR binding은 있어도 immunogenicity claim은 불가",
        "",
        "중요: 여기 rows는 전부 SIMULATED_NOT_REAL_ASSAY_LABEL입니다. 실제 label TSV와 분리되어 있습니다.",
    ]
    (OUT / "ASSAY_FEEDBACK_SCENARIO_ONE_PAGE_KR.md").write_text("\n".join(kr) + "\n")


def write_html(outcomes: pd.DataFrame, transition: pd.DataFrame, summary: dict) -> None:
    figures = "\n".join(
        f"""
        <article>
          <img src="assets/cross_neo_md_audit/{name}.png" alt="{title}">
          <h3>{title}</h3>
          <p>{desc}</p>
          <a href="assets/cross_neo_md_audit/{name}.png">PNG</a>
          <a href="assets/cross_neo_md_audit/{name}.pdf">PDF</a>
        </article>
        """
        for name, title, desc in [
            ("fig_md60_assay_scenario_posterior_heatmap", "Scenario Posterior Heatmap", "How flagship candidates move under different assay outcomes."),
            ("fig_md61_assay_scenario_claim_transitions", "Claim Transitions", "Allowed and blocked claim states under each scenario."),
            ("fig_md62_scenario_posterior_delta", "Posterior Delta", "Which scenarios raise or lower readiness the most."),
            ("fig_md63_assay_scenario_simulator_flow", "Simulator Flow", "Planning-only path from what-if assay result to next move."),
        ]
    )
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Assay Scenario Simulator</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5;text-decoration:none}} header{{padding:40px 30px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1360px;margin:0 auto;padding:24px 30px 70px}} h1{{font-family:Georgia,serif;font-size:44px;margin:0;color:#fff8e8}}
h2{{font-family:Georgia,serif;color:#fff2d0}} .lead,.muted{{color:#aebdd1}} section{{border-bottom:1px solid #26364d;padding:22px 0}}
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin-top:18px}} .stat{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:12px}}
.stat b{{display:block;color:#f2c46d;font-size:26px}} table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{border:1px solid #293b55;padding:7px;vertical-align:top}} th{{background:#13243a;color:#f2c46d}}
.grid2{{display:grid;grid-template-columns:repeat(2,1fr);gap:16px}} article{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:13px}} img{{width:100%;background:white;border-radius:6px}}
.warn{{border-left:4px solid #ff7b72;background:#111d2e;padding:12px 14px}} .links{{display:flex;flex-wrap:wrap;gap:10px}} .pill{{border:1px solid #395170;border-radius:999px;padding:7px 10px;background:#101d30}}
@media(max-width:900px){{.grid2{{grid-template-columns:1fr}} h1{{font-size:32px}}}}
</style></head><body>
<header><div class="wrap">
<h1>CROSS-Neo Assay Scenario Simulator</h1>
<p class="lead">Planning-only what-if layer for mutant/WT/decoy assay outcomes. It shows how posterior readiness, claim state, and next experiment would change.</p>
<div class="stats">
<div class="stat"><b>{summary['n_candidates']}</b><span>candidate rows</span></div>
<div class="stat"><b>{summary['n_scenarios']}</b><span>what-if scenarios</span></div>
<div class="stat"><b>{summary['n_outcome_rows']}</b><span>scenario outcomes</span></div>
<div class="stat"><b>{summary['n_simulated_assay_rows']}</b><span>simulated assay rows</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Boundary</h2>
<div class="warn">Every row here is SIMULATED_NOT_REAL_ASSAY_LABEL. Use this for planning and reviewer explanation only; do not train or claim from it.</div>
<div class="links">
<a class="pill" href="assets/cross_neo_md_audit/ASSAY_FEEDBACK_SCENARIO_SIMULATOR_REPORT.md">Full report</a>
<a class="pill" href="assets/cross_neo_md_audit/ASSAY_FEEDBACK_SCENARIO_ONE_PAGE_KR.md">Korean brief</a>
<a class="pill" href="cross_neo_assay_feedback_loop.html">Real assay feedback loop</a>
<a class="pill" href="cross_neo_physics_gate.html">Physics gate</a>
<a class="pill" href="assets/cross_neo_md_audit/assay_feedback_scenario_outcomes.tsv">scenario outcomes TSV</a>
</div></section>
<section><h2>Figures</h2><div class="grid2">{figures}</div></section>
<section><h2>Flagship Scenario Outcomes</h2>{table_html(outcomes[outcomes['peptide'].isin(['GADGVGKSAL', 'HMTEVVRHC'])], ['peptide','hla_4digit','scenario','prior_cross_neo_i_readiness','scenario_posterior','posterior_delta','claim_state','allowed_claim_if_real','next_step_if_real'], 40)}</section>
<section><h2>Claim Transitions</h2>{table_html(transition, ['scenario','claim_state','allowed_claim_if_real','n_candidates','mean_posterior','mean_delta'], 50)}</section>
</main></body></html>"""
    PAGE.write_text(html_text)


def deploy(paths: list[Path]) -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            safe_copy(path, ASSET / path.name)
            safe_copy(path, WEB_ASSET / path.name)
    for fig in [
        "fig_md60_assay_scenario_posterior_heatmap",
        "fig_md61_assay_scenario_claim_transitions",
        "fig_md62_scenario_posterior_delta",
        "fig_md63_assay_scenario_simulator_flow",
    ]:
        for ext in [".png", ".pdf"]:
            path = FIG / f"{fig}{ext}"
            if path.exists():
                safe_copy(path, ASSET / path.name)
                safe_copy(path, WEB_ASSET / path.name)
    WEB.mkdir(parents=True, exist_ok=True)
    safe_copy(PAGE, WEB_PAGE)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    outcomes, simulated_assay, transition = build_scenarios()
    outcomes.to_csv(OUT / "assay_feedback_scenario_outcomes.tsv", sep="\t", index=False)
    simulated_assay.to_csv(OUT / "simulated_assay_feedback_rows.tsv", sep="\t", index=False)
    transition.to_csv(OUT / "assay_feedback_claim_transitions.tsv", sep="\t", index=False)
    make_figures(outcomes, transition)
    write_reports(outcomes, simulated_assay, transition)
    summary = {
        "n_candidates": int(outcomes["row_id"].nunique()),
        "n_scenarios": int(outcomes["scenario"].nunique()),
        "n_outcome_rows": int(len(outcomes)),
        "n_simulated_assay_rows": int(len(simulated_assay)),
        "figures": [
            "fig_md60_assay_scenario_posterior_heatmap",
            "fig_md61_assay_scenario_claim_transitions",
            "fig_md62_scenario_posterior_delta",
            "fig_md63_assay_scenario_simulator_flow",
        ],
        "boundary": "simulated what-if scenarios only; not observed assay labels",
    }
    (OUT / "assay_feedback_scenario_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_html(outcomes, transition, summary)
    deploy(
        [
            OUT / "assay_feedback_scenario_outcomes.tsv",
            OUT / "simulated_assay_feedback_rows.tsv",
            OUT / "assay_feedback_claim_transitions.tsv",
            OUT / "ASSAY_FEEDBACK_SCENARIO_SIMULATOR_REPORT.md",
            OUT / "ASSAY_FEEDBACK_SCENARIO_ONE_PAGE_KR.md",
            OUT / "assay_feedback_scenario_summary.json",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[assay-scenarios] wrote {OUT}")
    print(f"[assay-scenarios] page {PAGE}")


if __name__ == "__main__":
    main()
