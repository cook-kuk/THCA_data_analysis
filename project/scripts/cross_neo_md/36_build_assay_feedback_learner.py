#!/usr/bin/env python3
"""Build a data-driven assay feedback learner for CROSS-Neo.

This is a closed-loop decision layer. It does not invent immunogenicity labels.
If real WT/decoy-controlled assay rows are present, it updates candidate
posteriors and claim states. If no assay rows exist yet, it emits templates,
priors, next-experiment utilities, and a web page that is ready for data entry.
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
PKG = MD_OUT / "high_impact_decision_package"
PRECLIN = MD_OUT / "preclinical_validation_protocol"
IMMUNO = MD_OUT / "immunogenicity_prediction_upgrade"
OUT = MD_OUT / "assay_feedback_learner"
FIG = MD_OUT / "figures"
DATA_IN = REPO / "project/data/cross_neo_assay_feedback/assay_feedback.tsv"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_assay_feedback_loop.html"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
WEB_PAGE = WEB / "cross_neo_assay_feedback_loop.html"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def nfloat(value: object, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default


def boolish(value: object) -> float:
    if pd.isna(value):
        return math.nan
    if isinstance(value, (int, float)):
        return 1.0 if float(value) > 0 else 0.0
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "pass", "positive", "pos"}:
        return 1.0
    if text in {"0", "false", "no", "n", "fail", "negative", "neg"}:
        return 0.0
    return math.nan


def esc(value: object) -> str:
    return html.escape("" if pd.isna(value) else str(value))


def fmt(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.{digits}f}"
    except Exception:
        return esc(value)


def table_html(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    if df.empty:
        return "<p class='muted'>No rows.</p>"
    keep = [c for c in cols if c in df.columns]
    rows = ["<table><thead><tr>"]
    rows.extend(f"<th>{esc(c.replace('_', ' '))}</th>" for c in keep)
    rows.append("</tr></thead><tbody>")
    for _, row in df.head(n).iterrows():
        rows.append("<tr>")
        for c in keep:
            value = row.get(c, "")
            rows.append(f"<td>{fmt(value) if isinstance(value, (int, float)) else esc(value)}</td>")
        rows.append("</tr>")
    rows.append("</tbody></table>")
    return "\n".join(rows)


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
        print(f"[assay-feedback] skip permission-denied copy {dst}")


def build_schema() -> pd.DataFrame:
    rows = [
        ("row_id", "string", "CROSS-Neo candidate row_id. Required when available."),
        ("peptide", "string", "Mutant peptide sequence."),
        ("hla_4digit", "string", "HLA allele, e.g. HLA-A*02:01."),
        ("wildtype_peptide", "string", "Matched WT peptide. Required for specificity claim."),
        ("decoy_peptide", "string", "Anchor-preserved or scrambled decoy peptide."),
        ("assay_date", "string", "YYYY-MM-DD if available."),
        ("batch_id", "string", "Assay batch or plate identifier."),
        ("assay_type", "enum", "hla_stability|multimer|tcr_reporter|elispot|ics|activation_marker|killing."),
        ("antigen_type", "enum", "mutant|wildtype|decoy|irrelevant|positive_control|background."),
        ("mutant_signal", "float", "Mutant condition signal in assay units."),
        ("wt_signal", "float", "Matched WT signal in same units."),
        ("decoy_signal", "float", "Decoy signal in same units."),
        ("irrelevant_signal", "float", "Irrelevant pHLA/peptide control signal."),
        ("background_signal", "float", "No peptide, vehicle, unstimulated, or assay background signal."),
        ("positive_control_signal", "float", "Positive control signal."),
        ("replicate_n", "integer", "Number of biological or technical replicates summarized."),
        ("signal_units", "string", "MFI, spots/1e5 cells, %positive, %lysis, etc."),
        ("pass_threshold", "float", "Pre-specified pass threshold in signal units or fold-change."),
        ("mutant_positive", "bool", "Whether mutant condition passes assay-specific threshold."),
        ("wt_positive", "bool", "Whether WT condition is positive. Positive WT blocks mutant-specific claim."),
        ("decoy_positive", "bool", "Whether decoy condition is positive. Positive decoy blocks specificity claim."),
        ("specificity_pass", "bool", "Mutant exceeds WT/decoy/irrelevant under matched conditions."),
        ("presentation_pass", "bool", "HLA binding/stability/MS presentation pass."),
        ("recognition_pass", "bool", "TCR/multimer/reporter recognition pass."),
        ("activation_pass", "bool", "ELISpot/ICS/CD137/cytokine activation pass."),
        ("killing_pass", "bool", "Target-cell killing pass."),
        ("notes", "string", "Free-text caveat. Do not use this as a model label."),
    ]
    return pd.DataFrame(rows, columns=["column", "type", "description"])


def load_inputs() -> dict[str, pd.DataFrame]:
    data = {
        "top13": read_tsv(PKG / "no_false_positive_top13_candidates.tsv"),
        "preclin": read_tsv(PRECLIN / "preclinical_readiness_scorecard.tsv"),
        "immuno": read_tsv(IMMUNO / "cross_neo_i_candidate_readiness.tsv"),
        "assay": read_tsv(DATA_IN),
    }
    if data["top13"].empty or data["immuno"].empty:
        raise FileNotFoundError("Missing Top-13 or CROSS-Neo-I readiness inputs")
    return data


def build_template(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    top = data["top13"].copy()
    pre = data["preclin"].copy()
    rows = []
    for _, row in top.sort_values("priority_order").head(13).iterrows():
        pre_row = pre[pre["source_row_id"].astype(str).eq(str(row.get("row_id", "")))]
        wt = row.get("tentative_wt_sequence", row.get("wildtype_peptide", ""))
        decoy = row.get("anchor_preserved_decoy", "")
        if not pre_row.empty:
            wt = pre_row.iloc[0].get("wildtype_peptide", wt)
            decoy = pre_row.iloc[0].get("decoy_peptide", decoy)
        rows.append(
            {
                "row_id": row.get("row_id", ""),
                "peptide": row.get("peptide", ""),
                "hla_4digit": row.get("hla_4digit", ""),
                "wildtype_peptide": wt,
                "decoy_peptide": decoy,
                "assay_date": "",
                "batch_id": "",
                "assay_type": "",
                "antigen_type": "mutant/wildtype/decoy/irrelevant/positive_control/background",
                "mutant_signal": "",
                "wt_signal": "",
                "decoy_signal": "",
                "irrelevant_signal": "",
                "background_signal": "",
                "positive_control_signal": "",
                "replicate_n": "",
                "signal_units": "",
                "pass_threshold": "",
                "mutant_positive": "",
                "wt_positive": "",
                "decoy_positive": "",
                "specificity_pass": "",
                "presentation_pass": "",
                "recognition_pass": "",
                "activation_pass": "",
                "killing_pass": "",
                "notes": "fill this row with real assay data; blank rows are not used as labels",
            }
        )
    return pd.DataFrame(rows)


def summarize_assays(assay: pd.DataFrame) -> pd.DataFrame:
    if assay.empty:
        return pd.DataFrame()
    assay = assay.copy()
    if "simulation_flag" in assay.columns:
        assay = assay[~assay["simulation_flag"].astype(str).str.upper().str.contains("SIMULATED", na=False)].copy()
    if "notes" in assay.columns:
        assay = assay[~assay["notes"].astype(str).str.contains("blank rows are not used as labels|fill this row", case=False, na=False)].copy()
    if "assay_type" in assay.columns:
        assay_type = assay["assay_type"].astype(str).str.strip()
        assay = assay[assay_type.ne("") & assay_type.str.lower().ne("nan")].copy()
    evidence_cols = [
        "mutant_signal",
        "wt_signal",
        "decoy_signal",
        "irrelevant_signal",
        "mutant_positive",
        "wt_positive",
        "decoy_positive",
        "specificity_pass",
        "presentation_pass",
        "recognition_pass",
        "activation_pass",
        "killing_pass",
    ]
    present_cols = [c for c in evidence_cols if c in assay.columns]
    if present_cols:
        has_evidence = pd.Series(False, index=assay.index)
        for col in present_cols:
            values = assay[col]
            has_evidence = has_evidence | (values.notna() & values.astype(str).str.strip().ne(""))
        assay = assay[has_evidence].copy()
    if assay.empty:
        return pd.DataFrame()
    for col in [
        "mutant_positive",
        "wt_positive",
        "decoy_positive",
        "specificity_pass",
        "presentation_pass",
        "recognition_pass",
        "activation_pass",
        "killing_pass",
    ]:
        if col in assay.columns:
            assay[col + "_num"] = assay[col].map(boolish)
    rows = []
    group_cols = [c for c in ["row_id", "peptide", "hla_4digit"] if c in assay.columns]
    for keys, sub in assay.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        item = dict(zip(group_cols, keys))
        signal_cols = ["mutant_signal", "wt_signal", "decoy_signal", "irrelevant_signal", "background_signal"]
        for col in signal_cols:
            if col in sub.columns:
                item[col + "_mean"] = pd.to_numeric(sub[col], errors="coerce").mean()
        for col in [
            "mutant_positive",
            "wt_positive",
            "decoy_positive",
            "specificity_pass",
            "presentation_pass",
            "recognition_pass",
            "activation_pass",
            "killing_pass",
        ]:
            ncol = col + "_num"
            if ncol in sub.columns:
                item[col + "_rate"] = sub[ncol].mean()
                item[col + "_observed_n"] = int(sub[ncol].notna().sum())
        item["assay_rows"] = len(sub)
        item["assay_types"] = "|".join(sorted(set(sub.get("assay_type", pd.Series(dtype=str)).dropna().astype(str))))
        rows.append(item)
    return pd.DataFrame(rows)


def posterior_from_assay(prior: float, summary: pd.Series | None) -> tuple[float, float, str, dict[str, float]]:
    alpha = 1.0 + prior * 6.0
    beta = 1.0 + (1.0 - prior) * 6.0
    evidence = {
        "presentation": math.nan,
        "recognition": math.nan,
        "activation": math.nan,
        "killing": math.nan,
        "specificity": math.nan,
        "wt_risk": math.nan,
        "decoy_risk": math.nan,
    }
    if summary is not None and not summary.empty:
        weights = {
            "presentation_pass_rate": 1.0,
            "recognition_pass_rate": 1.5,
            "activation_pass_rate": 3.0,
            "killing_pass_rate": 4.0,
            "specificity_pass_rate": 2.5,
        }
        names = {
            "presentation_pass_rate": "presentation",
            "recognition_pass_rate": "recognition",
            "activation_pass_rate": "activation",
            "killing_pass_rate": "killing",
            "specificity_pass_rate": "specificity",
        }
        for col, weight in weights.items():
            value = nfloat(summary.get(col), math.nan)
            if not math.isnan(value):
                evidence[names[col]] = value
                alpha += value * weight
                beta += (1.0 - value) * weight
        wt_risk = nfloat(summary.get("wt_positive_rate"), math.nan)
        decoy_risk = nfloat(summary.get("decoy_positive_rate"), math.nan)
        evidence["wt_risk"] = wt_risk
        evidence["decoy_risk"] = decoy_risk
        if not math.isnan(wt_risk):
            beta += wt_risk * 3.0
        if not math.isnan(decoy_risk):
            beta += decoy_risk * 3.0
    posterior = alpha / (alpha + beta)
    uncertainty = posterior * (1.0 - posterior)
    if summary is None or summary.empty:
        state = "PRIOR_ONLY_NO_ASSAY_LABEL"
    elif not math.isnan(evidence["wt_risk"]) and evidence["wt_risk"] > 0:
        state = "SPECIFICITY_RISK_HOLD_WT_POSITIVE"
    elif not math.isnan(evidence["decoy_risk"]) and evidence["decoy_risk"] > 0:
        state = "SPECIFICITY_RISK_HOLD_DECOY_POSITIVE"
    elif evidence["killing"] == 1.0 and evidence["specificity"] == 1.0:
        state = "FUNCTIONAL_SPECIFICITY_SUPPORTED"
    elif evidence["activation"] == 1.0 and evidence["specificity"] == 1.0:
        state = "ASSAY_SPECIFIC_IMMUNOGENICITY_SUPPORTED"
    elif evidence["recognition"] == 1.0 and evidence["specificity"] == 1.0:
        state = "RECOGNITION_PLAUSIBLE"
    elif evidence["presentation"] == 1.0:
        state = "PRESENTATION_SUPPORTED"
    else:
        state = "ASSAY_INCONCLUSIVE_OR_NEGATIVE"
    return posterior, uncertainty, state, evidence


def next_experiment(row: pd.Series, evidence: dict[str, float], has_assay: bool) -> tuple[str, float, str]:
    readiness = nfloat(row.get("cross_neo_i_readiness_score"), 0.0)
    p = nfloat(row.get("posterior_immunogenicity_readiness"), readiness)
    uncertainty = p * (1 - p)
    if not has_assay or math.isnan(evidence["presentation"]):
        return "WT/decoy-controlled HLA stability or pHLA binding", 0.25 + readiness + uncertainty, "cheapest gate before T-cell assays"
    if evidence["presentation"] < 1.0:
        return "hold or repeat presentation with reagent QC", 0.15 + uncertainty, "presentation failed or is unclear"
    if math.isnan(evidence["recognition"]):
        return "mutant/WT/decoy multimer or TCR reporter", 0.35 + readiness + uncertainty, "tests recognition without claiming activation"
    if evidence["recognition"] < 1.0:
        return "downgrade or search alternate TCR context", 0.20 + uncertainty, "recognition gate failed"
    if math.isnan(evidence["activation"]):
        return "ELISpot/ICS/CD137 activation assay", 0.50 + readiness + uncertainty, "first assay that can support immunogenicity claim"
    if evidence["activation"] < 1.0:
        return "dose-response repeat or downgrade to recognition-only", 0.25 + uncertainty, "activation gate failed"
    if math.isnan(evidence["killing"]):
        return "matched HLA/mutation target-cell killing assay", 0.45 + readiness + uncertainty, "functional follow-up after activation"
    return "independent replicate or external heldout validation", 0.30 + uncertainty, "strengthens reproducibility and manuscript claim"


def build_posteriors(data: dict[str, pd.DataFrame], assay_summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = data["immuno"].copy()
    if base.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    summaries = {}
    if not assay_summary.empty:
        for _, row in assay_summary.iterrows():
            key = str(row.get("row_id", "")) or f"{row.get('peptide','')}|{row.get('hla_4digit','')}"
            summaries[key] = row
    rows = []
    claims = []
    active = []
    for _, row in base.iterrows():
        key = str(row.get("row_id", "")) or f"{row.get('peptide','')}|{row.get('hla_4digit','')}"
        summary = summaries.get(key)
        has_assay = summary is not None
        prior = nfloat(row.get("cross_neo_i_readiness_score"), 0.0)
        posterior, uncertainty, state, evidence = posterior_from_assay(prior, summary)
        exp, utility, reason = next_experiment(pd.Series({**row.to_dict(), "posterior_immunogenicity_readiness": posterior}), evidence, has_assay)
        update = {
            "row_id": row.get("row_id", ""),
            "peptide": row.get("peptide", ""),
            "hla_4digit": row.get("hla_4digit", ""),
            "source_dataset": row.get("source_dataset", ""),
            "prior_cross_neo_i_readiness": prior,
            "posterior_immunogenicity_readiness": posterior,
            "posterior_uncertainty": uncertainty,
            "assay_rows_observed": int(summary.get("assay_rows", 0)) if has_assay else 0,
            "claim_state": state,
            "presentation_evidence": evidence["presentation"],
            "recognition_evidence": evidence["recognition"],
            "activation_evidence": evidence["activation"],
            "killing_evidence": evidence["killing"],
            "specificity_evidence": evidence["specificity"],
            "wt_positive_risk": evidence["wt_risk"],
            "decoy_positive_risk": evidence["decoy_risk"],
            "next_experiment": exp,
            "active_learning_utility": utility,
            "utility_reason": reason,
            "claim_boundary": "posterior is decision-support readiness, not clinical or universal immunogenicity probability",
        }
        rows.append(update)
        claims.append(
            {
                "row_id": row.get("row_id", ""),
                "peptide": row.get("peptide", ""),
                "hla_4digit": row.get("hla_4digit", ""),
                "claim_state": state,
                "allowed_claim_now": claim_allowed(state),
                "blocked_claim": "clinical utility; universal immunogenicity; MD-only proof",
                "required_next_gate": exp,
            }
        )
        active.append(update)
    post = pd.DataFrame(rows).sort_values(["active_learning_utility", "posterior_immunogenicity_readiness"], ascending=False)
    claim_df = pd.DataFrame(claims)
    active_df = pd.DataFrame(active).sort_values("active_learning_utility", ascending=False)
    return post, claim_df, active_df


def claim_allowed(state: str) -> str:
    return {
        "PRIOR_ONLY_NO_ASSAY_LABEL": "computational prioritization only",
        "PRESENTATION_SUPPORTED": "presentation plausibility in tested assay",
        "RECOGNITION_PLAUSIBLE": "recognition plausibility in tested TCR/cell context",
        "ASSAY_SPECIFIC_IMMUNOGENICITY_SUPPORTED": "assay-specific immunogenicity in tested context",
        "FUNCTIONAL_SPECIFICITY_SUPPORTED": "functional specificity in tested HLA/mutation system",
        "SPECIFICITY_RISK_HOLD_WT_POSITIVE": "no mutant-specific claim; WT cross-reactivity risk",
        "SPECIFICITY_RISK_HOLD_DECOY_POSITIVE": "no specificity claim; decoy/nonspecific signal risk",
        "ASSAY_INCONCLUSIVE_OR_NEGATIVE": "hold or downgrade until repeat/orthogonal assay",
    }.get(state, "claim requires manual review")


def build_rules() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ("mutant positive, WT negative, decoy negative", "increase posterior; allow next gate", "mutant-specific signal"),
            ("mutant positive, WT positive", "penalize posterior; hold specificity claim", "WT cross-reactivity risk"),
            ("mutant positive, decoy positive", "penalize posterior; repeat with orthogonal decoy", "nonspecific or motif artifact"),
            ("presentation positive, activation negative", "downgrade to presentation-only", "presentation is not recognition/activation"),
            ("recognition positive, activation negative", "mark as binding without productive activation", "TCR binding may be nonproductive"),
            ("activation positive, killing missing", "promote to killing assay", "activation is not functional killing"),
            ("killing positive with clean controls", "promote to functional specificity claim in tested context", "strongest preclinical claim"),
        ],
        columns=["observed_pattern", "model_update", "interpretation"],
    )


def make_figures(post: pd.DataFrame, active: pd.DataFrame, claim_df: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    fig, ax = plt.subplots(figsize=(11.5, 4.5))
    ax.axis("off")
    nodes = [
        ("CROSS-Neo\nprior", "DL/TCR/MD\nreadiness"),
        ("Assay TSV", "mutant/WT/decoy\nsignals"),
        ("Bayesian\nupdate", "posterior +\nuncertainty"),
        ("Claim\nstate", "allowed vs\nblocked"),
        ("Next\nexperiment", "active-learning\nqueue"),
    ]
    xs = np.linspace(0.08, 0.92, len(nodes))
    for i, (title, body) in enumerate(nodes):
        ax.add_patch(plt.Rectangle((xs[i] - 0.075, 0.48), 0.15, 0.22, fc="#10243a", ec="#edf5ff", lw=1.2))
        ax.text(xs[i], 0.63, title, ha="center", va="center", color="#f2c46d", fontsize=10, weight="bold")
        ax.text(xs[i], 0.53, body, ha="center", va="center", color="white", fontsize=8.5)
        if i < len(nodes) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.085, 0.59), xytext=(xs[i] + 0.085, 0.59), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(0.5, 0.84, "Closed-loop assay feedback learner", ha="center", fontsize=15, weight="bold")
    ax.text(0.5, 0.22, "Real assay rows update the model. Blank templates never become labels.", ha="center", color="#9b1c31", fontsize=10)
    fig.savefig(FIG / "fig_md56_assay_feedback_closed_loop.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md56_assay_feedback_closed_loop.pdf", bbox_inches="tight")
    plt.close(fig)

    if not post.empty:
        show = post.sort_values("posterior_immunogenicity_readiness", ascending=False).head(10)
        y = np.arange(len(show))
        fig, ax = plt.subplots(figsize=(10.8, 5.2))
        ax.barh(y + 0.18, show["prior_cross_neo_i_readiness"], height=0.34, label="prior", color="#4e79a7")
        ax.barh(y - 0.18, show["posterior_immunogenicity_readiness"], height=0.34, label="posterior", color="#59a14f")
        ax.set_yticks(y, show["peptide"] + "\n" + show["hla_4digit"].astype(str))
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("decision-support readiness")
        ax.set_title("Prior-to-posterior candidate update")
        ax.legend()
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(FIG / "fig_md57_prior_posterior_update.png", dpi=220)
        fig.savefig(FIG / "fig_md57_prior_posterior_update.pdf")
        plt.close(fig)

    if not active.empty:
        show = active.head(10)
        fig, ax = plt.subplots(figsize=(11, 5.4))
        ax.barh(np.arange(len(show)), show["active_learning_utility"], color="#f28e2b")
        ax.set_yticks(np.arange(len(show)), show["peptide"] + "\n" + show["next_experiment"].str.slice(0, 36))
        ax.invert_yaxis()
        ax.set_xlabel("active-learning utility")
        ax.set_title("Next experiment queue")
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(FIG / "fig_md58_active_learning_queue.png", dpi=220)
        fig.savefig(FIG / "fig_md58_active_learning_queue.pdf")
        plt.close(fig)

    if not claim_df.empty:
        counts = claim_df["claim_state"].value_counts()
        fig, ax = plt.subplots(figsize=(10.5, 4.8))
        ax.barh(np.arange(len(counts)), counts.values, color="#76b7b2")
        ax.set_yticks(np.arange(len(counts)), counts.index)
        ax.invert_yaxis()
        ax.set_xlabel("candidate count")
        ax.set_title("Claim-state distribution")
        fig.tight_layout()
        fig.savefig(FIG / "fig_md59_claim_state_distribution.png", dpi=220)
        fig.savefig(FIG / "fig_md59_claim_state_distribution.pdf")
        plt.close(fig)


def write_report(
    schema: pd.DataFrame,
    template: pd.DataFrame,
    assay_summary: pd.DataFrame,
    post: pd.DataFrame,
    claim_df: pd.DataFrame,
    active: pd.DataFrame,
    rules: pd.DataFrame,
) -> None:
    assay_present = not assay_summary.empty
    lines = [
        "# CROSS-Neo Assay Feedback Learner",
        "",
        "## Decision",
        "",
        "This package makes the neoantigen workflow data-driven. CROSS-Neo produces priors; WT/decoy-controlled assay rows update posterior readiness, claim state, and the next experiment queue.",
        "",
        "## Current Assay Data Status",
        "",
        f"- Assay feedback file: `{DATA_IN}`",
        f"- Real assay summaries detected: {len(assay_summary)}",
        f"- Mode: {'posterior update from assay data' if assay_present else 'prior-only template and active-learning plan'}",
        "",
        "## Candidate Posterior Updates",
        "",
        md_table(
            post,
            [
                "row_id",
                "peptide",
                "hla_4digit",
                "prior_cross_neo_i_readiness",
                "posterior_immunogenicity_readiness",
                "posterior_uncertainty",
                "assay_rows_observed",
                "claim_state",
                "next_experiment",
                "active_learning_utility",
            ],
            20,
        ),
        "",
        "## Claim State Updates",
        "",
        md_table(claim_df, ["row_id", "peptide", "hla_4digit", "claim_state", "allowed_claim_now", "required_next_gate"], 20),
        "",
        "## Next Experiment Queue",
        "",
        md_table(active, ["peptide", "hla_4digit", "next_experiment", "active_learning_utility", "utility_reason"], 20),
        "",
        "## Feedback Rules",
        "",
        md_table(rules, ["observed_pattern", "model_update", "interpretation"], 20),
        "",
        "## Required Input Schema",
        "",
        md_table(schema, ["column", "type", "description"], 40),
        "",
        "## Boundary",
        "",
        "The posterior is a decision-support readiness score for choosing the next experiment. It is not a clinical probability and not universal immunogenicity proof.",
    ]
    (OUT / "ASSAY_FEEDBACK_LEARNER_REPORT.md").write_text("\n".join(lines) + "\n")

    kr = [
        "# CROSS-Neo assay feedback learner 요약",
        "",
        "이제 구조가 실제로 data-driven closed loop가 됩니다.",
        "",
        "1. CROSS-Neo가 후보 prior를 냅니다.",
        "2. mutant/WT/decoy assay TSV를 넣습니다.",
        "3. posterior readiness, claim state, 다음 실험 우선순위가 자동으로 바뀝니다.",
        "4. WT positive 또는 decoy positive이면 mutant-specific claim이 자동으로 막힙니다.",
        "5. activation positive + WT/decoy negative이면 assay-specific immunogenicity claim으로 올라갈 수 있습니다.",
        "",
        f"현재 모드: {'실제 assay label 반영 중' if assay_present else '아직 실제 assay label 없음; 템플릿과 다음 실험 큐 생성'}",
        "",
        "중요: 이 점수는 실험 우선순위용 posterior readiness입니다. 임상 확률이나 면역원성 확정이 아닙니다.",
    ]
    (OUT / "ASSAY_FEEDBACK_ONE_PAGE_KR.md").write_text("\n".join(kr) + "\n")


def write_html(post: pd.DataFrame, claim_df: pd.DataFrame, active: pd.DataFrame, rules: pd.DataFrame, summary: dict) -> None:
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
            ("fig_md56_assay_feedback_closed_loop", "Closed Loop", "How CROSS-Neo priors, assay rows, posterior updates, and next experiments connect."),
            ("fig_md57_prior_posterior_update", "Prior-To-Posterior", "Candidate readiness before and after available assay evidence."),
            ("fig_md58_active_learning_queue", "Next Experiment Queue", "Which assay should be run next for each candidate."),
            ("fig_md59_claim_state_distribution", "Claim State", "Current allowed claim state across candidates."),
        ]
    )
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Assay Feedback Loop</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5;text-decoration:none}} header{{padding:40px 30px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1360px;margin:0 auto;padding:24px 30px 70px}} h1{{font-family:Georgia,serif;font-size:44px;margin:0;color:#fff8e8}}
h2{{font-family:Georgia,serif;color:#fff2d0}} .lead,.muted{{color:#aebdd1}} section{{border-bottom:1px solid #26364d;padding:22px 0}}
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin-top:18px}} .stat{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:12px}}
.stat b{{display:block;color:#f2c46d;font-size:26px}} table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{border:1px solid #293b55;padding:7px;vertical-align:top}} th{{background:#13243a;color:#f2c46d}}
.grid2{{display:grid;grid-template-columns:repeat(2,1fr);gap:16px}} article{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:13px}} img{{width:100%;background:white;border-radius:6px}}
.ok{{border-left:4px solid #39d4b5;background:#111d2e;padding:12px 14px}} .warn{{border-left:4px solid #ff7b72;background:#111d2e;padding:12px 14px}}
.links{{display:flex;flex-wrap:wrap;gap:10px}} .pill{{border:1px solid #395170;border-radius:999px;padding:7px 10px;background:#101d30}}
@media(max-width:900px){{.grid2{{grid-template-columns:1fr}} h1{{font-size:32px}}}}
</style></head><body>
<header><div class="wrap">
<h1>CROSS-Neo Assay Feedback Loop</h1>
<p class="lead">Closed-loop data-driven neoantigen triage: model prior -> WT/decoy assay rows -> posterior readiness -> claim state -> next experiment.</p>
<div class="stats">
<div class="stat"><b>{summary['n_candidates']}</b><span>candidate rows</span></div>
<div class="stat"><b>{summary['n_assay_summaries']}</b><span>real assay summaries</span></div>
<div class="stat"><b>{summary['n_claim_states']}</b><span>claim states</span></div>
<div class="stat"><b>{summary['n_next_actions']}</b><span>next experiment rows</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Decision</h2>
<div class="ok">The loop is now operational. Add real rows to <code>{esc(DATA_IN)}</code>, rerun the script, and the posterior/claim/next-experiment tables update automatically.</div>
<div class="warn">Blank templates are never treated as labels. The posterior is a decision-support readiness score, not a clinical or universal immunogenicity probability.</div>
<div class="links">
<a class="pill" href="assets/cross_neo_md_audit/ASSAY_FEEDBACK_LEARNER_REPORT.md">Full report</a>
<a class="pill" href="assets/cross_neo_md_audit/ASSAY_FEEDBACK_ONE_PAGE_KR.md">Korean brief</a>
<a class="pill" href="assets/cross_neo_md_audit/assay_feedback_template.tsv">assay feedback template</a>
<a class="pill" href="cross_neo_assay_feedback_scenarios.html">what-if scenario simulator</a>
<a class="pill" href="cross_neo_physics_gate.html">physics gate</a>
<a class="pill" href="cross_neo_high_impact_decision.html">high-impact decision page</a>
</div></section>
<section><h2>Figures</h2><div class="grid2">{figures}</div></section>
<section><h2>Posterior Updates</h2>{table_html(post, ['row_id','peptide','hla_4digit','prior_cross_neo_i_readiness','posterior_immunogenicity_readiness','posterior_uncertainty','assay_rows_observed','claim_state','next_experiment','active_learning_utility'], 20)}</section>
<section><h2>Claim States</h2>{table_html(claim_df, ['row_id','peptide','hla_4digit','claim_state','allowed_claim_now','required_next_gate'], 20)}</section>
<section><h2>Next Experiment Queue</h2>{table_html(active, ['peptide','hla_4digit','next_experiment','active_learning_utility','utility_reason'], 20)}</section>
<section><h2>Feedback Rules</h2>{table_html(rules, ['observed_pattern','model_update','interpretation'], 20)}</section>
</main></body></html>"""
    PAGE.write_text(html_text)


def deploy(files: list[Path]) -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for path in files:
        if path.exists():
            safe_copy(path, ASSET / path.name)
            safe_copy(path, WEB_ASSET / path.name)
    for fig in [
        "fig_md56_assay_feedback_closed_loop",
        "fig_md57_prior_posterior_update",
        "fig_md58_active_learning_queue",
        "fig_md59_claim_state_distribution",
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
    (REPO / "project/data/cross_neo_assay_feedback").mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    schema = build_schema()
    template = build_template(data)
    assay_summary = summarize_assays(data["assay"])
    post, claim_df, active = build_posteriors(data, assay_summary)
    rules = build_rules()

    schema.to_csv(OUT / "assay_feedback_schema.tsv", sep="\t", index=False)
    template.to_csv(OUT / "assay_feedback_template.tsv", sep="\t", index=False)
    if not DATA_IN.exists():
        template.to_csv(DATA_IN, sep="\t", index=False)
    assay_summary.to_csv(OUT / "assay_feedback_summary.tsv", sep="\t", index=False)
    post.to_csv(OUT / "candidate_posterior_updates.tsv", sep="\t", index=False)
    claim_df.to_csv(OUT / "claim_state_updates.tsv", sep="\t", index=False)
    active.to_csv(OUT / "active_learning_next_experiments.tsv", sep="\t", index=False)
    rules.to_csv(OUT / "feedback_integration_rules.tsv", sep="\t", index=False)

    make_figures(post, active, claim_df)
    write_report(schema, template, assay_summary, post, claim_df, active, rules)

    summary = {
        "n_candidates": int(len(post)),
        "n_assay_summaries": int(len(assay_summary)),
        "n_claim_states": int(claim_df["claim_state"].nunique()) if not claim_df.empty else 0,
        "n_next_actions": int(len(active)),
        "assay_feedback_input": str(DATA_IN),
        "mode": "posterior_update_from_real_assay_rows" if not assay_summary.empty else "prior_only_template_ready",
        "figures": [
            "fig_md56_assay_feedback_closed_loop",
            "fig_md57_prior_posterior_update",
            "fig_md58_active_learning_queue",
            "fig_md59_claim_state_distribution",
        ],
        "boundary": "closed-loop assay learning, not immunogenicity proof without real WT/decoy assay labels",
    }
    (OUT / "assay_feedback_learner_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_html(post, claim_df, active, rules, summary)
    deploy(
        [
            OUT / "assay_feedback_schema.tsv",
            OUT / "assay_feedback_template.tsv",
            OUT / "assay_feedback_summary.tsv",
            OUT / "candidate_posterior_updates.tsv",
            OUT / "claim_state_updates.tsv",
            OUT / "active_learning_next_experiments.tsv",
            OUT / "feedback_integration_rules.tsv",
            OUT / "ASSAY_FEEDBACK_LEARNER_REPORT.md",
            OUT / "ASSAY_FEEDBACK_ONE_PAGE_KR.md",
            OUT / "assay_feedback_learner_summary.json",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[assay-feedback] wrote {OUT}")
    print(f"[assay-feedback] page {PAGE}")


if __name__ == "__main__":
    main()
