#!/usr/bin/env python3
"""Small-dataset uncertainty funnel for CROSS-Neo/TCR/MD candidates.

Goal: do not over-claim positives. Use ensembles, Bayesian shrinkage,
dropout/weight perturbation, TCR evidence, and MD/control evidence to remove
obviously weak or unstable candidates before expensive wetlab.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from scipy.stats import beta as beta_dist
except Exception:  # pragma: no cover
    beta_dist = None


REPO = Path(__file__).resolve().parents[3]
CROSS = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
OUT = MD_OUT / "small_dataset_uncertainty"
FIG = MD_OUT / "figures"

PREDICTION_FILES = [
    CROSS / "predictions/all_predictions.tsv",
    CROSS / "predictions/selective_ensemble_v2_1_predictions.tsv",
    CROSS / "predictions/fast_esm2_qk_gate_predictions.tsv",
    CROSS / "predictions/esm2_35m_fast_predictions.tsv",
    CROSS / "predictions/esm2_150m_fast_predictions.tsv",
    CROSS / "predictions/esm2_650m_fast_predictions.tsv",
]


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs) if path.exists() else pd.DataFrame()


def norm_hla(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    return s if s.startswith("HLA-") else f"HLA-{s}"


def add_key(df: pd.DataFrame, pep: str = "peptide", hla: str = "hla_4digit") -> pd.DataFrame:
    out = df.copy()
    out["peptide_norm"] = out[pep].astype(str).str.upper().str.strip()
    out["hla_norm"] = out[hla].map(norm_hla)
    out["pmhc_key"] = out["peptide_norm"] + "|" + out["hla_norm"]
    return out


def beta_interval(alpha: float, beta: float, q: float) -> float:
    if beta_dist is not None:
        return float(beta_dist.ppf(q, alpha, beta))
    mean = alpha / (alpha + beta)
    var = alpha * beta / (((alpha + beta) ** 2) * (alpha + beta + 1))
    z = 1.645 if q in {0.05, 0.95} else 1.96
    return float(np.clip(mean + (z if q > 0.5 else -z) * np.sqrt(var), 0, 1))


def load_candidate_frame() -> pd.DataFrame:
    wet = read_tsv(CROSS / "tcr_extension/tcr_wetlab_candidate_prioritization_unique_pmhc.tsv")
    if wet.empty:
        raise SystemExit("Missing wetlab candidate table")
    wet = add_key(wet)
    ultra = read_tsv(MD_OUT / "ultra_priority/ultra_wetlab_priority_candidates.tsv")
    if not ultra.empty:
        ultra = add_key(ultra)
        ucols = [
            "pmhc_key",
            "recommendation_tier",
            "ultra_priority_score",
            "confidence_score",
            "model_evidence_score",
            "tcr_resource_score",
            "md_structural_score",
            "control_readiness_score",
            "claim_risk_penalty",
            "md_label",
            "md_score",
            "live_completion_fraction",
            "tentative_wt_sequence",
            "wt_status",
            "anchor_preserved_decoy",
            "why",
        ]
        wet = wet.merge(ultra[[c for c in ucols if c in ultra.columns]], on="pmhc_key", how="left", suffixes=("", "_ultra"))

    registry = read_tsv(CROSS / "canonical_registry.tsv")
    if not registry.empty:
        registry = add_key(registry)
        reg_agg = registry.groupby("pmhc_key", as_index=False).agg(
            registry_row_count=("row_id", "nunique"),
            registry_positive_count=("label_binary", lambda s: int(pd.to_numeric(s, errors="coerce").fillna(0).sum())),
            registry_sources=("source_dataset", lambda s: "|".join(sorted(set(map(str, s.dropna()))))[:500]),
            public_overlap_flags_registry=("public_overlap_flags", lambda s: "|".join(sorted(set(map(str, s.dropna()))))[:500]),
        )
        wet = wet.merge(reg_agg, on="pmhc_key", how="left")
    return wet


def load_prediction_rows(candidate_rows: set[str]) -> pd.DataFrame:
    frames = []
    for path in PREDICTION_FILES:
        if not path.exists():
            continue
        usecols = None
        df = pd.read_csv(path, sep="\t", usecols=usecols)
        if "row_id" not in df.columns or "score" not in df.columns:
            continue
        df = df[df["row_id"].astype(str).isin(candidate_rows)].copy()
        if df.empty:
            continue
        df["prediction_file"] = path.name
        for col in ["split_name", "fold_id", "model_name", "model_family", "claim_status"]:
            if col not in df.columns:
                df[col] = ""
        df["score"] = pd.to_numeric(df["score"], errors="coerce").clip(0, 1)
        df = df[df["score"].notna()].copy()
        frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def aggregate_predictions(pred: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if pred.empty:
        return pd.DataFrame(), pd.DataFrame()
    pred["model_key"] = (
        pred["prediction_file"].astype(str)
        + "|"
        + pred["split_name"].astype(str)
        + "|"
        + pred["fold_id"].astype(str)
        + "|"
        + pred["model_name"].astype(str)
    )
    pred["family_key"] = pred["model_family"].fillna(pred["model_name"]).astype(str)
    agg = pred.groupby("row_id").agg(
        n_model_scores=("score", "count"),
        n_unique_models=("model_key", "nunique"),
        n_model_families=("family_key", "nunique"),
        ensemble_mean=("score", "mean"),
        ensemble_median=("score", "median"),
        ensemble_std=("score", "std"),
        ensemble_min=("score", "min"),
        ensemble_max=("score", "max"),
        ensemble_q05=("score", lambda s: float(np.quantile(s, 0.05))),
        ensemble_q25=("score", lambda s: float(np.quantile(s, 0.25))),
        ensemble_q75=("score", lambda s: float(np.quantile(s, 0.75))),
        ensemble_q95=("score", lambda s: float(np.quantile(s, 0.95))),
        vote_gt_050=("score", lambda s: float((s >= 0.50).mean())),
        vote_gt_080=("score", lambda s: float((s >= 0.80).mean())),
        rank_pct_median=("rank_pct", "median"),
        rank_pct_q25=("rank_pct", lambda s: float(np.nanquantile(pd.to_numeric(s, errors="coerce"), 0.25)) if pd.to_numeric(s, errors="coerce").notna().any() else np.nan),
    ).reset_index()
    agg["ensemble_std"] = agg["ensemble_std"].fillna(0)
    agg["ensemble_width_90"] = agg["ensemble_q95"] - agg["ensemble_q05"]
    return agg, pred


def bayesian_scores(pred: pd.DataFrame, cand: pd.DataFrame) -> pd.DataFrame:
    rows = []
    pred_group = {k: v for k, v in pred.groupby("row_id")} if not pred.empty else {}
    global_prev = float(pd.to_numeric(cand.get("label_binary", pd.Series([0.1])), errors="coerce").fillna(0).mean())
    global_prev = float(np.clip(global_prev, 0.02, 0.98))
    prior_strength = 4.0
    for row_id in cand["row_id"].astype(str).unique():
        scores = pred_group.get(row_id, pd.DataFrame()).get("score", pd.Series(dtype=float)).dropna().astype(float).clip(0, 1)
        n = len(scores)
        if n:
            # Effective sample size shrinks correlated model/fold predictions.
            eff_n = min(float(n), max(3.0, np.sqrt(n) * 2.0))
            scaled_sum = float(scores.mean() * eff_n)
        else:
            eff_n = 0.0
            scaled_sum = 0.0
        alpha = 1.0 + prior_strength * global_prev + scaled_sum
        beta = 1.0 + prior_strength * (1 - global_prev) + (eff_n - scaled_sum)
        rows.append(
            {
                "row_id": row_id,
                "bayes_effective_n": eff_n,
                "bayes_alpha": alpha,
                "bayes_beta": beta,
                "bayes_mean": alpha / (alpha + beta),
                "bayes_q05": beta_interval(alpha, beta, 0.05),
                "bayes_q25": beta_interval(alpha, beta, 0.25),
                "bayes_q75": beta_interval(alpha, beta, 0.75),
                "bayes_q95": beta_interval(alpha, beta, 0.95),
            }
        )
    out = pd.DataFrame(rows)
    out["bayes_width_90"] = out["bayes_q95"] - out["bayes_q05"]
    return out


def perturbation_scores(pred: pd.DataFrame, seed: int = 17) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    if pred.empty:
        return pd.DataFrame()
    rows = []
    dropout_rates = [0.0, 0.1, 0.25, 0.5, 0.75]
    for row_id, grp in pred.groupby("row_id"):
        scores = grp["score"].dropna().astype(float).clip(1e-4, 1 - 1e-4).to_numpy()
        if len(scores) == 0:
            continue
        samples = []
        by_dropout = {}
        for dr in dropout_rates:
            vals = []
            for _ in range(400):
                keep = rng.random(len(scores)) >= dr
                if not keep.any():
                    keep[rng.integers(0, len(scores))] = True
                s = scores[keep]
                # Dirichlet weight perturbation plus simple logit-temperature noise.
                weights = rng.dirichlet(np.ones(len(s)))
                logit = np.log(s / (1 - s))
                temp = rng.uniform(0.75, 1.35)
                noisy = 1 / (1 + np.exp(-(logit / temp + rng.normal(0, 0.18, size=len(s)))))
                vals.append(float(np.sum(weights * noisy)))
            by_dropout[dr] = np.array(vals)
            samples.extend(vals)
        samples = np.array(samples)
        row = {
            "row_id": row_id,
            "perturb_mean": float(samples.mean()),
            "perturb_median": float(np.median(samples)),
            "perturb_q05": float(np.quantile(samples, 0.05)),
            "perturb_q25": float(np.quantile(samples, 0.25)),
            "perturb_q75": float(np.quantile(samples, 0.75)),
            "perturb_q95": float(np.quantile(samples, 0.95)),
            "perturb_width_90": float(np.quantile(samples, 0.95) - np.quantile(samples, 0.05)),
            "perturb_prob_gt_050": float((samples >= 0.50).mean()),
            "perturb_prob_gt_080": float((samples >= 0.80).mean()),
        }
        for dr, vals in by_dropout.items():
            row[f"dropout_{int(dr*100):02d}_median"] = float(np.median(vals))
            row[f"dropout_{int(dr*100):02d}_width90"] = float(np.quantile(vals, 0.95) - np.quantile(vals, 0.05))
        row["dropout_sensitivity"] = max(row[f"dropout_{int(dr*100):02d}_median"] for dr in dropout_rates) - min(
            row[f"dropout_{int(dr*100):02d}_median"] for dr in dropout_rates
        )
        rows.append(row)
    return pd.DataFrame(rows)


def decide(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in [
        "ensemble_mean",
        "ensemble_q05",
        "ensemble_q95",
        "vote_gt_050",
        "bayes_mean",
        "bayes_q05",
        "bayes_q95",
        "perturb_median",
        "perturb_q05",
        "perturb_q95",
        "perturb_prob_gt_050",
        "perturb_width_90",
        "dropout_sensitivity",
        "tcr_evidence_count",
        "paired_tcr_evidence_count",
        "md_structural_score",
        "control_readiness_score",
        "ultra_priority_score",
        "claim_risk_penalty",
    ]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    out["has_tcr_evidence"] = out.get("tcr_evidence_count", 0).fillna(0) > 0
    out["has_paired_tcr"] = out.get("paired_tcr_evidence_count", 0).fillna(0) > 0
    out["has_md_support"] = out.get("md_structural_score", 0).fillna(0) >= 0.45
    out["active_md_pending"] = out.get("recommendation_tier", "").fillna("").astype(str).str.contains("PENDING_MD", na=False)
    out["unstable_model"] = (out.get("perturb_width_90", 0).fillna(0) >= 0.45) | (out.get("dropout_sensitivity", 0).fillna(0) >= 0.18)
    out["low_posterior"] = out.get("bayes_q95", 0).fillna(0) < 0.40
    out["weak_lower_bound"] = out.get("bayes_q05", 0).fillna(0) < 0.25
    out["high_model_support"] = (out.get("bayes_mean", 0).fillna(0) >= 0.55) & (out.get("perturb_prob_gt_050", 0).fillna(0) >= 0.70)

    decisions = []
    reasons = []
    for _, r in out.iterrows():
        tier = str(r.get("recommendation_tier", ""))
        if "EXPERIMENT_NOW" in tier:
            decisions.append("SURVIVE_EXPERIMENT_NOW_WITH_CONTROLS")
            reasons.append("tier-A integrated TCR/MD/control evidence; still not immunogenicity proof")
        elif "PENDING_MD" in tier:
            decisions.append("SURVIVE_PENDING_MD_COMPLETION")
            reasons.append("strong TCR evidence and active MD; wait for full trajectory before final wetlab order")
        elif bool(r.get("low_posterior")) and not bool(r.get("has_tcr_evidence")) and not bool(r.get("has_md_support")):
            decisions.append("CULL_LOW_POSTERIOR_NO_TCR_NO_MD")
            reasons.append("Bayesian upper bound low and no independent TCR/MD support")
        elif bool(r.get("unstable_model")) and not bool(r.get("has_paired_tcr")) and not bool(r.get("has_md_support")):
            decisions.append("CULL_DROPOUT_UNSTABLE_NO_ORTHOGONAL_SUPPORT")
            reasons.append("model probability changes strongly under perturbation/dropout and lacks paired TCR/MD support")
        elif bool(r.get("high_model_support")) and bool(r.get("has_tcr_evidence")):
            decisions.append("REVIEW_FOR_TCR_WT_CURATION")
            reasons.append("model/TCR signal survives but needs WT, structure, or paired TCR curation")
        elif bool(r.get("has_tcr_evidence")) or float(r.get("ultra_priority_score", 0) or 0) >= 0.15:
            decisions.append("HOLD_FOR_MORE_EVIDENCE")
            reasons.append("some evidence exists but not enough for wetlab escalation")
        else:
            decisions.append("CULL_LOW_ACTIONABILITY")
            reasons.append("insufficient model/TCR/MD/control evidence for this small dataset")
    out["funnel_decision"] = decisions
    out["funnel_reason"] = reasons
    out["funnel_priority_score"] = (
        0.30 * out.get("bayes_mean", 0).fillna(0)
        + 0.20 * out.get("perturb_prob_gt_050", 0).fillna(0)
        + 0.20 * out.get("md_structural_score", 0).fillna(0)
        + 0.20 * out.get("tcr_resource_score", 0).fillna(0)
        + 0.10 * out.get("control_readiness_score", 0).fillna(0)
        - 0.15 * out.get("claim_risk_penalty", 0).fillna(0)
    ).clip(0, 1)
    order = {
        "SURVIVE_EXPERIMENT_NOW_WITH_CONTROLS": 0,
        "SURVIVE_PENDING_MD_COMPLETION": 1,
        "REVIEW_FOR_TCR_WT_CURATION": 2,
        "HOLD_FOR_MORE_EVIDENCE": 3,
        "CULL_DROPOUT_UNSTABLE_NO_ORTHOGONAL_SUPPORT": 4,
        "CULL_LOW_POSTERIOR_NO_TCR_NO_MD": 5,
        "CULL_LOW_ACTIONABILITY": 6,
    }
    out["funnel_order"] = out["funnel_decision"].map(order).fillna(9).astype(int)
    return out.sort_values(["funnel_order", "funnel_priority_score"], ascending=[True, False])


def plot_outputs(final: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    counts = final["funnel_decision"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(10.8, 5.6))
    ax.barh(counts.index, counts.values, color="#4d96ff")
    ax.set_xlabel("Candidate count")
    ax.set_title("Small-dataset uncertainty funnel: keep only actionable candidates")
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(FIG / "fig_md16_small_dataset_uncertainty_funnel.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md16_small_dataset_uncertainty_funnel.pdf", bbox_inches="tight")
    plt.close(fig)

    show = final.head(20).copy()
    label = show["peptide"].astype(str) + "\n" + show["hla_4digit"].astype(str)
    fig, ax = plt.subplots(figsize=(11.5, 6))
    ax.errorbar(
        show["bayes_mean"],
        label,
        xerr=[show["bayes_mean"] - show["bayes_q05"], show["bayes_q95"] - show["bayes_mean"]],
        fmt="o",
        color="#2fbf71",
        ecolor="#8fd8ad",
        label="Bayesian posterior 90%",
    )
    ax.scatter(show["perturb_median"], label, color="#e76f51", marker="x", label="dropout/weight perturb median")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Score")
    ax.set_title("Bayesian and dropout uncertainty for highest-action candidates")
    ax.legend(frameon=False)
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(FIG / "fig_md17_bayesian_dropout_uncertainty.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md17_bayesian_dropout_uncertainty.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(final: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cols = [
        "row_id",
        "peptide",
        "hla_4digit",
        "funnel_decision",
        "funnel_priority_score",
        "funnel_reason",
        "bayes_mean",
        "bayes_q05",
        "bayes_q95",
        "perturb_median",
        "perturb_q05",
        "perturb_q95",
        "perturb_width_90",
        "dropout_sensitivity",
        "ensemble_mean",
        "ensemble_std",
        "vote_gt_050",
        "recommendation_tier",
        "md_label",
        "md_structural_score",
        "tcr_evidence_count",
        "paired_tcr_evidence_count",
        "control_readiness_score",
        "claim_risk_penalty",
    ]
    final[[c for c in cols if c in final.columns]].to_csv(OUT / "small_dataset_uncertainty_funnel.tsv", sep="\t", index=False)
    final[final["funnel_decision"].str.startswith("SURVIVE") | final["funnel_decision"].str.startswith("REVIEW")].to_csv(
        OUT / "small_dataset_wetlab_survivors.tsv", sep="\t", index=False
    )
    final[final["funnel_decision"].str.startswith("CULL")].to_csv(OUT / "small_dataset_culled_candidates.tsv", sep="\t", index=False)
    counts = final["funnel_decision"].value_counts().rename_axis("decision").reset_index(name="n")
    counts.to_csv(OUT / "culling_decision_counts.tsv", sep="\t", index=False)
    top = final[[c for c in cols if c in final.columns]].head(20)
    lines = [
        "# Small-Dataset Bayesian/Dropout Uncertainty Funnel",
        "",
        "## Philosophy",
        "",
        "This is a culling strategy, not a positive-calling oracle. It keeps candidates only when model evidence survives perturbation or when independent TCR/MD evidence supports further work.",
        "",
        "## Decision Counts",
        "",
        counts.to_markdown(index=False),
        "",
        "## Top Actionable Candidates",
        "",
        top.to_markdown(index=False),
        "",
        "## Rules",
        "",
        "- `SURVIVE_EXPERIMENT_NOW_WITH_CONTROLS`: integrated TCR/MD/control evidence supports wetlab prioritization.",
        "- `SURVIVE_PENDING_MD_COMPLETION`: TCR/model evidence is strong but full MD must finish.",
        "- `REVIEW_FOR_TCR_WT_CURATION`: model/TCR survives enough to curate WT, paired TCR, or structure before wetlab.",
        "- `CULL_*`: not discarded biologically forever, but removed from expensive wetlab priority for this small dataset.",
        "",
        "## Claim Boundary",
        "",
        "Bayesian posterior and dropout perturbation scores are ranking/uncertainty summaries. They are not calibrated immunogenicity probabilities unless validated on strict external wetlab datasets.",
    ]
    (OUT / "small_dataset_uncertainty_funnel_report.md").write_text("\n".join(lines) + "\n")
    (OUT / "small_dataset_uncertainty_summary.json").write_text(
        json.dumps(
            {
                "n_candidates": int(len(final)),
                "decision_counts": counts.set_index("decision")["n"].to_dict(),
                "top_rows": final[["row_id", "peptide", "hla_4digit", "funnel_decision", "funnel_priority_score"]].head(5).to_dict("records"),
            },
            indent=2,
        )
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cand = load_candidate_frame()
    pred = load_prediction_rows(set(cand["row_id"].astype(str)))
    pred_agg, pred_long = aggregate_predictions(pred)
    bayes = bayesian_scores(pred_long, cand)
    perturb = perturbation_scores(pred_long)
    pred_agg.to_csv(OUT / "model_score_aggregation.tsv", sep="\t", index=False)
    bayes.to_csv(OUT / "bayesian_posterior_scores.tsv", sep="\t", index=False)
    perturb.to_csv(OUT / "dropout_perturbation_summary.tsv", sep="\t", index=False)

    final = cand.merge(pred_agg, on="row_id", how="left").merge(bayes, on="row_id", how="left").merge(perturb, on="row_id", how="left")
    final = decide(final)
    write_report(final)
    plot_outputs(final)
    print("[small-dataset-funnel]", (OUT / "small_dataset_uncertainty_summary.json").read_text())


if __name__ == "__main__":
    main()
