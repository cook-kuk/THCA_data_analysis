#!/usr/bin/env python3
"""Search threshold presets for the DL-first neoantigen funnel.

This turns the slider dashboard into decision intelligence:
find operating points that trade off precision, recall, false positives,
and expensive MD/wetlab escalation size using the available labels.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
SCRIPT_DIR = REPO / "project/scripts/cross_neo_md"
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
OUT = MD_OUT / "dl_threshold_optimization"
FIG = MD_OUT / "figures"

CALL_RULES = ["md_escalation", "wetlab_shortlist", "structure_md_supported", "tcr_supported", "robust_dl", "non_culled"]


def load_ui_module():
    path = SCRIPT_DIR / "24_build_dl_funnel_slider_dashboard.py"
    spec = importlib.util.spec_from_file_location("dl_slider_ui", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def prepare_frame() -> pd.DataFrame:
    mod = load_ui_module()
    dl = mod.load_dl_module()
    df = mod.add_baker_scores(dl.score_and_gate(dl.load_inputs()))
    return df


def arr(df: pd.DataFrame, col: str, default: float = 0.0) -> np.ndarray:
    if col not in df.columns:
        return np.full(len(df), default, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default).to_numpy(float)


def bool_arr(df: pd.DataFrame, col: str) -> np.ndarray:
    if col not in df.columns:
        return np.zeros(len(df), dtype=bool)
    s = df[col]
    if s.dtype == bool:
        return s.to_numpy(bool)
    return s.fillna("").astype(str).str.lower().isin(["true", "1", "yes"]).to_numpy(bool)


def prepare_arrays(df: pd.DataFrame) -> dict[str, np.ndarray]:
    hla = df["hla_4digit"].fillna("").astype(str).str.upper()
    md_score = np.maximum(arr(df, "md_structural_score"), arr(df, "md_score"))
    return {
        "label": arr(df, "label_binary", 0).astype(int),
        "pep_len": arr(df, "peptide_length", 0),
        "is_class_i": hla.str.match(r"^HLA-[ABC]\*").to_numpy(bool),
        "dropout_sensitivity": arr(df, "dropout_sensitivity"),
        "perturb_width_90": arr(df, "perturb_width_90"),
        "main_dl_score": arr(df, "main_dl_score"),
        "ensemble_q95": arr(df, "ensemble_q95"),
        "pmhc_score_mean": arr(df, "pmhc_score_mean"),
        "bayes_mean": arr(df, "bayes_mean"),
        "perturb_prob_gt_050": arr(df, "perturb_prob_gt_050"),
        "bayes_q95": arr(df, "bayes_q95"),
        "tcr_evidence_count": arr(df, "tcr_evidence_count"),
        "cancer_context_evidence_count": arr(df, "cancer_context_evidence_count"),
        "pathogen_context_evidence_count": arr(df, "pathogen_context_evidence_count"),
        "paired_tcr_evidence_count": arr(df, "paired_tcr_evidence_count"),
        "tcr_augmented_score_mean": arr(df, "tcr_augmented_score_mean"),
        "best_tcr_augmented_score": arr(df, "best_tcr_augmented_score"),
        "baker_structural_score": arr(df, "baker_structural_score"),
        "md_score": md_score,
        "live_completion_fraction": arr(df, "live_completion_fraction"),
        "wt_or_decoy_ready": bool_arr(df, "wt_or_decoy_ready"),
    }


def evaluate(data: dict[str, np.ndarray], cfg: dict[str, float], call_rule: str) -> dict[str, float | int | str]:
    label = data["label"]
    pep_len = data["pep_len"]
    is_class_i = data["is_class_i"]
    supp = np.where(
        is_class_i,
        (pep_len >= cfg["class_i_min_len"]) & (pep_len <= cfg["class_i_max_len"]),
        (pep_len >= cfg["other_min_len"]) & (pep_len <= cfg["other_max_len"]),
    )
    dropout_unstable = (data["dropout_sensitivity"] >= cfg["dropout_sens"]) | (data["perturb_width_90"] >= cfg["perturb_width"])
    main_dl = data["main_dl_score"]
    stage1 = supp & ((main_dl >= cfg["main_dl"]) | (data["ensemble_q95"] >= cfg["ensemble_upper"]) | (data["pmhc_score_mean"] >= cfg["pmhc"]))
    stage2 = stage1 & (
        ((data["bayes_mean"] >= cfg["bayes_mean"]) & (data["perturb_prob_gt_050"] >= cfg["perturb_prob"]) & (~dropout_unstable))
        | (data["bayes_q95"] >= cfg["bayes_upper"])
        | (main_dl >= cfg["robust_main"])
    )
    stage3 = stage2 & (~dropout_unstable) & (
        (data["bayes_mean"] >= cfg["robust_bayes"])
        | (data["perturb_prob_gt_050"] >= cfg["robust_perturb"])
        | (main_dl >= cfg["robust_main"])
    )
    source_tcr = (
        (data["tcr_evidence_count"] >= cfg["min_tcr_evidence"])
        & (data["cancer_context_evidence_count"] >= cfg["min_cancer_context"])
        & (data["pathogen_context_evidence_count"] <= cfg["max_pathogen_context"])
    )
    paired_tcr = source_tcr & (data["paired_tcr_evidence_count"] >= cfg["min_paired_tcr"])
    tcr_support = (data["tcr_augmented_score_mean"] >= cfg["tcr_branch"]) | (data["best_tcr_augmented_score"] >= cfg["tcr_rescue"])
    tcr_rescue = (
        supp
        & paired_tcr
        & (data["tcr_augmented_score_mean"] >= cfg["tcr_rescue"])
        & (data["best_tcr_augmented_score"] >= cfg["tcr_rescue"])
        & (data["bayes_q95"] >= cfg["tcr_rescue_bayes_upper"])
    )
    baker_pass = data["baker_structural_score"] >= cfg["min_baker_structural"]
    md_pass = (
        (data["md_score"] >= cfg["min_md_structural"])
        | (data["live_completion_fraction"] >= cfg["min_live_completion"])
    )

    primary_paired = stage3 & paired_tcr & tcr_support
    tcr_context = stage3 & source_tcr & tcr_support
    discordant = stage3 & source_tcr & (~tcr_support)
    cull = (~supp) | (~stage1) | (stage1 & (~stage2) & (~tcr_rescue)) | (dropout_unstable & (~paired_tcr))

    if call_rule == "md_escalation":
        pred = primary_paired | tcr_rescue
    elif call_rule == "wetlab_shortlist":
        pred = primary_paired & data["wt_or_decoy_ready"]
    elif call_rule == "structure_md_supported":
        pred = stage3 & source_tcr & tcr_support & baker_pass & md_pass
    elif call_rule == "tcr_supported":
        pred = tcr_context
    elif call_rule == "robust_dl":
        pred = stage3
    elif call_rule == "non_culled":
        pred = (~cull) & (~discordant)
    else:
        pred = primary_paired | tcr_rescue

    tp = int((pred & (label == 1)).sum())
    fp = int((pred & (label == 0)).sum())
    tn = int(((~pred) & (label == 0)).sum())
    fn = int(((~pred) & (label == 1)).sum())
    called = int(pred.sum())
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    specificity = tn / (tn + fp) if tn + fp else 0.0
    fpr = fp / (fp + tn) if fp + tn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    enrichment = precision / max(float((label == 1).mean()), 1e-9) if called else 0.0

    row: dict[str, float | int | str] = {
        "call_rule": call_rule,
        "called_positive": called,
        "TP": tp,
        "TN": tn,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "FPR": fpr,
        "F1": f1,
        "enrichment_over_prevalence": enrichment,
        "stage1_n": int(stage1.sum()),
        "stage2_n": int(stage2.sum()),
        "stage3_n": int(stage3.sum()),
        "source_tcr_n": int((stage3 & source_tcr).sum()),
        "paired_tcr_n": int((stage3 & paired_tcr).sum()),
        "tcr_supported_n": int((stage3 & source_tcr & tcr_support).sum()),
        "baker_supported_n": int((stage3 & source_tcr & tcr_support & baker_pass).sum()),
        "md_supported_n": int((stage3 & source_tcr & tcr_support & md_pass).sum()),
    }
    row.update(cfg)
    return row


def sample_configs(n: int, seed: int = 19) -> list[dict[str, float]]:
    rng = np.random.default_rng(seed)
    base_mod = load_ui_module()
    base = dict(base_mod.DEFAULTS)
    configs = [base]
    for _ in range(n - 1):
        class_min = int(rng.choice([8, 8, 8, 9]))
        class_max = int(rng.choice([9, 10, 11, 11, 12, 14]))
        if class_max < class_min:
            class_max = class_min
        other_min = int(rng.choice([8, 9, 10, 12]))
        other_max = int(rng.choice([15, 20, 25, 30, 35]))
        if other_max < other_min:
            other_max = other_min
        cfg = {
            "class_i_min_len": class_min,
            "class_i_max_len": class_max,
            "other_min_len": other_min,
            "other_max_len": other_max,
            "main_dl": float(rng.uniform(0.25, 0.85)),
            "ensemble_upper": float(rng.uniform(0.45, 0.98)),
            "pmhc": float(rng.uniform(0.30, 0.90)),
            "bayes_mean": float(rng.uniform(0.25, 0.75)),
            "perturb_prob": float(rng.uniform(0.05, 0.95)),
            "bayes_upper": float(rng.uniform(0.40, 0.90)),
            "robust_bayes": float(rng.uniform(0.35, 0.80)),
            "robust_perturb": float(rng.uniform(0.30, 0.98)),
            "robust_main": float(rng.uniform(0.45, 0.95)),
            "dropout_sens": float(rng.uniform(0.02, 0.30)),
            "perturb_width": float(rng.uniform(0.10, 0.70)),
            "min_tcr_evidence": int(rng.choice([0, 1, 1, 2, 5, 10, 20, 50])),
            "min_paired_tcr": int(rng.choice([0, 1, 1, 2, 3, 5, 10])),
            "min_cancer_context": int(rng.choice([0, 1, 1, 2, 5, 10, 20, 50])),
            "max_pathogen_context": int(rng.choice([0, 0, 1, 5, 20, 100, 500])),
            "tcr_branch": float(rng.uniform(0.20, 0.99)),
            "tcr_rescue": float(rng.uniform(0.40, 0.99)),
            "tcr_rescue_bayes_upper": float(rng.uniform(0.30, 0.80)),
            "min_baker_structural": float(rng.uniform(0.00, 0.95)),
            "min_md_structural": float(rng.uniform(0.00, 0.90)),
            "min_live_completion": float(rng.uniform(0.00, 0.99)),
        }
        configs.append(cfg)
    return configs


def pareto_frontier(df: pd.DataFrame) -> pd.DataFrame:
    # Higher TP, higher precision, lower FP, and fewer called candidates are desirable.
    cand = df[df["called_positive"] > 0].copy()
    if cand.empty:
        return cand
    cand = cand.sort_values(["FP", "called_positive", "precision", "TP"], ascending=[True, True, False, False])
    frontier = []
    best_tp = -1
    best_precision = -1.0
    for _, r in cand.iterrows():
        if r["TP"] > best_tp or r["precision"] > best_precision:
            frontier.append(r)
            best_tp = max(best_tp, r["TP"])
            best_precision = max(best_precision, r["precision"])
    return pd.DataFrame(frontier)


def choose_presets(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    nonzero = results[results["called_positive"] > 0].copy()
    if nonzero.empty:
        return pd.DataFrame()

    def pick(name: str, subset: pd.DataFrame, sort_cols: list[str], asc: list[bool]) -> None:
        if subset.empty:
            return
        r = subset.sort_values(sort_cols, ascending=asc).iloc[0].copy()
        r["preset_name"] = name
        rows.append(r)

    pick("NO_FALSE_POSITIVE_MAX_TP", nonzero[nonzero["FP"] == 0], ["TP", "called_positive", "recall"], [False, True, False])
    pick("HIGH_PRECISION_MIN_FP", nonzero[nonzero["called_positive"] <= 20], ["precision", "TP", "FP", "called_positive"], [False, False, True, True])
    pick("BALANCED_F1", nonzero, ["F1", "precision", "recall"], [False, False, False])
    pick("RECALL_PRESERVING", nonzero[nonzero["FP"] <= 50], ["recall", "precision", "FP"], [False, False, True])
    pick("WETLAB_ULTRA_STRICT", nonzero[nonzero["call_rule"] == "wetlab_shortlist"], ["precision", "TP", "FP"], [False, False, True])
    pick("STRUCTURE_MD_STRICT", nonzero[nonzero["call_rule"] == "structure_md_supported"], ["precision", "TP", "FP"], [False, False, True])
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.drop_duplicates(subset=["preset_name"]).reset_index(drop=True)


def plot_results(results: pd.DataFrame, presets: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    nonzero = results[results["called_positive"] > 0].copy()
    if nonzero.empty:
        return
    fig, ax = plt.subplots(figsize=(8.4, 6.2))
    for rule, sub in nonzero.groupby("call_rule"):
        ax.scatter(sub["recall"], sub["precision"], s=np.clip(sub["called_positive"], 5, 80), alpha=0.35, label=rule)
    if not presets.empty:
        ax.scatter(presets["recall"], presets["precision"], s=160, marker="*", color="#f2c46d", edgecolor="#111", label="recommended presets")
        for _, r in presets.iterrows():
            ax.annotate(str(r["preset_name"]).replace("_", "\n"), (r["recall"], r["precision"]), fontsize=7)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Threshold search: precision-recall operating points")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=7)
    fig.savefig(FIG / "fig_md22_threshold_precision_recall_frontier.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md22_threshold_precision_recall_frontier.pdf", bbox_inches="tight")
    plt.close(fig)

    if not presets.empty:
        labels = presets["preset_name"].astype(str).str.replace("_", "\n")
        x = np.arange(len(presets))
        fig, ax = plt.subplots(figsize=(10.5, 5.6))
        ax.bar(x - 0.3, presets["TP"], width=0.2, label="TP", color="#2fbf71")
        ax.bar(x - 0.1, presets["FP"], width=0.2, label="FP", color="#ff7b72")
        ax.bar(x + 0.1, presets["FN"], width=0.2, label="FN", color="#ffcc66")
        ax.bar(x + 0.3, presets["called_positive"], width=0.2, label="called", color="#4d96ff")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("Count")
        ax.set_title("Recommended threshold presets: confusion counts")
        ax.legend(frameon=False)
        ax.grid(axis="y", alpha=0.25)
        fig.savefig(FIG / "fig_md23_threshold_preset_confusion.png", dpi=220, bbox_inches="tight")
        fig.savefig(FIG / "fig_md23_threshold_preset_confusion.pdf", bbox_inches="tight")
        plt.close(fig)


def write_outputs(results: pd.DataFrame, frontier: pd.DataFrame, presets: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    results.to_csv(OUT / "threshold_random_search_results.tsv", sep="\t", index=False)
    frontier.to_csv(OUT / "threshold_pareto_frontier.tsv", sep="\t", index=False)
    presets.to_csv(OUT / "threshold_recommended_presets.tsv", sep="\t", index=False)
    threshold_keys = list(load_ui_module().DEFAULTS.keys())
    preset_payload = []
    for _, r in presets.iterrows():
        thresholds = {k: float(r[k]) for k in threshold_keys if k in r}
        for k in ["class_i_min_len", "class_i_max_len", "other_min_len", "other_max_len", "min_tcr_evidence", "min_paired_tcr", "min_cancer_context", "max_pathogen_context"]:
            if k in thresholds:
                thresholds[k] = int(round(thresholds[k]))
        preset_payload.append(
            {
                "name": r["preset_name"],
                "call_rule": r["call_rule"],
                "thresholds": thresholds,
                "metrics": {
                    "called_positive": int(r["called_positive"]),
                    "TP": int(r["TP"]),
                    "TN": int(r["TN"]),
                    "FP": int(r["FP"]),
                    "FN": int(r["FN"]),
                    "precision": float(r["precision"]),
                    "recall": float(r["recall"]),
                    "F1": float(r["F1"]),
                    "specificity": float(r["specificity"]),
                    "FPR": float(r["FPR"]),
                },
            }
        )
    (OUT / "threshold_recommended_presets.json").write_text(json.dumps(preset_payload, indent=2))
    lines = [
        "# DL Funnel Threshold Optimization",
        "",
        "## Purpose",
        "",
        "This searches slider thresholds against the available labels to create operating presets instead of stopping at manual threshold tuning.",
        "",
        "## Recommended Presets",
        "",
        presets[
            [
                "preset_name",
                "call_rule",
                "called_positive",
                "TP",
                "TN",
                "FP",
                "FN",
                "precision",
                "recall",
                "F1",
                "specificity",
                "FPR",
            ]
        ].to_markdown(index=False)
        if not presets.empty
        else "No presets found.",
        "",
        "## Claim Boundary",
        "",
        "These presets are optimized on the current labeled candidate table and must be treated as decision-support settings, not externally validated clinical thresholds.",
    ]
    (OUT / "threshold_optimization_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = prepare_frame()
    data = prepare_arrays(df)
    configs = sample_configs(20000)
    rows = []
    for cfg in configs:
        for rule in CALL_RULES:
            rows.append(evaluate(data, cfg, rule))
    results = pd.DataFrame(rows)
    frontier = pareto_frontier(results)
    presets = choose_presets(results)
    plot_results(results, presets)
    write_outputs(results, frontier, presets)
    print("[threshold-opt]", json.dumps({
        "n_results": int(len(results)),
        "n_frontier": int(len(frontier)),
        "presets": presets[["preset_name", "call_rule", "called_positive", "TP", "FP", "precision", "recall", "F1"]].to_dict("records") if not presets.empty else [],
    }, indent=2))


if __name__ == "__main__":
    main()
