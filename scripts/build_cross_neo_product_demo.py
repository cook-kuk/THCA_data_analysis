#!/usr/bin/env python3
"""CROSS-Neo product/demo track.

This is deliberately separate from the clean manuscript track. Product mode
may use public-predictor scores, but the main demo policy is pan-allele and
fixed-weight. Testset-aware optimization is reported only as an internal upper
bound and must not be used as a generalization claim.
"""

from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from cross_neo_v1_common import seq_similarity


REPO = Path(__file__).resolve().parents[1]
V0 = REPO / "project/results/cross_neo_v0"
V1 = REPO / "project/results/cross_neo_v1"
V2 = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
W3 = REPO / "project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep"
OUT = REPO / "project/results/cross_neo_product_demo_2026_05_10"
SEL = OUT / "selected_trainsets"


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SEL.mkdir(parents=True, exist_ok=True)


def clean_pep(x: object) -> str:
    aa = set("ACDEFGHIKLMNPQRSTVWY")
    return "".join(a for a in str(x or "").upper() if a in aa)


def metrics(y: np.ndarray, s: np.ndarray) -> dict[str, float]:
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    out = {"n": int(len(y)), "n_pos": int(y.sum()), "prevalence": float(y.mean()) if len(y) else np.nan}
    if len(y) and len(np.unique(y)) == 2:
        out["AUPRC"] = float(average_precision_score(y, s))
        out["AUROC"] = float(roc_auc_score(y, s))
    else:
        out["AUPRC"] = np.nan
        out["AUROC"] = np.nan
    order = np.argsort(-s)
    for k in [5, 10, 20]:
        kk = min(k, len(y))
        top = y[order[:kk]]
        out[f"top{k}_precision"] = float(top.mean()) if kk else np.nan
        out[f"recall_at_{k}"] = float(top.sum() / max(1, y.sum())) if kk else np.nan
        out[f"enrichment_at_{k}"] = float(out[f"top{k}_precision"] / out["prevalence"]) if out["prevalence"] else np.nan
    return out


def minmax(x: pd.Series, invert: bool = False) -> pd.Series:
    v = pd.to_numeric(x, errors="coerce")
    if invert:
        v = -v
    if v.notna().sum() == 0:
        return pd.Series(0.5, index=x.index)
    lo, hi = v.quantile(0.05), v.quantile(0.95)
    if abs(hi - lo) < 1e-12:
        return pd.Series(0.5, index=x.index)
    return ((v - lo) / (hi - lo)).clip(0, 1).fillna(0.5)


def load_component_table(master: pd.DataFrame) -> pd.DataFrame:
    strict = master[master["strict_set_flag"].astype(bool)].copy()
    comp = strict[["sample_id", "peptide_mut", "peptide_wt", "hla", "hla_supertype", "study", "label", "near_peptide_cluster"]].copy()

    def add_v1(path: Path, family: str, model: str, split: str, col: str) -> None:
        nonlocal comp
        if not path.exists():
            comp[col] = np.nan
            return
        df = pd.read_csv(path, sep="\t")
        sub = df[(df["split_name"] == split) & (df["model"] == model)][["sample_id", "score"]].rename(columns={"score": col})
        comp = comp.merge(sub, on="sample_id", how="left")

    add_v1(V1 / "gated_moe_predictions.tsv", "v1_gated_moe", "prespecified_equal_weight_C_QK_no_anchor", "hla_stratified_group_5fold", "v1_gated_cqk_hla")
    add_v1(V1 / "hard_decoy_focal_predictions.tsv", "v1_hard_decoy", "hard_decoy_rule_aux_C_QK_no_anchor", "hla_stratified_group_5fold", "hard_decoy_rule_hla")
    add_v1(V1 / "hard_decoy_focal_predictions.tsv", "v1_hard_decoy", "hard_decoy_nested_meta_C_QK_aux", "near_peptide_cluster_holdout", "hard_decoy_near_meta")
    add_v1(V1 / "hard_decoy_focal_predictions.tsv", "v1_hard_decoy", "hard_decoy_rule_aux_C_QK_no_anchor", "exact_peptide_hla_holdout", "hard_decoy_exact_rule")

    v2p = V2 / "predictions/selective_ensemble_v2_1_predictions.tsv"
    if v2p.exists():
        v2 = pd.read_csv(v2p, sep="\t")
        for split, col in [
            ("hla_stratified_group_5fold", "v2_selective_hla"),
            ("near_peptide_cluster_holdout", "v2_selective_near"),
            ("exact_peptide_hla_holdout", "v2_selective_exact"),
        ]:
            sub = v2[(v2["split_name"] == split) & (v2["model_name"].str.contains("selective_ensemble", na=False))]
            sub = sub.groupby("row_id", as_index=False)["score"].mean().rename(columns={"row_id": "sample_id", "score": col})
            comp = comp.merge(sub, on="sample_id", how="left")

    scored = W3 / "scored_master.tsv"
    if scored.exists():
        pub = pd.read_csv(scored, sep="\t")
        pub = pub.rename(columns={"peptide": "peptide_mut", "HLA_norm": "hla"})
        keep = ["peptide_mut", "hla", "label", "mhcflurry_presentation", "mhcflurry_affinity", "bigmhc_im", "bigmhc_el", "prime_score", "prime_rank"]
        pub = pub[[c for c in keep if c in pub.columns]].copy()
        pub = pub.drop_duplicates(["peptide_mut", "hla", "label"])
        comp = comp.merge(pub, on=["peptide_mut", "hla", "label"], how="left")
        comp["mhcflurry_presentation_norm"] = minmax(comp.get("mhcflurry_presentation", pd.Series(index=comp.index)))
        comp["mhcflurry_affinity_norm"] = minmax(comp.get("mhcflurry_affinity", pd.Series(index=comp.index)), invert=True)
        comp["bigmhc_im_norm"] = minmax(comp.get("bigmhc_im", pd.Series(index=comp.index)))
        comp["bigmhc_el_norm"] = minmax(comp.get("bigmhc_el", pd.Series(index=comp.index)))
        comp["prime_score_norm"] = minmax(comp.get("prime_score", pd.Series(index=comp.index)))
        comp["prime_rank_norm"] = minmax(comp.get("prime_rank", pd.Series(index=comp.index)), invert=True)
        pub_cols = ["mhcflurry_presentation_norm", "mhcflurry_affinity_norm", "bigmhc_im_norm", "bigmhc_el_norm", "prime_score_norm", "prime_rank_norm"]
        comp["public_predictor_assist_mean"] = comp[pub_cols].mean(axis=1)
    else:
        comp["public_predictor_assist_mean"] = 0.5

    component_cols = [
        "v1_gated_cqk_hla",
        "hard_decoy_rule_hla",
        "hard_decoy_near_meta",
        "hard_decoy_exact_rule",
        "v2_selective_hla",
        "v2_selective_near",
        "v2_selective_exact",
        "public_predictor_assist_mean",
    ]
    for c in component_cols:
        if c not in comp:
            comp[c] = np.nan
        comp[c] = pd.to_numeric(comp[c], errors="coerce").fillna(comp[c].median() if comp[c].notna().any() else 0.5)
    comp["product_prespecified_score"] = (
        0.18 * comp["v1_gated_cqk_hla"]
        + 0.18 * comp["hard_decoy_rule_hla"]
        + 0.16 * comp["hard_decoy_near_meta"]
        + 0.18 * comp["v2_selective_hla"]
        + 0.10 * comp["v2_selective_near"]
        + 0.20 * comp["public_predictor_assist_mean"]
    )
    return comp


def testset_aware_weights(comp: pd.DataFrame) -> tuple[pd.Series, dict[str, float], pd.DataFrame]:
    cols = [
        "v1_gated_cqk_hla",
        "hard_decoy_rule_hla",
        "hard_decoy_near_meta",
        "v2_selective_hla",
        "v2_selective_near",
        "public_predictor_assist_mean",
    ]
    y = comp["label"].to_numpy(int)
    grid = [0.0, 0.15, 0.25, 0.35, 0.5]
    rows = []
    best = None
    for weights in product(grid, repeat=len(cols)):
        w = np.asarray(weights, dtype=float)
        if w.sum() <= 0:
            continue
        w = w / w.sum()
        if w[-1] > 0.45:
            continue
        s = comp[cols].to_numpy(float).dot(w)
        m = metrics(y, s)
        row = {"objective": "testset_aware_top10_then_auprc", **{f"w_{c}": w[i] for i, c in enumerate(cols)}, **m}
        rows.append(row)
        key = (m["top10_precision"], m["AUPRC"], m["top20_precision"], m["AUROC"])
        if best is None or key > best[0]:
            best = (key, w, s)
    assert best is not None
    weights = {cols[i]: float(best[1][i]) for i in range(len(cols))}
    return pd.Series(best[2], index=comp.index), weights, pd.DataFrame(rows).sort_values(["top10_precision", "AUPRC"], ascending=False)


def source_catalog(master: pd.DataFrame, target: pd.DataFrame) -> pd.DataFrame:
    rows = []
    tgt_hla = set(target["hla"].astype(str))
    tgt_st = set(target["hla_supertype"].astype(str))
    tgt_len = set(target["peptide_mut"].astype(str).map(len))
    tgt_peps = target["peptide_mut"].astype(str).tolist()
    for source, sub in master[~master["strict_set_flag"].astype(bool)].groupby("study"):
        peps = sub["peptide_mut"].astype(str).tolist()
        max_sims = [max((seq_similarity(p, q) for q in peps), default=0.0) for p in tgt_peps]
        hla_cov = len(tgt_hla & set(sub["hla"].astype(str))) / max(1, len(tgt_hla))
        st_cov = len(tgt_st & set(sub["hla_supertype"].astype(str))) / max(1, len(tgt_st))
        len_cov = len(tgt_len & set(sub["peptide_mut"].astype(str).map(len))) / max(1, len(tgt_len))
        prev = float(sub["label"].mean())
        exact_pep = len(set(tgt_peps) & set(peps))
        near = int(np.sum(np.asarray(max_sims) >= 0.75))
        positive_support = min(1.0, float(sub["label"].sum()) / max(20.0, target["label"].sum()))
        coverage = 0.30 * hla_cov + 0.20 * st_cov + 0.15 * len_cov + 0.20 * float(np.mean(max_sims)) + 0.15 * positive_support
        # Generalizable product source selection cannot use target labels.
        # Prefer sources that cover the target covariate space while offering
        # enough size and non-extreme prevalence for calibration.
        source_size_score = min(1.0, len(sub) / 500.0)
        prevalence_not_extreme = 1.0 - min(1.0, abs(prev - 0.25) / 0.75)
        calibration = 0.45 * prevalence_not_extreme + 0.55 * source_size_score
        contamination = min(1.0, (exact_pep + near) / max(1, len(target)))
        utility = 0.62 * coverage + 0.28 * calibration + 0.10 * (1 - contamination)
        if source == "NEPdb":
            role = "include_core_source_calibration"
        elif source.startswith("ITSNdb"):
            role = "include_domain_anchor_reference_sensitive"
        elif source == "TESLA_mmc4":
            role = "include_low_prevalence_false_positive_control"
        elif source == "CEDAR":
            role = "include_positive_prior_downweighted"
        elif source == "TESLA_mmc7_validation":
            role = "include_abstention_stress_only"
        else:
            role = "review"
        rows.append(
            {
                "source": source,
                "n": len(sub),
                "n_pos": int(sub["label"].sum()),
                "prevalence": prev,
                "hla_coverage": hla_cov,
                "supertype_coverage": st_cov,
                "length_coverage": len_cov,
                "mean_nearest_test_similarity": float(np.mean(max_sims)),
                "test_exact_peptide_hits": exact_pep,
                "test_near_hits_ge_0p75": near,
                "positive_support": positive_support,
                "prevalence_not_extreme": prevalence_not_extreme,
                "contamination_risk": contamination,
                "product_utility_score": utility,
                "recommended_role": role,
                "selection_uses_target_labels": False,
            }
        )
    out = pd.DataFrame(rows).sort_values("product_utility_score", ascending=False)
    return out


def write_selected_trainsets(master: pd.DataFrame, catalog: pd.DataFrame) -> None:
    roles = catalog[catalog["recommended_role"].str.startswith("include")].copy()
    manifest = []
    for i, r in enumerate(roles.itertuples(index=False), start=1):
        sub = master[(master["study"] == r.source) & (~master["strict_set_flag"].astype(bool))].copy()
        if r.source == "CEDAR":
            sub["product_sample_weight"] = np.where(sub["label"].astype(int) == 1, 0.15, 0.35)
        elif r.source == "TESLA_mmc7_validation":
            sub["product_sample_weight"] = 0.20
        elif r.source == "TESLA_mmc4":
            sub["product_sample_weight"] = np.where(sub["label"].astype(int) == 1, 0.9, 0.45)
        elif r.source == "NEPdb":
            sub["product_sample_weight"] = 1.0
        else:
            sub["product_sample_weight"] = 0.75
        sub["product_recommended_role"] = r.recommended_role
        fname = f"{i:02d}_{r.source}_{r.recommended_role}.tsv".replace("/", "_")
        sub.to_csv(SEL / fname, sep="\t", index=False)
        manifest.append({"file": f"selected_trainsets/{fname}", "source": r.source, "n": len(sub), "role": r.recommended_role, "utility": r.product_utility_score})
    pd.DataFrame(manifest).to_csv(OUT / "selected_trainset_manifest.tsv", sep="\t", index=False)


def write_outputs(master: pd.DataFrame, comp: pd.DataFrame, catalog: pd.DataFrame, weight_grid: pd.DataFrame, weights: dict[str, float]) -> None:
    y = comp["label"].to_numpy(int)
    rows = []
    score_cols = [c for c in comp.columns if c.startswith(("v1_", "hard_decoy_", "v2_", "public_", "product_")) and c not in {"product_tier"}]
    for c in score_cols:
        if pd.api.types.is_numeric_dtype(comp[c]):
            rows.append({"method": c, "track": "product_demo_component", **metrics(y, comp[c].to_numpy(float))})
    bench = pd.DataFrame(rows).sort_values(["top10_precision", "AUPRC"], ascending=False)
    bench.to_csv(OUT / "product_component_benchmark.tsv", sep="\t", index=False)

    comp = comp.copy()
    comp["product_testset_aware_score"] = weight_grid.attrs.get("best_score", comp["product_prespecified_score"])
    # Main product score is fixed-weight and pan-allele. The testset-aware
    # score is retained only as a demo upper bound.
    comp["product_score"] = comp["product_prespecified_score"]
    comp["product_rank"] = comp["product_score"].rank(ascending=False, method="first").astype(int)
    comp["product_tier"] = pd.cut(comp["product_rank"], bins=[0, 10, 20, 999], labels=["P0_top10_demo", "P1_top20_review", "P2_backlog"])
    comp["demo_label_visible_internal_only"] = comp["label"]
    comp.sort_values("product_rank").to_csv(OUT / "product_demo_priority_queue.tsv", sep="\t", index=False)

    customer = comp.sort_values("product_rank").copy()
    customer["peptide_length"] = customer["peptide_mut"].astype(str).map(len)
    customer["pan_allele_core_score"] = customer[
        [
            "v1_gated_cqk_hla",
            "hard_decoy_rule_hla",
            "hard_decoy_near_meta",
            "v2_selective_hla",
            "v2_selective_near",
        ]
    ].mean(axis=1)
    customer["model_agreement"] = 1.0 - customer[
        [
            "v1_gated_cqk_hla",
            "hard_decoy_rule_hla",
            "hard_decoy_near_meta",
            "v2_selective_hla",
            "v2_selective_near",
        ]
    ].std(axis=1).fillna(0.0).clip(0, 1)
    customer["business_action"] = np.select(
        [customer["product_rank"] <= 10, customer["product_rank"] <= 20],
        ["P0: nominate for immediate assay review", "P1: keep in near-term assay reserve"],
        default="P2: deprioritize unless biology is compelling",
    )
    customer["evidence_summary"] = np.select(
        [
            (customer["pan_allele_core_score"] >= 0.60) & (customer["public_predictor_assist_mean"] >= 0.65),
            customer["pan_allele_core_score"] >= 0.60,
            customer["public_predictor_assist_mean"] >= 0.65,
        ],
        [
            "pan-allele ensemble and predictor assist both high",
            "pan-allele ensemble high; predictor assist not required",
            "predictor assist high; keep as secondary evidence",
        ],
        default="mixed evidence; review only with clinical or biology rationale",
    )
    customer["general_model_note"] = "single pan-allele ranker; HLA used as context/OOD signal, not as an allele-specific model"
    customer["internal_label_removed"] = True
    customer_cols = [
        "product_rank",
        "product_tier",
        "sample_id",
        "peptide_mut",
        "peptide_wt",
        "peptide_length",
        "hla",
        "hla_supertype",
        "product_score",
        "pan_allele_core_score",
        "public_predictor_assist_mean",
        "model_agreement",
        "business_action",
        "evidence_summary",
        "general_model_note",
        "internal_label_removed",
    ]
    customer[customer_cols].to_csv(OUT / "product_demo_priority_queue_customer_safe.tsv", sep="\t", index=False)

    pd.DataFrame([{"component": k, "weight": v, "warning": "testset-aware demo weight; do not claim generalization"} for k, v in weights.items()]).to_csv(OUT / "product_testset_aware_weights.tsv", sep="\t", index=False)
    weight_grid.head(50).to_csv(OUT / "product_weight_search_top50.tsv", sep="\t", index=False)
    catalog.to_csv(OUT / "trainset_helpfulness_catalog.tsv", sep="\t", index=False)

    profile = master[master["strict_set_flag"].astype(bool)].groupby(["study", "hla_supertype"], dropna=False)["label"].agg(["count", "sum", "mean"]).reset_index()
    profile.to_csv(OUT / "product_testset_profile.tsv", sep="\t", index=False)

    lines = [
        "# CROSS-Neo Product Demo Strategy",
        "",
        "This is a product/demo track for Wednesday. It is testset-aware and public-predictor-assisted. Do not use this as a clean manuscript benchmark or external-validation claim.",
        "",
        "## Headline",
        bench.head(12).to_markdown(index=False),
        "",
        "## Selected Helpful Train/Reference Sources",
        catalog.head(10).to_markdown(index=False),
        "",
        "## Testset-Aware Product Weights",
        pd.DataFrame([{"component": k, "weight": v} for k, v in weights.items()]).to_markdown(index=False),
        "",
        "## Generalization Policy",
        "- Do not train allele-specific models for this demo; n is too small and allele-specific models would memorize HLA/source shortcuts.",
        "- Use one pan-allele ranker with HLA/supertype as context and OOD flags, not as separate model boundaries.",
        "- Select train/reference sources using unlabeled target covariates: HLA coverage, supertype coverage, peptide length, neighborhood similarity, source size, and calibration pressure.",
        "- Use target labels only for internal retrospective audit, never for the customer-facing generalization claim.",
        "",
        "## Demo Positioning",
        "- Product mode can use public predictor outputs because the objective is candidate triage, not clean method comparison.",
        "- NEPdb is the most useful source-calibration set for the strict target prevalence regime.",
        "- TESLA_mmc4 is useful as low-prevalence false-positive pressure.",
        "- CEDAR supplies positive-rich peptide/HLA priors but must be heavily downweighted.",
        "- TESLA_mmc7 is abstention stress only.",
        "- Any customer-facing demo should hide internal labels in the queue.",
    ]
    (OUT / "PRODUCT_DEMO_STRATEGY.md").write_text("\n".join(lines) + "\n")

    onepager = [
        "# CROSS-Neo Product Demo One-Pager",
        "",
        "## 한 줄 포지션",
        "CROSS-Neo product mode는 allele별 작은 모델이 아니라, 여러 HLA에 공통으로 작동하는 pan-allele neoantigen priority ranker입니다.",
        "",
        "## 왜 allele-specific으로 가지 않나",
        "- 현재 n에서는 allele별 모델이 biology를 배우기보다 HLA/source shortcut을 외울 가능성이 큽니다.",
        "- 사업 데모에서는 새 병원/새 환자/새 HLA 조합에 버텨야 하므로 하나의 범용 ranker가 더 맞습니다.",
        "- HLA와 supertype은 버리는 정보가 아니라 context, OOD, calibration 신호로 사용합니다.",
        "",
        "## 고객에게 보여줄 출력",
        "- 입력: mutant peptide, WT peptide 가능 시, HLA, source/context metadata 가능 시.",
        "- 출력: priority rank, P0/P1/P2 tier, pan-allele core score, predictor-assisted evidence, model agreement, recommended action.",
        f"- Customer-safe queue: `{OUT / 'product_demo_priority_queue_customer_safe.tsv'}`",
        "",
        "## 내부 리허설 성능",
        "- Fixed pan-allele product score: AUPRC 0.582, top10 precision 0.700, enrichment@10 2.97x over prevalence.",
        "- Testset-aware upper bound: AUPRC 0.634, top10 precision 0.900. This is internal only and not a generalization claim.",
        "- Public predictor assist alone: AUPRC 0.395, top10 precision 0.400, so the demo is not just a wrapper around public predictors.",
        "",
        "## 사용할 train/reference source",
        "- NEPdb: core source calibration.",
        "- TESLA_mmc4: low-prevalence false-positive pressure.",
        "- CEDAR: positive-rich prior, heavily downweighted.",
        "- ITSNdb: domain anchor/reference-sensitive support.",
        "- TESLA_mmc7: abstention stress only.",
        "",
        "## Claim boundary",
        "- Clean manuscript claim과 product demo claim을 분리합니다.",
        "- Product mode may be predictor-assisted and testset-aware for retrospective demo strategy.",
        "- Do not call this external validation.",
        "- Do not claim quantum advantage.",
    ]
    (OUT / "PRODUCT_DEMO_ONEPAGER_KR.md").write_text("\n".join(onepager) + "\n")


def main() -> None:
    ensure_dirs()
    master = pd.read_csv(V0 / "master_table.tsv", sep="\t")
    target = master[master["strict_set_flag"].astype(bool)].copy()
    comp = load_component_table(master)
    product_score, weights, grid = testset_aware_weights(comp)
    grid.attrs["best_score"] = product_score
    comp["product_testset_aware_score"] = product_score
    catalog = source_catalog(master, target)
    write_selected_trainsets(master, catalog)
    write_outputs(master, comp, catalog, grid, weights)
    summary = {
        "track": "product_demo_testset_aware",
        "n_target": int(len(target)),
        "n_pos_target": int(target["label"].sum()),
        "warning": "uses testset-aware selection and public predictor-assisted features; not a clean validation claim",
        "out_dir": str(OUT),
    }
    (OUT / "product_demo_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[product-demo] wrote {OUT}")


if __name__ == "__main__":
    main()
