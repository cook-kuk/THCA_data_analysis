#!/usr/bin/env python3
"""Build a stronger competitor gauntlet for CROSS-Neo product/demo mode.

This is intentionally a benchmark/positioning artifact, not a clean training
pipeline. Public predictor scores are allowed here only as competitors or
predictor-assisted product signals. They remain forbidden as clean manuscript
features.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


REPO = Path(__file__).resolve().parents[1]
V0 = REPO / "project/results/cross_neo_v0"
V1 = REPO / "project/results/cross_neo_v1"
V2 = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
NEO = REPO / "project/results/p_neo_bayesian_2026_05_09"
PRODUCT = REPO / "project/results/cross_neo_product_demo_2026_05_10"
OUT = PRODUCT / "strong_competitors"


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def normalize_score(s: pd.Series, invert: bool = False) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    if invert:
        x = -x
    if x.notna().sum() == 0:
        return x
    lo = x.quantile(0.02)
    hi = x.quantile(0.98)
    if not np.isfinite(lo) or not np.isfinite(hi) or abs(hi - lo) < 1e-12:
        return pd.Series(0.5, index=x.index)
    return ((x - lo) / (hi - lo)).clip(0, 1)


def metrics(y: pd.Series, score: pd.Series) -> dict[str, float]:
    df = pd.DataFrame({"label": y, "score": score}).dropna()
    yy = df["label"].astype(int).to_numpy()
    ss = df["score"].astype(float).to_numpy()
    out: dict[str, float] = {
        "n_scored": int(len(df)),
        "n_pos": int(yy.sum()) if len(yy) else 0,
        "prevalence": float(yy.mean()) if len(yy) else np.nan,
    }
    if len(yy) and len(np.unique(yy)) == 2:
        out["AUPRC"] = float(average_precision_score(yy, ss))
        out["AUROC"] = float(roc_auc_score(yy, ss))
    else:
        out["AUPRC"] = np.nan
        out["AUROC"] = np.nan
    order = np.argsort(-ss)
    for k in [5, 10, 20]:
        kk = min(k, len(yy))
        top = yy[order[:kk]]
        out[f"top{k}_precision"] = float(top.mean()) if kk else np.nan
        out[f"recall_at_{k}"] = float(top.sum() / max(1, yy.sum())) if kk else np.nan
        out[f"enrichment_at_{k}"] = (
            float(out[f"top{k}_precision"] / out["prevalence"])
            if out.get("prevalence") and np.isfinite(out.get("prevalence", np.nan))
            else np.nan
        )
    return out


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def target_table() -> pd.DataFrame:
    master = read_tsv(V0 / "master_table.tsv")
    strict = master[master["strict_set_flag"].astype(bool)].copy()
    strict = strict.rename(columns={"peptide_mut": "peptide"})
    strict["row_key"] = (
        strict["peptide"].astype(str) + "|" + strict["hla"].astype(str) + "|" + strict["label"].astype(str)
    )
    return strict[["sample_id", "row_key", "peptide", "hla", "hla_supertype", "study", "label"]]


def merge_score(base: pd.DataFrame, path: Path, score_col: str, out_col: str, invert: bool = False) -> pd.DataFrame:
    df = read_tsv(path)
    if df.empty or score_col not in df.columns:
        base[out_col] = np.nan
        return base
    if "HLA_norm" in df.columns and "hla" not in df.columns:
        df = df.rename(columns={"HLA_norm": "hla"})
    if "peptide_mut" in df.columns and "peptide" not in df.columns:
        df = df.rename(columns={"peptide_mut": "peptide"})
    need = ["peptide", "hla", "label", score_col]
    if not set(need).issubset(df.columns):
        base[out_col] = np.nan
        return base
    sub = df[need].copy()
    sub["row_key"] = sub["peptide"].astype(str) + "|" + sub["hla"].astype(str) + "|" + sub["label"].astype(str)
    sub[out_col] = normalize_score(sub[score_col], invert=invert)
    sub = sub.groupby("row_key", as_index=False)[out_col].mean()
    return base.merge(sub, on="row_key", how="left")


def add_public_competitors(base: pd.DataFrame) -> pd.DataFrame:
    scored = NEO / "wave3_algorithm_sweep/scored_master.tsv"
    base = merge_score(base, scored, "mhcflurry_presentation", "MHCflurry_presentation")
    base = merge_score(base, scored, "mhcflurry_affinity", "MHCflurry_affinity_inverse", invert=True)
    base = merge_score(base, scored, "bigmhc_im", "BigMHC_IM")
    base = merge_score(base, scored, "bigmhc_el", "BigMHC_EL")
    base = merge_score(base, scored, "prime_score", "PRIME_score")
    base = merge_score(base, scored, "prime_rank", "PRIME_rank_inverse", invert=True)

    wave11 = NEO / "wave11/wave11_predictions"
    base = merge_score(base, wave11 / "DeepImmuno__itsndb.tsv", "score_deepimmuno", "DeepImmuno")
    base = merge_score(base, wave11 / "TransPHLA__itsndb.tsv", "score_transphla", "TransPHLA")
    base = merge_score(base, wave11 / "NetMHCpan_4.1__itsndb.tsv", "score_netmhcpan", "NetMHCpan_4_1")
    base = merge_score(base, wave11 / "MHCflurry__itsndb.tsv", "score_mhcflurry", "MHCflurry_wave11")
    base = merge_score(base, wave11 / "PRIME__itsndb.tsv", "score_prime", "PRIME_wave11")
    base = merge_score(base, wave11 / "BigMHC_IM__itsndb.tsv", "score_bigmhc_im", "BigMHC_IM_wave11")

    base = merge_score(base, NEO / "wave9/predictions_tscape.tsv", "score_tscape", "TSCAPE")
    base = merge_score(base, NEO / "wave9/predictions_mhcnuggets.tsv", "score_mhcnuggets", "MHCnuggets")
    base = merge_score(base, NEO / "wave9/predictions_netmhcstabpan.tsv", "score_stabpan", "NetMHCstabpan")
    base = merge_score(base, NEO / "wave3_netmhcpan/predictions.tsv", "score_netmhcpan", "NetMHCpan_wave3")
    base = merge_score(base, NEO / "wave3_deepimmuno/predictions.tsv", "score_deepimmuno", "DeepImmuno_wave3")
    base = merge_score(base, NEO / "wave3_transphla/predictions.tsv", "score_transphla", "TransPHLA_wave3")

    public_sets = {
        "Public_binding_mean": ["MHCflurry_presentation", "MHCflurry_affinity_inverse", "BigMHC_EL", "NetMHCpan_4_1"],
        "Public_immunogenicity_mean": ["BigMHC_IM", "DeepImmuno", "PRIME_score", "TransPHLA", "TSCAPE"],
        "Public_all_mean": [
            "MHCflurry_presentation",
            "MHCflurry_affinity_inverse",
            "BigMHC_IM",
            "BigMHC_EL",
            "PRIME_score",
            "PRIME_rank_inverse",
            "DeepImmuno",
            "TransPHLA",
            "TSCAPE",
            "MHCnuggets",
            "NetMHCstabpan",
        ],
    }
    for name, cols in public_sets.items():
        have = [c for c in cols if c in base.columns]
        base[name] = base[have].mean(axis=1) if have else np.nan
    return base


def add_cross_neo_scores(base: pd.DataFrame) -> pd.DataFrame:
    product = read_tsv(PRODUCT / "product_demo_priority_queue.tsv")
    for c in [
        "product_prespecified_score",
        "product_testset_aware_score",
        "v1_gated_cqk_hla",
        "hard_decoy_rule_hla",
        "hard_decoy_near_meta",
        "v2_selective_hla",
        "v2_selective_near",
        "v2_selective_exact",
    ]:
        if not product.empty and c in product.columns:
            sub = product[["sample_id", c]].copy()
            base = base.merge(sub, on="sample_id", how="left")
    return base


def v2_locked_rows() -> pd.DataFrame:
    paths = [
        V2 / "metrics/selective_ensemble_v2_1_comparison.tsv",
        V2 / "metrics/fast_esm2_qk_gate_primary_comparison.tsv",
    ]
    rows = []
    for path in paths:
        df = read_tsv(path)
        if df.empty:
            continue
        if "pre_policy_best" in df.columns:
            for r in df.itertuples(index=False):
                rows.append(
                    {
                        "benchmark_scope": "locked_split_best_available",
                        "method": f"v2_pre_policy_best::{r.split_name}::{r.pre_policy_best}",
                        "method_family": "CROSS_Neo_v2_internal",
                        "n_scored": np.nan,
                        "n_pos": np.nan,
                        "prevalence": np.nan,
                        "AUPRC": r.pre_policy_AUPRC,
                        "AUROC": np.nan,
                        "top5_precision": np.nan,
                        "top10_precision": r.pre_policy_top10,
                        "top20_precision": np.nan,
                        "coverage": np.nan,
                        "claim_boundary": "internal locked split; not external validation",
                    }
                )
        elif "combined_best" in df.columns:
            for r in df.itertuples(index=False):
                rows.append(
                    {
                        "benchmark_scope": "locked_split_best_available",
                        "method": f"fast_esm2_qk_gate::{r.split_name}::{r.combined_best}",
                        "method_family": "CROSS_Neo_v2_internal",
                        "n_scored": np.nan,
                        "n_pos": np.nan,
                        "prevalence": np.nan,
                        "AUPRC": r.combined_AUPRC,
                        "AUROC": np.nan,
                        "top5_precision": np.nan,
                        "top10_precision": r.combined_top10,
                        "top20_precision": np.nan,
                        "coverage": np.nan,
                        "claim_boundary": "internal locked split; not external validation",
                    }
                )
    return pd.DataFrame(rows)


def wave11_context_rows() -> pd.DataFrame:
    rows = []
    df = read_tsv(NEO / "wave11/wave11_mega_table_expanded.tsv")
    if not df.empty:
        for r in df.itertuples(index=False):
            rows.append(
                {
                    "benchmark_scope": f"wave11_public_test::{r.test_bundle}",
                    "method": r.algorithm,
                    "method_family": "public_or_prior_internal_competitor",
                    "n_scored": r.n,
                    "n_pos": r.n_pos,
                    "prevalence": np.nan,
                    "AUPRC": np.nan,
                    "AUROC": r.AUROC,
                    "top5_precision": np.nan,
                    "top10_precision": np.nan,
                    "top20_precision": np.nan,
                    "coverage": np.nan,
                    "claim_boundary": "context only; AUROC-only public benchmark, possible overlap/training-status caveats",
                }
            )
    return pd.DataFrame(rows)


def strict_gauntlet_rows(score_matrix: pd.DataFrame) -> pd.DataFrame:
    competitors = {
        "Product_fixed_pan_allele_score": ("product_prespecified_score", "CROSS_Neo_product"),
        "Product_testset_aware_upper_bound": ("product_testset_aware_score", "CROSS_Neo_product_upper_bound"),
        "v2_selective_exact": ("v2_selective_exact", "CROSS_Neo_v2_internal"),
        "v2_selective_hla": ("v2_selective_hla", "CROSS_Neo_v2_internal"),
        "v2_selective_near": ("v2_selective_near", "CROSS_Neo_v2_internal"),
        "hard_decoy_rule_hla": ("hard_decoy_rule_hla", "CROSS_Neo_v1_internal"),
        "hard_decoy_near_meta": ("hard_decoy_near_meta", "CROSS_Neo_v1_internal"),
        "v1_gated_cqk_hla": ("v1_gated_cqk_hla", "CROSS_Neo_v1_internal"),
        "Public_all_mean": ("Public_all_mean", "public_predictor_panel"),
        "Public_immunogenicity_mean": ("Public_immunogenicity_mean", "public_predictor_panel"),
        "Public_binding_mean": ("Public_binding_mean", "public_predictor_panel"),
        "BigMHC_IM": ("BigMHC_IM", "public_predictor"),
        "BigMHC_EL": ("BigMHC_EL", "public_predictor"),
        "BigMHC_IM_wave11": ("BigMHC_IM_wave11", "public_predictor"),
        "DeepImmuno": ("DeepImmuno", "public_predictor"),
        "DeepImmuno_wave3": ("DeepImmuno_wave3", "public_predictor"),
        "PRIME_score": ("PRIME_score", "public_predictor"),
        "PRIME_rank_inverse": ("PRIME_rank_inverse", "public_predictor"),
        "PRIME_wave11": ("PRIME_wave11", "public_predictor"),
        "TransPHLA": ("TransPHLA", "public_predictor"),
        "TransPHLA_wave3": ("TransPHLA_wave3", "public_predictor"),
        "MHCflurry_presentation": ("MHCflurry_presentation", "public_predictor"),
        "MHCflurry_affinity_inverse": ("MHCflurry_affinity_inverse", "public_predictor"),
        "MHCflurry_wave11": ("MHCflurry_wave11", "public_predictor"),
        "NetMHCpan_4_1": ("NetMHCpan_4_1", "public_predictor"),
        "NetMHCpan_wave3": ("NetMHCpan_wave3", "public_predictor"),
        "NetMHCstabpan": ("NetMHCstabpan", "public_predictor"),
        "MHCnuggets": ("MHCnuggets", "public_predictor"),
        "TSCAPE": ("TSCAPE", "public_predictor"),
    }
    rows = []
    y = score_matrix["label"]
    for method, (col, family) in competitors.items():
        if col not in score_matrix.columns:
            continue
        m = metrics(y, score_matrix[col])
        if m["n_scored"] == 0:
            continue
        rows.append(
            {
                "benchmark_scope": "strict_product_demo_set_n89",
                "method": method,
                "method_family": family,
                **m,
                "coverage": float(m["n_scored"] / len(score_matrix)) if len(score_matrix) else np.nan,
                "claim_boundary": (
                    "customer-demo main score; internal retrospective only"
                    if method == "Product_fixed_pan_allele_score"
                    else "testset-aware upper bound; do not claim generalization"
                    if method == "Product_testset_aware_upper_bound"
                    else "competitor benchmark only; public predictor scores are not clean CROSS-Neo training features"
                    if "public" in family
                    else "internal locked/product component; not external validation"
                ),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values(["benchmark_scope", "top10_precision", "AUPRC"], ascending=[True, False, False])


def write_report(all_rows: pd.DataFrame, score_matrix: pd.DataFrame) -> None:
    strict = all_rows[all_rows["benchmark_scope"] == "strict_product_demo_set_n89"].copy()
    product = strict[strict["method"] == "Product_fixed_pan_allele_score"].iloc[0]
    public_best = strict[strict["method_family"].str.contains("public", na=False)].sort_values(
        ["top10_precision", "AUPRC"], ascending=False
    ).head(1)
    public_best_row = public_best.iloc[0] if len(public_best) else None
    top = strict.head(15)

    lines = [
        "# Strong Competitor Gauntlet",
        "",
        "This is a business/product benchmark artifact. Public predictor scores are used only as competitors or predictor-assisted context, not as clean manuscript training features.",
        "",
        "## Strict Product Demo Set",
        top[
            [
                "method",
                "method_family",
                "n_scored",
                "AUPRC",
                "AUROC",
                "top5_precision",
                "top10_precision",
                "top20_precision",
                "coverage",
            ]
        ].to_markdown(index=False),
        "",
        "## Bottom Line",
        f"- Fixed pan-allele product score: AUPRC {product.AUPRC:.3f}, top10 precision {product.top10_precision:.3f}.",
    ]
    if public_best_row is not None:
        lines.append(
            f"- Best public-only competitor in this gauntlet: {public_best_row.method}, "
            f"AUPRC {public_best_row.AUPRC:.3f}, top10 precision {public_best_row.top10_precision:.3f}."
        )
    lines += [
        "- Product_testset_aware_upper_bound is internal rehearsal only and should not be used as a generalization claim.",
        "- The right customer-facing claim is stronger practical prioritization on an internal retrospective demo set, not external validation.",
        "",
        "## Strong Competitors Included",
        "- BigMHC, DeepImmuno, PRIME, TransPHLA, MHCflurry, NetMHCpan, NetMHCstabpan, MHCnuggets, T-SCAPE.",
        "- CROSS-Neo v1/v2 internal challengers, including hard-decoy, gated C+QK, selective ensemble, and ESM2/QK gate summaries.",
        "",
        "## Strong Competitors Still Missing Locally",
        "- IMPROVE: strong broad-scale neoepitope immunogenicity model, but no local score file was found.",
        "- MixMHCpred standalone: indirectly represented through PRIME-style scoring, but no standalone local output was found for the strict set.",
        "- DeepHLApan / CIImm / GraphMHC: useful literature comparators, but no local runnable score files were found.",
        "",
        "## Files",
        f"- Score matrix: `{OUT / 'strong_competitor_score_matrix.tsv'}`",
        f"- Gauntlet metrics: `{OUT / 'strong_competitor_gauntlet.tsv'}`",
        f"- Context public benchmark: `{OUT / 'strong_competitor_context_public_benchmarks.tsv'}`",
        f"- Missing competitor manifest: `{OUT / 'missing_strong_competitor_manifest.tsv'}`",
        f"- Figure: `{OUT / 'figure_strong_competitor_gauntlet.png'}`",
    ]
    (OUT / "STRONG_COMPETITOR_GAUNTLET.md").write_text("\n".join(lines) + "\n")

    summary = {
        "strict_product_fixed_AUPRC": float(product.AUPRC),
        "strict_product_fixed_top10": float(product.top10_precision),
        "best_public_method": None if public_best_row is None else str(public_best_row.method),
        "best_public_AUPRC": None if public_best_row is None else float(public_best_row.AUPRC),
        "best_public_top10": None if public_best_row is None else float(public_best_row.top10_precision),
        "n_strict": int(len(score_matrix)),
    }
    (OUT / "strong_competitor_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def write_missing_manifest() -> None:
    rows = [
        {
            "competitor": "IMPROVE",
            "priority": "high",
            "why_it_matters": "large broad-scale T-cell-recognition neoepitope model; strong reviewer-facing competitor",
            "local_status": "not_found",
            "needed_input": "peptide, HLA, mutation context/features or runnable model-specific feature set",
            "url_or_search": "IMPROVE neoepitope immunogenicity 17500 candidates 467 recognized PubMed 38633261",
        },
        {
            "competitor": "MixMHCpred standalone",
            "priority": "high",
            "why_it_matters": "widely used HLA-I ligand presentation predictor and PRIME dependency",
            "local_status": "not_found_as_standalone_strict_scores",
            "needed_input": "peptide, HLA allele",
            "url_or_search": "MixMHCpred HLA-I ligand predictor",
        },
        {
            "competitor": "DeepHLApan",
            "priority": "medium",
            "why_it_matters": "deep pan-allele immunogenicity/presentation comparator in neoantigen reviews",
            "local_status": "not_found",
            "needed_input": "peptide, HLA allele or pseudo-sequence",
            "url_or_search": "DeepHLApan neoantigen immunogenicity predictor",
        },
        {
            "competitor": "GraphMHC",
            "priority": "medium",
            "why_it_matters": "structure/GNN-style neoantigen comparator",
            "local_status": "not_found",
            "needed_input": "peptide, HLA, structure/graph features or runnable package",
            "url_or_search": "GraphMHC neoantigen prediction graph neural network structure PubMed 38536842",
        },
        {
            "competitor": "CIImm",
            "priority": "medium",
            "why_it_matters": "class-I pMHC immunogenicity predictor mentioned in comparative reviews",
            "local_status": "not_found",
            "needed_input": "peptide, HLA allele",
            "url_or_search": "CIImm T cell class I pMHC immunogenicity predictor",
        },
    ]
    pd.DataFrame(rows).to_csv(OUT / "missing_strong_competitor_manifest.tsv", sep="\t", index=False)


def write_figure(strict_rows: pd.DataFrame) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    plot = strict_rows.copy()
    plot = plot.sort_values(["top10_precision", "AUPRC"], ascending=False).head(16)
    labels = plot["method"].str.replace("_", " ", regex=False)
    y = np.arange(len(plot))
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = np.where(plot["method_family"].str.contains("CROSS", na=False), "#2f6f9f", "#b36b2c")
    ax.barh(y, plot["top10_precision"], color=colors, alpha=0.9, label="top10 precision")
    ax.scatter(plot["AUPRC"], y, color="#111111", s=35, label="AUPRC")
    ax.axvline(plot["prevalence"].iloc[0], color="#777777", linestyle="--", linewidth=1, label="prevalence")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Metric value")
    ax.set_title("Internal strict product demo: stronger competitor gauntlet")
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "figure_strong_competitor_gauntlet.png", dpi=220)
    plt.close(fig)


def main() -> None:
    ensure_dirs()
    score_matrix = target_table()
    score_matrix = add_public_competitors(score_matrix)
    score_matrix = add_cross_neo_scores(score_matrix)
    score_matrix.to_csv(OUT / "strong_competitor_score_matrix.tsv", sep="\t", index=False)

    strict_rows = strict_gauntlet_rows(score_matrix)
    context = pd.concat([v2_locked_rows(), wave11_context_rows()], ignore_index=True)
    all_rows = pd.concat([strict_rows, context], ignore_index=True)
    all_rows.to_csv(OUT / "strong_competitor_gauntlet.tsv", sep="\t", index=False)
    context.to_csv(OUT / "strong_competitor_context_public_benchmarks.tsv", sep="\t", index=False)
    write_missing_manifest()
    write_figure(strict_rows)
    write_report(all_rows, score_matrix)
    print(f"[strong-competitors] wrote {OUT}")


if __name__ == "__main__":
    main()
