#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import (
    ensure_run_dir,
    metric_binary,
    patient_topn_metrics,
    peptide_features,
    precision_recall_at_k,
    safe_read_table,
    write_md,
    write_tsv,
)


NMI_ALLOWED_BRANCH_TOKENS = [
    "Structure_LR",
    "Wave8",
    "TCR",
    "SelfSim",
    "ESM2_Bayesian",
    "GP_quantum",
    "VQC",
    "W7A",
    "W7B",
    "quantum",
    "qk_clean",
    "qk_no_anchor",
    "qk_quantum",
    "nested_learned_gate_C_QK_structure",
    "Stack_mean",
    "Stack_median",
    "Stack_LR",
]

NMI_BLOCKED_TOKENS = [
    "KG_GA",
    "BAR-Neo",
    "BigMHC",
    "MHCflurry",
    "NetMHC",
    "PRIME",
    "DeepImmuno",
    "public",
    "product",
    "impact_portfolio",
]

PUBLIC_COMPARATOR_TOKENS = ["BigMHC", "MHCflurry", "NetMHC", "PRIME"]


def is_nmi_allowed_branch(name: object) -> bool:
    s = str(name)
    low = s.lower()
    if any(tok.lower() in low for tok in NMI_BLOCKED_TOKENS):
        return False
    return any(tok.lower() in low for tok in NMI_ALLOWED_BRANCH_TOKENS)


def clean_numeric_context(d: pd.DataFrame) -> pd.DataFrame:
    out = peptide_features(d)
    for c in [
        "tcr_motif_score",
        "self_exact_match",
        "self_hamming1_count_log",
        "self_hamming2_count_log",
        "self_blosum_max",
        "foreignness_score",
        "self_similarity_score",
        "tcr_recognition_score",
        "expression_tpm",
        "vaf",
        "apm_score",
        "tap1_expr",
        "tap2_expr",
    ]:
        if c in d.columns:
            out[c] = pd.to_numeric(d[c], errors="coerce")
    pep = d.get("peptide_mut", pd.Series("", index=d.index)).fillna("").astype(str)
    for aa in list("ACDEFGHIKLMNPQRSTVWY"):
        out[f"aa_frac_{aa}"] = pep.map(lambda x, aa=aa: x.count(aa) / len(x) if x else np.nan)
    if "hla_allele" in d:
        locus = d["hla_allele"].fillna("").astype(str).str.extract(r"HLA-?([A-Z]+)\*", expand=False).fillna("UNK")
        out = pd.concat([out, pd.get_dummies(locus, prefix="hla_locus", dtype=float)], axis=1)
    if "hla_supertype" in d:
        supertype = d["hla_supertype"].fillna("UNK").astype(str)
        top = supertype.value_counts().head(18).index
        supertype = supertype.where(supertype.isin(top), "OTHER")
        out = pd.concat([out, pd.get_dummies(supertype, prefix="hla_supertype", dtype=float)], axis=1)
    return out


def make_branch_features(local: pd.DataFrame, canon: pd.DataFrame) -> pd.DataFrame:
    if local.empty:
        return pd.DataFrame({"candidate_id": canon["candidate_id"]})
    d = local[local["local_model_name"].map(is_nmi_allowed_branch)].copy()
    if d.empty:
        return pd.DataFrame({"candidate_id": canon["candidate_id"]})
    d["normalized_score"] = pd.to_numeric(d["normalized_score"], errors="coerce")
    wide = d.pivot_table(index="candidate_id", columns="local_model_name", values="normalized_score", aggfunc="mean").reset_index()
    return wide


def fit_oof(X: pd.DataFrame, y: pd.Series, groups: pd.Series, kind: str) -> pd.Series:
    from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold, StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.loc[:, X.notna().any(axis=0)].copy()
    y = pd.to_numeric(y, errors="coerce").fillna(0).astype(int)
    pred = pd.Series(np.nan, index=X.index, dtype=float)
    groups = groups.fillna("missing").astype(str)
    if groups.nunique() >= 3:
        splitter = GroupKFold(n_splits=min(5, groups.nunique()))
        splits = splitter.split(X, y, groups)
    else:
        n = min(5, int(y.value_counts().min())) if y.nunique() > 1 else 2
        splitter = StratifiedKFold(n_splits=max(2, n), shuffle=True, random_state=29)
        splits = splitter.split(X, y)
    for tr, te in splits:
        if y.iloc[tr].nunique() < 2:
            pred.iloc[te] = float(y.iloc[tr].mean())
            continue
        if kind == "lr":
            model = make_pipeline(
                SimpleImputer(strategy="median", add_indicator=True),
                StandardScaler(with_mean=False),
                LogisticRegression(max_iter=2000, class_weight="balanced", random_state=29),
            )
        elif kind == "rf":
            model = make_pipeline(
                SimpleImputer(strategy="median", add_indicator=True),
                RandomForestClassifier(
                    n_estimators=120,
                    max_features="sqrt",
                    min_samples_leaf=6,
                    class_weight="balanced_subsample",
                    n_jobs=-1,
                    random_state=29,
                ),
            )
        else:
            model = make_pipeline(
                SimpleImputer(strategy="median", add_indicator=True),
                ExtraTreesClassifier(
                    n_estimators=120,
                    max_features="sqrt",
                    min_samples_leaf=6,
                    class_weight="balanced",
                    n_jobs=-1,
                    random_state=29,
                ),
            )
        model.fit(X.iloc[tr], y.iloc[tr])
        pred.iloc[te] = model.predict_proba(X.iloc[te])[:, 1]
    return pred.clip(0, 1)


def transparent_nmi_score(feat: pd.DataFrame, branch_cols: list[str]) -> pd.Series:
    parts = []
    weights = []
    sequence_cols = [c for c in feat.columns if c.startswith(("aa_frac_", "hydrophobic", "charged", "aromatic", "polar", "mut_wt", "peptide_length"))]
    tcr_cols = [c for c in feat.columns if any(tok in c.lower() for tok in ["tcr", "selfsim", "self_", "waves"])]
    structure_cols = [c for c in branch_cols if any(tok in c.lower() for tok in ["structure", "esm2"])]
    quantum_cols = [c for c in branch_cols if any(tok in c.lower() for tok in ["quantum", "qk", "vqc", "gp_"])]
    stack_cols = [c for c in branch_cols if any(tok in c.lower() for tok in ["w7a", "w7b", "stack"])]
    for cols, w in [
        (sequence_cols, 0.18),
        (tcr_cols, 0.28),
        (structure_cols, 0.20),
        (quantum_cols, 0.17),
        (stack_cols, 0.17),
    ]:
        cols = [c for c in cols if c in feat.columns]
        if cols:
            x = feat[cols].apply(pd.to_numeric, errors="coerce")
            x = x.rank(pct=True).mean(axis=1, skipna=True)
            parts.append(x)
            weights.append(w)
    if not parts:
        return pd.Series(0.5, index=feat.index)
    score = sum(p.fillna(0.5) * w for p, w in zip(parts, weights)) / sum(weights)
    return score.clip(0, 1)


def add_metrics(rows: list[dict], df: pd.DataFrame, score_col: str, label_col: str, model: str, split: str, features: str, clean_allowed: bool = True) -> None:
    eval_df = df[df[label_col].notna()].copy()
    eval_df[score_col] = pd.to_numeric(eval_df[score_col], errors="coerce")
    eval_df = eval_df[eval_df[score_col].notna()]
    if eval_df.empty:
        return
    m = metric_binary(eval_df[label_col], eval_df[score_col])
    for k in [10, 20, 34, 50]:
        m.update(precision_recall_at_k(eval_df, score_col, label_col, k))
    for n in [20, 34]:
        m.update(patient_topn_metrics(eval_df, score_col, label_col, n))
    m.update(
        {
            "track": "nmi_clean_method_track",
            "model": model,
            "split": split,
            "features": features,
            "score_col": score_col,
            "clean_track_allowed": clean_allowed,
        }
    )
    rows.append(m)


def write_figures(outdir: Path, leaderboard: pd.DataFrame) -> None:
    try:
        import matplotlib.pyplot as plt

        figdir = outdir / "figures"
        figdir.mkdir(parents=True, exist_ok=True)
        top = leaderboard.sort_values("AUPRC", ascending=False).head(12).copy()
        fig, ax = plt.subplots(figsize=(12, 6), facecolor="#080b12")
        ax.set_facecolor("#111827")
        colors = ["#26d9ff" if str(x).startswith("NMI") else "#a78bfa" for x in top["model"]]
        ax.barh(top["model"], top["AUPRC"], color=colors)
        ax.invert_yaxis()
        ax.set_xlabel("AUPRC")
        ax.set_title("NMI clean method vs clean/public comparators", color="white", fontsize=15)
        ax.tick_params(colors="white")
        ax.xaxis.label.set_color("white")
        for spine in ax.spines.values():
            spine.set_color("#263348")
        ax.grid(axis="x", color="#263348", alpha=0.6)
        fig.tight_layout()
        fig.savefig(figdir / "nmi_clean_method_leaderboard.png", dpi=220)
        plt.close(fig)
    except Exception:
        pass


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    label_col = "label_immunogenicity"

    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    canon[label_col] = pd.to_numeric(canon[label_col], errors="coerce")
    local = safe_read_table(outdir / "local_scores" / "all_local_scores.tsv.gz")
    branch_wide = make_branch_features(local, canon)
    d = canon.merge(branch_wide, on="candidate_id", how="left")

    context = clean_numeric_context(d)
    branch_cols = [c for c in branch_wide.columns if c != "candidate_id" and c in d.columns]
    branch_feat = d[branch_cols].apply(pd.to_numeric, errors="coerce") if branch_cols else pd.DataFrame(index=d.index)
    branch_present = branch_feat.notna().astype(float).add_suffix("__available") if not branch_feat.empty else pd.DataFrame(index=d.index)
    core_feat = context.copy()
    fusion_feat = pd.concat([context, branch_feat, branch_present], axis=1)

    use = d[label_col].notna()
    groups = d.loc[use, "source_study"] if "source_study" in d.columns else d.loc[use, "dataset_source"]
    d.loc[use, "NMI_core_LR_source_oof"] = fit_oof(core_feat.loc[use], d.loc[use, label_col], groups, "lr")
    d.loc[use, "NMI_branch_LR_source_oof"] = fit_oof(fusion_feat.loc[use], d.loc[use, label_col], groups, "lr")
    d.loc[use, "NMI_branch_RF_source_oof"] = fit_oof(fusion_feat.loc[use], d.loc[use, label_col], groups, "rf")
    d.loc[use, "NMI_branch_ET_source_oof"] = fit_oof(fusion_feat.loc[use], d.loc[use, label_col], groups, "et")
    d["NMI_transparent_score"] = transparent_nmi_score(pd.concat([context, branch_feat], axis=1), branch_cols)
    for col in ["NMI_core_LR_source_oof", "NMI_branch_LR_source_oof", "NMI_branch_RF_source_oof", "NMI_branch_ET_source_oof"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d["NMI_clean_consensus"] = d[["NMI_branch_LR_source_oof", "NMI_branch_RF_source_oof", "NMI_branch_ET_source_oof", "NMI_transparent_score"]].mean(axis=1, skipna=True).clip(0, 1)
    def mean_existing(cols: list[str]) -> pd.Series:
        have = [c for c in cols if c in d.columns]
        if not have:
            return pd.Series(np.nan, index=d.index)
        return d[have].apply(pd.to_numeric, errors="coerce").mean(axis=1, skipna=True)

    d["NMI_all_branch_mean"] = mean_existing(
        [
            "W7B_stacked",
            "ESM2_Bayesian",
            "W7A_QK_only",
            "W7A_full",
            "Structure_LR",
            "Wave8_TCR_SelfSim_no_exact",
            "Wave8_TCR_motif_only",
            "Wave8_TCR_SelfSim_full",
            "GP_quantum",
            "qk_clean_no_tcr_gamma1",
        ]
    )
    d["NMI_tcr_qk_structure"] = (
        0.35 * pd.to_numeric(d.get("Wave8_TCR_SelfSim_no_exact"), errors="coerce")
        + 0.25 * pd.to_numeric(d.get("W7A_QK_only"), errors="coerce")
        + 0.20 * pd.to_numeric(d.get("Structure_LR"), errors="coerce")
        + 0.20 * pd.to_numeric(d.get("ESM2_Bayesian"), errors="coerce")
    )
    d["NMI_equal_top5"] = mean_existing(["W7B_stacked", "ESM2_Bayesian", "W7A_QK_only", "W7A_full", "Structure_LR"])
    d["NMI_w7b_esm2_qk"] = (
        0.45 * pd.to_numeric(d.get("W7B_stacked"), errors="coerce")
        + 0.35 * pd.to_numeric(d.get("ESM2_Bayesian"), errors="coerce")
        + 0.20 * pd.to_numeric(d.get("W7A_QK_only"), errors="coerce")
    )

    flags_path = outdir / "metrics" / "leakage_candidate_flags.tsv"
    if flags_path.exists():
        flags = safe_read_table(flags_path)[["candidate_id", "any_existing_overlap_flag"]]
        d = d.merge(flags, on="candidate_id", how="left")
        d["any_existing_overlap_flag"] = d["any_existing_overlap_flag"].fillna(False).astype(bool)
    else:
        d["any_existing_overlap_flag"] = False

    rows: list[dict] = []
    for col, model, split, features in [
        ("NMI_all_branch_mean", "NMI all-branch mean", "nmi_eligible_locked_formula", "locked clean local branch formula on NMI-eligible rows"),
        ("NMI_tcr_qk_structure", "NMI TCR-QK-Structure", "nmi_eligible_locked_formula", "locked TCR/quantum/structure/ESM2 formula on NMI-eligible rows"),
        ("NMI_equal_top5", "NMI equal top5", "nmi_eligible_locked_formula", "W7B/ESM2/QK/full/Structure equal-weight formula"),
        ("NMI_w7b_esm2_qk", "NMI W7B-ESM2-QK", "nmi_eligible_locked_formula", "W7B/ESM2/QK weighted formula"),
        ("NMI_clean_consensus", "NMI-clean consensus", "source_grouped_oof_fallback", "peptide/HLA/self/TCR/tumor context + own branch fusion"),
        ("NMI_branch_RF_source_oof", "NMI-branch RF", "source_grouped_oof_fallback", "own branch fusion, source-grouped OOF"),
        ("NMI_branch_ET_source_oof", "NMI-branch ExtraTrees", "source_grouped_oof_fallback", "own branch fusion, source-grouped OOF"),
        ("NMI_branch_LR_source_oof", "NMI-branch LR", "source_grouped_oof_fallback", "transparent own branch fusion, source-grouped OOF"),
        ("NMI_core_LR_source_oof", "NMI-core LR", "source_grouped_oof_fallback", "peptide/HLA/self/TCR/tumor context only"),
        ("NMI_transparent_score", "NMI-transparent formula", "source_grouped_oof_fallback", "fixed clean branch-weight formula"),
    ]:
        add_metrics(rows, d, col, label_col, model, split, features)

    strict = d[~d["any_existing_overlap_flag"]].copy()
    for col, model, split, features in [
        ("NMI_all_branch_mean", "NMI all-branch mean", "strict_no_existing_overlap_nmi_formula", "locked clean local branch formula on NMI-eligible rows"),
        ("NMI_tcr_qk_structure", "NMI TCR-QK-Structure", "strict_no_existing_overlap_nmi_formula", "locked TCR/quantum/structure/ESM2 formula on NMI-eligible rows"),
        ("NMI_equal_top5", "NMI equal top5", "strict_no_existing_overlap_nmi_formula", "W7B/ESM2/QK/full/Structure equal-weight formula"),
        ("NMI_w7b_esm2_qk", "NMI W7B-ESM2-QK", "strict_no_existing_overlap_nmi_formula", "W7B/ESM2/QK weighted formula"),
        ("NMI_clean_consensus", "NMI-clean consensus", "strict_no_existing_overlap_fallback", "peptide/HLA/self/TCR/tumor context + own branch fusion"),
        ("NMI_branch_RF_source_oof", "NMI-branch RF", "strict_no_existing_overlap_fallback", "own branch fusion, source-grouped OOF"),
        ("NMI_branch_ET_source_oof", "NMI-branch ExtraTrees", "strict_no_existing_overlap_fallback", "own branch fusion, source-grouped OOF"),
        ("NMI_branch_LR_source_oof", "NMI-branch LR", "strict_no_existing_overlap_fallback", "transparent own branch fusion, source-grouped OOF"),
        ("NMI_core_LR_source_oof", "NMI-core LR", "strict_no_existing_overlap_fallback", "peptide/HLA/self/TCR/tumor context only"),
        ("NMI_transparent_score", "NMI-transparent formula", "strict_no_existing_overlap_fallback", "fixed clean branch-weight formula"),
    ]:
        add_metrics(rows, strict, col, label_col, model, split, features)

    ext_path = outdir / "external_scores" / "all_external_scores.tsv.gz"
    if ext_path.exists():
        ext = safe_read_table(ext_path)
        ext = ext[ext["model_name"].astype(str).str.contains("|".join(PUBLIC_COMPARATOR_TOKENS), case=False, na=False)].copy()
        if not ext.empty:
            ext["normalized_score"] = pd.to_numeric(ext["normalized_score"], errors="coerce")
            ext_wide = ext.pivot_table(index="candidate_id", columns="model_name", values="normalized_score", aggfunc="mean").reset_index()
            comp = d[["candidate_id", "patient_id", label_col, "any_existing_overlap_flag"]].merge(ext_wide, on="candidate_id", how="left")
            for c in [x for x in ext_wide.columns if x != "candidate_id"]:
                add_metrics(rows, comp, c, label_col, f"public comparator: {c}", "retrospective_public_comparator", "frozen public predictor score", clean_allowed=False)
            eligible_ids = set(d.loc[d["NMI_all_branch_mean"].notna(), "candidate_id"].astype(str))
            comp_eligible = comp[comp["candidate_id"].astype(str).isin(eligible_ids)].copy()
            for c in [x for x in ext_wide.columns if x != "candidate_id"]:
                add_metrics(rows, comp_eligible, c, label_col, f"same-set public comparator: {c}", "nmi_eligible_same_rows", "frozen public predictor on NMI-eligible rows", clean_allowed=False)

    lb = pd.DataFrame(rows)
    ordered = [
        "track",
        "model",
        "split",
        "features",
        "clean_track_allowed",
        "n",
        "positives",
        "AUROC",
        "AUPRC",
        "Precision@10",
        "Precision@20",
        "Precision@34",
        "Recall@20",
        "Recall@34",
        "patient_hit_rate@20",
        "patient_recall@20",
        "patient_hit_rate@34",
        "patient_recall@34",
        "patients_evaluated",
        "calibration_brier",
    ]
    for c in ordered:
        if c not in lb:
            lb[c] = np.nan
    lb = lb[ordered + [c for c in lb.columns if c not in ordered]].sort_values(["split", "AUPRC", "patient_hit_rate@34"], ascending=[True, False, False])
    write_tsv(lb, outdir / "metrics" / "nmi_clean_method_leaderboard.tsv")

    pred_cols = [
        "candidate_id",
        "patient_id",
        "dataset_source",
        "source_study",
        "tumor_type",
        "gene",
        "peptide_mut",
        "peptide_wt",
        "hla_allele",
        label_col,
        "NMI_clean_consensus",
        "NMI_all_branch_mean",
        "NMI_tcr_qk_structure",
        "NMI_equal_top5",
        "NMI_w7b_esm2_qk",
        "NMI_branch_RF_source_oof",
        "NMI_branch_LR_source_oof",
        "NMI_core_LR_source_oof",
        "NMI_transparent_score",
        "any_existing_overlap_flag",
    ]
    for c in pred_cols:
        if c not in d.columns:
            d[c] = np.nan
    write_tsv(d[pred_cols].sort_values("NMI_clean_consensus", ascending=False), outdir / "predictions" / "nmi_clean_method_predictions.tsv")
    top34 = d.copy()
    top34["NMI_method_score"] = top34["NMI_all_branch_mean"].combine_first(top34["NMI_clean_consensus"])
    top34["nmi_rank_within_patient"] = top34.groupby("patient_id")["NMI_method_score"].rank(method="first", ascending=False)
    write_tsv(top34[top34["nmi_rank_within_patient"] <= 34][pred_cols + ["nmi_rank_within_patient"]].sort_values(["patient_id", "nmi_rank_within_patient"]), outdir / "predictions" / "patient_top34_candidates_nmi_clean.tsv")

    branch_manifest = pd.DataFrame(
        {
            "branch_feature": branch_cols,
            "n_available": [int(d[c].notna().sum()) for c in branch_cols],
            "allowed_in_nmi": True,
            "reason": "local/own branch score; public predictor tokens blocked",
        }
    )
    write_tsv(branch_manifest, outdir / "reports" / "nmi_branch_feature_manifest.tsv")
    write_figures(outdir, lb)

    best = lb[lb["model"].astype(str).str.startswith("NMI")].sort_values("AUPRC", ascending=False).head(1)
    best_line = "_No NMI metrics available._"
    if not best.empty:
        r = best.iloc[0]
        best_line = f"`{r['model']}` / `{r['split']}`: AUPRC={r['AUPRC']:.3f}, AUROC={r['AUROC']:.3f}, patient_hit@34={r.get('patient_hit_rate@34', np.nan):.3f}."
    md = [
        "# NMI Method Paper Strategy",
        "",
        "**NMI = Neoantigen Multimodal Immunogenicity.**",
        "",
        "## Decision",
        "Use NMI as the method-paper core, not the public-predictor-containing GA/RL production controller.",
        "",
        "## Why this is cleaner",
        "- No NetMHCpan, MHCflurry, BigMHC, PRIME, HLApollo, BAR-Neo, KG-GA, or product-impact scores are used as NMI training features.",
        "- NMI uses peptide/HLA/self/TCR/tumor context plus local branches only: Structure, ESM2, Wave8/TCR-self-similarity, quantum/QK, and W7 stacks.",
        "- Public predictors are reported only as frozen comparators.",
        "",
        "## Best current NMI row",
        f"- {best_line}",
        "",
        "## Leaderboard",
        lb.head(18).to_markdown(index=False),
        "",
        "## Branch feature manifest",
        branch_manifest.head(40).to_markdown(index=False) if not branch_manifest.empty else "_No branch features found._",
        "",
        "## Method-paper claim boundary",
        "- Allowed: NMI is a leakage-aware clean multimodal immunogenicity ranker evaluated retrospectively on public labels.",
        "- Allowed: NMI separates clean-method evidence from production stacking and public predictor comparators.",
        "- Forbidden: clinical vaccine efficacy, verified antigen presentation without MS, immunogenicity without T-cell assay labels, or universal superiority across all cohorts.",
        "",
        "## Paper title candidates",
        "1. NMI: a leakage-aware multimodal immunogenicity ranker for patient-level neoantigen prioritization",
        "2. CLEAN-NMI: separating TCR-visible immunogenicity from HLA presentation in neoantigen ranking",
        "3. Neoantigen Multimodal Immunogenicity modeling under leakage-controlled patient-level evaluation",
        "",
        "## Next method-paper gate",
        "Run NMI on a locked external patient-level set with pre-registered splits: leave-study-out, leave-HLA-supertype-out, and hospital-heldout when available.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "NMI_METHOD_PAPER_STRATEGY.md")
    print(outdir / "reports" / "NMI_METHOD_PAPER_STRATEGY.md")


if __name__ == "__main__":
    main()
