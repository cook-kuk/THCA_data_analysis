#!/usr/bin/env python3
"""Paper 9 Sprint 2 metabolic vulnerability prioritization.

Reframes Paper 9 as a research-use metabolic vulnerability prioritization map.
No protected data, new large downloads, manuscript edits, or clinical claims.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, spearmanr


ROOT = Path(__file__).resolve().parents[2]
P9 = ROOT / "project/results/p3_p9_full_execution/paper9"
P9S1 = ROOT / "project/results/p3_p9_full_execution/paper9_sprint1_gls"
RAW = P9 / "raw"
OUT = ROOT / "project/results/p3_p9_full_execution/paper9_sprint2_metabolic"
REPORT = ROOT / "project/reports/p3_p9_full_execution"


MODULES = {
    "GLS": {"type": "gene", "genes": ["GLS"], "class": "glutaminase"},
    "GLUD1": {"type": "gene", "genes": ["GLUD1"], "class": "glutamate_dehydrogenase"},
    "GLUD2": {"type": "gene", "genes": ["GLUD2"], "class": "glutamate_dehydrogenase"},
    "SLC1A5": {"type": "gene", "genes": ["SLC1A5"], "class": "glutamine_transport"},
    "SLC7A5": {"type": "gene", "genes": ["SLC7A5"], "class": "amino_acid_transport"},
    "ASNS": {"type": "gene", "genes": ["ASNS"], "class": "amino_acid_synthesis"},
    "GOT1": {"type": "gene", "genes": ["GOT1"], "class": "transaminase"},
    "GOT2": {"type": "gene", "genes": ["GOT2"], "class": "transaminase"},
    "GPT2": {"type": "gene", "genes": ["GPT2"], "class": "transaminase"},
    "MYC_glutamine_addiction_module": {
        "type": "module",
        "genes": ["MYC", "GLS", "SLC1A5", "SLC7A5", "ASNS", "LDHA"],
        "class": "myc_glutamine_addiction",
    },
    "LDHA_glycolysis_comparator": {
        "type": "module",
        "genes": ["LDHA", "SLC2A1", "HK2", "PFKP", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM"],
        "class": "glycolysis_comparator",
    },
    "OXPHOS_comparator": {
        "type": "module",
        "genes": ["NDUFS1", "NDUFS2", "NDUFA9", "NDUFB8", "SDHA", "UQCRC1", "COX4I1", "COX5A", "ATP5F1A", "ATP5F1B"],
        "class": "oxphos_comparator",
    },
    "NAMPT_NAD_salvage_comparator": {
        "type": "module",
        "genes": ["NAMPT", "NAPRT", "NADSYN1", "NMNAT1", "NMNAT2", "NMNAT3", "NMRK1", "NMRK2", "QPRT"],
        "class": "nad_salvage_comparator",
    },
}

DRUG_REGEX = {
    "GLS": r"glutamin|gls|bptes|cb-839|telaglenastat|compound 968",
    "GLUD1": r"glud|glutamate dehydrogenase|egcg|r162",
    "GLUD2": r"glud|glutamate dehydrogenase|egcg|r162",
    "SLC1A5": r"slc1a5|asct2|v-9302|benzylserine",
    "SLC7A5": r"slc7a5|lat1|jph203|bch",
    "ASNS": r"asns|asparaginase",
    "GOT1": r"got1|aminooxyacetate|aoa|transaminase",
    "GOT2": r"got2|aminooxyacetate|aoa|transaminase",
    "GPT2": r"gpt2|alanine transaminase|aminooxyacetate|aoa",
    "MYC_glutamine_addiction_module": r"\bmyc\b|brd4|\bjq1\b|glutamin|bptes|cb-839|telaglenastat|proteasome|bortezomib",
    "LDHA_glycolysis_comparator": r"ldha|lactate dehydrogenase|fx11|gsk2837808|glycolysis|2-deoxyglucose",
    "OXPHOS_comparator": r"oxphos|mitochondrial|complex i|phenformin|metformin|rotenone|oligomycin|antimycin",
    "NAMPT_NAD_salvage_comparator": r"nampt|fk866|ap866|gmx1778|kpt-9274|ot-82|nad",
}


def ensure() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    REPORT.mkdir(parents=True, exist_ok=True)


def clean_gene(x: str) -> str:
    return re.sub(r"\s*\(\d+\)$", "", str(x).strip()).upper()


def zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(skipna=True)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean(skipna=True)) / sd


def corr(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    tmp = pd.concat([pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")], axis=1).dropna()
    if len(tmp) < 4 or tmp.iloc[:, 0].nunique() < 2 or tmp.iloc[:, 1].nunique() < 2:
        return np.nan, np.nan, len(tmp)
    stat = spearmanr(tmp.iloc[:, 0], tmp.iloc[:, 1])
    return float(stat.statistic), float(stat.pvalue), len(tmp)


def read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "ModelID" in df.columns:
        ids = df["ModelID"].astype(str)
        meta = {"Unnamed: 0", "SequencingID", "ModelConditionID", "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel"}
        df = df.drop(columns=[c for c in meta if c in df.columns])
        df.index = ids
    else:
        first = df.columns[0]
        df = df.rename(columns={first: "ModelID"}).set_index("ModelID")
    df.columns = [clean_gene(c) for c in df.columns]
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.groupby(level=0).first()
    return df.T.groupby(level=0).mean(numeric_only=True).T


def strict_thyroid_mask(df: pd.DataFrame) -> pd.Series:
    cols = [c for c in ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode"] if c in df.columns]
    text = df[cols].fillna("").astype(str).agg(" ".join, axis=1).str.lower() if cols else pd.Series("", index=df.index)
    thyroid = text.str.contains("thyroid|thap|thfo|thpa|thme|thpd", regex=True)
    non_thyroid = text.str.contains("head and neck squamous|sarcoma, nos|hnsc|sarcnos", regex=True)
    return thyroid & ~non_thyroid


def module_values(matrix: pd.DataFrame, genes: list[str], dependency_strength: bool = False) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in matrix.columns]
    if not present:
        return pd.Series(index=matrix.index, dtype=float), []
    if len(present) == 1:
        val = pd.to_numeric(matrix[present[0]], errors="coerce")
    else:
        val = matrix[present].apply(zscore, axis=0).mean(axis=1)
    if dependency_strength:
        val = -val
    return val, present


def centered_corr(df: pd.DataFrame, x_col: str, y_col: str, group_col: str = "cancer_type") -> tuple[float, float, int]:
    tmp = df[[x_col, y_col, group_col]].copy()
    tmp[x_col] = pd.to_numeric(tmp[x_col], errors="coerce")
    tmp[y_col] = pd.to_numeric(tmp[y_col], errors="coerce")
    tmp = tmp.dropna()
    if len(tmp) < 4:
        return np.nan, np.nan, len(tmp)
    tmp["x_center"] = tmp[x_col] - tmp.groupby(group_col)[x_col].transform("mean")
    tmp["y_center"] = tmp[y_col] - tmp.groupby(group_col)[y_col].transform("mean")
    return corr(tmp["x_center"], tmp["y_center"])


def condition_map(path: Path) -> dict[str, str]:
    df = pd.read_csv(path)
    label_cols = [c for c in df.columns if re.search("name|compound|target|moa|mechanism|id", c, re.I)]
    mapping = {}
    for _, r in df.iterrows():
        label = " ".join(str(r[c]) for c in label_cols if pd.notna(r[c]))
        for id_col in ["CompoundID", "broad_id", "profile_id"]:
            if id_col in df.columns and pd.notna(r[id_col]):
                mapping[str(r[id_col])] = label
    return mapping


def drug_summary(lineage: pd.DataFrame) -> pd.DataFrame:
    base = lineage[["ModelID", "lineage_silenced_score", "strict_thyroid_model"]].copy()
    rows = []
    matrices = [
        ("GDSC1_AUC", RAW / "GDSC1AUCMatrix.csv", RAW / "GDSC1Log2ViabilityCollapsedConditions.csv", "lower_more_sensitive"),
        ("GDSC2_AUC", RAW / "GDSC2AUCMatrix.csv", RAW / "GDSC2Log2ViabilityCollapsedConditions.csv", "lower_more_sensitive"),
        ("CTRP_AUC", RAW / "CTRPAUCMatrix.csv", RAW / "CTRPLog2ViabilityCollapsedConditions.csv", "lower_more_sensitive"),
    ]
    for source, matrix_path, meta_path, direction in matrices:
        labels = condition_map(meta_path)
        mat = pd.read_csv(matrix_path).rename(columns={"Unnamed: 0": "ModelID"})
        for candidate, rx_text in DRUG_REGEX.items():
            rx = re.compile(rx_text, re.I)
            hit_cols = [k for k, v in labels.items() if k in mat.columns and rx.search(f"{k} {v}")]
            for col in hit_cols:
                joined = base.merge(mat[["ModelID", col]], on="ModelID", how="inner")
                sens = -pd.to_numeric(joined[col], errors="coerce") if direction == "lower_more_sensitive" else pd.to_numeric(joined[col], errors="coerce")
                r, p, n = corr(joined["lineage_silenced_score"], sens)
                th = joined["strict_thyroid_model"] == True
                rt, pt, nt = corr(joined.loc[th, "lineage_silenced_score"], sens[th])
                rows.append(
                    {
                        "candidate": candidate,
                        "source": source,
                        "compound_id": col,
                        "compound_label": labels.get(col, ""),
                        "spearman_lineage_sensitivity": r,
                        "p_lineage_sensitivity": p,
                        "n_models": n,
                        "strict_thyroid_n": nt,
                        "strict_thyroid_spearman": rt,
                        "strict_thyroid_p": pt,
                        "direction_note": "positive sensitivity correlation supports higher sensitivity with lineage silencing",
                    }
                )

    prism_meta = pd.read_csv(RAW / "Repurposing_Public_24Q2_Treatment_Meta_Data.csv")
    prism_text = prism_meta.astype(str).agg(" ".join, axis=1)
    for candidate, rx_text in DRUG_REGEX.items():
        hit = prism_meta[prism_text.str.contains(rx_text, case=False, regex=True, na=False)]
        if hit.empty:
            rows.append(
                {
                    "candidate": candidate,
                    "source": "PRISM24Q2_LFC_COLLAPSED",
                    "compound_id": "none_found",
                    "compound_label": "No metadata match by Sprint 2 regex",
                    "spearman_lineage_sensitivity": np.nan,
                    "p_lineage_sensitivity": np.nan,
                    "n_models": 0,
                    "strict_thyroid_n": 0,
                    "strict_thyroid_spearman": np.nan,
                    "strict_thyroid_p": np.nan,
                    "direction_note": "not testable",
                }
            )
    drug = pd.DataFrame(rows)
    drug.to_csv(OUT / "paper9_drug_concordance_summary.tsv", sep="\t", index=False)
    return drug


def main() -> None:
    ensure()
    lineage = pd.read_csv(P9 / "p9_cellline_lineage_scores.tsv", sep="\t")
    lineage["strict_thyroid_model"] = strict_thyroid_mask(lineage).values
    expr = read_matrix(RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
    dep = read_matrix(RAW / "CRISPRGeneEffect.csv")
    common = set(pd.read_csv(RAW / "AchillesCommonEssentialControls.csv").iloc[:, 0].dropna().astype(str).map(clean_gene))
    nonessential = set(pd.read_csv(RAW / "AchillesNonessentialControls.csv").iloc[:, 0].dropna().astype(str).map(clean_gene))

    common_present = sorted([g for g in common if g in dep.columns])
    nonessential_present = sorted([g for g in nonessential if g in dep.columns])
    model_base = lineage[["ModelID", "lineage_silenced_score", "cancer_type", "strict_thyroid_model"]].copy()
    model_base["common_essential_proxy"] = dep[common_present].mean(axis=1).reindex(model_base["ModelID"]).values if common_present else np.nan
    model_base["nonessential_proxy"] = dep[nonessential_present].mean(axis=1).reindex(model_base["ModelID"]).values if nonessential_present else np.nan

    assoc_rows = []
    artifact_rows = []
    model_scores = model_base.copy()
    for candidate, spec in MODULES.items():
        expr_val, expr_present = module_values(expr, spec["genes"], dependency_strength=False)
        dep_strength, dep_present = module_values(dep, spec["genes"], dependency_strength=True)
        raw_dep, raw_dep_present = module_values(dep, spec["genes"], dependency_strength=False)
        model_scores[f"{candidate}__expression_score"] = expr_val.reindex(model_scores["ModelID"]).values
        model_scores[f"{candidate}__dependency_strength"] = dep_strength.reindex(model_scores["ModelID"]).values
        model_scores[f"{candidate}__raw_gene_effect_or_module"] = raw_dep.reindex(model_scores["ModelID"]).values

        tmp = model_base.copy()
        tmp["expression_score"] = expr_val.reindex(tmp["ModelID"]).values
        tmp["dependency_strength"] = dep_strength.reindex(tmp["ModelID"]).values
        tmp["raw_dependency"] = raw_dep.reindex(tmp["ModelID"]).values

        dep_r, dep_p, dep_n = corr(tmp["lineage_silenced_score"], tmp["dependency_strength"])
        expr_r, expr_p, expr_n = corr(tmp["lineage_silenced_score"], tmp["expression_score"])
        cdep_r, cdep_p, cdep_n = centered_corr(tmp, "lineage_silenced_score", "dependency_strength")
        cexpr_r, cexpr_p, cexpr_n = centered_corr(tmp, "lineage_silenced_score", "expression_score")
        th = tmp["strict_thyroid_model"] == True
        th_dep_r, th_dep_p, th_dep_n = corr(tmp.loc[th, "lineage_silenced_score"], tmp.loc[th, "dependency_strength"])
        th_expr_r, th_expr_p, th_expr_n = corr(tmp.loc[th, "lineage_silenced_score"], tmp.loc[th, "expression_score"])
        assoc_rows.append(
            {
                "candidate": candidate,
                "candidate_type": spec["type"],
                "class": spec["class"],
                "genes_requested": ";".join(spec["genes"]),
                "dependency_genes_present": ";".join(dep_present),
                "expression_genes_present": ";".join(expr_present),
                "dependency_direction": "positive dependency_strength means stronger dependency with lineage silencing",
                "dependency_spearman": dep_r,
                "dependency_p": dep_p,
                "dependency_n": dep_n,
                "expression_spearman": expr_r,
                "expression_p": expr_p,
                "expression_n": expr_n,
                "pan_cancer_corrected_dependency_spearman": cdep_r,
                "pan_cancer_corrected_dependency_p": cdep_p,
                "pan_cancer_corrected_dependency_n": cdep_n,
                "pan_cancer_corrected_expression_spearman": cexpr_r,
                "pan_cancer_corrected_expression_p": cexpr_p,
                "strict_thyroid_dependency_spearman": th_dep_r,
                "strict_thyroid_dependency_p": th_dep_p,
                "strict_thyroid_dependency_n": th_dep_n,
                "strict_thyroid_expression_spearman": th_expr_r,
                "strict_thyroid_expression_p": th_expr_p,
                "strict_thyroid_expression_n": th_expr_n,
            }
        )

        ce_r, ce_p, ce_n = corr(tmp["dependency_strength"], tmp["common_essential_proxy"])
        ne_r, ne_p, ne_n = corr(tmp["dependency_strength"], tmp["nonessential_proxy"])
        lineage_labels = tmp.dropna(subset=["dependency_strength", "cancer_type"])
        groups = [g["dependency_strength"].values for _, g in lineage_labels.groupby("cancer_type") if len(g) >= 5]
        kw_p = kruskal(*groups).pvalue if len(groups) >= 3 else np.nan
        artifact_rows.append(
            {
                "candidate": candidate,
                "common_essential_proxy_spearman": ce_r,
                "common_essential_proxy_p": ce_p,
                "common_essential_proxy_n": ce_n,
                "nonessential_proxy_spearman": ne_r,
                "nonessential_proxy_p": ne_p,
                "nonessential_proxy_n": ne_n,
                "tissue_lineage_kruskal_p": kw_p,
                "common_essential_overlap_genes": ";".join([g for g in spec["genes"] if g in common]),
                "artifact_risk": artifact_risk(spec["genes"], ce_r, kw_p),
            }
        )

    assoc = pd.DataFrame(assoc_rows)
    assoc.to_csv(OUT / "paper9_module_dependency_associations.tsv", sep="\t", index=False)
    artifact = pd.DataFrame(artifact_rows)
    artifact.to_csv(OUT / "paper9_artifact_control_summary.tsv", sep="\t", index=False)
    model_scores.to_csv(OUT / "paper9_metabolic_model_scores.tsv", sep="\t", index=False)

    drug = drug_summary(lineage)
    ranking = build_ranking(assoc, artifact, drug)
    ranking.to_csv(OUT / "paper9_metabolic_candidate_ranking.tsv", sep="\t", index=False)
    write_docs(ranking, assoc, artifact, drug)


def artifact_risk(genes: list[str], common_proxy_r: float, tissue_p: float) -> str:
    overlap = {"MYC", "LDHA", "NAMPT"}.intersection(genes)
    if overlap or (pd.notna(common_proxy_r) and abs(common_proxy_r) >= 0.25):
        return "moderate_high_proliferation_or_common_essential_risk"
    if pd.notna(tissue_p) and tissue_p < 1e-5:
        return "moderate_tissue_context_risk"
    return "low_to_moderate"


def summarize_drug(candidate: str, drug: pd.DataFrame) -> tuple[str, int, float, float, str]:
    hit = drug[(drug["candidate"] == candidate) & (drug["compound_id"] != "none_found")].copy()
    if hit.empty:
        return "not_testable", 0, np.nan, np.nan, "no matched processed drug metadata"
    hit["abs_r"] = hit["spearman_lineage_sensitivity"].abs()
    best = hit.sort_values(["p_lineage_sensitivity", "abs_r"], ascending=[True, False], na_position="last").iloc[0]
    support = "supportive" if pd.notna(best["spearman_lineage_sensitivity"]) and best["spearman_lineage_sensitivity"] > 0 and best["p_lineage_sensitivity"] < 0.05 else "weak_or_discordant"
    label = f"{best['source']} {best['compound_label']} r={best['spearman_lineage_sensitivity']:.3g}, p={best['p_lineage_sensitivity']:.3g}"
    return support, int(len(hit)), float(best["spearman_lineage_sensitivity"]), float(best["p_lineage_sensitivity"]), label


def build_ranking(assoc: pd.DataFrame, artifact: pd.DataFrame, drug: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, a in assoc.iterrows():
        cand = a["candidate"]
        art = artifact[artifact["candidate"] == cand].iloc[0]
        drug_support, drug_n, drug_r, drug_p, drug_label = summarize_drug(cand, drug)
        dependency_support = "supportive" if pd.notna(a["dependency_spearman"]) and a["dependency_spearman"] > 0 and a["dependency_p"] < 0.05 else "weak_or_negative"
        expression_support = "supportive" if pd.notna(a["expression_spearman"]) and abs(a["expression_spearman"]) >= 0.1 and a["expression_p"] < 0.05 else "weak_or_none"
        corrected_support = "supportive" if pd.notna(a["pan_cancer_corrected_dependency_spearman"]) and a["pan_cancer_corrected_dependency_spearman"] > 0 and a["pan_cancer_corrected_dependency_p"] < 0.1 else "weak_or_none"
        thyroid_coverage = f"strict dependency n={int(a['strict_thyroid_dependency_n'])}; strict expression n={int(a['strict_thyroid_expression_n'])}; underpowered"
        risk = art["artifact_risk"]
        claim = "allowed_candidate_hypothesis"
        tier = "LOW"
        verdict = "WATCHLIST"
        if cand == "GLS" and dependency_support == "supportive":
            tier = "MEDIUM_MINUS"
            verdict = "KEEP_AS_CANDIDATE"
        elif dependency_support == "supportive" and corrected_support == "supportive" and not str(risk).startswith("moderate_high"):
            tier = "MEDIUM"
            verdict = "PRIORITIZE_FOR_WETLAB"
        elif dependency_support == "supportive":
            tier = "HYPOTHESIS_SUPPORT"
            verdict = "KEEP_AS_AXIS_OR_SECONDARY"
        if str(risk).startswith("moderate_high"):
            tier = "LOW_ARTIFACT_GUARDED" if tier == "LOW" else f"{tier}_ARTIFACT_GUARDED"
            claim = "weak_artifact_guard_required"
        if "comparator" in cand.lower():
            claim = "weak_comparator_only"
            verdict = "COMPARATOR_ONLY" if dependency_support != "supportive" else verdict
        rows.append(
            {
                "candidate": cand,
                "candidate_type": a["candidate_type"],
                "class": a["class"],
                "dependency_evidence": f"{dependency_support}; r={a['dependency_spearman']:.3g}, p={a['dependency_p']:.3g}, n={int(a['dependency_n'])}",
                "expression_evidence": f"{expression_support}; r={a['expression_spearman']:.3g}, p={a['expression_p']:.3g}, n={int(a['expression_n'])}",
                "drug_response_evidence": f"{drug_support}; matches={drug_n}; {drug_label}",
                "artifact_risk": risk,
                "thyroid_specific_coverage": thyroid_coverage,
                "pan_cancer_corrected_association": f"{corrected_support}; r={a['pan_cancer_corrected_dependency_spearman']:.3g}, p={a['pan_cancer_corrected_dependency_p']:.3g}",
                "claim_status": claim,
                "claim_boundary": "allowed: metabolic candidate prioritization; weak: thyroid-specific inference; forbidden: validated synthetic lethality or treatment recommendation",
                "evidence_tier": tier,
                "verdict": verdict,
            }
        )
    priority = {"MEDIUM": 0, "MEDIUM_MINUS": 1, "MEDIUM_ARTIFACT_GUARDED": 2, "HYPOTHESIS_SUPPORT": 3, "HYPOTHESIS_SUPPORT_ARTIFACT_GUARDED": 4, "LOW_ARTIFACT_GUARDED": 5, "LOW": 6}
    out = pd.DataFrame(rows)
    out["_priority"] = out["evidence_tier"].map(priority).fillna(9)
    out = out.sort_values(["_priority", "candidate"]).drop(columns=["_priority"])
    return out


def write_docs(ranking: pd.DataFrame, assoc: pd.DataFrame, artifact: pd.DataFrame, drug: pd.DataFrame) -> None:
    (REPORT / "paper9_sprint2_plan.md").write_text(sprint2_plan_text(ranking) + "\n")
    (REPORT / "paper9_professor_update_kr.md").write_text(professor_update_kr(ranking) + "\n")
    (REPORT / "paper9_go_no_go_criteria.md").write_text(go_no_go_text() + "\n")


def sprint2_plan_text(ranking: pd.DataFrame) -> str:
    top = ranking.head(6)
    lines = [
        "# Paper 9 Sprint 2 Metabolic Vulnerability Plan",
        "",
        "## Reframed Goal",
        "Paper 9 is framed as a metabolic vulnerability prioritization map for lineage-silenced thyroid cancer, not as validated synthetic lethality.",
        "",
        "## Current Strongest Claim",
        "Lineage-silenced thyroid cancer may exhibit glutamine-axis vulnerability candidates requiring experimental validation.",
        "",
        "## Sprint 2 Modules",
        "- Glutamine-axis genes: GLS, GLUD1, GLUD2, SLC1A5, SLC7A5, ASNS, GOT1, GOT2, GPT2.",
        "- MYC glutamine-addiction module.",
        "- LDHA/glycolysis comparator.",
        "- OXPHOS comparator.",
        "- NAMPT/NAD salvage comparator.",
        "",
        "## Analysis Tasks",
        "- Use existing DepMap/CCLE CRISPR and expression matrices only.",
        "- Score individual genes and modules for dependency strength and expression.",
        "- Estimate pan-cancer and cancer-type-centered associations with the lineage-silenced score.",
        "- Audit common-essential, nonessential, tissue-context, MYC/LDHA/NAMPT artifact risks.",
        "- Summarize processed PRISM/GDSC/CTRP drug evidence when metadata matches are available.",
        "",
        "## Current Prioritization Snapshot",
    ]
    for _, r in top.iterrows():
        lines.append(f"- {r['candidate']}: {r['verdict']} ({r['evidence_tier']}); claim={r['claim_status']}.")
    lines.extend(
        [
            "",
            "## Decision Boundary",
            "- Allowed: candidate metabolic vulnerability and research-use prioritization.",
            "- Weak: thyroid-specific claims, because strict thyroid CRISPR overlap remains small.",
            "- Forbidden: validated synthetic lethality, patient selection, clinical treatment recommendation.",
        ]
    )
    return "\n".join(lines)


def professor_update_kr(ranking: pd.DataFrame) -> str:
    gls = ranking[ranking["candidate"] == "GLS"].iloc[0]
    lines = [
        "# Paper 9 교수님 업데이트: 대사 취약성 우선순위 지도",
        "",
        "## 핵심 결론",
        "Paper 9는 GLS 단일 타깃의 synthetic lethality 논문이 아니라, lineage-silenced thyroid cancer에서 대사 취약성 후보를 우선순위화하는 지도로 재정의하는 것이 가장 안전합니다.",
        "",
        "## 현재 GLS 상태",
        f"- GLS verdict: {gls['verdict']} ({gls['evidence_tier']}).",
        "- pan-cancer dependency 방향은 가설과 맞지만 효과 크기는 작고, thyroid-only 분석은 모델 수가 부족합니다.",
        "- GDSC BPTES drug concordance는 현재 지지적이지 않습니다.",
        "",
        "## 확장 방향",
        "- GLS 단독이 아니라 glutamine transport/metabolism module, MYC glutamine-addiction module, glycolysis/OXPHOS/NAD salvage comparator를 함께 평가합니다.",
        "- 이 프레임은 negative/weak drug evidence를 숨기지 않고, wet-lab 검증 후보를 선별하는 구조입니다.",
        "",
        "## 사용할 수 있는 표현",
        "- lineage-silenced thyroid cancer may exhibit glutamine-axis vulnerability candidates requiring experimental validation.",
        "- metabolic vulnerability prioritization map.",
        "",
        "## 피해야 할 표현",
        "- validated synthetic lethality.",
        "- clinical treatment recommendation.",
        "- patient selection biomarker.",
    ]
    return "\n".join(lines)


def go_no_go_text() -> str:
    return "\n".join(
        [
            "# Paper 9 Go/No-Go Criteria",
            "",
            "## Go",
            "- At least one glutamine-axis gene/module has supportive pan-cancer dependency association after cancer-type centering.",
            "- Artifact risk is low-to-moderate or explicitly controllable with common-essential/proliferation controls.",
            "- Strict thyroid model coverage is sufficient for a descriptive audit, even if underpowered for inference.",
            "- Wet-lab assays are feasible for DM1-like and non-DM1-like thyroid models.",
            "",
            "## Conditional Go",
            "- Evidence is pan-cancer supportive but thyroid-only underpowered.",
            "- Drug concordance is weak or absent, but dependency/expression evidence supports a mechanistic wet-lab screen.",
            "- Module-level signal is stronger than any single target.",
            "",
            "## No-Go / Demote",
            "- Signal disappears after cancer-type centering.",
            "- Candidate tracks with common-essential/proliferation proxies.",
            "- Drug concordance is consistently opposite and dependency evidence is weak.",
            "- No feasible thyroid wet-lab model exists.",
            "",
            "## Claim Boundary",
            "- Go does not mean clinical actionability.",
            "- Go does not mean validated synthetic lethality.",
            "- Go means candidate prioritization for experimental validation.",
        ]
    )


if __name__ == "__main__":
    main()
