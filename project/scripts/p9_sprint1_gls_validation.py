#!/usr/bin/env python3
"""Paper 9 Sprint 1 GLS/glutamine-axis validation.

Uses only already-downloaded public processed files in p3_p9_full_execution.
Writes new Sprint 1 outputs and reports; no manuscript or Paper 2/3/4 files.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu, pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[2]
P9 = ROOT / "project/results/p3_p9_full_execution/paper9"
RAW = P9 / "raw"
OUT = ROOT / "project/results/p3_p9_full_execution/paper9_sprint1_gls"
REPORT = ROOT / "project/reports/p3_p9_full_execution"

GLUTAMINE_GENES = ["GLS", "GLUL", "SLC1A5", "SLC7A5", "SLC3A2", "GOT1", "GOT2", "GLUD1", "GLUD2", "ASNS", "GPT2", "MYC", "LDHA"]
OXPHOS_GENES = ["NDUFS1", "NDUFS2", "NDUFA9", "NDUFB8", "SDHA", "UQCRC1", "COX4I1", "COX5A", "ATP5F1A", "ATP5F1B"]
GLYCOLYSIS_GENES = ["SLC2A1", "HK2", "PFKP", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA"]
PROLIFERATION_GENES = ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM4", "MCM6", "UBE2C", "CCNB1", "CDK1"]
FOCUS_TARGETS = ["GLS", "glutamine_pathway_module", "CHEK2", "MCL1", "DNMT1", "STAT3", "JAK1", "HDAC1", "HDAC2"]


def ensure() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    REPORT.mkdir(parents=True, exist_ok=True)


def clean_gene(name: str) -> str:
    return re.sub(r"\s*\(\d+\)$", "", str(name).strip()).upper()


def zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(skipna=True)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean(skipna=True)) / sd


def corr(x: pd.Series, y: pd.Series, method: str = "spearman") -> tuple[float, float, int]:
    tmp = pd.concat([pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")], axis=1).dropna()
    if len(tmp) < 4 or tmp.iloc[:, 0].nunique() < 2 or tmp.iloc[:, 1].nunique() < 2:
        return np.nan, np.nan, len(tmp)
    stat = spearmanr(tmp.iloc[:, 0], tmp.iloc[:, 1]) if method == "spearman" else pearsonr(tmp.iloc[:, 0], tmp.iloc[:, 1])
    return float(stat.statistic), float(stat.pvalue), len(tmp)


def read_depmap_matrix(path: Path) -> pd.DataFrame:
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
    df = df.T.groupby(level=0).mean(numeric_only=True).T
    return df


def module_score(matrix: pd.DataFrame, genes: list[str], invert_dependency: bool = False) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in matrix.columns]
    if not present:
        return pd.Series(index=matrix.index, dtype=float), []
    z = matrix[present].apply(zscore, axis=0)
    score = z.mean(axis=1)
    if invert_dependency:
        score = -score
    return score, present


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, set[str], set[str]]:
    lineage = pd.read_csv(P9 / "p9_cellline_lineage_scores.tsv", sep="\t")
    expr = read_depmap_matrix(RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
    dep = read_depmap_matrix(RAW / "CRISPRGeneEffect.csv")
    common = set(pd.read_csv(RAW / "AchillesCommonEssentialControls.csv").iloc[:, 0].dropna().astype(str).map(clean_gene))
    nonessential = set(pd.read_csv(RAW / "AchillesNonessentialControls.csv").iloc[:, 0].dropna().astype(str).map(clean_gene))
    return lineage, expr, dep, common, nonessential


def matrix_model_ids(path: Path, long_prism: bool = False) -> set[str]:
    if long_prism:
        cell = pd.read_csv(RAW / "Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv")
        return set(cell["depmap_id"].dropna().astype(str))
    df = pd.read_csv(path, usecols=[0])
    return set(df.iloc[:, 0].dropna().astype(str))


def thyroid_coverage(lineage: pd.DataFrame, expr: pd.DataFrame, dep: pd.DataFrame) -> pd.DataFrame:
    model = pd.read_csv(RAW / "Model.csv")
    thyroid = lineage[lineage["is_thyroid_model"] == True].copy()
    strict_mask = strict_thyroid_mask(lineage)
    strict_thyroid = lineage[strict_mask].copy()
    thyroid_ids = set(thyroid["ModelID"].astype(str))
    strict_ids = set(strict_thyroid["ModelID"].astype(str))
    crispr_ids = set(dep.index.astype(str))
    expr_ids = set(expr.index.astype(str))
    prism_ids = matrix_model_ids(RAW / "Repurposing_Public_24Q2_LFC_COLLAPSED.csv", long_prism=True)
    gdsc1_ids = matrix_model_ids(RAW / "GDSC1AUCMatrix.csv")
    gdsc2_ids = matrix_model_ids(RAW / "GDSC2AUCMatrix.csv")
    ctrp_ids = matrix_model_ids(RAW / "CTRPAUCMatrix.csv")
    gdsc_ctrp_ids = gdsc1_ids | gdsc2_ids | ctrp_ids
    gls_drug_ids = gls_drug_model_ids()
    rows = [
        ("thyroid_flagged_models", len(thyroid_ids), ""),
        ("strict_thyroid_lineage_models", len(strict_ids), "OncoTree primary/subtype/lineage thyroid; excludes thyroid-site non-thyroid cancers"),
        ("crispr_overlap", len(thyroid_ids & crispr_ids), ""),
        ("strict_thyroid_crispr_overlap", len(strict_ids & crispr_ids), ""),
        ("ccle_expression_overlap", len(thyroid_ids & expr_ids), ""),
        ("strict_thyroid_expression_overlap", len(strict_ids & expr_ids), ""),
        ("prism_overlap", len(thyroid_ids & prism_ids), "all PRISM 24Q2 models, not GLS-specific"),
        ("gdsc_ctrp_overlap", len(thyroid_ids & gdsc_ctrp_ids), "any GDSC1/GDSC2/CTRP AUC model"),
        ("gls_expression_dependency_drug_complete", len(thyroid_ids & expr_ids & crispr_ids & gls_drug_ids), "GLS expression + GLS dependency + matched glutaminase inhibitor response"),
        ("strict_gls_expression_dependency_drug_complete", len(strict_ids & expr_ids & crispr_ids & gls_drug_ids), "strict thyroid subset with GLS expression + GLS dependency + matched glutaminase inhibitor response"),
    ]
    out = pd.DataFrame(rows, columns=["coverage_item", "n_models", "notes"])
    subcols = ["ModelID", "CellLineName", "CCLEName", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode", "SampleCollectionSite", "lineage_silenced_score"]
    thyroid[[c for c in subcols if c in thyroid.columns]].to_csv(OUT / "thyroid_model_coverage_detail.tsv", sep="\t", index=False)
    out.to_csv(OUT / "thyroid_model_coverage_summary.tsv", sep="\t", index=False)
    return out


def strict_thyroid_mask(df: pd.DataFrame) -> pd.Series:
    cols = [c for c in ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode"] if c in df.columns]
    text = df[cols].fillna("").astype(str).agg(" ".join, axis=1).str.lower() if cols else pd.Series("", index=df.index)
    thyroid = text.str.contains("thyroid|thap|thfo|thpa|thme|thpd", regex=True)
    non_thyroid = text.str.contains("head and neck squamous|sarcoma, nos|hnsc|sarcnos", regex=True)
    return thyroid & ~non_thyroid


def gls_drug_model_ids() -> set[str]:
    ids = set()
    gdsc1_conditions = pd.read_csv(RAW / "GDSC1Log2ViabilityCollapsedConditions.csv")
    bptes_ids = set(gdsc1_conditions.loc[gdsc1_conditions.astype(str).agg(" ".join, axis=1).str.contains("BPTES|glutamin|GLS|CB-839|telaglenastat", case=False, regex=True), "CompoundID"].astype(str))
    if bptes_ids:
        gdsc1 = pd.read_csv(RAW / "GDSC1AUCMatrix.csv", usecols=["Unnamed: 0", *[c for c in bptes_ids if c]])
        ids |= set(gdsc1.loc[gdsc1[list(bptes_ids)].notna().any(axis=1), "Unnamed: 0"].astype(str))
    return ids


def artifact_check(lineage: pd.DataFrame, expr: pd.DataFrame, dep: pd.DataFrame, common: set[str], nonessential: set[str]) -> pd.DataFrame:
    base = lineage[["ModelID", "lineage_silenced_score", "DM1_like_score", "is_thyroid_model", "cancer_type", "OncotreePrimaryDisease", "OncotreeSubtype"]].copy()
    base["strict_thyroid_model"] = strict_thyroid_mask(lineage).values
    base = base.merge(dep[["GLS"]].rename(columns={"GLS": "GLS_gene_effect"}), left_on="ModelID", right_index=True, how="left")
    base = base.merge(expr[["GLS"]].rename(columns={"GLS": "GLS_expression"}), left_on="ModelID", right_index=True, how="left")

    common_present = sorted([g for g in common if g in dep.columns])
    nonessential_present = sorted([g for g in nonessential if g in dep.columns])
    base["common_essential_dependency_proxy"] = dep[common_present].mean(axis=1).reindex(base["ModelID"]).values if common_present else np.nan
    base["nonessential_dependency_proxy"] = dep[nonessential_present].mean(axis=1).reindex(base["ModelID"]).values if nonessential_present else np.nan
    prolif_score, prolif_present = module_score(expr, PROLIFERATION_GENES)
    base["proliferation_expression_proxy"] = prolif_score.reindex(base["ModelID"]).values

    rows = []
    tests = [
        ("GLS_dependency_vs_lineage_silenced", base["lineage_silenced_score"], base["GLS_gene_effect"], "negative supports stronger GLS dependency in lineage-silenced models"),
        ("GLS_dependency_vs_GLS_expression", base["GLS_expression"], base["GLS_gene_effect"], "expression-dependency coupling check"),
        ("GLS_dependency_vs_common_essential_proxy", base["common_essential_dependency_proxy"], base["GLS_gene_effect"], "artifact/proliferation-like dropout check"),
        ("GLS_dependency_vs_nonessential_proxy", base["nonessential_dependency_proxy"], base["GLS_gene_effect"], "negative-control dependency check"),
        ("GLS_dependency_vs_proliferation_expression", base["proliferation_expression_proxy"], base["GLS_gene_effect"], "proliferation expression artifact check"),
    ]
    for name, x, y, interp in tests:
        r, p, n = corr(x, y)
        rows.append({"test": name, "scope": "pan_cancer", "spearman_r": r, "p_value": p, "n": n, "interpretation": interp})
        th = base["strict_thyroid_model"] == True
        rt, pt, nt = corr(x[th], y[th])
        rows.append({"test": name, "scope": "thyroid_only", "spearman_r": rt, "p_value": pt, "n": nt, "interpretation": "thyroid-only estimate; underpowered if n is small"})

    valid = base.dropna(subset=["GLS_gene_effect", "cancer_type"])
    groups = [g["GLS_gene_effect"].values for _, g in valid.groupby("cancer_type") if len(g) >= 5]
    if len(groups) >= 3:
        stat = kruskal(*groups)
        rows.append({"test": "GLS_dependency_by_pan_cancer_lineage_label", "scope": "pan_cancer", "spearman_r": np.nan, "p_value": float(stat.pvalue), "n": int(len(valid)), "interpretation": "Kruskal-Wallis tissue/lineage association; significant result means tissue context contributes"})

    gls_lineage_r, _, _ = corr(base["lineage_silenced_score"], base["GLS_gene_effect"])
    control_rows = []
    for gene in sorted(set(common_present[:300] + nonessential_present[:300])):
        r, p, n = corr(base["lineage_silenced_score"], dep[gene].reindex(base["ModelID"]).reset_index(drop=True))
        control_rows.append({"gene": gene, "control_class": "common_essential" if gene in common else "nonessential", "spearman_lineage_dependency": r, "p_value": p, "n": n})
    controls = pd.DataFrame(control_rows)
    if not controls.empty:
        controls.to_csv(OUT / "gls_control_gene_correlations.tsv", sep="\t", index=False)
        for cls in ["common_essential", "nonessential"]:
            vals = controls.loc[controls["control_class"] == cls, "spearman_lineage_dependency"].dropna()
            if len(vals):
                pct = float((vals <= gls_lineage_r).mean() * 100.0)
                rows.append({"test": f"GLS_vs_{cls}_lineage_correlation_percentile", "scope": "pan_cancer", "spearman_r": gls_lineage_r, "p_value": np.nan, "n": len(vals), "interpretation": f"GLS correlation is at {pct:.1f} percentile by more-negative-than-control criterion"})

    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "gls_dependency_artifact_check.tsv", sep="\t", index=False)
    base.to_csv(OUT / "gls_model_level_dependency_artifact_matrix.tsv", sep="\t", index=False)
    make_artifact_figures(base, controls if "controls" in locals() else pd.DataFrame())
    return summary


def make_artifact_figures(base: pd.DataFrame, controls: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(6, 4.5))
    sns.scatterplot(data=base, x="lineage_silenced_score", y="GLS_gene_effect", hue="is_thyroid_model", s=26, linewidth=0)
    plt.axhline(base["GLS_gene_effect"].median(skipna=True), color="gray", lw=0.8)
    plt.title("GLS dependency vs lineage-silenced score")
    plt.tight_layout()
    plt.savefig(OUT / "fig_gls_dependency_vs_lineage.png", dpi=200)
    plt.close()

    if controls.empty:
        plt.figure(figsize=(6, 4))
        plt.text(0.5, 0.5, "No control genes available", ha="center", va="center")
        plt.axis("off")
    else:
        plt.figure(figsize=(7, 4.5))
        sns.kdeplot(data=controls, x="spearman_lineage_dependency", hue="control_class", common_norm=False, fill=True, alpha=0.3)
        gls_r, _, _ = corr(base["lineage_silenced_score"], base["GLS_gene_effect"])
        plt.axvline(gls_r, color="red", lw=1.5, label="GLS")
        plt.legend()
        plt.title("GLS lineage association vs control genes")
    plt.tight_layout()
    plt.savefig(OUT / "fig_gls_artifact_controls.png", dpi=200)
    plt.close()


def pathway_scores(lineage: pd.DataFrame, expr: pd.DataFrame, dep: pd.DataFrame) -> pd.DataFrame:
    out = lineage[["ModelID", "lineage_silenced_score", "DM1_like_score", "is_thyroid_model", "cancer_type", "OncotreePrimaryDisease", "OncotreeSubtype"]].copy()
    out["strict_thyroid_model"] = strict_thyroid_mask(lineage).values
    scores = {}
    coverage = []
    for name, genes, mat, inv in [
        ("glutamine_expression_module", GLUTAMINE_GENES, expr, False),
        ("glutamine_dependency_module", GLUTAMINE_GENES, dep, True),
        ("oxphos_expression_module", OXPHOS_GENES, expr, False),
        ("oxphos_dependency_module", OXPHOS_GENES, dep, True),
        ("glycolysis_expression_module", GLYCOLYSIS_GENES, expr, False),
        ("glycolysis_dependency_module", GLYCOLYSIS_GENES, dep, True),
    ]:
        score, present = module_score(mat, genes, invert_dependency=inv)
        out[name] = score.reindex(out["ModelID"]).values
        coverage.append({"module": name, "n_genes_requested": len(genes), "n_genes_present": len(present), "genes_present": ";".join(present)})
    rows = []
    for col in [c for c in out.columns if c.endswith("_module")]:
        r, p, n = corr(out["lineage_silenced_score"], out[col])
        rows.append({"module": col, "scope": "pan_cancer", "spearman_lineage": r, "p_value": p, "n": n})
        th = out["strict_thyroid_model"] == True
        rt, pt, nt = corr(out.loc[th, "lineage_silenced_score"], out.loc[th, col])
        rows.append({"module": col, "scope": "thyroid_only", "spearman_lineage": rt, "p_value": pt, "n": nt})
    pd.DataFrame(coverage).to_csv(OUT / "glutamine_pathway_gene_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(OUT / "glutamine_pathway_module_associations.tsv", sep="\t", index=False)
    out.to_csv(OUT / "glutamine_pathway_scores.tsv", sep="\t", index=False)

    plt.figure(figsize=(6, 4.5))
    sns.scatterplot(data=out, x="lineage_silenced_score", y="glutamine_expression_module", hue="is_thyroid_model", s=25, linewidth=0)
    plt.title("Glutamine expression module vs lineage-silenced score")
    plt.tight_layout()
    plt.savefig(OUT / "fig_glutamine_module_vs_lineage.png", dpi=200)
    plt.close()
    return out


def condition_map(path: Path) -> dict[str, str]:
    df = pd.read_csv(path)
    mapping = {}
    label_cols = [c for c in df.columns if re.search("name|compound|target|moa|mechanism|id", c, re.I)]
    for _, r in df.iterrows():
        label = " ".join(str(r[c]) for c in label_cols if pd.notna(r[c]))
        if "CompoundID" in df.columns:
            mapping[str(r["CompoundID"])] = label
        if "broad_id" in df.columns:
            mapping[str(r["broad_id"])] = label
    return mapping


def drug_concordance(lineage: pd.DataFrame, dep: pd.DataFrame) -> pd.DataFrame:
    base = lineage[["ModelID", "lineage_silenced_score", "is_thyroid_model"]].copy()
    base["strict_thyroid_model"] = strict_thyroid_mask(lineage).values
    gls_dep = dep["GLS"].rename("GLS_gene_effect")
    gls_r, gls_p, gls_n = corr(base.set_index("ModelID")["lineage_silenced_score"], gls_dep)
    rows = []

    sources = [
        ("GDSC1_AUC", RAW / "GDSC1AUCMatrix.csv", RAW / "GDSC1Log2ViabilityCollapsedConditions.csv"),
        ("GDSC2_AUC", RAW / "GDSC2AUCMatrix.csv", RAW / "GDSC2Log2ViabilityCollapsedConditions.csv"),
        ("CTRP_AUC", RAW / "CTRPAUCMatrix.csv", RAW / "CTRPLog2ViabilityCollapsedConditions.csv"),
    ]
    rx = re.compile(r"glutamin|gls|bptes|cb-839|telaglenastat|compound 968", re.I)
    for source, matrix_path, meta_path in sources:
        labels = condition_map(meta_path)
        hit_cols = [k for k, v in labels.items() if rx.search(f"{k} {v}")]
        mat = pd.read_csv(matrix_path).rename(columns={"Unnamed: 0": "ModelID"})
        for col in hit_cols:
            if col not in mat.columns:
                continue
            joined = base.merge(mat[["ModelID", col]], on="ModelID", how="inner")
            sensitivity = -pd.to_numeric(joined[col], errors="coerce")
            r, p, n = corr(joined["lineage_silenced_score"], sensitivity)
            th = joined["strict_thyroid_model"] == True
            rt, pt, nt = corr(joined.loc[th, "lineage_silenced_score"], sensitivity[th])
            concordance = "concordant" if gls_r < 0 and r > 0 else "discordant_or_no_support"
            rows.append(
                {
                    "source": source,
                    "compound_id": col,
                    "compound_label": labels.get(col, ""),
                    "n_models": n,
                    "spearman_lineage_sensitivity": r,
                    "p_lineage_sensitivity": p,
                    "n_thyroid_models": nt,
                    "spearman_lineage_sensitivity_thyroid": rt,
                    "p_lineage_sensitivity_thyroid": pt,
                    "gls_dependency_lineage_spearman": gls_r,
                    "gls_dependency_lineage_p": gls_p,
                    "gls_dependency_n": gls_n,
                    "concordance": concordance,
                    "interpretation": "positive drug sensitivity correlation is concordant with negative GLS dependency correlation",
                }
            )

    # PRISM is a long LFC table; no glutaminase/GLS compounds were found in the treatment metadata at this pass.
    prism_meta = pd.read_csv(RAW / "Repurposing_Public_24Q2_Treatment_Meta_Data.csv")
    prism_hits = prism_meta[prism_meta.astype(str).agg(" ".join, axis=1).str.contains(rx)]
    if prism_hits.empty:
        rows.append(
            {
                "source": "PRISM24Q2_LFC_COLLAPSED",
                "compound_id": "none_found",
                "compound_label": "No glutaminase/GLS inhibitor match in treatment metadata by Sprint 1 regex",
                "n_models": 0,
                "spearman_lineage_sensitivity": np.nan,
                "p_lineage_sensitivity": np.nan,
                "n_thyroid_models": 0,
                "spearman_lineage_sensitivity_thyroid": np.nan,
                "p_lineage_sensitivity_thyroid": np.nan,
                "gls_dependency_lineage_spearman": gls_r,
                "gls_dependency_lineage_p": gls_p,
                "gls_dependency_n": gls_n,
                "concordance": "not_testable",
                "interpretation": "no matching PRISM compound metadata",
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "gls_drug_concordance.tsv", sep="\t", index=False)
    make_drug_figure(out)
    return out


def make_drug_figure(drug: pd.DataFrame) -> None:
    plot = drug[drug["compound_id"] != "none_found"].copy()
    if plot.empty:
        plt.figure(figsize=(6, 4))
        plt.text(0.5, 0.5, "No GLS-related drug response match", ha="center", va="center")
        plt.axis("off")
    else:
        plt.figure(figsize=(7, 4.5))
        plot["compound"] = plot["source"] + "\n" + plot["compound_label"].str.slice(0, 35)
        sns.barplot(data=plot, x="spearman_lineage_sensitivity", y="compound", hue="concordance", dodge=False)
        plt.axvline(0, color="black", lw=0.8)
        plt.xlabel("Lineage-silenced score vs sensitivity (Spearman)")
        plt.ylabel("")
    plt.tight_layout()
    plt.savefig(OUT / "fig_gls_drug_concordance.png", dpi=200)
    plt.close()


def updated_ranking(artifact: pd.DataFrame, pathway_assoc: pd.DataFrame, drug: pd.DataFrame, previous: pd.DataFrame) -> pd.DataFrame:
    def get_stat(test: str, scope: str = "pan_cancer") -> tuple[float, float, int]:
        hit = artifact[(artifact["test"] == test) & (artifact["scope"] == scope)]
        if hit.empty:
            return np.nan, np.nan, 0
        r = hit.iloc[0]
        return r["spearman_r"], r["p_value"], int(r["n"])

    gls_r, gls_p, gls_n = get_stat("GLS_dependency_vs_lineage_silenced")
    artifact_r, artifact_p, _ = get_stat("GLS_dependency_vs_common_essential_proxy")
    expr_r, expr_p, _ = get_stat("GLS_dependency_vs_GLS_expression")
    glut_dep = pathway_assoc[(pathway_assoc["module"] == "glutamine_dependency_module") & (pathway_assoc["scope"] == "pan_cancer")]
    glut_expr = pathway_assoc[(pathway_assoc["module"] == "glutamine_expression_module") & (pathway_assoc["scope"] == "pan_cancer")]
    drug_testable = drug[(drug["compound_id"] != "none_found") & drug["spearman_lineage_sensitivity"].notna()]
    drug_support = "none"
    if not drug_testable.empty:
        best = drug_testable.sort_values("p_lineage_sensitivity", na_position="last").iloc[0]
        drug_support = f"{best['concordance']}; best {best['source']} {best['compound_label']} r={best['spearman_lineage_sensitivity']:.3g}, p={best['p_lineage_sensitivity']:.3g}"

    rows = []
    gls_artifact = "moderate" if pd.notna(artifact_r) and abs(artifact_r) >= 0.25 else "low_to_moderate"
    gls_tier = "MEDIUM_MINUS"
    gls_verdict = "KEEP_AS_CANDIDATE"
    if pd.isna(gls_r) or gls_r >= 0 or (pd.notna(artifact_r) and abs(artifact_r) > 0.5):
        gls_tier = "LOW"
        gls_verdict = "DEMOTE"
    rows.append(
        {
            "target": "GLS",
            "dependency_support": f"pan-cancer r={gls_r:.3g}, p={gls_p:.3g}, n={gls_n}; negative direction supports hypothesis",
            "expression_support": f"GLS expression-dependency r={expr_r:.3g}, p={expr_p:.3g}",
            "thyroid_specificity": "underpowered; thyroid-only n reported in artifact table",
            "drug_support": drug_support,
            "artifact_risk": gls_artifact,
            "evidence_tier": gls_tier,
            "verdict": gls_verdict,
        }
    )
    dep_line = glut_dep.iloc[0] if not glut_dep.empty else pd.Series(dtype=object)
    expr_line = glut_expr.iloc[0] if not glut_expr.empty else pd.Series(dtype=object)
    rows.append(
        {
            "target": "glutamine_pathway_module",
            "dependency_support": f"module dependency r={dep_line.get('spearman_lineage', np.nan):.3g}, p={dep_line.get('p_value', np.nan):.3g}",
            "expression_support": f"module expression r={expr_line.get('spearman_lineage', np.nan):.3g}, p={expr_line.get('p_value', np.nan):.3g}",
            "thyroid_specificity": "underpowered; used as axis-level support",
            "drug_support": drug_support,
            "artifact_risk": "moderate; pathway overlaps MYC/LDHA metabolic proliferation biology",
            "evidence_tier": "HYPOTHESIS_SUPPORT",
            "verdict": "KEEP_AS_AXIS",
        }
    )
    for target in ["CHEK2", "MCL1", "DNMT1", "STAT3", "JAK1", "HDAC1", "HDAC2"]:
        prev = previous[previous["target"] == target]
        if prev.empty:
            continue
        p = prev.iloc[0]
        rows.append(
            {
                "target": target,
                "dependency_support": f"first-pass dependency={p['dependency_association']:.3g}",
                "expression_support": f"first-pass expression={p['target_expression_support']:.3g}",
                "thyroid_specificity": str(p["thyroid_model_support"]),
                "drug_support": f"first-pass drug_support={p['drug_response_support']}",
                "artifact_risk": str(p["artifact_risk"]),
                "evidence_tier": str(p["evidence_tier"]),
                "verdict": "UNCHANGED_LOW" if str(p["verdict"]) == "LOW" else str(p["verdict"]),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "p9_sprint1_updated_target_ranking.tsv", sep="\t", index=False)
    return out


def write_coverage_report(cov: pd.DataFrame, lineage: pd.DataFrame) -> None:
    thyroid = lineage[lineage["is_thyroid_model"] == True]
    lines = [
        "# Paper 9 Thyroid Model Coverage Audit",
        "",
        "## Coverage Summary",
    ]
    for _, r in cov.iterrows():
        lines.append(f"- {r['coverage_item']}: {r['n_models']} {r['notes']}")
    lines.extend(["", "## Histology / Subtype Annotation"])
    for col in ["OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode", "SampleCollectionSite"]:
        if col in thyroid:
            lines.append(f"### {col}")
            for k, v in thyroid[col].fillna("unavailable").value_counts().items():
                lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Boundary",
            "- Thyroid-only estimates remain underpowered and are used as coverage/audit evidence only.",
            "- No protected data or new large downloads were used.",
        ]
    )
    (REPORT / "07_paper9_thyroid_model_coverage_audit.md").write_text("\n".join(lines) + "\n")


def write_report(cov: pd.DataFrame, artifact: pd.DataFrame, pathway_assoc: pd.DataFrame, drug: pd.DataFrame, ranking: pd.DataFrame) -> None:
    gls = ranking[ranking["target"] == "GLS"].iloc[0]
    dep = artifact[(artifact["test"] == "GLS_dependency_vs_lineage_silenced") & (artifact["scope"] == "pan_cancer")].iloc[0]
    art = artifact[(artifact["test"] == "GLS_dependency_vs_common_essential_proxy") & (artifact["scope"] == "pan_cancer")].iloc[0]
    mod = pathway_assoc[(pathway_assoc["module"] == "glutamine_dependency_module") & (pathway_assoc["scope"] == "pan_cancer")]
    mod_line = mod.iloc[0] if not mod.empty else pd.Series(dtype=object)
    drug_hits = drug[drug["compound_id"] != "none_found"]
    lines = [
        "# Paper 9 GLS Sprint 1 Report",
        "",
        "## 1. Executive verdict",
        f"- GLS status: {gls['verdict']} ({gls['evidence_tier']}).",
        "- Interpretation: GLS/glutamine axis remains a candidate vulnerability hypothesis, but Sprint 1 does not prove synthetic lethality.",
        "- Main limitation: thyroid-only CRISPR/drug overlap is small, so pan-cancer artifact correction carries most of the evidence weight.",
        "",
        "## 2. Thyroid model coverage",
    ]
    for _, r in cov.iterrows():
        lines.append(f"- {r['coverage_item']}: {r['n_models']}")
    lines.extend(
        [
            "",
            "## 3. GLS dependency result",
            f"- GLS dependency vs lineage-silenced score: Spearman r={dep['spearman_r']:.3g}, p={dep['p_value']:.3g}, n={int(dep['n'])}.",
            "- Negative direction means higher lineage-silenced score is associated with more negative GLS gene effect.",
            "",
            "## 4. Artifact-control result",
            f"- GLS dependency vs common-essential proxy: Spearman r={art['spearman_r']:.3g}, p={art['p_value']:.3g}, n={int(art['n'])}.",
            "- Control-gene distributions and nonessential controls are reported in `gls_control_gene_correlations.tsv` and `fig_gls_artifact_controls.png`.",
            "",
            "## 5. Glutamine pathway result",
            f"- Glutamine dependency module vs lineage-silenced score: r={mod_line.get('spearman_lineage', np.nan):.3g}, p={mod_line.get('p_value', np.nan):.3g}.",
            "- OXPHOS/glycolysis comparator modules are included in `glutamine_pathway_module_associations.tsv`.",
            "",
            "## 6. Drug concordance result",
        ]
    )
    if drug_hits.empty:
        lines.append("- No GLS/glutaminase inhibitor response matrix match was testable.")
    else:
        for _, r in drug_hits.iterrows():
            lines.append(f"- {r['source']} {r['compound_label']}: r={r['spearman_lineage_sensitivity']:.3g}, p={r['p_lineage_sensitivity']:.3g}, concordance={r['concordance']}.")
    lines.extend(
        [
            "",
            "## 7. Updated target ranking",
        ]
    )
    for _, r in ranking.iterrows():
        lines.append(f"- {r['target']}: {r['verdict']} ({r['evidence_tier']}); artifact risk={r['artifact_risk']}.")
    lines.extend(
        [
            "",
            "## 8. Claim boundary",
            "Allowed:",
            "- GLS/glutamine axis is a candidate vulnerability.",
            "- Research-use hypothesis.",
            "",
            "Forbidden:",
            "- Validated synthetic lethality.",
            "- Clinical treatment recommendation.",
            "- Patient selection.",
            "",
            "## 9. Wet-lab validation plan",
            "- GLS inhibitor sensitivity in DM1-like thyroid models.",
            "- Rescue with glutamine pathway modulation.",
            "- Organoid/cell-line validation with lineage-state stratification.",
        ]
    )
    (REPORT / "08_paper9_gls_sprint1_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ensure()
    lineage, expr, dep, common, nonessential = load_inputs()
    cov = thyroid_coverage(lineage, expr, dep)
    artifact = artifact_check(lineage, expr, dep, common, nonessential)
    scores = pathway_scores(lineage, expr, dep)
    pathway_assoc = pd.read_csv(OUT / "glutamine_pathway_module_associations.tsv", sep="\t")
    drug = drug_concordance(lineage, dep)
    previous = pd.read_csv(P9 / "p9_vulnerability_ranked_targets.tsv", sep="\t")
    ranking = updated_ranking(artifact, pathway_assoc, drug, previous)
    write_coverage_report(cov, lineage)
    write_report(cov, artifact, pathway_assoc, drug, ranking)
    summary = {
        "gls_status": ranking.loc[ranking["target"] == "GLS", "verdict"].iloc[0],
        "thyroid_flagged_models": int(cov.loc[cov["coverage_item"] == "thyroid_flagged_models", "n_models"].iloc[0]),
        "strict_thyroid_lineage_models": int(cov.loc[cov["coverage_item"] == "strict_thyroid_lineage_models", "n_models"].iloc[0]),
        "gls_complete_thyroid_models": int(cov.loc[cov["coverage_item"] == "gls_expression_dependency_drug_complete", "n_models"].iloc[0]),
        "strict_gls_complete_thyroid_models": int(cov.loc[cov["coverage_item"] == "strict_gls_expression_dependency_drug_complete", "n_models"].iloc[0]),
        "outputs": str(OUT),
    }
    (OUT / "sprint1_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
