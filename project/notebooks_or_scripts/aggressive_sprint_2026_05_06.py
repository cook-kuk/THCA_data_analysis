#!/usr/bin/env python3
"""
Rapid feasibility tests for aggressive BRAF/RAS-negative thyroid cancer paper ideas.

Outputs only into project/results/aggressive_sprint_2026_05_06/.
Does not edit Paper 1 or manuscript assets.
"""

from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path

import GEOparse
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/aggressive_sprint_2026_05_06"
GEO_DIR = OUT / "geo"
OUT.mkdir(parents=True, exist_ok=True)
GEO_DIR.mkdir(parents=True, exist_ok=True)


RAI6 = ["SLC5A5", "TPO", "TSHR", "TG", "PAX8", "NKX2-1"]
RAI8 = ["SLC5A5", "TPO", "TSHR", "TG", "PAX8", "NKX2-1", "FOXE1", "DIO1"]


def cohen_d(a: pd.Series, b: pd.Series) -> float:
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    s = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / s) if s else float("nan")


def mw_p(a: pd.Series, b: pd.Series) -> float:
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 1 or len(b) < 1:
        return float("nan")
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def auc_pair(y: pd.Series, score: pd.Series) -> tuple[float, float]:
    df = pd.DataFrame({"y": y, "score": score}).dropna()
    if df["y"].nunique() != 2:
        return float("nan"), float("nan")
    raw = float(roc_auc_score(df["y"], df["score"]))
    return raw, max(raw, 1.0 - raw)


def find_exact_gene_probes(gpl: GEOparse.GEOTypes.GPL, genes: list[str]) -> dict[str, list[str]]:
    table = gpl.table.copy()
    probes = {}
    for gene in genes:
        hits = []
        for _, row in table.iterrows():
            annot = str(row.get("SPOT_ID.1", "")) + " " + str(row.get("SPOT_ID", ""))
            gene_pat = re.escape(gene)
            exact = (
                f"({gene})" in annot
                or re.search(rf"\b{gene_pat}\s*\[gene_biotype", annot) is not None
                or re.search(rf"\b{gene_pat}\s*\[Source:HGNC Symbol", annot) is not None
            )
            if exact:
                hits.append(str(row["ID"]))
        probes[gene] = sorted(set(hits))
    return probes


def parse_gsm_metadata(gsm) -> dict[str, str]:
    out = {
        "sample_id": gsm.name,
        "title": gsm.metadata.get("title", [""])[0],
        "source": gsm.metadata.get("source_name_ch1", [""])[0],
    }
    for item in gsm.metadata.get("characteristics_ch1", []):
        if ":" in item:
            k, v = item.split(":", 1)
            out[k.strip().lower().replace(" ", "_").replace("/", "_")] = v.strip()
    return out


def run_gse151179() -> dict:
    gse = GEOparse.get_GEO("GSE151179", destdir=str(GEO_DIR), silent=True)
    gpl = GEOparse.get_GEO("GPL23159", destdir=str(GEO_DIR), silent=True)
    probes = find_exact_gene_probes(gpl, RAI8)

    expr = {}
    meta = []
    for gsm_id, gsm in gse.gsms.items():
        meta.append(parse_gsm_metadata(gsm))
        s = gsm.table.set_index("ID_REF")["VALUE"].astype(float)
        expr[gsm_id] = s
    expr = pd.DataFrame(expr)
    meta = pd.DataFrame(meta).set_index("sample_id")

    gene_expr = {}
    for gene, ids in probes.items():
        present = [p for p in ids if p in expr.index]
        if present:
            gene_expr[gene] = expr.loc[present].mean(axis=0)
        else:
            gene_expr[gene] = pd.Series(index=expr.columns, dtype=float)
    gene_expr = pd.DataFrame(gene_expr)

    z = gene_expr.apply(lambda col: (col - col.mean()) / col.std(ddof=0), axis=0)
    scores = gene_expr.copy()
    scores["rai6_raw_mean"] = gene_expr[RAI6].mean(axis=1)
    scores["rai8_raw_mean"] = gene_expr[RAI8].mean(axis=1)
    scores["rai6_z_mean"] = z[RAI6].mean(axis=1)
    scores["rai8_z_mean"] = z[RAI8].mean(axis=1)
    scores = scores.join(meta, how="left")

    scores["is_tumor"] = scores["source"].str.contains("carcinoma", case=False, na=False)
    scores["is_primary"] = scores.get("tissue_type", "").astype(str).str.contains("Primary tumor", case=False, na=False)
    scores["is_pre_rai"] = scores.get("collection_before_after_rai", "").astype(str).str.lower().eq("before")
    scores["response_refractory"] = scores.get("patient_rai_responce", "").astype(str).str.lower().eq("refractory").astype(int)
    uptake = scores.get("rai_uptake_at_the_metastatic_site", "").astype(str).str.lower()
    scores["uptake_no"] = uptake.eq("no").astype(int)

    rows = []
    subsets = {
        "tumor_all": scores[scores["is_tumor"]],
        "primary_only": scores[scores["is_tumor"] & scores["is_primary"]],
        "pre_rai_tumor": scores[scores["is_tumor"] & scores["is_pre_rai"]],
        "metastatic_tumor": scores[scores["is_tumor"] & ~scores["is_primary"]],
    }
    score_cols = ["rai6_raw_mean", "rai8_raw_mean", "rai6_z_mean", "rai8_z_mean"] + RAI8
    for subset_name, df in subsets.items():
        for endpoint, ycol in [("response_refractory", "response_refractory"), ("uptake_no", "uptake_no")]:
            if ycol not in df:
                continue
            for score_col in score_cols:
                raw_auc, oriented_auc = auc_pair(df[ycol], df[score_col])
                pos = df.loc[df[ycol] == 1, score_col]
                neg = df.loc[df[ycol] == 0, score_col]
                rows.append(
                    {
                        "subset": subset_name,
                        "endpoint": endpoint,
                        "score": score_col,
                        "n": int(df[[ycol, score_col]].dropna().shape[0]),
                        "n_pos": int((df[ycol] == 1).sum()),
                        "n_neg": int((df[ycol] == 0).sum()),
                        "mean_pos": float(pos.mean()) if len(pos) else np.nan,
                        "mean_neg": float(neg.mean()) if len(neg) else np.nan,
                        "cohens_d_pos_vs_neg": cohen_d(pos, neg),
                        "mw_p": mw_p(pos, neg),
                        "auc_pos_high": raw_auc,
                        "auc_best_direction": oriented_auc,
                        "lower_score_predicts_pos": bool(raw_auc < 0.5) if not np.isnan(raw_auc) else None,
                    }
                )
    summary = pd.DataFrame(rows).sort_values(["endpoint", "auc_best_direction"], ascending=[True, False])
    scores.to_csv(OUT / "gse151179_rai_scores.tsv", sep="\t", index=True)
    summary.to_csv(OUT / "gse151179_rai_validation_summary.tsv", sep="\t", index=False)
    probe_rows = [{"gene": g, "n_probes": len(v), "probes": ",".join(v)} for g, v in probes.items()]
    pd.DataFrame(probe_rows).to_csv(OUT / "gse151179_probe_map.tsv", sep="\t", index=False)

    best_response = summary[summary["endpoint"] == "response_refractory"].head(10).to_dict("records")
    best_uptake = summary[summary["endpoint"] == "uptake_no"].head(10).to_dict("records")
    return {
        "n_samples": int(scores.shape[0]),
        "n_tumor": int(scores["is_tumor"].sum()),
        "n_primary": int((scores["is_tumor"] & scores["is_primary"]).sum()),
        "probe_map": {g: len(v) for g, v in probes.items()},
        "best_response_tests": best_response,
        "best_uptake_tests": best_uptake,
    }


def run_local_mechanism_tests() -> dict:
    master = pd.read_csv(ROOT / "project/results/v17/tables/sample_master_v17_tert.tsv", sep="\t")
    master["patient_id"] = master["sample_id"].str[:12]
    master_tumor = master[master["sample_id"].str.contains("-01", na=False)].drop_duplicates("patient_id")

    p2d = pd.read_csv(ROOT / "project/results/dark_matter_phase2/p2d_per_sample_classification.tsv", sep="\t")
    p2d = p2d.merge(master_tumor[["patient_id", "histology_subtype", "stage", "age", "rai_score_v17", "tds16_score_v17"]], left_on="tcga_short", right_on="patient_id", how="left")

    methyl = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
    methyl["patient_id"] = methyl["sample_short"].str[:12]
    p2d = p2d.merge(methyl[["patient_id", "mean_8g_beta"]], on="patient_id", how="left")

    driver_summary = []
    for cls, df in p2d.groupby("driver_class", dropna=False):
        driver_summary.append(
            {
                "driver_class": cls,
                "n": int(len(df)),
                "DM1_n": int((df["v17_dark_cluster"] == "DM1").sum()),
                "DM2_n": int((df["v17_dark_cluster"] == "DM2").sum()),
                "PFI_events": int(df["PFI"].fillna(0).sum()),
                "PFI_rate": float(df["PFI"].mean()),
                "rai_mean": float(df["rai_score_v17"].mean()),
                "tds16_mean": float(df["tds16_score_v17"].mean()),
                "methyl8_mean": float(df["mean_8g_beta"].mean()),
                "age_mean": float(df["age"].mean()),
            }
        )
    pd.DataFrame(driver_summary).to_csv(OUT / "driver_class_score_summary.tsv", sep="\t", index=False)
    p2d.to_csv(OUT / "tcga_driver_class_mechanism_matrix.tsv", sep="\t", index=False)

    # DM1 vs DM2 methylation/expression correlation in BRAF/RAS-negative cases.
    dm = p2d[p2d["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
    corr_rows = []
    for x in ["mean_8g_beta"]:
        for y in ["rai_score_v17", "tds16_score_v17"]:
            sub = dm[[x, y]].dropna()
            if len(sub) > 2:
                r, p = stats.spearmanr(sub[x], sub[y])
                corr_rows.append({"x": x, "y": y, "n": int(len(sub)), "spearman_r": float(r), "p": float(p)})
    pd.DataFrame(corr_rows).to_csv(OUT / "methylation_expression_correlations.tsv", sep="\t", index=False)

    # Arm-level CNV enrichment in class6 true driver-negative vs all others.
    arms = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round8/r8_1_arm_cnv_long.tsv", sep="\t")
    arms["patient_id"] = arms["sample_short"].str[:12]
    arms = arms.merge(p2d[["tcga_short", "driver_class", "v17_dark_cluster"]], left_on="patient_id", right_on="tcga_short", how="inner")
    cnv_rows = []
    for arm, d in arms.groupby("arm"):
        for call in ["Gain", "Loss"]:
            tab = pd.crosstab(d["driver_class"].eq("Class6_True_driver_neg"), d["call"].eq(call))
            for val in [False, True]:
                if val not in tab.index:
                    tab.loc[val] = [0, 0]
            for val in [False, True]:
                if val not in tab.columns:
                    tab[val] = 0
            table = [[tab.loc[True, True], tab.loc[True, False]], [tab.loc[False, True], tab.loc[False, False]]]
            odds, p = stats.fisher_exact(table)
            cnv_rows.append(
                {
                    "arm": arm,
                    "call": call,
                    "class6_call_pct": float(d.loc[d["driver_class"].eq("Class6_True_driver_neg"), "call"].eq(call).mean() * 100),
                    "other_call_pct": float(d.loc[~d["driver_class"].eq("Class6_True_driver_neg"), "call"].eq(call).mean() * 100),
                    "odds_ratio": float(odds) if np.isfinite(odds) else np.nan,
                    "fisher_p": float(p),
                }
            )
    cnv = pd.DataFrame(cnv_rows).sort_values("fisher_p")
    cnv["bh_q"] = np.minimum(1, cnv["fisher_p"] * len(cnv) / (np.arange(len(cnv)) + 1))
    cnv.to_csv(OUT / "class6_arm_cnv_enrichment.tsv", sep="\t", index=False)

    # Pull existing high-signal local results into a concise scorecard.
    methyl_json = json.loads((ROOT / "project/results/audit_2026_04_30/round5/r5_2_methylation_DM.json").read_text())
    f4_json = json.loads((ROOT / "project/results/audit_2026_04_30/round3/f4_fusion_dm_analysis.json").read_text())
    k2_json = json.loads((ROOT / "project/results/dark_matter_phase2/k2_mutation_summary.json").read_text())
    p2d_summary = json.loads((ROOT / "project/results/dark_matter_phase2/p2d_summary.json").read_text())

    return {
        "driver_class_summary": driver_summary,
        "methylation_dm_existing": methyl_json,
        "methyl_expr_corr": corr_rows,
        "top_class6_cnv": cnv.head(10).to_dict("records"),
        "fusion_existing_f4": f4_json,
        "k2_mutation_summary": k2_json,
        "p2d_summary": p2d_summary,
    }


def write_markdown(gse: dict, mech: dict) -> None:
    best_response = pd.DataFrame(gse["best_response_tests"])
    best_uptake = pd.DataFrame(gse["best_uptake_tests"])
    driver_summary = pd.DataFrame(mech["driver_class_summary"])
    top_cnv = pd.DataFrame(mech["top_class6_cnv"])
    corr = pd.DataFrame(mech["methyl_expr_corr"])

    def md_table(df: pd.DataFrame) -> str:
        if df is None or len(df) == 0:
            return "No rows."
        d = df.copy()
        for col in d.columns:
            if pd.api.types.is_float_dtype(d[col]):
                d[col] = d[col].map(lambda x: "" if pd.isna(x) else f"{x:.4g}")
            else:
                d[col] = d[col].map(lambda x: "" if pd.isna(x) else str(x))
        header = "| " + " | ".join(d.columns) + " |"
        sep = "| " + " | ".join(["---"] * len(d.columns)) + " |"
        rows = ["| " + " | ".join(map(str, row)) + " |" for row in d.to_numpy()]
        return "\n".join([header, sep] + rows)

    # Decision ranking intentionally penalizes overlap with Paper 1 and endpoint weakness.
    lines = []
    lines.append("# Aggressive-topic feasibility sprint\n")
    lines.append("Date: 2026-05-06\n")
    lines.append("Scope: rapid local + public-download tests; no Paper 1 edits.\n")
    lines.append("## Bottom line\n")
    lines.append("The strongest aggressive new angle is **not** generic dark-matter taxonomy. The best candidate is **DM1/true driver-negative hidden fusion or regulatory mechanism**, but it still needs independent validation because local fusion evidence conflicts by source. The direct RAI-refractory paper is currently weak after re-scoring GSE151179.\n")

    lines.append("## 1. Direct RAI-refractory validation: weak / no-go as headline\n")
    lines.append(f"GSE151179 parsed from GEO: {gse['n_samples']} samples, {gse['n_tumor']} tumors, {gse['n_primary']} primary tumors. Probe recovery: {gse['probe_map']}.\n")
    lines.append("Best refractory-response tests by oriented AUC:\n")
    lines.append(md_table(best_response[["subset", "score", "n", "n_pos", "n_neg", "mean_pos", "mean_neg", "cohens_d_pos_vs_neg", "mw_p", "auc_pos_high", "auc_best_direction", "lower_score_predicts_pos"]].head(8)))
    lines.append("\n\nBest metastatic RAI-uptake tests by oriented AUC:\n")
    lines.append(md_table(best_uptake[["subset", "score", "n", "n_pos", "n_neg", "mean_pos", "mean_neg", "cohens_d_pos_vs_neg", "mw_p", "auc_pos_high", "auc_best_direction", "lower_score_predicts_pos"]].head(8)))
    lines.append("\n\nInterpretation: if lower thyroid-lineage/RAI score predicts refractory/no-uptake, the signal is at best modest and sample-imbalanced. This can support a limitation or pilot validation, not an aggressive standalone RAI-refractory classifier.\n")

    lines.append("## 2. True driver-negative / hidden mechanism: most interesting, but validation-gated\n")
    lines.append("Driver-class score summary:\n")
    lines.append(md_table(driver_summary))
    lines.append("\n\nMethylation-expression correlations in DM1/DM2:\n")
    lines.append(md_table(corr) if len(corr) else "No correlations computed.")
    lines.append("\n\nExisting methylation result: DM1 is hypermethylated across the 8 thyroid-lineage/RAI genes versus DM2; TPO shows delta beta 0.415, Cohen d 2.299, p=1.88e-18. This is a real mechanism signal, but it overlaps heavily with the Paper 1 differentiation axis.\n")
    lines.append("Class6 arm-level CNV enrichment top rows:\n")
    lines.append(md_table(top_cnv[["arm", "call", "class6_call_pct", "other_call_pct", "odds_ratio", "fisher_p", "bh_q"]]))
    lines.append("\n\nInterpretation: class6 CNV arm hits are not strong after multiple testing. CNV alone is not an aggressive paper.\n")

    lines.append("## 3. Fusion angle: high upside, current evidence inconsistent\n")
    f4 = mech["fusion_existing_f4"]
    lines.append(f"PanCancer/SV-style local fusion result reports DM1 fusion-positive {f4['fusion_rate_by_dm']['fusion_pct'].get('DM1')}%, DM2 {f4['fusion_rate_by_dm']['fusion_pct'].get('DM2')}%, not-DM {f4['fusion_rate_by_dm']['fusion_pct'].get('not_DM')}%. RET/NTRK/ALK/BRAF fusion classes are concentrated in DM1 in that table.\n")
    lines.append("But v17 curated fusion overlay is sparse: DM1 has ALK 1, NTRK 1, RET 1; DM2 has 0. Therefore this is **the most aggressive testable new topic**, but it is not publishable until raw RNA/SV calls are reconciled and externally validated.\n")

    lines.append("## 4. DICER1/EIF1AX class: real but small\n")
    k2 = mech["k2_mutation_summary"]
    p2d = mech["p2d_summary"]
    lines.append(f"TCGA class counts: {p2d['class_counts']}. K2/Yoo: BRAF/RAS-negative {k2['dark_matter']['n']}/{k2['n_total']} = {k2['dark_matter']['pct']}%; DICER1+EIF1AX in K2 dark matter = {k2['dark_matter_alt_driver']['DICER1'] + k2['dark_matter_alt_driver']['EIF1AX']} cases.\n")
    lines.append("Interpretation: biologically clean, replicated frequency, but n=7 in TCGA and n=7 in K2 dark-matter. Good as a class inside a broader paper, weak as a standalone aggressive paper.\n")

    lines.append("## 5. TERT-only class: too small right now\n")
    lines.append("Local class map: TERT-only n=5, PFI events 2/5. This is provocative but too small for a standalone paper without external TERT-rich cohort access.\n")

    lines.append("## Ranked recommendation\n")
    rows = [
        ["1", "DM1 hidden fusion/regulatory driver", "Highest novelty; least Paper 1 overlap", "Needs raw fusion/SV reconciliation + external validation", "Attack next"],
        ["2", "True driver-negative epigenetic silencing", "Strong methylation signal", "Overlaps Paper 1 lineage axis", "Only if framed as mechanism, not biomarker"],
        ["3", "DICER1/EIF1AX FVPTC-like class", "Replicated rare-driver enrichment", "Small n", "Use as section, not standalone"],
        ["4", "RAI-refractory classifier", "Direct public labels available", "GSE151179 signal weak/imbalanced", "No-go as headline today"],
        ["5", "TERT-only triple-negative", "Provocative event rate", "n=5", "Hold until external cohort"],
    ]
    lines.append(md_table(pd.DataFrame(rows, columns=["rank", "topic", "why", "blocker", "decision"])))
    lines.append("\n## Immediate next work if continuing today\n")
    lines.append("1. Reconcile fusion evidence: compare cBioPortal PanCancer structural variants, v17 fusion overlay, and raw GSE184362/GSE232237 metadata.  \n2. If raw FASTQ or STAR-Fusion-compatible files exist for K2/GSE184362, run fusion calling on BRAF/RAS-negative or DM1-like samples.  \n3. For epigenetic mechanism, test whether promoter methylation explains expression after adjusting for histology and BRAF/RAS class.  \n")

    (OUT / "aggressive_topic_sprint_report.md").write_text("\n".join(lines))


def main() -> None:
    gse = run_gse151179()
    mech = run_local_mechanism_tests()
    combined = {"gse151179": gse, "mechanism": mech}
    (OUT / "aggressive_sprint_summary.json").write_text(json.dumps(combined, indent=2, ensure_ascii=False))
    write_markdown(gse, mech)
    print(f"Wrote sprint outputs to {OUT}")


if __name__ == "__main__":
    main()
