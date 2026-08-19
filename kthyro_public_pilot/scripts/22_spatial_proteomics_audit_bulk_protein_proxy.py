#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr

from pilot_utils import fdr_bh, pilot_root, setup_logging, zscore_df


PROTEIN_MODULES = {
    "RAI_differentiation_protein": ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "FOXE1", "DIO1"],
    "HLA_I_APM_protein": ["B2M", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    "HLA_II_APC_protein": ["CD74"],
    "CD8_cytotoxic_protein": ["CD8A", "PRF1"],
    "IFN_APM_activation_protein": ["STAT1", "IRF1", "PSMB8", "PSMB9"],
    "checkpoint_suppression_protein": ["TOX", "LAG3", "CD274", "PDCD1LG2"],
    "myeloid_suppressive_protein": ["ITGAM", "MRC1", "CD163", "TGFB1"],
    "TLS_B_cell_protein": ["CXCL13", "CCL21", "JCHAIN", "MS4A1", "CD79A", "MZB1"],
    "TGFB_barrier_single_protein": ["TGFB1"],
}

GROUP_ORDER = {"PTC": 0, "PTC_P": 0.5, "PDTC": 1, "ATC_P": 1.5, "ATC": 2, "T_Ca_m": 2}
CONTRASTS = [("ATC", "PTC"), ("PDTC", "PTC"), ("ATC", "PDTC"), ("ATC_P", "PTC_P")]


def cohens_d(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / pooled) if np.isfinite(pooled) and pooled > 0 else np.nan


def load_protein(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t")
    expr = raw.set_index("gene").T
    expr.index.name = "sample_id"
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr["group"] = [str(x).split("|")[-1] for x in expr.index]
    return expr.reset_index()


def score_modules(expr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = [c for c in expr.columns if c not in {"sample_id", "group"}]
    z = zscore_df(expr[genes])
    out = expr[["sample_id", "group"]].copy()
    coverage = []
    for name, module_genes in PROTEIN_MODULES.items():
        found = [g for g in module_genes if g in z.columns]
        out[name] = z[found].mean(axis=1) if found else np.nan
        coverage.append(
            {
                "protein_module": name,
                "n_requested": len(module_genes),
                "n_found": len(found),
                "found_proteins": ",".join(found),
                "claim_boundary": "Bulk proteomics availability; not spatial proteomics.",
            }
        )
    out["protein_aggressive_visibility_score"] = -out["RAI_differentiation_protein"] + out["HLA_I_APM_protein"] + out["myeloid_suppressive_protein"]
    return out, pd.DataFrame(coverage)


def contrasts(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    modules = [c for c in scores.columns if c not in {"sample_id", "group"}]
    for a_group, b_group in CONTRASTS:
        a_df = scores[scores["group"].eq(a_group)]
        b_df = scores[scores["group"].eq(b_group)]
        if len(a_df) < 2 or len(b_df) < 2:
            continue
        for mod in modules:
            a = pd.to_numeric(a_df[mod], errors="coerce")
            b = pd.to_numeric(b_df[mod], errors="coerce")
            p = mannwhitneyu(a.dropna(), b.dropna(), alternative="two-sided").pvalue if len(a.dropna()) >= 2 and len(b.dropna()) >= 2 else np.nan
            rows.append(
                {
                    "contrast": f"{a_group}_vs_{b_group}",
                    "group_a": a_group,
                    "group_b": b_group,
                    "module": mod,
                    "n_group_a": len(a.dropna()),
                    "n_group_b": len(b.dropna()),
                    "mean_group_a": a.mean(),
                    "mean_group_b": b.mean(),
                    "cohens_d_group_a_minus_b": cohens_d(a, b),
                    "mannwhitney_p": p,
                    "claim_boundary": "Bulk proteomics contrast validates protein-level directionality, not spatial localization.",
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_q"] = fdr_bh(out["mannwhitney_p"])
    return out


def trends(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    scores = scores.copy()
    scores["dediff_order"] = scores["group"].map(GROUP_ORDER)
    for mod in [c for c in scores.columns if c not in {"sample_id", "group", "dediff_order"}]:
        frame = scores[["dediff_order", mod]].replace([np.inf, -np.inf], np.nan).dropna()
        rho, p = (np.nan, np.nan)
        if len(frame) >= 5 and frame["dediff_order"].nunique() > 1 and frame[mod].nunique() > 1:
            rho, p = spearmanr(frame["dediff_order"], frame[mod])
        rows.append(
            {
                "module": mod,
                "n_samples": len(frame),
                "spearman_dediff_order_rho": rho,
                "spearman_p": p,
                "mean_PTC": scores.loc[scores["group"].eq("PTC"), mod].mean(),
                "mean_PDTC": scores.loc[scores["group"].eq("PDTC"), mod].mean(),
                "mean_ATC": scores.loc[scores["group"].eq("ATC"), mod].mean(),
                "claim_boundary": "Bulk proteomics dedifferentiation trend; not spatial and not treatment-response proof.",
            }
        )
    out = pd.DataFrame(rows)
    out["fdr_q"] = fdr_bh(out["spearman_p"])
    return out.sort_values("fdr_q")


def audit_table(root: Path) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "resource": "GSE301163",
                "modality": "GeoMx DSP spatial transcriptomics",
                "thyroid_specific": "yes",
                "spatial": "yes_ROI",
                "protein_level": "no_detected_as_RNA_matrix",
                "public_matrix_available": "yes",
                "local_path": str(root.parent / "project/data/external_geomx_GSE301163/GSE301163_Normalized_ST_matrix.txt.gz"),
                "usable_for_this_question": "spatial_RNA_only",
                "claim_boundary": "Good spatial ROI transcriptomics; not spatial proteomics.",
            },
            {
                "resource": "Donati et al. Virchows Archiv 2025",
                "modality": "GeoMx DSP pathologist demonstration",
                "thyroid_specific": "includes thyroid anaplastic carcinoma example",
                "spatial": "yes_ROI",
                "protein_level": "platform can measure protein/RNA but public count matrix not located",
                "public_matrix_available": "not_found",
                "local_path": "",
                "usable_for_this_question": "literature_feasibility_only",
                "claim_boundary": "Supports feasibility of DSP in thyroid ATC morphology; no reusable public protein matrix found.",
            },
            {
                "resource": "Haq/Bychkov/Mete/Jeon/Jung Endocrine Pathology 2025",
                "modality": "Spatial transcriptomics plus IHC validation",
                "thyroid_specific": "yes_ATC",
                "spatial": "yes_ROI_for_RNA",
                "protein_level": "IHC validation of selected markers",
                "public_matrix_available": "not_found_for_spatial_protein",
                "local_path": "",
                "usable_for_this_question": "literature_marker_support",
                "claim_boundary": "Supports marker prioritization, not a raw public spatial proteomics dataset.",
            },
            {
                "resource": "Local proteogenomic_v1 Mun/CellRepMed-style protein table",
                "modality": "Bulk proteomics",
                "thyroid_specific": "yes",
                "spatial": "no",
                "protein_level": "yes",
                "public_matrix_available": "local_processed",
                "local_path": str(root.parent / "project/results/proteogenomic_v1/paper3_mun2025_dediff_layer/protein_target_genes_S1C.tsv"),
                "usable_for_this_question": "bulk_protein_proxy",
                "claim_boundary": "Useful protein-level orthogonal validation but cannot prove spatial protein niche localization.",
            },
        ]
    )


def marker_support_table() -> pd.DataFrame:
    rows = []
    for marker, direction, axis, source in [
        ("COL7A1", "ATC protein/IHC up in literature", "ECM/invasion barrier", "Endocrine Pathology 2025 ATC spatial transcriptomic + IHC study"),
        ("LAMC2", "ATC protein/IHC up in literature", "ECM/invasion barrier", "Endocrine Pathology 2025 ATC spatial transcriptomic + IHC study"),
        ("SPHK1", "ATC protein/IHC up in literature", "aggressive/metabolic stress", "Endocrine Pathology 2025 ATC spatial transcriptomic + IHC study"),
        ("SRPX2", "ATC protein/IHC up in literature", "EMT/invasion", "Endocrine Pathology 2025 ATC spatial transcriptomic + IHC study"),
        ("CD24", "ATC protein/IHC down in literature", "differentiation/luminal marker loss", "Endocrine Pathology 2025 ATC spatial transcriptomic + IHC study"),
        ("EPHX1", "ATC protein/IHC down in literature", "differentiation/metabolic marker loss", "Endocrine Pathology 2025 ATC spatial transcriptomic + IHC study"),
        ("GPX3", "ATC protein/IHC down in literature", "oxidative/differentiation marker loss", "Endocrine Pathology 2025 ATC spatial transcriptomic + IHC study"),
        ("RBM47", "ATC protein/IHC down in literature", "epithelial differentiation loss", "Endocrine Pathology 2025 ATC spatial transcriptomic + IHC study"),
    ]:
        rows.append(
            {
                "marker": marker,
                "reported_direction": direction,
                "kthyro_axis_relevance": axis,
                "source": source,
                "proposal_use": "extended_validation_panel_candidate",
                "claim_boundary": "Literature/IHC marker support; not reanalyzed raw public spatial proteomics.",
            }
        )
    return pd.DataFrame(rows)


def plot_scores(scores: pd.DataFrame, fig_dir: Path) -> None:
    modules = [
        "RAI_differentiation_protein",
        "HLA_I_APM_protein",
        "CD8_cytotoxic_protein",
        "myeloid_suppressive_protein",
        "checkpoint_suppression_protein",
        "protein_aggressive_visibility_score",
    ]
    long = scores.melt(id_vars=["sample_id", "group"], value_vars=[m for m in modules if m in scores.columns], var_name="module", value_name="score")
    order = ["PTC", "PDTC", "ATC", "PTC_P", "ATC_P", "T_Ca_m"]
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.8), squeeze=False)
    for ax, mod in zip(axes.ravel(), modules):
        sub = long[long["module"].eq(mod)]
        sns.boxplot(data=sub, x="group", y="score", order=[g for g in order if g in set(sub["group"])], color="#9ecae1", fliersize=1.5, ax=ax)
        ax.set_title(mod.replace("_", " "), fontsize=10, weight="bold")
        ax.tick_params(axis="x", rotation=35)
        ax.set_xlabel("")
        ax.set_ylabel("protein module mean-z")
    fig.suptitle("Bulk proteomics validates K-Thyro protein-level axes (not spatial)", fontsize=15, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(fig_dir / "bulk_proteomics_kthyro_axis_boxplots.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "bulk_proteomics_kthyro_axis_boxplots.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_contrasts(contrast_df: pd.DataFrame, fig_dir: Path) -> None:
    sub = contrast_df[contrast_df["contrast"].isin(["ATC_vs_PTC", "PDTC_vs_PTC"])].copy()
    if sub.empty:
        return
    pivot = sub.pivot(index="module", columns="contrast", values="cohens_d_group_a_minus_b")
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    sns.heatmap(pivot, cmap="vlag", center=0, annot=True, fmt=".2f", linewidths=0.35, linecolor="#222", ax=ax)
    ax.set_title("Protein-level ATC/PDTC vs PTC effects", weight="bold")
    fig.tight_layout()
    fig.savefig(fig_dir / "bulk_proteomics_kthyro_contrast_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "bulk_proteomics_kthyro_contrast_heatmap.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=pilot_root())
    parser.add_argument("--protein-path", type=Path, default=pilot_root().parent / "project/results/proteogenomic_v1/paper3_mun2025_dediff_layer/protein_target_genes_S1C.tsv")
    args = parser.parse_args()
    logger = setup_logging("22_spatial_proteomics_audit_bulk_protein_proxy")
    out = args.root / "results" / "extra_analyses"
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    audit = audit_table(args.root)
    audit.to_csv(tables / "spatial_proteomics_public_resource_audit.tsv", sep="\t", index=False)
    marker_support_table().to_csv(tables / "literature_iHC_spatial_marker_support_candidates.tsv", sep="\t", index=False)

    if args.protein_path.exists():
        protein = load_protein(args.protein_path)
        scores, coverage = score_modules(protein)
        contrast_df = contrasts(scores)
        trend_df = trends(scores)
        scores.to_csv(tables / "bulk_proteomics_kthyro_module_scores.tsv", sep="\t", index=False)
        coverage.to_csv(tables / "bulk_proteomics_kthyro_module_coverage.tsv", sep="\t", index=False)
        contrast_df.to_csv(tables / "bulk_proteomics_kthyro_module_contrasts.tsv", sep="\t", index=False)
        trend_df.to_csv(tables / "bulk_proteomics_kthyro_dediff_trends.tsv", sep="\t", index=False)
        plot_scores(scores, figs)
        plot_contrasts(contrast_df, figs)
        logger.info("Wrote protein scores for %d samples.", len(scores))
    else:
        pd.DataFrame([{"status": "not_run", "reason": f"missing {args.protein_path}"}]).to_csv(tables / "bulk_proteomics_kthyro_module_scores.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
