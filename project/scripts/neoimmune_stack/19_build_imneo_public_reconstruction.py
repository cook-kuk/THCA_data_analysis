#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import (
    ensure_run_dir,
    metric_binary,
    patient_topn_metrics,
    precision_recall_at_k,
    safe_read_table,
    write_md,
    write_tsv,
)


LABEL_COL = "label_immunogenicity"

PATENT_DISCLOSED_SOURCES = [
    {
        "source_name": "NEPdb",
        "patent_role": "training-like public source",
        "canonical_match_rule": "source_study contains NEPdb",
        "status": "present",
        "claim_boundary": "public-source reconstruction only; not CG Invites proprietary peptide list",
    },
    {
        "source_name": "dbPepNeo",
        "patent_role": "training-like public source",
        "canonical_match_rule": "source_study contains dbPepNeo",
        "status": "present as dbPepNeo2",
        "claim_boundary": "public-source reconstruction only; counts need not match patent training split",
    },
    {
        "source_name": "PRIME",
        "patent_role": "training-like public source / public predictor comparator",
        "canonical_match_rule": "PRIME score exists as frozen comparator; no PRIME source table",
        "status": "predictor score present, source table absent",
        "claim_boundary": "can compare PRIME score, cannot reconstruct PRIME training rows",
    },
    {
        "source_name": "INeo-Epp",
        "patent_role": "training-like public source",
        "canonical_match_rule": "source_study contains INeo or Epp",
        "status": "absent",
        "claim_boundary": "missing from current local canonical table",
    },
    {
        "source_name": "TESLA",
        "patent_role": "independent-like public benchmark",
        "canonical_match_rule": "source_study contains TESLA",
        "status": "present",
        "claim_boundary": "patient-matched public benchmark; not the proprietary imNEO selected set",
    },
    {
        "source_name": "IEDB t-cell DB",
        "patent_role": "independent-like public benchmark",
        "canonical_match_rule": "source_study contains IEDB",
        "status": "absent as explicit source",
        "claim_boundary": "do not silently relabel CEDAR as IEDB",
    },
    {
        "source_name": "McPAS-TCR",
        "patent_role": "independent-like public benchmark",
        "canonical_match_rule": "source_study contains McPAS",
        "status": "absent",
        "claim_boundary": "requires additional ingestion",
    },
    {
        "source_name": "VDJdb",
        "patent_role": "independent-like public benchmark",
        "canonical_match_rule": "source_study contains VDJ",
        "status": "absent",
        "claim_boundary": "requires additional ingestion",
    },
]

NMI_MODELS = {
    "NMI_all_branch_mean": "NMI all-branch mean",
    "NMI_tcr_qk_structure": "NMI TCR-QK-Structure",
    "NMI_equal_top5": "NMI equal top5",
    "NMI_w7b_esm2_qk": "NMI W7B-ESM2-QK",
}

LOCAL_ABLATIONS = {
    "Wave8_TCR_SelfSim_full": "Wave8 TCR/SelfSim",
    "Wave8_TCR_SelfSim_no_exact": "Wave8 TCR/SelfSim no-exact",
    "Structure_LR": "Structure LR",
    "ESM2_Bayesian": "ESM2 Bayesian",
    "W7A_QK_only": "W7A QK only",
    "W7B_stacked": "W7B stacked",
    "GP_quantum": "GP quantum",
}

PUBLIC_COMPARATORS = {
    "BigMHC_IM": "BigMHC-IM",
    "MHCflurry_2.0_presentation": "MHCflurry presentation",
    "PRIME": "PRIME",
    "NetMHCpan_4.1_EL": "NetMHCpan EL",
}


def classify_source(source: object) -> str:
    s = "" if pd.isna(source) else str(source).lower()
    if "nepdb" in s:
        return "NEPdb"
    if "dbpepneo" in s:
        return "dbPepNeo"
    if "tesla" in s:
        return "TESLA"
    if "cedar" in s or "improve" in s:
        return "CEDAR/IMPROVE context"
    if "itsndb" in s:
        return "ITSNdb context"
    if "iedb" in s:
        return "IEDB"
    if "mcpas" in s:
        return "McPAS-TCR"
    if "vdj" in s:
        return "VDJdb"
    if "ineo" in s or "epp" in s:
        return "INeo-Epp"
    return "other"


def pivot_scores(scores: pd.DataFrame, name_col: str, value_col: str) -> pd.DataFrame:
    if scores.empty or "candidate_id" not in scores:
        return pd.DataFrame(columns=["candidate_id"])
    d = scores.copy()
    d[value_col] = pd.to_numeric(d[value_col], errors="coerce")
    return d.pivot_table(index="candidate_id", columns=name_col, values=value_col, aggfunc="mean").reset_index()


def add_metric(
    rows: list[dict],
    df: pd.DataFrame,
    score_col: str,
    model_label: str,
    subset_name: str,
    track: str,
    clean_allowed: bool,
    source_family: str,
    claim_boundary: str,
) -> None:
    if score_col not in df:
        return
    use = df[df[LABEL_COL].notna()].copy()
    use[score_col] = pd.to_numeric(use[score_col], errors="coerce")
    use = use[use[score_col].notna()]
    if use.empty:
        return
    m = metric_binary(use[LABEL_COL], use[score_col])
    for k in [10, 20, 34, 50]:
        m.update(precision_recall_at_k(use, score_col, LABEL_COL, k))
    for n in [20, 34]:
        m.update(patient_topn_metrics(use, score_col, LABEL_COL, n))
    patient_ids = use["patient_id"].fillna("unknown_patient").astype(str) if "patient_id" in use else pd.Series([], dtype=str)
    real_patients = patient_ids[~patient_ids.isin(["", "unknown", "unknown_patient", "nan", "None"])].nunique()
    patient_level_claimable = bool(real_patients >= 2)
    if not patient_level_claimable:
        for n in [20, 34]:
            m[f"patient_hit_rate@{n}"] = np.nan
            m[f"patient_recall@{n}"] = np.nan
        m["patients_evaluated"] = 0.0
    m.update(
        {
            "subset": subset_name,
            "source_family": source_family,
            "track": track,
            "model": model_label,
            "score_col": score_col,
            "clean_track_allowed": clean_allowed,
            "claim_boundary": claim_boundary,
            "unique_patients": float(patient_ids.nunique()) if len(patient_ids) else np.nan,
            "real_patients": float(real_patients),
            "patient_level_claimable": patient_level_claimable,
            "unique_sources": ", ".join(sorted(use["source_family"].dropna().astype(str).unique())) if "source_family" in use else "",
        }
    )
    rows.append(m)


def write_figures(outdir: Path, leaderboard: pd.DataFrame, source_inventory: pd.DataFrame) -> None:
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    try:
        import cv2
        import matplotlib.pyplot as plt

        head = leaderboard[
            leaderboard["subset"].isin(["imNEO_core_head_to_head_common", "imNEO_core_nmi_locked_rows"])
        ].copy()
        if head.empty:
            head = leaderboard.copy()
        head = head.sort_values(["subset", "AUPRC"], ascending=[True, False]).head(16)

        fig, ax = plt.subplots(figsize=(12, 7), facecolor="#080b12")
        ax.set_facecolor("#111827")
        colors = [
            "#26d9ff" if str(x).startswith("NMI") else "#5ee6a8" if "Wave8" in str(x) else "#a78bfa"
            for x in head["model"]
        ]
        labels = head["model"] + "\n" + head["subset"].str.replace("imNEO_", "", regex=False)
        ax.barh(labels, head["AUPRC"], color=colors)
        ax.invert_yaxis()
        ax.set_xlabel("AUPRC")
        ax.set_title("imNEO-public reconstruction: confident clean algorithms vs public comparators", color="white", fontsize=14)
        ax.tick_params(colors="white", labelsize=9)
        ax.xaxis.label.set_color("white")
        for spine in ax.spines.values():
            spine.set_color("#263348")
        ax.grid(axis="x", color="#263348", alpha=0.55)
        fig.tight_layout()
        raw = figdir / "imneo_public_reconstruction_leaderboard_raw.png"
        fig.savefig(raw, dpi=220)
        plt.close(fig)

        img = cv2.imread(str(raw))
        if img is not None:
            band_h = 92
            band = np.zeros((band_h, img.shape[1], 3), dtype=np.uint8)
            band[:, :] = (18, 11, 8)
            title = "NMI is compared only on public-source rows visible from patent/presentation context"
            cv2.putText(band, title, (32, 38), cv2.FONT_HERSHEY_SIMPLEX, 0.78, (255, 217, 38), 2, cv2.LINE_AA)
            cv2.putText(
                band,
                "Claim boundary: not proprietary imNEO peptides; retrospective labels; same-row head-to-head where possible.",
                (32, 70),
                cv2.FONT_HERSHEY_SIMPLEX,
                0.52,
                (199, 176, 159),
                1,
                cv2.LINE_AA,
            )
            img2 = np.vstack([band, img])
            cv2.imwrite(str(figdir / "imneo_public_reconstruction_leaderboard.png"), img2)

        inv = source_inventory.copy()
        inv["present_rows"] = pd.to_numeric(inv["present_rows"], errors="coerce").fillna(0)
        fig, ax = plt.subplots(figsize=(11, 5), facecolor="#080b12")
        ax.set_facecolor("#111827")
        ax.bar(inv["source_name"], inv["present_rows"], color="#f6c768")
        ax.set_yscale("symlog")
        ax.set_ylabel("canonical rows (symlog)")
        ax.set_title("Patent-disclosed public sources: local reconstruction coverage", color="white", fontsize=14)
        ax.tick_params(colors="white", labelrotation=25)
        ax.yaxis.label.set_color("white")
        for spine in ax.spines.values():
            spine.set_color("#263348")
        ax.grid(axis="y", color="#263348", alpha=0.55)
        fig.tight_layout()
        fig.savefig(figdir / "imneo_public_source_coverage.png", dpi=220)
        plt.close(fig)
    except Exception:
        pass


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)

    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    canon[LABEL_COL] = pd.to_numeric(canon[LABEL_COL], errors="coerce")
    canon["source_family"] = canon["source_study"].fillna(canon["dataset_source"]).map(classify_source)

    nmi_path = outdir / "predictions" / "nmi_clean_method_predictions.tsv"
    if not nmi_path.exists():
        raise FileNotFoundError(f"Run 18_train_nmi_clean_method.py first: {nmi_path}")
    nmi = safe_read_table(nmi_path)
    nmi_cols = [c for c in NMI_MODELS if c in nmi.columns]

    local = safe_read_table(outdir / "local_scores" / "all_local_scores.tsv.gz")
    local_wide = pivot_scores(local[local["local_model_name"].isin(LOCAL_ABLATIONS)], "local_model_name", "normalized_score")

    ext = safe_read_table(outdir / "external_scores" / "all_external_scores.tsv.gz")
    ext_wide = pivot_scores(ext[ext["model_name"].isin(PUBLIC_COMPARATORS)], "model_name", "normalized_score")

    base_cols = [
        "candidate_id",
        "patient_id",
        "sample_id",
        "dataset_source",
        "source_study",
        "source_family",
        "tumor_type",
        "gene",
        "mutation_id",
        "peptide_mut",
        "peptide_wt",
        "peptide_length",
        "hla_allele",
        "hla_class",
        "assay_type",
        "validation_level",
        LABEL_COL,
        "any_existing_overlap_flag",
        "strict_set_flag",
    ]
    for c in base_cols:
        if c not in canon:
            canon[c] = np.nan
    d = canon[base_cols].merge(nmi[["candidate_id"] + nmi_cols], on="candidate_id", how="left")
    d = d.merge(local_wide, on="candidate_id", how="left").merge(ext_wide, on="candidate_id", how="left")

    core = d[d["source_family"].isin(["NEPdb", "dbPepNeo", "TESLA"])].copy()
    locked_rows = core[core["NMI_all_branch_mean"].notna()].copy()
    common_cols = ["NMI_all_branch_mean", "BigMHC_IM", "MHCflurry_2.0_presentation", "PRIME"]
    common = core.dropna(subset=[c for c in common_cols if c in core.columns]).copy()
    strict_common = common[~common["any_existing_overlap_flag"].fillna(False).astype(bool)].copy()
    nepdb = core[core["source_family"].eq("NEPdb")].copy()
    tesla = core[core["source_family"].eq("TESLA")].copy()

    subsets = [
        ("imNEO_core_all_available", core, "NEPdb + dbPepNeo + TESLA"),
        ("imNEO_core_nmi_locked_rows", locked_rows, "NMI-locked rows inside patent-like core"),
        ("imNEO_core_head_to_head_common", common, "NMI + BigMHC + MHCflurry + PRIME common rows"),
        ("imNEO_core_strict_no_existing_overlap_common", strict_common, "common rows without existing-overlap flag"),
        ("NEPdb_all_available", nepdb, "NEPdb only"),
        ("TESLA_all_available", tesla, "TESLA only"),
    ]

    rows: list[dict] = []
    for subset_name, sub, source_desc in subsets:
        if sub.empty:
            continue
        for col, label in NMI_MODELS.items():
            add_metric(
                rows,
                sub,
                col,
                label,
                subset_name,
                "confident_clean_NMI",
                True,
                source_desc,
                "own clean algorithm; public-source reconstruction, not proprietary imNEO peptide list",
            )
        for col, label in LOCAL_ABLATIONS.items():
            add_metric(
                rows,
                sub,
                col,
                label,
                subset_name,
                "local_clean_ablation",
                True,
                source_desc,
                "local branch ablation; not a standalone clinical predictor",
            )
        for col, label in PUBLIC_COMPARATORS.items():
            add_metric(
                rows,
                sub,
                col,
                label,
                subset_name,
                "frozen_public_comparator",
                False,
                source_desc,
                "frozen public comparator; not used as NMI training feature",
            )

    leaderboard = pd.DataFrame(rows)
    if not leaderboard.empty:
        order = [
            "subset",
            "source_family",
            "track",
            "model",
            "clean_track_allowed",
            "n",
            "positives",
            "unique_patients",
            "real_patients",
            "patient_level_claimable",
            "AUROC",
            "AUPRC",
            "Precision@10",
            "Precision@20",
            "Precision@34",
            "Recall@20",
            "Recall@34",
            "patient_hit_rate@20",
            "patient_hit_rate@34",
            "patients_evaluated",
            "claim_boundary",
        ]
        for c in order:
            if c not in leaderboard:
                leaderboard[c] = np.nan
        leaderboard = leaderboard[order + [c for c in leaderboard.columns if c not in order]]
        leaderboard = leaderboard.sort_values(["subset", "AUPRC", "AUROC"], ascending=[True, False, False])
    write_tsv(leaderboard, outdir / "metrics" / "imneo_public_reconstruction_leaderboard.tsv")

    inv_rows = []
    for item in PATENT_DISCLOSED_SOURCES:
        name = item["source_name"]
        fam = classify_source(name)
        if name == "dbPepNeo":
            mask = d["source_family"].eq("dbPepNeo")
        elif name == "IEDB t-cell DB":
            mask = d["source_family"].eq("IEDB")
        elif name == "PRIME":
            mask = pd.Series(False, index=d.index)
        else:
            mask = d["source_family"].eq(fam)
        present = d[mask]
        inv_rows.append(
            {
                **item,
                "present_rows": int(len(present)),
                "positive_labels": int(pd.to_numeric(present[LABEL_COL], errors="coerce").fillna(0).sum()) if len(present) else 0,
                "nmi_locked_rows": int(present["NMI_all_branch_mean"].notna().sum()) if len(present) and "NMI_all_branch_mean" in present else 0,
                "public_comparator_rows": int(
                    present[[c for c in PUBLIC_COMPARATORS if c in present]].notna().any(axis=1).sum()
                )
                if len(present)
                else 0,
            }
        )
    source_inventory = pd.DataFrame(inv_rows)
    write_tsv(source_inventory, outdir / "reports" / "imneo_public_source_inventory.tsv")

    ranked_cols = [
        "candidate_id",
        "patient_id",
        "dataset_source",
        "source_family",
        "tumor_type",
        "gene",
        "peptide_mut",
        "peptide_wt",
        "hla_allele",
        LABEL_COL,
        "NMI_all_branch_mean",
        "NMI_tcr_qk_structure",
        "Wave8_TCR_SelfSim_full",
        "Structure_LR",
        "ESM2_Bayesian",
        "W7A_QK_only",
        "BigMHC_IM",
        "MHCflurry_2.0_presentation",
        "PRIME",
        "NetMHCpan_4.1_EL",
        "any_existing_overlap_flag",
    ]
    for c in ranked_cols:
        if c not in d:
            d[c] = np.nan
    common_ranked = common.copy()
    common_ranked["NMI_primary_score"] = common_ranked["NMI_all_branch_mean"].combine_first(common_ranked["NMI_tcr_qk_structure"])
    write_tsv(
        common_ranked[ranked_cols + ["NMI_primary_score"]].sort_values("NMI_primary_score", ascending=False),
        outdir / "predictions" / "imneo_public_reconstruction_ranked_candidates.tsv",
    )

    write_figures(outdir, leaderboard, source_inventory)

    best = leaderboard[
        leaderboard["subset"].eq("imNEO_core_head_to_head_common") & leaderboard["model"].astype(str).str.startswith("NMI")
    ].sort_values("AUPRC", ascending=False)
    best_public = leaderboard[
        leaderboard["subset"].eq("imNEO_core_head_to_head_common") & leaderboard["track"].eq("frozen_public_comparator")
    ].sort_values("AUPRC", ascending=False)
    best_line = "No NMI head-to-head metric available."
    if not best.empty:
        r = best.iloc[0]
        best_line = f"{r['model']}: n={int(r['n'])}, positives={int(r['positives'])}, AUPRC={r['AUPRC']:.3f}, AUROC={r['AUROC']:.3f}."
    public_line = "No public comparator head-to-head metric available."
    if not best_public.empty:
        r = best_public.iloc[0]
        public_line = f"{r['model']}: n={int(r['n'])}, positives={int(r['positives'])}, AUPRC={r['AUPRC']:.3f}, AUROC={r['AUROC']:.3f}."

    md = [
        "# imNEO Public-Source Reconstruction Benchmark",
        "",
        "## Decision",
        "Compare only the public-source space disclosed by CG Invites/imNEO patent/presentation context against our confident clean algorithm family.",
        "",
        "## Brutal Claim Boundary",
        "- This is **not** CG Invites' proprietary peptide list.",
        "- This is **not** a reverse-engineered imNEO model.",
        "- This is a public-source reconstruction using local canonical rows that match the disclosed data-source names.",
        "- The fair head-to-head set is the common rows where NMI and frozen public comparators all have scores.",
        "- Patient-level top-N is **not claimable** in this reconstruction when public source rows collapse to unknown/synthetic patient IDs.",
        "",
        "## Why This Is The Right Narrow Comparison",
        "- It avoids pretending we have undisclosed company-selected peptides.",
        "- It directly tests our clean NMI branch against the same public-source universe visible from the patent/presentation trail.",
        "- It preserves the clean-method boundary: BigMHC/MHCflurry/PRIME/NetMHCpan are comparators, not NMI training features.",
        "",
        "## Coverage",
        source_inventory.to_markdown(index=False),
        "",
        "## Main Head-to-Head Result",
        f"- Best NMI: `{best_line}`",
        f"- Best frozen public comparator: `{public_line}`",
        "",
        "## Leaderboard",
        leaderboard.head(40).to_markdown(index=False) if not leaderboard.empty else "_No metrics generated._",
        "",
        "## Interpretation",
        "- If NMI wins on `imNEO_core_head_to_head_common`, the defensible claim is: our clean multimodal immunogenicity branch outperforms frozen public predictors on the currently reconstructed patent-like public-source overlap.",
        "- If NMI does not win, the defensible claim is: NMI coverage and source transfer need extension before competitor-facing claims.",
        "- NEPdb coverage is currently the cleanest bridge because it contains the NMI locked rows and matches one disclosed public source.",
        "- TESLA/dbPepNeo are present in canonical but do not yet have NMI locked-branch coverage; this is the next ingestion/feature-extension target.",
        "- Patient-level hit-rate columns are intentionally blank/NA unless at least two real patient IDs are present.",
        "",
        "## Output Files",
        f"- `metrics/imneo_public_reconstruction_leaderboard.tsv`",
        f"- `reports/imneo_public_source_inventory.tsv`",
        f"- `predictions/imneo_public_reconstruction_ranked_candidates.tsv`",
        f"- `figures/imneo_public_reconstruction_leaderboard.png`",
        f"- `figures/imneo_public_source_coverage.png`",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "IMNEO_PUBLIC_RECONSTRUCTION_BENCHMARK.md")
    print(outdir / "reports" / "IMNEO_PUBLIC_RECONSTRUCTION_BENCHMARK.md")


if __name__ == "__main__":
    main()
