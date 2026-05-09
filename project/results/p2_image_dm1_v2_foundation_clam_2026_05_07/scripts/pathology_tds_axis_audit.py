"""
Pathology scout 2: does the UNI image-DM1 score read thyroid differentiation?

The image-DM1 classifier itself is not used as a new DM1 claim here. Instead,
this audit asks whether the center-held-out UNI image score behaves as a
morphology readout of RNA-derived thyroid differentiation / dedifferentiation.

Critical confound check:
  Residualize both image score and dedifferentiation proxy against DM label,
  TSS, histology, molecular subtype, sex, and age. If the residual association
  survives, this is a candidate pathology finding. If it dies, it is just label
  leakage / cohort structure.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict


REPO = Path(os.environ.get("THCA_REPO_ROOT", Path.cwd()))
P2_ROOT = Path(
    os.environ.get(
        "THCA_P2_ROOT",
        REPO / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07",
    )
)
OUT_DIR = P2_ROOT / "analysis_supp/pathology_tds_axis_2026_05_09"
MANIFEST = P2_ROOT / "phase2_tcga_clam/slide_manifest.tsv"
UNI_LOTO = P2_ROOT / "analysis_supp/audit_uni_loto/uni_loto_oof_predictions.tsv"
UNI_OOF = P2_ROOT / "phase2_tcga_clam_UNI/clam_per_slide_predictions.tsv"
META = REPO / "project/metadata/sample_master_v3.tsv"


def safe_spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 4 or np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return float("nan"), float("nan")
    rho, p = spearmanr(x[mask], y[mask])
    return float(rho), float(p)


def case_barcode(sample_id: str) -> str:
    return str(sample_id)[:12]


def load_table() -> pd.DataFrame:
    man = pd.read_csv(MANIFEST, sep="\t")
    man = man[man["dm"].isin(["DM1", "DM2"])].copy()
    man["case_barcode"] = man["submitter_id"].str[:12]
    man["tss"] = man["submitter_id"].str.split("-").str[1]
    loto = pd.read_csv(UNI_LOTO, sep="\t").rename(columns={"slide": "file_id"})
    oof = pd.read_csv(UNI_OOF, sep="\t").rename(columns={"slide": "file_id", "prob_DM1": "prob_DM1_uni_oof"})
    meta = pd.read_csv(META, sep="\t")
    meta["case_barcode"] = meta["sample_id"].map(case_barcode)
    meta = meta[(meta["dataset"] == "TCGA-THCA") & (meta["normal_vs_tumor"] == "tumor")].copy()
    meta = meta.drop_duplicates("case_barcode", keep="first")
    cols = [
        "case_barcode",
        "sample_id",
        "histology_subtype",
        "driver_anchor",
        "molecular_subtype",
        "tds_group",
        "ajcc_stage_group",
        "sex",
        "age",
        "tds_score",
        "dedifferentiation_proxy_score",
        "v3_anchor",
    ]
    df = man.merge(loto[["file_id", "prob_DM1_uni_loto"]], on="file_id", how="inner")
    df = df.merge(oof[["file_id", "fold", "label", "prob_DM1_uni_oof"]], on="file_id", how="left")
    df = df.merge(meta[cols], on="case_barcode", how="left")
    df["dm_label"] = (df["dm"] == "DM1").astype(int)
    for col in ["prob_DM1_uni_loto", "prob_DM1_uni_oof", "tds_score", "dedifferentiation_proxy_score", "age"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def design_matrix(df: pd.DataFrame, covars: list[str]) -> pd.DataFrame:
    parts = []
    for cov in covars:
        if cov == "age":
            vals = pd.to_numeric(df["age"], errors="coerce")
            vals = vals.fillna(vals.median())
            parts.append(pd.DataFrame({"age": vals.astype(float)}, index=df.index))
        elif cov == "dm_label":
            parts.append(pd.DataFrame({"dm_label": df["dm_label"].astype(float)}, index=df.index))
        else:
            d = pd.get_dummies(df[cov].fillna("missing").astype(str), prefix=cov, drop_first=True)
            if d.shape[1]:
                parts.append(d)
    if not parts:
        return pd.DataFrame(index=df.index)
    return pd.concat(parts, axis=1)


def residualize(df: pd.DataFrame, value_col: str, covars: list[str]) -> np.ndarray:
    y = pd.to_numeric(df[value_col], errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(y)
    x = design_matrix(df, covars)
    resid = np.full(len(df), np.nan)
    if x.shape[1] == 0:
        resid[mask] = y[mask] - np.nanmean(y[mask])
        return resid
    if mask.sum() <= x.shape[1] + 3:
        return resid
    lr = LinearRegression()
    lr.fit(x.loc[mask].to_numpy(dtype=float), y[mask])
    resid[mask] = y[mask] - lr.predict(x.loc[mask].to_numpy(dtype=float))
    return resid


def spearman_rows(df: pd.DataFrame) -> pd.DataFrame:
    target = "dedifferentiation_proxy_score"
    score_cols = ["prob_DM1_uni_loto", "prob_DM1_uni_oof"]
    covar_sets = {
        "none": [],
        "dm_label": ["dm_label"],
        "tss": ["tss"],
        "tss_dm_label": ["tss", "dm_label"],
        "histology_molecular_sex_age": ["histology_subtype", "molecular_subtype", "sex", "age"],
        "full": ["dm_label", "tss", "histology_subtype", "molecular_subtype", "sex", "age"],
    }
    rows = []
    for score in score_cols:
        for label, covars in covar_sets.items():
            sub = df[[score, target] + [c for c in covars if c in df.columns or c == "dm_label"]].copy()
            sub = sub.dropna(subset=[score, target])
            if label == "none":
                xs = sub[score].to_numpy(dtype=float)
                ys = sub[target].to_numpy(dtype=float)
            else:
                xs = residualize(sub.assign(dm_label=df.loc[sub.index, "dm_label"]), score, covars)
                ys = residualize(sub.assign(dm_label=df.loc[sub.index, "dm_label"]), target, covars)
            rho, p = safe_spearman(xs, ys)
            rows.append(
                {
                    "score": score,
                    "target": target,
                    "adjustment": label,
                    "n": int(np.isfinite(xs).sum() if label != "none" else len(sub)),
                    "spearman_rho": rho,
                    "spearman_p": p,
                }
            )
    for score in score_cols:
        for group_col, group_val in [
            ("dm", "DM1"),
            ("dm", "DM2"),
            ("molecular_subtype", "BRAF_like"),
            ("molecular_subtype", "RAS_like"),
            ("histology_subtype", "cPTC"),
            ("histology_subtype", "FVPTC"),
            ("sex", "Female"),
            ("sex", "Male"),
        ]:
            sub = df[(df[group_col] == group_val) & df[score].notna() & df[target].notna()].copy()
            if len(sub) < 8:
                continue
            rho, p = safe_spearman(sub[score], sub[target])
            rows.append(
                {
                    "score": score,
                    "target": target,
                    "adjustment": f"within_{group_col}_{group_val}",
                    "n": int(len(sub)),
                    "spearman_rho": rho,
                    "spearman_p": p,
                }
            )
    return pd.DataFrame(rows)


def auc_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    valid = df[df["dedifferentiation_proxy_score"].notna()].copy()
    q66 = valid["dedifferentiation_proxy_score"].quantile(2 / 3)
    q33 = valid["dedifferentiation_proxy_score"].quantile(1 / 3)
    valid["dediff_top_tertile"] = (valid["dedifferentiation_proxy_score"] >= q66).astype(int)
    valid["dediff_bottom_tertile"] = (valid["dedifferentiation_proxy_score"] <= q33).astype(int)
    for score in ["prob_DM1_uni_loto", "prob_DM1_uni_oof"]:
        sub = valid[valid[score].notna()].copy()
        y = sub["dediff_top_tertile"].to_numpy()
        s = sub[score].to_numpy(dtype=float)
        rows.append(
            {
                "score": score,
                "target": "dediff_top_tertile_vs_rest",
                "n": int(len(sub)),
                "n_pos": int(y.sum()),
                "auc_pos_high": float(roc_auc_score(y, s)),
                "mannwhitney_p": float(
                    mannwhitneyu(s[y == 1], s[y == 0], alternative="two-sided").pvalue
                ),
            }
        )
        extreme = sub[(sub["dediff_top_tertile"] == 1) | (sub["dediff_bottom_tertile"] == 1)].copy()
        y2 = extreme["dediff_top_tertile"].to_numpy()
        s2 = extreme[score].to_numpy(dtype=float)
        rows.append(
            {
                "score": score,
                "target": "dediff_top_vs_bottom_tertile",
                "n": int(len(extreme)),
                "n_pos": int(y2.sum()),
                "auc_pos_high": float(roc_auc_score(y2, s2)),
                "mannwhitney_p": float(
                    mannwhitneyu(s2[y2 == 1], s2[y2 == 0], alternative="two-sided").pvalue
                ),
            }
        )
    # Clinical metadata-only baseline for top-tertile dediff, to benchmark image score.
    sub = valid.dropna(subset=["prob_DM1_uni_loto"]).copy()
    y = sub["dediff_top_tertile"].to_numpy()
    if len(np.unique(y)) == 2 and len(sub) >= 20:
        x_meta = design_matrix(sub, ["dm_label", "tss", "histology_subtype", "molecular_subtype", "sex", "age"])
        x_img = pd.concat([x_meta, sub[["prob_DM1_uni_loto"]].reset_index(drop=True)], axis=1)
        cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=9)
        for name, x in [("metadata_only", x_meta), ("metadata_plus_image", x_img)]:
            try:
                lr = LogisticRegression(max_iter=1000, solver="liblinear")
                prob = cross_val_predict(lr, x.to_numpy(dtype=float), y, cv=cv, method="predict_proba")[:, 1]
                rows.append(
                    {
                        "score": name,
                        "target": "dediff_top_tertile_oof_lr",
                        "n": int(len(sub)),
                        "n_pos": int(y.sum()),
                        "auc_pos_high": float(roc_auc_score(y, prob)),
                        "mannwhitney_p": "NA",
                    }
                )
            except Exception as exc:
                rows.append(
                    {
                        "score": name,
                        "target": "dediff_top_tertile_oof_lr",
                        "n": int(len(sub)),
                        "n_pos": int(y.sum()),
                        "auc_pos_high": float("nan"),
                        "mannwhitney_p": f"failed: {exc}",
                    }
                )
    return pd.DataFrame(rows)


def r2_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for value in ["prob_DM1_uni_loto", "dedifferentiation_proxy_score"]:
        sub = df[df[value].notna()].copy()
        y = sub[value].to_numpy(dtype=float)
        for covar_set, covars in {
            "dm_label": ["dm_label"],
            "tss": ["tss"],
            "histology_molecular": ["histology_subtype", "molecular_subtype"],
            "full_no_image": ["dm_label", "tss", "histology_subtype", "molecular_subtype", "sex", "age"],
        }.items():
            x = design_matrix(sub, covars)
            if x.shape[1] == 0 or len(sub) <= x.shape[1] + 3:
                continue
            lr = LinearRegression().fit(x.to_numpy(dtype=float), y)
            rows.append({"value": value, "explained_by": covar_set, "n": int(len(sub)), "r2": float(lr.score(x, y))})
    return pd.DataFrame(rows)


def plot(df: pd.DataFrame, spearman: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    ax = axes[0]
    colors = {"BRAF_like": "#2f6f8f", "RAS_like": "#7b1f2a", "unknown": "#777777"}
    for subtype, sub in df.groupby("molecular_subtype"):
        ax.scatter(
            sub["dedifferentiation_proxy_score"],
            sub["prob_DM1_uni_loto"],
            s=42,
            alpha=0.85,
            label=subtype,
            color=colors.get(subtype, "#999999"),
        )
    raw = spearman[
        (spearman["score"] == "prob_DM1_uni_loto")
        & (spearman["target"] == "dedifferentiation_proxy_score")
        & (spearman["adjustment"] == "none")
    ].iloc[0]
    ax.set_title(f"UNI LOTO vs RNA dediff proxy: rho={raw.spearman_rho:.2f}, p={raw.spearman_p:.1e}")
    ax.set_xlabel("RNA dedifferentiation proxy score")
    ax.set_ylabel("UNI LOTO prob_DM1")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    order = ["low", "mid", "high"]
    data = [df.loc[df["tds_group"] == g, "prob_DM1_uni_loto"].dropna() for g in order]
    ax.boxplot(data, labels=order, patch_artist=True)
    for i, vals in enumerate(data, start=1):
        x = np.full(len(vals), i) + np.linspace(-0.10, 0.10, len(vals))
        ax.scatter(x, vals, s=34, alpha=0.8, color="#2f6f8f")
    ax.set_title("UNI LOTO image score by TDS group")
    ax.set_xlabel("TDS group")
    ax.set_ylabel("UNI LOTO prob_DM1")

    ax = axes[2]
    sub = spearman[spearman["score"] == "prob_DM1_uni_loto"].copy()
    keep = ["none", "dm_label", "tss_dm_label", "histology_molecular_sex_age", "full"]
    sub = sub[sub["adjustment"].isin(keep)]
    ax.barh(sub["adjustment"], sub["spearman_rho"], color="#7b1f2a")
    ax.axvline(0, color="#999", lw=0.8)
    ax.set_title("Partial Spearman after residualization")
    ax.set_xlabel("rho with dediff proxy")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_pathology_tds_axis.png", dpi=180)
    fig.savefig(out_dir / "fig_pathology_tds_axis.pdf")
    plt.close(fig)


def write_report(summary: dict, out_dir: Path) -> None:
    lines = [
        "# Pathology scout 2 - UNI image score as thyroid dedifferentiation morphology axis",
        "",
        "## Verdict",
        "",
        f"- Verdict: **{summary['verdict']}**",
        f"- Raw Spearman with RNA dedifferentiation proxy: {summary['raw_rho']:.3f} (p={summary['raw_p']:.3g})",
        f"- After DM-label residualization: {summary['dm_label_resid_rho']:.3f} (p={summary['dm_label_resid_p']:.3g})",
        f"- After full residualization: {summary['full_resid_rho']:.3f} (p={summary['full_resid_p']:.3g})",
        f"- Top-vs-bottom dediff tertile AUC: {summary['top_bottom_auc']:.3f}",
        "",
        "## Interpretation",
        "",
    ]
    if summary["verdict"] == "CANDIDATE_PASS":
        lines.append(
            "The center-held-out UNI image score tracks RNA-derived dedifferentiation "
            "beyond the binary DM1 label and major metadata covariates. This is a "
            "candidate pathology finding: H&E morphology may carry a continuous "
            "thyroid-differentiation signal aligned with the Paper 1 DM1 axis."
        )
    elif summary["verdict"] == "RAW_ONLY":
        lines.append(
            "The raw image-score association with dedifferentiation is real, but it "
            "does not survive conservative residualization. Treat it as explanation "
            "for the image-DM1 model, not as a new independent finding."
        )
    else:
        lines.append(
            "No robust image-score/dedifferentiation association was found."
        )
    lines += [
        "",
        "## Output files",
        "",
        "- `pathology_tds_axis_merged.tsv`",
        "- `pathology_tds_axis_spearman.tsv`",
        "- `pathology_tds_axis_auc.tsv`",
        "- `pathology_tds_axis_covariate_r2.tsv`",
        "- `fig_pathology_tds_axis.png`",
    ]
    (out_dir / "PATHOLOGY_TDS_AXIS_REPORT.md").write_text("\n".join(lines) + "\n")


def pick_row(df: pd.DataFrame, score: str, adjustment: str) -> pd.Series:
    return df[(df["score"] == score) & (df["adjustment"] == adjustment)].iloc[0]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_table()
    spearman = spearman_rows(df)
    aucs = auc_rows(df)
    r2 = r2_rows(df)
    plot(df, spearman, OUT_DIR)

    df.to_csv(OUT_DIR / "pathology_tds_axis_merged.tsv", sep="\t", index=False, na_rep="NA")
    spearman.to_csv(OUT_DIR / "pathology_tds_axis_spearman.tsv", sep="\t", index=False, na_rep="NA")
    aucs.to_csv(OUT_DIR / "pathology_tds_axis_auc.tsv", sep="\t", index=False, na_rep="NA")
    r2.to_csv(OUT_DIR / "pathology_tds_axis_covariate_r2.tsv", sep="\t", index=False, na_rep="NA")

    raw = pick_row(spearman, "prob_DM1_uni_loto", "none")
    dm_resid = pick_row(spearman, "prob_DM1_uni_loto", "dm_label")
    full_resid = pick_row(spearman, "prob_DM1_uni_loto", "full")
    tb = aucs[
        (aucs["score"] == "prob_DM1_uni_loto")
        & (aucs["target"] == "dediff_top_vs_bottom_tertile")
    ].iloc[0]
    if raw.spearman_rho >= 0.45 and dm_resid.spearman_rho >= 0.30 and full_resid.spearman_rho >= 0.25:
        verdict = "CANDIDATE_PASS"
    elif raw.spearman_rho >= 0.35:
        verdict = "RAW_ONLY"
    else:
        verdict = "NO_GO"
    summary = {
        "verdict": verdict,
        "n": int(len(df)),
        "raw_rho": float(raw.spearman_rho),
        "raw_p": float(raw.spearman_p),
        "dm_label_resid_rho": float(dm_resid.spearman_rho),
        "dm_label_resid_p": float(dm_resid.spearman_p),
        "full_resid_rho": float(full_resid.spearman_rho),
        "full_resid_p": float(full_resid.spearman_p),
        "top_bottom_auc": float(tb.auc_pos_high),
        "top_bottom_p": float(tb.mannwhitney_p),
    }
    (OUT_DIR / "PATHOLOGY_TDS_AXIS_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_report(summary, OUT_DIR)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
