"""v15: K2 panel-only profile overlay.

K2/PRJEB11591 has an 8-gene kallisto mini-index quantification but no MAPK
genes, so it cannot be added to the v14 MAPK x Panel-8 forest. This script
adds a conservative score-level overlay:

  * fit the existing scale-invariant, within-sample-centered 8-gene TCGA model;
  * score TCGA, Lee/GSE213647, and K2 on the same centered profile axis;
  * report K2 as a panel-only reserve row, not a MAPK mechanism row.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/p_deconv_2026_05_08"
HUB = ROOT / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/paper1"

PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]


def center_rows(x: np.ndarray) -> np.ndarray:
    return x - x.mean(axis=1, keepdims=True)


def gene_z_mean(frame: pd.DataFrame, genes: list[str]) -> pd.Series:
    sub = frame.loc[genes].T.astype(float)
    z = (sub - sub.mean(axis=0)) / (sub.std(axis=0, ddof=0) + 1e-9)
    return z.mean(axis=1)


def tag_k2_sample(alias: object) -> tuple[str, str]:
    if not isinstance(alias, str):
        return "unknown", "unknown"
    if "-N" in alias:
        return "Normal", "normal"
    if "FA" in alias:
        return "FA", "tumor"
    if "FT" in alias:
        return "FTC", "tumor"
    if "FV" in alias:
        return "FVPTC", "tumor"
    if "PT" in alias:
        return "PTC", "tumor"
    return "unknown", "unknown"


def score_matrix(
    mat: pd.DataFrame,
    cohort: str,
    group: pd.Series,
    model: LogisticRegression,
    scaler: StandardScaler,
    source: str,
    extra: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """mat is genes x samples, log2-like expression."""
    have = [g for g in PANEL_8 if g in mat.index]
    if len(have) != len(PANEL_8):
        missing = sorted(set(PANEL_8) - set(have))
        raise RuntimeError(f"{cohort}: missing panel genes: {missing}")
    x = mat.loc[PANEL_8].T.astype(float)
    x_centered = center_rows(x.values)
    p_dm2 = model.predict_proba(scaler.transform(x_centered))[:, 1]
    out = pd.DataFrame(
        {
            "cohort": cohort,
            "sample_id": x.index,
            "group": group.reindex(x.index).fillna("unknown").values,
            "p_DM2_centered": p_dm2,
            "DM_call_centered": np.where(p_dm2 >= 0.5, "DM2", "DM1"),
            "panel_z_mean": gene_z_mean(mat.loc[PANEL_8], PANEL_8).reindex(x.index).values,
            "source": source,
        }
    )
    if extra is not None:
        extra_df = extra.reset_index()
        if "sample_id" not in extra_df.columns:
            extra_df = extra_df.rename(columns={extra_df.columns[0]: "sample_id"})
        out = out.merge(extra_df, on="sample_id", how="left")
    return out


def read_lee_symbol_log2() -> tuple[pd.DataFrame, pd.DataFrame]:
    lee = pd.read_csv(
        "/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_log2.tsv",
        sep="\t",
        index_col=0,
    )
    gmap = pd.read_csv("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv", sep="\t")
    ens_to_sym = dict(zip(gmap["ensembl"], gmap["symbol"]))
    ens_target = [e for e, s in ens_to_sym.items() if s in PANEL_8]
    sub = lee.loc[lee.index.intersection(ens_target)].copy()
    sub.index = [ens_to_sym[e] for e in sub.index]
    sub = sub.groupby(sub.index).mean()
    meta = pd.read_csv("/data/thca/v17_korean/GSE213647/sample_sheet_clinical.tsv", sep="\t")
    meta = meta.set_index("gsm")
    meta_small = meta[["histology", "tissue_type", "cell_type", "cell_subtype"]].copy()
    meta_small["tumor_normal"] = np.where(meta_small["tissue_type"].str.lower().eq("normal"), "normal", "tumor")
    return sub, meta_small


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cohort, group), sub in scores.groupby(["cohort", "group"], dropna=False):
        vals = sub["p_DM2_centered"].dropna()
        rows.append(
            {
                "cohort": cohort,
                "group": group,
                "n": int(len(vals)),
                "mean_p_DM2_centered": float(vals.mean()),
                "median_p_DM2_centered": float(vals.median()),
                "sd_p_DM2_centered": float(vals.std(ddof=1)) if len(vals) > 1 else np.nan,
                "q25_p_DM2_centered": float(vals.quantile(0.25)),
                "q75_p_DM2_centered": float(vals.quantile(0.75)),
                "dm2_call_fraction": float((vals >= 0.5).mean()),
                "mean_panel_z": float(sub["panel_z_mean"].mean()),
                "median_panel_z": float(sub["panel_z_mean"].median()),
            }
        )
    return pd.DataFrame(rows)


def mann_whitney(a: pd.Series, b: pd.Series) -> dict[str, float]:
    a = a.dropna()
    b = b.dropna()
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "mannwhitney_u": float(u),
        "mannwhitney_p": float(p),
    }


def plot_outputs(scores: pd.DataFrame, summary: pd.DataFrame, metrics: dict[str, object]) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.22,
        }
    )
    fig = plt.figure(figsize=(13.4, 9.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    order = [
        ("TCGA-THCA", "TCGA_DM1"),
        ("TCGA-THCA", "TCGA_DM2"),
        ("Lee_GSE213647", "Lee_normal"),
        ("Lee_GSE213647", "Lee_PTC_tumor"),
        ("Lee_GSE213647", "Lee_ATC_PD_tumor"),
        ("K2_PRJEB11591", "K2_tumor"),
        ("K2_PRJEB11591", "K2_normal"),
    ]
    labels = ["TCGA\nDM1", "TCGA\nDM2", "Lee\nnormal", "Lee\nPTC", "Lee\nATC/PD", "K2\ntumor", "K2\nnormal"]
    data = [
        scores[(scores["cohort"] == c) & (scores["group"] == g)]["p_DM2_centered"].dropna().values
        for c, g in order
    ]
    colors = ["#c44e52", "#4c72b0", "#55a868", "#8172b3", "#dd8452", "#ccb974", "#64b5cd"]
    bp = ax1.boxplot(data, patch_artist=True, showfliers=False, widths=0.64)
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.55)
        patch.set_edgecolor("#2c3440")
    for i, vals in enumerate(data, start=1):
        if len(vals) == 0:
            continue
        rng = np.random.default_rng(100 + i)
        x = rng.normal(i, 0.045, size=len(vals))
        ax1.scatter(x, vals, s=8, color=colors[i - 1], alpha=0.36, linewidths=0)
    ax1.axhline(0.5, color="#222", lw=1, ls="--")
    ax1.set_ylim(-0.03, 1.03)
    ax1.set_xticks(range(1, len(labels) + 1))
    ax1.set_xticklabels(labels)
    ax1.set_ylabel("profile-centered p(DM2)")
    ax1.set_title("A. Panel-only score overlay; K2 cannot enter MAPK forest", loc="left", fontweight="bold")

    k2 = scores[scores["cohort"] == "K2_PRJEB11591"].copy()
    for group, col in [("K2_tumor", "#ccb974"), ("K2_normal", "#64b5cd")]:
        vals = k2[k2["group"] == group]["p_DM2_centered"].dropna()
        ax2.hist(vals, bins=np.linspace(0, 1, 21), alpha=0.65, color=col, label=f"{group} n={len(vals)}")
    ax2.axvline(0.5, color="#222", lw=1, ls="--")
    ax2.set_xlabel("profile-centered p(DM2)")
    ax2.set_ylabel("sample count")
    ax2.legend(frameon=False)
    ax2.set_title("B. K2 mini-index distribution", loc="left", fontweight="bold")

    heat_groups = [
        ("TCGA_DM1", scores[(scores["cohort"] == "TCGA-THCA") & (scores["group"] == "TCGA_DM1")]),
        ("TCGA_DM2", scores[(scores["cohort"] == "TCGA-THCA") & (scores["group"] == "TCGA_DM2")]),
        ("Lee_PTC", scores[(scores["cohort"] == "Lee_GSE213647") & (scores["group"] == "Lee_PTC_tumor")]),
        ("K2_tumor", scores[(scores["cohort"] == "K2_PRJEB11591") & (scores["group"] == "K2_tumor")]),
    ]
    # Rebuild per-gene centered medians from the stored raw score files.
    tcga_raw = pd.read_csv(OUT / "v15_k2_panel_overlay_centered_gene_profiles.tsv", sep="\t")
    hm = []
    for group, _sub in heat_groups:
        temp = tcga_raw[tcga_raw["group"] == group]
        hm.append([temp[g].median() for g in PANEL_8])
    hm = np.array(hm)
    im = ax3.imshow(hm, aspect="auto", cmap="RdBu_r", vmin=-2.2, vmax=2.2)
    ax3.set_xticks(range(len(PANEL_8)))
    ax3.set_xticklabels(PANEL_8, rotation=45, ha="right")
    ax3.set_yticks(range(len(heat_groups)))
    ax3.set_yticklabels([x[0] for x in heat_groups])
    ax3.set_title("C. Median centered 8-gene profile", loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax3, shrink=0.8, label="TCGA-scaled centered feature")

    bar_df = summary.set_index(["cohort", "group"]).reindex(order).reset_index()
    ax4.bar(range(len(bar_df)), bar_df["dm2_call_fraction"], color=colors, alpha=0.72)
    ax4.set_ylim(0, 1.02)
    ax4.set_xticks(range(len(labels)))
    ax4.set_xticklabels(labels)
    ax4.set_ylabel("fraction called DM2")
    ax4.set_title("D. Threshold summary (p>=0.5)", loc="left", fontweight="bold")
    for i, row in bar_df.iterrows():
        if pd.isna(row["dm2_call_fraction"]):
            continue
        ax4.text(i, row["dm2_call_fraction"] + 0.025, f"n={int(row['n'])}", ha="center", va="bottom", fontsize=7)

    fig.suptitle(
        "Supp Fig SX_v15 -- K2 PRJEB11591 8-gene mini-index overlay (panel-only, MAPK unavailable)",
        fontsize=14,
        fontweight="bold",
    )
    for ext in ["png", "pdf"]:
        fig.savefig(OUT / f"Fig_SX_v15_K2_panel_overlay.{ext}", dpi=220 if ext == "png" else None)
        fig.savefig(ASSET / f"Fig_SX_v15_K2_panel_overlay.{ext}", dpi=220 if ext == "png" else None)
    plt.close(fig)


def write_html(summary: pd.DataFrame, metrics: dict[str, object]) -> None:
    rows = "\n".join(
        f"<tr><td>{r.cohort}</td><td>{r.group}</td><td>{int(r.n)}</td>"
        f"<td>{r.mean_p_DM2_centered:.3f}</td><td>{r.median_p_DM2_centered:.3f}</td>"
        f"<td>{r.dm2_call_fraction:.3f}</td><td>{r.mean_panel_z:.3f}</td></tr>"
        for r in summary.itertuples(index=False)
    )
    html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Paper 1 v15 K2 Panel Overlay</title>
<style>
body{{margin:0;background:#081019;color:#edf5ff;font-family:Inter,"Noto Sans KR",Arial,sans-serif;line-height:1.55}}
a{{color:#44d3ad;text-decoration:none}}a:hover{{text-decoration:underline}}code,.path{{font-family:"JetBrains Mono",monospace}}
.hero{{padding:54px 32px 34px;background:linear-gradient(135deg,#0a1320 0%,#16263b 62%,#241317 100%);border-bottom:1px solid #26364d}}
.inner{{max-width:1180px;margin:0 auto}}.k{{font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.18em;color:#ffd28a;text-transform:uppercase;font-weight:800}}
h1{{font-family:"Cormorant Garamond",Georgia,serif;font-size:48px;line-height:1;margin:10px 0;color:#fff6dd}}.lead{{font-family:Georgia,serif;font-size:18px;color:#d7e1ee;max-width:940px}}
.stats{{display:grid;grid-template-columns:repeat(6,1fr);gap:10px;margin-top:20px}}.stat{{border:1px solid rgba(255,210,138,.3);background:rgba(255,255,255,.045);border-radius:8px;padding:12px}}.stat b{{display:block;font-family:"Cormorant Garamond",Georgia,serif;font-size:28px;color:#ffd28a;line-height:1}}.stat span{{font-family:"JetBrains Mono",monospace;font-size:9px;color:#9fb0c5;text-transform:uppercase;letter-spacing:.06em}}
main{{max-width:1180px;margin:0 auto;padding:28px 32px 60px}}section{{border-bottom:1px solid #26364d;padding:24px 0}}h2{{font-family:"Cormorant Garamond",Georgia,serif;font-size:31px;color:#fff6dd}}.box{{border-left:4px solid #ffd28a;background:#101a28;border-radius:0 8px 8px 0;padding:14px 16px;margin:12px 0}}.warn{{border-left-color:#ff8a6b}}table{{border-collapse:collapse;width:100%;font-size:12.5px}}th,td{{border:1px solid #26364d;padding:7px 9px;text-align:left}}th{{background:#17253d;color:#ffd28a;font-family:"JetBrains Mono",monospace;font-size:10px;text-transform:uppercase}}td{{color:#d7e1ed}}img{{max-width:100%;background:white;border:1px solid #26364d;border-radius:8px}}.path{{background:#0a111d;border:1px solid #26364d;border-radius:4px;padding:1px 5px;color:#44d3ad;font-size:11px}}
@media(max-width:850px){{.stats{{grid-template-columns:repeat(2,1fr)}}h1{{font-size:36px}}main,.hero{{padding-left:18px;padding-right:18px}}}}
</style></head><body>
<header class="hero"><div class="inner"><div class="k">Paper 1 mechanism v15 · K2 panel-only overlay · 2026-05-09</div>
<h1>K2 is a panel-only calibration stress test, not a MAPK validation cohort.</h1>
<p class="lead">PRJEB11591/K2 has local 8-gene mini-index quantification for 260 runs. K2 tumor runs score DM2-like, but K2 normal runs also score high, so this is a cross-platform score overlay and calibration audit only. Because MAPK-output genes are absent, no causal, tumor-normal, prevalence, or MAPK-correlation claim is made from K2.</p>
<div class="stats">
<div class="stat"><b>{metrics['k2_n_total']}</b><span>K2 scored runs</span></div>
<div class="stat"><b>{metrics['k2_n_tumor']}</b><span>K2 tumor runs</span></div>
<div class="stat"><b>{metrics['k2_tumor_dm2_fraction']:.2f}</b><span>K2 tumor DM2 fraction</span></div>
<div class="stat"><b>{metrics['k2_normal_median_p_DM2']:.3f}</b><span>K2 normal median p(DM2)</span></div>
<div class="stat"><b>{metrics['tcga_cv_auc_mean']:.3f}</b><span>TCGA 5-fold AUC</span></div>
<div class="stat"><b>0</b><span>MAPK genes in K2</span></div>
</div></div></header>
<main>
<section><h2>Disposition</h2>
<div class="box"><b>Use:</b> reviewer-reserve Korean 8-gene profile overlay and calibration audit.</div>
<div class="box warn"><b>Do not use:</b> as a sixth MAPK x Panel-8 forest cohort. K2 mini-index contains only the eight thyroid panel genes.</div>
<div class="box warn"><b>Do not infer prevalence or tumor-normal discrimination:</b> K2 normal samples also have high p(DM2), indicating platform/mini-index baseline shift under this projection.</div>
</section>
<section><h2>Figure</h2>
<img src="assets/paper1/Fig_SX_v15_K2_panel_overlay.png" alt="K2 panel overlay figure"></section>
<section><h2>Group Summary</h2>
<table><thead><tr><th>Cohort</th><th>Group</th><th>n</th><th>mean p(DM2)</th><th>median p(DM2)</th><th>DM2 fraction</th><th>mean panel z</th></tr></thead><tbody>
{rows}
</tbody></table></section>
<section><h2>Source Paths</h2>
<p><span class="path">project/results/p_deconv_2026_05_08/run_v15_k2_panel_overlay.py</span></p>
<p><span class="path">project/results/p_deconv_2026_05_08/v15_k2_panel_overlay_scores.tsv</span></p>
<p><span class="path">project/results/p_deconv_2026_05_08/v15_k2_panel_overlay_summary.tsv</span></p>
<p><span class="path">/data/thca/repo_results/v17_korean/K2_8gene_tpm_matrix_v4.tsv</span></p>
</section></main></body></html>"""
    (HUB / "paper1_k2_panel_overlay_v15.html").write_text(html)


def main() -> None:
    ASSET.mkdir(parents=True, exist_ok=True)

    # TCGA training matrix and labels.
    tcga = pd.read_csv(
        "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
        sep="\t",
        index_col=0,
    )
    labels = pd.read_csv(ROOT / "project/results/v17_realfix/R1A_cluster_labels.tsv", sep="\t")
    labels["dm_label"] = np.where(labels["cluster"].str.startswith("DM2"), "TCGA_DM2", "TCGA_DM1")
    labels["y_dm2"] = labels["cluster"].str.startswith("DM2").astype(int)
    labels = labels.set_index("sample_id")
    common = [s for s in labels.index if s in tcga.columns]
    x_train = tcga.loc[PANEL_8, common].T.astype(float)
    x_train_centered = center_rows(x_train.values)
    y = labels.loc[common, "y_dm2"].values
    scaler = StandardScaler().fit(x_train_centered)
    model = LogisticRegression(C=1.0, max_iter=2000, random_state=42).fit(scaler.transform(x_train_centered), y)
    cv = cross_val_score(
        LogisticRegression(C=1.0, max_iter=2000, random_state=42),
        scaler.transform(x_train_centered),
        y,
        cv=StratifiedKFold(5, shuffle=True, random_state=42),
        scoring="roc_auc",
    )
    train_auc = roc_auc_score(y, model.predict_proba(scaler.transform(x_train_centered))[:, 1])

    tcga_group = labels.loc[tcga.columns.intersection(labels.index), "dm_label"]
    tcga_scores = score_matrix(
        tcga.loc[PANEL_8, tcga_group.index],
        "TCGA-THCA",
        tcga_group,
        model,
        scaler,
        "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
        extra=labels.loc[tcga_group.index, ["cluster"]].rename(columns={"cluster": "histology"}),
    )
    tcga_scores["tumor_normal"] = "tumor"

    lee, lee_meta = read_lee_symbol_log2()
    lee_group = pd.Series("Lee_other", index=lee.columns)
    lee_group.loc[lee_meta.index[lee_meta["tumor_normal"].eq("normal")].intersection(lee.columns)] = "Lee_normal"
    lee_group.loc[
        lee_meta.index[(lee_meta["tumor_normal"].eq("tumor")) & (lee_meta["histology"].eq("PTC"))].intersection(lee.columns)
    ] = "Lee_PTC_tumor"
    lee_group.loc[
        lee_meta.index[
            (lee_meta["tumor_normal"].eq("tumor")) & (lee_meta["histology"].isin(["ATC", "PD"]))
        ].intersection(lee.columns)
    ] = "Lee_ATC_PD_tumor"
    lee_scores = score_matrix(
        lee.loc[PANEL_8],
        "Lee_GSE213647",
        lee_group,
        model,
        scaler,
        "/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_log2.tsv",
        extra=lee_meta,
    )

    k2 = pd.read_csv("/data/thca/repo_results/v17_korean/K2_8gene_tpm_matrix_v4.tsv", sep="\t").set_index("run")
    k2_meta = pd.read_csv("/data/thca/repo_results/v17_korean/K1A_prjeb11591_runs.tsv", sep="\t").set_index("run_accession")
    tags = k2_meta["sample_alias"].apply(tag_k2_sample)
    k2_meta["histology"] = [x[0] for x in tags]
    k2_meta["tumor_normal"] = [x[1] for x in tags]
    k2_extra = k2_meta[["sample_alias", "histology", "tumor_normal", "read_count", "base_count"]]
    k2_group = pd.Series("K2_unknown", index=k2.index)
    k2_group.loc[k2_extra.index[k2_extra["tumor_normal"].eq("tumor")].intersection(k2.index)] = "K2_tumor"
    k2_group.loc[k2_extra.index[k2_extra["tumor_normal"].eq("normal")].intersection(k2.index)] = "K2_normal"
    k2_log = np.log2(k2[PANEL_8].T + 1.0)
    k2_scores = score_matrix(
        k2_log,
        "K2_PRJEB11591",
        k2_group,
        model,
        scaler,
        "/data/thca/repo_results/v17_korean/K2_8gene_tpm_matrix_v4.tsv",
        extra=k2_extra,
    )

    scores = pd.concat([tcga_scores, lee_scores, k2_scores], ignore_index=True)
    scores = scores[[c for c in scores.columns if c != "source"] + ["source"]]
    scores.to_csv(OUT / "v15_k2_panel_overlay_scores.tsv", sep="\t", index=False)

    # Store TCGA-scaled centered gene profiles for heatmap audit.
    profile_rows = []
    for source_scores, mat in [
        (tcga_scores, tcga.loc[PANEL_8, tcga_scores["sample_id"]]),
        (lee_scores, lee.loc[PANEL_8, lee_scores["sample_id"]]),
        (k2_scores, k2_log.loc[PANEL_8, k2_scores["sample_id"]]),
    ]:
        x = mat.T.astype(float)
        centered_scaled = scaler.transform(center_rows(x.values))
        temp = pd.DataFrame(centered_scaled, columns=PANEL_8)
        temp.insert(0, "sample_id", x.index.values)
        temp = temp.merge(source_scores[["sample_id", "cohort", "group"]], on="sample_id", how="left")
        profile_rows.append(temp)
    profiles = pd.concat(profile_rows, ignore_index=True)
    profiles.to_csv(OUT / "v15_k2_panel_overlay_centered_gene_profiles.tsv", sep="\t", index=False)

    summary = summarize(scores)
    summary.to_csv(OUT / "v15_k2_panel_overlay_summary.tsv", sep="\t", index=False)

    k2_tumor = scores[(scores["cohort"].eq("K2_PRJEB11591")) & (scores["group"].eq("K2_tumor"))]
    k2_normal = scores[(scores["cohort"].eq("K2_PRJEB11591")) & (scores["group"].eq("K2_normal"))]
    tcga_dm2 = scores[(scores["cohort"].eq("TCGA-THCA")) & (scores["group"].eq("TCGA_DM2"))]
    tcga_dm1 = scores[(scores["cohort"].eq("TCGA-THCA")) & (scores["group"].eq("TCGA_DM1"))]
    metrics: dict[str, object] = {
        "method": "TCGA-trained scale-invariant within-sample-centered 8-gene LogisticRegression; p_DM2 is score-level overlay only for K2.",
        "panel_genes": PANEL_8,
        "tcga_training_n": int(len(common)),
        "tcga_train_auc_overfit_biased": float(train_auc),
        "tcga_cv_auc_mean": float(cv.mean()),
        "tcga_cv_auc_sd": float(cv.std()),
        "tcga_cv_auc_min": float(cv.min()),
        "tcga_cv_auc_max": float(cv.max()),
        "k2_n_total": int(len(k2_scores)),
        "k2_n_tumor": int(len(k2_tumor)),
        "k2_n_normal": int(len(k2_normal)),
        "k2_tumor_dm2_fraction": float((k2_tumor["p_DM2_centered"] >= 0.5).mean()),
        "k2_normal_dm2_fraction": float((k2_normal["p_DM2_centered"] >= 0.5).mean()) if len(k2_normal) else None,
        "k2_tumor_median_p_DM2": float(k2_tumor["p_DM2_centered"].median()),
        "k2_tumor_mean_p_DM2": float(k2_tumor["p_DM2_centered"].mean()),
        "k2_normal_median_p_DM2": float(k2_normal["p_DM2_centered"].median()) if len(k2_normal) else None,
        "k2_vs_tcga_dm2": mann_whitney(k2_tumor["p_DM2_centered"], tcga_dm2["p_DM2_centered"]),
        "k2_vs_tcga_dm1": mann_whitney(k2_tumor["p_DM2_centered"], tcga_dm1["p_DM2_centered"]),
        "limitation": "K2 mini-index contains only Panel-8 genes; MAPK-output genes are unavailable, so K2 is not a v14 forest cohort.",
    }
    with open(OUT / "v15_k2_panel_overlay_metrics.json", "w") as fh:
        json.dump(metrics, fh, indent=2)

    plot_outputs(scores, summary, metrics)
    write_html(summary, metrics)

    print(json.dumps(metrics, indent=2))
    print(f"Wrote {OUT/'v15_k2_panel_overlay_scores.tsv'}")
    print(f"Wrote {OUT/'Fig_SX_v15_K2_panel_overlay.png'}")
    print(f"Wrote {HUB/'paper1_k2_panel_overlay_v15.html'}")


if __name__ == "__main__":
    main()
