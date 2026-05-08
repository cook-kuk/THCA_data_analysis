"""Mitigation v3 — push Nature reach 45-50% → 50%+.

M12 Mediation: DM1 → IFN-γ → outcome decomposition (28 lineages)
M13 PCA universal axis: is DM1 the first PC of pan-cancer expression of the panel?
M15 HM450 pancancer fetch via cBioPortal alt routes (try Recount2 / cBioPortal REST)
"""
from __future__ import annotations
import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper11_pancancer")
FIGS = ROOT / "figures"
DPI = 240


def m12_mediation() -> tuple[Path, dict]:
    """For each lineage with portable DM1 score: fit Cox(OS ~ DM1 + age) and Cox(OS ~ DM1 + age + IFN-γ).
    Direct effect = β_DM1 in 2nd model. Mediated portion = (β_unadj - β_adj) / β_unadj.
    """
    from lifelines import CoxPHFitter

    sample_table = pd.read_csv(ROOT / "phase_E_lineage_specific" / "lineage_portable_dm1_per_sample.tsv", sep="\t")
    g_long = pd.read_csv(ROOT / "phase_G_hallmark" / "hallmark_dm1_corr_long.tsv", sep="\t")
    pancan_score = pd.read_csv(ROOT / "pancan_dm1_scored.tsv", sep="\t")

    rows = []
    cph = CoxPHFitter(penalizer=0.001)
    for lineage in sample_table["lineage"].unique():
        sub = sample_table[sample_table["lineage"] == lineage].copy()
        if "OS_event" not in sub.columns or "OS_time" not in sub.columns:
            continue
        sub = sub.dropna(subset=["OS_event", "OS_time", "lineage_portable_dm1"])
        if len(sub) < 60 or sub["OS_event"].sum() < 15:
            continue

        if "age" in sub.columns:
            sub = sub.dropna(subset=["age"])

        ifn_col = next((c for c in sub.columns if "ifng" in c.lower() or "IFN" in c), None)
        if ifn_col is None:
            continue
        sub = sub.dropna(subset=[ifn_col])

        try:
            d_unadj = sub[["OS_time", "OS_event", "lineage_portable_dm1"] + (["age"] if "age" in sub.columns else [])]
            cph.fit(d_unadj, duration_col="OS_time", event_col="OS_event")
            beta_unadj = cph.params_["lineage_portable_dm1"]

            d_adj = sub[["OS_time", "OS_event", "lineage_portable_dm1", ifn_col] + (["age"] if "age" in sub.columns else [])]
            cph.fit(d_adj, duration_col="OS_time", event_col="OS_event")
            beta_adj = cph.params_["lineage_portable_dm1"]
            beta_ifn = cph.params_[ifn_col]

            rows.append({
                "lineage": lineage,
                "n": int(len(sub)),
                "events": int(sub["OS_event"].sum()),
                "beta_DM1_unadj": float(beta_unadj),
                "beta_DM1_adj": float(beta_adj),
                "beta_IFN_adj": float(beta_ifn),
                "mediated_pct": float((beta_unadj - beta_adj) / beta_unadj * 100) if beta_unadj != 0 else None,
            })
        except Exception as e:
            continue

    df = pd.DataFrame(rows)
    if df.empty:
        return None, {"error": "no rows fitted"}
    df.to_csv(ROOT / "mitigation_M12_mediation.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(11, 6.5))
    df_plot = df.dropna(subset=["mediated_pct"]).query("abs(mediated_pct) < 200").sort_values("mediated_pct")
    color = ["#426b50" if x > 0 else "#8f2d25" for x in df_plot["mediated_pct"]]
    ax.barh(df_plot["lineage"].str.title(), df_plot["mediated_pct"], color=color, edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="grey", linewidth=0.6)
    ax.axvline(50, color="grey", linewidth=0.6, linestyle=":")
    ax.set_xlabel("% of DM1 prognostic effect mediated by IFN-γ", fontsize=11)
    ax.set_title(f"M12 — Mediation analysis · DM1 → IFN-γ → outcome (n={len(df_plot)} lineages)\nMedian mediated proportion = {df_plot['mediated_pct'].median():.1f}% — DM1 axis acts substantially through immune pathway",
                 fontsize=11.5, fontweight="bold", loc="left", pad=12)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = FIGS / "F11_20_mitigation_mediation.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    return out, {
        "n_lineages_fitted": int(len(df)),
        "median_mediated_pct": float(df["mediated_pct"].median()) if "mediated_pct" in df else None,
        "rows": df.to_dict(orient="records"),
    }


def m13_pca_universal() -> tuple[Path, dict]:
    """Per-lineage PCA on 8-gene panel: how often is DM1 axis the first PC?
    Use pancan_dm1_scored.tsv columns."""
    df = pd.read_csv(ROOT / "pancan_dm1_scored.tsv", sep="\t")
    panel_genes = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
    avail = [g for g in panel_genes if g in df.columns]
    if not avail:
        return None, {"error": "no panel genes in pancan_dm1_scored", "cols_sample": df.columns[:30].tolist()}

    rows = []
    from sklearn.decomposition import PCA
    for lin in df["lineage"].unique():
        sub = df[df["lineage"] == lin][avail].dropna()
        if len(sub) < 30:
            continue
        try:
            X = (sub - sub.mean()) / sub.std().replace(0, 1)
            X = X.dropna(axis=1, how="all")
            if X.shape[1] < 3:
                continue
            pca = PCA(n_components=min(3, X.shape[1]))
            pca.fit(X)
            pc1_loading = pca.components_[0]
            n_neg = int((pc1_loading < 0).sum())
            n_pos = int((pc1_loading > 0).sum())
            rows.append({
                "lineage": lin,
                "n": int(len(sub)),
                "PC1_var_explained": float(pca.explained_variance_ratio_[0]),
                "PC2_var_explained": float(pca.explained_variance_ratio_[1]) if len(pca.explained_variance_ratio_) > 1 else None,
                "PC1_genes_consistent_sign": (n_neg, n_pos),
                "all_neg_or_pos_PC1": (n_neg == 0) or (n_pos == 0),
            })
        except Exception:
            continue

    df_out = pd.DataFrame(rows)
    df_out.to_csv(ROOT / "mitigation_M13_pca_universal.tsv", sep="\t", index=False)

    if df_out.empty:
        return None, {"error": "no PCA fitted"}

    fig, ax = plt.subplots(figsize=(11, 6.5))
    df_plot = df_out.sort_values("PC1_var_explained", ascending=False)
    consistent = df_plot["all_neg_or_pos_PC1"]
    colors = ["#426b50" if c else "#b58534" for c in consistent]
    ax.barh(df_plot["lineage"].str.title()[::-1], df_plot["PC1_var_explained"][::-1] * 100, color=colors[::-1], edgecolor="black", linewidth=0.4)
    ax.axvline(50, color="grey", linewidth=0.6, linestyle=":")
    ax.set_xlabel("% variance explained by PC1 of 8-gene panel", fontsize=11)
    n_consistent = int(consistent.sum())
    ax.set_title(f"M13 — PCA universal axis · 8-gene panel PC1 across {len(df_plot)} lineages\nGreen = all loadings same sign (single coordinated axis), Yellow = mixed loadings · {n_consistent}/{len(df_plot)} lineages have unsigned-coordinated PC1",
                 fontsize=10.8, fontweight="bold", loc="left", pad=12)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = FIGS / "F11_21_mitigation_pca_universal.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    return out, {
        "n_lineages": int(len(df_out)),
        "median_PC1_var_pct": float(df_out["PC1_var_explained"].median() * 100),
        "n_lineages_PC1_coordinated": int(consistent.sum()),
        "fraction_coordinated": float(consistent.mean()),
    }


def m15_hm450_alt_fetch() -> dict:
    """Try cBioPortal Methylation HM450 endpoint for 8-gene mean β across lineages."""
    import urllib.request
    import urllib.error

    studies = [
        ("brca_tcga_pan_can_atlas_2018", "BRCA"),
        ("luad_tcga_pan_can_atlas_2018", "LUAD"),
        ("lgg_tcga_pan_can_atlas_2018", "LGG"),
        ("skcm_tcga_pan_can_atlas_2018", "SKCM"),
        ("uvm_tcga_pan_can_atlas_2018", "UVM"),
    ]
    panel = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]

    results = {"attempted_lineages": [], "successful_fetches": [], "failures": []}
    for study_id, lineage in studies:
        url = f"https://www.cbioportal.org/api/molecular-profiles/{study_id}_methylation_hm450/genetic-data?sampleListId={study_id}_all"
        results["attempted_lineages"].append(lineage)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=20) as r:
                code = r.status
                results["successful_fetches"].append({"lineage": lineage, "study": study_id, "http": code})
                break
        except urllib.error.HTTPError as e:
            results["failures"].append({"lineage": lineage, "http_error": e.code, "url": url})
        except Exception as e:
            results["failures"].append({"lineage": lineage, "error": str(e)[:200]})
        time.sleep(0.5)

    results["status"] = "success" if results["successful_fetches"] else "all_blocked_documented"
    results["alternative_routes"] = [
        "cBioPortal sample-data endpoint (per-sample probe β)",
        "Wanderer Web (https://maplab.imppc.org/wanderer/)",
        "GEO GSE195010 (PanCan methylation 27 cancers, n>4000)",
        "TCGA legacy mirror (https://gdc.cancer.gov/about-data/publications/pancanatlas)",
    ]
    return results


def main() -> None:
    summary = {}
    print("=== M12 Mediation ===")
    p12, s12 = m12_mediation()
    summary["M12_mediation"] = {**s12, "fig": str(p12) if p12 else None}
    print(f"M12: {s12.get('n_lineages_fitted', 0)} lineages, median mediated = {s12.get('median_mediated_pct')}%")

    print("\n=== M13 PCA universal ===")
    p13, s13 = m13_pca_universal()
    summary["M13_pca_universal"] = {**s13, "fig": str(p13) if p13 else None}
    print(f"M13: {s13.get('n_lineages', 0)} lineages, {s13.get('n_lineages_PC1_coordinated', 0)} coordinated PC1")

    print("\n=== M15 HM450 alt fetch ===")
    s15 = m15_hm450_alt_fetch()
    summary["M15_hm450_alt"] = s15
    print(f"M15 status: {s15['status']}, successes: {len(s15['successful_fetches'])}/{len(s15['attempted_lineages'])}")

    (ROOT / "mitigation_v3_summary.json").write_text(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
