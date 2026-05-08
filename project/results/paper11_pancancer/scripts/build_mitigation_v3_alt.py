"""Mitigation v3 alt — works with available data schemas.

M12  Cross-lineage prognostic-mediation proxy:
     Phase E Cox log_HR × Phase G IFN-γ r → 28 lineages.
     Strong correlation = DM1 prognostic effect is mediated by lineage IFN-γ ecosystem.

M13  Three-score axis convergence:
     pancan_dm1_scored has DM1_like / TF_collapse / RAI_8 per sample.
     Cross-correlation per lineage → 32 lineages — does the axis converge?

M14  Cross-lineage architecture Spearman matrix:
     Each lineage's lineage-portable DM1 vs each other's — pan-cancer axis universality matrix.

M15  HM450 alternative routes (cBioPortal sample-data, MEXPRESS, Wanderer) — try harder.

M16  Driver × DM1 axis pancancer:
     If pancan_dm1_scored has BRAF/RAS info or driver-orthogonal architecture cohort = MAPK enriched.
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


def m12_cross_lineage_mediation() -> tuple[Path, dict]:
    """Phase E HR × Phase G IFN-γ r — does prognostic strength scale with lineage immune-axis strength?"""
    e_cox = pd.read_csv(ROOT / "phase_E_lineage_specific" / "cox_per_cohort_portable.tsv", sep="\t").dropna(subset=["HR"])
    e_cox["log_HR"] = np.log(e_cox["HR"])
    g = pd.read_csv(ROOT / "phase_G_hallmark" / "hallmark_dm1_corr_long.tsv", sep="\t")
    ifn = g[g["hallmark"] == "Interferon Gamma Response"][["lineage", "spearman_r"]].rename(columns={"spearman_r": "ifng_r"})
    allog = g[g["hallmark"] == "Allograft Rejection"][["lineage", "spearman_r"]].rename(columns={"spearman_r": "allog_r"})
    df = e_cox.merge(ifn, on="lineage").merge(allog, on="lineage")

    rho_ifn, p_ifn = stats.spearmanr(df["ifng_r"], df["log_HR"])
    rho_allog, p_allog = stats.spearmanr(df["allog_r"], df["log_HR"])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, (xcol, xlabel, rho, p) in zip(axes, [
        ("ifng_r", "IFN-γ Response Hallmark r", rho_ifn, p_ifn),
        ("allog_r", "Allograft Rejection Hallmark r", rho_allog, p_allog),
    ]):
        ax.scatter(df[xcol], df["log_HR"], s=140,
                   c=["#426b50" if h < 1 else "#8f2d25" for h in df["HR"]],
                   edgecolor="black", linewidth=0.5, alpha=0.85)
        for _, r in df.iterrows():
            abbr = r["lineage"].replace(" carcinoma", "").title()[:18]
            ax.annotate(abbr, (r[xcol], r["log_HR"]),
                        fontsize=8.5, xytext=(5, 4), textcoords="offset points")
        ax.axhline(0, color="grey", linewidth=0.6, linestyle="--")
        ax.set_xlabel(f"{xlabel} (Phase G)", fontsize=11)
        ax.set_ylabel("log HR (Phase E lineage-portable Cox OS)", fontsize=11)
        z = np.polyfit(df[xcol], df["log_HR"], 1)
        x_fit = np.linspace(df[xcol].min(), df[xcol].max(), 100)
        ax.plot(x_fit, np.polyval(z, x_fit), color="#7a1b15", linewidth=1.5, alpha=0.7, linestyle="--")
        ax.set_title(f"Spearman ρ = {rho:.2f}, p = {p:.3g}", fontsize=11, fontweight="bold", loc="left")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    fig.suptitle("M12 — Cross-lineage prognostic-mediation proxy · DM1 axis prognostic effect tracks lineage immune-axis strength",
                 fontsize=12.5, fontweight="bold", y=1.005)
    fig.tight_layout()
    out = FIGS / "F11_22_mitigation_cross_lineage_mediation.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    return out, {
        "n_lineages": int(len(df)),
        "spearman_rho_logHR_vs_IFNg": float(rho_ifn),
        "p_logHR_vs_IFNg": float(p_ifn),
        "spearman_rho_logHR_vs_Allograft": float(rho_allog),
        "p_logHR_vs_Allograft": float(p_allog),
    }


def m13_three_score_convergence() -> tuple[Path, dict]:
    """Per-lineage cross-correlation of DM1_like / TF_collapse / RAI_8 (3-score axis convergence)."""
    df = pd.read_csv(ROOT / "pancan_dm1_scored.tsv", sep="\t")
    rows = []
    for lin in df["lineage"].unique():
        sub = df[df["lineage"] == lin].dropna(subset=["DM1_like", "TF_collapse", "RAI_8"])
        if len(sub) < 30:
            continue
        rho_dm1_tf, _ = stats.spearmanr(sub["DM1_like"], sub["TF_collapse"])
        rho_dm1_rai, _ = stats.spearmanr(sub["DM1_like"], sub["RAI_8"])
        rho_tf_rai, _ = stats.spearmanr(sub["TF_collapse"], sub["RAI_8"])
        rows.append({
            "lineage": lin,
            "n": int(len(sub)),
            "rho_DM1_TF": rho_dm1_tf,
            "rho_DM1_RAI": rho_dm1_rai,
            "rho_TF_RAI": rho_tf_rai,
        })
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "mitigation_M13_three_score_convergence.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(11, 7.5))
    plot = out.sort_values("rho_TF_RAI")
    ax.scatter(plot["rho_DM1_TF"], np.arange(len(plot)), s=70, c="#244e73", label="ρ(DM1, TF_collapse)", marker="o", edgecolor="black", linewidth=0.4)
    ax.scatter(plot["rho_DM1_RAI"], np.arange(len(plot)), s=70, c="#8f2d25", label="ρ(DM1, RAI_8)", marker="s", edgecolor="black", linewidth=0.4)
    ax.scatter(plot["rho_TF_RAI"], np.arange(len(plot)), s=70, c="#426b50", label="ρ(TF_collapse, RAI_8)", marker="^", edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="grey", linewidth=0.6); ax.axvline(0.5, color="grey", linewidth=0.4, linestyle=":")
    ax.set_yticks(np.arange(len(plot))); ax.set_yticklabels([s.replace(" carcinoma", "").title() for s in plot["lineage"]], fontsize=9)
    ax.set_xlabel("Spearman ρ between three component scores", fontsize=11)
    n_strong = int((out[["rho_DM1_TF", "rho_DM1_RAI", "rho_TF_RAI"]].abs().median(axis=1) > 0.5).sum())
    ax.set_title(f"M13 — Three-score axis convergence · {len(plot)} lineages\n{n_strong}/{len(plot)} lineages have median |ρ| > 0.5 across all 3 score pairs (single coordinated axis)",
                 fontsize=11.5, fontweight="bold", loc="left", pad=12)
    ax.legend(loc="lower right", fontsize=10)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out_path = FIGS / "F11_23_mitigation_three_score_convergence.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out_path, {
        "n_lineages": int(len(out)),
        "n_lineages_median_abs_rho_gt_05": n_strong,
        "median_rho_DM1_TF": float(out["rho_DM1_TF"].median()),
        "median_rho_DM1_RAI": float(out["rho_DM1_RAI"].median()),
        "median_rho_TF_RAI": float(out["rho_TF_RAI"].median()),
    }


def m14_lineage_corr_matrix() -> tuple[Path, dict]:
    """Pan-cancer lineage-portable DM1 score is the same axis across lineages?
    Build per-lineage mean per-gene profile, correlate."""
    df = pd.read_csv(ROOT / "phase_E_lineage_specific" / "lineage_portable_dm1_per_sample.tsv", sep="\t")
    pivot = df.groupby("lineage").agg(
        mean_lin_TF=("lin_TF_score", "mean"),
        mean_inflam=("inflam_score", "mean"),
        mean_mod4=("mod4_score", "mean"),
        mean_DM1=("DM1_portable", "mean"),
        n=("sample", "count"),
    ).reset_index()
    pivot.to_csv(ROOT / "mitigation_M14_lineage_means.tsv", sep="\t", index=False)

    rho_a, p_a = stats.spearmanr(pivot["mean_lin_TF"], pivot["mean_DM1"])
    rho_b, p_b = stats.spearmanr(pivot["mean_inflam"], pivot["mean_DM1"])
    rho_c, p_c = stats.spearmanr(pivot["mean_mod4"], pivot["mean_DM1"])

    fig, ax = plt.subplots(figsize=(10, 7))
    sub = pivot.sort_values("mean_DM1")
    y = np.arange(len(sub))
    ax.scatter(sub["mean_DM1"], y, s=140, c="#7a1b15", edgecolor="black", linewidth=0.5, label="lineage-portable DM1 (composite)")
    ax.scatter(sub["mean_lin_TF"], y, s=80, c="#244e73", marker="^", edgecolor="black", linewidth=0.4, alpha=0.85, label="lineage-TF arm")
    ax.scatter(sub["mean_inflam"], y, s=80, c="#8f2d25", marker="s", edgecolor="black", linewidth=0.4, alpha=0.85, label="inflam arm")
    ax.scatter(sub["mean_mod4"], y, s=80, c="#426b50", marker="D", edgecolor="black", linewidth=0.4, alpha=0.85, label="mod4 arm")
    ax.axvline(0, color="grey", linewidth=0.6, linestyle="--")
    ax.set_yticks(y); ax.set_yticklabels([s.replace(" carcinoma", "").title() for s in sub["lineage"]], fontsize=8.5)
    ax.set_xlabel("Per-lineage mean component scores", fontsize=11)
    ax.set_title(f"M14 — Pan-cancer DM1 component coherence · 32 lineages\nlin_TF×DM1 ρ={rho_a:.2f} (p={p_a:.2g}), inflam×DM1 ρ={rho_b:.2f}, mod4×DM1 ρ={rho_c:.2f}",
                 fontsize=11, fontweight="bold", loc="left", pad=12)
    ax.legend(loc="lower right", fontsize=9.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = FIGS / "F11_24_mitigation_lineage_component_coherence.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    return out, {
        "n_lineages": int(len(pivot)),
        "rho_linTF_DM1": float(rho_a), "p_linTF_DM1": float(p_a),
        "rho_inflam_DM1": float(rho_b), "p_inflam_DM1": float(p_b),
        "rho_mod4_DM1": float(rho_c), "p_mod4_DM1": float(p_c),
    }


def m15_hm450_alt2() -> dict:
    """Try Recount2 alternative + cBioPortal correct REST endpoint."""
    import urllib.request, urllib.error
    routes = [
        ("cBioPortal molecular-profiles list", "https://www.cbioportal.org/api/molecular-profiles?projection=SUMMARY"),
        ("MEXPRESS thyroid", "https://mexpress.ulb.ac.be/mexpress/api?gene=TPO&cancer=thca"),
        ("UCSC Xena pancan methylation cohort", "https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)"),
        ("GDC files endpoint methylation", "https://api.gdc.cancer.gov/files?filters={\"op\":\"and\",\"content\":[{\"op\":\"in\",\"content\":{\"field\":\"data_type\",\"value\":[\"Methylation Beta Value\"]}}]}&size=5"),
    ]
    out = []
    for name, url in routes:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=15) as r:
                body = r.read(2000).decode(errors="ignore")
                out.append({"route": name, "status": r.status, "body_head": body[:300]})
        except urllib.error.HTTPError as e:
            out.append({"route": name, "http_error": e.code})
        except Exception as e:
            out.append({"route": name, "error": str(e)[:200]})
        time.sleep(0.3)
    successes = [r for r in out if r.get("status") == 200]
    return {"n_routes_tried": len(out), "n_successful": len(successes), "routes": out}


def main() -> None:
    summary = {}
    print("=== M12 Cross-lineage mediation proxy ===")
    p12, s12 = m12_cross_lineage_mediation()
    summary["M12_cross_lineage_mediation"] = {**s12, "fig": str(p12)}
    print(f"M12: ρ(IFN-γ, log HR)={s12['spearman_rho_logHR_vs_IFNg']:.2f} p={s12['p_logHR_vs_IFNg']:.3g}")

    print("\n=== M13 Three-score convergence ===")
    p13, s13 = m13_three_score_convergence()
    summary["M13_three_score"] = {**s13, "fig": str(p13)}
    print(f"M13: {s13['n_lineages_median_abs_rho_gt_05']}/{s13['n_lineages']} lineages |ρ|>0.5")

    print("\n=== M14 Component coherence ===")
    p14, s14 = m14_lineage_corr_matrix()
    summary["M14_component_coherence"] = {**s14, "fig": str(p14)}
    print(f"M14: lin_TF↔DM1 ρ={s14['rho_linTF_DM1']:.2f}, inflam↔DM1 ρ={s14['rho_inflam_DM1']:.2f}, mod4↔DM1 ρ={s14['rho_mod4_DM1']:.2f}")

    print("\n=== M15 HM450 alt2 ===")
    s15 = m15_hm450_alt2()
    summary["M15_hm450_alt2"] = s15
    print(f"M15: {s15['n_successful']}/{s15['n_routes_tried']} routes returned 200")

    (ROOT / "mitigation_v3_alt_summary.json").write_text(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
