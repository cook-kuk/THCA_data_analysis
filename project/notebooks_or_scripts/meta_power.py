#!/usr/bin/env python3
"""
Meta-analysis across TCGA / GSE27155 / GSE126698 DE results
plus forest & funnel plots.

Idempotent: rerunning overwrites outputs with the same content.
Does NOT modify any existing pipeline artifact.

Inputs
------
results/tables/biomarker_de_full.tsv

Outputs
-------
results/tables/meta_analysis.tsv
reports/html/figs_interactive/meta_forest_top20.html
reports/html/figs_interactive/meta_funnel.html
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT = Path("/home/seungho/personal/THCA_data_analysis/project")
DE_PATH = PROJECT / "results/tables/biomarker_de_full.tsv"
OUT_META = PROJECT / "results/tables/meta_analysis.tsv"
OUT_FOREST = PROJECT / "reports/html/figs_interactive/meta_forest_top20.html"
OUT_FUNNEL = PROJECT / "reports/html/figs_interactive/meta_funnel.html"
VENDORED_PLOTLY = "../assets/vendor/plotly/plotly-2.35.2.min.js"

# Cohort sizes (two-group setup: BRAF_like vs RAS_like-ish)
# TCGA: n_braf and n_ras per gene are in the table.
# Externals: we only have pval/log2FC, but class splits are stable per study.
# Use the pipeline-assumed totals and split 2:1 to match TCGA-class proportions.
COHORT_N = {
    "tcga":     {"n1_col": "n_braf_tcga", "n2_col": "n_ras_tcga"},
    "27155":    {"n1_fixed": 48, "n2_fixed": 24},   # 72 total, 2:1 heuristic
    "126698":   {"n1_fixed": 8,  "n2_fixed": 4},    # 12 total, 2:1 heuristic
}


def log2fc_se(log2fc: float, d: float, n1: int, n2: int) -> float | None:
    """Approximate SE of log2FC from Cohen's d and sample sizes.

    Variance of group-mean difference on the log2 scale:
        SE(diff) = sqrt( ((n1+n2)/(n1*n2)) + (d**2 / (2*(n1+n2))) ) * pooled_sd
    We approximate pooled_sd from d: d = diff / pooled_sd => pooled_sd = diff / d.
    If d ~ 0 we fall back to a conservative SE using log2fc itself.
    """
    if not (np.isfinite(log2fc) and np.isfinite(d) and n1 and n2):
        return None
    if abs(d) < 1e-6:
        # cannot invert — use a large SE proxy
        return abs(log2fc) + 1.0 if np.isfinite(log2fc) else None
    pooled_sd = log2fc / d
    if not np.isfinite(pooled_sd) or pooled_sd == 0:
        return None
    base = (n1 + n2) / (n1 * n2) + (d ** 2) / (2 * (n1 + n2))
    se = math.sqrt(max(base, 1e-12)) * abs(pooled_sd)
    return se if np.isfinite(se) and se > 0 else None


def p_to_z(p: float, effect_sign: float) -> float | None:
    if p is None or not np.isfinite(p) or p <= 0 or p >= 1:
        return None
    z = stats.norm.isf(p / 2.0)
    return z * (1.0 if effect_sign >= 0 else -1.0)


def cohort_se_from_p(log2fc: float, p: float) -> float | None:
    """If Cohen's d missing (external cohorts), derive SE from two-sided p.
       log2fc / SE = z -> SE = log2fc / z.
    """
    if not (np.isfinite(log2fc) and np.isfinite(p)) or p <= 0 or p >= 1:
        return None
    z = stats.norm.isf(p / 2.0)
    if z <= 0 or not np.isfinite(z):
        return None
    if log2fc == 0:
        return 1.0
    return abs(log2fc) / z


def fixed_effect(effs: np.ndarray, ses: np.ndarray) -> tuple[float, float, float, float]:
    """Return (mean, se, z, p) under fixed-effect inverse-variance."""
    w = 1.0 / (ses ** 2)
    mu = float((w * effs).sum() / w.sum())
    se = float(math.sqrt(1.0 / w.sum()))
    z = mu / se if se else 0.0
    p = 2.0 * stats.norm.sf(abs(z))
    return mu, se, z, p


def heterogeneity(effs: np.ndarray, ses: np.ndarray, fe_mu: float) -> tuple[float, float, float, float]:
    """Cochran's Q, df, I² (percent), tau² via DerSimonian-Laird."""
    w = 1.0 / (ses ** 2)
    Q = float((w * (effs - fe_mu) ** 2).sum())
    df = len(effs) - 1
    if df <= 0:
        return Q, float(df), 0.0, 0.0
    C = float(w.sum() - (w ** 2).sum() / w.sum())
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    I2 = max(0.0, (Q - df) / Q) * 100.0 if Q > 0 else 0.0
    return Q, float(df), float(I2), float(tau2)


def random_effect(effs: np.ndarray, ses: np.ndarray, tau2: float) -> tuple[float, float, float, float]:
    w = 1.0 / (ses ** 2 + tau2)
    mu = float((w * effs).sum() / w.sum())
    se = float(math.sqrt(1.0 / w.sum()))
    z = mu / se if se else 0.0
    p = 2.0 * stats.norm.sf(abs(z))
    return mu, se, z, p


def fishers_combine(pvals: list[float]) -> float:
    clean = [p for p in pvals if p is not None and np.isfinite(p) and 0 < p <= 1]
    if len(clean) < 2:
        return float("nan")
    X = -2.0 * sum(math.log(p) for p in clean)
    df = 2 * len(clean)
    return float(stats.chi2.sf(X, df))


# ---------------------------------------------------------------- main
def run_meta() -> pd.DataFrame:
    df = pd.read_csv(DE_PATH, sep="\t")
    records = []
    for _, row in df.iterrows():
        cohort_effs: list[tuple[str, float, float, float]] = []  # (cohort, eff, se, p)
        # TCGA
        if row.get("present_tcga", False) and np.isfinite(row.get("log2FC_tcga", np.nan)):
            se_t = log2fc_se(row["log2FC_tcga"], row.get("cohens_d_tcga", 0.0),
                             int(row.get("n_braf_tcga") or 0), int(row.get("n_ras_tcga") or 0))
            if se_t:
                cohort_effs.append(("TCGA", float(row["log2FC_tcga"]), float(se_t), float(row["pval_tcga"])))
        # GSE27155
        if row.get("present_27155", False) and np.isfinite(row.get("log2FC_27155", np.nan)):
            se_g = cohort_se_from_p(row["log2FC_27155"], row.get("pval_27155", np.nan))
            if se_g:
                cohort_effs.append(("GSE27155", float(row["log2FC_27155"]), float(se_g), float(row["pval_27155"])))
        # GSE126698
        if row.get("present_126698", False) and np.isfinite(row.get("log2FC_126698", np.nan)):
            se_g = cohort_se_from_p(row["log2FC_126698"], row.get("pval_126698", np.nan))
            if se_g:
                cohort_effs.append(("GSE126698", float(row["log2FC_126698"]), float(se_g), float(row["pval_126698"])))

        if len(cohort_effs) < 2:
            continue
        effs = np.array([c[1] for c in cohort_effs], dtype=float)
        ses = np.array([c[2] for c in cohort_effs], dtype=float)
        pvs = [c[3] for c in cohort_effs]
        fe_mu, fe_se, fe_z, fe_p = fixed_effect(effs, ses)
        Q, df_Q, I2, tau2 = heterogeneity(effs, ses, fe_mu)
        re_mu, re_se, re_z, re_p = random_effect(effs, ses, tau2)
        fisher_p = fishers_combine(pvs)
        signs = np.sign(effs)
        concord = bool(len(set(signs.tolist())) == 1 and 0 not in signs)

        records.append({
            "gene": row["gene"],
            "n_cohorts": len(cohort_effs),
            "cohorts": ",".join(c[0] for c in cohort_effs),
            "log2fc_tcga": row.get("log2FC_tcga"),
            "log2fc_27155": row.get("log2FC_27155"),
            "log2fc_126698": row.get("log2FC_126698"),
            "fe_log2fc": fe_mu, "fe_se": fe_se, "fe_z": fe_z, "fe_p": fe_p,
            "meta_log2fc": re_mu, "meta_se": re_se, "meta_z": re_z, "meta_p": re_p,
            "fisher_p": fisher_p,
            "Q": Q, "Q_df": df_Q, "I2": I2, "tau2": tau2,
            "direction_concordance": concord,
            "is_known": bool(row.get("is_known", False)),
            "is_novel_validated": bool(row.get("is_novel_validated", False)),
            "cohort_effs_json": str([(c[0], c[1], c[2]) for c in cohort_effs]),
        })

    meta = pd.DataFrame.from_records(records)
    if len(meta):
        # FDR (BH) across genes on random-effect p-value
        pvals = meta["meta_p"].values
        order = np.argsort(pvals)
        ranked = pvals[order]
        m = len(ranked)
        q = ranked * m / (np.arange(m) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q_out = np.empty(m); q_out[order] = np.clip(q, 0, 1)
        meta["meta_q"] = q_out
    OUT_META.parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(OUT_META, sep="\t", index=False)
    print(f"[meta] wrote {OUT_META} rows={len(meta)}")
    return meta


# ---------------------------------------------------------------- plots
def _html_shell(title: str, fig_json: str, extra_js: str = "") -> str:
    return f"""<!doctype html>
<html lang="ko"><head><meta charset="utf-8"><title>{title}</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>html,body{{margin:0;background:#07132b;color:#e2e8f0;font-family:Inter,system-ui,sans-serif}}#fig{{width:100%;height:100vh}}</style>
</head><body>
<div id="fig"></div>
<script src="{VENDORED_PLOTLY}"></script>
<script>
const SPEC = {fig_json};
Plotly.newPlot('fig', SPEC.data, SPEC.layout, {{responsive:true, displaylogo:false}});
{extra_js}
</script>
</body></html>
"""


def make_forest(meta: pd.DataFrame) -> None:
    import json
    # Rank by meta_q (lowest) among genes with ≥2 cohorts; require direction concordance first
    top = meta.sort_values(["direction_concordance", "meta_q"], ascending=[False, True]).head(20).copy()
    top = top.iloc[::-1]  # invert so best gene plots at top

    # Consistent per-gene cohort spacing
    traces = []
    colors = {"TCGA": "#5eead4", "GSE27155": "#fca5a5", "GSE126698": "#fcd34d"}
    y_labels = []
    y_positions = []
    diamond_y = []
    diamond_x = []
    diamond_lo = []
    diamond_hi = []
    diamond_text = []

    y_cursor = 0
    yticks = []
    yticktext = []
    for _, r in top.iterrows():
        # parse cohort_effs_json (stored as Python literal)
        try:
            cohort_effs = eval(r["cohort_effs_json"])  # noqa: S307 - trusted self-made
        except Exception:
            continue
        label = f"{r['gene']}"
        yticks.append(y_cursor + (len(cohort_effs) + 1) / 2)
        yticktext.append(label)
        for (cname, eff, se) in cohort_effs:
            lo, hi = eff - 1.96 * se, eff + 1.96 * se
            traces.append({
                "type": "scatter", "mode": "markers",
                "x": [eff], "y": [y_cursor],
                "error_x": {"type": "data", "array": [1.96 * se], "thickness": 1.2, "width": 4, "color": colors.get(cname, "#94a3b8")},
                "marker": {"color": colors.get(cname, "#94a3b8"), "size": 7},
                "name": cname, "legendgroup": cname, "showlegend": False,
                "hovertemplate": f"<b>{r['gene']}</b><br>{cname}: log2FC={eff:.3f} (95%CI {lo:.2f}, {hi:.2f})<extra></extra>",
            })
            y_cursor += 1
        # diamond for summary (random-effect)
        diamond_y.append(y_cursor)
        diamond_x.append(r["meta_log2fc"])
        diamond_lo.append(r["meta_log2fc"] - 1.96 * r["meta_se"])
        diamond_hi.append(r["meta_log2fc"] + 1.96 * r["meta_se"])
        diamond_text.append(f"<b>{r['gene']}</b> meta log2FC={r['meta_log2fc']:.3f} q={r['meta_q']:.2e} I²={r['I2']:.0f}%")
        y_cursor += 2

    # summary diamonds as single trace
    traces.append({
        "type": "scatter", "mode": "markers",
        "x": diamond_x, "y": diamond_y,
        "error_x": {"type": "data", "symmetric": False,
                    "array": [hi - x for hi, x in zip(diamond_hi, diamond_x)],
                    "arrayminus": [x - lo for x, lo in zip(diamond_x, diamond_lo)],
                    "thickness": 2.4, "width": 8, "color": "#c084fc"},
        "marker": {"color": "#c084fc", "size": 12, "symbol": "diamond"},
        "name": "Random-effect summary",
        "hovertext": diamond_text, "hoverinfo": "text",
    })
    # Legend proxies for the 4 categories
    for cname, col in colors.items():
        traces.append({
            "type": "scatter", "mode": "markers",
            "x": [None], "y": [None],
            "marker": {"color": col, "size": 7},
            "name": cname, "showlegend": True,
        })
    traces.append({
        "type": "scatter", "mode": "markers",
        "x": [None], "y": [None],
        "marker": {"color": "#c084fc", "size": 12, "symbol": "diamond"},
        "name": "Meta (RE)",
    })

    layout = {
        "title": {"text": "메타분석 Forest Plot · 상위 20개 유전자 (direction-concordant)", "x": 0.01, "font": {"color": "#e2e8f0"}},
        "paper_bgcolor": "#07132b", "plot_bgcolor": "#0b1a3a",
        "font": {"color": "#cbd5e1"},
        "xaxis": {"title": "log2 fold change (BRAF_like vs RAS_like)", "zeroline": True, "zerolinecolor": "#64748b", "gridcolor": "rgba(148,163,184,.18)"},
        "yaxis": {"tickvals": yticks, "ticktext": yticktext, "autorange": "reversed", "gridcolor": "rgba(148,163,184,.08)"},
        "height": max(600, 32 * len(top) + 120),
        "margin": {"l": 90, "r": 30, "t": 60, "b": 60},
        "legend": {"orientation": "h", "x": 0, "y": 1.05, "bgcolor": "rgba(0,0,0,0)"},
        "hoverlabel": {"bgcolor": "#0b1a3a", "bordercolor": "#5eead4"},
    }
    spec = {"data": traces, "layout": layout}
    OUT_FOREST.parent.mkdir(parents=True, exist_ok=True)
    OUT_FOREST.write_text(_html_shell("Meta Forest (Top 20)", json.dumps(spec)))
    print(f"[meta] forest → {OUT_FOREST}")


def make_funnel(meta: pd.DataFrame) -> None:
    import json
    # funnel: x=meta_log2fc (random-effect), y=1/SE (inverse) or SE (flipped)
    m = meta.dropna(subset=["meta_log2fc", "meta_se"]).copy()
    if not len(m):
        return
    # Sub-sample for browser stability
    if len(m) > 5000:
        m = m.sample(5000, random_state=42)
    colors = np.where(m["direction_concordance"], "#5eead4", "#f87171")
    trace = {
        "type": "scattergl", "mode": "markers",
        "x": m["meta_log2fc"].tolist(), "y": m["meta_se"].tolist(),
        "marker": {"color": colors.tolist(), "size": 4, "opacity": 0.55,
                   "line": {"width": 0}},
        "text": m["gene"].tolist(),
        "hovertemplate": "<b>%{text}</b><br>log2FC=%{x:.2f}<br>SE=%{y:.3f}<extra></extra>",
        "name": "genes",
    }
    # Pseudo-confidence funnel lines (95%)
    max_se = float(m["meta_se"].quantile(0.99))
    ses = np.linspace(0.001, max_se, 50).tolist()
    xs_hi = [1.96 * s for s in ses]
    xs_lo = [-1.96 * s for s in ses]
    layout = {
        "title": {"text": "Funnel Plot · 출판 편향 점검 (random-effect)", "x": 0.01, "font": {"color": "#e2e8f0"}},
        "paper_bgcolor": "#07132b", "plot_bgcolor": "#0b1a3a",
        "font": {"color": "#cbd5e1"},
        "xaxis": {"title": "meta log2 fold change", "gridcolor": "rgba(148,163,184,.18)"},
        "yaxis": {"title": "Standard Error", "autorange": "reversed", "gridcolor": "rgba(148,163,184,.18)"},
        "height": 640,
        "margin": {"l": 80, "r": 30, "t": 60, "b": 60},
        "shapes": [],
    }
    funnel_hi = {"type": "scatter", "mode": "lines", "x": xs_hi, "y": ses,
                 "line": {"color": "rgba(148,163,184,.6)", "dash": "dash", "width": 1},
                 "name": "95% pseudo-CI"}
    funnel_lo = {"type": "scatter", "mode": "lines", "x": xs_lo, "y": ses,
                 "line": {"color": "rgba(148,163,184,.6)", "dash": "dash", "width": 1},
                 "name": "95% pseudo-CI", "showlegend": False}
    spec = {"data": [trace, funnel_hi, funnel_lo], "layout": layout}
    OUT_FUNNEL.parent.mkdir(parents=True, exist_ok=True)
    OUT_FUNNEL.write_text(_html_shell("Meta Funnel", json.dumps(spec)))
    print(f"[meta] funnel → {OUT_FUNNEL}")


# ---------------------------------------------------------------- CLI
if __name__ == "__main__":
    print("=== THCA meta-analysis + power bundle ===")
    meta = run_meta()
    if len(meta):
        make_forest(meta)
        make_funnel(meta)
    # summary stats → stdout
    print("\n[top-10 most-replicated, concordant direction]")
    sel = meta[meta["direction_concordance"]].sort_values("meta_q").head(10)
    for _, r in sel.iterrows():
        print(f"  {r['gene']:<12s}  meta log2FC={r['meta_log2fc']:+.3f}  q={r['meta_q']:.2e}  I²={r['I2']:.0f}%  cohorts={r['n_cohorts']}")
    print("\nDone.")
