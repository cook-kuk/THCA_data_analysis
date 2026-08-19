#!/usr/bin/env python3
"""Create reviewer-safe CLEAN-NeoBench visual dashboard assets."""

from __future__ import annotations

import argparse
import html
import json
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import ensure_dir, update_manifest


PAGE_NAME = "clean_neobench_visual_dashboard_2026_05_10.html"
ASSET_DIR = "clean_neobench_barneo_visuals"
BG = "#0d1117"
PANEL = "#101820"
INK = "#e6edf3"
MUTED = "#9aa7b4"
LINE = "#2a3441"
GOLD = "#e3b341"
CYAN = "#5eead4"
GREEN = "#86efac"
RED = "#f87171"
VIOLET = "#b58cff"
BLUE = "#60a5fa"
ROLE_COLORS = {
    "anchor": GOLD,
    "internal_candidate": CYAN,
    "bounded_fallback": VIOLET,
    "caveated_public_comparator": RED,
    "uncertainty_only": BLUE,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    p.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub root")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def nnum(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def safe_int(value: Any) -> int:
    try:
        return int(value)
    except Exception:
        return 0


def prep_ax(ax: plt.Axes, title: str, subtitle: str = "") -> None:
    ax.set_facecolor(PANEL)
    ax.grid(axis="x", color=LINE, alpha=0.45, linewidth=0.8)
    ax.tick_params(colors=MUTED, labelsize=8)
    for spine in ax.spines.values():
        spine.set_color(LINE)
    ax.set_title("")
    ax.text(0, 1.085, title, transform=ax.transAxes, color=INK, fontsize=13, fontweight="bold", va="bottom", clip_on=False)
    if subtitle:
        ax.text(0, 1.045, subtitle, transform=ax.transAxes, color=MUTED, fontsize=8, va="bottom", clip_on=False)


def save(fig: plt.Figure, fig_dir: Path, asset_dir: Path, name: str) -> str:
    fig.patch.set_facecolor(BG)
    out = fig_dir / name
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    shutil.copy2(out, asset_dir / name)
    return name


def fig_leaderboard(lb: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(11, 7))
    if lb.empty:
        ax.text(0.5, 0.5, "No leaderboard", color=INK, ha="center")
    else:
        df = lb.head(18).iloc[::-1].copy()
        colors = [ROLE_COLORS.get(str(x), MUTED) for x in df["method_role"]]
        ax.barh(df["method_name"], nnum(df, "mean_AUPRC"), color=colors, alpha=0.9)
        ax.set_xlabel("Mean AUPRC", color=MUTED)
    prep_ax(ax, "Leaderboard by AUPRC", "Role colors expose anchor/internal/fallback/public caveats.")
    return save(fig, fig_dir, asset_dir, "fig01_leaderboard_auprc.png")


def fig_source(master: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(11, 5.5))
    if master.empty:
        ax.text(0.5, 0.5, "No master table", color=INK, ha="center")
    else:
        g = master.groupby("source_name").agg(n=("candidate_id", "size"), prev=("label", "mean")).reset_index()
        g = g.sort_values("n", ascending=False).head(12)
        x = np.arange(len(g))
        ax.bar(x, g["n"], color=CYAN, alpha=0.82)
        ax.set_xticks(x)
        ax.set_xticklabels(g["source_name"].astype(str), rotation=30, ha="right")
        ax.set_ylabel("Candidates", color=MUTED)
        ax2 = ax.twinx()
        ax2.plot(x, g["prev"], color=GOLD, marker="o", linewidth=2.2)
        ax2.set_ylim(0, 1.05)
        ax2.set_ylabel("Positive prevalence", color=GOLD)
        ax2.tick_params(colors=GOLD, labelsize=8)
        for spine in ax2.spines.values():
            spine.set_color(LINE)
    prep_ax(ax, "Source prevalence shift", "Positive prevalence differs strongly across CEDAR/TESLA/ITSNdb sources.")
    return save(fig, fig_dir, asset_dir, "fig02_source_prevalence_shift.png")


def fig_funnel(master: pd.DataFrame, barneo: pd.DataFrame, bma: pd.DataFrame, ctx: pd.DataFrame, fail: pd.DataFrame, xscore: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    vals = [
        len(master),
        int((~barneo["abstain"].astype(bool)).sum()) if "abstain" in barneo else 0,
        int((~bma["bma_abstain"].astype(bool)).sum()) if "bma_abstain" in bma else 0,
        int(ctx["contextual_clean_claim_allowed"].astype(bool).sum()) if "contextual_clean_claim_allowed" in ctx else 0,
        int((~fail["failure_aware_abstain"].astype(bool)).sum()) if "failure_aware_abstain" in fail else 0,
        int((xscore.get("barneo_x_primary_action", pd.Series(dtype=str)) == "priority_review_candidate").sum()) if not xscore.empty else 0,
    ]
    labels = ["all", "BAR-Neo\nnon-abstain", "BMA\nnon-abstain", "contextual\nclean", "failure-aware\nnon-abstain", "X priority"]
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    bars = ax.bar(np.arange(len(vals)), vals, color=[CYAN, BLUE, VIOLET, GREEN, GOLD, RED], alpha=0.9)
    ax.set_xticks(np.arange(len(vals)))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Rows", color=MUTED)
    ymax = max(vals + [1]) * 1.18
    ax.set_ylim(0, ymax)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + ymax * 0.025, f"{v:,}", color=INK, ha="center", fontsize=9)
    prep_ax(ax, "Abstention funnel", "Discovery candidates narrow sharply to reviewer-safe prompts.")
    return save(fig, fig_dir, asset_dir, "fig03_abstention_funnel.png")


def fig_vulnerability(vuln: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(11, 6.2))
    cols = [
        "source_score_prevalence_corr",
        "low_prevalence_high_ranked_negative_rate",
        "rare_hla_missed_positive_rate",
        "external_missed_positive_rate",
        "distribution_vulnerability_score",
    ]
    if vuln.empty:
        ax.text(0.5, 0.5, "No vulnerability audit", color=INK, ha="center")
    else:
        df = vuln.head(14).copy()
        use = [c for c in cols if c in df]
        data = df[use].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy()
        im = ax.imshow(data, aspect="auto", cmap="magma", vmin=0)
        ax.set_yticks(np.arange(len(df)))
        ax.set_yticklabels(df["method_name"].astype(str))
        ax.set_xticks(np.arange(len(use)))
        ax.set_xticklabels([x.replace("_", "\n") for x in use])
        cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        cb.ax.tick_params(colors=MUTED, labelsize=8)
        cb.outline.set_edgecolor(LINE)
    prep_ax(ax, "Distribution vulnerability heatmap", "Higher cells indicate source/HLA/external fragility.")
    return save(fig, fig_dir, asset_dir, "fig04_method_vulnerability_heatmap.png")


def fig_challenges(summary: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(11, 6.2))
    if summary.empty:
        ax.text(0.5, 0.5, "No challenge pack", color=INK, ha="center")
    else:
        df = summary.sort_values("n_unique_candidates", ascending=True)
        y = np.arange(len(df))
        ax.barh(y, df["n_unique_candidates"], color=CYAN, alpha=0.86)
        ax.set_yticks(y)
        ax.set_yticklabels(df["challenge_axis"].astype(str))
        ax.set_xlabel("Unique candidates", color=MUTED)
        prev = pd.to_numeric(df["positive_prevalence"], errors="coerce").fillna(0)
        for i, (n, p) in enumerate(zip(df["n_unique_candidates"], prev)):
            ax.text(float(n) + 2, i, f"prev={p:.2f}", color=GOLD, va="center", fontsize=8)
    prep_ax(ax, "Challenge pack axes", "Stress tests and manual-review queues, not clinical selections.")
    return save(fig, fig_dir, asset_dir, "fig05_challenge_axes.png")


def fig_context(ctx: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(8.4, 7.0))
    if ctx.empty:
        ax.text(0.5, 0.5, "No contextual BMA", color=INK, ha="center")
    else:
        label = pd.to_numeric(ctx.get("label", pd.Series(np.nan, index=ctx.index)), errors="coerce")
        color = np.where(label.eq(1), GREEN, RED)
        high_leak = ctx.get("leakage_risk_level", pd.Series("", index=ctx.index)).astype(str).str.lower().eq("high")
        ax.scatter(nnum(ctx, "contextual_bma_score"), nnum(ctx, "contextual_confidence_score"), c=color, s=np.where(high_leak, 18, 10), alpha=0.36, edgecolors="none")
        ax.axhline(0.45, color=GOLD, linestyle="--", linewidth=1.1)
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlabel("Contextual BMA score", color=MUTED)
        ax.set_ylabel("Contextual confidence", color=MUTED)
    prep_ax(ax, "Contextual score vs confidence", "Green=positive, red=negative; larger points are high leakage risk.")
    return save(fig, fig_dir, asset_dir, "fig06_contextual_score_confidence.png")


def fig_patient(req: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    if req.empty:
        ax.text(0.5, 0.5, "No patient-gate requirements", color=INK, ha="center")
    else:
        p = req.groupby(["gate_group", "metadata_status"]).size().reset_index(name="n")
        pivot = p.pivot(index="gate_group", columns="metadata_status", values="n").fillna(0)
        bottom = np.zeros(len(pivot))
        x = np.arange(len(pivot))
        for col, color in zip(pivot.columns, [RED, GOLD, GREEN, CYAN, VIOLET, BLUE]):
            ax.bar(x, pivot[col], bottom=bottom, label=str(col), color=color, alpha=0.88)
            bottom += pivot[col].to_numpy()
        ax.set_xticks(x)
        ax.set_xticklabels(pivot.index.astype(str), rotation=25, ha="right")
        ax.legend(frameon=False, fontsize=8, labelcolor=INK)
        ax.set_ylabel("Required fields", color=MUTED)
    prep_ax(ax, "PAAD/THCA metadata blockers", "Missing patient context keeps patient-gated output demo-only.")
    return save(fig, fig_dir, asset_dir, "fig07_patient_metadata_blockers.png")


def fig_public(audit: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    if audit.empty:
        ax.text(0.5, 0.5, "No public audit", color=INK, ha="center")
    else:
        df = audit.copy()
        if "uses_public_pretraining" in df:
            df = df[df["uses_public_pretraining"].astype(str).str.lower().isin(["true", "1"])]
        allowed = df.get("clean_comparator_allowed_after_audit", pd.Series(False, index=df.index)).astype(str).str.lower().isin(["true", "1"])
        counts = pd.Series(np.where(allowed, "clean allowed", "caveated/unresolved")).value_counts()
        ax.pie(counts.values, labels=counts.index, autopct="%1.0f%%", colors=[GREEN, RED], textprops={"color": INK, "fontsize": 9})
    ax.set_title("Public pretrained comparator status", color=INK, fontsize=13, fontweight="bold")
    return save(fig, fig_dir, asset_dir, "fig08_public_tool_caveat.png")


def fig_split_heatmap(metrics: pd.DataFrame, leaderboard: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(12, 8))
    if metrics.empty or leaderboard.empty:
        ax.text(0.5, 0.5, "No split metrics", color=INK, ha="center")
    else:
        focus = leaderboard.head(16)["method_name"].astype(str).tolist()
        sub = metrics[metrics["method_name"].isin(focus)].copy()
        grouped = (
            sub.groupby(["method_name", "split_contract"])["AUPRC"]
            .median()
            .reset_index()
        )
        pivot = grouped.pivot(index="method_name", columns="split_contract", values="AUPRC").reindex(focus)
        pivot = pivot.dropna(axis=1, how="all")
        data = pivot.apply(pd.to_numeric, errors="coerce").fillna(np.nan).to_numpy()
        im = ax.imshow(data, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        ax.set_yticks(np.arange(len(pivot.index)))
        ax.set_yticklabels(pivot.index.astype(str))
        ax.set_xticks(np.arange(len(pivot.columns)))
        ax.set_xticklabels([c.replace("_", "\n") for c in pivot.columns], rotation=0)
        cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        cb.ax.tick_params(colors=MUTED, labelsize=8)
        cb.outline.set_edgecolor(LINE)
    prep_ax(ax, "Split-contract AUPRC heatmap", "Top methods by leaderboard; missing/single-class splits stay blank.")
    return save(fig, fig_dir, asset_dir, "fig09_split_contract_auprc_heatmap.png")


def fig_delta_anchor(metrics: pd.DataFrame, leaderboard: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(11, 7))
    if metrics.empty or "Structure_LR" not in set(metrics.get("method_name", [])):
        ax.text(0.5, 0.5, "Structure_LR anchor metrics unavailable", color=INK, ha="center")
    else:
        focus = [m for m in leaderboard.head(14)["method_name"].astype(str).tolist() if m != "Structure_LR"]
        sub = metrics[metrics["method_name"].isin(focus + ["Structure_LR"])].copy()
        key_cols = ["split_contract", "split_group"]
        anchor = sub[sub["method_name"].eq("Structure_LR")][key_cols + ["AUPRC"]].rename(columns={"AUPRC": "anchor_AUPRC"})
        merged = sub[~sub["method_name"].eq("Structure_LR")].merge(anchor, on=key_cols, how="inner")
        merged["delta_AUPRC"] = pd.to_numeric(merged["AUPRC"], errors="coerce") - pd.to_numeric(merged["anchor_AUPRC"], errors="coerce")
        agg = merged.groupby("method_name")["delta_AUPRC"].median().reindex(focus).dropna().sort_values()
        colors = np.where(agg >= 0, GREEN, RED)
        ax.barh(agg.index, agg.values, color=colors, alpha=0.9)
        ax.axvline(0, color=GOLD, linewidth=1.2)
        ax.set_xlabel("Median AUPRC delta vs Structure_LR", color=MUTED)
    prep_ax(ax, "Win/loss versus Structure_LR anchor", "Positive bars beat the honest local anchor across matched split rows.")
    return save(fig, fig_dir, asset_dir, "fig10_delta_vs_structure_lr.png")


def fig_context_weight_heatmap(weights: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(12, 7.5))
    if weights.empty:
        ax.text(0.5, 0.5, "No contextual weights", color=INK, ha="center")
    else:
        contexts = ["global", "low_prevalence", "rare_hla", "external_or_holdout", "high_leakage_review_only", "korean_hla_focus"]
        top_methods = (
            weights[weights["context_label"].isin(contexts)]
            .groupby("method_name")["contextual_weight"]
            .max()
            .sort_values(ascending=False)
            .head(18)
            .index
        )
        sub = weights[weights["context_label"].isin(contexts) & weights["method_name"].isin(top_methods)]
        pivot = sub.pivot_table(index="method_name", columns="context_label", values="contextual_weight", aggfunc="max").reindex(top_methods)
        pivot = pivot.reindex(columns=[c for c in contexts if c in pivot.columns])
        data = pivot.fillna(0).to_numpy()
        im = ax.imshow(data, aspect="auto", cmap="cividis")
        ax.set_yticks(np.arange(len(pivot.index)))
        ax.set_yticklabels(pivot.index.astype(str))
        ax.set_xticks(np.arange(len(pivot.columns)))
        ax.set_xticklabels([c.replace("_", "\n") for c in pivot.columns])
        cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        cb.ax.tick_params(colors=MUTED, labelsize=8)
        cb.outline.set_edgecolor(LINE)
    prep_ax(ax, "Contextual method posterior weights", "Weights change by source/HLA/stress context instead of one global ensemble.")
    return save(fig, fig_dir, asset_dir, "fig11_contextual_method_weights.png")


def fig_claim_safe_top(explain: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(11, 6.5))
    if explain.empty:
        ax.text(0.5, 0.5, "No BAR-Neo-X explanations", color=INK, ha="center")
    else:
        df = explain.sort_values("barneo_x_claim_safe_rank").head(15).iloc[::-1].copy()
        y = np.arange(len(df))
        ax.barh(y, pd.to_numeric(df["barneo_x_discovery_score"], errors="coerce"), color=BLUE, alpha=0.45, label="discovery")
        ax.barh(y, pd.to_numeric(df["barneo_x_claim_safe_score"], errors="coerce"), color=GREEN, alpha=0.88, label="claim-safe")
        labels = df["candidate_id"].astype(str) + " · " + df["hla_allele_4digit"].astype(str)
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_xlim(0, 1.0)
        ax.legend(frameon=False, fontsize=8, labelcolor=INK)
        ax.set_xlabel("Score", color=MUTED)
    prep_ax(ax, "BAR-Neo-X claim-safe top rows", "Claim-safe score is deliberately lower when metadata/leakage penalties apply.")
    return save(fig, fig_dir, asset_dir, "fig12_barneo_x_claim_safe_top.png")


def fig_stress_actions(stress_scores: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    if stress_scores.empty or "stress_guarded_action" not in stress_scores:
        ax.text(0.5, 0.5, "No stress-guarded scores", color=INK, ha="center")
    else:
        counts = stress_scores["stress_guarded_action"].value_counts().sort_values()
        colors = [GREEN if "claim_safe" in x else GOLD if "priority" in x else CYAN if "support" in x else RED if "abstain" in x else BLUE for x in counts.index]
        ax.barh(counts.index.astype(str), counts.values, color=colors, alpha=0.9)
        ax.set_xlabel("Candidates", color=MUTED)
        for i, v in enumerate(counts.values):
            ax.text(v + max(counts.values) * 0.015, i, f"{int(v):,}", color=INK, va="center", fontsize=9)
    prep_ax(ax, "Stress-guarded BAR-Neo actions", "Action bins after source/HLA/leakage-aware method downweighting.")
    return save(fig, fig_dir, asset_dir, "fig13_stress_guarded_actions.png")


def fig_stress_top(stress_scores: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(11, 6.5))
    if stress_scores.empty:
        ax.text(0.5, 0.5, "No stress-guarded scores", color=INK, ha="center")
    else:
        df = stress_scores.sort_values("stress_guarded_rank_global").head(18).iloc[::-1].copy()
        labels = df["candidate_id"].astype(str) + " · " + df["hla_allele_4digit"].astype(str)
        y = np.arange(len(df))
        ax.barh(y, nnum(df, "stress_guarded_discovery_score"), color=BLUE, alpha=0.40, label="discovery")
        ax.barh(y, nnum(df, "stress_guarded_claim_safe_score"), color=GREEN, alpha=0.82, label="claim-safe")
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_xlim(0, 1.0)
        ax.legend(frameon=False, fontsize=8, labelcolor=INK)
        ax.set_xlabel("Stress-guarded score", color=MUTED)
    prep_ax(ax, "Stress-guarded top candidates", "Top rows are manual-audit research triage candidates, not clinical selections.")
    return save(fig, fig_dir, asset_dir, "fig14_stress_guarded_top_candidates.png")


def fig_stress_method_weights(method_weights: pd.DataFrame, fig_dir: Path, asset_dir: Path) -> str:
    fig, ax = plt.subplots(figsize=(11, 7))
    if method_weights.empty:
        ax.text(0.5, 0.5, "No stress-guarded method weights", color=INK, ha="center")
    else:
        df = method_weights.sort_values("stress_guarded_base_weight", ascending=False).head(18).iloc[::-1].copy()
        colors = [ROLE_COLORS.get(str(x), MUTED) for x in df["method_role"]]
        ax.barh(df["method_name"], nnum(df, "stress_guarded_base_weight"), color=colors, alpha=0.9)
        ax.set_xlabel("Stress-guarded base weight", color=MUTED)
    prep_ax(ax, "Stress-guarded method weights", "Weights reward methods that survive source/HLA stress and cap public/fallback support.")
    return save(fig, fig_dir, asset_dir, "fig15_stress_guarded_method_weights.png")


def write_html(hub_root: Path, output_root: Path, figures: list[str], stats: dict[str, int]) -> None:
    cards = "\n".join(f'<section><img src="assets/{ASSET_DIR}/{html.escape(f)}"><p>{html.escape(f)}</p></section>' for f in figures)
    stat_html = "".join(f"<div><b>{v:,}</b><span>{html.escape(k)}</span></div>" for k, v in stats.items())
    text = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CLEAN-NeoBench Visual Dashboard</title>
<style>
body{{margin:0;background:{BG};color:{INK};font-family:JetBrains Mono,ui-monospace,Menlo,monospace}} header{{padding:42px 34px 24px;border-bottom:1px solid {LINE}}}
h1{{font-family:Georgia,serif;font-size:48px;margin:0 0 12px}} .lead{{color:#c8d1dc;max-width:980px}}
.stats{{display:grid;grid-template-columns:repeat(6,minmax(120px,1fr));gap:10px;margin-top:22px}} .stats div{{border:1px solid {LINE};background:{PANEL};padding:13px}}
.stats b{{display:block;color:{CYAN};font-size:24px}} .stats span{{display:block;color:{MUTED};font-size:12px}}
main{{max-width:1320px;margin:auto;padding:26px;display:grid;gap:22px}} section{{border:1px solid {LINE};background:{PANEL};padding:14px}} img{{width:100%;height:auto;display:block}} p{{color:{MUTED};font-size:12px}}
@media(max-width:900px){{.stats{{grid-template-columns:1fr 1fr}}h1{{font-size:34px}}main{{padding:14px}}}}
</style></head><body><header><h1>CLEAN-NeoBench Visual Dashboard</h1>
<p class="lead">Reviewer-safe visuals for leakage-aware neoantigen benchmarking, contextual BAR-Neo-BMA, failure modes, challenge-pack design, public-tool caveats, and PAAD/THCA metadata blockers. No SOTA, clinical selection, or quantum-advantage claim is made.</p>
<div class="stats">{stat_html}</div></header><main>{cards}</main>
<p style="padding:0 26px 40px;color:{MUTED}">Output root: {html.escape(str(output_root))}</p></body></html>"""
    (hub_root / PAGE_NAME).write_text(text)


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    hub_root = Path(args.hub_root)
    fig_dir = output_root / "figures"
    asset_dir = hub_root / "assets" / ASSET_DIR
    ensure_dir(fig_dir)
    ensure_dir(asset_dir)
    ensure_dir(hub_root)

    lb = read_tsv(output_root / "clean_neobench_leaderboard.tsv")
    master = read_tsv(output_root / "clean_neobench_master.tsv")
    barneo = read_tsv(output_root / "barneo_candidate_scores.tsv")
    bma = read_tsv(output_root / "barneo_bma_candidate_scores.tsv")
    ctx = read_tsv(output_root / "barneo_contextual_bma_candidate_scores.tsv")
    fail = read_tsv(output_root / "barneo_failure_aware_candidate_scores.tsv")
    xscore = read_tsv(output_root / "barneo_x_candidate_scores.tsv")
    vuln = read_tsv(output_root / "clean_neobench_method_distribution_vulnerability.tsv")
    challenge = read_tsv(output_root / "clean_neobench_challenge_axis_summary.tsv")
    patient = read_tsv(output_root / "patient_gated_clean_neo_metadata_requirements.tsv")
    public = read_tsv(output_root / "clean_neobench_public_tool_overlap_audit.tsv")
    metrics = read_tsv(output_root / "clean_neobench_split_metrics.tsv")
    ctx_weights = read_tsv(output_root / "barneo_contextual_bma_method_weights.tsv")
    x_explain = read_tsv(output_root / "barneo_x_candidate_explanations.tsv")
    stress_scores = read_tsv(output_root / "barneo_stress_guarded_candidate_scores.tsv")
    stress_weights = read_tsv(output_root / "barneo_stress_guarded_method_weights.tsv")
    high_impact = read_tsv(output_root / "barneo_high_impact_lead_candidates.tsv")
    reviewer_kill = read_tsv(output_root / "barneo_high_impact_reviewer_kill_audit.tsv")

    figures = [
        fig_leaderboard(lb, fig_dir, asset_dir),
        fig_source(master, fig_dir, asset_dir),
        fig_funnel(master, barneo, bma, ctx, fail, xscore, fig_dir, asset_dir),
        fig_vulnerability(vuln, fig_dir, asset_dir),
        fig_challenges(challenge, fig_dir, asset_dir),
        fig_context(ctx, fig_dir, asset_dir),
        fig_patient(patient, fig_dir, asset_dir),
        fig_public(public, fig_dir, asset_dir),
        fig_split_heatmap(metrics, lb, fig_dir, asset_dir),
        fig_delta_anchor(metrics, lb, fig_dir, asset_dir),
        fig_context_weight_heatmap(ctx_weights, fig_dir, asset_dir),
        fig_claim_safe_top(x_explain, fig_dir, asset_dir),
        fig_stress_actions(stress_scores, fig_dir, asset_dir),
        fig_stress_top(stress_scores, fig_dir, asset_dir),
        fig_stress_method_weights(stress_weights, fig_dir, asset_dir),
    ]

    if not public.empty and "uses_public_pretraining" in public:
        public_rows = public[public["uses_public_pretraining"].astype(str).str.lower().isin(["true", "1"])]
    else:
        public_rows = public.iloc[0:0].copy()
    public_clean_allowed = (
        int(public_rows.get("clean_comparator_allowed_after_audit", pd.Series(dtype=bool)).astype(str).str.lower().isin(["true", "1"]).sum())
        if not public_rows.empty
        else 0
    )
    stats = {
        "Candidates": len(master),
        "Methods": safe_int(lb["method_name"].nunique()) if "method_name" in lb else 0,
        "Challenge Rows": safe_int(challenge["n_unique_candidates"].sum()) if "n_unique_candidates" in challenge else 0,
        "Contextual Clean Claims": int(ctx["contextual_clean_claim_allowed"].astype(bool).sum()) if "contextual_clean_claim_allowed" in ctx else 0,
        "X Priority": int((xscore.get("barneo_x_primary_action", pd.Series(dtype=str)) == "priority_review_candidate").sum()) if not xscore.empty else 0,
        "Stress Claim-Audit": int((stress_scores.get("stress_guarded_action", pd.Series(dtype=str)) == "claim_safe_candidate_after_manual_audit").sum()) if not stress_scores.empty else 0,
        "HI Leads": len(high_impact),
        "Kill Pass": int((reviewer_kill.get("reviewer_kill_disposition", pd.Series(dtype=str)) == "passes_current_reviewer_kill_audit").sum()) if not reviewer_kill.empty else 0,
        "Public Clean Allowed": public_clean_allowed,
    }
    write_html(hub_root, output_root, figures, stats)
    (output_root / "CLEAN_NEOBENCH_VISUAL_DASHBOARD.md").write_text(
        "# CLEAN-NeoBench Visual Dashboard\n\n"
        + "\n".join(f"- `figures/{f}`" for f in figures)
        + f"\n\nHTML: `{hub_root / PAGE_NAME}`\n"
    )

    outputs = ["CLEAN_NEOBENCH_VISUAL_DASHBOARD.md", PAGE_NAME] + [f"figures/{f}" for f in figures]
    update_manifest(output_root, "clean_neobench_visual_dashboard", {"outputs": outputs, "n_figures": len(figures), "warnings": ["Descriptive visual dashboard; not validation proof."]})
    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for name in outputs:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})["n_visual_dashboard_figures"] = len(figures)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"html": str(hub_root / PAGE_NAME), "figures": figures, "n_figures": len(figures)}, indent=2))


if __name__ == "__main__":
    main()
