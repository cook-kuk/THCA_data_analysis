"""Fig 10 — chr7 arm-level CNA forest + I-131 dose by DM + within-z K2 calibration."""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
R9 = ROOT / "project/results/audit_2026_04_30/round9"
FIG = ROOT / "project/reports/html/figs_interactive/v17"


def main() -> None:
    # data
    boot = json.loads((R9 / "r9_3_summary.json").read_text())["bootstrap"]
    i131 = pd.DataFrame(json.loads((R9 / "r9_1_summary.json").read_text())["i131_by_dm"])
    k2_z = pd.read_csv(R9 / "r9_4_k2_within_z.tsv", sep="\t")
    rtk = pd.read_csv(R9 / "r9_2_chr7_rtk_summary.tsv", sep="\t")

    fig = make_subplots(
        rows=1,
        cols=4,
        subplot_titles=(
            "A. chr arm Gain Δ (DM1 − not_DM) bootstrap CI",
            "B. chr7 RTK gene log2CNA Cohen's d",
            "C. I-131 total dose by DM (n with dose / pct)",
            "D. K2 within-sample-z calibration (TPM-fix)",
        ),
        column_widths=[0.27, 0.23, 0.25, 0.25],
    )

    # A — forest plot
    arms = list(boot.keys())
    means = [boot[a]["mean_diff"] for a in arms]
    lows = [boot[a]["ci_low"] for a in arms]
    highs = [boot[a]["ci_high"] for a in arms]
    err_minus = [m - lo for m, lo in zip(means, lows)]
    err_plus = [hi - m for m, hi in zip(means, highs)]
    fig.add_trace(
        go.Scatter(
            x=means,
            y=arms,
            mode="markers",
            error_x=dict(type="data", array=err_plus, arrayminus=err_minus, color="#444"),
            marker=dict(size=14, color=["#d62728" if a in ("7p", "7q") else "#888" for a in arms]),
            showlegend=False,
        ),
        row=1, col=1,
    )
    fig.add_vline(x=0, line=dict(color="grey", dash="dot"), row=1, col=1)
    fig.update_xaxes(title="Δ rate (gain%)", row=1, col=1)

    # B — RTK Cohen's d bar
    rtk_sorted = rtk.sort_values("cohens_d_dm1_vs_notDM", ascending=True)
    fig.add_trace(
        go.Bar(
            y=rtk_sorted["gene"],
            x=rtk_sorted["cohens_d_dm1_vs_notDM"],
            orientation="h",
            marker_color=[
                "#d62728" if g in ("BRAF", "EGFR", "MET") else ("#888" if g != "RET" else "#1f77b4")
                for g in rtk_sorted["gene"]
            ],
            text=[f"d={d:.2f}" for d in rtk_sorted["cohens_d_dm1_vs_notDM"]],
            textposition="auto",
            showlegend=False,
        ),
        row=1, col=2,
    )
    fig.add_vline(x=0, line=dict(color="grey", dash="dot"), row=1, col=2)
    fig.update_xaxes(title="Cohen's d (DM1 vs not_DM)", row=1, col=2)

    # C — I-131 dose by DM
    fig.add_trace(
        go.Bar(
            x=i131["cluster"],
            y=i131["pct_received_i131"],
            marker_color=["#d62728", "#ff7f0e", "#888"],
            text=[
                f"{r['n_with_i131']}/{r['n_total']} ({r['pct_received_i131']}%)<br>median {r['median_mCi'] or 'NA'} mCi"
                for _, r in i131.iterrows()
            ],
            textposition="auto",
            showlegend=False,
        ),
        row=1, col=3,
    )
    fig.update_yaxes(title="% with I-131 dose recorded", row=1, col=3)

    # D — K2 within-z calibration: scatter rai_within_z_top4 vs p_DM2
    fig.add_trace(
        go.Scatter(
            x=k2_z["rai_within_z_top4"],
            y=k2_z["p_DM2"],
            mode="markers",
            marker=dict(
                color=k2_z["rai_within_z_top4"], colorscale="RdBu_r",
                size=6, showscale=True,
                colorbar=dict(title="rai z (top4)", x=1.02),
            ),
            showlegend=False,
            text=k2_z["DM_call"],
        ),
        row=1, col=4,
    )
    fig.update_xaxes(title="K2 within-sample-z (top4 RAI genes)", row=1, col=4)
    fig.update_yaxes(title="p_DM2 (original absolute-form)", row=1, col=4)

    fig.update_layout(
        height=560, width=1700,
        title="Fig 10. chr7 arm-level CNA forest + chr7 RTK Cohen's d + TCGA I-131 dose by DM + K2 within-z calibration",
        template="plotly_white", margin=dict(l=70, r=120, t=80, b=70),
    )
    out_html = FIG / "v17_fig10_R9_chr7_RAI_K2.html"
    fig.write_html(str(out_html), include_plotlyjs="cdn")
    print(f"wrote {out_html}")


if __name__ == "__main__":
    main()
