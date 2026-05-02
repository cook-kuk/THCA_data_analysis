"""Fig 11 — R10 integrated: TCGA Cabrita TLS by DM, 3-way survival, chr7 × fusion within DM1, MSK chr7."""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
R10 = ROOT / "project/results/audit_2026_04_30/round10"
FIG = ROOT / "project/reports/html/figs_interactive/v17"


def main() -> None:
    tls_summary = json.loads((R10 / "r10_5_summary.json").read_text())
    tls_per = pd.read_csv(R10 / "r10_5_tcga_tls_per_sample.tsv", sep="\t")
    surv = pd.read_csv(R10 / "r10_2_3way_combo_summary.tsv", sep="\t")
    chr7 = pd.read_csv(R10 / "r10_3_chr7_x_fusion_4cell.tsv", sep="\t")
    msk = json.loads((R10 / "r10_1_msk_chr7_replication.json").read_text())

    fig = make_subplots(
        rows=1, cols=4,
        subplot_titles=(
            f"A. TCGA Cabrita 12-gene TLS by DM<br>(DM1 vs DM2 d={tls_summary['cohens_d_DM1_vs_DM2']}, p={tls_summary['mw_p_DM1_vs_DM2']:.1e})",
            "B. 3-way DM × fusion × TERT survival<br>(median OS days)",
            "C. chr7 gain × RET-fusion within DM1<br>(co-occurrence)",
            "D. MSK 2016 GISTIC discrete<br>(chr7 RTK, % Gain/Amp)",
        ),
        column_widths=[0.27, 0.28, 0.23, 0.22],
    )

    # A: TLS violin
    color_map = {"DM1": "#d62728", "DM2": "#ff7f0e", "not_DM": "#888"}
    master = pd.read_csv(
        ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv", sep="\t"
    )
    master["sample_short"] = master["sample_id"].str.slice(0, 15)
    master["dm"] = master["v17_dark_cluster"].fillna("not_DM").replace("", "not_DM")
    merged = tls_per.merge(master[["sample_short", "dm"]], on="sample_short", how="inner")
    for cluster in ["DM1", "DM2", "not_DM"]:
        sub = merged[merged["dm"] == cluster]["TLS_score"].dropna()
        fig.add_trace(
            go.Violin(
                y=sub, name=f"{cluster} (n={len(sub)})",
                line_color=color_map[cluster], fillcolor=color_map[cluster],
                opacity=0.5, box_visible=True, meanline_visible=True,
            ),
            row=1, col=1,
        )
    fig.update_yaxes(title="TLS z-score (mean of 12 genes)", row=1, col=1)

    # B: bar chart median OS by 8 combos
    surv_sorted = surv.sort_values("median_os")
    combo_cols = []
    for c in surv_sorted["combo"]:
        if c.startswith("DM1_fus1"):
            combo_cols.append("#d62728")
        elif c.startswith("DM1"):
            combo_cols.append("#ff7f0e")
        elif c.startswith("notDM1_fus1"):
            combo_cols.append("#1f77b4")
        else:
            combo_cols.append("#888")
    fig.add_trace(
        go.Bar(
            y=surv_sorted["combo"],
            x=surv_sorted["median_os"],
            orientation="h",
            marker_color=combo_cols,
            text=[f"n={n} ({int(e)} ev)" for n, e in zip(surv_sorted["n"], surv_sorted["events"])],
            textposition="auto",
            showlegend=False,
        ),
        row=1, col=2,
    )
    fig.update_xaxes(title="Median OS (days)", row=1, col=2)

    # C: chr7 × fusion within DM1 (4-cell)
    dm1_only = chr7[chr7["cohort"] == "DM1"].copy()
    dm1_only["label"] = dm1_only.apply(lambda r: f"chr7={r['chr7_gain']}|fus={r['fusion']}", axis=1)
    fig.add_trace(
        go.Bar(
            y=dm1_only["label"],
            x=dm1_only["n"],
            orientation="h",
            marker_color=["#d62728" if (r["chr7_gain"] == 1 and r["fusion"] == 1) else "#888"
                          for _, r in dm1_only.iterrows()],
            text=[f"n={r['n']} | rai med={r['rai_med'] if r['rai_med'] else 'NA'}" for _, r in dm1_only.iterrows()],
            textposition="auto",
            showlegend=False,
        ),
        row=1, col=3,
    )
    fig.update_xaxes(title="DM1 sample count", row=1, col=3)

    # D: MSK + GATCI chr7 RTK GISTIC %
    msk_rows = []
    for study, payload in msk.items():
        if "per_gene" not in payload:
            continue
        for g in payload["per_gene"]:
            msk_rows.append({"study": study, **g})
    msk_df = pd.DataFrame(msk_rows)
    if not msk_df.empty:
        for study in msk_df["study"].unique():
            sub = msk_df[msk_df["study"] == study]
            fig.add_trace(
                go.Bar(
                    x=sub["gene"],
                    y=sub["gain_amp_pct"],
                    name=f"{study} (n={int(sub['n'].iloc[0])})",
                ),
                row=1, col=4,
            )
        fig.update_xaxes(row=1, col=4)
        fig.update_yaxes(title="% Gain/Amp", row=1, col=4)

    fig.update_layout(
        height=560, width=1800, barmode="group",
        title="Fig 11. R10 — TCGA Cabrita TLS replication + 3-way survival + chr7 × fusion DM1 4-cell + MSK GISTIC chr7 RTK",
        template="plotly_white", margin=dict(l=70, r=20, t=110, b=60),
    )
    out = FIG / "v17_fig11_R10_TLS_3way_chr7_MSK.html"
    fig.write_html(str(out), include_plotlyjs="cdn")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
