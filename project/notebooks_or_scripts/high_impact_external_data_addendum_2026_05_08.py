#!/usr/bin/env python3
"""Summarize extra public data pulled during the 11-topic pilot sprint."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08"
DATA = OUT / "data_search"
REPORT = ROOT / "project/reports/2026_05_08_high_impact_extra_data_addendum.md"


def fmt(x) -> str:
    if pd.isna(x):
        return "NA"
    if isinstance(x, (float, np.floating)):
        if abs(x) < 1e-4 and x != 0:
            return f"{x:.1e}"
        return f"{x:.3g}"
    return str(x)


def main() -> None:
    cbio_arm = pd.read_csv(DATA / "cbio_class6_dm2_dm1_armdriver_validation.tsv", sep="\t")
    gene_tcga = pd.read_csv(DATA / "cbio_class6_dm2_dm1_selected_gene_gistic.tsv", sep="\t")

    summaries = []
    for study in ["thyroid_gatci_2024", "thyroid_mskcc_2016", "thpa_tcga_gdc"]:
        g = pd.read_csv(DATA / f"{study}_selected_gene_gistic.tsv", sep="\t")
        gm = pd.read_csv(DATA / "cbio_gene_map.tsv", sep="\t")[["entrezGeneId", "hugoGeneSymbol"]]
        g = g.merge(gm, on="entrezGeneId", how="left")
        g["gain_or_amp"] = g["alteration"].astype(float) >= 1
        g["loss_or_del"] = g["alteration"].astype(float) <= -1
        piv = (
            g.groupby("hugoGeneSymbol")
            .agg(
                n=("sampleId", "nunique"),
                gain_rate=("gain_or_amp", "mean"),
                loss_rate=("loss_or_del", "mean"),
                mean_alt=("alteration", "mean"),
            )
            .reset_index()
            .sort_values("gain_rate", ascending=False)
        )
        piv.insert(0, "study", study)
        summaries.append(piv)
    cross = pd.concat(summaries, ignore_index=True)
    cross.to_csv(OUT / "external_cna_selected_gene_summary.tsv", sep="\t", index=False)

    topic_update = pd.DataFrame(
        [
            {
                "topic_id": "T00",
                "old_grade": "B+",
                "new_grade": "A-",
                "reason": "cBioPortal ARMDRIVER_CN independently supports 7-arm Class6-DM2 enrichment: 35.4% vs 1.3%, Fisher p=1.6e-7.",
            },
            {
                "topic_id": "T01",
                "old_grade": "B+",
                "new_grade": "B+",
                "reason": "GSE301163 already included; extra cBio data adds context but not independent DICER1/DGCR8 class validation.",
            },
            {
                "topic_id": "T11",
                "old_grade": "B-",
                "new_grade": "B",
                "reason": "GATCI 2024 ATC GISTIC shows frequent EGFR/MET/TERT/CDK6/BRAF gains, supporting an advanced-CNA progression layer, but not the same PTC residual-class test.",
            },
        ]
    )
    topic_update.to_csv(OUT / "topic_grade_updates_after_extra_data.tsv", sep="\t", index=False)

    lines = []
    lines.append("# 추가 public data addendum — 2026-05-08\n")
    lines.append("사용자 요청에 따라 11-topic pilot 중 추가 public 데이터를 더 가져와서 바로 붙일 수 있는 검증 레이어를 확인했다.\n")

    lines.append("## 새로 가져온 데이터\n")
    lines.append("| Source | What | Local file | Use |")
    lines.append("|---|---|---|---|")
    lines.append("| cBioPortal TCGA-THCA Cell 2014 | sample clinical, ARMDRIVER_CN, ARM_SCNA_CLUSTER, selected-gene GISTIC | `data_search/cbio_thca_tcga_pub_*` | T00 CNV residual class independent validation |")
    lines.append("| cBioPortal TCGA-GDC 2025 | selected-gene GISTIC, hg38 GDC mirror | `data_search/thpa_tcga_gdc_*` | T00 consistency check |")
    lines.append("| cBioPortal GATCI 2024 ATC | n=190 ATC GISTIC + clinical | `data_search/thyroid_gatci_2024_*` | T11 advanced CNA progression context |")
    lines.append("| cBioPortal MSK 2016 PDTC/ATC | n=117 targeted-panel GISTIC + clinical | `data_search/thyroid_mskcc_2016_*` | negative/sparse targeted-panel context |")
    lines.append("| Web snapshots | JCI Insight 2026 atlas + GSE301163 OmicsDI page | `data_search/external_source_snapshots.json` | prior-art guardrail |")

    lines.append("\n## T00 CNV residual class — cBioPortal validation\n")
    lines.append("| Metric | DM2 rate | DM1 rate | OR | p |")
    lines.append("|---|---:|---:|---:|---:|")
    for row in cbio_arm.head(10).itertuples(index=False):
        lines.append(f"| {row.metric} | {fmt(row.DM2_rate)} | {fmt(row.DM1_rate)} | {fmt(row.OR)} | {fmt(row.p)} |")
    lines.append("\n해석: local arm-CNV table뿐 아니라 cBioPortal의 TCGA-THCA sample-level `ARMDRIVER_CN`에서도 같은 7-arm signature가 DM2에 강하게 몰린다. T00은 `B+`에서 `A-`로 상향 가능하다.\n")

    lines.append("## T00 selected-gene GISTIC support within Class6\n")
    lines.append("| Gene | Event | DM2 rate | DM1 rate | OR | p |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for row in gene_tcga.head(18).itertuples(index=False):
        lines.append(f"| {row.gene} | {row.event} | {fmt(row.DM2_rate)} | {fmt(row.DM1_rate)} | {fmt(row.OR)} | {fmt(row.p)} |")
    lines.append("\n해석: 7q/7p representative genes (MET, CDK6, BRAF, EGFR)와 12q/16p 관련 selected genes가 Class6-DM2에서 gene-level GISTIC gain으로도 반복된다. 단 gene-level은 arm-level의 proxy라서 main claim은 arm-CNV로 유지.\n")

    lines.append("## Advanced thyroid cancer external context\n")
    lines.append("| Study | Top selected CNA gains | Interpretation |")
    lines.append("|---|---|---|")
    for study in ["thyroid_gatci_2024", "thyroid_mskcc_2016", "thpa_tcga_gdc"]:
        top = cross[cross["study"] == study].sort_values("gain_rate", ascending=False).head(8)
        top_txt = ", ".join([f"{r.hugoGeneSymbol} {r.gain_rate:.1%}" for r in top.itertuples(index=False)])
        if study == "thyroid_gatci_2024":
            interp = "ATC에서 EGFR/MET/TERT/CDK6/BRAF gains가 흔해 advanced-CNA progression context를 제공."
        elif study == "thyroid_mskcc_2016":
            interp = "Targeted-panel GISTIC가 매우 sparse; negative control/coverage caveat로만 사용."
        else:
            interp = "TCGA-GDC mirror에서 TCGA selected-gene gain frequency가 비슷한 방향으로 확인."
        lines.append(f"| {study} | {top_txt} | {interp} |")

    lines.append("\n## 반영된 판정\n")
    lines.append("| Topic | Old | New | Reason |")
    lines.append("|---|---:|---:|---|")
    for row in topic_update.itertuples(index=False):
        lines.append(f"| {row.topic_id} | {row.old_grade} | {row.new_grade} | {row.reason} |")

    lines.append("\n## Files\n")
    lines.append(f"- `{OUT / 'external_cna_selected_gene_summary.tsv'}`\n")
    lines.append(f"- `{OUT / 'topic_grade_updates_after_extra_data.tsv'}`\n")
    lines.append(f"- `{REPORT}`\n")
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"[write] {OUT / 'external_cna_selected_gene_summary.tsv'}")
    print(f"[write] {OUT / 'topic_grade_updates_after_extra_data.tsv'}")
    print(f"[write] {REPORT}")


if __name__ == "__main__":
    main()
