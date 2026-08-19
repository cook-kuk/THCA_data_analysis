#!/usr/bin/env python3
from __future__ import annotations

import html
import math
import shutil
import textwrap
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch
import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PROP = ROOT / "kthyro_ctc_emt_samsung_proposal"
NOYU = PROP / "outputs/in_silico_paper_no_yu"
SAMSUNG = PROP / "outputs"
PILOT = ROOT / "kthyro_public_pilot/results"
OUT = PROP / "outputs/paper_flow_master_dossier"
FIG = OUT / "figures"
TAB = OUT / "tables"
REP = OUT / "reports"
SRC = OUT / "source_assets"

HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = "kthyro_ctc_emt_paper_flow"
HUB_ASSET = HUB / "assets" / ASSET
LIVE_ASSET = LIVE / "assets" / ASSET
HTML = HUB / "kthyro_ctc_emt_paper_flow_master_dossier.html"
LIVE_HTML = LIVE / "kthyro_ctc_emt_paper_flow_master_dossier.html"
LOCAL_HTML = REP / "kthyro_ctc_emt_paper_flow_master_dossier.html"
ZIP_NAME = "kthyro_ctc_emt_paper_flow_packet_2026_05_09.zip"
ZIP_PATH = OUT / ZIP_NAME

BG = "#080d18"
PANEL = "#101a2c"
PANEL2 = "#16213a"
LINE = "#263954"
TEXT = "#e9f0fb"
MUTED = "#9aa9bd"
GOLD = "#ffd28a"
CYAN = "#5bdcff"
TEAL = "#35d39d"
PURPLE = "#b38cff"
ORANGE = "#ff9b66"
WARN = "#ff8066"
GREEN = "#64e0aa"

KR_FONT = Path("/home/seungho/.local/share/fonts/NanumGothic-Regular.ttf")
KR_FONT_BOLD = Path("/home/seungho/.local/share/fonts/NanumGothic-Bold.ttf")
if KR_FONT.exists():
    fm.fontManager.addfont(str(KR_FONT))
    if KR_FONT_BOLD.exists():
        fm.fontManager.addfont(str(KR_FONT_BOLD))
    plt.rcParams["font.family"] = "NanumGothic"


def mkdirs() -> None:
    for d in [OUT, FIG, TAB, REP, SRC, HUB_ASSET, LIVE_ASSET]:
        d.mkdir(parents=True, exist_ok=True)


def esc(x) -> str:
    return html.escape("" if x is None else str(x), quote=True)


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def fmt(x, nd=3) -> str:
    try:
        if pd.isna(x):
            return ""
        v = float(x)
        if abs(v) >= 1000:
            return f"{v:,.0f}"
        if abs(v) >= 100:
            return f"{v:.0f}"
        return f"{v:.{nd}g}"
    except Exception:
        return str(x)


def df_table(df: pd.DataFrame, cols: list[str] | None = None, max_rows: int = 18) -> str:
    if df.empty:
        return "<p class='muted'>No table available.</p>"
    view = df.copy()
    if cols:
        cols = [c for c in cols if c in view.columns]
        view = view[cols]
    view = view.head(max_rows)
    out = ["<table class='t'><thead><tr>"]
    for c in view.columns:
        out.append(f"<th>{esc(c)}</th>")
    out.append("</tr></thead><tbody>")
    for _, row in view.iterrows():
        out.append("<tr>")
        for c in view.columns:
            v = row[c]
            numeric = isinstance(v, (int, float)) and not isinstance(v, bool)
            out.append(f"<td{' class=num' if numeric else ''}>{esc(fmt(v) if numeric else v)}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def setup(name: str, title: str, subtitle: str = "", w: float = 15.5, h: float = 8.2):
    fig, ax = plt.subplots(figsize=(w, h), dpi=180)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.035, 0.925, title, color=TEXT, fontsize=25, fontweight="bold", transform=ax.transAxes)
    if subtitle:
        ax.text(0.037, 0.872, subtitle, color=MUTED, fontsize=11.5, transform=ax.transAxes)
    return fig, ax, FIG / name


def save(fig, path: Path) -> None:
    fig.savefig(path, facecolor=BG, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def box(ax, x, y, w, h, title, sub="", color=CYAN, fs=13, align="center", fill=PANEL):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.014",
        linewidth=1.4,
        edgecolor=color,
        facecolor=fill,
        transform=ax.transAxes,
    )
    ax.add_patch(patch)
    tx = x + (w / 2 if align == "center" else 0.018)
    ha = align
    ax.text(tx, y + h * 0.67, title, color=TEXT, fontsize=fs, fontweight="bold", ha=ha, va="center", transform=ax.transAxes)
    if sub:
        wrapped = "\n".join(textwrap.wrap(sub, width=max(18, int(w * 82))))
        ax.text(tx, y + h * 0.30, wrapped, color=MUTED, fontsize=max(fs - 5, 7.5), ha=ha, va="center", transform=ax.transAxes)


def arrow(ax, x1, y1, x2, y2, color=GOLD):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        xycoords=ax.transAxes,
        textcoords=ax.transAxes,
        arrowprops=dict(arrowstyle="-|>", color=color, lw=2, shrinkA=5, shrinkB=5),
    )


def collect_metrics() -> dict:
    braf = read_tsv(PILOT / "extra_analyses/tables/tcga_braf_only_label_distribution.tsv")
    var = read_tsv(PILOT / "extra_analyses/tables/tcga_axis_variance_models.tsv")
    spatial = read_tsv(PILOT / "tables/spatial_coherence_statistics.tsv")
    scrna = read_tsv(PILOT / "extra_analyses/tables/scrna_module_celltype_attribution.tsv")
    prot = read_tsv(PILOT / "extra_analyses/tables/bulk_proteomics_kthyro_module_contrasts.tsv")
    metrics: dict[str, object] = {}
    if not braf.empty:
        metrics["braf_n"] = int(braf["n"].sum())
        metrics["braf_labels"] = int(braf.shape[0])
        metrics["braf_largest_pct"] = float(braf["fraction_of_braf"].max() * 100)
        metrics["braf_ambiguity_pct"] = float((1 - braf["fraction_of_braf"].max()) * 100)
    else:
        metrics.update(braf_n="NA", braf_labels="NA", braf_largest_pct=0, braf_ambiguity_pct=0)
    if not spatial.empty:
        niche = spatial[spatial["label_set"].eq("niche_label")].dropna(subset=["z"])
        metrics["spatial_slides"] = int(niche["sample_id"].nunique())
        metrics["spatial_z_gt2"] = int((niche["z"] > 2).sum())
        metrics["spatial_z_median"] = float(niche["z"].median())
    else:
        metrics.update(spatial_slides="NA", spatial_z_gt2="NA", spatial_z_median=0)
    spot_file = PILOT / "tables/spatial_spot_vulnerability_scores.tsv"
    metrics["spatial_spots"] = 57144
    if spot_file.exists():
        try:
            metrics["spatial_spots"] = sum(1 for _ in spot_file.open()) - 1
        except OSError:
            pass
    metrics["scrna_modules"] = int(scrna.shape[0]) if not scrna.empty else "NA"
    metrics["proteomics_rows"] = int(prot.shape[0]) if not prot.empty else "NA"
    if not var.empty and "r_squared" in var.columns:
        metrics["driver_model_count"] = int(var.shape[0])
        metrics["driver_only_r2_median_pct"] = float(var.loc[var["model"].eq("driver_only"), "r_squared"].median() * 100)
    else:
        metrics["driver_model_count"] = "NA"
        metrics["driver_only_r2_median_pct"] = 0
    return metrics


def write_master_tables(metrics: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    data_rows = [
        {
            "dataset_id": "TCGA-THCA",
            "type": "bulk tumor RNA/mutation/clinical public cohort",
            "local_input": "kthyro_public_pilot/results/tables + extra_analyses/tables",
            "n_or_scale": "505 primary tumors in pilot context; BRAF-mutant subset 274",
            "used_for": "mutation is not enough; tissue-state heterogeneity; driver/state axis support",
            "key_variables": "driver group, vulnerability/state label, RAI/HLA/immune/CD8 axes, covariates",
            "allowed_claim": "public PTC tissue states split beyond mutation group",
            "forbidden_claim": "does not prove CTC shedding or post-thyroidectomy transition",
        },
        {
            "dataset_id": "GSE250521 spatial transcriptomics",
            "type": "thyroid spatial transcriptomics",
            "local_input": "spatial_coherence_statistics.tsv; spatial maps/adjacency/hotspot tables",
            "n_or_scale": f"{metrics.get('spatial_slides')} slides; {metrics.get('spatial_spots')} spots",
            "used_for": "spatial tissue-state organization and ROI/marker logic",
            "key_variables": "slide, spot, niche label, module score, same-neighbor coherence z",
            "allowed_claim": "tissue states are spatially organized in public thyroid tissue data",
            "forbidden_claim": "spatial tissue map is not CTC origin proof",
        },
        {
            "dataset_id": "PTC scRNA marker context",
            "type": "single-cell RNA reference summarized in pilot tables",
            "local_input": "scrna_module_celltype_attribution.tsv; scRNA marker figures",
            "n_or_scale": f"{metrics.get('scrna_modules')} module rows in current summary",
            "used_for": "cell-type context for epithelial, mesenchymal, thyroid-lineage, stress modules",
            "key_variables": "module, top cell type, top-minus-second attribution, cell count",
            "allowed_claim": "marker modules have interpretable tissue/cell context",
            "forbidden_claim": "does not itself define circulating CTC state",
        },
        {
            "dataset_id": "K-THYRO public bulk proteomics proxy",
            "type": "bulk proteomics thyroid cancer state proxy",
            "local_input": "bulk_proteomics_kthyro_module_contrasts.tsv; dediff_trends.tsv",
            "n_or_scale": "461 thyroid samples in local proxy summary",
            "used_for": "orthogonal protein-direction support for RAI loss and myeloid/TGFB gain",
            "key_variables": "module, contrast, Cohen's d, FDR, dedifferentiation-order rho",
            "allowed_claim": "protein-level direction is consistent with tissue-state prior",
            "forbidden_claim": "not spatial proteomics and not blood/CTC measurement",
        },
        {
            "dataset_id": "Yu et al. 2024 CTC study",
            "type": "published prospective PTC serial CTC anchor",
            "local_input": "literature anchor summarized in Samsung pitch pack; raw patient table not local",
            "n_or_scale": "62 prospective PTC patients; pre-op, 2 weeks, 3 months; 87% CTC detection reported",
            "used_for": "feasibility anchor that PTC CTC-EMT states are measurable around surgery",
            "key_variables": "timepoint, CTC count, epithelial/hybrid/mesenchymal phenotype",
            "allowed_claim": "published CTC-EMT measurement around thyroidectomy is feasible",
            "forbidden_claim": "we cannot reanalyze patient-level transition without raw Yu/hospital data",
        },
        {
            "dataset_id": "Future Samsung/Yu prospective cohort",
            "type": "planned clinical research dataset",
            "local_input": "not yet collected",
            "n_or_scale": "serial blood T0/T1/T2/T3 plus tissue-normal NGS/cfDNA/CTC-enriched NGS",
            "used_for": "test primary biological endpoint",
            "key_variables": "CTC E/EM/M, tissue clone, cfDNA persistence, Tg/US/pathology/outcome",
            "allowed_claim": "proposal will test thyroidectomy as perturbation",
            "forbidden_claim": "no completed validation or clinical diagnostic claim yet",
        },
        {
            "dataset_id": "Inocras/vendor workflow",
            "type": "planned sequencing/assay service inputs",
            "local_input": "vendor questionnaire and budget reports",
            "n_or_scale": "not a novelty source; execution infrastructure only",
            "used_for": "assay feasibility, QC, sample logistics, low-input CTC/cfDNA testing",
            "key_variables": "input amount, QC pass rate, panel coverage, turnaround, failure handling",
            "allowed_claim": "vendor can support implementation if QC thresholds pass",
            "forbidden_claim": "vendor technology is not the scientific novelty",
        },
    ]
    data = pd.DataFrame(data_rows)
    data.to_csv(TAB / "data_dictionary_master.tsv", sep="\t", index=False)

    fig_rows = [
        ["Fig. 1", "Clinical gap and perturbation model", "PF01_manuscript_flow.png; slide01/02 figures", "Yu 2024 literature + proposal logic", "Introduction / Rationale", "Thyroidectomy is a controlled human perturbation for serial CTC-EMT monitoring", "not a clinical diagnostic"],
        ["Fig. 2", "Public data layer/provenance map", "PF02_data_provenance_map.png; F02_public_data_layer_map.png", "TCGA, spatial, scRNA, proteomics, Yu literature", "Results 1 / Methods overview", "public data can construct a prior map", "does not prove CTC biology"],
        ["Fig. 3", "Mutation alone is insufficient", "F03_tcga_braf_state_split.png; tcga_braf_only_vulnerability_axis_heatmap.png", "TCGA-THCA BRAF-mutant subset", "Results 2", "BRAF tumors split across tissue states", "no recurrence/CTC claim"],
        ["Fig. 4", "Spatial tissue-state organization", "F05_spatial_tissue_state_organization.png; spatial_coherence_summary.png", "GSE250521 spatial transcriptomics", "Results 3", "tissue states are spatially coherent", "not a shedding origin proof"],
        ["Fig. 5", "Single-cell/proteomics marker support", "F06_scrna_marker_context.png; F07_proteomics_dediff_direction.png", "scRNA module table + bulk proteomics proxy", "Results 4", "marker modules have tissue/cell/protein rationale", "not CTC validation"],
        ["Fig. 6", "CTC-EMT-NGS prior panel", "F08_ctc_emt_ngs_prior_panel.png; PF03_marker_to_ngs_map.png", "marker prior table + thyroid cancer drivers", "Results 5 / Proposal bridge", "define research-grade phenotype/genotype panel", "not a locked clinical assay"],
        ["Fig. 7", "Dataset gap and prospective design", "PF04_dataset_gap_bridge.png; slide08/09 figures", "public dataset gap table + Samsung study design", "Discussion / Future validation", "integrated serial CTC-EMT-NGS dataset is missing and needed", "not solved by public data"],
        ["Fig. 8", "Claim boundary and kill criteria", "PF05_claim_boundary_matrix.png; F09_publishability_decision.png", "claim boundary reports and reviewer defense", "Discussion / Reviewer defense", "explicit boundaries make proposal credible", "no overclaiming"],
        ["Table 1", "Data dictionary", "data_dictionary_master.tsv", "all local/public/planned datasets", "Methods / Supplement", "each data layer has defined role and limits", "no untracked source"],
        ["Table 2", "Figure/table manifest", "figure_table_manifest.tsv", "all selected outputs", "Supplement", "every figure maps to a source and claim", "no floating figure"],
        ["Table 3", "Manuscript flow outline", "manuscript_flow_outline.tsv", "proposal/no-Yu outputs", "Writing plan", "result sequence is fixed", "not final prose"],
        ["Table 4", "Claim boundary matrix", "claim_boundary_matrix.tsv", "claim-boundary reports", "Reviewer defense", "allowed/forbidden statements are separated", "not defensive handwaving"],
    ]
    manifest = pd.DataFrame(fig_rows, columns=["paper_item", "title", "asset_or_table", "data_used", "manuscript_position", "allowed_takeaway", "claim_boundary"])
    manifest.to_csv(TAB / "figure_table_manifest.tsv", sep="\t", index=False)

    flow_rows = [
        ["Title/Abstract", "Name the paper as a framework/prior-map paper", "No-Yu data can support a translational rationale, not discovery", "Fig. 1 summary schematic", "Do not say postoperative transition was discovered"],
        ["Introduction paragraph 1", "PTC is usually curable but residual-risk biology is hard to read after surgery", "creates clinical gap", "Fig. 1A", "avoid implying current tools are useless"],
        ["Introduction paragraph 2", "Yu 2024 shows CTC-EMT measurement around thyroidectomy is feasible", "published anchor", "literature anchor box", "do not present Yu as our dataset"],
        ["Introduction paragraph 3", "No public dataset integrates serial CTC-EMT, matched NGS, cfDNA, tissue context, and outcome", "main gap", "Fig. 7 / Table gap", "avoid saying no one has ever studied CTCs"],
        ["Results 1", "Public multi-omics can build a prior but not test the transition", "data provenance first", "Fig. 2", "no CTC conclusion"],
        ["Results 2", "BRAF/mutation alone does not resolve tissue-state heterogeneity", "TCGA BRAF split and variance models", "Fig. 3", "no causality"],
        ["Results 3", "Tissue states are spatially organized", "spatial coherence and adjacency", "Fig. 4", "not CTC origin"],
        ["Results 4", "Marker modules have cell and protein context", "scRNA/proteomics support", "Fig. 5", "not validated blood assay"],
        ["Results 5", "Construct CTC-EMT-NGS prior map", "panel logic and clone anchor", "Fig. 6", "research-grade, not diagnostic"],
        ["Discussion", "The real test requires serial hospital CTC/cfDNA/tissue NGS", "bridge to Samsung Science proposal", "Fig. 7", "future validation needed"],
        ["Reviewer defense", "Use hard claim boundaries and kill criteria", "credibility", "Fig. 8 / Table 4", "no recurrence prediction"],
    ]
    flow = pd.DataFrame(flow_rows, columns=["section", "story_step", "why_here", "figure_or_table", "boundary_sentence"])
    flow.to_csv(TAB / "manuscript_flow_outline.tsv", sep="\t", index=False)

    claim_rows = [
        ["Allowed", "Yu 2024 supports feasibility of serial CTC-EMT measurement in PTC", "Slide 3; literature anchor row", "published summary, not local raw reanalysis"],
        ["Allowed", "Public omics show thyroid cancer tissue-state heterogeneity beyond mutation", "Fig. 3; TCGA BRAF tables", "tissue-state support only"],
        ["Allowed", "Spatial data support organized tissue-state context", "Fig. 4; spatial coherence table", "not CTC origin proof"],
        ["Allowed", "scRNA/proteomics support marker-prior selection", "Fig. 5; marker table", "not blood validation"],
        ["Allowed", "Matched tissue-normal NGS is needed to interpret postoperative blood signals", "Fig. 6-7; NGS map", "proposal logic"],
        ["Forbidden", "Our cohort proves thyroidectomy-induced CTC-EMT transition", "none", "requires Yu/hospital raw serial CTC data"],
        ["Forbidden", "Public pilot proves CTC shedding", "none", "public pilot is tissue/bulk/spatial only"],
        ["Forbidden", "CTC/cfDNA NGS will work in every early PTC patient", "none", "must be stage-gated by QC"],
        ["Forbidden", "This is a ready clinical diagnostic", "none", "research-grade monitoring framework only"],
        ["Forbidden", "Inocras/vendor method is the novelty", "none", "vendor is infrastructure"],
    ]
    claims = pd.DataFrame(claim_rows, columns=["status", "statement", "support", "boundary"])
    claims.to_csv(TAB / "claim_boundary_matrix.tsv", sep="\t", index=False)
    return data, manifest, flow, claims


def figure_explanation_rows() -> list[dict[str, str]]:
    return [
        {
            "figure": "PF01",
            "image": "PF01_manuscript_flow.png",
            "title": "논문 전체 논리 흐름",
            "role": "첫 장에서 독자에게 이 논문이 CTC discovery가 아니라 prior-map/framework 논문임을 고정한다.",
            "data_used": "Yu 2024 published CTC anchor, Samsung Science proposal logic, public-pilot evidence layers.",
            "how_to_read": "왼쪽 clinical gap에서 시작해 Yu 2024 feasibility, public in silico prior, future prospective test로 이동한다. 아래 네 result box는 본문 Results의 순서다.",
            "paper_use": "Introduction 마지막 또는 Results overview 첫 그림.",
            "narration": "이 그림은 '왜 지금 공개 데이터만으로도 쓸 말이 있는가'와 '왜 병원 serial CTC 데이터가 여전히 필요한가'를 동시에 보여준다.",
            "boundary": "이 그림만으로 thyroidectomy 후 CTC-EMT transition이 증명됐다고 말하지 않는다.",
        },
        {
            "figure": "PF02",
            "image": "PF02_data_provenance_map.png",
            "title": "데이터 출처와 역할 정의",
            "role": "각 데이터가 무슨 일을 하는지, 그리고 어디까지 못 하는지를 한 화면에 고정한다.",
            "data_used": "TCGA-THCA, GSE250521 spatial, PTC scRNA, K-THYRO proteomics proxy, Yu 2024, future cohort, vendor workflow.",
            "how_to_read": "각 박스의 첫 줄은 데이터 레이어, 둘째 줄은 사용 목적, 마지막 limit 문장은 금지 claim이다.",
            "paper_use": "Methods/Data overview와 Reviewer defense에 사용.",
            "narration": "공개 데이터는 CTC를 직접 보여주는 게 아니라 tissue-state와 marker/ROI/genotype prior를 만든다.",
            "boundary": "public pilot이 CTC shedding을 증명한다는 문장은 금지.",
        },
        {
            "figure": "PF03",
            "image": "PF03_marker_to_ngs_map.png",
            "title": "Marker-to-NGS 해석 구조",
            "role": "CTC phenotype만으로 과장하지 않고, clone anchoring을 요구하는 과학적 안전장치.",
            "data_used": "CTC-EMT marker prior table, thyroid-lineage markers, thyroid cancer driver/private variant logic.",
            "how_to_read": "Epithelial/hybrid/mesenchymal/thyroid-lineage/survival modules가 CTC fraction state로 모이고, 오른쪽 tumor-informed NGS가 clone identity를 확인한다.",
            "paper_use": "Results 5 또는 prospective study design의 핵심 schematic.",
            "narration": "표현형은 후보 신호이고, NGS는 그 신호가 진짜 tumor-derived인지 묻는 장치다.",
            "boundary": "marker co-expression을 viable metastasis나 recurrence predictor로 말하지 않는다.",
        },
        {
            "figure": "PF04",
            "image": "PF04_dataset_gap_bridge.png",
            "title": "빠진 데이터셋과 제안서 필요성",
            "role": "왜 Samsung/Yu prospective cohort가 필요한지 가장 직접적으로 설명한다.",
            "data_used": "serial CTC-NGS dataset gap map, public dataset inventory, Samsung study design.",
            "how_to_read": "위의 다섯 구성요소가 public dataset에서 동시에 빠져 있다. 아래는 public data가 할 수 있는 일과 prospective cohort가 해야 할 일을 나눈다.",
            "paper_use": "Discussion 또는 proposal justification.",
            "narration": "공개 데이터는 prior를 만들지만, thyroidectomy라는 perturbation test는 새 cohort 없이는 불가능하다.",
            "boundary": "public data만으로 postoperative outcome association을 주장하지 않는다.",
        },
        {
            "figure": "PF05",
            "image": "PF05_claim_boundary_matrix.png",
            "title": "허용 claim / 금지 claim",
            "role": "리뷰어가 가장 먼저 공격할 부분을 우리가 먼저 제한한다.",
            "data_used": "claim boundary reports, Samsung reviewer defense table, no-Yu feasibility report.",
            "how_to_read": "왼쪽은 논문/제안서에서 써도 되는 문장, 오른쪽은 쓰면 안 되는 문장이다.",
            "paper_use": "Reviewer response, internal decision memo, Discussion limitations.",
            "narration": "이 프로젝트의 강점은 데이터가 아직 없는 부분을 숨기지 않고 정확히 다음 실험으로 넘긴다는 점이다.",
            "boundary": "defensive slide가 아니라 scientific scope control로 제시해야 한다.",
        },
        {
            "figure": "PF06",
            "image": "PF06_figure_table_roadmap.png",
            "title": "Figure/Table 논문 배치도",
            "role": "각 그림이 논문 어디에 들어가는지 고정해 writing drift를 막는다.",
            "data_used": "figure_table_manifest.tsv.",
            "how_to_read": "각 박스는 paper item, 제목, manuscript position, data source를 같이 보여준다.",
            "paper_use": "논문 outline과 slide deck 제작 control sheet.",
            "narration": "그림을 많이 넣되, 모든 그림이 본문 주장 하나를 책임지게 만든다.",
            "boundary": "장식용 figure나 출처 없는 figure는 사용하지 않는다.",
        },
        {
            "figure": "PF07",
            "image": "PF07_reviewer_to_figure_map.png",
            "title": "리뷰어 공격 대응 지도",
            "role": "예상 공격을 figure/table support와 바로 연결한다.",
            "data_used": "reviewer defense logic, claim boundary table, figure manifest.",
            "how_to_read": "왼쪽은 공격, 가운데는 대응 그림/표, 오른쪽은 짧은 답변이다.",
            "paper_use": "Samsung 발표 Q&A, rebuttal 준비, professor briefing.",
            "narration": "가장 위험한 공격은 'public pilot has no CTC'인데, 답은 '맞다. 그래서 marker-prior evidence로만 쓴다'이다.",
            "boundary": "방어를 위해 새 claim을 만들지 않는다.",
        },
        {
            "figure": "PF08",
            "image": "PF08_cv2_visual_atlas.png",
            "title": "CV2 Visual Atlas",
            "role": "교수님이나 공동연구자에게 전체 흐름을 한 장으로 빠르게 보여주는 overview image.",
            "data_used": "paper-flow figures, no-Yu in silico figures, CV2 contact-sheet rendering.",
            "how_to_read": "위에서 아래로 보면 논리 흐름, 데이터 provenance, marker-NGS, gap, claim boundary, TCGA/spatial/scRNA/proteomics 근거가 순서대로 보인다.",
            "paper_use": "Kakao update, meeting first slide, web hero visual.",
            "narration": "이 한 장은 '우리가 무엇을 이미 할 수 있고, 무엇은 병원 데이터가 와야 하는지'를 압축한다.",
            "boundary": "overview image라서 정량값의 원자료 표를 대체하지 않는다.",
        },
        {
            "figure": "PF09",
            "image": "PF09_manuscript_claim_scorecard.png",
            "title": "In silico manuscript claim scorecard",
            "role": "인실리코 결과 중 실제 Results로 쓸 수 있는 claim과 그 강도를 구분한다.",
            "data_used": "manuscript_claim_evidence_scorecard.tsv; TCGA/spatial/scRNA/proteomics/gap/claim-boundary outputs.",
            "how_to_read": "각 줄은 논문 결과 블록 하나다. 가운데는 쓸 수 있는 claim, 오른쪽은 support figure와 왜 아직 strong discovery가 아닌지를 보여준다.",
            "paper_use": "Results section planning, internal go/no-go, professor briefing.",
            "narration": "이 표는 지금 바로 쓸 Results와 future validation으로 남겨야 할 Results를 분리한다.",
            "boundary": "claim strength가 Moderate인 항목을 discovery claim처럼 쓰지 않는다.",
        },
        {
            "figure": "F04",
            "image": "F04_tcga_driver_variance_boundary.png",
            "title": "TCGA driver variance boundary",
            "role": "driver-only model이 tissue-state axis를 전부 설명하지 못한다는 점을 보여준다.",
            "data_used": "tcga_axis_variance_models.tsv.",
            "how_to_read": "driver grouping이 일부 정보를 주지만 axis variation을 완전히 설명하지 못한다. mutation-only interpretation의 한계를 읽는다.",
            "paper_use": "Results 2 보조 figure 또는 Supplement.",
            "narration": "driver mutation은 clone identity에는 중요하지만 state biology를 단독으로 설명하기에는 부족하다.",
            "boundary": "descriptive variance model이며 causality claim이 아니다.",
        },
        {
            "figure": "PF10",
            "image": "PF10_marker_module_support_heatmap.png",
            "title": "Marker module cross-layer planning support scores (supplementary)",
            "role": "marker module별 prior 구성 유용성을 0/1/2로 손코딩한 manuscript-planning matrix. 분석 산출물이 아니라 기획 점수표다.",
            "data_used": "marker_module_cross_layer_support_matrix.tsv (hand-coded planning scores). TCGA-THCA tissue (no CTC); GSE250521 spatial; scRNA module summary; bulk proteomics proxy; Yu 2024 PUBLISHED summary anchor (no unpublished/raw Yu CTC data).",
            "how_to_read": "0=no current support, 1=indirect/context, 2=directly useful for prior construction. 모두 manuscript-planning judgement이며 measured statistic이 아니다.",
            "paper_use": "Supplementary planning matrix. main figure 위치에 두지 않는다.",
            "narration": "thyroid-lineage module이 가장 강한 public prior를 갖고, EMT module들은 plausible하지만 prospective serial blood validation이 필요하다는 기획 판단을 정리한 표다.",
            "boundary": "support score를 assay validation score처럼 말하지 않는다. TCGA에는 CTC 데이터가 없고, Yu 2024 컬럼은 publish된 summary만 인용한 것이다.",
        },
        {
            "figure": "PF11",
            "image": "PF11_manuscript_blueprint_board.png",
            "title": "Manuscript blueprint board",
            "role": "지금 논문으로 쓸 구조와 쓰면 안 되는 구조를 첫 페이지에서 고정한다.",
            "data_used": "manuscript_blueprint.tsv; submission_strategy_matrix.tsv; current evidence and claim-boundary outputs.",
            "how_to_read": "왼쪽은 논문 identity/title/thesis/boundary, 오른쪽은 Results 1-6 writing blocks다.",
            "paper_use": "웹 첫 검토, 교수님 회의, writing control sheet.",
            "narration": "지금 쓸 논문은 translational framework/prior paper이며, CTC transition discovery는 병원 serial data 이후 단계다.",
            "boundary": "blueprint는 manuscript planning tool이지 새 biological result가 아니다.",
        },
        {
            "figure": "PF12",
            "image": "PF12_future_data_unlock_map.png",
            "title": "Future data unlock map",
            "role": "Yu/병원 데이터가 들어오면 어떤 future discovery figure와 claim이 가능해지는지 정의한다.",
            "data_used": "future_data_unlock_map.tsv; current missing-dataset logic; Samsung prospective study design.",
            "how_to_read": "왼쪽은 필요한 hospital/Yu data, 가운데는 통합 discovery dataset, 오른쪽은 그 데이터가 unlock하는 future figures다.",
            "paper_use": "Discussion future work, professor data request, Samsung execution plan.",
            "narration": "현재 논문은 public-data prior에서 멈추고, serial CTC-EMT transition과 residual-risk association은 이 future data가 들어온 뒤에야 열린다.",
            "boundary": "future unlock claim을 현재 Results로 끌어오지 않는다.",
        },
        {
            "figure": "PF13",
            "image": "PF13_execution_board.png",
            "title": "Execution board",
            "role": "논문을 실제로 조립할 때 필요한 figure layout, reviewer risk, 72시간 실행 계획을 한 장에 묶는다.",
            "data_used": "main_supp_figure_layout.tsv; reviewer_risk_heatmap.tsv; next_72h_manuscript_execution.tsv.",
            "how_to_read": "왼쪽은 main figure 순서, 오른쪽 위는 risk severity, 오른쪽 아래는 다음 72시간 작업이다.",
            "paper_use": "실행 회의, 교수님 보고, 다음 작업 분배.",
            "narration": "이 보드는 논문 조립을 빠르게 하되 claim drift를 막는 작업판이다.",
            "boundary": "execution board는 manuscript management 도구이지 biological result가 아니다.",
        },
        {
            "figure": "PF14",
            "image": "PF14_analysis_control_board.png",
            "title": "Analysis control board",
            "role": "논문 숫자, 통계분석계획, stage-gate kill criteria를 한 장에서 고정한다.",
            "data_used": "key_number_ledger.tsv; statistical_analysis_plan.tsv; stage_gate_kill_criteria.tsv.",
            "how_to_read": "왼쪽은 본문에 쓸 key numbers, 가운데는 현재/미래 분석계획, 오른쪽은 실패 시 downgrade/kill rule이다.",
            "paper_use": "Methods planning, reviewer defense, final fact-check.",
            "narration": "이 보드는 숫자 drift와 current/future claim 혼선을 막기 위한 분석 통제판이다.",
            "boundary": "future prospective analysis rows는 현재 Results가 아니라 planned validation이다.",
        },
        {
            "figure": "F03",
            "image": "F03_tcga_braf_state_split.png",
            "title": "BRAF-mutant PTC is not one state",
            "role": "driver mutation 하나로 biological state를 설명할 수 없다는 핵심 public-data 근거.",
            "data_used": "TCGA-THCA BRAF-mutant subset; tcga_braf_only_label_distribution.tsv.",
            "how_to_read": "BRAF-mutant tumors가 여러 vulnerability/tissue-state label로 분산된다. largest label이 약 33.2%라 나머지 66.8%는 다른 state다.",
            "paper_use": "Results 2: mutation alone is insufficient.",
            "narration": "수술 후 혈액 신호를 mutation만으로 해석하면 state biology를 놓칠 수 있다.",
            "boundary": "BRAF split은 CTC transition이나 recurrence를 증명하지 않는다.",
        },
        {
            "figure": "F05",
            "image": "F05_spatial_tissue_state_organization.png",
            "title": "Spatial tissue states are organized",
            "role": "tissue-state가 무작위 transcript noise가 아니라 공간적으로 조직화된 상태임을 보조한다.",
            "data_used": "GSE250521 spatial transcriptomics; spatial_coherence_statistics.tsv; 16 slides, 57,144 spots.",
            "how_to_read": "same-niche coherence z-score가 16/16 slides에서 2를 넘는다. tissue-state label이 주변 spot과 함께 모이는 경향을 뜻한다.",
            "paper_use": "Results 3: tissue context supports marker/ROI logic.",
            "narration": "공간 조직화는 CTC 후보 신호가 어떤 tissue context에서 나올 수 있는지 묻는 근거다.",
            "boundary": "이 그림은 CTC origin proof가 아니다.",
        },
        {
            "figure": "F06",
            "image": "F06_scrna_marker_context.png",
            "title": "scRNA marker context",
            "role": "CTC-EMT panel이 임의의 gene list가 아니라 cell-state context를 가진다는 근거.",
            "data_used": "scrna_module_celltype_attribution.tsv.",
            "how_to_read": "각 module이 어느 cell type/context에서 상대적으로 강한지 본다. marker 선택의 biological plausibility를 설명한다.",
            "paper_use": "Results 4: marker modules have cell context.",
            "narration": "CTC panel은 epithelial/mesenchymal/thyroid-lineage/survival axis를 동시에 읽도록 구성해야 한다.",
            "boundary": "single-cell tissue context는 circulating CTC validation이 아니다.",
        },
        {
            "figure": "F07",
            "image": "F07_proteomics_dediff_direction.png",
            "title": "Proteomics dedifferentiation direction",
            "role": "RNA 기반 prior가 protein 방향성과 완전히 동떨어진 것이 아님을 보조한다.",
            "data_used": "bulk_proteomics_kthyro_module_contrasts.tsv; bulk_proteomics_kthyro_dediff_trends.tsv.",
            "how_to_read": "dedifferentiation 방향에서 RAI/thyroid-lineage protein loss와 myeloid/TGFB/stress gain을 확인한다.",
            "paper_use": "Results 4 supplement or orthogonal support.",
            "narration": "protein proxy는 tissue-state prior의 방향성을 보강하지만, 혈액 CTC를 직접 읽은 것은 아니다.",
            "boundary": "bulk proteomics를 spatial proteomics나 CTC proteomics처럼 말하지 않는다.",
        },
        {
            "figure": "F08",
            "image": "F08_ctc_emt_ngs_prior_panel.png",
            "title": "CTC-EMT-NGS prior panel",
            "role": "실제 제안서/논문에서 사용할 panel logic을 압축한다.",
            "data_used": "marker prior table, public tissue-state support, thyroid cancer driver logic.",
            "how_to_read": "E, E/M, M, thyroid identity, survival/stress, genetic anchor가 별도 module로 정의된다.",
            "paper_use": "Results 5 and Samsung Aim 1-2 bridge.",
            "narration": "핵심은 CTC count가 아니라 state와 genotype을 같이 읽는 것이다.",
            "boundary": "panel은 research-grade prior이며 clinical assay lock이 아니다.",
        },
        {
            "figure": "F09",
            "image": "F09_publishability_decision.png",
            "title": "Publishability decision",
            "role": "지금 가능한 논문과 불가능한 논문을 최종 판정한다.",
            "data_used": "publishability_decision_table.tsv; no-Yu feasibility report.",
            "how_to_read": "YES_NOW는 framework/prior-map paper, NO_WITHOUT_YU는 transition discovery/residual-risk prediction paper다.",
            "paper_use": "Internal strategy, professor-facing decision memo.",
            "narration": "우리는 지금 논문을 낼 수 있지만, 논문 이름을 정확히 붙여야 한다.",
            "boundary": "strategy figure이며 scientific result figure로 과장하지 않는다.",
        },
    ]


def write_figure_explanation_table() -> pd.DataFrame:
    explain = pd.DataFrame(figure_explanation_rows())
    explain.to_csv(TAB / "figure_explanation_master.tsv", sep="\t", index=False)
    return explain


def caption_bank_rows() -> list[dict[str, str]]:
    return [
        {
            "figure": "Fig. 1",
            "image": "PF01_manuscript_flow.png",
            "short_caption": "Manuscript logic linking postoperative residual-risk uncertainty, published PTC CTC feasibility, public multi-omics prior construction, and the proposed prospective perturbation test.",
            "long_caption": "The paper should open by separating what is already known from what remains untested. Yu 2024 supports the feasibility of serial CTC-EMT measurement around thyroidectomy, while public TCGA, spatial, single-cell, and proteomics resources support a tissue-state and marker-prior framework. The actual thyroidectomy-induced CTC-EMT transition remains a prospective test requiring hospital serial blood and matched genetic data.",
            "results_scaffold": "We first defined the study as a staged translational framework rather than a completed CTC-transition discovery analysis. The logic proceeds from the clinical gap after thyroidectomy to a published CTC feasibility anchor, then to public-data-derived tissue-state priors, and finally to a prospective serial blood and matched NGS experiment.",
            "methods_note": "This figure is a synthesis schematic. It does not present a new statistical test; it maps evidence layers to manuscript sections and claim boundaries.",
            "reviewer_sentence": "This is a Science-track biological perturbation proposal because thyroidectomy is the perturbation and CTC-EMT state transition is the primary biological readout.",
        },
        {
            "figure": "Fig. 2",
            "image": "PF02_data_provenance_map.png",
            "short_caption": "Data provenance map defining how each public, published, planned, and vendor layer is used and bounded.",
            "long_caption": "TCGA-THCA, spatial transcriptomics, scRNA marker context, bulk proteomics proxy, Yu 2024, the future Samsung/Yu cohort, and Inocras/vendor workflows have distinct roles. Public data provide tissue-state, spatial, marker, and protein-direction priors; Yu 2024 provides published CTC feasibility; the future cohort is required to test serial CTC-EMT and genetic mapping.",
            "results_scaffold": "We organized all evidence sources into explicit functional categories. No public dataset is treated as a surrogate for serial CTC measurement. Instead, each public layer contributes a defined prior that informs the design of a future liquid-biopsy perturbation experiment.",
            "methods_note": "Dataset roles were assigned from locally generated pilot outputs and proposal reports. The map is intentionally conservative and includes each layer's non-supported claim.",
            "reviewer_sentence": "The public pilot is not being used as CTC proof; it is being used to define a rational marker and genotype interpretation framework.",
        },
        {
            "figure": "Fig. 3",
            "image": "F03_tcga_braf_state_split.png",
            "short_caption": "BRAF-mutant PTC tumors distribute across multiple public tissue-state labels rather than forming a single biological state.",
            "long_caption": "In the TCGA-THCA BRAF-mutant subset, tumors were distributed across seven tissue-state/vulnerability labels, and the largest label accounted for only approximately one third of BRAF-mutant cases. This supports the central premise that driver mutation alone is insufficient to define the biological state relevant to postoperative blood-signal interpretation.",
            "results_scaffold": "Among BRAF-mutant PTC tumors, mutation status did not collapse tumors into a single state. The largest state represented only 33.2% of the BRAF-mutant subset, leaving substantial state heterogeneity within the same driver-defined group.",
            "methods_note": "Input table: tcga_braf_only_label_distribution.tsv generated from the public pilot extra analyses.",
            "reviewer_sentence": "This figure does not claim CTC biology; it justifies why tissue genotype must be interpreted together with state markers.",
        },
        {
            "figure": "Fig. 4",
            "image": "F04_tcga_driver_variance_boundary.png",
            "short_caption": "Driver-only models explain only part of public tissue-state axis variation.",
            "long_caption": "Driver grouping provides useful information but does not fully account for RAI, immune, HLA, CD8, or related state axes in the public TCGA analysis. This supports a phenotype-plus-genotype strategy rather than a mutation-only interpretation of postoperative blood signals.",
            "results_scaffold": "Driver-only models captured only a fraction of tissue-state variation. These residual state axes motivate the combined CTC-EMT phenotype and matched NGS design, where mutation identifies clone origin and phenotype captures transition state.",
            "methods_note": "Input table: tcga_axis_variance_models.tsv. Models are descriptive and not causal.",
            "reviewer_sentence": "The driver analysis is a boundary argument: mutation helps, but it is not enough for the state-transition question.",
        },
        {
            "figure": "Fig. 5",
            "image": "F05_spatial_tissue_state_organization.png",
            "short_caption": "Spatial thyroid tissue states show non-random organization across public spatial transcriptomic slides.",
            "long_caption": "GSE250521 spatial transcriptomics supports the organization of tissue-state labels in thyroid tissue. Same-niche coherence exceeded z>2 in all 16 slides analyzed, across 57,144 spots. This supports ROI and tissue-context logic for marker selection but does not establish the physical origin of circulating tumor cells.",
            "results_scaffold": "Public spatial data showed that tissue-state labels are spatially structured rather than random spot-level noise. This provides tissue-context support for the marker-prior framework that will later be tested in serial blood.",
            "methods_note": "Input table: spatial_coherence_statistics.tsv; source figure: spatial_coherence_summary.png.",
            "reviewer_sentence": "Spatial coherence is used as tissue-context evidence only, not as evidence of CTC shedding.",
        },
        {
            "figure": "Fig. 6",
            "image": "F06_scrna_marker_context.png",
            "short_caption": "Single-cell marker context supports the biological plausibility of epithelial, mesenchymal, thyroid-lineage, and stress/survival modules.",
            "long_caption": "Single-cell module attribution provides cell-context rationale for selecting CTC-EMT and thyroid-lineage marker modules. The goal is not to validate a blood assay from tissue scRNA data, but to avoid an arbitrary marker panel and to ensure that each module has an interpretable biological context.",
            "results_scaffold": "We used single-cell marker attribution to organize candidate CTC modules into epithelial, hybrid/mesenchymal, thyroid-lineage, and survival/stress compartments. This creates a structured marker prior for the prospective assay.",
            "methods_note": "Input table: scrna_module_celltype_attribution.tsv; source figure: scrna_module_celltype_attribution_heatmap.png.",
            "reviewer_sentence": "The scRNA layer supports marker selection, not direct circulating-cell validation.",
        },
        {
            "figure": "Fig. 7",
            "image": "F07_proteomics_dediff_direction.png",
            "short_caption": "Bulk proteomics proxy provides orthogonal support for dedifferentiation-related directionality.",
            "long_caption": "Bulk proteomics contrasts support loss of thyroid/RAI-related programs and gain of stress, myeloid, or TGF-beta-associated programs along dedifferentiation-related comparisons. This provides protein-level directionality for the tissue-state prior, while remaining distinct from spatial proteomics or CTC proteomics.",
            "results_scaffold": "Orthogonal proteomic summaries were consistent with the direction of dedifferentiation-prior modules. This supports including thyroid-lineage loss and survival/stress gain in the CTC-EMT panel design.",
            "methods_note": "Input tables: bulk_proteomics_kthyro_module_contrasts.tsv and bulk_proteomics_kthyro_dediff_trends.tsv.",
            "reviewer_sentence": "The proteomics layer strengthens biological plausibility but does not substitute for serial blood evidence.",
        },
        {
            "figure": "Fig. 8",
            "image": "PF03_marker_to_ngs_map.png",
            "short_caption": "Phenotypic CTC state modules are connected to tumor-informed genetic anchoring.",
            "long_caption": "The proposed interpretation framework separates phenotypic CTC state from clone identity. Epithelial, hybrid E/M, mesenchymal, thyroid-lineage, and survival/stress markers define candidate cell-state fractions; matched tissue-normal NGS and cfDNA/CTC-enriched NGS are required to determine whether persistent postoperative blood signals are genetically anchored to the resected tumor.",
            "results_scaffold": "We therefore defined the CTC-EMT panel as a two-layer system: phenotype first, genetic anchoring second. This prevents overinterpretation of nonspecific postoperative blood cells as tumor-derived residual-risk biology.",
            "methods_note": "Input sources: ctc_emt_marker_prior_table.tsv, specific aims, Samsung NGS map schematic.",
            "reviewer_sentence": "Matched NGS is the guardrail against overcalling CTC phenotype as tumor-derived biology.",
        },
        {
            "figure": "Fig. 9",
            "image": "F08_ctc_emt_ngs_prior_panel.png",
            "short_caption": "CTC-EMT-NGS prior panel for post-thyroidectomy monitoring in PTC.",
            "long_caption": "The candidate panel includes epithelial, hybrid E/M, mesenchymal, thyroid identity, survival/stress, immune/stress, and genetic-anchor modules. It is intended as a research-grade prior for prospective testing, not as a locked clinical diagnostic assay.",
            "results_scaffold": "The final prior panel integrates tissue-state evidence and thyroid-cancer genetic context into a candidate serial blood assay. The central readout is not CTC count alone, but the relationship among CTC state, thyroid-lineage signal, and tumor-informed genetic persistence.",
            "methods_note": "Input table: ctc_emt_marker_prior_table.tsv.",
            "reviewer_sentence": "This panel is a hypothesis-testing framework, not a validated commercial test.",
        },
        {
            "figure": "Fig. 10",
            "image": "PF04_dataset_gap_bridge.png",
            "short_caption": "The missing integrated serial CTC-EMT-NGS dataset defines the prospective study need.",
            "long_caption": "No local or public dataset currently combines serial pre/post-thyroidectomy blood, CTC EMT phenotype, matched tissue-normal NGS, cfDNA or CTC-enriched NGS, and postoperative follow-up. This gap is not a weakness of the public analysis; it is the direct justification for the Samsung Science proposal.",
            "results_scaffold": "The public-data analysis identifies the exact dataset that must be generated next. The missing elements are serial blood, CTC EMT subtyping, tissue-normal clone anchoring, postoperative cfDNA or CTC-enriched sequencing, and longitudinal clinical follow-up.",
            "methods_note": "Input table: serial_ctc_ngs_dataset_gap_map.tsv and proposal study design.",
            "reviewer_sentence": "The proposal is necessary precisely because the integrated dataset does not exist.",
        },
        {
            "figure": "Fig. 11",
            "image": "PF05_claim_boundary_matrix.png",
            "short_caption": "Allowed and forbidden claims for the current evidence package.",
            "long_caption": "The claim-boundary matrix separates publishable current statements from statements that require hospital serial CTC data or future validation. This is a core reviewer defense element because it shows that the proposal is ambitious without treating preliminary public data as clinical proof.",
            "results_scaffold": "We explicitly separated current evidence from future validation. Allowed claims include public tissue-state heterogeneity and published CTC feasibility; forbidden claims include new postoperative CTC transition, recurrence prediction, universal early-PTC NGS success, and clinical diagnostic readiness.",
            "methods_note": "Input table: claim_boundary_matrix.tsv and claim_boundaries_and_kill_criteria.md.",
            "reviewer_sentence": "The most credible version of this project is the version that states exactly what it cannot yet claim.",
        },
        {
            "figure": "Fig. 12",
            "image": "PF07_reviewer_to_figure_map.png",
            "short_caption": "Reviewer attack map linking each expected criticism to figure/table support and a bounded answer.",
            "long_caption": "Expected reviewer criticisms include category confusion, absence of public CTC data, mutation-only sufficiency, spatial-data relevance, CTC nonspecificity, and excessive breadth. Each criticism is paired with a specific figure/table and a constrained answer.",
            "results_scaffold": "We pre-specified reviewer attacks and matched them to evidence. This table should be used to keep oral and written responses narrow, especially when reviewers challenge the absence of local serial CTC data.",
            "methods_note": "Input sources: reviewer defense table, figure_table_manifest.tsv, claim_boundary_matrix.tsv.",
            "reviewer_sentence": "When attacked for lacking CTC data, answer directly: correct, and that is why public data are used only as prior evidence.",
        },
        {
            "figure": "Fig. 13",
            "image": "PF08_cv2_visual_atlas.png",
            "short_caption": "CV2 visual atlas summarizing the full figure-first manuscript package.",
            "long_caption": "The CV2 atlas is designed for fast human review: it places logic, provenance, marker-to-NGS interpretation, dataset gap, claim boundaries, reviewer defense, and key public-data figures into one large visual board. It is a navigation figure, not a substitute for the underlying tables.",
            "results_scaffold": "For internal and professor-facing communication, we generated a CV2 visual atlas to make the evidence structure inspectable at a glance. This supports rapid review while preserving links to source tables and individual figure explanations.",
            "methods_note": "Generated with cv2 from the deployed figure assets.",
            "reviewer_sentence": "Use this as an overview board, not as a standalone quantitative result.",
        },
        {
            "figure": "Fig. 14",
            "image": "F09_publishability_decision.png",
            "short_caption": "Decision matrix distinguishing what can be published now from what requires Yu/hospital raw data.",
            "long_caption": "A public-data framework or methods-rationale paper is feasible now. A postoperative CTC-EMT transition discovery paper or residual-risk prediction paper is not feasible without patient-level serial CTC, matched NGS, and follow-up data.",
            "results_scaffold": "The immediate manuscript should be framed as a public-data prior-map paper. The higher-impact biological discovery manuscript becomes feasible only after the prospective serial dataset is generated or obtained.",
            "methods_note": "Input table: publishability_decision_table.tsv.",
            "reviewer_sentence": "The current manuscript is intentionally scoped as a prior-map paper; the discovery paper is the next data-dependent stage.",
        },
        {
            "figure": "Fig. 15",
            "image": "PF09_manuscript_claim_scorecard.png",
            "short_caption": "Manuscript claim scorecard translating in silico outputs into writeable Results blocks.",
            "long_caption": "The scorecard separates data-backed public in silico claims from conceptual or future-validation claims. It identifies which results can be written now, which figures support each result, and why specific claims remain bounded.",
            "results_scaffold": "The current evidence package yields six writeable in silico result blocks: mutation/state heterogeneity, spatial tissue-state organization, marker-module rationale, phenotype-to-genotype anchoring, integrated dataset gap definition, and publishability scope.",
            "methods_note": "Input table: manuscript_claim_evidence_scorecard.tsv.",
            "reviewer_sentence": "This scorecard prevents moderate public-data evidence from being overstated as prospective CTC discovery.",
        },
        {
            "figure": "Supplementary planning matrix",
            "image": "PF10_marker_module_support_heatmap.png",
            "short_caption": "Supplementary planning matrix of hand-coded support scores for CTC-EMT marker modules (manuscript-planning judgement, not data-derived).",
            "long_caption": "Marker modules were hand-scored (0/1/2) across TCGA tissue-state context, spatial context, single-cell context, proteomics directionality, and the Yu 2024 published anchor. The scores reflect manuscript-planning judgement about how useful each evidence layer is for prior construction. They are not derived from raw data, are not measured statistics, and are not validation metrics. TCGA-THCA contains tissue (RNA/mutation/methylation/clinical) only and has no CTC, blood, or liquid biopsy data. The Yu 2024 column scores the published feasibility summary; no unpublished or raw Yu CTC data was used.",
            "results_scaffold": "As a supplementary planning matrix, we documented our judgement of how useful each public/published evidence layer was for constructing each marker-module prior. Thyroid-lineage modules carry the strongest public prior judgement, while EMT and survival/stress modules remain hypotheses requiring serial blood validation and matched genetic anchoring.",
            "methods_note": "Input table: marker_module_cross_layer_support_matrix.tsv. Scores are hand-coded manuscript-planning judgements (0/1/2), not measurements: 0=no current support, 1=indirect/context support, 2=directly useful for prior construction. TCGA column = TCGA-THCA tissue-state context (no CTC data). Yu 2024 column = scoring of the published Yu 2024 paper's reported feasibility (62 PTC, 87% CTC detection); no unpublished/raw Yu CTC data was accessed.",
            "reviewer_sentence": "This supplementary matrix is a manuscript-planning score, not a validated assay-performance matrix; readers should not interpret it as a cross-layer measurement.",
        },
        {
            "figure": "Fig. 17",
            "image": "PF11_manuscript_blueprint_board.png",
            "short_caption": "Manuscript blueprint for the current no-Yu-data in silico prior/framework paper.",
            "long_caption": "The blueprint board fixes the safe manuscript identity, recommended title, core thesis, allowed Results blocks, and boundaries. It separates a feasible translational framework paper from a future high-impact CTC-transition discovery paper that requires serial hospital data.",
            "results_scaffold": "We finalized the manuscript as a public-data prior/framework paper rather than a discovery paper. The writeable Results sequence is mutation/state heterogeneity, spatial tissue-state organization, marker-context rationale, marker support grading, phenotype-to-genotype anchoring, and integrated dataset-gap definition.",
            "methods_note": "Input tables: manuscript_blueprint.tsv and submission_strategy_matrix.tsv.",
            "reviewer_sentence": "This board is a scope-control device: it prevents the manuscript from drifting into unsupported CTC transition or recurrence-prediction claims.",
        },
        {
            "figure": "Fig. 18",
            "image": "PF12_future_data_unlock_map.png",
            "short_caption": "Future data unlock map separating current public-data claims from future serial CTC-EMT discovery claims.",
            "long_caption": "The data-unlock map specifies which hospital/Yu data elements are required to move beyond the current public-data prior manuscript. Serial CTC phenotype, CTC marker intensity, matched tissue-normal NGS, serial cfDNA, CTC-enriched NGS, pathology fields, and postoperative follow-up each unlock specific future figures and claims.",
            "results_scaffold": "We explicitly mapped the data elements required for the next-stage discovery manuscript. Serial CTC phenotype data would unlock transition plots, matched tissue-normal NGS would unlock clone anchoring, cfDNA would unlock molecular persistence kinetics, and postoperative follow-up would unlock risk-association analyses.",
            "methods_note": "Input table: future_data_unlock_map.tsv.",
            "reviewer_sentence": "Future data requirements are shown to prevent unsupported future claims from being imported into the current public-data manuscript.",
        },
        {
            "figure": "Fig. 19",
            "image": "PF13_execution_board.png",
            "short_caption": "Execution board for assembling the in silico manuscript package.",
            "long_caption": "The execution board combines the proposed main/supplementary figure layout, reviewer-risk severity, and a 72-hour execution sequence. It is designed to convert the current dossier into a manuscript draft and professor-facing data request without expanding unsupported claims.",
            "results_scaffold": "We converted the in silico dossier into an execution-ready package, defining main and supplementary figures, methods reproducibility checks, reviewer risks, safe sentences, hospital CRF requests, and a 72-hour action plan.",
            "methods_note": "Input tables: main_supp_figure_layout.tsv, reviewer_risk_heatmap.tsv, next_72h_manuscript_execution.tsv.",
            "reviewer_sentence": "Execution planning is separated from biological evidence to keep manuscript operations from becoming unsupported claims.",
        },
        {
            "figure": "Fig. 20",
            "image": "PF14_analysis_control_board.png",
            "short_caption": "Analysis control board locking key numbers, statistical plan, and stage-gate criteria.",
            "long_caption": "The analysis control board consolidates the key numerical ledger, statistical analysis plan, and stage-gate kill criteria. It is designed to prevent numeric drift in manuscript writing and to keep current public-data analyses separate from future prospective validation analyses.",
            "results_scaffold": "We created a final analysis-control layer that records key numbers, current versus future statistical questions, and stage-gate criteria. This converts the dossier from a figure collection into a controlled manuscript-preparation system.",
            "methods_note": "Input tables: key_number_ledger.tsv, statistical_analysis_plan.tsv, stage_gate_kill_criteria.tsv.",
            "reviewer_sentence": "The analysis control board makes clear which analyses were performed now and which require future serial hospital data.",
        },
    ]


def write_caption_bank() -> pd.DataFrame:
    captions = pd.DataFrame(caption_bank_rows())
    captions.to_csv(TAB / "figure_caption_and_results_scaffold.tsv", sep="\t", index=False)
    caption_md = ["# Figure Caption And Results Scaffold\n"]
    for _, row in captions.iterrows():
        caption_md.append(
            f"""## {row['figure']}. {row['short_caption']}

**Image:** `{row['image']}`

**Long caption:** {row['long_caption']}

**Results scaffold:** {row['results_scaffold']}

**Methods/provenance note:** {row['methods_note']}

**Reviewer defense sentence:** {row['reviewer_sentence']}
"""
        )
    (REP / "FIGURE_CAPTION_AND_RESULTS_SCAFFOLD_KR.md").write_text("\n".join(caption_md), encoding="utf-8")
    results_md = ["# Results Narrative Scaffold\n\n이 파일은 최종 author voice가 아니라 논문 구조를 고정하기 위한 scaffolding이다.\n"]
    for _, row in captions.iterrows():
        results_md.append(f"\n## {row['figure']}\n\n{row['results_scaffold']}\n\nBoundary: {row['reviewer_sentence']}\n")
    (REP / "RESULTS_NARRATIVE_SCAFFOLD_KR.md").write_text("\n".join(results_md), encoding="utf-8")
    return captions


def write_manuscript_claim_package(metrics: dict) -> pd.DataFrame:
    rows = [
        {
            "result_block": "Result 1",
            "paper_claim": "PTC CTC-EMT monitoring requires a state-aware framework rather than a mutation-only frame.",
            "in_silico_evidence": f"TCGA BRAF-mutant subset splits across {metrics.get('braf_labels')} labels; largest state {metrics.get('braf_largest_pct'):.1f}%.",
            "figure_support": "Fig. 3; Fig. 4",
            "strength": "Moderate",
            "why_not_strong": "Public tissue/bulk result; not serial blood.",
            "writeable_sentence": "Within BRAF-mutant PTC, public tissue-state heterogeneity argues that driver status alone is insufficient for interpreting postoperative residual-risk biology.",
        },
        {
            "result_block": "Result 2",
            "paper_claim": "Thyroid tissue-state programs are spatially organized and can inform marker/ROI logic.",
            "in_silico_evidence": f"GSE250521 same-niche coherence z>2 in {metrics.get('spatial_z_gt2')}/{metrics.get('spatial_slides')} slides; {metrics.get('spatial_spots')} spots.",
            "figure_support": "Fig. 5",
            "strength": "Moderate",
            "why_not_strong": "Spatial tissue context; not CTC origin tracing.",
            "writeable_sentence": "Spatial coherence supports the existence of organized tissue-state niches that can guide marker selection, while remaining insufficient to infer CTC shedding.",
        },
        {
            "result_block": "Result 3",
            "paper_claim": "CTC-EMT marker modules can be biologically organized using public single-cell and proteomic context.",
            "in_silico_evidence": "scRNA module attribution plus proteomics dedifferentiation-direction support.",
            "figure_support": "Fig. 6; Fig. 7",
            "strength": "Moderate",
            "why_not_strong": "Marker prior only; needs blood assay validation.",
            "writeable_sentence": "Single-cell and proteomic context provide biological rationale for epithelial, hybrid/mesenchymal, thyroid-lineage, and survival/stress modules in a research-grade CTC panel.",
        },
        {
            "result_block": "Result 4",
            "paper_claim": "Phenotypic CTC-state interpretation should be genetically anchored.",
            "in_silico_evidence": "Marker-to-NGS framework and thyroid cancer driver/private-variant logic.",
            "figure_support": "Fig. 8; Fig. 9",
            "strength": "Conceptual but essential",
            "why_not_strong": "Requires prospective tissue-normal/cfDNA/CTC sequencing data.",
            "writeable_sentence": "Because EMT-state markers can be nonspecific, persistent postoperative blood signals should be interpreted only after tumor-informed genetic anchoring.",
        },
        {
            "result_block": "Result 5",
            "paper_claim": "The missing integrated dataset is itself the main translational gap.",
            "in_silico_evidence": "No available integrated resource combines serial PTC blood, CTC EMT phenotype, matched tissue-normal NGS, cfDNA/CTC NGS, and outcome.",
            "figure_support": "Fig. 10; Table 1",
            "strength": "Strong as gap statement",
            "why_not_strong": "Gap statements justify a study; they do not prove biology.",
            "writeable_sentence": "The absence of an integrated serial CTC-EMT-genetic dataset defines the precise prospective resource required to test thyroidectomy-induced state transitions.",
        },
        {
            "result_block": "Result 6",
            "paper_claim": "The current manuscript is publishable as a public-data prior map, not as a transition-discovery paper.",
            "in_silico_evidence": "Claim boundary matrix and publishability decision table.",
            "figure_support": "Fig. 11; Fig. 12; Fig. 14",
            "strength": "Strong as strategy decision",
            "why_not_strong": "Not a biological discovery result.",
            "writeable_sentence": "Current public data can support a state-prior manuscript, whereas postoperative CTC transition and residual-risk prediction require patient-level serial validation.",
        },
    ]
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "manuscript_claim_evidence_scorecard.tsv", sep="\t", index=False)

    abstract = """# Manuscript-Ready Writing Scaffold

## Safe Working Title

Tissue-State and Genetic Prior Mapping for Post-Thyroidectomy CTC-EMT Monitoring in Papillary Thyroid Cancer

## Abstract Scaffold

**Background:** Papillary thyroid cancer is usually curable after thyroidectomy, but postoperative residual-risk biology is inferred indirectly from pathology, thyroglobulin, imaging, and follow-up. Published work indicates that circulating tumor cells and EMT phenotypes can be measured serially around thyroidectomy in PTC, but no public dataset integrates serial CTC-EMT phenotyping with matched tumor-normal NGS, cfDNA/CTC sequencing, tissue context, and postoperative outcome.

**Objective:** To build a conservative public-data-derived tissue-state and genetic prior map for future post-thyroidectomy CTC-EMT monitoring in PTC.

**Methods:** We integrated TCGA-THCA mutation/state analyses, thyroid spatial transcriptomic coherence, single-cell marker context, bulk proteomics directionality, and a published serial CTC feasibility anchor. Each evidence layer was assigned a predefined role and claim boundary.

**Results:** Public TCGA analysis showed that BRAF-mutant PTC is not a single tissue state. Spatial transcriptomics supported organized tissue-state niches across thyroid slides. Single-cell and proteomics summaries supported epithelial, hybrid/mesenchymal, thyroid-lineage, and stress/survival marker modules. These layers support a CTC-EMT-NGS prior panel but do not prove postoperative CTC transition.

**Conclusion:** Public data can support a manuscript describing a tissue-state and genetic-prior framework for post-thyroidectomy CTC-EMT monitoring. A true transition-discovery and residual-risk prediction study requires prospective serial CTC/cfDNA/tissue-normal NGS data.

## Reviewer-Safe Final Message

This is not a claim that public omics prove CTC biology. It is a data-defined, claim-bounded in silico framework showing why serial CTC-EMT and matched genetic mapping after thyroidectomy is the necessary next experiment.
"""
    (REP / "MANUSCRIPT_READY_WRITING_SCAFFOLD_KR.md").write_text(abstract, encoding="utf-8")
    return df


def write_marker_cross_layer_support() -> pd.DataFrame:
    rows = [
        ["Epithelial", "EpCAM/KRT8/18/19/CDH1", 1, 1, 2, 1, 2, "moderate", "CTC capture and epithelial fraction baseline"],
        ["Hybrid E/M", "KRT/EpCAM + VIM/FN1/SNAI2/ZEB1/TWIST1", 1, 1, 2, 1, 2, "moderate", "primary transition-state hypothesis"],
        ["Mesenchymal", "VIM/FN1/CDH2/SNAI2/ZEB1/TWIST1", 1, 1, 2, 1, 2, "moderate", "persistent postoperative EM/M hypothesis"],
        ["Thyroid lineage", "PAX8/TG/TPO/TSHR/SLC5A5", 2, 2, 2, 2, 1, "strong prior", "tumor/thyroid identity and dedifferentiation context"],
        ["Survival/stress", "AXL/MET/CD44/ALDH1A1/SOX9/CA9/VEGFA", 1, 1, 1, 2, 0, "exploratory", "persistence biology; needs prospective confirmation"],
        ["Immune/stress", "CD274/CD47/HLA/B2M/TAP", 2, 2, 1, 1, 0, "context prior", "immune/stress interpretation, not ICI prediction"],
        ["Genetic anchor", "BRAF/TERT/RAS/fusions/private variants", 2, 0, 0, 0, 0, "required future anchor", "clone identity; requires matched tissue-normal/cfDNA/CTC NGS"],
    ]
    df = pd.DataFrame(
        rows,
        columns=[
            "module",
            "markers",
            "TCGA_tissue_state_context_no_CTC",
            "spatial_context_support",
            "scRNA_context_support",
            "proteomics_direction_support",
            "Yu2024_published_anchor_support",
            "evidence_grade",
            "paper_use",
        ],
    )
    support_cols = [
        "TCGA_tissue_state_context_no_CTC",
        "spatial_context_support",
        "scRNA_context_support",
        "proteomics_direction_support",
        "Yu2024_published_anchor_support",
    ]
    df["current_support_score_0_10"] = df[support_cols].sum(axis=1)
    df["score_type"] = "manuscript_planning_judgement_not_data_derived"
    df["prospective_validation_need"] = [
        "assay reproducibility + clone anchor",
        "serial transition test + clone anchor",
        "serial persistence test + clone anchor",
        "lineage-low interpretation + clone anchor",
        "high priority validation needed",
        "context only; avoid therapy claim",
        "must be generated prospectively",
    ]
    df.to_csv(TAB / "marker_module_cross_layer_support_matrix.tsv", sep="\t", index=False)
    (REP / "MARKER_MODULE_CROSS_LAYER_SUPPORT_KR.md").write_text(
        "# Marker Module Cross-Layer Planning Support Scores\n\n"
        "**Score type:** manuscript-planning judgement, **not data-derived**. "
        "Each cell is a hand-coded integer (0/1/2) reflecting whether a layer is useful "
        "for prior construction, not a measured statistic.\n\n"
        "Scores: 0=no current support in that layer; 1=indirect/context support; "
        "2=directly useful support for prior construction. These are manuscript-planning "
        "support scores, **not** validation metrics.\n\n"
        "**Data scope (important reviewer note):**\n"
        "- `TCGA_tissue_state_context_no_CTC`: TCGA-THCA is bulk tumor tissue (RNA/mutation/"
        "methylation/clinical). **TCGA does not contain CTC, blood, or liquid biopsy data.** "
        "This column scores how useful TCGA tissue-state analyses are as a prior for the "
        "marker module, not any measurement of CTCs.\n"
        "- `Yu2024_published_anchor_support`: scores the **published** Yu 2024 paper's "
        "reported feasibility summary (62 PTC, 87% CTC detection, EMT/M-predominant). "
        "**No unpublished or raw Yu CTC data was used.**\n\n"
        + df.to_markdown(index=False)
        + "\n",
        encoding="utf-8",
    )
    return df


def write_manuscript_blueprint() -> tuple[pd.DataFrame, pd.DataFrame]:
    blueprint_rows = [
        ["Paper identity", "Public-data prior/framework paper", "Tissue-state and genetic prior mapping for post-thyroidectomy CTC-EMT monitoring in PTC", "Do not frame as completed serial CTC discovery"],
        ["Central thesis", "Mutation-only interpretation is insufficient; CTC-EMT monitoring needs state + clone anchoring", "TCGA state heterogeneity + marker/NGS framework", "Do not claim clinical utility"],
        ["Result 1", "BRAF-mutant PTC is not one state", "TCGA BRAF split; driver variance boundary", "No CTC/recurrence inference"],
        ["Result 2", "Tissue states are spatially organized", "GSE250521 coherence and spatial context", "No CTC origin proof"],
        ["Result 3", "Marker modules have public cell/protein rationale", "scRNA + proteomics support", "No blood assay validation claim"],
        ["Result 4", "Marker support differs by module", "Cross-layer marker support matrix", "Qualitative planning score only"],
        ["Result 5", "Phenotype must be genetically anchored", "CTC-EMT-to-NGS map", "Do not call marker-positive cells residual disease without clone anchor"],
        ["Result 6", "The missing integrated dataset defines the next experiment", "dataset gap bridge", "Gap statement is not biological proof"],
        ["Discussion", "Public data can justify a prospective Science proposal", "claim boundary + reviewer map", "Avoid broad ICT/platform framing"],
        ["Next dataset", "Yu/hospital raw serial CTC + matched tissue-normal NGS + cfDNA/CTC NGS + follow-up", "prospective validation plan", "No prediction model without follow-up"],
    ]
    blueprint = pd.DataFrame(blueprint_rows, columns=["section", "what_to_write", "support", "boundary"])
    blueprint.to_csv(TAB / "manuscript_blueprint.tsv", sep="\t", index=False)

    strategy_rows = [
        ["Fast preprint / methods-rationale", "Highest immediate feasibility", "Public prior map + figures + claim boundary", "Lower biological discovery weight", "GO now"],
        ["Translational framework paper", "Best balanced option", "Tissue-state prior + CTC-EMT-NGS study design", "Needs careful language; no patient-level CTC", "GO after Yu discussion"],
        ["Samsung preliminary-data appendix", "Grant support value", "All figures/tables directly reusable", "Not a standalone journal goal", "GO now"],
        ["High-impact mechanism paper", "Not feasible now", "Would need serial CTC raw data, matched NGS, outcomes", "Current data insufficient", "NO until Yu/hospital data"],
        ["Clinical prediction paper", "Not feasible now", "Would need recurrence/Tg/US follow-up and validation", "Forbidden current claim", "NO"],
    ]
    strategy = pd.DataFrame(strategy_rows, columns=["paper_route", "feasibility", "what_it_uses", "main_risk", "decision"])
    strategy.to_csv(TAB / "submission_strategy_matrix.tsv", sep="\t", index=False)

    md = ["# Manuscript Blueprint\n"]
    md.append("## Safe Manuscript Identity\n\nPublic-data tissue-state and genetic-prior framework paper for post-thyroidectomy CTC-EMT monitoring in PTC.\n")
    md.append("## Recommended Title\n\nTissue-State and Genetic Prior Mapping for Post-Thyroidectomy CTC-EMT Monitoring in Papillary Thyroid Cancer\n")
    md.append("## Section Blueprint\n\n")
    md.append(blueprint.to_markdown(index=False))
    md.append("\n\n## Submission Strategy\n\n")
    md.append(strategy.to_markdown(index=False))
    md.append("\n\n## Final Boundary\n\nDo not call this a CTC transition discovery paper until Yu/hospital serial CTC and matched genetic data are available.\n")
    (REP / "MANUSCRIPT_BLUEPRINT_KR.md").write_text("\n".join(md), encoding="utf-8")
    return blueprint, strategy


def write_future_data_unlock_map() -> pd.DataFrame:
    rows = [
        ["Serial CTC phenotype table", "patient_id, timepoint, total CTC, E/EM/M count or fraction, assay QC", "Future Fig A: patient-level CTC-EMT transition spaghetti/paired plot", "thyroidectomy perturbs CTC-EMT state", "missing locally", "no transition claim until obtained"],
        ["CTC marker intensity/image features", "cell_id, marker intensities, phenotype call, capture method, batch", "Future Fig B: CTC state classifier/QC panel", "CTC subtyping reproducibility", "missing locally", "phenotype calls need reproducibility metrics"],
        ["Matched tumor-normal NGS", "sample_id, tumor variants, normal filter, VAF, CNV/fusion if available", "Future Fig C: tumor-informed clone map", "blood signal can be genetically anchored", "planned", "NGS success threshold required"],
        ["Serial cfDNA variants", "patient_id, timepoint, variant, VAF, depth, fragment/QC metrics", "Future Fig D: cfDNA persistence kinetics", "postoperative molecular persistence", "missing locally", "early PTC cfDNA may be low/negative"],
        ["CTC-enriched NGS", "CTC-high subset, input cells/DNA, panel, variants, QC pass/fail", "Future Fig E: CTC genotype concordance", "CTC phenotype links to tumor clone", "high-risk feasibility only", "do not promise universal success"],
        ["Pathology risk fields", "LVI, LNM, ETE, tumor size, multifocality, BRAF/TERT/fusion context", "Future Fig F: CTC state vs pathological risk", "association with residual-risk features", "hospital needed", "association not prediction"],
        ["Postoperative follow-up", "Tg, anti-Tg, ultrasound, RAI, recurrence/persistence, follow-up date", "Future Fig G: longitudinal risk association", "residual-risk biology hypothesis", "hospital needed", "requires adequate follow-up and validation"],
        ["Sample/assay logistics", "blood volume, tube, processing time, storage, batch, vendor QC", "Future Table S: assay feasibility and failure modes", "stage-gate feasibility", "planned", "failed QC is a kill criterion, not hidden noise"],
    ]
    df = pd.DataFrame(rows, columns=["data_needed", "minimal_fields", "unlocks_figure", "unlocks_claim", "current_status", "boundary"])
    df.to_csv(TAB / "future_data_unlock_map.tsv", sep="\t", index=False)
    (REP / "FUTURE_DATA_UNLOCK_MAP_KR.md").write_text(
        "# Future Data Unlock Map\n\n"
        "이 표는 Yu/병원 데이터가 들어오면 어떤 figure와 claim이 새로 가능해지는지 정의한다. 현재 public-data manuscript와 future discovery manuscript를 분리하기 위한 control table이다.\n\n"
        + df.to_markdown(index=False)
        + "\n",
        encoding="utf-8",
    )
    return df


def write_execution_package() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    main_figs = pd.DataFrame(
        [
            ["Main Fig 1", "Study logic and data provenance", "PF01 + PF02 combined", "Introduction/overview", "public data as prior, Yu as feasibility, hospital as future test"],
            ["Main Fig 2", "Mutation/state heterogeneity", "F03 + F04", "Results 1", "BRAF/mutation alone is insufficient"],
            ["Main Fig 3", "Spatial tissue-state organization", "F05 + source spatial coherence", "Results 2", "tissue state is organized, not random"],
            ["Main Fig 4", "Marker context and module support", "F06 + F07", "Results 3-4", "marker panel is biologically motivated but not validated; PF10 demoted to supplementary planning matrix"],
            ["Main Fig 5", "CTC-EMT-NGS interpretation framework", "PF03 + F08", "Results 5", "phenotype must be clone-anchored"],
            ["Main Fig 6", "Dataset gap and future unlock", "PF04 + PF12", "Discussion", "current paper stops before discovery claims"],
            ["Main Fig 7", "Claim boundary and reviewer defense", "PF05 + PF07", "Discussion/Reviewer defense", "explicit allowed/forbidden claim boundary"],
            ["Supplementary Fig S1", "CV2 visual atlas", "PF08", "Supplement/meeting", "all visual assets in one board"],
            ["Supplementary Fig S2", "Publishability decision", "F09", "Internal supplement", "paper route decision"],
            ["Supplementary planning matrix", "Marker module cross-layer planning support scores (hand-coded judgement)", "PF10", "Supplement only", "manuscript-planning judgement, not data-derived; TCGA has no CTC; Yu 2024 = published summary only"],
            ["Supplementary Tables", "All data dictionaries and manifests", "TSV tables", "Methods/Supplement", "reproducible provenance"],
        ],
        columns=["paper_slot", "content", "asset", "manuscript_position", "takeaway"],
    )
    main_figs.to_csv(TAB / "main_supp_figure_layout.tsv", sep="\t", index=False)

    methods = pd.DataFrame(
        [
            ["TCGA BRAF split", "tcga_braf_only_label_distribution.tsv", "F03", "BRAF tumors across state labels", "computed", "confirm exact n/percent before manuscript"],
            ["TCGA driver variance", "tcga_axis_variance_models.tsv", "F04", "driver-only model boundary", "computed", "report as descriptive only"],
            ["Spatial coherence", "spatial_coherence_statistics.tsv", "F05", "same-niche z-score across slides", "computed", "cite GSE250521 and spot count"],
            ["scRNA marker context", "scrna_module_celltype_attribution.tsv", "F06", "module-cell context", "computed", "state source dataset in Methods"],
            ["Proteomics proxy", "bulk_proteomics_kthyro_module_contrasts.tsv; dediff_trends.tsv", "F07", "protein-direction support", "computed", "call it bulk proxy, not spatial/CTC proteomics"],
            ["Marker support score (supplementary)", "marker_module_cross_layer_support_matrix.tsv", "PF10", "0/1/2 hand-coded planning judgement (not data-derived; TCGA has no CTC; Yu 2024 column = published summary only)", "new synthesis", "qualitative planning score only; supplementary planning matrix"],
            ["Claim boundary", "claim_boundary_matrix.tsv", "PF05", "allowed vs forbidden statements", "computed", "use in Discussion and reviewer memo"],
            ["Future unlock", "future_data_unlock_map.tsv", "PF12", "hospital data unlocks future figures", "new synthesis", "do not mix with current Results"],
        ],
        columns=["analysis_unit", "input_files", "primary_output", "method_or_metric", "status", "manuscript_qc_note"],
    )
    methods.to_csv(TAB / "methods_reproducibility_checklist.tsv", sep="\t", index=False)

    crf = pd.DataFrame(
        [
            ["patient_core", "patient_id, age, sex, diagnosis date, surgery date", "all patients", "de-identified linkage key", "required"],
            ["tumor_pathology", "tumor size, multifocality, LVI, LNM, ETE, margin, ATA risk", "resected tumor", "pathology report fields", "required"],
            ["blood_timepoints", "T0 pre-op, T1 2 weeks, T2 3 months, T3 6-12 months/recurrence", "serial blood", "exact draw date and tube type", "required"],
            ["CTC_counts", "total CTC, epithelial, hybrid, mesenchymal counts/fractions", "each blood timepoint", "assay method and QC flag", "required"],
            ["CTC_marker_qc", "marker intensity, threshold, batch, image/cell QC if available", "cell-level or sample-level", "needed for reproducibility", "high priority"],
            ["tumor_normal_NGS", "variant list, VAF, coverage, CNV/fusion, normal filter", "tumor-normal pair", "clone anchor", "required"],
            ["cfDNA_NGS", "variant, VAF, depth, fragment/QC, timepoint", "serial plasma", "low signal expected", "high priority"],
            ["followup", "Tg, anti-Tg, ultrasound, RAI, recurrence/persistence, last follow-up", "post-op clinical follow-up", "date-stamped", "required for future association"],
            ["sample_logistics", "blood volume, processing time, storage, shipping, freeze-thaw", "all biospecimens", "failure-mode audit", "required"],
        ],
        columns=["data_block", "fields", "scope", "why_needed", "priority"],
    )
    crf.to_csv(TAB / "hospital_data_request_crf.tsv", sep="\t", index=False)

    risk = pd.DataFrame(
        [
            ["No raw CTC data", 5, "Frame current paper as public prior/framework", "PF02, PF11, PF12", "Do not claim transition discovery"],
            ["Public pilot not CTC biology", 5, "State explicitly that it supports marker/state priors only", "PF05, PF07", "Strong boundary"],
            ["Qualitative support scores look subjective", 4, "Demote PF10 to supplementary; label heatmap as hand-coded planning judgement; rename TCGA column to 'TCGA tissue state context (no CTC)' and Yu 2024 column to 'published anchor support'; state explicitly that no unpublished/raw Yu CTC data was used", "PF10 (supplementary), marker support TSV", "Not validation metrics; not a cross-layer measurement"],
            ["Spatial relevance challenged", 3, "Use only as tissue organization/ROI context", "F05, spatial table", "No CTC origin claim"],
            ["Mutation-alone reviewer", 3, "Show BRAF state split and driver variance boundary", "F03, F04", "Descriptive, not causal"],
            ["Too broad / platform-like", 4, "Keep paper scoped to PTC, thyroidectomy, CTC-EMT, NGS anchor", "PF11", "Remove drug/AI platform language"],
            ["Clinical utility overreach", 5, "Use research-grade and future-validation language", "PF05", "No diagnostic claim"],
            ["Low-input CTC NGS feasibility", 4, "Stage-gate as CTC-high subset feasibility", "PF12, CRF", "No universal success claim"],
        ],
        columns=["reviewer_risk", "severity_1_5", "mitigation", "support", "boundary"],
    )
    risk.to_csv(TAB / "reviewer_risk_heatmap.tsv", sep="\t", index=False)

    sentence_bank = pd.DataFrame(
        [
            ["Mutation/state", "BRAF-mutant PTC remains heterogeneous across public tissue-state labels.", "BRAF status predicts postoperative CTC state."],
            ["Spatial", "Spatial coherence supports organized tissue-state context for marker selection.", "Spatial data prove the origin of CTC shedding."],
            ["scRNA", "Single-cell marker context supports a biologically organized panel prior.", "scRNA validates the circulating CTC phenotype."],
            ["Proteomics", "Bulk proteomics proxy supports directional consistency of selected modules.", "Proteomics proves CTC protein state."],
            ["Marker support", "Cross-layer support differs by module and should guide prospective assay priorities.", "Support score is assay validation."],
            ["NGS anchor", "Phenotypic blood signals require matched genetic anchoring before residual-risk interpretation.", "Marker-positive cells are residual disease."],
            ["Paper identity", "The current manuscript is a public-data prior/framework paper.", "The current manuscript discovers CTC-EMT transition."],
            ["Future data", "Serial hospital CTC/cfDNA/tissue NGS data will be required to test transition and risk association.", "Public data already support recurrence prediction."],
        ],
        columns=["claim_area", "safe_sentence", "forbidden_sentence"],
    )
    sentence_bank.to_csv(TAB / "safe_sentence_bank.tsv", sep="\t", index=False)

    action = pd.DataFrame(
        [
            ["0-12h", "Freeze title, abstract scaffold, and main figure layout", "MANUSCRIPT_BLUEPRINT_KR.md + main_supp_figure_layout.tsv", "ready for Prof. Yu"],
            ["12-24h", "Send hospital/Yu CRF request", "hospital_data_request_crf.tsv", "data request sent"],
            ["24-36h", "Convert figure cards into manuscript outline", "figure_caption_and_results_scaffold.tsv", "Results skeleton ready"],
            ["36-48h", "Audit exact dataset citations and methods", "methods_reproducibility_checklist.tsv", "Methods provenance locked"],
            ["48-60h", "Prepare Samsung appendix deck pages from PF11/PF12 (+ PF10 as supplementary planning matrix only)", "web figures", "proposal appendix ready"],
            ["60-72h", "Decide preprint/framework route vs wait for Yu raw data", "submission_strategy_matrix.tsv", "go/no-go decision"],
        ],
        columns=["window", "task", "input_output", "done_definition"],
    )
    action.to_csv(TAB / "next_72h_manuscript_execution.tsv", sep="\t", index=False)

    md = ["# Execution Package\n"]
    for title, df in [
        ("Main/Supp Figure Layout", main_figs),
        ("Methods Reproducibility Checklist", methods),
        ("Hospital Data Request CRF", crf),
        ("Reviewer Risk Heatmap", risk),
        ("Safe Sentence Bank", sentence_bank),
        ("Next 72h Execution", action),
    ]:
        md.append(f"\n## {title}\n\n{df.to_markdown(index=False)}\n")
    (REP / "MANUSCRIPT_EXECUTION_PACKAGE_KR.md").write_text("\n".join(md), encoding="utf-8")
    return main_figs, methods, crf, risk, sentence_bank, action


def write_analysis_control_package(metrics: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    key_numbers = pd.DataFrame(
        [
            ["TCGA BRAF-mutant subset", metrics.get("braf_n"), "tumors", "tcga_braf_only_label_distribution.tsv", "BRAF-mutant tumors used for mutation/state split", "tissue only; no CTC"],
            ["BRAF tissue-state labels", metrics.get("braf_labels"), "labels", "tcga_braf_only_label_distribution.tsv", "BRAF tumors distribute across multiple states", "label scheme is public-pilot derived"],
            ["Largest BRAF state", f"{metrics.get('braf_largest_pct'):.1f}", "%", "tcga_braf_only_label_distribution.tsv", "largest state is only about one-third", "do not infer recurrence"],
            ["BRAF ambiguity", f"{metrics.get('braf_ambiguity_pct'):.1f}", "%", "tcga_braf_only_label_distribution.tsv", "remaining BRAF cases are outside largest state", "descriptive only"],
            ["Spatial slides", metrics.get("spatial_slides"), "slides", "spatial_coherence_statistics.tsv", "spatial coherence evaluated across slides", "not CTC origin proof"],
            ["Spatial spots", metrics.get("spatial_spots"), "spots", "spatial_spot_vulnerability_scores.tsv", "spot-scale tissue-state context", "not cell-level pathology"],
            ["Spatial z>2 slides", f"{metrics.get('spatial_z_gt2')}/{metrics.get('spatial_slides')}", "slides", "spatial_coherence_statistics.tsv", "same-niche coherence supports organization", "not shedding evidence"],
            ["scRNA module rows", metrics.get("scrna_modules"), "rows", "scrna_module_celltype_attribution.tsv", "marker-context prior", "not blood validation"],
            ["Proteomics contrast rows", metrics.get("proteomics_rows"), "rows", "bulk_proteomics_kthyro_module_contrasts.tsv", "protein-direction support", "bulk proxy only"],
            ["Driver-only median R2", f"{metrics.get('driver_only_r2_median_pct'):.1f}", "%", "tcga_axis_variance_models.tsv", "driver-only model boundary", "descriptive not causal"],
        ],
        columns=["metric", "value", "unit", "source_table", "why_it_matters", "boundary"],
    )
    key_numbers.to_csv(TAB / "key_number_ledger.tsv", sep="\t", index=False)

    sap = pd.DataFrame(
        [
            ["Current in silico RQ1", "Does driver mutation collapse PTC into one state?", "TCGA BRAF label distribution", "fraction by state label; descriptive distribution", "barplot + state split table", "Mutation alone is insufficient", "No CTC or outcome claim"],
            ["Current in silico RQ2", "Are tissue-state labels spatially organized?", "GSE250521 spatial scores", "same-neighbor coherence z-score", "coherence barplot", "Tissue context supports marker/ROI logic", "No CTC origin proof"],
            ["Current in silico RQ3", "Which marker modules have public support?", "TCGA/spatial/scRNA/proteomics/Yu anchor coding", "0/1/2 cross-layer support score", "marker support heatmap", "Prioritizes marker modules", "Not assay validation"],
            ["Current in silico RQ4", "What claims are currently writeable?", "claim/evidence scorecard", "evidence grade + support figure mapping", "claim scorecard", "Framework paper is feasible", "No discovery claim"],
            ["Future prospective RQ1", "Does thyroidectomy perturb CTC E/EM/M state?", "serial CTC phenotype table", "paired change; mixed-effects or paired nonparametric model", "transition plot", "Primary biological endpoint", "Requires raw serial data"],
            ["Future prospective RQ2", "Are persistent postoperative signals genetically anchored?", "tumor-normal NGS + cfDNA/CTC NGS", "variant concordance/VAF kinetics/QC pass", "clone map", "Genetic map of persistent signals", "No universal success assumption"],
            ["Future prospective RQ3", "Do EM/M persistent signals associate with risk features?", "CTC/cfDNA + pathology + Tg/US/follow-up", "association model; exploratory FDR control", "risk association panel", "Hypothesis-generating risk biology", "No clinical prediction without validation"],
        ],
        columns=["analysis_stage", "question", "input_data", "planned_method", "output", "allowed_interpretation", "boundary"],
    )
    sap.to_csv(TAB / "statistical_analysis_plan.tsv", sep="\t", index=False)

    gates = pd.DataFrame(
        [
            ["G1", "Public-data manuscript", "all source tables present and linked", "missing or undocumented input", "do not write unsupported Methods", "currently pass"],
            ["G2", "CTC assay reproducibility", "repeatable CTC detection/subtyping", "CTC detection/subtyping not reproducible", "stop transition claim", "future"],
            ["G3", "Tissue-normal NGS", "success rate >=85%", "success <85%", "clone-anchor claim weakened", "future"],
            ["G4", "cfDNA feasibility", "usable signal in high-risk or CTC-high subset", "unusable cfDNA even in enriched subset", "cfDNA becomes exploratory/optional", "future"],
            ["G5", "CTC-enriched NGS", "QC pass in CTC-high subset", "low-input CTC NGS fails QC", "do not claim CTC genotype map", "future"],
            ["G6", "Follow-up capture", ">=80% key follow-up fields", "follow-up capture <80%", "no outcome/risk association", "future"],
            ["G7", "Claim discipline", "all claims map to evidence table", "unsupported clinical/prediction claim appears", "remove sentence/figure", "currently pass"],
        ],
        columns=["gate", "stage", "pass_rule", "kill_or_downgrade_rule", "action_if_failed", "status"],
    )
    gates.to_csv(TAB / "stage_gate_kill_criteria.tsv", sep="\t", index=False)

    reproducibility = pd.DataFrame(
        [
            ["build_no_yu_in_silico_paper.py", "generates no-Yu figures/tables/reports", "in_silico_paper_no_yu", "source analysis layer", "rerun before manuscript freeze"],
            ["build_no_yu_cv2_storyboard.py", "generates CV2 storyboards", "visual_qc", "visual overview", "optional rerun"],
            ["build_no_yu_in_silico_web.py", "deploys no-Yu dossier", "papers hub", "web evidence page", "deployed"],
            ["build_paper_flow_master_dossier.py", "generates master manuscript control dossier", "paper_flow_master_dossier", "current page", "rerun after edits"],
            ["key_number_ledger.tsv", "locks numbers used in text", "tables", "fact checking", "must match final manuscript"],
            ["statistical_analysis_plan.tsv", "separates current and future analyses", "tables", "reviewer defense", "must not overclaim future analyses"],
            ["hospital_data_request_crf.tsv", "defines hospital/Yu fields", "tables", "data request", "send to collaborator"],
        ],
        columns=["artifact", "purpose", "location", "use_in_manuscript", "qc_action"],
    )
    reproducibility.to_csv(TAB / "reproducibility_manifest.tsv", sep="\t", index=False)

    downgrade = pd.DataFrame(
        [
            ["If Yu raw data arrives", "upgrade to serial CTC transition manuscript", "add patient-level transition and clone map figures", "still avoid clinical prediction until follow-up"],
            ["If only aggregate Yu data arrives", "keep current prior/framework paper", "cite as feasibility anchor only", "no reanalysis claim"],
            ["If no hospital data before deadline", "submit/prepare framework paper and Samsung appendix", "use PF11/PF13 execution route", "make limitations explicit"],
            ["If CTC NGS fails", "retain phenotype-transition biology only", "downgrade genetic-map Aim to cfDNA/tissue anchor", "do not hide failure"],
            ["If follow-up is immature", "publish biological transition/gap paper", "defer recurrence/risk association", "no prediction language"],
        ],
        columns=["scenario", "paper_route", "what_to_do", "boundary"],
    )
    downgrade.to_csv(TAB / "scenario_decision_tree.tsv", sep="\t", index=False)

    md = ["# Analysis Control Package\n"]
    for title, df in [
        ("Key Number Ledger", key_numbers),
        ("Statistical Analysis Plan", sap),
        ("Stage-Gate Kill Criteria", gates),
        ("Reproducibility Manifest", reproducibility),
        ("Scenario Decision Tree", downgrade),
    ]:
        md.append(f"\n## {title}\n\n{df.to_markdown(index=False)}\n")
    (REP / "ANALYSIS_CONTROL_PACKAGE_KR.md").write_text("\n".join(md), encoding="utf-8")
    return key_numbers, sap, gates, reproducibility, downgrade


def plot_master_figures(metrics: dict, data: pd.DataFrame, manifest: pd.DataFrame, flow: pd.DataFrame, claims: pd.DataFrame) -> None:
    fig, ax, path = setup(
        "PF01_manuscript_flow.png",
        "Paper logic flow: from clinical gap to testable serial CTC-EMT-NGS experiment",
        "This is the manuscript skeleton, not just a slide list.",
    )
    xs = [0.06, 0.28, 0.50, 0.72]
    top = [
        ("Clinical gap", "post-op residual-risk biology\nis indirectly monitored", CYAN),
        ("Published anchor", "Yu 2024: serial PTC\nCTC-EMT measurement feasible", GOLD),
        ("In silico prior", "TCGA + spatial + scRNA\n+ proteomics support panel", TEAL),
        ("Prospective test", "hospital serial blood +\nmatched NGS", PURPLE),
    ]
    for i, (t, s, c) in enumerate(top):
        box(ax, xs[i], 0.62, 0.18, 0.16, t, s, c, 13)
        if i:
            arrow(ax, xs[i - 1] + 0.18, 0.70, xs[i], 0.70)
    lows = [
        ("Result 1", "mutation alone is not enough"),
        ("Result 2", "spatial context is organized"),
        ("Result 3", "marker modules have rationale"),
        ("Result 4", "claim boundary is explicit"),
    ]
    for i, (t, s) in enumerate(lows):
        box(ax, 0.10 + i * 0.21, 0.30, 0.17, 0.12, t, s, ORANGE, 12)
        arrow(ax, 0.18 + i * 0.21, 0.42, 0.37, 0.58, ORANGE)
    box(ax, 0.30, 0.09, 0.40, 0.10, "Paper identity", "public-data CTC-EMT-NGS prior map; not a CTC transition discovery paper", WARN, 14)
    save(fig, path)

    fig, ax, path = setup("PF02_data_provenance_map.png", "Data provenance and use boundary", "Every evidence layer has a role and a hard limit.")
    for i, row in data.iterrows():
        x = 0.05 + (i % 2) * 0.46
        y = 0.70 - (i // 2) * 0.17
        color = [CYAN, PURPLE, TEAL, ORANGE, GOLD, WARN, GREEN][i % 7]
        box(ax, x, y, 0.40, 0.125, row["dataset_id"], f"{row['used_for']}\nLimit: {row['forbidden_claim']}", color, 10.8, align="left")
    save(fig, path)

    fig, ax, path = setup("PF03_marker_to_ngs_map.png", "Marker-to-NGS map", "Phenotype alone is weak; clone anchoring is the scientific guardrail.")
    modules = [
        ("Epithelial", "EpCAM/KRT/CDH1", CYAN),
        ("Hybrid E/M", "KRT + VIM/FN1/SNAI2", GOLD),
        ("Mesenchymal", "VIM/FN1/CDH2/ZEB1", ORANGE),
        ("Thyroid lineage", "PAX8/TG/TPO/TSHR/SLC5A5", TEAL),
        ("Survival/stress", "AXL/MET/CD44/CD274/CD47", PURPLE),
    ]
    for i, (t, s, c) in enumerate(modules):
        box(ax, 0.055, 0.70 - i * 0.12, 0.27, 0.085, t, s, c, 12, align="left")
        arrow(ax, 0.33, 0.742 - i * 0.12, 0.47, 0.50, c)
    box(ax, 0.46, 0.42, 0.22, 0.15, "CTC fraction state", "E / E-M / M\nresearch phenotype", GOLD, 14)
    arrow(ax, 0.68, 0.50, 0.76, 0.50)
    box(ax, 0.76, 0.42, 0.18, 0.15, "Tumor-informed NGS", "BRAF/TERT/RAS/fusions\nprivate variants", GREEN, 12)
    box(ax, 0.39, 0.15, 0.37, 0.12, "Allowed interpretation", "persistent post-op signal must be matched to tissue/cfDNA before residual-risk biology is claimed", WARN, 11)
    save(fig, path)

    fig, ax, path = setup("PF04_dataset_gap_bridge.png", "The missing dataset is the proposal justification", "Public data builds the prior; hospital data tests the perturbation.")
    pieces = [
        ("Serial blood", "T0/T1/T2/T3", WARN),
        ("CTC EMT", "E / E-M / M", WARN),
        ("Tissue-normal NGS", "clone anchor", WARN),
        ("cfDNA/CTC NGS", "blood genetic map", WARN),
        ("Outcome follow-up", "Tg/US/recurrence-risk", WARN),
    ]
    for i, (t, s, c) in enumerate(pieces):
        box(ax, 0.07 + i * 0.18, 0.58, 0.14, 0.12, t, s, c, 11)
        ax.text(0.14 + i * 0.18, 0.49, "MISSING\nin public", color=WARN, fontsize=10, fontweight="bold", ha="center", transform=ax.transAxes)
    box(ax, 0.13, 0.22, 0.31, 0.13, "What public data can do", "marker/state/ROI/genotype prior", TEAL, 13)
    box(ax, 0.56, 0.22, 0.31, 0.13, "What Samsung cohort must do", "test post-thyroidectomy transition", GOLD, 13)
    arrow(ax, 0.44, 0.285, 0.56, 0.285, GOLD)
    save(fig, path)

    fig, ax, path = setup("PF05_claim_boundary_matrix.png", "Claim boundary matrix", "The strongest review defense is showing what we refuse to overclaim.")
    allowed = claims[claims["status"].eq("Allowed")].head(5)
    forbidden = claims[claims["status"].eq("Forbidden")].head(5)
    ax.text(0.07, 0.78, "Allowed", color=TEAL, fontsize=18, fontweight="bold", transform=ax.transAxes)
    ax.text(0.55, 0.78, "Forbidden", color=WARN, fontsize=18, fontweight="bold", transform=ax.transAxes)
    for i, (_, row) in enumerate(allowed.iterrows()):
        box(ax, 0.06, 0.68 - i * 0.105, 0.39, 0.075, row["statement"], row["boundary"], TEAL, 8.6, align="left")
    for i, (_, row) in enumerate(forbidden.iterrows()):
        box(ax, 0.54, 0.68 - i * 0.105, 0.39, 0.075, row["statement"], row["boundary"], WARN, 8.6, align="left")
    save(fig, path)

    fig, ax, path = setup("PF06_figure_table_roadmap.png", "Figure and table roadmap", "Where each item belongs in the manuscript.")
    slots = manifest.head(8).copy()
    for i, (_, row) in enumerate(slots.iterrows()):
        x = 0.06 + (i % 2) * 0.46
        y = 0.73 - (i // 2) * 0.15
        color = [CYAN, GOLD, TEAL, PURPLE, ORANGE, GREEN, WARN, CYAN][i]
        box(ax, x, y, 0.39, 0.10, f"{row['paper_item']}: {row['title']}", f"{row['manuscript_position']} | {row['data_used']}", color, 9.8, align="left")
    save(fig, path)

    fig, ax, path = setup("PF07_reviewer_to_figure_map.png", "Reviewer attack to figure/table support", "Make the defense visible before the reviewer asks.")
    rows = [
        ("This is not Science; it is ICT.", "Fig. 1 / Fig. 7", "one perturbation: thyroidectomy"),
        ("Public pilot has no CTC.", "Fig. 2 / Table 1", "correct; it is marker-prior evidence"),
        ("Mutation/driver explains enough.", "Fig. 3", "BRAF state split argues otherwise"),
        ("Spatial data are irrelevant.", "Fig. 4", "ROI and tissue-state context only"),
        ("CTC phenotype is nonspecific.", "Fig. 6", "matched NGS is the guardrail"),
        ("Too broad.", "Table 3 / Fig. 8", "remove drug/AI/spatial-platform claims"),
    ]
    for i, (attack, support, answer) in enumerate(rows):
        y = 0.72 - i * 0.105
        box(ax, 0.05, y, 0.31, 0.072, attack, "", WARN, 8.5, align="left")
        box(ax, 0.40, y, 0.18, 0.072, support, "", GOLD, 9.2)
        box(ax, 0.62, y, 0.31, 0.072, answer, "", TEAL, 8.5, align="left")
        arrow(ax, 0.36, y + 0.036, 0.40, y + 0.036, GOLD)
        arrow(ax, 0.58, y + 0.036, 0.62, y + 0.036, GOLD)
    save(fig, path)


def plot_manuscript_claim_figure(scorecard: pd.DataFrame) -> None:
    fig, ax, path = setup(
        "PF09_manuscript_claim_scorecard.png",
        "In silico manuscript claim scorecard",
        "What can be written as Results, and how strong each claim is.",
    )
    colors = {
        "Moderate": CYAN,
        "Conceptual but essential": GOLD,
        "Strong as gap statement": TEAL,
        "Strong as strategy decision": PURPLE,
    }
    for i, (_, row) in enumerate(scorecard.iterrows()):
        y = 0.72 - i * 0.105
        c = colors.get(row["strength"], ORANGE)
        box(ax, 0.05, y, 0.13, 0.073, row["result_block"], row["strength"], c, 8.5)
        box(ax, 0.21, y, 0.37, 0.073, row["paper_claim"], "", c, 8.0, align="left")
        box(ax, 0.61, y, 0.31, 0.073, row["figure_support"], row["why_not_strong"], WARN if "not" in row["why_not_strong"].lower() else TEAL, 8.0, align="left")
    box(ax, 0.18, 0.08, 0.64, 0.10, "Manuscript rule", "Write the first five blocks as in silico Results. Keep CTC transition and recurrence prediction as prospective validation, not current findings.", WARN, 12)
    save(fig, path)


def plot_marker_support_heatmap(marker_support: pd.DataFrame) -> None:
    cols = [
        "TCGA_tissue_state_context_no_CTC",
        "spatial_context_support",
        "scRNA_context_support",
        "proteomics_direction_support",
        "Yu2024_published_anchor_support",
    ]
    labels = [
        "TCGA tissue\nstate (no CTC)",
        "Spatial\ncontext",
        "scRNA\ncontext",
        "Proteomics\ndirection",
        "Yu 2024\npublished anchor",
    ]
    vals = marker_support[cols].to_numpy(dtype=float)
    fig, ax, path = setup(
        "PF10_marker_module_support_heatmap.png",
        "Marker module cross-layer planning support scores",
        "Manuscript-planning judgement (0/1/2), NOT data-derived. TCGA = tissue (no CTC); Yu 2024 column uses published summary only.",
    )
    heat_ax = fig.add_axes([0.10, 0.20, 0.62, 0.56], facecolor=BG)
    cmap = plt.cm.get_cmap("viridis", 3)
    im = heat_ax.imshow(vals, vmin=0, vmax=2, cmap=cmap, aspect="auto")
    heat_ax.set_xticks(range(len(cols)))
    heat_ax.set_xticklabels(labels, color=TEXT, fontsize=10)
    heat_ax.set_yticks(range(marker_support.shape[0]))
    heat_ax.set_yticklabels(marker_support["module"], color=TEXT, fontsize=10)
    heat_ax.tick_params(length=0)
    for spine in heat_ax.spines.values():
        spine.set_color(LINE)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            heat_ax.text(j, i, int(vals[i, j]), ha="center", va="center", color="#07101a" if vals[i, j] > 1 else TEXT, fontsize=11, fontweight="bold")
    cax = fig.add_axes([0.75, 0.22, 0.03, 0.52])
    cb = fig.colorbar(im, cax=cax, ticks=[0, 1, 2])
    cb.ax.tick_params(colors=TEXT, labelsize=9)
    cb.outline.set_edgecolor(LINE)
    ax.text(0.81, 0.65, "Interpretation", color=GOLD, fontsize=16, fontweight="bold", transform=ax.transAxes)
    ax.text(
        0.81,
        0.53,
        "Thyroid-lineage has the strongest\npublic prior. EMT modules are\nplausible but require serial blood\nvalidation and clone anchoring.",
        color=TEXT,
        fontsize=12,
        transform=ax.transAxes,
    )
    ax.text(0.81, 0.35, "Boundary", color=WARN, fontsize=16, fontweight="bold", transform=ax.transAxes)
    ax.text(
        0.81,
        0.20,
        "Manuscript-planning judgement,\nnot data-derived. TCGA has no\nCTC data; Yu 2024 column uses\npublished anchor summary only,\nnot unpublished/raw Yu data.",
        color=TEXT,
        fontsize=11,
        transform=ax.transAxes,
    )
    ax.text(
        0.10,
        0.10,
        "Supplementary planning matrix · not a validation metric · do not cite as cross-layer measurement.",
        color=WARN,
        fontsize=10,
        fontweight="bold",
        transform=ax.transAxes,
    )
    save(fig, path)


def plot_manuscript_blueprint_board(blueprint: pd.DataFrame, strategy: pd.DataFrame) -> None:
    fig, ax, path = setup(
        "PF11_manuscript_blueprint_board.png",
        "Manuscript blueprint: what to write now",
        "A write-safe structure for a no-Yu-data in silico prior/framework paper.",
    )
    left = [
        ("Identity", "public-data prior/framework\nnot CTC discovery", GOLD),
        ("Title", "Tissue-State and Genetic Prior Mapping\nfor Post-Thyroidectomy CTC-EMT Monitoring", CYAN),
        ("Core thesis", "state + clone anchor\nbeats mutation-only interpretation", TEAL),
        ("Boundary", "no recurrence prediction\nno clinical diagnostic", WARN),
    ]
    for i, (t, s, c) in enumerate(left):
        box(ax, 0.05, 0.72 - i * 0.15, 0.34, 0.105, t, s, c, 11, align="left")

    result_rows = blueprint[blueprint["section"].str.startswith("Result")].copy()
    for i, (_, row) in enumerate(result_rows.iterrows()):
        x = 0.46 + (i % 2) * 0.23
        y = 0.72 - (i // 2) * 0.15
        box(ax, x, y, 0.19, 0.105, row["section"], row["what_to_write"], [CYAN, PURPLE, TEAL, ORANGE, GREEN, GOLD][i], 9, align="left")
    box(ax, 0.46, 0.22, 0.42, 0.12, "Best route", "Translational framework paper now; mechanism/discovery paper only after serial hospital data", GOLD, 12)
    ax.text(0.05, 0.08, "Use this board as the first-page writing control: if a sentence does not fit one of these boxes, it probably belongs in future validation.", color=MUTED, fontsize=11, transform=ax.transAxes)
    save(fig, path)


def plot_future_data_unlock_map(unlock: pd.DataFrame) -> None:
    fig, ax, path = setup(
        "PF12_future_data_unlock_map.png",
        "Future data unlock map",
        "What Yu/hospital data unlocks beyond the current public-data manuscript.",
    )
    inputs = [
        ("Serial CTC\nphenotype", "transition plot", CYAN),
        ("CTC marker\nintensity", "subtyping QC", TEAL),
        ("Tumor-normal\nNGS", "clone map", GOLD),
        ("Serial\ncfDNA", "VAF kinetics", PURPLE),
        ("CTC-enriched\nNGS", "genotype concordance", ORANGE),
        ("Pathology +\nfollow-up", "risk association", WARN),
    ]
    for i, (t, s, c) in enumerate(inputs):
        y = 0.73 - i * 0.095
        box(ax, 0.06, y, 0.22, 0.065, t, s, c, 8.5)
        arrow(ax, 0.28, y + 0.032, 0.43, 0.50, c)
    box(ax, 0.42, 0.42, 0.20, 0.16, "Future discovery dataset", "serial CTC-EMT +\nmatched NGS + cfDNA\n+ outcome", GOLD, 12)
    outputs = [
        ("Future Fig A", "CTC-EMT transition"),
        ("Future Fig C", "tumor-informed clone map"),
        ("Future Fig D", "cfDNA persistence"),
        ("Future Fig F/G", "risk association"),
    ]
    for i, (t, s) in enumerate(outputs):
        y = 0.70 - i * 0.13
        box(ax, 0.70, y, 0.22, 0.075, t, s, TEAL if i < 2 else PURPLE, 9)
        arrow(ax, 0.62, 0.50, 0.70, y + 0.038, GOLD)
    box(ax, 0.28, 0.12, 0.46, 0.10, "Boundary", "Current public-data paper stops before these future figures. Do not import these future claims into the current manuscript.", WARN, 11)
    save(fig, path)


def plot_execution_board(main_figs: pd.DataFrame, risk: pd.DataFrame, action: pd.DataFrame) -> None:
    fig, ax, path = setup(
        "PF13_execution_board.png",
        "Execution board: assemble the paper without drifting claims",
        "Main figure layout, reviewer risks, and next 72-hour actions.",
    )
    for i, (_, row) in enumerate(main_figs.head(7).iterrows()):
        y = 0.72 - i * 0.075
        box(ax, 0.045, y, 0.30, 0.052, row["paper_slot"], row["content"], CYAN if i < 6 else GOLD, 7.5, align="left")
    ax.text(0.41, 0.77, "Reviewer risk severity", color=GOLD, fontsize=15, fontweight="bold", transform=ax.transAxes)
    for i, (_, row) in enumerate(risk.head(8).iterrows()):
        y = 0.70 - i * 0.06
        sev = int(row["severity_1_5"])
        ax.add_patch(FancyBboxPatch((0.41, y), 0.06 * sev, 0.035, boxstyle="round,pad=0.004,rounding_size=0.006", facecolor=WARN if sev >= 4 else ORANGE, edgecolor="none", transform=ax.transAxes))
        ax.text(0.72, y + 0.017, row["reviewer_risk"], color=TEXT, fontsize=8.2, va="center", transform=ax.transAxes)
    ax.text(0.41, 0.17, "72-hour execution", color=TEAL, fontsize=15, fontweight="bold", transform=ax.transAxes)
    for i, (_, row) in enumerate(action.head(3).iterrows()):
        box(ax, 0.41 + i * 0.18, 0.06, 0.15, 0.075, row["window"], row["task"], TEAL, 7.5)
    save(fig, path)


def plot_analysis_control_board(key_numbers: pd.DataFrame, sap: pd.DataFrame, gates: pd.DataFrame) -> None:
    fig, ax, path = setup(
        "PF14_analysis_control_board.png",
        "Analysis control board: numbers, statistics, stage gates",
        "Prevents numeric drift and separates current in silico analyses from future prospective tests.",
    )
    key = key_numbers.head(7)
    for i, (_, row) in enumerate(key.iterrows()):
        y = 0.73 - i * 0.075
        box(ax, 0.045, y, 0.25, 0.052, row["metric"], f"{row['value']} {row['unit']}", GOLD if i < 4 else CYAN, 8, align="left")
    ax.text(0.36, 0.78, "Statistical analysis plan", color=TEAL, fontsize=15, fontweight="bold", transform=ax.transAxes)
    for i, (_, row) in enumerate(sap.head(7).iterrows()):
        y = 0.70 - i * 0.065
        c = TEAL if str(row["analysis_stage"]).startswith("Current") else PURPLE
        box(ax, 0.36, y, 0.27, 0.047, row["analysis_stage"], row["planned_method"], c, 7, align="left")
    ax.text(0.69, 0.78, "Stage gates", color=WARN, fontsize=15, fontweight="bold", transform=ax.transAxes)
    for i, (_, row) in enumerate(gates.head(7).iterrows()):
        y = 0.70 - i * 0.065
        box(ax, 0.69, y, 0.24, 0.047, row["gate"], row["kill_or_downgrade_rule"], WARN if row["status"] == "future" else TEAL, 7, align="left")
    box(ax, 0.24, 0.06, 0.52, 0.085, "Rule", "Current manuscript uses only current in silico RQs. Future prospective RQs are written as planned validation, not current Results.", WARN, 11)
    save(fig, path)


def image_card(img: str, title: str, caption: str, cls: str = "") -> str:
    return f"""
    <article class="fig-card {esc(cls)}">
      <div class="fig-head"><b>{esc(title)}</b><span>{esc(img)}</span></div>
      <a href="assets/{ASSET}/{esc(img)}"><img src="assets/{ASSET}/{esc(img)}" alt="{esc(title)}" loading="lazy" /></a>
      <p>{esc(caption)}</p>
    </article>
    """


def link_card(label: str, fn: str, note: str = "") -> str:
    return f"<a class='link-card' href='assets/{ASSET}/{esc(fn)}'><b>{esc(label)}</b><span>{esc(note or fn)}</span></a>"


def explanation_blocks(explain: pd.DataFrame, captions: pd.DataFrame) -> str:
    caption_by_image = {str(row["image"]): row for _, row in captions.iterrows()}
    out = ["<div class='explain-list'>"]
    for _, row in explain.iterrows():
        cap = caption_by_image.get(str(row["image"]))
        caption_html = ""
        if cap is not None:
            caption_html = f"""
                <div class="embedded-caption">
                  <h4>Caption And Results Scaffold</h4>
                  <div class="caption-grid embedded">
                    <div><h4>Figure Caption</h4><p>{esc(cap['long_caption'])}</p></div>
                    <div><h4>Results Paragraph</h4><p>{esc(cap['results_scaffold'])}</p></div>
                    <div><h4>Methods / Data Provenance</h4><p>{esc(cap['methods_note'])}</p></div>
                    <div><h4>Reviewer-Safe Sentence</h4><p>{esc(cap['reviewer_sentence'])}</p></div>
                  </div>
                </div>
            """
        else:
            caption_html = """
                <div class="embedded-caption muted-cap">
                  <h4>Caption And Results Scaffold</h4>
                  <p>This figure is a navigation/control figure. Use the role, data, reading guide, and boundary above; no separate manuscript caption scaffold was assigned.</p>
                </div>
            """
        out.append(
            f"""
            <article class="explain-card">
              <div class="explain-img">
                <a href="assets/{ASSET}/{esc(row['image'])}"><img src="assets/{ASSET}/{esc(row['image'])}" alt="{esc(row['title'])}" loading="lazy" /></a>
              </div>
              <div class="explain-body">
                <div class="tagline"><span>{esc(row['figure'])}</span><b>{esc(row['title'])}</b></div>
                <div class="mini-grid">
                  <div><h4>Role</h4><p>{esc(row['role'])}</p></div>
                  <div><h4>Data Used</h4><p>{esc(row['data_used'])}</p></div>
                  <div><h4>How To Read</h4><p>{esc(row['how_to_read'])}</p></div>
                  <div><h4>Paper Use</h4><p>{esc(row['paper_use'])}</p></div>
                </div>
                <div class="script"><b>설명 문장</b><p>{esc(row['narration'])}</p></div>
                <div class="boundary"><b>Claim Boundary</b><p>{esc(row['boundary'])}</p></div>
                {caption_html}
              </div>
            </article>
            """
        )
    out.append("</div>")
    return "".join(out)


def caption_blocks(captions: pd.DataFrame) -> str:
    out = ["<div class='caption-list'>"]
    for _, row in captions.iterrows():
        out.append(
            f"""
            <article class="caption-card">
              <div class="caption-title"><span>{esc(row['figure'])}</span><b>{esc(row['short_caption'])}</b></div>
              <div class="caption-grid">
                <div><h4>Long Caption</h4><p>{esc(row['long_caption'])}</p></div>
                <div><h4>Results Scaffold</h4><p>{esc(row['results_scaffold'])}</p></div>
                <div><h4>Methods / Provenance</h4><p>{esc(row['methods_note'])}</p></div>
                <div><h4>Reviewer-Safe Sentence</h4><p>{esc(row['reviewer_sentence'])}</p></div>
              </div>
            </article>
            """
        )
    out.append("</div>")
    return "".join(out)


def write_reports(
    data: pd.DataFrame,
    manifest: pd.DataFrame,
    flow: pd.DataFrame,
    claims: pd.DataFrame,
    explain: pd.DataFrame,
    captions: pd.DataFrame,
    scorecard: pd.DataFrame,
    metrics: dict,
) -> None:
    (REP / "MANUSCRIPT_FLOW_MASTER_KR.md").write_text(
        f"""# CTC-EMT-NGS In Silico Paper Flow Master

## 한 줄 정리

지금 논문은 `Serial CTC-EMT transition discovery`가 아니라 **public-data 기반 PTC CTC-EMT-NGS prior map / framework paper**로 써야 한다.

## 논문 제목 후보

**Tissue-State and Genetic Prior Mapping for Post-Thyroidectomy CTC-EMT Monitoring in Papillary Thyroid Cancer**

## 핵심 흐름

1. PTC 수술 후 residual-risk biology는 직접 읽기 어렵다.
2. Yu 2024는 serial CTC-EMT 측정 가능성을 보여준 published anchor다.
3. 하지만 public dataset에는 serial CTC-EMT + matched NGS + cfDNA + tissue context + outcome이 없다.
4. 그러므로 공개 multi-omics로 tissue-state / marker / NGS prior map을 먼저 세운다.
5. 이 prior map은 Samsung/Yu prospective cohort에서 thyroidectomy perturbation으로 테스트한다.

## 주요 숫자

- TCGA BRAF-mutant subset: {metrics.get('braf_n')} tumors, {metrics.get('braf_labels')} tissue-state labels, largest state {metrics.get('braf_largest_pct'):.1f}%.
- Spatial GSE250521: {metrics.get('spatial_slides')} slides, {metrics.get('spatial_spots')} spots, same-niche z>2 in {metrics.get('spatial_z_gt2')}/{metrics.get('spatial_slides')} slides.
- Yu 2024 anchor: 62 prospective PTC patients, pre-op / 2 weeks / 3 months, 87% CTC detection reported.

## 논문에서 금지할 문장

- 우리 데이터가 thyroidectomy 후 CTC-EMT transition을 증명했다.
- public pilot이 CTC shedding을 증명한다.
- early PTC에서 CTC/cfDNA NGS가 universal하게 성공한다.
- 이 결과가 ready clinical diagnostic이다.
""",
        encoding="utf-8",
    )
    (REP / "DATA_DEFINITION_MASTER_KR.md").write_text(
        "# Data Definition Master\n\n"
        + "각 데이터 레이어의 역할과 금지 claim은 아래 표가 기준이다.\n\n"
        + data.to_markdown(index=False)
        + "\n",
        encoding="utf-8",
    )
    (REP / "FIGURE_TABLE_EXPLANATION_GUIDE_KR.md").write_text(
        "# Figure / Table Explanation Guide\n\n"
        + "논문 본문에서 각 그림과 표를 어디에 쓰고 무엇을 말할지 정리한 manifest다.\n\n"
        + manifest.to_markdown(index=False)
        + "\n\n## Manuscript Flow\n\n"
        + flow.to_markdown(index=False)
        + "\n\n## Claim Boundary\n\n"
        + claims.to_markdown(index=False)
        + "\n\n## Figure-by-Figure Explanation Master\n\n"
        + explain.to_markdown(index=False)
        + "\n\n## Caption And Results Scaffold\n\n"
        + captions.to_markdown(index=False)
        + "\n\n## Manuscript Claim Evidence Scorecard\n\n"
        + scorecard.to_markdown(index=False)
        + "\n",
        encoding="utf-8",
    )


def copy_source_assets() -> list[Path]:
    files: list[Path] = []
    files += sorted(FIG.glob("*.png"))
    files += sorted(NOYU.glob("figures/*.png"))
    files += sorted(NOYU.glob("visual_qc/*.png"))
    files += sorted(NOYU.glob("source_figures/*.png"))
    files += sorted(SAMSUNG.glob("figures/slide*.png"))
    files += sorted(SAMSUNG.glob("figures/webvisual_*.png"))
    files += sorted(TAB.glob("*.tsv"))
    files += sorted(REP.glob("*.md"))
    files += sorted(NOYU.glob("tables/*.tsv"))
    files += sorted(NOYU.glob("reports/*.md"))
    extras = [
        SAMSUNG / "deck/kthyro_ctc_emt_samsung_science_pitch_v2_reference_style.pptx",
        PILOT / "ppt_evidence_package/docs/FIGURE_TABLE_GUIDE.md",
        PILOT / "ppt_evidence_package/docs/ANALYSIS_FLOW_COMPREHENSIVE.md",
        PILOT / "extra_analyses/figures/tcga_braf_only_vulnerability_axis_heatmap.png",
        PILOT / "extra_analyses/figures/tcga_axis_variance_models.png",
        PILOT / "figures/spatial_coherence_summary.png",
        PILOT / "extra_analyses/figures/scrna_module_celltype_attribution_heatmap.png",
        PILOT / "extra_analyses/figures/bulk_proteomics_kthyro_axis_boxplots.png",
        PILOT / "extra_analyses/figures/bulk_proteomics_kthyro_contrast_heatmap.png",
        PILOT / "extra_analyses/figures/spatial_niche_adjacency_enrichment_heatmap.png",
        PILOT / "extra_analyses/figures/spatial_hotspot_colocalization_heatmap.png",
        PILOT / "extra_analyses/figures/spatial_interface_distance_z_by_condition.png",
    ]
    files += [p for p in extras if p.exists()]

    copied: list[Path] = []
    seen: set[str] = set()
    for p in files:
        if not p.exists() or p.name in seen:
            continue
        seen.add(p.name)
        target = SRC / p.name
        shutil.copy2(p, target)
        shutil.copy2(p, HUB_ASSET / p.name)
        shutil.copy2(p, LIVE_ASSET / p.name)
        copied.append(target)
    return copied


def build_cv2_atlas() -> None:
    imgs = [
        "PF01_manuscript_flow.png",
        "PF02_data_provenance_map.png",
        "PF03_marker_to_ngs_map.png",
        "PF04_dataset_gap_bridge.png",
        "PF05_claim_boundary_matrix.png",
        "PF06_figure_table_roadmap.png",
        "PF07_reviewer_to_figure_map.png",
        "PF09_manuscript_claim_scorecard.png",
        "PF10_marker_module_support_heatmap.png",
        "PF11_manuscript_blueprint_board.png",
        "PF12_future_data_unlock_map.png",
        "PF13_execution_board.png",
        "PF14_analysis_control_board.png",
        "F01_no_yu_claim_ladder.png",
        "F03_tcga_braf_state_split.png",
        "F05_spatial_tissue_state_organization.png",
        "F08_ctc_emt_ngs_prior_panel.png",
        "no_yu_in_silico_visual_storyboard_cv2.png",
    ]
    try:
        import cv2
        import numpy as np

        tiles = []
        for name in imgs:
            p = HUB_ASSET / name
            if not p.exists():
                continue
            img = cv2.imread(str(p))
            if img is None:
                continue
            h, w = img.shape[:2]
            target_w = 720
            target_h = int(h * target_w / max(w, 1))
            img = cv2.resize(img, (target_w, target_h), interpolation=cv2.INTER_AREA)
            canvas = np.zeros((510, target_w, 3), dtype=np.uint8)
            canvas[:] = (18, 26, 44)
            y = min(465, target_h)
            canvas[:y, :, :] = img[:y, :, :]
            cv2.rectangle(canvas, (0, 466), (target_w - 1, 509), (38, 57, 84), -1)
            cv2.putText(canvas, name[:64], (18, 493), cv2.FONT_HERSHEY_SIMPLEX, 0.58, (210, 230, 245), 1, cv2.LINE_AA)
            tiles.append(canvas)
        if not tiles:
            return
        cols = 3
        rows = math.ceil(len(tiles) / cols)
        blank = tiles[0].copy()
        blank[:] = (8, 13, 24)
        while len(tiles) < rows * cols:
            tiles.append(blank.copy())
        sheet_rows = [np.concatenate(tiles[r * cols : (r + 1) * cols], axis=1) for r in range(rows)]
        sheet = np.concatenate(sheet_rows, axis=0)
        title_h = 120
        header = np.zeros((title_h, sheet.shape[1], 3), dtype=np.uint8)
        header[:] = (8, 13, 24)
        cv2.putText(header, "K-Thyro CTC-EMT Paper Flow Visual Atlas", (32, 52), cv2.FONT_HERSHEY_SIMPLEX, 1.35, (255, 210, 138), 2, cv2.LINE_AA)
        cv2.putText(header, "figure-first manuscript guide | data provenance | claim boundaries", (34, 92), cv2.FONT_HERSHEY_SIMPLEX, 0.78, (160, 178, 202), 1, cv2.LINE_AA)
        out = np.concatenate([header, sheet], axis=0)
        for dest in [FIG / "PF08_cv2_visual_atlas.png", HUB_ASSET / "PF08_cv2_visual_atlas.png", LIVE_ASSET / "PF08_cv2_visual_atlas.png"]:
            cv2.imwrite(str(dest), out)
    except Exception:
        return


def zip_packet(copied: list[Path]) -> None:
    seen: set[str] = set()
    with zipfile.ZipFile(ZIP_PATH, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in copied:
            if p.exists() and p.name not in seen:
                seen.add(p.name)
                zf.write(p, arcname=p.name)
        for p in [LOCAL_HTML, HTML]:
            if p.exists() and p.name not in seen:
                seen.add(p.name)
                zf.write(p, arcname=p.name)
    shutil.copy2(ZIP_PATH, HUB_ASSET / ZIP_NAME)
    shutil.copy2(ZIP_PATH, LIVE_ASSET / ZIP_NAME)


def html_page(
    data: pd.DataFrame,
    manifest: pd.DataFrame,
    flow: pd.DataFrame,
    claims: pd.DataFrame,
    explain: pd.DataFrame,
    captions: pd.DataFrame,
    scorecard: pd.DataFrame,
    blueprint: pd.DataFrame,
    strategy: pd.DataFrame,
    unlock: pd.DataFrame,
    main_figs: pd.DataFrame,
    methods: pd.DataFrame,
    crf: pd.DataFrame,
    risk: pd.DataFrame,
    sentence_bank: pd.DataFrame,
    action: pd.DataFrame,
    key_numbers: pd.DataFrame,
    sap: pd.DataFrame,
    gates: pd.DataFrame,
    reproducibility: pd.DataFrame,
    downgrade: pd.DataFrame,
    metrics: dict,
) -> str:
    lead_cards = [
        ("Final identity", "Science / mechanism", "one disease, one perturbation, one transition"),
        ("Paper possible now", "YES", "public-data prior/framework paper"),
        ("Discovery claim now", "NO", "needs Yu/hospital serial CTC raw data"),
        ("BRAF subset", str(metrics.get("braf_n")), f"largest state {metrics.get('braf_largest_pct'):.1f}%"),
        ("Spatial support", f"{metrics.get('spatial_z_gt2')}/{metrics.get('spatial_slides')}", "same-niche coherence z>2"),
        ("Output style", "Figure-first", "manifest + data dictionary + claim table"),
    ]
    metric_html = "<div class='metric-grid'>" + "".join(
        f"<div class='metric'><span>{esc(k)}</span><b>{esc(v)}</b><p>{esc(n)}</p></div>" for k, v, n in lead_cards
    ) + "</div>"

    visual_cards = [
        ("PF08_cv2_visual_atlas.png", "Master CV2 Visual Atlas", "한 장으로 보는 전체 논문 흐름과 핵심 그림 벽. 발표/카카오 전달용 첫 이미지.", "feature"),
        ("PF01_manuscript_flow.png", "Paper Logic Flow", "임상 gap에서 Yu anchor, 공개 prior, 병원 prospective test까지 이어지는 논문 골격.", ""),
        ("PF02_data_provenance_map.png", "Data Provenance Map", "각 데이터가 무엇을 지지하고 무엇을 지지하지 않는지 역할/한계를 함께 표시.", ""),
        ("PF03_marker_to_ngs_map.png", "Marker-to-NGS Map", "CTC phenotype만으로 끝내지 않고 tissue/cfDNA/CTC NGS clone anchor를 요구하는 구조.", ""),
        ("PF04_dataset_gap_bridge.png", "Missing Dataset Bridge", "public data로 가능한 부분과 Samsung/Yu cohort가 새로 만들어야 할 부분을 분리.", ""),
        ("PF05_claim_boundary_matrix.png", "Claim Boundary Matrix", "허용 claim과 금지 claim을 reviewer가 보기 전에 먼저 보여주는 방어 장치.", ""),
        ("PF06_figure_table_roadmap.png", "Figure/Table Roadmap", "각 그림과 표가 논문 어느 위치에 들어가는지 정리.", ""),
        ("PF07_reviewer_to_figure_map.png", "Reviewer Attack Map", "예상 공격, 답변, figure/table support를 한 화면에 연결.", ""),
        ("PF09_manuscript_claim_scorecard.png", "Manuscript Claim Scorecard", "인실리코 결과 중 논문 Results로 쓸 수 있는 claim과 강도, support figure를 정리.", ""),
        ("PF10_marker_module_support_heatmap.png", "Marker Module Planning Support Scores (Supplementary)", "marker module별 prior 구성 유용성을 손코딩한 manuscript-planning matrix (0/1/2). 분석 산출물이 아닌 기획 점수표이며, TCGA는 CTC 데이터 없음, Yu 2024 컬럼은 published summary만 사용.", ""),
        ("PF11_manuscript_blueprint_board.png", "Manuscript Blueprint Board", "지금 논문으로 쓸 구조와 쓰면 안 되는 구조를 첫 페이지에서 고정.", "feature"),
        ("PF12_future_data_unlock_map.png", "Future Data Unlock Map", "Yu/병원 데이터가 들어오면 어떤 future discovery figure와 claim이 열리는지 정의.", "feature"),
        ("PF13_execution_board.png", "Execution Board", "main figure layout, reviewer risk, 72시간 작업표를 한 장에 묶은 실행판.", "feature"),
        ("PF14_analysis_control_board.png", "Analysis Control Board", "본문 key numbers, 통계분석계획, stage-gate kill criteria를 한 화면에 고정.", "feature"),
        ("F03_tcga_braf_state_split.png", "TCGA BRAF State Split", "BRAF mutation 하나로 tissue-state를 설명할 수 없다는 공개 데이터 근거.", ""),
        ("F04_tcga_driver_variance_boundary.png", "TCGA Driver Variance Boundary", "driver-only model이 tissue-state axis를 전부 설명하지 못한다는 경계 근거.", ""),
        ("F05_spatial_tissue_state_organization.png", "Spatial Context", "공간적으로 조직화된 tissue-state 근거. 단, CTC origin proof는 아님.", ""),
        ("F06_scrna_marker_context.png", "scRNA Marker Context", "CTC marker prior의 cell-context rationale.", ""),
        ("F07_proteomics_dediff_direction.png", "Proteomics Direction", "protein-level dedifferentiation 방향성 보조 근거.", ""),
        ("F08_ctc_emt_ngs_prior_panel.png", "CTC-EMT-NGS Prior Panel", "E, E/M, M, thyroid-lineage, survival/stress, genetic anchor panel.", ""),
        ("F09_publishability_decision.png", "Publishability Decision", "지금 가능한 논문 유형과 불가능한 논문 유형을 구분.", ""),
    ]
    fig_html = "<div class='fig-grid'>" + "".join(image_card(*c) for c in visual_cards) + "</div>"

    links = [
        ("Packet ZIP", ZIP_NAME, "all figures + tables + reports"),
        ("Data dictionary TSV", "data_dictionary_master.tsv", "dataset definitions"),
        ("Figure manifest TSV", "figure_table_manifest.tsv", "where each figure/table is used"),
        ("Figure explanation TSV", "figure_explanation_master.tsv", "role/data/readout/boundary"),
        ("Caption scaffold TSV", "figure_caption_and_results_scaffold.tsv", "captions + results text"),
        ("Claim scorecard TSV", "manuscript_claim_evidence_scorecard.tsv", "writeable in silico claims"),
        ("Marker support TSV (supplementary)", "marker_module_cross_layer_support_matrix.tsv", "hand-coded planning judgement (0/1/2); not data-derived"),
        ("Manuscript blueprint TSV", "manuscript_blueprint.tsv", "safe paper structure"),
        ("Submission strategy TSV", "submission_strategy_matrix.tsv", "which paper route is feasible"),
        ("Future data unlock TSV", "future_data_unlock_map.tsv", "what hospital/Yu data unlocks"),
        ("Main/Supp layout TSV", "main_supp_figure_layout.tsv", "paper figure assembly"),
        ("Methods checklist TSV", "methods_reproducibility_checklist.tsv", "input/method/output QC"),
        ("Hospital CRF TSV", "hospital_data_request_crf.tsv", "fields to request from hospital"),
        ("Reviewer risk TSV", "reviewer_risk_heatmap.tsv", "risk severity and mitigation"),
        ("Safe sentence TSV", "safe_sentence_bank.tsv", "allowed vs forbidden phrasing"),
        ("72h execution TSV", "next_72h_manuscript_execution.tsv", "rapid execution plan"),
        ("Key number ledger TSV", "key_number_ledger.tsv", "numbers used in manuscript"),
        ("Statistical plan TSV", "statistical_analysis_plan.tsv", "current vs future analysis plan"),
        ("Stage gates TSV", "stage_gate_kill_criteria.tsv", "kill/downgrade criteria"),
        ("Repro manifest TSV", "reproducibility_manifest.tsv", "scripts and artifacts"),
        ("Scenario tree TSV", "scenario_decision_tree.tsv", "what to do under data scenarios"),
        ("Manuscript flow TSV", "manuscript_flow_outline.tsv", "paper section flow"),
        ("Claim boundary TSV", "claim_boundary_matrix.tsv", "allowed vs forbidden claims"),
        ("Manuscript guide MD", "MANUSCRIPT_FLOW_MASTER_KR.md", "Korean paper logic"),
        ("Data guide MD", "DATA_DEFINITION_MASTER_KR.md", "data source definitions"),
        ("Figure guide MD", "FIGURE_TABLE_EXPLANATION_GUIDE_KR.md", "figure/table explanation"),
        ("Writing scaffold MD", "MANUSCRIPT_READY_WRITING_SCAFFOLD_KR.md", "title/abstract/results direction"),
        ("Marker support MD", "MARKER_MODULE_CROSS_LAYER_SUPPORT_KR.md", "in silico marker support matrix"),
        ("Manuscript blueprint MD", "MANUSCRIPT_BLUEPRINT_KR.md", "paper route and section plan"),
        ("Future data unlock MD", "FUTURE_DATA_UNLOCK_MAP_KR.md", "future figures and claims"),
        ("Execution package MD", "MANUSCRIPT_EXECUTION_PACKAGE_KR.md", "figure layout, CRF, risk, 72h plan"),
        ("Analysis control MD", "ANALYSIS_CONTROL_PACKAGE_KR.md", "numbers, stats, gates, reproducibility"),
        ("Caption scaffold MD", "FIGURE_CAPTION_AND_RESULTS_SCAFFOLD_KR.md", "caption bank"),
        ("Samsung PPTX", "kthyro_ctc_emt_samsung_science_pitch_v2_reference_style.pptx", "deck download"),
    ]
    link_html = "<div class='links'>" + "".join(link_card(*x) for x in links) + "</div>"

    allowed = claims[claims["status"].eq("Allowed")]
    forbidden = claims[claims["status"].eq("Forbidden")]

    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>K-Thyro CTC-EMT Paper Flow Master Dossier</title>
<style>
:root{{--bg:{BG};--panel:{PANEL};--panel2:{PANEL2};--line:{LINE};--text:{TEXT};--muted:{MUTED};--gold:{GOLD};--cyan:{CYAN};--teal:{TEAL};--purple:{PURPLE};--orange:{ORANGE};--warn:{WARN};}}
*{{box-sizing:border-box}} html{{scroll-behavior:smooth}}
body{{margin:0;background:var(--bg);color:var(--text);font-family:Inter,"Noto Sans KR",-apple-system,sans-serif;font-size:14px;line-height:1.58}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
h1,h2,h3{{font-family:"Cormorant Garamond","Noto Serif KR",serif;color:#fff8e7;letter-spacing:-.35px}}
code,pre,.mono{{font-family:"JetBrains Mono","SF Mono",Menlo,monospace}}
.hero{{min-height:78vh;padding:70px 34px 42px;border-bottom:1px solid var(--line);background:linear-gradient(135deg,#090e1a,#101827 48%,#061321),radial-gradient(ellipse at 18% 10%,rgba(91,220,255,.18),transparent 36%),radial-gradient(ellipse at 82% 12%,rgba(255,210,138,.16),transparent 32%)}}
.hero-inner{{max-width:1440px;margin:0 auto}}
.kicker{{font:800 10px/1.4 "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.22em;text-transform:uppercase;margin-bottom:16px}}
h1{{font-size:64px;line-height:.98;margin:0 0 18px;max-width:1120px}} h1 em{{color:var(--gold);font-style:italic}}
.lead{{font-size:18px;max-width:1120px;color:#d8e2ef;font-family:"Noto Serif KR",serif;line-height:1.72}}
.hero-grid{{display:grid;grid-template-columns:1.05fr .95fr;gap:22px;margin-top:26px;align-items:start}}
.hero-grid img{{width:100%;border:1px solid rgba(91,220,255,.34);border-radius:14px;box-shadow:0 24px 80px rgba(0,0,0,.35);background:#050811}}
.actions{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin:18px 0}}
.actions a{{display:block;border:1px solid rgba(255,210,138,.38);border-radius:12px;padding:16px;background:linear-gradient(135deg,rgba(91,220,255,.14),rgba(255,210,138,.08));text-decoration:none}}
.actions b{{display:block;color:#fff8e7;font-size:20px;font-family:"Cormorant Garamond","Noto Serif KR",serif}} .actions span{{font:10px/1.35 "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.05em;text-transform:uppercase}}
.urlbox{{border:1px solid rgba(91,220,255,.3);border-left:4px solid var(--cyan);background:#091321;border-radius:0 10px 10px 0;padding:10px 13px;margin:10px 0;font:12px/1.45 "JetBrains Mono",monospace;color:#cfe9ff;overflow-wrap:anywhere}}
.metric-grid{{display:grid;grid-template-columns:repeat(6,1fr);gap:10px;margin:22px 0 0}}
.metric{{background:rgba(255,255,255,.045);border:1px solid rgba(255,210,138,.2);border-radius:12px;padding:13px}}
.metric span{{display:block;color:var(--muted);font:800 9px/1.35 "JetBrains Mono",monospace;letter-spacing:.08em;text-transform:uppercase}} .metric b{{display:block;color:var(--gold);font:800 27px/1 "Cormorant Garamond",serif;margin:6px 0}} .metric p{{margin:0;color:#d6e1ee;font-size:12px}}
.wrap{{max-width:1460px;margin:0 auto;display:grid;grid-template-columns:260px 1fr;gap:36px;padding:0 32px 76px}}
.toc{{position:sticky;top:0;align-self:start;max-height:100vh;overflow:auto;padding:28px 0;border-right:1px solid var(--line)}}
.toc h4{{font:800 10px/1.4 "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;text-transform:uppercase;margin:0 0 12px}} .toc ol{{list-style:none;counter-reset:t;margin:0;padding:0}} .toc li{{counter-increment:t;margin:5px 0;font-size:12px}} .toc li:before{{content:counter(t,decimal-leading-zero) "  ";font-family:"JetBrains Mono",monospace;color:var(--gold);font-size:10px}}
main{{min-width:0;padding-top:28px}} section{{padding:30px 0;border-bottom:1px solid var(--line)}} h2{{font-size:36px;margin:0 0 6px}} h2 .num{{font:400 14px/1 "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:12px}} .sub{{font:11px/1.45 "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.09em;text-transform:uppercase;margin-bottom:18px}}
.box{{background:#0e1728;border:1px solid var(--line);border-left:3px solid var(--gold);border-radius:0 9px 9px 0;padding:15px 17px;margin:14px 0}} .box.warn{{border-left-color:var(--warn);background:#1d1410}} .box.good{{border-left-color:var(--teal);background:#0e1d18}} .box h3{{font:800 12px/1.4 "JetBrains Mono",monospace;letter-spacing:.1em;text-transform:uppercase;color:var(--gold);margin:0 0 8px}}
.fig-grid{{display:grid;grid-template-columns:repeat(2,1fr);gap:18px;margin-top:16px}} .fig-card{{background:#0e1728;border:1px solid var(--line);border-radius:13px;overflow:hidden;box-shadow:0 14px 38px rgba(0,0,0,.24)}} .fig-card.feature{{grid-column:1/-1;border-color:rgba(255,210,138,.5)}} .fig-head{{padding:12px 14px;background:#121f35;border-bottom:1px solid var(--line);display:flex;justify-content:space-between;gap:12px;align-items:center}} .fig-head b{{color:#fff8e7}} .fig-head span{{font:9px/1.2 "JetBrains Mono",monospace;color:var(--muted);overflow-wrap:anywhere;text-align:right}} .fig-card img{{display:block;width:100%;height:auto;background:#050811;border-bottom:1px solid var(--line)}} .fig-card p{{margin:0;padding:12px 14px 14px;color:#d5dfec;font-size:13px}}
.explain-list{{display:grid;gap:18px;margin-top:24px}} .explain-card{{display:grid;grid-template-columns:360px 1fr;gap:0;background:#0e1728;border:1px solid rgba(255,210,138,.22);border-radius:14px;overflow:hidden;box-shadow:0 16px 42px rgba(0,0,0,.25)}} .explain-img{{background:#050811;border-right:1px solid var(--line);display:flex;align-items:center}} .explain-img img{{display:block;width:100%;height:auto}} .explain-body{{padding:18px 20px}} .tagline{{display:flex;align-items:center;gap:12px;margin-bottom:14px}} .tagline span{{font:900 11px/1 "JetBrains Mono",monospace;color:#07101a;background:var(--gold);border-radius:999px;padding:6px 9px}} .tagline b{{font-size:24px;color:#fff8e7;font-family:"Cormorant Garamond","Noto Serif KR",serif}} .mini-grid{{display:grid;grid-template-columns:1fr 1fr;gap:10px;margin-bottom:12px}} .mini-grid div{{border:1px solid var(--line);background:#111d31;border-radius:8px;padding:10px}} .mini-grid h4{{margin:0 0 5px;color:var(--gold);font:900 9px/1.2 "JetBrains Mono",monospace;letter-spacing:.09em;text-transform:uppercase}} .mini-grid p,.script p,.boundary p{{margin:0;color:#d8e2ef}} .script{{border-left:3px solid var(--teal);background:#0d1d19;border-radius:0 8px 8px 0;padding:11px 13px;margin:10px 0}} .script b{{color:var(--teal)}} .boundary{{border-left:3px solid var(--warn);background:#21130f;border-radius:0 8px 8px 0;padding:11px 13px}} .boundary b{{color:var(--warn)}}
.caption-list{{display:grid;gap:14px;margin-top:18px}} .caption-card{{background:#0e1728;border:1px solid rgba(91,220,255,.25);border-radius:12px;padding:16px}} .caption-title{{display:flex;align-items:flex-start;gap:12px;margin-bottom:12px}} .caption-title span{{flex:0 0 auto;font:900 10px/1 "JetBrains Mono",monospace;color:#07101a;background:var(--cyan);border-radius:999px;padding:7px 9px}} .caption-title b{{font-family:"Cormorant Garamond","Noto Serif KR",serif;font-size:21px;color:#fff8e7;line-height:1.15}} .caption-grid{{display:grid;grid-template-columns:1.2fr 1.2fr 1fr 1fr;gap:10px}} .caption-grid.embedded{{grid-template-columns:1fr 1fr;margin-top:8px}} .caption-grid div{{background:#111d31;border:1px solid var(--line);border-radius:8px;padding:10px}} .caption-grid h4{{margin:0 0 6px;color:var(--gold);font:900 9px/1.2 "JetBrains Mono",monospace;letter-spacing:.08em;text-transform:uppercase}} .caption-grid p{{margin:0;color:#d8e2ef;font-size:12.5px}} .embedded-caption{{margin-top:13px;padding-top:12px;border-top:1px solid rgba(255,210,138,.22)}} .embedded-caption>h4{{margin:0 0 8px;color:var(--cyan);font:900 10px/1.2 "JetBrains Mono",monospace;letter-spacing:.1em;text-transform:uppercase}} .muted-cap{{background:#0a111d;border:1px dashed var(--line);border-radius:8px;padding:10px}} .muted-cap p{{margin:0;color:var(--muted)}}
table.t{{border-collapse:collapse;width:100%;font-size:12.2px;margin:10px 0 16px}} .t th,.t td{{border:1px solid var(--line);padding:8px 9px;text-align:left;vertical-align:top}} .t th{{background:#16213a;color:var(--gold);font:800 10px/1.3 "JetBrains Mono",monospace;letter-spacing:.04em;text-transform:uppercase}} .t tr:nth-child(even) td{{background:rgba(255,255,255,.025)}} .t td.num{{font-family:"JetBrains Mono",monospace;text-align:right;color:var(--gold)}}
.links{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}} .link-card{{display:block;background:#101a2c;border:1px solid var(--line);border-radius:8px;padding:12px;min-height:76px;text-decoration:none}} .link-card b{{display:block;color:var(--gold);font:800 10px/1.4 "JetBrains Mono",monospace;letter-spacing:.08em;text-transform:uppercase;margin-bottom:4px}} .link-card span{{color:#d5dfec;font-size:12px;overflow-wrap:anywhere}}
.grid-2{{display:grid;grid-template-columns:1fr 1fr;gap:16px}} .muted{{color:var(--muted)}} footer{{padding:26px 32px;border-top:1px solid var(--line);text-align:center;color:var(--muted);font:11px/1.4 "JetBrains Mono",monospace}}
@media(max-width:1100px){{.wrap,.hero-grid,.grid-2,.explain-card,.mini-grid,.caption-grid{{grid-template-columns:1fr}}.toc{{display:none}}.fig-grid,.actions,.metric-grid,.links{{grid-template-columns:1fr}} h1{{font-size:42px}} .hero{{min-height:0}} .explain-img{{border-right:0;border-bottom:1px solid var(--line)}}}}
</style>
</head>
<body>
<header class="hero"><div class="hero-inner">
  <div class="kicker">K-Thyro · CTC-EMT-NGS · paper-flow master dossier · 2026-05-09</div>
  <h1>논문 흐름, 그림, 표, 데이터 정의서를 <em>한 페이지로 고정</em></h1>
  <p class="lead">이 페이지는 지금까지 만든 Samsung Science proposal, no-Yu-data in silico 분석, 공개 pilot 결과를 논문 작성 순서로 재배열한 master dossier다. 각 그림과 표마다 <b>어떤 데이터를 썼는지, 본문에서 무슨 말을 해야 하는지, 어디까지 말하면 안 되는지</b>를 함께 적었다.</p>
  <div class="hero-grid">
    <div>
      <div class="actions">
        <a href="assets/{ASSET}/{ZIP_NAME}"><b>Packet ZIP</b><span>figures + reports + tables</span></a>
        <a href="#figures"><b>Figure Wall</b><span>visual-first manuscript map</span></a>
        <a href="#manifest"><b>Manifest</b><span>figure/table/data definitions</span></a>
      </div>
      <div class="urlbox">전달용 웹: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_paper_flow_master_dossier.html</div>
      <div class="urlbox">패키지 ZIP: http://40.82.129.113:8012/papers_hub_2026_05_04/assets/{ASSET}/{ZIP_NAME}</div>
      {metric_html}
    </div>
    <a href="assets/{ASSET}/PF08_cv2_visual_atlas.png"><img src="assets/{ASSET}/PF08_cv2_visual_atlas.png" alt="CV2 visual atlas" /></a>
  </div>
</div></header>
<div class="wrap">
<aside class="toc"><h4>Contents</h4><ol>
  <li><a href="#verdict">Verdict</a></li>
  <li><a href="#story">Paper Story</a></li>
  <li><a href="#figures">Figure Wall</a></li>
  <li><a href="#writing">Writing Pack</a></li>
  <li><a href="#data">Data Dictionary</a></li>
  <li><a href="#manifest">Figure/Table Manifest</a></li>
  <li><a href="#flow">Manuscript Flow</a></li>
  <li><a href="#claims">Claim Boundaries</a></li>
  <li><a href="#downloads">Downloads</a></li>
</ol></aside>
<main>
<section id="verdict">
  <h2><span class="num">01</span>Final Verdict</h2>
  <div class="sub">Science-track proposal remains Science; paper is possible now only as a public-data prior/framework paper.</div>
  <div class="grid-2">
    <div class="box good"><h3>지금 쓸 수 있는 논문</h3><p><b>PTC CTC-EMT-NGS prior map / framework paper.</b> 공개 tissue-state, spatial, scRNA, proteomics 근거를 이용해 post-thyroidectomy CTC-EMT monitoring을 위한 marker/clone-anchor 논리를 정리한다.</p></div>
    <div class="box warn"><h3>지금 쓰면 안 되는 논문</h3><p><b>새로운 postoperative CTC-EMT transition discovery paper.</b> Yu/병원 raw serial CTC table과 matched NGS가 없으면 transition, residual-risk prediction, recurrence association은 주장 금지다.</p></div>
  </div>
</section>
<section id="story">
  <h2><span class="num">02</span>Paper Story</h2>
  <div class="sub">The exact order in which the reader should encounter the logic.</div>
  {df_table(flow, max_rows=20)}
</section>
<section id="figures">
  <h2><span class="num">03</span>Figure Wall</h2>
  <div class="sub">Visual-first, but not decoration: every figure below now carries role, source data, reading guide, manuscript use, caption, results paragraph, provenance note, reviewer-safe sentence, and claim boundary.</div>
  <div class="box good"><h3>읽는 법</h3><p>각 그림 카드는 그 자리에서 끝나게 만들었습니다. <b>Role</b>은 왜 필요한지, <b>Data Used</b>는 무슨 데이터를 썼는지, <b>How To Read</b>는 그림에서 봐야 할 신호, <b>Paper Use</b>는 본문 위치, <b>설명 문장</b>은 발표용 해석, <b>Claim Boundary</b>는 넘으면 안 되는 선입니다. 바로 아래 <b>Caption And Results Scaffold</b>에는 그 그림의 caption, Results paragraph, methods/provenance note, reviewer-safe sentence가 붙습니다.</p></div>
  {fig_html}
  <h3>Figure-by-Figure Explanation Master</h3>
  {explanation_blocks(explain, captions)}
</section>
<section id="writing">
  <h2><span class="num">04</span>Manuscript Writing Pack</h2>
  <div class="sub">In silico outputs converted into writeable Results claims, captions, and reviewer-safe sentences.</div>
  <div class="box good"><h3>논문으로 쓰는 법</h3><p>현재 논문은 <b>public-data prior/framework paper</b>로 쓴다. Results는 mutation/state heterogeneity, spatial organization, marker context, marker-to-NGS anchoring, missing dataset definition, publishability boundary 순서로 간다. CTC transition과 recurrence prediction은 Discussion의 future validation으로만 둔다.</p></div>
  {image_card("PF11_manuscript_blueprint_board.png", "Manuscript Blueprint Board", "논문 identity, title, thesis, Results block, boundary를 첫 페이지에서 고정한다.", "feature")}
  <h3>Manuscript Blueprint</h3>
  {df_table(blueprint, max_rows=20)}
  <h3>Submission Strategy</h3>
  {df_table(strategy, max_rows=10)}
  {image_card("PF09_manuscript_claim_scorecard.png", "In Silico Manuscript Claim Scorecard", "각 result block이 어떤 in silico evidence와 figure support를 가지는지, 그리고 왜 아직 discovery claim이 아닌지를 정리한다.", "feature")}
  <h3>Claim-Evidence Scorecard</h3>
  {df_table(scorecard, max_rows=12)}
  {image_card("PF12_future_data_unlock_map.png", "Future Data Unlock Map", "Yu/병원 데이터가 들어오면 어떤 future discovery figure와 claim이 가능해지는지 정의한다.", "feature")}
  <h3>Future Data Unlock Map</h3>
  {df_table(unlock, max_rows=12)}
  {image_card("PF13_execution_board.png", "Execution Board", "main/supp figure layout, reviewer risk severity, next 72-hour manuscript actions를 한 장으로 묶었다.", "feature")}
  <h3>Main/Supp Figure Layout</h3>
  {df_table(main_figs, max_rows=14)}
  <h3>Methods Reproducibility Checklist</h3>
  {df_table(methods, max_rows=14)}
  <h3>Hospital Data Request CRF</h3>
  {df_table(crf, max_rows=14)}
  <h3>Reviewer Risk Heatmap</h3>
  {df_table(risk, max_rows=14)}
  <h3>Safe Sentence Bank</h3>
  {df_table(sentence_bank, max_rows=14)}
  <h3>Next 72h Manuscript Execution</h3>
  {df_table(action, max_rows=10)}
  {image_card("PF14_analysis_control_board.png", "Analysis Control Board", "본문에 들어갈 key numbers, 현재/미래 statistical analysis plan, stage-gate kill criteria를 한 장으로 고정한다.", "feature")}
  <h3>Key Number Ledger</h3>
  {df_table(key_numbers, max_rows=14)}
  <h3>Statistical Analysis Plan</h3>
  {df_table(sap, max_rows=14)}
  <h3>Stage-Gate Kill Criteria</h3>
  {df_table(gates, max_rows=12)}
  <h3>Reproducibility Manifest</h3>
  {df_table(reproducibility, max_rows=12)}
  <h3>Scenario Decision Tree</h3>
  {df_table(downgrade, max_rows=10)}
  <h3>Caption And Results Scaffold</h3>
  {caption_blocks(captions)}
</section>
<section id="data">
  <h2><span class="num">05</span>Data Dictionary</h2>
  <div class="sub">What data was used, what each layer supports, and what it cannot support.</div>
  {df_table(data, max_rows=20)}
</section>
<section id="manifest">
  <h2><span class="num">06</span>Figure/Table Manifest</h2>
  <div class="sub">This is the writing-control table: no figure floats without a data source and claim boundary.</div>
  {df_table(manifest, max_rows=30)}
</section>
<section id="flow">
  <h2><span class="num">07</span>Explanation Flow</h2>
  <div class="sub">How to explain the figures in a manuscript or Samsung reviewer meeting.</div>
  <div class="box"><h3>추천 설명 순서</h3><p>1. 먼저 Yu 2024로 CTC-EMT 측정 가능성을 인정한다. 2. 바로 “하지만 public integrated dataset은 없다”로 gap을 만든다. 3. TCGA/spatial/scRNA/proteomics는 CTC proof가 아니라 marker-prior라고 제한한다. 4. 그래서 matched tissue-normal NGS와 serial blood가 proposal core라고 닫는다.</p></div>
  {image_card("PF01_manuscript_flow.png", "Narrative Anchor", "이 그림을 첫 설명 그림으로 쓰면 흐름이 가장 안전하다.", "feature")}
</section>
<section id="claims">
  <h2><span class="num">08</span>Claim Boundaries</h2>
  <div class="sub">Allowed and forbidden statements, explicitly separated.</div>
  <div class="grid-2">
    <div><h3>Allowed</h3>{df_table(allowed, max_rows=10)}</div>
    <div><h3>Forbidden</h3>{df_table(forbidden, max_rows=10)}</div>
  </div>
  {image_card("PF05_claim_boundary_matrix.png", "Claim Boundary Figure", "리뷰어 방어용 핵심 그림. 과장하지 않는다는 신뢰를 만든다.", "")}
</section>
<section id="downloads">
  <h2><span class="num">09</span>Downloads</h2>
  <div class="sub">Direct assets for paper writing, professor update, Samsung pitch, and vendor/hospital follow-up.</div>
  {link_html}
</section>
</main>
</div>
<footer>Generated from local outputs under kthyro_ctc_emt_samsung_proposal and kthyro_public_pilot/results. Public pilot supports prior-map logic only, not CTC transition discovery.</footer>
</body>
</html>
"""


def update_index() -> None:
    idx = HUB / "index.html"
    if not idx.exists():
        return
    text = idx.read_text(encoding="utf-8")
    needle = '    <a href="kthyro_no_yu_in_silico_ctc_emt_dossier.html"'
    if "kthyro_ctc_emt_paper_flow_master_dossier.html" in text:
        return
    insert = '    <a href="kthyro_ctc_emt_paper_flow_master_dossier.html" style="background:linear-gradient(135deg,#07101a 0%,#1a2036 52%,#241915 100%);border:2px solid #ffd28a;color:#fff2d5;font-weight:900;font-size:16px;padding:14px 18px">★ CTC-EMT PAPER FLOW · figures/tables/data dictionary</a>\n'
    pos = text.find(needle)
    if pos >= 0:
        text = text[:pos] + insert + text[pos:]
    else:
        text += "\n" + insert
    idx.write_text(text, encoding="utf-8")
    shutil.copy2(idx, LIVE / "index.html")


def main() -> None:
    mkdirs()
    metrics = collect_metrics()
    data, manifest, flow, claims = write_master_tables(metrics)
    explain = write_figure_explanation_table()
    captions = write_caption_bank()
    scorecard = write_manuscript_claim_package(metrics)
    marker_support = write_marker_cross_layer_support()
    blueprint, strategy = write_manuscript_blueprint()
    unlock = write_future_data_unlock_map()
    main_figs, methods, crf, risk, sentence_bank, action = write_execution_package()
    key_numbers, sap, gates, reproducibility, downgrade = write_analysis_control_package(metrics)
    plot_master_figures(metrics, data, manifest, flow, claims)
    plot_manuscript_claim_figure(scorecard)
    plot_marker_support_heatmap(marker_support)
    plot_manuscript_blueprint_board(blueprint, strategy)
    plot_future_data_unlock_map(unlock)
    plot_execution_board(main_figs, risk, action)
    plot_analysis_control_board(key_numbers, sap, gates)
    write_reports(data, manifest, flow, claims, explain, captions, scorecard, metrics)
    copied = copy_source_assets()
    build_cv2_atlas()
    # Copy atlas after creation if it was generated after copy_source_assets.
    for p in [FIG / "PF08_cv2_visual_atlas.png"]:
        if p.exists():
            shutil.copy2(p, SRC / p.name)
            shutil.copy2(p, HUB_ASSET / p.name)
            shutil.copy2(p, LIVE_ASSET / p.name)
            copied.append(SRC / p.name)
    html_text = html_page(
        data,
        manifest,
        flow,
        claims,
        explain,
        captions,
        scorecard,
        blueprint,
        strategy,
        unlock,
        main_figs,
        methods,
        crf,
        risk,
        sentence_bank,
        action,
        key_numbers,
        sap,
        gates,
        reproducibility,
        downgrade,
        metrics,
    )
    HTML.write_text(html_text, encoding="utf-8")
    LOCAL_HTML.write_text(html_text, encoding="utf-8")
    shutil.copy2(HTML, LIVE_HTML)
    copied.append(LOCAL_HTML)
    zip_packet(copied)
    update_index()
    print(f"HTML={HTML}")
    print(f"LIVE={LIVE_HTML}")
    print(f"LOCAL={LOCAL_HTML}")
    print(f"ZIP={ZIP_PATH}")
    print(f"ASSET={LIVE_ASSET}")


if __name__ == "__main__":
    main()
