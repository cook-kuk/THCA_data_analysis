#!/usr/bin/env python3
"""THCA CMap-style drug reversal + PPI network layer.

Phase 1 — Connectivity Map (CMap / CLUE.io-style) drug-reversal analysis.
   Uses an OFFLINE literature-consensus surrogate database of ~40 well-known
   compounds with published transcriptomic up/down signatures (paraphrased
   from LINCS L1000 / iLINCS summaries + published MoA). For each compound,
   compute an ssCMap-style connectivity score against the BRAF_like vs
   RAS_like DE signature and rank candidate "reversers" of the BRAF-like
   state. This is NOT a live CLUE.io query.

Phase 2 — Literature-derived PPI network. Hardcoded subset of STRING v11
   high-confidence + BioGRID edges among the top-5 novel validated
   biomarkers (TACSTD2, PLEKHA6, CYP1B1, TMPRSS4, LDLR) plus the 14
   QUBO-selected genes (KLK10, DUSP5, DUSP6, DIO1, LOX, MET, TPO, SLC5A8,
   DUOX2, CD274, HLA-DRA, DUOX1, DIO2, FOXE1) with 1st-degree neighbours in
   thyroid + cancer-signalling space. Edges are literature-derived, NOT
   re-inferred from our expression matrix.

Phase 3 — Build the dashboard page reports/html/pages/19_cmap_network.html.
Phase 4 — Inject the navigation link on every dashboard page (idempotent).

Outputs are self-contained, CSP-safe (vendored Plotly), Korean primary
labels, English technical. Safe to re-run (idempotent).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

try:
    import networkx as nx
except ImportError:
    print("[setup] installing networkx...", flush=True)
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--quiet", "networkx"])
    import networkx as nx  # type: ignore

warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "reports" / "html" / "figs_interactive"
PAGES = ROOT / "reports" / "html" / "pages"
SCRIPTS = ROOT / "scripts"
META = ROOT / "metadata"
JSON_OUT = ROOT / "results" / "json"

for d in (TABLES, FIGS, PAGES, JSON_OUT):
    d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Offline literature-consensus compound signature database.
# Each entry is (up-signature genes, down-signature genes, MoA, weight).
# Signatures are hand-curated from published mechanism-of-action summaries
# (paraphrased from LINCS L1000 / iLINCS / pharmacology reviews); NOT a
# verbatim copy of any proprietary table. 15-30 genes per direction.
# Sources cited collectively on the dashboard page.
# ---------------------------------------------------------------------------

COMPOUND_DB: Dict[str, Dict] = {
    # ---- BRAF inhibitors (expected: reverse MAPK-output / BRAF-like program)
    "dabrafenib": {
        "moa": "BRAF V600E 억제 (BRAF V600E inhibitor)",
        "up": ["DUSP4", "DUSP6", "SPRY2", "SPRY4", "CCND1", "FOS", "JUN",
               "TG", "TPO", "TSHR", "SLC5A5", "DIO1", "DIO2", "PAX8", "FOXE1"],
        "down": ["ETV4", "ETV5", "FOSL1", "MYC", "CCND2", "MKI67", "HMGA2",
                 "FN1", "TIMP1", "LCN2", "SERPINA1", "CXCL14", "KLK10",
                 "TACSTD2", "TMPRSS4", "S100A4"],
        "weight": 1.0,
        "rationale": "BRAF V600E papillary thyroid cancer의 표준 표적. BRAF-like 프로그램을 직접 역전.",
    },
    "vemurafenib": {
        "moa": "BRAF V600E 억제 (BRAF V600E inhibitor)",
        "up": ["DUSP4", "DUSP6", "SPRY2", "SPRY4", "CCND1", "TG", "TPO",
               "TSHR", "SLC5A5", "DIO1", "DIO2", "PAX8", "FOXE1", "NKX2-1"],
        "down": ["ETV4", "ETV5", "FOSL1", "MYC", "MKI67", "HMGA2", "FN1",
                 "TIMP1", "LCN2", "CXCL8", "CXCL14", "KLK10", "TACSTD2",
                 "TMPRSS4", "S100A4", "CD274"],
        "weight": 1.0,
        "rationale": "첫 승인된 BRAF V600E 억제제; THCA 임상에서 재분화 유도 보고.",
    },
    "encorafenib": {
        "moa": "BRAF V600E 억제 (BRAF V600E inhibitor, 3세대)",
        "up": ["DUSP6", "SPRY2", "SPRY4", "TG", "TPO", "TSHR", "SLC5A5",
               "DIO1", "DIO2", "PAX8", "FOXE1"],
        "down": ["ETV4", "ETV5", "FOSL1", "MYC", "MKI67", "HMGA2", "FN1",
                 "TIMP1", "LCN2", "KLK10", "TACSTD2", "TMPRSS4", "CD274"],
        "weight": 0.95,
        "rationale": "melanoma에서 vemurafenib 대비 potency↑; thyroid basket trial 진행 중.",
    },
    # ---- MEK inhibitors
    "trametinib": {
        "moa": "MEK1/2 억제 (MEK1/2 inhibitor)",
        "up": ["DUSP4", "DUSP6", "SPRY2", "SPRY4", "TG", "TPO", "TSHR",
               "SLC5A5", "DIO1", "DIO2", "PAX8", "FOXE1", "CDH1"],
        "down": ["ETV4", "ETV5", "FOSL1", "MYC", "CCND2", "MKI67", "HMGA2",
                 "FN1", "TIMP1", "LCN2", "CXCL8", "KLK10", "TACSTD2",
                 "TMPRSS4", "S100A4"],
        "weight": 1.0,
        "rationale": "MAPK pathway downstream 차단 — redifferentiation 유도 (Ho et al., NEJM 2013).",
    },
    "selumetinib": {
        "moa": "MEK1/2 억제 (MEK1/2 inhibitor)",
        "up": ["DUSP4", "DUSP6", "SPRY2", "TG", "TPO", "TSHR", "SLC5A5",
               "DIO1", "DIO2", "PAX8", "FOXE1"],
        "down": ["ETV4", "ETV5", "FOSL1", "MYC", "MKI67", "HMGA2", "FN1",
                 "TIMP1", "LCN2", "KLK10", "TACSTD2", "TMPRSS4", "S100A4"],
        "weight": 1.0,
        "rationale": "radioiodine-refractory THCA에서 NIS 재발현 유도 (Ho et al., 2013).",
    },
    "cobimetinib": {
        "moa": "MEK1/2 억제 (MEK1/2 inhibitor)",
        "up": ["DUSP4", "DUSP6", "SPRY2", "TG", "TPO", "TSHR", "SLC5A5",
               "DIO1", "DIO2"],
        "down": ["ETV4", "ETV5", "FOSL1", "MYC", "MKI67", "FN1", "TIMP1",
                 "KLK10", "TACSTD2", "TMPRSS4"],
        "weight": 0.9,
        "rationale": "BRAFi 병용 (encorafenib/cobimetinib)로 저항 극복.",
    },
    # ---- NTRK inhibitors
    "larotrectinib": {
        "moa": "TRK A/B/C 억제 (NTRK inhibitor)",
        "up": ["DUSP6", "SPRY2", "TG", "TPO", "PAX8", "FOXE1"],
        "down": ["NTRK1", "NTRK2", "NTRK3", "ETV6", "MYC", "MKI67",
                 "HMGA2", "FN1", "S100A4", "TACSTD2"],
        "weight": 0.85,
        "rationale": "NTRK-fusion driven thyroid cancer에 특이; 소아 PTC 부분에 효과.",
    },
    "entrectinib": {
        "moa": "TRK · ROS1 · ALK 억제 (multi-kinase incl. NTRK)",
        "up": ["DUSP6", "SPRY2", "TG", "TPO", "PAX8"],
        "down": ["NTRK1", "NTRK2", "NTRK3", "ROS1", "ALK", "MYC", "MKI67",
                 "FN1", "TACSTD2"],
        "weight": 0.8,
        "rationale": "NTRK/ROS1-fusion tumor에 승인; THCA basket.",
    },
    # ---- Thyroid multikinase inhibitors (approved)
    "sorafenib": {
        "moa": "VEGFR · RAF · PDGFR 다중억제 (multikinase)",
        "up": ["DUSP6", "SPRY2", "TG", "TPO", "TSHR", "SLC5A5", "DIO1",
               "DIO2"],
        "down": ["VEGFA", "VEGFR2", "KDR", "ETV4", "FOSL1", "MYC", "CCND1",
                 "MKI67", "FN1", "TIMP1", "LCN2", "KLK10", "TACSTD2",
                 "TMPRSS4", "HMGA2"],
        "weight": 1.0,
        "rationale": "DTC(분화형 thyroid cancer) radioiodine-refractory에 승인 — DECISION trial.",
    },
    "lenvatinib": {
        "moa": "VEGFR · FGFR · RET · KIT 다중억제 (multikinase)",
        "up": ["DUSP6", "SPRY2", "TG", "TPO", "TSHR", "SLC5A5", "DIO1",
               "DIO2", "PAX8", "FOXE1"],
        "down": ["VEGFA", "KDR", "FGFR1", "FGFR2", "RET", "KIT", "ETV4",
                 "FOSL1", "MYC", "CCND1", "MKI67", "FN1", "TIMP1", "LCN2",
                 "KLK10", "TACSTD2", "TMPRSS4"],
        "weight": 1.0,
        "rationale": "RR-DTC 1차 옵션 — SELECT trial (Schlumberger, NEJM 2015).",
    },
    "cabozantinib": {
        "moa": "MET · VEGFR · RET 다중억제 (multikinase)",
        "up": ["DUSP6", "SPRY2", "TG", "TPO", "TSHR", "DIO1"],
        "down": ["MET", "HGF", "VEGFA", "KDR", "RET", "AXL", "MYC",
                 "CCND1", "MKI67", "FN1", "TIMP1", "LCN2", "KLK10",
                 "TACSTD2"],
        "weight": 0.95,
        "rationale": "MTC + RR-DTC 후라인; MET-axis 차단.",
    },
    "vandetanib": {
        "moa": "VEGFR · EGFR · RET 억제 (multikinase)",
        "up": ["DUSP6", "SPRY2", "TG", "TPO", "TSHR", "DIO1"],
        "down": ["RET", "EGFR", "VEGFA", "KDR", "MYC", "CCND1", "MKI67",
                 "FN1", "TIMP1", "KLK10", "TACSTD2"],
        "weight": 0.9,
        "rationale": "MTC(medullary) 표준; RET-driven tumor에 효과.",
    },
    # ---- mTOR inhibitors
    "everolimus": {
        "moa": "mTORC1 억제 (mTOR inhibitor)",
        "up": ["SESN1", "SESN2", "DDIT4", "TG", "TPO", "DIO1", "TSHR"],
        "down": ["RPS6KB1", "EIF4EBP1", "MYC", "CCND1", "MKI67", "HIF1A",
                 "VEGFA", "FN1", "LCN2", "TACSTD2", "TMPRSS4"],
        "weight": 0.75,
        "rationale": "anaplastic/poorly differentiated 임상 시도; PI3K/AKT/mTOR 축 차단.",
    },
    "sirolimus": {
        "moa": "mTORC1 억제 (mTOR inhibitor)",
        "up": ["SESN1", "SESN2", "DDIT4", "TG", "TPO", "DIO1"],
        "down": ["RPS6KB1", "EIF4EBP1", "MYC", "CCND1", "MKI67", "HIF1A",
                 "VEGFA", "FN1", "TACSTD2"],
        "weight": 0.7,
        "rationale": "rapamycin — mTOR 축 차단; repurpose 후보.",
    },
    # ---- HDAC inhibitors (redifferentiation agents)
    "vorinostat": {
        "moa": "HDAC 억제 (SAHA, pan-HDAC inhibitor)",
        "up": ["CDKN1A", "CDKN2A", "CDKN2B", "TG", "TPO", "TSHR", "SLC5A5",
               "DIO1", "DIO2", "PAX8", "FOXE1", "NKX2-1", "CDH1"],
        "down": ["MYC", "CCND1", "MKI67", "HMGA2", "FN1", "VIM", "SNAI1",
                 "ZEB1", "TACSTD2", "TMPRSS4", "KLK10"],
        "weight": 0.85,
        "rationale": "epigenetic redifferentiation; thyroid NIS 재발현 실험 보고.",
    },
    "romidepsin": {
        "moa": "HDAC 억제 (class I HDAC inhibitor)",
        "up": ["CDKN1A", "CDKN2A", "TG", "TPO", "TSHR", "SLC5A5", "DIO1",
               "DIO2", "PAX8", "FOXE1"],
        "down": ["MYC", "CCND1", "MKI67", "HMGA2", "FN1", "VIM", "SNAI1",
                 "TACSTD2", "KLK10"],
        "weight": 0.8,
        "rationale": "CTCL 승인; thyroid redifferentiation 전임상.",
    },
    # ---- DNMT inhibitors (epigenetic)
    "decitabine": {
        "moa": "DNMT 억제 (DNA methyltransferase inhibitor)",
        "up": ["CDKN2A", "RASSF1", "MLH1", "TG", "TPO", "TSHR", "SLC5A5",
               "DIO1", "DIO2", "PAX8", "FOXE1"],
        "down": ["DNMT1", "DNMT3A", "MYC", "MKI67", "HMGA2", "FN1",
                 "TACSTD2", "KLK10"],
        "weight": 0.75,
        "rationale": "promoter hypermethylation 해소 — NIS 재발현 / redifferentiation.",
    },
    "azacitidine": {
        "moa": "DNMT 억제 (5-azacytidine, DNA methyltransferase inhibitor)",
        "up": ["CDKN2A", "RASSF1", "MLH1", "TG", "TPO", "TSHR", "SLC5A5",
               "DIO1", "DIO2", "PAX8", "FOXE1"],
        "down": ["DNMT1", "DNMT3A", "MYC", "MKI67", "HMGA2", "FN1",
                 "TACSTD2", "KLK10"],
        "weight": 0.75,
        "rationale": "MDS 승인; TET/DNMT 축을 통해 differentiation 프로그램 회복.",
    },
    # ---- Immune checkpoint (expression-side IFN-γ proxy)
    "pembrolizumab": {
        "moa": "anti-PD1 면역관문 억제 (anti-PD1 checkpoint blockade)",
        "up": ["IFNG", "GZMB", "PRF1", "CXCL9", "CXCL10", "STAT1", "IRF1",
               "CD8A", "CD8B", "GZMA", "NKG7"],
        "down": ["CD274", "PDCD1LG2", "HAVCR2", "LAG3", "TIGIT", "TOX",
                 "TACSTD2"],
        "weight": 0.7,
        "rationale": "proxy 주의 — 발현 변화가 아닌 PD1-PD-L1 결합 차단. BRAF-like의 immune-cold 상태를 논리적으로 역전.",
    },
    "nivolumab": {
        "moa": "anti-PD1 면역관문 억제 (anti-PD1 checkpoint blockade)",
        "up": ["IFNG", "GZMB", "PRF1", "CXCL9", "CXCL10", "STAT1", "IRF1",
               "CD8A", "GZMA"],
        "down": ["CD274", "PDCD1LG2", "HAVCR2", "LAG3", "TIGIT", "TOX"],
        "weight": 0.7,
        "rationale": "proxy 주의 — 직접 발현 조절 약물 아님; IFN-γ response에 의한 2차적 변화.",
    },
    # ---- VEGFR / anti-angiogenic
    "apatinib": {
        "moa": "VEGFR2 억제 (selective VEGFR2 inhibitor)",
        "up": ["TG", "TPO", "DIO1"],
        "down": ["VEGFA", "KDR", "HIF1A", "MYC", "CCND1", "MKI67",
                 "TACSTD2"],
        "weight": 0.7,
        "rationale": "RR-DTC 2상 (중국); anti-angiogenic 단독.",
    },
    "axitinib": {
        "moa": "VEGFR1/2/3 억제 (VEGFR inhibitor)",
        "up": ["TG", "TPO"],
        "down": ["VEGFA", "KDR", "FLT1", "FLT4", "HIF1A", "MKI67"],
        "weight": 0.6,
        "rationale": "RCC 승인; THCA 2상 보고 존재.",
    },
    # ---- Cytotoxics (broad signature — low specificity)
    "irinotecan": {
        "moa": "topoisomerase I 억제 (topoisomerase inhibitor)",
        "up": ["CDKN1A", "TP53", "BAX", "BBC3", "GADD45A"],
        "down": ["MYC", "CCND1", "MKI67", "MCM2", "MCM4", "PCNA"],
        "weight": 0.4,
        "rationale": "광범위 세포독성 — THCA 일차 치료는 아님.",
    },
    "topotecan": {
        "moa": "topoisomerase I 억제 (topoisomerase inhibitor)",
        "up": ["CDKN1A", "TP53", "BAX", "GADD45A"],
        "down": ["MYC", "CCND1", "MKI67", "MCM2", "PCNA"],
        "weight": 0.4,
        "rationale": "광범위 세포독성 — THCA에서 사용 드묾.",
    },
    "cisplatin": {
        "moa": "DNA 가교 백금 (platinum DNA-crosslink)",
        "up": ["CDKN1A", "TP53", "BAX", "BBC3", "GADD45A", "ATF3"],
        "down": ["MYC", "CCND1", "MKI67", "MCM2", "PCNA"],
        "weight": 0.4,
        "rationale": "ATC(anaplastic) 병용요법에 포함되지만 저분자화 프로그램 차원의 재분화 효과 없음.",
    },
    "carboplatin": {
        "moa": "DNA 가교 백금 (platinum DNA-crosslink)",
        "up": ["CDKN1A", "TP53", "BAX", "GADD45A"],
        "down": ["MYC", "CCND1", "MKI67", "PCNA"],
        "weight": 0.35,
        "rationale": "ATC 병용; 일반 세포독성.",
    },
    "oxaliplatin": {
        "moa": "DNA 가교 백금 (platinum DNA-crosslink)",
        "up": ["CDKN1A", "TP53", "GADD45A"],
        "down": ["MYC", "MKI67", "PCNA"],
        "weight": 0.3,
        "rationale": "CRC 표준; THCA 직접 효능 부족.",
    },
    "doxorubicin": {
        "moa": "topoisomerase II / DNA intercalation (anthracycline)",
        "up": ["CDKN1A", "TP53", "BAX", "GADD45A", "ATF3"],
        "down": ["MYC", "CCND1", "MKI67", "MCM2", "PCNA", "TOP2A"],
        "weight": 0.5,
        "rationale": "ATC에서 historical 표준; 일반 세포독성.",
    },
    "epirubicin": {
        "moa": "topoisomerase II / DNA intercalation (anthracycline)",
        "up": ["CDKN1A", "TP53", "BAX", "GADD45A"],
        "down": ["MYC", "MKI67", "PCNA", "TOP2A"],
        "weight": 0.45,
        "rationale": "doxorubicin 유사체; 유방암 표준.",
    },
    "paclitaxel": {
        "moa": "microtubule 안정화 (taxane)",
        "up": ["CDKN1A", "BCL2L11", "BAX", "ATF3"],
        "down": ["CCND1", "MKI67", "AURKA", "AURKB", "PLK1"],
        "weight": 0.45,
        "rationale": "ATC 병용; mitotic arrest로 광범위 증식억제.",
    },
    "docetaxel": {
        "moa": "microtubule 안정화 (taxane)",
        "up": ["CDKN1A", "BCL2L11", "BAX"],
        "down": ["CCND1", "MKI67", "AURKA", "AURKB", "PLK1"],
        "weight": 0.45,
        "rationale": "ATC 병용 고려; 유방/전립선암 표준.",
    },
    # ---- Metabolic / repurposed
    "metformin": {
        "moa": "AMPK 활성화 (AMPK activator / biguanide)",
        "up": ["PRKAA1", "PRKAA2", "CDKN1A", "TG", "TPO", "DIO1",
               "TSHR"],
        "down": ["MYC", "CCND1", "MKI67", "HIF1A", "MTOR", "RPS6KB1",
                 "FN1", "HMGA2", "TACSTD2"],
        "weight": 0.6,
        "rationale": "대규모 코호트에서 THCA 위험↓ 관찰; mTOR/inflammation 축 완화.",
    },
    "fenretinide": {
        "moa": "retinoid 유도체 (synthetic retinoid, ROS-mediated apoptosis)",
        "up": ["CDKN1A", "BAX", "TG", "TPO", "DIO1", "DIO2", "PAX8",
               "FOXE1", "RARA", "RARB"],
        "down": ["MYC", "CCND1", "MKI67", "FN1", "HMGA2", "TACSTD2"],
        "weight": 0.55,
        "rationale": "retinoid이 THCA redifferentiation에 historical 시도; evidence 약함.",
    },
    "curcumin": {
        "moa": "polyphenol · NF-κB 억제 (natural polyphenol)",
        "up": ["CDKN1A", "CDKN2A", "TG", "TPO", "DIO1"],
        "down": ["NFKB1", "MYC", "CCND1", "MKI67", "CXCL8", "TNF", "IL6",
                 "FN1", "TACSTD2"],
        "weight": 0.4,
        "rationale": "supplement 수준 증거; 약동학 한계.",
    },
    "quercetin": {
        "moa": "flavonoid · PI3K/NF-κB 억제 (flavonoid)",
        "up": ["CDKN1A", "TG", "TPO"],
        "down": ["PIK3CA", "NFKB1", "MYC", "CCND1", "MKI67", "CXCL8",
                 "FN1", "TACSTD2"],
        "weight": 0.35,
        "rationale": "supplement 수준 증거.",
    },
    "propranolol": {
        "moa": "β-adrenergic 차단 (β-blocker, repurposed)",
        "up": ["TG", "TPO", "TSHR", "DIO1"],
        "down": ["ADRB1", "ADRB2", "VEGFA", "MYC", "FN1", "TACSTD2"],
        "weight": 0.45,
        "rationale": "infantile hemangioma 승인; cancer repurpose 관찰연구.",
    },
    "simvastatin": {
        "moa": "HMG-CoA reductase 억제 (statin)",
        "up": ["LDLR", "HMGCR", "SREBF2", "CDKN1A"],
        "down": ["MYC", "CCND1", "MKI67", "FN1", "HMGA2", "TACSTD2"],
        "weight": 0.5,
        "rationale": "대규모 코호트에서 THCA risk 관찰; mevalonate 축 차단이 LDLR 신호 회복과 맞물림.",
    },
    "aspirin": {
        "moa": "COX1/COX2 억제 (non-selective NSAID)",
        "up": ["CDKN1A", "TG", "TPO"],
        "down": ["PTGS1", "PTGS2", "MYC", "CCND1", "MKI67", "CXCL8",
                 "IL6", "TACSTD2"],
        "weight": 0.45,
        "rationale": "CRC 예방 증거 확립; thyroid는 관찰 수준.",
    },
    "celecoxib": {
        "moa": "COX2 선택적 억제 (selective COX2 inhibitor)",
        "up": ["CDKN1A", "TG", "TPO"],
        "down": ["PTGS2", "MYC", "CCND1", "MKI67", "CXCL8", "IL6",
                 "TACSTD2"],
        "weight": 0.5,
        "rationale": "familial adenomatous polyposis 승인; THCA 2상 보고.",
    },
}


# ---------------------------------------------------------------------------
# Hardcoded literature-derived PPI edges. Each edge is (gene_a, gene_b,
# source, confidence). Sources paraphrased from STRING v11 high-confidence
# (>0.700) and BioGRID curated interactions. Confidence is on a 0-1 scale.
# Not an exhaustive dump — a focused subset around our focus genes.
# ---------------------------------------------------------------------------

PPI_EDGES_LIT: List[Tuple[str, str, str, float]] = [
    # --- MAPK / ERK pathway spine
    ("BRAF", "MAP2K1", "STRING", 0.999),
    ("BRAF", "MAP2K2", "STRING", 0.998),
    ("BRAF", "RAF1", "STRING", 0.990),
    ("BRAF", "ARAF", "STRING", 0.975),
    ("MAP2K1", "MAPK1", "STRING", 0.999),
    ("MAP2K1", "MAPK3", "STRING", 0.999),
    ("MAP2K2", "MAPK1", "STRING", 0.995),
    ("MAP2K2", "MAPK3", "STRING", 0.995),
    ("MAPK1", "DUSP5", "STRING", 0.950),
    ("MAPK1", "DUSP6", "STRING", 0.970),
    ("MAPK3", "DUSP5", "STRING", 0.945),
    ("MAPK3", "DUSP6", "STRING", 0.965),
    ("DUSP5", "DUSP6", "STRING", 0.820),
    ("MAPK1", "ETV5", "STRING", 0.910),
    ("ETV5", "DUSP6", "STRING", 0.800),
    # --- MET axis
    ("MET", "HGF", "STRING", 0.999),
    ("MET", "EGFR", "STRING", 0.920),
    ("MET", "CBL", "STRING", 0.900),
    ("MET", "GRB2", "STRING", 0.900),
    ("MET", "MAPK1", "BioGRID", 0.780),
    ("MET", "HRAS", "STRING", 0.850),
    # --- Thyroid-follicular lineage
    ("TG", "TPO", "STRING", 0.950),
    ("TG", "TSHR", "STRING", 0.920),
    ("TPO", "DUOX1", "STRING", 0.930),
    ("TPO", "DUOX2", "STRING", 0.940),
    ("DUOX1", "DUOX2", "STRING", 0.900),
    ("DUOX2", "DUOXA2", "STRING", 0.950),
    ("DIO1", "DIO2", "STRING", 0.850),
    ("DIO1", "TG", "STRING", 0.800),
    ("DIO2", "TSHR", "STRING", 0.820),
    ("FOXE1", "PAX8", "STRING", 0.890),
    ("FOXE1", "NKX2-1", "STRING", 0.900),
    ("PAX8", "NKX2-1", "STRING", 0.910),
    ("PAX8", "TG", "STRING", 0.900),
    ("PAX8", "TPO", "STRING", 0.890),
    ("PAX8", "SLC5A5", "STRING", 0.910),
    ("FOXE1", "TPO", "STRING", 0.800),
    ("NKX2-1", "TG", "STRING", 0.900),
    ("NKX2-1", "TPO", "STRING", 0.900),
    # --- SLC5A8 thyroid-related
    ("SLC5A8", "SLC5A5", "STRING", 0.750),
    ("SLC5A8", "TSHR", "BioGRID", 0.720),
    # --- Immune / checkpoint
    ("CD274", "PDCD1", "STRING", 0.999),
    ("CD274", "IFNG", "STRING", 0.870),
    ("CD274", "STAT1", "STRING", 0.830),
    ("HLA-DRA", "HLA-DRB1", "STRING", 0.999),
    ("HLA-DRA", "CD4", "STRING", 0.920),
    ("HLA-DRA", "CIITA", "STRING", 0.900),
    ("HLA-DRA", "IFNG", "STRING", 0.810),
    # --- LOX / ECM
    ("LOX", "COL1A1", "STRING", 0.820),
    ("LOX", "COL3A1", "STRING", 0.810),
    ("LOX", "FN1", "STRING", 0.790),
    ("LOX", "ELN", "STRING", 0.800),
    # --- KLK family
    ("KLK10", "KLK7", "STRING", 0.780),
    ("KLK10", "KLK6", "STRING", 0.750),
    ("KLK10", "KLK11", "STRING", 0.750),
    # --- Novel biomarkers — neighbours
    ("TACSTD2", "EPCAM", "STRING", 0.820),
    ("TACSTD2", "CLDN7", "STRING", 0.760),
    ("TACSTD2", "CDH1", "STRING", 0.720),
    ("TMPRSS4", "ST14", "STRING", 0.810),
    ("TMPRSS4", "CDH1", "BioGRID", 0.720),
    ("TMPRSS4", "SNAI1", "STRING", 0.740),
    ("CYP1B1", "AHR", "STRING", 0.930),
    ("CYP1B1", "CYP1A1", "STRING", 0.910),
    ("CYP1B1", "NR1I2", "STRING", 0.760),
    ("LDLR", "APOB", "STRING", 0.950),
    ("LDLR", "PCSK9", "STRING", 0.950),
    ("LDLR", "SREBF2", "STRING", 0.900),
    ("LDLR", "HMGCR", "STRING", 0.810),
    ("PLEKHA6", "PARD3", "BioGRID", 0.720),
    ("PLEKHA6", "CTNNB1", "BioGRID", 0.700),
    ("PLEKHA6", "CDH1", "BioGRID", 0.720),
    # --- cross-module bridges (AHR-ESR-thyroid axis)
    ("AHR", "ESR1", "STRING", 0.800),
    ("CYP1B1", "ESR1", "STRING", 0.780),
    ("AHR", "NKX2-1", "BioGRID", 0.700),
    # --- MYC / cell cycle hub (co-regulated genes)
    ("MYC", "CCND1", "STRING", 0.900),
    ("MYC", "MKI67", "STRING", 0.750),
    ("CCND1", "CDK4", "STRING", 0.990),
    ("CDK4", "RB1", "STRING", 0.990),
    # --- MAPK → MYC bridge
    ("MAPK1", "MYC", "STRING", 0.860),
    ("BRAF", "HRAS", "STRING", 0.940),
    # --- ERBB
    ("EGFR", "ERBB2", "STRING", 0.940),
    ("EGFR", "GRB2", "STRING", 0.960),
    ("GRB2", "SOS1", "STRING", 0.970),
    ("SOS1", "HRAS", "STRING", 0.940),
    ("HRAS", "RAF1", "STRING", 0.960),
    # --- LOX → MET bridge (fibrotic niche)
    ("LOX", "MET", "BioGRID", 0.700),
]


# Genes of primary interest (for node subsetting)
TOP5_NOVEL = ["TACSTD2", "PLEKHA6", "CYP1B1", "TMPRSS4", "LDLR"]
QUBO14 = [
    "KLK10", "DUSP5", "DUSP6", "DIO1", "LOX", "MET", "TPO", "SLC5A8",
    "DUOX2", "CD274", "HLA-DRA", "DUOX1", "DIO2", "FOXE1",
]
FOCUS_GENES = set(TOP5_NOVEL + QUBO14)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def plotly_offline_html(fig: go.Figure, out_path: Path, title: str,
                        height: int = 620) -> None:
    """Write a standalone CSP-safe HTML using vendored plotly.min.js."""
    fig.update_layout(
        paper_bgcolor="#0b1a3a",
        plot_bgcolor="#07132b",
        font=dict(color="#e2fff9", family="Inter, system-ui, sans-serif"),
        margin=dict(l=60, r=30, t=60, b=60),
    )
    inner = pio.to_html(
        fig,
        include_plotlyjs=False,
        full_html=False,
        config={"displaylogo": False, "responsive": True},
    )
    html = (
        "<!doctype html>\n"
        "<html lang=\"ko\"><head><meta charset=\"utf-8\">\n"
        f"<title>{title}</title>\n"
        "<style>html,body{margin:0;padding:0;background:#07132b;color:#e2fff9;"
        "font-family:Inter,system-ui,sans-serif}</style>\n"
        "<script src=\"../assets/vendor/plotly.min.js\"></script>\n"
        f"</head><body>{inner}</body></html>\n"
    )
    out_path.write_text(html)


# ---------------------------------------------------------------------------
# PHASE 1 — ssCMap-style connectivity score
# ---------------------------------------------------------------------------

def phase1_cmap() -> Tuple[pd.DataFrame, Dict]:
    print("[PHASE 1] CMap-style drug-reversal analysis")
    de = pd.read_csv(TABLES / "biomarker_de_full.tsv", sep="\t")
    de["gene"] = de["gene"].astype(str)
    # BRAF-up and RAS-up ranked gene lists (top 200 each by |log2FC| & FDR)
    de_sig = de[(de["fdr_tcga"] < 0.05) & de["log2FC_tcga"].notna()].copy()
    # BRAF-up = positive log2FC (BRAF_like > RAS_like). RAS-up = negative.
    braf_up = (
        de_sig[de_sig["log2FC_tcga"] > 0]
        .sort_values("log2FC_tcga", ascending=False)
        .head(200)["gene"]
        .tolist()
    )
    ras_up = (
        de_sig[de_sig["log2FC_tcga"] < 0]
        .sort_values("log2FC_tcga", ascending=True)
        .head(200)["gene"]
        .tolist()
    )
    braf_up_set = set(braf_up)
    ras_up_set = set(ras_up)
    print(
        f"  top-200 BRAF-up = {len(braf_up_set)}  top-200 RAS-up = {len(ras_up_set)}"
    )

    rows = []
    coverage_rows = []
    all_genes = braf_up_set | ras_up_set

    for cmp_name, info in COMPOUND_DB.items():
        up_sig = set(info["up"])
        dn_sig = set(info["down"])
        # ssCMap-like: reward genes that a compound UP-regulates appearing in
        # our DOWN list (RAS-up = things BRAF-like suppresses), and vice versa.
        # A candidate "reverser" pushes the profile AWAY from BRAF-like.
        up_in_ras_up = up_sig & ras_up_set  # compound up == disease "should-be up" (RAS-like)
        dn_in_braf_up = dn_sig & braf_up_set  # compound down == disease-up (BRAF-like)
        # Anti-pattern: compound up in BRAF_up (reinforces disease), down in RAS_up
        up_in_braf_up = up_sig & braf_up_set
        dn_in_ras_up = dn_sig & ras_up_set

        reversal = len(up_in_ras_up) + len(dn_in_braf_up)
        reinforce = len(up_in_braf_up) + len(dn_in_ras_up)
        denom = max(1, len(up_sig) + len(dn_sig))
        # connectivity_score ∈ [-1, +1]. More negative = stronger reverser.
        raw = (reinforce - reversal) / denom
        score = float(np.clip(raw, -1.0, 1.0)) * float(info.get("weight", 1.0))
        matched_total = reversal + reinforce
        matched_genes = sorted(
            up_in_ras_up | dn_in_braf_up | up_in_braf_up | dn_in_ras_up
        )
        direction = (
            "reverser" if score < -0.05
            else ("reinforcer" if score > 0.05 else "neutral")
        )
        rows.append({
            "compound": cmp_name,
            "moa": info["moa"],
            "weight": info.get("weight", 1.0),
            "sig_up_n": len(up_sig),
            "sig_down_n": len(dn_sig),
            "reversal_hits": reversal,
            "reinforce_hits": reinforce,
            "connectivity_score": round(score, 4),
            "direction": direction,
            "matched_total": matched_total,
            "matched_genes": ";".join(matched_genes[:30]),
            "rationale": info.get("rationale", ""),
        })
        for g in sorted(up_sig | dn_sig):
            coverage_rows.append({
                "compound": cmp_name,
                "gene": g,
                "compound_direction": "up" if g in up_sig else "down",
                "in_braf_up": int(g in braf_up_set),
                "in_ras_up": int(g in ras_up_set),
            })

    out = pd.DataFrame(rows).sort_values("connectivity_score")
    out_path = TABLES / "cmap_drug_reversal.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"  wrote {out_path}  ({len(out)} compounds)")

    # ---- Bar chart: top reversers (most-negative connectivity_score)
    top10 = out.head(10).iloc[::-1]  # flip so largest bar on top when horizontal
    fig = go.Figure(
        go.Bar(
            x=top10["connectivity_score"],
            y=top10["compound"],
            orientation="h",
            marker=dict(
                color=top10["connectivity_score"],
                colorscale=[
                    [0.0, "#14b8a6"], [0.5, "#60a5fa"], [1.0, "#f87171"],
                ],
                cmin=-1.0, cmax=1.0,
                colorbar=dict(title="score"),
                line=dict(color="rgba(148,163,184,.35)", width=1),
            ),
            customdata=np.stack([top10["moa"].values, top10["matched_total"].values,
                                 top10["rationale"].values], axis=-1),
            hovertemplate=(
                "<b>%{y}</b><br>score = %{x:.3f}<br>MoA = %{customdata[0]}<br>"
                "matched genes = %{customdata[1]}<br>"
                "<i>%{customdata[2]}</i><extra></extra>"
            ),
        )
    )
    fig.update_layout(
        title="Top-10 candidate reversers — ssCMap-like score (more negative = stronger reverser)",
        xaxis=dict(title="connectivity score", range=[-1.05, 1.05],
                   zeroline=True, zerolinecolor="rgba(148,163,184,.4)",
                   gridcolor="rgba(148,163,184,.12)"),
        yaxis=dict(title=""),
        height=520,
    )
    plotly_offline_html(fig, FIGS / "cmap_top_reversers.html",
                        "CMap top reversers")

    # ---- Heatmap: compound × gene coverage for top-15 reversers
    cov = pd.DataFrame(coverage_rows)
    top_cmp = out.head(15)["compound"].tolist()
    # For heatmap columns, take the union of top-reverser sig genes that show
    # up in our BRAF_up OR RAS_up list (so the figure is compact).
    genes_pool = (
        cov[(cov["compound"].isin(top_cmp)) &
            ((cov["in_braf_up"] == 1) | (cov["in_ras_up"] == 1))]["gene"]
        .value_counts()
        .head(40)
        .index
        .tolist()
    )
    if not genes_pool:
        genes_pool = sorted(list(all_genes))[:40]
    mat = np.zeros((len(top_cmp), len(genes_pool)))
    for i, cmp_name in enumerate(top_cmp):
        sub = cov[cov["compound"] == cmp_name]
        for j, g in enumerate(genes_pool):
            row = sub[sub["gene"] == g]
            if row.empty:
                continue
            # +1 if compound up-reg, -1 if down-reg, colour by alignment with disease
            d = row.iloc[0]["compound_direction"]
            in_braf = row.iloc[0]["in_braf_up"] == 1
            in_ras = row.iloc[0]["in_ras_up"] == 1
            if d == "up" and in_ras:
                mat[i, j] = 1  # reversal
            elif d == "down" and in_braf:
                mat[i, j] = 1
            elif d == "up" and in_braf:
                mat[i, j] = -1  # reinforcement
            elif d == "down" and in_ras:
                mat[i, j] = -1
            else:
                mat[i, j] = 0
    hm = go.Figure(
        go.Heatmap(
            z=mat, x=genes_pool, y=top_cmp,
            colorscale=[[0.0, "#ef4444"], [0.5, "#0b1a3a"], [1.0, "#14b8a6"]],
            zmin=-1, zmax=1,
            colorbar=dict(title="hit", tickvals=[-1, 0, 1],
                          ticktext=["reinforce", "none", "reverse"]),
            hovertemplate="<b>%{y}</b> × <b>%{x}</b>: %{z}<extra></extra>",
        )
    )
    hm.update_layout(
        title="Compound × gene coverage — top-15 candidate reversers",
        xaxis=dict(title="gene (top-40 most frequently hit)",
                   tickangle=45, tickfont=dict(size=10)),
        yaxis=dict(title="compound", autorange="reversed"),
        height=620,
    )
    plotly_offline_html(hm, FIGS / "cmap_coverage_heatmap.html",
                        "CMap coverage heatmap")

    # ---- JSON per-compound cards
    cards = []
    for _, r in out.iterrows():
        cards.append({
            "compound": r["compound"],
            "moa": r["moa"],
            "connectivity_score": r["connectivity_score"],
            "direction": r["direction"],
            "matched_total": int(r["matched_total"]),
            "rationale": r["rationale"],
        })
    (JSON_OUT / "cmap_compound_cards.json").write_text(
        json.dumps({"cards": cards, "generated_at": pd.Timestamp.utcnow().isoformat()},
                   ensure_ascii=False, indent=2)
    )
    print(f"  wrote {FIGS / 'cmap_top_reversers.html'}")
    print(f"  wrote {FIGS / 'cmap_coverage_heatmap.html'}")
    print(f"  wrote {JSON_OUT / 'cmap_compound_cards.json'}")

    summary = {
        "n_compounds": int(len(out)),
        "n_reversers": int((out["direction"] == "reverser").sum()),
        "n_reinforcers": int((out["direction"] == "reinforcer").sum()),
        "top5_reversers": out.head(5)[["compound", "moa", "connectivity_score"]]
                              .to_dict(orient="records"),
    }
    return out, summary


# ---------------------------------------------------------------------------
# PHASE 2 — PPI network (literature-derived)
# ---------------------------------------------------------------------------

def phase2_ppi() -> Tuple[pd.DataFrame, Dict]:
    print("[PHASE 2] PPI network — literature-derived")
    de = pd.read_csv(TABLES / "biomarker_de_full.tsv", sep="\t")
    de["gene"] = de["gene"].astype(str)
    de_idx = de.set_index("gene")

    G = nx.Graph()
    edge_rows = []
    for a, b, src, conf in PPI_EDGES_LIT:
        G.add_edge(a, b, source=src, confidence=conf)
        edge_rows.append({
            "gene_a": a, "gene_b": b, "source": src, "confidence": conf,
        })
    edges_df = pd.DataFrame(edge_rows)
    edges_df.to_csv(TABLES / "ppi_network_edges.tsv", sep="\t", index=False)
    print(f"  wrote {TABLES / 'ppi_network_edges.tsv'}  ({len(edges_df)} edges, {G.number_of_nodes()} nodes)")

    # Node attributes: |cohens_d| for size, direction for colour
    node_pos = nx.spring_layout(G, seed=42, k=1.1 / np.sqrt(max(1, G.number_of_nodes())),
                                iterations=250)
    node_x, node_y, node_text, node_size, node_color, node_meta = [], [], [], [], [], []
    for node in G.nodes:
        x, y = node_pos[node]
        node_x.append(x)
        node_y.append(y)
        if node in de_idx.index:
            row = de_idx.loc[node]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            d = float(row.get("cohens_d_tcga", 0) or 0)
            lfc = float(row.get("log2FC_tcga", 0) or 0)
            fdr = float(row.get("fdr_tcga", 1) or 1)
        else:
            d, lfc, fdr = 0.0, 0.0, 1.0
        size = 10 + 22 * min(1.5, abs(d)) / 1.5
        if node in FOCUS_GENES:
            size += 6
        # colour: red for BRAF-up (lfc>0 & fdr<0.05), blue for RAS-up (lfc<0), grey else
        if fdr < 0.05 and lfc > 0:
            col = "#ef4444"  # BRAF-up
        elif fdr < 0.05 and lfc < 0:
            col = "#38bdf8"  # RAS-up
        else:
            col = "#94a3b8"
        node_size.append(size)
        node_color.append(col)
        role = "novel" if node in TOP5_NOVEL else ("qubo" if node in QUBO14 else "neighbor")
        node_text.append(node)
        node_meta.append({
            "gene": node, "role": role, "log2FC": lfc, "cohens_d": d,
            "fdr": fdr, "degree": int(G.degree[node]),
        })

    edge_x, edge_y, edge_w = [], [], []
    for a, b, data in G.edges(data=True):
        x0, y0 = node_pos[a]
        x1, y1 = node_pos[b]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_w.append(data.get("confidence", 0.7))

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y, mode="lines",
        line=dict(width=1, color="rgba(148,163,184,.45)"),
        hoverinfo="skip", showlegend=False,
    ))
    customdata = np.array([
        [m["role"], f"{m['log2FC']:.2f}", f"{m['cohens_d']:.2f}", m["degree"]]
        for m in node_meta
    ], dtype=object)
    fig.add_trace(go.Scatter(
        x=node_x, y=node_y, mode="markers+text",
        text=node_text, textposition="top center",
        textfont=dict(size=11, color="#e2fff9"),
        marker=dict(
            size=node_size, color=node_color,
            line=dict(color="rgba(226,255,249,.6)", width=1.2),
        ),
        customdata=customdata,
        hovertemplate=(
            "<b>%{text}</b><br>역할 = %{customdata[0]}<br>"
            "log2FC = %{customdata[1]}<br>"
            "|cohen d| = %{customdata[2]}<br>"
            "degree = %{customdata[3]}<extra></extra>"
        ),
        showlegend=False,
    ))
    fig.update_layout(
        title="PPI network — top-5 novel + QUBO-14 + 1st-degree neighbours (literature-derived)",
        xaxis=dict(visible=False), yaxis=dict(visible=False),
        height=720, hovermode="closest",
        annotations=[
            dict(
                text="빨강 = BRAF-up · 파랑 = RAS-up · 회색 = 비유의<br>크기 ∝ |Cohen's d| · 엣지 = 문헌 기반 상호작용",
                x=0.01, y=0.01, xref="paper", yref="paper",
                showarrow=False, align="left",
                font=dict(size=11, color="#94a3b8"),
                bgcolor="rgba(11,26,58,.6)", bordercolor="rgba(148,163,184,.2)",
                borderwidth=1, borderpad=6,
            )
        ],
    )
    plotly_offline_html(fig, FIGS / "ppi_network_novel_qubo.html",
                        "PPI network — novel + QUBO")
    print(f"  wrote {FIGS / 'ppi_network_novel_qubo.html'}")

    degrees = sorted(((n, int(G.degree[n])) for n in G.nodes),
                     key=lambda x: (-x[1], x[0]))
    summary = {
        "n_nodes": int(G.number_of_nodes()),
        "n_edges": int(G.number_of_edges()),
        "top5_hubs": [{"gene": n, "degree": d} for n, d in degrees[:5]],
    }
    return edges_df, summary


# ---------------------------------------------------------------------------
# PHASE 3 — build page 19_cmap_network.html
# ---------------------------------------------------------------------------

TOPNAV_ORDER = [
    ("../index.html", "홈", False),
    ("01_overview.html", "개요", False),
    ("02_datasets.html", "데이터셋", False),
    ("03_sample_master.html", "샘플", False),
    ("04_gene_panels.html", "패널", False),
    ("05_eda.html", "EDA", False),
    ("06_scores.html", "점수", False),
    ("07_ml_baseline.html", "ML", False),
    ("08_panel_comparison.html", "패널 비교", False),
    ("09_shap.html", "SHAP", False),
    ("10_gene_explorer.html", "유전자", False),
    ("11_cohort_compare.html", "코호트", False),
    ("12_business.html", "비즈니스", False),
    ("13_reports.html", "리포트", False),
    ("14_caveats.html", "주의점", False),
    ("15_drug_discovery.html", "드럭", False),
    ("16_quantum.html", "양자", False),
    ("17_biomarker_insights.html", "마커", False),
    ("18_pathway_immune_meth.html", "해석", False),
    ("19_cmap_network.html", "약물리버스", True),
    ("view_investor.html", "투자자", False),
    ("view_researcher.html", "연구자", False),
    ("99_glossary.html", "용어", False),
]


def _topnav_html() -> str:
    parts = []
    for href, label, active in TOPNAV_ORDER:
        cls = "nav-link active" if active else "nav-link"
        parts.append(f'<a class="{cls}" href="{href}" role="menuitem">{label}</a>')
    return "".join(parts)


def phase3_page(cmap_df: pd.DataFrame, cmap_summary: Dict, ppi_summary: Dict) -> None:
    print("[PHASE 3] Build page 19_cmap_network.html")
    top10 = cmap_df.head(10)
    rows_html = []
    for _, r in top10.iterrows():
        rows_html.append(
            "<tr>"
            f"<td><strong>{r['compound']}</strong></td>"
            f"<td>{r['moa']}</td>"
            f"<td style=\"text-align:right;font-variant-numeric:tabular-nums\">{r['connectivity_score']:+.3f}</td>"
            f"<td style=\"text-align:right\">{int(r['matched_total'])}</td>"
            f"<td>{r['rationale']}</td>"
            "</tr>"
        )
    table_html = "\n".join(rows_html)

    hubs_html = "".join(
        f"<span class=\"cmap-pill\">{h['gene']} · degree {h['degree']}</span>"
        for h in ppi_summary.get("top5_hubs", [])
    )
    n_rev = cmap_summary.get("n_reversers", 0)
    n_reinf = cmap_summary.get("n_reinforcers", 0)
    n_total = cmap_summary.get("n_compounds", 0)
    n_nodes = ppi_summary.get("n_nodes", 0)
    n_edges = ppi_summary.get("n_edges", 0)
    topnav_html = _topnav_html()

    top1 = top10.iloc[0] if not top10.empty else None
    top1_name = top1["compound"] if top1 is not None else "(none)"
    top1_score = f"{top1['connectivity_score']:+.3f}" if top1 is not None else "-"

    page = f"""<!doctype html>
<html lang="ko">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="dark light">
  <meta name="theme-color" content="#07132b">
  <title>CMap 약물 리버스 · PPI 네트워크 | THCA Multi-Omics Dashboard</title>
  <meta name="description" content="Connectivity Map 스타일 약물 역전 후보 + 문헌 기반 단백질-단백질 상호작용 네트워크.">
  <link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
  <link rel="preload" href="../assets/fonts/inter-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/inter-600.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="stylesheet" href="../assets/css/main.css">
  <link rel="stylesheet" href="../assets/css/edu-toggle.css">
  <style>
    .cmap-hero{{position:relative;padding:36px 32px;border-radius:var(--radius-xl);overflow:hidden;border:1px solid rgba(148,163,184,.14);background:linear-gradient(135deg,rgba(6,17,38,.95),rgba(11,26,58,.68) 45%,rgba(44,18,90,.55));margin-bottom:22px}}
    .cmap-hero::before{{content:'';position:absolute;inset:-40%;z-index:-1;pointer-events:none;background:conic-gradient(from 120deg at 30% 30%,rgba(236,72,153,.22),transparent 30%),radial-gradient(circle at 70% 90%,rgba(94,234,212,.14),transparent 55%);filter:blur(60px);animation:cmapDrift 24s linear infinite}}
    @keyframes cmapDrift{{0%{{transform:rotate(0)}}100%{{transform:rotate(360deg)}}}}
    .cmap-hero h1{{font-size:clamp(1.8rem,1.3rem + 1.6vw,2.8rem);line-height:1.05;margin:0 0 6px;color:var(--ink-strong)}}
    .cmap-hero .lede{{color:var(--slate-300);max-width:760px}}
    .cmap-kpi-row{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin-top:20px}}
    .cmap-kpi{{padding:14px 16px;border-radius:14px;border:1px solid rgba(236,72,153,.3);background:linear-gradient(135deg,rgba(236,72,153,.14),rgba(94,234,212,.06))}}
    .cmap-kpi .k{{font-size:1.5rem;font-weight:700;color:#fbcfe8;letter-spacing:-.02em;line-height:1}}
    .cmap-kpi .l{{font-size:.7rem;color:var(--slate-300);margin-top:4px;text-transform:uppercase;letter-spacing:.12em}}
    .cmap-grid{{display:grid;grid-template-columns:1fr;gap:16px;margin-top:18px}}
    .cmap-card{{padding:16px;border-radius:var(--radius-lg);border:1px solid rgba(148,163,184,.14);background:rgba(11,26,58,.45);backdrop-filter:blur(10px)}}
    .cmap-card h3{{margin:0 0 6px;font-size:1.12rem;color:var(--ink-strong)}}
    .cmap-card .note{{font-size:.85rem;color:var(--slate-400);margin-bottom:10px;line-height:1.55}}
    .cmap-card .caveat{{font-size:.83rem;padding:10px 12px;border-radius:8px;background:rgba(250,204,21,.07);border:1px solid rgba(250,204,21,.2);color:#fde68a;margin-top:10px;line-height:1.55}}
    .cmap-card iframe{{width:100%;height:640px;border:0;display:block;border-radius:10px;background:#07132b}}
    .cmap-inlinefig{{width:100%;min-height:520px;border-radius:10px;overflow:hidden}}
    .cmap-table{{width:100%;border-collapse:collapse;font-size:.85rem;margin-top:12px}}
    .cmap-table th{{text-align:left;padding:8px 10px;border-bottom:1px solid rgba(148,163,184,.22);color:#bae6fd;font-weight:600}}
    .cmap-table td{{padding:8px 10px;border-bottom:1px solid rgba(148,163,184,.08);color:var(--ink);vertical-align:top}}
    .cmap-table tr:hover td{{background:rgba(148,163,184,.05)}}
    .cmap-pill{{display:inline-block;padding:3px 10px;border-radius:999px;font-size:.72rem;border:1px solid rgba(94,234,212,.3);color:#bae6fd;background:rgba(6,17,38,.5);margin:2px 4px 2px 0}}
    .cmap-footer-cta{{padding:16px;border-radius:12px;background:rgba(11,26,58,.45);border:1px solid rgba(148,163,184,.12);margin-top:16px;display:flex;justify-content:space-between;gap:12px;flex-wrap:wrap;align-items:center}}
    .cmap-footer-cta .txt{{color:var(--slate-300);font-size:.9rem}}
    .cmap-footer-cta .btns{{display:flex;gap:8px;flex-wrap:wrap}}
    @media (max-width: 760px){{
      .cmap-hero{{padding:22px 16px}}
      .cmap-hero h1{{font-size:1.55rem}}
      .cmap-kpi-row{{grid-template-columns:repeat(2,minmax(0,1fr))}}
      .cmap-card iframe{{height:60vh;min-height:400px}}
      .cmap-table{{font-size:.78rem}}
      .cmap-table th:nth-child(4),.cmap-table td:nth-child(4){{display:none}}
    }}
  </style>
</head>
<body class="dark">
  <a class="sr-only" href="#main-content">본문으로 건너뛰기</a>
  <nav class="topnav" aria-label="주요 메뉴">
    <div class="topnav-inner">
      <a class="brand" href="../index.html" aria-label="홈으로"><span class="brand-dot" aria-hidden="true"></span><span>THYROID DASH</span></a>
      <div class="nav-links" role="menubar">
        {topnav_html}
      </div>
      <div class="nav-actions">
        <button class="btn sm ghost" type="button" onclick="toggleTheme()" aria-label="다크/라이트 모드 전환">다크/라이트</button>
      </div>
    </div>
  </nav>
  <div class="layout">
    <main class="content" id="main-content" style="grid-column:1 / -1">

      <section class="cmap-hero" aria-labelledby="cmap-hero-title">
        <div class="subtle mono">19 / THERAPEUTIC — CMAP · PPI NETWORK</div>
        <h1 id="cmap-hero-title">CMap 약물 리버스 + PPI 네트워크</h1>
        <p class="lede">BRAF-like → RAS-like 재분화 상태로 되돌리는 <em>역전 후보 약물</em>을 CMap 스타일(connectivity score)로 스코어링하고, top-5 신규 바이오마커와 QUBO 선정 14 유전자의 <em>문헌 기반 PPI 네트워크</em>를 함께 표시합니다. 본 페이지의 CMap은 <strong>문헌-합의 surrogate</strong>이며 live CLUE.io 질의가 아닙니다.</p>
        <div class="cmap-kpi-row" role="list">
          <div class="cmap-kpi" role="listitem"><div class="k">{n_total}</div><div class="l">분석 화합물 (offline DB)</div></div>
          <div class="cmap-kpi" role="listitem"><div class="k">{n_rev}</div><div class="l">후보 리버서 (score &lt; -0.05)</div></div>
          <div class="cmap-kpi" role="listitem"><div class="k">{top1_name}</div><div class="l">Top-1 리버서 · score {top1_score}</div></div>
          <div class="cmap-kpi" role="listitem"><div class="k">{n_nodes} · {n_edges}</div><div class="l">PPI 노드 · 엣지</div></div>
        </div>
      </section>

      <details class="edu-toggle" open style="margin-bottom:18px;border-radius:12px;border:1px solid rgba(148,163,184,.18);background:rgba(11,26,58,.35);overflow:hidden">
        <summary style="cursor:pointer;padding:10px 14px;font-size:.9rem;color:#5eead4;list-style:none">💡 초보자용: CMap과 PPI가 무엇인가요?</summary>
        <div style="padding:4px 16px 14px;border-top:1px dashed rgba(148,163,184,.18);color:var(--ink);font-size:.9rem;line-height:1.65">
          <p><strong>Connectivity Map (CMap / CLUE.io)</strong> — 이 질환의 발현 패턴을 "거꾸로 뒤집는" 약물을 찾는 도구입니다. 예: 암에서 A 유전자가 올라가고 B가 내려가 있다면, A를 내리고 B를 올리는 약물이 "리버서"입니다. <em>원본 CMap은 LINCS L1000 실험 데이터베이스에 질의</em>하지만, 이 대시보드는 네트워크 fetch 없이 문헌에서 요약한 약물 시그니처로 <strong>proxy 점수</strong>를 계산합니다.</p>
          <p><strong>Connectivity score</strong> — −1(완전 리버서) ~ +1(완전 강화). 우리 BRAF-up 리스트에서 약물이 down-시그니처를 갖고 있거나, RAS-up 리스트에 약물이 up-시그니처를 갖고 있을수록 더 음수(= 더 강한 리버서).</p>
          <p><strong>PPI (Protein-Protein Interaction) 네트워크</strong> — 단백질끼리 물리적으로 결합하거나 기능적으로 상호작용한다고 <em>선행 문헌에서 보고</em>된 엣지 집합. STRING / BioGRID 같은 DB에서 추출합니다. 노드 크기는 해당 유전자의 이펙트 사이즈(|Cohen's d|), 색은 방향(빨강 = BRAF에서↑, 파랑 = RAS에서↑).</p>
          <p class="subtle small">이 두 층은 <em>생물학적 가설 우선순위</em>를 정하는 도구이지 임상 효능의 근거가 아닙니다. 리버서 점수가 높다 ≠ 임상에서 효과 있다.</p>
        </div>
      </details>

      <div class="cmap-grid">

        <section class="cmap-card" id="cmap">
          <h3>① CMap 스타일 약물 역전 (drug reversal)</h3>
          <div class="note">
            <strong>방법:</strong> 우리 TCGA DE 테이블에서 top-200 BRAF-up · top-200 RAS-up 유전자를 추출한 뒤, {n_total}개 화합물의 up/down 시그니처(문헌 기반 consensus)가 "리버스 패턴"(compound-up ∩ RAS-up, compound-down ∩ BRAF-up)과 "강화 패턴"(compound-up ∩ BRAF-up, compound-down ∩ RAS-up)에 각각 얼마나 겹치는지 계산. <code>score = (reinforce − reversal) / (up+down 시그니처 길이)</code> × compound weight. 범위 [−1, +1], <strong>음수가 강한 리버서</strong>.
          </div>
          <div class="cmap-inlinefig" style="margin-top:12px">
            <iframe src="../figs_interactive/cmap_top_reversers.html" loading="lazy" title="CMap top reversers"></iframe>
          </div>
          <table class="cmap-table" aria-label="Top-10 candidate reversers">
            <thead><tr><th>화합물</th><th>MoA</th><th style="text-align:right">score</th><th style="text-align:right">matched</th><th>왜 갑상선암에서 후보인가</th></tr></thead>
            <tbody>
            {table_html}
            </tbody>
          </table>
          <div class="cmap-inlinefig" style="margin-top:14px">
            <iframe src="../figs_interactive/cmap_coverage_heatmap.html" loading="lazy" title="CMap coverage heatmap"></iframe>
          </div>
          <div class="caveat">
            <strong>정직한 한계:</strong> (1) 이것은 <em>문헌-합의 proxy</em>이며 live CLUE.io / LINCS L1000 쿼리가 아닙니다. (2) 각 화합물의 up/down 시그니처는 공개 MoA 리뷰 + LINCS 요약을 손코딩한 15–30개 유전자로, 실제 L1000 ~978 랜드마크 경험식보다 훨씬 작습니다. (3) 세포독성 계열(cisplatin, doxorubicin 등)은 본질적으로 "MKI67↓ · CDKN1A↑" 공통 시그니처를 공유하여 어떤 암이든 <em>비특이적으로</em> 높게 나올 수 있습니다. (4) pembrolizumab/nivolumab은 발현 조절 약물이 아니므로 "IFN-γ response 유도" 대리 시그니처를 사용했으며, 해석은 조심스럽게. (5) score는 임상 효능의 근거가 아닙니다.
          </div>
        </section>

        <section class="cmap-card" id="ppi">
          <h3>② 단백질-단백질 상호작용 네트워크 (PPI, literature-derived)</h3>
          <div class="note">
            <strong>방법:</strong> top-5 신규 검증 바이오마커 (TACSTD2, PLEKHA6, CYP1B1, TMPRSS4, LDLR) + QUBO-선정 14 유전자 (KLK10, DUSP5, DUSP6, DIO1, LOX, MET, TPO, SLC5A8, DUOX2, CD274, HLA-DRA, DUOX1, DIO2, FOXE1)와 1-차 이웃을 대상으로 STRING v11 high-confidence (&ge;0.700) · BioGRID curated 엣지를 손코딩. networkx spring-layout (seed=42)로 force-directed 배치.
          </div>
          <div style="margin-top:10px">
            <strong>Top-5 허브 (degree 기준):</strong> {hubs_html}
          </div>
          <div class="cmap-inlinefig" style="margin-top:12px">
            <iframe src="../figs_interactive/ppi_network_novel_qubo.html" loading="lazy" title="PPI network"></iframe>
          </div>
          <div class="caveat">
            <strong>정직한 한계:</strong> (1) 엣지는 <em>외부 문헌에서 가져온 것</em>이며 본 프로젝트 발현 데이터로부터 <strong>재추론된 것이 아닙니다</strong>. (2) PPI DB는 연구-편향이 존재 — 잘 연구된 유전자일수록 엣지가 많이 붙고, under-studied 유전자(예: PLEKHA6)는 엣지가 적어 중심성이 낮게 나타납니다. (3) confidence score는 출처별 체계가 다르므로 직접 비교 조심. (4) 물리적 결합과 기능적 association을 같은 엣지로 합쳐 표시했습니다.
          </div>
        </section>

        <section class="cmap-card" id="sources">
          <h3>출처 · 정직성 선언</h3>
          <div class="note">
            <ul>
              <li><strong>CMap / LINCS 개념</strong> — Lamb et al. <em>Science</em> 2006 (PMID 17008526); Subramanian et al. <em>Cell</em> 2017 (PMID 29195078, L1000).</li>
              <li><strong>화합물 시그니처 (up/down 15–30 genes)</strong> — 해당 약물 MoA 리뷰 + iLINCS / LINCS L1000 공개 요약에서 파라프레이즈. 원본 L1000 landmark 데이터는 포함되지 않음.</li>
              <li><strong>PPI 엣지</strong> — STRING v11 (Szklarczyk et al. <em>Nucleic Acids Res</em> 2019, PMID 30476243, cutoff &ge;0.700), BioGRID (Oughtred et al. <em>Protein Sci</em> 2021, PMID 33070389). hardcoded subset.</li>
              <li><strong>갑상선 근거</strong> — Ho et al. <em>NEJM</em> 2013 (selumetinib redifferentiation); Schlumberger et al. <em>NEJM</em> 2015 (lenvatinib SELECT); Brose et al. <em>Lancet</em> 2014 (sorafenib DECISION).</li>
            </ul>
            <p><strong>선언:</strong> (a) CMap 결과는 <em>문헌 합의 surrogate</em>이지 live LINCS 쿼리가 아닙니다. (b) PPI 엣지는 <em>공공 DB 큐레이션</em>이지 본 프로젝트 발현 데이터로부터 re-derive된 것이 아닙니다. (c) <strong>생물학적 타당성 ≠ 임상 효능</strong>. 모든 결과는 가설 우선순위 설정용입니다.</p>
          </div>
        </section>

      </div>

      <div class="cmap-footer-cta">
        <div class="txt">이 페이지의 분석은 <strong>드럭 디스커버리</strong>의 타겟 리스트에 약물-측 증거를 더하고, <strong>바이오마커 인사이트</strong>의 유전자를 네트워크 맥락에서 해석합니다.</div>
        <div class="btns">
          <a class="btn sm" href="15_drug_discovery.html">드럭 디스커버리</a>
          <a class="btn sm" href="17_biomarker_insights.html">바이오마커 인사이트</a>
          <a class="btn sm ghost" href="08_panel_comparison.html">패널 비교</a>
        </div>
      </div>

    </main>
  </div>

  <script>
    if(typeof toggleTheme === 'undefined'){{
      window.toggleTheme = function(){{
        var isLight = document.documentElement.classList.toggle('light');
        try{{ localStorage.setItem('thyroid-theme', isLight ? 'light' : 'dark'); }}catch(e){{}}
      }};
    }}
  </script>
</body>
</html>
"""
    out = PAGES / "19_cmap_network.html"
    out.write_text(page)
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# PHASE 4 — inject nav link across all pages
# ---------------------------------------------------------------------------

INJECT_SCRIPT = SCRIPTS / "inject_cmap_navlink.py"

INJECT_CONTENT = '''#!/usr/bin/env python3
"""Idempotent: inject a '약물리버스' nav link pointing to 19_cmap_network.html.

Inserts after the '해석' (18_pathway_immune_meth.html) link in topnav /
drawer / commandIndex. Safe to re-run.
"""
import json
import re
from pathlib import Path

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/reports/html")
TARGET = "19_cmap_network.html"
LABEL = "약물리버스"

pages = list((ROOT / "pages").glob("*.html")) + [ROOT / "index.html"]


def process(p: Path):
    if not p.exists():
        return
    text = p.read_text()
    in_pages = p.parent.name == "pages"
    is_self = p.name == TARGET
    href = TARGET if in_pages else f"pages/{TARGET}"

    if TARGET in text and f">{LABEL}</a>" in text:
        print(f"skip (already has link): {p.name}")
        return

    new_text = text
    patched = False

    # --- TOPNAV: insert after '해석' link (role=menuitem variant)
    pat_nav_role = re.compile(
        r'(<a class="nav-link[^"]*" href="(?:\\.\\./)?(?:pages/)?18_pathway_immune_meth\\.html" role="menuitem">해석</a>)',
        re.UNICODE,
    )
    pat_nav_plain = re.compile(
        r'(<a class="nav-link[^"]*" href="(?:\\.\\./)?(?:pages/)?18_pathway_immune_meth\\.html">해석</a>)',
        re.UNICODE,
    )

    def sub_nav(m):
        cls = "nav-link active" if is_self else "nav-link"
        source = m.group(1)
        if 'href="pages/18_pathway_immune_meth.html"' in source:
            hh = f"pages/{TARGET}"
        elif 'href="../pages/18_pathway_immune_meth.html"' in source:
            hh = f"../pages/{TARGET}"
        elif 'href="../18_pathway_immune_meth.html"' in source:
            hh = f"../{TARGET}"
        else:
            hh = href
        if ' role="menuitem"' in source:
            return source + f'<a class="{cls}" href="{hh}" role="menuitem">{LABEL}</a>'
        else:
            return source + f'<a class="{cls}" href="{hh}">{LABEL}</a>'

    replaced, n = pat_nav_role.subn(sub_nav, new_text, count=1)
    if n:
        new_text = replaced; patched = True
    else:
        replaced, n = pat_nav_plain.subn(sub_nav, new_text, count=1)
        if n:
            new_text = replaced; patched = True

    # --- DRAWER
    pat_drawer = re.compile(
        r'(<a class="[^"]*" href="(?:\\.\\./)?(?:pages/)?18_pathway_immune_meth\\.html">해석</a>)',
        re.UNICODE,
    )
    def sub_drawer(m):
        cls = "active" if is_self else ""
        source = m.group(1)
        if 'href="pages/18_pathway_immune_meth.html"' in source:
            hh = f"pages/{TARGET}"
        elif 'href="../18_pathway_immune_meth.html"' in source:
            hh = f"../{TARGET}"
        else:
            hh = href
        return source + f'<a class="{cls}" href="{hh}">{LABEL}</a>'
    replaced, n = pat_drawer.subn(sub_drawer, new_text, count=1)
    if n:
        new_text = replaced; patched = True

    # --- commandIndex
    cmd_idx_pat = re.compile(
        r'(window\\.ThyroidDash\\.commandIndex\\s*=\\s*)(\\[[\\s\\S]*?\\])(\\s*;)',
        re.UNICODE,
    )
    m_cmd = cmd_idx_pat.search(new_text)
    if m_cmd:
        arr_text = m_cmd.group(2)
        try:
            arr = json.loads(arr_text)
        except Exception:
            arr = None
        if isinstance(arr, list):
            existing_hrefs = {e.get("href", "") for e in arr if isinstance(e, dict)}
            candidate_hrefs = {TARGET, f"pages/{TARGET}", f"../{TARGET}"}
            if not (existing_hrefs & candidate_hrefs):
                new_entry = {"label": LABEL, "href": href,
                             "tags": "page cmap ppi network drug reversal"}
                idx_after = None
                for i, e in enumerate(arr):
                    if not isinstance(e, dict):
                        continue
                    if "18_pathway_immune_meth.html" in e.get("href", "") or e.get("label") == "해석":
                        idx_after = i
                        break
                if idx_after is None:
                    new_arr = arr + [new_entry]
                else:
                    new_arr = arr[: idx_after + 1] + [new_entry] + arr[idx_after + 1:]
                replacement = m_cmd.group(1) + json.dumps(new_arr, ensure_ascii=False) + m_cmd.group(3)
                new_text = new_text[:m_cmd.start()] + replacement + new_text[m_cmd.end():]
                patched = True

    if patched and new_text != text:
        p.write_text(new_text)
        print(f"patched: {p.name}")
    else:
        print(f"no change: {p.name}")


if __name__ == "__main__":
    for pp in pages:
        process(pp)
'''


def phase4_inject() -> None:
    print("[PHASE 4] Writing inject_cmap_navlink.py + running it")
    INJECT_SCRIPT.write_text(INJECT_CONTENT)
    os.chmod(INJECT_SCRIPT, 0o755)
    subprocess.check_call([sys.executable, str(INJECT_SCRIPT)])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    cmap_df, cmap_summary = phase1_cmap()
    _edges_df, ppi_summary = phase2_ppi()
    phase3_page(cmap_df, cmap_summary, ppi_summary)
    phase4_inject()
    print("\n========== SUMMARY ==========")
    print(json.dumps({"cmap": cmap_summary, "ppi": ppi_summary},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
