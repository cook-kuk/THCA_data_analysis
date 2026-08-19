#!/usr/bin/env python3
"""Scout public patent/company sources for an industrial neoantigen benchmark.

The output is an evidence inventory, not a training set. Public industrial
sources should be used only after freezing the candidate algorithm, as an
out-of-distribution stress test. Confidential, leaked, paywalled or
terms-violating sources are intentionally excluded.
"""

from __future__ import annotations

import hashlib
import html
import json
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from urllib.parse import urlparse

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "project/results/public_industrial_neoantigen_benchmark_scout_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")

AA_RE = re.compile(r"(?<![A-Z])([ACDEFGHIKLMNPQRSTVWY]{8,35})(?![A-Z])")
HLA_CLASSI_RE = re.compile(
    r"(?:HLA[-\s]?)?([ABC])\*?(\d{2})[:]?(\d{2})|(?:HLA[-\s]?)?([ABC])(\d{4})",
    re.IGNORECASE,
)
HLA_CLASSII_RE = re.compile(
    r"(?:HLA[-\s]?)?((?:DRB[1-5]|DQA1|DQB1|DPA1|DPB1))\*?(\d{2,3})[:]?(\d{2})|"
    r"(?:HLA[-\s]?)?((?:DRB[1-5]|DQA1|DQB1|DPA1|DPB1))[-_\s]?(\d{2,3})[-_:]?(\d{2})",
    re.IGNORECASE,
)
ASSAY_RE = re.compile(
    r"ELISPOT|IFN[\s-]*γ|IFN-gamma|tetramer|multimer|T[-\s]?cell|CD8|CD4|response|immunogenic|killing|DTH",
    re.IGNORECASE,
)
LABEL_RE = re.compile(r"\+\+\+|\+\+|\+|-\s|No response|RESPONSE|responder|non[-\s]?responder", re.IGNORECASE)
METHOD_RE = re.compile(
    r"HLA binding|T cell recognition|cleavage|presentation|flanking|pseudo-sequence|expression|RNA|WES|RNA-seq|"
    r"synthetic long peptide|SLP|mRNA|viral vector|machine[-\s]?learn|deep learning|AI|mutation position|rare HLA",
    re.IGNORECASE,
)
AA_STOPWORDS = {
    "RESEARCH",
    "PRACTICE",
    "CLINICAL",
    "PATIENTS",
    "ADAPTIVE",
    "RESPONSE",
    "RESPONSES",
    "TREATMENT",
    "THERAPY",
    "CANDIDATE",
    "CANDIDATES",
    "SCREENING",
    "MUTATION",
    "MUTATIONS",
    "SEQUENCE",
    "SEQUENCES",
    "PRESENTED",
    "PRESENTING",
    "ACTIVITY",
    "ANALYSIS",
    "IMMUNITY",
    "MELANOMA",
    "MEDICINE",
    "PRECISION",
    "SCIENTIFIC",
}


@dataclass(frozen=True)
class Source:
    source_id: str
    organization: str
    source_class: str
    expected_tier: str
    url: str
    rationale: str


SOURCES = [
    Source(
        "THERAGEN_DEEPOMICS_NEO_EN",
        "Theragen Bio",
        "company_platform",
        "C",
        "https://www.theragenbio.com/eng/rnd/deep.php",
        "DEEPOMICS NEO public platform: HLA I/II, rare HLA, T-cell activity and SLP design claims.",
    ),
    Source(
        "THERAGEN_DEEPOMICS_NEO_PR",
        "Theragen Bio",
        "company_press",
        "C",
        "https://theragenbio.com/eng/promotion/data.php/1000?category=66&code=data_eng&idx=1591&ptype=view",
        "Public patent announcement describing SLP immunogenicity, HLA binding, T-cell recognition and cleavage features.",
    ),
    Source(
        "THERAGENETEX_NEO_PATENT_KR",
        "Theragen Bio",
        "company_press",
        "C",
        "https://www.theragenetex.com/ko/news_detail/278/990",
        "Korean public announcement of neoantigen prediction patent using peptide and HLA allele sequences.",
    ),
    Source(
        "THERAGEN_JUSTIA_SLP_20240071564",
        "Theragen Bio",
        "patent",
        "C",
        "https://patents.justia.com/patent/20240071564",
        "Patent text for SLP immunogenicity prediction, cleavage features, binding and T-cell activity.",
    ),
    Source(
        "US12567480_GAZETTE_CLASSI_IM",
        "Unknown/Patent assignee",
        "patent",
        "C",
        "https://patentsgazette.uspto.gov/week09/OG/html/1544-1/US12567480-20260303.html",
        "Patent claim for peptide flanking region + HLA pseudo-sequence presentation/binding model.",
    ),
    Source(
        "PGV001_PIPELINE_FRONTIERS",
        "Mount Sinai PGV-001",
        "trial_pipeline",
        "B",
        "https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2017.01807/full",
        "Personalized vaccine computational pipeline with WES/RNA and ranked SLP vaccine contents.",
    ),
    Source(
        "EVAXION_AACR2026_EVX01",
        "Evaxion",
        "company_poster_page",
        "B",
        "https://evaxion.ai/scientific-posters/aacr-2026-evx-01/",
        "Final immunogenicity readout for AI-selected EVX-01 neoantigen vaccine.",
    ),
    Source(
        "EVAXION_EVX01_ONCOIMMUNOLOGY_PMC",
        "Evaxion",
        "clinical_publication",
        "A",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC11116868/",
        "EVX-01 clinical paper with immunogenicity, pMHC multimers and example vaccine epitope sequence.",
    ),
    Source(
        "EVAXION_ASCO2024_POSTER_PDF",
        "Evaxion",
        "company_poster_pdf",
        "A",
        "https://www.evaxion-biotech.com/media/nzblm4qr/asco2024-poster-evx01p2_final.pdf",
        "ASCO 2024 poster for AI-designed EVX-01 immunogenicity readout.",
    ),
    Source(
        "EVAXION_ASCO2024_POSTER_PDF_ALT",
        "Evaxion",
        "company_poster_pdf",
        "A",
        "https://evaxion.ai/wp-content/uploads/2024/04/asco2024-poster-evx01p2_final.pdf",
        "Alternate public ASCO 2024 poster URL for AI-designed EVX-01 immunogenicity readout.",
    ),
    Source(
        "BIONTECH_PD_INDEX_NATURE2023",
        "BioNTech/Genentech",
        "clinical_publication",
        "A",
        "https://www.nature.com/articles/s41586-023-06063-y",
        "Autogene cevumeran PDAC study with ELISpot-identified immunodominant neoantigens and TCR validation.",
    ),
    Source(
        "BIONTECH_ADVANCED_SOLID_NATUREMED2024",
        "BioNTech/Genentech",
        "clinical_publication",
        "A",
        "https://www.nature.com/articles/s41591-024-03334-7",
        "Autogene cevumeran phase 1 trial with ex vivo/IVS ELISpot and pMHC multimer validation.",
    ),
    Source(
        "RCC_NEOANTIGEN_NATURE2024",
        "Academic/Personalized vaccine",
        "clinical_publication",
        "A",
        "https://www.nature.com/articles/s41586-024-08507-5",
        "Renal cell carcinoma neoantigen vaccine study with deconvoluted IFNγ ELISpot assays.",
    ),
    Source(
        "MODERNA_MERCK_V940_PATSNAP",
        "Moderna/Merck",
        "industry_analysis",
        "C",
        "https://www.patsnap.com/resources/blog/articles/mrna-cancer-vaccine-pipeline-v940-bnt111-neoantigens/",
        "Public industrial summary of V940/mRNA-4157 and patent landscape.",
    ),
    Source(
        "NYKODE_VB10_NEO_PATENT_PR",
        "Nykode",
        "company_press_pdf",
        "C",
        "https://nykode.com/wp-content/uploads/2024/08/240814-Nykode-PR-VB10.NEO-patent_FINAL.pdf",
        "Public announcement of VB10.NEO individualized neoantigen vaccine patent.",
    ),
    Source(
        "NYKODE_VB10_NEO_PLATFORM",
        "Nykode",
        "company_platform",
        "C",
        "https://nykode.com/vb10-neo-a-clinically-validated-individualized-vaccine-platform/",
        "VB10.NEO clinical platform page.",
    ),
    Source(
        "NYKODE_VB10_JOURNEY",
        "Nykode",
        "company_platform",
        "C",
        "https://nykode.com/vb10-neo-the-journey-of-a-clinically-validatedindividualized-cancer-therapy-from-tumorbiopsy-to-vaccine-administration/",
        "VB10.NEO process page from biopsy to vaccine administration.",
    ),
    Source(
        "NYKODE_VB10_ASCO2025_PDF",
        "Nykode",
        "company_poster_pdf",
        "B",
        "https://nykode.com/wp-content/uploads/2025/06/N02_ASCO_poster_2025_v4.pdf",
        "ASCO 2025 VB10.NEO poster with IVS/ex vivo ELISpot immunogenicity rates.",
    ),
    Source(
        "TRANSGENE_TG4050_PRODUCT",
        "Transgene/NEC",
        "company_platform",
        "B",
        "https://www.transgene.com/tg4050-1/",
        "TG4050 product page: NEC prediction system and myvac individualized vaccine.",
    ),
    Source(
        "TRANSGENE_NEC_PRESS",
        "Transgene/NEC",
        "company_press",
        "B",
        "https://neclab.eu/about-us/press-releases/detail/transgene-and-nec-demonstrate-high-accuracy-of-ai-based-neoantigen-prediction-for-the-design-of-individualized-cancer-vaccine-tg4050",
        "Press release for high accuracy AI-based neoantigen prediction for TG4050.",
    ),
    Source(
        "TRANSGENE_TG4050_ASCO2023_PDF",
        "Transgene/NEC",
        "company_poster_pdf",
        "A",
        "https://www.transgene.com/wp-content/uploads/Transgene_ASCO2023_Safety_Immunogenicity_TG4050.pdf",
        "ASCO 2023 TG4050 poster with peptide sequences, HLA prediction and ELISpot/tetramer response snippets.",
    ),
    Source(
        "TRANSGENE_TG4050_ASCO2022_PDF",
        "Transgene/NEC",
        "company_poster_pdf",
        "B",
        "https://www.transgene.com/wp-content/uploads/Transgene_NEC_TG4050_ASCO2022_poster.pdf",
        "ASCO 2022 TG4050 phase 1 poster for ovarian/head-neck personalized vaccine.",
    ),
    Source(
        "TRANSGENE_TG4050_AACR2020_PDF",
        "Transgene/NEC",
        "company_poster_pdf",
        "A",
        "https://www.transgene.com/wp-content/uploads/posters/AACR2020-myvac.pdf",
        "TG4050 prediction performance poster with deconvoluted peptide immunogenicity claims.",
    ),
    Source(
        "TRANSGENE_TG4050_AACR2022_PDF",
        "Transgene/NEC",
        "company_poster_pdf",
        "B",
        "https://www.transgene.com/wp-content/uploads/Transgene_TG4050_AACR2022_poster.pdf",
        "TG4050 phase I poster with peptide-level ELISpot context.",
    ),
    Source(
        "TRANSGENE_TG4050_MEDRXIV2026",
        "Transgene/NEC",
        "preprint",
        "A",
        "https://www.medrxiv.org/content/10.64898/2026.01.06.25342687v1.full",
        "TG4050 randomized phase I preprint with immunogenicity and TCR response details.",
    ),
    Source(
        "TRANSGENE_TG4050_SITC2024_PDF",
        "Transgene/NEC",
        "company_poster_pdf",
        "A",
        "https://www.transgene.fr/wp-content/uploads/20241108_TG4050_SITC2024_Poster.pdf",
        "SITC 2024 TG4050 poster with sustained neoantigen-specific IFN-gamma ELISpot responses.",
    ),
    Source(
        "TRANSGENE_TG4050_SITC2025_PDF",
        "Transgene/NEC",
        "company_poster_pdf",
        "A",
        "https://www.transgene.com/wp-content/uploads/20251107_TG4050_SITC2025_Poster.pdf",
        "SITC 2025 TG4050 poster with ELISpot/tetramer response summary.",
    ),
    Source(
        "GENEOS_GNOSPV02_ASCO2023_TCR_PDF",
        "Geneos",
        "company_poster_pdf",
        "B",
        "https://www.geneostx.com/wp-content/uploads/2024/04/8.6.2.5.1.1_2023_ASCO_Poster_TCR_work_24pts_FINAL.pdf",
        "ASCO 2023 GNOS-PV02 poster with neoantigen-specific CD8/CD4 T-cell response context.",
    ),
    Source(
        "GENEOS_GNOSPV02_AACR2024_GT30_PDF",
        "Geneos",
        "company_poster_pdf",
        "B",
        "https://www.geneostx.com/wp-content/uploads/2024/04/AACR_2024-GT30-trial-poster.pdf",
        "AACR 2024 GNOS-PV02 poster with IFNg ELISpot response rates.",
    ),
    Source(
        "GEN009_JITC2020_ABSTRACT",
        "Genocea",
        "conference_abstract",
        "B",
        "https://jitc.bmj.com/content/8/Suppl_3/A251.1",
        "SITC/JITC 2020 GEN-009 abstract with ex vivo/in vitro fluorospot immunogenicity.",
    ),
    Source(
        "GEN009_JITC2021_ABSTRACT",
        "Genocea",
        "conference_abstract",
        "B",
        "https://jitc.bmj.com/content/9/Suppl_2/A551",
        "SITC/JITC 2021 GEN-009 abstract with durable immune responses.",
    ),
    Source(
        "ELICIO_SITC2024_KRAS_PDF",
        "Elicio",
        "company_poster_pdf",
        "C",
        "https://elicio.com/wp-content/uploads/2024/11/2024-SITC-poster-FINAL.pdf",
        "SITC 2024 ELI-002 7P KRAS vaccine poster; off-the-shelf mutant KRAS, useful as public antigen control rather than personalized neoantigen.",
    ),
    Source(
        "GRITSTONE_EP4576103_PUBCHEM",
        "Gritstone Bio",
        "patent",
        "C",
        "https://pubchem.ncbi.nlm.nih.gov/patent/EP-4576103-A2",
        "Gritstone neoantigen hotspot/presentation likelihood patent summary.",
    ),
    Source(
        "GRITSTONE_VIRAL_DELIVERY_USPTOREPORT",
        "Gritstone Bio",
        "patent",
        "C",
        "https://uspto.report/patent/app/20200010849",
        "Viral delivery of neoantigens patent application with selection/presentation module.",
    ),
    Source(
        "GRITSTONE_ALPHAVIRUS_PUBCHEM",
        "Gritstone Bio",
        "patent",
        "C",
        "https://pubchem.ncbi.nlm.nih.gov/patent/US-11504421-B2",
        "Alphavirus neoantigen vectors patent summary.",
    ),
    Source(
        "GRITSTONE_EP3389630_GREYB",
        "Gritstone Bio",
        "patent_analysis",
        "C",
        "https://ipverse.greyb.com/patents/EP3389630",
        "Neoantigen identification/manufacture/use patent analysis emphasizing MS and presentation modeling.",
    ),
    Source(
        "JUSTIA_NEOANTIGENS_METHODS_20230405102",
        "Unknown/Patent assignee",
        "patent",
        "B",
        "https://patents.justia.com/patent/20230405102",
        "Patent with many disclosed neoantigenic peptide/protein sequence examples.",
    ),
    Source(
        "GOOGLE_PATENTS_WO2023089203",
        "Unknown/Patent assignee",
        "patent",
        "C",
        "https://patents.google.com/patent/WO2023089203A1/en",
        "Neoantigen immunogenicity prediction patent with 25mer mut-seq, 8-12mer neo-peps and feature tools.",
    ),
    Source(
        "MODERNA_WO2017020026_CONCATEMERIC_RNA",
        "Moderna",
        "patent",
        "C",
        "https://patents.google.com/patent/WO2017020026A1/en",
        "Moderna concatemeric peptide epitope RNA patent cited in neoepitope vaccine families.",
    ),
    Source(
        "BIONTECH_HRP20201671_NEOANTIGEN_USEFULNESS",
        "BioNTech/TRON",
        "patent",
        "C",
        "https://patents.google.com/patent/HRP20201671T1/en",
        "BioNTech/TRON method for predicting usefulness of neoantigens for immunotherapy.",
    ),
    Source(
        "BIONTECH_WO2020252039_NEOANTIGEN_COMPOSITIONS",
        "BioNTech",
        "patent",
        "C",
        "https://patents.google.com/patent/WO2020252039A1/en",
        "BioNTech neoantigen compositions and immunotherapeutic polypeptides.",
    ),
    Source(
        "BIONTECH_WO2017118702_NEOEPITOPE_RNA",
        "BioNTech",
        "patent",
        "C",
        "https://patents.google.com/patent/WO2017118702A1/en",
        "BioNTech neoepitope RNA cancer vaccine patent family.",
    ),
    Source(
        "GRITSTONE_IL273030_NEOANTIGEN_T_CELL_THERAPY",
        "Gritstone Bio",
        "patent",
        "C",
        "https://patents.google.com/patent/IL273030A/en",
        "Gritstone neoantigen identification for T-cell therapy patent.",
    ),
    Source(
        "GRITSTONE_AU2016369519_MANUFACTURE_USE",
        "Gritstone Bio",
        "patent",
        "C",
        "https://patents.google.com/patent/AU2016369519A1/en",
        "Gritstone neoantigen identification, manufacture and use patent.",
    ),
    Source(
        "GRITSTONE_EP3694532_HOTSPOTS",
        "Gritstone Bio",
        "patent",
        "C",
        "https://patents.google.com/patent/EP3694532A4/en",
        "Gritstone neoantigen identification using hotspots patent.",
    ),
    Source(
        "GRITSTONE_WO2018208856_ALPHAVIRUS",
        "Gritstone Bio",
        "patent",
        "C",
        "https://patents.google.com/patent/WO2018208856A1/en",
        "Gritstone alphavirus neoantigen vectors patent.",
    ),
    Source(
        "NYKODE_EP3856957_SELECTING_NEOEPITOPES",
        "Nykode",
        "patent",
        "A",
        "https://patents.google.com/patent/EP3856957A1/en",
        "Nykode method for selecting neoepitopes; includes ELISPOT immunogenicity examples and response rates.",
    ),
    Source(
        "NYKODE_US20220370579_THERAPEUTIC_NEOEPITOPE",
        "Nykode",
        "patent",
        "C",
        "https://patents.google.com/patent/US20220370579A1/en",
        "Nykode therapeutic anticancer neoepitope vaccine patent.",
    ),
    Source(
        "NYKODE_US20240350601_THERAPEUTIC_NEOEPITOPE",
        "Nykode",
        "patent",
        "C",
        "https://patents.google.com/patent/US20240350601A1/en",
        "Nykode therapeutic anticancer neoepitope vaccine continuation with sequence examples.",
    ),
    Source(
        "EVAXION_WO2021123232_NEOEPITOPE_CONSTRUCTS",
        "Evaxion",
        "patent",
        "C",
        "https://patents.google.com/patent/WO2021123232A1/en",
        "Evaxion nucleic acid vaccination using neo-epitope encoding constructs.",
    ),
    Source(
        "GENEOS_EP4522752_WNT_GNOSPV02",
        "Geneos",
        "patent",
        "A",
        "https://patents.google.com/patent/EP4522752A1/en",
        "Geneos WNT-related cancer vaccine patent with GNOS-PV02 ELISpot/TCR response examples.",
    ),
    Source(
        "NOUSCOM_NZ759940_MSI_SHARED_NEOANTIGENS",
        "Nouscom",
        "patent",
        "C",
        "https://patents.google.com/patent/NZ759940A/en",
        "Nouscom universal vaccine based on shared MSI tumor neoantigens.",
    ),
    Source(
        "NOUSCOM_WO2019012091_NEOANTIGEN_COMPOSITION",
        "Nouscom",
        "patent",
        "C",
        "https://patents.google.com/patent/WO2019012091A1/en",
        "Nouscom neoantigen-based vaccine composition patent.",
    ),
    Source(
        "NOUSCOM_AU2018300051_NEOANTIGEN_COMPOSITION",
        "Nouscom",
        "patent",
        "C",
        "https://patents.google.com/patent/AU2018300051B2/en",
        "Nouscom neoantigen vaccine composition patent with sequence-heavy claims.",
    ),
    Source(
        "NOUSCOM_COMPANY_PLATFORM",
        "Nouscom",
        "company_platform",
        "C",
        "https://nouscom.com/",
        "Nouscom platform page: viral vectors encoding strings of tumor neoantigens.",
    ),
    Source(
        "NOUSCOM_TECHNOLOGY",
        "Nouscom",
        "company_platform",
        "C",
        "https://nouscom.com/technology/",
        "Nouscom technology page for NOUS-209/NOUS-PEV platforms.",
    ),
    Source(
        "NOUSCOM_PRECLINICAL_NEWS",
        "Nouscom",
        "company_press",
        "B",
        "https://nouscom.com/2017/11/13/nouscoms-neoantigen-based-vaccine-synergizes-with-nktr-214-to-cure-established-tumors-in-preclinical-model/",
        "Preclinical neoantigen vaccine response context.",
    ),
    Source(
        "NIH_NEOANTIGEN_TCR_PATENT",
        "NIH",
        "patent",
        "C",
        "https://www.techtransfer.nih.gov/patent/e-067-2017-0-us-02",
        "Public patent page for isolating neoantigen-specific TCR sequences.",
    ),
    Source(
        "BAYVAX_PLATFORM",
        "BayVax",
        "company_platform",
        "D",
        "https://www.bayvaxbio.com/",
        "Company platform page for AI-powered immune programming and personalized neoantigen vaccines.",
    ),
]


def source_to_filename(source: Source) -> str:
    digest = hashlib.md5(source.url.encode()).hexdigest()[:10]
    suffix = ".pdf" if source.url.lower().split("?")[0].endswith(".pdf") else ".html"
    return f"{source.source_id}_{digest}{suffix}"


def fetch_source(source: Source, raw_dir: Path) -> tuple[bool, str, str, str]:
    raw_path = raw_dir / source_to_filename(source)
    headers = {
        "User-Agent": "Mozilla/5.0 public-benchmark-scout/0.1 (public sources only)",
        "Accept": "text/html,application/pdf,*/*",
    }
    try:
        resp = requests.get(source.url, headers=headers, timeout=30)
        resp.raise_for_status()
        raw_path.write_bytes(resp.content)
    except Exception as exc:
        return False, "", "", f"fetch_failed: {exc}"

    content_type = resp.headers.get("content-type", "")
    is_pdf = raw_path.suffix.lower() == ".pdf" or "pdf" in content_type.lower()
    try:
        if is_pdf:
            if shutil.which("pdftotext") is None:
                return True, str(raw_path), "", "pdftotext_missing"
            txt = subprocess.check_output(
                ["pdftotext", "-layout", str(raw_path), "-"],
                stderr=subprocess.DEVNULL,
                timeout=60,
            ).decode("utf-8", errors="ignore")
        else:
            soup = BeautifulSoup(resp.text, "html.parser")
            for tag in soup(["script", "style", "noscript"]):
                tag.decompose()
            txt = soup.get_text("\n")
        txt = re.sub(r"\n{3,}", "\n\n", txt)
        return True, str(raw_path), txt, ""
    except Exception as exc:
        return True, str(raw_path), "", f"text_extract_failed: {exc}"


def normalize_hla_class_i(match: re.Match[str]) -> str:
    if match.group(1):
        return f"HLA-{match.group(1).upper()}*{match.group(2)}:{match.group(3)}"
    return f"HLA-{match.group(4).upper()}*{match.group(5)[:2]}:{match.group(5)[2:]}"


def normalize_hla_class_ii(match: re.Match[str]) -> str:
    if match.group(1):
        return f"HLA-{match.group(1).upper()}*{match.group(2)}:{match.group(3)}"
    return f"HLA-{match.group(4).upper()}*{match.group(5)}:{match.group(6)}"


def line_context(text: str, token: str, width: int = 160) -> str:
    idx = text.find(token)
    if idx < 0:
        return ""
    start = max(0, idx - width)
    end = min(len(text), idx + len(token) + width)
    return re.sub(r"\s+", " ", text[start:end]).strip()


def extract_entities(source: Source, text: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    entities: list[dict[str, object]] = []
    evidence_lines: list[dict[str, object]] = []
    seen_seq: set[str] = set()
    for m in AA_RE.finditer(text):
        seq = m.group(1)
        # Avoid obvious all-one-letter or very low-complexity tokens.
        if len(set(seq)) < 4:
            continue
        if seq in seen_seq:
            continue
        seen_seq.add(seq)
        ctx = line_context(text, seq)
        entities.append(
            {
                "source_id": source.source_id,
                "organization": source.organization,
                "entity_type": "aa_sequence",
                "value": seq,
                "length": len(seq),
                "context": ctx[:500],
            }
        )
    seen_hla_i: set[str] = set()
    for m in HLA_CLASSI_RE.finditer(text):
        hla = normalize_hla_class_i(m)
        if hla in seen_hla_i:
            continue
        seen_hla_i.add(hla)
        entities.append(
            {
                "source_id": source.source_id,
                "organization": source.organization,
                "entity_type": "hla_allele",
                "mhc_class": "I",
                "value": hla,
                "length": np.nan,
                "context": line_context(text, m.group(0))[:500],
            }
        )
    seen_hla_ii: set[str] = set()
    for m in HLA_CLASSII_RE.finditer(text):
        hla = normalize_hla_class_ii(m)
        if hla in seen_hla_ii:
            continue
        seen_hla_ii.add(hla)
        entities.append(
            {
                "source_id": source.source_id,
                "organization": source.organization,
                "entity_type": "hla_classII_allele",
                "mhc_class": "II",
                "value": hla,
                "length": np.nan,
                "context": line_context(text, m.group(0))[:500],
            }
        )

    for line in text.splitlines():
        clean = re.sub(r"\s+", " ", line).strip()
        if len(clean) < 20:
            continue
        if ASSAY_RE.search(clean) or METHOD_RE.search(clean):
            evidence_lines.append(
                {
                    "source_id": source.source_id,
                    "organization": source.organization,
                    "line_type": "assay_or_method",
                    "has_assay_term": bool(ASSAY_RE.search(clean)),
                    "has_label_term": bool(LABEL_RE.search(clean)),
                    "has_method_term": bool(METHOD_RE.search(clean)),
                    "line": clean[:1000],
                }
            )
    return entities, evidence_lines[:250]


def tier_from_counts(expected: str, n_seq: int, n_hla: int, n_assay_lines: int, n_label_lines: int) -> str:
    if n_seq >= 5 and n_assay_lines >= 3 and n_label_lines >= 2:
        return "A_candidate_metric"
    if n_seq >= 1 and n_assay_lines >= 2:
        return "B_candidate_enrichment"
    if n_hla >= 1 or n_assay_lines >= 1 or expected in {"B", "C"}:
        return "C_method_or_context"
    return "D_context_only"


def process_source(source: Source, raw_dir: Path, text_dir: Path) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    ok, raw_path, text, warning = fetch_source(source, raw_dir)
    if text:
        (text_dir / f"{source.source_id}.txt").write_text(text, encoding="utf-8", errors="ignore")
    entities, evidence = extract_entities(source, text) if text else ([], [])
    n_seq = sum(1 for e in entities if e["entity_type"] == "aa_sequence")
    n_hla_i = sum(1 for e in entities if e["entity_type"] == "hla_allele")
    n_hla_ii = sum(1 for e in entities if e["entity_type"] == "hla_classII_allele")
    n_hla = n_hla_i + n_hla_ii
    n_assay = sum(1 for e in evidence if e["has_assay_term"])
    n_label = sum(1 for e in evidence if e["has_label_term"])
    inferred = tier_from_counts(source.expected_tier, n_seq, n_hla, n_assay, n_label)
    inventory_row = {
        "source_id": source.source_id,
        "organization": source.organization,
        "source_class": source.source_class,
        "expected_tier": source.expected_tier,
        "inferred_tier": inferred,
        "fetch_ok": ok,
        "url": source.url,
        "raw_path": raw_path,
        "text_chars": len(text),
        "n_aa_sequences": n_seq,
        "n_hla_alleles": n_hla,
        "n_hla_classI_alleles": n_hla_i,
        "n_hla_classII_alleles": n_hla_ii,
        "n_assay_method_lines": n_assay,
        "n_label_like_lines": n_label,
        "rationale": source.rationale,
        "warning": warning,
    }
    return inventory_row, entities, evidence


def build_inventory(workers: int = 12) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw_dir = OUT / "raw_public_sources"
    raw_dir.mkdir(parents=True, exist_ok=True)
    inventory_rows = []
    entity_rows = []
    evidence_rows = []
    text_dir = OUT / "source_text"
    text_dir.mkdir(parents=True, exist_ok=True)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(process_source, source, raw_dir, text_dir): source for source in SOURCES}
        for fut in as_completed(futures):
            try:
                inventory_row, entities, evidence = fut.result()
            except Exception as exc:
                source = futures[fut]
                inventory_row = {
                    "source_id": source.source_id,
                    "organization": source.organization,
                    "source_class": source.source_class,
                    "expected_tier": source.expected_tier,
                    "inferred_tier": "D_context_only",
                    "fetch_ok": False,
                    "url": source.url,
                    "raw_path": "",
                    "text_chars": 0,
                    "n_aa_sequences": 0,
                    "n_hla_alleles": 0,
                    "n_hla_classI_alleles": 0,
                    "n_hla_classII_alleles": 0,
                    "n_assay_method_lines": 0,
                    "n_label_like_lines": 0,
                    "rationale": source.rationale,
                    "warning": f"worker_failed: {exc}",
                }
                entities, evidence = [], []
            inventory_rows.append(inventory_row)
            entity_rows.extend(entities)
            evidence_rows.extend(evidence)
    inv = pd.DataFrame(inventory_rows)
    ent = pd.DataFrame(entity_rows)
    ev = pd.DataFrame(evidence_rows)
    return inv, ent, ev


def build_candidate_pairs(inv: pd.DataFrame, ent: pd.DataFrame, ev: pd.DataFrame) -> pd.DataFrame:
    if ent.empty:
        return pd.DataFrame()
    seqs = ent[ent["entity_type"].eq("aa_sequence")].copy()
    hlas_i = ent[ent["entity_type"].eq("hla_allele")].copy()
    hlas_ii = ent[ent["entity_type"].eq("hla_classII_allele")].copy()
    rows = []
    for _, seq in seqs.iterrows():
        source_hlas_i = hlas_i[hlas_i["source_id"].eq(seq["source_id"])]["value"].drop_duplicates().tolist()
        source_hlas_ii = hlas_ii[hlas_ii["source_id"].eq(seq["source_id"])]["value"].drop_duplicates().tolist()
        context = str(seq.get("context", ""))
        assay_context = ASSAY_RE.search(context) is not None
        label_context = LABEL_RE.search(context) is not None
        if source_hlas_i or source_hlas_ii:
            for hla in source_hlas_i[:8]:
                rows.append(
                    {
                        "source_id": seq["source_id"],
                        "organization": seq["organization"],
                        "peptide_or_sequence": seq["value"],
                        "length": seq["length"],
                        "hla_allele": hla,
                        "mhc_class": "I",
                        "candidate_label_tier": "A_like" if assay_context and label_context else "B_or_C_like",
                        "assay_context": assay_context,
                        "label_context": label_context,
                        "context": context,
                    }
                )
            for hla in source_hlas_ii[:8]:
                rows.append(
                    {
                        "source_id": seq["source_id"],
                        "organization": seq["organization"],
                        "peptide_or_sequence": seq["value"],
                        "length": seq["length"],
                        "hla_allele": hla,
                        "mhc_class": "II",
                        "candidate_label_tier": "A_like" if assay_context and label_context else "B_or_C_like",
                        "assay_context": assay_context,
                        "label_context": label_context,
                        "context": context,
                    }
                )
        else:
            rows.append(
                {
                    "source_id": seq["source_id"],
                    "organization": seq["organization"],
                    "peptide_or_sequence": seq["value"],
                    "length": seq["length"],
                    "hla_allele": "",
                    "mhc_class": "",
                    "candidate_label_tier": "A_like" if assay_context and label_context else "B_or_C_like",
                    "assay_context": assay_context,
                    "label_context": label_context,
                    "context": context,
                }
            )
    pairs = pd.DataFrame(rows)
    if pairs.empty:
        return pairs
    # Prioritize epitope-like rows for direct predictor scoring.
    pairs["direct_classI_score_ready"] = (
        pairs["mhc_class"].eq("I") & pairs["length"].between(8, 15) & pairs["hla_allele"].astype(str).ne("")
    )
    pairs["direct_classII_score_ready"] = (
        pairs["mhc_class"].eq("II") & pairs["length"].between(13, 25) & pairs["hla_allele"].astype(str).ne("")
    )
    pairs["pan_vaccine_score_ready"] = pairs["direct_classI_score_ready"] | pairs["direct_classII_score_ready"]
    pairs["slp_or_long_sequence"] = pairs["length"].astype(float).gt(15)
    pairs = pairs.sort_values(
        ["pan_vaccine_score_ready", "direct_classI_score_ready", "direct_classII_score_ready", "assay_context", "label_context", "organization", "source_id"],
        ascending=[False, False, False, False, False, True, True],
    )
    return pairs


def build_high_confidence_pairs() -> pd.DataFrame:
    """Extract stricter peptide/HLA/label rows from the same line of public text."""
    text_dir = OUT / "source_text"
    rows = []
    if not text_dir.exists():
        return pd.DataFrame()
    for path in sorted(text_dir.glob("*.txt")):
        source_id = path.stem
        text = path.read_text(encoding="utf-8", errors="ignore")
        for lineno, line in enumerate(text.splitlines(), start=1):
            clean = re.sub(r"\s+", " ", line).strip()
            if len(clean) < 20:
                continue
            seqs = []
            for m in AA_RE.finditer(clean):
                seq = m.group(1)
                if len(seq) < 8 or len(seq) > 15:
                    continue
                if seq in AA_STOPWORDS:
                    continue
                if len(set(seq)) < 4:
                    continue
                seqs.append(seq)
            if not seqs:
                continue
            hlas = [(normalize_hla_class_i(m), "I") for m in HLA_CLASSI_RE.finditer(clean)]
            hlas.extend((normalize_hla_class_ii(m), "II") for m in HLA_CLASSII_RE.finditer(clean))
            if not hlas:
                continue
            label = np.nan
            endpoint = "unknown"
            raw = ""
            if re.search(r"\bNo response\b", clean, re.IGNORECASE):
                label = 0
                endpoint = "tetramer_or_response"
                raw = "No response"
            elif re.search(r"\bRESPONSE\b", clean):
                label = 1
                endpoint = "tetramer_or_response"
                raw = "RESPONSE"
            elif re.search(r"\+\+\+|\+\+|\+", clean):
                label = 1
                endpoint = "ELISPOT_like"
                raw = "plus_sign"
            elif re.search(r"\s-\s", f" {clean} "):
                label = 0
                endpoint = "ELISPOT_like"
                raw = "minus_sign"
            for seq in sorted(set(seqs)):
                for hla, mhc_class in sorted(set(hlas)):
                    if mhc_class == "I" and not (8 <= len(seq) <= 15):
                        continue
                    if mhc_class == "II" and not (13 <= len(seq) <= 25):
                        continue
                    rows.append(
                        {
                            "source_id": source_id,
                            "line_number": lineno,
                            "peptide": seq,
                            "peptide_length": len(seq),
                            "hla_allele": hla,
                            "mhc_class": mhc_class,
                            "response_label": label,
                            "assay_endpoint_inferred": endpoint,
                            "label_raw": raw,
                            "high_confidence_metric_ready": pd.notna(label),
                            "context_line": clean[:1000],
                        }
                    )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.drop_duplicates(["source_id", "peptide", "hla_allele", "mhc_class", "response_label", "context_line"]).sort_values(
        ["high_confidence_metric_ready", "mhc_class", "source_id", "peptide"], ascending=[False, True, True, True]
    )


def build_theragen_map(inv: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "proposal_claim": "Industrial-style blind benchmark, not another academic benchmark",
            "why_theragen_should_care": "DEEPOMICS NEO already emphasizes HLA binding, T-cell activity, cleavage and SLP design; the benchmark tests those same industrial questions with public external evidence.",
            "supporting_sources": "Theragen platform/PR; TG4050; EVX-01; autogene cevumeran; patent corpus",
            "deliverable": "Freeze KG-GA/CROSS_claimsafe, score public industrial sources, compare against BigMHC and public predictors.",
        },
        {
            "proposal_claim": "Detect suspicious overfitting instead of hiding it",
            "why_theragen_should_care": "If one algorithm is too good on patent-derived data, we flag overlap/proxy/source leakage before any partnership claim.",
            "supporting_sources": "label tiering; source overlap audit; peptide/HLA duplicate audit",
            "deliverable": "A/B/C/D source tiers, overlap report, hard negative and low-confidence label separation.",
        },
        {
            "proposal_claim": "Theragen-specific differentiator: East Asian rare HLA and SLP positioning",
            "why_theragen_should_care": "Their public materials explicitly claim rare HLA and SLP mutation-position design; BioDarwin-RL can turn these into explicit method primitives and validation endpoints.",
            "supporting_sources": "Theragen DEEPOMICS NEO public page and patent announcement",
            "deliverable": "Rare-HLA and SLP-cleavage benchmark lanes; not just 9-mer class-I ranking.",
        },
        {
            "proposal_claim": "Wetlab-ready closing experiment",
            "why_theragen_should_care": "Public industrial benchmark can choose cases where KG-GA, DEEPOMICS-like logic and BigMHC disagree, then test disagreement in 96-well assays.",
            "supporting_sources": "TG4050 peptide/HLA/ELISpot examples; EVX-01 immunogenicity; autogene cevumeran ELISpot",
            "deliverable": "Disagreement plate: KG-GA-high/BigMHC-low, BigMHC-high/KG-GA-low, both-high, both-low controls.",
        },
    ]
    return pd.DataFrame(rows)


def write_outputs(inv: pd.DataFrame, ent: pd.DataFrame, ev: pd.DataFrame, pairs: pd.DataFrame, proposal: pd.DataFrame) -> None:
    inv.to_csv(OUT / "public_industrial_source_inventory.tsv", sep="\t", index=False)
    ent.to_csv(OUT / "public_industrial_extracted_entities.tsv", sep="\t", index=False)
    ev.to_csv(OUT / "public_industrial_evidence_lines.tsv", sep="\t", index=False)
    pairs.to_csv(OUT / "public_industrial_candidate_pairs.tsv", sep="\t", index=False)
    proposal.to_csv(OUT / "theragen_proposal_map.tsv", sep="\t", index=False)


def write_locked_testset_inputs(high_conf: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if high_conf.empty:
        locked = pd.DataFrame()
        locked_ii = pd.DataFrame()
        locked_pan = pd.DataFrame()
        bigmhc = pd.DataFrame()
    else:
        metric_ready = high_conf[high_conf["high_confidence_metric_ready"]].copy()
        locked = metric_ready[
            metric_ready["mhc_class"].eq("I")
            & metric_ready["peptide_length"].between(8, 15)
            & metric_ready["hla_allele"].astype(str).ne("")
        ].copy()
        locked = locked.reset_index(drop=True)
        locked.insert(0, "industrial_candidate_id", [f"INDPUB_{i+1:04d}" for i in range(len(locked))])
        locked["label"] = pd.to_numeric(locked["response_label"], errors="coerce").astype(int)
        locked["benchmark_use"] = "locked_public_industrial_classI_metric_after_manual_QA"
        locked["manual_QA_required"] = True
        locked["do_not_train"] = True
        locked_ii = metric_ready[
            metric_ready["mhc_class"].eq("II")
            & metric_ready["peptide_length"].between(13, 25)
            & metric_ready["hla_allele"].astype(str).ne("")
        ].copy()
        locked_ii = locked_ii.reset_index(drop=True)
        locked_ii.insert(0, "industrial_candidate_id", [f"INDPUBII_{i+1:04d}" for i in range(len(locked_ii))])
        if not locked_ii.empty:
            locked_ii["label"] = pd.to_numeric(locked_ii["response_label"], errors="coerce").astype(int)
        else:
            locked_ii["label"] = pd.Series(dtype=int)
        locked_ii["benchmark_use"] = "locked_public_industrial_classII_metric_after_manual_QA"
        locked_ii["manual_QA_required"] = True
        locked_ii["do_not_train"] = True
        locked_pan = pd.concat([locked, locked_ii], ignore_index=True, sort=False)
        if not locked_pan.empty:
            locked_pan["pan_vaccine_axis"] = np.where(locked_pan["mhc_class"].eq("II"), "CD4_helper", "CD8_cytotoxic")
        bigmhc = locked[
            ["industrial_candidate_id", "hla_allele", "peptide", "label", "source_id"]
        ].rename(
            columns={
                "industrial_candidate_id": "candidate_id",
                "hla_allele": "mhc",
                "peptide": "pep",
                "source_id": "source_name",
            }
        )
    locked.to_csv(OUT / "public_industrial_locked_testset_v0.tsv", sep="\t", index=False)
    locked_ii.to_csv(OUT / "public_industrial_locked_classII_testset_v0.tsv", sep="\t", index=False)
    locked_pan.to_csv(OUT / "public_industrial_locked_pan_vaccine_testset_v0.tsv", sep="\t", index=False)
    bigmhc.to_csv(OUT / "public_industrial_locked_testset_bigmhc_input_v0.csv", index=False)
    return locked, bigmhc, locked_ii, locked_pan


def plot_outputs(inv: pd.DataFrame, pairs: pd.DataFrame, ev: pd.DataFrame) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 4.8))
    counts = inv["inferred_tier"].value_counts().sort_index()
    ax.bar(counts.index, counts.values, color="#2f6f73")
    ax.set_ylabel("Public sources")
    ax.set_title("Industrial public source tiers")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_source_tiers.png", dpi=190)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5.2))
    org = inv.groupby("organization")[["n_aa_sequences", "n_hla_alleles", "n_assay_method_lines"]].sum().sort_values("n_assay_method_lines")
    org.plot(kind="barh", ax=ax, color=["#69a7a2", "#c59b3b", "#8f9cc2"])
    ax.set_xlabel("Extracted public evidence count")
    ax.set_title("Evidence density by organization/source family")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_evidence_density.png", dpi=190)
    plt.close(fig)

    if not pairs.empty:
        fig, ax = plt.subplots(figsize=(8.8, 4.8))
        vals = pd.Series(
            {
                "Class-I score ready": int(pairs["direct_classI_score_ready"].sum()),
                "Class-II score ready": int(pairs["direct_classII_score_ready"].sum()) if "direct_classII_score_ready" in pairs else 0,
                "Needs parsing/SLP": int((~pairs["pan_vaccine_score_ready"]).sum()) if "pan_vaccine_score_ready" in pairs else int((~pairs["direct_classI_score_ready"]).sum()),
            }
        )
        vals = vals[vals.gt(0)]
        ax.pie(vals.values, labels=vals.index, autopct="%1.0f%%", colors=["#2f6f73", "#d16b5f"])
        ax.set_title("Direct scoring readiness of extracted sequence/HLA candidates")
        fig.tight_layout()
        fig.savefig(fig_dir / "fig3_scoring_readiness.png", dpi=190)
        plt.close(fig)

    graph = nx.DiGraph()
    graph.add_node("Public Industrial OOD Benchmark", node_type="benchmark")
    for _, row in inv.iterrows():
        graph.add_node(row["organization"], node_type="organization")
        graph.add_node(row["source_id"], node_type=row["source_class"], tier=row["inferred_tier"])
        graph.add_edge(row["organization"], row["source_id"], relation="publishes")
        graph.add_edge(row["source_id"], "Public Industrial OOD Benchmark", relation="candidate_source")
    fig, ax = plt.subplots(figsize=(12, 8))
    pos = nx.spring_layout(graph, seed=21, k=0.9)
    color_map = {"benchmark": "#c59b3b", "organization": "#69a7a2"}
    colors = [color_map.get(graph.nodes[n].get("node_type", ""), "#d7dde2") for n in graph.nodes]
    sizes = [900 if n == "Public Industrial OOD Benchmark" else 280 if graph.nodes[n].get("node_type") == "organization" else 110 for n in graph.nodes]
    nx.draw_networkx_nodes(graph, pos, node_color=colors, node_size=sizes, edgecolors="#27323a", linewidths=0.5, ax=ax)
    nx.draw_networkx_edges(graph, pos, arrows=False, alpha=0.28, edge_color="#9aa7af", ax=ax)
    labels = {n: n for n in graph.nodes if graph.nodes[n].get("node_type") in {"benchmark", "organization"}}
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=7, ax=ax)
    ax.set_title("Public industrial neoantigen benchmark source graph")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig4_source_graph.png", dpi=220)
    plt.close(fig)


def table_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "<p>No rows.</p>"
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_report(
    inv: pd.DataFrame,
    ent: pd.DataFrame,
    ev: pd.DataFrame,
    pairs: pd.DataFrame,
    proposal: pd.DataFrame,
    high_conf: pd.DataFrame,
) -> None:
    tier_counts = inv["inferred_tier"].value_counts().to_dict()
    ready_n = int(pairs["direct_classI_score_ready"].sum()) if not pairs.empty else 0
    ready_ii_n = int(pairs["direct_classII_score_ready"].sum()) if not pairs.empty and "direct_classII_score_ready" in pairs else 0
    high_conf_n = int(high_conf["high_confidence_metric_ready"].sum()) if not high_conf.empty else 0
    lines = [
        "# Public industrial neoantigen benchmark scout",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Bottom line",
        "",
        "- This is a public-only industrial benchmark scout for a Theragen-facing proposal.",
        "- Sources are patents, company pages, posters and clinical publications, not confidential materials.",
        "- Use these data as a locked out-of-distribution testset after freezing KG-GA/CROSS_claimsafe; do not train on them.",
        f"- Source inventory: {len(inv)} public sources.",
        f"- Extracted entities: {len(ent)} sequence/HLA entities.",
        f"- Candidate peptide/HLA rows: {len(pairs)}, direct class-I score-ready rows: {ready_n}, direct class-II score-ready rows: {ready_ii_n}.",
        f"- High-confidence same-line peptide/HLA/label rows: {high_conf_n}.",
        f"- Tier counts: {tier_counts}.",
        "",
        "## Source inventory",
        "",
        inv.sort_values(["inferred_tier", "organization"]).to_markdown(index=False),
        "",
        "## Candidate peptide/HLA extraction preview",
        "",
        pairs.head(80).to_markdown(index=False) if not pairs.empty else "No candidate pairs extracted.",
        "",
        "## High-confidence same-line candidate rows",
        "",
        high_conf.head(120).to_markdown(index=False) if not high_conf.empty else "No high-confidence candidate rows extracted.",
        "",
        "## Theragen proposal map",
        "",
        proposal.to_markdown(index=False),
        "",
        "## Guardrails",
        "",
        "- Public sources only.",
        "- No confidential, leaked or terms-violating company material.",
        "- No training on these sources before benchmark evaluation.",
        "- Treat patent/company selected peptides without true negatives as enrichment tests, not supervised metrics.",
        "- Run overlap audit against CEDAR/IEDB/TESLA/NEPdb before scoring.",
        "- If a model is suspiciously perfect on one company source, flag possible sequence/source overlap rather than claiming victory.",
        "- Wetlab use of patented sequences requires FTO/legal review; computational benchmarking is the intended first use.",
        "",
    ]
    (OUT / "PUBLIC_INDUSTRIAL_BENCHMARK_SCOUT_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def write_html(
    inv: pd.DataFrame,
    ent: pd.DataFrame,
    ev: pd.DataFrame,
    pairs: pd.DataFrame,
    proposal: pd.DataFrame,
    high_conf: pd.DataFrame,
) -> dict[str, object]:
    asset_dir = HUB / "assets/public_industrial_neoantigen_benchmark_scout"
    asset_dir.mkdir(parents=True, exist_ok=True)
    for fig in (OUT / "figures").glob("*.png"):
        shutil.copy2(fig, asset_dir / fig.name)
    ready_n = int(pairs["direct_classI_score_ready"].sum()) if not pairs.empty else 0
    ready_ii_n = int(pairs["direct_classII_score_ready"].sum()) if not pairs.empty and "direct_classII_score_ready" in pairs else 0
    high_conf_n = int(high_conf["high_confidence_metric_ready"].sum()) if not high_conf.empty else 0
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Public Industrial Neoantigen Benchmark Scout</title>
  <style>
    :root {{ --bg:#101418; --panel:#151c22; --ink:#e9edf0; --muted:#9aa7af; --gold:#c59b3b; --line:#29343c; --teal:#69a7a2; }}
    body {{ margin:0; background:var(--bg); color:var(--ink); font:15px/1.55 system-ui, -apple-system, Segoe UI, sans-serif; }}
    header {{ padding:42px 5vw 30px; background:#0c1115; border-bottom:1px solid var(--line); }}
    .kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.12em; font-weight:700; font-size:12px; }}
    h1 {{ margin:.3rem 0 .65rem; font-size:clamp(30px,4vw,54px); line-height:1.05; }}
    .lead {{ color:#d0d8dd; max-width:1040px; font-size:18px; }}
    .stats {{ display:grid; grid-template-columns:repeat(5,minmax(130px,1fr)); gap:10px; margin-top:22px; }}
    .stat {{ border:1px solid var(--line); background:var(--panel); padding:13px 14px; border-radius:8px; }}
    .stat b {{ display:block; font-size:23px; }}
    .stat span {{ color:var(--muted); font-size:12px; }}
    main {{ padding:30px 5vw 60px; display:grid; grid-template-columns:270px minmax(0,1fr); gap:28px; }}
    nav {{ position:sticky; top:0; align-self:start; max-height:100vh; overflow:auto; padding:14px; border:1px solid var(--line); border-radius:8px; background:#111820; }}
    nav a {{ display:block; color:#d6dee2; text-decoration:none; padding:7px 0; border-bottom:1px solid #202b32; }}
    h2 {{ border-bottom:1px solid var(--line); padding-bottom:8px; }}
    .num {{ color:var(--gold); margin-right:8px; }}
    .note {{ border-left:3px solid var(--gold); background:#171f26; padding:10px 14px; color:#dce3e7; }}
    table.data {{ width:100%; border-collapse:collapse; font-size:13px; margin:12px 0 28px; }}
    table.data th, table.data td {{ border-bottom:1px solid var(--line); padding:7px 8px; text-align:left; vertical-align:top; }}
    table.data th {{ color:#f4d891; background:#141b21; }}
    .grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:18px; }}
    img {{ max-width:100%; border:1px solid var(--line); border-radius:8px; background:#fff; }}
    code {{ color:#f4d891; }}
    @media(max-width:900px) {{ main {{ grid-template-columns:1fr; }} nav {{ position:static; }} .stats {{ grid-template-columns:repeat(2,1fr); }} .grid {{ grid-template-columns:1fr; }} }}
  </style>
</head>
<body>
<header>
  <div class="kicker">Public-only industrial OOD benchmark · Theragen proposal scaffold</div>
  <h1>Patent and company-source neoantigen benchmark scout</h1>
  <p class="lead">This page inventories public industrial sources that can become a locked out-of-distribution benchmark for KG-GA, CROSS_claimsafe and public comparators. It is not training data and excludes confidential material.</p>
  <div class="stats">
    <div class="stat"><b>{len(inv)}</b><span>public sources</span></div>
    <div class="stat"><b>{len(ent)}</b><span>sequence/HLA entities</span></div>
    <div class="stat"><b>{len(pairs)}</b><span>candidate rows</span></div>
    <div class="stat"><b>{ready_n}</b><span>class-I score-ready rows</span></div>
    <div class="stat"><b>{ready_ii_n}</b><span>class-II score-ready rows</span></div>
    <div class="stat"><b>{high_conf_n}</b><span>same-line labeled rows</span></div>
    <div class="stat"><b>{inv['organization'].nunique()}</b><span>organization families</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#summary">01 Summary</a>
  <a href="#figures">02 Figures</a>
  <a href="#sources">03 Sources</a>
  <a href="#candidates">04 Candidate Rows</a>
  <a href="#highconf">05 High-Confidence Rows</a>
  <a href="#proposal">06 Theragen Map</a>
  <a href="#evidence">07 Evidence Lines</a>
  <a href="#guardrails">08 Guardrails</a>
</nav>
<article>
<section id="summary">
  <h2><span class="num">01</span>Summary</h2>
  <p class="note">Use this as a scout layer. The next locked step is source QA, overlap removal, manual curation of Tier A labels, then frozen scoring by KG-GA, BigMHC and CROSS_claimsafe.</p>
</section>
<section id="figures">
  <h2><span class="num">02</span>Figures</h2>
  <div class="grid">
    <img src="assets/public_industrial_neoantigen_benchmark_scout/fig1_source_tiers.png" alt="source tiers">
    <img src="assets/public_industrial_neoantigen_benchmark_scout/fig2_evidence_density.png" alt="evidence density">
  </div>
  <div class="grid">
    <img src="assets/public_industrial_neoantigen_benchmark_scout/fig3_scoring_readiness.png" alt="scoring readiness">
    <img src="assets/public_industrial_neoantigen_benchmark_scout/fig4_source_graph.png" alt="source graph">
  </div>
</section>
<section id="sources">
  <h2><span class="num">03</span>Source Inventory</h2>
  {table_html(inv.sort_values(['inferred_tier','organization']), 60)}
</section>
<section id="candidates">
  <h2><span class="num">04</span>Extracted Candidate Rows</h2>
  {table_html(pairs, 80)}
</section>
<section id="highconf">
  <h2><span class="num">05</span>High-Confidence Same-Line Rows</h2>
  {table_html(high_conf, 100)}
</section>
<section id="proposal">
  <h2><span class="num">06</span>Theragen Proposal Map</h2>
  {table_html(proposal, 20)}
</section>
<section id="evidence">
  <h2><span class="num">07</span>Evidence Lines</h2>
  {table_html(ev.sort_values(['has_label_term','has_assay_term'], ascending=[False, False]), 80)}
</section>
<section id="guardrails">
  <h2><span class="num">08</span>Guardrails</h2>
  <p>Public sources only. No confidential material. No training on this scout set before benchmark evaluation. Patent/company selected peptides without true negatives are enrichment tests, not supervised metrics. Any suspiciously perfect model result triggers overlap/proxy audit.</p>
</section>
</article>
</main>
</body>
</html>
"""
    html_path = HUB / "public_industrial_neoantigen_benchmark_scout.html"
    html_path.write_text(html_text, encoding="utf-8")
    live_ok = False
    warnings: list[str] = []
    try:
        live_asset_dir = LIVE_HUB / "assets/public_industrial_neoantigen_benchmark_scout"
        live_asset_dir.mkdir(parents=True, exist_ok=True)
        for fig in asset_dir.glob("*.png"):
            shutil.copy2(fig, live_asset_dir / fig.name)
        shutil.copy2(html_path, LIVE_HUB / html_path.name)
        live_ok = True
    except Exception as exc:
        warnings.append(str(exc))
    return {
        "html_path": str(html_path),
        "live_html_path": str(LIVE_HUB / html_path.name),
        "live_deploy_ok": live_ok,
        "live_deploy_warnings": warnings,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    inv, ent, ev = build_inventory()
    pairs = build_candidate_pairs(inv, ent, ev)
    high_conf = build_high_confidence_pairs()
    locked, bigmhc_input, locked_ii, locked_pan = write_locked_testset_inputs(high_conf)
    proposal = build_theragen_map(inv)
    write_outputs(inv, ent, ev, pairs, proposal)
    high_conf.to_csv(OUT / "public_industrial_high_confidence_candidate_pairs.tsv", sep="\t", index=False)
    plot_outputs(inv, pairs, ev)
    write_report(inv, ent, ev, pairs, proposal, high_conf)
    html_info = write_html(inv, ent, ev, pairs, proposal, high_conf)
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "out_dir": str(OUT),
        "n_public_sources": int(len(inv)),
        "n_fetch_ok": int(inv["fetch_ok"].sum()),
        "n_organizations": int(inv["organization"].nunique()),
        "n_extracted_entities": int(len(ent)),
        "n_candidate_pairs": int(len(pairs)),
        "n_direct_classI_score_ready": int(pairs["direct_classI_score_ready"].sum()) if not pairs.empty else 0,
        "n_direct_classII_score_ready": int(pairs["direct_classII_score_ready"].sum()) if not pairs.empty and "direct_classII_score_ready" in pairs else 0,
        "n_high_confidence_candidate_pairs": int(len(high_conf)),
        "n_high_confidence_metric_ready": int(high_conf["high_confidence_metric_ready"].sum()) if not high_conf.empty else 0,
        "n_locked_classI_testset_v0": int(len(locked)),
        "n_locked_classII_testset_v0": int(len(locked_ii)),
        "n_locked_pan_vaccine_testset_v0": int(len(locked_pan)),
        "n_locked_testset_v0": int(len(locked)),
        "tier_counts": inv["inferred_tier"].value_counts().to_dict(),
        "claim_boundary": "public-only scout; use as locked OOD benchmark after algorithm freeze, not as training data",
        **html_info,
        "output_files": [
            str(OUT / "public_industrial_source_inventory.tsv"),
            str(OUT / "public_industrial_extracted_entities.tsv"),
            str(OUT / "public_industrial_evidence_lines.tsv"),
            str(OUT / "public_industrial_candidate_pairs.tsv"),
            str(OUT / "public_industrial_high_confidence_candidate_pairs.tsv"),
            str(OUT / "public_industrial_locked_testset_v0.tsv"),
            str(OUT / "public_industrial_locked_classII_testset_v0.tsv"),
            str(OUT / "public_industrial_locked_pan_vaccine_testset_v0.tsv"),
            str(OUT / "public_industrial_locked_testset_bigmhc_input_v0.csv"),
            str(OUT / "theragen_proposal_map.tsv"),
            str(OUT / "PUBLIC_INDUSTRIAL_BENCHMARK_SCOUT_REPORT.md"),
        ],
    }
    (OUT / "public_industrial_benchmark_scout_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
