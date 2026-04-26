"""Inject beginner-friendly education toggles into every THCA dashboard page.

Design:
  - Idempotent: re-running does not duplicate toggles. Marker comment
    ``<!-- edu:<section_id> -->`` anchors the per-section insertion; marker
    ``<!-- edu-css -->`` anchors the stylesheet link; marker
    ``<!-- edu-glossary-nav -->`` anchors the topnav "용어" link injection.
  - Data-driven: content is loaded from
    ``reports/html/assets/data/education_content.json``. If that file is
    missing, a built-in default is written out first so operators can edit it
    and re-run.
  - Non-destructive: we do NOT modify templates, the build pipeline, or any
    Python scripts. Edits are localized inline replacements of existing tags.
  - Markdown is pre-rendered to HTML at injection time (small, deterministic
    subset: headings, bold, inline code, lists, paragraphs). This avoids a
    client-side render race.

Run:
    python3 scripts/inject_education_toggles.py

The module is safe to import; calling ``main()`` performs the injection.
"""
from __future__ import annotations

import html
import json
import re
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
HTML_ROOT = ROOT / "reports" / "html"
PAGES_DIR = HTML_ROOT / "pages"
INDEX_PATH = HTML_ROOT / "index.html"
CONTENT_PATH = HTML_ROOT / "assets" / "data" / "education_content.json"
CSS_PATH = HTML_ROOT / "assets" / "css" / "edu-toggle.css"
GLOSSARY_PATH = PAGES_DIR / "99_glossary.html"


# ---------------------------------------------------------------------------
# Tiny markdown -> HTML renderer (sufficient for our education bodies).
# ---------------------------------------------------------------------------

_INLINE_CODE = re.compile(r"`([^`]+)`")
_BOLD = re.compile(r"\*\*([^*]+)\*\*")
_ITALIC = re.compile(r"(?<!\*)\*([^*]+)\*(?!\*)")


def _render_inline(line: str) -> str:
    # Escape first, then re-inject formatting markers.
    out = html.escape(line, quote=False)
    out = _INLINE_CODE.sub(lambda m: f"<code>{m.group(1)}</code>", out)
    out = _BOLD.sub(lambda m: f"<strong>{m.group(1)}</strong>", out)
    out = _ITALIC.sub(lambda m: f"<em>{m.group(1)}</em>", out)
    return out


def markdown_to_html(text: str) -> str:
    """Minimal markdown renderer: #, ##, ###, -, 1. / paragraphs / inline."""
    lines = text.strip().splitlines()
    out: list[str] = []
    i = 0
    n = len(lines)
    while i < n:
        raw = lines[i].rstrip()
        stripped = raw.strip()
        if not stripped:
            i += 1
            continue
        # Headings
        m = re.match(r"^(#{1,4})\s+(.*)$", stripped)
        if m:
            level = len(m.group(1))
            content = _render_inline(m.group(2).strip())
            # Cap visual heading at h4 (inside a details block)
            tag = f"h{min(level + 2, 5)}"  # h1 -> h3, h2 -> h4
            out.append(f"<{tag}>{content}</{tag}>")
            i += 1
            continue
        # Unordered list
        if stripped.startswith("- ") or stripped.startswith("* "):
            items: list[str] = []
            while i < n:
                s = lines[i].strip()
                if s.startswith("- ") or s.startswith("* "):
                    items.append(_render_inline(s[2:].strip()))
                    i += 1
                elif not s:
                    i += 1
                    break
                else:
                    break
            out.append("<ul>" + "".join(f"<li>{it}</li>" for it in items) + "</ul>")
            continue
        # Ordered list
        if re.match(r"^\d+\.\s+", stripped):
            items = []
            while i < n:
                s = lines[i].strip()
                m2 = re.match(r"^\d+\.\s+(.*)$", s)
                if m2:
                    items.append(_render_inline(m2.group(1).strip()))
                    i += 1
                elif not s:
                    i += 1
                    break
                else:
                    break
            out.append("<ol>" + "".join(f"<li>{it}</li>" for it in items) + "</ol>")
            continue
        # Paragraph: collect until blank or structural line
        buf: list[str] = [stripped]
        i += 1
        while i < n:
            s = lines[i].strip()
            if not s:
                break
            if s.startswith(("#", "- ", "* ")) or re.match(r"^\d+\.\s+", s):
                break
            buf.append(s)
            i += 1
        out.append("<p>" + _render_inline(" ".join(buf)) + "</p>")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Education content: (page_basename, section_id) -> {title, summary, body_md}
# Page basename includes the ``.html`` suffix. ``index.html`` uses full name.
# ---------------------------------------------------------------------------


def _build_default_content() -> dict:
    """Seed content dictionary. Covers every major section across pages."""

    generic_template = {
        "title": "이 섹션 읽는 법",
        "summary": "섹션을 처음 보는 분을 위한 안내",
        "body_md": (
            "## 무엇을 보여주나\n"
            "- 이 섹션은 페이지 상단의 주제를 상세히 풀어놓은 블록입니다.\n"
            "- 표·그림·지표가 섞여 있을 수 있으며, 각 요소는 바로 위 제목과 같은 질문에 답합니다.\n\n"
            "## 초보자용 체크리스트\n"
            "- 먼저 페이지의 hero 영역 문장을 읽고 \"무엇을 비교하는가\"를 확인하세요.\n"
            "- 각 카드/지표의 `subtle mono` 캡션을 읽으면 단위를 알 수 있습니다.\n"
            "- 본문에 등장하는 약어(AUC, BRS, TDS 등)가 낯설면 상단 내비게이션의 `용어` 페이지를 참고하세요.\n\n"
            "## 주의할 점\n"
            "- 표/그림의 숫자는 한 번의 빌드 결과이며, 다른 코호트에서는 달라질 수 있습니다.\n"
            "- 외부 데이터셋과 내부 교차검증 숫자를 혼동하지 마세요."
        ),
    }

    content: dict[str, dict] = {}

    def add(key: str, title: str, summary: str, body: str) -> None:
        content[key] = {"title": title, "summary": summary, "body_md": body}

    # ---- index.html ----
    add(
        "index.html:findings",
        "Key Findings 한눈에 읽기",
        "상단 요약 카드가 말하는 것",
        "## 이 섹션이 답하는 질문\n"
        "- 이 프로젝트에서 지금 가장 핵심이 되는 숫자 세 가지가 무엇인가?\n"
        "- 내부 데이터(TCGA)와 외부 데이터(GEO)의 성능 차이를 어떻게 보나?\n\n"
        "## 어떻게 읽는가\n"
        "- `Best internal AUC 0.935`는 **TCGA-THCA** 데이터에서 5-fold 교차검증 중 최고점입니다.\n"
        "- `현실적 외부 AUC 0.980`은 GSE27155 데이터로 확인한 `sanity-check`이며, 조직학 기반 proxy라 실제 molecular subtype 예측보다 과대 추정될 수 있습니다.\n"
        "- `라벨 수정 632`는 `GSE213647`에서 라벨 오류를 잡아 다시 정렬한 샘플 수입니다.\n\n"
        "## 주의할 점\n"
        "- internal AUC 단독으로는 성능을 판단하기 어렵습니다. 외부 AUC, calibration, threshold metric까지 함께 봐야 합니다.\n"
        "- histology-proxy AUC는 **실제 BRAF/RAS 분자 라벨이 아닌** 조직학 라벨을 기준으로 한 값이라는 점을 기억하세요.",
    )
    add(
        "index.html:explorer",
        "Interactive Explorer 사용법",
        "14개 페이지를 타일로 탐색",
        "## 무엇을 보여주나\n"
        "- 대시보드 전체 페이지를 카드 그리드로 보여주는 허브입니다.\n"
        "- 각 타일의 제목을 클릭하면 해당 페이지로 이동합니다.\n\n"
        "## 어떻게 읽는가\n"
        "- 카드의 `subtle mono` 번호는 페이지 순서이며, 분석 파이프라인 순서와 일치합니다.\n"
        "- 분석 흐름상 `01 → 04`(데이터와 패널) → `05/06`(구조와 점수) → `07/08`(모델) → `09/10/11`(해석과 비교) → `12`(사업성)의 순서로 읽을 때 가장 이해가 빠릅니다.\n\n"
        "## 단축키 팁\n"
        "- `⌘K` 또는 `Ctrl+K`로 검색 팔레트를 열어 섹션 제목만 입력해도 바로 점프할 수 있습니다.",
    )
    add(
        "index.html:latest",
        "Latest Build 의미",
        "이번 빌드가 사용한 데이터 스냅샷",
        "## 이 숫자들의 뜻\n"
        "- `Build <타임스탬프>`는 **이 페이지가 언제 렌더링됐는지**입니다.\n"
        "- `7 datasets`, `1509 samples`, `51711 genes`는 실제로 사용된 데이터 크기입니다.\n\n"
        "## 아래 그림 읽는 법\n"
        "- `panel perf curve`는 16 / 32 / 67 / 112 유전자 패널에서 AUC가 어떻게 변하는지 비교한 곡선입니다.\n"
        "- 유전자 수가 늘어날수록 AUC는 보통 오르지만, 일정 지점부터는 상승 폭이 작아집니다(Occam의 면도날).\n\n"
        "## 주의할 점\n"
        "- 빌드 해시가 `no-git`이면 로컬 실험 결과이며, 동료와 공유할 때는 git revision을 재확인해야 합니다.",
    )
    add(
        "index.html:highlights",
        "Visual Highlights 읽는 법",
        "시각적 하이라이트 미리보기",
        "## 무엇을 보여주나\n"
        "- 전체 분석에서 시각적으로 가장 설득력이 있는 Figure 몇 장을 preview 형태로 보여줍니다.\n\n"
        "## 어떻게 읽는가\n"
        "- 각 Figure는 PNG preview → 확대해서 보기(Plotly 인터랙티브) → 새 탭(원본 HTML) 순으로 접근 가능합니다.\n"
        "- `Cite` 버튼을 누르면 인용 가능한 경로/타이틀이 복사됩니다.",
    )

    # ---- 01 overview ----
    add(
        "01_overview.html:hero",
        "연구 개요와 핵심 구조",
        "THCA와 BRAF/RAS 축을 30초에 이해하기",
        "## 왜 이 연구를 하는가\n"
        "- THCA(Thyroid Carcinoma, 갑상선암)는 상당수가 PTC(Papillary Thyroid Carcinoma)로, 임상적으로 **BRAF-like**와 **RAS-like** 두 분자 subtype으로 크게 갈립니다.\n"
        "- BRAF-like와 RAS-like는 예후·치료 반응·dedifferentiation 위험이 다르기 때문에 어느 쪽인지를 빠르고 저렴하게 판단하는 것이 임상 가치가 있습니다.\n\n"
        "## molecular subtype이 뭔가\n"
        "- 형태(histology)만으로는 구분이 어려운 환자들을 **유전자 발현 패턴**으로 분류한 라벨입니다.\n"
        "- TCGA-THCA 논문(2014, Cell)에서 제안된 `BRS(BRAF-RAS Score)` 축이 사실상 산업 표준입니다.\n\n"
        "## 이 대시보드가 답하는 질문\n"
        "- 공개 데이터만으로 BRAF-like vs RAS-like 분류가 가능한가?\n"
        "- 16~112개 유전자 정도의 **compact panel**이면 충분한가, 아니면 전체 전사체가 필요한가?\n"
        "- 내부(TCGA) 성능이 외부(GEO) 코호트에서 재현되는가?",
    )
    add(
        "01_overview.html:interpretation",
        "이 페이지에서 확인할 것",
        "어떤 순서로 대시보드를 둘러볼지",
        "## 권장 읽기 순서\n"
        "1. 먼저 이 개요 페이지에서 **샘플 수/데이터 출처/주된 축**을 확인합니다.\n"
        "2. `02 데이터셋`에서 각 코호트의 플랫폼(RNA-seq vs microarray)과 샘플 크기를 확인합니다.\n"
        "3. `03 샘플`에서 BRAF_like 1,075 / RAS_like 159 / dedifferentiated 82가 어떻게 나뉘었는지 봅니다.\n"
        "4. `04 패널` → `07 ML` → `08 패널 비교`로 모델 성능을 확인합니다.\n\n"
        "## 용어 미리보기\n"
        "- **cPTC**: classical PTC (가장 흔한 유두암 아형)\n"
        "- **FVPTC**: follicular variant PTC (여포성 변종, RAS-like에 가까움)\n"
        "- **FTC/ATC/PDTC/MTC**: 여포암/미분화/저분화/수질성 (각각 드물지만 임상적으로 중요)\n\n"
        "## 주의할 점\n"
        "- mutation anchor(예: BRAF V600E 검출)는 보조 컬럼이며 최종 라벨은 histology + molecular annotation 일치 시에만 확정합니다.",
    )

    # ---- 02 datasets ----
    add(
        "02_datasets.html:dataset-summary",
        "데이터셋 요약: TCGA와 GEO 이해하기",
        "공개 데이터 출처와 기술 플랫폼",
        "## 사용한 공개 데이터\n"
        "- **TCGA (The Cancer Genome Atlas)**: 미국 NCI/NHGRI 주도 대규모 암 유전체 프로젝트. THCA 코호트는 약 500명 PTC 환자의 RNA-seq / methylation / mutation을 제공합니다.\n"
        "- **GEO (Gene Expression Omnibus)**: NCBI가 운영하는 공개 유전자 발현 저장소. 각 연구실이 직접 업로드한 microarray·RNA-seq 데이터가 모여 있습니다.\n\n"
        "## 기술 플랫폼이 뭔가\n"
        "- **bulk RNA-seq**: 조직 전체를 갈아 mRNA를 대량 시퀀싱. 가장 정밀하지만 비쌉니다.\n"
        "- **microarray**: 미리 설계된 probe에 hybridization. 저렴하고 오래된 방식이며 동적 범위가 좁습니다.\n"
        "- **methylation (450K / EPIC array)**: DNA의 CpG 부위 메틸화 정도를 측정. 이 대시보드에서는 QC 전용입니다.\n\n"
        "## 왜 혼합하지 않나\n"
        "- RNA-seq과 microarray는 값의 스케일과 noise 구조가 달라 섞어 학습하면 **batch effect**가 성능을 흐립니다.",
    )
    add(
        "02_datasets.html:dataset-table",
        "데이터셋 테이블 읽는 법",
        "표 컬럼이 의미하는 것",
        "## 각 컬럼 의미\n"
        "- `sample count`: 해당 데이터셋에서 실제 사용한 샘플 수 (QC 통과 기준).\n"
        "- `platform`: Illumina HiSeq / Affymetrix U133 / Agilent 등 측정 장비.\n"
        "- `normal/tumor`: 정상 조직 샘플과 종양 샘플 비율. 정상 비율이 너무 낮으면 배경이 되는 발현 범위를 가늠하기 어렵습니다.\n"
        "- `A/B grade`: 우리가 매긴 품질 등급(메타데이터 완전성, sample size, reproducibility).\n\n"
        "## 활용 팁\n"
        "- 검색 박스에 `TCGA`, `GSE`, `PTC` 등 키워드를 넣어 즉시 필터링할 수 있습니다.\n"
        "- 정렬 가능한 컬럼(▲/▼)은 sample size 비교에 유용합니다.",
    )
    add(
        "02_datasets.html:dataset-viz",
        "전역 뷰 그래프",
        "데이터셋을 한 화면에 시각화",
        "## 이 그래프가 보여주는 것\n"
        "- 각 데이터셋의 **상대적 크기**와 **플랫폼 조합**을 한눈에 비교할 수 있게 합니다.\n\n"
        "## 읽을 때 주의\n"
        "- bar 길이는 샘플 수이지 품질이 아닙니다. 품질은 `A/B grade` 표기를 함께 봐야 합니다.\n"
        "- 정상 샘플이 거의 없는 데이터셋은 종양-정상 비교에 사용할 수 없습니다.",
    )
    add(
        "02_datasets.html:dataset-supplement",
        "추가 데이터셋 Figure Wall",
        "보조 그림 모음",
        "## 언제 참고하나\n"
        "- 메인 분석에 직접 들어가지는 않지만 각 데이터셋의 특이한 분포, 결측 패턴을 확인하고 싶을 때 봅니다.\n"
        "- 각 Figure는 `PNG preview → 확대 → 새 탭` 순으로 접근 가능합니다.",
    )

    # ---- 03 sample master ----
    add(
        "03_sample_master.html:sample-summary",
        "샘플 레벨 메타데이터 요약",
        "sample_master.tsv가 뭔가",
        "## sample_master.tsv\n"
        "- 모든 데이터셋을 통합한 **샘플 단위** 메타데이터 테이블입니다.\n"
        "- 하나의 행 = 하나의 조직 샘플, 열에는 `dataset`, `histology`, `molecular_subtype`, `platform`, `confidence` 등이 들어갑니다.\n\n"
        "## 조직학 아형 한 줄 설명\n"
        "- **cPTC**: Classical Papillary Thyroid Carcinoma. 가장 흔한 유두암, BRAF V600E가 많음.\n"
        "- **FVPTC**: Follicular Variant PTC. 여포 구조를 가진 PTC, RAS-like에 가까움.\n"
        "- **FTC**: Follicular Thyroid Carcinoma. 여포암, PTC와 구분되는 follicle 구조.\n"
        "- **ATC**: Anaplastic Thyroid Carcinoma. 미분화암, 매우 공격적.\n"
        "- **PDTC**: Poorly Differentiated Thyroid Carcinoma. 저분화암, PTC와 ATC 중간.\n"
        "- **MTC**: Medullary Thyroid Carcinoma. 수질성 갑상선암, C-cell 기원이라 follicular-derived와 섞으면 안 됩니다.\n\n"
        "## 주의할 점\n"
        "- `confidence` 컬럼이 낮은 샘플은 라벨이 확실하지 않다는 뜻이며, 핵심 모델 학습에서는 제외됩니다.",
    )
    add(
        "03_sample_master.html:sample-grid",
        "검색형 테이블 사용법",
        "필터와 검색 팁",
        "## 필터링\n"
        "- 상단 검색창에 gene symbol, sample ID, dataset ID 중 아무 키워드나 넣으면 즉시 좁혀집니다.\n"
        "- 컬럼 헤더를 클릭하면 그 컬럼 기준으로 정렬됩니다(내림/오름차순 전환).\n\n"
        "## 내보내기\n"
        "- 필터된 결과는 상단 `TSV` 버튼으로 그대로 다운로드할 수 있습니다.\n\n"
        "## 주의\n"
        "- 표 상단의 총 샘플 수 카운터는 **현재 필터된 상태**의 수입니다. 전체 수(1,509)와 혼동하지 마세요.",
    )
    add(
        "03_sample_master.html:sample-figs",
        "요약 도표 읽는 법",
        "histology × dataset 교차 뷰",
        "## 무엇을 보여주나\n"
        "- 각 데이터셋에서 어떤 조직학 아형이 얼마나 많은지 교차 집계한 도표입니다.\n"
        "- `confidence` 분포 히스토그램으로 라벨 신뢰도 전반을 확인할 수 있습니다.\n\n"
        "## 어떻게 읽는가\n"
        "- 특정 아형이 하나의 데이터셋에만 쏠려 있다면, 그 아형의 성능 일반화는 조심해야 합니다.",
    )

    # ---- 04 gene panels ----
    add(
        "04_gene_panels.html:panel-summary",
        "TDS16 / TierA67 / BRS71 / ThyroSeq112 패널 비교",
        "세 가지 연구 기반 유전자 패널",
        "## 각 패널의 정체\n"
        "- **TDS16 (Thyroid Differentiation Score, 16 genes)**: TCGA-THCA(Cell 2014) 논문에서 제안. 갑상선 세포가 얼마나 \"정상처럼 분화되었는지\"를 정량화합니다.\n"
        "- **BRS71 (BRAF-RAS Score, 71 genes)**: 같은 논문에서 BRAF vs RAS 축을 계산하는 proxy. 축 값이 양수면 BRAF-like, 음수면 RAS-like.\n"
        "- **TierA67**: 후속 연구에서 greedy feature selection으로 뽑힌 67개 유전자. 외부 재현성이 좋은 upper tier.\n"
        "- **ThyroSeq-112 surrogate**: 상용 `ThyroSeq v3` 패널(112 gene)을 공개 데이터에서 재현한 surrogate. 실제 상용 패널과 완전히 동일하지 않으므로 `surrogate`라고 부릅니다.\n\n"
        "## 왜 이걸 비교하나\n"
        "- 패널 크기가 클수록 성능은 조금씩 올라가지만 **측정 비용과 복잡도**도 커집니다. 의료현장에서 실효성은 ROC-AUC 외에도 비용 최적점에서 결정됩니다.",
    )
    add(
        "04_gene_panels.html:panel-viz",
        "패널 overlap / coverage 의미",
        "유전자 겹침과 커버리지 해석",
        "## coverage란\n"
        "- 해당 데이터셋의 플랫폼에서 패널 유전자가 **측정 가능**한 비율입니다.\n"
        "- 예: 112 gene 패널 중 microarray 플랫폼에서 실제 probe가 있는 유전자가 108개라면 coverage는 `108/112 = 96.4%`.\n\n"
        "## 왜 missing이 생기나\n"
        "- 오래된 microarray는 일부 유전자의 probe가 설계되지 않았습니다.\n"
        "- RNA-seq에서도 매우 짧거나 expression이 낮은 유전자는 드문드문 검출되지 않을 수 있습니다.\n\n"
        "## overlap (Venn/UpSet)\n"
        "- 여러 패널이 공통으로 포함하는 \"핵심\" 유전자를 시각화합니다. 세 패널 모두에 등장하는 유전자일수록 \"안정적인 분류 신호\"로 해석됩니다.",
    )
    add(
        "04_gene_panels.html:panel-benchmark",
        "패널 성능 벤치마크 미리보기",
        "성능 비교 요약",
        "## 빠르게 읽기\n"
        "- 16 gene ≈ BRS axis를 근사한 최소 패널, 112 gene ≈ 상업용 테스트 규모.\n"
        "- 두 극단 사이에서 AUC가 얼마나 늘고 얼마나 줄어드는지를 보여줍니다.\n\n"
        "## 주의할 점\n"
        "- 같은 AUC라도 외부 데이터셋에서 재현되는지는 `07 ML`, `08 패널 비교`에서 따로 확인해야 합니다.",
    )

    # ---- 05 EDA ----
    add(
        "05_eda.html:eda-guide",
        "PCA / UMAP / t-SNE를 신뢰할 때와 말 때",
        "차원 축소 기법 3총사",
        "## 차원 축소란\n"
        "- 유전자 수만 명(예: 2만)을 2차원으로 압축해서 샘플 간 유사성을 보는 작업입니다.\n\n"
        "## 각각 언제 쓰나\n"
        "- **PCA**: 선형 투영. 분산이 큰 축부터 고르므로 **전체적 구조**를 믿을 만하게 보여줍니다. Batch effect 진단에 제일 먼저 써야 합니다.\n"
        "- **UMAP**: 비선형. 국소 이웃 구조를 잘 보존하지만 축의 거리는 의미가 약합니다. 클러스터 개수/모양을 과신하지 마세요.\n"
        "- **t-SNE**: 비선형. perplexity에 민감하며, 멀리 떨어진 클러스터 간 거리는 믿지 말아야 합니다.\n\n"
        "## 언제 의심해야 하나\n"
        "- 같은 데이터를 다른 random seed로 돌렸을 때 UMAP 모양이 크게 바뀌면 그 구조는 해석하지 않는 게 안전합니다.\n"
        "- PCA에서 첫 축 분산이 60% 이상이면 그 축은 `batch` 또는 `platform`일 가능성이 높습니다.",
    )
    add(
        "05_eda.html:embeddings",
        "임베딩 공간 해석법",
        "클러스터를 어떻게 읽을까",
        "## 색깔이 의미하는 것\n"
        "- 점 하나 = 샘플 하나. 색깔은 `histology`, `dataset`, `molecular_subtype` 중 하나를 표현합니다.\n\n"
        "## 읽기 팁\n"
        "- `histology`로 칠했을 때 cPTC와 FVPTC가 살짝 섞이는 것은 정상입니다. 두 그룹은 연속체에 가깝습니다.\n"
        "- `dataset`으로 칠했을 때 데이터셋별로 완전히 분리된다면 **batch effect가 강하다는 신호**이므로 후속 harmonization이 필요합니다.\n\n"
        "## 주의\n"
        "- 낮은 차원 투영에서 \"떨어져 있으면 다르다\"는 참이지만 \"가까이 있으면 같다\"는 꼭 참이 아닙니다. ROC/AUC와 함께 판단하세요.",
    )

    # ---- 06 scores ----
    add(
        "06_scores.html:scores-intro",
        "점수 해석 프레임: TDS/BRS/dediff",
        "세 점수의 생물학적 의미",
        "## TDS (Thyroid Differentiation Score)\n"
        "- 정상 갑상선 세포가 가지는 고유 기능 유전자(TG, TPO, SLC5A5 등) 발현 강도의 평균.\n"
        "- **높을수록 분화 잘 되어 있음 = 정상에 가까움**. 낮을수록 dedifferentiation 진행.\n\n"
        "## BRS (BRAF-RAS Score)\n"
        "- 71개 유전자 집합에서 BRAF-like vs RAS-like 축에 투영한 값.\n"
        "- **양수(+)이면 BRAF-like, 음수(−)이면 RAS-like**.\n\n"
        "## Dedifferentiation score\n"
        "- 저분화/미분화 특징(예: EMT, 줄기세포성)이 얼마나 나타나는지를 보여주는 값.\n"
        "- 높을수록 ATC/PDTC에 가깝고 예후가 나쁜 경향.\n\n"
        "## 한 줄 요약\n"
        "- TDS ↓ + dediff ↑ 이면 공격적, BRS 부호가 subtype의 방향을 알려줍니다.",
    )
    add(
        "06_scores.html:scores-a",
        "분포(violin/box) 해석",
        "violin plot 처음 보는 경우",
        "## violin plot이 뭔가\n"
        "- box plot(중앙값/사분위)에 **확률밀도 모양**을 덧씌운 플롯. 한쪽으로 치우친 분포, 이중 모드 분포 등을 한눈에 보여줍니다.\n\n"
        "## 읽는 법\n"
        "- 허리가 두꺼운 부분 = 그 값 근처에 샘플이 많다.\n"
        "- 좌우 폭이 모두 길면 variance가 큼 = 서브그룹이 혼재했을 가능성.\n\n"
        "## 이 섹션에서\n"
        "- subtype 축(BRAF_like / RAS_like / normal / dediff)으로 비교하므로 서로 다른 subtype의 violin이 얼마나 겹치는지가 \"얼마나 구분하기 어려운가\"를 말해줍니다.",
    )
    add(
        "06_scores.html:scores-b",
        "2D / Radar 차트 읽기",
        "복합 지표 비교법",
        "## 2D scatter\n"
        "- 두 점수(예: TDS × BRS)를 축으로 놓아 각 샘플을 점으로 찍습니다.\n"
        "- 사분면(+,+ / +,− 등)별 분포를 보면 subtype의 생물학적 축을 이해할 수 있습니다.\n\n"
        "## Radar chart\n"
        "- 여러 점수(TDS, BRS, dediff, ...)를 축으로 돌려놓고 한 샘플/그룹의 **다차원 프로파일**을 polygon으로 그립니다.\n"
        "- 폴리곤의 크기/모양 비대칭성으로 \"어떤 축이 강한가\"를 빠르게 비교합니다.\n\n"
        "## 주의\n"
        "- Radar는 축 순서가 바뀌면 모양이 달라지므로 축 배열이 고정돼 있어야 공정한 비교가 됩니다.",
    )
    add(
        "06_scores.html:scores-linkout",
        "Embedding 연동 보기",
        "점수와 임베딩을 같이 보는 이유",
        "## 왜 연동하나\n"
        "- 점수 분포만 보면 \"수치는 비슷한데 실제 샘플은 전혀 다른 클러스터\"에 있는 경우를 놓칠 수 있습니다.\n"
        "- 2D 임베딩 위에 점수를 색으로 덧씌우면 \"어느 영역에서 BRAF-like 확신이 세진다\"를 눈으로 확인할 수 있습니다.\n\n"
        "## 팁\n"
        "- UMAP 좌표 위에 BRS를 heatmap처럼 칠했을 때 grad가 부드럽게 흐르면 그 축이 실제로 연속적 signal임을 뒷받침합니다.",
    )

    # ---- 07 ML baseline ----
    add(
        "07_ml_baseline.html:ml-nav",
        "섹션 빠른 이동",
        "ML 페이지 구성",
        "## 이 페이지 구성\n"
        "1. **Task 정의**: 무엇을 예측하는지.\n"
        "2. **곡선 기반 평가**: ROC/PR 곡선.\n"
        "3. **요약 성능**: 메트릭 표.\n"
        "4. **Threshold metrics**: 결정 임계값에 따른 성능.\n"
        "5. **Interactive explorer**: 슬라이더로 threshold를 직접 움직여보기.\n\n"
        "## 권장 순서\n"
        "- 이 페이지는 위에서 아래로 읽기만 해도 하나의 완결된 보고서입니다.",
    )
    add(
        "07_ml_baseline.html:ml-task",
        "이번 baseline의 정확한 정의",
        "무엇을 예측하는 문제인가",
        "## Task\n"
        "- 입력: 샘플 하나의 유전자 발현 프로파일(패널 내 유전자만).\n"
        "- 출력: **BRAF_like vs RAS_like** 이진 분류 확률.\n\n"
        "## 라벨은 어떻게 주어지나\n"
        "- TCGA-THCA는 공식 molecular subtype 주석이 있음 → 그대로 사용.\n"
        "- 외부 GEO는 논문/supplementary의 annotation 또는 histology proxy로 부여.\n"
        "- `unknown` / `normal` / `dediff` 샘플은 이번 baseline에서는 **제외**합니다.\n\n"
        "## 모델\n"
        "- 로지스틱 회귀(`LogReg_l2`) + 표준화. Baseline으로 의도적으로 단순화했습니다.\n"
        "- 5-fold stratified CV, stratification 기준은 dataset × label.\n\n"
        "## 왜 이렇게 단순한가\n"
        "- 단순한 모델이 **재현하기 쉽고, 외부 일반화 성능의 바닥을 정직하게 보여줍니다**. 복잡한 모델은 이후 단계에서 비교합니다.",
    )
    add(
        "07_ml_baseline.html:ml-curves",
        "ROC와 PR 곡선 완전 정복",
        "두 곡선이 보여주는 다른 진실",
        "## ROC curve\n"
        "- 가로축 **FPR (False Positive Rate)** = 음성을 잘못 양성이라고 부른 비율.\n"
        "- 세로축 **TPR (True Positive Rate, Recall)** = 진짜 양성을 양성이라고 부른 비율.\n"
        "- 임계값(threshold)을 0→1로 움직이며 점을 찍어 곡선을 만듭니다.\n\n"
        "## AUC 해석\n"
        "- `AUC = 0.5` → 완전 무작위.\n"
        "- `AUC = 1.0` → 완벽 분리.\n"
        "- `0.8~0.9` → 실용적, `>0.95` → 외부 검증에서 무너질 수 있으니 재현성 확인 필수.\n\n"
        "## PR curve (Precision-Recall)\n"
        "- 가로 **Recall**, 세로 **Precision**.\n"
        "- 양성 클래스가 **희소**할 때 ROC보다 더 정직한 그림을 줍니다.\n"
        "- PR-AUC = average precision.\n\n"
        "## 한 줄 규칙\n"
        "- 클래스 balance가 맞을 땐 ROC, 한쪽이 10% 이하로 드물 땐 PR을 우선 보세요.",
    )
    add(
        "07_ml_baseline.html:ml-summary",
        "요약 성능 표 한 줄 한 줄",
        "메트릭 정의와 읽는 순서",
        "## 각 메트릭\n"
        "- **ROC-AUC**: threshold 독립적 랭킹 품질.\n"
        "- **PR-AUC (AP)**: 희소 양성에서 랭킹 품질.\n"
        "- **Balanced accuracy**: `(TPR + TNR) / 2`. 클래스 불균형이 있어도 공정.\n"
        "- **F1**: precision과 recall의 조화 평균, 0~1.\n"
        "- **MCC (Matthews Correlation Coefficient)**: −1~1, **클래스 불균형에 가장 강건**한 메트릭.\n"
        "- **Brier score**: `mean((p − y)^2)`. **낮을수록** 잘 calibrated. 0.25 = 무의미한 예측, 0.0 = 완벽.\n\n"
        "## 읽는 순서\n"
        "1. ROC-AUC로 전반 순위 확인.\n"
        "2. Balanced accuracy/MCC로 threshold 부근 품질 확인.\n"
        "3. Brier로 확률 자체가 믿을 만한지 확인.\n\n"
        "## 주의\n"
        "- 단일 fold의 값에 과하게 의존하지 말고 CV 평균 ± SD를 같이 보세요.",
    )
    add(
        "07_ml_baseline.html:threshold-table",
        "Threshold metrics 왜 중요한가",
        "결정 임계값이 성능을 결정합니다",
        "## threshold란\n"
        "- 모델이 출력하는 확률 `p ∈ [0,1]`을 \"양성\"으로 부르는 기준값입니다.\n"
        "- 기본값은 0.5이지만 **이게 최적이라는 보장은 없습니다**.\n\n"
        "## Youden J\n"
        "- `J = TPR − FPR`이 최대가 되는 threshold.\n"
        "- 민감도와 특이도의 균형점으로, 의료 스크리닝에서 자주 쓰입니다.\n\n"
        "## Balanced accuracy가 threshold에 따라 반전하는 이유\n"
        "- threshold가 너무 낮으면 모두를 양성으로 불러 TPR=1이지만 TNR=0.\n"
        "- 너무 높으면 반대로 TNR=1이지만 TPR=0.\n"
        "- 최적값은 보통 0.3~0.6 구간에서 중간 어딘가에 있습니다.\n\n"
        "## 표 읽는 팁\n"
        "- 같은 AUC라도 임상적 의사결정은 threshold를 어떻게 잡느냐에 따라 크게 달라진다는 점이 핵심입니다.",
    )
    add(
        "07_ml_baseline.html:threshold-explorer",
        "Interactive threshold explorer 쓰는 법",
        "슬라이더로 직접 확인하기",
        "## 단계별 사용법\n"
        "1. 모델과 데이터셋을 드롭다운에서 선택합니다.\n"
        "2. **슬라이더를 움직이면** 실시간으로 다음 값이 갱신됩니다: TPR, FPR, Precision, F1, Balanced accuracy, MCC.\n"
        "3. confusion matrix 타일이 색을 바꿔 잘/틀림 분포를 보여줍니다.\n\n"
        "## 무엇을 확인하면 좋나\n"
        "- threshold를 조금만 밀어도 F1/MCC가 크게 흔들린다면, 모델 확률이 calibrated되지 않았을 가능성이 있습니다.\n"
        "- 임상 목적이 \"양성 놓치지 않기\"라면 TPR 우선, \"헛 양성 최소화\"라면 Precision 우선으로 threshold를 골라야 합니다.\n\n"
        "## 팁\n"
        "- 슬라이더 아래 숫자 입력란에 값을 직접 넣을 수도 있습니다. 0.5가 기본이지만 Youden J에서 제안한 값을 넣어 비교해 보세요.",
    )

    # ---- 08 panel comparison ----
    add(
        "08_panel_comparison.html:perf",
        "성능 vs 패널 크기 tradeoff",
        "작은 패널 vs 큰 패널",
        "## 기본 경향\n"
        "- 패널 크기를 16 → 32 → 67 → 112로 키우면 **AUC가 단조 증가**하는 게 일반적입니다.\n"
        "- 그러나 증가폭은 빠르게 줄어듭니다(**diminishing returns**).\n\n"
        "## Occam의 면도날\n"
        "- 같은 성능이면 유전자 수가 적은 패널이 임상에 더 좋습니다(비용↓, QC 용이, 실패율↓).\n"
        "- 따라서 \"AUC 1% 더\"가 유전자 50개 더를 정당화하는지 따져봐야 합니다.\n\n"
        "## 읽는 팁\n"
        "- x축=패널 크기, y축=AUC인 곡선에서 **기울기가 급격히 완만해지는 knee point**가 실무적 최적점입니다.",
    )
    add(
        "08_panel_comparison.html:gain",
        "증분 정보량 (incremental gain)",
        "유전자 하나를 더하면 얼마나 얻나",
        "## 개념\n"
        "- 패널에 유전자를 하나씩 더할 때마다 **AUC가 얼마나 오르는지** 측정한 값.\n"
        "- 수학적으로는 information gain / mutual information과 유사합니다.\n\n"
        "## 읽기 팁\n"
        "- 상위에 나오는 유전자는 **서로 상관이 낮고, 라벨과는 상관이 높은** 유전자입니다.\n"
        "- 뒤로 갈수록 gain이 0에 수렴하는 이유는 새 유전자가 기존 유전자의 신호와 겹치기 때문입니다(redundancy).\n\n"
        "## 주의\n"
        "- gain 순위는 feature selection 알고리즘에 의존합니다. greedy vs L1 vs 상관 기반은 다른 순위를 줄 수 있습니다.",
    )

    # ---- 09 SHAP ----
    add(
        "09_shap.html:shap-summary",
        "SHAP 값 직관적 이해",
        "이 feature가 예측을 얼마나 밀었나",
        "## SHAP이 뭔가\n"
        "- **SHapley Additive exPlanations**. 게임이론의 Shapley value를 ML 해석에 적용.\n"
        "- 직관: \"이 feature가 이 예측을 base probability로부터 위로/아래로 얼마나 밀었는가\"를 **항상 합하면 실제 예측과 정확히 일치**하도록 나눠줍니다.\n\n"
        "## 전역 vs 국소\n"
        "- **전역(global)**: 모든 샘플에서 평균적으로 중요한 feature 순위.\n"
        "- **국소(local)**: 특정 샘플 한 개의 예측이 왜 그렇게 나왔는지 분해.\n\n"
        "## 주의\n"
        "- SHAP은 \"모델이 그렇게 판단했다\"는 설명이지 **생물학적 인과관계의 증거가 아닙니다**.\n"
        "- 상관관계가 높은 feature들끼리는 SHAP이 임의로 몫을 나누므로 개별 유전자 순위에 과하게 의존하지 마세요.",
    )
    add(
        "09_shap.html:shap-global",
        "전역 중요도 bar chart 읽기",
        "상위 feature bar plot",
        "## 그림 형식\n"
        "- 세로축 feature 이름, 가로축 평균 |SHAP| 값.\n"
        "- 막대가 길수록 평균적으로 예측을 많이 밀고 당긴 feature입니다.\n\n"
        "## 생물학적 감 잡기\n"
        "- 상위권에 THYROID DIFFERENTIATION 관련 유전자(TG, TPO, SLC5A5)가 올라오면 TDS 신호를 주로 보고 있다는 뜻입니다.\n"
        "- BRAF 하위 경로 유전자(DUSP, FOSL 등)가 상위에 있으면 BRAF 축을 반영.\n\n"
        "## 주의\n"
        "- 순위는 패널 구성에 따라 크게 바뀝니다. panel 간 비교는 `08` 페이지를 함께 보세요.",
    )
    add(
        "09_shap.html:shap-waterfall",
        "샘플별 워터폴 플롯",
        "한 샘플의 예측이 왜 그렇게 나왔나",
        "## 읽는 법\n"
        "- 왼쪽: base value (모든 샘플의 평균 예측 로짓).\n"
        "- 가운데: 각 feature가 **위(+)**로 또는 **아래(−)**로 얼마나 예측을 밀었는지 단계별 누적.\n"
        "- 오른쪽: 최종 예측 확률.\n\n"
        "## 예시\n"
        "- `BRAF downstream gene ↑ → +0.3 push`, `TG ↓ → +0.1 push` 이런 식으로 합해져서 최종 p=0.87 같은 값에 도달.\n\n"
        "## 활용\n"
        "- 개별 환자에게 \"왜 BRAF-like로 분류됐나\"를 설명할 때 쓸 수 있습니다(단, 최종 임상 판단의 근거로 단독 사용은 금물).",
    )

    # ---- 10 gene explorer ----
    add(
        "10_gene_explorer.html:gene-intro",
        "유전자 발현 violin plot 읽기",
        "개별 유전자의 cohort별 분포",
        "## 무엇을 보여주나\n"
        "- 원하는 유전자의 발현을 **조직학 아형 × 코호트**로 쪼개서 분포 모양으로 비교합니다.\n\n"
        "## 어떻게 읽는가\n"
        "- 같은 유전자라도 데이터셋에 따라 절대값이 달라질 수 있습니다(플랫폼 차이).\n"
        "- 따라서 **같은 플롯 내 상대 비교**는 의미 있고, **플롯 간 절대값 비교**는 조심해야 합니다.\n\n"
        "## 드러나는 것\n"
        "- 특정 유전자가 BRAF_like에서 일관되게 높고 RAS_like에서 낮다면, 이 유전자는 `08 패널 비교`의 상위권에 있을 가능성이 높습니다.",
    )
    add(
        "10_gene_explorer.html:gene-search",
        "유전자 검색 사용법",
        "심볼과 별칭으로 찾기",
        "## 사용법\n"
        "- 검색창에 `TG`, `TPO`, `BRAF` 같은 공식 gene symbol을 입력하면 자동완성됩니다.\n"
        "- 별칭(alias)도 지원. 예: `NIS → SLC5A5`로 자동 매핑.\n\n"
        "## 팁\n"
        "- 한 번에 여러 유전자를 쉼표로 구분해 넣으면 작은 multi-panel 비교를 볼 수 있습니다.",
    )
    add(
        "10_gene_explorer.html:gene-default",
        "기본 유전자 Figure",
        "추천 유전자 미리 준비",
        "## 왜 기본값이 있나\n"
        "- 처음 방문한 사람이 \"어떤 유전자부터 봐야 하지?\" 고민하지 않도록, TDS/BRS 핵심 유전자 몇 개를 미리 보여줍니다.\n\n"
        "## 포함된 유전자 예시\n"
        "- **TG (Thyroglobulin)**: 갑상선 분화 표지자의 대표. 정상>PTC>ATC 순으로 감소.\n"
        "- **TPO**: 갑상선 과산화효소. 정상 조직 특이적.\n"
        "- **BRAF**: 자체 발현이 아니라 downstream 경로 유전자와 함께 봐야 subtype 신호가 보입니다.",
    )
    add(
        "10_gene_explorer.html:gene-note",
        "섹션 설명",
        "이 페이지의 한계",
        "## 한계\n"
        "- 이 페이지는 **단일 유전자 시각화** 도구이며, 다변량 모델의 결정은 `07 ML`과 `09 SHAP`을 함께 봐야 합니다.\n"
        "- 발현이 높다고 \"그 유전자의 돌연변이가 있다\"는 뜻이 아닙니다. 발현과 돌연변이는 별개 층의 데이터입니다.",
    )

    # ---- 11 cohort compare ----
    add(
        "11_cohort_compare.html:cohort-text",
        "코호트 비교 원칙: batch effect와 harmonization",
        "여러 데이터셋을 어떻게 공정하게 비교하나",
        "## batch effect란\n"
        "- 실험실, 장비, 시약, 일자에 따라 같은 생물학적 샘플이라도 측정값이 달라지는 현상.\n"
        "- 이걸 잡지 않으면 \"코호트가 다르다\"가 \"subtype이 다르다\"로 잘못 읽힙니다.\n\n"
        "## harmonization\n"
        "- batch effect를 통계적으로 제거하는 과정.\n"
        "- 대표 기법: **ComBat**(empirical Bayes 기반), **limma removeBatchEffect**, **quantile normalization**.\n\n"
        "## 이 대시보드에서\n"
        "- 주된 모델은 **cross-dataset 학습 없이** 각 데이터셋 내부에서만 CV를 돌립니다. harmonization 과정에서 신호가 깎일 위험을 피하려는 의도입니다.",
    )
    add(
        "11_cohort_compare.html:cohort-toggle",
        "분포 축 전환 (by dataset / histology / label)",
        "같은 데이터 다르게 보기",
        "## 세 가지 축\n"
        "- **by dataset**: 코호트별 분포. 기술적 차이(플랫폼/batch)를 본다.\n"
        "- **by histology**: cPTC/FVPTC/FTC/ATC 등 조직학 차이. 임상 의미가 큰 차이.\n"
        "- **by label**: BRAF_like vs RAS_like. 모델이 배우는 축.\n\n"
        "## 읽기 팁\n"
        "- by histology에서는 선명한 차이가 있는데 by dataset으로 바꾸면 무너진다 → **일반화 가능한 신호**.\n"
        "- by dataset만으로도 선명한 차이가 있다 → batch-driven 가능성, 모델이 이걸 학습하면 외부에서 실패합니다.",
    )
    add(
        "11_cohort_compare.html:cohort-figs",
        "코호트 비교 뷰",
        "side-by-side 보는 법",
        "## 보는 순서\n"
        "1. 왼쪽 코호트와 오른쪽 코호트의 subtype 비율을 비교.\n"
        "2. 같은 subtype 내에서도 발현 분포가 겹치는지 확인.\n"
        "3. 겹치면 모델이 cross-dataset에서도 잘 돌아갈 가능성↑.",
    )
    add(
        "11_cohort_compare.html:cohort-support",
        "보조 구조 도표",
        "추가 맥락 그림",
        "## 언제 참고하나\n"
        "- 메인 비교에서 이상한 점을 발견했을 때 배경 분포, 결측 패턴을 확인하는 용도입니다.",
    )

    # ---- 12 business ----
    add(
        "12_business.html:biz-context",
        "TAM / SAM / SOM과 NPV 기본",
        "사업성 용어 15분 정리",
        "## 용어\n"
        "- **TAM (Total Addressable Market)**: 이론상 시장 전체 크기.\n"
        "- **SAM (Serviceable Addressable Market)**: 우리 제품이 실제로 닿을 수 있는 부분.\n"
        "- **SOM (Serviceable Obtainable Market)**: 단기간(3~5년)에 확보 가능한 부분.\n"
        "- **NPV (Net Present Value)**: 미래 현금흐름을 현재가치로 할인한 합. 양수면 투자 가치가 있다고 봅니다.\n\n"
        "## pricing NPV\n"
        "- 테스트당 가격 × 연간 검사 수 × 채택률 − 운영비 를 연도별로 할인(보통 8~12%)해 합산.\n\n"
        "## 주의\n"
        "- 이 페이지 숫자는 **가정에 민감**합니다. 가격이 20% 움직이면 NPV는 훨씬 더 크게 움직입니다.\n"
        "- 최종 사업 결정의 근거가 아니라 \"상대적 우위/열위\" 비교용으로 봐주세요.",
    )
    add(
        "12_business.html:biz-figs",
        "비즈니스 프레임 Figure",
        "가정을 시각화",
        "## 무엇을 보여주나\n"
        "- 가격·채택률·검사량의 서로 다른 가정 조합에서 3년 누적 매출·NPV가 어떻게 달라지는지 시나리오별로 보여줍니다.\n\n"
        "## 읽는 팁\n"
        "- 한 가정에 모든 결과가 몰려 있다면 그 가정이 가장 민감한 레버입니다.\n"
        "- 시나리오 간 차이가 크다면 **가정의 근거가 불확실**하다는 뜻이므로 추가 근거가 필요합니다.",
    )
    add(
        "12_business.html:biz-table",
        "경쟁 구도 표",
        "시장 내 비교",
        "## 비교 축\n"
        "- 테스트당 가격, 운영 friction, 필요 검체량, turnaround, evidence level 등.\n\n"
        "## 주의\n"
        "- 경쟁사 숫자는 공개 자료 기반이라 최신이 아닐 수 있습니다. 의사결정 직전에 반드시 최신 가격표/임상 근거를 재확인하세요.",
    )
    add(
        "12_business.html:biz-anchor",
        "과학적 anchor",
        "사업성 주장을 과학적으로 묶기",
        "## 왜 anchor가 필요한가\n"
        "- 사업성 주장은 \"AUC 0.93\"처럼 모델 성능에 묶여야 신뢰가 갑니다.\n"
        "- 이 섹션은 `07 ML` 결과와 이 페이지의 NPV를 한 줄에 연결합니다.",
    )

    # ---- 13 reports ----
    add(
        "13_reports.html:reports-intro",
        "문서 읽기 가이드",
        "각 markdown 문서가 언제 쓰이나",
        "## 이 페이지 구성\n"
        "- 분석 파이프라인에서 생성된 **원문 markdown 리포트**들을 브라우저에서 바로 렌더링합니다.\n\n"
        "## 추천 순서\n"
        "1. `analysis_summary.md` → 전체 방법론 요약.\n"
        "2. `ml_report_v2.md` → 모델 성능 상세.\n"
        "3. `biomarker_analysis.md` → 생물학적 해석.\n"
        "4. `next_steps.md` → 다음 할 일.\n\n"
        "## 팁\n"
        "- markdown 내 링크는 상대 경로이므로 이 페이지에서 바로 클릭 가능합니다.",
    )
    add(
        "13_reports.html:reports-md",
        "Rendered markdown",
        "랜더링된 문서",
        "## 렌더링 엔진\n"
        "- `marked.js`(클라이언트)로 렌더링합니다.\n"
        "- 수식은 지원하지 않으므로 `$...$` 표기는 원본 그대로 보일 수 있습니다.",
    )

    # ---- 14 caveats ----
    add(
        "14_caveats.html:risk",
        "리스크 레지스터 읽는 법",
        "각 리스크 항목이 뭐야",
        "## 리스크 등급\n"
        "- **High**: 지금 결론을 뒤집을 수 있는 리스크.\n"
        "- **Medium**: 성능 숫자를 조정해야 하는 리스크.\n"
        "- **Low**: 문서화 개선으로 커버되는 리스크.\n\n"
        "## 대표 항목\n"
        "- *batch effect* — 데이터셋별 측정 편차.\n"
        "- *label leakage* — 학습 셋에 외부 정보가 섞였을 가능성.\n"
        "- *survivorship bias* — QC에서 걸러진 샘플이 한쪽 subtype에 몰려 있을 가능성.\n\n"
        "## 주의\n"
        "- 리스크를 읽고도 **완화 계획**을 확인하지 않으면 단순 공포 조성에 그칩니다. 각 항목의 mitigation 컬럼을 함께 보세요.",
    )
    add(
        "14_caveats.html:failure",
        "Known failure modes",
        "이미 발견한 실패 시나리오",
        "## 예시\n"
        "- *microarray 전용 모델을 RNA-seq 데이터에 적용* → 스케일 불일치로 AUC 급락.\n"
        "- *normal 샘플 비중이 극단적으로 낮은 코호트* → baseline 분포가 왜곡되어 threshold가 엉뚱한 위치로.\n"
        "- *MTC 샘플을 follicular-derived와 섞어 학습* → MTC가 outlier로 가중치를 잡아먹음.\n\n"
        "## 활용\n"
        "- 새 데이터셋을 추가할 때 이 목록을 점검표처럼 쓰면 같은 함정을 피할 수 있습니다.",
    )
    add(
        "14_caveats.html:caveats-context",
        "한계 해석용 구조 도표",
        "배경 시각화",
        "## 무엇을 보여주나\n"
        "- 리스크/실패 모드의 배경이 되는 데이터 구조를 한 번에 보여주는 보조 그림입니다.\n\n"
        "## 읽기 팁\n"
        "- 이 그림 자체가 결론이 아니라 **앞서 읽은 리스크 항목을 시각적으로 확인**하는 용도입니다.",
    )

    # ---- 15 drug discovery ----
    add(
        "15_drug_discovery.html:hero",
        "Target-to-Drug 파이프라인 5단계",
        "신약 개발 단계 한눈 정리",
        "## 5단계가 뭔가\n"
        "1. **Biomarker**: 질병 표지자/표적 후보 발굴. 이 대시보드의 바이오마커 분석이 여기에 해당.\n"
        "2. **Target validation**: 그 표적을 건드리면 정말 질병이 좋아지는지 생물학적으로 검증(세포/동물 실험).\n"
        "3. **Hit discovery**: 표적에 붙는 화합물을 수만~수백만 개 스크리닝해 `hit` 선별.\n"
        "4. **Lead optimization**: hit의 약물성(흡수/분포/대사/배설/독성)을 개선해 `lead`로 발전.\n"
        "5. **Preclinical**: 동물 모델에서 안전성/효능 확인, 임상 진입 준비.\n\n"
        "## 이 페이지의 범위\n"
        "- 위 단계 중 **1~3의 초반**(in-silico triage)에 해당. 실험 단계는 포함하지 않습니다.\n\n"
        "## 주의\n"
        "- 컴퓨터 기반 pre-screening 결과는 가설이지 증거가 아닙니다. 실제 약효는 실험으로 확인해야 합니다.",
    )
    add(
        "15_drug_discovery.html:row1",
        "타겟 분류: Novel / Emerging / Validated",
        "왼쪽 필터 칩이 의미하는 것",
        "## 분류 기준\n"
        "- **Novel**: 갑상선암에서 표적 근거가 거의 없는 신규. 발견 가치 높지만 실패 위험 큼.\n"
        "- **Emerging**: 초기 임상/비교 연구 존재. 중간 위험.\n"
        "- **Validated**: 승인약 또는 후기 임상 존재. 상대적으로 안정.\n\n"
        "## 오른쪽 metrics\n"
        "- `compounds`: ChEMBL에 연결된 화합물 수.\n"
        "- `best pChEMBL`: 가장 강한 활성값 (`9 = 1 nM`, `7 = 100 nM`).\n"
        "- `thyroid papers`: 이 표적을 갑상선 맥락에서 언급한 PubMed 논문 수.\n"
        "- `novelty score`: literature + mutation data 기반 신규성 점수.\n\n"
        "## SMILES 미리보기\n"
        "- 오른쪽 canvas에 SMILES 문자열을 2D 구조로 렌더링합니다. 화합물을 하나 클릭하면 해당 구조가 나타납니다.",
    )
    add(
        "15_drug_discovery.html:row-network",
        "Target ↔ Compound 관계망 읽기",
        "네트워크 그래프 해석",
        "## 노드 색\n"
        "- teal(청록) 큰 원 = 표적 유전자.\n"
        "- cyan(하늘) 원 = 대표 화합물.\n\n"
        "## 엣지(연결)\n"
        "- 두 노드를 잇는 선은 ChEMBL에 등록된 **bioactivity 근거**를 의미.\n"
        "- 원 크기는 pChEMBL에 비례. 큰 원일수록 강한 활성.\n\n"
        "## 읽는 팁\n"
        "- 한 화합물이 여러 표적에 연결되면 `polypharmacology` 후보이며, 부작용 위험도 함께 봐야 합니다.",
    )
    add(
        "15_drug_discovery.html:row-table",
        "후보 화합물 전체 표",
        "테이블 컬럼 의미",
        "## 컬럼\n"
        "- `ChEMBL ID`: ChEMBL 데이터베이스의 고유 식별자.\n"
        "- `pChEMBL`: 활성값 로그 스케일. 7 이상이 일반적 기대선.\n"
        "- `targets`: 이 화합물이 붙는 표적 수(많으면 다중성).\n"
        "- `phase`: 최고 임상 단계(0=preclinical, 4=approved).\n\n"
        "## 정렬\n"
        "- 컬럼 헤더를 클릭하면 정렬. `pChEMBL` 내림차순으로 보면 강력한 후보부터 나타납니다.",
    )
    add(
        "15_drug_discovery.html:row-lit",
        "문헌 근거 스트립",
        "PubMed 카드 읽는 법",
        "## 각 카드\n"
        "- 연도, 제목, 근거 문장, PMID 링크.\n\n"
        "## 체크포인트\n"
        "- 최근 5년 논문 비중 → 주제 활성도.\n"
        "- 제목에 표적 + 갑상선 키워드 동시 포함 → 관련성 높음.\n"
        "- 리뷰와 원저 논문의 비율 → 원저가 적으면 실증이 아직 부족.",
    )

    return content, generic_template


# ---------------------------------------------------------------------------
# Glossary page content: term -> {summary, body_md}
# ---------------------------------------------------------------------------


def _build_glossary() -> list[dict]:
    """Each entry: {id, term, summary, body_md, category}."""
    entries: list[dict] = []

    def E(cat: str, term: str, summary: str, body: str) -> None:
        entries.append(
            {
                "id": "term-" + re.sub(r"[^a-z0-9]+", "-", term.lower()).strip("-"),
                "term": term,
                "summary": summary,
                "body_md": body,
                "category": cat,
            }
        )

    # --- 모델 평가 ---
    E(
        "모델 평가",
        "AUC",
        "ROC 곡선 아래 면적",
        "## 정의\n"
        "- ROC curve 아래의 면적(0~1).\n"
        "- `0.5`는 무작위, `1.0`은 완벽 분리.\n\n"
        "## 실무 해석\n"
        "- `0.7~0.8` 쓸만함, `0.8~0.9` 좋음, `>0.95` 외부 재현성 점검 필수.\n\n"
        "## 주의\n"
        "- 클래스 불균형이 심하면 AUC가 과대평가될 수 있어 PR-AUC를 같이 봐야 합니다.",
    )
    E(
        "모델 평가",
        "ROC Curve",
        "FPR vs TPR 곡선",
        "## 축\n"
        "- x축: **FPR** (False Positive Rate).\n"
        "- y축: **TPR** (True Positive Rate = Recall).\n\n"
        "## 곡선 그리는 방법\n"
        "- threshold를 0→1로 움직이며 각 점의 (FPR, TPR)을 찍습니다.\n\n"
        "## 해석\n"
        "- 좌상단에 가까울수록 좋음. 대각선은 무작위.",
    )
    E(
        "모델 평가",
        "PR Curve / PR-AUC",
        "Precision-Recall 곡선과 AP",
        "## 정의\n"
        "- 가로 Recall, 세로 Precision인 곡선.\n"
        "- **PR-AUC = Average Precision (AP)**.\n\n"
        "## 언제 쓰나\n"
        "- 양성 클래스가 드물 때(예: 10% 이하). ROC보다 정직한 그림을 줍니다.",
    )
    E(
        "모델 평가",
        "Confidence Interval (95%)",
        "95% 신뢰구간",
        "## 정의\n"
        "- 동일 실험을 무한 반복했을 때 **95%의 경우 이 구간 안에 참값이 들어간다**는 의미(빈도주의적 해석).\n\n"
        "## 실무\n"
        "- AUC 0.92 (95% CI 0.88–0.95)처럼 표기.\n"
        "- 구간이 넓으면 샘플 수가 적거나 분산이 크다는 신호.",
    )
    E(
        "모델 평가",
        "Bootstrap",
        "재표집으로 CI 추정",
        "## 방법\n"
        "- 원래 데이터에서 **복원 추출**로 N개 표본을 만든다. 이걸 1000번 반복해 통계량 분포를 구한다.\n\n"
        "## 장점\n"
        "- 분포 가정 없이 CI를 구할 수 있다.\n\n"
        "## 주의\n"
        "- 표본 자체가 편향됐다면 bootstrap CI도 편향됩니다.",
    )
    E(
        "모델 평가",
        "Youden J",
        "민감도+특이도-1이 최대인 threshold",
        "## 정의\n"
        "- `J = TPR − FPR` = 민감도 + 특이도 − 1.\n\n"
        "## 사용\n"
        "- `J`가 최대가 되는 threshold를 선택하면 민감도와 특이도 균형이 잡힙니다.\n\n"
        "## 의료 맥락\n"
        "- 스크리닝 목적에서 기본 임계값으로 자주 쓰입니다.",
    )
    E(
        "모델 평가",
        "Brier Score",
        "확률 예측의 제곱오차",
        "## 정의\n"
        "- `Brier = mean((p_i − y_i)^2)`, y∈{0,1}.\n\n"
        "## 해석\n"
        "- **낮을수록 좋음**. 0=완벽, 0.25=무의미 예측(모두 0.5 예측).\n\n"
        "## 왜 중요한가\n"
        "- AUC만 높으면 순위만 맞다는 뜻. Brier가 낮아야 확률값 자체가 의미 있습니다(calibration).",
    )
    E(
        "모델 평가",
        "Cohen's d",
        "두 그룹 간 효과 크기",
        "## 정의\n"
        "- `d = (μ1 − μ2) / pooled_SD`.\n\n"
        "## 해석 규칙\n"
        "- `0.2` 작음, `0.5` 중간, `0.8` 큼.\n\n"
        "## 활용\n"
        "- violin plot에서 subtype 간 차이의 크기를 한 숫자로 요약할 때.",
    )
    E(
        "모델 평가",
        "FDR / BH correction",
        "다중 비교 보정",
        "## 문제\n"
        "- 유전자 2만 개를 t-test로 각각 비교하면 우연히 유의한 것이 많이 나온다(`p<0.05`에서도 5%가 false positive).\n\n"
        "## Benjamini-Hochberg\n"
        "- p-value를 정렬한 뒤 rank 기반으로 보정해 **기대 FDR(False Discovery Rate)**를 제어.\n\n"
        "## 실무\n"
        "- 보통 `FDR < 0.05` 또는 `0.10` 기준으로 필터링합니다.",
    )
    E(
        "모델 평가",
        "Cross-Validation (k-fold, stratified, nested)",
        "교차검증 유형",
        "## k-fold\n"
        "- 데이터를 k개로 나눠 (k-1)에서 학습, 1에서 평가. 이걸 k번 돌려 평균.\n\n"
        "## stratified\n"
        "- fold별로 클래스 비율을 원본과 같게 유지. 불균형 데이터에서 필수.\n\n"
        "## nested\n"
        "- 하이퍼파라미터 탐색(inner fold)과 최종 평가(outer fold)를 분리. **정직한 성능 추정**을 위해 필요.",
    )
    E(
        "모델 평가",
        "Class Imbalance / Balanced Accuracy",
        "불균형과 균형 정확도",
        "## class imbalance\n"
        "- 한 클래스가 다른 클래스보다 훨씬 많은 상황(예: 85% vs 15%).\n"
        "- 일반 accuracy는 다수 클래스를 외워만 해도 0.85가 돼 무의미.\n\n"
        "## balanced accuracy\n"
        "- `(TPR + TNR) / 2`. 두 클래스에 동등 가중치.\n\n"
        "## 추가 팁\n"
        "- MCC, F1, ROC-AUC, PR-AUC와 함께 해석하는 것이 안전합니다.",
    )
    E(
        "모델 평가",
        "Decision Threshold",
        "양성 판정 기준값",
        "## 정의\n"
        "- 모델 확률 `p`에 대해 `p >= τ`이면 양성이라고 부르는 `τ` 값.\n\n"
        "## 기본값 0.5의 함정\n"
        "- 클래스 불균형이나 비용 비대칭이 있으면 0.5는 최적이 아닙니다.\n\n"
        "## 실무\n"
        "- Youden J, 비용 함수 최소화, 또는 정책 기반으로 결정.",
    )
    E(
        "모델 평가",
        "Feature Selection",
        "유전자(특징) 선택 기법",
        "## 종류\n"
        "- **filter**: t-test, mutual information 등으로 각 유전자 점수 매김.\n"
        "- **wrapper**: greedy forward/backward로 실제 모델 성능 기준 탐색.\n"
        "- **embedded**: L1 regression, tree importance.\n\n"
        "## 주의\n"
        "- feature selection을 CV 바깥에서 하면 **leakage**가 생깁니다. nested CV 내부에서 수행해야 정직합니다.",
    )

    # --- 분자 용어 ---
    E(
        "생물학 (분자)",
        "BRAF V600E",
        "PTC에서 가장 흔한 돌연변이",
        "## 정의\n"
        "- BRAF 유전자 600번 아미노산이 Valine(V)에서 Glutamate(E)로 바뀌는 점 돌연변이.\n\n"
        "## 임상\n"
        "- cPTC의 40~60%에서 발견. BRAF-like subtype의 핵심 driver.\n"
        "- vemurafenib, dabrafenib 같은 BRAF inhibitor의 표적.\n\n"
        "## 주의\n"
        "- V600E가 있어도 모두 공격적이지는 않습니다. TDS, dediff score와 함께 판단.",
    )
    E(
        "생물학 (분자)",
        "NRAS / HRAS / KRAS",
        "RAS family 세 멤버",
        "## 정의\n"
        "- 세포 분열 신호전달을 시작하는 작은 GTPase 세 종류.\n\n"
        "## 갑상선암에서\n"
        "- NRAS Q61R이 대표적. FVPTC, FTC에서 흔함.\n"
        "- RAS-like subtype의 driver.\n\n"
        "## 치료\n"
        "- RAS 직접 표적은 오래 어려웠으나 최근 KRAS G12C inhibitor가 등장.",
    )
    E(
        "생물학 (분자)",
        "TDS (Thyroid Differentiation Score)",
        "갑상선 분화 점수",
        "## 정의\n"
        "- 정상 갑상선 기능 유전자(TG, TPO, SLC5A5 등) 발현 평균.\n\n"
        "## 해석\n"
        "- **높으면 정상에 가까움**, 낮으면 dedifferentiation.\n\n"
        "## 출처\n"
        "- TCGA-THCA 논문(Cell, 2014).",
    )
    E(
        "생물학 (분자)",
        "BRS (BRAF-RAS Score)",
        "BRAF vs RAS 축 점수",
        "## 정의\n"
        "- 71개 유전자에서 BRAF-like와 RAS-like 중심의 차이를 반영한 축.\n\n"
        "## 해석\n"
        "- **양수(+) → BRAF-like**, **음수(−) → RAS-like**.\n\n"
        "## 용도\n"
        "- molecular subtype 연속 지표로 쓰입니다. 이 대시보드 전반의 핵심 축.",
    )
    E(
        "생물학 (조직학)",
        "PTC / cPTC / FVPTC",
        "유두암과 주요 변종",
        "## PTC (Papillary Thyroid Carcinoma)\n"
        "- 갑상선암 중 가장 흔함(~80%).\n\n"
        "## cPTC (Classical PTC)\n"
        "- PTC의 classical 변종, 유두 구조 명확. BRAF V600E가 많음.\n\n"
        "## FVPTC (Follicular Variant PTC)\n"
        "- 유두 구조 대신 여포 구조. RAS-like에 가깝고 예후는 cPTC와 비슷하거나 더 좋음.",
    )
    E(
        "생물학 (조직학)",
        "FTC / ATC / PDTC / MTC",
        "기타 갑상선암 종류",
        "## FTC\n"
        "- Follicular Thyroid Carcinoma. 여포암.\n\n"
        "## ATC\n"
        "- Anaplastic Thyroid Carcinoma. 미분화암. 극히 공격적, 생존율 낮음.\n\n"
        "## PDTC\n"
        "- Poorly Differentiated. 저분화, PTC와 ATC 중간.\n\n"
        "## MTC\n"
        "- Medullary Thyroid Carcinoma. 수질성. **C-cell 기원**이라 follicular-derived와 별도.",
    )
    E(
        "생물학 (일반)",
        "Gene Expression",
        "유전자 발현",
        "## 정의\n"
        "- DNA에서 mRNA가 전사되고 단백질로 번역되는 과정의 강도.\n"
        "- 실무에서는 **mRNA 양**을 측정한 값(RNA-seq count, microarray intensity)을 말합니다.\n\n"
        "## 단위\n"
        "- TPM, FPKM (RNA-seq), log2 intensity (microarray).",
    )
    E(
        "생물학 (일반)",
        "mRNA / RNA-seq / Microarray",
        "전사체 측정 기술",
        "## mRNA\n"
        "- 단백질을 만들기 위해 DNA에서 복사된 중간 분자.\n\n"
        "## RNA-seq\n"
        "- mRNA를 cDNA로 바꿔 시퀀싱. 넓은 동적 범위, probe 없음.\n\n"
        "## Microarray\n"
        "- 미리 설계된 probe에 mRNA가 붙으면 형광 강도 측정. 저비용, 오래됐음.",
    )
    E(
        "생물학 (일반)",
        "Methylation / Beta Value",
        "DNA 메틸화 측정",
        "## methylation\n"
        "- CpG 부위에 붙는 메틸기(-CH3)의 정도. 유전자 발현을 끄는 대표 조절.\n\n"
        "## beta value\n"
        "- `methylated / (methylated + unmethylated + offset)`, 0~1.\n"
        "- 0이면 완전 비메틸, 1이면 완전 메틸.",
    )
    E(
        "생물학 (일반)",
        "Mutation Types",
        "돌연변이 종류",
        "## missense\n"
        "- 아미노산 하나가 바뀜 (예: V600E).\n\n"
        "## nonsense\n"
        "- 중간에 stop codon 생성 → 단백질 잘림.\n\n"
        "## frameshift\n"
        "- 삽입/결실로 reading frame이 밀려 이후가 엉망.\n\n"
        "## splice site\n"
        "- intron/exon 경계 파괴, splicing 이상.",
    )
    E(
        "생물학 (일반)",
        "Log2 Fold Change",
        "발현 변화의 로그 비율",
        "## 정의\n"
        "- `log2(A / B)`. A가 B의 2배이면 `+1`, 절반이면 `−1`.\n\n"
        "## 왜 로그인가\n"
        "- 발현은 수십 배 차이가 흔하므로 로그로 펼쳐야 시각적 비교가 쉽고 통계도 안정.",
    )
    E(
        "생물학 (일반)",
        "Volcano Plot",
        "변화 크기와 유의성 동시 보기",
        "## 축\n"
        "- x: log2 fold change, y: `-log10(p-value)`.\n\n"
        "## 해석\n"
        "- 우상단·좌상단 점: 크게 변했고 통계적으로 유의함 → 관심 유전자.\n"
        "- 색으로 FDR 통과 여부를 구분하기도 합니다.",
    )
    E(
        "생물학 (일반)",
        "Heatmap / Clustering",
        "히트맵과 군집화",
        "## heatmap\n"
        "- 행=유전자, 열=샘플, 색=발현. 패턴을 눈으로 확인.\n\n"
        "## clustering\n"
        "- hierarchical(dendrogram) 또는 k-means로 비슷한 유전자/샘플을 묶음.\n\n"
        "## 주의\n"
        "- 정규화(row z-score 등) 방식에 따라 보이는 패턴이 달라집니다.",
    )
    E(
        "생물학 (일반)",
        "Batch Effect / ComBat",
        "기술 잡음과 보정",
        "## batch effect\n"
        "- 실험 일자·장비·시약 차이로 생기는 측정 편차.\n\n"
        "## ComBat\n"
        "- empirical Bayes 기반 보정 기법. batch 변수만 알면 적용 가능.\n\n"
        "## 주의\n"
        "- batch와 관심 변수가 너무 얽혀 있으면 보정이 관심 신호까지 깎을 수 있습니다.",
    )
    E(
        "모델 평가",
        "SHAP",
        "모델 해석 값",
        "## 정의\n"
        "- Shapley value(게임이론)를 ML 예측에 적용한 해석 도구.\n\n"
        "## 해석\n"
        "- 각 feature가 예측을 base에서 위/아래로 얼마나 밀었는지.\n"
        "- 모든 feature의 SHAP 합 + base = 실제 예측(덧셈 가능).\n\n"
        "## 주의\n"
        "- 상관 feature끼리는 몫이 임의로 나뉘므로 개별 순위 신뢰는 과하게 하지 마세요.",
    )

    # --- 데이터 출처 ---
    E(
        "데이터 소스",
        "TCGA / GEO / GDC",
        "공개 유전체 데이터 저장소",
        "## TCGA\n"
        "- The Cancer Genome Atlas. NCI/NHGRI 공동. 33개 암종의 multi-omics 공개.\n\n"
        "## GEO\n"
        "- NCBI Gene Expression Omnibus. 각 연구팀이 업로드한 공개 전사체 데이터 저장소.\n\n"
        "## GDC\n"
        "- Genomic Data Commons. TCGA 후속 공식 배포처. 표준화된 API로 raw/processed 데이터 접근.",
    )
    E(
        "데이터 소스",
        "PubChem / ChEMBL / DrugBank",
        "화합물/약물 데이터베이스",
        "## PubChem\n"
        "- NIH 산하. 화합물 구조, 생물활성 데이터 광범위.\n\n"
        "## ChEMBL\n"
        "- EMBL-EBI. 약물성 화합물 위주, 표적과 활성값(pChEMBL) 풍부.\n\n"
        "## DrugBank\n"
        "- 승인약/실험약 정보. 임상/약동학 중심.",
    )
    E(
        "드럭 디스커버리",
        "pChEMBL Value",
        "활성값 로그 스케일",
        "## 정의\n"
        "- `pChEMBL = -log10(IC50 또는 Ki [mol/L])`.\n\n"
        "## 대응표\n"
        "- `pChEMBL 5` → IC50 10 μM\n"
        "- `pChEMBL 6` → 1 μM\n"
        "- `pChEMBL 7` → 100 nM\n"
        "- `pChEMBL 8` → 10 nM\n"
        "- `pChEMBL 9` → 1 nM\n\n"
        "## 실무\n"
        "- 신약 후보로 보통 **pChEMBL ≥ 7** 수준을 기대합니다.",
    )
    E(
        "드럭 디스커버리",
        "SMILES",
        "화합물 문자열 표현",
        "## 정의\n"
        "- Simplified Molecular Input Line Entry System.\n"
        "- 원자·결합을 ASCII 문자열로 표현해 화합물 구조를 한 줄로 저장.\n\n"
        "## 예시\n"
        "- aspirin: `CC(=O)OC1=CC=CC=C1C(=O)O`.\n\n"
        "## 활용\n"
        "- 머신러닝 입력, 2D 구조 렌더링, 검색 키.",
    )
    E(
        "드럭 디스커버리",
        "PDB / Structure-Based Drug Design",
        "단백질 구조 기반 설계",
        "## PDB\n"
        "- Protein Data Bank. 단백질 3D 구조 저장소(X-ray, cryo-EM, NMR).\n\n"
        "## structure-based drug design\n"
        "- 표적의 결합 포켓 구조를 알면 **도킹·자유에너지 계산**으로 후보를 스크리닝.\n\n"
        "## 한계\n"
        "- 구조가 있어도 동적 배치, 용매, 알로스테릭 효과는 포착하기 어렵습니다.",
    )
    E(
        "드럭 디스커버리",
        "Druggability",
        "표적이 약으로 공략 가능한가",
        "## 정의\n"
        "- 작은 분자(<500 Da)가 표적에 충분히 강하게 붙어 생물학적 효과를 낼 수 있는 정도.\n\n"
        "## 판단 요소\n"
        "- 결합 포켓의 깊이/모양, 표면 소수성, 기존 ligand 유무.\n\n"
        "## 점수\n"
        "- fpocket, SiteMap 같은 도구가 수치화합니다.",
    )
    E(
        "드럭 디스커버리",
        "ADC (Antibody-Drug Conjugate)",
        "항체-약물 접합체",
        "## 정의\n"
        "- 표적 세포 표면 단백질에 특이적인 항체에 **세포독성 payload**를 화학적으로 연결한 약.\n\n"
        "## 장점\n"
        "- 종양 세포에 선택적으로 독성을 전달해 off-target을 줄임.\n\n"
        "## 대표약\n"
        "- trastuzumab-emtansine (HER2+ 유방암), enfortumab vedotin.",
    )
    E(
        "드럭 디스커버리",
        "Kinase / Protease / Phosphatase",
        "효소 세 가족",
        "## kinase\n"
        "- 단백질/분자에 **인산기를 붙이는** 효소. BRAF, MAPK 등. 신약 표적으로 가장 인기.\n\n"
        "## protease\n"
        "- 단백질을 **자르는** 효소. HIV protease inhibitor가 대표 성공사례.\n\n"
        "## phosphatase\n"
        "- 인산기를 **떼는** 효소. kinase의 반대 역할.",
    )
    E(
        "드럭 디스커버리",
        "Target Classification (Novel / Emerging / Validated)",
        "표적 성숙도 분류",
        "## novel_target\n"
        "- 신규 발견, 공개 문헌/임상 근거 거의 없음. 발견 가치는 크지만 실패 위험도 큼.\n\n"
        "## emerging\n"
        "- 초기 증거(일부 논문, 초기 임상) 축적. 우리 데이터로 재확인이 의미 있음.\n\n"
        "## validated\n"
        "- 승인약 또는 후기 임상이 존재. 재현 중심.",
    )
    E(
        "드럭 디스커버리",
        "Literature (PubMed)",
        "문헌 근거 읽는 법",
        "## 체크할 점\n"
        "- **recency**: 최근 5년 내 논문 비중이 높을수록 활발한 주제.\n"
        "- **sample size**: n=1000 이상의 대규모 연구에 가중치.\n"
        "- **relevance**: 갑상선 / THCA / PTC / FTC 등 키워드 매칭.\n\n"
        "## 주의\n"
        "- 제목/초록만으로 판단하면 overclaim 논문에 속기 쉽습니다. 가급적 원문 method와 결론을 확인.",
    )

    # --- 양자 ML ---
    E(
        "양자 ML",
        "VQC (Variational Quantum Classifier)",
        "파라미터화 양자 분류기",
        "## 개념\n"
        "- 데이터를 양자 상태로 인코딩 → 학습 가능한 회로 U(θ) 적용 → 측정 → 분류.\n"
        "- 고전 optimizer가 θ를 업데이트.\n\n"
        "## 언제 도움 될까\n"
        "- 표현하기 어려운 커널을 양자가 더 적은 자원으로 표현 가능할 때(이론).\n\n"
        "## 언제 안 되나\n"
        "- noisy qubit으로는 얕은 회로만 가능. 대부분의 tabular 데이터는 고전이 더 잘 합니다.",
    )
    E(
        "양자 ML",
        "QSVM (Quantum SVM)",
        "양자 커널 기반 SVM",
        "## 개념\n"
        "- 커널 `K(x,x')`를 양자 회로에서 계산한 뒤 고전 SVM에 공급.\n\n"
        "## 기대\n"
        "- 고차원 Hilbert space에서 classical로 계산 비싼 커널을 저렴하게.\n\n"
        "## 현실\n"
        "- 현재 하드웨어에서는 샘플 수가 커지면 측정 반복 횟수가 폭증해 실용성이 제한됩니다.",
    )
    E(
        "양자 ML",
        "QAOA",
        "조합 최적화용 양자 알고리즘",
        "## 정의\n"
        "- Quantum Approximate Optimization Algorithm.\n"
        "- MaxCut 같은 QUBO 문제를 얕은 회로로 근사.\n\n"
        "## 쓰임\n"
        "- feature selection(유전자 조합 최적화) 등 조합 문제에 이론적 후보.\n\n"
        "## 한계\n"
        "- p(회로 깊이)가 작으면 고전 heuristic과 경쟁하기 어렵고, 크면 noise에 무너집니다.",
    )
    E(
        "양자 ML",
        "QNN (Quantum Neural Network)",
        "양자 회로로 만든 신경망",
        "## 개념\n"
        "- 파라미터화된 양자 회로(PQC)를 신경망 레이어처럼 쌓고 gradient descent로 학습.\n\n"
        "## 특징\n"
        "- trainability 문제(barren plateau)가 심각합니다.\n\n"
        "## 실용성\n"
        "- 현재는 연구 수준. 작은 toy 문제에서만 classical baseline과 비교 가능.",
    )
    E(
        "양자 ML",
        "qPCA (Quantum PCA)",
        "양자 버전 주성분 분석",
        "## 아이디어\n"
        "- 양자 상태의 밀도행렬에서 주요 고유값을 지수적으로 빠르게 추정(이론).\n\n"
        "## 조건\n"
        "- 데이터가 이미 양자 상태로 준비되어 있어야 합니다(QRAM 가정).\n\n"
        "## 현실\n"
        "- QRAM이 없으므로 실 데이터에서 빠른 속도를 보장하기 어렵습니다.",
    )
    E(
        "양자 ML",
        "Grover",
        "비정렬 DB 검색 √N 알고리즘",
        "## 개념\n"
        "- N개 중 표시된 항목을 찾는 데 `O(√N)`.\n\n"
        "## ML 접목\n"
        "- amplitude amplification으로 rare event 탐색, sampling 가속.\n\n"
        "## 한계\n"
        "- oracle 구현 비용을 고려하면 실 데이터 전체에서 속도 이득을 보기 어렵습니다.",
    )
    E(
        "양자 ML",
        "QBoost",
        "양자 어닐러 기반 앙상블",
        "## 개념\n"
        "- 약분류기 선택을 QUBO로 정식화해 D-Wave 같은 annealer에 맡김.\n\n"
        "## 기대\n"
        "- 가중치 조합 공간이 클 때 대안적 탐색.\n\n"
        "## 현실\n"
        "- 문제 크기/임베딩 제약이 크고, 고전 XGBoost를 이기는 사례는 제한적.",
    )
    E(
        "양자 ML",
        "QKernel",
        "양자 커널 함수",
        "## 개념\n"
        "- 두 데이터 포인트를 양자 상태 `|φ(x)⟩`, `|φ(x')⟩`로 인코딩 후 `|⟨φ(x)|φ(x')⟩|^2`를 커널로.\n\n"
        "## 용도\n"
        "- SVM, kernel regression 등 classical kernel method의 드롭인 대체.\n\n"
        "## 주의\n"
        "- 커널 평가당 여러 번의 회로 실행이 필요합니다.",
    )
    E(
        "양자 ML",
        "QUBO",
        "이진 제곱 최적화 형태",
        "## 정의\n"
        "- Quadratic Unconstrained Binary Optimization. `minimize x^T Q x`, `x∈{0,1}^n`.\n\n"
        "## 쓰임\n"
        "- QAOA(게이트형)와 D-Wave(어닐러)가 공통으로 받는 표준 입력 형태.\n\n"
        "## ML 접목\n"
        "- feature selection, clustering, portfolio를 QUBO로 재표현해 양자 기기에 전달.",
    )
    E(
        "양자 ML",
        "D-Wave Annealing",
        "양자 어닐링 머신",
        "## 개념\n"
        "- 이징/QUBO 문제의 에너지 최저 상태를 quantum tunneling으로 찾는 방식.\n\n"
        "## 장점\n"
        "- 수천 qubit 급 문제 입력 가능.\n\n"
        "## 한계\n"
        "- 완전 일반 목적이 아니고, 문제→그래프 임베딩 비용이 크며, 고전 heuristic을 일관되게 이기지 못하는 사례가 많습니다.",
    )

    return entries


# ---------------------------------------------------------------------------
# HTML building
# ---------------------------------------------------------------------------


def build_details_block(section_id: str, entry: dict) -> str:
    body_html = markdown_to_html(entry["body_md"])
    title = html.escape(entry.get("title", "").strip())
    summary_text = html.escape(entry.get("summary", "").strip() or title)
    return (
        f'<details class="edu-toggle" data-edu-id="{section_id}">'
        f'<!-- edu:{section_id} -->'
        f'<summary>{summary_text}</summary>'
        f'<div class="edu-body">'
        f'<h3>{title}</h3>'
        f'{body_html}'
        f'</div>'
        f'</details>'
    )


# ---------------------------------------------------------------------------
# Injection logic
# ---------------------------------------------------------------------------


_SECTION_H2_RE = re.compile(
    r'(<h2\s+class="section-title"\s+id="sec-(?P<sid>[^"]+)"[^>]*>.*?</h2>)',
    re.IGNORECASE | re.DOTALL,
)

# Fallback pattern: <section id="xxx"> ... <h2 ...>label</h2>
# Covers two variants observed in the repo:
#   (a) index.html: <section id="x"> ... <h2 class="section-title">...</h2>
#   (b) 15_drug_discovery.html: <section ... id="x"> ... <h2>...</h2>
# We bound the span with a lazy match to the first closing </h2>.
_SECTION_PARENT_RE = re.compile(
    r'(<section\b[^>]*\bid="(?P<sid>[^"]+)"[^>]*>(?:(?!</section>).){0,2000}?<h2\b[^>]*>[^<]*</h2>)',
    re.IGNORECASE | re.DOTALL,
)


_DETAILS_BLOCK_RE = re.compile(
    r'<details class="edu-toggle" data-edu-id="(?P<sid>[^"]+)">.*?</details>\s*',
    re.DOTALL,
)


def cleanup_orphan_toggles(html_text: str, page_name: str, content: dict) -> tuple[str, int]:
    """Remove toggle blocks injected by a prior pass that are not in ``content``
    AND whose sid refers to a chrome/utility section (no matching h2 id on the
    page). This protects against accidental generic-fills from early runs.
    """
    removed = 0

    def _remove(m: re.Match) -> str:
        nonlocal removed
        sid = m.group("sid")
        key = f"{page_name}:{sid}"
        if key in content:
            return m.group(0)
        # If there's a matching <h2 ... id="sec-<sid>"> then the block is legit
        # (generic fallback for a curated-but-missing key, which we allow).
        if f'id="sec-{sid}"' in html_text:
            return m.group(0)
        removed += 1
        return ""

    new_html = _DETAILS_BLOCK_RE.sub(_remove, html_text)
    return new_html, removed


def inject_into_page(html_text: str, page_name: str, content: dict, generic: dict) -> tuple[str, int]:
    """Insert a <details> block after every matching h2 on the page.

    Two strategies, run in order:
      1. h2 that has ``id="sec-<sid>"`` directly.
      2. Parent ``<section id="<sid>">`` that contains a ``<h2 class="section-title">``.
    Idempotency is guaranteed by the ``<!-- edu:<sid> -->`` marker comment.

    Returns (new_html, toggles_inserted).
    """
    inserted = 0

    def _make_block(sid: str) -> str:
        key = f"{page_name}:{sid}"
        entry = content.get(key)
        if entry is None:
            entry = dict(generic)
            entry["title"] = f"{sid} — 섹션 안내"
        return build_details_block(sid, entry)

    def replace_direct(match: re.Match) -> str:
        nonlocal inserted
        sid = match.group("sid")
        full = match.group(1)
        if f"<!-- edu:{sid} -->" in html_text:
            return full
        inserted += 1
        return full + "\n" + _make_block(sid)

    new_html = _SECTION_H2_RE.sub(replace_direct, html_text)

    # Fallback pass on the updated html.
    # Skip framework/chrome ids that aren't user-facing sections.
    skip_ids = {
        "command-palette",
        "figure-modal",
        "td-drawer",
        "td-drawer-backdrop",
        "main-content",
    }

    def replace_parent(match: re.Match) -> str:
        nonlocal inserted, new_html
        sid = match.group("sid")
        full = match.group(1)
        if sid in skip_ids:
            return full
        if f"<!-- edu:{sid} -->" in new_html:
            return full
        # Skip if the section already contains an h2 with id="sec-..." that we handle.
        if 'id="sec-' in full:
            return full
        # Only inject when we have a curated entry — do not generic-fill for chrome sections.
        key = f"{page_name}:{sid}"
        if key not in content:
            return full
        inserted += 1
        return full + "\n" + _make_block(sid)

    new_html = _SECTION_PARENT_RE.sub(replace_parent, new_html)
    return new_html, inserted


def inject_css_link(html_text: str, href: str) -> str:
    """Inject stylesheet link once, right before </head>."""
    if "<!-- edu-css -->" in html_text:
        return html_text
    link = f'<!-- edu-css --><link rel="stylesheet" href="{href}">'
    return html_text.replace("</head>", link + "\n</head>", 1)


def inject_glossary_navlink(html_text: str, href: str) -> str:
    """Add a '용어' link to the topnav and drawer. Idempotent via marker."""
    if "<!-- edu-glossary-nav -->" in html_text:
        return html_text
    # Insert right before closing </div> of .nav-links (first occurrence).
    # Pattern: find the block that opens with class="nav-links" and ends with </div>.
    nav_pattern = re.compile(r'(<div class="nav-links"[^>]*>)(.*?)(</div>)', re.DOTALL)
    first_done = {"n": 0}

    def _replace(m: re.Match) -> str:
        if first_done["n"] > 0:
            return m.group(0)
        first_done["n"] += 1
        inside = m.group(2)
        link_html = (
            '<!-- edu-glossary-nav --><a class="nav-link" href="'
            + href
            + '" role="menuitem">용어</a>'
        )
        return m.group(1) + inside + link_html + m.group(3)

    new_html = nav_pattern.sub(_replace, html_text, count=1)
    # Also add to drawer (<aside id="td-drawer"> contains bare <a> links).
    drawer_re = re.compile(r'(<aside id="td-drawer"[^>]*>.*?)(</aside>)', re.DOTALL)

    def _drawer(m: re.Match) -> str:
        inner = m.group(1)
        if 'href="' + href + '"' in inner:
            return m.group(0)
        drawer_link = f'<a class="" href="{href}">용어</a>'
        return inner + drawer_link + m.group(2)

    new_html = drawer_re.sub(_drawer, new_html, count=1)
    return new_html


# ---------------------------------------------------------------------------
# Glossary page builder
# ---------------------------------------------------------------------------


def build_glossary_page(entries: list[dict]) -> str:
    # Group by category preserving insertion order.
    cats: dict[str, list[dict]] = {}
    for e in entries:
        cats.setdefault(e["category"], []).append(e)

    toc = "".join(
        f'<a href="#cat-{re.sub(r"[^a-z0-9]+", "-", cat.lower()).strip("-")}">{html.escape(cat)}</a>'
        for cat in cats
    )
    sections_html: list[str] = []
    for cat, items in cats.items():
        cat_id = "cat-" + re.sub(r"[^a-z0-9]+", "-", cat.lower()).strip("-")
        cards = []
        for e in items:
            block = build_details_block(e["id"], {
                "title": e["term"],
                "summary": e["summary"],
                "body_md": e["body_md"],
            })
            card = (
                f'<div class="glossary-card" id="{e["id"]}">'
                f'<div class="glossary-term">{html.escape(e["term"])}</div>'
                f'{block}'
                f'</div>'
            )
            cards.append(card)
        sections_html.append(
            f'<section id="{cat_id}" class="panel glass" aria-labelledby="hd-{cat_id}">'
            f'<h2 class="section-title" id="hd-{cat_id}">{html.escape(cat)}</h2>'
            f'<div class="glossary-grid">{"".join(cards)}</div>'
            f'</section>'
        )

    nav_links_html = (
        '<a class="nav-link" href="../index.html">홈</a>'
        '<a class="nav-link" href="01_overview.html">개요</a>'
        '<a class="nav-link" href="02_datasets.html">데이터셋</a>'
        '<a class="nav-link" href="03_sample_master.html">샘플</a>'
        '<a class="nav-link" href="04_gene_panels.html">패널</a>'
        '<a class="nav-link" href="05_eda.html">EDA</a>'
        '<a class="nav-link" href="06_scores.html">점수</a>'
        '<a class="nav-link" href="07_ml_baseline.html">ML</a>'
        '<a class="nav-link" href="08_panel_comparison.html">패널 비교</a>'
        '<a class="nav-link" href="09_shap.html">SHAP</a>'
        '<a class="nav-link" href="10_gene_explorer.html">유전자</a>'
        '<a class="nav-link" href="11_cohort_compare.html">코호트</a>'
        '<a class="nav-link" href="12_business.html">비즈니스</a>'
        '<a class="nav-link" href="13_reports.html">리포트</a>'
        '<a class="nav-link" href="14_caveats.html">주의점</a>'
        '<a class="nav-link" href="15_drug_discovery.html">드럭</a>'
        '<!-- edu-glossary-nav --><a class="nav-link active" href="99_glossary.html">용어</a>'
    )

    return f"""<!doctype html>
<html lang="ko">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="dark light">
<meta name="theme-color" content="#07132b">
<title>용어집 | THCA Multi-Omics Dashboard</title>
<meta name="description" content="생물정보/ML/양자ML 용어 한 장 가이드">
<link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
<link rel="stylesheet" href="../assets/css/main.css">
<!-- edu-css --><link rel="stylesheet" href="../assets/css/edu-toggle.css">
<style>
  .glossary-grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:14px;margin-top:12px}}
  .glossary-card{{padding:12px 14px;border-radius:14px;border:1px solid rgba(148,163,184,.18);background:rgba(11,26,58,.35)}}
  .glossary-term{{font-weight:700;font-size:1.05rem;color:#e7fff8;letter-spacing:-.01em;margin-bottom:4px}}
  .glossary-card .edu-toggle{{margin:6px 0 0}}
</style>
</head>
<body class="dark">
<nav class="topnav" aria-label="주요 메뉴">
  <div class="topnav-inner">
    <a class="brand" href="../index.html" aria-label="홈으로"><span class="brand-dot" aria-hidden="true"></span><span>THYROID DASH</span></a>
    <div class="nav-links" role="menubar">{nav_links_html}</div>
    <div class="nav-actions">
      <button class="btn sm ghost" type="button" onclick="toggleTheme()" aria-label="다크/라이트 모드 전환">다크/라이트</button>
    </div>
  </div>
</nav>
<aside id="td-drawer" class="drawer" aria-hidden="true" aria-label="모바일 메뉴">
  <div class="drawer-header"><strong>메뉴</strong></div>
  <a href="../index.html">홈</a><a href="01_overview.html">개요</a><a href="02_datasets.html">데이터셋</a><a href="03_sample_master.html">샘플</a><a href="04_gene_panels.html">패널</a><a href="05_eda.html">EDA</a><a href="06_scores.html">점수</a><a href="07_ml_baseline.html">ML</a><a href="08_panel_comparison.html">패널 비교</a><a href="09_shap.html">SHAP</a><a href="10_gene_explorer.html">유전자</a><a href="11_cohort_compare.html">코호트</a><a href="12_business.html">비즈니스</a><a href="13_reports.html">리포트</a><a href="14_caveats.html">주의점</a><a href="15_drug_discovery.html">드럭</a><a class="active" href="99_glossary.html">용어</a>
</aside>
<div class="layout">
  <aside class="sidebar glass" aria-label="용어집 목차">
    <h3>카테고리</h3>
    {toc}
  </aside>
  <main class="content" id="main-content">
    <section class="hero" aria-labelledby="glossary-title">
      <div class="subtle mono">99 / GLOSSARY</div>
      <h1 id="glossary-title">용어집</h1>
      <p class="lede">생물정보, 머신러닝, 양자 ML까지 이 대시보드에 등장하는 모든 기술 용어를 한 페이지에 모았습니다. 항목을 펼치면 200자 안팎의 초보자용 설명이 나옵니다.</p>
    </section>
    {''.join(sections_html)}
    <footer class="footer glass" role="contentinfo">
      <div><span class="mono">THCA Dashboard</span> · 용어집 · 새 항목이 필요하면 analysis team에 요청하세요.</div>
    </footer>
  </main>
</div>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# CSS
# ---------------------------------------------------------------------------


CSS_BODY = """.edu-toggle{margin:12px 0 8px;border-radius:12px;border:1px solid rgba(148,163,184,.18);background:rgba(11,26,58,.35);overflow:hidden}
.edu-toggle>summary{cursor:pointer;padding:10px 14px;font-size:.85rem;color:var(--teal-300, #5eead4);list-style:none;position:relative}
.edu-toggle>summary::before{content:"\\1F4A1 \\0020 \\CC98\\C74C \\0020 \\BCF4\\B294 \\0020 \\ACBD\\C6B0 \\0020 \\2014 \\0020 ";font-weight:500}
.edu-toggle[open]>summary::before{content:"\\1F4D6 \\0020 \\BC30\\ACBD \\0020 \\C124\\BA85 \\0020 \\2014 \\0020 "}
.edu-toggle>summary::-webkit-details-marker{display:none}
.edu-toggle>summary::after{content:"\\25B8";float:right;transition:transform .2s}
.edu-toggle[open]>summary::after{transform:rotate(90deg)}
.edu-toggle .edu-body{padding:4px 16px 14px;border-top:1px dashed rgba(148,163,184,.18);color:var(--ink, #e2e8f0);font-size:.9rem;line-height:1.6}
.edu-toggle .edu-body h3{font-size:.95rem;margin:10px 0 4px;color:var(--teal-300, #5eead4)}
.edu-toggle .edu-body h4{font-size:.85rem;margin:8px 0 3px;color:var(--teal-300, #5eead4);text-transform:uppercase;letter-spacing:.06em}
.edu-toggle .edu-body h5{font-size:.8rem;margin:6px 0 2px;color:var(--teal-300, #5eead4)}
.edu-toggle .edu-body ul,.edu-toggle .edu-body ol{padding-left:20px;margin:4px 0 8px}
.edu-toggle .edu-body li{margin:2px 0}
.edu-toggle .edu-body p{margin:6px 0}
.edu-toggle .edu-body code{padding:1px 6px;border-radius:4px;background:rgba(20,184,166,.12);color:var(--teal-300, #5eead4);font-size:.85rem}
.edu-toggle .edu-body strong{color:#e7fff8}
"""


# Use literal strings for CSS content — avoid hex escape subtleties by using actual unicode.
CSS_BODY = (
    ".edu-toggle{margin:12px 0 8px;border-radius:12px;border:1px solid rgba(148,163,184,.18);background:rgba(11,26,58,.35);overflow:hidden}\n"
    ".edu-toggle>summary{cursor:pointer;padding:10px 14px;font-size:.85rem;color:var(--teal-300, #5eead4);list-style:none;position:relative}\n"
    ".edu-toggle>summary::before{content:\"\U0001F4A1 \\00A0 \\CC98\\C74C \\00A0 \\BCF4\\B294 \\00A0 \\ACBD\\C6B0 \\00A0 \\2014 \\00A0 \";font-weight:500}\n"
    ".edu-toggle[open]>summary::before{content:\"\U0001F4D6 \\00A0 \\BC30\\ACBD \\00A0 \\C124\\BA85 \\00A0 \\2014 \\00A0 \"}\n"
    ".edu-toggle>summary::-webkit-details-marker{display:none}\n"
    ".edu-toggle>summary::after{content:\"\\25B8\";float:right;transition:transform .2s}\n"
    ".edu-toggle[open]>summary::after{transform:rotate(90deg)}\n"
    ".edu-toggle .edu-body{padding:4px 16px 14px;border-top:1px dashed rgba(148,163,184,.18);color:var(--ink, #e2e8f0);font-size:.9rem;line-height:1.6}\n"
    ".edu-toggle .edu-body h3{font-size:.95rem;margin:10px 0 4px;color:var(--teal-300, #5eead4)}\n"
    ".edu-toggle .edu-body h4{font-size:.85rem;margin:8px 0 3px;color:var(--teal-300, #5eead4);text-transform:uppercase;letter-spacing:.06em}\n"
    ".edu-toggle .edu-body h5{font-size:.8rem;margin:6px 0 2px;color:var(--teal-300, #5eead4)}\n"
    ".edu-toggle .edu-body ul,.edu-toggle .edu-body ol{padding-left:20px;margin:4px 0 8px}\n"
    ".edu-toggle .edu-body li{margin:2px 0}\n"
    ".edu-toggle .edu-body p{margin:6px 0}\n"
    ".edu-toggle .edu-body code{padding:1px 6px;border-radius:4px;background:rgba(20,184,166,.12);color:var(--teal-300, #5eead4);font-size:.85rem}\n"
    ".edu-toggle .edu-body strong{color:#e7fff8}\n"
)


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------


def main() -> dict:
    stats: dict = {"pages": 0, "toggles": 0, "per_page": {}, "glossary_terms": 0, "json_size": 0}

    # 1) Ensure content JSON exists.
    CONTENT_PATH.parent.mkdir(parents=True, exist_ok=True)
    content, generic = _build_default_content()
    if CONTENT_PATH.exists():
        try:
            existing = json.loads(CONTENT_PATH.read_text(encoding="utf-8"))
            if isinstance(existing, dict) and existing.get("sections"):
                # Merge: existing file wins (operator edits), missing keys filled in.
                existing_sec = existing.get("sections", {})
                for k, v in content.items():
                    existing_sec.setdefault(k, v)
                payload = {"sections": existing_sec, "generic": existing.get("generic", generic)}
                CONTENT_PATH.write_text(
                    json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
                )
                content = existing_sec
                generic = payload["generic"]
            else:
                payload = {"sections": content, "generic": generic}
                CONTENT_PATH.write_text(
                    json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
                )
        except Exception:
            payload = {"sections": content, "generic": generic}
            CONTENT_PATH.write_text(
                json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
            )
    else:
        payload = {"sections": content, "generic": generic}
        CONTENT_PATH.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
        )
    stats["json_size"] = CONTENT_PATH.stat().st_size

    # 2) Write CSS.
    CSS_PATH.parent.mkdir(parents=True, exist_ok=True)
    CSS_PATH.write_text(CSS_BODY, encoding="utf-8")

    # 3) Build glossary page.
    glossary_entries = _build_glossary()
    GLOSSARY_PATH.write_text(build_glossary_page(glossary_entries), encoding="utf-8")
    stats["glossary_terms"] = len(glossary_entries)

    # 4) Walk all pages.
    targets: list[Path] = [INDEX_PATH]
    targets.extend(sorted(p for p in PAGES_DIR.glob("*.html") if p.name != "99_glossary.html"))

    for page in targets:
        original = page.read_text(encoding="utf-8")
        page_name = page.name
        # CSS path is relative depending on directory depth.
        css_href = (
            "assets/css/edu-toggle.css" if page.parent.name == "html" else "../assets/css/edu-toggle.css"
        )
        glossary_href = (
            "pages/99_glossary.html" if page.parent.name == "html" else "99_glossary.html"
        )
        text = original
        text = inject_css_link(text, css_href)
        text = inject_glossary_navlink(text, glossary_href)
        text, _removed = cleanup_orphan_toggles(text, page_name, content)
        text, inserted = inject_into_page(text, page_name, content, generic)
        if text != original:
            page.write_text(text, encoding="utf-8")
        stats["per_page"][page_name] = inserted
        stats["toggles"] += inserted
        stats["pages"] += 1

    # Regenerate glossary too (in case content changed).
    GLOSSARY_PATH.write_text(build_glossary_page(glossary_entries), encoding="utf-8")
    return stats


if __name__ == "__main__":
    s = main()
    import pprint

    pprint.pprint(s)
