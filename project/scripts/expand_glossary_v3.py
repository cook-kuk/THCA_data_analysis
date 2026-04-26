"""Extend the THCA dashboard glossary (99_glossary.html) with v3 additions.

Adds a new section titled "v3 추가 용어 (v3 additions)" with beginner-friendly
<details> entries for terms used by the v3 pages (honesty audit, panel size
curve, three-class, survival, bethesda, multimodal).

Idempotent:
  - Marker comment <!-- v3-glossary-additions --> wraps the new section.
  - If present, re-running is a no-op.
  - Skips any term whose id (term-<slug>) already exists in the file.

Do NOT modify any Python pipeline script or serve_secure.py.
"""
from __future__ import annotations

import html
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GLOSSARY_PATH = ROOT / "reports" / "html" / "pages" / "99_glossary.html"
MARKER = "<!-- v3-glossary-additions -->"


# ---------------------------------------------------------------------------
# Tiny markdown -> HTML (mirrors inject_education_toggles.markdown_to_html
# so the visual style matches other glossary bodies exactly).
# ---------------------------------------------------------------------------

_INLINE_CODE = re.compile(r"`([^`]+)`")
_BOLD = re.compile(r"\*\*([^*]+)\*\*")
_ITALIC = re.compile(r"(?<!\*)\*([^*]+)\*(?!\*)")


def _render_inline(line: str) -> str:
    out = html.escape(line, quote=False)
    out = _INLINE_CODE.sub(lambda m: f"<code>{m.group(1)}</code>", out)
    out = _BOLD.sub(lambda m: f"<strong>{m.group(1)}</strong>", out)
    out = _ITALIC.sub(lambda m: f"<em>{m.group(1)}</em>", out)
    return out


def markdown_to_html(text: str) -> str:
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
        m = re.match(r"^(#{1,4})\s+(.*)$", stripped)
        if m:
            level = len(m.group(1))
            content = _render_inline(m.group(2).strip())
            tag = f"h{min(level + 2, 5)}"
            out.append(f"<{tag}>{content}</{tag}>")
            i += 1
            continue
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


def _slug(term: str) -> str:
    """Mirror the slug convention used by existing entries (term-<slug>)."""
    s = term.lower()
    s = re.sub(r"[^a-z0-9]+", "-", s)
    s = re.sub(r"-+", "-", s).strip("-")
    return f"term-{s}"


# ---------------------------------------------------------------------------
# New v3 glossary entries.
# Each entry: (term display, short summary, markdown body 60-150 words).
# Terms already in glossary (Brier, Hazard ratio) are intentionally skipped.
# ---------------------------------------------------------------------------

ENTRIES: list[tuple[str, str, str]] = [
    (
        "Leakage curve (누출 곡선)",
        "train leak이 성능을 얼마나 부풀리는가",
        "## 무엇을 측정하나\n"
        "- 학습 데이터와 평가 데이터 사이의 **겹침(sample leak)**이 AUC를 얼마나 끌어올리는지 시뮬레이션한 곡선입니다.\n"
        "- x축은 인위적으로 넣은 leak 비율(0→100%), y축은 그에 따른 내부 AUC입니다.\n\n"
        "## 읽는 법\n"
        "- **leak=0%**에서의 AUC가 진짜 성능에 가깝습니다.\n"
        "- leak을 조금만 넣어도 AUC가 급격히 오르면, 원래 실험도 *우연한* 겹침에 취약하다는 신호입니다.\n"
        "- 곡선이 완만하면 평가 설계가 leak에 robust.\n\n"
        "## 주의\n"
        "- leakage curve는 **시뮬레이션**이며 실제 leak을 탐지하지는 않습니다. nested CV, group split과 함께 봐야 합니다.",
    ),
    (
        "Permutation null (순열 null 분포)",
        "라벨을 섞어서 얻은 empirical p",
        "## 정의\n"
        "- 라벨을 무작위로 <strong>섞은(permute)</strong> 뒤 같은 분류기를 학습/평가해 얻는 **null AUC 분포**입니다.\n\n"
        "## 왜 쓰나\n"
        "- 이론적 p-value 가정(정규성, 독립성)이 흔들리는 고차원 `n≪p` 상황에서 **분포-free**로 유의성을 판정.\n"
        "- 실제 AUC가 null 분포의 상위 1%에 들어가면 `p<0.01`이라 말할 수 있습니다.\n\n"
        "## 읽는 법\n"
        "- 히스토그램 중심은 0.5 근처여야 정상. 중심이 0.6 이상이면 파이프라인 내부에 leak 의심.\n\n"
        "## 계산 비용\n"
        "- 보통 1,000~10,000회 permutation이 필요합니다. GPU 없이도 tabular면 충분합니다.",
    ),
    (
        "LODO (Leave-One-Dataset-Out)",
        "한 코호트 통째로 빼고 평가",
        "## 정의\n"
        "- **Leave-One-Dataset-Out**: 여러 코호트 중 하나를 통째로 테스트 셋으로 빼고, 나머지로만 학습합니다.\n\n"
        "## 왜 중요한가\n"
        "- 같은 코호트 내 random CV는 **batch-specific 신호**까지 공유해 AUC가 과대 평가됩니다.\n"
        "- LODO는 실제 임상 배포 상황(새 병원, 새 플랫폼)에 가장 가깝습니다.\n\n"
        "## 해석\n"
        "- LODO AUC < 내부 CV AUC면 <em>generalization gap</em>이 존재한다는 증거.\n"
        "- 여러 held-out 데이터셋의 평균·최소값을 함께 보고해야 정직합니다.\n\n"
        "## 실무\n"
        "- 최소 3개 이상의 코호트가 있어야 LODO가 의미 있습니다.",
    ),
    (
        "Batch effect / Dataset identifiability",
        "데이터셋을 맞히면 경고음이 울린다",
        "## 정의\n"
        "- **Batch effect**: 실험 일자·장비·시약·파이프라인이 달라서 생기는 샘플 간 조직적 편차.\n"
        "- **Dataset identifiability**: 라벨(BRAF/RAS) 대신 **데이터셋 이름**을 타깃으로 같은 feature로 분류기를 학습했을 때의 AUC.\n\n"
        "## 왜 경고 신호인가\n"
        "- dataset-prediction AUC가 0.99 같으면 feature가 생물학보다 **어느 lab에서 찍었는지**를 학습하고 있다는 뜻입니다.\n"
        "- 이 경우 BRAF 분류 AUC가 높아도 재현성이 낮습니다.\n\n"
        "## 완화책\n"
        "- ComBat, harmony, RUV 등 보정 기법, 또는 robust feature(정규화 강한 유전자)만 선택.",
    ),
    (
        "Isotonic calibration / Platt scaling",
        "확률 재보정 두 가지 방법",
        "## 정의\n"
        "- 분류기의 출력 점수를 **잘 맞는 확률**로 다시 매핑하는 기법입니다.\n\n"
        "## Platt scaling\n"
        "- 점수에 `sigmoid(A·score + B)` 로지스틱을 적합. `A, B` 두 파라미터만 학습해 **작은 데이터**에서 안전.\n\n"
        "## Isotonic regression\n"
        "- 단조 증가 비모수 회귀. 더 유연하지만 **샘플 수가 많아야** 과적합 없음.\n\n"
        "## 언제 쓰나\n"
        "- SVM, tree ensemble의 raw 점수는 확률처럼 사용하면 잘 안 맞습니다(Brier 큼).\n"
        "- reliability diagram 대각선에서 벗어나면 calibration 필요.\n\n"
        "## 주의\n"
        "- calibration은 <strong>validation fold</strong>에서 적합해야지 training에서 적합하면 over-confident.",
    ),
    (
        "Decision curve analysis (DCA) / Net benefit",
        "AUC 너머의 임상 유용성",
        "## 정의\n"
        "- **Net Benefit (NB)** = `TP/n − FP/n × (p_t / (1−p_t))`.\n"
        "- 여기서 `p_t`는 **threshold probability**(환자·의사가 감내할 위험).\n\n"
        "## 읽는 법\n"
        "- x축: 임계확률 `p_t`, y축: net benefit.\n"
        "- 모델 곡선이 **Treat-all / Treat-none** 두 기준선 위에 있는 구간에서 모델이 의사결정에 도움이 됩니다.\n\n"
        "## 왜 AUC만으론 부족한가\n"
        "- AUC는 순위만 보지만 실제 임상은 어느 threshold에서 surgery·biopsy를 결정하는지가 중요합니다.\n"
        "- DCA는 오분류 비용 비대칭을 직접 반영합니다 (literature-reviewed: Vickers & Elkin 2006).",
    ),
    (
        "NPV / PPV at threshold",
        "임상 운영점 지표",
        "## 정의\n"
        "- **PPV (Positive Predictive Value)** = `TP / (TP + FP)`: 양성이라고 말한 환자 중 진짜 양성 비율.\n"
        "- **NPV (Negative Predictive Value)** = `TN / (TN + FN)`: 음성이라고 말한 환자 중 진짜 음성 비율.\n\n"
        "## 왜 threshold가 필요한가\n"
        "- 확률 모델은 연속값을 뱉으므로, 실제 임상에서는 특정 cutoff에서 양·음을 가릅니다.\n\n"
        "## 유병률 의존성\n"
        "- PPV/NPV는 **prevalence**가 바뀌면 같이 바뀝니다. 스크리닝(유병률 낮음) 환경에서 PPV가 뚝 떨어지는 이유.\n"
        "- sensitivity, specificity와 함께 reporting해야 오해가 없습니다.",
    ),
    (
        "Bethesda III/IV triage",
        "indeterminate FNA 분류",
        "## Bethesda System\n"
        "- 갑상선 세침흡인(FNA) 세포검사 결과를 6단계로 표준화한 체계입니다 (literature-reviewed).\n"
        "- **III (AUS/FLUS)**: atypia of undetermined significance. 악성 위험 ~10–30%.\n"
        "- **IV (FN/SFN)**: follicular neoplasm. 악성 위험 ~25–40%.\n\n"
        "## 왜 문제인가\n"
        "- 수많은 환자가 III/IV로 분류돼 **진단 수술**을 받지만, 실제로는 양성인 경우가 대부분입니다.\n"
        "- molecular triage가 NPV를 높이면 **불필요한 갑상선절제 수술**을 줄일 수 있습니다.\n\n"
        "## 지표\n"
        "- Bethesda III/IV만 모은 sub-cohort의 **surgery avoidance %**가 핵심 KPI.",
    ),
    (
        "Log-rank test",
        "KM 곡선 차이의 검정",
        "## 정의\n"
        "- 두(또는 그 이상) 생존 곡선이 **동일한지**를 검정하는 비모수 방법.\n"
        "- 각 event 시점에서 관측/기대 event 수 차이를 카이제곱 형태로 합칩니다.\n\n"
        "## 해석\n"
        "- `p < 0.05`면 두 군의 생존 분포가 통계적으로 다름.\n"
        "- effect size는 말해주지 않으므로 **HR**과 함께 reporting합니다.\n\n"
        "## 가정\n"
        "- proportional hazards: hazard 비율이 시간에 걸쳐 일정해야 합니다. 곡선이 교차하면 log-rank의 검정력이 떨어집니다.",
    ),
    (
        "Concordance index (C-index)",
        "생존 모델 ranking 지표",
        "## 정의\n"
        "- 무작위로 뽑은 두 환자 pair 중 **더 위험한 환자가 더 일찍 event**를 겪는 비율.\n"
        "- 분류의 AUC를 생존 데이터로 일반화한 지표로 볼 수 있습니다.\n\n"
        "## 해석\n"
        "- `0.5` = 무작위, `0.7–0.8` = 쓸만함, `>0.8` = 좋음.\n"
        "- censoring 있는 pair는 비교에서 제외됩니다 (Harrell 1982, literature-reviewed).\n\n"
        "## 주의\n"
        "- C-index는 <strong>calibration</strong>은 보지 않고 **순위**만 봅니다. 절대 위험이 맞는지는 Brier, calibration plot으로 별도 확인.",
    ),
    (
        "Cox proportional hazards",
        "semi-parametric 생존 회귀",
        "## 정의\n"
        "- `h(t | x) = h0(t) × exp(β·x)`.\n"
        "- baseline hazard `h0(t)`는 비모수로 두고, **공변량 효과만 로그-선형**으로 모델링.\n\n"
        "## 장점\n"
        "- baseline 분포(Weibull 등)를 가정할 필요 없음.\n"
        "- `exp(β)`가 곧 **Hazard Ratio (HR)** → 해석 쉬움.\n\n"
        "## 가정\n"
        "- **proportional hazards**: 공변량 효과가 시간에 걸쳐 일정. Schoenfeld residual로 검증.\n"
        "- 위반 시 time-varying coefficient 또는 stratified Cox로 보완.",
    ),
    (
        "Fixed vs random effect meta-analysis",
        "이종성 처리 방식",
        "## Fixed-effect\n"
        "- 모든 연구가 **동일한 진짜 효과 크기**를 공유한다고 가정.\n"
        "- 연구 간 차이는 sampling error 뿐.\n\n"
        "## Random-effect\n"
        "- 각 연구가 자신의 **진짜 효과**를 가지며, 그 효과들이 분포(`τ²`)를 이룬다고 가정.\n"
        "- 코호트·플랫폼이 다르면 보통 random-effect가 정직합니다 (DerSimonian-Laird, literature-reviewed).\n\n"
        "## 선택 기준\n"
        "- 연구가 근본적으로 비슷(동일 프로토콜, 동일 집단) → fixed.\n"
        "- 이질적 데이터(우리 THCA 상황처럼 GEO + TCGA 혼합) → **random-effect** 기본.",
    ),
    (
        "I² statistic",
        "이종성 비율 정량화",
        "## 정의\n"
        "- `I² = max(0, (Q − df) / Q) × 100%`.\n"
        "- 전체 분산 중 **연구 간(between-study) 분산** 비율(%)입니다.\n\n"
        "## 해석\n"
        "- `I² ≈ 25%`: 낮음, `50%`: 중간, `75%+`: 높음 (Higgins 가이드라인, literature-reviewed).\n\n"
        "## 왜 중요한가\n"
        "- I²가 높으면 단순 평균이 아닌 **random-effect pooled estimate**를 써야 합니다.\n"
        "- subgroup 분석이나 meta-regression으로 이종성 원인을 찾는 것이 다음 단계.\n\n"
        "## 주의\n"
        "- 연구 수가 적으면 I² 추정 자체가 불안정합니다. CI를 같이 보고하세요.",
    ),
    (
        "Cochran Q test",
        "이종성의 유의성 검정",
        "## 정의\n"
        "- 여러 연구 효과 크기가 **같은가**를 검정. 귀무가설 `H0: 모든 θ_i 동일`.\n"
        "- 통계량 `Q = Σ w_i (θ_i − θ̂)²`는 자유도 `k−1`의 카이제곱 분포를 따릅니다.\n\n"
        "## 해석\n"
        "- `p < 0.1`이 관례적 이종성 기준(일반 `0.05`보다 느슨).\n"
        "- 유의하면 random-effect 모델로 전환.\n\n"
        "## 한계\n"
        "- 연구 수가 적으면 **검정력 부족**, 많으면 사소한 차이도 유의해집니다. 그래서 **I²와 세트로** 해석합니다.",
    ),
    (
        "Statistical power (통계적 검정력)",
        "1 − β의 의미",
        "## 정의\n"
        "- `power = 1 − β = P(유의 결과 | 효과가 실제로 존재)`.\n"
        "- **진짜 차이가 있을 때 발견할 확률**입니다.\n\n"
        "## 관례\n"
        "- `0.80` (=80%)이 표준 목표, `0.90`이 엄격 목표.\n\n"
        "## 결정 요소\n"
        "- 효과 크기 (Cohen’s d), 표본 수 n, 유의수준 `α`, 양측/단측 여부.\n\n"
        "## 왜 미리 계산해야 하나\n"
        "- post-hoc power는 의미가 적습니다. **사전 power analysis**로 n을 맞추지 않으면 음성 결과를 \"효과 없음\"으로 잘못 해석할 수 있습니다.\n"
        "- meta-analysis power는 study 수와 `τ²`에 의해 결정됩니다.",
    ),
    (
        "ssCMap / connectivity score",
        "화합물 reversal 점수",
        "## 정의\n"
        "- **ssCMap (single-sample Connectivity Map)**: 한 환자의 disease signature와 수천 개 화합물의 유전자 발현 지문을 비교해 얻는 상관 점수.\n"
        "- 점수가 **음수(−)**일수록 질병 시그니처를 **뒤집는** 화합물(=후보 치료제).\n\n"
        "## 해석\n"
        "- 양수(+)는 질병을 <em>강화</em>하는 방향 → off-target 위험 신호.\n"
        "- 절대값이 0.3 이상이면 신호로 간주하는 것이 일반적.\n\n"
        "## 근거\n"
        "- CMap(Lamb 2006), L1000 landmark 978 유전자 기반 확장판(literature-reviewed).\n\n"
        "## 한계\n"
        "- in vitro(HA1E, A549 등) 세포주 데이터 기반이라 갑상선 특이도는 직접 실험 검증 필요.",
    ),
    (
        "PPI (protein-protein interaction)",
        "단백질 상호작용 그래프",
        "## 정의\n"
        "- 단백질 A와 B가 물리적으로 결합하거나 신호전달에서 상호작용하는 관계를 **간선(edge)**으로 표현한 네트워크.\n\n"
        "## 노드와 간선\n"
        "- **노드**: 단백질(=유전자).\n"
        "- **간선**: STRING, BioGRID 등에서 수집한 상호작용. 점수(`confidence`)가 0–1로 표시됩니다.\n\n"
        "## 쓰임\n"
        "- driver 주변의 **hub** 단백질을 찾거나, 신약 표적의 downstream 네트워크를 이해하는 데 씁니다.\n\n"
        "## 주의\n"
        "- 문헌 편향이 큽니다: 연구가 많이 된 단백질일수록 간선이 많아 보입니다. 신규 단백질이 **hub가 아니라는 보장은 없음**.",
    ),
    (
        "Dedifferentiation axis / dediff proxy score",
        "갑상선암 진행의 생물학 축",
        "## 개념\n"
        "- 갑상선 세포는 정상에서 PTC → 저분화(PDTC) → 미분화(ATC)로 가면서 <strong>갑상선 기능 유전자</strong>(TG, TPO, SLC5A5, TSHR 등)가 꺼집니다.\n"
        "- 이 상실 정도를 한 축으로 표현한 것이 **dedifferentiation axis**입니다.\n\n"
        "## dediff proxy score\n"
        "- TDS(Thyroid Differentiation Score)의 역방향 또는 ATC/PDTC 시그니처 기반 점수로 근사.\n"
        "- 높을수록 정상에서 멀어져 더 공격적인 표현형.\n\n"
        "## 임상 함의\n"
        "- RAI (방사성요오드) 치료 반응은 분화가 유지돼야 좋습니다. dediff 높으면 **RAI-refractory** 가능성.\n"
        "- BRAF V600E + 높은 dediff 조합은 예후 감시 대상.",
    ),
]


# ---------------------------------------------------------------------------
# Build the section HTML (single-line style, matching existing glossary).
# ---------------------------------------------------------------------------

def build_section(entries: list[tuple[str, str, str]]) -> str:
    cards: list[str] = []
    for term, summary, body_md in entries:
        term_id = _slug(term)
        body_html = markdown_to_html(body_md)
        term_esc = html.escape(term, quote=False)
        summary_esc = html.escape(summary, quote=False)
        card = (
            f'<div class="glossary-card" id="{term_id}">'
            f'<div class="glossary-term">{term_esc}</div>'
            f'<details class="edu-toggle" data-edu-id="{term_id}">'
            f'<!-- edu:{term_id} -->'
            f'<summary>{summary_esc}</summary>'
            f'<div class="edu-body">{body_html}</div>'
            f'</details>'
            f'</div>'
        )
        cards.append(card)
    return (
        f'{MARKER}'
        f'<section id="cat-v3" class="panel glass" aria-labelledby="hd-cat-v3">'
        f'<h2 class="section-title" id="hd-cat-v3">v3 추가 용어 (v3 additions)</h2>'
        f'<div class="glossary-grid">'
        f'{"".join(cards)}'
        f'</div></section>'
    )


def main() -> int:
    if not GLOSSARY_PATH.exists():
        print(f"[expand_glossary_v3] missing: {GLOSSARY_PATH}")
        return 1
    text = GLOSSARY_PATH.read_text(encoding="utf-8")

    # Idempotency: if marker present, rebuild (replace) the section in-place.
    already_had_marker = MARKER in text

    # Per-term idempotency: drop entries whose id is already present elsewhere.
    filtered: list[tuple[str, str, str]] = []
    skipped: list[str] = []
    for term, summary, body in ENTRIES:
        tid = _slug(term)
        id_attr = f'id="{tid}"'
        # If the id already exists AND it's outside our v3 section, skip it.
        if id_attr in text and MARKER not in text:
            skipped.append(term)
            continue
        filtered.append((term, summary, body))

    new_section = build_section(filtered)

    if already_had_marker:
        # Replace existing v3 section entirely.
        pattern = re.compile(
            re.escape(MARKER) + r"<section[^>]*id=\"cat-v3\".*?</section>",
            re.DOTALL,
        )
        new_text, n = pattern.subn(new_section, text, count=1)
        if n == 0:
            # Marker present but malformed; fall back to insertion.
            new_text = text.replace(MARKER, "") + new_section
    else:
        # Insert before the footer.
        footer_idx = text.find("<footer")
        if footer_idx == -1:
            new_text = text + new_section
        else:
            new_text = text[:footer_idx] + new_section + text[footer_idx:]

    if new_text != text:
        GLOSSARY_PATH.write_text(new_text, encoding="utf-8")
        print(
            f"[expand_glossary_v3] added={len(filtered)} skipped={len(skipped)} "
            f"rebuilt_marker={already_had_marker}"
        )
        if skipped:
            print(f"  skipped-already-present: {skipped}")
    else:
        print("[expand_glossary_v3] no change")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
