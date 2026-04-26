"""Inject beginner-friendly education toggles into v3 pages 19–24.

Target pages (created by upstream v3 build):
  19_honesty_audit.html
  20_panel_size_curve.html
  21_three_class.html
  22_survival_v3.html
  23_bethesda.html
  24_multimodal.html

Each page gets a <details class="edu-toggle"> block inserted **after the
first <section>** (hero) using an idempotency marker:
    <!-- edu:v3-<pageN> -->

Notes:
  - Idempotent: re-running does not duplicate the toggle.
  - Ensures the edu-toggle.css stylesheet link is present; if not, appends
    it before </head>.
  - Skips pages that are not (yet) on disk — logs them so the caller can
    report.
  - Does NOT modify serve_secure.py, any Python pipeline, index.html, or
    pages 1–18. Only writes the six v3 pages listed above.
"""
from __future__ import annotations

import html
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGES_DIR = ROOT / "reports" / "html" / "pages"
CSS_LINK = '<link rel="stylesheet" href="../assets/css/edu-toggle.css">'
CSS_MARKER = "<!-- edu-css -->"


# ---------------------------------------------------------------------------
# Markdown -> HTML (same renderer as inject_education_toggles, kept local
# to avoid import-time side effects from that script's __main__ block).
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
    i, n = 0, len(lines)
    while i < n:
        stripped = lines[i].strip()
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
                if s.startswith(("- ", "* ")):
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


# ---------------------------------------------------------------------------
# Per-page content.
# Keys map to target filenames. Each entry has: page_tag (for marker),
# summary (short line shown in collapsed state), title, body_md.
# ---------------------------------------------------------------------------

# The upstream agent ships the pages under slightly different filenames than
# the task spec listed. Map the canonical spec keys to any alias filenames we
# have observed on disk so the injector works for both naming conventions.
FILENAME_ALIASES: dict[str, list[str]] = {
    "20_panel_size_curve.html": ["20_panel_size_curve.html", "20_panel_size.html"],
    "22_survival_v3.html": ["22_survival_v3.html", "22_survival.html"],
    "23_bethesda.html": ["23_bethesda.html", "23_bethesda_sim.html"],
    # 19, 21, 24 used in both spec and actual build.
    "19_honesty_audit.html": ["19_honesty_audit.html"],
    "21_three_class.html": ["21_three_class.html"],
    "24_multimodal.html": ["24_multimodal.html"],
}


PAGE_CONTENT: dict[str, dict] = {
    "19_honesty_audit.html": {
        "page_tag": "v3-19",
        "summary": "왜 honesty audit이 필요한가 (처음 오신 분을 위한 안내)",
        "title": "Honesty audit 읽는 법",
        "body_md": (
            "## 왜 이 페이지가 존재하는가\n"
            "- 단일 코호트에서 계산한 완벽에 가까운 AUC(예: `0.99`)는 **진짜 실력**이 아니라 leakage·batch·label noise의 합작품일 수 있습니다.\n"
            "- honesty audit은 \"우리 숫자가 부풀려졌는가?\"를 체계적으로 자가 점검하는 섹션입니다.\n\n"
            "## 이 페이지의 구성 요소\n"
            "- **Leakage curve (누출 곡선)**: train/test가 겹칠 때 AUC가 얼마나 상승하는지 보여줍니다. leak=0에서의 값이 정직한 추정치입니다.\n"
            "- **Permutation null (순열 null 분포)**: 라벨을 무작위로 섞어 학습했을 때의 AUC 히스토그램. 실제 AUC가 이 분포의 꼬리에 있으면 유의, 중앙 근처면 의심.\n"
            "- **Dataset identifiability**: 라벨 대신 **데이터셋 이름**을 맞히는 분류기 성능. `AUC>0.9`면 모델이 생물학이 아닌 batch를 학습 중일 가능성.\n"
            "- **Decision curve (DCA)**: ROC 너머의 **임상 순이익(net benefit)**. 같은 AUC라도 threshold·비용 비대칭에 따라 실제 유용성이 달라집니다.\n\n"
            "## 읽는 순서 제안\n"
            "1. 먼저 permutation null에서 실제 AUC의 p-value를 확인.\n"
            "2. leakage curve의 slope가 가파르면 group-wise split을 강화해야 합니다.\n"
            "3. dataset identifiability가 높으면 ComBat·harmony 보정 재검토.\n"
            "4. DCA로 임상 threshold 구간에서 모델이 treat-all/treat-none보다 **위에 있는지** 확인.\n\n"
            "## 용어\n"
            "- 모든 약어는 상단 `용어` 페이지의 *v3 추가 용어* 섹션에서 찾을 수 있습니다."
        ),
    },
    "20_panel_size_curve.html": {
        "page_tag": "v3-20",
        "summary": "panel 크기 vs 성능 — 얼마나 작은 패널이 충분한가",
        "title": "Panel size curve 읽는 법",
        "body_md": (
            "## 이 페이지의 질문\n"
            "- BRAF-like vs RAS-like 분류에 정말로 전체 전사체(~20,000 유전자)가 필요한가?\n"
            "- 16, 32, 67, 112 유전자 패널 중 **어디에서 성능 곡선이 평탄해지는가**?\n\n"
            "## Occam의 면도날 관점\n"
            "- 동일한 성능이면 **더 단순한 모델**이 임상적으로 우월합니다.\n"
            "- 유전자 수가 적을수록 qPCR, NanoString 같은 저비용 플랫폼에서 임상화가 쉽습니다.\n\n"
            "## Plateau 탐지\n"
            "- 곡선의 기울기가 작아지기 시작하는 **팔꿈치(elbow)** 지점이 실용적 최소 패널 크기입니다.\n"
            "- 숫자적으로는 패널을 늘렸을 때 AUC 상승 폭이 `<0.005` 수준이면 추가 이득이 미미하다고 봅니다.\n\n"
            "## 임상 롤아웃 함의\n"
            "- 작은 패널은 **비용·턴어라운드·reproducibility** 모두에서 유리합니다.\n"
            "- 단, 외부 코호트에서도 같은 plateau가 재현되는지(LODO 검증)가 전제입니다.\n\n"
            "## 주의\n"
            "- 같은 AUC라도 **calibration**과 **DCA 순이익**까지 함께 봐야 진짜 임상 가치를 판단할 수 있습니다."
        ),
    },
    "21_three_class.html": {
        "page_tag": "v3-21",
        "summary": "왜 3-class(BRAF / RAS / dediff)로 나누는가",
        "title": "Three-class 분류 읽는 법",
        "body_md": (
            "## 왜 2-class가 아니라 3-class인가\n"
            "- **BRAF-like vs RAS-like** 이분법은 많은 PTC를 잘 설명하지만, **dedifferentiated (PDTC/ATC-like)** 샘플은 별도 축에 있습니다.\n"
            "- dediff 축은 예후가 현저히 나쁘고 RAI(방사성요오드) 저항성과 연관되어, 임상에서 분리해 보는 것이 유용합니다.\n\n"
            "## 생물학적 의미\n"
            "- dedifferentiation은 갑상선 기능 유전자(TG, TPO, TSHR, SLC5A5)의 소실을 동반합니다.\n"
            "- BRAF-like 중에서도 dediff로 진행한 샘플은 별도 치료 전략이 필요합니다.\n\n"
            "## Softmax vs One-vs-Rest (OvR)\n"
            "- **Softmax**: 세 클래스 확률의 합 = 1. 서로 **배타적**이라고 보는 모델.\n"
            "- **OvR**: 각 클래스 vs 나머지로 별도 이진 분류기를 학습. 클래스가 **중첩**될 때(BRAF-like + dediff) 더 자연스럽습니다.\n"
            "- 우리 데이터에서는 BRAF-like이면서 dediff인 샘플이 존재하므로 OvR 결과도 함께 보고합니다.\n\n"
            "## 읽을 때 주의\n"
            "- 세 클래스 AUC를 단순 평균하면 **macro-AUC**, 샘플 비례 가중이면 **micro-AUC**. 소수 클래스(dediff) 성능이 중요하므로 macro를 우선 확인."
        ),
    },
    "22_survival_v3.html": {
        "page_tag": "v3-22",
        "summary": "생존 분석 입문 — KM, log-rank, Cox",
        "title": "Survival curve 읽는 법",
        "body_md": (
            "## Kaplan-Meier (KM) 곡선\n"
            "- x축: 시간(개월·년), y축: **생존 확률** `S(t)`.\n"
            "- 계단이 떨어질 때마다 event(사망 또는 재발)가 발생한 것이며, **작은 틱 마크**는 censored(관찰 종료) 환자입니다.\n\n"
            "## Log-rank 검정\n"
            "- 두 곡선이 통계적으로 다른지 판정합니다. `p<0.05`이면 생존 분포가 다르다는 증거.\n"
            "- **효과 크기**는 말해주지 않으므로 HR과 함께 reporting합니다.\n\n"
            "## Cox proportional hazards\n"
            "- `h(t|x) = h0(t) × exp(β·x)`. `exp(β)` = **Hazard Ratio**.\n"
            "- HR `1.5`이면 해당 군의 순간 위험률이 1.5배.\n"
            "- proportional hazards 가정(시간에 걸쳐 효과 일정)이 깨지면 Schoenfeld residual로 감지합니다.\n\n"
            "## THCA 특성: 저사건율(censoring)\n"
            "- PTC는 전반적으로 **5년 생존율 >98%**로, event가 매우 적습니다.\n"
            "- 이 때문에 KM 곡선 간 차이가 작고, **power**가 부족해 보일 수 있습니다.\n"
            "- 해결책: DFS(disease-free survival) 같은 더 빈번한 endpoint, 더 긴 추적, pooled 코호트.\n\n"
            "## 읽을 때 주의\n"
            "- 곡선 말단은 샘플 수가 급감해 **신뢰구간이 매우 넓어집니다**. 끝부분의 차이는 과해석 금물."
        ),
    },
    "23_bethesda.html": {
        "page_tag": "v3-23",
        "summary": "Bethesda FNA 분류와 NPV/PPV 트레이드오프",
        "title": "Bethesda triage 읽는 법",
        "body_md": (
            "## Bethesda System이란\n"
            "- 갑상선 세침흡인(FNA) 세포검사 결과를 **6단계**로 표준화한 체계 (literature-reviewed).\n"
            "- 단계별 악성 위험이 다르며, III/IV는 **indeterminate**(불확실)로 분류됩니다.\n\n"
            "## 왜 III/IV가 임상 문제인가\n"
            "- **III (AUS/FLUS)** 악성 위험 ~10–30%, **IV (FN/SFN)** ~25–40%.\n"
            "- 최종 확진을 위해 많은 환자가 **진단 수술(diagnostic lobectomy)**을 받지만, 대부분 양성 판정으로 밝혀집니다.\n"
            "- 즉, **불필요한 수술률**이 핵심 KPI입니다.\n\n"
            "## NPV vs PPV 트레이드오프\n"
            "- **NPV (음성 예측값)**: 모델이 \"양성 아님\"이라 한 환자 중 실제 양성이 아닐 비율. `NPV>95%`이면 수술을 안전하게 피할 수 있는 **rule-out** 도구.\n"
            "- **PPV (양성 예측값)**: 모델이 \"양성\"이라 한 환자 중 진짜 양성 비율. **rule-in** 용도.\n"
            "- Bethesda III/IV triage의 1차 목표는 보통 **높은 NPV**입니다(수술을 피하는 용도).\n\n"
            "## 유병률 의존성\n"
            "- NPV/PPV는 prevalence에 따라 크게 달라집니다. 외부 코호트의 prevalence가 다르면 숫자가 바뀝니다. sensitivity/specificity와 함께 reporting하세요.\n\n"
            "## 읽을 때 주의\n"
            "- 모델이 cPTC까지 전부 음성으로 보내지 않도록 **confusion matrix**와 sub-category별 민감도를 확인해야 합니다."
        ),
    },
    "24_multimodal.html": {
        "page_tag": "v3-24",
        "summary": "멀티오믹스 fusion — 어떤 modality가 어떤 이야기를 하는가",
        "title": "Multimodal fusion 읽는 법",
        "body_md": (
            "## 모달리티별로 무엇을 보는가\n"
            "- **RNA (전사체)**: 현재 세포가 **실제로 만들고 있는** mRNA 수준. BRAF/RAS subtype 분리의 핵심 축.\n"
            "- **Methylation (DNA 메틸화)**: **epigenetic memory**. 정상 조직과의 거리, 분화 상실 정도를 장기 관점으로 보여줍니다.\n"
            "- **SCNA (somatic copy-number alteration)**: 염색체 구간 복제 수 변화. 대규모 유전체 불안정성(ATC 방향)을 포착.\n"
            "- **Fusion / mutation**: BRAF V600E, RET 융합 등 **driver 이벤트**. 치료 표적과 직결.\n\n"
            "## 왜 fusion 모델이 필요한가\n"
            "- 한 모달리티가 놓치는 환자를 다른 모달리티가 포착합니다(예: methylation-only dediff).\n"
            "- late fusion(각 모달 확률의 가중합)과 early fusion(특징 concat 후 학습)이 각각 장단점.\n\n"
            "## 교차 코호트(LODO) 검증이 왜 gold standard인가\n"
            "- 한 코호트 내에서 random CV를 돌리면 **batch-specific 신호**를 학습해 AUC가 과대 평가됩니다.\n"
            "- **LODO (Leave-One-Dataset-Out)**는 실제 임상 배포(새 병원, 새 플랫폼)와 가장 가까운 시나리오.\n"
            "- LODO AUC가 내부 CV AUC보다 낮은 것이 정상입니다. 두 값을 나란히 보고해야 정직합니다.\n\n"
            "## 읽을 때 주의\n"
            "- 모달리티 추가가 항상 도움 되는 것은 아닙니다. **marginal gain**이 작고 샘플이 줄면 오히려 손해일 수 있습니다(missing-by-modality 패널티).\n"
            "- feature importance는 fusion 모델에서 **모달 간 상관** 때문에 해석이 까다롭습니다. per-modality 단독 AUC도 reporting하세요."
        ),
    },
}


# ---------------------------------------------------------------------------
# Injection helpers.
# ---------------------------------------------------------------------------

def _ensure_css(text: str) -> str:
    """Ensure the edu-toggle.css link is present in <head>."""
    if 'edu-toggle.css' in text:
        return text
    # Prefer inserting via existing marker convention used elsewhere.
    if CSS_MARKER in text:
        return text.replace(CSS_MARKER, f"{CSS_MARKER}{CSS_LINK}")
    # Fallback: inject right before </head>.
    head_close = text.find("</head>")
    if head_close == -1:
        return text
    tagged = f"<!-- edu-css -->{CSS_LINK}"
    return text[:head_close] + tagged + text[head_close:]


def _build_toggle_block(page_tag: str, title: str, summary: str, body_md: str) -> str:
    body_html = markdown_to_html(body_md)
    summary_esc = html.escape(summary, quote=False)
    title_esc = html.escape(title, quote=False)
    return (
        f'<details class="edu-toggle" data-edu-id="{page_tag}">'
        f'<!-- edu:{page_tag} -->'
        f'<summary>{summary_esc}</summary>'
        f'<div class="edu-body">'
        f'<h3>{title_esc}</h3>'
        f'{body_html}'
        f'</div>'
        f'</details>'
    )


_FIRST_SECTION_CLOSE_RE = re.compile(r"</section>", re.IGNORECASE)


def _inject_after_first_section(text: str, block: str) -> str | None:
    """Insert ``block`` right after the first </section> occurrence.

    Returns the new text, or None if no <section> was found.
    """
    m = _FIRST_SECTION_CLOSE_RE.search(text)
    if not m:
        return None
    insert_at = m.end()
    return text[:insert_at] + block + text[insert_at:]


def _already_injected(text: str, page_tag: str) -> bool:
    return f"<!-- edu:{page_tag} -->" in text


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------

def main() -> int:
    injected: list[str] = []
    already: list[str] = []
    missing: list[str] = []
    no_section: list[str] = []

    for canonical_fname, info in PAGE_CONTENT.items():
        # Resolve which alias actually exists on disk.
        aliases = FILENAME_ALIASES.get(canonical_fname, [canonical_fname])
        resolved: Path | None = None
        resolved_name: str = canonical_fname
        for alias in aliases:
            candidate = PAGES_DIR / alias
            if candidate.exists():
                resolved = candidate
                resolved_name = alias
                break
        if resolved is None:
            missing.append(canonical_fname)
            continue

        text = resolved.read_text(encoding="utf-8")
        original = text

        tag = info["page_tag"]
        if _already_injected(text, tag):
            # Still ensure CSS link present (cheap, idempotent).
            text = _ensure_css(text)
            if text != original:
                resolved.write_text(text, encoding="utf-8")
            already.append(resolved_name)
            continue

        text = _ensure_css(text)
        block = _build_toggle_block(
            page_tag=tag,
            title=info["title"],
            summary=info["summary"],
            body_md=info["body_md"],
        )
        new_text = _inject_after_first_section(text, block)
        if new_text is None:
            no_section.append(resolved_name)
            continue

        resolved.write_text(new_text, encoding="utf-8")
        injected.append(resolved_name)

    print(f"[inject_v3_edutoggles] injected={len(injected)} "
          f"already={len(already)} missing={len(missing)} "
          f"no_section={len(no_section)}")
    if injected:
        print(f"  injected: {injected}")
    if already:
        print(f"  already-had-marker: {already}")
    if missing:
        print(f"  missing-from-disk: {missing}")
    if no_section:
        print(f"  no-<section>-tag: {no_section}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
