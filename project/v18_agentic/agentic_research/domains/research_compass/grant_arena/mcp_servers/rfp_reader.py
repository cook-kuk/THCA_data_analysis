"""rfp_reader_mcp — parse a funding RFP into structured criteria.

In offline mode, accepts free-text and applies a few keyword heuristics
to populate the RFP dataclass. Live mode (v0.3) would use an LLM tool.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass
class RFP:
    program: str = ""
    track: str = ""               # "Science" | "Technology" | ...
    budget_scale: str = ""        # e.g. "30억"
    duration_years: int = 0
    must_have: list[str] = field(default_factory=list)
    preferred_keywords: list[str] = field(default_factory=list)
    forbidden_or_risky: list[str] = field(default_factory=list)
    raw_text: str = ""


_BUDGET_RE = re.compile(r"(\d+)\s*억")
_YEARS_RE = re.compile(r"(\d+)\s*년")
_KOREAN_TECH_KEYWORDS = (
    "Advanced AI", "차세대 로봇", "디지털 헬스", "정밀의료",
    "유전자 치료제", "단백질 치료제", "바이오 공정", "바이오 융합",
    "양자", "Foundation", "파운데이션", "원천기술",
)


def parse_rfp(text: str, *, program: str = "", track: str = "") -> RFP:
    rfp = RFP(program=program, track=track, raw_text=text)
    if (m := _BUDGET_RE.search(text)):
        rfp.budget_scale = f"{m.group(1)}억"
    if (m := _YEARS_RE.search(text)):
        rfp.duration_years = int(m.group(1))
    rfp.preferred_keywords = [kw for kw in _KOREAN_TECH_KEYWORDS if kw in text]

    must_have_hints = ["창의성", "원천성", "도전성", "차별성", "실행가능성", "산업 임팩트"]
    rfp.must_have = [h for h in must_have_hints if h in text]

    risky_phrases = ["RCT", "임상시험 자동화", "FDA 승인", "임상 의사결정 자동화"]
    rfp.forbidden_or_risky = [r for r in risky_phrases if r.lower() in text.lower()]

    return rfp
