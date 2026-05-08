"""Stage 1 — CV → ResearchProfile.

Parses a researcher CV (markdown or plain text) into structured signals the
ranker can join against journal feeds, method papers, and dataset metadata.

Design: keep extraction *lossy on purpose*. We want a few high-confidence
buckets (skills, prior topics, methods used, datasets touched), not a full
ontology. Ambiguity is fine — the ranker tolerates noisy keywords.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ResearchProfile:
    """Structured fingerprint of a researcher."""

    name: str
    skills: list[str] = field(default_factory=list)            # e.g. "RNA-seq", "deep learning"
    prior_topics: list[str] = field(default_factory=list)      # e.g. "thyroid cancer", "HLA"
    methods_used: list[str] = field(default_factory=list)      # e.g. "ComBat", "LODO CV"
    datasets_used: list[str] = field(default_factory=list)     # e.g. "TCGA-THCA", "GSE286332"
    keywords: list[str] = field(default_factory=list)          # union, deduped, lowercased

    def all_terms(self) -> set[str]:
        """Tokenize every field down to atomic terms (hyphen-preserving) and union.

        CV entries are usually full phrases ("RNA-seq analysis (bulk + single-cell)");
        journal-item keywords are atomic ("RNA-seq"). Without tokenization the
        Jaccard overlap is always 0. Hyphens are preserved so "RNA-seq" survives
        as a single token.
        """
        out: set[str] = set()
        for s in self.skills + self.prior_topics + self.methods_used + self.datasets_used + self.keywords:
            if not s:
                continue
            for tok in re.split(r"[^A-Za-z0-9\-_+]+", s):
                tok = tok.strip("-_").lower()
                if len(tok) >= 3:
                    out.add(tok)
        return out


_SECTION_HEADS = {
    "skills":   re.compile(r"^#+\s*(skills|technical skills|expertise|기술|역량)", re.I | re.M),
    "topics":   re.compile(r"^#+\s*(research|publications|topics|주제|연구)", re.I | re.M),
    "methods":  re.compile(r"^#+\s*(methods|techniques|방법)", re.I | re.M),
    "datasets": re.compile(r"^#+\s*(datasets|data|cohorts|코호트|데이터)", re.I | re.M),
}

_BULLET = re.compile(r"^[\s]*[-*•]\s+(.+)$", re.M)
_GENE_LIKE = re.compile(r"\b(GSE\d+|TCGA[-_][A-Z]+|UKB|GeneBass|GTEx|cBioPortal|PRJ[A-Z]+\d+|GEO|SRA|EGA)\b")


def _extract_section(text: str, head_re: re.Pattern) -> str:
    m = head_re.search(text)
    if not m:
        return ""
    start = m.end()
    next_head = re.search(r"^#+\s+\S", text[start:], re.M)
    return text[start:start + next_head.start()] if next_head else text[start:]


def _bullets(section: str) -> list[str]:
    return [b.strip() for b in _BULLET.findall(section) if b.strip()]


def load_cv(path_or_text: str | Path, name: str | None = None) -> ResearchProfile:
    """Load CV from path or inline text. Returns ResearchProfile.

    Section heads recognized (English + Korean): Skills/기술, Research/연구,
    Methods/방법, Datasets/데이터. Anything outside falls into ``keywords``
    via a simple GSE/TCGA/UKB regex sweep.
    """
    if isinstance(path_or_text, Path) or (isinstance(path_or_text, str) and len(path_or_text) < 256 and Path(path_or_text).exists()):
        text = Path(path_or_text).read_text()
        derived_name = Path(path_or_text).stem
    else:
        text = str(path_or_text)
        derived_name = "anonymous"

    profile = ResearchProfile(name=name or derived_name)
    profile.skills = _bullets(_extract_section(text, _SECTION_HEADS["skills"]))
    profile.prior_topics = _bullets(_extract_section(text, _SECTION_HEADS["topics"]))
    profile.methods_used = _bullets(_extract_section(text, _SECTION_HEADS["methods"]))
    profile.datasets_used = _bullets(_extract_section(text, _SECTION_HEADS["datasets"]))
    profile.keywords = sorted(set(_GENE_LIKE.findall(text)))

    if not any([profile.skills, profile.prior_topics, profile.methods_used, profile.datasets_used]):
        profile.keywords = list(set(profile.keywords) | {w.lower() for w in re.findall(r"\b[A-Z][a-z]{4,}\b", text)[:20]})

    return profile
