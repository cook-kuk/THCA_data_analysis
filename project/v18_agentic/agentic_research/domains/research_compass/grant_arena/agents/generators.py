"""Generator agents — each proposes grant topics from a different angle.

Offline scaffold: no LLM calls. Each agent uses simple deterministic rules
over the inputs (CV terms, RFP keywords, datasets, BRIC items, prior art)
to compose plausible candidate topics. The output shape is stable so v0.3
can swap individual generators for LLM-backed implementations.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Protocol, Sequence

from ...cv_profile import ResearchProfile
from ...data_registry import DatasetEntry
from ...method_index import MethodEntry, MethodIndex
from ..mcp_servers.bric_trend import BRICItem
from ..mcp_servers.rfp_reader import RFP


@dataclass
class GrantTopic:
    title: str
    proposed_by: str                         # generator name
    keywords: list[str] = field(default_factory=list)
    suggested_dataset: DatasetEntry | None = None
    suggested_method: MethodEntry | None = None
    rationale: str = ""
    factor_inputs: dict[str, float] = field(default_factory=dict)  # for scoring
    elo: float = 1200.0


class GeneratorAgent(Protocol):
    name: str

    def generate(
        self,
        *,
        profile: ResearchProfile,
        rfp: RFP,
        datasets: Sequence[DatasetEntry],
        method_index: MethodIndex,
        bric_items: Sequence[BRICItem],
        n: int = 5,
    ) -> list[GrantTopic]: ...


def _seed_factors(novelty: float = 0.5, feasibility: float = 0.5,
                  pi_fit: float = 0.5, funder_fit: float = 0.5,
                  publishability: float = 0.5, patentability: float = 0.4,
                  clinical_impact: float = 0.5, scalability: float = 0.5,
                  pilot_speed: float = 0.5, patent_collision_risk: float = 0.1,
                  overclaim_risk: float = 0.1) -> dict[str, float]:
    return {
        "novelty": novelty, "feasibility": feasibility, "pi_fit": pi_fit,
        "funder_fit": funder_fit, "publishability": publishability,
        "patentability": patentability, "clinical_impact": clinical_impact,
        "scalability": scalability, "pilot_speed": pilot_speed,
        "patent_collision_risk": patent_collision_risk,
        "overclaim_risk": overclaim_risk,
    }


class _Base:
    name: str = "base"

    def _topic(self, title: str, *, keywords: Iterable[str],
               dataset: DatasetEntry | None = None,
               method: MethodEntry | None = None,
               rationale: str = "", **factor_overrides: float) -> GrantTopic:
        factors = _seed_factors()
        factors.update(factor_overrides)
        return GrantTopic(
            title=title, proposed_by=self.name,
            keywords=[k for k in keywords if k],
            suggested_dataset=dataset, suggested_method=method,
            rationale=rationale, factor_inputs=factors,
        )


class ClinicalNeedAgent(_Base):
    name = "clinical_need"

    def generate(self, *, profile, rfp, datasets, method_index, bric_items, n=5):
        topics = []
        for topic in profile.prior_topics[:n]:
            topics.append(self._topic(
                title=f"Unmet clinical-need study in {topic}",
                keywords=[topic, "clinical decision", *profile.skills[:2]],
                rationale=f"Anchored to PI prior topic '{topic}'",
                pi_fit=0.85, clinical_impact=0.8, novelty=0.45,
            ))
        return topics[:n]


class FrontierTechAgent(_Base):
    name = "frontier_tech"

    def generate(self, *, profile, rfp, datasets, method_index, bric_items, n=5):
        hits = method_index.search(profile.methods_used + rfp.preferred_keywords, top_k=n) if method_index else []
        topics = []
        for entry, score in hits:
            topics.append(self._topic(
                title=f"{entry.name}-driven study in {profile.prior_topics[0] if profile.prior_topics else 'oncology'}",
                keywords=[entry.name, *entry.domains, *entry.keywords[:2]],
                method=entry,
                rationale=f"Method-fit {score:.2f} from method_index",
                novelty=0.7, funder_fit=0.7, publishability=0.7,
            ))
        return topics[:n]


class CancerBiologyAgent(_Base):
    name = "cancer_biology"

    def generate(self, *, profile, rfp, datasets, method_index, bric_items, n=5):
        cancer_kw = [t for t in profile.prior_topics if "cancer" in t.lower() or "thyroid" in t.lower() or "ptc" in t.lower()]
        out = []
        for kw in cancer_kw[:n]:
            ds = next((d for d in datasets if any(c in d.keywords for c in ["thyroid", "PTC", "cancer"])), None)
            out.append(self._topic(
                title=f"Mechanistic biology of {kw}",
                keywords=[kw, "EMT", "dedifferentiation", "TME"],
                dataset=ds,
                rationale="Cancer biology angle on PI's domain",
                novelty=0.6, clinical_impact=0.75, publishability=0.75,
            ))
        return out[:n]


class BioDataReuseAgent(_Base):
    name = "biodata_reuse"

    def generate(self, *, profile, rfp, datasets, method_index, bric_items, n=5):
        open_ds = [d for d in datasets if d.access == "open"][:n]
        return [
            self._topic(
                title=f"Reusing {d.accession} for {profile.prior_topics[0] if profile.prior_topics else 'novel inference'}",
                keywords=[d.accession, *d.keywords],
                dataset=d,
                rationale=f"High-reuse open dataset ({d.source})",
                feasibility=0.85, pilot_speed=0.8, novelty=0.4,
            )
            for d in open_ds
        ]


class KoreanTrendAgent(_Base):
    name = "korean_trend"

    def generate(self, *, profile, rfp, datasets, method_index, bric_items, n=5):
        return [
            self._topic(
                title=f"Korean replication of: {item.title}",
                keywords=[*item.keywords, "Korean cohort"],
                rationale=f"BRIC trend ({item.korean_pi}, {item.journal})",
                funder_fit=0.7, novelty=0.45, feasibility=0.7,
            )
            for item in (bric_items or [])[:n]
        ]


class PatentWhiteSpaceAgent(_Base):
    name = "patent_white_space"

    def generate(self, *, profile, rfp, datasets, method_index, bric_items, n=5):
        return [self._topic(
            title=f"Patent-differentiated angle on {profile.prior_topics[0] if profile.prior_topics else 'AI in oncology'}",
            keywords=[*profile.methods_used[:2], "claim differentiation"],
            rationale="White-space framing to reduce prior-art collision",
            patentability=0.8, patent_collision_risk=0.15, novelty=0.6,
        )]


class FundingStrategyAgent(_Base):
    name = "funding_strategy"

    def generate(self, *, profile, rfp, datasets, method_index, bric_items, n=5):
        kws = rfp.preferred_keywords[:3]
        if not kws:
            return []
        return [self._topic(
            title=f"{' × '.join(kws)} platform for {profile.prior_topics[0] if profile.prior_topics else 'biomedical research'}",
            keywords=kws,
            rationale=f"Aligned to RFP keywords {kws}",
            funder_fit=0.9, scalability=0.7, novelty=0.55,
        )]


class CrazyIdeaAgent(_Base):
    name = "crazy_idea"

    def generate(self, *, profile, rfp, datasets, method_index, bric_items, n=5):
        return [self._topic(
            title=f"Counterfactual digital twin × {profile.skills[0] if profile.skills else 'AI'} × {profile.prior_topics[0] if profile.prior_topics else 'cancer'}",
            keywords=["counterfactual", "digital twin", *profile.skills[:1]],
            rationale="Out-of-distribution combinatorial proposal",
            novelty=0.9, feasibility=0.4, overclaim_risk=0.25,
        )]
