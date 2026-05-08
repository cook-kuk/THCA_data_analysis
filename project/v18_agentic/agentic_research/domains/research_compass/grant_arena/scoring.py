"""10-factor composite scoring for grant topics.

Extends research_compass v0's 4-factor weighted sum (competence/novelty/
feasibility/method_fit) to a 10-factor scheme tuned for mega-grant
discovery. Weights mirror the brainstorm spec; tweak via ScoringWeights.

    score = + 0.18 novelty
            + 0.15 feasibility
            + 0.15 pi_fit
            + 0.12 funder_fit          (e.g. Samsung Future Tech)
            + 0.10 publishability
            + 0.10 patentability
            + 0.10 clinical_impact
            + 0.05 scalability
            + 0.05 pilot_speed
            - 0.15 patent_collision_risk
            - 0.10 overclaim_risk
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ScoringWeights:
    novelty: float = 0.18
    feasibility: float = 0.15
    pi_fit: float = 0.15
    funder_fit: float = 0.12
    publishability: float = 0.10
    patentability: float = 0.10
    clinical_impact: float = 0.10
    scalability: float = 0.05
    pilot_speed: float = 0.05
    patent_collision_risk: float = 0.15  # subtracted
    overclaim_risk: float = 0.10         # subtracted


@dataclass
class OpportunityScore:
    total: float
    components: dict[str, float] = field(default_factory=dict)
    tier: str = "LOW"  # HIGH ≥0.55, MEDIUM ≥0.30, LOW <0.30

    def explain(self) -> str:
        parts = ", ".join(f"{k}={v:.2f}" for k, v in self.components.items())
        return f"{self.tier} (total={self.total:.3f}) | {parts}"


def _tier_for(score: float) -> str:
    if score >= 0.55:
        return "HIGH"
    if score >= 0.30:
        return "MEDIUM"
    return "LOW"


def score_topic(
    *,
    novelty: float,
    feasibility: float,
    pi_fit: float,
    funder_fit: float,
    publishability: float,
    patentability: float,
    clinical_impact: float,
    scalability: float = 0.5,
    pilot_speed: float = 0.5,
    patent_collision_risk: float = 0.0,
    overclaim_risk: float = 0.0,
    weights: ScoringWeights | None = None,
) -> OpportunityScore:
    """All factor inputs are in [0, 1]. Risk factors are subtracted."""
    w = weights or ScoringWeights()
    components = {
        "novelty": novelty,
        "feasibility": feasibility,
        "pi_fit": pi_fit,
        "funder_fit": funder_fit,
        "publishability": publishability,
        "patentability": patentability,
        "clinical_impact": clinical_impact,
        "scalability": scalability,
        "pilot_speed": pilot_speed,
        "patent_collision_risk": patent_collision_risk,
        "overclaim_risk": overclaim_risk,
    }
    total = (
        w.novelty * novelty
        + w.feasibility * feasibility
        + w.pi_fit * pi_fit
        + w.funder_fit * funder_fit
        + w.publishability * publishability
        + w.patentability * patentability
        + w.clinical_impact * clinical_impact
        + w.scalability * scalability
        + w.pilot_speed * pilot_speed
        - w.patent_collision_risk * patent_collision_risk
        - w.overclaim_risk * overclaim_risk
    )
    return OpportunityScore(
        total=round(total, 3),
        components={k: round(v, 3) for k, v in components.items()},
        tier=_tier_for(total),
    )
