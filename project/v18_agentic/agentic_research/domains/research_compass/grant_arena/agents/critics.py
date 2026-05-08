"""Critic agents — attack proposed topics and adjust factor inputs.

Each critic returns a Critique with a textual attack and a dict of factor
deltas applied to the topic's `factor_inputs` (clamped to [0, 1]).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

from ..mcp_servers.patent_priorart import screen_collision_risk
from .generators import GrantTopic


@dataclass
class Critique:
    by: str                                  # critic name
    attack: str
    factor_deltas: dict[str, float] = field(default_factory=dict)


class CriticAgent(Protocol):
    name: str

    def critique(self, topic: GrantTopic) -> Critique: ...


def _apply(topic: GrantTopic, deltas: dict[str, float]) -> None:
    for k, dv in deltas.items():
        topic.factor_inputs[k] = max(0.0, min(1.0, topic.factor_inputs.get(k, 0.5) + dv))


class FeasibilityCritic:
    name = "feasibility"

    def critique(self, topic: GrantTopic) -> Critique:
        feas = topic.factor_inputs.get("feasibility", 0.5)
        if topic.suggested_dataset is None:
            attack = "No anchored dataset — feasibility unclear."
            deltas = {"feasibility": -0.1}
        elif feas > 0.8:
            attack = "Feasibility claim very high — verify IRB/data access path."
            deltas = {"feasibility": -0.05}
        else:
            attack = "Feasibility plausible at face value."
            deltas = {}
        c = Critique(by=self.name, attack=attack, factor_deltas=deltas)
        _apply(topic, deltas)
        return c


class NoveltyCritic:
    name = "novelty"

    def critique(self, topic: GrantTopic) -> Critique:
        nov = topic.factor_inputs.get("novelty", 0.5)
        if nov < 0.5:
            attack = "Topic is incremental over published reuse studies."
            deltas = {"novelty": -0.1, "publishability": -0.05}
        elif nov > 0.85:
            attack = "Very high novelty — risk of being viewed as speculative."
            deltas = {"overclaim_risk": +0.1}
        else:
            attack = "Novelty in the defensible range."
            deltas = {}
        c = Critique(by=self.name, attack=attack, factor_deltas=deltas)
        _apply(topic, deltas)
        return c


class PatentCritic:
    name = "patent"

    def critique(self, topic: GrantTopic) -> Critique:
        hits = screen_collision_risk(topic.title, topic.keywords)
        if not hits:
            attack = "No representative prior-art collision found in stub corpus."
            deltas = {"patent_collision_risk": -0.05, "patentability": +0.05}
        else:
            top = hits[0]
            attack = f"Possible prior-art overlap: {top.accession} ({top.title}, sim={top.similarity:.2f})"
            deltas = {"patent_collision_risk": +min(0.4, top.similarity), "patentability": -0.05}
        c = Critique(by=self.name, attack=attack, factor_deltas=deltas)
        _apply(topic, deltas)
        return c


class Reviewer2Agent:
    name = "reviewer_2"

    def critique(self, topic: GrantTopic) -> Critique:
        attacks = []
        deltas: dict[str, float] = {}
        if topic.factor_inputs.get("clinical_impact", 0) < 0.5:
            attacks.append("Clinical impact path unclear")
            deltas["clinical_impact"] = -0.05
        if topic.factor_inputs.get("funder_fit", 0) < 0.5:
            attacks.append("Funder-fit weak — would be deprioritized")
            deltas["funder_fit"] = -0.05
        if topic.factor_inputs.get("novelty", 0) > 0.8 and topic.factor_inputs.get("feasibility", 1) < 0.5:
            attacks.append("Ambition outpaces feasibility — overclaim risk")
            deltas["overclaim_risk"] = +0.1
        if not attacks:
            attacks.append("No major reviewer-2 red flags raised on first pass.")
        c = Critique(by=self.name, attack=" | ".join(attacks), factor_deltas=deltas)
        _apply(topic, deltas)
        return c


class BudgetCritic:
    name = "budget"

    def critique(self, topic: GrantTopic) -> Critique:
        scal = topic.factor_inputs.get("scalability", 0.5)
        if scal < 0.4:
            attack = "Scope too small to justify mega-grant budget."
            deltas = {"funder_fit": -0.05}
        elif scal > 0.85:
            attack = "Scope sprawling — risks looking unfocused."
            deltas = {"feasibility": -0.05}
        else:
            attack = "Budget-justification scope plausible."
            deltas = {}
        c = Critique(by=self.name, attack=attack, factor_deltas=deltas)
        _apply(topic, deltas)
        return c
