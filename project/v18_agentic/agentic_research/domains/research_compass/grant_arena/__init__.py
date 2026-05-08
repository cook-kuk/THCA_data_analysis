"""grant_arena — competitive multi-agent extension of research_compass.

Distilled from a long brainstorm with another LLM ("Grant Arena MCP" /
"Research Future Arena"). Layered on top of research_compass v0:

    research_compass (v0)        grant_arena (v0.2)
    ----------------------       --------------------------------
    cv_profile             →     cv_ontology (capability/asset/evidence split)
    journal_feed (BRIC)    →     bric_trend MCP (Korean PI/topic radar)
    method_index           →     ai_topvenue MCP (NeurIPS/ICML/ICLR/AAAI/CVPR filter)
    data_registry          →     biodata_reuse MCP (TCGA/GEO/DepMap... reuse score)
    topic_ranker (4-factor)→     scoring (10-factor) + arena (Elo + debate + evolve)
                           +     rfp_reader MCP (RFP → evaluation criteria)
                           +     jcr_paper MCP (JCR Q1/top-10% filter)
                           +     patent_priorart MCP (KIPRIS/WIPO/Google Patents screening)

The arena instantiates competing GeneratorAgents (Clinical-Need, Frontier-Tech,
Cancer-Biology, BioData-Reuse, Korean-Trend, Patent-White-Space, Funding-Strategy,
Crazy-Idea), runs pairwise Elo battles, exposes the survivors to CriticAgents
(Feasibility, Novelty, Patent, Reviewer-2, Budget), then EvolutionAgent mutates
the top survivors before final judging.

Marathon-mode: offline-first, no network, no LLM calls. Live wiring is v0.3.
"""
from __future__ import annotations

from .scoring import OpportunityScore, score_topic, ScoringWeights
from .arena import GrantArena, ArenaResult, BattleLog
from .mcp_servers.rfp_reader import RFP, parse_rfp
from .mcp_servers.biodata_reuse import DatasetReuseScore, score_dataset_reuse, default_registry
from .mcp_servers.jcr_paper import is_jcr_q1, journal_tier, JCRRecord
from .mcp_servers.ai_topvenue import is_top_ai_venue, TOP_AI_VENUES
from .mcp_servers.bric_trend import BRICItem, search_bric_hanbitsa
from .mcp_servers.patent_priorart import PriorArtHit, screen_collision_risk

__all__ = [
    "OpportunityScore",
    "score_topic",
    "ScoringWeights",
    "GrantArena",
    "ArenaResult",
    "BattleLog",
    "RFP",
    "parse_rfp",
    "DatasetReuseScore",
    "score_dataset_reuse",
    "default_registry",
    "is_jcr_q1",
    "journal_tier",
    "JCRRecord",
    "is_top_ai_venue",
    "TOP_AI_VENUES",
    "BRICItem",
    "search_bric_hanbitsa",
    "PriorArtHit",
    "screen_collision_risk",
]
