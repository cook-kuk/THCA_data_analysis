"""thyroid_replay.py -- 200-line reproduction of the v17 thyroid sprint shape.

This is the live-demo script. Running it does NOT re-run the actual analysis;
it walks an example agent through the same five-pattern shape (fix-amplify,
tiered DAG, dual dump, reviewer loop, failure catalog) using mock handlers.
The point is to show the seams.

Run:
    python -m agentic_research.examples.thyroid_replay --max-iter 3
"""
from __future__ import annotations

import argparse
import asyncio
import random
from pathlib import Path

from agentic_research import (
    DualDumpRecorder,
    FailureCatalog,
    Fact,
    FixAmplifyAgent,
    Goal,
    ProjectMemory,
    ReviewerLoop,
    Task,
    TieredParallelOrchestrator,
    Verifier,
)
from agentic_research.patterns.tiered_dag import Tier
from agentic_research.patterns.reviewer_loop import Critique, DefenseEntry
from agentic_research.patterns.failure_catalog import FailureRecord


THYROID_PRIOR = {
    "data_paths": ["project/data_processed/"],
    "panels": ["TierA67_clean", "BRS71", "TDS16", "RAI8"],
    "cohorts": ["TCGA-THCA", "GSE27155", "GSE76039", "GSE126698", "GSE213647"],
    "axes": ["DM1", "DM2"],
    "drivers": ["BRAF", "RAS", "TERT", "TripleNeg"],
    "target_venue": ["npj Precision Oncology", "Genome Medicine"],
}


# --------------------------------------------------------------------------
# Mock handlers for FixAmplifyAgent. Real handlers would call into real code.
# --------------------------------------------------------------------------
async def fix_handler(goal: Goal, last: dict) -> dict:
    print(f"  [FIX] iteration target: {goal.description[:80]}")
    await asyncio.sleep(0.05)
    return {
        "fix_target": "DM1/DM2 cluster definition",
        "auc_estimate": round(0.80 + 0.05 * random.random(), 3),
    }


async def amplify_handler(goal: Goal, last: dict) -> dict:
    print("  [AMPLIFY] extending to external cohorts")
    await asyncio.sleep(0.05)
    return {
        "external_aucs": {
            "GSE76039": round(0.92 + 0.05 * random.random(), 3),
            "GSE126698": round(0.85 + 0.04 * random.random(), 3),
        },
        "tert_recovered_n": 36,
    }


async def synthesis_handler(goal: Goal, last: dict) -> dict:
    print("  [SYNTHESIS] preparing audit + paper section")
    await asyncio.sleep(0.05)
    return {
        "manuscript_section": "Results 4.2 -- 4-group survival",
        "logrank_p": 4.92e-06,
    }


# --------------------------------------------------------------------------
# Tiered DAG -- mimics v17p35 Phase B
# --------------------------------------------------------------------------
async def t1_consensus_clustering():
    await asyncio.sleep(0.05)
    return {"clusters": ["DM1", "DM2"], "n_DM1": 109, "n_DM2": 69}


async def t1_panel_screen():
    await asyncio.sleep(0.05)
    return {"chosen_panel": "TierA67_clean", "auc_cv": 0.954}


async def t2_external_validation():
    await asyncio.sleep(0.05)
    return {"GSE76039_auc": 0.974}


async def t2_survival():
    await asyncio.sleep(0.05)
    return {"logrank_p_4group": 3.78e-05}


async def t3_paper_draft():
    await asyncio.sleep(0.05)
    return {"sections_drafted": 5}


async def t3_reviewer_defense():
    await asyncio.sleep(0.05)
    return {"attacks_defended": 17}


def build_dag() -> TieredParallelOrchestrator:
    return TieredParallelOrchestrator(
        [
            Tier(
                name="discovery",
                tasks=[
                    Task(name="consensus_clustering", func=t1_consensus_clustering),
                    Task(name="panel_screen", func=t1_panel_screen),
                ],
            ),
            Tier(
                name="validation",
                tasks=[
                    Task(name="external_validation", func=t2_external_validation),
                    Task(name="survival", func=t2_survival),
                ],
            ),
            Tier(
                name="paper",
                tasks=[
                    Task(name="paper_draft", func=t3_paper_draft),
                    Task(name="reviewer_defense", func=t3_reviewer_defense),
                ],
            ),
        ]
    )


# --------------------------------------------------------------------------
# Reviewer loop -- mock external reviewer
# --------------------------------------------------------------------------
async def mock_reviewer(artifacts: dict) -> Critique:
    return Critique(
        reviewer_id="external_LLM_v1",
        attack_vectors=[
            {"id": "A1", "title": "TERT promoter calls -- provenance?", "severity": 8},
            {"id": "A2", "title": "DM1/DM2 reproducibility on GSE27155", "severity": 6},
            {"id": "A3", "title": "Confounding by stage", "severity": 5},
        ],
    )


async def mock_defense_builder(critique: Critique, artifacts: dict) -> list[DefenseEntry]:
    out: list[DefenseEntry] = []
    for attack in critique.top_attacks(3):
        out.append(
            DefenseEntry(
                attack_id=attack["id"],
                defense_text=f"Defense for {attack['title']}: see audit + figure.",
            )
        )
    return out


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
async def main(max_iter: int):
    workspace = Path("/tmp/v18_thyroid_replay")
    workspace.mkdir(parents=True, exist_ok=True)

    print("==========================================================")
    print("  v18 thyroid_replay -- agentic_research demo")
    print("==========================================================")

    # Memory + dual-dump + failure catalog
    memory = ProjectMemory(workspace / "memory")
    recorder = DualDumpRecorder(workspace / "dumps")
    catalog = FailureCatalog(workspace / "failures.jsonl")

    catalog.log(
        FailureRecord(
            failure_id="V1_TERT_RECOVERY",
            title="v1 TERT recovery missed thca_tcga_pub",
            what_was_tried="cBioPortal: thca_tcga + thca_tcga_pan_can_atlas_2018",
            why_it_failed="Those two studies use WES MAF; TERT promoter not captured",
            next_action="Iterate ALL cBio thyroid studies, including thca_tcga_pub",
            severity="HIGH",
        )
    )

    # PATTERN 1: fix-amplify-synthesis cycle
    print("\n[PATTERN 1] FixAmplifyAgent")
    goal = Goal(
        description="Reproduce DM1/DM2 -> 8-gene RAI panel -> 4-group survival narrative",
        domain_prior=THYROID_PRIOR,
        success_criteria=["AUC > 0.95", "logrank p < 1e-4"],
        max_iterations=max_iter,
    )
    agent = FixAmplifyAgent(
        goal=goal,
        fix_handler=fix_handler,
        amplify_handler=amplify_handler,
        synthesis_handler=synthesis_handler,
        memory=memory,
        verifier=Verifier(),
    )
    result = await agent.run()
    print(f"  -> iterations: {result.iterations}, success: {result.success}")
    memory.fact(Fact(topic="DM1_DM2", claim="DM1/DM2 reproduces", confidence="HIGH"))

    # PATTERN 2: Tiered parallel DAG
    print("\n[PATTERN 2] TieredParallelOrchestrator")
    dag_summary = await build_dag().run()
    print(f"  -> {dag_summary['ok']}/{dag_summary['n_tasks']} OK across {len(dag_summary['by_tier'])} tiers")

    # PATTERN 3: Dual dump
    print("\n[PATTERN 3] DualDumpRecorder")
    paths = recorder.dump(
        name="phase_B_summary",
        raw_payload=dag_summary,
        compressed_md=(
            "# Phase B summary\n\n"
            f"- Tasks: {dag_summary['n_tasks']}\n"
            f"- OK: {dag_summary['ok']}\n"
            f"- Wall time: {dag_summary['wall_seconds']:.2f}s\n"
        ),
    )
    print(f"  -> raw: {paths.raw.name}, compressed: {paths.compressed.name}")

    # PATTERN 4: Reviewer loop
    print("\n[PATTERN 4] ReviewerLoop")
    loop = ReviewerLoop(reviewers=[mock_reviewer], defense_builder=mock_defense_builder)
    critiques, defenses = await loop.review(result.artifacts)
    print(f"  -> {sum(len(c.attack_vectors) for c in critiques)} attacks, {len(defenses)} defenses")

    # PATTERN 5: Failure catalog
    print("\n[PATTERN 5] FailureCatalog")
    print(f"  -> {len(catalog.list())} entries")

    print("\n==========================================================")
    print(f"  Demo complete. Workspace: {workspace}")
    print("==========================================================")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-iter", type=int, default=3)
    args = parser.parse_args()
    asyncio.run(main(args.max_iter))
