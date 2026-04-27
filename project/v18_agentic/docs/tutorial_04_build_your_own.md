# Tutorial 4 — Build your own agent in 30 minutes

You'll build a small agent that grades a markdown manuscript draft against
a checklist. By the end you'll have used Goal, FixAmplifyAgent, Verifier,
ProjectMemory, and ReviewerLoop.

## Setup (1 min)

```bash
pip install -e .
```

## Step 1 — Define the goal (3 min)

```python
from agentic_research import Goal

goal = Goal(
    description="Manuscript draft passes the npj checklist",
    domain_prior={
        "checklist": [
            "abstract <= 150 words",
            "primary cohort declared in section 2",
            "external validation in section 3",
            "limitations declared",
            "data availability statement present",
        ],
    },
    success_criteria=["all checklist items present"],
    max_iterations=4,
)
```

## Step 2 — Write a verifier (5 min)

```python
from agentic_research import Verifier

def has_section(name):
    def check(artifacts):
        text = artifacts.get("draft_text", "")
        ok = name.lower() in text.lower()
        return ok, f"section '{name}' present" if ok else f"section '{name}' missing"
    return check

verifier = Verifier(checks=[
    has_section("Abstract"),
    has_section("External validation"),
    has_section("Limitations"),
    has_section("Data availability"),
])
```

## Step 3 — Wire up FixAmplifyAgent (10 min)

```python
from agentic_research import FixAmplifyAgent, ProjectMemory

DRAFT = "# Title\n\n## Abstract\n...\n## Methods\n"

async def fix(goal, last):
    return {"draft_text": DRAFT, "step": "fix"}

async def amplify(goal, last):
    text = last.get("draft_text", DRAFT) + "\n## External validation\n...\n"
    return {"draft_text": text, "step": "amplify"}

async def synthesis(goal, last):
    text = last.get("draft_text", DRAFT) + "\n## Limitations\n## Data availability\n"
    return {"draft_text": text, "step": "synthesis"}

agent = FixAmplifyAgent(
    goal=goal,
    fix_handler=fix, amplify_handler=amplify, synthesis_handler=synthesis,
    memory=ProjectMemory("./manuscript_workspace"),
    verifier=verifier,
)
```

## Step 4 — Run it (1 min)

```python
import asyncio
result = asyncio.run(agent.run())
print(result)
```

## Step 5 — Inspect the audit trail (5 min)

```bash
ls manuscript_workspace/audit/
# 1729...0_iter_0.json  1729...0_iter_1.json  1729...0_iter_2.json
```

Each iteration's artifacts and verifier verdict are recorded. The verifier
notes per iteration are how the agent decides whether to continue or stop.

## Step 6 — Add a reviewer loop (5 min)

```python
from agentic_research import ReviewerLoop
from agentic_research.patterns.reviewer_loop import Critique, DefenseEntry

async def style_reviewer(artifacts):
    text = artifacts.get("draft_text", "")
    attacks = []
    if len(text.split()) > 4500:
        attacks.append({"id": "WORDCOUNT", "title": "draft over 4500 words", "severity": 5})
    return Critique(reviewer_id="style_v1", attack_vectors=attacks)

async def defense_builder(critique, artifacts):
    return [DefenseEntry(attack_id=a["id"], defense_text="Will trim.") for a in critique.top_attacks(3)]

loop = ReviewerLoop(reviewers=[style_reviewer], defense_builder=defense_builder)
critiques, defenses = asyncio.run(loop.review(result.artifacts))
```

## What you've built

A 50-line agent that:
- Knows what "done" means (Goal + Verifier).
- Iterates fix → amplify → synthesis until done (FixAmplifyAgent).
- Records every iteration's audit (ProjectMemory).
- Runs an external reviewer over the final artifact (ReviewerLoop).

That's the same shape the v17 thyroid sprint used at multi-month scale —
just smaller.
