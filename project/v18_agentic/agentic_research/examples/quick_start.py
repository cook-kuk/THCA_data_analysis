"""quick_start.py -- 5-minute hello world."""
from __future__ import annotations

import asyncio

from agentic_research import Goal, FixAmplifyAgent, Verifier


async def fix(goal, last):
    return {"step": "fix", "claim": "found a defect"}


async def amplify(goal, last):
    return {"step": "amplify", "claim": "extended to 3 cohorts"}


async def synthesis(goal, last):
    return {"step": "synthesis", "claim": "wrote the section"}


async def main():
    agent = FixAmplifyAgent(
        goal=Goal(description="Demo run", max_iterations=3),
        fix_handler=fix,
        amplify_handler=amplify,
        synthesis_handler=synthesis,
        verifier=Verifier(),
    )
    res = await agent.run()
    print(res)


if __name__ == "__main__":
    asyncio.run(main())
