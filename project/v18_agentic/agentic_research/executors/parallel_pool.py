"""ParallelPoolExecutor -- multiprocessing pool for CPU-bound work."""
from __future__ import annotations

import asyncio
from concurrent.futures import ProcessPoolExecutor
from typing import Any, Callable


class ParallelPoolExecutor:
    """Wraps ProcessPoolExecutor in an async-friendly run_many."""

    def __init__(self, max_workers: int | None = None):
        self.max_workers = max_workers

    async def run_many(self, funcs: list[Callable[[], Any]]) -> list[Any]:
        loop = asyncio.get_running_loop()
        with ProcessPoolExecutor(max_workers=self.max_workers) as pool:
            return await asyncio.gather(*[loop.run_in_executor(pool, f) for f in funcs])
