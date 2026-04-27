"""AsyncFetchExecutor -- aiohttp-backed multi-URL fetcher with audit trail.

Distilled from the v17 TERT recovery v2 sprint, where 9 sources x ~80 URLs
ran concurrently with per-URL status / bytes / timing logged.
"""
from __future__ import annotations

import asyncio
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

try:
    import aiohttp  # noqa: F401
    HAVE_AIOHTTP = True
except ImportError:
    HAVE_AIOHTTP = False


@dataclass
class FetchResult:
    url: str
    ok: bool
    status: int | None = None
    bytes_: int = 0
    saved_to: str | None = None
    elapsed_seconds: float = 0.0
    error: str | None = None


class AsyncFetchExecutor:
    def __init__(self, save_root: str | Path, ua: str = "agentic-research/0.1"):
        if not HAVE_AIOHTTP:
            raise ImportError("aiohttp required for AsyncFetchExecutor")
        self.save_root = Path(save_root)
        self.save_root.mkdir(parents=True, exist_ok=True)
        self.ua = ua

    async def fetch_one(self, session, url: str, save_name: str | None = None, timeout: int = 30) -> FetchResult:
        import aiohttp
        t0 = time.time()
        try:
            async with session.get(
                url,
                timeout=aiohttp.ClientTimeout(total=timeout),
                headers={"User-Agent": self.ua},
                allow_redirects=True,
            ) as resp:
                data = await resp.read()
                saved = None
                if resp.status == 200 and save_name:
                    p = self.save_root / save_name
                    p.write_bytes(data)
                    saved = str(p)
                return FetchResult(
                    url=url,
                    ok=resp.status == 200,
                    status=resp.status,
                    bytes_=len(data),
                    saved_to=saved,
                    elapsed_seconds=time.time() - t0,
                )
        except asyncio.TimeoutError:
            return FetchResult(url=url, ok=False, error="timeout", elapsed_seconds=time.time() - t0)
        except Exception as e:  # noqa: BLE001
            return FetchResult(url=url, ok=False, error=f"{type(e).__name__}: {e}", elapsed_seconds=time.time() - t0)

    async def run_many(self, urls: list[tuple[str, str | None]], timeout: int = 30) -> list[FetchResult]:
        import aiohttp
        connector = aiohttp.TCPConnector(limit=20, ssl=False)
        async with aiohttp.ClientSession(connector=connector, trust_env=True) as session:
            return await asyncio.gather(
                *[self.fetch_one(session, url, name, timeout) for url, name in urls]
            )
