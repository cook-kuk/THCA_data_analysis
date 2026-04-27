"""Shared helpers for v17 TERT recovery v2 (9 sources)."""
from __future__ import annotations

import asyncio
import json
import logging
import re
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Iterable

import aiohttp

REPO = Path(__file__).resolve().parents[1].parent
V2_DIR = REPO / "project" / "results" / "v17_tert_recovery" / "v2"
RAW = V2_DIR / "raw"
PARSED = V2_DIR / "parsed"
LOGS = V2_DIR / "logs"
for d in (V2_DIR, RAW, PARSED, LOGS):
    d.mkdir(parents=True, exist_ok=True)

# TERT promoter coordinates (chr5)
PROMOTER_GRCH37 = {"C228T": 1295228, "C250T": 1295250}
PROMOTER_GRCH38 = {"C228T": 1295113, "C250T": 1295135}
PROMOTER_REGION_BUFFER = 250  # ±bp tolerance for liftover/annotation drift

UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/124.0 Safari/537.36 v17-tert-recovery"
)


def configure_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(LOGS / f"{name}.log", mode="w")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(fh)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter(f"[{name}] %(message)s"))
    logger.addHandler(sh)
    return logger


@dataclass
class Attempt:
    url: str
    ok: bool = False
    status: int | None = None
    bytes_: int = 0
    saved_to: str | None = None
    error: str | None = None
    elapsed: float = 0.0


@dataclass
class SourceResult:
    source_id: str
    label: str
    n_records: int = 0
    n_promoter_mutations: int = 0
    n_tcga_matched: int = 0
    attempts: list[Attempt] = field(default_factory=list)
    output_files: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    success: bool = False

    def to_json(self) -> dict:
        d = asdict(self)
        return d


async def fetch_url(
    session: aiohttp.ClientSession,
    url: str,
    save_path: Path | None = None,
    timeout: int = 30,
    expect_binary: bool | None = None,
    headers: dict | None = None,
) -> Attempt:
    """Fetch URL, optionally save to disk. Returns Attempt with metadata."""
    t0 = time.time()
    h = {"User-Agent": UA}
    if headers:
        h.update(headers)
    try:
        async with session.get(url, timeout=aiohttp.ClientTimeout(total=timeout), headers=h, allow_redirects=True) as resp:
            data = await resp.read()
            elapsed = time.time() - t0
            ct = resp.headers.get("Content-Type", "")
            saved = None
            if resp.status == 200 and save_path is not None:
                save_path.parent.mkdir(parents=True, exist_ok=True)
                save_path.write_bytes(data)
                saved = str(save_path.relative_to(REPO))
            return Attempt(
                url=url,
                ok=resp.status == 200,
                status=resp.status,
                bytes_=len(data),
                saved_to=saved,
                elapsed=elapsed,
            )
    except asyncio.TimeoutError:
        return Attempt(url=url, ok=False, error="timeout", elapsed=time.time() - t0)
    except Exception as e:  # noqa: BLE001
        return Attempt(url=url, ok=False, error=f"{type(e).__name__}: {e}", elapsed=time.time() - t0)


def is_promoter_position(chrom: str | int, pos: int) -> bool:
    """Loose match around C228T / C250T (GRCh37 or GRCh38)."""
    try:
        p = int(pos)
    except (ValueError, TypeError):
        return False
    chrom_s = str(chrom).lower().lstrip("chr")
    if chrom_s != "5":
        return False
    targets = list(PROMOTER_GRCH37.values()) + list(PROMOTER_GRCH38.values())
    return any(abs(p - t) <= PROMOTER_REGION_BUFFER for t in targets)


TCGA_BARCODE_RE = re.compile(r"TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}(?:-[0-9A-Z]{2,3}[A-Z]?)?", re.I)


def extract_tcga_barcodes(text: str) -> list[str]:
    return sorted({m.group(0).upper() for m in TCGA_BARCODE_RE.finditer(text)})


def write_attempt_log(source_id: str, attempts: list[Attempt]) -> Path:
    p = LOGS / f"{source_id}_attempts.json"
    p.write_text(json.dumps([asdict(a) for a in attempts], indent=2))
    return p


def write_result(result: SourceResult) -> Path:
    p = V2_DIR / f"{result.source_id}_result.json"
    p.write_text(json.dumps(result.to_json(), indent=2, default=str))
    return p


def make_session_kwargs() -> dict:
    return dict(
        connector=aiohttp.TCPConnector(limit=20, ssl=False),
        trust_env=True,
    )
