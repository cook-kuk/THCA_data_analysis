"""Orchestrator — run all 9 TERT recovery sources concurrently."""
from __future__ import annotations

import asyncio
import json
import time
from dataclasses import asdict
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))

from v17_tert_v2_common import LOGS, V2_DIR, configure_logger

LOG = configure_logger("orchestrator")

SOURCE_MODULES = [
    "v17_tert_v2_S1_liu2017",
    "v17_tert_v2_S2_landa2016",
    "v17_tert_v2_S3_yoo2019",
    "v17_tert_v2_S4_pozdeyev2024",
    "v17_tert_v2_S5_cosmic",
    "v17_tert_v2_S6_cbioportal_all",
    "v17_tert_v2_S7_gdc_controlled",
    "v17_tert_v2_S8_preprint",
    "v17_tert_v2_S9_publication_mining",
]


async def run_one(module_name: str, timeout_s: int = 1800) -> dict:
    """Import the module and call its async run()."""
    t0 = time.time()
    try:
        mod = __import__(module_name)
        result = await asyncio.wait_for(mod.run(), timeout=timeout_s)
        return {
            "source": module_name,
            "status": "ok",
            "elapsed": time.time() - t0,
            "result": result.to_json(),
        }
    except asyncio.TimeoutError:
        return {"source": module_name, "status": "timeout", "elapsed": time.time() - t0}
    except Exception as e:  # noqa: BLE001
        import traceback

        return {
            "source": module_name,
            "status": "error",
            "error": f"{type(e).__name__}: {e}",
            "trace": traceback.format_exc(),
            "elapsed": time.time() - t0,
        }


async def main() -> None:
    LOG.info("Starting %d sources", len(SOURCE_MODULES))
    t0 = time.time()
    tasks = [run_one(m) for m in SOURCE_MODULES]
    results = await asyncio.gather(*tasks, return_exceptions=False)
    elapsed = time.time() - t0
    LOG.info("All sources finished in %.1fs", elapsed)

    summary = {
        "wall_seconds": elapsed,
        "n_sources": len(SOURCE_MODULES),
        "n_ok": sum(1 for r in results if r["status"] == "ok"),
        "n_error": sum(1 for r in results if r["status"] == "error"),
        "n_timeout": sum(1 for r in results if r["status"] == "timeout"),
        "results": results,
    }
    out = V2_DIR / "orchestrator_summary.json"
    out.write_text(json.dumps(summary, indent=2, default=str))
    LOG.info("Summary written to %s", out)

    # Quick stdout summary
    print("\n=== ORCHESTRATOR SUMMARY ===")
    for r in results:
        src = r["source"].replace("v17_tert_v2_", "")
        if r["status"] == "ok":
            res = r.get("result", {})
            print(
                f"  {src}: OK  records={res.get('n_records', 0)} "
                f"promoter={res.get('n_promoter_mutations', 0)} "
                f"tcga_match={res.get('n_tcga_matched', 0)} "
                f"({r['elapsed']:.1f}s)"
            )
        else:
            print(f"  {src}: {r['status'].upper()} {r.get('error', '')} ({r['elapsed']:.1f}s)")
    print("============================")


if __name__ == "__main__":
    asyncio.run(main())
