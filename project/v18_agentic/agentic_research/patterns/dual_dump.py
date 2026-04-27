"""Pattern 3 -- Dual Dump (raw + compressed).

For every audit-eligible artifact, write two files:
  * raw      -- the full, machine-readable payload (json/tsv/log) for replay
  * compressed -- a markdown summary the user can read in 60 seconds

This pattern was the single most useful safety net in v3 -> v17p35: when a
phase looked successful but the audit failed, the compressed dump told us
which assumption broke without re-running anything.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class DumpPaths:
    raw: Path
    compressed: Path


class DualDumpRecorder:
    def __init__(self, root: str | Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def dump(
        self,
        name: str,
        raw_payload: Any,
        compressed_md: str,
        raw_format: str = "json",
    ) -> DumpPaths:
        raw_path = self.root / f"{name}.{raw_format}"
        if raw_format == "json":
            raw_path.write_text(json.dumps(raw_payload, indent=2, default=str))
        elif raw_format in {"tsv", "txt", "log"}:
            raw_path.write_text(str(raw_payload))
        else:
            raise ValueError(f"unsupported raw_format: {raw_format}")
        compressed_path = self.root / f"{name}.md"
        compressed_path.write_text(compressed_md, encoding="utf-8")
        return DumpPaths(raw=raw_path, compressed=compressed_path)
