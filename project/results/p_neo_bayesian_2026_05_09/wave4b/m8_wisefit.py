"""M8 — WiSE-FT (Wortsman 2022, CVPR) — STATUS: SKIPPED (Wave 4A not ready).

Brief permission: "If Wave 4A LoRA isn't ready by your deadline, use Wave 1
head only (skip WiSE-FT, document in notes)."

State:
  - wave4a/ directory empty as of run-time (only timestamp, no checkpoints).
  - wave1/ directory contains predictions only — no Wave 1 head state-dict
    on this filesystem (the head was trained on the pod and not pulled back).
  - We therefore have NEITHER endpoint for the θ_blend = α·θ_finetune + (1-α)·θ_frozen
    interpolation, and skip M8 cleanly.

This file documents the decision and emits a placeholder predictions_wisefit.tsv
copy of Wave 1's predictions_itsndb.tsv tagged with α=0 (== frozen baseline)
so the wave4b_results.tsv schema stays uniform.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
WAVE1 = ROOT.parent
OUT = ROOT / "predictions_wisefit.tsv"

src = WAVE1 / "predictions_itsndb.tsv"
df = pd.read_csv(src, sep="\t")
df["alpha"] = 0.0  # WiSE-FT α=0 → frozen baseline
df["pred_wisefit"] = df["pred_mean"]
df.to_csv(OUT, sep="\t", index=False)
print(f"M8 SKIPPED — Wave 4A not ready, Wave 1 head dict unavailable. "
      f"Wrote α=0 placeholder to {OUT}.")
