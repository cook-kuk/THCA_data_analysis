#!/usr/bin/env python3
"""CROSS-Neo 2.0 shared utilities."""

from __future__ import annotations

import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


SEED = 20260509
REPO = Path(__file__).resolve().parents[3]
V0 = REPO / "project/results/cross_neo_v0"
V1 = REPO / "project/results/cross_neo_v1_lockdown"
OUT = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
SCRIPTS = REPO / "project/scripts/cross_neo_v2"

AA = "ACDEFGHIKLMNPQRSTVWY"
KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"E":-3.5,"Q":-3.5,"G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
CHG = {a: 0.0 for a in AA} | {"R":1.0,"K":1.0,"H":0.1,"D":-1.0,"E":-1.0}
HEL = {"A":1.42,"R":0.98,"N":0.67,"D":1.01,"C":0.70,"E":1.51,"Q":1.11,"G":0.57,"H":1.00,"I":1.08,"L":1.21,"K":1.16,"M":1.45,"F":1.13,"P":0.57,"S":0.77,"T":0.83,"W":1.08,"Y":0.69,"V":1.06}


def ensure_dirs() -> None:
    for p in [OUT, OUT / "features", OUT / "splits", OUT / "metrics", OUT / "predictions", OUT / "figures", OUT / "output", OUT / "corpora", OUT / "structure_jobs", OUT / "logs"]:
        p.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def safe_parquet(df: pd.DataFrame, path: Path) -> None:
    try:
        df.to_parquet(path, index=False)
    except Exception:
        df.to_pickle(path.with_suffix(path.suffix + ".pkl"))
        path.with_suffix(path.suffix + ".UNAVAILABLE.txt").write_text("parquet engine unavailable; see .pkl fallback\n")


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet" and path.exists():
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t")


def normalize_hla(hla: object) -> str:
    h = str(hla or "").strip().upper().replace("HLA-", "")
    if h in {"", "NAN", "NONE"}:
        return ""
    h = h.replace("_", "*", 1) if "*" not in h and len(h) >= 5 else h
    if "*" not in h and h[0] in "ABC" and len(h) >= 5:
        h = f"{h[0]}*{h[1:]}"
    if ":" not in h and "*" in h:
        locus, rest = h.split("*", 1)
        rest = rest.replace(":", "")
        if len(rest) >= 4:
            h = f"{locus}*{rest[:2]}:{rest[2:4]}"
    return "HLA-" + h if h and not h.startswith("HLA-") else h


def hla_gene(hla: object) -> str:
    h = normalize_hla(hla)
    return h.split("*", 1)[0].replace("HLA-", "") if "*" in h else ""


def hla_supertype(hla: object) -> str:
    h = normalize_hla(hla)
    if "*" not in h:
        return ""
    gene = hla_gene(h)
    fam = h.split("*", 1)[1].split(":", 1)[0]
    return f"{gene}{fam}"


def kmer_set(seq: object, k: int) -> set[str]:
    s = str(seq or "")
    if len(s) < k:
        return {s} if s else set()
    return {s[i:i+k] for i in range(len(s)-k+1)}


def seq_similarity(a: object, b: object) -> float:
    a, b = str(a or ""), str(b or "")
    if not a or not b:
        return 0.0
    ident = sum(x == y for x, y in zip(a, b)) / max(len(a), len(b))
    vals = [ident]
    for k in (3, 4, 5):
        ka, kb = kmer_set(a, k), kmer_set(b, k)
        vals.append(0.0 if not ka or not kb else len(ka & kb) / len(ka | kb))
    return float(np.mean(vals))


def levenshtein(a: object, b: object, max_cutoff: int = 3) -> int:
    a, b = str(a or ""), str(b or "")
    if abs(len(a) - len(b)) > max_cutoff:
        return max_cutoff + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        row_min = i
        for j, cb in enumerate(b, 1):
            val = min(prev[j] + 1, cur[j-1] + 1, prev[j-1] + (ca != cb))
            cur.append(val)
            row_min = min(row_min, val)
        if row_min > max_cutoff:
            return max_cutoff + 1
        prev = cur
    return int(prev[-1])


def peptide_cluster(seq: object) -> str:
    s = str(seq or "")
    if len(s) < 4:
        return f"L{len(s)}_{s}"
    return f"L{len(s)}_{s[:2]}_{s[-2:]}"


def stable_hash(text: str, mod: int) -> int:
    return int(hashlib.md5(text.encode()).hexdigest()[:8], 16) % mod


def aa_physchem(seq: object, prefix: str) -> dict[str, float]:
    s = str(seq or "").upper()
    out = {f"{prefix}_len": float(len(s)), f"{prefix}_missing": float(not bool(s))}
    for a in AA:
        out[f"{prefix}_frac_{a}"] = s.count(a) / max(1, len(s))
    for name, tab in [("kd", KD), ("chg", CHG), ("hel", HEL)]:
        arr = np.array([tab.get(a, 0.0) for a in s], dtype=float)
        out[f"{prefix}_{name}_mean"] = float(arr.mean()) if len(arr) else 0.0
        out[f"{prefix}_{name}_std"] = float(arr.std()) if len(arr) else 0.0
        out[f"{prefix}_{name}_min"] = float(arr.min()) if len(arr) else 0.0
        out[f"{prefix}_{name}_max"] = float(arr.max()) if len(arr) else 0.0
    bins = np.zeros(64, dtype=float)
    for k in (2, 3, 4):
        for km in kmer_set(s, k):
            bins[stable_hash(f"{k}:{km}", len(bins))] += 1
    bins = bins / max(1.0, bins.sum())
    for i, v in enumerate(bins):
        out[f"{prefix}_hash_{i:02d}"] = float(v)
    return out


def metrics(y: np.ndarray, score: np.ndarray, ks=(5, 10, 20)) -> dict[str, float]:
    y = np.asarray(y, dtype=int)
    s = np.asarray(score, dtype=float)
    mask = np.isfinite(s)
    y, s = y[mask], s[mask]
    out = {"n": int(len(y)), "n_pos": int(y.sum()) if len(y) else 0, "prevalence": float(y.mean()) if len(y) else np.nan}
    if len(y) < 2 or len(np.unique(y)) < 2:
        out |= {"AUPRC": np.nan, "AUROC": np.nan, "Brier": np.nan, "ECE": np.nan}
    else:
        out |= {"AUPRC": float(average_precision_score(y, s)), "AUROC": float(roc_auc_score(y, s)), "Brier": float(brier_score_loss(y, np.clip(s, 0, 1))), "ECE": ece(y, s)}
    order = np.argsort(-s)
    for k in ks:
        kk = min(k, len(y))
        top = y[order[:kk]] if kk else np.array([])
        prec = float(top.mean()) if kk else np.nan
        out[f"top{k}_precision"] = prec
        out[f"recall_at_{k}"] = float(top.sum() / max(1, y.sum())) if kk else np.nan
        out[f"enrichment_at_{k}"] = float(prec / out["prevalence"]) if kk and out["prevalence"] > 0 else np.nan
    return out


def ece(y: np.ndarray, score: np.ndarray, n_bins: int = 10) -> float:
    y = np.asarray(y, dtype=int)
    s = np.clip(np.asarray(score, dtype=float), 0, 1)
    val = 0.0
    for lo, hi in zip(np.linspace(0, 1, n_bins+1)[:-1], np.linspace(0, 1, n_bins+1)[1:]):
        m = (s >= lo) & (s < hi if hi < 1 else s <= hi)
        if m.sum():
            val += float(m.mean() * abs(s[m].mean() - y[m].mean()))
    return val


def env_report() -> dict[str, object]:
    return {"python": platform.python_version(), "platform": platform.platform(), "seed": SEED, "repo": str(REPO)}
