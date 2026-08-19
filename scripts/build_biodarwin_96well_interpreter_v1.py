#!/usr/bin/env python3
"""Interpret the BioDarwin 96-well router plate with all practical local algorithms.

This is a wetlab-facing scorebook, not a claim generator.
It adds:
  * BigMHC IM / EL on valid class-I rows
  * RF_biophys trained on the master benchmark train split
  * BioDarwin public-anchor v3
  * RF-anchor heuristic and regime router
  * GA public-feature adapter v0
  * BioDarwin mode-bank v4

The output is a single plate-wide table plus a locked-control metric table.
"""
from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PLATE = ROOT / "project/results/biodarwin_96well_router_plate_2026_05_11/biodarwin_96well_router_plate_v1.tsv"
OUT = ROOT / "project/results/biodarwin_96well_interpreter_2026_05_11"
OUT.mkdir(parents=True, exist_ok=True)
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
BIGMHC = Path("/tmp/bigmhc")
WAVE = ROOT / "project/results/p_neo_bayesian_2026_05_09"
GA = ROOT / "project/results/public_industrial_neoantigen_benchmark_scout_2026_05_10/ga_public_feature_adapter_v0"
RF_COMMON = ROOT / "project/results/cancer_vaccine_robustness_2026_05_09"

sys.path.insert(0, str(RF_COMMON))
from _common import biophys, build_hla_onehot_factory, fit_rf, normalize_hla  # noqa: E402


AA = set("ACDEFGHIKLMNPQRSTVWY")
FEATURES = [
    "bigmhc_im_norm",
    "bigmhc_el_norm",
    "peptide_quality_proxy",
    "length_pref",
    "diversity",
    "hydrophobic",
    "aromatic",
    "charged_balance",
    "hla_A",
    "hla_B",
    "hla_C",
    "im_x_quality",
    "el_x_quality",
    "im_x_el",
]
TRANSFORMS = ["identity", "sqrt", "square", "sigmoid"]


def minmax_array(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    finite = np.isfinite(x)
    out = np.zeros_like(x, dtype=float)
    if not finite.any():
        return out
    lo = float(np.min(x[finite]))
    hi = float(np.max(x[finite]))
    if hi == lo:
        out[finite] = 0.5
    else:
        out[finite] = (x[finite] - lo) / (hi - lo)
    return out


def minmax_series(s: pd.Series) -> pd.Series:
    return pd.Series(minmax_array(pd.to_numeric(s, errors="coerce").to_numpy()), index=s.index)


def transform(x: np.ndarray, name: str) -> np.ndarray:
    x = np.nan_to_num(x, nan=0.0, posinf=1.0, neginf=0.0)
    x = np.clip(x, 0.0, 1.0)
    if name == "identity":
        return x
    if name == "sqrt":
        return np.sqrt(x)
    if name == "square":
        return x * x
    if name == "sigmoid":
        return 1.0 / (1.0 + np.exp(-8.0 * (x - 0.5)))
    raise ValueError(name)


def valid_peptide(x: object) -> bool:
    if pd.isna(x):
        return False
    s = str(x).strip().upper()
    return 8 <= len(s) <= 15 and set(s).issubset(AA)


def peptide_features(peptide: object) -> dict[str, float]:
    pep = str(peptide or "").upper()
    pep = "".join(aa for aa in pep if aa in AA)
    if not pep:
        return {
            "seq_cytotoxic_prior": 0.0,
            "seq_helper_prior": 0.0,
            "helper_promiscuity_proxy": 0.0,
            "slp_processability_proxy": 0.0,
            "seq_diversity": 0.0,
            "hydrophobic": 0.0,
            "aromatic": 0.0,
            "charged_balance": 0.0,
            "length_pref": 0.0,
            "diversity": 0.0,
            "peptide_quality_proxy": 0.0,
        }
    n = len(pep)
    diversity = min(1.0, len(set(pep)) / min(12.0, max(1.0, float(n))))
    hydrophobic = sum(aa in "AILMVPFWY" for aa in pep) / n
    aromatic = sum(aa in "FWY" for aa in pep) / n
    charged = sum(aa in "KRDE" for aa in pep) / n
    charged_balance = 1.0 - min(1.0, abs(charged - 0.22) / 0.35)
    class_i_len = 1.0 if n in {9, 10, 11} else 0.72 if n in {8, 12, 13, 14, 15} else 0.20
    class_ii_len = math.exp(-abs(n - 15.0) / 5.0) if 8 <= n <= 35 else 0.0
    helper_promiscuity = 0.38 * diversity + 0.23 * charged_balance + 0.20 * hydrophobic + 0.19 * class_ii_len
    slp_processability = 0.35 * class_ii_len + 0.25 * diversity + 0.20 * (1.0 - abs(hydrophobic - 0.45)) + 0.20 * charged_balance
    cytotoxic = 0.40 * class_i_len + 0.20 * diversity + 0.17 * hydrophobic + 0.13 * aromatic + 0.10 * charged_balance
    helper = 0.32 * class_ii_len + 0.24 * helper_promiscuity + 0.20 * slp_processability + 0.14 * diversity + 0.10 * charged_balance
    quality = 0.42 * class_i_len + 0.22 * diversity + 0.16 * hydrophobic + 0.10 * aromatic + 0.10 * charged_balance
    length_pref = 1.0 if n in {9, 10, 11} else 0.75 if n in {8, 12, 13, 14, 15} else 0.0
    return {
        "seq_cytotoxic_prior": float(np.clip(cytotoxic, 0, 1)),
        "seq_helper_prior": float(np.clip(helper, 0, 1)),
        "helper_promiscuity_proxy": float(np.clip(helper_promiscuity, 0, 1)),
        "slp_processability_proxy": float(np.clip(slp_processability, 0, 1)),
        "seq_diversity": float(np.clip(diversity, 0, 1)),
        "hydrophobic": float(np.clip(hydrophobic, 0, 1)),
        "aromatic": float(np.clip(aromatic, 0, 1)),
        "charged_balance": float(np.clip(charged_balance, 0, 1)),
        "length_pref": float(length_pref),
        "diversity": float(min(1.0, len(set(pep)) / 8.0)),
        "peptide_quality_proxy": float(np.clip(quality, 0, 1)),
    }


def hla_features(hla: pd.Series) -> pd.DataFrame:
    h = hla.astype(str).str.upper().str.replace("HLA-", "", regex=False)
    return pd.DataFrame(
        {
            "hla_A": h.str.startswith("A").astype(float),
            "hla_B": h.str.startswith("B").astype(float),
            "hla_C": h.str.startswith("C").astype(float),
        },
        index=hla.index,
    )


def build_common_features(df: pd.DataFrame, peptide_col: str, hla_col: str) -> pd.DataFrame:
    out = df.copy()
    pf = pd.DataFrame([peptide_features(x) for x in out[peptide_col]], index=out.index)
    hf = hla_features(out[hla_col])
    for col in pf.columns:
        out[col] = pf[col].astype(float)
    for col in hf.columns:
        out[col] = hf[col].astype(float)
    out["im_x_quality"] = 0.0
    out["el_x_quality"] = 0.0
    out["im_x_el"] = 0.0
    return out


def run_bigmhc(bigmhc_dir: Path, input_csv: Path, model: str, out_csv: Path, force: bool = False) -> None:
    if out_csv.exists() and not force:
        return
    script = bigmhc_dir / "src/predict.py"
    if not script.exists():
        raise FileNotFoundError(f"BigMHC predict.py not found: {script}")
    cmd = [
        "python3",
        str(script),
        "-i",
        str(input_csv),
        "-m",
        model,
        "-a",
        "1",
        "-p",
        "2",
        "-c",
        "1",
        "-d",
        "cpu",
        "-b",
        "256",
        "-o",
        str(out_csv),
    ]
    subprocess.run(cmd, check=True, cwd=str(bigmhc_dir))


def add_bigmhc_scores(df: pd.DataFrame, out_dir: Path, force: bool = False) -> pd.DataFrame:
    out = df.copy().reset_index(drop=True)
    eligible = out[out["is_bigmhc_eligible"]].copy()
    eligible["mhc"] = eligible["hla_allele"].astype(str).map(lambda x: str(x).upper().replace("HLA-", ""))
    eligible["pep"] = eligible["peptide"].astype(str).str.upper()
    input_csv = out_dir / "plate_bigmhc_input.csv"
    eligible[["candidate_id", "mhc", "pep", "plate_well", "source_id"]].to_csv(input_csv, index=False)

    im_out = out_dir / "plate_bigmhc_im.csv"
    el_out = out_dir / "plate_bigmhc_el.csv"
    run_bigmhc(BIGMHC, input_csv, "im", im_out, force=force)
    run_bigmhc(BIGMHC, input_csv, "el", el_out, force=force)

    im = pd.read_csv(im_out)
    el = pd.read_csv(el_out)
    im = im.rename(columns={"BigMHC_IM": "BigMHC_IM"})
    el = el.rename(columns={"BigMHC_EL": "BigMHC_EL"})
    merged = eligible[["candidate_id"]].copy()
    for frame, col in [(im, "BigMHC_IM"), (el, "BigMHC_EL")]:
        if "candidate_id" not in frame.columns:
            if len(frame) != len(eligible):
                raise ValueError(f"BigMHC output row mismatch for {col}: {len(frame)} vs {len(eligible)}")
            frame = frame.copy()
            frame.insert(0, "candidate_id", eligible["candidate_id"].to_numpy())
        merged = merged.merge(frame[["candidate_id", col]], on="candidate_id", how="left")
    out = out.merge(merged, on="candidate_id", how="left")
    return out


def score_rf(train_df: pd.DataFrame, eval_df: pd.DataFrame) -> tuple[object, np.ndarray]:
    train_df = train_df.copy()
    eval_df = eval_df.copy()
    train_df["hla_norm"] = train_df["hla_allele_4digit"].apply(normalize_hla)
    eval_df["hla_norm"] = eval_df["hla_allele"].apply(normalize_hla).fillna("UNKNOWN")
    train_df = train_df.dropna(subset=["hla_norm"]).reset_index(drop=True)
    encode, _, _ = build_hla_onehot_factory(train_df["hla_norm"].tolist())

    def feats(df: pd.DataFrame) -> np.ndarray:
        bp = biophys(df["peptide"].tolist())
        hla = encode(df["hla_norm"].tolist())
        return np.hstack([bp, hla]).astype(np.float32)

    Xtr = feats(train_df)
    ytr = train_df["label"].astype(int).to_numpy()
    rf = fit_rf(Xtr, ytr)
    return rf, feats(eval_df)


def transform_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["BigMHC_IM_norm"] = minmax_series(out["BigMHC_IM"]).fillna(0.0)
    out["BigMHC_EL_norm"] = minmax_series(out["BigMHC_EL"]).fillna(0.0)
    out["bigmhc_im_norm"] = out["BigMHC_IM_norm"]
    out["bigmhc_el_norm"] = out["BigMHC_EL_norm"]
    out["rf_biophys_norm"] = minmax_series(out["RF_biophys_score"]).fillna(0.0)
    out["RF_biophys_norm"] = out["rf_biophys_norm"]
    out["rf_bigmhc_gap"] = np.abs(out["rf_biophys_norm"] - out["BigMHC_IM_norm"])
    out["im_x_quality"] = out["BigMHC_IM_norm"] * out["peptide_quality_proxy"]
    out["el_x_quality"] = out["BigMHC_EL_norm"] * out["peptide_quality_proxy"]
    out["im_x_el"] = out["BigMHC_IM_norm"] * out["BigMHC_EL_norm"]
    return out


def add_ga_public_adapter(df: pd.DataFrame) -> pd.DataFrame:
    weights = pd.read_csv(GA / "ga_public_feature_best_weights.tsv", sep="\t")
    gates = pd.read_csv(GA / "ga_public_feature_best_gates.tsv", sep="\t")
    weight_map = {r.feature: (float(r.weight), str(r.transform)) for r in weights.itertuples(index=False)}
    g = {r.gate: float(r.value) for r in gates.itertuples(index=False)}
    out = df.copy()
    raw = np.zeros(len(out), dtype=float)
    for feature, (weight, transform_name) in weight_map.items():
        if feature not in out.columns:
            raise KeyError(f"Missing GA feature: {feature}")
        raw += weight * transform(pd.to_numeric(out[feature], errors="coerce").fillna(0.0).to_numpy(dtype=float), transform_name)
    penalty = np.ones(len(out), dtype=float)
    penalty[out["BigMHC_IM_norm"].to_numpy(dtype=float) < g["im_floor"]] *= 1.0 - g["floor_penalty"]
    penalty[out["BigMHC_EL_norm"].to_numpy(dtype=float) < g["el_floor"]] *= 1.0 - 0.5 * g["floor_penalty"]
    complexity = np.array([1.0 if x else 0.0 for x in out["is_complexity_row"]], dtype=float)
    raw = raw * np.clip(penalty, 0.05, 1.0) - g["complexity_penalty"] * complexity
    out["GA_public_feature_adapter_v0_raw"] = raw
    out["GA_public_feature_adapter_v0"] = minmax_array(raw)
    out["GA_public_feature_adapter_v0_formula"] = "GA frozen weights + floors"
    return out


def add_biodarwin_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["BioDarwin_public_anchor_v3_raw"] = (
        0.95 * out["BigMHC_IM_norm"].to_numpy(dtype=float)
        - 0.05 * out["BigMHC_EL_norm"].to_numpy(dtype=float)
        + 0.10 * out["seq_cytotoxic_prior"].to_numpy(dtype=float)
        + 0.02 * out["seq_helper_prior"].to_numpy(dtype=float)
    )
    out["BioDarwin_public_anchor_v3"] = minmax_array(out["BioDarwin_public_anchor_v3_raw"].to_numpy(dtype=float))
    out["rf_public_anchor_gap"] = np.abs(out["rf_biophys_norm"] - out["BioDarwin_public_anchor_v3"])
    out["BioDarwin_mode_bank_v4"] = minmax_array(
        0.80 * out["BioDarwin_public_anchor_v3"].to_numpy(dtype=float)
        + 0.20 * out["BigMHC_IM_norm"].to_numpy(dtype=float)
    )
    out["BioDarwin_mode_bank_mode"] = "public_anchor"
    out["BioDarwin_RF_anchor_heuristic_v1"] = minmax_array(
        0.45 * out["rf_biophys_norm"].to_numpy(dtype=float)
        + 0.25 * out["BigMHC_IM_norm"].to_numpy(dtype=float)
        + 0.20 * out["BioDarwin_public_anchor_v3"].to_numpy(dtype=float)
        + 0.10 * (1.0 - out["rf_bigmhc_gap"].to_numpy(dtype=float))
    )
    out["BioDarwin_regime_router_v1"] = minmax_array(
        0.70 * out["BioDarwin_public_anchor_v3"].to_numpy(dtype=float)
        + 0.30 * out["BigMHC_IM_norm"].to_numpy(dtype=float)
        - 0.10 * out["rf_bigmhc_gap"].to_numpy(dtype=float)
    )
    out["score_consensus"] = out[[
        "BigMHC_IM_norm",
        "rf_biophys_norm",
        "BioDarwin_public_anchor_v3",
        "GA_public_feature_adapter_v0",
        "BioDarwin_mode_bank_v4",
        "BioDarwin_regime_router_v1",
    ]].mean(axis=1)
    out["score_disagreement"] = out[[
        "BigMHC_IM_norm",
        "rf_biophys_norm",
        "BioDarwin_public_anchor_v3",
        "GA_public_feature_adapter_v0",
        "BioDarwin_mode_bank_v4",
        "BioDarwin_regime_router_v1",
    ]].std(axis=1).fillna(0.0)
    out["topline_vote_count"] = (
        (out["BigMHC_IM_norm"] >= out["BigMHC_IM_norm"].median()).astype(int)
        + (out["rf_biophys_norm"] >= out["rf_biophys_norm"].median()).astype(int)
        + (out["BioDarwin_public_anchor_v3"] >= out["BioDarwin_public_anchor_v3"].median()).astype(int)
        + (out["GA_public_feature_adapter_v0"] >= out["GA_public_feature_adapter_v0"].median()).astype(int)
        + (out["BioDarwin_mode_bank_v4"] >= out["BioDarwin_mode_bank_v4"].median()).astype(int)
        + (out["BioDarwin_regime_router_v1"] >= out["BioDarwin_regime_router_v1"].median()).astype(int)
    )
    return out


def metric_row(df: pd.DataFrame, score_col: str, algorithm: str) -> dict[str, object]:
    d = df[["label", score_col]].copy()
    d["label"] = pd.to_numeric(d["label"], errors="coerce")
    d[score_col] = pd.to_numeric(d[score_col], errors="coerce")
    d = d.dropna()
    row: dict[str, object] = {
        "algorithm": algorithm,
        "score_col": score_col,
        "n": int(len(d)),
        "positives": int(d["label"].sum()),
    }
    if d["label"].nunique() == 2 and len(d) > 0:
        row["AUPRC"] = float(average_precision_score(d["label"], d[score_col]))
        row["AUROC"] = float(roc_auc_score(d["label"], d[score_col]))
    else:
        row["AUPRC"] = np.nan
        row["AUROC"] = np.nan
    ordered = d.sort_values(score_col, ascending=False)
    for k in [5, 10, 24]:
        top = ordered.head(min(k, len(ordered)))
        row[f"top{k}_hits"] = int(top["label"].sum()) if len(top) else 0
        row[f"top{k}_precision"] = float(top["label"].mean()) if len(top) else np.nan
    return row


def assign_interpretation(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["wetlab_priority"] = np.select(
        [
            out["mhc_class"].astype(str).eq("II"),
            out["analysis_role"].eq("locked_metric_control") & out["label"].eq(1),
            out["analysis_role"].eq("locked_metric_control") & out["label"].eq(0),
            out["analysis_role"].eq("negative_control"),
            out["score_consensus"] >= 0.72,
            out["score_consensus"] <= 0.35,
        ],
        [
            "helper_lane",
            "hit_reference",
            "fail_reference",
            "fail_reference",
            "hit_candidate",
            "fail_candidate",
        ],
        default="probe_candidate",
    )
    out["wetlab_reason"] = np.select(
        [
            out["mhc_class"].astype(str).eq("II"),
            out["analysis_role"].eq("locked_metric_control") & out["label"].eq(1),
            out["analysis_role"].eq("locked_metric_control") & out["label"].eq(0),
            out["analysis_role"].eq("negative_control"),
            out["score_disagreement"] >= 0.18,
            out["score_consensus"] >= 0.72,
            out["score_consensus"] <= 0.35,
        ],
        [
            "class-II helper lane",
            "known positive control",
            "known negative control",
            "designed fail control",
            "router disagreement probe",
            "consensus hit candidate",
            "consensus fail candidate",
        ],
        default="balanced probe",
    )
    return out


def safe_metrics(df: pd.DataFrame, score_col: str) -> dict[str, object]:
    d = df[df["analysis_role"].eq("locked_metric_control")][["label", score_col]].copy()
    d["label"] = pd.to_numeric(d["label"], errors="coerce")
    d[score_col] = pd.to_numeric(d[score_col], errors="coerce")
    d = d.dropna()
    out = {
        "score_col": score_col,
        "n": int(len(d)),
        "positives": int(d["label"].sum()) if len(d) else 0,
    }
    if d["label"].nunique() == 2 and len(d) > 0:
        out["AUPRC"] = float(average_precision_score(d["label"], d[score_col]))
        out["AUROC"] = float(roc_auc_score(d["label"], d[score_col]))
    else:
        out["AUPRC"] = np.nan
        out["AUROC"] = np.nan
    ordered = d.sort_values(score_col, ascending=False)
    for k in [5, 10, 24]:
        top = ordered.head(min(k, len(ordered)))
        out[f"top{k}_hits"] = int(top["label"].sum()) if len(top) else 0
        out[f"top{k}_precision"] = float(top["label"].mean()) if len(top) else np.nan
    return out


def table_html(df: pd.DataFrame, max_rows: int = 40) -> str:
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_html(metrics: pd.DataFrame, plate: pd.DataFrame, summary: dict[str, object]) -> None:
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BioDarwin 96-well interpreter v1</title>
<style>
body {{ margin:0; background:#0f1418; color:#e8edf0; font:14px/1.55 system-ui,-apple-system,Segoe UI,sans-serif; }}
header {{ padding:34px 5vw 20px; background:#0c1014; border-bottom:1px solid #27313a; }}
h1 {{ margin:.2rem 0 .5rem; font-size:34px; line-height:1.06; }}
.lead {{ color:#cfd8dd; max-width:1100px; }}
main {{ padding:24px 5vw 56px; }}
.note {{ background:#171e24; border-left:3px solid #c59b3b; padding:10px 14px; margin:0 0 18px; }}
table.data {{ width:100%; border-collapse:collapse; font-size:13px; margin:10px 0 28px; }}
table.data th, table.data td {{ border-bottom:1px solid #2a3640; padding:7px 8px; text-align:left; vertical-align:top; }}
table.data th {{ color:#f1d58b; background:#141b21; position:sticky; top:0; }}
code {{ color:#f1d58b; }}
</style>
</head>
<body>
<header>
  <div style="color:#c59b3b;text-transform:uppercase;letter-spacing:.12em;font-weight:700;font-size:12px;">BioDarwin 96-well interpreter</div>
  <h1>All practical local algorithms on one wetlab plate</h1>
  <p class="lead">BigMHC IM/EL, RF_biophys, GA public adapter, BioDarwin public-anchor v3, mode bank v4, RF-anchor heuristic, and regime router v1 are scored on the same 96-well board.</p>
</header>
<main>
  <div class="note">Locked-control metrics are computed only on the 25 labeled rows. Class-II and invalid-peptide rows are carried through as helper/probe lanes, and BigMHC is only called on valid class-I peptides.</div>
  <h2>Summary</h2>
  <pre>{json.dumps(summary, indent=2)}</pre>
  <h2>Locked-control metrics</h2>
  {table_html(metrics, 40)}
  <h2>Top plate rows by consensus</h2>
  {table_html(plate.sort_values("score_consensus", ascending=False)[[
    "plate_well","source_id","peptide","hla_allele","mhc_class","analysis_role","wetlab_priority","wetlab_reason",
    "BigMHC_IM","BigMHC_EL","RF_biophys_score","GA_public_feature_adapter_v0","BioDarwin_public_anchor_v3",
    "BioDarwin_mode_bank_v4","BioDarwin_regime_router_v1","score_consensus","score_disagreement","topline_vote_count"
  ]], 40)}
</main>
</body>
</html>
"""
    path = HUB / "biodarwin_96well_interpreter_v1.html"
    path.write_text(html, encoding="utf-8")
    try:
        live = LIVE_HUB / path.name
        shutil.copy2(path, live)
    except Exception as exc:
        print(f"deploy warning: {exc}")


def main() -> None:
    t0 = time.time()
    plate = pd.read_csv(PLATE, sep="\t")
    plate = plate.copy().reset_index(drop=True)
    plate["candidate_id"] = plate["plate_well"].astype(str)
    plate["label"] = pd.to_numeric(plate["response_label"], errors="coerce")
    plate["hla_norm"] = plate["hla_allele"].astype(str).map(normalize_hla)
    plate["is_bigmhc_eligible"] = (
        plate["mhc_class"].astype(str).str.upper().eq("I")
        & plate["hla_norm"].notna()
        & plate["peptide"].map(valid_peptide)
    )
    plate["is_complexity_row"] = plate["analysis_role"].astype(str).eq("exploratory_probe") | plate["analysis_role"].astype(str).eq("negative_control")
    plate = build_common_features(plate, "peptide", "hla_allele")

    # BigMHC on valid class-I rows only.
    plate = add_bigmhc_scores(plate, OUT, force=False)

    # RF_biophys trained on the master train split.
    train_pool = pd.read_csv(WAVE / "bundle.tsv", sep="\t")
    train_pool = train_pool[train_pool["split"].eq("train")].dropna(subset=["HLA_norm"]).copy()
    train_pool["hla_allele_4digit"] = train_pool["HLA_norm"]
    train_pool["peptide"] = train_pool["peptide"].astype(str).str.upper()
    rf, X_plate = score_rf(train_pool, plate)
    plate["RF_biophys_score"] = rf.predict_proba(X_plate)[:, 1]

    # Normalize raw scores and build public-anchored layers.
    plate = transform_features(plate)
    plate = add_ga_public_adapter(plate)
    plate = add_biodarwin_scores(plate)
    plate = assign_interpretation(plate)

    # Metrics only on labeled locked-control rows.
    metrics_rows = []
    for algo, col in [
        ("BigMHC_IM", "BigMHC_IM"),
        ("BigMHC_EL", "BigMHC_EL"),
        ("RF_biophys", "RF_biophys_norm"),
        ("GA_public_feature_adapter_v0", "GA_public_feature_adapter_v0"),
        ("BioDarwin_public_anchor_v3", "BioDarwin_public_anchor_v3"),
        ("BioDarwin_mode_bank_v4", "BioDarwin_mode_bank_v4"),
        ("BioDarwin_RF_anchor_heuristic_v1", "BioDarwin_RF_anchor_heuristic_v1"),
        ("BioDarwin_regime_router_v1", "BioDarwin_regime_router_v1"),
    ]:
        metrics_rows.append({"algorithm": algo, **safe_metrics(plate, col)})
    metrics = pd.DataFrame(metrics_rows).sort_values(["AUPRC", "AUROC"], ascending=False)

    # Compact outputs.
    ordered_cols = [
        "plate_well",
        "candidate_id",
        "source_id",
        "line_number",
        "peptide",
        "hla_allele",
        "mhc_class",
        "candidate_label_tier",
        "response_label",
        "label",
        "label_raw",
        "plate_bucket",
        "analysis_role",
        "hit_fail_interpretation",
        "wetlab_priority",
        "wetlab_reason",
        "BigMHC_IM",
        "BigMHC_EL",
        "RF_biophys_score",
        "GA_public_feature_adapter_v0",
        "BioDarwin_public_anchor_v3",
        "BioDarwin_mode_bank_v4",
        "BioDarwin_RF_anchor_heuristic_v1",
        "BioDarwin_regime_router_v1",
        "score_consensus",
        "score_disagreement",
        "topline_vote_count",
    ]
    plate_out = plate[ordered_cols].copy()

    plate_out.to_csv(OUT / "biodarwin_96well_interpreter_v1.tsv", sep="\t", index=False)
    metrics.to_csv(OUT / "biodarwin_96well_interpreter_v1_metrics.tsv", sep="\t", index=False)

    summary = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "plate_rows": int(len(plate_out)),
        "bigmhc_scored_rows": int(plate["BigMHC_IM"].notna().sum()),
        "bigmhc_eligible_rows": int(plate["is_bigmhc_eligible"].sum()),
        "locked_metric_rows": int((plate["analysis_role"] == "locked_metric_control").sum()),
        "positive_locked_rows": int((plate["analysis_role"].eq("locked_metric_control") & plate["label"].eq(1)).sum()),
        "negative_locked_rows": int((plate["analysis_role"].eq("locked_metric_control") & plate["label"].eq(0)).sum()),
        "helper_lane_rows": int((plate["wetlab_priority"] == "helper_lane").sum()),
        "hit_candidates": int((plate["wetlab_priority"] == "hit_candidate").sum()),
        "fail_candidates": int((plate["wetlab_priority"] == "fail_candidate").sum()),
        "probe_candidates": int((plate["wetlab_priority"] == "probe_candidate").sum()),
        "best_locked_control": metrics.head(1).to_dict(orient="records"),
        "outputs": {
            "plate": str(OUT / "biodarwin_96well_interpreter_v1.tsv"),
            "metrics": str(OUT / "biodarwin_96well_interpreter_v1_metrics.tsv"),
            "html": str(HUB / "biodarwin_96well_interpreter_v1.html"),
        },
        "elapsed_s": round(time.time() - t0, 2),
    }
    (OUT / "biodarwin_96well_interpreter_v1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    report = [
        "# BioDarwin 96-well interpreter v1",
        "",
        f"Generated: {summary['generated_at']}",
        "",
        "## Locked-control metrics",
        metrics.to_markdown(index=False),
        "",
        "## Summary",
        pd.DataFrame([summary]).to_markdown(index=False),
        "",
        "## Decision rule",
        "- `hit_candidate`: consensus score >= 0.72",
        "- `fail_candidate`: consensus score <= 0.35",
        "- `helper_lane`: class-II rows",
        "- `probe_candidate`: everything else",
    ]
    (OUT / "BIODARWIN_96WELL_INTERPRETER_V1.md").write_text("\n".join(report), encoding="utf-8")
    write_html(metrics, plate_out, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
