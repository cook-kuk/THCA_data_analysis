#!/usr/bin/env python3
"""Merge split RunPod WSI embeddings and recompute 50-WSI DM1/DM2 LOSO."""

from __future__ import annotations

import json
import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut


ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = Path(
    os.environ.get(
        "RUNPOD_WSI_RUN_DIR",
        ROOT / "project/data/processed/TCGA-THCA-WSI-DM/runpod_split_2026_05_09",
    )
)
MANIFEST = ROOT / "project/data/manifests/tcga_thca_wsi_dm_balanced.tsv"
OUT = RUN_DIR / "merged"


def load_embeddings() -> tuple[dict[str, np.ndarray], dict[str, str], list[dict[str, object]]]:
    merged: dict[str, np.ndarray] = {}
    backbones: dict[str, str] = {}
    rows: list[dict[str, object]] = []
    for pod_dir in sorted(RUN_DIR.glob("pod_*")):
        pkl = pod_dir / "embeddings/slide_embeddings.pkl"
        if not pkl.exists():
            continue
        with pkl.open("rb") as f:
            pkg = pickle.load(f)
        emb = pkg.get("embeddings", {})
        backbone = pkg.get("backbone", "unknown")
        backbones[pod_dir.name] = backbone
        for file_id, vec in emb.items():
            merged[file_id] = np.asarray(vec)
            rows.append(
                {
                    "pod": pod_dir.name,
                    "file_id": file_id,
                    "embedding_dim": int(np.asarray(vec).shape[0]),
                    "backbone": backbone,
                }
            )
    return merged, backbones, rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    embeddings, backbones, emb_rows = load_embeddings()

    resolved = []
    for pod_dir in sorted(RUN_DIR.glob("pod_*")):
        path = pod_dir / "wsi_resolved_manifest.tsv"
        if path.exists():
            df = pd.read_csv(path, sep="\t")
            df["pod"] = pod_dir.name
            resolved.append(df)
    resolved_df = pd.concat(resolved, ignore_index=True).drop_duplicates("file_id")
    manifest = pd.read_csv(MANIFEST, sep="\t")
    meta = resolved_df.merge(manifest, on="case_id", how="left")

    X, y, groups, out_rows = [], [], [], []
    for _, row in meta.iterrows():
        file_id = row["file_id"]
        dm = row.get("dm")
        if file_id not in embeddings or dm not in {"DM1", "DM2"}:
            continue
        X.append(embeddings[file_id])
        y.append(0 if dm == "DM1" else 1)
        groups.append(row["case_id"])
        out_rows.append(row.to_dict())

    X_arr = np.asarray(X)
    y_arr = np.asarray(y)
    groups_arr = np.asarray(groups)
    preds = np.full(len(y_arr), np.nan)

    for train_idx, test_idx in LeaveOneGroupOut().split(X_arr, y_arr, groups_arr):
        if len(np.unique(y_arr[train_idx])) < 2:
            continue
        class_counts = np.bincount(y_arr[train_idx])
        inner_cv = int(min(3, class_counts[class_counts > 0].min()))
        if inner_cv < 2:
            inner_cv = 2
        model = LogisticRegressionCV(
            Cs=[0.01, 0.1, 1, 10],
            cv=inner_cv,
            max_iter=5000,
        ).fit(X_arr[train_idx], y_arr[train_idx])
        preds[test_idx] = model.predict_proba(X_arr[test_idx])[:, 1]

    ok = ~np.isnan(preds)
    auc = float(roc_auc_score(y_arr[ok], preds[ok])) if len(np.unique(y_arr[ok])) == 2 else None
    acc = float(accuracy_score(y_arr[ok], (preds[ok] > 0.5).astype(int))) if ok.any() else None

    pred_df = pd.DataFrame(out_rows)
    pred_df["dm2_pred_prob"] = preds
    pred_df["dm_true"] = np.where(y_arr == 0, "DM1", "DM2")
    pred_df.to_csv(OUT / "merged_dm1_vs_dm2_loso_predictions.tsv", sep="\t", index=False)

    emb_meta = pd.DataFrame(emb_rows)
    emb_meta.to_csv(OUT / "embedding_coverage.tsv", sep="\t", index=False)

    missing = meta.loc[~meta["file_id"].isin(embeddings), ["pod", "file_id", "case_id", "dm", "file_name"]]
    missing.to_csv(OUT / "missing_embedding_slides.tsv", sep="\t", index=False)

    with (OUT / "merged_slide_embeddings.pkl").open("wb") as f:
        pickle.dump({"backbones": backbones, "embeddings": embeddings}, f)

    metrics = {
        "n_manifest_slides": int(len(meta)),
        "n_embedded_slides": int(len(pred_df)),
        "n_missing_embeddings": int(len(missing)),
        "n_cases": int(pd.Series(groups_arr).nunique()),
        "n_DM1": int((y_arr == 0).sum()),
        "n_DM2": int((y_arr == 1).sum()),
        "auc_dm2_probability": auc,
        "accuracy_at_0_5": acc,
        "backbones": backbones,
        "note": "DINOv2/UNI split RunPod embeddings merged locally; LOSO grouped by case_id.",
    }
    (OUT / "merged_loso_metrics.json").write_text(json.dumps(metrics, indent=2))

    summary = [
        "# RunPod split WSI merge summary",
        "",
        f"- Manifest slides: {metrics['n_manifest_slides']}",
        f"- Embedded slides used: {metrics['n_embedded_slides']}",
        f"- Missing embeddings: {metrics['n_missing_embeddings']}",
        f"- Class counts: DM1={metrics['n_DM1']}, DM2={metrics['n_DM2']}",
        f"- Merged LOSO AUC (DM2 probability): {auc:.3f}" if auc is not None else "- Merged LOSO AUC: NA",
        f"- Accuracy at 0.5: {acc:.3f}" if acc is not None else "- Accuracy at 0.5: NA",
        "",
        "Outputs:",
        "- `merged_dm1_vs_dm2_loso_predictions.tsv`",
        "- `merged_loso_metrics.json`",
        "- `embedding_coverage.tsv`",
        "- `missing_embedding_slides.tsv`",
    ]
    (OUT / "SUMMARY.md").write_text("\n".join(summary) + "\n")

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
