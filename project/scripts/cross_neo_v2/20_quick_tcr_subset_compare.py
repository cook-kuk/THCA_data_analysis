#!/usr/bin/env python3
"""Quick subset comparisons for the CROSS-Neo-TCR extension.

This is a diagnostic pilot, not a final TCR-aware benchmark. It answers two
near-term questions quickly:

1. On the main CROSS-Neo registry, do local TCR evidence features add signal
   over existing pMHC features on a small set of strict splits?
2. On paired public TCR-pMHC rows, can simple sequence/cross features separate
   cognate examples from shuffled TCR decoys?
"""

from __future__ import annotations

import argparse
import math
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_extraction import DictVectorizer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from common import AA, OUT, SEED, hla_supertype, kmer_set, peptide_cluster, safe_parquet, seq_similarity

warnings.filterwarnings("ignore")


TCR_OUT = OUT / "tcr_extension"
PILOT_OUT = TCR_OUT / "pilot_subset_compare"


DEFAULT_SPLITS = [
    "repeated_stratified_5x5_internal",
    "exact_peptide_hla_holdout",
    "near_peptide_cluster_holdout",
    "source_heldout_NEPdb",
]


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    view = df.copy()
    view = view.where(view.notna(), "NA").replace("", "NA")
    view.to_csv(path, sep="\t", index=False, na_rep="NA")


def read_feature(path: Path) -> pd.DataFrame:
    if path.exists():
        try:
            return pd.read_parquet(path)
        except Exception:
            pass
    pkl = path.with_suffix(path.suffix + ".pkl")
    if pkl.exists():
        return pd.read_pickle(pkl)
    tsv = path.with_suffix(".tsv")
    if tsv.exists():
        return pd.read_csv(tsv, sep="\t", low_memory=False)
    raise FileNotFoundError(path)


def numeric_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        if c != "row_id":
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.fillna(0)


def merge_features(frames: list[pd.DataFrame]) -> pd.DataFrame:
    x = frames[0].copy()
    for f in frames[1:]:
        x = x.merge(f, on="row_id", how="left", suffixes=("", "_dup"))
    dup = [c for c in x.columns if c.endswith("_dup")]
    if dup:
        x = x.drop(columns=dup)
    return numeric_features(x)


def tcr_evidence_features(linked: pd.DataFrame, include_public_label_counts: bool = True) -> pd.DataFrame:
    cols = [
        "tcr_evidence_count",
        "paired_tcr_evidence_count",
        "beta_only_tcr_evidence_count",
        "cancer_context_evidence_count",
        "pathogen_context_evidence_count",
        "structure_evidence_count",
    ]
    if include_public_label_counts:
        cols.extend([
            "peptide_hla_tcr_label_count",
            "tcr_binding_positive_count",
            "tcr_binding_negative_count",
        ])
    out = linked[["row_id"]].copy()
    for col in cols:
        vals = pd.to_numeric(linked.get(col, 0), errors="coerce").fillna(0)
        out[col] = vals
        out[f"log1p_{col}"] = np.log1p(vals)
    out["has_any_tcr_evidence"] = (out["tcr_evidence_count"] > 0).astype(int)
    out["has_paired_tcr_evidence"] = (out["paired_tcr_evidence_count"] > 0).astype(int)
    out["has_structure_tcr_evidence"] = (out["structure_evidence_count"] > 0).astype(int)
    if include_public_label_counts:
        out["tcr_positive_minus_negative"] = out["tcr_binding_positive_count"] - out["tcr_binding_negative_count"]
        out["tcr_positive_fraction"] = out["tcr_binding_positive_count"] / (
            out["tcr_binding_positive_count"] + out["tcr_binding_negative_count"] + 1e-6
        )
    cats = pd.get_dummies(linked["tcr_link_category"].fillna("missing"), prefix="tcrcat", dtype=int)
    return pd.concat([out, cats], axis=1)


def summarize_tcr_hits(hits: pd.DataFrame) -> dict[str, float]:
    if hits.empty:
        return {
            "tcr_evidence_count": 0.0,
            "paired_tcr_evidence_count": 0.0,
            "beta_only_tcr_evidence_count": 0.0,
            "peptide_hla_tcr_label_count": 0.0,
            "cancer_context_evidence_count": 0.0,
            "pathogen_context_evidence_count": 0.0,
            "tcr_binding_positive_count": 0.0,
            "tcr_binding_negative_count": 0.0,
            "structure_evidence_count": 0.0,
        }
    labels = pd.to_numeric(hits.get("binding_label_binary", pd.Series([], dtype=float)), errors="coerce")
    return {
        "tcr_evidence_count": float(len(hits)),
        "paired_tcr_evidence_count": float(hits.get("paired_tcr_available", False).fillna(False).astype(bool).sum()),
        "beta_only_tcr_evidence_count": float(hits.get("beta_only_tcr_available", False).fillna(False).astype(bool).sum()),
        "peptide_hla_tcr_label_count": float((hits.get("peptide_hla_available", False).fillna(False).astype(bool) & labels.notna()).sum()),
        "cancer_context_evidence_count": float(hits.get("is_cancer_context", False).fillna(False).astype(bool).sum()),
        "pathogen_context_evidence_count": float(hits.get("is_pathogen_context", False).fillna(False).astype(bool).sum()),
        "tcr_binding_positive_count": float((labels == 1).sum()),
        "tcr_binding_negative_count": float((labels == 0).sum()),
        "structure_evidence_count": float(hits.get("structure_pdb_id", "").fillna("").astype(str).ne("").sum()),
    }


def source_excluded_tcr_evidence_features(
    neo: pd.DataFrame,
    include_public_label_counts: bool = True,
) -> pd.DataFrame:
    """Recompute TCR evidence after excluding rows from the same source dataset."""
    tcr_cols = [
        "source_dataset",
        "peptide",
        "hla_4digit",
        "hla_supertype",
        "paired_tcr_available",
        "beta_only_tcr_available",
        "peptide_hla_available",
        "binding_label_binary",
        "is_cancer_context",
        "is_pathogen_context",
        "structure_pdb_id",
    ]
    tcr = pd.read_parquet(TCR_OUT / "tcr_registry.parquet", columns=tcr_cols).copy()
    tcr["peptide"] = tcr["peptide"].fillna("").astype(str)
    tcr["hla_4digit"] = tcr["hla_4digit"].fillna("").astype(str)
    tcr["exact_key"] = tcr["peptide"] + "|" + tcr["hla_4digit"]
    tcr["peptide_key"] = tcr["peptide"]
    tcr["near_key"] = tcr["peptide"].map(peptide_cluster).astype(str) + "|" + tcr["hla_supertype"].fillna("").astype(str)

    exact_idx: dict[str, list[int]] = defaultdict(list)
    peptide_idx: dict[str, list[int]] = defaultdict(list)
    near_idx: dict[str, list[int]] = defaultdict(list)
    for idx, row in tcr.iterrows():
        pep = row["peptide"]
        if not pep:
            continue
        exact_idx[row["exact_key"]].append(idx)
        peptide_idx[row["peptide_key"]].append(idx)
        near_idx[row["near_key"]].append(idx)

    rows: list[dict[str, object]] = []
    for _, row in neo.iterrows():
        pep = str(row.get("peptide", "") or "")
        hla = str(row.get("hla_4digit", "") or "")
        supertype = hla_supertype(hla)
        exact_key = f"{pep}|{hla}"
        near_key = f"{peptide_cluster(pep)}|{supertype}"
        source = str(row.get("source_dataset", "") or "")

        def get_hits(index: dict[str, list[int]], key: str) -> pd.DataFrame:
            ids = index.get(key, [])
            if not ids:
                return pd.DataFrame(columns=tcr.columns)
            hits = tcr.loc[ids]
            return hits[~hits["source_dataset"].astype(str).eq(source)].copy()

        exact_hits = get_hits(exact_idx, exact_key)
        peptide_hits = get_hits(peptide_idx, pep)
        near_hits = get_hits(near_idx, near_key)
        if len(exact_hits) and exact_hits.get("paired_tcr_available", False).fillna(False).astype(bool).any():
            category = "exact_tcr_pmhc_match"
            chosen = exact_hits
        elif len(exact_hits):
            category = "peptide_hla_match"
            chosen = exact_hits
        elif len(peptide_hits):
            category = "peptide_only_match"
            chosen = peptide_hits
        elif len(near_hits):
            category = "near_peptide_match"
            chosen = near_hits
        else:
            category = "no_tcr_match"
            chosen = pd.DataFrame(columns=tcr.columns)
        rec = {"row_id": row["row_id"], "tcr_link_category": category}
        rec.update(summarize_tcr_hits(chosen))
        rows.append(rec)
    linked = pd.DataFrame(rows)
    out = tcr_evidence_features(linked, include_public_label_counts=include_public_label_counts)
    rename = {c: f"source_excluded_{c}" for c in out.columns if c != "row_id"}
    return out.rename(columns=rename)


def make_lr() -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    C=0.5,
                    class_weight="balanced",
                    max_iter=3000,
                    solver="liblinear",
                    random_state=SEED,
                ),
            ),
        ]
    )


def score_metrics(y: np.ndarray, score: np.ndarray, prefix: str = "") -> dict[str, float]:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    out: dict[str, float] = {
        f"{prefix}n": float(len(y)),
        f"{prefix}n_pos": float(y.sum()),
        f"{prefix}prevalence": float(y.mean()) if len(y) else np.nan,
    }
    if len(np.unique(y)) == 2:
        out[f"{prefix}auprc"] = float(average_precision_score(y, score))
        out[f"{prefix}auroc"] = float(roc_auc_score(y, score))
    else:
        out[f"{prefix}auprc"] = np.nan
        out[f"{prefix}auroc"] = np.nan
    for k in [5, 10, 20]:
        kk = min(k, len(y))
        if kk:
            # Deterministic tiny jitter prevents arbitrary input-order wins when
            # a control model assigns identical scores to many rows.
            jitter = np.random.default_rng(SEED + k + len(y)).normal(0, 1e-12, len(score))
            idx = np.argsort(-(score + jitter))[:kk]
            precision = float(y[idx].mean())
        else:
            precision = np.nan
        out[f"{prefix}precision_at_{k}"] = precision
        prev = out[f"{prefix}prevalence"]
        out[f"{prefix}enrichment_at_{k}"] = precision / prev if prev and not np.isnan(prev) else np.nan
    return out


def evaluate_main_registry(
    split_names: list[str],
    max_folds_per_split: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    reg = pd.read_csv(OUT / "canonical_registry.tsv", sep="\t", low_memory=False)
    reg["label_binary"] = pd.to_numeric(reg["label_binary"], errors="coerce").fillna(0).astype(int)
    linked = pd.read_parquet(TCR_OUT / "tcr_neo_linked_registry.parquet")

    cf = read_feature(OUT / "features/counterfactual_features.parquet")
    plm = read_feature(OUT / "features/plm_embeddings.parquet")
    struct = read_feature(OUT / "features/structure_features.parquet")
    proc = read_feature(OUT / "features/processing_context_features.parquet")
    tcr = tcr_evidence_features(linked, include_public_label_counts=True)
    tcr_nolabel = tcr_evidence_features(linked, include_public_label_counts=False)
    tcr_exsource = source_excluded_tcr_evidence_features(reg, include_public_label_counts=True)
    tcr_exsource_nolabel = source_excluded_tcr_evidence_features(reg, include_public_label_counts=False)

    feature_sets = {
        "pmhc_counterfactual": numeric_features(cf),
        "pmhc_multimodal": merge_features([cf, plm, struct, proc]),
        "tcr_evidence_only": numeric_features(tcr),
        "tcr_evidence_nolabel_only": numeric_features(tcr_nolabel),
        "tcr_evidence_exsource_only": numeric_features(tcr_exsource),
        "tcr_evidence_exsource_nolabel_only": numeric_features(tcr_exsource_nolabel),
        "pmhc_counterfactual_plus_tcr": merge_features([cf, tcr]),
        "pmhc_counterfactual_plus_tcr_nolabel": merge_features([cf, tcr_nolabel]),
        "pmhc_counterfactual_plus_tcr_exsource": merge_features([cf, tcr_exsource]),
        "pmhc_counterfactual_plus_tcr_exsource_nolabel": merge_features([cf, tcr_exsource_nolabel]),
        "pmhc_multimodal_plus_tcr": merge_features([cf, plm, struct, proc, tcr]),
        "pmhc_multimodal_plus_tcr_nolabel": merge_features([cf, plm, struct, proc, tcr_nolabel]),
        "pmhc_multimodal_plus_tcr_exsource": merge_features([cf, plm, struct, proc, tcr_exsource]),
        "pmhc_multimodal_plus_tcr_exsource_nolabel": merge_features([cf, plm, struct, proc, tcr_exsource_nolabel]),
    }

    meta_cols = ["row_id", "label_binary"]
    meta = reg[meta_cols].merge(
        linked[["row_id", "tcr_link_category", "tcr_evidence_count", "paired_tcr_evidence_count"]],
        on="row_id",
        how="left",
    )
    pred_rows: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []

    for split_name in split_names:
        split_path = OUT / "splits" / f"{split_name}.tsv"
        if not split_path.exists():
            continue
        split = pd.read_csv(split_path, sep="\t", low_memory=False)
        for fold_i, (fold_id, fold) in enumerate(split.groupby("fold_id", sort=True)):
            if fold_i >= max_folds_per_split:
                break
            train_ids = fold.loc[fold["role"].eq("train"), "row_id"].astype(str)
            test_ids = fold.loc[fold["role"].eq("test"), "row_id"].astype(str)
            train_meta = meta[meta["row_id"].astype(str).isin(train_ids)].copy()
            test_meta = meta[meta["row_id"].astype(str).isin(test_ids)].copy()
            if train_meta["label_binary"].nunique() < 2 or test_meta.empty:
                continue
            for model_name, feats in feature_sets.items():
                train = train_meta[["row_id", "label_binary"]].merge(feats, on="row_id", how="left")
                test = test_meta[["row_id", "label_binary"]].merge(feats, on="row_id", how="left")
                x_cols = [c for c in train.columns if c not in {"row_id", "label_binary"}]
                model = make_lr()
                model.fit(train[x_cols], train["label_binary"].values)
                score = model.predict_proba(test[x_cols])[:, 1]
                pred = test_meta.copy()
                pred["split_name"] = split_name
                pred["fold_id"] = fold_id
                pred["model_name"] = model_name
                pred["score"] = score
                pred_rows.append(pred)

                base = {
                    "pilot": "cross_neo_label_tcr_evidence",
                    "split_name": split_name,
                    "fold_id": fold_id,
                    "model_name": model_name,
                    "subset": "all_rows",
                }
                base.update(score_metrics(test["label_binary"].values, score))
                metric_rows.append(base)

                evidence_mask = test_meta["tcr_evidence_count"].fillna(0).astype(float).gt(0).values
                if evidence_mask.sum() >= 20:
                    sub = {
                        "pilot": "cross_neo_label_tcr_evidence",
                        "split_name": split_name,
                        "fold_id": fold_id,
                        "model_name": model_name,
                        "subset": "rows_with_any_tcr_evidence",
                    }
                    sub.update(score_metrics(test["label_binary"].values[evidence_mask], score[evidence_mask]))
                    metric_rows.append(sub)

                exact_mask = test_meta["tcr_link_category"].fillna("").eq("exact_tcr_pmhc_match").values
                if exact_mask.sum() >= 10:
                    sub = {
                        "pilot": "cross_neo_label_tcr_evidence",
                        "split_name": split_name,
                        "fold_id": fold_id,
                        "model_name": model_name,
                        "subset": "exact_tcr_pmhc_match",
                    }
                    sub.update(score_metrics(test["label_binary"].values[exact_mask], score[exact_mask]))
                    metric_rows.append(sub)

    predictions = pd.concat(pred_rows, ignore_index=True) if pred_rows else pd.DataFrame()
    metrics = pd.DataFrame(metric_rows)
    return predictions, metrics


def aa_comp(seq: str, prefix: str) -> dict[str, float]:
    seq = str(seq or "")
    n = max(1, len(seq))
    out = {f"{prefix}_len": float(len(seq))}
    for aa in AA:
        out[f"{prefix}_frac_{aa}"] = seq.count(aa) / n
    return out


def jaccard_kmer(a: str, b: str, k: int) -> float:
    ka = kmer_set(a, k)
    kb = kmer_set(b, k)
    return 0.0 if not ka or not kb else len(ka & kb) / len(ka | kb)


def decoy_feature_dict(row: pd.Series, mode: str) -> dict[str, object]:
    pep = str(row.get("peptide", "") or "")
    ca = str(row.get("cdr3_alpha", "") or "")
    cb = str(row.get("cdr3_beta", "") or "")
    feats: dict[str, object] = {}

    if mode in {"pmhc_only", "pmhc_tcr", "pmhc_tcr_cross"}:
        feats.update(aa_comp(pep, "pep"))
        feats["hla"] = f"hla={row.get('hla_4digit', '')}"
        feats["hla_supertype"] = f"hla_supertype={row.get('hla_supertype', '')}"
        feats["mhc_gene"] = f"mhc_gene={row.get('mhc_gene', '')}"
        feats["mhc_class"] = f"mhc_class={row.get('mhc_class', '')}"

    if mode in {"tcr_only", "pmhc_tcr", "pmhc_tcr_cross"}:
        feats.update(aa_comp(ca, "cdr3a"))
        feats.update(aa_comp(cb, "cdr3b"))
        feats["trav"] = f"trav={row.get('tcr_alpha_v', '')}"
        feats["traj"] = f"traj={row.get('tcr_alpha_j', '')}"
        feats["trbv"] = f"trbv={row.get('tcr_beta_v', '')}"
        feats["trbj"] = f"trbj={row.get('tcr_beta_j', '')}"
        feats["paired"] = float(bool(row.get("paired_tcr_available", True)))

    if mode == "pmhc_tcr_cross":
        feats["pep_cdr3a_len_delta"] = abs(len(pep) - len(ca))
        feats["pep_cdr3b_len_delta"] = abs(len(pep) - len(cb))
        feats["pep_cdr3a_seq_similarity"] = seq_similarity(pep, ca)
        feats["pep_cdr3b_seq_similarity"] = seq_similarity(pep, cb)
        for k in [2, 3]:
            feats[f"pep_cdr3a_k{k}_jaccard"] = jaccard_kmer(pep, ca, k)
            feats[f"pep_cdr3b_k{k}_jaccard"] = jaccard_kmer(pep, cb, k)
        for aa in AA:
            feats[f"pep_x_cdr3a_frac_{aa}"] = feats.get(f"pep_frac_{aa}", 0.0) * feats.get(f"cdr3a_frac_{aa}", 0.0)
            feats[f"pep_x_cdr3b_frac_{aa}"] = feats.get(f"pep_frac_{aa}", 0.0) * feats.get(f"cdr3b_frac_{aa}", 0.0)
    return feats


def make_dict_lr() -> Pipeline:
    return Pipeline(
        [
            ("vec", DictVectorizer(sparse=True)),
            ("scaler", StandardScaler(with_mean=False)),
            (
                "clf",
                LogisticRegression(
                    C=0.5,
                    class_weight="balanced",
                    max_iter=3000,
                    solver="liblinear",
                    random_state=SEED,
                ),
            ),
        ]
    )


def load_decoy_positives(max_positive: int) -> pd.DataFrame:
    cols = [
        "row_id",
        "source_dataset",
        "peptide",
        "hla_4digit",
        "hla_supertype",
        "mhc_class",
        "mhc_gene",
        "tcr_alpha_v",
        "tcr_alpha_j",
        "cdr3_alpha",
        "tcr_beta_v",
        "tcr_beta_j",
        "cdr3_beta",
        "paired_tcr_available",
        "binding_label_binary",
        "exact_peptide_hla_key",
    ]
    reg = pd.read_parquet(TCR_OUT / "tcr_registry.parquet", columns=cols)
    reg = reg[
        reg["paired_tcr_available"].astype(bool)
        & reg["peptide"].astype(str).ne("")
        & reg["hla_4digit"].astype(str).ne("")
        & pd.to_numeric(reg["binding_label_binary"], errors="coerce").eq(1)
        & reg["mhc_class"].astype(str).eq("I")
        & reg["mhc_gene"].astype(str).isin(["A", "B", "C"])
        & reg["peptide"].astype(str).str.len().between(8, 11)
    ].copy()
    if len(reg) > max_positive:
        reg = reg.sample(n=max_positive, random_state=SEED)
    reg = reg.reset_index(drop=True)
    reg["decoy_type"] = "observed_positive"
    reg["label"] = 1
    return reg


def make_shuffled_decoys(pos: pd.DataFrame, seed: int) -> pd.DataFrame:
    pos = pos.reset_index(drop=True)
    donor = pos.sample(frac=1, random_state=seed).reset_index(drop=True)
    for _ in range(10):
        same = donor["exact_peptide_hla_key"].astype(str).values == pos["exact_peptide_hla_key"].astype(str).values
        if not same.any():
            break
        donor.loc[same] = donor.loc[same].sample(frac=1, random_state=seed + 31).values
    if (donor["exact_peptide_hla_key"].astype(str).values == pos["exact_peptide_hla_key"].astype(str).values).any():
        donor = donor.shift(1).fillna(donor.iloc[-1])

    decoy = pos.copy()
    for col in ["tcr_alpha_v", "tcr_alpha_j", "cdr3_alpha", "tcr_beta_v", "tcr_beta_j", "cdr3_beta"]:
        decoy[col] = donor[col].values
    decoy["row_id"] = decoy["row_id"].astype(str) + "_SHUFFLED_TCR_DECOY"
    decoy["decoy_type"] = "same_pmhc_shuffled_tcr"
    decoy["label"] = 0
    return decoy


def evaluate_decoy_pilot(max_positive: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    positives = load_decoy_positives(max_positive=max_positive)
    modes = ["pmhc_only", "tcr_only", "pmhc_tcr", "pmhc_tcr_cross"]
    groups = positives["exact_peptide_hla_key"].astype(str).values
    y_pos = positives["label"].astype(int).values
    n_splits = min(5, max(2, len(np.unique(groups))))
    gkf = GroupKFold(n_splits=n_splits)
    metric_rows: list[dict[str, object]] = []
    pred_rows: list[pd.DataFrame] = []

    for fold_id, (train_idx, test_idx) in enumerate(gkf.split(positives, y_pos, groups=groups)):
        train_pos = positives.iloc[train_idx].copy()
        test_pos = positives.iloc[test_idx].copy()
        train_decoy = make_shuffled_decoys(train_pos, seed=SEED + 1000 + fold_id)
        test_decoy = make_shuffled_decoys(test_pos, seed=SEED + 2000 + fold_id)
        train_data = pd.concat([train_pos, train_decoy], ignore_index=True)
        test_data = pd.concat([test_pos, test_decoy], ignore_index=True)
        y_train = train_data["label"].astype(int).values
        y_test = test_data["label"].astype(int).values
        train_features = {mode: [decoy_feature_dict(row, mode) for _, row in train_data.iterrows()] for mode in modes}
        test_features = {mode: [decoy_feature_dict(row, mode) for _, row in test_data.iterrows()] for mode in modes}
        for mode in modes:
            model = make_dict_lr()
            train_dicts = train_features[mode]
            test_dicts = test_features[mode]
            model.fit(train_dicts, y_train)
            score = model.predict_proba(test_dicts)[:, 1]
            pred = test_data[["row_id", "source_dataset", "peptide", "hla_4digit", "cdr3_alpha", "cdr3_beta", "label", "decoy_type"]].copy()
            pred["fold_id"] = fold_id
            pred["model_name"] = mode
            pred["score"] = score
            pred_rows.append(pred)
            rec = {
                "pilot": "paired_tcr_pmhc_shuffled_decoy",
                "split_name": "grouped_by_peptide_hla",
                "fold_id": fold_id,
                "model_name": mode,
                "subset": "paired_class_i_positive_vs_shuffled_tcr_decoy",
            }
            rec.update(score_metrics(y_test, score))
            rec["n_unique_test_pmhc"] = float(test_data["exact_peptide_hla_key"].nunique())
            metric_rows.append(rec)
    return pd.concat(pred_rows, ignore_index=True), pd.DataFrame(metric_rows)


def aggregate_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return metrics
    metric_cols = [
        c
        for c in metrics.columns
        if c
        not in {
            "pilot",
            "split_name",
            "fold_id",
            "model_name",
            "subset",
        }
        and pd.api.types.is_numeric_dtype(metrics[c])
    ]
    agg = (
        metrics.groupby(["pilot", "split_name", "model_name", "subset"], dropna=False)[metric_cols]
        .agg(["mean", "std"])
        .reset_index()
    )
    agg.columns = ["_".join([str(x) for x in col if str(x)]) for col in agg.columns]
    return agg


def key_deltas(aggregate: pd.DataFrame) -> pd.DataFrame:
    main = aggregate[
        aggregate["pilot"].eq("cross_neo_label_tcr_evidence")
        & aggregate["subset"].eq("all_rows")
    ].copy()
    pairs = [
        ("pmhc_counterfactual", "pmhc_counterfactual_plus_tcr"),
        ("pmhc_counterfactual", "pmhc_counterfactual_plus_tcr_nolabel"),
        ("pmhc_counterfactual", "pmhc_counterfactual_plus_tcr_exsource"),
        ("pmhc_counterfactual", "pmhc_counterfactual_plus_tcr_exsource_nolabel"),
        ("pmhc_multimodal", "pmhc_multimodal_plus_tcr"),
        ("pmhc_multimodal", "pmhc_multimodal_plus_tcr_nolabel"),
        ("pmhc_multimodal", "pmhc_multimodal_plus_tcr_exsource"),
        ("pmhc_multimodal", "pmhc_multimodal_plus_tcr_exsource_nolabel"),
        ("pmhc_counterfactual", "tcr_evidence_only"),
        ("pmhc_counterfactual", "tcr_evidence_nolabel_only"),
        ("pmhc_counterfactual", "tcr_evidence_exsource_only"),
        ("pmhc_counterfactual", "tcr_evidence_exsource_nolabel_only"),
    ]
    rows: list[dict[str, object]] = []
    for split_name in sorted(main["split_name"].dropna().unique()):
        for baseline, comparison in pairs:
            a = main[main["split_name"].eq(split_name) & main["model_name"].eq(baseline)]
            b = main[main["split_name"].eq(split_name) & main["model_name"].eq(comparison)]
            if a.empty or b.empty:
                continue
            aa = a.iloc[0]
            bb = b.iloc[0]
            rows.append({
                "split_name": split_name,
                "baseline_model": baseline,
                "comparison_model": comparison,
                "baseline_auprc": aa.get("auprc_mean", np.nan),
                "comparison_auprc": bb.get("auprc_mean", np.nan),
                "delta_auprc": bb.get("auprc_mean", np.nan) - aa.get("auprc_mean", np.nan),
                "baseline_auroc": aa.get("auroc_mean", np.nan),
                "comparison_auroc": bb.get("auroc_mean", np.nan),
                "delta_auroc": bb.get("auroc_mean", np.nan) - aa.get("auroc_mean", np.nan),
                "baseline_p20": aa.get("precision_at_20_mean", np.nan),
                "comparison_p20": bb.get("precision_at_20_mean", np.nan),
                "delta_p20": bb.get("precision_at_20_mean", np.nan) - aa.get("precision_at_20_mean", np.nan),
            })
    return pd.DataFrame(rows)


def write_report(main_metrics: pd.DataFrame, decoy_metrics: pd.DataFrame, aggregate: pd.DataFrame) -> None:
    PILOT_OUT.mkdir(parents=True, exist_ok=True)

    def best_table(pilot: str, subset: str, split: str | None = None) -> str:
        df = aggregate[(aggregate["pilot"].eq(pilot)) & (aggregate["subset"].eq(subset))].copy()
        if split is not None:
            df = df[df["split_name"].eq(split)]
        keep = [
            "split_name",
            "model_name",
            "n_mean",
            "n_pos_mean",
            "prevalence_mean",
            "auprc_mean",
            "auroc_mean",
            "precision_at_20_mean",
            "enrichment_at_20_mean",
        ]
        keep = [c for c in keep if c in df.columns]
        if df.empty:
            return "_No rows._"
        return df[keep].sort_values([c for c in ["split_name", "auprc_mean"] if c in keep], ascending=[True, False] if "split_name" in keep else [False]).to_markdown(index=False)

    lines = [
        "# CROSS-Neo-TCR Quick Subset Comparison",
        "",
        "## Claim Boundary",
        "",
        "This is a fast diagnostic pilot. It is not a final TCR-aware SOTA benchmark and does not transfer public TCR labels into the main CROSS-Neo task.",
        "",
        "## Pilot A: Main CROSS-Neo Labels With TCR Evidence Features",
        "",
        "Question: on existing CROSS-Neo rows, does adding TCR registry evidence improve ranking compared with pMHC-only features?",
        "",
        best_table("cross_neo_label_tcr_evidence", "all_rows"),
        "",
        "### Rows With Any TCR Evidence",
        "",
        best_table("cross_neo_label_tcr_evidence", "rows_with_any_tcr_evidence"),
        "",
        "## Pilot B: Paired TCR-pMHC Positives Versus Shuffled TCR Decoys",
        "",
        "Question: on paired class-I public positives, can simple sequence and peptide-TCR cross features distinguish observed TCR-pMHC pairs from same-pMHC shuffled-TCR decoys?",
        "",
        best_table("paired_tcr_pmhc_shuffled_decoy", "paired_class_i_positive_vs_shuffled_tcr_decoy"),
        "",
        "## Interpretation",
        "",
        "- A TCR gain in Pilot A is diagnostic evidence that local TCR evidence overlaps with useful ranking signal, not proof of universal TCR recognition prediction.",
        "- Pilot B is a decoy stress test. It checks whether the current simple TCR sequence features are non-random before running expensive structure jobs.",
        "- Strong claims still require paired TCR alpha/beta data, strict TCR/peptide-HLA/source-heldout splits, and external validation.",
        "",
        "## Outputs",
        "",
        "- `cross_neo_tcr_evidence_pilot_metrics.tsv`",
        "- `cross_neo_tcr_evidence_pilot_predictions.parquet`",
        "- `paired_tcr_decoy_pilot_metrics.tsv`",
        "- `paired_tcr_decoy_pilot_predictions.parquet`",
        "- `quick_tcr_subset_compare_aggregate.tsv`",
        "- `quick_tcr_subset_compare_key_deltas.tsv`",
    ]
    (PILOT_OUT / "quick_tcr_subset_compare_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--splits", nargs="*", default=DEFAULT_SPLITS)
    parser.add_argument("--max-folds-per-split", type=int, default=5)
    parser.add_argument("--max-decoy-positives", type=int, default=5000)
    args = parser.parse_args()

    PILOT_OUT.mkdir(parents=True, exist_ok=True)
    main_pred, main_metrics = evaluate_main_registry(args.splits, args.max_folds_per_split)
    decoy_pred, decoy_metrics = evaluate_decoy_pilot(args.max_decoy_positives)
    aggregate = aggregate_metrics(pd.concat([main_metrics, decoy_metrics], ignore_index=True))
    deltas = key_deltas(aggregate)

    write_tsv(main_metrics, PILOT_OUT / "cross_neo_tcr_evidence_pilot_metrics.tsv")
    write_tsv(decoy_metrics, PILOT_OUT / "paired_tcr_decoy_pilot_metrics.tsv")
    write_tsv(aggregate, PILOT_OUT / "quick_tcr_subset_compare_aggregate.tsv")
    write_tsv(deltas, PILOT_OUT / "quick_tcr_subset_compare_key_deltas.tsv")
    safe_parquet(main_pred, PILOT_OUT / "cross_neo_tcr_evidence_pilot_predictions.parquet")
    safe_parquet(decoy_pred, PILOT_OUT / "paired_tcr_decoy_pilot_predictions.parquet")
    write_report(main_metrics, decoy_metrics, aggregate)

    main_best = aggregate[
        aggregate["pilot"].eq("cross_neo_label_tcr_evidence")
        & aggregate["subset"].eq("all_rows")
        & aggregate["split_name"].eq("exact_peptide_hla_holdout")
    ].sort_values("auprc_mean", ascending=False).head(1)
    decoy_best = aggregate[
        aggregate["pilot"].eq("paired_tcr_pmhc_shuffled_decoy")
    ].sort_values("auprc_mean", ascending=False).head(1)
    print(f"[tcr-pilot] out={PILOT_OUT}")
    if not main_best.empty:
        r = main_best.iloc[0]
        print(f"[tcr-pilot] main exact-holdout best={r['model_name']} auprc={r['auprc_mean']:.4f} auroc={r['auroc_mean']:.4f}")
    if not decoy_best.empty:
        r = decoy_best.iloc[0]
        print(f"[tcr-pilot] decoy best={r['model_name']} auprc={r['auprc_mean']:.4f} auroc={r['auroc_mean']:.4f}")


if __name__ == "__main__":
    main()
