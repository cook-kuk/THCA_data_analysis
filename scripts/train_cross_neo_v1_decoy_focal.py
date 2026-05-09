#!/usr/bin/env python3
"""Decoy-negative + focal-style positive-pattern ranker for CROSS-Neo v1.

This branch is intentionally small and split-safe:
- positive pattern statistics are computed from the train fold only;
- decoys are generated from train-fold positives only;
- labeled negatives are downweighted as PU/ambiguous;
- focal-style weights are fit only on train + generated decoys.
"""

from __future__ import annotations

import json
import warnings
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler

from cross_neo_v1_common import V1, ensure_v1_dirs, get_fold_ids, load_folds, load_master, metrics, seq_similarity


AA = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_I = {a: i for i, a in enumerate(AA)}
ANCHOR_POS = {1, -1}  # 0-based p2 and C-terminal anchor proxy.
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning)


def clean_pep(x: object) -> str:
    return "".join(a for a in str(x or "").upper() if a in AA)


def make_decoy(pep: str, bg: list[str], rng: np.random.Generator) -> str:
    if not pep:
        return pep
    chars = list(pep)
    mutable = [i for i in range(len(chars)) if i not in {1, len(chars) - 1}]
    if not mutable:
        mutable = list(range(len(chars)))
    n_mut = max(1, min(3, int(round(len(chars) * 0.25))))
    for i in rng.choice(mutable, size=min(n_mut, len(mutable)), replace=False):
        old = chars[i]
        choices = [a for a in bg if a != old] or [a for a in AA if a != old]
        chars[i] = str(rng.choice(choices))
    decoy = "".join(chars)
    if decoy == pep and len(chars):
        i = int(rng.choice(mutable))
        choices = [a for a in AA if a != chars[i]]
        chars[i] = str(rng.choice(choices))
        decoy = "".join(chars)
    return decoy


class PositivePattern:
    def __init__(self, train: pd.DataFrame, smooth: float = 0.5):
        self.smooth = smooth
        self.pos = train[train["label"].astype(int) == 1].copy()
        self.bg = train.copy()
        self.global_pwm = self._pwm(self.pos)
        self.bg_pwm = self._pwm(self.bg)
        self.hla_pwm = {}
        for st, sub in self.pos.groupby("hla_supertype"):
            if len(sub) >= 3:
                self.hla_pwm[str(st)] = self._pwm(sub)
        self.pos_peps = [clean_pep(p) for p in self.pos["peptide_mut"]]
        self.pos_by_len = defaultdict(list)
        for p in self.pos_peps:
            self.pos_by_len[len(p)].append(p)
        self.bg_aa = []
        for p in self.bg["peptide_mut"]:
            self.bg_aa.extend(list(clean_pep(p)))
        if not self.bg_aa:
            self.bg_aa = list(AA)

    def _pwm(self, df: pd.DataFrame) -> dict[int, np.ndarray]:
        counts: dict[int, np.ndarray] = {}
        for pep in df["peptide_mut"]:
            p = clean_pep(pep)
            for i, aa in enumerate(p):
                if aa not in AA_TO_I:
                    continue
                counts.setdefault(i, np.zeros(len(AA), dtype=float))[AA_TO_I[aa]] += 1
        for i in range(15):
            arr = counts.setdefault(i, np.zeros(len(AA), dtype=float))
            arr += self.smooth
            counts[i] = arr / arr.sum()
        return counts

    def logodds(self, pep: str, hla_supertype: str = "") -> tuple[float, float]:
        p = clean_pep(pep)
        g = 0.0
        h = 0.0
        hpwm = self.hla_pwm.get(str(hla_supertype), self.global_pwm)
        for i, aa in enumerate(p[:15]):
            if aa not in AA_TO_I:
                continue
            idx = AA_TO_I[aa]
            bg = max(self.bg_pwm.get(i, np.ones(len(AA)) / len(AA))[idx], 1e-6)
            gp = max(self.global_pwm.get(i, np.ones(len(AA)) / len(AA))[idx], 1e-6)
            hp = max(hpwm.get(i, np.ones(len(AA)) / len(AA))[idx], 1e-6)
            g += float(np.log(gp / bg))
            h += float(np.log(hp / bg))
        return g, h

    def nearest_pos_similarity(self, pep: str) -> float:
        p = clean_pep(pep)
        refs = self.pos_by_len.get(len(p), self.pos_peps)
        return max((seq_similarity(p, q) for q in refs), default=0.0)

    def row_features(self, df: pd.DataFrame) -> pd.DataFrame:
        rows = []
        for _, r in df.iterrows():
            pep = clean_pep(r.get("peptide_mut", ""))
            st = str(r.get("hla_supertype", ""))
            glo, hlo = self.logodds(pep, st)
            aa_counts = Counter(pep)
            vals = {
                "sample_id": r.get("sample_id", ""),
                "label": int(r.get("label", 0)),
                "peptide_len": len(pep),
                "pattern_global_logodds": glo,
                "pattern_hla_logodds": hlo,
                "nearest_train_positive_similarity": self.nearest_pos_similarity(pep),
                "anchor_aromatic": float(sum(pep[i] in "FWY" for i in [1, len(pep) - 1] if 0 <= i < len(pep))),
                "anchor_hydrophobic": float(sum(pep[i] in "AILMFWVY" for i in [1, len(pep) - 1] if 0 <= i < len(pep))),
                "nonanchor_aromatic_frac": float(sum(a in "FWY" for i, a in enumerate(pep) if i not in {1, len(pep) - 1}) / max(1, len(pep) - 2)),
                "nonanchor_charged_frac": float(sum(a in "DEKRH" for i, a in enumerate(pep) if i not in {1, len(pep) - 1}) / max(1, len(pep) - 2)),
                "hla_seen_in_positive_train": int(str(r.get("hla", "")) in set(self.pos["hla"].astype(str))),
                "hla_supertype_seen_in_positive_train": int(st in set(self.pos["hla_supertype"].astype(str))),
            }
            for aa in AA:
                vals[f"aa_frac_{aa}"] = aa_counts.get(aa, 0) / max(1, len(pep))
            for pos in range(9):
                aa = pep[pos] if pos < len(pep) else "X"
                for sym in AA:
                    vals[f"pos{pos+1:02d}_{sym}"] = int(aa == sym)
            rows.append(vals)
        return pd.DataFrame(rows)

    def positive_pattern_table(self) -> pd.DataFrame:
        rows = []
        for i in range(15):
            gp = self.global_pwm.get(i, np.ones(len(AA)) / len(AA))
            bg = self.bg_pwm.get(i, np.ones(len(AA)) / len(AA))
            enrich = np.log(np.clip(gp, 1e-6, 1) / np.clip(bg, 1e-6, 1))
            top = np.argsort(-enrich)[:5]
            rows.append(
                {
                    "position_1based": i + 1,
                    "top_positive_enriched_residues": ",".join(AA[j] for j in top),
                    "top_logodds": ",".join(f"{enrich[j]:.3f}" for j in top),
                    "positive_entropy": float(-(gp * np.log2(gp)).sum()),
                    "background_entropy": float(-(bg * np.log2(bg)).sum()),
                }
            )
        return pd.DataFrame(rows)


def build_train_with_decoys(train: pd.DataFrame, pattern: PositivePattern, n_decoys_per_pos: int = 4) -> pd.DataFrame:
    rng = np.random.default_rng(20260509 + len(train))
    pieces = [train.copy()]
    decoys = []
    pos = train[train["label"].astype(int) == 1]
    existing = set(train["peptide_mut"].astype(str))
    for _, r in pos.iterrows():
        pep = clean_pep(r["peptide_mut"])
        made = 0
        tries = 0
        while made < n_decoys_per_pos and tries < n_decoys_per_pos * 8:
            tries += 1
            d = make_decoy(pep, pattern.bg_aa, rng)
            if not d or d in existing:
                continue
            rr = r.copy()
            rr["sample_id"] = f"{r['sample_id']}__decoy{made:02d}"
            rr["peptide_mut"] = d
            rr["label"] = 0
            rr["is_decoy"] = 1
            decoys.append(rr)
            made += 1
    if decoys:
        pieces.append(pd.DataFrame(decoys))
    out = pd.concat(pieces, ignore_index=True)
    out["is_decoy"] = out.get("is_decoy", 0)
    out["is_decoy"] = out["is_decoy"].fillna(0).astype(int)
    return out


def focal_fit_predict(train_feat: pd.DataFrame, test_feat: pd.DataFrame, cols: list[str], gamma: float = 2.0, alpha: float = 0.75) -> np.ndarray:
    y = train_feat["label"].to_numpy(int)
    if len(np.unique(y)) < 2:
        return np.full(len(test_feat), y.mean() if len(y) else 0.5)
    xtr = train_feat[cols].fillna(0).to_numpy(float)
    xte = test_feat[cols].fillna(0).to_numpy(float)
    scaler = StandardScaler()
    xtr = scaler.fit_transform(xtr)
    xte = scaler.transform(xte)
    base_w = np.ones(len(y), dtype=float)
    base_w[(y == 0) & (train_feat.get("is_decoy", pd.Series(0, index=train_feat.index)).to_numpy(int) == 0)] = 0.25
    base_w[(y == 0) & (train_feat.get("is_decoy", pd.Series(0, index=train_feat.index)).to_numpy(int) == 1)] = 1.0
    base_w[y == 1] = 1.5
    sample_w = base_w.copy()
    clf = LogisticRegression(max_iter=400, C=0.35, solver="liblinear")
    for _ in range(2):
        clf.fit(xtr, y, sample_weight=sample_w)
        p = clf.predict_proba(xtr)[:, 1]
        pt = np.where(y == 1, p, 1 - p)
        a = np.where(y == 1, alpha, 1 - alpha)
        sample_w = base_w * a * np.power(1 - np.clip(pt, 1e-4, 1 - 1e-4), gamma)
        sample_w = sample_w * len(sample_w) / max(sample_w.sum(), 1e-9)
    return clf.predict_proba(xte)[:, 1]


def evaluate_locked(master: pd.DataFrame, folds: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pred_rows = []
    pattern_rows = []
    decoy_rows = []
    for (split, fold), _ in folds.groupby(["split_name", "fold_id"], sort=False):
        train_ids, test_ids = get_fold_ids(master, folds, split, fold)
        train = master[master["sample_id"].isin(train_ids)].copy()
        test = master[master["sample_id"].isin(test_ids)].copy()
        if len(train) < 12 or len(test) < 2 or train["label"].nunique() < 2 or test["label"].nunique() < 2:
            continue
        pattern = PositivePattern(train)
        aug = build_train_with_decoys(train, pattern)
        trf = pattern.row_features(aug)
        tef = pattern.row_features(test)
        trf["is_decoy"] = aug["is_decoy"].to_numpy(int)
        cols = [c for c in trf.columns if c not in {"sample_id", "label", "is_decoy"}]
        score = focal_fit_predict(trf, tef, cols)
        for sid, y, p in zip(test["sample_id"], test["label"], score):
            pred_rows.append({"split_name": split, "fold_id": fold, "model": "decoy_focal_positive_pattern", "sample_id": sid, "label": int(y), "score": float(p)})
        pt = pattern.positive_pattern_table()
        pt["split_name"] = split
        pt["fold_id"] = fold
        pattern_rows.append(pt)
        decoy_rows.append({"split_name": split, "fold_id": fold, "train_n": len(train), "train_pos": int(train["label"].sum()), "decoy_n": int(aug["is_decoy"].sum()), "features_n": len(cols)})
        print(f"[v1-decoy-focal] split={split} fold={fold} train={len(train)} decoys={int(aug['is_decoy'].sum())}", flush=True)
    return pd.DataFrame(pred_rows), pd.concat(pattern_rows, ignore_index=True), pd.DataFrame(decoy_rows)


def evaluate_source(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    pool = master[master["study"].isin(["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"])].copy()
    for heldout in ["NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation", "CEDAR"]:
        train = pool[pool["study"] != heldout].copy()
        test = pool[pool["study"] == heldout].copy()
        if len(test) < 10 or train["label"].nunique() < 2 or test["label"].nunique() < 2:
            continue
        pattern = PositivePattern(train)
        aug = build_train_with_decoys(train, pattern)
        trf = pattern.row_features(aug)
        tef = pattern.row_features(test)
        trf["is_decoy"] = aug["is_decoy"].to_numpy(int)
        cols = [c for c in trf.columns if c not in {"sample_id", "label", "is_decoy"}]
        score = focal_fit_predict(trf, tef, cols)
        for sid, y, p in zip(test["sample_id"], test["label"], score):
            rows.append({"split_name": f"source_heldout_{heldout}", "fold_id": heldout, "model": "decoy_focal_positive_pattern", "sample_id": sid, "label": int(y), "score": float(p)})
        print(f"[v1-decoy-focal] source_heldout={heldout} train={len(train)} decoys={int(aug['is_decoy'].sum())}", flush=True)
    return pd.DataFrame(rows)


def main() -> None:
    ensure_v1_dirs()
    master = load_master().copy()
    folds = load_folds()
    locked_pred, pattern_table, decoy_manifest = evaluate_locked(master, folds)
    source_pred = evaluate_source(master)
    pred = pd.concat([locked_pred, source_pred], ignore_index=True)
    pred.to_csv(V1 / "decoy_focal_predictions.tsv", sep="\t", index=False)
    pattern_table.to_csv(V1 / "positive_pattern_pwm.tsv", sep="\t", index=False)
    decoy_manifest.to_csv(V1 / "decoy_generation_manifest.tsv", sep="\t", index=False)
    rows = []
    for (split, model), sub in pred.groupby(["split_name", "model"]):
        rows.append({"split_name": split, "model": model, **metrics(sub["label"].to_numpy(), sub["score"].to_numpy())})
    met = pd.DataFrame(rows).sort_values(["split_name", "AUPRC"], ascending=[True, False])
    met.to_csv(V1 / "decoy_focal_metrics.tsv", sep="\t", index=False)
    # Descriptive all-strict positive pattern, not used for fold training.
    strict = master[master["strict_set_flag"].astype(bool)].copy()
    desc = PositivePattern(strict).positive_pattern_table()
    lines = [
        "# CROSS-Neo v1 Decoy/Focal Positive-Pattern Report",
        "",
        "This branch implements the earlier decoy-negative idea in a fold-safe way. Decoys and positive PWM/log-odds features are generated from each train fold only.",
        "",
        "## Metrics",
        met.sort_values("AUPRC", ascending=False).head(20).to_markdown(index=False),
        "",
        "## Decoy Manifest Summary",
        decoy_manifest.describe(include="all").to_markdown(),
        "",
        "## Descriptive Strict Positive Pattern",
        desc.to_markdown(index=False),
        "",
        "Interpretation: this is a focal-style reweighted logistic ranker over positive-pattern and decoy features; public predictor scores are not used.",
    ]
    (V1 / "decoy_focal_positive_pattern_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-decoy-focal] predictions={len(pred)} metrics={len(met)}")


if __name__ == "__main__":
    main()
