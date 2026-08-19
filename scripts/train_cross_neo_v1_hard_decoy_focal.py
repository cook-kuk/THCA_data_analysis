#!/usr/bin/env python3
"""Hard-decoy focal repair for CROSS-Neo v1.

This is intentionally conservative. Hard decoys are generated inside each
train fold and train only an auxiliary sequence-pattern model. The locked
prediction stays anchored to existing fold-safe C_counterfactual + QK scores.
"""

from __future__ import annotations

import math
import warnings
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from cross_neo_v1_common import (
    SEED,
    V1,
    ensure_v1_dirs,
    expert_wide_for_fold,
    get_fold_ids,
    load_folds,
    load_master,
    metrics,
    read_v0_expert_predictions,
    seq_similarity,
)

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_I = {a: i for i, a in enumerate(AA)}
HYDRO = set("AILMFWVY")
AROM = set("FWY")
CHARGE = set("DEKRH")

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


def clean_pep(x: object) -> str:
    return "".join(a for a in str(x or "").upper() if a in AA)


def anchor_pos(n: int) -> set[int]:
    out = {n - 1} if n > 0 else set()
    if n > 2:
        out.add(1)
    return out


def ent(p: np.ndarray) -> float:
    q = np.clip(np.asarray(p, dtype=float), 1e-12, 1.0)
    return float(-(q * np.log2(q)).sum())


class PatternProfile:
    def __init__(self, train: pd.DataFrame, smooth: float = 0.75):
        self.train = train.copy()
        self.pos = train[train["label"].astype(int) == 1].copy()
        self.smooth = smooth
        self.pos_pwm = self._pwm(self.pos)
        self.bg_pwm = self._pwm(self.train)
        self.st_pwm = {str(k): self._pwm(v) for k, v in self.pos.groupby("hla_supertype") if len(v) >= 3}
        self.pos_peps = [clean_pep(p) for p in self.pos["peptide_mut"]]
        self.pos_by_len = defaultdict(list)
        for p in self.pos_peps:
            self.pos_by_len[len(p)].append(p)
        self.bg_residues = self._bg_residues()
        self.reliability = self._reliability()

    def _pwm(self, rows: pd.DataFrame) -> dict[int, np.ndarray]:
        counts: dict[int, np.ndarray] = {}
        for pep in rows["peptide_mut"].astype(str):
            for i, aa in enumerate(clean_pep(pep)[:15]):
                counts.setdefault(i, np.zeros(len(AA), dtype=float))[AA_TO_I[aa]] += 1
        for i in range(15):
            arr = counts.setdefault(i, np.zeros(len(AA), dtype=float))
            arr += self.smooth
            counts[i] = arr / arr.sum()
        return counts

    def _bg_residues(self) -> dict[tuple[str, int, int], list[str]]:
        out: dict[tuple[str, int, int], list[str]] = defaultdict(list)
        for _, r in self.train.iterrows():
            pep = clean_pep(r.get("peptide_mut", ""))
            st = str(r.get("hla_supertype", ""))
            for i, aa in enumerate(pep):
                out[(st, len(pep), i)].append(aa)
                out[("*", len(pep), i)].append(aa)
                out[("*", -1, i)].append(aa)
        return out

    def _reliability(self) -> float:
        abs_lods = []
        entropy_delta = []
        for i in range(15):
            lod = np.log(np.clip(self.pos_pwm[i], 1e-8, 1) / np.clip(self.bg_pwm[i], 1e-8, 1))
            abs_lods.append(float(np.mean(np.abs(lod))))
            entropy_delta.append(ent(self.bg_pwm[i]) - ent(self.pos_pwm[i]))
        x = (len(self.pos) - 10) / 7 + (np.mean(abs_lods) - 0.08) * 4 + np.mean(entropy_delta) * 0.4
        return float(np.clip(1 / (1 + math.exp(-max(-40, min(40, x)))), 0, 0.75))

    def lod(self, pep: str, st: str = "") -> tuple[float, float]:
        p = clean_pep(pep)
        st_pwm = self.st_pwm.get(str(st), self.pos_pwm)
        g = s = 0.0
        for i, aa in enumerate(p[:15]):
            j = AA_TO_I[aa]
            bg = max(self.bg_pwm[i][j], 1e-8)
            g += math.log(max(self.pos_pwm[i][j], 1e-8) / bg)
            s += math.log(max(st_pwm[i][j], 1e-8) / bg)
        return float(g), float(s)

    def nearest_pos(self, pep: str) -> float:
        p = clean_pep(pep)
        refs = self.pos_by_len.get(len(p), self.pos_peps)
        return max((seq_similarity(p, q) for q in refs), default=0.0)

    def sample_aa(self, st: str, n: int, i: int, old: str, rng: np.random.Generator) -> str:
        pool = self.bg_residues.get((st, n, i)) or self.bg_residues.get(("*", n, i)) or self.bg_residues.get(("*", -1, i)) or list(AA)
        choices = [a for a in pool if a != old] or [a for a in AA if a != old]
        return str(rng.choice(choices))

    def features(self, rows: pd.DataFrame) -> pd.DataFrame:
        pos_hlas = set(self.pos["hla"].astype(str))
        pos_sts = set(self.pos["hla_supertype"].astype(str))
        out = []
        for _, r in rows.iterrows():
            pep = clean_pep(r.get("peptide_mut", ""))
            st = str(r.get("hla_supertype", ""))
            glo, slo = self.lod(pep, st)
            anc = anchor_pos(len(pep))
            non = [a for i, a in enumerate(pep) if i not in anc]
            cnt = Counter(pep)
            vals = {
                "sample_id": r.get("sample_id", ""),
                "label": int(r.get("label", 0)),
                "peptide_len": len(pep),
                "pattern_global_lod": glo,
                "pattern_supertype_lod": slo,
                "nearest_train_positive_similarity": self.nearest_pos(pep),
                "anchor_hydrophobic": float(sum(pep[i] in HYDRO for i in anc if i < len(pep))),
                "anchor_aromatic": float(sum(pep[i] in AROM for i in anc if i < len(pep))),
                "nonanchor_hydrophobic_frac": float(sum(a in HYDRO for a in non) / max(1, len(non))),
                "nonanchor_aromatic_frac": float(sum(a in AROM for a in non) / max(1, len(non))),
                "nonanchor_charged_frac": float(sum(a in CHARGE for a in non) / max(1, len(non))),
                "hla_seen_pos": int(str(r.get("hla", "")) in pos_hlas),
                "supertype_seen_pos": int(st in pos_sts),
            }
            for aa in AA:
                vals[f"aa_frac_{aa}"] = cnt.get(aa, 0) / max(1, len(pep))
            for i in range(9):
                here = pep[i] if i < len(pep) else "X"
                for aa in AA:
                    vals[f"pos{i+1:02d}_{aa}"] = int(here == aa)
            out.append(vals)
        return pd.DataFrame(out)

    def pwm_report(self) -> pd.DataFrame:
        rows = []
        for i in range(15):
            lod = np.log(np.clip(self.pos_pwm[i], 1e-8, 1) / np.clip(self.bg_pwm[i], 1e-8, 1))
            top = np.argsort(-lod)[:5]
            rows.append(
                {
                    "position_1based": i + 1,
                    "top_positive_enriched_residues": ",".join(AA[j] for j in top),
                    "top_logodds": ",".join(f"{lod[j]:.3f}" for j in top),
                    "positive_entropy": ent(self.pos_pwm[i]),
                    "background_entropy": ent(self.bg_pwm[i]),
                }
            )
        return pd.DataFrame(rows)


def make_hard_decoys(train: pd.DataFrame, prof: PatternProfile, n_keep: int = 8) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + len(train) + int(train["label"].sum()) * 17)
    existing = set(train["peptide_mut"].astype(str))
    decoys = []
    for _, r in train[train["label"].astype(int) == 1].iterrows():
        pep = clean_pep(r.get("peptide_mut", ""))
        mutable = [i for i in range(len(pep)) if i not in anchor_pos(len(pep))]
        if len(mutable) < 2:
            continue
        cand = []
        for _ in range(72):
            chars = list(pep)
            n_mut = min(len(mutable), int(rng.integers(2, min(4, len(mutable)) + 1)))
            for i in rng.choice(mutable, size=n_mut, replace=False):
                chars[int(i)] = prof.sample_aa(str(r.get("hla_supertype", "")), len(pep), int(i), chars[int(i)], rng)
            d = "".join(chars)
            if d == pep or d in existing:
                continue
            sim = seq_similarity(pep, d)
            # seq_similarity averages identity with 3/4/5-mer Jaccard scores;
            # two substitutions in a 9-mer are biologically near but often fall
            # below 0.55. Keep a wider hard-decoy band and rank by similarity.
            if sim < 0.30:
                continue
            g, s = prof.lod(d, str(r.get("hla_supertype", "")))
            cand.append((sim + 0.08 * g + 0.08 * s, d))
        if not cand:
            # Last-resort hard decoys: one non-anchor edit. This prevents
            # degenerate zero-decoy folds while still avoiding anchor-only
            # artifacts.
            for i in mutable[: min(4, len(mutable))]:
                chars = list(pep)
                chars[i] = prof.sample_aa(str(r.get("hla_supertype", "")), len(pep), int(i), chars[i], rng)
                d = "".join(chars)
                if d != pep and d not in existing:
                    cand.append((seq_similarity(pep, d), d))
        for j, (_, d) in enumerate(sorted(cand, reverse=True)[:n_keep]):
            rr = r.copy()
            rr["sample_id"] = f"{r['sample_id']}__harddecoy{j:02d}"
            rr["peptide_mut"] = d
            rr["label"] = 0
            rr["is_hard_decoy"] = 1
            decoys.append(rr)
            existing.add(d)
    return pd.DataFrame(decoys)


def focal_predict(xtr: pd.DataFrame, y: np.ndarray, xte: pd.DataFrame, sample_weight: np.ndarray | None = None) -> np.ndarray:
    y = np.asarray(y, dtype=int)
    if len(np.unique(y)) < 2:
        return np.full(len(xte), float(y.mean()) if len(y) else 0.5)
    scaler = StandardScaler()
    xt = scaler.fit_transform(xtr.fillna(0).to_numpy(float))
    xv = scaler.transform(xte.fillna(0).to_numpy(float))
    base = np.ones(len(y)) if sample_weight is None else np.asarray(sample_weight, dtype=float)
    sw = base.copy()
    clf = LogisticRegression(max_iter=500, C=0.35, solver="liblinear")
    for _ in range(3):
        clf.fit(xt, y, sample_weight=sw)
        p = clf.predict_proba(xt)[:, 1]
        pt = np.where(y == 1, p, 1 - p)
        alpha = np.where(y == 1, 0.70, 0.30)
        sw = base * alpha * np.power(1 - np.clip(pt, 1e-4, 1 - 1e-4), 2.0)
        sw *= len(sw) / max(sw.sum(), 1e-9)
    return clf.predict_proba(xv)[:, 1]


def aux_score(train: pd.DataFrame, test: pd.DataFrame) -> tuple[np.ndarray, dict[str, float], pd.DataFrame]:
    prof = PatternProfile(train)
    decoys = make_hard_decoys(train, prof)
    pos = train[train["label"].astype(int) == 1].copy()
    pos["is_hard_decoy"] = 0
    aug = pd.concat([pos, decoys], ignore_index=True, sort=False) if len(decoys) else pos
    feat_tr = prof.features(aug)
    feat_te = prof.features(test)
    cols = [c for c in feat_tr.columns if c not in {"sample_id", "label"}]
    sw = np.where(feat_tr["label"].to_numpy(int) == 1, 1.35, 1.0)
    pred = focal_predict(feat_tr[cols], feat_tr["label"].to_numpy(int), feat_te[cols], sw)
    info = {
        "train_n": len(train),
        "train_pos": int(train["label"].sum()),
        "hard_decoy_n": int(len(decoys)),
        "pattern_reliability": float(prof.reliability),
    }
    return pred, info, prof.pwm_report()


def inner_oof_aux(train: pd.DataFrame) -> pd.DataFrame:
    y = train["label"].to_numpy(int)
    if len(train) < 24 or len(np.unique(y)) < 2 or min(np.bincount(y)) < 3:
        p, info, _ = aux_score(train, train)
        return pd.DataFrame({"sample_id": train["sample_id"].values, "decoy_aux": p, "decoy_reliability": info["pattern_reliability"]})
    n_splits = min(3, int(min(np.bincount(y))))
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    parts = []
    train = train.reset_index(drop=True)
    for tr_idx, va_idx in skf.split(train, y):
        tr = train.iloc[tr_idx].copy()
        va = train.iloc[va_idx].copy()
        p, info, _ = aux_score(tr, va)
        parts.append(pd.DataFrame({"sample_id": va["sample_id"].values, "decoy_aux": p, "decoy_reliability": info["pattern_reliability"]}))
    return pd.concat(parts, ignore_index=True)


def build_meta(rows: pd.DataFrame, wide: pd.DataFrame, aux: pd.DataFrame) -> pd.DataFrame:
    out = rows[["sample_id", "label", "peptide_mut", "hla", "hla_supertype", "study"]].copy()
    out["peptide_len"] = out["peptide_mut"].map(lambda x: len(clean_pep(x)))
    out = out.merge(wide, on=["sample_id", "label"], how="left")
    out = out.merge(aux, on="sample_id", how="left")
    for c in ["C_counterfactual_rf", "qk_no_anchor_gamma1", "qk_quantum_only_gamma1", "D_structure_geometry_rf"]:
        if c not in out.columns:
            out[c] = 0.5
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.5)
    out["decoy_aux"] = pd.to_numeric(out["decoy_aux"], errors="coerce").fillna(0.5)
    out["decoy_reliability"] = pd.to_numeric(out["decoy_reliability"], errors="coerce").fillna(0.0)
    out["base_C_QK"] = 0.5 * out["C_counterfactual_rf"] + 0.5 * out["qk_no_anchor_gamma1"]
    out["base_C_QK_quantum"] = 0.5 * out["C_counterfactual_rf"] + 0.5 * out["qk_quantum_only_gamma1"]
    out["C_QK_disagreement"] = (out["C_counterfactual_rf"] - out["qk_no_anchor_gamma1"]).abs()
    out["decoy_agreement"] = (1 - (out["decoy_aux"] - out["base_C_QK"]).abs()).clip(0, 1)
    return out


META_COLS = [
    "C_counterfactual_rf",
    "qk_no_anchor_gamma1",
    "qk_quantum_only_gamma1",
    "D_structure_geometry_rf",
    "base_C_QK",
    "base_C_QK_quantum",
    "C_QK_disagreement",
    "decoy_aux",
    "decoy_reliability",
    "decoy_agreement",
    "peptide_len",
]


def run_locked() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    master = load_master().copy()
    folds = load_folds()
    experts = read_v0_expert_predictions()
    pred_rows = []
    weight_rows = []
    manifest_rows = []
    pattern_rows = []
    for (split, fold), _ in folds.groupby(["split_name", "fold_id"], sort=False):
        train_ids, test_ids = get_fold_ids(master, folds, split, fold)
        train = master[master["sample_id"].isin(train_ids)].copy()
        test = master[master["sample_id"].isin(test_ids)].copy()
        if len(test) < 2 or train["label"].nunique() < 2 or test["label"].nunique() < 2:
            continue
        p_aux, info, pattern = aux_score(train, test)
        aux_test = pd.DataFrame({"sample_id": test["sample_id"].values, "decoy_aux": p_aux, "decoy_reliability": info["pattern_reliability"]})
        aux_train = inner_oof_aux(train)
        train_wide = expert_wide_for_fold(experts, split, fold, train_ids, exclude_fold=True)
        test_wide = expert_wide_for_fold(experts, split, fold, test_ids, exclude_fold=False)
        train_meta = build_meta(train, train_wide, aux_train)
        test_meta = build_meta(test, test_wide, aux_test)

        base = test_meta["base_C_QK"].to_numpy(float)
        w_aux = np.clip(0.16 * test_meta["decoy_reliability"].to_numpy(float) * test_meta["decoy_agreement"].to_numpy(float), 0, 0.12)
        rule = np.clip((1 - w_aux) * base + w_aux * test_meta["decoy_aux"].to_numpy(float), 0, 1)
        for sid, y, s in zip(test["sample_id"], test["label"], rule):
            pred_rows.append({"split_name": split, "fold_id": fold, "model": "hard_decoy_rule_aux_C_QK_no_anchor", "sample_id": sid, "label": int(y), "score": float(s)})
        weight_rows.append(pd.DataFrame({"split_name": split, "fold_id": fold, "model": "hard_decoy_rule_aux_C_QK_no_anchor", "sample_id": test["sample_id"].values, "w_base_C_QK": 1 - w_aux, "w_hard_decoy_aux": w_aux, "pattern_reliability": test_meta["decoy_reliability"].to_numpy(float), "decoy_agreement": test_meta["decoy_agreement"].to_numpy(float)}))

        if len(train_meta) >= 24 and train_meta["label"].nunique() == 2:
            y = train_meta["label"].to_numpy(int)
            sw = np.where(y == 1, 1.25, 0.65)
            meta = focal_predict(train_meta[META_COLS], y, test_meta[META_COLS], sw)
            for sid, yy, ss in zip(test["sample_id"], test["label"], meta):
                pred_rows.append({"split_name": split, "fold_id": fold, "model": "hard_decoy_nested_meta_C_QK_aux", "sample_id": sid, "label": int(yy), "score": float(ss)})

        manifest_rows.append({"split_name": split, "fold_id": fold, **info})
        pattern["split_name"] = split
        pattern["fold_id"] = fold
        pattern_rows.append(pattern)
        print(f"[v1-hard-decoy] split={split} fold={fold} train={len(train)} test={len(test)} decoys={info['hard_decoy_n']} rel={info['pattern_reliability']:.3f}", flush=True)
    return pd.DataFrame(pred_rows), pd.concat(weight_rows, ignore_index=True), pd.DataFrame(manifest_rows), pd.concat(pattern_rows, ignore_index=True)


def run_source_stress() -> tuple[pd.DataFrame, pd.DataFrame]:
    master = load_master().copy()
    pool = master[master["study"].isin(["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"])].copy()
    rows = []
    manifests = []
    for heldout in ["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"]:
        train = pool[pool["study"] != heldout].copy()
        test = pool[pool["study"] == heldout].copy()
        if len(test) < 10 or train["label"].nunique() < 2 or test["label"].nunique() < 2:
            continue
        aux, info, _ = aux_score(train, test)
        w_aux = min(0.25, 0.18 * float(info["pattern_reliability"]))
        score = np.clip((1 - w_aux) * 0.5 + w_aux * aux, 0, 1)
        for sid, y, s in zip(test["sample_id"], test["label"], score):
            rows.append({"split_name": f"source_heldout_{heldout}", "fold_id": heldout, "model": "hard_decoy_sequence_only_source_stress", "sample_id": sid, "label": int(y), "score": float(s)})
        manifests.append({"split_name": f"source_heldout_{heldout}", "fold_id": heldout, **info})
        print(f"[v1-hard-decoy] source_heldout={heldout} train={len(train)} test={len(test)} decoys={info['hard_decoy_n']} rel={info['pattern_reliability']:.3f}", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(manifests)


def write_report(met: pd.DataFrame, manifest: pd.DataFrame, pattern: pd.DataFrame) -> None:
    hla = met[met["split_name"] == "hla_stratified_group_5fold"].sort_values("AUPRC", ascending=False)
    near = met[met["split_name"] == "near_peptide_cluster_holdout"].sort_values("AUPRC", ascending=False)
    src = met[met["split_name"].str.startswith("source_heldout", na=False)].sort_values("AUPRC", ascending=False)
    lines = [
        "# CROSS-Neo v1 Hard-Decoy/Focal Repair",
        "",
        "Hard decoys are now a bounded auxiliary contrastive confidence signal, not a standalone motif model. Public predictor scores are not used.",
        "",
        "## HLA-Stratified Locked",
        hla.to_markdown(index=False) if len(hla) else "Not available.",
        "",
        "## Near-Peptide Locked",
        near.to_markdown(index=False) if len(near) else "Not available.",
        "",
        "## Source-Stress",
        src.to_markdown(index=False) if len(src) else "Not available.",
        "",
        "## Decoy Manifest",
        manifest.head(80).to_markdown(index=False) if len(manifest) else "Not available.",
        "",
        "## Interpretation",
        "- Use the rule auxiliary model if it improves locked AUPRC/top-k without hurting near-peptide holdout.",
        "- Use sequence-only source stress only as a business triage diagnostic, not external validation.",
        "- Fold safety: decoys and motif statistics are created from train rows only.",
    ]
    if len(pattern):
        avg = pattern.groupby("position_1based", as_index=False).agg(top_positive_enriched_residues=("top_positive_enriched_residues", "first"), mean_positive_entropy=("positive_entropy", "mean"), mean_background_entropy=("background_entropy", "mean"))
        lines.extend(["", "## Fold-Averaged Positive Pattern", avg.to_markdown(index=False)])
    (V1 / "hard_decoy_focal_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ensure_v1_dirs()
    pred, weights, manifest, pattern = run_locked()
    src_pred, src_manifest = run_source_stress()
    pred = pd.concat([pred, src_pred], ignore_index=True, sort=False)
    manifest = pd.concat([manifest, src_manifest], ignore_index=True, sort=False)
    pred.to_csv(V1 / "hard_decoy_focal_predictions.tsv", sep="\t", index=False)
    weights.to_csv(V1 / "hard_decoy_focal_fold_weights.tsv", sep="\t", index=False)
    manifest.to_csv(V1 / "hard_decoy_generation_manifest.tsv", sep="\t", index=False)
    pattern.to_csv(V1 / "hard_decoy_positive_pattern_pwm.tsv", sep="\t", index=False)
    rows = []
    for (split, model), sub in pred.groupby(["split_name", "model"]):
        rows.append({"split_name": split, "model": model, **metrics(sub["label"].to_numpy(), sub["score"].to_numpy())})
    met = pd.DataFrame(rows).sort_values(["split_name", "AUPRC"], ascending=[True, False])
    met.to_csv(V1 / "hard_decoy_focal_metrics.tsv", sep="\t", index=False)
    write_report(met, manifest, pattern)
    print(f"[v1-hard-decoy] predictions={len(pred)} metrics={len(met)}")


if __name__ == "__main__":
    main()
