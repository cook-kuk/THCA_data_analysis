#!/usr/bin/env python3
"""Fixed OOF late-fusion candidates for CROSS-Neo v0."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from cross_neo_v0_common import OUT, ensure_dirs


def main() -> None:
    ensure_dirs()
    pred = pd.read_csv(OUT / "oof_predictions.tsv", sep="\t")
    qk = pd.read_csv(OUT / "qk_fallback_predictions.tsv", sep="\t")
    rows = []
    # Prespecified equal-weight fusion is the audit-friendly candidate. The
    # non-0.5 weights are marked exploratory and must not be promoted.
    weights = [(0.5, "prespecified_equal_weight"), (0.25, "exploratory"), (0.75, "exploratory")]
    for split in sorted(set(pred.split_name) & set(qk.split_name)):
        base = pred[
            (pred.split_name == split)
            & (pred.feature_group == "C_counterfactual")
            & (pred.model == "rf_secondary")
        ][["split_name", "fold_id", "sample_id", "label", "score"]].rename(columns={"score": "p_counterfactual"})
        for branch in ["qk_no_anchor_gamma1", "qk_quantum_only_gamma1", "qk_clean_no_tcr_gamma1"]:
            q = qk[(qk.split_name == split) & (qk.branch == branch)][
                ["split_name", "fold_id", "sample_id", "score"]
            ].rename(columns={"score": "p_qk"})
            m = base.merge(q, on=["split_name", "fold_id", "sample_id"])
            if len(m) == 0 or m.label.nunique() < 2:
                continue
            for w, status in weights:
                score = w * m.p_counterfactual + (1 - w) * m.p_qk
                y = m.label.to_numpy(int)
                order = np.argsort(-score.to_numpy(float))
                rows.append(
                    {
                        "split_name": split,
                        "fusion": f"counterfactual_rf_plus_{branch}_w{w}",
                        "status": status,
                        "n": len(m),
                        "n_pos": int(y.sum()),
                        "prevalence": float(y.mean()),
                        "AUPRC": float(average_precision_score(y, score)),
                        "AUROC": float(roc_auc_score(y, score)),
                        "Brier": float(brier_score_loss(y, np.clip(score, 0, 1))),
                        "top5_precision": float(y[order[:5]].mean()),
                        "top10_precision": float(y[order[:10]].mean()),
                    }
                )
    out = pd.DataFrame(rows).sort_values(["AUPRC", "top10_precision"], ascending=False)
    out.to_csv(OUT / "late_fusion_oof_metrics.tsv", sep="\t", index=False)
    print(f"[late-fusion] rows={len(out)} best_AUPRC={out['AUPRC'].max():.3f}")


if __name__ == "__main__":
    main()
