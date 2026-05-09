"""Reconstruct wave4a_results.tsv from saved predictions_<method>.tsv files.

Avoids re-running already-completed methods. Computes AUROC + bootstrap CI
fresh from prediction files.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss


def ece(y, p, n_bins=15):
    y, p = np.asarray(y), np.asarray(p)
    bins = np.linspace(0, 1, n_bins + 1)
    ece_v = 0.0
    n = len(y)
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (p >= lo) & (p < hi if hi < 1 else p <= hi)
        if mask.sum() == 0:
            continue
        ece_v += (mask.sum() / n) * abs(p[mask].mean() - y[mask].mean())
    return float(ece_v)


def auroc_ci(y, p, n_boot=1000, seed=0):
    y = np.asarray(y); p = np.asarray(p)
    if len(set(y)) < 2 or len(y) < 5:
        return None, None, None
    auc = float(roc_auc_score(y, p))
    rng = np.random.default_rng(seed)
    n = len(y)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(set(y[idx])) < 2:
            continue
        boots.append(roc_auc_score(y[idx], p[idx]))
    if not boots:
        return auc, None, None
    return auc, float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))


def metrics_block(y, p):
    y = np.asarray(y); p = np.asarray(p)
    if len(set(y)) < 2 or len(y) < 5:
        return None
    auc, lo, hi = auroc_ci(y, p)
    return {
        "AUROC": auc,
        "AUROC_lo95": lo,
        "AUROC_hi95": hi,
        "AUPRC": float(average_precision_score(y, p)),
        "Brier": float(brier_score_loss(y, p)),
        "ECE": ece(y, p),
        "n": int(len(y)),
        "n_pos": int(np.sum(y)),
    }


def process_method(method, pred_path, bundle_path):
    """Read predictions for a method and reconstruct AUROC for each testset."""
    if not Path(pred_path).exists():
        print(f"  no predictions for {method} at {pred_path}")
        return []
    pred = pd.read_csv(pred_path, sep="\t")
    bundle = pd.read_csv(bundle_path, sep="\t")
    bundle["HLA_norm"] = bundle["HLA_norm"].fillna("")
    if "in_master" not in bundle.columns:
        bundle["in_master"] = False
    in_master_lookup = bundle.set_index(["peptide", "HLA_norm"])["in_master"].to_dict()

    rows = []
    # In-domain 5-fold
    for fold in range(5):
        sub = pred[pred["split"] == f"in_domain_fold{fold}"]
        if len(sub) == 0:
            continue
        m = metrics_block(sub["label"].values, sub["pred_mean"].values)
        if m:
            rows.append({"method": method, "testset": "in_domain_5fold", "fold": fold,
                         "in_master": "NA", **m})
    # LOSO
    for split_val in pred["split"].unique():
        if not isinstance(split_val, str) or not split_val.startswith("loso_"):
            continue
        sub = pred[pred["split"] == split_val]
        if len(sub) == 0:
            continue
        m = metrics_block(sub["label"].values, sub["pred_mean"].values)
        if m:
            rows.append({"method": method, "testset": split_val, "fold": "NA",
                         "in_master": "NA", **m})
    # ITSNdb (split = ext_itsndb_main / ext_itsndb_val)
    itsn = pred[pred["split"].astype(str).str.startswith("ext_itsndb")]
    if len(itsn) > 0:
        # Add in_master if missing
        if "in_master" not in itsn.columns or itsn["in_master"].isna().any():
            itsn = itsn.copy()
            itsn["in_master"] = itsn.apply(
                lambda r: bool(in_master_lookup.get((r["peptide"], r["HLA_norm"]), False)), axis=1)
        else:
            itsn = itsn.copy()
            itsn["in_master"] = itsn["in_master"].astype(bool)
        for sub_name, mask, im in [
            ("ITSNdb_main", itsn["split"] == "ext_itsndb_main", "mixed"),
            ("ITSNdb_Val",  itsn["split"] == "ext_itsndb_val", "mixed"),
            ("ITSNdb_combined", np.ones(len(itsn), dtype=bool), "mixed"),
            ("ITSNdb_no_overlap", ~itsn["in_master"].values, False),
            ("ITSNdb_in_master", itsn["in_master"].values, True),
        ]:
            sub = itsn[mask]
            m = metrics_block(sub["label"].values, sub["pred_mean"].values)
            if m:
                rows.append({"method": method, "testset": sub_name, "fold": "NA",
                             "in_master": im, **m})
    # Venus  — note: predictions_*.tsv stores per-(window, HLA) rows for venus,
    # NOT aggregated. We need to aggregate top10_mean per protein here.
    venus = pred[pred["split"].astype(str).str.startswith("ext_venus")]
    if len(venus) > 0:
        # Group by source (which holds split) and protein-id; we don't have
        # protein_id in the prediction TSV. Skip Venus reconstruction unless
        # prediction file has protein_id.
        if "protein_id" in venus.columns:
            for split_name in ["ext_venus_test", "ext_venus_valid"]:
                sub = venus[venus["split"] == split_name]
                if len(sub) == 0:
                    continue
                agg_rows = []
                for (pid, label), g in sub.groupby(["protein_id", "label"]):
                    sc = g["pred_mean"].values
                    agg_rows.append({"protein_id": pid, "label": int(label),
                                     "top10_mean": float(np.sort(sc)[::-1][:10].mean())})
                vdf = pd.DataFrame(agg_rows)
                if len(vdf) >= 5 and vdf["label"].nunique() == 2:
                    m = metrics_block(vdf["label"].values, vdf["top10_mean"].values)
                    if m:
                        rows.append({"method": method, "testset": f"venus_{split_name}_top10",
                                     "fold": "NA", "in_master": "NA", **m})
    return rows


def parse_log_venus(log_path, method):
    """Fallback: scrape Venus AUROC from run.log if prediction file lacks it.

    Parses lines like:
      [groupdro]     ext_venus_test top10_mean: n=78 AUROC=0.802 [0.693,0.893]
    """
    rows = []
    if not Path(log_path).exists():
        return rows
    text = Path(log_path).read_text()
    for line in text.splitlines():
        if f"[{method}]" not in line or "venus" not in line.lower():
            continue
        # Extract split, AUROC, lo, hi, n
        # crude regex via str.find
        try:
            sp = "ext_venus_test" if "ext_venus_test" in line else (
                "ext_venus_valid" if "ext_venus_valid" in line else None)
            if not sp:
                continue
            n = int(line.split("n=")[1].split()[0])
            auc = float(line.split("AUROC=")[1].split()[0])
            # Find the LAST [...] (the CI bracket — the first [ is [method])
            ci = line.rsplit("[", 1)[1].split("]")[0]
            lo, hi = [float(x) for x in ci.split(",")]
            rows.append({"method": method, "testset": f"venus_{sp}_top10", "fold": "NA",
                         "in_master": "NA", "AUROC": auc, "AUROC_lo95": lo, "AUROC_hi95": hi,
                         "AUPRC": None, "Brier": None, "ECE": None, "n": n, "n_pos": None})
        except Exception:
            pass
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred_dir", default=".")
    ap.add_argument("--bundle", default="../bundle.tsv")
    ap.add_argument("--log", default="run.log")
    ap.add_argument("--out", default="wave4a_results.tsv")
    ap.add_argument("--methods", default="groupdro,miro,mole,lora")
    args = ap.parse_args()

    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    all_rows = []
    for m in methods:
        pred_path = Path(args.pred_dir) / f"predictions_{m}.tsv"
        rows = process_method(m, pred_path, args.bundle)
        # Add Venus from log as fallback
        existing_testsets = {r["testset"] for r in rows}
        for r in parse_log_venus(args.log, m):
            if r["testset"] not in existing_testsets:
                rows.append(r)
        all_rows.extend(rows)
        print(f"  {m}: {len(rows)} rows reconstructed")

    df = pd.DataFrame(all_rows)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"saved {len(df)} rows to {args.out}")
    # Print headline
    print("\n--- ITSNdb_no_overlap headline ---")
    sub = df[df["testset"] == "ITSNdb_no_overlap"]
    print(sub[["method", "AUROC", "AUROC_lo95", "AUROC_hi95", "n", "n_pos"]].to_string(index=False))


if __name__ == "__main__":
    main()
