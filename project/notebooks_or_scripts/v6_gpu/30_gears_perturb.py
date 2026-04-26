#!/usr/bin/env python3
"""v6 Wave 2 — GEARS in-silico KO of 8 druggable targets + 28 pairs + BRAF-like-B rescue.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
from itertools import combinations

DRUGGABLE = ["CYP1B1","TACSTD2","TMPRSS4","PLEKHA6","LDLR","GABRB2","B3GNT3","PTPRE"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", default="/opt/thyroid-dash/project/results/v6_scrna")
    ap.add_argument("--combos", action="store_true", default=True)
    ap.add_argument("--pretrained", default="gears-norman")
    ap.add_argument("--model_cache", default="/data/thca/scrna/models/gears-norman")
    args = ap.parse_args()
    out = Path(args.output_dir); (out/"perturbation").mkdir(parents=True, exist_ok=True)

    try:
        import torch, scanpy as sc, numpy as np, pandas as pd
        from gears import GEARS, PertData
    except Exception as e:
        print(f"FATAL: gears install required (pip install git+https://github.com/snap-stanford/GEARS.git)\n{e}",
              file=sys.stderr); sys.exit(2)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[gears] device={device}")

    pert_data = PertData(args.model_cache)
    pert_data.load(data_path=args.model_cache)
    pert_data.prepare_split(split="simulation", seed=1)
    model = GEARS(pert_data, device=device)
    model.load_pretrained(args.model_cache)

    # Single-target KO
    rows = []
    for t in DRUGGABLE:
        try:
            pred = model.predict([[t]])
            mean_shift = float(np.mean(np.abs(pred)))
            rows.append({"target":t, "n_genes_changed": int((np.abs(pred)>0.1).sum()),
                         "mean_abs_shift": mean_shift})
        except Exception as e:
            rows.append({"target":t, "n_genes_changed":0, "mean_abs_shift":None, "error":str(e)})
    pd.DataFrame(rows).to_csv(out/"perturbation/gears_single_ko.tsv", sep="\t", index=False)

    # Combinations
    if args.combos:
        rows2 = []
        for a, b in combinations(DRUGGABLE, 2):
            try:
                p_combo = model.predict([[a, b]])
                p_a = model.predict([[a]])
                p_b = model.predict([[b]])
                synergy = float(np.mean(np.abs(p_combo - (p_a + p_b))))
                rows2.append({"target_A":a, "target_B":b, "synergy":synergy})
            except Exception as e:
                rows2.append({"target_A":a, "target_B":b, "synergy":None, "error":str(e)})
        pd.DataFrame(rows2).to_csv(out/"perturbation/gears_combo_ko.tsv", sep="\t", index=False)

    # BRAF-like-B rescue placeholder — depends on cell-state annotation
    rescue_rows = [{"target":t, "rescue_score":None,
                    "note":"compute after cell-state annotation; cosine shift toward normal-thyrocyte centroid"}
                   for t in DRUGGABLE]
    pd.DataFrame(rescue_rows).to_csv(out/"perturbation/gears_braflike_b_rescue.tsv", sep="\t", index=False)

    Path(f"{args.output_dir}/.stamp_30_gears_perturb").touch()
    print("[gears] DONE")

if __name__ == "__main__":
    main()
