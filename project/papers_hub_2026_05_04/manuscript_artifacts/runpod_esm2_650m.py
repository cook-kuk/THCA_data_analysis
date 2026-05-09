#!/usr/bin/env python3
"""RunPod-deploy script: ESM2-650M GPU embedding + LOSO benchmark.

Run on RunPod (RTX A4000 or better):
  python3 runpod_esm2_650m.py --benchmark benchmark_clean.tsv --max-rows 4500
"""
import argparse, time, json, warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np
import pandas as pd

KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,
      "H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,
      "T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
AA = "ACDEFGHIKLMNPQRSTVWY"

def biophys(peps):
    X = np.zeros((len(peps), 24), dtype=np.float32)
    for i, p in enumerate(peps):
        X[i,0] = sum(KD.get(a,0) for a in p)/max(len(p),1)
        X[i,1] = p.count("R")+p.count("K")+p.count("H")*0.1 - p.count("D") - p.count("E")
        X[i,2] = (p.count("F")+p.count("W")+p.count("Y"))/max(len(p),1)
        X[i,3] = len(p)
        for j, aa in enumerate(AA):
            X[i,4+j] = p.count(aa)/max(len(p),1)
    return X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", required=True)
    ap.add_argument("--max-rows", type=int, default=4500)
    ap.add_argument("--out", default="option_b_650m_results.tsv")
    ap.add_argument("--cache-emb", default="esm2_650m_embeddings.npy")
    args = ap.parse_args()

    df = pd.read_csv(args.benchmark, sep="\t", low_memory=False)
    df = df[df["safety"].isin(["HELD_OUT","EXTERNAL_TEST","CLINICAL_TEST"])].reset_index(drop=True)
    df = df[df["HLA"].notna() & (df["peptide"].str.len() <= 15)].reset_index(drop=True)
    if len(df) > args.max_rows:
        df = df.groupby("source", group_keys=False).apply(
            lambda g: g.sample(min(len(g), int(args.max_rows * len(g) / len(df))), random_state=42)
        ).reset_index(drop=True)
    y = df["label"].values
    print(f"benchmark: {len(df):,}")

    # === Run ESM2-650M on GPU ===
    print("loading ESM2-650M...")
    import torch, esm
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"  device: {device}")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()
    print(f"  params: {sum(p.numel() for p in model.parameters())/1e6:.1f}M")

    peps = df["peptide"].tolist()
    embed_dim = 1280  # ESM2-650M
    E = np.zeros((len(peps), embed_dim), dtype=np.float32)
    BATCH = 32
    t0 = time.time()
    with torch.no_grad():
        for i in range(0, len(peps), BATCH):
            batch = [(f"p{k}", peps[i+k]) for k in range(min(BATCH, len(peps)-i))]
            _, _, tokens = batch_converter(batch)
            tokens = tokens.to(device)
            out = model(tokens, repr_layers=[33])
            reps = out["representations"][33]
            for k in range(len(batch)):
                pep_len = len(peps[i+k])
                E[i+k] = reps[k, 1:pep_len+1].mean(dim=0).cpu().numpy()
            if (i // BATCH) % 20 == 0:
                done = i + len(batch)
                el = time.time() - t0
                eta = el * (len(peps) - done) / max(done, 1)
                print(f"  {done}/{len(peps)}  el={el:.0f}s  ETA={eta:.0f}s")
    print(f"  done in {time.time()-t0:.0f}s")
    np.save(args.cache_emb, E)
    print(f"  saved {args.cache_emb}")

    # === LOSO benchmark ===
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import roc_auc_score

    bp = biophys(peps)
    feature_sets = {
        "biophys": bp,
        "esm2_650m": E,
        "esm2_650m_plus_biophys": np.hstack([E, bp]),
    }

    rows = []
    for fname, X in feature_sets.items():
        for mt in ["logreg", "rf", "mlp_2layer"]:
            scale = mt in ("logreg", "mlp_2layer")
            aurocs = []; per = []
            for s in pd.Series(df["source"].values).unique():
                ti = np.where(df["source"].values != s)[0]
                vi = np.where(df["source"].values == s)[0]
                if len(vi) < 50 or len(ti) < 200: continue
                if y[ti].sum() < 30 or (1-y[ti]).sum() < 30: continue
                Xtr, Xte = X[ti], X[vi]
                if scale:
                    sc = StandardScaler().fit(Xtr)
                    Xtr = sc.transform(Xtr); Xte = sc.transform(Xte)
                if mt == "logreg":
                    clf = LogisticRegression(C=1.0, max_iter=300, class_weight="balanced", n_jobs=-1)
                elif mt == "rf":
                    clf = RandomForestClassifier(n_estimators=100, max_depth=8, n_jobs=-1, random_state=42, class_weight="balanced")
                elif mt == "mlp_2layer":
                    clf = MLPClassifier(hidden_layer_sizes=(128, 64), max_iter=200, random_state=42, early_stopping=True)
                clf.fit(Xtr, y[ti])
                au = roc_auc_score(y[vi], clf.predict_proba(Xte)[:,1])
                aurocs.append(au)
                per.append({"source":s, "AUROC":float(au)})
            if aurocs:
                rows.append({"feature_set":fname, "model":mt,
                             "AUROC_mean":float(np.mean(aurocs)),
                             "AUROC_std":float(np.std(aurocs)),
                             "n_sources":len(per),
                             "per_source":json.dumps(per)})
                print(f"  {fname:25s} {mt:14s} LOSO AUROC={np.mean(aurocs):.3f} ± {np.std(aurocs):.3f}")

    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"saved {args.out}")
    print("\n=== HEADLINE ===")
    if rows:
        best = max(rows, key=lambda r: r["AUROC_mean"])
        print(f"Best: {best['feature_set']} + {best['model']} = {best['AUROC_mean']:.3f}")


if __name__ == "__main__":
    main()
