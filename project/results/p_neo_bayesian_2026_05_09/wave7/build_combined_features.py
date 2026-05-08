"""Wave 7 — Build combined feature matrix combining proven strengths.

Sources:
  1. wave2/structure_features.tsv  — 18-dim PWM + BLOSUM + physico (covers train + ITSNdb only)
  2. wave3_mhcflurry/predictions.tsv — presentation_score (only ITSNdb here; we re-run on full bundle)
  3. wave25/pwm_features.tsv — pwm_raw, pwm_pct (PWM scores)
  4. hla_pseudo.tsv — 34aa pseudo-seq → one-hot 680d → PCA(16)

The structure_features.tsv only covers ITSNdb (no Venus). We re-run the structure
feature pipeline on the FULL bundle so we have features for cross-source LOSO and Venus.

We also re-run MHCflurry on the FULL bundle to get presentation_score for everyone.

Output:
  combined_features.tsv  (peptide, HLA_norm, label, source, split, in_master, [features...])
"""
from __future__ import annotations
import os, sys, time, json, math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE7 = ROOT / "wave7"
BUNDLE = ROOT / "bundle.tsv"
HLA_PSEUDO = ROOT / "hla_pseudo.tsv"
WAVE2_STRUCT = ROOT / "wave2" / "structure_features.tsv"
WAVE25_PWM = ROOT / "wave25" / "pwm_features.tsv"
WAVE3_MHC = ROOT / "wave3_mhcflurry" / "predictions.tsv"
OUT = WAVE7 / "combined_features.tsv"
META = WAVE7 / "feature_meta.json"

t0 = time.time()
print("[load] bundle...", flush=True)
df = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})
df["pep_len"] = df["peptide"].str.len()

AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA_set = set(AA20)
mask_aa = df["peptide"].apply(lambda s: set(s).issubset(AA_set))
df = df[mask_aa].reset_index(drop=True).copy()
print(f"[load] {len(df)} rows post-AA filter", flush=True)

# ---------------- Structure feature recomputation (Wave 2 logic, adapted) ----------------
# Kyte-Doolittle hydrophobicity
KD = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}
CHG = {'D': -1, 'E': -1, 'K': 1, 'R': 1, 'H': 0.5}
# Chou-Fasman helix/sheet propensities
HEL = {'A':1.42,'R':0.98,'N':0.67,'D':1.01,'C':0.70,'Q':1.11,'E':1.51,'G':0.57,'H':1.00,'I':1.08,'L':1.21,'K':1.16,'M':1.45,'F':1.13,'P':0.57,'S':0.77,'T':0.83,'W':1.08,'Y':0.69,'V':1.06}
SHE = {'A':0.83,'R':0.93,'N':0.89,'D':0.54,'C':1.19,'Q':1.10,'E':0.37,'G':0.75,'H':0.87,'I':1.60,'L':1.30,'K':0.74,'M':1.05,'F':1.38,'P':0.55,'S':0.75,'T':1.19,'W':1.37,'Y':1.47,'V':1.70}
AROMATIC = set("FWY")

# BLOSUM62 self-similarity (we use match-to-self for max BLOSUM, plus mean)
# Compact lookup for max-BLOSUM-to-known-binder: we approximate by using
# BLOSUM62 self-matches to a mini consensus (just for feature consistency w/ Wave 2).
# Wave 2's structure_features.tsv defined blosum_max as max BLOSUM62 substitution
# score across positions (proxy for "binder-likeness"). We replicate using
# self-substitution score per position (i.e. diagonal). This is identical to
# Wave 2 since BLOSUM62 diagonal is the per-AA self-score.
BLOSUM_SELF = {'A':4,'R':5,'N':6,'D':6,'C':9,'Q':5,'E':5,'G':6,'H':8,'I':4,'L':4,'K':5,'M':5,'F':6,'P':7,'S':4,'T':5,'W':11,'Y':7,'V':4}

# Build a per-(HLA, position) PWM from training-positive peptides
print("[pwm] building PWM from train positives...", flush=True)
train_pos = df[(df["split"] == "train") & (df["label"] == 1)].copy()
print(f"[pwm] {len(train_pos)} train+ peptides for PWM", flush=True)

# Per-HLA per-position frequency (length-normalized: P1, P2, ..., POmega = last)
def build_pwm(group: pd.DataFrame, n_positions: int = 9) -> np.ndarray:
    """Return a (n_positions, 20) PWM log-odds matrix."""
    counts = np.zeros((n_positions, 20), dtype=float)
    aa2i = {a: i for i, a in enumerate(AA20)}
    for pep in group["peptide"]:
        L = len(pep)
        # Use first 4 positions (anchor-rich) and last 5 positions (POmega-1, POmega, ...)
        for i, aa in enumerate(pep):
            # Map to canonical 9-position frame:
            # positions 0..3 → 0..3
            # positions L-5..L-1 → 4..8
            if i < 4:
                pos = i
            elif i >= L - 5:
                pos = 4 + (i - (L - 5))
            else:
                continue
            if aa in aa2i:
                counts[pos, aa2i[aa]] += 1
    # Add pseudocount and convert to log-odds vs uniform
    counts += 1.0  # smoothing
    freqs = counts / counts.sum(axis=1, keepdims=True)
    bg = 1.0 / 20.0
    pwm = np.log(freqs / bg)
    return pwm

pwm_by_hla: dict[str, np.ndarray] = {}
for hla, grp in train_pos.groupby("HLA_norm"):
    if len(grp) < 5:
        continue
    pwm_by_hla[hla] = build_pwm(grp, n_positions=9)
print(f"[pwm] built {len(pwm_by_hla)} HLA-specific PWMs", flush=True)

# Global PWM fallback for unseen HLAs
GLOBAL_PWM = build_pwm(train_pos, n_positions=9)

def pep_to_pwm_score(pep: str, pwm: np.ndarray) -> tuple[float, float, float, float]:
    """Return (p2, pOmega, mean, sum) of PWM log-odds for the canonical 9-position framing."""
    L = len(pep)
    aa2i = {a: i for i, a in enumerate(AA20)}
    scores = []
    p2 = pOmega = 0.0
    for i, aa in enumerate(pep):
        if i < 4:
            pos = i
        elif i >= L - 5:
            pos = 4 + (i - (L - 5))
        else:
            continue
        if aa not in aa2i:
            continue
        s = pwm[pos, aa2i[aa]]
        scores.append(s)
        if pos == 1:
            p2 = s
        if pos == 8:
            pOmega = s
    if not scores:
        return (0.0, 0.0, 0.0, 0.0)
    return (p2, pOmega, float(np.mean(scores)), float(np.sum(scores)))

def physico(pep: str) -> dict:
    kd = [KD.get(a, 0.0) for a in pep]
    chg = [CHG.get(a, 0.0) for a in pep]
    hel = [HEL.get(a, 1.0) for a in pep]
    she = [SHE.get(a, 1.0) for a in pep]
    blo = [BLOSUM_SELF.get(a, 0) for a in pep]
    arom = sum(1 for a in pep if a in AROMATIC) / max(1, len(pep))
    return {
        'kd_mean': float(np.mean(kd)),
        'kd_std': float(np.std(kd)),
        'kd_max': float(np.max(kd)),
        'kd_min': float(np.min(kd)),
        'chg_mean': float(np.mean(chg)),
        'chg_sum': float(np.sum(chg)),
        'hel_mean': float(np.mean(hel)),
        'hel_std': float(np.std(hel)),
        'she_mean': float(np.mean(she)),
        'she_std': float(np.std(she)),
        'blosum_max': float(np.max(blo)),
        'blosum_mean': float(np.mean(blo)),
        'pep_len': float(len(pep)),
        'aromatic_frac': float(arom),
    }

print("[features] computing structural features for full bundle...", flush=True)
struct_rows = []
for _, r in df.iterrows():
    pwm = pwm_by_hla.get(r["HLA_norm"], GLOBAL_PWM)
    p2, pOm, pmean, psum = pep_to_pwm_score(r["peptide"], pwm)
    feat = {
        'pwm_p2': p2, 'pwm_pOmega': pOm, 'pwm_mean': pmean, 'pwm_sum': psum,
    }
    feat.update(physico(r["peptide"]))
    struct_rows.append(feat)
struct_df = pd.DataFrame(struct_rows)
struct_df.index = df.index
print(f"[features] structure {struct_df.shape}", flush=True)
STRUCT_COLS = list(struct_df.columns)

# ---------------- HLA pseudo-seq → one-hot → PCA 16d ----------------
print("[hla_pca] one-hot pseudo-seq + PCA(16)...", flush=True)
pseudo = pd.read_csv(HLA_PSEUDO, sep="\t").set_index("HLA_norm")
def onehot_pseudo(seq: str) -> np.ndarray:
    aa2i = {a: i for i, a in enumerate(AA20)}
    L = len(seq)
    M = np.zeros((L, 20), dtype=float)
    for i, a in enumerate(seq):
        if a in aa2i:
            M[i, aa2i[a]] = 1
    return M.flatten()

hla_unique = df["HLA_norm"].unique()
n_kept = 0
oh_rows = []
hla_to_idx = {}
for h in hla_unique:
    if h in pseudo.index:
        seq = pseudo.loc[h, "pseudo_seq"]
        if isinstance(seq, str) and len(seq) > 0:
            oh_rows.append(onehot_pseudo(seq))
            hla_to_idx[h] = n_kept
            n_kept += 1
oh_arr = np.array(oh_rows)  # (n_hla, 680)

from sklearn.decomposition import PCA
n_comp = min(16, oh_arr.shape[0] - 1, oh_arr.shape[1])
pca = PCA(n_components=n_comp, random_state=0)
oh_pca = pca.fit_transform(oh_arr)  # (n_hla, n_comp)
print(f"[hla_pca] {n_kept} HLAs → PCA dim {oh_pca.shape[1]} (var explained {pca.explained_variance_ratio_.sum():.3f})", flush=True)

hla_feat = np.zeros((len(df), n_comp), dtype=float)
for i, h in enumerate(df["HLA_norm"]):
    if h in hla_to_idx:
        hla_feat[i] = oh_pca[hla_to_idx[h]]
HLA_COLS = [f"hla_pc{i+1}" for i in range(n_comp)]

# ---------------- Wave25 PWM features (pwm_raw, pwm_pct) merge ----------------
print("[wave25] merging PWM features (raw+pct)...", flush=True)
w25 = pd.read_csv(WAVE25_PWM, sep="\t")
w25 = w25.rename(columns={"hla": "HLA_norm"})
w25 = w25[["peptide", "HLA_norm", "pwm_raw", "pwm_pct"]].drop_duplicates(["peptide", "HLA_norm"])

# ---------------- Re-run MHCflurry on full bundle (presentation_score) ----------------
print("[mhcflurry] loading predictor for full bundle...", flush=True)
df["mhc_filter_ok"] = df["pep_len"].between(8, 11)
n_drop_len = (~df["mhc_filter_ok"]).sum()
print(f"[mhcflurry] {n_drop_len} rows outside 8-11mer (will get NaN for MHCflurry)", flush=True)

from mhcflurry import Class1PresentationPredictor
predictor = Class1PresentationPredictor.load()
supported = set(predictor.supported_alleles)

mhc_score = np.full(len(df), np.nan, dtype=float)
mhc_aff = np.full(len(df), np.nan, dtype=float)
mhc_proc = np.full(len(df), np.nan, dtype=float)
mhc_pct = np.full(len(df), np.nan, dtype=float)

mask = df["mhc_filter_ok"] & df["HLA_norm"].isin(supported)
df_score = df[mask].copy().reset_index()  # keep original idx
print(f"[mhcflurry] scoring {len(df_score)} eligible rows...", flush=True)

t1 = time.time()
sample_names = [f"s{i}" for i in range(len(df_score))]
alleles_dict = {sn: [hla] for sn, hla in zip(sample_names, df_score["HLA_norm"].tolist())}
preds = predictor.predict(
    peptides=df_score["peptide"].tolist(),
    alleles=alleles_dict,
    sample_names=sample_names,
    verbose=0,
)
print(f"[mhcflurry] cols: {list(preds.columns)}", flush=True)
preds = preds.sort_values("peptide_num").reset_index(drop=True)

# Map back to df indices
orig_indices = df_score["index"].to_numpy()
mhc_score[orig_indices] = preds["presentation_score"].to_numpy()
if "processing_score" in preds.columns:
    mhc_proc[orig_indices] = preds["processing_score"].to_numpy()
if "affinity" in preds.columns:
    mhc_aff[orig_indices] = preds["affinity"].to_numpy()
if "affinity_percentile" in preds.columns:
    mhc_pct[orig_indices] = preds["affinity_percentile"].to_numpy()

print(f"[mhcflurry] done in {time.time()-t1:.1f}s", flush=True)

# Use logs for affinity (nM) where available; cap percentile rank at 100
mhc_aff_log = np.where(np.isfinite(mhc_aff), np.log1p(mhc_aff), np.nan)
mhc_pct_log = np.where(np.isfinite(mhc_pct), np.log1p(mhc_pct), np.nan)

MHC_COLS = ["mhc_present", "mhc_process", "mhc_aff_log", "mhc_pct_log"]

# ---------------- Assemble final feature matrix ----------------
out = df[["peptide", "HLA_norm", "label", "source", "split", "in_master", "pep_len"]].copy()
for col in STRUCT_COLS:
    out[col] = struct_df[col].to_numpy()
for j, col in enumerate(HLA_COLS):
    out[col] = hla_feat[:, j]
# wave25 pwm
out = out.merge(w25, on=["peptide", "HLA_norm"], how="left")
# MHCflurry
out["mhc_present"] = mhc_score
out["mhc_process"] = mhc_proc
out["mhc_aff_log"] = mhc_aff_log
out["mhc_pct_log"] = mhc_pct_log

# Fill NaNs in MHCflurry/wave25 columns with column median (per-source) so missing rows stay usable
for col in ["pwm_raw", "pwm_pct"] + MHC_COLS:
    if col in out.columns:
        med = out[col].median(skipna=True)
        n_nan = out[col].isna().sum()
        out[col] = out[col].fillna(med)
        if n_nan > 0:
            print(f"[fillna] {col}: {n_nan} NaN → median={med:.4f}", flush=True)

print(f"[done] combined feature matrix: {out.shape}", flush=True)
print(f"[done] feature columns: {list(out.columns)}", flush=True)

WAVE7.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, sep="\t", index=False)
print(f"[write] {OUT}", flush=True)

meta = {
    "n_rows": len(out),
    "structure_cols": STRUCT_COLS,
    "hla_cols": HLA_COLS,
    "mhc_cols": MHC_COLS,
    "wave25_cols": ["pwm_raw", "pwm_pct"],
    "feature_total_dim": len(STRUCT_COLS) + len(HLA_COLS) + len(MHC_COLS) + 2,
    "n_pwm_per_hla": len(pwm_by_hla),
    "pca_var_explained_sum": float(pca.explained_variance_ratio_.sum()),
    "elapsed_sec": time.time() - t0,
}
META.write_text(json.dumps(meta, indent=2))
print(f"[meta] feature_total_dim={meta['feature_total_dim']}", flush=True)
print(f"[done] elapsed {time.time()-t0:.1f}s", flush=True)
