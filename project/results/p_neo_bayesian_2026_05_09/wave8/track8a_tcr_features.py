"""Wave 8 — Track 8A: TCR engagement proxy features.

Builds a TCR-engagement reference set from local VDJdb + McPAS-TCR (cancer-tagged + all-human),
then computes for each test peptide x HLA:
  - tcr_motif_score: max k-mer Jaccard (k=3,4,5 averaged) similarity to any reference peptide for the SAME HLA
                     (falls back to global if HLA missing)
  - tcr_motif_count: # reference peptides within Hamming distance <=2 (length-matched), same-HLA
  - tcr_class_score: per HLA-supertype (HLA-A02, HLA-A24, HLA-B07, ...) reference engagement frequency
                     (log1p of count of references for that supertype, normalized 0..1)

Reference data (local):
  - /data/neoantigen_vaccine_hub/data_raw/tcr/vdjdb.slim.txt   (~75K rows)
  - /data/neoantigen_vaccine_hub/data_raw/tcr/mcpas/McPAS-TCR.csv

Output: wave8/tcr_features.tsv
"""
from __future__ import annotations
import os, sys, time, json, math, re
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE8 = ROOT / "wave8"
BUNDLE = ROOT / "bundle.tsv"
VDJDB = Path("/data/neoantigen_vaccine_hub/data_raw/tcr/vdjdb.slim.txt")
MCPAS = Path("/data/neoantigen_vaccine_hub/data_raw/tcr/mcpas/McPAS-TCR.csv")

OUT = WAVE8 / "tcr_features.tsv"
META = WAVE8 / "tcr_features_meta.json"

AA20 = set("ACDEFGHIKLMNPQRSTVWY")

def norm_hla(h: str) -> str:
    if not isinstance(h, str):
        return ""
    h = h.strip().upper().replace(" ", "")
    # Common patterns: HLA-A*02:01, HLA-A02:01, HLA-A2, A*02:01, HLA-A0201
    m = re.match(r"^(?:HLA-)?([ABC])\*?(\d{1,2}):?(\d{2})", h)
    if m:
        return f"HLA-{m.group(1)}*{int(m.group(2)):02d}:{m.group(3)}"
    m2 = re.match(r"^(?:HLA-)?([ABC])\*?(\d{1,2})$", h)  # supertype only
    if m2:
        return f"HLA-{m2.group(1)}*{int(m2.group(2)):02d}:01"
    # E.g. HLA-A02:01 (no asterisk)
    m3 = re.match(r"^HLA-([ABC])(\d{2}):(\d{2})$", h)
    if m3:
        return f"HLA-{m3.group(1)}*{m3.group(2)}:{m3.group(3)}"
    return h

def supertype(hla: str) -> str:
    """Return supertype prefix: 'HLA-A02', 'HLA-B07' etc."""
    m = re.match(r"^HLA-([ABC])\*(\d{2}):", hla)
    if m:
        return f"HLA-{m.group(1)}{m.group(2)}"
    return ""

def kmers(pep: str, k: int):
    if len(pep) < k:
        return set()
    return {pep[i:i+k] for i in range(len(pep) - k + 1)}

def jaccard(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    inter = len(a & b)
    uni = len(a | b)
    return inter / uni if uni > 0 else 0.0

def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return 999
    return sum(x != y for x, y in zip(a, b))


def main():
    t0 = time.time()
    WAVE8.mkdir(parents=True, exist_ok=True)

    # ---- Load reference TCR-engaged peptides ----
    print("[ref] loading VDJdb...", flush=True)
    vdj = pd.read_csv(VDJDB, sep="\t", low_memory=False)
    # Use human, MHC-I only, length 8-12
    vdj_h = vdj[(vdj["species"] == "HomoSapiens") & (vdj["mhc.class"] == "MHCI")].copy()
    vdj_h["peptide"] = vdj_h["antigen.epitope"].astype(str)
    # Compose HLA from mhc.a (HLA-A*02:01 already-formatted in VDJdb)
    vdj_h["HLA_norm"] = vdj_h["mhc.a"].astype(str).map(norm_hla)
    vdj_h = vdj_h[
        vdj_h["peptide"].apply(lambda s: 8 <= len(s) <= 12 and set(s).issubset(AA20))
    ].copy()
    vdj_ref = vdj_h[["peptide", "HLA_norm", "antigen.species"]].drop_duplicates()
    print(f"[ref] VDJdb: {len(vdj_ref)} unique (peptide, HLA) human MHC-I", flush=True)

    print("[ref] loading McPAS-TCR...", flush=True)
    mc = pd.read_csv(MCPAS, low_memory=False, encoding_errors="replace")
    mc = mc[mc["Species"].astype(str).str.contains("Human", na=False, case=False)].copy()
    mc["peptide"] = mc["Epitope.peptide"].astype(str)
    mc["HLA_norm"] = mc["MHC"].astype(str).map(norm_hla)
    mc = mc[
        mc["peptide"].apply(lambda s: 8 <= len(s) <= 12 and set(s).issubset(AA20))
    ].copy()
    mc_ref = mc[["peptide", "HLA_norm", "Pathology"]].rename(columns={"Pathology": "antigen.species"})
    mc_ref = mc_ref.drop_duplicates()
    print(f"[ref] McPAS: {len(mc_ref)} unique (peptide, HLA) human MHC-I", flush=True)

    ref = pd.concat([vdj_ref, mc_ref], ignore_index=True).drop_duplicates(["peptide", "HLA_norm"])
    print(f"[ref] combined unique reference: {len(ref)}", flush=True)

    # Per-HLA reference set
    ref_by_hla: dict[str, list[str]] = defaultdict(list)
    for pep, hla in zip(ref["peptide"], ref["HLA_norm"]):
        ref_by_hla[hla].append(pep)
    # Global reference
    ref_all = list(set(ref["peptide"].tolist()))
    print(f"[ref] {len(ref_by_hla)} unique reference HLAs; global pool {len(ref_all)} peptides", flush=True)

    # Per-supertype counts (for tcr_class_score)
    super_counts: Counter = Counter()
    for hla in ref_by_hla:
        st = supertype(hla)
        if st:
            super_counts[st] += len(ref_by_hla[hla])
    max_super = max(super_counts.values()) if super_counts else 1
    super_score = {st: math.log1p(n) / math.log1p(max_super) for st, n in super_counts.items()}
    print(f"[ref] supertype score range: {min(super_score.values()):.3f}..{max(super_score.values()):.3f}", flush=True)

    # Pre-build kmer sets per reference peptide for the union of k=3,4,5 (Jaccard average)
    def multi_k(pep: str) -> tuple[set, set, set]:
        return kmers(pep, 3), kmers(pep, 4), kmers(pep, 5)

    print("[ref] precomputing kmer sets...", flush=True)
    ref_kmers_by_hla: dict[str, list[tuple[set, set, set]]] = {}
    for hla, peps in ref_by_hla.items():
        ref_kmers_by_hla[hla] = [multi_k(p) for p in peps]
    ref_kmers_global = [multi_k(p) for p in ref_all]

    # ---- Score test bundle ----
    print("[score] loading bundle...", flush=True)
    df = pd.read_csv(BUNDLE, sep="\t", low_memory=False)
    df = df[df["peptide"].apply(lambda s: isinstance(s, str) and set(s).issubset(AA20))].copy()
    df["HLA_norm"] = df["HLA_norm"].fillna("").astype(str)
    print(f"[score] {len(df)} bundle rows to score", flush=True)

    out_motif_score = np.zeros(len(df), dtype=float)
    out_motif_count = np.zeros(len(df), dtype=int)
    out_class_score = np.zeros(len(df), dtype=float)

    # cache: per (HLA_norm, peptide_len) → list of length-matched ref peptides for Hamming check
    ref_by_hla_len: dict[tuple, list[str]] = {}
    for hla, peps in ref_by_hla.items():
        for p in peps:
            ref_by_hla_len.setdefault((hla, len(p)), []).append(p)
    # Same for global
    ref_by_len: dict[int, list[str]] = defaultdict(list)
    for p in ref_all:
        ref_by_len[len(p)].append(p)

    t_score = time.time()
    n_total = len(df)
    for i, (pep, hla) in enumerate(zip(df["peptide"].tolist(), df["HLA_norm"].tolist())):
        if i % 500 == 0:
            print(f"[score] {i}/{n_total} elapsed {time.time()-t_score:.1f}s", flush=True)
        pk3, pk4, pk5 = multi_k(pep)
        # Try same-HLA ref pool; fallback to global
        ref_kmer_pool = ref_kmers_by_hla.get(hla)
        ref_pep_pool_len = ref_by_hla_len.get((hla, len(pep)))
        used_global = False
        if not ref_kmer_pool:
            ref_kmer_pool = ref_kmers_global
            ref_pep_pool_len = ref_by_len.get(len(pep), [])
            used_global = True

        # Motif score: average jaccard at k=3,4,5, max over reference set
        max_avg = 0.0
        for rk3, rk4, rk5 in ref_kmer_pool:
            j3 = jaccard(pk3, rk3)
            j4 = jaccard(pk4, rk4)
            j5 = jaccard(pk5, rk5)
            avg = (j3 + j4 + j5) / 3.0
            if avg > max_avg:
                max_avg = avg
        # If we fell back to global, multiply by 0.5 to discount cross-HLA evidence
        out_motif_score[i] = max_avg * (0.5 if used_global else 1.0)

        # Hamming<=2 count, length-matched
        cnt = 0
        if ref_pep_pool_len:
            for rp in ref_pep_pool_len:
                if hamming(pep, rp) <= 2:
                    cnt += 1
        out_motif_count[i] = cnt

        # Supertype class score
        st = supertype(hla)
        out_class_score[i] = super_score.get(st, 0.0)

    print(f"[score] done in {time.time()-t_score:.1f}s", flush=True)

    # Log-transform count
    out_motif_count_log = np.log1p(out_motif_count.astype(float))

    out_df = df[["peptide", "HLA_norm"]].copy()
    out_df = out_df.rename(columns={"HLA_norm": "hla"})
    out_df["tcr_motif_score"] = out_motif_score
    out_df["tcr_motif_count_log"] = out_motif_count_log
    out_df["tcr_class_score"] = out_class_score
    out_df.to_csv(OUT, sep="\t", index=False)

    meta = {
        "n_rows": int(len(out_df)),
        "n_ref_unique": int(len(ref)),
        "n_ref_hla": int(len(ref_by_hla)),
        "n_supertypes": int(len(super_score)),
        "supertype_top": dict(super_counts.most_common(8)),
        "fraction_motif_score_gt_0": float((out_motif_score > 0).mean()),
        "elapsed_sec": time.time() - t0,
    }
    META.write_text(json.dumps(meta, indent=2))
    print(f"[done] wrote {OUT}", flush=True)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
