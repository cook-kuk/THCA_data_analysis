"""Wave 8 — Track 8C: Self-similarity to UniProt human SwissProt proteome.

Memory-efficient version: builds plain set of human k-mers per length;
for each query peptide, enumerates its Hamming-ball-2 neighbors and looks up
in the set (set membership is O(1)).

For L=12: |Hamming-1| = 19*L = 228 ; |Hamming-2| = C(L,2)*19^2 = 66*361 ≈ 24K.
2871 peptides * 24K = ~70M lookups, fast.

Features:
  - self_exact_match : 1/0 if peptide is in human k-mer set
  - self_hamming1_count_log : log1p(# human k-mers within Hamming distance 1)
  - self_hamming2_count_log : log1p(# within Hamming distance 2; includes h0,h1,h2)
  - self_blosum_max : max BLOSUM62 substitution-sum over Hamming-ball neighbors
                      (falls back to self-blosum if no neighbor found within ball-2).

Restrict to peptides of length 8..12 (longer rows, e.g. Venus full-protein sequences,
get NaN-default features that downstream eval fills with median).
"""
from __future__ import annotations
import os, sys, time, json, math, gzip, urllib.request, itertools
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE8 = ROOT / "wave8"
BUNDLE = ROOT / "bundle.tsv"
FASTA = WAVE8 / "human_swissprot.fa"
OUT = WAVE8 / "self_similarity_features.tsv"
META = WAVE8 / "self_similarity_meta.json"

AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA20_set = set(AA20)

# BLOSUM62
BLOSUM62_AA = "ARNDCQEGHILKMFPSTWYV"
BLOSUM62_RAW = """4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4"""
_blosum_arr = np.array([[int(x) for x in row.split()] for row in BLOSUM62_RAW.strip().split("\n")])
_blosum_idx = {a: i for i, a in enumerate(BLOSUM62_AA)}

def blosum_score(a: str, b: str) -> int:
    s = 0
    for x, y in zip(a, b):
        ix = _blosum_idx.get(x); iy = _blosum_idx.get(y)
        if ix is None or iy is None:
            continue
        s += int(_blosum_arr[ix, iy])
    return s


def fetch_human_proteome():
    if FASTA.exists() and FASTA.stat().st_size > 5_000_000:
        print(f"[fetch] using cached {FASTA} ({FASTA.stat().st_size/1e6:.1f} MB)", flush=True)
        return
    url = ("https://rest.uniprot.org/uniprotkb/stream"
           "?compressed=true&format=fasta"
           "&query=organism_id%3A9606+AND+reviewed%3Atrue")
    print(f"[fetch] downloading {url}...", flush=True)
    tmp_gz = FASTA.with_suffix(".fa.gz")
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp_gz, "wb") as fh:
        while True:
            chunk = r.read(1 << 16)
            if not chunk:
                break
            fh.write(chunk)
    print("[fetch] decompressing...", flush=True)
    with gzip.open(tmp_gz, "rb") as gz, open(FASTA, "wb") as out:
        while True:
            chunk = gz.read(1 << 16)
            if not chunk:
                break
            out.write(chunk)
    tmp_gz.unlink(missing_ok=True)
    print(f"[fetch] wrote {FASTA} ({FASTA.stat().st_size/1e6:.1f} MB)", flush=True)


def parse_proteome():
    seqs = []
    cur = []
    with open(FASTA) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))
    cleaned = []
    for s in seqs:
        s2 = "".join(c if c in AA20_set else "X" for c in s.upper())
        cleaned.append(s2)
    return cleaned


def build_kmer_sets(proteins: list, lengths: list):
    kmer: dict = {L: set() for L in lengths}
    for prot in proteins:
        n = len(prot)
        for L in lengths:
            if n < L:
                continue
            for i in range(n - L + 1):
                w = prot[i:i+L]
                if "X" in w:
                    continue
                kmer[L].add(w)
    return kmer


def hamming_ball2(query: str, kmer_set: set):
    """Enumerate Hamming-ball-2 neighbors of query and look up in kmer_set.
    Return dict {neighbor: hamming_distance}.
    Memory: O(|hits|), not O(neighbors).
    """
    L = len(query)
    hits = {}
    # Hamming 0
    if query in kmer_set:
        hits[query] = 0
    # Hamming 1
    q = list(query)
    for i in range(L):
        orig = q[i]
        for aa in AA20:
            if aa == orig:
                continue
            q[i] = aa
            cand = "".join(q)
            if cand in kmer_set and cand not in hits:
                hits[cand] = 1
        q[i] = orig
    # Hamming 2
    for i in range(L):
        orig_i = q[i]
        for aa_i in AA20:
            if aa_i == orig_i:
                continue
            q[i] = aa_i
            for j in range(i + 1, L):
                orig_j = q[j]
                for aa_j in AA20:
                    if aa_j == orig_j:
                        continue
                    q[j] = aa_j
                    cand = "".join(q)
                    if cand in kmer_set and cand not in hits:
                        hits[cand] = 2
                q[j] = orig_j
            q[i] = orig_i
    return hits


def main():
    t0 = time.time()
    WAVE8.mkdir(parents=True, exist_ok=True)

    fetch_human_proteome()
    print("[parse] reading proteome...", flush=True)
    proteins = parse_proteome()
    print(f"[parse] {len(proteins)} proteins", flush=True)

    df = pd.read_csv(BUNDLE, sep="\t", low_memory=False)
    df = df[df["peptide"].apply(lambda s: isinstance(s, str))].copy()
    df["HLA_norm"] = df["HLA_norm"].fillna("").astype(str)
    df["pep_len"] = df["peptide"].str.len()

    LENS_ALL = sorted(df["pep_len"].unique().tolist())
    LENS = [L for L in LENS_ALL if 8 <= L <= 12]
    print(f"[bundle] indexing lengths {LENS} (excluding {[L for L in LENS_ALL if L not in LENS]})", flush=True)

    print("[idx] building human k-mer sets...", flush=True)
    kmer = build_kmer_sets(proteins, LENS)
    for L in LENS:
        print(f"[idx] L={L}: {len(kmer[L])} unique human k-mers", flush=True)

    n = len(df)
    out_exact = np.zeros(n, dtype=int)
    out_h1 = np.full(n, np.nan, dtype=float)
    out_h2 = np.full(n, np.nan, dtype=float)
    out_blosum = np.full(n, np.nan, dtype=float)

    t_score = time.time()
    print("[score] scoring peptides...", flush=True)
    for i, (pep, L) in enumerate(zip(df["peptide"].tolist(), df["pep_len"].tolist())):
        if i % 250 == 0:
            print(f"[score] {i}/{n} elapsed {time.time()-t_score:.1f}s", flush=True)
        if L not in kmer or set(pep) - AA20_set:
            continue
        hits = hamming_ball2(pep, kmer[L])
        if pep in kmer[L]:
            out_exact[i] = 1
        h1 = sum(1 for v in hits.values() if v == 1)
        h2 = sum(1 for v in hits.values() if v == 2)
        out_h1[i] = h1
        out_h2[i] = h2
        if hits:
            best = max(blosum_score(pep, w) for w in hits.keys())
            out_blosum[i] = best
        else:
            out_blosum[i] = blosum_score(pep, pep)

    print(f"[score] done in {time.time()-t_score:.1f}s", flush=True)

    # Log-transform counts; keep NaN for unsupported peptides (eval will fillna by median)
    out_h1_log = np.where(np.isfinite(out_h1), np.log1p(out_h1), np.nan)
    out_h2_log = np.where(np.isfinite(out_h2), np.log1p(out_h2), np.nan)

    out_df = df[["peptide", "HLA_norm"]].rename(columns={"HLA_norm": "hla"}).copy()
    out_df["self_exact_match"] = out_exact
    out_df["self_hamming1_count_log"] = out_h1_log
    out_df["self_hamming2_count_log"] = out_h2_log
    out_df["self_blosum_max"] = out_blosum
    out_df.to_csv(OUT, sep="\t", index=False)

    valid_mask = np.isfinite(out_h1)
    meta = {
        "n_rows": int(len(out_df)),
        "n_proteins": int(len(proteins)),
        "kmer_counts": {int(L): int(len(kmer[L])) for L in LENS},
        "n_scored": int(valid_mask.sum()),
        "exact_match_rate_scored": float(out_exact[valid_mask].mean()) if valid_mask.any() else 0.0,
        "mean_h1_count_scored": float(out_h1[valid_mask].mean()) if valid_mask.any() else 0.0,
        "mean_h2_count_scored": float(out_h2[valid_mask].mean()) if valid_mask.any() else 0.0,
        "elapsed_sec": time.time() - t0,
    }
    META.write_text(json.dumps(meta, indent=2))
    print(f"[done] wrote {OUT}", flush=True)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
