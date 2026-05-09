#!/usr/bin/env python3
"""Score ITSNdb scoring set with PRIME 2.1 (per-allele, parallelized)."""
import subprocess
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = Path(__file__).resolve().parent
PRIME = "/data/thca/_repo_offload/wave3_tools/PRIME/PRIME"
MIX = "/data/thca/_repo_offload/wave3_tools/MixMHCpred/MixMHCpred"

df = pd.read_csv(ROOT / "scoring_set_itsndb.tsv", sep="\t")

def prime_allele(hla):
    return hla.replace("HLA-", "").replace("*", "").replace(":", "")

df["prime_allele"] = df["HLA_norm"].apply(prime_allele)

prime_dir = ROOT / "prime_inputs"
out_dir = ROOT / "prime_outputs"
out_dir.mkdir(exist_ok=True)

allele_groups = list(df.groupby("prime_allele"))
print(f"PRIME: {len(allele_groups)} alleles to score")

def run_one(allele, grp):
    inp = prime_dir / f"peps_{allele}.txt"
    out = out_dir / f"out_{allele}.txt"
    cmd = [PRIME, "-i", str(inp), "-o", str(out), "-a", allele, "-mix", MIX]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    return allele, grp, out, r

def parse_output(allele, grp, out, r):
    if r.returncode != 0 or not Path(out).exists():
        return None, {"allele": allele, "n": len(grp),
                      "stderr": r.stderr[-500:], "stdout": r.stdout[-500:]}
    with open(out) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    if not lines:
        return None, {"allele": allele, "n": len(grp), "reason": "empty output"}
    header = lines[0].rstrip("\n").split("\t")
    rows = [ln.rstrip("\n").split("\t") for ln in lines[1:] if ln.strip()]
    out_df = pd.DataFrame(rows, columns=header)
    score_col = f"Score_{allele}"
    rank_col = f"%Rank_{allele}"
    if score_col not in out_df.columns:
        score_col = "Score_bestAllele"
        rank_col = "%Rank_bestAllele"
    sub = out_df[["Peptide", score_col, rank_col]].copy()
    sub.columns = ["peptide", "prime_score", "prime_rank"]
    sub["prime_allele"] = allele
    sub["prime_score"] = pd.to_numeric(sub["prime_score"], errors="coerce")
    sub["prime_rank"] = pd.to_numeric(sub["prime_rank"], errors="coerce")
    return sub, None

results = []
errors = []
done = 0
with ThreadPoolExecutor(max_workers=4) as ex:
    futures = {ex.submit(run_one, a, g): a for a, g in allele_groups}
    for fut in as_completed(futures):
        allele = futures[fut]
        try:
            allele, grp, out, r = fut.result()
        except Exception as e:
            errors.append({"allele": allele, "n": "?", "reason": f"exception {e!r}"})
            print(f"EXC {allele}: {e!r}")
            continue
        sub, err = parse_output(allele, grp, out, r)
        if err is not None:
            errors.append(err)
            print(f"FAILED {allele}")
            continue
        results.append(sub)
        done += 1
        if done % 5 == 0 or done == len(allele_groups):
            print(f"  done {done}/{len(allele_groups)}")

if results:
    all_res = pd.concat(results, ignore_index=True)
    print(f"PRIME scored: {len(all_res)} rows across {len(results)} alleles")
    df_m = df.merge(all_res, on=["peptide", "prime_allele"], how="left")
    df_m.to_csv(ROOT / "prime_scored.tsv", sep="\t", index=False)
    print(f"saved prime_scored.tsv; with score: {df_m['prime_score'].notna().sum()}/{len(df_m)}")
else:
    print("NO PRIME results")

if errors:
    pd.DataFrame(errors).to_csv(ROOT / "prime_errors.tsv", sep="\t", index=False)
    print(f"PRIME failures: {len(errors)} alleles")
    for e in errors:
        print(f"  {e['allele']} n={e['n']}: {e.get('reason', 'runtime')}")
