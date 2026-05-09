#!/usr/bin/env python3
"""RunPod-deploy script for AlphaFold-Multimer pMHC + TCR-pMHC complex prediction.

Designed to run on RunPod RTX A6000 / RTX 3090 / H100 GPU instance.
Uses ColabFold-batch for fast MSA-free / templates-free pipeline.

Usage on the pod (after upload):
  python3 deploy_alphafold_runpod.py --inputs candidates.fasta --gpu 0

Expected runtime per pMHC complex (HLA + β2m + peptide ~ 380 aa):
  - RTX A4000: 8-15 min
  - RTX 3090:  4-8 min
  - RTX A6000: 3-6 min
  - H100:      2-4 min

Cost estimate (RunPod spot):
  - 30 complexes on RTX 3090 spot $0.20/hr × 4hr = ~$0.80
  - 30 complexes on RTX A6000 spot $0.50/hr × 3hr = ~$1.50
"""
import argparse, json, os, time, subprocess
from pathlib import Path

# ====== Cost / time estimates ======
GPU_TIMING = {
    "RTX A4000": {"per_complex_min": 12, "$_per_hr": 0.20},
    "RTX 3090":  {"per_complex_min": 6,  "$_per_hr": 0.20},
    "RTX A6000": {"per_complex_min": 5,  "$_per_hr": 0.50},
    "H100":      {"per_complex_min": 3,  "$_per_hr": 2.00},
}


def install_colabfold():
    """One-time install on the pod."""
    print("Installing ColabFold (~30 sec on warm pod)...")
    subprocess.run(["pip", "install", "-q", "colabfold[cuda]"], check=True)
    print("ColabFold installed")


def build_pmhc_fasta(peptides, hla_alleles, out_fasta):
    """Build a multi-chain FASTA for pMHC complexes.
    Format: each complex is HLA heavy chain + β2m + peptide, separated by ':'.
    """
    # HLA-A*02:01 heavy chain (residues 1-275, alpha1+alpha2+alpha3 + TM stub)
    HLA_A0201_SEQ = (
        "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAH"
        "SQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAA"
        "QTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYP"
        "AEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE"
    )
    # β2-microglobulin (98 aa)
    B2M_SEQ = (
        "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTE"
        "FTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
    )

    with open(out_fasta, "w") as f:
        for pep, allele in zip(peptides, hla_alleles):
            name = f"{allele.replace('*','_').replace(':','_')}_{pep}"
            # Multimer format: chains separated by ':'
            seq = f"{HLA_A0201_SEQ}:{B2M_SEQ}:{pep}"
            f.write(f">{name}\n{seq}\n")
    print(f"  wrote {out_fasta} with {len(peptides)} pMHC complexes")


def run_colabfold(fasta_path, output_dir, gpu=0, num_recycle=3, num_models=1):
    """Run ColabFold-batch."""
    cmd = [
        "colabfold_batch",
        "--num-recycle", str(num_recycle),
        "--num-models", str(num_models),
        "--model-type", "alphafold2_multimer_v3",
        "--msa-mode", "single_sequence",  # No MSA = much faster but lower accuracy
        str(fasta_path),
        str(output_dir),
    ]
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    print(f"  command: {' '.join(cmd)}")
    t0 = time.time()
    subprocess.run(cmd, env=env, check=True)
    print(f"  ColabFold done in {time.time()-t0:.0f}s")


def parse_outputs(output_dir):
    """Pull pLDDT + pTM scores from ColabFold JSON outputs."""
    out_dir = Path(output_dir)
    results = []
    for json_file in out_dir.glob("*_scores_rank_*.json"):
        try:
            d = json.loads(json_file.read_text())
            results.append({
                "name": json_file.stem,
                "pLDDT_mean": d.get("plddt", [0])[0] if isinstance(d.get("plddt"), list) else d.get("plddt"),
                "pTM": d.get("ptm"),
                "ipTM": d.get("iptm"),  # inter-protein TM score
                "ranking_confidence": d.get("ranking_confidence"),
            })
        except Exception:
            pass
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-tsv", default="kras_g12_predictions.tsv",
                    help="TSV with peptide + HLA columns (top KRAS candidates)")
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--output-dir", default="/workspace/af_pmhc_results")
    ap.add_argument("--gpu", type=int, default=0)
    ap.add_argument("--skip-install", action="store_true")
    args = ap.parse_args()

    if not args.skip_install:
        install_colabfold()

    # Build top-N pMHC FASTA
    import pandas as pd
    df = pd.read_csv(args.input_tsv, sep="\t").sort_values("pred_proba", ascending=False).head(args.top_n)
    print(f"Using top-{len(df)} candidates from {args.input_tsv}")

    # Pair with HLA-A*11:01 for KRAS G12D peptides (Wells 2020 anchor)
    peptides = df["peptide"].tolist()
    hlas = ["HLA-A*11:01"] * len(peptides)
    fasta = Path(args.output_dir).parent / "candidates.fasta"
    fasta.parent.mkdir(parents=True, exist_ok=True)
    build_pmhc_fasta(peptides, hlas, fasta)

    Path(args.output_dir).mkdir(exist_ok=True, parents=True)
    run_colabfold(fasta, args.output_dir, gpu=args.gpu)

    results = parse_outputs(args.output_dir)
    print(f"\nResults summary:")
    for r in results:
        print(f"  {r['name'][:50]:50s}  pLDDT={r.get('pLDDT_mean','?')}  pTM={r.get('pTM','?')}  ipTM={r.get('ipTM','?')}")

    # Save scores
    pd.DataFrame(results).to_csv(Path(args.output_dir) / "scores_summary.tsv", sep="\t", index=False)
    print(f"\nDone. PDB files + scores in {args.output_dir}")


if __name__ == "__main__":
    main()
