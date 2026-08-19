# RunPod GPU deployment — ESM2-650M + NetMHCpan strict-split benchmark

**Goal**: scale beyond ESM2-150M (CPU verified) to test whether further PLM scaling → continued LOSO gains.

## Cost estimate
- RunPod RTX A4000 (16 GB VRAM): $0.20-0.30/hr
- Total job: ~30 min for ESM2-650M embedding + 30 min for benchmark = ~$0.30
- RunPod RTX 3090 (24 GB) or A6000: $0.40-0.80/hr
- Best fit: A4000 spot @ $0.20/hr · ~$0.20 total

## Step 1 — Spin up RunPod instance

```bash
# RunPod CLI install (or use web UI at https://www.runpod.io/)
pip install runpod
export RUNPOD_API_KEY=...

# Pick template: PyTorch 2.x + CUDA 12.x (Ubuntu 22.04)
runpod pod create \
  --gpu-type "RTX A4000" \
  --image-name "runpod/pytorch:2.4.0-cuda12.1-py3.11" \
  --container-disk 30 \
  --volume-disk 0 \
  --name "lumenix-esm650"
```

Or via web UI:
1. https://www.runpod.io/console/deploy
2. Pick "RTX A4000" or higher
3. Image: `runpod/pytorch:2.4.0-cuda12.1-py3.11`
4. 30 GB disk (ESM2-650M weights ~3.5 GB)

## Step 2 — Upload required files (from this machine)

Files to copy to the pod:

```
/data/neoantigen_vaccine_hub/experiments/what_matters_neoantigen/cache/benchmark_clean.tsv
/data/neoantigen_vaccine_hub/data_processed/esm2_embeddings_full_meta.tsv
/data/neoantigen_vaccine_hub/scripts/runpod_esm2_650m.py  (created below)
```

Use scp / rsync to RunPod's IP:port (provided after pod start):
```bash
scp -P <pod_port> -i ~/.ssh/runpod_id \
  /data/neoantigen_vaccine_hub/experiments/what_matters_neoantigen/cache/benchmark_clean.tsv \
  /data/neoantigen_vaccine_hub/scripts/runpod_esm2_650m.py \
  root@<pod_ip>:/workspace/
```

## Step 3 — Run on the pod

```bash
ssh -p <pod_port> root@<pod_ip>

cd /workspace
pip install fair-esm scikit-learn pandas xgboost matplotlib

python3 runpod_esm2_650m.py \
  --benchmark benchmark_clean.tsv \
  --max-rows 4500 \
  --models esm2_t33_650M_UR50D
```

## Step 4 — Pull results back

```bash
scp -P <pod_port> -i ~/.ssh/runpod_id \
  root@<pod_ip>:/workspace/option_b_650m_results.tsv \
  root@<pod_ip>:/workspace/esm2_650m_embeddings.npy \
  /data/neoantigen_vaccine_hub/data_processed/

ssh -p <pod_port> root@<pod_ip> "shutdown -h now"   # stop billing
```

## Expected runtime
- ESM2-650M embedding 4,500 peptides on A4000: ~3 min
- LOSO benchmark with 5 models: ~2 min
- Total wall: ~5-10 min
- Total cost: $0.05-0.10

## Expected results

If the systematic-benchmark hypothesis holds:
- ESM2-650M LOSO AUROC ≥ ESM2-150M (currently 0.75 RF)
- NEPdb-out should approach 0.90+
- Within-source CV not expected to change much

If ESM2-650M does NOT improve LOSO:
- Confirms representation has saturated for this corpus
- Suggests task-aligned binding features (MHCflurry/NetMHCpan/PRIME) are the next-priority lever
- Strengthens the manuscript's "what matters" framing

## NetMHCpan note

NetMHCpan binary requires a license (free for academic) from DTU.
Alternative: MHCflurry (open-source, pip installable, included in this codebase).

**MHCflurry models already downloaded** at `~/.local/share/mhcflurry/4/2.2.0/models_class1_presentation/` (135 MB).

The mhcflurry script (`scripts/option_b_run.py`) had a minor API issue (allele list 6-element cap); fix:

```python
# Use Class1AffinityPredictor directly (per-row prediction):
from mhcflurry import Class1AffinityPredictor
predictor = Class1AffinityPredictor.load()
res = predictor.predict_to_dataframe(
    peptides=df["peptide"].tolist(),
    alleles=df["HLA"].tolist(),  # parallel list, not genotype
    n_flanking_length=0,
)
# Output columns: prediction (nM), prediction_percentile
```

This fix is straightforward (~15 min) and can be done locally on CPU; not blocking.

## Self-contained script for the pod

See `runpod_esm2_650m.py` (also in this directory).
