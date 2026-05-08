#!/bin/bash
# Pod-side runner. Assumes /workspace/neo_bayes contains bundle.tsv, hla_pseudo.tsv,
# embed_esm2.py, train_bayesian_neo.py.
set -e
cd /workspace/neo_bayes

# Faster Hugging Face download path
export HF_HUB_ENABLE_HF_TRANSFER=0
export TRANSFORMERS_OFFLINE=0

echo "=== Step 1: install missing deps ==="
pip install -q transformers==4.40.0 accelerate scikit-learn pandas numpy matplotlib seaborn safetensors huggingface_hub 2>&1 | tail -5

echo "=== Step 2: embeddings ==="
python3 embed_esm2.py 2>&1 | tail -40

echo "=== Step 3: training + eval (5 seeds, 20 epochs each) ==="
python3 train_bayesian_neo.py --bundle bundle.tsv --embeddings embeddings.pt \
    --out_dir . --epochs 20 --n_seeds 5 --T 30 --mixup 0.2 2>&1 | tail -200

echo "=== Step 4: ablation no-DANN ==="
mkdir -p ablation_no_dann
python3 train_bayesian_neo.py --bundle bundle.tsv --embeddings embeddings.pt \
    --out_dir ablation_no_dann --epochs 20 --n_seeds 5 --T 30 --mixup 0.2 --no_dann 2>&1 | tail -50

echo "=== DONE ==="
ls -la
