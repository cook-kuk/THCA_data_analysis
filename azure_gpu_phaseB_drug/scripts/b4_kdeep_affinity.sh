#!/bin/bash
# B4 — KDeep / DeepPurpose binding affinity prediction
# Time: 4-6h on A100 40GB · Cost: $25-40
set -euo pipefail
cd /workspace/azure_gpu_phaseB_drug

OUT=results/b4_kdeep
mkdir -p $OUT

python3 << 'PYEOF'
import pandas as pd, torch
from DeepPurpose import models, utils
from pathlib import Path
from collections import OrderedDict

# Load 8 target sequences
import requests
TARGETS = OrderedDict([
    ('TACSTD2', 'P09758'), ('TMPRSS4', 'Q9NRS4'), ('PLEKHA6', 'Q9Y2H5'),
    ('CYP1B1', 'Q16678'), ('LDLR', 'P01130'), ('GABRB2', 'P47870'),
    ('B3GNT3', 'Q9Y2A9'), ('PTPRE', 'P23469')
])
seqs = {}
for t, up in TARGETS.items():
    fa = requests.get(f'https://www.uniprot.org/uniprot/{up}.fasta').text
    seqs[t] = fa.split('\n', 1)[1].replace('\n', '')

# Load compounds
df = pd.read_csv('/workspace/v13_drug_discovery/prioritization/ranked_candidates.tsv', sep='\t')

# DeepPurpose model — use pre-trained MPNN_CNN (DAVIS)
config = utils.generate_config(drug_encoding='MPNN', target_encoding='CNN', cls_hidden_dims=[1024,1024,512], train_epoch=50, LR=0.001, batch_size=128)
model = models.model_initialize(**config)
model.load_pretrained('/workspace/DeepPurpose/save_folder/MPNN_CNN_DAVIS')

# Predict
out_rows = []
for t, seq in seqs.items():
    sub = df[df.get('target', '') == t].head(50)
    if sub.empty:
        sub = df.head(50)  # fallback
    smiles = sub['smiles'].dropna().tolist()
    if not smiles:
        continue
    Y = model.predict(utils.data_process(smiles, [seq] * len(smiles), [0]*len(smiles), drug_encoding='MPNN', target_encoding='CNN', split_method='no_split'))
    for s, y in zip(smiles, Y):
        out_rows.append({'target': t, 'smiles': s, 'pkd_pred': y})

pd.DataFrame(out_rows).to_csv('results/b4_kdeep/pkd_predictions.tsv', sep='\t', index=False)
print(f'[B4] {len(out_rows)} predictions')
PYEOF

echo "[B4] complete · pKd predictions"
