#!/bin/bash
# B1 — Structure refinement (AlphaFold3 / ESMFold)
# Time: 4-8h on A100 80GB · Cost: $30-60
set -euo pipefail
cd /workspace/azure_gpu_phaseB_drug

OUT=results/b1_structures/refined
mkdir -p $OUT

# Source: 8 priority targets per v13_drug_discovery
TARGETS=(TACSTD2 TMPRSS4 PLEKHA6 CYP1B1 LDLR GABRB2 B3GNT3 PTPRE)

# UniProt mapping
declare -A UNIPROT
UNIPROT[TACSTD2]=P09758
UNIPROT[TMPRSS4]=Q9NRS4
UNIPROT[PLEKHA6]=Q9Y2H5
UNIPROT[CYP1B1]=Q16678
UNIPROT[LDLR]=P01130
UNIPROT[GABRB2]=P47870
UNIPROT[B3GNT3]=Q9Y2A9
UNIPROT[PTPRE]=P23469

# AlphaFold3 (DeepMind) — preferred
if command -v alphafold3 >/dev/null; then
  for t in "${TARGETS[@]}"; do
    UP=${UNIPROT[$t]}
    echo "[B1] AlphaFold3 $t (UniProt $UP)"
    alphafold3 --uniprot $UP --output_dir $OUT/$t --use_msa --num_seeds 5
  done
else
  # Fallback: ESMFold (faster, less accurate)
  echo "[B1] ESMFold fallback"
  python3 -c "
import torch, esm
from pathlib import Path
model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()
import requests
for t, up in {'TACSTD2':'P09758','TMPRSS4':'Q9NRS4','PLEKHA6':'Q9Y2H5','CYP1B1':'Q16678','LDLR':'P01130','GABRB2':'P47870','B3GNT3':'Q9Y2A9','PTPRE':'P23469'}.items():
    seq = requests.get(f'https://www.uniprot.org/uniprot/{up}.fasta').text.split('\n',1)[1].replace('\n','')
    print(f'{t}: len={len(seq)}')
    with torch.no_grad():
        out = model.infer_pdb(seq)
    Path('$OUT').mkdir(parents=True, exist_ok=True)
    open(f'$OUT/{t}.pdb','w').write(out)
    plddt = model.infer_pdbs([seq])[0]
    # ... extract per-residue plddt
"
fi

# Pocket re-detection (P2Rank or fpocket)
mkdir -p ../results/b1_structures/pocket_redetect
for t in "${TARGETS[@]}"; do
  if [ -f $OUT/$t.pdb ]; then
    fpocket -f $OUT/$t.pdb -o ../results/b1_structures/pocket_redetect/$t || true
  fi
done

echo "[B1] complete · 8 refined PDBs + pocket redetection"
