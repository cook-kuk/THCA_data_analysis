#!/bin/bash
# Pod C — AlphaFold HLA-DPB1*05:01 + TSHR/Tg/TPO autoreactive peptide structure
# Run inside RunPod A100 80GB pod after SSH connect
#
# Strategy: ColabFold (LocalColabFold) + 5-10 HLA-peptide complex predictions
# Expected runtime: 8-12 hr
# Expected cost: $12-15
set -e

WORK=/workspace/alphafold
mkdir -p $WORK/inputs $WORK/outputs $WORK/databases
cd $WORK

# ============================================================
# Phase 1 — Install LocalColabFold
# ============================================================
echo "=== [Phase 1] LocalColabFold install ==="
DEBIAN_FRONTEND=noninteractive apt-get update -qq
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq curl wget git bzip2

# Install LocalColabFold (~5 min)
wget -q https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
bash install_colabbatch_linux.sh 2>&1 | tail -5

export PATH="$WORK/localcolabfold/colabfold-conda/bin:$PATH"
which colabfold_batch
colabfold_batch --help 2>&1 | head -5

# ============================================================
# Phase 2 — Prepare HLA-peptide input sequences
# ============================================================
echo "=== [Phase 2] HLA-peptide complex inputs ==="
cd $WORK/inputs

# HLA-DPB1*05:01 alpha + beta chains + TSHR/Tg/TPO autoreactive 9-mer peptides
# (Using known T-cell epitopes from IEDB for thyroid autoimmunity)

cat > hla_dpb1_0501_tshr_p1.fasta <<EOF
>HLA-DPA1*02:02
IKADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGRAFSFEAQGGLANIAILNNNLNTLIQRSNHTQATNDPPEVTVFPKEPVELGQPNTLICHIDKFFPPVLNVTWLCNGELVTEGVAESLFLPRTDYSFHKFHYLPFLPSAEDFYDCKVEHWGLDQPLLKHWE
>HLA-DPB1*05:01
RATPENYLFQGRQECYAFNGTQRFLERYIYNREEFARFDSDVGEFRAVTELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVGESFTVQRRVEPKVTVYPSKTQPLQHHNLLVCSVSGFYPGSIEVRWFRNGQEEKAGVVSTGLIQNGDWTFQTLVMLETVPRSGEVYTCQVEHPSV
>TSHR_p1_residues_252-261
KLEFLNVCFC
EOF

cat > hla_dpb1_0501_tg_p1.fasta <<EOF
>HLA-DPA1*02:02
IKADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGRAFSFEAQGGLANIAILNNNLNTLIQRSNHTQATNDPPEVTVFPKEPVELGQPNTLICHIDKFFPPVLNVTWLCNGELVTEGVAESLFLPRTDYSFHKFHYLPFLPSAEDFYDCKVEHWGLDQPLLKHWE
>HLA-DPB1*05:01
RATPENYLFQGRQECYAFNGTQRFLERYIYNREEFARFDSDVGEFRAVTELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVGESFTVQRRVEPKVTVYPSKTQPLQHHNLLVCSVSGFYPGSIEVRWFRNGQEEKAGVVSTGLIQNGDWTFQTLVMLETVPRSGEVYTCQVEHPSV
>Thyroglobulin_p1_residues_2549-2570
DAKFKAFFNAYVDDLPKKMS
EOF

cat > hla_dpb1_0501_tpo_p1.fasta <<EOF
>HLA-DPA1*02:02
IKADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGRAFSFEAQGGLANIAILNNNLNTLIQRSNHTQATNDPPEVTVFPKEPVELGQPNTLICHIDKFFPPVLNVTWLCNGELVTEGVAESLFLPRTDYSFHKFHYLPFLPSAEDFYDCKVEHWGLDQPLLKHWE
>HLA-DPB1*05:01
RATPENYLFQGRQECYAFNGTQRFLERYIYNREEFARFDSDVGEFRAVTELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVGESFTVQRRVEPKVTVYPSKTQPLQHHNLLVCSVSGFYPGSIEVRWFRNGQEEKAGVVSTGLIQNGDWTFQTLVMLETVPRSGEVYTCQVEHPSV
>TPO_p1_residues_535-552
GMPFITPSPNTLPEEAKE
EOF

# B*46:01 + TSHR (Asian-specific HLA-I context)
cat > hla_b4601_tshr_p1.fasta <<EOF
>HLA-B*46:01
GSHSMRYFYTAMSRPGRGEPRFITVGYVDDTLFVRFDSDAASPRTEPRAPWIEQEGPEYWDRETQISKTNTQTYRESLRNLRGYYNQSEDGSHTLQRMYGCDLGPDGRLLRGYDQYAYDGKDYIALNEDLSSWTAADTAAQITQRKWEAARVAEQLRAYLEGTCVEWLRRYLENGKETLQRA
>TSHR_p1_residues_55-63
ALTDLELHV
EOF

# ============================================================
# Phase 3 — Run AlphaFold (4 structures)
# ============================================================
echo "=== [Phase 3] Run AlphaFold predictions ==="
cd $WORK
for f in inputs/hla_dpb1_0501_tshr_p1.fasta \
         inputs/hla_dpb1_0501_tg_p1.fasta \
         inputs/hla_dpb1_0501_tpo_p1.fasta \
         inputs/hla_b4601_tshr_p1.fasta; do
    name=$(basename $f .fasta)
    echo "--- $name ---"
    colabfold_batch $f $WORK/outputs/$name \
                     --num-models 3 \
                     --num-recycle 3 \
                     --rank plddt 2>&1 | tail -10
done

ls -la $WORK/outputs/

# ============================================================
# Phase 4 — Extract pLDDT + interface pTM scores
# ============================================================
echo "=== [Phase 4] Score extraction ==="
python3 << 'EOF'
import os, json, glob
results = []
for outdir in glob.glob("/workspace/alphafold/outputs/*"):
    name = os.path.basename(outdir)
    for json_f in glob.glob(f"{outdir}/*_scores_rank_001*.json"):
        with open(json_f) as f:
            data = json.load(f)
        results.append({
            "structure": name,
            "plddt_mean": sum(data.get("plddt", [0])) / len(data.get("plddt", [1])),
            "ptm": data.get("ptm", None),
            "iptm": data.get("iptm", None),
        })
import json
print(json.dumps(results, indent=2))
with open("/workspace/alphafold/outputs/SUMMARY.json", "w") as f:
    json.dump(results, f, indent=2)
EOF

ls -la $WORK/outputs/SUMMARY.json
echo "=== ALL DONE ==="
