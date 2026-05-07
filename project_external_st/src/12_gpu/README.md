# Phase B GPU package — Cancer Cell push (G1 + G2 + G3)

## Quick start (RunPod)

```bash
# 1. unpack
tar xzf phaseB_gpu_pkg.tar.gz
cd phaseB_gpu_pkg

# 2. install deps
bash setup.sh

# 3. G1 — pathology AI POC (5,600 H&E tiles → DM1 prediction LOSO)
python g1_embed.py --model resnet50 --meta data/all_tile_metadata.tsv.gz --out results/embeddings_resnet50.npz
python g1_ridge_loso.py --embed results/embeddings_resnet50.npz --meta data/all_tile_metadata.tsv.gz --target DM1_like_score_resid --out-dir results/g1_resnet50

# (optional) UNI foundation model — requires HuggingFace gated access
export HF_TOKEN=<your_token>
python g1_embed.py --model uni --out results/embeddings_uni.npz
python g1_ridge_loso.py --embed results/embeddings_uni.npz --target DM1_like_score_resid --out-dir results/g1_uni

# 4. G2 — TCGA-THCA digital pathology (n=500 H&E + bulk RNA-seq)
bash g2_tcga_he_pipeline.sh data/g2_tcga_he      # builds GDC manifest
gdc-client download -m data/g2_tcga_he/manifest_gdc.txt -d data/g2_tcga_he/slides/   # 100-200 GB
python g2_tile_embed.py --manifest data/g2_tcga_he/manifest.tsv --slides-dir data/g2_tcga_he/slides --model resnet50 --delete-after
python g2_slide_regress.py

# 5. G3 — cell-type deconvolution (paired scRNA → ST)
python g3_celldart.py --scrna-dir data/scrna_gse250521 --st-meta data/all_tile_metadata.tsv.gz
```

## Hardware requirements

| Sprint | GPU | Disk | Time |
|---|---|---|---|
| G1 ResNet50 | any (T4 OK) | 1 GB | 10-30 min |
| G1 UNI | T4+ (8 GB VRAM) | 2 GB | 30-60 min |
| G2 download | — | **100-200 GB** | 4-12 h bandwidth |
| G2 tile + embed | T4+ | rolls disk if `--delete-after` | 4-8 h GPU |
| G2 regression | any | 1 GB | 5 min |
| G3 CellDART | T4+ | 1 GB | 30-60 min |

**Recommended RunPod**: A100 / 4090, 200 GB disk, 32 GB RAM.

## Foundation model access

| Model | License | Setup |
|---|---|---|
| ResNet50 (default) | none | already in torchvision |
| UNI | gated, MahmoodLab HF | accept license + HF_TOKEN env var |
| CONCH | gated, MahmoodLab HF | accept license + HF_TOKEN env var |
| Virchow2 | gated, paige-ai HF | accept license + HF_TOKEN env var |

If gating is a barrier, **start with ResNet50** — POC published with ImageNet weights, then re-run with UNI when access cleared.

## Expected output (success criteria)

| Output | Path | Cancer Cell threshold |
|---|---|---|
| G1 LOSO Spearman | `results/g1_resnet50/summary.tsv` | **r > 0.40** |
| G1 LOSO AUROC | same | **AUC > 0.70** |
| G2 5-fold CV Spearman | `results/g2_regression/fold_metrics.tsv` | **r > 0.50** |
| G3 cell-type panel | `results/g3_celldart/celltype_panel_expression.tsv` | DM1-high cells = thyrocyte (epithelial) class |

## scp commands (from RunPod back to source VM)

```bash
# bring back results
scp -r results/ seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project_external_st/results/extra/phaseB_gpu/
```
