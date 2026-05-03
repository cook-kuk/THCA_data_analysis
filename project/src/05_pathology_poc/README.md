# Pathology POC — H&E tile → DM1_like_score regression

Goal: prove that an off-the-shelf foundation-model embedding of H&E tiles can predict
the ST-derived DM1_like_score (= -RAI_8 z-score, within-sample). If yes, Paper 1 reach
gains a "tissue → molecular" inference angle without external IHC.

## What `extract_tiles.py` produces

- `data/processed/GSE250521/tiles/{sample_id}/{spot_id}_size{tile_size}.png`
  spot-centered tile from the hires H&E (Visium tissue_hires_image), default 224×224
  in hires-image pixel space (= ~110 µm at default scalefactor; tune via `--tile-size`)
- `results/03_pathology_poc/tile_metadata.tsv.gz`
  one row per tile: spot_id, sample_id, stage, tile_path, tile_size, x_hires, y_hires,
  RAI_8_score, DM1_like_score, TDS_like_score

## TODO (skeleton, not yet implemented)

```python
# 1. Embed every tile with a frozen foundation model
#    Candidates (require manual access acceptance):
#    - UNI       MahmoodLab/UNI                 (HF, gated, non-commercial)
#    - CONCH     MahmoodLab/CONCH               (HF, gated, non-commercial)
#    - Virchow2  paige-ai/Virchow2              (HF, gated)
#    - Prov-GigaPath prov-gigapath/prov-gigapath (HF, gated)
#    Fast fallback: torchvision ResNet50 ImageNet weights — weak but unblocked.

import torch, torchvision.models as M
import torchvision.transforms as T
from PIL import Image

# choose backbone
backbone = M.resnet50(weights=M.ResNet50_Weights.IMAGENET1K_V2)
backbone.fc = torch.nn.Identity()
backbone.eval().cuda()

xform = T.Compose([T.Resize(224), T.CenterCrop(224), T.ToTensor(),
                   T.Normalize([0.485,0.456,0.406],[0.229,0.224,0.225])])

# 2. Forward pass, save embeddings to .npz keyed by spot_id
# 3. Train Ridge with leave-one-slide-out CV
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import roc_auc_score

# meta = pd.read_csv(...tile_metadata.tsv.gz, sep='\t')
# embeddings = np.load(...embeddings.npz)
# for held_out_sample in meta['sample_id'].unique():
#     train = meta[meta.sample_id != held_out_sample]
#     test  = meta[meta.sample_id == held_out_sample]
#     X_tr = embeddings[train.spot_id.values]; y_tr = train.DM1_like_score.values
#     X_te = embeddings[test.spot_id.values];  y_te = test.DM1_like_score.values
#     mdl = Ridge(alpha=1.0).fit(X_tr, y_tr); pred = mdl.predict(X_te)
#     # metrics: spearman(y_te,pred), pearson, AUROC top25 vs bottom25
```

## Why leave-one-slide-out, not random split

Spots within a slide share batch effects (sectioning, staining, sequencing depth).
Random tile split silently inflates performance because the model memorizes per-slide
artifacts. LOSO is the only honest validation for a 16-slide POC.

## Expected POC metrics (rough, tune after run)

- Spearman r ≥ 0.30 across LOSO folds → POC works, justifies foundation-model upgrade
- AUROC (DM1 top 25% vs bottom 25%) ≥ 0.70 → tile-level discrimination is real
- below those → ResNet50 too weak; gate UNI/CONCH access and retry

## Scope guard (Paper isolation)

- This POC is Paper 1 (cancer / 8-gene driver-excluded sub-stratifier) territory.
- Do NOT extend to HT/TLS/BCR clonality (Paper 2) or GD HLA (Paper 3) here.
- DM1_like_score is the only target. Do NOT add Hashimoto signature regression in this script.
