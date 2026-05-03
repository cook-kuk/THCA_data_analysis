# Azure T4 — Phase A GPU pipeline (H&E → DM1_like_score_resid)

**Goal**: Verify whether ST-aligned H&E tiles can predict the depth-residualized DM1_like score (Paper 1's 8-gene RAI/DM1 axis) using ResNet50 ImageNet embeddings + Ridge LOSO. Result drives go/no-go for Paper 2 viability and IP feasibility.

**Hard rules** (encoded in `run_phaseA_gpu.sh`)
- random tile split **BANNED** (leave-one-slide-out only)
- target = `DM1_like_score_resid` (depth-residualized) only — no raw score
- no RET fusion claim
- no TCGA WSI download in this phase
- T4 only; **no A100** until T4 Phase A passes

---

## 0. Files in this package

```
phaseA_gpu_pkg/
├── README_AZURE_GPU_PHASEA.md            ← this file
├── requirements_pathology_gpu.txt
├── run_phaseA_gpu.sh                     ← entrypoint
└── thca_phaseA/                          ← THCA_ROOT after untar
    └── project/
        ├── src/05_pathology_poc/         ← 6 .py scripts
        ├── data/processed/GSE250521/
        │   ├── tiles/  (3,200 × 224 px PNG, ~320 MB)
        │   └── GSM*/*.scored.h5ad  (16 files, ~329 MB, for neg ctrl)
        ├── results/01_spatial_score/all_spots_scored.tsv.gz
        ├── results/03_pathology_poc/tile_metadata*.tsv.gz
        └── reports/                      ← verdict written here
```

Tarball size ≈ **~660 MB**. SCP over Azure peering ~5–10 min from a typical home connection.

---

## 1. Provision the Azure T4 VM

```bash
# Local (Azure CLI). Replace RG, NAME, REGION as appropriate.
RG="thca-pathA-rg"
NAME="thca-pathA-t4"
REGION="koreacentral"           # or eastus2 if T4 quota easier there

az group create -n "$RG" -l "$REGION"

az vm create \
  --resource-group "$RG" \
  --name "$NAME" \
  --image Canonical:0001-com-ubuntu-server-jammy:22_04-lts-gen2:latest \
  --size Standard_NC4as_T4_v3 \
  --admin-username azureuser \
  --ssh-key-values ~/.ssh/id_rsa.pub \
  --os-disk-size-gb 128 \
  --public-ip-sku Standard \
  --priority Spot --eviction-policy Deallocate --max-price -1
# (Spot saves ~70%. Drop the spot flags if you need a guaranteed run.)

# Open SSH only (default). Get the public IP:
az vm show -d -g "$RG" -n "$NAME" --query publicIps -o tsv
```

If T4 quota is denied in your subscription, request a quota increase or pick another region. **Do not auto-fall back to A100** — the spec disallows it.

---

## 2. NVIDIA driver + CUDA torch (one-time on the VM)

SSH in:

```bash
ssh azureuser@<GPU_VM_IP>
```

Then on the VM:

```bash
# NVIDIA driver (Ubuntu 22.04 default repo includes nvidia-driver-535 which works for T4)
sudo apt-get update
sudo apt-get install -y build-essential python3-venv python3-pip nvidia-driver-535
sudo reboot
# wait ~60 s, ssh back in
nvidia-smi    # should show "Tesla T4" with driver 535+

# Python env
python3 -m venv ~/venv-pathA
source ~/venv-pathA/bin/activate
pip install --upgrade pip
pip install torch==2.4.0 torchvision==0.19.0 --index-url https://download.pytorch.org/whl/cu121
# requirements ship inside the package but you can pre-install:
# pip install numpy pandas scipy scikit-learn statsmodels matplotlib Pillow anndata h5py tqdm
```

Verify GPU is visible to PyTorch:

```bash
python3 -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
# expect: True Tesla T4
```

---

## 3. Build + transfer the package (run on the local VM)

```bash
# Local (this VM, not the GPU box):
bash /home/seungho/personal/THCA_data_analysis/azure_gpu_phaseA/package_for_gpu.sh
# → /home/seungho/personal/THCA_data_analysis/azure_gpu_phaseA/phaseA_gpu_pkg.tar.gz (~660 MB)

# Send to GPU VM:
scp /home/seungho/personal/THCA_data_analysis/azure_gpu_phaseA/phaseA_gpu_pkg.tar.gz \
    azureuser@<GPU_VM_IP>:~/

# Or rsync if resumable transfer is preferable:
rsync -avhP /home/seungho/personal/THCA_data_analysis/azure_gpu_phaseA/phaseA_gpu_pkg.tar.gz \
    azureuser@<GPU_VM_IP>:~/
```

---

## 4. Run Phase A on the GPU VM

```bash
ssh azureuser@<GPU_VM_IP>
source ~/venv-pathA/bin/activate
tar xzf ~/phaseA_gpu_pkg.tar.gz
cd ~/phaseA_gpu_pkg
pip install -r requirements_pathology_gpu.txt
bash run_phaseA_gpu.sh 2>&1 | tee run_phaseA_$(date +%Y%m%d_%H%M%S).log
```

Expected wall time on T4 ≈ **5–10 min total** (3,200 tile embed ≈ 30 s on T4 + LOSO + 100 random panels).

Steps the script runs (each prints its section header):
0. `nvidia-smi` + torch CUDA check (aborts if no GPU)
1. `compute_resid_labels.py` — re-residualize labels (idempotent)
2. `embed_resnet50.py --tile-size 224 --device cuda` → `embeddings_resnet50_224.npz`
3. `train_loso_ridge.py` — Ridge + ElasticNet, 3 targets, 16-fold LOSO
4. `negative_controls.py --n-random 100` — housekeeping + 100 random 8-gene panels
5. `plot_pred_vs_obs.py` — pooled scatter + per-slide bar
6. **Auto-verdict**: writes `pathology_dm1_phaseA_gpu_verdict_<date>.md` with GO/BORDERLINE/NO-GO

---

## 5. Send results back to local

```bash
# On GPU VM, pack results:
tar czf ~/phaseA_results.tar.gz \
    -C ~/phaseA_gpu_pkg/thca_phaseA/project \
    results/03_pathology_poc \
    reports

# Local:
scp azureuser@<GPU_VM_IP>:~/phaseA_results.tar.gz \
    /home/seungho/personal/THCA_data_analysis/azure_gpu_phaseA/

# Untar in place to merge back into the project tree (only path that matters
# is project/reports/* and project/results/03_pathology_poc/*):
cd /home/seungho/personal/THCA_data_analysis
tar xzf azure_gpu_phaseA/phaseA_results.tar.gz -C project/
```

---

## 6. Auto-gate decisions

After `run_phaseA_gpu.sh` finishes, read the verdict file:

| Verdict | Meaning | Next |
|---|---|---|
| **GO** | Spearman ≥ 0.30 AND AUROC ≥ 0.70 AND controls weak | Build & ship `phaseB_gpu_pkg.tar.gz` (multi-cohort LODO). Same VM. |
| **BORDERLINE** | Spearman ∈ [0.20, 0.30) OR AUROC ∈ [0.65, 0.70) | Apply for UNI/CONCH/Virchow2 HF access. Re-run with foundation embedder. |
| **NO-GO** | Spearman < 0.20 AND AUROC < 0.65 | Stop. Reframe concept before more investment. |

---

## 7. Cost / cleanup

`Standard_NC4as_T4_v3` Spot ≈ $0.12–0.18 / hr (region-dependent). 10-min Phase A ≈ **<$0.10**.

After verdict + results pulled back, **deallocate** the VM (don't delete RG yet — keep for Phase B re-use):

```bash
az vm deallocate -g thca-pathA-rg -n thca-pathA-t4
# (or full delete after Phase B is done:)
# az group delete -n thca-pathA-rg --yes
```

---

## 8. Failure modes & mitigations

| Symptom | Likely cause | Fix |
|---|---|---|
| `nvidia-smi: command not found` | driver not installed or no reboot after install | `sudo apt install nvidia-driver-535 && sudo reboot` |
| `torch.cuda.is_available() == False` | wrong torch wheel (cpu-only) | reinstall with `--index-url https://download.pytorch.org/whl/cu121` |
| OOM during embed | should not happen on T4 with batch=128 + 224 px ResNet50 | drop `--batch` to 64 |
| `permission denied: ssh-key` for VM create | wrong key path | `--ssh-key-values $(cat ~/.ssh/id_rsa.pub)` literal |
| spot eviction mid-run | spot price > max | re-create with `--priority Regular` (~3× cost) |
| `extract_tiles.py` referenced in errors | not needed in Phase A 224-only run; ignore | n/a |

---

## 9. What this run does NOT do (deliberately)

- ❌ TCGA WSI anything
- ❌ A100 / H100
- ❌ UNI / CONCH / Virchow2 (gated; only after BORDERLINE/GO)
- ❌ 448 / 672 px multi-resolution (deferred to Phase B if A passes)
- ❌ Phase B / C auto-trigger (verdict prints next command, you ship next package)

---

## 10. Generated locally

- `package_for_gpu.sh` — tarball builder
- `run_phaseA_gpu.sh` — bundled GPU runner
- `requirements_pathology_gpu.txt` — pip deps
- 6 Python scripts under `project/src/05_pathology_poc/`

All scripts honor `THCA_ROOT` env var (defaults to local path); GPU runner sets it to `$HOME/phaseA_gpu_pkg/thca_phaseA`.
