# v6 Wave 2 — GPU runner bundle

Single-cell foundation-model pipeline for the THCA DIAL project. Six scripts (scGPT / Geneformer / GEARS / CellRank 2 / CellOracle / VEGA) plus a setup script and a master runner.

## Hardware

- **GPU**: NVIDIA, ≥ 16 GB VRAM (T4 / L4 / A10 / A100). Smaller GPUs (8 GB) work for everything except scGPT batch_size=32 — drop to 8 if needed.
- **RAM**: ≥ 32 GB system RAM. Some scripts (CellOracle, GEARS) hold the GRN / pretrained graph in RAM.
- **Disk**: 30 GB for downloaded models + intermediate h5ad checkpoints.

## Wall time

| GPU class | Total expected |
|---|---|
| T4 (Colab Pro+ free) | 4–6 h |
| L4 / A10 (Paperspace, AWS g5) | 2–3 h |
| A100 80 GB | 1.5–2 h |

Per-step ETAs are logged at the top of each script.

## Provisioning a GPU box

**Colab Pro+** (cheapest path, ≈ $50/mo):
1. Pro+ subscription, request L4 / A100 from runtime menu.
2. Mount Drive. Upload this `v6_gpu_scripts.zip`.
3. `!unzip -q v6_gpu_scripts.zip -d /content/v6_gpu && cd /content/v6_gpu && bash 00_setup.sh`
4. Stage the input h5ad `classical_baseline.h5ad` to `/data/thca/scrna/processed/` (or pass `INPUT=/path/to/your.h5ad`).
5. `bash run_all.sh`
6. After completion, rsync results back to your project: `rsync -av /opt/thyroid-dash/project/results/v6_scrna/ user@host:/opt/thyroid-dash/project/results/v6_scrna/`

**Paperspace Gradient** (≈ $0.5–1.5/h, no time cap):
1. Spin up an L4 or A10 instance with the PyTorch image.
2. `scp v6_gpu_scripts.zip paperspace:/notebooks/`
3. Same as Colab steps 3–6.

**Lambda Cloud / RunPod / Vast.ai**: any CUDA 12 image with PyTorch ≥ 2.2 works.

## Step order and inputs

| Step | Script | Inputs | Outputs |
|---|---|---|---|
| 0 | `00_setup.sh` | none | venv, pinned packages, `logs/v6_gpu_setup.log` |
| 1 | `10_scgpt_embedding.py` | `--input` h5ad | `scgpt_cell_embeddings.npy`, attention table, DIAL table, cell-type table |
| 2 | `20_geneformer.py` | h5ad | `geneformer_ranking.tsv` (gene × attention) |
| 3 | `30_gears_perturb.py` | h5ad | single + combo KO, BRAF-like-B rescue tables |
| 4 | `40_cellrank.py` | h5ad with spliced/unspliced layers | macrostate h5ad, driver gene table, interactive UMAP |
| 5 | `50_celloracle_tf.py` | h5ad | TF screening table, top-10 candidates md |
| 6 | `60_vega.py` | h5ad + Hallmark gmt | pathway activation table, per-cell pathway DIAL |

## CellRank velocity caveat

`40_cellrank.py` requires `spliced` and `unspliced` layers in the AnnData. If your input lacks them, re-process from raw fastq with velocyto or skylight, or use a pre-computed Pu et al. 2021 velocyto checkpoint (search for `GSE184362_velocyto.h5ad`). Without velocity layers, this step exits 2 with a clear error.

## Continuation on error

`run_all.sh` defaults to halting on first failure. Set `CONTINUE_ON_ERROR=1 bash run_all.sh` to push through and collect every step's exit code.

## Stamp files

Each script writes `$RES/.stamp_<scriptname>` on success. Quick progress view:
```
ls -la /opt/thyroid-dash/project/results/v6_scrna/.stamp_*
```

## Re-rendering the v6 dashboard after a Wave 2 run

The dashboard at `/reports/html/pages/v6_scrna_foundation.html` reads the multi-scale DIAL TSV and the figure HTMLs. After Wave 2 finishes:
1. Rebuild the multi-scale DIAL bar chart by appending scGPT and VEGA rows to `results/v6_scrna/dial_scrna/multi_resolution_dial.tsv`.
2. Re-run the small chart-building snippet inlined in the dashboard build script (or just re-run the chart cell in `notebooks_or_scripts/v6_build_dashboard.py` if a sibling rebuild script exists).
3. The dashboard refreshes on next load (no rebuild needed for iframes).

## Troubleshooting

- `scgpt` import fails → install with `pip install git+https://github.com/bowang-lab/scGPT.git` (the PyPI package is stale).
- `vega` import fails → upstream package was renamed; try both `vega-tools` and the LucasESBS git repo.
- CUDA OOM in scGPT → drop `--batch_size` to 8 or 16.
- Geneformer slow → cap cells with the in-script `n_cells = min(adata_b.n_obs, 5000)` line.
- CellOracle base GRN download blocked → pre-stage the file from `celloracle.data.load_human_promoter_base_GRN()` on a machine with internet, then copy `~/.celloracle_data/` to the GPU box.

## Citing

Each step's pre-trained model and method has its own paper — see `reports/v6/v6_scrna_foundation_section.md` for the canonical citation list.
