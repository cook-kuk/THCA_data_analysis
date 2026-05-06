#!/bin/bash
# G2 — TCGA-THCA H&E digital pathology pipeline.
# 1. Download diagnostic FFPE SVS slides from GDC (~500, 100-200 GB)
# 2. Tile each slide at 256×256 patches with overlap, filter background
# 3. Embed with same foundation model used in G1
# 4. Aggregate to slide-level → predict bulk DM1 score (S TCGA result)
# 5. Cross-validate
set -euo pipefail

DATA_DIR="${1:-data/g2_tcga_he}"
SCRIPT_DIR="$(dirname "$0")"
mkdir -p "$DATA_DIR"

echo "[G2] Step 1 — manifest from GDC API for TCGA-THCA diagnostic slides"
python3 - <<EOF
import json, urllib.request, sys, pathlib
filt = {"op":"and","content":[
    {"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-THCA"]}},
    {"op":"in","content":{"field":"data_type","value":["Slide Image"]}},
    {"op":"in","content":{"field":"experimental_strategy","value":["Diagnostic Slide"]}},
]}
url = "https://api.gdc.cancer.gov/files"
import urllib.parse
params = urllib.parse.urlencode({
    "filters": json.dumps(filt),
    "size": "1000",
    "format": "JSON",
    "fields": "file_name,file_id,file_size,cases.submitter_id"
})
with urllib.request.urlopen(f"{url}?{params}", timeout=60) as r:
    data = json.load(r)
hits = data["data"]["hits"]
print(f"GDC returned {len(hits)} TCGA-THCA diagnostic slides")
out_path = pathlib.Path("$DATA_DIR/manifest.tsv")
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w") as f:
    f.write("file_id\tfile_name\tfile_size\tcase_id\n")
    for h in hits:
        case = h.get("cases",[{}])[0].get("submitter_id","")
        f.write(f"{h['file_id']}\t{h['file_name']}\t{h['file_size']}\t{case}\n")
print(f"manifest → {out_path}")
EOF

echo "[G2] Step 2 — install gdc-client + openslide-python on RunPod"
echo "  (run manually if missing):"
echo "    pip install openslide-python tiatoolbox"
echo "    sudo apt-get install -y openslide-tools"
echo "    wget https://gdc.cancer.gov/files/public/file/gdc-client_2.3.0_Ubuntu_x64-py3.8-ubuntu-20.04.zip"

echo "[G2] Step 3 — download diagnostic slides (use gdc-client manifest)"
echo "  gdc-client download -m $DATA_DIR/manifest_gdc.txt -d $DATA_DIR/slides/"
echo "  Estimated: 100-200 GB."

echo "[G2] Step 4 — tile + embed: see g2_tile_embed.py"
echo "[G2] Step 5 — train slide-level regression: see g2_slide_regress.py"

echo
echo "Note: G2 download requires significant bandwidth + disk. RunPod 200 GB+ disk recommended."
echo "If disk-constrained, use g2_tile_embed.py with --stream-mode (download → tile → embed → delete in pipeline)"
