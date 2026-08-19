"""Crop each main figure composite into 7 individual panel PNG thumbnails.

Reads composite PNGs and slices them into A-G sub-images per the 2x4 grid layout used in main figure scripts:
   row1 cols (0,1,2,3) = panels A B C D
   row2 cols (0,1,2,3) = panels E F G + (often G spans cols 2-3)

For each fig we use known panel layouts to determine bounding boxes.
"""
import os, numpy as np
from PIL import Image

ROOT_OUT = "/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main"
ASSETS   = "/home/seungho/personal/THCA_data_analysis/project/papers_hub_2026_05_04/assets/paper1_nc"
os.makedirs(ASSETS, exist_ok=True)

# Panel grid specifications: figure_no -> list of (label, (col_start, col_end, row_start, row_end)) on a 4-col 2-row grid
panel_specs = {
    1: [("A",(0,1,0,1)), ("B",(1,2,0,1)), ("C",(2,3,0,1)), ("D",(3,4,0,1)),
        ("E",(0,1,1,2)), ("F",(1,2,1,2)), ("G",(2,4,1,2))],
    2: [("A",(0,1,0,1)), ("B",(1,2,0,1)), ("C",(2,3,0,1)), ("D",(3,4,0,1)),
        ("E",(0,1,1,2)), ("F",(1,2,1,2)), ("G",(2,4,1,2))],
    3: [("A",(0,2,0,1)), ("B",(2,3,0,1)), ("C",(3,4,0,1)),
        ("D",(0,1,1,2)), ("E",(1,2,1,2)), ("F",(2,3,1,2)), ("G",(3,4,1,2))],
    4: [("A",(0,1,0,1)), ("B",(1,2,0,1)), ("C",(2,3,0,1)), ("D",(3,4,0,1)),
        ("E",(0,1,1,2)), ("F",(1,3,1,2)), ("G",(3,4,1,2))],
    5: [("A",(0,1,0,1)), ("B",(1,2,0,1)), ("C",(2,3,0,1)), ("D_msk",(3,4,0,1)),
        ("D_port",(0,1,1,2)), ("E_ffpe",(1,2,1,2)), ("F_roc",(2,3,1,2)), ("G_calib",(3,4,1,2))],
    6: [("A",(0,1,0,1)), ("B",(1,2,0,1)), ("C",(2,3,0,1)), ("D",(3,4,0,1)),
        ("E",(0,1,1,2)), ("F",(1,3,1,2)), ("G",(3,4,1,2))],
}
files = {
    1: "Fig1_discovery_axis.png",
    2: "Fig2_fusion_mechanism.png",
    3: "Fig3_epigenetic.png",
    4: "Fig4_sc_validation.png",
    5: "Fig5_survival_portability.png",
    6: "Fig6_reflex_translation.png",
}

for fig_no, fname in files.items():
    img = Image.open(os.path.join(ROOT_OUT, fname))
    W, H = img.size
    # account for top suptitle margin (~7%) and tight bbox already trimmed
    title_h = int(H * 0.07)
    grid_h = H - title_h
    col_w = W / 4
    row_h = grid_h / 2
    for label,(c0,c1,r0,r1) in panel_specs[fig_no]:
        left = int(c0 * col_w)
        right= int(c1 * col_w)
        top  = title_h + int(r0 * row_h)
        bot  = title_h + int(r1 * row_h)
        crop = img.crop((max(0,left-4), max(0,top-4), min(W,right+4), min(H,bot+4)))
        out = f"{ASSETS}/Fig{fig_no}_panel_{label}.png"
        # downsize to smaller thumbnail
        crop.thumbnail((900, 900))
        crop.save(out, optimize=True)
print("Panel thumbnails built.")
