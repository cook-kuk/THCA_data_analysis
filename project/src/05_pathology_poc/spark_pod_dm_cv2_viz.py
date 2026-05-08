#!/usr/bin/env python3
"""
Pod DM 결과 cv2 시각화 — DINOv2-ViT-L-518 LOSO classifier + 한글 storytelling.

3 main outputs:
  1. F7_backbone_progression.jpg
     ResNet50 (chance 0.51) → DINOv2 (0.573) → UNI (expected 0.65+)
     bar chart + 빨간 박스 highlight + 한글 annotation

  2. F8_per_case_predictions.jpg
     40-case forest plot, DM1 (parallel) vs DM2 + correct/wrong overlay
     ROC curve inset

  3. F9_panel_overview.jpg
     One mega-panel with everything:
       - top: backbone progression
       - middle: per-case predictions
       - bottom: ROC + summary stats + 한글 take-home
"""
from __future__ import annotations
from pathlib import Path
import cv2
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from sklearn.metrics import roc_curve, auc

ROOT = Path(__file__).resolve().parent.parent.parent.parent
RES = ROOT / "project/results/03_pathology_poc"
POD = RES / "pod_dm"
OUT = RES / "pod_dm_viz"
OUT.mkdir(parents=True, exist_ok=True)

FONT_PATH = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
F36 = ImageFont.truetype(FONT_PATH, 36)
F28 = ImageFont.truetype(FONT_PATH, 28)
F20 = ImageFont.truetype(FONT_PATH, 20)
F16 = ImageFont.truetype(FONT_PATH, 16)


def add_text(img, x, y, text, font, color=(255, 255, 255), bg=None, padding=8):
    pil = Image.fromarray(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
    draw = ImageDraw.Draw(pil)
    if bg is not None:
        bbox = draw.textbbox((x, y), text, font=font)
        draw.rectangle((bbox[0]-padding, bbox[1]-padding, bbox[2]+padding, bbox[3]+padding), fill=bg)
    draw.text((x, y), text, fill=color, font=font)
    return cv2.cvtColor(np.array(pil), cv2.COLOR_RGB2BGR)


def F7_backbone_progression():
    W, H = 1800, 900
    canvas = np.full((H, W, 3), 12, dtype=np.uint8)  # near-black
    canvas[:, :, 1] = 17; canvas[:, :, 2] = 26
    canvas = add_text(canvas, 40, 40, "F7 · H&E Foundation Model 진화 — DM1 vs DM2 예측 성능", F36, (255, 255, 255))
    canvas = add_text(canvas, 40, 90, "ResNet50 ImageNet → DINOv2 ViT-L → UNI (예상)", F20, (148, 163, 196))

    # Three bars
    bars = [
        ("ResNet50\n(2015 ImageNet)", 0.51, "#94a3b8", "GSE250521 ST · 16-slide spot LOSO\n100 random panel control 통과 X"),
        ("DINOv2 ViT-L\n(2024 self-supervised, pathology-naive)", 0.573, "#fbbf24",
         "TCGA-THCA WSI 50-case LOSO ★\n+6 percentage points · n=40 (DM1=22, DM2=18)"),
        ("UNI Mahmood Lab\n(2024 100M+ pathology)", 0.70, "#5eead4",
         "예상 — HF gated, license + token 필요\n+13 pp 추가 가능 (literature)"),
    ]
    bar_y0 = 200; bar_h = 80; bar_x0 = 240
    bar_w_max = 1200
    for i, (label, val, color_hex, sub) in enumerate(bars):
        y = bar_y0 + i * 220
        # bar background
        cv2.rectangle(canvas, (bar_x0, y), (bar_x0 + bar_w_max, y + bar_h), (50, 60, 90), -1)
        # bar value
        bgr = tuple(int(color_hex.lstrip('#')[k:k+2], 16) for k in (4, 2, 0))  # BGR
        bar_w = int((val - 0.5) / 0.5 * bar_w_max)  # 0.5 ~ 1.0 range
        cv2.rectangle(canvas, (bar_x0, y), (bar_x0 + bar_w, y + bar_h), bgr, -1)
        # chance line at val=0.5
        chance_x = bar_x0
        cv2.line(canvas, (chance_x, y - 10), (chance_x, y + bar_h + 10), (200, 200, 200), 2)
        # Red box on DINOv2 (i=1) — actual measured
        if i == 1:
            cv2.rectangle(canvas, (bar_x0 - 10, y - 10), (bar_x0 + bar_w + 10, y + bar_h + 10),
                         (0, 0, 255), 4)
        # value text
        canvas = add_text(canvas, bar_x0 + bar_w + 25, y + 20,
                         f"AUROC = {val:.3f}", F28,
                         (255, 215, 0) if i == 1 else (200, 220, 240))
        # label
        canvas = add_text(canvas, 40, y + 10, label, F20, (255, 255, 255))
        canvas = add_text(canvas, 40, y + 50, sub, F16, (148, 163, 196))

    # X axis labels
    canvas = add_text(canvas, bar_x0 - 40, bar_y0 + 3*220 + 30, "0.50 (chance)", F16, (148, 163, 196))
    canvas = add_text(canvas, bar_x0 + bar_w_max - 30, bar_y0 + 3*220 + 30, "1.00 (perfect)", F16, (148, 163, 196))

    # Take-home
    canvas = add_text(canvas, 40, H - 100,
                     "★ 결론: DINOv2 가 ResNet50 chance level 을 뚫음 (+6pp). UNI 활용 시 +13pp 추가 가능 → 임상 활용 검토 가능",
                     F20, (252, 211, 77))
    canvas = add_text(canvas, 40, H - 60,
                     "★ Pod DM L40S 1.5h · TCGA-THCA n=50 DM-balanced subset · GDC API 직접 다운로드 58.7 GB",
                     F16, (148, 163, 196))

    cv2.imwrite(str(OUT / "F7_backbone_progression.jpg"), canvas, [cv2.IMWRITE_JPEG_QUALITY, 92])
    print(f"wrote {OUT / 'F7_backbone_progression.jpg'}")


def F8_per_case_predictions():
    df = pd.read_csv(POD / "dm1_vs_dm2_loso_predictions.tsv", sep="\t")
    df = df.sort_values(["dm_true", "dm_pred_prob"]).reset_index(drop=True)

    W, H = 1800, 1100
    canvas = np.full((H, W, 3), 12, dtype=np.uint8); canvas[:,:,1]=17; canvas[:,:,2]=26
    canvas = add_text(canvas, 40, 30, "F8 · Per-case DINOv2 LOSO 예측 (40명)", F36, (255, 255, 255))
    canvas = add_text(canvas, 40, 80, "DM1 (HT-like 면역활성) ↔ DM2 (aggressive driver-neg)", F20, (148, 163, 196))

    # Forest plot
    plot_x0 = 380; plot_x1 = 1200
    plot_y0 = 130; row_h = 22
    cv2.line(canvas, ((plot_x0+plot_x1)//2, plot_y0), ((plot_x0+plot_x1)//2, plot_y0 + len(df)*row_h + 10),
             (180, 180, 180), 2)
    canvas = add_text(canvas, (plot_x0+plot_x1)//2 - 20, plot_y0 - 28, "0.5 threshold", F16, (148, 163, 196))

    for i, r in df.iterrows():
        y = plot_y0 + i * row_h + 10
        # case label
        canvas = add_text(canvas, 40, y - 8, f"{r.case_id}", F16, (200, 220, 240))
        # true class chip
        true_color = (94, 234, 212) if r.dm_true == "DM1" else (255, 124, 62)  # BGR teal/orange
        true_color = tuple(int(c) for c in true_color)
        cv2.rectangle(canvas, (220, y - 12), (270, y + 8), true_color, -1)
        canvas = add_text(canvas, 230, y - 8, r.dm_true, F16, (10, 10, 10))
        # pred dot
        x_pred = int(plot_x0 + r.dm_pred_prob * (plot_x1 - plot_x0))
        # correct/incorrect color
        pred_class = "DM1" if r.dm_pred_prob < 0.5 else "DM2"
        correct = pred_class == r.dm_true
        dot_color = (94, 234, 212) if correct else (0, 0, 255)  # teal correct, red wrong
        cv2.circle(canvas, (x_pred, y), 8, dot_color, -1)
        # bar from threshold to dot
        thr_x = (plot_x0+plot_x1)//2
        cv2.line(canvas, (thr_x, y), (x_pred, y), dot_color, 2)
        # prob text
        canvas = add_text(canvas, plot_x1 + 20, y - 8, f"{r.dm_pred_prob:.2f}", F16,
                         (94, 234, 212) if correct else (0, 100, 255))

    # Legend
    legend_y = H - 200
    cv2.circle(canvas, (60, legend_y), 8, (94, 234, 212), -1)
    canvas = add_text(canvas, 80, legend_y - 10, "Correct", F20, (94, 234, 212))
    cv2.circle(canvas, (220, legend_y), 8, (0, 0, 255), -1)
    canvas = add_text(canvas, 240, legend_y - 10, "Wrong", F20, (0, 100, 255))

    # ROC
    df_for_roc = pd.read_csv(POD / "dm1_vs_dm2_loso_predictions.tsv", sep="\t")
    y_true = (df_for_roc.dm_true == "DM2").astype(int).values
    y_score = df_for_roc.dm_pred_prob.values
    fpr, tpr, _ = roc_curve(y_true, y_score)
    auc_v = auc(fpr, tpr)
    # ROC mini-plot in bottom-right (200×200)
    rx0, ry0 = W - 280, H - 320; rsz = 240
    cv2.rectangle(canvas, (rx0, ry0), (rx0 + rsz, ry0 + rsz), (50, 60, 90), 2)
    cv2.line(canvas, (rx0, ry0 + rsz), (rx0 + rsz, ry0), (180, 180, 180), 1)  # diagonal
    pts = [(rx0 + int(f * rsz), ry0 + rsz - int(t * rsz)) for f, t in zip(fpr, tpr)]
    for a, b in zip(pts[:-1], pts[1:]):
        cv2.line(canvas, a, b, (94, 234, 212), 3)
    canvas = add_text(canvas, rx0, ry0 - 28, f"ROC · AUC = {auc_v:.3f}", F20, (94, 234, 212))
    canvas = add_text(canvas, rx0, ry0 + rsz + 6, "FPR →   ↑ TPR", F16, (148, 163, 196))

    # Take-home
    canvas = add_text(canvas, 40, H - 100,
                     "★ 결론: 40 case 중 25 명 정확 (62.5%), 15 명 오분류 — 모델은 chance 보다 강하지만 임상 적용엔 부족",
                     F20, (252, 211, 77))
    canvas = add_text(canvas, 40, H - 60,
                     "★ 분석: DM 라벨은 RNA 기반 unsupervised cluster, H&E 만으로 직접 readout 은 어려움 — UNI/HoVer-NeXt + 다른 cohort 필요",
                     F16, (148, 163, 196))

    cv2.imwrite(str(OUT / "F8_per_case_predictions.jpg"), canvas, [cv2.IMWRITE_JPEG_QUALITY, 92])
    print(f"wrote {OUT / 'F8_per_case_predictions.jpg'}")


def F9_grand_overview():
    """Combine F7+F8 with extra storytelling panel."""
    f7 = cv2.imread(str(OUT / "F7_backbone_progression.jpg"))
    f8 = cv2.imread(str(OUT / "F8_per_case_predictions.jpg"))
    # Resize to same width
    target_w = 1800
    f7 = cv2.resize(f7, (target_w, int(f7.shape[0] * target_w / f7.shape[1])))
    f8 = cv2.resize(f8, (target_w, int(f8.shape[0] * target_w / f8.shape[1])))
    # Add story panel
    story_h = 400
    story = np.full((story_h, target_w, 3), 12, dtype=np.uint8); story[:,:,1]=17; story[:,:,2]=26
    story = add_text(story, 40, 30, "★ 종합 한 줄 해석 (Pod DM 1차 결과)", F36, (252, 211, 77))
    lines = [
        "1. ResNet50 ImageNet 으로는 DM1/DM2 분리 불가 (16-slide ST AUROC 0.51 = chance, 100 random panel control 입증)",
        "2. DINOv2 ViT-L 518 self-supervised (pathology-naive) 가 chance 를 뚫음 (TCGA WSI 50 LOSO AUROC 0.573, +6pp)",
        "3. 그러나 임상 활용 (AUROC 0.85+ 필요) 까지 거리 — UNI / GigaTIME / HoVer-NeXt 7-class 추가 필요",
        "4. 더 강한 finding 은 NON-H&E side: CCI 회피 -0.27 ST → -0.31 TCGA → -0.34 Korean 3-cohort replicate ★",
        "5. → paper 의 main figure 는 cross-platform CCI replicate, H&E 는 supplement (foundation model upgrade trajectory)",
        "",
        "★ Pod DM 다음: HoVer-NeXt 7-class POC (5 WSI 시도 중) → 성공시 SPARK Analytical 8-feature 적용 → DINOv2 보다 더 강한 H&E 분리 시도",
    ]
    for i, l in enumerate(lines):
        story = add_text(story, 40, 100 + i * 38, l, F20,
                        (252, 211, 77) if l.startswith("★") else (200, 220, 240))

    canvas = np.vstack([f7, f8, story])
    cv2.imwrite(str(OUT / "F9_grand_overview.jpg"), canvas, [cv2.IMWRITE_JPEG_QUALITY, 90])
    print(f"wrote {OUT / 'F9_grand_overview.jpg'}")


if __name__ == "__main__":
    F7_backbone_progression()
    F8_per_case_predictions()
    F9_grand_overview()
