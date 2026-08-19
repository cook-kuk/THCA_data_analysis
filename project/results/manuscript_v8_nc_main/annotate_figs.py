"""Annotate Fig 1-6 + atlas with Korean callouts using PIL + OpenCV.

Each figure gets:
  - 우측 또는 하단에 한글 strip (figure-wide 결론)
  - 핵심 panel 위에 색상 highlight box
  - 핵심 stat 옆 화살표 + 한글 콜아웃
"""
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageFilter
import cv2

ROOT = "/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main"
FONT_PATH = "/home/seungho/.local/share/fonts/NanumGothic-Regular.ttf"

# Palette aligned with figure style
DM1_RED   = (185, 28, 28)
AMBER     = (180, 83, 9)
INK       = (15, 23, 42)
LIGHT_BG  = (250, 250, 248)
HIGHLIGHT = (255, 235, 167)   # soft yellow
ACCENT    = (217, 119, 6)     # vivid amber
CRIMSON   = (220, 38, 38)
ROYAL     = (37, 99, 235)
WHITE     = (255, 255, 255)

def font(size, bold=False):
    # Use NanumGothic for Korean text; fall back default
    try:
        return ImageFont.truetype(FONT_PATH, size)
    except Exception:
        return ImageFont.load_default()

def callout(draw, xy, text, size=20, fc=AMBER, tc=WHITE, pad=10, max_width=None):
    """Draw a filled rounded rectangle with Korean text inside."""
    f = font(size, bold=True)
    # Compute text bbox
    bbox = draw.textbbox((0,0), text, font=f)
    tw, th = bbox[2]-bbox[0], bbox[3]-bbox[1]
    x, y = xy
    box = [x, y, x + tw + pad*2, y + th + pad*2]
    draw.rounded_rectangle(box, radius=8, fill=fc, outline=None)
    draw.text((x+pad, y+pad-bbox[1]), text, font=f, fill=tc)
    return box

def arrow(img, p1, p2, color=ACCENT, thickness=4):
    """Draw an arrow on cv2 array (BGR)."""
    cv2.arrowedLine(img, p1, p2, color[::-1], thickness, tipLength=0.18, line_type=cv2.LINE_AA)

def highlight_rect(img, p1, p2, color=ACCENT, thickness=4, alpha=0.0):
    """Draw a rectangle (no fill, colored border) on cv2 array."""
    cv2.rectangle(img, p1, p2, color[::-1], thickness, lineType=cv2.LINE_AA)

def add_bottom_strip(pil_img, lines, strip_height=110, side_label="결론"):
    """Append a colored strip at the bottom with Korean conclusion lines."""
    W, H = pil_img.size
    new_h = H + strip_height + 8
    canvas = Image.new("RGB", (W, new_h), LIGHT_BG)
    canvas.paste(pil_img, (0, 0))
    draw = ImageDraw.Draw(canvas)
    # Top divider
    draw.line([(0, H+4), (W, H+4)], fill=(203, 213, 225), width=2)
    # Strip background
    strip_box = [0, H+8, W, new_h]
    draw.rectangle(strip_box, fill=(254, 252, 232))
    # Side label badge
    f_label = font(18, bold=True)
    label_box = draw.textbbox((0,0), side_label, font=f_label)
    pad = 12
    badge = [16, H+22, 16 + (label_box[2]-label_box[0]) + pad*2, H+22 + (label_box[3]-label_box[1]) + pad*2]
    draw.rounded_rectangle(badge, radius=8, fill=DM1_RED)
    draw.text((badge[0]+pad, badge[1]+pad-label_box[1]), side_label, font=f_label, fill=WHITE)
    # Body lines
    f_body = font(20, bold=False)
    y_text = H + 28
    x_text = badge[2] + 18
    max_w = W - x_text - 24
    for ln in lines:
        # word-wrap (Korean — break at spaces)
        words = ln.split(" ")
        cur = ""
        for w in words:
            test = (cur + " " + w).strip()
            bw = draw.textbbox((0,0), test, font=f_body)[2]
            if bw > max_w and cur:
                draw.text((x_text, y_text), cur, font=f_body, fill=INK)
                y_text += 28
                cur = w
            else:
                cur = test
        if cur:
            draw.text((x_text, y_text), cur, font=f_body, fill=INK)
            y_text += 28
    return canvas

def annotate(src, dst, panel_highlights=None, callouts=None, bottom_lines=None, side_label="결론"):
    """One-stop helper: load src, draw highlights + callouts via cv2, append bottom strip via PIL, save dst."""
    bgr = cv2.imread(src)
    if bgr is None:
        print(f"⨯ cannot read {src}")
        return
    # 1. cv2 highlights
    if panel_highlights:
        for (p1, p2, color) in panel_highlights:
            highlight_rect(bgr, p1, p2, color=color, thickness=5)
    # 2. cv2 arrows
    if callouts:
        for c in callouts:
            if "arrow" in c:
                arrow(bgr, c["arrow"][0], c["arrow"][1], color=c.get("color", ACCENT), thickness=c.get("thickness", 4))
    # cv2 → PIL
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    pil = Image.fromarray(rgb)
    # 3. PIL callout text boxes
    if callouts:
        draw = ImageDraw.Draw(pil)
        for c in callouts:
            if "text" in c and "at" in c:
                callout(draw, c["at"], c["text"], size=c.get("size", 22),
                        fc=c.get("box_color", AMBER), tc=c.get("text_color", WHITE), pad=c.get("pad", 10))
    # 4. bottom Korean strip
    if bottom_lines:
        h = 30 + len(bottom_lines)*32
        pil = add_bottom_strip(pil, bottom_lines, strip_height=h, side_label=side_label)
    pil.save(dst, quality=92, optimize=True)
    print(f"✓ {os.path.basename(dst)}")


# ===== Read figure dimensions to make coordinates portable =====
def dims(path):
    img = cv2.imread(path)
    return img.shape[1], img.shape[0]  # W, H

# ============ FIG 1 (Discovery axis) ============
W, H = dims(f"{ROOT}/Fig1_discovery_axis.png")
# Estimated panel positions in pixel coords (figsize=(18, 12.6), dpi=180 -> ~3240×2268)
# Layout: top row 4 cols (a,b,c,d) at roughly y=160..1100; bottom row 3 cols at y=1230..2100
# We'll do high-level callouts pointing at key regions
annotate(
    f"{ROOT}/Fig1_discovery_axis.png",
    f"{ROOT}/Fig1_annotated.png",
    panel_highlights=[
        # KM panel (g) - bottom-right roughly
        ((int(W*0.50), int(H*0.55)), (int(W*0.985), int(H*0.92)), ACCENT),
        # UMAP panel (c) - top mid
        ((int(W*0.49), int(H*0.07)), (int(W*0.735), int(H*0.45)), CRIMSON),
    ],
    callouts=[
        {"at": (int(W*0.50), int(H*0.95)), "text": "★ 핵심: DM1 (28.4 %) / DM2 군 명확 분리 + 생존 곡선 가파른 분리",
         "size": 22, "box_color": DM1_RED},
    ],
    bottom_lines=[
        "그림 1 — DM1/DM2 발견 축.  TCGA n=504 종양에서 8-유전자 RAI 패널 (TSHR · TPO · TG · PAX8 · NKX2-1 · FOXE1 · DIO1 · SLC5A5) 에",
        "KMeans (k=2) 적용 시 두 군이 PCA 공간에서 명확히 분리 (DM1 28.4 %, DM2 71.6 %).  핵심: BRAF · TERT · KRAS · NRAS · HRAS 모든 driver 의",
        "단일 feature AUC ≈ 0.5 — DM1 축은 어떤 driver 변이와도 평행이 아닌 새로운 분자 축.  Pan-genome top-5000 ARI = 0.92 ≈ TIERA67 (0.90) 로",
        "패널이 우연 발견이 아닌 더 큰 transcriptional axis 의 의미 있는 lens 임을 입증.  KM 곡선의 가파른 분리가 임상 의의의 직접 증거.",
    ],
    side_label="핵심 메시지"
)

# ============ FIG 2 (Fusion mechanism) ============
W, H = dims(f"{ROOT}/Fig2_fusion_mechanism.png")
annotate(
    f"{ROOT}/Fig2_fusion_mechanism.png",
    f"{ROOT}/Fig2_annotated.png",
    panel_highlights=[
        # Panel c (fusion enrichment) top-mid
        ((int(W*0.49), int(H*0.07)), (int(W*0.735), int(H*0.45)), CRIMSON),
        # Panel f (FVPTC mosaic) bottom-mid
        ((int(W*0.27), int(H*0.55)), (int(W*0.50), int(H*0.92)), ACCENT),
    ],
    callouts=[
        {"at": (int(W*0.50), int(H*0.95)), "text": "★ 핵심: DM1 76.8 % 융합+ (OR 7.41) + FVPTC 농축 (OR 17.9)",
         "size": 22, "box_color": DM1_RED},
    ],
    bottom_lines=[
        "그림 2 — DM1 융합 기전.  Xing 2014 BRAF−/TERT− dark matter 180 명 중 131 (73 %) 이 DM1/DM2 로 재분류.  핵심 발견:",
        "DM1 종양의 76.8 % (63/82) 가 kinase 융합 양성 (RET 33 · NTRK 10 · ALK 4 · BRAF 5).  Fisher OR = 7.41 [4.38, 12.55], p = 1.9 × 10⁻¹³.",
        "TCGA RET-양성 33 종양 중 27 (81.8 %) 가 DM1 으로 호출 — RNA-first reflex 알고리즘의 직접 정당화.  FVPTC 조직형 OR = 17.9 (p = 3 × 10⁻³¹) —",
        "패널은 RAI 생물학 기반으로 선정되었으나 조직형 substructure 까지 회복.  MSK-IMPACT 117 명 진행성 코호트에서 동일한 fusion 분포가 재현.",
    ],
    side_label="핵심 메시지"
)

# ============ FIG 3 (Epigenetic silencing) ============
W, H = dims(f"{ROOT}/Fig3_epigenetic.png")
annotate(
    f"{ROOT}/Fig3_epigenetic.png",
    f"{ROOT}/Fig3_annotated.png",
    panel_highlights=[
        # Heatmap panel (a) - top-left spans 2 cols
        ((int(W*0.04), int(H*0.06)), (int(W*0.475), int(H*0.46)), CRIMSON),
        # Landa convergence panel (g) - bottom right
        ((int(W*0.755), int(H*0.55)), (int(W*0.985), int(H*0.92)), ACCENT),
    ],
    callouts=[
        {"at": (int(W*0.50), int(H*0.95)), "text": "★ 핵심: TPO d=2.30 (p=1.9×10⁻¹⁸) + Landa 5/8 overlap = reverse-causality lock",
         "size": 22, "box_color": DM1_RED},
    ],
    bottom_lines=[
        "그림 3 — DM1 후성유전 침묵.  HM450 promoter β heatmap (TCGA n=503) 에서 DM1 종양의 8-유전자가 일관되게 메틸화.",
        "TPO Cohen's d = 2.30 (p = 1.9 × 10⁻¹⁸) — 갑상선암 epigenetics 에서 보기 드문 effect size.  DM1 mean β = 0.385 vs DM2 0.253 (+52 %).",
        "5 코호트 1,287 명 MAPK output × Panel-8 pooled Spearman ρ = −0.327 [−0.376, −0.278] — MAPK 활성이 높을수록 panel 발현 감소.",
        "Landa 2016 ATC silenced gene list 와 우리 8-유전자 중 5 개 (TG · TSHR · TPO · PAX8 · DIO1) 직접 1 : 1 overlap = convergent biology, 독립 설계 검증.",
    ],
    side_label="핵심 메시지"
)

# ============ FIG 4 (sc validation) ============
W, H = dims(f"{ROOT}/Fig4_sc_validation.png")
annotate(
    f"{ROOT}/Fig4_sc_validation.png",
    f"{ROOT}/Fig4_annotated.png",
    panel_highlights=[
        # UMAP panel (a)
        ((int(W*0.04), int(H*0.07)), (int(W*0.255), int(H*0.46)), CRIMSON),
        # Pseudotime panel (f) - bottom mid spans 2
        ((int(W*0.27), int(H*0.55)), (int(W*0.735), int(H*0.92)), ACCENT),
    ],
    callouts=[
        {"at": (int(W*0.50), int(H*0.95)), "text": "★ 핵심: thyrocyte 본질적 신호 + 4 cell type 단조 변화 |ρ|≥0.95",
         "size": 22, "box_color": DM1_RED},
    ],
    bottom_lines=[
        "그림 4 — 단일세포 외부 검증.  Lu 2023 GSE193581 thyrocyte UMAP (KRT8∩KRT19∩EPCAM 필터, n=14,624 cells) 에서 DM score gradient",
        "가 stromal/immune 이 아닌 thyrocyte 자체에서 명확.  Pu 2021 paired 6 환자 모두 per-patient r = 0.798–0.886 (Bonferroni p < 10⁻¹⁰).",
        "10-decile composition pseudotime 에서 4 compartment 단조 변화 — Malignant ↓ Epithelial ↑ Myeloid ↓ Endothelial ↑, 모두 |ρ| ≥ 0.95.",
        "TCGA n=513 + Lee n=632 양쪽에서 동일하게 재현.  Reviewer 의 1차 공격 (\"bulk RNA 신호가 stromal/immune 혼동 아니냐\") 직접 차단.",
    ],
    side_label="핵심 메시지"
)

# ============ FIG 5 (Survival + portability) ============
W, H = dims(f"{ROOT}/Fig5_survival_portability.png")
annotate(
    f"{ROOT}/Fig5_survival_portability.png",
    f"{ROOT}/Fig5_annotated.png",
    panel_highlights=[
        # Forest panel (a)
        ((int(W*0.04), int(H*0.07)), (int(W*0.255), int(H*0.46)), CRIMSON),
        # FFPE vs FF (f)
        ((int(W*0.27), int(H*0.55)), (int(W*0.495), int(H*0.92)), ACCENT),
    ],
    callouts=[
        {"at": (int(W*0.50), int(H*0.95)), "text": "★ 핵심: pooled HR=2.53 [1.31, 4.89] I²=0% + FFPE 호환 (KS p=0.44)",
         "size": 22, "box_color": DM1_RED},
    ],
    bottom_lines=[
        "그림 5 — 통합 생존 + 임상 portability.  TCGA HR 2.30 + MSK HR 2.67 의 DerSimonian-Laird pooled HR = 2.53 [1.31, 4.89], I² = 0 %.",
        "Age, stage III/IV, TERT promoter mut 보정 후에도 DM1 hazard 독립 유지.  FFPE (Lee n=632) vs FF (TCGA n=504) Kolmogorov-Smirnov p = 0.44 —",
        "임상실 일상 표본인 paraffin-embedded 표본으로도 검사 가능.  Time-dep ROC 1y AUC = 0.83 · 3y = 0.78 · 5y = 0.72 (DM1 + age + stage, IPCW).",
        "5-year calibration 의 Hosmer-Lemeshow p = 0.31 (NS) 로 모델이 잘 calibrate 됨.  단순 분자 분류가 아닌 임상 deploy 가능 axis.",
    ],
    side_label="핵심 메시지"
)

# ============ FIG 6 (Reflex pathway) ============
W, H = dims(f"{ROOT}/Fig6_reflex_translation.png")
annotate(
    f"{ROOT}/Fig6_reflex_translation.png",
    f"{ROOT}/Fig6_annotated.png",
    panel_highlights=[
        # Reflex flowchart (a)
        ((int(W*0.04), int(H*0.07)), (int(W*0.255), int(H*0.46)), CRIMSON),
        # post-RAI box (d)
        ((int(W*0.755), int(H*0.07)), (int(W*0.985), int(H*0.46)), ACCENT),
    ],
    callouts=[
        {"at": (int(W*0.50), int(H*0.95)), "text": "★ 핵심: 1,000 PTC 당 48 selpercatinib 후보 + post-RAI 종양 = DM1 state",
         "size": 22, "box_color": DM1_RED},
    ],
    bottom_lines=[
        "그림 6 — 임상 reflex 경로 + translation.  RNA → 8-유전자 패널 → DM1 호출 → fusion test (RET/NTRK/ALK/BRAF) → selpercatinib (Wirth 2020 NEJM",
        "LIBRETTO-001 ORR 79 %) / larotrectinib 후보군 식별.  추정 yield: 1,000 PTC 환자 당 48 명의 selpercatinib 후보.  ATA 2015/2025 risk tier mosaic 에서",
        "DM1 은 중간 위험군 (RAI 결정의 가장 큰 불확실 zone) 에 과대표현.  GSE151179 post-RAI 종양은 DM1 transcriptional state 와 일치 (Cohen's d = −1.01,",
        "MW p = 1 × 10⁻⁴) — DM1 = pre-existing RAI-resistant phenotype.  Decision-curve 에서 DM1-stratified strategy 가 treat-all/none 을 dominate.",
    ],
    side_label="핵심 메시지"
)

# ============ External Validation Atlas (overlay) ============
W, H = dims(f"{ROOT}/External_Validation_Atlas.png")
annotate(
    f"{ROOT}/External_Validation_Atlas.png",
    f"{ROOT}/External_Validation_Atlas_annotated.png",
    panel_highlights=None,
    callouts=[
        {"at": (int(W*0.36), int(H*0.965)), "text": "★ 14/15 외부 코호트 direction-consistent · 1 attenuation (GSE286332) = 모델 예측 HT route",
         "size": 22, "box_color": DM1_RED},
    ],
    bottom_lines=[
        "External Validation Atlas — 15 + TCGA = 16 코호트 카드.  6 modality (bulk RNA · methylation · single-cell · protein · drug · pan-cancer) ·",
        "4 인종 (한국 · 동아시아 · 유럽 · 다인종) · 3 질병 단계 (primary PTC · 진행성 PDTC/ATC · post-RAI refractory).  핵심: 14 / 15 외부 코호트",
        "direction-consistent.  1 개 attenuation (GSE286332 PTC+HT) 은 우연이 아닌 사전 가설 (\"HT route 에서는 MAPK 가 driver 가 아니므로 attenuated\")",
        "의 직접 검증.  Reviewer 의 가장 흔한 공격 \"TCGA cherry-pick 아니냐\" 의 1차 차단선이며 NC-tier 진입의 핵심 강점.",
    ],
    side_label="핵심 메시지"
)

# ============ EV master forest (canonical) ============
W, H = dims(f"{ROOT}/../../papers_hub_2026_05_04/assets/paper1_nc/EV_master_forest.png")
src = "/home/seungho/personal/THCA_data_analysis/project/papers_hub_2026_05_04/assets/paper1_nc/EV_master_forest.png"
annotate(
    src,
    f"{ROOT}/EV_master_forest_annotated.png",
    panel_highlights=None,
    callouts=[
        {"at": (int(W*0.32), int(H*0.94)), "text": "★ 평균 Cohen's d = 2.81 · median 2.37 · 14 entry 모두 ≥ 1.56",
         "size": 22, "box_color": DM1_RED},
    ],
    bottom_lines=[
        "Master Cross-Cohort Forest — 14 entries × 11 distinct cohorts.  Lee 2024 within-cohort KMeans d = 5.93 (all) / 5.48 (tumour-only),",
        "Landa 2016 PDTC+ATC d = 4.08, GSE33630 ATC vs PTC d = 3.91, Lee 2024 ATC vs Normal d = 3.06, Mun 2025 protein d = 2.86, GSE65144 d = 2.79,",
        "K2 within-cohort d = 1.94, Pu 2021 sc d = 1.56.  평균 d = 2.81, median = 2.37 — 갑상선암 분자 아형 문헌에서 보기 드문 폭과 깊이.",
        "Modality 분포: Korean RNA 5 · Western array 4 · Western protein 2 · Western RNA 3 — 인종 / platform / modality 균등 배분 cross-validation.",
    ],
    side_label="핵심 메시지"
)

# ============ EV per-gene matrix ============
src = "/home/seungho/personal/THCA_data_analysis/project/papers_hub_2026_05_04/assets/paper1_nc/EV_pergene_matrix.png"
W, H = dims(src)
annotate(
    src,
    f"{ROOT}/EV_pergene_matrix_annotated.png",
    panel_highlights=None,
    callouts=[
        {"at": (int(W*0.30), int(H*0.94)), "text": "★ 80 / 80 cell direction-consistent — 단일 유전자 의존성 없음",
         "size": 22, "box_color": DM1_RED},
    ],
    bottom_lines=[
        "Per-Gene × Cohort × Contrast Cohen's d Matrix — 8 panel 유전자 × 10 contrast × 4 cohort = 80 cell.",
        "단 하나의 cell 도 방향이 뒤집히지 않음.  TPO 가 최대 d = +3.53 (Lee 2024 ATC vs Normal); 가장 약한 SLC5A5 도 모든 cohort 에서",
        "동일 방향 유지.  핵심 차단점: \"8-유전자 중 일부만 작동\" 또는 \"단일 cohort 효과\" 가설을 동시에 reject.  진정한 multi-gene",
        "multi-cohort 분자 axis.  Reviewer 의 panel-selection-bias 공격에 대한 가장 강한 답.",
    ],
    side_label="핵심 메시지"
)

print("\n=== Annotation complete ===")
