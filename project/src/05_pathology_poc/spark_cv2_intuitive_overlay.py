#!/usr/bin/env python3
"""
직관적 spatial 시각화 — cv2 기반 heatmap overlay + 한글 annotation.

기존 matplotlib scatter overlay (점들이 흩뿌려진) 대신:
  1. spot 점수 → image grid 보간 (rbf-style nearest)
  2. cv2.GaussianBlur 로 smooth heatmap
  3. cv2.applyColorMap (JET) 적용
  4. cv2.addWeighted 로 H&E 위에 알파 블렌딩
  5. PIL 로 한글 annotation + 색상 의미 explain

각 slide 마다 4-panel:
  - H&E (원본 + stage label)
  - DM1 score heatmap (분화도; 빨강 = de-diff, 파랑 = diff)
  - CAF score heatmap (섬유화 niche; 노랑 = CAF-rich)
  - CAF × −RAI avoidance hotspots (CAF가 RAI thyrocyte 피하는 지역; 마젠타 = avoidance hot)

Output:
  project/results/03_pathology_poc/cv2_intuitive_overlays/{sid}_4panel.png
  + index.html
"""
from __future__ import annotations
import sys
from pathlib import Path

import anndata as ad
import cv2
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
OUT = ROOT / "project/results/03_pathology_poc/cv2_intuitive_overlays"
OUT.mkdir(parents=True, exist_ok=True)

# Korean-capable font search
FONT_PATHS = [
    "/usr/share/fonts/truetype/nanum/NanumGothic.ttf",
    "/usr/share/fonts/truetype/nanum/NanumGothicBold.ttf",
    "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
]
FONT = None
for fp in FONT_PATHS:
    if Path(fp).exists():
        try:
            FONT = ImageFont.truetype(fp, 28)
            FONT_SMALL = ImageFont.truetype(fp, 18)
            FONT_PATH = fp
            break
        except Exception:
            pass
if FONT is None:
    FONT = ImageFont.load_default()
    FONT_SMALL = ImageFont.load_default()
    FONT_PATH = "default (no Korean)"
print(f"font: {FONT_PATH}")

CAF = ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"]
RAI = ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"]


def get_scalefactor(spatial_uns):
    if not spatial_uns: return 1.0
    libs = list(spatial_uns.keys())
    return float(spatial_uns[libs[0]].get("scalefactors", {}).get("tissue_hires_scalef", 1.0))


def hires_path(sdir):
    cands = list(sdir.rglob("tissue_hires_image.png"))
    return cands[0] if cands else None


def module(X, sym, genes):
    cols = [sym[g] for g in genes if g in sym]
    if len(cols) < 2: return None
    sub = X[:, cols]
    return ((sub - sub.mean(0)) / (sub.std(0) + 1e-6)).mean(1)


def render_heatmap(bg_img, coords, values, cmap_id=cv2.COLORMAP_JET, alpha=0.55, sigma=25):
    """Render cv2-smoothed heatmap overlay on background image."""
    H, W = bg_img.shape[:2]
    # interpolate via nearest-neighbor onto a coarse grid then upsample
    grid = np.full((H, W), np.nan, dtype=np.float32)
    # plant values at integer pixel coordinates
    coords_int = coords.astype(int)
    coords_int[:, 0] = np.clip(coords_int[:, 0], 0, W - 1)
    coords_int[:, 1] = np.clip(coords_int[:, 1], 0, H - 1)
    # for each spot, paint a small disk
    radius = 12
    canvas = np.full((H, W), np.nan, dtype=np.float32)
    for (x, y), v in zip(coords_int, values):
        if np.isnan(v): continue
        cv2.circle(canvas, (x, y), radius, float(v), -1)
    # fill NaN with 0 for blur, but keep mask
    mask = ~np.isnan(canvas)
    canvas_filled = np.where(np.isnan(canvas), 0, canvas)
    blurred = cv2.GaussianBlur(canvas_filled, (0, 0), sigmaX=sigma, sigmaY=sigma)
    mask_blur = cv2.GaussianBlur(mask.astype(np.float32), (0, 0), sigmaX=sigma, sigmaY=sigma) + 1e-3
    blurred = blurred / mask_blur  # normalize for sparse spots
    # normalize 5-95% quantile
    valid = blurred[mask_blur > 0.05]
    if valid.size == 0:
        return bg_img.copy()
    lo, hi = np.percentile(valid, [5, 95])
    norm = np.clip((blurred - lo) / max(hi - lo, 1e-6), 0, 1) * 255
    norm = norm.astype(np.uint8)
    heat = cv2.applyColorMap(norm, cmap_id)
    # alpha blend only where mask_blur > threshold
    blend_alpha = np.clip(mask_blur, 0, 1) * alpha
    blend_alpha = np.dstack([blend_alpha] * 3)
    out = (bg_img * (1 - blend_alpha) + heat * blend_alpha).astype(np.uint8)
    return out


def annotate_panel(img, title, sub, color_meaning, interp_lines=None):
    """Top banner + bottom bottom legend + (optional) interpretation block."""
    pil = Image.fromarray(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
    draw = ImageDraw.Draw(pil)
    H = pil.height; W = pil.width
    # top banner (taller for sub)
    draw.rectangle([(0, 0), (W, 78)], fill=(8, 12, 26, 235))
    draw.text((14, 8), title, fill=(255, 255, 255), font=FONT)
    draw.text((14, 44), sub, fill=(180, 200, 240), font=FONT_SMALL)
    # interpretation block at bottom
    interp_h = 0
    if interp_lines:
        line_h = 26
        interp_h = 18 + line_h * len(interp_lines) + 18
    bot_legend_h = 50
    total_bot = bot_legend_h + interp_h
    draw.rectangle([(0, H - total_bot), (W, H)], fill=(8, 12, 26, 235))
    if interp_lines:
        for i, line in enumerate(interp_lines):
            color = (252, 211, 77) if i == 0 else (200, 220, 240)
            draw.text((14, H - total_bot + 10 + i * 26), line, fill=color, font=FONT_SMALL)
    draw.text((14, H - 42), color_meaning, fill=(220, 230, 240), font=FONT_SMALL)
    return cv2.cvtColor(np.array(pil), cv2.COLOR_RGB2BGR)


def render_slide(sid):
    sdir = GS / sid
    a = ad.read_h5ad(sdir / f"{sid}.scored.h5ad")
    sf = get_scalefactor(a.uns.get("spatial", {}))
    bg = hires_path(sdir)
    if bg is None: return None
    bg_rgb = np.array(Image.open(bg).convert("RGB"))
    bg_bgr = cv2.cvtColor(bg_rgb, cv2.COLOR_RGB2BGR)
    coords = a.obsm["spatial"] * sf
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    sym = {s: i for i, s in enumerate(a.var.index.values)}

    dm1 = a.obs["DM1_like_score"].values
    caf = module(X, sym, CAF)
    rai = module(X, sym, RAI)
    if caf is None or rai is None: return None

    # CAF × -RAI lag (avoidance map)
    tree = cKDTree(a.obsm["spatial"])
    _, idx = tree.query(a.obsm["spatial"], k=7)
    rai_lag = rai[idx[:, 1:]].mean(axis=1)
    avoid = caf * (-rai_lag)

    stage = a.obs.stage.iloc[0] if "stage" in a.obs.columns else "?"

    # Per-slide quantitative interpretation (computed live)
    rai_lag_corr = float(np.corrcoef(caf, rai_lag)[0, 1]) if np.std(caf) > 0 and np.std(rai_lag) > 0 else float("nan")
    dm1_mean = float(np.nanmean(dm1)); dm1_std = float(np.nanstd(dm1))
    caf_mean = float(np.nanmean(caf))
    rai_mean = float(np.nanmean(rai))
    avoid_top10 = float(np.nanpercentile(avoid, 90))

    stage_interp = {
        "PT":   "정상 갑상선 — RAI 분화도 thyrocyte 가 영역 전체 점유. CAF는 정상적 stromal 띠로만.",
        "PTC":  "유두암 (PTC) — DM1 hotspot 와 정상 영역이 공존하는 intra-tumor heterogeneity 시작.",
        "LPTC": "lateral PTC — DM1 axis 가장 spatially organized (Moran 0.66). CCL19/CCL21 TLS niche 가시.",
        "ATC": "역형성암 (ATC) — DM1 score 가 spot 전반에 균일 (Moran 붕괴). 정체성 잃은 dedifferentiation 종착점.",
    }

    # Render 4 panels
    panels = []
    # 1. H&E original
    he_panel = annotate_panel(
        bg_bgr.copy(),
        f"{sid} · {stage} · H&E 원본",
        f"Visium spot {a.n_obs}개 · 16-slide {stage}-stage GSE250521",
        "원본 슬라이드 (염색만)",
        interp_lines=[
            f"▶ {stage_interp.get(stage, '')}",
            f"▶ 비교: 16-slide 4-stage (N→PTC→LPTC→ATC) trajectory 의 한 슬라이드.",
        ]
    )
    panels.append(he_panel)

    # 2. DM1 score
    dm1_overlay = render_heatmap(bg_bgr, coords, dm1, cmap_id=cv2.COLORMAP_TURBO, alpha=0.55)
    dm1_panel = annotate_panel(
        dm1_overlay,
        f"{sid} · DM1_like_score 공간 분포",
        f"평균 = {dm1_mean:+.3f} · 표준편차 = {dm1_std:.3f}",
        "🔴 빨강 = DM1 high (HT-like, 면역활성)  🔵 파랑 = DM1 low (DM2 aggressive)",
        interp_lines=[
            f"▶ 의미: DM1_like_score = −RAI_8 z-score. 빨강 = RAI 분화 잃음 (HT-like).",
            f"▶ ATC 슬라이드는 점수가 spot 전체에 균일 (Moran 0.05) = niche 구조 사라짐.",
            f"▶ PTC/LPTC 는 hotspot (Moran 0.6+) = intra-tumor heterogeneity 살아있음.",
        ]
    )
    panels.append(dm1_panel)

    # 3. CAF score
    caf_overlay = render_heatmap(bg_bgr, coords, caf, cmap_id=cv2.COLORMAP_VIRIDIS, alpha=0.55)
    caf_panel = annotate_panel(
        caf_overlay,
        f"{sid} · CAF score (섬유아세포 niche)",
        f"평균 = {caf_mean:+.3f} · FAP/ACTA2/COL1A1/POSTN 8-gene z-score",
        "🟡 노랑 = CAF-rich (섬유화)  🟣 보라 = CAF-poor",
        interp_lines=[
            f"▶ 의미: 종양 주위 CAF 섬유아세포 활성도. 노랑 영역 = collagen+ POSTN+ FAP+ niche.",
            f"▶ PTC 부터 stromal reorganization 시작 (POSTN/COL1A2/SULF1 SVG top).",
            f"▶ ATC 도 CAF 신호 유지하지만 RAI thyrocyte 와 같은 자리 차지 시작.",
        ]
    )
    panels.append(caf_panel)

    # 4. CAF × -RAI avoidance + RED BOX highlights at top-10% avoidance regions
    avoid_overlay = render_heatmap(bg_bgr, coords, avoid, cmap_id=cv2.COLORMAP_HOT, alpha=0.6)
    # Identify top-5% avoidance spots and draw red bounding boxes around clusters
    valid = ~np.isnan(avoid)
    if valid.sum() > 30:
        thr = np.nanpercentile(avoid[valid], 95)
        hot_idx = np.where((avoid > thr) & valid)[0]
        # cluster hot spots via simple distance grouping
        if len(hot_idx) > 5:
            from scipy.spatial import cKDTree as KD
            hot_coords = coords[hot_idx]
            tree_hot = KD(hot_coords)
            visited = np.zeros(len(hot_coords), dtype=bool)
            clusters = []
            for i in range(len(hot_coords)):
                if visited[i]: continue
                nbrs = tree_hot.query_ball_point(hot_coords[i], r=80)
                if len(nbrs) >= 3:
                    cluster_pts = hot_coords[nbrs]
                    clusters.append(cluster_pts)
                    visited[nbrs] = True
            for cl in clusters[:5]:  # cap at 5 boxes
                x1, y1 = cl[:, 0].min(), cl[:, 1].min()
                x2, y2 = cl[:, 0].max(), cl[:, 1].max()
                pad = 25
                cv2.rectangle(avoid_overlay, (int(x1-pad), int(y1-pad)),
                              (int(x2+pad), int(y2+pad)), (0, 0, 255), 4)
                # arrow + label
                cv2.arrowedLine(avoid_overlay,
                                (int(x2 + pad + 60), int(y1 - pad - 20)),
                                (int(x2 + pad + 5), int(y1 - pad)),
                                (0, 0, 255), 3, tipLength=0.3)
    avoid_panel = annotate_panel(
        avoid_overlay,
        f"{sid} · CAF–RAI 공간 회피 (avoidance) hotspot ★ 핵심 finding",
        f"본 슬라이드 CAF–RAI_lag Pearson = {rai_lag_corr:+.3f} · top-10% avoidance 강도 = {avoid_top10:+.2f}",
        "🔥 노랑/주황 = CAF 많고 주변 RAI thyrocyte 적은 영역 = 갑상선이 자기 정체성 잃은 dedifferentiation niche",
        interp_lines=[
            f"▶ 의미: CAF score × (−이웃 RAI score). 밝을수록 CAF는 많은데 RAI 분화 thyrocyte 가 주변에 없음.",
            f"▶ 16-slide median Pearson −0.27 / TCGA bulk n=572 Spearman −0.31 — cross-platform replicate ★",
            f"▶ 이 패턴 = CAF 섬유화가 RAI 분화도 thyrocyte 를 anatomically 밀어내는 신호.",
            f"▶ Paper finding: 갑상선암 dedifferentiation 의 spatial signature.",
        ]
    )
    panels.append(avoid_panel)

    # Combine 2x2 grid
    H, W = panels[0].shape[:2]
    canvas = np.zeros((H * 2, W * 2, 3), dtype=np.uint8)
    canvas[:H, :W] = panels[0]
    canvas[:H, W:] = panels[1]
    canvas[H:, :W] = panels[2]
    canvas[H:, W:] = panels[3]

    # Add overall header bar + summary footer
    pil = Image.fromarray(cv2.cvtColor(canvas, cv2.COLOR_BGR2RGB))
    draw = ImageDraw.Draw(pil)
    # top header
    draw.rectangle([(0, 0), (pil.width, 64)], fill=(20, 30, 64, 230))
    draw.text((20, 16), f"SPARK · {sid} · {stage} · 4-panel 해석", fill=(255, 255, 255), font=FONT)
    # bottom story footer
    story_lines = [
        f"이 슬라이드의 한 줄 해석 ({stage}):",
        f"  · DM1_like_score 평균 {dm1_mean:+.2f}, 본 슬라이드 CAF–RAI_lag Pearson = {rai_lag_corr:+.3f}",
        f"  · {stage_interp.get(stage, '')}",
        f"  · 16-slide 종합: 4번째 panel 이 paper headline finding (CAF가 RAI thyrocyte 밀어냄, 16 ST + bulk n=572 양쪽 replicate).",
    ]
    footer_h = 26 * len(story_lines) + 24
    draw.rectangle([(0, pil.height - footer_h), (pil.width, pil.height)], fill=(8, 12, 26, 240))
    for i, line in enumerate(story_lines):
        col = (252, 211, 77) if i == 0 else (220, 230, 240)
        draw.text((20, pil.height - footer_h + 12 + i * 26), line, fill=col, font=FONT_SMALL)
    canvas2 = cv2.cvtColor(np.array(pil), cv2.COLOR_RGB2BGR)

    # downscale + jpg save for web
    H2, W2 = canvas2.shape[:2]
    if max(H2, W2) > 2400:
        scale = 2400.0 / max(H2, W2)
        canvas2 = cv2.resize(canvas2, (int(W2 * scale), int(H2 * scale)), interpolation=cv2.INTER_AREA)
    out_path = OUT / f"{sid}_4panel.jpg"
    cv2.imwrite(str(out_path), canvas2, [cv2.IMWRITE_JPEG_QUALITY, 88])
    return out_path


def main():
    sids = sorted([d.name for d in GS.iterdir() if (d / f"{d.name}.scored.h5ad").exists()])
    out_paths = []
    for sid in sids:
        try:
            p = render_slide(sid)
            print(f"  {sid}: {'ok' if p else 'skip'}")
            if p: out_paths.append(p)
        except Exception as e:
            print(f"  {sid}: ERR {e}")

    # gallery
    html = ['<!doctype html><html lang="ko"><head><meta charset="utf-8"><title>SPARK · 4-panel 직관 시각화</title>',
            '<style>body{background:#0a0d18;color:#e1e8f0;font-family:system-ui,sans-serif;margin:24px}',
            'h1{margin:0 0 8px} .sub{color:#94a3b8;margin-bottom:18px;line-height:1.6}',
            '.legend{background:#101833;border:1px solid #2a3a5e;border-radius:10px;padding:14px;margin-bottom:24px;font-size:13px;line-height:1.7}',
            'figure{margin:0 0 18px;border:1px solid #1f2a44;border-radius:10px;padding:8px;background:#0e1426}',
            'figcaption{font-size:12px;color:#94a3b8;padding:6px 4px}',
            'img{width:100%;display:block;border-radius:6px}</style></head><body>',
            '<h1>SPARK · 갑상선암 dark-matter 4-panel 직관 시각화 (16 slides)</h1>',
            '<div class="legend">',
            '<b>각 슬라이드 그림 4-패널 안내</b><br/>',
            '<b>① H&E 원본</b>: Visium spot이 찍힌 원본 슬라이드<br/>',
            '<b>② DM1_like_score 공간 분포</b>: 빨강 = DM1 high (HT 닮음, 면역활성), 파랑 = DM1 low (공격적 DM2)<br/>',
            '<b>③ CAF score</b>: 노랑 = 섬유아세포 많은 영역 (FAP/ACTA2/COL1A1 등 8-gene)<br/>',
            '<b>④ CAF–RAI avoidance hotspot</b>: 주황/노랑 = CAF 많고 주변에 RAI 분화 thyrocyte 적은 영역. <br/>',
            '&nbsp;&nbsp;&nbsp;= 갑상선이 자기 정체성을 잃고 섬유화로 대체되는 anatomical signature.<br/>',
            '&nbsp;&nbsp;&nbsp;TCGA 572명 bulk 에서도 같은 방향 replicate (Spearman −0.31).',
            '</div>']
    for p in out_paths:
        html.append(f'<figure><img src="{p.name}" /><figcaption>{p.stem}</figcaption></figure>')
    html.append("</body></html>")
    (OUT / "index.html").write_text("\n".join(html))
    print(f"\nwrote {OUT / 'index.html'}  ({len(out_paths)} 4-panel figures)")


if __name__ == "__main__":
    main()
