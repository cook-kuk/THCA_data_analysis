#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from textwrap import dedent

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "project" / "results" / "hla_two_paper_synthesis_2026_05_09"
OUT = BASE / "execution_board"
FIG_OUT = OUT / "figures"
HUB = ROOT / "project" / "papers_hub_2026_05_04"


PAPER2_PANELS = [
    (
        "A. HT-overlap bulk",
        "GSE286332: HLA-II d=3.65; TLS d=3.01",
        BASE / "extra_external_validation/figures/F01_gse286332_ht_overlap_hla_tls_effects.png",
    ),
    (
        "B. Exact small-n robustness",
        "HLA-II exact p=0.000062; AP/TLS AUC=0.96",
        BASE / "level_up_wave2/figures/F03_gse286332_exact_permutation_bootstrap.png",
    ),
    (
        "C. Korean bulk validation",
        "GSE213647 n=632: HLA-II/AP d=1.39, FDR=4.03e-40",
        BASE / "gse213647_korean_bulk_validation/figures/F02_gse213647_ptc_effect_sizes.png",
    ),
    (
        "D. Thyroid Visium progression",
        "GSE250521: 16 slides, 55,873 spots; CD74/MIF rho=0.76",
        BASE / "gse250521_thyroid_visium_hla_validation/figures/F02_gse250521_stage_trend_rho.png",
    ),
]

PAPER4_PANELS = [
    (
        "A. Allele triangulation",
        "C*01:02 specificity-led; DPB1*05:01 AITD-broad",
        BASE / "level_up_wave2/figures/F01_paper4_allele_triangulation_map.png",
    ),
    (
        "B. Endpoint robustness",
        "Strict GD separated from broad hyperthyroid/HT endpoints",
        BASE / "level_up_wave2/figures/F02_paper4_endpoint_robustness_heatmap.png",
    ),
    (
        "C. AITD spatial tissue biology",
        "GSE248205: CD74/MIF d=3.62; AP/TLS d=2.26",
        BASE / "gse248205_aitd_spatial_validation/figures/F03_gse248205_spot_burden_heatmap.png",
    ),
    (
        "D. Older-array corroboration",
        "GSE29315 HT: HLA-II/AP d=5.45, FDR=0.0023",
        BASE / "gse29315_aitd_array_validation/figures/F02_gse29315_ht_effect_sizes.png",
    ),
]


def font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def fit_image(path: Path, width: int, height: int) -> Image.Image:
    img = Image.open(path).convert("RGB")
    img.thumbnail((width, height), Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", (width, height), "white")
    x = (width - img.width) // 2
    y = (height - img.height) // 2
    canvas.paste(img, (x, y))
    return canvas


def draw_panel(
    canvas: Image.Image,
    xy: tuple[int, int],
    panel_size: tuple[int, int],
    title: str,
    subtitle: str,
    image_path: Path,
) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y = xy
    w, h = panel_size
    header_h = 78
    draw.rounded_rectangle([x, y, x + w, y + h], radius=18, fill="#ffffff", outline="#cad4dc", width=2)
    draw.rounded_rectangle([x, y, x + w, y + header_h], radius=18, fill="#0d1823")
    draw.rectangle([x, y + header_h - 18, x + w, y + header_h], fill="#0d1823")
    draw.text((x + 18, y + 13), title, font=font(28, True), fill="#ffd28a")
    draw.text((x + 18, y + 47), subtitle, font=font(18), fill="#d7e4ea")
    img = fit_image(image_path, w - 26, h - header_h - 26)
    canvas.paste(img, (x + 13, y + header_h + 13))


def make_montage(panels: list[tuple[str, str, Path]], title: str, subtitle: str, out_path: Path) -> None:
    missing = [str(path) for _, _, path in panels if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing montage inputs:\n" + "\n".join(missing))

    panel_w, panel_h = 760, 560
    margin, gap = 36, 24
    title_h = 132
    width = margin * 2 + panel_w * 2 + gap
    height = margin * 2 + title_h + panel_h * 2 + gap
    canvas = Image.new("RGB", (width, height), "#071016")
    draw = ImageDraw.Draw(canvas)
    draw.text((margin, 28), title, font=font(42, True), fill="#ffffff")
    draw.text((margin, 82), subtitle, font=font(22), fill="#d7e4ea")

    positions = [
        (margin, margin + title_h),
        (margin + panel_w + gap, margin + title_h),
        (margin, margin + title_h + panel_h + gap),
        (margin + panel_w + gap, margin + title_h + panel_h + gap),
    ]
    for pos, panel in zip(positions, panels):
        draw_panel(canvas, pos, (panel_w, panel_h), *panel)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_path, optimize=True, quality=95)


def write_execution_md(md_path: Path) -> None:
    md_path.write_text(
        dedent(
            """\
            # HLA Execution Board — 2026-05-09

            ## Verdict

            HLA stays as two separated papers:

            - Paper 2: HT-overlap PTC antigen-presentation / TLS / IFN-gamma expression ecology. This is the stronger near-term NComm challenge.
            - Paper 4: Pan-Asian Graves/AITD HLA genetics. This is strong, but NComm-level upgrade depends on Korean adult GD germline HLA validation or equivalent source-level genetics.

            ## Do Not Cross The Streams

            - Paper 2 must not claim HLA allele-based thyroid cancer risk.
            - Paper 4 must not use AITD tissue-expression data as genotype replication.
            - DRB1*04:05 is a prospective HT-PTC candidate only.
            - DPB1*05:01 is not an HT-PTC allele finding.

            ## Immediate Execution

            1. Freeze Paper 2 as expression ecology, not allele association.
            2. Use the Paper 2 montage as the advisor-facing evidence wall.
            3. Keep Paper 4 as a separate GD-HLA genetics track.
            4. Gate Paper 4 NComm claims on Korean adult GD HLA NGS / imputation access.
            """
        ),
        encoding="utf-8",
    )


def write_html(html_path: Path) -> None:
    paper2_img = "../results/hla_two_paper_synthesis_2026_05_09/execution_board/figures/HLA_PAPER2_CORE_MONTAGE.png"
    paper4_img = "../results/hla_two_paper_synthesis_2026_05_09/execution_board/figures/HLA_PAPER4_CORE_MONTAGE.png"
    html_path.write_text(
        dedent(
            f"""\
            <!DOCTYPE html>
            <html lang="ko">
            <head>
            <meta charset="UTF-8" />
            <meta name="viewport" content="width=device-width, initial-scale=1.0" />
            <title>HLA Execution Board | 2026-05-09</title>
            <style>
            :root{{--bg:#071016;--panel:#0d1823;--ink:#edf6fb;--muted:#9fb1bd;--line:#273a4a;--gold:#ffd28a;--green:#93e2aa;--red:#ff8a6b;--blue:#8bbcff;--teal:#65dcc8}}
            *{{box-sizing:border-box}} html,body{{margin:0;background:var(--bg);color:var(--ink);font-family:Inter,"Noto Sans KR",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;line-height:1.56}}
            a{{color:var(--teal);text-decoration:none}} a:hover{{text-decoration:underline}}
            .hero{{min-height:68vh;padding:64px 30px 36px;background:linear-gradient(135deg,#071016 0%,#102234 62%,#05090d 100%);border-bottom:1px solid var(--line)}}
            .inner{{max-width:1280px;margin:0 auto}} .kicker{{font-family:"JetBrains Mono","SF Mono",monospace;color:var(--gold);letter-spacing:.18em;text-transform:uppercase;font-size:11px;font-weight:900;margin:0 0 16px}}
            h1{{font-family:Georgia,"Times New Roman","Noto Serif KR",serif;font-size:68px;line-height:.95;margin:0 0 18px;color:#fff;letter-spacing:0}}
            .lead{{font-family:Georgia,"Times New Roman","Noto Serif KR",serif;font-size:23px;color:#d8e5ec;max-width:1080px;margin:0}} .lead b{{color:var(--gold)}}
            .stats{{display:grid;grid-template-columns:repeat(6,1fr);gap:10px;margin-top:28px}} .stat{{border:1px solid rgba(255,210,138,.24);background:rgba(255,255,255,.04);border-radius:8px;min-height:116px;padding:14px}}
            .stat b{{display:block;font-family:Georgia,serif;font-size:34px;line-height:1;color:var(--gold)}} .stat span{{display:block;font-weight:900;margin-top:5px}} .stat small{{display:block;color:var(--muted);font-family:"JetBrains Mono","SF Mono",monospace;font-size:10px;margin-top:4px}}
            .actions{{display:flex;flex-wrap:wrap;gap:10px;margin-top:24px}} .actions a{{border:1px solid var(--line);border-radius:6px;background:#0f1c28;color:var(--ink);padding:10px 13px;font-family:"JetBrains Mono","SF Mono",monospace;font-size:12px;font-weight:900}} .actions a.primary{{border-color:var(--gold);color:var(--gold)}}
            main{{max-width:1280px;margin:0 auto;padding:30px}} section{{border-bottom:1px solid var(--line);padding:32px 0}} h2{{font-family:Georgia,"Times New Roman","Noto Serif KR",serif;font-size:42px;margin:0 0 8px;color:#fff}} h3{{margin:0 0 8px;color:#fff;font-size:19px}}
            .sub{{font-family:"JetBrains Mono","SF Mono",monospace;color:var(--muted);font-size:11px;letter-spacing:.12em;text-transform:uppercase;margin:0 0 18px}}
            .grid{{display:grid;gap:14px}} .g2{{grid-template-columns:repeat(2,minmax(0,1fr))}} .g3{{grid-template-columns:repeat(3,minmax(0,1fr))}}
            .box{{border:1px solid var(--line);background:rgba(13,24,35,.94);border-radius:8px;padding:18px}} .box p{{margin:0;color:#d1dde4}} .box ul{{margin:8px 0 0;padding-left:18px;color:#d1dde4}} .box li{{margin:6px 0}}
            .green{{border-top:5px solid var(--green)}} .gold{{border-top:5px solid var(--gold)}} .red{{border-top:5px solid var(--red)}} .blue{{border-top:5px solid var(--blue)}}
            table{{width:100%;border-collapse:collapse;margin:12px 0 18px;font-size:14px}} th,td{{border:1px solid var(--line);padding:10px;text-align:left;vertical-align:top}} th{{background:#142333;color:var(--gold);font-family:"JetBrains Mono","SF Mono",monospace;font-size:11px;text-transform:uppercase;letter-spacing:.06em}} td{{background:rgba(255,255,255,.025);color:#d1dde4}}
            .good{{color:var(--green);font-weight:900}} .warn{{color:var(--gold);font-weight:900}} .bad{{color:var(--red);font-weight:900}}
            figure{{margin:18px 0;background:#fff;border:1px solid #d6e0e6;border-radius:10px;overflow:hidden}} figure img{{width:100%;display:block}} figcaption{{background:#0f1c28;color:#d2e2ea;border-top:1px solid var(--line);padding:12px 14px;font-size:13px}} figcaption b{{color:var(--gold)}}
            .footer{{font-family:"JetBrains Mono","SF Mono",monospace;color:var(--muted);font-size:11px;margin-top:18px}}
            code{{font-family:"JetBrains Mono","SF Mono",monospace;color:#eaf7ff}}
            @media(max-width:900px){{h1{{font-size:44px}}.hero{{padding:46px 18px 30px}}main{{padding:22px 18px 70px}}.stats{{grid-template-columns:repeat(2,1fr)}}.g2,.g3{{grid-template-columns:1fr}}}}
            </style>
            </head>
            <body>
            <header class="hero">
              <div class="inner">
                <p class="kicker">HLA execution board · 2026-05-09</p>
                <h1>HLA는 두 논문으로 계속 간다. 섞지 말고, 둘 다 키운다.</h1>
                <p class="lead">Paper 2는 <b>HT-overlap PTC HLA-II/AP/TLS expression ecology</b>로 NComm 도전권이 있고, Paper 4는 <b>Pan-Asian Graves/AITD HLA genetics</b>로 별도 트랙이다. 핵심은 allele claim과 expression claim을 절대 섞지 않는 것이다.</p>
                <div class="stats">
                  <div class="stat"><b>2</b><span>separate papers</span><small>PTC ecology / GD genetics</small></div>
                  <div class="stat"><b>0.96</b><span>AP/TLS AUC</span><small>GSE286332 exact robustness</small></div>
                  <div class="stat"><b>632</b><span>Korean bulk</span><small>GSE213647</small></div>
                  <div class="stat"><b>55,873</b><span>Visium spots</span><small>GSE250521</small></div>
                  <div class="stat"><b>5.45</b><span>AITD HLA-II/AP d</span><small>GSE29315 HT array</small></div>
                  <div class="stat"><b>K2/GD</b><span>upgrade gates</span><small>FFPE / adult GD HLA</small></div>
                </div>
                <div class="actions">
                  <a class="primary" href="#paper2">Paper 2 montage</a>
                  <a href="#paper4">Paper 4 montage</a>
                  <a href="hla_ncomm_upgrade_2026_05_09.html">Full HLA NComm dossier</a>
                  <a href="../results/hla_two_paper_synthesis_2026_05_09/execution_board/HLA_EXECUTION_BOARD.md">MD execution board</a>
                  <a href="index.html">Portfolio hub</a>
                </div>
              </div>
            </header>
            <main>
              <section id="verdict">
                <h2>Decision</h2>
                <p class="sub">what changed after the HLA level-up</p>
                <div class="grid g3">
                  <div class="box green"><h3>Paper 2 is the near-term push</h3><p>Use the expression ecology spine: HT-overlap, HLA-II/AP, IFN-gamma, B/TLS, scRNA/spatial and negative stress-tests. This is the strongest current HLA story.</p></div>
                  <div class="box gold"><h3>Paper 4 stays separate</h3><p>GD/AITD HLA genetics has a real Pan-Asian architecture, but NComm-level claim needs Korean adult GD HLA validation or equivalent genetics access.</p></div>
                  <div class="box red"><h3>Boundary is non-negotiable</h3><p>Do not claim PTC HLA allele risk. Do not call tissue expression genotype replication. Do not call DRB1*04:05 validated.</p></div>
                </div>
              </section>
              <section id="paper2">
                <h2>Paper 2 Evidence Wall</h2>
                <p class="sub">HT-overlap PTC antigen-presentation/TLS ecology</p>
                <figure>
                  <img src="{paper2_img}" alt="Paper 2 HLA core evidence montage" />
                  <figcaption><b>Readout.</b> This makes Paper 2 an expression-ecology paper: HT-overlap bulk signal, exact small-n robustness, Korean bulk generalization and thyroid Visium progression all point toward HLA/AP-TLS biology.</figcaption>
                </figure>
                <table>
                  <thead><tr><th>Dataset</th><th>Key result</th><th>Use</th></tr></thead>
                  <tbody>
                    <tr><td>GSE286332</td><td>HLA-II d=3.65, exact p=0.000062; AP/TLS composite d=3.08, AUC=0.96.</td><td class="good">Main HT-overlap validation.</td></tr>
                    <tr><td>GSE213647</td><td>Korean bulk n=632; PTC tumor vs normal HLA-II/AP d=1.39, FDR=4.03e-40.</td><td class="good">Korean generalization.</td></tr>
                    <tr><td>GSE250521</td><td>16 Visium slides, 55,873 spots; CD74/MIF rho=0.76 q=0.0027; HLA-II/AP rho=0.61 q=0.018.</td><td class="good">Spatial progression.</td></tr>
                    <tr><td>GSE6004</td><td>Invasion vs normal negative: HLA-II/AP d=-0.42, AP/TLS d=-0.59, FDR=0.829.</td><td class="warn">Specificity stress-test.</td></tr>
                  </tbody>
                </table>
              </section>
              <section id="paper4">
                <h2>Paper 4 Evidence Wall</h2>
                <p class="sub">Pan-Asian Graves/AITD HLA genetics plus tissue mechanism</p>
                <figure>
                  <img src="{paper4_img}" alt="Paper 4 HLA core evidence montage" />
                  <figcaption><b>Readout.</b> This keeps Paper 4 as a genetics paper: allele triangulation and endpoint robustness define the HLA architecture; AITD spatial/array data support tissue antigen-presentation biology without being genotype replication.</figcaption>
                </figure>
                <table>
                  <thead><tr><th>Axis</th><th>Key result</th><th>Use</th></tr></thead>
                  <tbody>
                    <tr><td>DPB1*05:01</td><td>Pan-Asian anchor, but AITD-broad rather than GD-specific.</td><td class="good">Main class-II anchor with caveat.</td></tr>
                    <tr><td>C*01:02</td><td>Highest GD-specificity signal in the current triangulation.</td><td class="good">Specificity-led GD candidate.</td></tr>
                    <tr><td>GSE248205</td><td>AITD spatial: CD74/MIF d=3.62; AP/TLS d=2.26.</td><td class="warn">Tissue mechanism only.</td></tr>
                    <tr><td>GSE29315</td><td>HT array: HLA-II/AP d=5.45, exact p=0.0007, FDR=0.0023.</td><td class="warn">Independent mechanism corroboration.</td></tr>
                  </tbody>
                </table>
              </section>
              <section id="queue">
                <h2>Next Execution Queue</h2>
                <p class="sub">paper-blocking only</p>
                <div class="grid g2">
                  <div class="box green"><h3>Now</h3><ul><li>Freeze Paper 2 around the expression-ecology claim.</li><li>Use the Paper 2 montage for advisor review.</li><li>Keep allele candidates in prospective-validation language only.</li></ul></div>
                  <div class="box gold"><h3>Gate</h3><ul><li>K2/FFPE H&E or Korean PTC+HT data upgrades Paper 2.</li><li>Korean adult GD HLA-NGS or imputed genetics upgrades Paper 4.</li><li>Do not spend marathon time on non-blocking HLA allele joins.</li></ul></div>
                </div>
                <p class="footer">Generated by <code>scripts/hla_deepdive_2026_05_08/build_hla_execution_board_2026_05_09.py</code>. Figures and MD live under <code>project/results/hla_two_paper_synthesis_2026_05_09/execution_board/</code>.</p>
              </section>
            </main>
            </body>
            </html>
            """
        ),
        encoding="utf-8",
    )


def main() -> None:
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    make_montage(
        PAPER2_PANELS,
        "Paper 2 core: HT-overlap PTC HLA/AP-TLS ecology",
        "Expression-module evidence only; allele association remains prospective.",
        FIG_OUT / "HLA_PAPER2_CORE_MONTAGE.png",
    )
    make_montage(
        PAPER4_PANELS,
        "Paper 4 core: Pan-Asian GD/AITD HLA architecture",
        "Allele/genetics track only; tissue expression supports mechanism, not genotype replication.",
        FIG_OUT / "HLA_PAPER4_CORE_MONTAGE.png",
    )
    write_execution_md(OUT / "HLA_EXECUTION_BOARD.md")
    write_html(HUB / "hla_execution_board_2026_05_09.html")


if __name__ == "__main__":
    main()
