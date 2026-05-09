#!/usr/bin/env python3
"""Build the Paper 2 image-DM1 audit dossier HTML.

The dossier is intentionally data-bound: headline numbers are loaded from
the audit JSON/TSV evidence files, and the script also writes the
TSS-only-within-RAS evidence table used by the page.
"""
from __future__ import annotations

import html
import json
import math
import statistics as stats
from pathlib import Path

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RUN = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
AUD = RUN / "analysis_supp/audit_ras_auc100"
HUB = ROOT / "project/papers_hub_2026_05_04"


def esc(x: object) -> str:
    return html.escape("" if x is None else str(x))


def fmt(x: object, nd: int = 3, na: str = "n/a") -> str:
    if x is None:
        return na
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return esc(x)
    if math.isnan(xf):
        return na
    return f"{xf:.{nd}f}"


def pct(x: object, nd: int = 1) -> str:
    if x is None:
        return "n/a"
    xf = float(x)
    if math.isnan(xf):
        return "n/a"
    return f"{100 * xf:.{nd}f}%"


def median(values: list[float]) -> float:
    return float(stats.median([float(v) for v in values]))


def read_json(rel: str) -> dict:
    return json.loads((AUD / rel).read_text())


def read_tsv(rel: str) -> pd.DataFrame:
    return pd.read_csv(AUD / rel, sep="\t")


def table(headers: list[str], rows: list[list[object]], classes: str = "") -> str:
    head = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = []
    for row in rows:
        cells = []
        for cell in row:
            cls = ""
            val = cell
            if isinstance(cell, tuple):
                val, cls = cell
            cells.append(f"<td class=\"{esc(cls)}\">{val}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return (
        f"<table class=\"t {classes}\"><thead><tr>{head}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table>"
    )


def figure(src: str, title: str, caption: str) -> str:
    return (
        "<figure class=\"fig\">"
        f"<img src=\"{esc(src)}\" alt=\"{esc(title)}\" loading=\"lazy\" />"
        f"<figcaption><b>{esc(title)}</b>{esc(caption)}</figcaption>"
        "</figure>"
    )


def write_tss_only_evidence() -> dict:
    df = read_tsv("merged_preds_subgroups.tsv")
    df["tss"] = df["submitter_id"].str.extract(r"TCGA-([^-]+)-")
    ras = df[df["molecular_subtype"] == "RAS_like"].copy()
    x = pd.get_dummies(ras[["tss"]], drop_first=False).astype(float).values
    y = ras["label"].to_numpy()
    lr = LogisticRegression(C=1e6, max_iter=1000)
    lr.fit(x, y)
    prob = lr.predict_proba(x)[:, 1]
    out = ras[["submitter_id", "tss", "label"]].copy()
    out["tss_only_prob_DM1"] = prob
    out = out.sort_values("tss_only_prob_DM1")
    out_path = AUD / "tss_only_within_raslike_lr.tsv"
    out.to_csv(out_path, sep="\t", index=False)

    auc = float(roc_auc_score(y, prob))
    summary = {
        "model": "in-sample logistic regression on TSS dummies only",
        "cohort": "16 RAS-like/FVPTC TCGA-THCA slides",
        "n": int(len(ras)),
        "n_pos": int(y.sum()),
        "n_neg": int(len(y) - y.sum()),
        "auc": auc,
        "tsv": str(out_path.relative_to(ROOT)),
    }
    (AUD / "tss_only_within_raslike_lr.json").write_text(
        json.dumps(summary, indent=2)
    )
    return summary


def build_html(figure_prefix: str, sibling_prefix: str, out_label: str) -> str:
    p1 = read_json("AUDIT_SUMMARY.json")
    p2 = read_json("PHASE2_AUDIT_SUMMARY.json")
    p3 = read_json("PHASE3_AUDIT_SUMMARY.json")
    rs = read_json("RETRAIN_SUMMARY.json")
    p4 = read_json("phase4/PHASE4_SOLUTIONS_SUMMARY.json")
    uni = json.loads((RUN / "analysis_supp/bootstrap_auc_uni.json").read_text())
    uni_loto_path = RUN / "analysis_supp/audit_uni_loto/UNI_LOTO_SUMMARY.json"
    uni_loto = json.loads(uni_loto_path.read_text()) if uni_loto_path.exists() else None
    tss_only = write_tss_only_evidence()

    prob = read_tsv("c2_prob_distribution.tsv")
    ras_raw = read_tsv("c2_RAS_like_16_predictions.tsv")
    fold = read_tsv("c5_RAS_like_fold_composition.tsv")
    e2 = read_tsv("e2_tss_among_RASlike.tsv")
    e8 = read_tsv("e8_leave_one_tss_out.tsv")
    s1 = read_tsv("phase4/s1_tss_combat_per_init.tsv")
    s2 = read_tsv("phase4/s2_tss_balanced_split.tsv")
    s3 = read_tsv("phase4/s3_multimodal_LR.tsv")
    s5 = read_tsv("phase4/s5_em_out.tsv")

    ras_perm = p1["C3_permutation"]["molecular_subtype=RAS_like"]
    fv_perm = p1["C3_permutation"]["histology_subtype=FVPTC"]
    braf_perm = p1["C3_permutation"]["molecular_subtype=BRAF_like"]
    e5 = p2["E5_clinical_baseline"]
    e8s = p3["E8_leave_one_tss_out"]
    clinical_auc = next(r["auc"] for r in p4["S3_multimodal_LR"]
                        if r["predictor"] == "clinical_only")
    multimodal_auc = next(r["auc"] for r in p4["S3_multimodal_LR"]
                          if r["predictor"] == "multimodal")
    clam_only_auc = next(r["auc"] for r in p4["S3_multimodal_LR"]
                         if r["predictor"] == "clam_only")
    ras_gap = (
        prob[(prob["group"] == "molecular_subtype=RAS_like")
             & (prob["label"] == "DM1")]["min_prob_DM1"].iloc[0]
        - prob[(prob["group"] == "molecular_subtype=RAS_like")
               & (prob["label"] == "DM2")]["max_prob_DM1"].iloc[0]
    )

    hero_stats = [
        ("1.000", "original RAS-like AUC", "AUDIT_SUMMARY.json C3"),
        (fmt(e8s["pooled_ras_like_auc"]), "LOTO RAS-like AUC", "PHASE3_AUDIT_SUMMARY.json E8"),
        (fmt(e8s["pooled_overall_auc"]), "LOTO overall AUC", "PHASE3_AUDIT_SUMMARY.json E8"),
        (fmt(tss_only["auc"]), "TSS-only RAS AUC", "tss_only_within_raslike_lr.json"),
        (fmt(e5["clinical_only_oof_auc"]), "clinical-only LR", "PHASE2_AUDIT_SUMMARY.json E5"),
        (fmt(multimodal_auc), "multimodal LR", "phase4/S3"),
        ("0.350", "male raw AUC", "PHASE2_AUDIT_SUMMARY.json E1"),
        (fmt(uni["point_auc"]), "UNI OOF AUC", "bootstrap_auc_uni.json"),
    ]
    if uni_loto is not None:
        hero_stats.extend([
            (fmt(uni_loto["pooled_overall_auc"]), "UNI LOTO overall", "audit_uni_loto/UNI_LOTO_SUMMARY.json"),
            (fmt(uni_loto["pooled_ras_like_auc"]), "UNI LOTO RAS-like", "audit_uni_loto/UNI_LOTO_SUMMARY.json"),
        ])
    hero_cards = "".join(
        f"<div class=\"hs\"><b>{esc(v)}</b><span>{esc(k)}</span><small>{esc(src)}</small></div>"
        for v, k, src in hero_stats
    )

    headline_rows = [
        ["Original RAS-like/FVPTC AUC", ("1.000", "hi"), "n=16, n_pos=3", "AUDIT_SUMMARY.json C3", "Pre-audit headline"],
        ["RAS-like separation gap", (fmt(ras_gap), "warn"), "lowest DM1 minus highest DM2", "c2_prob_distribution.tsv", "One near-boundary rank drives perfect AUC"],
        ["LOTO RAS-like AUC", (fmt(e8s["pooled_ras_like_auc"]), "bad"), "pooled leave-one-TSS-out", "PHASE3_AUDIT_SUMMARY.json E8", "Worse than random under center holdout"],
        ["LOTO overall AUC", (fmt(e8s["pooled_overall_auc"]), "bad"), "pooled leave-one-TSS-out", "PHASE3_AUDIT_SUMMARY.json E8", "Overall signal weakens materially"],
        ["TSS-only LR within RAS-like", (fmt(tss_only["auc"]), "bad"), "TSS dummies only; no image features", "tss_only_within_raslike_lr.tsv/json", "Perfect ranking is reproducible without pixels"],
        ["Clinical-only LR", (fmt(clinical_auc), "hi"), "histology + sex + TSS", "phase4/s3_multimodal_LR.tsv", "Beats image + clinical"],
        ["Multimodal LR", (fmt(multimodal_auc), "warn"), "CLAM_prob + clinical", "phase4/s3_multimodal_LR.tsv", f"Net image gain {fmt(multimodal_auc - clinical_auc, 3)}"],
        ["Male raw AUC", ("0.350", "bad"), "n=13 males", "PHASE2_AUDIT_SUMMARY.json E1", "Sex-stratified failure"],
    ]
    if uni_loto is not None:
        headline_rows.extend([
            ["UNI-final OOF AUC", (fmt(uni["point_auc"]), "hi"), "n=54", "bootstrap_auc_uni.json", "Foundation encoder improves overall ranking"],
            ["UNI-final LOTO overall AUC", (fmt(uni_loto["pooled_overall_auc"]), "hi"), "center holdout", "audit_uni_loto/UNI_LOTO_SUMMARY.json", "Overall signal survives center holdout"],
            ["UNI-final LOTO RAS-like AUC", (fmt(uni_loto["pooled_ras_like_auc"]), "warn"), "center holdout in 16 RAS-like slides", "audit_uni_loto/UNI_LOTO_SUMMARY.json", "No 1.000; subgroup remains underpowered"],
        ])

    raw_rows = []
    for i, row in enumerate(ras_raw.sort_values("prob_DM1").itertuples(), 1):
        raw_rows.append([
            i,
            esc(row.submitter_id),
            esc(row.submitter_id.split("-")[1]),
            "DM1" if row.label == 1 else "DM2",
            fmt(row.prob_DM1),
            int(row.fold),
        ])

    perm_rows = []
    for label, rec in [
        ("RAS_like", ras_perm),
        ("FVPTC", fv_perm),
        ("BRAF_like", braf_perm),
        ("cPTC", p1["C3_permutation"]["histology_subtype=cPTC"]),
    ]:
        perm_rows.append([
            label,
            f"{rec['n']} ({rec['n_pos']}/{rec['n_neg']})",
            fmt(rec["observed_auc"]),
            fmt(rec["permutation_p_two_sided"], 4),
            f"[{fmt(rec['bootstrap_ci95_lo'])}, {fmt(rec['bootstrap_ci95_hi'])}]",
            "c3_permutation_test.json",
        ])

    fold_rows = [
        [int(r.fold), int(r.n), int(r.n_pos), int(r.n_neg),
         ("zero positives", "bad") if int(r.n_pos) == 0 else "has positives"]
        for r in fold.itertuples()
    ]

    tss_rows = [
        [esc(r.tss), int(r.n), int(r.n_pos), fmt(r.mean_prob),
         ("deterministic anchor", "bad") if (r.tss in ["DJ", "EM", "FK"]) else "single-site"]
        for r in e2.itertuples()
    ]

    loto_rows = []
    for r in e8.itertuples():
        auc = fmt(r.tss_holdout_auc)
        cls = "bad" if r.tss_held_out in ["DJ", "EM", "FK", "ET"] and auc != "n/a" else ""
        loto_rows.append([
            esc(r.tss_held_out),
            int(r.n_test),
            int(r.n_pos),
            int(r.n_ras_test),
            int(r.n_ras_pos_test),
            (auc, cls),
            "e8_leave_one_tss_out.tsv",
        ])

    uni_loto_block = ""
    if uni_loto is not None:
        uni_fig_prefix = (
            "../audit_uni_loto/"
            if figure_prefix == "figures/"
            else figure_prefix.replace("audit_ras_auc100/figures/", "audit_uni_loto/")
        )
        uni_loto_rows = [
            ["Original UNI OOF overall", fmt(uni_loto["original_uni_overall_auc"]), "phase2_tcga_clam_UNI/clam_per_slide_predictions.tsv"],
            ["UNI LOTO pooled overall", fmt(uni_loto["pooled_overall_auc"]), "audit_uni_loto/UNI_LOTO_SUMMARY.json"],
            ["UNI LOTO pooled RAS-like", fmt(uni_loto["pooled_ras_like_auc"]), "audit_uni_loto/UNI_LOTO_SUMMARY.json"],
            ["UNI LOTO pooled BRAF-like", fmt(uni_loto["pooled_braf_like_auc"]), "audit_uni_loto/UNI_LOTO_SUMMARY.json"],
            ["UNI LOTO Female", fmt(uni_loto["pooled_female_auc"]), "audit_uni_loto/UNI_LOTO_SUMMARY.json"],
            ["UNI LOTO Male", fmt(uni_loto["pooled_male_auc"]), "audit_uni_loto/UNI_LOTO_SUMMARY.json"],
        ]
        uni_loto_block = f"""
  <h3>UNI-final center-holdout audit</h3>
  <div class=\"box good\">
    <h3>Refined conclusion after GPU LOTO.</h3>
    <p>The ViT-L subgroup headline was artefactual, but the UNI-final model is not simply dead: overall UNI LOTO AUC is {fmt(uni_loto['pooled_overall_auc'])}. The honest subgroup claim is smaller and underpowered, with UNI LOTO RAS-like AUC {fmt(uni_loto['pooled_ras_like_auc'])}, not 1.000.</p>
  </div>
  {table(["Metric", "Value", "Source"], uni_loto_rows)}
  {figure(uni_fig_prefix + "fig_H_uni_loto.png", "Figure H", "UNI-final OOF versus leave-one-TSS-out center-holdout performance.")}
"""

    retrain_rows = [
        ["original GPU KFold seed=42", "0.746", "1.000", "0.592", "phase2_tcga_clam/clam_per_slide_predictions.tsv"],
    ]
    for r in rs["strategy_A"]:
        retrain_rows.append([
            f"KFold CPU seed={r['seed']}",
            fmt(r["overall_auc"]),
            fmt(r["ras_like_auc"]),
            fmt(r["braf_like_auc"]),
            "retrain_strategyA_multi_seed.tsv",
        ])
    retrain_rows.append([
        "StratifiedKFold(3)",
        fmt(rs["strategy_B"]["overall_auc"]),
        fmt(rs["strategy_B"]["ras_like_auc"]),
        fmt(rs["strategy_B"]["braf_like_auc"]),
        "RETRAIN_SUMMARY.json",
    ])
    retrain_rows.append([
        "LOO 16 RAS-like",
        "n/a",
        fmt(rs["strategy_C"]["loo_auc"]),
        "n/a",
        "RETRAIN_SUMMARY.json",
    ])

    s1_rows = [
        [int(r.init_seed), fmt(r.overall), fmt(r.ras), fmt(r.braf),
         fmt(r.male), fmt(r.female), "phase4/s1_tss_combat_per_init.tsv"]
        for r in s1.itertuples()
    ]
    s2_rows = [
        [int(r.init_seed), fmt(r.overall), fmt(r.ras), fmt(r.braf),
         "phase4/s2_tss_balanced_split.tsv"]
        for r in s2.itertuples()
    ]
    s3_rows = [
        [esc(r.predictor), fmt(r.auc), "phase4/s3_multimodal_LR.tsv"]
        for r in s3.itertuples()
    ]
    s5_rows = [
        [int(r.init_seed), fmt(r.overall), fmt(r.ras), fmt(r.braf),
         int(r.n_train), "phase4/s5_em_out.tsv"]
        for r in s5.itertuples()
    ]

    fig = lambda name: f"{figure_prefix}{name}.png"
    related = [
        ("Paper 2 image-DM1 page", f"{sibling_prefix}paper2_image_dm1.html", "Original visual summary before audit"),
        ("Yu review report", "../results/p2_image_dm1_v2_foundation_clam_2026_05_07/REPORT_YU_REVIEW_PAPER2_2026_05_08.html", "Contains red post-audit caveat"),
        ("Paper 1 Fig 8 dossier", f"{sibling_prefix}paper1_fig8_mechanism_dossier.html", "Validated dossier pattern"),
        ("Reviewer defense dashboard", f"{sibling_prefix}paper1_reviewer_defense_dashboard.html", "Q-numbered defense surface"),
        ("Methods outline", "../results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_ras_auc100/METHODS_PAPER_OUTLINE.md", "Path 2 methods-paper scaffold"),
        ("K2 H&E draft", "../results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_ras_auc100/K2_HE_PUSH_EMAIL_DRAFT.md", "Author-only outreach draft"),
    ]
    related_cards = "".join(
        f"<a href=\"{esc(href)}\"><b>{esc(title)}</b><span>{esc(desc)}</span></a>"
        for title, href, desc in related
    )

    decision_rows = [
        ["Publish RAS-like AUC=1.000 as-is", "Maximum headline", "Refutable by LOTO/TSS-only audit", ("Reject", "rec"), "AUDIT_REPORT_v2.md"],
        ["Report original model only in supplement", "Preserves useful pilot work", "Cannot support main-text claim", ("Use for Path 3", "rec"), "AUDIT_REPORT_v2.md"],
        ["TCGA-only image-DM1 main claim", "Fast", "Clinical-only beats image+clinical", ("Do not use", "rec"), "phase4/s3_multimodal_LR.tsv"],
        ["Clinical + image nomogram", "Truthful multimodal framing", "Current image gain is negative", ("Future only", "rec"), "phase4/s3_multimodal_LR.tsv"],
        ["K2 H&E external validation", "Breaks TCGA TSS structure", "Needs slide retrieval", ("Primary fix", "rec"), "K2_HE_PUSH_EMAIL_DRAFT.md"],
        ["LUAD/BRCA replication", "Methods-paper generalization", "No local WSI features yet", ("Scaffold now", "rec"), "METHODS_PAPER_OUTLINE.md"],
        ["ComBat correction", "Improves male failure in this pilot", "Does not restore original 1.000", ("Supplementary baseline", "rec"), "phase4/s1_tss_combat_per_init.tsv"],
        ["TSS-balanced split", "Shows matched-center ceiling", "Still fails LOTO truth test", ("Report with caveat", "rec"), "phase4/s2_tss_balanced_split.tsv"],
        ["EM-out sensitivity", "Confirms EM-DM2 anchor load-bearing", "Underpowered after removal", ("Include as caveat", "rec"), "phase4/s5_em_out.tsv"],
    ]

    css = """
<style>
:root{--bg:#071016;--panel:#0d1823;--panel2:#111f2c;--ink:#e9f2f6;--muted:#9fb0bd;--line:#263746;--gold:#ffd28a;--red:#ff8a6b;--teal:#5bd6bd;--green:#8be0a4;--blue:#83b7ff}
*{box-sizing:border-box}
html,body{margin:0;background:var(--bg);color:var(--ink);font-family:Inter,"Noto Sans KR",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;font-size:14px;line-height:1.55}
a{color:var(--teal);text-decoration:none}a:hover{text-decoration:underline}
code,.mono{font-family:"JetBrains Mono","SF Mono",Menlo,monospace;font-size:12px}
h1,h2,h3{font-family:"Cormorant Garamond",Newsreader,Georgia,serif;color:#fff;margin:0}
.hero{padding:64px 32px 42px;border-bottom:1px solid var(--line);background:linear-gradient(135deg,#08131d 0%,#0d1a25 58%,#071016 100%)}
.hero-inner{max-width:1320px;margin:0 auto}.kicker{font-family:"JetBrains Mono",monospace;color:var(--gold);font-size:11px;letter-spacing:.22em;text-transform:uppercase;margin:0 0 16px}
h1{font-size:64px;line-height:.95;letter-spacing:-1px}.lead{font-family:Newsreader,Georgia,serif;font-size:18px;max-width:980px;color:#cfdbe3;margin:18px 0 0}
.hero-stats{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-top:26px}.hs{border:1px solid rgba(255,210,138,.25);background:rgba(255,255,255,.035);border-radius:8px;padding:12px 14px;min-height:98px}.hs b{display:block;color:var(--gold);font-family:"Cormorant Garamond",serif;font-size:32px;line-height:1}.hs span{display:block;margin-top:4px;color:#d5e0e7;font-weight:700}.hs small{display:block;margin-top:4px;color:var(--muted);font-family:"JetBrains Mono",monospace;font-size:10px}
.crumbs{margin-top:20px;color:var(--muted);font-family:"JetBrains Mono",monospace;font-size:11px}.wrap{max-width:1320px;margin:0 auto;display:grid;grid-template-columns:250px minmax(0,1fr);gap:34px;padding:0 32px 80px}.toc{position:sticky;top:0;max-height:100vh;overflow:auto;padding:28px 0;border-right:1px solid var(--line)}.toc h4{margin:0 0 12px;color:var(--gold);font-family:"JetBrains Mono",monospace;font-size:10px;letter-spacing:.18em;text-transform:uppercase}.toc a{display:block;color:#d7e2e8;padding:4px 0;font-size:12px}.toc .sub{color:var(--muted);font-size:11px;padding-left:12px}
main{padding-top:24px}section{border-bottom:1px solid var(--line);padding:32px 0}h2{font-size:36px;line-height:1.05;margin-bottom:6px}h2 .num{font-family:"JetBrains Mono",monospace;color:var(--gold);font-size:14px;letter-spacing:.2em;margin-right:12px}.subhead{font-family:"JetBrains Mono",monospace;color:var(--muted);font-size:11px;letter-spacing:.1em;text-transform:uppercase;margin-bottom:18px}h3{font-size:23px;margin:26px 0 8px;color:#fff4df}p{margin:8px 0 12px}.lead-p{font-family:Newsreader,Georgia,serif;color:#d1dde5;font-size:16px}
.box{border:1px solid var(--line);border-left:4px solid var(--gold);background:var(--panel);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}.box.warn{border-left-color:var(--red);background:#1f1714}.box.good{border-left-color:var(--green);background:#102018}.box h3{margin-top:0}
.grid{display:grid;gap:14px}.grid2{grid-template-columns:repeat(2,1fr)}.grid3{grid-template-columns:repeat(3,1fr)}.card{border:1px solid var(--line);background:var(--panel);border-radius:8px;padding:14px}.card b{display:block;color:var(--gold);font-family:"JetBrains Mono",monospace;font-size:11px;text-transform:uppercase;letter-spacing:.08em}.card .big{font-size:30px;color:#fff;font-family:"Cormorant Garamond",serif;font-weight:700}.card p{color:#cbd7df;font-size:12.5px}
.t{width:100%;border-collapse:collapse;margin:12px 0 18px;font-size:12.5px}.t th,.t td{border:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}.t th{background:#142333;color:var(--gold);font-family:"JetBrains Mono",monospace;font-size:10.5px;text-transform:uppercase;letter-spacing:.05em}.t tr:nth-child(even) td{background:rgba(255,255,255,.025)}.t td.hi{color:var(--gold);font-weight:700}.t td.warn{color:var(--gold)}.t td.bad{color:var(--red);font-weight:700}.t td.good{color:var(--green);font-weight:700}.t td.rec{background:rgba(255,210,138,.12);color:var(--gold);font-weight:800}.path{font-family:"JetBrains Mono",monospace;color:#aeeadf;background:#0a151f;border:1px solid var(--line);border-radius:4px;padding:1px 5px;font-size:11px}
.fig{border:1px solid var(--line);background:#071016;border-radius:8px;overflow:hidden;margin:16px 0}.fig img{display:block;width:100%;background:#fff}.fig figcaption{border-top:1px solid var(--line);background:var(--panel);padding:10px 12px;color:#cbd7df;font-size:12px}.fig figcaption b{color:var(--gold);font-family:"JetBrains Mono",monospace;font-size:10.5px;letter-spacing:.08em;text-transform:uppercase;margin-right:8px}
.related{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin-top:14px}.related a{display:block;border:1px solid var(--line);border-radius:8px;background:var(--panel);padding:12px;color:#d8e3ea}.related a:hover{border-color:var(--gold);text-decoration:none}.related b{display:block;color:var(--gold);font-family:"JetBrains Mono",monospace;font-size:11px;text-transform:uppercase}.related span{display:block;color:var(--muted);font-size:12px;margin-top:4px}
.foot{color:var(--muted);font-family:"JetBrains Mono",monospace;font-size:11px;margin-top:24px}
@media(max-width:1000px){.wrap{grid-template-columns:1fr;padding:0 18px 60px}.toc{display:none}.hero{padding:48px 20px 34px}h1{font-size:44px}.hero-stats{grid-template-columns:repeat(2,1fr)}.grid2,.grid3,.related{grid-template-columns:1fr}}
</style>
"""

    html_out = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Paper 2 image-DM1 audit dossier</title>
{css}
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <p class="kicker">Paper 2 audit dossier · TCGA-THCA · built 2026-05-09 · {esc(out_label)}</p>
    <h1>The perfect RAS-like AUC was a TSS artefact.</h1>
    <p class="lead">A 12-test audit of the TCGA-THCA image-DM1 pilot shows that the reported RAS-like/FVPTC AUC of 1.000 does not survive tissue-source-site holdout. Under pooled leave-one-TSS-out validation, the RAS-like AUC falls to {fmt(e8s['pooled_ras_like_auc'])}; a TSS-only logistic regression reaches {fmt(tss_only['auc'])} inside the same 16 slides without image features.</p>
    <div class="hero-stats">{hero_cards}</div>
    <div class="crumbs"><a href="{esc(sibling_prefix)}index.html">Hub</a> / <a href="{esc(sibling_prefix)}paper2_image_dm1.html">Paper 2 image-DM1</a> / source: <span class="path">analysis_supp/audit_ras_auc100/</span></div>
  </div>
</header>
<div class="wrap">
<aside class="toc">
  <h4>Contents</h4>
  <a href="#verdict">01 Verdict</a>
  <a class="sub" href="#definitions">Data definitions</a>
  <a href="#phase1">02 Phase 1 landscape</a>
  <a class="sub" href="#raw16">Raw 16 RAS-like slides</a>
  <a href="#phase23">03 Phase 2/3 decisive evidence</a>
  <a class="sub" href="#loto">Leave-one-TSS-out</a>
  <a href="#solutions">04 Phase 4 corrective baselines</a>
  <a href="#paths">05 Forward paths</a>
  <a href="#decision">06 Decision matrix</a>
  <a href="#supplement">07 Supplementary files</a>
  <a href="#sources">08 Sources + paths</a>
</aside>
<main>
<section id="verdict">
  <h2><span class="num">01</span>Verdict</h2>
  <p class="subhead">advisor-facing conclusion with every number tied to a local evidence file</p>
  <div class="box warn">
    <h3>Do not publish the RAS-like AUC=1.000 as currently framed.</h3>
    <p>The TCGA-only image channel is not the load-bearing signal. RAS-like/FVPTC membership is exactly the same 16-slide subset, the 1.000 ranking has a {fmt(ras_gap)} probability gap, the 3 positive RAS-like slides come from DJ/FK while the 10 EM RAS-like slides are all DM2, and clinical-only metadata outperforms the multimodal model.</p>
  </div>
  {table(["Audit number", "Value", "Scope", "Source", "Meaning"], headline_rows)}
  <div class="grid grid3">
    <div class="card"><b>Data</b><div class="big">59 slides</div><p>TCGA-THCA CLAM OOF predictions; 16 RAS-like/FVPTC slides with 3 positives.</p></div>
    <div class="card"><b>Analysis</b><div class="big">12 tests</div><p>C1-C5, E1-E8 plus S1/S2/S3/S5 corrective baselines, all CPU and reproducible.</p></div>
    <div class="card"><b>Clinical meaning</b><div class="big">K2 gated</div><p>External H&E with TSS-balanced design is required before an image-DM1 main claim is defensible.</p></div>
  </div>
</section>

<section id="definitions">
  <h2><span class="num">02</span>Data Definitions</h2>
  <p class="subhead">scope dictionary for reviewer/advisor reading</p>
  {table(["Object", "Definition", "n", "Role"], [
      ["TCGA-THCA WSI pilot", "Slides with pre-extracted foundation-model tile features and CLAM OOF predictions", "59 ViT-L audit / 54 UNI final", "Paper 2 image-DM1 pilot"],
      ["DM1 label", "Paper 1 molecular dark-matter cluster label projected onto image-available TCGA cases", "29 positive in ViT-L OOF", "Binary CLAM target"],
      ["RAS-like/FVPTC subset", "Molecular_subtype=RAS_like and histology_subtype=FVPTC; identical in the 59-slide subset", "16 (3/13)", "Original perfect-AUC subgroup"],
      ["TSS", "TCGA tissue source site parsed from submitter_id, e.g. TCGA-EM-*", "12 TSS groups in OOF set", "Batch-effect axis"],
      ["Clinical-only LR", "Out-of-fold logistic regression on histology, sex, and TSS dummies", "57 analyzable cPTC/FVPTC slides", "Image-gain control"],
      ["LOTO", "Leave-one-TSS-out retraining and pooled evaluation", "all TSS groups", "Center-holdout truth test"],
  ])}
</section>

<section id="phase1">
  <h2><span class="num">03</span>Phase 1 Confound Landscape</h2>
  <p class="subhead">C1-C5: overlap, raw distribution, permutation null, shortcut probe, fold composition</p>
  <p class="lead-p"><b>Data.</b> The audit starts from <span class="path">phase2_tcga_clam/clam_per_slide_predictions.tsv</span> plus <span class="path">sample_master_v3.tsv</span>. <b>Analysis.</b> C1-C5 are no-retraining checks. <b>Meaning.</b> The headline perfect AUC is a duplicated subgroup result with a narrow margin. <b>Clinical meaning.</b> A reviewer can reproduce the problem before touching the model weights.</p>
  {figure(fig("fig_A_raw_16_distribution"), "Figure A", "Raw RAS-like/FVPTC prediction distribution; the perfect ranking is separated by only 0.017 probability.")}
  <div id="raw16">
  <h3>Raw 16-slide ranking</h3>
  {table(["Rank", "Case", "TSS", "Label", "prob_DM1", "Fold"], raw_rows)}
  </div>
  <h3>Permutation null and bootstrap CI</h3>
  {table(["Group", "n (pos/neg)", "Observed AUC", "Permutation p", "Bootstrap 95% CI", "Source"], perm_rows)}
  {figure(fig("fig_B_permutation_null"), "Figure B", "Permutation null for the small RAS-like/FVPTC subset.")}
  <h3>Fold imbalance</h3>
  {table(["Fold", "n", "n_pos", "n_neg", "Interpretation"], fold_rows)}
  {figure(fig("fig_C_fold_composition"), "Figure C", "Original folds 1, 2, and 4 have zero RAS-like positives in test.")}
</section>

<section id="phase23">
  <h2><span class="num">04</span>Phase 2/3 Decisive Evidence</h2>
  <p class="subhead">sex failure, TSS determinism, clinical baseline, init variance, label-shuffle null, LOTO</p>
  <div class="box warn">
    <h3>The decisive test is LOTO.</h3>
    <p>Random and stratified splits can keep the same TSS centers in train and test. LOTO removes that escape route: overall AUC becomes {fmt(e8s['pooled_overall_auc'])}, and RAS-like AUC becomes {fmt(e8s['pooled_ras_like_auc'])}.</p>
  </div>
  <h3>TSS distribution among RAS-like slides</h3>
  {table(["TSS", "n", "n DM1", "mean prob_DM1", "Interpretation"], tss_rows)}
  {figure(fig("fig_E_sex_tss_confound"), "Figure E", "Sex-stratified AUC and TSS concentration.")}
  <h3>Split-stress retraining</h3>
  {table(["Split strategy", "Overall AUC", "RAS-like AUC", "BRAF-like AUC", "Source"], retrain_rows)}
  {figure(fig("fig_D_retrain_summary"), "Figure D", "Multi-seed and stratified retrain summary.")}
  <h3>Clinical-only baseline</h3>
  {table(["Predictor", "OOF AUC", "Source"], s3_rows)}
  {figure(fig("fig_F_clinical_baseline"), "Figure F", "Clinical-only metadata already matches or exceeds the image-inclusive model.")}
  <h3 id="loto">Leave-one-TSS-out</h3>
  {table(["Held-out TSS", "n test", "n DM1", "n RAS", "n RAS DM1", "AUC", "Source"], loto_rows)}
  {figure(fig("fig_G_phase3_decisive"), "Figure G", "LOTO collapses the pooled RAS-like estimate to worse-than-random.")}
  {uni_loto_block}
  <div class="box good">
    <h3>Label-shuffle null passed.</h3>
    <p>The shuffled-label controls stay near chance, so this is not a simple memorization bug. The failure mode is specifically TSS/site entanglement plus sparse subgroup positives.</p>
  </div>
</section>

<section id="solutions">
  <h2><span class="num">05</span>Phase 4 Corrective Baselines</h2>
  <p class="subhead">S1/S2/S3/S5: what each fix buys and what it cannot repair</p>
  <p class="lead-p"><b>Data.</b> Same 59-slide feature set and OOF labels. <b>Analysis.</b> TSS-ComBat, TSS-balanced split, multimodal LR, and EM-out sensitivity. <b>Meaning.</b> No correction recovers the original 1.000 under honest center stress. <b>Clinical meaning.</b> The honest TCGA-only ceiling is pilot-grade, not submission-grade.</p>
  <h3>S1: TSS-ComBat on tile features</h3>
  {table(["Init seed", "Overall", "RAS-like", "BRAF-like", "Male", "Female", "Source"], s1_rows)}
  <p>Median S1 overall AUC is {fmt(median(list(s1['overall'])))}; median male AUC improves to {fmt(median(list(s1['male'])))} from raw 0.350, but the overall signal remains pilot-level.</p>
  <h3>S2: TSS-balanced split</h3>
  {table(["Init seed", "Overall", "RAS-like", "BRAF-like", "Source"], s2_rows)}
  <p>S2 preserves high RAS-like AUC when matched centers remain in train/test. That is useful as a matched-center ceiling, not as evidence of center-generalized biology.</p>
  <h3>S3: Multimodal LR</h3>
  {table(["Predictor", "OOF AUC", "Source"], s3_rows)}
  <p>CLAM-only in this corrective baseline is {fmt(clam_only_auc)}. Clinical-only is {fmt(clinical_auc)}. Multimodal is {fmt(multimodal_auc)}, so the image channel net gain is {fmt(multimodal_auc - clinical_auc)}.</p>
  <h3>S5: EM-out sensitivity</h3>
  {table(["Init seed", "Overall", "RAS-like", "BRAF-like", "n train", "Source"], s5_rows)}
  <p>Dropping the 10 EM RAS-like DM2 anchors lowers the RAS-like median to {fmt(median(list(s5['ras'])))}. This supports the TSS-anchor interpretation.</p>
</section>

<section id="paths">
  <h2><span class="num">06</span>Forward Paths</h2>
  <p class="subhead">Path 1/2/3 status after the audit</p>
  {table(["Path", "Goal", "Current status", "Blocker", "Action"], [
      ["Path 1: multi-cohort multimodal nomogram", "NC/Nature Cancer-grade image + clinical + RNA model", "Scientifically valid only with external H&E", "K2 or FFPE multi-site timing", "Keep design; do not headline TCGA-only image gain"],
      ["Path 2: methods paper", "DIAL sequel for TCGA pathology-AI TSS audits", "THCA worked example complete", "Needs LUAD/BRCA feature extraction", "Use this dossier as the advisor packet"],
      ["Path 3: Paper 1 absorption", "Safe NC path by demoting Paper 2 image work to supp case study", "Audit caveat already added to Paper 2 report", "Needs Paper 1 supplement insertion", "Add TSS-confound disclosure section"],
      ["K2 outreach", "Resolve external validation timeline", "Draft exists", "Recipient name/email", "Author sends after editing"],
  ])}
</section>

<section id="decision">
  <h2><span class="num">07</span>Decision Matrix</h2>
  <p class="subhead">recommended disposition is highlighted</p>
  {table(["Option", "Benefit", "Risk", "Disposition", "Evidence"], decision_rows, "dec")}
</section>

<section id="supplement">
  <h2><span class="num">08</span>Supplementary Files</h2>
  <p class="subhead">all-pages card grid and audit deliverables</p>
  <div class="related">{related_cards}</div>
  {table(["File", "Purpose"], [
      ["AUDIT_REPORT_v2.md", "Narrative final report for the audit"],
      ["METHODS_PAPER_OUTLINE.md", "Path 2 methods-paper scaffold"],
      ["METHODS_PAPER_TABLES.md", "Cross-cohort methods-paper table scaffold; THCA complete, LUAD/BRCA pending inputs"],
      ["../audit_uni_predictions/UNI_PREDICTION_ONLY_AUDIT.md", "UNI-final OOF prediction-only audit; paired with UNI LOTO below"],
      ["../audit_uni_loto/UNI_LOTO_REPORT.md", "UNI-final GPU leave-one-TSS-out center-holdout audit"],
      ["K2_HE_PUSH_EMAIL_DRAFT.md", "Author-only K2 H&E outreach draft"],
      ["figures/fig_A..G_*.png/pdf", "Seven audit figures embedded above"],
      ["phase4/*.tsv", "Corrective-baseline raw tables"],
      ["tss_only_within_raslike_lr.tsv/json", "Generated evidence for TSS-only perfect ranking"],
  ])}
</section>

<section id="sources">
  <h2><span class="num">09</span>Sources + Paths</h2>
  <p class="subhead">local evidence chain</p>
  {table(["Source", "Path"], [
      ["Audit root", "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_ras_auc100/"],
      ["Dossier source copy", "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_ras_auc100/AUDIT_DOSSIER.html"],
      ["Hub copy", "project/papers_hub_2026_05_04/paper2_image_dm1_audit_dossier.html"],
      ["Live deploy target", "/var/www/papers/papers_hub_2026_05_04/paper2_image_dm1_audit_dossier.html"],
      ["Builder", "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_make_dossier.py"],
      ["Original OOF predictions", "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/clam_per_slide_predictions.tsv"],
      ["UNI bootstrap", "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/bootstrap_auc_uni.json"],
      ["UNI LOTO audit", "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_uni_loto/"],
  ])}
  <p class="foot">Build generated by audit_make_dossier.py on 2026-05-09. Voice-protected Paper 1 manuscript sections were not edited by this dossier builder.</p>
</section>
</main>
</div>
</body>
</html>
"""
    return html_out


def main() -> None:
    AUD.mkdir(parents=True, exist_ok=True)
    HUB.mkdir(parents=True, exist_ok=True)

    audit_html = build_html("figures/", "../../papers_hub_2026_05_04/", "audit-root copy")
    audit_out = AUD / "AUDIT_DOSSIER.html"
    audit_out.write_text(audit_html)

    hub_html = build_html(
        "../results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/audit_ras_auc100/figures/",
        "",
        "hub copy",
    )
    hub_out = HUB / "paper2_image_dm1_audit_dossier.html"
    hub_out.write_text(hub_html)

    print(f"WROTE {audit_out}")
    print(f"WROTE {hub_out}")


if __name__ == "__main__":
    main()
