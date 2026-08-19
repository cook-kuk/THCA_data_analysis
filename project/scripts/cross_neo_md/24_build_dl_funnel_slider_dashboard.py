#!/usr/bin/env python3
"""Build an interactive DL-first funnel dashboard with threshold sliders."""

from __future__ import annotations

import html
import importlib.util
import json
import shutil
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[3]
SCRIPT_DIR = REPO / "project/scripts/cross_neo_md"
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
OUT = MD_OUT / "dl_first_funnel_ui"
HUB = REPO / "project/papers_hub_2026_05_04"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
PAGE = HUB / "cross_neo_dl_first_slider_dashboard.html"
WEB_PAGE = WEB / "cross_neo_dl_first_slider_dashboard.html"


DEFAULTS = {
    "class_i_min_len": 8,
    "class_i_max_len": 11,
    "other_min_len": 9,
    "other_max_len": 25,
    "main_dl": 0.45,
    "ensemble_upper": 0.75,
    "pmhc": 0.65,
    "bayes_mean": 0.50,
    "perturb_prob": 0.50,
    "bayes_upper": 0.70,
    "robust_bayes": 0.55,
    "robust_perturb": 0.70,
    "robust_main": 0.75,
    "dropout_sens": 0.12,
    "perturb_width": 0.45,
    "min_tcr_evidence": 1,
    "min_paired_tcr": 1,
    "min_cancer_context": 1,
    "max_pathogen_context": 0,
    "tcr_branch": 0.65,
    "tcr_rescue": 0.85,
    "tcr_rescue_bayes_upper": 0.50,
    "min_baker_structural": 0.70,
    "min_md_structural": 0.45,
    "min_live_completion": 0.90,
}

SLIDER_META = {
    "class_i_min_len": {"min": 5, "max": 15, "step": 1},
    "class_i_max_len": {"min": 6, "max": 25, "step": 1},
    "other_min_len": {"min": 5, "max": 30, "step": 1},
    "other_max_len": {"min": 8, "max": 40, "step": 1},
    "min_tcr_evidence": {"min": 0, "max": 120, "step": 1},
    "min_paired_tcr": {"min": 0, "max": 30, "step": 1},
    "min_cancer_context": {"min": 0, "max": 120, "step": 1},
    "max_pathogen_context": {"min": 0, "max": 500, "step": 1},
}


def load_dl_module():
    path = SCRIPT_DIR / "22_dl_first_candidate_funnel.py"
    spec = importlib.util.spec_from_file_location("dl_first_candidate_funnel", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def boolify(x: object) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"true", "1", "yes"}


def add_baker_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    score_path = MD_OUT / "baker_rosetta_filter/baker_rosetta_fallback_interface_scores.tsv"
    if score_path.exists():
        scores = pd.read_csv(score_path, sep="\t")
        if {"row_id", "fallback_structural_score"}.issubset(scores.columns):
            agg = (
                scores.groupby("row_id", as_index=False)
                .agg(baker_structural_score=("fallback_structural_score", "max"))
            )
            out = out.merge(agg, on="row_id", how="left")
    if "baker_structural_score" not in out.columns:
        out["baker_structural_score"] = 0.0
    out["baker_structural_score"] = pd.to_numeric(out["baker_structural_score"], errors="coerce").fillna(0.0)
    return out


def stage_label_counts(df: pd.DataFrame) -> pd.DataFrame:
    baker_pass = pd.to_numeric(df.get("baker_structural_score", pd.Series(0, index=df.index)), errors="coerce").fillna(0) >= DEFAULTS["min_baker_structural"]
    md_score = pd.to_numeric(df.get("md_structural_score", pd.Series(0, index=df.index)), errors="coerce").fillna(0)
    md_score = md_score.combine(pd.to_numeric(df.get("md_score", pd.Series(0, index=df.index)), errors="coerce").fillna(0), max)
    live_completion = pd.to_numeric(df.get("live_completion_fraction", pd.Series(0, index=df.index)), errors="coerce").fillna(0)
    md_pass = (md_score >= DEFAULTS["min_md_structural"]) | (live_completion >= DEFAULTS["min_live_completion"])
    stages = [
        ("00_all_candidates", pd.Series(True, index=df.index)),
        ("01_supported_peptide_length", df["supported_length_main_dl"].map(boolify)),
        ("02_main_dl_broad_pass", df["stage1_main_dl_broad_pass"].map(boolify)),
        ("03_bayesian_dropout_pass", df["stage2_uncertainty_pass"].map(boolify)),
        ("04_robust_main_dl_pass", df["stage3_robust_dl_pass"].map(boolify)),
        ("05_source_compatible_tcr_after_dl", df["stage3_robust_dl_pass"].map(boolify) & df["source_compatible_tcr"].map(boolify)),
        ("06_paired_tcr_after_dl", df["stage3_robust_dl_pass"].map(boolify) & df["paired_source_compatible_tcr"].map(boolify)),
        (
            "07_tcr_branch_supported_after_dl",
            df["stage3_robust_dl_pass"].map(boolify) & df["source_compatible_tcr"].map(boolify) & df["tcr_branch_support"].map(boolify),
        ),
        (
            "08_baker_structure_supported",
            df["stage3_robust_dl_pass"].map(boolify) & df["source_compatible_tcr"].map(boolify) & df["tcr_branch_support"].map(boolify) & baker_pass,
        ),
        (
            "09_md_supported_or_complete",
            df["stage3_robust_dl_pass"].map(boolify) & df["source_compatible_tcr"].map(boolify) & df["tcr_branch_support"].map(boolify) & md_pass,
        ),
        ("10_tcr_rescue_reviews", df["dl_first_decision"].eq("TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL")),
        (
            "11_md_escalation_after_cheap_gates",
            df["dl_first_decision"].isin(["DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR", "TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL"]),
        ),
        ("12_wetlab_shortlist_after_dl", df["dl_first_decision"].eq("DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR") & df["wt_or_decoy_ready"].map(boolify)),
    ]
    rows = []
    labels = pd.to_numeric(df["label_binary"], errors="coerce").fillna(-1).astype(int)
    for stage, mask in stages:
        sub = df[mask]
        lab = labels[mask]
        pos = int((lab == 1).sum())
        neg = int((lab == 0).sum())
        unknown = int((lab < 0).sum())
        rows.append(
            {
                "stage": stage,
                "n": int(len(sub)),
                "positive": pos,
                "negative": neg,
                "unknown": unknown,
                "positive_fraction": pos / len(sub) if len(sub) else 0.0,
            }
        )
    return pd.DataFrame(rows)


def decision_label_counts(df: pd.DataFrame) -> pd.DataFrame:
    labels = pd.to_numeric(df["label_binary"], errors="coerce").fillna(-1).astype(int)
    rows = []
    for decision, idx in df.groupby("dl_first_decision").groups.items():
        lab = labels.loc[idx]
        rows.append(
            {
                "decision": decision,
                "n": int(len(idx)),
                "positive": int((lab == 1).sum()),
                "negative": int((lab == 0).sum()),
                "unknown": int((lab < 0).sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(["n", "decision"], ascending=[False, True])


def js_rows(df: pd.DataFrame) -> list[dict[str, object]]:
    keep = [
        "row_id",
        "peptide",
        "hla_4digit",
        "label_binary",
        "peptide_length",
        "main_dl_score",
        "ensemble_q95",
        "pmhc_score_mean",
        "bayes_mean",
        "bayes_q05",
        "bayes_q95",
        "perturb_prob_gt_050",
        "dropout_sensitivity",
        "perturb_width_90",
        "tcr_augmented_score_mean",
        "best_tcr_augmented_score",
        "tcr_evidence_count",
        "paired_tcr_evidence_count",
        "cancer_context_evidence_count",
        "pathogen_context_evidence_count",
        "md_score",
        "md_structural_score",
        "live_completion_fraction",
        "control_readiness_score",
        "baker_structural_score",
        "source_compatible_tcr",
        "paired_source_compatible_tcr",
        "wt_or_decoy_ready",
        "md_label",
    ]
    present = [c for c in keep if c in df.columns]
    small = df[present].copy()
    for c in small.columns:
        if small[c].dtype == object:
            small[c] = small[c].fillna("")
        else:
            small[c] = pd.to_numeric(small[c], errors="coerce").fillna(0)
    small["source_compatible_tcr"] = small.get("source_compatible_tcr", False).map(boolify)
    small["paired_source_compatible_tcr"] = small.get("paired_source_compatible_tcr", False).map(boolify)
    small["wt_or_decoy_ready"] = small.get("wt_or_decoy_ready", False).map(boolify)
    return small.to_dict("records")


def esc(x: object) -> str:
    return html.escape(str(x))


def build_html(rows: list[dict[str, object]], stage_counts: pd.DataFrame, decision_counts: pd.DataFrame) -> str:
    data_json = json.dumps(rows, ensure_ascii=False)
    defaults_json = json.dumps(DEFAULTS)
    slider_meta_json = json.dumps(SLIDER_META)
    preset_path = MD_OUT / "dl_threshold_optimization/threshold_recommended_presets.json"
    presets_json = preset_path.read_text() if preset_path.exists() else "[]"
    stage_table = stage_counts.to_html(index=False, classes="data", border=0)
    decision_table = decision_counts.to_html(index=False, classes="data", border=0)
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo DL-First Funnel Slider</title>
<style>
:root{{--bg:#07111f;--panel:#0e1b2d;--line:#29405f;--text:#edf5ff;--muted:#9fb0c7;--gold:#f2c46d;--blue:#4d96ff;--red:#ff7b72;--green:#2fbf71}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--text);font-family:Inter,Arial,sans-serif;line-height:1.45}}
header{{padding:28px 30px;border-bottom:1px solid var(--line);background:#0f2035}} h1{{margin:0;font-family:Georgia,serif;font-size:34px}}
main{{display:grid;grid-template-columns:360px 1fr;gap:18px;max-width:1500px;margin:0 auto;padding:18px}}
aside,section{{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:14px}}
.controls{{position:sticky;top:12px;max-height:calc(100vh - 24px);overflow:auto}} label{{display:block;margin:12px 0 4px;color:var(--gold);font-size:12px}}
input[type=range]{{width:100%}} .value{{float:right;color:var(--text)}} .stats{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-bottom:14px}}
select{{width:100%;background:#0a1626;color:var(--text);border:1px solid var(--line);border-radius:6px;padding:8px;margin-top:4px}}
.presetGrid{{display:grid;grid-template-columns:1fr;gap:7px;margin:10px 0 14px}} .preset{{text-align:left;background:#13243a;border:1px solid #3e5c82;border-radius:7px;padding:8px;cursor:pointer;color:var(--text)}} .preset b{{color:var(--gold)}} .preset span{{display:block;color:var(--muted);font-size:11px}}
.stat{{background:#0a1626;border:1px solid var(--line);border-radius:8px;padding:10px}} .stat b{{display:block;color:var(--gold);font-size:24px}}
.confusion{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:12px 0}} .cell{{border:1px solid var(--line);border-radius:8px;padding:12px;background:#0a1626}} .cell b{{display:block;font-size:28px}}
.tp b{{color:var(--green)}} .tn b{{color:#7cb7ff}} .fp b{{color:var(--red)}} .fn b{{color:#ffcc66}}
.badge{{display:inline-block;border-radius:999px;padding:2px 7px;font-size:11px;border:1px solid var(--line);white-space:nowrap}} .label-pos,.pred-pos,.out-tp{{background:#123d2a;color:#8ff0b6}} .label-neg,.pred-neg,.out-tn{{background:#102848;color:#9dccff}} .out-fp{{background:#4a1d24;color:#ffb0b9}} .out-fn{{background:#473514;color:#ffe19a}}
.bars{{display:grid;gap:8px}} .barrow{{display:grid;grid-template-columns:260px 1fr 80px;gap:10px;align-items:center}}
.bar{{height:26px;background:#101d30;border:1px solid var(--line);border-radius:6px;overflow:hidden;display:flex}} .pos{{background:var(--green)}} .neg{{background:var(--red)}} .unk{{background:#777}}
table{{border-collapse:collapse;width:100%;font-size:12px;margin-top:10px}} th,td{{border:1px solid var(--line);padding:6px;vertical-align:top}} th{{background:#13243a;color:var(--gold)}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-top:14px}} .muted{{color:var(--muted)}} a{{color:#39d4b5}}
button{{background:#18314f;color:var(--text);border:1px solid #3e5c82;border-radius:6px;padding:8px 10px;cursor:pointer;margin-top:10px}}
@media(max-width:1000px){{main{{grid-template-columns:1fr}}.controls{{position:relative}}.stats,.grid{{grid-template-columns:1fr}}}}
</style></head><body>
<header><h1>CROSS-Neo DL-First Funnel Slider</h1>
<p class="muted">Actual current run counts are shown first. Move thresholds to see candidate, positive, and negative counts change before structure/MD escalation.</p>
<p><a href="cross_neo_md_audit_dossier.html">Back to MD dossier</a></p></header>
<main>
<aside class="controls">
<h2>Thresholds</h2>
<label for="callRule">Positive call rule <span class="value">for TP/TN/FP/FN</span></label>
<select id="callRule">
<option value="md_escalation" selected>MD escalation candidates</option>
<option value="wetlab_shortlist">Wetlab shortlist</option>
<option value="structure_md_supported">Structure/MD supported</option>
<option value="tcr_supported">TCR-supported robust DL</option>
<option value="robust_dl">Robust DL pass</option>
<option value="non_culled">All non-culled/review rows</option>
</select>
<label>Optimized presets <span class="value">label-aware</span></label>
<div class="presetGrid" id="presetButtons"></div>
<div id="sliders"></div>
<button id="reset">Reset current run thresholds</button>
<p class="muted">MD is not used to rescue failed cheap gates. It remains a late audit layer.</p>
</aside>
<div>
<section>
<h2>Live Funnel</h2>
<div class="stats" id="stats"></div>
<div class="bars" id="bars"></div>
</section>
<section>
<h2>Real-Time Label Audit</h2>
<p class="muted">TP/TN/FP/FN update when thresholds or the positive-call rule change. Labels are the dataset labels; model-positive means the selected rule calls the row positive.</p>
<div class="confusion" id="confusion"></div>
<div id="metrics" class="muted"></div>
</section>
<div class="grid">
<section><h2>Current Run Stage Counts</h2>{stage_table}</section>
<section><h2>Current Run Decision Counts</h2>{decision_table}</section>
</div>
<section>
<h2>Top Remaining Rows</h2>
<table id="top"><thead><tr><th>row</th><th>peptide</th><th>HLA</th><th>actual label</th><th>model call</th><th>outcome</th><th>decision</th><th>main DL</th><th>Bayes</th><th>TCR branch</th><th>Baker</th><th>MD</th></tr></thead><tbody></tbody></table>
</section>
</div>
</main>
<script>
const rows = {data_json};
const defaults = {defaults_json};
const sliderMeta = {slider_meta_json};
const presets = {presets_json};
const labels = {{
  class_i_min_len: "class-I peptide min length",
  class_i_max_len: "class-I peptide max length",
  other_min_len: "other peptide min length",
  other_max_len: "other peptide max length",
  main_dl: "main DL broad threshold",
  ensemble_upper: "ensemble upper rescue",
  pmhc: "pMHC score rescue",
  bayes_mean: "Bayesian mean",
  perturb_prob: "perturb prob > 0.50",
  bayes_upper: "Bayesian upper rescue",
  robust_bayes: "robust Bayesian mean",
  robust_perturb: "robust perturb prob",
  robust_main: "robust main DL",
  dropout_sens: "max dropout sensitivity",
  perturb_width: "max perturb width",
  min_tcr_evidence: "min TCR evidence rows",
  min_paired_tcr: "min paired TCR rows",
  min_cancer_context: "min cancer-context TCR rows",
  max_pathogen_context: "max pathogen-context TCR rows",
  tcr_branch: "TCR branch support",
  tcr_rescue: "TCR rescue threshold",
  tcr_rescue_bayes_upper: "TCR rescue Bayesian upper",
  min_baker_structural: "min Baker/static structure score",
  min_md_structural: "min MD structural score",
  min_live_completion: "min live MD completion"
}};
const sliderBox = document.getElementById("sliders");
const vals = {{...defaults}};
let callRule = "md_escalation";
document.getElementById("callRule").addEventListener("change", e => {{ callRule = e.target.value; render(); }});
const presetBox = document.getElementById("presetButtons");
for (const p of presets) {{
  const b = document.createElement("button");
  b.className = "preset";
  const m = p.metrics || {{}};
  b.innerHTML = `<b>${{String(p.name||"preset").replaceAll("_"," ")}}</b><span>${{p.call_rule}} · TP ${{m.TP??0}} / FP ${{m.FP??0}} · precision ${{Number(m.precision||0).toFixed(3)}} · recall ${{Number(m.recall||0).toFixed(3)}}</span>`;
  b.onclick = () => {{
    callRule = p.call_rule || callRule;
    document.getElementById("callRule").value = callRule;
    for (const [k,v] of Object.entries(p.thresholds || {{}})) {{
      if (!(k in vals) || !document.getElementById(k)) continue;
      vals[k] = Number(v);
      document.getElementById(k).value = vals[k];
      const meta = sliderMeta[k] || {{step:0.01}};
      document.getElementById("v_"+k).textContent = meta.step >= 1 ? Number(vals[k]).toFixed(0) : Number(vals[k]).toFixed(2);
    }}
    render();
  }};
  presetBox.appendChild(b);
}}
for (const [k,v] of Object.entries(defaults)) {{
  const wrap = document.createElement("div");
  const meta = sliderMeta[k] || {{min:0, max:1, step:0.01}};
  const fmt = x => meta.step >= 1 ? Number(x).toFixed(0) : Number(x).toFixed(2);
  wrap.innerHTML = `<label>${{labels[k]}} <span class="value" id="v_${{k}}">${{fmt(v)}}</span></label><input id="${{k}}" type="range" min="${{meta.min}}" max="${{meta.max}}" step="${{meta.step}}" value="${{v}}">`;
  sliderBox.appendChild(wrap);
  wrap.querySelector("input").addEventListener("input", e => {{ vals[k] = Number(e.target.value); document.getElementById("v_"+k).textContent = fmt(vals[k]); render(); }});
}}
document.getElementById("reset").onclick = () => {{
  for (const [k,v] of Object.entries(defaults)) {{
    const meta = sliderMeta[k] || {{step:0.01}};
    vals[k]=v; document.getElementById(k).value=v; document.getElementById("v_"+k).textContent=meta.step >= 1 ? Number(v).toFixed(0) : Number(v).toFixed(2);
  }}
  render();
}};
function num(x) {{ const y=Number(x); return Number.isFinite(y)?y:0; }}
function supported(r) {{
  const h = String(r.hla_4digit||"").toUpperCase();
  const l = num(r.peptide_length);
  return /^HLA-[ABC]\\*/.test(h) ? (l>=vals.class_i_min_len && l<=vals.class_i_max_len) : (l>=vals.other_min_len && l<=vals.other_max_len);
}}
function classify(r) {{
  const supp = supported(r);
  const dropoutUnstable = num(r.dropout_sensitivity) >= vals.dropout_sens || num(r.perturb_width_90) >= vals.perturb_width;
  const stage1 = supp && (num(r.main_dl_score)>=vals.main_dl || num(r.ensemble_q95)>=vals.ensemble_upper || num(r.pmhc_score_mean)>=vals.pmhc);
  const stage2 = stage1 && (((num(r.bayes_mean)>=vals.bayes_mean) && (num(r.perturb_prob_gt_050)>=vals.perturb_prob) && !dropoutUnstable) || num(r.bayes_q95)>=vals.bayes_upper || num(r.main_dl_score)>=vals.robust_main);
  const stage3 = stage2 && !dropoutUnstable && (num(r.bayes_mean)>=vals.robust_bayes || num(r.perturb_prob_gt_050)>=vals.robust_perturb || num(r.main_dl_score)>=vals.robust_main);
  const sourceTcr = num(r.tcr_evidence_count) >= vals.min_tcr_evidence && num(r.cancer_context_evidence_count) >= vals.min_cancer_context && num(r.pathogen_context_evidence_count) <= vals.max_pathogen_context;
  const pairedTcr = sourceTcr && num(r.paired_tcr_evidence_count) >= vals.min_paired_tcr;
  const tcrSupport = num(r.tcr_augmented_score_mean)>=vals.tcr_branch || num(r.best_tcr_augmented_score)>=vals.tcr_rescue;
  const tcrRescue = supp && pairedTcr && num(r.tcr_augmented_score_mean)>=vals.tcr_rescue && num(r.best_tcr_augmented_score)>=vals.tcr_rescue && num(r.bayes_q95)>=vals.tcr_rescue_bayes_upper;
  const bakerPass = num(r.baker_structural_score) >= vals.min_baker_structural;
  const mdPass = num(r.md_structural_score) >= vals.min_md_structural || num(r.md_score) >= vals.min_md_structural || num(r.live_completion_fraction) >= vals.min_live_completion;
  let decision = "";
  if (!supp) decision = "DL_CULL_UNSUPPORTED_PEPTIDE_LENGTH";
  else if (!stage1) decision = "DL_CULL_LOW_MAIN_MODEL_SCORE";
  else if (!stage2) decision = tcrRescue ? "TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL" : "DL_CULL_UNCERTAIN_OR_LOW_POSTERIOR";
  else if (dropoutUnstable && !pairedTcr) decision = "DL_CULL_DROPOUT_UNSTABLE";
  else if (stage3 && sourceTcr && !tcrSupport) decision = "DL_TCR_DISCORDANT_AUDIT_NOT_POSITIVE";
  else if (stage3 && pairedTcr && tcrSupport) decision = "DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR";
  else if (stage3 && sourceTcr && tcrSupport) decision = "DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY";
  else if (stage3) decision = "DL_PRIMARY_SURVIVOR_MODEL_ONLY";
  else if (tcrRescue) decision = "TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL";
  else decision = "DL_HOLD_MODEL_SIGNAL_NOT_ROBUST";
  return {{supp,stage1,stage2,stage3,sourceTcr,pairedTcr,tcrSupport,tcrRescue,bakerPass,mdPass,decision,dropoutUnstable}};
}}
function isPredPositive(r,c) {{
  if (callRule === "non_culled") return !c.decision.startsWith("DL_CULL") && c.decision !== "DL_TCR_DISCORDANT_AUDIT_NOT_POSITIVE";
  if (callRule === "robust_dl") return c.stage3;
  if (callRule === "tcr_supported") return c.stage3 && c.sourceTcr && c.tcrSupport;
  if (callRule === "structure_md_supported") return c.stage3 && c.sourceTcr && c.tcrSupport && c.bakerPass && c.mdPass;
  if (callRule === "wetlab_shortlist") return c.decision === "DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR" && !!r.wt_or_decoy_ready;
  return ["DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR","TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL"].includes(c.decision);
}}
function labelBadge(r) {{
  return num(r.label_binary) === 1 ? `<span class="badge label-pos">actual positive</span>` : `<span class="badge label-neg">actual negative</span>`;
}}
function predBadge(pred) {{
  return pred ? `<span class="badge pred-pos">model positive</span>` : `<span class="badge pred-neg">model negative</span>`;
}}
function outcomeName(r,c) {{
  const actual = num(r.label_binary) === 1;
  const pred = isPredPositive(r,c);
  if (pred && actual) return "TP";
  if (!pred && !actual) return "TN";
  if (pred && !actual) return "FP";
  return "FN";
}}
function outcomeBadge(outcome) {{
  const label = {{TP:"true positive", TN:"true negative", FP:"false positive", FN:"false negative"}}[outcome];
  return `<span class="badge out-${{outcome.toLowerCase()}}">${{label}}</span>`;
}}
const stageDefs = [
  ["00_all_candidates", r=>true],
  ["01_supported_peptide_length", r=>classify(r).supp],
  ["02_main_dl_broad_pass", r=>classify(r).stage1],
  ["03_bayesian_dropout_pass", r=>classify(r).stage2],
  ["04_robust_main_dl_pass", r=>classify(r).stage3],
  ["05_source_compatible_tcr_after_dl", r=>classify(r).stage3 && classify(r).sourceTcr],
  ["06_paired_tcr_after_dl", r=>classify(r).stage3 && classify(r).pairedTcr],
  ["07_tcr_branch_supported_after_dl", r=>classify(r).stage3 && classify(r).sourceTcr && classify(r).tcrSupport],
  ["08_baker_structure_supported", r=>classify(r).stage3 && classify(r).sourceTcr && classify(r).tcrSupport && classify(r).bakerPass],
  ["09_md_supported_or_complete", r=>classify(r).stage3 && classify(r).sourceTcr && classify(r).tcrSupport && classify(r).mdPass],
  ["10_tcr_rescue_reviews", r=>classify(r).decision==="TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL"],
  ["11_md_escalation_after_cheap_gates", r=>["DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR","TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL"].includes(classify(r).decision)],
  ["12_wetlab_shortlist_after_dl", r=>classify(r).decision==="DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR" && !!r.wt_or_decoy_ready],
];
function countsFor(maskFn) {{
  let n=0,pos=0,neg=0,unk=0;
  for (const r of rows) if (maskFn(r)) {{ n++; const lab=num(r.label_binary); if(lab===1) pos++; else if(lab===0) neg++; else unk++; }}
  return {{n,pos,neg,unk}};
}}
function render() {{
  const stageCounts = stageDefs.map(([name,fn]) => [name, countsFor(fn)]);
  const classified = rows.map(r => [r, classify(r)]);
  let tp=0, tn=0, fp=0, fn=0;
  for (const [r,c] of classified) {{
    const out = outcomeName(r,c);
    if (out === "TP") tp++;
    else if (out === "TN") tn++;
    else if (out === "FP") fp++;
    else fn++;
  }}
  const precision = tp + fp ? tp / (tp + fp) : 0;
  const recall = tp + fn ? tp / (tp + fn) : 0;
  const specificity = tn + fp ? tn / (tn + fp) : 0;
  const fpr = fp + tn ? fp / (fp + tn) : 0;
  const last = stageCounts[stageCounts.length-1][1];
  const md = stageCounts[stageCounts.length-2][1];
  const all = stageCounts[0][1];
  document.getElementById("stats").innerHTML = `
    <div class="stat"><b>${{all.n}}</b><span>start</span></div>
    <div class="stat"><b>${{stageCounts[2][1].n}}</b><span>main DL pass</span></div>
    <div class="stat"><b>${{stageCounts[4][1].n}}</b><span>robust DL</span></div>
    <div class="stat"><b>${{md.n}}</b><span>MD escalation</span></div>
    <div class="stat"><b>${{last.n}}</b><span>wetlab shortlist</span></div>`;
  const maxN = all.n || 1;
  document.getElementById("bars").innerHTML = stageCounts.map(([name,c]) => {{
    const pw = 100*c.pos/maxN, nw = 100*c.neg/maxN, uw = 100*c.unk/maxN;
    return `<div class="barrow"><div>${{name}}</div><div class="bar"><div class="pos" style="width:${{pw}}%"></div><div class="neg" style="width:${{nw}}%"></div><div class="unk" style="width:${{uw}}%"></div></div><div>${{c.n}} (+${{c.pos}}/-${{c.neg}})</div></div>`;
  }}).join("");
  document.getElementById("confusion").innerHTML = `
    <div class="cell tp"><span>true positive</span><b>${{tp}}</b></div>
    <div class="cell tn"><span>true negative</span><b>${{tn}}</b></div>
    <div class="cell fp"><span>false positive</span><b>${{fp}}</b></div>
    <div class="cell fn"><span>false negative</span><b>${{fn}}</b></div>`;
  document.getElementById("metrics").textContent = `precision=${{precision.toFixed(3)}} · recall=${{recall.toFixed(3)}} · specificity=${{specificity.toFixed(3)}} · FPR=${{fpr.toFixed(3)}}`;
  const ranked = classified
    .filter(x => !x[1].decision.startsWith("DL_CULL"))
    .sort((a,b) => Number(isPredPositive(b[0],b[1])) - Number(isPredPositive(a[0],a[1])) || num(b[0].main_dl_score)+num(b[0].tcr_augmented_score_mean)-num(a[0].main_dl_score)-num(a[0].tcr_augmented_score_mean))
    .slice(0,40);
  document.querySelector("#top tbody").innerHTML = ranked.map(([r,c]) => {{
    const pred = isPredPositive(r,c);
    const out = outcomeName(r,c);
    return `<tr><td>${{r.row_id}}</td><td>${{r.peptide}}</td><td>${{r.hla_4digit}}</td><td>${{labelBadge(r)}}</td><td>${{predBadge(pred)}}</td><td>${{outcomeBadge(out)}}</td><td>${{c.decision}}</td><td>${{num(r.main_dl_score).toFixed(3)}}</td><td>${{num(r.bayes_mean).toFixed(3)}}</td><td>${{num(r.tcr_augmented_score_mean).toFixed(3)}}</td><td>${{num(r.baker_structural_score).toFixed(3)}}</td><td>${{r.md_label||""}}</td></tr>`;
  }}).join("");
}}
render();
</script></body></html>"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mod = load_dl_module()
    df = add_baker_scores(mod.score_and_gate(mod.load_inputs()))
    stage_counts = stage_label_counts(df)
    decision_counts = decision_label_counts(df)
    stage_counts.to_csv(OUT / "dl_first_stage_label_counts.tsv", sep="\t", index=False)
    decision_counts.to_csv(OUT / "dl_first_decision_label_counts.tsv", sep="\t", index=False)
    payload_rows = js_rows(df)
    (OUT / "dl_first_slider_payload.json").write_text(json.dumps(payload_rows, ensure_ascii=False, indent=2))
    html_text = build_html(payload_rows, stage_counts, decision_counts)
    PAGE.write_text(html_text)
    WEB.mkdir(parents=True, exist_ok=True)
    shutil.copy2(PAGE, WEB_PAGE)
    print(f"[dl-funnel-ui] wrote {PAGE}")
    print(f"[dl-funnel-ui] deployed {WEB_PAGE}")
    print(stage_counts.to_string(index=False))


if __name__ == "__main__":
    main()
