#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import html
from pathlib import Path

import pandas as pd

from neoimmune_common import REPO_ROOT, ensure_run_dir, safe_read_table


def fmt(x, digits=3):
    try:
        if pd.isna(x):
            return "NA"
        return f"{float(x):.{digits}f}"
    except Exception:
        return html.escape(str(x))


def table_html(df: pd.DataFrame, cols: list[str], n: int = 10) -> str:
    if df.empty:
        return "<p class='muted'>No rows available.</p>"
    use = [c for c in cols if c in df.columns]
    d = df[use].head(n).copy()
    head = "".join(f"<th>{html.escape(c)}</th>" for c in use)
    rows = []
    for _, r in d.iterrows():
        tds = "".join(f"<td>{fmt(r[c]) if isinstance(r[c], float) else html.escape(str(r[c]))}</td>" for c in use)
        rows.append(f"<tr>{tds}</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def data_uri(path: Path) -> str:
    try:
        raw = path.read_bytes()
        return "data:image/png;base64," + base64.b64encode(raw).decode("ascii")
    except Exception:
        return ""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--html", default=str(REPO_ROOT / "project/papers_hub_2026_05_04/neoimmune_stack_clean_neo_plus_plus.html"))
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    ext = safe_read_table(outdir / "external_scores" / "all_external_scores.tsv.gz")
    local = safe_read_table(outdir / "local_scores" / "all_local_scores.tsv.gz")
    clean = safe_read_table(outdir / "metrics" / "leaderboard_clean_track_strict_no_overlap.tsv")
    prod = safe_read_table(outdir / "metrics" / "leaderboard_production_track_strict_no_overlap.tsv")
    leak = safe_read_table(outdir / "metrics" / "leakage_summary.tsv")
    top20 = safe_read_table(outdir / "predictions" / "patient_top20_candidates.tsv")
    alg_matrix = safe_read_table(outdir / "reports" / "algorithm_integration_matrix.tsv") if (outdir / "reports" / "algorithm_integration_matrix.tsv").exists() else pd.DataFrame()
    data_req = safe_read_table(outdir / "reports" / "hospital_data_request_sheet.tsv") if (outdir / "reports" / "hospital_data_request_sheet.tsv").exists() else pd.DataFrame()
    fig_manifest = safe_read_table(outdir / "reports" / "paper_figure_table_manifest_kr.tsv") if (outdir / "reports" / "paper_figure_table_manifest_kr.tsv").exists() else pd.DataFrame()
    public_inv = safe_read_table(outdir / "reports" / "public_patient_level_data_inventory.tsv") if (outdir / "reports" / "public_patient_level_data_inventory.tsv").exists() else pd.DataFrame()
    rescue_summary = safe_read_table(outdir / "metrics" / "blockbuster_rescue_summary.tsv") if (outdir / "metrics" / "blockbuster_rescue_summary.tsv").exists() else pd.DataFrame()
    top34_metrics = safe_read_table(outdir / "metrics" / "patient_top34_moderna_like_metrics.tsv") if (outdir / "metrics" / "patient_top34_moderna_like_metrics.tsv").exists() else pd.DataFrame()
    ga_metrics = safe_read_table(outdir / "metrics" / "ga_rl_algorithm_metrics.tsv") if (outdir / "metrics" / "ga_rl_algorithm_metrics.tsv").exists() else pd.DataFrame()
    ga_weights = safe_read_table(outdir / "reports" / "kg_ga_feature_weights.tsv") if (outdir / "reports" / "kg_ga_feature_weights.tsv").exists() else pd.DataFrame()
    nmi_metrics = safe_read_table(outdir / "metrics" / "nmi_clean_method_leaderboard.tsv") if (outdir / "metrics" / "nmi_clean_method_leaderboard.tsv").exists() else pd.DataFrame()
    imneo_metrics = safe_read_table(outdir / "metrics" / "imneo_public_reconstruction_leaderboard.tsv") if (outdir / "metrics" / "imneo_public_reconstruction_leaderboard.tsv").exists() else pd.DataFrame()
    imneo_sources = safe_read_table(outdir / "reports" / "imneo_public_source_inventory.tsv") if (outdir / "reports" / "imneo_public_source_inventory.tsv").exists() else pd.DataFrame()
    fig = outdir / "figures"
    fig_system = data_uri(fig / "neoimmune_stack_system_map.png")
    fig_clean = data_uri(fig / "strict_clean_leaderboard.png")
    fig_prod = data_uri(fig / "strict_production_leaderboard.png")
    fig_gap = data_uri(fig / "patient_data_gap.png")
    fig_claim = data_uri(fig / "claim_ladder.png")
    fig_rescue = data_uri(fig / "blockbuster_rescue_engine.png")
    fig_wetlab = data_uri(fig / "wetlab_bridge_to_publishable_evidence.png")
    fig_market = data_uri(fig / "business_market_wedge.png")
    fig_nmi = data_uri(fig / "nmi_clean_method_leaderboard.png")
    fig_imneo = data_uri(fig / "imneo_public_reconstruction_leaderboard.png")
    fig_imneo_sources = data_uri(fig / "imneo_public_source_coverage.png")
    fig_source_gap = data_uri(fig / "source_transfer_gap.png")
    fig_source_rescue = data_uri(fig / "source_rescue_queue.png")
    fig_collab = data_uri(fig / "collaborator_positive_handoff.png")
    fig_deck11 = data_uri(fig / "deck_closeout_slide11.png")
    fig_deck12 = data_uri(fig / "deck_closeout_slide12.png")
    real_patients = int(canon["patient_id"].fillna("unknown_patient").astype(str).ne("unknown_patient").sum()) if "patient_id" in canon else 0
    real_patient_unique = int(canon.loc[canon["patient_id"].fillna("unknown_patient").astype(str).ne("unknown_patient"), "patient_id"].astype(str).nunique()) if "patient_id" in canon else 0
    imneo_head = imneo_metrics[imneo_metrics["subset"].astype(str).eq("imNEO_core_head_to_head_common")].copy() if not imneo_metrics.empty and "subset" in imneo_metrics else pd.DataFrame()
    best_clean = clean.iloc[0] if not clean.empty else {}
    best_prod = prod.iloc[0] if not prod.empty else {}
    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>NeoImmune-Stack / CLEAN-Neo++</title>
<style>
:root {{
  --bg:#080b12; --panel:#111827; --panel2:#0d1320; --text:#eef5ff; --muted:#9fb0c7;
  --cyan:#26d9ff; --violet:#a78bfa; --gold:#f6c768; --red:#fb7185; --green:#5ee6a8; --line:#263348;
}}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:radial-gradient(circle at 30% 0%, #16223b 0, #080b12 44%); color:var(--text); font:15px/1.55 Inter, system-ui, -apple-system, Segoe UI, sans-serif; }}
a {{ color:var(--cyan); text-decoration:none; }}
.wrap {{ max-width:1320px; margin:0 auto; padding:28px; }}
.hero {{ min-height:420px; display:grid; grid-template-columns:1.25fr .75fr; gap:28px; align-items:center; border-bottom:1px solid var(--line); }}
.kicker {{ color:var(--cyan); text-transform:uppercase; letter-spacing:.12em; font-weight:700; font-size:12px; }}
h1 {{ font-size:56px; line-height:1.0; margin:12px 0 18px; letter-spacing:0; }}
h2 {{ font-size:24px; margin:0 0 14px; }}
h3 {{ font-size:18px; margin:0 0 10px; }}
.lead {{ font-size:21px; color:#dbeafe; max-width:920px; }}
.verdict {{ border:1px solid var(--line); background:linear-gradient(160deg,#111827,#111827 60%,#1b1730); padding:22px; border-radius:8px; }}
.verdict b {{ color:var(--gold); }}
.stats {{ display:grid; grid-template-columns:repeat(4,1fr); gap:12px; margin-top:24px; }}
.stat {{ background:#0d1320; border:1px solid var(--line); padding:16px; border-radius:8px; min-height:94px; }}
.stat .num {{ font-size:30px; font-weight:800; color:var(--cyan); }}
.stat .lab {{ color:var(--muted); font-size:12px; text-transform:uppercase; letter-spacing:.08em; }}
.grid {{ display:grid; grid-template-columns:1fr 1fr; gap:18px; margin:22px 0; }}
.panel {{ background:rgba(17,24,39,.92); border:1px solid var(--line); border-radius:8px; padding:20px; overflow:auto; }}
.wide {{ grid-column:1 / -1; }}
.muted {{ color:var(--muted); }}
.danger {{ color:var(--red); font-weight:700; }}
.ok {{ color:var(--green); font-weight:700; }}
.gold {{ color:var(--gold); font-weight:700; }}
table {{ width:100%; border-collapse:collapse; font-size:13px; }}
th,td {{ padding:8px 9px; border-bottom:1px solid #263348; text-align:left; vertical-align:top; }}
th {{ color:#b6c7e1; font-size:11px; text-transform:uppercase; letter-spacing:.06em; }}
.flow {{ display:grid; grid-template-columns:repeat(5,1fr); gap:12px; align-items:stretch; }}
.node {{ border:1px solid var(--line); background:#0d1320; border-radius:8px; padding:15px; min-height:132px; position:relative; }}
.node:after {{ content:"→"; position:absolute; right:-14px; top:45%; color:var(--cyan); font-size:24px; }}
.node:last-child:after {{ content:""; }}
.node .t {{ color:var(--cyan); font-weight:800; margin-bottom:8px; }}
.pill {{ display:inline-block; padding:3px 8px; border:1px solid var(--line); border-radius:999px; margin:3px 4px 3px 0; color:#dbeafe; background:#0b1220; font-size:12px; }}
.bar {{ height:9px; border-radius:99px; background:#1f2937; overflow:hidden; }}
.fill {{ height:100%; background:linear-gradient(90deg,var(--cyan),var(--violet)); }}
.boundary {{ display:grid; grid-template-columns:1fr 1fr; gap:14px; }}
.boundary ul {{ margin:8px 0 0 18px; padding:0; }}
.footer {{ color:var(--muted); border-top:1px solid var(--line); padding:22px 0 40px; margin-top:24px; }}
img.figure {{ width:100%; border:1px solid var(--line); border-radius:8px; display:block; background:#080b12; }}
.figcap {{ color:var(--muted); font-size:13px; margin-top:8px; }}
@media (max-width:900px) {{ .hero,.grid,.flow,.boundary {{ grid-template-columns:1fr; }} h1 {{ font-size:40px; }} .stats {{ grid-template-columns:1fr 1fr; }} .node:after {{ content:""; }} }}
</style>
</head>
<body>
<div class="wrap">
  <section class="hero">
    <div>
      <div class="kicker">NeoImmune-Stack / CLEAN-Neo++</div>
      <h1>Leakage-aware cancer vaccine candidate ranking layer</h1>
      <p class="lead">내 알고리즘과 공개 HLA/pMHC/immunogenicity 도구를 합친 통합 레이어입니다. 핵심은 clinical efficacy claim이 아니라, patient-context-aware top-N ranking과 leakage-aware immunogenicity strategy입니다.</p>
      <div class="stats">
        <div class="stat"><div class="num">{len(canon):,}</div><div class="lab">canonical candidates</div></div>
        <div class="stat"><div class="num">{len(local):,}</div><div class="lab">local model scores</div></div>
        <div class="stat"><div class="num">{len(ext):,}</div><div class="lab">external score rows</div></div>
        <div class="stat"><div class="num">{real_patient_unique:,}</div><div class="lab">real patients</div></div>
      </div>
    </div>
    <div class="verdict">
      <h2>Brutal Verdict</h2>
      <p><b>대박 가능성:</b> product/integration layer로는 강합니다. 기존 BAR-Neo, Wave8/TCR, Structure, ESM2, quantum branch와 public comparators가 한 구조로 묶였습니다.</p>
      <p><span class="ok">업그레이드:</span> 공개 통합 테이블에서 {real_patient_unique:,}명 / {real_patients:,} patient-ID rows가 살아났습니다. Retrospective public patient-level ranking scaffold는 이제 주장 가능합니다.</p>
      <p><span class="danger">아직 금지:</span> prospective hospital utility, clinical vaccine efficacy, MS 없는 presentation-confirmed claim, T-cell assay 없는 immunogenicity-confirmed claim.</p>
    </div>
  </section>

  <section class="grid">
    <div class="panel wide">
      <h2>Master System Figure</h2>
      <img class="figure" src="{fig_system}" alt="NeoImmune-Stack system map">
      <div class="figcap">Role: shows how local algorithms and public predictors are separated into clean science and production tracks. Claim boundary: architecture and evaluation scaffold, not clinical validation.</div>
    </div>

    <div class="panel wide">
      <h2>Blockbuster Rescue Engine</h2>
      <img class="figure" src="{fig_rescue}" alt="Blockbuster rescue engine">
      <div class="figcap">Role: shows why local immunogenicity branches matter beyond public binding/presentation predictors. Claim boundary: hypothesis-generation audit, not vaccine efficacy.</div>
      {table_html(rescue_summary, ["metric","value"], 10)}
    </div>

    <div class="panel wide">
      <h2>Moderna/V940 Public Translation</h2>
      <p class="muted">Public V940 evidence supports a top34 patient-specific neoantigen cap, adjuvant PD-1 combination framing, and T-cell response endpoints. It does not expose Moderna's patient-level raw neoantigen tables or proprietary selection algorithm.</p>
      {table_html(top34_metrics, ["track","score_col","N","patient_hit_rate@34","patient_recall@34","patients_evaluated","rationale"], 6)}
    </div>

    <div class="panel wide">
      <h2>GA/RL-Discovered Controller Applied</h2>
      <p class="muted">KG-GA evolutionary architecture search is now integrated as a production/experiment-priority branch. It is not clean-track allowed because the evolved controller includes public-predictor and product-value components.</p>
      {table_html(ga_metrics, ["algorithm","score_column","n","positives","AUPRC","AUROC","Precision@34","Recall@34","patient_hit_rate@34","patient_recall@34","clean_track_allowed"], 10)}
      <h3>Top evolved weights</h3>
      {table_html(ga_weights, ["feature","concept","role","weight","transform"], 12)}
    </div>

    <div class="panel wide">
      <h2>NMI Clean Method Track</h2>
      <p class="muted">NMI = Neoantigen Multimodal Immunogenicity. This is the method-paper branch: no NetMHCpan, MHCflurry, BigMHC, PRIME, BAR-Neo, KG-GA, or product-impact scores are used as NMI training features. Public predictors are comparators only.</p>
      <img class="figure" src="{fig_nmi}" alt="NMI clean method leaderboard">
      <div class="figcap">Role: shows whether the clean in-house multimodal branch beats frozen public predictors on the same NMI-eligible rows. Claim boundary: retrospective public labels; not clinical vaccine efficacy.</div>
      {table_html(nmi_metrics, ["model","split","clean_track_allowed","n","positives","AUPRC","AUROC","Precision@34","Recall@34","patient_hit_rate@34","patient_recall@34"], 16)}
    </div>

    <div class="panel wide">
      <h2>imNEO Public-Source Reconstruction</h2>
      <p class="muted">This section compares only the public-source space visible from CG Invites/imNEO patent/presentation context. It is not the proprietary imNEO peptide list and not a reverse-engineered company model. The confident comparison is NMI locked/common rows vs frozen public comparators.</p>
      <img class="figure" src="{fig_imneo}" alt="imNEO public reconstruction leaderboard">
      <div class="figcap">Role: competitor-facing but claim-safe comparison. Claim boundary: public-source reconstruction; retrospective immunogenicity labels; patient-level top-N is not claimable here because rows collapse to unknown/synthetic patient IDs.</div>
      {table_html(imneo_head, ["model","track","clean_track_allowed","n","positives","real_patients","patient_level_claimable","AUPRC","AUROC","Precision@34","Recall@34"], 14)}
      <h3>Public-source coverage</h3>
      <img class="figure" src="{fig_imneo_sources}" alt="imNEO public source coverage">
      {table_html(imneo_sources, ["source_name","patent_role","status","present_rows","positive_labels","nmi_locked_rows","public_comparator_rows","claim_boundary"], 10)}
    </div>

    <div class="panel wide">
      <h2>Source Transfer Gap</h2>
      <p class="muted">This is the direct answer to where the clean algorithm already beats frozen public comparators and where coverage is still thin. It is the next reviewer-facing slide after the imNEO public-source reconstruction.</p>
      <img class="figure" src="{fig_source_gap}" alt="source transfer gap">
      <div class="figcap">Role: source-family gap audit. Claim boundary: source-family level only; not a clinical validation or patient-level benchmark.</div>
      {table_html(safe_read_table(outdir / "metrics" / "source_transfer_gap_summary.tsv") if (outdir / "metrics" / "source_transfer_gap_summary.tsv").exists() else pd.DataFrame(), ["source_family","rows","positives","nmi_best_model","nmi_best_auprc","public_best_model","public_best_auprc","delta_nmi_minus_public","patient_claimable"], 12)}
    </div>

    <div class="panel wide">
      <h2>Source Rescue Queue</h2>
      <p class="muted">This is the practical handoff list for TESLA, dbPepNeo, CEDAR, and IMPROVE families. Positive high-priority rows are the ones to bring to collaborators first; negative rows remain failure controls.</p>
      <img class="figure" src="{fig_source_rescue}" alt="source rescue queue">
      <div class="figcap">Role: candidate-level rescue queue. Claim boundary: retrospective prioritization only; not a validated vaccine target set.</div>
      {table_html(safe_read_table(outdir / "metrics" / "source_rescue_queue_summary.tsv") if (outdir / "metrics" / "source_rescue_queue_summary.tsv").exists() else pd.DataFrame(), ["dataset_source","rows","positives","top_priority","median_priority","mean_clean","mean_public_best","mean_delta"], 12)}
    </div>

    <div class="panel wide">
      <h2>Collaborator Handoff</h2>
      <p class="muted">This is the minimal wet-lab-facing shortlist. The goal is not to maximize the table length; it is to make the first discussion with collaborators obvious.</p>
      <img class="figure" src="{fig_collab}" alt="collaborator positive handoff">
      <div class="figcap">Role: top positive rows by source family. Claim boundary: retrospective handoff only; no candidate is validated without orthogonal assay support.</div>
      {table_html(safe_read_table(outdir / "predictions" / "collaborator_positive_handoff.tsv") if (outdir / "predictions" / "collaborator_positive_handoff.tsv").exists() else pd.DataFrame(), ["candidate_id","dataset_source","peptide_mut","hla_allele","label_immunogenicity","public_best","clean_science_score","production_stack_score","local_rescue_delta_vs_public","model_disagreement_score","rescue_priority"], 18)}
      <h3>Failure controls</h3>
      {table_html(safe_read_table(outdir / "predictions" / "collaborator_control_handoff.tsv") if (outdir / "predictions" / "collaborator_control_handoff.tsv").exists() else pd.DataFrame(), ["candidate_id","dataset_source","peptide_mut","hla_allele","label_immunogenicity","public_best","clean_science_score","production_stack_score","local_rescue_delta_vs_public","model_disagreement_score","rescue_priority"], 12)}
    </div>

    <div class="panel wide">
      <h2>Deck Closeout</h2>
      <p class="muted">These are the final two slides for the collaborator deck. Slide 11 is the positive handoff. Slide 12 is the control and claim-boundary slide.</p>
      <div class="grid" style="grid-template-columns:1fr 1fr;">
        <div><img class="figure" src="{fig_deck11}" alt="deck closeout slide 11"><div class="figcap">Slide 11: hand off the shortlist.</div></div>
        <div><img class="figure" src="{fig_deck12}" alt="deck closeout slide 12"><div class="figcap">Slide 12: keep controls and claims visible.</div></div>
      </div>
      <p class="muted">Report: <code>{html.escape(str(outdir / "reports/DECK_CLOSEOUT.md"))}</code></p>
    </div>

    <div class="panel wide">
      <h2>System Flow</h2>
      <div class="flow">
        <div class="node"><div class="t">1. Candidate Sources</div><span class="pill">CLEAN-Neobench</span><span class="pill">cross_neo_v0</span><span class="pill">neoantigen hub</span><p class="muted">Presentation and immunogenicity labels are kept separate.</p></div>
        <div class="node"><div class="t">2. Local Algorithms</div><span class="pill">Structure_LR</span><span class="pill">Wave8_TCR</span><span class="pill">ESM2</span><span class="pill">Quantum</span><p class="muted">Allowed in clean science track.</p></div>
        <div class="node"><div class="t">3. Public Predictors</div><span class="pill">MHCflurry</span><span class="pill">BigMHC</span><span class="pill">PRIME</span><span class="pill">NetMHCpan</span><p class="muted">Frozen features only; production track only.</p></div>
        <div class="node"><div class="t">4. Leakage Audit</div><span class="pill">peptide</span><span class="pill">pHLA</span><span class="pill">study</span><span class="pill">patient</span><p class="muted">Strict no-existing-overlap views generated.</p></div>
        <div class="node"><div class="t">5. Top-N Handoff</div><span class="pill">top20</span><span class="pill">failure audit</span><span class="pill">wet-lab queue</span><p class="muted">Not validated vaccine candidates.</p></div>
      </div>
    </div>

    <div class="panel">
      <h2>Strict Clean Track</h2>
      <p class="muted">External public predictor scores excluded. Best row is artifact-level unless source-heldout/no-overlap is explicit.</p>
      <img class="figure" src="{fig_clean}" alt="Strict clean leaderboard">
      <div class="figcap">Role: emphasizes your local algorithms under the clean track. Claim boundary: strict no-existing-overlap view, not prospective patient validation.</div>
      {table_html(clean, ["model","split","n","positives","AUPRC","AUROC","Precision@20","Recall@20"], 10)}
    </div>

    <div class="panel">
      <h2>Strict Production Track</h2>
      <p class="muted">Frozen public scores allowed. This is practical prioritization, not a clean novel-model claim.</p>
      <img class="figure" src="{fig_prod}" alt="Strict production leaderboard">
      <div class="figcap">Role: shows practical stack/comparator behavior. Claim boundary: production ranking, not clean algorithm novelty.</div>
      {table_html(prod, ["model","split","n","positives","AUPRC","AUROC","Precision@20","Recall@20"], 10)}
    </div>

    <div class="panel">
      <h2>Leakage Controls</h2>
      {table_html(leak, ["leakage_control","n_flagged","denominator","fraction"], 8)}
    </div>

    <div class="panel">
      <h2>Top-20 Handoff Warning</h2>
      <img class="figure" src="{fig_gap}" alt="Patient data gap">
      <div class="figcap">Role: makes the next required hospital data obvious. Claim boundary: patient-level endpoint is not claimable from public scaffold alone.</div>
      <p class="ok">Public patient IDs are now recovered enough for retrospective top-N scaffold reporting. This is still not prospective hospital validation.</p>
      {table_html(top20, ["patient_id","candidate_id","dataset_source","peptide_mut","hla_allele","production_stack_score","clean_science_score","model_disagreement_score"], 8)}
    </div>

    <div class="panel wide">
      <h2>Wet-lab Bridge</h2>
      <img class="figure" src="{fig_wetlab}" alt="Wet-lab bridge to publishable evidence">
      <div class="figcap">Role: converts top-N ranking into assay-confirmed evidence. Claim boundary: presentation requires MS; immunogenicity requires T-cell assay.</div>
    </div>

    <div class="panel wide">
      <h2>Business Wedge</h2>
      <img class="figure" src="{fig_market}" alt="Business market wedge">
      <div class="figcap">Role: positions NeoImmune-Stack as a vendor-agnostic operating layer, not another single predictor.</div>
    </div>

    <div class="panel wide">
      <h2>Claim Boundary Matrix</h2>
      <img class="figure" src="{fig_claim}" alt="Claim ladder">
      <div class="figcap">Role: converts reviewer risk into a claim ladder. Claim boundary: separates what is strong now from what needs wet-lab or patient data.</div>
      <div class="boundary">
        <div>
          <h3 class="ok">Allowed</h3>
          <ul>
            <li>Leakage-aware integration of existing local algorithms and public frozen comparators.</li>
            <li>Clean-science track excludes NetMHCpan/MHCflurry/BigMHC/PRIME features.</li>
            <li>Production stack is a practical candidate prioritization layer.</li>
            <li>Wave8/TCR-self-similarity is a biologically important immunogenicity branch.</li>
          </ul>
        </div>
        <div>
          <h3 class="danger">Forbidden</h3>
          <ul>
            <li>No clinical vaccine efficacy claim.</li>
            <li>No antigen presentation claim without MS or equivalent evidence.</li>
            <li>No immunogenicity claim without T-cell assay labels.</li>
            <li>No patient-level endpoint claim until true patient IDs exist.</li>
          </ul>
        </div>
      </div>
    </div>

    <div class="panel wide">
      <h2>Algorithm Integration Matrix</h2>
      <p class="muted">This is the core answer to “내 알고리즘이랑 합치는거야?” Yes: local algorithms are the clean-science core; public predictors are frozen production comparators/features.</p>
      {table_html(alg_matrix, ["component","source_type","role","clean_track_feature_allowed","production_track_feature_allowed","claim_boundary"], 12)}
    </div>

    <div class="panel wide">
      <h2>Hospital Data Contract</h2>
      <p class="muted">This is the gate that turns the scaffold into a real patient-level top-N paper/business package.</p>
      {table_html(data_req, ["field","priority","why_needed","acceptable_format","risk_if_missing"], 15)}
    </div>

    <div class="panel wide">
      <h2>Open Patient-Level Data Hunt</h2>
      <p class="muted">There are open/public patient-level data sources. The strongest immediate local source is NEPdb; TESLA and IMPROVE are the next recovery targets for richer original patient-level fields.</p>
      {table_html(public_inv, ["dataset","public_status","patient_level","rows","unique_patients","best_use","source_url"], 10)}
    </div>

    <div class="panel wide">
      <h2>Paper Figure/Table Manifest</h2>
      <p class="muted">This turns the analysis package into a manuscript-ready visual story.</p>
      {table_html(fig_manifest, ["item","title","file","role_in_story","claim_boundary"], 12)}
    </div>

    <div class="panel wide">
      <h2>Outputs</h2>
      <p><span class="gold">Result root:</span> <code>{html.escape(str(outdir))}</code></p>
      <p><span class="gold">Strategy report:</span> <code>{html.escape(str(outdir / "reports/NEOIMMUNE_STACK_STRATEGY_REPORT.md"))}</code></p>
      <p><span class="gold">Impact memo:</span> <code>{html.escape(str(outdir / "reports/CLEAN_NEO_PLUS_PLUS_IMPACT_MEMO_KR.md"))}</code></p>
      <p><span class="gold">Reviewer defense:</span> <code>{html.escape(str(outdir / "reports/reviewer_attack_defense_table_kr.md"))}</code></p>
      <p><span class="gold">Hospital data request:</span> <code>{html.escape(str(outdir / "reports/hospital_data_request_sheet.tsv"))}</code></p>
      <p><span class="gold">Manuscript scaffold:</span> <code>{html.escape(str(outdir / "reports/manuscript_results_scaffold_kr.md"))}</code></p>
      <p><span class="gold">Business brief:</span> <code>{html.escape(str(outdir / "reports/business_product_brief_kr.md"))}</code></p>
      <p><span class="gold">Blockbuster memo:</span> <code>{html.escape(str(outdir / "reports/BLOCKBUSTER_UPGRADE_MEMO_KR.md"))}</code></p>
      <p><span class="gold">Public data hunt:</span> <code>{html.escape(str(outdir / "reports/public_patient_level_data_hunt.md"))}</code></p>
      <p><span class="gold">Moderna evidence memo:</span> <code>{html.escape(str(outdir / "reports/moderna_public_data_and_combo_evidence_kr.md"))}</code></p>
      <p><span class="gold">Moderna top34 candidates:</span> <code>{html.escape(str(outdir / "predictions/patient_top34_candidates_moderna_like.tsv"))}</code></p>
      <p><span class="gold">GA/RL applied report:</span> <code>{html.escape(str(outdir / "reports/ga_rl_algorithm_applied_report.md"))}</code></p>
      <p><span class="gold">GA/RL top34 candidates:</span> <code>{html.escape(str(outdir / "predictions/patient_top34_candidates_kg_ga.tsv"))}</code></p>
      <p><span class="gold">NMI method report:</span> <code>{html.escape(str(outdir / "reports/NMI_METHOD_PAPER_STRATEGY.md"))}</code></p>
      <p><span class="gold">NMI leaderboard:</span> <code>{html.escape(str(outdir / "metrics/nmi_clean_method_leaderboard.tsv"))}</code></p>
      <p><span class="gold">NMI top34 candidates:</span> <code>{html.escape(str(outdir / "predictions/patient_top34_candidates_nmi_clean.tsv"))}</code></p>
      <p><span class="gold">imNEO public reconstruction report:</span> <code>{html.escape(str(outdir / "reports/IMNEO_PUBLIC_RECONSTRUCTION_BENCHMARK.md"))}</code></p>
      <p><span class="gold">imNEO public reconstruction leaderboard:</span> <code>{html.escape(str(outdir / "metrics/imneo_public_reconstruction_leaderboard.tsv"))}</code></p>
      <p><span class="gold">imNEO public reconstruction ranked candidates:</span> <code>{html.escape(str(outdir / "predictions/imneo_public_reconstruction_ranked_candidates.tsv"))}</code></p>
      <p><span class="gold">Wet-lab queue scaffold:</span> <code>{html.escape(str(outdir / "predictions/wetlab_candidate_queue_scaffold.tsv"))}</code></p>
      <p><span class="gold">Canonical table:</span> <code>{html.escape(str(outdir / "data/canonical_candidates.tsv.gz"))}</code> and <code>canonical_candidates.parquet</code></p>
      <p><span class="gold">Runner:</span> <code>bash project/scripts/neoimmune_stack/run_neoimmune_stack.sh project/results/neoimmune_stack_2026_05_10</code></p>
    </div>
  </section>
  <div class="footer">Generated by NeoImmune-Stack integration scripts. Every number should be interpreted with the stated leakage and label boundaries.</div>
</div>
</body>
</html>
"""
    path = Path(args.html)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html_text, encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
