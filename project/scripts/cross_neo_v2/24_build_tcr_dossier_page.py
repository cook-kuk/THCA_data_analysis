#!/usr/bin/env python3
"""Build and deploy a one-page CROSS-Neo-TCR dossier for papers_hub."""

from __future__ import annotations

import html
import shutil
from pathlib import Path

import pandas as pd

from common import OUT, REPO


TCR_OUT = OUT / "tcr_extension"
HUB = REPO / "project/papers_hub_2026_05_04"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/cross_neo_tcr_extension"
WEB_ASSET_DIR = WEB / "assets/cross_neo_tcr_extension"
PAGE = HUB / "cross_neo_tcr_extension_dossier.html"
WEB_PAGE = WEB / "cross_neo_tcr_extension_dossier.html"
PACKET_NAME = "cross_neo_tcr_reviewer_packet_2026_05_09.zip"
OPENMM_PILOT_PACKAGE = "openmm_p0_10ns_pilot_package_2026_05_09.tar.gz"


def esc(x: object) -> str:
    return html.escape("" if pd.isna(x) else str(x))


def fmt_num(x: object, digits: int = 3) -> str:
    if pd.isna(x):
        return "NA"
    try:
        val = float(x)
    except Exception:
        return esc(x)
    if abs(val) >= 1000:
        return f"{val:,.0f}"
    return f"{val:.{digits}f}"


def table(df: pd.DataFrame, cols: list[str], max_rows: int = 20, num_cols: set[str] | None = None) -> str:
    num_cols = num_cols or set()
    rows = ["<table class='t'><thead><tr>"]
    for c in cols:
        rows.append(f"<th>{esc(c.replace('_', ' '))}</th>")
    rows.append("</tr></thead><tbody>")
    for _, r in df.head(max_rows).iterrows():
        rows.append("<tr>")
        for c in cols:
            cls = " class='num'" if c in num_cols else ""
            val = fmt_num(r[c], 4) if c in num_cols else esc(r[c])
            rows.append(f"<td{cls}>{val}</td>")
        rows.append("</tr>")
    rows.append("</tbody></table>")
    return "\n".join(rows)


def copy_figures() -> list[str]:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    WEB_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    names = []
    for p in sorted((TCR_OUT / "figures").glob("fig_tcr*.png")):
        dest = ASSET_DIR / p.name
        shutil.copy2(p, dest)
        shutil.copy2(p, WEB_ASSET_DIR / p.name)
        names.append(p.stem)
    return names


def main() -> None:
    fig_names = copy_figures()
    counts = pd.read_csv(TCR_OUT / "tcr_registry_source_counts.tsv", sep="\t")
    linkage = pd.read_csv(TCR_OUT / "tcr_neo_linkage_summary.tsv", sep="\t")
    metrics = pd.read_csv(TCR_OUT / "model_discovery/external_tcr_expert_benchmark/external_tcr_expert_benchmark_metrics.tsv", sep="\t")
    source_metrics = pd.read_csv(TCR_OUT / "model_discovery/external_tcr_expert_benchmark/external_tcr_expert_source_heldout_metrics.tsv", sep="\t")
    wetlab = pd.read_csv(TCR_OUT / "model_discovery/wetlab_external_tcr_expert_scores/wetlab_candidates_external_tcr_expert_ranked.tsv", sep="\t")
    md_queue = pd.read_csv(TCR_OUT / "md_escalation/md_escalation_queue_top20.tsv", sep="\t")
    p0_pilots = pd.read_csv(TCR_OUT / "md_escalation/p0_structures/p0_md_pilot_ready_complexes.tsv", sep="\t")
    min_qc = pd.read_csv(TCR_OUT / "md_escalation/p0_structures/openmm_minimized/openmm_p0_minimization_qc.tsv", sep="\t")
    smoke_qc = pd.read_csv(TCR_OUT / "md_escalation/p0_structures/openmm_pilot_10ns_package/p0_openmm_smoke_qc.tsv", sep="\t")
    model_matrix = pd.read_csv(TCR_OUT / "model_discovery/tcr_model_discovery_matrix.tsv", sep="\t")
    jobs = pd.read_csv(TCR_OUT / "structure_jobs/tcr_structure_job_summary.tsv", sep="\t")

    pooled = metrics[metrics["panel"].eq("pooled")].copy().sort_values("auprc", ascending=False)
    total_rows = int(counts["n"].sum())
    paired = int(counts["paired_tcr"].sum())
    any_tcr = int(counts["any_tcr"].sum())
    peptide_hla = int(counts["peptide_hla"].sum())
    structures = int(counts["structures"].sum())
    exact = int(linkage.loc[linkage["tcr_link_category"].eq("exact_tcr_pmhc_match"), "n_neo_rows"].sum())
    overlap = int(linkage.loc[~linkage["tcr_link_category"].eq("no_tcr_match"), "n_neo_rows"].sum())
    no_match = int(linkage.loc[linkage["tcr_link_category"].eq("no_tcr_match"), "n_neo_rows"].sum())
    ensemble = pooled[pooled["model"].eq("mean_pMTnet_TEPCAM")].iloc[0]
    pmtnet = pooled[pooled["model"].eq("pMTnet")].iloc[0]
    top_wetlab = wetlab.iloc[0]
    top_md = md_queue.iloc[0]
    p0_md = int(md_queue["md_tier"].eq("P0_MD_TCR_pMHC").sum())
    p0_pilot_ready = int(len(p0_pilots))
    openmm_min_ok = int(min_qc["status"].eq("ok").sum())
    smoke_ok = int(smoke_qc["status"].eq("ok").sum())

    source_cols = ["source_dataset", "n", "paired_tcr", "any_tcr", "peptide_hla", "cancer_context", "pathogen_context", "structures"]
    linkage_cols = ["tcr_link_category", "n_neo_rows", "mean_tcr_evidence_count", "paired_tcr_evidence_count", "cancer_context_evidence_count"]
    perf_cols = ["model", "n", "auprc", "auroc", "positive_mean", "negative_mean"]
    source_perf_cols = ["heldout_source", "n", "auprc", "auroc", "positive_mean", "negative_mean"]
    wet_cols = ["peptide", "hla_4digit", "external_tcr_pairs_total", "external_tcr_expert_mean", "external_tcr_expert_max", "wetlab_priority_score_external_adjusted"]
    md_cols = ["peptide", "hla_4digit", "md_tier", "md_escalation_score", "prediction_fragility_score", "recommended_md_protocol"]
    pilot_cols = ["target_peptide", "hla_4digit", "pdb_id", "peptide_chain", "selected_chains", "contacts_tcr_peptide", "pilot_status"]
    min_cols = ["target_peptide", "hla_4digit", "pdb_id", "status", "atom_count_after_repair", "energy_initial_kj_mol", "energy_minimized_kj_mol"]
    smoke_cols = ["smoke_id", "status", "steps", "peptide_rmsd_final_nm", "peptide_mhc_contacts_final", "peptide_tcr_contacts_final"]
    model_cols = ["model", "priority", "input", "local_status", "claim_use", "blocker"]
    source_num_cols = set(source_cols) - {"source_dataset"}
    linkage_num_cols = set(linkage_cols) - {"tcr_link_category"}
    perf_num_cols = {"n", "auprc", "auroc", "positive_mean", "negative_mean"}
    wet_num_cols = {
        "external_tcr_pairs_total",
        "external_tcr_expert_mean",
        "external_tcr_expert_max",
        "wetlab_priority_score_external_adjusted",
    }
    md_num_cols = {"md_escalation_score", "prediction_fragility_score"}
    pilot_num_cols = {"contacts_tcr_peptide"}
    min_num_cols = {"atom_count_after_repair", "energy_initial_kj_mol", "energy_minimized_kj_mol"}
    smoke_num_cols = {"steps", "peptide_rmsd_final_nm", "peptide_mhc_contacts_final", "peptide_tcr_contacts_final"}
    job_num_cols = {c for c in jobs.columns if c != "priority"}

    fig_cards = []
    fig_titles = {
        "fig_tcr1_recognition_gap": "Recognition Gap",
        "fig_tcr2_architecture": "Architecture",
        "fig_tcr3_structure_pipeline": "Structure Pipeline",
        "fig_tcr4_mutant_wt_interface_delta": "Mutant-WT Delta",
        "fig_tcr5_tcr_available_subset_performance": "Performance",
        "fig_tcr6_case_studies": "Case Studies",
        "fig_tcr7_claim_boundary": "Claim Boundary",
    }
    for name in fig_names:
        fig_cards.append(
            f"""
            <article class="fig-card">
              <img src="assets/cross_neo_tcr_extension/{name}.png" alt="{esc(fig_titles.get(name, name))}">
              <h4>{esc(fig_titles.get(name, name))}</h4>
              <a href="assets/cross_neo_tcr_extension/{name}.png">PNG</a>
            </article>
            """
        )

    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo-TCR Extension Dossier</title>
<style>
:root{{--bg:#07111f;--panel:#0d1a2b;--panel2:#101f33;--ink:#eef5ff;--muted:#9dafc4;--line:#22324a;--gold:#f2c46d;--teal:#38d6b1;--blue:#70a7ff;--red:#ff7d70;--green:#67d59a;}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font-family:Inter,"Noto Sans KR",-apple-system,BlinkMacSystemFont,sans-serif;font-size:14px;line-height:1.55}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7e3}}
code,.mono{{font-family:"JetBrains Mono","SF Mono",monospace;font-size:12px}}
.hero{{padding:70px 34px 44px;background:linear-gradient(135deg,#101f33,#07111f 70%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1320px;margin:0 auto}}
.kicker{{font-family:"JetBrains Mono",monospace;color:var(--gold);letter-spacing:.28em;font-size:11px;text-transform:uppercase;font-weight:700}}
h1{{font-size:72px;line-height:.95;margin:18px 0 14px}}
.lead{{max-width:980px;color:#d7e2f0;font-size:18px;font-family:"Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(8,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.045);border:1px solid rgba(242,196,109,.22);border-radius:10px;padding:13px}}
.stat b{{display:block;color:var(--gold);font-size:26px;font-family:"Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font-family:"JetBrains Mono",monospace;font-size:9px;text-transform:uppercase;letter-spacing:.08em;margin-top:5px}}
.crumbs{{margin-top:22px;color:var(--muted);font-family:"JetBrains Mono",monospace;font-size:11px}}
.wrap{{max-width:1320px;margin:0 auto;display:grid;grid-template-columns:240px 1fr;gap:34px;padding:0 30px 80px}}
.toc{{position:sticky;top:0;align-self:start;max-height:100vh;overflow:auto;padding:28px 0;border-right:1px solid var(--line)}}
.toc h4{{font-family:"JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;font-size:10px;text-transform:uppercase}}
.toc a{{display:block;color:#cfdbeb;padding:5px 4px;font-size:12px}}
main{{min-width:0;padding-top:26px}}
section{{border-bottom:1px solid var(--line);padding:30px 0}}
h2{{font-size:36px;margin:0 0 8px}} h2 .num{{font-family:"JetBrains Mono",monospace;color:var(--gold);font-size:13px;margin-right:12px;letter-spacing:.16em}}
.sub{{font-family:"JetBrains Mono",monospace;color:var(--muted);font-size:11px;letter-spacing:.08em;text-transform:uppercase;margin-bottom:14px}}
.grid{{display:grid;gap:14px}} .grid2{{grid-template-columns:1fr 1fr}} .grid3{{grid-template-columns:repeat(3,1fr)}} .grid4{{grid-template-columns:repeat(4,1fr)}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:15px}}
.card h3,.card h4{{margin:0 0 8px;color:var(--gold)}} .card .big{{font-size:34px;color:#fff7e3;font-family:"Cormorant Garamond",serif;font-weight:700}}
.callout{{border-left:3px solid var(--gold);background:#101c2e;border-top:1px solid var(--line);border-right:1px solid var(--line);border-bottom:1px solid var(--line);padding:14px 16px;margin:14px 0;border-radius:0 8px 8px 0}}
.good{{border-left-color:var(--green)}} .warn{{border-left-color:var(--red)}}
table.t{{border-collapse:collapse;width:100%;font-size:12px;margin:10px 0 16px}} .t th,.t td{{border:1px solid var(--line);padding:7px 8px;vertical-align:top}} .t th{{background:#14243a;color:var(--gold);font-family:"JetBrains Mono",monospace;font-size:10px;text-transform:uppercase;letter-spacing:.04em}} .t tr:nth-child(even) td{{background:rgba(255,255,255,.025)}} .t td.num{{text-align:right;font-family:"JetBrains Mono",monospace;color:#e8f0fb}}
.tag{{display:inline-block;border:1px solid rgba(242,196,109,.35);background:rgba(242,196,109,.12);color:var(--gold);border-radius:99px;padding:2px 8px;font-family:"JetBrains Mono",monospace;font-size:10px;text-transform:uppercase;letter-spacing:.05em;margin-right:5px}}
.fig-grid{{display:grid;grid-template-columns:repeat(2,1fr);gap:16px}} .fig-card{{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:12px}} .fig-card img{{width:100%;border-radius:7px;border:1px solid var(--line);background:#fff}} .fig-card h4{{margin:8px 0 2px;color:var(--gold)}}
.path{{font-family:"JetBrains Mono",monospace;color:#c7d6ea;font-size:11px;word-break:break-all}}
footer{{border-top:1px solid var(--line);padding:28px;text-align:center;color:var(--muted);font-family:"JetBrains Mono",monospace;font-size:11px}}
@media(max-width:1050px){{.stats{{grid-template-columns:repeat(4,1fr)}}.wrap{{grid-template-columns:1fr}}.toc{{display:none}}.fig-grid,.grid2,.grid3,.grid4{{grid-template-columns:1fr}}h1{{font-size:48px}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">CROSS-Neo 2.0 · TCR-aware extension · 2026-05-09</div>
    <h1>CROSS-Neo-TCR<br>Extension Dossier</h1>
    <p class="lead">A recognition-aware diagnostic branch for TCR-available neoantigen cases. It does not replace the main pMHC ranker; it adds optional pMTnet + TEPCAM expert evidence, structure-readiness auditing, and wetlab prioritization.</p>
    <div class="stats">
      <div class="stat"><b>{total_rows:,}</b><span>TCR registry rows</span></div>
      <div class="stat"><b>{paired:,}</b><span>paired alpha/beta</span></div>
      <div class="stat"><b>{peptide_hla:,}</b><span>peptide-HLA rows</span></div>
      <div class="stat"><b>{overlap:,}</b><span>CROSS-Neo overlap</span></div>
      <div class="stat"><b>{exact:,}</b><span>exact TCR-pMHC</span></div>
      <div class="stat"><b>{ensemble.auprc:.3f}</b><span>ensemble AUPRC</span></div>
      <div class="stat"><b>{ensemble.auroc:.3f}</b><span>ensemble AUROC</span></div>
      <div class="stat"><b>{esc(top_wetlab.peptide)}</b><span>top wetlab peptide</span></div>
    </div>
    <div class="crumbs"><a href="index.html">Hub index</a> · <a href="assets/cross_neo_tcr_extension/{PACKET_NAME}">Reviewer packet ZIP</a> · <a href="cancer_vaccine_full_dossier.html">Cancer vaccine dossier</a> · <a href="paper1_fig8_mechanism_dossier.html">Paper 1 mechanism dossier</a></div>
  </div>
</header>
<div class="wrap">
<nav class="toc">
  <h4>Contents</h4>
  <a href="#decision">01 Decision</a>
  <a href="#registry">02 Registry</a>
  <a href="#linkage">03 CROSS-Neo linkage</a>
  <a href="#experts">04 External experts</a>
  <a href="#wetlab">05 Wetlab candidates</a>
  <a href="#md">06 MD escalation</a>
  <a href="#structure">07 Structure status</a>
  <a href="#figures">08 Figures</a>
  <a href="#claims">09 Claim boundary</a>
  <a href="#sources">10 Sources</a>
</nav>
<main>
<section id="decision">
  <h2><span class="num">01</span>Decision</h2>
  <div class="sub">What this module is allowed to claim now</div>
  <div class="grid grid3">
    <div class="card"><h4>Promote now</h4><p><span class="tag">Diagnostic</span><span class="tag">Wetlab</span></p><p>Use as optional expert for TCR-available cases and prioritization.</p></div>
    <div class="card"><h4>Hold</h4><p><span class="tag">Main method</span></p><p>No universal TCR-aware neoantigen claim until strict external paired-TCR benchmarks support it.</p></div>
    <div class="card"><h4>Best practical score</h4><div class="big">{ensemble.auprc:.3f}</div><p>pMTnet + TEPCAM mean ensemble AUPRC; AUROC {ensemble.auroc:.3f}.</p></div>
  </div>
  <div class="callout good"><b>Bottom line:</b> keep CROSS-Neo pMHC as the primary ranker. Add pMTnet + TEPCAM as a TCR-aware optional recognition layer when modelable TCR evidence exists.</div>
</section>

<section id="registry">
  <h2><span class="num">02</span>TCR Data Registry</h2>
  <div class="sub">Canonical local TCR-pMHC inventory</div>
  <div class="grid grid4">
    <div class="card"><div class="big">{any_tcr:,}</div><p>rows with any TCR sequence</p></div>
    <div class="card"><div class="big">{paired:,}</div><p>paired alpha/beta rows</p></div>
    <div class="card"><div class="big">{int(counts['cancer_context'].sum()):,}</div><p>cancer-context rows</p></div>
    <div class="card"><div class="big">{structures:,}</div><p>structure/PDB evidence rows</p></div>
  </div>
  {table(counts.sort_values("n", ascending=False), source_cols, max_rows=10, num_cols=source_num_cols)}
</section>

<section id="linkage">
  <h2><span class="num">03</span>CROSS-Neo Linkage</h2>
  <div class="sub">TCR evidence linked to main neoantigen registry</div>
  <p>{overlap:,}/2,715 CROSS-Neo rows have some TCR-resource overlap; {no_match:,} have no TCR match. Peptide-only and near-peptide links are diagnostic evidence only.</p>
  {table(linkage.sort_values("n_neo_rows", ascending=False), linkage_cols, max_rows=10, num_cols=linkage_num_cols)}
</section>

<section id="experts">
  <h2><span class="num">04</span>External TCR Experts</h2>
  <div class="sub">pMTnet, TEPCAM, ERGO-II, PanPep, NetTCR and structure tools</div>
  <div class="callout"><b>Benchmark:</b> six 100-positive challenge panels with same-panel shuffled-TCR decoys. Decoys are synthetic, so this is model selection evidence, not external SOTA validation.</div>
  {table(pooled, perf_cols, max_rows=10, num_cols=perf_num_cols)}
  <h3>Source-heldout ensemble</h3>
  {table(source_metrics, source_perf_cols, max_rows=10, num_cols=perf_num_cols)}
  <h3>Runtime model inventory</h3>
  {table(model_matrix[model_cols].head(12), model_cols, max_rows=12)}
</section>

<section id="wetlab">
  <h2><span class="num">05</span>Wetlab Candidate Rescoring</h2>
  <div class="sub">External expert scores over exact peptide-HLA TCR rows</div>
  <div class="grid grid2">
    <div class="card"><h4>Priority 1</h4><div class="big">{esc(wetlab.iloc[0].peptide)}</div><p>{esc(wetlab.iloc[0].hla_4digit)} · external mean {fmt_num(wetlab.iloc[0].external_tcr_expert_mean, 3)} · adjusted priority {fmt_num(wetlab.iloc[0].wetlab_priority_score_external_adjusted, 3)}</p></div>
    <div class="card"><h4>Priority 2</h4><div class="big">{esc(wetlab.iloc[1].peptide)}</div><p>{esc(wetlab.iloc[1].hla_4digit)} · external mean {fmt_num(wetlab.iloc[1].external_tcr_expert_mean, 3)} · adjusted priority {fmt_num(wetlab.iloc[1].wetlab_priority_score_external_adjusted, 3)}</p></div>
  </div>
  {table(wetlab, wet_cols, max_rows=15, num_cols=wet_num_cols)}
</section>

<section id="md">
  <h2><span class="num">06</span>MD Escalation Queue</h2>
  <div class="sub">Prediction-fragile, high-value candidates for expensive simulation</div>
  <div class="grid grid3">
    <div class="card"><h4>P0 MD cases</h4><div class="big">{p0_md}</div><p>exact paired TCR evidence plus external expert support</p></div>
    <div class="card"><h4>Top MD case</h4><div class="big">{esc(top_md.peptide)}</div><p>{esc(top_md.hla_4digit)} · score {fmt_num(top_md.md_escalation_score, 3)}</p></div>
    <div class="card"><h4>PDB pilots</h4><div class="big">{p0_pilot_ready}</div><p>existing TCR-pMHC pilot complexes; minimized {openmm_min_ok}/2; smoke dynamics {smoke_ok}/2</p></div>
  </div>
  <div class="callout warn"><b>Runtime reality:</b> OpenMM/MDTraj/PDBFixer are available locally, but only CPU/Reference OpenMM platforms are present. Production MD should run on a CUDA-enabled GPU pod. MD is not treated as proof of immunogenicity.</div>
  {table(md_queue, md_cols, max_rows=12, num_cols=md_num_cols)}
  <h3>P0 solved-structure pilots</h3>
  {table(p0_pilots, pilot_cols, max_rows=12, num_cols=pilot_num_cols)}
  <h3>OpenMM repair/minimization sanity check</h3>
  {table(min_qc, min_cols, max_rows=4, num_cols=min_num_cols)}
  <h3>OpenMM local smoke dynamics</h3>
  {table(smoke_qc, smoke_cols, max_rows=4, num_cols=smoke_num_cols)}
  <div class="callout good"><b>RunPod package:</b> <a href="assets/cross_neo_tcr_extension/{OPENMM_PILOT_PACKAGE}">{OPENMM_PILOT_PACKAGE}</a></div>
</section>

<section id="structure">
  <h2><span class="num">07</span>Structure Branch Status</h2>
  <div class="sub">Useful, but currently sequence-limited</div>
  <div class="callout warn"><b>Current boundary:</b> structure features are parser/QC and readiness diagnostics. No mutant-WT interface-delta claim yet because full TCR/MHC chains are missing for most rows.</div>
  {table(jobs, list(jobs.columns), max_rows=12, num_cols=job_num_cols)}
</section>

<section id="figures">
  <h2><span class="num">08</span>Figure Package</h2>
  <div class="sub">PNG/PDF figures generated from local result tables</div>
  <div class="fig-grid">
    {"".join(fig_cards)}
  </div>
</section>

<section id="claims">
  <h2><span class="num">09</span>Claim Boundary</h2>
  <div class="grid grid2">
    <div class="card"><h4>Allowed</h4><p>Optional TCR expert, diagnostic case interpretation, wetlab prioritization, and missingness reality.</p></div>
    <div class="card"><h4>Forbidden for now</h4><p>Clinical utility, universal TCR-aware prediction, structure-is-correct claims, and pathogen-to-cancer label transfer.</p></div>
  </div>
</section>

<section id="sources">
  <h2><span class="num">10</span>Sources And Paths</h2>
  <div class="sub">Every number on this page comes from these local artifacts</div>
  <p class="path">{TCR_OUT / "tcr_registry_source_counts.tsv"}</p>
  <p class="path">{TCR_OUT / "tcr_neo_linkage_summary.tsv"}</p>
  <p class="path">{TCR_OUT / "model_discovery/external_tcr_expert_benchmark/external_tcr_expert_benchmark_metrics.tsv"}</p>
  <p class="path">{TCR_OUT / "model_discovery/wetlab_external_tcr_expert_scores/wetlab_candidates_external_tcr_expert_ranked.tsv"}</p>
  <p class="path">{TCR_OUT / "md_escalation/md_escalation_queue_top20.tsv"}</p>
  <p class="path">{TCR_OUT / "TCR_EXTENSION_REVIEWER_QA.md"}</p>
  <p class="path">{TCR_OUT / "figures"}</p>
  <div class="callout good"><b>Reviewer packet:</b> <a href="assets/cross_neo_tcr_extension/{PACKET_NAME}">{PACKET_NAME}</a></div>
</section>
</main>
</div>
<footer>Generated by project/scripts/cross_neo_v2/24_build_tcr_dossier_page.py · CROSS-Neo-TCR diagnostic branch · claim-limited</footer>
</body>
</html>
"""
    PAGE.write_text(html_text)
    shutil.copy2(PAGE, WEB_PAGE)
    print(f"[tcr-dossier] wrote {PAGE}")
    print(f"[tcr-dossier] deployed {WEB_PAGE}")


if __name__ == "__main__":
    main()
