"""Generate all 10 schematic + data figures for the one-page audit, fast."""
from __future__ import annotations
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

# Register Korean fonts (NanumGothic + Noto Sans CJK KR) so 한글 doesn't break
for _p in ["/home/seungho/.local/share/fonts/NanumGothic-Regular.ttf",
           "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
           "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"]:
    if os.path.exists(_p):
        try: fm.fontManager.addfont(_p)
        except Exception: pass

OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)

# ---------- common style ----------
plt.rcParams.update({
    "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": ["NanumGothic", "Noto Sans CJK JP", "DejaVu Sans"],
    "axes.unicode_minus": False,
})
RED, BLUE, GREY, ORANGE, GREEN, GOLD = "#c0392b", "#2c5e9c", "#7f8fa6", "#e67e22", "#3C6B4F", "#d4a017"

# ============================================================
# fig01 — clinical workflow current vs future
# ============================================================
fig, ax = plt.subplots(figsize=(13, 5.5))
ax.set_xlim(0, 14); ax.set_ylim(0, 6); ax.axis("off")
def box(x, y, w, h, text, color, fontsize=9):
    r = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04", linewidth=1.2,
                       edgecolor="black", facecolor=color)
    ax.add_patch(r)
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fontsize, wrap=True)
def arrow(x1, y1, x2, y2):
    a = FancyArrowPatch((x1,y1),(x2,y2), arrowstyle="-|>", mutation_scale=18, color="black", lw=1.5)
    ax.add_patch(a)

# row 1 — current
ax.text(0.2, 5.4, "CURRENT WORKFLOW (current standard)", fontsize=12, fontweight="bold")
box(0.3, 4.0, 2.2, 1.0, "Surgery\n(thyroidectomy)", "#FFE5D9")
box(2.8, 4.0, 2.2, 1.0, "RAI attempt\n(empirical)", "#FFE5D9")
box(5.3, 4.0, 2.4, 1.0, "RAI failure\nrecognition delay\n(months–years)", "#FFCCCC")
box(8.0, 4.0, 2.4, 1.0, "Hemato-oncology /\nsystemic review", "#FFE5D9")
box(10.7, 4.0, 2.6, 1.0, "Aggressive disease\noften progressed by then", "#FFAAAA", fontsize=9)
for x1, x2 in [(2.5,2.8),(5.0,5.3),(7.7,8.0),(10.4,10.7)]: arrow(x1,4.5,x2,4.5)

# row 2 — proposed (hypothesis only)
ax.text(0.2, 3.0, "PROPOSED HYPOTHESIS — early triage scaffold (NOT validated; NOT for current clinical use)", fontsize=12, fontweight="bold", color=BLUE)
box(0.3, 1.4, 2.2, 1.0, "Surgery", "#E8F4FD")
box(2.8, 1.4, 2.4, 1.0, "Compact RNA readout\n(8-gene RAI lineage)", "#E8F4FD")
box(5.5, 1.4, 2.5, 1.0, "Early flag of\nRAI-lineage failure biology", "#FFF4D6")
box(8.3, 1.4, 2.4, 1.0, "Earlier review of\nsystemic / molecular options", "#E8F4FD")
box(11.0, 1.4, 2.4, 1.0, "Reduced delay\n(hypothesis only)", "#D4F4DD", fontsize=9)
for x1, x2 in [(2.5,2.8),(5.2,5.5),(8.0,8.3),(10.7,11.0)]: arrow(x1,1.9,x2,1.9)

ax.text(0.2, 0.4, "★ Triage scaffold only. NOT a validated RAI response predictor. NOT a treatment-selection tool. Prospective trial required.",
        fontsize=9, color="#7B1F2A", style="italic")
plt.tight_layout()
plt.savefig(OUT/"fig01_clinical_workflow_current_vs_future.png", dpi=150, bbox_inches="tight"); plt.close()


# ============================================================
# fig02 — dataset suitability matrix
# ============================================================
data = [
    # (name, scope, role, n, badge)
    ("TCGA-THCA",          "primary mixed",  "discovery",         504,  "main"),
    ("TCGA HM450",         "primary mixed",  "mechanism",         503,  "main"),
    ("cBioPortal SV/TERT", "primary mixed",  "discovery",         "504/36", "main"),
    ("MSK-IMPACT (Landa)", "advanced only",  "aggressive replication", 117, "main"),
    ("GSE76039",           "advanced only",  "external valid.",   37,   "main"),
    ("GSE33630",           "ATC+PTC+N",      "external valid.",   105,  "main"),
    ("GSE65144",           "ATC+N",          "external valid.",   25,   "main"),
    ("GSE29265",           "ATC+PTC+N",      "external valid.",   49,   "supp"),
    ("GSE53157",           "PDTC+PTC+FVPTC+FTC+N", "sensitivity",  26,   "supp"),
    ("GSE213647 (Lee 24)", "Korean PTC",     "external valid.",   632,  "main"),
    ("K2/PRJEB11591",      "Korean PTC",     "external valid.",   260,  "main"),
    ("GSE184362 (Pu 21)",  "scRNA primary",  "mechanism",         "6/158K", "supp"),
    ("GSE193581 (Lu 23)",  "scRNA + cell line", "mechanism",      23,   "supp"),
    ("GSE241184 (Phase 1)","scRNA",          "author-indep.",     1,    "supp"),
    ("GSE250521 (spatial)","ST PT/PTC/LPTC/ATC", "spatial supp",  16,   "supp"),
    ("GSE126698",          "RNA-seq metadata-only", "DROP",       28,   "drop"),
    ("GSE60542",           "PTC nodal mets", "DROP (lymphoid)",   92,   "drop"),
    ("GSE286332",          "PTC + PTC+HT",   "Paper 2 territory", 18,   "drop"),
    ("GSE151179",          "post-RAI",       "Paper 3 territory", 52,   "drop"),
    ("Pozdeyev 2018 supp", "advanced 779",   "supp cite only",    779,  "supp"),
    ("Bundang SNUH",       "outreach",       "future",            0,    "drop"),
]
fig, ax = plt.subplots(figsize=(13.5, 9))
ax.axis("off")
ax.set_title("Paper 1 — dataset suitability matrix (main / supp / drop)", fontsize=13, fontweight="bold")
header = ["Dataset", "Scope", "Role", "n", "Badge"]
for i, h in enumerate(header):
    ax.text([0.02, 0.28, 0.50, 0.78, 0.92][i], 0.96, h, fontweight="bold", fontsize=10, transform=ax.transAxes)
badge_color = {"main": GOLD, "supp": "#cfd8dc", "drop": "#bf3434"}
for i, (name, scope, role, n, badge) in enumerate(data):
    y = 0.92 - i*0.043
    ax.text(0.02, y, name, fontsize=9, transform=ax.transAxes, fontweight="bold")
    ax.text(0.28, y, scope, fontsize=9, transform=ax.transAxes)
    ax.text(0.50, y, role, fontsize=9, transform=ax.transAxes)
    ax.text(0.78, y, str(n), fontsize=9, transform=ax.transAxes)
    ax.add_patch(Rectangle((0.91, y-0.012), 0.07, 0.028, transform=ax.transAxes,
                           facecolor=badge_color[badge], edgecolor="black", linewidth=0.5))
    ax.text(0.945, y+0.002, badge.upper(), fontsize=8, transform=ax.transAxes, ha="center",
            color="white" if badge=="drop" else "black", fontweight="bold")
plt.tight_layout()
plt.savefig(OUT/"fig02_dataset_suitability_matrix.png", dpi=150, bbox_inches="tight"); plt.close()


# ============================================================
# fig03 — gene set hierarchy
# ============================================================
fig, ax = plt.subplots(figsize=(13, 7))
ax.set_xlim(0, 14); ax.set_ylim(0, 9); ax.axis("off")
ax.set_title("Paper 1 — gene-set hierarchy (literature-curated discovery → deployment readout)",
             fontsize=12, fontweight="bold")
# top — pan-genome (data driven)
box(4.5, 7.6, 5, 1.0, "Pan-genome top-5000 (MAD)\nARI = 0.92  (data-driven; literature-free)", "#E8F4FD", 10)
# middle — TIERA67 (literature)
box(4.5, 5.8, 5, 1.0, "TIERA67  (literature-curated 67-gene candidate pool)\nARI = 0.90", "#FFF4D6", 10)
# branches
ax.text(0.5, 4.4, "7 sub-categories (Supp Table S1):", fontsize=10, fontweight="bold")
cats = [("TDS_core 16",    "Yoo 2016 PLOS Genet"),
        ("MAPK_output 10", "ERK feedback"),
        ("Driver_anchor 12","TCGA 2014 + COSMIC"),
        ("Aggressive 10",  "p53/PTEN/etc"),
        ("Dediff/EMT 10",  "Kalluri 2009 EMT"),
        ("Immune-light 5", "checkpoint/Treg"),
        ("Lineage extra 4","HHEX/GLIS3/IYD/MET/KLK10")]
for i, (n, src) in enumerate(cats):
    box(0.3+i*1.92, 3.4, 1.85, 0.85, f"{n}\n{src}", "#F0F0F0", 7.5)
# bottom — RAI_8
box(0.5, 1.6, 4.5, 1.2, "RAI_8 (8-gene)  ★ MAIN DEPLOYABLE READOUT\nSLC5A5/TPO/TG/TSHR/PAX8/NKX2-1/FOXE1/DIO1\nARI 0.49 (modest by design)  •  AUC 0.962", RED, 9)
box(5.3, 1.6, 4.0, 1.2, "THYROID_NONOVERLAP (8)\nSLC26A4/IYD/DUOX1/DUOX2/TFF3/HHEX/GLIS3/DIO2\nzero gene overlap with RAI_8", BLUE, 9)
box(9.6, 1.6, 4.0, 1.2, "TF_collapse (4 lineage TFs)\nFOXE1/NKX2-1/PAX8/HHEX\nDamante 4-TF backbone", "#cfd8dc", 9)
# arrows
arrow(7.0, 7.6, 7.0, 6.8)
arrow(7.0, 5.8, 7.0, 4.7)
for i in range(7):
    arrow(0.3+i*1.92+0.92, 3.4, 0.3+i*1.92+0.92, 2.85)
# bottom — drop
box(0.5, 0.1, 6, 1.0, "DROP / DEMOTED:  TROP2 (no main); 10-gene q5 (TRPS1 ad-hoc); 55-gene driver-excluded (superseded by TIERA67)", "#bf3434", 9)
box(7.0, 0.1, 6, 1.0, "TDS-16 comparator: AUC 0.975  (ΔAUC vs 8-gene = 0.013, NS)", "#cfd8dc", 9)
plt.tight_layout()
plt.savefig(OUT/"fig03_gene_set_hierarchy.png", dpi=150, bbox_inches="tight"); plt.close()


# ============================================================
# fig04 — main evidence chain
# ============================================================
fig, ax = plt.subplots(figsize=(14, 4.5))
ax.set_xlim(0, 14); ax.set_ylim(0, 4); ax.axis("off")
ax.set_title("Paper 1 — main evidence chain  (deployment-vs-discovery 2-layer architecture)",
             fontsize=12, fontweight="bold")
chain = [
    ("a) Driver-orthogonal\nlineage axis", "BRAF V600E vs WT mRNA\nd = -0.04 (NS)", "#FFE5D9"),
    ("b) Compact 8-gene\nreadout", "AUC 0.962 (5-fold CV)\nΔ vs 16 = 0.013 NS", "#FFF4D6"),
    ("c) External validation\n(GPL570 + Lee 2024)", "ρ -0.84 to -0.94\nDM1_like vs NONOVERLAP", "#E8F4FD"),
    ("d) Aggressive replication\n(MSK + Landa)", "Pooled OS HR 2.53\n(switch to PFI primary)", "#D4F4DD"),
    ("e) Mechanism support\n(HM450 + TF collapse)", "TPO d 2.30; mean β 0.385\nvs DM2 0.253", GOLD),
]
for i, (h, sub, c) in enumerate(chain):
    x = 0.3 + i*2.75
    box(x, 1.5, 2.5, 2.0, f"{h}\n\n{sub}", c, 9)
    if i < 4: arrow(x+2.55, 2.5, x+2.75, 2.5)
ax.text(0.3, 0.6, "★ 8-gene = deployable readout / TIERA67 + pan-genome top-5000 = discovery axis (2-layer architecture)",
        fontsize=10, color=RED, fontweight="bold")
plt.tight_layout()
plt.savefig(OUT/"fig04_main_evidence_chain.png", dpi=150, bbox_inches="tight"); plt.close()


# ============================================================
# fig05 — external validation summary chart
# ============================================================
cohorts = ["GSE33630\n(105)", "GSE29265\n(49)", "GSE65144\n(25)", "GSE53157\n(26)"]
spearman_dm1_nonoverlap = [-0.93, -0.85, -0.94, -0.84]
fig, ax = plt.subplots(figsize=(10, 5))
bars = ax.bar(cohorts, spearman_dm1_nonoverlap, color=[BLUE]*4, edgecolor="black")
ax.set_ylim(-1.0, 0.0)
ax.axhline(0, color="black", lw=0.8)
ax.axhline(-0.7, color=RED, lw=1, linestyle="--", label="strong-anticorrelation threshold (-0.7)")
ax.set_ylabel("Spearman ρ  (DM1_like_score vs THYROID_NONOVERLAP_score)")
ax.set_title("External GPL570 validation — DM1 axis vs zero-overlap orthogonal lineage panel\n(every cohort independently replicates ρ ≤ -0.84)", fontsize=11)
for b, v in zip(bars, spearman_dm1_nonoverlap):
    ax.text(b.get_x()+b.get_width()/2, v-0.04, f"{v}", ha="center", fontsize=10, fontweight="bold", color="white")
ax.legend()
plt.tight_layout()
plt.savefig(OUT/"fig05_external_validation_summary.png", dpi=150, bbox_inches="tight"); plt.close()


# ============================================================
# fig06 — aggressive cohort context chart
# ============================================================
fig, ax = plt.subplots(figsize=(11, 5))
labels = ["TCGA-THCA\nprimary (504)", "MSK-IMPACT\nadvanced (117)", "GSE76039 Landa\nadvanced (37)", "GSE53157\nPDTC arm (5)"]
adv_pct = [10, 100, 100, 100]   # advanced % of cohort
event_pct_pfi = [11.6, None, None, None]  # PFI events
colors = ["#94d4a4", RED, RED, ORANGE]
ax.barh(range(len(labels)), adv_pct, color=colors, edgecolor="black")
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels)
ax.set_xlabel("% advanced disease (PDTC + ATC) within cohort")
ax.set_xlim(0, 110)
ax.set_title("Aggressive cohort context — explicit bias disclosure required for Landa / MSK-IMPACT\n(disclosure: SA6 MSK-IMPACT bias panel + Methods M1)",
             fontsize=11)
for i, p in enumerate(adv_pct):
    ax.text(p+1, i, f"{p}%", va="center", fontsize=10, fontweight="bold")
ax.axvline(50, color="black", lw=0.5, linestyle=":")
ax.text(52, -0.7, "advanced-enriched →", fontsize=9, color=RED)
plt.tight_layout()
plt.savefig(OUT/"fig06_aggressive_cohort_context.png", dpi=150, bbox_inches="tight"); plt.close()


# ============================================================
# fig07 — main vs supplement vs drop figure map
# ============================================================
items = [
    ("Driver-orthogonal landscape (Fig 2)", "main"),
    ("8-gene + ARI ladder (Fig 1)", "main"),
    ("External GPL570 4-cohort (Fig 3 incl. scatter grid)", "main"),
    ("PFI + age-adjusted Cox (Fig 4 reframed; OS sparse)", "main"),
    ("HM450 mechanism (Fig 5/8)", "main"),
    ("Korean GSE213647 + K2 (Fig 6)", "main"),
    ("DM1 sub-A vs sub-B (with versioning caveat)", "main"),
    ("Single-cell GSE184362/193581/241184", "supp"),
    ("Spatial GSE250521 (PT→ATC marginal)", "supp"),
    ("MSK-IMPACT bias panel (SA6)", "supp"),
    ("8-gene self + 16-gene zero-overlap matrix", "main"),
    ("Pan-genome top-5000 ARI ladder", "supp"),
    ("TIERA67 candidate pool table (S1)", "supp"),
    ("PDTC GSE53157 sensitivity", "supp"),
    ("Landa 2016 cite-save heatmap (SD1)", "supp"),
    ("10-gene q5_thyroid_tf_network (TRPS1 ad-hoc)", "drop"),
    ("H&E / WSI pathology DM1 inference", "drop"),
    ("TROP2 main / title / spatial claim", "drop"),
    ("GSE126698 Series Matrix metadata-only", "drop"),
    ("GSE60542 lymphoid confounding", "drop"),
]
counts = {"main": 0, "supp": 0, "drop": 0}
for _, c in items: counts[c] += 1
fig, ax = plt.subplots(figsize=(13, 8))
ax.axis("off")
ax.set_title(f"Paper 1 — figure registry  (main {counts['main']} / supp {counts['supp']} / drop {counts['drop']})",
             fontsize=12, fontweight="bold")
col_map = {"main": "#94d4a4", "supp": "#f9c863", "drop": "#bf3434"}
for i, (n, c) in enumerate(items):
    y = 0.95 - i*0.044
    ax.add_patch(Rectangle((0.02, y-0.014), 0.06, 0.030, transform=ax.transAxes,
                           facecolor=col_map[c], edgecolor="black"))
    ax.text(0.05, y+0.001, c.upper(), fontsize=8, ha="center", transform=ax.transAxes,
            color="white" if c=="drop" else "black", fontweight="bold")
    ax.text(0.10, y+0.001, n, fontsize=10, transform=ax.transAxes)
plt.tight_layout()
plt.savefig(OUT/"fig07_main_vs_supplement_map.png", dpi=150, bbox_inches="tight"); plt.close()


# ============================================================
# fig08 — claim boundary board (allowed vs forbidden)
# ============================================================
allowed = [
    "compact RAI-lineage readout",
    "direction-consistent lineage silencing",
    "advanced-disease replication",
    "zero-overlap module validation",
    "driver-orthogonal axis",
    "candidate triage scaffold (hypothesis)",
    "supportive mechanism arm",
    "consistent with / supports",
    "PFI primary; OS secondary",
    "East-Asian generalizability",
]
forbidden = [
    "TROP2 vulnerability (title)",
    "TROP2-targetable (headline)",
    "validated RAI response predictor",
    "treatment selection tool",
    "clinical utility proven",
    "progression proven / survival validation",
    "fusion validation (causal)",
    "H&E inferable / pathology DM1",
    "demonstrates causality",
    "all risks resolved",
]
fig, ax = plt.subplots(figsize=(13, 7))
ax.axis("off")
ax.set_title("Paper 1 — claim boundary  (allowed wording ↔ forbidden wording)", fontsize=12, fontweight="bold")
ax.text(0.05, 0.93, "[OK]  ALLOWED", fontsize=14, fontweight="bold", color=GREEN, transform=ax.transAxes)
ax.text(0.55, 0.93, "[NO]  FORBIDDEN", fontsize=14, fontweight="bold", color=RED, transform=ax.transAxes)
for i, w in enumerate(allowed):
    ax.text(0.05, 0.85-i*0.07, f"+  {w}", fontsize=10, transform=ax.transAxes, color=GREEN)
for i, w in enumerate(forbidden):
    ax.text(0.55, 0.85-i*0.07, f"-  {w}", fontsize=10, transform=ax.transAxes, color=RED)
plt.tight_layout()
plt.savefig(OUT/"fig08_claim_boundary_board.png", dpi=150, bbox_inches="tight"); plt.close()


# ============================================================
# fig09 — reviewer risk heatmap
# ============================================================
risks = [
    ("Why 8 genes not TDS-16",                "medium"),
    ("Why not 67",                            "medium"),
    ("Just dedifferentiation relabeled",      "high"),
    ("External cohorts aggressive-biased",    "medium"),
    ("Predicts RAI response",                 "very_high"),
    ("Clinical utility overclaim",            "very_high"),
    ("TROP2 demotion",                        "very_high"),
    ("H&E / pathology branch missing",        "very_high"),
    ("NONOVERLAP novelty enough",             "low"),
    ("Cheaper TDS surrogate",                 "low"),
    ("OS swap / age confound",                "high"),
    ("Stage doesn't separate",                "high"),
    ("Sub-A/sub-B versioning gap",            "high"),
    ("Spatial signal weak",                   "high"),
    ("DM1 sub-B = NBNR claim",                "medium"),
    ("TIERA67 acronym (project-internal)",    "low"),
    ("Mechanism causality",                   "medium"),
    ("Korean cohort applicability",           "low"),
    ("Decitabine + I-131 trial framing",      "medium"),
    ("Selpercatinib reflex framing",          "medium"),
]
risk_color = {"low":"#94d4a4","medium":"#f9c863","high":ORANGE,"very_high":RED}
fig, ax = plt.subplots(figsize=(13, 7.5))
ax.axis("off")
ax.set_title("Reviewer-risk heatmap (top 20)", fontsize=12, fontweight="bold")
for i, (q, r) in enumerate(risks):
    y = 0.96 - i*0.046
    ax.add_patch(Rectangle((0.02, y-0.013), 0.16, 0.030, transform=ax.transAxes,
                           facecolor=risk_color[r], edgecolor="black"))
    ax.text(0.10, y+0.002, r.upper().replace("_"," "), fontsize=8, ha="center",
            color="white" if r in {"high","very_high"} else "black", transform=ax.transAxes,
            fontweight="bold")
    ax.text(0.20, y+0.002, q, fontsize=10, transform=ax.transAxes)
plt.tight_layout()
plt.savefig(OUT/"fig09_reviewer_risk_heatmap.png", dpi=150, bbox_inches="tight"); plt.close()


# ============================================================
# fig10 — keep / move / drop figure map (compressed registry view)
# ============================================================
fig, ax = plt.subplots(figsize=(12, 6))
keep = ["F1 clinical workflow","F2 axis robustness + 8-gene","F3 GPL570 external","F4 PFI + Cox","F5 mechanism HM450",
        "8/16-gene self+cross matrix","Korean GSE213647 + K2","Landa cite-save SD1"]
supp = ["Single-cell sc 3-cohort","Spatial GSE250521 (marginal)","TIERA67 ARI ladder",
        "PDTC GSE53157 sensitivity","Sub-A/sub-B versioning panel","MSK-IMPACT bias panel SA6"]
drop = ["10-gene q5 (TRPS1 ad-hoc)","H&E pathology DM1","TROP2 main/title","GSE126698 metadata-only","GSE60542 lymphoid"]
ax.set_title("Keep / Supp / Drop figure registry", fontsize=12, fontweight="bold")
ax.axis("off")
for i, (group, items, col) in enumerate([("KEEP MAIN", keep, GREEN), ("SUPP", supp, "#d4a017"), ("DROP", drop, RED)]):
    x = 0.05 + i*0.32
    ax.add_patch(Rectangle((x-0.005, 0.92-0.01), 0.30, 0.06, transform=ax.transAxes,
                           facecolor=col, edgecolor="black"))
    ax.text(x+0.14, 0.94, group, fontsize=11, ha="center", fontweight="bold",
            color="white" if group=="DROP" else "black", transform=ax.transAxes)
    for j, it in enumerate(items):
        ax.text(x, 0.86-j*0.06, "• "+it, fontsize=9, transform=ax.transAxes)
plt.tight_layout()
plt.savefig(OUT/"fig10_keep_move_drop_figure_registry.png", dpi=150, bbox_inches="tight"); plt.close()

print("[done] 10 figures written to", OUT)
