"""
v13 Task 4 — Multi-modality feasibility matrix (8 targets × 10 modalities).

Each (target, modality) cell carries:
  - evidence class: approved | clinical | preclinical | hypothesis
  - count of agents (best-effort)
  - short rationale
  - a curated representative agent (if any)

Curated from biology + current landscape, then augmented with ClinicalTrials.gov
v2 API queries for (drug × thyroid) where a named agent exists.

Output:
  $RES/modalities/target_modality_matrix.tsv
  $RES/modalities/target_modality_matrix.json
  $RES/modalities/ct_gov_hits.tsv
"""
from __future__ import annotations
import json, time
from pathlib import Path
import pandas as pd
import requests

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v13_drug_discovery"
MOD = RES / "modalities"
MOD.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "thca-v13 (kukshomr@gmail.com)", "Accept": "application/json"}

MODALITIES = ["SM", "mAb", "ADC", "bispecific", "siRNA", "PROTAC", "CAR_T", "radioligand", "mRNA", "ASO"]

# Curated landscape. evidence: approved > clinical > preclinical > hypothesis > none.
# count is best-effort existing agents (not necessarily thyroid-indicated).
LANDSCAPE = {
    "CYP1B1": {
        "SM":          ("clinical",    12, "flavonoid inhibitors (resveratrol Ph3, quercetin Ph3, luteolin Ph2); TMS preclinical; α-naphthoflavone tool compound", "resveratrol"),
        "mAb":         ("none",         0, "intracellular ER-membrane enzyme — antibodies non-viable", ""),
        "ADC":         ("none",         0, "intracellular target — not ADC-accessible", ""),
        "bispecific":  ("none",         0, "intracellular target", ""),
        "siRNA":       ("preclinical",  2, "CYP1B1 siRNA shown to sensitize tumors in vitro; no clinical program", ""),
        "PROTAC":      ("hypothesis",   0, "CYP1B1 PROTAC design feasible via heme-pocket ligand + CRBN warhead; no published program", ""),
        "CAR_T":       ("none",         0, "no surface epitope", ""),
        "radioligand": ("none",         0, "intracellular; no high-affinity radiopharmaceutical", ""),
        "mRNA":        ("none",         0, "no therapeutic rationale to replace/deliver", ""),
        "ASO":         ("preclinical",  1, "Gen-2 ASO knockdown reported in breast cancer models", ""),
    },
    "LDLR": {
        "SM":          ("approved",    30, "statins up-regulate LDLR via SREBP axis (atorvastatin, rosuvastatin, etc.); ezetimibe complementary", "atorvastatin"),
        "mAb":         ("approved",     3, "PCSK9 mAbs stabilize surface LDLR: evolocumab, alirocumab (approved); tafolecimab (clinical CN)", "evolocumab"),
        "ADC":         ("none",         0, "LDLR is widely expressed on hepatocytes — ADC would be hepatotoxic", ""),
        "bispecific":  ("hypothesis",   0, "PCSK9 × LDLR bispecifics described conceptually; no clinical program", ""),
        "siRNA":       ("approved",     2, "inclisiran (PCSK9 siRNA, Novartis) → raises LDLR indirectly; lepodisiran Ph3", "inclisiran"),
        "PROTAC":      ("none",         0, "LDLR up-regulation is therapeutic goal — degradation counter-productive", ""),
        "CAR_T":       ("none",         0, "no specific tumor-restricted expression", ""),
        "radioligand": ("preclinical",  1, "Cu-64 LDL particle imaging explored; not therapeutic", ""),
        "mRNA":        ("preclinical",  2, "LDLR mRNA replacement for homozygous FH under investigation (Moderna)", ""),
        "ASO":         ("approved",     1, "mipomersen (ApoB ASO) reduces LDL indirectly; volanesorsen (ApoC3) approved EU", "mipomersen"),
    },
    "TACSTD2": {
        "SM":          ("none",         0, "extracellular adhesion molecule — no small-molecule chemistry", ""),
        "mAb":         ("clinical",     5, "RN927C, datopotamab, sacituzumab all have naked-mAb PK data", ""),
        "ADC":         ("approved",     4, "🔥 sacituzumab govitecan (Trodelvy) approved breast/urothelial; datopotamab deruxtecan (Dato-DXd) Ph3 lung; SKB264 Ph3 CN; BAT8003 preclinical", "sacituzumab_govitecan"),
        "bispecific":  ("clinical",     2, "TROP2 × CD3 bispecifics entering Ph1 (Chugai CMG901-like, IGM)", ""),
        "siRNA":       ("preclinical",  1, "TROP2 siRNA anti-proliferative in vitro; no clinical", ""),
        "PROTAC":      ("none",         0, "extracellular target — not amenable to cytosolic E3 machinery", ""),
        "CAR_T":       ("clinical",     3, "TROP2 CAR-T Ph1 (NCT04151186 gastric; Legend Biotech)", ""),
        "radioligand": ("preclinical",  2, "⁶⁴Cu- and ¹⁷⁷Lu-labelled anti-TROP2 preclinical; Actinium Pharma exploratory", ""),
        "mRNA":        ("none",         0, "not a replacement target", ""),
        "ASO":         ("none",         0, "less competitive vs. ADC for a surface-expressed oncogene", ""),
    },
    "TMPRSS4": {
        "SM":          ("clinical",     6, "nafamostat (approved JP pancreatitis/DIC; tested as host-protease inhibitor); camostat (approved JP); gabexate; sivelestat", "nafamostat"),
        "mAb":         ("preclinical",  1, "anti-TMPRSS4 mAb against extracellular protease domain reported (Kim 2020)", ""),
        "ADC":         ("hypothesis",   0, "extracellular protease domain internalizing — conceptually ADC-amenable", ""),
        "bispecific":  ("none",         0, "no rationale beyond ADC", ""),
        "siRNA":       ("preclinical",  3, "TMPRSS4 siRNA reduces invasion in colon/thyroid models (Kim lab)", ""),
        "PROTAC":      ("hypothesis",   0, "plasma-membrane PROTAC uncommon; one exception family (LyTACs)", ""),
        "CAR_T":       ("none",         0, "low broad expression for CAR-T", ""),
        "radioligand": ("none",         0, "not developed", ""),
        "mRNA":        ("none",         0, "", ""),
        "ASO":         ("preclinical",  1, "ASO prototypes reported academically", ""),
    },
    "GABRB2": {
        "SM":          ("approved",    20, "benzodiazepines (diazepam, midazolam, clonazepam); Z-drugs (zolpidem); anaesthetics (etomidate, propofol); flumazenil antagonist", "diazepam"),
        "mAb":         ("none",         0, "ion-channel intracellular/membrane-embedded – not an antibody target", ""),
        "ADC":         ("none",         0, "", ""),
        "bispecific":  ("none",         0, "", ""),
        "siRNA":       ("preclinical",  1, "GABRB2 siRNA explored in neuropsychiatric models", ""),
        "PROTAC":      ("hypothesis",   0, "ion channel PROTAC (e.g., GABAAR) remains exploratory", ""),
        "CAR_T":       ("none",         0, "", ""),
        "radioligand": ("approved",     1, "[¹¹C]flumazenil PET tracer (benzodiazepine site) — diagnostic", ""),
        "mRNA":        ("none",         0, "", ""),
        "ASO":         ("preclinical",  1, "antisense modulation described", ""),
    },
    "PLEKHA6": {
        "SM":          ("hypothesis",   0, "scaffold PH domain — challenging pocket; no chemical probes", ""),
        "mAb":         ("none",         0, "intracellular junctional scaffold", ""),
        "ADC":         ("none",         0, "", ""),
        "bispecific":  ("none",         0, "", ""),
        "siRNA":       ("preclinical",  1, "PLEKHA6 siRNA reduces copper-ATP7A junction retention in vitro", ""),
        "PROTAC":      ("hypothesis",   0, "PH-domain PROTACs feasible (cf. AKT PH-domain PROTACs); de novo", ""),
        "CAR_T":       ("none",         0, "", ""),
        "radioligand": ("none",         0, "", ""),
        "mRNA":        ("none",         0, "", ""),
        "ASO":         ("preclinical",  1, "exploratory", ""),
    },
    "PTPRE": {
        "SM":          ("preclinical",  5, "broad PTP-family inhibitors (NSC-87877, suramin, ertiprotafib (discontinued)); no PTPRE-selective drug", ""),
        "mAb":         ("hypothesis",   0, "receptor PTP with extracellular domain — mAb conceptually feasible", ""),
        "ADC":         ("hypothesis",   0, "", ""),
        "bispecific":  ("none",         0, "", ""),
        "siRNA":       ("preclinical",  2, "PTPRE siRNA decreases proliferation in breast/lymphoma models", ""),
        "PROTAC":      ("hypothesis",   0, "PTP PROTACs described (e.g., SHP2 degraders) — extensible", ""),
        "CAR_T":       ("none",         0, "", ""),
        "radioligand": ("none",         0, "", ""),
        "mRNA":        ("none",         0, "", ""),
        "ASO":         ("preclinical",  1, "exploratory", ""),
    },
    "B3GNT3": {
        "SM":          ("preclinical",  2, "no selective B3GNT3 small molecule; pan-glycosyltransferase inhibitors non-selective", ""),
        "mAb":         ("none",         0, "Golgi-resident enzyme — inaccessible to mAbs", ""),
        "ADC":         ("none",         0, "", ""),
        "bispecific":  ("none",         0, "", ""),
        "siRNA":       ("preclinical",  2, "B3GNT3 knockdown disrupts PD-L1 glycosylation and sensitises to anti-PD-1 (Cell 2018)", ""),
        "PROTAC":      ("hypothesis",   0, "membrane-bound Golgi enzyme; PROTAC feasibility unclear", ""),
        "CAR_T":       ("none",         0, "", ""),
        "radioligand": ("none",         0, "", ""),
        "mRNA":        ("none",         0, "", ""),
        "ASO":         ("preclinical",  1, "exploratory", ""),
    },
}

# ---- ClinicalTrials.gov v2 lookups for flagship drugs ----
CTGOV_QUERIES = [
    ("sacituzumab", "thyroid"),
    ("sacituzumab", "cancer"),
    ("datopotamab", "thyroid"),
    ("trop2",       "thyroid"),
    ("inclisiran",  "cardiovascular"),
    ("evolocumab",  "ldlr"),
    ("nafamostat",  "cancer"),
    ("camostat",    "cancer"),
    ("resveratrol", "thyroid"),
    ("quercetin",   "thyroid"),
    ("cannabidiol", "thyroid"),
]


def ctgov_query(intr: str, cond: str | None = None):
    base = "https://clinicaltrials.gov/api/v2/studies"
    params = {"query.intr": intr, "pageSize": 20}
    if cond:
        params["query.cond"] = cond
    try:
        r = requests.get(base, params=params, headers=UA, timeout=30)
        if r.status_code != 200:
            return []
        js = r.json()
        rows = []
        for s in js.get("studies", []):
            p = s.get("protocolSection", {})
            ident = p.get("identificationModule", {})
            status = p.get("statusModule", {})
            design = p.get("designModule", {})
            conds = (p.get("conditionsModule", {}).get("conditions") or [])
            rows.append({
                "intr": intr,
                "cond": cond,
                "nct": ident.get("nctId"),
                "title": (ident.get("briefTitle") or "")[:120],
                "phase": ",".join(design.get("phases", []) or []),
                "status": status.get("overallStatus"),
                "start": (status.get("startDateStruct", {}) or {}).get("date"),
                "conditions": "; ".join(conds)[:200],
            })
        return rows
    except Exception as e:
        print(f"  CT.gov fail {intr}/{cond}: {e}")
        return []


def main():
    # Long-form TSV
    rows = []
    for gene, mods in LANDSCAPE.items():
        for m in MODALITIES:
            ev, ct, rat, rep = mods.get(m, ("none", 0, "", ""))
            rows.append({
                "gene": gene,
                "modality": m,
                "evidence_class": ev,
                "agent_count": ct,
                "representative_agent": rep,
                "rationale": rat,
            })
    df = pd.DataFrame(rows)
    df.to_csv(MOD / "target_modality_matrix.tsv", sep="\t", index=False)

    # Wide JSON (evidence_class) for dashboard heatmap
    ev_score = {"approved": 4, "clinical": 3, "preclinical": 2, "hypothesis": 1, "none": 0}
    wide = {}
    for gene, mods in LANDSCAPE.items():
        wide[gene] = {m: {
            "evidence": mods[m][0],
            "score": ev_score[mods[m][0]],
            "count": mods[m][1],
            "agent": mods[m][3],
            "rationale": mods[m][2],
        } for m in MODALITIES}
    (MOD / "target_modality_matrix.json").write_text(json.dumps(wide, indent=2))

    # CT.gov hits
    ct_rows = []
    for intr, cond in CTGOV_QUERIES:
        hits = ctgov_query(intr, cond)
        ct_rows.extend(hits)
        print(f"  CT.gov {intr:>15} × {cond or '-':<15} → {len(hits)} studies")
        time.sleep(0.5)
    if ct_rows:
        pd.DataFrame(ct_rows).to_csv(MOD / "ct_gov_hits.tsv", sep="\t", index=False)
        print(f"Total CT.gov rows: {len(ct_rows)} → {MOD/'ct_gov_hits.tsv'}")

    # Pretty print heatmap
    print("\n=== Modality evidence matrix (score 0=none 4=approved) ===")
    print(f"{'gene':<9} " + " ".join(f"{m:>10}" for m in MODALITIES))
    for g in LANDSCAPE:
        print(f"{g:<9} " + " ".join(f"{ev_score[LANDSCAPE[g][m][0]]:>10}" for m in MODALITIES))


if __name__ == "__main__":
    main()
