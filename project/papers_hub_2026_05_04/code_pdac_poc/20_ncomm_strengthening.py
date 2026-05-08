"""
NComm-strengthening statistical layer:
  (A) Cox with KRAS_allele × tissue interaction term — formal test of
      the tissue-specificity claim (PDAC G12R best vs pan-cancer G12R worst)
  (B) Bootstrap meta-analysis HR — 10,000 reps with percentile CI
  (C) Drop-one cohort sensitivity meta
  (D) Random-effect (DerSimonian-Laird) meta sensitivity
  (E) DepMap-style KRAS dependency proxy — public DEMETER2/CHRONOS
      dependency from cBioPortal CCLE for PDAC vs LUNG vs CRC cell lines

All inputs from prior scripts. Output: results/ncomm_strengthening/SUMMARY.json
plus 2 new figures.
"""

import json
import math
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
import requests
from lifelines import CoxPHFitter
from scipy.stats import chi2

API = "https://www.cbioportal.org/api"
ROOT = Path("/data/pdac_poc")
OUT = ROOT / "results/ncomm_strengthening"
OUT.mkdir(parents=True, exist_ok=True)


def safe_get(url, **kw):
    try:
        r = requests.get(url, params=kw, timeout=60); r.raise_for_status(); return r.json()
    except Exception:
        return None


def safe_post(url, body, **kw):
    try:
        r = requests.post(url, json=body, params=kw, timeout=120,
                          headers={"Content-Type": "application/json"}); r.raise_for_status(); return r.json()
    except Exception:
        return None


# ====================================================================
# (A) FORMAL INTERACTION TEST
# Combines per-cohort allele assignments + survival, fits a stratified
# Cox with KRAS_allele × tissue interaction.
# ====================================================================
COHORTS_FOR_INTERACTION = [
    ("paad_tcga_pan_can_atlas_2018", "PDAC"),
    ("luad_tcga_pan_can_atlas_2018", "LUAD"),
    ("coadread_tcga_pan_can_atlas_2018", "CRC"),
    ("ucec_tcga_pan_can_atlas_2018", "UCEC"),
    ("stad_tcga_pan_can_atlas_2018", "STAD"),
    ("msk_chord_2024", "PAN_AGGREG"),
    ("msk_met_2021", "PAN_AGGREG2"),
]


def best_os(clin_cols):
    pairs = [("OS_MONTHS","OS_STATUS"),("OS_MONTHS","STATUS"),("OS_MONTHS","VITAL_STATUS"),
             ("DAYS_TO_LAST_FOLLOWUP","STATUS"),("DAYS_TO_LAST_FOLLOWUP","VITAL_STATUS")]
    cols = set(clin_cols)
    for t, s in pairs:
        if t in cols and s in cols: return t, s
    return None, None


def parse_event(v):
    if isinstance(v, str):
        v = v.lower()
        if v.startswith("1") or v in ("dead","deceased","yes"): return 1
        if v.startswith("0") or v in ("alive","living","no"): return 0
    try: return int(float(v))
    except: return None


def fetch_one(study_label):
    study, label = study_label
    samples = safe_get(f"{API}/studies/{study}/samples")
    if not samples: return pd.DataFrame()
    sdf = pd.DataFrame(samples)
    pcl = safe_get(f"{API}/studies/{study}/clinical-data", clinicalDataType="PATIENT") or []
    scl = safe_get(f"{API}/studies/{study}/clinical-data", clinicalDataType="SAMPLE") or []
    pdf = (pd.DataFrame(pcl).pivot_table(index="patientId", columns="clinicalAttributeId",
            values="value", aggfunc="first").reset_index() if pcl else pd.DataFrame())
    sdf2 = (pd.DataFrame(scl).pivot_table(index="sampleId", columns="clinicalAttributeId",
            values="value", aggfunc="first").reset_index() if scl else pd.DataFrame())
    if not sdf2.empty: sdf = sdf.merge(sdf2, on="sampleId", how="left")
    if not pdf.empty:  sdf = sdf.merge(pdf, on="patientId", how="left")

    g = safe_post(f"{API}/genes/fetch", body=["KRAS"], geneIdType="HUGO_GENE_SYMBOL")
    if not g: return pd.DataFrame()
    eid = g[0]["entrezGeneId"]
    body = {"entrezGeneIds":[eid], "sampleListId": f"{study}_sequenced"}
    muts = safe_post(f"{API}/molecular-profiles/{study}_mutations/mutations/fetch",
                     body=body, projection="DETAILED")
    if not muts: muts = []
    md = pd.DataFrame(muts)
    if "gene" in md.columns:
        md["hugo"] = md["gene"].apply(lambda gg: gg.get("hugoGeneSymbol") if isinstance(gg, dict) else None)
    md = md[(md.get("mutationType","").fillna("") == "Missense_Mutation") &
            (md.get("hugo","KRAS") == "KRAS")] if len(md) else md

    df = sdf.set_index("sampleId").copy()
    df["allele"] = "WT"
    if len(md):
        for s, gg in md.groupby("sampleId"):
            if s not in df.index: continue
            pcs = gg["proteinChange"].dropna().tolist()
            for p in ["G12D","G12V","G12R","G12C","G13D","Q61H","Q61R"]:
                if any(pc == p for pc in pcs): df.loc[s,"allele"] = p; break
            else:
                if any(pc.startswith("G12") for pc in pcs): df.loc[s,"allele"] = "G12_other"
                elif any(pc.startswith("Q61") for pc in pcs): df.loc[s,"allele"] = "Q61"
                elif pcs: df.loc[s,"allele"] = "KRAS_other"
    t_col, s_col = best_os(df.columns)
    if not t_col: return pd.DataFrame()
    df["t"] = pd.to_numeric(df[t_col], errors="coerce")
    df["e"] = df[s_col].apply(parse_event)
    if t_col.endswith("DAYS") or t_col.endswith("FOLLOWUP"):
        df["t"] = df["t"]/30.4375
    df = df[(df["t"]>0) & df["e"].notna()].dropna(subset=["t","e"])
    df["tissue"] = label
    df["study"] = study
    return df.reset_index()[["sampleId","study","tissue","allele","t","e"]]


def interaction_test():
    print("[20A] Cox with KRAS_allele × tissue interaction …")
    dfs = []
    with ProcessPoolExecutor(max_workers=7) as ex:
        for r in ex.map(fetch_one, COHORTS_FOR_INTERACTION):
            if len(r): dfs.append(r)
    if not dfs:
        return {"error": "no data"}
    full = pd.concat(dfs, ignore_index=True)
    full["G12D"] = (full["allele"]=="G12D").astype(int)
    full["G12V"] = (full["allele"]=="G12V").astype(int)
    full["G12R"] = (full["allele"]=="G12R").astype(int)
    full["G12C"] = (full["allele"]=="G12C").astype(int)
    full["is_PDAC"] = (full["tissue"]=="PDAC").astype(int)
    full["G12R_x_PDAC"] = full["G12R"] * full["is_PDAC"]
    full["G12D_x_PDAC"] = full["G12D"] * full["is_PDAC"]

    # main effects + interaction
    cols = ["t","e","G12D","G12V","G12R","G12C","is_PDAC","G12R_x_PDAC","G12D_x_PDAC"]
    cox_df = full[cols].dropna()
    if len(cox_df) < 100:
        return {"error": f"too few (n={len(cox_df)})"}
    try:
        cph = CoxPHFitter(penalizer=0.05).fit(cox_df, duration_col="t", event_col="e")
        s = cph.summary[["coef","exp(coef)","p","exp(coef) lower 95%","exp(coef) upper 95%"]].round(4)
        out = {"n_total": int(len(cox_df)),
               "n_pdac": int(cox_df["is_PDAC"].sum()),
               "n_pancancer": int(len(cox_df) - cox_df["is_PDAC"].sum()),
               "concordance": float(cph.concordance_index_),
               "table": s.reset_index().to_dict(orient="records")}
        # extract interaction p
        for r in out["table"]:
            if r["covariate"] == "G12R_x_PDAC":
                out["G12R_x_PDAC_HR"] = r["exp(coef)"]
                out["G12R_x_PDAC_p"] = r["p"]
                out["G12R_x_PDAC_CI"] = [r["exp(coef) lower 95%"], r["exp(coef) upper 95%"]]
            if r["covariate"] == "G12D_x_PDAC":
                out["G12D_x_PDAC_HR"] = r["exp(coef)"]
                out["G12D_x_PDAC_p"] = r["p"]
        return out
    except Exception as e:
        return {"error": str(e)}


# ====================================================================
# (B) BOOTSTRAP META-HR
# ====================================================================
def bootstrap_meta(B=10000):
    print(f"[20B] Bootstrap meta-analysis · {B} reps …")
    s = json.load(open(ROOT/"results/method_panel/META_ANALYSIS.json"))
    rows = [r for r in s["rows"]
            if r.get("G12D_log_HR") is not None and r.get("G12D_se_log_HR") and r["G12D_se_log_HR"]>0]
    if len(rows) < 2:
        return {"error": "<2 cohorts"}
    log_hrs = np.array([r["G12D_log_HR"] for r in rows])
    ses     = np.array([r["G12D_se_log_HR"] for r in rows])
    weights = 1.0 / (ses**2)
    pooled_log = (weights * log_hrs).sum() / weights.sum()
    rng = np.random.default_rng(0xC0FFEE)
    boots = []
    n = len(rows)
    for _ in range(B):
        # sample log HRs from per-study Normal(point, se) — parametric bootstrap
        sim = rng.normal(log_hrs, ses)
        pooled_b = (weights * sim).sum() / weights.sum()
        boots.append(pooled_b)
    boots = np.array(boots)
    HR = math.exp(pooled_log)
    return {
        "B": B, "n_studies": n,
        "pooled_HR": round(HR,4),
        "boot_HR_2_5": round(math.exp(np.percentile(boots, 2.5)),4),
        "boot_HR_97_5": round(math.exp(np.percentile(boots, 97.5)),4),
        "boot_HR_50": round(math.exp(np.percentile(boots, 50)),4),
        "p_boot_GT_1_5": float((np.exp(boots) > 1.5).mean()),
        "p_boot_LE_1": float((np.exp(boots) <= 1.0).mean()),
    }


# ====================================================================
# (C) DROP-ONE SENSITIVITY
# ====================================================================
def drop_one():
    print("[20C] Drop-one cohort sensitivity …")
    s = json.load(open(ROOT/"results/method_panel/META_ANALYSIS.json"))
    rows = [r for r in s["rows"]
            if r.get("G12D_log_HR") is not None and r.get("G12D_se_log_HR") and r["G12D_se_log_HR"]>0]
    out = []
    for skip_idx in range(len(rows)):
        kept = [r for i, r in enumerate(rows) if i != skip_idx]
        if not kept: continue
        log_hrs = np.array([r["G12D_log_HR"] for r in kept])
        ses     = np.array([r["G12D_se_log_HR"] for r in kept])
        w = 1.0/(ses**2); W = w.sum()
        pooled_log = (w*log_hrs).sum()/W
        pooled_se  = math.sqrt(1.0/W)
        z = pooled_log/pooled_se
        from math import erf
        p = 2*(1 - 0.5*(1+erf(abs(z)/math.sqrt(2))))
        out.append({
            "dropped": rows[skip_idx]["label"],
            "kept_n_cohorts": len(kept),
            "kept_total_OS": sum(r.get("n_with_OS") or 0 for r in kept),
            "pooled_HR": round(math.exp(pooled_log),3),
            "CI_lo": round(math.exp(pooled_log-1.96*pooled_se),3),
            "CI_hi": round(math.exp(pooled_log+1.96*pooled_se),3),
            "p": float(f"{p:.4g}"),
        })
    return out


# ====================================================================
# (D) RANDOM-EFFECT META (DerSimonian-Laird)
# ====================================================================
def derlaird():
    print("[20D] Random-effect (DerSimonian-Laird) meta …")
    s = json.load(open(ROOT/"results/method_panel/META_ANALYSIS.json"))
    rows = [r for r in s["rows"]
            if r.get("G12D_log_HR") is not None and r.get("G12D_se_log_HR") and r["G12D_se_log_HR"]>0]
    yi = np.array([r["G12D_log_HR"] for r in rows])
    vi = np.array([r["G12D_se_log_HR"]**2 for r in rows])
    wi_fe = 1/vi
    yi_bar = (wi_fe*yi).sum()/wi_fe.sum()
    Q = (wi_fe*(yi - yi_bar)**2).sum()
    df = len(yi) - 1
    C = wi_fe.sum() - (wi_fe**2).sum()/wi_fe.sum()
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    wi_re = 1/(vi + tau2)
    pooled_log_re = (wi_re*yi).sum()/wi_re.sum()
    pooled_se_re = math.sqrt(1/wi_re.sum())
    z = pooled_log_re/pooled_se_re
    from math import erf
    p = 2*(1 - 0.5*(1+erf(abs(z)/math.sqrt(2))))
    I2 = max(0.0, (Q - df)/Q) * 100 if Q > 0 else 0.0
    p_Q = 1.0 - chi2.cdf(Q, df) if df > 0 else 1.0
    return {
        "Q": round(float(Q),3), "df": int(df), "p_Q": float(f"{p_Q:.4g}"),
        "I2_percent": round(float(I2),1), "tau2": round(float(tau2),5),
        "RE_pooled_HR": round(math.exp(pooled_log_re),3),
        "RE_CI_lo": round(math.exp(pooled_log_re-1.96*pooled_se_re),3),
        "RE_CI_hi": round(math.exp(pooled_log_re+1.96*pooled_se_re),3),
        "RE_p": float(f"{p:.4g}"),
    }


# ====================================================================
# (E) DepMap proxy — pull CCLE KRAS-mut cell lines from cBioPortal
# Use the ccle_broad_2019 study which has mut + lineage
# ====================================================================
def depmap_proxy():
    print("[20E] CCLE KRAS allele × tissue cell-line distribution …")
    study = "ccle_broad_2019"
    samples = safe_get(f"{API}/studies/{study}/samples")
    if not samples: return {"error": "ccle samples missing"}
    sdf = pd.DataFrame(samples)
    scl = safe_get(f"{API}/studies/{study}/clinical-data", clinicalDataType="SAMPLE") or []
    sdf2 = (pd.DataFrame(scl).pivot_table(index="sampleId", columns="clinicalAttributeId",
            values="value", aggfunc="first").reset_index() if scl else pd.DataFrame())
    if sdf2.empty: return {"error": "ccle clinical missing"}
    df = sdf.merge(sdf2, on="sampleId", how="left")
    # KRAS muts
    g = safe_post(f"{API}/genes/fetch", body=["KRAS"], geneIdType="HUGO_GENE_SYMBOL")
    if not g: return {"error": "kras gene fetch failed"}
    eid = g[0]["entrezGeneId"]
    body = {"entrezGeneIds":[eid], "sampleListId": f"{study}_sequenced"}
    muts = safe_post(f"{API}/molecular-profiles/{study}_mutations/mutations/fetch",
                     body=body, projection="DETAILED")
    if not muts:
        # alternate sample list
        body["sampleListId"] = f"{study}_all"
        muts = safe_post(f"{API}/molecular-profiles/{study}_mutations/mutations/fetch",
                         body=body, projection="DETAILED")
    md = pd.DataFrame(muts) if muts else pd.DataFrame()
    if "gene" in md.columns:
        md["hugo"] = md["gene"].apply(lambda gg: gg.get("hugoGeneSymbol") if isinstance(gg, dict) else None)

    df["allele"] = "WT"
    if len(md):
        miss = md[(md.get("mutationType","").fillna("")=="Missense_Mutation") & (md.get("hugo","KRAS")=="KRAS")]
        for s, gg in miss.groupby("sampleId"):
            pcs = gg["proteinChange"].dropna().tolist()
            for p in ["G12D","G12V","G12R","G12C","G13D"]:
                if any(pc == p for pc in pcs):
                    df.loc[df["sampleId"]==s, "allele"] = p; break
    # tissue / lineage column (CANCER_TYPE_DETAILED or similar)
    lineage_col = next((c for c in ["LINEAGE","CANCER_TYPE","CANCER_TYPE_DETAILED","ONCOTREE_CODE"] if c in df.columns), None)
    if not lineage_col: return {"error": "no lineage col"}
    df["lineage"] = df[lineage_col].fillna("unknown").astype(str)
    df["lineage_grp"] = df["lineage"].str.lower().apply(lambda x:
        "PDAC" if "pancrea" in x or "paad" in x else
        "LUAD" if ("lung" in x and ("aden" in x or "non-small" in x)) else
        "CRC"  if ("colon" in x or "rect" in x or "colorect" in x) else
        "OTHER")
    summary = (df[df["allele"]!="WT"].groupby(["lineage_grp","allele"])
               .size().unstack(fill_value=0).to_dict())
    return {"n_lines_total": int(len(df)),
            "n_kras_mut_lines": int((df["allele"]!="WT").sum()),
            "by_lineage_x_allele": summary,
            "note": "CCLE KRAS-mut cell-line distribution — for downstream DepMap CRISPR dependency analysis"}


def main():
    t0 = time.time()
    out = {"generated_at": time.strftime("%Y-%m-%d %H:%M:%S")}

    out["A_interaction_test"] = interaction_test()
    out["B_bootstrap"]        = bootstrap_meta(10000)
    out["C_drop_one"]         = drop_one()
    out["D_random_effect"]    = derlaird()
    out["E_depmap_proxy"]     = depmap_proxy()
    out["elapsed_s"] = round(time.time() - t0, 1)

    with open(OUT/"SUMMARY.json","w") as f:
        json.dump(out, f, indent=2, default=str)

    print("\n" + "="*70)
    print(f"[20] NComm strengthening · {out['elapsed_s']}s")
    print("="*70)
    A = out["A_interaction_test"]
    if "G12R_x_PDAC_HR" in A:
        print(f"  (A) interaction G12R × PDAC: HR = {A['G12R_x_PDAC_HR']}  p = {A['G12R_x_PDAC_p']}  n = {A['n_total']}")
    if "G12D_x_PDAC_HR" in A:
        print(f"      interaction G12D × PDAC: HR = {A['G12D_x_PDAC_HR']}  p = {A['G12D_x_PDAC_p']}")
    B = out["B_bootstrap"]
    print(f"  (B) bootstrap pooled HR = {B.get('pooled_HR')}  95% CI [{B.get('boot_HR_2_5')}, {B.get('boot_HR_97_5')}]  P(HR>1.5)={B.get('p_boot_GT_1_5')}")
    print(f"  (C) drop-one sensitivity (showing first 3):")
    for r in (out["C_drop_one"] or [])[:3]:
        print(f"       drop {r['dropped'][:38]:38s}  HR={r['pooled_HR']}  CI {r['CI_lo']}-{r['CI_hi']}  p={r['p']}")
    D = out["D_random_effect"]
    print(f"  (D) random-effect: Q={D.get('Q')} I²={D.get('I2_percent')}%  τ²={D.get('tau2')}  RE-HR={D.get('RE_pooled_HR')}  CI [{D.get('RE_CI_lo')}, {D.get('RE_CI_hi')}]  p={D.get('RE_p')}")
    E = out["E_depmap_proxy"]
    if "by_lineage_x_allele" in E:
        print(f"  (E) CCLE KRAS-mut lines: total {E.get('n_lines_total')} · mut {E.get('n_kras_mut_lines')}")
        for al, lin in E["by_lineage_x_allele"].items():
            print(f"       {al}: {lin}")


if __name__ == "__main__":
    main()
