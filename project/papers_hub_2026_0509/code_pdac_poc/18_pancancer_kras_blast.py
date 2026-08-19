"""
Pan-cancer KRAS-allele blast — use all 8 CPUs.

Strategy:
  1. Discover ALL cBioPortal studies with KRAS mutations
  2. ProcessPoolExecutor — per study, fetch + Cox + per-allele HR
  3. Aggregate by cancer type → meta-analysis per cancer
  4. Cross-cancer KRAS allele biology comparison

Hits:
  - PDAC (pancreatic) — already done, use as anchor
  - CRC (coadread)
  - LUAD (lung adeno)
  - LUSC (lung squamous)
  - GIST + endometrial + bile duct (lower-frequency KRAS)

Output: /data/pdac_poc/results/pancancer_kras/
"""

import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import pandas as pd
import numpy as np
import requests
from lifelines import CoxPHFitter

API = "https://www.cbioportal.org/api"
OUT = Path("/data/pdac_poc/results/pancancer_kras")
OUT.mkdir(parents=True, exist_ok=True)

# Curated study list — TCGA Pan-Cancer Atlas + key MSK
STUDIES = [
    # PDAC anchor
    ("paad_tcga_pan_can_atlas_2018", "PDAC"),
    # CRC
    ("coadread_tcga_pan_can_atlas_2018", "CRC"),
    ("crc_msk_2018", "CRC"),
    # NSCLC
    ("luad_tcga_pan_can_atlas_2018", "LUAD"),
    ("lusc_tcga_pan_can_atlas_2018", "LUSC"),
    ("lung_msk_2017", "LUNG"),
    ("nsclc_pd1_msk_2018", "LUNG_ICI"),
    # Bile duct
    ("chol_tcga_pan_can_atlas_2018", "CHOL"),
    # Other KRAS-relevant
    ("ucec_tcga_pan_can_atlas_2018", "UCEC"),
    ("cesc_tcga_pan_can_atlas_2018", "CESC"),
    ("stad_tcga_pan_can_atlas_2018", "STAD"),
    # Aggregate MSK
    ("msk_impact_2017", "PAN_MSK_IMPACT"),
    # MSK-CHORD
    ("msk_chord_2024", "PAN_MSK_CHORD"),
    # MSK pan-cancer met tropism
    ("msk_met_2021", "PAN_MET"),
]


def safe_get(url, **kw):
    try:
        r = requests.get(url, params=kw, timeout=60); r.raise_for_status()
        return r.json()
    except Exception:
        return None


def safe_post(url, body, **kw):
    try:
        r = requests.post(url, json=body, params=kw, timeout=120,
                          headers={"Content-Type":"application/json"}); r.raise_for_status()
        return r.json()
    except Exception:
        return None


def best_os_cols(clin_cols):
    pairs = [("OS_MONTHS","OS_STATUS"),("OS_MONTHS","STATUS"),("OS_MONTHS","VITAL_STATUS"),
             ("DAYS_TO_LAST_FOLLOWUP","STATUS"),("DAYS_TO_LAST_FOLLOWUP","VITAL_STATUS"),
             ("OVERALL_SURVIVAL_MONTHS","OVERALL_SURVIVAL_STATUS")]
    cols = set(clin_cols)
    for t,s in pairs:
        if t in cols and s in cols: return t, s
    return None, None


def parse_event(v):
    if isinstance(v,str):
        v=v.lower()
        if v.startswith("1") or v in ("dead","deceased","yes"): return 1
        if v.startswith("0") or v in ("alive","living","no"): return 0
    try: return int(float(v))
    except: return None


def per_study(study_label):
    study, label = study_label
    out = {"study": study, "cancer": label}

    samples = safe_get(f"{API}/studies/{study}/samples")
    if not samples:
        out["status"] = "no samples"; return out
    sdf = pd.DataFrame(samples)
    sample_ids = sdf["sampleId"].tolist()
    out["n_samples"] = len(sample_ids)

    pcl = safe_get(f"{API}/studies/{study}/clinical-data", clinicalDataType="PATIENT")
    scl = safe_get(f"{API}/studies/{study}/clinical-data", clinicalDataType="SAMPLE")
    pdf = (pd.DataFrame(pcl).pivot_table(index="patientId", columns="clinicalAttributeId",
            values="value", aggfunc="first").reset_index() if pcl else pd.DataFrame())
    sdf2 = (pd.DataFrame(scl).pivot_table(index="sampleId", columns="clinicalAttributeId",
            values="value", aggfunc="first").reset_index() if scl else pd.DataFrame())
    if not sdf2.empty: sdf = sdf.merge(sdf2, on="sampleId", how="left")
    if not pdf.empty:  sdf = sdf.merge(pdf, on="patientId", how="left")
    clin = sdf

    # KRAS muts
    g = safe_post(f"{API}/genes/fetch", body=["KRAS"], geneIdType="HUGO_GENE_SYMBOL")
    if not g:
        out["status"] = "kras gene fetch failed"; return out
    eid = g[0]["entrezGeneId"]
    body = {"entrezGeneIds":[eid], "sampleListId": f"{study}_sequenced"}
    muts = safe_post(f"{API}/molecular-profiles/{study}_mutations/mutations/fetch",
                     body=body, projection="DETAILED")
    if not muts:
        # try alternative profile id naming
        for alt in [f"{study}_kras_mut", f"{study}_mutation"]:
            muts = safe_post(f"{API}/molecular-profiles/{alt}/mutations/fetch",
                              body=body, projection="DETAILED")
            if muts: break
    if not muts:
        out["status"] = "no kras mut profile"; return out

    # allele assign
    md = pd.DataFrame(muts)
    if "gene" in md.columns:
        md["hugo"] = md["gene"].apply(lambda gg: gg.get("hugoGeneSymbol") if isinstance(gg, dict) else None)
    md = md[(md.get("mutationType","").fillna("") == "Missense_Mutation") &
            (md.get("hugo","KRAS") == "KRAS")]
    allele = pd.Series("WT", index=sample_ids)
    for s, gg in md.groupby("sampleId"):
        if s not in allele.index: continue
        pcs = gg["proteinChange"].dropna().tolist()
        for p in ["G12D","G12V","G12R","G12C","G13D","Q61H","Q61R","Q61L"]:
            if any(pc == p for pc in pcs): allele.loc[s] = p; break
        else:
            if any(pc.startswith("G12") for pc in pcs): allele.loc[s] = "G12_other"
            elif any(pc.startswith("Q61") for pc in pcs): allele.loc[s] = "Q61"
            elif pcs: allele.loc[s] = "KRAS_other"
    out["allele_counts"] = allele.value_counts().to_dict()

    # Cox
    t_col, s_col = best_os_cols(clin.columns)
    if not t_col:
        out["status"] = "no OS"; return out
    df = clin.set_index("sampleId").copy()
    df["allele"] = df.index.map(allele)
    df["t"] = pd.to_numeric(df[t_col], errors="coerce")
    df["e"] = df[s_col].apply(parse_event)
    if t_col.endswith("DAYS") or t_col.endswith("FOLLOWUP"):
        df["t"] = df["t"]/30.4375
    df = df[(df["t"]>0) & df["e"].notna()].dropna(subset=["t","e"])
    if len(df) < 30:
        out["status"] = f"n_OS={len(df)} <30"; return out
    out["n_with_OS"] = len(df)

    for al in ["G12D","G12V","G12R","G12C","G13D","Q61H","Q61R"]:
        df[al] = (df["allele"]==al).astype(int)
    cov_cols = [al for al in ["G12D","G12V","G12R","G12C","G13D"]
                if df[al].sum() >= 5 and df[al].sum() <= len(df)-5]
    if not cov_cols:
        out["status"] = "no allele with n>=5"; return out
    try:
        cph = CoxPHFitter(penalizer=0.05).fit(df[["t","e"]+cov_cols],
                                               duration_col="t", event_col="e")
        s = cph.summary[["coef","exp(coef)","p","exp(coef) lower 95%","exp(coef) upper 95%"]].round(4)
        out["cox"] = s.reset_index().to_dict(orient="records")
        out["concordance"] = float(cph.concordance_index_)
        for r in out["cox"]:
            al = r["covariate"]
            out[f"{al}_HR"] = r["exp(coef)"]
            out[f"{al}_p"]  = r["p"]
            if r["exp(coef) upper 95%"]>0 and r["exp(coef) lower 95%"]>0:
                out[f"{al}_se_log"] = (math.log(r["exp(coef) upper 95%"])
                                        - math.log(r["exp(coef) lower 95%"]))/3.92
                out[f"{al}_log_HR"] = math.log(r["exp(coef)"])
    except Exception as e:
        out["cox_error"] = str(e)
    return out


def fixed_meta(rows, allele):
    obs = [(r[f"{allele}_log_HR"], r[f"{allele}_se_log"], r["cancer"], r["n_with_OS"])
           for r in rows
           if r.get(f"{allele}_log_HR") is not None
           and r.get(f"{allele}_se_log") and r[f"{allele}_se_log"]>0]
    if not obs: return None
    w = [1/(s*s) for _,s,*_ in obs]; W = sum(w)
    pooled_log = sum(wi*li for wi,(li,*_) in zip(w,obs))/W
    pooled_se = math.sqrt(1/W)
    HR = math.exp(pooled_log)
    z = pooled_log/pooled_se
    from math import erf
    p = 2*(1 - 0.5*(1+erf(abs(z)/math.sqrt(2))))
    return {"allele": allele, "pooled_HR": round(HR,3),
            "CI_lo": round(math.exp(pooled_log-1.96*pooled_se),3),
            "CI_hi": round(math.exp(pooled_log+1.96*pooled_se),3),
            "p_value": float(f"{p:.4g}"), "n_studies": len(obs),
            "n_total_OS": sum(o[3] for o in obs),
            "studies": [o[2] for o in obs]}


def main():
    t0 = time.time()
    print(f"[18] launching {len(STUDIES)} cohorts on {os.cpu_count()} CPUs")
    rows = []
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        futs = {ex.submit(per_study, s): s for s in STUDIES}
        for fut in as_completed(futs):
            s = futs[fut]
            try:
                r = fut.result()
                rows.append(r)
                hrs = []
                for al in ["G12D","G12V","G12R","G12C","G13D"]:
                    if r.get(f"{al}_HR") is not None:
                        hrs.append(f"{al} HR={r[f'{al}_HR']:.2f} p={r[f'{al}_p']:.3g} n_alleled={r['allele_counts'].get(al,0)}")
                print(f"   ✓ {s[1]:8s} {s[0]:42s} n={r.get('n_samples')} OS={r.get('n_with_OS','—')} :: {' · '.join(hrs) if hrs else r.get('status','-')}")
            except Exception as e:
                rows.append({"study":s[0],"cancer":s[1],"error":str(e)})
                print(f"   ✗ {s[1]:8s} {s[0]:42s}  {e}")

    # cancer-stratified meta
    metas = {}
    rows_with_data = [r for r in rows if r.get("n_with_OS")]
    cancers = sorted({r["cancer"] for r in rows_with_data})
    for c in cancers:
        sub = [r for r in rows_with_data if r["cancer"]==c]
        for al in ["G12D","G12V","G12R","G12C","G13D"]:
            m = fixed_meta(sub, al)
            if m: metas[f"{c}_{al}"] = m

    # global pan-cancer meta (excluding controls)
    global_pan = {}
    for al in ["G12D","G12V","G12R","G12C","G13D"]:
        m = fixed_meta([r for r in rows_with_data if r["cancer"] not in ("CHOL",)], al)
        if m: global_pan[al] = m

    out = {"generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
           "elapsed_s": round(time.time()-t0,1),
           "n_studies_attempted": len(STUDIES),
           "n_studies_with_OS": len(rows_with_data),
           "per_cancer_meta": metas,
           "pan_cancer_meta": global_pan,
           "rows": rows}
    with open(OUT/"PANCANCER_META.json","w") as f:
        json.dump(out,f,indent=2,default=str)

    print("\n" + "="*70)
    print(f"[18] PAN-CANCER KRAS BLAST  {round(time.time()-t0,1)}s on {os.cpu_count()} cores")
    print("="*70)
    print("\n  cancer-stratified KRAS allele meta-HRs (CI · p · k · n):")
    for k, m in sorted(metas.items()):
        print(f"   {k:25s}  HR={m['pooled_HR']}  CI {m['CI_lo']}-{m['CI_hi']}  p={m['p_value']}  k={m['n_studies']} n={m['n_total_OS']}")
    print("\n  PAN-CANCER pooled (excluding controls):")
    for al, m in global_pan.items():
        print(f"   {al:8s}  HR={m['pooled_HR']}  CI {m['CI_lo']}-{m['CI_hi']}  p={m['p_value']}  k={m['n_studies']} n={m['n_total_OS']}")


if __name__ == "__main__":
    main()
