"""
v13 Task 6 — ADMET prediction for top virtual-screen hits.

admet-ai not installed in this environment, so we use RDKit-native proxies:
  - Lipinski rule of 5 + Veber filter (oral bioavailability proxy)
  - QED (quantitative drug-likeness)
  - ESOL logS estimate (aqueous solubility)
  - BBB likelihood heuristic (MW/logP/TPSA thresholds)
  - hERG liability heuristic (basic amine + logP + ring count)
  - PAINS / Brenk structural alerts (via RDKit FilterCatalog)
  - Hepatotoxicity proxy (reactive groups + lipophilicity)

Assembles a composite "ADMET favorability" score for ranking.

Output: $RES/admet/top_candidates_admet.tsv
"""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, FilterCatalog

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v13_drug_discovery"
DOCK = RES / "docking"
ADMET = RES / "admet"
ADMET.mkdir(parents=True, exist_ok=True)


# Pre-build catalogs once
_PAINS = FilterCatalog.FilterCatalogParams()
_PAINS.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
PAINS = FilterCatalog.FilterCatalog(_PAINS)

_BRENK = FilterCatalog.FilterCatalogParams()
_BRENK.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
BRENK = FilterCatalog.FilterCatalog(_BRENK)


# SMARTS for common reactive/tox substructures as a last-resort hepatotox proxy
HEP_FLAGS = [
    ("nitro_arene",   "[N+](=O)[O-]"),
    ("thiophene",     "c1ccsc1"),
    ("furan",         "c1ccoc1"),
    ("aniline_free",  "[NX3H2][c]"),
    ("hydrazine",     "NN"),
    ("acyl_halide",   "[CX3](=O)[F,Cl,Br,I]"),
    ("epoxide",       "C1OC1"),
    ("aldehyde",      "[CX3H1](=O)"),
]
HEP_PAT = [(n, Chem.MolFromSmarts(s)) for n, s in HEP_FLAGS]


def esol_logS(m):
    """Simplified ESOL (Delaney 2004) — MW + logP + RotB + aromatic prop."""
    logp = Descriptors.MolLogP(m)
    mw = Descriptors.MolWt(m)
    rotb = Descriptors.NumRotatableBonds(m)
    ar = sum(1 for a in m.GetAtoms() if a.GetIsAromatic())
    ap = ar / max(1, m.GetNumHeavyAtoms())
    return 0.16 - 0.63*logp - 0.0062*mw + 0.066*rotb - 0.74*ap


def lipinski(m):
    return dict(
        mw_ok   = Descriptors.MolWt(m) <= 500,
        logp_ok = Descriptors.MolLogP(m) <= 5,
        hba_ok  = Descriptors.NumHAcceptors(m) <= 10,
        hbd_ok  = Descriptors.NumHDonors(m) <= 5,
    )


def veber(m):
    return dict(
        rotb_ok = Descriptors.NumRotatableBonds(m) <= 10,
        tpsa_ok = Descriptors.TPSA(m) <= 140,
    )


def bbb_score(m):
    """Heuristic: MW 150-450, logP 1-4, TPSA < 90, HBD < 3 → BBB-likely.
       Returns 0..1."""
    mw = Descriptors.MolWt(m)
    logp = Descriptors.MolLogP(m)
    tpsa = Descriptors.TPSA(m)
    hbd = Descriptors.NumHDonors(m)
    s = 0
    if 150 <= mw <= 450: s += 0.3
    if 1 <= logp <= 4:   s += 0.3
    if tpsa <= 90:       s += 0.25
    if hbd <= 3:         s += 0.15
    return round(s, 3)


def herg_risk(m):
    """Heuristic: basic amine + logP > 3 + aromatic rings >= 2 → hERG risk."""
    logp = Descriptors.MolLogP(m)
    bn = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N" and a.GetTotalNumHs() <= 2 and a.GetDegree() >= 2)
    has_basic = any(
        a.GetSymbol() == "N" and a.GetFormalCharge() == 0 and a.GetTotalNumHs() >= 0
        and not a.GetIsAromatic() and a.GetDegree() <= 3
        for a in m.GetAtoms()
    )
    n_aro_rings = Chem.rdMolDescriptors.CalcNumAromaticRings(m)
    s = 0
    if has_basic:       s += 0.4
    if logp > 3:        s += 0.3
    if n_aro_rings >= 2: s += 0.3
    return round(min(1.0, s), 3)


def hepatotox_flags(m):
    hits = [n for n, p in HEP_PAT if m.HasSubstructMatch(p)]
    return hits


def structural_alerts(m):
    pains_hits = [e.GetDescription() for e in PAINS.GetMatches(m)]
    brenk_hits = [e.GetDescription() for e in BRENK.GetMatches(m)]
    return pains_hits, brenk_hits


def compute_admet(smi: str) -> dict | None:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    try:
        out = {
            "mw":       round(Descriptors.MolWt(m), 2),
            "logp":     round(Descriptors.MolLogP(m), 2),
            "tpsa":     round(Descriptors.TPSA(m), 2),
            "hba":      Descriptors.NumHAcceptors(m),
            "hbd":      Descriptors.NumHDonors(m),
            "rotb":     Descriptors.NumRotatableBonds(m),
            "ar_rings": Chem.rdMolDescriptors.CalcNumAromaticRings(m),
            "rings":    Descriptors.RingCount(m),
            "qed":      round(Descriptors.qed(m), 3),
            "fsp3":     round(Descriptors.FractionCSP3(m), 3),
            "logS":     round(esol_logS(m), 2),
        }
        lip = lipinski(m); ver = veber(m)
        out["lipinski_pass"] = all(lip.values())
        out["veber_pass"] = all(ver.values())
        out["bbb_score"] = bbb_score(m)
        out["herg_risk"] = herg_risk(m)
        hep = hepatotox_flags(m)
        out["hepatotox_flags"] = ";".join(hep)
        pains, brenk = structural_alerts(m)
        out["pains_hits"] = ";".join(pains[:3])
        out["brenk_hits"] = ";".join(brenk[:3])
        out["pains_count"] = len(pains)
        out["brenk_count"] = len(brenk)

        # Composite ADMET favorability (0..1, higher better)
        fav = 0.0
        fav += 0.25 * out["qed"]
        fav += 0.15 * (1.0 if out["lipinski_pass"] else 0.5)
        fav += 0.10 * (1.0 if out["veber_pass"] else 0.5)
        fav += 0.15 * (1.0 - out["herg_risk"])
        fav += 0.10 * (1.0 if len(hep) == 0 else max(0.0, 1.0 - 0.3*len(hep)))
        fav += 0.10 * (1.0 if out["pains_count"] == 0 else 0.0)
        fav += 0.08 * (1.0 if out["brenk_count"] <= 2 else max(0.0, 1.0 - 0.1*out["brenk_count"]))
        logs_norm = max(0.0, min(1.0, (out["logS"] + 6) / 6))  # -6..0 → 0..1
        fav += 0.07 * logs_norm
        out["admet_favorability"] = round(min(1.0, fav), 3)
        return out
    except Exception as e:
        return None


def main():
    # gather all target top-50 files
    frames = []
    for f in sorted(DOCK.glob("*_top50.tsv")):
        df = pd.read_csv(f, sep="\t")
        frames.append(df)
    if not frames:
        print("No docking outputs found")
        return
    candidates = pd.concat(frames, ignore_index=True)
    print(f"Candidates: {len(candidates)} rows across targets")

    admet_rows = []
    for _, r in candidates.iterrows():
        smi = r.get("smiles")
        if not isinstance(smi, str) or not smi:
            continue
        d = compute_admet(smi)
        if d is None:
            continue
        d["chembl_id"] = r.get("chembl_id")
        d["name"] = r.get("name")
        d["target"] = r.get("target")
        d["smiles"] = smi
        d["max_sim_to_ref"] = r.get("max_sim_to_ref")
        d["max_phase"] = r.get("max_phase")
        d["first_approval"] = r.get("first_approval")
        d["oral"] = r.get("oral")
        d["parenteral"] = r.get("parenteral")
        d["withdrawn"] = r.get("withdrawn")
        admet_rows.append(d)

    out = pd.DataFrame(admet_rows)
    # drop duplicates by (target, chembl_id)
    out = out.drop_duplicates(subset=["target", "chembl_id"])
    out.to_csv(ADMET / "top_candidates_admet.tsv", sep="\t", index=False)
    print(f"Wrote {ADMET/'top_candidates_admet.tsv'}  ({len(out)} rows)")

    # Summary
    print("\nBy target — mean ADMET favorability (top 50):")
    print(out.groupby("target")["admet_favorability"].agg(["count", "mean", "max"]).round(3))

    # Flag concerning
    concerning = out[(out["herg_risk"] > 0.6) | (out["pains_count"] > 0) | (out["hepatotox_flags"] != "")]
    print(f"\nCandidates with at least one red flag: {len(concerning)}/{len(out)}")


if __name__ == "__main__":
    main()
