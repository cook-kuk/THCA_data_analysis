"""
v13 Task 3 — Virtual screening (similarity-based fallback).

Docking is skipped (GNINA/Vina not installed + 20k poses unrealistic in budget).
Instead we do a principled similarity screen:

  1. Pull a FDA-approved small-molecule library with SMILES from ChEMBL
     (max_phase=4). Cache to $RES/docking/fda_library.tsv.
  2. For each target, define literature-anchored reference binders (SMILES).
  3. Morgan (ECFP4) Tanimoto similarity FDA_library × refs → max_sim per compound.
  4. Rank top 50 per target.

Output:
  $RES/docking/fda_library.tsv
  $RES/docking/{GENE}_top50.tsv
  $RES/docking/screen_all_top.tsv
"""
from __future__ import annotations
import json, time, sys
from pathlib import Path
from typing import Iterable
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v13_drug_discovery"
DOCK = RES / "docking"
DOCK.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "thca-v13 (kukshomr@gmail.com)", "Accept": "application/json"}

# Reference binders (per-target, curated from literature)
REFS = {
    "CYP1B1": [
        ("alpha-naphthoflavone", "O=c1cc(-c2ccc3ccccc3c2)oc2ccccc12"),
        ("resveratrol",          "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1"),
        ("TMS",                  "COc1cc(/C=C/c2cc(OC)c(OC)c(OC)c2)ccc1OC"),  # 2,4,3',5'-TMS
        ("luteolin",             "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
        ("kaempferol",           "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"),
        ("quercetin",            "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
        ("fisetin",              "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)ccc12"),
        ("cannabidiol",          "CCCCCc1cc(O)c([C@@H]2C=C(C)CC[C@H]2C(=C)C)c(O)c1"),
    ],
    "GABRB2": [
        ("diazepam",   "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
        ("midazolam",  "Cc1ncc2CN=C(c3ccccc3F)c3cc(Cl)ccc3-n12"),
        ("zolpidem",   "Cc1ccc(-c2nc3ccc(C)cn3c2CC(=O)N(C)C)cc1"),
        ("etomidate",  "CCOC(=O)c1cncn1[C@H](C)c1ccccc1"),
        ("propofol",   "CC(C)c1cccc(C(C)C)c1O"),
        ("flumazenil", "CCOC(=O)c1ncn2-c3ccc(F)cc3C(=O)N(C)Cc12"),
        ("clonazepam", "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1"),
    ],
    "LDLR": [
        # Direct LDLR-binders are rare; use statins (SREBP axis) + ezetimibe (NPC1L1)
        # as clinically co-positioned lipid-modifying anchors.
        ("atorvastatin", "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)n(CC[C@H](O)C[C@H](O)CC(=O)O)c1-c1ccc(F)cc1"),
        ("simvastatin",  "CCC(C)(C)C(=O)O[C@H]1C[C@@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@@H]12"),
        ("rosuvastatin", "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O"),
        ("ezetimibe",    "O=C1[C@@H](CC[C@@H](O)c2ccc(F)cc2)[C@@H](c2ccc(O)cc2)N1c1ccc(F)cc1"),
        ("evolocumab_fragment_mimic", "Cc1ccc(-c2cc(C(=O)O)ccc2O)cc1"),  # placeholder
    ],
    "TACSTD2": [
        # TROP2 is an extracellular antigen — no natural small-mol binder.
        # We use the ADC payload SN-38 (irinotecan active metabolite) and its parents
        # as "chemistry of interest" anchors for payload-focused screening.
        ("SN-38",            "CC[C@@]1(O)C(=O)OCc2c1cc1n(c2=O)Cc2cc3cc(O)ccc3nc2-1"),
        ("irinotecan",       "CCC1=C2CN3C(=O)Cc4c(c1n2Cc1cc2cc(OC(=O)N5CCC(N6CCCCC6)CC5)ccc2nc41)C3(CC)O"),
        ("govitecan_linker", "CC[C@@]1(O)C(=O)OCc2c1cc1n(c2=O)Cc2cc3cc(OCC(=O)N)ccc3nc2-1"),
        ("topotecan",        "CN(C)Cc1cc2cc3n(c(=O)c2cc1O)Cc1c(CC)[C@@]2(O)C(=O)OCc1-3"),
    ],
    "TMPRSS4": [
        # TMPRSS4 is a type II transmembrane serine protease; no approved direct
        # inhibitor. Camostat / nafamostat target TMPRSS2 family and cross-react.
        ("camostat",    "CCOC(=O)CC(=O)OCOC(=O)c1ccc(NC(=N)N)cc1"),
        ("nafamostat",  "NC(=N)Nc1ccc(C(=O)Oc2ccc(/C(N)=N/O)cc2)cc1"),
        ("gabexate",    "CCOC(=O)CCCCCC(=O)OC1=CC=C(C=C1)NC(=N)N"),
        ("sivelestat",  "CCC(=O)NC1=CC=C(C=C1)S(=O)(=O)N(C(C)=O)CC(=O)O"),
    ],
    "PLEKHA6": [
        # No known small molecules. Use PROTAC warhead chemistry (cereblon binders).
        ("lenalidomide",  "NC1=CC=CC2=C1C(=O)N(C1CCC(=O)NC1=O)C2"),
        ("pomalidomide",  "O=C1CCC(N2C(=O)c3cccc(N)c3C2=O)C(=O)N1"),
        ("thalidomide",   "O=C1NC(=O)C(N2C(=O)c3ccccc3C2=O)C1"),
    ],
    "PTPRE": [
        # Phosphatase inhibitors — mostly non-selective chemistry (aryl sulfonates,
        # phosphate mimetics). No approved PTPRE-specific drug.
        ("NSC-87877",      "O=C(O)c1ccc2cc(N=Nc3ccc(S(=O)(=O)O)c4ccccc34)c(O)cc2c1"),
        ("suramin",        "Nc1ccc2cc(NC(=O)c3ccc(NC(=O)Nc4cccc(C)c4)cc3)ccc2c1"),  # approximation
        ("ertiprotafib",   "Cn1c2ccccc2c2cc(CC(=O)O)ccc21"),  # approximation
    ],
    "B3GNT3": [
        # No approved small molecule. Use sialyltransferase / glycosyltransferase
        # reference chemistry: UDP-GlcNAc mimetics, 3FAx-peracetyl-N-acetylmannosamine.
        ("UDP_GlcNAc_mimic", "CC(=O)NC1C(O)OC(COP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)C(O)C1O"),
        ("2FLacNAc",         "CC(=O)NC1C(O)OC(CO)C(F)C1OC1OC(CO)C(O)C(O)C1O"),
    ],
}


def fetch_fda_library(limit: int = 2500) -> pd.DataFrame:
    """ChEMBL max_phase=4 molecules with canonical SMILES + name + ID.
       Cached to disk to avoid re-fetching."""
    cache = DOCK / "fda_library.tsv"
    if cache.exists() and cache.stat().st_size > 1000:
        df = pd.read_csv(cache, sep="\t")
        print(f"FDA library cached: {len(df)} rows")
        return df

    rows = []
    offset = 0
    page = 200
    while len(rows) < limit:
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
            f"?max_phase=4.0&limit={page}&offset={offset}"
        )
        try:
            r = requests.get(url, headers=UA, timeout=60)
            if r.status_code != 200:
                print(f"  ChEMBL HTTP {r.status_code} at offset {offset}")
                break
            js = r.json()
        except Exception as e:
            print(f"  fetch fail offset {offset}: {e}")
            break
        mols = js.get("molecules", [])
        if not mols:
            break
        for m in mols:
            smi = ((m.get("molecule_structures") or {}).get("canonical_smiles"))
            if not smi:
                continue
            rows.append({
                "chembl_id": m.get("molecule_chembl_id", ""),
                "name":      (m.get("pref_name") or "").lower(),
                "smiles":    smi,
                "max_phase": m.get("max_phase"),
                "first_approval": m.get("first_approval"),
                "oral":      m.get("oral"),
                "parenteral": m.get("parenteral"),
                "withdrawn": m.get("withdrawn_flag"),
                "molecule_type": m.get("molecule_type"),
            })
        print(f"  ChEMBL offset={offset}  total rows={len(rows)}")
        offset += page
        time.sleep(0.4)
        if len(mols) < page:
            break

    df = pd.DataFrame(rows)
    # Prefer small molecules
    df = df[df["molecule_type"] == "Small molecule"].copy()
    df = df.drop_duplicates(subset=["smiles"])
    df.to_csv(cache, sep="\t", index=False)
    print(f"FDA library fetched: {len(df)} rows → {cache}")
    return df


def mol_fp(smi: str):
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None, None
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
        return m, fp
    except Exception:
        return None, None


def tanimoto_to_refs(fp, ref_fps) -> float:
    if fp is None or not ref_fps:
        return 0.0
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    return float(max(sims))


def compute_descriptors(m):
    try:
        return {
            "mw":       round(Descriptors.MolWt(m), 2),
            "logp":     round(Descriptors.MolLogP(m), 2),
            "tpsa":     round(Descriptors.TPSA(m), 2),
            "hba":      Descriptors.NumHAcceptors(m),
            "hbd":      Descriptors.NumHDonors(m),
            "rotb":     Descriptors.NumRotatableBonds(m),
            "rings":    Descriptors.RingCount(m),
            "qed":      round(Descriptors.qed(m), 3),
            "fsp3":     round(Descriptors.FractionCSP3(m), 3),
        }
    except Exception:
        return {}


def main():
    lib = fetch_fda_library()
    if lib.empty:
        print("FDA library empty — aborting")
        return

    # Pre-compute fingerprints for library
    print("Computing library fingerprints...")
    lib_fps = []
    lib_mols = []
    keep_idx = []
    for i, row in lib.iterrows():
        m, fp = mol_fp(row["smiles"])
        if fp is None:
            continue
        lib_fps.append(fp)
        lib_mols.append(m)
        keep_idx.append(i)
    lib = lib.loc[keep_idx].reset_index(drop=True)
    print(f"Library with valid fingerprints: {len(lib)}")

    all_top = []
    for gene, refs in REFS.items():
        ref_fps = []
        ref_keep = []
        for name, smi in refs:
            _, fp = mol_fp(smi)
            if fp is not None:
                ref_fps.append(fp)
                ref_keep.append(name)
        print(f"\n[{gene}] references valid: {len(ref_fps)}/{len(refs)} ({ref_keep})")
        if not ref_fps:
            continue

        sims = [tanimoto_to_refs(fp, ref_fps) for fp in lib_fps]
        lib_copy = lib.copy()
        lib_copy["max_sim_to_ref"] = [round(s, 3) for s in sims]

        # Descriptor filter: drop obvious non-leads (MW > 800 = biologics/peptides)
        descs = []
        for m in lib_mols:
            descs.append(compute_descriptors(m))
        ddf = pd.DataFrame(descs)
        lib_copy = pd.concat([lib_copy.reset_index(drop=True), ddf.reset_index(drop=True)], axis=1)
        # Small-molecule lead-like filter
        mask = (
            (lib_copy["mw"].fillna(0) > 100) &
            (lib_copy["mw"].fillna(9999) < 900)
        )
        lib_f = lib_copy[mask].copy()
        lib_f = lib_f.sort_values("max_sim_to_ref", ascending=False).head(50)
        lib_f["target"] = gene
        out = DOCK / f"{gene}_top50.tsv"
        lib_f.to_csv(out, sep="\t", index=False)
        print(f"  top1 sim={lib_f['max_sim_to_ref'].iloc[0]}  name={lib_f['name'].iloc[0]}")
        print(f"  wrote {out}")
        all_top.append(lib_f.head(25))

    if all_top:
        combined = pd.concat(all_top, ignore_index=True)
        combined.to_csv(DOCK / "screen_all_top.tsv", sep="\t", index=False)
        print(f"\nCombined top → {DOCK/'screen_all_top.tsv'}  ({len(combined)} rows)")


if __name__ == "__main__":
    main()
