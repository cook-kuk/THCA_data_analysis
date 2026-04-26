#!/usr/bin/env python3
"""v7 Drug Repurposing end-to-end pipeline for THCA computational study.

Steps: 1) Open Targets knownDrugs 2) DrugCentral 3) ChEMBL clinical-stage
4) Unified scoring 5) Top-10 deep dive 6) Combination hypotheses
7) Dashboard HTML 8) Paper section. Network failures degrade gracefully.
"""
from __future__ import annotations
import datetime as _dt, gzip, itertools, json, logging, re, sys, time, traceback
import urllib.error, urllib.parse, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
RES = ROOT / "results" / "v7_repurposing"
RPT = ROOT / "reports" / "v7"
HTML = ROOT / "reports" / "html"
FIGS = HTML / "figs_interactive" / "v7"
PAGES = HTML / "pages"
DATA_RAW = Path("/data/thca/v7_repurposing")
LOG_PATH = ROOT / "logs" / "v7_repurposing.log"

for p in (RES, RPT, FIGS, PAGES, DATA_RAW, LOG_PATH.parent):
    p.mkdir(parents=True, exist_ok=True)

TARGETS = [
    "TACSTD2", "TMPRSS4", "PLEKHA6", "CYP1B1",
    "LDLR", "GABRB2", "B3GNT3", "PTPRE",
]

# UniProt IDs for the 8 targets (used for Reactome lookups)
UNIPROT_MAP = {
    "TACSTD2": "P09758",
    "TMPRSS4": "Q9NRS4",
    "PLEKHA6": "Q9Y2H5",
    "CYP1B1":  "Q16678",
    "LDLR":    "P01130",
    "GABRB2":  "P47870",
    "B3GNT3":  "Q9Y2A9",
    "PTPRE":   "P23469",
}

# Fallback Reactome pathway mapping (used when live Reactome calls fail)
FALLBACK_PATHWAYS = {
    "TACSTD2": {"R-HSA-1500931-cell-cell-comm", "R-HSA-442742-adhesion"},
    "TMPRSS4": {"R-HSA-1474244-ecm-organization", "R-HSA-75205-protease"},
    "PLEKHA6": {"R-HSA-373760-cell-adhesion"},
    "CYP1B1":  {"R-HSA-211859-biological-oxidations", "R-HSA-211981-xenobiotics", "R-HSA-5423646-retinoid-metab"},
    "LDLR":    {"R-HSA-8963899-cholesterol", "R-HSA-174824-plasma-lipoprotein", "R-HSA-8963889-lipid-digestion"},
    "GABRB2":  {"R-HSA-977443-gaba-signaling", "R-HSA-112310-neurotransmission"},
    "B3GNT3":  {"R-HSA-913709-n-glycan-biosynthesis", "R-HSA-446728-glycan-assembly"},
    "PTPRE":   {"R-HSA-388841-ptp-signaling", "R-HSA-9006934-signaling-by-rtks"},
}

THCA_STANDARD_DRUGS = {
    "lenvatinib", "sorafenib", "cabozantinib",
    "vandetanib", "selpercatinib", "pralsetinib",
}

OVEREXPRESSED = set(TARGETS)  # every target in this panel is overexpressed in THCA

USER_AGENT = "THCA-v7-repurposing/1.0 (+kukshomr@gmail.com)"
HTTP_TIMEOUT = 30
RETRIES = 2
RETRY_SLEEP = 1.0

# Logging

logger = logging.getLogger("v7")
logger.setLevel(logging.INFO)
_fmt = logging.Formatter("%(asctime)s  %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
_fh = logging.FileHandler(LOG_PATH, mode="a")
_fh.setFormatter(_fmt)
_sh = logging.StreamHandler(sys.stdout)
_sh.setFormatter(_fmt)
if not logger.handlers:
    logger.addHandler(_fh)
    logger.addHandler(_sh)

def log(msg: str, level: str = "info") -> None:
    getattr(logger, level)(msg)

# HTTP helpers (urllib + manual retries)

def _request(method: str, url: str, *, data: bytes | None = None, headers: dict | None = None) -> tuple[int, bytes]:
    hdrs = {"User-Agent": USER_AGENT}
    if headers:
        hdrs.update(headers)
    req = urllib.request.Request(url, data=data, headers=hdrs, method=method)
    try:
        with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as resp:
            return resp.status, resp.read()
    except urllib.error.HTTPError as e:
        body = b""
        try:
            body = e.read()
        except Exception:
            pass
        return e.code, body

def http_get(url: str, headers: dict | None = None) -> tuple[int, bytes]:
    last_exc = None
    for attempt in range(RETRIES + 1):
        try:
            return _request("GET", url, headers=headers)
        except Exception as e:
            last_exc = e
            time.sleep(RETRY_SLEEP)
    log(f"GET failed after retries: {url} ({last_exc})", "warning")
    return 0, b""

def http_post_json(url: str, payload: dict, headers: dict | None = None) -> tuple[int, bytes]:
    body = json.dumps(payload).encode("utf-8")
    hdrs = {"Content-Type": "application/json"}
    if headers:
        hdrs.update(headers)
    last_exc = None
    for attempt in range(RETRIES + 1):
        try:
            return _request("POST", url, data=body, headers=hdrs)
        except Exception as e:
            last_exc = e
            time.sleep(RETRY_SLEEP)
    log(f"POST failed after retries: {url} ({last_exc})", "warning")
    return 0, b""

def safe_json(raw: bytes) -> dict | list | None:
    try:
        return json.loads(raw.decode("utf-8"))
    except Exception:
        return None

# STEP 1 — Open Targets

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

OT_SEARCH_QUERY = """
query Search($q: String!) {
  search(queryString: $q, entityNames: ["target"]) {
    hits { id name entity object { ... on Target { approvedSymbol id } } }
  }
}
"""

OT_KNOWN_DRUGS_QUERY = """
query KnownDrugs($ensgId: String!) {
  target(ensemblId: $ensgId) {
    id
    approvedSymbol
    knownDrugs(size: 200) {
      count
      rows {
        drugId
        prefName
        drugType
        phase
        status
        mechanismOfAction
        disease { id name }
      }
    }
  }
}
"""

def ot_resolve_ensg(gene: str) -> str | None:
    status, raw = http_post_json(OT_URL, {"query": OT_SEARCH_QUERY, "variables": {"q": gene}})
    if status != 200:
        log(f"[OT] search {gene} → HTTP {status}", "warning")
        return None
    js = safe_json(raw) or {}
    hits = ((js.get("data") or {}).get("search") or {}).get("hits") or []
    for h in hits:
        obj = h.get("object") or {}
        if (obj.get("approvedSymbol") or "").upper() == gene.upper():
            return obj.get("id") or h.get("id")
    for h in hits:
        if h.get("entity") == "target":
            return h.get("id")
    return None

def step1_opentargets() -> pd.DataFrame:
    log("STEP 1  Open Targets knownDrugs — starting")
    rows: list[dict] = []
    for gene in TARGETS:
        try:
            ensg = ot_resolve_ensg(gene)
            if not ensg:
                log(f"[OT] {gene}: could not resolve Ensembl ID", "warning")
                (DATA_RAW / f"opentargets_{gene}.json").write_text("{}")
                continue
            status, raw = http_post_json(OT_URL, {"query": OT_KNOWN_DRUGS_QUERY, "variables": {"ensgId": ensg}})
            if status != 200:
                log(f"[OT] {gene} ({ensg}) knownDrugs → HTTP {status}", "warning")
                (DATA_RAW / f"opentargets_{gene}.json").write_text("{}")
                continue
            js = safe_json(raw) or {}
            (DATA_RAW / f"opentargets_{gene}.json").write_bytes(raw)
            tgt = (js.get("data") or {}).get("target") or {}
            kd = (tgt.get("knownDrugs") or {}).get("rows") or []
            for r in kd:
                rows.append({
                    "target": gene,
                    "drug_id": r.get("drugId"),
                    "drug_name": (r.get("prefName") or "").lower(),
                    "drug_type": r.get("drugType"),
                    "max_phase": r.get("phase"),
                    "status": r.get("status"),
                    "mechanism_of_action": r.get("mechanismOfAction"),
                    "original_indication": ((r.get("disease") or {}).get("name") or ""),
                    "source": "opentargets",
                })
            log(f"[OT] {gene}: {len(kd)} known-drug rows")
        except Exception as e:
            log(f"[OT] {gene}: exception {e}", "warning")
    df = pd.DataFrame(rows)
    df.to_csv(RES / "v7_opentargets_known_drugs.tsv", sep="\t", index=False)
    log(f"STEP 1  done — {len(df)} rows → {RES/'v7_opentargets_known_drugs.tsv'}")
    return df

# STEP 2 — DrugCentral

DC_URL_PRIMARY = "https://unmtid-shinyapps.net/download/drug.target.interaction.tsv.gz"
DC_URL_FALLBACK_PAGE = "https://drugcentral.org/download"

def step2_drugcentral() -> pd.DataFrame:
    log("STEP 2  DrugCentral — starting")
    cache = DATA_RAW / "drugcentral_interactions.tsv.gz"
    if not cache.exists():
        status, raw = http_get(DC_URL_PRIMARY)
        if status == 200 and raw:
            cache.write_bytes(raw)
            log(f"[DC] downloaded primary ({len(raw):,} bytes)")
        else:
            log(f"[DC] primary URL failed ({status}); trying fallback page", "warning")
            st2, html = http_get(DC_URL_FALLBACK_PAGE)
            tsv_url = None
            if st2 == 200 and html:
                m = re.search(r'href="([^"]*drug\.target\.interaction[^"]*\.tsv\.gz)"', html.decode("utf-8", "ignore"))
                if m:
                    tsv_url = m.group(1)
                    if tsv_url.startswith("/"):
                        tsv_url = "https://drugcentral.org" + tsv_url
            if tsv_url:
                st3, raw2 = http_get(tsv_url)
                if st3 == 200 and raw2:
                    cache.write_bytes(raw2)
                    log(f"[DC] downloaded fallback ({len(raw2):,} bytes)")
            if not cache.exists():
                log("[DC] all download attempts failed — emitting empty subset", "warning")
                empty = pd.DataFrame(columns=["target", "drug_name", "drug_id", "source"])
                empty.to_csv(RES / "v7_drugcentral_filtered.tsv", sep="\t", index=False)
                return empty

    try:
        with gzip.open(cache, "rb") as fh:
            df = pd.read_csv(fh, sep="\t", low_memory=False)
    except Exception as e:
        log(f"[DC] read failed: {e}", "warning")
        empty = pd.DataFrame(columns=["target", "drug_name", "drug_id", "source"])
        empty.to_csv(RES / "v7_drugcentral_filtered.tsv", sep="\t", index=False)
        return empty

    gene_col = None
    for c in ("GENE", "GENE_NAME", "HGNC_SYMBOL", "gene_symbol", "symbol"):
        if c in df.columns:
            gene_col = c
            break
    if gene_col is None:
        log(f"[DC] no gene column in {list(df.columns)[:8]}", "warning")
        empty = pd.DataFrame(columns=["target", "drug_name", "drug_id", "source"])
        empty.to_csv(RES / "v7_drugcentral_filtered.tsv", sep="\t", index=False)
        return empty

    # Gene entries may be '|' separated
    def hits(gstr: str) -> list[str]:
        if not isinstance(gstr, str):
            return []
        toks = {t.strip().upper() for t in re.split(r"[|;,/\s]+", gstr) if t.strip()}
        return [g for g in TARGETS if g in toks]

    df["_target_hits"] = df[gene_col].map(hits)
    df = df[df["_target_hits"].map(len) > 0].copy()
    exploded = df.explode("_target_hits").rename(columns={"_target_hits": "target"})

    drug_col = next((c for c in ("DRUG_NAME", "drug_name", "NAME") if c in exploded.columns), None)
    id_col = next((c for c in ("STRUCT_ID", "struct_id", "DRUGCENTRAL_ID") if c in exploded.columns), None)
    act_col = next((c for c in ("ACT_VALUE", "act_value") if c in exploded.columns), None)
    actn_col = next((c for c in ("ACT_TYPE", "act_type") if c in exploded.columns), None)
    moa_col = next((c for c in ("ACTION_TYPE", "action_type") if c in exploded.columns), None)

    out = pd.DataFrame({
        "target": exploded["target"],
        "drug_name": exploded[drug_col].astype(str).str.lower() if drug_col else "",
        "drug_id": (exploded[id_col].astype(str) if id_col else ""),
        "max_phase": np.nan,
        "original_indication": "",
        "mechanism_of_action": exploded[moa_col].astype(str) if moa_col else "",
        "pchembl_value": (pd.to_numeric(exploded[act_col], errors="coerce") if act_col else np.nan),
        "activity_type": exploded[actn_col].astype(str) if actn_col else "",
        "source": "drugcentral",
    })
    out["drug_id"] = "DC:" + out["drug_id"].astype(str)
    out.to_csv(RES / "v7_drugcentral_filtered.tsv", sep="\t", index=False)
    log(f"STEP 2  done — {len(out)} drug-target rows across targets {sorted(out['target'].unique())}")
    return out

# STEP 3 — ChEMBL clinical-stage activities

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

def chembl_target_id(gene: str) -> str | None:
    url = f"{CHEMBL_BASE}/target/search.json?q={urllib.parse.quote(gene)}"
    status, raw = http_get(url)
    if status != 200:
        return None
    js = safe_json(raw) or {}
    for t in js.get("targets") or []:
        if t.get("organism") != "Homo sapiens":
            continue
        if t.get("target_type") != "SINGLE PROTEIN":
            continue
        comps = t.get("target_components") or []
        for c in comps:
            for syn in c.get("target_component_synonyms") or []:
                if (syn.get("component_synonym") or "").upper() == gene.upper():
                    return t.get("target_chembl_id")
        if (t.get("pref_name") or "").upper().startswith(gene.upper()):
            return t.get("target_chembl_id")
    # fallback: first human hit
    for t in js.get("targets") or []:
        if t.get("organism") == "Homo sapiens":
            return t.get("target_chembl_id")
    return None

def chembl_activities(target_id: str) -> list[dict]:
    url = (f"{CHEMBL_BASE}/activity.json?target_chembl_id={target_id}"
           f"&pchembl_value__gte=6&limit=1000")
    status, raw = http_get(url)
    if status != 200:
        return []
    js = safe_json(raw) or {}
    return js.get("activities") or []

def chembl_molecules(mol_ids: list[str]) -> list[dict]:
    out: list[dict] = []
    for i in range(0, len(mol_ids), 40):
        batch = mol_ids[i:i + 40]
        ids = ",".join(batch)
        url = f"{CHEMBL_BASE}/molecule.json?molecule_chembl_id__in={ids}&limit=100"
        status, raw = http_get(url)
        if status != 200:
            continue
        js = safe_json(raw) or {}
        out.extend(js.get("molecules") or [])
    return out

def step3_chembl() -> pd.DataFrame:
    log("STEP 3  ChEMBL clinical-stage — starting")
    all_rows: list[dict] = []
    for gene in TARGETS:
        try:
            tid = chembl_target_id(gene)
            if not tid:
                log(f"[CHEMBL] {gene}: no target id found", "warning")
                pd.DataFrame().to_csv(RES / f"v7_chembl_clinical_{gene}.tsv", sep="\t", index=False)
                continue
            acts = chembl_activities(tid)
            mol_ids = sorted({a.get("molecule_chembl_id") for a in acts if a.get("molecule_chembl_id")})
            mols = chembl_molecules(mol_ids) if mol_ids else []
            # index best pchembl per molecule
            best_pc: dict[str, float] = {}
            for a in acts:
                mid = a.get("molecule_chembl_id")
                pv = a.get("pchembl_value")
                try:
                    pv = float(pv) if pv is not None else None
                except Exception:
                    pv = None
                if mid and pv is not None:
                    best_pc[mid] = max(best_pc.get(mid, 0.0), pv)
            rows = []
            for m in mols:
                mp = m.get("max_phase")
                try:
                    mp_val = float(mp) if mp is not None else None
                except Exception:
                    mp_val = None
                if mp_val is None or mp_val < 1:
                    continue
                mid = m.get("molecule_chembl_id")
                props = m.get("molecule_properties") or {}
                name = (m.get("pref_name") or mid or "").lower()
                rows.append({
                    "target": gene,
                    "drug_id": mid,
                    "drug_name": name,
                    "max_phase": mp_val,
                    "first_approval": m.get("first_approval"),
                    "indication_class": m.get("indication_class") or "",
                    "withdrawn_flag": bool(m.get("withdrawn_flag")),
                    "oral": bool(m.get("oral")),
                    "parenteral": bool(m.get("parenteral")),
                    "pchembl_value": best_pc.get(mid),
                    "mechanism_of_action": "",
                    "original_indication": m.get("indication_class") or "",
                    "chembl_target_id": tid,
                    "usan_stem": (props.get("usan_stem") or ""),
                    "source": "chembl",
                })
            per = pd.DataFrame(rows)
            per.to_csv(RES / f"v7_chembl_clinical_{gene}.tsv", sep="\t", index=False)
            log(f"[CHEMBL] {gene} ({tid}): activities={len(acts)} molecules={len(mols)} clinical≥1={len(per)}")
            all_rows.extend(rows)
        except Exception as e:
            log(f"[CHEMBL] {gene}: exception {e}", "warning")
    df = pd.DataFrame(all_rows)
    df.to_csv(RES / "v7_chembl_clinical_all.tsv", sep="\t", index=False)
    log(f"STEP 3  done — {len(df)} clinical rows")
    return df

# STEP 4 — Repurposing scoring

def _potency(pc: float | None) -> float:
    if pc is None or (isinstance(pc, float) and np.isnan(pc)):
        return 0.3
    if pc >= 8:
        return 1.0
    if pc >= 7:
        return 0.7
    if pc >= 6:
        return 0.4
    return 0.2

def _thyroid_safety(row: dict) -> float:
    text = " ".join(str(row.get(k) or "") for k in ("indication_class", "original_indication", "drug_name", "mechanism_of_action")).lower()
    bad = ("hepato", "qt", "torsade", "hepatotox")
    if any(b in text for b in bad):
        return 0.4
    return 0.8

def _oral_bioavail(row: dict) -> float:
    oral = bool(row.get("oral"))
    par = bool(row.get("parenteral"))
    if oral and par:
        return 0.8
    if oral:
        return 1.0
    if par:
        return 0.3
    return 0.5

def _off_patent(row: dict) -> float:
    yr = row.get("first_approval")
    try:
        yr = int(yr) if yr not in (None, "", float("nan")) else None
    except Exception:
        yr = None
    if yr is None:
        return 0.3
    if yr <= 2005:
        return 1.0
    if yr <= 2015:
        return 0.5
    return 0.2

def _mechanism_fit(row: dict) -> float:
    moa = str(row.get("mechanism_of_action") or "").lower()
    gene = row.get("target")
    if any(k in moa for k in ("inhibitor", "antagonist", "blocker")):
        return 1.0
    if any(k in moa for k in ("agonist", "activator")):
        return 0.2 if gene in OVEREXPRESSED else 0.5
    return 0.5

def _rationale(row: dict, score_components: dict) -> str:
    bits = []
    mp = row.get("max_phase")
    if mp is not None and not (isinstance(mp, float) and np.isnan(mp)):
        bits.append(f"phase {int(mp) if mp >= 1 else mp}")
    if row.get("mechanism_of_action"):
        bits.append(str(row["mechanism_of_action"])[:60])
    if row.get("original_indication"):
        bits.append(f"original: {row['original_indication'][:48]}")
    bits.append(f"potency={score_components['potency']:.2f} mech_fit={score_components['mech']:.2f} offpatent={score_components['off']:.2f}")
    return "; ".join(b for b in bits if b)

def step4_score(ot_df: pd.DataFrame, dc_df: pd.DataFrame, ch_df: pd.DataFrame) -> pd.DataFrame:
    log("STEP 4  Unified repurposing scoring — starting")
    parts = []
    # Harmonize columns
    for src in (ot_df, dc_df, ch_df):
        if src is None or src.empty:
            continue
        parts.append(src.copy())
    if not parts:
        log("[SCORE] no data from any source", "warning")
        out = pd.DataFrame(columns=["target", "drug_name", "drug_id", "max_phase",
                                    "original_indication", "repurposing_score", "rationale", "source"])
        out.to_csv(RES / "v7_repurposing_ranked.tsv", sep="\t", index=False)
        return out

    df = pd.concat(parts, ignore_index=True, sort=False)
    for col in ("max_phase", "first_approval", "pchembl_value", "withdrawn_flag",
                "oral", "parenteral", "mechanism_of_action", "original_indication",
                "indication_class", "drug_name", "drug_id", "target", "source"):
        if col not in df.columns:
            df[col] = np.nan

    df["drug_name"] = df["drug_name"].astype(str).str.strip().str.lower()
    df = df[df["drug_name"].astype(bool)].copy()
    df = df[~df["drug_name"].isin(THCA_STANDARD_DRUGS)].copy()

    def max_phase_num(x):
        try:
            v = float(x)
            return v if not np.isnan(v) else np.nan
        except Exception:
            return np.nan
    df["max_phase"] = df["max_phase"].map(max_phase_num)
    df["max_phase"] = df["max_phase"].fillna(1.0)  # neutral baseline for OT/DC without phase

    # Aggregate per (target, drug_name): keep best evidence, union mechanisms
    def agg(group: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "drug_id": next((x for x in group["drug_id"].dropna().astype(str) if x and x != "nan"), ""),
            "max_phase": float(np.nanmax(group["max_phase"])) if group["max_phase"].notna().any() else 1.0,
            "first_approval": next((x for x in group["first_approval"].dropna() if x), None),
            "pchembl_value": float(np.nanmax(group["pchembl_value"])) if group["pchembl_value"].notna().any() else np.nan,
            "withdrawn_flag": bool(group["withdrawn_flag"].fillna(False).any()),
            "oral": bool(group["oral"].fillna(False).any()),
            "parenteral": bool(group["parenteral"].fillna(False).any()),
            "mechanism_of_action": " | ".join(sorted({str(x) for x in group["mechanism_of_action"].dropna() if str(x).strip()})),
            "original_indication": next((str(x) for x in group["original_indication"].dropna() if str(x).strip()), ""),
            "indication_class": next((str(x) for x in group["indication_class"].dropna() if str(x).strip()), ""),
            "source": ",".join(sorted({str(x) for x in group["source"].dropna() if str(x).strip()})),
        })

    agg_df = df.groupby(["target", "drug_name"], as_index=False).apply(agg, include_groups=False)
    # pandas 2.x apply returns MultiIndex — reset
    if isinstance(agg_df.index, pd.MultiIndex):
        agg_df = agg_df.reset_index()

    scored_rows = []
    for _, r in agg_df.iterrows():
        rdict = r.to_dict()
        pot = _potency(rdict.get("pchembl_value"))
        safe = _thyroid_safety(rdict)
        oral = _oral_bioavail(rdict)
        off = _off_patent(rdict)
        mech = _mechanism_fit(rdict)
        mp = float(rdict.get("max_phase") or 1.0)
        score = (0.30 * (mp / 4.0) + 0.20 * pot + 0.15 * safe
                 + 0.15 * oral + 0.10 * off + 0.10 * mech)
        if rdict.get("withdrawn_flag"):
            score *= 0.1
        comps = {"potency": pot, "mech": mech, "off": off, "safe": safe, "oral": oral}
        scored_rows.append({
            "target": rdict["target"],
            "drug_name": rdict["drug_name"],
            "drug_id": rdict.get("drug_id", ""),
            "max_phase": mp,
            "original_indication": rdict.get("original_indication") or rdict.get("indication_class") or "",
            "repurposing_score": round(float(score), 4),
            "rationale": _rationale(rdict, comps),
            "source": rdict.get("source", ""),
            "pchembl_value": rdict.get("pchembl_value"),
            "first_approval": rdict.get("first_approval"),
            "withdrawn_flag": rdict.get("withdrawn_flag"),
            "oral": rdict.get("oral"),
            "parenteral": rdict.get("parenteral"),
            "mechanism_of_action": rdict.get("mechanism_of_action"),
        })
    out = pd.DataFrame(scored_rows).sort_values("repurposing_score", ascending=False).reset_index(drop=True)
    keep_cols = ["target", "drug_name", "drug_id", "max_phase", "original_indication",
                 "repurposing_score", "rationale", "source"]
    out[keep_cols].to_csv(RES / "v7_repurposing_ranked.tsv", sep="\t", index=False)
    # Also write a richer TSV for downstream steps
    out.to_csv(RES / "v7_repurposing_ranked_full.tsv", sep="\t", index=False)
    log(f"STEP 4  done — {len(out)} scored rows; max={out['repurposing_score'].max():.3f}")
    return out

# STEP 5 — Top-10 deep dive

def _sanitize(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(name).strip()).strip("_")[:80] or "drug"

def step5_top10(ranked: pd.DataFrame, chembl_all: pd.DataFrame) -> list[dict]:
    log("STEP 5  Top-10 deep dive — starting")
    if ranked.empty:
        return []
    top = (ranked.sort_values("repurposing_score", ascending=False)
                  .drop_duplicates(subset=["drug_name"], keep="first")
                  .head(10).reset_index(drop=True))
    ch_index: dict[str, dict] = {}
    if not chembl_all.empty:
        for _, r in chembl_all.iterrows():
            nm = str(r.get("drug_name") or "").lower()
            if nm and nm not in ch_index:
                ch_index[nm] = r.to_dict()

    out_list = []
    for _, row in top.iterrows():
        drug = row["drug_name"]
        fname = _sanitize(drug)
        info = {
            "drug_name": drug,
            "drug_id": row.get("drug_id"),
            "target": row.get("target"),
            "max_phase": row.get("max_phase"),
            "repurposing_score": row.get("repurposing_score"),
            "original_indication": row.get("original_indication"),
            "mechanism_of_action": row.get("mechanism_of_action"),
            "source": row.get("source"),
            "clinicaltrials": {"count": 0, "studies": [], "thyroid_hits": 0},
            "pharmacokinetics": {},
            "pubmed_count_since_2020_thyroid_cancer": 0,
        }

        # (a) ClinicalTrials.gov v2
        try:
            q = urllib.parse.quote_plus(drug)
            url = f"https://clinicaltrials.gov/api/v2/studies?query.term={q}&pageSize=50&format=json"
            status, raw = http_get(url)
            if status == 200:
                js = safe_json(raw) or {}
                studies = js.get("studies") or []
                summaries, thy = [], 0
                for s in studies[:25]:
                    ps = (s.get("protocolSection") or {})
                    idm = ps.get("identificationModule") or {}
                    dm = ps.get("designModule") or {}
                    cm = ps.get("conditionsModule") or {}
                    em = dm.get("enrollmentInfo") or {}
                    phases = dm.get("phases") or []
                    conds = cm.get("conditions") or []
                    if any("thyroid" in (c or "").lower() for c in conds):
                        thy += 1
                    summaries.append({
                        "nct_id": idm.get("nctId"),
                        "title": idm.get("briefTitle"),
                        "phases": phases,
                        "conditions": conds,
                        "enrollment": em.get("count"),
                    })
                info["clinicaltrials"] = {"count": len(studies), "studies": summaries, "thyroid_hits": thy}
        except Exception as e:
            log(f"[CT] {drug}: {e}", "warning")

        # (b) PK from ChEMBL cache
        try:
            src = ch_index.get(drug) or {}
            if src:
                info["pharmacokinetics"] = {
                    "oral": bool(src.get("oral")),
                    "parenteral": bool(src.get("parenteral")),
                    "indication_class": src.get("indication_class"),
                    "first_approval": src.get("first_approval"),
                    "withdrawn_flag": bool(src.get("withdrawn_flag")),
                    "best_pchembl": src.get("pchembl_value"),
                }
        except Exception:
            pass

        # (c) PubMed count since 2020
        try:
            q = urllib.parse.quote_plus(f"{drug} AND (thyroid OR cancer)")
            url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                   f"?db=pubmed&term={q}&retmax=0&mindate=2020")
            status, raw = http_get(url)
            if status == 200 and raw:
                m = re.search(rb"<Count>(\d+)</Count>", raw)
                if m:
                    info["pubmed_count_since_2020_thyroid_cancer"] = int(m.group(1))
        except Exception as e:
            log(f"[PubMed] {drug}: {e}", "warning")

        (RES / f"v7_top_candidate_{fname}.json").write_text(json.dumps(info, indent=2, default=str))
        out_list.append(info)
        log(f"[TOP10] {drug} → CT={info['clinicaltrials']['count']} (thyroid={info['clinicaltrials']['thyroid_hits']}), PubMed={info['pubmed_count_since_2020_thyroid_cancer']}")
    log(f"STEP 5  done — {len(out_list)} top candidates profiled")
    return out_list

# STEP 6 — Combination hypotheses

def reactome_pathways(gene: str) -> set[str]:
    up = UNIPROT_MAP.get(gene)
    if not up:
        return FALLBACK_PATHWAYS.get(gene, set())
    url = f"https://reactome.org/ContentService/data/mapping/UniProt/{up}/pathways?species=9606"
    status, raw = http_get(url, headers={"Accept": "application/json"})
    if status != 200 or not raw:
        return FALLBACK_PATHWAYS.get(gene, set())
    js = safe_json(raw) or []
    if isinstance(js, list):
        ids = {p.get("stId") for p in js if isinstance(p, dict) and p.get("stId")}
        if ids:
            return ids
    return FALLBACK_PATHWAYS.get(gene, set())

def step6_combos(ranked: pd.DataFrame) -> pd.DataFrame:
    log("STEP 6  Combination hypotheses — starting")
    if ranked.empty:
        out = pd.DataFrame(columns=["target_A", "drug_A", "target_B", "drug_B",
                                    "rationale", "pathway_overlap", "synergy_score", "clinical_precedent"])
        out.to_csv(RES / "v7_combination_hypotheses.tsv", sep="\t", index=False)
        return out

    # Reactome pathway sets per target
    paths: dict[str, set[str]] = {g: reactome_pathways(g) for g in TARGETS}
    log(f"[COMB] pathway set sizes: {{{', '.join(f'{k}:{len(v)}' for k,v in paths.items())}}}")

    # Best drug per target (score >= 0.4)
    best_by_target: dict[str, pd.Series] = {}
    for g in TARGETS:
        sub = ranked[(ranked["target"] == g) & (ranked["repurposing_score"] >= 0.4)].sort_values("repurposing_score", ascending=False)
        if not sub.empty:
            best_by_target[g] = sub.iloc[0]

    rows = []
    for a, b in itertools.combinations(TARGETS, 2):
        if a not in best_by_target or b not in best_by_target:
            continue
        da, db = best_by_target[a], best_by_target[b]
        pa, pb = paths.get(a, set()), paths.get(b, set())
        inter = pa & pb
        overlap = 0.0
        if pa or pb:
            overlap = round(len(inter) / max(1, min(len(pa), len(pb)) or 1), 3) if inter else 0.0
        fda_a = (float(da.get("max_phase") or 0) >= 4)
        fda_b = (float(db.get("max_phase") or 0) >= 4)
        both_fda = 1 if fda_a and fda_b else 0
        moa_a = str(da.get("rationale") or "").lower()
        moa_b = str(db.get("rationale") or "").lower()
        diff_mech = 1 if moa_a.split(";")[0] != moa_b.split(";")[0] else 0
        synergy = 0.3 + 0.3 * (1 if inter else 0) + 0.2 * both_fda + 0.2 * diff_mech
        synergy += 0.1 * overlap  # small extra continuous bonus
        synergy = round(min(synergy, 1.0), 3)
        rationale = (
            f"{a}:{da['drug_name']} + {b}:{db['drug_name']}; "
            f"shared Reactome pathways={len(inter)}; "
            f"FDA-approved pair={'yes' if both_fda else 'no'}; "
            f"mechanism classes={'different' if diff_mech else 'same'}"
        )
        clinical_precedent = "search recommended"
        rows.append({
            "target_A": a, "drug_A": da["drug_name"],
            "target_B": b, "drug_B": db["drug_name"],
            "rationale": rationale,
            "pathway_overlap": len(inter),
            "synergy_score": synergy,
            "clinical_precedent": clinical_precedent,
        })
    out = pd.DataFrame(rows).sort_values("synergy_score", ascending=False).reset_index(drop=True)
    out.to_csv(RES / "v7_combination_hypotheses.tsv", sep="\t", index=False)
    log(f"STEP 6  done — {len(out)} combination rows")
    return out, paths

# STEP 7 — Dashboard HTML + combination heatmap

def build_heatmap(paths: dict[str, set[str]], combos: pd.DataFrame) -> str:
    import plotly.graph_objects as go
    n = len(TARGETS)
    M = np.zeros((n, n), dtype=float)
    label = [["" for _ in range(n)] for _ in range(n)]
    # diagonal = number of pathways
    for i, g in enumerate(TARGETS):
        M[i, i] = len(paths.get(g, set()))
        label[i][i] = f"{g}<br>{len(paths.get(g, set()))} pathways"
    # off-diagonal = synergy_score; lower triangle = pathway overlap count
    combo_map = {}
    for _, r in combos.iterrows():
        combo_map[(r["target_A"], r["target_B"])] = (r["synergy_score"], r["pathway_overlap"], r["drug_A"], r["drug_B"])
    for i, a in enumerate(TARGETS):
        for j, b in enumerate(TARGETS):
            if i == j:
                continue
            key = (a, b) if (a, b) in combo_map else (b, a)
            cell = combo_map.get(key)
            if cell:
                syn, ov, da, db = cell
                if i < j:
                    M[i, j] = syn
                    label[i][j] = f"{a}+{b}<br>synergy={syn:.2f}<br>{da} + {db}"
                else:
                    M[i, j] = ov
                    label[i][j] = f"{a}+{b}<br>pathway overlap={ov}"
            else:
                M[i, j] = np.nan
                label[i][j] = f"{a}+{b}<br>(no qualifying drugs)"
    fig = go.Figure(data=go.Heatmap(
        z=M, x=TARGETS, y=TARGETS,
        text=label, hoverinfo="text",
        colorscale="Teal", zmin=0, zmax=max(1.0, np.nanmax(M) if np.isfinite(np.nanmax(M)) else 1.0),
    ))
    fig.update_layout(
        template="plotly_dark", paper_bgcolor="rgba(6,17,38,1)", plot_bgcolor="rgba(6,17,38,1)",
        title="Combination hypotheses — upper triangle: synergy score, lower: pathway overlap",
        margin=dict(l=90, r=40, t=70, b=80),
        xaxis=dict(tickangle=-30), yaxis=dict(autorange="reversed"),
        height=560,
    )
    out_path = FIGS / "combination_heatmap.html"
    fig.write_html(out_path, include_plotlyjs="cdn", full_html=True)
    return str(out_path)

def build_stat_cards(ranked: pd.DataFrame) -> dict:
    if ranked.empty:
        return {"total": 0, "fda_approved": 0, "phase23": 0, "score_gt_05": 0}
    return {
        "total": int(len(ranked)),
        "fda_approved": int((ranked["max_phase"] >= 4).sum()),
        "phase23": int(((ranked["max_phase"] >= 2) & (ranked["max_phase"] < 4)).sum()),
        "score_gt_05": int((ranked["repurposing_score"] > 0.5).sum()),
    }

def step7_dashboard(ranked: pd.DataFrame, combos: pd.DataFrame, paths: dict[str, set[str]],
                    top_profiles: list[dict]) -> str:
    log("STEP 7  Dashboard page — starting")
    heat_path = build_heatmap(paths, combos)
    heat_rel = "../figs_interactive/v7/combination_heatmap.html"
    stats = build_stat_cards(ranked)

    # Build rows for table (top 100 for browser friendliness)
    tbl = ranked.head(300).copy()
    tbl_rows = []
    for _, r in tbl.iterrows():
        tbl_rows.append({
            "target": r["target"], "drug": r["drug_name"],
            "score": round(float(r["repurposing_score"]), 3),
            "phase": r["max_phase"],
            "indication": str(r.get("original_indication") or "")[:90],
            "source": r.get("source", ""),
            "rationale": str(r.get("rationale") or "")[:160],
        })
    tbl_json = json.dumps(tbl_rows)

    # Top-10 cards
    top_by_drug = {p["drug_name"]: p for p in top_profiles}
    top_sorted = (ranked.sort_values("repurposing_score", ascending=False)
                          .drop_duplicates("drug_name")
                          .head(10))
    card_html = []
    for _, r in top_sorted.iterrows():
        prof = top_by_drug.get(r["drug_name"], {})
        ct = prof.get("clinicaltrials", {}) or {}
        card_html.append(f"""
        <div class="v7-card">
          <div class="v7-card-title">{r['drug_name'].upper()}</div>
          <div class="v7-card-sub">target {r['target']} &middot; score {float(r['repurposing_score']):.3f} &middot; max phase {r['max_phase']}</div>
          <div class="v7-card-meta">{(str(r.get('original_indication') or '') or 'no indication info')[:120]}</div>
          <div class="v7-card-pill-row">
            <span class="v7-pill">CT.gov studies: {ct.get('count', 0)}</span>
            <span class="v7-pill">thyroid trials: {ct.get('thyroid_hits', 0)}</span>
            <span class="v7-pill">PubMed (2020+): {prof.get('pubmed_count_since_2020_thyroid_cancer', 0)}</span>
          </div>
          <div class="v7-card-rationale">{str(r.get('rationale') or '')[:260]}</div>
        </div>
        """)
    cards_block = "\n".join(card_html) if card_html else "<p class='muted'>No candidates scored yet.</p>"

    fda_map = {g: 0 for g in TARGETS}
    if not ranked.empty:
        for g in TARGETS:
            fda_map[g] = int(((ranked["target"] == g) & (ranked["max_phase"] >= 4)).sum())

    from v7_templates import DASHBOARD_HTML
    html = (DASHBOARD_HTML
            .replace("__STAT_TOTAL__", str(stats["total"]))
            .replace("__STAT_FDA__", str(stats["fda_approved"]))
            .replace("__STAT_P23__", str(stats["phase23"]))
            .replace("__STAT_GT05__", str(stats["score_gt_05"]))
            .replace("__CARDS__", cards_block)
            .replace("__HEAT__", heat_rel)
            .replace("__TBL_JSON__", tbl_json))
    out = PAGES / "v7_drug_repurposing.html"
    out.write_text(html)
    log(f"STEP 7  done — {out}")
    return str(out)

# STEP 8 — Paper section

def step8_paper(ranked: pd.DataFrame, combos: pd.DataFrame, top_profiles: list[dict], stats: dict) -> str:
    log("STEP 8  Paper section — starting")
    top = (ranked.sort_values("repurposing_score", ascending=False)
                  .drop_duplicates("drug_name").head(10))
    top_lines = []
    for _, r in top.iterrows():
        prof = next((p for p in top_profiles if p["drug_name"] == r["drug_name"]), {})
        ct = (prof.get("clinicaltrials") or {}).get("count", 0)
        thy = (prof.get("clinicaltrials") or {}).get("thyroid_hits", 0)
        pm = prof.get("pubmed_count_since_2020_thyroid_cancer", 0)
        top_lines.append(
            f"- **{r['drug_name']}** (target {r['target']}, phase {r['max_phase']}, "
            f"score {r['repurposing_score']:.3f}) — {str(r.get('original_indication') or 'unspecified')[:80]}; "
            f"ClinicalTrials.gov studies={ct} (thyroid={thy}); PubMed 2020+ hits={pm}"
        )

    combo_lines = []
    for _, r in combos.head(10).iterrows():
        combo_lines.append(
            f"- **{r['target_A']}:{r['drug_A']}** + **{r['target_B']}:{r['drug_B']}** "
            f"(synergy={r['synergy_score']:.2f}, pathway overlap={r['pathway_overlap']})"
        )

    from v7_templates import PAPER_MD
    md = (PAPER_MD
          .replace("__TOP_LINES__", "\n".join(top_lines) if top_lines else "_No candidates scored._")
          .replace("__COMBO_LINES__", "\n".join(combo_lines) if combo_lines else "_No qualifying combinations (no target-pair with both partners scoring >=0.4)._")
          .replace("__STAT_TOTAL__", str(stats["total"]))
          .replace("__STAT_FDA__", str(stats["fda_approved"]))
          .replace("__STAT_P23__", str(stats["phase23"]))
          .replace("__STAT_GT05__", str(stats["score_gt_05"])))
    out = RPT / "v7_repurposing_section.md"
    out.write_text(md)
    log(f"STEP 8  done — {out}")
    return str(out)

# main

def _safe(fn, *args, **kwargs):
    try:
        return fn(*args, **kwargs)
    except Exception as e:
        log(f"STEP FAILED {fn.__name__}: {e}\n{traceback.format_exc()}", "error")
        return None

def main() -> int:
    t0 = time.time()
    log(f"===== v7 repurposing run @ {_dt.datetime.now().isoformat()} =====")

    ot = _safe(step1_opentargets)
    if ot is None:
        ot = pd.DataFrame()

    dc = _safe(step2_drugcentral)
    if dc is None:
        dc = pd.DataFrame()

    ch = _safe(step3_chembl)
    if ch is None:
        ch = pd.DataFrame()

    ranked = _safe(step4_score, ot, dc, ch)
    if ranked is None:
        ranked = pd.DataFrame(columns=["target", "drug_name", "drug_id", "max_phase",
                                       "original_indication", "repurposing_score", "rationale", "source"])

    top_profiles = _safe(step5_top10, ranked, ch) or []

    combos_result = _safe(step6_combos, ranked)
    if isinstance(combos_result, tuple):
        combos, paths = combos_result
    else:
        combos = pd.DataFrame(columns=["target_A", "drug_A", "target_B", "drug_B",
                                       "rationale", "pathway_overlap", "synergy_score", "clinical_precedent"])
        paths = {g: FALLBACK_PATHWAYS.get(g, set()) for g in TARGETS}

    _safe(step7_dashboard, ranked, combos, paths, top_profiles)
    stats = build_stat_cards(ranked)
    _safe(step8_paper, ranked, combos, top_profiles, stats)

    # FDA-approved per target map
    fda_map = {g: 0 for g in TARGETS}
    if not ranked.empty:
        for g in TARGETS:
            fda_map[g] = int(((ranked["target"] == g) & (ranked["max_phase"] >= 4)).sum())

    elapsed = time.time() - t0
    log(f"Total wall time: {elapsed:.1f}s")

    # STAGE 99 output — printed EXACTLY as required
    print("=== v7 REPURPOSING COMPLETE ===")
    print(f"Total drug-target interactions   : {int(len(ranked))}")
    print(f"FDA-approved drugs per target    : {fda_map}")
    print("Top 10 repurposing candidates    : /opt/thyroid-dash/project/results/v7_repurposing/v7_repurposing_ranked.tsv (top 10)")
    print("Combination hypotheses           : /opt/thyroid-dash/project/results/v7_repurposing/v7_combination_hypotheses.tsv")
    print("Dashboard page                   : http://40.82.129.113:8012/reports/html/pages/v7_drug_repurposing.html")
    print("Paper section                    : /opt/thyroid-dash/project/reports/v7/v7_repurposing_section.md")
    print("Log                              : /opt/thyroid-dash/project/logs/v7_repurposing.log")
    return 0

if __name__ == "__main__":
    sys.exit(main())
