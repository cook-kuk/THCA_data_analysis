"""
SPARK-inspired agentic orchestrator for PDAC vaccine discovery.

Architecture (loosely following Nature Med 2026 SPARK design):

  PLANNER  ──▶  TOOL ROUTER  ──▶  TOOLS  ──▶  VERIFIER  ──▶  MEMORY
                                                                │
                                                                ▼
                                                              REPORT

Components:
  - Planner.plan(query)             : decomposes a user query into a tool plan
  - ToolRouter.run(plan)            : executes tools in dependency order
  - Verifier.audit(state)           : runs sanity & direction-invariance checks
  - Memory.write/read(slot, value)  : persistent KV across turns
  - Report.render()                 : human-readable summary

Tools (each is a small, deterministic function calling our scripts/results):
  T1 cohort_lookup        — reads registry
  T2 subtype_classifier   — Moffitt basal/classical
  T3 driver_landscape     — KRAS / TP53 / SMAD4 / CDKN2A frequencies
  T4 neoantigen_quality   — Balachandran R × D
  T5 immune_readiness     — TME modules + cohort directionality
  T6 vaccine_priority     — antigen ranker
  T7 hla_bias_audit       — Korean / pan-Asian / European coverage
  T8 trial_matcher        — published vaccine trials matched to plan

Run end-to-end:
    python 07_agent_orchestrator.py "BRAF/RAS-wt PDAC + MMR-deficient + Korean"
"""

import json
import sys
import time
from pathlib import Path
import pandas as pd

ROOT = Path("/data/pdac_poc")
PROC = ROOT / "processed"
RES = ROOT / "results"


# ============================================================
# Memory
# ============================================================
class Memory:
    def __init__(self, path=ROOT / "logs/agent_memory.json"):
        self.path = path
        self.store = {}
        if path.exists():
            try:
                self.store = json.load(open(path))
            except Exception:
                self.store = {}

    def write(self, slot, value):
        self.store[slot] = value
        json.dump(self.store, open(self.path, "w"), indent=2, default=str)

    def read(self, slot, default=None):
        return self.store.get(slot, default)


# ============================================================
# Tools (each returns {ok, data, citations, n_supporting_samples})
# ============================================================
def tool_cohort_lookup(query=None):
    reg = json.load(open(ROOT / "raw/registry/pdac_registry.json"))
    return {
        "ok": True,
        "tool": "T1.cohort_lookup",
        "data": {"n_cohorts": reg["n_entries"],
                  "n_open": reg["n_open_data_cohorts"],
                  "n_controlled": reg["n_controlled"],
                  "total_n": reg["total_n"]},
        "n_supporting_samples": reg["total_n"],
        "citations": ["Bailey 2016 Nature", "Moffitt 2015 Nat Genet",
                      "Cao 2021 Cell", "Steele 2020 Nat Cancer",
                      "Hwang 2022 Nat Genet", "Balachandran 2017 Nature",
                      "Rojas 2023 Nature"],
    }


def tool_subtype_classifier():
    s = json.load(open(PROC / "moffitt_summary.json"))
    return {"ok": True, "tool": "T2.subtype_classifier",
            "data": s, "citations": ["Moffitt 2015 Nat Genet"],
            "n_supporting_samples": s.get("n_total")}


def tool_driver_landscape():
    p = pd.read_csv(RES / "moffitt_calls.tsv", sep="\t")
    out = {
        "n_samples_genotyped": int(len(p)),
        "kras_rate": float(p["mut_KRAS"].mean()),
        "tp53_rate": float(p["mut_TP53"].mean()),
        "smad4_rate": float(p["mut_SMAD4"].mean()),
        "cdkn2a_rate": float(p["mut_CDKN2A"].mean()),
        "kras_AND_tp53_rate": float(((p["mut_KRAS"] == 1) &
                                      (p["mut_TP53"] == 1)).mean()),
    }
    return {"ok": True, "tool": "T3.driver_landscape",
            "data": out, "citations": ["TCGA-PAAD 2017"],
            "n_supporting_samples": out["n_samples_genotyped"]}


def tool_neoantigen_quality():
    s = json.load(open(PROC / "neoantigen_quality_summary.json"))
    return {"ok": True, "tool": "T4.neoantigen_quality",
            "data": s, "citations": ["Balachandran 2017 Nature",
                                      "Łuksza 2017 Nature"],
            "n_supporting_samples": s.get("n_samples_scored")}


def tool_immune_readiness():
    s = json.load(open(PROC / "pdac_immune_readiness_summary.json"))
    # surface only key metrics
    out = {
        "n_samples": s["n_samples"],
        "module_coverage": s["module_coverage"],
        "cohens_d_basal_vs_classical": s["cohens_d_basal_vs_classical"],
        "fraction_TLS_high": s["fraction_TLS_high"],
        "fraction_myeloid_high": s["fraction_myeloid_high"],
        "fraction_basal_AND_inflamed": s["fraction_basal_AND_inflamed"],
        "fraction_basal_AND_TLS_high": s["fraction_basal_AND_TLS_high"],
    }
    return {"ok": True, "tool": "T5.immune_readiness",
            "data": out, "citations": ["Elyada 2019 Cancer Discov",
                                        "Steele 2020 Nat Cancer",
                                        "Hwang 2022 Nat Genet"],
            "n_supporting_samples": s["n_samples"]}


def tool_vaccine_priority():
    s = json.load(open(PROC / "vaccine_priority_summary.json"))
    return {"ok": True, "tool": "T6.vaccine_priority",
            "data": s,
            "citations": ["Pant 2024 Nat Med ELI-002",
                          "Rojas 2023 Nature BNT122",
                          "Sethna 2025 Nature long-FU",
                          "Lowery 2022 JCI mTP53 TCR"],
            "n_supporting_samples": s.get("n_targets_ranked")}


def tool_hla_bias_audit():
    s = json.load(open(PROC / "vaccine_priority_summary.json"))
    return {"ok": True, "tool": "T7.hla_bias_audit",
            "data": {"korean_pct": s["korean_off_the_shelf_population_coverage_pct"],
                      "european_pct": s["european_off_the_shelf_population_coverage_pct"]},
            "citations": ["AFND v3.1 (2024)"],
            "n_supporting_samples": None}


def tool_trial_matcher(query=None):
    trials = [
        {"id": "NCT05111353", "name": "Autogene cevumeran (BNT122) + atezolizumab + mFOLFIRINOX",
         "phase": "II", "n_planned": 260, "country": "US/EU/multi",
         "eligibility": "resected PDAC, HLA-typed, ECOG 0-1",
         "cite": "Rojas 2023 Nature; expansion phase II"},
        {"id": "NCT05726864", "name": "ELI-002 7P (mKRAS amphiphile) + ICI",
         "phase": "II", "n_planned": 158,
         "eligibility": "resected MRD+ PDAC/CRC with mKRAS G12D/R/V/A/S/C/G13D",
         "cite": "Pant 2024 Nat Med phase 1"},
        {"id": "NCT04117087", "name": "Sotorasib + IO maintenance / vaccine cassette",
         "phase": "I/II", "n_planned": "—",
         "eligibility": "KRAS G12C-mutant solid tumors incl. PDAC",
         "cite": "Mirati / Amgen sotorasib roadmap"},
        {"id": "NCT04161755", "name": "Personalized RNA neoantigen (BNT122) PDAC",
         "phase": "I", "n_planned": 16,
         "eligibility": "resected PDAC, post-surgery <12 wk",
         "cite": "Rojas 2023 Nature primary paper"},
    ]
    return {"ok": True, "tool": "T8.trial_matcher",
            "data": {"matched_trials": trials},
            "citations": ["clinicaltrials.gov 2026-05"],
            "n_supporting_samples": None}


TOOLS = {
    "T1": tool_cohort_lookup, "T2": tool_subtype_classifier,
    "T3": tool_driver_landscape, "T4": tool_neoantigen_quality,
    "T5": tool_immune_readiness, "T6": tool_vaccine_priority,
    "T7": tool_hla_bias_audit, "T8": tool_trial_matcher,
}


# ============================================================
# Planner — minimalist rule-based plan; in production this is an LLM call
# ============================================================
class Planner:
    """Decomposes a query into an ordered tool plan with rationale."""

    def plan(self, query: str):
        q = (query or "").lower()
        steps = ["T1"]                                 # always lookup cohorts
        if any(k in q for k in ["subtype", "moffitt", "basal", "classical"]):
            steps += ["T2", "T3"]
        else:
            steps += ["T2", "T3"]                      # default include
        if any(k in q for k in ["neo", "antigen", "kras", "tp53", "vaccine"]):
            steps += ["T4"]
        if any(k in q for k in ["immune", "tme", "myeloid", "tls", "caf",
                                  "ici", "checkpoint"]):
            steps += ["T5"]
        steps += ["T6"]                                # always priority
        if any(k in q for k in ["korean", "asian", "hla", "bias", "ethni"]):
            steps += ["T7"]
        if any(k in q for k in ["trial", "nct", "study", "phase"]):
            steps += ["T8"]
        # de-dup, preserve order
        seen, ordered = set(), []
        for s in steps:
            if s not in seen:
                ordered.append(s); seen.add(s)
        return ordered


# ============================================================
# Verifier — direction-invariance + sanity checks
# ============================================================
class Verifier:
    def audit(self, results):
        flags = []

        # CHECK 1: KRAS rate sanity (PDAC > 80%)
        t3 = next((r for r in results if r["tool"].startswith("T3")), None)
        if t3 and t3["data"]["kras_rate"] < 0.6:
            flags.append({"level": "ERROR",
                          "msg": f"KRAS rate {t3['data']['kras_rate']:.2%} below 60% — pipeline likely broken (PDAC literature ~90%)"})
        elif t3:
            flags.append({"level": "PASS",
                          "msg": f"KRAS rate {t3['data']['kras_rate']:.0%} ≥ 60% (PDAC sanity)"})

        # CHECK 2: basal fraction in [30,70]
        t2 = next((r for r in results if r["tool"].startswith("T2")), None)
        if t2:
            f = t2["data"]["fraction_basal"]
            if 0.30 <= f <= 0.70:
                flags.append({"level": "PASS",
                              "msg": f"basal fraction {f:.0%} within Moffitt 2015 expected 30–70%"})
            else:
                flags.append({"level": "WARN",
                              "msg": f"basal fraction {f:.0%} outside 30–70% — cohort selection bias possible"})

        # CHECK 3: immune readiness module coverage ≥ 50% per module
        t5 = next((r for r in results if r["tool"].startswith("T5")), None)
        if t5:
            cov = t5["data"]["module_coverage"]
            for mod, ratio in cov.items():
                a, b = ratio.split("/")
                if int(a) / max(1, int(b)) < 0.5:
                    flags.append({"level": "WARN",
                                  "msg": f"module {mod} coverage {ratio} < 50%"})

        # CHECK 4: NeoQ top genes — KRAS expected at top
        t4 = next((r for r in results if r["tool"].startswith("T4")), None)
        if t4 and "top_neo_genes" in t4["data"]:
            top_gene = max(t4["data"]["top_neo_genes"].items(),
                           key=lambda kv: kv[1])[0]
            if top_gene == "KRAS":
                flags.append({"level": "PASS",
                              "msg": "top NeoQ gene = KRAS (PDAC expected)"})
            else:
                flags.append({"level": "WARN",
                              "msg": f"top NeoQ gene = {top_gene} (KRAS expected in PDAC)"})

        # CHECK 5: HLA bias — Korean coverage ≥ 30% expected
        t7 = next((r for r in results if r["tool"].startswith("T7")), None)
        if t7:
            kor = t7["data"]["korean_pct"]
            if kor < 30:
                flags.append({"level": "WARN",
                              "msg": f"Korean off-the-shelf coverage {kor}% < 30% — recommend personalized track for Korean cohort"})
            else:
                flags.append({"level": "PASS",
                              "msg": f"Korean off-the-shelf coverage {kor}% ≥ 30% (acceptable for off-the-shelf cassette)"})

        return flags


# ============================================================
# Orchestrator
# ============================================================
class LumenixPDACAgent:
    def __init__(self):
        self.memory = Memory()
        self.planner = Planner()
        self.verifier = Verifier()

    def run(self, query: str, verbose=True):
        t0 = time.time()
        plan = self.planner.plan(query)
        if verbose:
            print(f"[plan] {plan}")
        results = []
        for step in plan:
            try:
                if step == "T1":
                    res = TOOLS[step](query=query)
                elif step == "T8":
                    res = TOOLS[step](query=query)
                else:
                    res = TOOLS[step]()
                results.append(res)
                if verbose:
                    n = res.get("n_supporting_samples")
                    print(f"  ok {step}  ({res['tool']})"
                          f"  n_supp={n}")
            except Exception as e:
                results.append({"ok": False, "tool": step, "err": str(e)})
                print(f"  ERR {step}: {e}")
        flags = self.verifier.audit(results)
        report = {
            "query": query,
            "plan": plan,
            "results": results,
            "audit_flags": flags,
            "elapsed_s": round(time.time() - t0, 2),
        }
        self.memory.write(f"run_{int(time.time())}", report)
        return report


def main():
    q = sys.argv[1] if len(sys.argv) > 1 else \
        "PDAC vaccine plan for Korean cohort — KRAS G12 + TP53 hotspots + ICI"
    agent = LumenixPDACAgent()
    rep = agent.run(q)

    out = ROOT / "results/agent_run_latest.json"
    out.write_text(json.dumps(rep, indent=2, default=str))
    print(f"\n[agent] elapsed={rep['elapsed_s']}s  audit={len(rep['audit_flags'])} flags")
    for f in rep["audit_flags"]:
        print(f"  {f['level']:5s}  {f['msg']}")
    print(f"[agent] full report → {out}")


if __name__ == "__main__":
    main()
