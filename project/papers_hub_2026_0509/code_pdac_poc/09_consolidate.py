"""
Consolidate everything into one JSON for the HTML to consume.
"""

import json
from pathlib import Path
import pandas as pd

ROOT = Path("/data/pdac_poc")
PROC = ROOT / "processed"
RES = ROOT / "results"

out = {}

# Manifest
out["data_manifest"] = json.load(open(PROC / "manifest_01_fetch.json"))

# Registry (subset useful for HTML)
reg = json.load(open(ROOT / "raw/registry/pdac_registry.json"))
out["registry_summary"] = {
    "n_entries": reg["n_entries"],
    "n_open": reg["n_open_data_cohorts"],
    "n_controlled": reg["n_controlled"],
    "total_n": reg["total_n"],
}
out["registry_top"] = sorted(
    [r for r in reg["registry"] if isinstance(r["n"], int)],
    key=lambda r: r["n"], reverse=True
)[:20]
out["registry_full"] = reg["registry"]

# Moffitt
out["moffitt"] = json.load(open(PROC / "moffitt_summary.json"))

# Neoantigen
out["neoantigen_quality"] = json.load(open(PROC / "neoantigen_quality_summary.json"))

# Immune readiness
out["immune_readiness"] = json.load(open(PROC / "pdac_immune_readiness_summary.json"))

# Vaccine priority
out["vaccine_priority"] = json.load(open(PROC / "vaccine_priority_summary.json"))
vt = pd.read_csv(RES / "vaccine_targets.tsv", sep="\t")
out["vaccine_targets_table"] = vt.to_dict(orient="records")

# Agent run
out["agent_latest_run"] = json.load(open(RES / "agent_run_latest.json"))

# Per-figure manifest
out["figures"] = sorted(p.name for p in (ROOT / "figures").glob("*.png"))

with open(RES / "pdac_poc_consolidated.json", "w") as f:
    json.dump(out, f, indent=2, default=str)

print("[09] consolidated →", RES / "pdac_poc_consolidated.json")
print(f"  manifest n_samples={out['data_manifest']['samples']}")
print(f"  moffitt basal={out['moffitt']['n_basal']} class={out['moffitt']['n_classical']}")
print(f"  neoQ samples={out['neoantigen_quality']['n_samples_scored']}  KRAS-G12={out['neoantigen_quality']['kras_g12_count']}")
print(f"  vaccine targets={out['vaccine_priority']['n_targets_ranked']}")
print(f"  Korean coverage={out['vaccine_priority']['korean_off_the_shelf_population_coverage_pct']}%")
print(f"  registry n={out['registry_summary']['n_entries']} open={out['registry_summary']['n_open']}")
