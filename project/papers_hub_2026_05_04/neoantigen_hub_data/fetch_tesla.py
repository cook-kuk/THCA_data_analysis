"""TESLA · Wells 2020 Cell supplementary tables direct download from PMC.

PMID 33038342 · PMC7652061 · 7 xlsx supplementary tables.
"""
import json, time, requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

OUT = Path("/data/neoantigen_vaccine_hub/data_raw/tesla/cell_supp_tables")
OUT.mkdir(parents=True, exist_ok=True)
LOG = Path("/data/neoantigen_vaccine_hub/logs"); LOG.mkdir(parents=True, exist_ok=True)

PMC_BASE = "https://pmc.ncbi.nlm.nih.gov/articles/instance/7652061/bin"
FILES = [
    f"NIHMS1634895-supplement-{i}.xlsx" for i in range(1,8)
] + ["NIHMS1634895-supplement-8.pdf"]

S = requests.Session()
S.headers.update({"User-Agent":"Mozilla/5.0 Lumenix-NV-Hub/1.0",
                   "Accept":"*/*"})


def download_one(fname):
    url = f"{PMC_BASE}/{fname}"
    dest = OUT/fname
    if dest.exists() and dest.stat().st_size > 1000:
        return {"file":fname,"status":"ALREADY_EXISTS","size_kb":dest.stat().st_size//1024}
    t0 = time.time()
    try:
        r = S.get(url, timeout=120, stream=True)
        if not r.ok:
            return {"file":fname,"status":"FAILED","http":r.status_code,"url":url}
        with open(dest,"wb") as f:
            for chunk in r.iter_content(8192):
                f.write(chunk)
        size = dest.stat().st_size
        if size < 500:
            return {"file":fname,"status":"FAILED","reason":"size too small","size":size}
        return {"file":fname,"status":"SUCCESS","size_kb":size//1024,"elapsed_s":round(time.time()-t0,1)}
    except Exception as e:
        return {"file":fname,"status":"FAILED","err":str(e)[:200]}


def main():
    print(f"[TESLA] downloading {len(FILES)} supplementary files from PMC7652061 …")
    results = []
    with ThreadPoolExecutor(max_workers=4) as ex:
        for r in ex.map(download_one, FILES):
            results.append(r)
            sz = r.get('size_kb','?')
            print(f"   {r['status']:14s}  {r['file']:38s}  {sz} KB  ({r.get('elapsed_s','?')}s)")
    n_ok = sum(1 for r in results if r["status"] in ("SUCCESS","ALREADY_EXISTS"))
    summary = {"source":"TESLA","n_files":len(FILES),"n_success":n_ok,
               "files":results,"fetched_at":time.strftime("%Y-%m-%d %H:%M:%S"),
               "url_base":PMC_BASE,"pmid":"33038342","doi":"10.1016/j.cell.2020.09.015"}
    json.dump(summary, open(LOG/"tesla_status.json","w"), indent=2, default=str)
    print(f"[TESLA] {n_ok}/{len(FILES)} fetched · status saved")


if __name__ == "__main__":
    main()
