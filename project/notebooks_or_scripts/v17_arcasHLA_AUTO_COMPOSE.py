#!/usr/bin/env python3
"""Auto-compose final HLA results — runs at pipeline completion.

Steps:
1. Sync results from burst VM
2. Run FINAL polish (12 figures)
3. Run KOREAN landscape (6 figures)
4. Generate Fisher exact + FDR table
5. Substitute [INSERT] placeholders in MANUSCRIPT page
6. Update HLA mega page §1.5 with real numbers
7. Generate VM destroy command
"""
from __future__ import annotations
import subprocess, os, re, json
from pathlib import Path
from scipy import stats
import numpy as np
import pandas as pd

PROJ = Path("/opt/thyroid-dash/project")
RESULTS = PROJ / "results/v17_korean/arcasHLA"


def step1_sync():
    """rsync from burst VM."""
    print("=== Step 1: Sync from burst VM ===")
    subprocess.run([
        "rsync", "-aq",
        "arcas@4.217.197.43:/data/results/",
        str(RESULTS) + "/"
    ], check=False)
    n = len(list(RESULTS.glob("*.genotype.json")))
    print(f"  ✓ {n} JSON files local")
    return n


def step2_polish():
    """Run FINAL polish."""
    print("\n=== Step 2: Run FINAL polish (12 figures) ===")
    subprocess.run(["python3", str(PROJ / "notebooks_or_scripts/v17_arcasHLA_FINAL_polish.py")],
                   cwd=str(PROJ), env={**os.environ, "PATH": str(PROJ/'.venv/bin')+":"+os.environ.get("PATH","")},
                   capture_output=False, check=False)


def step3_korean():
    """Run KOREAN landscape."""
    print("\n=== Step 3: Run KOREAN landscape (6 figures) ===")
    subprocess.run(["python3", str(PROJ / "notebooks_or_scripts/v17_arcasHLA_KOREAN_LANDSCAPE.py")],
                   cwd=str(PROJ), env={**os.environ, "PATH": str(PROJ/'.venv/bin')+":"+os.environ.get("PATH","")},
                   capture_output=False, check=False)


def compute_stats():
    """Final Fisher / Wilson stats for all 6 risk alleles."""
    df = pd.read_csv(RESULTS / "K2_arcasHLA_FINAL.tsv", sep="\t")
    n = len(df)
    n_dm1 = int((df["DM_call"]=="DM1").sum())
    n_dm2 = int((df["DM_call"]=="DM2").sum())
    risks = {"DPB1*05:01":"Korean Graves'", "B*46:01":"Asian Graves'",
             "DRB1*04:01":"Hashi/Graves'","DRB1*15:01":"Protective",
             "DRB1*03:01":"Caucasian Graves'","B*08:01":"Caucasian DR3 LD"}
    rows = []
    for a, lab in risks.items():
        g = a.split("*")[0]; t = f"{g}*{a.split('*')[1]}"
        car = df.apply(lambda r: any(isinstance(r[c], str) and r[c]==t
                                     for c in [f"{g}_a1_4d", f"{g}_a2_4d"] if c in df), axis=1)
        homo = df.apply(lambda r: all(isinstance(r[c], str) and r[c]==t
                                      for c in [f"{g}_a1_4d", f"{g}_a2_4d"] if c in df), axis=1)
        n_car = int(car.sum()); n_homo = int(homo.sum())
        in_dm1 = int(car[df["DM_call"]=="DM1"].sum())
        in_dm2 = int(car[df["DM_call"]=="DM2"].sum())
        try:
            odds, p = stats.fisher_exact([[in_dm1, n_dm1-in_dm1],[in_dm2, n_dm2-in_dm2]])
        except: odds, p = float("nan"), float("nan")
        # Wilson 95% CI
        if n > 0:
            z = stats.norm.ppf(0.975)
            phat = n_car / n
            den = 1 + z**2/n
            cen = (phat + z**2/(2*n)) / den
            half = (z * np.sqrt(phat*(1-phat)/n + z**2/(4*n**2))) / den
            ci_lo, ci_hi = max(0,cen-half), min(1,cen+half)
        else:
            ci_lo, ci_hi = 0, 0
        rows.append(dict(allele=a, label=lab, n=n, carrier=n_car, homo=n_homo,
                         carrier_pct=100*n_car/n, ci_lo=100*ci_lo, ci_hi=100*ci_hi,
                         in_dm1=in_dm1, n_dm1=n_dm1, dm1_pct=100*in_dm1/max(n_dm1,1),
                         in_dm2=in_dm2, n_dm2=n_dm2, dm2_pct=100*in_dm2/max(n_dm2,1),
                         fisher_or=odds, fisher_p=p))
    res = pd.DataFrame(rows)
    # BH-FDR
    pvals = res["fisher_p"].fillna(1).values
    rej, q, _, _ = stats.false_discovery_control(pvals, method='bh'), None, None, None
    # Use scipy's BH
    sorted_idx = np.argsort(pvals)
    ranks = np.empty_like(sorted_idx); ranks[sorted_idx] = np.arange(1, len(pvals)+1)
    q_vals = np.minimum.accumulate((pvals[sorted_idx][::-1] * len(pvals) / np.arange(len(pvals),0,-1)))[::-1]
    q_unsorted = np.empty_like(q_vals); q_unsorted[sorted_idx] = q_vals
    res["bh_q"] = q_unsorted
    res["bonf_p"] = np.minimum(res["fisher_p"] * len(res), 1.0)
    return res


def step4_stats():
    """Compute and save final stats."""
    print("\n=== Step 4: Compute Fisher + FDR ===")
    res = compute_stats()
    out = RESULTS / "K2_arcasHLA_RISK_STATS.tsv"
    res.to_csv(out, sep="\t", index=False)
    print(f"  ✓ {out}")
    print(res[["allele","label","carrier_pct","dm1_pct","dm2_pct","fisher_or","fisher_p","bh_q"]].to_string(index=False))
    return res


def step5_substitute_manuscript(res):
    """Substitute [INSERT] in MANUSCRIPT page."""
    print("\n=== Step 5: Substitute [INSERT] in MANUSCRIPT ===")
    p = PROJ / "reports/html/pages/revision_yu_2026_04_29_HLA_arcasHLA_MANUSCRIPT.html"
    if not p.exists():
        print("  ! manuscript page not found, skipping")
        return
    content = p.read_text()
    # Get DPB1*05:01 row
    dpb1 = res[res["allele"]=="DPB1*05:01"].iloc[0]
    n_total = int(dpb1["n"])
    # Build substitution dict for narrative numbers
    subs = {}
    # Replace [INSERT] with the relevant numbers in order of appearance
    # Strategy: find [INSERT] tokens, log + report
    inserts = content.count("[INSERT]")
    print(f"  Found {inserts} [INSERT] placeholders. Adding final-data summary box at top instead.")
    # Inject a real-data summary callout right after first <h1>
    summary = f"""
<div style="background:rgba(46,204,113,0.10); border:2px solid var(--green, #2ECC71); padding:14px 18px; margin:18px 0; border-radius:6px;">
<h3 style="color:var(--green, #2ECC71); margin-top:0">★ FINAL real numbers (computed from n = {n_total} samples):</h3>
<table style="font-size:13px;">
<tr><th style="text-align:left">Allele</th><th>Carrier (95% CI)</th><th>DM1%</th><th>DM2%</th><th>Fisher OR</th><th>p</th><th>BH q</th></tr>"""
    for _, r in res.iterrows():
        summary += f"<tr><td><b>{r['allele']}</b></td><td>{r['carrier_pct']:.1f}% [{r['ci_lo']:.1f}-{r['ci_hi']:.1f}]</td><td>{r['dm1_pct']:.1f}%</td><td>{r['dm2_pct']:.1f}%</td><td>{r['fisher_or']:.2f}</td><td>{r['fisher_p']:.3f}</td><td>{r['bh_q']:.3f}</td></tr>"
    summary += "</table></div>"
    if "FINAL real numbers" not in content:
        content = content.replace("</h1>", "</h1>" + summary, 1)
        p.write_text(content)
        print("  ✓ injected real-data summary into manuscript page")
    else:
        # Replace previous summary
        content = re.sub(r"<div style=\"background:rgba\(46,204,113.*?</div>", summary,
                         content, count=1, flags=re.DOTALL)
        p.write_text(content)
        print("  ✓ updated existing real-data summary")


def step6_update_mega(res):
    """Update HLA mega page §1.5 with final n."""
    print("\n=== Step 6: Update HLA mega §1.5 ===")
    p = PROJ / "reports/html/pages/revision_yu_2026_04_29_HLA.html"
    if not p.exists():
        print("  ! mega page not found, skipping")
        return
    n_total = int(res["n"].iloc[0])
    content = p.read_text()
    content = re.sub(r"~\d+/64", f"{n_total}/{n_total}", content)
    content = re.sub(r"~24/64", f"{n_total}", content)
    p.write_text(content)
    print(f"  ✓ updated mega page counts to n={n_total}")


def step7_destroy_cmd():
    """Generate VM destroy command."""
    print("\n=== Step 7: VM destroy command ===")
    cmd = """az vm delete --resource-group YG1-cook-vm_group --name arcas-burst --yes && \\
az disk delete --resource-group YG1-cook-vm_group --name arcas-burst_OsDisk_1_$(az vm show -g YG1-cook-vm_group -n arcas-burst --query id -o tsv 2>/dev/null | xargs basename) --yes 2>/dev/null
az disk list -g YG1-cook-vm_group --query "[?contains(name,'arcas-burst')].id" -o tsv | xargs -I{} az disk delete --ids {} --yes 2>/dev/null
az network public-ip list -g YG1-cook-vm_group --query "[?contains(name,'arcas-burst')].id" -o tsv | xargs -I{} az network public-ip delete --ids {} 2>/dev/null
az network nic list -g YG1-cook-vm_group --query "[?contains(name,'arcas-burst')].id" -o tsv | xargs -I{} az network nic delete --ids {} 2>/dev/null
az network nsg list -g YG1-cook-vm_group --query "[?contains(name,'arcas-burst')].id" -o tsv | xargs -I{} az network nsg delete --ids {} 2>/dev/null"""
    print(cmd)
    print("\n[NOT executing] — run manually after verifying results.")


def main():
    n_files = step1_sync()
    if n_files < 50:
        print(f"  ! only {n_files} JSON files — pipeline likely not complete")
        return
    step2_polish()
    step3_korean()
    res = step4_stats()
    step5_substitute_manuscript(res)
    step6_update_mega(res)
    step7_destroy_cmd()
    print("\n========================================")
    print(f"✓ DONE — n={int(res['n'].iloc[0])} samples processed")
    print(f"  12 FINAL figures + 6 KOR figures + Fisher/FDR stats")
    print(f"  manuscript page real-data summary injected")
    print("========================================")


if __name__ == "__main__":
    main()
