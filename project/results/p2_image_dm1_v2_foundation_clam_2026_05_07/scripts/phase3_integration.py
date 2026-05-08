"""Phase 3 — Local integration of Phase 1 + Phase 2 results."""
from __future__ import annotations
import argparse, re
from pathlib import Path
import pandas as pd

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path,
                   default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"))
    args = p.parse_args()
    p1 = args.root / "phase1_gse250521"; p2 = args.root / "phase2_tcga_clam"
    out = args.root / "phase3_integration"; out.mkdir(parents=True, exist_ok=True)

    # Phase 1 verdict
    p1_corr = p1 / "uni_dm1_correlation_per_slide.tsv"
    if p1_corr.exists():
        df1 = pd.read_csv(p1_corr, sep="\t"); pc1 = df1[df1["pc"]==1]
        n_slides = len(pc1); pass_n = (pc1["rho_DM1"].abs()>0.3).sum()
        pass_rate = pass_n / max(n_slides,1); max_rho = pc1["rho_DM1"].abs().max()
        p1_verdict = "PASS" if pass_rate>=0.5 else ("MARGINAL" if pass_rate>=0.25 else "FAIL")
    else:
        p1_verdict="MISSING"; n_slides=0; pass_n=0; pass_rate=0; max_rho=float("nan")

    # Phase 2 verdict
    p2_report = p2 / "PHASE2_REPORT.md"; overall_auc=float("nan"); p2_verdict="MISSING"
    if p2_report.exists():
        text = p2_report.read_text()
        m = re.search(r"Pooled cross-fold AUC:\s*([\d.]+)", text)
        if m: overall_auc = float(m.group(1))
        if overall_auc>0.70: p2_verdict="PASS"
        elif overall_auc>0.60: p2_verdict="MARGINAL"
        elif overall_auc>0: p2_verdict="FAIL"

    # Combined
    if p1_verdict=="PASS" and p2_verdict=="PASS": combined="PASS_LAUNCH_PAPER2"
    elif p1_verdict in {"PASS","MARGINAL"} and p2_verdict in {"PASS","MARGINAL"}: combined="MARGINAL_PAPER1_SUPP"
    elif p1_verdict=="PASS": combined="PARTIAL_P1_ONLY"
    elif p1_verdict=="FAIL": combined="CLOSURE_REAFFIRM_V2"
    else: combined="INCONCLUSIVE"

    rep = [
        f"# Sprint result — Foundation-model + CLAM v2  (2026-05-07)\n",
        f"## Verdicts\n",
        f"- Phase 1 (GSE250521 ↔ DM1): **{p1_verdict}**",
        f"  - n_slides={n_slides}; |ρ|>0.3 in PC1: {pass_n}/{n_slides} ({pass_rate:.0%}); max |ρ|={max_rho:.3f}",
        f"- Phase 2 (TCGA CLAM AUC): **{p2_verdict}**  AUC={overall_auc:.3f}",
        f"\n## Combined verdict\n\n**{combined}**\n",
    ]
    if combined=="PASS_LAUNCH_PAPER2":
        rep.append("→ Paper 1 spatial supp augment + Paper 2 image-DM1 launch unlocked.")
    elif combined=="MARGINAL_PAPER1_SUPP":
        rep.append("→ Paper 1 spatial supp augment; Paper 2 launch deferred.")
    elif combined=="PARTIAL_P1_ONLY":
        rep.append("→ Paper 1 supp only; Phase 2 did not pass.")
    elif combined=="CLOSURE_REAFFIRM_V2":
        rep.append("→ closure NO-GO re-confirmed at v2 architecture.")
    else:
        rep.append("→ inconclusive (one or more phases missing).")
    (out/"SPRINT_RESULT.md").write_text("\n".join(rep)); print("\n".join(rep))

if __name__ == "__main__":
    main()
