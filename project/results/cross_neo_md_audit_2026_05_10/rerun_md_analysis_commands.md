# Rerun Commands

```bash
# Update live status and regenerate all currently possible MD outputs.
python project/scripts/cross_neo_md/00_parse_openmm_status.py
python project/scripts/cross_neo_md/01_annotate_complex.py
python project/scripts/cross_neo_md/02_md_qc.py
python project/scripts/cross_neo_md/03_analyze_pmhc_contacts.py
python project/scripts/cross_neo_md/04_analyze_tcr_contacts.py
python project/scripts/cross_neo_md/05_counterfactual_md_analysis.py
python project/scripts/cross_neo_md/06_replicate_consistency.py
python project/scripts/cross_neo_md/07_md_evidence_score.py
python project/scripts/cross_neo_md/08_integrate_md_with_cross_neo.py
python project/scripts/cross_neo_md/09_generate_md_figures.py
python project/scripts/cross_neo_md/11_generate_md_extra_visuals.py
python project/scripts/cross_neo_md/10_build_md_visual_dossier.py
python project/scripts/cross_neo_md/12_write_md_decision_report.py

# While 6VRN is running, sync state only.
OUT=project/results/cross_neo_md_audit_2026_05_10/remote_sync/pod1_thca_neo_bayesian_aux
mkdir -p "$OUT"
ssh -i ~/.runpod/ssh/RunPod-Key-Go -p 20878 root@135.84.176.142 \
  'cd /workspace/openmm_pilot_10ns_package && tar --exclude=trajectory.dcd --exclude=final.chk -czf - prod_10ns_6VRN_1fs300K' \
  | tar -xzf - -C "$OUT"

# After 6VRN finishes, sync full trajectory and rerun the analysis stack above.
rsync -avP -e 'ssh -i ~/.runpod/ssh/RunPod-Key-Go -p 20878' \
  root@135.84.176.142:/workspace/openmm_pilot_10ns_package/prod_10ns_6VRN_1fs300K/ \
  project/results/cross_neo_md_audit_2026_05_10/remote_sync/pod1_thca_neo_bayesian_aux/prod_10ns_6VRN_1fs300K/
```
