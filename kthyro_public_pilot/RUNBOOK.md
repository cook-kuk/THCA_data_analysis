# Runbook

1. Audit environment and local files:

```bash
python scripts/00_audit_environment.py
python scripts/01_audit_local_data.py
```

2. Link/download TCGA-THCA and score bulk axes:

```bash
python scripts/02_download_or_link_tcga_thca.py
python scripts/03_tcga_bulk_vulnerability.py
```

3. Link/download spatial thyroid data and score spatial axes:

```bash
python scripts/04_download_or_link_spatial_thyroid.py
python scripts/05_spatial_vulnerability_mapping.py
```

4. Run optional reference/drug/integration/report stages:

```bash
python scripts/06_scrna_reference_audit.py
python scripts/07_depm_prism_drug_pilot.py
python scripts/08_integrate_evidence.py
python scripts/09_make_proposal_figures.py
python scripts/10_write_report.py
```

5. End-to-end:

```bash
bash scripts/run_all.sh
```
