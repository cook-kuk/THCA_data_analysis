# Therapy Spectrum Triage

Date: 2026-05-11

This folder extends the GEO therapy-spectrum triage into a larger modality map and adds query templates for automation.

Files:
- `therapy_spectrum_triage.tsv` - expanded modality triage table
- `therapy_spectrum_query_templates.tsv` - GEO query templates for the agent screener

Snapshot:
- Rows: 16
- GO / WATCH / RESERVE: 10 / 5 / 1
- Top total-score lanes: neoantigen / mRNA vaccine / ICI combo; resistance / combo prioritization / response modeling; CAR-T / CAR-NK / CAR-NKT / TIL / TCR-T

Operating rule:
- Frame the platform as therapy triage + response prediction + biomarker selection + potency/QC + combo prioritization.
- Promote only studies with intervention + response readout + deployable output.
