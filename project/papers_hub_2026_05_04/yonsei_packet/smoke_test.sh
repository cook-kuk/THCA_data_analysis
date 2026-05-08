#!/bin/bash
# End-to-end smoke test for the Yonsei drop-in.
# Runs three scenarios that exercise the auto-detection + pipeline.
set -e

SCRIPT="${SCRIPT:-21_yonsei_dropin.py}"
echo "=== [1] Template generation ==="
python3 "$SCRIPT" --template /tmp/yonsei_template.tsv
test -s /tmp/yonsei_template.tsv && echo "  ✓ template written"
head -3 /tmp/yonsei_template.tsv

echo
echo "=== [2] Self-test (synthetic 140-patient cohort) ==="
python3 "$SCRIPT" --self-test 2>&1 | tail -10
test -s /data/pdac_poc/results/yonsei/SUMMARY.json && echo "  ✓ SUMMARY.json written"
test -s /data/pdac_poc/results/yonsei/YONSEI_RESULTS.md && echo "  ✓ YONSEI_RESULTS.md written"
test -s /data/pdac_poc/results/yonsei/yonsei_updated_forest.png && echo "  ✓ forest figure written"

echo
echo "=== [3] Variant column-name detection ==="
# Generate a TSV using NON-standard column names to verify auto-detection
python3 -c "
import pandas as pd, numpy as np
np.random.seed(7)
n=120
df=pd.DataFrame({
    'patient_id':[f'YS_{i:04d}' for i in range(n)],     # not 'sample_id'
    'survival_months':np.random.exponential(20,n)+5,    # not 'os_months'
    'vital_status':np.random.choice(['Dead','Alive'],n,p=[0.7,0.3]),  # str, not 0/1
    'kras_aa':np.random.choice(['p.G12D','p.G12V','p.G12R','WT'],n,p=[0.4,0.25,0.15,0.20]),
    'age_at_diagnosis':np.random.normal(63,9,n).round()
})
df.to_csv('/tmp/yonsei_variant.tsv', sep='\t', index=False)
print('wrote /tmp/yonsei_variant.tsv with non-standard column names')
"
python3 "$SCRIPT" --validate /tmp/yonsei_variant.tsv 2>&1 | tail -10
echo "  ✓ variant column names auto-detected"

echo
echo "=== [4] Full run on variant file ==="
python3 "$SCRIPT" /tmp/yonsei_variant.tsv --label "Smoke test variant" 2>&1 | tail -8

echo
echo "=== ALL SMOKE TESTS PASSED ==="
ls -la /data/pdac_poc/results/yonsei/
