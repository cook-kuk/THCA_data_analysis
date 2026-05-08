#!/bin/bash
# Phase 1+2 결과 받은 후 로컬에서 한방에 통합:
#   1. tar 받아서 풀기
#   2. phase3_integration.py
#   3. phase3_audit_page_update.py
#   4. phase3_make_supp_figure.py
#   5. memory_update_proposal.md 출력 + 사용자 결정 대기

set -e
cd /home/seungho/personal/THCA_data_analysis
ROOT=project/results/p2_image_dm1_v2_foundation_clam_2026_05_07

if [ "$#" -gt 0 ]; then
  TAR_PATH="$1"
  echo "=== extract $TAR_PATH ==="
  tar xzf "$TAR_PATH" -C "$ROOT"
fi

echo ""
echo "=== Phase 3 integration ==="
.venv/bin/python "$ROOT/scripts/phase3_integration.py" --root "$ROOT"

echo ""
echo "=== Phase 3 supp figure build ==="
.venv/bin/python "$ROOT/scripts/phase3_make_supp_figure.py" --root "$ROOT"

echo ""
echo "=== Phase 3 audit page update (Paper 1 one-page audit) ==="
.venv/bin/python "$ROOT/scripts/phase3_audit_page_update.py" \
  --audit_html project/manuscript_v8/p1_onepage_audit.html \
  --sprint_root "$ROOT"

echo ""
echo "=== outputs ==="
ls -lh "$ROOT/phase3_integration/"

echo ""
echo "Next: review SPRINT_RESULT.md + supp_figure_image_dm1_v2.png + audit page"
