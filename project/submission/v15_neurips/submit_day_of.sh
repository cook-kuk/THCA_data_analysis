#!/usr/bin/env bash
# v15 NeurIPS — submission day-of preflight script
# Usage:  bash submit_day_of.sh
# Run from /opt/thyroid-dash/project/.
# Performs final compile + 8 sanity gates + checksum.
# Exit code 0 = ready to submit; non-zero = block.

set -e

PROJECT="/opt/thyroid-dash/project"
SUBMIT="$PROJECT/reports/v15_neurips/submit"
PKG="$PROJECT/submission/v15_neurips"
PDF_FINAL="$PKG/v15_neurips_SUBMIT_FINAL.pdf"

cd "$SUBMIT" || { echo "ERR: cannot cd to $SUBMIT"; exit 1; }

echo "============================================================"
echo "v15 NeurIPS submission day-of preflight"
echo "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo "============================================================"

# ---------- gate 1: anonymous repo URL is filled (not placeholder) ----
if grep -q 'v15-dial-XXXXXX' v15_missing_drafts_v2.tex; then
  echo "❌ FAIL gate 1: anonymous URL still placeholder. Replace 'v15-dial-XXXXXX' in v15_missing_drafts_v2.tex with the real anonymous.4open.science URL."
  exit 1
fi
echo "✓ gate 1: anonymous URL filled"

# ---------- gate 2: clean compile ------------------------------------
echo "→ recompiling (pdflatex; bibtex; pdflatex; pdflatex)..."
pdflatex -interaction=nonstopmode v15_neurips_SUBMIT.tex >/tmp/g2_p1.log 2>&1
bibtex v15_neurips_SUBMIT >/tmp/g2_b.log 2>&1 || true  # bibtex returns nonzero on warnings
pdflatex -interaction=nonstopmode v15_neurips_SUBMIT.tex >/tmp/g2_p2.log 2>&1
pdflatex -interaction=nonstopmode v15_neurips_SUBMIT.tex >/tmp/g2_p3.log 2>&1

if grep -q '! LaTeX Error' v15_neurips_SUBMIT.log; then
  echo "❌ FAIL gate 2: LaTeX Error in compile. Check /tmp/g2_*.log."
  grep -A2 'LaTeX Error' v15_neurips_SUBMIT.log | head -10
  exit 1
fi
if grep -qE 'Reference .* undefined|Citation .* undefined' v15_neurips_SUBMIT.log; then
  echo "❌ FAIL gate 2: undefined reference / citation."
  grep -E 'Reference .* undefined|Citation .* undefined' v15_neurips_SUBMIT.log | head -5
  exit 1
fi
echo "✓ gate 2: compile clean (no errors, no undefined refs)"

# ---------- gate 3: \TODO macros gone --------------------------------
TODO_COUNT=$(grep -c '\\TODO' v15_neurips_SUBMIT.tex theorem2_proof_v2.tex v15_missing_drafts_v2.tex checklist_v15_filled.tex 2>/dev/null | awk -F: '{s+=$2} END {print s}')
if [ "${TODO_COUNT:-0}" -gt 0 ]; then
  echo "❌ FAIL gate 3: $TODO_COUNT \\TODO markers remain"
  exit 1
fi
echo "✓ gate 3: 0 \\TODO markers"

# ---------- gate 4: anonymisation grep --------------------------------
python3 - <<'PYEOF' || { echo "❌ FAIL gate 4: author identifier leaked into PDF"; exit 1; }
from pdfminer.high_level import extract_text
import sys
txt = extract_text("v15_neurips_SUBMIT.pdf")
leaks = []
for m in ['Seungho','kukshomr','Cornerstone','SNUBH','seungho']:
    if m in txt: leaks.append(m)
if leaks:
    print("LEAKS:", leaks)
    sys.exit(1)
PYEOF
echo "✓ gate 4: anonymisation clean (text extract)"

# ---------- gate 5: PDF metadata anonymous ---------------------------
META=$(strings -n 5 v15_neurips_SUBMIT.pdf | grep -iE 'Author|Title|Subject|Keywords' | grep -vE '\(\)' | head -3)
if [ -n "$META" ]; then
  echo "⚠️ gate 5 WARN: PDF metadata may have non-empty fields:"
  echo "$META"
  echo "Consider running:  exiftool -all= -overwrite_original v15_neurips_SUBMIT.pdf"
fi
echo "✓ gate 5: PDF metadata check (warn-only)"

# ---------- gate 6: page count within limit --------------------------
PAGES=$(grep -E 'Output written.*pages' v15_neurips_SUBMIT.log | tail -1 | grep -oE '[0-9]+ pages' | awk '{print $1}')
echo "→ PDF has $PAGES pages total (NeurIPS 2026: main ≤10 pages + checklist + unlimited appendix)"
if [ "${PAGES:-0}" -gt 25 ]; then
  echo "⚠️ gate 6 WARN: $PAGES pages is unusually high; verify main body still ≤10"
fi
echo "✓ gate 6: page count $PAGES (informational)"

# ---------- gate 7: NeurIPS Paper Checklist present ------------------
python3 - <<'PYEOF' || { echo "❌ FAIL gate 7: checklist section missing from PDF"; exit 1; }
from pdfminer.high_level import extract_text
import sys
txt = extract_text("v15_neurips_SUBMIT.pdf")
if "NeurIPS Paper Checklist" not in txt:
    print("Checklist section NOT FOUND in PDF text — desk-reject risk!")
    sys.exit(1)
PYEOF
echo "✓ gate 7: NeurIPS Paper Checklist present in PDF"

# ---------- gate 8: copy + checksum ----------------------------------
cp -f v15_neurips_SUBMIT.pdf "$PDF_FINAL"
SHA=$(sha256sum v15_neurips_SUBMIT.pdf | awk '{print $1}')
SIZE=$(stat -c '%s' v15_neurips_SUBMIT.pdf)
echo "✓ gate 8: copied to $PDF_FINAL"
echo "  sha256: $SHA"
echo "  bytes : $SIZE"

# ---------- gate 9: anonymous_code.zip exists ------------------------
if [ ! -f "$PKG/supplementary/anonymous_code.zip" ]; then
  echo "❌ FAIL gate 9: $PKG/supplementary/anonymous_code.zip missing"
  exit 1
fi
ZIP_SHA=$(sha256sum "$PKG/supplementary/anonymous_code.zip" | awk '{print $1}')
ZIP_SIZE=$(stat -c '%s' "$PKG/supplementary/anonymous_code.zip")
echo "✓ gate 9: supplementary zip exists ($ZIP_SIZE bytes, sha256 ${ZIP_SHA:0:16}...)"

# ---------- final summary --------------------------------------------
echo ""
echo "============================================================"
echo "🟢 ALL GATES PASS — READY TO SUBMIT TO OPENREVIEW"
echo "============================================================"
echo ""
echo "Submit to https://openreview.net (NeurIPS 2026 venue)"
echo "Phase 1 abstract deadline: May 4, 2026 AOE"
echo "Phase 2 paper deadline:    May 6, 2026 AOE"
echo ""
echo "Files to upload:"
echo "  Main PDF:        $PDF_FINAL"
echo "  Supplementary:   $PKG/supplementary/anonymous_code.zip"
echo ""
echo "Form-field text:   $PKG/openreview_form_fields.md"
echo ""
echo "After submission, take screenshot of submission ID + email."
