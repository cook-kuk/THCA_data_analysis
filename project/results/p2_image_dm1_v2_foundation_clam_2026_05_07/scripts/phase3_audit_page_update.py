"""Phase 3b — Update Paper 1 one-page audit HTML with image-DM1 v2 results.

Reads SPRINT_RESULT.md + Phase 1 report + Phase 2 report.
Inserts a new section "★★★ Image-DM1 v2 sprint results" into p1_onepage_audit.html.
"""
from __future__ import annotations
import argparse, re
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--audit_html", type=Path,
                   default=Path("project/manuscript_v8/p1_onepage_audit.html"))
    p.add_argument("--sprint_root", type=Path,
                   default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"))
    args = p.parse_args()

    sprint_result = (args.sprint_root/"phase3_integration"/"SPRINT_RESULT.md")
    if not sprint_result.exists():
        print(f"[skip] {sprint_result} not found"); return
    body = sprint_result.read_text()

    # Build HTML section
    html_section = f"""
<!-- IMG-DM1-V2-SECTION-START -->
<section id="image-dm1-v2" style="margin:38px 0;padding:18px 22px;border:2px solid #c0392b;border-radius:8px;background:#fff8f0">
<h2 style="color:#7B1F2A;border-bottom:3px solid #7B1F2A;padding-bottom:6px">★★★ Image-DM1 v2 sprint results — foundation model + CLAM</h2>
<div class="lede">2026-05-07 sprint: closure NO-GO re-opened with foundation-model (UNI/ViT-L) + CLAM gated-attention MIL. ResNet50 closure baseline ~0.55.</div>
<pre style="background:#fafafa;border-left:4px solid #7B1F2A;padding:10px 14px;font-family:'JetBrains Mono',monospace;font-size:13px;white-space:pre-wrap">{body}</pre>
</section>
<!-- IMG-DM1-V2-SECTION-END -->
"""

    audit_html = args.audit_html.read_text(encoding="utf-8")
    # remove old block
    audit_html = re.sub(r"<!-- IMG-DM1-V2-SECTION-START -->.*?<!-- IMG-DM1-V2-SECTION-END -->\n?", "", audit_html, flags=re.S)
    # insert after section id="hero"
    m = re.search(r"</section>\s*<div class=\"divider-band\"></div>", audit_html)
    if m:
        idx = m.end()
        audit_html = audit_html[:idx] + "\n" + html_section + audit_html[idx:]
    args.audit_html.write_text(audit_html, encoding="utf-8")
    print(f"[done] inserted v2 section into {args.audit_html}")

if __name__ == "__main__":
    main()
