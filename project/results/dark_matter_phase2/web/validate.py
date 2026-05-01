"""Validate v2.html — check img paths, data paths, basic HTML structure."""
import re
from pathlib import Path

HTML = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase2/web/v2.html")
WEB = HTML.parent

content = HTML.read_text()
print(f"HTML size: {len(content):,} chars, {content.count(chr(10))} lines")
print(f"<img> tags: {content.count('<img ')}")
print(f"<table> tags: {content.count('<table')}")
print(f"<section> tags: {content.count('<section ')}")
print(f"<h2> tags: {content.count('<h2')}")
print(f"<h3> tags: {content.count('<h3')}")

# Tag balance check
for tag in ['div', 'section', 'table', 'tbody', 'thead', 'tr', 'script', 'style']:
    opens = len(re.findall(rf'<{tag}[\s>]', content))
    closes = len(re.findall(rf'</{tag}>', content))
    bal = "✓" if opens == closes else f"⚠ MISMATCH ({opens}-{closes}={opens-closes})"
    print(f"  <{tag}> open={opens}  close={closes}  {bal}")

# Check img src paths
print("\n=== Image asset check ===")
imgs = re.findall(r'<img\s+src="([^"]+)"', content)
missing = []
for src in imgs:
    if src.startswith("http"):
        continue
    fp = WEB / src
    if not fp.exists():
        missing.append(src)
        print(f"  ✗ MISSING: {src}")
    else:
        print(f"  ✓ {src}  ({fp.stat().st_size:,} bytes)")

print(f"\n{len(imgs)} img tags total, {len(missing)} missing")

# Check data href downloads
print("\n=== Data download asset check ===")
data_links = re.findall(r'href="(data/[^"]+)"', content)
for d in set(data_links):
    fp = WEB / d
    if fp.exists():
        print(f"  ✓ {d}  ({fp.stat().st_size:,} bytes)")
    else:
        print(f"  ✗ MISSING: {d}")

# Check fetch() data paths in JS
print("\n=== JS fetch() data paths ===")
fetches = re.findall(r"fetch\s*\(\s*['\"]([^'\"]+)['\"]", content)
for f in set(fetches):
    fp = WEB / f
    if fp.exists():
        print(f"  ✓ {f}  ({fp.stat().st_size:,} bytes)")
    else:
        print(f"  ✗ MISSING: {f}")

# Check all anchor #targets exist
print("\n=== Anchor target check ===")
anchors = re.findall(r'href="#([^"]+)"', content)
ids = set(re.findall(r'id="([^"]+)"', content))
broken = [a for a in set(anchors) if a not in ids and a != "top"]
if broken:
    print(f"  ⚠ Broken anchors: {broken}")
else:
    print(f"  ✓ All {len(set(anchors))} anchors point to existing IDs")

print("\n=== Final ===")
print(f"  Sections: {content.count('<section ')}")
print(f"  Tables: {content.count('<table')}")
print(f"  Figures (img): {len(imgs)}")
print(f"  Glossary terms (dt): {content.count('<dt>')}")
print(f"  References (ref-list li): {len(re.findall(r'<li[^>]*>', content))}")
