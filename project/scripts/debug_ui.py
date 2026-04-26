from pathlib import Path
import json,re
root=Path('project/reports/html')
files=[root/'index.html', *sorted((root/'pages').glob('*.html'))]
out=[]
for p in files:
    t=p.read_text(encoding='utf-8')
    out.append({'file':str(p),'has_viewport_fit':'viewport-fit=cover' in t,'has_nav_js':'assets/js/nav.js' in t,'has_mobile_js':'assets/js/mobile.js' in t,'has_filters_js':'assets/js/filters.js' in t,'has_downloads_js':'assets/js/downloads.js' in t,'uses_local_plotly':'assets/vendor/plotly/' in t})
print(json.dumps(out,ensure_ascii=False,indent=2))
