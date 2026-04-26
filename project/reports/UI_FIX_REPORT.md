# UI_FIX_REPORT

## detected root cause(s)
- Root cause 1: 모바일 메뉴 버튼 자체가 렌더링되지 않았음.
- Root cause 2: 홈페이지 상단 네비게이션의 상대 경로가 일부 잘못 생성됨.
- Root cause 3: 모바일 전용 drawer / backdrop / scroll lock / ESC close 로직이 없었음.
- Root cause 4: CDN 의존성이 강했고 로컬 fallback 사용이 일관되지 않았음.

## exact files changed
- project/reports/html/index.html
- project/reports/html/pages/*.html
- project/reports/html/assets/css/main.css
- project/reports/html/assets/js/main.js
- project/reports/html/assets/js/nav.js
- project/reports/html/assets/js/mobile.js
- project/reports/html/assets/js/filters.js
- project/reports/html/assets/js/downloads.js

## before/after summary
- Before: desktop 위주 정적 네비게이션, 모바일 메뉴 없음, 경로 취약.
- After: 모바일 햄버거 drawer, filter drawer, body scroll lock, share/theme/filter action, chart action buttons 강화.

## menu fix summary
- topnav 우선순위 z-index 상향
- decorative layer pointer-events 차단
- 모바일 메뉴 버튼 동적 주입
- backdrop click / ESC close / aria-expanded 반영
- 메뉴와 필터 drawer 동시 open 방지

## mobile compatibility summary
- 360px~1440px 기준 단일 컬럼 / 2열 / 3열 responsive 조정
- chart action row 모바일에서 inline overlay 대신 block row 전환
- touch target 44px 이상 확보

## remaining minor issues
- Playwright 브라우저 설치 실패 가능성 있음
- index 원문 HTML의 일부 잘못된 href는 JS에서 런타임 정규화함
- Kaleido PNG는 Chrome 미설치로 placeholder 파일이 남아 있음
