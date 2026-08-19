# DM1 한글 스토리 웹 (Vite + React + TS)

DM1 갑상선암 발견을 5–7 분에 설명하는 한글 스크롤 프레젠테이션.

## 실행

```bash
cd /home/seungho/personal/THCA_data_analysis/project/dm1_story_web
npm install                 # 최초 1회만
npm run dev                 # 0.0.0.0:5173 로 바인딩
```

`vite.config.ts` 의 `server.host = "0.0.0.0"` · `port = 5173` 으로 외부 IP 접근 가능.

## 접속

- 내부: <http://localhost:5173>
- 외부: `http://<VM-EXTERNAL-IP>:5173`
- 현재 VM 외부 IP: `40.82.129.113`  ·  최종 URL: <http://40.82.129.113:5173>

## 방화벽 체크리스트

- [ ] Azure NSG 에 inbound TCP 5173 허용 (CLAUDE.md 의 `Allow-Dev-Ports` 규칙에 이미 포함).
- [ ] `ufw status` — Ubuntu 방화벽이 5173 차단하지 않는지 확인.
- [ ] 프로세스가 `0.0.0.0:5173` 에 바인딩됐는지: `ss -tlnp | grep 5173`.
- [ ] 동일 VPN / 사설망에서 접속하는지 (public IP 차단 환경이면).

## 그림 자산

`public/figures/` 에 다음 파일이 있어야 합니다 (이미 복사됨):
- `Fig1_annotated.png` (Discovery composite)
- `Fig1_panel_D.png` (P_DM1 sorted heatmap)
- `Fig3_annotated.png` (Epigenetic silencing composite)
- `Fig5_annotated.png` (Pooled survival composite)
- `Fig6_annotated.png` (Reflex pathway composite)

원본 위치: `project/papers_hub_2026_05_04/assets/paper1_nc/`

새 그림을 교체하려면 같은 파일명으로 덮어쓰기만 하면 됨.

## 구조

```
src/
  main.tsx
  App.tsx                  ← 6 섹션 한글 스토리
  styles.css               ← clean white theme · Noto Serif KR
  components/
    StorySection.tsx       ← 섹션 래퍼 (앵커 + 헤드라인 + lead + body)
    FigureCard.tsx         ← 그림 카드 (lightbox 클릭)
    AxisDiagram.tsx        ← DM2 ↔ DM1 분화 축 SVG
    WorkflowDiagram.tsx    ← surgery → FFPE → 8-gene → DM call SVG
    DarkMatterDiagram.tsx  ← canonical drivers → dark matter → DM1 axis SVG
    DriverLayerDiagram.tsx ← driver 변이 × DM 상태 two-layer SVG
public/
  figures/
    *.png                  ← 실제 figure 자산
```

## 빌드 (정적 배포용)

```bash
npm run build               # dist/ 생성
sudo cp -r dist/* /var/www/papers/dm1_story_web/
```

`/var/www/papers/dm1_story_web/` 경로로 nginx 정식 배포 시 port 80 으로도 접근 가능.

## 핵심 보장 사항

- **모든 수치는 manuscript audit-locked 분석 기준** (현재 분석 기준 라벨).
- DM1 은 **risk stratification axis** 로만 표현 — treatment-selection biomarker 로 과장하지 않음.
- Limitation box 에 prospective validation 전 단계임을 명시.
- 외부 공유 전 PI 검토 필요.
