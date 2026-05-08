# 수동으로 받아주셔야 할 데이터 목록

**Lumenix · Neoantigen Vaccine Hub · 2026-05-07**
**연락처**: Seungho Cook · kukshomr@gmail.com

---

## 📊 현재 상태 (자동 수집 + 사용자 업로드 완료분)

| Source | 상태 | n records | License |
|---|---|---|---|
| ✅ **TESLA mmc.zip** (사용자 업로드) | 완료 | 918 (608 mmc4 + 310 mmc7) | Cell Press · derived CC-BY |
| ✅ CEDAR API | 자동 완료 | 5,000 | CC-BY |
| ✅ IEDB T-cell full v3 | 자동 완료 (1.3 GB CSV) | 100,000 (570k subset) | CC-BY 4.0 |
| ✅ **TSNAdb v2 validated** | 자동 완료 | 1,856 (L4 HELD_OUT) | academic |
| ✅ TSNAdb v2 predicted (전체 cancer) | 자동 완료 | 31,552 (L0) | academic |
| ✅ NeoRanking Gartner | 자동 완료 | 46 (L5 patient-matched) | academic |
| ✅ BigMHC GitHub examples | partial | 9 | CC-BY 4.0 |
| ✅ VDJdb 2023 | 자동 완료 | 62,178 TCR-pMHC | CC-BY-NC 4.0 |
| ✅ PRIME / MixMHCpred GitHub | 자동 완료 (binary 제외) | code only | academic |
| ✅ Clinical vaccine evidence | curated | 3 (L6) | derived CC-BY |
| 🔄 **Neodb** (Zenodo 16892216) | 다운 중 (884 MB) | TBD | check |

**현재 master**: 139,384 records · L5+ patient-matched 964 · L4+ immunogenic-credible 102,857 · leakage-free benchmark 2,823

---

## ❌ 못 받은 것 — **직접 받아주세요**

### 🔴 1순위 · IMPROVE (Borch 2024 Frontiers Immunology) — 가장 큰 임팩트

**왜 필요한가**: 17,500 tested neoepitope candidates + **467 T-cell-recognized** + 70 patients. 큰 patient-matched T-cell screen — TESLA 다음으로 큰 L5 cohort.

**검색 방법**:
- Google Scholar: `"IMPROVE" Borch 2024 neoantigen frontiers immunology`
- 또는 PubMed: `Borch IMPROVE neoantigen T-cell`
- DOI 추정: `10.3389/fimmu.2024.XXXXXXX` (정확한 번호 모름)

**다운 절차**:
1. 정확한 paper page 찾기 (frontiersin.org)
2. "Supplementary Material" / "Data Sheet 1" / "Data Sheet 2" 등 모든 첨부파일 다운로드
3. 보통 .xlsx 또는 .zip 형식

**드롭 위치**: `/data/neoantigen_vaccine_hub/data_raw/improve/`

---

### 🔴 2순위 · NEPdb (positive + negative labels) — class balance 결정적

**왜 필요한가**: 우리 master에 **negative label이 매우 부족**. TSNAdb validated는 거의 다 positive. NEPdb는 명시적 positive AND negative T-cell-validated neoepitopes 둘 다 가짐. class balance 정상화에 critical.

**다운 절차**:
1. https://nep.whu.edu.cn/ 접속 (사이트 reachable, 한국어/영어 다 됨)
2. 네비게이션에서 "Download" / "Bulk Data" / "Statistics" 메뉴 찾기
3. 등록 form (이름 + 이메일 + 연구목적) 작성
4. TSV/CSV 받기

**원하는 파일명** (예시):
- `nepdb_all.tsv` 또는 `nepdb_positive.tsv` + `nepdb_negative.tsv`
- 컬럼 예상: peptide, hla, mhc_class, tumor_type, gene, mutation, t_cell_response (positive/negative), reference

**드롭 위치**: `/data/neoantigen_vaccine_hub/data_raw/nepdb/`

---

### 🟡 3순위 · McPAS-TCR — TCR-pMHC pathology layer

**왜 필요한가**: pathology-associated TCR + autoreactivity (Hashimoto/Graves) cross-reactivity 검사. 갑상선 vaccine 안전성 검증 layer.

**다운 절차**:
1. https://friedmanlab.weizmann.ac.il/McPAS-TCR/ 접속
2. "Download data" 클릭
3. 이메일 + 이름 + 연구목적 form 작성
4. CSV 받기

**License 주의**: academic use only — 외부 redistribute 금지.

**드롭 위치**: `/data/neoantigen_vaccine_hub/data_raw/tcr/mcpas/`

---

### 🟡 4순위 · dbPepNeo (거의 dead site) — 시도만

**현 상태**: 원래 URL `www.biostatistics.online/dbPepNeo/` DNS 끊어짐. Wayback Machine은 index 페이지만 있고 zip 파일 archive는 없음.

**시도 옵션**:
- GitHub 검색: https://github.com/search?q=dbPepNeo+neoantigen (mirror 있을 수 있음)
- 논문 corresponding author 이메일 (dbPepNeo paper에서 contact info)
- 또는 **TSNAdb validated에 dbPepNeo merged sources 이미 포함** — TSNAdb로 부분 대체 됨

**원하는 파일** (찾으시면):
- `dbPepNeo_LC.zip` (LC = MS only)
- `dbPepNeo_MC.zip` (MC = MS + WES/WGS)
- `dbPepNeo_HC.zip` (HC = T-cell-validated) ← 가장 중요

**드롭 위치**: `/data/neoantigen_vaccine_hub/data_raw/dbpepneo/`

---

### 🟢 5순위 · BigMHC Mendeley bulk — 매우 큰 파일 (선택)

**왜 필요한가**: BigMHC paper의 진짜 training data. 수백만 rows MS-eluted ligand + immunogenicity + MANAFEST.

**다운 절차**:
1. https://data.mendeley.com/datasets/dvmz6pkzvb/4 접속
2. "Download all files" 또는 개별 파일 다운로드
3. CSV 파일들:
   - `el_train.csv` · `el_val.csv` · `el_test.csv` (MS eluted ligand)
   - `im_train.csv` · `im_val.csv` · `im_test.csv` (immunogenicity)
   - `manafest.csv` (T-cell screen) ← 가장 중요
   - `iedb.csv` · `pseudoseqs.csv` · `summary.csv`

**License**: CC-BY 4.0 (Mendeley)
**드롭 위치**: `/data/neoantigen_vaccine_hub/data_raw/bigmhc/`

---

## 📤 받으신 후 보내주는 방법 (편한 것 하나)

### (A) Claude Desktop drag-and-drop
파일을 채팅에 직접 첨부 (xlsx, csv, zip 다 가능)

### (B) base64 paste (작은 파일)
```bash
cd ~/Downloads
zip -r my_data.zip improve_suppl.xlsx nepdb_all.csv mcpas.csv
base64 my_data.zip | head -c 200000
# 채팅에 paste
```

### (C) scp / sftp (서버 접근 권한 있으시면)
```bash
# 1순위: IMPROVE
scp ~/Downloads/improve_suppl.xlsx \
    seungho@40.82.129.113:/data/neoantigen_vaccine_hub/data_raw/improve/

# 2순위: NEPdb
scp ~/Downloads/nepdb_all.tsv \
    seungho@40.82.129.113:/data/neoantigen_vaccine_hub/data_raw/nepdb/

# 3순위: McPAS-TCR
scp ~/Downloads/McPAS-TCR.csv \
    seungho@40.82.129.113:/data/neoantigen_vaccine_hub/data_raw/tcr/mcpas/

# 5순위: BigMHC bulk
scp ~/Downloads/manafest.csv \
    seungho@40.82.129.113:/data/neoantigen_vaccine_hub/data_raw/bigmhc/
```

### (D) 또는 그냥 같은 위치 (project 폴더 root)에 zip 던져두세요
```bash
scp my_data.zip seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/
```

저번에 TESLA 처럼 `mmc.zip`을 거기 두셔서 발견했던 것처럼 — 그냥 `/home/seungho/personal/THCA_data_analysis/` 또는 `/data/neoantigen_vaccine_hub/_inbox/` 에 던져주세요.

---

## ⏱ 예상 시간

| Source | 받는 시간 | 기대 추가 records |
|---|---|---|
| IMPROVE | 5-10분 (논문 찾기 + suppl 다운) | +17,967 (17,500 + 467 L5) |
| NEPdb | 10분 (form + download) | +수천 (positive + negative) |
| McPAS-TCR | 5분 (form + download) | +수천 TCR-pMHC |
| dbPepNeo | 시도 (성공 가능성 낮음) | 옵션 |
| BigMHC bulk | 5분 (큰 파일이라 download만 시간) | 수백만 (training/test) |

**3개 (IMPROVE + NEPdb + McPAS-TCR) 만 받으셔도** master 139k → 160k+ · leakage-free benchmark 2,823 → 5,500+ · NComm-급 cohort 완성.

---

## 📋 자동 ingest 명령어 (받자마자 실행)

받으시는 즉시 (또는 사용자가 알려주시면):
```bash
# Master 재빌드
python3 /data/neoantigen_vaccine_hub/scripts/build_master_dataset.py
# Leakage audit 재계산
python3 /data/neoantigen_vaccine_hub/scripts/assign_test_set_safety.py
# Benchmark 재실행
python3 /data/neoantigen_vaccine_hub/scripts/benchmark_models.py --quick
# Reports 재생성
python3 /data/neoantigen_vaccine_hub/scripts/generate_report.py
# 외부 IP 배포
rsync -a /home/seungho/personal/THCA_data_analysis/project/papers_hub_2026_05_04/neoantigen_hub_data/ \
    /var/www/papers/papers_hub_2026_05_04/neoantigen_hub_data/
```

---

## 🔗 외부 IP 접속 (받은 즉시 업데이트됨)

```
http://40.82.129.113/papers_hub_2026_05_04/neoantigen_hub.html
http://40.82.129.113/papers_hub_2026_05_04/neoantigen_hub_data/leakage_audit.json
http://40.82.129.113/papers_hub_2026_05_04/neoantigen_hub_data/source_x_test_set_safety.csv
```

---

## ❓ 문제 생기면

- 컬럼 이름 다를 때 → 그냥 보내주시면 자동 detect 시도
- 형식이 xlsx/csv/tsv/zip 어느거든 OK
- 인코딩 (UTF-8 vs CP949) 자동 처리

---

**Seungho Cook · Lumenix · 2026-05-07**
