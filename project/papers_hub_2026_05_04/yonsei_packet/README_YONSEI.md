# 연세대 PDAC 코호트 통합 안내서 (One-pager)

이 안내서대로 따라하시면 **n ≥ 100 환자 데이터를 받자마자 1분 안에** 메타분석 + 논문용 그림 + 원고 본문이 자동 생성됩니다.

---

## TL;DR · 가장 빠른 길

```bash
# 1. 패킷 받으시고 압축 풀어주세요 (lumenix_yonsei_packet.tar.gz)
tar xvf lumenix_yonsei_packet.tar.gz
cd lumenix_yonsei_packet

# 2. Python 환경 (3.10 이상)
pip install -r requirements.txt

# 3. (선택) 빈 템플릿 받기
python 21_yonsei_dropin.py --template yonsei_template.tsv

# 4. 데이터 채우신 다음 (또는 실데이터 파일 그대로 쓰셔도 됩니다 — 컬럼명 자동인식)
python 21_yonsei_dropin.py yonsei_real.tsv --label "Yonsei PDAC 2026"

# 끝. 1분 안에 다음이 생성됩니다:
#   results/yonsei/SUMMARY.json
#   results/yonsei/YONSEI_RESULTS.md     ← 원고 한 단락 그대로
#   results/yonsei/yonsei_updated_forest.png/svg     ← 메타-forest 그림
```

---

## 1. 필수 5컬럼 (대소문자·언더바 무관)

| 컬럼 | 자료형 | 설명 | 자동인식되는 다른 이름들 |
|---|---|---|---|
| `sample_id` | 문자열 | 환자 ID (익명화됨) | `sample`, `patient_id`, `case_id`, `tcga_barcode` … |
| `os_months` | 실수 | 전체생존기간 (개월) | `os_time`, `survival_months`, `overall_survival`, `follow_up_months` … |
| `os_event` | 정수 | 1 = 사망, 0 = 검열 | `vital_status`, `death`, `dead_alive`, `status` … (문자열도 자동변환) |
| `kras_allele` | 문자열 | `G12D`, `G12V`, `G12R`, `G12C`, `G13D`, `Q61H`, `Q61R`, `WT` 중 하나 | `kras`, `kras_aa`, `kras_mutation`, `mutation`, `protein_change` …<br>(`p.G12D` 형식, `KRAS p.G12D` 형식, `G12D` 형식 모두 자동변환) |
| `age` | 실수 | 진단 시 나이 (년) | `age_at_diagnosis`, `age_dx`, `diagnosis_age` … |

## 2. 선택 컬럼 (있으면 분석이 풍부해짐)

| 컬럼 | 풀어주는 분석 |
|---|---|
| `moffitt_call` | basal/classical 분류 → 4-cohort paradox 검증 |
| `ancestry` | Korean / Asian / European → HLA-bias 맥락 |
| `hla_a`, `hla_b`, `hla_drb1` 등 | 환자별 mut × HLA actionability 직접 검증 |

## 3. 입력 파일 포맷

| 포맷 | 자동지원? |
|---|---|
| TSV (`.tsv`) | ✅ 권장 |
| CSV (`.csv`) | ✅ |
| Excel (`.xlsx` / `.xls`) | ✅ openpyxl 사용 |
| MAF | ✅ `--maf-mode` 자동변환 (KRAS missense만 추출) |

## 4. 사전 점검 (실데이터 받기 전 확인)

```bash
# 빈 템플릿 받아서 5행 sample 그대로 한번 돌려보세요
python 21_yonsei_dropin.py --template /tmp/check.tsv
python 21_yonsei_dropin.py --self-test    # 합성 140명으로 end-to-end 검증
```

## 5. 실데이터 받자마자 1분 안에

```bash
# (a) 형식만 점검
python 21_yonsei_dropin.py --validate yonsei.tsv

# (b) 통과하면 풀-파이프라인 실행
python 21_yonsei_dropin.py yonsei.tsv --label "Yonsei PDAC 2026"
```

산출물 (5개):

```
results/yonsei/
├── SUMMARY.json                  ← 기계판독용 (Cox table + meta-pool 결과)
├── YONSEI_RESULTS.md             ← 논문 원고 한 단락 그대로 복붙용
├── yonsei_updated_forest.png     ← Yonsei row 보라색 + pooled diamond
├── yonsei_updated_forest.svg     ← 출판용 vector
└── yonsei_input_validated.tsv    ← QA용 canonicalized 입력
```

## 6. Docker로 한방 (Python 환경 안 만지고)

```bash
docker build -t lumenix/pdac-yonsei .
docker run --rm -v $(pwd):/work lumenix/pdac-yonsei /work/yonsei.tsv --label "Yonsei PDAC 2026"
```

## 7. Privacy / IRB

- **개인식별정보는 필요 없음** — `sample_id`는 익명화된 임의 문자열이면 됩니다.
- 원하시면 RNA z-score는 **연세 내부에서 사전계산** 후 per-sample 요약통계만 공유 가능 (raw FASTQ/BAM 외부 반출 불필요).
- Docker 모드는 완전 오프라인 — 외부 네트워크 접근 0회.

## 8. 자주 묻는 질문

**Q. 컬럼명이 우리 LIS 시스템 표준이라 다른데?**
A. 위 표의 "자동인식되는 다른 이름들" 리스트에 들어 있으면 그대로 인식됩니다. 그래도 안 되면 `--validate` 모드로 어떤 이름이 인식되는지 확인 후 헤더만 수정해주세요.

**Q. KRAS 변이가 `c.35G>A`처럼 cDNA 표기로 와 있어요.**
A. 사전에 protein change(`G12D`)로 변환해주세요. cDNA → protein 변환은 변이 콜러 (Strelka2/Mutect2) 출력에 보통 같이 나옵니다.

**Q. n이 100 미만이면?**
A. 50명까지는 Cox 적합 가능 (각 대립유전자 카운트가 ≥ 5만 되면). 30명 미만은 logrank만 출력합니다.

**Q. OS가 일수(days) 단위인데?**
A. `days_to_last_followup` 같은 일수 컬럼을 자동인식하면 자동으로 30.4375로 나누어 개월로 변환합니다.

**Q. `os_event`가 "Dead"/"Alive" 같은 문자열인데?**
A. 자동변환됩니다 (`Dead`, `dead`, `1` → 1, `Alive`, `alive`, `0`, `Censored` → 0).

**Q. 실패하면?**
A. 첫줄 에러 메시지 + 인식된 컬럼 리스트 + 어떤 컬럼이 빠졌는지 정확히 알려줍니다. 그대로 캡처해서 보내주세요.

## 9. 결과 해석

`YONSEI_RESULTS.md` 안에 다음이 자동 작성됩니다:

- 컬럼 자동인식 로그 (어떤 이름 → 어떤 표준컬럼으로 매핑되었는지)
- Yonsei 코호트 per-allele Cox 표 (HR · 95% CI · p)
- Updated meta-analysis (k=5 → k=6, n=3,113 → n=3,213+)
- **논문 원고용 한 단락** (그대로 복붙 가능)

## 10. 연락처

문제 / 질문 / 컬럼 매핑 충돌 시:
**Seungho Cook · Lumenix · kukshomr@gmail.com · 2026-05-07**

Slack/이메일로 첫 5행 (PII 가려서) 보내주시면 30분 안에 정확한 입력 명령어 회신드립니다.
