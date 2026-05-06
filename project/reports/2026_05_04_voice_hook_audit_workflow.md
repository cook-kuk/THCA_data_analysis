# Paper 1 Hook ¶1 — Writing Guide + Audit Workflow (본인 키보드용)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Mode:** Marathon scaffolding — voice-protected Hook ¶1 작성용 가이드 + 검토 workflow.
**Companion:** `2026_05_04_voice_hook_fact_brief.md` (fact brief, 동일 디렉터리).
**Hook ¶1 자체 작성 금지** — Claude는 본 문서로 가이드 제공만, 사용자가 paste 후 audit 실행.

---

## [1] Hook writing guide — bullets only

- 길이: ~120–180 words / 4 sentences (Cell Rep Med Intro ¶1 norm)
- Voice: 본인 (Korean academic English; OUTLINE ¶1 line 37 토글 "no mechanistic compass" vs "without molecular guidance" 본인 결정)
- Tense: present for clinical/biological facts; past for cohort/published findings
- Frame order (사용자 추천): clinical unmet need → BRAF/RAS limit → dark-matter gap → transcriptional/RAI-lineage axis 필요성
- Hook ¶1 boundary: clinical + driver-paradigm limit only. 다음 단락으로 hold —
  - ¶2 = existing molecular landscape (TCGA 2014 / Yoo 2016 / Landa 2016 / Pan 2025)
  - ¶3 = Korean/Asian HLA bridge
  - ¶4 = Aim ("Here we …")
- Result spoilers ¶1 금지: DM1/DM2 axis, sub-B NBNR, HT-overlap mechanism, HLA cohort 숫자, 5-pillar overview → 모두 ¶4 / Results 영역
- Cite envelope:
  - ATA 2015 (mandatory)
  - TCGA 2014 (Cancer Network) optional
  - **Landa 2016 JCI 126(3):1052–1066** preferred over Krishnamoorthy 2025 (`v17_landa2016_cite_save` 메모리: 4 fact 모두 misattribution 위험)
- Density cap: ≤3 numerical citations in Hook ¶1 (else reads like Results)
- Verb floor: `supports` / `consistent with` / `stratifies` / `complements` / `does not capture` — 금지: `proves` / `establishes` / `definitively` / `first` / `novel` (without literature anchor)
- Pre-write 체크: `references.bib`에 Landa 2016 / ATA 2015 / TCGA 2014 entries 존재 확인 → Discussion §3.1 framing alignment 후 paste

---

## [2] Recommended 4-sentence structure — sentence labels only

| 문장 | label | 내용 | 숫자 cite |
|---|---|---|---|
| **S1** | Clinical magnitude | DTC recurrence 15–35% + ATA 2015 anatomic stratification 3–5%↔50–75% range | ATA 2015 |
| **S2** | Driver paradigm utility | BRAF/RAS/TERT mutation-based stratification clinical role description (no critique yet) | (optional) TCGA 2014 |
| **S3** | Specific limitation | driver mutation status does not propagate to transcript-level expression OR partition the BRAF-neg ∩ RAS-neg subset OR predict RAI-lineage trajectory — pick **one** axis (not all three) | 1 numerical cite max from §3 list |
| **S4** | Bridge / motivation | transcriptional differentiation axis = the missing layer (no result, no Aim verb). "Here we" 금지 — that opens ¶4 | none |

---

## [3] Safe numbers (Hook only — pick ≤3 total)

| domain | 값 | 사용 권장 |
|---|---|---|
| DTC recurrence | **15–35%** | S1 |
| ATA 2015 risk range | **3–5%** (low) ↔ **50–75%** (high) | S1 |
| BRAF V600E vs WT transcript | **Cohen d = −0.044, p = 0.567** | S3 (cite n only if room) |
| Driver single-feature AUC ceiling | **BRAF 0.602, TERT 0.578** | S3 (BRAF d 사용 시 redundant) |
| Driver_anchor 12-gene cluster vs DM1/DM2 | **ARI = −0.007** | S3 (powerful, high return per word) |
| BRAF-neg ∩ RAS-neg subset size | **n = 156** (TCGA-THCA) | S3 |

**Density discipline:** S1 = 1 cite (ATA range), S3 = 1 numerical cite (pick one of d / AUC / ARI / n=156). S2/S4 may be number-free.

---

## [4] Avoid in Hook

| category | items |
|---|---|
| **Result spoilers** (R/¶4 territory) | DM1/DM2 axis, 8-gene panel composition, sub-A/sub-B, 94.6% mut-neg, Hashimoto-like 18–30%, DPB1\*05:01 53.2%, n=874, GSE286332 10,380 DEGs, IFN-γ FDR=2e-4, HLA-II d=+3.65, TIERA67 ARI=0.903, PFI HR=2.04 |
| **Closure-NO-GO territory** | H&E, WSI, image, morphology, pathology-AI, tile, spot-level, foundation model, image-DM1 closure metrics (r=0.022, AUROC=0.511) |
| **Out-of-paper (Paper 2/3/4)** | PTC+HT mechanism, BCR/TLS, cookHLA, arcasHLA, Bundang FFPE, GD / TSAb / exophthalmos, ICI / neoantigen / HLA loss, Pan-Asian HLA bridge framing |
| **Therapeutic / wet-lab** | TROP2, sacituzumab, IHC, organoid, drug overlay, actionability |
| **Speculative / forbidden** | Cancer Cell %, Nat Cancer reach, "all risks resolved", "first", "novel", "definitively", "establishes", venue-probability |
| **Cite risk** | Krishnamoorthy 2025 (4 fact misattribution risk) → substitute Landa 2016 JCI 126(3):1052–1066 |

---

## [5] After-user-draft audit checklist (붙여넣기 후 Claude 실행)

### A. Factual accuracy
- **A1:** 각 숫자 claim → §3 safe-numbers 리스트와 1:1 매치 (sign / decimal / percent 기호)
- **A2:** BRAF d = −0.044 사용 시 음수 부호 + 소수 3자리 보존 확인
- **A3:** ATA 2015 range 인용 시 3–5% / 50–75% 정확
- **A4:** 인용된 모든 reference → `project/reports/2026_05_03_references.bib` entry 존재 확인
- **A5:** §4 avoid-list 숫자 누락 침투 grep

### B. Overclaim risk
- **B1:** verb scan — `proves` / `establishes` / `definitively` / `demonstrates the framework` 등장 시 FLAG
- **B2:** "all risks resolved" / "모든 risk 해소" / "complete picture" / "fully explains" FLAG
- **B3:** `first` / `novel` 사용 시 prior literature anchor 없으면 FLAG
- **B4:** venue-probability / "ready for X journal" 표현 FLAG
- **B5:** TROP2 / sacituzumab / wet-lab / actionability 단어 등장 FLAG

### C. Scope contamination
- **C1:** H&E / WSI / pathology / image / morphology / tile / spot — 단어 검출 시 FLAG (Paper 2A NO-GO 영역)
- **C2:** PTC+HT mechanism / Hashimoto signature / TLS / BCR / IGHV — 등장 시 ¶1로는 부적합 FLAG (Paper 2B 또는 R5 territory)
- **C3:** ICI / neoantigen / HLA loss / immune ecotype / ICI-readiness — Paper 3 frozen territory FLAG
- **C4:** GD / Graves / TSAb / exophthalmos / Pan-Asian / cookHLA — Paper 4 backlog FLAG
- **C5:** Bundang / arcasHLA / RunPod / Azure / Foundation model — infra/outreach 누출 FLAG

### D. Voice / structure
- **D1:** 문장 수 ≤5; 단어 수 120–200 범위 확인
- **D2:** 숫자 인용 ≤3 (S1 ATA + 1 numerical = 권장 2; 절대 상한 3)
- **D3:** "Here we" / "We hypothesize" / "We propose" / "In this study" 시작 → ¶4 Aim 영역 침범 FLAG
- **D4:** Krishnamoorthy 2025 인용 검출 시 → Landa 2016 substitution 권고
- **D5:** ATA 2015 인용 누락 시 → 추가 권고

### E. Cross-section consistency
- **E1:** Discussion §3.1 framing tone과 모순 여부 (`v17_landa2016_cite_save` 메모리 + OUTLINE Discussion D1 영역과 정합)
- **E2:** Limitations §3.4 negative-disclosure 방향성과 모순 여부 (Hook이 "all resolved" 톤이면 §3.4와 충돌)
- **E3:** ¶4 Aim 문장과 중복되는 result claim 검출

### F. Audit output 형식
- 결과 코드: `PASS` / `FLAG-FACT` / `FLAG-OVERCLAIM` / `FLAG-SCOPE` / `FLAG-VOICE` / `FLAG-XSEC` (multi-flag 가능)
- 위반 항목별 line-anchored 노트 ("S2: 'X' phrase — suggested factual replacement value: Y / source: Z")
- Claude는 paragraph 자체 rewrite 금지. **suggested factual delta** + **flag** 만 제공.
- PASS 시: bib entry 검증 명령 (`grep -c LandaJCI2016 references.bib`) 추가 권고.

---

## Voice-protected reminder

- Hook ¶1 = 본인 키보드 strict. Default Claude 권한 = scaffolding/infra만 (`v17_sprint_vs_marathon_violation` 메모리).
- "고고" / "다 해줘" / "faster" 명령은 voice-protected sprint generate 권한 NOT 부여.
- Hook 작성 후 paste 시 본 문서 §[5] audit checklist 자동 실행 요청 가능.
- Audit 결과는 flag + suggested factual delta만; Claude가 Hook을 rewrite하지 않음.

---

Ready for user-written Hook paragraph.
