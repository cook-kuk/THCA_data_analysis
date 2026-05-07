# Decisions Pending — Seungho 결정 필요 항목 일람

**작성일:** 2026-05-07
**Marathon deadline:** 2026-06-13 (Paper 1 bioRxiv 제출)
**Source:** Paper 1 submission-readiness audit + 2026-05-07 hub deep reconciliation 후 잔여
**모든 voice-protected는 본인 키보드만** (per `v17_sprint_vs_marathon_violation` 메모리 룰)

---

## 🔴 P0 — must-decide before bioRxiv (paper-blocking)

### D-A · Voice-protected 5개 작성 (Seungho 키보드만)
**Marathon mode rule:** "고고/faster/다 해줘"는 voice-protected sprint generate 권한 NOT.

| # | 파일 | 라인 | spec |
|---|---|---|---|
| **A1** Hook (≤45 words) | `project/manuscript_v8/04_intro_1_1_hook.md` + `03_introduction.md:15` | placeholder | 3 fact anchors lock(>98%/~20%/~23%) · "dark matter"/"DM1"/"76.8%" preview 금지 · 두 seed 옵션 (Alt B refined / Hybrid "compass") `02_outline.md:87-91`에 있음 |
| **A2** Discussion §3.1 opening (~200 words) | `06_discussion.md:13` | placeholder | 3-layer pathology framing + Landa 2016 dial-back convergence ("identifies upstream signature consistent with") · 본인 voice |
| **A3** Limitations §3.4 (250-350 words / 4 paragraphs) | `06_discussion.md:43` | empty | **8 specific items 필수**: ① TCGA OS 16 events ② MSK selection bias ③ Korean survival absent ④ methylation external absent ⑤ RAI outcome absent (panel name vs data gap) ⑥ within-DM1 fusion-independence n=19 underpowered ⑦ EA fusion validation absent ⑧ sub-B mechanism deferred |
| **A4** Cover letter paragraph 1 (~120-150 words / 3 sentences) | `08_cover_letter.md:19` | placeholder | S1 unmet need (BRAF/RAS-neg ~23% / ATA gap) · S2 contribution (DM1 76.8% fusion + TPO d=2.30 + HR 2.53) · S3 venue fit (Cell Rep Med = mechanism + multi-cohort + transl) |
| **A5** Reviewer Q9 (mechanism of silencing) | `09_reviewer_qa.md:48` | placeholder | "promoter methylation correlates without proving causation" + 후보 upstream (DNMT, SETDB1 등 본인 데이터 따라) + hold the line on "suggested not proven" |

**추정 시간:** 2-4시간 focused author time.

---

### D-B · 결정 단답 (factual, 입력 후 Claude 컴파일 가능)

| # | 결정 항목 | 옵션 | 추천 | 파일 |
|---|---|---|---|---|
| **B1** | Affiliations (Seungho + Yu Hyeong-won) | 본인 직위/소속 적기 | "TBD" 자리 fill | `01_abstract.md` author block + `08_cover_letter.md:42-48` |
| **B2** | 교신저자 email (Yu Hyeong-won) | 사용 email 입력 | "TBD" 자리 fill | `08_cover_letter.md:48` |

**추정 시간:** 15분.

---

## 🟡 P1 — high-value before bioRxiv (substantive but not paper-blocking)

### D-C · Aim sentence split (voice-adjacent)

| # | 결정 | spec |
|---|---|---|
| **C1** | `03_introduction.md:37` 145-word single sentence를 2개로 split | Sentence 1: "We applied an 8-gene RAI-responsiveness panel — independently selected from canonical thyroid differentiation biology (Yoo et al., 2016) before access to the Landa 2016 ATC-silenced gene list — across TCGA-THCA (n=504), MSK-IMPACT thyroid (n=117), Korean cohorts (n=865), and external single-cell datasets." (95w)<br/>Sentence 2: "We show that this panel resolves the BRAF/RAS-negative compartment into a DM1/DM2 axis that captures most tyrosine-kinase-fusion-positive tumors, harbors fusion-independent promoter hypermethylation of thyroid differentiation genes, and supports a clinically interpretable framework for reflex fusion testing and prospective evaluation of epigenetic-targeted RAI re-induction." (~50w) |

**추정 시간:** 15분 (구조는 정해져 있으므로 voice tweak만).

### D-D · Figure 정리

| # | 결정 | 추천 | 파일 |
|---|---|---|---|
| **D1 (P1-5)** | Fig 8C (within-DM1 fusion+ vs −, n=19, NS) main → supp 이동 여부 | 권장: 이동. small-n NS test as main panel은 reviewer 공격 표적 | `05_figure_captions.md` + `12_fig8_epigenetic.py` |
| **D2 (P1-7)** | Fig 7D (sub-A vs sub-B teaser) main → supp 이동해 main figure 8→7 | 권장: 이동 (Cell Press 표준 6-7) | `05_figure_captions.md` |

**Note:** Fig 8C/Fig 7D 이동을 결정하면 Claude가 즉시 figure caption 재배치 + supp 번호 부여 가능.

---

## 🟢 P2 — polish (post-bioRxiv OK)

| # | 항목 | 추천 |
|---|---|---|
| **E1** | STAR Methods cohort 수치 cross-check vs abstract | Claude 자동 가능 |
| **E2** | Repository URL + DOI 보관 language | bioRxiv DOI 받은 후 채움 |
| **E3** | 5명 → 6-7명 Suggested reviewers 확장 검토 (East Asian 1-2명 추가) | optional |
| **E4** | ATA 2025 (Ringel 2025) cite 검증 (PMID 40844370) | bib note 이미 있음 |

---

## 🟣 Strategic / cross-paper

### F · Paper 1 vs Paper 5 동시 submission 회피 (이미 진행중)

| # | 결정 | 현 상태 | 추천 |
|---|---|---|---|
| **F1** | Paper 5 (npj 8-gene RAI biomarker) reactivation 시점 | Ship-ready 2026-04-27, paused. Manuscript v6 + cover letter v6 + 49 figures + 12 supp tables 완료 | Paper 1 (Cell Rep Med) submission 후 결과 보고 결정. Paper 1 reach 실패시 Paper 5 npj fallback. 동시 submission ethics 위반. |
| **F2** | Paper 5 → npj submission portal click 누가 | "user sends + clicks submit" per `v17_npj_ship_status` | Seungho 본인 click. 4 outreach drafts 미발송 상태 — 별도 trigger. |

### G · Hub framing reconciliation 결과 검증

| # | 결정 | 적용 후 review |
|---|---|---|
| **G1** | paper1.html hub deep reconciliation (manuscript v8 framing 정렬) — **이번 세션 적용 완료**. Hub의 §3-§8 historical record 처리 | http://40.82.129.113:8012/papers_hub_2026_05_04/paper1.html 에서 직접 확인 후 ✓/× |
| **G2** | §0 driver framing fix (pool 포함, panel input 분리, leak-control comparator) — **적용 완료** | manuscript Q1과 일치 검증 → ✓ |
| **G3** | §8 TROP2 → "Historical note" 재라벨 — **적용 완료** | 페이지 review 후 ✓/× |

### H · GSE286332 drop 확인

| # | 결정 | 현 상태 |
|---|---|---|
| **H1** | Korean cohorts n=874 → n=865 (K2 235 + Lee 630, GSE286332-PTC(9) drop) — **이번 세션 적용 완료** | 4 파일 (`01_abstract.md`, `04_results.md`, `07_star_methods.md`, `02_outline.md`) 정합. Paper 2 boundary 보호. **확인 필요: 본인이 이 결정에 동의하는가?** 만약 keep 원하면 4 파일 revert. |

---

## 🟠 Operational / git / infrastructure

### I · Git state cleanup

| # | 결정 | 옵션 | 추천 |
|---|---|---|---|
| **I1** | 1782 pre-existing deleted tracked files (rsync 전부터 missing) | (a) `git rm --cached <files>` 일괄 정리 (b) 일부 복원 일부 cached 제거 (c) 그대로 두고 다음 commit에서 자연스럽게 | **(a)** — 어차피 disk에 없는 파일들. Commit message에 "session 시작 전 stale tracked deletion 일괄 정리" 명시. |
| **I2** | 145 untracked files 처리 | (i) Track B-lite 결과 commit + bak/.md gitignore (ii) Marathon 끝까지 보류 | **(i)** 부분 적용 — 명확한 artifact는 commit, bak 파일 삭제 |
| **I3** | 이번 세션 변경사항 commit 시점 | 즉시 / voice-protected 들어간 후 / bioRxiv 직전 | **voice-protected 5개 들어간 후** 단일 commit 권장: "session 2026-05-07 disk migration + 12-paper section 0 + CLAUDE.md + manuscript v8 hub reconciliation + voice-protected" |

### J · Pod / data infrastructure (이미 처리됨)

| # | 결정 | 현 상태 |
|---|---|---|
| **J1** | RunPod pod 사용 시점 | EXITED $0/hr 또는 새 A6000 pod 사용 — Marathon writing은 GPU 불필요. 6/13 ship 후 Track B 시작 시 재셋업. |
| **J2** | arcasHLA 32G `_repo_offload` 데이터 운명 | 그대로 보존 (디스크 358G 여유). 6/13 ship 후 정리 검토. |

---

## 🔵 결정 후 즉시 Claude 가능 (mechanical)

본인 D-A (voice-protected 5개) + D-B (affiliations) + D-C (Aim split) + D-D (Fig 정리) 결정 후 Claude는 다음 작업 즉시 가능:

1. `10_full_manuscript_compiled.md` regen
2. PDF build (기존 pipeline)
3. bioRxiv supplementary tables 최종 정렬
4. Cover letter formatting (date, 발신/수신 채움)
5. Reviewer Q&A 12개 final consistency check
6. STAR Methods cohort 수치 cross-check
7. (선택) git commit 작성
8. (선택) Hub paper1.html에 manuscript-final-aligned 마지막 polish

---

## 🎯 Quick-action template (체크박스)

```
P0 — bioRxiv-blocking
[ ] A1 Hook (~45w)             — 04_intro_1_1_hook.md
[ ] A2 Discussion 3.1 opening   — 06_discussion.md:13
[ ] A3 Limitations 3.4 (8 items) — 06_discussion.md:43
[ ] A4 Cover letter para 1     — 08_cover_letter.md:19
[ ] A5 Reviewer Q9             — 09_reviewer_qa.md:48
[ ] B1 Affiliations input      — 01_abstract.md + 08_cover_letter.md
[ ] B2 Corresponding email     — 08_cover_letter.md:48

P1 — high-value
[ ] C1 Aim §1.4 sentence split — 03_introduction.md:37
[ ] D1 Fig 8C demote decision  — main → supp 여부
[ ] D2 Fig 7D demote decision  — main 8→7 여부

Strategic verification (already applied — verify only)
[ ] G1 paper1.html hub review (http://40.82.129.113:8012/...) — OK / 추가 변경 요청
[ ] H1 GSE286332 drop OK 확인 — n=865 keep / revert to n=874

Operational
[ ] I1 Git stale 1782 cleanup 시점 — 즉시 / 후순위
[ ] I2 Untracked 145 처리 정책   — commit / gitignore / 보류
[ ] I3 본 세션 commit timing     — voice-protected 후 / 즉시 / 6/13 직전

Post-decision (Claude 자동)
[ ] Compile full_manuscript_compiled.md regen
[ ] PDF build
[ ] bioRxiv submission preparation
```

---

## 🕒 추정 잔여 시간 (6/13 deadline까지)

| Owner | Item | 시간 |
|---|---|---|
| Seungho | A1-A5 voice-protected 5개 | 2-4시간 |
| Seungho | B1-B2 affiliations | 15분 |
| Seungho | C1 Aim split | 15분 |
| Seungho | D1-D2 figure demote 결정 | 15분 |
| Seungho | G1/H1 hub + GSE286332 verify | 30분 |
| Claude | Compile + PDF + supp polish | 1시간 |
| Seungho | bioRxiv upload click | 30분 |

**총: Seungho 약 4-6시간 + Claude 1시간 = 5-7시간 잔여 작업.**
**6/13 deadline까지 5주+ 여유** — 매우 충분.

---

## 📚 참조

- 원 audit: `/home/seungho/personal/THCA_data_analysis/SESSION_FULL_2026_05_07_DISK_POD_HUB.md` D1-D8 항목
- Paper 1 submission audit: 본 세션 대화 — "Paper 1 Submission-Readiness Audit (2026-05-07)" + "Edit Plan §B"
- 메모리: `v17_marathon_mode_post_pillar1` · `v17_sprint_vs_marathon_violation` · `v17_npj_ship_status` · `v18_paper2_HT_isolated`
- Hub: http://40.82.129.113:8012/papers_hub_2026_05_04/paper1.html
- CLAUDE.md: `/home/seungho/personal/THCA_data_analysis/CLAUDE.md`

---

*Generated 2026-05-07 by Claude Opus 4.7 (1M) — Paper 1 audit consolidation.*
*All voice-protected items remain author keyboard.*
