# v17 Methods Detailed


> ⚠️ **v5.2 retraction notice (2026-04-25):** some DIA-AUC and DIAL numbers below were computed under v5.1 leaky-ComBat protocol (full-pooled-data ComBat before LODO split, information leak). Under proper per-fold ComBat, the numbers shift. Current submission-ready figures are at `reports/html/pages/v17_npj_robustness.html` (manuscript v3 scenario-B reframe). See `reports/v5p2/v5p2_critical_assessment.md`.


## 1. 전체 설계 원리

v17의 목적은 갑상선암을 `BRAF-like vs RAS-like` 라는 2분법으로만 요약하지 않고, 그 바깥에 남는 표본을 다시 구조적으로 설명하는 데 있다. 데이터 흐름은 다음 다섯 단계로 나뉜다.

1. **Dark Matter 정의**: `TCGA-THCA` tumor 중 `driver_anchor ∉ {BRAF, RAS}` 인 표본을 추출한다.
2. **Dark Matter 내 subtype discovery**: TierA67 expression만으로 내부 구조를 찾는다.
3. **v17 driver map**: 기존 `driver_anchor` 위에 로컬 mutation summary와 fusion proxy를 덧씌워 확장 분류를 만든다.
4. **cross-cohort audit**: TCGA에서 학습한 Dark Matter subtype signature가 external cohort로 transfer되는지 본다.
5. **trajectory**: 전체 tumor를 저차원 manifold에 놓고 pseudotime 축에서 분화 상실과 고위험 쪽 이동을 본다.

핵심 원리는 단순하다. **driver classification**, **expression subtype**, **batch audit**, **progression geometry**를 별개로 보지 않고 한 표본 좌표계로 통합한다.

---

## 2. 입력 데이터와 표본 정의

기준 metadata는 `project/metadata/sample_master_v3.tsv` 이고, 전체 1,509개 sample을 포함한다. v17에서 직접 사용한 주요 subset은 다음과 같다.

- TCGA-THCA tumor: 513명
- TCGA Dark Matter: 178명
- trajectory usable set: 645 tumor samples

expression matrix는 각 코호트별 log2 scale matrix를 사용했다. `read_expr()`는 파일의 첫 열을 `gene_symbol`로 해석하고, 중복 gene symbol은 첫 occurrence만 유지한다. 이후 matrix를 gene-by-sample에서 sample-by-gene로 transpose했다.

---

## 3. Signature scoring

RAI score, TDS16 score는 모두 **단순 평균 signature score**로 계산했다.

표본 \(i\) 와 gene set \(G\) 에 대해 score는

\[
S_i(G) = \frac{1}{|G_i|} \sum_{g \in G_i} x_{ig}
\]

로 정의했다. 여기서

- \(x_{ig}\): sample \(i\) 의 gene \(g\) log-expression
- \(G_i\): 해당 sample matrix에 실제로 존재하는 gene subset

즉 누락 gene은 0으로 대체하지 않고 **관측된 gene만 평균**했다. 이건 platform 간 gene coverage 차이가 있는 external cohort에서 편향을 줄이기 위한 보수적 선택이다.

---

## 4. Task 1: Dark Matter consensus clustering

### 4.1 Feature space

Dark Matter clustering은 `TierA67_clean` gene panel 위에서 수행했다. 실제 matrix coverage는 코호트 gene 존재성에 의해 결정된다.

각 sample vector \(x_i \in \mathbb{R}^p\) 에 대해 먼저 `StandardScaler` 를 적용하여

\[
z_{ig} = \frac{x_{ig} - \mu_g}{\sigma_g}
\]

로 gene별 평균 0, 분산 1 스케일로 정규화했다.

그 다음 PCA를 적용해 상위 10개 주성분만 사용했다. PCA는 공분산 행렬

\[
\Sigma = \frac{1}{n-1} Z^\top Z
\]

의 고유분해로 얻은 loading vector를 사용하며, sample의 저차원 좌표는

\[
t_i = Z_i W_{1:10}
\]

이다.

이 단계의 이유는 두 가지다.

1. 54-55 gene panel에서도 gene 간 상관이 강하므로, PCA가 노이즈보다 공통 latent axis를 강조한다.
2. consensus clustering에서 반복 seed를 돌릴 때 원 feature보다 안정적이다.

### 4.2 Leiden / fallback clustering

사용 가능할 때는 `igraph + leidenalg` 로 k-nearest-neighbor graph를 만들고 Leiden partition을 사용했다. Leiden은 modularity 또는 CPM 계열 objective를 최적화하면서 disconnected community를 줄이는 graph clustering이다.

Leiden이 unavailable이면 `AgglomerativeClustering(linkage='ward')` 로 fallback 하도록 만들었다. Ward linkage는 매 병합 단계에서 cluster 내 제곱합 증가량

\[
\Delta(A,B) = \frac{|A||B|}{|A|+|B|}\|\bar{x}_A - \bar{x}_B\|^2
\]

이 최소인 pair를 합친다.

### 4.3 Consensus matrix

`K=2..10` 과 seed 30개를 돌렸다. 각 반복 \(r\) 에서 두 sample \(i,j\) 가 같은 cluster에 배정되면 \(M_{ij}^{(r)}=1\), 아니면 0으로 둔다. consensus matrix는

\[
C_{ij} = \frac{1}{R}\sum_{r=1}^{R} M_{ij}^{(r)}
\]

이다.

cluster stability는 같은 cluster 안 sample pair의 평균 consensus로 계산했다.

\[
\mathrm{Stability}(k) = \frac{1}{\sum_c |P_c|}\sum_c \sum_{(i,j)\in P_c} C_{ij}
\]

여기서 \(P_c\) 는 cluster \(c\) 내부 pair 집합이다. 최종 \(K\) 는 stability 최대값을 주는 값으로 선택했다.

이번 실행에서는

- best \(K=2\)
- stability \(= 0.9788\)

로 매우 높은 내부 일관성을 보였다.

### 4.4 Marker detection

각 cluster에 대해 one-vs-rest 방식으로 gene별 차이를 계산했다. gene \(g\) 에 대해

\[
t_g = \frac{\bar{x}_{g,1} - \bar{x}_{g,0}}
{\sqrt{s_{g,1}^2/n_1 + s_{g,0}^2/n_0}}
\]

형태의 Welch-type two-sample \(t\)-statistic을 사용했다. 이는 두 그룹 분산이 동일하다는 가정을 두지 않는다는 점에서 보수적이다.

p-value는 Benjamini-Hochberg FDR로 보정했다. \(m\)개의 p-value를 오름차순 \(p_{(1)} \le \cdots \le p_{(m)}\) 로 두면,

\[
q_{(i)} = \min_{j \ge i} \frac{m}{j}p_{(j)}
\]

로 계산한다.

log2 fold-change는 실질적으로

\[
\Delta_g = \bar{x}_{g,\mathrm{cluster}} - \bar{x}_{g,\mathrm{rest}}
\]

형태의 mean difference로 해석하면 된다. 원 matrix 자체가 이미 log2 scale이기 때문이다.

---

## 5. Task 2: Quad-group analysis

### 5.1 왜 로컬 MAF로 전환했는가

원래 설계는 cBioPortal `thca_tcga_pan_can_atlas_2018` mutation table을 직접 fetch해서 TERT promoter를 찾는 것이었다. 하지만 실행 환경에서 공개 archive fetch가 `403` 또는 LFS pointer 문제로 깨졌다. 그래서 로컬에 이미 내려받아 둔 GDC MAF로 대체했다.

### 5.2 로컬 mutation summary

각 MAF 파일에서 다음 gene만 골라 summary를 만들었다.

- `BRAF`, `HRAS`, `KRAS`, `NRAS`, `TERT`
- fusion 관련 gene: `RET`, `NTRK1`, `NTRK3`, `ALK`, `PAX8`, `PPARG`, `THADA`
- extra driver: `DICER1`, `EIF1AX`, `PPM1D`, `EZH1`, `EZH2`, `ATM`, `CHEK2`, `TP53`, `PIK3CA`, `AKT1`

추가로 세 개의 binary flag를 만들었다.

1. `braf_v600e`
2. `ras_hotspot`
3. `tert_promoter`

여기서 `ras_hotspot` 은 단백질 annotation에 `Q61`, `G12`, `G13` 이 들어가는지를 기준으로 잡았고, `tert_promoter` 는 canonical promoter hotspot 위치 또는 promoter-like classification을 찾도록 구현했다.

### 5.3 중요한 한계: TERT promoter 비관측

실행 결과 `tert_promoter` 는 0건이었다. 이건 THCA biology가 TERT promoter mutation이 없다는 뜻이 아니라, **현재 로컬 WXS MAF가 promoter hotspot을 포괄하지 못한다는 뜻**으로 해석해야 한다. 따라서 이번 v17의 TERT axis는 완전한 negative result가 아니라 **missing modality problem**이다.

### 5.4 임상 비교

임상 공변량은 `tcga_thca_clinical_extended.tsv` 에서 가져왔다.

- `stage`
- `age_at_diagnosis`
- `gender`
- `os_days`
- `os_event`

group별 평균 비교는 우선 ANOVA proxy를 사용했고, 사후비교는 Tukey HSD를 저장했다. 엄밀한 Welch ANOVA formula를 직접 구현하지는 않았고, 이번 스프린트 목적상 exploratory 비교로 유지했다.

### 5.5 Kaplan-Meier

생존함수는 Kaplan-Meier estimator를 직접 계산했다.

event time \(t_j\) 마다 at-risk 수를 \(n_j\), event 수를 \(d_j\) 라고 하면,

\[
\hat{S}(t) = \prod_{t_j \le t}\left(1 - \frac{d_j}{n_j}\right)
\]

이다.

group 간 비교는 log-rank test를 썼다. 각 event time마다 기대 event 수와 분산을 계산해

\[
\chi^2 = \frac{(O_1 - E_1)^2}{V_1}
\]

형태의 통계량을 사용했다.

이번 실행에서는 TERT group이 비어 있으므로 OS KM은 사실상 `BRAF only`, `RAS only`, `triple-negative` 3군 비교다.

---

## 6. Task 3: v17 driver landscape

### 6.1 우선순위 규칙

driver assignment는 우선순위 기반 deterministic rule로 구현했다.

1. `BRAF V600E`
2. `RAS hotspot`
3. `TERT`
4. `RET fusion`
5. `NTRK fusion`
6. `ALK fusion`
7. `PAX8PPARG`
8. `DICER1/EIF1AX/PPM1D`
9. `TP53`
10. fallback = existing `driver_anchor`

수학적으로는 multi-label mutation 상태 \(m_i\) 를 입력으로 받아, 위 우선순위 predicate를 차례대로 통과하는 **ordered decision function**

\[
f(i) = \arg\min_r \{ r : \phi_r(m_i)=1 \}
\]

로 볼 수 있다.

### 6.2 fusion proxy

이번 실행에서는 직접적인 TCGA fusion flat file 대신 `v3_fusion_anchor_tcga.tsv` 를 이용했다. 따라서 RET/NTRK/ALK/PAX8-PPARG 분류는 **완전한 raw fusion caller 기반 재산출이 아니라, 기존 internal anchor의 proxy mapping** 이다.

이는 figure를 빠르게 복구하기 위한 pragmatic choice였고, methods section에는 반드시 이 한계를 명시해야 한다.

### 6.3 confusion matrix

`molecular_subtype` (v14 라벨) 과 `driver_anchor_v17` 의 교차표를 만들어,

\[
N_{ab} = \#\{ i : \text{v14}_i=a,\ \text{v17}_i=b \}
\]

를 heatmap으로 그렸다. 이 confusion matrix는 “expression-defined subtype” 과 “driver-defined taxonomy” 사이의 불일치를 시각화한다.

---

## 7. Task 4: DIAL cross-cohort audit

### 7.1 철학

이 단계는 “DM1/DM2가 TCGA 안에서만 보이는 artifact인가?” 를 묻는 단계다. 하지만 external cohort에는 정답 cluster label이 없다. 그래서 **pseudo-label transfer** 구조를 썼다.

### 7.2 classifier

각 cluster에 대해 top marker 20개를 signature로 잡고, TCGA에서는 `cluster vs rest` logistic regression classifier를 적합했다.

로지스틱 회귀는

\[
\Pr(Y_i=1\mid x_i) = \sigma(\beta_0 + x_i^\top \beta)
\]

이며,

\[
\sigma(z) = \frac{1}{1+e^{-z}}
\]

이다. 구현은 `StandardScaler + LogisticRegression(max_iter=2000)` pipeline이다.

### 7.3 pseudo-label construction

external sample \(x\) 에 대해 TCGA positive centroid \(\mu_1\), negative centroid \(\mu_0\) 를 만들고,

\[
d_1(x) = \|x - \mu_1\|^2,\quad d_0(x)=\|x-\mu_0\|^2
\]

를 계산한다. \(d_1 < d_0\) 이면 pseudo-positive 로 본다.

이 pseudo-label을 기준으로 external predicted probability의 ROC AUC를 계산했다. true biological truth는 아니지만, **TCGA-defined geometry가 외부 코호트에서도 유지되는지** 보는 transfer proxy다.

### 7.4 direction-invariant AUC

direction-invariant AUC는

\[
\mathrm{DIA\text{-}AUC} = \max(\mathrm{AUC}, 1-\mathrm{AUC})
\]

로 계산했다. label 방향이 뒤집혀도 separation strength 자체는 유지되므로, transfer robustness를 볼 때는 sign보다 magnitude가 중요하다는 철학이다.

### 7.5 identifiability

cohort identifiability는 TCGA vs EXT batch label을 logistic regression으로 다시 맞출 수 있는지를 보는 값이다. 이 값이 1에 가까우면 cohort separation이 명확하고, 0.5에 가까우면 cohort indistinguishable 하다.

즉 identifiability는

\[
\mathrm{Identifiability} = \mathrm{AUC}\big(\hat{p}(B=\mathrm{EXT}\mid x), B\big)
\]

로 볼 수 있다.

이번 결과는 매우 인상적이다.

- transfer DIA-AUC는 0.969-1.000으로 높다.
- 그러나 ComBat 적용 후 identifiability는 1.000/0.9998에서 0.2577/0.1573으로 크게 떨어진다.

이건 “분리되는 biology” 와 “사라지는 cohort signature” 사이 긴장을 한 그림 안에 보여준다.

---

## 8. Task 5: Trajectory

### 8.1 공통 gene 교집합

여러 코호트 expression matrix를 단순 세로 결합하면 gene coverage 차이 때문에 NaN이 생긴다. 그래서 trajectory에서는 각 cohort에서 실제로 존재하는 TierA67 gene의 **교집합 50개**만 남겼다.

이는 batch correction과 manifold learning에서 필수다. 결측 유전자를 0으로 채우면 platform artifact가 pseudotime보다 더 강하게 작동할 수 있다.

### 8.2 ComBat

trajectory는 먼저 cohort batch를 보정하기 위해 `pycombat_norm` 을 시도했다. ComBat은 gene \(g\), sample \(i\), batch \(b(i)\) 에 대해 대략

\[
x_{ig} = \alpha_g + \mathbf{c}_i^\top \beta_g + \gamma_{g,b(i)} + \delta_{g,b(i)}\varepsilon_{ig}
\]

형태의 additive + multiplicative batch effect model을 가정한다.

여기서

- \(\alpha_g\): gene baseline
- \(\mathbf{c}_i\): biological covariate
- \(\gamma_{g,b}\): batch-specific location shift
- \(\delta_{g,b}\): batch-specific scale

를 empirical Bayes 방식으로 추정해 교정한다.

이번 trajectory 실행에서는 3-cohort usable subset에서 ComBat이 성공했다.

### 8.3 Manifold learning

Scanpy로 다음 순서를 탔다.

1. scaling
2. PCA
3. kNN graph
4. UMAP
5. diffusion map
6. DPT (diffusion pseudotime)

Diffusion map은 local transition graph를 Markov process로 보고, 확산 연산자의 leading eigenvector들로 sample manifold를 매개한다. DPT는 선택한 root cell에서 다른 sample까지의 diffusion distance를 이용한 진행축이다.

수학적으로는 transition matrix \(P\) 의 spectral decomposition 기반 저차원 embedding에서,

\[
\Psi_t(x_i) = \left(\lambda_1^t \psi_1(i), \lambda_2^t \psi_2(i), \dots\right)
\]

를 사용한다고 이해하면 된다.

root는 TDS16 score가 가장 높은 sample로 잡았다. 즉 **가장 thyroid-differentiated 상태를 시작점** 으로 가정하고, 그로부터 멀어지는 방향을 dedifferentiation / progression 쪽으로 해석했다.

### 8.4 dynamic genes

각 gene \(g\) 에 대해 pseudotime \(t_i\) 와 expression \(x_{ig}\) 의 Spearman correlation \(\rho_g\) 를 계산했다.

\[
\rho_g = \mathrm{corr}_{\mathrm{rank}}(x_{\cdot g}, t)
\]

Spearman을 쓴 이유는 pseudotime relation이 반드시 선형일 필요가 없고, monotone trend만 봐도 생물학 해석력이 충분하기 때문이다.

상위 음의 상관 유전자 `TPO`, `SLC26A4`, `DIO1`, `TG`, `DIO2`, `PAX8` 는 classic thyroid differentiation axis다. 반대로 양의 상관 `DUSP5`, `MET`, `LOX` 는 MAPK feedback/invasive remodeling 쪽 신호로 해석 가능하다.

---

## 9. 재현성과 난수

- Python: `project/.venv/bin/python`
- random seed: `42`
- sklearn random_state: `42`
- UMAP / CV split / PCA-related stochastic step도 가능하면 `42` 로 고정

이렇게 고정한 이유는 v17이 “정확한 cluster ID 재현”보다 “핵심 서사의 안정성”을 더 중요하게 보기 때문이다. 즉 seed가 달라도 **K=2 안정성**, **DIA-AUC 고값**, **differentiation gene의 pseudotime 음의 상관**이 유지되어야 한다.

---

## 10. 해석상 주의점

1. **TERT promoter 부재는 biological zero가 아니라 assay zero일 수 있다.**
2. **fusion은 raw caller 기반이 아니라 existing proxy 기반이다.**
3. **trajectory는 현재 5-cohort full이 아니라 3-cohort usable subset 결과다.**
4. **DIAL audit의 external pseudo-label은 진짜 ground truth가 아니라 centroid-based transfer proxy다.**

그럼에도 이번 v17의 강점은 분명하다. **Dark Matter를 “남는 찌꺼기”가 아니라 구조화된 subtype space로 승격**했고, 그것이 external transfer와 progression geometry 안에서도 일관된 패턴을 보였다는 점이다.
