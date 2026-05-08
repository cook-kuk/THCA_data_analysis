# HF_TOKEN — UNI/CONCH foundation model 액세스 가이드

**Status:** UNI gated; v2 sprint 현재 ImageNet ViT-L fallback 으로 진행 중. UNI 또는 CONCH 액세스 받으면 *진정한 foundation model* 결과로 Phase 1+2 재실행 가능.

## 1. UNI 액세스 (5분, 학술용 즉시 승인)

1. https://huggingface.co/MahmoodLab/UNI 방문
2. "Access this model" 클릭 → academic email 입력 + agree to license
3. 보통 즉시 ~1시간 내 승인 (수동 review)
4. https://huggingface.co/settings/tokens → "New token" → "read" scope → 복사
5. Pod 에 export:
   ```bash
   export HF_TOKEN=hf_xxx_your_token
   ```
6. Phase 1 재실행:
   ```bash
   cd /workspace/project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/
   bash runpod_dispatch/run_phase1.sh
   ```

## 2. CONCH 액세스 (alternative)

- https://huggingface.co/MahmoodLab/CONCH
- 동일 절차

## 3. ViT-L fallback (현재 상태) 의 실제 의미

- **closure baseline 은 ResNet50 (ImageNet)** — 전혀 못 잡음 (~0.55)
- **현재 v2 fallback 은 ViT-L (ImageNet)** — 더 큰 / 더 강한 generic vision encoder
  - ImageNet 만 학습; pathology 도메인 모름
  - 그러나 generic visual texture / morphology 는 캡처
  - Closure 의 ResNet50 보다 약간 좋을 가능성 있음
- **UNI 는 Mass-100K pathology WSI 학습** — pathology 도메인 specific
  - generic ImageNet 훈련 ViT-L 와 *완전히 다른 generalization*
- **결론:** ViT-L fallback 으로 *0.55 < AUC < 0.70* 나오면 → "v2 ViT-L 가 closure 를 부분 능가" 를 paper 에 적을 수 있음
  - 그래도 UNI 를 받으면 더 강한 결과 가능 — *future paper / supplement* 에 추가
  - CONCH (multimodal) 는 다른 특성 → exploratory 용도

## 4. 만약 UNI 액세스 안 받으면

- 본 sprint 은 ViT-L baseline 으로 종결
- 결과 framing: "with ImageNet ViT-L (foundation-model not used due to gating); UNI / CONCH 사용시 더 강한 결과 가능"
- closure 와의 comparison 가능성: closure ResNet50 AUC ~0.55 vs v2 ViT-L AUC = X
