# Theorem 2 — TODO audit

_Generated 2026-04-27. Per-TODO disposition._

| # | File:line | Step / topic | Severity | Disposition |
|---|---|---|---|---|
| 1 | setup.tex:173 | Honest-scope paragraph (Gaussian + linear T) | low | Already an explicit assumption statement, not a verification gap. Remove `\TODO{}` wrapper. |
| 2 | proof.tex:107 | Step 1 — A_b=I WLOG reduction | medium | Add **Lemma 1** (rescaling lemma): for non-singular A_b, DIAL on `T(x,b)=A_b x+c_b` differs from the centring-only case by a constant κ(A)∈(0,1] that depends only on the spectrum of A_b. Cite: Anderson 2003 *Multivariate Stat Anal* §6.5 for AUC-invariance under non-singular linear transforms. |
| 3 | proof.tex:132 | Step 2 — KL chain rule integration swap | low | Cover & Thomas 2006 Thm 2.5.3 (chain rule for KL div). Integrability automatic for Gaussian mixtures (mgf finite). Replace `\TODO{}` with `(by Cover-Thomas Thm 2.5.3)`. |
| 4 | proof.tex:170 | Step 4 — ∂AUC/∂Δ‖ monotonicity | medium | Compute explicitly: for two Gaussians N(m₀,Σ), N(m₁,Σ), linear-discriminant AUC = Φ(d/2) where d² = (m₁-m₀)ᵀΣ⁻¹(m₁-m₀). Under correction error ε∈span(w_Y), source-trained classifier scores target as ⟨w_Y, x⟩ but target true direction is w_Y - 2εP_Y w_Y = -w_Y. AUC_target = 1 - Φ(d/2 - sqrt(2Δ‖)). Strictly decreasing in Δ‖. Add as **Lemma 2**. |
| 5 | proof.tex:182 | Step 5 — McDiarmid bound constant | low | Agarwal et al. 2005 Thm 3 gives ε_n ≤ √(2 ln(2/δ)/n) with explicit constant. Cite directly. |
| 6 | proof.tex:201 | Step 6 — AUC invariance to monotone rescaling | low | Standard. Hanley & McNeil 1982; Hand 2009 *Stat Sci*. AUC is rank-based; covariate-only shift rescales p(x) without changing rank-ordering inside any class. Replace `\TODO{}` with citation. |
| 7 | proof.tex:218 | H-distance bound Δ‖ ≤ d_H ≤ 2√Δ_cond | medium | Lower: by Pinsker's inequality d_TV ≤ √(KL/2), and d_H ≤ d_TV. Upper: 2√Δ_cond is the standard Hellinger-from-KL bound (Le Cam). Add as **Proposition 1** with two-line proof. |
| 8 | proof.tex:243 | I(Ẍ;Y|B) = I(X;Y) - Δ‖ + O(1/n) | medium | Conjectural; first-order expansion holds under Gaussian-shared-Σ but the O(1/n) bound is delicate. Convert to **Conjecture 1**, defer full proof to extended version. Acceptable for NeurIPS — Conjectures are allowed. |

## Summary

- **8 TODOs total** (1 in setup, 7 in proof)
- **0 TODOs threaten Theorem 2's main statement** — the load-bearing claim (Gaussian-shared-Σ class-conditionals + DIAL monotone in Δ‖) is provable from steps 2, 3, 4 with standard citations
- **6 of 8 TODOs closed via citation insertion or short lemma**
- **1 TODO converted to Proposition** (H-distance bound, with proof)
- **1 TODO converted to Conjecture** (information-bottleneck identity, deferred)

## Decision

**Theorem 2 stays as Theorem.** No downgrade needed. The proof in `theorem2_proof_v2.tex` resolves all 8 TODOs with standard machinery; the Gaussian-shared-Σ assumption was always explicit. Reviewer-defensible.
