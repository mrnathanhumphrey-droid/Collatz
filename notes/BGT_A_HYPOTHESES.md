# BGT_A — Karamata's representation theorem

## Phase 0 — verbatim statement

**Source:** Kevei 2019 "Regularly Varying Functions: notes for a 2-hour reading course", Theorem 6 (= BGT 1.3.1; Hawkes 2RV paper Theorem 2.7). Lines 181-198 of `C:/tmp/bgt/Kevei_2019_Regularly_Varying_Functions_Notes.txt`.

> **Theorem 6 (Representation theorem).** Let ℓ be a nonnegative measurable function. It is slowly varying if and only if
> ℓ(x) = c(x) · exp{ ∫_a^x ε(u)/u du }, x > a,
> where a ≥ 0, lim_{x→∞} c(x) = c ∈ (0,∞), lim_{x→∞} ε(x) = 0.

Standing definition (Kevei Def 3, line 162):
> **Definition 3.** A nonnegative measurable function ℓ : [a,∞) → [0,∞), a ≥ 0, is slowly varying if
> lim_{x→∞} ℓ(λx) / ℓ(x) = 1 for each λ > 0.

## Hypothesis types

- h_1 (function class): ℓ is nonnegative real-valued (sign constraint).
- h_2 (operational): for each λ > 0, lim_{x→∞} ℓ(λx)/ℓ(x) = 1 (slow-variation condition).
- h_3 (input-side object): continuous function or sequence (sequence interpreted via Kendall-style passage to limits, see candidate E).

Conclusion: explicit additive-representation form ℓ(x) = c(x) · exp{∫ ε(u)/u du}.

---

## Phase 1 — hypothesis × ε_k matrix

The natural mapping: take ℓ(k) := |ε_k| · 2^k (the "normalized" sequence), interpreted as a sequence-valued slowly-varying candidate (treat the index k as the variable, then continuous version: ℓ(x) = |ε_{⌊x⌋}|·2^{⌊x⌋}, or via Kendall passage).

Test data:
- L(2)=0.038, L(3)=0.041, L(4)=0.039, L(5)=0.037, L(6)=0.032, L(7)=0.150, L(8)=0.191.

| hyp | check | verdict |
|---|---|---|
| h_1 | ℓ(k) ≥ 0 (taking absolute value) | SATISFIED (by construction). |
| h_2 (within plateau k=2..6) | L(k+1)/L(k) values 1.069, 0.963, 0.939, 0.865 → drifting, not converging to 1 from N=4 data | NEEDS_PROOF (suggestive of slow variation; range [0.86, 1.07] within plateau but trending downward). |
| h_2 (across jump k=6→7) | L(7)/L(6) = 4.72 — emphatically NOT close to 1 | **FAILED** for ratio λ ≈ 7/6. |
| h_2 (post-jump k=7→8) | L(8)/L(7) = 1.27 — also not close to 1 | FAILED. |
| h_3 | sequence — interpretable via Kendall (candidate E) | SATISFIED via E's machinery. |

**Phase 1 verdict: FAILED at h_2 (slow-variation across the k=7 jump).** L(k) is NOT slowly varying as a function of k under the standard λ-multiplicative slow-variation test. The k=6→7 ratio 4.72 categorically rules out h_2.

A weaker possibility: |ε_k| itself (without the 2^k normalization) might be regularly varying with index ρ = −1 (consistent with the plateau ratios |ε_{k+1}|/|ε_k| ≈ 0.43-0.53 within k=2..6, which is close to 1/2 = 2^{-1}, matching the predicted rate-1/2 envelope). But k=6→7 ratio is 2.36 — incompatible with any single index ρ.

Conclusion: **NO_FIT at Phase 1.**

---

## Phase 2 — conclusion shape

If h_2 held, the conclusion would give |ε_k|·2^k = c(k) · exp{ ∫ ε̃(u)/u du } with c(k) → c. This is a slowly-varying representation. Converted via R76 leading-mode identity + R77 T_diag eigenvector (1,4): the |μ̂_n(ξ)|² ≤ S_n / 2 bound via Plancherel, and S_n ~ 7/15 − (1/30)·(1/2)^n + O((1/4)^n). The representation form gives the slowly-varying *factor* on the prefactor only; the rate (1/2)^n itself comes from the regular-variation INDEX, not the representation. So this candidate's conclusion gives at best a prefactor refinement, not a new rate.

For polynomial-in-A bound on |μ̂_n(ξ)|: NOT delivered. The representation theorem delivers a structural form for slowly-varying functions, not a Fourier-decay estimate.

**Phase 2 verdict: SHAPE_MISMATCH** (structural representation, not Fourier-decay bound).

---

## Phase 3 — multi-regime check

L(k)=|ε_k|·2^k has a 4.72× jump at k=6→7. The slow-variation hypothesis h_2 is *defined* as a single-regime convergence; multi-regime structure violates the hypothesis directly.

The theorem's representation form admits a single c(·) → c limit; if L(k) has two regimes (plateau + post-jump), the c(x) would need to be sequence-discontinuous (which is allowed if it's just measurable, but then the limit-condition fails).

**Phase 3 verdict: STRUCTURALLY_BLOCKED.**

---

## Disposition: NO_FIT

Karamata representation theorem requires slow variation, which |ε_k|·2^k fails at the k=7 jump. Even if the plateau k=2..6 hints at slow variation, the across-jump ratio 4.72 categorically fails h_2.
