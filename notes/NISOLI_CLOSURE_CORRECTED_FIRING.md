# NISOLI_CLOSURE_CORRECTED_FIRING — firing cells of the corrected-rate Nisoli closure

**Date:** 2026-05-12. Phase 3 of Route A.

---

## 1. Firing cells, under `C_A = 1` parameterisation

The closure inequality `|K_bil(r)| · K^{−A} · M_3'' < 1` (M_3'' = 24.43, p = 3) fires at the following minimum-(A, K) per r:

### r = 3 (strict bilinear 2√N, N=9, |K_bil| = 6, threshold 146.6):

| A | K_min | product at K_min | smallest plausibly-feasible cell |
|---|---|---|---|
| 1 | 147 | ~0.998 | (r=3, A=1, K=147) — fires marginally |
| 2 | 13 | 0.867 | **(r=3, A=2, K=13)** ← smallest firing cell at r=3 |
| 3 | 6 | 0.679 | **(r=3, A=3, K=6)** ← smallest K firing cell |
| 5 | 3 | 0.602 | (r=3, A=5, K=3) ← if A is large, K can be tiny |
| 10 | 2 | 0.143 | (r=3, A=10, K=2) |

### r = 4 (polylog-free 2√3·√N, N=27, |K_bil| = 18, threshold 439.7):

| A | K_min | product | smallest cell |
|---|---|---|---|
| 1 | 440 | ~0.999 | (r=4, A=1, K=440) |
| 2 | 21 | 0.998 | **(r=4, A=2, K=21)** |
| 3 | 8 | 0.859 | **(r=4, A=3, K=8)** |
| 5 | 4 | 0.430 | (r=4, A=5, K=4) |

### r = 5 (polylog-free, N=81, |K_bil| ≈ 31.18, threshold 761.7):

| A | K_min | product | smallest cell |
|---|---|---|---|
| 2 | 28 | 0.972 | **(r=5, A=2, K=28)** |
| 3 | 10 | 0.762 | **(r=5, A=3, K=10)** |
| 5 | 4 | 0.736 | (r=5, A=5, K=4) |

### r = 6 (polylog-free, N=243, |K_bil| ≈ 54, threshold 1319):

| A | K_min | product | smallest cell |
|---|---|---|---|
| 2 | 37 | 0.964 | (r=6, A=2, K=37) |
| 3 | 11 | 0.992 | (r=6, A=3, K=11) |
| 5 | 5 | 0.378 | (r=6, A=5, K=5) |

---

## 2. The smallest absolute firing cell

Across all (r, A, K) with A ≤ 10 and K ≤ 20:

> **(r=3, A=3, K=6): product = |K_bil(3)| · M_3'' · K^{−A} = 6 · 24.43 · 6^{−3} = 146.6/216 = 0.679 < 1** ✓

This is a **single-digit-A, single-digit-K cell where the closure inequality fires** under the optimistic `C_A = 1` parameterisation.

Comparable in r=4: **(r=4, A=3, K=8)** with product 0.859. r=5: **(r=5, A=3, K=10)** with product 0.762.

These are mathematically the firing cells if Tao Prop 1.17 delivered the `C_A = 1` algebraic decay.

---

## 3. The actual Tao Prop 1.17 firing question

Tao Prop 1.17's effective form (with C_A absorbed):

> **`|K_bil(r)| · C_A · n^{−A} · M_3'' < 1`** (where K and n are conflated — Tao's "n" is the level of the Syracuse modulus 3^n; here we identify n = K for closure parameterisation)

Required: `n^A > |K_bil(r)| · C_A · M_3''`.

Per BOOKKEEPING_PHASE1_DISPOSITION:
- C_A grows at least like `A^{O(A)}` (optimistic floor) up to `exp(exp(A²))` (faithful bookkeeping).
- At A = 3, K = 6: optimistic floor `C_3 ≥ 3^{c·3} ≈ 3^3 = 27` to `3^9 = 19683`. Faithful bookkeeping ≥ `exp(exp(9))` ≈ astronomical.

Substituting `C_A = A^A` (very optimistic):

| r | A | K | K^A | required (|K_bil|·C_A·M_3'') |
|---|---|---|---|---|
| 3 | 3 | 6 | 216 | 6 · 3³ · 24.43 = 6 · 27 · 24.43 = 3,958 | NO (216 < 3,958) |
| 3 | 2 | 13 | 169 | 6 · 2² · 24.43 = 6 · 4 · 24.43 = 587 | NO (169 < 587) |
| 3 | 5 | 3 | 243 | 6 · 5⁵ · 24.43 = 6 · 3125 · 24.43 = 458,062 | NO |
| 3 | 10 | 2 | 1024 | 6 · 10¹⁰ · 24.43 = 1.47 × 10¹² | NO |

**Under the optimistic-floor `C_A = A^A`, NO cell with reasonable (A, K, r) fires.**

Under faithful bookkeeping `C_A ≥ exp(exp(A²))`, the situation is dramatically worse — at A = 2, C_2 ≥ exp(exp(4)) ≈ exp(54.6) ≈ 5 × 10²³, fully eclipsing any K^A gain.

---

## 4. Conditional disposition

> **The Nisoli closure inequality fires at small-(A, K, r) under the optimistic `C_A = 1` parameterisation, but does NOT fire under any realistic Tao Prop 1.17 bookkeeping of `C_A`.**

This routes to:

- **H_NISOLI_CLOSES_AT_R3** (CONDITIONALLY, given `C_A = 1`): firing cell (r=3, A=3, K=6) — well-defined, sub-1 product, structurally the smallest closure. **NOT a rigorous closure** because the implicit `C_A = 1` is not what Tao delivers.
- **H_A_EXTRACTION_HARD** (UNCONDITIONALLY, given BOOKKEEPING_PHASE1's finding): A cannot be extracted from Tao Prop 1.17 in any form that makes the inequality fire with rigorous constants. The corrected-rate Nisoli closure is **not closable via this route**, mirroring the rate-1/2 case (R77.3 falsified + BOOKKEEPING_PHASE1 INFEASIBLE).

The disposition for the corrected-rate probe lands at **H_A_EXTRACTION_HARD**, not at H_NISOLI_CLOSES_AT_R3, because the Tao constant is the load-bearing obstruction.

---

## 5. What would unlock H_NISOLI_CLOSES_AT_R3 (rigorously)

Three independent unlock paths:

### (a) A new derivation of `|μ̂_n(ξ)| ≤ C · n^{−A}` with polynomial `C` in A

Outside of Tao's renewal-process method. This is the "different Fourier-analytic bound on the Syrac MGF with polynomial-in-A constant" alternative noted in BOOKKEEPING_PHASE1_DISPOSITION §"Alternative routes." **Status: not in the project; would require novel technique.**

### (b) The cross-freq machinery delivering A directly

Cross-freq derives Off_lin's weights `W_+(g) = 2^{−g+1}/15` for g ∈ {2, 4, 6, ...}. The W_+ weights are **geometric in g** — they decay like `(1/2)^g` × constant. This gives an `A_cross_freq = ∞` in some sense (geometric decay beats any polynomial). However:

- This is the decay of Off_lin's contributions to T_lead, not the decay of `|μ̂_n(ξ)|`.
- These are different quantities: cross-freq weights are matrix-entry decays of T_lead; Tao's A is the per-frequency Fourier decay of the Syracuse measure.

So cross-freq's geometric decay is **structurally a different object** from Tao's A — it can't substitute directly. **The cross-freq machinery does NOT directly deliver Tao's A.**

### (c) T_lead's eigenvalue structure gives A = ∞ semantically

T_lead has spectrum {43/45, 0}. The eigenvalue 0 (on (9, −7)) means certain modes are killed immediately by T_lead. One might argue this corresponds to "infinite A" for those modes. But:

- T_lead is a within-level operator; Nisoli's A applies to the asymptotic-in-n decay of `|μ̂_n(ξ)|`, not the within-level rank-1 structure of T_lead.
- The two A's measure different things.

So **T_lead's spectrum gives a positive structural finding (43/45) but does NOT replace Tao Prop 1.17's A in the Nisoli closure inequality.**

---

## 6. Summary

| Question | Answer |
|---|---|
| Does the closure inequality fire under `C_A = 1`? | **YES** at (r=3, A=3, K=6); product 0.679. |
| Does it fire under Tao's actual bookkeeping `C_A`? | **NO** at any plausible (A, K) under either optimistic-floor or faithful-bookkeeping `C_A`. |
| Can A be extracted from cross-freq, T_lead spectrum, or other project machinery? | **NO** — these are structurally different objects from Tao's A. |
| What would unlock? | A novel polynomial-in-A bound on the Syracuse Fourier coefficient, **outside Tao's method**. |

> **Phase 3 disposition: closure inequality is conditional-firing at (r=3, A=3, K=6) but unconditionally non-firing under realistic C_A. The Phase 4 honest read documents A-extractability.**

---

## 7. Files

- NISOLI_CLOSURE_CORRECTED_FIRING.md (this file)
- BOOKKEEPING_PHASE1_DISPOSITION.md (Tao C_A INFEASIBLE)
- NISOLI_CLOSURE_CORRECTED_HONEST_READ.md — Phase 4
