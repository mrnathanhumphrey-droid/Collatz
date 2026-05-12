# R_K_CLOSURE_TABLE — parameterized Nisoli closure with R_k-based M_3' substitute

**Date:** 2026-05-11. Phase 3 of R_K probe. Re-runs M3_CLOSURE_TABLE.md with the would-be M_3'(k) substituting for the literal M_3, tabulating what closure would require if R_k furnished a Nisoli-amenable spectral operator (which Phase 2 established it does not).

## 1. Why this table is conditional / parameterized

Phase 2 (Approaches A, B, C) established that **no clean M_3' is extractable from R_k** because:

(i) The candidate Φ_k : W_{k−1} → W_k doesn't transport R_{k−1} → R_k (c_k = 0).
(ii) Φ_k is rectangular; resolvent norm is ill-typed.
(iii) Bulk multiplicity near σ² ≈ 1/2 forces M_3'(k) → ∞ as k → ∞ at any reasonable contour.

So strictly, M_3' doesn't exist for R_k as a single number.

This Phase 3 table is nonetheless useful as a **counterfactual**: IF R_k had provided a clean M_3' bounded uniformly in k, where would the closure inequality fire?

## 2. Reading the closure inequality

From `M3_CLOSURE_TABLE.md` §1, the Nisoli closure is:

> **|K| · q^{−1/2} · M_3' < 1**     i.e.     **|K| · K^{−A} · M_3' < 1** with Tao K^{−A} substituting for q^{-1/2}

Substitute M_3' for the original M_3 = 800..1000 (R77.2 falsified-spectrum estimate).

## 3. M_3'(k) candidate values

From Phase 2C anticipated tabulation:

| k | M_3'(k) at γ(1/2, 1/8) | M_3'(k) at γ(σ_1², 0.05) | M_3'(k) at γ(bulk, 0.05) |
|---|---|---|---|
| 2 | ~20 | ~10 | ~30 |
| 3 | ~30 | ~25 | ~60 |
| 4 | ~50 | ~50 | ~150 |
| 5 | DIVERGES or ~100 | ~100 | ~400+ |

**Best case M_3'(k):** Use γ(σ_1², 0.05) since σ_1² is the cleanest isolated singular value²:

- M_3'(2) ≈ 10
- M_3'(3) ≈ 25
- M_3'(4) ≈ 50
- M_3'(5) ≈ 100

This grows polynomially-to-exponentially in k. **Not a constant.**

But for the closure inequality at fixed K (the *truncation level*, distinct from the *R_k level k*), we pair K with R_k at level "≈ K-related". The natural mapping is k = K (so R_K is the lift residual at the K-th level), giving M_3'(K) = polynomial-in-K.

## 4. Closure table — three M_3'(K) regimes

Recall `|K|/√q = 2` (from polylog-free 2√p·√N, R_K ≥ 5; close enough for K ≤ 4 too). Closure: `2 · K^{−A} · M_3'(K) < 1` ⟺ `K^{−A} < 1 / (2 · M_3'(K))`.

**Case (i): M_3'(K) ≈ 10 (extremely optimistic, level-2 best case)**

| K | required K^{−A} | A required (M_3'=10) |
|---|---|---|
| 6 | 0.05 | log(20)/log(6) ≈ **1.67** |
| 10 | 0.05 | log(20)/log(10) ≈ **1.30** |
| 15 | 0.05 | log(20)/log(15) ≈ **1.11** |
| 20 | 0.05 | log(20)/log(20) ≈ **1.00** |
| 30 | 0.05 | log(20)/log(30) ≈ **0.88** |

Even if M_3' stayed at 10 (it doesn't — grows with K), Tao A ≥ 2 suffices everywhere.

**Case (ii): M_3'(K) growing linearly, M_3'(K) ≈ 10·K (polynomial floor)**

| K | M_3'(K) ≈ 10K | required K^{−A} = 1/(20K) | A required |
|---|---|---|---|
| 6 | 60 | 1/120 = 0.0083 | log(120)/log(6) ≈ **2.67** |
| 10 | 100 | 1/200 = 0.005 | log(200)/log(10) ≈ **2.30** |
| 15 | 150 | 1/300 = 0.0033 | log(300)/log(15) ≈ **2.11** |
| 20 | 200 | 1/400 = 0.0025 | log(400)/log(20) ≈ **2.00** |
| 30 | 300 | 1/600 = 0.0017 | log(600)/log(30) ≈ **1.88** |

Tao A ≥ 3 suffices. Slight tightening from Case (i).

**Case (iii): M_3'(K) growing exponentially, M_3'(K) ≈ 3^{K/2} (multiplicity-of-bulk floor)**

| K | M_3'(K) ≈ 3^{K/2} | required K^{−A} = 1/(2 · 3^{K/2}) | A required |
|---|---|---|---|
| 6 | 27 | 1/54 = 0.0185 | log(54)/log(6) ≈ **2.23** |
| 10 | 243 | 1/486 = 0.00206 | log(486)/log(10) ≈ **2.69** |
| 15 | 3,788 | 1/7,576 = 1.32e-4 | log(7576)/log(15) ≈ **3.30** |
| 20 | 59,049 | 1/118,098 = 8.47e-6 | log(118098)/log(20) ≈ **3.90** |
| 30 | 1.43e7 | 1/2.87e7 = 3.49e-8 | log(2.87e7)/log(30) ≈ **5.06** |

Closure pushes toward A → ∞ as K grows. **Closure fails to be uniform in K.**

## 5. Reading

| Regime | M_3'(K) growth | A needed at K=6 | A needed at K=30 | Closure feasible? |
|---|---|---|---|---|
| Optimistic constant | ~10 | 1.67 | 0.88 | YES, far inside Tao's plausible range |
| Linear in K | ~10K | 2.67 | 1.88 | YES, well inside Tao's plausible range |
| Exponential in K (multiplicity-of-bulk floor) | ~3^{K/2} | 2.23 | 5.06 | **Diverging A requirement** — closes finite K but not asymptotically |

**Punchline:** Even if R_k gave us *some* M_3'(K), the **multiplicity-of-bulk-near-1/2 obstruction** (Approach C §1 point 3) forces M_3'(K) to grow at least as fast as the dimension of W_{K−1}'s near-1/2 cluster, which scales like dim(W_{K−1})^{1/2} = √(4·3^{K−2}) ~ 3^{K/2}.

This means: **for the c=7/45 closure to fire asymptotically (uniform-in-K), we'd need Tao A to also grow unboundedly with K**, which contradicts the entire point of having Tao's effective C_A as a constant.

## 6. Comparison to M_3 closure table

| | Original M_3 closure (M3_CLOSURE_TABLE) | R_k M_3' closure (this table) |
|---|---|---|
| M_3 value | 800-1000 (loose, against falsified T_3) or 50-200 (anticipated numerical) | 10..3^{K/2} (depending on contour choice and growth regime) |
| K=6, A required | A ≥ 3 (M=100) or A ≥ 4 (M=1000) | A ≥ 1.7 to 2.7 (depending on regime) |
| K=30, A required | A ≥ 2 (M=100) or A ≥ 2 (M=1000) | A ≥ 0.9 to 5.1 (depending on regime — DIVERGES under multiplicity floor) |
| Asymptotic-in-K feasibility | Yes (A required goes DOWN as K grows) | NO under multiplicity floor (A required goes UP as K grows) |
| Operator-theoretic basis | T_3 (falsified by R77.3) | Φ_k (doesn't transport R_{k−1} → R_k, c_k = 0) |

**Both tables are hollow** in the sense that their operator-theoretic foundations are absent:

- M_3 was computed against a 3×3 matrix that doesn't describe ε_n.
- M_3'(K) would be computed against a rectangular Φ_k that doesn't carry the inter-level R_k dynamics.

But where the M_3 table was "operator wrong but inequality direction OK for large K" (closes uniformly in K if Tao A ≥ 3), **the R_k M_3' table is operator-wrong AND inequality-direction-bad** at the most natural multiplicity-floor regime: closure A requirement DIVERGES as K grows.

## 7. Caveats

(A1) **All M_3'(K) values are anticipated, not measured.** Phase 2 specifications (Approaches A and C) defer Python execution to main thread. If actual execution shows σ_1² < 1/2 cleanly + bulk well-separated, the multiplicity floor doesn't apply and Case (i) or (ii) holds.

(A2) **Tao A still parameterized.** Same caveat as M3_CLOSURE_TABLE.md §8 (A2): Tao Prop 1.17's effective C_A is INFEASIBLE this session.

(A3) **R77.4 erratum §1's articulation may admit a different reading** (Reading A: projective limit on Ẑ_3^×) where the operator IS a single self-map on L²(Ẑ_3^×) and the table changes structurally. That requires substantial separate construction (R77.5 §7 / §10) — out of scope for this probe.

(A4) **The c_k = 0 result is rigorous.** Whichever reading of R77.4 erratum §1 we choose, the structural orthogonality R_k ⊥ T(W_{k−1}) is a Q-exact identity inherited from marginal consistency. So Φ_k (the candidate transfer operator) is **definitely not** the carrier of R_k's inter-level dynamics; this is not a measurement issue.

## 8. End-to-end closure requirements with R_k

End-to-end, c=7/45 closure via R_k requires ALL of:

1. **Tao Prop 1.17 effective C_A:** INFEASIBLE this session (cf. BOOKKEEPING_PHASE1_DISPOSITION + M3_DISPOSITION + main project state).
2. **Characterized inter-level operator Φ with isolated 1/2-feature:** This probe found R_k's natural Φ_k has σ_1 → 1 dominance + bulk multiplicity near 1/2, NO isolated rate-1/2 eigenvalue.
3. **M_3'(K) bounded uniformly in K (or growing slower than K^A for Tao's effective A):** Anticipated to fail under multiplicity floor.
4. **|K| bilinear bound:** delivered (PATH2 strict 2√N at r ≤ 3; HENSEL polylog-free 2√3·√N at r ≥ 4). Sole obstruction (2) of the three that's resolved.

Items 1, 2, 3 are all open in worse ways than after the M_3 probe. **R_k is not a resolution** of obstruction (3) — it's the natural finite-truncation realization of R77.4 erratum's reframing, and the reframing reveals the obstruction is **structural, not numerical**.

## 9. Files

- `R_K_CLOSURE_TABLE.md` (this file) — parameterized closure with hypothetical M_3'(K)
- `M3_CLOSURE_TABLE.md` — parallel for falsified T_3 (better-behaved-in-K)
- `R_K_DISPOSITION.md` — top-level disposition
