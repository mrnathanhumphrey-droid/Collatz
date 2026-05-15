# TAUBERIAN_RESCOPE_B_HYPOTHESIS_CHECK (Chevalier 1.16 × inputs (1)-(4))

**Date:** 2026-05-13.

For each candidate target sequence b_n (from §"Notational mapping" in B_HYPOTHESES), check each hypothesis against each input.

---

## Empirical input (1) pattern, refined

Normalized magnitudes |ε_k|·2^k:

| k | ε_k | |ε_k|·2^k | ratio |ε_{k+1}/ε_k| |
|---|---|---|---|
| 1 | +0.20000 | 0.4000 | — |
| 2 | +0.00952 | 0.0381 | 0.048 |
| 3 | −0.00509 | 0.0407 | 0.535 |
| 4 | −0.00245 | 0.0392 | 0.482 |
| 5 | −0.00115 | 0.0369 | 0.470 |
| 6 | −0.000498 | 0.0319 | 0.432 |
| 7 | −0.00118 | 0.1504 | 2.360 ← JUMP |
| 8 | −0.000746 | 0.1909 | 0.634 |

For k=2..6, ratio |ε_{k+1}/ε_k| ∈ [0.43, 0.53] suggesting |ε_k| ~ C · 2^{-k}. **At k=6→7 the ratio jumps to 2.36**; |ε_7| is *larger* than |ε_6|. The geometric-decay hypothesis with rate 1/2 is FALSIFIED at k=7. (Sign pattern: +, +, −, −, −, −, −, − — also informational.)

---

## Target sequence T1: b_n = ε_n (signed exact rationals)

| Hypothesis | Status | Reason |
|---|---|---|
| h_1: g continuous on D̄ | UNVERIFIABLE | We only have 8 coefficients; need analytic-continuation argument from the *full* ε_k structure to extend g across |z|=1. We do not have such an argument from inputs (1) alone. |
| h_2: g holomorphic on D | UNVERIFIABLE | Same: need radius of convergence ≥ 1 (i.e. lim sup |ε_n|^{1/n} ≤ 1) — observed 8-term data is consistent with radius > 1 (numerically lim sup of |ε_k|^{1/k} for k=1..8 gives roughly 0.5; ratio test → ~0.5 for k≤6) but k=7 jump muddies this. UNVERIFIABLE without convergence proof on the full sequence. |
| h_3: Σ ε_n z^n has positive radius | LIKELY TRUE but FRAGILE | First 8 terms give finite estimates consistent with radius ≥ 2, but the k=7,8 pattern jump is unmodeled. |
| h_4: ∃ meromorphic h_p on nbhd of D(1,1)^{1/2} with single pole at 0 of mult M ≥ 1 | UNVERIFIABLE | The substitution g(z) = h_p(√(1−z)) is the load-bearing claim. To find h_p we'd need to know g(z) globally; with only 8 coefficients we cannot identify h_p. |
| h_5: M ≥ 1 finite, integer | UNVERIFIABLE → BLOCKER | The whole exercise is to *find* M from the asymptotic. We don't even have evidence the sequence ε_n admits an n^{-(3/2 − M)} asymptotic — the k=7 jump *contradicts* a clean n^{α} asymptotic. |
| h_6: g(z) = h_p(√(1−z)) for all z ∈ D | UNVERIFIABLE | Same as h_4. |

**T1 verdict:** **BLOCKER** — input (1) gives 8 numerical coefficients with an anomaly at k=7, insufficient to verify *any* of the analytic-continuation hypotheses. Even ignoring the k=7 anomaly, h_4 and h_6 require knowing g globally, not just its first 8 Taylor coefficients.

---

## Target sequence T2: b_n = |ε_n| · 2^n (normalized magnitudes)

Same analysis — replacing ε_n with |ε_n|·2^n only changes the radius of convergence (now ≥ 1 since |ε_k|·2^k ≤ 1 empirically for k ≤ 8). h_4, h_5, h_6 still UNVERIFIABLE for the same reason.

The empirical *non-monotonicity* (jump at k=7) of |ε_k|·2^k is **active evidence AGAINST h_6**: an asymptotic of the form b_n ~ D n^{M − 3/2} (1 + d_1/n + ...) would predict b_n monotone in n for large n, of consistent sign. The k=7 jump and sign-mixing makes the b_n ~ D n^{M-3/2} fit *poorly* on the observed range.

**T2 verdict:** **PARTIAL→NO_FIT** — h_1, h_2 plausible; h_4, h_5, h_6 unverifiable AND empirically disfavored.

---

## Target sequence T3: b_n = |μ̂_n(ξ)|² (the closure target)

| Hypothesis | Status | Reason |
|---|---|---|
| h_1, h_2, h_3 | UNVERIFIABLE | Need μ̂_n(ξ) asymptotic — that's the closure target itself. **Mode H circularity.** |
| h_4, h_5, h_6 | UNVERIFIABLE / Mode H | Same. |

**T3 verdict:** **BLOCKER (Mode H circular)** — the target sequence is the closure target; using Chevalier 1.16 to bound it requires already knowing its analytic structure.

---

## Target sequence T4: b_n = |ε_n|² (squared)

Same issues: 8 coefficients, k=7 jump.

**T4 verdict:** **BLOCKER**.

---

## Phase 1 (h × I) matrix — Chevalier 1.16

|  | (1) ε_k k=1..8 | (2) C1 renewal-walk | (3) C2 BMP F_1 support | (4) BT archimedean place |
|---|---|---|---|---|
| h_1 (g continuous on D̄) | UNVERIFIABLE | N/A | N/A | N/A |
| h_2 (g holomorphic on D) | UNVERIFIABLE | N/A | N/A | N/A |
| h_3 (power series at 0) | UNVERIFIABLE / partial | N/A | N/A | N/A |
| h_4 (∃ meromorphic h_p) | UNVERIFIABLE | N/A | N/A | N/A |
| h_5 (M ≥ 1 finite integer) | UNVERIFIABLE → empirically suggestive of NO clean M | N/A | N/A | N/A |
| h_6 (g(z) = h_p(√(1-z))) | UNVERIFIABLE | N/A | N/A | N/A |

Inputs (2), (3), (4) are **structural** descriptions (renewal-walk / model-set / archimedean-place) that do not directly enter Theorem 1.16's hypothesis slots — they characterize the *substrate* on which μ_n lives but Theorem 1.16's hypotheses are about the *generating function* on the unit disc. So they bear on h_3-h_6 only indirectly through what they say about ε_n's analytic structure.

**Aggregate disposition for B: BLOCKER.**

Reason: Theorem 1.16's hypotheses h_4, h_5, h_6 require *a priori* knowledge of the generating function's square-root profile h_p (with explicit pole at 0 of multiplicity M). Inputs (1)-(4) do NOT supply h_p:
- (1) gives only 8 numerical coefficients with an empirical anomaly at k=7 that *defeats* a clean (1−z)^{1/2 − M} asymptotic.
- (2) renewal-walk gives the recursion structure but not the generating-series analytic continuation.
- (3) BMP gives the support diffraction (a *different* analytic object — Fourier on ℝ, not a generating function on D̄).
- (4) BT archimedean-place finding tells us the c=7/45 phenomenon is global/adelic, not encoded in any single p-adic generating series.

**The M parameter is therefore unknown.** This is the BLOCKER case.

See M_PARAMETER file for further analysis.
