# TAUBERIAN_RESCOPE_B_M_PARAMETER (Chevalier 1.16 M-value computation)

**Date:** 2026-05-13.

Theorem 1.16's conclusion: b_n = D · n^{M − 3/2} (1 + d_1/n + … + d_K/n^K + O(1/n^{K+1})).

The M parameter is the multiplicity of the pole of h_p at 0. To check if "M is realizable", we need to back out M from the empirical b_n behavior:

  n^{M − 3/2} fits the rate of decay/growth.

## Empirical M backout — input (1) candidates

### T2: b_n = |ε_n|·2^n, n = 1..8

If b_n ~ D · n^{M − 3/2}, then log b_n ~ const + (M − 3/2) log n.

| n | b_n | log n | log b_n | implied M − 3/2 from n=2 ratio |
|---|---|---|---|---|
| 1 | 0.4000 | 0.0 | −0.916 | — |
| 2 | 0.0381 | 0.693 | −3.267 | (−3.267 + 0.916)/0.693 = −3.39, so M ≈ −1.89 |
| 3 | 0.0407 | 1.099 | −3.200 | — |
| 4 | 0.0392 | 1.386 | −3.238 | — |
| 5 | 0.0369 | 1.609 | −3.299 | — |
| 6 | 0.0319 | 1.792 | −3.444 | — |
| 7 | 0.1504 | 1.946 | −1.894 | ← jump |
| 8 | 0.1909 | 2.079 | −1.656 | — |

The pattern is **not** a clean power law. From n=2 to n=6, b_n is approximately constant (~0.035 ± small wobble) — implying M − 3/2 ≈ 0, i.e. **M ≈ 3/2**, **but M is required to be a positive integer**. M = 1 gives n^{-1/2} = decreasing; M = 2 gives n^{1/2} = increasing; observed b_n is approximately constant for n=2..6, which is neither.

Then at n=7, b_n jumps to 0.15, breaking even the "approximately constant" pattern.

**No integer M ≥ 1 fits the observed b_n behavior.** Theorem 1.16's M parameter is **NOT realizable** from this data.

### T1: b_n = ε_n (signed)

b_n alternates in sign with sign pattern +, +, −, −, −, −, −, − — Theorem 1.16's asymptotic b_n = D n^{M − 3/2} (1 + d_1/n + …) predicts a *fixed sign* eventually (sign of D). The actual b_n changes sign at k=3 once and then is consistently negative. So the sign pattern is *consistent with* a single-sign asymptotic for n ≥ 3 with D < 0, BUT the magnitudes don't follow a clean n^{M-3/2}:

| n | |ε_n| | log n | log |ε_n| |
|---|---|---|---|
| 3 | 5.09e-3 | 1.099 | −5.28 |
| 4 | 2.45e-3 | 1.386 | −6.01 |
| 5 | 1.15e-3 | 1.609 | −6.77 |
| 6 | 4.98e-4 | 1.792 | −7.60 |
| 7 | 1.18e-3 | 1.946 | −6.75 (jump) |
| 8 | 7.46e-4 | 2.079 | −7.20 |

Slope from n=3..6: (−7.60 − (−5.28))/(1.792 − 1.099) = −3.35. So if |ε_n| ~ n^α this gives α ≈ −3.35, i.e. M − 3/2 ≈ −3.35, i.e. **M ≈ −1.85**.

But Theorem 1.16 requires **M ≥ 1 integer**. **M < 0 is OUT OF RANGE.** Theorem 1.16's hypothesis violated.

(Note: M < 0 would mean h_p has no pole at 0, only a zero — but then it would be holomorphic and we'd be in Theorem 1.14 not 1.16. Theorem 1.14 gives b_n = (1/n^{3/2}) (C + …), i.e. fixed exponent −3/2. The observed n^{-3.35} decay is *steeper than n^{-3/2}*, which Theorem 1.14 cannot deliver either.)

---

## What if we try a multi-singularity / multiple-square-root version?

Theorem 1.16 is for a *single* pole at 0 of h_p. The k=7 jump is suspiciously consistent with a *second* singularity contributing — i.e., the generating function has multiple dominant singularities and the asymptotic is a sum of two terms with different α exponents. **But this is Flajolet-Sedgewick Theorem VI.5 (multiple singularities), not Chevalier 1.16.** Chevalier 1.16 cannot accommodate the multi-singularity case in its single-pole-at-0 hypothesis.

---

## M-parameter verdict

**M is NOT realizable as a positive integer ≥ 1 from inputs (1)-(4).**

Empirical slope on n=3..6 gives M ≈ −1.85 (out of theorem's range M ≥ 1).
At n=7,8 the slope changes sign (b_n grows), contradicting any monotone-decay asymptotic n^{M−3/2} with M ≤ 3/2.

Theorem 1.16 disposition: **BLOCKER (parameter unrealizable + insufficient data + Mode H trap on h_p)**, which combined with the HYPOTHESIS_CHECK BLOCKER for h_4-h_6 gives final disposition **BLOCKER** for B.

---

## Note: what data would unblock B?

To realize Chevalier 1.16 we would need:
(a) ε_9, ε_10, …, ε_K (K ≥ ~15-20) computed exactly to extend the asymptotic fitting.
(b) An *a priori* analytic-continuation argument for g(z) = Σ ε_n z^n past |z| = 1.
(c) An *a priori* identification of h_p (or at least its pole multiplicity M) coming from the dynamics, not from numerical fit.

(a) alone is insufficient: even with 20 terms the fitting still numerical-only would not satisfy h_4 (existence of meromorphic h_p) verbatim.

Compute is free (per brief), but extending to ε_9..ε_20 requires the v2_vec_pool solver — a non-trivial multi-week compute task, and even completion of (a) does not address (b) or (c). The k=7 anomaly is **structurally informative** even before more terms — it indicates a SINGLE-singularity n^{M-3/2} asymptotic does not fit.
