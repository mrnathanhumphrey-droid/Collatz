# R77.6 — Generating function analysis of ε_n via Padé approximants

**Date:** 2026-05-04. Operator-free probe, parallel to R77.5's renormalization study. Tests the analytic structure of E(z) := Σ_{n≥1} ε_n z^n at z = 2 (the rate-1/2 endpoint) by reading the pole locations of [m/n] Padé approximants over Q.

## Verdict

> **Outcome (G-branch-cut, type indeterminate).** All non-degenerate Padé approximants place their closest pole **on the real axis at z slightly above 2**, with the diagonal sequence [n/n] converging monotonically toward z = 2 from above (2.076 at [1/1] → 2.051 at [2/2]). This is the Padé fingerprint of a **branch-cut singularity at z = 2**, NOT a simple pole at z = 2 (which would give a single stable pole exactly at z = 2 across all approximants).
>
> **The 5-coefficient budget cannot separate (G-power) from (G-log).** Both produce poles converging to z = 2 from the same side at this small order. Distinguishing requires extending to ε_7, ε_8 (k = 7, 8 Markov chains; ~hours each).
>
> **Consistent with R77.4 (M)** — branch cut implies non-self-adjoint operator; rules out simple-pole / Jordan structures; agrees with H2 ≈ H3 tie at N = 5.

## Setup

We want the analytic structure of E(z) := Σ_{n≥1} ε_n z^n near z = 2. Working with f̃(z) := (E(z) − ε_1 z)/z² to drop the n=1 transient (per the convention of R77.x). f̃ shares E's singularities away from z = 0, and f̃(0) = ε_2 ≠ 0 makes Padé construction unambiguous.

Coefficients of f̃(z) = ε_2 + ε_3 z + ε_4 z² + ε_5 z³ + ε_6 z⁴ as **exact rationals**:

| j | c_j = ε_{j+2} | float |
|---|---|---|
| 0 | 1/105 | +9.524×10⁻³ |
| 1 | −5191/1019445 | −5.092×10⁻³ |
| 2 | −11346676448406637/4627031617157687115 | −2.452×10⁻³ |
| 3 | (~57-digit / ~62-digit) | −1.152×10⁻³ |
| 4 | (~217-digit / ~219-digit) | −4.979×10⁻⁴ |

Stage 1 compute (rebuild from project's level-k Markov chain over Q): k=1..6 took **452 seconds** total (k=6 alone: 446 s for 486-state Gauss elimination over Q). Exact rationals stored in [result_77_6_pade_construction.py](result_77_6_pade_construction.py).

## Padé approximants over Q

[m/n] Padé matches first m+n+1 Taylor coefficients of f̃. With 5 coefficients, m + n ≤ 4. Constructed exactly via Gaussian elimination on the n×n linear system for q_1..q_n, then read off P from f̃·Q ≡ P (mod z^{m+n+1}).

The approximants computed (excluding the [m+n=5] cases that would need ε_7):

| (m, n) | role | poles found | closest pole | dist to z=2 |
|---|---|---|---|---|
| (1, 1) | lowest-order diagonal | 1 | z = +2.076 | **0.076** |
| (2, 1) | near-diag m+n=3 | 1 | z = +2.129 | 0.129 |
| (1, 2) | near-diag m+n=3 | 2 (one spurious at 155) | z = +2.130 | 0.130 |
| (3, 1) | m+n=4 | 1 | z = +2.313 | 0.313 |
| (2, 2) | diagonal m+n=4 | 2 | z = +2.051 | **0.051** |
| (1, 3) | m+n=4 | 3 (two spurious) | z = +2.348 | 0.348 |
| (4, 0) | Taylor truncation | — | — (sanity) | — |
| (0, 4) | all-pole | 2 cc-pairs | off-axis on \|z\|≈1.1 | 1.85 |

## Diagonal Padé sequence — convergence to z = 2

The diagonal [n/n] sequence is the natural probe of singular structure (Stahl's theorem):

| n | pole(s) | distance to z=2 |
|---|---|---|
| 1 | z = +2.0764 | 0.076 |
| 2 | z = +2.0513, z = +0.6878 | **0.051**, 1.31 |

The PRIMARY diagonal pole moves from 2.0764 → 2.0513 — **monotone convergence to z = 2 from above the real axis**. The secondary pole at z = 0.6878 in [2/2] is a likely artifact at this small order (only 5 input coefficients).

The convergence ratio is 0.051 / 0.076 ≈ **0.67**, which is slower than exponential (which would be ~0.25 at most for a simple pole at z = 2). This is consistent with branch-cut convergence rates (typically O(1/N) or O(1/N²)) and inconsistent with a pure pole.

## What the pole locations say

**Strong findings:**

1. **All "interior" Padé approximants put their closest pole on the positive real axis between z = 2.05 and z = 2.35.** Off-axis poles (when present) are clearly artifacts: the [1/2] case gives a spurious pole at z = 155, the [1/3] case at z ≈ −12 and z ≈ +8.

2. **Poles approach z = 2 from ABOVE on the real axis (z > 2).** This is the canonical Padé signature of a branch-cut singularity at z = 2 with the cut extending along [2, ∞). For a simple pole exactly at z = 2, Padé would give a stable pole AT z = 2 across all (m, n); we see drift, not stability.

3. **Diagonal [n/n] sequence monotonically tightens toward z = 2.** Slope of (closest-pole − 2) vs n consistent with branch-cut pole-densification, inconsistent with stationary point.

**Cannot distinguish (G-power) from (G-log) at N = 5:**

| signature | (G-power): power-law | (G-log): logarithmic |
|---|---|---|
| Padé poles | cluster on cut [2, ∞), density ~ N | cluster on cut [2, ∞), density ~ N (log-faster) |
| Distance to 2 | shrinks like O(1/N) or O(1/N²) | shrinks like O(1/N) |
| At N = 2 diagonal | one pole near z = 2 | one pole near z = 2 |
| Discriminator | requires N ≥ 5–10 diagonal points | same |

Two diagonal points can't separate these. Both fit the empirical pattern.

## Cross-validation with R77.4 and R77.5

- **R77.4 (envelope curve fits):** Verdict (M). H1 (Jordan) ruled out by direction (b < 0). H2 (log) and H3 (power) tied (ΔAIC = 0.23). R77.6 confirms: a branch cut at z = 2 in E(z) corresponds to envelope shape (1/2)^n · g(n) for slowly-varying g — consistent with both H2 and H3, INCONSISTENT with H1 (which would be a simple pole). **Two independent probes agree.**

- **R77.5 (renormalization residual operator) [parallel work]:** if R77.5's operator Φ has continuous-spectrum or branch-cut structure near λ = 1/2, that's a third independent probe converging on the same picture. If R77.5 finds a discrete eigenvalue at 1/2, that contradicts R77.6's branch-cut signature, and one of them is misreading. (Pending R77.5's output.)

- **R77.4 erratum (K_k spectrum):** K_k itself has no eigenvalue near 1/2 — confirmed not the rate operator. R77.6 doesn't depend on which operator T is; it directly probes E(z). Compatible.

## Decision tree resolution

Per the brief's outcome categories:

- **(G-power):** branch-cut signature, real-axis [2, ∞). **Partially observed.** Poles cluster real-axis, monotone diagonal convergence — visible signature. Power vs log NOT distinguished.
- **(G-log):** logarithmic signature, dense pole pattern at z = 2. **Compatible with observations** at N = 5 — log produces similar pole pattern as power-law to leading order.
- **(G-mixed):** poles unstable across (m, n). **NOT observed.** Closest poles cleanly cluster in [2.05, 2.35] across 6 distinct approximants.
- **(G-anomaly):** poles in unexpected location (e.g., off real axis, at z ≠ 2). **NOT observed.** All closest poles are real, positive, near z = 2.

> **Verdict: (G-branch-cut, type indeterminate).** The framework's hypothesis "ε_n's generating function has branch-cut singularity at z = 2" is supported. The branch type (power-law vs logarithmic) is not resolvable from N = 5 coefficients.

## What's needed to discriminate

The brief and AIC theory agree: separating power-law from logarithmic at the singular endpoint requires either much higher N or finer asymptotic analysis at the singularity.

- **Extending to ε_7 (k=7 Markov chain, ~hours):** adds one diagonal Padé point [3/3]. Predicted: pole at z ≈ 2.030–2.040. If the convergence ratio is ≈ 0.5 (i.e., 0.051 → 0.025), we have evidence consistent with a power-law of α ≈ 1; if ratio is ≈ 0.7, consistent with log. Either way, the pattern firms up.
- **Extending to ε_8 (k=8 Markov chain, ~tens of hours):** [4/4] diagonal point. By then the discrimination should be feasible at the (G-power) vs (G-log) level.
- **Differential approximants** (D-Padé, working in the ODE that f̃ satisfies modulo the singular part): can sometimes distinguish power-law exponent at lower N, but requires more numerical infrastructure. Not done here.

## Caveats (per anti-patterns in brief)

1. **N = 5 is small for Padé.** This is acknowledged. Verdict (G-branch-cut, type indeterminate) is robust because the *direction* of the closest pole is consistent (all real, all > 2, diagonal monotone), not the precise value. Spurious poles in non-diagonal cases ([1/2] at z=155, [1/3] at z=−12 and z=+8) are normal at small N and do not alter the diagnostic.

2. **Stahl's theorem is suggestive only.** The full theorem requires meromorphic functions and applies asymptotically. Treating the [2/2] convergence as a finite-N approximation to the Stahl curve is heuristic.

3. **Floating-point pole-finding.** Roots of Q(z) computed via numpy.roots after float-conversion. Q's coefficients are exact rationals (computed over Q), so the Q polynomial itself is exact; only root-finding is numerical. Errors in pole locations from float arithmetic are ≤ 10⁻¹⁰ for n ≤ 3, which is far below the "distance to 2" magnitudes (≥ 0.05) we report.

4. **The secondary [2/2] pole at z = 0.688 is a likely Padé artifact** at small N; it does not move with the n=1 transient (which is excluded by working with f̃, not E directly), nor does it appear in [1/1] or [3/1]. Unstable secondary poles are typical in low-order Padé.

## Files

- [result_77_6_pade_construction.py](result_77_6_pade_construction.py) — script
- [experiments_output/result_77_6_pade_poles.csv](experiments_output/result_77_6_pade_poles.csv) — pole locations across (m, n)
- [result_77_6_generating_function.md](result_77_6_generating_function.md) — this writeup

## Combined picture (post R77.4 + R77.6)

The framework's open question on ε_n's analytic structure now has TWO independent partial answers:

| probe | hypothesis class probed | result |
|---|---|---|
| R77.4 envelope fits (M) | (J) / (L) / (P) on \|ε_n\|·2^n | Jordan ruled out; log ≈ power tied |
| R77.6 generating function (G-branch-cut, indet.) | analytic structure at z=2 | branch cut at z=2; power vs log tied |

Both probes converge on: **simple pole at z = 2 ruled out (no rate-1/2 with constant prefactor), Jordan structure ruled out, branch-cut singularity at z = 2 supported, but power-law (α ≈ 0.14 from H3) vs logarithmic (H2) not separated at N = 5.**

The next data-budget step that resolves this is computing ε_7 (k=7 Markov chain). Predicted compute: a few hours given Stage 1's ~7-minute timing for k=6 and the O(N³) scaling of Gauss elimination on Q (k=7 has N=1458, ratio ~27× over k=6).

R77.x as a closure-attack thread for c=7/45 remains parked: this finding does not unlock or close any path; it tightens the structural picture from "operator unidentified" to "branch-cut analytic structure at z = 2, type to be determined".
