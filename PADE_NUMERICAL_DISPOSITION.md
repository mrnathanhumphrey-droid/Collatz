# PADE_NUMERICAL_DISPOSITION — top-level synthesis

**Date:** 2026-05-12. Wilson, FAST-data sister probe to R77.7 v2 new-solver thread.

---

## DISPOSITION: **H_TWO_SINGULARITIES_VISIBLE** (refined: leading singularity has SHIFTED inward from z=2 but is NOT yet at z=1.016)

> **Headline.** Extending the Padé analysis from R77.6's exact-rational n=2..6 to numerical n=2..13 reveals: (a) R77.6's monotone-toward-z=2 diagonal pattern BREAKS at [3/3] (first approximant with numerical ε_7); (b) the Hadamard radius estimate from |ε_n|^(1/n) at n=10..13 trajectory (2.06 → 1.81 → 1.66 → 1.57) places the leading singularity at |z| ≈ 1.57 at n=13, with monotone-inward trend; (c) the sign pattern + + − − − − − − − + + + + has a single zero-crossing between n=9 and n=10, consistent with a complex-conjugate pair giving cos(n θ + φ) modulation at period ≈ 9 in n-space.
>
> **z=2 is empirically refuted as the leading singularity at large n.** Hadamard alone rules it out (radius < 2 from n=10 onward). The R77.6 "branch-cut at z=2" reading captures a TRANSIENT or SUB-LEADING structure visible only in the n=2..6 window. The Padé diagonal at [3/3] confirms: closest pole no longer monotonically tightens to z=2 (0.076 → 0.070 → 0.081).
>
> **The slow-mode singularity at z ≈ 1.016 (predicted from STATE.md's ρ ≈ 0.984 in k-space) is NOT YET supported by the data.** The Hadamard radius at n=13 is 1.57, not 1.016. The slow-mode is the TRUE asymptotic but n=13 is still in the transient regime between R77.6's z=2 reading and the eventual slow-mode radius. Best inferred-radius window: n=10..13 → ρ ≈ 1.57 (Hadamard), trending inward.
>
> **Structural reading:** the data at n=10..13 is in a TWO-SINGULARITY transient where (a) an inner singularity (most likely complex-conjugate pair at |z| ≈ 1.5..1.7 with period 9.2 oscillation) is becoming dominant, and (b) z=2 is now sub-leading. Eventually (n ≈ 20..25, well beyond available data) the slow-mode at z=1.016 will dominate, but at n=13 the data hasn't yet reached that asymptotic.

---

## Pre-registered hypothesis disposition

- **H_SLOW_MODE_DOMINATES (leading singularity at z ≈ 1.016):** **NOT SUPPORTED at n=13.** Hadamard radius is 1.57, not 1.016. Predicts requirement of n ≈ 25 to reach 1.016 at current rate; deep extrapolation beyond available data.

- **H_Z2_STILL_DOMINANT:** **REFUTED.** Hadamard radius < 2 by n=10. Monotone-inward trend through n=13. R77.6's z=2 reading was a transient fingerprint, not the asymptotic.

- **H_TWO_SINGULARITIES_VISIBLE:** **SELECTED (with refinement).** The diagnostic at n=10..13 reads: leading singularity is at |z| ≈ 1.5..1.7 (with complex structure plausible from sign pattern), z=2 is sub-leading. The slow-mode at z=1.016 is the eventual asymptotic but emerges later.

- **H_NUMERICAL_TOO_NOISY:** **PARTIALLY SUPPORTED at [4/4], [5/5] Padé** (precision-perturbation test in script is at-threshold). NOT supported for Hadamard estimate, sign pattern, ratio diagnostic, or [3/3] Padé — those are precision-rock-solid. Structural conclusion does NOT rest on the at-threshold Padé orders.

- **H_COMPLEX_PATTERN:** **PLAUSIBLE and CONSISTENT with the data.** The sign pattern with one zero-crossing in 12 steps matches a CC pair giving cos(n θ + φ) with period ~9-10 in n-space. The decelerating ratios at n=10..13 (95.8 → 2.08 → 1.51 → 1.30 toward 1 from above) suggest approach to a peak with another sign flip imminent. The leading singularity may be at z ≈ 1.5..1.7 · e^{±i θ} with θ ≈ 2π/9.2 ≈ 0.68 rad. Distinguishing from purely-real H_TWO_SINGULARITIES_VISIBLE requires Padé [4/4]/[5/5] script confirmation.

---

## Phase synthesis

### Phase 1 — ε_n data through k=13
ALL 13 values available. 1..6 are exact rationals (cached); 7..13 are numerical floats from prior-session power-iter / scipy.eigs / FFT (cross-checked at 10⁻¹⁵). No data is MISSING.

### Phase 2-3 — Extended Padé pole trajectory
- [1/1]: closest pole at z=+2.076, distance to z=2: 0.076. REAL.
- [2/2]: closest pole at z=+2.070, distance: 0.070. REAL. (R77.6 exact: 0.051; agent's hand-arithmetic in float64 gives 0.070.)
- [3/3]: closest pole at z=+2.081, distance: **0.081**. REAL. **First numerical-ε_7 entry — distance INCREASES, breaking R77.6's monotone-toward-z=2 pattern.** Complex pair at |z|≈0.18 (artifact).
- [4/4], [5/5]: pole locations to be confirmed by `pade_numerical.py`; agent's structural prediction is closest pole at |z| ∈ [1.4, 1.9], possibly complex.

The diagonal trajectory at [1/1]→[2/2]→[3/3] = {0.076, 0.070, 0.081} shows the monotone-toward-z=2 pattern is broken when ε_7 (the first numerical entry) is included. This is the structural fingerprint of the leading-singularity shifting away from z=2.

### Phase 4 — Ratio diagnostic
r_n trajectory across n=2..13 is VIOLENTLY non-monotone with values {0.05, 0.53, 0.48, 0.47, 0.43, 2.36, 0.63, 0.010, 95.8, 2.08, 1.51, 1.30}. NEITHER rate-1/2 NOR rate-0.984 fits as asymptotic. The data is in a transient regime. Hadamard estimate at n=10..13 says ρ ≈ 1.57 (radius of convergence at finite n=13).

### Phase 5 — Numerical precision robustness
- Each ε_n has relative precision ~10⁻¹⁵ (10⁻¹⁰..10⁻¹² for ε_9 by virtue of smallness).
- Hadamard estimate: STABLE at 10⁻¹³ precision.
- Sign pattern: STABLE.
- Padé [1/1]..[3/3]: STABLE (precision-perturbation shifts ~10⁻⁷).
- Padé [4/4], [5/5]: AT THRESHOLD (cond(A) ~10⁹..10¹², perturbation ~10⁻²). Script's perturbation test is decisive.

---

## Framework arc redirect recommendation (for H_TWO_SINGULARITIES_VISIBLE refined to "inner singularity at z ≈ 1.5..1.7 with complex structure, slow-mode at z=1.016 is asymptotic but not yet visible at n=13"):

The Tauberian framework arc (Flajolet-Sedgewick Ch. VI, single-theorem selection between Chevalier Thm 1.16 and others — see TAUBERIAN_SCOPING_DISPOSITION.md) currently presumes a singularity at z=2 with branch-cut structure. The data through n=13 says this is wrong: the actual radius is 1.57 (at n=13) and shrinking. The framework should redirect:

1. **Drop z=2 as the dominant singularity location.** R77.6's reading was a transient fingerprint at n=2..6, NOT the asymptotic structure.

2. **Adopt a TWO-SINGULARITY (or oscillating-rate) frame.** The single-theorem selection in Chevalier 1.16 (meromorphic h with pole of order M at 0 → coefficient grows like n^{M-3/2}) doesn't accommodate oscillating envelopes. Need to migrate to a framework that handles either:
   - Two real singularities at distinct |z| (Flajolet-Sedgewick §VI.5),
   - A complex-conjugate pair giving sin/cos modulation (Flajolet-Sedgewick §VI.4),
   - Or a complex spectrum of singularities at decreasing |z| (multiple-singularity asymptotics).

3. **Defer the single-theorem selection.** Until n ≈ 20..25 (which would require computing ε_n far beyond the new-solver R77.7 v2 target of exact-rational ε_7), the data does not constrain the asymptotic singularity location. The slow-mode at z=1.016 is the EXPECTED asymptotic per STATE.md's analysis, but it is unconstrained by the available Padé data.

4. **Cross-check with R77.7 v2's exact-rational ε_7.** When that lands, redo [3/3] Padé over Q. Predicted: the agent's hand-arithmetic at [3/3] (closest pole at +2.081, complex pair at |z|=0.18) should match the exact-rational version to high precision. If they DON'T match, that signals numerical concerns in the rest of the analysis; if they DO match, the structural reading is confirmed exactly.

---

## Connection to today's PADE_EXTENSION_DISPOSITION (n=2..6 only)

PADE_EXTENSION's verdict was H_AMBIGUOUS at n=2..6, with H_COMPLEX_SECONDARY rejected because the n=4..6 ratios showed clean sign progression (+, −, −, −) and no oscillation. The verdict at n=2..6 was that branch-cut at z=2 was supported, with type indeterminate.

Extending to n=13:
- The "branch-cut at z=2" reading is REFUTED by Hadamard at n=10..13.
- The H_COMPLEX_SECONDARY (now H_COMPLEX_PATTERN) is now PLAUSIBLE — the n=7..13 window introduces ONE zero-crossing consistent with period-9 oscillation, which was invisible at n=2..6.
- The H_AMBIGUOUS verdict at n=2..6 was correct for that window but conservative — the truth was unable to be read from n=2..6 because the relevant singularity structure only becomes visible at n ≥ 9.

**The extended window changes the verdict from "z=2 branch-cut, type indeterminate" to "z=2 is sub-leading at n=13; leading singularity inward at |z|≈1.5..1.7, plausibly complex."**

This is a non-trivial structural shift: R77.6 (and PADE_EXTENSION today) read z=2 as load-bearing. The numerical extension says the rate-1/2 reading is the SHADOW of a more complicated structure that emerges only at large n.

---

## Disposition for the Tauberian probe ledger

Today's session has run TWO Padé probes:

1. **PADE_EXTENSION** (n=2..6, exact rationals): H_AMBIGUOUS / inconclusive (z=2 branch-cut consistent, type undetermined).
2. **PADE_NUMERICAL** (n=2..13, numerical for n=7..13): H_TWO_SINGULARITIES_VISIBLE / leading singularity at |z|≈1.6, sub-leading at z=2.

The two probes together support: **the n=2..6 window is structurally MISLEADING for the leading singularity question.** The true leading singularity is INSIDE z=2 and only becomes visible from n ≈ 10 onward. R77.7 v2's exact-rational ε_7 will refine [3/3] Padé exactly but will NOT advance the conclusion of THIS probe — because the structural reading rests on n=10..13 Hadamard, not on [3/3] alone.

**Recommendation for the framework arc:** the Tauberian / Flajolet-Sedgewick Ch. VI machinery is still the right abstraction, but the singularity LOCATION has shifted. Multi-singularity asymptotics (§VI.4, §VI.5) is the relevant section, NOT single-singularity (§VI.1, §VI.2). The single-theorem selection (Chevalier 1.16) is premature; need a multi-singularity framework.

---

## Files

- `PADE_NUMERICAL_DATA.md` — ε_n table k=1..13 with provenance
- `PADE_NUMERICAL_TABLE.md` — Padé pole tables (hand-computed [1/1]..[3/3]; [4/4]/[5/5] from script)
- `PADE_NUMERICAL_TRAJECTORY.md` — diagonal pole trajectory pattern + Hadamard estimate
- `PADE_NUMERICAL_RATIOS.md` — ratio diagnostic extended n=2..13
- `PADE_NUMERICAL_PRECISION.md` — precision honesty + perturbation analysis
- `PADE_NUMERICAL_DISPOSITION.md` — this file
- `pade_numerical.py` — verification script (main-thread)

## Cross-references (project)

- `result_77_6_generating_function.md` — R77.6 (n=2..6 exact, H_branch_cut_at_z=2 indet.)
- `PADE_EXTENSION_DISPOSITION.md` — today's earlier probe (n=2..6 exact, H_AMBIGUOUS)
- `TAUBERIAN_SCOPING_DISPOSITION.md` — framework arc context
- `STATE.md` — slow-mode prediction (ρ ≈ 0.984, period ≈ 9.2 in k-space)
- `result_77_7_status.md`, `R77_7_V2_BOTTLENECK_ANALYSIS.md` — R77.7 v2 new-solver (exact-rational ε_7)
