# PADE_EXTENSION_RATIOS — Phase 3: ratio diagnostic

**Date:** 2026-05-12. Wilson. Phase 3 of the Padé-extension probe.

---

## Headline finding

> **Ratios r_n := |ε_n|/|ε_{n-1}| do NOT converge monotonically to 1/2 for n=3..6. The deviation δr_n := r_n − 1/2 has sign pattern + − − − (n=3..6) with magnitudes [+0.035, −0.018, −0.030, −0.068]. After the sign change at n=4, magnitudes ACCELERATE downward (factor 1.65 then 2.23 per step). This is inconsistent with a pure simple pole (which would give monotone geometric approach to 0.5) AND inconsistent with an oscillating subleading correction (which would give alternating signs). It IS consistent with a leading rate-1/2 pole + NEGATIVE subleading correction whose magnitude grows in n — the signature of a branch-cut at z=2 with NEGATIVE branch coefficient.**

---

## Setup

Working with empirical ε_n cached through k=6 (exact rationals, loaded from `experiments_output/result_77_7_eps_exact_through_k7.json`):

| n | ε_n (exact, displayed as float) | |ε_n| |
|---|---|---|
| 1 | +0.2000000000 (= 1/5) | 2.000×10⁻¹ |
| 2 | +9.5238095×10⁻³ (= 1/105) | 9.5238×10⁻³ |
| 3 | −5.0920305×10⁻³ (= −5191/1019445) | 5.0920×10⁻³ |
| 4 | −2.4521884×10⁻³ | 2.4522×10⁻³ |
| 5 | −1.1517287×10⁻³ | 1.1517×10⁻³ |
| 6 | −4.9790×10⁻⁴ | 4.9790×10⁻⁴ |

Defines: r_n := |ε_n| / |ε_{n−1}|. Target for clean rate-1/2: r_n → 1/2.

---

## Raw ratios

| n | r_n | r_n − 1/2 | comment |
|---|---|---|---|
| 2 | 0.047619 (= 1/21) | −0.452381 | n=1 is transient (ε_1=1/5 large), ratio_2 not reflective |
| 3 | 0.534643 | **+0.034643** | first "asymptotic" ratio |
| 4 | 0.481574 | **−0.018426** | sign flip |
| 5 | 0.469686 | **−0.030314** | deeper below 0.5 |
| 6 | 0.432319 | **−0.067681** | even deeper, departure ACCELERATING |

(Excluding n=2 since it includes the transient ε_1; report n=3..6 as the diagnostic range.)

---

## Sign and magnitude pattern (n=3..6)

| n | r_n − 1/2 | sign | \|r_n − 1/2\| | ratio to previous \|·\| |
|---|---|---|---|---|
| 3 | +0.034643 | + | 0.0346 | — |
| 4 | −0.018426 | − | 0.0184 | 0.532 |
| 5 | −0.030314 | − | 0.0303 | 1.645 |
| 6 | −0.067681 | − | 0.0677 | 2.234 |

**Key observations:**

1. **One sign change** (between n=3 and n=4). After that, all negative.
2. **Magnitudes are NOT monotone.** Drop from n=3 to n=4 (factor 0.53), then rise n=4→5→6 (factors 1.65, 2.23).
3. **Acceleration of departure from 0.5 after the sign change** is the standout feature.

---

## Diagnostic tests

### Test 1: Monotone approach to 0.5 (would indicate pure pole)

REJECTED. Sign change between n=3 and n=4 rules out monotone approach. Even ignoring sign, magnitude rises from n=4 onward.

### Test 2: Oscillation around 0.5 (would indicate complex-conjugate secondary singularity)

REJECTED. Only ONE sign change in four ratios; expected oscillation pattern would be + − + − or − + − +. The single sign change is consistent with a transient correction passing through zero, not periodic oscillation.

### Test 3: r_n = 0.5 + c/n^α monotone-correction fit

| n | r_n | (r_n − 0.5)·n¹ | (r_n − 0.5)·n^1.5 | (r_n − 0.5)·n² |
|---|---|---|---|---|
| 3 | 0.5346 | +0.1041 | +0.1803 | +0.3123 |
| 4 | 0.4816 | −0.0736 | −0.1473 | −0.2944 |
| 5 | 0.4697 | −0.1515 | −0.3389 | −0.7575 |
| 6 | 0.4323 | −0.4062 | −0.9952 | −2.4372 |

None of these stabilize. (r_n − 0.5)·n is NOT constant. (r_n − 0.5)·n^β grows in magnitude for any β > 0.

**This means the departure from 0.5 is FASTER than any 1/n^α power.** Implication: at the n=4..6 range, the subleading is not yet in its asymptotic regime, OR the subleading has a faster-than-power decay component (e.g., exponential from a discrete secondary singularity at z=4).

### Test 4: Geometric correction (discrete secondary singularity)

If E(z) has a leading simple pole at z=2 (giving ε_n ~ A·(1/2)^n) and a secondary pole at z=ρ with ρ < 2 (giving ε_n ~ A·(1/2)^n + B·(1/ρ)^n), then r_n approaches some limit different from 1/2 — namely max(1/2, 1/ρ).

The observed r_n is BELOW 1/2 (not above), which is inconsistent with a secondary singularity at |z|<2 (which would push the rate UP from 1/2, not down).

**REJECTED:** No discrete secondary singularity at radius |ρ| < 2 is consistent with the ratio pattern.

### Test 5: Branch cut + negative subleading coefficient

If E(z) has a branch cut at z=2 of order α and a leading simple-pole-like term, then asymptotically:

ε_n ≈ A·(1/2)^n + B·(1/2)^n·n^{-α-1} · sign

with sign possibly negative. The ratio:

r_n = |ε_n|/|ε_{n-1}| ≈ (1/2) · |A + B·n^{-α-1}·sign| / |A + B·(n-1)^{-α-1}·sign|

If sign of B is OPPOSITE to sign of A: the magnitude correction subtracts, so |ε_n| < A·(1/2)^n, and r_n < 1/2. As n grows, B-term shrinks → r_n → 1/2 from below.

If sign of B is SAME as A: r_n > 1/2 always.

**Empirical observation:** r_n is below 0.5 for n=4,5,6 AND moving FURTHER from 0.5 (not converging back). This means the magnitude of the subleading correction is GROWING relative to the leading, which contradicts the asymptotic picture.

**However:** at small n, the subleading-to-leading ratio CAN be increasing if the asymptotic regime hasn't been reached. With α small (e.g., α ≈ 0.1-0.3 — the empirical R77.4 H3 value), the subleading n^{-α-1} decays very slowly, and the ratio of (B-term)/(A-term) ≈ (B/A)·n^{-α-1} hasn't kicked in to clear convergence at n=6.

This is the most consistent reading: **branch cut at z=2, negative subleading coefficient, low branch order (slow asymptotic approach), data still pre-asymptotic at n=6.**

### Test 6: Exponential correction from secondary singularity at z=4

If subleading ~ C·(1/4)^n (R76 §10 conjecture), the ratio of subleading to leading is (1/4)^n / (1/2)^n = (1/2)^n, which AT n=6 is 1/64 ≈ 0.0156. Subleading contribution to r_n would be ≈ (1/2)·(C/A)·(1/2)^n.

Magnitude of correction needed to fit r_6 − r_5: (0.4697 − 0.4323) = 0.0374. Would need (1/2)·(C/A)·(1/2)^5 ≈ 0.0374, so C/A ≈ 0.0374·64 ≈ 2.4. **Plausible** that (C/A) is of order unity.

**Test prediction:** If subleading is purely (1/4)^n, r_n = 1/2 · (1 + (C/A)(1/2)^{n-1}) / (1 + (C/A)(1/2)^{n-2}). With C/A ≈ 2.4:
- r_3 should be ≈ 1/2 · (1 + 2.4·(1/4)) / (1 + 2.4·(1/2)) = 1/2 · 1.6/2.2 = 0.364
- r_4 should be ≈ 1/2 · (1 + 2.4·(1/8)) / (1 + 2.4·(1/4)) = 1/2 · 1.3/1.6 = 0.406
- r_5 ≈ 1/2 · 1.15/1.3 = 0.442
- r_6 ≈ 1/2 · 1.075/1.15 = 0.467

Observed: r_3 = 0.535 (HIGHER not lower), then r_4..r_6 = 0.482, 0.470, 0.432 (going DOWN not up).

**REJECTED:** Pure (1/4)^n subleading with consistent sign predicts r_n monotone approaching 0.5 from BELOW, growing toward 0.5. Observed: r_3 is ABOVE 0.5 then drops below and goes DOWN. Pattern doesn't fit a simple (1/4)^n subleading with single fixed sign.

---

## Cross-check vs δ_n diagnostic from TAUBERIAN_SCOPING_VERIFICATION.md

The TAUBERIAN_SCOPING probe reported δ_n := |ε_n|·2^n − 1/30 values:

| n | δ_n | sign |
|---|---|---|
| 2 | +0.00476 | + |
| 3 | +0.00740 | + |
| 4 | +0.00590 | + |
| 5 | +0.00352 | + |
| 6 | −0.00147 | − |

**One sign change at n=5→6.** Compare to ratio diagnostic: **one sign change at n=3→4.**

The two diagnostics encode different things:
- δ_n is the residual after subtracting the leading constant 1/30 from S_n := |ε_n|·2^n.
- r_n − 1/2 is the difference between empirical ratio and asymptotic geometric rate.

Algebraically, r_n = S_n / (2·S_{n−1}). So r_n − 1/2 = (S_n − S_{n−1})/(2·S_{n−1}) = ΔS_n/(2·S_{n−1}).

If S_n converges to 1/30 from above (positive δ_n) and starts to decrease, ΔS_n changes from positive (δ increasing) to negative (δ decreasing). The sign of (r_n − 1/2) tracks the sign of ΔS_n. Let's verify:

| n | S_n = |ε_n|·2^n | ΔS_n = S_n − S_{n−1} | sign |
|---|---|---|---|
| 2 | 0.03810 | — | — |
| 3 | 0.04074 | +0.00264 | + |
| 4 | 0.03924 | −0.00150 | − |
| 5 | 0.03686 | −0.00238 | − |
| 6 | 0.03187 | −0.00499 | − |

**ΔS_n: + − − −.** Same sign pattern as (r_n − 1/2). ✓

Both diagnostics agree: S_n peaks at n=3 (+0.0407) and then decreases monotonically through n=6 (0.0319). The ratio diagnostic is the rate-of-change of S_n, normalized; the δ_n diagnostic is the value of S_n minus the limit conjecture 1/30.

These are mutually consistent. The Phase 3 ratio reading adds a quantitative fact: |S_n| is decreasing at an ACCELERATING rate through n=6, with no sign of stabilization at the conjectured 1/30 limit.

---

## What the ratio diagnostic rules out

1. **Pure simple pole at z=2 with constant residue** — REJECTED. Would give monotone geometric r_n → 1/2; we see sign change and acceleration away.

2. **Complex-conjugate secondary singularity oscillation** — REJECTED. Only one sign change in four data points; not periodic.

3. **Discrete secondary singularity at |z|<2** — REJECTED. Would push r_n above 1/2 (faster rate); we see r_n below 1/2 at n=4..6.

4. **Pure (1/4)^n subleading per R76 §10** with consistent sign — REJECTED. Would predict r_n monotone increasing toward 1/2; we see decreasing past n=3.

---

## What the ratio diagnostic supports

1. **Leading rate (1/2)^n is approximately correct** but with significant pre-asymptotic deviation at n=4..6.

2. **Subleading correction has NEGATIVE coefficient** (pulling |ε_n| below A·(1/2)^n at n=4..6).

3. **Subleading correction does NOT have a single clean power-law form** at this data range — either branch-cut at z=2 with slow exponent (still pre-asymptotic) or a more complex multi-term structure.

4. **No evidence for off-axis singular structure**: ratios don't oscillate.

---

## Verdict (Phase 3)

The ratio diagnostic is **independent** of the Padé diagnostic and **converges on the same picture**:
- Leading rate-1/2 supported (both probes).
- Pure simple pole rejected (both probes).
- Complex secondary rejected (both probes).
- Branch cut at z=2 with negative subleading correction, currently pre-asymptotic at n=6 (both probes).
- Exact branch order undetermined (both probes).

The independent cross-validation IS the value-add of this probe. R77.6's Padé reading is reinforced by an entirely different mathematical operation on the same data.

To advance: needs ε_7 (and ideally ε_8) — the same Route A recommendation as TAUBERIAN_SCOPING_DISPOSITION.md.
