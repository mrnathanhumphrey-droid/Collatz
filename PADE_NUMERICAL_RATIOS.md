# PADE_NUMERICAL_RATIOS — Phase 4: ratio diagnostic extended through n=13

**Date:** 2026-05-12. Wilson. Extends PADE_EXTENSION_RATIOS to numerical n=7..13.

## r_n = |ε_n|/|ε_{n-1}|

| n | ε_n | |ε_n| | r_n = \|ε_n\|/\|ε_{n-1}\| | r_n − 1/2 | r_n − 0.984 |
|---|---|---|---|---|---|
| 2 | +9.524e-3 | 9.524e-3 | 0.04762 | -0.4524 | -0.9364 |
| 3 | −5.092e-3 | 5.092e-3 | 0.5347 | +0.0347 | -0.4493 |
| 4 | −2.452e-3 | 2.452e-3 | 0.4816 | -0.0184 | -0.5024 |
| 5 | −1.152e-3 | 1.152e-3 | 0.4697 | -0.0303 | -0.5143 |
| 6 | −4.979e-4 | 4.979e-4 | 0.4323 | -0.0677 | -0.5517 |
| 7 | −1.175e-3 | 1.175e-3 | **2.3604** | +1.8604 | +1.3764 |
| 8 | −7.455e-4 | 7.455e-4 | 0.6344 | +0.1344 | -0.3496 |
| 9 | −7.520e-6 | 7.520e-6 | **0.01009** | -0.4899 | -0.9739 |
| 10 | +7.208e-4 | 7.208e-4 | **95.841** | +95.341 | +94.857 |
| 11 | +1.502e-3 | 1.502e-3 | 2.0839 | +1.5839 | +1.0999 |
| 12 | +2.275e-3 | 2.275e-3 | 1.5145 | +1.0145 | +0.5305 |
| 13 | +2.948e-3 | 2.948e-3 | 1.2969 | +0.7969 | +0.3129 |

## Reading

The ratio diagnostic is **violently non-monotone** across n=2..13. The PADE_EXTENSION_DISPOSITION's reading "ratios LEAVING 0.5 in absolute terms (departure 0.018 → 0.030 → 0.068)" at n=4..6 was correct but only captures the pre-transition phase.

Three distinct regimes:

### Regime A (n=3..6): apparent approach to 1/2, branch-cut-consistent

r_3..r_6 = 0.535, 0.482, 0.470, 0.432. Monotonically declining from above to below 0.5. Departure |r_n − 0.5| grows from 0.035 to 0.068. This is what R77.6 read as "branch-cut at z=2 with finite-N corrections."

### Regime B (n=7..9): catastrophic oscillation, no fixed point

r_7 = 2.360 (jump up), r_8 = 0.634 (jump down), r_9 = 0.0101 (catastrophic collapse — ε_9 is near-zero node).

Neither rate-1/2 nor rate-0.984 captures this. The ratio diagnostic is incompatible with a SINGLE-MODE asymptotic at n=7..9. This is the structural signature of **at least two competing singularities** (or one complex-conjugate pair giving oscillation).

### Regime C (n=10..13): post-zero-crossing growth, decelerating ratios

r_10 = 95.84 (huge — coming up from the near-zero ε_9), r_11 = 2.084, r_12 = 1.515, r_13 = 1.297.

The ratios are MONOTONE DECREASING from r_10 onward, approaching 1 from above (but with another sign-flip plausibly imminent before they converge).

|r_n − 1/2| at n=10..13: 95.34, 1.58, 1.01, 0.80. Still decreasing, but very far from 0.5.
|r_n − 0.984| at n=10..13: 94.86, 1.10, 0.53, 0.31. Decreasing FASTER, suggesting (with very rough extrapolation) ratios may approach rate ~1 from above, not 0.5 or 0.984.

## limsup |ε_n|^(1/n) diagnostic (Hadamard)

| n | \|ε_n\|^(1/n) | radius ρ = 1/\|ε_n\|^(1/n) |
|---|---|---|
| 2 | 0.0976 | 10.25 |
| 3 | 0.1716 | 5.83 |
| 4 | 0.2229 | 4.49 |
| 5 | 0.2585 | 3.87 |
| 6 | 0.2870 | 3.48 |
| 7 | 0.3676 | 2.72 |
| 8 | 0.3925 | 2.55 |
| 9 | 0.2787 | 3.59 (artifact of near-zero ε_9) |
| 10 | 0.4852 | 2.06 |
| 11 | 0.5536 | 1.81 |
| 12 | 0.6022 | 1.66 |
| 13 | 0.6388 | **1.57** |

|ε_n|^(1/n) is MONOTONE INCREASING from n=2..8, dips at n=9 (near-zero artifact), then resumes monotone INCREASING from n=10..13. The implied radius of convergence is monotone DECREASING (except for the n=9 dip).

**The Hadamard estimate at n=13 says the leading singularity is at ρ ≈ 1.57**, not at z=2. The trend is INWARD: from ρ=2.06 at n=10 to ρ=1.57 at n=13. The slow-mode prediction ρ=1.016 would require continued inward shrinking; at n=13 we are 50% above that target.

For the slow-mode prediction to be the asymptotic, |ε_n|^(1/n) must continue to increase from 0.64 toward 0.984. The current rate of change is ~0.035 per step (n=10 to 13). To reach 0.984 requires ~10 more steps. Suggests the slow-mode asymptotic emerges around n ≈ 25, well beyond available data.

## Sign pattern: ++−−−−−−−++++

Single zero-crossing between n=9 and n=10. The ε_9 near-zero is the crossing point.

A single zero-crossing across n=2..13 is consistent with:
- A single complex-conjugate pair giving cos(n θ + φ) modulation with period 2π/θ. The half-period (sign change) span is 8..9 steps (depending on where in the cycle we entered). This matches the period-9.2 prediction from STATE.md.
- NOT consistent with a single real positive singularity giving monotonic sign.
- NOT consistent with two real singularities (would not produce oscillation, only sign-stable decays).

**Sign pattern strongly supports H_COMPLEX_PATTERN over H_TWO_SINGULARITIES_VISIBLE with two REAL singularities.**

## Sub-rate fit: phase + envelope

If ε_n = A · ρ^n · cos(n θ + φ), then |ε_n|·ρ^{-n} oscillates with period 2π/θ. Trying ρ = 1.57 (Hadamard at n=13):

|ε_n|·1.57^n at n=2..13:
- n=2: 0.0095·2.46 = 0.0235
- n=10: 7.21e-4·83.0 = 0.0599
- n=13: 2.95e-3·1037 = 3.06

This is NOT oscillating with fixed envelope — it's still growing. So ρ=1.57 is NOT the asymptotic rate; the actual rate is < 1.57 (= radius < 1.57). The Hadamard estimate at finite n is an OVER-estimate of the true radius.

Trying ρ = 1.016 (slow-mode):
|ε_n|·1.016^n at n=13: 2.95e-3·1.227 = 3.62e-3
|ε_n|·1.016^n at n=10: 7.21e-4·1.173 = 8.45e-4
Envelope at n=10: 8.45e-4; at n=13: 3.62e-3. RATIO 4.3 across 3 steps. So 1.016 is too large a radius (envelope still growing).

Trying ρ at radius such that the n=10..13 envelope is approximately stable: solve for r where |ε_13|·r^{-13}/|ε_10|·r^{-10} = 1 → (|ε_13|/|ε_10|)·r^{-3} = 1 → r^3 = 2.95e-3/7.21e-4 = 4.09 → r = 1.60. So ρ ≈ 1/1.60 = 0.625, radius 1.60. Matches Hadamard at n=13. So in the WINDOW n=10..13, the dominant rate is ≈ 1/1.60.

This is consistent with the closest-singularity-at-z≈1.6 reading. **NOT consistent with slow-mode-at-1.016 prediction at this n.**

## Comparison to slow-mode model (STATE.md prediction)

STATE.md says the slow oscillating mode emerges at k=7+ with ρ ≈ 0.984 per k-step and period ≈ 9.2 in k-space. In our notation, this is the radius of convergence prediction (radius = 1/0.984 ≈ 1.016).

The data at n=10..13 supports oscillation (sign pattern, decelerating ratios) BUT at an envelope rate of ρ ≈ 0.625, NOT 0.984. The data is in a TRANSIENT regime BETWEEN the rate-1/2 (n=2..6) asymptotic and the eventual slow-mode asymptotic.

**Implication:** the period 9.2 might already be visible (single sign-crossing in 8-step window), but the envelope rate has not yet converged to the slow-mode value. Either:
- (a) STATE.md's ρ=0.984 is correct as TRUE asymptotic but n=13 is still pre-asymptotic.
- (b) The slow-mode rate is actually steeper (~0.625 to 0.7 at this finite n), and STATE.md's 0.984 is itself an under-estimate from a fitting window that included transient data.

The probe doesn't resolve this; either way, the slow-mode is the live structural object, NOT z=2.

## Cross-link to PADE_EXTENSION_DISPOSITION

PADE_EXTENSION_DISPOSITION concluded H_AMBIGUOUS within n=2..6, with H_PURE_SIMPLE_POLE and H_COMPLEX_SECONDARY both REJECTED. The n=7..13 extension makes H_COMPLEX_SECONDARY (or H_COMPLEX_PATTERN — the present-tense relabeling) the FAVORED hypothesis:

- The sign pattern with ONE zero-crossing in 12-step window matches cos(n θ + φ) with period ~9-10. ✓
- The envelope is decelerating non-monotonically. ✓ (rate-1/2 → spike → decay → spike → growth)
- Padé [3/3] develops complex pair (at small |z|, likely artifact — but is the structural signature).

The diagonal Padé at [4/4]/[5/5] should make this concrete — if those approximants produce CC pairs at |z| ≈ 1.5..1.7, then H_COMPLEX_PATTERN is structurally confirmed.
