# Probe G4 — the partner at L=4 (Track E). Lambda A10, within-block power iteration.

**Date:** 2026-07-18  Partner = ρ(M_tower) = dominant eigenvalue of the γ≠0 principal submatrix (G0).
Method: plain power iteration WITHIN the tower block (233,280 states, 225.3M nnz), pure CPU/scipy, two
independent random starts. Sanctioned instrument (within-block runner-up ~0.29, ratio ~0.87 nominal; observed
effective ratio ~0.97 → ~1156 iters to res 1e-10). No ARPACK, no shift-invert. Instance terminated. ≈$1.

## R1 — ρ(M_tower, 4)
**ρ_4 = 0.33349990132**  (start1 0.3334999013221966, start2 0.33349990132218055 — **agree to 1.6e-14**;
res 9.9e-11 / 9.6e-11; 1156 / 1155 iterations). Clean, 8+ significant digits, two-start plateau agreement.

## R2 — Δ_4 = 1/3 − ρ_4  (braid point 3; c₀(4)=1/3+(2/3)2⁻⁵⁴, correction 3.7e-17 < double eps ⟹ 1/3 is the reference)
**Δ_4 = −1.66568e-4  (sign NEGATIVE).**
| L | 2 | 3 | 4 |
|---|---|---|---|
| Δ = 1/3 − ρ | −2.911e-3 | +9.958e-5 | **−1.66568e-4** |
| ρ vs 1/3 | above | below | **above** |

**The partner BRAIDS around 1/3.** Sign sequence −, +, − : ρ crosses 1/3 **twice** (between L2–L3 and again
L3–L4). Magnitudes 2.911e-3 → 0.996e-4 → 1.666e-4 — **non-monotone** (shrinks hard, then grows back up).
The sign is reported as a finding (no pre-registered pick, standing stance).

## R3 — the 27⁻ᴸ shot (zero evidential weight until DERIVED; NO fit under any outcome)
Predicted |Δ_4| ≈ 2·27⁻⁴ = **3.763e-6**.  Measured |Δ_4| = **1.666e-4**.  **ratio meas/pred = 44.3×.**
**The 27⁻ᴸ pattern DIES at L=4**, at the cost of one number — exactly the pre-registered falsification.
Post-mortem (not a fit): 2·27⁻² = 2.74e-3 (vs 2.911e-3) and 2·27⁻³ = 1.02e-4 (vs 9.96e-5, near-perfect) had
looked like a law — but L4 shows the L3 near-hit was **the oscillation passing near a zero-crossing**, not a
geometric decay. The braid, not a monotone cascade, governs Δ.

## R4 — γ-level mass profile of the converged partner eigenvector (both starts identical to ~1e-13)
| v₃(γ) level | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| measured mass frac | 0.666559 | 0.223220 | 0.076751 | 0.033470 |
| pre-reg [2/3,2/9,2/27,1/27] | 0.666667 | 0.222222 | 0.074074 | 0.037037 |
Pre-registered SHAPE class (geometric cascade + truncation tail) **CONFIRMED**: dominant 2/3, 2/9, then two
truncated deeper levels. The profile has **extended one level deeper than L=3** (L3 anchor [0.67,0.22,0.11,0]
cut off at v₃=2; L4 now populates v₃=3 at 0.0335), consistent with truncation at v₃ = L−1. Quantitative
deviations from pure 2·3⁻⁽ʲ⁺¹⁾ at the deep levels are ~3e-3 (levels 2,3), reported as deviations.

## R5 — SKIPPED (reported, per spec)
The F2-4 coupling needs the c₀ **right** eigenvector r₀; at L=4 r₀ is near-degenerate with the partner in the
FULL operator (ρ_partner/c₀ ≈ 0.9997) → power iteration cannot separate the two (F2-4's own recorded caveat).
Faithful L=4 coupling needs r₀'s closed form or a deflation method beyond the sanctioned within-block power
iteration. Skipped; g_4 / defectiveness-meter deferred to a method that resolves the near-degenerate pair.

## Headline
Partner nailed to 11 digits at L=4. The two load-bearing readouts: **(1) Δ braids — sign −,+,− with two
crossings of 1/3, non-monotone magnitude**; **(2) the 27⁻ᴸ magnitude law is falsified at L=4 (44× miss).**
The partner's approach to 1/3 is an **oscillation around it, not a geometric descent** — and the L3 point that
looked law-abiding was the braid near a crossing. Three braid points (−2.911e-3, +9.958e-5, −1.66568e-4) now
sit on the table for Track D (Wilson's cap-sector coupling ledger) to derive; they are not fitted.

Probe `probes/l4_partner_g4.py`; profile dump `outputs/g4_profile_L4.json`; log `logs/l4_partner_g4_log.txt`.
