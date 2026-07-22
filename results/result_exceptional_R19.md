# Probe R19 — typical or exceptional — **the decider: the equipartition/regularity route SURVIVES R18-D**

> **⚠️ CORRECTIONS (R20, Wilson).** Two framings below are retracted — the **result (R19-B: the fixed low-harmonic
> coefficients are typical) stands**, but: (1) **R19-A's "pinned to ⟨2⟩ / 2-adic resonance" is VACUOUS** — 2 is a
> primitive root mod 3^r, so ⟨2⟩ is the *entire* unit group and "the argmax is a power of 2" is content-free. The
> real content is R20-C: the additive max migrates to the **trivial character** at dist (2/3)^r. (2) **R19-D's
> "A-side dominant spike m*=3^{r−1}, ≈0.233, not 7/30" is VACATED** — R20-D proves A_r(3^{r−1}) = −S_r/2 exactly
> (→ −7/30), a definitional DC-split artifact, not a spike. R19-D's *conclusion* (the additive spike doesn't touch
> the b_r route) survives on R19-B + R20-A/C, not on the A-side argmax.

**Date:** 2026-07-21  Reuses R7/R9/R10; exact rationals on the A-side. Probe `probes/probe_exceptional_R19.py`.
Per-r shape statistics only (no cross-r rates, no periods). Disambiguates the **additive** spike (R18-D) from the
**A-side** coefficients A_r(m) = γ_r(τ_m) − γ_{r−1}(τ_m) = C_{r+1}(m)/3 (real, exact; the dlog-angular moments,
R12-B; A_r(0)=S_r). **The verdict: they are different spikes on orthogonal orbits, and the additive max-death does
not touch the coefficients b_r depends on.**

## R19-A — THE ARGMAX (measurement): the additive spike is **pinned to the ⟨2⟩ orbit** (group-coherent, angle-wandering)
argmax_{3∤ξ}|μ̂_r(ξ)| and its conjugate N−ξ*:

| r | ξ* | = | ξ*/N | v₃(ξ*) | ξ* ∈ ⟨2⟩? | conj N−ξ* ∈ ⟨2⟩? |
|---|---|---|---|---|---|---|
| 2 | 4 | 2² | 0.444 | 0 | ✓ | ✓ (5) |
| 3 | 19 | 2¹² mod 27 | 0.704 | 0 | ✓ | ✓ (8=2³) |
| 4 | 16 | 2⁴ | 0.198 | 0 | ✓ | ✓ (65) |
| 5 | 32 | 2⁵ | 0.132 | 0 | ✓ | ✓ (211) |
| 6 | 64 | 2⁶ | 0.088 | 0 | ✓ | ✓ (665) |
| 7 | 256 | 2⁸ | 0.117 | 0 | ✓ | ✓ (1931) |

**Every argmax is a power of 2 mod 3^r** (the conjugate pair too), always in the **finest stratum v₃=0**. The naive
mod-3^{r−1} congruence test is False (ξ*/N wanders: 0.44, 0.70, 0.20, 0.13, 0.09, 0.12), but that is because
2^a mod 3^r **equidistributes as an angle** — the *group* location is coherent even when the *angle* is not. So the
answer is neither "fixed frequency" nor "wandering": **the additive spike is locked to the ⟨2⟩ direction — the
shift generator 2^{−v} of the dynamics itself** (the R17 slow mode, now identified as the 2-adic resonance).

## R19-B — TYPICAL-OR-EXCEPTIONAL (measurement; **THE decider**): the fixed low-harmonic coefficients are **typical, not exceptional**
|A_r(m)|² / stratum-typical, typical = (‖δ‖²_A/r)/N_{v₃(m)}, N_j = 2·3^{r−1−j}:

| r | m=1 (v₃0) | m=2 (v₃0) | m=3 (v₃1) | m=4 (v₃0) | m=9 (v₃2) | m=27 (v₃3) |
|---|---|---|---|---|---|---|
| 3 | 0.019 | 0.007 | 0.676 | 0.894 | 1.359 | — |
| 4 | 0.015 | 0.002 | 0.031 | 0.256 | 0.103 | 1.874 |
| 5 | 0.038 | 0.002 | 1.783 | 0.005 | 0.220 | 0.002 |
| 6 | 0.105 | 1.142 | 3.291 | 0.652 | 4.799 | 0.975 |
| 7 | 0.183 | 1.055 | 0.090 | 0.131 | 5.885 | 2.170 |

**The fixed small-m ratios are O(1) and oscillating — NOT growing.** m=1 is systematically **sub-typical** (0.02–0.18,
depleted — the R13-B depletion at the trivial-character end); m=2,3,4 hover O(1) (≤3.3, oscillating, no trend);
m=9 reaches 4.8–5.9 at r=6,7 but non-monotonically (1.36, 0.10, 0.22, 4.8, 5.9). **The low-harmonic coefficients
that b_r = ⟨δ_r, Re w⟩ couples to (small |m|, because Re w is smooth) are typical-or-depleted, not exceptional.**
⟹ **the additive max-spike (R18-D) does not reproduce as a growing spike on these coefficients; the
equipartition/regularity route survives.**

## R19-C — WITHIN-STRATUM DISTRIBUTION (measurement, NO fit): the A-side **does** spike — on isolated members, not the small m
Within-stratum max/typical of |A_r(m)|², and top-3 mass fraction:

| r | j | #members | max/typ | top-3 frac |
|---|---|---|---|---|
| 5 | 0 | 162 | 12.5 | 0.212 |
| 5 | 1 | 54 | 7.5 | 0.338 |
| 5 | 2 | 18 | 4.0 | 0.573 |
| 6 | 0 | 486 | 21.4 | 0.113 |
| 6 | 1 | 162 | 9.7 | 0.160 |
| 6 | 2 | 54 | 9.0 | 0.398 |
| 7 | 0 | 1458 | **40.1** | 0.076 |
| 7 | 1 | 486 | 12.4 | 0.074 |
| 7 | 2 | 162 | 12.3 | 0.198 |

**max/typ grows within every stratum (j=0: 12.5 → 40.1), strongest in the finest** — so the additive spiking (R18-D)
**does** reproduce on the A-side. But the **top-3 fraction shrinks** (0.21 → 0.08 at j=0): the spike is one isolated
exceptional member, not a concentration in the low harmonics. **This mildly refines R18-B: the equipartition is a
stratum-average phenomenon; member-by-member there are exceptional spikes** — but (R19-D) they land on the
resonance orbits, not the fixed small m.

## R19-D — SPIKE LOCATION vs FIXED m: the two spikes are on **orthogonal orbits in opposite strata**
Additive spike ξ* ∈ ⟨2⟩ (v₃=0, finest; transported A-stratum = v₃(ξ*) = 0 by the U support law). A-side argmax over
all m:

| r | additive ξ* (v₃) | A-side m* | = | v₃(m*) | \|A_r(m*)\| | m* in {1,2,3,4,9,27}? |
|---|---|---|---|---|---|---|
| 3 | 19 (0) | 9 | 3² | 2 | 0.2308 | yes |
| 4 | 16 (0) | 27 | 3³ | 3 | 0.2321 | yes |
| 5 | 32 (0) | 81 | 3⁴ | 4 | 0.2328 | no |
| 6 | 64 (0) | 243 | 3⁵ | 5 | 0.2331 | no |
| 7 | 256 (0) | 729 | 3⁶ | 6 | 0.2327 | no |

**The A-side dominant spike is m* = 3^{r−1} — a power of 3, the coarsest order-3 mode (v₃=r−1)** — with
|A_r(3^{r−1})| ≈ **0.233 constant** (exact rationals −15685/67963, −71597824829728245/…, …; converging near 0.2327,
**not** 7/30=0.2333). So:
- **Additive spike:** ⟨2⟩ orbit, **finest** stratum (v₃=0), magnitude grows relative to typical (R18-D).
- **A-side spike:** ⟨3⟩ orbit, **coarsest** stratum (v₃=r−1), magnitude **constant** ≈0.233.
- **Different generators (2 vs 3), opposite strata (0 vs r−1).** The transported additive stratum v₃(ξ*)=0 never
  equals v₃(m*)=r−1. And m* is a fixed small element only at r=3,4; for r≥5 it has marched out of {1,2,3,4,9,27}.

**Verdict:** the additive max-spike lands **away** from the fixed small-m coefficients. Moreover both resonances
(⟨2⟩ at wandering angle, ⟨3⟩ at harmonic index 3^{r−1}) pair with **negligible Re w weight** — Re w is smooth, so
(Re w)^(m) is large only at low harmonic index |m|, exactly where R19-B shows the coefficients are typical. **So
R18-D's "max route dead" is an additive-⟨2⟩ statement that does not touch the equipartition/regularity route for
b_r. Two distinct obstructions, not one.**

## R19-E — max_ξ|μ̂_r(ξ)| — six points + ratios (NOT a rate)
max|μ̂_r| is algebraic, not rational; float, with the successive-ratio column:

| r | max\|μ̂_r\| | ratios |
|---|---|---|
| 2 | 0.37792363 | — |
| 3 | 0.25223684 | 0.66743 |
| 4 | 0.17699888 | 0.70172 |
| 5 | 0.12927357 | 0.73036 |
| 6 | 0.09610639 | 0.74343 |
| 7 | 0.07587004 | 0.78944 |

Six points reported; no exponent named, no extrapolation.

## Status
**R19: the decider answers in favor of the regularity route.** **A** — the additive spike is **pinned to the ⟨2⟩
orbit** (the shift generator; finest stratum; group-coherent, angle-wandering — the R17 slow mode identified as the
2-adic resonance). **B (decider)** — the fixed low-harmonic coefficients A_r(m), m∈{1,2,3,4}, are **typical or
depleted (O(1), oscillating, not growing)**; m=1 sub-typical (R13-B). **C** — within-stratum there **is** a growing
max/typ spike (40× at r=7, j=0) on isolated members, top-3 fraction shrinking — refines R18-B's equipartition to a
stratum-average statement. **D** — the additive spike (⟨2⟩, v₃=0, growing) and the A-side dominant spike
(m*=3^{r−1}∈⟨3⟩, v₃=r−1, constant ≈0.233) are on **orthogonal orbits in opposite strata**, both away from the
low-harmonic small m and both paired with negligible Re w weight. **E** — six-point ratio column, no rate named.

**Consequence for the crux (owed to the pen):** **R18-D's max-coefficient death does NOT kill the equipartition /
low-frequency-decay route.** The additive max-spike is the ⟨2⟩ resonance (the slow mode); the theorem-relevant
object b_r = ⟨δ_r, Re w⟩ couples only to the low-harmonic A_r(m), which R19-B shows are typical. So the owed estimate
remains **summable cancellation of the low-harmonic C-table coefficients A_r(m)** (R18 crux, R12-B) — now known to be
**uncontaminated by the additive spike**, which lives on a different orbit. No fitting; exact A-side rationals,
labeled numeric additive magnitudes; the decider reported as a decider.
