# Probe D-1 — the fork: H1 (finite Jordan cluster) vs H2 (essential-curve finite sections)

**Date:** 2026-07-18  Direct/from-data (L=2,3 spectra from D1-A dense; L=4 top modes from D1-C block-SpMV);
F4 crossing rider = small local power iteration. No fits. **Verdict: H2 favored** (essential-curve finite
sections condensing on a symbol curve through 1/3), with F1's binary count the sole H1-reading discriminator —
and it reads H1 only because the accumulation has not yet reached the small-δ balls at L≤4.

## F1 — mode census N(δ, L) = #{λ : |λ − 1/3| < δ}
| δ | L=2 | L=3 | L=4 |
|---|---|---|---|
| 0.05 | 1 | 1 | 1 |
| 0.02 | 1 | 1 | 1 |
| 0.01 | 0 | 1 | 1 |
**Binary count: BOUNDED (=1) at accessible (δ, L) → reads H1.** BUT the leading complex pair's distance to 1/3
**contracts geometrically**: 0.3899 (L2) → 0.2065 (L3) → 0.0763 (L4), ratios 0.530, 0.370 → toward 1/3. The
modes ARE accumulating onto 1/3 (the essential-curve signature); they simply have not entered the δ≤0.05 balls
by L=4 (extrapolating ~×0.4/level, the pair enters δ=0.05 near L=5). **Distance trend → H2.**

## F2 — phase arithmetic (leading complex pair, raw)
arg(pair): 1.47615 (L2) → 0.65630 (L3) → 0.23064 (L4). Contraction ratios **0.4446, 0.3514 → 1/3.** The phase
quantizes and contracts toward the 3^{−L} rate — the finite-section angle-quantization signature. **→ H2**
(H1 would give a fixed fan or irregular phases).

## F3 — doublet anatomy
The two adjacent complex pairs and their splitting:
- L=2: leading pairs |λ| = 0.2362, 0.2041 — split 0.032, **no tight doublet yet**.
- L=3: 0.237640+0.183030j, 0.234999+0.183155j — **splitting 2.64e-3**.
- L=4: 0.320423+0.075242j, 0.320223+0.075252j — **splitting 2.00e-4**.
Splitting **shrinks 2.64e-3 → 2.00e-4 (ratio 0.076) with L**. A stable broken degeneracy (H1) would hold; a
splitting that collapses with L is **adjacent finite-section modes → H2**.

## F4 — the crossing rider (decisive; c₀ vs partner at λ = 0.4, 0.5, 0.6)
Pre-registered SHAPE: c₀ ~ (1−λ)/(1+λ) + fold; the partner stays within fluctuation scale of 1/3. **Confirmed:**

| λ | c₀ = Σw² (L=3) | (1−λ)/(1+λ) | partner ρ(M_tower) (L=3) | \|ρ − 1/3\| |
|---|---|---|---|---|
| 0.4 | 0.428571 | 0.428571 | 0.333946 | 6e-4 |
| 0.5 | 0.333336 | 0.333333 | 0.333236 | 1e-4 |
| 0.6 | 0.250051 | 0.250000 | 0.333076 | 3e-4 |

**c₀ moves as (1−λ)/(1+λ)** (fold-corrections visible only at L=2, D=6: c₀=0.432 vs 0.4286); **the partner stays
pinned at 1/3 for every λ** (|ρ−1/3| ≤ 6e-4, fluctuation scale). So **λ=½ is exactly where c₀ crosses the
partner's fixed 1/3 line** — off-resonance c₀ departs while the partner sits at 1/3, so the two SEPARATE and the
coalescence is specific to λ=½. The partner is pinned at the **critically-weighted central value 1/3** (weight-
free, consistent with G-D0's D0.3 ρ(S)=1/3 ∀λ); c₀ is the moving crosser.

## Verdict → names D-2
**H2 (essential-curve finite sections).** Three trend discriminators (pair-distance contraction, phase
quantization → 3^{−L}, doublet splitting collapse) and the crossing rider all point to finite sections of an
infinite-tape object whose spectrum condenses onto a **symbol curve through 1/3**, with 1/3 as the curve's
critically-weighted central value (the partner) and the complex pairs accumulating onto it as L→∞. The only
H1-reading discriminator is F1's raw count, bounded solely because the accumulation is still approaching the
small-δ balls at L≤4 — not evidence of a finite cluster.

⟹ **D-2 = THE SYMBOL THEOREM**: compute the spectrum of the infinite-tape fiber symbol (D-FORM + shell
telescoping); target = a curve through 1/3 with the mean-field identity (1/27)[[5,4],[4,5]]-core value 1/3 as
its central point — "the tower is critically weighted at (q=3, λ=½)" as an exact statement. The λ-crossing (F4)
and the phase/doublet contraction (F2/F3) are the finite-truncation data D-2's symbol must reproduce.

Probe `probes/probe_trackD1_fork.py`.
