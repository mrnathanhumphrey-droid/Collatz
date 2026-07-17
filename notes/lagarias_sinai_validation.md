# Result 68 — Lagarias-Sinai heuristic validation: outcome (γ), structural deviations of order 0.5%-25%

**Date:** 2026-05-03. Tests whether v_t = v_2(3m_t + 1) along Syracuse orbits
follows Geom(½) at high empirical precision. Sample: 250M+ v values across
4 orbit-scale points N ∈ {2²⁸, 2³⁰, 2³², 2³⁴}. Code: `lagarias_sinai_validation.py`.

**Verdict (γ):** the heuristic holds at marginal level (E[v] ≈ 2.0 within 0.3%,
P(v=1) ≈ 0.5 within 0.03%) but has STRUCTURAL deviations at specific j values
that **do not vanish with sample size**. Implications for downstream
closures:

| Closure | Required precision | Actual precision | Status |
|---|---|---|---|
| Tao K_h = 3/log(4/3) | E[v] = 2.000 exactly | E[v] = 1.995 (-0.5%) | ε ≈ 0.026 in K_h |
| R64 \|μ̂(1/3)\|² | P(v even) = 1/3 exactly | P(v even) = 0.337 (+1.0%) | propagates into closed form |
| R65 4^(-k) decay | Geom tail at all j | 0.5%-25% deviations at specific j | decay law not exact |

## Methodology

Walked Collatz orbits at N ∈ {2²⁸, 2³⁰, 2³², 2³⁴} with parallel numba walker.
Recorded v_t per Syracuse step, m_t mod 32 at step entry, total σ, log(start).
Sample sizes:

| N | orbits | v samples |
|---|---|---|
| 2²⁸ | 2.0M | 130.7M |
| 2³⁰ | 1.0M | 70.2M |
| 2³² | 0.5M | 37.5M |
| 2³⁴ | 0.2M | 16.0M |

## Direct empirical distribution

```
   j     2^-j     N=2^28      N=2^30      N=2^32      N=2^34   ratio @ 2^34
   1   0.5000    0.50005    0.50002    0.50012    0.50014    1.0003
   2   0.2500    0.24258    0.24319    0.24350    0.24387    0.9755
   3   0.1250    0.12440    0.12443    0.12445    0.12450    0.9960
   4   0.0625    0.07803    0.07696    0.07602    0.07521    1.2033  ← +20%
   5   0.0313    0.03009    0.03014    0.03027    0.03031    0.9699
   6   0.0156    0.01266    0.01286    0.01304    0.01318    0.8434  ← -16%
   7   0.0078    0.00615    0.00626    0.00636    0.00643    0.8235  ← -18%
   8   0.0039    0.00287    0.00295    0.00299    0.00309    0.7900
   9   0.0020    0.00138    0.00141    0.00145    0.00147    0.7525
  10   0.0010    0.00122    0.00120    0.00119    0.00118    1.2106  ← +21%
```

**Pattern:** systematic deviations at j=4 (+20%) and j=10 (+21%); deficits
at j=6,7,8,9 of -16% to -25%. Pattern persists and shifts only mildly across
N — these are structural, not noise.

**Chi-squared GoF at N=2³⁴:** χ²/dof = 4373 (massively rejects Geom(½)).
For comparison, χ²/dof ≈ 1 would indicate Geom(½) holds within sampling
fluctuation. The factor of 4373 means the deviations are ~66 σ each.

## v-parity test (R64 foundation)

```
  N         P̂(v even)    Geom(1/2)    deviation
  2^28      0.337558     0.333333     +0.004224
  2^30      0.337357     0.333333     +0.004024
  2^32      0.336956     0.333333     +0.003623
  2^34      0.336731     0.333333     +0.003398
```

P(v even) is consistently **above** 1/3 by ~1%. The deviation slowly
shrinks with N (0.42% → 0.34% over a factor of 64 in N), but does not
extrapolate cleanly to zero.

## Higher moments

```
  N        E[v]      Var[v]     E[(v-2)³]   excess_kurt
  Geom(½)  2.0000    2.0000     6.0000      6.0000
  2^28     1.9946    1.8728     5.0176      5.6404
  2^30     1.9948    1.8806     5.0818      5.7150
  2^32     1.9951    1.8885     5.1400      5.7680
  2^34     1.9955    1.8962     5.2019      5.8391
```

- E[v] is below 2.0 by 0.005 (-0.25%)
- Var[v] is below 2.0 by ~0.11 (-5.5% — biggest deviation)
- Skewness × σ³ deviates by -0.8 to -1.0 (-13% to -17%)
- Excess kurtosis below 6.0 by 0.16-0.36 (-3% to -6%)

**Variance deficit -0.11 is the largest moment deviation** — the heuristic
is "too dispersed" relative to the empirical reality.

## σ-band conditional REVEALS more deviation, not less

The naive expectation: high-σ orbits (long survivors) provide cleanest
Geom(½) statistics due to mixing. Empirical reality at N=2³²:

```
  band       n_v        P(v=1)     P(v=2)     P(v=3)     P(v=4)
  Geom(½)               0.5000     0.2500     0.1250     0.0625
  0-25      5.2M        0.4201     0.2383     0.1509     0.0996
  25-50     7.9M        0.4782     0.2432     0.1324     0.0831
  50-75    10.2M        0.5072     0.2455     0.1217     0.0744
  75-95    10.6M        0.5314     0.2449     0.1141     0.0669
  95-100    3.6M        0.5550     0.2434     0.1078     0.0583

  Per-band χ²/dof:  22980, 4994, 3869, 8159, 5595
```

**The bands have RADICALLY different v distributions.** Band 0-25 has
P(v=1) = 0.420 (5σ-band orbits with FAST descent — many high-v steps);
band 95-100 has P(v=1) = 0.555 (long-tail orbits with slow descent —
many low-v steps).

The unconditional P(v=1) ≈ 0.500 is an **arithmetic average across bands
that happens to land near the Geom(½) value**. The "Geom(½) holds for
the trajectory measure" framing was a coincidence of band-averaging,
not a structural truth.

## Convergence rate

Deviations |P̂(v=j) - 2^(-j)| at N ∈ {2²⁸, 2³⁰, 2³², 2³⁴} for j=2:
0.0074, 0.0068, 0.0065, 0.0061. Reduction across 4× more samples (more
than 8x N): just 18%. Deviation is NOT 1/√N stochastic — it asymptotes
to a structural floor.

```
  N         1/√(n_v)         avg|dev j=1..5|
  2^28      8.7e-05           4.95e-03    ← 57× larger than CLT prediction
  2^30      1.2e-04           4.59e-03    ← 39× larger
  2^32      1.6e-04           4.33e-03    ← 27× larger
  2^34      2.5e-04           4.09e-03    ← 16× larger
```

The deviations exceed 1/√N by 16-57×. Confirms structural origin.

## Per-residue determinism (Step 2 quick check)

P(v=1 | r mod 32) at N=2³²:

```
  r ≡ 1 mod 4 (r ∈ {1, 5, 9, 13, 17, 21, 25, 29}): P(v=1) = 0.000
  r ≡ 3 mod 4 (r ∈ {3, 7, 11, 15, 19, 23, 27, 31}): P(v=1) = 1.000
```

**v=1 is fully deterministic given m mod 4.** Half of all residues force
v=1; other half force v≥2. The "P(v=1) ≈ 0.5" marginal is exactly the
50/50 split across r mod 4 (matching uniform residue distribution mod 4),
NOT a stochastic Geom(½) draw. This is **arithmetic determinism**, not
randomness.

This pattern repeats at deeper levels: v=1,2 deterministic given m mod 8;
v=1,2,3 deterministic given m mod 16; etc. The v-distribution at m mod 2^k
is fully deterministic; the "randomness" is only the marginal over m mod 2^k.

## Why deviations exist

The trajectory measure on Z is NOT uniform mod 2^k for k ≥ 3 (R63/R65
established mod-3 asymmetry; an analogous mod-2^k asymmetry follows by
inverse-Syracuse predecessor structure mod 2^k).

Specifically: the trajectory measure has weight ratios at each m mod 2^k
that propagate from inverse-tree branching. Those non-uniform weights
project to specific v-deviations from Geom(½):

- m residues mod 16 force specific v=4: P(v=4 | m ≡ 5 mod 16, …) = 1
- Trajectory measure over-weights certain mod-16 residues
- Net effect: P(v=4) is excess by 20% over Geom(½)

This is the same mechanism producing R63's mod-3 mass asymmetry, just
at the v-distribution level.

## Implications for downstream closures

### Tao K_h = 3/log(4/3)

Bridge equation `s_mean(r) = α_det(r) + K_h log N + ε` verified at slope
1.000 ± 0.005. Sensitivity: dK_h / dE[v] ≈ −K_h / E[v] = −5.21.

Empirical E[v] = 1.9955 (deviation −0.0045). Implied K_h correction:
+0.024. So actual empirical "K_h" ≈ 10.452 vs Tao's 3/log(4/3) = 10.428.

**The 0.5% verification at slope 1.000 is at ~0.005 precision — entirely
consistent with the 0.0023-magnitude correction implied by E[v] deviation.**
Tao K_h holds AT THE PRECISION OF THE VERIFICATION but not exactly.

### R64 closed form for |μ̂(1/3)|²

R64's closed form uses P(v even) = 1/3 exactly. Empirical P(v even) =
0.337 (+1% deviation). The closed form produces 0.306; this propagates
linearly into the predicted resonance.

Empirical |μ̂(1/3)|² (R63) = 0.306. R64's path-counting closed form gives
0.306. The match is at the 4-decimal level by design — the path-counting
analytics use the SAME population mass-fractions that produce the empirical
value. So the "closure" is somewhat tautological; it's not a derivation
of (a, b, c) from first-principles Geom(½), it's a regrouping of empirical
mass-fractions through the closed-form formula.

### R65 conjectured 4^(-k) decay

The decay law |μ̂(a/3^k)|² ≈ 0.31 × 4^(-(k-1)) was conjectured from k=1,2,3
empirical values. Geom(½) tail predictions show 16-25% structural
deviations at specific j. So the 4^(-k) law is at best **approximate**;
exact form requires the structural mod-2^k asymmetry corrections.

## Verdict

| Outcome | Status |
|---|---|
| (α) heuristic holds to ε ≤ 0.001 | **REJECTED** — deviations 0.5%-25% |
| (β) heuristic holds approximately, structural corrections | **PARTIAL** — applies at marginal level (E[v], P(v=1)) but not at full distribution |
| (γ) heuristic fails in ways affecting downstream | **PRIMARY** — full distribution has structural deviations that do not vanish with N |

The clean version: **the Lagarias-Sinai heuristic v ~ Geom(½) is a marginal
approximation that holds for E[v] and P(v=1) within 1% but has structural
deviations of 5-25% at specific j due to the trajectory measure's
non-uniform distribution mod 2^k.** Geom(½) is not the underlying
distribution — it's the average over arithmetic-deterministic v values
across residue classes mod 2^k.

## What this means for the body of work

**Substantive content unchanged:**
- R63 |μ̂(1/3)|² = 0.306 stands (computed from empirical mass-fractions)
- R65 3-adic specificity stands (mechanism is exact, decay constants approximate)
- Tao K_h verification at slope 1.000 ± 0.005 stands (at this precision)

**Reframings recommended:**
- Tao K_h: "verified at precision 0.005 in slope, consistent with E[v] = 2.0 ± 0.005"
  → not "exact verification of Tao's heuristic"
- R64 closed form: "regrouping of empirical mass-fractions through path-counting"
  → not "derivation of (a, b, c) from first principles"
- R65 4^(-k) decay: "conjectured approximate decay" → not "exact decay law"
- Lagarias-Sinai heuristic: "approximate at marginal moments, structural
  corrections at full distribution"

**Important positive finding:** v IS arithmetic-deterministic given m mod 2^k.
The "randomness" is the trajectory measure's distribution over m mod 2^k.
The Lagarias-Sinai heuristic's TRUE form is:

> **v_t given m_t mod 2^k is exactly determined by arithmetic; the
>  Geom(½)-like marginal arises from the trajectory measure being
>  uniform-on-m mod 2 but non-uniform mod 2^k for k ≥ 3.**

This is a STRONGER, more precise statement than the original heuristic.
Pinning down the trajectory measure's distribution mod 2^k for k → ∞
would give exact closed forms for every v-distribution moment.

## Files

- `lagarias_sinai_validation.py` — script
- `lagarias_sinai_validation_log.txt` — full output log
- `v_distribution_by_N.csv` — P̂(v=j) at each N
- `v_higher_moments.csv` — moments table
- `v_conditional_distributions.csv` — r mod 32 + σ-band conditionals
- `convergence_rate.csv` — N-scaling of deviations

## Concrete next moves

1. **Trajectory measure mod 2^k profile**: compute m mod 2^k distribution
   on the trajectory measure for k ∈ {3, 4, 5, 6, ...}. This gives the
   exact corrections to Geom(½).
2. **Path-counting with corrected mass-fractions**: re-derive R64's closed
   form using empirical mod-2^k mass-fractions instead of Geom(½) assumption.
   This should give exact agreement with empirical |μ̂(1/3)|² to all decimals.
3. **K_h refined**: derive K_h from empirical E[v] (= 1.9955) instead of
   assumed E[v] = 2.0. This gives "empirical K_h" ≈ 10.452.
4. **Honesty in framework synthesis chapter**: flag every appearance of
   "Geom(½)" or "K_h = 3/log(4/3)" with the precision qualifier from this
   validation.
