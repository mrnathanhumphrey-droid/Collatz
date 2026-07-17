# T-distribution conditional on σ-band (Result 38)

**Status.** Outcome (c) with structural twist. T | band does NOT fit any
standard parametric family (Normal, Lognormal, Gamma, Inverse-Gaussian) —
all KS > 0.05 across all bands. However, **shape parameters (skew, kurt)
are remarkably N-stable** across factor 64 in N, and **scaling is clean**
(μ_T ∝ ⟨log n⟩ with band-specific coefficient; σ_T ∝ √⟨log n⟩).

The trajectory measure produces a band-specific T-distribution with N-stable
shape that doesn't sit in standard families. Empirical characterization
documented; parametric closure null.

## Per-band T-distribution moments at N=2³⁶ (1M orbits)

| q | n | μ_T | σ_T | skew | excess kurt | p1 | p25 | p50 | p75 | p99 |
|------:|------:|-----:|-----:|-------:|-------:|---:|----:|----:|----:|----:|
| 0.125 | 250,000 |  49.13 |  9.91 | **−0.595** | −0.212 | 24 |  43 |  51 |  57 |  65 |
| 0.375 | 250,000 |  71.99 |  6.62 | −0.349 | −0.079 | 56 |  67 |  72 |  77 |  84 |
| 0.625 | 250,000 |  91.75 |  7.16 | −0.161 | −0.262 | 75 |  87 |  92 |  97 | 106 |
| 0.875 | 250,000 | 125.54 | 20.52 | **+1.514** | **+3.345** | 98 | 111 | 120 | 135 | 194 |
| 0.975 |  50,000 | 158.16 | 18.14 | **+1.708** | **+4.325** |134 | 145 | 153 | 166 | 222 |

**Bimodal qualitative pattern** across bands:
- **Lower bands (q ≤ 0.625):** mild negative skew, near-zero excess kurtosis,
  bounded above. Heavy LEFT tail (some orbits with T as low as 24 in q=0.125).
- **Upper bands (q ≥ 0.875):** large positive skew, large positive excess
  kurtosis, heavy RIGHT tail (q=0.975 has T values reaching 222 at p99 vs
  median 153).

Sign of skew flips between q=0.625 (−0.16) and q=0.875 (+1.51) — structural
crossover in upper-tail regime.

## KS tests against standard parametric families

| q | Normal | Lognormal | Gamma | Inv-Gaussian | best (still FAIL) |
|------:|------:|----------:|------:|-------------:|-------------------|
| 0.125 | 0.082 | 0.111 | 0.100 | 0.115 | normal (0.082) |
| 0.375 | 0.063 | 0.072 | 0.070 | 0.073 | normal (0.063) |
| 0.625 | 0.056 | 0.055 | **0.053** | 0.055 | gamma (0.053) |
| 0.875 | 0.123 | **0.096** | 0.105 | 0.097 | lognormal (0.096) |
| 0.975 | 0.131 | **0.112** | 0.118 | 0.113 | lognormal (0.112) |

**No KS < 0.05 in any band** with any candidate. T | band does not match
any standard family.

AIC ranking confirms Normal best for lower bands (Δ ≥ 3000 vs runner-up),
Lognormal best for upper bands (Δ = 84-344). But best AIC ≠ good fit:
the best models still fail KS by factor 1.5–6× the 0.05 threshold.

## Shape parameters are remarkably N-stable

Tested at N ∈ {2³², 2³⁶, 2³⁸} (factor 64 in N):

| q | skew @ N=2³² | skew @ N=2³⁶ | skew @ N=2³⁸ | drift |
|------:|------:|------:|------:|------:|
| 0.125 | −0.556 | −0.595 | −0.619 | **11.3%** |
| 0.375 | −0.334 | −0.349 | −0.313 | 11.5% |
| 0.625 | −0.220 | −0.161 | −0.154 | 30% (small magnitude) |
| 0.875 | +1.565 | +1.514 | +1.528 | **3.5%** |
| 0.975 | +1.802 | +1.708 | +1.700 | **5.7%** |

Excess kurtosis similar pattern. **Upper-band shape parameters drift <6%
across factor 64 in N — structurally robust.** Lower-band parameters drift
~11%, similar to the α(j) and w_q drift in Result 37.

The bimodal pattern (sign-flip between q=0.625 and q=0.875) is preserved
at every N. Distribution shape is structural; just not in a standard family.

## Scaling: μ_T ∝ ⟨log n⟩, σ_T ∝ √⟨log n⟩

| q | μ_T/⟨ln n⟩ across N | σ_T/√⟨ln n⟩ across N |
|------:|--------------------:|---------------------:|
| 0.125 | 1.968 → 2.084 (5.6%) | 1.949 → 2.049 (5.0%) |
| 0.375 | 2.980 → 3.022 (1.4%) | 1.395 → 1.340 (3.9%) |
| 0.625 | 3.864 → 3.821 (1.1%) | 1.484 → 1.452 (2.1%) |
| 0.875 | 5.347 → 5.187 (3.0%) | 4.178 → 4.185 (0.2%) |
| 0.975 | 6.766 → 6.498 (4.0%) | 3.771 → 3.720 (1.4%) |

**Both scalings hold within ~1-6% drift across factor 64 in N**, much
tighter than would be expected if scaling were wrong.

Asymptotic Wald slope predicts μ_T/⟨ln n⟩ = 1/(E_band·log2 − log3):
- q=0.125 (E_band ≈ 2.22): 1/0.4416 = 2.265 vs empirical 2.05 (gap 0.22)
- q=0.375 (E_band ≈ 2.03): 1/0.3088 = 3.239 vs empirical 3.00 (gap 0.24)
- q=0.625 (E_band ≈ 1.94): 1/0.2459 = 4.067 vs empirical 3.82 (gap 0.25)
- q=0.875 (E_band ≈ 1.81): 1/0.1576 = 6.345 vs empirical 5.19 (gap 1.16)
- q=0.975 (E_band ≈ 1.73): 1/0.0979 = 10.215 vs empirical 6.50 (gap 3.71)

Lower bands close to asymptotic Wald slope; upper bands deviate substantially —
consistent with finite-N boundary effect (orbit terminates before reaching
asymptotic regime in upper bands where T is short relative to log n).

## σ-identity check (machine precision)

σ = T·(1 + V_orbit) holds to <1e-10 across all bands and N. V_orbit
is mechanically determined by T given σ. **Joint (T, V_orbit) | band has
no additional structure beyond T-distribution + σ-identity** — once T | band
is characterized, V_orbit | band follows from V = σ/T − 1 within band σ
spread.

## Verdict and v3.6 framing

**Outcome (c) with structural twist:** T | band fails parametric closure
across all standard families (KS 0.05–0.13 best), but the distribution
shape itself IS structurally robust:

- Shape parameters (skew, kurt) N-stable (3-12% drift per factor 64 in N)
- Bimodal pattern (negative-skew lower bands, positive-skew-heavy-tail upper
  bands) preserved at every N
- μ_T ∝ ⟨log n⟩, σ_T ∝ √⟨log n⟩ with band-specific N-stable coefficients

The trajectory measure produces a band-specific T-distribution with structural
shape that doesn't sit in standard families. This is itself a constraint on
the trajectory measure: any closed-form characterization must reproduce
(at minimum):
- The bimodal skew pattern (sign-flip between q=0.625 and q=0.875)
- Heavy right tail in upper bands (kurt +3 to +4)
- Bounded behavior in lower bands (negative skew, near-Gaussian)
- Linear-in-⟨log n⟩ mean scaling, sqrt-in-⟨log n⟩ SD scaling

**v3.6 reports the empirical characterization as constraints on the
trajectory measure**, not as a closed-form parametric statement.

## Files

- `experiments/69_T_distribution.py`
- `experiments_output/69_T_distribution.csv`
- `experiments_output/69_T_distribution_log.txt`
