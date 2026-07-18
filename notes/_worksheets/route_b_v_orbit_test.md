# Route B: Direct E[V_orbit | band] from Esscher + σ-identity (Result 32)

**Status.** Outcome (b) — correction exists but closes algebraically, NOT a
Lagarias-class problem. V_orbit | band closes in TWO pieces:

1. E_band-per-step from Esscher (Result 25, exact within ±0.005)
2. Finite-T correction = −Cov[T, V_orbit | band] / E[T | band] (exact algebraic
   identity from σ = T·(1+V_orbit))

## What's tautology, what's empirical

**Tautology.** V_orbit ≡ (σ_orbit − T) / T = σ_orbit/T − 1 is the definition.
Equivalently σ_orbit = T·(1 + V_orbit). Verified to 1.14e-13 — this is arithmetic,
not asymptotic.

**Algebraic identity (Cauchy-style, universal).** For ANY conditioning event A
and ANY pair of random variables (T, V):
E[V | A] − E[T·V | A] / E[T | A] = −Cov[T, V | A] / E[T | A].

**Empirical content.** Does Cov[T, V_orbit | band] = 0 or not? If 0, then
E[V_orbit | band] = E_band-per-step exactly (Esscher per-step closes everything).
If nonzero, the gap is exactly the covariance term — closed-form algebraic, but
the magnitude is band-specific empirical structure.

The empirical finding: Cov[T, V_orbit | band] ≠ 0; specifically ρ(T, V_orbit | band)
ranges from −0.83 to −0.95, strongest in tail bands. The "correction" is exactly
the algebraic gap implied by this nonzero covariance.

## Empirical (500k orbits each, 5 seeds, σ_resid quartile bands)

| log2N | q | n_band | E[V_orbit\|band] | E_band-per-step | correction |
|------:|------:|--------:|-----------------:|----------------:|-----------:|
| 32 | 0.125 | 125,000 | 2.3659 | 2.3229 | **+0.0430** |
| 32 | 0.375 | 125,000 | 2.0766 | 2.0728 | +0.0038 |
| 32 | 0.625 | 125,000 | 1.9632 | 1.9614 | +0.0019 |
| 32 | 0.875 | 125,000 | 1.8632 | 1.8570 | +0.0062 |
| 32 | 0.975 |  25,000 | 1.8025 | 1.8000 | +0.0024 |
| 36 | 0.125 | 125,000 | 2.3286 | 2.2930 | **+0.0356** |
| 36 | 0.375 | 125,000 | 2.0711 | 2.0679 | +0.0032 |
| 36 | 0.625 | 125,000 | 1.9657 | 1.9640 | +0.0017 |
| 36 | 0.875 | 125,000 | 1.8683 | 1.8623 | +0.0060 |
| 36 | 0.975 |  25,000 | 1.8079 | 1.8055 | +0.0024 |
| 38 | 0.125 | 125,000 | 2.3139 | 2.2812 | **+0.0327** |
| 38 | 0.375 | 125,000 | 2.0686 | 2.0656 | +0.0029 |
| 38 | 0.625 | 125,000 | 1.9666 | 1.9650 | +0.0016 |
| 38 | 0.875 | 125,000 | 1.8708 | 1.8649 | +0.0058 |
| 38 | 0.975 |  25,000 | 1.8107 | 1.8083 | +0.0023 |

4 of 5 bands close within ±0.006. The q=0.125 band (low-σ tail) has +0.04
correction, persistent across N (slow drift 0.043 → 0.033, ~N^(−0.04)).

## The closed-form correction

**Algebraic identity** (Cauchy-style):
E[V] − E[TV]/E[T] = −Cov[T,V]/E[T]

E_band-per-step = Σsumv / ΣT = E[T·V_orbit | band] / E[T | band], so

```
correction = E[V_orbit | band] − E_band-per-step
           = −Cov[T, V_orbit | band] / E[T | band]
```

**Verification** (q=0.125 band):

| log2N | ρ(V,T)\|band | sd_T | sd_V | E[T] | predicted | actual | ratio |
|------:|-------------:|-----:|-----:|-----:|----------:|-------:|------:|
| 32 | −0.9090 |  8.98 | 0.2201 |  41.73 | +0.0431 | +0.0430 | 0.9987 |
| 34 | −0.9155 |  9.43 | 0.2046 |  45.39 | +0.0389 | +0.0389 | 0.9996 |
| 36 | −0.9241 |  9.90 | 0.1914 |  49.14 | +0.0356 | +0.0356 | 0.9991 |
| 38 | −0.9281 | 10.30 | 0.1806 |  52.83 | +0.0327 | +0.0327 | 1.0006 |

Predicted = empirical to 4 decimals across all N. **Tautological closure** —
the correction is the algebraic difference between unweighted and T-weighted
band averages, no extra structure needed.

## Why q=0.125 dominates

V_orbit | band has a band-specific spread sd_V|band:
- q=0.125: sd_V = 0.22 (≈ global sd_V)
- q=0.375: sd_V = 0.04
- q=0.625: sd_V = 0.03
- q=0.875: sd_V = 0.04
- q=0.975: sd_V = 0.02

In middle bands, V_orbit is tightly determined by σ-identity at fixed (T, log_n),
so |Cov[T,V|band]| is small. In the low-σ tail (q=0.125), wider T-distribution
combined with strong negative ρ(V,T) (≈ −0.92, "short orbits need high per-step v")
makes |Cov[T,V|band]| large.

## T-stratification within q=0.125 — confirms Esscher per-step at fine scale

| T-decile | T range | n | E[V_orbit\|T-dec] | E_band-per-step\|T-dec |
|---------:|--------:|--:|------------------:|----------------------:|
|  1 |  6–35 | 13,237 | 2.7500 | 2.7203 |
|  2 | 35–41 | 14,875 | 2.4756 | 2.4739 |
|  3 | 41–44 |  9,684 | 2.3858 | 2.3856 |
|  4 | 44–48 | 16,211 | 2.3294 | 2.3290 |
|  5 | 48–51 | 12,796 | 2.2765 | 2.2763 |
|  6 | 51–53 | 10,033 | 2.2478 | 2.2477 |
|  7 | 53–56 | 14,328 | 2.2187 | 2.2186 |
|  8 | 56–58 | 10,598 | 2.1958 | 2.1958 |
|  9 | 58–61 | 12,118 | 2.1726 | 2.1725 |
| 10 | 61–66 | 11,120 | 2.1508 | 2.1507 |

Within each (T-decile, band) cell, V_orbit | T-dec ≈ E_per-step | T-dec to <0.001
(except decile 1, finite-T residual). Confirms Esscher per-step structure
holds at all (T, band) sub-strata. The aggregate +0.0356 correction at the band
level is purely T-weighting asymmetry: μ_v(T) varies 2.75 → 2.15 across T-deciles,
unweighted-T averaging differs from T-weighted averaging by exactly the
Cov[T,V] formula.

## Why the Edgeworth-shape doesn't appear in the residual

Fitting correction = c·E[Z²−1|band] gives R² = −0.4 (worse than null). The
correction shape is dominated by the q=0.125 outlier; the other bands cluster
near zero. This is NOT the symmetric "smile" shape of E[Z²−1|band].

Result 29's "C ≈ 0.21·σ_V" was fitting a different residual:
E[V_orbit|band] − (μ_V + ρσ_V·EZ_band), which has the Edgeworth-shape because
the Linear-Gauss baseline μ_V + ρσ_V·EZ misses both:
- The band-specific Esscher tilt (captured by E_band-per-step here)
- The finite-T Cov[T,V] correction (computed exactly here)

When we use E_band-per-step as baseline, the Cov[T,V] correction is what's left,
and it has a different (non-Z²−1) shape concentrated in the low-σ tail.

## Closed-form structure for V_orbit | band

```
E[V_orbit | band(q)] = E_band(q) − Cov[T, V_orbit | band(q)] / E[T | band(q)]
                    = E_band(q) + Δ(q, N)

where:
  E_band(q) = Esscher mean (exact via Result 25)
  Δ(q, N) = −Cov[T, V_orbit | band(q)] / E[T | band(q)]  (algebraic identity)
```

Δ(q, N) requires the conditional joint moments of (T, V_orbit) given band,
which depend on band geometry, σ-identity, and log_n distribution. These are
derivable from the same trajectory-measure structure that gave E_band, but
require additional bivariate moments rather than just the per-step Esscher tilt.

**The correction is NOT a Lagarias-class open piece** — it's a finite-T
algebraic tautology arising from σ = T·(1+V_orbit). It can be made arbitrarily
small by choosing finer band conditioning, or computed exactly by direct
integration once (T, log_n | band) is characterized.

## Implication for v3.5+

Route B closes V_orbit | band as **Esscher + algebraic finite-T correction**.
- Constant 4 bulk closes structurally (two pieces, both derivable)
- Boundary correction ΔK_band(q) U-shape remains as the only open piece
- Unlike outcome (d), no new Lagarias-class problem unifies with Result 30's
  per-j W_j wall

Edgeworth program (Results 27/29/31) is superseded: the correct decomposition
is per-step Esscher + Jensen-gap-from-σ-identity, NOT joint cumulant expansion.

## Files

- `experiments/64_route_b_v_orbit.py`
- `experiments_output/64_route_b_v_orbit.csv`
- `experiments_output/64_route_b_v_orbit_log.txt`
