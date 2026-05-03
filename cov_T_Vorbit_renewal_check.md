# Cov[T, V_orbit | band] within renewal framework — outcome (a) confirmed (Result 53)

**Status.** Decisive. Renewal decomposition (Results 47/49) and σ-identity
Cauchy correction (Result 32) are **mutually consistent**. The renewal model
gives the asymptotic per-step v average; the σ-identity correction Cov/E[T]
is the exact finite-T Cauchy gap between per-step and per-orbit averages.

## Two decompositions, one quantity

For E[V_orbit | band]:

**Form (1)** — σ-identity (Result 32, exact algebraic tautology):
> E[V_orbit | B] = E_band-per-step(B) − Cov[T, V_orbit | B] / E[T | B]

**Form (2)** — renewal model (Results 47/49, asymptotic in T):
> E[V_orbit | B] ≈ μ_d(B) + λ(B)·(E_V(B) − μ_d(B))

where:
- μ_d(B) = avg v at non-cylinder steps given band
- E_V(B) = avg v at cylinder visits given band
- λ(B) = K/T = visit rate given band

## Empirical verification at N=2³⁶, 500K orbits

| band | ⟨T⟩ | direct ⟨V_orbit⟩ | (1) σ-id | (2) renewal | (1)−direct | (2)−direct |
|------|---:|---:|---:|---:|---:|---:|
| 0–25 |  49.14 | 2.3286 | 2.3286 | **2.2930** | +0.000000 | **−0.03564** |
| 25–50 | 71.95 | 2.0711 | 2.0711 | 2.0679 | +0.000000 | −0.00320 |
| 50–75 | 91.70 | 1.9657 | 1.9657 | 1.9640 | +0.000000 | −0.00170 |
| 75–95 | 117.38 | 1.8833 | 1.8833 | 1.8814 | +0.000000 | −0.00197 |
| 95–100 | 158.14 | 1.8079 | 1.8079 | 1.8055 | +0.000000 | −0.00237 |

Form (1) matches direct ⟨V_orbit⟩ to **machine precision** (Cauchy tautology).

Form (2) matches within 0.0017–0.0357 — the gap is exactly **−Cov[T,V|B]/E[T|B]**:

| band | Cov[T,V\|B] | E[T\|B] | −Cov/E[T] (predicted gap) | actual gap |
|------|---:|---:|---:|---:|
| 0–25 | −1.7514 | 49.14 | **+0.03564** | +0.03564 |
| 25–50 | −0.2304 | 71.95 | +0.00320 | +0.00320 |
| 50–75 | −0.1558 | 91.70 | +0.00170 | +0.00170 |
| 75–95 | −0.2314 | 117.38 | +0.00197 | +0.00197 |
| 95–100 | −0.3754 | 158.14 | +0.00237 | +0.00237 |

**Predicted gap = actual gap to 5 decimals across all bands.**

## Structural interpretation

Form (2) IS the per-step v average pooled across orbit steps in band:

```
E_per_step(B) = (sum of v across all steps in band orbits) / (total steps)
              = [K·E_V + (T−K)·μ_d] / T (orbit-summed)
              = λ·E_V + (1−λ)·μ_d
              = μ_d + λ·(E_V − μ_d)   ← Form (2)
```

So **Form (2) ≡ E_band-per-step** by construction. The renewal-model
asymptotic prediction IS the per-step v average.

Form (1) corrects for the finite-T per-orbit vs per-step averaging difference.
The Cauchy identity:
```
E[V_orbit | B] = E[T·V_orbit] / E[T] − Cov[T, V_orbit | B] / E[T | B]
              = E_per_step(B) − Cov[T, V_orbit | B] / E[T | B]
              ≡ Form (1)
```

So **Form (1) − Form (2) = −Cov/E[T]** by tautology. Both decompositions
are EXACTLY consistent.

## Why q=0.125 has largest gap

q=0.125 (low-σ tail): ⟨T⟩ = 49 (shortest orbits), |Cov[T,V|B]| = 1.75 (largest).
The product Cov/E[T] = 0.036 is the largest finite-T correction. As ⟨T⟩ grows
(higher σ-bands), the correction shrinks: 0.0017 at q=0.625.

In the asymptotic limit T → ∞, Form (2) becomes exact. For finite T at any
band, Form (1)'s Cov/E[T] correction quantifies the deviation.

## Cross-check with Result 33

Result 33 established that ΔK_band U-shape (boundary correction for K_h) was
largely a Result 32 baseline artifact. With per-orbit E_V_orbit (correct
baseline) instead of per-step E_band:

ΔK_band shrunk from up to +8.96 (using per-step) to ±0.72 (using per-orbit).

This is consistent with the Form (1) vs Form (2) distinction here:
- Per-step average corresponds to renewal asymptotic (Form 2)
- Per-orbit average corresponds to σ-identity exact (Form 1)
- Result 33's "correction" of K_bulk used the per-orbit baseline

The picture is internally consistent across Results 32, 33, 47, 49, and now 53.

## Renewal-model parameters per band

| band | E_V (visit avg) | μ_d (non-cyl avg) | λ (visit rate) |
|------|---:|---:|---:|
| 0–25 | 6.199 | 1.901 | 0.0912 |
| 25–50 | 5.963 | 1.802 | 0.0639 |
| 50–75 | 5.841 | 1.748 | 0.0527 |
| 75–95 | 5.744 | 1.704 | 0.0440 |
| 95–100 | 5.671 | 1.659 | 0.0365 |

Trends:
- E_V at visits decreases with band (high σ → orbits visit r=53 mod 64 more, V=5 deterministic)
- μ_d at non-cylinder steps decreases (high-σ orbits traverse lower-v residues)
- λ decreases (high σ → fewer visits per step)

All band-conditional parameters consistent with Result 50/52 visit-measure
analysis.

## Verdict — outcome (a)

**The renewal decomposition is internally consistent with Result 32's σ-identity
closure.** Form (1) is the algebraic tautology; Form (2) is the renewal
asymptotic; their difference is the Cauchy correction Cov/E[T], verified to
5 decimals across all bands.

The bridge equation framework is structurally tight at this level. Constants
1, 2 closed; Constant 4 bulk closes via two equivalent decompositions
(σ-identity = exact Cauchy, renewal = asymptotic per-step) that match;
Constant 4 boundary closes within ±0.7 (Result 33); Constant 3 reduces to
the visit-measure on r=21 cylinder (Results 42, 50/52).

## For v3.6 / Chang correspondence

> The σ-identity correction term (Result 32) and the renewal decomposition
> (Results 47/49) provide two equivalent decompositions of E[V_orbit | band].
> Renewal gives the asymptotic per-step v average μ_d + λ·(E_V − μ_d). The
> σ-identity adds the finite-T Cauchy correction Cov[T, V_orbit | B] / E[T | B],
> verified to 5 decimals across all bands. For long-T bands (high σ), the
> correction is <0.003; for short-T bands (low σ), up to 0.036. The two
> frameworks agree exactly in the asymptotic limit and differ only by the
> finite-T Cauchy gap.

## Files

- `experiments/77_renewal_cov_check.py`
- `experiments_output/77_renewal_cov_check_log.txt`
- `experiments_output/77_renewal_cov_check.csv`

Compute: 1.1s (500K orbits, walking + analysis).
