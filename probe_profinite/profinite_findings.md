# Profinite-Limit Analysis of π_∞ on Z_3

**Date:** 2026-05-06.  Levels analyzed: k = 5 … 12.  K=12 used as π_∞ proxy.

## Phase 1: Inverse-limit consistency

CONSISTENCY OK: max relative L∞ residual = 1.522e-15 (machine precision). π_k tower is a valid inverse system; π_∞ well-defined as the projective limit.

Per-level residuals ||P_k π_{k+1} - π_k|| (k → k+1):

| k | L1 | L∞ | TV | rel L∞ |
|---|---|---|---|---|
| 5 | 1.520e-16 | 1.388e-17 | 7.598e-17 | 3.072e-16 |
| 6 | 1.883e-16 | 6.939e-18 | 9.416e-17 | 3.061e-16 |
| 7 | 3.721e-16 | 6.939e-18 | 1.861e-16 | 6.107e-16 |
| 8 | 4.425e-16 | 6.939e-18 | 2.213e-16 | 1.219e-15 |
| 9 | 5.084e-16 | 4.337e-18 | 2.542e-16 | 1.522e-15 |
| 10 | 4.880e-16 | 1.084e-18 | 2.440e-16 | 7.606e-16 |
| 11 | 4.931e-16 | 8.674e-19 | 2.465e-16 | 1.216e-15 |

## Phase 2: Moments, entropy, S_∞

Coordinate: x = s / 3^k ∈ [0, 1) for s ∈ (Z/3^k)*.

| k | n | E[x] | Var[x] | H(π_k) | H_uniform | H - H_uni | ε_k |
|---|---|---|---|---|---|---|---|
| 5 | 162 | 0.5418010226 | 8.2749e-02 | 4.6072 | 5.0876 | -4.8040e-01 | +0.0000e+00 |
| 6 | 486 | 0.5226087662 | 8.3408e-02 | 5.6566 | 6.1862 | -5.2959e-01 | -4.9791e-04 |
| 7 | 1,458 | 0.5121235232 | 8.5758e-02 | 6.7151 | 7.2848 | -5.6972e-01 | -1.1752e-03 |
| 8 | 4,374 | 0.5005342298 | 8.3591e-02 | 7.7802 | 8.3834 | -6.0324e-01 | -7.4555e-04 |
| 9 | 13,122 | 0.5010586565 | 8.3312e-02 | 8.8504 | 9.4820 | -6.3163e-01 | -7.5203e-06 |
| 10 | 39,366 | 0.5002843739 | 8.3528e-02 | 9.9247 | 10.5807 | -6.5599e-01 | +7.2075e-04 |
| 11 | 118,098 | 0.5003590063 | 8.3356e-02 | 11.0022 | 11.6793 | -6.7704e-01 | +1.5020e-03 |
| 12 | 354,294 | 0.5000574513 | 8.3419e-02 | 12.0825 | 12.7779 | -6.9536e-01 | +2.2747e-03 |

**S_∞ at the c=7/45 character group:** approximated by S_K = 7/15 + ε_K = 0.4667 + +2.2747e-03 = 0.4689413804. 

|ε_12| = 2.2747e-03.  Pre-registered S_∞ = 7/15 = 0.4666666667 exactly. Whether ε_k → 0 or saturates is the open extrapolation question (separate analysis on the order-3 recurrence).

## Phase 3: Convergence rates of π_k → π_∞

Lifted norm: ||lift_k^12(π_k) - π_12||_p where lift sends π_k uniformly across each fiber of (Z/3^12)* → (Z/3^k)*.

| k | L1 | L2 | L∞ | TV |
|---|---|---|---|---|
| 5 | 4.7130e-01 | 2.4801e-03 | 3.3610e-04 | 2.3565e-01 |
| 6 | 4.1386e-01 | 2.2964e-03 | 3.2566e-04 | 2.0693e-01 |
| 7 | 3.6189e-01 | 2.0971e-03 | 3.1000e-04 | 1.8094e-01 |
| 8 | 3.1135e-01 | 1.8764e-03 | 2.8650e-04 | 1.5568e-01 |
| 9 | 2.6114e-01 | 1.6257e-03 | 2.5124e-04 | 1.3057e-01 |
| 10 | 2.0986e-01 | 1.3279e-03 | 1.9836e-04 | 1.0493e-01 |
| 11 | 1.5134e-01 | 9.3936e-04 | 1.4743e-04 | 7.5671e-02 |

Geometric decay fits y_k = A · ρ^k (OLS in log-space):

| Norm | ρ | A | R² | n |
|---|---|---|---|---|
| L1 | 0.833711 | 1.2523e+00 | 0.970272 | 7 |
| L2 | 0.858785 | 5.8305e-03 | 0.926689 | 7 |
| Linf | 0.877038 | 7.3037e-04 | 0.874451 | 7 |
| TV | 0.833711 | 6.2613e-01 | 0.970272 | 7 |

**Comparison to ρ_slow ≈ 0.83:** the order-3 recurrence on ε_2..ε_12 has dominant real root ≈ 0.83. 
Mean rate across norms: ρ = 0.8508  (differs from 0.83 by 2.5%).

**Outcome A** — slow rate is the genuine inverse-limit convergence rate, agreeing with ε_k fit to within 5% across multiple norms.

## Files

- `result_profinite_consistency.csv`
- `result_profinite_moments.csv`
- `result_convergence_rates.csv`
- `pi_infinity_cylinder_representation.npy`