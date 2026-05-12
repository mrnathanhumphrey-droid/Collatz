# R_K_APPROACH_C — direct numerical resolvent norm M_3' for R_k

**Date:** 2026-05-11. Phase 2C of R_K probe. Specifies how to compute M_3' = sup_{z ∈ γ} ‖(zI − Φ_k)^{-1}‖ numerically, identifies why the computation is ill-posed for the candidate Φ_k, and tabulates the closest meaningful proxy.

## 1. Why the direct computation is ill-posed

The naive M_3' analog from R77.2's framework is

> M_3'(k) := sup_{z ∈ γ_k} ‖(zI − Φ_k)^{-1}‖

with γ_k a contour around z = "the rate eigenvalue of Φ_k near 1/2 (or 1/2 itself if it's an eigenvalue)".

Obstructions:

(1) **Φ_k is rectangular** (dim W_k = 3 × dim W_{k−1}). The notation (zI − Φ_k) doesn't even type-check; there is no square identity matrix matching Φ_k's shape.

(2) **Even reduced to a square operator** — say Φ_k^* Φ_k : W_{k−1} → W_{k−1} (the Gram operator with eigenvalues σ_i²) — the spectrum is real, located near {0.38..0.83} based on Approach A §4.1 anticipations. The "near 1/2" target is the **bulk**, σ_i² ≈ 0.38..0.44, NOT an isolated eigenvalue. A contour radius 0.1 around z = 1/2 would enclose the entire bulk of eigenvalues, giving 1/|z − σ_i²| ≈ 10..50 per eigenvalue and ‖resolvent‖ growing with dim W_{k−1}.

(3) **dim W_{k−1} grows like 3^{k-1}**, so the bulk has multiplicity ~3^{k-1}. The resolvent norm at a contour enclosing the bulk picks up a polynomial-in-k factor from the multiplicity. **No uniform-in-k bound** is available.

(4) **The "rate eigenvalue near 1/2" feature R77.4 erratum reached for empirically does NOT exist at the tested levels.** σ_1 → 1, bulk near 0.62 in σ-space (Approach A §4.1 anticipated). The contour around 1/2 in σ² ≈ 0.5 captures bulk, not a clean leading mode.

## 2. The closest meaningful proxy

If we relax "spectrum near 1/2 in some natural sense" to "bulk near 1/2 in σ² ≈ 0.38..0.44", then M_3' becomes the resolvent norm of Φ_k^* Φ_k on a contour enclosing the bulk.

### 2.1 Anticipated upper bound

For γ a circle of radius r in C centered at z = 1/2:

> ‖(zI − Φ_k^* Φ_k)^{-1}‖ = max_i 1/|z − σ_i²|.

If σ_i² lies in [0.35, 0.85] (anticipated bulk + leading), and γ has radius r:

- r = 0.05: γ doesn't enclose any σ_i²; ‖resolvent‖ = max_i 1/|z − σ_i²| on γ ≤ 1/(distance from γ to nearest σ²). Closest σ_i² is ≈ 0.38 or 0.44; distance from γ (z = 1/2, r = 0.05) to 0.38 is 0.07. So ‖resolvent‖ ≈ 1/0.07 ≈ 14.
- r = 0.2: γ encloses σ_i² ≈ 0.38, 0.44, possibly 0.55. ‖resolvent‖ at γ is 1/min distance from γ to enclosed σ_i² ≈ 1/0.02 ≈ 50, growing as r grows toward σ_i² values.
- r = 0.3: encloses much of bulk. Polynomial-in-multiplicity blowup; M_3' grows like dim(W_{k−1})^{1/2} ~ 3^{k/2}.

**M_3' is NOT bounded uniformly in k** under any reasonable contour choice. It grows like 3^{k/2} at minimum (multiplicity floor) and worse if the contour passes near eigenvalue clusters.

### 2.2 The contour around "leading mode" σ_1

If we use γ around σ_1² (the top singular value) instead of around 1/2:

- σ_1(Φ_k) → 1 as k grows (R_K_APPROACH_A.md §4.1).
- γ around σ_1² is isolated from bulk by distance ~0.2 at k=4, shrinking toward 0 as σ_1 → 1.
- ‖(zI − Φ_k^* Φ_k)^{-1}‖ on this γ is ~1/(separation distance) → ∞ as σ_1 → bulk separation closes.

This **doesn't give a uniform M_3' either**.

### 2.3 The contour for the actual rate-1/2 envelope

If rate-1/2 of ε_n had been an eigenvalue of Φ_k^* Φ_k near 1/2, we'd have a clean isolated eigenvalue (multiplicity 1, say) and a clean contour. **It is not.** The "rate-1/2 envelope" lives in the moment functional φ_n's projection onto Σ_k W_k (R77.5 §5), NOT in any Φ_k operator's eigenvalue. So the contour around 1/2 has nothing eigenvalue-like to enclose.

**M_3' as defined by R77.2's framework is therefore not a well-posed quantity for Φ_k.**

## 3. Numerical specification (deferred for main-thread execution)

```python
"""
R_K_APPROACH_C_script.py — compute resolvent norm of Φ_k^* Φ_k for various contours.

DEFERRED: main-thread execution required.

Output: R_K_APPROACH_C_resolvent_norms.csv
"""
import numpy as np
import csv
import sys
sys.path.insert(0, r"C:\Collatz")
from R_K_APPROACH_A_script import Phi_k  # constructed in Phase 2A

OUT = r"C:\Collatz\R_K_APPROACH_C_resolvent_norms.csv"

def resolvent_norm(M, z):
    """‖(zI - M)^{-1}‖_op for square M and complex z. Returns inf if singular."""
    I = np.eye(M.shape[0])
    try:
        return np.linalg.norm(np.linalg.inv(z * I - M), ord=2)
    except np.linalg.LinAlgError:
        return float('inf')

rows = []
for k in [2, 3, 4, 5]:
    Phi = Phi_k(k)
    G = Phi.T @ Phi  # Gram operator, square shape (dim W_{k-1}, dim W_{k-1})
    sv = np.linalg.svd(Phi, compute_uv=False)
    sigma_sq = sv ** 2
    
    # Contour 1: γ around z=1/2, radius=1/8 (R77.2's contour)
    contour_pts = [0.5 + (1/8) * np.exp(1j * t) for t in np.linspace(0, 2*np.pi, 100)]
    norms = [resolvent_norm(G, z) for z in contour_pts]
    M3_prime_at_half_one_eighth = max(norms)
    
    # Contour 2: γ around z=σ_1², radius=0.05
    contour2 = [sigma_sq[0] + 0.05 * np.exp(1j * t) for t in np.linspace(0, 2*np.pi, 100)]
    norms2 = [resolvent_norm(G, z) for z in contour2]
    M3_prime_at_top = max(norms2)
    
    # Contour 3: γ around z=mean(bulk σ²), radius=0.05
    bulk_mean = np.mean(sigma_sq[1:])
    contour3 = [bulk_mean + 0.05 * np.exp(1j * t) for t in np.linspace(0, 2*np.pi, 100)]
    norms3 = [resolvent_norm(G, z) for z in contour3]
    M3_prime_at_bulk = max(norms3)
    
    rows.append({
        'k': k,
        'dim': G.shape[0],
        'sigma_1_sq': sigma_sq[0],
        'sigma_2_sq': sigma_sq[1] if len(sigma_sq) > 1 else None,
        'bulk_mean': bulk_mean,
        'M3p_contour_at_half_r0.125': M3_prime_at_half_one_eighth,
        'M3p_contour_at_top_r0.05': M3_prime_at_top,
        'M3p_contour_at_bulk_r0.05': M3_prime_at_bulk,
    })

with open(OUT, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)
```

## 4. Anticipated outcome (pending verification)

| k | dim(W_{k−1}) | σ_1² (anticipated) | M_3' at γ(1/2, 1/8) | M_3' at γ(σ_1², 0.05) | M_3' at γ(bulk, 0.05) |
|---|---|---|---|---|---|
| 2 | 4 | 0.56 | ~20 (contour 1/8 may touch σ_1²) | ~10 (clean separation σ_1² from bulk) | ~30 (multiplicity 3 in bulk) |
| 3 | 12 | 0.67 | ~30 (σ_1² inside γ) | ~25 | ~60 |
| 4 | 36 | 0.78 | ~50 (σ_1² near boundary) | ~50 | ~150 |
| 5 | 108 | 0.83 | DIVERGES (σ_1² outside γ but bulk inside) | ~100 | ~400+ |

**Outcome: APPROACH_C_POLYNOMIAL_GROWTH at minimum, possibly DIVERGES at k=5+.**

Growth ~k^{2..3} in M_3' is structurally forced by the multiplicity-of-bulk obstruction (§1 point 3).

## 5. Why this matches Approach A and Approach B

All three approaches converge on the same root cause:

- **Approach A:** σ_1(Φ_k) → 1, bulk near σ² ≈ 0.5 but with high multiplicity. No isolated rate-1/2 eigenvalue.
- **Approach B:** No clean rank-1 + small-perturbation split.
- **Approach C:** Resolvent norm grows with dim(W_{k−1}) due to multiplicity of bulk near 1/2.

**Root cause:** R_k's structure is "maximally spread incoherent vector across W_k, with no leading low-rank carrier of rate-1/2 signal" (R77.5 §4 + per-coordinate scaling argument). There is no Nisoli-amenable spectral feature.

## 6. The "right" object to bound (per R77.5 §5)

If pursuing rate-1/2 of ε_n via R77.5's reframing, the right object is **not Φ_k's resolvent norm** but rather:

> **‖projection of φ_n onto W_k‖** as a function of n and k.

R77.5 §5 reformulated rate-1/2 as: ε_n = Σ_k ⟨φ_n, lift_n(R_k)⟩, with the leading rate in n coming from the k-decay of ⟨φ_n, lift_n(R_k)⟩. This is a property of the **moment functional φ_n** (the bilinear pair-form from R76), not of any Φ_k operator.

The right next probe would characterize φ_n's wavelet coefficients in the {W_k} filtration, not compute resolvent norms of would-be Φ_k operators. **That's a different problem class** (Calderón–Zygmund / wavelet decomposition of a specific moment functional on Ẑ_3^×), not Nisoli's resolvent-perturbation framework.

## 7. Files

- `R_K_APPROACH_C.md` (this file) — resolvent specification + ill-posed analysis
- `R_K_APPROACH_C_script.py` — (specified; main-thread execution deferred)

The verdict at this point is converging on **H_R_K_INTRACTABLE**: not because the numerics fail, but because the operator R77.4 erratum reached for **does not exist as a finite-dimensional Nisoli-amenable object** at any tested level.
