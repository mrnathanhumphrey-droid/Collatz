# M3_APPROACH_C — numerical computation

**Date:** 2026-05-11. Phase 2C of M3 probe. Specifies the numerical computation of M_3 = sup_γ ‖R(z, T_3)‖_op for the literal R77.2 T_3, plus its k-dependence at higher modulus 3^k (k = 2..5) if T_3 generalizes.

## C.1 The script (written, not executed — harness denies python)

```python
# M3_approach_c.py — numerical sup_γ ‖R(z, T_3)‖_op
import numpy as np

def T3():
    """R77.2 companion matrix of recursion (7/8, −7/32, 1/64)."""
    return np.array([[7/8, -7/32, 1/64],
                     [1.0,  0.0,   0.0 ],
                     [0.0,  1.0,   0.0 ]], dtype=np.complex128)

def resolvent_norm_on_contour(T, center=0.5, radius=0.125, n_pts=1000):
    """Compute sup_{z on circle} ||(z*I - T)^{-1}||_op."""
    I = np.eye(T.shape[0], dtype=np.complex128)
    theta = np.linspace(0, 2*np.pi, n_pts, endpoint=False)
    zs = center + radius * np.exp(1j * theta)
    norms = []
    for z in zs:
        R = np.linalg.inv(z * I - T)
        # operator 2-norm = largest singular value
        norms.append(np.linalg.norm(R, ord=2))
    return max(norms), np.array(norms)

T = T3()
sup_R, all_norms = resolvent_norm_on_contour(T)
print(f"M_3 (numerical) = {sup_R:.4f}")
print(f"Mean ||R||      = {all_norms.mean():.4f}")
print(f"Min  ||R||      = {all_norms.min():.4f}")

# Spectrum verification
eigs = np.linalg.eigvals(T)
print(f"Eigenvalues: {sorted(eigs.real)}")
# Expected: [0.125, 0.25, 0.5]

# Condition number of eigenvector basis
eigs, V = np.linalg.eig(T)
kappa = np.linalg.cond(V, p=2)
print(f"kappa(V) (2-norm cond): {kappa:.4f}")
```

## C.2 Expected output (analytical prediction)

From Approach A:

- `M_3 ≤ 944` (κ(V)-based bound, factor-of-118 from spectral radius).
- `M_3 ≥ 8` (spectral radius).
- True value will sit somewhere in [8, 944].

Analytical estimate: for a non-normal 3×3 with this conditioning, the operator norm is typically `(spectral-radius) × √(κ(V))` or thereabouts — `8 × √118 ≈ 87` as a rough midpoint. R77.2's "sharper 800–1000" estimate uses crude bounds; the true numerical value is likely **in the 50–200 range**.

> **Anticipated M_3 ≈ 50–200** based on analytical prediction. R77.2's quoted 800–1000 is an upper bound, not a tight estimate.

## C.3 k-dependence (extension to higher modulus)

The user task mentions "construct T_3's matrix at modulus 3^k for k=2..5". This is **not directly defined in project artifacts** — R77.2's T_3 is a fixed 3×3 matrix (the companion of a recursion conjectured to hold for ε_n at all n), not an operator parameterized by modulus.

Two possible interpretations:

### C.3.1 Interpretation 1: T_3 = T_3 fixed, ε_n at higher k

R77.2's T_3 is the operator that *would* govern ε_n for all n ≥ 1 if the 3-mode recursion held. The matrix is fixed; only ε_n's are sampled at different n. M_3 doesn't depend on k under this interpretation.

### C.3.2 Interpretation 2: T_k = companion of (rank-k) recursion fit to ε_n at higher k

A natural generalization is to fit ε_n to a k-mode recursion with k+1 coefficients, building a k×k companion matrix. R77.3 §7 (Stage A.4) tested this for k=4 (4-mode `(1/2, 1/4, 1/8, 1/16)`), found 5–10% residual at n=5,6 — also FALSIFIED.

So **T_k for higher k is also a falsified ansatz**. The candidate matrix would exist, but its spectrum doesn't match ε_n's actual dynamics.

### C.3.3 What numerical M_k would tell us (if pursued)

If one were to compute the resolvent norm of these falsified higher-rank T_k operators on contours around their conjectured "1/2" eigenvalue, one would get a sequence M_3, M_4, M_5, ... that tracks the **mathematical** behavior of the companion-matrix family, but **says nothing about Nisoli closure for c=7/45** because none of these T_k describe ε_n.

The stabilization-vs-polynomial-vs-divergence question is mathematically interesting (companion matrices of recursions with closely-spaced eigenvalues have predictable conditioning), but operationally moot.

## C.4 Outcome

> **APPROACH_C_DEFERRED.** Numerical computation specified; not executed (harness denies python in this task per the operational setup). Anticipated value `M_3 ∈ [50, 200]` based on Approach A's bounds and the typical non-normal-matrix ratio between operator norm and spectral radius. R77.2's quoted "800–1000" is the upper-bound estimate, not the true value.
>
> Even if executed, the value is **operationally moot** because R77.3 falsified the underlying T_3 as a description of ε_n's dynamics.

## C.5 Honest scope

(1) The script as written would execute in milliseconds on any standard machine; the bound is mathematically clean. The blocker is python denial in this task, plus the deeper blocker that the result doesn't translate to Nisoli closure under R77.3.

(2) Numerical M_3 ≈ 50–200 would be a tighter upper bound for downstream closure-table parameterization, but doesn't change the order-of-magnitude takeaway: closure at K=6 needs `|K| < √q / M_3 ≈ 47/100 ≈ 0.5` (with tight M_3=100) vs empirical `|K_max| ≈ 16.6` — still 33× too large.

(3) Even closing the M_3 = 100 vs M_3 = 1000 gap doesn't fix the closure inequality; the bilinear bound `|K|` is the load-bearing problem at K=6, and the spectrum-of-T issue dominates the operator-norm question.
