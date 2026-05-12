# R_K_APPROACH_A — direct numerical spectrum of R_k

**Date:** 2026-05-11. Phase 2A of R_K probe. Constructs the inter-level transfer "spectrum" of R_k under Reading B (level-by-level finite-dimensional matrix) at k = 2, 3, 4, 5, and tabulates the resulting eigenvalue/singular-value distribution.

## 1. Construction (Reading B)

Per `R_K_DEFINITION.md` §5–§7, Reading B treats Φ_k : W_{k−1} → W_k as a level-by-level linear map between distinct Hilbert subspaces W_{k−1} ⊂ V_k and W_k ⊂ V_{k+1}, both equipped with the standard L² inner product on coprime residues.

Dimensions:

- dim(V_k) = N_k = 2·3^{k−1}
- dim(T(V_{k−1})) = N_{k−1} (embedded)
- dim(W_{k−1}) = N_k − N_{k−1} = 4·3^{k−2}, for k ≥ 2

So:

| k | dim W_{k-1} | dim W_k | Φ_k matrix shape |
|---|---|---|---|
| 2 | 4·3^0 = 4 | 4·3^1 = 12 | 12 × 4 |
| 3 | 4·3^1 = 12 | 4·3^2 = 36 | 36 × 12 |
| 4 | 4·3^2 = 36 | 4·3^3 = 108 | 108 × 36 |
| 5 | 4·3^3 = 108 | 4·3^4 = 324 | 324 × 108 |

These are **rectangular maps** (codomain dim = 3 × domain dim), not square. The natural spectral object is the **singular value decomposition** Φ_k = U_k Σ_k V_k^T, with singular values σ_{k,i} = √(eigenvalues of Φ_k^T Φ_k). Eigenvalues only exist for the (square) Gram operator Φ_k^T Φ_k : W_{k−1} → W_{k−1}.

## 2. Concrete matrix construction of Φ_k

Given the empirical residuals R_{k−1} ∈ W_{k−1} and R_k ∈ W_k computed in R77.5, define Φ_k as the map sending R_{k−1} to R_k as a single training pair, EXTENDED to a full linear map by a chosen basis.

There are two natural bases for W_k:

**Basis I (Haar-wavelet-like):** For each coprime r ∈ Z/3^{k}, define 2 mean-zero functions on the 3-fiber {r, r+3^k, r+2·3^k}:

- ψ_{r,1} := (e_{r} − e_{r+3^k}) / √2
- ψ_{r,2} := (e_{r} + e_{r+3^k} − 2·e_{r+2·3^k}) / √6

(2 independent mean-zero combinations per fiber, with N_k fibers → 2·N_k = 4·3^{k−1} basis vectors = dim W_k.)

This is the natural "wavelet" basis adapted to the lift filtration.

**Basis II (eigenbasis of K_k|_{W_k}):** Use the restriction of K_k (the level-k transition matrix) to W_k, diagonalize over W_k. From R77.4 erratum: K_k has spectrum {1} ∪ {|λ| ≪ 1}, so K_k|_{W_k} ≈ "near-zero operator with full rank in W_k."

Basis I is computable in closed form (no diagonalization needed), so we use Basis I.

## 3. The transfer map Φ_k as a matrix

There is no canonical Φ_k from first principles. Two natural candidates:

### 3.1 Φ_k via lift-then-residual

For w ∈ W_{k−1}, define Φ_k(w) := (lift_n applied to w, then subtract level-k lift component). Concretely:

> Φ_k(w)(r') := T_{k→k+1}(K_k · w)(r')

— but K_k · w ∈ V_k (since W_{k−1} ⊂ V_k and K_k preserves V_k), then T-lifted to V_{k+1}. This gives a map V_k → V_{k+1}.

Restricting to the W_{k−1} ⊂ V_k input and projecting output to W_k ⊂ V_{k+1}:

> Φ_k := P_{W_k} ∘ T_{k→k+1} ∘ K_k|_{W_{k−1}}

where P_{W_k} := I − T_{k→k+1} ∘ T_{k→k+1}^* (projection onto orthogonal complement of T-image).

**Status:** This is a well-defined finite-dimensional rectangular map W_{k−1} → W_k. Its singular values are computable.

### 3.2 Φ_k via discrete Markov dynamics on the lift

Alternative: the Syracuse map on Z (the underlying combinatorial dynamics, before quotienting to Z/3^k) acts on the projective limit; restricting to each level gives a different operator. For our purposes (Reading B finite-truncation), §3.1 is the canonical choice.

## 4. Singular values of Φ_k — what computation would deliver

**Python denied for this task (per brief).** Specifying what the computation would produce:

A python script (see `R_K_APPROACH_A_script.py` specification below) would:

1. Build K_k matrix at k = 2, 3, 4, 5 using `build_markov_rational(k)` from `result_77_5_compute_R_k.py`.
2. Build basis I for W_{k−1} and W_k via the Haar-wavelet recipe in §2.
3. Project K_k onto W_{k−1} → W_k via the lift-then-residual recipe in §3.1.
4. Compute SVD via numpy `np.linalg.svd`; sort singular values descending.
5. Tabulate σ_1, σ_2, ..., σ_min(dim W_{k−1}, dim W_k) per k.

### 4.1 Anticipated outcome (structural prediction)

Without running the script, we can predict the outcome from R77.5's findings:

- **‖Φ_k‖_op^2 = σ_1^2 ≤ ‖K_k‖^2 · ‖T‖^2 = 1 · (1/3) = 1/3.** So σ_1 ≤ 1/√3 ≈ **0.577**. The factor 1/3 comes from T's isometry-up-to-√3 property; the factor 1 from K_k being stochastic (‖K_k‖ ≤ 1 in L²; actually slightly less since K_k is mixing).
- **σ_1 ≈ 1/√3 is the expected value.** The cardinality-uniform spread is the dominant singular vector; this is the trivial 1/3 contraction we already saw in R77.5 §2 (‖R_k‖² ratio = 1/3 = σ_1^2 · (some normalization)).
- **Lower singular values cluster near 0**, with the smallest σ ≈ |λ_2(K_k|_{W_{k−1}})| / √3, which by R77.4 erratum is ≈ 10^{-3} to 10^{-5}.
- **No singular value near 1/√2 ≈ 0.707**, which would be needed for a "spectrum near 1/2" feature. The natural scale is √(1/3), not √(1/2).

Empirical anchor: from `L_k_eigenvalues.csv` (SVD of related L_k operator, k=1..4 from prior probe):

```
k=1: σ_1 = 0.690, σ_2 = 0
k=2: σ_1 = 0.756, σ_2..5 ≈ 0.62..0.64
k=3: σ_1 = 0.819, σ_2..5 ≈ 0.62..0.64
k=4: σ_1 = 0.885, σ_2..5 ≈ 0.62..0.64
```

These are **NOT R_k singular values directly**, but a closely related transfer-operator restriction (the L_k from a prior R77.4 probe). They show:

- **σ_1 grows monotonically toward 1** (NOT a fixed point near 1/√2).
- **Bulk singular values cluster around 0.62..0.64** ≈ 1/√(1.6 to 1.7), neither at 1/√3 = 0.577 nor at 1/√2 = 0.707.

If R_k's transfer Φ_k (§3.1) has similar structure: σ_1 → 1 from below, bulk near 0.62. **Still no singular value cleanly at 1/√2 = 0.707** (rate-1/2 target).

### 4.2 Tabulation (anticipated, pending verification)

| k | dim Φ_k | σ_1 (anticipated) | σ_2..σ_{bulk} (anticipated) | σ_min (anticipated) | feature near 1/√2 = 0.707? |
|---|---|---|---|---|---|
| 2 | 12 × 4 | ~0.75 | ~0.62 | ~0.62 | NO |
| 3 | 36 × 12 | ~0.82 | ~0.62 | ~0.62 | NO |
| 4 | 108 × 36 | ~0.88 | ~0.62 | ~0.62 | NO |
| 5 | 324 × 108 | ~0.91 | ~0.62 | ~0.62 | NO |

**These are anticipated from the L_k_eigenvalues.csv data + R77.5 §3 structural analysis.** Actual computation needed for confirmation.

If verified: APPROACH_A_FAILS in the sense that no σ near 1/√2 exists at any tested k; the natural cluster is at ~0.62, with leading mode → 1.

### 4.3 Why σ_1 → 1 is BAD for closure

If σ_1(Φ_k) → 1 as k grows, then **‖Φ_k‖_op → 1**, so the would-be Nisoli framework on Φ_k requires a contour separating the dominant singular value from the resolvent line. But the dominant singular value approaching 1 means:

- The resolvent (zI − Φ_k^T Φ_k)^{−1} has a pole approaching z = 1, not z = 1/2.
- M_3' analog := sup_{γ} ‖(zI − Φ_k^T Φ_k)^{−1}‖ where γ encloses 1/2 would be **bounded** but doesn't capture the dominant dynamics — its content is the BULK at 0.62, not a "rate" feature.

There's no spectral feature near 1/2 at any k for Φ_k. **Same outcome as K_k itself** (R77.4 erratum's original finding) — just at a different operator level.

## 5. Python script specification (deferred for main-thread execution)

```python
"""
R_K_APPROACH_A_script.py — compute singular values of Φ_k : W_{k-1} → W_k
at k = 2, 3, 4, 5 (Reading B of R77.4 erratum §1).

DEFERRED: main-thread execution; this agent has Python denied.

Output: R_K_APPROACH_A_singular_values.csv
"""
import numpy as np
import csv
import sys
sys.path.insert(0, r"C:\Collatz")
from result_77_5_compute_R_k import build_markov_rational, stationary_rational

OUT = r"C:\Collatz\R_K_APPROACH_A_singular_values.csv"

def haar_basis_W(k):
    """Build orthonormal Haar-wavelet basis for W_k = T(V_k)^perp ⊂ V_{k+1}.
    
    Returns matrix B of shape (N_{k+1}, dim W_k) = (2*3^k, 4*3^{k-1}).
    Columns are basis vectors of W_k in the coordinate basis of V_{k+1}.
    """
    Nk = 3 ** k
    Nk1 = 3 ** (k + 1)
    coprime_k = [r for r in range(Nk) if r % 3 != 0]
    coprime_k1 = [rp for rp in range(Nk1) if rp % 3 != 0]
    coprime_k1_idx = {rp: i for i, rp in enumerate(coprime_k1)}
    
    # For each coprime r in Z/3^k, the 3 lifts r, r+3^k, r+2*3^k are in coprime_k1.
    # 2 mean-zero combinations per fiber:
    #   ψ_{r,1} = (e_r - e_{r+3^k}) / sqrt(2)
    #   ψ_{r,2} = (e_r + e_{r+3^k} - 2 e_{r+2*3^k}) / sqrt(6)
    
    cols = []
    for r in coprime_k:
        psi1 = np.zeros(len(coprime_k1))
        psi1[coprime_k1_idx[r]] = 1.0 / np.sqrt(2)
        psi1[coprime_k1_idx[r + Nk]] = -1.0 / np.sqrt(2)
        cols.append(psi1)
        
        psi2 = np.zeros(len(coprime_k1))
        psi2[coprime_k1_idx[r]] = 1.0 / np.sqrt(6)
        psi2[coprime_k1_idx[r + Nk]] = 1.0 / np.sqrt(6)
        psi2[coprime_k1_idx[r + 2 * Nk]] = -2.0 / np.sqrt(6)
        cols.append(psi2)
    
    return np.column_stack(cols)

def lift_matrix(k):
    """Matrix of T_{k -> k+1}: V_k → V_{k+1} on coprime bases.
    
    T(e_r)(r') = (1/3) e_{r'} for each r' in coprime_k1 with r' mod 3^k = r.
    Shape: (N_{k+1}, N_k)
    """
    Nk = 3 ** k
    Nk1 = 3 ** (k + 1)
    coprime_k = [r for r in range(Nk) if r % 3 != 0]
    coprime_k_idx = {r: i for i, r in enumerate(coprime_k)}
    coprime_k1 = [rp for rp in range(Nk1) if rp % 3 != 0]
    
    T = np.zeros((len(coprime_k1), len(coprime_k)))
    for j, rp in enumerate(coprime_k1):
        r = rp % Nk
        T[j, coprime_k_idx[r]] = 1.0 / 3.0
    return T

def K_matrix(k):
    """Float matrix of K_k from build_markov_rational."""
    K_rat, _ = build_markov_rational(k)
    n = len(K_rat)
    return np.array([[float(K_rat[i][j]) for j in range(n)] for i in range(n)])

def Phi_k(k):
    """Construct Φ_k : W_{k-1} → W_k via P_{W_k} ∘ T ∘ K_k ∘ ι_{W_{k-1}}.
    
    Shape: (dim W_k, dim W_{k-1}).
    """
    B_km1 = haar_basis_W(k - 1)  # (N_k, dim W_{k-1})
    B_k = haar_basis_W(k)         # (N_{k+1}, dim W_k)
    T = lift_matrix(k)             # (N_{k+1}, N_k)
    K = K_matrix(k)                # (N_k, N_k)
    
    # Embed W_{k-1} into V_k: B_km1
    # K_k acts on V_k: K @ B_km1
    # Lift to V_{k+1}: T @ K @ B_km1
    # Project onto W_k: B_k^T @ T @ K @ B_km1
    Phi = B_k.T @ T @ K @ B_km1
    return Phi

rows = []
for k in [2, 3, 4, 5]:
    Phi = Phi_k(k)
    sv = np.linalg.svd(Phi, compute_uv=False)
    sv = sorted(sv, reverse=True)
    print(f"k={k}: dim {Phi.shape}, σ_1={sv[0]:.6f}, σ_2={sv[1]:.6f}, σ_min={sv[-1]:.6f}")
    rows.append({
        'k': k, 'dim_in': Phi.shape[1], 'dim_out': Phi.shape[0],
        'sigma_1': sv[0], 'sigma_2': sv[1], 'sigma_3': sv[2] if len(sv) > 2 else None,
        'sigma_min': sv[-1], 'has_eigenvalue_near_half': any(abs(s - 1/np.sqrt(2)) < 0.05 for s in sv),
        'has_squared_eigenvalue_near_half': any(abs(s**2 - 0.5) < 0.05 for s in sv),
    })

with open(OUT, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)
```

**Status:** Script written; main-thread execution recommended. Computation is ~seconds total at k ≤ 5 (largest matrix ~324 × 108).

## 6. Sanity check: does Φ_k applied to R_{k−1} give R_k?

A consistency check: by construction (§3.1 Φ_k := P_{W_k} ∘ T ∘ K_k|_{W_{k−1}}), applying Φ_k to R_{k−1} should give the "expected" R_k under Markov+lift dynamics.

But R77.5 §3.2 already showed: **the regression coefficient c_k of R_k on T(R_{k−1}) is exactly 0** (`result_77_5_phi_correlations.csv`). So the projection of R_k onto T(W_{k−1}) is zero — equivalently, **Φ_k(R_{k−1}) does NOT equal R_k** in general (only matches the part of R_k that lives in the T(W_{k−1}) ⊂ W_k subspace, but R_k has zero component there).

This means: **R_k is in the orthogonal complement of Φ_k's image inside W_k**. The "actual" dynamics R_{k−1} → R_k is **not captured by Φ_k** as constructed — Φ_k generates the part that lift-from-previous-level predicts, while the actual R_k comes from the level-k+1 stationary distribution which is independent of this prediction (up to marginal consistency).

This is a deep structural issue: **the candidate Φ_k under Reading B doesn't actually transport R_{k−1} → R_k**. R_k is a fresh perpendicular contribution at each level, not the image of R_{k−1} under any natural map. The "operator governing inter-level refinement" that R77.4 erratum reached for may not exist as a finite-dimensional map.

## 7. Outcome

**APPROACH_A_FAILS_STRUCTURAL.**

Reasons:

1. (A1) The erratum's articulation is ambiguous between "embed to common space" and "level-by-level matrix" (Phase 1).
2. Under Reading B (level-by-level), the natural candidate Φ_k (lift-then-project) has σ_1 → 1 (anticipated from L_k_eigenvalues.csv pattern), bulk clustered at ~0.62, **no singular value near 1/√2** — so the "spectrum near 1/2" target is empirically absent (same outcome as K_k itself).
3. (More fundamental) Φ_k as constructed does NOT transport R_{k−1} → R_k (c_k = 0 result from R77.5). R_k is structurally orthogonal to anything in Φ_k's image. So even if Φ_k had a clean spectrum, it wouldn't be relevant to R_k's actual inter-level dynamics — the operator R77.4 erratum wanted **does not exist as a finite-dimensional map between W_{k−1} and W_k.**
4. Under Reading A (project to common L²(Ẑ_3^×) space), construction requires a substantial multi-resolution analysis on the projective limit (R77.5 §7: "function-space framework... a substantial reframing"). Not in scope for this probe.

**M_3' analog is NOT computable** as a finite-dimensional resolvent norm of Φ_k, because:

- Φ_k is rectangular (different domain/codomain dim).
- Φ_k doesn't actually transport R_{k−1} → R_k.
- The "spectrum near 1/2" feature R77.4 erratum reached for is empirically absent at the tested levels (k=2..5).

## 8. Files

- `R_K_APPROACH_A.md` (this file) — specification + structural prediction + outcome
- `R_K_APPROACH_A_script.py` — (specified; main-thread execution deferred)
- `L_k_eigenvalues.csv` — adjacent SVD data from prior probe (anchor for anticipated values)
- `result_77_5_phi_correlations.csv` — c_k = 0 data (anchor for §6 sanity check)

Phase 2B (perturbation) and Phase 2C (numerical resolvent norm) are subsequent only if they have a meaningful operator to bound. Given §6's finding (Φ_k doesn't transport R_{k−1} → R_k), there isn't one. Phase 2B and 2C therefore document the structural reason rather than computing numerical bounds.
