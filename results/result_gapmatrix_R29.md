# Probe R29 — the gap matrix — **the matrix method FAILS the sanity gate; |λ₂|≈½ survives independently (R29-D)**

**Date:** 2026-07-22  λ=½. Probe `probes/probe_gapmatrix_R29.py` (reuses R28 ν̂/R machinery). Tests whether the
gap-index operator `M_r[m,m']=3·P(2(m'−m))·κ_{2m}(r)` (a fixed lattice `d∈2ℤ`, so buildable unlike R26's growing
spaces) diagonalizes to give `|λ₂|` directly. **It does not: the construction fails its own sanity gate.**

## R29-C — does κ stabilize? **κ_0=1 exact, but κ_2 DRIFTS** (the obstruction, run first)
`κ_{2m}(r) = R_r(2m) / [3·Σ_{m'} P(2(m'−m))·R_{r−1}(2m')]`, D=10:

| r | κ_0 | κ_2 (m=1) | κ_4 (m=2) | κ_6 (m=3) | chan1/S_r | chan2/S_r | chan3/S_r |
|---|---|---|---|---|---|---|---|
| 2 | 1.00000 | 0.48571 | 0.34286 | 0.66667 | 0.700 | 0.175 | 0.109 |
| 5 | 1.00000 | 0.40469 | 0.38655 | 0.70024 | 0.759 | 0.133 | 0.089 |
| 7 | 1.00000 | **0.36048** | 0.36353 | 0.70320 | 0.767 | 0.130 | 0.084 |

`κ_0=1` exactly (Diagonal Flatness). But **`κ_2` drifts down monotonically (0.486→0.360, still declining at r=7)** —
it does *not* stabilize. `κ_4, κ_6` are stabler. The *S_r-normalized channel contributions* stabilize (chan1≈0.77,
chan2≈0.13, chan3≈0.086), but κ — what the matrix is built from — does not. **The matrix has no clean r→∞ limit.**

## R29-A — build & diagonalize: **FAILS the sanity gate** (leading ≠ 1)
| D | r=7 leading | \|λ₂\| |
|---|---|---|
| 4 | 1.0859 | 0.925 |
| 6 | 1.0834 | 0.972 |
| 8 | 1.0832 | 0.980 |
| 10 | 1.0830 | 1.005 |

**The leading eigenvalue is 1.08, not 1** — the pre-registered sanity gate (`leading→1=ρ at criticality`) **fails**,
so by the pre-registration the construction is wrong. And `|λ₂|` does *not* approach ½ — it drifts toward 1. The
structural reason: `3·diag(κ)·K` has the `3K` convolution symbol peaking at 5/3, only damped to ~1.08 by κ — a
**spurious mode**. The physical R-vector grows *linearly* (R(0)=3X_r, eigenvalue 1 with a Jordan block), so the
physical rate is a *subdominant, non-leading* eigenvalue of M, contaminated by the spurious 1.08. Building M from a
single `R_{r−1}→R_r` pair via κ reproduces that one vector but not the operator's spectrum.

## R29-B — truncation convergence: **does NOT converge** ⟹ A void
`|λ₂|` = 0.925 (D=4) → 0.972 → 0.980 → 1.005 (D=10) — **increasing toward 1, D-steps not shrinking** (4.7e−2, 7.8e−3,
2.5e−2). By the pre-registered criterion ("D=8,D=10 disagree ⟹ truncation illegitimate, A void"), the truncation is
illegitimate. The eigenvalues are truncation artifacts, not a converged spectrum.

## R29-D — relative deflation residual: **|λ₂| ≈ ½, confirmed independently of the matrix**
Exact Λ local ratios bracket ½: `Λ₄/Λ₃ = 0.492669`, `Λ₅/Λ₄ = 0.502757`. `μ_r/|Λ_r|` for `μ_r=Λ_{r+1}−z·Λ_r`:

| z | r=3 | r=4 |
|---|---|---|
| 0.49 | +2.7e−3 | +1.3e−2 |
| 0.50 | −7.3e−3 | +2.8e−3 |
| 0.503 | −1.0e−2 | −2.4e−4 |

No single z zeros both residuals — the two local ratios **straddle ½ symmetrically** (0.4927, 0.5028; mean 0.4977),
consistent with **`|λ₂| = ½` exactly** (the deviation is the secondary mode, not a shift of the real mode). This is
the cleanest independent read, agreeing with R18-A (0.493/0.503) and R27 ({3,4}=0.503).

## Status
**R29: the gap-matrix eigenproblem FAILS; the answer survives.** **C** — `κ_2` drifts (does not stabilize), so the
matrix has no clean limit. **A** — diagonalizing `3·diag(κ)·K` gives **leading eigenvalue 1.08 ≠ 1**, failing the
pre-registered sanity gate: the `diag(κ)K` ansatz carries a spurious mode (the `3K` symbol peaks at 5/3) and does
**not** reproduce the physical spectrum (R grows *linearly*, rate 1, a defective eigenvalue the matrix misses).
**B** — `|λ₂|` does not converge in D (drifts to 1), so the truncation is illegitimate and A is void. **D** — `|λ₂|≈½`
is confirmed *independently* (exact Λ ratios straddle ½), consistent with 1/2 exactly.

**Consequence for the crux (owed to the pen):** the "unblocked eigenproblem" is **not** unblocked by this
construction — the growing-spaces obstruction survives in disguise as spurious modes of the single-pair `diag(κ)K`
matrix (leading 1.08, non-convergent in D, κ drifting). So R29 does **not** deliver `|λ₂|` or the period from a
finite eigenproblem. **What does survive:** `|λ₂| ≈ ½` (R29-D, exact Λ ratios; + R18-A + R27), so **gap-survival is
unchanged** — it continues to rest on the converging numerical lines, not a rigorous finite matrix. The correct
finite operator (if one exists) is *not* `3·diag(κ)·K` — the stable object is the *S_r-normalized channel
structure* (chan_m/S_r, which does stabilize), not κ; reformulating the operator on that basis is the pen's item.
The period remains unresolved. No fitting; exact κ_0=1 and Λ-ratio gates, labeled spectra with the sanity-gate
failure and D-non-convergence reported plainly as a failure, not smoothed.
