# CANDIDATE_A_PHI_CONSTRUCTION — explicit construction of φ_n on V_n

**Date:** 2026-05-12. Wilson (analyst) reporting to Nathan. Phase 1 of the Reading A scoping probe.

---

## Goal

Construct φ_n explicitly as a vector in V_n (function space on coprime residues mod 3^n) such that

  ε_n = S_n − 7/15 = ⟨φ_n, π_n⟩ − 7/15

so that ε_n can be expressed via the bilinear pair-form structure of R76 and projected through the W_k filtration.

---

## R76 articulation, gap, and resolution

R76 defines

  S_n = M_n(1) = Σ_{ξ ∈ Z/3^n, 3∤ξ} |μ̂_n(ξ)|²

where μ̂_n(ξ) := Σ_r π_n(r) e^{−2πi r ξ / 3^n} is the characteristic function of π_n. R76 does **not** explicitly write φ_n as a vector in V_n; the project's existing R76 code (`bilinear_pair_operator.py`) computes M_n via complex-valued FFT-style sums. To construct φ_n as a real-valued (in fact, integer-valued up to π_n weighting) functional on V_n, expand S_n directly:

  S_n = Σ_{ξ: 3∤ξ} |μ̂_n(ξ)|²
      = Σ_{ξ: 3∤ξ} Σ_r Σ_s π_n(r) π_n(s) e^{−2πi (r − s) ξ / 3^n}
      = Σ_r Σ_s π_n(r) π_n(s) · K_n(r − s)

where

  **K_n(d) := Σ_{ξ ∈ Z/3^n, 3∤ξ} e^{−2πi d ξ / 3^n}.**

This is a real, integer-valued kernel on Z/3^n. Computing it:

  K_n(d) = (Σ_{ξ ∈ Z/3^n} e^{−2πi d ξ / 3^n}) − (Σ_{ξ ∈ Z/3^n, 3|ξ} e^{−2πi d ξ / 3^n})

Both sums are characters; the first equals 3^n if d ≡ 0 mod 3^n else 0; the second (substitute ξ = 3η, η ∈ Z/3^{n−1}) equals 3^{n−1} if d ≡ 0 mod 3^{n−1} else 0. Hence:

  **K_n(d) = 3^n · 1[d ≡ 0 mod 3^n] − 3^{n−1} · 1[d ≡ 0 mod 3^{n−1}].**

Explicitly on d ∈ Z/3^n:

  K_n(0)              = 3^n − 3^{n−1} = 2 · 3^{n−1}
  K_n(3^{n−1})        = − 3^{n−1}
  K_n(2 · 3^{n−1})    = − 3^{n−1}
  K_n(d)              = 0   for all other d ∈ Z/3^n.

K_n is a 3-valued integer function on Z/3^n, supported on the **5 cosets of 3^{n−1} Z/3^n** (i.e., d ≡ 0 mod 3^{n−1}).

### Definition of φ_n

Define φ_n ∈ V_n by

  **φ_n(r) := Σ_{s coprime in Z/3^n} π_n(s) · K_n(r − s).**

This makes φ_n a function on coprime states in Z/3^n, with values in Q (rational, since K_n is integer-valued and π_n ∈ Q^{N_n}). Then by construction:

  ⟨φ_n, π_n⟩ = Σ_r Σ_s π_n(r) π_n(s) K_n(r − s) = S_n.

### Linearization caveat (load-bearing)

φ_n is *not* a fixed functional independent of π_n — it depends on π_n itself (K_n * π_n is a convolution against the **same** π_n we're contracting against). This is the standard linearization-of-quadratic-form trick: at each n, define φ_n by the level-n π_n, then the question becomes how φ_n projects onto the W_k filtration. The pre-registration treats this construction as the natural lift of R76's bilinear pair-form moment to the W_k decomposition.

---

## Verification at n = 1, 2, 3 (and extended through n = 6)

`candidate_a_compute.py` computes φ_n exactly over Q via `fractions.Fraction` and verifies S_n = ⟨φ_n, π_n⟩ against R76's known values:

| n | S_n = ⟨φ_n, π_n⟩ (exact) | R76 known | match |
|---|--------------------------|-----------|-------|
| 1 | 2/3                      | 2/3       | ✓     |
| 2 | 10/21                    | 10/21     | ✓     |
| 3 | 31370/67963              | 31370/67963 | ✓   |
| 4 | 143195649659456490 / 308468774477179141 | (matches R77.6 ε_4 = S_4 − 7/15 = −11346676448406637/4627031617157687115) | ✓ |
| 5 | (large fraction, matches R77.6) | ✓ | ✓ |
| 6 | (very large fraction, matches R77.6) | ✓ | ✓ |

Floats: S_1 = 0.6667, S_2 = 0.4762, S_3 = 0.4616, S_4 = 0.4642, S_5 = 0.4655, S_6 = 0.4662 — monotonically approaching 7/15 ≈ 0.4667 from below for n ≥ 3, matching R76 / R77.x.

ε_n = S_n − 7/15 reproduced exactly:
- ε_1 = +1/5 ≈ +0.2000 (above limit; n=1 transient)
- ε_2 = +1/105 ≈ +9.524e-3 (above limit)
- ε_3 = −5191/1019445 ≈ −5.092e-3 (crossed below)
- ε_4 ≈ −2.452e-3
- ε_5 ≈ −1.152e-3
- ε_6 ≈ −4.979e-4

**Verdict:** φ_n construction is exact, reproduces S_n, reproduces ε_n. Phase 1 unblocked.

---

## Structural property of φ_n (load-bearing for Phase 4)

K_n(d) is supported on d ≡ 0 mod 3^{n−1}. So φ_n(r) depends only on the **3-fiber-anti-mean** behavior of π_n at scale 3^{n−1}:

  φ_n(r) = 2 · 3^{n−1} · π_n(r) − 3^{n−1} · (π_n(r + 3^{n−1}) + π_n(r + 2 · 3^{n−1}))
         = 3^{n−1} · (3 · π_n(r) − (π_n(r) + π_n(r + 3^{n−1}) + π_n(r + 2 · 3^{n−1})))

(where indices are mod 3^n on coprime states). The second term in the last line is 3 × (3-fiber average of π_n at r). So

  **φ_n(r) = 3^n · (π_n(r) − π̄_n(r))**, where π̄_n(r) is the 3-fiber-average at the level-(n−1) fiber containing r.

This means **φ_n has zero 3-fiber-mean at level 3^{n−1}** — i.e., φ_n ∈ W_{n−1} ⊂ V_n by definition of W_{n−1} (R77.5 §3).

This is **the critical structural fact that determines Phase 4's outcome**: φ_n lives entirely in the finest-scale W_{n−1} subspace, not spread across W_0 ⊕ ... ⊕ W_{n−1}.

---

## Files

- `candidate_a_compute.py` — full computation
- `candidate_a_phi_n_verify.csv` — exact-rational S_n / ε_n / sum-c table per level
