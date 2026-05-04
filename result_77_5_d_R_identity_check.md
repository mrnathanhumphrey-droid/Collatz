# R77.5 follow-up — Identity check: ‖R_k‖² = ‖d_{k+1}‖² over Q

**Date:** 2026-05-04. Tests whether R77.5's lift-residual norm ‖R_k‖² and R74's deviation norm ‖d_{k+1}‖² are the same exact rational, equivalently whether ‖R_k‖² · 3^k = S_{k+1}/3 → 7/45 is a rigorous identity inherited from R74 or an empirical proximity at the 0.001 level.

## Verdict in one line

> **Outcome (IDENTITY).** ‖R_k‖² = ‖d_{k+1}‖² and ‖R_k‖² · 3^k = S_{k+1}/3 hold as **exact rational equalities at every k = 1, 2, 3, 4, 5** (5/5 PASS, 5/5 PASS). The identity is algebraic; the empirical proximity 0.155 ≈ 7/45 is the consequence, not a coincidence.

**Files (this result):**
- `result_77_5_d_R_identity_check.md` (this writeup)
- `result_77_5_d_R_identity_check.py` (verification script)
- `result_77_5_d_R_norms.csv` (exact-rational table with both norms and the test booleans)

---

## 1. The two definitions

**R74 deviation norm** (from `c_seven_forty_fifth_derivation.py:compute_d_squared_exact`):

> ‖d_{k+1}‖²_R74 := Σ_{r' coprime in Z/3^{k+1}} π_{k+1}(r')² − (1/3) · Σ_{r coprime in Z/3^k} π_k(r)²

This is a scalar (sum), not a vector — the "deviation" is implicit.

**R77.5 lift residual norm:**

> R_k(r') := π_{k+1}(r') − T(π_k)(r'),    T(π_k)(r') := π_k(r' mod 3^k) / 3
>
> ‖R_k‖² := Σ_{r' coprime in Z/3^{k+1}} R_k(r')²

This is a vector L² norm.

## 2. Algebraic identity ‖R_k‖² = ‖d_{k+1}‖²

Expand the L² norm:

  ‖R_k‖² = Σ π_{k+1}(r')² − 2·Σ π_{k+1}(r')·T(π_k)(r') + Σ T(π_k)(r')².

For each term:

**Cross term.** Σ_{r'} π_{k+1}(r') · T(π_k)(r') = Σ_{r'} π_{k+1}(r') · π_k(r' mod 3^k) / 3. Group by r = r' mod 3^k:
  = (1/3) · Σ_{r coprime in Z/3^k} π_k(r) · [Σ_{r' lifts of r} π_{k+1}(r')]
  = (1/3) · Σ_r π_k(r) · π_k(r)         (by **marginal consistency** of the projective Markov system)
  = (1/3) · Σ π_k².

**Lift self-norm.** Σ_{r'} T(π_k)(r')² = Σ_{r'} (π_k(r' mod 3^k)/3)² = (1/9) · 3 · Σ_r π_k(r)² = Σ π_k² / 3.

**Combine:**

  ‖R_k‖² = Σ π_{k+1}² − 2·(1/3)·Σ π_k² + (1/3)·Σ π_k² = **Σ π_{k+1}² − (1/3)·Σ π_k² = ‖d_{k+1}‖²_R74.** □

Note: the marginal consistency Σ_{r' lifts of r} π_{k+1}(r') = π_k(r) is the only non-trivial input. It holds because the Syracuse Markov chain is coherent under reduction mod 3^k — verified via standard projective-system arguments and confirmed numerically at every k.

## 3. Empirical verification at k = 1..5 (all over Q via fractions.Fraction)

Tests run by `result_77_5_d_R_identity_check.py`:

**Test (A): ‖R_k‖² = ‖d_{k+1}‖²?**

**Test (B): ‖R_k‖² · 3^k = S_{k+1}/3?** (where S_{k+1} := X_{k+1} − X_k, X_j := 3^j · Σ π_j² — R74's pre-S formulation)

| k | ‖R_k‖² (decimal) | ‖d_{k+1}‖² (decimal) | (A) | ‖R_k‖²·3^k | S_{k+1}/3 | (B) |
|---|---|---|-----|---|---|-----|
| 1 | 5.291×10⁻² | 5.291×10⁻² | **PASS** | 0.1587302 | 0.1587302 | **PASS** |
| 2 | 1.710×10⁻² | 1.710×10⁻² | **PASS** | 0.1538582 | 0.1538582 | **PASS** |
| 3 | 5.731×10⁻³ | 5.731×10⁻³ | **PASS** | 0.1547381 | 0.1547381 | **PASS** |
| 4 | 1.916×10⁻³ | 1.916×10⁻³ | **PASS** | 0.1551716 | 0.1551716 | **PASS** |
| 5 | 6.395×10⁻⁴ | 6.395×10⁻⁴ | **PASS** | 0.1553896 | 0.1553896 | **PASS** |

5/5 PASS for (A), 5/5 PASS for (B). All comparisons over Q via fractions.Fraction; not float comparisons.

Convergence to 7/45 = 0.15555…:

| k | ‖R_k‖² · 3^k | 7/45 − ‖R_k‖²·3^k |
|---|---|---|
| 1 | 0.158730 | −0.003175 |
| 2 | 0.153858 | +0.001697 |
| 3 | 0.154738 | +0.000817 |
| 4 | 0.155172 | +0.000384 |
| 5 | 0.155390 | +0.000166 |

Monotone-from-below for k ≥ 2; the deviation 7/45 − ‖R_k‖²·3^k = (S_∞ − S_{k+1})/3 = −ε_{k+1}/3, matching the empirical |ε_n|·2^n ≈ 0.04 envelope (so deviation ∼ 0.04/3 · 2^{−k}).

## 4. Strategic implications for v3.7.1

**R77.5's contribution is geometric reframing of R74, not new theorem content.**

- The vector R_k and the R74 quantity d_{k+1} are **literally the same vector** in the same basis (both are π_{k+1} − T(π_k)).
- The L² norms agree as exact rationals, not just empirically.
- The "multi-resolution decomposition V_{k+1} = T(V_k) ⊕ W_k with R_k ∈ W_k" framing in R77.5 is a clean geometric language for the same object R74 was already studying.
- The 7/45 limit is **rigorous identity** modulo the open S_∞ = 7/15 conjecture. It is not a coincidence at the 0.001 level.

**What this means for R77.2's displacement claim:**

R77.5's writeup argued that "rate-1/2 lives in the moment functional φ_n's projection onto Σ_k W_k, not in any single operator's spectrum." That conclusion stands — the multi-resolution / wavelet-like geometry is real, structurally tied to marginal consistency, and a legitimate framework for further analysis. But:

- R77.5 does NOT introduce a new mathematical object beyond R74's d_{k+1}.
- It introduces a new **language** (multi-resolution / projective-limit / wavelet) for R74.
- The transfer-operator analysis on Ẑ_3^× that R77.5 recommended is still the right next step — the projective-limit framing IS load-bearing — but the underlying object is unchanged from R74.

**For v3.7.1's Path I:** R77.5 is a "framing chapter," not a "result chapter." The substantive open problem remains rate-1/2 closure, which now has the clean multi-resolution language for stating the problem precisely:

> **Conjecture (rate-1/2 in multi-resolution form).** The bilinear pair-form moment functional φ_n decomposes against the orthogonal-complement filtration W_k as Σ_k ⟨φ_n, lift_n(R_k)⟩, and the leading rate of this sum is exactly 1/2 in n.

This is mathematically equivalent to the original |ε_n|·2^n bounded conjecture, but stated in a basis where the structure is geometric.

## 5. What this does NOT close

- **S_∞ = 7/15 itself remains open.** The 7/45 limit of ‖R_k‖²·3^k follows from S_∞ = 7/15, which is itself a separate conjecture (empirically certified but not proven).
- **Rate-1/2 of ε_n remains open.** The identity |ε_n|·2^n ≈ 0.04 (bounded) is the substantive analytic claim; R77.5's reframing doesn't prove it.
- **Tao Prop 1.17 effective C_A** remains the gate to a Nisoli-style rigorous closure of rate-1/2, if that route is pursued. R77.5's reframing suggests an alternate route via transfer-operator analysis on Ẑ_3^×, but that's a substantial separate construction.

## 6. Anti-pattern audit

- **Did not** claim identity from float agreement. Tested over Q via fractions.Fraction at k=1..5; all 10 tests (5 × 2 identities) PASS as rational equalities.
- **Did derive** the algebraic identity (§2) so the test is checking against a proof, not a hope.
- **Did identify** the structural source of the identity (marginal consistency of the projective Markov system). This is the same fact that caused the c_k = 0 result in R77.5 Stage 2.

## 7. Files

- `result_77_5_d_R_identity_check.py` — verification script (runs cleanly via `python ...`)
- `result_77_5_d_R_norms.csv` — exact-rational table with numerator/denominator pairs and both test booleans
