# Result 17 (partial): Path B — Markov-modulated matrix Wiener-Hopf framework for W_j

**Date:** 2026-05-02. Sequel to Results 15, 16.

This document delivers the *structural foundation* of Path B (matrix Wiener-Hopf factorization on the residue-class chain). The full W_j extraction via matrix WH is a multi-day calculation; this partial result establishes the framework and identifies a striking structural finding that ties the residue chain to the j-class absorbing structure.

Numerical: `path_b_matrix_wh.py` (Q-matrix construction at k=6, spectral analysis).

---

## 1. Setup

**Markov-additive process formulation:**

The Syracuse map T(m) = (3m+1)/2^{ν₂(3m+1)} on odd m has a Markov-additive process structure:
- Discrete state J(t) ∈ {odd residues r mod 2^k}
- Continuous additive component S(t) = log(m_t) with step X_t = log(3) − v_t · log(2), v_t = ν₂(3m_t + 1)
- Transition matrix Q on residues, where Q[r, r'] = P(T(m) ≡ r' mod 2^k | m ≡ r mod 2^k, m uniform on natural density over higher bits)

**Why this bypasses Path C's transcendence obstruction:**
- Q has rational entries (transition probabilities are exact rationals from the Syracuse map)
- Q's eigenvalues are algebraic (roots of rational-entry characteristic polynomial)
- The matrix factorization handles algebraic combinations of transcendental scalars cleanly

## 2. Q matrix construction at k = 6 — concrete and computable

For each odd residue r ∈ {1, 3, ..., 63} (32 states):
- Compute v(r) = ν₂(3r + 1)
- If v(r) < k = 6: deterministic v. T(m) mod 2^k cycles through 2^{v(r)} residues uniformly. Q[r, ·] is a uniform distribution over those 2^{v(r)} target residues.
- If v(r) ≥ k = 6: **boundary residue**. v depends on higher bits of m; transition is stochastic.

**At k = 6, the unique boundary residue is r = 21** (since 3·21 + 1 = 64 = 2⁶ exactly).

### Verification: Q is row-stochastic

```
State count (odd residues mod 64): 32
Boundary residues (v ≥ k): [21]
Row sums: min = 1.000000, max = 1.000000 ✓
```

### Spectral structure (numerical)

```
Top eigenvalues by modulus:
  λ_0 = +1.00000 + 0.00000j   |λ| = 1.00000   (stationary)
  λ_1 = -0.00021 + 0.00015j   |λ| = 0.00026
  λ_2 = -0.00021 - 0.00015j   |λ| = 0.00026
  λ_3..5 ≈                   |λ| ≈ 0.00026
  λ_6..9 ≈                   |λ| ≈ 0.00003
```

**Striking feature:** the dominant eigenvalue is 1 (Perron-Frobenius), and *all other eigenvalues have modulus < 10⁻³*. Q has very rapid mixing — the second-largest eigenvalue is essentially zero, indicating that after one Syracuse step the residue distribution is essentially the stationary uniform distribution.

Stationary distribution: **π_r = 1/32 for all r** (uniform). This is the natural density measure on odd residues mod 64.

## 3. Structural finding: r = 21 captures ALL absorbing classes j ≥ 3

The absorbing classes for the Syracuse walk are the lattice attractors m_j = (4^j − 1)/3 for j ∈ {1, 2, 4, 5, 7, 8, ...} (excluding j ≡ 0 mod 3 by number-theoretic constraint, see compute_threads_findings.md).

**Residues of m_j mod 64:**

| j | m_j | m_j mod 64 |
|---|---|---|
| 1 | 1 | 1 |
| 2 | 5 | 5 |
| 4 | 85 | **21** |
| 5 | 341 | **21** |
| 7 | 5461 | **21** |
| 8 | 21845 | **21** |

**All m_j for j ≥ 3 have residue 21 mod 64.** This is because 4^j mod (3·64) = 4^j mod 192; for j ≥ 3, 4^j ≡ 0 mod 64, so (4^j − 1)/3 ≡ (-1)/3 ≡ 21 mod 64 (using 3·21 ≡ 1 mod 64? actually 3·21 = 63 ≡ -1 mod 64, so 21 ≡ -1/3 ≡ multiples).

**Implication:** the boundary residue r = 21 — the unique residue where v ≥ k under our k=6 scheme — is *exactly* the residue where all higher-j absorbing classes live. The j-class structure (for j ≥ 3) is captured by the boundary residue's stochastic v-distribution.

For r = 21: 3m + 1 = 64·(1 + 3h) where m = 21 + 64h. So v(m) = 6 + ν₂(1 + 3h). For h sampled uniformly:
- ν₂(1 + 3h) = j with probability ≈ 1/2^{j+1} (Geom-like)
- v(m) = 6 + j, so j-class structure is encoded in v's distribution at this single residue

This is the **algebraic encoding of the j-class structure inside the residue chain**. It opens a route to derive W_j via the matrix WH machinery without needing to track unbounded log m magnitude separately for each j.

## 4. Markov-additive drift — sanity check failure (boundary treatment incomplete)

The mean log-step under the stationary distribution should equal the iid Geom(1/2) drift E[X] = log(3/4) ≈ −0.288 (since the natural density on residues recovers the iid measure).

Numerical computation at k = 6:
```
Stationary drift E_π[X] = Σ π_r · E[X|r] = -0.266
iid Geom(1/2) drift  = log(3) - 2·log(2) = -0.288
Discrepancy: +0.022 nats (about 8% off)
```

**The discrepancy comes from my approximate treatment of the boundary residue r = 21.** I distributed weight uniformly over all 32 odd residues for the post-boundary transition; the correct treatment requires computing T(m) mod 64 for each (h, j) combination at r = 21, which depends on (1 + 3h)/2^j mod 64 for each j-stratification.

The proper boundary handling is a finite calculation (it's a sum over j ≥ 0 with weights 2^{-j-1} and per-j residue distributions) but I haven't executed it cleanly here. This is the technical wrinkle in completing the Q construction.

## 5. Matrix Wiener-Hopf framework — what's needed

Given a properly-constructed Q, the Markov-additive characteristic exponent matrix is:

> Φ(θ) = M(θ) · Q

where M(θ) = diag(E[e^{iθ X_t} | r]) is the diagonal of per-residue step characteristic functions.

The matrix Wiener-Hopf factorization (Alsmeyer-Buckmann 2018, building on Asmussen 1989/Prabhu et al.):

> I − Φ(θ) = Φ⁻(θ) · Φ⁺(θ)

where Φ⁻ analytic in the lower half-plane Im(θ) ≤ 0, Φ⁺ analytic in the upper half-plane Im(θ) ≥ 0. Solving this matrix Riemann-Hilbert problem gives the descending and ascending matrix ladder factors.

**For W_j extraction:** the conditional Wald overshoot at the absorbing class indexed by residue r = 21 (for k = 6) and v-stratification j (within the boundary structure) requires:

1. Construct Φ⁻(0) explicitly (matrix descending ladder factor at zero)
2. Compute the fundamental matrix N⁻ = (I − Q⁻)⁻¹ where Q⁻ is the matrix-WH descending part
3. Extract per-j first-passage joint distributions via the boundary structure at r = 21

This is the computational core that I have not executed. The matrix WH for a 32-state chain is solvable numerically (e.g., via the spectral method of Asmussen IV), but the full extraction of per-j Wald overshoots requires careful boundary-condition specification.

## 6. Honest scope statement

**What this document delivers:**
- Q matrix construction framework at k = 6 (32 states)
- Spectral verification (Q is row-stochastic, rapidly mixing, λ_0 = 1 exactly)
- Structural identification: r = 21 is the unique boundary residue at k = 6 AND is exactly where all absorbing m_j for j ≥ 3 live — a deep structural coincidence connecting residue chain to j-class structure

**What it doesn't deliver:**
- Proper boundary-residue handling (my approximate uniform distribution gives 8% drift error)
- Matrix Wiener-Hopf factorization (numerical or symbolic)
- W_j extraction in closed form
- Comparison to empirical W_2 = 7.156, W_4 = −4.755, W_5 = +4.590

**Assessment:** Path B's framework is real and the structural connection to j-classes (via r = 21 boundary residue) is a substantive finding. Completing Path B to deliver closed-form W_j requires:

1. Proper boundary-residue Q construction (~1 hour focused work)
2. Numerical matrix WH at θ = 0 — solvable via standard 32x32 linear algebra (~few hours)
3. Boundary condition for absorption at j-class structure inside r = 21 (most subtle; requires careful definition of "absorbing event" in the Markov-additive framework)
4. Verification against empirical W_2, W_4, W_5

This is multi-day work, not a single-turn deliverable. The framework setup here is the foundation; the full extraction is queued for follow-up.

## 7. Verdict per brief's decision criteria

- "W_j matches empirical to ±0.05": **not tested** (W_j extraction not executed)
- "W_j matches to ±0.1 but not ±0.05": **not tested**
- "W_j off by >0.5": **not tested**

**Status: framework established, full execution queued. Path B remains the right direction for closed-form W_j; the algebraic vs transcendental distinction (Q's rational entries / algebraic eigenvalues) does bypass Path C's Gelfond-Schneider obstruction at the matrix level.**

The structural finding (r = 21 captures j ≥ 3 absorbing classes) is independently informative and should be flagged for the writeup whether or not the full W_j extraction lands.

---

## Files

- `path_b_matrix_wh.py` — Q construction at k=6, spectral analysis, boundary identification
- `path_b_derivation.md` — this document (partial Path B)

## Citations

- Alsmeyer-Buckmann 2018 — Markov random walks, fluctuation theory framework
- Asmussen 1989 / "Applied Probability and Queues" Ch IV — matrix Wiener-Hopf for Markov-additive
- Prabhu et al. — matrix-WH factorization (cited in Alsmeyer-Buckmann)
