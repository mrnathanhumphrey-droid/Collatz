# Result 19: Path B continuation — matrix Wiener-Hopf framework completed; per-j W_j requires lattice landing

**Date:** 2026-05-02. Sequel to Results 15, 16, 17.

This document delivers the *completed* matrix Wiener-Hopf framework on the residue-class chain at k=6 (and k=10), including the boundary-residue fix that closes Result 17's 8% drift error. The matrix WH machinery yields a clean **universal** prediction for W_j across all j, equal to the descending ladder mean. The empirical per-j variation (W_2 = +7.16, W_4 = -4.76, W_5 = +4.59 — including the sign flip between j=4 and j=5) is **not** reproduced by the residue chain. The framework's gap is identified precisely: it is integer-lattice landing structure that lives below the residue-chain resolution.

Code: `path_b_matrix_wh_v2.py` (drift exact at 1e-13, k=6 and k=10 sweeps, asymptotic-regime W_j extraction).

---

## 1. What's fixed from Result 17 (partial)

The Q matrix construction at k=6 had a subtle bookkeeping error in v1: for the boundary residue r=21, only the j=0 stratum's log-step was retained per (r,r') pair, giving an 8% drift discrepancy (E_π[X] = -0.266 vs target -0.288 nats).

**Fix:** track the per-residue *step distribution* P(v|r) directly, then compute E[X|r] = log(3) - E[v|r]·log(2). For boundary r=21:

```
P(v=6+j | r=21) = 2^(-j-1)  for j = 0, 1, 2, ...
E[v | r=21] = 6 + Σ j·2^(-j-1) = 6 + 1 = 7
E[X | r=21] = log(3) - 7·log(2) = -3.7534 nats
```

Stationary distribution remains uniform π_r = 1/32. Drift now matches log(3/4) to **machine precision (1e-13)**:

```
E_π[X] = -0.28768207  (target: log(3) - 2 log(2) = -0.28768207)
```

This closes Result 17's drift error. The Q construction is now exact.

## 2. The Markov-additive characteristic exponent factorizes cleanly

For the residue chain on odd residues mod 2^k, the Markov-additive characteristic function is:

> Φ(θ) = D(θ) · Q

where D(θ) is *diagonal* with per-state step CFs:

- For non-boundary r: D(θ)[r,r] = 3^(iθ) · 2^(-iθ·v(r)) [deterministic v]
- For boundary r=21 at k=6: D(θ)[21,21] = 3^(iθ) · 2^(1 − 6iθ) / (2^(1+iθ) − 1) [closed form]

The boundary-residue CF closed form follows from the geometric sum:

```
Σ_{j=0}^∞ 2^(-j-1) · 2^(-iθ(6+j)) = 2^(-6iθ-1) · 1/(1 - 2^(-1-iθ)) = 2^(1-6iθ)/(2^(1+iθ) - 1)
```

This is a **rational function in 2^(iθ)**, the same algebraic class as the iid Path C symbol (Result 15). The per-residue CFs are all rational in 2^(iθ) modulo the universal 3^(iθ) factor. The Markov-additive matrix WH inherits this algebraic structure.

## 3. Spectral structure of Q is stable across k

Eigenvalues at k=6 (32-state):
```
λ_0 = +1.000           (stationary, Perron-Frobenius)
λ_1..31: |λ| ≈ 2.6×10⁻⁴ (rapid mixing)
```

At k=10 (512-state):
```
λ_0 = +1.000
non-stationary eigenvalues all |λ| < 10⁻³
```

The chain mixes essentially in **one step**. Practical consequence: the residue at first descent below any level -L (for L ≫ 0) is approximately uniform over the 32 (or 512) odd residues, regardless of starting state. Conditional residue probabilities don't carry strong j-dependent structure.

## 4. Where the m_j absorbing residues live

| k | Modulus | m_2=5 | m_3=21 | m_4=85 | m_5=341 | m_7=5461 | m_8=21845 |
|---|---|---|---|---|---|---|---|
| 6 | 64 | 5 | 21 (B) | 21 (B) | 21 (B) | 21 (B) | 21 (B) |
| 10 | 1024 | 5 | 21 | 85 | 341 (B) | 341 (B) | 341 (B) |
| 12 | 4096 | 5 | 21 | 85 | 341 | 1365 (B) | 1365 (B) |

(B) = boundary residue at that k (where v_k(3r+1) ≥ k).

**Pattern:** m_j separates from m_{j-1} when 2^k > m_j, i.e., k > 2j - log_2(3) ≈ 2j - 1.585. To distinguish j ∈ {2, 3, 4, 5} requires k ≥ 9, so we use k=10. The boundary residue at each k is *exactly* the smallest m_j with 2^k = 4^j (k=6 → r=21=m_3; k=10 → r=341=m_5; k=12 → r=1365=m_8).

Remarkably: **the boundary residue at each k is itself an attractor m_j**. The j-class structure is encoded *into* the residue chain's boundary structure.

## 5. Numerical W_j extraction via Markov-additive simulation

We run 500K-orbit Markov-additive sims at k=6, walking from S=0 until first descent below -L = -(log(2^36) - log(m_j)), conditioning on end residue at the first crossing. This is the asymptotic regime (L ≈ 19-23 nats ≫ typical step size 0.3 nats).

Convention: W_j (Syracuse-step units) = E[overshoot magnitude] / log(4/3).

| j | m_j | L (nats) | E[overshoot] (nats) | W_j_pred (Syracuse steps) | W_j_emp |
|---|---|---|---|---|---|
| 2 | 5 | 23.34 | 0.989 ± 0.008 | **3.44 ± 0.03** | +7.156 |
| 4 | 85 | 20.51 | 0.981 ± 0.008 | **3.41 ± 0.03** | −4.755 |
| 5 | 341 | 19.12 | 0.991 ± 0.008 | **3.44 ± 0.03** | +4.590 |

**The MAP simulation gives a UNIVERSAL value W_j_pred ≈ 3.44 across all j**, independent of:
- Which target lattice point (m_j)
- Which absorbing residue (r_j)
- Whether the residue is boundary (j=4, 5 at k=6) or non-boundary (j=2 at k=6, all distinct at k=10)

Identification of the universal value: 3.44 ≈ E[L⁻] / log(4/3) = 1.0046 / 0.288, the **strict descending ladder mean** of the Markov-additive walk, in Syracuse-step units. This matches the iid value (Path A baseline) since the residue chain mixes in essentially one step — the Markov correction to E[L⁻] is sub-percent.

## 6. The empirical per-j variation is structurally outside the residue chain

Empirical (50M-orbit, N=2^36):
```
ΔW_j = W_j_emp - 3.44   (deviation from MAP universal prediction)
ΔW_2 = +3.72   (j=2)
ΔW_4 = -8.20   (j=4, sign flipped)
ΔW_5 = +1.15   (j=5)
```

The deviations:
- Are not random (50M-orbit precision pins them at ±0.01)
- Vary by 12 step units across j (more than any natural Markov correction at the residue level)
- Include a **sign flip** between adjacent j (W_4 < 0 < W_5)

**Why the residue chain misses it.** Empirical W_j is conditioned on the orbit landing **exactly on the integer m_j** (i.e., m_τ = m_j). In the continuous MAP framework, this is a measure-zero event. The MAP can condition on (residue at first descent below -L, S at first descent below -L) but cannot condition on the exact integer landing m_t = m_j without coupling to the integer lattice.

**Algebraic identification of the missing structure.** When an orbit at residue r=21 (k=6) has current m at residue 21, the next step's v stratum is v_2(1+3h) where h = (m-21)/64. For absorption at attractor m_j with j ≥ 3:
- m = m_j requires h = (m_j - 21)/64 = (4^(j-3) - 1)/3 = m_(j-3)
- Then v_2(1 + 3·m_(j-3)) = v_2(4^(j-3)) = 2(j-3)
- So absorption at m_j corresponds to the v=6+2(j-3) sub-stratum within the boundary residue's distribution

Each j ≥ 3 absorbs at a different sub-stratum of the boundary residue — but in the continuous MAP, all sub-strata get summed into the boundary-residue marginal, losing the per-j distinction.

**Concrete prediction for closure path.** A higher-resolution analysis tracking (residue mod 2^k, sub-stratum j_inner = v_2(1+3h)) jointly would distinguish per-j W_j. This is a 32×N_strata state space (effectively a deeper Markov-additive process with finer alphabet). Each sub-stratum corresponds to a distinct integer-lattice landing pattern.

## 7. Honest scope statement

**What this delivers:**
- Q construction at k=6 corrected to machine-precision drift (Result 17 fix landed)
- Markov-additive characteristic function Φ(θ) = D(θ)·Q with closed-form per-residue CFs
- Boundary-residue CF in closed form (rational in 2^(iθ))
- Spectral verification (rapid mixing at all k tested)
- Numerical W_j extraction via MAP simulation in asymptotic regime
- Universal value W_j_pred ≈ 3.44 = E[L⁻]/log(4/3) across all j
- Identification of the precise gap to per-j empirical W_j

**What it does NOT deliver:**
- Per-j W_j matching empirical to ±0.05 (the brief's success target). Best gap is +1.15 step units for j=5; W_4 has SIGN FLIP at -8.2 step units gap. The framework's prediction is universal across j; per-j variation is below the residue chain's resolution.
- Closed-form symbolic W_j expressions
- Coupling to integer-lattice landing (the missing piece)

**Verdict per brief's decision criteria:**
- "W_j matches empirical to ±0.05": **NO**
- "W_j matches to ±0.1 but not ±0.05": **NO**
- "W_j off by >0.5": **YES** (gap of 1-8 step units across j)

Per the brief: "If off by >0.5: framework as formulated misses structure. Document the gap." This is that documentation.

## 8. Why this is still substantive — what the framework closes

The matrix WH framework on the residue chain **does** close several things even though it doesn't close per-j W_j:

1. **Confirms the iid Wald-Lorden baseline.** The MAP sim gives E[L⁻] = 1.005 nats, matching the iid Path A baseline (Result 15) to sub-percent. The residue-chain Markov correction to ladder mean is ≤ 0.5%, consistent with rapid mixing.

2. **Identifies the boundary-residue / attractor coincidence as algebraic.** The boundary residue at each k is *exactly* an m_j attractor; this is a structural finding (Result 17) confirmed and given an algebraic reason here (Section 6).

3. **Closes the drift identity.** E_π[X] = -log(4/3) is now provable via the matrix construction, not just observable empirically.

4. **Gives the asymptotic-regime baseline for the bridge equation.** The W_j universal value 3.44 step units is the leading-order term in:

> ε_S = Σ_j P(j) · [W_j − log(m_j)/log(4/3) + 1]

If we replace empirical W_j with the universal MAP prediction 3.44, ε_S_pred = 3.44 - Σ_j P(j) · log(m_j)/log(4/3) + 1. With P(j=2) = 0.938 dominating: ε_S_pred ≈ 3.44 - 0.938·5.59 - 0.024·15.43 - 0.038·20.27 + 1 ≈ 3.44 - 5.245 - 0.370 - 0.770 + 1 = -1.95 step units. Empirical ε_S ≈ 1.375. So the universal W_j gives WRONG ε_S; per-j W_j is essential.

This is consistent: the per-j variation in W_j is what builds ε_S. The framework gives the baseline; the lattice-landing residual is the per-j correction that nets out to ε_S.

## 9. Where the program actually closes — and what's left

**Closed (now structurally confirmed):**
- ⟨α_det⟩ = log(6)/log(4/3) (Result 1)
- P(j) entry distribution via absorbing-Markov chain on integers (compute_threads_findings.md addendum)
- Drift identity E_π[X] = -log(4/3) on residue chain (Result 19)
- Algebraic boundary-residue / m_j coincidence (Results 17, 19)

**Open (named gap from prior session, NOT closed by Path B alone):**
- Per-j W_j: requires either:
  (a) Discrete renewal theory on integer chain — fails by truncation as shown in compute_threads addendum
  (b) Higher-order Markov-additive analysis tracking sub-strata of boundary residue
  (c) Edgeworth-style lattice expansion of the integer-landing distribution

The named open theoretical problem (asymptotic W_j requires absorbing-Markov first-passage at small-target boundary) remains. The matrix WH framework on residues is **necessary but not sufficient**: it gives the universal Lorden-style baseline, not per-j corrections.

## 10. Files

- `path_b_matrix_wh_v2.py` — corrected Q + matrix Φ(θ) construction + asymptotic W_j extraction (k=6 and k=10)
- `path_b_matrix_wh.py` — original (Result 17 partial); kept for reference
- `path_b_derivation.md` — Result 17 partial (preserves boundary-residue / attractor structural finding)
- `path_b_continuation.md` — this document (Result 19)
- `closed_form_findings.md` — entry below

## 11. Citations

- Asmussen 1989 / "Applied Probability and Queues" Ch IV — matrix Wiener-Hopf for Markov-additive
- Alsmeyer-Buckmann 2018 — Markov random walks, fluctuation theory
- Lorden 1970 — asymptotic stationary excess for renewal overshoot
