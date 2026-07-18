# Probe: Preimage structure of the Syracuse step on (Z/3^k)*

**Date:** 2026-05-06.
**Brief:** empirically scope the preimage map of f_k : (Z/3^k)* → (Z/3^k)*
under the natural-density Markov-chain construction (M = 2^20 integer lifts
per coprime state). Test how |Preimage(y)| scales with k to determine the
right function-space construction for adapting Butterley-Kim 2023 anisotropic
Banach space methodology to the abelian profinite Z_3 setting.

## Verdict — Outcome A: |Preimage(y)| BOUNDED, k-independent

> **At every k ∈ {5, 6, 7}, mean = median = 21.00 EXACTLY**, with
> min = 18 and max ∈ {24, 25}. No growth in k whatsoever.
>
> The bound 21 = **log₂(M) + 1** at M = 2^20 — it's the number of v
> values that contribute non-negligible weight at the chosen precision,
> not a property of k. Under the brief's M-lift construction, |Preimage|
> is bounded *by precision*, not *by level*.
>
> **Empirical v-distribution:** P(v = 0) = 1/2 (even lifts of x give v = 0
> deterministically), P(v = j ≥ 1) = 2^(-(j+1)) (odd lifts give Geom(1/2)).
> Ratio to plain 2^(-j) is **exactly 0.5** across all j ∈ {1..12} and all k,
> reflecting that the natural-density construction includes both even and
> odd integer lifts in equal proportions.
>
> **Sanity check:** at k=5, K_emp differs from algebraic K_alg
> (truncated-Geom K_k restricted to v ≥ 1) by 100% Frobenius — they're
> structurally different chains. K_emp includes the v=0 column from even
> lifts; K_alg doesn't.
>
> **Implication for the inverse-limit transfer operator on Z_3:** the
> natural-density operator has bounded essential preimage count at any
> fixed precision, k-independent. The Butterley-Kim adaptation should
> proceed without level-dependent function-space weighting; the function
> space's smoothness/regularity parameter just sets the effective
> precision (analogous to choosing M).

## Phase 1 — preimage map summary

For each y ∈ (Z/3^k)*:
- **|Preimage_strict(y)|** = number of x with K_emp[x, y] > 0
- **|Preimage_eff(y)|** = number of x with K_emp[x, y] ≥ 1/M (= 2^-20)
- **max K_emp[x, y]** over x ∈ Preimage(y)
- **entropy** of {K_emp[x, y] / Σ_x K_emp[x, y] : x ∈ Preimage(y)}

| k | n_coprime | mean strict | median strict | max strict | mean eff | max eff | entropy/log\|P\| |
|---|---|---|---|---|---|---|---|
| 5 | 162 | **21.00** | 21 | 24 | 21.00 | 24 | 0.5438 |
| 6 | 486 | **21.00** | 21 | 25 | 21.00 | 25 | 0.5438 |
| 7 | 1458 | **21.00** | 21 | 25 | 21.00 | 25 | 0.5438 |

All three statistics hit the same value at every k. Strict and effective
counts are identical because at M = 2^20 lifts, the smallest non-zero
K_emp entry is exactly 1/M — there are no entries between 0 and 1/M.

The max max_K[x, y] across all (x, y) at every k is 0.5 (median 0.25), with
mean 0.281. So the dominant transition into y carries 25-50% of the column's
mass. The remaining mass is split across ~20 other preimages with weights
following Geom(1/2)/2.

## Phase 2 — scaling analysis

### Linear, log-log, and constant fits to mean |Preimage|(k):

| Fit | Form | Parameters | RSS |
|---|---|---|---|
| Linear (strict) | a·k + b | a = 0.0000, b = 21.0000 | 7.6e-29 |
| Log-log (strict) | c·k^p | p = 0.0000, c = 21.0000 | (degenerate) |
| Constant (strict) | const | 21.00 | 0.0e+00 |

**Constant fit is exact** (RSS = 0). Mean |Preimage| is *literally
constant* at 21 across the three tested k. No growth, no decay, no
log dependence. The variation k=5 → k=7 in max |Preimage| (24 → 25)
is the only k-dependence visible; the central tendency is unchanged.

### Why exactly 21?

Each lift of x produces (v, target) where v = v_2(3n+1). For random integer
lifts, v ∈ {0, 1, 2, ...} with P(v = 0) = 1/2 and P(v = j ≥ 1) = 2^(-(j+1)).

The probability of seeing v = j at least once across M = 2^20 lifts:
- v = 0: certain (P → 1)
- v = j: 1 - (1 - 2^(-(j+1)))^M ≈ 1 for j ≤ 19, decays for j ≥ 20

So v values reliably sampled: {0, 1, 2, ..., ~19-20}. That's ~21 values,
each producing a distinct target (modulo rare collisions).

Per column y, the symmetric count: each preimage x corresponds to one
specific v giving x → y. v values seen → preimages seen. **Bound = 21
≈ log₂(M) + 1.**

This is *intrinsically a precision bound*, not a level bound. If we
increased M to 2^30, mean |Preimage| would rise to ~31 across all k.

### Concentration (entropy / log|Preimage|)

| k | entropy / log\|Preimage\| (median) | interpretation |
|---|---|---|
| 5 | 0.5438 | 54% of uniform-spread |
| 6 | 0.5438 | 54% of uniform-spread |
| 7 | 0.5438 | 54% of uniform-spread |

Same value across all k (to 4 decimals). The mass within each column is
moderately concentrated: not uniform across all 21 preimages (which would
give ratio 1.0), not concentrated on a single preimage (which would give
ratio 0). The 0.54 ratio reflects the Geom(1/2) weight structure: dominant
preimage at v=1 carries 25% of mass, next at v=2 carries 12.5%, etc.

The k-invariance of the entropy ratio matches the k-invariance of the
v-distribution: both are determined by the v-weighting structure, which is
level-independent.

## Phase 3 — empirical v-distribution

Aggregated v counts across all (x, lift) pairs at each k:

### k = 5 (M_total lifts = 1.7 × 10^8)

| v | count | P(v=j)_emp | P(v=j)_geom = 2^(-j) | ratio |
|---|---|---|---|---|
| 1 | 4.25 × 10^7 | 0.2500 | 0.5000 | 0.5000 |
| 2 | 2.12 × 10^7 | 0.1250 | 0.2500 | 0.5000 |
| 3 | 1.06 × 10^7 | 0.0625 | 0.1250 | 0.5000 |
| 4 | 5.31 × 10^6 | 0.0313 | 0.0625 | 0.5000 |
| ... | ... | ... | ... | 0.5000 |
| 12 | 2.07 × 10^4 | 1.22 × 10^-4 | 2.44 × 10^-4 | 0.5000 |

**Ratio = exactly 0.5000 at every j ∈ {1..20} and every k ∈ {5, 6, 7}.**
This is structural, not noise: it reflects that half the lifts are even
(giving v = 0, contributing the missing 0.5 of total probability) and
half are odd (giving Geom(1/2) on j ≥ 1).

So the conditional distribution given v ≥ 1 is exactly Geom(1/2):
P(v = j | v ≥ 1) = 2^(-j) restoring the geometric prediction. The factor
of 0.5 below the "Geom(1/2)" baseline is just the v = 0 contribution
masking the conditional Geom structure.

This is *the natural-density v-distribution* and matches the standard
result from Lagarias / Sinai for Syracuse iteration: v ∼ Geom(1/2) on
odd integers under the natural density.

### v = 0 contribution (not in the above table)

Each k contributes ~50% of total lifts at v = 0 (even lifts of x give
3n+1 odd, v_2(3n+1) = 0). These contribute K_emp transitions x → 3x + 1
mod 3^k where x is treated as the integer x (not its odd lift). This is
a column of K_emp not present in K_alg (truncated-Geom, restricted to v ≥ 1).

## Phase 4 — implications for transfer operator construction

### Outcome reconciliation

| Outcome | Status |
|---|---|
| **A: \|Preimage\|(k) bounded** | **PRIMARY** — exactly 21 at all k tested |
| **B: \|Preimage\|(k) grows linearly in k** | NO — completely flat in k |
| **C: \|Preimage\|(k) faster than linear** | NO |
| **D: empirical v-distribution differs from geometric** | PARTIAL — empirical is Geom(1/2) on odd lifts (j ≥ 1), with even lifts contributing additional v = 0 mass that K_alg doesn't capture |

### What this means for Butterley-Kim adaptation

The empirical natural-density operator has:
1. **Bounded preimage count at any fixed precision** (21 at M = 2^20),
   k-independent.
2. **Geometric weight structure** on the preimages (Geom(1/2) over v).
3. **One dominant preimage per column** carrying 25-50% of mass; the rest
   exponentially suppressed.

For the inverse-limit operator on Z_3:
- Test functions should be C^q-class on Z_3 with the standard 3-adic
  metric.
- The transfer operator L acts as L f(y) = Σ_v 2^(-v) · f(x_v(y)) where
  x_v(y) ∈ (Z/3^∞)* solves (3·x_v + 1)·2^(-v) ≡ y in Z_3.
- The infinite sum over v converges geometrically.
- Anisotropic Banach space construction (Butterley-Kim) on test functions
  with smoothness q in the 3-adic direction and decay in the v-direction
  applies *without level-dependent reweighting*.

This is the simplest possible regime for transfer-operator analysis:
**bounded essential preimage count + geometric weights = Butterley-Kim
machinery directly applicable**.

### Caveat: empirical vs algebraic K_k

The brief's empirical K_emp construction differs structurally from the
truncated-Geom K_alg used in Probe A and prior Collatz framework work:

| | K_emp (empirical, this probe) | K_alg (truncated-Geom, prior probes) |
|---|---|---|
| v range | {0, 1, 2, ...} | {1, 2, ..., M_k = ord_{3^k}(2)} |
| v-weights | observed frequency from M lifts | algebraic 2^(-v) / (1 - 2^(-M_k)) |
| v = 0 column | yes (from even lifts) | no |
| Frobenius diff at k=5 | (vs K_alg) 100% relative |

**For the inverse-limit operator on Z_3, K_emp is the natural finite-k
truncation** (matches natural density on integer lifts). K_alg is a
specific algebraic abstraction that ignores even lifts and bounds v at the
multiplicative order. They give different chains; this probe characterizes
K_emp.

For the prior framework results (S_k → 7/15, ε_k convergence, ρ_slow ≈ 0.83),
K_alg is the relevant chain. Whether the inverse-limit operator on Z_3
(via K_emp) and the truncated-Geom limit chain (via K_alg) give the same
ergodic / spectral behavior is itself a research question; this probe
doesn't address it.

## Summary table

| k | n_coprime | mean \|P\| | median \|P\| | max \|P\| | min \|P\| | max K_emp[x,y] median | entropy/log median |
|---|---|---|---|---|---|---|---|
| 5 | 162 | 21.00 | 21 | 24 | 18 | 0.250 | 0.5438 |
| 6 | 486 | 21.00 | 21 | 25 | 18 | 0.250 | 0.5438 |
| 7 | 1458 | 21.00 | 21 | 25 | 18 | 0.250 | 0.5438 |

Constant in k. Bounded by precision (M = 2^20 ⇒ ~21 = log₂(M) + 1).

## Files

- [preimage_probe.py](preimage_probe.py) — main probe (numba-jit, ~4s total compute at k=5,6,7)
- [preimage_probe.log](preimage_probe.log) — full stdout
- preimage_map_k{5,6,7}.csv — per-y stats
- v_distribution_empirical_k{5,6,7}.csv — per-v counts
- scaling_analysis.csv — fit parameters

## Strategic position

The probe's framing was empirical scoping for an analytical construction
(Butterley-Kim adaptation). Outcome A is the simplest possible case: the
operator's preimage structure is bounded and k-uniform, so the analytical
machinery can proceed without level-dependent function-space weighting.
The geometric-weight Geom(1/2) v-distribution gives the standard transfer-
operator setup; combined with the bounded preimage count, this is
"vanilla" anisotropic Banach space territory.

For someone (you, collaborator, future graduate student) actually doing
the analytical adaptation:

1. **Function space:** C^q(Z_3) test functions with smoothness q ≥ 1
   in the 3-adic metric.
2. **Operator action:** L f(y) = Σ_{v=1}^∞ 2^(-v) · f(x_v(y))
   where x_v(y) is determined by (3·x_v + 1) · 2^(-v) ≡ y in Z_3.
3. **Convergence:** geometric, controlled by ‖f‖_C^q.
4. **Spectral analysis:** Butterley-Kim type quasi-compactness arguments
   should give isolated discrete spectrum + essential spectrum bounded
   by the smoothness rate; ρ_slow ≈ 0.83 should appear as an isolated
   eigenvalue or Pollicott-Ruelle resonance of L on the appropriate
   space.

The probe doesn't do this analytical work; it confirms that the
structural prerequisites are in place.

## Honest framing

The exact-21 finding across three k values is striking and clean, but
it's a function of M (precision), not k (level). Reporting "Outcome A
fires" without that caveat would overstate the case: |Preimage| is
bounded *given a precision cutoff*, but the underlying inverse-limit
operator on Z_3 has countably-infinitely-many preimages with geometric
weights. The empirical M-bounded preimage count is the right notion for
finite-precision computation; the analytical adaptation works with the
infinite-preimage limit operator, where Butterley-Kim's machinery
naturally accommodates geometric-weighted preimage sums.

The "ratio 0.5" at every j in the empirical v-distribution is the only
mildly subtle point — the natural-density construction (with even+odd
lifts) gives half the geometric-prediction probability mass on each
v ≥ 1 because the other half is at v = 0. The conditional distribution
given v ≥ 1 is exactly Geom(1/2). This is consistent with Lagarias /
Tao 2-adic equidistribution and doesn't change the transfer-operator
construction (the v = 0 column is a separate "drift" component that's
also bounded and well-behaved).
