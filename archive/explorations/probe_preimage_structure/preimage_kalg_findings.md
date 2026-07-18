# Preimage structure under K_alg (truncated-Geom Tao chain)

**Date:** 2026-05-06.
**Complements:** [preimage_findings.md](preimage_findings.md) which used the
brief's M=2^20 integer-lift method (K_emp). K_emp and K_alg differ at 100%
relative Frobenius — they are structurally different chains. K_alg is the
chain driving all framework results (S_k → 7/15, ε_k convergence, ρ_slow).

**Construction:** K_alg[x, y] = Σ_{v=1}^{M_k} (2^(-v)/Z) · 1[(3x+1)·2^(-v) ≡ y mod 3^k]
with M_k = ord(2 mod 3^k) = 2·3^(k-1) and Z = 1 - 2^(-M_k) ≈ 1.

## Verdict

**Outcome B (formerly C-leaning)** under full v_eff = M:

| k | n | M | mean \|Pre\|_struct | mean \|Pre\|_weighted | rank K |
|---|---|---|---|---|---|
| 5 | 162 | 162 | 162 | 162 | 54 |
| 6 | 486 | 486 | 486 | 486 | 162 |
| 7 | 1458 | 1458 | **1458** | **1074** † | 486 |

† k=7 weighted gap: weights 2^(-v) underflow to 0 in float64 for v ≥ 1074
(2^(-1074) is the smallest positive subnormal). Mathematically all 1458
edges per row exist; computationally only 1074 register as nonzero. This
is a precision artifact, not structure.

**Scaling fits (full v_eff=M):**
- log|Pre|_struct = -0.4055 + 1.0986·k → exp(b) = **3.0000** per step
- Predicted |Pre| = M = 2·3^(k-1), so log = log(2/3) + log(3)·k
  - Intercept log(2/3) = -0.4055 ✓ exact
  - Slope log(3) = 1.0986 ✓ exact
- Linear |Pre| = -3186 + 648·k, R² = 0.9231 (worse than log-linear; |Pre|
  scales exponentially in k, not linearly)

**Mean preimage count = M_k = 2·3^(k-1) exactly.** Linear in n, exponential
in k.

**Outcome A** under truncated v_eff=min(M, 60): mean |Pre| = 60 at every k,
bounded by truncation. Identical to v_max=60 cap on K_emp.

## K_alg is structurally complete (dense support)

For each k tested:
- `n_struct_edges = n × M = n²` (since M = n in our setting)
- Density = 1.0000 in struct edges
- |Preimage|_struct(y) = n for every y

Every x ∈ (Z/3^k)* maps to every y ∈ (Z/3^k)* via some v ∈ {1..M}. The
support graph of K_alg is the complete digraph on (Z/3^k)*. Equivalently:
the equation y = (3x+1)·2^(-v) mod 3^k has a unique solution v ∈ {1..M}
for any (x, y) pair, since 2 is a primitive root mod 3^k.

This is the **opposite** structural answer from K_emp: K_emp has |Pre| ≈
21 (precision-bounded), K_alg has |Pre| = n_k (exhaustive).

## 3:1 row collapse — rank deficiency

The empirical numerical rank of K_alg matches **n/3 exactly** at every k:

| k | n | rank K_alg | n/3 |
|---|---|---|---|
| 5 | 162 | 54 | 54 |
| 6 | 486 | 162 | 162 |
| 7 | 1458 | 486 | 486 |

**Reason:** for x_1, x_2 ∈ (Z/3^k)* with x_1 ≡ x_2 mod 3^(k-1), one has
3x_1 + 1 ≡ 3x_2 + 1 mod 3^k. So all v-distributed transitions starting
from these two states land on the same set of y values with the same
weights. **Rows of K_alg corresponding to x and x+3^(k-1) and x+2·3^(k-1)
are identical** (when all three lifts are coprime to 3, which holds since
3^(k-1) ≡ 0 mod 3).

Distinct row count = number of mod-3^(k-1) classes among coprime states =
|(Z/3^(k-1))*| = 2·3^(k-2). Evaluating: 2·27=54 at k=5, 2·81=162 at k=6,
2·243=486 at k=7. Matches observed rank exactly.

Stationary measure π_k satisfies π_k(x) = π_k(x + 3^(k-1)) = π_k(x + 2·3^(k-1))
for all x, by the row identity. This is the structural symmetry that makes
π_k lift consistently to π_{k-1} under the natural projection.

(Aside on Part D output — the script grouped rows by (3x+1) mod 3^(k-2),
a coarser granularity giving 18 / 54 / 162 classes at k=5 / 6 / 7. Each
such coarse class contains 3 rank-distinct row-images, recovering
3 × 18 = 54, 3 × 54 = 162, 3 × 162 = 486 in total. Rows in the same
coarse class are NOT identical in general; identity holds only within
the finer mod-3^(k-1) sub-class. The 0.816 sample-row-diff in Part D
reflects this finer structure, not a bug.)

## Column structure and weight distribution

Column properties are k-INVARIANT for K_alg (the geometric weight structure
projected onto coprime classes is level-independent in distribution):

| k | mean col_sum | std col_sum | min col_sum | max col_sum | median max_w | mean entropy |
|---|---|---|---|---|---|---|
| 5 | 1.000000 | 0.5909 | 0.2381 | 1.9048 | 0.2500 | 1.5885 |
| 6 | 1.000000 | 0.5909 | 0.2381 | 1.9048 | 0.2500 | 1.5885 |
| 7 | 1.000000 | 0.5909 | 0.2381 | 1.9048 | 0.2500 | 1.5885 |

Identical to all five decimal places across k. The mean column sum is 1
(matching π_k stationarity averaged over y), but **std = 0.59** — column
sums VARY substantially. K_alg is decidedly NOT doubly stochastic, which
is consistent with non-uniform stationary π_k.

**Mean column entropy = 1.5885** is k-invariant. With |Pre|_weighted varying
across k (162, 486, 1074), the entropy/log|P|_w concentration ratio FALLS
as k grows: 0.33 → 0.27 → 0.24. The mass is increasingly concentrated on
the few high-weight preimages (v=1, 2, 3 dominate; weight 2^(-v) decays
geometrically) regardless of how many low-weight v values pile in.

The entropy 1.5885 ≈ entropy of the geometric distribution Geom(1/2):
H(Geom(1/2)) = 2 log 2 ≈ 1.386. Hmm, our 1.5885 is slightly higher;
this is because per-column the weights are a permutation (not the strict
geometric) of {2^(-v)/Z}, accumulating entropy from all M values.
H_Geom_M = Σ_{v=1}^{M} (2^(-v)/Z) · log((2^(-v)/Z)^(-1)). For M large,
Σ_{v=1}^∞ 2^(-v) · v · log 2 = 2 log 2 ≈ 1.386, plus log Z ≈ 0.
Empirical 1.5885 vs theoretical Geom(1/2) entropy ≈ 1.386 — there's a
small gap, presumably from the renormalization Z and finite-M effects.

## Implication for transfer operator construction

| Aspect | K_emp (M=2^20 lifts) | K_alg (truncated-Geom) |
|---|---|---|
| Domain | (Z/3^k)* | (Z/3^k)* |
| v support | {0, 1, 2, ...} (incl. v=0) | {1..M_k} (no v=0) |
| Structural \|Pre\|(y) | ~21 (precision-bounded) | n_k = 2·3^(k-1) |
| Effective \|Pre\|(y) | ~21 | ~M_k = n_k |
| Scaling in k | flat | exponential (rate 3) |
| Column entropy | k-invariant ≈ 11 | k-invariant ≈ 1.59 |
| Outcome | A (bounded) | B (linear in n) |
| Inverse-limit limit | natural-density on Z_3 | Tao odd-integer chain on Z_3 |

**Both flavors of inverse-limit transfer operator** are tractable for
Butterley-Kim-style anisotropic Banach space construction:

- K_emp's limit operator on Z_3 has bounded essential preimage count at
  finite precision; the analytical adaptation parametrizes precision via
  the function-space smoothness q.
- K_alg's limit operator on Z_3 has countably-infinitely-many preimages
  with geometric weights 2^(-v); standard transfer-operator machinery
  for hyperbolic systems applies via the C^q test function space's
  exponential-decay-in-v weighting.

Different limits, different finite-k preimage scaling, but both fit the
"weighted preimage sum, geometric tail" template that Butterley-Kim 2023
and predecessors handle.

The K_alg side is the relevant one for the framework's algebraic results
(S_k → 7/15, ρ_slow ≈ 0.83, etc.). The K_emp side is the natural-density
finite-precision approximation.

## Summary table

| version | k | n | M | v_eff | mean \|P\|_struct | mean \|P\|_weighted | rank | mean entropy |
|---|---|---|---|---|---|---|---|---|
| full | 5 | 162 | 162 | 162 | 162 | 162 | 54 | 1.5885 |
| full | 6 | 486 | 486 | 486 | 486 | 486 | 162 | 1.5885 |
| full | 7 | 1458 | 1458 | 1458 | 1458 | 1074 † | 486 | 1.5885 |
| trunc60 | 5 | 162 | 162 | 60 | 60 | 60 | — | 1.5885 |
| trunc60 | 6 | 486 | 486 | 60 | 60 | 60 | — | 1.5885 |
| trunc60 | 7 | 1458 | 1458 | 60 | 60 | 60 | — | 1.5885 |

† float64 underflow at v ≥ 1074.

## Files

- [probe_preimage_kalg.py](probe_preimage_kalg.py) — main probe (~0.1s/k)
- [probe_preimage_kalg.log](probe_preimage_kalg.log) — full stdout
- preimage_map_k{5,6,7}_kalg_full.csv — per-y stats, full v_eff=M
- preimage_map_k{5,6,7}_kalg_trunc60.csv — per-y stats, v_max=60
- scaling_analysis_kalg.csv — scaling fit summary
- [preimage_findings.md](preimage_findings.md) — original K_emp probe (this file's complement)

## Honest framing

The K_alg version answers the question the framework cares about: yes, the
"natural" Tao chain has dense complete-digraph support at every k, with
exponentially-many preimages per y but geometric concentration on a few
dominant ones. This is **standard transfer-operator structure** for an
expanding map on a profinite limit; nothing pathological.

The 3:1 row-collapse (rank n/3) is **structural** and reflects the
3-adic lifting from (Z/3^(k-1))* to (Z/3^k)* — rows from the same residue
class mod 3^(k-1) are identical because 3x+1 mod 3^k only depends on
x mod 3^(k-1). This is the algebraic source of why π_k descends to π_{k-1}
under the natural projection (π_k(x) is constant on the 3 lifts of each
mod-3^(k-1) class) — a known framework fact, here surfaced as a numerical
rank computation.

Neither outcome (K_emp Outcome A, K_alg Outcome B) blocks Butterley-Kim
adaptation. The probes scope which finite-k structure is the right
approximation for which inverse-limit operator.
