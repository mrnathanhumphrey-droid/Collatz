# Edgeworth Standardization Test (Result 30)

**Question.** Does choice of standardization (raw σ vs σ−K_h·log_n) explain the
N-growing discrepancy between empirical C and σ_V·κ_111/2 in the bivariate
Edgeworth conditional-mean correction for E_band(q)?

**Answer. NO. Outcome (c).** Standardization is not the issue. The discrepancy
is a structural mismatch between the Gram-Charlier coefficient σ_V·κ_111/2
(which vanishes as N→∞) and the empirical correction C (which plateaus at
~0.21·σ_V).

## Method

Three standardizations of σ tested at N ∈ {2³², 2³⁴, 2³⁶, 2³⁸}, 500k orbits each:

- **A**: Z_σ = (σ − E[σ])/SD[σ] using raw σ
- **B**: Z_σ = (σ − K_h·log n − c)/SD[σ_resid], **theoretical** K_h = 3/log(4/3) ≈ 10.4282
- **B'**: Z_σ = (σ − β_emp·log n − c)/SD[σ_resid], **empirical** β regressed on data

For each, compute:
- ρ = corr(V, Z_σ)
- κ_111 = E[Z_V · Z_σ²]
- C_pred = σ_V · κ_111 / 2 (Gram-Charlier)
- C_emp = OLS slope of corr_emp = E[V|band]_emp − (μ_V + ρσ_V·E[Z|band]) regressed on E[Z²−1|band]

## Why standardization makes no difference

| log2N | Var[log n] | Var[σ_resid] | drift contribution |
|-------|-----------:|-------------:|-------------------:|
| 32    | 1.003      | 5,425        | 0.019%             |
| 34    | 1.004      | 5,823        | 0.018%             |
| 36    | 1.004      | 6,202        | 0.017%             |
| 38    | 1.004      | 6,569        | 0.017%             |

For uniform-on-[1,N] sampling, Var[log n] → 1 (asymptotic). The K_h·log n drift
contributes <0.02% of σ-variance. Standardization A and B are numerically
indistinguishable.

## Results — three standardizations across N

| log2N | std | ρ | κ_111 | C_pred | C_emp | ratio | R² |
|-----:|:--:|------:|--------:|---------:|---------:|------:|------:|
| 32 | A  | −0.8458 | 0.05408 | 0.00596 | 0.04730 |  7.94 | 0.869 |
| 32 | B  | −0.8516 | 0.05268 | 0.00580 | 0.04808 |  8.28 | 0.863 |
| 32 | B' | −0.8516 | 0.05268 | 0.00580 | 0.04812 |  8.29 | 0.863 |
| 34 | A  | −0.8542 | 0.04172 | 0.00436 | 0.04375 | 10.04 | 0.863 |
| 34 | B  | −0.8598 | 0.03980 | 0.00416 | 0.04486 | 10.79 | 0.848 |
| 34 | B' | −0.8597 | 0.03978 | 0.00415 | 0.04481 | 10.79 | 0.848 |
| 36 | A  | −0.8623 | 0.03252 | 0.00323 | 0.04104 | 12.69 | 0.869 |
| 36 | B  | −0.8679 | 0.03150 | 0.00313 | 0.04166 | 13.30 | 0.866 |
| 36 | B' | −0.8679 | 0.03150 | 0.00313 | 0.04167 | 13.30 | 0.866 |
| 38 | A  | −0.8690 | 0.02785 | 0.00265 | 0.03880 | 14.63 | 0.869 |
| 38 | B  | −0.8743 | 0.02598 | 0.00247 | 0.03927 | 15.88 | 0.857 |
| 38 | B' | −0.8743 | 0.02598 | 0.00247 | 0.03928 | 15.88 | 0.857 |

CV(ratio_B across N) = 0.234 — **N-varying, not N-stable.** Outcome (c).

ratio_B ≈ 1.265·log2N − 32.20 (R² ≈ 1.00). Linear divergence in log2N.

## What's actually wrong

| quantity | log-log slope vs N |
|----------|-------------------:|
| C_emp     | N^(−0.142) |
| σ_V       | N^(−0.102) |
| skew[V]   | N^(−0.109) |
| κ_111     | **N^(−0.490)** |
| σ_V·κ_111/2 | N^(−0.592) |

κ_111 decays **3.5× faster** than C_emp. Asymptotically:
- κ_111 → 0 (joint distribution becomes Gaussian in standardized third order)
- C_emp → ~0.21·σ_V_∞ > 0 (persistent)

The Gram-Charlier coefficient σ_V·κ_111/2 captures a **vanishing** subleading
piece. The persistent piece is driven by the persistent marginal V skew
(skew[V] ≈ 1.75–2.04, decays at the same rate as σ_V).

## Closer-to-stable rescalings

| log2N | C_emp/σ_V | C_emp/(σ_V·skew_V/2) |
|------:|----------:|---------------------:|
| 32 | 0.2181 | 0.2138 |
| 34 | 0.2147 | 0.2203 |
| 36 | 0.2095 | 0.2302 |
| 38 | 0.2062 | 0.2353 |

C_emp/σ_V drifts only 5.5% across factor 64 in N (vs ratio drifting 92%).
**C ≈ 0.21·σ_V** is the empirical structural pattern, not the Gram-Charlier form.

## Why the bivariate Edgeworth fails here

ρ ≈ −0.87 puts the joint distribution near the rank-degenerate limit
(1−ρ² = 0.24). Standard bivariate Gram-Charlier expansions converge poorly
here: each correction term carries factors of (1−ρ²)^(−k) for various k, and
the truncation order needed grows.

Marginal V skew is large (≈ 2) and persistent. Conditional mean E[V|σ-band]
inherits asymmetry from V's marginal, not just from joint cross-cumulants
captured by κ_111.

In the asymptotic limit:
- σ standardization → Gaussian (CLT)
- κ_111 (mixed standardized cross-cumulant) → 0
- Marginal V remains non-Gaussian (Geom-skew from per-step v_t structure)
- E[V|σ-band] retains a persistent O(σ_V) shape correction

The shape function E[Z²−1|band] is the right basis (R² = 0.87) but the
coefficient must come from V's marginal structure, not from a vanishing
mixed cumulant.

## Decisive negative result for v3.5+

C does **not** close as σ_V·κ_111_B/2 under any standardization (A, B, or B').
Outcome (c) confirmed: full Stuart-Ord bivariate derivation with multiple
cumulants, OR a non-Gaussian-base expansion (saddlepoint / direct
Esscher-tilted conditional), is required.

The next-cheapest test for closing C: try σ_V · skew[V] · g(ρ) / k for some
combination, or derive the conditional mean directly from the joint Geom
structure (since E[v|band] is now well-characterized via per-band Esscher
inversion in Result 25).

## Files

- `experiments/63_edgeworth_standardization_test.py`
- `experiments_output/63_edgeworth_standardization_test.csv`
- `experiments_output/63_edgeworth_standardization_test_log.txt`
