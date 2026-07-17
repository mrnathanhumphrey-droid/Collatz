# Hateley Framework — Duality with the Syracuse Residue Framework Eliminated

**Date:** 2026-05-05.

## One-line result

> The framework's Syracuse residue-Markov structure on Z/3^k coprime classes (carrying S_k → 7/15) and Hateley's transfer-operator framework on N (carrying block-decay c_j ~ C·6^(-j)) describe **non-overlapping structural information** about Collatz dynamics. No natural transformation in either direction maps one framework's content onto the other's. Duality is eliminated.

## Setup

Two frameworks are compared:

- **Framework (this project).** Syracuse Markov chain on coprime classes of Z/3^k, with stationary distribution π_k. The framework's central object is the Plancherel mass S_k of π_k at non-trivial mod-3 characters; conjectured limit S_∞ = 7/15, and the rate-1/2 envelope of |ε_n| := |S_n − 7/15| is the load-bearing analytic content.
- **Hateley 2026** ("The Collatz Conjecture and the Spectral Calculus for Arithmetic Dynamics", preprint posted 2026-01-28). A backward transfer-operator framework on N with a multiscale Banach space and a predicted block-average asymptotic c_j = (1/|I_j|) · Σ_{n ∈ I_j} h(n) ~ C·6^(-j) on the dyadic-triadic blocks I_j = [6^j, 2·6^j).

The duality question: is one framework's structural content a corollary of the other's via a natural map between their spaces?

## Two-direction empirical test

**Direction 1: Framework → block decay (lift via 1/n weight).**
Define h_lift_k(n) = π_k(n mod 3^k) / n on coprime n, zero otherwise. Compute c_j(h_lift_k) on N = 6^8, j = 0..7. Fit log(c_j) vs j.

| lift | slope | (slope − target)/|target| | R² | C_fit | C_predicted = ln(2)/3^k | C_fit/C_pred |
|---|---:|---:|---:|---:|---:|---:|
| uniform baseline | -1.856 | -3.57% | 0.998 | 0.636 | 0.462 | 1.376 |
| k=1 | -1.821 | -1.61% | 0.999 | 0.266 | 0.231 | 1.153 |
| k=2 | -1.839 | -2.61% | 0.999 | 0.0975 | 0.0770 | 1.266 |
| k=3 | -1.800 | -0.43% | 0.998 | 0.0264 | 0.0257 | 1.027 |
| k=4 | -1.806 | -0.78% | 0.999 | 0.00909 | 0.00857 | 1.061 |
| k=5 | -1.803 | -0.63% | 0.998 | 0.00294 | 0.00286 | 1.030 |

The −ln(6) slope is reproduced under any 1/n-weighted lift, including the *uniform-residue baseline* (no π_k information at all). The 6^(−j) decay is a generic property of harmonic-sum integration over dyadic-triadic blocks, not a structural consequence of π_k. The prefactor C_lift_k = ln(2)/3^k follows from Σ_r π_k(r) = 1 alone — the fine structure of π_k does not enter.

**Direction 2: block decay → framework (project to coprime classes mod 3^k).**
Define Π_k(h)(r) = (1/Z) · Σ_{n ≡ r mod 3^k, coprime} h(n) · w(n) for natural weights w. For any block-uniform asymptotic h ~ C/n and any weight w(n) depending only on n's value, the projection onto Z/3^k coprime classes is asymptotically uniform (each class has density 1/(c·3^k); harmonic-sum contributions equalize across r), giving:

| k | S~_k (uniform projection) | S_k (framework) |
|---:|---:|---:|
| 1 | 1/2 = 0.500 | 2/3 = 0.667 |
| 2 | 0 | 10/21 = 0.476 |
| 3 | 0 | 0.4621 |
| 4 | 0 | 0.4639 |
| 5 | 0 | 0.4673 |

(For k ≥ 2 the Plancherel mass at non-trivial mod-3 characters of a uniform measure on coprime classes vanishes by the Ramanujan-sum identity μ(3^k) = 0 for k ≥ 2.) Tested with the Syracuse-acceleration weight w(n) = 1 + v_2(3n+1) as well — same conclusion, since v_2(3n+1) depends on n mod 2-power, which is independent of n mod 3^k by CRT.

## Reading

The two frameworks are **structurally orthogonal**:

- The framework's information content lives in residue distributions on Z/3^k coprime classes — the non-trivial mod-3 character mass, S_k.
- Hateley's framework's information content lives in block averages on the 6-adic block decomposition I_j.

Each framework's natural map onto the other's space erases its own distinctive content. The framework's S_k structure is *not visible* in any block average. Hateley's block decay is *not constrained* by any residue-class measure on Z/3^k.

## Conclusion

The Syracuse residue framework's c = 7/45 work and Hateley 2026's framework address different structural questions about Collatz dynamics. Neither implies the other; the two are independent contributions. The framework's rate-1/2 envelope and 7/45 closure are its own results, with no dual derivation through Hateley's spectral apparatus.

## Files

- `result_hateley_duality_eliminated.md` — this writeup (the only Hateley-related artifact in the repository).
