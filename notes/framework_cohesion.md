# Framework Cohesion: K_k in the Diaconis-Graham / Ayyer-Singla / Diaconis-Saloff-Coste / Eberhard-Varjú Lineage

**Author:** Nathan Humphrey
**Date:** May 6, 2026
**Project:** github.com/mrnathanhumphrey-droid/Collatz

---

## Purpose

This document records the empirically-confirmed placement of the framework's central operator K_k on (Z/3^k)* within the established analytical lineage of Markov chain spectral analysis on finite abelian groups. It documents which structural features of that lineage transfer cleanly to K_k, which quantitative tools require modification, and which closure direction the combined empirical evidence supports.

The document is dated and authored. The empirical results referenced are recorded in the framework's probe outputs at the cited paths in the Collatz repository. Together they establish the framework characterization below as publicly documented original work as of the date above.

---

## The intellectual lineage

The framework's empirical work on K_k operates within a research program with four foundational papers:

**Diaconis & Graham 1992**, "An affine walk on the hypercube" (J. Comp. Appl. Math. 41, 215-235), establishes the Fourier-analytic framework for X_{n+1} = a_n·X_n + ε_n on finite abelian groups. The key technical structure: Fourier transform of the n-fold iterated distribution Q̂_n(χ) factors as a product of per-step Fourier factors when the dynamics is iid across steps. The technique is laid out in Section 2.1.

**Ayyer & Singla 2019** (arXiv:1605.05089), "Random Motion on Finite Rings, I: Commutative Rings," extends DG into the multiplicative-with-random-coefficient case. Theorem 2.3 gives the per-character single-step factor as λ_χ = Σ_x χ(x) β_x for finite chain rings, and Corollary 4.14 establishes Jordan block structure of size k - e where e is the character conductor.

**Diaconis & Saloff-Coste 1994**, "Moderate growth and random walk on finite groups" (GAFA 4(1), 1-36), establishes sharp-rate theorems for abelian groups satisfying the moderate growth condition V(n) ≥ A·(n/γ)^d. The combination of moderate growth with eigenvalue bounds produces sharp mixing-time bounds in the form t_mix ≤ C(A,d)·γ²·log|G|.

**Eberhard & Varjú 2020** (arXiv:2003.08117), "Mixing Time of the Chung-Diaconis-Graham Random Process," extends DG into the fixed-multiplier random-shift case. The asymptotic mixing rate is governed by entropy H(a, μ) of the associated self-similar Cantor-like measure on R, formalizing the Bernoulli-convolution interpretation of CDG dynamics.

These four papers together define the analytical lineage to which the framework's methodology belongs. The empirical work documented in the Collatz repository is best understood as a specific instance of this lineage applied to the Syracuse residue dynamics.

---

## Empirically-confirmed structural fits

Three structural features of the lineage transfer cleanly to K_k, with empirical confirmation across k = 5, 6, 7.

### 1. Ayyer-Singla Jordan block prediction (probe_framework_predictions Phase 2)

AS Corollary 4.14 predicts that for a Markov chain on a finite chain ring R = Z/p^k Z, eigenvalues have algebraic multiplicity k - e where e is the character conductor, with geometric multiplicity 1. The dominant non-trivial character (e = 0) produces a Jordan block of size exactly k, capturing essentially the entire dimension.

Direct measurement on K_k:

| k | Dimension n_k | Block size at top |λ| | Dimension captured | Fraction |
|---|---|---|---|---|
| 5 | 162 | 5 | 157 | 96.9% |
| 6 | 486 | 6 | 481 | 99.0% |
| 7 | 1458 | 7 | 1444 | 99.0% |

The prediction matches at three k-values without parameter tuning. Maximum Jordan block size equals k at every level; algebraic multiplicity captures 96-99% of the operator's dimension. (Z/3^k)* is the unit group of the chain ring Z/3^k, not the chain ring itself, but the AS prediction shape transfers directly.

### 2. Diaconis-Saloff-Coste moderate growth precondition (probe_framework_test Phase 3)

DS-C 1994 requires the carrier group to satisfy V(n) ≥ A·(n/γ)^d for the sharp-rate framework to apply. Direct verification on (Z/3^k)*:

- Group structure: cyclic, order n_k = 2·3^(k-1)
- Generating set: {2, 2^(-1)} (since 2 is a primitive root mod 3^k for all k ≥ 1)
- Cayley diameter γ_k = n_k / 2
- Volume growth: d = 1, polynomial in n
- Moderate growth precondition: holds

The structural prerequisite for sharp-rate analysis is met. (Z/3^k)* is in the moderate-growth class of finite abelian groups that DS-C's framework addresses.

### 3. Cesàro-averaged framework compatibility (probe_framework_test Phase 2)

DG's iid product structure Q̂_n(χ) = F(χ)^n requires F(χ) to be approximately constant in n. Direct measurement on K_k reveals a structured pattern:

| k | Within-trajectory CoV (median) | Across-trajectory F̄(χ) variation |
|---|---|---|
| 5 | 0.284 | 0.056 |
| 6 | 0.286 | 0.073 |
| 7 | 0.275 | 0.068 |

Within-trajectory variation (~28%) far exceeds the iid product-structure tolerance. But trajectory-averaged F̄(χ) varies only ~7% across initial conditions. The Cesàro-averaged per-character factor is a meaningful invariant of K_k that is approximately initial-condition-independent.

Interpretation: K_k's near-defective Jordan structure prevents step-by-step product structure (characters aren't eigenvectors because the operator is non-diagonalizable), but ergodic averaging projects onto the dominant mode regardless of basis defect. The framework needs the iid assumption replaced by Cesàro averaging rather than abandoned.

---

## Empirically-walked-back quantitative claims

Three quantitative tools from the lineage do not transfer cleanly to K_k. The walks-back are decisive.

### 1. Ayyer-Singla character-sum eigenvalue formula (probe_ayyer_singla Outcome C)

AS gives explicit eigenvalue formula λ_χ = Σ_x χ(x) β_x for the multiplicative-only operator B_R on a finite chain ring. Applied to K_k's multiplicative-only piece B_k on (Z/3^k)*, the formula produces eigenvalues on the circle |λ - 2/3| = 1/3, so |λ_χ| ∈ [1/3, 1].

K_k's measured non-trivial eigenvalues have |λ| ~ 10^(-3) at k = 5, 6, 7. The character-sum prediction and the measured spectrum live in non-overlapping magnitude shells, separated by roughly three orders of magnitude.

The +1 additive shift in the Syracuse step is not perturbative. It suppresses the spectrum from the |λ - 2/3| = 1/3 circle to |λ| ~ 10^(-3) in a way that cannot be recovered by treating the additive shift as a small correction to the multiplicative-only operator. The character-sum formula does not predict K_k's spectrum; it predicts B_k's spectrum, which is structurally different.

### 2. Diaconis-Saloff-Coste sharp-rate Cayley-diameter bound (probe_framework_predictions Phase 3)

DS-C Theorem 1.2 gives the bound t_mix(ε) ≤ C(A,d)·γ²·log(1/ε) for moderate-growth groups. Applied to (Z/3^k)*: γ_k = n_k/2, d = 1, predicts t_mix on the order of γ_k².

Direct measurement of K_k's mixing time:

| k | n_k | γ_k = n_k/2 | DS-C bound order | Measured t_mix(1/4) | Looseness factor |
|---|---|---|---|---|---|
| 5 | 162 | 81 | ~6,500 | 5.5 | ~1,200× |
| 6 | 486 | 243 | ~59,000 | 6.5 | ~9,000× |
| 7 | 1458 | 729 | ~530,000 | 7.5 | ~70,000× |

Mixing time = Jordan block size = k + 0.5 steps at every k tested. The DS-C bound is loose by 2,000× to 100,000× and the looseness grows with k.

The bound is loose because DS-C's framework is derived for diagonalizable group walks where mixing is governed by spectral gap. K_k is Jordan-block-dominated, mixing through generalized-eigenspace structure rather than via spectral gap on a clean eigenbasis. The mixing time is set by Jordan block size, not by Cayley diameter and spectral gap.

### 3. Diaconis-Graham step-by-step product structure (probe_framework_test Phase 2)

DG's Section 2.1 establishes Q̂_n(χ) = F(χ)^n for iid dynamics. Direct measurement on K_k shows within-trajectory CoV of 27-29% in F_n(χ) magnitude across n = 2..10, with phase standard deviation of 0.75-0.94 rad.

The step-by-step product structure does not hold for K_k. Individual iterates oscillate in Fourier coordinates because characters are not eigenvectors of K_k (the operator is essentially non-diagonalizable, with eigenvector matrix condition number growing from 2.5×10^14 at k=5 to 2.4×10^17 at k=7). The product structure is recoverable only under Cesàro averaging.

---

## Synthesis: where ρ_slow lives

The framework's central convergence claim is S_k = Σ |π̂_k(χ)|^2 over a specific character orbit converges to 7/15. Direct matrix-free power iteration through k = 13 confirms the convergence at machine precision per level (ε_k tracked to ~10^(-15) absolute). The geometric **rate** ρ_slow ≈ 0.83 is identified empirically through **inverse-limit** measurement: L^1 and total-variation distances of finite-k lifts to π_∞ decay at rate ρ = 0.834 (probe_profinite, R² = 0.97 across k = 5..11).

**Walk-back (2026-05-06, ε_13 measurement):** an earlier framing identified ρ_slow with the dominant root of an order-3 linear recurrence fitted on ε_k. That fit is **not stable**: the dominant root drifts wildly with training window (0.577 at k=2..10, 0.740 at k=2..11, 1.030 at k=2..12, 1.115 at k=2..13), and the order-3 recurrence prediction for ε_13 is 9.2% off the measured value (predicted +2.677e-3, measured +2.948e-3). The 0.83-from-recurrence claim was a finite-k coincidence from one specific window before the post-zero-crossing growth (ε_10..ε_13 all positive and growing in magnitude) destabilized the fit. The reliable identification of ρ_slow is the inverse-limit L^1/TV decay rate; the recurrence root is not a structural feature of the operator. See `probe_epsilon_13/epsilon_13_findings.md` for the full diagnostic.

The combined empirical evidence places ρ_slow in a specific location relative to the lineage's machinery:

**ρ_slow is not in K_k's L^2 eigenspectrum at any tested k.** K_k's |λ_2| is approximately 10^(-3); ρ_slow ≈ 0.83 is in a non-overlapping magnitude shell. Confirmed in probe_eigenvalue_spectrum and reinforced by probe_ayyer_singla Outcome C bimodality.

**ρ_slow is not in the multiplicative-only character-sum spectrum.** That spectrum populates the circle |λ - 2/3| = 1/3 densely; many λ_χ values happen to land near magnitude 0.83 by circle artifact, but no specific character carries ρ_slow as a structural identifier. Confirmed in probe_ayyer_singla Outcome D not firing.

**ρ_slow is not in the singular spectrum of R_k.** R_k's singular values are band-supported on [0.49, 0.67], k-stable to 5 significant figures across k = 5, 6, 7. The dominant singular value σ_1 ≈ 0.6706 does not match ρ_slow. Confirmed in probe_R_operator and probe_mode_amplitudes_v2.

**ρ_slow is in the inverse-limit convergence rate.** L^1 and total-variation distances of finite-k lifts to π_∞ decay at rate 0.834. Confirmed in probe_profinite (R² = 0.97 across k = 5..11). This is the standalone reliable identification of ρ_slow; the previously-claimed "matches the order-3 recurrence's dominant root within 0.5%" framing has been walked back (see Synthesis section) — the recurrence root is window-unstable and not a structural feature of the operator.

The leading hypothesis the empirical evidence supports: ρ_slow is a Pollicott-Ruelle-type resonance of the inverse-limit transfer operator on Z_3, manifesting at finite k as slow leak through the Jordan generalized-eigenspace tower. The L^2-character basis sees only a Cesàro-smoothed shadow of this structure because K_k has no clean character-diagonalization at any finite k (cond(V) ~ 10^14 to 10^17 across k = 5, 6, 7).

This hypothesis is consistent with three independent empirical findings: K_k's essential non-diagonalizability (probe_ayyer_singla Outcome E), the Cesàro-recoverable framework compatibility (probe_framework_test Outcome B), and the mixing-time-equals-Jordan-block-size structure (probe_framework_predictions Outcome B Phase 3). It is also consistent with the entropy-deficit per-step decay rate ρ_Δ ≈ 0.85 (probe_self_similarity Phase 4) being in the same neighborhood as ρ_slow but structurally distinct, pending higher-k resolution.

The hypothesis is not proven. Constructing the inverse-limit transfer operator on Z_3 in the function-space framework where its Pollicott-Ruelle resonances become a tractable spectral object is research mathematics that requires adapting Liverani-style anisotropic Banach space machinery to the abelian profinite setting. That work is beyond the scope of empirical probes and represents the natural next analytical phase of the framework.

---

## Closure direction the evidence supports

The empirical evidence justifies a specific closure direction without proving it:

The closure of c = 7/45 lies in characterizing ρ_slow as a Pollicott-Ruelle resonance of the inverse-limit transfer operator on Z_3, accessed through Jordan-block-extended versions of the AS / DG / DS-C analytical machinery applied to a Cesàro-averaged or anisotropic-Banach-space refinement of the natural L^2 framework. The character-theoretic / Fourier-analytic / sharp-rate machinery from the foundational papers is the right intellectual home; the specific quantitative tools require modification to handle K_k's defective structure and the inverse-limit transition.

Concretely, the analytical reading direction this implies:

- Anisotropic Banach space construction for transfer operators (Liverani 1995, Faure-Sjöstrand 2011, Baladi-Tsujii 2007) adapted to discrete profinite settings
- Spectral analysis of random walks on profinite groups (Varjú 2013, Eberhard-Varjú 2020, Hermon-Olesker-Taylor 2021) with abelian-specific extensions
- Direct construction of inverse-limit operators on Z_p (Bendikov-Bobikau-Pittet and successors)
- Jordan-block-extended versions of Diaconis-Graham 1992's Fourier framework, of Ayyer-Singla 2019's chain-ring eigenvalue analysis, and of Diaconis-Saloff-Coste 1994's sharp-rate machinery

This is months-to-years of analytical work. It is not closure. But the empirical evidence has narrowed the closure question from "where does ρ_slow come from" to "what specific function space and operator construction makes ρ_slow appear as an isolated resonance of the inverse-limit transfer operator," which is a substantially more constrained question than the framework had at its empirical inception.

---

## What this document is not

This document does not claim closure of c = 7/45. The constant 7/15 (equivalently 7/45 in the framework's normalization) has been characterized empirically and structurally but not derived analytically from first principles.

This document does not claim full unification of the AS / DG / DS-C / EV machinery for K_k. The framework's structural placement is empirically confirmed; specific quantitative tools are walked back. The unification is qualified rather than complete.

This document does not assert priority claims over the cited foundational works. Diaconis, Graham, Ayyer, Singla, Saloff-Coste, Eberhard, and Varjú established the analytical lineage in which the framework's empirical work operates. The framework's contribution is the empirical placement of K_k within that lineage and the specific qualifications on which tools transfer.

---

## Provenance

Empirical results referenced are documented in the following probe directories within the Collatz repository:

- probe_eigenvalue_spectrum: K_k eigenvalue measurements at k = 5, 6, 7
- probe_R_operator: R_k singular spectrum, rank, k-stability
- probe_mode_amplitudes_v2: δ_k decomposition on R_k singular vectors
- probe_epsilon_12: order-3 recurrence verification at k = 12, residual 7.4×10^(-16)
- probe_offset_sweep: conjugation theorem K_k(c) = σ_c · K_k(1) · σ_c^(-1) at machine precision
- probe_profinite: ρ_slow as inverse-limit convergence rate, L^1 and TV match within 0.5%
- probe_smoothness: 1619 club orbital family mass concentration
- probe_ayyer_singla_test: Outcome C (character-sum spectrum decoupled from K_k spectrum) and Outcome E (cond(V) growing exponentially in k)
- probe_self_similarity: Phase 1 IFS correction, Phase 4 entropy-deficit ρ_Δ ≈ 0.85 ambiguous between same-mode and distinct-modes pending higher-k
- probe_framework_test: Outcome B (Cesàro-averaged framework compatibility, moderate growth confirmed)
- probe_framework_predictions: Outcome B (Jordan block size = k captures 99% dimension; mixing time = block size; per-character bimodal partition)

---

— Nathan Humphrey, Annapolis, MD
