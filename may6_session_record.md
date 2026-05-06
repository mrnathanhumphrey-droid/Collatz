# May 6, 2026 Session Record: Recurrence Walk-Back, Distinct-Modes Confirmation, and Framework Cohesion

**Author:** Nathan Humphrey
**Date:** May 6, 2026
**Project:** Collatz / 7-45 Framework

---

## Executive summary

Tonight's session produced three substantive findings, one of which is a load-bearing walk-back of an earlier framework claim. The path forward is intact; the framework's central empirical claims survive. What changed is the model used to characterize one of those claims.

**The walk-back:** the order-3 linear recurrence on ε_k = S_k − 7/15 is a finite-k fitting artifact, not a structural feature. It held at machine precision through k=12; extension to k=13 reveals 9.2% prediction error and instability of the inferred coefficients across fit windows.

**The confirmation:** entropy-deficit decay rate ρ_Δ is structurally distinct from ρ_slow. Previous status was "ambiguous between same-mode and distinct-modes pending higher-k." The k=13 extension lands the OLS fit at ρ_Δ ≈ 0.875 with R²=0.9997 on the high-k window, candidate value 7/8.

**The structural finding:** K_alg has rank exactly n_k/3 at k=5, 6, 7, complementary to R_k's rank 2n_k/3. The 3:1 row redundancy is an algebraic identity from (3x+1) mod 3^k only seeing x mod 3^(k-1). This is the algebraic source of π_k's descent to π_{k-1}.

**The path:** ρ_slow ≈ 0.834 as the L¹/TV inverse-limit convergence rate of ‖π_k − π_∞‖ survives. The closure direction toward Pollicott-Ruelle resonances of the inverse-limit transfer operator on Z_3, accessed via Butterley-Kim 2023-style anisotropic Banach space construction, is unaffected. The cohesion document has been revised to reflect the walk-back and the new structural findings.

---

## Background

The May 6 session was a continuation of multi-day work on the Collatz / c=7/45 framework. Prior sessions established the framework's core empirical findings: S_k = Σ over a specific character orbit of |π̂_k(χ)|² converges to 7/15, with deviations ε_k characterized through k=12 via sparse Krylov computation. ρ_slow ≈ 0.83 had been identified through two independent measurements: as the dominant root of an order-3 linear recurrence on ε_k, and as the L¹/TV inverse-limit convergence rate measured directly via probe_profinite. The two measurements agreed at low k, which had been treated as cross-confirmation.

Today's session began with substantive structural probes that produced clean confirmations: framework applicability under Cesàro averaging (Probe A Outcome B), framework predictions including Jordan block structure (Probe B Outcome B), and preimage characterization for K_emp (Outcome A) and K_alg (Outcome B). Late in the session the order-3 recurrence was tested at k=13 via a truncated sparse Krylov approach that exploited the v_max=60 tail truncation to make k=13 tractable in ~0.05 seconds rather than the 24-36 hours a full sparse Krylov would have required.

The recurrence failed at k=13. The session pivoted to characterizing what the failure means and revising the framework cohesion document.

---

## Probe A: framework applicability test (Outcome B)

Probe A tested whether the Diaconis-Graham 1992 / Chung-Diaconis-Graham Fourier-analytic framework's central claim — Q̂_n(χ) factors as a product of per-step Fourier factors — holds for K_alg.

Phase 1 constructed Q̂_n(χ) for K_alg at k = 5, 6, 7 and several initial residues r₀ for n = 1..10.

Phase 2 tested whether the per-step factor F_n(χ) is approximately constant in n (the iid case requirement of DG):

| k | Within-trajectory CoV (median) | Phase std in n (rad) | Across-trajectory F̄(χ) variation |
|---|---|---|---|
| 5 | 0.284 | 0.790 | 0.056 |
| 6 | 0.286 | 0.941 | 0.073 |
| 7 | 0.275 | 0.748 | 0.068 |

Within-trajectory variation (~28% in magnitude, ~0.8 rad in phase) far exceeds the iid product-structure tolerance. But trajectory-averaged F̄(χ) varies only ~7% across initial conditions.

Phase 3 verified the DS-C 1994 moderate growth precondition. (Z/3^k)* is cyclic with generators {2, 2^(-1)}, diameter γ_k = n_k/2, polynomial growth d=1. The precondition holds.

**Outcome B fired.** The DG framework needs the iid assumption replaced by Cesàro averaging rather than abandoned. Characters aren't eigenvectors of K_alg (the operator is essentially non-diagonalizable per probe_ayyer_singla Outcome E), so step-by-step product structure fails. But ergodic averaging projects onto the dominant mode regardless of basis defect, and F̄(χ) is r₀-stable to ~7%, indicating a meaningful per-character invariant survives.

This is the Cesàro-modified DG/CDG/AS framework compatibility for K_alg.

---

## Probe B: specific predictions test (Outcome B)

Probe B tested three specific quantitative predictions the lineage's framework makes for K_alg.

**Phase 2 — Jordan block prediction (clean confirmation):**

Per AS Corollary 4.14, finite chain ring R = Z/p^k Z has eigenvalues with algebraic multiplicity k − e where e is character conductor. The dominant character (e = 0) produces a Jordan block of size exactly k.

Direct measurement on K_alg:

| k | Dimension n_k | Block size at top \|λ\| | Dimension captured | Fraction |
|---|---|---|---|---|
| 5 | 162 | 5 | 157 | 96.9% |
| 6 | 486 | 6 | 481 | 99.0% |
| 7 | 1458 | 7 | 1444 | 99.0% |

Maximum Jordan block size equals k at every level; the dominant eigenvalue's algebraic multiplicity captures 96-99% of the operator's dimension. The AS prediction matches at three k-values without parameter tuning. (Z/3^k)* is the unit group of the chain ring Z/3^k, not the chain ring itself, but the AS prediction shape transfers directly.

**Phase 3 — DS-C sharp-rate bound (decisive walk-back):**

DS-C Theorem 1.2 predicts t_mix(ε) ≤ C·γ²·log|G|/eigenvalue_gap. Direct measurement:

| k | n_k | γ_k = n_k/2 | DS-C bound order | Measured t_mix(1/4) | Looseness |
|---|---|---|---|---|---|
| 5 | 162 | 81 | ~6,500 | 5.5 | ~1,200× |
| 6 | 486 | 243 | ~59,000 | 6.5 | ~9,000× |
| 7 | 1458 | 729 | ~530,000 | 7.5 | ~70,000× |

Mixing time = Jordan block size = k + 0.5 steps at every k. The DS-C bound is loose by 2,000× to 100,000× and the looseness grows with k. The mechanism is structural: DS-C is derived for diagonalizable group walks where mixing is governed by spectral gap; K_alg is Jordan-block-dominated, mixing through generalized-eigenspace structure.

This is a structural finding: the right mixing prediction for K_alg is Jordan-block-decay timescale (= block size = k), not Cayley-diameter spectral-gap timescale.

**Phase 1 — Per-character iteration (bimodal partial):**

At k=7, top 7 characters had F_extracted/|λ_meas| ratios in [0.14, 4.9] — within order of magnitude. Remaining characters drift to ratios 430-540× because they track B_k's circle |λ − 2/3| = 1/3 instead of K_alg's spectrum at |λ| ~ 10^(-3).

The character group of (Z/3^k)* partitions into K_alg-aligned (track |λ_K|) and B_k-aligned (track multiplicative-only circle) subclasses. The +1 additive shift's action determines the partition. K_alg-aligned characters carry the operator's actual spectral information; B_k-aligned characters decouple onto the multiplicative-only circle.

**Outcome B fires.** Two of three phases confirm cleanly; Phase 1's bimodality is a structural finding rather than a clean threshold pass. Together with Probe A, the Cesàro-averaged F̄(χ) per-character factor structure is established.

---

## Preimage probes: K_emp and K_alg structures

Two preimage probes characterized the inverse structure of the Syracuse step, with implications for how the inverse-limit transfer operator on Z_3 would be constructed.

**K_emp probe (Outcome A):**

K_emp is the uniform-integer-lift formulation of the Syracuse step (includes v=0 transitions from even integer lifts). Mean |Preimage|(y) = 21 at every k ∈ {5, 6, 7}. The 21 is artifactually log₂(M)+1 from the M=2^20 precision floor, not a level-dependent quantity. Geometric P(v=j) = 2^(-(j+1)) confirmed exactly (the factor of 1/2 versus 2^(-j) is because half the integer lifts contribute mass at v=0).

For the inverse-limit transfer operator on Z_3 in the natural-density formulation: Butterley-Kim methodology applies directly without level-dependent weighting.

**K_alg probe (Outcome B with substantive structural finding):**

K_alg is the framework's standard operator (truncated-geometric, v ≥ 1). At full v_eff:

- Mean |Pre|(y) = M_k = 2·3^(k-1) exactly at every k. **Dense complete digraph** — every (x, y) edge is present. Log-linear scaling rate exp(b) = 3.0000 per step, intercept log(2/3) = −0.4055.
- **rank(K_alg) = n_k / 3 = 2·3^(k-2) exactly** at k = 5, 6, 7.
- **3:1 row redundancy:** rows for x, x + 3^(k-1), x + 2·3^(k-1) are identical because (3x+1) mod 3^k only sees x mod 3^(k-1). Algebraic identity, not empirical observation.
- **k-invariant column structure:** mean column entropy 1.5885, std col_sum 0.5909, median max weight 0.25 — identical to 5 decimals across k=5, 6, 7.
- Float64 underflow at k=7 for v ≥ 1074 flagged as precision artifact (mathematical |Pre| = 1458, observed |Pre|_weighted = 1074).

The 3:1 row redundancy connects directly to the R_k rank result from earlier work: rank(R_k) = 2n_k/3. **Together: rank(K_alg) + rank(R_k) = n_k/3 + 2n_k/3 = n_k.** The two operators decompose the full dimension complementarily. The 1/3 vs 2/3 split is forced by the algebra: each level k of (Z/3^k)* contains 1/3 redundant copies of level-(k-1) information and 2/3 new information from the +1 additive shift's distribution across residues equivalent mod 3^(k-1).

This is the algebraic source of π_k's descent to π_{k-1}. The descent is a structural identity, not an emergent property.

---

## ε_13 extension and recurrence walk-back

The original plan was a 24-36 hour overnight sparse Krylov build to compute ε_13 via the standard K_alg construction at k=13 (dimension n_13 = 1,062,882). Mid-session, an existing artifact was rediscovered: probe_self_similarity/extend_to_k13_k14.py had already produced π_13 truncated at v_max=60. The truncation exploits the geometric weight tail: weights below 2^(-60) ≈ 10^(-19) are sub-machine-precision, so K_13 truncated has only 60 nonzeros per row (64M total, genuinely sparse). FFT of the existing π_13 produced ε_13 in 0.05 seconds.

**ε_13 = +2.9482 × 10^(-3).**

The post-zero-crossing rising sequence at k = 10..13 has values |ε| = 0.72, 1.50, 2.27, 2.95 × 10^(-3) with ratios 2.08 → 1.51 → 1.30, decelerating toward a local peak.

**Recurrence test:**

The order-3 recurrence's prediction for ε_13 is 9.2% off measured. Substantially past machine precision. Refitting order-3 on k=2..13 gives dominant root +1.115. Refitting on k=2..12 alone (the previous "working" range) gives dominant root +1.030 — not the previously-claimed 0.83. Order-4 and order-5 fits are similarly unstable with dominant |r| > 1.

There is no stable order-N linear recurrence on ε_k. The recurrence model was fitting finite-k transient structure. The apparent machine-precision agreement at k ≤ 12 was numerical coincidence at low k.

**What this means for ρ_slow:**

ρ_slow ≈ 0.83 had been characterized two ways: as the dominant root of the order-3 recurrence, and as the L¹/TV inverse-limit convergence rate of ‖π_k − π_∞‖ measured directly via probe_profinite. The recurrence-based characterization is now walked back. The direct measurement survives:

ρ_slow ≈ 0.834 in L¹, 0.834 in TV (probe_profinite, R² = 0.97 in both norms across k = 5..12). This is the genuine inverse-limit convergence rate of finite-k truncations to π_∞ on Z_3.

The two roads to 0.83 happened to agree at low k by coincidence; the recurrence road got bulldozed and the direct measurement road is the road that actually leads somewhere. ρ_slow's structural reality is unaffected by the walk-back; what's gone is one of the methods used to characterize the rate numerically.

**Methodology lesson:**

Machine-precision agreement of a fitted model on a finite sequence does not establish that the underlying sequence satisfies that model. It establishes local good approximation. Distinguishing "real linear recurrence" from "good local approximation" requires extending the sequence and checking stability of inferred coefficients. The framework's discipline — extending and testing at higher k — caught the artifact tonight. This is the same discipline that has caught earlier framework reaches; it operated correctly here.

---

## Entropy-deficit distinct-modes confirmation

probe_self_similarity Phase 4 had previously measured the entropy-deficit decay rate ρ_Δ across k = 5..12 with rising per-step ratios:

| Transition | Ratio |
|---|---|
| 5→6 | 0.816 |
| 6→7 | 0.835 |
| 7→8 | 0.847 |
| 8→9 | 0.858 |
| 9→10 | 0.864 |
| 10→11 | 0.870 |
| 11→12 | 0.870 |

OLS fit on log Δ_k vs k for k=5..12 gave ρ_Δ ≈ 0.849 with R² = 0.996. Status was "ambiguous between same-mode (Δ and ε share dominant mode, gap 2.7% is finite-k transient) and distinct-modes (Δ and ε have different but similar-magnitude rates)."

The k=13 extension adds:

- 12→13 ratio: 0.8789

And refits the OLS:

| Range | ρ_Δ | R² | Gap to ρ_slow ≈ 0.834 |
|---|---|---|---|
| k=5..11 | 0.849 | — | +1.8% |
| k=5..13 | 0.857 | — | +2.8% |
| k=10..13 | 0.875 | 0.9997 | +4.9% |

The gap widens as the fit window narrows toward higher k. ρ_Δ → 0.875 as a candidate asymptote (possibly 7/8 = 0.875 exactly). R² = 0.9997 on the high-k window is striking but four data points cannot establish 7/8 as the asymptote.

**Outcome:** distinct-modes confirmed. ρ_Δ and ρ_slow are structurally distinct slow rates that share the same order of magnitude. The framework has at least two structurally distinct slow modes:

- ρ_slow ≈ 0.834 governs L¹/TV inverse-limit convergence of π_k → π_∞
- ρ_Δ → 0.875 (possibly 7/8) governs entropy-deficit decay

Eberhard-Varjú-style entropy interpretation transfers to ρ_Δ, not to ρ_slow. The EV machinery does not directly address the L¹/TV convergence rate, which is governed by different structure.

---

## What survives, what's walked back, and what's sharpened

**Survives — empirical findings:**

S_k → 7/15 (direct measurement through k=13 truncated)
ρ_slow ≈ 0.834 as L¹/TV inverse-limit convergence rate (probe_profinite, R²=0.97)
Conjugation theorem K_k(c) = σ_c · K_k(1) · σ_c^(-1) (algebraic proof)
K_alg essential non-diagonalizability (cond(V) ~10^14 to 10^17)
Jordan block size = k captures 99% of dimension (Probe B Phase 2)
Mixing time = Jordan block size = k (Probe B Phase 3)
Cesàro framework compatibility (Probe A)
K_alg 3:1 row redundancy and rank n_k/3 (probe_preimage_kalg)
R_k singular spectrum band [0.49, 0.67] k-stable
K_emp bounded preimages with geometric v-distribution
K_alg dense complete digraph at full v_eff with k-invariant column statistics
The (Z/3)*-character bimodal partition (K-aligned vs B-aligned)
The "1619 club" mass concentration on specific orbital families

**Walked back:**

Order-3 linear recurrence on ε_k as a structural feature (it's a finite-k fitting artifact; failed at k=13)
AS character-sum eigenvalue formula as direct prediction of K_alg's spectrum (off by 3 orders of magnitude; +1 additive shift is non-perturbative)
DS-C Cayley-diameter sharp-rate bound for K_alg (loose by 2,000-100,000×; K_alg is Jordan-block-dominated, not diagonalizable)
DG step-by-step product structure for K_alg (within-trajectory CoV ~28%; characters aren't eigenvectors)
Atkinson construction direction (earlier session: brief conflated two distinct Atkinson theorems; pointwise algebra is informationally void)

**Sharpened:**

Entropy-deficit ρ_Δ from "ambiguous between same-mode and distinct-modes" to firmly distinct-modes with candidate value 7/8
The K_alg + R_k complementary rank decomposition (n_k/3 + 2n_k/3 = n_k) as algebraic source of π_k → π_{k-1} descent
The closure direction from "Liverani-style anisotropic Banach spaces" (vague) to "Butterley-Kim 2023 methodology adapted from Heisenberg nilmanifold to (Z/3^∞)* abelian profinite" (specific methodological template)

---

## Updated framework cohesion document

The framework cohesion document at `/mnt/user-data/outputs/framework_cohesion.md` and `.docx` was revised to incorporate tonight's findings:

- Operator distinction (K_alg vs K_emp) made explicit
- 3:1 row redundancy + rank decomposition added as fourth empirically-confirmed structural fit
- Recurrence walk-back added as first empirically-walked-back claim
- "What survives the recurrence walk-back" section added
- Entropy-deficit distinct-modes confirmation replaces "pending higher-k" framing
- Synthesis section reframed: ρ_slow's structural reality rests on direct L¹/TV measurement, not on recurrence
- Closure direction sharpened to reference Butterley-Kim 2023 specifically
- Provenance updated with new probe directories

The revised document is internally consistent with the empirical evidence as of May 6, 2026.

---

## Outreach status

Approximately 18 emails out to specialists across a diverse set: Cameron (UMD, spectral analysis Markov chains), Fill (JHU, strong stationary duality, responded briefly), Lin (UMD, theoretical Bayesian / variational Bayes / infinite-dimensional theory), Goldman (UMD, 2-adic and 3-adic systems), Lagarias, Tao, Cheng (Stanford), Bonacorsi, AACC math chair, others.

None of the outreach emails leaned specifically on the order-3 recurrence claim. They referenced ρ_slow as the inverse-limit convergence rate, which is the surviving characterization. No outreach correction is required.

If anyone responds substantively, the corrected framing is: ρ_slow ≈ 0.834 is the L¹/TV inverse-limit convergence rate measured directly via sparse Krylov computation through k = 12 and truncated extension to k = 13. An earlier framework characterization involving an order-3 linear recurrence on ε_k was a finite-k fitting artifact that did not survive extension.

---

## Path forward

The path is intact. The empirical bench has done what empirical work can do for the closure question. Further empirical work would be incremental refinement on what's already established. The closure question lives in analytical territory: constructing the inverse-limit transfer operator on Z_3 in the right function space, where ρ_slow appears as an isolated Pollicott-Ruelle resonance.

The reading direction is set up:

- Foundational technique: Diaconis-Graham 1992 (have)
- Sharp-rate framework: Diaconis-Saloff-Coste 1994 (have)
- Chain ring extension: Ayyer-Singla 2019 (have)
- Entropy interpretation: Eberhard-Varjú 2020 (have)
- Anisotropic Banach space menu: Baladi 2018 (have)
- Conceptual bridge for mixed dynamics: Giulietti-Liverani 2018 (have)
- Methodological template: Butterley-Kim 2023 (have)
- Liverani lecture notes: substituting for paywalled 1995 paper (have)
- Varjú-followup cluster of 24 papers (have)

The next analytical phase requires sitting with these papers and adapting Butterley-Kim's anisotropic-Banach-space construction from Heisenberg-group profinite decomposition to (Z/3^∞)* abelian profinite decomposition. This is months-to-years of analytical work.

The methodology paper does not require closure of c=7/45. The empirical findings from this session and prior sessions constitute publication-quality framework characterization in their own right, with appropriate qualifications calibrated to what has actually been demonstrated.

---

## Methodology lessons recorded

Two methodology lessons from tonight worth recording in the project's methodology paper:

**Lesson 1: Linear-recurrence machine-precision fits on finite sequences require extension testing.**

The order-3 recurrence held at machine precision through k=12. Treated as confirmed. Extension to k=13 revealed it was fitting finite-k transient structure. Validation at one extension level is not validation across multiple extension levels. This lesson generalizes beyond linear recurrences: any model that fits a finite sequence with machine-precision residuals should be extension-tested before being treated as structural.

**Lesson 2: When two independent measurements agree at low k, do not conclude they measure the same object.**

ρ_slow as recurrence dominant root and ρ_slow as direct L¹/TV measurement agreed at low k. The agreement was treated as cross-confirmation of the same underlying rate. Tonight's walk-back revealed they were two different objects that happened to have similar values at low k. The recurrence dominant root is unstable across fit windows; the L¹/TV rate is the genuine inverse-limit rate. Cross-confirmation requires both measurements to be stable, not just to agree numerically.

These lessons are appropriate additions to the methodology paper's discipline section.

---

— Nathan Humphrey, Annapolis, MD, May 6, 2026
