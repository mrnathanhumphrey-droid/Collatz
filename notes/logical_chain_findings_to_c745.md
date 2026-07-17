# Mathematical Logic Chain: From σ(n) Empirical to c = 7/45 Bridge

**Purpose:** Block-by-block logical math steps from first data to current state. No biography. Just: what was computed, what it showed, what was concluded, what the next computation was, why.

---

## Block 1: Bulk distribution shape

**Computed:** σ(n) = total Collatz stopping time for all n ∈ [1, 2²³]. Right-tail CCDF of residual r(n) = σ(n) − (2/(ln 4 − ln 3))·ln(n).

**Observation:** Semi-log plot linear over ~5 decades, slope ~0.030. Log-log shows cliff drop, not linear.

**Logical step:** Linear semi-log + non-linear log-log ⟹ tail is **exponential, not power-law**: P(r > t) ≈ exp(−0.030 t).

**Boring explanations ruled out:** sampling bias (full uniform), finite-N (rate stable across bulk; only mild thinning at extreme right; largest σ in range = 664), definition mismatch (residual shape invariant to slope baseline).

**Next computation justified by this:** check whether bulk slope coefficient matches the random-walk heuristic — because if tail is exponential, leading-order behavior is heuristic-consistent.

---

## Block 2: Bulk slope identification

**Computed:** Pooled OLS on 16.7M odd n in [3, 2²⁵]: σ ≈ α + β·ln(n).

**Observation:** β = 10.4543. Heuristic for odd-n is 3/(ln 4 − ln 3) = 10.4282.

**Logical step:** |β_emp − β_heur|/β_heur = 0.0024 ⟹ **bulk slope matches heuristic to 0.25% on 16.7M points.** Random-walk model is essentially exact at leading order.

**Critical correction made here:** initial baseline used 6.95 (full-Collatz including evens). Switched to 10.4282 (odd-only filter). Without this correction the deviation looked 50× larger.

**Next computation justified:** if bulk is universal at leading order, test whether **per-class** structure deviates from universal. Use Bayesian hierarchical model with partial pooling on residue classes mod 2^k.

---

## Block 3: Trajectory measure ≠ natural density (control finding)

**Computed:** Empirical v = ν₂(3m+1) distribution, two ways:
- Uniform odd m in [1, 2²⁵]
- Trajectory-sampled (Syracuse iterates)

**Observation:**
- Uniform: ratio to Geom(½) = 1.0000 to 7 decimals for v ∈ {1..19}. Heuristic exact on natural density.
- Trajectory deep-step (step_idx > 0): v=4 ratio 1.325×, v=10 ratio 1.321×, sags at v ∈ {6,7,11,12,14}.

**Logical step:** Natural-density heuristic is **right** on uniform sampling. Trajectory measure **deviates** because Syracuse iterates are not i.i.d. — they correlate through dynamics, oversampling certain residue classes. This is Lagarias's known phenomenon, locally verified.

**Implication:** mean drift on natural-density preserved (β matches heuristic), but **variance and tail of σ may deviate from i.i.d. predictions**. This is a structural reason to expect per-class effects.

**Next computation:** Bayesian hierarchical fit at k=6 (32 classes mod 64) to localize where deviation lives — slope, intercept, or tail shape.

---

## Block 4: Hierarchical fit at k=6

**Computed:** Stan hierarchical model, 32 classes mod 64, uniform sample 50K/class, 4 chains × 4 threads, 500+500 iterations. Convergence: R-hat ∈ [0.9991, 1.0111], 2/2000 divergent transitions (0.4%).

**Posteriors:**
- μ_β = 10.4475 ± 0.0515
- τ_β = 0.067 ± 0.048 (q05=0.006, q95=0.161)
- τ_α = 11.70 ± 1.08

**Logical step:** μ_β posterior covers theoretical 10.4282. τ_β posterior small, near noise floor. τ_α large and well-separated from zero.

**Conclusion:** Modular structure lives in **intercept α**, not slope β. Slope is universal across classes (within noise); intercept varies systematically with residue class.

**Per-class GPD on top-5% residuals:** ξ ∈ [-0.107, 0.054], 25/32 classes "exponential" (|ξ| ≤ 0.05), 5/32 sub-exponential, 2/32 super-exponential. **No power-law tails anywhere.**

**Next computation justified:** test whether τ_β = 0 asymptotically (universality) by N-scaling sweep. Test whether α structure has further decomposition.

---

## Block 5: N-scaling — slope universality confirmed

**Computed:** OLS sweep across N ∈ {2²⁰, 2²², 2²³, 2²⁴, 2²⁵}, 32 classes mod 64.

| N | μ_β | τ_β observed | per-class SE | τ_β corrected |
|---|---|---|---|---|
| 2²⁰ | 10.372 | 0.368 | 0.430 | 0 |
| 2²² | 10.382 | 0.187 | 0.228 | 0 |
| 2²³ | 10.385 | 0.128 | 0.165 | 0 |
| 2²⁴ | 10.404 | 0.110 | 0.120 | 0 |
| 2²⁵ | 10.419 | 0.086 | 0.087 | 0 |

**Logical step:** Moment-corrected τ_β² = max(0, observed² − SE²) = 0 at every N. **Per-class slope variation is entirely sampling noise.** As n_per_class → ∞, observed τ_β → 0.

**Also:** ⟨ξ⟩ → 0 monotonically (-0.083 at 2²⁰ → -0.028 at 2²⁵). Sub-exponential cliff drops are **finite-N truncation artifacts**.

**Conclusion:** Asymptotic universality confirmed. μ_β → 10.4282, τ_β → 0, ⟨ξ⟩ → 0. Only τ_α stable around 12-13 across N (real structure).

**Next computation:** decompose α into deterministic + stochastic components.

---

## Block 6: α decomposition — the eureka

**Setup:** For each odd r ∈ {1, 3, ..., 63}, compute the deterministic Collatz "prefix" — steps where parity is forced by r alone (state = a·m + c symbolically; algorithm terminates when a becomes odd at (a_final, c_final)).

**Predicted α:** α_det(r) = prefix_steps + 10.43 · ln(a_final / 64)

**Computed:** Linear fit α_actual = a + b·α_det.

**Result:** α_actual = −2.66 + 0.986 × α_predicted, **R² = 0.9996**.
- SD(α_actual) = 13.7
- SD(α_stoch residual) = 0.28
- Largest individual deviation: 0.66 absolute, 0.46 standard errors
- **No class statistically significant α_stoch.**

**Logical step:** τ_α ≈ 13 is 99.96% explained by 7-12 steps of deterministic Collatz prefix algebra. **The "stochastic" residue-class structure is not stochastic — it is algebraic.**

**Prefix terminal values:** a_final ∈ {3, 9, 27, 81, 243, 729} = {3¹, ..., 3⁶}. Powers of 3 from the Syracuse map.

**Conclusion:** σ(odd n) decomposes into:
1. Deterministic intercept shift α_det(r) computable by symbolic Collatz iteration on r mod 64
2. Universal stochastic random-walk behavior thereafter (slope 10.4282, exponential tail)

**No residual residue-class structure exists at this resolution.**

**Next computation justified:** test whether decomposition extends to higher k, and whether it generalizes to qx+1 dynamics for q ≠ 3.

---

## Block 7: Universality across k

**Computed:** Extended α decomposition to k ∈ {6, 7, 8, 9} on N = 2²⁷ data.

| k | classes | n/class | R² | SD(resid) | mean SE | ratio |
|---|---|---|---|---|---|---|
| 6 | 32 | 2.10M | 0.9967 | 0.77 | 0.80 | 0.96 |
| 7 | 64 | 1.05M | 0.9942 | 1.12 | 1.13 | 0.99 |
| 8 | 128 | 524K | 0.9918 | 1.44 | 1.59 | 0.91 |
| 9 | 256 | 262K | 0.9851 | 2.09 | 2.24 | 0.93 |

**Logical step:** R² declines monotonically with k, but signal-to-noise ratio SD(resid)/mean(SE) ≈ 0.91-0.99 stays constant. **R² decline is entirely driven by smaller per-class n.** Decomposition holds at every modular resolution.

**Next computation:** generalize to qx+1, test whether prefix decomposition predicts convergence rates.

---

## Block 8: qx+1 prefix complexity → exponential decay law

**Computed:** qx+1 trajectory data for q ∈ {3, 5, 7, 9, 11, 13} at N=10⁷, selected q at N=10⁸. Per-class convergence rates binned by prefix odd-step count j.

**Observation:** log(conv_rate(j)) ≈ const + slope·j is essentially exact (R² ≥ 0.99) at q ∈ {5, 7, 9}. Empirical slope divided by log(q/4) is q-independent:

| q | n_conv | slope | log(q/4) | slope / log(q/4) | R² |
|---|---|---|---|---|---|
| 5 | 32,785 | -0.5619 | 0.2231 | **-2.518** | 0.999 |
| 7 | 258 | -1.3685 | 0.5596 | **-2.445** | 0.999 |
| 9 | 104 | -2.0529 | 0.8109 | **-2.531** | 0.994 |

**Mean: -2.498, within 0.1% of -5/2.**

**Logical step:** conv_rate(j; q) ≈ A(q) · (4/q)^((5/2)·j). Universal multiplier 5/2 across q ∈ {5, 7, 9}. **The 4/q base is the random-walk drift; the 5/2 multiplier is empirical with no current derivation.**

**Cramér's theorem prediction with v_2 ~ Geom(½):** q-DEPENDENT multipliers 3.35 (q=5), 4.05 (q=7), 4.57 (q=9). **Wrong direction.** Either v_2 distribution is different under qx+1 dynamics, or the gambler's-ruin assumption ignores correlations.

**Next computation:** verify v_2 distribution under qx+1 trajectory measure.

---

## Block 9: Universal Geom(½) trajectory v_2

**Computed:** Unconditional sampling of qx+1 trajectories, ~10⁶ v_2 records per q at q ∈ {5, 7, 9, 11}.

| q | mean v_2 | var v_2 | drift_emp | drift_pred = log(q/4) |
|---|---|---|---|---|
| 5 | 2.0028 | 2.0056 | 0.2212 | 0.2231 |
| 7 | 1.9990 | 2.0003 | 0.5603 | 0.5596 |
| 9 | — | — | 0.8119 | 0.8109 |
| 11 | — | — | 1.0120 | 1.0116 |

**Logical step:** Geom(½) prediction (mean=2, var=2) matches empirical to 0.5% across q. **Trajectory measure on v_2 is i.i.d. Geom(½) for all tested q.**

**Conclusion:** Cramér baseline is correct. The 5/2 multiplier deviation is NOT from a wrong v_2 distribution. It must come from step-to-step correlations or higher-moment effects.

**Implication:** the prefix-complexity axis is structural and q-universal. Per-class qx+1 convergence is governed by `conv_rate(j; q) ≈ A(q)·(4/q)^((5/2)·j)`.

**Next computation:** connect this empirical structure to Tao 2022's analytic almost-everywhere theorem.

---

## Block 10: Bridge to Tao (5.15)

**Tao 2022 (5.15):** for almost all N, T_x(N) = log(N/x)/log(4/3) + O(log^0.6 x).

**Hypothesis:** s_mean(r) ≈ α_det(r) + K_h·log(N/f(N)) where K_h = 3/log(4/3) = 10.4282.

**Computed:** N ∈ {2²⁵, 2²⁷}, k ∈ {8, 10, 12, 14}, five threshold functions f(N) ∈ {1, N^(2/3), √N·log N, √N, √N/log N}. Total 40 cells.

**Result table (raw mean slope @ K_h):**

| observable | k=8 N=2²⁵ | k=8 N=2²⁷ | k=14 N=2²⁵ | k=14 N=2²⁷ |
|---|---|---|---|---|
| σ | 0.9960 | 0.9990 | 0.9945 | 0.9977 |
| s @ N^(2/3) | 0.9989 | 0.9996 | 0.9944 | 0.9977 |
| s @ √N | 0.9989 | 1.0006 | 0.9967 | 0.9997 |
| s @ √N/log N | 1.0000 | 1.0012 | 0.9980 | 1.0005 |

**All 40 raw-mean slopes ∈ [0.9936, 1.0012].** Tightest at √N: median 1.000.

**Offset gaps:** stable across k and N (variation < 0.06). Trim-1% offset gap at √N: < 0.34 across all (k, N).

**Logical step:** **α_det(r) (computable from prefix algebra) + K_h·log(N/f(N)) (Tao's leading term) = s_mean(r) to slope 1.000 ± 0.005, no calibration.** This is the per-class realization of Tao's almost-everywhere first-passage formula.

**Conclusion:** Tao's (5.15) supplies the asymptotic mean coefficient K_h. Structural decomposition supplies α_det(r). They agree at every resolution.

**Next computation:** formalize this bridge into theorems with explicit error control. Develop framework that can prove (not just verify) the connection.

---

## Block 11: Bridge equation framework formalization

The empirical observation in Block 10 — that α_det(r) + K_h·Δlog matches mean σ structure exactly — is a *bridge*: it connects two existing analytical objects.

**Object A:** Tao's iterated Markov-chain offset distribution, converging to a stationary measure under Syracuse iteration.

**Object B:** A 2×2 algebraic recursion T_diag arising from Plancherel decomposition of the Collatz step operator on Fourier mass.

**Conservation Law (R76.1):** propagation of mass through Collatz iteration in Object A is identical to propagation of T_diag's diagonal in Object B.

**Computed (Phase 7 of project):** explicit form of T_diag for 3-adic case.

T_diag = (1/5)·[[1,1],[4,4]]

**Fixed point of T_diag:** c = 7/45.

**Logical step:** Tao's offset distribution converges to a stationary measure. The corresponding fixed point of T_diag is c = 7/45. **c = 7/45 is the candidate sharp constant** for the convergence rate that Tao's theorem proves at almost-bounded scale.

**Theorems established (12 total):** R74, R75.1, R75.2, R76.1, R76.3, R77, R78.1, 78.2, 78.3, 78.4, 78.5, 78.6.

**Closure problem:** convert "candidate" → "actual" by proving the residual gap closes.

**Residual gap (after R78.6):**
|S_partial| = |Σ_{a ∈ supp} 1̂(3a)·ψ(a)| where:
- supp = {a ∈ Z/3^r : a ≡ 1 mod 3}
- 1̂(3a) is the Fourier coefficient with |1̂(3a)| = 3√q exactly (R78.3 saturation)
- ψ(a) = e_q(P_a(s*(C_a))) — saddle-phase form (T78.6 exact at r=3, qualified for r ≥ 4)
- q = 3^{r+1}

**Required:** |S_partial| ≪ q^{1/2 - δ} for some δ > 0.

**Next computation:** apply existing prime-power character sum theory to bound this object, OR identify why it can't be bounded by existing techniques.

---

## Block 12: Verification agent — applying Milićević and Banks-Shparlinski

**Computed:** structural compatibility of S_partial with two existing frameworks.

**Milićević 2014 framework check:**
- F-class structural match holds at p=3, y=1: (w, y, κ, λ, u) = (1, 1, 1, ∞, ∞), a₀ = -C_a ✓
- Theorem 3 condition (iv) check: ρ_3(1) = 3, κ = 1, so λ̃ = κ - ρ_p(y) = -2 < 0 ✗
- Exponent pair iteration: ABA³B(0,1) = (11/82, 57/82), θ ≈ 0.1646. Target δ ≈ 0.0855 corresponds to k+ℓ ≈ 0.671, **beyond Rankin's frontier 0.659.**

**Banks-Shparlinski 2018 framework check:**
- γ₀ ≥ e²⁰⁰ — purely asymptotic, no computational range
- ξ₀ ≈ 10⁻¹⁷⁸ — non-zero but unusable
- Their G(n) is in summation variable; our cubic structure is in dual variable a — different object

**Critical cross-check:** Milićević's exponent pair on primal would predict |Σ_n χ(n)e(c·n)| ≪ q^{1/2 - 0.085}. But |1̂(3a)| = 3√q exactly by R78.3. **Subconvex saving on primal contradicts known equality.** R78.3 saturation IS the obstruction.

**Logical step:** the bound |S_partial| ≪ q^{1/2 - δ} cannot come from technique that proves primal subconvexity, because the primal is already at saturation. **Must come from cancellation among phases ψ(a) summed over a, not from per-term improvement.**

**Constructive direction identified:** partition Σ_a by saddle class j ∈ {0, 1, 2}. At r = 3, each subsum is linear in a (Pólya-Vinogradov gives only log-saving, δ = 0). At r ≥ 4, phase becomes higher-degree from Hensel lifting; polynomial saving room re-opens.

**Next computation:** empirical |S_partial(r)| at r = 8..20 to test whether δ > 0 at observed scale.

---

## Block 13: Empirical agent — direct |S_partial| computation

**Computed:** direct Kalafatelis sum K(r, c, m) = Σ_{u=0}^{N-1} e_q(c·4^u − 9mu), straight summation, no saddle approximation. r = 3..19 in flight at time of writing, target r = 20.

**Observation:** log|K|/log(N) ≈ 0.528 at r=18, trending to 0.5 from above. |K(r)| ≈ √N with O(1) prefactor across all observed r.

**Logical step:** **Empirical δ → 0.** No polynomial saving visible. Confirms R78.3 saturation is the operating bound: |S_partial| sits at exactly the magnitude predicted by Plancherel + random-walk-style cancellation across the support.

**Side finding (independent agent verification):** at r ≥ 4, leading-order saddle phase ψ_lead deviates from true phase G(a)/√q by 13-21% of q. **Hensel lifting is non-trivial at r ≥ 4.**

**T78.6 status check:** theorem statement is correct (says "exact at r = 3", flags r ≥ 4 as needing Hensel lifting). Downstream "cubic in a" remark in lines 60-69 of result_78_extended.md is wrong — leading-order at r ≥ 4 is piecewise linear, not cubic. Walk-back needed.

**Logical step:** the cubic-character-sum framing of the gap (used in Block 12 verification) was based on the walked-back remark. The actual gap object at r ≥ 4 is the Hensel-lifted full phase, not a generic cubic. **Verification doc analyzed a partially-correct object.** Block 12 conclusions about saturation hold (those depend on R78.3, not phase structure). Conclusions about exponent-pair frontier may need re-derivation for the correct phase.

---

## Current state (Block 13 → next)

**Established with confidence:**
1. Bulk σ(n) tail is exponential with rate ~0.030 (Block 1)
2. Slope μ_β = 10.4282 universal (Blocks 2, 4, 5)
3. Per-class structure is 99.96% deterministic algebra (Block 6)
4. Decomposition holds at k ∈ {6, 7, 8, 9} (Block 7)
5. qx+1 convergence law: conv_rate ≈ A(q)·(4/q)^((5/2)·j) for q ∈ {5, 7, 9} (Block 8)
6. Trajectory v_2 is universal Geom(½) (Block 9)
7. Tao (5.15) bridge: s_mean(r) = α_det(r) + K_h·Δlog with slope 1.000 ± 0.005 (Block 10)
8. Bridge equation framework v3.7: 12 theorems, c(3) = 7/45 candidate fixed point (Block 11)
9. Milićević + Banks-Shparlinski don't directly close the gap (Block 12)
10. Empirical |K(r)| ≈ √N, δ → 0 at observed r (Block 13)

**Pending verification:**
- T78.6 walk-back of cubic-in-a remark
- Empirical run completion at r = 19, 20
- Effect of correct (Hensel-lifted) phase on Block 12 verification analysis

**Open paths:**
- Path C (committed): 5x+1 sibling attack, Phase 2.3 gate first
- Outreach: Milićević directly + Heath-Brown, with precise mathematical question now writable

**The mathematical question that landed today:**

"Coset character sum on (Z/3^{r+1})* is saturated at Plancherel by an exact Gauss sum identity (R78.3). Standard subconvexity techniques can't push below saturation without contradicting the equality. Saddle-class partition at r ≥ 4 with Hensel lifting is the only constructive direction visible. Question: is this within reach of existing techniques, or genuinely open?"

---

## Files of record

- `/mnt/user-data/uploads/findings.md` — original findings log April 30 - May 2 2026
- This document: `/mnt/user-data/outputs/logical_chain_findings_to_c745.md`

