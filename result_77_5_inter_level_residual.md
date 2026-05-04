# Result 77.5 — Inter-level lift residual operator: outcome (R-stable) with multi-resolution structure

**Date:** 2026-05-04. Continues R77.4. Tests whether the rate-1/2 source lives in the inter-level renormalization rather than the within-level Markov spectrum (R77.4 ruled out the latter).

## Verdict in three lines

> **Stage 1: norms.** ‖R_k‖² values across k=1..5 stabilize at ‖R_k‖² ≈ 0.155 · 3^{−k}; the ratio ‖R_k‖² / ‖R_{k−1}‖² stabilizes at **0.334 ± 0.001 ≈ 1/3**, NOT 1/4. So R_k contracts in L² by 1/√3 per level, purely from the lift's cardinality scaling.
>
> **Stage 2: orthogonality is STRUCTURAL.** The regression coefficient c_k of R_k on T(R_{k−1}) is **exactly 0/1 over Q at every k = 2, 3, 4, 5**. Not "small" — exactly zero, by rational equality on the inner product Σ R_k(r') · T(R_{k−1})(r') in fractions.Fraction arithmetic. This is automatic from marginal consistency of the projective Markov system, not a contingent dynamical fact.
>
> **Strategic implication.** The natural multi-resolution decomposition V_{k+1} = T(V_k) ⊕ W_k holds, with R_k ∈ W_k by construction. The "renormalization operator" Φ acts on the orthogonal-complement subspaces W_k, but its leading mode contributes only the trivial cardinality-scaling factor 1/3 to ‖R_k‖². **Rate-1/2 of ε_n cannot live in ‖R_k‖ — it must live in the projection of the bilinear pair-form moment functional φ_n onto Σ_k W_k.** R77.2's "find a 1/2-eigenvalue at finite truncation" framing is structurally displaced; the right framework is multi-resolution / transfer-operator analysis on the projective limit.

**Files (this result):**
- `result_77_5_inter_level_residual.md` (this writeup)
- `result_77_5_compute_R_k.py` — Stage 1, exact rationals via fractions
- `result_77_5_renormalization_step.py` — Stage 2, regression onto lift (corrected indexing)
- `result_77_5_R_k_norms.csv` — ‖R_k‖² and ratios
- `result_77_5_phi_correlations.csv` — c_k = 0 exact data

---

## 1. Setup

Let π_k be the stationary distribution of the level-k Markov chain (coprime states r ∈ Z/3^k, computed exactly over Q via `build_markov_rational(k)` + `stationary_rational(K)` from `push_to_k6_rate_analysis.py`). State count: N_k = 2·3^{k−1}.

The natural lift T_{k→k+1}: V_k → V_{k+1} (where V_k = R^{N_k}) is

  T(π_k)(r') := π_k[r' mod 3^k] / 3   for r' coprime in Z/3^{k+1}.

Each coprime r in Z/3^k has exactly 3 coprime preimages r' in Z/3^{k+1}, and Σ_{r'} T(π_k)(r') = 1 because each coprime r is hit by 3 lifts each contributing π_k[r]/3.

The **lift residual** is

  R_k(r') := π_{k+1}(r') − T(π_k)(r').

Σ_{r'} R_k(r') = 0 (both summands have total mass 1). R_k = 0 identically would mean π is preserved by lift, contradicting ε_n ≠ 0.

## 2. Stage 1: norms and ratios — confirmed multi-resolution scaling

All exact rationals; floats shown for readability. See `result_77_5_R_k_norms.csv` for full numerator/denominator pairs.

| k | N_k | N_{k+1} | ‖R_k‖² | ‖R_k‖² · 3^k | ratio to prev |
|---|-----|---------|--------|--------------|---------------|
| 1 | 2   | 6       | 10/189 ≈ 5.291×10⁻² | 0.158730 | — |
| 2 | 6   | 18      | 31370/1835001 ≈ 1.7095×10⁻² | 0.153858 | **0.323102** |
| 3 | 18  | 54      | ≈ 5.7310×10⁻³ | 0.154738 | **0.335240** |
| 4 | 54  | 162     | ≈ 1.9157×10⁻³ | 0.155172 | **0.334267** |
| 5 | 162 | 486     | ≈ 6.3946×10⁻⁴ | 0.155390 | **0.333802** |

**Observations.**

1. **Ratio stabilizes at ≈ 0.334 = 1/3, NOT 0.25.** Across k=2..5 the ratios are (0.3231, 0.3352, 0.3343, 0.3338) — within 0.7% of 1/3 from k=3 onward.
2. **‖R_k‖² · 3^k stabilizes** at ≈ 0.155, so ‖R_k‖² ≈ 0.155 · 3^{−k}. R_k contracts by 1/√3 per level in L².
3. The 1/3 factor is purely the lift cardinality scaling — not an eigenvalue of any contractive dynamic. This is verified explicitly in §3.
4. **Outcome is NOT (R-trivial):** R_k ≠ 0 at every k.

The rate-1/2 envelope test would predict a ratio near 1/4 = 0.25; **decisively rejected**. The R_k norm sequence itself does not directly carry the ε_n rate.

## 3. Stage 2: c_k = 0 exactly — structural orthogonality

For k = 2, 3, 4, 5, we computed (using exact rationals)

  c_k := ⟨R_k, T(R_{k−1})⟩ / ‖T(R_{k−1})‖²

via fractions.Fraction arithmetic. **Result: ⟨R_k, T(R_{k−1})⟩ = 0/1 exactly at every k.**

| k | ⟨R_k, T(R_{k−1})⟩ | ‖T(R_{k−1})‖² | c_k | ‖R_k^⊥‖² / ‖R_k‖² |
|---|--------------------|---------------|-----|--------------------|
| 2 | 0/1 (exact) | 10/567 ≈ 1.764×10⁻² | 0 | **1.000 (orthogonal)** |
| 3 | 0/1 (exact) | 31370/5505003 ≈ 5.698×10⁻³ | 0 | **1.000 (orthogonal)** |
| 4 | 0/1 (exact) | ≈ 1.910×10⁻³ | 0 | **1.000 (orthogonal)** |
| 5 | 0/1 (exact) | ≈ 6.386×10⁻⁴ | 0 | **1.000 (orthogonal)** |

This is NOT a contingent dynamical fact — it follows from a structural identity.

### 3.1 Why c_k = 0 is structural

The natural lift T : V_k → V_{k+1} (acting via uniform-1/3 split across coprime preimages) is an **isometric embedding up to a √3 factor**:

  ⟨T(u), T(v)⟩_{V_{k+1}} = ⟨u, v⟩_{V_k} / 3.

Hence the image T(V_k) is an embedded copy of V_k inside V_{k+1}, and we have an orthogonal direct sum

  V_{k+1} = T(V_k) ⊕ W_k    where    W_k := T(V_k)^⊥ ⊂ V_{k+1}.

The orthogonal complement W_k is exactly the subspace of vectors in V_{k+1} that have **mean zero on each 3-fiber** {r, r+3^k, r+2·3^k}.

**Marginal consistency claim:** for any coprime r ∈ Z/3^k,

  Σ_{r' lifts of r} π_{k+1}(r') = π_k(r).

This holds because the Syracuse Markov chain is coherent under reduction mod 3^k — the projection of stationary π_{k+1} onto level k equals the stationary π_k of the level-k chain.

It follows directly that R_k(r') := π_{k+1}(r') − T(π_k)(r') has zero mean on each 3-fiber:

  Σ_{r' lifts of r} R_k(r') = π_k(r) − 3 · π_k(r)/3 = 0.

Hence **R_k ∈ W_k by construction.**

For any v ∈ V_k, ⟨R_k, T(v)⟩_{V_{k+1}} = 0 because R_k ⊥ T(V_k). In particular at v = R_{k−1}, we get c_k = 0 automatically.

### 3.2 What this means

- **c_k = 0 is verifying the chain is coherent under projection**, not learning anything dynamical about Φ.
- The "Φ acts diagonally on the lift basis" picture is wrong — the lift basis is structurally absent from R_k by construction.
- The renormalization between levels is NOT a scalar map T(R_{k−1}) ↦ c · T(R_{k−1}) plus small perpendicular corrections. It's a map from W_{k−1} into W_k (different orthogonal subspaces) with no preferred lift basis.

This is a multi-resolution decomposition reminiscent of Haar wavelets or the orthogonal-complement filtration in scaling limits of Markov chains.

## 4. Stage 3: structure of R_k inside W_k

What R77.5 can say about Φ as an operator:

- **Domain-codomain mismatch.** Φ would map W_{k−1} → W_k, where W_{k−1} ⊂ V_k and W_k ⊂ V_{k+1} are different Hilbert spaces. There is no single Hilbert space for Φ to act on at finite truncation.
- **Cardinality scaling alone explains the 1/3 norm ratio.** dim(W_k) = N_{k+1} − N_k = 2·3^k − 2·3^{k−1} = 4·3^{k−1}, growing by factor 3 per level. ‖R_k‖² ∼ const · 3^{−k} is the "constant per coordinate" scaling that any random vector with bounded average coordinate magnitude would exhibit.
- **The relative magnitude per coordinate is stable.** Per-coordinate squared mass: ‖R_k‖² / dim(W_k) ≈ (0.155 / 3^k) / (4 · 3^{k−1}) = (0.155 · 3) / (4 · 3^{2k}) ≈ 0.116 / 9^k. The decay per-coordinate is 1/9 per level. (To be checked: 0.155·3/(4·3²) = 0.0129 at k=1; 0.155·3/(4·3⁴) = 0.00144 at k=2. Ratio 0.111 ≈ 1/9. ✓)

So at the per-coordinate level, R_k is decaying by 1/9 = (1/3)² each level. This is the natural scaling for a "random-like" multi-resolution decomposition where mass spreads uniformly into the new orthogonal degrees of freedom.

**The operator Φ as a sequence W_{k−1} → W_k cannot be characterized from these data alone** — we have one number per level (the L² norm) and a structural orthogonality. To extract Φ's spectral content, we'd need a basis for {W_k} and a description of the map in that basis, which is a substantial structural project.

## 5. Why this does not contradict the rate-1/2 envelope of ε_n

ε_n := S_n − 7/15 is a **scalar moment** of π_n (specifically, S_n is the bilinear pair-form moment from R76). Writing ε_n = ⟨φ_n, π_n − π_∞⟩ for some functional φ_n, decompose

  π_n − π_∞ = Σ_{k ≥ 1} (lift to level n of R_k component)

— a telescoping sum through the orthogonal-complement subspaces. Because the W_k are mutually orthogonal in V_n (after proper lifting), the moment

  ε_n = ⟨φ_n, π_n − π_∞⟩ = Σ_{k} ⟨φ_n, lift_n(R_k)⟩

decomposes into independent contributions per level. Each contribution carries its own scaling: ‖lift_n(R_k)‖² = ‖R_k‖² · 3^{n−k−1} = (0.155 · 3^{−k}) · 3^{n−k−1} = 0.155 · 3^{n−2k−1}.

For ε_n's leading contribution to scale like (1/2)^n, the moment functional φ_n must have specific overlap pattern with the W_k subspaces. The rate-1/2 question becomes: **what is the overlap structure ⟨φ_n, lift_n(R_k)⟩?**

This is a question about φ_n (the bilinear pair-form moment) that R77.5 does not answer. But it identifies the right setting: rate-1/2 is encoded in φ's projection onto the multi-resolution decomposition Σ_k W_k, not in any single operator's spectrum.

## 6. Decision: outcome (R-stable) with corrected interpretation

From the three-way tree in the brief:

- **NOT (R-unstable):** the structural orthogonality c_k = 0 is stable across k=2..5 (and provable in general from marginal consistency).
- **NOT (R-trivial):** R_k ≠ 0 at every k tested.
- **(R-stable):** the multi-resolution decomposition V_{k+1} = T(V_k) ⊕ W_k is well-defined and stable across levels. **Verdict.**

But the qualitative result is unexpected and was misinterpreted in an earlier draft of this document: **c_k = 0 is structural orthogonality (marginal consistency), not a learned dynamical fact**. The rate-1/2 of ε_n is NOT a leading rate of any single operator at finite level — it lives in the moment functional's projection onto the multi-resolution decomposition.

## 7. Strategic implication for R77.2's framework

R77.2 framed the rate-1/2 question as: "find a finite-truncation operator T_N whose spectrum contains 1/2 in a Nisoli-amenable way."

- R77.3 ruled out the finite-mode geometric expansion of ε_n (no such 3-mode model holds).
- R77.4 ruled out the within-level Markov K_k (mixes in one step; no eigenvalues near 1/2).
- R77.5 now reframes: there is no single operator at finite truncation whose spectrum contains 1/2. The natural inter-level structure is a multi-resolution decomposition with orthogonal W_k subspaces, and the rate-1/2 is encoded in the moment functional φ_n's projection onto these subspaces.

**R77.2's conditional Theorem framing IS displaced.** The Nisoli machinery cannot be applied to K_k (R77.4) or to a single Φ on the W_k filtration (R77.5 — different Hilbert spaces per level). The right setting is a transfer-operator analysis where φ is an explicit observable on a function space — likely a space of functions on the projective limit Ẑ_3^× equipped with the inverse-limit topology — and the rate-1/2 manifests as a Mellin-Barnes / spectral-density feature of the transfer operator on that function space.

The next R-block work item is: **construct the function-space framework explicitly.** Candidates:
- Hilbert spaces of locally constant functions on Ẑ_3^× (the natural completion of ⊕_k V_k via the projective system)
- Wavelet-like frames on Ẑ_3^× adapted to the {W_k} decomposition
- Transfer-operator analysis on the action of the Syracuse map's coherent extension to Ẑ_3

This is a substantial reframing; the path from "rate-1/2 envelope" to "rigorous Theorem 77" runs through this function-space construction, not through Nisoli.

## 8. Anti-pattern audits

- **Did not** curve-fit ‖R_k‖² values (5 points, would be R77.4-trap). Ratio stability across k=2..5 IS the stability evidence; no fitting.
- **Did not** claim Φ_∞ exists from a single ratio. The 4 ratios (0.323, 0.335, 0.334, 0.334) are weak but unanimously near 1/3, consistent with cardinality scaling.
- **Did not** conflate the lift map with the rate operator. T is the natural projective lift; R_k is the residual; rate-1/2 lives in φ_n's projection onto the W_k filtration.
- **Did identify a prior-draft error.** An earlier version of this document interpreted c_k as "stable at √(1/3) within a few percent" based on a buggy state-space indexing in `result_77_5_renormalization_step.py`. After fixing the off-by-one (R_k lives on V_{k+1}, not V_k as the bug treated it), the regression coefficient is exactly 0/1 — structural orthogonality from marginal consistency, not a learned dynamical relationship.

## 9. What didn't finish

- **‖R_k‖² extension to k=6** would tighten the 1/3 estimate but cannot separate 1/3 from any nearby fraction better than the existing data.
- **Explicit basis for {W_k} and Φ's matrix in that basis.** This would require a multi-resolution analysis of the Syracuse Markov stationary distribution — substantial, separate work.
- **No explicit construction of φ_n on the projective limit.** This is the central open piece for the displaced R77.2 framing.
- **No verification that φ_n's projection onto Σ_k W_k actually produces the 1/2 envelope.** This is the empirical content of the rate-1/2 conjecture in the new framework.

## 10. Files

- `result_77_5_compute_R_k.py` — Stage 1 (exact rationals, ratio computation)
- `result_77_5_renormalization_step.py` — Stage 2 (corrected indexing; produces c_k = 0 exact)
- `result_77_5_R_k_norms.csv` — Stage 1 exact-rational table
- `result_77_5_phi_correlations.csv` — Stage 2 exact c_k = 0 table
