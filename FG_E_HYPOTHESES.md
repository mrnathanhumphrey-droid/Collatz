# FG Candidate E — Bourgain discretized sum-product Fourier decay on R^n (Li 2018)

**PDF:** Discretized_SumProduct_Fourier_Decay_Rn.pdf (Li 2018, generalizing Bourgain 2010).
**Extracted text:** `C:/tmp/fg/bourgain_sumproduct_rn.txt`.

---

## Li 2018 Theorem 1.1 (VERBATIM, p. 1)

> "Given κ_0 > 0, there exist ǫ, ǫ_1 > 0 and k ∈ N such that the following holds for δ > 0 small enough. Let μ be a probability measure on [1/2, 1]^n ⊂ R^n which satisfies (δ, κ_0, ǫ) projective non concentration assumption, that is,
>     ∀ρ ≥ δ, sup_{a ∈ R, v ∈ S^{n-1}} (π_v)_* μ(B_R(a, ρ)) = sup_{a,v} μ{x | ⟨v, x⟩ ∈ B_R(a, ρ)} ≤ δ^{-ǫ} ρ^{κ_0}. (1.1)
> Then for all ξ ∈ R^n with ||ξ|| ∈ [δ^{-1/2}, δ^{-1}],
>     |μ̂_k(ξ)| = | ∫ exp(2iπ ⟨ξ, x_1 … x_k⟩) dμ(x_1) … dμ(x_k) | ≤ δ^{ǫ_1}. (1.2)"

### Hypotheses (typed):

- h_E.1.1.ambient: μ probability measure on [1/2, 1]^n ⊂ R^n. The space is **R^n with coordinate-multiplication** product structure. [TYPE (i)]
- h_E.1.1.nonconc: (δ, κ_0, ǫ) **projective non-concentration**: μ has no concentration on hyperplanes (one-dim linear subspaces) at any scale ρ ≥ δ. [TYPE (iii) — regularity of μ]
- h_E.1.1.frequency: ||ξ|| ∈ [δ^{-1/2}, δ^{-1}] (intermediate frequency window). [TYPE (iv) — output regime]

### Conclusion C_E.1.1:

- |μ̂_k(ξ)| ≤ δ^{ǫ_1}: polynomial Fourier decay of the k-fold multiplicative convolution. **Quantitative**, polynomial in δ.

The conclusion is essentially polynomial-in-A on the Fourier side (with the size parameter being δ, and the frequency being polynomial in 1/δ).

---

## Phase 1 — hypothesis × input matrix

| Hypothesis | (1) Tao | (2) C1 | (3) R75/76/77 | (4) eps_k |
|---|---|---|---|---|
| h_E.ambient: μ on [1/2,1]^n ⊂ R^n with coord. mult. | Syracuse μ_n lives on **(Z/3^n)***, a finite cyclic group, NOT on R^n. The chain step **does** use multiplication (x → x · 2^{-v} mod 3^n), but the modulus is a discrete prime power, not [1/2,1]^n ⊂ R^n. FAILED. | FAILED | FAILED | FAILED |
| h_E.nonconc: projective non-concentration of μ | The Syracuse stationary π_n on (Z/3^n)* has class-asymptotic-fractions (1/3, 2/3) for classes c=1, c=2 (R64.B); within each class it is approximately uniform (Plancherel μ̂_n decay → 0 fast for high ξ). Non-concentration on hyperplanes is vacuous on a finite cyclic group (no "hyperplane" structure). **NEEDS_PROOF / N/A**. | N/A | NEEDS_PROOF | NEEDS_PROOF |
| h_E.frequency window ||ξ|| ∈ [δ^{-1/2}, δ^{-1}] | Syracuse's frequency variable ξ ∈ Z/3^n with 3 ∤ ξ; size of ξ is bounded by 3^n. The "δ" analog would be δ ≈ 3^{-2n} (to put 3^n into [δ^{-1/2}, δ^{-1}]); but the conclusion δ^{ǫ_1} = 3^{-2ǫ_1 n} translates to exponential decay in n, which IS polynomial-in-A in the natural Syracuse parameter. SATISFIED if dictionary translation is admissible. | SATISFIED | SATISFIED | SATISFIED |
| h_E.k iterations of μ | Syracuse has n iterations of the chain step; matches "k-fold multiplicative convolution" but k is **fixed** in Li's theorem (k depends on κ_0) while n is growing. NEEDS_PROOF (whether the bound is uniform in n is the question). | NEEDS_PROOF | NEEDS_PROOF | NEEDS_PROOF |

**Phase 1 disposition: NEEDS_PROOF dominant.** The ambient mismatch is the categorically harder issue: R^n with coord. mult. vs (Z/3^n)* with the Tao chain.

---

## Phase 2 — conclusion shape

If hypotheses were satisfied (modulo translation): the bound |μ̂_k(ξ)| ≤ δ^{ǫ_1} translates to |μ̂_n(ξ)| ≤ 3^{-2ǫ_1 n} = (3^n)^{-2ǫ_1}.

With ξ ranging over Z/3^n in the natural window, this is polynomial-in-3^n decay, i.e., polynomial-in-A (with A = 2ǫ_1 a fixed constant depending on chain regularity, NOT free over A).

**Conclusion shape: PARTIAL match.** The shape is polynomial decay but with a fixed exponent ǫ_1, not "any A > 0". Sharper than Tao's qualitative S_χ(n) ≪_A n^{-A} only in the case where ǫ_1 can be made large; the standard Bourgain sum-product gives ǫ_1 small (typically ~10^{-3}). This delivers polynomial decay, *not* polynomial-in-A.

Note: the relationship "polynomial decay with fixed exponent" ↔ "any-A decay" requires bootstrapping (iterate the bound on convolutions) — this is exactly how Tao's recursion converts |f(x,3)|-bound at single scale to S_χ(n) ≪_A n^{-A}. So **Li 2018 / Bourgain delivers the single-scale ingredient**, and Tao's iteration scheme converts that to any-A. They are complementary, not equivalent.

---

## Phase 3 — profinite extension

Can Li/Bourgain extend to (Z/3^n)* with chain ξ → ξ · 2^{-v} mod 3^n?

The key ingredient of Li/Bourgain is the **discretized sum-product theorem on R^n** (Li Thm 1.4, p. 2): if A, X satisfy non-concentration and dimension assumptions, then N_δ(X+X) + sup_{a ∈ A} N_δ(X + aX) ≥ δ^{-ǫ} N_δ(X).

The analog on (Z/3^n)*: **Bourgain-Glibichuk-Konyagin sum-product on Z/p^n.** This exists and is studied. The discretized version on multiplicative groups of prime-power moduli was a goal of the BGK / shparlinski-school work scoped in cluster C1 (closed-negative on Syracuse for *different* reasons — Cochrane 1.2 degree-blowup; HB phase incompatibility).

But: the C1 disposition closed Bourgain-school *complete sum* machinery on Syracuse on **category-of-object** grounds — Syracuse is on tuple-space, not Z/p^n. The same applies here: Li Thm 1.1 requires μ supported on R^n with sup-norm regularity, whereas Syracuse step is **on tuple-space N+1 via 2^{-v}**, mapping into (Z/3^n)* via the chain. The non-concentration assumption (1.1) on the *step distribution* (2-adic Geom(2)) — interpreted as concentration on hyperplanes of R^n — is **vacuously satisfied** (Geom(2) is supported on the orbit of 2 ∈ Z_3^*, a 1-dim object; the "hyperplane" non-concentration is vacuous).

But the *iterated multiplicative convolution μ_k* on Z_3^* via 2^{-v_1} · 2^{-v_2} · … · 2^{-v_k} = 2^{-(v_1+…+v_k)} — this is **still a power of 2**. The Syracuse chain is NOT μ_k of Geom(2), because Tao's recursion is not a simple multiplicative convolution: it's a *3-step Markov chain* with multiplication by 2^{-v} combined with **modulus reduction from 3^{n+1} to 3^n** (the 3-adic descent). The "multiplicative convolution k times" interpretation of Li/Bourgain would give a measure on ⟨2⟩ ⊂ Z_3^* concentrated at 2^{-(v_1+…+v_k)} — a 1-dim object, not the full stationary.

**Phase 3 disposition: STRUCTURALLY_BLOCKED at the iteration mapping.** The k-fold multiplicative convolution doesn't reconstruct the Syracuse stationary measure π_∞ on Z_3^*; π_∞ has full support on Z_3^* with class-fractions (1/3, 2/3), not just on ⟨2⟩.

---

## Disposition E: **NO_FIT** (categorical at iteration-mapping; partial conclusion-shape match).

- Phase 1: ambient mismatch (R^n vs Z/3^n) + iteration-mapping mismatch.
- Phase 2: conclusion shape is the right *kind* (polynomial decay) but with fixed exponent, requiring Tao-iteration bootstrap to reach polynomial-in-A.
- Phase 3: the multiplicative-convolution interpretation doesn't reconstruct the Syracuse stationary measure — k-fold μ on Z_3^* only fills ⟨2⟩ ⊂ Z_3^*, not the full stationary support.

**Caveat**: Li/Bourgain Theorem 1.1 IS the single-step ingredient that Tao Section 7 (Lemma 7.2 ⊕ pair-grouping) reproduces for Syracuse. Tao's *proof* of S_χ(n) ≪_A n^{-A} **already uses** a discretized-sum-product-style argument at the single-scale level (the |f(x,3)| ≤ cos(πθ) bound is essentially a single-step Fourier coefficient bound). So Li 2018 doesn't add new content beyond what Tao already implements; the rigorous c=7/45 rate (1/30 coefficient and (1+3^n)-fine-frequency identification) requires sharper machinery on the **iteration / pair-grouping** side, not the single-step side.
