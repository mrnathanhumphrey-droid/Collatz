# ADELIC_H — Bourgain-Gamburd-Sarnak generalized Selberg 3/16 + affine sieve

**Source:** C:/tmp/adelic/Bourgain_Gamburd_Sarnak_Selberg_316_Affine_Sieve.txt. Main results Theorems 1.1, 1.2, 1.3.

## Verbatim Theorem 1.1

> "Theorem 1.1. Let Λ be a finitely generated subgroup of SL(2, ℤ) with δ(Λ) > 1/2. For q ≥ 1 let Λ(q) be the 'congruence' subgroup {x ∈ Λ : x ≡ I mod q}. There is ε = ε(Λ) > 0 such that
>    λ_1(Λ(q)) ≥ λ_0(Λ(q)) + ε,
> for all square-free q ≥ 1 (note that λ_0(Λ(q)) = λ_0(Λ))."

## Verbatim Theorem 1.2

> "Theorem 1.2. Let Λ = ⟨S⟩ be a finitely generated subgroup of SL(2, ℝ) with δ(Λ) > 1/2. Let {N_i} be a family of finite index normal subgroups of Λ. Then the following are equivalent
> (1) The Cayley graphs G(Λ/N_i, S) form a family of expanders.
> (2) There is ε = ε(Λ) > 0 such that λ_1(Λ/N_i) ≥ λ_0(Λ/N_i) + ε."

## Verbatim Theorem 1.3 (lattice point count)

> "Theorem 1.3. Let Λ be a finitely generated subgroup of SL(2, ℤ) with δ(Λ) > 1/2. Assume that q is square free and (q, q_0) = 1, where q_0 is provided by the strong approximation theorem [19]. There is ε_1 > 0 depending on Λ such that for any g ∈ SL_2(q) we have
>    |{γ ∈ Λ | ∥γ∥ ≤ T and γ ≡ g mod q}| = c_Λ T^{2δ} / |SL_2(q)| + O(q^3 T^{2δ − ε_1})."

## Hypotheses isolated

- **h1 (GROUP):** Λ a finitely generated subgroup of SL(2, ℤ) (or SL(2, ℝ)).
- **h2 (DIMENSION):** δ(Λ) > 1/2, where δ is the Hausdorff dimension of the limit set of Λ acting on ℍ² = hyperbolic plane.
- **h3 (ORBIT):** ⟨S⟩ generates infinite-volume hyperbolic surface X = Λ\ℍ²; Cayley graphs of congruence quotients exist.
- **h4 (CONGRUENCE STRUCTURE):** Λ(q) := {x ∈ Λ : x ≡ I mod q} is well-defined; q square-free.
- **CONCLUSION:** λ_1(Λ(q)) ≥ λ_0(Λ(q)) + ε uniformly in q; this is a *spectral gap* on infinite-volume congruence covers, generalizing Selberg's 3/16.

## Hypothesis × input check

| Hyp | Syracuse setting |
|---|---|
| h1 (Λ ⊂ SL(2, ℤ)) | The Tao recursion family {M_v = [[3, 1], [0, 2^v]] : v ≥ 1} is a subset of SL(2, ℚ) (with appropriate normalization) but NOT of SL(2, ℤ): det(M_v) = 6·2^{v-1}, not 1. After projectivizing to PGL(2, ℚ), the elements live in PGL_2(ℚ), but BGS works in SL(2, ℤ). Categorical: the integer condition fails. |
| h2 (δ(Λ) > 1/2) | The limit set of {M_v} in ∂T_3 = ℙ¹(ℚ_3) is a single point (∞) — the boundary fixed point of the pencil per BT_CANDIDATE_CONSTRUCTIONS. Hausdorff dimension of a single point is 0, not > 1/2. Also wrong ambient space (BGS works in ℍ², Tao acts on T_3). |
| h3 (infinite-vol hyperbolic surface) | Tao doesn't produce a hyperbolic surface. Different categorical setting. |
| h4 (congruence Λ(q)) | Not applicable — Λ isn't in SL(2, ℤ). |

## Disposition for H

**NO_FIT (categorical).**

BGS's framework requires:
- A finitely generated subgroup of SL(2, ℤ) acting on ℍ² (archimedean hyperbolic plane);
- Hausdorff dimension > 1/2 of the limit set;
- Congruence covers with spectral gap.

Tao's family is in PGL(2, ℚ_3) (3-adic), not SL(2, ℤ); its limit set in ℙ¹(ℚ_3) is a single point; there's no hyperbolic surface. Different categories.

**Adelic factorization tag:** **GLOBAL_BUT_PLACE_BLIND** — BGS uses archimedean ordering ∥γ∥ ≤ T but produces a spectral gap statement on the *archimedean* Laplacian Δ on ℍ². No place-by-place factorization in the Tate sense. The "archimedean" here is hyperbolic-symmetric-space archimedean, not Syracuse's "trajectory is bounded archimedean-ly" archimedean. Different uses of "archimedean".

## Note on speculative re-interpretation

If one COULD find an alternative algebraic completion of Tao's group structure that lands in SL(2, ℤ) (not PGL(2, ℚ_3)), then BGS might apply. The polynomial-in-A bound BGS gives (via Theorem 1.3, with O(q^3 T^{2δ − ε_1}) error) would be exactly the kind of "explicit polynomial-in-A bound" required by the decision rule for SELECTED. But the alternative algebraic completion is itself an open problem — would require its own theorem. Mode H circular.
