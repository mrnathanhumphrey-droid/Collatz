# C2_THEOREM_5_9_HYPOTHESES — verbatim hypothesis statements

**Date:** 2026-05-12. Cluster 2 cut-and-project probe, Phase 1.

Sources: Richard–Strungaru 2017 "Pure point diffraction and Poisson summation" (arxiv 1512.00912) and "Short guide to pure point diffraction in cut-and-project sets" (arxiv 1606.08831).

There are TWO Theorem 5.9 statements relevant here, which I quote separately.

---

## (a) Cut-and-project scheme definition

**Source: 1606.08831 §5 (page 9):**

> "Let G = R^d and let H be a compactly generated LCA group. We write π_G : G × H → G, π_H : G × H → H for the canonical projections. Given a lattice L in G × H, the triple (G, H, L) is called a **cut-and-project scheme** if (i) π_G is one-to-one on L, and if (ii) π_H(L) is dense in H. Let us call these two conditions the projection assumptions."

**Source: 1512.00912 §5.3 (page 19):**

> "Let (G, H, L) be a cut-and-project scheme with σ-compact G…"

i.e. 1512.00912 generalises G from R^d to **σ-compact LCA G** (load-bearing for ℤ_3, which is compact hence σ-compact).

**Remark 5.1 (1606.08831, p. 9):**

> "(i) Given (G, H, L), we have that (Ĝ, Ĥ, L^0) is also a cut-and-project scheme. Indeed, π_Ĝ is one-to-one on L^0 if and only if π_H(L) is dense in H, and π_Ĥ(L^0) is dense in Ĥ if and only if π_G is one-to-one on L. This is a consequence of Pontryagin duality, see e.g. [45, Sec. 5].
> (ii) In the Euclidean setting G = R^d and H = R^n, assume that L is a rotated scaled copy of Z^{d+n}. Then either of the two projection assumptions implies the other by duality, compare Figure 2."

---

## (b) Window W and regular model set

**Source: 1606.08831 Definition 5.2 (page 9):**

> "Let a cut-and-project scheme (G, H, L) and a **regular window** W ⊆ H be given, i.e., a relatively compact and measurable set with non-empty interior such that θ_H(∂W) = 0. Then ⋏(W) = π_G(L ∩ (G × W)) is called a **regular model set**. If W is relatively compact and measurable, then ⋏(W) is called a **weak model set**."

So:
- **regular window** = relatively compact + measurable + non-empty interior + boundary Haar-measure zero.
- **regular model set** uses a regular window.
- **weak model set** drops the interior and boundary conditions; only relative compactness + measurability remain.

**Remark 5.3(ii) (1606.08831, p. 10):**

> "Since we may pass from H to π_H(L), even if the second projection assumption does not hold we can assume that π_H(L) is dense in H without loss of generality. A prominent example where the second projection assumption is violated is the Penrose point set, when projected from G × H = R^{2+3} and L a rotated copy of Z^5, compare [18] and [4, Rem. 7.8]."

---

## (c) Three load-bearing conditions

Combining Definition 5.2 and §5 (1606.08831):

1. **π_G injective on L** (projection assumption (i))
2. **π_H(L) dense in H** (projection assumption (ii))
3. **W regular**: relatively compact, measurable, non-empty interior, θ_H(∂W) = 0

Strict regular-model-set Theorem 5.9 requires all three.

---

## (d) Theorem 5.9 — regular case (1606.08831)

**Verbatim, page 12:**

> "**Theorem 5.9** (diffraction formula for regular model sets). Consider the Dirac comb ω_{1_W} for some regular window W ⊆ H in some cut-and-project scheme (G, H, L). Then ω_{1_W} has autocorrelation γ and diffraction γ̂ given by
>
>   γ = dens(L) · ω_{1_W * 1̃_W} ,    γ̂ = dens(L)² · ω_{|1̂_W|²} .   □"

**Remark 5.10(ii):** "The diffraction formula might no longer be valid if the second projection assumption is violated."

Mass at character χ ∈ π_Ĝ(L^0): intensity I(χ) = dens(L)² · |1̂_W(χ⋆)|², where (χ, χ⋆) ∈ L^0.

---

## (e) Theorem 5.9 — weak-model-set extension (1512.00912)

**Verbatim, page 21:**

> "**Theorem 5.9.** [3, Thm. 7] Let (G, H, L) be a cut-and-project scheme with σ-compact G and let ω_h be the Dirac comb of a weak model set, i.e., h = 1_W for some relatively compact measurable W ⊂ H. Assume that there exists a van Hove sequence (A_n)_{n∈N} in G such that ω_h has **maximal density** with respect to (A_n)_n, i.e.,
>
>   lim_{n→∞} (1/θ_G(A_n)) ω_h(A_n) = dens(L) · θ_H(W) .
>
> Then, with respect to the given van Hove sequence (A_n)_n, the weak model set ω_h has autocorrelation γ and diffraction γ̂ given by
>
>   γ = dens(L) · ω_g * g̃ ,    γ̂ = dens(L)² · ω_{|ĝ|²} ,
>
> where g = 1_W. In particular, ω_h has pure point diffraction. □"

**Remark 5.11 (1606.08831, p. 12) — same content, short form:**

> "Any weak model set ⋏(W) satisfies the inequality dens(⋏(W)) ≤ dens(L) · θ_H(W), which can be proved by approximating W from above using regular windows [29, Prop. 3.4]. We say that a weak model set has maximal density if we have equality in the above expression, for some fixed averaging sequence. Regular model sets are of maximal density by Theorem 5.7. The approximation argument can be used to show that Theorem 5.9 even holds for weak model sets of maximal density [7, Cor. 6], compare also [32, Sec. 3.3.2]."

Crucial: σ-compact-G Theorem 5.9 ONLY requires (a) lattice L in G × H, (b) some W relatively-compact and measurable, (c) maximal-density realization along some van Hove sequence. Strict π_H-density and W-regularity are NOT needed when maximal density holds.

---

## (f) Theorem 5.4 (1606.08831) — weighted-comb PSF analog

**Verbatim, page 10:**

> "**Theorem 5.4.** Let (G, H, L) be a cut-and-project scheme and h ∈ P_K(H). Then ω̂_h ∈ M(Ĝ) is a positive measure, and we have
>
>   ⟨ω_h, g⟩ = dens(L) · ⟨ω_ĥ, ĝ⟩
>
> for all g ∈ P_K(G). **This result holds without the projection assumptions on (G, H, L).**" (emphasis mine)

P_K(G) = positive-definite continuous compactly supported functions.

The **weighted Dirac comb** (page 10, Eqn. before Thm 5.4):

> "For h: H → C bounded and compactly supported, consider the weighted Dirac comb ω_h ∈ M_∞(G) defined by
>
>   ω_h = Σ_{(x,y)∈L} h(y) δ_x .
>
> If π_G is one-to-one on L, we may identify ⋏(W) with ω_{1_W}."

---

## (g) Theorem 5.7 (1606.08831) — density formula

> "**Theorem 5.7** (Density formula for weighted model sets). Let (G, H, L) be a cut-and-project scheme and h: H → C Riemann integrable. Then M(ω_h) = dens(L) · ∫_H h(y) dθ_H(y). **This result uses only the second projection assumption.** In particular if W ⊆ H is a regular window, we have M(⋏(W)) = dens(L) · θ_H(W)."

---

## (h) BMP super-singular model set (1996/2000) — adelic G × H

**Source: bmp_visible.txt, pages 34–36, "Further connections and directions":**

> "Above, we have emphasized that the sets of visible lattice points V_Γ and the set of k-th-power-free numbers F_k differ from any regular model set … by a set of positive density, suggesting that they cannot be obtained from the cut-and-project construction in any natural way. **However, there is a way of obtaining these sets by cut-and-project using the rational adeles instead of Euclidean space as the hyperspace** … From this point of view, V_Γ and F_k are 'super-singular' model sets. This comes about because these sets are the result of sieving over primes: F_k, for example, is what remains of Z after removing the zero residue class mod p^k for each p."

The BMP scheme for F_k (page 35):

>   R  ←π←  A_Q  →π_int→  ∏'_p Q_p
>            ∪
>            Q
>
> "with componentwise projections and ∏'(Z_p) denoting the restricted product. The image π(Q) is, of course, dense in R and the denseness of π_int(Q) in ∏ Q_p is equivalent to the **Strong Approximation Theorem** … If we choose for the window the closed and open set Ω := ∏_p Z_p we obtain the model set Λ(Ω) = Z ∈ R. … If instead we choose Ω := ∏_p (Z_p \ p^k Z_p) with k ≥ 2 we obtain the 'thin' model set Λ(Ω) = F_k. We say 'thin' because, being not relatively dense, F_k is not a regular model set. The reason that Schlottmann's result does not apply here is that **Ω has empty interior**: Ω contains no basic open set because it projects to a proper subset of Z_p in every non-Archimedean component."

Key takeaway: BMP's framework is **G = R, H = ∏'_p Q_p, L = Q diagonal**, and the "super-singular" window has empty interior so is *not* a regular window. F_k is a *weak* model set with positive density (the Möbius-mass calculation), and it IS pure-point diffractive — this is what fired in BMP Theorems 4–5. The σ-compact-G Theorem 5.9 of 1512.00912 covers exactly this configuration when maximal density holds.

---

## Bochner foundation (Theorem B.3, 1606.08831 Appendix B, paraphrased)

A bounded continuous positive-definite function f on G is the Fourier transform of a finite positive measure on Ĝ. The proof of Theorem 5.9 traces through Proposition 5.8 (autocorrelation is positive-definite by construction γ = ω_{1_W} ⊛ ω̃_{1_W}) and then applies Bochner.

---

## Summary — what Syracuse μ_n must satisfy

For the strict regular-model-set Theorem 5.9 (1606.08831 form):

- **G = LCA group where μ_n / Λ lives**
- **H = compactly generated LCA group** (the BMP "internal space")
- **L ⊂ G × H** a discrete cocompact subgroup
- **π_G injective on L**
- **π_H(L) dense in H**
- **W ⊂ H** regular window: rel. compact, measurable, non-empty interior, ∂W has Haar measure zero
- Set Λ identified with π_G(L ∩ (G × W)); Dirac comb ω_Λ

For the weak-model-set version (1512.00912 Theorem 5.9):

- G σ-compact LCA
- L lattice in G × H
- W ⊂ H rel. compact measurable (no regularity)
- van Hove (A_n) on G witnessing **maximal density**: lim ω_{1_W}(A_n)/θ_G(A_n) = dens(L) · θ_H(W)

For the weighted-comb version (Theorem 5.4):

- (G, H, L) cut-and-project (NO projection assumptions needed for Thm 5.4 itself)
- h ∈ P_K(H), i.e., positive-definite continuous compactly supported
- Then ω_h = Σ_{(x,y)∈L} h(y) δ_x has ω̂_h = dens(L) · ω_{ĥ}

These three forms are nested progressively weaker.
