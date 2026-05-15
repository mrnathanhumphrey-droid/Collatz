# FG Candidate C — Bourgain-Furman-Lindenstrauss-Mozes (BFLM) Theorem A

**PDFs:** BFLM_2007 (Comptes Rendus announcement), BFLM_2011 (JAMS full version).
**Extracted text:** `C:/tmp/fg/bflm_2007.txt`, `C:/tmp/fg/bflm_2011.txt`.

---

## BFLM 2011 Theorem A (VERBATIM, JAMS 2011, p. 232)

> "Let Γ < SL_d(R) satisfy (Γ-1) and (Γ-2) above, and let ν be a probability measure supported on a set of generators of Γ satisfying (1.1). Then for any 0 < λ < λ_1(ν) there is a constant C = C(ν, λ) so that if for a point x ∈ T^d the measure μ_n = ν^{*n} * δ_x satisfies that for some a ∈ Z^d \\ {0}
>     |μ̂_n(a)| > t > 0, with n > C · log(2 ||a|| / t),
> then x admits a rational approximation p/q for p ∈ Z^d and q ∈ Z+ satisfying
>     |x − p/q| < e^{−λ n} and |q| < (2||a||/t)^C."

### Standing assumptions (1.1), (Γ-1), (Γ-2):

> "(1.1)  ∫ ||g||^η dν(g) < ∞ for some η > 0" [finite exponential moment]
> "(Γ-1)  Γ is strongly irreducible in R^d (no finite union of proper subspaces is Γ-invariant)"
> "(Γ-2)  Γ contains a proximal element, i.e. an element with a dominant eigenvalue which is a simple root of its characteristic polynomial"

### Lyapunov definition (from p. 231):
> λ_1(ν) = lim_{n→∞} (1/n) log ||g_1 g_2 … g_n||, ν^{Z+}-a.s. (Furstenberg-Kesten Lyapunov exponent)

### Hypotheses (typed):

- h_C.A.group: Γ ⊂ **SL_d(R)** acting on T^d = R^d/Z^d. [TYPE (i) — real Lie group; semigroup of integer matrices]
- h_C.A.Γ-1: Γ strongly irreducible on R^d. [TYPE (i)]
- h_C.A.Γ-2: Γ contains proximal element with simple dominant eigenvalue. [TYPE (i)]
- h_C.A.walk_moment: ν on Γ with **finite exponential moment** ∫ ||g||^η dν < ∞. [TYPE (ii)]
- h_C.A.walk_gen: supp(ν) generates Γ. [TYPE (ii)]
- h_C.A.μ_n: μ_n = ν^{*n} * δ_x (n-step random walk on T^d starting from x). [TYPE (iii)]

### Conclusion C_C.A:

- Dichotomy on x: either |μ̂_n(a)| ≤ t, OR x has Diophantine bad approximation by rationals with controlled denominator.
- Quantitative parameters: C = C(ν, λ), λ < λ_1, n > C log(2||a||/t).

This delivers (via Cor C) a **polynomial-in-||b||/||t|| decay**: for x Diophantine generic, max_{0<||b||<B} |μ̂_n(b)| < B · e^{-c_2 n/M}.

---

## Phase 1 — hypothesis × input matrix

| Hypothesis | (1) Tao | (2) C1 | (3) R75/76/77 | (4) eps_k |
|---|---|---|---|---|
| h_C.A.group: Γ ⊂ SL_d(R) acting on T^d | **FAILED**: Syracuse chain has no SL_d(R) action; ambient is the abelian profinite (Z/3^n)*, not the torus T^d under a non-abelian matrix group action | FAILED | FAILED | FAILED |
| h_C.A.Γ-1: strongly irreducible | N/A (no R^d) | N/A | N/A | N/A |
| h_C.A.Γ-2: proximal element | N/A (abelian → no proximal element; all elements act by scalar, no contracting projective direction) | N/A | N/A | N/A |
| h_C.A.walk: ν on Γ, exp moment | N/A | N/A | N/A | N/A |
| h_C.A.μ_n: ν^{*n} * δ_x | The Syracuse μ_n is the n-step stationary measure on (Z/3^n)*, NOT a convolution-of-fixed-ν on a fixed torus. The chain *changes its modulus* with n. SATISFIED only if we re-interpret μ_n as the stationary distribution of a chain on the inverse limit, but BFLM requires fixed-d torus. FAILED. | FAILED | FAILED | FAILED |
| h_C.A.fixed dimension d | The Syracuse chain lives on (Z/3^n)* with **growing modulus 3^n**; BFLM's T^d is fixed-d. The category mismatch is the same as TAUBERIAN/F (Singha Roy LSD). FAILED. | FAILED | FAILED | FAILED |

**Phase 1: NO_FIT** on multiple hypotheses (group, action, fixed dimension).

---

## Phase 2 — conclusion shape

If hypotheses were satisfied: BFLM Theorem A gives a Diophantine dichotomy. Corollary C(1) (p. 232) translates this into a polynomial-in-||b||/||t|| bound:

> "Assume x is Diophantine generic in the sense that for some M and Q, |x − p/q| > q^{-M} for all q ≥ Q, p ∈ Z^d. Then for n > c_1 log Q, max_{b ∈ Z^d, 0 < ||b|| < B} |μ̂_n(b)| < B · e^{-c_2 n/M}."

This **IS** polynomial-in-A in spirit: the bound is exponentially decaying in n (the walk-step count), uniformly over Fourier modes in a B-window. With c_2 depending on the spectral gap and M on Diophantine class of x.

**Conclusion shape match: STRONG** — exactly the polynomial-in-A decay we seek. The obstruction is purely at the hypothesis side.

---

## Phase 3 — profinite extension

Could BFLM Theorem A extend to "Γ ⊂ SL_d(Z_3) acting on (Z_3/3^n Z_3)*"?

Steps in the proof (BFLM 2011, §3-§5):
1. **Diophantine analysis on T^d.** Uses the *continuous* torus structure and *real* Lyapunov exponent. The 3-adic analog would replace T^d with Z_3^d / lattice, and Lyapunov with 3-adic valuation. Some recent work (e.g., Benoist-Quint p-adic; Li 2018 hints, see FG_F) suggests this is technically possible but requires re-doing all geometric arguments in the p-adic setting.

2. **Sum-product on R / R^n.** BFLM uses Bourgain's discretized sum-product on R (Bou10) at multiple steps. The p-adic / non-archimedean discretized sum-product is studied (e.g., Bourgain-Gamburd-Sarnak in finite-fields setting), but the BFLM proof's specific use of *additive structure on R / R^n* doesn't directly transfer to (Z/3^n)*. NEEDS_TECHNICAL_WORK substantial.

3. **Crucial: Γ proximal on R^d.** A proximal element has a dominant eigenvalue. In the abelian Syracuse case (Z/3^n)*, every element acts by multiplication by a unit — every "eigenvalue" is just the unit itself, never dominant in any contractive sense. **There is no proximal element in the abelian setting.**

   This is the same category-of-object barrier hit in BT_DISPOSITION (Bruhat-Tits adelic): the chain (Z/3^n)* lacks the *non-commutative dynamics* that powers proximality.

4. **Single-element acting on whole torus.** BFLM's action g · x = g·x mod Z^d uses a single g ∈ SL_d(Z) for each step, with x ranging over T^d. The Syracuse step changes (n+1, x) → (n, x') by ξ → ξ · 2^{-v} mod 3^n, with v ~ Geom(2). The "v" being random with infinite support is structurally different from BFLM's "g ∈ Γ randomly with finite-support measure on Γ".

**Phase 3 disposition for BFLM**: STRUCTURALLY_BLOCKED. The proximality + non-commutative + fixed-dimension hypotheses are all categorically incompatible with the abelian profinite chain. The Diophantine-approximation argument (the *engine* of the proof) depends on the embedding x ∈ T^d, which has no analog in (Z/3^n)*.

---

## Disposition C: **NO_FIT** (categorical, with deep extension obstruction at Phase 3).

- Phase 1: multiple hypothesis failures (group, action, fixed-d, proximality).
- Phase 2: conclusion shape matches the target IF hypotheses were satisfied.
- Phase 3: extension structurally blocked by proximality + fixed-dimension.

**Surprise:** the Tao bound S_χ(n) ≪_A n^{-A} that Tao **proved** for Syracuse (1909.03562 Prop 7.1) **is** essentially the BFLM-style conclusion. But Tao's proof does NOT go through BFLM — it goes through a custom pair-grouping + Plancherel argument specific to the 2-adic chain. The category-of-object mismatch is real: BFLM's hypotheses give Tao-style decay on real toral actions, but Tao's proof on the profinite chain uses **different machinery**. The closure of c = 7/45's *rate* (not the decay itself) requires that machinery to be sharpened to give polynomial-in-A error, not just the qualitative decay.
