# C4_REPROBE_V2_HYPOTHESIS_CHECK

**Date:** 2026-05-14. Phase 2 hypothesis matching: do Cébron 2013 §2.7 / §4.4 + Goldsheid-Margulis 1989 §1, §5, §6 close the C4 (Ch 11 renewal-Egorov) gap for Syracuse's iterated Tao recursion?

Inputs: `C4_REPROBE_V2_CEBRON_HYPOTHESES.md`, `C4_REPROBE_V2_GOLDSHEID_MARGULIS_HYPOTHESES.md`, `C1_TAO_RECURSION_FORM.md`, `C4_REPROBE_TAO_RMT_DISPOSITION.md`, `result_77_T_lead_spectrum.md`.

---

## The C4 target — Syracuse's iterated Tao recursion

From `C1_TAO_RECURSION_FORM.md` (Tao 1909.03562 §7.1, eq 7.2 + eq 7.5):

> S_χ(n) := E χ( 2^{-a_1} + 3·2^{-a_{[1,2]}} + ⋯ + 3^{n-1}·2^{-a_{[1,n]}} )
>        ≤ E ∏_{j ∈ [n/2]} | f( 3^{2j-2} 2^{-b_{[1,j]}}, b_j ) |        (after pair-grouping)

where (a_1, …, a_n) ≡ Geom(2)^n iid, χ is the **additive character** Z[1/2] → C, x ↦ e^{−2πi ξ (x mod 3^n) / 3^n}, and the inner f is itself an expectation conditional on a pair-sum.

The project-internal one-step-peel form (R75/R77):

> μ̂_{n+1}(ξ) = Σ_v 2^{-v} A_v(ξ) μ̂_n(2^{-v} ξ mod 3^n)

The renewal-Egorov question (Ch 11): does there exist a composition formula at the OPERATOR level that controls the iterate T^n, with **off-diagonal v ≠ v' cross-frequency bilinear terms** treated explicitly?

---

## H1. Does Cébron §2.7 (Theorem 2.13) apply to Syracuse's iterated recursion?

**Theorem 2.13 statement:** τ(P(A B)|B) = (e^{D_A} P)(B) for A, B free in (A, τ), τ(A_i), τ(B_i) ≠ 0, P ∈ C{X_i : i ∈ I}.

### Match table — Theorem 2.13 hypotheses vs. Syracuse object

| Cébron hypothesis | Syracuse iterated Tao recursion | Match? |
|---|---|---|
| (A, τ) a W*-probability space | The Syracuse "step" is the operator 2^{-v} A_v(ξ) Δ_v on L²((Z/3^n)*), NOT a non-commutative random variable in a W* algebra. The trace τ on the natural von Neumann algebra of operators on L²((Z/3^n)*) IS well-defined (finite-dim case, just normalized matrix trace), but it does not have the "free probability" structure that L_κ free log-cumulants require. | NO (structural — see below) |
| A, B FREE in (A, τ) | The step operators at iteration level j and j+1 are NOT free. They share the multiplicative-unit action of (Z/3^n)*, which couples them. In Tao's eq 7.5, the inner-factor f at j and at j+1 share the variable b_{[1,j]} (and j+1 multiplies by 3 the j-th cyclic level). This is **explicit non-freeness via shared 2-adic exponent**. | NO — **freeness is exactly what fails**. R77's v ≠ v' bilinear coupling = non-freeness measure. |
| Non-zero traces τ(A_i), τ(B_i) ≠ 0 | The single-step operator 2^{-v} A_v(ξ) has trace ≠ 0 generically; not the obstruction. | YES |
| P ∈ C{X_i : i ∈ I} (formal polynomial in non-commuting indeterminates with traces) | Syracuse's iterated map is NOT a polynomial in formal indeterminates with trace — it is a chain of conditional expectations over Geom(2)-distributed integer exponents. The structural object is an **expectation along a renewal walk**, NOT a polynomial expression evaluated at iid free random variables. | NO — category mismatch. The trace-polynomial calculus C{X_i} captures the moment structure of free variables; the Syracuse renewal walk indexes by integers, not by non-commutative random variables. |

### Verdict on H1: Theorem 2.13 does NOT apply directly

The composition formula τ(P(AB)|B) = (e^{D_A} P)(B) is **explicitly conditional on A and B being free**. In the Syracuse iteration, "A" and "B" would have to be the j-th and j+1-th step operators, and **these are NOT free** — they are arithmetically coupled through the common 2-adic exponent structure. This is the **same non-freeness obstruction** that the v1 Tao-RMT probe identified.

Cébron supplies a multiplicative composition formula **only for free A, B**. It does not supply a substitute formula for **non-free A, B**.

---

## H2. Does Cébron's large-N limit (§4.4 Theorem 4.6) cover the Syracuse N = 3^n → ∞ regime?

**Theorem 4.6 statement:** E[tr(P_0(G^{(N)}_t)) ⋯ tr(P_n(G^{(N)}_t))] = τ(P_0(G_t)) ⋯ τ(P_n(G_t)) + O(1/N²).

### Match table

| Theorem 4.6 hypothesis | Syracuse N = 3^n setting | Match? |
|---|---|---|
| G^{(N)}_t = right Brownian motion on **smooth Lie group GL_N(C)** | Syracuse's iterated operator acts on **finite cyclic group (Z/3^n)\***, a **profinite (discrete) group**. There is no smooth Lie structure; (Z/3^n)\* is a finite cyclic group of order 2·3^{n-1}. | NO — categorical mismatch (smooth Lie ↔ profinite/discrete) |
| Asymptotic parameter: matrix size N → ∞ | Syracuse asymptotic parameter: cyclic level n → ∞, i.e. N = 3^n. **The role of N is the same in a numerical sense** (a dimension going to infinity), but the underlying geometry is different. | PARTIAL — N → ∞ matches the LIMIT FORM but not the underlying object |
| Polynomial in (X, X*, X⁻¹, X*⁻¹) | The Syracuse iterate is a sum of phase-factor expectations, **not** a polynomial in operator-valued indeterminates. After one-step-peel R75, the form is Σ_v 2^{-v} A_v Δ_v, which is a **shift sum**, not a polynomial in non-commuting variables. | NO |
| Rate: O(1/N²) | For Syracuse this would be O(1/3^{2n}), which is exponentially small. There is no a-priori reason to expect this rate. The rate-½ R77 conjecture predicts spectral radius 1/√3 = O(3^{-n/2}), much weaker than O(3^{-2n}). | NO — rate scales are incompatible |

### Verdict on H2: Large-N limit does NOT transport directly

The Cébron / Biane large-N convergence is a **smooth-Lie-group-to-W*-algebra approximation**. Its driving mechanism (Theorem 4.6's proof, p. 50-51) is the **smoothness of the exponential map on End(C^d{...})**, which collapses the (1/N²) correction Δ̃_{GL} on polynomials. This **explicitly uses the smooth structure** of GL_N(C).

Syracuse's (Z/3^n)* has **no smooth structure**: it is a discrete cyclic group. The "exponential map" of a Lie algebra is replaced by the discrete-log structure of the multiplicative group, which is **not** a smooth perturbation of any limit object.

---

## H3. Does Goldsheid-Margulis Theorem 1.2 / 5.4 give the leading Lyapunov exponent for Syracuse?

**Theorem 1.2 statement:** Under ln^+‖A‖ ∈ L_1, the iid product A(n, ω) = A_n ⋯ A_1 satisfies lim n^{-1} ln μ_i(n, ω) = χ_i(ω) almost surely.

### Match table

| GM Thm 1.2 hypothesis | Syracuse iterated Tao recursion | Match? |
|---|---|---|
| A : Ω → GL(m, R) measurable, **finite-dimensional matrices** | The Syracuse step at level n is an operator on L²((Z/3^n)*), which IS finite-dimensional (dim 2·3^{n-1}). So at fixed n it's an iid random matrix on R^{2·3^{n-1}}. But **the level n changes with the iteration**: each step takes a level-n object to a level-(n+1) object via multiplication by 3 on the cyclic structure. This is **NOT iid in a fixed GL(m, R)**. | PARTIAL — works in fixed-level approximation; fails for the level-changing iteration |
| ln^+ ‖A(ω)‖ ∈ L_1 | The Geom(2)-weighted single-step operator has log-norm bounded by O(n) (since the cyclic group order is 3^n); under the geometric weighting the integrability is fine. | YES at fixed level n |
| iid in fixed GL(m, R) | At fixed level n, the step operators **are** iid in a fixed GL(2·3^{n-1}, R) (the Geom(2) randomness is iid). But the relevant limit is n → ∞, NOT iteration-count k → ∞ at fixed n. | PARTIAL — same as above |

### What GM does close

GM **gives the leading scalar Lyapunov index** in the fixed-level n setting: at fixed n, the iterated operator T^k on L²((Z/3^n)*) has a leading Lyapunov exponent (depending on n) given by lim k^{-1} ln ‖T^k x‖ a.s.

In the n → ∞ limit (which is the Syracuse limit), this Lyapunov exponent should converge to **λ_Syracuse = log 3 − 2 log 2 ≈ −0.288** (the known scalar leading rate, per the project's standing result).

### Verdict on H3: GM gives the SCALAR leading-order exponent; nothing more

GM closes the **scalar Lyapunov-index** question for the fixed-level n problem. It does **NOT** supply:
- An operator-valued composition formula (it gives a scalar exponent, not the spectral measure of the iterated operator)
- Cross-frequency v ≠ v' coupling control (Oseledets handles iid products, and the conclusion is about norms of A(n)x for fixed vectors, NOT about the spectral measure of A(n))
- The n → ∞ (cyclic level → ∞) limit — only the k → ∞ (iteration count → ∞) limit at fixed n

This is the **structural ceiling of GM**: leading scalar exponent only. **It is upstream of Ch 11, not at the level of Ch 11.**

---

## H4. Does Cébron's framework provide a "freeness escape" via operator-valued / amalgamated free probability?

Cébron works in **scalar-valued** (A, τ) — a single W*-probability space with single trace τ. The paper does NOT develop:
- **Operator-valued** free probability (Voiculescu's amalgamated free probability over a subalgebra B ⊂ A)
- **B-free independence** (the natural generalization for non-free A, B that are free over a common subalgebra)
- Any other generalization that would relax freeness

These exist in Voiculescu's later work (Voiculescu 1995 "Operations on Certain Non-Commutative Operator-Valued Random Variables", or Speicher 1998 Memoirs AMS), but they are **NOT in Cébron 2013**.

### Verdict on H4: NO. Cébron stays scalar-valued; the freeness escape is OUTSIDE this paper.

The freeness obstruction identified in v1 (Tao-RMT) remains. Cébron supplies a strictly cleaner multiplicative composition formula than Tao RMT (Theorem 2.13 is the explicit statement Tao RMT punts to Speicher), but the FREENESS hypothesis is still load-bearing.

---

## H5. Does Cébron's framework port from GL_N(C) to (Z/3^n)*?

(Detailed in `C4_REPROBE_V2_PROFINITE_PORT.md`.)

Short answer: the load-bearing steps that USE the smooth Lie group structure of GL_N(C) are:
- The Brownian motion construction (§3.3, §3.5): requires a **Laplacian on a smooth manifold**.
- The Itô calculus identities (eq 3.9, eq 3.15, eq 4.6, eq 4.7): require **smooth differential operators**.
- The large-N convergence proof (§4.4 Thm 4.6): uses the **smoothness of exp on End(C^d{X_i})** to drop O(1/N²) corrections.

None of these has a direct profinite analog. The W*-probability framework (§1.1, §2.1) is **abstract** and might be portable; but the **specific composition formula e^{D_A}** is built from free log-cumulants on Hopf algebra Y(k), which presupposes the **scalar-valued free-probability** structure. There is no a-priori obstruction to defining an analogous Hopf algebra over a profinite base, but the construction is original (research-monograph level).

---

## Summary match table

| C4 / Ch 11 requirement | Cébron supplies? | GM supplies? |
|---|---|---|
| Composition formula for **products** of i.i.d. operators | YES (Thm 2.13) — but **only for FREE A, B** | NO (gives scalar exponent, not composition) |
| Renewal-product iterate E[∏_j ...] | YES in operator-valued form (Thm 2.13 + Prop 2.12 chain), but freeness required at every step | NO |
| n → ∞ (cyclic level) iteration limit | NO — Cébron's "N → ∞" is matrix size, not cyclic level | NO — GM iterates count k → ∞ at fixed level |
| Cross-frequency v ≠ v' off-diagonal bilinear coupling | NO — Thm 2.13 assumes freeness, which is exactly the absence of this coupling | NO — scalar exponent absorbs cross-frequency into its limit |
| Profinite / discrete base | NO — built on smooth GL_N(C) | NO — built on GL(m, R) ↪ ergodic dynamical system |
| Operator-valued / amalgamated extension to handle non-free A, B | NO — paper is scalar-valued only | NO |
| Leading-order Lyapunov scalar | NO (Cébron doesn't address) | YES — Thm 1.2 gives λ_Syracuse = log 3 − 2 log 2 |

---

## Disposition (preview; full statement in C4_REPROBE_V2_DISPOSITION.md)

**Cébron Theorem 2.13 is the EXACT multiplicative composition formula Tao RMT punted to Speicher.** It is the cleanest scalar-valued multiplicative-free-convolution theorem in the corpus. But it **requires freeness**, which Syracuse's R77 cross-frequency coupling explicitly violates.

**Cébron's large-N limit (Theorem 4.6) is N = matrix size → ∞**, with the convergence driven by **smoothness of GL_N(C)**. This is NOT the Syracuse "N = 3^n → ∞" cyclic-level regime, which has no smooth structure.

**Goldsheid-Margulis Theorem 1.2 gives the scalar leading Lyapunov exponent** under iid product structure on GL(m, R). This is upstream of Ch 11: it closes the SCALAR leading-order question (already known: λ_Syracuse = log 3 − 2 log 2), not the operator-valued composition.

**Net status:** the corpus delivers a **sharper version of the v1 finding** (multiplicative composition formula is now explicit, not punted), but **the freeness obstruction + profinite portability ceiling remain**. C4 status is **C4_PARTIAL_OUTSIDE**: the needed next step is **operator-valued / amalgamated free probability** (Voiculescu 1995 / Speicher 1998), which is **not in this corpus**.
