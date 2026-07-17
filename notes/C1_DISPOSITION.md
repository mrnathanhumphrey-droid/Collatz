# C1 Disposition — Syracuse μ̂_n(ξ) expressibility against Cluster 1

**Date:** 2026-05-12
**Probe:** Cluster 1 (Cochrane / BC / Heath-Brown / Konyagin / Kowalski) expressibility test for Syracuse μ̂_n(ξ).
**Pre-registered priors:** BC_DIRECT 15%, COCHRANE_MIXED 35%, HEATH_BROWN 10%, ABEL_REDUCTION 25%, DOESNT_FIT 15%.

---

## Disposition: **H_C1_DOESNT_FIT**

All four sub-cases (BC 4.7, Cochrane 1.2, Heath-Brown Heilbronn, Abel reduction) fail on Syracuse via independent structural mismatches that trace to a unifying meta-failure.

This was the 15%-prior outcome. Empirical posterior settled on it after Phase 1 made the tuple-space-vs-Z/3^n domain mismatch concrete and Phase 3b.iv exposed the degree-blowup mechanism for Cochrane 1.2.

---

## Load-bearing one-line summary per phase

**Phase 0 (PDF):** All 8 load-bearing UTF-8 transcripts already extracted at `C:/tmp/crystal/`; Tao 1909.03562 newly extracted via pypdf 6.10.2 to `C:/tmp/crystal/tao_1909.txt` (58 pp). No access blockers.

**Phase 1 (Tao recursion form):** Syracuse μ̂_n(ξ) = E χ(linear combination of 2^{-a_{[1,j]}}) is **an expectation over a Geom(2)^n tuple-space**, expressible (Tao §7.1 eq 7.5) as a **product over a 2D renewal walk of conditional expectations of an additive character**. The character χ is purely *additive*, no multiplicative χ enters; the summation domain is *tuple-space*, not Z/3^n or a multiplicative subgroup; the phase function is *2-adic exponential* in tuple variables.

**Phase 2 (Theorem hypotheses verbatim):**
- BC 4.7 needs q with "few prime factors" + p_i > q^ε ⇒ α_i < 1/ε; H a multiplicative subgroup of Z_q^*. For q = 3^n: α_1 = n → ∞ violates p_1 > q^ε for any fixed ε.
- Cochrane 1.2 / Cor 1.1: complete sum on Z/p^m with multiplicative-additive factored integrand, polynomial degrees d_1, d_2 *bounded uniformly in m*.
- HB Heilbronn Thm 1: phase n^p mod p^2 — base of the power equals the modulus prime; Fermat-quotient lift load-bearing.
- Shparlinski Problem 1 (multiplicative-character BGK): OPEN; doesn't bite directly because Syracuse uses additive.

**Phase 3 (sub-fits):**
- 3a BC 4.7: NO. (i) q = 3^n violates few-prime-factors-with-p_i>q^ε; (ii) Syracuse domain is tuple-space, no multiplicative-subgroup structure on (Z/3^n)^*.
- 3b Cochrane 1.2: NO. (i) Full μ̂_n: weight 2^{-Σa_i} not unimodular multiplicative; (ii) Tao's f(x,3): 2-term sub-trivial sum; (iii) R78-style polynomial identification: d ≍ n, bound degenerates to trivial.
- 3c Heath-Brown: NO. Phase 2^{-a} mod 3^n (coprime base 2 vs modulus prime 3) ≠ n^p mod p^2 (matching base); no Fermat-quotient lift.
- 3d Abel reduction: NO. Partial summation reshuffles within tuple-space; cannot change domain category from tuple-space to Z/p^m complete sum.

**Phase 4 (R78 D=0 disambiguation):**
- R78 tested Cochrane Theorem 2 (coefficient extraction, complete-sum vanishing) on Kalafatelis (1+3)^u family. D=0 obstruction was real and specific to Theorem 2.
- Cochrane Theorem 1.2 fails on Syracuse via a *different* mechanism: degree blowup d_1+d_2 ≍ n makes the (1.17) bound degenerate to trivial.
- Both obstructions trace to the algebraic root **4 = 1+3 in Z_3**: makes natural polynomial expansion p-adically trivial at the top (R78 D=0) AND p-adically degree-large at the lift (Theorem 1.2 degeneration).
- Audit of "Mode E inherited-claim trap": confirmed Theorem 2's D=0 doesn't apply wholesale to Theorem 1.2, but Theorem 1.2 fails on a different ground.

**Phase 5 (adversarial check):**
- (A1) **polynomial-in-A vs poly-in-p^m:** Cochrane 1.2 bound is poly in p^m with exponent 1-1/(d_1+d_2). For poly-in-A, need d_1, d_2 bounded in A; Syracuse forces d ≍ n. **Fails.**
- (A2) **uniformly in a:** HB's uniformity in a coprime to p preserved in Tao's bookkeeping only conditional on HB applying — which it doesn't.
- (A3) **R78 inherited-claim audit:** done in Phase 4. R78 D=0 specific to Theorem 2; Theorem 1.2 failure is structurally distinct but same algebraic root.
- (A4) **Shparlinski Problem 1:** doesn't bite — Syracuse is additive-character, closed by BC. The BC route's failure is the few-prime-factors hypothesis on q = 3^n.
- (A5) **A_v vs μ̂_n:** project-internal A_v (R75/R77 one-step peel) is a *single phase coefficient* at fixed v, not an exp-sum itself. The exp-sum character lives in the outer Σ_v which is a 1D finite (geometric-tail) sum. Bounds on A_v are trivial (|A_v| ≤ 1); the load is on the outer-Σ_v telescoping, which is Tao §7.4's actual battleground, not a Cluster 1 setting.

---

## Failure mode: category-of-object mismatch (arithmetic version)

Same meta-pattern as the five-probe Fourier-decay arc, restated for Cluster 1 arithmetic methods:

- **The five-probe arc closed** modern continuous-Fourier-decay machinery (BKS, Furstenberg, ARHW, drift) on category-of-object grounds — those frameworks target continuous/smooth-dynamical or transient-mixing objects.
- **This Cluster 1 probe closes** modern discrete-arithmetic exponential-sum machinery (BGK / BC / Cochrane / HB) on a parallel category-of-object grounds — those theorems target *single complete sums on Z/(p^m) or its multiplicative subgroups with polynomial / rational / Fermat-quotient phase functions*.
- **Syracuse μ̂_n(ξ) is structurally** a *product over a 2D renewal walk* indexed by a Geom(2)^n tuple-space, with phase function 2-adic exponential in tuple variables. **Neither category fits.**

Mode H (target-object selection) — even within Cluster 1 the three candidate target-objects (full μ̂_n / project-internal A_v / Tao's f(x,b) conditional expectation) **all** fail to fit Cluster 1 by independent structural reasons. There is no target-object choice that rescues Cluster 1 directly.

---

## Routing forward

Cluster 1 direct attack closed-negative. Per the brief and `POLYNOMIAL_IN_A_LANDSCAPE.md`:

- **Tauberian arc (primary):** Flajolet-Sedgewick / Chevalier 2507.15394 Thm 1.16 candidates. Operates on generating series E(z) = Σ ε_n z^n / Plancherel-trace Σ_k |π̂_k(ξ)|² — **different object**, escapes the category trap. Gated on ε_7 → ε_8 exact-rational compute (R77.7 v2 solver redesign).
- **Cluster 2 (BMP / PSF / cut-and-project) (secondary, parallel to Tauberian):** scope separately. The "different category" route — discrete-arithmetic Poisson summation / cut-and-project, native to Z/3^n structure.
- **Bourgain-Konyagin discrete sum-product on the orbit {2^v mod 3^n}** stays pre-cleared categorically. The *orbit* of 2 is the multiplicative subgroup ⟨2⟩ = (Z/3^n)^* (size 2·3^{n-1}), but the sum we need is not Σ_{x ∈ ⟨2⟩} e_q(ξ x) — that's BC 4.7's target which fails (Phase 3a). The right setting on the orbit is multiplicative-energy / sum-product on the *additive translates* of ⟨2^d⟩ subgroups, which is a separate scoping question.

**Recommendation:** consolidate this probe + the five-probe arc into a paper-section "**Category-of-object barriers to Syracuse Fourier decay**" — the negative result is itself paper-worthy: modern Fourier-decay AND modern exponential-sum machinery both fail on Syracuse for *parallel* category reasons. The closure path is then forced through Tauberian + Cluster 2.

---

## Deliverables in this probe

- `C:/Collatz/C1_TAO_RECURSION_FORM.md` (Phase 1 verbatim)
- `C:/Collatz/C1_THEOREM_HYPOTHESES.md` (Phase 2 verbatim)
- `C:/Collatz/C1_FIT_ATTEMPTS.md` (Phase 3, four sub-fits)
- `C:/Collatz/C1_R78_DISAMBIGUATION.md` (Phase 4 inherited-claim audit)
- `C:/Collatz/C1_DISPOSITION.md` (this file)

Tao 1909.03562 newly extracted to `C:/tmp/crystal/tao_1909.txt`.

No git operations performed.
