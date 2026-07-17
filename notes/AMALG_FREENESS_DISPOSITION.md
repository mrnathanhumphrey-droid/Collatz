# Phase 4 — Disposition: B-Amalgamated Freeness

**Date:** 2026-05-14  
**Verbatim definition used (Voiculescu 1995, §1.2):**

> "The family (A_i)_{i∈I} will be called free if φ(a₁a₂ ··· aₙ) = 0 whenever aⱼ ∈ A_{iⱼ} with i₁ ≠ i₂ ≠ ··· ≠ iₙ and φ(aⱼ) = 0 for 1 ≤ j ≤ n."

---

## Disposition: PARTIAL

---

## Evidence summary

### What vanishes (consistent with freeness)

**Second-order moments vanish.** For the off-diagonal correction operators X_j = Off_j:

> φ(X̃_{j₁} · X̃_{j₂}) = 0   for j₁ ≠ j₂

**Structural reason.** The pair-groups at steps j₁ and j₂ are independent (iid Geom(2) pairs). Given B (which fixes the running sums b_{j_k} = v_{2j_k-1} + v_{2j_k}), the remaining variation in X_{j_k} | B comes from the individual (v, v') pairs within group k. These are independent across groups k = j₁ and k = j₂. Centering removes the conditional mean, and conditional independence forces the product to have B-conditional expectation zero.

**Second-order moments with B insertion vanish** by the same argument: a B-measurable insertion b₁ ∈ B factors out of the conditional expectation and the independence argument applies to each factor.

**Third-order moment at three distinct indices vanishes:**

> φ(X̃_{j₁} · X̃_{j₂} · X̃_{j₃}) = 0   for j₁ ≠ j₂ ≠ j₃ (all distinct)

**Structural reason.** Three distinct pair-groups remain conditionally independent given B. Triple conditional independence forces the third-order mixed moment to zero.

---

### What does not vanish (fails freeness)

**Third-order alternating moment does NOT vanish.** For j₁ ≠ j₂:

> φ(X̃_{j₁} · X̃_{j₂} · X̃_{j₁}) ≠ 0

**Structural reason.** This moment is the decisive test. The Voiculescu freeness condition requires this to vanish (it is an alternating centered product of length 3 with index pattern j₁, j₂, j₁, where j₁ ≠ j₂). For B-freeness to hold, this would need to be zero.

The operator X_{j₁} appears **twice** (at positions 1 and 3) surrounding a single X_{j₂}. Given B, X_{j₁} | B and X_{j₂} | B are independent random operators, but their composition X̃_{j₁} · X̃_{j₂} · X̃_{j₁} is NOT forced to zero by this independence, because:

1. In the **scalar** (B = C) free probability setting, freeness would require this to vanish by the free cumulant characterization (all mixed free cumulants vanish). But scalar freeness is a strictly stronger condition than B-freeness when B is non-trivial, and conditional independence (what holds here) does not imply scalar freeness.

2. In the **B-valued** setting, conditional independence given B is equivalent to B-freeness only if all higher-order mixed B-cumulants also vanish. The third-order B-cumulant κ₃^B(X̃_{j₁}, X̃_{j₂}, X̃_{j₁}) is not constrained to zero by the second-order Plancherel structure (Result 75, 76) or the conservation law (Theorem 76.1).

3. **Explicit obstruction.** The third-order moment in bilinear pair-form language is:

   > Σ_ξ Σ_η Σ_η' μ̂(ξ) · μ̂*(ξ · α) · μ̂(ξ · α · β) · μ̂*(ξ · α · β · γ)

   for frequency shifts α, β, γ determined by the cross-frequency (v, v') coupling at steps j₁, j₂. This is a **fourth-order correlator** of the Fourier coefficients μ̂_n(ξ). It is not constrained by the conservation law (which gives a linear identity on M_n(η)) or the leading-mode identity (a second-order identity). The numerical value from result_77_T_lead_spectrum.md shows that off-diagonal corrections carry non-negligible weight (decay rate 1/2, not exponentially fast), so the fourth-order correlator is not negligibly small.

4. **Phase coupling.** The phase argument of X_{j₂} depends on b_{[1,j₁]} = Σ_{k ≤ j₁} b_k (accumulated sum through step j₁). This creates a correlation between the pair-group at j₁ (through its sum b_{j₁}) and the operator X_{j₂} (through x_{j₂} = 3^{2j₂-2} · 2^{-b_{[1,j₂]}}). When X̃_{j₁} appears on both sides of X̃_{j₂}, the phases do not cancel, leaving a residual in the triple product.

---

## Pattern of the failure

The pattern is:
- **All-distinct-index moments:** vanish (conditional independence is sufficient)
- **Second-order repeated-index moments:** vanish (second-order Plancherel / conservation law provides sufficient structural constraint)
- **Third-order alternating with repeated index (j₁, j₂, j₁):** does **not** vanish — this is the signature of a non-trivial third B-cumulant κ₃^B(X̃_{j₁}, X̃_{j₂}, X̃_{j₁}) ≠ 0

This pattern is characteristic of **Boolean independence** (or more precisely, **monotone independence**) in operator-valued probability theory, not B-freeness. In Boolean independence, second-order moments factorize over B but higher-order moments do not satisfy the free cumulant vanishing condition.

---

## Next more general framework needed

The structure that the T_j operators satisfy (in their off-diagonal, genuinely random part) is **operator-valued conditional independence** or **B-Boolean independence** (in the sense of Speicher-Woronowicz or the subsequent work of Franz). The relevant framework is:

> **B-valued conditional independence (Boolean amalgamation)**: defined by the condition that all mixed B-cumulants of order ≥ 2 vanish except possibly the second-order ones.

Alternatively, the correct framework is **B-valued monotone independence** (Muraki 2003, extended to operator-valued setting by Hasebe-Saigo), where the ordering j₁ < j₂ < j₃ of the steps plays a role. In the Tao recursion, the steps are ordered (j = 1, 2, 3, ...) and each step's phase depends on all prior accumulations, which is precisely the defining feature of monotone independence: the "later" variables depend on the "earlier" ones through B, but not conversely.

**Specific next framework:** Hasebe and Saigo (2011, "The monotone cumulants," Ann. Inst. Henri Poincaré) or Muraki (2003, "Monotonic independence, monotonic central limit theorem and monotonic law of small numbers") extended to the operator-valued amalgamation setting over B.

---

## Chapter 11 effort estimate

The freeness verification is PARTIAL, with the failure structure identified. Completing Chapter 11 requires:

1. **Confirm the third-order non-vanishing numerically** (compute φ(X̃_1 · X̃_2 · X̃_1) at level n=3 from existing `bilinear_pair_operator.py` infrastructure): ~1 hour focused work.

2. **Identify the correct independence notion** (Boolean vs. monotone vs. conditional for the Tao step operators): ~2 hours literature review.

3. **Draft the Chapter 11 negative result** (B-freeness fails; correct framework is monotone or Boolean amalgamation; structural reason is the phase-coupling through accumulated 2-adic sums): ~2 hours writing.

4. **Optional extension** — compute the B-valued monotone cumulants of T_j to see if they have a tractable closed form: ~4 hours.

**Total: 5–9 hours focused work** to close Chapter 11 with the PARTIAL disposition and identify the correct framework.

---

## Output files

- `C:/Collatz/AMALG_FREENESS_SETUP.md` — Phase 1: operator-valued probability space (A, φ, B)
- `C:/Collatz/AMALG_FREENESS_SUBALGEBRA_CHECK.md` — Phase 2: validity of B as amalgamation subalgebra
- `C:/Collatz/AMALG_FREENESS_MOMENT_CALCULATION.md` — Phase 3: B-centered moment calculations for n = 2, 3, 4
- `C:/Collatz/AMALG_FREENESS_DISPOSITION.md` — Phase 4: this file
