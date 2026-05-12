# SL2_FRAMEWORK_TRANSFER — does the Furstenberg framework transfer to μ_n?

**Date:** 2026-05-12. Phase 3 of the SL_2(ℝ)-embedding structural-compatibility probe.

This file evaluates whether, under the only Phase-2 candidate that passes G1 (Candidate A: 2-atom
μ on SL_2(K) with E[M] = T_lead), the resulting Furstenberg measure ν on P^1(ℝ) bears a
relationship to the Syracuse stationary measure μ_n on ℤ_3 in a way that transfers polynomial
Fourier decay.

---

## 1. What Candidate A delivers, under the Furstenberg framework

Assume (for the sake of argument, given the Phase 2 caveats) that Candidate A's 2-atom
construction satisfies G1-G5:

  μ = (1/2) δ_{M_+} + (1/2) δ_{M_-} on SL_2(K), K = ℚ(√7)
  M_+, M_- ∈ SL_2(K), both elliptic (trace 43/45 each) or hyperbolic (per choice of Δ)
  E[M] = T_lead = (1/45)·[[7, 9], [28, 36]]

The Furstenberg framework then says there is a unique stationary ν on P^1(ℝ) satisfying

    ν = (1/2) (M_+)_* ν + (1/2) (M_-)_* ν.

If hypotheses G1-G4 hold and (Li 2020) G5 holds (Zariski-density), then

    |ν̂(k)| ≪ |k|^{−β}        for some β = β(μ) > 0.

The Frostman-dimension paper (2601.14061) gives an effective form of β if the Hochman-Solomyak
algebraic-data conditions are met sharply.

---

## 2. The transfer question: does this ν project onto μ_n?

This is the load-bearing question. Two structural mismatches:

### 2.1 Mismatch in fixed-point equations

ν satisfies the Furstenberg fixed-point equation:

    ν(B) = (1/2) ν(M_+^{-1} B) + (1/2) ν(M_-^{-1} B)           for Borel B ⊂ P^1.

μ_n (Syracuse stationary on ℤ_3) satisfies a DIFFERENT fixed-point equation derived from Tao's
recursion. Specifically, μ_n is the stationary measure of the **inverse Syracuse Markov chain
on ℤ_3** (per project context, R58 onwards). The transition kernel is:

    K_Syracuse(r → r') = ½ if (r, r') consistent with one parity branch of (3r+1)/2 mod 3^n; 0 else.

This is a Markov chain on ℤ_3 (a discrete profinite set), NOT a projective action on P^1(ℝ).

Even if we identify P^1(ℝ) with ℝ ∪ {∞} via stereographic projection, ν is a continuous measure
on ℝ, while μ_n is a measure on a discrete or profinite set ℤ_3. **The two measures live on
fundamentally different topological spaces.**

### 2.2 Mismatch in operator structure

T_lead = E[M] is the **linearized first-moment operator** of the cross-frequency-derived 2D
projection of Tao's bilinear recursion. The corresponding ν is the stationary measure of the
PROJECTIVE RANDOM WALK driven by μ.

But Tao's recursion is NOT the projective action of μ on P^1. Tao's recursion is a bilinear
operator on FOURIER COEFFICIENTS (the μ̂_n(ξ) for ξ ∈ ℤ/3^n), with the bilinear coupling coming
from convolution in physical space.

The cross-frequency closure (CROSS_FREQ_DISPOSITION) found that the closed family is V_M, not
the simpler 2D space (P_+, P_-). The T_lead in 2D is the (1, 4)-projection of T_V on V_M.
**ν** would be the Furstenberg measure of the SL_2 random walk generating T_lead — which lives
on P^1(ℝ). **This is not even on the same scale as μ_n.**

### 2.3 Even if a projection-style transfer existed

Suppose ν on P^1(ℝ) is related to μ_n via some Borel map π: P^1 → ℤ_3 (e.g., via 3-adic
expansion of x ∈ ℝ identified with (0, ∞) ⊂ P^1). Then π_* ν is a measure on ℤ_3, comparable
to μ_n in principle. **But:**

  - π is generically singular/non-Lipschitz from ℝ to ℤ_3 (different topologies; ℝ is
    archimedean, ℤ_3 is non-archimedean). Maps preserving Fourier decay between these are
    SPECIFIC algebraic-arithmetic constructions, not generic.

  - Even if a π exists, the pushforward π_* ν has Fourier decay |π_* ν̂(k)| ≪ ... bounded by the
    pullback under π:  π pulls back ℤ_3-characters to functions on P^1, which are NOT P^1-
    characters in any natural sense. The polynomial decay |ν̂(k)| ≪ |k|^{-β} on archimedean
    frequencies k ∈ ℝ has no direct translation to 3-adic frequencies of π_* ν on ℤ_3.

  - Even ignoring topological mismatches, the **convolution structure** of Furstenberg vs.
    Tao-bilinear differs. Furstenberg's ν satisfies a one-step convolution equation; Tao's μ_n
    satisfies a bilinear recursion. The two fixed-point structures don't naturally project onto
    each other.

---

## 3. Could T_lead be derived from ν via a different route?

A potentially cleaner direction: instead of asking "does ν project to μ_n?", ask "is T_lead's
spectral structure recoverable from ν?"

  - The leading Lyapunov exponent χ_1 of μ on SL_2(K) is related to ‖M_+‖, ‖M_-‖ growth. For
    the construction in §A.5 of Phase 2 with both atoms having trace 43/45 (elliptic), we have
    ‖M‖² ≈ 1 + |off-diagonal|² (elliptic rotation amplitude), so log ‖M_n ⋯ M_1‖ grows
    sub-linearly (often logarithmically) — i.e., χ_1 = 0 or χ_1 is small.

  - **χ_1 = 0 is the boundary case**: Furstenberg's positivity theorem REQUIRES non-elementary
    + finite moment ⟹ χ_1 > 0. Elliptic-generator setups fall on the negation of "non-
    elementary" (they're inside a compact subgroup conjugate to SO(2)). So **the construction in
    §A.5 of Phase 2 with elliptic atoms FAILS G2 (non-elementary) outright** and gives χ_1 = 0,
    no Furstenberg measure.

  - For the hyperbolic-atom variant (which Phase 2 §A.6 found requires entries in higher
    algebraic extensions and may not have rational lifts), the framework would apply, χ_1 > 0,
    and ν exists on P^1(ℝ). But this ν has no relation to μ_n.

  - Hochman-Solomyak's dimension formula dim ν = min{1, h_RW / 2χ_1} gives a number between 0
    and 1 for ν's Hausdorff dimension. T_lead's eigenvalue 43/45 is unrelated to χ_1 (T_lead is
    E[M], not the Lyapunov-exponent-controlling product). **T_lead's spectrum is not recoverable
    from ν's dimensional invariants.**

---

## 4. What polynomial-in-A bound WOULD result if the transfer worked?

Hypothetically, if Candidate A's construction gave a ν on P^1 with polynomial Fourier decay
|ν̂(k)| ≪ |k|^{-β}, AND there were a Lipschitz pushforward π: P^1 → ℤ_3 with positive Jacobian,
then π_* ν would inherit polynomial decay |π_* ν̂(k)| ≪ |k|^{-β} on 3-adic Fourier modes.

  - The "A" in Tao Prop 1.17 is the level-mass parameter (number of high-frequency Fourier
    coefficients summed over), connected to the iteration depth n through n = O(log A) or
    similar.

  - β would be the Furstenberg-measure decay exponent — typically irrational, dependent on
    the entries of M_+ and M_-, and **not polynomial in A**.

  - Even under best-case assumptions, β is a CONSTANT depending on μ, not a polynomial in A.
    The bound would be |μ̂_n(ξ)| ≪ |ξ|^{-β}, where |ξ| ≪ 3^n. Translating to A via the level
    counting: this gives `|μ̂_n(ξ)| ≪ A^{-β}` for A = number-of-frequencies-summed ≪ 3^n.

  - **A-polynomial decay |μ̂_n(ξ)| ≪ A^{-γ} for γ rational/polynomial in A** is the form
    Tao Prop 1.17's effective C_A wants. The Furstenberg-framework β IS polynomial-in-A (in the
    sense that |ξ|^{-β} is polynomial decay), but the EXPONENT β is a fixed irrational number,
    not improvable as a function of A.

This is actually still useful for Tao Prop 1.17 IF β > 0 with explicit value. But the issue is
that there's no recipe to compute β from T_lead alone — β depends on the FULL measure μ on
SL_2(K), not on E[M] = T_lead.

---

## 5. Honest summary of Phase 3

Even ASSUMING Candidate A's construction gives a valid Furstenberg-framework setup (which
requires non-elliptic atoms — a separate construction question), the resulting ν on P^1(ℝ) has
the following relations to μ_n:

| Question | Answer |
|---|---|
| Is ν = μ_n? | No — different measures, different spaces. |
| Is ν a projection of μ_n? | No natural projection exists. |
| Is μ_n a projection of ν? | No natural projection P^1(ℝ) → ℤ_3 preserving Fourier decay. |
| Is T_lead's spectrum recoverable from ν? | No — T_lead = E[M], ν = stationary measure of projective walk; different invariants. |
| Does polynomial decay of ν̂ imply polynomial decay of μ̂_n? | Only if the transfer-map T1 exists with positive Jacobian; no such map known. |

**Transfer gate T1 fails for the only Phase-2 candidate that passed G1.**

---

## 6. Comparison to He-de Saxcé alternative

For completeness, what would happen under He-de Saxcé (torus equidistribution)? The Syracuse
update (3x + 1)/2^v involves multiplication by 3 (acts on 3-adic side) and division by 2^v
(acts on 2-adic side). The joint action is **not in SL_2(ℤ)** — the determinant is 3 / 2^v,
which is not 1. So even the torus alternative fails its prerequisite.

---

## 7. Phase 3 verdict

**No transfer.** Even ASSUMING the best possible reading of Candidate A as a viable SL_2-
embedding satisfying the Furstenberg framework's hypotheses, the resulting Furstenberg measure
ν on P^1(ℝ) has no natural pushforward / pullback / projection relationship to the Syracuse
stationary measure μ_n on ℤ_3. The two measures live on different topological spaces (continuous
ℝ vs. profinite ℤ_3), satisfy different fixed-point equations (one-step projective convolution
vs. bilinear Fourier recursion), and have no shared invariants accessible via T_lead alone.

The Furstenberg framework's hypotheses might be satisfiable for a 2-atom μ extending T_lead,
but the FRAMEWORK CONCLUSION (polynomial decay of ν̂) does not transfer to μ̂_n via any
identifiable mechanism.
